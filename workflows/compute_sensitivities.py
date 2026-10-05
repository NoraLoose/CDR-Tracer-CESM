#!/usr/bin/env python

import xarray as xr
import s3fs
import numpy as np
import cftime
import glob
import os
import re
from multiprocessing import Pool
from datetime import timedelta
import argparse

from carbonate_sensitivities import carbonate_sensitivity

# -----------------------------
# Parse arguments
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument(
    "--config",
    choices=["monthly", "daily", "daily-monthly"],
    required=True,
    help="monthly: SMYLE-FOSI monthly tseries.  daily: cdr-atlas-v0 control daily "
         "history files, beta/eta evaluated per day.  daily-monthly: the SAME daily "
         "control run but first averaged to monthly means, so beta/eta come from "
         "monthly-mean ALK/DIC/T/S/nutrients (f(mean), the same way "
         "compute_truth_sensitivities.py builds the truth _altco2 baseline).",
)
args = parser.parse_args()

CONFIG = args.config

# daily-monthly reuses the daily loader, then resamples to monthly means
RESAMPLE_TO_MONTHLY = CONFIG == "daily-monthly"

# -----------------------------
# Configurations
# -----------------------------
if CONFIG == "monthly":

    file_suffix = "monthly"
    USE_S3 = True
    TIME_CHUNK = None  # process all at once

    S3_DATA_FILE = (
        "s3://us-west-2.opendata.source.coop/cworthy/"
        "oae-efficiency-atlas/data/control/"
        "g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.SMYLE.005.pop.h.TEMP.030601-036812.nc"
    )

    selectors = {
        "time": slice(480, 492 + 12 * 22),
        "z_t": 0,
    }

    VAR_MAP = {
        "temp": "TEMP",
        "salt": "SALT",
        "ALK": "ALK",
        "DIC": "DIC",
        "PO4": "PO4",
        "SiO3": "SiO3",
    }

else:

    file_suffix = "monthly_from_daily" if RESAMPLE_TO_MONTHLY else "daily"
    USE_S3 = False
    # PyCO2SYS batch size (time steps).  daily: 30 daily frames per chunk.
    # daily-monthly: the monthly means are already in memory by this point, so
    # the chunk just bounds PyCO2SYS's own working set.
    TIME_CHUNK = 24 if RESAMPLE_TO_MONTHLY else 30

    LOCAL_PREFIX = (
        "/global/cfs/projectdirs/m4746/Users/nora/Ocean-CDR-Atlas-v0/"
        "data/archive/smyle.cdr-atlas-v0.control.001_backup/ocn/hist/"
        "smyle.cdr-atlas-v0.control.001.pop.h.nday1"
    )

    LOCAL_SUFFIX = "nc"

    selectors = {
        "z_t": 0,
    }

    VAR_MAP = {
        "temp": "SST",
        "salt": "SSS",
        "ALK": "ALK",
        "DIC": "DIC",
        "PO4": "PO4",
        "SiO3": "SiO3",
    }


# -----------------------------
# Load dataset
# -----------------------------
def load_data():

    ds = xr.Dataset()

    if USE_S3:

        fs = s3fs.S3FileSystem(anon=True)

        for key, var in VAR_MAP.items():

            s3_path = S3_DATA_FILE.replace("TEMP", var)

            with fs.open(s3_path, "rb") as f:
                ds0 = xr.open_dataset(f, decode_timedelta=False)
                ds0 = ds0.isel(**selectors)
                ds0.load()
                ds[key] = ds0[var]

    else:

        path = f"{LOCAL_PREFIX}.*.{LOCAL_SUFFIX}"

        ds0 = xr.open_mfdataset(
            path,
            decode_times=True,
            decode_timedelta=True,
        )

        ds0 = ds0.isel(**selectors)

        for key, var in VAR_MAP.items():
            ds[key] = ds0[var]

    ds["time_bound"] = ds0["time_bound"]

    return ds


# -----------------------------
# Fix time axis
# -----------------------------
def move_time_to_middle(ds):

    time_bounds = ds["time_bound"].values

    mid_times = [
        start + timedelta(days=(end - start).days / 2)
        for start, end in time_bounds
    ]

    ds = ds.assign_coords(time=("time", mid_times))

    return ds


def _mid_month(year, month):
    """Interval midpoint of a calendar month (matches move_time_to_middle)."""
    start = cftime.DatetimeNoLeap(year, month, 1)
    ey, em = (year + 1, 1) if month == 12 else (year, month + 1)
    end = cftime.DatetimeNoLeap(ey, em, 1)
    return start + timedelta(days=(end - start).days / 2)


def _ym_from_filename(path):
    """(year, month) from '...pop.h.nday1.YYYY-MM-01.nc'."""
    m = re.search(r"\.nday1\.(\d{4})-(\d{2})", os.path.basename(path))
    return int(m.group(1)), int(m.group(2))


def _monthly_mean_one_file(path):
    """Time-mean of the surface VAR_MAP fields in one daily history file.

    The 3-D fields (ALK/DIC/PO4/SiO3) are stored contiguously, so reading the
    z_t=0 hyperslab is the bottleneck (~7 s/file cold); the driver runs this
    over a process pool.
    """
    d = xr.open_dataset(path, decode_times=False, decode_timedelta=False)
    if "z_t" in d.dims:
        d = d.isel(z_t=0)
    out = {k: d[v].mean("time").values.astype("float32") for k, v in VAR_MAP.items()}
    d.close()
    y, m = _ym_from_filename(path)
    return _mid_month(y, m), out


def load_monthly_from_daily(nproc=None):
    """Monthly-mean surface fields from the daily control history files.

    Each ``.nday1.*`` file is exactly one calendar month of daily data, so the
    monthly mean is just that file's time mean - no dask resample (which is very
    slow here).  beta/eta are then computed from these monthly-mean inputs, i.e.
    f(mean), the same construction compute_truth_sensitivities.py uses for the
    truth _altco2 baseline, so the two are directly comparable.

    Returns an in-memory Dataset (time, nlat, nlon) on a mid-month axis, with
    TLONG / TLAT, ready for compute_carbonate_sensitivity / the chunk loop.
    """
    files = sorted(glob.glob(f"{LOCAL_PREFIX}.*.{LOCAL_SUFFIX}"))
    if not files:
        raise FileNotFoundError(f"no daily history files: {LOCAL_PREFIX}.*.{LOCAL_SUFFIX}")

    if nproc is None:
        nproc = min(24, len(files),
                    int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 8)))
    print(f"  {len(files)} monthly files, {_ym_from_filename(files[0])} .. "
          f"{_ym_from_filename(files[-1])};  {nproc} workers")

    mids, cols = [], {k: [] for k in VAR_MAP}
    with Pool(nproc) as pool:
        for i, (mid, out) in enumerate(pool.imap(_monthly_mean_one_file, files), 1):
            mids.append(mid)
            for k in VAR_MAP:
                cols[k].append(out[k])
            if i % 24 == 0 or i == len(files):
                print(f"    {i}/{len(files)}")

    g = xr.open_dataset(files[0], decode_times=False, decode_timedelta=False)
    mesh = {"TLONG": g["TLONG"].load(), "TLAT": g["TLAT"].load()}
    g.close()

    return xr.Dataset(
        {k: (("time", "nlat", "nlon"), np.stack(cols[k])) for k in VAR_MAP},
        coords={"time": ("time", np.array(mids)), **mesh},
    )


# -----------------------------
# Carbonate sensitivities
# -----------------------------
def compute_carbonate_sensitivity(ds):
    """beta (dDICdCO2) and eta (dDICdALK) from the surface carbonate state.

    Thin wrapper around carbonate_sensitivities.carbonate_sensitivity that maps
    this script's field names (temp/salt from VAR_MAP) onto its arguments.
    """
    return carbonate_sensitivity(
        ds.ALK, ds.DIC, ds.salt, ds.temp, ds.PO4, ds.SiO3
    )


# -----------------------------
# Run pipeline
# -----------------------------
if RESAMPLE_TO_MONTHLY:
    print("Averaging daily control-run files to monthly means...")
    ds = load_monthly_from_daily()
else:
    print("Loading data...")
    ds = load_data()
    print("Fixing time axis...")
    ds = move_time_to_middle(ds)

print("Computing carbonate sensitivities...")
if TIME_CHUNK is None:
    beta, eta = compute_carbonate_sensitivity(ds)
else:
    n_time = ds.sizes["time"]
    betas, etas = [], []
    for i in range(0, n_time, TIME_CHUNK):
        print(f"  chunk {i}–{min(i + TIME_CHUNK, n_time) - 1} / {n_time - 1}")
        ds_chunk = ds.isel(time=slice(i, i + TIME_CHUNK))
        ds_chunk.load()
        beta_chunk, eta_chunk = compute_carbonate_sensitivity(ds_chunk)
        betas.append(beta_chunk)
        etas.append(eta_chunk)
    beta = np.concatenate(betas, axis=0)
    eta = np.concatenate(etas, axis=0)

print("Building dataset...")
ds_out = xr.Dataset(
    data_vars={
        "dDICdCO2": (["time", "nlat", "nlon"], beta),
        "dDICdALK": (["time", "nlat", "nlon"], eta),
    },
    coords={
        "time": ds.time,
        "TLAT": ds.TLAT,
        "TLONG": ds.TLONG,
    },
)

OUTPUT_FILE = (
    f"/global/cfs/projectdirs/m4746/Projects/OAE-Efficiency-Map/data/"
    f"carbonate-sensitivities/carbonate_sensitivity_{file_suffix}.nc"
)

print(f"Writing {OUTPUT_FILE}")
ds_out.to_netcdf(
    OUTPUT_FILE,
    encoding={
        "time": {
            "units": "days since 0001-01-01 00:00:00",
            "calendar": "noleap",
        }
    },
)

print("Done.")
