#!/usr/bin/env python
"""
Download and cache surface fields from OAE CESM-MARBL truth experiments (S3).

Reads per-month history files from the public S3 archive and writes one file
per polygon into that polygon's truth analysis subdirectory:

    {analysis_base}/{truth_subdir}/surface.nc
        ALK         [mmol m⁻³, nlat × nlon]  — surface alkalinity
        ALK_ALT_CO2 [mmol m⁻³, nlat × nlon]  — surface ALK (alt CO2)
        PH          [—,         nlat × nlon]  — surface pH
        PH_ALT_CO2  [—,         nlat × nlon]  — surface pH (alt CO2)

Coordinates:
    elapsed_months, year, month

These cached fields are the input for compute_dalk_timeseries.py (truth dALK)
and can be used directly for pH analysis.

Usage:
    python compute_truth_surface.py --polygons 0,27,140,437,651
    python compute_truth_surface.py --polygons 0 --overwrite
"""

import argparse
import os
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import h5py
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import (
    INVERSE_POLYGON_MAP,
    analysis_base as ANALYSIS_BASE,
    oae_experiment_paths,
    _get_s3fs,
)

warnings.filterwarnings("ignore", category=xr.SerializationWarning)

OUTPUT_NAME = "surface.nc"
SURFACE_VARS = ["ALK", "ALK_ALT_CO2", "PH", "PH_ALT_CO2"]
_N_WORKERS = 32


def _read_one(s3_path):
    """Read surface fields from one S3 HDF5 file.

    Uses h5py with block-caching to minimise HTTP round-trips.
    CESM POP monthly files have shape (1, nlat, nlon) for 2-D fields
    and (1, z_t, nlat, nlon) for 3-D fields; we take index [0] or [0,0].
    """
    of = _get_s3fs().open(s3_path, cache_type="blockcache", block_size=2**22)
    result = {}
    try:
        with h5py.File(of, mode="r") as hf:
            for v in SURFACE_VARS:
                raw = hf[v]
                arr = raw[0, 0] if raw.ndim == 4 else raw[0]
                result[v] = np.array(arr, dtype=np.float32)
    finally:
        of.close()
    return result


def run(pid, realization, intervention_year, intervention_month, overwrite):
    b, p    = INVERSE_POLYGON_MAP[pid]
    pid_str = f"{pid:03d}"
    label   = f"OAE {pid:03d}  ({b}, {p})"

    truth_subdir, file_glob = oae_experiment_paths(
        pid_str,
        intervention_month=intervention_month,
        intervention_year=str(intervention_year),
        realization=realization,
    )

    out_path = os.path.join(ANALYSIS_BASE, truth_subdir, OUTPUT_NAME)
    if not overwrite and os.path.exists(out_path):
        print(f"  {label}: exists, skipping")
        return

    fs    = _get_s3fs()
    files = sorted(fs.glob(file_glob))
    if not files:
        print(f"  {label}: no files found ({file_glob}), skipping")
        return
    n = len(files)
    print(f"  {label}: {n} months")

    results = {}
    completed = 0
    with ThreadPoolExecutor(max_workers=_N_WORKERS) as pool:
        future_to_idx = {pool.submit(_read_one, f"s3://{p}"): i
                         for i, p in enumerate(files)}
        for future in as_completed(future_to_idx):
            t_idx = future_to_idx[future]
            results[t_idx] = future.result()
            completed += 1
            if completed % 12 == 0:
                print(f"    {completed}/{n} months done …")

    frames = {v: [results[i][v] for i in range(n)] for v in SURFACE_VARS}

    int_m  = int(intervention_month)
    months = np.arange(n, dtype=int)
    coords = {
        "elapsed_months": ("time", months,
                           {"long_name": "elapsed months since injection start"}),
        "year":  ("time", intervention_year + (int_m - 1 + months) // 12),
        "month": ("time", (int_m - 1 + months) % 12 + 1),
    }
    _enc = {"zlib": True, "complevel": 4}

    _attrs = {
        "ALK":         {"long_name": "surface alkalinity",          "units": "mmol m-3"},
        "ALK_ALT_CO2": {"long_name": "surface alkalinity (alt CO2)","units": "mmol m-3"},
        "PH":          {"long_name": "surface pH",                  "units": "—"},
        "PH_ALT_CO2":  {"long_name": "surface pH (alt CO2)",        "units": "—"},
    }

    data_vars = {
        v: xr.DataArray(np.stack(frames[v], axis=0),
                        dims=["time", "nlat", "nlon"], attrs=_attrs[v])
        for v in SURFACE_VARS
    }

    ds_out = xr.Dataset(data_vars, coords=coords,
                        attrs=dict(mode="oae", polygon=pid, basin=b, basin_polygon=p,
                                   realization=realization,
                                   intervention_year=intervention_year,
                                   intervention_month=int_m,
                                   experiment=truth_subdir))
    os.makedirs(os.path.join(ANALYSIS_BASE, truth_subdir), exist_ok=True)
    ds_out.to_netcdf(out_path, encoding={v: _enc for v in SURFACE_VARS})
    print(f"    wrote {out_path}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--polygons", required=True,
                    help="comma-separated polygon ids, e.g. 0,27,140,437,651")
    ap.add_argument("--realization", default="001")
    ap.add_argument("--intervention-year",  type=int, default=1999)
    ap.add_argument("--intervention-month", default="01")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    polygons = sorted({int(x) for x in args.polygons.split(",") if x.strip()})
    for pid in polygons:
        run(pid, args.realization, args.intervention_year,
            args.intervention_month, args.overwrite)
    print("\nDone.")


if __name__ == "__main__":
    main()
