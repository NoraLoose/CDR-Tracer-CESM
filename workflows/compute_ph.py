#!/usr/bin/env python
"""
Precompute surface pH 2D field timeseries for the CDR-tracer (antitracer) run.

Uses two file sources (control run for the antitracer base state):

  1. Control monthly means  (pop.h.YYYY-MM.nc)
       PH  → ph_ctrl  (2D surface field, read directly from MARBL)

  2. Control nday1 snapshots  (pop.h.nday1.YYYY-MM-01.nc)
       ALK, DIC, SST, SSS, PO4, SiO3  →  base carbonate state

  3. CDR-tracer antitracer monthly output  (pop.h.YYYY-MM.nc)
       DELTAALK{pid}, DELTADIC{pid}  →  perturbation on top of control

CDR pH (2D field):
    ph_cdr(t,i,j) = pH( ALK_ctrl + DELTAALK{pid},
                        DIC_ctrl + DELTADIC{pid},
                        SST_ctrl, SSS_ctrl, PO4_ctrl, SiO3_ctrl )

Output:
    {analysis_base}/{antitracer_subdir}/ph_{pid}.nc
        ph_ctrl  [nlat × nlon]  — control surface pH (MARBL)
        ph_cdr   [nlat × nlon]  — CDR-tracer surface pH (PyCO2SYS)

A separate script computes max(δpH) from these fields alongside truth surface.nc.

Usage:
    python compute_ph.py --mode oae --polygons 27
    python compute_ph.py --mode oae --polygons 0,27,437 --overwrite
"""

import argparse
import glob
import os
import warnings

import numpy as np
import pop_tools
import PyCO2SYS as pyco2
import xarray as xr

from analysis import (
    INVERSE_POLYGON_MAP,
    analysis_base as ANALYSIS_BASE,
    deficit_tracer_experiment_paths,
)

warnings.filterwarnings("ignore", category=xr.SerializationWarning)

# ---------------------------------------------------------------------------
# Control run paths
# ---------------------------------------------------------------------------
CONTROL_HIST = (
    "/global/cfs/projectdirs/m4746/Users/nora/Ocean-CDR-Atlas-v0/data/archive/"
    "smyle.cdr-atlas-v0.control.001_backup/ocn/hist/"
)
# monthly means: PH, PH_ALT_CO2, ALK_ALT_CO2, DIC_ALT_CO2
CTRL_MONTHLY_GLOB = CONTROL_HIST + "smyle.cdr-atlas-v0.control.001.pop.h.[0-9][0-9][0-9][0-9]-[0-9][0-9].nc"
# nday1 snapshots: ALK, DIC, SST, SSS, PO4, SiO3
CTRL_NDAY1_GLOB   = CONTROL_HIST + "smyle.cdr-atlas-v0.control.001.pop.h.nday1.*.nc"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SEAWATER_DENSITY = 1025.0
_TO_UMOL_KG = 1000.0 / SEAWATER_DENSITY   # mmol m⁻³ → umol kg⁻¹

# ---------------------------------------------------------------------------
# Grid (lazy init)
# ---------------------------------------------------------------------------
_TAREA_M2   = None
_OCEAN_MASK = None

def _init_grid():
    global _TAREA_M2, _OCEAN_MASK
    if _TAREA_M2 is not None:
        return
    pg = pop_tools.get_grid("POP_gx1v7")
    _TAREA_M2   = pg.TAREA.values * 1e-4   # cm² → m²
    _OCEAN_MASK = pg.KMT.values > 0

# ---------------------------------------------------------------------------
# File index helpers
# ---------------------------------------------------------------------------
def _ym_from_monthly(path):
    """(year, month) from '...pop.h.YYYY-MM.nc'."""
    stem = path.split(".pop.h.")[-1].replace(".nc", "")   # '0347-01'
    y, m = (int(x) for x in stem.split("-"))
    return y, m

def _ym_from_nday1(path):
    """(year, month) from '...pop.h.nday1.YYYY-MM-DD.nc'."""
    stem = path.split(".pop.h.nday1.")[-1].replace(".nc", "")   # '0347-01-01'
    parts = stem.split("-")
    return int(parts[0]), int(parts[1])

def _build_index(files, key_fn):
    """dict (year, month) → path."""
    return {key_fn(f): f for f in files}

def _elapsed_months(anti_files):
    y0, m0 = _ym_from_monthly(anti_files[0])
    ry = y0 if m0 > 1 else y0 - 1
    rm = m0 - 1  if m0 > 1 else 12
    return np.array(
        [(_ym_from_monthly(f)[0] - ry) * 12 + (_ym_from_monthly(f)[1] - rm)
         for f in anti_files],
        dtype=int,
    )

# ---------------------------------------------------------------------------
# pH computation
# ---------------------------------------------------------------------------
def _compute_ph_field(alk_mmol, dic_mmol, sst, sss, po4_mmol, sio3_mmol):
    """Return 2D surface pH field (nlat, nlon), masked to ocean cells."""
    csys = pyco2.sys(
        par1=alk_mmol  * _TO_UMOL_KG,
        par2=dic_mmol  * _TO_UMOL_KG,
        par1_type=1,
        par2_type=2,
        salinity=sss,
        temperature=sst,
        total_phosphate=po4_mmol  * _TO_UMOL_KG,
        total_silicate=sio3_mmol  * _TO_UMOL_KG,
    )
    ph = csys["pH"].astype(np.float32)
    ph[~_OCEAN_MASK] = np.nan
    return ph

# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------
def run(mode, pid, suffix, intervention_month="01", intervention_year="1999",
        realization="001", overwrite=False):
    _init_grid()

    pid_str = f"{pid:03d}"
    b, p = INVERSE_POLYGON_MAP[pid]
    print(f"\n=== polygon {pid_str}  ({b}, polygon {p})  mode={mode} ===")

    anti_subdir, anti_glob = deficit_tracer_experiment_paths(
        suffix,
        intervention_month=intervention_month,
        intervention_year=intervention_year,
        realization=realization,
    )

    out_path = os.path.join(ANALYSIS_BASE, anti_subdir, f"ph_{pid_str}.nc")
    if not overwrite and os.path.exists(out_path):
        print(f"  already done, skipping")
        return

    anti_files   = sorted(glob.glob(anti_glob))
    ctrl_monthly = _build_index(sorted(glob.glob(CTRL_MONTHLY_GLOB)), _ym_from_monthly)
    ctrl_nday1   = _build_index(sorted(glob.glob(CTRL_NDAY1_GLOB)),   _ym_from_nday1)

    if not anti_files:
        print(f"  [WARN] no antitracer files: {anti_glob}")
        return

    nt = len(anti_files)
    print(f"  {nt} monthly time steps")

    elapsed = _elapsed_months(anti_files)
    years   = np.array([_ym_from_monthly(f)[0] for f in anti_files], dtype=int)
    months  = np.array([_ym_from_monthly(f)[1] for f in anti_files], dtype=int)

    ph_ctrl_frames = []
    ph_cdr_frames  = []

    alk_var = f"DELTAALK{pid_str}"
    dic_var = f"DELTADIC{pid_str}"

    for i, af in enumerate(anti_files):
        ym = _ym_from_monthly(af)

        # ── control PH: read directly from MARBL monthly output ──────────────
        ctrl_ph = np.full((384, 320), np.nan, dtype=np.float32)
        if ym in ctrl_monthly:
            ds_cm = xr.open_dataset(ctrl_monthly[ym], decode_timedelta=False).squeeze(drop=True)
            ctrl_ph = ds_cm["PH"].values.astype(np.float32)
            ctrl_ph[~_OCEAN_MASK] = np.nan
            ds_cm.close()
        ph_ctrl_frames.append(ctrl_ph)

        # ── CDR pH: PyCO2SYS with ctrl state + antitracer deltas ─────────────
        cdr_ph = np.full((384, 320), np.nan, dtype=np.float32)
        if ym in ctrl_nday1:
            ds_nd = xr.open_dataset(ctrl_nday1[ym], decode_timedelta=False).squeeze(drop=True)
            alk  = ds_nd["ALK"].isel(z_t=0).values
            dic  = ds_nd["DIC"].isel(z_t=0).values
            sst  = ds_nd["SST"].values
            sss  = ds_nd["SSS"].values
            po4  = ds_nd["PO4"].isel(z_t=0).values
            sio3 = ds_nd["SiO3"].isel(z_t=0).values
            ds_nd.close()

            ds_a = xr.open_dataset(af, decode_timedelta=False).squeeze(drop=True)
            d_alk = ds_a[alk_var].isel(z_t=0).values if alk_var in ds_a else np.zeros_like(alk)
            d_dic = ds_a[dic_var].isel(z_t=0).values if dic_var in ds_a else np.zeros_like(dic)
            ds_a.close()

            cdr_ph = _compute_ph_field(alk + d_alk, dic + d_dic, sst, sss, po4, sio3)
        else:
            print(f"  [WARN] no nday1 file for {ym}, ph_cdr will be NaN for month {i}")
        ph_cdr_frames.append(cdr_ph)

        if (i + 1) % 12 == 0:
            print(f"  {i+1}/{nt}")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    coords = dict(
        elapsed_months=("time", elapsed),
        year=("time", years),
        month=("time", months),
    )
    _enc = {"zlib": True, "complevel": 4}
    ds_out = xr.Dataset(
        {
            "ph_ctrl": xr.DataArray(
                np.stack(ph_ctrl_frames, axis=0),
                dims=["time", "nlat", "nlon"],
                attrs={"long_name": "control surface pH (MARBL)", "units": "—"}),
            "ph_cdr": xr.DataArray(
                np.stack(ph_cdr_frames, axis=0),
                dims=["time", "nlat", "nlon"],
                attrs={"long_name": "CDR-tracer surface pH",
                       "units": "—",
                       "source": "PyCO2SYS(ALK_ctrl + DELTAALK, DIC_ctrl + DELTADIC)"}),
        },
        coords=coords,
        attrs={"polygon_id": pid, "basin": b, "polygon": p, "mode": mode, "suffix": suffix},
    )
    ds_out.to_netcdf(out_path, encoding={v: _enc for v in ds_out.data_vars})
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Precompute surface pH timeseries for the CDR-tracer run."
    )
    parser.add_argument("--mode", default="oae", choices=["oae", "dor"])
    parser.add_argument("--polygons", required=True,
                        help="Polygon ID or comma-separated list, e.g. 27 or 0,27,437")
    parser.add_argument("--suffix", default="all-oae-daily")
    parser.add_argument("--intervention-month", default="01")
    parser.add_argument("--intervention-year",  default="1999")
    parser.add_argument("--realization", default="001")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    for pid in [int(p) for p in args.polygons.split(",")]:
        run(
            mode=args.mode,
            pid=pid,
            suffix=args.suffix,
            intervention_month=args.intervention_month,
            intervention_year=args.intervention_year,
            realization=args.realization,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
