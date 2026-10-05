#!/usr/bin/env python
"""
Compute max(δpH) maps for truth and CDR-tracer experiments.

For each polygon, takes the per-cell temporal maximum of the pH perturbation:

    Truth:      δpH(t,i,j) = PH(t,i,j) − PH_ALT_CO2(t,i,j)
    CDR tracer: δpH(t,i,j) = ph_cdr(t,i,j) − ph_ctrl(t,i,j)

Reads from:
    {truth_subdir}/surface.nc          (written by compute_truth_surface.py)
    {antitracer_subdir}/ph_{pid}.nc    (written by compute_ph.py)

Writes:
    {truth_subdir}/dph_max.nc
        dph_max        [nlat × nlon]  — max δpH over all timesteps
        dph_max_month  [nlat × nlon]  — elapsed month when max occurs

    {antitracer_subdir}/dph_max_{pid}.nc
        dph_max        [nlat × nlon]  — max δpH over all timesteps
        dph_max_month  [nlat × nlon]  — elapsed month when max occurs

Usage:
    python compute_max_dph.py --polygons 0,27,140,437,651 --suffix all-oae-daily
    python compute_max_dph.py --polygons 0 --overwrite
"""

import argparse
import os
import sys
import warnings

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import (
    INVERSE_POLYGON_MAP,
    MODES,
    analysis_base as ANALYSIS_BASE,
    deficit_tracer_experiment_paths,
)

warnings.filterwarnings("ignore", category=xr.SerializationWarning)

_ENC = {"zlib": True, "complevel": 4, "dtype": "float32"}


def _max_dph(dph_arr, elapsed_months):
    """
    Given dph_arr (time, nlat, nlon), return (max_val, month_of_max) both (nlat, nlon).
    NaN where no valid data.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        idx = np.nanargmax(dph_arr, axis=0)          # (nlat, nlon)
        val = np.nanmax(dph_arr,    axis=0)           # (nlat, nlon)
    month_of_max = elapsed_months[idx].astype(np.float32)
    all_nan = np.all(~np.isfinite(dph_arr), axis=0)
    val[all_nan]          = np.nan
    month_of_max[all_nan] = np.nan
    return val.astype(np.float32), month_of_max


def _write(out_path, dph_max, dph_max_month, attrs):
    ds = xr.Dataset(
        {
            "dph_max": xr.DataArray(dph_max, dims=["nlat", "nlon"],
                attrs={"long_name": "maximum surface δpH", "units": "—"}),
            "dph_max_month": xr.DataArray(dph_max_month, dims=["nlat", "nlon"],
                attrs={"long_name": "elapsed month of maximum δpH", "units": "months"}),
        },
        attrs=attrs,
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    ds.to_netcdf(out_path, encoding={v: _ENC for v in ds.data_vars})
    print(f"    wrote {out_path}")


def run(pid, suffix, realization, intervention_year, intervention_month, overwrite):
    b, p    = INVERSE_POLYGON_MAP[pid]
    pid_str = f"{pid:03d}"
    label   = f"OAE {pid:03d}  ({b}, {p})"

    truth_subdir, _ = MODES["oae"](pid_str, realization=realization)
    anti_subdir, _  = deficit_tracer_experiment_paths(
        suffix,
        intervention_month=intervention_month,
        intervention_year=str(intervention_year),
        realization=realization,
    )

    out_truth = os.path.join(ANALYSIS_BASE, truth_subdir, "dph_max.nc")
    out_anti  = os.path.join(ANALYSIS_BASE, anti_subdir,  f"dph_max_{pid_str}.nc")

    truth_done = not overwrite and os.path.exists(out_truth)
    anti_done  = not overwrite and os.path.exists(out_anti)
    if truth_done and anti_done:
        print(f"  {label}: outputs exist, skipping")
        return

    _base = dict(polygon=pid, basin=b, basin_polygon=p,
                 realization=realization,
                 intervention_year=intervention_year,
                 intervention_month=int(intervention_month))

    # ── Truth ──────────────────────────────────────────────────────────────────
    if not truth_done:
        surface_nc = os.path.join(ANALYSIS_BASE, truth_subdir, "surface.nc")
        if not os.path.exists(surface_nc):
            print(f"  {label}: surface.nc missing, run compute_truth_surface.py first")
        else:
            ds = xr.open_dataset(surface_nc)
            dph = (ds["PH"] - ds["PH_ALT_CO2"]).values        # (time, nlat, nlon)
            elapsed = ds["elapsed_months"].values
            ds.close()
            val, mon = _max_dph(dph, elapsed)
            _write(out_truth, val, mon, {**_base, "source": "surface.nc"})

    # ── CDR tracer ─────────────────────────────────────────────────────────────
    if not anti_done:
        ph_nc = os.path.join(ANALYSIS_BASE, anti_subdir, f"ph_{pid_str}.nc")
        if not os.path.exists(ph_nc):
            print(f"  {label}: ph_{pid_str}.nc missing, run compute_ph.py first")
        else:
            ds = xr.open_dataset(ph_nc)
            dph = (ds["ph_cdr"] - ds["ph_ctrl"]).values        # (time, nlat, nlon)
            elapsed = ds["elapsed_months"].values
            ds.close()
            val, mon = _max_dph(dph, elapsed)
            _write(out_anti, val, mon, {**_base, "suffix": suffix,
                                        "source": f"ph_{pid_str}.nc"})


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--polygons", required=True,
                    help="comma-separated polygon ids, e.g. 0,27,140,437,651")
    ap.add_argument("--suffix", default="all-oae-daily")
    ap.add_argument("--realization", default="001")
    ap.add_argument("--intervention-year",  type=int, default=1999)
    ap.add_argument("--intervention-month", default="01")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    polygons = sorted({int(x) for x in args.polygons.split(",") if x.strip()})
    for pid in polygons:
        run(pid, args.suffix, args.realization,
            args.intervention_year, args.intervention_month, args.overwrite)
    print("\nDone.")


if __name__ == "__main__":
    main()
