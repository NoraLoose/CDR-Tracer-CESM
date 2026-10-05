#!/usr/bin/env python
"""
Precompute surface dALK timeseries for truth and CDR-tracer (antitracer) experiments.

Reads from locally cached surface.nc files (written by compute_truth_surface.py)
for truth, and from antitracer monthly history files for the CDR tracer.

Output files:

    {analysis_base}/{truth_subdir}/dalk_timeseries.nc
        surf_true  [mmol m⁻¹]  — Σ_{i,j} (ALK − ALK_ALT_CO2)(z=0) · TAREA

    {analysis_base}/{antitracer_subdir}/dalk_timeseries_{pid}.nc
        surf_anti  [mmol m⁻¹]  — Σ_{i,j} ΔALK(z=0) · TAREA

Coordinates:
    elapsed_months, year, month

Run compute_truth_surface.py first to cache surface.nc for each polygon.

Usage:
    python compute_dalk_timeseries.py --polygons 0,27,140,437,651 --suffix all-oae-daily
    python compute_dalk_timeseries.py --polygons 0 --overwrite
"""

import argparse
import glob
import os
import sys
import warnings

import numpy as np
import pop_tools
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import (
    INVERSE_POLYGON_MAP,
    MODES,
    analysis_base as ANALYSIS_BASE,
    deficit_tracer_experiment_paths,
    open_deficit_tracer_experiment,
)

warnings.filterwarnings("ignore", category=xr.SerializationWarning)

# ── Grid ──────────────────────────────────────────────────────────────────────

_TAREA_M2 = None
_OCEAN_DA  = None


def _init_grid():
    global _TAREA_M2, _OCEAN_DA
    if _TAREA_M2 is not None:
        return
    print("Loading POP gx1v7 grid …")
    pg        = pop_tools.get_grid("POP_gx1v7")
    _TAREA_M2 = xr.DataArray(pg.TAREA.values * 1e-4, dims=["nlat", "nlon"])
    _OCEAN_DA  = xr.DataArray(pg.KMT.values > 0,      dims=["nlat", "nlon"])
    print(f"  ocean cells: {int(_OCEAN_DA.values.sum()):,}")


def _surface_integral(da):
    """Σ_{i,j} da · TAREA  [same units as da · m²]."""
    return float((da * _TAREA_M2).where(_OCEAN_DA).sum().values)


# ── Workers ───────────────────────────────────────────────────────────────────

def run_truth(mode, pid, realization, intervention_year, intervention_month, overwrite):
    b, p    = INVERSE_POLYGON_MAP[pid]
    pid_str = f"{pid:03d}"
    label   = f"{mode.upper()} {pid:03d}  ({b}, {p}) [truth]"

    truth_subdir, _ = MODES[mode](pid_str, realization=realization)
    surface_nc = os.path.join(ANALYSIS_BASE, truth_subdir, "surface.nc")
    out_path   = os.path.join(ANALYSIS_BASE, truth_subdir, "dalk_timeseries.nc")

    if not overwrite and os.path.exists(out_path):
        print(f"  {label}: exists, skipping")
        return

    if not os.path.exists(surface_nc):
        print(f"  {label}: surface.nc not found, run compute_truth_surface.py first")
        return

    ds   = xr.open_dataset(surface_nc)
    dalk = (ds["ALK"] - ds["ALK_ALT_CO2"])
    n    = len(ds["time"])

    surf_true = np.array([_surface_integral(dalk.isel(time=t)) for t in range(n)])

    coords = {k: ds.coords[k] for k in ("elapsed_months", "year", "month")}
    ds_out = xr.Dataset(
        {"surf_true": xr.DataArray(surf_true, dims=["time"],
            attrs={"long_name": "truth surface δALK area integral",
                   "units": "mmol m-1",
                   "description": "Σ_{i,j} (ALK−ALK_ALT_CO2)(z=0,i,j)·TAREA(i,j)"})},
        coords=coords,
        attrs=dict(mode=mode, polygon=pid, basin=b, basin_polygon=p,
                   realization=realization,
                   intervention_year=intervention_year,
                   intervention_month=int(intervention_month)),
    )
    ds_out.to_netcdf(out_path, encoding={"surf_true": {"zlib": True, "complevel": 1}})
    print(f"    wrote {out_path}")
    ds.close()


def run_anti(pid, suffix, realization, intervention_year, intervention_month, overwrite):
    b, p    = INVERSE_POLYGON_MAP[pid]
    pid_str = f"{pid:03d}"
    label   = f"polygon {pid:03d}  ({b}, {p}) [antitracer]"

    anti_subdir, anti_glob = deficit_tracer_experiment_paths(
        suffix,
        intervention_month=intervention_month,
        intervention_year=str(intervention_year),
        realization=realization,
    )

    out_path = os.path.join(ANALYSIS_BASE, anti_subdir, f"dalk_timeseries_{pid_str}.nc")
    if not overwrite and os.path.exists(out_path):
        print(f"  {label}: exists, skipping")
        return

    files = sorted(glob.glob(anti_glob))
    if not files:
        print(f"  {label}: no antitracer files ({anti_glob}), skipping")
        return
    n = len(files)
    print(f"  {label}: {n} months")

    alk_var    = f"DELTAALK{pid_str}"
    surf_anti  = np.full(n, np.nan)
    _noted_missing = False

    for t_idx in range(n):
        if t_idx % 12 == 0:
            yr = intervention_year + t_idx // 12
            print(f"    month {t_idx:3d}  (year {yr}) …")
        try:
            ds_a = open_deficit_tracer_experiment(
                suffix,
                intervention_month=intervention_month,
                intervention_year=str(intervention_year),
                realization=realization,
                file_index=t_idx,
            )
            if alk_var in ds_a.data_vars:
                surf_anti[t_idx] = _surface_integral(
                    ds_a[alk_var].squeeze(drop=True).isel(z_t=0))
            elif not _noted_missing:
                print(f"    Note: '{alk_var}' not in {suffix}; will be NaN.")
                _noted_missing = True
        except (FileNotFoundError, IndexError, OSError):
            pass

    int_m  = int(intervention_month)
    months = np.arange(n, dtype=int)
    coords = {
        "elapsed_months": ("time", months,
                           {"long_name": "elapsed months since injection start"}),
        "year":  ("time", intervention_year + (int_m - 1 + months) // 12),
        "month": ("time", (int_m - 1 + months) % 12 + 1),
    }
    ds_out = xr.Dataset(
        {"surf_anti": xr.DataArray(surf_anti, dims=["time"],
            attrs={"long_name": "CDR-tracer surface δALK area integral",
                   "units": "mmol m-1",
                   "description": "Σ_{i,j} ΔALK(z=0,i,j)·TAREA(i,j)"})},
        coords=coords,
        attrs=dict(polygon=pid, basin=b, basin_polygon=p, suffix=suffix,
                   realization=realization,
                   intervention_year=intervention_year,
                   intervention_month=int_m),
    )
    os.makedirs(os.path.join(ANALYSIS_BASE, anti_subdir), exist_ok=True)
    ds_out.to_netcdf(out_path, encoding={"surf_anti": {"zlib": True, "complevel": 1}})
    print(f"    wrote {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--polygons", required=True,
                    help="comma-separated polygon ids, e.g. 0,27,140,437,651")
    ap.add_argument("--mode", choices=["oae", "dor", "both"], default="oae")
    ap.add_argument("--suffix", default="all-oae-daily")
    ap.add_argument("--realization", default="001")
    ap.add_argument("--intervention-year",  type=int, default=1999)
    ap.add_argument("--intervention-month", default="01")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    polygons = sorted({int(x) for x in args.polygons.split(",") if x.strip()})
    modes    = ["oae", "dor"] if args.mode == "both" else [args.mode]

    _init_grid()

    for pid in polygons:
        for mode in modes:
            run_truth(mode, pid, args.realization,
                      args.intervention_year, args.intervention_month, args.overwrite)
        run_anti(pid, args.suffix, args.realization,
                 args.intervention_year, args.intervention_month, args.overwrite)

    print("\nDone.")


if __name__ == "__main__":
    main()
