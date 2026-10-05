import os
import argparse
import xarray as xr

from analysis import compute_air_sea_flux, finalize, _time_dim, deficit_tracer_experiment_paths


def main(
    suffix: str,
    mode: str,
    polygon_ids,
    intervention_month="01",
    intervention_year="1999",
):
    """
    Precompute the temporally-integrated (cumulative over the full
    experiment, i.e. 15 years) approximate air-sea CO2 flux map for
    selected polygons in a given deficit-tracer suffix. Opens the (large,
    multi-polygon) experiment dataset once and reuses it across all
    requested polygons, saving one NetCDF (2-D spatial field) per polygon.
    """
    subdir, path = deficit_tracer_experiment_paths(suffix, intervention_month, intervention_year)

    base_path = "/global/cfs/cdirs/m4746/Users/nora/Ocean-CDR-Atlas-v0/data/analysis/"

    pending = []
    for pid in polygon_ids:
        output_file = f"{base_path}/{subdir}/flux_map_{pid}.nc"
        if os.path.exists(output_file):
            print(f"Skipping polygon {pid}, file already exists: {output_file}")
            continue
        pending.append(pid)

    if not pending:
        return

    # These experiment files bundle every polygon's tracer/forcing variables
    # into each history file (thousands of variables total). Restricting each
    # per-file dataset to only the variables we need *before* xarray combines
    # them across the ~180 files avoids aligning/concatenating the thousands
    # of unused polygon variables, which otherwise dominates open_mfdataset's
    # runtime even though the actual data volume needed is tiny.
    needed_vars = {"KMT", "time_bound"}
    for pid in pending:
        needed_vars.update([f"STF_DELTADIC{pid}", f"DELTADIC{pid}_FORCING"])

    def _keep_needed(ds):
        return ds[[v for v in needed_vars if v in ds.variables]]

    ds = xr.open_mfdataset(path, decode_timedelta=False, preprocess=_keep_needed)

    for pid in pending:
        output_file = f"{base_path}/{subdir}/flux_map_{pid}.nc"

        print(f"Processing polygon {pid}")

        flux = compute_air_sea_flux(
            ds, f"STF_DELTADIC{pid}", f"DELTADIC{pid}_FORCING", cumulative=True
        )  # (total surface flux - external dic forcing) to get fco2
        if mode == "dor":
            flux = -flux  # minus sign because positive DIC forcing was applied (from alk forcing files)

        flux = finalize(flux.assign_coords(polygon_id=pid), ds)
        flux = flux.isel({_time_dim(flux): -1})

        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        flux.to_dataset(name="flux").to_netcdf(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Precompute temporally-integrated approximate air-sea "
        "CO2 flux maps for selected polygons in a deficit-tracer experiment."
    )
    parser.add_argument("suffix", type=str, help="Suffix for the approximate experiment (e.g., 'all', 'test3')")
    parser.add_argument("mode", type=str, choices=["dor", "oae"], help="Mode ('dor', 'oae')")
    parser.add_argument(
        "--polygons",
        type=str,
        nargs="+",
        required=True,
        help="Polygon IDs to compute (e.g. 000 001 653)",
    )

    args = parser.parse_args()
    main(args.suffix, args.mode, polygon_ids=args.polygons)
