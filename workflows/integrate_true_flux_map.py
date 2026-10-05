import os
import argparse
import xarray as xr

from analysis import compute_air_sea_flux, finalize, _time_dim, dor_experiment_paths, oae_experiment_paths

MODES = {
    "dor": dor_experiment_paths,
    "oae": oae_experiment_paths,
}


def main(
    mode,
    polygon_ids,
    intervention_month="01",
    intervention_year="1999",
):
    """
    Precompute the temporally-integrated (cumulative over the full
    experiment, i.e. 15 years) air-sea CO2 flux map for selected polygons
    in the true CESM-MARBL experiments. Saves one NetCDF (2-D spatial
    field) per polygon, so plots like `compare_air_sea_flux` can load the
    result instantly instead of re-reading and cumsum-ing 180 monthly
    files on every call.
    """
    analysis_base = (
        "/global/cfs/projectdirs/m4746/Users/nora/"
        "Ocean-CDR-Atlas-v0/data/analysis"
    )

    experiment_paths = MODES[mode]

    for pid in polygon_ids:
        subdir, path = experiment_paths(pid, intervention_month, intervention_year)
        output_file = f"{analysis_base}/{subdir}/flux_map.nc"

        if os.path.exists(output_file):
            print(f"Skipping polygon {pid}, file exists: {output_file}")
            continue

        print(f"Processing polygon {pid}")

        ds = xr.open_mfdataset(path, decode_timedelta=False)

        flux = compute_air_sea_flux(ds, "FG_CO2", "FG_ALT_CO2", cumulative=True)
        flux = finalize(flux, ds)
        flux = flux.isel({_time_dim(flux): -1})

        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        flux.to_dataset(name="flux").to_netcdf(output_file)

        ds.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Precompute temporally-integrated air-sea CO2 flux maps "
        "for selected polygons in the true CESM-MARBL experiments."
    )
    parser.add_argument(
        "mode",
        choices=["dor", "oae"],
        help="Experiment type: 'dor' or 'oae'",
    )
    parser.add_argument(
        "--polygons",
        type=str,
        nargs="+",
        required=True,
        help="Polygon IDs to compute (e.g. 000 001 653)",
    )

    args = parser.parse_args()

    main(mode=args.mode, polygon_ids=args.polygons)
