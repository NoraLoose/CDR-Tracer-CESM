# CDR-Tracer-CESM

Workflow to reproduce the global CESM/MARBL CDR tracer ("antitracer") atlas from Loose et al.
(in review, JAMES): ocean alkalinity enhancement (OAE) and direct ocean removal (DOR)
efficiency curves for all 690 polygons of the
[OAE Efficiency Atlas](https://github.com/CWorthy-ocean/Ocean-CDR-Atlas-v0), computed from a
single CESM run with passive CDR tracers instead of one fully coupled CESM/MARBL simulation per
polygon.

The ROMS counterpart is [CDR-Tracer-ROMS](https://github.com/NoraLoose/CDR-Tracer-ROMS). The
notebooks that make the paper's figures are in
[CDR-Tracer-Paper](https://github.com/NoraLoose/CDR-Tracer-Paper).

This repository started as a fork of
[CWorthy-ocean/Ocean-CDR-Atlas-v0](https://github.com/CWorthy-ocean/Ocean-CDR-Atlas-v0)
(Matt Long, Anderson Banihirwe, and contributors). `atlas.py`, `cesm.py`, `machine.py` and
`config.py` derive from that repository. The truth experiments that the CDR tracer atlas is
validated against are the atlas simulations produced there. The OAE truth output is public on
[Source Cooperative](https://source.coop/cworthy/oae-efficiency-atlas), and the scripts read it
from there directly (`oae_experiment_paths` in `analysis.py`). The DOR truth output is read from
the NERSC archive `/global/cfs/projectdirs/m4746/Projects/Ocean-CDR-Atlas-v0/data/archive`
(`dor_experiment_paths`).

## Setup

```bash
conda env create -f environment.yml
conda activate cworthy
```

All Perlmutter paths (scratch, data, cases, codes) and NERSC account info live in
`workflows/config.py`. To run this yourself, edit `account` and `dir_project_root` there. Several
analysis scripts (`analysis.py`, and steps 05–09) also hard-code the NERSC paths where the
published atlas data were written. Search for `/global/cfs` to find them.

### CESM antitracer branch

The antitracer capability lives on the `anti-tracer` branch of
[CWorthy-ocean/cesm2.2.0](https://github.com/CWorthy-ocean/cesm2.2.0/tree/anti-tracer). The
paper's atlas used commit `0c1fdae`. You don't need to clone it by hand: the first time you
`import config` (which `atlas.py`/`cesm.py` do), `workflows/config.py` clones that repo into
`codes/cesm2.2.0` and checks out `anti-tracer` automatically if it isn't already there.

## Running CDR tracer experiments (`workflows/01` – `09`)

Run these in order. Notebooks are run top-to-bottom in Jupyter; `.sh` scripts are submitted with
`sbatch` from inside `workflows/`.

1. **`01-run-control.ipynb`** — Builds and submits the baseline SMYLE control case
   (`smyle.cdr-atlas-v0.control.001`) with daily POP output, via
   `atlas.global_irf_map("DOR", vintage)`. This control run is the ocean state everything else
   branches from, and its daily output feeds step 2.

   Daily output isn't on by default in stock POP2 — it comes from the SourceMods under
   `workflows/input/cesm2.2.0/SourceMods/`. `ocn.base.tavg.csh` there adds a second (`nday1`)
   tavg stream alongside the default monthly one, and `gx1v7_tavg_contents`/`ecosys_diagnostics`
   mark which fields (e.g. `SST`, `SSS`, `DIC`, `ALK`) get written to it. `cesm.create_smyle_clone`
   (called under the hood by `build()`) copies these into the case's `SourceMods/src.pop/`
   automatically based on the `curtail_output` flag, which `atlas.py` sets per case: the control
   case built here has `curtail_output=False` hard-coded (`atlas.py:643`), so it's
   **`gx1v7/src.pop/`** that actually gets used for this run — that's the directory to edit if you
   need to change which fields the control run writes daily. The antitracer cases built in step 4
   are hard-coded to `curtail_output=True` instead, so they use `curtail-output-gx1v7/src.pop/`
   (smaller ecosys output) rather than this one.

2. **`02-compute_sensitivities.sh`** — `sbatch 02-compute_sensitivities.sh`. Runs
   `compute_sensitivities.py --config daily`, which reads the control run's daily
   TEMP/SALT/ALK/DIC/PO4/SiO3 output and uses PyCO2SYS to compute the local carbonate-system
   sensitivities β = d(DIC)/d(pCO2) and η = d(DIC)/d(ALK) at daily resolution, writing
   `carbonate_sensitivity_daily.nc`. These are the β/η forcing fields that let the CDR tracer
   linearize the true nonlinear carbonate response. Requires step 1's daily output.

3. **`03-make-yamls.ipynb`** — Writes the YAML experiment configs into `antitracer-configs/`
   that describe each CDR tracer run: which basin/polygon combinations get forcing, whether it's
   a DOR (single DIC-deficit tracer) or OAE (paired ALK/DIC tracer) experiment, where the
   per-polygon alkalinity forcing files live, and where to find the β/η file from step 2. Edit
   the basin/suffix cells at the top to add new experiment sets.

4. **`04-atlas-global-antitracer-map-build-run.ipynb`** — Reads a YAML from step 3 via
   `atlas.global_irf_map("ANTITRACER", vintage, antitracer_config=<yaml path>)`, and builds a
    single CESM case from that, and submits
   (`calc.build(...)`, `calc.compute(...)`). Edit the `antitracer_config` path near the top to
   pick which YAML to run.

5. **`05-integrate.sh`** — `sbatch 05-integrate.sh <SUFFIX> <MODE>`, e.g.
   `sbatch 05-integrate.sh all-dor-daily dor`, where `SUFFIX` matches the suffix used in steps 3/4
   and `MODE` is `dor` or `oae`. Runs `integrate.py`, which spatially integrates each polygon's
   antitracer output into a single uptake/forcing time series (`integrated_<polygon>.nc` under
   `data/analysis/`). Requires step 4's case output.

6. **`06-integrate-true.sh`** — `sbatch 06-integrate-true.sh`. Runs `integrate_true.py`
   (edit the script to select `dor` or `oae`), which does the equivalent spatial integration for
   the "true" full CESM-MARBL DOR/OAE reference experiments (archived separately by the
   companion Ocean-CDR-Atlas project (https://github.com/CWorthy-ocean/Ocean-CDR-Atlas-v0),
   not built by steps 1–5) to produce the ground-truth curves the CDR tracer results are validated
   against.

7. **`07-validation-integral.ipynb`** — Loads the integrated curves from steps 5 and 6
   (`open_deficit_tracer_curve` / `open_true_curve` in `analysis.py`) and plots polygon maps,
   normalized uptake-efficiency curves, air-sea flux, and forcing, comparing the antitracer
   approximation against the true simulation. Run once 5 and 6 have output for the
   polygons/suffix you want to check.

8. **`08-integrate-flux-map.sh`** — `sbatch 08-integrate-flux-map.sh`. Runs
   `integrate_flux_map.py`, which time-integrates the CDR tracer air-sea CO2 flux anomaly into
   maps for selected polygons (`flux_map_<polygon>.nc`). These are the maps behind the paper's
   polygon flux figures (Fig. 12, S38–S40). Requires step 4's case output.

9. **`09-integrate-true-flux-map.sh`** — `sbatch 09-integrate-true-flux-map.sh`. Runs
   `integrate_true_flux_map.py`, the same for the truth experiments (`flux_map.nc` per truth
   experiment).

In steps 8 and 9, the polygon list and DOR/OAE mode are set by commenting/uncommenting lines in
the `.sh` file.

Dependency chain: `01 → 02 → 03 → 04 → 05 → 07`, with `06` running in parallel off separately
archived "true" experiments rather than off step 4's output. `08` needs `04`; `09` runs off the
archived truth experiments like `06`.

## Running your own experiment

Steps 1–9 as shipped only cover the basins/polygons/dates already baked into the config. To force
a new location, or a period outside what's already been computed, you need to do some of the
following before/inside step 3:

1. **Make a forcing file.** Each experiment in the YAML has its own `forcing_file`/`varname`
   pair, written straight into a per-tracer namelist entry (`cesm.py` ~line 458) — there's no
   requirement that every polygon/tracer point at a distinct file, or that the variable be named
   `alk_forcing`. That's just the convention `03-make-yamls.ipynb`'s `alk_forcing_file()` and
   `consolidate_forcing.py` happen to follow (one file per basin/polygon, named
   `alk-forcing-<basin>.<polygon>-<date>.nc`, variable `alk_forcing`, plus grid fields `KMT`/
   `TAREA`). Pointing several experiments at the same file/varname, or splitting one polygon's
   forcing across a few files, should work too (untested) — just set `forcing_file`/`varname`
   accordingly for each experiment entry in the YAML. `is_alk` controls how a tracer is applied
   (OAE alkalinity source vs. DOR DIC sink), independent of the file/variable naming. Generating
   forcing data for a new basin/polygon isn't done anywhere in this repo — that's produced
   upstream, in the companion OAE-Efficiency-Map forcing pipeline.

2. **Make sure your β/η (carbonate-sensitivity) file covers your run period.** The YAML's
   `beta_forcing`/`eta_forcing` blocks point at `carbonate_sensitivity_daily.nc` from step 2,
   which is only valid over the `year_first`–`year_last` window baked into it at the time it was
   computed (currently model years 346–368, i.e. however much daily control output existed when
   you last ran step 2). If your intervention period falls outside that window:
   - re-run the control case (step 1) far enough to cover the new period — extend `nyear_baseline`
     in `atlas.py`'s baseline-case block (`atlas.py:596`), or point `refdate` at a different model
     start;
   - re-run `02-compute_sensitivities.sh` against that extended daily output to regenerate
     `carbonate_sensitivity_daily.nc`;
   - update `year_first`/`year_last`/`year_align` in your YAML's `beta_forcing`/`eta_forcing`
     blocks (step 3) to match.

3. **Start dates are hard-coded.** `global_irf_map.__init__` only recognizes four start dates,
   each mapped to a fixed model ref date (`atlas.py:579–580`):
   ```python
   start_dates = ["1999-01", "1999-04", "1999-07", "1999-10"]
   ref_dates   = ["0347-01-01", "0347-04-01", "0347-07-01", "0347-10-01"]
   ```
   If your YAML's `date` isn't one of those four, step 4 raises `Invalid ANTITRACER start_date`.
   To add a new one, add a matching pair to both lists — the ref_date is whichever model calendar
   date in the control run corresponds to that real calendar date (the control run's `0347-01-01`
   origin stands in for real `1999-01-01`), and it has to fall inside the daily output range
   from step 1.

4. **New basins/polygons need their index registered in multiple places.** The `basins`/
   `npolygon` maps that assign each (basin, polygon) pair its global master index are duplicated
   in `atlas.py`, `analysis.py`, `consolidate_forcing.py`, and `03-make-yamls.ipynb` — a new
   location has to be added consistently to all of them, or the master index used for antitracer
   tracer bookkeeping and output filenames will disagree between build, integrate, and validation
   steps.

## Other files in `workflows/`

- **`atlas.py`** — Core orchestration library: the `global_irf_map` class that builds/submits CESM
  cases (control, DOR, OAE, or ANTITRACER forcing), tracks case status, and bundles/submits Slurm
  jobs. For DOR/OAE/ERW it builds one case per basin/polygon; for ANTITRACER it builds a single
  case per YAML, packing every basin/polygon experiment listed in it into that one simulation as
  separate tracers running simultaneously. Used by notebooks 01 and 04.
- **`cesm.py`** — CESM case-management helpers (cloning/configuring a case from the SMYLE
  reference case, setting up antitracer/CDR forcing namelists, submitting builds/runs). Used by
  `atlas.py`.
- **`machine.py`** — Perlmutter/Slurm helpers: polling `squeue` for running/building jobs
  (`building_jobids`, `queue_info`) and spinning up an ad hoc Dask cluster (`dask_cluster`). Used
  by `atlas.py` and `cesm.py`. Two functions in this file, `running_jobids` and `JobState`, aren't
  called anywhere else in the repo.
- **`config.py`** — Central config: project/account names, Perlmutter paths (`scratch`, `data`,
  `cases`, `codes`), and the CESM antitracer-branch clone/checkout logic described above. Imported
  by `atlas.py`, `cesm.py`, and `machine.py`.
- **`analysis.py`** — Shared analysis library: basin/polygon index maps, time-axis utilities, path
  helpers for locating true (`dor_experiment_paths`, `oae_experiment_paths`) and antitracer
  (`deficit_tracer_experiment_paths`) experiment output, and the spatial-integration/flux
  functions used by `integrate.py`, `integrate_true.py`, and notebook 07.
- **`compute_sensitivities.py`** — Script invoked by `02-compute_sensitivities.sh`; also supports
  a `--config monthly` mode (reads control output from a public S3 bucket) that isn't currently
  used by the 01–09 pipeline (which uses `--config daily`).
- **`integrate.py`** — Script invoked by `05-integrate.sh`.
- **`integrate_true.py`** — Script invoked by `06-integrate-true.sh`.
- **`consolidate_forcing.py`** / **`submit_consolidate.sh`** — One-off utility
  (`sbatch submit_consolidate.sh`) that consolidates the individual per-basin/polygon
  `alk-forcing-*.nc` files into a single `deficit_tracer_forcing_1999-01.nc` with one variable per
  polygon. Not part of the main 01–09 pipeline.
- **`interactive_job_perlmutter.sh`** — `salloc` one-liner for grabbing an interactive Perlmutter
  compute node (4 hr, `cpu` constraint, account `m4632`) to run notebooks/scripts outside batch.
- **`carbonate_sensitivities.py`** — β/η computation with PyCO2SYS, used by
  `compute_sensitivities.py`.
- **`atlas_cesm_overview.md`** — Short overview of how `atlas.py` and `cesm.py` divide the work.
- **`antitracer-configs/`** — YAML experiment configs, mostly generated by `03-make-yamls.ipynb`
  (`*-all-dor-daily`, `*-all-oae-daily`, `*-test-dor(-daily)`, `*-test-oae(-daily)`). The paper
  uses `all-dor-daily` and `all-oae-daily`.
- **`input/cesm2.2.0/`** — SourceMods and reference-case files copied into each CESM case.
- **`logs/`** — `sbatch`/Slurm stdout and stderr logs from the scripts above.

## Known issues

### Manual `./case.submit` resubmit can SIGSEGV (NetCDF library mismatch)

If a running case crashes and you resubmit it by hand from its case directory
(`cd cases/<case>; ./case.submit`), it can die ~20s into the restart with:

```
rank 1: forrtl: severe (174): SIGSEGV, segmentation fault
  libnetcdff_intel.so → nf_mod.F90:903 (nf90_inq_varid)
                      → seq_io_read_mod.F90:272
                      → seq_infodata_mod.F90:877
```

right after opening the coupler restart file. **Cause:** Perlmutter CPE upgrades (e.g. 23.12→25.09)
change what `CRAY_LD_LIBRARY_PATH` prepends at `srun` time, so the wrong (newer) `libnetcdff_intel.so`
gets loaded instead of the Intel 2023.2 NetCDF libs the executable was actually built against. The
bundle-submit path that `atlas.py` generates for a fresh build/run (`submit_bundle`/`submit_build`,
`atlas.py` ~lines 100–138) already exports the correct `LD_LIBRARY_PATH` and unsets
`CRAY_LD_LIBRARY_PATH` before `srun`, so cases submitted through steps 1/4 don't hit this. But
CIME's auto-generated `.case.run` — used for every RESUBMIT and for a manual `./case.submit` — did
not, until it was fixed at the source.

**Fix:** the same library pin was added directly to
`codes/cesm2.2.0/cime/config/cesm/machines/template.case.run` (commit `8df28ca`, "Pin Intel 2023.2
NetCDF/compiler libraries before any srun call", on the `anti-tracer` branch). Every case's
`.case.run` is regenerated from this template, so the fix applies automatically to future
auto-resubmits and manual submits alike — **except** CIME only regenerates `.case.run` when
`env_batch.xml`'s timestamp differs from its locked copy in `LockedFiles/env_batch.xml`. If you hit
this crash on an existing case whose `.case.run` predates the template fix, `touch env_batch.xml`
in the case directory before running `./case.submit` to force regeneration.
