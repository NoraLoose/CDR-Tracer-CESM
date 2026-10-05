#!/bin/bash
#SBATCH -A m4632
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH -J truth_surface
#SBATCH -o logs/truth_surface.%j.out
#SBATCH -e logs/truth_surface.%j.err

# Submit with: sbatch 10-compute-truth-surface.sh
# Downloads surface fields (ALK, ALK_ALT_CO2, PH, PH_ALT_CO2) from S3
# for all 690 polygons. Reads ~180 files per polygon from S3 — expect
# ~5-10 min/polygon. Run this before compute_dalk_timeseries.sh and
# compute_max_dph.sh.
# Add --overwrite to force redownload of existing files.

module load conda
conda activate cworthy

cd "$SLURM_SUBMIT_DIR"

for PID in $(seq 0 689); do
    echo "polygon=$PID"
    python compute_truth_surface.py \
        --polygons "$PID"
done
