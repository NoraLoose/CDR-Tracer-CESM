#!/bin/bash
#SBATCH -A m4632
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH -J dalk_ts
#SBATCH -o logs/dalk_ts.%j.out
#SBATCH -e logs/dalk_ts.%j.err

# Submit with: sbatch 11-compute-dalk-timeseries.sh
# Requires 10-compute-truth-surface.sh to have run first (reads surface.nc).
# Loops over all 690 polygons serially in a single job.
# Add --overwrite to force recompute of existing files.

MODE="oae"
SUFFIX="all-oae-daily"

module load conda
conda activate cworthy

cd "$SLURM_SUBMIT_DIR"

for PID in $(seq 0 689); do
    echo "polygon=$PID  mode=$MODE  suffix=$SUFFIX"
    python compute_dalk_timeseries.py \
        --mode "$MODE" \
        --polygons "$PID" \
        --suffix "$SUFFIX"
done
