#!/bin/bash
#SBATCH -A m4632
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 24:00:00
#SBATCH -J ph_ts
#SBATCH -o logs/ph_ts.%j.out
#SBATCH -e logs/ph_ts.%j.err

# Submit with: sbatch 12-compute-ph.sh
# Can run in parallel with 11-compute-dalk-timeseries.sh (independent inputs).
# PyCO2SYS is called once per surface grid cell per month — expect ~1-3 min/polygon.
# Add --overwrite to force recompute of existing files.

MODE="oae"
SUFFIX="all-oae-daily"

module load conda
conda activate cworthy

cd "$SLURM_SUBMIT_DIR"

for PID in $(seq 0 689); do
    echo "polygon=$PID  mode=$MODE  suffix=$SUFFIX"
    python compute_ph.py \
        --mode "$MODE" \
        --polygons "$PID" \
        --suffix "$SUFFIX"
done
