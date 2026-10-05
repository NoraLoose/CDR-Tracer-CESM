#!/bin/bash
#SBATCH -A m4632
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 4:00:00
#SBATCH -J max_dph
#SBATCH -o logs/max_dph.%j.out
#SBATCH -e logs/max_dph.%j.err

# Submit with: sbatch 13-compute-max-dph.sh
# Requires compute_truth_surface.py and compute_ph.py to have run first.
# Reads surface.nc (truth) and ph_{pid}.nc (CDR tracer) and writes
# dph_max.nc / dph_max_{pid}.nc. Pure NumPy — fast (~1 min/polygon).
# Add --overwrite to force recompute of existing files.

SUFFIX="all-oae-daily"

module load conda
conda activate cworthy

cd "$SLURM_SUBMIT_DIR"

for PID in $(seq 0 689); do
    echo "polygon=$PID  suffix=$SUFFIX"
    python compute_max_dph.py \
        --polygons "$PID" \
        --suffix "$SUFFIX"
done
