#!/bin/bash
#SBATCH -A m4632
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH -J integrate_true_flux_map
#SBATCH -o logs/integrate_true_flux_map.%j.out
#SBATCH -e logs/integrate_true_flux_map.%j.err
#SBATCH --cpus-per-task=1

# Load conda module
module load conda
conda activate cworthy

#srun -n 1 python integrate_true_flux_map.py dor --polygons 000 001 002 150 151 152 350 351 352 651 652 653
#srun -n 1 python integrate_true_flux_map.py oae --polygons 000 001 002 150 151 152 350 351 352 651 652 653
srun -n 1 python integrate_true_flux_map.py dor --polygons 437 376 027 091 674 653 164 678 011 140 060 127
#srun -n 1 python integrate_true_flux_map.py oae --polygons 437 376 027 091 674 653 164 678 011 140 060 127

