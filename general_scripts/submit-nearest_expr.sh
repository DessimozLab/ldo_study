#!/usr/bin/env bash
#SBATCH --job-name=expr
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

echo "Running on `hostname`"
echo "Start `date`"
# note: call this for each species with e.g., ANOCA.npz
python nearest_expr.py raw/$1 nearest_neighbour/$1

echo "Finished `date`"
