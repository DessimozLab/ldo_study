#!/usr/bin/env bash
#$ -l tmem=15G
#$ -l h_vmem=15G
#$-l h_rt=12:00:0
#$ -R y
#$ -cwd -V
#$ -notify
#$ -j y -o /SAN/ugi/SSSoc/alex-2019/family_scaling/log
#$ -S /bin/bash
#$ -N pairwise-exprcalc
#SBATCH --job-name=exprcalc
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

echo "Running on `hostname`"
echo "Start `date`"
# note: call this for each species with e.g., ANOCA
python ./scripts/expression_calc.py $1 ./results/expression ./results/pairwise/$1.tsv ./results/pairwise/$1.res.tsv

echo "Finished `date`"
