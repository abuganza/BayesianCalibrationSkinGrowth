#!/bin/sh -l
#SBATCH -A abuganza
#SBATCH -N 1 -n 24 -t 48:00:00
#SBATCH --job-name=testjob
module load abaqus/2020
module load anaconda
unset SLURM_GTIDS
python parastudy_growth.py
