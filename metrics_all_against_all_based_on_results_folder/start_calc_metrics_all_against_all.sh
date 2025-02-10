#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=1-0
#SBATCH --error=/nfs/home/students/l.spindler/systems_biomed/slurm/err/%j.err
#SBATCH --output=/nfs/home/students/l.spindler/systems_biomed/slurm/out/%j.out
#SBATCH --job-name=metrics_all

ml load r/4.4.2

srun Rscript /nfs/home/students/l.spindler/calc_metrics_all_against_all.R