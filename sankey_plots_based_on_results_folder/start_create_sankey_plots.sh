#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=70G
#SBATCH --time=1-0
#SBATCH --error=/nfs/home/students/l.spindler/systems_biomed/slurm/err/%j.err
#SBATCH --output=/nfs/home/students/l.spindler/systems_biomed/slurm/out/%j.out
#SBATCH --job-name=sankey

ml load r/4.4.2

srun Rscript /nfs/home/students/l.spindler/create_sankey_plots.R