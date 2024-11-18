#!/bin/sh
#SBATCH --job-name E3PM_vs_E3CM
#SBATCH --mem 75G
#SBATCH --tasks 20
#SBATCH --output logs/%x.%j.stdout
#SBATCH --error logs/%x.%j.stderr
#SBATCH --partition cpu-short
#SBATCH --time 08:00:00 ## HH:MM:SS
#SBATCH --propagate=NONE

# source settings
source $HOME/.bash_profile

# run R script
Rscript 04a_E3PM_vs_E3CM.R
