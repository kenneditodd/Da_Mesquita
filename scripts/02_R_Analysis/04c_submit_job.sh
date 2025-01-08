#!/bin/sh
#SBATCH --job-name E4PM_vs_E4CM
#SBATCH --mem 50G
#SBATCH --tasks 20
#SBATCH --output logs/%x.%j.stdout
#SBATCH --error logs/%x.%j.stderr
#SBATCH --partition cpu-short
#SBATCH --time 08:00:00 ## HH:MM:SS
#SBATCH --propagate=NONE

# source settings
source $HOME/.bash_profile

# run R script
Rscript 05a_E4PM_vs_E4CM.R
