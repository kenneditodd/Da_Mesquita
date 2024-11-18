#!/bin/sh
#SBATCH --job-name presto_find_markers
#SBATCH --mem 100G
#SBATCH --tasks 10
#SBATCH --output logs/%x.%j.stdout
#SBATCH --error logs/%x.%j.stderr
#SBATCH --partition cpu-med
#SBATCH --time 24:00:00 ## HH:MM:SS
#SBATCH --propagate=NONE

# source settings
source $HOME/.bash_profile

# run script
Rscript 03a_presto_find_markers.R