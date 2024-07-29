#!/bin/sh
#SBATCH --job-name n10x_count
#SBATCH --mem 75G
#SBATCH --tasks 20
#SBATCH --output logs/%x.%j.stdout
#SBATCH --error logs/%x.%j.stderr
#SBATCH --partition cpu-med
#SBATCH --time 08:00:00 ## HH:MM:SS
#SBATCH --propagate=NONE

# change directory to your desired output folder
cd /research/labs/neurology/fryer/m214960/Da_Mesquita/E3_E4_pilot/counts

# source settings
source $HOME/.bash_profile

# get cellranger version
cellranger -V

# print sample
sample=$1
echo "sample: $sample"

# run cellranger
cellranger count \
	--id=$sample \
	--sample=$sample \
	--fastqs=/research/labs/neurology/fryer/projects/Da_Mesquita/E3_E4_pilot \
	--transcriptome=/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-mm10-2020-A \
	--localcores=16 \
	--localmem=50 


