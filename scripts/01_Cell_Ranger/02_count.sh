#!/bin/sh
#SBATCH --job-name n10x_count       # Name of the job
#SBATCH --mem 50G                   # Amount of memory allocated for the job
#SBATCH --tasks 32                  # Number of tasks (or CPU cores)
#SBATCH --output logs/%x.%j.stdout  # File for standard output
#SBATCH --error logs/%x.%j.stderr   # File for standard error output
#SBATCH --partition cpu-med         # Specifies the partition (queue)
#SBATCH --time 08:00:00             # Maximum time the job is allowed to run, HH:MM:SS
#SBATCH --propagate=NONE            # Controls which signals are propagated to child processes

# source settings and environment variables
source $HOME/.bash_profile
source ../../refs/.env

# change directory to your desired output folder
cd $COUNTS_DIR

# get cellranger version
cellranger -V

# print sample passed from 03_sample_loop.sh script
SAMPLE=$1
echo "sample: $SAMPLE"

# run cellranger
cellranger count \
	--id=$SAMPLE \
	--sample=$SAMPLE \
	--fastqs=$FASTQ_DIR \
	--transcriptome=$ANNOTATION_REFS/mouse/refdata-gex-mm10-2020-A \
	--localcores=$SLURM_NTASKS \
	--localmem=$(($SLURM_MEM_PER_NODE / 1024))

# key:
# --id, output folder named after the sample
# --sample, sample name matching the FASTQ files
# --fastqs, directory with FASTQ files
# --transcriptome, path to reference transcriptome
# --localcores, number of CPU cores to use (32 in this case)
# --localmem, memory allocation in GB (50 GB, converting MB to GB)