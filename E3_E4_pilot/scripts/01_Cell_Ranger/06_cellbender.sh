#!/bin/sh
#SBATCH --job-name cellbender       # Name of the job
#SBATCH --mem 50G                   # Amount of memory allocated for the job
#SBATCH --output logs/%x.%j.stdout  # File for standard output
#SBATCH --error logs/%x.%j.stderr   # File for standard error output
#SBATCH --partition gpu             # Specifies the partition (queue)
#SBATCH --time 08:00:00             # Maximum time the job is allowed to run, HH:MM:SS
#SBATCH --gres=gpu:1                # Request 1 GPU

# change to top level of project directory
cd ../../

# source settings and environment variables
source $HOME/.bash_profile
source refs/.env

# cellbender version
cellbender --version

# print sample passed from 03_sample_loop.sh script
SAMPLE=$1
echo "sample: $SAMPLE"

# check cuda
python -c "import torch; print(torch.cuda.is_available())"

# run cellbender remove-background
cellbender remove-background \
      --cuda \
      --input counts/$SAMPLE/outs/raw_feature_bc_matrix.h5 \
      --output cellbender/$SAMPLE_cellbender.h5 \
      --expected-cells 5000
      
# Key
# --input un-filtered data file, .h5 or matrix directory from 10x
# --output output file location containing .h5 extension
# --cude will run on reference GPU