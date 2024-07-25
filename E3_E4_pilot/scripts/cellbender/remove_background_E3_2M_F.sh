#!/bin/bash
#$ -cwd  
#$ -N cellbender_E3_2M_F
#$ -q 4-day,lg-mem  
#$ -M todd.kennedi@mayo.edu  
#$ -m abe  
#$ -pe threaded 24
#$ -l h_vmem=12G  
#$ -notify  
#$ -j y  

# source settings
source $HOME/.bash_profile

# activate conda env
source $HOME/.bash_profile
module load python
conda activate cellbender

# change directory to output folder
cd /research/labs/neurology/fryer/m214960/Ferreira_Da_Mesquita/cellbender

# run cellbender with cellranger count raw_feature_bc_matrix.h5 output from each sample
cellbender remove-background \
			--input /research/labs/neurology/fryer/m214960/Ferreira_Da_Mesquita/count/E3_2M_F/outs/raw_feature_bc_matrix.h5 \
			--output pass2_E3_2M_F_fpr_0.05.h5 \
			--expected-cells 10000 \
			--total-droplets-included 50000 \
			--fpr 0.05 \
			--epochs 150

# KEY
# --cuda: flag if using GPU
# --expected-cells: number of cells expected a priori from the experimental design
# --total-droplets-included: Choose a number that goes a few thousand barcodes into the empty droplet plateau. Include some droplets that you think are surely empty. But be aware that the larger this number, the longer the algorithm takes to run (linear).
# --fpr: false positive rate, a value of 0.01 is generally quite good, but you can generate a few output count matrices and compare them by choosing a few values: 0.01 0.05 0.1
