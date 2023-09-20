#!/bin/bash
#$ -cwd  
#$ -N n10x_count  
#$ -q 1-day  
#$ -M todd.kennedi@mayo.edu  
#$ -m abe  
#$ -pe threaded 16
#$ -l h_vmem=16G  
#$ -notify  
#$ -j y  

# change directory to your desired output folder
cd /research/labs/neurology/fryer/m214960/Da_Mesquita/count

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
	--fastqs=/research/labs/neurology/fryer/projects/Da_Mesquita/mouse/scRNA \
	--transcriptome=/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-mm10-2020-A \
	--localcores=16 \
	--localmem=50 


