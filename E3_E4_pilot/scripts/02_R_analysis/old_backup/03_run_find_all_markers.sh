#!/bin/sh
#$ -cwd
#$ -N find_all_markers_no_integration
#$ -m abe
#$ -M todd.kennedi@mayo.edu
#$ -l h_vmem=30G
#$ -q 1-day,4-day,lg-mem
#$ -notify

# load environment
source $HOME/.bash_profile
source ~/.bash_mayobiomods7

# change directory
CWD="/research/labs/neurology/fryer/m214960/Da_Mesquita/scripts/02_R"
cd $CWD

# run script
Rscript  02_find_all_markers.R