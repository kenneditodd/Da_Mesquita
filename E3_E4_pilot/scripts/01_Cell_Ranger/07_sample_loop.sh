#!/bin/bash
for sample in E3_14M_F E3_14M_M E4_14M_F E4_14M_M
do
	sbatch 06_cellbender.sh $sample
done
