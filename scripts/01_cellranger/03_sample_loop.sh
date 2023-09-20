#!/bin/bash
for sample in E3_2M_F E3_14M_F E3_2M_M E3_14M_M E4_2M_F E4_14M_F E4_2M_M E4_14M_M
do
	qsub 02_count.sh $sample
done

