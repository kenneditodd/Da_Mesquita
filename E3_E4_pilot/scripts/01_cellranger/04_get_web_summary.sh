#!/bin/bash

# extract webs summaries
for sample in E3_2M_F E3_14M_F E3_2M_M E3_14M_M E4_2M_F E4_14M_F E4_2M_M E4_14M_M
do
	cd /research/labs/neurology/fryer/m214960/Da_Mesquita/count/$sample/outs
	cp web_summary.html ../../web_summaries/"$sample"_web_summary.html
done
