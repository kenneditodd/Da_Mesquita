#!/bin/bash
# Loop through cellranger count output folders to copy and rename html files to 
# one desired output folder

# source environment variables
source ../../refs/.env

# extract webs summaries
for sample in E3_14M_F E3_14M_M E4_14M_F E4_14M_M
do
	cd $COUNTS_DIR/$sample/outs
	cp web_summary.html ../../web_summaries/"$sample"_web_summary.html
done
