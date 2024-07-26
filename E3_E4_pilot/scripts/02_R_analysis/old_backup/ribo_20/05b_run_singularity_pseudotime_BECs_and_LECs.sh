#!/bin/sh

# make sure to run on mforgehn and not mforgers
# run script in container
singularity exec --nv -B /research/labs/neurology/fryer/m214960 /research/bsi/tools/biotools/monocle3-R/monocle3-leafcutter.img Rscript 05a_pseudotime_BECs_and_LECs.R
								 
