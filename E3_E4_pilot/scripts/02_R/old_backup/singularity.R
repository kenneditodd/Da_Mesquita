# launch singularity
# singularity exec --nv /research/bsi/tools/biotools/monocle3-R/monocle3-leafcutter.img R

# load libraries
library(Seurat)

# read data
endothelial.annotated <- readRDS("/research/labs/neurology/fryer/m214960/Da_Mesquita/rObjects/mouse_endothelial_annotated.rds")