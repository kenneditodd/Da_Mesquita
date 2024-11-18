# load packages
setwd(".")
library(Seurat)      # Read10X_h5()
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(SeuratWrappers)

# read data
mouse.unannotated <- readRDS("../../rObjects/pass1_unannotated.rds")
Idents(mouse.unannotated) <- "seurat_clusters"
DefaultAssay(mouse.unannotated) <- "RNA"
mouse.unannotated <- NormalizeData(mouse.unannotated)

# Find markers for each cluster
markers <- SeuratWrappers::RunPrestoAll(
  object = mouse.unannotated,
  assay = "RNA",
  slot = "counts",
  only.pos = TRUE
)

# save
saveRDS(markers, "../../rObjects/pass1_unannotated_cluster_markers.rds")