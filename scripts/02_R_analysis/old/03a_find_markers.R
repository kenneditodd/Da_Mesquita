# load packages
setwd(".")
#library(Seurat, lib.loc = "/home/mayo/m214960/R/x86_64-pc-linux-gnu-library/4.2") 
library(SeuratObject, lib.loc = "/usr/local/biotools/rpackages/R-4.2.2-2023-02-01")
library(Seurat, lib.loc = "/usr/local/biotools/rpackages/R-4.2.2-2023-02-01")
options(Seurat.object.assay.version = "v4")

# read data
options(mc.cores = 10)
mouse.unannotated <- readRDS("../../rObjects/pass1_unannotated.rds")
Idents(mouse.unannotated) <- "seurat_clusters"
DefaultAssay(mouse.unannotated) <- "RNA"
mouse.unannotated <- NormalizeData(mouse.unannotated)

# Find markers for each cluster
# Compares single cluster vs all other clusters
# Default is logfc.threshold = 0.25, min.pct = 0.5
Idents(mouse.unannotated) <- "seurat_clusters"
all.markers <- FindAllMarkers(object = mouse.unannotated,
                              assay = "RNA",
                              test.use = "MAST",
                              verbose = TRUE,
                              densify = TRUE)

# add column
all.markers$delta_pct <- abs(all.markers$pct.1 - all.markers$pct.2)

# rename columns and rows
rownames(all.markers) <- 1:nrow(all.markers)
all.markers <- all.markers[,c(6,7,1,5,2:4,8)]
colnames(all.markers)[c(6,7)] <- c("pct_1","pct_2")

# save
saveRDS(all.markers, "../../rObjects/pass1_unannotated_cluster_markers.rds")