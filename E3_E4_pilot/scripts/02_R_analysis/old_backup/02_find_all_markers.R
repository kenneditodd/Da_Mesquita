# set working directory
setwd(".")

# libraries
require(Seurat)

# read object
mouse.unannotated <- readRDS("../../rObjects/no_integration_unannotated.rds")
Idents(mouse.unannotated) <- "seurat_clusters"
DefaultAssay(mouse.unannotated) <- "RNA"
mouse.unannotated <- NormalizeData(mouse.unannotated)

# Find markers for each cluster
# Compares single cluster vs all other clusters
# Default is logfc.threshold = 0.25, min.pct = 0.5
Idents(mouse.unannotated) <- "SCT_snn_res.0.5"
all.markers <- FindAllMarkers(object = mouse.unannotated,
                              assay = "RNA",
                              test.use = "MAST",
                              verbose = TRUE)

# add column
all.markers$delta_pct <- abs(all.markers$pct.1 - all.markers$pct.2)

# rename columns and rows
rownames(all.markers) <- 1:nrow(all.markers)
all.markers <- all.markers[,c(6,7,1,5,2:4,8)]
colnames(all.markers)[c(6,7)] <- c("pct_1","pct_2")

# save
saveRDS(all.markers, "../../rObjects/no_integration_unannotated_cluster_markers.rds")
