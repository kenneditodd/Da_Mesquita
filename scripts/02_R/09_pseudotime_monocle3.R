# Kennedi Todd
# September 12, 2022

# Set working directory
setwd(".")

# Load libraries
library(ggplot2)
library(monocle3)
library(Seurat)
library(SeuratWrappers)

# Set seed
set.seed(8)

# Read preprocessed SeuratObject
endothelial.annotated <- readRDS("mouse_endothelial_annotated.rds")
DefaultAssay(endothelial.annotated) <- "RNA"
Idents(endothelial.annotated) <- "individual_clusters"

# Plot UMAP
DimPlot(endothelial.annotated,
        cols = c("firebrick1","gold","green","cyan","blue","darkorchid1", "lightgray"))

# Convert SeuratObject to cell_data_set
cds <- as.cell_data_set(endothelial.annotated, default.reduction = "UMAP")

# Preview cell metadata
head(colData(cds))

# Set gene metadata
head(fData(cds)) # empty
head(rownames(fData(cds))) # short gene names
fData(cds)$gene_short_name <- rownames(fData(cds)) # assign
head(fData(cds)) # check

# Preview counts
head(counts(cds))

# Transfer cluster info
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters$UMAP$partitions <- recreate.partition
cds@clusters$UMAP$clusters <- endothelial.annotated$individual_clusters
cds@int_colData@listData$reducedDims$UMAP <-endothelial.annotated@reductions$umap@cell.embeddings

# Plot UMAP of cell_data_set object
p1 <- plot_cells(cds,
                 color_cells_by = "individual_clusters",
                 label_groups_by_cluster = FALSE,
                 group_label_size = 5) +
  scale_color_manual(values = c("firebrick1","gold","green","cyan","blue","darkorchid1", "lightgray")) +
  theme(legend.position = "right")
p1

# Learn trajectory
cds <- learn_graph(cds, use_partition = TRUE)

# Plot trajectory
p2 <- plot_cells(cds,
                 color_cells_by = "individual_clusters",
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 group_label_size = 5) +
  scale_color_manual(values = c("firebrick1","gold","green","cyan","blue","darkorchid1", "lightgray")) +
  theme(legend.position = "right")
p2

# Order cells in pseudotime
# Must have prior knowledge of root cells - will use LECs
starting.cells <- cds@colData$individual_clusters == "LECs"
cds <- order_cells(cds, reduction_method = "UMAP", 
                   root_cells = colnames(cds[,clusters(cds) == "LECs"]))

# Plot pseudotime
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 group_label_size = 5)
p3
# split pseudotime plot
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 group_label_size = 5) + facet_wrap(~sex)
p3



