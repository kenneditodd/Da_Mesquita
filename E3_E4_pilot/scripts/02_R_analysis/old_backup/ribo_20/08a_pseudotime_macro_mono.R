# Kennedi Todd
# September 12, 2022

# Set working directory
setwd("/research/labs/neurology/fryer/m214960/Da_Mesquita/scripts/R")

# Load libraries
library(ggplot2)
library(monocle3)
library(Seurat)
library(SeuratWrappers)

# Set seed
set.seed(8)

# save function
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

# checkpoint
print("Setup ready")

# Read preprocessed SeuratObject
macro.annotated <- readRDS("../../rObjects/mouse_macro_annotated.rds")
DefaultAssay(macro.annotated) <- "RNA"
Idents(macro.annotated) <- "merged_clusters"

# Add new labels
sex.isoform <- paste0(macro.annotated$sex," ",macro.annotated$isoform)
macro.annotated$sex_isoform <- sex.isoform

# Plot UMAP
path <- "../../results/ribo_20/reclusterMacro/pseudotime/UMAP_merged_clusters.pdf"
pdf(path, width = 8, height = 6)
DimPlot(macro.annotated,
        cols = c("firebrick1","gold","green","forestgreen","cyan","blue","gray40"))
dev.off()
print("Normal umap saved")

# Convert SeuratObject to cell_data_set
cds <- as.cell_data_set(macro.annotated, default.reduction = "UMAP")

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
cds@clusters$UMAP$clusters <- macro.annotated$merged_clusters
cds@int_colData@listData$reducedDims$UMAP <-macro.annotated@reductions$umap@cell.embeddings

# Plot UMAP of cell_data_set object
plot_cells(cds,
           color_cells_by = "merged_clusters",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c("firebrick1","gold","green","cyan","blue","darkorchid1", "lightgray")) +
  theme(legend.position = "right")
dev.off()

# Learn trajectory
cds <- learn_graph(cds, use_partition = TRUE)

# Plot trajectory
path <- "../../results/ribo_20/reclusterMacro/pseudotime/UMAP_merged_cluster_trajectory.pdf"
pdf(path, height = 6, width = 8)
plot_cells(cds,
           color_cells_by = "merged_clusters",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c("firebrick1","gold","green","cyan","blue","darkorchid1", "lightgray")) +
  theme(legend.position = "right")
dev.off()
print("Trajectory saved")

# Order cells in pseudotime
# Must have prior knowledge of root cells - will use LECs
starting.cells <- cds@colData$merged_clusters == "Monocytes"
cds <- order_cells(cds, reduction_method = "UMAP", 
                   root_cells = colnames(cds[,clusters(cds) == "Monocytes"]))

# Plot pseudotime
path <- "../../results/ribo_20/reclusterMacro/pseudotime/UMAP_merged_cluster_pseudotime.pdf"
pdf(path, height = 6, width = 8)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)
dev.off()
# split by sex
path <- "../../results/ribo_20/reclusterMacro/pseudotime/UMAP_individual_cluster_pseudotime_split_sex.pdf"
pdf(path, height = 6, width = 8)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5) + facet_wrap(~sex)
dev.off()
# split by age
path <- "../../results/ribo_20/reclusterMacro/pseudotime/UMAP_individual_cluster_pseudotime_split_age.pdf"
pdf(path, height = 6, width = 8)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5) + facet_wrap(~age)
dev.off()
# split by isoform
path <- "../../results/ribo_20/reclusterMacro/pseudotime/UMAP_individual_cluster_pseudotime_split_isoform.pdf"
pdf(path, height = 6, width = 8)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5) + facet_wrap(~isoform)
dev.off()
# split by sex and isoform
path <- "../../results/ribo_20/reclusterMacro/pseudotime/UMAP_individual_cluster_pseudotime_split_sex&isoform.pdf"
pdf(path, height = 6, width = 8)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5) + facet_wrap(~sex_isoform)
dev.off()

# checkpoint
print("Pseudotime saved")

