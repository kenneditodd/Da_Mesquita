# Kennedi Todd
# Jan 2, 2024

# set wd
setwd(".")

# load packages
library(ComplexUpset) 
library(dplyr)        # ungroup()
library(ggrepel)      # geom_text_repel()
library(ggtree)       # BuildClusterTree()
library(gridExtra)    # grid.arrange()
library(gtools)       # smartbind()
library(parallel)     # detectCores()
library(plotly)       # plot_ly()
library(Seurat)       # DimPlot()
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(tidyr)        # %>%
library(UpSetR)       # fromList()

# work in parallel
options(mc.cores = 20)

# load object
mouse.annotated <- readRDS("../../rObjects/pass1_annotated.rds")

############################### DE: E4PM vs E4CM ###############################

# print comparison
print("E4PM vs E4CM")

# intitialize variables
cell_types <- levels(mouse.annotated$annotated_clusters)
e3.male.df <- data.frame()

# loop through clusters
for (i in cell_types) {
  print(i)
  cluster <- subset(mouse.annotated, annotated_clusters == i)
  Idents(cluster) <- cluster$group
  markers <- FindMarkers(object = cluster,
                         ident.1 = "E4PM",
                         ident.2 = "E4CM",
                         only.pos = FALSE, # default
                         min.pct = 0.10,   # default
                         test.use = "MAST",
                         verbose = TRUE,
                         assay = "RNA")
  if(nrow(markers) == 0) {
    next
  }
  markers$cluster <- i
  markers$gene <- rownames(markers)
  e3.male.df <- rbind(e3.male.df, markers)
}

# reformat table
colnames(e3.male.df)[c(3,4)] <- c("percent_E4PM","percent_E4CM")
rownames(e3.male.df) <- 1:nrow(e3.male.df)
e3.male.df$percent_difference <- abs(e3.male.df$percent_E4PM - e3.male.df$percent_E4CM)
e3.male.df <- e3.male.df[,c(6,7,1,5,2,3,4,8)]

# write table
write.table(e3.male.df, "../../results/all_clusters_pass1/DEGs/DEG_tables/E4PM_vs_E4CM_DEGs.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
