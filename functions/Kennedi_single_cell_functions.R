# Kennedi Todd
# July 24, 2024

#' cell_cycle_QC
#'
#' @param obj A Seurat object.
#' @param species A character with value mouse or human. Necessary for cell cycle phase marker list.
#' @param outDir Path to output directory where QC plots will be saved.
#' @param verbose Display messages
#' @return A character vector containing cell phases.
#' @export
#' @examples
#' #' mySeuratObject[["phase"]] <- cell_cycle_QC(mySeuratObject, "mouse")
cell_cycle_QC <- function(obj, species, markersPath, sampleCol = "sample", outDir = NULL, verbose = TRUE) {
  
  # required packages
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) {stop("obj argument is not a Seurat object")}
  if(!species %in% c("mouse","human")){stop("species argument is not a valid option")}
  if(!is(markersPath,"character")){stop("markers is not a character")}
  if(!is(sampleCol,"character")){stop("sampleCol is not a character")}
  if(!is(verbose,"logical")) { stop("verbose argument is not a logical")}
  
  # set output dir
  output <- "./"
  if(!is.null(outDir)) {output <- outDir}
  if(!endsWith(output,"/")) { output <- paste0(output,"/")}
  
  # Raw counts are not comparable between cells. Each cell has a different nCount_RNA. 
  # The log normalization, ((nCount_RNA / nFeature_RNA) * log1p, is taken in order to explore variation
  if(verbose){print("Normalizing data")}
  obj <- NormalizeData(obj)
  
  # load cell cycle markers
  phase.markers <- read.delim(markersPath, header = TRUE, sep = "\t")
  phase.markers <- subset(phase.markers, species == species)
  
  # save cycle marker list
  write.table(phase.markers, 
              paste0(output, species, "_cell_cycle_phase_markers.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # subset cell cycle markers
  g2m <- phase.markers[phase.markers$phase == "G2/M", "gene_name"]
  s <- phase.markers[phase.markers$phase == "S", "gene_name"]
  
  # score cells
  if(verbose){print("Cell cycle scoring")}
  obj <- CellCycleScoring(obj,
                          g2m.features = g2m,
                          s.features = s,
                          set.ident = TRUE)
  
  # find variable features
  if(verbose) {print("Finding variable features")}
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  
  # scale
  if(verbose){print("Scaling data")}
  obj <- ScaleData(obj)
  
  # run PCA
  if(verbose){print("Running PCA")}
  obj <- RunPCA(obj)
  
  # plot and save PCA plot
  if(verbose){print("Outputing plots")}
  pdf(paste0(output, "cell_cycle_pca.pdf"), height = 4, width = 6)
  pca <- DimPlot(obj,
                 reduction = "pca",
                 group.by = "Phase",
                 split.by = "Phase")
  print(pca)
  dev.off()
  
  # update metadata with sample column name given
  meta <- dplyr::rename(obj@meta.data, cell_cycle_sample_split=sampleCol)
  obj@meta.data <- meta
  
  # percent cells per phase plot
  pdf(paste0(output, "cell_cycle_percent_phase_per_sample.pdf"), height = 4, width = 6)
  percent.phase <- obj@meta.data %>%
    group_by(cell_cycle_sample_split, Phase) %>%
    dplyr::count() %>%
    group_by(cell_cycle_sample_split) %>%
    dplyr::mutate(percent = 100*n/sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = cell_cycle_sample_split, y = percent, fill = Phase)) +
    geom_col() +
    ggtitle("Percentage of phase per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Sample")
  print(percent.phase)
  dev.off()
  
  # return phase vector
  return(obj$Phase)
}


#' mitochondria_QC
#'
#' @param obj 
#' @param outDir 
#' @param verbose 
#' @param mitoCol
#' @return
#' @export
#'
#' @examples
mitochondria_QC <- function(obj, outDir = NULL, mitoCol = NULL, sampleCol = "sample", verbose = TRUE) {
  
  # required packages
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) {stop("obj argument is not a Seurat object.")}
  if(!is(outDir,"character")) {stop("outDir argument is not a character.")}
  if(!is(mitoCol,"character") && !is(mitoCol, "NULL")) { stop("mitoCol argument is not a character")}
  if(!is(sampleCol,"character")){stop("sampleCol is not a character")}
  if(!is(verbose,"logical")) { stop("verbose argument is not a logical")}
  
  # set output directory
  output <- "./"
  if(!is.null(outDir)) {output <- outDir}
  if(!endsWith(output,"/")) { output <- paste0(output,"/")}
  
  # Raw counts are not comparable between cells. Each cell has a different nCount_RNA. 
  # The log normalization, ((nCount_RNA / nFeature_RNA) * log1p, is taken in order to explore variation
  if(verbose){print("Normalizing data")}
  obj <- NormalizeData(obj)
  
  # find variable features
  if(verbose) {print("Finding variable features")}
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  
  # scale
  if(verbose){print("Scaling data")}
  obj <- ScaleData(obj)
  
  # run PCA
  if(verbose){print("Running PCA")}
  obj <- RunPCA(obj)
  
  # set percent.mt column name if named alternatively
  mito_column <- "percent.mt"
  if(!is.null(mitoCol)) {mito_column <- mitoCol}
  
  # set quartile values
  meta <- obj@meta.data
  first <- as.numeric(summary(meta[,mito_column])[2])
  mean <- as.numeric(summary(meta[,mito_column])[4])
  third <- as.numeric(summary(meta[,mito_column])[5])
  
  # turn percent.mt into factor based on quartile value
  if(verbose){print("Factoring percent.mt")}
  obj[["mito_factor"]] <- cut(meta[,mito_column],
                              breaks = c(-Inf, first, mean, third, Inf),
                              labels = c("Low", "Medium", "Medium high", "High"))
  
  # plot and save PCA
  if(verbose){print("Outputing plots")}
  pdf(paste0(output, "mitocondria_mito_factor_pca.pdf"), height = 4, width = 6)
  pca <- DimPlot(obj,
                 reduction = "pca",
                 group.by = "mito_factor",
                 split.by = "mito_factor")
  print(pca)
  dev.off()
  
  # update metadata with sample column name given
  meta <- dplyr::rename(obj@meta.data, mito_sample_split=sampleCol)
  obj@meta.data <- meta
  
  # percent cells per mito.factor per sample plot
  pdf(paste0(output, "mitochondria_percent_mito_factor_per_sample.pdf"), height = 4, width = 6)
  percent.mito.factor <- obj@meta.data %>%
    group_by(mito_sample_split, mito_factor) %>%
    dplyr::count() %>%
    group_by(mito_sample_split) %>%
    dplyr::mutate(percent = 100*n/sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = mito_sample_split, y = percent, fill = mito_factor)) +
    geom_col() +
    ggtitle("Percentage of mito_factor per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Sample")
  print(percent.mito.factor)
  dev.off()
  
  # return phase vector
  return(obj$mito_factor)
}

#' DE_within_each_cluster
#'
#' @param obj 
#' @param outDir 
#' @param clusterCol 
#' @param groupCol
#' @return
#' @export
#'
#' @examples
DE_within_each_cluster <- function(obj, outDir, clusterCol = "annotated_clusters", groupCol = "group", group1, group2) {
  
  # add column
  obj$clusterCol <- obj[[clusterCol]]
  obj$groupCol <- obj[[groupCol]]
  
  # initialize variables
  cell_types <- levels(obj$clusterCol)
  master.df <- data.frame()
  
  # loop through clusters
  for (i in cell_types) {
    
    # print cell type you're on
    print(i)
    
    # subset object by cell type
    cluster <- subset(obj, clusterCol == i)
    
    # set idents based on group
    Idents(cluster) <- cluster$groupCol
    
    # differential expression
    markers <- FindMarkers(object = cluster,
                           ident.1 = group1,
                           ident.2 = group2,
                           only.pos = FALSE, # default
                           min.pct = 0.10,   # default
                           test.use = "MAST",
                           verbose = TRUE,
                           assay = "RNA")
    
    # check if no DEGs
    if(nrow(markers) == 0) {
      next
    }
    
    # add markers to master table
    markers$cluster <- i
    markers$gene <- rownames(markers)
    master.df <- rbind(master.df, markers)
  }
  
  # rename columns
  master.df <- master.df %>%
    rename(
      !!paste0("percent_", group1) := "pct.1",
      !!paste0("percent_", group2) := "pct.2"
    )
  
  # set row names
  rownames(master.df) <- 1:nrow(master.df)
  
  # calculate difference
  master.df$percent_difference <- 
    abs(master.df[[paste0("percent_", group1)]] - master.df[[paste0("percent_", group2)]])
  
  
  # reorder columns
  master.df <- master.df[,c(6,7,1,5,2,3,4,8)]
  
  # check if the outDir ends with a slash
  if(!stringr::str_ends(outDir, "/")) {
    outDir <- paste0(outDir, "/")
  }
  
  # write table
  write.table(x = master.df, 
              file = paste0(outDir, group1, "_vs_", group2, "_DEGs.tsv"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
}
