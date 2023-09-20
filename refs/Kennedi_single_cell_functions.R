# Kennedi Todd
# April 19, 2023

#' find_min_pc
#'
#' @param obj A Seurat object
#' @param reduct By default find_min_pc() will use harmony dimensional reduction if present. Otherwise, pca reduction will be used. Use this argument to force and specify which reduction you wish to use.
#' @param plot By default find_min_pc() will not plot the minimum principal component cutoff. Set argument equal to TRUE if you want plotted.
#' @return An integer with the minimum PC you should use for UMAP and clustering.
#' @export
#' @examples
#' min.pc <- find_min_pc(mySeuratObject)
find_min_pc <- function(obj, reduct = NULL, plot = FALSE){
  
  # libraries
  #require(ggplot2)
  
  # initialize variables
  reduction <- "pca"
  
  # check if input is a seurat object and if they have a pca/harmony reduction
  if(!is(obj,"Seurat")){
    stop("Input is not a Seurat object.")
  }
  if(!"pca" %in% names(obj@reductions)){
    stop("Input does not have pca reduction.")
  }
  if("harmony" %in% names(obj@reductions)){
    print("Using harmony dimensional reduction over pca.")
    reduction <- "harmony"
  }
  if(!is.null(reduct)) {
    reduction <- reduct
  }
  
  # Determine percent of variation associated with each PC
  stdv <- obj[[reduction]]@stdev
  sum.stdv <- sum(obj[[reduction]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  
  # Calculate cumulative percents for each PC
  cumulative <- cumsum(percent.stdv)
  
  # Determine which PC exhibits cumulative percent greater than 90% and
  # and % variation associated with the PC as less than 5
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  
  # Determine the difference between variation of PC and subsequent PC
  # last point where change of % of variation is more than 0.1%
  co2 <- sort(which(
    (percent.stdv[1:length(percent.stdv) - 1] - 
       percent.stdv[2:length(percent.stdv)]) > 0.1), 
    decreasing = T)[1] + 1
  
  # min.pc
  min.pc <- min(co1,co2)
  
  if(plot){
    # Create a dataframe with values
    plot_df <- data.frame(pct = percent.stdv, 
                          cumu = cumulative, 
                          rank = 1:length(percent.stdv))
    
    # Elbow plot to visualize 
    p <- ggplot(plot_df, aes(cumulative, percent.stdv, label = rank, color = rank > min.pc)) + 
      geom_text() + 
      geom_vline(xintercept = 90, color = "grey") + 
      geom_hline(yintercept = min(percent.stdv[percent.stdv > 5]), color = "grey") +
      theme_bw()
    print(p)
  }
  
  # min.pc
  return(min.pc)
}


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
cell_cycle_QC <- function(obj, species, outDir = NULL, verbose = TRUE) {
  
  # required packages
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) {stop("obj argument is not a Seurat object")}
  if(!species %in% c("mouse","human")){stop("species argument is not a valid option")}
  if(!is(verbose,"logical")) { stop("verbose argument is not a logical")}
  
  # set output directory
  output <- "./"
  if(!is.null(outDir)) {output <- outDir}
  if(!endsWith(output,"/")) { output <- paste0(output,"/")}
  
  # Raw counts are not comparable between cells. Each cell has a different nCount_RNA. 
  # The log normalization, ((nCount_RNA / nFeature_RNA) * log1p, is taken in order to explore variation
  if(verbose){print("Normalizing data")}
  obj <- NormalizeData(obj)
  
  # load cell cycle markers
  path <- "/research/labs/neurology/fryer/m214960/Da_Mesquita/refs/cell_cycle_markers.tsv"
  phase.markers <- read.delim(path, header = TRUE, sep = "\t")
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
  
  # percent cells per phase plot
  pdf(paste0(output, "cell_cycle_percent_phase_per_sample.pdf"), height = 4, width = 6)
  percent.phase <- obj@meta.data %>%
    group_by(sample, Phase) %>%
    dplyr::count() %>%
    group_by(sample) %>%
    dplyr::mutate(percent = 100*n/sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = sample, y = percent, fill = Phase)) +
    geom_col() +
    ggtitle("Percentage of phase per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
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
mitochondria_QC <- function(obj, outDir = NULL, mitoCol = NULL, verbose = TRUE) {
  
  # required packages
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) {stop("obj argument is not a Seurat object.")}
  if(!is(outDir,"character")) {stop("outDir argument is not a character.")}
  if(!is(mitoCol,"character") && !is(mitoCol, "NULL")) { stop("mitoCol argument is not a character")}
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
  
  # percent cells per mito.factor per sample plot
  pdf(paste0(output, "mitochondria_percent_mito_factor_per_sample.pdf"), height = 4, width = 6)
  percent.mito.factor <- obj@meta.data %>%
    group_by(sample, mito_factor) %>%
    dplyr::count() %>%
    group_by(sample) %>%
    dplyr::mutate(percent = 100*n/sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = sample, y = percent, fill = mito_factor)) +
    geom_col() +
    ggtitle("Percentage of mito_factor per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(percent.mito.factor)
  dev.off()
  
  # return phase vector
  return(obj$mito_factor)
}


#' subset_recluster
#'
#' @param obj 
#' @param cellIdent 
#' @param cellType 
#' @param projectID
#' @param plot logical indicating wether or not to display plots
#' @param integrate logical indicating wether integration should be performed
#' @return
#' @export
#'
#' @examples
subset_recluster <- function(obj, cellIdent, cellType, projectID, keepCol, plot = TRUE, integrate = FALSE, verbose = TRUE) {
  
  # load packages
  require(scCustomize)
  require(sctransform)
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) { stop("obj argument is not a Seurat object")}
  if(!is(cellIdent,"character")) { stop("cellIdent argument is not a character")}
  if(!is(cellType,"character")) { stop("cellType argument is not a character")}
  if(!is(projectID,"character")) { stop("projectID argument is not a character")}
  if(!is(keepCol,"numeric")) { stop("keepCol argument is not a numeric")}
  if(!is(plot,"logical")) {stop("plot argument is not a logical")}
  if(!is(integrate,"logical")) {stop("integrate argument is not a logical")}
  if(!is(verbose,"logical")) {stop("verbose argument is not a logical")}
  
  # FIX LATER
  if(integrate) { stop("Integration not yet supported")}
  # FIX LATER
  # make sure to throw error if even 1 cellType not in cellIdent
  
  # check if cellIdent and cellType exist within the seurat object
  if(!cellIdent %in% colnames(obj@meta.data)) {stop("cellIdent was not found in obj metadata")}
  if(!all(cellType %in% obj@meta.data[,cellIdent])) {stop("cellType was not found within cellIdent")}
  
  # subset and create new obj based on cellType
  if(verbose) {print("Subsetting data and creating new seurat object")}
  keep <- which(obj@meta.data[,cellIdent] %in% cellType)
  meta <- obj@meta.data[keep,keepCol]
  counts <- obj@assays$RNA@counts[,keep]
  obj <- CreateSeuratObject(counts, project = projectID, meta.data = meta) # overwrite
  
  # Split and SCTransform
  if(verbose) {print("Splitting object by sample and running sctransform")}
  Idents(obj) <- obj$sample
  obj <- SplitObject(obj, split.by = "sample")
  obj <- lapply(obj, SCTransform)
  var.features <- SelectIntegrationFeatures(obj)
  obj.merged <- Merge_Seurat_List(obj)
  VariableFeatures(obj.merged) <- var.features
  remove(obj)
  
  # Run PCA on the merged object
  if(verbose) {print("Running PCA")}
  obj.merged <- RunPCA(object = obj.merged, assay = "SCT")
  min.pc <- find_min_pc(obj.merged)
  
  # Run UMAP and clustering
  if(verbose) {print("Running UMAP and clustering")}
  obj.merged <- RunUMAP(obj.merged,
                        dims = 1:min.pc,
                        assay = "SCT",
                        reduction = "pca",
                        n.components = 3) # set to 3 to use with VR
  obj.merged <- FindNeighbors(object = obj.merged,
                              assay = "SCT",
                              reduction = "pca",
                              dims = 1:min.pc)
  obj.merged <- FindClusters(object = obj.merged,
                             algorithm = 1, # 1 = Louvain
                             resolution = seq(0.1, 0.5,by = 0.1))
  obj.merged$seurat_clusters <- obj.merged$SCT_snn_res.0.5
  
  # reset assay and idents
  DefaultAssay(obj.merged) <- "RNA"
  Idents(obj.merged) <- "seurat_clusters"
  obj.merged <- NormalizeData(obj.merged)
  
  return(obj.merged)
  gc()
}

#' find_cluster_markers
#'
#' @param obj 
#' @param cellIdent 
#'
#' @return
#' @export
#'
#' @examples
find_cluster_markers <- function(obj,cellIdent) {
  
  # load packages
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) { stop("obj argument is not a Seurat object")}
  if(!is(cellIdent,"character")) { stop("cellIdent argument is not a character")}
  
  # find markers
  Idents(obj) <- cellIdent
  all.markers <- FindAllMarkers(object = obj,
                                assay = "RNA",
                                test.use = "MAST",
                                verbose = TRUE)
  
  # add column
  all.markers$delta_pct <- abs(all.markers$pct.1 - all.markers$pct.2)
  
  # rename columns and rows
  rownames(all.markers) <- 1:nrow(all.markers)
  all.markers <- all.markers[,c(6,7,1,5,2:4,8)]
  colnames(all.markers)[c(6,7)] <- c("pct_1","pct_2")
  
  # return data frame
  return(all.markers)
}


#' Title
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
saveToPDF <- function(...){
  d = dev.copy(pdf,...)
  dev.off(d)
}