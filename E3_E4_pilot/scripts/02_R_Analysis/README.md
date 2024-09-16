# R Analysis Scripts
- *01_processing.Rmd* will read counts data, perform prefiltering QC, quality filtering, postfiltering QC, and clustering.
- *02_pass1_annotation.Rmd* annotates the clusters created from 01_processing.Rmd.
- *03_pass1_quality_reclustering.Rmd* reclusters data into two groups. This was to further inspect problematic clusters (i.e. ambient RNA clusters) identified in 02_pass_1_annotation.
- *04_pass2_recluster_and_annotation.Rmd* removes problematic clusters identified in 02_pass1_annotation.Rmd and 03_pass1_quality_reclustering.Rmd and then reclusters the remaining cells. This creates a second, final, cleaned pass of the data with all cell types.
- *05_pass2_QC_and_DE.Rmd* takes peforms QC and differential expression on the final annotated model create in 04_pass2_recluster_and_annotation.Rmd.
- *06_pass2_subcluster.Rmd* will recluster specific cell types to identify subpopulations.
- *07_pass2_endothelial_QC_and_DE.Rmd* performs QC and differential expression on the reclustered endothelial cells created in 06_pass2_subcluster.Rmd. 
- *08_pass2_macrophages_monocytes_QC_and_DE.Rmd* performs QC and differential expression on the reclustered macrophages & monocytes created in 06_pass2_subcluster.Rmd. 
