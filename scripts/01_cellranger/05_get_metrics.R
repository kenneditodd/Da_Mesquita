# load libraries
library(stringr)

# file locations
locations <- c("/Da_Mesquita/count/E3_2M_F/outs/metrics_summary.csv",
               "/Da_Mesquita/count/E3_2M_M/outs/metrics_summary.csv",
               "/Da_Mesquita/count/E3_14M_F/outs/metrics_summary.csv",
               "/Da_Mesquita/count/E3_14M_M/outs/metrics_summary.csv",
               "/Da_Mesquita/count/E4_2M_F/outs/metrics_summary.csv",
               "/Da_Mesquita/count/E4_2M_M/outs/metrics_summary.csv",
               "/Da_Mesquita/count/E4_14M_F/outs/metrics_summary.csv",
               "/Da_Mesquita/count/E4_14M_M/outs/metrics_summary.csv")

# sample names
names <- str_match(locations, "/Da_Mesquita/count/(.+)/outs/metrics_summary.csv")[,2]

# initialize df and loop through files
df <- data.frame()
for (i in 1:length(locations)) {
  if (i == 1) {
    df <- read.csv(locations[i])
  } else {
    row <- read.csv(locations[i])[1,]
    df <- rbind(df,row)
  }
}

rownames(df) <- names
c.names <- c("estimated_cells", "mean_reads", "median_genes", "number_reads",
                  "valid_barcodes", "sequencing_saturation", "Q30_bases_barcode",
                  "Q30_bases_read", "Q30_bases_UMI", "reads_mapped_genome", "confident_reads_mapped_genome",
                  "confident_intergenic_reads_mapped", "confident_intronic_reads_mapped",
                  "confident_exonic_reads_mapped", "confident_reads_mapped_transcriptome",
                  "reads_mapped_antisense", "fraction_reads", "total_genes", "median_UMI")
colnames(df) <- c.names

write.table(df, 
            "/Da_Mesquita/count/web_summaries/overall_metrics.tsv",
            sep = "\t",
            quote = FALSE)


