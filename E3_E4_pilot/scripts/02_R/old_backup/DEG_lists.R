setwd("/research/labs/neurology/fryer/m214960/Da_Mesquita/scripts/R")

isoform.df <- read.delim("../../results/DEGs/E4_vs_E3_DEGs.tsv",
                         sep = "\t")
isoform.df <- isoform.df[isoform.df$cluster == "Dendritic cells",]
isoform.df <- isoform.df[isoform.df$p_val_adj < 0.05,]
isoform.up <- isoform.df[isoform.df$avg_log2FC > 0,]
isoform.up <- isoform.up[,2]

write.table(isoform.up,
            "../../results/DEGs/DEG_lists/isoform_dendritic_cells_up.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)