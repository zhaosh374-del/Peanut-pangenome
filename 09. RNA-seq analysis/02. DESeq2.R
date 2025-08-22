library(DESeq2)

sample_data <- read.csv(sample_info_file, header = TRUE)
sample_data$group <- factor(sample_data$group)
sample_data$group <- relevel(sample_data$group, ref = reference_level)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_data,
  design = ~ group
)

dds_analyzed <- DESeq(dds)

results <- results(dds_analyzed)

results_clean <- na.omit(results)
downregulated <- results_clean[results_clean$log2FoldChange <= -lfc_threshold & 
                               results_clean$padj <= padj_threshold, ]
upregulated <- results_clean[results_clean$log2FoldChange >= lfc_threshold & 
                             results_clean$padj <= padj_threshold, ]

significant_genes <- rbind(upregulated, downregulated)
write.csv(significant_genes, file = paste0(output_dir, sample_name, ".sig.csv"))
write.csv(results, file = paste0(output_dir, sample_name, ".csv"))
