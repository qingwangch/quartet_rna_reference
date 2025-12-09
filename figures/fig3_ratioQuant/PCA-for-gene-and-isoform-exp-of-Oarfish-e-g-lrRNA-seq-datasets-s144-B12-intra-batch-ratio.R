#' PCA-for-gene-and-isoform-exp-of-Oarfish-e-g-lrRNA-seq-datasets-s144-B12-intra-batch-ratio
#' Qingwang Chen
#' 2025-07-02
#' Description: This script was used for PCA analysis of gene and isoform
#' quantification results from Oarfish

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())

# library import
# library(tidyverse)
library(readr)
library(dplyr)
library(data.table)
library(edgeR)
library(pvca)
library(Biobase)
library(ggplot2)

# Source custom functions
source("/vast/projects/quartet_rna_refdata/analysis/R/figures/functions/calSNR_pca.R")
source("/vast/projects/quartet_rna_refdata/analysis/R/figures/functions/plot_snr.R")

# data import
tx_count_oarfish_G_d12_s144 <- readRDS("./isoforms/Rdata/quartet-LO-minimap2-oarfish-g-txq-count-s144t252913.rds")
# expr.mat.tx <- tx_count_oarfish_G_d12_s144[,-grep("HCC",colnames(tx_count_oarfish_G_d12_s144))]
expr.mat.tx <- cpm(tx_count_oarfish_G_d12_s144) %>% as.data.frame()

# Obtain code information
# Get metadata
metadata <- data.frame(
  sample_id = colnames(expr.mat.tx),
  lib = sapply(strsplit(colnames(expr.mat.tx), "_"), function(x) x[1]),
  tech = sapply(strsplit(colnames(expr.mat.tx), "_"), function(x) x[2]),
  lab_code = sapply(strsplit(colnames(expr.mat.tx), "_"), function(x) x[3]),
  batch = sapply(strsplit(colnames(expr.mat.tx), "_"), function(x) x[4]),
  sample = sapply(strsplit(colnames(expr.mat.tx), "_"), function(x) x[5]),
  rep = sapply(strsplit(colnames(expr.mat.tx), "_"), function(x) x[6])
)

metadata$batch_id <- paste0(metadata$lib, "_", metadata$tech, "_",
                            metadata$lab_code, "_", metadata$batch)

metadata$code <- metadata$sample_id
# Define group function
group <- function(x) {
  group <- rep(NA, length(x))
  group[grepl("D5", x)] <- "D5"
  group[grepl("D6", x)] <- "D6"
  group[grepl("F7", x)] <- "F7"
  group[grepl("M8", x)] <- "M8"
  group[grepl("HCC1395", x)] <- "HCC1395"
  group[grepl("HCC1395BL", x)] <- "HCC1395BL"
  group
}
metadata$group <- group(metadata$code)

#### isoform SNR per batch (All)
metadata$Batch <- metadata$batch_id
ubatch <- unique(metadata$batch_id)

# Ratio--------
# Helper function to calculate ratios based on batch
calculate_ratios <- function(log_expr_mat, metadata, batch_col) {
  ubatch <- unique(as.character(metadata[[batch_col]]))
  ratio_D6 <- matrix(0, ncol = ncol(log_expr_mat), nrow = nrow(log_expr_mat))
  rownames(ratio_D6) <- rownames(log_expr_mat)
  colnames(ratio_D6) <- colnames(log_expr_mat)
  
  for (i in 1:length(ubatch)) {
    name <- as.character(metadata$sample_id[metadata[[batch_col]] == ubatch[i]])
    logexpr_batch <- log_expr_mat[, colnames(log_expr_mat) %in% name]
    m <- rowMeans(logexpr_batch[, grep("D6", colnames(logexpr_batch))])
    mat <- apply(logexpr_batch, 2, function(x) { x - m })
    ratio_D6[rownames(mat), colnames(mat)] <- mat
  }
  
  ratio_D6 <- ratio_D6[apply(ratio_D6, 1, var) != 0, ]  # Remove zero variance rows
  return(ratio_D6)
}

# Loop through each batch to calculate SNR and plot
plots_isoform_q <- list()
snr_isoform_q <- list()
plots_isoform_q_r <- list()
snr_isoform_q_r <- list()

for (batch in ubatch) {
  print(batch)
  # Subset the metadata and expression data for the current batch
  metadata_batch <- metadata[metadata$batch_id == batch,]
  
  # Assuming you are dealing with log2 transformed data adjusted with a pseudocount of 0.01
  expr_mat_cleaned <- na.omit(expr.mat.tx[metadata_batch$sample_id])
  # Set threshold, e.g., retain genes with CPM > 1
  expr_mat_cleaned <- expr_mat_cleaned[rowMeans(expr_mat_cleaned > 1) > 0.5, ]
  
  log_expr_mat <- log2(expr_mat_cleaned + 1)
  # Remove features which variance is zero to ensure the PCA
  log_expr_mat <- log_expr_mat[which(apply(log_expr_mat, 1, var) != 0), ]
  ratio_tx_D6 <- calculate_ratios(log_expr_mat, metadata_batch, "batch_id")
  
  # Calculate SNR
  #metadata_batch$sample <- sapply(metadata_batch$code, sampleq) # Ensure 'sampleq' function is defined
  SNR <- calSNR_pca(log_expr_mat, metadata_batch$sample) # Ensure 'calSNR_pca' function is defined
  SNR_r <- calSNR_pca(ratio_tx_D6, metadata_batch$sample)
  
  # save SNR
  snr_data <- data.frame(batch_id = batch, SNR = SNR$signoise_db)
  snr_isoform_q[[batch]] <- snr_data
  
  snr_data_r <- data.frame(batch_id = batch, SNR = SNR_r$signoise_db)
  snr_isoform_q_r[[batch]] <- snr_data_r
  
  # Plot SNR; ensure 'plot_snr_quartet' function is defined and works as expected
  plot_data <- plot_snr_quartet(SNR, metadata_batch, log_expr_mat, paste(batch, " isoforms s12"), 6, 6)
  p <- plot_data$p# + geom_text_repel(aes(label = sample_id), nudge_x = 0.05, nudge_y = 0.05, size = 2.5, max.overlaps = 100)
  plot_data_r <- plot_snr_quartet(SNR_r, metadata_batch, ratio_tx_D6, paste(batch, " ratio isoforms s12"), 6, 6)
  p_r <- plot_data_r$p
    
  # Store the plot in a list
  plots_isoform_q[[batch]] <- p
  plots_isoform_q_r[[batch]] <- p_r
}

# 合并所有批次的SNR数据为一个数据框
snr_isoform_quartet <- do.call(rbind, snr_isoform_q)
# snr_isoform_all$group <- "all"
snr_isoform_quartet$group <- "quartet"
# 导出SNR数据到CSV文件
write.csv(snr_isoform_quartet, "./results/tables/lrRNAseq_SNR_isoform_quartet_oarfish_s144.csv", row.names = FALSE)

# 查看图片
combined_plot_isoform_q <- plot_grid(
  plots_isoform_q$P_ONT_LN_B2,
  plots_isoform_q$P_ONT_LG_B2,
  plots_isoform_q$D_ONT_LW_B1,
  plots_isoform_q$M_PAB_LG_B2,
  plots_isoform_q$D_ONT_LG_B1,
  plots_isoform_q$P_ONT_LG_B1,
  plots_isoform_q$M_PAB_LN_B1,
  plots_isoform_q$P_ONT_LB_B1,
  plots_isoform_q$I_PAB_LN_B1,
  plots_isoform_q$I_PAB_LN_B2, 
  plots_isoform_q$M_PAB_LG_B1,
  plots_isoform_q$P_ONT_LN_B1,
  nrow = 3, # 设置列数为 1，即垂直排列
  align = 'v' # 垂直对齐图形
)
combined_plot_isoform_q

ggsave("./results/plots/merged_PCA_oarfish_g4_tx_s144_B12.pdf", combined_plot_isoform_q, width = 12, height = 9, dpi = 300)
ggsave("./results/plots/merged_PCA_oarfish_g4_tx_s144_B12.png", combined_plot_isoform_q, width = 12, height = 9, dpi = 300)

# 合并所有批次的SNR数据为一个数据框
snr_isoform_quartet_r <- do.call(rbind, snr_isoform_q_r)
# snr_isoform_all$group <- "all"
snr_isoform_quartet_r$group <- "quartet"
# 导出SNR数据到CSV文件
write.csv(snr_isoform_quartet_r, "./results/tables/lrRNAseq_SNR_isoform_quartet_oarfish_ratio_s144.csv", row.names = FALSE)


# 查看图片
combined_plot_isoform_q_r <- plot_grid(
  plots_isoform_q_r$P_ONT_LN_B2,
  plots_isoform_q_r$P_ONT_LG_B2,
  plots_isoform_q_r$D_ONT_LW_B1,
  plots_isoform_q_r$M_PAB_LG_B2,
  plots_isoform_q_r$D_ONT_LG_B1,
  plots_isoform_q_r$P_ONT_LG_B1,
  plots_isoform_q_r$M_PAB_LN_B1,
  plots_isoform_q_r$P_ONT_LB_B1,
  plots_isoform_q_r$I_PAB_LN_B1,
  plots_isoform_q_r$I_PAB_LN_B2, 
  plots_isoform_q_r$M_PAB_LG_B1,
  plots_isoform_q_r$P_ONT_LN_B1,
  nrow = 3, # 设置列数为 1，即垂直排列
  align = 'v' # 垂直对齐图形
)
combined_plot_isoform_q_r
ggsave("./results/plots/merged_PCA_oarfish_g4_tx_ratio_s144_B12.pdf", combined_plot_isoform_q_r, width = 12, height = 9, dpi = 300)
ggsave("./results/plots/merged_PCA_oarfish_g4_tx_ratio_s144_B12.png", combined_plot_isoform_q_r, width = 12, height = 9, dpi = 300)

# combined_plot_isoform_qc <- plot_grid(
#   plots_isoform_q$P_ONT_LN_B2,
#   plots_isoform_q$D_ONT_LG_B1,
#   plots_isoform_q$P_ONT_LG_B2,
#   plots_isoform_q$P_ONT_LG_B1,
#   plots_isoform_q$M_PAB_LN_B1,
#   plots_isoform_q$M_PAB_LG_B2,
#   nrow = 2, # 设置列数为 1，即垂直排列
#   align = 'v' # 垂直对齐图形
# )
# combined_plot_isoform_qc
# ggsave("./results/plots/merged_PCA_oarfish_tx_s72_B6.pdf", combined_plot_isoform_qc, width = 9, height = 6, dpi = 300)
# ggsave("./results/plots/merged_PCA_oarfish_tx_s72_B6.png", combined_plot_isoform_qc, width = 9, height = 6, dpi = 300)

##### SNR11
for (batch in ubatch) {
  print(batch)
  # Subset the metadata and expression data for the current batch
  metadata_batch <- metadata[metadata$batch_id == batch,]
  
  # Assuming you are dealing with log2 transformed data adjusted with a pseudocount of 0.01
  expr_mat_cleaned <- na.omit(expr.mat.tx[metadata_batch$sample_id])
  # Set threshold, e.g., retain genes with CPM > 1
  expr_mat_cleaned <- expr_mat_cleaned[rowMeans(expr_mat_cleaned > 1) > 0.5, ]
  
  log_expr_mat <- log2(expr_mat_cleaned + 1)
  # Remove features with zero variance to ensure PCA works
  log_expr_mat <- log_expr_mat[which(apply(log_expr_mat, 1, var) != 0), ]
  
  # Initialize list to store SNR values for each sample removal
  snr_values <- c()
  
  # Loop over each sample in the batch, removing one sample at a time
  for (sample_id in metadata_batch$sample_id) {
    # Remove the current sample and subset the data accordingly
    expr_mat_subset <- log_expr_mat[, !colnames(log_expr_mat) %in% sample_id]
    metadata_subset <- metadata_batch[metadata_batch$sample_id != sample_id,]
    print(sample_id)
    # Calculate SNR for the remaining samples
    SNR <- calSNR_pca(expr_mat_subset, metadata_subset$sample) # Ensure 'calSNR_pca' function is defined
    
    # Store the SNR value along with the batch and removed sample
    snr_values <- c(snr_values, SNR$signoise_db)
  }
  
  # Create a data frame for the SNR values of this batch
  snr_data <- data.frame(batch_id = batch, SNR11 = snr_values, sample_removed = metadata_batch$sample_id)
  
  # Store the SNR values for this batch
  snr_isoform_q[[batch]] <- snr_data
}

# Combine all batch SNR data
snr_isoform_q_combined <- do.call(rbind, snr_isoform_q)
# Save the combined SNR data to CSV
write.csv(snr_isoform_q_combined, "./results/tables/lrRNAseq_oarfish_SNR_isoform_quartet_s144_SNR12.csv", row.names = FALSE)

# plot SNR11
snr_isoform_q_combined <- merge(snr_isoform_q_combined,snr_isoform_quartet)
snr_isoform_q_combined$SNR <- round(snr_isoform_q_combined$SNR,2)
p1 <- ggplot(snr_isoform_q_combined, aes(x = reorder(batch_id,-SNR),fill=batch_id)) +
  geom_bar(aes(y=SNR),stat = "identity",color="black", position = position_dodge(width = 0.3),
           linewidth=0.7 ,width = 0.8) +
  geom_point(aes(y =SNR11), size = 3,shape=21) +
  geom_text(aes(label = scales::comma(SNR), y = SNR),  # 添加这一行来显示标签
            position = position_dodge(width = 0.3),
            size = 4,
            vjust = -2.5) +  # 标签的位置微调
  geom_hline(yintercept = mean(snr_isoform_q_combined$SNR), linetype = "dashed", color = "#FF2D00") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold",size = 18),
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", color = "black",size = 14),
        legend.title = element_text(color = "black", face = "bold"),
        legend.text = element_text(color = "black", face = "bold"),
        legend.position = "none",
        legend.background = element_blank(),
        panel.grid = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        #   panel.background = element_rect(fill = c("gray","blue"), colour = NA),
        strip.background.x = element_rect(fill = "white", colour = "white"),
        strip.text.x = element_text(colour = "black", face = "bold",size = 16),
        strip.placement = "inside",
        strip.switch.pad.grid = unit(1, "inch")) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Batch", y = "SNR", title = "SNR11 isoforms");p1
ggsave("./results/plots/barplot_oarfish_SNR11_s144_B12.png", p1, width = 12, height = 5, dpi = 300)
ggsave("./results/plots/barplot_oarfish_SNR11_s144_B12.pdf", p1, width = 12, height = 5, dpi = 300)

# remove all
rm(list=ls())
gc()
