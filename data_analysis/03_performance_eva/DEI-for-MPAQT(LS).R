#' DEI-for-MPAQT(LS)
#' Qingwang Chen
#' 2025-04-18

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# 指定新的库路径
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# 修改 .libPaths，使其优先使用自定义路径
.libPaths(c(custom_lib_path, .libPaths()))

# 检查当前的库路径
print(.libPaths())

# library import
library(limma)
library(edgeR)
library(dplyr)
library(stringr)
library(data.table)

# 输入路径
tpm_file <- "/vast/projects/quartet_rna_refdata/analysis/mpaqt/projects/quartet/transcript_tpm_LR_SR_s144.csv"
output_dir <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform_tpm"
pipeline <- "MPAQT"
group1 <- "D6"
group2_list <- c("D5", "F7", "M8")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 加载 TPM 表达矩阵
tpm <- read.csv(tpm_file)
rownames(tpm) <- tpm$transcript_id
tpm$transcript_id <- NULL

# 样本名列表
sample_names <- colnames(tpm)
batch_ids <- unique(sapply(strsplit(sample_names, "_"), function(x) paste(x[1:4], collapse = "_")))

# 初始化结果列表
all_results <- list()

# 主函数
for (batch_id in batch_ids) {
  samples_in_batch <- sample_names[grepl(batch_id, sample_names)]
  batch_expr <- tpm[, samples_in_batch]
  
  group_labels <- sapply(strsplit(samples_in_batch, "_"), function(x) x[5])
  
  for (group2 in group2_list) {
    compare_samples <- which(group_labels %in% c(group1, group2))
    expr_filtered <- batch_expr[, compare_samples]
    
    # 明确 factor 水平顺序
    group <- factor(group_labels[compare_samples], levels = c(group1, group2))  # D6 为 baseline，对应 group1
    
    # 构建设计矩阵
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    
    # 构造对比组
    contrast <- makeContrasts(contrasts = paste0(group2, "-", group1), levels = design)
    
    # voom + limma
    v <- voom(expr_filtered, design, normalize.method = "quantile")
    fit <- lmFit(v, design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    
    # 提取结果
    res <- topTable(fit2, number = Inf, sort.by = "P")
    res$comparison <- paste0(group2, "vs", group1)
    res$feature <- rownames(res)
    res$batch_id <- batch_id
    res$pipeline <- pipeline
    res$feature_type <- "non-DE"
    res$feature_type[which(res$adj.P.Val < 0.05 & res$logFC >= 1)] <- "up-regulate"
    res$feature_type[which(res$adj.P.Val < 0.05 & res$logFC <= -1)] <- "down-regulate"
    
    # 保存单个文件
    output_file <- file.path(output_dir, paste0("DE_results_", pipeline, "_", batch_id, "_", group2, "vs", group1, ".txt"))
    write.table(res, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Saved:", output_file, "\n")
    
    # 存入总列表
    all_results[[paste(batch_id, group2, sep = "_")]] <- res
  }
}

# 合并所有结果并保存
merged_results <- dplyr::bind_rows(all_results)
test_results <- fread("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_Salmon_isoform.txt") %>% as.data.frame()
# 统一列名
colnames(merged_results)[colnames(merged_results) == "P.Value"] <- "PValue"
colnames(merged_results)[colnames(merged_results) == "adj.P.Val"] <- "FDR"
merged_results$analysis_type <- "isoform"
head(merged_results)
merged_output_file <- file.path(output_dir, paste0("DE_results_", pipeline, "_merged_all.txt"))
write.table(merged_results, merged_output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("All merged results saved to:", merged_output_file, "\n")

# remove all
rm(list=ls())
gc()
