#' PCC-absolute-relative-long-read-from-oarfish-s144
#' Qingwang Chen
#' 2025-07-02
#' 

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())

# library import
library(readr)
library(dplyr)
library(data.table)
library(pheatmap)
library(edgeR)
library(pvca)
library(Biobase)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
# library(ggVennDiagram)
library(cowplot)
library(tidysdm)
library(patchwork)

# data import
expr.mat.tx <- readRDS("./isoforms/Rdata/quartet-LO-minimap2-oarfish-g-txq-count-s144t252913.rds")

# combined_expression_data <- cpm(expr.mat.tx)
# get metadata 
columns <- colnames(expr.mat.tx)

metadata <- data.frame(
  sample_id = columns,
  lib = sapply(strsplit(columns, "_"), `[`, 1),
  platform = sapply(strsplit(columns, "_"), `[`, 2),
  lab = sapply(strsplit(columns, "_"), `[`, 3),
  batch = sapply(strsplit(columns, "_"), `[`, 4),
  sample = sapply(strsplit(columns, "_"), `[`, 5),
  rep = sapply(strsplit(columns, "_"), `[`, 6)
)
metadata$batch_id <- paste(metadata$lib,metadata$platform,metadata$lab,metadata$batch,sep = "_")
metadata$sample_rep <- paste(metadata$sample,metadata$rep,sep = "_")
metadata$type <- metadata$sample
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
metadata$batch <- metadata$batch_id
metadata$Batch <- metadata$batch_id

### feature selection
# generate_profiles <- function(metaMat, counts, denominator, threshold_pct = 0.5) {
#   ###### First perform feature selection before obtaining absolute and ratio profiles
#   ## Otherwise, cbind will fill with many NA values and the PCA plot will appear as two lines
#   filtered_genes <- list()
#   for (batch_id in unique(metaMat$batch)) {
#     print(paste("Processing batch:", batch_id))
#     # Extract the columns and group information for the current batch
#     cols <- metaMat %>% filter(batch == batch_id) %>% pull(sample_id)
#     group <- metaMat %>% filter(batch == batch_id) %>% pull(sample)
#     subcounts <- counts[, cols]
#     # Denominator stable expression filter: at least 2 samples with counts >= 2
#     denom_cols <- cols[group == denominator]
#     denom_filter <- rowSums(subcounts[, denom_cols, drop = FALSE] >= 2) >= 2
#     # Overall expression filter: expressed in more than 50% of the samples
#     expressed_filter <- rowSums(subcounts > 0) >= ncol(subcounts) * threshold_pct
#     # Store the gene row names for the current batch
#     filtered_genes[[batch_id]] <- rownames(subcounts)[denom_filter & expressed_filter]
#   }
#   # unique(unlist(filtered_genes))  # Union: approximately 52k genes
#   # Reduce(intersect, filtered_genes)  # Intersection: approximately 14k genes
#   # Get genes that satisfy the condition in at least 50% of batches, approximately 28k genes
#   gene_counts <- table(unlist(filtered_genes))
#   num_batches <- length(filtered_genes)
#   threshold <- num_batches * threshold_pct
#   filtered_genes_50pct <- names(gene_counts[gene_counts >= threshold])
#   return(filtered_genes_50pct)
# }

## expressed isoforms
# expr.mat.tx.expressed <- generate_profiles(metadata, expr.mat.tx, "D6", 0.5)

# read features
rt <- fread("/vast/projects/quartet_rna_refdata/analysis/oarfish_quant/results/tables/rt_tx.2.csv") %>% as.data.frame()
expr.mat.tx.expressed <- rt$intial_rn[rt$Pobject==1]
enst_ids <- vapply(
  strsplit(expr.mat.tx.expressed, "\\|", fixed = FALSE),
  `[`, character(1), 1
)
# expr.mat.tx.expressed_cpm <- cpm(expr.mat.tx[expr.mat.tx.expressed,])
expr.mat.tx.expressed_cpm <- cpm(expr.mat.tx[enst_ids,])
log_expr.mat.tx.expressed_cpm <- apply(expr.mat.tx.expressed_cpm, 2, function(x) {log2(x + 1)})

# cor
# retain the features whose count >=3 above 90% samples
# retained_features <- rownames(count_merge[rowMeans(count_merge >= 3) > 0.9, ])

# 计算样本之间的相关性（使用皮尔逊相关系数）
cor_matrix_r <- cor(log_expr.mat.tx.expressed_cpm, method = "pearson")
# cor_matrix_r <- cor(ratio_tx_D6, method = "pearson")
head(cor_matrix_r)

# 选择颜色映射
colormap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 从 metadata 中提取需要的注释列
annotations <- metadata[, c("lab", "platform", "lib", "sample")]
rownames(annotations) <- metadata$sample_id

# 创建注释颜色列表
annotation_colors <- list(
  lab = c("LB" = "#d35400",  # Burnt orange instead of red
          "LG" = "#28a745",
          "LN" = "#8e44ad",
          "LW" = "#2980b9"), # Moderate blue instead of blue
  platform = c("ONT" = "#6c9eb4",
               "PAB" = "#e01a97"),
  lib = c("P" = "#FFDD44",
          "I" = "#16a085",  # Teal instead of green
          "M" = "#A9DFBF",
          "D" = "#85C1AE"),
  sample = c("D5" = "#4CC3D9",
             "D6" = "#7BC8A4",
             "F7" = "#FFC65D",
             "M8" = "#F16745")
)


# 绘制热图并添加注释
p1 <- pheatmap(mat = cor_matrix_r,                    # 数据矩阵
               color = colormap,  # 颜色映射
               cluster_rows = TRUE,                 # 行聚类
               cluster_cols = TRUE,                 # 列聚类
               # scale = "row",                       # 按行标准化
               show_rownames = F,                # 显示行名
               show_colnames = F,                # 显示列名
               fontsize = 12,                       # 字体大小
               fontsize_row = 10,                   # 行名字体大小
               fontsize_col = 10,                   # 列名字体大小
               border_color = "white",              # 单元格边框颜色
               legend = TRUE,                       # 显示色条
               annotation_col = annotations,        # 添加列注释
               annotation_colors = annotation_colors, # 注释颜色
               main = "PCC of absolute expression LR s144") # 标题
p1
# 保存为 PDF 文件
ggsave("./results/plots/cor-absolute-isoform-expression-oarfish-lr-s144.png", plot = p1, width = 9, height = 8)

## ratio
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

ratio_tx_D6 <- calculate_ratios(log_expr.mat.tx.expressed_cpm, metadata, "batch_id")

# 计算样本之间的相关性（使用皮尔逊相关系数）
cor_matrix_rr <- cor(ratio_tx_D6[,-grep("D6",colnames(ratio_tx_D6))], method = "pearson")
# cor_matrix_r <- cor(ratio_tx_D6, method = "pearson")
head(cor_matrix_rr)

# 选择颜色映射
colormap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 从 metadata 中提取需要的注释列
annotations <- metadata[, c("lab", "platform", "lib", "sample")]
rownames(annotations) <- metadata$sample_id


# 绘制热图并添加注释
p2 <- pheatmap(mat = cor_matrix_rr,                    # 数据矩阵
               color = colormap,  # 颜色映射
               cluster_rows = TRUE,                 # 行聚类
               cluster_cols = TRUE,                 # 列聚类
               # scale = "row",                       # 按行标准化
               show_rownames = F,                # 显示行名
               show_colnames = F,                # 显示列名
               fontsize = 12,                       # 字体大小
               fontsize_row = 10,                   # 行名字体大小
               fontsize_col = 10,                   # 列名字体大小
               border_color = "white",              # 单元格边框颜色
               legend = TRUE,                       # 显示色条
               annotation_col = annotations,        # 添加列注释
               annotation_colors = annotation_colors, # 注释颜色
               main = "PCC of relative expression LR s144") # 标题
p2
# 保存为 PDF 文件
ggsave("./results/plots/cor-relative-isoform-expression-oarfish-lr-s144.png", plot = p2, width = 9, height = 8)

# 使用 cowplot 共享图例
# 提取 pheatmap 的 grob
p1_grob <- ggplotify::as.ggplot(p1$gtable)
p2_grob <- ggplotify::as.ggplot(p2$gtable)
combined_plot <- plot_grid(
  p1_grob + theme(legend.position = "none"),  # 移除 p1 的图例
  p2_grob + theme(legend.position = "none"),  # 移除 p2 的图例
  align = "h",  
  ncol = 2,     # 单列布局
  labels = c("")  # 添加子图标签
)
combined_plot

ggsave("./results/plots/abs-cor-relative-isoform-expression-lr-s144.pdf", plot = combined_plot, width = 10, height = 5)

# correlation by batch id
# Prepare data for the absolute PCC (cor_matrix)
cor_df_absolute <- as.data.frame(as.table(cor_matrix_r))
colnames(cor_df_absolute) <- c("Sample1", "Sample2", "PCC")
test1 <- cor_df_absolute[cor_df_absolute$PCC<0,]
# 删除 Sample1 和 Sample2 完全相同的行
cor_df_absolute <- cor_df_absolute %>%
  filter(Sample1 != Sample2)

# 更新分类逻辑
cor_df_absolute <- cor_df_absolute %>%
  mutate(
    # 提取批次信息，假设批次信息在“_B1”前
    batch_1 = sapply(strsplit(as.character(Sample1), "_"), function(x) paste(x[1:4], collapse = "_")),
    batch_2 = sapply(strsplit(as.character(Sample2), "_"), function(x) paste(x[1:4], collapse = "_")),
    
    # 提取样本类型信息，假设样本类型为后缀（例如，D5, D6, F7等）
    sample_type_1 = gsub(".*_(D[0-9]+|F[0-9]+|M[0-9]+).*", "\\1", Sample1),
    sample_type_2 = gsub(".*_(D[0-9]+|F[0-9]+|M[0-9]+).*", "\\1", Sample2),
    
    # 提取protocol信息
    protocol1 = sapply(strsplit(as.character(Sample1), "_"), function(x) x[1]),
    protocol2 = sapply(strsplit(as.character(Sample2), "_"), function(x) x[1]),
    
    # 提取platform信息
    platform1 = sapply(strsplit(as.character(Sample1), "_"), function(x) x[2]),
    platform2 = sapply(strsplit(as.character(Sample2), "_"), function(x) x[2]),
    
    # 提取lab信息
    lab1 = sapply(strsplit(as.character(Sample1), "_"), function(x) x[3]),
    lab2 = sapply(strsplit(as.character(Sample2), "_"), function(x) x[3]),
    
    # Class逻辑
    Class = case_when(
      # Cross-platform
      platform1 != platform2 ~ "Cross-platform",
      # Cross-protocol
      protocol1 != protocol2 ~ "Cross-protocol",
      # Cross-lab
      lab1 != lab2 ~ "Cross-lab",
      # Intra-batch
      batch_1 == batch_2 ~ "Intra-batch",
      # Cross-batch
      batch_1 != batch_2 ~ "Cross-time",
      # Others
      TRUE ~ "Others"
    ),
    
    # Group 逻辑
    Group = case_when(
      # Cross-platform 分为 intra-sample 和 cross-sample
      Class == "Cross-platform" & sample_type_1 == sample_type_2 ~ "Cross-platform intra-sample",
      Class == "Cross-platform" & sample_type_1 != sample_type_2 ~ "Cross-platform cross-sample",
      
      # Cross-protocol 分为 intra-sample 和 cross-sample
      Class == "Cross-protocol" & sample_type_1 == sample_type_2 ~ "Cross-protocol intra-sample",
      Class == "Cross-protocol" & sample_type_1 != sample_type_2 ~ "Cross-protocol cross-sample",
      
      # Cross-lab 分为 intra-sample 和 cross-sample
      Class == "Cross-lab" & sample_type_1 == sample_type_2 ~ "Cross-lab intra-sample",
      Class == "Cross-lab" & sample_type_1 != sample_type_2 ~ "Cross-lab cross-sample",
      
      # Intra-batch 分为 intra-sample 和 cross-sample
      Class == "Intra-batch" & sample_type_1 == sample_type_2 ~ "Intra-batch intra-sample",
      Class == "Intra-batch" & sample_type_1 != sample_type_2 ~ "Intra-batch cross-sample",
      
      # Cross-batch 分为 intra-sample 和 cross-sample
      Class == "Cross-time" & sample_type_1 == sample_type_2 ~ "Cross-time intra-sample",
      Class == "Cross-time" & sample_type_1 != sample_type_2 ~ "Cross-time cross-sample",
      
      # 其他
      TRUE ~ "Others"
    )
  )

# Statistic
cor_df_absolute %>%
  group_by(Class) %>%
  summarize(Count = n())

cor_df_absolute %>%
  group_by(Group) %>%
  summarize(Count = n())

# 创建新的列：GroupType 用于区分 intra-sample 和 cross-sample
cor_df_absolute$GroupType <- ifelse(grepl("intra-sample", cor_df_absolute$Group), "intra-sample", "cross-sample")

# Create the split violin plot
# Calculate the median PCC for each Class
class_order <- cor_df_absolute %>%
  group_by(Class) %>%
  summarize(median_PCC = median(PCC, na.rm = TRUE)) %>%
  arrange(median_PCC) %>%
  pull(Class)

# Convert 'Class' to a factor with levels ordered by median PCC
cor_df_absolute <- cor_df_absolute %>%
  mutate(Class = factor(Class, levels = class_order))

library(ggplot2)
library(scales)

# 1. 按 PCC 中位数对 Class 排序
cor_df_absolute$Class <- reorder(cor_df_absolute$Class, cor_df_absolute$PCC, FUN = median)

# 2. 定义 GroupType 顺序（图例中先 cross-sample，再 intra-sample）
group_order <- c("cross-sample", "intra-sample")

# 3. 创建合并因子 ClassGroup，因子水平顺序对应“Class × GroupType”
cor_df_absolute$ClassGroup <- factor(
  paste(cor_df_absolute$Class, cor_df_absolute$GroupType, sep = " - "),
  levels = unlist(lapply(levels(cor_df_absolute$Class), function(cl) {
    paste(cl, group_order, sep = " - ")
  }))
)

# 4. 为每个 Class 指定一个基础颜色（这里给出了 5 种，可自行增减）
unique_classes <- levels(cor_df_absolute$Class)
base_colors <- c("#FB9A99", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")[seq_along(unique_classes)]

# 5. 生成颜色映射：同一 Class 的 cross-sample 半透明、intra-sample 不透明
color_mapping <- sapply(levels(cor_df_absolute$ClassGroup), function(cg) {
  parts <- strsplit(cg, " - ")[[1]]
  cl <- parts[1]
  gt <- parts[2]
  col <- base_colors[which(unique_classes == cl)]
  if(gt == "cross-sample") alpha(col, 0.3) else col
})
names(color_mapping) <- levels(cor_df_absolute$ClassGroup)

# 6. 绘图：左右小提琴 + 中间箱线图
p3 <- ggplot(cor_df_absolute, aes(x = Class, y = PCC)) +
  geom_split_violin(aes(fill = ClassGroup), trim = TRUE, position = position_dodge(width = 0.2)) +
  scale_fill_manual(values = color_mapping) +
  labs(
    title = "Before correction",
    x = "Class",
    y = "Pearson correlation coefficient",
    fill = "Class - Group Type"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

print(p3)

# 保存为 PDF 文件
ggsave("./results/plots/pcc-boxplot-absolute-isoform-expression-lr-s144.png", plot = p3, width = 7, height = 6)
ggsave("./results/plots/pcc-boxplot-absolute-isoform-expression-lr-s144.pdf", plot = p3, width = 7, height = 6)

# 若 cor_df_absolute 不是 data.table，可先转换
setDT(cor_df_absolute)

# 计算每个 ClassGroup 的 PCC 中位数
pcc_median <- cor_df_absolute[
  , .(median_PCC = median(PCC, na.rm = TRUE)),
  by = ClassGroup
][order(ClassGroup)]

print(pcc_median)

# Prepare data for the ratio PCC (cor_matrix_rr)
cor_df_relative <- as.data.frame(as.table(cor_matrix_rr))
colnames(cor_df_relative) <- c("Sample1", "Sample2", "PCC")
test1 <- cor_df_relative[cor_df_relative$PCC<0,]
# 删除 Sample1 和 Sample2 完全相同的行
cor_df_relative <- cor_df_relative %>%
  filter(Sample1 != Sample2)

# 更新分类逻辑
cor_df_relative <- cor_df_relative %>%
  mutate(
    # 提取批次信息，假设批次信息在“_B1”前
    batch_1 = sapply(strsplit(as.character(Sample1), "_"), function(x) paste(x[1:4], collapse = "_")),
    batch_2 = sapply(strsplit(as.character(Sample2), "_"), function(x) paste(x[1:4], collapse = "_")),
    
    # 提取样本类型信息，假设样本类型为后缀（例如，D5, D6, F7等）
    sample_type_1 = gsub(".*_(D[0-9]+|F[0-9]+|M[0-9]+).*", "\\1", Sample1),
    sample_type_2 = gsub(".*_(D[0-9]+|F[0-9]+|M[0-9]+).*", "\\1", Sample2),
    
    # 提取protocol信息
    protocol1 = sapply(strsplit(as.character(Sample1), "_"), function(x) x[1]),
    protocol2 = sapply(strsplit(as.character(Sample2), "_"), function(x) x[1]),
    
    # 提取platform信息
    platform1 = sapply(strsplit(as.character(Sample1), "_"), function(x) x[2]),
    platform2 = sapply(strsplit(as.character(Sample2), "_"), function(x) x[2]),
    
    # 提取lab信息
    lab1 = sapply(strsplit(as.character(Sample1), "_"), function(x) x[3]),
    lab2 = sapply(strsplit(as.character(Sample2), "_"), function(x) x[3]),
    
    # Class逻辑
    Class = case_when(
      # Cross-platform
      platform1 != platform2 ~ "Cross-platform",
      # Cross-protocol
      protocol1 != protocol2 ~ "Cross-protocol",
      # Cross-lab
      lab1 != lab2 ~ "Cross-lab",
      # Intra-batch
      batch_1 == batch_2 ~ "Intra-batch",
      # Cross-batch
      batch_1 != batch_2 ~ "Cross-time",
      # Others
      TRUE ~ "Others"
    ),
    
    # Group 逻辑
    Group = case_when(
      # Cross-platform 分为 intra-sample 和 cross-sample
      Class == "Cross-platform" & sample_type_1 == sample_type_2 ~ "Cross-platform intra-sample",
      Class == "Cross-platform" & sample_type_1 != sample_type_2 ~ "Cross-platform cross-sample",
      
      # Cross-protocol 分为 intra-sample 和 cross-sample
      Class == "Cross-protocol" & sample_type_1 == sample_type_2 ~ "Cross-protocol intra-sample",
      Class == "Cross-protocol" & sample_type_1 != sample_type_2 ~ "Cross-protocol cross-sample",
      
      # Cross-lab 分为 intra-sample 和 cross-sample
      Class == "Cross-lab" & sample_type_1 == sample_type_2 ~ "Cross-lab intra-sample",
      Class == "Cross-lab" & sample_type_1 != sample_type_2 ~ "Cross-lab cross-sample",
      
      # Intra-batch 分为 intra-sample 和 cross-sample
      Class == "Intra-batch" & sample_type_1 == sample_type_2 ~ "Intra-batch intra-sample",
      Class == "Intra-batch" & sample_type_1 != sample_type_2 ~ "Intra-batch cross-sample",
      
      # Cross-batch 分为 intra-sample 和 cross-sample
      Class == "Cross-time" & sample_type_1 == sample_type_2 ~ "Cross-time intra-sample",
      Class == "Cross-time" & sample_type_1 != sample_type_2 ~ "Cross-time cross-sample",
      
      # 其他
      TRUE ~ "Others"
    )
  )

# Statistic
cor_df_relative %>%
  group_by(Class) %>%
  summarize(Count = n())

cor_df_relative %>%
  group_by(Group) %>%
  summarize(Count = n())

# 创建新的列：GroupType 用于区分 intra-sample 和 cross-sample
cor_df_relative$GroupType <- ifelse(grepl("intra-sample", cor_df_relative$Group), "intra-sample", "cross-sample")

# Create the split violin plot
# Calculate the median PCC for each Class
# class_order <- cor_df_relative %>%
#   group_by(Class) %>%
#   summarize(median_PCC = median(PCC, na.rm = TRUE)) %>%
#   arrange(median_PCC) %>%
#   pull(Class)

# Convert 'Class' to a factor with levels ordered by median PCC
cor_df_relative <- cor_df_relative %>%
  mutate(Class = factor(Class, levels = class_order))

# cor_df_relative$Class <- reorder(cor_df_relative$Class, cor_df_relative$PCC, FUN = median)

# 先获取图3最终使用的 Class 顺序
order_in_figure3 <- levels(cor_df_absolute$Class)

# 让图4的数据 cor_df_relative$Class 拥有与图3相同的因子水平
cor_df_relative$Class <- factor(cor_df_relative$Class, levels = order_in_figure3)

# 2. 定义 GroupType 顺序（图例中先 cross-sample，再 intra-sample）
group_order <- c("cross-sample", "intra-sample")

# 3. 创建合并因子 ClassGroup，因子水平顺序对应“Class × GroupType”
cor_df_relative$ClassGroup <- factor(
  paste(cor_df_relative$Class, cor_df_relative$GroupType, sep = " - "),
  levels = unlist(lapply(levels(cor_df_relative$Class), function(cl) {
    paste(cl, group_order, sep = " - ")
  }))
)

# 4. 为每个 Class 指定一个基础颜色（这里给出了 5 种，可自行增减）
unique_classes <- levels(cor_df_relative$Class)
base_colors <- c("#FB9A99", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")[seq_along(unique_classes)]

# 5. 生成颜色映射：同一 Class 的 cross-sample 半透明、intra-sample 不透明
color_mapping <- sapply(levels(cor_df_relative$ClassGroup), function(cg) {
  parts <- strsplit(cg, " - ")[[1]]
  cl <- parts[1]
  gt <- parts[2]
  col <- base_colors[which(unique_classes == cl)]
  if(gt == "cross-sample") alpha(col, 0.3) else col
})
names(color_mapping) <- levels(cor_df_relative$ClassGroup)

# 6. 绘图：左右小提琴 + 中间箱线图
p4 <- ggplot(cor_df_relative, aes(x = Class, y = PCC)) +
  geom_split_violin(aes(fill = ClassGroup), trim = TRUE, position = position_dodge(width = 0.2)) +
  scale_fill_manual(values = color_mapping) +
  labs(
    title = "After correction",
    x = "Class",
    y = "Pearson correlation coefficient",
    fill = "Class - Group Type"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

print(p4)

# 保存为 PDF 文件
ggsave("./results/plots/pcc-boxplot-relative-isoform-expression-lr-s144.png", plot = p3, width = 7, height = 6)

# 若 cor_df_absolute 不是 data.table，可先转换
setDT(cor_df_relative)

# 计算每个 ClassGroup 的 PCC 中位数
pcc_median <- cor_df_relative[
  , .(median_PCC = median(PCC, na.rm = TRUE)),
  by = ClassGroup
][order(ClassGroup)]

print(pcc_median)

# merge plot
p3_n <- p3+theme(legend.position = "none")
p4_n <- p4+theme(legend.position = "none")


pm <- p3_n+p4_n
pm
# 保存为 PDF 文件
ggsave("./results/plots/pcc-boxplot-abs-relative-isoform-expression-lr-s144.png", plot = pm, width = 10, height = 6)
ggsave("./results/plots/pcc-boxplot-abs-relative-isoform-expression-lr-s144.pdf", plot = pm, width = 10, height = 6)

## pcc per batch
head(cor_df_absolute)
cor_df_absolute_intraB <- cor_df_absolute[cor_df_absolute$Class=="Intra-batch",]
head(cor_df_relative)
cor_df_relative_intraB <- cor_df_relative[cor_df_relative$Class=="Intra-batch",]

# Create the plot
# 1. 如果 batch_1 不是因子，则转换（并保持原有顺序）
cor_df_absolute_intraB$batch_1 <- factor(cor_df_absolute_intraB$batch_1)

# 2. 定义 GroupType 顺序（先 cross-sample，再 intra-sample）
group_order <- c("cross-sample", "intra-sample")

# 3. 创建合并变量 BatchGroup，并设定因子水平顺序：每个 batch_1 下依次排列 group_order
unique_batches <- levels(cor_df_absolute_intraB$batch_1)
desired_levels <- unlist(lapply(unique_batches, function(b) {
  paste(b, group_order, sep = " - ")
}))
cor_df_absolute_intraB$BatchGroup <- factor(
  paste(cor_df_absolute_intraB$batch_1, cor_df_absolute_intraB$GroupType, sep = " - "),
  levels = desired_levels
)

# 4. 基于 Set2 调色板为每个 batch 指定基础颜色，再根据 GroupType 调整透明度：
base_colors <- brewer.pal(n = length(unique_batches), name = "Set2")
color_mapping <- sapply(desired_levels, function(bg) {
  parts <- strsplit(bg, " - ")[[1]]
  batch <- parts[1]
  gt <- parts[2]
  base <- base_colors[which(unique_batches == batch)]
  if(gt == "cross-sample") alpha(base, 0.5) else alpha(base, 1)
})
names(color_mapping) <- desired_levels

# 5. 绘图：violins 使用合并变量 BatchGroup（生成合并图例），箱线图只用 batch_1 分组且不显示图例
p5 <- ggplot(cor_df_absolute_intraB, aes(x = batch_1, y = PCC)) +
  geom_split_violin(aes(fill = BatchGroup), trim = TRUE, position = position_dodge(width = 0.2)) +
  scale_fill_manual(values = color_mapping) +
  labs(
    title = "Before correction",
    x = "Batch",
    y = "Pearson correlation coefficient",
    fill = "Batch - Group Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )
print(p5)

# 1. 确保 batch_1 为因子（保持原有顺序）
cor_df_relative_intraB$batch_1 <- factor(cor_df_relative_intraB$batch_1)

# 2. 定义 GroupType 的顺序
group_order <- c("cross-sample", "intra-sample")

# 3. 创建合并变量 BatchGroup，并设置因子水平顺序
unique_batches <- levels(cor_df_relative_intraB$batch_1)
desired_levels <- unlist(lapply(unique_batches, function(b) {
  paste(b, group_order, sep = " - ")
}))
cor_df_relative_intraB$BatchGroup <- factor(
  paste(cor_df_relative_intraB$batch_1, cor_df_relative_intraB$GroupType, sep = " - "),
  levels = desired_levels
)

# 4. 根据 Set2 调色板生成颜色映射：同一批次中，cross-sample 用半透明，intra-sample 用不透明
base_colors <- brewer.pal(n = length(unique_batches), name = "Set2")
color_mapping <- sapply(desired_levels, function(bg) {
  parts <- strsplit(bg, " - ")[[1]]
  batch <- parts[1]
  gt <- parts[2]
  base <- base_colors[which(unique_batches == batch)]
  if (gt == "cross-sample") alpha(base, 0.5) else alpha(base, 1)
})
names(color_mapping) <- desired_levels

# 5. 绘图：violins 使用合并变量 BatchGroup（生成合并图例），箱线图按 batch_1 分组显示单个箱线
p6 <- ggplot(cor_df_relative_intraB, aes(x = batch_1, y = PCC)) +
  geom_split_violin(aes(fill = BatchGroup), trim = TRUE, position = position_dodge(width = 0.2)) +
  scale_fill_manual(values = color_mapping) +
  labs(
    title = "Ratio correction",
    x = "Batch",
    y = "Pearson correlation coefficient",
    fill = "Batch - Group Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

print(p6)

# merge plot
p5_n <- p5+theme(legend.position = "none")
p6_n <- p6+theme(legend.position = "none")


pm2 <- p5_n+p6_n
pm2

ggsave("./results/plots/pcc-boxplot-perBacth-abs-relative-isoform-expression-lr-s144.png", plot = pm2, width = 10, height = 6)
ggsave("./results/plots/pcc-boxplot-perBacth-abs-relative-isoform-expression-lr-s144.pdf", plot = pm2, width = 10, height = 6)

# remove all
rm(list=ls())
gc()

