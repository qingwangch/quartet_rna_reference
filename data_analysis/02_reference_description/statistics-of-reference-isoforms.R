#' statistics-of-reference-isoforms
#' Qingwang Chen
#' 2025-04-09
#' 

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())


# library import
library(edgeR)
library(arrow)
library(ggplot2)
library(reshape2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

# data import
## LO
############ refFC
refFC_final_p4_lo <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b6_p4_s72_u_20250403.csv")
table(refFC_final_p4_lo$compare)
############ refDEI
refDEI_final_p4_lo <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250407.csv")
refDEI_final_p4_lo <- refDEI_final_p4_lo %>% filter(refDEI_final_p4_lo$Final != "non-DEI")
table(refDEI_final_p4_lo$compare)
## SO
############ refFC
refFC_final_p4_so <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/ref_expr_b13_p4_s156_u_20250407.csv")
table(refFC_final_p4_so$compare)
############ refDEI
refDEI_final_p4_so <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/RefData_DEIs_all_isoforms_classified_u_20250407.csv")
refDEI_final_p4_so <- refDEI_final_p4_so %>% filter(refDEI_final_p4_so$Final != "non-DEI")
table(refDEI_final_p4_so$compare)
## LS
############ refFC
refFC_final_p4_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/ref_expr_LS_t2_u_20250407.csv")
table(refFC_final_p4_ls$compare)
############ refDEI
refDEI_final_p4_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/RefData_DEIs_LS_t2_u_20250407.csv")
refDEI_final_p4_ls <- refDEI_final_p4_ls %>% filter(refDEI_final_p4_ls$Final != "non-DEI")
table(refDEI_final_p4_ls$compare)

## data process
# 函数：读取并统计 compare 数量
read_and_count <- function(file, type, dataset) {
  df <- read.csv(file)
  
  # 仅对 refDEI 进行 Final != "non-DEI" 过滤
  if (type == "refDEI") {
    df <- df %>% filter(Final != "non-DEI")
  }
  
  df %>%
    count(compare, name = "count") %>%
    mutate(type = type, dataset = dataset)
}

# 加载各数据表
fc_lo  <- read_and_count("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b6_p4_s72_u_20250403.csv", "refFC", "LO")
dei_lo <- read_and_count("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250407.csv", "refDEI", "LO")

fc_so  <- read_and_count("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/ref_expr_b13_p4_s156_u_20250407.csv", "refFC", "SO")
dei_so <- read_and_count("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/RefData_DEIs_all_isoforms_classified_u_20250407.csv", "refDEI", "SO")

fc_ls  <- read_and_count("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/ref_expr_LS_t2_u_20250407.csv", "refFC", "LS")
dei_ls <- read_and_count("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/RefData_DEIs_LS_t2_u_20250407.csv", "refDEI", "LS")

# 合并所有表格
all_data <- bind_rows(fc_lo, dei_lo, fc_so, dei_so, fc_ls, dei_ls)

# boxplot (fig 4b)
# 设置 factor 顺序（确保横轴有序）
all_data$compare <- factor(all_data$compare, levels = c("D5/D6", "F7/D6", "M8/D6"))
all_data$dataset <- factor(all_data$dataset, levels = c("LO", "SO", "LS"))
all_data$type <- factor(all_data$type, levels = c("refFC", "refDEI"))

# 自定义颜色（NBT常用色系：蓝-橙-绿）
# 假设 all_data 里包含 count, dataset, compare, type
# 计算缩放因子：根据 max 值之间的比例
scale_factor <- max(all_data$count[all_data$type == "refFC"]) / max(all_data$count[all_data$type == "refDEI"])

# 创建新列，refDEI 数据进行缩放
all_data_scaled <- all_data %>%
  mutate(scaled_count = ifelse(type == "refDEI", count * scale_factor, count))

# 绘图
p1 <- ggplot(all_data_scaled, aes(x = dataset, y = scaled_count, fill = compare)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~type, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c("D5/D6" = "#1b9e77", "F7/D6" = "#d95f02", "M8/D6" = "#7570b3")) +
  
  # 添加次坐标轴
  scale_y_continuous(
    name = "refFC Isoform Count",
    sec.axis = sec_axis(~ . / scale_factor, name = "refDEI Isoform Count")
  ) +
  
  labs(
    x = "Datasets",
    fill = "Compare Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 14),
    axis.title.y.left = element_text(color = "#333333", size = 13, face = "bold"),
    axis.title.y.right = element_text(color = "#555555", size = 13, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA)
  )

# 打印图
print(p1)

# 保存高分辨率图像
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/compare_isoform_count_quadrant_by_dataset.png", plot = p1, width = 5, height = 5, dpi = 300)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/compare_isoform_count_quadrant_by_dataset.pdf", plot = p1, width = 5, height = 5, dpi = 300)

# scatter plot
## LO vs SO
# 根据 isoform_compare 合并两个数据框
merged_data <- merge(
  refFC_final_p4_so,
  refFC_final_p4_lo,
  by = "isoform_compare",
  suffixes = c("_so", "_lo")
)

# 检查合并后的数据框
head(merged_data)
# 计算相关系数
correlation <- cor(merged_data$log2FC_so, merged_data$log2FC_lo, use = "complete.obs")
cat("Correlation coefficient (log2FC):", correlation, "\n")

# 绘制散点图
# 自定义颜色映射，与前一图一致
nbt_colors <- c(
  "D5/D6" = "#1b9e77",
  "F7/D6" = "#d95f02",
  "M8/D6" = "#7570b3"
)

p2 <- ggplot(merged_data, aes(x = log2FC_so, y = log2FC_lo, color = compare_so)) +
  geom_point(alpha = 0.4, size = 2) +  # 更小更透明的点
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  scale_color_manual(values = nbt_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    subtitle = paste0("N = ", nrow(merged_data)),
    x = "log2FC (SO reference)",
    y = "log2FC (LO reference)",
    color = "Comparison"
  ) +
  annotate("text", x = min(merged_data$log2FC_so, na.rm = TRUE),
           y = max(merged_data$log2FC_lo, na.rm = TRUE),
           label = paste0("r = ", round(correlation, 2)),
           hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p2)

# 保存图像为 PNG 文件
ggsave(filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/Correlation_of_log2FC_between_SO_and_LO_reference_datasets.png", plot = p2, width = 6, height = 5, dpi = 300)
ggsave(filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/Correlation_of_log2FC_between_SO_and_LO_reference_datasets.pdf", plot = p2, width = 6, height = 5, dpi = 300)

## LS vs qPCR
qpcr_results <- fread("/vast/projects/quartet_rna_refdata/analysis/isoforms/qpcr/log2fc_20250409_U.csv")

# 根据 isoform_compare 合并两个数据框
merged_data <- merge(
  qpcr_results,
  refFC_final_p4_ls,
  by = "isoform_compare",
  suffixes = c("_pcr", "_ref")
)

# 检查合并后的数据框
head(merged_data)
# 计算相关系数
correlation <- cor(merged_data$log2FC_pcr, merged_data$log2FC_ref, use = "complete.obs")
cat("Correlation coefficient (log2FC):", correlation, "\n")

# 绘制散点图
table(merged_data$source)
p3 <- ggplot(merged_data, aes(x = log2FC_ref, y = log2FC_pcr, color = compare_pcr)) +
  geom_point(alpha = 0.4, size = 2) +  # 更小更透明的点
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  scale_color_manual(values = nbt_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    subtitle = paste0("N = ", nrow(merged_data)),
    x = "log2FC (LS reference)",
    y = "log2FC (qPCR)",
    color = "Comparison"
  ) +
  annotate("text", x = min(merged_data$log2FC_pcr, na.rm = TRUE),
           y = max(merged_data$log2FC_ref, na.rm = TRUE),
           label = paste0("r = ", round(correlation, 2)),
           hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p3)

# 保存图像为 PNG 文件
ggsave(filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/Correlation_of_log2FC_between_qPCR_and_LS_reference_datasets.png", plot = p3, width = 6, height = 5, dpi = 300)
ggsave(filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/Correlation_of_log2FC_between_qPCR_and_LS_reference_datasets.pdf", plot = p3, width = 6, height = 5, dpi = 300)

## LS vs dRNA
# will be updated
p4 <- NULL

# merge figure
p1_n <- p1+theme(legend.position = "none")
p2_n <- p2+theme(legend.position = "none")
p3_n <- p3+theme(legend.position = "none")

combined_plot <- plot_grid(
  plot_grid(
    p1_n, p2_n, p3_n, p4,
    align = "h", ncol = 2, labels = c("b", "c","d","e")  # 横向排列子图，并添加标签 A 和 B
  )
)
combined_plot

# 保存图像为 PNG 文件
ggsave(filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/fig4bcde.png", plot = combined_plot, width = 7, height = 8, dpi = 300)
ggsave(filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/fig4bcde.pdf", plot = combined_plot, width = 7, height = 8, dpi = 300)


## remove all
rm(list=ls())
gc()


