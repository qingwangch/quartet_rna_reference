#' junctions_per_read
#' Qingwang Chen
#' 2025-07-09

# setwd
getwd()
setwd("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/minimap2")

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())


# library import
library(megadepth)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# define the function
process_bam_file <- function(test_bam) {
  # 运行 bam_to_junctions()
  example_jxs <- bam_to_junctions(test_bam, all_junctions = TRUE, overwrite = TRUE, long_reads = TRUE)
  
  # 读取数据
  example_jxs <- read_junction_table(example_jxs)
  
  # 使用 data.table 来统计每个 read_name 的 junction 数量
  # .N 是一个特殊的符号，表示当前组的大小
  example_jxs_dt <- as.data.table(example_jxs)
  junction_count <- example_jxs_dt[, .(junctions_covered = sum(unique > 0)), by = read_name]
  
  # 计算每个 read_name 对应的 junctions_covered 占比
  junction_count_percentage <- junction_count[, .(n = .N), by = junctions_covered]
  junction_count_percentage[, percentage := n / sum(n) * 100]
  
  # 提取路径中的文件夹名称
  sample_id <- str_split(basename(test_bam), "\\.")[[1]][1]
  
  # 添加样本 ID 列
  junction_count_percentage$sample_id <- sample_id
  
  return(junction_count_percentage)
}

# 假设有多个样本路径，按批次处理
# sample_bam_paths <- list.files("/vast/projects/quartet_rna_refdata/analysis/minimap2/", 
#                                pattern = "*.sorted.only_primary.bam$", recursive = T)
sample_bam_paths <- list.files("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/minimap2/", 
                               pattern = "*.sorted.only_primary.bam$", recursive = T)

# 只选择四个批次
sample_bam_paths <- sample_bam_paths[grep("D_ONT_LG_B1|D_ONT_LW_B1|M_PAB_LN_B1|M_PAB_LG_B2|P_ONT_LG_B1|P_ONT_LG_B2|P_ONT_LN_B2",sample_bam_paths)]
sample_bam_paths <- sample_bam_paths[-grep("HCC1395|SIRV|D_ONT_LW_B1_D6_1.minimap2.sorted.only_primary.bam",sample_bam_paths)]

# 分批处理，每次处理 10 个文件
batch_size <- 10
n_batches <- ceiling(length(sample_bam_paths) / batch_size)

all_results <- list()

for (i in 1:n_batches) {
  # 计算当前批次的文件范围
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(sample_bam_paths))
  
  # 获取当前批次的文件
  batch_files <- sample_bam_paths[start_idx:end_idx]
  print(batch_files)
  
  # 处理当前批次的文件
  batch_results <- lapply(batch_files, process_bam_file)
  
  # 合并当前批次的结果
  all_results[[i]] <- rbindlist(batch_results)
}

# 合并所有批次的结果
final_result <- rbindlist(all_results)
final_result$data_type <- "long read"

# 查看最终结果
head(final_result)

# save
fwrite(final_result,"/vast/projects/quartet_rna_refdata/analysis/minimap2/junctions/quartet_LR_junctions_s84.csv")

# 1. 对 junctions_covered 进行分组，并计算每组的占比
# 根据 junctions_covered 分组
final_result_grouped <- final_result %>%
  mutate(junction_group = case_when(
    junctions_covered == 1 ~ "1",
    junctions_covered == 2 ~ "2",
    junctions_covered == 3 ~ "3",
    junctions_covered == 4 ~ "4",
    junctions_covered == 5 ~ "5",
    junctions_covered == 6 ~ "6",
    junctions_covered >= 7 & junctions_covered <= 105 ~ "7-105",  # 这里将 7-105 作为一个组
    TRUE ~ as.character(junctions_covered)  # 其他情况使用原始值
  ))

# 对 7-105 组进行合并（按 sample_id 汇总 percentage），其他组保持原样
final_result_grouped_summarised <- final_result_grouped %>%
  group_by(sample_id, junction_group) %>%
  summarise(
    percentage = ifelse(junction_group == "7-105", sum(percentage), first(percentage)),
    .groups = "drop"
  ) %>% distinct()

# 绘制箱型图
ggplot(final_result_grouped_summarised, aes(x = junction_group, y = percentage, fill = junction_group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, outlier.colour = "#FF5733", # 设置离群点的颜色、形状和大小
  ) +  # 设置箱子的宽度
  scale_fill_manual(values = rep("#bdd3fd", length(unique(final_result_grouped_summarised$junction_group)))) + # 使用统一的颜色填充
  theme_minimal(base_size = 16) +  # 基础字体大小
  labs(
    x = "Number of Junctions Covered", 
    y = "Percentage (%)", 
    title = "Percentage of Junctions Covered by Reads"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 旋转x轴标签并设置大小
    axis.text.y = element_text(size = 12),  # 设置y轴标签大小
    axis.title = element_text(size = 14),  # 设置轴标题的字体大小
    plot.title = element_text(size = 16, hjust = 0.5),  # 设置标题的字体大小并居中
    legend.position = "none",  # 不显示图例
    panel.grid.major = element_line(color = "grey90", size = 0.5),  # 设置主网格线的颜色和大小
    panel.grid.minor = element_line(color = "grey95", size = 0.25),  # 设置次网格线的颜色和大小
    panel.border = element_blank()  # 去除边框
  ) 


# remove all
rm(list=ls())
gc()

