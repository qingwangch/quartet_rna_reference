#' fig2_comprison_LR84_SR60_junctions_per_read
#' Qingwang Chen
#' 2025-07-10
#' 

# library import
# library(megadepth)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# data import
quartet_LR_junctions_s84 <- fread("/vast/projects/quartet_rna_refdata/analysis/minimap2/junctions/quartet_LR_junctions_s84.csv") %>% as.data.frame()
# quartet_LR_junctions_s48$data_type <- "long read"
quartet_SR_junctions_s60 <- fread("/vast/projects/quartet_rna_refdata/analysis/star_bam/junction/quartet_SR_junctions_s60.csv") %>% as.data.frame()
quartet_SR_junctions_s60$sample_id <- gsub("_downsampled_5000000_reads.bam","",quartet_SR_junctions_s60$sample_id )

# data merge
quartet_LS_junctions_s144 <- rbind(quartet_LR_junctions_s84,quartet_SR_junctions_s60)

# 1. 对 junctions_covered 进行分组，并计算每组的占比
# 根据 junctions_covered 分组
final_result_grouped <- quartet_LS_junctions_s144 %>%
  mutate(junction_group = case_when(
    junctions_covered == 1 ~ "1",
    junctions_covered == 2 ~ "2",
    junctions_covered == 3 ~ "3",
    junctions_covered == 4 ~ "4",
    junctions_covered == 5 ~ "5",
    junctions_covered == 6 ~ "6",
    junctions_covered >= 7 & junctions_covered <= 114 ~ "7-114",  # 这里将 7-105 作为一个组
    TRUE ~ as.character(junctions_covered)  # 其他情况使用原始值
  ))

# 对 7-105 组进行合并（按 sample_id 汇总 percentage），其他组保持原样
final_result_grouped_summarised <- final_result_grouped %>%
  group_by(sample_id, junction_group) %>%
  summarise(
    percentage = ifelse(junction_group == "7-105", sum(percentage), dplyr::first(percentage)),
    .groups = "drop"
  ) %>% distinct()

final_result_grouped_summarised <- merge(final_result_grouped_summarised,quartet_LS_junctions_s144[,c("sample_id","data_type")]) %>% distinct()

# Calculate the median percentage for each group
medians <- final_result_grouped_summarised %>%
  group_by(junction_group, data_type) %>%
  summarise(median_percentage = median(percentage))

# Create the boxplot and add lines connecting the medians
p1 <- ggplot(final_result_grouped_summarised, aes(x = junction_group, y = percentage, fill = data_type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, outlier.colour = "white") +
  geom_point(data = medians, aes(x = junction_group, y = median_percentage, group = data_type), 
             position = position_dodge(width = 0.75), size = 3, shape = 21, fill = "white") +
  geom_line(data = medians, aes(x = junction_group, y = median_percentage, group = data_type, color = data_type), 
            position = position_dodge(width = 0.75), size = 1) +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +  # 调整箱线图填充颜色
  scale_color_manual(values = c("#4E79A7", "#F28E2B")) +  # 调整线条颜色
  labs(
    x = "Number of junctions per read", 
    y = "Percentage of reads (%)", 
    title = NULL
  ) +
  theme_minimal(base_size = 14) +  # 设置整体最小主题样式
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # X轴标签旋转，大小优化
    axis.text.y = element_text(size = 12),  # Y轴标签大小
    axis.title.x = element_text(size = 14, face = "bold"),  # X轴标题字体加粗
    axis.title.y = element_text(size = 14, face = "bold"),  # Y轴标题字体加粗
    legend.position = c(0.8, 0.85),  # 将图例放在图内，右上角
    legend.title = element_blank(),  # 图例标题移除
    panel.grid = element_blank(),  # 移除背景网格线
    panel.background = element_blank(),  # 移除背景
    axis.line = element_line(size = 0.8, colour = "black"),  # 添加坐标轴线
    axis.ticks = element_line(size = 0.5, color = "black"),  # 添加轴须须
    axis.ticks.length = unit(0.2, "cm"),  # 设置轴须须的长度
    plot.margin = margin(10, 10, 10, 10)  # 调整边距
  )
p1


# 保存图像为 PNG 格式
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/junctions_per_read_LS_sr60_lr84.png", plot = p1, width = 6, height = 6, dpi = 300)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/junctions_per_read_LS_sr60_lr84.pdf", plot = p1, width = 6, height = 6, dpi = 300)

# ----------------transcripts per reads------------
read_to_transcript_lr <- readRDS("/vast/projects/quartet_rna_refdata/analysis/bambu/quartet_lr_s84_gtf/read_to_transcript.rds")
read_to_transcript_sr <- readRDS("/vast/projects/quartet_rna_refdata/analysis/bambu/quartet_sr_s60_ds/read_to_transcript.rds")

process_read_to_transcript <- function(data, data_type_label, output_file = NULL) {
  # ----------------------------
  # 计算每条 read 的覆盖转录本总数
  # ----------------------------
  data <- data %>%
    mutate(
      equalMatchCount = sapply(equalMatches, function(x) ifelse(is.null(x), 0, length(x))),
      compatibleMatchCount = sapply(compatibleMatches, function(x) ifelse(is.null(x), 0, length(x))),
      totalMatchCount = equalMatchCount + compatibleMatchCount
    )
  
  # ----------------------------
  # 统计 totalMatchCount 的频数分布
  # ----------------------------
  count_table <- table(data$totalMatchCount)
  numeric_names <- as.numeric(names(count_table))
  
  # 合并 TotalMatchCount >= 11 的频数
  merged_count_table <- c(
    count_table[numeric_names < 11],
    "11+" = sum(count_table[numeric_names >= 11])
  )
  
  # ----------------------------
  # 计算合并后频数的比例
  # ----------------------------
  merged_proportion_table <- prop.table(merged_count_table)
  
  # ----------------------------
  # 转换结果为数据框格式
  # ----------------------------
  result_df <- data.frame(
    TotalMatchCount = names(merged_count_table),
    Frequency = as.vector(merged_count_table),
    Proportion = as.vector(merged_proportion_table)
  )
  
  # 设置 TotalMatchCount 为因子变量，并控制水平顺序
  result_df$TotalMatchCount <- factor(
    result_df$TotalMatchCount,
    levels = c(as.character(0:10), "11+"),
    ordered = TRUE
  )
  
  # 添加数据类型标签
  result_df$data_type <- data_type_label
  
  # 返回处理结果
  return(result_df)
}

# Calculate the median percentage for each group
result_df_lr <- process_read_to_transcript(read_to_transcript_lr,"long read")
result_df_sr <- process_read_to_transcript(read_to_transcript_sr,"short read")

# merged
result_df <- rbind(result_df_lr,result_df_sr)
medians <- result_df %>%
  group_by(TotalMatchCount, data_type) %>%
  summarise(median_proportion = median(Proportion))


# plot
p2 <- ggplot(result_df, aes(x = TotalMatchCount, y = Proportion, fill = data_type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, outlier.colour = "white") +
  geom_point(data = medians, aes(x = TotalMatchCount, y = median_proportion, group = data_type), 
             position = position_dodge(width = 0.75), size = 3, shape = 21, fill = "white") +
  geom_line(data = medians, aes(x = TotalMatchCount, y = median_proportion, group = data_type, color = data_type), 
            position = position_dodge(width = 0.75), size = 1) +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +  # 调整箱线图填充颜色
  scale_color_manual(values = c("#4E79A7", "#F28E2B")) +  # 调整线条颜色
  labs(
    x = "Transcript assignments per read", 
    y = "Proportion of reads (%)", 
    title = NULL
  ) +
  theme_minimal(base_size = 14) +  # 设置整体最小主题样式
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # X轴标签旋转，大小优化
    axis.text.y = element_text(size = 12),  # Y轴标签大小
    axis.title.x = element_text(size = 14, face = "bold"),  # X轴标题字体加粗
    axis.title.y = element_text(size = 14, face = "bold"),  # Y轴标题字体加粗
    legend.position = c(0.8, 0.85),  # 将图例放在图内，右上角
    legend.title = element_blank(),  # 图例标题移除
    panel.grid = element_blank(),  # 移除背景网格线
    panel.background = element_blank(),  # 移除背景
    axis.line = element_line(size = 0.8, colour = "black"),  # 添加坐标轴线
    axis.ticks = element_line(size = 0.5, color = "black"),  # 添加轴须须
    plot.margin = margin(10, 10, 10, 10)  # 调整边距
  )
p2

# 保存图像为 PNG 格式
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/transcript_assignments_per_read_LR84_SR60.png", plot = p2, width = 6, height = 6, dpi = 300)

# combined_plot <- plot_grid(
#   p1 ,  
#   p2 ,  
#   align = "h",  
#   ncol = 2,     # 单列布局
#   labels = c("")  # 添加子图标签
# )
library(cowplot)
combined_plot <- plot_grid(
  p1 ,  
  p2 ,  
  p3 ,
  align = "h",  
  ncol = 3,     # 单列布局
  labels = c("")  # 添加子图标签
)
combined_plot

ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/junctions_transcript_assignments_per_read_overdispersion_LR84_SR60.pdf", plot = combined_plot, width = 13, height = 5)


# remove all
rm(list=ls())
gc()
