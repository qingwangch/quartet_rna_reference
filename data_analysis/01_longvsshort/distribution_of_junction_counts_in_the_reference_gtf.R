#' distribution_of_junction_counts_in_the_reference_gtf
#' Qingwang Chen
#' 2025-04-02
#' Local

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())

# Load libraries
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(rtracklayer)
library(ggpubr)
library(data.table)

# 读取本地GTF文件
gtf_file <- "./data/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"  # 替换为实际GTF文件路径
gtf_data <- rtracklayer::import(gtf_file)

# 提取转录本和 exon 信息
gtf_exons <- as.data.frame(gtf_data) %>%
  filter(type == "exon") %>%  # 仅保留 exon 信息
  dplyr::select(transcript_id, exon_id, start, end)  # 提取需要的列

# Step 2: 统计每个转录本的 exon 数量
transcript_exon_counts <- gtf_exons %>%
  group_by(transcript_id) %>%
  summarise(
    exon_count = n(),  # 统计 exon 的数量
    junction_count = exon_count - 1  # 计算 junction 数量
  )

# Step 3: 统计不同 junction 数目的转录本数量及比例
junction_distribution <- transcript_exon_counts %>%
  group_by(junction_count) %>%
  summarise(
    transcript_count = n(),  # 每种 junction 数对应的转录本数量
    percentage = (transcript_count / nrow(transcript_exon_counts)) * 100  # 计算比例
  ) %>%
  arrange(junction_count)  # 按 junction 数目升序排列

# 查看结果
print(junction_distribution)

# 合并 junction_count >= 7 的数据
junction_distribution <- junction_distribution %>%
  mutate(
    junction_group = ifelse(junction_count >= 7, ">=7", as.character(junction_count))
  ) %>%
  group_by(junction_group) %>%
  summarise(
    transcript_count = sum(transcript_count),
    percentage = sum(percentage),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(gsub(">=7", "7", junction_group)))  # 确保排序正确



# 可视化
junction_distribution$junction_group <- factor(junction_distribution$junction_group,
                                               levels = c("0","1","2","3","4","5","6",">=7"))
p1 <- ggplot(junction_distribution, aes(x = junction_group, y = percentage)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = round(percentage, 2)), vjust = -0.5, size = 4) +
  labs(
    title = "Distribution of Junction Counts \nReference Transcriptome",
    x = "The Number of Junctions",
    y = "Percentage of Transcripts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),  # 加粗 X 轴字体
    axis.text.y = element_text(face = "bold", size = 12),  # 加粗 Y 轴字体
    axis.title.x = element_text(face = "bold", size = 14),  # 加粗 X 轴标题
    axis.title.y = element_text(face = "bold", size = 14),  # 加粗 Y 轴标题
    axis.line = element_line(color = "black"),  # 添加 X 和 Y 轴线
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # 标题加粗
  )
p1
# Save the plot in high resolution for publication
ggsave("./results/plot/distribution_of_junction_counts_in_the_reference_transcriptome_bars.png", plot = p1, width = 6, height = 6, dpi = 300)

# remove all
rm(list=ls())
gc()
