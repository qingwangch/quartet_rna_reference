#' overdispersion-short-read-from-salmon-and-long-read-from-oarfish-sr60-lr84
#' Qingwang Chen
#' 2025-04-02
#' 

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())


# library import
library(edgeR)
source("/vast/projects/quartet_rna_refdata/analysis/R/figures/functions/catchOarfish.R")
library(arrow)
library(ggplot2)
library(reshape2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(tidyr) 

## salmon
salmon_quants <- dir("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/salmon_quant/")[grep("P_ILM_L8_B1|R_ILM_L8_B1|P_BGI_L3_B1|R_BGI_L3_B1|R_ELE_LH_B1", dir("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/salmon_quant/"))]
#[grep("P_", dir("./data/salmon_quant/"))]
salmon_quants <- paste0("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/salmon_quant/",salmon_quants)
catch <- catchSalmon(salmon_quants,verbose=T)
# scaled.counts <- catch$counts/catch$annotation$Overdispersion
annn_salmon <- catch$annotation
count_salmon <- catch$counts
colnames(count_salmon) <- gsub("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/salmon_quant/","",colnames(count_salmon))
head(count_salmon)
annn_salmon$source <- "salmon"
head(annn_salmon)
annn_salmon$Transcript <- rownames(annn_salmon)
# Save both data frames into a single RDS file
saveRDS(list(annn_salmon = annn_salmon, count_salmon = count_salmon), file = "/vast/projects/quartet_rna_refdata/analysis/R/figures/Rdata/salmon_quartet_s60.rds")

## oarfish
oarfish_quants <- dir("/vast/projects/quartet_rna_refdata/analysis/oarfish_quant")[grep("P_ONT_LN_B2|P_ONT_LG_B1|P_ONT_LG_B2|D_ONT_LG_B1|D_ONT_LW_B1|M_PAB_LN_B1|M_PAB_LG_B2", dir("/vast/projects/quartet_rna_refdata/analysis/oarfish_quant"))]
oarfish_quants <- oarfish_quants[-grep("HCC|D_ONT_LW_B1_D6_1$",oarfish_quants)]
oarfish_quants <- paste0("/vast/projects/quartet_rna_refdata/analysis/oarfish_quant/",oarfish_quants)
s <- catchOarfish(oarfish_quants,verbose=T)
annn_oarfish <- s$annotation
count_oarfish <- s$counts
colnames(count_oarfish) <- gsub("/vast/projects/quartet_rna_refdata/analysis/oarfish_quant/","",colnames(count_oarfish))
annn_oarfish$source <- "oarfish"
# Save both data frames into a single RDS file
saveRDS(list(annn_oarfish = annn_oarfish, count_oarfish = count_oarfish), file = "/vast/projects/quartet_rna_refdata/analysis/R/figures/Rdata/oarfish_quartet_s84.rds")


annn_oarfish$Transcript <- sapply(rownames(annn_oarfish), function(x) strsplit(x, "\\|")[[1]][1])

# 拆分注释信息
annotation_df <- do.call(rbind, strsplit(rownames(annn_oarfish), "\\|")) %>% as.data.frame(stringsAsFactors = FALSE)
colnames(annotation_df) <- c(
  "Transcript_ID", "Gene_ID", "Additional1", "Additional2",
  "Transcript_Name", "Gene_Name", "Length", "Gene_Type"
)
geneid2txid <- annotation_df[,c("Transcript_ID","Gene_ID")]

# 检查是否有重复的 Transcript
duplicates <- annn_oarfish$Transcript[duplicated(annn_oarfish$Transcript)]

# 查看重复值
if (length(duplicates) > 0) {
  print("重复的 Transcript:")
  print(duplicates)
} else {
  print("没有重复的 Transcript。")
}

# 合并annotation
annn_merge <- merge(annn_salmon,annn_oarfish,by="Transcript")

annn_merge <- merge(geneid2txid,annn_merge,by.x="Transcript_ID",by.y = "Transcript")

head(annn_merge)

# 统计每个基因的转录本数量，并计算不同来源的 Ambiguity 平均值
gene_transcript_count <- annn_merge %>%
  group_by(Gene_ID) %>%
  summarise(
    Transcript_Count = n(),  # 统计每个基因的转录本数量
    Ambiguity_SR = mean(Overdispersion.x[source.x == "salmon"], na.rm = TRUE),  # SR 的平均 Overdispersion
    Ambiguity_LR = mean(Overdispersion.y[source.y == "oarfish"], na.rm = TRUE)  # LR 的平均 Overdispersion
  )

# 将数据从宽格式转换为长格式
gene_transcript_count_long <- gene_transcript_count %>%
  pivot_longer(
    cols = starts_with("Ambiguity"),
    names_to = "Source",
    values_to = "Ambiguity"
  ) %>%
  mutate(
    Source = case_when(
      Source == "Ambiguity_SR" ~ "SR",
      Source == "Ambiguity_LR" ~ "LR"
    )
  )

# 添加类别分组（用于模仿图中 [2, 6], [8, 11] 等分组）
gene_transcript_count_long <- gene_transcript_count_long %>%
  mutate(
    Group = cut(
      Transcript_Count,
      breaks = c(0, 1, 5, 10, 20, 50, 296), 
      labels = c("1", "[2, 6)", "[6, 11)", "[11, 21)", "[21, 51)", "[51, 296]")
    )
  )

# 计算每个组的总转录本数量（统计的是转录本数，而非基因数）
group_counts <- gene_transcript_count_long %>%
  filter(Source == "SR") %>%  # 仅计算 SR 作为代表
  group_by(Group) %>%
  summarise(Transcript_Count = sum(Transcript_Count), .groups = "drop")  # 统计每个 Group 内的转录本总数

# 绘图
gene_transcript_count_long$Source <- factor(gene_transcript_count_long$Source,levels = c("SR","LR"))
p3 <- ggplot(gene_transcript_count_long, aes(x = Group, y = Ambiguity, fill = Source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.6) +  # 设置箱线图宽度
  geom_text(
    data = group_counts, 
    aes(x = Group, y = max(gene_transcript_count_long$Ambiguity, na.rm = TRUE) * 1.1, 
        label = Transcript_Count),  # 显示转录本数量
    inherit.aes = FALSE,  # 不继承全局 aes
    size = 4,
    color = "#F28E2B"
  ) +  # 每个类别只显示一个转录本数量
  scale_fill_manual(values = c("SR" = "#F28E2B", "LR" = "#4E79A7")) +  # 手动设置箱线图填充颜色
  labs(
    x = "Number of transcripts per gene", 
    y = "Assignment ambiguity",
    title = NULL
  ) +
  scale_y_log10(limits = c(1, NA)) +  # 对数刻度，并设置 y 轴从 1 开始
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # X轴标签样式
    axis.text.y = element_text(size = 12),  # Y轴标签样式
    axis.title.x = element_text(size = 14, face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(size = 14, face = "bold"),  # Y轴标题加粗
    legend.position = c(0.8, 0.85),  # 图例位置
    legend.title = element_blank(),  # 图例标题移除
    panel.grid = element_blank(),  # 移除背景网格线
    panel.background = element_blank(),  # 移除背景色
    axis.line = element_line(size = 0.8, colour = "black"),  # 添加坐标轴线
    axis.ticks = element_line(size = 0.5, color = "black"),  # 添加轴须须
    axis.ticks.length = unit(0.2, "cm"),  # 设置轴须须长度
    plot.margin = margin(10, 10, 10, 10)  # 调整边距
  ) +
  annotate(
    "text", x = 3, y = max(gene_transcript_count_long$Ambiguity, na.rm = TRUE) * 2, 
    label = "",
    hjust = 0.5, size = 4, color = "black", fontface = "italic"
  )  # 添加说明

# 显示图形
print(p3)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/overdispersion_sr60_lr84_byNum.pdf", plot = p3, width = 5, height = 5, dpi = 300)

### by length
annn_merge_plot <- log2(annn_merge[,c("Overdispersion.x","Overdispersion.y")])

# 进一步画一下根据Length分组的Overdispersion
## 去掉Overdispersion超过100的行
annn_merge <- annn_merge[annn_merge$Overdispersion.x<30,]

# 按 Length 分组
annn_merge$Length_Group <- cut(
  annn_merge$Length.x, 
  breaks = c(0, 1000, 2000, 3000, 4000, 5000, Inf), 
  labels = c("0-1000", "1001-2000", "2001-3000", "3001-4000", "4001-5000", ">5000")
)

# 转换数据为长格式
annn_merge_long <- annn_merge %>%
  pivot_longer(
    cols = c(Overdispersion.x, Overdispersion.y), 
    names_to = "Source", 
    values_to = "Overdispersion"
  )

# 绘制boxplot
p2 <- ggplot(annn_merge_long, aes(x = Length_Group, y = Overdispersion, fill = Source)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏离群点
  scale_fill_manual(
    values = c("Overdispersion.x" = "#F28E2B", "Overdispersion.y" = "#4E79A7"),
    labels = c("Salmon (SR)", "Oarfish (LR)")
  ) +
  labs(
    x = "Transcript Length Group",
    y = "Overdispersion",
    fill = "Source",
    title = "Overdispersion by Transcript Length Group"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.position = "top",
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    axis.line = element_line(color = "black", size = 0.8)  # 显示 x 和 y 轴线
  )
p2

ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/overdispersion_LS_sr60_lr84_byLength.png", plot = p2, width = 6, height = 6, dpi = 300)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/overdispersion_LS_sr60_lr84_byLength.pdf", plot = p2, width = 5, height = 5, dpi = 300)


## raw data bind
count_salmon <- as.data.frame(count_salmon)
count_oarfish <- as.data.frame(count_oarfish)

# count_salmon$Transcript <- rownames(count_salmon)
# count_oarfish$Transcript <- sapply(rownames(count_oarfish), function(x) strsplit(x, "\\|")[[1]][1])

count_merge <- merge(count_salmon,count_salmon,by="Transcript")
rownames(count_merge) <- count_merge$Transcript
count_merge$Transcript <- NULL

saveRDS(count_merge, file = "/vast/projects/quartet_rna_refdata/analysis/R/figures/Rdata/salmon_oarfish_quartet_s120.rds")

# remove all
rm(list = ls())
gc()
