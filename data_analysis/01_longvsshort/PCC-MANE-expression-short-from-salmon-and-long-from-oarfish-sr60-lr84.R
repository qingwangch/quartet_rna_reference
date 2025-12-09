#' PCC-MANE-expression-short-from-salmon-and-long-from-oarfish-sr60-lr84
#' Qingwang Chen
#' 2025-06-12
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
library(tibble)  
library(edgeR)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

library(arrow)
library(ggpubr)
library(reshape2)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
# library(ggVennDiagram)
# library(tidysdm)
library(cowplot)

# data import
# count_merge <- readRDS("/vast/projects/quartet_rna_refdata/analysis/R/figures/Rdata/salmon_oarfish_quartet_s120.rds")
count_lr <- readRDS("./isoforms/Rdata/quartet-LO-minimap2-oarfish-g-txq-count-s144t252913.rds")
rownames(count_lr) <- sub("\\|.*", "", rownames(count_lr))
count_sr <- readRDS("/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/isoform_level_count_salmon_G_b22_s264.rds")

count_merge <- full_join(
  rownames_to_column(as.data.frame(count_sr), "tx"),
  rownames_to_column(as.data.frame(count_lr), "tx"),
  by = "tx"
) %>% column_to_rownames("tx")
count_merge[is.na(count_merge)] <- 0        # NA → 0

# Process metadata
process_metadata <- function(expr_mat) {
  metadata <- data.frame(
    sample_id = colnames(expr_mat),
    lib = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[1]),
    platform = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[2]),
    lab = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[3]),
    Batch = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[4]),
    sample = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[5]),
    rep = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[6])
  )
  metadata$group <- ifelse(grepl("D5|D6|F7|M8", metadata$sample_id), "Quartet", ifelse(grepl("HCC1395", metadata$sample_id), "HCC1395", "HCC1395BL"))
  metadata$batch_id <- paste(metadata$lib, metadata$platform, metadata$lab, metadata$Batch, sep = "_")
  metadata$Batch <- metadata$batch_id
  metadata$code <- metadata$sample_id
  # write.csv(metadata, "./results/tables/metadata4LS_LG_LN_s48.csv", quote = FALSE)
  return(metadata)
}

metadata <- process_metadata(count_merge)
metadata$tech <- ifelse(metadata$platform %in% c("BGI","ILM","ELE"), "short", "long")

# process batch
batch_of_compare <- c("P_BGI_L3_B1","R_BGI_L3_B1","P_ILM_L8_B1","R_ILM_L8_B1","R_ELE_LH_B1",
                      "D_ONT_LG_B1","D_ONT_LW_B1","M_PAB_LN_B1","M_PAB_LG_B2","P_ONT_LN_B2","P_ONT_LG_B1","P_ONT_LG_B2")
keep_samples <- metadata$sample_id[metadata$batch_id %in% batch_of_compare]
count_merge  <- count_merge[ , keep_samples]
metadata     <- metadata   %>% filter(sample_id %in% keep_samples)

# process isoform
MANE_isoforms <- fread("/vast/projects/quartet_rna_refdata/analysis/MANE/MANE.GRCh38.v1.4.summary.txt") %>% as.data.frame()
mane_isoforms <- MANE_isoforms$Ensembl_nuc
mani_keep    <- intersect(MANE_isoforms$Ensembl_nuc, rownames(count_merge))
count_merge  <- count_merge[mani_keep, ]
## 选定 MANE isoform 行后保存 ----
out_path <- "/vast/projects/quartet_rna_refdata/analysis/R/figures/Rdata/count_merge_MANE_filtered_s144_lr84sr60.rds"

saveRDS(count_merge, out_path)
cat("✅ Filtered count matrix saved to:\n   ", out_path, "\n")

### jaccard index
# ── 二值化表达：1 = expressed, 0 = not expressed ────────────────
expr_bin <- (count_merge > 0) * 1    
binary <- as.matrix(expr_bin) > 0
n <- ncol(binary)
jaccard_mat <- matrix(0, n, n, dimnames = list(colnames(binary), colnames(binary)))

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    inter <- sum(binary[,i] & binary[,j])
    union <- sum(binary[,i] | binary[,j])
    jacc  <- inter / union
    jaccard_mat[i, j] <- jaccard_mat[j, i] <- jacc
  }
  jaccard_mat[i, i] <- 1
}
jaccard_mat[n, n] <- 1

library(reshape2)
jacc_df <- melt(jaccard_mat,
                varnames = c("sample1", "sample2"),
                value.name = "Jaccard") %>%
  filter(as.character(sample1) < as.character(sample2))

# 添加 tech / category
tech_vec    <- setNames(metadata$tech, metadata$sample_id)
sample_vec  <- setNames(metadata$sample, metadata$sample_id)

jacc_df <- jacc_df %>%
  mutate(
    tech1 = tech_vec[sample1],
    tech2 = tech_vec[sample2],
    category = case_when(
      tech1 == "long"  & tech2 == "long"  ~ "long-vs-long",
      tech1 == "short" & tech2 == "short" ~ "short-vs-short",
      TRUE                                ~ "long-vs-short"
    ),
    sample_cat = if_else(sample_vec[sample1] == sample_vec[sample2],
                         paste0("intra-",sample_vec[sample1]),
                         "inter-sample")
  )

## ────────────────────────────────────────────────────────────────────
## 1  因子顺序
## ────────────────────────────────────────────────────────────────────
jacc_df$category <- factor(jacc_df$category,
                           levels = c("long-vs-long",
                                      "short-vs-short",
                                      "long-vs-short"))

## ────────────────────────────────────────────────────────────────────
## 2  样本分组颜色 & 透明度
## ────────────────────────────────────────────────────────────────────
fill_cols <- c("intra-D5"      = "#4CC3D9",
               "intra-D6"      = "#7BC8A4",
               "intra-F7"      = "#FFC65D",
               "intra-M8"      = "#F16745",
               "inter-sample"  = "#A58CD9")

alpha_cols <- c("long-vs-long"   = 1.00,
                "short-vs-short" = 0.50,
                "long-vs-short"  = 0.20)

## ────────────────────────────────────────────────────────────────────
## 3  每 facet×category 的 n
## ────────────────────────────────────────────────────────────────────
n_df <- jacc_df %>%
  group_by(sample_cat, category) %>%
  summarise(n   = n(),
            y_n = max(Jaccard) + 0.02,
            .groups = "drop")

## ────────────────────────────────────────────────────────────────────
## 4  仅比较 “short-vs-short” vs “long-vs-long”
## ────────────────────────────────────────────────────────────────────
j_sub <- jacc_df %>%
  filter(category %in% c("short-vs-short", "long-vs-long"))

p_tbl <- compare_means(
  Jaccard ~ category,
  group.by = "sample_cat",
  data     = j_sub,
  method   = "wilcox.test",
  p.adjust.method   = "BH"
) %>%
  mutate(
    category   = "short-vs-short",
    y.position = max(j_sub$Jaccard) + 0.03,
    xmin       = "short-vs-short",
    xmax       = "long-vs-long"
  ) %>%
  mutate(
    p.adj.signif = case_when(          # ★ 用校正后 p 值重新打星
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj <= 0.05   ~ "*",
      TRUE            ~ "ns"
    )
  )


## ────────────────────────────────────────────────────────────────────
## 5  绘图
## ────────────────────────────────────────────────────────────────────
p_jacc <- ggplot(jacc_df,
            aes(x = category, y = Jaccard,
                fill  = sample_cat,
                alpha = category)) +
  
  geom_boxplot(width = .48, colour = "grey25", outlier.shape = NA) +
  
  geom_text(data = n_df,
            aes(x = category, y = y_n,
                label = paste0("n=", n)),
            vjust = 0, size = 3.3, fontface = "italic",
            inherit.aes = FALSE) +
  
  facet_wrap(~ sample_cat, nrow = 1, scales = "free_x") +
  
  scale_fill_manual(values = fill_cols, name = "Sample group") +
  scale_alpha_manual(values = alpha_cols, name = "Tech category") +
  
  scale_y_continuous(limits = c(0.6, 1.0),
                     breaks = seq(0.6, 1, 0.1)) +
  
  labs(title = "Jaccard concordance on MANE isoforms",
       x = NULL, y = "Jaccard index") +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        legend.key.height = unit(0.4, "cm"),
        strip.background  = element_blank(),
        strip.text        = element_text(face = "bold", size = 11),
        axis.text.x       = element_blank(),
        axis.ticks.x      = element_line(size = 0.55),
        plot.title        = element_text(face = "bold", hjust = .5, size = 15)) +
  
  # 手动显著性括号 (与 facet 变量同名列：sample_cat)
  stat_pvalue_manual(
    p_tbl,
    label        = "p.adj.signif",
    xmin         = "xmin",
    xmax         = "xmax",
    y.position   = "y.position",
    bracket.size = 0.35,
    tip.length   = 0.015,
    colour       = "black",   # ← 改成黑色
    size         = 3.2,       # 可顺便调字大小
    inherit.aes  = FALSE
  )
p_jacc

p_jacc+canvas(width = 12,height = 5)

# 保存为 PNG 格式
ggsave("./results/plots/jacc-absolute-isoform-expression-ls-s144.png", plot = p_jacc, width = 12, height = 5, dpi = 300)
ggsave("./results/plots/jacc-absolute-isoform-expression-ls-s144.pdf", plot = p_jacc, width = 12, height = 5, dpi = 300)

## -------- optional: summary table ------------------------------------
jacc_df %>%
  group_by(category,sample_cat) %>%
  summarise(mean_Jaccard = mean(Jaccard),
            median_Jaccard = median(Jaccard),
            n_pairs = n())

# cor
# retain the features whose count >=3 above 90% samples
# 使用 cpm() 函数计算 CPM 矩阵
## ---- 在调用 cpm() 之前加这一段 -----------------------------------
lib_size <- colSums(count_merge, na.rm = TRUE)
# 把 library size = 0 的列剔除，并同步更新 metadata
zero_lib <- names(lib_size)[lib_size == 0]

cpm_df <- cpm(count_merge) %>% as.data.frame()

# 计算样本之间的相关性（使用皮尔逊相关系数）
# 过滤低表达 isoform（以 CPM 为例）
keep <- rowSums(cpm_df >= 1) >= ceiling(0.5 * ncol(cpm_df))# 50% 样本 CPM≥1

cpm_df_filt   <- cpm_df[keep, ]
log2cpm_df    <- log2(cpm_df_filt + 0.01)
cor_matrix <- cor(log2cpm_df, method = "pearson")
# cor_matrix <- cor(log2cpm_df[,-grep("D6",colnames(log2cpm_df))], method = "pearson")

head(cor_matrix)

## -------- 1. reshape correlation matrix to long format ---------------
corr_df <- melt(cor_matrix, varnames = c("sample1","sample2"),
                value.name = "PCC") %>%
  filter(as.character(sample1) < as.character(sample2))

## ───────── 2. 补充技术 & 样本标签 ─────────
tech_vec   <- setNames(metadata$tech,    metadata$sample_id)
sample_vec <- setNames(metadata$sample,  metadata$sample_id)

corr_df <- corr_df %>%
  mutate(
    tech1   = tech_vec[sample1],
    tech2   = tech_vec[sample2],
    tech_cat = case_when(
      tech1 == "long"  & tech2 == "long"  ~ "long-vs-long",
      tech1 == "short" & tech2 == "short" ~ "short-vs-short",
      TRUE                                ~ "long-vs-short"
    ),
    sample_tag1 = sample_vec[sample1],
    sample_tag2 = sample_vec[sample2],
    sample_cat = if_else(sample_vec[sample1] == sample_vec[sample2],
                         paste0("intra-",sample_vec[sample1]),
                         "inter-sample")
  )

corr_df$tech_cat <- factor(corr_df$tech_cat,
                           levels = c("long-vs-long",
                                      "short-vs-short",
                                      "long-vs-short"))


## ───────── 3. 配色 & 透明度 ─────────
fill_cols <- c("intra-D5"      = "#4CC3D9",
               "intra-D6"      = "#7BC8A4",
               "intra-F7"      = "#FFC65D",
               "intra-M8"      = "#F16745",
               "inter-sample"  = "#A58CD9")

alpha_cols <- c("long-vs-long"   = 1.00,
                "short-vs-short" = 0.50,
                "long-vs-short"  = 0.20)

## ───────── 4. n 注释表 ─────────
n_df <- corr_df %>%
  group_by(sample_cat, tech_cat) %>%
  summarise(n   = n(),
            y_n = max(PCC) + 0.02,
            .groups = "drop")

## ───────── 5. 显著性表（long-long vs short-short）─────────
p_tbl <- corr_df %>%                         # ← 确保是 data.frame/tibble
  filter(tech_cat %in% c("long-vs-long", "short-vs-short")) %>%
  compare_means(
    formula  = PCC ~ tech_cat,
    data     = .,                # 关键：把当前管道数据作为 data
    group.by = "sample_cat",
    method   = "wilcox.test"
  ) %>%
  mutate(
    y.position = max(corr_df$PCC) + 0.04,
    xmin       = "short-vs-short",
    xmax       = "long-vs-long"
  )

## ───────── 6. 画图 ─────────
p_cor <- ggplot(corr_df,
            aes(x = tech_cat, y = PCC,
                fill  = sample_cat,
                alpha = tech_cat)) +
  geom_boxplot(width = .55, colour = "grey30", outlier.shape = NA) +
  
  geom_text(data = n_df,
            aes(x = tech_cat, y = y_n, label = paste0("n=", n)),
            vjust = 0, size = 3.2, fontface = "italic", inherit.aes = FALSE) +
  
  facet_wrap(~ sample_cat, nrow = 1, scales = "free_x") +
  
  scale_fill_manual(values = fill_cols, name = "Sample group") +
  scale_alpha_manual(values = alpha_cols, name = "Tech category") +
  
  scale_y_continuous(limits = c(0.4, 1.10),
                     breaks = seq(0.4, 1.0, 0.1)) +
  
  stat_pvalue_manual(
    p_tbl,
    label        = "p.signif",
    xmin         = "xmin", xmax = "xmax", y.position = "y.position",
    tip.length   = .015, bracket.size = .35,
    colour       = "black",       # ← 黑色显著性标记
    inherit.aes  = FALSE
  ) +
  
  labs(title = "Pearson correlation on MANE isoforms",
       x = NULL, y = "PCC (log₂CPM)") +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_line(size = .5),
        plot.title       = element_text(face = "bold", hjust = .5))

print(p_cor)
p_cor+canvas(width = 12,height = 5)

# 保存为 PNG 格式
ggsave("./results/plots/pcc-absolute-isoform-expression-ls-s144.png", plot = p_jacc, width = 12, height = 5, dpi = 300)
ggsave("./results/plots/pcc-absolute-isoform-expression-ls-s144.pdf", plot = p_jacc, width = 12, height = 5, dpi = 300)

## -------- optional: summary table ------------------------------------
corr_df %>%
  group_by(tech_cat,sample_cat) %>%
  summarise(mean_PCC = mean(PCC),
            median_PCC = median(PCC),
            n_pairs = n())

##
pm1 <- p_jacc/p_cor
pm1

#### pcc-heatmap

# 选择颜色映射
colormap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 从 metadata 中提取需要的注释列
annotations <- metadata[, c("tech", "platform", "lib", "sample")]
rownames(annotations) <- metadata$sample_id

# 创建注释颜色列表
annotation_colors <- list(
  tech = c("short" = "#f38d26", "long" = "#4f79a8"),        # tech 注释颜色
  platform = c("BGI" = "#16559a", "ILM" = "#f89e35", "ELE" = "red", "ONT" = "#6c9eb4", "PAB" = "#e01a97"), # platform 注释颜色
  lib = c("P" = "#FFDD44", "R" = "#F781BF", "M" = "#A9DFBF", "D" = "#85C1AE"),            # lib 注释颜色
  sample = c("D5" = "#4CC3D9", "D6" = "#7BC8A4", "F7" = "#FFC65D", "M8" = "#F16745")       # sample 注释颜色
)

# 绘制热图并添加注释
p1 <- pheatmap(mat = cor_matrix,                    # 数据矩阵
               color = colormap,  # 颜色映射
               breaks = seq(0, 1, length.out = length(colormap) + 1),  # Set value range from 0 to 1
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
               main = "Correlation of absolute isoform expression") # 标题
p1
p1+canvas(width = 6,height = 5)

# 保存为 PDF 文件
ggsave("./results/plots/cor-absolute-isoform-expression-ls-s144.png", plot = p1, width = 6, height = 5)
ggsave("./results/plots/cor-absolute-isoform-expression-ls-s144.pdf", plot = p1, width = 6, height = 5)

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

ratio_tx_D6 <- calculate_ratios(log2cpm_df, metadata, "batch_id")

# 计算样本之间的相关性（使用皮尔逊相关系数）
cor_matrix_r <- cor(ratio_tx_D6[,-grep("D6",colnames(ratio_tx_D6))], method = "pearson")
# cor_matrix_r <- cor(ratio_tx_D6, method = "pearson")
head(cor_matrix_r)

# 选择颜色映射
colormap <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 从 metadata 中提取需要的注释列
annotations <- metadata[, c("tech", "platform", "lib", "sample")]
rownames(annotations) <- metadata$sample_id

# 创建注释颜色列表
annotation_colors <- list(
  tech = c("short" = "#f38d26", "long" = "#4f79a8"),        # tech 注释颜色
  platform = c("BGI" = "#16559a", "ILM" = "#f89e35", "ELE" = "red", "ONT" = "#6c9eb4", "PAB" = "#e01a97"), # platform 注释颜色
  lib = c("P" = "#FFDD44", "R" = "#F781BF", "M" = "#A9DFBF", "D" = "#85C1AE"),            # lib 注释颜色
  sample = c("D5" = "#4CC3D9", "D6" = "#7BC8A4", "F7" = "#FFC65D", "M8" = "#F16745")       # sample 注释颜色
)

# 绘制热图并添加注释
p2 <- pheatmap(mat = cor_matrix_r,                    # 数据矩阵
               color = colormap,  # 颜色映射
               breaks = seq(0, 1, length.out = length(colormap) + 1),  # Set value range from 0 to 1
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
               main = "PCC of relative expression LS s144") # 标题
p2
# 保存为 PDF 文件
ggsave("./results/plots/cor-relative-isoform-expression-ls-s108.png", plot = p2, width = 6, height = 5)
ggsave("./results/plots/cor-relative-isoform-expression-ls-s108.pdf", plot = p2, width = 6, height = 5)

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

ggsave("./results/plots/abs-ratio-cor-isoform-expression-lrsr-s144.pdf", plot = combined_plot, width = 12, height = 8)

#### count-heatmap
##──────────────────────── 1. 选取要展示的转录本 ────────────────────────
##   • 只在 ≥50% 样本里 CPM≥1   • 取方差 Top-2000 以便图更清晰
cpm_df  <- cpm(count_merge) |> as.data.frame()
keep    <- rowSums(cpm_df >= 1) >= 0.5 * ncol(cpm_df)
expr_mat <- log2(cpm_df[keep, ] + 0.01)

## 方差排序并截前 2 000（如想全量可跳过此步）
# var_top <- head(order(matrixStats::rowVars(as.matrix(expr_mat)), decreasing = TRUE), 2000)
# expr_top <- expr_mat[var_top, ]         # log₂-CPM 矩阵
expr_top <- expr_mat

##──────────── 2. ratio-to-D6（仍用 log₂ 比值）──────────────
ratio_top <- calculate_ratios(expr_top, metadata, "batch_id")   # 复用你之前的函数

##──────────── 3. 颜色 / 注释一致化──────────────
seq_col <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)   # 连续型

ann_col <- metadata[, c("tech","platform","lib","sample")]
rownames(ann_col) <- metadata$sample_id
ann_cols <- list(
  tech     = c(short = "#f38d26", long = "#4f79a8"),
  platform = c(BGI="#16559a", ILM="#f89e35", ELE="red",
               ONT="#6c9eb4", PAB="#e01a97"),
  lib      = c(P="#FFDD44", R="#F781BF", M="#A9DFBF", D="#85C1AE"),
  sample   = c(D5="#4CC3D9", D6="#7BC8A4",
               F7="#FFC65D", M8="#F16745")
)

##──────────── 4. heatmap：绝对表达 ──────────────
expr_top <- expr_top[apply(expr_top, 1, sd, na.rm = TRUE) != 0, ]
ht_abs <- pheatmap(expr_top,
                   color = seq_col,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   show_rownames = FALSE, show_colnames = FALSE,
                   border_color = NA, legend = TRUE,
                   scale = "row",                       # 每行 0-mean / sd
                   annotation_col = ann_col,
                   annotation_colors = ann_cols,
                   fontsize = 11,
                   main = "log2CPM (absolute)")
ht_abs

##──────────── 5. heatmap：ratio 表达──────────────
div_col <- colorRampPalette(rev(brewer.pal(11, "PuOr")))(100)   # 发散型调色

ht_ratio <- pheatmap(ratio_top,
                     color = div_col,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = FALSE, show_colnames = FALSE,
                     border_color = NA, legend = TRUE,
                     scale = "row",
                     annotation_col = ann_col,
                     annotation_colors = ann_cols,
                     fontsize = 11,
                     main = "log2CPM (ratio to D6)")
ht_ratio

##──────────── 6. 并排输出（cowplot / patchwork）──────────────
library(ggplotify); library(cowplot)

abs_g  <- as.ggplot(ht_abs$gtable)
ratio_g<- as.ggplot(ht_ratio$gtable)

pm2 <- plot_grid(abs_g + theme(legend.position="none"),
          ratio_g + theme(legend.position="none"),
          ncol = 2, align = "h", rel_widths = c(1,1),
          labels = c("",""))
pm2
ggsave("./results/plots/Heatmap_abs_ratio_MANE_s144.pdf",pm2, width = 12, height = 8)
ggsave("./results/plots/Heatmap_abs_ratio_MANE_s144.png",pm2, width = 12, height = 8)







## boxlplot-PCC

# Prepare data for the absolute PCC (cor_matrix)
cor_df_absolute <- as.data.frame(as.table(cor_matrix))
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
    
    # 确定技术类型
    tech1 = case_when(
      grepl("BGI|ILM", Sample1) ~ "short",
      grepl("ONT|PAB", Sample1) ~ "long",
      TRUE ~ "unknown"
    ),
    tech2 = case_when(
      grepl("BGI|ILM", Sample2) ~ "short",
      grepl("ONT|PAB", Sample2) ~ "long",
      TRUE ~ "unknown"
    ),
    
    # Class逻辑
    Class = case_when(
      # Cross-tech
      tech1 != tech2 ~ "Cross-tech",
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
      # Intra-batch 分为 intra-sample 和 cross-sample
      Class == "Intra-batch" & sample_type_1 == sample_type_2 ~ "Intra-batch intra-sample",
      Class == "Intra-batch" & sample_type_1 != sample_type_2 ~ "Intra-batch cross-sample",
      
      # Cross-tech 分为 intra-sample 和 cross-sample
      Class == "Cross-tech" & sample_type_1 == sample_type_2 ~ "Cross-tech intra-sample",
      Class == "Cross-tech" & sample_type_1 != sample_type_2 ~ "Cross-tech cross-sample",
      
      # Cross-platform 分为 intra-sample 和 cross-sample
      Class == "Cross-platform" & sample_type_1 == sample_type_2 ~ "Cross-platform intra-sample",
      Class == "Cross-platform" & sample_type_1 != sample_type_2 ~ "Cross-platform cross-sample",
      
      # Cross-protocol 分为 intra-sample 和 cross-sample
      Class == "Cross-protocol" & sample_type_1 == sample_type_2 ~ "Cross-protocol intra-sample",
      Class == "Cross-protocol" & sample_type_1 != sample_type_2 ~ "Cross-protocol cross-sample",
      
      # Cross-lab 分为 intra-sample 和 cross-sample
      Class == "Cross-lab" & sample_type_1 == sample_type_2 ~ "Cross-lab intra-sample",
      Class == "Cross-lab" & sample_type_1 != sample_type_2 ~ "Cross-lab cross-sample",
      
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
base_colors <- c("#FB9A99", "#e53531", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")[seq_along(unique_classes)]

# 5. 生成颜色映射：同一 Class 的 cross-sample 半透明、intra-sample 不透明
color_mapping <- sapply(levels(cor_df_absolute$ClassGroup), function(cg) {
  parts <- strsplit(cg, " - ")[[1]]
  cl <- parts[1]
  gt <- parts[2]
  col <- base_colors[which(unique_classes == cl)]
  if(gt == "cross-sample") alpha(col, 0.5) else col
})
names(color_mapping) <- levels(cor_df_absolute$ClassGroup)

# 6. 绘图：左右小提琴 + 中间箱线图
p3 <- ggplot(cor_df_absolute, aes(x = Class, y = PCC)) +
  geom_split_violin(aes(fill = ClassGroup), trim = TRUE, position = position_dodge(width = 0.2)) +
  geom_boxplot(aes(fill = ClassGroup, group = Class), width = 0.1, outlier.shape = NA, color = "black") +
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
ggsave("./results/plots/pcc-boxplot-absolute-isoform-expression-ls-s120.png", plot = p3, width = 7, height = 6)

# Prepare data for the relative PCC (cor_matrix_r)
cor_df_relative <- as.data.frame(as.table(cor_matrix_r))
colnames(cor_df_relative) <- c("Sample1", "Sample2", "PCC")
cor_df_relative <- cor_df_relative %>%
  filter(Sample1 != Sample2)
test <- cor_df_relative[cor_df_relative$PCC<0,]
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
    
    # 确定技术类型
    tech1 = case_when(
      grepl("BGI|ILM", Sample1) ~ "short",
      grepl("ONT|PAB", Sample1) ~ "long",
      TRUE ~ "unknown"
    ),
    tech2 = case_when(
      grepl("BGI|ILM", Sample2) ~ "short",
      grepl("ONT|PAB", Sample2) ~ "long",
      TRUE ~ "unknown"
    ),
    
    # Class逻辑
    Class = case_when(
      # Cross-tech
      tech1 != tech2 ~ "Cross-tech",
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
      # Intra-batch 分为 intra-sample 和 cross-sample
      Class == "Intra-batch" & sample_type_1 == sample_type_2 ~ "Intra-batch intra-sample",
      Class == "Intra-batch" & sample_type_1 != sample_type_2 ~ "Intra-batch cross-sample",
      
      # Cross-tech 分为 intra-sample 和 cross-sample
      Class == "Cross-tech" & sample_type_1 == sample_type_2 ~ "Cross-tech intra-sample",
      Class == "Cross-tech" & sample_type_1 != sample_type_2 ~ "Cross-tech cross-sample",
      
      # Cross-platform 分为 intra-sample 和 cross-sample
      Class == "Cross-platform" & sample_type_1 == sample_type_2 ~ "Cross-platform intra-sample",
      Class == "Cross-platform" & sample_type_1 != sample_type_2 ~ "Cross-platform cross-sample",
      
      # Cross-protocol 分为 intra-sample 和 cross-sample
      Class == "Cross-protocol" & sample_type_1 == sample_type_2 ~ "Cross-protocol intra-sample",
      Class == "Cross-protocol" & sample_type_1 != sample_type_2 ~ "Cross-protocol cross-sample",
      
      # Cross-lab 分为 intra-sample 和 cross-sample
      Class == "Cross-lab" & sample_type_1 == sample_type_2 ~ "Cross-lab intra-sample",
      Class == "Cross-lab" & sample_type_1 != sample_type_2 ~ "Cross-lab cross-sample",
      
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
# 将 cor_df_relative 的 Class 按照 p3 的顺序设定
cor_df_relative$Class <- factor(cor_df_relative$Class, levels = levels(cor_df_absolute$Class))

# 定义 GroupType 顺序（先 cross-sample，再 intra-sample）
group_order <- c("cross-sample", "intra-sample")

# 创建合并因子 ClassGroup，因子水平顺序与 p3 保持一致
cor_df_relative$ClassGroup <- factor(
  paste(cor_df_relative$Class, cor_df_relative$GroupType, sep = " - "),
  levels = unlist(lapply(levels(cor_df_absolute$Class), function(cl) {
    paste(cl, group_order, sep = " - ")
  }))
)

# 绘制图形：左右小提琴 + 中间箱线图，图例显示合并后的 ClassGroup
p4 <- ggplot(cor_df_relative, aes(x = Class, y = PCC)) +
  geom_split_violin(aes(fill = ClassGroup), trim = TRUE, position = position_dodge(width = 0.2)) +
  geom_boxplot(aes(fill = ClassGroup, group = Class), 
               width = 0.1, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = color_mapping) +
  labs(
    title = "After correction",
    x = "Class",
    y = "Pearson correlation coefficient",
    fill = "Class - Group Type"
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

print(p4)
# 保存为 PDF 文件
ggsave("./results/plots/pcc-boxplot-relative-isoform-expression-ls-s120.png", plot = p4, width = 7, height = 6)

# merge plot
p3_n <- p3+theme(legend.position = "none")
p4_n <- p4+theme(legend.position = "none")


pm <- p3_n+p4_n
pm
# 保存为 PDF 文件
ggsave("./results/plots/pcc-boxplot-abs-relative-isoform-expression-lr-s120.png", plot = pm, width = 10, height = 5)

# remove all
rm(list=ls())
gc()
