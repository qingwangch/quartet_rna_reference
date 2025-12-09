#' fig4-statistics-of-reference-isoforms-ASEs
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

# install.packages("Rcpp", type = "source", lib=custom_lib_path)
# library(Rcpp)
# library import
library(edgeR)
library(arrow)
library(ggplot2)
library(reshape2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(grid) 

# data import
# 通用读取函数
read_and_count <- function(file, type, dataset, skip_filter = FALSE) {
  df <- read.csv(file)
  colnames(df) <- tolower(colnames(df))  # 标准化列名
  
  if (!"compare" %in% colnames(df)) stop("Missing 'compare' column in: ", file)
  
  if (!skip_filter && type %in% c("refDEI", "refDAS")) {
    if (!"final" %in% colnames(df)) stop("Missing 'final' column in: ", file)
    df <- df %>% filter(final != ifelse(type == "refDEI", "non-DEI", "non-DAS"))
  }
  
  df %>%
    count(compare, name = "count") %>%
    mutate(type = type, dataset = dataset)
}

# 加载所有数据
get_combined_data <- function(mode = c("isoform", "ase")) {
  mode <- match.arg(mode)
  
  if (mode == "isoform") {
    files <- list(
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b6_p4_s72_u_20250403.csv", "refFC", "LO", TRUE),
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250407.csv", "refDEI", "LO", FALSE),
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/ref_expr_b13_p4_s156_u_20250407.csv", "refFC", "SO", TRUE),
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/RefData_DEIs_all_isoforms_classified_u_20250407.csv", "refDEI", "SO", FALSE),
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/ref_expr_LS_t2_u_20250407.csv", "refFC", "LS", TRUE),
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/RefData_DEIs_LS_t2_u_20250407.csv", "refDEI", "LS", FALSE)
    )
  } else {
    files <- list(
      list("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/ref_expr_as_b6_p4_s72_u_20250425.csv", "refPSI", "LO", TRUE),
      list("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/RefData_AS_all_DAS_classified_u_0425.csv", "refDAS", "LO", FALSE),
      list("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/so-e-g/ref_expr_as_b13_p4_s72_u_20250425.csv", "refPSI", "SO", TRUE),
      list("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/so-e-g/RefData_AS_all_DAS_classified_u_0425.csv", "refDAS", "SO", FALSE),
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/ref_expr_ase_LS_t2_u_20250425.csv", "refPSI", "LS", TRUE),
      list("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/RefData_DAS_LS_t2_u_20250425.csv", "refDAS", "LS", FALSE)
    )
  }
  
  data_list <- lapply(files, function(x) {
    read_and_count(file = x[[1]], type = x[[2]], dataset = x[[3]], skip_filter = x[[4]])
  })
  
  bind_rows(data_list)
}

# 绘图函数
# 你原来的绘图函数，稍作修改，在函数开始处保存数据
plot_count_comparison_nbt <- function(data, type_fc, type_dei, ylab_fc, ylab_dei,
                                      out_prefix = "./results/plots/count_comparison") {
  # --- 1. 预处理 ---
  data <- data %>%
    mutate(
      compare = factor(compare, levels = c("D5/D6", "F7/D6", "M8/D6")),
      dataset = factor(dataset, levels = c("LO", "SO", "LS")),
      type    = factor(type,    levels = c(type_fc, type_dei))
    )
  
  # 先把原始整理好的 data 保存下来
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  write.csv(data,
            file = paste0(out_prefix, ".csv"),
            row.names = FALSE)
  
  # --- 2. 缩放 second y-axis ---
  scale_factor <- max(data$count[data$type == type_fc]) /
    max(data$count[data$type == type_dei])
  
  data_scaled <- data %>%
    mutate(scaled_count = ifelse(type == type_dei,
                                 count * scale_factor,
                                 count))
  
  # 再把用于绘图的 data_scaled 保存下来
  write.csv(data_scaled,
            file = paste0(out_prefix, "_scaled.csv"),
            row.names = FALSE)
  
  # --- 3. 画图 ---
  nbt_palette <- c(
    "D5/D6" = "#1b9e77",
    "F7/D6" = "#d95f02",
    "M8/D6" = "#7570b3"
  )
  facet_labels <- c(type_fc = ylab_fc, type_dei = ylab_dei)
  
  p <- ggplot(data_scaled, aes(x = dataset, y = scaled_count, fill = compare)) +
    geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.6) +
    facet_wrap(~ type, nrow = 1, scales = "free_x",
               labeller = labeller(type = facet_labels)) +
    scale_fill_manual(values = nbt_palette) +
    scale_y_continuous(
      name     = ylab_fc,
      sec.axis = sec_axis(~ . / scale_factor, name = ylab_dei)
    ) +
    labs(x = "Dataset", fill = "Sample pair") +
    theme_minimal(base_size = 14) +
    theme(
      text                = element_text(family = "Helvetica", color = "#222222"),
      strip.text          = element_text(face = "bold", size = 14, color = "#FFFFFF"),
      strip.background    = element_rect(fill = "#333333", color = NA),
      axis.title.x        = element_text(face = "bold", size = 14, margin = margin(t = 8)),
      axis.title.y.left   = element_text(face = "bold", size = 14, margin = margin(r = 8)),
      axis.title.y.right  = element_text(face = "bold", size = 14, margin = margin(l = 8)),
      axis.text.x         = element_text(size = 12, face = "bold"),
      axis.text.y         = element_text(size = 12),
      axis.line           = element_line(color = "#222222", size = 0.6),
      axis.ticks          = element_line(color = "#222222", size = 0.6),
      axis.ticks.length   = unit(0.25, "cm"),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.background    = element_rect(fill = "white", color = NA),
      legend.position     = "top",
      legend.direction    = "horizontal",
      legend.key.size     = unit(0.6, "cm"),
      legend.key.width    = unit(1, "cm"),
      legend.background   = element_blank(),
      legend.title        = element_text(face = "bold", size = 12),
      legend.text         = element_text(size = 11),
      plot.margin         = margin(10, 10, 10, 10)
    )
  
  return(p)
}

# === 主流程 ===
# Isoform
isoform_data <- get_combined_data("isoform")
p1 <- plot_count_comparison_nbt(
  data     = isoform_data, 
  type_fc  = "refFC", 
  type_dei = "refDEI", 
  ylab_fc  = "refFC Isoform Count", 
  ylab_dei = "refDEI Count",
  out_prefix= "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/count_ref_Isoform"
)

p1

# ASE
ase_data <- get_combined_data("ase")
p2 <- plot_count_comparison_nbt(
  ase_data,
  type_fc  = "refPSI",
  type_dei = "refDAS",
  ylab_fc  = "refPSI ASE Count",
  ylab_dei = "refDAS Count",
  out_prefix= "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/count_ref_ASE"
)
p2

# 拼图输出
p1 <- p1 +
  theme(
    axis.title.x = element_blank()  # 把上面子图的 x 轴标题去掉
  )
pb <- (p1 + p2) +
  plot_layout(ncol = 1,    # 垂直排列
              guides = "collect") &  # 收集（合并）所有子图的图例
  theme(legend.position = "top")  # 只在最上面保留一次图例
pb
# 保存
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/fig4b.png", pb, width = 8, height = 9, dpi = 300)

# merge data
count_ref_Isoform <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/count_ref_Isoform.csv") %>% as.data.frame()
count_ref_ASE <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/count_ref_ASE.csv") %>% as.data.frame()
count_refData <- rbind(count_ref_Isoform,count_ref_ASE)
fwrite(count_refData,"/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/data/fig4b_count_refdata.csv")

############### Fig4c & 4d
# NBT 配色
nbt_colors <- c(
  "D5/D6" = "#1b9e77",
  "F7/D6" = "#d95f02",
  "M8/D6" = "#7570b3"
)

# 通用：绘制相关性散点图
# df_x, df_y: 两个 data.frame
# by: 合并键
# suffixes: 合并后来自 x,y 的列名后缀
# x_var, y_var: 合并前在 df_x/df_y 中的测度列名
# group_var: 用于着色的分组列名（同样在 df_x/df_y 中）
# lab_prefix: subtitle 中的前缀（如 "Isoform" 或 "ASE"）
# xlab, ylab: 坐标轴标题
plot_corr <- function(df_x, df_y, by = "isoform_compare",
                      suffixes = c("_so","_lo"),
                      x_var, y_var, group_var,
                      lab_prefix, xlab, ylab) {
  merged <- merge(df_x, df_y, by = by, suffixes = suffixes)
  
  # 如果是 ΔPSI，需要先重命名列
  # （假如原列叫 mean_delta_psi_mean）
  if (x_var == "delta_psi") {
    names(merged)[names(merged) == paste0("mean_delta_psi_mean", suffixes[1])] <- "delta_psi_so"
    x_var_full <- "delta_psi_so"
  } else {
    x_var_full <- paste0(x_var, suffixes[1])
  }
  if (y_var == "delta_psi") {
    names(merged)[names(merged) == paste0("mean_delta_psi_mean", suffixes[2])] <- "delta_psi_lo"
    y_var_full <- "delta_psi_lo"
  } else {
    y_var_full <- paste0(y_var, suffixes[2])
  }
  
  # 对 ΔPSI 图做过滤：只保留 |ΔPSI|>0.05
  if (x_var == "delta_psi") {
    merged <- merged %>%
      filter(abs(.data[[x_var_full]]) > 0.05,
             abs(.data[[y_var_full]]) > 0.05)
  }
  
  n <- nrow(merged)
  r <- cor(merged[[x_var_full]], merged[[y_var_full]], use = "complete.obs")
  group_full <- paste0(group_var, suffixes[1])
  
  ggplot(merged, aes_string(x = x_var_full, y = y_var_full, color = group_full)) +
    geom_point(alpha = 0.4, size = 2) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
    scale_color_manual(values = nbt_colors) +
    labs(
      subtitle = paste0(lab_prefix, ", N = ", n),
      x = xlab, y = ylab, color = "Comparison"
    ) +
    annotate("text",
             x = min(merged[[x_var_full]], na.rm = TRUE),
             y = max(merged[[y_var_full]], na.rm = TRUE),
             label = paste0("r = ", round(r, 2)),
             hjust = 0, vjust = 1, fontface = "bold", size = 4) +
    theme_minimal(base_size = 14) +
    theme(
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 12),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm")
    )
}

# ——————————————
# 1) Isoform log2FC（Fig4c）
# ——————————————
fc_lo_csv   <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b6_p4_s72_u_20250403.csv"
fc_so_csv   <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/ref_expr_b13_p4_s156_u_20250407.csv"

iso_lo <- read.csv(fc_lo_csv, stringsAsFactors = FALSE)
iso_so <- read.csv(fc_so_csv, stringsAsFactors = FALSE)

p4c <- plot_corr(
  df_x = iso_so, df_y = iso_lo,
  x_var = "log2FC", y_var = "log2FC", group_var = "compare",
  lab_prefix = "Isoform", 
  xlab = "log2FC (SO reference)", ylab = "log2FC (LO reference)"
)
p4c
# ——————————————
# 2) ASE ΔPSI（Fig4d）
# ——————————————
das_lo_csv <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/RefData_AS_all_DAS_classified_u_0425.csv"
das_so_csv <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/so-e-g/RefData_AS_all_DAS_classified_u_0425.csv"

das_lo <- read.csv(das_lo_csv) %>% filter(abs(mean_delta_psi_mean) > 0.05)
das_so <- read.csv(das_so_csv) %>% filter(abs(mean_delta_psi_mean) > 0.05)

# 注意：这两个表中 ΔPSI 列叫 mean_delta_psi_mean
# 我们在绘图函数里会重命名
p4d <- plot_corr(
  df_x = das_so, df_y = das_lo,
  x_var = "delta_psi", y_var = "delta_psi", group_var = "compare",
  lab_prefix = "ASE", 
  xlab = expression(Delta*"PSI (SO reference)"),
  ylab = expression(Delta*"PSI (LO reference)")
)
p4d

# ——————————————
# 3) 合并并输出
# ——————————————
pcd <- (p4c + p4d) +
  plot_layout(ncol = 1,    # 垂直排列
              guides = "collect") &  # 收集（合并）所有子图的图例
  theme(legend.position = "top")  # 只在最上面保留一次图例
pcd
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/fig4cd.png", pcd, width = 4, height = 9, dpi = 300)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/fig4cd.pdf", pcd, width = 4, height = 9, dpi = 300)

## remove all
rm(list=ls())
gc()
