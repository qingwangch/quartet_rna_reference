# validation-RC-RMSE-cal.R
# 2025-06-13
# Qingwang Chen

# 工作目录 & R 库路径
setwd("/vast/projects/quartet_rna_refdata/analysis/")
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"
.libPaths(c(custom_lib_path, .libPaths()))

# 加载必要的包
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(purrr)

# ——— 通用函数 ———

# 针对 isoform 的 RC+RMSE 计算
calc_metrics_isoform <- function(ref_df, test_df) {
  inner_join(
    ref_df %>% select(isoform_compare, compare, source, log2FC_ref = log2FC),
    test_df %>% rename(log2FC_val = log2FC) %>% select(isoform_compare, compare, log2FC_val),
    by = c("isoform_compare", "compare")
  ) %>%
    group_by(compare, source) %>%
    summarise(
      RC       = cor(log2FC_ref, log2FC_val, use = "complete.obs"),
      RMSE     = sqrt(mean((log2FC_val - log2FC_ref)^2, na.rm = TRUE)),
      n_points = n(),
      .groups  = "drop"
    )
}

# 针对可变剪接（ASE）的 RC+RMSE 计算
calc_metrics_ase <- function(ref_df, test_df) {
  inner_join(
    ref_df %>% rename(psi_ref = mean_delta_psi_mean) %>% select(ASE_compare, compare, source, psi_ref),
    test_df %>% rename(psi_val = delta_PSI) %>% select(ASE_compare, psi_val),
    by = c("ASE_compare")
  ) %>%
    group_by(compare, source) %>%
    summarise(
      RC       = cor(psi_ref, psi_val, use = "complete.obs"),
      RMSE     = sqrt(mean((psi_val - psi_ref)^2, na.rm = TRUE)),
      n_points = n(),
      .groups  = "drop"
    )
}

# ——— 数据读取 ———

# Isoform 参考集
ref_iso_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv")
ref_dei_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv") %>%
  filter(Final != "non-DEI")

# ASE 参考集
ref_ase_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv")
ref_das_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/RefData_AS_all_DAS_classified_u_0520.csv") %>%
  filter(feature_type != "non-DAS")%>%
  mutate(ASE_compare = paste(ASE, compare))

# qPCR 验证数据（Isoform + ASE 共用）
qpcr_iso <- fread("isoforms/qpcr/log2fc_20250409_U.csv")   %>% as.data.frame()
qpcr_ase <- fread("isoforms/qpcr/deltapsi_20250513_U.csv") %>% as.data.frame()

# dRNA 验证数据
drna_iso_raw <- fread("isoforms/de/isoform/D_ONT_LW_B1/DE_results_Oarfish_isoform.txt") %>% as.data.frame()
drna_iso <- drna_iso_raw %>%
  rename(log2FC = logFC) %>%
  mutate(
    isoform         = sub("\\|.*$", "", feature),
    compare         = gsub("vs", "/", comparison),
    isoform_compare = paste(isoform, compare),
    PValue          = PValue
  ) %>%
  select(isoform_compare, compare, log2FC)

drna_ase_raw <- fread("suppa2/oarfish/suppa2_diffsplice/D_ONT_LW_B1/leg_DAS_drna.txt") %>% as.data.frame()
drna_ase <- drna_ase_raw %>%
  rename(delta_PSI = delta_psi, compare = compare) %>%
  filter(
    !is.na(delta_PSI),
    abs(delta_PSI) >= 0.05
  ) %>%
  mutate(ASE_compare = paste(ase, compare) ) %>% 
  select(ASE_compare, compare, delta_PSI)

# ——— 批量计算并合并 ———
# 1) 定义参考集和验证集
isoform_refs <- list(
  RefFC  = ref_iso_ls,
  RefDEI = ref_dei_ls
)
validations_iso <- list(
  qPCR  = qpcr_iso,
  dRNA  = drna_iso
)

ase_refs <- list(
  RefPSI = ref_ase_ls,
  RefDAS = ref_das_ls
)
validations_ase <- list(
  qPCR  = qpcr_ase,
  dRNA  = drna_ase
)

# 2) 用基础 for 循环跑 Isoform 部分
iso_metrics <- list()
cnt <- 1
for(rtype in names(isoform_refs)) {
  ref_df <- isoform_refs[[rtype]]
  for(vname in names(validations_iso)) {
    test_df <- validations_iso[[vname]]
    # 计算 RC/RMSE
    metr <- calc_metrics_isoform(ref_df, test_df)
    # 打上类型标签
    metr$type <- paste0(rtype, "_vs_", vname)
    iso_metrics[[cnt]] <- metr
    cnt <- cnt + 1
  }
}
all_metrics_iso <- bind_rows(iso_metrics)

# 3) 用同样的方法跑 ASE 部分
ase_metrics <- list()
cnt <- 1
for(rtype in names(ase_refs)) {
  ref_df <- ase_refs[[rtype]]
  for(vname in names(validations_ase)) {
    test_df <- validations_ase[[vname]]
    metr <- calc_metrics_ase(ref_df, test_df)
    metr$type <- paste0(rtype, "_vs_", vname)
    ase_metrics[[cnt]] <- metr
    cnt <- cnt + 1
  }
}
all_metrics_ase <- bind_rows(ase_metrics)

# 4) 合并两部分并写出
all_metrics <- bind_rows(
  all_metrics_iso %>% mutate(metric = "isoform"),
  all_metrics_ase %>% mutate(metric = "ASE")
)

# 写出结果
fwrite(all_metrics, "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/data/fig4e_validated_metrics_RefData.csv")

# 可视化验证结果
# 1) 拉长，把 RC/ RMSE 变成一列
metrics_long <- all_metrics %>%
  pivot_longer(
    cols      = c(RC, RMSE),
    names_to  = "measure",
    values_to = "value"
  ) %>%
  mutate(
    compare = factor(compare, levels = c("D5/D6","F7/D6","M8/D6")),
    type    = factor(type,   levels = c(
      "RefFC_vs_qPCR","RefFC_vs_dRNA",
      "RefDEI_vs_qPCR","RefDEI_vs_dRNA",
      "RefPSI_vs_qPCR","RefPSI_vs_dRNA",
      "RefDAS_vs_qPCR","RefDAS_vs_dRNA"
    )),
    measure = factor(measure, levels = c("RC","RMSE"))
  )

# 2) 分面标签
facet_type_labels <- c(
  RefFC_vs_qPCR   = "Isoform: qPCR",
  RefFC_vs_dRNA   = "Isoform: dRNA",
  RefDEI_vs_qPCR  = "DEI: qPCR",
  RefDEI_vs_dRNA  = "DEI: dRNA",
  RefPSI_vs_qPCR  = "ASE: qPCR",
  RefPSI_vs_dRNA  = "ASE: dRNA",
  RefDAS_vs_qPCR  = "DAS: qPCR",
  RefDAS_vs_dRNA  = "DAS: dRNA"
)
facet_measure_labels <- c(RC = "RC", RMSE = "RMSE")

# 3) 配色
nbt_palette <- c("D5/D6"="#1b9e77","F7/D6"="#d95f02","M8/D6"="#7570b3")

# 4) 绘图：行是指标，列是 type，柱子按 sample pair 填色
p <- ggplot(metrics_long, aes(x = compare, y = value, fill = source)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  facet_grid(
    rows    = vars(measure),
    cols    = vars(type),
    scales  = "free_y",
    labeller = labeller(
      type    = facet_type_labels,
      measure = facet_measure_labels
    )
  ) +
  scale_fill_manual(
    name   = "Strategy",
    values = c(LO = "#1f77b4", LS = "#ff7f0e", SO = "#2ca02c")  # 这里替换成你的 nbt_palette 中对应 source 的配色
  ) +
  labs(
    x    = "Sample pair",
    y    = "Value"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major    = element_line(color = "grey80", size = 0.5),
    panel.grid.minor    = element_line(color = "grey90", size = 0.25),
    axis.line           = element_line(color = "#222222", size = 0.6),
    axis.ticks          = element_line(color = "#222222", size = 0.6),
    axis.ticks.length   = unit(0.2, "cm"),
    strip.text          = element_text(face = "bold", size = 12, colour = "#FFFFFF"),
    strip.background    = element_rect(fill = "#333333", colour = NA),
    axis.title.x        = element_text(face = "bold", size = 13, margin = margin(t = 8)),
    axis.title.y        = element_text(face = "bold", size = 13, margin = margin(r = 8)),
    axis.text.x         = element_text(size = 11, face = "bold"),
    axis.text.y         = element_text(size = 11),
    legend.position     = "bottom",
    legend.title        = element_text(face = "bold", size = 12),
    legend.text         = element_text(size = 11),
    panel.background    = element_rect(fill = "white", colour = NA),
    plot.margin         = margin(5, 5, 5, 5)
  )

print(p)

# 假设你的图对象叫 p
ggsave(
  filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/iso_as_validation_qpcr_drna.png",
  plot     = p,
  width    = 15,    # 宽度（英寸）
  height   = 8,    # 高度（英寸）
  dpi      = 300   # 分辨率
)

ggsave(
  filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/iso_as_validation_qpcr_drna.pdf",
  plot     = p,
  width    = 15,    # 宽度（英寸）
  height   = 8,    # 高度（英寸）
  dpi      = 300   # 分辨率
)

# remove all
rm(list=ls())
gc()
