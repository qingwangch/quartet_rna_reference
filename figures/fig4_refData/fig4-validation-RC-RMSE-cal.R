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

## ── 1. 通用函数 ──────────────────────────────────────────────────────────────
calc_metrics_isoform <- function(ref_df, test_df){
  inner_join(
    ref_df  %>% select(isoform_compare, compare, source, log2FC_ref = log2FC),
    test_df %>% select(isoform_compare, compare, log2FC_val = log2FC),
    by = c("isoform_compare","compare")
  ) %>%
    group_by(compare, source) %>%
    summarise(
      RC   = cor(log2FC_ref, log2FC_val, use = "complete.obs"),
      RMSE = sqrt(mean((log2FC_val - log2FC_ref)^2, na.rm = TRUE)),
      n    = n(), .groups = "drop"
    )
}

calc_metrics_ase <- function(ref_df, test_df) {
  clamp <- function(x) pmin(pmax(x, 0.001), 0.999)   # 防止 ±Inf
  inner_join(
    ref_df  |> select(ASE_compare, compare, source,
                      psi_ref = mean_delta_psi_mean),
    test_df |> select(ASE_compare, 
                      psi_val = delta_PSI),
    by = "ASE_compare"
  ) |>
    mutate(
      psi_ref_logit = log(clamp(psi_ref) / (1 - clamp(psi_ref))),
      psi_val_logit = log(clamp(psi_val) / (1 - clamp(psi_val)))
    ) |>
    group_by(compare, source) |>
    summarise(
      RC        = cor(psi_ref_logit, psi_val_logit,
                      use = "complete.obs"),
      RMSE = sqrt(mean((psi_val_logit - psi_ref_logit)^2,
                            na.rm = TRUE)),
      n         = n(),
      .groups   = "drop"
    )
}


## ── 2. 读取参考集 ────────────────────────────────────────────────────────────
# Isoform (refFC / refDEI)
ref_fc  <- fread("isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv") |> mutate(source = "LO")
ref_dei <- fread("isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv") |>
  filter(Final != "non-DEI") |> mutate(source = "LO")

# ASE (refPSI / refDAS)
ref_psi <- fread("suppa2/ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv") |> 
  mutate(
  ASE_compare = paste(ASE, compare),   # ★补这一列
  source      = "LO"
)
ref_das <- fread("suppa2/ref_data_construction/lo-e-g/RefData_AS_all_DAS_classified_u_0520.csv") |>
  filter(feature_type != "non-DAS") |>
  mutate(ASE_compare = paste(ASE, compare)) |> mutate(source = "LO")

## ── 3. 读取 qPCR 验证集 ─────────────────────────────────────────────────────
qpcr_iso <- fread("isoforms/qpcr/log2fc_20250409_U.csv")   |> 
  mutate(source = "qPCR")                        # ← 添加 source
qpcr_ase <- fread("isoforms/qpcr/deltapsi_20250513_U.csv") |> 
  mutate(source = "qPCR")

## ── 4. 计算 RC / RMSE ───────────────────────────────────────────────────────
iso_metrics <- bind_rows(
  calc_metrics_isoform(ref_fc , qpcr_iso) |> mutate(type = "RefFC_vs_qPCR"),
  calc_metrics_isoform(ref_dei, qpcr_iso) |> mutate(type = "RefDEI_vs_qPCR")
)

ase_metrics <- bind_rows(
  calc_metrics_ase(ref_psi, qpcr_ase) |> mutate(type = "RefPSI_vs_qPCR"),
  calc_metrics_ase(ref_das, qpcr_ase) |> mutate(type = "RefDAS_vs_qPCR")
)

all_metrics <- bind_rows(
  iso_metrics |> mutate(metric = "Isoform"),
  ase_metrics |> mutate(metric = "ASE")
)

# 写出结果
fwrite(all_metrics, "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/data/fig4b_validated_metrics_RefData.csv")

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

p+canvas(width = 8,height = 4)

# 假设你的图对象叫 p
ggsave(
  filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/iso_as_validation_qpcr.png",
  plot     = p,
  width    = 15,    # 宽度（英寸）
  height   = 8,    # 高度（英寸）
  dpi      = 300   # 分辨率
)

ggsave(
  filename = "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/iso_as_validation_qpcr.pdf",
  plot     = p,
  width    = 15,    # 宽度（英寸）
  height   = 8,    # 高度（英寸）
  dpi      = 300   # 分辨率
)

# remove all
rm(list=ls())
gc()
