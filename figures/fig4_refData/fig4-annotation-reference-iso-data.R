# annotation-reference-iso-data.R
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
library(tidytext) 

## ── 2. 读取参考集 ────────────────────────────────────────────────────────────
# Isoform (refFC / refDEI)
ref_fc  <- fread("isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv") |> mutate(source = "LO")
ref_dei <- fread("isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv") |>
  filter(Final != "non-DEI") |> mutate(source = "LO")
count <- dplyr::count
## ------------- 1. 统计 Top-8 transcript type & “Other” ----------------
## 1. 重新统计 —— 只按 transcript_type 排序 -------------------------
top_n <- 8
top_types <- ref_fc %>%
  count(transcript_type, sort = TRUE) %>%
  slice_head(n = top_n) %>%
  pull(transcript_type)

iso_dot_df <- ref_fc %>%
  mutate(
    type_raw = if_else(transcript_type %in% top_types,
                       transcript_type, "Other"),
    compare  = factor(compare, levels = c("D5/D6", "F7/D6", "M8/D6"))
  ) %>%
  count(compare, type_raw, name = "n") %>%
  ## 现在 x 轴只用 type_raw，本身就唯一
  mutate(type_raw = fct_reorder(type_raw, n, .fun = sum, .desc = TRUE))

## 2. 新配色 ----------------------------------------------------------
type_cols <- c(
  "protein_coding"             = "#ffd92f",
  "lncRNA"                     = "#4DAF4A",   # ← 替换
  "retained_intron"            = "#E41A1C",   # ← 替换
  "nonsense_mediated_decay"    = "#984EA3",   # ← 替换
  "processed_transcript"       = "#e7298a",
  "protein_coding_CDS_not_def" = "#66a61e",
  "IG_C_gene"                  = "#e6ab02",
  "TEC"                        = "#a6761d",
  "Other"                      = "#666666"
)

## 3. 画 dot-plot -----------------------------------------------------
p_dot <- ggplot(iso_dot_df,
                aes(x = type_raw, y = compare,
                    size = n, fill = type_raw, colour = type_raw)) +
  geom_point(shape = 21, stroke = .35) +
  geom_text(aes(label = n), colour = "white",
            size = 2.8, fontface = "bold") +
  scale_size(range = c(4, 16), guide = "none") +
  scale_fill_manual(values = type_cols, guide = "none") +
  scale_colour_manual(values = type_cols, guide = "none") +
  labs(title = "Transcript-type composition (refFC isoforms, LO)",
       x = "Transcript type (Top 8 + Other)",
       y = "Sample pair") +
  theme_classic(base_size = 12.5) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 9.5),
    axis.text.y  = element_text(size = 11, face = "bold"),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    plot.title   = element_text(face = "bold", size = 14, hjust = .5)
  )

print(p_dot)
p_dot+canvas(width = 8,height = 4)

ggsave("figures/fig4/isoform_type_dotplot_top8_colour.png",
       p_dot, width = 8, height = 4, dpi = 300)
ggsave("figures/fig4/isoform_type_dotplot_top8_colour.pdf",
       p_dot, width = 8, height = 4, dpi = 300)

# remove all
rm(list=ls())
gc()
