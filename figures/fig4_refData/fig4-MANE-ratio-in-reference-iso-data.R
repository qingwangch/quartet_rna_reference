# MANE-ratio-in-reference-iso-data.R
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

# MANE
MANE_summary <- fread("/vast/projects/quartet_rna_refdata/analysis/MANE/MANE.GRCh38.v1.4.summary.txt") %>% as.data.frame()
MANE_isoforms <- MANE_summary$Ensembl_nuc

## ── 3. 计算每个 sample-pair 的 MANE 覆盖率 ────────────────────────────
# ① 标记 refFC 里哪些 isoform 属于 MANE
ref_fc_mane <- ref_fc %>%
  mutate(
    is_MANE = isoform %in% MANE_isoforms,          # TRUE / FALSE
    compare = factor(compare,                      # 固定显示顺序
                     levels = c("D5/D6","F7/D6","M8/D6"))
  )

# ② 统计：总数、MANE 数、百分比
mane_stat <- ref_fc_mane %>%
  group_by(compare) %>%
  summarise(
    n_total = n(),
    n_MANE  = sum(is_MANE),
    pct     = 100 * n_MANE / n_total,
    .groups = "drop"
  )

print(mane_stat)
#> # A tibble: 3 × 4    compare n_total n_MANE   pct
#>   <fct>      <int>   <int> <dbl>
#> 1 D5/D6        ...     ...  ...
#> 2 F7/D6        ...     ...  ...
#> 3 M8/D6        ...     ...  ...

## ── 4. 绘图  ─────────────────────────────────────────────────────────
nbt_cols <- c("D5/D6" = "#1b9e77",
              "F7/D6" = "#d95f02",
              "M8/D6" = "#7570b3")

p_mane <- ggplot(mane_stat,
                 aes(x = compare, y = pct, fill = compare)) +
  geom_col(width = .55, colour = "black", alpha = .85) +
  geom_text(aes(label = paste0(round(pct,1), "%\n(",
                               n_MANE, "/", n_total, ")")),
            vjust = -0.25, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = nbt_cols) +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, .05))) +
  labs(title = "MANE isoform validation in refFC (LO)",
       subtitle = "Percentage of MANE transcripts within each sample pair",
       x = "Sample pair", y = "MANE validation rate (%)") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.ticks.x = element_line(size = 0.55),
    plot.title   = element_text(face = "bold", size = 14, hjust = .5),
    plot.subtitle= element_text(hjust = .5)
  )

print(p_mane)
p_mane+canvas(width = 5,height = 4)

## ── 5. 保存  ────────────────────────────────────────────────────────
ggsave("figures/fig4/refFC_MANE_rate_bar.png",
       p_mane, width = 5, height = 4, dpi = 300)
ggsave("figures/fig4/refFC_MANE_rate_bar.pdf",
       p_mane, width = 5, height = 4, dpi = 300)

# remove all
rm(list=ls())
gc()


