# num_quartet_assembly_isoform_as.R
# 2025-07-10
# Qingwang Chen

# 工作目录 & R 库路径
setwd("/vast/projects/quartet_rna_refdata/analysis/")
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"
.libPaths(c(custom_lib_path, .libPaths()))

# 加载必要的包
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(purrr)
library(tidytext) 
library(scales)
library(stringr) 
library(ggview) 
library(patchwork)

# data import
isoform_num <- fread("/vast/projects/quartet_rna_refdata/slurm/annotation/gtf_counts.csv")

as_num <- fread("/vast/projects/quartet_rna_refdata/slurm/suppa2/events_generation/ioe_event_stats.csv")

## ---------- 1. panel-e  Genes / Isoforms ---------------------------------
ds_order <- c("GENCODE", "GENCODE+Quartets", "GENCODE+Quartets_filter")

## ❶ 手动色板：深色＝Annotated；同色 25% 透明度＝Novel
pal_ds <- c("GENCODE+Quartets"        = "#4E79A7",
            "GENCODE+Quartets_filter" = "#9F5FE6",
            "GENCODE"                 = "#E15759")

## ---- 整理长表 ----
iso_long <- isoform_num %>% 
  ## ① 先把 File 重编码成目标字符串，再设定 factor 顺序
  mutate(dataset = recode(File,
                          "GENCODE_plus_Quartets"        = "GENCODE+Quartets",
                          "GENCODE_plus_Quartets_filter" = "GENCODE+Quartets_filter",
                          "GENCODE"                      = "GENCODE"),
         dataset = factor(dataset, levels = ds_order)) %>% 
  ## ② 只保留 Ann_/Novel_ 列
  select(dataset, matches("^(Ann|Novel)_")) %>%                   
  pivot_longer(-dataset,
               names_to      = c("state", "feature"),
               names_sep     = "_",
               values_to     = "value") %>% 
  mutate(feature = recode(feature,
                          Gene = "Genes",
                          Tx   = "Isoforms"),
         state   = factor(state,
                          levels = c("Ann","Novel"),
                          labels = c("Annotated","Novel")),
         alpha   = ifelse(state == "Annotated", 1, .25))

## ----- 画图：同数据集纵向堆叠，数据集之间横向 dodge -----
## state 的因子顺序决定堆叠先后：Annotated→Novel
iso_long <- iso_long %>%                     # 保持因子顺序
  mutate(state = factor(state,
                        levels = c("Annotated", "Novel")))

pe <- ggplot(iso_long,
             aes(x = dataset, y = value,
                 fill  = dataset,
                 alpha = state)) +
  geom_col(width = .55, colour = "black", linewidth = .25,
           position = position_stack(reverse = TRUE)) +  # ← 反转
  
  facet_grid(. ~ feature, scales = "free_x", space = "free_x") +
  
  scale_fill_manual(values = pal_ds, name = NULL) +
  scale_alpha_manual(values = c(Annotated = 1, Novel = .25),
                     name   = NULL,
                     labels = c("annotated", "novel")) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, .05))) +
  
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.ticks   = element_line(colour = "black"),
    axis.line    = element_line(colour = "black"),
    axis.text.x  = element_blank(),
    strip.text   = element_blank(),
    strip.background = element_blank(),
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.box.just      = "left",
    legend.margin        = margin(2, 2, 2, 2),
    legend.text          = element_text(size = 10),
    plot.margin          = margin(5.5, 14, 5.5, 5.5)
  ) +
  guides(
    alpha = guide_legend(order = 2,
                         override.aes = list(fill = "grey40", colour = "black")),
    fill  = guide_legend(order = 1)
  )

pe
pe+canvas(width = 3,height = 4)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4e_num_iso.pdf", plot = pe, width = 3, height = 4)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4e_num_iso.png", plot = pe, width = 3, height = 4, dpi = 300)

## ---------- 2. panel-g  Splicing events ----------------------------------
evt_keep <- c("A3","A5","AF","AL","MX","SE","RI")
evt_lab  <- c(A3 = "alternate 3' splice",
              A5 = "alternate 5' splice",
              AF = "alternate first exon",
              AL = "alternate last exon",
              MX = "mutually exclusive exons",
              SE = "exon skipping",
              RI = "intron retention")

as_long <- as_num %>%
  mutate(dataset = recode(file,
                          "GENCODE_plus_Quartets"        = "GENCODE+Quartets",
                          "GENCODE_plus_Quartets_filter" = "GENCODE+Quartets_filter",
                          "GENCODE"                      = "GENCODE"),
         dataset = factor(dataset, levels = ds_order)) %>%
  select(dataset, all_of(evt_keep)) %>%
  pivot_longer(-dataset,
               names_to = "event", values_to = "count") %>%
  mutate(event = factor(event, levels = evt_keep, labels = evt_lab))

pg <- ggplot(as_long,
             aes(event, count,
                 fill = dataset)) +
  geom_col(position = position_dodge(width = .7),
           width = .6, colour = "black", linewidth = .25) +
  scale_fill_manual(values = pal_ds, name = NULL) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = NULL, title = "Splicing events") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x  = element_text(angle = 30, hjust = 1),
    axis.ticks   = element_line(colour = "black"),
    axis.line    = element_line(colour = "black"),
    legend.position = "none",
    plot.title   = element_text(face = "bold")
  )
pg
pg+canvas(width = 3,height = 4)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4g_num_as.pdf", plot = pg, width = 3, height = 4)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4g_num_as.png", plot = pg, width = 3, height = 4, dpi = 300)

# ---------------------------------------------------------------
# ‼ 请确保 panel_ABCD、pe、p、pg 三个对象仍在环境中
#   • panel_ABCD : 你 earlier 拼好的 a–d
#   • pe         : Genes / Isoforms 叠柱图   → e
#   • p          : novel-isoform 分类柱图   → f
#   • pg         : splicing-event 柱图      → g
# ---------------------------------------------------------------

## ---- e  f  g  横排 --------------------------------------------------------
efg_row <- (pe | p | pg) +                             # 3 列
  plot_layout(widths = c(1, 1.35, 1.35)) +             # 不再 guides='collect'
  plot_annotation(tag_levels = list(letters[5:7]),
                  theme = theme(plot.tag = element_text(face = "bold",
                                                        size = 12,
                                                        vjust = 1.2))) &
  theme(legend.position = "none")      # ★ 这一行将子图全部隐藏图例
efg_row
efg_row+canvas(width = 12,height = 4)

## ===============================================================
## 7.  保存：180 mm (=7.1 in) 宽，160 mm 高；300 dpi
## ===============================================================
out_dir <- "/vast/projects/quartet_rna_refdata/analysis/figures/fig4/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "Fig4efg_final.pdf"),
       plot   = efg_row,
       width  = 12,   # ≈ 180 mm
       height = 4,   # NBT 单栏建议高度
       dpi = 300)

ggsave(file.path(out_dir, "Fig4efg_final.png"),
       plot   = efg_row,
       width  = 12,
       height = 4,
       dpi    = 300)

message("✔ Fig. 4 saved to: ", out_dir, "Fig4efg_final.[pdf|png]")

# remove all
rm(list=ls())
gc()
