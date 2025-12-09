#' Metrics-vis_p10
#' 2024-06-23
#' Qingwang Chen

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())

# Import necessary libraries
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(viridis)

# data import
all_metrics <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv")


# ───────────── 1. 仅取 isoform-level 指标 ───────────────────────────
iso_df <- all_metrics %>%
  as_tibble() %>%
  filter(type == "isoform")

# ───────────── 2. 设定因子顺序  ─────────────────────────────────────
## (1) group 三个等级按固定顺序
iso_df <- iso_df %>%
  mutate(group = factor(group, levels = c("LO", "SO", "LS"))) %>%
  mutate(Metric = factor(Metric, levels = c("RC", "RMSE", "MCC", "FNR")))%>%
  mutate(pipeline = factor(pipeline, levels = c("RSEM", "Salmon", "kallisto", "StringTie2-SO", "MPAQT", 
                                                "miniQuant", "Bambu", "IsoQuant", "Oarfish", "StringTie2-LO")))

## group 箱线颜色
box_fills <- c(LO = "#E0E0E0",
               SO = "#9ECAE1",
               LS = "#3182BD")

## 2️⃣  手动配色向量必须用同样的名字/顺序
pl_cols <- c(
  "RSEM"          = "#1f77b4",
  "Salmon"        = "#ff7f0e",
  "kallisto"      = "#2ca02c",
  "StringTie2-SO" = "#d62728",
  "MPAQT"         = "#9467bd",
  "miniQuant"     = "#8c564b",
  "Bambu"         = "#e377c2",
  "IsoQuant"      = "#7f7f7f",
  "Oarfish"       = "#bcbd22",
  "StringTie2-LO" = "#17becf"
)

## 3️⃣  画图时用 breaks / values 保持次序
p_box <- ggplot(iso_df,
                aes(x = group, y = value,
                    fill = group, colour = pipeline)) +
  
  geom_boxplot(width = .55, colour = "black",
               linewidth = .25, outlier.shape = NA, alpha = .7) +
  
  geom_jitter(width = .12, height = 0, size = 1.3, alpha = .85) +
  
  scale_fill_manual(values = box_fills, name = "Group") +
  scale_colour_manual(values = pl_cols,
                      breaks = names(pl_cols),   # ← 指定顺序
                      name   = "Pipeline") +
  
  facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
  
  labs(x = NULL, y = "Metric value",
       title = "Isoform-level performance") +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85"),
    axis.line          = element_line(colour = "black", linewidth = .55),
    axis.ticks         = element_line(colour = "black", linewidth = .55),
    strip.text         = element_text(face = "bold", size = 11),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 9, face = "bold")
  ) +
  guides(
    fill   = guide_legend(order = 1),
    colour = guide_legend(order = 2, nrow = 2, byrow = TRUE,
                          override.aes = list(size = 3))
  )

print(p_box)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig6/Fig5abcd_performance_eveluation.pdf",
       plot = p_box, width = 12, height = 4, dpi = 300)

## ───────────── 3. 改用 protocols_platforms 作为颜色分组 ─────────────
##  (1) 设定 factor 顺序（例：按手动想要的展示次序排）
proto_levels <- c("D_ONT", "M_PAB", "I_PAB",
                  "P_ONT", "P_BGI", "P_ILM",
                  "R_BGI", "R_ILM", "R_ELE")   # ← 按你自己的需要排列
iso_df <- iso_df %>%
  mutate(protocols_platforms = factor(protocols_platforms,
                                      levels = proto_levels))

##  (2) 给每个平台手动或自动配色
# 若想自动配色，可用 viridis 调色：
library(viridis)
proto_cols <- setNames(viridis(length(proto_levels), option = "plasma"),
                       proto_levels)

# 也可以手动指定：
# proto_cols <- c(
#   "P_BGI" = "#1b9e77", "P_NVG" = "#d95f02", "P_GRM" = "#7570b3",
#   "D_BGI" = "#e7298a", "D_NVG" = "#66a61e", "D_GRM" = "#e6ab02",
#   "M_BGI" = "#a6761d", "M_NVG" = "#666666", "M_GRM" = "#1f78b4"
# )

##  (3) 绘图
p_box_proto <- ggplot(iso_df,
                      aes(x = protocols_platforms, y = value,
                          fill = group, colour = protocols_platforms)) +
  
  geom_boxplot(width = .55, colour = "black",
               linewidth = .25, outlier.shape = NA, alpha = .7) +
  # geom_jitter(width = .12, height = 0, size = 1.3, alpha = .85) +
  
  scale_fill_manual(values = box_fills, name = "Group") +
  scale_colour_manual(values = proto_cols,
                      breaks = names(proto_cols),
                      name   = "Protocol / Platform") +
  
  facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
  labs(x = "Protocol / Platform", y = "Metric value",
       title = "Isoform-level performance by protocol / platform") +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85"),
    axis.line          = element_line(colour = "black", linewidth = .55),
    axis.ticks         = element_line(colour = "black", linewidth = .55),
    strip.text         = element_text(face = "bold", size = 11),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 9, face = "bold"),
    axis.text.x        = element_text(angle = 30, hjust = 1)
  ) +
  guides(
    fill   = guide_legend(order = 1),
    colour = guide_legend(order = 2, nrow = 2, byrow = TRUE,
                          override.aes = list(size = 3))
  )

print(p_box_proto)

## ——— 保存 ————————————————————————————
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig6/Fig5abcd_performance_by_protocol.pdf",
       plot = p_box_proto, width = 12, height = 4, dpi = 300)


# remove all
rm(list=ls())
gc()
