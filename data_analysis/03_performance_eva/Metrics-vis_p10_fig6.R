#!/usr/bin/env Rscript
## -----------------------------------------------------------------
##  Minimal NBT-style plot: colour = Protocol, box fill = Group
## -----------------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(tidyr)
library(scico)
library(scales) 
library(colorspace)

## ─── 1. 设定 ------------------------------------------------------------------
metric_bg   <- c(RC="#f7fcf5", RMSE="#fee0d2", MCC="#eff3ff", FNR="#fff7bc")
metrics_vec <- names(metric_bg)

fill_group <- c(LO="#4C72B0", LS="#00A087", SO="#DD8452")
proto_levels <- c("D_ONT","M_PAB","P_ONT","I_PAB","P_ILM","P_BGI","R_ILM","R_BGI","R_ELE")
proto_shapes <- setNames(c(16,17,15,18,0,1,2,3,4), proto_levels)

## ─── 2. 读数据 & 整理 ---------------------------------------------------------
dat <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv")

plot_df <- dat |>
  filter(Metric %in% metrics_vec) |>
  mutate(
    group     = factor(group, levels = names(fill_group)),
    Metric    = factor(Metric, levels = metrics_vec),
    type      = factor(type, c("isoform","AS"), c("Isoform","AS")),
    protocols_platforms = factor(protocols_platforms, levels = proto_levels)
  )

stat_tbl <- plot_df %>%                          # 原数据框
  group_by(Metric, type, group) %>%                     # 按 Metric 和 type 分组
  summarise(
    mean = mean(value, na.rm = TRUE),            # 计算均值
    sd   = sd(value,   na.rm = TRUE),            # 计算标准差
    n    = n(),                                  # 可选：观测值个数
    .groups = "drop"
  )

stat_tbl

## ─── 3. 单面板函数 (strip 开关) ---------------------------------------------
build_panel <- function(metric, show_strip = TRUE){
  df <- filter(plot_df, Metric == metric)
  
  ggplot(df, aes(group, value)) +
    geom_boxplot(aes(fill = group),
                 width = .55, outlier.shape = NA,
                 colour = "black", linewidth = .35, alpha = .85) +
    geom_jitter(aes(shape = protocols_platforms),
                width = .15, size = 1.3, alpha = .8, colour = "grey20") +
    scale_fill_manual(values = fill_group,  name = "Group") +
    scale_shape_manual(values = proto_shapes, name = "Protocol",
                       guide  = guide_legend(nrow = 1, byrow = TRUE,
                                             override.aes = list(size = 2))) +
    
    ## 关键：strip 放到 y 轴外侧
    facet_grid(rows = vars(type),
               switch = "y",             # strip 在左
               scales = "free_y") +
    labs(x = NULL, y = NULL, title = metric) +
    
    theme_cowplot(10) %+replace%
    theme(
      panel.background   = element_blank(),
      panel.border       = element_rect(colour="black", fill=NA, linewidth=.4),
      axis.text.x        = element_text(size = 7),
      ## strip 放到最外侧
      strip.placement    = "outside",
      strip.switch.pad.grid = unit(1.5, "mm"),    # strip 与面板间距
      strip.text.y.left  = if (show_strip)        # 仅左图保留文字
        element_text(size = 8, face = "bold",
                     angle = 0, hjust = 1, vjust = 0.5,
                     margin = margin(r = 4))
      else
        element_blank(),
      axis.text.y        = element_text(size = 7),
      legend.position    = "bottom",
      legend.box         = "horizontal",
      legend.title       = element_text(size = 8),
      legend.text        = element_text(size = 7)
    )
}

## ─── 4. 生成 3 幅子图 ---------------------------------------------------------
p_a <- build_panel("RC",   show_strip = TRUE)   # 保留 strip
p_b <- build_panel("RMSE", show_strip = FALSE)  # 去除 strip
p_c <- build_panel("MCC",  show_strip = FALSE)  # 去除 strip

## ─── 5. 拼图 & 收集图例 ------------------------------------------------------
fig6abc <- (p_a | p_b | p_c) +
  plot_layout(guides = "collect", widths = c(2,2,2)) +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = "bottom")

print(fig6abc)

## 5. 输出 ----------------------------------------------------------
ggsave("figures/fig6/Fig_LO_SO_LS_box_protocolDots.pdf",
       fig6abc, width = 14, height = 12, dpi = 600, useDingbats = FALSE)

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table);  library(dplyr);   library(tidyr)
  library(forcats);     library(ggplot2); library(cowplot)
  library(RColorBrewer);library(colorspace); library(patchwork)
})

## ── 0. 读取 & 指标同向化 + z-score ────────────────────────────────
df <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv") |>
  filter(Metric %in% c("RC","RMSE","MCC")) |>
  mutate(
    value_adj = ifelse(Metric == "RMSE", -value, value),
    pipeline  = sub("-SUPPA2$", "", pipeline),
    ## 提取 protocol：若已有列用它；否则从 batch_id 解析
    protocol  = sub("^([^_]+_[^_]+).*", "\\1", batch_id)
  ) |>
  group_by(Metric) |>
  mutate(z = scale(value_adj)[, 1]) |>
  ungroup()

## ── 1. composite-z：type × group × protocol × pipeline ───────────
comp <- df |>
  group_by(type, group, protocol, pipeline) |>
  summarise(Score = mean(z, na.rm = TRUE), .groups = "drop") |>
  mutate(
    type  = factor(type,  c("isoform","AS"), c("Isoform","AS")),
    group = factor(group, c("LO","LS","SO")),
    protocol = factor(protocol)                 # 排序稍后再定
  ) |>
  group_by(pipeline) |>
  filter(any(!is.na(Score))) |>                # 去掉全 NA 管线
  ungroup()

## ── 2. y 轴顺序：LO+LS 共排，SO 单排 (按均值降序) ────────────────
ord_proto_lo_ls <- comp |>
  filter(group %in% c("LO","LS")) |>
  group_by(protocol) |>
  summarise(avg = mean(Score, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(avg)) |>
  pull(protocol)

ord_proto_so <- comp |>
  filter(group == "SO") |>
  group_by(protocol) |>
  summarise(avg = mean(Score, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(avg)) |>
  pull(protocol)

comp_lo_ls <- comp |> filter(group %in% c("LO","LS")) |>
  mutate(protocol = factor(protocol, ord_proto_lo_ls))
comp_so    <- comp |> filter(group == "SO") |>
  mutate(protocol = factor(protocol, ord_proto_so))

## ── 3. 每 facet Top-3（由 protocol 内 Score 排名） ───────────────
get_top <- function(dat){
  dat |> filter(!is.na(Score)) |>
    group_by(type, group) |>
    slice_max(Score, n = 3, with_ties = FALSE) |>
    mutate(rank = factor(row_number())) |>
    ungroup()
}
top_lo_ls <- get_top(comp_lo_ls)
top_so    <- get_top(comp_so)

rank_shapes <- c("1"=8,"2"=19,"3"=17)          # ★ ● ▲

## ── 4. 公共配色 / 主题 ───────────────────────────────────────────
cols_soft <- desaturate(brewer.pal(11,"RdBu"), .12)
fill_grad <- scale_fill_gradientn(
  colours = cols_soft, limits = c(-1.5, 1.5),
  values  = scales::rescale(c(-1.5, 0, 1.5)),
  oob     = scales::squish, name = "Composite z"
)
base_thm <- theme_cowplot(10) %+replace% theme(
  panel.border     = element_rect(colour="black", fill = NA, linewidth = .4),
  strip.placement  = "outside",
  panel.spacing    = unit(0.9, "mm"),
  axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
  axis.text.y      = element_text(size = 7)
)

## ── 5. 作图函数：x = pipeline, y = protocol ─────────────────────
bubble_plot <- function(dat, topdf,
                        show_y = TRUE, title_lab = NULL){
  ggplot(dat |> filter(!is.na(Score)),
         aes(pipeline, protocol)) +
    geom_point(aes(fill = Score, size = abs(Score)),
               shape = 21, colour = "grey25", stroke = .25) +
    geom_point(data = topdf, aes(shape = rank),
               size = 3, stroke = .6, colour = "black", fill = NA) +
    scale_size(range = c(1.2, 5), guide = "none") +
    scale_shape_manual(values = rank_shapes, name = "Top") +
    fill_grad +
    facet_grid(type ~ ., switch = "y",
               scales = "free_y", space = "free_y") +
    scale_x_discrete(drop = TRUE) +
    scale_y_discrete(drop = TRUE) +
    labs(title = title_lab, x = NULL, y = NULL) +
    base_thm +
    theme(
      plot.title   = element_text(face = "bold", hjust = .5, size = 9),
      axis.text.y  = if (show_y) element_text(size = 7) else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      legend.position = "none"
    )
}

## ── 6. 三张子图：LO/LS 共 y，SO 独立 y ─────────────────────────
p_lo <- bubble_plot(comp_lo_ls |> filter(group == "LO"),
                    top_lo_ls |> filter(group == "LO"),
                    show_y = TRUE,  title_lab = "LO")

p_ls <- bubble_plot(comp_lo_ls |> filter(group == "LS"),
                    top_lo_ls |> filter(group == "LS"),
                    show_y = FALSE, title_lab = "LS") +
  theme(strip.text.y = element_blank())

p_so <- bubble_plot(comp_so,
                    top_so,
                    show_y = TRUE,  title_lab = "SO") +
  theme(strip.text.y = element_blank())

## ── 7. 拼图 & 收集图例 ──────────────────────────────────────────
fig6def <- (p_lo | p_ls | p_so) +
  plot_layout(guides = "collect", widths = c(1.5, 1, 1.5)) &
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7)
  ) &
  guides(
    fill  = guide_colourbar(barwidth = 10, barheight = .4,
                            title.position = "top", title.hjust = .5),
    shape = guide_legend(nrow = 1, override.aes = list(size = 3))
  )

print(fig6def)

final_plot <- fig6abc / fig6def
final_plot
ggsave("figures/fig6/fig6_metrics_pipelins_Heatmap_batch_pipeline_composite.pdf",
       final_plot, width = 9, height = 12, dpi = 300, useDingbats = FALSE)

# remove all
rm(list=ls())
gc()
