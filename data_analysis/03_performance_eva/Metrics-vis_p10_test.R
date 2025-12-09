#!/usr/bin/env Rscript
## Metrics-vis_p10  •  2025-07-04
## Qingwang Chen
## ---------------------------------------------------------------------

## 0. 运行环境 ----------------------------------------------------------
setwd("/vast/projects/quartet_rna_refdata/analysis/")

custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"
.libPaths(c(custom_lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(scales)
})

## 1. 数据读取 ----------------------------------------------------------
dat <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv")

## ---------------- 1. 公共配色 & 因子水平 ------------------------
grp_levels   <- c("LO","SO","LS")
pipe_levels  <- unique(dat$pipeline)             # 或手动排序
proto_levels <- unique(dat$protocols_platforms)  # 或手动排序

box_fills <- c(LO="#E0E0E0",SO="#9ECAE1",LS="#3182BD")
pipe_cols  <- setNames(viridis::viridis(length(pipe_levels), option="C"),
                       pipe_levels)
proto_cols <- setNames(viridis::viridis(length(proto_levels), option="D"),
                       proto_levels)

dat <- dat %>%
  mutate(
    group    = factor(group,    levels = grp_levels),
    pipeline = factor(pipeline, levels = pipe_levels),
    protocols_platforms = factor(protocols_platforms, levels = proto_levels)
  )

## ---------------- 2. Metric 正/反向定义 -------------------------
## TRUE = 越大越好, FALSE = 越小越好
metric_direction <- c(
  IM   = FALSE,
  CM   = TRUE,
  RE   = TRUE,
  MRD  = FALSE,
  SCC  = TRUE,
  PET  = TRUE
)
## 若你的 Metric 名称不同，可用
# metric_direction <- setNames(rep(TRUE, length(unique(dat$Metric))),
#                              unique(dat$Metric)); metric_direction["IM"] <- FALSE

## ---------------- 3. 辅助函数 -----------------------------------
## 3.1 通用 boxplot
make_box <- function(df, metrics, facet_title, x_col, col_vals){
  df_use <- df %>% filter(Metric %in% metrics) %>%
    mutate(Metric = factor(Metric, levels = metrics))
  
  ggplot(df_use,
         aes_string(x = x_col, y = "value",
                    fill = "group", colour = x_col)) +
    geom_boxplot(width=.5, colour="black",
                 linewidth=.25, outlier.shape=NA, alpha=.7) +
    geom_jitter(width=.15, height=0, size=1.2, alpha=.85) +
    scale_fill_manual(values = box_fills, name = "Group") +
    scale_colour_manual(values = col_vals, name = NULL) +
    facet_wrap(~Metric, nrow=length(metrics), scales="free_y") +
    labs(x = NULL, y = "Value", title = facet_title) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour="grey85"),
      axis.text.x = element_text(angle=30, hjust=1, size=7),
      strip.text  = element_text(face="bold"),
      legend.position="none"
    )
}

## 3.2 barplot ±SD (cell mix / SRIV-set4 / sim)
make_bar <- function(df, metric, x_col, col_vals, ylab){
  df_sum <- df %>%
    filter(Metric==metric) %>%
    group_by(.data[[x_col]]) %>%
    summarise(mean = mean(value, na.rm=TRUE),
              sd   = sd(value,   na.rm=TRUE),
              .groups="drop")
  
  ggplot(df_sum, aes_string(x=x_col, y="mean",
                            fill=x_col)) +
    geom_col(width=.65, colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                  width=.2) +
    scale_fill_manual(values=col_vals, guide="none") +
    labs(x=NULL, y=ylab, title=metric) +
    theme_minimal(base_size=8) +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=7),
          panel.grid.major.x = element_blank())
}

## 3.3 Top-3 bubble
make_bubble <- function(df, by_col){
  ## 计算平均得分
  score_df <- df %>%
    group_by(.data[[by_col]], Metric) %>%
    summarise(score = mean(value, na.rm=TRUE), .groups="drop")
  
  ## 根据正/反向把高低转成 “好” 的方向
  score_df <- score_df %>%
    rowwise() %>%
    mutate(adj = ifelse(metric_direction[Metric], score, -score)) %>%
    ungroup()
  
  ## 排名 & 取 Top3
  score_df <- score_df %>%
    group_by(Metric) %>%
    mutate(rank = rank(-adj, ties.method = "min")) %>%
    ungroup() %>%
    filter(rank <= 3)
  
  ## 画气泡
  ggplot(score_df,
         aes_string(x="Metric", y=by_col, size="abs(score)",
                    fill="score")) +
    geom_point(shape=21, colour="black") +
    scale_fill_gradient2(low="lightgrey", mid="sandybrown",
                         high="sienna", midpoint=median(score_df$score),
                         guide=guide_colourbar(title="Score")) +
    scale_size(range=c2(3,8), guide="none") +
    geom_text(aes(label = ifelse(rank==1, "★",
                                 ifelse(rank==2,"✚","●"))),
              vjust=.5,hjust=.5,size=3) +
    theme_minimal(base_size=8) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
}

## ---------------- 4. 面板 b  (Real data, boxplot) -------------------
iso_metrics <- c("IM","CM","RE")   # adjust if needed
as_metrics  <- c("MRD","SCC","PET")

real_df <- dat %>% filter(sample_type == "Real")

p_b1 <- make_box(real_df %>% filter(type=="isoform"),
                 iso_metrics, "Quant tools", "pipeline", pipe_cols)
p_b2 <- make_box(real_df %>% filter(type=="isoform"),
                 iso_metrics, "Prot/Plat", "protocols_platforms", proto_cols)

## ---------------- 5. 面板 c (Cell mix, bar ±SD) ---------------------
cellmix_df <- dat %>% filter(sample_type == "CellMix")

p_c1 <- make_bar(cellmix_df, "MRD","pipeline",pipe_cols,"Mean")
p_c2 <- make_bar(cellmix_df, "SCC","pipeline",pipe_cols,"Mean")

## ---------------- 6. 面板 d/e (SRIV-set4 & Sim) ---------------------
set4_df <- dat %>% filter(sample_type == "SRIVset4")
sim_df  <- dat %>% filter(sample_type == "Simulation")

p_d1 <- make_bar(set4_df,"MRD","pipeline",pipe_cols,"Mean")
p_d2 <- make_bar(set4_df,"SCC","pipeline",pipe_cols,"Mean")
p_e1 <- make_bar(sim_df ,"MRD","pipeline",pipe_cols,"Mean")
p_e2 <- make_bar(sim_df ,"SCC","pipeline",pipe_cols,"Mean")

## ---------------- 7. 面板 f/g  (Top-3 bubble) -----------------------
p_f <- make_bubble(real_df %>% filter(type=="isoform"), "pipeline")
p_g <- make_bubble(real_df %>% filter(type=="isoform"), "protocols_platforms")

## ---------------- 8. patchwork 组装 -------------------------------
fig_full <- (p_b1 | p_b2) /                # row1 boxplots
  (p_c1 | p_c2) /                # row2 cell-mix
  (p_d1 | p_d2) /                # row3 set4
  (p_e1 | p_e2) /                # row4 sim
  (p_f  | p_g) +                 # row5 bubbles
  plot_layout(guides="collect") &
  theme(legend.position = "bottom")

ggsave("figures/fig6/FigX_full_composite.pdf",
       fig_full, width = 16, height = 18, dpi = 300)