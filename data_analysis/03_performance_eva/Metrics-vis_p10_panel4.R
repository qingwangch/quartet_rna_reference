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

## 1. 基础设定 ------------------------------------------------------
metric_bg  <- c(RC="#f7fcf5", RMSE="#fee0d2",
                MCC="#eff3ff", FNR="#fff7bc")
metrics_vec <- names(metric_bg)

fill_group <- c(LO="#9ECAE1", LS="#08519c", SO="#3182BD")   # 箱体填充
grp_levels <- names(fill_group)

proto_levels <- c("D_ONT","M_PAB","P_ONT","I_PAB",
                  "P_ILM","P_BGI","R_ILM","R_BGI","R_ELE")
proto_cols   <- setNames(brewer.pal(9, "Set1"), proto_levels)

## 2. 读数据 --------------------------------------------------------
dat <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv")

plot_df <- dat %>%
  filter(Metric %in% metrics_vec) %>%
  mutate(
    group     = factor(group, levels = grp_levels),
    Metric    = factor(Metric, levels = metrics_vec),
    type      = factor(type, levels = c("isoform","AS"),
                       labels = c("Isoform","AS")),
    protocols_platforms = factor(protocols_platforms, levels = proto_levels)
  )

## 3. 单面板函数 ----------------------------------------------------
build_panel <- function(metric){
  
  df <- plot_df %>% filter(Metric == metric)
  
  ggplot(df, aes(x = group, y = value)) +
    
    # 箱体
    geom_boxplot(aes(fill = group),
                 width = .55, outlier.shape = NA,
                 colour = "black", size = .35, alpha = .8) +
    
    # 散点
    geom_jitter(aes(colour = protocols_platforms),
                width = .15, size = 1.2, alpha = .85, stroke = 0) +
    
    scale_fill_manual(values = fill_group,  name = "Group") +
    scale_colour_manual(values = proto_cols, name = "Protocol") +
    
    facet_grid(rows = vars(type), switch = "y", scales = "free_y") +
    labs(x = NULL, y = NULL, title = metric) +
    
    theme_cowplot(9) %+replace%
    theme(
      panel.background = element_rect(fill = metric_bg[[metric]], colour = NA),
      panel.border     = element_rect(colour="black", fill=NA, linewidth=.4),
      axis.text.x      = element_text(size = 7),
      strip.text.y     = element_text(face = "bold", size = 8),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7)
    )
}

## 4. 生成四幅子图 --------------------------------------------------
p_a <- build_panel("RC")
p_b <- build_panel("RMSE")
p_c <- build_panel("MCC")

fig6abc <- (p_a | p_b | p_c ) +
  plot_annotation(tag_levels='a') +     # 自动 a–f 标面板
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7)
  ) &
  guides(
    fill   = guide_legend(title = "Group",
                          nrow = 1, byrow = TRUE, order = 1),
    colour = guide_legend(title = "Protocol",
                          nrow = 1, byrow = TRUE, order = 2)
  )
print(fig6abc)

## 5. 输出 ----------------------------------------------------------
ggsave("figures/fig6/Fig_LO_SO_LS_box_protocolDots.pdf",
       fig6abc, width = 14, height = 12, dpi = 600, useDingbats = FALSE)


# best practice
## ---------- 0. 读 & 方向一致 + z-score --------------------------
df <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv")

df <- df %>% 
  filter(Metric %in% c("RC","RMSE","MCC")) %>%
  mutate(value_adj = ifelse(Metric=="RMSE", -value, value)) %>%
  group_by(Metric) %>% 
  mutate(z = (value_adj - mean(value_adj, na.rm=TRUE)) / sd(value_adj, na.rm=TRUE)) %>% 
  ungroup()

## ---------- 1. 综合得分 -----------------------------------------
comp <- df %>% 
  group_by(type, group, batch_id, pipeline) %>% 
  summarise(Score = mean(z, na.rm=TRUE), .groups="drop") %>%
  mutate(
    group = factor(group, levels=c("LO","SO","LS")),
    type  = factor(type, levels=c("isoform","AS"), labels=c("Isoform","AS"))
  )

pipe_levels  <- sort(unique(comp$pipeline))
batch_levels <- sort(unique(comp$batch_id))

comp <- comp %>%
  mutate(pipeline = factor(pipeline, levels=pipe_levels),
         batch_id = factor(batch_id, levels=batch_levels))

## ---------- 2. 单面板函数 ----------------------------------------
rank_shapes <- c("1"=8,"2"=19,"3"=17)    # ★ ● ▲
## --- 自定义柔和 RdBu -------------------------------------------------
base_cols <- RColorBrewer::brewer.pal(11, "RdBu")

# 在两端加一点灰度，让过渡更柔和
soft_cols <- desaturate(base_cols, amount = 0.1)

## 全局极值
rng <- range(comp$Score, na.rm = TRUE)  # comp 是综合得分表

# p_fill <- scale_fill_gradientn(
#   colours   = soft_cols,
#   values    = scales::rescale(c(rng[1], 0, rng[2])),   # 保证 0 在中点
#   limits    = rng,
#   name      = "Composite z",
#   na.value  = "white",
#   guide     = guide_colorbar(barwidth = unit(0.35, "npc"),
#                              barheight = 0.4,
#                              title.position = "top",
#                              title.hjust = 0.5)
# )

p_fill <- scale_fill_gradientn(
  colours = soft_cols,
  limits  = c(-1.5, 1.5),
  oob     = scales::squish,
  values  = scales::rescale(c(-1.5, 0, 1.5)),
  name    = "Composite z"
)
# cols <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
# 
# p_fill <- scale_fill_stepsn(colours = cols, limits = rng, n.breaks = 9,
#                             name = "Composite z")
                           
make_heat <- function(df_one){
  df_one <- df_one %>% droplevels() %>% complete(batch_id, pipeline, fill=list(Score=NA))
  
  top3 <- df_one %>% filter(!is.na(Score)) %>% 
    arrange(desc(Score)) %>% slice_head(n=3) %>% 
    mutate(rank=factor(row_number()))
  
  ggplot(df_one, aes(pipeline, batch_id, fill=Score)) +
    geom_tile(colour="grey90") +
    p_fill +
    geom_point(data=top3, aes(shape=rank), colour="black", size=2.6, stroke=.35) +
    scale_shape_manual(values=rank_shapes, name="Top",
                       labels=c("1st","2nd","3rd"), breaks=levels(top3$rank)) +
    theme_cowplot(9) %+replace%
    theme(
      axis.text.x  = element_text(angle=45, hjust=1, vjust=1, size=6.8),
      axis.text.y  = element_text(size=6.8),
      panel.border = element_blank()
    ) +
    labs(x=NULL, y=NULL)
}

## ---------- 3. 生成 2×3 网格 ------------------------------------
# helper to drop/keep x-axis labels
drop_x <- theme(axis.text.x = element_blank())

iso_lo <- make_heat(comp %>% filter(type=="Isoform", group=="LO")) + drop_x
iso_so <- make_heat(comp %>% filter(type=="Isoform", group=="SO")) + drop_x
iso_ls <- make_heat(comp %>% filter(type=="Isoform", group=="LS")) + drop_x

as_lo  <- make_heat(comp %>% filter(type=="AS",      group=="LO"))
as_so  <- make_heat(comp %>% filter(type=="AS",      group=="SO"))
as_ls  <- make_heat(comp %>% filter(type=="AS",      group=="LS")) +
  labs(x="Pipeline")        # 仅右下角显示轴标题

row_iso <- iso_lo | iso_ls | iso_so
row_as  <- as_lo  | as_ls  | as_so

fig6defg <- (row_iso / row_as) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box      = "vertical",
    legend.margin   = margin(t = 0, unit = "pt"),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7)
  ) &
  guides(
    fill   = guide_colourbar(barwidth = 10, barheight = .4, order=1,
                             title.position="top", title.hjust=.5),
    shape  = guide_legend(nrow = 1, order=2),
    colour = guide_legend(nrow = 1, order=3),
    fill_new  = "collect"  # patchwork 1.2+
  )
fig6defg
ggsave("figures/fig6/Heatmap_batch_pipeline_composite.pdf",
       fig6defg, width = 9, height = 7, dpi = 600, useDingbats = FALSE)

final_plot <- fig6abc / fig6defg
final_plot
ggsave("figures/fig6/fig6_metrics_pipelins_Heatmap_batch_pipeline_composite.pdf",
       final_plot, width = 9, height = 12, dpi = 300, useDingbats = FALSE)

# remove all
rm(list=ls())
gc()
