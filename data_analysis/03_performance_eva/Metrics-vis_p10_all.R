#!/usr/bin/env Rscript
## ---------------------------------------------------------------
##  Panel: fill = protocol  • colour = pipeline  • shape = compare
## ---------------------------------------------------------------
library(data.table); library(dplyr)
library(ggplot2);    library(cowplot); library(patchwork)

## 1. 颜色与顺序 -----------------------------------------------------
metric_bg   <- c(RC="#f7fcf5", RMSE="#fee0d2",
                 MCC="#eff3ff", FNR="#fff7bc")

grp_levels  <- c("LO","SO","LS")
metrics_vec <- names(metric_bg)

#── Pipeline ∶ 用 Okabe-Ito (色盲友好，8 色)，不足再叠加灰调
pipe_levels <- c("miniQuant","MPAQT","Bambu","IsoQuant","Oarfish",
                 "StringTie2-LO","RSEM","Salmon","kallisto","StringTie2-SO",
                 "miniQuant-SUPPA2","MPAQT-SUPPA2","Bambu-SUPPA2","IsoQuant-SUPPA2",
                 "Oarfish-SUPPA2","StringTie2-LO-SUPPA2","RSEM-SUPPA2",
                 "Salmon-SUPPA2","kallisto-SUPPA2","StringTie2-SO-SUPPA2")

okabe8 <- c("#E69F00","#56B4E9","#009E73","#F0E442",
            "#0072B2","#D55E00","#CC79A7","#000000")
pipe_cols <- setNames(c(okabe8, scales::alpha(okabe8, .6))[1:length(pipe_levels)],
                      pipe_levels)

#── Protocol ∶ 9 种；箱体填充
proto_levels <- c("D_ONT","M_PAB","P_ONT","I_PAB",
                  "P_ILM","P_BGI","R_ILM","R_BGI","R_ELE")
proto_fills  <- setNames(RColorBrewer::brewer.pal(9,"Set3")[c(2,3,4,5,6,7,8,9,1)],
                         proto_levels)

#── Compare ∶ 点形状
compare_shapes <- c("D5/D6" = 21, "F7/D6" = 22, "M8/D6" = 24)  # 按实际水平扩充

## 2. 读数据 ---------------------------------------------------------
dat <- fread("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv")

plot_df <- dat %>%
  filter(Metric %in% metrics_vec) %>%
  mutate(
    group     = factor(group, levels = grp_levels),
    Metric    = factor(Metric, levels = metrics_vec),
    type      = factor(type, levels = c("isoform","AS"), labels = c("Isoform","AS")),
    pipeline  = factor(pipeline,           levels = pipe_levels),
    protocols_platforms = factor(protocols_platforms, levels = proto_levels),
    compare   = factor(compare,            levels = names(compare_shapes))
  )

## 3. 单图函数 --------------------------------------------------------
panel_fun <- function(metric){
  
  df <- plot_df %>% filter(Metric == metric)
  
  ggplot(df, aes(x = group, y = value,
                 fill   = protocols_platforms,
                 colour = pipeline,
                 shape  = compare)) +
    
    geom_boxplot(width=.55, size=.35, outlier.shape = NA) +
    geom_jitter(size=1.2, width=.15, stroke=.25, alpha=.85) +
    
    scale_fill_manual(values = proto_fills,  name = "Protocol") +
    scale_colour_manual(values = pipe_cols,  name = "Pipeline") +
    scale_shape_manual(values = compare_shapes, name = "Compare") +
    
    facet_grid(rows = vars(type), switch = "y", scales = "free_y") +
    labs(x = NULL, y = NULL, title = metric) +
    
    theme_cowplot(9) %+replace%
    theme(
      panel.background = element_rect(fill = metric_bg[[metric]], colour = NA),
      panel.border     = element_rect(colour="black", fill=NA, linewidth=.4),
      axis.text.x      = element_text(size=7),
      strip.text.y     = element_text(face="bold", size=8),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.title     = element_text(size=8),
      legend.text      = element_text(size=7)
    )
}

## 4. 拼四图 ----------------------------------------------------------
p_a <- panel_fun("RC");  p_b <- panel_fun("RMSE")
p_c <- panel_fun("MCC"); p_d <- panel_fun("FNR")

final_plot <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect")  &
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7))
final_plot


ggsave("figures/fig6/Fig_metrics_fillProtocol_colourPipeline_shapeCompare.pdf",
       final_plot, width = 14, height = 12, dpi = 600, useDingbats = FALSE)

# remove all
rm(list=ls())
gc()
