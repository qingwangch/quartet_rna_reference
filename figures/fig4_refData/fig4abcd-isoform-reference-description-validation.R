## ============================================================================
##  Fig.4  —  Isoform reference (LO)   │   Qingwang Chen · 2025-06-13
##  ─ Panel A  content  : total count + type-composition
##  ─ Panel B  reliability : LO vs qPCR scatter (r + RMSE + N)
##  ─ Panel C  clinical relevance : MANE coverage %
## ============================================================================

## 0. 环境 --------------------------------------------------------------------
setwd("/vast/projects/quartet_rna_refdata/analysis/")
.libPaths(c("/vast/projects/quartet_rna_refdata/my_r_packages", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(data.table); library(ggplot2); library(cowplot)
  library(forcats);  library(purrr);  library(patchwork)
})

nbt_cols <- c("D5/D6"="#1b9e77","F7/D6"="#d95f02","M8/D6"="#7570b3")

## 1. 读取 refFC（LO） ---------------------------------------------------------
ref_fc  <- fread("isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv")
ref_fc  <- mutate(ref_fc, compare = factor(compare, levels=c("D5/D6","F7/D6","M8/D6")))
count <- dplyr::count
transcrpt_type <- ref_fc[,c("isoform","transcript_type")] %>% distinct()
## 1A. 总量 + DEI  -------------------------------------------------------------
## refDEI (LO) 仅用于计数
ref_dei <- fread("isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv") |>
  filter(Final!="non-DEI") |>
  mutate(compare=factor(compare,levels=levels(ref_fc$compare)))
ref_dei <- merge(ref_dei,transcrpt_type)

## ─── A. Total-count panel  --------------------------------------------------
total_df <- bind_rows(
  dplyr::count(ref_fc , compare,  name = "n") |> mutate(cat = "refFC"),
  dplyr::count(ref_dei, compare,  name = "n") |> mutate(cat = "refDEI")
) |>
  ## ❶ 设定 facet 顺序：refFC ➜ refDEI
  mutate(cat = factor(cat, levels = c("refFC", "refDEI")))

p_total <- ggplot(total_df,
                  aes(y = compare,          # 横向柱形
                      x = n,
                      fill = compare)) +
  geom_col(width = .55, colour = "black", alpha = .85,
           show.legend = FALSE) +
  scale_fill_manual(values = nbt_cols) +
  ## ❷ 纵向排列：ncol = 1  (而不是 nrow = 1)
  facet_wrap(~ cat, ncol = 1,
             labeller = labeller(
               cat = c(refFC  = "refFC isoform",
                       refDEI = "refDEI"))) +
  scale_x_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = "", y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = .8),
    axis.text.y  = element_text(size = 11, face = "bold"),
    axis.ticks.y = element_line(size = .55),
    axis.ticks.x = element_line(size = .55),
    axis.line    = element_line(size = .6)
  )

p_total

## -----------------------------------------------------------
## 1B. transcript-type 组成（Top-8 + Other, facetted by cat）
## -----------------------------------------------------------

## 「Top-6」仍按 refFC 频次来定
top8 <- ref_fc %>%
  count(transcript_type, sort = TRUE) %>%
  slice_head(n = 6) %>% pull(transcript_type)

## ❶ 先把 refFC / refDEI 合并，加上 cat 标签 -----------------
type_df <- bind_rows(
  ref_fc  %>% mutate(cat = "refFC isoform"),
  ref_dei %>% mutate(cat = "refDEI")
) %>%
  mutate(
    cat  = factor(cat, levels = c("refFC isoform", "refDEI")),
    compare = factor(compare,
                     levels = c("D5/D6","F7/D6","M8/D6")),
    type_raw = fct_other(transcript_type,
                         keep = top8,
                         other_level = "Other")
  ) %>%
  count(cat, compare, type_raw, name = "n") %>%
  ## 为了让「x 轴」在两个 facet 中顺序一致
  group_by(type_raw) %>%               # 全局总和排序
  mutate(total_n = sum(n, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(type_raw = fct_reorder(type_raw, total_n, .desc = TRUE))

## ❷ 颜色字典（同之前，只要覆盖 top8+Other 即可） ----------
type_cols <- c(
  "protein_coding"               = "#FFD92F",
  "lncRNA"                       = "#4DAF4A",
  "retained_intron"              = "#E41A1C",
  "nonsense_mediated_decay"      = "#984EA3",
  "processed_transcript"         = "#E7298A",
  "protein_coding_CDS_not_defined"= "#66A61E",
  "IG_C_gene"                    = "#E6AB02",
  "TEC"                          = "#A6761D",
  "Other"                        = "#666666"
)

## ❸ dot-plot ------------------------------------------------
# p_type <- ggplot(type_df,
#                  aes(x = type_raw, y = compare,
#                      size = n, fill = type_raw, colour = type_raw)) +
#   geom_point(shape = 21, stroke = .35, colour = "black") +
#   geom_text(aes(label = n),
#             colour = "white", size = 2.7, fontface = "bold") +
#   scale_size(range = c(3.8, 15), guide = "none") +
#   scale_fill_manual(values = type_cols, guide = "none") +
#   scale_colour_manual(values = type_cols, guide = "none") +
#   facet_wrap(~ cat, ncol = 1,                # ↙ 跟 p_total 一样纵向 facet
#              labeller = labeller(
#                cat = c("refFC isoform" = "refFC isoform",
#                        "refDEI"       = "refDEI")
#              )) +
#   labs(x = "Transcript type (Top-6 + Other)",
#        y = NULL) +
#   theme_classic(base_size = 12) +
#   theme(
#     strip.text      = element_text(face = "bold", size = 11),
#     strip.background = element_rect(fill = "white", colour = "black", linewidth = .8),
#     axis.text.x     = element_text(angle = 20, hjust = 1, size = 9.5),
#     axis.text.y     = element_blank(),
#     axis.ticks.y    = element_blank(),
#     axis.ticks.x    = element_line(),
#     axis.line.x     = element_line()
#   )
# p_type

## 连续型配色 —— 这里用蓝紫色渐变，可改成 scale_fill_gradient(low = "...", high = "...")
hm_pal <- scales::seq_gradient_pal(low = "#d6e8fa", high = "#1a1e91", space = "Lab")

p_heat <- ggplot(type_df,
                 aes(x = type_raw,           # 横轴：转录本类型
                     y = compare,            # 纵轴：比较类别
                     fill = n)) +            # 颜色映射到计数 / 频率
  geom_tile(colour = "grey80", linewidth = .3) +   # 画色块
  geom_text(aes(label = n),                    # 叠加数字
            colour = "black", size = 3, fontface = "bold") +
  scale_fill_gradientn(
    colours = hm_pal(seq(0, 1, length.out = 100)),
    limits  = c(0, max(type_df$n, na.rm = TRUE)),  # 颜色范围
    name    = "Frequency (%)"
  ) +
  facet_wrap(~ cat, ncol = 1,
             labeller = labeller(
               cat = c("refFC isoform" = "refFC isoform",
                       "refDEI"       = "refDEI")
             )) +
  labs(x = "", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = .8),
    axis.text.x     = element_text(angle = 25, hjust = 1, size = 9.5),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.line.x     = element_line(),
    legend.position = "right"
  )+
  scale_fill_gradientn(
    colours = c("#e0f3ff", "#64b5f6", "#00796B"),  # 浅蓝 → 天蓝 → 深绿
    limits  = c(0, max(type_df$n)),
    name    = "Frequency"
  )

print(p_heat)

## -----------------------------------------------------------
## Panel AB : p_total (左)  +  p_type (右)
## -----------------------------------------------------------
panel_AB <- plot_grid(
  # 左：a
  p_total + theme(plot.margin = margin(4, 4, 4, 4)),
  # 右：b
  p_heat  + theme(plot.margin = margin(4, 4, 4, 4)),
  nrow       = 1,              # 并排
  rel_widths = c(1.2, 2),      # 适当调左右宽比
  align      = "h",            # y 轴对齐
  axis       = "tb",           # 同步上下须线
  labels     = c("a", "b"),    # 直接打 a / b
  label_size = 14,
  label_fontface = "bold"
)

# 预览
panel_AB


## 2. Panel C  —  LO vs qPCR scatter ------------------------------------------
qpcr_iso <- fread("isoforms/qpcr/log2fc_20250409_U.csv")

setkey(ref_fc,  isoform_compare, compare)
setkey(qpcr_iso, isoform_compare, compare)

scat_dt <- merge(
  ref_fc[, .(isoform_compare, compare, fc_ref = log2FC)],
  qpcr_iso[,.(isoform_compare, compare, fc_val = log2FC)]
)

## 若想保持原来的比较顺序
scat_dt[, compare := factor(compare, levels = unique(ref_fc$compare))]

##------------ 2. 统计表 --------------------------------------------------------
stat_dt <- scat_dt[, .(
  N    = .N,
  r    = cor(fc_ref, fc_val, use = "complete.obs"),
  RMSE = sqrt(mean((fc_val - fc_ref)^2, na.rm = TRUE))
), by = compare]

##------------ 3. 颜色 & 坐标范围设定 -----------------------------------------
nbt_cols <- c("D5/D6"="#1b9e77","F7/D6"="#d95f02","M8/D6"="#7570b3")

rng <- range(unlist(scat_dt[, .(fc_ref, fc_val)]), na.rm = TRUE)
lim <- range(rng) + c(-1, 1) * diff(rng) * 0.05   # 5% padding

##------------ 4. 绘图 ----------------------------------------------------------
p_scatter <- ggplot(scat_dt, aes(fc_ref, fc_val, colour = compare)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = .4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80", linewidth = .3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80", linewidth = .3) +
  geom_point(alpha = .35, size = 1.6, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = .45) +
  facet_wrap(~compare, nrow = 1) +
  geom_text(
    data = stat_dt,
    aes(x = -Inf, y =  Inf,
        label = sprintf("r = %.2f\nRMSE = %.2f\nN = %d", r, RMSE, N)),
    hjust = -0.05, vjust = 1.10, size = 3.4, fontface = "bold"
  ) +
  scale_colour_manual(values = nbt_cols) +
  labs(
    x = expression(log[2]*"FC in reference"),
    y = expression(log[2]*"FC by qPCR")
  ) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = .6),
    axis.ticks.length = unit(0.18, "cm"),
    axis.ticks        = element_line(colour = "black", linewidth = .5),
    strip.text        = element_text(size = 12, face = "bold", margin = margin(b = 4)),
    strip.background  = element_blank(),
    legend.position = "none"
  )

##------------ 5. Patchwork（如需后续拼图） ------------------------------------
panel_C <- p_scatter + plot_layout(guides = "collect") &
  theme(plot.margin = margin(4, 4, 4, 4))

panel_C          # 输出

## 3. Panel D  —  MANE 覆盖率 ---------------------------------------------------
mane_vec <- fread("MANE/MANE.GRCh38.v1.4.summary.txt")$Ensembl_nuc

mane_stat <- ref_fc |>
  mutate(compare=factor(compare,levels=levels(ref_fc$compare)),
         is_MANE=isoform %in% mane_vec) |>
  group_by(compare) |>
  summarise(n_total=n(),
            n_MANE=sum(is_MANE),
            pct=100*n_MANE/n_total,.groups="drop")

p_mane <- ggplot(mane_stat,
                 aes(compare,pct,fill=compare))+
  geom_col(width=.55,color="black",alpha=.9,show.legend=FALSE)+
  geom_text(aes(label=sprintf("%.1f%%\n(%d/%d)",pct,n_MANE,n_total)),
            vjust=-0.25,size=3.4,fontface="bold")+
  scale_fill_manual(values=nbt_cols)+
  scale_y_continuous(limits=c(0,105),expand=expansion(mult=c(0,.05)))+
  labs(x=NULL,y="MANE validation rate (%)")+
  theme_classic(base_size=13)+
  theme(axis.text.x=element_text(size=11,face="bold"),
        axis.ticks.x=element_blank())

panel_D <- p_mane + theme(plot.margin=margin(4,4,4,4))
panel_D

## 4. 组合整图  ---------------------------------------------------------------
## ── 1. 先各自去掉内部标签（如果之前加过 cowplot 标签的话） -------------
p_total_nolab   <- p_total   + labs(tag = NULL)
# p_type_nolab    <- p_type    + labs(tag = NULL)
p_heat
p_scatter_nolab <- p_scatter + labs(tag = NULL)   # 你的对象名若是 p_iso_scatter，请替换
p_mane_nolab    <- p_mane    + labs(tag = NULL)

## ── 2. 拼图：第一行 a:b=1:3，第二行 c:d=3:1 ----------------------------
# fig4_ab   <- p_total_nolab   + p_type_nolab    + plot_layout(widths = c(1, 2))
fig4_ab   <- p_total_nolab   + p_heat    + plot_layout(widths = c(1, 1.7))
fig4_ab
fig4_cd   <- p_scatter_nolab + p_mane_nolab    + plot_layout(widths = c(2, 1))

fig4_final <- fig4_ab /
  fig4_cd +
  plot_layout(heights  = c(1.5, 1)) +
  plot_annotation(tag_levels = "a") &          # 自动 a,b,c,d
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(4, 4, 4, 4))

## ── 3. 预览 & 保存 -----------------------------------------------------------
fig4_final    # 查看

## 5. 保存 ---------------------------------------------------------------------
dir.create("figures/fig4",showWarnings=FALSE,recursive=TRUE)
ggsave("figures/fig4/fig4abcd_final.png",fig4_final,width=10,height=8,dpi=300)
ggsave("figures/fig4/fig4abcd_final.pdf",fig4_final,width=10,height=8,dpi=300)

# remove all
rm(list =ls())
gc()
