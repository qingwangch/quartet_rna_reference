# as-reference-description-validation.R
# 2025-06-13
# Qingwang Chen

# 工作目录 & R 库路径
setwd("/vast/projects/quartet_rna_refdata/analysis/")
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"
.libPaths(c(custom_lib_path, .libPaths()))

# 加载必要的包
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)   # 拼图
library(scales)      # comma()
library(ComplexUpset)
library(eulerr)
library(cowplot)

## ---------- 0. 读数据 ------------------------------------------------------
as_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv"
as_dt   <- fread(as_file)

## ---------- 1. 事件码 → 可读类别 ------------------------------------------
evt_lut <- c(
  SE = "Exon skipping",
  AF = "Alternative first exon",
  AL = "Alternative last exon",
  RI = "Intron retention",
  A3 = "Alternative 3′ end",
  A5 = "Alternative 5′ end",
  MX = "Mutually exclusive exons"
)

as_dt <- as_dt %>%
  mutate(evt_code = str_extract(ASE, ";([A-Z0-9]+):") |> str_remove_all("[:;]"),
         evt_name = recode(evt_code, !!!evt_lut, .default = NA_character_)) %>%
  filter(!is.na(evt_name))

evt_order <- as_dt %>% 
  as_tibble() %>%  # 确保是 tibble 类型
  dplyr::count(evt_name, name = "tot") %>% 
  arrange(tot) %>% 
  pull(evt_name)                       # e.g.  "Exon skipping" → …

panelA_df <- as_dt %>%
  distinct(compare, isoform_compare, evt_name) %>%
  dplyr::count(compare, evt_name, name = "n") %>%
  mutate(evt_name = factor(evt_name, levels = evt_order))   # ❷ 重新设定顺序


## ---------- panel a：sample × event ---------------------------------------
pA <- ggplot(panelA_df,
             aes(x = n, y = evt_name)) +
  geom_col(fill = "#BAD7FF", width = .6) +
  facet_grid(. ~ compare, scales = "free_x", switch = "y") +
  scale_x_continuous(labels = comma) +
  labs(x = "Number of alternative splicing events", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    # —— 网格 / 轴线 —— #
    panel.grid.major.x = element_line(colour = "grey85", linewidth = .3),
    panel.grid.major.y = element_blank(),
    axis.line.x        = element_line(colour = "black", linewidth = .5),
    axis.line.y        = element_line(colour = "black", linewidth = .5),
    axis.ticks.x       = element_line(colour = "black", linewidth = .5),
    axis.ticks.y       = element_line(colour = "black", linewidth = .5),
    axis.text.x        = element_text(margin = margin(t = 1.5)),
    axis.text.y        = element_text(margin = margin(r = 2)),
    
    # —— Facet strip —— #
    strip.placement = "outside",
    strip.text.x    = element_text(face = "bold", size = 11, margin = margin(b = 4)),
    strip.text.y    = element_text(face = "bold"),
    
    # —— 面板/图边缘 —— #
    panel.spacing.x = unit(0.7, "lines"),
    plot.margin     = margin(5.5, 5.5, 5.5, 6.5)
  )
pA

# -------------------------------------------------------------------
#  PANEL b  –  distribution of AS-event counts per gene – stacked % –
# -------------------------------------------------------------------
gene_bins <- as_dt %>% 
  mutate(gene_id = sub(";.*$", "", ASE)) %>% 
  distinct(compare, gene_id, ASE) %>% 
  dplyr::count(compare, gene_id, name = "as_n") %>% 
  mutate(bin = cut(as_n, breaks = c(0,1,2,3,4,5, Inf),
                   labels = c("1","2","3","4","5","≥6"), right = TRUE))

comp_dist <- gene_bins %>% 
  dplyr::count(compare, bin, name = "gene_n") %>% 
  group_by(compare) %>% 
  mutate(pct = gene_n / sum(gene_n)) %>% 
  ungroup()

# totals for “n = …”
tot_lbl <- gene_bins %>% 
  dplyr::count(compare, name = "total") %>% 
  mutate(label = paste0("n = ", format(total, big.mark = ",")))

bin_lvls <- c("1","2","3","4","5","≥6")
bin_cols <- setNames(colorRampPalette(c("#E1EEFF", "#3367C6"))(length(bin_lvls)),
                     bin_lvls)

pB <- ggplot(comp_dist, aes(compare, pct, fill = bin)) +
  geom_col(width = .75, colour = "black", linewidth = .25) +
  geom_text(data = tot_lbl, inherit.aes = FALSE,
            aes(compare, 1, label = label),
            vjust = -0.4, size = 3.1, fontface = "bold") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, .07))) +
  scale_fill_manual(values = bin_cols, name = "# AS events\nper gene") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Genes (fraction within compare)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85"),
    axis.line   = element_line(colour = "black", linewidth = .55),
    axis.ticks  = element_line(colour = "black", linewidth = .55),
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    plot.margin     = margin(5.5, 8, 5.5, 5.5)
  )
pB

# -------------------------------------------------------------------
#  COMBINE  –  final figure (a | b)  –––––––––––––––––––––––––––––––––
# -------------------------------------------------------------------
final_5ab <- (pA | pB) + 
  plot_layout(widths = c(2, 1), guides = "collect") &
  theme(axis.title.x = element_blank())

# preview
print(final_5ab)

ggsave("figures/fig5/ab_refAS_descript_bar.png",
       final_5ab, width = 12, height = 5, dpi = 300)
ggsave("figures/fig5/ab_refAS_descript_bar.pdf",
       final_5ab, width = 12, height = 5, dpi = 300)

gene_as_counts <- as_dt %>% 
  mutate(
    gene_id = sub(";.*$", "", ASE)            # keep ENSG… (strip everything after 1st ‘;’)
  ) %>% 
  distinct(compare, gene_id, ASE, .keep_all = FALSE) %>%  # one row = one event
  dplyr::count(compare, gene_id, name = "as_n")      # how many events for that gene in this compare

summary_tbl <- gene_as_counts %>% 
  group_by(compare) %>% 
  summarise(
    genes_total  = n(),                 # how many genes in that compare
    events_total = sum(as_n),           # total AS events for the compare
    mean_events  = mean(as_n),
    median_events = median(as_n),
    max_events   = max(as_n)
  )

print(summary_tbl, n = Inf)


# -------------------------------------------------------------------
#  PANEL c  – validation of AS events – scatter plots % –
# -------------------------------------------------------------------
qpcr_ase <- fread("isoforms/qpcr/deltapsi_20250513_U.csv") 
ref_ase <- as_dt
## ───────────────────────────────────────────────────────────────
## 2. 合并 & 整理  ------------------------------------------------
##  如果 ref_ase 中没有 compare，可根据 isoform_compare 衍生
##  例如 isoform_compare 末尾自带 “D5/D6”；如已有可跳过
if (!"compare" %in% names(ref_ase)) {
  ref_ase[, compare := sub(".*\\s", "", isoform_compare)]   # 仅示例
}

##  inner merge
setkey(ref_ase,   isoform_compare)
setkey(qpcr_ase,  ASE_compare)

ase_dt <- merge(
  ref_ase[, .(isoform_compare, compare, fc_ref = mean_delta_psi_mean)],
  qpcr_ase[,.(isoform_compare = ASE_compare, fc_val = delta_PSI)]
)

##  保持比较顺序
ase_dt[, compare := factor(compare,
                           levels = c("D5/D6", "F7/D6", "M8/D6"))]

## ───────────────────────────────────────────────────────────────
## 3. 统计表  -----------------------------------------------------
stat_dt <- ase_dt[, .(
  N    = .N,
  r    = cor(fc_ref, fc_val, use = "complete.obs"),
  RMSE = sqrt(mean((fc_val - fc_ref)^2, na.rm = TRUE))
), by = compare]

## ───────────────────────────────────────────────────────────────
## 4. 颜色 & 坐标范围  -------------------------------------------
nbt_cols <- c("D5/D6" = "#1b9e77",
              "F7/D6" = "#d95f02",
              "M8/D6" = "#7570b3")

rng <- range(unlist(ase_dt[, .(fc_ref, fc_val)]), na.rm = TRUE)
lim <- range(rng) + c(-1, 1) * diff(rng) * 0.05   # 5 % padding

## ───────────────────────────────────────────────────────────────
## 5. 绘图  -------------------------------------------------------
p_scatter_ase <- ggplot(ase_dt, aes(fc_ref, fc_val, colour = compare)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = .4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey80", linewidth = .3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey80", linewidth = .3) +
  geom_point(alpha = .35, size = 1.6, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "black", linewidth = .45) +
  facet_wrap(~compare, nrow = 1) +
  geom_text(data = stat_dt,
            aes(x = -Inf, y =  Inf,
                label = sprintf("r = %.2f\nRMSE = %.2f\nN = %d",
                                r, RMSE, N)),
            hjust = -0.05, vjust = 1.10,
            size = 3.4, fontface = "bold") +
  scale_colour_manual(values = nbt_cols) +
  labs(
    x = expression(Delta*"PSI in reference (LO)"),
    y = expression(Delta*"PSI by qPCR")
  ) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = .6),
    axis.ticks.length = unit(0.18, "cm"),
    axis.ticks        = element_line(colour = "black", linewidth = .5),
    strip.text        = element_text(size = 12, face = "bold",
                                     margin = margin(b = 4)),
    strip.background  = element_blank(),
    legend.position   = "none"
  )
p_scatter_ase
panel_C <- p_scatter_ase + plot_layout(guides="collect") &
  theme(plot.margin=margin(4,4,4,4))
panel_C

# -------------------------------------------------------------------
#  PANEL d  – correlation of AS events and alternative isoforms – scatter plots % –
# -------------------------------------------------------------------

############ refFC
## AS
############ refFC
refFC_as_final_p4_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv")
table(refFC_as_final_p4_ls$compare)

############ refDAS
refDAS_final_p4_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/RefData_AS_all_DAS_classified_u_0520.csv")
refDAS_final_p4_ls <- refDAS_final_p4_ls %>% filter(feature_type != "non-DAS")
table(refDAS_final_p4_ls$compare)

## isoform
############ refFC
refFC_final_p4_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv")
table(refFC_final_p4_ls$compare)
############ refDEI
refDEI_final_p4_ls <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv")
refDEI_final_p4_ls <- refDEI_final_p4_ls %>% filter(refDEI_final_p4_ls$Final != "non-DEI")
table(refDEI_final_p4_ls$compare)

### AS 
transcript_ids <- refDAS_final_p4_ls$alternative_transcripts
# transcript_ids <- refFC_as_final_p4_ls$alternative_transcripts
# 将每行文本按逗号分割，并将结果展平为一个向量
transcript_list <- unlist(strsplit(as.character(transcript_ids), split = ","))

# 统计唯一基因和转录本的数量
# 统计唯一转录本的数量
num_transcripts <- length(unique(transcript_list))
print(num_transcripts)

# setdiff(transcript_list,refDEI_final_p4_ls$isoform)
setdiff(transcript_list,refFC_final_p4_ls$isoform)

# common_isoforms <- intersect(transcript_list,refDEI_final_p4_ls$isoform)
common_isoforms <- intersect(transcript_list,refFC_final_p4_ls$isoform)

# 构建组合键 (isoform + compare) for DEIs
# deis <- refDEI_final_p4_ls %>%
#   mutate(
#     isoform = trimws(isoform),
#     compare = trimws(compare),
#     key = paste(isoform, compare, sep = "_")
#   ) %>%
#   pull(key)
deis <- refFC_final_p4_ls %>%
  mutate(
    isoform = trimws(isoform),
    compare = trimws(compare),
    key = paste(isoform, compare, sep = "_")
  ) %>%
  pull(key)
# 构建组合键 (isoform + compare) for DASEs
dases <- refDAS_final_p4_ls %>%
  select(compare, alternative_transcripts, mean_delta_psi_mean) %>%
  rowwise() %>%
  mutate(
    isoform = strsplit(as.character(alternative_transcripts), ",")
  ) %>%
  unnest(isoform) %>%
  mutate(
    isoform = trimws(isoform),
    compare = trimws(compare),
    key = paste(isoform, compare, sep = "_")
  ) %>%
  pull(key)

# 构建集合计数
venn_counts <- euler(c(
  `DEIs` = length(setdiff(deis, dases)),
  `DASEs` = length(setdiff(dases, deis)),
  `DEIs&DASEs` = length(intersect(deis, dases))
))

# 绘图
venn_plot <- plot(
  venn_counts,
  fills = list(fill = c("#999999", "#1F78B4"), alpha = 0.7),
  labels = list(font = 2, cex = 1.2),
  edges = TRUE,
  quantities = TRUE,
  legend = FALSE
)
venn_plot
# ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig5/Fig5h_proportional_venn_DEIs_DASEs.png",
#        plot = venn_plot, width = 4.5, height = 4.5, dpi = 300)

# 1. 提取 refDEI 中 log2FC
# refDEI_common <- refDEI_final_p4_ls %>%
#   # filter(Final != "non-DEI") %>%
#   mutate(
#     isoform = trimws(isoform),
#     compare = trimws(compare),
#     key = paste(isoform, compare, sep = "_")
#   ) %>%
#   select(key, isoform, compare, log2FC) %>%
#   distinct(key, .keep_all = TRUE)
refDEI_common <- refFC_final_p4_ls %>%
  mutate(
    isoform = trimws(isoform),
    compare = trimws(compare),
    key = paste(isoform, compare, sep = "_")
  ) %>%
  select(key, isoform, compare, log2FC) %>%
  distinct(key, .keep_all = TRUE)
# 2. 提取 refDAS 中 delta_psi
refDAS_common <- refDAS_final_p4_ls %>%
  select(compare, alternative_transcripts, mean_delta_psi_mean) %>%
  rowwise() %>%
  mutate(isoform = strsplit(as.character(alternative_transcripts), ",")) %>%
  unnest(isoform) %>%
  mutate(
    isoform = trimws(isoform),
    compare = trimws(compare),
    key = paste(isoform, compare, sep = "_")
  ) %>%
  select(key, isoform, compare, delta_psi = mean_delta_psi_mean) %>%
  distinct(key, .keep_all = TRUE)

# 3. 合并
merged_common <- inner_join(refDAS_common, refDEI_common, by = "key")
setdiff(merged_common$key, intersect(unique(deis), unique(dases)))
# fwrite(merged_common, "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/data/fig4f_merged_isoform_AS_ref.csv")


# 4. 计算每个 compare 分组的相关系数
cor_stats <- merged_common %>%
  group_by(compare.x) %>%
  summarise(r = round(cor(delta_psi, log2FC, use = "complete.obs"), 2)) %>%
  rename(compare = compare.x)



# 5. 绘图
nbt_colors <- c("D5/D6" = "#1b9e77", "F7/D6" = "#d95f02", "M8/D6" = "#7570b3")

# 创建散点图
p1 <- ggplot(merged_common, aes(x = delta_psi, y = log2FC, color = compare.x)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  scale_color_manual(values = nbt_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = expression(Delta*"PSI (DASEs)"),
    y = "log2FC (Ref DEIs)",
    color = "Comparison"
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  # 添加 r 值注释
  annotate("text",
           x = min(merged_common$delta_psi, na.rm = TRUE),
           y = max(merged_common$log2FC, na.rm = TRUE),
           label = paste0("r = ", round(cor(merged_common$delta_psi, merged_common$log2FC, use = "complete.obs"), 2),"\n N = ", nrow(merged_common)),
           hjust = 0, vjust = 2.2,
           fontface = "bold", size = 5)
print(p1)

# 拼图
final_combined <- plot_grid(pA,pB,p_scatter_ase, p1, labels = c("e", "f","g", "h"), ncol = 2, rel_widths = c(2, 1),align = "hv")

final_combined
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/Fig5efgh_correlation_combined.pdf",
       plot = final_combined, width = 10, height = 7, dpi = 300)

# remove all
rm(list=ls())
gc()
