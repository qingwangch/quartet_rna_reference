# novel_isoforms_classification.R
# 2025-06-15
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
library(scales)
library(stringr) 
library(ggview) 
# ---------------------------------------------------------------
#  extract_novel_isoforms_classification.R
#  ------------------------
#  Read a SQANTI3 "classification.txt", extract novel isoforms,
#  write them to a TSV, and print basic summary statistics.
#
# ---------------------------------------------------------------


# ── 1. 读取 SQANTI3 classification 表 ────────────────────────────────────────
sq3_cls <- fread(
  "/vast/projects/quartet_rna_refdata/analysis/sqanti3/novel/quartet_s84_c_gencode_v43_filter/quartet_s84_c_gencode_v43_filter.SQ3_classification.txt"
)
cat("Rows read in classification:", nrow(sq3_cls), "\n")

# ── 2. 读取 GTF（仅 transcript 行）并提取 transcript_id ──────────────────────
gtf_path <- "/vast/projects/quartet_rna_refdata/analysis/sqanti3/novel/quartet_s84_c_gencode_v43_filter/filter_ud/quartet_s84_c_gencode_v43_filter_filter_sqanti.cds.sorted.gene.gtf"

gtf_dt <- fread(
  gtf_path,
  sep = "\t",
  header = FALSE,
  skip = "#",           # 跳过 header 行
  col.names = c("chr","source","type","start","end","score",
                "strand","phase","attr")
)[type == "transcript"] |>
  # 提取 transcript_id
  mutate(
    transcript_id = str_match(attr, 'transcript_id "([^"]+)"')[,2],
    gene_id       = str_match(attr, 'gene_id "([^"]+)"')[,2]
  )

cat("Transcript rows in GTF:", nrow(gtf_dt), "\n")

# ── 3. 识别 novel transcript（不以 ENST 开头即可视为 novel） ────────────────
novel_tx <- gtf_dt |>
  filter(!str_detect(transcript_id, "^ENST")) |>
  pull(transcript_id)

cat("Novel transcript IDs:", length(novel_tx), "\n")

# ── 4. 对齐到 SQANTI3 classification 表 ───────────────────────────────────────
# SQANTI3 默认 isoform 列就是 transcript id
novel_cls <- sq3_cls[ isoform %in% novel_tx ]

cat("Matched rows in classification:", nrow(novel_cls), "\n")

# ── 5. 保存 novel_cls ────────────────────────────────────────────────────────
out_dir  <- "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/data/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "st3_quartet_novel_sqanti3_classification.tsv")
fwrite(novel_cls, out_file, sep = "\t")         # 如需 CSV 改成 `sep=","`

cat("✓ Saved", nrow(novel_cls), "novel transcripts to:", out_file, "\n")

outfile <- "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig4/novel_isoform_classification"            # no extension

# ---------------- 2.  CODING STATUS ---------------------------
## derive a 3-level factor
## 1. 数据整理 ---------------------------------------------------------------
sq3 <- novel_cls %>%                   # novel_cls 来自你前面代码
  mutate(coding_status = case_when(
    coding == "coding" & !predicted_NMD ~ "Protein-coding",
    coding == "coding" &  predicted_NMD ~ "NMD-sensitive",
    TRUE                                ~ "Non-coding"
  ),
  coding_status = factor(coding_status,
                         levels = c("Protein-coding",
                                    "NMD-sensitive",
                                    "Non-coding")))

class_order <- c("novel_in_catalog","novel_not_in_catalog",
                 "fusion", "incomplete-splice_match",
                 "antisense", "intergenic","genic")

plot_df <- sq3 %>%
  filter(structural_category %in% class_order) %>%
  mutate(structural_category = factor(structural_category,
                                      levels = class_order,
                                      labels = c("Novel in catalog","Novel not in catalog",
                                                 "Fusion", "Incomplete splice match",
                                                 "Antisense", "Intergenic","Genic"))) %>%
  dplyr::count(structural_category, coding_status, name = "n") %>%
  group_by(structural_category) %>%
  mutate(total        = sum(n),
         pc_pct       = n[coding_status == "Protein-coding"] / total,
         pc_pct_label = ifelse(is.na(pc_pct), "(0%)",
                               sprintf("(%s%%)", round(pc_pct*100)))) %>%
  ungroup()

## ── 顶端标签（两行） ─────────────────────────────────────────────
top_lbl <- plot_df %>%
  group_by(structural_category) %>%
  summarise(total = unique(total),
            pc_pct = unique(pc_pct[!is.na(pc_pct)]),
            .groups = "drop") %>%
  mutate(label = sprintf("%s\n(%s%%)",
                         scales::comma(total),
                         round(pc_pct*100)),
         y_pos = total)            # 放在柱顶

## ---- 4. plot ---------------------------------------------------------------
cols <- c(
  "Protein-coding" = "#FF7043",  # 珊瑚橘
  "NMD-sensitive"  = "#2BA3C6",  # 海蓝
  "Non-coding"     = "#756BB1"   # 淡靛紫
)

p <- ggplot(plot_df,
            aes(structural_category, n, fill = coding_status)) +
  geom_col(width = .75, colour = "black", linewidth = .25) +
  
  ## 顶端两行文字
  geom_text(data = top_lbl, inherit.aes = FALSE,
            aes(structural_category, y_pos, label = label),
            vjust = -0.35, size = 3.4, fontface = "bold") +
  
  scale_fill_manual(values = cols, name = "Coding status") +     # 换行让内部图例更紧凑
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, .12))) +
  
  labs(x = "Isoform classification\n(SQANTI)", y = NULL,
       title    = sprintf("Novel isoforms (total = %s)",
                          scales::comma(sum(plot_df$n))),
       subtitle = "percentage protein-coding") +
  
  coord_cartesian(clip = "off") +
  
  theme_minimal(base_size = 13) +
  theme(
    ## 网格去掉，轴线 / 须线保留
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = .3),
    axis.line.x  = element_line(colour = "black", linewidth = .6),
    axis.line.y  = element_line(colour = "black", linewidth = .6),
    axis.ticks.x = element_line(colour = "black", linewidth = .5),
    axis.ticks.y = element_line(colour = "black", linewidth = .5),
    
    ## x 轴不显示数字（分类标签也不要）——如果想留文字，把 element_blank 换回 element_text
    axis.text.x  = element_text(angle = 35, hjust = 1),
    
    ## 图例放图内右上角
    legend.position   = c(.82, .82),        # 相对坐标 (0-1, 0-1)
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", colour = NA, linewidth = .3),
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9),
    
    plot.title      = element_text(face = "bold", hjust = .5, size = 14),
    plot.subtitle   = element_text(hjust = .5, margin = margin(b = 8)),
    plot.margin     = margin(5.5, 14, 5.5, 5.5)
  )

print(p)
p+canvas(width = 8,height = 4)

# ---------------- 5.  OUTPUT ----------------------------------
ggsave(paste0(outfile, ".pdf"), plot = p, width = 8, height = 4)
ggsave(paste0(outfile, ".png"), plot = p, width = 8, height = 4, dpi = 300)

# also save the tidy count table
fwrite(plot_df, file = paste0(outfile, "_counts.tsv"), sep = "\t")

message("✔ Figure saved to ",
        paste0(outfile, ".pdf / .png"),
        "\n✔ Counts table  : ", outfile, "_counts.tsv")

# remove all
rm(list = ls())
gc()
