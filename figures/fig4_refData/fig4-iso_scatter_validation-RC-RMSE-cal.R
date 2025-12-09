# iso_scatter_validation-RC-RMSE-cal.R
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

## ── 2. 读取参考集 ────────────────────────────────────────────────────────────
# Isoform (refFC / refDEI)
ref_fc  <- fread("isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv") |> mutate(source = "LO")
ref_dei <- fread("isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv") |>
  filter(Final != "non-DEI") |> mutate(source = "LO")

qpcr_iso <- fread("isoforms/qpcr/log2fc_20250409_U.csv")   |> 
  mutate(source = "qPCR")                        # ← 添加 source
qpcr_ase <- fread("isoforms/qpcr/deltapsi_20250513_U.csv") |> 
  mutate(source = "qPCR")

## ── Isoform scatter plot : RefFC (LO)  vs  qPCR ───────────────

# 1. 准备数据 ------------------------------------------------------------------
iso_scatter_df <- inner_join(
  ref_fc  %>%                       # come from your earlier code
    select(isoform_compare, compare, isoform, log2FC_ref = log2FC),
  qpcr_iso %>%                      # qPCR validation file
    select(isoform_compare, compare, log2FC_val = log2FC),
  by = c("isoform_compare", "compare")
)
fwrite(iso_scatter_df,"figures/fig4/data/isoform_scatter_refFC_vs_qPCR.csv")

# annotation
annotation <- fread("/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf")
## ---- 1. 仅保留 transcript 这一 feature 行 -------------------------------
anno_tx <- annotation[V3 == "transcript"]

## ---- 2. 定义要抽取的属性关键字 ------------------------------------------
attrs <- c("gene_id", "transcript_id",
           "gene_type", "gene_name",
           "transcript_type", "transcript_name")

## ---- 3. 写一个辅助函数：从 V9 中抓取属性值 -------------------------------
get_attr <- function(x, key) {
  m <- str_match(x, paste0(key, ' "([^"]+)"'))
  m[, 2]               # 若匹配不到返回 NA
}

## ---- 4. 逐列解析并生成新的 data.table ------------------------------------
res <- anno_tx[, lapply(attrs, \(k) get_attr(V9, k))]
setnames(res, attrs)

## ---- 5. 去重（同一转录本只保留一行） ------------------------------------
res <- unique(res, by = "transcript_id")
length(unique(res$gene_id))
fwrite(res, file = "figures/fig4/data/gene_tx_mapping.tsv", sep = "\t")

# merge annotaion
iso_scatter_ann <- iso_scatter_df %>%
  left_join(res, by = c("isoform" = "transcript_id")) %>%
  relocate(isoform, .before = everything()) 

fwrite(iso_scatter_ann, "figures/fig4/data/iso_scatter_annotated.csv")

# （可选）把三组 sample pair 设定 factor 顺序，便于 Facet 排序
iso_scatter_df$compare <- factor(
  iso_scatter_df$compare,
  levels = c("D5/D6", "F7/D6", "M8/D6")
)

# 2. 计算每个 sample pair 的 Pearson r，用来标到图里 ----------------------------
stat_tbl <- iso_scatter_df |>
  group_by(compare) |>
  summarise(
    N     = n(),
    r     = cor(log2FC_ref, log2FC_val, use = "complete.obs"),
    RMSE  = sqrt(mean((log2FC_val - log2FC_ref)^2, na.rm = TRUE)),
    .groups = "drop"
  )

# 3. 画图 ----------------------------------------------------------------------
p_iso_scatter <- ggplot(iso_scatter_df,
                        aes(x = log2FC_ref,
                            y = log2FC_val,
                            colour = compare)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey40") +
  geom_point(alpha = .35, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "black", linewidth = .4) +
  facet_wrap(~ compare, nrow = 1) +
  scale_colour_manual(values = c("D5/D6"="#1b9e77",
                                 "F7/D6"="#d95f02",
                                 "M8/D6"="#7570b3")) +
  labs(title    = "Isoform reference (LO) vs qPCR validation",
       x = "log2FC in reference (LO)",
       y = "log2FC by qPCR",
       colour = "Sample pair") +
  theme_classic(base_size = 13) +
  theme(strip.text      = element_text(face = "bold", size = 12),
        strip.background  = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(0.18, "cm")) +
  
  ## ←← 新的三行注释：r, RMSE 和 N
  geom_text(
    data = stat_tbl,
    aes(x = -Inf, y = Inf,
        label = sprintf(" r = %.2f\nRMSE = %.2f\n N = %d", r, RMSE, N)),
    hjust = -0.05, vjust = 1.1,
    size  = 3.3, fontface = "bold"
  )

# 5. 输出 ----------------------------------------------------------------------
print(p_iso_scatter)
p_iso_scatter+canvas(width = 8,height = 4)

ggsave(
  filename = "figures/fig4/isoform_scatter_refFC_vs_qPCR.png",
  plot     = p_iso_scatter,
  width    = 8,  # inch
  height   = 4,  # 一行 facets
  dpi      = 300
)
ggsave(
  filename = "figures/fig4/isoform_scatter_refFC_vs_qPCR.pdf",
  plot     = p_iso_scatter,
  width    = 8,
  height   = 4
)

# remove all
rm(list=ls())
gc()
