#' heatmap_of_m6A_sites
#' Qingwang Chen
#' 2025-07-11

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# 指定新的库路径
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# 修改 .libPaths，使其优先使用自定义路径
.libPaths(c(custom_lib_path, .libPaths()))

print(.libPaths())


library(dplyr)
library(tidyr)
library(pheatmap)
library(purrr)
library(edgeR)
library(stringr)

coverage_threshold <- 100  # 最低 Nvalid_cov
mod_threshold      <- 90  # percent_modified 过滤暂不使用

read_m6A_bed <- function(bed_file) {
  cols <- c(
    "chrom","start","end","mod","score","strand",
    "start2","end2","color","Nvalid_cov","percent_modified",
    "Nmod","Ncanonical","Nother_mod","Ndelete",
    "Nfail","Ndiff","Nnocall"
  )
  df <- read.table(bed_file, header = FALSE, sep = "\t", quote = "")
  colnames(df) <- cols
  
  df %>%
    filter(
      Nvalid_cov >= coverage_threshold,
      percent_modified >= mod_threshold,
      Nmod        >= 10                 # 过滤掉 Nmod 小于 10 的位点
    ) %>%
    transmute(
      site_id = paste(chrom, start, end, strand, sep = "_"),
      Nmod    = as.integer(Nmod)
    )
}

bed_files_1 <- c(
  "/vast/scratch/users/xu.ya/quartet_rna_project/huge_data_result/D_ONT_LG_B1_F7_1/pileup.bed", "/vast/scratch/users/xu.ya/quartet_rna_project/huge_data_result/D_ONT_LG_B1_F7_2/pileup.bed", "/vast/scratch/users/xu.ya/quartet_rna_project/huge_data_result/D_ONT_LG_B1_F7_3/pileup.bed",
  "/vast/scratch/users/xu.ya/huge_data/m6a/D_ONT_LG_B1_M8_1/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/m6a/D_ONT_LG_B1_M8_2/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/m6a/D_ONT_LG_B1_M8_3/pileup.bed",
  "/vast/scratch/users/xu.ya/quartet_rna_project/ONT-DR-GRM-R01/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/m6a/D_ONT_LG_B1_D5_2/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/m6a/D_ONT_LG_B1_D5_3/pileup.bed",
  "/vast/scratch/users/xu.ya/quartet_rna_project/pod5_pass/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/m6a/D_ONT_LG_B1_D6_2/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/m6a/D_ONT_LG_B1_D6_3/pileup.bed"
)

bed_files_2 <- c(
  "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_F7_1/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_F7_2/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_F7_3/pileup.bed",
  "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_M8_1/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_M8_2/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_M8_3/pileup.bed",
  "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_D5_1/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_D5_2/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_D5_3/pileup.bed",
  "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_D6_1/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_D6_2/pileup.bed", "/vast/scratch/users/xu.ya/huge_data/modkit/D_ONT_LW_B1_D6_3/pileup.bed"
)

sample_names <- c(
  paste0("D_ONT_LG_B1_", c("F7_1","F7_2","F7_3","M8_1","M8_2","M8_3",
                      "D5_1","D5_2","D5_3","D6_1","D6_2","D6_3")),
  paste0("D_ONT_LW_B1_", c("F7_1","F7_2","F7_3","M8_1","M8_2","M8_3",
                      "D5_1","D5_2","D5_3","D6_1","D6_2","D6_3"))
)

bed_files <- setNames(c(bed_files_1, bed_files_2), sample_names)
bed_list <- map(bed_files, read_m6A_bed)
saveRDS(bed_list,"/vast/projects/quartet_rna_refdata/analysis/R/supp/m6a/bed_list_Nmod_s24.rds")

# bed_list <- readRDS("/vast/projects/quartet_rna_refdata/analysis/R/supp/m6a/bed_list.rds")
bed_list <- readRDS("/vast/projects/quartet_rna_refdata/analysis/R/supp/m6a/bed_list_Nmod50_s24.rds")

# 从 bed_list 构建 long_df
long_df <- map2_df(bed_list, names(bed_list), ~ {
  .x %>% mutate(sample = .y)
})

# pivot 成宽矩阵 rows = site_id, cols = sample, values = Nmod
wide_df <- long_df %>%
  pivot_wider(names_from  = sample,
              values_from = Nmod,
              values_fill = 0) %>% 
  as.data.frame()

# 计算 log2CPM
# lib_sizes <- colSums(mat0)
# cpm_mat    <- sweep(mat0, 2, lib_sizes/1e6, FUN = "/")
mat0 <- wide_df
rownames(mat0) <- mat0$site_id
mat0$site_id <- NULL
head(mat0)

cpm_mat <- cpm(mat0)

# FILTER SITES: 只保留 ≥100% 样本 CPM ≥ 10 的位点
cpm_df <- as.data.frame(cpm_mat)
# keep   <- rowSums(cpm_df >= 10) >= 1 * ncol(cpm_df)
keep   <- rowSums(mat0 >= 50) >= 1 * ncol(mat0)
length(keep)
# mat0     <- mat0[keep, , drop=FALSE]
cpm_mat  <- cpm_mat[keep, , drop=FALSE]
log2cpm  <- log2(cpm_mat + 0.01)

meta_df <- data.frame(
  sample_id = colnames(log2cpm),
  batch_id  = rep(c("Batch1","Batch2"), each = 12)
)

# 定义 calculate_ratios() 
calculate_ratios <- function(log_expr_mat, metadata, batch_col) {
  ubatch   <- unique(metadata[[batch_col]])
  ratio_D6 <- matrix(NA,
                     nrow = nrow(log_expr_mat),
                     ncol = ncol(log_expr_mat),
                     dimnames = dimnames(log_expr_mat))
  for (b in ubatch) {
    samps        <- metadata$sample_id[metadata[[batch_col]] == b]
    mat_batch    <- log_expr_mat[, samps, drop=FALSE]
    ref_mean     <- rowMeans(mat_batch[, grep("D6", samps), drop=FALSE])
    ratio_D6[, samps] <- sweep(mat_batch, 1, ref_mean, FUN = "-")
  }
  ratio_D6[apply(ratio_D6, 1, var, na.rm=TRUE) != 0, , drop=FALSE]
}
ratio_mat <- calculate_ratios(log2cpm, meta_df, "batch_id")

samps <- colnames(log2cpm)
batch <- sapply(strsplit(samps, "_"), function(x) paste(x[1:4], collapse = "_"))
group <- sapply(strsplit(samps, "_"), function(x) paste(x[5], collapse = "_"))

col_ann <- data.frame(
  Batch = factor(batch, levels=c("D_ONT_LG_B1","D_ONT_LW_B1")),
  Group = factor(group, levels=c("F7","M8","D5","D6"))
)
rownames(col_ann) <- samps
col_colors <- list(
  Batch = c(D_ONT_LG_B1="#1f78b4", D_ONT_LW_B1="#33a02c"),
  Group = c(F7="#FFC65D", M8="#F16745", D5="#4CC3D9", D6="#7BC8A4")
)
#  heatmap
cor_mat <- cor(ratio_mat, 
               use    = "pairwise.complete.obs", 
               method = "pearson")

pheatmap(
  cor_mat,
  main              = "Sample–Sample Correlation of log2CPM–D6",
  color             = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
  annotation_row    = col_ann,
  annotation_col    = col_ann,
  annotation_colors = col_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  border_color      = NA
)

