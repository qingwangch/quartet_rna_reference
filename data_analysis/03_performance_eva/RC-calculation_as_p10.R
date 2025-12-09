#' RC-calculation_as_p10
#' 2024-07-04
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
library(data.table)

# Define a function to calculate RC for each batch (with consistent sample pair key)
calculate_batch_rc <- function(ref_data, test_data) {
  # 构建 isoform_compare 键用于 merge
  test_data$ase_compare <- paste(test_data$ase, test_data$compare)
  ref_data$ase_compare <- ref_data$isoform_compare
  
  # 合并
  merged_data <- merge(
    ref_data[, c("ase_compare", "mean_delta_psi_mean", "compare")],
    test_data[, c("ase_compare", "delta_psi", "batch")],
    by = "ase_compare"
  )
  
  # 只保留同时不为 NA 的行
  merged_data <- merged_data %>%
    filter(!is.na(mean_delta_psi_mean), !is.na(delta_psi))
  
  # 按 batch_id 和 compare 分组，计算 Pearson 相关系数（RC）和点数
  rc_results <- merged_data %>%
    group_by(batch, compare) %>%
    summarise(
      RC = cor(mean_delta_psi_mean, delta_psi, use = "complete.obs"),
      n_points = n(),
      .groups = "drop"
    )
  
  return(rc_results)
}

# Example usage
ref_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv"
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/salmon_seg_DAS_s264.txt"

# Load data
ref_data <- fread(ref_file) %>% as.data.frame()
test_data <- fread(test_file) %>% as.data.frame()

# Calculate RC for each batch
rc_results <- calculate_batch_rc(ref_data, test_data)

# View results
print(rc_results)

# SR
# batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
#                            "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
#                            "R_ILM_L6_B1", "P_BGI_L6_B1")
## Salmon
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/salmon_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_salmon <- rc_results
rc_salmon$pipeline <- "Salmon-SUPPA2"

## RSEM
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/rsem_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_rsem <- rc_results
rc_rsem$pipeline <- "RSEM-SUPPA2"

## kallisto
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/kallisto_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_kallisto <- rc_results
rc_kallisto$pipeline <- "kallisto-SUPPA2"

## StringTie2-SO
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/stringtie2_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_stringtie2_so <- rc_results
rc_stringtie2_so$pipeline <- "StringTie2-SO-SUPPA2"

## Merge results
rc_results_p4_so <- rbind(rc_salmon, rc_rsem, rc_kallisto, rc_stringtie2_so)
rc_results_p4_so$group <- "SO"

# LR
## IsoQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/isoquant/isoquant_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
batch_for_construction <- unique(test_data$batch)
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_isoquant <- rc_results
rc_isoquant$pipeline <- "IsoQuant-SUPPA2"

## Bambu
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/bambu/bambu_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_bambu <- rc_results
rc_bambu$pipeline <- "Bambu-SUPPA2"

## Oarfish
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/oarfish/oarfish_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_oarfish <- rc_results
rc_oarfish$pipeline <- "Oarfish-SUPPA2"

## StringTie2-LO
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/stringtie2_lo/stlo_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
test_data <- test_data[test_data$batch %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_stringtie2_lo <- rc_results
rc_stringtie2_lo$pipeline <- "StringTie2-LO-SUPPA2"

## Merge results
rc_results_p4_lo <- rbind(rc_isoquant, rc_bambu, rc_oarfish, rc_stringtie2_lo)
rc_results_p4_lo$group <- "LO"

# LS
## MPAQT
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/mpaqt/mpaqt_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_mpaqt <- rc_results
rc_mpaqt$pipeline <- "MPAQT-SUPPA2"
rc_mpaqt$group <- "LS"

## miniQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/miniquant/miniquant_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_miniquant <- rc_results
rc_miniquant$pipeline <- "miniQuant-SUPPA2"
rc_miniquant$group <- "LS"

## merge
rc_results_p8 <- rbind(rc_results_p4_so,rc_results_p4_lo)
rc_results_p10 <- rbind(rc_results_p8,rc_mpaqt,rc_miniquant)

# Save the results to a CSV file
write.csv(rc_results_p10, file = "/vast/projects/quartet_rna_refdata/analysis/performance/performance_assesement_AS_RC_results_LS_p10.csv", row.names = FALSE)

# summary
rc_summary <- rc_results_p10 %>%
  group_by(group) %>%
  summarise(
    mean_RC = mean(RC),
    sd_RC = sd(RC),
    .groups = "drop"
  )

print(rc_summary)

# remove all
rm(list=ls())
gc()
