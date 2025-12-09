#' RC-calculation_p10
#' 2024-07-03
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
  # 统一 comparison 格式（如 D5vsD6 -> D5/D6）
  test_data$compare <- gsub("vs", "/", test_data$comparison)
  
  # 构建 isoform_compare 键用于 merge
  test_data$isoform_compare <- paste(test_data$feature, test_data$compare)
  ref_data$isoform_compare <- paste(ref_data$isoform, ref_data$compare)
  
  # 合并
  merged_data <- merge(
    ref_data[, c("isoform_compare", "log2FC", "compare")],
    test_data[, c("isoform_compare", "logFC", "batch_id")],
    by = "isoform_compare"
  )
  
  # 按 batch_id 和 compare 分组，计算 Pearson 相关系数（RC）和点数
  rc_results <- merged_data %>%
    group_by(batch_id, compare) %>%
    summarise(
      RC = cor(log2FC, logFC, use = "complete.obs"),
      n_points = n(),
      .groups = "drop"
    )
  
  return(rc_results)
}

# Example usage
ref_file <- "isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv"
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_Salmon_isoform.txt"

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
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_Salmon_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_salmon <- rc_results
rc_salmon$pipeline <- "Salmon"

## RSEM
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_RSEM_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_rsem <- rc_results
rc_rsem$pipeline <- "RSEM"

## kallisto
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_kallisto_id_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_kallisto <- rc_results
rc_kallisto$pipeline <- "kallisto"

## StringTie2-SO
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_StringTie2_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_stringtie2_so <- rc_results
rc_stringtie2_so$pipeline <- "StringTie2-SO"

## Merge results
rc_results_p4_so <- rbind(rc_salmon, rc_rsem, rc_kallisto, rc_stringtie2_so)
rc_results_p4_so$group <- "SO"

# LR
## IsoQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_IsoQuant_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_isoquant <- rc_results
rc_isoquant$pipeline <- "IsoQuant"

## Bambu
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_Bambu_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_bambu <- rc_results
rc_bambu$pipeline <- "Bambu"

## Oarfish
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_Oarfish_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_oarfish <- rc_results
rc_oarfish$pipeline <- "Oarfish"

## StringTie2-LO
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_StringTie2_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_stringtie2_lo <- rc_results
rc_stringtie2_lo$pipeline <- "StringTie2-LO"

## Merge results
rc_results_p4_lo <- rbind(rc_isoquant, rc_bambu, rc_oarfish, rc_stringtie2_lo)
rc_results_p4_lo$group <- "LO"

# LS
## MPAQT
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform_tpm/DE_results_MPAQT_merged_all.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_mpaqt <- rc_results
rc_mpaqt$pipeline <- "MPAQT"
rc_mpaqt$group <- "LS"

## miniquant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_miniQuant_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rc_results <- calculate_batch_rc(ref_data, test_data)
rc_miniquant <- rc_results
rc_miniquant$pipeline <- "miniQuant"
rc_miniquant$group <- "LS"

## merge
rc_results_p8 <- rbind(rc_results_p4_so,rc_results_p4_lo)
rc_results_p10 <- rbind(rc_results_p8,rc_mpaqt,rc_miniquant)

# Save the results to a CSV file
write.csv(rc_results_p10, file = "/vast/projects/quartet_rna_refdata/analysis/performance/performance_assesement_RC_results_LO_p10.csv", row.names = FALSE)

# summary
rc_summary <- rc_results_p10 %>%
  group_by(pipeline) %>%
  summarise(
    mean_RC = mean(RC),
    sd_RC = sd(RC),
    .groups = "drop"
  )

print(rc_summary)
 
# remove all
rm(list=ls())
gc()
