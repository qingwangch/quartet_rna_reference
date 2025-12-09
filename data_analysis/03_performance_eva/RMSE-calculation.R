#' RMSE-calculation
#' Qingwang Chen
#' 2025-06-23
#' 

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
library(Metrics) # For the rmse function

# Define a function to calculate RMSE for each batch and comparison
calculate_batch_rmse <- function(ref_data, test_data) {
  # 统一 comparison 格式（如 D5vsD6 -> D5/D6）
  test_data$compare <- gsub("vs", "/", test_data$comparison)
  
  # 构建 isoform_compare 键
  test_data$isoform_compare <- paste(test_data$feature, test_data$compare)
  ref_data$isoform_compare <- paste(ref_data$isoform, ref_data$compare)
  
  # 合并
  merged_data <- merge(
    ref_data[, c("isoform_compare", "log2FC", "compare")], 
    test_data[, c("isoform_compare", "logFC", "batch_id")],
    by = "isoform_compare"
  )
  
  # Group by batch_id 和 compare，计算 RMSE 和点数
  rmse_results <- merged_data %>%
    group_by(batch_id, compare) %>%
    summarise(
      RMSE = sqrt(mean((log2FC - logFC)^2, na.rm = TRUE)),
      n_points = n(),
      .groups = "drop"
    )
  
  return(rmse_results)
}

# Example usage
ref_file <- "isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv"

# Load reference data
ref_data <- fread(ref_file) %>% as.data.frame()

# Calculate RMSE for each pipeline
### SO
# batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
#                            "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
#                            "R_ILM_L6_B1", "P_BGI_L6_B1")

## Salmon
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_Salmon_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_salmon <- rmse_results
rmse_salmon$pipeline <- "Salmon"

## RSEM
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_RSEM_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_rsem <- rmse_results
rmse_rsem$pipeline <- "RSEM"

## kallisto
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_kallisto_id_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_kallisto <- rmse_results
rmse_kallisto$pipeline <- "kallisto"

## StringTie2
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_StringTie2_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch_id %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_stringtie2 <- rmse_results
rmse_stringtie2$pipeline <- "StringTie2-SO"

## Merge results
rmse_results_p4 <- rbind(rmse_salmon, rmse_rsem, rmse_kallisto, rmse_stringtie2)
rmse_results_p4$group <- "SO"

# Save the results to a CSV file
# write.csv(rmse_results_p4, file = "./results/tables/performance_assesement_RMSE_results_SR_p4.csv", row.names = FALSE)

# LR
## IsoQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_IsoQuant_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_isoquant <- rmse_results
rmse_isoquant$pipeline <- "IsoQuant"

## Bambu
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_Bambu_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_bambu <- rmse_results
rmse_bambu$pipeline <- "Bambu"

## Oarfish
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_Oarfish_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_oarfish<- rmse_results
rmse_oarfish$pipeline <- "Oarfish"

## StringTie2
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_StringTie2_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_stringtie2_lo <- rmse_results
rmse_stringtie2_lo$pipeline <- "StringTie2-LO"

## Merge results
rmse_results_p4_lo <- rbind(rmse_isoquant, rmse_bambu, rmse_oarfish, rmse_stringtie2_lo)
rmse_results_p4_lo$group <- "LO"

# LS
## MPAQT
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform_tpm/DE_results_MPAQT_merged_all.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_mpaqt <- rmse_results
rmse_mpaqt$pipeline <- "MPAQT"
rmse_mpaqt$group <- "LS"

## miniquant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_miniQuant_isoform.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_miniquant <- rmse_results
rmse_miniquant$pipeline <- "miniQuant"
rmse_miniquant$group <- "LS"

## Merge all results
rmse_results_p8 <- rbind(rmse_results_p4,rmse_results_p4_lo)
rmse_results_p10 <- rbind(rmse_results_p8,rmse_mpaqt,rmse_miniquant)

# Save the results to a CSV file
write.csv(rmse_results_p10, file = "/vast/projects/quartet_rna_refdata/analysis/performance/performance_assesement_RMSE_results_LO_p10.csv", row.names = FALSE)

# summary
rmse_summary <- rmse_results_p10 %>%
  group_by(pipeline) %>%
  summarise(
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    .groups = "drop"
  )

print(rmse_summary)

# Clean up the environment
rm(list = ls())
gc()
