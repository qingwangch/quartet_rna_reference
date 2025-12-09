#' RMSE-calculation-as
#' Qingwang Chen
#' 2025-04-23
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
  # 1) 构建合并用的 key
  test_data$ase_compare <- paste(test_data$ase, test_data$compare)
  ref_data$ase_compare  <- ref_data$isoform_compare
  
  # 2) 合并
  merged_data <- merge(
    ref_data[, c("ase_compare", "mean_delta_psi_mean", "compare")],
    test_data[, c("ase_compare", "delta_psi", "batch")],
    by = "ase_compare",
    all = FALSE
  )
  
  # 3) 去掉任何含 NA 的行
  merged_data <- merged_data %>%
    filter(!is.na(mean_delta_psi_mean), !is.na(delta_psi))
  
  # 4) 按 batch 和 compare 分组，计算 RMSE 和有效点数
  rmse_results <- merged_data %>%
    group_by(batch, compare) %>%
    summarise(
      RMSE     = sqrt(mean((mean_delta_psi_mean - delta_psi)^2)),
      n_points = n(),
      .groups  = "drop"
    )
  
  return(rmse_results)
}


# Ref
ref_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv"

# Load reference data
ref_data <- fread(ref_file) %>% as.data.frame()

# Calculate RMSE for each pipeline
### SO
# batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
#                            "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
#                            "R_ILM_L6_B1", "P_BGI_L6_B1")

## Salmon
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/salmon_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_salmon <- rmse_results
rmse_salmon$pipeline <- "Salmon-SUPPA2"

## RSEM
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/rsem_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_rsem <- rmse_results
rmse_rsem$pipeline <- "RSEM-SUPPA2"

## kallisto
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/kallisto_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_kallisto <- rmse_results
rmse_kallisto$pipeline <- "kallisto-SUPPA2"

## StringTie2
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/stringtie2_seg_DAS_s264.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_stringtie2 <- rmse_results
rmse_stringtie2$pipeline <- "StringTie2-SO-SUPPA2"

## Merge results
rmse_results_p4 <- rbind(rmse_salmon, rmse_rsem, rmse_kallisto, rmse_stringtie2)
rmse_results_p4$group <- "SO"

# Save the results to a CSV file
# write.csv(rmse_results_p4, file = "./results/tables/performance_assesement_RMSE_results_SR_p4.csv", row.names = FALSE)

# LR
## IsoQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/isoquant/isoquant_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
batch_for_construction <- unique(test_data$batch)
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_isoquant <- rmse_results
rmse_isoquant$pipeline <- "IsoQuant-SUPPA2"

## Bambu
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/bambu/bambu_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_bambu <- rmse_results
rmse_bambu$pipeline <- "Bambu-SUPPA2"

## Oarfish
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/oarfish/oarfish_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_oarfish<- rmse_results
rmse_oarfish$pipeline <- "Oarfish-SUPPA2"

## StringTie2
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/stringtie2_lo/stlo_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
# test_data <- test_data[test_data$batch %in% batch_for_construction,]
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_stringtie2_lo <- rmse_results
rmse_stringtie2_lo$pipeline <- "StringTie2-LO-SUPPA2"

## Merge results
rmse_results_p4_lo <- rbind(rmse_isoquant, rmse_bambu, rmse_oarfish, rmse_stringtie2_lo)
rmse_results_p4_lo$group <- "LO"

# LS
## MPAQT
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/mpaqt/mpaqt_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_mpaqt <- rmse_results
rmse_mpaqt$pipeline <- "MPAQT-SUPPA2"
rmse_mpaqt$group <- "LS"

## miniQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/miniquant/miniquant_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
rmse_results <- calculate_batch_rmse(ref_data, test_data)
rmse_miniquant <- rmse_results
rmse_miniquant$pipeline <- "miniQuant-SUPPA2"
rmse_miniquant$group <- "LS"

## Merge all results
rmse_results_p8 <- rbind(rmse_results_p4,rmse_results_p4_lo)
rmse_results_p10 <- rbind(rmse_results_p8,rmse_mpaqt,rmse_miniquant)

# Save the results to a CSV file
write.csv(rmse_results_p10, file = "/vast/projects/quartet_rna_refdata/analysis/performance/performance_assesement_AS_RMSE_results_LS_p10.csv", row.names = FALSE)

# summary
rmse_summary <- rmse_results_p10 %>%
  group_by(group) %>%
  summarise(
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    .groups = "drop"
  )

print(rmse_summary)

# Clean up the environment
rm(list = ls())
gc()
