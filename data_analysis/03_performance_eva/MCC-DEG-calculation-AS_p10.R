#' MCC-DEG-calculation-AS_p10
#' 2025-07-04
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

# 计算 Matthews Correlation Coefficient
calculate_mcc <- function(tp, tn, fp, fn) {
  # 分子
  num <- as.numeric(tp * tn - fp * fn)
  # 分母
  den <- sqrt(
    as.numeric(tp + fp) * 
      as.numeric(tp + fn) * 
      as.numeric(tn + fp) * 
      as.numeric(tn + fn)
  )
  # 防止除零
  if (den == 0) {
    return(NA_real_)
  }
  num / den
}

# Define a function to calculate MCC
calculate_batch_performance <- function(ref_deis, test_deis) {
  # 生成 isoform_compare 列，格式一致
  ref_deis <- ref_deis %>%
    mutate(ase_compare = paste(ASE, compare))  # e.g., ENST000004 / D5/D6
  
  test_deis <- test_deis %>%
    mutate(feature_type = case_when(
      p_value < 0.05 & delta_psi >= 0.05 ~ "up-regulate",
      p_value < 0.05 & delta_psi <= -0.05 ~ "down-regulate",
      TRUE ~ "non-DAS"
    )) %>%
    mutate(
      ase_compare = paste(ase, compare)
    )
  
  # 提取并合并感兴趣的列
  test_deis_f <- test_deis %>% select(ase_compare, feature_type, batch, compare)
  ref_deis_f <- ref_deis %>% select(ase_compare, feature_type_ref)
  
  # 合并
  merged_data <- merge(ref_deis_f, test_deis_f, by = "ase_compare") %>%
    mutate(
      ref_label = case_when(
        feature_type_ref == "up-regulate" ~ "Upregulate",
        feature_type_ref == "down-regulate" ~ "Downregulate",
        TRUE ~ "NonDAS"
      ),
      test_label = case_when(
        feature_type == "up-regulate" ~ "Upregulate",
        feature_type == "down-regulate" ~ "Downregulate",
        TRUE ~ "NonDAS"
      )
    )
  
  # 汇总指标，按 batch + compare 分组
  results <- merged_data %>%
    group_by(batch, compare) %>%
    summarise(
      tp = sum((ref_label == "Upregulate" & test_label == "Upregulate") |
                 (ref_label == "Downregulate" & test_label == "Downregulate"), na.rm = TRUE),
      tn = sum(ref_label == "NonDAS" & test_label == "NonDAS", na.rm = TRUE),
      fp = sum((ref_label == "NonDAS" & test_label %in% c("Upregulate", "Downregulate")) |
                 (ref_label == "Upregulate" & test_label == "Downregulate") |
                 (ref_label == "Downregulate" & test_label == "Upregulate"), na.rm = TRUE),
      fn = sum((ref_label %in% c("Upregulate", "Downregulate") & test_label == "NonDAS") |
                 (ref_label == "Upregulate" & test_label == "Downregulate") |
                 (ref_label == "Downregulate" & test_label == "Upregulate"), na.rm = TRUE),
      MCC = calculate_mcc(tp, tn, fp, fn),
      Precision = ifelse((tp + fp) == 0, NA, tp / (tp + fp)),
      Recall = ifelse((tp + fn) == 0, NA, tp / (tp + fn)),
      FNR = ifelse((tp + fn) == 0, NA, fn / (tp + fn)),
      F1 = ifelse((2 * tp + fp + fn) == 0, NA, 2 * tp / (2 * tp + fp + fn)),
      n_points = n(),
      .groups = "drop"
    )
  
  return(results)
}


# Example usage
# ref_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/RefData_DAS_LS_t2_u_20250408.csv"
ref_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/RefData_AS_all_DAS_classified_u_0520.csv"

test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/salmon_seg_DAS_s264.txt"

# Load data
ref_dases <- fread(ref_file) %>% as.data.frame()
# Add missing columns to ref_dases
ref_dases <- ref_dases %>%
  mutate(
    feature_type_ref = feature_type,
  )

test_dases <- fread(test_file) %>% as.data.frame()

# Calculate MCC for each batch
dase_results <- calculate_batch_performance(ref_dases, test_dases)

# View results
print(dase_results)

# SO
# batch_for_construction <- c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
#                             "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
#                             "R_ILM_L6_B1", "P_BGI_L6_B1")
## Salmon
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/salmon_seg_DAS_s264.txt"
test_dases <- fread(test_file) %>% as.data.frame()
# test_dases <- test_dases[test_dases$batch %in% batch_for_construction,]
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_salmon <- dase_results
dase_salmon$pipeline <- "Salmon-SUPPA2"

## RSEM
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/rsem_seg_DAS_s264.txt"
test_dases <- fread(test_file) %>% as.data.frame()
# test_dases <- test_dases[test_dases$batch %in% batch_for_construction,]
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_rsem <- dase_results
dase_rsem$pipeline <- "RSEM-SUPPA2"

## kallisto
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/kallisto_seg_DAS_s264.txt"
test_dases <- fread(test_file) %>% as.data.frame()
# test_dases <- test_dases[test_dases$batch %in% batch_for_construction,]
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_kallisto <- dase_results
dase_kallisto$pipeline <- "kallisto-SUPPA2"

## StringTie2-SO
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/stringtie2_seg_DAS_s264.txt"
test_dases <- fread(test_file) %>% as.data.frame()
# test_dases <- test_dases[test_dases$batch %in% batch_for_construction,]
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_stringtie2_so <- dase_results
dase_stringtie2_so$pipeline <- "StringTie2-SO-SUPPA2"

## Merge results
dase_results_p4 <- rbind(dase_salmon, dase_rsem, dase_kallisto, dase_stringtie2_so)
dase_results_p4$group <- "SO"

mcc_summary <- dase_results_p4 %>%
  group_by(pipeline) %>%
  summarise(
    mean_MCC = mean(MCC),
    sd_MCC = sd(MCC),
    .groups = "drop"
  )

print(mcc_summary)

# LO
## IsoQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/isoquant/isoquant_leg_DAS_s144.txt"
test_dases <- fread(test_file) %>% as.data.frame()
batch_for_construction <- unique(test_data$batch)
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_isoquant <- dase_results
dase_isoquant$pipeline <- "IsoQuant-SUPPA2"

## Bambu
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/bambu/bambu_leg_DAS_s144.txt"
test_dases <- fread(test_file) %>% as.data.frame()
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_bambu <- dase_results
dase_bambu$pipeline <- "Bambu-SUPPA2"

## Oarfish
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/oarfish/oarfish_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
test_data <- test_data[test_data$batch %in% batch_for_construction,]
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_oarfish <- dase_results
dase_oarfish$pipeline <- "Oarfish-SUPPA2"

## StringTie2
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/stringtie2_lo/stlo_leg_DAS_s144.txt"
test_data <- fread(test_file) %>% as.data.frame()
test_data <- test_data[test_data$batch %in% batch_for_construction,]
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_stringtie2_lo <- dase_results
dase_stringtie2_lo$pipeline <- "StringTie2-LO-SUPPA2"

## Merge results
dase_results_p4_lo <- rbind(dase_isoquant, dase_bambu, dase_oarfish, dase_stringtie2_lo)
dase_results_p4_lo$group <- "LO"

## MPAQT
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/mpaqt/mpaqt_leg_DAS_s144.txt"
test_dases <- fread(test_file) %>% as.data.frame()
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_mpaqt <- dase_results
dase_mpaqt$pipeline <- "MPAQT-SUPPA2"
dase_mpaqt$group <- "LS"

## miniQuant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/suppa2/miniquant/miniquant_leg_DAS_s144.txt"
test_dases <- fread(test_file) %>% as.data.frame()
dase_results <- calculate_batch_performance(ref_dases, test_dases)
dase_miniquant <- dase_results
dase_miniquant$pipeline <- "miniQuant-SUPPA2"
dase_miniquant$group <- "LS"

## Merge results
dase_results_p8 <- rbind(dase_results_p4, dase_results_p4_lo)
dase_results_p10 <- rbind(dase_results_p8, dase_mpaqt, dase_miniquant)

# Save the results to a CSV file
write.csv(dase_results_p10, file = "/vast/projects/quartet_rna_refdata/analysis/performance/performance_assesement_DAS_MCC_results_LS_p10.csv", row.names = FALSE)

mcc_summary <- dase_results_p10 %>%
  group_by(group) %>%
  summarise(
    mean_MCC = mean(MCC, na.rm = TRUE),
    sd_MCC = sd(MCC, na.rm = TRUE),
    .groups = "drop"
  )

print(mcc_summary)

# Clean up the environment
rm(list = ls())
gc()
