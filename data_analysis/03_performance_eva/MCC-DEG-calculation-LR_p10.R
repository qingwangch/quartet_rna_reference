#' MCC-DEG-calculation-LR_p8
#' 2025-06-23
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
    mutate(isoform_compare = paste(isoform, compare))  # e.g., ENST000004 / D5/D6
  
  test_deis <- test_deis %>%
    mutate(feature_type = case_when(
      PValue < 0.05 & logFC >= 1 ~ "up-regulate",
      PValue < 0.05 & logFC <= -1 ~ "down-regulate",
      TRUE ~ "non-DEG"
    )) %>%
    mutate(
      compare = gsub("vs", "/", comparison),  # D5vsD6 -> D5/D6
      isoform_compare = paste(feature, compare)
    )
  
  # 提取并合并感兴趣的列
  test_deis_f <- test_deis %>% select(isoform_compare, feature_type, batch_id, compare)
  ref_deis_f <- ref_deis %>% select(isoform_compare, Final)
  
  # 合并
  merged_data <- merge(ref_deis_f, test_deis_f, by = "isoform_compare") %>%
    mutate(
      ref_label = case_when(
        Final == "up-regulate" ~ "Upregulate",
        Final == "down-regulate" ~ "Downregulate",
        TRUE ~ "NonDEG"
      ),
      test_label = case_when(
        feature_type == "up-regulate" ~ "Upregulate",
        feature_type == "down-regulate" ~ "Downregulate",
        TRUE ~ "NonDEG"
      )
    )
  
  # 汇总指标，按 batch_id + compare 分组
  results <- merged_data %>%
    group_by(batch_id, compare) %>%
    summarise(
      tp = sum((ref_label == "Upregulate" & test_label == "Upregulate") |
                 (ref_label == "Downregulate" & test_label == "Downregulate"), na.rm = TRUE),
      tn = sum(ref_label == "NonDEG" & test_label == "NonDEG", na.rm = TRUE),
      fp = sum((ref_label == "NonDEG" & test_label %in% c("Upregulate", "Downregulate")) |
                 (ref_label == "Upregulate" & test_label == "Downregulate") |
                 (ref_label == "Downregulate" & test_label == "Upregulate"), na.rm = TRUE),
      fn = sum((ref_label %in% c("Upregulate", "Downregulate") & test_label == "NonDEG") |
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
ref_file <- "isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv"

test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_Salmon_isoform.txt"

# Load data
ref_deis <- fread(ref_file) %>% as.data.frame()
# Add missing columns to ref_deis
ref_deis <- ref_deis %>%
  mutate(
    feature_type = Final,
  )

test_deis <- fread(test_file) %>% as.data.frame()

# Calculate MCC for each batch
dei_results <- calculate_batch_performance(ref_deis, test_deis)

# View results
print(dei_results)

# SO
# batch_for_construction <- c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
#                             "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
#                             "R_ILM_L6_B1", "P_BGI_L6_B1")
## Salmon
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_Salmon_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
# test_deis <- test_deis[test_deis$batch_id %in% batch_for_construction,]
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_salmon <- dei_results
dei_salmon$pipeline <- "Salmon"

## RSEM
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_RSEM_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
# test_deis <- test_deis[test_deis$batch_id %in% batch_for_construction,]
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_rsem <- dei_results
dei_rsem$pipeline <- "RSEM"

## kallisto
test_file <- "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_kallisto_id_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
# test_deis <- test_deis[test_deis$batch_id %in% batch_for_construction,]
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_kallisto <- dei_results
dei_kallisto$pipeline <- "kallisto"

## StringTie2-SO
test_file <-"/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_StringTie2_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
# test_deis <- test_deis[test_deis$batch_id %in% batch_for_construction,]
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_stringtie2_so <- dei_results
dei_stringtie2_so$pipeline <- "StringTie2-SO"

## Merge results
dei_results_p4 <- rbind(dei_salmon, dei_rsem, dei_kallisto, dei_stringtie2_so)
dei_results_p4$group <- "SO"

# LO
## IsoQuant
test_file <-  "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_IsoQuant_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_isoquant <- dei_results
dei_isoquant$pipeline <- "IsoQuant"

## Bambu
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_Bambu_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_bambu <- dei_results
dei_bambu$pipeline <- "Bambu"

## Oarfish
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_Oarfish_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
test_deis <- test_deis %>% 
  mutate(feature = sub("^((ENST[^|]+)).*", "\\1", feature))

dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_oarfish <- dei_results
dei_oarfish$pipeline <- "Oarfish"

## StringTie2
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_StringTie2_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_stringtie2_lo <- dei_results
dei_stringtie2_lo$pipeline <- "StringTie2-LO"

## Merge results
dei_results_p4_lo <- rbind(dei_isoquant, dei_bambu, dei_oarfish, dei_stringtie2_lo)
dei_results_p4_lo$group <- "LO"

## MPAQT
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform_tpm/DE_results_MPAQT_merged_all.txt"
test_deis <- fread(test_file) %>% as.data.frame()
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_mpaqt <- dei_results
dei_mpaqt$pipeline <- "MPAQT"
dei_mpaqt$group <- "LS"

## miniquant
test_file <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/DE_results_miniQuant_isoform.txt"
test_deis <- fread(test_file) %>% as.data.frame()
dei_results <- calculate_batch_performance(ref_deis, test_deis)
dei_miniquant <- dei_results
dei_miniquant$pipeline <- "miniQuant"
dei_miniquant$group <- "LS"

## Merge results
dei_results_p8 <- rbind(dei_results_p4, dei_results_p4_lo)
dei_results_p10 <- rbind(dei_results_p8, dei_mpaqt, dei_miniquant)

# Save the results to a CSV file
write.csv(dei_results_p10, file = "/vast/projects/quartet_rna_refdata/analysis/performance/performance_assesement_DEG_MCC_results_LO_p10.csv", row.names = FALSE)

mcc_summary <- dei_results_p10 %>%
  group_by(group) %>%
  summarise(
    mean_MCC = mean(MCC),
    sd_MCC = sd(MCC),
    .groups = "drop"
  )

print(mcc_summary)

# Clean up the environment
rm(list = ls())
gc()
