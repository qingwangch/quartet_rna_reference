#' DEI-ref-datasets-d7-p4-lo
#' Multi-input Support
#' 2025-05-20

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# Specify a new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths so that the custom path has priority
.libPaths(c(custom_lib_path, .libPaths()))

# Check the current library paths
print(.libPaths())

# library import
library(magrittr)
library(purrr)
library(dplyr)
library(tidyr)
library(data.table)

# Define the calculation function
## Relative Standard Error (RSE)
rse <- function(x) {
  sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) * sqrt(length(x[!is.na(x)])))
}

# Calculate the metric for fold change (FC)
calculate_fc_metrics <- function(var1, refFC_list) {
  fc_values <- c()
  p_vals <- c()
  
  for (refFC in refFC_list) {
    matching_row <- refFC[refFC$Var1 == var1, ]
    if (nrow(matching_row) > 0) {
      fc_values <- c(fc_values, matching_row$fc)
      p_vals   <- c(p_vals, matching_row$medianp)
    }
  }
  
  list(
    fc_mean = ifelse(length(fc_values) > 0, mean(fc_values, na.rm = TRUE), NA),
    log2FC = ifelse(length(fc_values) > 0, mean(log2(fc_values), na.rm = TRUE), NA),
    medianp  = ifelse(length(p_vals) > 0, median(p_vals,  na.rm = TRUE) , NA),
    uncertainty = ifelse(length(fc_values) > 0, rse(fc_values) * 100, NA)
  )
}

# Automatically load the reference dataset
ref_files <- list(
  bambu = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/bambu_RefData_DEIs_s84.csv",
  stringtie = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/stringtie2_RefData_DEIs_s84.csv",
  oarfish = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/oarfish_RefData_DEIs_s84.csv",
  isoquant = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/isoquant_RefData_DEIs_s84.csv"
)

# Load all reference files and store them in refDEI_list
refDEI_list <- lapply(ref_files, function(file) {
  fread(file) %>%
    as.data.frame()
})

# Merge all reference datasets and retain the isoforms in the Final column
ref_data_combined <- bind_rows(lapply(names(ref_files), function(source) {
  fread(ref_files[[source]]) %>%
    as_tibble() %>%  # 将 data.table 转换为 tibble
    mutate(source = source) %>%
    dplyr::select(Var1, Final, source)  # 明确指定 dplyr::select()
})) %>% unique()

# Frequency statistics: count how many times each Var1–Final combination appears across different sources
var1_final_counts <- ref_data_combined %>%
  group_by(Var1, Final) %>%
  summarise(Freq = n(), .groups = 'drop') %>%
  separate(Var1, into = c("isoform", "compare"), sep = " ")
var1_final_counts$isoform_compare <- paste(var1_final_counts$isoform,var1_final_counts$compare)

# Define a function to perform feature_type classification
classify_feature_type <- function(df) {
  # Group by isoform_compare
  classified_df <- df %>%
    group_by(isoform_compare) %>%
    summarise(
      Final_unique = n_distinct(Final),  
      Freq_min = min(Freq),              
      Freq_max = max(Freq),              
      Final_consensus = ifelse(Final_unique == 1 & Freq_min >= 2, unique(Final), NA), 
      .groups = "drop"
    ) %>%
    mutate(
      feature_type = case_when(
        !is.na(Final_consensus) ~ Final_consensus,  # If Final is consistent and the support frequency is ≥ 2
        TRUE ~ "Conflict"                           
      )
    )
  
  # Merge the feature_type results back into the original data
  final_result <- df %>%
    left_join(classified_df %>% dplyr::select(isoform_compare, feature_type), by = "isoform_compare")
  
  return(final_result)
}

# Example application
# Call the classification function
final_classified_df <- classify_feature_type(var1_final_counts)
final_classified_df <- final_classified_df[final_classified_df$feature_type != "Conflict",]
final_classified_df <- final_classified_df[final_classified_df$Final != "conflicting",]

# Calculate fc, log2FC, and uncertainty
fc_metrics_list <- lapply(final_classified_df$isoform_compare, calculate_fc_metrics, refFC_list = refDEI_list)

# Merge the calculation results into the var1_final_counts data frame
fc_metrics_df <- bind_rows(fc_metrics_list)
result_df <- bind_cols(final_classified_df, fc_metrics_df)
head(result_df)
# Calculate the quantiles of uncertainty
uncertainty_quantiles <- quantile(result_df$uncertainty, c(0.25, 0.50, 0.75))

# Add a classification column: classify all isoforms with Final
result_df <- result_df %>%
  mutate(
    classification = case_when(
      Freq >= 3 & uncertainty <= uncertainty_quantiles[2] ~ "High Confidence",
      Freq >= 2 & uncertainty > uncertainty_quantiles[2] & uncertainty <= uncertainty_quantiles[3] ~ "Medium Confidence",
      TRUE ~ "Low Confidence"
    )
  )

# print
print(table(result_df$classification))

# save
write.csv(result_df, "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv", row.names = FALSE)

# remove all 
rm(list=ls())
gc()
