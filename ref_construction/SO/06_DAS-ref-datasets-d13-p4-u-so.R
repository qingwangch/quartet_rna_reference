#' DAS-ref-datasets-d13-p4-u-so
#' Multi-input Support
#' 2025-07-09

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/suppa2/")

# Specify a new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths so that the custom path has priority
.libPaths(c(custom_lib_path, .libPaths()))

# Check the current library paths
print(.libPaths())


# Import necessary libraries
library(tidyr)
library(data.table)

# Define the calculation function
## Relative Standard Error (RSE)
rse <- function(x) {
  if (length(x[!is.na(x)]) <= 1) {
    return(0)  # If there is only one non-NA data point, return 0 to avoid NA
  }
  sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) * sqrt(length(x[!is.na(x)])))
}


# Calculate fold-change metrics (ΔPSI)
calculate_mean_delta_psi_metrics <- function(var1, ref_mean_delta_psi_list) {
  mean_delta_psi_values <- c()
  
  for (ref_mean_delta_psi in ref_mean_delta_psi_list) {
    matching_row <- ref_mean_delta_psi[ref_mean_delta_psi$var1 == var1, ]
    if (nrow(matching_row) > 0) {
      mean_delta_psi_values <- c(mean_delta_psi_values, matching_row$mean_delta_psi)
    }
  }
  
  list(
    mean_delta_psi_mean = ifelse(length(mean_delta_psi_values) > 0, 
                                 mean(mean_delta_psi_values, na.rm = TRUE), NA),
    uncertainty = ifelse(length(mean_delta_psi_values) > 0, 
                         rse(mean_delta_psi_values) * 100, NA)
  )
}


# Automatically load reference datasets
ref_files <- list(
  kallisto = "./ref_data_construction/kallisto_s156_RefData_DAS.csv",
  stringtie = "./ref_data_construction/stso_s156_RefData_DAS.csv",
  salmon = "./ref_data_construction/salmon_s156_RefData_DAS.csv",
  rsem = "./ref_data_construction/rsem_s156_RefData_DAS.csv"
)

# Load all reference files and store them in `refDEI_list`
refDEI_list <- lapply(ref_files, function(file) {
  fread(file) %>%
    as.data.frame()
})

# Combine all reference datasets, keeping every isoform in the `final` column
ref_data_combined <- bind_rows(lapply(names(ref_files), function(source) {
  fread(ref_files[[source]]) %>%
    as_tibble() %>%  # convert data.table to tibble
    mutate(source = source) %>%
    dplyr::select(var1, final, source)  # explicitly use dplyr::select()
})) %>% unique()

# Frequency statistics: count how many times each Var1–Final pair appears across sources
var1_final_counts <- ref_data_combined %>%
  group_by(var1, final) %>%
  summarise(Freq = n(), .groups = 'drop') %>%
  separate(var1, into = c("isoform", "compare"), sep = " ")
var1_final_counts$isoform_compare <- paste(var1_final_counts$isoform,var1_final_counts$compare)

# Define a function to classify `feature_type`
classify_feature_type <- function(df) {
  # Group by isoform_compare
  classified_df <- df %>%
    group_by(isoform_compare) %>%
    summarise(
      Final_unique = n_distinct(final),  
      Freq_min = min(Freq),              
      Freq_max = max(Freq),              
      Final_consensus = ifelse(Final_unique == 1 & Freq_min >= 2, unique(final), NA), 
      .groups = "drop"
    ) %>%
    mutate(
      feature_type = case_when(
        !is.na(Final_consensus) ~ Final_consensus,  
        TRUE ~ "Conflict"                           
      )
    )
  
  # Merge the `feature_type` results back into the original data
  final_result <- df %>%
    left_join(classified_df %>% dplyr::select(isoform_compare, feature_type), by = "isoform_compare")
  
  return(final_result)
}

# Example usage
# Apply the classification function
final_classified_df <- classify_feature_type(var1_final_counts)
final_classified_df <- final_classified_df[final_classified_df$feature_type != "Conflict",]
final_classified_df <- final_classified_df[final_classified_df$final != "conflicting",]

# Compute mean ΔPSI and uncertainty
fc_metrics_list <- lapply(final_classified_df$isoform_compare, calculate_mean_delta_psi_metrics, refDEI_list)

# Merge the calculated metrics into final_classified_df
fc_metrics_df <- bind_rows(fc_metrics_list)
result_df <- bind_cols(final_classified_df, fc_metrics_df)

# Add annotations
# overall
# Read the local GTF-derived IOE file
gtf_file <- fread("/vast/projects/quartet_rna_refdata/reference/suppa2/gencode_v43.all_AS.ioe") %>% as.data.frame()
gtf_file <- gtf_file[,c("event_id","alternative_transcripts","total_transcripts")]
result_df <- merge(result_df,gtf_file,by.x = "isoform", by.y = "event_id")


# Compute uncertainty quantiles
uncertainty_quantiles <- quantile(result_df$uncertainty, c(0.25, 0.50, 0.75))

# Add a classification column: classify all isoforms by Final
result_df <- result_df %>%
  mutate(
    classification = case_when(
      Freq >= 4 ~ "High Confidence",
      Freq >= 3 ~ "Medium Confidence",
      TRUE ~ "Low Confidence"
    )
  )

# Inspect the classification results
print(table(result_df$classification))
table(result_df$compare)
colnames(result_df)[1] <- "ASE"
# Save the results
write.csv(result_df, "./ref_data_construction/so-e-g/RefData_AS_all_DAS_classified_u_0910.csv", row.names = FALSE)

# filter DAS
das_result_df <- result_df[result_df$final != "non-DAS", ]
# Inspect filtered classification results
print(table(das_result_df$classification))
print(table(das_result_df$compare))

# remove all 
rm(list=ls())
gc()
