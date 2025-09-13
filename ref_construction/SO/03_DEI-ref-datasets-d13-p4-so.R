#' DEI-ref-datasets-d13-p4-so
#' Multi-input Support
#' 2025-07-09
#' 

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


# DEI reference datasets
refDEI_kallisto <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/kallisto_RefData_DEIs.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refDEI_stringite <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/stringtie2_RefData_DEIs.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refDEI_salmon <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/salmon_RefData_DEIs.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refDEI_rsem <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/rsem_RefData_DEIs.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)

# Add a source column to indicate the origin of each data frame
refDEI_kallisto$source <- "kallisto"
refDEI_stringite$source <- "stringite"
refDEI_salmon$source <- "salmon"
refDEI_rsem$source <- "rsem"

# Combine the four data frames
ref_data_combined <- rbind(refDEI_kallisto[, c("Var1", "Final", "source")],
                           refDEI_stringite[, c("Var1", "Final", "source")],
                           refDEI_salmon[, c("Var1", "Final", "source")],
                           refDEI_rsem[, c("Var1", "Final", "source")])

# Remove duplicate rows
ref_data_combined <- unique(ref_data_combined)

# Frequency statistics: count how many times each Var1–Final pair appears across sources
var1_final_counts <- ref_data_combined %>%
  group_by(Var1, Final) %>%
  summarise(Freq = n(), .groups = 'drop') %>%
  separate(Var1, into = c("isoform", "compare"), sep = " ")
var1_final_counts$isoform_compare <- paste(var1_final_counts$isoform,var1_final_counts$compare)

# Define a function to classify `feature_type`
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
        !is.na(Final_consensus) ~ Final_consensus,  
        TRUE ~ "conflicting"                           
      )
    )
  
  final_result <- df %>%
    left_join(classified_df %>% dplyr::select(isoform_compare, feature_type), by = "isoform_compare")
  
  return(final_result)
}

# Example usage
# Apply the classification function
final_classified_df <- classify_feature_type(var1_final_counts)
final_classified_df <- final_classified_df[final_classified_df$feature_type != "conflicting",]
final_classified_df <- final_classified_df[final_classified_df$Final != "conflicting",]

# Define the calculation function
## Relative Standard Error (RSE)
rse <- function(x) {
  sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) * sqrt(length(x[!is.na(x)])))
}

# Load all reference files
refDEI_list <- list(refDEI_kallisto, refDEI_stringite, refDEI_salmon, refDEI_rsem)
# Combine refDEI_list into one data frame, adding a source tag (optional)
combined_refDEI <- bind_rows(
  lapply(seq_along(refDEI_list), function(i) {
    refDEI <- refDEI_list[[i]]
    refDEI$source <- paste0("ref_", i)  # 添加来源信息（可选）
    refDEI
  })
)

# Convert to data.table
setDT(combined_refDEI)
setDT(final_classified_df)

# Set a key for combined_refDEI
setkey(combined_refDEI, Var1)

# Pre-compute fc_mean, log2FC, and uncertainty
metrics <- combined_refDEI[, .(
  fc_mean = mean(fc, na.rm = TRUE),
  log2FC = mean(log2(fc), na.rm = TRUE),
  uncertainty = rse(fc) * 100
), by = Var1]

# Merge the metrics into final_classified_df
final_classified_df <- merge(final_classified_df, metrics, by.x = "isoform_compare", by.y = "Var1", all.x = TRUE)

# Inspect the results
print(head(final_classified_df))
table(final_classified_df$compare)

# Compute uncertainty quantiles (25th, 50th, 75th)
uncertainty_quantiles <- quantile(final_classified_df$uncertainty, c(0.25, 0.50, 0.75))

# Add a classification column: classify all isoforms by Final
final_classified_df <- final_classified_df %>%
  mutate(
    classification = case_when(
      Freq >= 3 & uncertainty <= uncertainty_quantiles[2] ~ "High Confidence",
      Freq >= 2 & uncertainty > uncertainty_quantiles[2] & uncertainty <= uncertainty_quantiles[3] ~ "Medium Confidence",
      TRUE ~ "Low Confidence"
    )
  )

# Check the classification results
print(table(final_classified_df$classification))

# Save the final data frame
write.csv(final_classified_df, "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/RefData_DEIs_all_isoforms_classified_u_20250910.csv", row.names = FALSE)

# remove all 
rm(list=ls())
gc()
