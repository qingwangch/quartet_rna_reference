#' isoform-ref-datasets-d13-p4-so
#' Multi-input Support
#' 2025-07-09

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# Specify a new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths so that the custom path has priority
.libPaths(c(custom_lib_path, .libPaths()))

# Check the current library paths
print(.libPaths())

# Import necessary libraries
library(edgeR)
library(magrittr)
library(purrr)
library(dplyr)
library(data.table)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(tidyr)

# Ratio-based reference datasets
refFC_final_kallisto <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/kallisto_ref_expr.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_stringite <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/stringtie2_ref_expr.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_salmon <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/salmon_ref_expr.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_rsem <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/rsem_ref_expr.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)

# Extract the unique values in column Var1 from each data frame
var1_kallisto <- unique(refFC_final_kallisto$Var1)
var1_stringite <- unique(refFC_final_stringite$Var1)
var1_salmon <- unique(refFC_final_salmon$Var1)
var1_rsem <- unique(refFC_final_rsem$Var1)

# Combine all unique values
all_var1 <- c(var1_kallisto, var1_stringite, var1_salmon, var1_rsem)

# Count how many times each Var1 appears
var1_counts <- table(all_var1)

# Keep Var1 values that appear at least twice
var1_at_least_two <- var1_counts[var1_counts >= 2]

# Create a new data frame containing these Var1 values and their counts
result_df <- data.frame(Var1 = names(var1_at_least_two), Freq = as.integer(var1_at_least_two))
# Split Var1 into isoform and compare columns
result_df <- separate(result_df, col = Var1, into = c("isoform", "compare"), sep = " ")
result_df$isoform_compare <- paste(result_df$isoform,result_df$compare)
table(result_df$compare)

## uchar
rse <- function(x) {
  sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) * sqrt(length(x[!is.na(x)])))
}

# Load all reference files
refFC_list <- list(refFC_final_kallisto, refFC_final_stringite, refFC_final_salmon, refFC_final_rsem)
# Combine refFC_list into one data frame, adding a source tag (optional)
combined_refFC <- bind_rows(
  lapply(seq_along(refFC_list), function(i) {
    refFC <- refFC_list[[i]]
    refFC$source <- paste0("ref_", i)  # 添加来源信息（可选）
    refFC
  })
)

# Convert to data.table
setDT(combined_refFC)
setDT(result_df)

# Set a key for combined_refFC
setkey(combined_refFC, Var1)

# Pre-compute fc_mean, log2FC, and uncertainty
metrics <- combined_refFC[, .(
  fc_mean = mean(fc, na.rm = TRUE),
  log2FC = mean(log2(fc), na.rm = TRUE),
  uncertainty = rse(fc) * 100
), by = Var1]

# Merge the metrics into result_df
result_df <- merge(result_df, metrics, by.x = "isoform_compare", by.y = "Var1", all.x = TRUE)

# Inspect the results
print(head(result_df))
table(result_df$compare)

# Add annotations
# overall
# Read local GTF-derived transcript info
transcript_info<- fread("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/transcript_info_reference.csv") %>% as.data.frame()
result_df <- merge(result_df,transcript_info,by.x = "isoform", by.y = "transcript_id")

# Calculate the 25th, 50th, and 75th percentiles of uncertainty
uncertainty_q25 <- quantile(result_df$uncertainty, 0.25)
uncertainty_q50 <- quantile(result_df$uncertainty, 0.50)
uncertainty_q75 <- quantile(result_df$uncertainty, 0.75)

# Print the quantiles
print(c(Q25 = uncertainty_q25, Q50 = uncertainty_q50, Q75 = uncertainty_q75))

# Define the classification function
classify_isoforms <- function(Freq, uncertainty, q50, q75) {
  if (Freq >= 3 && uncertainty <= q50) {
    return("High Confidence")
  } else if (Freq >= 2 && uncertainty > q50 && uncertainty <= q75) {
    return("Medium Confidence")
  } else {
    return("Low Confidence")
  }
}

# Add a classification column to result_df
result_df <- result_df %>%
  mutate(
    classification = mapply(
      classify_isoforms, 
      Freq = Freq, 
      uncertainty = uncertainty, 
      q50 = uncertainty_q50, 
      q75 = uncertainty_q75
    )
  )

# Check the classification results
table(result_df$classification)

# Save the final data frame
refFC_final_p4 <- result_df
write.csv(refFC_final_p4, "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/ref_expr_b13_p4_s156_u_20250910.csv", row.names = FALSE)

# summary the reference datasets
# Define the classification criteria
refFC_final_p4 <- refFC_final_p4 %>%
  mutate(
    Condition = case_when(
      classification == "High Confidence" ~ paste("Freq >= 3 & uncertainty <=", round(uncertainty_q50, 2)),
      classification == "Medium Confidence" ~ paste("Freq >= 2 &", round(uncertainty_q50, 2), "< uncertainty <=", round(uncertainty_q75, 2)),
      classification == "Low Confidence" ~ "Others"
    )
  )

# Count the number of isoforms in each class
classification_counts <- refFC_final_p4[ , .N, by = .(classification, Condition)]
classification_counts[ , Percentage := N / sum(N) * 100]

# Count isoforms with uncertainty < 30
num_isoforms_below_30 <- sum(refFC_final_p4$uncertainty < 30, na.rm = TRUE)

# Total number of isoforms
total_isoforms <- nrow(refFC_final_p4)

# Proportion of isoforms with uncertainty < 30
percentage_below_30 <- (num_isoforms_below_30 / total_isoforms) * 100

# Output the result
cat(sprintf("uncertainty 小于 30%% 的 isoform 占比: %.2f%%\n", percentage_below_30))

# remove all
rm(list = ls())
gc()
