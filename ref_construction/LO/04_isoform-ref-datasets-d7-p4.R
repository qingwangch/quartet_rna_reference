#' isoform-ref-datasets-d7-p4
#' Multi-input Support
#' 2024-05-16

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

# Ratio-based reference datasets
refFC_final_bambu <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/bambu_ref_expr_s84.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_stringite <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/stringtie2_ref_expr_s84.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_oarfish <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/oarfish_ref_expr_s84.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_isoquant <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/isoquant_ref_expr_s84.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)

# Extract the unique values from the Var1 column in each data frame
var1_bambu <- unique(refFC_final_bambu$Var1)
var1_stringite <- unique(refFC_final_stringite$Var1)
var1_oarfish <- unique(refFC_final_oarfish$Var1)
var1_isoquant <- unique(refFC_final_isoquant$Var1)

# Combine all unique values
all_var1 <- c(var1_bambu, var1_stringite, var1_oarfish, var1_isoquant)

# Count the number of occurrences of each Var1
var1_counts <- table(all_var1)

# Filter Var1 values that appear at least twice
var1_at_least_two <- var1_counts[var1_counts >= 2]

# Create a new data frame that contains these Var1 values and their occurrence counts
result_df <- data.frame(Var1 = names(var1_at_least_two), Freq = as.integer(var1_at_least_two))
# Use the separate() function to split the Var1 column into two columns: isoform and compare
result_df <- separate(result_df, col = Var1, into = c("isoform", "compare"), sep = " ")
table(result_df$compare)

# Define the calculation function
## uchar
rse <- function(x) {
  if (length(x[!is.na(x)]) <= 1) {
    return(0)  # If there is only one data point, return 0 directly to avoid NA
  }
  sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) * sqrt(length(x[!is.na(x)])))
}

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

# Load all reference files
refFC_list <- list(refFC_final_bambu, refFC_final_stringite, refFC_final_oarfish, refFC_final_isoquant)

# Add fc, log2FC, and uncertainty columns to result_df
result_df <- result_df %>%
  rowwise() %>%
  mutate(
    fc = calculate_fc_metrics(paste(isoform, compare), refFC_list)$fc_mean,
    log2FC = calculate_fc_metrics(paste(isoform, compare), refFC_list)$log2FC,
    medianp = calculate_fc_metrics(paste(isoform, compare), refFC_list)$medianp,
    uncertainty = calculate_fc_metrics(paste(isoform, compare), refFC_list)$uncertainty,
    isoform_compare = paste(isoform, compare)
  ) %>%
  ungroup()

# check
print(head(result_df))

# overall annotation
transcript_info<- fread("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/transcript_info_reference.csv") %>% as.data.frame()
result_df <- merge(result_df,transcript_info,by.x = "isoform", by.y = "transcript_id")

# Calculate the 25th and 50th percentiles of uncertainty
uncertainty_q25 <- quantile(result_df$uncertainty, 0.25)
uncertainty_q50 <- quantile(result_df$uncertainty, 0.50)
uncertainty_q75 <- quantile(result_df$uncertainty, 0.75)

# check
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

# check
table(result_df$classification)

# save
refFC_final_p4 <- result_df
write.csv(refFC_final_p4, "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv", row.names = FALSE)

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

# Calculate the count for each category
classification_counts <- refFC_final_p4 %>%
  count(classification, Condition) %>%
  ungroup() %>%
  mutate(Percentage = n / sum(n) * 100)

# Create a descriptive table
description_table <- classification_counts %>%
  select(Classification = classification, Condition, Count = n, Percentage)

# check
print(description_table)

# Count the number of isoforms with uncertainty less than 30
num_isoforms_below_30 <- sum(refFC_final_p4$uncertainty < 30, na.rm = TRUE)

# Calculate the total number of isoforms
total_isoforms <- nrow(refFC_final_p4)

# Compute the proportion
percentage_below_30 <- (num_isoforms_below_30 / total_isoforms) * 100

# output
cat(sprintf("uncertainty 小于 30%% 的 isoform 占比: %.2f%%\n", percentage_below_30))

# remove all
rm(list = ls())
gc()
