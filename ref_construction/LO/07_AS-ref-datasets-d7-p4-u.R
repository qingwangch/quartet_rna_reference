#' AS-ref-datasets-d7-p4-u
#' Multi-input Support
#' 2025-05-20

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/suppa2/")

# Specify a new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths so that the custom path has priority
.libPaths(c(custom_lib_path, .libPaths()))

# Check the current library paths
print(.libPaths())

# library import
library(data.table)
library(tidyr) 

# Ratio-based reference datasets
refFC_final_bambu <- read.csv("./ref_data_construction/bambu_s84_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_stringite <- read.csv("./ref_data_construction/stringtie2_s84_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_oarfish <- read.csv("./ref_data_construction/oarfish_s84_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_isoquant <- read.csv("./ref_data_construction/isoquant_s84_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)

# Extract the unique values in column Var1 from each data frame
var1_bambu <- unique(refFC_final_bambu$Var1)
var1_stringite <- unique(refFC_final_stringite$Var1)
var1_oarfish <- unique(refFC_final_oarfish$Var1)
var1_isoquant <- unique(refFC_final_isoquant$Var1)

# Combine all unique values
all_var1 <- c(var1_bambu, var1_stringite, var1_oarfish, var1_isoquant)

# Count how many times each Var1 appears
var1_counts <- table(all_var1)

# Keep Var1 values that appear at least twice
var1_at_least_two <- var1_counts[var1_counts >= 2]

# Create a new data frame containing these Var1 values and their counts
result_df <- data.frame(Var1 = names(var1_at_least_two), Freq = as.integer(var1_at_least_two))
# 使用 separate() 函数将 Var1 列拆分为 isoform 和 compare 两列
result_df <- separate(result_df, col = Var1, into = c("isoform", "compare"), sep = " ")
table(result_df$compare)

# Define the calculation function
## uchar
rse <- function(x) {
  if (length(x[!is.na(x)]) <= 1) {
    return(0)  # 如果数据点只有 1 个，直接返回 0，避免 NA
  }
  sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) * sqrt(length(x[!is.na(x)])))
}

# Calculate the metric for fold change (delta PSI)
calculate_mean_delta_psi_metrics <- function(var1, ref_mean_delta_psi_list) {
  mean_delta_psi_values <- c()
  p_vals <- c()
  
  for (ref_mean_delta_psi in ref_mean_delta_psi_list) {
    matching_row <- ref_mean_delta_psi[ref_mean_delta_psi$Var1 == var1, ]
    if (nrow(matching_row) > 0) {
      mean_delta_psi_values <- c(mean_delta_psi_values, matching_row$mean_delta_psi)
      p_vals   <- c(p_vals, matching_row$medianp)
    }
  }
  
  list(
    mean_delta_psi_mean = ifelse(length(mean_delta_psi_values) > 0, 
                                 mean(mean_delta_psi_values, na.rm = TRUE), NA),
    medianp  = ifelse(length(p_vals) > 0, median(p_vals,  na.rm = TRUE) , NA),
    uncertainty = ifelse(length(mean_delta_psi_values[!is.na(mean_delta_psi_values)]) > 0, 
                         rse(mean_delta_psi_values) * 100, NA)
  )
}

# Load all reference files
refFC_list <- list(refFC_final_bambu, refFC_final_stringite, refFC_final_oarfish, refFC_final_isoquant)

# Add metric columns (mean ΔPSI, median p-value, and uncertainty) to result_df
result_df <- result_df %>%
  rowwise() %>%
  mutate(
    mean_delta_psi_mean = calculate_mean_delta_psi_metrics(paste(isoform, compare), refFC_list)$mean_delta_psi_mean,
    medianp = calculate_mean_delta_psi_metrics(paste(isoform, compare), refFC_list)$medianp,
    uncertainty = calculate_mean_delta_psi_metrics(paste(isoform, compare), refFC_list)$uncertainty,
    isoform_compare = paste(isoform, compare)
  ) %>%
  ungroup()

# Inspect the results
print(head(result_df))
table(result_df$compare)

# Add annotations
# overall
# Read the local GTF file
gtf_file <- fread("/vast/projects/quartet_rna_refdata/reference/suppa2/gencode_v43.all_AS.ioe") %>% as.data.frame()
gtf_file <- gtf_file[,c("event_id","alternative_transcripts","total_transcripts")]
result_df <- merge(result_df,gtf_file,by.x = "isoform", by.y = "event_id")

# #uncertainty  25%、50% 
# uncertainty_q25 <- quantile(result_df$uncertainty, 0.25)
# uncertainty_q50 <- quantile(result_df$uncertainty, 0.50)
# uncertainty_q75 <- quantile(result_df$uncertainty, 0.75)
# 
# # PRINT
# print(c(Q25 = uncertainty_q25, Q50 = uncertainty_q50, Q75 = uncertainty_q75))

# define
# classify_isoforms <- function(Freq, uncertainty, q50, q75) {
#   if (Freq >= 3 && uncertainty <= q50) {
#     return("High Confidence")
#   } else if (Freq >= 2 && uncertainty > q50 && uncertainty <= q75) {
#     return("Medium Confidence")
#   } else {
#     return("Low Confidence")
#   }
# }
classify_ase <- function(Freq) {
  if (Freq >= 4) {
    return("High Confidence")
  } else if (Freq >= 3) {
    return("Medium Confidence")
  } else {
    return("Low Confidence")
  }
}

# Add a classification column to result_df
result_df <- result_df %>%
  mutate(
    classification = mapply(
      classify_ase, 
      Freq = Freq
    )
  )

# Inspect the classification results
table(result_df$classification)

# Save to another data frame
refFC_final_p4 <- result_df
colnames(refFC_final_p4)[1] <- "ASE"
write.csv(refFC_final_p4, "./ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv", row.names = FALSE)

# summary the reference datasets
# Define classification criteria
refFC_final_p4 <- refFC_final_p4 %>%
  mutate(
    Condition = case_when(
      classification == "High Confidence" ~ paste("Freq >= 4"),
      classification == "Medium Confidence" ~ paste("Freq >= 3"),
      classification == "Low Confidence" ~ "Others"
    )
  )

# Count the number of each category
classification_counts <- refFC_final_p4 %>%
  count(classification, Condition) %>%
  ungroup() %>%
  mutate(Percentage = n / sum(n) * 100)

# Create a descriptive table
description_table <- classification_counts %>%
  select(Classification = classification, Condition, Count = n, Percentage)

# view the results
print(description_table)

# Count the number of isoforms with uncertainty < 30
num_isoforms_below_30 <- sum(refFC_final_p4$uncertainty < 30, na.rm = TRUE)

# Calculate the total number of isoforms
total_isoforms <- nrow(refFC_final_p4)

# Calculate the proportion
percentage_below_30 <- (num_isoforms_below_30 / total_isoforms) * 100

# output
cat(sprintf("uncertainty 小于 30%% 的 isoform 占比: %.2f%%\n", percentage_below_30))

# remove all
rm(list = ls())
gc()
