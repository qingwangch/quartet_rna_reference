#' AS-ref-datasets-d13-p4-u-so
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


# library import
library(data.table)
library(tidyr) 

# Ratio-based reference datasets
refFC_final_kallisto <- read.csv("./ref_data_construction/kallisto_s156_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_stringite <- read.csv("./ref_data_construction/stso_s156_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_salmon <- read.csv("./ref_data_construction/salmon_s156_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)
refFC_final_rsem <- read.csv("./ref_data_construction/rsem_s156_ref_expr_das.csv",header=T,stringsAsFactors=F,row.names=1,check.names=F)

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
table(result_df$compare)

# Define the calculation function
## Relative Standard Error (RSE)
rse <- function(x) {
  if (length(x[!is.na(x)]) <= 1) {
    return(0)  # If there is only one non-NA data point, return 0 to avoid NA
  }
  sd(x, na.rm = TRUE) / (mean(x, na.rm = TRUE) * sqrt(length(x[!is.na(x)])))
}

# Calculate the metric for fold change (delta PSI)
calculate_mean_delta_psi_metrics <- function(var1, ref_mean_delta_psi_list) {
  mean_delta_psi_values <- c()
  
  for (ref_mean_delta_psi in ref_mean_delta_psi_list) {
    matching_row <- ref_mean_delta_psi[ref_mean_delta_psi$Var1 == var1, ]
    if (nrow(matching_row) > 0) {
      mean_delta_psi_values <- c(mean_delta_psi_values, matching_row$mean_delta_psi)
    }
  }
  
  list(
    mean_delta_psi_mean = ifelse(length(mean_delta_psi_values) > 0, 
                                 mean(mean_delta_psi_values, na.rm = TRUE), NA),
    uncertainty = ifelse(length(mean_delta_psi_values[!is.na(mean_delta_psi_values)]) > 0, 
                         rse(mean_delta_psi_values) * 100, NA)
  )
}

# Load all reference files
refFC_list <- list(refFC_final_kallisto, refFC_final_stringite, refFC_final_salmon, refFC_final_rsem)

# Add fc, log2FC, and uncertainty columns to result_df
result_df <- result_df %>%
  rowwise() %>%
  mutate(
    mean_delta_psi_mean = calculate_mean_delta_psi_metrics(paste(isoform, compare), refFC_list)$mean_delta_psi_mean,
    uncertainty = calculate_mean_delta_psi_metrics(paste(isoform, compare), refFC_list)$uncertainty,
    isoform_compare = paste(isoform, compare)
  ) %>%
  ungroup()

# Inspect the results
print(head(result_df))
table(result_df$compare)

# Add annotations
# overall
# Read the local GTF-derived IOE file
gtf_file <- fread("/vast/projects/quartet_rna_refdata/reference/suppa2/gencode_v43.all_AS.ioe") %>% as.data.frame()
gtf_file <- gtf_file[,c("event_id","alternative_transcripts","total_transcripts")]
result_df <- merge(result_df,gtf_file,by.x = "isoform", by.y = "event_id")

classify_ase <- function(Freq) {
  if (Freq >= 4) {
    return("High Confidence")
  } else if (Freq >= 3) {
    return("Medium Confidence")
  } else {
    return("Low Confidence")
  }
}

# Merge transcript information into result_df
result_df <- result_df %>%
  mutate(
    classification = mapply(
      classify_ase, 
      Freq = Freq
    )
  )

# Check the classification results
table(result_df$classification)

# Save the data frame to another object and CSV file
refFC_final_p4 <- result_df
colnames(refFC_final_p4)[1] <- "ASE"
write.csv(refFC_final_p4, "./ref_data_construction/so-e-g/ref_expr_as_b13_p4_s156_u_20250910.csv", row.names = FALSE)

# Summarize the reference data sets
# Define the textual condition for each class
refFC_final_p4 <- refFC_final_p4 %>%
  mutate(
    Condition = case_when(
      classification == "High Confidence" ~ paste("Freq >= 4"),
      classification == "Medium Confidence" ~ paste("Freq >= 3"),
      classification == "Low Confidence" ~ "Others"
    )
  )

# Count the number of isoforms in each class
refFC_final_p4 <- as.data.table(refFC_final_p4)

classification_counts <- refFC_final_p4[ , .N, by = .(classification, Condition)]
classification_counts[ , Percentage := N / sum(N) * 100]
classification_counts

classification_counts <- refFC_final_p4 %>%
  count(classification, Condition) %>%
  ungroup() %>%
  mutate(Percentage = n / sum(n) * 100)


# Count the number of isoforms with uncertainty < 30
num_isoforms_below_30 <- sum(refFC_final_p4$uncertainty < 30, na.rm = TRUE)

# Total number of isoforms
total_isoforms <- nrow(refFC_final_p4)

# Calculate the proportion
percentage_below_30 <- (num_isoforms_below_30 / total_isoforms) * 100

# Output the result
cat(sprintf("Proportion of isoforms with uncertainty < 30%%: %.2f%%\n", percentage_below_30))

# remove all
rm(list = ls())
gc()
