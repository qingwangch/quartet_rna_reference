#' isoform-ref-datasets-lo-qpcr-validation
#' Qingwang Chen
#' 2025-05-19

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
library(ComplexUpset)

# data import
############ refFC
refFC_final_p4 <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv")

### Count how many references come from SR and how many from LO
table(refFC_final_p4$source)

############ qPCR
ct_qpcr <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/qpcr/raw_ct_20250326_U.csv")

# Replace "Undetermined" with NA and convert to numeric
ct_qpcr <- ct_qpcr %>%
  mutate(across(starts_with("Value"), ~ ifelse(. == "Undetermined", NA, as.numeric(.))))

# check
str(ct_qpcr)
head(ct_qpcr)

calculate_log2FC_pval <- function(data, control_group = "D6") {
  
  # Get all experimental groups (excluding D6)
  experimental_groups <- unique(data$Sample)
  experimental_groups <- experimental_groups[experimental_groups != control_group]
  
  # Store results
  results <- list()
  
  for (exp_group in experimental_groups) {
    
    # Select data for control group (D6) and experimental group (e.g., D5, F7, M8)
    control_data <- data %>% 
      filter(Sample == control_group) %>% 
      select(Ensembl_Id, starts_with("Value")) 
    
    exp_data <- data %>% 
      filter(Sample == exp_group) %>% 
      select(Ensembl_Id, starts_with("Value"))
    
    # Ensure gene IDs match
    merged_data <- inner_join(
      control_data, exp_data,
      by = "Ensembl_Id", suffix = c("_control", "_exp")
    )
    
    # Check if merged_data is empty
    if (nrow(merged_data) == 0) {
      warning(
        paste0(
          "No matching Ensembl_Id found between control (",
          control_group, ") and experiment (", exp_group, "). Skipping..."
        )
      )
      next  # skip current experimental group
    }
    
    # Select `Value` columns for control and experimental groups
    control_values <- merged_data %>% select(contains("control"))
    exp_values     <- merged_data %>% select(contains("exp"))
    
    # Step 1: remove transcripts where all values are NA
    valid_rows <- rowSums(!is.na(control_values)) > 0 &
                  rowSums(!is.na(exp_values))   > 0
    merged_data    <- merged_data[valid_rows, ]
    control_values <- control_values[valid_rows, ]
    exp_values     <- exp_values[valid_rows, ]
    
    # Check if any data remain after filtering valid_rows
    if (nrow(merged_data) == 0) {
      warning(
        paste0("All values were filtered out for experiment (",
               exp_group, "). Skipping...")
      )
      next
    }
    
    # Step 2: compute log2FC (adjusted for qPCR Ct values)
    mean_Ct_control <- rowMeans(control_values, na.rm = TRUE)
    mean_Ct_exp     <- rowMeans(exp_values,   na.rm = TRUE)
    
    log2FC_values <- -(mean_Ct_exp - mean_Ct_control)  # -ΔCt for qPCR
    
    # Step 3: compute p-value (t-test), ensuring sufficient data
    p_values <- mapply(function(ctrl, exp) {
      # ensure at least two non-NA values
      if (length(na.omit(ctrl)) > 1 && length(na.omit(exp)) > 1) {
        t.test(ctrl, exp, paired = FALSE, var.equal = TRUE)$p.value
      } else {
        NA  # return NA if data are insufficient
      }
    }, 
    split(control_values, 1:nrow(control_values)), 
    split(exp_values,     1:nrow(exp_values)))
    
    # Ensure lengths of log2FC and p_value match merged_data
    if (length(log2FC_values) != nrow(merged_data)) {
      warning(paste0("Mismatch in log2FC length for ", exp_group, ". Adjusting..."))
      log2FC_values <- rep(NA, nrow(merged_data))
    }
    if (length(p_values) != nrow(merged_data)) {
      warning(paste0("Mismatch in p-value length for ", exp_group, ". Adjusting..."))
      p_values <- rep(NA, nrow(merged_data))
    }
    
    # Step 4: merge results
    result_df <- data.frame(
      isoform  = merged_data$Ensembl_Id,
      log2FC   = log2FC_values,
      p_value  = p_values,
      compare  = paste0(exp_group, "/", control_group)
    )
    
    results[[exp_group]] <- result_df
  }
  
  # Combine results from all experimental groups
  final_results <- bind_rows(results)
  return(final_results)
}

# Calculate log2FC and p-value
final_results <- calculate_log2FC_pval(ct_qpcr)
final_results$isoform_compare <- paste(final_results$isoform,final_results$compare)
final_results <- final_results %>%
  mutate(feature_type = case_when(
    log2FC > 0.5 & p_value < 0.05  ~ "up-regulate",
    log2FC < -0.5 & p_value < 0.05 ~ "down-regulate",
    TRUE ~ "non-DEI"
  ))

# view
head(final_results)
colnames(final_results)
# fwrite(final_results,"/vast/projects/quartet_rna_refdata/analysis/isoforms/qpcr/log2fc_20250409_U.csv")

# Merge the two data frames by isoform_compare
merged_data <- merge(
  final_results,
  refFC_final_p4,
  by = "isoform_compare",
  suffixes = c("_pcr", "_ref")
)

# Absolute difference between qPCR and reference log2FC
merged_data$delta <- abs(merged_data$log2FC_pcr-merged_data$log2FC_ref)
  
# Preview the merged data
head(merged_data)
# Compute correlation coefficient
correlation <- cor(merged_data$log2FC_pcr, merged_data$log2FC_ref, use = "complete.obs")
cat("Correlation coefficient (log2FC):", correlation, "\n")

# Scatter plot comparing log2FC values
table(merged_data$source)
p1 <- ggplot(merged_data, aes(x = log2FC_ref, y = log2FC_pcr, color = compare_pcr)) +
  geom_point(alpha = 0.7, size = 3) +                              # point transparency and size
  geom_smooth(method = "lm", color = "black", se = FALSE,
              linetype = "dashed") +                               # add linear regression line
  scale_color_viridis_d() +                                        # Viridis palette
  theme_minimal(base_size = 14) +                                  # minimal theme
  labs(
    title    = "Correlation of log2FC (Ref isoforms)",
    subtitle = paste0("N = ", nrow(merged_data)),
    x        = "LO reference",
    y        = "qPCR",
    color    = "Comparison"
  ) +
  annotate(
    "text",
    x      = min(merged_data$log2FC_ref, na.rm = TRUE),
    y      = max(merged_data$log2FC_ref, na.rm = TRUE),
    label  = paste0("r = ", round(correlation, 2)),
    hjust  = 0, vjust = 1, size = 5, fontface = "bold", color = "black"  # display correlation coefficient
  ) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, size = 14),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = c(0.85, 0.15),                                 # place legend inside the plot (upper-right)
    legend.background = element_rect(fill = "white",
                                     color = "black", size = 0.5),   # add legend background
    panel.border    = element_rect(color = "black", fill = NA, size = 1),  # add border
    panel.grid.major = element_blank(),                             # remove major grid lines
    panel.grid.minor = element_blank(),                             # remove minor grid lines
    axis.line       = element_line(color = "black"),                # add axis lines
    axis.ticks      = element_line(color = "black", size = 0.8),    # add axis ticks
    axis.ticks.length = unit(0.2, "cm")                             # set tick length
  )

# Print the plot
print(p1)

### RefDEI
refDEI <-read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_b7_p4_s84_u_20250516.csv")

# Inner-join on isoform_compare; add suffixes _pcr and _ref
common_df <- inner_join(final_results, refDEI, by = "isoform_compare", suffix = c("_pcr", "_ref"))

# Add a consistency flag
common_df <- common_df %>%
  mutate(consistency = if_else(feature_type == Final, "consistent", "inconsistent"))

# Count consistent vs. inconsistent cases
consistency_summary <- common_df %>%
  count(consistency)

print(consistency_summary)

# Keep rows where feature_type is consistent and not “non-DEI”
common_consistent <- common_df %>% filter(consistency == "consistent" & Final != "non-DEI")

# Correlation of log2FC values
log2fc_cor <- cor(common_consistent$log2FC_pcr, common_consistent$log2FC_ref, use = "complete.obs")
print(paste("Correlation coefficient:", round(log2fc_cor, 3)))

# Scatter plot comparing log2FC values
table(common_consistent$source)
p2 <- ggplot(
        common_consistent,
        aes(x = log2FC_ref, y = log2FC_pcr, color = compare_pcr)
      ) +
  geom_point(alpha = 0.7, size = 3) +                          # point transparency and size
  geom_smooth(method = "lm", color = "black", se = FALSE,
              linetype = "dashed") +                           # add linear regression line
  scale_color_viridis_d() +                                    # Viridis colour palette
  theme_minimal(base_size = 14) +                              # minimal theme
  labs(
    title    = "Correlation of log2FC (Ref DEIs)",
    subtitle = paste0("N = ", nrow(common_consistent)),
    x        = "LO reference",
    y        = "qPCR",
    color    = "Comparison"
  ) +
  annotate(
    "text",
    x      = min(common_consistent$log2FC_ref, na.rm = TRUE),
    y      = max(common_consistent$log2FC_ref, na.rm = TRUE),
    label  = paste0("r = ", round(log2fc_cor, 2)),
    hjust  = 0, vjust = 1, size = 5, fontface = "bold", color = "black"  # show correlation coefficient
  ) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, size = 14),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = c(0.85, 0.15),                           # place legend inside the plot (upper-right)
    legend.background = element_rect(fill = "white",
                                     color = "black", size = 0.5),
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line       = element_line(color = "black"),
    axis.ticks      = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")                        # tick length
  )

# Display the plot
print(p2)

# 3. Combine the two sub-plots (p1 and p2) with a shared legend
combined_plot <- plot_grid(
  plot_grid(
    p1, p2,
    align  = "h",
    ncol   = 2,
    labels = c("", "")      # arrange horizontally; adjust labels if you want "A", "B"
  )
)

# 4. Display the combined figure
print(combined_plot)

# save
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/validation_lo_ref_qpcr_isoforms_DEIs.pdf", plot = combined_plot, width = 10, height = 5)

# remove all
rm(list=ls())
gc()

