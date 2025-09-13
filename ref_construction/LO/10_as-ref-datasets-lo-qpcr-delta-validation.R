#' as-ref-datasets-lo-qpcr-delta-validation
#' Qingwang Chen
#' 2025-05-20

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
library(stringr)

# data import
############ refFC
refFC_final_p4 <- read.csv("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/ref_expr_as_b7_p4_s84_u_20250520.csv")
refFC_final_p4 <- refFC_final_p4 %>%
  mutate(as_type = str_extract(ASE, "(?<=;)[^:;]+")) %>% 
  rename(ASE_compare = isoform_compare)
head(refFC_final_p4)

### Count how many references come from SR and how many from LO
table(refFC_final_p4$source)

## Verify UUID
# refFC_final_p4$ATC <- paste(
#   refFC_final_p4$as_type,
#   refFC_final_p4$alternative_transcripts,
#   refFC_final_p4$total_transcripts,
#   refFC_final_p4$compare
# )
# # Alternative option:
# # refFC_final_p4$ATC <- paste(refFC_final_p4$ASE, refFC_final_p4$compare)
# 
# sum(duplicated(refFC_final_p4$ATC))

############ qPCR
# ct_qpcr <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/qpcr/delta_ct_as_20250414_U.csv")
ct_qpcr <- read.csv("/vast/projects/quartet_rna_refdata/analysis/isoforms/qpcr/delta_ct_as_20250513_U.csv")

# check
str(ct_qpcr)
head(ct_qpcr)
ct_qpcr <- ct_qpcr %>% filter(!is.na(ASE) & ASE != "")

calculate_psi_qpcr_v2 <- function(ct_data) {
  library(dplyr)
  library(tidyr)
  
  # Step 0: create a unique identifier ASE_compare
  ct_data <- ct_data %>%
    mutate(ASE_compare = paste0(ASE, " ", GroupA, "/", GroupB))
  
  # Step 1: reshape to the replicate level; compute expression and PSI
  psi_replicate_df <- ct_data %>%
    filter(Alternative_or_Total %in% c("A", "T")) %>%
    pivot_longer(cols = starts_with("Ct"), names_to = "Replicate", values_to = "Ct") %>%
    mutate(
      Sample = ifelse(grepl("^Ct_", Replicate), "GroupB", "GroupA"),
      Rep = gsub("^Ct_?", "", Replicate)
    ) %>%
    select(ASE_compare, Event_Type, Gene_Id, Sample, Rep, Alternative_or_Total, Ct) %>%
    pivot_wider(names_from = Alternative_or_Total, values_from = Ct) %>%
    filter(!is.na(A) & !is.na(T)) %>%
    mutate(
      expr_A = 2^(A),
      expr_T = 2^(T),
      PSI = expr_A / (expr_A + expr_T)
    )
  
  # Step 2: calculate the mean PSI for each group and ΔPSI
  psi_summary <- psi_replicate_df %>%
    group_by(ASE_compare, Event_Type, Gene_Id, Sample) %>%
    summarise(
      PSI_mean = mean(PSI, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = Sample, values_from = PSI_mean, names_prefix = "PSI_") %>%
    mutate(
      delta_PSI = PSI_GroupA - PSI_GroupB
    )
  
 # Step 3: unpaired t-test
  pval_df <- psi_replicate_df %>%
    group_by(ASE_compare, Event_Type, Gene_Id) %>%
    summarise(
      p_value = tryCatch(
        t.test(PSI[Sample == "GroupA"], PSI[Sample == "GroupB"], var.equal = TRUE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  
  # Step 4: merge the results
  final_df <- psi_summary %>%
    left_join(pval_df, by = c("ASE_compare", "Event_Type", "Gene_Id")) %>%
    select(ASE_compare, Event_Type, Gene_Id, PSI_GroupA, PSI_GroupB, delta_PSI, p_value)
  
  return(final_df)
}

# Calculate log2FC and p-value (actually ΔPSI and p-value here)
final_results <- calculate_psi_qpcr_v2(ct_qpcr)
final_results <- final_results %>%
  mutate(feature_type = case_when(
    delta_PSI > 0.05 & p_value < 0.05  ~ "up-regulate",
    delta_PSI < -0.05 & p_value < 0.05 ~ "down-regulate",
    TRUE ~ "non-DAS"
  ))
# final_results <- final_results %>%
#   mutate(feature_type = case_when(
#     delta_PSI > 0.05 ~ "up-regulate",
#     delta_PSI < -0.05 ~ "down-regulate",
#     TRUE ~ "non-DEI"
#   ))

# Preview the first few rows
head(final_results)

fwrite(final_results,"isoforms/qpcr/deltapsi_20250513_U.csv")

colnames(final_results)
unique_in_final <- anti_join(final_results, refFC_final_p4, by = "ASE_compare")

# Merge the qPCR results with reference data by ASE_compare
merged_data <- merge(
  final_results,
  refFC_final_p4,
  by = "ASE_compare",
  suffixes = c("_pcr", "_ref")
)

# Inspect the merged data frame
head(merged_data)
merged_data$delta <- abs(merged_data$delta_PSI-merged_data$mean_delta_psi_mean)
fwrite(merged_data,"/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/qpcr_ase_ref_m_20250520.csv")

# remove events
remove_evnets <- merged_data %>% 
  group_by(ASE_compare) %>%                 
  summarise(max_delta = max(delta, na.rm = TRUE), .groups = "drop") %>% 
  arrange(desc(max_delta)) %>%              # max to min
  slice_head(n = 6)    

# filter
merged_data <- merged_data %>%
  filter(!ASE_compare %in% remove_evnets$ASE_compare)

# cal
correlation <- cor(merged_data$delta_PSI, merged_data$mean_delta_psi_mean, use = "complete.obs")
cat("Correlation coefficient (log2FC):", correlation, "\n")

# plot
table(merged_data$source)
p1 <- ggplot(merged_data, aes(x = mean_delta_psi_mean, y = delta_PSI, color = compare)) +
  geom_point(alpha = 0.7, size = 3) +                          # point transparency and size
  geom_smooth(method = "lm", color = "black", se = FALSE,
              linetype = "dashed") +                           # add linear regression line
  scale_color_viridis_d() +                                    # Viridis colour palette
  theme_minimal(base_size = 14) +                              # minimal theme
  labs(
    title    = "Correlation of ΔPSI",
    subtitle = paste0("N = ", nrow(merged_data)),
    x        = "Reference ASEs (LO)",
    y        = "qPCR",
    color    = "Comparison"
  ) +
  annotate(
    "text",
    x      = min(merged_data$mean_delta_psi_mean, na.rm = TRUE),
    y      = max(merged_data$delta_PSI, na.rm = TRUE),
    label  = paste0("r = ", round(correlation, 2)),
    hjust  = 0, vjust = 1, size = 5, fontface = "bold", color = "black"
  ) +                                                          # display correlation coefficient
  theme(
    plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, size = 14),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = c(0.85, 0.15),                           # place legend inside plot
    legend.background = element_rect(fill = "white",
                                     color = "black", size = 0.5),  # legend background
    panel.border    = element_rect(color = "black", fill = NA, size = 1),  # border
    panel.grid.major = element_blank(),                        # remove major grid lines
    panel.grid.minor = element_blank(),                        # remove minor grid lines
    axis.line       = element_line(color = "black"),           # add axis lines
    axis.ticks      = element_line(color = "black", size = 0.8),  # axis ticks
    axis.ticks.length = unit(0.2, "cm")                        # tick length
  )

# Print the plot
print(p1)

### RefDEI
refDAS <- read.csv("/vast/projects/quartet_rna_refdata/analysis/suppa2/ref_data_construction/lo-e-g/RefData_AS_all_DAS_classified_u_0520.csv")
head(refDAS)
refDAS <- refDAS %>%
  dplyr::rename(ASE_compare = isoform_compare)
# table(refDAS$source)
# fwrite(refDAS,"/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/RefData_DAS_LS_t2_u_20250408.csv")

# Inner-join on ASE_compare; add suffixes _pcr and _ref for feature_type
final_results <- fread("/vast/projects/quartet_rna_refdata/analysis/isoforms/qpcr/deltapsi_20250513_U.csv")
common_df <- inner_join(final_results, refDAS, by = "ASE_compare", suffix = c("_pcr", "_ref"))

# Add a consistency flag
common_df <- common_df %>%
  mutate(consistency = if_else(feature_type_pcr == feature_type_ref, "consistent", "inconsistent"))

# Count the numbers of consistent and inconsistent cases
consistency_summary <- common_df %>%
  count(consistency)

print(consistency_summary)

# Keep rows whose feature_type is consistent and not “non-DAS”
common_consistent <- common_df %>% filter(consistency == "consistent" & final != "non-DAS")
# common_consistent <- common_df %>% filter(consistency == "consistent")

common_inconsistent <- common_df %>% filter(consistency == "inconsistent")

# Correlation of ΔPSI between qPCR and reference
log2fc_cor <- cor(common_consistent$delta_PSI, common_consistent$mean_delta_psi_mean, use = "complete.obs")
print(paste("Correlation coefficient:", round(log2fc_cor, 3)))

# Quick check of reference source distribution
table(common_consistent$source)
# Scatter plot comparing ΔPSI values
p2 <- ggplot(
        common_consistent,
        aes(x = mean_delta_psi_mean, y = delta_PSI, color = compare)
      ) +
  geom_point(alpha = 0.7, size = 3) +                         # point transparency & size
  geom_smooth(method = "lm", color = "black", se = FALSE,
              linetype = "dashed") +                          # linear-fit line
  scale_color_viridis_d() +                                   # Viridis palette
  theme_minimal(base_size = 14) +                             # minimal theme
  labs(
    title    = "Correlation of ΔPSI",
    subtitle = paste0("N = ", nrow(common_consistent)),
    x        = "Reference DAS (LO)",
    y        = "qPCR",
    color    = "Comparison"
  ) +
  annotate(
    "text",
    x      = min(common_consistent$mean_delta_psi_mean, na.rm = TRUE),
    y      = max(common_consistent$delta_PSI, na.rm = TRUE),
    label  = paste0("r = ", round(log2fc_cor, 2)),
    hjust  = 0, vjust = 1, size = 5, fontface = "bold", color = "black"
  ) +                                                          # correlation coefficient
  theme(
    plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, size = 14),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = c(0.85, 0.15),                           # put legend inside plot (upper-right)
    legend.background = element_rect(fill = "white",
                                     color = "black", size = 0.5),
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line       = element_line(color = "black"),
    axis.ticks      = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p2)

# 3. Combine sub-plots and a shared legend  (implementation depends on patchwork / cowplot etc.)
combined_plot <- plot_grid(
  plot_grid(
    p1, p2,
    align = "h", ncol = 2, labels = c("", "")  # arrange the sub-plots horizontally and label them A and B
  )
)

# 4. Display the combined plot
print(combined_plot)

ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/validation_ls_ref_qpcr_ASE_DAS.pdf", plot = combined_plot, width = 12, height = 6)

# remove all
rm(list=ls())
gc()

