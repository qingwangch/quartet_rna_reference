#' Metrics-statistic_p10
#' 2024-06-23
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
library(tidyr)
library(data.table)
library(purrr)

# data import
## isoform
files <- list(
  RC = "performance/performance_assesement_RC_results_LO_p10.csv",
  RMSE = "performance/performance_assesement_RMSE_results_LO_p10.csv",
  MCC = "performance/performance_assesement_DEG_MCC_results_LO_p10.csv",
  FNR = "performance/performance_assesement_DEG_MCC_results_LO_p10.csv"
)

# 整合数据
combined_data <- map_dfr(names(files), function(metric) {
  data <- read.csv(files[[metric]]) %>%
    select(batch_id, compare, pipeline, group, value = !!sym(metric)) %>%
    mutate(Metric = metric)
  return(data)
})

# 使用 strsplit 和 sapply 更新 Protocols_platforms 列
combined_data$protocols_platforms <- sapply(strsplit(combined_data$batch_id, "_"), function(x) {
  paste0(x[1], "_", x[2])
})

combined_data$labs <- sapply(strsplit(combined_data$batch_id, "_"), function(x) {
  x[3]
})

isoform_metrics <- combined_data
table(isoform_metrics$Metric)
isoform_metrics$type <- "isoform"
fwrite(isoform_metrics, "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_isoform.csv")

## AS
files <- list(
  RC = "performance/performance_assesement_AS_RC_results_LS_p10.csv",
  RMSE = "performance/performance_assesement_AS_RMSE_results_LS_p10.csv",
  MCC = "performance/performance_assesement_DAS_MCC_results_LS_p10.csv",
  FNR = "performance/performance_assesement_DAS_MCC_results_LS_p10.csv"
)

# 整合数据
combined_data <- map_dfr(names(files), function(metric) {
  data <- read.csv(files[[metric]]) %>%
    select(batch, compare, pipeline, group, value = !!sym(metric)) %>%
    mutate(Metric = metric)
  return(data)
})

# 使用 strsplit 和 sapply 更新 Protocols_platforms 列
combined_data$protocols_platforms <- sapply(strsplit(combined_data$batch, "_"), function(x) {
  paste0(x[1], "_", x[2])
})

combined_data$labs <- sapply(strsplit(combined_data$batch, "_"), function(x) {
  x[3]
})

ase_metrics <- combined_data
table(ase_metrics$Metric)
ase_metrics$type <- "AS"
colnames(ase_metrics)[1] <- "batch_id"

# merge all metrcis
all_metrics <- rbind(isoform_metrics,ase_metrics)
fwrite(all_metrics, "/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv")

metric_summary <- all_metrics %>%
  group_by(Metric, group, type) %>%     # 按 Metric、group、type 分组
  summarise(
    mean_value = mean(value, na.rm = TRUE),  # 对 value 列求均值
    sd_value   = sd(value, na.rm = TRUE),    # 对 value 列求标准差
    n          = n(),                        # 样本量
    .groups    = "drop"
  )

print(metric_summary)

metric_summary %>%
  filter(Metric == "RC",type == "isoform")

metric_summary %>%
  filter(Metric == "RMSE",type == "isoform")

metric_summary %>%
  filter(Metric == "MCC",type == "isoform")

metric_summary %>%
  filter(Metric == "FNR",type == "isoform")

metric_summary %>%
  filter(Metric == "RC",type == "AS")

metric_summary %>%
  filter(Metric == "RMSE",type == "AS")

metric_summary %>%
  filter(Metric == "MCC",type == "AS")

metric_summary %>%
  filter(Metric == "FNR",type == "AS")


# remove all
rm(list=ls())
gc()

