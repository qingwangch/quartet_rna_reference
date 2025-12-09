#' junctions_per_read_SR_multi
#' Qingwang Chen
#' 2025-07-09

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/star_bam/")

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())


# library import
library(megadepth)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# define the function
process_bam_file <- function(test_bam) {
  # 运行 bam_to_junctions()
  example_jxs <- bam_to_junctions(test_bam, all_junctions = TRUE, overwrite = TRUE)
  
  # 读取数据
  example_jxs <- read_junction_table(example_jxs)
  
  # 使用 data.table 来统计每个 read_name 的 junction 数量
  # .N 是一个特殊的符号，表示当前组的大小
  example_jxs_dt <- as.data.table(example_jxs)
  example_jxs_dt <- example_jxs_dt[unique > 0]
  junction_count <- example_jxs_dt[, .(junctions_covered = .N), by = read_name]
  
  # 计算每个 read_name 对应的 junctions_covered 占比
  junction_count_percentage <- junction_count[, .(n = .N), by = junctions_covered]
  junction_count_percentage[, percentage := n / sum(n) * 100]
  
  # 提取路径中的文件夹名称
  sample_id <- gsub(".downsampled_5000000_reads.bam","",test_bam)
  
  # 添加样本 ID 列
  junction_count_percentage$sample_id <- sample_id
  
  return(junction_count_percentage)
}

# 假设有多个样本路径，按批次处理
sample_bam_paths <- list.files("/vast/projects/quartet_rna_refdata/analysis/star_bam/", 
                               pattern = "*downsampled_5000000_reads.bam$", recursive = T)

# 分批处理，每次处理 10 个文件
batch_size <- 10
n_batches <- ceiling(length(sample_bam_paths) / batch_size)

all_results <- list()

for (i in 1:n_batches) {
  # 计算当前批次的文件范围
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(sample_bam_paths))
  
  # 获取当前批次的文件
  batch_files <- sample_bam_paths[start_idx:end_idx]
  print(batch_files)
  
  # 处理当前批次的文件
  batch_results <- lapply(batch_files, process_bam_file)
  
  # 合并当前批次的结果
  all_results[[i]] <- rbindlist(batch_results)
}

# 合并所有批次的结果
final_result <- rbindlist(all_results)
final_result$data_type <- "short read"

# 查看最终结果
head(final_result)

# save
fwrite(final_result,"./junction/quartet_SR_junctions_s60.csv")

# remove all
rm(list=ls())
gc()
