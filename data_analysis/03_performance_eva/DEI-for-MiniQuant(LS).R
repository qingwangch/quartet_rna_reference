#' DEI-for-MiniQuant(LS)
#' Qingwang Chen
#' 2025-06-23

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# 指定新的库路径
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# 修改 .libPaths，使其优先使用自定义路径
.libPaths(c(custom_lib_path, .libPaths()))

# 检查当前的库路径
print(.libPaths())

# Import libraries
library(edgeR)
library(dplyr)
library(tibble)

# Define function: DE analysis
run_DE_analysis <- function(counts, group1, group2, batch_id, pipeline, analysis_type = "isoform") {
  # Parameters:
  # counts: Expression count matrix
  # group1: Control group name (e.g., "D6")
  # group2: Experimental group name (e.g., "D5")
  # batch_id: Current batch identifier
  # pipeline: Pipeline name (e.g., "Bambu", "StringTie")
  # analysis_type: Type of analysis ("isoform" or "gene")
  
  # Extract group information
  group_labels <- sapply(strsplit(colnames(counts), "_"), function(x) x[5])
  group <- factor(group_labels)
  
  # Filter samples for group1 and group2
  selected_samples <- which(group %in% c(group1, group2))
  counts_filtered <- counts[, selected_samples]
  group_filtered <- factor(group[selected_samples])
  
  # Create DGEList object
  y <- DGEList(counts = counts_filtered, group = group_filtered)
  
  # Filter low-expression genes/isoforms
  keep <- filterByExpr(y, min.count = 1) # For lrRNA-seq
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # Normalize and estimate dispersion
  y <- calcNormFactors(y)
  design <- model.matrix(~0 + group_filtered)
  colnames(design) <- levels(group_filtered)  # Name columns by group levels
  y <- estimateDisp(y, design)
  
  # Define contrasts
  contrast_matrix <- makeContrasts(
    contrasts = paste0(group2, "-", group1),
    levels = design
  )
  
  # Print the contrast for verification
  cat("Contrast being tested:", paste0(group2, "-", group1), "\n")
  print(contrast_matrix)
  
  # Fit GLM model and perform differential expression testing
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, contrast = contrast_matrix[, 1])  # Compare group2 against group1
  
  # Retrieve results
  results <- topTags(qlf, n = Inf)
  results_df <- results$table
  results_df$batch_id <- batch_id
  results_df$comparison <- paste0(group2, "vs", group1)
  results_df$feature <- rownames(results_df)
  results_df$pipeline <- pipeline
  
  # Filter significant results (FDR < 0.05)
  significant <- results_df[results_df$FDR < 0.05, ]
  
  # Print summary
  cat("Pipeline:", pipeline, "\n")
  cat("Batch ID:", batch_id, "\n")
  cat("Analysis Type:", analysis_type, "\n")
  cat("Comparison:", group2, "vs", group1, "\n")
  cat("Number of significantly differentially expressed", analysis_type, ":", nrow(significant), "\n")
  
  # Return results
  return(results_df)
}

# Function to save results
save_results <- function(results, output_file) {
  write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Results saved to:", output_file, "\n")
}

# Function to preprocess counts and extract batch IDs
preprocess_counts <- function(file_path) {
  counts <- readRDS(file_path) %>% as.data.frame()
  counts[is.na(counts)] <- 0
  batch_ids <- unique(sapply(strsplit(colnames(counts), "_"), function(x) paste(x[1:4], collapse = "_")))
  return(list(counts = counts, batch_ids = batch_ids))
}

# Main function to perform DE analysis across batches and groups
perform_DE_analysis <- function(input_file, group1, group2_list, pipeline, output_dir, analysis_type = "isoform") {
  # Preprocess counts and extract batch IDs
  data <- preprocess_counts(input_file)
  counts <- data$counts
  batch_ids <- data$batch_ids
  
  # Initialize result container
  all_results <- data.frame()
  
  # Loop through each batch and perform DE analysis
  for (batch_id in batch_ids) {
    # Filter counts for current batch
    batch_counts <- counts[, grepl(batch_id, colnames(counts), fixed = TRUE)]
    
    # Loop through experimental groups
    for (group2 in group2_list) {
      cat("Processing pipeline:", pipeline, "Batch:", batch_id, "Comparison:", group2, "vs", group1, "\n")
      
      # Run DE analysis
      results <- run_DE_analysis(
        counts = batch_counts, 
        group1 = group1, 
        group2 = group2, 
        batch_id = batch_id, 
        pipeline = pipeline,
        analysis_type = analysis_type
      )
      
      # Combine results
      all_results <- bind_rows(all_results, results)
    }
  }
  
  # Add metadata columns
  all_results$analysis_type <- analysis_type
  all_results$feature_type <- "non-DE"
  all_results$feature_type[intersect(which(all_results$FDR < 0.05), which(all_results$logFC >= 1))] <- "up-regulate"
  all_results$feature_type[intersect(which(all_results$FDR < 0.05), which(all_results$logFC <= -1))] <- "down-regulate"
  
  # Save final results
  output_file <- file.path(output_dir, paste0("DE_results_", pipeline, "_", analysis_type, ".txt"))
  save_results(all_results, output_file)
}

miniquant_s144 <- fread("/vast/projects/quartet_rna_refdata/slurm/miniQuant/gencode/quartet_s144_longRead_matrix.tsv") %>% as.data.frame()
bambu_s84 <- readRDS("/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-minimap2-bambu-g-txq-count-s84t274031.rds")
samples <- colnames(bambu_s84)
miniquant_s84 <- miniquant_s144[,c("Transcript_id",samples)]

miniquant_s84 <- miniquant_s84 %>% 
  mutate(ENST_id = sub("^(ENST[0-9]+\\.[0-9]+).*", "\\1", Transcript_id))

# ---- 2. 按 ENST_id 合并，取每列（除了 ENST_id）最大值 ----
miniquant_s84_agg <- miniquant_s84 %>% 
  select(-Transcript_id) %>%                 # 去掉原 Transcript_id
  group_by(ENST_id) %>%                      # 以 ENST_id 分组
  summarise(across(everything(), max), .groups = "drop")  # 对其余列取最大值

# ---- 3. 将 ENST_id 设为行名（若需要） ----
miniquant_s84_mat <- miniquant_s84_agg %>% 
  column_to_rownames("ENST_id")              # 现在行名唯一，不会报错
dim(miniquant_s84_mat)
saveRDS(miniquant_s84_mat,"/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-miniquant-g-txq-count-s84t252739.rds")

# isoform
input_files <- list(
  list(
    file_path = "/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-miniquant-g-txq-count-s84t252739.rds",
    pipeline = "MiniQuant"
  )
)

group1 <- "D6"
group2_list <- c("D5", "F7", "M8")
output_dir <- "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/D7"
analysis_type <- "isoform"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Batch processing
for (input in input_files) {
  input_file <- input$file_path
  pipeline <- input$pipeline
  perform_DE_analysis(input_file, group1, group2_list, pipeline, output_dir, analysis_type)
}

# Clean up
rm(list = ls())
gc()

