#' generalized-isoform-reference-dataset-constrcution-short-read
#' Qingwang Chen
#' Date: 2025-07-09

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

# Process metadata
process_metadata <- function(expr_mat) {
  metadata <- data.frame(
    sample_id = colnames(expr_mat),
    lib = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[1]),
    platform = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[2]),
    lab = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[3]),
    Batch = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[4]),
    sample = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[5]),
    rep = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[6]),
    caller = sapply(strsplit(colnames(expr_mat), "_"), function(x) x[7])
  )
  metadata$batch_id <- paste(metadata$lib, metadata$platform, metadata$lab, metadata$Batch, sep = "_")
  metadata$Batch <- metadata$batch_id
  metadata$code <- metadata$sample_id
  return(metadata)
}

# Function to generate reference datasets
generate_reference_dataset <- function(expr_file, batch_for_construction,output_prefix, dataset_support_count, tool_name, dei_file = NULL) {
  # Load expression matrix
  expr.mat.tx <- readRDS(expr_file)
  expr.mat.tx[is.na(expr.mat.tx)] <- 0
  
  # Load metadata
  metadata <- process_metadata(expr.mat.tx)
  metadata$code <- metadata$sample_id
  metadata$Batch <- metadata$batch_id
  # batch_for_construction <- c("D_ONT_LG_B1", "P_ONT_LG_B1", "P_ONT_LN_B2", "M_PAB_LN_B1", "M_PAB_LG_B2")
  metadata <- metadata[metadata$batch_id %in% batch_for_construction,]
  
  # Filter expression data
  expr.mat.tx <- expr.mat.tx[, metadata$sample_id]
  
  # Process metadata
  metadata$batch_sample <- paste(metadata$batch_id, metadata$sample, sep = "_")
  meta_1 <- metadata
  ###################----------Char_detect-----------##############
  # Melt expression matrix
  exprMat_count_M <- reshape2::melt(as.matrix(expr.mat.tx))
  colnames(exprMat_count_M)[1:3] <- c("isoform", "library", "count")
  
  # Merge metadata
  exprMat_count_M$batch_sample <- meta_1$batch_sample[match(exprMat_count_M$library, meta_1$sample_id)]
  
  # Filter isoforms with counts >= 10
  expd_M <- exprMat_count_M[exprMat_count_M$count >= 10,]
  
  # Calculate detected isoforms per sample
  expd_pg <- data.frame(tapply(expd_M$isoform, as.factor(paste(expd_M$isoform, expd_M$batch_sample)), length))
  colnames(expd_pg) <- c("Num")
  expd_pg$isoform <- sapply(strsplit(as.character(rownames(expd_pg)), " "), function(x) x[1])
  expd_pg$batch_sample <- sapply(strsplit(as.character(rownames(expd_pg)), " "), function(x) x[2])
  expd_pg$batch <- meta_1$batch_id[match(expd_pg$batch_sample, meta_1$batch_sample)]
  expd_pg$sample <- meta_1$sample[match(expd_pg$batch_sample, meta_1$batch_sample)]
  
  # Filter isoforms detected in multiple samples
  expd_pg_f <- expd_pg[expd_pg$Num > 1,]
  table(expd_pg_f$batch,expd_pg_f$sample)
  # Save detected genes across samples
  write.csv(expd_pg_f, paste0(output_prefix, "_detect_genelist_acrossSample.csv"))
  col_means <- apply(table(expd_pg_f$batch,expd_pg_f$sample), 2, mean)
  print("Detected isoforms in any batch:")
  print(col_means)
  
  # Detect isoforms in all batches
  expd_refs <- lapply(unique(meta_1$sample), function(sample) {
    subset <- expd_pg_f[expd_pg_f$sample == sample,]
    expd_sample <- data.frame(table(subset$isoform))
    colnames(expd_sample) <- c("isoform", "freq")
    expd_sample$sample <- sample
    return(expd_sample)
  }) %>%
    bind_rows()
  
  # Filter isoforms detected in all batches
  expd_refs_f <- expd_refs[expd_refs$freq == length(unique(meta_1$batch_id)),]
  print("Detected isoforms by consensus:")
  print(table(expd_refs_f$sample))
  # Save results for isoforms detected across batches
  write.csv(expd_refs_f, paste0(output_prefix, "_detect_genelist_acrossBatch.csv"))
  # Calculate the number of genes detected in each sample
  alldet <- names(which(table(expd_refs_f$isoform) == 4));length(alldet)
  
  # DEI Analysis (if DEI file is provided)
  if (!is.null(dei_file) && file.exists(dei_file)) {
    DEIs <- fread(dei_file) %>% as.data.frame()
    DEIs <- DEIs[DEIs$batch_id %in% batch_for_construction,]
    DEIs <- DEIs %>%
      mutate(
        groupa = sapply(strsplit(comparison, "vs"), function(x) x[2]),
        groupb = sapply(strsplit(comparison, "vs"), function(x) x[1]),
        compare = paste(groupb, "/", groupa, sep = "")
      )
    
    meta <- metadata
    sample_combn <- data.frame(
      sampleA=rep("D6",3),
      sampleB=c("D5","F7","M8")
    )
    
    sample_combn$compare<-paste(sample_combn$sampleB,"/",sample_combn$sampleA,sep="")
    
    # define sample pair
    detect_isoforms_pairs <- data.frame(
      isoform = character(),
      compare = character(),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:nrow(sample_combn)) {
      s <- intersect(
        expd_refs_f$isoform[expd_refs_f$sample == sample_combn$sampleA[i]],
        expd_refs_f$isoform[expd_refs_f$sample == sample_combn$sampleB[i]]
      )
      
      # 如果没有交集，也可以跳过或记录空
      if (length(s) > 0) {
        detect_isoforms_pairs <- rbind(
          detect_isoforms_pairs,
          data.frame(
            isoform = s,
            compare = as.character(sample_combn$compare[i]),
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    # colnames(detect_isoforms_pairs)<-c("isoform","compare")
    detect_isoforms_pairs<-data.frame(detect_isoforms_pairs)
    
    detect_isoforms_pairs$isoform_compare<-paste(detect_isoforms_pairs$isoform,detect_isoforms_pairs$compare)
    print("Detected in both groups in each sample pair:")
    print(table(detect_isoforms_pairs$compare))
    
    # FC and P
    DEIs_p <- DEIs
    # DEIs_p$protocol <- sapply(strsplit(as.character(DEIs_p$batch),"_"),function(x){x[1]})
    DEIs_p$isoform_compare <- paste(DEIs_p$feature,DEIs_p$compare)
    
    DEIs_p_dec <- DEIs_p[DEIs_p$isoform_compare %in% detect_isoforms_pairs$isoform_compare,]
    
    # print(table(DEIs_p_dec$batch,DEIs_p_dec$compare))
    
    # #filter PValue<0.05
    DEIs_p_f <- DEIs_p_dec[DEIs_p_dec$PValue<0.05,]
    # DEIs_p_f <- DEIs_p_dec[abs(DEIs_p_dec$logfc)>1,]
    print("Detected DEIs in each batch:")
    print(table(DEIs_p_f$batch,DEIs_p_f$compare))
    
    DEIs_p_cal <- data.frame(table(paste(DEIs_p_f$feature,DEIs_p_f$compare)))
    DEIs_p_cal$isoform <- sapply(strsplit(as.character(DEIs_p_cal$Var1)," "),function(x){x[1]})
    DEIs_p_cal$compare <- sapply(strsplit(as.character(DEIs_p_cal$Var1)," "),function(x){x[2]})
    
    DEIs_p_cal_f <- DEIs_p_cal[DEIs_p_cal$Freq >= dataset_support_count,]
    DEIs_p_cal_f_fd <- DEIs_p_cal_f[DEIs_p_cal_f$Var1 %in% detect_isoforms_pairs$isoform_compare,]
    print("Differentially expressed by consensus:")
    print(table(DEIs_p_cal_f_fd$compare))
    table(DEIs_p_cal_f_fd$compare)/nrow(expr.mat.tx) * 100
    
    DEIs_p_cal_f_fd <- DEIs_p_cal_f_fd[order(DEIs_p_cal_f_fd$compare),]
    DEIs_p_cal_f_fd$isoform_compare <- paste(DEIs_p_cal_f_fd$isoform,DEIs_p_cal_f_fd$compare)
    
    ####meanFC medianp,seFC
    DEIs_p$isoform_compare <- paste(DEIs_p$feature,DEIs_p$compare)
    
    #function #
    isoform_compare=levels(as.factor(DEIs_p$isoform_compare))
    
    mm.DEG<-data.frame(
      meanlogFC=tapply(DEIs_p$logFC,as.factor(DEIs_p$isoform_compare),mean),
      medianp=tapply(DEIs_p$PValue,as.factor(DEIs_p$isoform_compare),median),
      isoform_compare=levels(as.factor(DEIs_p$isoform_compare))
    )
    
    # ref_FC_f2$fc <- 2^(mm.DEG$meanlogFC[match(ref_FC_f2$isoform_compare,mm.DEG$isoform_compare)])
    # ref_FC_f2$medianp <- mm.DEG$medianp[match(ref_FC_f2$isoform_compare,mm.DEG$isoform_compare)]
    DEIs_p_cal_f_fd$fc <- 2^(mm.DEG$meanlogFC[match(DEIs_p_cal_f_fd$isoform_compare,mm.DEG$isoform_compare)])
    DEIs_p_cal_f_fd$medianp <- mm.DEG$medianp[match(DEIs_p_cal_f_fd$isoform_compare,mm.DEG$isoform_compare)]
    
    # Output results to CSV file
    # output_file4 <- paste0(opt$out_dir, "/ref_expr_", Sys.Date(), ".csv")
    # write.csv(ref_FC_f2, output_file4)
    write.csv(DEIs_p_cal_f_fd, paste0(output_prefix, "_ref_expr.csv"))
    
    ###################----------char_refDEI-----------##############
    ## preprocess
    DEIs_p$protocol<-sapply(strsplit(as.character(DEIs_p$batch),"_"),function(x){x[1]})
    DEIs_p$isoform_compare<-paste(DEIs_p$feature,DEIs_p$compare)
    
    DEIs_p$DEI_type<-"non-DEI"
    DEIs_p$DEI_type[intersect(which(DEIs_p$PValue<0.05),which(DEIs_p$logFC>=1))]<-"up-regulate"
    DEIs_p$DEI_type[intersect(which(DEIs_p$PValue<0.05),which(DEIs_p$logFC<=(-1)))]<-"down-regulate"
    
    ########DEIs_cal 
    DEI_ref_cal<-data.frame(
      N_up=tapply(DEIs_p$DEI_type,as.factor(DEIs_p$isoform_compare),function(x){length(which(x=="up-regulate"))}),
      N_non=tapply(DEIs_p$DEI_type,as.factor(DEIs_p$isoform_compare),function(x){length(which(x=="non-DEI"))}),
      N_down=tapply(DEIs_p$DEI_type,as.factor(DEIs_p$isoform_compare),function(x){length(which(x=="down-regulate"))}))
    
    DEI_ref_cal$Final<-"non-DEI"
    DEI_ref_cal$Final[intersect(which(DEI_ref_cal$N_up>=1),which(DEI_ref_cal$N_down>=1))]<-"conflicting"
    
    DEI_ref_cal$Final[intersect(which(DEI_ref_cal$N_up>=dataset_support_count),which(DEI_ref_cal$N_down==0))]<-"up-regulate"
    DEI_ref_cal$Final[intersect(which(DEI_ref_cal$N_down>=dataset_support_count),which(DEI_ref_cal$N_up==0))]<-"down-regulate"
    DEI_ref_cal$isoform<-sapply(strsplit(rownames(DEI_ref_cal)," "),function(x){x[1]})
    DEI_ref_cal$compare<-sapply(strsplit(rownames(DEI_ref_cal)," "),function(x){x[2]})
    
    DEI_ref_cal$Var1<-rownames(DEI_ref_cal)
    
    ######ref_FC
    refFC_final <- DEIs_p_cal_f_fd
    DEI_ref_cal_f<-DEI_ref_cal[DEI_ref_cal$Var1 %in% refFC_final$isoform_compare, ]
    
    refFC_final$N_up<-DEI_ref_cal$N_up[match(refFC_final$isoform_compare,DEI_ref_cal$Var1)]
    refFC_final$N_non<-DEI_ref_cal$N_non[match(refFC_final$isoform_compare,DEI_ref_cal$Var1)]
    refFC_final$N_down<-DEI_ref_cal$N_down[match(refFC_final$isoform_compare,DEI_ref_cal$Var1)]
    refFC_final$Final<-DEI_ref_cal$Final[match(refFC_final$isoform_compare,DEI_ref_cal$Var1)]
    
    print(table(refFC_final$Final,refFC_final$compare))
    
    ## filter
    DEI_ref<-refFC_final[grep("regulate",refFC_final$Final),]
    print("High-confidence differentially expressed isoforms:")
    print(table(DEI_ref$Final,DEI_ref$compare))
    table(DEI_ref$Final,DEI_ref$compare)/nrow(expr.mat.tx) * 100
    # Ensure column names are lowercase with underscores
    # colnames(refFC_final) <- tolower(gsub(".", "_", colnames(refFC_final), fixed = TRUE))
    write.csv(refFC_final, paste0(output_prefix, "_RefData_DEIs.csv"))
  }
  
  message("Processing completed for tool: ", tool_name)
}

# Specify inputs for different tools
tools <- list(
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/isoform_level_count_kallisto_G_id_b22_s264.rds",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/kallisto",
    dataset_support_count = 7,
    tool_name = "kallisto",
    dei_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_kallisto_id_isoform.txt"
  ),
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/isoform_level_count_salmon_G_b22_s264.rds",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/salmon",
    dataset_support_count = 7,
    tool_name = "salmon",
    dei_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_Salmon_isoform.txt"
  ),
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/isoform_level_count_rsem_G_b22_s264.rds",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/rsem",
    dataset_support_count = 7,
    tool_name = "rsem",
    dei_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_RSEM_isoform.txt"
  ),
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/isoform_level_count_stringtie2_G_b22_s264.rds",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/stringtie2",
    dataset_support_count = 7,
    tool_name = "stringtie2",
    dei_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_StringTie2_isoform.txt"
  )
  # Add more tools here if necessary
)

# tools <- list(
#   list(
#     expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/isoform_level_count_stringtie2_G_b21_s252.rds",
#     batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
#                                "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
#                                "R_ILM_L6_B1", "P_BGI_L6_B1"),
#     output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/stringtie2",
#     dataset_support_count = 4,
#     tool_name = "stringtie2",
#     dei_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/tximport/srrnaseq/de/isoform/DE_results_StringTie2_isoform.txt"
#   )
# )

# Loop through tools and process datasets
for (tool in tools) {
  generate_reference_dataset(
    expr_file = tool$expr_file,
    batch_for_construction = tool$batch_for_construction,
    output_prefix = tool$output_prefix,
    dataset_support_count = tool$dataset_support_count,
    tool_name = tool$tool_name,
    dei_file = tool$dei_file
  )
}

# remove all
rm(list=ls())
gc()

