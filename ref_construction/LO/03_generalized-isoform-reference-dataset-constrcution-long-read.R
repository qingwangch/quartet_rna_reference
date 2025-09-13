#' generalized-isoform-reference-dataset-constrcution-long-read
#' Qingwang Chen
#' 2025-05-16

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
# library(tidyverse)
library(magrittr)
library(purrr)
library(dplyr)
library(data.table)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(reshape2)

# Function to generate reference datasets
generate_reference_dataset <- function(expr_file, metadata_file, batch_for_construction,output_prefix, dataset_support_count, tool_name, dei_file = NULL) {
  # Load expression matrix
  expr.mat.tx <- readRDS(expr_file)
  expr.mat.tx[is.na(expr.mat.tx)] <- 0
  
  # ===== Load Metadata =====
  metadata <- fread(metadata_file) %>% as.data.frame()
  metadata$batch_sample <- paste(metadata$batch_id, metadata$sample, sep = "_")
  metadata <- metadata[metadata$batch_id %in% batch_for_construction, ]
  metadata <- metadata[metadata$group=="Quartet", ]
  
  # ===== Filter Expression Data =====
  expr.mat.tx <- expr.mat.tx[, metadata$sample_id]
  expr.mat.tx[is.na(expr.mat.tx)] <- 0
  
  # Process metadata
  metadata$batch_sample <- paste(metadata$batch_id, metadata$sample, sep = "_")
  meta_1 <- metadata
  ###################----------Char_detect-----------##############
  # Melt expression matrix
  exprMat_count_M <-reshape2::melt(as.matrix(expr.mat.tx))
  colnames(exprMat_count_M)[1:3] <- c("isoform", "library", "count")
  
  # Merge metadata
  exprMat_count_M$batch_sample <- meta_1$batch_sample[match(exprMat_count_M$library, meta_1$sample_id)]
  
  # Filter genes with counts >= 1
  expd_M <- exprMat_count_M[exprMat_count_M$count >= 1,]
  
  # Calculate detected genes per sample
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
    detect_isoforms_pairs<-c()
    for (i in 1:nrow(sample_combn)){
      s<-intersect(expd_refs_f$isoform[expd_refs_f$sample==sample_combn$sampleA[i]],expd_refs_f$isoform[expd_refs_f$sample==sample_combn$sampleB[i]])
      detect_isoforms_pairs<-rbind(detect_isoforms_pairs,cbind(as.character(s),as.character(sample_combn$compare[i])))
    }
    colnames(detect_isoforms_pairs)<-c("isoform","compare")
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
    write.csv(DEIs_p_cal_f_fd, paste0(output_prefix, "_ref_expr_s84.csv"))
    
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
    
    DEI_ref_cal$Final[intersect(which(DEI_ref_cal$N_up>=ceiling(length(batch_for_construction)/2)),which(DEI_ref_cal$N_down==0))]<-"up-regulate"
    DEI_ref_cal$Final[intersect(which(DEI_ref_cal$N_down>=ceiling(length(batch_for_construction)/2)),which(DEI_ref_cal$N_up==0))]<-"down-regulate"
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
    write.csv(refFC_final, paste0(output_prefix, "_RefData_DEIs_s84.csv"))
  }
  
  message("Processing completed for tool: ", tool_name)
}

# Specify inputs for different tools
tools <- list(
  list(
    expr_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-minimap2-stringtie2-g-txq-count-s84t274031.rds",
    metadata_file = "/vast/projects/quartet_rna_refdata/analysis/suppa2/metadata_s192.csv",
    batch_for_construction = c("P_ONT_LG_B1", "P_ONT_LG_B2", "P_ONT_LN_B2", "M_PAB_LN_B1", "D_ONT_LG_B1", "D_ONT_LW_B1", "M_PAB_LG_B2"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/stringtie2",
    dataset_support_count = 3,
    tool_name = "StringTie2",
    dei_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/D7/DE_results_StringTie2_isoform.txt"
  ),
  list(
    expr_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-minimap2-bambu-g-txq-count-s84t274031.rds",
    metadata_file = "/vast/projects/quartet_rna_refdata/analysis/suppa2/metadata_s192.csv",
    batch_for_construction = c("P_ONT_LG_B1", "P_ONT_LG_B2", "P_ONT_LN_B2", "M_PAB_LN_B1", "D_ONT_LG_B1", "D_ONT_LW_B1", "M_PAB_LG_B2"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/bambu",
    dataset_support_count = 3,
    tool_name = "Bambu",
    dei_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/D7/DE_results_Bambu_isoform.txt"
  ),
  list(
    expr_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-minimap2-isoquant-g-txq-count-s84t252913.rds",
    metadata_file = "/vast/projects/quartet_rna_refdata/analysis/suppa2/metadata_s192.csv",
    batch_for_construction = c("P_ONT_LG_B1", "P_ONT_LG_B2", "P_ONT_LN_B2", "M_PAB_LN_B1", "D_ONT_LG_B1", "D_ONT_LW_B1", "M_PAB_LG_B2"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/isoquant",
    dataset_support_count = 3,
    tool_name = "IsoQuant",
    dei_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/D7/DE_results_IsoQuant_isoform.txt"
  ),
  list(
    expr_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-minimap2-oarfish-g-txq-count-s84t252913.rds",
    metadata_file = "/vast/projects/quartet_rna_refdata/analysis/suppa2/metadata_s192.csv",
    batch_for_construction = c("P_ONT_LG_B1", "P_ONT_LG_B2", "P_ONT_LN_B2", "M_PAB_LN_B1", "D_ONT_LG_B1", "D_ONT_LW_B1", "M_PAB_LG_B2"),
    output_prefix = "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/D7/oarfish",
    dataset_support_count = 3,
    tool_name = "Oarfish",
    dei_file = "/vast/projects/quartet_rna_refdata/analysis/isoforms/de/isoform/D7/DE_results_Oarfish_isoform.txt"
  )
  # Add more tools here if necessary
)

# Loop through tools and process datasets
for (tool in tools) {
  generate_reference_dataset(
    expr_file = tool$expr_file,
    metadata_file = tool$metadata_file,
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
