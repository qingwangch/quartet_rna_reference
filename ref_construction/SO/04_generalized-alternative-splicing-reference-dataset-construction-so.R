#' generalized-alternative-splicing-reference-dataset-construction-so
#' Qingwang Chen
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


# ===== Library Import =====
# library(tidyverse)
library(reshape2)
library(data.table)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(future.apply)
library(progressr)

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


# ===== Function: Generate Reference Dataset for Alternative Splicing =====
generate_as_reference_dataset <- function(expr_file, batch_for_construction, output_prefix, tool_name, das_file = NULL) {
  
  # ===== Load Expression Matrix (PSI values) =====
  expr.mat.tx <- fread(expr_file) %>% as.data.frame()
  rownames(expr.mat.tx) <- expr.mat.tx$V1
  expr.mat.tx$V1 <- NULL
  
  # ===== Load Metadata =====
  metadata <- process_metadata(expr.mat.tx)
  metadata$code <- metadata$sample_id
  metadata$Batch <- metadata$batch_id
  metadata$group <- ifelse(grepl("HCC", metadata$sample_id),"HCC1395","Quartet")
  # batch_for_construction <- c("D_ONT_LG_B1", "P_ONT_LG_B1", "P_ONT_LN_B2", "M_PAB_LN_B1", "M_PAB_LG_B2")
  metadata <- metadata[metadata$batch_id %in% batch_for_construction,]
  metadata <- metadata[metadata$group=="Quartet", ]
  
  # Filter expression data
  expr.mat.tx <- expr.mat.tx[, metadata$sample_id]
  
  # Process metadata
  metadata$batch_sample <- paste(metadata$batch_id, metadata$sample, sep = "_")
  meta_1 <- metadata
  
  # ===== Filter Expression Data =====
  expr.mat.tx <- expr.mat.tx[, metadata$sample_id]
  expr.mat.tx[is.na(expr.mat.tx)] <- 0
  
  # Remove AS events that are always 0 or 1
  filtered_expr.mat.tx <- expr.mat.tx[!apply(expr.mat.tx, 1, function(row) all(row %in% c(0, 1))), ]
  
  # ===== Prepare Long-Format Data =====
  exprMat_psi_M <- reshape2::melt(as.matrix(filtered_expr.mat.tx))
  colnames(exprMat_psi_M)[1:3] <- c("ASE", "library", "psi")
  exprMat_psi_M$batch_sample <- metadata$batch_sample[match(exprMat_psi_M$library, metadata$sample_id)]
  
  # ===== Filter PSI Values (0.05 ≤ PSI ≤ 0.95) =====
  expd_M <- exprMat_psi_M %>%
    group_by(ASE) %>%
    filter(psi >= 0.05 & psi <= 0.95) %>%
    ungroup()
  expd_M <- expd_M[expd_M$library %in% metadata$sample_id,]
  
  # ===== Count AS Events Per Sample =====
  expd_pg <- data.frame(tapply(expd_M$ASE, as.factor(paste(expd_M$ASE, expd_M$batch_sample)), length))
  colnames(expd_pg) <- "Num"
  expd_pg$ASE <- sapply(strsplit(as.character(rownames(expd_pg)), " "), function(x) x[1])
  expd_pg$batch_sample <- sapply(strsplit(as.character(rownames(expd_pg)), " "), function(x) x[2])
  expd_pg$batch <- metadata$batch_id[match(expd_pg$batch_sample, metadata$batch_sample)]
  expd_pg$sample <- metadata$sample[match(expd_pg$batch_sample, metadata$batch_sample)]
  
  # ===== Filter ASEs Detected in Multiple Samples =====
  expd_pg_f <- expd_pg[expd_pg$Num > 1, ]
  print(table(expd_pg_f$batch,expd_pg_f$sample))
  col_means <- apply(table(expd_pg_f$batch,expd_pg_f$sample), 2, mean);print(col_means)
  write.csv(expd_pg_f, paste0(output_prefix, "_detect_aselist_acrossSample.csv"))
  
  # ===== Consensus ASEs Across Batches =====
  # Ensure metadata columns are character type
  meta <- metadata
  for (i in 1:ncol(meta)) {
    meta[,i] <- as.character(meta[,i])
  }
  
  # Read the list of detected genes
  expd_ase1_p <- expd_pg_f
  
  # Sample combinations
  sample_combn <- data.frame(
    sampleA=rep("D6", 3),
    sampleB=c("D5", "F7", "M8")
  )
  sample_combn$compare <- paste(sample_combn$sampleB, "/", sample_combn$sampleA, sep="")
  
  # Process samples
  usample <- unique(as.character(meta$sample))
  
  expd_refs <- c()
  for (i in 1:length(usample)) {
    x <- expd_ase1_p[expd_ase1_p$sample == usample[i],]
    expd_sample <- data.frame(table(x$ASE))
    expd_sample <- cbind(expd_sample, usample[i])
    expd_refs <- rbind(expd_refs, expd_sample)
  }
  
  expd_refs <- data.frame(expd_refs)
  colnames(expd_refs) <- c("ASE", "freq", "sample")
  expd_refs$ASE <- as.character(expd_refs$ASE)
  expd_refs$sample <- as.character(expd_refs$sample)
  
  # Filter ASEs detected in all batches
  expd_refs_f <- expd_refs[expd_refs$freq == length(unique(metadata$batch_id)), ]
  write.csv(expd_refs_f, paste0(output_prefix, "_detect_aselist_acrossBatch.csv"))
  print(table(expd_refs_f$sample))
  
  # ===== Differential AS Analysis (DAS) =====
  if (!is.null(das_file) && file.exists(das_file)) {
    
    DAS_p <- fread(das_file) %>% as.data.frame()
    # ===== Filter DAS Data =====
    DAS_p <- DAS_p[DAS_p$batch %in% batch_for_construction,]
    
    sample_combn <- data.frame(sampleA = c("D5", "F7", "M8"), sampleB = "D6")
    sample_combn$compare <- paste(sample_combn$sampleA, "/", sample_combn$sampleB, sep = "")
    sample_combn[] <- lapply(sample_combn, as.character)
    
    # define sample pair
    detect_das_pairs<-c()
    for (i in 1:nrow(sample_combn)){
      s<-intersect(expd_refs_f$ASE[expd_refs_f$sample==sample_combn$sampleA[i]],expd_refs_f$ASE[expd_refs_f$sample==sample_combn$sampleB[i]])
      detect_das_pairs<-rbind(detect_das_pairs,cbind(s,as.character(sample_combn$compare[i])))
    }
    colnames(detect_das_pairs)<-c("das","compare")
    detect_das_pairs<-data.frame(detect_das_pairs)
    
    detect_das_pairs$gene_compare<-paste(detect_das_pairs$das,detect_das_pairs$compare)
    
    print(table(detect_das_pairs$compare))
    
    # delta_PSI and P
    DAS_p$gene_compare <- paste(DAS_p$ase,DAS_p$compare)
    
    DAS_p_dec <- DAS_p[DAS_p$gene_compare %in% detect_das_pairs$gene_compare,]
    print(table(DAS_p_dec$batch,DAS_p_dec$compare))
    
    # Filter DAS events with ΔPSI >= 0.05
    # DAS_p_f <- DAS_p_dec[DAS_p_dec$p_value<0.05,]
    DAS_p_f <- DAS_p_dec[abs(DAS_p_dec$delta_psi)>=0.05,]
    # DEGs_p_f <- DEGs_p_dec[abs(DEGs_p_dec$logfc)>1,]
    print(table(DAS_p_f$batch,DAS_p_f$compare))
    
    DAS_p_cal <- data.frame(table(paste(DAS_p_f$ase,DAS_p_f$compare)))
    DAS_p_cal$ase <- sapply(strsplit(as.character(DAS_p_cal$Var1)," "),function(x){x[1]})
    DAS_p_cal$compare <- sapply(strsplit(as.character(DAS_p_cal$Var1)," "),function(x){x[2]})
    
    DAS_p_cal_f <- DAS_p_cal[DAS_p_cal$Freq>=2,]
    DAS_p_cal_f_fd <- DAS_p_cal_f[DAS_p_cal_f$Var1 %in% detect_das_pairs$gene_compare,]
    print(table(DAS_p_cal_f_fd$compare))
    
    DAS_p_cal_f_fd <- DAS_p_cal_f_fd[order(DAS_p_cal_f_fd$compare),]
    DAS_p_cal_f_fd$gene_compare <- paste(DAS_p_cal_f_fd$ase,DAS_p_cal_f_fd$compare)
    
    ####meanFC medianp,seFC
    DAS_p$gene_compare <- paste(DAS_p$ase,DAS_p$compare)
    
    #function #
    gene_compare=levels(as.factor(DAS_p$gene_compare))
    
    mm.DAS<-data.frame(
      mean_delta_psi=tapply(DAS_p$delta_psi,as.factor(DAS_p$gene_compare),mean),
      medianp=tapply(DAS_p$p_value,as.factor(DAS_p$gene_compare),median),
      gene_compare=levels(as.factor(DAS_p$gene_compare))
    )
    
    DAS_p_cal_f_fd$mean_delta_psi <- mm.DAS$mean_delta_psi[match(DAS_p_cal_f_fd$gene_compare,mm.DAS$gene_compare)]
    DAS_p_cal_f_fd$medianp <- mm.DAS$medianp[match(DAS_p_cal_f_fd$gene_compare,mm.DAS$gene_compare)]
    
    # Output results to CSV file
    write.csv(DAS_p_cal_f_fd, paste0(output_prefix, "_ref_expr_das.csv"))
    
    # ===== Reference DAS (Consensus Across Batches) =====
    # preprocess
    DAS_p$protocol<-sapply(strsplit(as.character(DAS_p$batch),"_"),function(x){x[1]})
    DAS_p$gene_compare<-paste(DAS_p$ase,DAS_p$compare)
    
    DAS_p$DAS_type<-"non-DAS"
    DAS_p$DAS_type[intersect(which(DAS_p$p_value<0.05),which(DAS_p$delta_psi>=0.05))]<-"up-regulate"
    DAS_p$DAS_type[intersect(which(DAS_p$p_value<0.05),which(DAS_p$delta_psi<=(-0.05)))]<-"down-regulate"
    
    ########DEGs_cal 
    DAS_ref_cal<-data.frame(
      N_up=tapply(DAS_p$DAS_type,as.factor(DAS_p$gene_compare),function(x){length(which(x=="up-regulate"))}),
      N_non=tapply(DAS_p$DAS_type,as.factor(DAS_p$gene_compare),function(x){length(which(x=="non-DAS"))}),
      N_down=tapply(DAS_p$DAS_type,as.factor(DAS_p$gene_compare),function(x){length(which(x=="down-regulate"))}))
    
    DAS_ref_cal$Final<-"non-DAS"
    DAS_ref_cal$Final[intersect(which(DAS_ref_cal$N_up>=1),which(DAS_ref_cal$N_down>=1))]<-"conflicting"
    
    DAS_ref_cal$Final[intersect(which(DAS_ref_cal$N_up>=ceiling(length(batch_for_construction)/2)),which(DAS_ref_cal$N_down==0))]<-"up-regulate"
    DAS_ref_cal$Final[intersect(which(DAS_ref_cal$N_down>=ceiling(length(batch_for_construction)/2)),which(DAS_ref_cal$N_up==0))]<-"down-regulate"
    DAS_ref_cal$ase<-sapply(strsplit(rownames(DAS_ref_cal)," "),function(x){x[1]})
    DAS_ref_cal$compare<-sapply(strsplit(rownames(DAS_ref_cal)," "),function(x){x[2]})
    
    DAS_ref_cal$Var1<-rownames(DAS_ref_cal)
    
    ######ref_FC
    refFC_final <- DAS_p_cal_f_fd
    DAS_ref_cal_f<-DAS_ref_cal[DAS_ref_cal$Var1 %in% refFC_final$gene_compare, ]
    
    refFC_final$N_up<-DAS_ref_cal$N_up[match(refFC_final$gene_compare,DAS_ref_cal$Var1)]
    refFC_final$N_non<-DAS_ref_cal$N_non[match(refFC_final$gene_compare,DAS_ref_cal$Var1)]
    refFC_final$N_down<-DAS_ref_cal$N_down[match(refFC_final$gene_compare,DAS_ref_cal$Var1)]
    refFC_final$Final<-DAS_ref_cal$Final[match(refFC_final$gene_compare,DAS_ref_cal$Var1)]
    
    print(table(refFC_final$Final,refFC_final$compare))
    
    ## filter
    DAS_ref<-refFC_final[grep("regulate",refFC_final$Final),]
    print(table(DAS_ref$Final,DAS_ref$compare))
    
    # Ensure column names are lowercase with underscores
    colnames(refFC_final) <- tolower(gsub(".", "_", colnames(refFC_final), fixed = TRUE))
    
    write.csv(refFC_final, paste0(output_prefix, "_RefData_DAS.csv"))
  }
  
  print(paste("Processing completed for tool:", tool_name))
}

# ===== Define Tools for Processing =====
tools <- list(
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/kallisto/quartet_as_events_psi_combined_s264_so_e_g.psi",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "./ref_data_construction/kallisto_s156",
    tool_name = "kallisto",
    das_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/kallisto_seg_DAS_s264.txt"
  ),
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/rsem/quartet_as_events_psi_combined_s264_so_e_g.psi",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "./ref_data_construction/rsem_s156",
    tool_name = "rsem",
    das_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/rsem_seg_DAS_s264.txt"
  ),
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/salmon/quartet_as_events_psi_combined_s264_so_e_g.psi",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "./ref_data_construction/salmon_s156",
    tool_name = "salmon",
    das_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/salmon_seg_DAS_s264.txt"
  ),
  list(
    expr_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/stringtie2/SO-E/quartet_as_events_psi_combined_s264_so_e_g.psi",
    batch_for_construction = c("P_ILM_L8_B1", "R_ILM_L4_B2", "R_BGI_L3_B1", "R_ILM_L4_B3", "R_ILM_L2_B2", "R_ILM_L8_B1",
                               "R_ILM_L1_B1", "P_ILM_L1_B1", "R_ILM_L5_B1", "P_BGI_L3_B1", "R_BGI_L6_B1",
                               "R_ILM_L6_B1", "P_BGI_L6_B1"),
    output_prefix = "./ref_data_construction/stso_s156",
    tool_name = "stringtie2",
    das_file = "/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/suppa2/diffsplice/stringtie2_seg_DAS_s264.txt"
  )
)

# ===== Loop Through Tools and Generate Reference Datasets =====
for (tool in tools) {
  generate_as_reference_dataset(
    expr_file = tool$expr_file,
    batch_for_construction = tool$batch_for_construction,
    output_prefix = tool$output_prefix,
    tool_name = tool$tool_name,
    das_file = tool$das_file
  )
}

# ===== Clean Up =====
rm(list = ls())
gc()

