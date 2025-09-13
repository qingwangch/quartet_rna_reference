#' isoform-data-preparation-for-reference-dataset-construction
#' Qingwang Chen
#' 2025-05-15

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# Specify a new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths so that the custom path has priority
.libPaths(c(custom_lib_path, .libPaths()))

# Check the current library paths
print(.libPaths())

# Load libraries
library(data.table)
library(magrittr)
library(purrr)
library(dplyr)

# Import data
## Isoform ExprData
### IsoQuant
expr.mat.tx.iq.d1 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_ont_d_s12/D_ONT_LG_B1/D_ONT_LG_B1.transcript_grouped_counts.tsv") %>% as.data.frame()
expr.mat.tx.iq.d2 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_ont_d_lw_s12/D_ONT_LW_B1/D_ONT_LW_B1.transcript_grouped_counts.tsv") %>% as.data.frame()

expr.mat.tx.i.p1 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_ont_s24/P_ONT_LG_B1/P_ONT_LG_B1.transcript_grouped_counts.tsv") %>% as.data.frame()
expr.mat.tx.i.p2 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_P_ONT_LG_B2/P_ONT_LG_B2/P_ONT_LG_B2.transcript_grouped_counts.tsv") %>% as.data.frame()
expr.mat.tx.i.p3 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_ont_s24/P_ONT_LN_B2/P_ONT_LN_B2.transcript_grouped_counts.tsv") %>% as.data.frame()
# expr.mat.tx.i.p3 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_ont_s24/P_ONT_LN_B2/P_ONT_LN_B2.transcript_grouped_counts.tsv") %>% as.data.frame()


expr.mat.tx.i.m1 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_pab_m_s12/M_PAB_LN_B1/M_PAB_LN_B1.transcript_grouped_counts.tsv") %>% as.data.frame()
expr.mat.tx.i.m2 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_pab_m_lg_s12/M_PAB_LG_B2/M_PAB_LG_B2.transcript_grouped_counts.tsv") %>% as.data.frame()
# expr.mat.tx.i.m2 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_pab_m_lg_s12/M_PAB_LG_B2/M_PAB_LG_B2.transcript_grouped_counts.tsv") %>% as.data.frame()
# expr.mat.tx.i.i1 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_I_PAB_LN_B1_s12/I_PAB_LN_B1/I_PAB_LN_B1.transcript_grouped_counts.tsv") %>% as.data.frame()
# expr.mat.tx.i.i2 <- fread("/vast/projects/quartet_rna_refdata/analysis/isoquant/quartet_I_PAB_LN_B2_s12/I_PAB_LN_B2/I_PAB_LN_B2.transcript_grouped_counts.tsv") %>% as.data.frame()

# Put all data frames into a list
df_list <- list(expr.mat.tx.iq.d1, expr.mat.tx.iq.d2, expr.mat.tx.i.p1, expr.mat.tx.i.p2, expr.mat.tx.i.p3, expr.mat.tx.i.m1, expr.mat.tx.i.m2)

# Merge all data frames by the "#feature_id" column using reduce + full_join
merged_df <- reduce(df_list, full_join, by = "#feature_id")

# Set row names and remove the feature ID column
rownames(merged_df) <- merged_df[,1]; merged_df[,1] <- NULL
dim(merged_df)
saveRDS(merged_df,"./isoforms/Rdata/quartet-LO-minimap2-isoquant-g-txq-count-s84t252913.rds")

### StringTie2
expr.mat.tx.st <- fread("/vast/projects/quartet_rna_refdata/analysis/stringtie2/LO-E/count/ref_data_s84/merged_tx_count_s84/merged_count_from_cov.txt") %>% as.data.frame()
rownames(expr.mat.tx.st) <- expr.mat.tx.st$transcript_id; expr.mat.tx.st$transcript_id <- NULL
dim(expr.mat.tx.st)
saveRDS(expr.mat.tx.st,"./isoforms/Rdata/quartet-LO-minimap2-stringtie2-g-txq-count-s84t274031.rds")

### StringTie3
# expr.mat.tx.st <- fread("/vast/projects/quartet_rna_refdata/analysis/stringtie3/LO-E/count/ref_data_s72/merged_tx_count_s72/merged_count_from_cov.txt") %>% as.data.frame()
# rownames(expr.mat.tx.st) <- expr.mat.tx.st$transcript_id; expr.mat.tx.st$transcript_id <- NULL
# dim(expr.mat.tx.st)
# saveRDS(expr.mat.tx.st,"./isoforms/Rdata/quartet-LO-minimap2-stringtie3-g-txq-count-s72t27407.rds")

### Bambu
expr.mat.tx.ba <- fread("/vast/projects/quartet_rna_refdata/analysis/bambu/quartet_lr_s72/Quartet_counts_transcript.txt") %>% as.data.frame()
expr.mat.tx.ba.lw <- fread("/vast/projects/quartet_rna_refdata/analysis/bambu/quartet_lr_s12_lw/Quartet_counts_transcript.txt") %>% as.data.frame()

# 1. Remove rows in TXNAME column that contain "bambu" (case-insensitive)
expr.mat.tx.ba <- expr.mat.tx.ba %>%
  filter(!grepl("bambu", TXNAME, ignore.case = TRUE))

expr.mat.tx.ba.lw <- expr.mat.tx.ba.lw %>%
  filter(!grepl("bambu", TXNAME, ignore.case = TRUE))

# 2. Set TXNAME column as row names
rownames(expr.mat.tx.ba) <- expr.mat.tx.ba$TXNAME
rownames(expr.mat.tx.ba.lw) <- expr.mat.tx.ba.lw$TXNAME

# 3. Remove GENEID and TXNAME columns (already used as row names)
expr.mat.tx.ba <- expr.mat.tx.ba %>%
  select(-GENEID, -TXNAME)
expr.mat.tx.ba.lw <- expr.mat.tx.ba.lw %>%
  select(-GENEID, -TXNAME)

# 4. Clean up sample names by removing ".minimap2.sorted.only_primary.bam" suffix
colnames(expr.mat.tx.ba) <- gsub("\\.minimap2\\.sorted\\.only_primary\\.bam$", "", colnames(expr.mat.tx.ba))
dim(expr.mat.tx.ba)
colnames(expr.mat.tx.ba.lw) <- gsub("\\.minimap2\\.sorted\\.only_primary\\.bam$", "", colnames(expr.mat.tx.ba.lw))
dim(expr.mat.tx.ba.lw)
expr.mat.tx.ba.m <- cbind(expr.mat.tx.ba,expr.mat.tx.ba.lw)
dim(expr.mat.tx.ba.m)
expr.mat.tx.ba.m <- expr.mat.tx.ba.m %>% 
  rename(D_ONT_LW_B1_D6_1 = D_ONT_LW_B1_D6_1_N)
saveRDS(expr.mat.tx.ba.m,"./isoforms/Rdata/quartet-LO-minimap2-bambu-g-txq-count-s84t274031.rds")

### Oarfish
count_oarfish <- readRDS("/vast/projects/quartet_rna_refdata/analysis/isoforms/Rdata/quartet-LO-minimap2-oarfish-g-txq-count-s144t252913.rds")
dim(count_oarfish)
metadata <- fread("/vast/projects/quartet_rna_refdata/analysis/suppa2/metadata_s192.csv") %>% as.data.frame()
batch_for_construction = c("P_ONT_LG_B1", "P_ONT_LG_B2", "P_ONT_LN_B2", "M_PAB_LN_B1", "D_ONT_LG_B1", "M_PAB_LG_B2","D_ONT_LW_B1")
metadata <- metadata[metadata$batch_id %in% batch_for_construction, ]
metadata <- metadata[metadata$group=="Quartet", ]

count_oarfish <- count_oarfish[, metadata$sample_id]
# Extract the first field before "|" from each row name and update the row names of count_oarfish
rownames(count_oarfish) <- sapply(strsplit(rownames(count_oarfish), "\\|"), `[`, 1)
head(rownames(count_oarfish))
dim(count_oarfish)
# Save both data frames into a single RDS file
saveRDS(count_oarfish, "./isoforms/Rdata/quartet-LO-minimap2-oarfish-g-txq-count-s84t252913.rds")

# Save gene count data

# Clear all variables
rm(list=ls())
gc()
