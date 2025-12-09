#!/bin/bash
#SBATCH --job-name=bambu_multi_sample  # Job name
#SBATCH --partition=long               # Use the long partition
#SBATCH --output=logs/bambu_output_%j.log # Output log
#SBATCH --time=96:00:00               # Maximum wall time (2 weeks)
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=64             # Number of CPUs per task
#SBATCH --mem=252G                     # Maximum memory allocation

# Load Conda
source /home/users/allstaff/chen.q/miniconda3/bin/activate /home/users/allstaff/chen.q/miniconda3/envs/bambu_env

# Define paths
BAM_LIST=/vast/projects/quartet_rna_refdata/slurm/bambu/quartet_lrbam_list_s12_lw.csv
OUTPUT_DIR=/vast/projects/quartet_rna_refdata/analysis/bambu/quartet_lr_s12_lw
GTF_FILE=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf
FASTA_FILE=/vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
THREADS=64  # Use maximum CPUs allocated

# Track start time
START_TIME=$(date +%s)
echo "Job started at: $(date)"

# Ensure BiocFileCache is installed
Rscript -e "if (!requireNamespace('BiocFileCache', quietly = TRUE)) BiocManager::install('BiocFileCache')"

# Run the R analysis
Rscript /vast/projects/quartet_rna_refdata/slurm/bambu/bambu_multiSample_LR.R \
  --bam_files $BAM_LIST \
  --output_dir $OUTPUT_DIR \
  --gtf $GTF_FILE \
  --fasta $FASTA_FILE 

# Track end time and print runtime
END_TIME=$(date +%s)
echo "Job finished at: $(date)"
RUNTIME=$((END_TIME - START_TIME))
echo "Total runtime: $RUNTIME seconds"
