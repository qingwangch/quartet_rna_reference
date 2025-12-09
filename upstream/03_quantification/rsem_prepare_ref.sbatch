#!/bin/bash
#SBATCH --job-name=rsem_prepare_ref          # Job name
#SBATCH --output=logs/rsem_prepare_ref_%j.log # Standard output and error log
#SBATCH --ntasks=1                           # Single task job
#SBATCH --cpus-per-task=8                    # Number of CPU cores
#SBATCH --mem=128G                            # Memory allocation
#SBATCH --time=48:00:00                      # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au       # Where to send email notifications

# Enable error detection and pipe failure detection
set -eo pipefail

# Load conda environment
source ~/miniconda3/bin/activate /home/users/allstaff/chen.q/miniconda3/envs/rsem-env

# Define paths
genome_ref=/vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
annotation_gtf=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf
workdir=/vast
workdir1=/vast/projects/quartet_rna_refdata/ref_genome/rsem

# Create necessary directories
mkdir -p ${workdir1}
mkdir -p ${workdir1}/ref/

# Set working directory
cd ${workdir}

# Log start time
date
echo "rsem_prepare_reference begins"
start_time=$(date +%s)

# Run rsem-prepare-reference with STAR 2.7.10b
rsem-prepare-reference \
    --gtf ${annotation_gtf} \
    --star \
    --star-path ~/miniconda3/envs/rsem-env/bin/ \
    ${genome_ref} \
    ${workdir1}/ref/human_ensembl

# Log end time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "rsem_prepare_reference ends"
echo "Execution time: $execution_time seconds"
