#!/bin/bash
#SBATCH --job-name=kallisto_index            # Job name
#SBATCH --output=logs/kallisto_index_%j.log  # Standard output and error log
#SBATCH --ntasks=1                          # Single task job
#SBATCH --cpus-per-task=8                   # Number of CPU cores
#SBATCH --mem=32G                           # Memory allocation
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au      # Where to send email notifications

# Enable error detection and pipe failure detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define paths
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast

# Transcriptome file
transcriptome_fa=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.transcripts.fa.gz

# Kallisto index output directory
workdir1=/vast/projects/quartet_rna_refdata/ref_transcriptome/kallisto
mkdir -p ${workdir1}

# Log start time
date
echo "kallisto_index begins"

start_time=$(date +%s)

# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Run Kallisto indexing with Singularity
singularity exec --bind $BIND_PATHS ${imgDIR}/rna_quantification_0.3.sif kallisto index \
    -k 31 \
    -i ${workdir1}/kallisto_trans \
    ${transcriptome_fa}

# Log end time and execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "kallisto_index ends"

echo "Execution time: ${execution_time} seconds"
