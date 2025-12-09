#!/bin/bash
#SBATCH --job-name=bam_to_sam_conversion         # Job name
#SBATCH --output=logs/bam_to_sam_%j.log          # Standard output and error log
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --cpus-per-task=8                        # Number of CPU cores
#SBATCH --mem=16G                                # Memory allocation
#SBATCH --time=12:00:00                          # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au           # Where to send email notifications

# Enable error detection and pipe failure detection
set -eo pipefail

# Load required modules
module load singularity
module load samtools/1.21

# Ensure logs directory exists
mkdir -p logs

# Set working directory
cd /vast/projects/quartet_rna_refdata/

# Log start time
date
echo "BAM to SAM conversion begins"
start_time=$(date +%s)

# Parse named parameters
BAM_FILE=""
OUTPUT_DIR=""
SAMPLE_ID=""

for arg in "$@"; do
    case $arg in
        --bam=*) BAM_FILE="${arg#*=}"; shift ;;
        --output-dir=*) OUTPUT_DIR="${arg#*=}"; shift ;;
        --sample-id=*) SAMPLE_ID="${arg#*=}"; shift ;;
        --help)
            echo "Usage: sbatch bam_to_sam.slurm --bam=<bam_file> --output-dir=<output_dir> --sample-id=<sample_id>"
            exit 1 ;;
        *) echo "Unknown parameter: $arg"; exit 1 ;;
    esac
done

# Validate required parameters
if [[ -z "$BAM_FILE" || -z "$OUTPUT_DIR" || -z "$SAMPLE_ID" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: sbatch bam_to_sam.slurm --bam=<bam_file> --output-dir=<output_dir> --sample-id=<sample_id>"
    exit 1
fi

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Define output SAM file
OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_ID}.sam"

# Print inputs and outputs
echo "Input BAM file: ${BAM_FILE}"
echo "Output SAM file: ${OUTPUT_SAM}"
echo "Sample ID: ${SAMPLE_ID}"

# Run samtools to convert BAM to SAM
samtools view -h -o ${OUTPUT_SAM} ${BAM_FILE}

# Check if conversion was successful
if [[ -f "${OUTPUT_SAM}" ]]; then
    echo "Conversion completed successfully: ${OUTPUT_SAM}"
else
    echo "Error: Conversion failed. Output file was not created."
    exit 1
fi

# Log end time and execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "BAM to SAM conversion ends"
echo "Execution time: ${execution_time} seconds"
