#!/bin/bash
#SBATCH --job-name=rsem_quant                # Job name
#SBATCH --output=logs/rsem_quant_%j.log     # Standard output and error log
#SBATCH --ntasks=1                          # Single task job
#SBATCH --cpus-per-task=16                  # Number of CPU cores
#SBATCH --mem=64G                           # Memory allocation
#SBATCH --time=48:00:00                     # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au      # Where to send email notifications

# Enable error detection and pipe failure detection
set -eo pipefail

# Load conda environment
source ~/miniconda3/bin/activate /home/users/allstaff/chen.q/miniconda3/envs/rsem-env

# Define paths
rsem_index=/vast/projects/quartet_rna_refdata/ref_genome/rsem/ref/human_ensembl
workdir=/vast
workdir1=/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/rsem

# Function to display usage
usage() {
    echo "Usage: sbatch rsem_quant_slurm.sh --fq_r1 <R1.fastq.gz> --fq_r2 <R2.fastq.gz> --sample_id <sample_name> --fwd_prob <probability>"
    echo "Example: sbatch rsem_quant_slurm.sh --fq_r1 sample_R1.fastq.gz --fq_r2 sample_R2.fastq.gz --sample_id sample_name --fwd_prob 0.5"
    exit 1
}

# Parse command-line arguments
fq_r1=""
fq_r2=""
sample_id=""
fwd_prob=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fq_r1) fq_r1="$2"; shift ;;
        --fq_r2) fq_r2="$2"; shift ;;
        --sample_id) sample_id="$2"; shift ;;
        --fwd_prob) fwd_prob="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Validate required input arguments
if [[ -z "$fq_r1" || -z "$fq_r2" || -z "$sample_id" || -z "$fwd_prob" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Ensure workdir2 exists
workdir2="${workdir1}/${sample_id}"
mkdir -p ${workdir2}

# Log start time
date
echo "rsem_quant begins for sample: ${sample_id}"
echo "Input R1: ${fq_r1}"
echo "Input R2: ${fq_r2}"
echo "Output Directory: ${workdir2}"
echo "Forward Probability: ${fwd_prob}"
start_time=$(date +%s)

# Run RSEM quantification
rsem-calculate-expression --paired-end \
    --star \
    --star-path ~/miniconda3/envs/rsem-env/bin/ \
    --star-gzipped-read-file \
    -p 16 \
    --forward-prob ${fwd_prob} \
    ${fq_r1} \
    ${fq_r2} \
    ${rsem_index} \
    ${workdir2}/${sample_id}_rsem

# Log end time and execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "rsem_quant ends for sample: ${sample_id}"
echo "Execution time: ${execution_time} seconds"
