#!/bin/bash
#SBATCH --job-name=kallisto_quant             # Job name
#SBATCH --output=logs/kallisto_quant_%j.log   # Standard output and error log
#SBATCH --ntasks=1                           # Single task job
#SBATCH --cpus-per-task=16                  # Number of CPU cores
#SBATCH --mem=32G                           # Memory allocation
#SBATCH --time=24:00:00                      # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au       # Where to send email notifications

# Enable error detection and pipe failure detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define paths
imgDIR=/vast/projects/quartet_rna_refdata/images
PathToKallistoIndex=/vast/projects/quartet_rna_refdata/ref_transcriptome/kallisto/kallisto_trans

# Set working directory
workdir=/vast
cd ${workdir}

# Function to display usage
usage() {
    echo "Usage: $0 --fq_r1 PATH_TO_READ1 --fq_r2 PATH_TO_READ2 --sample_id SAMPLE_ID --strand-specific true|false"
    exit 1
}

# Parse input arguments
fq_r1=""
fq_r2=""
sample_id=""
strand_specific=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fq_r1) fq_r1="$2"; shift ;;
        --fq_r2) fq_r2="$2"; shift ;;
        --sample_id) sample_id="$2"; shift ;;
        --strand-specific) strand_specific="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Validate arguments
if [[ -z "$fq_r1" || -z "$fq_r2" || -z "$sample_id" || -z "$strand_specific" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if input files exist
if [[ ! -f "$fq_r1" || ! -f "$fq_r2" ]]; then
    echo "Error: Input FASTQ files not found: $fq_r1, $fq_r2"
    exit 1
fi

# Check if Kallisto index exists
if [[ ! -f "$PathToKallistoIndex" ]]; then
    echo "Error: Kallisto index not found: $PathToKallistoIndex"
    exit 1
fi

# Validate strand-specific parameter
if [[ "$strand_specific" != "true" && "$strand_specific" != "false" ]]; then
    echo "Error: Invalid value for --strand-specific. Use 'true' or 'false'."
    exit 1
fi

# Log start time
date
echo "kallisto_quant begins"
echo "Strand specificity: $strand_specific"

start_time=$(date +%s)

# Output directory for the current sample
workdir1=/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/kallisto/${sample_id}
mkdir -p "$workdir1"

# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Run Kallisto based on strand-specific value
if [[ "$strand_specific" == "true" ]]; then
    echo "Running Kallisto quantification with --rf-stranded..."
    singularity exec --bind $BIND_PATHS ${imgDIR}/rna_quantification_0.3.sif kallisto quant \
        -i ${PathToKallistoIndex} \
        -o ${workdir1} \
        --bootstrap-samples=50 \
        -t 16 \
        --rf-stranded ${fq_r1} ${fq_r2}
else
    echo "Running Kallisto quantification without strand-specificity..."
    singularity exec --bind $BIND_PATHS ${imgDIR}/rna_quantification_0.3.sif kallisto quant \
        -i ${PathToKallistoIndex} \
        -o ${workdir1} \
        --bootstrap-samples=50 \
        -t 16 \
        ${fq_r1} ${fq_r2}
fi

# Log end time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "kallisto_quant ends"
echo "Execution time: ${execution_time} seconds"
