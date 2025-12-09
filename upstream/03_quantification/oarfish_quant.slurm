#!/bin/bash
#SBATCH --job-name=oarfish_quant           # Job name
#SBATCH --output=logs/oarfish_quant_%j.log # Standard output and error log
#SBATCH --ntasks=1                         # Single task job
#SBATCH --cpus-per-task=8                  # Number of CPU cores
#SBATCH --mem=24G                          # Memory per node
#SBATCH --time=24:00:00                    # Time limit hrs:min:sec

# Enable error detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define paths
imgDIR=/vast/projects/quartet_rna_refdata/images/oarfish_0.6.2--hc9368f3_0.sif
outputDIR=/vast/projects/quartet_rna_refdata/analysis/oarfish_quant
reference=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.transcripts.fa  # Reference transcriptome

# Bind necessary directories for Singularity
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Parse input arguments
reads=""
sample_id=""
seq_tech="ont-cdna"  # Default to ONT cDNA sequencing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --reads) reads="$2"; shift ;;
        --sample-id) sample_id="$2"; shift ;;
        --seq-tech) seq_tech="$2"; shift ;;
        --help)
            echo "Usage: $0 --reads READS --sample-id SAMPLE_ID [--seq-tech SEQ_TECH]"
            exit 1 ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Validate inputs
if [[ -z "$reads" || -z "$sample_id" ]]; then
    echo "Error: --reads and --sample-id are required."
    echo "Usage: $0 --reads READS --sample-id SAMPLE_ID [--seq-tech SEQ_TECH]"
    exit 1
fi

# Ensure output directory exists
mkdir -p ${outputDIR}/${sample_id}
cd /vast/projects/quartet_rna_refdata/

# Log the start time
start_time=$(date +%s)
date
echo "Oarfish quantification begins for sample: $sample_id"

# Run Oarfish quantification with Singularity
singularity exec --bind $BIND_PATHS $imgDIR \
    oarfish --reads $reads \
            --reference $reference \
            --seq-tech $seq_tech \
            --model-coverage \
            --num-bootstraps 50\
            -o ${outputDIR}/${sample_id}/${sample_id} \
            -j 8

# Log the end time and execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "Oarfish quantification ends for sample: $sample_id"
echo "Execution time: $execution_time seconds"
