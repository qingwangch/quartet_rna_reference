#!/bin/bash
#SBATCH --job-name=geneBody_coverage         # Job name
#SBATCH --ntasks=1                           # Run on a single CPU
#SBATCH --cpus-per-task=4                    # Number of CPU cores
#SBATCH --mem=16G                            # Memory allocation
#SBATCH --time=48:00:00                      # Time limit hrs:min:sec
#SBATCH --output=logs/geneBody_coverage_%j.log    # Standard output and error log

# Enable error detection and pipe failure detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define paths
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast/projects/quartet_rna_refdata/
bed=/vast/projects/quartet_rna_refdata/ref_transcriptome/hg38_GENCODE_V43_Comprehensive.bed

# Parse named parameters
for arg in "$@"
do
    case $arg in
        --sample_id=*)
        sample_id="${arg#*=}"
        shift
        ;;
        --bam_list=*)
        bam_list="${arg#*=}"
        shift
        ;;
        *)
        echo "Invalid argument: $arg"
        exit 1
        ;;
    esac
done

# Check if both parameters are provided
if [ -z "$sample_id" ] || [ -z "$bam_list" ]; then
    echo "Error: Missing required parameters."
    echo "Usage: sbatch read_distribution.slurm --sample_id=<sample_id> --bam_list=<bam_list>"
    exit 1
fi

# Change to the working directory
cd ${workdir}

date
echo "read_distribution begins"

start_time=$(date +%s)

# Output directory
output=/vast/projects/quartet_rna_refdata/analysis/rseqc/${sample_id}

# Check if the output directory exists, create it if not
if [ ! -d "$output" ]; then
    mkdir -p "$output"
    echo "Directory created: $output"
else
    echo "Directory already exists: $output"
fi

# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Run junction_annotation.py using Singularity
singularity exec --bind $BIND_PATHS $imgDIR/rseqc_v1.sif /usr/local/bin/read_distribution.py -r ${bed} -i ${bam_list} > ${output}/${sample_id}_read_distribution.txt

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "read_distribution ends"

echo "Execution time: $execution_time seconds"
