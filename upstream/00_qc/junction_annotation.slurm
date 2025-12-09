#!/bin/bash
#SBATCH --job-name=junction_annotation      # Job name
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --cpus-per-task=4                   # Number of CPU cores
#SBATCH --mem=64G                           # Memory allocation
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec
#SBATCH --output=logs/junction_annotation_%j.log # Standard output and error log
#SBATCH --mail-type=FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au      # Where to send email notifications

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
        --bam_file=*)
        bam_file="${arg#*=}"
        shift
        ;;
        *)
        echo "Invalid argument: $arg"
        exit 1
        ;;
    esac
done

# Check if both parameters are provided
if [ -z "$sample_id" ] || [ -z "$bam_file" ]; then
    echo "Error: Missing required parameters."
    echo "Usage: sbatch junction_annotation.slurm --sample_id=<sample_id> --bam_file=<bam_file>"
    exit 1
fi

# Change to the working directory
cd ${workdir}

date
echo "junction_annotation begins"

start_time=$(date +%s)

# Output directory
output_dir=/vast/projects/quartet_rna_refdata/analysis/rseqc/${sample_id}

# Check if the output directory exists, create it if not
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
    echo "Directory created: $output_dir"
else
    echo "Directory already exists: $output_dir"
fi

# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Run junction_annotation.py using Singularity
singularity exec --bind $BIND_PATHS $imgDIR/rseqc_v1.sif /usr/local/bin/junction_annotation.py -i ${bam_file} -r ${bed} -o ${output_dir}/${sample_id}.junction_annotation_

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "junction_annotation ends"

echo "Execution time: $execution_time seconds"
