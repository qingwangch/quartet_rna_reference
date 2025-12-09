#!/bin/bash
#SBATCH --job-name=bam_flagstat          # Job name
#SBATCH --ntasks=1                       # Run on a single CPU
#SBATCH --cpus-per-task=4                # Number of CPU cores per task
#SBATCH --mem=16G                        # Job memory request
#SBATCH --time=01:00:00                  # Time limit hrs:min:sec
#SBATCH --output=logs/bam_flagstat_%j.log     # Standard output and error log

# Enable error detection
set -eo pipefail

# Load Singularity
module load singularity
imgDIR=/vast/projects/quartet_rna_refdata/images

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam) bam="$2"; shift ;;
        --outdir) outdir="$2"; shift ;;
        --sample_id) sample_id="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure required arguments are provided
if [ -z "$bam" ] || [ -z "$outdir" ] || [ -z "$sample_id" ]; then
    echo "Usage: $0 --bam <bam_file> --outdir <output_directory> --sample_id <sample_identifier>"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${outdir}"

# Define the working directory
workdir=/vast/projects/quartet_rna_refdata

# Change to the working directory
cd ${workdir}

date
echo "samtools flagstat begins"

# Run samtools flagstat using Singularity and save output with sample_id in specified directory
output_file="${outdir}/${sample_id}_flagstat.txt"
singularity exec $imgDIR/flair_2.0.0.sif samtools flagstat ${bam} > ${output_file}

echo "samtools flagstat completed"
date
echo "Output saved to ${output_file}"
