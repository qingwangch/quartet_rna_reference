#!/bin/bash
#SBATCH --job-name=stringtie_lo_e        # Job name
#SBATCH --ntasks=1                      # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem=32G                       # Job memory request
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=logs/stringtie_lo_e_%j.log    # Standard output and error log

# Enable error detection
set -eo pipefail

# Load the necessary Singularity module
module load singularity/

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_id) sample_id="$2"; shift ;;
        --bam) bam="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure both arguments are provided
if [ -z "$sample_id" ] || [ -z "$bam" ]; then
    echo "Usage: $0 --sample_id <sample_id> --bam <bam_file>"
    exit 1
fi

# Define the image directory and the working directory
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast/projects/quartet_rna_refdata
cd ${workdir}

# Start time and process information
date
echo "stringtie_lo_e begins"
start_time=$(date +%s)

# Define the working directory for the sample
workdir1=${workdir}/analysis/stringtie2/LO-E/${sample_id}

# Check if the directory exists, and create it if it doesn't
if [ ! -d "$workdir1" ]; then
    mkdir -p "$workdir1"
    echo "Directory created: $workdir1"
else
    echo "Directory already exists: $workdir1"
fi

# Path to the GTF file
gtf=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf

# Run StringTie with Singularity
singularity exec $imgDIR/stringtie2_latest.sif stringtie -G ${gtf} -e -B -p 8 -L --conservative -l ${sample_id} -o ${workdir1}/${sample_id}.gtf -A ${workdir1}/${sample_id}.gene.abundance.txt ${bam}

# End time and execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "stringtie_lo_e ends"
echo "Program execution time: $execution_time seconds"