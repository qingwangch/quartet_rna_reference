#!/bin/bash
#SBATCH --job-name=fastp                          # Job name
#SBATCH --output=logs/fastp_%j.log                # Standard output and error log (%j = JobID)
#SBATCH --error=logs/fastp_%j.log                 # Merge stderr and stdout
#SBATCH --ntasks=1                                # Single task
#SBATCH --cpus-per-task=8                         # Number of CPU cores
#SBATCH --mem=16G                                 # Memory per node
#SBATCH --time=24:00:00                           # Time limit hrs:min:sec
#SBATCH --mail-user=chen.q@wehi.edu.au            # Email address for notifications
#SBATCH --mail-type=FAIL                          # Send email on failure

# Enable error detection
set -eo pipefail

# Load Singularity (version 3.5.2)
module load singularity

# Image directory
imgDIR=/vast/projects/quartet_rna_refdata/images/fastp_0.24.1.sif

# Working directory
workdir=/vast/projects/quartet_rna_refdata
cd "${workdir}"

# Log start
date
echo "fastp begins for sample: ${sample_id}"

start_time=$(date +%s)

# Function to create directory if it doesn't exist
create_dir_if_not_exists() {
  if [ ! -d "$1" ]; then
    mkdir -p "$1"
    echo "Directory created: $1"
  else
    echo "Directory already exists: $1"
  fi
}

# Define output directories
workdir1="/vast/projects/quartet_rna_refdata/analysis/fastp/${sample_id}"
workdir2="/vast/projects/quartet_rna_refdata/analysis/trim_fastq/${sample_id}"

create_dir_if_not_exists "${workdir1}"
create_dir_if_not_exists "${workdir2}"

# Quality filtering
singularity exec --bind "${workdir},/vast/scratch/users/chen.q" "${imgDIR}" \
  fastp \
    --thread 8 \
    --trim_front1 "${trim_front1}" \
    --trim_tail1  "${trim_tail1}" \
    --max_len1    "${max_len1}" \
    --trim_front2 "${trim_front2}" \
    --trim_tail2  "${trim_tail2}" \
    --max_len2    "${max_len2}" \
    -i "${read1}" \
    -I "${read2}" \
    -o "${workdir2}/${sample_id}_R1.fastq.tmp1.gz" \
    -O "${workdir2}/${sample_id}_R2.fastq.tmp1.gz" \
    -j "${workdir1}/${sample_id}.qf.json" \
    -h "${workdir1}/${sample_id}.qf.html"

# Adapter trimming
singularity exec --bind "${workdir},/vast/scratch/users/chen.q" "${imgDIR}" \
  fastp \
    --thread 8 \
    -l "${length_required}" \
    -q "${qualified_quality_phred}" \
    -u "${length_required1}" \
    --adapter_sequence       "${adapter_sequence}" \
    --adapter_sequence_r2    "${adapter_sequence_r2}" \
    --detect_adapter_for_pe  \
    --trim_front1            "${trim_front1}" \
    --trim_tail1             "${trim_tail1}" \
    --max_len1               "${max_len1}" \
    --trim_front2            "${trim_front2}" \
    --trim_tail2             "${trim_tail2}" \
    --max_len2               "${max_len2}" \
    -i "${workdir2}/${sample_id}_R1.fastq.tmp1.gz" \
    -I "${workdir2}/${sample_id}_R2.fastq.tmp1.gz" \
    -o "${workdir2}/${sample_id}_R1.fastq.gz" \
    -O "${workdir2}/${sample_id}_R2.fastq.gz" \
    -j "${workdir1}/${sample_id}.at.json" \
    -h "${workdir1}/${sample_id}.at.html"

# Clean up temporary files
rm "${workdir2}/${sample_id}_R1.fastq.tmp1.gz"
rm "${workdir2}/${sample_id}_R2.fastq.tmp1.gz"

# Log end
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "fastp ends for sample: ${sample_id}"
echo "Execution time: ${execution_time} seconds"
