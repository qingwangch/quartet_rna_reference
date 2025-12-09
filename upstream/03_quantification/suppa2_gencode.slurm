#!/bin/bash

#SBATCH --job-name=suppa2_gencode           # Job name
#SBATCH --nodes=1                           # Use 1 node
#SBATCH --ntasks=16                          # Allocate 16 CPU cores
#SBATCH --mem=16G  # Request 16GB of memory (adjust as needed)
#SBATCH --output=logs/suppa2_gencode_%j.log # Output log file with job ID
#SBATCH --time=24:00:00                     # Maximum runtime (24 hours)

# Exit the script on any error
set -eo pipefail

# Load Singularity environment and clear previous modules
module purge
module load singularity

# Define the location of the Singularity image
imgDIR=/vast/projects/quartet_rna_refdata/images

# Change to the working directory
cd /vast/projects/quartet_rna_refdata/

date
echo "Process begins: suppa2 generate events"

start_time=$(date +%s)

# Define subdirectory paths
workdir1=/vast/projects/quartet_rna_refdata/reference/suppa2

# Check or create necessary directories
for dir in "$workdir1"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Directory created: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done

# Extract TPM values from the GTF file
gtf_file="/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"

if [ ! -f "$gtf_file" ]; then
    echo "Error: $gtf_file not found!"
    exit 1
fi

# Execute SUPPA2 tasks using the Singularity image

## Generate AS events
singularity exec "$imgDIR/suppa_2.3.sif" \
    suppa.py generateEvents -t 16 \
                            -i "$gtf_file" \
                            -o "${workdir1}/gencode_v43.AS" \
                            -f ioe \
                            -e SE SS MX RI FL

# Merge AS event files
cat "${workdir1}/gencode_v43.AS_"*_strict.ioe | \
    awk 'NR==1 || (NR!=1 && $0!~/^seqname/)' > "${workdir1}/gencode_v43.all_AS.ioe"


# Record the end time and calculate execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "Process ends: suppa2 generate events"
echo "Total execution time: ${execution_time} seconds"
