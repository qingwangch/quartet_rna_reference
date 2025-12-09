#!/bin/bash

#SBATCH --job-name=salmon_suppa2                # Job name
#SBATCH --nodes=1                               # Use 1 node
#SBATCH --ntasks=8                              # Allocate 8 CPU cores
#SBATCH --mem=16G                               # Request 16GB of memory (adjust as needed)
#SBATCH --output=logs/salmon_suppa2_%j.log      # Output log file with job ID
#SBATCH --time=24:00:00                         # Maximum runtime (24 hours)

# Exit the script on any error
set -eo pipefail

# Load Singularity environment and clear previous modules
module purge
module load singularity

# Define the location of the Singularity image
imgDIR=/vast/projects/quartet_rna_refdata/images

# Parse command-line arguments (named parameters)
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_id) sample_id="$2"; shift ;;
        --workdir) workdir="$2"; shift ;;
        --quant_file) quant_file="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure required parameters are provided
if [ -z "$sample_id" ] || [ -z "$workdir" ] || [ -z "$quant_file" ]; then
    echo "Error: --sample_id, --workdir, and --quant_file are required parameters."
    echo "Usage: sbatch salmon_suppa2.slurm --sample_id <ID> --workdir <path> --quant_file <path_to_quant.sf>"
    exit 1
fi

# Change to the working directory
cd /vast/scratch/users/chen.q/quartet_rna_refdata/

date
echo "Process begins: salmon-suppa2"

start_time=$(date +%s)

# Define subdirectory paths
workdir1=${workdir}/analysis/suppa2/salmon/${sample_id}

# Check or create necessary directories
if [ ! -d "$workdir1" ]; then
    mkdir -p "$workdir1"
    echo "Directory created: $workdir1"
else
    echo "Directory already exists: $workdir1"
fi

# Check if the quant.sf file exists
if [ ! -f "$quant_file" ]; then
    echo "Error: $quant_file not found!"
    exit 1
fi

# Extract the Name and TPM columns, set the header to the sample_id, and save to a new formatted file
formatted_file="${workdir1}/${sample_id}_formatted_tpm.txt"
{
    echo -e "${sample_id}"
    awk 'NR > 1 {print $1 "\t" $4}' "$quant_file"
} > "$formatted_file"

# Define paths to reference files
gtf_file="/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"
ioe_file="/vast/projects/quartet_rna_refdata/reference/suppa2/gencode_v43.all_AS.ioe"

# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Check if the GTF file exists
if [ ! -f "$gtf_file" ]; then
    echo "Error: $gtf_file not found!"
    exit 1
fi

# Execute SUPPA2 tasks using the Singularity image
## Perform PSI calculation for events
singularity exec --bind $BIND_PATHS "$imgDIR/suppa_2.3.sif" \
    suppa.py psiPerEvent -i "${ioe_file}" \
                         -e "${formatted_file}" \
                         -o "${workdir1}/${sample_id}_events"

## Perform PSI calculation for isoforms
singularity exec --bind $BIND_PATHS "$imgDIR/suppa_2.3.sif" \
    suppa.py psiPerIsoform -g "$gtf_file" \
                           -e "${formatted_file}" \
                           -o "${workdir1}/${sample_id}"

# Record the end time and calculate execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "Process ends: salmon-suppa2"
echo "Total execution time: ${execution_time} seconds"
