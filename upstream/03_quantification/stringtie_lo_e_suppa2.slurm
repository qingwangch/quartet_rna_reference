#!/bin/bash

#SBATCH --job-name=stringtie_lo_e_suppa2           # Job name
#SBATCH --nodes=1                                # Use 1 node
#SBATCH --ntasks=8                               # Allocate 8 CPU cores
#SBATCH --mem=16G  # Request 16GB of memory (adjust as needed)
#SBATCH --output=logs/stringtie_lo_e_suppa2_%j.log # Output log file with job ID
#SBATCH --time=24:00:00                          # Maximum runtime (24 hours)

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
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure required parameters are provided
if [ -z "$sample_id" ] || [ -z "$workdir" ]; then
    echo "Error: --sample_id and --workdir are required parameters."
    echo "Usage: sbatch stringtie_lo_e_suppa2.slurm --sample_id <ID> --workdir <path>"
    exit 1
fi

# Change to the working directory
cd /vast/projects/quartet_rna_refdata/

date
echo "Process begins: getTPM_LO and stringtie_lo_e-suppa"

start_time=$(date +%s)

# Define subdirectory paths
workdir1=${workdir}/analysis/suppa2/stringtie2/LO-E-G/${sample_id}

# Check or create necessary directories
for dir in "$workdir1" "$workdir1"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Directory created: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done

# Extract TPM values from the GTF file
gtf_file="/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"
ioe_file="/vast/projects/quartet_rna_refdata/reference/suppa2/gencode_v43.all_AS.ioe"
output_file="${workdir1}/${sample_id}_transcript_tpm.csv"

if [ ! -f "$gtf_file" ]; then
    echo "Error: $gtf_file not found!"
    exit 1
fi

# Extract TPM values and save them to a CSV file
workdir2=/vast/projects/quartet_rna_refdata/analysis/stringtie2/LO-E/${sample_id}
sample_gtf="${workdir2}/${sample_id}.gtf"
awk '/transcript_id/ && /TPM/ {gsub(/;/, "", $12); gsub(/;/, "", $NF); \
     split($12, transcript_id, "\""); split($NF, tpm, "\""); \
     print transcript_id[2] "\t" tpm[2]}' "$sample_gtf" > "$output_file"

echo "TPM extraction completed: $output_file"

# Merge and format the TPM file for SUPPA2
formatted_file="${workdir1}/${sample_id}.stringtie_lo.formatted.tpm.tsv"
echo "${sample_id}" > "${formatted_file}"
cat "${output_file}" >> "${formatted_file}"

# Execute SUPPA2 tasks using the Singularity image
## Perform PSI calculation
singularity exec "$imgDIR/suppa_2.3.sif" \
    suppa.py psiPerEvent -i "${ioe_file}" \
                         -e "${formatted_file}" \
                         -o "${workdir1}/${sample_id}_events"

singularity exec "$imgDIR/suppa_2.3.sif" \
    suppa.py psiPerIsoform -g "$gtf_file" \
                           -e "${formatted_file}" \
                           -o "${workdir1}/${sample_id}"

# Record the end time and calculate execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "Process ends: getTPM_LO-E and stringtie_lo_e-suppa"
echo "Total execution time: ${execution_time} seconds"
