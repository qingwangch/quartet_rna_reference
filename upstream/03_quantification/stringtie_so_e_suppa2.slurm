#!/bin/bash

#SBATCH --job-name=stringtie_so_e_suppa2           # Job name
#SBATCH --nodes=1                                # Use 1 node
#SBATCH --ntasks=8                               # Allocate 8 CPU cores
#SBATCH --mem=16G                                # Request 16GB of memory
#SBATCH --output=logs/stringtie_so_e_suppa2_%j.log # Output log file with job ID
#SBATCH --time=24:00:00                          # Maximum runtime
#SBATCH --mail-type=FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au       # Where to send email notifications

set -eo pipefail  # Exit script on any error

module purge
module load singularity

# Define the location of the Singularity image
imgDIR=/vast/projects/quartet_rna_refdata/images

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_id) sample_id="$2"; shift ;;
        --workdir) workdir="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z "$sample_id" ] || [ -z "$workdir" ]; then
    echo "Error: --sample_id and --workdir are required parameters."
    exit 1
fi

cd /vast/scratch/users/chen.q/quartet_rna_refdata/

date
echo "Process begins: getTPM_SO and stringtie_so_e-suppa"

start_time=$(date +%s)

# Define subdirectory paths
workdir1=${workdir}/analysis/suppa2/stringtie2/SO-E/${sample_id}

# Check or create necessary directory
if [ ! -d "$workdir1" ]; then
    mkdir -p "$workdir1"
    echo "Directory created: $workdir1"
else
    echo "Directory already exists: $workdir1"
fi

# Extract GTF file path dynamically
sample_gtf=$(find /vast/scratch/users/chen.q/quartet_rna_refdata/analysis/stringtie2/SO/ -type f -name "${sample_id}.gtf" | head -n 1)

if [ -z "$sample_gtf" ]; then
    echo "Error: GTF file for sample ID ${sample_id} not found!"
    exit 1
fi

echo "Found GTF file: $sample_gtf"

# Extract TPM values
output_file="${workdir1}/${sample_id}_transcript_tpm.csv"
awk '/transcript_id/ && /TPM/ {gsub(/;/, "", $12); gsub(/;/, "", $NF); \
     split($12, transcript_id, "\""); split($NF, tpm, "\""); \
     print transcript_id[2] "\t" tpm[2]}' "$sample_gtf" > "$output_file"
echo "TPM extraction completed: $output_file"

# Merge and format the TPM file for SUPPA2
formatted_file="${workdir1}/${sample_id}.stringtie_so.formatted.tpm.tsv"
echo "${sample_id}" > "${formatted_file}"
cat "${output_file}" >> "${formatted_file}"

# Execute SUPPA2 tasks
ioe_file="/vast/projects/quartet_rna_refdata/reference/suppa2/gencode_v43.all_AS.ioe"
gtf_file="/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"
# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

singularity exec --bind $BIND_PATHS "$imgDIR/suppa_2.3.sif" \
    suppa.py psiPerEvent -i "${ioe_file}" \
                         -e "${formatted_file}" \
                         -o "${workdir1}/${sample_id}_events"

singularity exec --bind $BIND_PATHS "$imgDIR/suppa_2.3.sif" \
    suppa.py psiPerIsoform -g "$gtf_file" \
                           -e "${formatted_file}" \
                           -o "${workdir1}/${sample_id}"

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "Process ends: getTPM_SO-E and stringtie_so_e-suppa"
echo "Total execution time: ${execution_time} seconds"
