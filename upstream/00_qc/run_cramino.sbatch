#!/bin/bash
#SBATCH --job-name=cramino_analysis      # Job name
#SBATCH --output=logs/only_primary/cramino_%j.log  # Combined output and error log (with job ID)
#SBATCH --time=12:00:00                   # Time limit
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=8                 # CPUs per task
#SBATCH --mem=4G                         # Memory allocation

# Load Singularity (or appropriate module)
module load singularity

# Parse BAM file and sample ID from command-line arguments
while getopts ":b:s:" opt; do
  case $opt in
    b) BAM_FILE="$OPTARG" ;;  # BAM file input
    s) SAMPLE_ID="$OPTARG" ;; # Sample ID input
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Ensure both BAM file and sample ID are provided
if [ -z "$BAM_FILE" ] || [ -z "$SAMPLE_ID" ]; then
  echo "Usage: sbatch run_cramino.sbatch -b <BAM_FILE> -s <SAMPLE_ID>"
  exit 1
fi

# Create a temporary working directory
cd /vast/projects/quartet_rna_refdata/

# Track the start time
START_TIME=$(date +%s)
echo "Task started at: $(date)"

# Run Cramino analysis inside Singularity
singularity exec -B /vast/projects/quartet_rna_refdata/ \
  /vast/projects/quartet_rna_refdata/images/cramino_0.14.5.sif cramino \
  --threads 8 \
  --reference /vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
  --hist \
  --spliced \
  "$BAM_FILE"

# # Move results to the designated output directory
# OUTPUT_DIR="/path/to/output_directory/${SAMPLE_ID}"
# mkdir -p "$OUTPUT_DIR"
# mv * "$OUTPUT_DIR"

# Track the end time
END_TIME=$(date +%s)
echo "Task finished at: $(date)"
RUNTIME=$((END_TIME - START_TIME))
echo "Total runtime: $RUNTIME seconds"
