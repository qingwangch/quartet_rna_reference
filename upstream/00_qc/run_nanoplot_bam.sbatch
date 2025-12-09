#!/bin/bash
#SBATCH --job-name=nanoplot_analysis      # Job name
#SBATCH --output=logs/nanoplot_%j.log     # Combined output and error log (with job ID)
#SBATCH --time=12:00:00                   # Time limit
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=8                 # CPUs per task
#SBATCH --mem=16G                         # Memory allocation

# Load Singularity
module load singularity/

# Parse command-line arguments for bam, output directory, and sample ID
while getopts ":b:o:s:" opt; do
  case $opt in
    b) BAM_FILE="$OPTARG" ;;     # Input BAM file
    o) OUTPUT_DIR="$OPTARG" ;;   # Output directory
    s) SAMPLE_ID="$OPTARG" ;;    # Sample ID
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Ensure all required parameters are provided
if [ -z "$BAM_FILE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$SAMPLE_ID" ]; then
  echo "Usage: sbatch run_nanoplot.sbatch -b <BAM_FILE> -o <OUTPUT_DIR> -s <SAMPLE_ID>"
  exit 1
fi

# Set the image directory
IMG_DIR="/vast/projects/quartet_rna_refdata/images/nanoplot_latest.sif"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Move to the working directory
cd /vast/

# Record the start time
echo "NanoPlot begins at: $(date)"
start_time=$(date +%s)

# Bind necessary directories for Singularity
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Run NanoPlot inside the Singularity container
singularity exec -B "$BIND_PATHS" "$IMG_DIR" \
  NanoPlot \
  -t 8 \
   --verbose \
  --bam "$BAM_FILE" \
  --maxlength 10000 \
  -f png \
  --plots kde dot \
  --dpi 300 \
  -o "$OUTPUT_DIR" \
  -p "${SAMPLE_ID}_bam_"

# Record the end time and calculate execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "NanoPlot ends at: $(date)"
echo "Execution time: $execution_time seconds"
