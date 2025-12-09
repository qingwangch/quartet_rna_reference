#!/bin/bash
#SBATCH --job-name=alignqc_analysis      # Job name
#SBATCH --output=logs/alignqc_output_%j.log # Output log (with job ID)
#SBATCH --error=logs/alignqc_error_%j.log   # Error log (with job ID)
#SBATCH --time=24:00:00                  # Time limit
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=8                # CPUs per task
#SBATCH --mem=64G                        # Memory allocation

# Load Singularity if necessary
module load singularity

# Parse command-line arguments for BAM file, sample ID, and output directory
while getopts ":b:s:o:" opt; do
  case $opt in
    b) BAM_FILE="$OPTARG" ;;   # BAM file input
    s) SAMPLE_ID="$OPTARG" ;;  # Sample ID input
    o) OUTPUT_DIR="$OPTARG" ;; # Output directory
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Ensure required arguments are provided
if [ -z "$BAM_FILE" ] || [ -z "$SAMPLE_ID" ] || [ -z "$OUTPUT_DIR" ]; then
  echo "Usage: sbatch run_alignqc.sbatch -b <BAM_FILE> -s <SAMPLE_ID> -o <OUTPUT_DIR>"
  exit 1
fi

# Ensure the output directory exists
#mkdir -p "$OUTPUT_DIR"

# Define output files
OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_report.xhtml"
PORTABLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE_ID}_report.portable.xhtml"

# Track the start time
START_TIME=$(date +%s)
echo "Task started at: $(date)"
cd /vast/projects/quartet_rna_refdata

# Run alignQC inside Singularity with locale settings and bind path
singularity exec \
  --env LC_ALL=C.UTF-8,LANG=C.UTF-8 \
  -B /vast/projects/quartet_rna_refdata \
  /vast/projects/quartet_rna_refdata/images/alignqc_2.0.5.sif alignqc analyze \
    "$BAM_FILE" \
    -g /vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    -t /vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz \
    -o "$OUTPUT_FILE" \
    --portable_output ${PORTABLE_OUTPUT} \
    --output_folder ${OUTPUT_DIR} \
    --tempdir /vast/projects/quartet_rna_refdata/tmpdir/alignQC

# Track the end time
END_TIME=$(date +%s)
echo "Task finished at: $(date)"
RUNTIME=$((END_TIME - START_TIME))
echo "Total runtime: $RUNTIME seconds"
