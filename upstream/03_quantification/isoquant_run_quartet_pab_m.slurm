#!/bin/bash
#SBATCH --job-name=isoquant_run_quartet_pab_m     # Job name
#SBATCH --ntasks=1                         # Run on a single CPU
#SBATCH --cpus-per-task=8                  # Number of CPU cores
#SBATCH --mem=128G                          # Memory allocation
#SBATCH --time=48:00:00                    # Time limit hrs:min:sec
#SBATCH --output=logs/iisoquant_run_quartet_pab_m_%j.log           # Standard output and error log

# Define paths
source /home/users/allstaff/chen.q/miniconda3/bin/activate isoquant

REFERENCE_GENOME=/vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa  # Path to reference genome in FASTA format
ANNOTATION_FILE=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.db    # Path to gene annotation in GTF format

# Parse named parameters
for arg in "$@"
do
    case $arg in
        --output_dir=*)
        OUTPUT_DIR="${arg#*=}"
        shift
        ;;
        --yaml_file=*)
        YAML_FILE="${arg#*=}"
        shift
        ;;
        --data_type=*)
        DATA_TYPE="${arg#*=}"
        shift
        ;;
        *)
        echo "Invalid argument: $arg"
        exit 1
        ;;
    esac
done

# Check if required parameters are provided
if [ -z "$OUTPUT_DIR" ] || [ -z "$YAML_FILE" ] || [ -z "$DATA_TYPE" ]; then
    echo "Error: Missing required parameters."
    echo "Usage: sbatch isoquant_run_quartet_ont_p.slurm --output_dir=<output_dir> --yaml_file=<yaml_file> --data_type=<data_type>"
    exit 1
fi

# Define other IsoQuant options
THREADS=8                                  # Number of threads to use

# Run IsoQuant
isoquant.py \
  --output $OUTPUT_DIR \
  --prefix "quartet" \
  --yaml $YAML_FILE \
  --reference $REFERENCE_GENOME \
  --genedb $ANNOTATION_FILE \
  --complete_genedb \
  --data_type $DATA_TYPE \
  --threads $THREADS \
  --sqanti_output \
  --count_exons

echo "IsoQuant analysis completed"
