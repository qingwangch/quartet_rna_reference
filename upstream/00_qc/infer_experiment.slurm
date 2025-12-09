#!/bin/bash
#SBATCH --job-name=infer_experiment          # Job name
#SBATCH --ntasks=1                           # Run on a single CPU
#SBATCH --cpus-per-task=4                    # Number of CPU cores
#SBATCH --mem=16G                            # Memory allocation
#SBATCH --time=04:00:00                      # Time limit hrs:min:sec
#SBATCH --output=infer_experiment_%j.log     # Standard output and error log

# Enable error detection and pipe failure detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define paths
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast/projects/quartet_rna_refdata/
bed=/vast/projects/quartet_rna_refdata/ref_transcriptome/hg38_GENCODE_V43_Comprehensive.bed

# Parse named parameters
for arg in "$@"
do
    case $arg in
        --sample_id=*)
        sample_id="${arg#*=}"
        shift
        ;;
        --bam_list=*)
        bam_list="${arg#*=}"
        shift
        ;;
        *)
        echo "Invalid argument: $arg"
        exit 1
        ;;
    esac
done

# Check if both parameters are provided
if [ -z "$sample_id" ] || [ -z "$bam_list" ]; then
    echo "Error: Missing required parameters."
    echo "Usage: sbatch infer_experiment_slurm.sh --sample_id=<sample_id> --bam_list=<bam_list>"
    exit 1
fi

# Change to the working directory
cd ${workdir}

date
echo "infer_experiment begins"

start_time=$(date +%s)

# Output directory
output=/vast/projects/quartet_rna_refdata/analysis/rseqc/${sample_id}

# Check if the output directory exists, create it if not
if [ ! -d "$output" ]; then
    mkdir -p "$output"
    echo "Directory created: $output"
else
    echo "Directory already exists: $output"
fi

# Run infer_experiment using Singularity
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

singularity exec -B "${BIND_PATHS}" $imgDIR/rseqc_v1.sif /usr/local/bin/infer_experiment.py -r ${bed} -i ${bam_list} > ${output}/${sample_id}_infer_experiment.txt

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "infer_experiment ends"

echo "Execution time: $execution_time seconds"
