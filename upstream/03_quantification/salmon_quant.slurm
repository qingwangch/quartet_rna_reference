#!/bin/bash
#SBATCH --job-name=salmon_quant              # Job name
#SBATCH --output=logs/salmon_quant_%j.log    # Standard output and error log
#SBATCH --ntasks=1                           # Single task job
#SBATCH --cpus-per-task=8                    # Number of CPU cores
#SBATCH --mem=24G                            # Memory per node
#SBATCH --time=24:00:00                      # Time limit hrs:min:sec

# Enable error detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define paths
imgDIR=/vast/projects/quartet_rna_refdata/images/salmon_1.10.3.sif
indexDIR=/vast/projects/quartet_rna_refdata/ref_transcriptome/salmon_index/
outputDIR=/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/salmon_quant

# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Parse input arguments
reads1=""
reads2=""
libtype=""
sample_id=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -1) reads1="$2"; shift ;;
        -2) reads2="$2"; shift ;;
        -l) libtype="$2"; shift ;;
        --sample-id) sample_id="$2"; shift ;;
        --help)
            echo "Usage: $0 -1 reads1.fq -2 reads2.fq -l LIBTYPE --sample-id SAMPLE_ID"
            exit 1 ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Validate inputs
if [[ -z "$reads1" || -z "$reads2" || -z "$libtype" || -z "$sample_id" ]]; then
    echo "Error: -1, -2, -l, and --sample-id parameters are required."
    echo "Usage: $0 -1 reads1.fq -2 reads2.fq -l LIBTYPE --sample-id SAMPLE_ID"
    exit 1
fi

# Ensure output directory exists
mkdir -p ${outputDIR}/${sample_id}
cd /vast/projects/quartet_rna_refdata

# Log the start time
start_time=$(date +%s)
date
echo "Salmon quantification begins for sample: $sample_id"

# Run Salmon quant with Singularity
singularity exec --bind $BIND_PATHS $imgDIR \
    salmon quant -i ${indexDIR} \
                 -l $libtype \
                 --gcBias \
                 -p 8 \
                 -1 $reads1 \
                 -2 $reads2 \
                 --numBootstraps 50\
                 --validateMappings \
                 -o ${outputDIR}/${sample_id}

# Log the end time and execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "Salmon quantification ends for sample: $sample_id"
echo "Execution time: $execution_time seconds"
