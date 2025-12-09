#!/bin/bash
#SBATCH --job-name=stringtie2_so             # Job name
#SBATCH --ntasks=1                           # Run on a single CPU
#SBATCH --cpus-per-task=8                    # Number of CPU cores
#SBATCH --mem=32G                            # Memory allocation
#SBATCH --time=12:00:00                      # Time limit hrs:min:sec
#SBATCH --output=logs/stringtie2_so_%j.log        # Standard output and error log
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au   # Where to send the email notifications

# Enable error detection and pipe failure detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define paths
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast
gtf=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf

# Parse named parameters (strand specificity and sample_id)
for arg in "$@"
do
    case $arg in
        --sample_id=*)
        sample_id="${arg#*=}"
        shift
        ;;
        --bam=*)
        bam="${arg#*=}"
        shift
        ;;
        --strand_specific=*)
        strand_specific="${arg#*=}"
        shift
        ;;
        *)
        echo "Invalid argument: $arg"
        exit 1
        ;;
    esac
done

# Ensure sample_id, bam, and strand_specific are provided
if [ -z "$sample_id" ] || [ -z "$bam" ] || [ -z "$strand_specific" ]; then
    echo "Error: Missing required arguments."
    echo "Usage: sbatch stringtie2_so_slurm.sh --sample_id=<sample_id> --bam=<bam_file> --strand_specific=<SO_NS|SO_SS>"
    exit 1
fi

# Change to the working directory
cd ${workdir}

date
echo "stringtie2 begins for sample: ${sample_id} in mode: ${strand_specific}"

start_time=$(date +%s)

workdir1=/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/stringtie2/SO/${sample_id}

# Check if the output directory exists, create it if not
if [ ! -d "$workdir1" ]; then
    mkdir -p "$workdir1"
    echo "Directory created: $workdir1"
else
    echo "Directory already exists: $workdir1"
fi

# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Run StringTie with or without strand specificity depending on the argument
if [ "$strand_specific" == "SO_SS" ]; then
    singularity exec --bind $BIND_PATHS $imgDIR/stringtie2_latest.sif stringtie -e -G ${gtf} -B -p 8 --rf -o ${workdir1}/${sample_id}.gtf -A ${workdir1}/${sample_id}.gene.abundance.txt ${bam}
else
    singularity exec --bind $BIND_PATHS $imgDIR/stringtie2_latest.sif stringtie -e -G ${gtf} -B -p 8 -o ${workdir1}/${sample_id}.gtf -A ${workdir1}/${sample_id}.gene.abundance.txt ${bam}
fi

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "stringtie2 ends for sample: ${sample_id}"

echo "Execution time: $execution_time seconds"
