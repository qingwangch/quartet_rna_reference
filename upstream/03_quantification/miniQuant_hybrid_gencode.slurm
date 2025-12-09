#!/bin/bash

#SBATCH --job-name=miniQuant_hybrid_gencode     # Job name
#SBATCH --output=logs/miniQuant_hybrid_gencode_%j.out    # Combined stdout and stderr
#SBATCH --partition=bigmem               # Use the bigmem partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128gb
#SBATCH --time=48:00:00

# Load Singularity
module load singularity

# Working directory and image
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast/projects/quartet_rna_refdata

# Parse named parameters passed from command line
for arg in "$@"
do
    case $arg in
        sample_id=*)
            sample_id="${arg#*=}"
            shift
            ;;
        long_reads=*)
            long_reads="${arg#*=}"
            shift
            ;;
        short_r1=*)
            short_r1="${arg#*=}"
            shift
            ;;
        short_r2=*)
            short_r2="${arg#*=}"
            shift
            ;;
        long_reads_library_prep=*)
            long_reads_library_prep="${arg#*=}"
            shift
            ;;
        short_reads_strandness=*)
            short_reads_strandness="${arg#*=}"
            shift
            ;;
        *)
            echo "Invalid argument: $arg"
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$sample_id" ] || [ -z "$long_reads" ] || [ -z "$short_r1" ] || [ -z "$short_r2" ] || [ -z "$long_reads_library_prep" ] || [ -z "$short_reads_strandness" ]; then
    echo "Error: Missing required parameters."
    echo "Usage: sbatch miniQuant_hybrid_gencode.slurm sample_id=<sample_id> long_reads=<long_reads.fa> short_r1=<SR_R1.fq> short_r2=<SR_R2.fq> long_reads_library_prep=<prep_type> short_reads_strandness=<strandness>"
    exit 1
fi

# Create output directory
output_dir=${workdir}/analysis/miniquant_gencode/${sample_id}
mkdir -p "$output_dir"

# Record start time
echo "miniQuant hybrid quantification begins"
date
start_time=$(date +%s)

# Define reference
ref_fa=${workdir}/ref_transcriptome/gencode.v43.transcripts.fa

# Run miniQuant (assuming singularity image is properly packed)
singularity exec -B /vast/projects/quartet_rna_refdata/,/vast/scratch/users/chen.q/quartet_rna_refdata/ ${imgDIR}/miniquant_v1.1.sif \
    /opt/miniQuant_linux/miniQuant quant \
    -r ${ref_fa} \
    -l ${long_reads} \
    --long_reads_library_prep ${long_reads_library_prep} \
    -1 ${short_r1} \
    -2 ${short_r2} \
    --short_reads_strandness ${short_reads_strandness} \
    -t 4 \
    -o ${output_dir}

rm -rf ${output_dir}/temp

# Record end time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "miniQuant hybrid quantification ends"
date
echo "Execution time: ${execution_time} seconds"
