#!/bin/bash
#SBATCH --job-name=star_mapping              # Job name
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --cpus-per-task=8                    # Number of CPU cores
#SBATCH --mem=80G                            # Memory allocation
#SBATCH --time=24:00:00                      # Time limit hrs:min:sec
#SBATCH --output=logs/star_mapping_%j.log         # Standard output and error log
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au   # Where to send the email notifications

# Enable error detection and pipe failure detection
set -eo pipefail

# Load Singularity module
module load singularity

# Define directories
imgDIR=/vast/projects/quartet_rna_refdata/images
genome_index=/vast/projects/quartet_rna_refdata/ref_genome/star_index
annotation_gtf=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf

# Parse named parameters
for arg in "$@"
do
    case $arg in
        --fq_r1=*)
        fq_r1="${arg#*=}"
        shift
        ;;
        --fq_r2=*)
        fq_r2="${arg#*=}"
        shift
        ;;
        --sample_id=*)
        sample_id="${arg#*=}"
        shift
        ;;
        *)
        echo "Invalid argument: $arg"
        exit 1
        ;;
    esac
done

# Ensure all required arguments are provided
if [ -z "$fq_r1" ] || [ -z "$fq_r2" ] || [ -z "$sample_id" ]; then
    echo "Error: Missing required arguments."
    echo "Usage: sbatch star_mapping_slurm.sh --fq_r1=<R1_file> --fq_r2=<R2_file> --sample_id=<sample_id>"
    exit 1
fi

# Set working directory
cd /vast/projects/quartet_rna_refdata/

date
echo "STAR mapping begins for sample: $sample_id"

start_time=$(date +%s)

# Create output directory for STAR mapping results
mkdir -p /vast/scratch/users/chen.q/quartet_rna_refdata/analysis/star/${sample_id}/

# Run STAR mapping with Singularity
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

singularity exec -B "${BIND_PATHS}" $imgDIR/star_2.7.10b.sif STAR --genomeDir ${genome_index} \
                                                --readFilesIn ${fq_r1} ${fq_r2} \
                                                --runThreadN 8 \
                                                --outFilterMultimapScoreRange 1 \
                                                --outFilterMultimapNmax 20 \
                                                --outFilterMismatchNmax 999 \
                                                --alignIntronMax 1000000 \
                                                --alignMatesGapMax 1000000 \
                                                --sjdbScore 1 \
                                                --alignSJDBoverhangMin 1 \
                                                --genomeLoad NoSharedMemory \
                                                --limitBAMsortRAM 80000000000 \
                                                --readFilesCommand zcat \
                                                --sjdbOverhang 149 \
                                                --outSAMstrandField intronMotif \
                                                --outSAMattributes NH HI NM MD AS jM jI XS \
                                                --sjdbGTFfile ${annotation_gtf} \
                                                --limitSjdbInsertNsj 2000000 \
                                                --outSAMunmapped Within \
                                                --outSAMtype BAM SortedByCoordinate \
                                                --outSAMheaderHD @HD VN:1.4 \
                                                --outSAMattrRGline ID::${sample_id} \
                                                --twopassMode Basic \
                                                --outSAMmultNmax 1 \
                                                --outFileNamePrefix /vast/scratch/users/chen.q/quartet_rna_refdata/analysis/star/${sample_id}/${sample_id}_ \
                                                --outTmpDir /vast/scratch/users/chen.q/quartet_rna_refdata/analysis/star/${sample_id}/tmp

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "STAR mapping ends for sample: $sample_id"

echo "Execution time: $execution_time seconds"
