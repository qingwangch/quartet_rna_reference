#!/bin/bash
#SBATCH --job-name=minimap2_pab         # Job name
#SBATCH --ntasks=1                      # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem=64G                       # Job memory request
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=logs/minimap2_%j.log        # Standard output and error log
#SBATCH --partition=bigmem 

# Enable error detection
set -eo pipefail

# Load Singularity
module load singularity/
imgDIR=/vast/projects/quartet_rna_refdata/images

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --workdir1) workdir1="$2"; shift ;;
        --sample_id) sample_id="$2"; shift ;;
        --fasta) fasta="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure all required arguments are provided
if [ -z "$workdir1" ] || [ -z "$sample_id" ] || [ -z "$fasta" ]; then
    echo "Usage: $0 --workdir1 <workdir1> --sample_id <sample_id> --fasta <fasta>"
    exit 1
fi

# Change to the working directory
cd /vast/projects/quartet_rna_refdata/

date
echo "minimap2_pab+samtools begins"

start_time=$(date +%s)

# Create directory if it doesn't exist
mkdir -p ${workdir1}

# Run minimap2
# Bind necessary directories
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

singularity exec --bind $BIND_PATHS $imgDIR/flair_2.0.0.sif minimap2 -t 8 \
    -ax splice:hq -uf \
    --junc-bed /vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.isoforms_junctions.bed \
    /vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    ${fasta} \
    -o ${workdir1}/${sample_id}.minimap2.aligned.all.sam

# Convert SAM to BAM
singularity exec --bind $BIND_PATHS $imgDIR/flair_2.0.0.sif \
    samtools view -hSb \
    ${workdir1}/${sample_id}.minimap2.aligned.all.sam > ${workdir1}/${sample_id}.minimap2.aligned.all.bam

# Filter BAM for quality >= 40 and primary alignments only
singularity exec --bind $BIND_PATHS $imgDIR/flair_2.0.0.sif \
    samtools view -q 40 -F 2304 \
    -hb ${workdir1}/${sample_id}.minimap2.aligned.all.bam > ${workdir1}/${sample_id}.minimap2.aligned.only_primary.bam

# Sort BAM and index
singularity exec --bind $BIND_PATHS $imgDIR/flair_2.0.0.sif \
    samtools sort ${workdir1}/${sample_id}.minimap2.aligned.only_primary.bam \
    -o ${workdir1}/${sample_id}.minimap2.sorted.only_primary.bam

singularity exec --bind $BIND_PATHS $imgDIR/flair_2.0.0.sif \
    samtools index ${workdir1}/${sample_id}.minimap2.sorted.only_primary.bam

rm ${workdir1}/${sample_id}.minimap2.aligned.all.sam ${workdir1}/${sample_id}.minimap2.aligned.only_primary.bam

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "minimap2_pab+samtools ends"

echo "Program execution time: $execution_time seconds"
