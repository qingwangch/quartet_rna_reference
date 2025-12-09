#!/bin/bash
#SBATCH --job-name=star_index              # Job name
#SBATCH --mem=100G
#SBATCH --ntasks=2                         # Run on a single CPU
#SBATCH --cpus-per-task=8                  # Number of CPU cores
#SBATCH --time=24:00:00                    # Time limit hrs:min:sec
#SBATCH --output=logs/star_index_%j.log         # Standard output and error log
#SBATCH --error=logs/star_index_%j.err          # Error log

# Enable error detection and pipe failure detection
set -eo pipefail

# Load Singularity module
module load singularity

# Directory where the Singularity image is stored
imgDIR=/vast/projects/quartet_rna_refdata/images

# Move to the working directory
cd /vast/projects/quartet_rna_refdata/

date
echo "STAR genome index generation begins"

start_time=$(date +%s)

# Final command to generate the STAR genome index
singularity exec $imgDIR/star_2.7.10b.sif STAR --runThreadN 8 \
                                               --runMode genomeGenerate \
                                               --genomeDir /vast/projects/quartet_rna_refdata/ref_genome/star_index \
                                               --genomeFastaFiles /vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
                                               --sjdbGTFfile /vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf \
                                               --sjdbOverhang 149

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "STAR genome index generation ends"
echo "Execution time: $execution_time seconds"
