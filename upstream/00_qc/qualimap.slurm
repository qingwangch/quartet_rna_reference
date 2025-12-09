#!/bin/bash

#SBATCH --job-name=qualimap                # Job name
#SBATCH --nodes=1                           # Run on a single node
#SBATCH --ntasks=1                          # Number of tasks
#SBATCH --cpus-per-task=16                    # Number of CPU cores
#SBATCH --mem=64G                            # Memory allocation
#SBATCH --time=48:00:00                      # Time limit hrs:min:sec
#SBATCH --output=logs/qualimap_%j.log       # Standard output log (Job ID included)
#SBATCH --error=logs/qualimap_%j.log        # Standard error log (Job ID included)

# Enable error detection and pipe failure detection
set -eo pipefail

# Read command-line named parameters
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_id) sample_id="$2"; shift ;;   # Sample ID
        --bam) bam="$2"; shift ;;               # BAM file path
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Set defaults for optional parameters if not provided
seq_protocol="non-strand-specific"
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast/projects/quartet_rna_refdata
ref_gtf=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf

# Load Singularity module
module load singularity

# Define working directory for this sample
cd $workdir

workdir1=/vast/projects/quartet_rna_refdata/analysis/qualimap/${sample_id}
mkdir -p "$workdir1"
echo "Working directory: $workdir1"

date
echo "qualimap begins"

start_time=$(date +%s)

# Run qualimapBAMqc
echo "qualimapBAMqc begins"
singularity exec $imgDIR/qualimap_v2.3.sif qualimap bamqc -bam ${bam} -c -outformat PDF:HTML -nt 16 -p ${seq_protocol} -outdir ${workdir1}/${sample_id}_bamqc --java-mem-size=32G

# Compress results
tar -zcvf ${workdir1}/${sample_id}_bamqc_qualimap.tar.gz ${workdir1}/${sample_id}_bamqc
echo "qualimapBAMqc ends"

# Run qualimapRNAseq
echo "qualimapRNAseq begins"
singularity exec $imgDIR/qualimap_v2.3.sif qualimap rnaseq -bam ${bam} -outformat HTML -outdir ${workdir1}/${sample_id}_RNAseq -gtf ${ref_gtf} --java-mem-size=32G

# Compress results
tar -zcvf ${workdir1}/${sample_id}_RNAseq_qualimap.tar.gz ${workdir1}/${sample_id}_RNAseq
echo "qualimapRNAseq ends"

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "qualimap ends"
echo "Execution time: $execution_time seconds"
