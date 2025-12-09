#!/bin/bash

#SBATCH --cpus-per-task=40
#SBATCH --mem=256G
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A30:2
#SBATCH --time=48:00:00
#SBATCH --job-name=modkit_call_mods
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=chen.q@wehi.edu.au
#SBATCH --output=logs/dorado_modkit_call_mods_%x_%j.log

# 1️⃣ Load necessary modules
module load dorado
module load singularity

# 2️⃣ Parse named command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_id) SAMPLE_ID="$2"; shift ;;
        *) echo "❌ Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# 3️⃣ Check required arguments
if [ -z "$SAMPLE_ID" ]; then
  echo "❌ ERROR: --sample_id is required."
  echo "Usage: $0 --sample_id <SAMPLE_ID>"
  exit 1
fi

# 4️⃣ Define paths and parameters
IMG_DIR=/vast/projects/quartet_rna_refdata/images
BASECALL_DIR=/vast/scratch/users/chen.q/quartet_rna_refdata/analysis/m6a/${SAMPLE_ID}
POD5_DIR=/vast/scratch/users/xu.ya/huge_data/${SAMPLE_ID}/pod5_pass
REF_FA=/vast/scratch/users/xu.ya/quartet_rna_project/reference/gencode.v43.transcripts.fa

# 5️⃣ Perform basecalling with m6A detection
dorado basecaller /vast/scratch/users/xu.ya/quartet_rna_project/ONT-DR-GRM-R01/rna004_130bps_sup@v5.1.0 \
  --modified-bases m6A_DRACH \
  --reference $REF_FA \
  --output-dir $BASECALL_DIR \
  $POD5_DIR

# 6️⃣ Run modkit pileup to compute base modification levels
singularity exec --bind /vast/projects/quartet_rna_refdata/,/vast/scratch/users/xu.ya/,/vast/scratch/users/chen.q/ \
  ${IMG_DIR}/RNA_modification.sif \
  modkit pileup ${BASECALL_DIR}/*.bam ${BASECALL_DIR}/pileup.bed \
  --log-filepath ${BASECALL_DIR}/pileup.log

# 7️⃣ (Optional) Use call-mods to generate a modBAM with 0/100% calls
# singularity exec ${IMG_DIR}/RNA_modification.sif modkit \
#   call-mods ${BASECALL_DIR}/*.bam -o ${BASECALL_DIR}/${SAMPLE_ID}_called_mods.bam

echo "✅ modkit pileup completed for $SAMPLE_ID!"
