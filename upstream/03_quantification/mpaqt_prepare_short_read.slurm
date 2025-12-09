#!/bin/bash
#SBATCH --job-name=mpaqt_short_read
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=logs/mpaqt_prepare_short_read_%j.log
#SBATCH --error=logs/mpaqt_prepare_short_read_%j.log

set -eo pipefail

# Load Singularity module
module load singularity

# ------------------------ Parse arguments ------------------------ #
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --project) project="$2"; shift ;;
    --sample) sample="$2"; shift ;;
    --fastq1) fastq1="$2"; shift ;;
    --fastq2) fastq2="$2"; shift ;;
    *) echo "❌ Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [[ -z "$project" || -z "$sample" || -z "$fastq1" || -z "$fastq2" ]]; then
  echo "Usage: sbatch $0 --project <project_name> --sample <sample_id> --fastq1 <read1.fq.gz> --fastq2 <read2.fq.gz>"
  exit 1
fi

# ------------------------ Setup ------------------------ #
imgDIR=/vast/projects/quartet_rna_refdata/images
mpaqt_root=/vast/projects/quartet_rna_refdata/analysis/mpaqt
bus_out="$mpaqt_root/projects/${project}/data/bus/${sample}"
kallisto_index="/vast/projects/quartet_rna_refdata/reference/mpaqt/gencode.v43.index.kallisto"

start_time=$(date +%s)
echo "🚀 [$sample] Processing starts at: $(date)"
echo "📁 Project: $project"
echo "🧪 Sample: $sample"
echo "📄 FASTQ1: $fastq1"
echo "📄 FASTQ2: $fastq2"

# ------------------------ Run kallisto bus ------------------------ #
echo "▶️ Running kallisto bus..."
mkdir -p "$bus_out"

singularity exec --bind $mpaqt_root/projects:/csglab/mpaqt/projects,$mpaqt_root/tmpdir:/csglab/mpaqt/tmp,/vast/projects/quartet_rna_refdata \
  $imgDIR/mpaqt_0.3.1.sif kallisto bus \
  --index $kallisto_index \
  --num \
  --paired \
  --technology bulk \
  --output-dir ${bus_out} \
  ${fastq1} \
  ${fastq2}

echo "✅ kallisto bus finished for $sample"

# ------------------------ Run MPAQT prepare ------------------------ #
echo "▶️ Running mpaqt prepare short-read..."

singularity exec --bind $mpaqt_root/projects:/csglab/mpaqt/projects, \
                   --bind $mpaqt_root/tmpdir:/csglab/mpaqt/tmp, \
                   --bind /vast/projects/quartet_rna_refdata:/vast/projects/quartet_rna_refdata \
    $imgDIR/mpaqt_0.3.1.sif /csglab/mpaqt/mpaqt prepare short-read \
  --project "$project" \
  --sample "$sample" \
  --bus "${bus_out}/output.bus" \
  --matrix_ec "${bus_out}/matrix.ec"

echo "✅ mpaqt short-read preparation complete for $sample"

# ------------------------ Wrap-up ------------------------ #
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "⏰ [$sample] Finished at: $(date)"
echo "🕒 Total runtime: ${duration} seconds"
