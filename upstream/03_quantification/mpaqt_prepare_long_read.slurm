#!/bin/bash
#SBATCH --job-name=mpaqt_long_read
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=logs/mpaqt_prepare_long_read_%j.log
#SBATCH --error=logs/mpaqt_prepare_long_read_%j.log

set -eo pipefail

# Load Singularity module
module load singularity

# ------------------------ Parse arguments ------------------------ #
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --project) project="$2"; shift ;;
    --sample) sample="$2"; shift ;;
    --counts) counts="$2"; shift ;;
    *) echo "❌ Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [[ -z "$project" || -z "$sample" || -z "$counts" ]]; then
  echo "Usage: sbatch $0 --project <project_name> --sample <sample_id> --counts <counts_file>"
  exit 1
fi

# ------------------------ Setup ------------------------ #
imgDIR=/vast/projects/quartet_rna_refdata/images
mpaqt_root=/vast/projects/quartet_rna_refdata/analysis/mpaqt

# 请确保 counts 文件的路径为宿主机上的绝对路径（如 /vast/projects/quartet_rna_refdata/analysis/mpaqt/projects/quartet/data/lr/...)
echo "🚀 [$sample] Long-read processing starts at: $(date)"
echo "📁 Project: $project"
echo "🧪 Sample: $sample"
echo "📄 Counts file: $counts"

# ------------------------ Prepare Long-read RNA-seq Data ------------------------ #
echo "▶️ Running mpaqt prepare long-read..."

# 绑定整个参考目录，mpaqt_root 和参考数据目录，以保证容器内能正确访问
singularity exec --bind $mpaqt_root/projects:/csglab/mpaqt/projects, \
                   --bind $mpaqt_root/tmpdir:/csglab/mpaqt/tmp, \
                   --bind /vast/projects/quartet_rna_refdata:/vast/projects/quartet_rna_refdata \
    $imgDIR/mpaqt_0.3.1.sif /csglab/mpaqt/mpaqt prepare long-read \
      --project "$project" \
      --sample "$sample" \
      --counts "$counts"

echo "✅ mpaqt long-read preparation complete for $sample"

# ------------------------ Wrap-up ------------------------ #
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "⏰ [$sample] Finished at: $(date)"
echo "🕒 Total runtime: ${duration} seconds"
