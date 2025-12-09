#!/bin/bash
#SBATCH --job-name=mpaqt_quant
#SBATCH --partition=regular               # Use the regular partition
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=logs/mpaqt_quant_%j.log
#SBATCH --error=logs/mpaqt_quant_%j.log

set -eo pipefail

# Load Singularity module
module load singularity

# ------------------------ Parse arguments ------------------------ #
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --project) project="$2"; shift ;;
    --sample) sample="$2"; shift ;;
    *) echo "❌ Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [[ -z "$project" || -z "$sample" ]]; then
  echo "Usage: sbatch $0 --project <project_name> --sample <sample_id>"
  exit 1
fi

# ------------------------ Setup ------------------------ #
# 定义镜像目录和 MPAQT 根目录（请根据实际路径修改）
imgDIR=/vast/projects/quartet_rna_refdata/images
mpaqt_root=/vast/projects/quartet_rna_refdata/analysis/mpaqt

start_time=$(date +%s)
echo "🚀 [$sample] Quantification starts at: $(date)"
echo "📁 Project: $project"
echo "🧪 Sample: $sample"

# ------------------------ Run MPAQT quant ------------------------ #
# 绑定目录：
#   - $mpaqt_root 绑定到容器内 /csglab 以保证容器能访问 MPAQT 相关数据
#   - /vast/projects/quartet_rna_refdata 挂载，使容器能访问参考数据和其它必需文件
singularity exec --bind $mpaqt_root/projects:/csglab/mpaqt/projects, \
                   --bind $mpaqt_root/tmpdir:/csglab/mpaqt/tmp, \
                   --bind /vast/projects/quartet_rna_refdata:/vast/projects/quartet_rna_refdata \
  $imgDIR/mpaqt_0.3.1.sif /csglab/mpaqt/mpaqt quant \
    --project "$project" \
    --sample "$sample"

echo "✅ mpaqt quantification complete for $sample"

# ------------------------ Wrap-up ------------------------ #
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "⏰ [$sample] Finished at: $(date)"
echo "🕒 Total runtime: ${duration} seconds"
