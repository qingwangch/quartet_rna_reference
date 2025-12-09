#!/bin/bash
#SBATCH --job-name=preseq
#SBATCH --partition=regular          # 按集群实际情况修改
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/preseq_%x_%j.log
#SBATCH --error=logs/preseq_%x_%j.log

# ---------------- 0. 命名参数解析 ----------------
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --sample_id) SAMPLE="$2";   shift 2 ;;
    --bam)       BAM="$2";      shift 2 ;;
    *) echo "Unknown argument $1"; exit 1 ;;
  esac
done

if [[ -z "${SAMPLE:-}" || -z "${BAM:-}" ]]; then
  echo "Usage: sbatch $0 --sample_id <ID> --bam <bam_path>"
  exit 1
fi
[[ -f "$BAM" ]] || { echo "BAM not found: $BAM"; exit 1; }

# # ---------------- 1. 路径设置 ----------------
# module load singularity
# SIF=/vast/projects/quartet_rna_refdata/images/preseq_3.2.0.sif             # 镜像
# Load Conda
source /home/users/allstaff/chen.q/miniconda3/bin/activate /home/users/allstaff/chen.q/miniconda3/envs/preseq3.2

OUTDIR=/vast/projects/quartet_rna_refdata/analysis/preseq/${SAMPLE}
mkdir -p "$OUTDIR"

echo "[$(date)]  Start Preseq for ${SAMPLE}"
cd /vast/projects/quartet_rna_refdata/

# ---------------- 2. c_curve ----------------
preseq c_curve -B -s 10000 -o "${OUTDIR}/${SAMPLE}.c_curve.txt" "$BAM"

echo "[$(date)]  Done   Preseq c_curve for ${SAMPLE}"

# ---------------- 3. lc_extrap ----------------
preseq lc_extrap -B -s 1000000 -o "${OUTDIR}/${SAMPLE}.lc_extrap.txt" "$BAM"

echo "[$(date)]  Done   Preseq lc_extrap for ${SAMPLE}"