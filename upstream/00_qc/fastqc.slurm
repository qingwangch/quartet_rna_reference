#!/bin/bash
###############################################################################
# Slurm job options
###############################################################################
#SBATCH --job-name=fastqc                    # 作业名称
#SBATCH --output=logs/fastqc_%j.log          # 标准输出和错误日志（合并）
#SBATCH --error=logs/fastqc_%j.log
#SBATCH --nodes=1                            # 节点数
#SBATCH --ntasks=1                           # 每个节点的任务数
#SBATCH --cpus-per-task=8                    # 每个任务的 CPU 核心数
#SBATCH --mem=8G                            # 内存
#SBATCH --time=24:00:00                      # 运行时限 (HH:MM:SS)
#SBATCH --mail-user=chen.q@wehi.edu.au       # 邮件通知地址
#SBATCH --mail-type=FAIL                     # 仅在失败时通知
###############################################################################

# 遇到错误立即退出，并让管道错误向上传递
set -eo pipefail

###############################################################################
# 软件环境
###############################################################################
module load singularity


###############################################################################
# 变量定义（可根据需要改为命令行参数）
###############################################################################
imgDIR=/vast/projects/quartet_rna_refdata/images
workdir=/vast/projects/quartet_rna_refdata
sample_id=${sample_id:-"SAMPLE"}   # 若未提前 export，可在提交脚本时 --export
read1=${read1:-"/path/to/read1.fastq.gz"}
read2=${read2:-"/path/to/read2.fastq.gz"}

###############################################################################
# 工作目录与日志
###############################################################################
cd "${workdir}"
echo "========== $(date)  fastqc begins =========="

start_time=$(date +%s)

result_dir=/vast/projects/quartet_rna_refdata/analysis/fastqc/${sample_id}
mkdir -p "${result_dir}"

###############################################################################
# 运行 FastQC（使用容器）
###############################################################################
BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

singularity exec -B "${BIND_PATHS}" \
    "${imgDIR}/fastqc_v0.12.0.sif" fastqc -t 8 -o "${result_dir}" "${read1}"

singularity exec -B "${BIND_PATHS}" \
    "${imgDIR}/fastqc_v0.12.0.sif" fastqc -t 8 -o "${result_dir}" "${read2}"

###############################################################################
# 结束与计时
###############################################################################
end_time=$(date +%s)
echo "========== $(date)  fastqc ends   =========="
echo "Execution time: $((end_time - start_time)) seconds"
