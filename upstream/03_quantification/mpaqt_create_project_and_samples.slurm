#!/bin/bash
#SBATCH --job-name=mpaqt_project_create
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=logs/mpaqt_create_project_%j.log
#SBATCH --error=logs/mpaqt_create_project_%j.log

set -eo pipefail
module load singularity

# ------------------------ Parse arguments ------------------------ #
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --index) index_path="$2"; shift ;;
    --project) project_name="$2"; shift ;;
    --samples) samples="$2"; shift ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
  shift
done

# Check required args
if [[ -z "$index_path" || -z "$project_name" || -z "$samples" ]]; then
  echo "Usage: sbatch $0 --index <index.mpaqt> --project <project_name> --samples \"sample1 sample2 ...\""
  exit 1
fi

# ------------------------ TMPDIR setup ------------------------ #
mpaqt_root=/vast/projects/quartet_rna_refdata/analysis/mpaqt
export TMPDIR="$mpaqt_root/tmpdir"   

# ------------------------ Project dir mapping ------------------------ #
PROJECT_DIR="/vast/projects/quartet_rna_refdata/analysis/mpaqt/projects"

imgDIR=/vast/projects/quartet_rna_refdata/images
cd /vast/projects/quartet_rna_refdata/

start_time=$(date +%s)
date
echo "Creating project $project_name with index $index_path"

# ------------------------ Run MPAQT: Create Project ------------------------ #
TMPDIR=$TMPDIR singularity exec \
  --bind "$mpaqt_root/tmpdir:/csglab/mpaqt/tmp" \
  --bind "$PROJECT_DIR:/csglab/mpaqt/projects" \
  "$imgDIR/mpaqt_0.3.1.sif" /csglab/mpaqt/mpaqt create project \
  --index "$index_path" "$project_name"

# ------------------------ Run MPAQT: Create Samples ------------------------ #
echo "Creating samples: $samples"

TMPDIR=$TMPDIR singularity exec \
  --bind "$mpaqt_root/tmpdir:/csglab/mpaqt/tmp" \
  --bind "$PROJECT_DIR:/csglab/mpaqt/projects" \
  "$imgDIR/mpaqt_0.3.1.sif" /csglab/mpaqt/mpaqt create sample \
  --project "$project_name" $samples

# ------------------------ Finish ------------------------ #
rm -rf "$TMPDIR"/*   # 仅删除 TMPDIR 下的所有文件和子目录，而不删除 TMPDIR 本身

end_time=$(date +%s)
execution_time=$((end_time - start_time))
date
echo "Project and sample creation complete. Time used: $execution_time seconds"
