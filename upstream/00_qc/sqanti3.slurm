#!/bin/bash
#SBATCH --job-name=sqanti3              # Job name
#SBATCH --ntasks=1                      # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem=32G                       # Job memory request
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=sqanti3_%j.log         # Standard output and error log

# Enable error detection
set -eo pipefail

# Load Singularity
module load singularity
# Define the image directory
imgDIR=/vast/projects/quartet_rna_refdata/images

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_id) sample_id="$2"; shift ;;
        --sample_gtf) sample_gtf="$2"; shift ;;
        --workdir1) workdir1="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure required arguments are provided
if [ -z "$sample_id" ] || [ -z "$sample_gtf" ] || [ -z "$workdir1" ]; then
    echo "Usage: $0 --sample_id <sample_id> --sample_gtf <sample_gtf> --workdir1 <workdir1>"
    exit 1
fi

# Change to the working directory
cd /vast/projects/quartet_rna_refdata/

date
echo "sqanti3 begins"

start_time=$(date +%s)

# Check if the output directory exists; create if it doesn't
if [ ! -d "$workdir1" ]; then
    mkdir -p "$workdir1"
    echo "Directory created: $workdir1"
else
    echo "Directory already exists: $workdir1"
fi

# Define the paths to the reference files
ref_gtf=/vast/projects/quartet_rna_refdata/ref_transcriptome/gencode.v43.chr_patch_hapl_scaff.annotation.gtf
ref_genome=/vast/projects/quartet_rna_refdata/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# Run SQANTI3 QC
awk '$7 != "." || NR < 3' ${sample_gtf} > ${workdir1}/${sample_id}_corrected.gtf

singularity exec $imgDIR/sqanti3_v5.2.sif bash -c "/opt2/sqanti3/5.2.0/SQANTI3-5.2/sqanti3_qc.py -t 16 -n 2 --aligner_choice minimap2 --genename --force_id_ignore --isoAnnotLite --saturation \
                                                   --report both \
                                                   --CAGE_peak /vast/projects/quartet_rna_refdata/ref_data/iso-seq/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
                                                   --polyA_motif_list /vast/projects/quartet_rna_refdata/ref_data/iso-seq/polyA.list.txt \
                                                   -c /vast/projects/quartet_rna_refdata/ref_data/iso-seq/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv \
                                                   -d ${workdir1} -o ${sample_id}.SQ3 ${workdir1}/${sample_id}_corrected.gtf ${ref_gtf} ${ref_genome}"

# Run SQANTI3 Filter (Rule-based filtering)
mkdir -p ${workdir1}/filter_rules
singularity exec $imgDIR/sqanti3_v5.2.sif bash -c "/opt2/sqanti3/5.2.0/SQANTI3-5.2/sqanti3_filter.py rules ${workdir1}/${sample_id}.SQ3_classification.txt \
                                                   -j /vast/projects/quartet_rna_refdata/slurm/sqanti3/filtering.json \
                                                   --isoforms ${workdir1}/${sample_id}.SQ3_corrected.fasta \
                                                   --gtf ${workdir1}/${sample_id}.SQ3_corrected.gtf\
                                                   -o ${sample_id}.SQ3.rules \
                                                   -d ${workdir1}/filter_rules"

# Run SQANTI3 Filter (Machine learning-based filtering)
mkdir -p ${workdir1}/filter_ml
singularity exec $imgDIR/sqanti3_v5.2.sif bash -c "/opt2/sqanti3/5.2.0/SQANTI3-5.2/sqanti3_filter.py ml ${workdir1}/filter_rules/${sample_id}.SQ3.rules_RulesFilter_result_classification.txt \
                                                   --isoforms ${workdir1}/filter_rules/${sample_id}.SQ3.rules.filtered.fasta \
                                                   --gtf ${workdir1}/filter_rules/${sample_id}.SQ3.rules.filtered.gtf \
                                                   --force_fsm_in TRUE \
                                                   -o ${sample_id}.SQ3.ml \
                                                   -d ${workdir1}/filter_ml"

# Run SQANTI3 Rescue (Rule-based rescue)
mkdir -p ${workdir1}/rescue_rules
singularity exec $imgDIR/sqanti3_v5.2.sif bash -c "sqanti3_rescue.py rules \
                                                --isoforms ${workdir1}/${sample_id}.SQ3_corrected.fasta \
                                                --gtf ${workdir1}/filter_ml/${sample_id}.SQ3.ml.filtered.gtf \
                                                -f ${ref_genome} \
                                                -g ${ref_gtf} \
                                                -k ${workdir1}/${sample_id}.SQ3_classification.txt \
                                                --mode automatic \
                                                -e all \
                                                -o ${sample_id}.SQ3.rescue_rules \
                                                -d ${workdir1}/rescue_rules \
                                                -j /vast/projects/quartet_rna_refdata/slurm/sqanti3/filtering.json \
                                                ${workdir1}/filter_ml/${sample_id}.SQ3.ml_MLresult_classification.txt"

# Run SQANTI3 Rescue (ML-based rescue)
mkdir -p ${workdir1}/rescue_ml
singularity exec $imgDIR/sqanti3_v5.2.sif bash -c "sqanti3_rescue.py ml \
                                                --isoforms ${workdir1}/${sample_id}.SQ3_corrected.fasta \
                                                --gtf ${workdir1}/filter_ml/${sample_id}.SQ3.ml.filtered.gtf \
                                                -f ${ref_genome} \
                                                -g ${ref_gtf} \
                                                -k ${workdir1}/${sample_id}.SQ3_classification.txt \
                                                --mode automatic \
                                                -e all \
                                                -o ${sample_id}.SQ3.rescue_ml \
                                                -d ${workdir1}/rescue_ml \
                                                -r ${workdir1}/filter_ml/randomforest.RData \
                                                -j 0.7 \
                                                ${workdir1}/filter_ml/${sample_id}.SQ3.ml_MLresult_classification.txt"



end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "sqanti3 ends"

echo "Program execution time: $execution_time seconds"
