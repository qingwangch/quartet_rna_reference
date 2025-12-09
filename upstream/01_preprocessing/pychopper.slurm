#!/bin/bash
#SBATCH --job-name=pychopper_trim         # Job name
#SBATCH --ntasks=1                        # Run on a single CPU
#SBATCH --cpus-per-task=8                 # Number of CPU cores per task
#SBATCH --mem=32G                         # Job memory request
#SBATCH --time=24:00:00                   # Time limit hrs:min:sec
#SBATCH --output=pychopper_trim_%j.log    # Standard output and error log
#SBATCH --mail-type=FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au       # Where to send email notifications

# Enable error detection
set -eo pipefail

# Load Singularity
module load singularity/
imgDIR=/vast/projects/quartet_rna_refdata/images

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --workdir1) workdir1="$2"; shift ;;
        --workdir2) workdir2="$2"; shift ;;
        --sample) sample="$2"; shift ;;
        --fastq) fastq="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure all required arguments are provided
if [ -z "$workdir1" ] || [ -z "$workdir2" ] || [ -z "$sample" ] || [ -z "$fastq" ]; then
    echo "Usage: $0 --workdir1 <workdir1> --workdir2 <workdir2> --sample <sample> --fastq <fastq>"
    exit 1
fi

# Change to the working directory
cd /vast/projects/quartet_rna_refdata

date
echo "pychopper+trim_polya begins"

start_time=$(date +%s)

# Create directories if they don't exist
mkdir -p ${workdir1}
mkdir -p ${workdir2}

# Run pychopper
echo "pychopper begins"
singularity exec $imgDIR/pychopper_latest.sif pychopper -t 8 -m edlib -z 100 -u ${workdir1}/${sample}.unclassified.fastq -r ${workdir1}/${sample}.cDNA_classifier_report.pdf -S ${workdir1}/${sample}.cDNA_classifier_report.tsv ${fastq} ${workdir1}/${sample}.full_length.fastq
echo "pychopper ends"

# Convert full_length.fastq to full_length.fasta
awk '{if(NR%4==1){print $1" "$NF}else if(NR%4==2){print $0}}' ${workdir1}/${sample}.full_length.fastq | sed 's/@/>/g' > ${workdir1}/${sample}.full_length.fasta

# Run trim_polya
echo "trim_polya begins"
singularity exec $imgDIR/trim_isoseq_polya_latest.sif /opt/conda/bin/trim_isoseq_polyA -G -t 8 -i ${workdir1}/${sample}.full_length.fasta > ${workdir2}/${sample}.full_length.trim_polyA.fasta
echo "trim_polya ends"

# Filter sequences shorter than 50 bp
awk 'NR%2{id=$0;next}NR%2==0{fa[id]=$0}END{for(i in fa){if(length(fa[i])>=50){print i"\n"fa[i]}}}' ${workdir2}/${sample}.full_length.trim_polyA.fasta > ${workdir2}/${sample}.full_length.trim_polyA.filter.fasta

# Clean up
rm ${workdir1}/${sample}.full_length.fastq

end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "pychopper+trim_polya ends"

echo "Program execution time: $execution_time seconds"