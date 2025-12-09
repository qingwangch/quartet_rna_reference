#!/bin/bash
#SBATCH --job-name=suppa2_diffsplice              # Job name
#SBATCH --nodes=1                                 # Use 1 node
#SBATCH --ntasks=8                                # Allocate 8 CPU cores
#SBATCH --mem=16G                                 # Request 16GB of memory
#SBATCH --output=logs/suppa2_diffsplice_%j.log    # Output log file
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mail-type=FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=chen.q@wehi.edu.au            # Email notifications

# Exit the script on any error
set -eo pipefail

# Load Singularity environment and clear previous modules
module purge
module load singularity
# Load Conda
source /home/users/allstaff/chen.q/miniconda3/bin/activate /home/users/allstaff/chen.q/miniconda3/envs/suppa_env

# Define the location of the Singularity image
imgDIR="/vast/projects/quartet_rna_refdata/images"

# Set default values
method="empirical"
area=1000
lower_bound=0.05
ioe_file="/vast/projects/quartet_rna_refdata/reference/suppa2/gencode_v43.all_AS.ioe"

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --cond1_psi) cond1_psi="$2"; shift ;;
        --cond2_psi) cond2_psi="$2"; shift ;;
        --cond1_tpm) cond1_tpm="$2"; shift ;;
        --cond2_tpm) cond2_tpm="$2"; shift ;;
        --ioe_file) ioe_file="$2"; shift ;;
        --output_prefix) output_prefix="$2"; shift ;;
        --method) method="$2"; shift ;;
        --area) area="$2"; shift ;;
        --lower_bound) lower_bound="$2"; shift ;;
        *) echo "❌ Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure required parameters are provided
if [ -z "$cond1_psi" ] || [ -z "$cond2_psi" ] || [ -z "$cond1_tpm" ] || [ -z "$cond2_tpm" ] || [ -z "$output_prefix" ]; then
    echo "❌ Error: Missing required parameters!"
    echo "Usage: sbatch suppa2_diffsplice.slurm --cond1_psi <file> --cond2_psi <file> --cond1_tpm <file> --cond2_tpm <file> --output_prefix <prefix> [--method <empirical>] [--area <1000>] [--lower_bound <0.05>]"
    exit 1
fi

# Create output directory if needed
output_dir=$(dirname "$output_prefix")
mkdir -p "$output_dir"

# Start execution
date
echo "🚀 Starting Differential Splicing Analysis with SUPPA2"
echo "🔹 Conditions: $cond1_psi vs. $cond2_psi"
echo "🔹 Expression Files: $cond1_tpm & $cond2_tpm"
echo "🔹 Method: $method | Area: $area | Lower Bound: $lower_bound"
echo "🔹 Using IOE File: $ioe_file"

start_time=$(date +%s)

# Bind necessary directories
#BIND_PATHS="/vast/scratch/users/chen.q,/vast/projects/quartet_rna_refdata"

# Execute SUPPA2 differential splicing analysis
#singularity exec --bind $BIND_PATHS "$imgDIR/suppa_2.3.sif" \
python /vast/projects/quartet_rna_refdata/softwares/SUPPA/SUPPA-2.4/suppa.py diffSplice --method "$method" \
                        --input "$ioe_file" \
                        --psi "$cond1_psi" "$cond2_psi" \
                        --tpm "$cond1_tpm" "$cond2_tpm" \
                        --area "$area" \
                        --lower-bound "$lower_bound" \
                        -gc \
                        -o "$output_prefix"

# Record the end time and calculate execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

date
echo "✅ Differential Splicing Analysis Completed"
echo "⏳ Total execution time: ${execution_time} seconds"
