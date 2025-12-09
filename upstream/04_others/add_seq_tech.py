# Define the mapping rules for seq-tech based on the sample_id prefix
def get_seq_tech(sample_id):
    if sample_id.startswith(("C_BGI", "P_ONT")):
        return "ont-cdna"
    elif sample_id.startswith("D_ONT"):
        return "ont-drna"
    elif sample_id.startswith("M_PAB"):
        return "pac-bio-hifi"
    else:
        return "unknown"  # Handle unexpected cases

# Process the input data
input_file = "fastq_s139.txt"  # Replace with your input filename
output_file = "fastq_s139_with_seq_tech.txt"

with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for line in f_in:
        if line.strip():  # Skip empty lines
            sample_id, path = line.strip().split(" ", 1)
            seq_tech = get_seq_tech(sample_id)
            f_out.write(f"{sample_id} {path} {seq_tech}\n")

print(f"Processed data saved to: {output_file}")
