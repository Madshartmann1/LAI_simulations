import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Modify FASTA file based on VCF")
parser.add_argument("-i", "--input", required=True, help="Input FASTA file path")
parser.add_argument("-v", "--vcf", required=True, help="VCF file path")
parser.add_argument("-o", "--output", required=True, help="Output FASTA file path")

# Parse arguments
args = parser.parse_args()

# File paths from arguments
fasta_path = args.input
vcf_path = args.vcf
output_path = args.output


# Read and store the entire FASTA sequence
with open(fasta_path, 'r') as fasta_file:
    # Skipping the header line
    header = fasta_file.readline()
    # Read the rest of the file and remove newlines
    fasta_sequence = ''.join(fasta_file.read().splitlines())

# Convert the sequence into a list for easy modification
fasta_list = list(fasta_sequence)

# Process the VCF file
with open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        if line.startswith('#'):
            continue  # Skip header lines
        parts = line.strip().split('\t')
        pos = int(parts[1]) - 1  # Convert to 0-based index
        ref = parts[3]  # Get the REF nucleotide

        # Replace the nucleotide in the FASTA sequence
        fasta_list[pos] = ref

# Join the list back into a string
modified_sequence = ''.join(fasta_list)

# Write the modified sequence to a new FASTA file
with open(output_path, 'w') as output_file:
    output_file.write(header)  # Write the original header
    output_file.write(modified_sequence)  # Write the modified sequence
    output_file.write('\n')  # Add a newline at the end
