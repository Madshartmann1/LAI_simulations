#!/bin/bash

# Input and output files
large_vcf="/home/projects/MAAG/msprime_deme/msprime/results/chrom1_500indi.vcf.gz"
output_vcf="/home/projects/MAAG/msprime_deme/msprime/results/chrom1_100indi.vcf"


# Function to select samples based on the segment range
select_samples_from_segment() {
    bcftools query -l "$large_vcf" | grep "^tsk_" | awk -F'_' -v start="$1" -v end="$2" '$2 >= start && $2 <= end' | shuf -n 100 | tr '\n' ',' | sed 's/,$//'
}
# Select 100 samples from each of the four populations
samples_pop_0_499=$(select_samples_from_segment 0 499)
samples_pop_500_999=$(select_samples_from_segment 500 999)
samples_pop_1000_1499=$(select_samples_from_segment 1000 1499)
samples_pop_1500_1999=$(select_samples_from_segment 1500 1999)


# Combine all selected samples
samples="$samples_pop_0_499,$samples_pop_500_999,$samples_pop_1000_1499,$samples_pop_1500_1999"

# Extract the selected samples and write them to a new VCF file
bcftools view -s "$samples" "$large_vcf" -Ov -o "$output_vcf"


