#!/bin/bash

# Assign input arguments to variables for better readability
large_vcf="$1"
small_vcf="$2"
output_merged_vcf="$3"

# Define function to select samples based on ranges for large VCF
select_samples_from_large_vcf() {
    bcftools query -l "$large_vcf" | grep 'tsk_' | awk -F'_' -v start="$1" -v end="$2" '$2 >= start && $2 <= end' | shuf -n 5 | tr '\n' ',' | sed 's/,$//'
}

# Select samples from each specified range
samples_0_499=$(select_samples_from_large_vcf 0 499)
samples_500_999=$(select_samples_from_large_vcf 500 999)
samples_1000_1499=$(select_samples_from_large_vcf 1000 1499)
samples_1499_1999=$(select_samples_from_large_vcf 1499 1999)


# Select random samples from small VCF
samples_small_vcf=$(bcftools query -l "$small_vcf" | shuf -n 5 | tr '\n' ',' | sed 's/,$//')


# Extract and merge the specified samples into a single VCF
temp_large_vcf_selected="large_vcf_selected.vcf.gz"
temp_small_vcf_selected="small_vcf_selected.vcf.gz"

bcftools view -s "$samples_0_499,$samples_500_999,$samples_1000_1499,$samples_1499_1999" "$large_vcf" -Oz -o "$temp_large_vcf_selected"
bcftools view -s "$samples_small_vcf" "$small_vcf" -Oz -o "$temp_small_vcf_selected"

# Index the compressed VCF files
tabix -p vcf "$temp_large_vcf_selected"
tabix -p vcf "$temp_small_vcf_selected"

# Merge the VCFs
bcftools merge "$temp_large_vcf_selected" "$temp_small_vcf_selected" -Ov -o "$output_merged_vcf"

# Clean up intermediate files
rm "$temp_large_vcf_selected" "$temp_small_vcf_selected"
rm "${temp_large_vcf_selected}.tbi" "${temp_small_vcf_selected}.tbi"
