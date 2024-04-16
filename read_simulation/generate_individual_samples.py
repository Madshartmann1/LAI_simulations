import argparse
from pysam import VariantFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def apply_variants_to_sequence(ref_sequence, variants, sample_idx, vcf):
    """
    Apply variants to the reference sequence to reconstruct the haplotype,
    selecting the correct allele based on the sample index (even: maternal, odd: paternal).
    """
    seq = list(ref_sequence)  # Convert to list for mutability
    individual_index = int(sample_idx[1:]) // 2  # Calculate individual's index based on sample_idx
    haplotype_suffix = "_M" if int(sample_idx[1:]) % 2 == 0 else "_P"  # Determine haplotype (Maternal or Paternal)

    # Assuming sample names are in the order they appear in the VCF header
    vcf_sample_names = list(vcf.header.samples)
    actual_sample_name = vcf_sample_names[individual_index]  # Get the actual sample name from VCF

    for variant in variants:
        # Assuming 'GT' is accessible via the sample name and that pysam returns a tuple
        sample_gt = variant.samples[actual_sample_name]['GT']
        haplotype_idx = 0 if haplotype_suffix == "_M" else 1  # 0 for maternal, 1 for paternal allele
        alleles = [variant.ref] + list(variant.alts)
        allele_idx = sample_gt[haplotype_idx]
        allele = alleles[allele_idx]
        pos = variant.pos - 1  # Adjust for 0-based indexing
        
        seq[pos] = allele if allele in alleles else ref_sequence[pos]

    return ''.join(seq)



def main(vcf_path, ref_path, sample_idx, output_path):
    # Load the reference sequence
    ref_seq_records = list(SeqIO.parse(ref_path, "fasta"))
    ref_seq = str(ref_seq_records[0].seq)
    
    # Open VCF file and fetch sample names
    vcf = VariantFile(vcf_path)
    sample_names = list(vcf.header.samples)
    
    # Calculate sample index (assuming n0 and n1 are the same individual, etc.)
    individual_idx = int(sample_idx[1:]) // 2  # Convert 'nX' to integer and halve it
    if individual_idx < len(sample_names):  # Check if it's within the phased sample range
        sample_name = sample_names[individual_idx]
    else:
        # Handle numerical samples differently if necessary
        print(f"Sample index {sample_idx} out of range for phased individuals.")
        return
    
    # Extract variants for the specified sample
    variants = [rec for rec in vcf.fetch() if sample_name in rec.samples]

    # Determine the haplotype suffix based on the sample index
    haplotype_suffix = "_M" if int(sample_idx[1:]) % 2 == 0 else "_P"

    # Adjust the sample name to include the haplotype suffix
    haplotype_sample_name = sample_name + haplotype_suffix
    
    # Assume function to apply variants to sequence
    haplotype_seq = apply_variants_to_sequence(ref_seq, variants, sample_idx, vcf)
    
    # Output handling 
    # Write haplotype sequence to output FASTA file
    haplotype_record = SeqRecord(Seq(haplotype_seq), id=haplotype_sample_name, description="")
    with open(output_path, "w") as output_file:
        SeqIO.write(haplotype_record, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate haplotype sequences from a phased VCF.")
    parser.add_argument("--vcf", required=True, help="Path to the VCF file.")
    parser.add_argument("--ref", required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("--sample", required=True, help="Sample index as nX.")
    parser.add_argument("--output", required=True, help="Output FASTA file path.")
    
    args = parser.parse_args()
    main(args.vcf, args.ref, args.sample, args.output)
