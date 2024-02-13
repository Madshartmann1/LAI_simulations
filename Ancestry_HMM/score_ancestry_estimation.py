import sys
import re
import gzip
from datetime import datetime


def print_with_timestamp(message):
    print(f"{datetime.now()}: {message}")

def read_intervals(intervals_file):
    """
    Parses the intervals file, organizing data by individual and chromosome type.
    Adjusted to output a nested dictionary for direct access to intervals by individual and chromosome type.

    Args:
        intervals_file (str): Path to the intervals file.

    Returns:
        dict: A nested dictionary {individual: {chrom: [(start_pos, end_pos, barcode), ...], ...}, ...}.
    """

    print_with_timestamp(f"Reading intervals from {intervals_file}...")

    intervals_dict = {}
    current_individual = None

    # Adjust to handle gzipped intervals file
    open_func = gzip.open if intervals_file.endswith('.gz') else open

    with open_func(intervals_file, 'rt') as file:  # 'rt' mode for reading as text
        for line in file:
            line = line.strip()
            # Check for individual designation
            if line.startswith("Individual"):
                current_individual = line.split()[1]
                intervals_dict[current_individual] = {"M": [], "P": []}
            else:
                parts = line.split()
                if parts and parts[0] in ("M", "P"):
                    current_chrom = parts[0]
                    # Process interval data on the same line as the chromosome designation
                    if len(parts) > 1:
                        interval_data = ' '.join(parts[1:]).split(', ')
                        if len(interval_data) == 3:
                            start_pos, end_pos, barcode = interval_data
                            intervals_dict[current_individual][current_chrom].append((int(start_pos), int(end_pos), barcode))
                elif parts:  # Processing interval lines without M or P designation
                    interval_data = line.split(', ')
                    if len(interval_data) == 3:
                        start_pos, end_pos, barcode = interval_data
                        intervals_dict[current_individual][current_chrom].append((int(start_pos), int(end_pos), barcode))

    print_with_timestamp("Finished reading intervals.")
    return intervals_dict


def score_estimates(individual_truth_data, best_posterior_path):
    positive_score = 0
    total_sites = 0
    with open(best_posterior_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            position, estimated_ancestry = line.strip().split('\t')
            position = int(position)
            total_sites += 1
            estimated_parts = sorted(estimated_ancestry.split('/'))  # Sort to ignore order in comparison

            # Collect the true ancestries for this position across M and P
            true_ancestries = []
            for parent_type, intervals in individual_truth_data.items():
                for start, end, true_ancestry in intervals:
                    if start <= position <= end:
                        true_ancestry_main_part = true_ancestry.split('#')[0]
                        true_ancestries.append(true_ancestry_main_part)

            # Sort and compare
            if len(true_ancestries) >= len(estimated_parts):
                true_ancestries_sorted = sorted(set(true_ancestries))  # Deduplicate and sort true ancestries
                if all(est_part in true_ancestries_sorted for est_part in estimated_parts):
                    positive_score += 1

    # Calculate TPCR
    tpcr = positive_score / total_sites if total_sites > 0 else 0
    return positive_score, total_sites, tpcr



def main(truth_file, best_posterior_file, output_file):
    all_truth_data = read_intervals(truth_file)
    individual_id = best_posterior_file.split('/')[-1].split('_')[0]  # Extract individual ID from filename
    individual_truth_data = all_truth_data.get(individual_id, {"M": [], "P": []})

    positive_score, total_sites, tpcr = score_estimates(individual_truth_data, best_posterior_file)
    
    with open(output_file, 'w') as out_file:
        out_file.write(f"Score for Individual {individual_id}: {positive_score}\n")
        out_file.write(f"Total Sites Evaluated: {total_sites}\n")
        out_file.write(f"TPCR (True Positive Call Rate): {tpcr:.4f}\n")

if __name__ == "__main__":
    truth_file_path = sys.argv[1]
    best_posterior_path = sys.argv[2]
    output_path = sys.argv[3]
    main(truth_file_path, best_posterior_path, output_path)
