import sys
import re
import gzip
from datetime import datetime
import numpy as np


def print_with_timestamp(message):
    print(f"{datetime.now()}: {message}")


def ancestry_key(ancestry):
    mapping = {'CEU': 'E', 'YRI': 'Y', 'CHB': 'C', 'KAR': 'K'}
    # Split the ancestry string, map to single characters, sort, and join
    parts = ancestry.split('-')
    ancestry_codes = [mapping.get(part, '') for part in parts]
    # If dealing with single ancestries, consider how they should be represented
    if len(ancestry_codes) == 1:
        ancestry_codes *= 2  # Repeat the code to match expected format
    return ''.join(sorted(ancestry_codes))


def initialize_ancestry_matrix_and_index():
    ancestries = ['E', 'Y', 'C', 'K']
    # Generate all unique combinations for single chromosomes, considering reversed pairs as the same
    combinations = sorted(set(["".join(sorted([a, b])) for a in ancestries for b in ancestries]))
    matrix_index = {combo: i for i, combo in enumerate(combinations)}  # Map combinations to indices
    matrix = np.zeros((len(combinations), len(combinations)))  # Initialize a matrix of zeros
    return matrix, matrix_index

def update_ancestry_matrix(matrix, matrix_index, true_ancestry, estimated_ancestry):
    true_key = ancestry_key(true_ancestry)
    est_key = ancestry_key(estimated_ancestry)

    try:
        i = matrix_index[true_key]
        j = matrix_index[est_key]
        matrix[i, j] += 1
    except KeyError as e:
        print(f"KeyError encountered with key: {e}. True key: '{true_key}', Estimated key: '{est_key}'")
        # Handle the error or debug further

def normalize_ancestry_matrix(matrix):
    total_sum = np.sum(matrix)
    if total_sum > 0:
        normalized_matrix = matrix / total_sum
    else:
        normalized_matrix = matrix  # If the matrix sum is 0, return the original matrix to avoid division by zero
    return normalized_matrix



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
    ancestry_matrix, matrix_index = initialize_ancestry_matrix_and_index()
    
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
            true_ancestries_sorted = sorted(set(true_ancestries))

            # Update score and matrix
            if len(true_ancestries) >= len(estimated_parts):
                true_str = '-'.join(true_ancestries_sorted)
                est_str = '-'.join(estimated_parts)
                update_ancestry_matrix(ancestry_matrix, matrix_index, true_str, est_str)
                if true_str == est_str:
                    positive_score += 1
                    
    # Calculate TPCR
    tpcr = positive_score / total_sites if total_sites > 0 else 0
    raw_matrix = np.copy(ancestry_matrix)  # Assuming 'matrix' is your populated matrix variable
    # Scale the matrix to 1 using total_sites
    # Normalize the ancestry matrix
    normalized_ancestry_matrix = normalize_ancestry_matrix(ancestry_matrix)

    return positive_score, total_sites, normalized_ancestry_matrix, matrix_index, tpcr, raw_matrix

def print_ancestry_matrix(matrix, matrix_index, out_file):
    ancestries = sorted(matrix_index.keys(), key=lambda x: matrix_index[x])
    header = ["True\\Estimated"] + ancestries  # Clarifying the matrix orientation

    # Write the column headers
    out_file.write('\t'.join(header) + '\n')
    
    # Iterate over the matrix to write each row
    for true_ancestry in ancestries:
        row = [true_ancestry] + [f"{matrix[matrix_index[true_ancestry], matrix_index[est_ancestry]]:.4f}" for est_ancestry in ancestries]
        out_file.write('\t'.join(row) + '\n')

def format_matrix_to_latex(matrix, title, matrix_type='raw', explanation="Rows represent true ancestry, and columns represent estimated ancestries."):
    """
    Generates a LaTeX table representation of the given matrix, dynamically formatting based on matrix type,
    and includes an explanation of the rows and columns.

    Args:
    - matrix: A 2D numpy array representing the matrix to be formatted.
    - title: A string representing the title of the matrix.
    - matrix_type: A string indicating the type of matrix ('raw' or 'normalized').
    - explanation: A string containing an explanation of what the rows and columns represent.
    """
    # Define ancestries and generate combinations
    ancestries = ['E', 'Y', 'C', 'K']
    combinations = sorted(set(["".join(sorted([a, b])) for a in ancestries for b in ancestries]))

    # Use combinations as both row and column names
    row_names = combinations
    col_names = combinations

    latex_str = f"\\begin{{table}}[h]\n\\centering\n\\caption{{{title}}}\n\\begin{{tabular}}{{{'|'.join(['c' for _ in range(len(col_names) + 1)])}}}\n\\hline\n"
     # Column headers
    latex_str += " & " + " & ".join(col_names) + " \\\\\n\\hline\n"
    # Rows with alternating background color
    for i, row in enumerate(matrix):
        if i % 2 == 0:  # Apply grey background to every other row
            latex_str += "\\rowcolor{Gray} "
        if matrix_type == 'raw':
            formatted_row = [f"{int(val)}" for val in row]
        elif matrix_type == 'normalized':
            formatted_row = [f"{val:.4f}".rstrip('0').rstrip('.') if val != 0 else '0' for val in row]
        latex_str += f"{row_names[i]} & " + " & ".join(formatted_row) + " \\\\\n"
    latex_str += "\\hline\n\\end{tabular}\n"
    if explanation:
        latex_str += f"\\caption*{{{explanation}}}\n"  # Add explanation as a caption note
    latex_str += "\\end{table}"
    
    return latex_str






def main(truth_file, best_posterior_file, output_file):
    all_truth_data = read_intervals(truth_file)
    individual_id = best_posterior_file.split('/')[-1].split('_')[0]  # Extract individual ID from filename
    individual_truth_data = all_truth_data.get(individual_id, {"M": [], "P": []})

    positive_score, total_sites, ancestry_matrix_scaled, matrix_index, tpcr, raw_matrix = score_estimates(individual_truth_data, best_posterior_file)
    
    with open(output_file, 'w') as out_file:
        out_file.write(f"Score for Individual {individual_id}: {positive_score}\n")
        out_file.write(f"Total Sites Evaluated: {total_sites}\n")
        out_file.write(f"TPCR (True Positive Call Rate): {tpcr:.4f}\n")
        print_ancestry_matrix(ancestry_matrix_scaled, matrix_index, out_file)
        out_file.write("\n")
        print_ancestry_matrix(raw_matrix, matrix_index, out_file)

    # Generate LaTeX for raw matrix (integers)
    latex_raw = format_matrix_to_latex(raw_matrix, "Raw Matrix", matrix_type='raw')
    # Generate LaTeX for normalized matrix (4 decimal places)
    latex_normalized = format_matrix_to_latex(ancestry_matrix_scaled, "Normalized Matrix", matrix_type='normalized')

    # Assuming output_file is something like "scoring_files/individualID_score.txt"
    latex_table_file = output_file.replace('_score.txt', '_score.table')

    # Now, latex_table_file would be something like "scoring_files/individualID_score.table"
    with open(latex_table_file, "w") as file:
        file.write(latex_raw + "\n\n" + latex_normalized)



if __name__ == "__main__":
    truth_file_path = sys.argv[1]
    best_posterior_path = sys.argv[2]
    output_path = sys.argv[3]
    main(truth_file_path, best_posterior_path, output_path)
