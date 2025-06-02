#!/opt/conda/bin/python

import pandas as pd
from typing import List
import os
from difflib import SequenceMatcher
import shutil


def find_common_prefix(strings: List[str]) -> str:
    """Extract the longest common part from multiple strings"""
    if not strings:
        return ""

    # Find common parts using SequenceMatcher
    matcher = SequenceMatcher(None, strings[0], strings[1])
    match = matcher.find_longest_match(0, len(strings[0]), 0, len(strings[1]))
    return strings[0][match.a:match.a + match.size]


def merge_transcripts(df: pd.DataFrame) -> pd.DataFrame:
    """Merge rows with non-Common Transcript values that have identical other columns"""

    # Group by columns other than Transcript
    cols_for_grouping = [col for col in df.columns if col != 'Transcript']

    # List to store results
    merged_rows = []

    # Process each group
    for _, group in df.groupby(cols_for_grouping):
        if len(group) == 1:  # For single rows, add as is
            merged_rows.append(group.iloc[0])
        else:
            # Extract rows with non-Common Transcript
            non_common = group[group['Transcript'] != 'Common']
            if len(non_common) > 1:  # Multiple non-Common rows
                # Extract common part from Transcript names
                common_transcript = find_common_prefix(non_common['Transcript'].tolist())
                if common_transcript:  # If common part found
                    new_row = non_common.iloc[0].copy()
                    new_row['Transcript'] = common_transcript
                    merged_rows.append(new_row)
            elif len(non_common) == 1:  # Single non-Common row
                merged_rows.append(non_common.iloc[0])
            else:  # All rows are Common
                merged_rows.append(group.iloc[0])

    return pd.DataFrame(merged_rows)


def main():
    """Main process"""
    # Set input/output directory
    input_dir = "./0610_list_variants"

    # Set input/output filenames
    input_filename = "variants_snp_mnp_processed_with_genome_pos.tsv"
    output_filename = "variants_snp_mnp_processed_with_genome_pos_merged.tsv"

    # Create full paths
    input_path = os.path.join(input_dir, input_filename)
    output_path = os.path.join(input_dir, output_filename)

    try:
        # Read TSV file
        df = pd.read_csv(input_path, sep='\t')

        # Process transcript merging
        merged_df = merge_transcripts(df)

        # Save results
        merged_df.to_csv(output_path, sep='\t', index=False)
        print(f"Processing completed. Output file: {output_path}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        raise


if __name__ == "__main__":
    main()