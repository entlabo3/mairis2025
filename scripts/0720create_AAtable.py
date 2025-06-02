#!/opt/conda/bin/python

import os
import pandas as pd
from Bio.Data import CodonTable
import shutil


def load_variant_table(file_path):
    """
    Load variant table
    """
    # Load without header
    df = pd.read_csv(file_path, sep="\t", header=None, skip_blank_lines=False)

    # Set first column to blank
    df.iloc[0, 0] = ""

    # Remove whitespace from string data
    df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    return df


def create_codon_dict():
    """
    Create codon to three-letter amino acid conversion dictionary from standard genetic code
    """
    # Dictionary for converting one-letter to three-letter amino acid codes
    one_to_three = {
        "A": "Ala",
        "C": "Cys",
        "D": "Asp",
        "E": "Glu",
        "F": "Phe",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "K": "Lys",
        "L": "Leu",
        "M": "Met",
        "N": "Asn",
        "P": "Pro",
        "Q": "Gln",
        "R": "Arg",
        "S": "Ser",
        "T": "Thr",
        "V": "Val",
        "W": "Trp",
        "Y": "Tyr",
        "*": "Stop",
    }

    # Get standard genetic code and create codon dictionary
    standard_table = CodonTable.standard_dna_table
    codon_dict = {}

    # Create dictionary for converting codons to three-letter amino acids (both lowercase and uppercase)
    for codon, aa in standard_table.forward_table.items():
        codon_dict[codon.lower()] = one_to_three[aa]
        codon_dict[codon.upper()] = one_to_three[aa]

    # Add stop codons (both lowercase and uppercase)
    for stop_codon in standard_table.stop_codons:
        codon_dict[stop_codon.lower()] = "Stop"
        codon_dict[stop_codon.upper()] = "Stop"

    return codon_dict


def process_value_vectorized(values, codon_dict):
    """
    Batch process each cell value
    Process based on value format:
    - "..." → output as is
    - Uppercase codon → convert to amino acid
    - Lowercase codon → "---" (synonymous substitution)
    """
    # Handle empty or missing values
    values = values.fillna("")

    # Return "..." as is
    mask_dots = values == "..."
    mask_upper = values.str.isupper() & (values.str.len() == 3)
    mask_lower = values.str.islower() & (values.str.len() == 3)

    # Process amino acid conversion and synonymous substitution
    result = values.copy()
    result[mask_upper] = values[mask_upper].map(lambda x: codon_dict.get(x, "Xxx"))
    result[mask_lower] = "---"
    result[~(mask_dots | mask_upper | mask_lower)] = "Xxx"

    return result


def process_variants_table(df, codon_dict):
    """
    Process variant table
    """
    # Keep first 7 rows as is, process row 8 and beyond
    result = df.copy()
    result.iloc[7:, 2:] = result.iloc[7:, 2:].apply(
        lambda col: process_value_vectorized(col, codon_dict)
    )
    return result


def main():
    """Main function"""
    input_folder = "./0710_variant_tables"
    input_file = os.path.join(input_folder, "variants_table_codons.tsv")
    output_folder = input_folder
    output_file = os.path.join(output_folder, "variants_table_aminoacids.tsv")

    try:
        # Create codon conversion dictionary
        codon_dict = create_codon_dict()

        # Load and process data
        df = load_variant_table(input_file)
        result_df = process_variants_table(df, codon_dict)

        # Save results (exclude index and column numbers)
        result_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="")
        print(f"Amino acid conversion table saved: {output_file}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 1

    return 0


if __name__ == "__main__":
    main()