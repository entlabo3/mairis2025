#!/opt/conda/bin/python

import pandas as pd
import os
import re
from collections import defaultdict
import shutil

def clean_output_directory(output_dir: str) -> None:
    """
    Clean up output directory

    Args:
        output_dir: Directory path to clean up
    """
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

def get_codon_display(sample_codon, ref_codon, sample_aa, ref_aa):
    """
    Function to determine the codon display format

    Args:
        sample_codon: Sample codon
        ref_codon: Reference codon
        sample_aa: Sample amino acid
        ref_aa: Reference amino acid

    Returns:
        str: Display codon information
    """
    if sample_codon == ref_codon:
        return "..."
    elif sample_aa == ref_aa:
        return sample_codon.lower()  # Synonymous substitution in lowercase
    else:
        return sample_codon.upper()  # Amino acid substitution in uppercase


def format_refaa_pos(value):
    """
    Function to convert MdAA_pos value to appropriate format
    """
    try:
        return str(int(value))
    except (ValueError, TypeError):
        return "-"


def create_position_data_dict(df, positions):
    """
    Function to create dictionary of position data
    """
    pos_data = {}
    for pos in positions:
        pos_rows = df[df["genome_pos1"] == pos]
        if not pos_rows.empty:
            pos_df = pos_rows.iloc[0]
            pos_data[pos] = {
                "genome_pos2": str(int(pos_df["genome_pos2"])),
                "genome_pos3": str(int(pos_df["genome_pos3"])),
                "md_aa": str(pos_df["Md_AA"]),
                "mdaa_pos": format_refaa_pos(pos_df["MdAA_pos"]),
                "ref_aa": str(pos_df["Ref_AA"]),
                "ref_codon": str(pos_df["Ref_codon"])
            }
    return pos_data


def create_sample_data_dict(df):
    """
    Function to create dictionary of sample data
    """
    sample_data = defaultdict(dict)
    for _, row in df.iterrows():
        if pd.notna(row["genome_pos1"]):
            pos = int(row["genome_pos1"])
            sample_data[(row["Sample"], pos)] = {
                "sample_codon": row["Sample_codon"],
                "sample_aa": row["Sample_AA"],
                "ref_aa": row["Ref_AA"],
                "ref_codon": row["Ref_codon"]
            }
    return sample_data


def create_codon_table(input_file, output_file):
    """
    Function to create codon table from variant data
    """
    # Load and preprocess data
    df = pd.read_csv(input_file, sep="\t")

    # Ensure numeric type
    df["genome_pos1"] = pd.to_numeric(df["genome_pos1"], errors="coerce")
    positions = sorted(df["genome_pos1"].dropna().unique().astype(int))

    # Pre-dictionary position data
    pos_data = create_position_data_dict(df, positions)

    # Pre-dictionary sample data
    sample_data = create_sample_data_dict(df)

    # Create header row
    headers = ["", "genome_pos1"] + [str(pos) for pos in positions]
    rows = [headers]

    # Create basic information rows
    basic_rows = [
        ["", "genome_pos2"] + [pos_data[pos]["genome_pos2"] for pos in positions],
        ["", "genome_pos3"] + [pos_data[pos]["genome_pos3"] for pos in positions],
        ["", "MdAA"] + [pos_data[pos]["md_aa"] for pos in positions],
        ["", "MdAA_pos"] + [pos_data[pos]["mdaa_pos"] for pos in positions],
        ["Reference", "AA"] + [pos_data[pos]["ref_aa"] for pos in positions],
        ["Reference", "Haplotype"] + [pos_data[pos]["ref_codon"] for pos in positions]
    ]
    rows.extend(basic_rows)

    # Process sample names
    samples = set()
    for sample in df["Sample"].unique():
        base_name = re.sub(r"_H[12]$", "", sample)
        samples.add(base_name)
    samples = sorted(list(samples))

    # Add codon information for each sample
    for base_sample in samples:
        for hap in ["H1", "H2"]:
            sample_name = f"{base_sample}_{hap}"
            sample_codons = []

            for pos in positions:
                key = (sample_name, pos)
                if key in sample_data:
                    data = sample_data[key]
                    display_codon = get_codon_display(
                        data["sample_codon"],
                        data["ref_codon"],
                        data["sample_aa"],
                        data["ref_aa"]
                    )
                else:
                    display_codon = "..."
                sample_codons.append(display_codon)

            rows.append([base_sample, hap] + sample_codons)

    # Convert to DataFrame and save
    output_df = pd.DataFrame(rows[1:], columns=rows[0])
    output_df.to_csv(output_file, sep="\t", index=False)


def main():
    """Main function"""
    input_folder = "./0610_list_variants"
    input_file = os.path.join(input_folder, "variants_with_refAA_pos.tsv")
    output_folder = "./0710_variant_tables"
    output_file = os.path.join(output_folder, "variants_table_codons.tsv")

    # Check output directory
    clean_output_directory(output_folder)

    create_codon_table(input_file, output_file)
    print(f"Created codon table: {output_file}")


if __name__ == "__main__":
    main()