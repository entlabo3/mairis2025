#!/opt/conda/bin/python

import pandas as pd
import os
from typing import Dict, Set


def get_transcripts_from_gff3(gff3_file: str) -> Set[str]:
    """
    Function to get transcript IDs from GFF3 file

    Args:
        gff3_file: Path to GFF3 file

    Returns:
        Set[str]: Set of transcript IDs
    """
    transcripts = set()

    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if fields[2] == 'transcript':
                attributes = dict(attr.split('=', 1) for attr in fields[8].split(';'))
                if 'ID' in attributes:
                    transcript_id = attributes['ID']
                    transcripts.add(transcript_id)

    if not transcripts:
        raise ValueError("Could not extract transcript IDs from GFF3 file")

    return transcripts


def load_variants_data(input_file: str) -> pd.DataFrame:
    """
    Function to load variant data

    Args:
        input_file: Input file path

    Returns:
        pd.DataFrame: Loaded DataFrame
    """
    return pd.read_csv(input_file, sep='\t')


def split_variants_by_type(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Function to split DataFrame by variant type

    Args:
        df: Input DataFrame

    Returns:
        Dict[str, pd.DataFrame]: Split DataFrames
    """
    snp_mnp_mask = df['Variant_type'].isin(['SNP', 'MNP'])
    ins_del_mask = df['Variant_type'].isin(['INS', 'DEL'])

    return {
        'snp_mnp': df[snp_mnp_mask].copy(),
        'ins_del': df[ins_del_mask].copy()
    }


def process_variants(df: pd.DataFrame, all_transcripts: Set[str]) -> pd.DataFrame:
    """
    Function to integrate only variants common to all transcripts

    Args:
        df: Input DataFrame
        all_transcripts: Set of all transcript IDs defined in GFF3 file

    Returns:
        pd.DataFrame: Processed DataFrame
    """
    if df.empty:
        return df

    # Use columns other than Transcript as grouping keys
    group_columns = [col for col in df.columns if col != 'Transcript']

    # List to store results
    result_rows = []

    # Process each group
    for _, group in df.groupby(group_columns):
        # Get transcripts in the group
        group_transcripts = set(group['Transcript'])

        # Only integrate if all transcripts exist
        if group_transcripts == all_transcripts:
            # Get representative row and change Transcript to 'Common'
            common_row = group.iloc[0].to_dict()
            common_row['Transcript'] = 'Common'
            result_rows.append(common_row)
        else:
            # Add individual rows if not all transcripts are included
            for _, row in group.iterrows():
                result_rows.append(row.to_dict())

    if not result_rows:
        # Return empty DataFrame if no results
        return pd.DataFrame(columns=df.columns)

    # Convert results to DataFrame
    result_df = pd.DataFrame(result_rows)

    # Keep same column order as original DataFrame
    result_df = result_df[df.columns]

    return result_df


def save_variants_data(df: pd.DataFrame, output_file: str) -> None:
    """
    Function to save processed data

    Args:
        df: DataFrame to save
        output_file: Output file path
    """
    if df.empty:
        print(f"Warning: {output_file} will not be created as DataFrame is empty")
        return

    df.to_csv(output_file, sep='\t', index=False)
    print(f"File saved: {output_file}")


def main():
    # Set input/output paths
    input_dir = "./0610_list_variants"
    input_file = os.path.join(input_dir, "variants.tsv")
    gff3_file = "./ref/ref_trimmed.gff3"

    # Set output filenames
    output_files = {
        'snp_mnp': os.path.join(input_dir, "variants_snp_mnp_processed.tsv"),
        'ins_del': os.path.join(input_dir, "variants_ins_del_processed.tsv")
    }

    # Check input file existence
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    if not os.path.exists(gff3_file):
        raise FileNotFoundError(f"GFF3 file not found: {gff3_file}")

    # Get transcript information from GFF3 file
    all_transcripts = get_transcripts_from_gff3(gff3_file)
    print(f"Detected transcript IDs from GFF3 file ({len(all_transcripts)} items):")
    for transcript in sorted(all_transcripts):
        print(f"  - {transcript}")

    # Load data
    df = load_variants_data(input_file)

    # Split by variant type
    split_dfs = split_variants_by_type(df)

    # Process and save each type of data
    for variant_type, variant_df in split_dfs.items():
        print(f"\nProcessing {variant_type} type variants...")

        # Process data
        processed_df = process_variants(variant_df, all_transcripts)

        # Display statistics
        n_total = len(variant_df)
        n_common = len(processed_df[processed_df['Transcript'] == 'Common']) if not processed_df.empty else 0
        print(f"Number of rows before processing: {n_total}")
        print(f"Number of Common rows: {n_common}")
        print(f"Number of rows after processing: {len(processed_df)}")

        # Save results
        output_file = output_files[variant_type]
        save_variants_data(processed_df, output_file)

    print("\nAll processing completed")


if __name__ == "__main__":
    main()