#!/opt/conda/bin/python

import os
import pandas as pd
import logging
from pathlib import Path


def setup_logger():
    """Set up the logger"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def parse_gff3(gff3_file: str) -> dict:
    """
    Parse GFF3 file and create mapping between genome positions and exon names

    Args:
        gff3_file: Path to the GFF3 file

    Returns:
        dict: Dictionary with genome positions as keys and exon names as values
    """
    exon_map = {}
    try:
        with open(gff3_file, 'r') as file:
            for line in file:
                if line.startswith("#"):
                    continue

                columns = line.strip().split("\t")
                if len(columns) < 9 or columns[2] != "CDS":
                    continue

                try:
                    start = int(columns[3])
                    end = int(columns[4])
                except ValueError:
                    continue

                # Extract exon name from attributes
                attributes = dict(
                    attr.split('=', 1) for attr in columns[8].split(";")
                    if '=' in attr
                )
                exon_name = attributes.get('Name')

                if exon_name:
                    for position in range(start, end + 1):
                        exon_map[position] = exon_name

        return exon_map

    except Exception as e:
        logging.error(f"Error parsing GFF3 file: {e}")
        raise


def process_file(input_file: str, output_file: str, exon_map: dict):
    """
    Process variant table file and add exon names

    Args:
        input_file: Path to input file
        output_file: Path to output file
        exon_map: Mapping between genome positions and exon names
    """
    try:
        # Load data
        data = pd.read_csv(input_file, sep="\t", header=None)

        # Generate exon name rows
        exon_rows = []
        for row_index in range(3):  # For top 3 rows
            exon_row = []
            # First two columns are empty
            exon_row.extend([""] * 2)

            # Process columns from the third onwards
            for col in range(2, data.shape[1]):
                cell_value = data.iloc[row_index, col]

                if pd.isna(cell_value) or not str(cell_value).isdigit():
                    exon_row.append("NA")
                else:
                    position = int(cell_value)
                    exon_row.append(exon_map.get(position, "NA"))

            exon_rows.append(exon_row)

        # Combine exon name rows with original data
        exon_rows_df = pd.DataFrame(exon_rows)
        processed_data = pd.concat([exon_rows_df, data], ignore_index=True)

        # Save results
        processed_data.to_csv(output_file, sep="\t", header=False, index=False)

    except Exception as e:
        logging.error(f"Error processing file: {e}")
        raise


def main():
    """Main processing function"""
    logger = setup_logger()

    # Set paths
    input_folder = "0710_variant_tables"
    gff_folder = "./ref"
    output_folder = "0710_variant_tables"

    files = {
        'aminoacids': {
            'input': os.path.join(input_folder, "variants_table_aminoacids.tsv"),
            'output': os.path.join(output_folder, "variants_table_aminoacids_with_exon_name.tsv")
        },
        'codons': {
            'input': os.path.join(input_folder, "variants_table_codons.tsv"),
            'output': os.path.join(output_folder, "variants_table_codons_with_exon_name.tsv")
        }
    }

    gff3_file = os.path.join(gff_folder, "ref_trimmed.gff3")

    try:
        # Check output directory
        Path(output_folder).mkdir(parents=True, exist_ok=True)

        # Check input files
        for file_type, paths in files.items():
            if not os.path.exists(paths['input']):
                raise FileNotFoundError(f"Input file for {file_type} not found: {paths['input']}")

        logger.info("Parsing GFF3 file...")
        exon_map = parse_gff3(gff3_file)
        logger.info(f"Created mapping for {len(exon_map)} genome positions")

        # Process each file
        for file_type, paths in files.items():
            logger.info(f"Processing {file_type} file...")
            process_file(paths['input'], paths['output'], exon_map)
            logger.info(f"Completed processing {file_type} file")

        logger.info("All processing completed successfully")
        logger.info("Created files with added exon names successfully!")

    except Exception as e:
        logger.error(f"An error occurred during processing: {e}")
        raise


if __name__ == "__main__":
    main()