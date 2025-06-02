#!/opt/conda/bin/python

import os
import csv
from pathlib import Path
import re
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

def parse_duplicate_metrics(metrics_file):
    """
    Parse Picard metrics file

    Args:
        metrics_file (str): Path to metrics file

    Returns:
        dict: Dictionary containing metrics information
    """
    try:
        with open(metrics_file, 'r') as f:
            content = f.read()

        # Extract metrics section using regex
        metrics_pattern = r"## METRICS CLASS\s+.*?\n(.*?)\n(.*?)\n"
        matches = re.search(metrics_pattern, content, re.DOTALL)

        if not matches:
            print(f"Warning: No metrics data found in {metrics_file}")
            return {}

        # Get header and value lines
        header_line = matches.group(1).strip()
        data_line = matches.group(2).strip()

        # Split by tab
        headers = header_line.split('\t')
        values = data_line.split('\t')

        print(f"Processing file: {metrics_file}")
        print(f"Headers: {headers}")
        print(f"Values: {values}")

        # Store metrics in dictionary
        metrics = {}
        for header, value in zip(headers, values):
            if header == 'LIBRARY':
                continue

            try:
                if header == 'PERCENT_DUPLICATION':
                    # Maintain 5 decimal places for percentage values
                    metrics[header] = round(float(value) * 100, 5)
                elif 'ESTIMATED_LIBRARY_SIZE' in header:
                    # Convert estimated library size to integer
                    metrics[header] = int(float(value)) if value != '' else 0
                elif any(x in header for x in ['EXAMINED', 'DUPLICATES', 'UNMAPPED']):
                    # Convert other numeric values to integer
                    metrics[header] = int(float(value)) if value != '' else 0
                else:
                    metrics[header] = value
            except ValueError:
                print(f"Warning: Cannot convert value '{value}' for {header}")
                metrics[header] = 0

        print(f"Analysis results: {metrics}")
        return metrics

    except Exception as e:
        print(f"Error: Problem occurred while analyzing metrics file {metrics_file}: {str(e)}")
        return {}


def save_metrics_summary(metrics_dict, output_file):
    """
    Save metrics information as TSV file

    Args:
        metrics_dict (dict): Metrics information for each sample
        output_file (str): Output file path
    """
    headers = [
        'Sample',
        'UNPAIRED_READS_EXAMINED',
        'READ_PAIRS_EXAMINED',
        'UNMAPPED_READS',
        'UNPAIRED_READ_DUPLICATES',
        'READ_PAIR_DUPLICATES',
        'PERCENT_DUPLICATION',
        'ESTIMATED_LIBRARY_SIZE'
    ]

    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(headers)

            for sample, metrics in sorted(metrics_dict.items()):
                row = [sample]
                for header in headers[1:]:
                    value = metrics.get(header, 0)
                    if header == 'PERCENT_DUPLICATION':
                        row.append(f"{value:.5f}")
                    else:
                        row.append(str(value))
                writer.writerow(row)

        print(f"Metrics summary saved: {output_file}")

    except Exception as e:
        print(f"Error: Problem occurred while saving metrics summary: {str(e)}")


def main():
    """Main process"""
    metrics_dir = "./0130_picard_metrics"
    output_dir = "./0140_metrics"
    clean_output_directory(output_dir)
    output_file = os.path.join(output_dir, "duplicate_metrics_summary.tsv")

    # Collect and analyze metrics files
    metrics_dict = {}
    metrics_files = list(Path(metrics_dir).glob("*_duplicate_metrics.tsv"))

    print(f"Number of files to process: {len(metrics_files)}")

    for metrics_file in metrics_files:
        sample_name = metrics_file.stem.replace("_duplicate_metrics.tsv", "")
        print(f"\nStarting process for sample {sample_name}")

        metrics = parse_duplicate_metrics(str(metrics_file))
        if metrics:
            metrics_dict[sample_name] = metrics

    # Create summary file
    if metrics_dict:
        save_metrics_summary(metrics_dict, output_file)
        print(f"\nNumber of samples analyzed: {len(metrics_dict)}")
    else:
        print("\nWarning: No valid metrics found")


if __name__ == "__main__":
    main()