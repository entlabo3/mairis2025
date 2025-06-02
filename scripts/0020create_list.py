#!/opt/conda/bin/python
import os
import re
from collections import defaultdict


def clean_name(name: str, pattern_to_remove) -> str:
    """
    Function to remove pattern_to_remove from sample names and file names

    Args:
        name (str): Original name

    Returns:
        str: Cleaned name
    """
    return re.sub(pattern_to_remove, '', name)


def get_fastq_files(data_dir):
    """Get FASTQ files from data directory and group them by R1 and R2"""
    files = defaultdict(dict)
    r1_files = []
    r2_files = []

    for filename in os.listdir(data_dir):
        if filename.endswith('.fastq.gz'):
            file_path = os.path.join(data_dir, filename)
            file_size = os.path.getsize(file_path)

            if '-R1' in filename or '_R1' in filename:
                r1_files.append((filename, file_size))
            elif '-R2' in filename or '_R2' in filename:
                r2_files.append((filename, file_size))
            else:
                # Process files without R1 or R2 individually
                base_name = filename.replace('.fastq.gz', '')
                files[base_name]['single'] = {'name': filename, 'size': file_size}

    # Process R1 files and find corresponding R2 files
    for r1_file, r1_size in r1_files:
        base_name = re.sub(r'[-_]R1', '', r1_file.replace('.fastq.gz', ''))
        r2_file = next((f for f, _ in r2_files if f.replace('R2', 'R1') == r1_file), None)

        if r2_file:
            r2_size = next(size for f, size in r2_files if f == r2_file)
            files[base_name]['R1'] = {'name': r1_file, 'size': r1_size}
            files[base_name]['R2'] = {'name': r2_file, 'size': r2_size}
        else:
            files[base_name]['unpaired'] = {'name': r1_file, 'size': r1_size}

    # Process remaining R2 files (those without corresponding R1)
    for r2_file, r2_size in r2_files:
        base_name = re.sub(r'[-_]R2', '', r2_file.replace('.fastq.gz', ''))
        if base_name not in files or 'R2' not in files[base_name]:
            files[base_name]['unpaired'] = {'name': r2_file, 'size': r2_size}

    return files


def write_paired_files(files, output_file, min_size):
    """Write paired files to output file if both files are present and exceed the minimum size"""
    with open(output_file, 'w') as pair_file:
        for sample, reads in files.items():
            if 'R1' in reads and 'R2' in reads:
                if reads['R1']['size'] >= min_size and reads['R2']['size'] >= min_size:
                    # Use cleaned sample name
                    cleaned_sample_name = clean_name(sample, r'_001')
                    pair_file.write(f"{cleaned_sample_name}\t{reads['R1']['name']}\t{reads['R2']['name']}\n")


def write_excluded_files(files, output_file, min_size):
    """Write unpaired files or files below minimum size to output file"""
    with open(output_file, 'w') as excluded_file:
        for sample, reads in files.items():
            if 'unpaired' in reads:
                excluded_file.write(f"{reads['unpaired']['name']} (Unpaired)\n")
            elif 'single' in reads:
                excluded_file.write(f"{reads['single']['name']} (No R1/R2 designation)\n")
            elif 'R1' in reads and 'R2' in reads:
                if reads['R1']['size'] < min_size or reads['R2']['size'] < min_size:
                    excluded_file.write(f"{reads['R1']['name']} (Size: {reads['R1']['size']} bytes)\n")
                    excluded_file.write(f"{reads['R2']['name']} (Size: {reads['R2']['size']} bytes)\n")


def process_fastq_files(data_dir, pair_output, excluded_output, min_size):
    """Process FASTQ files and output results"""
    files = get_fastq_files(data_dir)
    write_paired_files(files, pair_output, min_size)
    write_excluded_files(files, excluded_output, min_size)
    print("Processing completed.")


def main():
    """Main function"""
    data_dir = '../fastaq'
    pair_output_filename = 'paired_files.tsv'
    excluded_output_filename = 'excluded_files.txt'
    pair_output = os.path.join(data_dir, pair_output_filename)
    excluded_output = os.path.join(data_dir, excluded_output_filename)

    min_size = 100 * 1024  # Exclude files smaller than ex. 1Mb = 1000 * 1024 bytes from analysis
    process_fastq_files(data_dir, pair_output, excluded_output, min_size)


if __name__ == "__main__":
    main()