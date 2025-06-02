#!/opt/conda/bin/python
import gzip
from pathlib import Path
import sys
from typing import List, Tuple


def parse_genotype(gt: str) -> Tuple[str, str]:
    """
    Split genotype string into H1 and H2

    Args:
        gt: Genotype string (e.g., "0|1")

    Returns:
        Tuple[str, str]: (H1 value, H2 value)
    """
    try:
        h1, h2 = gt.split('|')
        return h1, h2
    except ValueError:
        print(f"Error: Invalid genotype format: {gt}")
        sys.exit(1)


def calculate_length_diff(ref: str, alt: str) -> int:
    """
    Calculate length difference between reference and alternative sequences

    Args:
        ref: Reference sequence
        alt: Alternative sequence

    Returns:
        int: Length difference (positive: insertion, negative: deletion)
    """
    return len(alt) - len(ref)


def has_length_variation(ref: str, alts: List[str]) -> bool:
    """
    Check if there are length differences between reference and alternative sequences

    Args:
        ref: Reference sequence
        alts: List of alternative sequences

    Returns:
        bool: True if there are length differences
    """
    ref_len = len(ref)
    return any(len(alt) != ref_len for alt in alts)


def process_vcf_file(input_path: Path) -> List[Tuple[str, str, int, str, str]]:
    """
    Process VCF file and return results list

    Args:
        input_path: Input VCF file path

    Returns:
        List[Tuple[str, str, int, str, str]]: Sample_Haplotype, position, length difference, REF, ALT list
    """
    results = []

    try:
        with gzip.open(input_path, 'rt') as f:
            # Process VCF file
            for line in f:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        # Get sample names (from 9th column onwards)
                        headers = line.strip().split('\t')
                        samples = headers[9:]
                    continue

                # Process data lines
                fields = line.strip().split('\t')
                pos = fields[1]
                ref = fields[3]
                alts = fields[4].split(',')

                # Skip if no length changes
                if not has_length_variation(ref, alts):
                    continue

                # Process each sample
                for i, sample in enumerate(samples):
                    genotype = fields[i + 9]
                    if genotype == './.':  # Skip missing values
                        continue

                    # Analyze genotype
                    h1, h2 = parse_genotype(genotype)

                    # Process H1
                    if h1 == '0':
                        length_diff_h1 = 0
                    else:
                        alt_idx = int(h1) - 1
                        alt_h1 = alts[alt_idx]
                        length_diff_h1 = calculate_length_diff(ref, alt_h1)
                        if length_diff_h1 != 0:
                            results.append((f"{sample}_H1", pos, length_diff_h1, ref, alt_h1))

                    # Process H2
                    if h2 == '0':
                        length_diff_h2 = 0
                    else:
                        alt_idx = int(h2) - 1
                        alt_h2 = alts[alt_idx]
                        length_diff_h2 = calculate_length_diff(ref, alt_h2)
                        if length_diff_h2 != 0:
                            results.append((f"{sample}_H2", pos, length_diff_h2, ref, alt_h2))

    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

    return results


def write_sorted_results(results: List[Tuple[str, str, int, str, str]], output_path: Path) -> None:
    """
    Sort results and write to output file

    Args:
        results: Results list
        output_path: Output file path
    """
    try:
        with open(output_path, 'w') as out:
            # Write header
            out.write("Sample_Haplotype\tPosition\tLength_diff\tREF\tALT\n")

            # Sort results by sample name
            sorted_results = sorted(results, key=lambda x: x[0])

            # Write results
            for sample_hap, pos, length_diff, ref, alt in sorted_results:
                out.write(f"{sample_hap}\t{pos}\t{length_diff:+d}\t{ref}\t{alt}\n")

    except Exception as e:
        print(f"Error writing results: {e}")
        sys.exit(1)


def main():
    """
    Main function
    """
    # Set input/output paths
    input_path = Path('./0310_beagle_output/beagle_result.vcf.gz')
    output_path = Path('./0330_indel_list/indel_summary.tsv')

    # Check input file existence
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)

    # Check output directory existence
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Process VCF file
    results = process_vcf_file(input_path)

    # Write results
    write_sorted_results(results, output_path)

    print(f"Processing complete. Results written to {output_path}")


if __name__ == "__main__":
    main()