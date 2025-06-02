#!/opt/conda/bin/python

import gzip
import os
import sys
import numpy as np
import shutil
from typing import List, Dict, Tuple


def clean_output_directory(output_dir: str) -> None:
    """Clean up output directory"""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")


class HaplotypeDistanceCalculator:
    """Class for calculating haplotype distances"""

    def __init__(self):
        self.samples = []  # List of sample names
        self.haplotypes = {}  # Dictionary to store haplotype data
        self.reference_haplotype = []  # Reference haplotype
        self.distance_matrix = None  # Distance matrix
        self.hap_names = None  # List of haplotype names

    def read_vcf_gz(self, vcf_path: str) -> None:
        """
        Read a gzipped VCF file

        Args:
            vcf_path: Path to input VCF file
        """
        try:
            with gzip.open(vcf_path, 'rt') as f:
                # Find header line
                for line in f:
                    if line.startswith('#CHROM'):
                        # Get sample names from 9th column onwards
                        self.samples = line.strip().split('\t')[9:]
                        break

                # Initialize haplotypes for each sample
                for sample in self.samples:
                    self.haplotypes[f"{sample}_H1"] = []
                    self.haplotypes[f"{sample}_H2"] = []

                # Load variant data
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    ref_allele = fields[3]  # Reference allele
                    alt_alleles = fields[4].split(',')  # Alternative alleles
                    genotypes = fields[9:]  # Sample genotype information

                    # Save reference allele as '0'
                    self.reference_haplotype.append('0')

                    # Process genotype for each sample
                    for sample_idx, gt in enumerate(genotypes):
                        sample = self.samples[sample_idx]

                        # Split genotype into haplotypes
                        try:
                            alleles = gt.split(':')[0].replace('|', '/')
                            if '/' in alleles:
                                hap1, hap2 = alleles.split('/')
                            else:
                                hap1 = hap2 = alleles

                            self.haplotypes[f"{sample}_H1"].append(hap1)
                            self.haplotypes[f"{sample}_H2"].append(hap2)
                        except IndexError:
                            print(f"Warning: Error processing genotype {gt} for sample {sample}")
                            continue

        except Exception as e:
            print(f"Error: An error occurred while reading VCF file: {str(e)}")
            sys.exit(1)

    def calculate_distances(self) -> None:
        """Calculate distances between all haplotype pairs"""
        # Create list of all haplotype names (including reference)
        self.hap_names = ["Reference"]
        for sample in self.samples:
            self.hap_names.extend([f"{sample}_H1", f"{sample}_H2"])

        n_haps = len(self.hap_names)
        self.distance_matrix = np.zeros((n_haps, n_haps))

        # Calculate distances with reference
        for i, hap_name in enumerate(self.hap_names[1:], 1):
            distance = self._calculate_pair_distance(
                self.reference_haplotype,
                self.haplotypes[hap_name]
            )
            self.distance_matrix[0, i] = distance
            self.distance_matrix[i, 0] = distance

        # Calculate distances between other haplotypes
        for i in range(1, n_haps):
            for j in range(i + 1, n_haps):
                distance = self._calculate_pair_distance(
                    self.haplotypes[self.hap_names[i]],
                    self.haplotypes[self.hap_names[j]]
                )
                self.distance_matrix[i, j] = distance
                self.distance_matrix[j, i] = distance

    def _calculate_pair_distance(self, hap1: List[str], hap2: List[str]) -> int:
        """
        Calculate distance between two haplotypes

        Args:
            hap1, hap2: Two haplotype arrays to compare

        Returns:
            int: Distance between two haplotypes
        """
        if len(hap1) != len(hap2):
            raise ValueError("Haplotype lengths do not match")

        return sum(1 for a, b in zip(hap1, hap2) if a != b)

    def save_matrix(self, output_path: str) -> None:
        """
        Save distance matrix to file
        Output sample names comma-separated in first line, numbers only from second line onwards

        Args:
            output_path: Output file path
        """
        try:
            with open(output_path, 'w') as f:
                # First line: Write haplotype names comma-separated
                header = ','.join(self.hap_names)
                f.write(f"{header}\n")

                # Write each row of distance matrix
                for row in self.distance_matrix:
                    row_str = ','.join(f"{val:.1f}" for val in row)
                    f.write(f"{row_str}\n")

            print(f"Distance matrix saved: {output_path}")

        except Exception as e:
            print(f"Error: An error occurred while saving distance matrix: {str(e)}")
            sys.exit(1)

    def print_matrix_info(self) -> None:
        """Display distance matrix information"""
        if self.distance_matrix is not None:
            shape = self.distance_matrix.shape
            print(f"Distance matrix size: {shape[0]} x {shape[1]}")
            print(f"Number of haplotypes: {len(self.hap_names)}")
            print(f"Number of samples: {len(self.samples)}")
        else:
            print("Distance matrix has not been calculated yet")


def main():
    """Main process"""
    # Set input/output paths
    input_path = './0310_beagle_output/beagle_result.vcf.gz'
    output_dir = './0510_distance_matrix'
    output_path = os.path.join(output_dir, "haplotype_distance_matrix.csv")

    # Create output directory
    clean_output_directory(output_dir)

    # Execute process
    calculator = HaplotypeDistanceCalculator()

    print("Starting VCF file reading...")
    calculator.read_vcf_gz(input_path)

    print("Starting distance matrix calculation...")
    calculator.calculate_distances()

    # Display matrix information
    calculator.print_matrix_info()

    print("Starting result saving...")
    calculator.save_matrix(output_path)

    print("Process completed")


if __name__ == "__main__":
    main()