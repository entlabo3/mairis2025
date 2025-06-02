#!/opt/conda/bin/python

import os
import subprocess
import shutil
import glob
import gzip
import re
from typing import TextIO, List


class VCFFixer:
    def __init__(self, input_vcf: str, output_vcf: str):
        """Initialize VCF file modification class"""
        self.input_vcf = input_vcf
        self.output_vcf = output_vcf
        self.sample_names: List[str] = []

    def process_header_line(self, line: str) -> str:
        """Process VCF header line"""
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            if len(fields) > 9:  # If sample columns exist
                self.sample_names = fields[9:]
                return line
        return line

    def fix_genotype_field(self, genotype_str: str) -> str:
        """Modify genotype field"""
        if genotype_str == ".":
            return "./.:.:.:.:.:.:.:."

        fields = genotype_str.split(':')
        if len(fields) != 8:
            return "./.:.:.:.:.:.:.:."

        gt = fields[0]
        if not re.match(r'^([0-9]+)[|/]([0-9]+)$', gt):
            fields[0] = "./."

        return ":".join(fields)

    def process_vcf_line(self, line: str) -> str:
        """Process one line of VCF"""
        if line.startswith('#'):
            return self.process_header_line(line)

        fields = line.strip().split('\t')
        if len(fields) < 10:
            return line

        for i in range(9, len(fields)):
            fields[i] = self.fix_genotype_field(fields[i])

        return '\t'.join(fields) + '\n'

    def process_file(self) -> None:
        """Process VCF file"""
        print(f"Starting modification of VCF file {self.input_vcf}")

        input_handle: TextIO
        if self.input_vcf.endswith('.gz'):
            input_handle = gzip.open(self.input_vcf, 'rt')
        else:
            input_handle = open(self.input_vcf, 'r')

        output_handle: TextIO
        if self.output_vcf.endswith('.gz'):
            output_handle = gzip.open(self.output_vcf, 'wt')
        else:
            output_handle = open(self.output_vcf, 'w')

        try:
            line_count = 0
            for line in input_handle:
                fixed_line = self.process_vcf_line(line)
                output_handle.write(fixed_line)
                line_count += 1
                if line_count % 1000 == 0:
                    print(f"Processed {line_count} lines")
        finally:
            input_handle.close()
            output_handle.close()
            if self.sample_names:
                print("Sample names:")
                for name in self.sample_names:
                    print(f"  - {name}")
            print(f"VCF file modification completed. Total {line_count} lines processed")


def clean_output_directory(output_dir):
    """Clean up output directory"""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(f"Output directory {output_dir} has been cleaned up.")


def run_beagle(input_file, output_file, memory="32g", nthreads=12):
    """Execute Beagle (imputation enabled)"""
    command = [
        "beagle",
        f"-Xmx{memory}",
        f"gt={input_file}",
        f"out={output_file}",
        f"nthreads={nthreads}",
        "impute=true"
    ]

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("Beagle execution successful (imputation enabled)")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Beagle execution error: {e}")
        print(f"Error output:\n{e.stderr}")
        raise


def create_index(vcf_file):
    """Create index for VCF file"""
    command = ["tabix", "-p", "vcf", vcf_file]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"Created index for {vcf_file}")
    except subprocess.CalledProcessError as e:
        print(f"Index creation error: {e}")
        print(f"Error output:\n{e.stderr}")
        raise


def compress_vcf(vcf_file):
    """Compress VCF file using bgzip"""
    try:
        subprocess.run(['bgzip', '-f', vcf_file], check=True)
        print(f"VCF file compressed: {vcf_file}.gz")
        return f"{vcf_file}.gz"
    except subprocess.CalledProcessError as e:
        print(f"VCF file compression error: {e}")
        print(f"Error output:\n{e.stderr}")
        raise


def main():
    # Path settings
    input_vcf = "./0210_freebayes_output/variants.vcf.gz"
    temp_dir = "./0310_beagle_output/_temp"
    output_dir = "./0310_beagle_output"

    # Prepare working directory
    clean_output_directory(output_dir)
    os.makedirs(temp_dir, exist_ok=True)

    try:
        # Modify VCF file
        fixed_vcf = os.path.join(temp_dir, "variants.fixed.vcf")
        fixer = VCFFixer(input_vcf, fixed_vcf)
        fixer.process_file()

        # Compress modified VCF file
        compressed_vcf = compress_vcf(fixed_vcf)

        # Create index for compressed VCF file
        create_index(compressed_vcf)

        # Execute Beagle
        output_prefix = os.path.join(output_dir, "beagle_result")
        run_beagle(compressed_vcf, output_prefix, memory="32g", nthreads=12)

        # Check Beagle output files and create index
        output_vcfs = glob.glob(os.path.join(output_dir, "*.vcf.gz"))
        if output_vcfs:
            create_index(output_vcfs[0])
        else:
            raise FileNotFoundError("Beagle output file not found")

        print("All processes completed successfully")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 1

    finally:
        # Delete temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    return 0


if __name__ == "__main__":
    exit(main())