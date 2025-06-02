#!/opt/conda/bin/python

import os
import subprocess
import shutil


def clean_output_directory(output_dir):
    """Clean up the output directory"""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(f"Created output directory {output_dir}")


def get_sample_names(vcf_file):
    """Get list of sample names from VCF file"""
    try:
        cmd = ["bcftools", "query", "-l", vcf_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        samples = result.stdout.strip().split('\n')
        print(f"Number of samples: {len(samples)}")
        return samples
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to get sample names: {e}")
        return []


def check_fasta_content(file_path):
    """Check FASTA file content and verify sequence existence"""
    try:
        with open(file_path, 'r') as f:
            content = f.read().strip()
            return len(content) > 0
    except Exception:
        return False


def create_consensus_sequence(vcf_file, ref_fasta, sample, output_dir):
    """
    Create consensus sequences for each haplotype of the specified sample as separate files
    """
    try:
        success = True

        # Process each haplotype
        for haplotype in [1, 2]:
            # Set output filename (including H1 or H2)
            output_file = os.path.join(
                output_dir,
                f"{sample}_H{haplotype}.fasta"
            )

            # Build bcftools command
            cmd = [
                "bcftools", "consensus",
                "-f", ref_fasta,
                "-s", sample,
                "-H", str(haplotype),
                "--haplotype", str(haplotype),
                "-o", output_file,
                vcf_file
            ]

            # Generate consensus sequence
            subprocess.run(cmd, check=True)

            # Verify sequence existence
            if not check_fasta_content(output_file):
                print(
                    f"Warning: Failed to generate sequence for sample {sample} "
                    f"haplotype {haplotype}"
                )
                success = False
                continue

            # Modify header
            with open(output_file, 'r') as f:
                content = f.read()

            # Replace first line starting with > with new header
            new_header = f">{sample}_H{haplotype}\n"
            modified_content = new_header + \
                               content[content.find('\n') + 1:]

            with open(output_file, 'w') as f:
                f.write(modified_content)

            print(
                f"Created consensus sequence for sample {sample} "
                f"haplotype {haplotype}"
            )

        return success

    except subprocess.CalledProcessError as e:
        print(f"Error: Error occurred while processing sample {sample}: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error: Problem occurred while processing sample {sample}: {e}")
        return False


def main():
    # Set input/output file paths
    input_dir = "./0310_beagle_output"
    output_dir = "./0320_bcftools_output_consensus"
    ref_fasta = "./ref/ref_trimmed.fasta"
    input_vcf_filename = "beagle_result.vcf.gz"

    # Create full path for input VCF file
    input_vcf = os.path.join(input_dir, input_vcf_filename)

    # Prepare output directory
    clean_output_directory(output_dir)

    # Check input file existence
    if not os.path.exists(input_vcf):
        print(f"Error: Input VCF file {input_vcf} not found")
        return
    if not os.path.exists(ref_fasta):
        print(f"Error: Reference FASTA file {ref_fasta} not found")
        return

    # Get sample names
    samples = get_sample_names(input_vcf)
    if not samples:
        print("Error: No samples found")
        return

    # Create consensus sequence for each sample
    success_count = 0
    for sample in samples:
        if create_consensus_sequence(input_vcf, ref_fasta, sample, output_dir):
            success_count += 1

    print(f"Processing complete: {success_count}/{len(samples)} samples processed")


if __name__ == "__main__":
    main()