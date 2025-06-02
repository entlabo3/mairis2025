#!/opt/conda/bin/python

import os
import subprocess
import shutil


def clean_output_directory(output_dir):
    """
    Function to directly delete contents of output directory without confirmation
    and create the directory if it doesn't exist
    """
    if not os.path.exists(output_dir):
        print(f"Output directory {output_dir} does not exist. Creating new directory.")
        os.makedirs(output_dir)
        return

    for item in os.listdir(output_dir):
        item_path = os.path.join(output_dir, item)
        if os.path.isfile(item_path) or os.path.islink(item_path):
            os.unlink(item_path)
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)
        print(f"Deleted {item_path}")

    print(f"All contents of output directory {output_dir} have been deleted.")


def create_bam_list(bam_folder, output_file):
    """
    Function to create a list of BAM files
    """
    bam_files = []
    for file in sorted(os.listdir(bam_folder)):
        if file.endswith('.bam'):
            bam_path = os.path.join(bam_folder, file)
            # Check index file
            bai_path = f"{bam_path}.bai"
            if not os.path.exists(bai_path):
                print(f"Creating index for {file}")
                subprocess.run(['samtools', 'index', bam_path], check=True)
            bam_files.append(bam_path)

    # Create BAM list file
    with open(output_file, 'w') as f:
        for bam in bam_files:
            f.write(f"{bam}\n")

    return len(bam_files)


def run_freebayes(reference_fasta, bam_list, bed_file, output_vcf, threads=4):
    """
    Function to run FreeBayes and generate VCF

    Minimum alternate allele count (--min-alternate-count 3)
        Minimum number of occurrences required for an alternate allele to be called as a variant
        In this case, a variant must be observed at least 3 times
        Helps filter out sequencing errors and random noise
        Too low values increase false positives, too high values may miss real low-frequency variants

    Minimum alternate allele fraction (--min-alternate-fraction 0.05)
        Minimum frequency required for an alternate allele to be called (5%)
        Ratio of reads with alternate allele to total read count at that position
        Example: With 100 total reads, at least 5 reads must show the variant
        Important parameter for detecting heterozygous and mosaic variants

    Minimum coverage (--min-coverage 10)
        Minimum read coverage required for variant calling
        At least 10 reads must be mapped at that position
        Prevents uncertain variant calls in low-coverage regions
        Too low values lead to unreliable calls, too high values lead to data loss

    Maximum number of alleles (--use-best-n-alleles 4)
        Maximum number of alleles to consider at one position
        In this case, considers up to 4 alleles including the reference allele
        Affects analysis of complex polymorphisms and multi-allelic variants
        Too high values increase computation time and potential false detections
    """

    command = [
        'freebayes',
        '-f', reference_fasta,  # Reference file
        '-L', bam_list,  # BAM file list
        '-t', bed_file,  # Target region
        '--min-alternate-count', '3',  # Minimum alternate allele count
        '--min-alternate-fraction', '0.05',  # Minimum alternate allele fraction
        '--min-coverage', '10',  # Minimum coverage
        '--pooled-continuous',  # Treat samples as pooled
        '-v', output_vcf  # Output VCF file
    ]

    try:
        print(f"Running FreeBayes: {' '.join(command)}")
        subprocess.run(command, check=True)
        print("FreeBayes execution completed")
        return True
    except subprocess.CalledProcessError as e:
        print(f"FreeBayes execution error: {str(e)}")
        return False


def bgzip_and_index_vcf(vcf_file):
    """
    Function to compress VCF file with bgzip and create index
    """
    # Compression with bgzip
    bgzip_output = f"{vcf_file}.gz"
    subprocess.run(['bgzip', '-f', vcf_file], check=True)

    # Create index with tabix
    subprocess.run(['tabix', '-p', 'vcf', bgzip_output], check=True)

    return bgzip_output


def main():
    # Set input/output paths
    reference_fasta = "./ref/ref_trimmed.fasta"
    bed_file = "./ref/ref_trimmed_cds.bed"
    bam_folder = "./0120_strobealign"
    output_folder = "./0210_freebayes_output"
    temp_folder = os.path.join(output_folder, "./_temp")

    # Prepare output directories
    for folder in [output_folder, temp_folder]:
        clean_output_directory(folder)

    # Create BAM list file
    bam_list = os.path.join(temp_folder, "bam_list.txt")
    num_bams = create_bam_list(bam_folder, bam_list)

    if num_bams == 0:
        print("No BAM files found")
        return

    # Generate VCF with FreeBayes
    raw_vcf = os.path.join(temp_folder, "raw.vcf")
    if not run_freebayes(reference_fasta, bam_list, bed_file, raw_vcf):
        print("FreeBayes execution failed")
        return

    # Compress and index VCF
    final_vcf = os.path.join(output_folder, "variants.vcf")
    shutil.copy(raw_vcf, final_vcf)

    try:
        compressed_vcf = bgzip_and_index_vcf(final_vcf)
        print(f"Processing completed: {compressed_vcf}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to compress and index VCF: {str(e)}")


if __name__ == "__main__":
    main()