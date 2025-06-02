#!/opt/conda/bin/python

import os
import re
import sys
import subprocess
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import logging
from typing import List, Tuple


def calculate_thread_allocation(total_cpu: int) -> tuple[int, int]:
    """
    Calculate thread allocation

    Args:
        total_cpu: Total available CPU cores

    Returns:
        tuple[int, int]: (Number of concurrent samples, Threads per sample)
    """
    # Reserve 2 threads for system
    available_threads = total_cpu - 2

    # Recommended threads per sample (8 threads)
    threads_per_sample = 8

    # Calculate number of concurrent samples
    max_workers = max(1, available_threads // threads_per_sample)

    return max_workers, threads_per_sample


def setup_logger(output_dir: str) -> logging.Logger:
    """
    Configure logger

    Args:
        output_dir: Output directory for log file

    Returns:
        logging.Logger: Configured logger
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Configure file handler
    log_file = os.path.join(output_dir, 'error.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Configure console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def clean_output_directory(output_dir: str) -> None:
    """
    Clean up output directory

    Args:
        output_dir: Directory path to clean up
    """
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)


def create_bam_index(bam_file: str, threads: int, logger: logging.Logger) -> None:
    """
    Create BAM file index

    Args:
        bam_file: BAM file path
        threads: Number of threads to use
        logger: Logger
    """
    try:
        cmd = ["samtools", "index", "-@", str(threads), bam_file]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Created BAM index for: {bam_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error creating BAM index: {e.stderr}")
        raise


def process_sample(
        sample_name: str,
        input_r1: str,
        input_r2: str,
        config: dict,
        threads: int,
        logger: logging.Logger
) -> str:
    """
    Execute alignment processing for one sample

    Args:
        sample_name: Sample name
        input_r1: R1 FASTQ file path
        input_r2: R2 FASTQ file path
        config: Configuration settings
        threads: Number of threads
        logger: Logger

    Returns:
        str: Output BAM file path
    """
    try:
        # Generate safe sample name
        safe_sample = re.sub(r'[^\w\-\.]', '_', sample_name)
        logger.info(f"Processing sample: {safe_sample}")

        # Output BAM file path
        output_bam = os.path.join(config['output_dir'], f"{safe_sample}.sorted.bam")

        # Build strobealign command and execute pipeline
        strobealign_cmd = f"""
        strobealign -t {threads} {config['reference']} {input_r1} {input_r2} |
        samtools addreplacerg -r '@RG\\tID:{safe_sample}\\tSM:{safe_sample}\\tPL:ILLUMINA' - |
        samtools sort -@ {threads} -m 2G - > {output_bam}
        """

        # Execute command
        subprocess.run(strobealign_cmd, shell=True, check=True)

        # Create BAM index
        create_bam_index(output_bam, threads, logger)

        logger.info(f"Successfully processed: {safe_sample}")
        return output_bam

    except Exception as e:
        logger.error(f"Error processing sample {sample_name}: {str(e)}")
        raise


def get_input_files(input_dir: str) -> List[Tuple[str, str, str]]:
    """
    Get paired-end FASTQ file information from input directory

    Args:
        input_dir: Input directory path

    Returns:
        List[Tuple[str, str, str]]: [(sample_name, R1_path, R2_path), ...]
    """
    samples = []
    r1_files = {}
    r2_files = {}

    # Collect files
    for file in os.listdir(input_dir):
        if not file.endswith('_paired.fastq.gz'):
            continue

        full_path = os.path.join(input_dir, file)
        if '_R1_' in file:
            sample_name = file.split('_R1_')[0]
            r1_files[sample_name] = full_path
        elif '_R2_' in file:
            sample_name = file.split('_R2_')[0]
            r2_files[sample_name] = full_path

    # Create pairs
    for sample in sorted(r1_files.keys()):
        if sample in r2_files:
            samples.append((sample, r1_files[sample], r2_files[sample]))

    return samples


def process_samples(
        samples: List[Tuple[str, str, str]],
        config: dict,
        max_workers: int,
        threads_per_sample: int,
        logger: logging.Logger
) -> bool:
    """
    Execute parallel processing of multiple samples

    Args:
        samples: Sample information list
        config: Configuration settings
        max_workers: Number of concurrent samples
        threads_per_sample: Threads per sample
        logger: Logger

    Returns:
        bool: Whether all samples were processed successfully
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for sample_name, r1_file, r2_file in samples:
            futures.append(
                executor.submit(
                    process_sample,
                    sample_name,
                    r1_file,
                    r2_file,
                    config,
                    threads_per_sample,
                    logger
                )
            )

        # Process results
        success = True
        for future in as_completed(futures):
            try:
                future.result()
            except Exception:
                success = False
                break

        return success


def main():
    """Main process"""
    # Configuration
    config = {
        'input_dir': "./0110_bbduk",
        'output_dir': "./0120_strobealign",
        'reference': "./ref/ref_trimmed.fasta"
    }

    # Prepare output directory
    clean_output_directory(config['output_dir'])

    # Configure logger
    logger = setup_logger(config['output_dir'])
    logger.info("Starting alignment process")

    try:
        # Get input files
        samples = get_input_files(config['input_dir'])
        logger.info(f"Found {len(samples)} samples to process")

        # Calculate thread allocation
        total_threads = multiprocessing.cpu_count()
        max_workers, threads_per_sample = calculate_thread_allocation(total_threads)

        logger.info(f"Total CPU threads: {total_threads}")
        logger.info(f"Reserved system threads: 2")
        logger.info(f"Available threads: {total_threads - 2}")
        logger.info(f"Threads per sample: {threads_per_sample}")
        logger.info(f"Max concurrent samples: {max_workers}")

        # Process samples
        if not process_samples(samples, config, max_workers, threads_per_sample, logger):
            logger.error("Processing was interrupted")
            sys.exit(1)

        logger.info("All sample processing completed")

    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()