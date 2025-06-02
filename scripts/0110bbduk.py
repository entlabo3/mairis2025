#!/opt/conda/bin/python

import os
import re
import sys
import subprocess
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import logging
from typing import Tuple, List


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

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    log_file = os.path.join(output_dir, 'error.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

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


def run_bbduk(
        sample_name: str,
        fastq1: str,
        fastq2: str,
        input_dir: str,
        output_dir: str,
        adapter_file: str,
        threads: int,
        logger: logging.Logger
) -> Tuple[str, str]:
    """
    Run BBDuk to perform read quality control

    Args:
        sample_name: Sample name
        fastq1: R1 FASTQ filename
        fastq2: R2 FASTQ filename
        input_dir: Input directory
        output_dir: Output directory
        adapter_file: Path to adapter file
        threads: Number of threads
        logger: Logger

    Returns:
        Tuple[str, str]: Paths to processed R1, R2 files
    """
    try:
        # Generate safe sample name
        safe_sample = re.sub(r'[^\w\-\.]', '_', sample_name)

        # Set output filenames
        output1_paired = os.path.join(output_dir, f"{safe_sample}_R1_paired.fastq.gz")
        output2_paired = os.path.join(output_dir, f"{safe_sample}_R2_paired.fastq.gz")
        stats_file = os.path.join(output_dir, f"{safe_sample}_stats.txt")

        # Build BBDuk command
        bbduk_cmd = [
            "bbduk.sh",
            f"in1={os.path.join(input_dir, fastq1)}",
            f"in2={os.path.join(input_dir, fastq2)}",
            f"out1={output1_paired}",
            f"out2={output2_paired}",
            f"ref={adapter_file}",
            f"threads={threads}",
            "ktrim=r",  # Trim adapters from right side
            "k=23",  # kmer size
            "mink=11",  # Minimum kmer size
            "hdist=1",  # Allowed mismatches
            "tpe",  # Trim both reads to same length
            "tbo",  # Trim if reads overlap
            "qtrim=rl",  # Trim low quality bases from both ends
            "trimq=20",  # Quality trimming threshold
            "minlen=36",  # Minimum read length
            "maq=20",  # Minimum average quality
            "ordered",  # Maintain input order
            "qin=33",  # Input quality encoding
            f"stats={stats_file}"  # Output statistics
        ]

        # Log processing start
        logger.info(f"Processing sample: {sample_name}")
        logger.info(f"Command: {' '.join(bbduk_cmd)}")

        # Execute BBDuk
        result = subprocess.run(
            bbduk_cmd,
            check=True,
            capture_output=True,
            text=True
        )

        # Log processing results
        if result.stdout:
            logger.info(f"BBDuk output for {sample_name}:\n{result.stdout}")

        logger.info(f"Successfully processed: {sample_name}")
        return (output1_paired, output2_paired)

    except subprocess.CalledProcessError as e:
        logger.error(f"BBDuk error for {sample_name}:\n{e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error processing {sample_name}: {str(e)}")
        raise


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
        samples: List of sample information [(sample_name, fastq1, fastq2), ...]
        config: Configuration dictionary
        max_workers: Number of concurrent samples
        threads_per_sample: Threads per sample
        logger: Logger

    Returns:
        bool: Whether all samples were processed successfully
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for sample_name, fastq1, fastq2 in samples:
            futures.append(
                executor.submit(
                    run_bbduk,
                    sample_name,
                    fastq1,
                    fastq2,
                    config['input_dir'],
                    config['output_dir'],
                    config['adapter_file'],
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
        'input_dir': "../fastaq",
        'output_dir': "./0110_bbduk",
        'adapter_file': "../scripts/ngs_adapter.fasta",  # Explicitly specify adapter file path
        'input_list': "../fastaq/paired_files.tsv"
    }

    # Prepare output directory
    clean_output_directory(config['output_dir'])

    # Check adapter file existence
    if not os.path.exists(config['adapter_file']):
        raise FileNotFoundError(f"Adapter file not found: {config['adapter_file']}")

    # Configure logger
    logger = setup_logger(config['output_dir'])
    logger.info("Starting BBDuk processing")
    logger.info(f"Using adapter file: {config['adapter_file']}")

    try:
        # Load sample information
        with open(config['input_list'], 'r') as f:
            samples = [line.strip().split('\t') for line in f]

        # Calculate thread allocation
        total_threads = multiprocessing.cpu_count()
        available_threads = total_threads - 2  # Reserve 2 threads for system
        max_workers = max(1, available_threads // 6)  # Allocate 6 threads per sample
        threads_per_sample = 6

        logger.info(f"Total CPU threads: {total_threads}")
        logger.info(f"Reserved system threads: 2")
        logger.info(f"Available threads for processing: {available_threads}")
        logger.info(f"Max concurrent samples: {max_workers}")
        logger.info(f"Threads per sample: {threads_per_sample}")

        # Process samples
        if not process_samples(
                samples,
                config,
                max_workers,
                threads_per_sample,
                logger
        ):
            logger.error("Processing aborted")
            sys.exit(1)

        logger.info("All sample processing completed")

    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()