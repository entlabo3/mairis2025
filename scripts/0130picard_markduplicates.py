#!/opt/conda/bin/python

import os
import subprocess
import shutil
import sys
import psutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from typing import List, Tuple


def setup_logger(output_dir: str) -> logging.Logger:
    """
    Configure logging settings

    Args:
        output_dir: Directory to save log file

    Returns:
        logging.Logger: Configured logger
    """
    logger = logging.getLogger('picard_processor')
    logger.setLevel(logging.INFO)

    # File handler configuration
    log_file = os.path.join(output_dir, 'processing.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    # Console handler configuration
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Formatter configuration
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers
    logger.addHandler(file_handler)
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


def get_memory_settings() -> dict:
    """
    Calculate memory settings for each tool based on system memory
    Allocate memory considering 4 concurrent processes

    Returns:
        dict: Memory settings
    """
    total_mem = psutil.virtual_memory().total / (1024 ** 3)  # Convert to GB
    mem_per_process = total_mem / 5  # 4 processes + system reserve

    return {
        'picard_mem': f"{int(mem_per_process * 0.7)}G",  # 70% memory per process
        'samtools_mem': f"{int(mem_per_process * 0.2)}G"  # 20% memory per process
    }


def convert_gff3_to_bed(gff3_file: str, bed_file: str, logger: logging.Logger) -> None:
    """
    Create BED file from GFF3 file for CDS regions

    Args:
        gff3_file: Input GFF3 file path
        bed_file: Output BED file path
        logger: Logger
    """
    try:
        with open(gff3_file, 'r') as gff, open(bed_file, 'w') as bed:
            for line in gff:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 8 and fields[2] == 'CDS':
                    chrom = fields[0]
                    start = str(int(fields[3]) - 1)  # BED is 0-based
                    end = fields[4]
                    name = "CDS"
                    score = "0"
                    strand = fields[6]
                    bed.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
        logger.info(f"Created BED file: {bed_file}")
    except Exception as e:
        logger.error(f"Failed to create BED file: {str(e)}")
        raise


def get_input_files(input_dir: str) -> List[Tuple[str, str]]:
    """
    Get list of BAM files from input directory

    Args:
        input_dir: Input directory path

    Returns:
        List[Tuple[str, str]]: List of (sample_name, BAM_file_path)
    """
    samples = []
    for file in sorted(os.listdir(input_dir)):
        if file.endswith('.bam'):
            sample_name = os.path.splitext(file)[0]
            bam_path = os.path.join(input_dir, file)
            samples.append((sample_name, bam_path))
    return samples


def process_sample(
        input_bam: str,
        sample_name: str,
        config: dict,
        mem_settings: dict,
        logger: logging.Logger
) -> str:
    """
    Process a single sample

    Args:
        input_bam: Input BAM file path
        sample_name: Sample name
        config: Configuration settings
        mem_settings: Memory settings
        logger: Logger

    Returns:
        str: Output BAM file path
    """
    try:
        logger.info(f"Processing sample: {sample_name}")

        # Temporary file for duplicate removal
        dedup_bam = os.path.join(config['temp_dir'], f"{sample_name}.dedup.bam")
        metrics_file = os.path.join(config['metrics_dir'], f"{sample_name}_duplicate_metrics.tsv")

        # Execute Picard command
        picard_cmd = [
            "picard",
            f"-Xmx{mem_settings['picard_mem']}",
            "MarkDuplicates",
            f"INPUT={input_bam}",
            f"OUTPUT={dedup_bam}",
            f"METRICS_FILE={metrics_file}",
            "REMOVE_DUPLICATES=true",
            "VALIDATION_STRINGENCY=LENIENT",
            "ASSUME_SORTED=true",
            "CREATE_INDEX=true",
            "USE_JDK_DEFLATER=true",
            "USE_JDK_INFLATER=true"
        ]

        subprocess.run(picard_cmd, check=True, capture_output=True, text=True)
        logger.info(f"Completed duplicate marking for: {sample_name}")

        # CDS region filtering and sorting
        # Change output filename to .sorted.markdup.bam
        final_bam = os.path.join(
            config['output_dir'],
            f"{os.path.splitext(os.path.basename(input_bam))[0].replace('.sorted', '')}.sorted.markdup.bam"
        )

        filter_cmd = f"""
        samtools view -@ 2 -b -L {config['cds_bed']} {dedup_bam} |
        samtools sort -@ 2 -m {mem_settings['samtools_mem']} - > {final_bam}
        """

        subprocess.run(filter_cmd, shell=True, check=True, capture_output=True, text=True)

        # Create index
        subprocess.run(
            ["samtools", "index", "-@", "2", final_bam],
            check=True,
            capture_output=True,
            text=True
        )

        # Remove temporary files
        os.remove(dedup_bam)
        dedup_bai = dedup_bam.replace('.bam', '.bai')
        if os.path.exists(dedup_bai):
            os.remove(dedup_bai)

        logger.info(f"Successfully processed: {sample_name}")
        return final_bam

    except subprocess.CalledProcessError as e:
        logger.error(f"Command error processing {sample_name}: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error processing sample {sample_name}: {str(e)}")
        raise


def process_samples(
        samples: List[Tuple[str, str]],
        config: dict,
        mem_settings: dict,
        logger: logging.Logger
) -> bool:
    """
    Process multiple samples in parallel (maximum 4 samples simultaneously)

    Args:
        samples: List of sample information [(sample_name, bam_file), ...]
        config: Configuration settings
        mem_settings: Memory settings
        logger: Logger

    Returns:
        bool: Whether all samples were processed successfully
    """
    failed_samples = []
    max_workers = 4  # Fixed to 4 concurrent processes

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_sample = {
            executor.submit(
                process_sample,
                bam_file,
                sample_name,
                config,
                mem_settings,
                logger
            ): sample_name
            for sample_name, bam_file in samples
        }

        for future in as_completed(future_to_sample):
            sample_name = future_to_sample[future]
            try:
                future.result()
            except Exception as e:
                failed_samples.append(sample_name)
                logger.error(f"Failed to process sample {sample_name}: {str(e)}")

    if failed_samples:
        logger.error(f"Failed samples: {', '.join(failed_samples)}")
        return False

    return True


def main():
    """Main process"""
    # Configuration
    config = {
        'input_dir': "./0120_strobealign",
        'output_dir': "./0130_picard_markduplicates",
        'metrics_dir': "./0130_picard_metrics",
        'temp_dir': "/opt/temp",
        'gff3_file': "./ref/ref_trimmed.gff3",
        'cds_bed': "./ref/ref_trimmed_cds.bed"
    }

    # Prepare output directories
    for dir_path in [config['output_dir'], config['temp_dir'], config['metrics_dir']]:
        clean_output_directory(dir_path)

    # Configure logger
    logger = setup_logger(config['output_dir'])
    logger.info("Starting duplicate marking and CDS filtering process")

    try:
        # Get memory settings
        mem_settings = get_memory_settings()
        logger.info(f"Memory settings: Picard: {mem_settings['picard_mem']}, "
                    f"Samtools: {mem_settings['samtools_mem']}")

        # Create BED file from GFF3
        convert_gff3_to_bed(config['gff3_file'], config['cds_bed'], logger)

        # Get input files
        samples = get_input_files(config['input_dir'])
        logger.info(f"Found {len(samples)} samples to process")
        logger.info("Processing with maximum 4 concurrent samples")

        # Process samples
        if not process_samples(samples, config, mem_settings, logger):
            logger.error("Failed to process some samples")
            sys.exit(1)

        logger.info("All sample processing completed")

    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()