#!/opt/conda/bin/python

import os
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from typing import Dict, List, Tuple


def setup_output_dir(output_dir: str) -> None:
    """Output directory setup"""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)


def read_reference(ref_fasta: str) -> SeqRecord:
    """Read reference sequence"""
    return next(SeqIO.parse(ref_fasta, "fasta"))


def read_indel_data(indel_file: str) -> Dict[str, List[dict]]:
    """Read and organize indel data"""
    if not os.path.exists(indel_file):
        return {}

    df = pd.read_csv(indel_file, sep='\t')
    indel_dict = {}

    for _, row in df.iterrows():
        sample = row['Sample_Haplotype']
        if sample not in indel_dict:
            indel_dict[sample] = []

        indel_dict[sample].append({
            'pos': row['Position'] - 1,  # Convert to 0-based
            'ref': row['REF'],
            'alt': row['ALT'],
            'length_diff': len(row['ALT']) - len(row['REF'])
        })

    # Sort by position
    for sample in indel_dict:
        indel_dict[sample].sort(key=lambda x: x['pos'])

    return indel_dict


def align_sequences(ref_seq: str, sample_seq: str, indels: List[dict]) -> Tuple[str, str]:
    """Process sequence alignment"""
    ref_list = list(ref_seq)
    sample_list = list(sample_seq)

    offset = 0  # Position adjustment value for insertions/deletions

    for indel in indels:
        pos = indel['pos']
        adj_pos = pos + offset  # Position considering offset
        length_diff = indel['length_diff']

        if length_diff > 0:  # For insertions
            # Insert gaps into reference
            gap_pos = adj_pos + len(indel['ref'])
            ref_list[gap_pos:gap_pos] = ['-'] * length_diff
        elif length_diff < 0:  # For deletions
            # Insert gaps into sample
            gap_pos = adj_pos + len(indel['alt'])
            sample_list[gap_pos:gap_pos] = ['-'] * abs(length_diff)

        offset += length_diff

    return ''.join(ref_list), ''.join(sample_list)


def process_sample(sample_path: str,
                   ref_record: SeqRecord,
                   indel_dict: Dict[str, List[dict]],
                   output_dir: str) -> None:
    """Process sample"""
    # Load sample sequence
    sample_record = next(SeqIO.parse(sample_path, "fasta"))

    # Get sample name (remove extension)
    sample_name = Path(sample_path).stem

    # Get sequences
    ref_seq = str(ref_record.seq)
    sample_seq = str(sample_record.seq)

    # Process if indel data exists
    if sample_name in indel_dict:
        ref_seq, sample_seq = align_sequences(
            ref_seq,
            sample_seq,
            indel_dict[sample_name]
        )

    # Create output records
    records = [
        SeqRecord(
            Seq(ref_seq),
            id=ref_record.id,
            description=ref_record.description
        ),
        SeqRecord(
            Seq(sample_seq),
            id=sample_record.id,
            description=sample_record.description
        )
    ]

    # Create output filename
    output_path = os.path.join(output_dir, f"{sample_name}_aligned.fasta")

    # Write sequences
    SeqIO.write(records, output_path, "fasta")


def main():
    # Parameter settings
    output_dir = "./0410_alignment_output"
    input_dir = "./0320_bcftools_output_consensus"
    ref_fasta = "./ref/ref_trimmed.fasta"
    indel_file = './0330_indel_list/indel_summary.tsv'

    # Setup output directory
    setup_output_dir(output_dir)

    # Read reference sequence
    ref_record = read_reference(ref_fasta)

    # Read indel data
    indel_dict = read_indel_data(indel_file)

    # Process input files
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta'):
            sample_path = os.path.join(input_dir, filename)
            process_sample(sample_path, ref_record, indel_dict, output_dir)

    processed_files = len([f for f in os.listdir(input_dir) if f.endswith('.fasta')])
    print(f"Alignment completed: {processed_files} sequences processed")
    print(f"All aligned sequences have been saved in: {output_dir}")


if __name__ == "__main__":
    main()