#!/opt/conda/bin/python

import os
from dataclasses import dataclass
from typing import Dict, List, Set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
from collections import defaultdict


def clean_output_directory(output_dir: str) -> None:
    """
    Clean up the output directory function

    Args:
        output_dir (str): Path to the output directory to clean up
    """
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

@dataclass
class ExonInfo:
    """Class to hold exon information"""
    name: str
    start: int
    end: int
    phase: int
    transcripts: Set[str]


class TranscriptBuilder:
    """Class to manage transcript construction"""

    def __init__(self, gff3_file: str):
        self.exons: Dict[str, ExonInfo] = {}  # name -> ExonInfo
        self.transcript_exons: Dict[str, List[str]] = defaultdict(list)  # transcript_id -> [exon_names]
        self.parse_gff3(gff3_file)

    def parse_gff3(self, gff3_file: str) -> None:
        """Extract exon information from GFF3 file"""
        try:
            with open(gff3_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue

                    fields = line.strip().split('\t')
                    if fields[2] != 'CDS':
                        continue

                    # Parse attributes
                    attrs = dict(attr.split('=') for attr in fields[8].split(';'))
                    name = attrs['Name']
                    transcripts = set(attrs['Parent'].split(','))

                    # Save exon information
                    self.exons[name] = ExonInfo(
                        name=name,
                        start=int(fields[3]),
                        end=int(fields[4]),
                        phase=int(fields[7]) if fields[7] != '.' else 0,
                        transcripts=transcripts
                    )

                    # Create exon list for each transcript
                    for transcript in transcripts:
                        self.transcript_exons[transcript].append(name)

            # Sort exon list by start position
            for transcript in self.transcript_exons:
                self.transcript_exons[transcript].sort(
                    key=lambda x: self.exons[x].start
                )

        except Exception as e:
            raise RuntimeError(f"Error occurred during GFF3 file analysis: {str(e)}")

    def get_exon_sequences_from_file(self, fasta_file: str) -> Dict[str, str]:
        """Extract exon sequences from FASTA file"""
        sequences = {}
        try:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # Extract exon name (e.g., get 'Exon01' from 'Fukuoka-3_S84_H1_Exon01')
                exon_name = record.id.split('_')[-1]
                sequences[exon_name] = str(record.seq)
        except Exception as e:
            raise RuntimeError(f"Error reading FASTA file {fasta_file}: {str(e)}")
        return sequences

    def build_transcripts(self, exon_sequences: Dict[str, str], sample_id: str) -> List[SeqRecord]:
        """Build transcripts"""
        transcripts = []
        for transcript_id, exon_names in self.transcript_exons.items():
            try:
                # Concatenate sequences
                sequence = ''
                for exon_name in exon_names:
                    if exon_name not in exon_sequences:
                        raise RuntimeError(f"Sequence not found for exon {exon_name}")
                    sequence += exon_sequences[exon_name]

                # Create SeqRecord
                record = SeqRecord(
                    Seq(sequence),
                    id=f"{sample_id}_{transcript_id}",
                    description=f"Transcript {transcript_id} from {sample_id}"
                )
                transcripts.append(record)

            except Exception as e:
                raise RuntimeError(f"Error during transcript {transcript_id} construction: {str(e)}")

        return transcripts


def process_fasta_file(builder: TranscriptBuilder,
                       input_file: str,
                       output_dir: str) -> None:
    """Process one FASTA file"""
    try:
        # Get sample ID (e.g., 'Fukuoka-3_S84_H1')
        sample_id = os.path.basename(input_file).replace('.fasta', '')

        # Load exon sequences
        exon_sequences = builder.get_exon_sequences_from_file(input_file)

        # Build transcripts
        transcripts = builder.build_transcripts(exon_sequences, sample_id)

        # Output results
        output_file = os.path.join(output_dir, f"{sample_id}_transcripts.fasta")
        SeqIO.write(transcripts, output_file, "fasta")

        print(f"{sample_id}: Generated {len(transcripts)} transcripts")

    except Exception as e:
        print(f"Error ({sample_id}): {str(e)}")


def main():
    """Main process"""
    # Set paths
    input_dir = "./0420_extract_exons"
    output_dir = "./0430_transcripts"
    gff3_file = "./ref/ref_trimmed.gff3"

    # Create output directory
    clean_output_directory(output_dir)

    try:
        # Initialize transcript builder
        builder = TranscriptBuilder(gff3_file)

        # Get list of input files
        input_files = sorted([
            os.path.join(input_dir, f) for f in os.listdir(input_dir)
            if f.endswith('.fasta')
        ])

        if not input_files:
            print(f"Warning: No FASTA files found in {input_dir}")
            return 1

        print(f"Number of files to process: {len(input_files)}")

        # Process each file
        for i, input_file in enumerate(input_files, 1):
            print(f"\nProcessing ({i}/{len(input_files)}): {os.path.basename(input_file)}")
            process_fasta_file(builder, input_file, output_dir)

        print(f"\nAll file processing completed. Results saved in {output_dir}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())