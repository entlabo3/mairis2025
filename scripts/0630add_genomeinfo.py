#!/opt/conda/bin/python

import sys
import csv
import os
from typing import Optional
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class CDSRegion:
    """Data class representing CDS region"""
    start: int  # Genomic start position
    end: int  # Genomic end position
    phase: int  # Frame position
    name: str  # Exon name
    exon_number: float  # Exon number (supports decimal)

    @property
    def length(self) -> int:
        """Returns the length of CDS region"""
        return self.end - self.start + 1


class TranscriptPositionConverter:
    """Class for converting transcript position information"""

    def __init__(self, gff3_file: str):
        """Load transcript information from GFF3 file"""
        self.transcripts = defaultdict(list)
        self.cds_cumulative_lengths = defaultdict(list)
        self._load_gff3(gff3_file)
        self._calculate_cumulative_lengths()

    def _load_gff3(self, gff3_file: str) -> None:
        """Load GFF3 file and store CDS information for each transcript"""
        try:
            with open(gff3_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.strip() == '':
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) != 9:
                        continue

                    if fields[2] != 'CDS':
                        continue

                    # Parse attributes
                    attrs = {}
                    for attr in fields[8].split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attrs[key] = value

                    if 'Parent' not in attrs:
                        continue

                    # Get transcript IDs from Parent attribute
                    parent_ids = attrs['Parent'].split(',')

                    try:
                        # Parse exon number as float
                        exon_number = float(attrs.get('exon_number', '0'))
                    except ValueError:
                        print(f"Warning: Invalid exon number: {attrs.get('exon_number', '0')}")
                        continue

                    cds_region = CDSRegion(
                        start=int(fields[3]),
                        end=int(fields[4]),
                        phase=int(fields[7]),
                        name=attrs.get('Name', ''),
                        exon_number=exon_number
                    )

                    # Associate CDS region with each transcript
                    for parent_id in parent_ids:
                        self.transcripts[parent_id].append(cds_region)

            # Sort CDS regions by exon number for each transcript (supporting decimals)
            for transcript in self.transcripts.values():
                transcript.sort(key=lambda x: x.exon_number)

        except Exception as e:
            raise RuntimeError(f"Failed to read GFF3 file: {str(e)}")

    def _calculate_cumulative_lengths(self):
        """Calculate cumulative CDS lengths for each transcript"""
        for transcript_name, cdss in self.transcripts.items():
            cumulative_length = 0
            for cds in cdss:
                self.cds_cumulative_lengths[transcript_name].append(
                    (cumulative_length + 1, cumulative_length + cds.length)
                )
                cumulative_length += cds.length

    def cds_to_genome_position(self, cds_pos: int, transcript_name: Optional[str] = None) -> int:
        """
        Convert CDS position to genomic position

        Args:
            cds_pos: Position in CDS (1-based)
            transcript_name: Transcript name (uses default transcript if not specified)

        Returns:
            Genomic position (1-based, integer)
        """
        if transcript_name is None:
            transcript_name = sorted(self.transcripts.keys())[0]

        if transcript_name not in self.transcripts:
            raise ValueError(f"Specified transcript name does not exist: {transcript_name}")

        if cds_pos < 1:
            raise ValueError("CDS position must be greater than or equal to 1")

        # Process cumulative length info and actual CDS regions simultaneously
        for cds, (start_cum, end_cum) in zip(
                self.transcripts[transcript_name],
                self.cds_cumulative_lengths[transcript_name]
        ):
            if start_cum <= cds_pos <= end_cum:
                relative_pos = cds_pos - start_cum
                return cds.start + relative_pos

        raise ValueError(f"Specified CDS position {cds_pos} exceeds transcript length")


def process_variants_file(input_file: str, output_file: str, converter: TranscriptPositionConverter):
    """
    Process variant information file and add genomic positions

    Args:
        input_file: Input TSV file path
        output_file: Output TSV file path
        converter: Position converter object
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    try:
        with open(input_file, 'r') as in_f, open(output_file, 'w', newline='') as out_f:
            # Read header from input file
            reader = csv.DictReader(in_f, delimiter='\t')

            # Create new header (original header + new fields)
            fieldnames = list(reader.fieldnames)
            if 'genome_pos1' not in fieldnames:
                fieldnames.extend(['genome_pos1', 'genome_pos2', 'genome_pos3'])

            writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            # Count rows
            row_count = 0

            # Process each row
            for row in reader:
                row_count += 1
                try:
                    transcript_name = row.get('Transcript', '')

                    # Process CDS positions (use 0 for empty or invalid values)
                    cds_positions = []
                    for pos_key in ['Ref_CDS_pos1', 'Ref_CDS_pos2', 'Ref_CDS_pos3']:
                        try:
                            pos = int(row.get(pos_key, '0').strip())
                            cds_positions.append(pos)
                        except ValueError:
                            cds_positions.append(0)

                    # Convert to genomic positions
                    for i, pos in enumerate(cds_positions, 1):
                        if pos > 0:
                            try:
                                # Don't specify transcript name for Common
                                if transcript_name == 'Common':
                                    genome_pos = converter.cds_to_genome_position(pos)
                                else:
                                    genome_pos = converter.cds_to_genome_position(pos, transcript_name)
                                row[f'genome_pos{i}'] = str(genome_pos)
                            except ValueError as e:
                                print(f"Warning: Failed to convert position {pos} for {transcript_name}: {str(e)}")
                                row[f'genome_pos{i}'] = ''
                        else:
                            row[f'genome_pos{i}'] = ''

                    # Write row
                    writer.writerow(row)

                except Exception as e:
                    print(f"Warning: Error occurred while processing row {row_count}: {str(e)}")
                    continue

            print(f"Processed rows: {row_count}")

    except Exception as e:
        raise RuntimeError(f"Error occurred during file processing: {str(e)}")


def main():
    """Main function"""
    # File path settings
    gff3_file = "./ref/ref_trimmed.gff3"
    input_dir = "./0610_list_variants"
    input_file = os.path.join(input_dir, "variants_snp_mnp_processed.tsv")
    output_file = os.path.join(input_dir, "variants_snp_mnp_processed_with_genome_pos.tsv")

    try:
        # Check input file existence
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        if not os.path.exists(gff3_file):
            raise FileNotFoundError(f"GFF3 file not found: {gff3_file}")

        # Create position conversion object
        converter = TranscriptPositionConverter(gff3_file)

        # Display available transcripts
        print("Available transcripts:")
        for transcript_id in sorted(converter.transcripts.keys()):
            print(f"- {transcript_id}")

        # Process variant information file
        process_variants_file(input_file, output_file, converter)

        print("Processing completed")

    except Exception as e:
        print(f"An error occurred: {str(e)}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())