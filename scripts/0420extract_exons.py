#!/opt/conda/bin/python

import os
from dataclasses import dataclass
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


@dataclass
class ExonRegion:
    """Data class representing exon regions"""
    start: int
    end: int
    name: str


def clean_output_directory(output_dir):
    """Clean up the output directory"""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(f"Created output directory {output_dir}")


def parse_gff3(gff3_file: str) -> list[ExonRegion]:
    """
    Extract CDS region information from GFF3 file

    Args:
        gff3_file: Path to GFF3 file

    Returns:
        List of CDS regions (sorted by start position)
    """
    cds_regions = []
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) != 9 or fields[2] != 'CDS':
                continue

            # Extract name from attributes
            attributes = dict(item.split('=') for item in fields[8].split(';'))
            name = attributes.get('Name', f"exon{len(cds_regions) + 1}")

            # GFF3 uses 1-based coordinates, so keep as is
            cds_regions.append(ExonRegion(
                start=int(fields[3]),
                end=int(fields[4]),
                name=name
            ))

    return sorted(cds_regions, key=lambda x: x.start)


def create_coordinate_map(ref_seq: str, query_seq: str) -> dict[int, int]:
    """
    Create coordinate conversion map from aligned sequences
    Returns in 1-based coordinate system

    Args:
        ref_seq: Reference sequence (with gaps)
        query_seq: Query sequence (with gaps)

    Returns:
        Conversion map from reference to query coordinates (1-based)
    """
    ref_pos = 1  # Reference sequence actual position (1-based)
    query_pos = 1  # Query sequence actual position (1-based)
    coord_map = {}

    for ref_char, query_char in zip(ref_seq, query_seq):
        if ref_char != '-':
            if query_char != '-':
                coord_map[ref_pos] = query_pos
            ref_pos += 1

        if query_char != '-':
            query_pos += 1

    return coord_map


def extract_exons(alignment_file: str, gff3_file: str, output_file: str):
    """
    Extract exon sequences from alignment results and save to FASTA file

    Args:
        alignment_file: Alignment result FASTA file
        gff3_file: Path to GFF3 file
        output_file: Path to output FASTA file
    """
    # Get exon regions from GFF3 (1-based coordinates)
    exon_regions = parse_gff3(gff3_file)

    # Load alignment results
    records = list(SeqIO.parse(alignment_file, 'fasta'))
    if len(records) != 2:
        raise ValueError(f"Alignment file must contain 2 sequences ({len(records)} sequences found)")

    # Get reference and sample sequences
    ref_seq = str(records[0].seq)
    query_seq = str(records[1].seq)
    sample_name = records[1].id

    # Create coordinate conversion map (1-based)
    coord_map = create_coordinate_map(ref_seq, query_seq)

    # Query sequence with gaps removed
    clean_query_seq = query_seq.replace('-', '')

    # Extract each exon sequence
    exon_records = []
    for region in exon_regions:
        # Convert reference coordinates to query coordinates (1-based)
        query_start = coord_map.get(region.start)
        query_end = coord_map.get(region.end)

        if query_start is None or query_end is None:
            print(f"Warning: Failed to convert coordinates {region.start}-{region.end} for exon {region.name}")
            continue

        # Extract exon sequence (convert to 0-based for slicing)
        exon_seq = clean_query_seq[query_start - 1:query_end]

        # Create FASTA record
        record = SeqRecord(
            Seq(exon_seq),
            id=f"{sample_name}_{region.name}",
            description=f"coordinates={query_start}-{query_end}"
        )
        exon_records.append(record)

    # Save results
    SeqIO.write(exon_records, output_file, 'fasta')


def process_sample_files(input_dir: str, output_dir: str, gff3_file: str):
    """
    Process all alignment files in the input directory

    Args:
        input_dir: Directory containing alignment files
        output_dir: Output directory
        gff3_file: Path to GFF3 file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Process alignment files
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta') or filename.endswith('.fa'):
            alignment_file = os.path.join(input_dir, filename)

            # Get sample name (ID of second sequence)
            with open(alignment_file) as f:
                records = list(SeqIO.parse(f, 'fasta'))
                if len(records) != 2:
                    print(f"Skip: {filename} - not 2 sequences")
                    continue
                sample_name = records[1].id

            # Set output filename
            output_file = os.path.join(output_dir, f"{sample_name}.fasta")

            try:
                extract_exons(alignment_file, gff3_file, output_file)
                print(f"Processing complete: {filename} -> {output_file}")
            except Exception as e:
                print(f"Error: Problem occurred while processing {filename}: {str(e)}")


def main():
    """Main process"""
    # Set input/output directories and GFF3 file
    input_dir = "./0410_alignment_output"
    output_dir = "./0420_extract_exons"
    gff3_file = "./ref/ref_trimmed.gff3"

    clean_output_directory(output_dir)

    try:
        process_sample_files(input_dir, output_dir, gff3_file)
        print("All processing completed")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


if __name__ == "__main__":
    main()