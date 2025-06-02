#!/opt/conda/bin/python
# improved_ref_files_processor_v5.py

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from pyfaidx import Faidx
import shutil
import subprocess

TRIM_SIZE = 999


def get_cds_info(gff3_file):
    cds_start, cds_end, cds_strand = float('inf'), 0, '+'
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#') or len(fields := line.strip().split('\t')) != 9 or fields[2] != 'CDS':
                continue
            start, end = int(fields[3]), int(fields[4])
            cds_start, cds_end = min(cds_start, start), max(cds_end, end)
            cds_strand = fields[6]
    return cds_start, cds_end, cds_strand


def trim_fasta(input_fasta, output_fasta, start, end):
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            trim_start = max(1, start - TRIM_SIZE)
            trim_end = min(len(record.seq), end + TRIM_SIZE)
            trimmed_seq = record.seq[trim_start - 1:trim_end]
            trimmed_record = SeqRecord(trimmed_seq, id=record.id,
                                       description=f"{record.description} (Trimmed {trim_start}-{trim_end})")
            SeqIO.write(trimmed_record, out_f, 'fasta')
    return trim_start


def update_gff3_coordinates_and_calculate_phase(input_gff3, output_gff3, offset, reverse=False, seq_length=None):
    """
    Update GFF3 file coordinates and recalculate phase.
    The phase of an exon start is the distance to the first base of the codon starting at that exon.
    Empty lines are excluded, only version header is retained.

    Args:
        input_gff3 (str): Input GFF3 file
        output_gff3 (str): Output GFF3 file
        offset (int): Coordinate adjustment value
        reverse (bool): Flag for reverse strand conversion
        seq_length (int): Sequence length (required when reverse=True)
    """
    version_line = None
    gene_lines = []
    mrna_lines = []
    cds_lines = []
    start_lines = []
    stop_lines = []

    # Read GFF3 file and classify each line appropriately
    with open(input_gff3, 'r') as in_f:
        for line in in_f:
            line = line.strip()
            if not line:  # Skip empty lines
                continue

            if line.startswith('##gff-version'):
                version_line = line
                continue

            if line.startswith('#'):  # Skip other comment lines
                continue

            fields = line.split('\t')
            if len(fields) != 9:  # Skip malformed lines
                continue

            if fields[2] == 'gene':
                gene_lines.append(fields)
            elif fields[2] == 'transcript':
                mrna_lines.append(fields)
            elif fields[2] == 'CDS':
                cds_lines.append(fields)
            elif fields[2] == 'start_codon':
                start_lines.append(fields)
            elif fields[2] == 'stop_codon':
                stop_lines.append(fields)

    # Sort CDS lines by exon number
    cds_lines.sort(key=lambda x: float(dict(item.split('=') for item in x[8].split(';'))['exon_number']))

    # Update coordinates and calculate phase
    cumulative_length = 0
    prev_exon_number = None

    for i, fields in enumerate(cds_lines):
        attributes = dict(item.split('=') for item in fields[8].split(';'))
        exon_number = float(attributes['exon_number'])
        start, end = int(fields[3]), int(fields[4])
        cds_length = end - start + 1

        # Update coordinates
        if reverse:
            new_start = seq_length - end + offset + 1
            new_end = seq_length - start + offset + 1
            fields[6] = '+' if fields[6] == '-' else '-'
        else:
            new_start = start - offset + 1
            new_end = end - offset + 1
        fields[3], fields[4] = str(max(1, new_start)), str(max(1, new_end))

        # Calculate phase
        if exon_number == 1:  # First exon
            phase = 0
            cumulative_length = cds_length
        elif exon_number == prev_exon_number:
            # Alternative splicing: use same phase as previous line
            phase = int(cds_lines[i - 1][7])
        else:
            # New exon: calculate new phase from cumulative length up to previous exon
            phase = (3 - (cumulative_length % 3)) % 3
            cumulative_length += cds_length

        # Update phase
        fields[7] = str(phase)
        prev_exon_number = exon_number

    # Final sort by POS (ascending)
    cds_lines.sort(key=lambda x: int(x[3]))

    # Update gene and mRNA line coordinates
    for lines in [gene_lines, mrna_lines]:
        for fields in lines:
            start, end = int(fields[3]), int(fields[4])
            if reverse:
                new_start = seq_length - end + offset + 1
                new_end = seq_length - start + offset + 1
                fields[6] = '+' if fields[6] == '-' else '-'
            else:
                new_start = start - offset + 1
                new_end = end - offset + 1
            fields[3], fields[4] = str(max(1, new_start)), str(max(1, new_end))

    # Write updated GFF3 file
    with open(output_gff3, 'w') as out_f:
        if version_line:
            out_f.write(version_line + '\n')
        for fields in gene_lines + mrna_lines + start_lines + cds_lines + stop_lines:
            out_f.write('\t'.join(fields) + '\n')


def reverse_complement_fasta(input_fasta, output_fasta):
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            rc_seq = record.seq.reverse_complement()
            rc_record = SeqRecord(rc_seq, id=record.id, description=f"{record.description} (Reverse Complement)")
            SeqIO.write(rc_record, out_f, 'fasta')


def get_fasta_length(input_fasta):
    return next(len(record.seq) for record in SeqIO.parse(input_fasta, 'fasta'))


def create_fai_index(fasta_file):
    Faidx(fasta_file)
    print(f"FAI index created for {fasta_file}")


def process_files(input_fasta, input_gff3, ref_dir):
    # Create ref directory if it doesn't exist
    os.makedirs(ref_dir, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(input_fasta))[0]
    base_path = os.path.join(ref_dir, base_name)

    cds_start, cds_end, cds_strand = get_cds_info(input_gff3)

    trimmed_fasta = f"{base_path}_trimmed.fasta"
    trim_offset = trim_fasta(input_fasta, trimmed_fasta, cds_start, cds_end)

    trimmed_gff3 = f"{base_path}_trimmed.gff3"
    update_gff3_coordinates_and_calculate_phase(input_gff3, trimmed_gff3, trim_offset)

    if cds_strand == '-':
        rc_fasta = f"{base_path}_trimmed_sense.fasta"
        reverse_complement_fasta(trimmed_fasta, rc_fasta)

        seq_length = get_fasta_length(trimmed_fasta)
        rc_gff3 = f"{base_path}_trimmed_sense.gff3"
        update_gff3_coordinates_and_calculate_phase(trimmed_gff3, rc_gff3, 0, reverse=True, seq_length=seq_length)

        create_fai_index(rc_fasta)
        print(f"Reverse complement files created: {rc_fasta}, {rc_gff3}")
        shutil.copy(rc_fasta, trimmed_fasta)
        shutil.copy(rc_gff3, trimmed_gff3)

    else:
        print("CDS strand is positive. No reverse complement needed.")

    create_fai_index(trimmed_fasta)

def extract_cds_fasta(gff3_file, fasta_file, output_file):
    """
    Creates a FASTA file containing only CDS regions from GFF3 and FASTA files.

    Args:
        gff3_file (str): Path to input GFF3 file
        fasta_file (str): Path to input FASTA file
        output_file (str): Path to output FASTA file
    """
    try:
        subprocess.run(["gffread", gff3_file, "-g", fasta_file, "-x", output_file], check=True)
        print(f"Created FASTA file with CDS regions: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to extract CDS regions - {e}")


def main():
    ref_dir = './ref'  # Define reference directory
    input_fasta = os.path.join(ref_dir, "ref.fasta")
    input_gff3 = os.path.join(ref_dir, "ref.gff3")
    process_files(input_fasta, input_gff3, ref_dir)

    ref_trimmed_gff3_file = os.path.join(ref_dir, "ref_trimmed.gff3")
    ref_trimmed_fasta_file = os.path.join(ref_dir, "ref_trimmed.fasta")
    ref_trimmed_output_file = os.path.join(ref_dir, "ref_trimmed_cds.fasta")

    extract_cds_fasta(ref_trimmed_gff3_file, ref_trimmed_fasta_file, ref_trimmed_output_file)




if __name__ == "__main__":
    main()