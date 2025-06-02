#!/opt/conda/bin/python

from Bio import SeqIO
import os
import shutil
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein


def reverse_complement_sequences(input_file, output_file):
    """
    Function to read sequences from a FASTA file, convert to reverse complement, and save to a new file
    """
    try:
        # Read sequences from input file
        sequences = list(SeqIO.parse(input_file, "fasta"))

        # Convert each sequence to reverse complement
        reversed_sequences = []
        for seq in sequences:
            reversed_seq = seq.reverse_complement(
                id=seq.id,
                description="Reverse complement of " + seq.description
            )
            reversed_sequences.append(reversed_seq)

        # Write converted sequences to output file
        SeqIO.write(reversed_sequences, output_file, "fasta")
        print(f"Converted genome sequences to reverse complement and saved to {output_file}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


def split_genome_seq(genome_seq, chunk_size=50):
    """
    Split a long string into 50-character substrings and return a list of tuples with starting positions
    Args:
        genome_seq: Target string for splitting
    Returns:
        List of tuples [(start_position, substring), ...]
    """
    # Get length of genome_seq
    seq_length = len(genome_seq)
    # Create tuples of 50-character substrings and their start positions
    result = [(i, genome_seq[i:i + chunk_size]) for i in range(seq_length - chunk_size)]
    return result


def find_exon_start(query, splitted_genome_seq_list):
    """
    Return the index of the substring with the smallest Levenshtein distance to the query string

    Args:
        query: Query string
        splitted_genome_seq_list: List of tuples [(start_position, substring), ...]

    Returns:
        Start position of the most similar substring (False if not found)
    """
    # Use rapidfuzz.process.extractOne to get the closest substring and its score
    best_match, score, index = process.extractOne(
        query,  # Query string
        [seq[1] for seq in splitted_genome_seq_list],  # List of substrings only
        score_cutoff=0  # Consider only scores above 0 (has some similarity)
    )

    # Return start position if match found
    if best_match:
        return splitted_genome_seq_list[index][0]  # Return start position
    else:
        return False  # Return False if no match found


def find_exon_end(query, genome_seq, exon_start, chunk_size=50):
    """
    Function to find exon end position
    """
    # Set starting position to 50
    pos = exon_start + chunk_size
    exon_end = -1

    # Process fragment containing stop codon
    # Compare substring with query
    if Levenshtein.distance(query, genome_seq[exon_start:exon_start + len(query)]) <= 5:
        return exon_start + len(query)

    # Process if stop codon not found
    while True:
        # Search for 'GT' pattern
        gt_index = genome_seq.find('GT', pos)

        if gt_index == -1:
            # Exit loop if 'GT' not found
            break

        # Get substring up to 'GT'
        substring = genome_seq[exon_start:gt_index - 1]

        # Compare substring with query
        if Levenshtein.distance(substring[:-10], query[:len(substring) - 10]) <= 5 \
                and substring[-10:] == query[len(substring) - 10:len(substring)]:
            exon_end = gt_index  # Match found
        else:
            break
        # Set next search position
        pos = gt_index + 2

    return exon_end  # Return no match


def main():
    # Set file paths
    genome_file = 'ref.fasta'
    temp_genome_file = 'temp_ref.fasta'  # Changed temporary filename
    cds_file = 'cds_full.fasta'
    output_file = 'exon_positions.txt'
    chunk_size = 50  # Length of sequence processing unit

    try:
        # Load genome and CDS sequences
        genome_record = next(SeqIO.parse(genome_file, 'fasta'))
        cds_record = next(SeqIO.parse(cds_file, 'fasta'))
        genome_seq = str(genome_record.seq).upper()
        cds_seq = str(cds_record.seq).upper()

        # Determine which genome sequence file to use
        if cds_seq[:50] in genome_seq:
            shutil.copy2(genome_file, temp_genome_file)
        else:
            reverse_complement_sequences(genome_file, temp_genome_file)

        # Load reference genome and CDS sequences
        genome_record = next(SeqIO.parse(temp_genome_file, 'fasta'))
        genome_seq = str(genome_record.seq).upper()

        # Split genome sequence into 50-character chunks
        splitted_genome_seq_list = split_genome_seq(genome_seq, chunk_size)
        ag_splitted_genome_seq_list = [
            (num, string) for num, string in splitted_genome_seq_list
            if string.startswith('AG')
        ]

        # Find Exon1 position and create lists for start position, end position, and exon sequence
        exon_start = [find_exon_start(cds_seq[0:chunk_size], splitted_genome_seq_list)]
        exon_end = [find_exon_end(cds_seq, genome_seq, exon_start[0])]
        exon_length = [exon_end[0] - exon_start[0]]
        exon_seq = [cds_seq[0:exon_length[0]]]

        num_exons = 0
        while not exon_seq[-1].endswith(('TTA', 'TAA', 'TAG')):
            num_exons += 1
            slice_cds_seq = cds_seq[sum(exon_length):]
            exon_start.append(
                find_exon_start('AG' + slice_cds_seq[:chunk_size - 2], ag_splitted_genome_seq_list) + 2
            )
            exon_end.append(find_exon_end(slice_cds_seq, genome_seq, exon_start[-1]))
            exon_length.append(exon_end[-1] - exon_start[-1])
            exon_seq.append(slice_cds_seq[:exon_length[-1]])
            if num_exons > 10:
                exon_seq.append('TTA')

        # Display results
        for i in range(num_exons + 1):
            print(f"Exon {i + 1}", exon_start[i], exon_end[i])

        # Write results to file (tab-separated on one line)
        with open(output_file, 'w') as f:
            for i in range(num_exons + 1):
                f.write(f"Exon{i + 1}\t{exon_start[i]}\t{exon_end[i] - 1}\n")

        print(f"Exon genomic position information has been written to {output_file}")

    finally:
        # Delete temporary file
        if os.path.exists(temp_genome_file):
            os.remove(temp_genome_file)
            print(f"Temporary file {temp_genome_file} has been deleted")


if __name__ == "__main__":
    main()