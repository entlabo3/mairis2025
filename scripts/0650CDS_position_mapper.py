#!/opt/conda/bin/python

import os
import subprocess
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import IUPACData


def one_to_three(aa):
    """
    Convert 1-letter amino acid code to 3-letter code

    Args:
        aa (str): 1-letter amino acid code

    Returns:
        str: 3-letter amino acid code
    """
    if aa is None or aa == "":
        return ""
    return IUPACData.protein_letters_1to3[aa].capitalize()


def translate_cds_to_protein(cds_seq):
    """
    Function to translate CDS to amino acid sequence

    Args:
        cds_seq (str): CDS sequence to translate

    Returns:
        str: Amino acid sequence
    """
    return str(Seq(cds_seq).translate())


def run_mafft_pairwise(ref_aa_file, query_aa_file, output_file):
    """
    Function to run pairwise alignment using MAFFT

    Args:
        ref_aa_file (str): Reference amino acid sequence file
        query_aa_file (str): Query amino acid sequence file
        output_file (str): Output file for alignment results
    """
    # Create input file
    combined_input = output_file + ".in"
    with open(combined_input, "w") as outf:
        # Copy reference sequence
        with open(ref_aa_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                SeqIO.write(record, outf, "fasta")
        # Copy query sequence
        with open(query_aa_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                SeqIO.write(record, outf, "fasta")

    # Execute MAFFT command
    cmd = ["mafft", "--quiet", "--auto", combined_input]

    try:
        # Write MAFFT execution results directly to output file
        with open(output_file, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)

        # Delete temporary file
        os.remove(combined_input)

    except subprocess.CalledProcessError as e:
        print(f"MAFFT execution error: {e.stderr.decode()}")
        raise


def create_position_mapping(alignment_file):
    """
    Function to create position mapping information from alignment file

    Args:
        alignment_file (str): MAFFT alignment result file

    Returns:
        dict: Mapping dictionary from CDS position to refAA position
        {cds_pos: (refaa_pos, aa)}
    """
    records = list(SeqIO.parse(alignment_file, "fasta"))
    ref_seq = str(records[0].seq)
    query_seq = str(records[1].seq)

    mapping = {}
    ref_pos = 0
    query_pos = 0

    for i, (ref_aa, query_aa) in enumerate(zip(ref_seq, query_seq)):
        if ref_aa != "-":
            ref_pos += 1
        if query_aa != "-":
            query_pos += 1
            if ref_aa != "-":
                mapping[query_pos] = (ref_pos, ref_aa)
            else:
                mapping[query_pos] = (None, None)

    return mapping


def process_cds_sequences(cds_fasta, ref_aa, temp_dir="./temp_alignment"):
    """
    Process CDS sequences and create mapping information for each transcript

    Args:
        cds_fasta (str): CDS FASTA file
        ref_aa (str): Reference amino acid sequence file
        temp_dir (str): Directory for temporary files

    Returns:
        dict: Mapping information for each transcript
    """
    # Prepare temporary directory
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)

    transcript_mappings = {}

    for record in SeqIO.parse(cds_fasta, "fasta"):
        # Translate CDS to amino acid sequence
        protein_seq = translate_cds_to_protein(str(record.seq))
        query_aa_file = os.path.join(temp_dir, f"{record.id}_aa.fasta")

        # Save translated amino acid sequence to temporary file
        with open(query_aa_file, "w") as f:
            SeqIO.write(
                SeqRecord(Seq(protein_seq), id=record.id, description=""), f, "fasta"
            )

        # Align with MAFFT
        aln_file = os.path.join(temp_dir, f"{record.id}_aln.fasta")
        run_mafft_pairwise(ref_aa, query_aa_file, aln_file)

        # Create mapping information
        transcript_mappings[record.id] = create_position_mapping(aln_file)

    return transcript_mappings


def main():
    # Set input/output files
    gff3_file = "./ref/ref_trimmed.gff3"
    cds_fasta_file = "./ref/ref_trimmed_cds.fasta"
    refAA_fasta_file = "./ref/RefAA.fasta"
    input_folder = "./0610_list_variants"
    input_file = os.path.join(input_folder, "variants_snp_mnp_processed_with_genome_pos_merged.tsv")
    output_file = os.path.join(input_folder, "variants_with_refaa_pos.tsv")
    temp_dir = os.path.join(input_folder, "temp_alignment")

    # Create mapping information
    transcript_mappings = process_cds_sequences(
        cds_fasta_file, refAA_fasta_file, temp_dir
    )

    # Process input file
    with open(input_file) as f_in, open(output_file, "w") as f_out:
        # Process header
        header = f_in.readline().strip()
        f_out.write(f"{header}\tMdAA_pos\tMd_AA\n")

        # Process each line
        for line in f_in:
            fields = line.strip().split("\t")
            transcript = fields[1]
            cds_pos = int(fields[3])  # Use Ref_CDS_pos1

            # Calculate amino acid position (1-based)
            aa_pos = (cds_pos - 1) // 3 + 1

            # Get mapping information
            if transcript == "Common":
                # Use information from any transcript where mapping is found
                for mapping in transcript_mappings.values():
                    if aa_pos in mapping:
                        refaa_pos, refaa = mapping[aa_pos]
                        if refaa_pos is not None:
                            fields.extend([str(refaa_pos), one_to_three(refaa)])
                            break
                else:
                    fields.extend(["", ""])
            else:
                # When specific transcript is specified
                transcript_key = None
                for key in transcript_mappings.keys():
                    if transcript in key:  # Search transcript by partial match
                        transcript_key = key
                        break

                if transcript_key and aa_pos in transcript_mappings[transcript_key]:
                    refaa_pos, refaa = transcript_mappings[transcript_key][aa_pos]
                    fields.extend(
                        [
                            str(refaa_pos) if refaa_pos is not None else "",
                            one_to_three(refaa) if refaa is not None else "",
                        ]
                    )
                else:
                    fields.extend(["", ""])

            f_out.write("\t".join(fields) + "\n")
    print(f"Processing completed. Output file: {output_file}")

    # Delete temporary directory
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


if __name__ == "__main__":
    main()