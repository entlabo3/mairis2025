#!/opt/conda/bin/python

import os
import tempfile
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from pathlib import Path


@dataclass
class Variant:
    """Class that holds mutation information"""
    sample: str
    transcript: str
    valid_translation: bool
    cds_pos1: int
    cds_pos2: int
    cds_pos3: int
    ref_codon: str
    sample_codon: str
    ref_pos1: int
    ref_pos2: int
    ref_pos3: int
    ref_aa: str
    sample_aa: str
    mutation_type: str
    variant_type: str


class SequenceAligner:
    """Class that performs sequence alignment using MAFFT"""

    @staticmethod
    def validate_alignment(seq1: str, seq2: str) -> bool:
        """Alignment result validation"""
        if len(seq1) != len(seq2):
            print(f"The sequence lengths do not match.: {len(seq1)} != {len(seq2)}")
            return False
        if len(seq1) == 0 or len(seq2) == 0:
            print("Contains an empty sequence")
            return False
        return True

    @staticmethod
    def parse_mafft_output(mafft_output: str) -> Optional[Tuple[str, str]]:
        """Parses the output of MAFFT and returns a pair of sequences"""
        sequences = []
        current_header = None
        current_seq = []

        for line in mafft_output.split('\n'):
            line = line.strip()
            if not line:  # skip blank line
                continue

            if line.startswith('>'):
                if current_header and current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
                current_header = line
            else:
                current_seq.append(line)

        # 最後の配列を追加
        if current_header and current_seq:
            sequences.append(''.join(current_seq))

        if len(sequences) != 2:
            print(f"Unexpected number of sequences: {len(sequences)}")
            return None

        if sequences:
            print(f"sequence length: {[len(seq) for seq in sequences]}")

        if not SequenceAligner.validate_alignment(*sequences):
            return None

        return tuple(sequences)

    @staticmethod
    def align_sequences(ref_seq: str, sample_seq: str) -> Tuple[str, str]:
        """Align the two sequences using MAFFT"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_file:
            temp_file.write(f">ref\n{ref_seq}\n>sample\n{sample_seq}\n")
            temp_fasta = temp_file.name

        try:
            cmd = ['mafft', '--quiet', '--auto', '--thread', '1', temp_fasta]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            aligned_seqs = SequenceAligner.parse_mafft_output(result.stdout)

            if aligned_seqs is None:
                raise ValueError("The alignment results are incorrect.")

            return aligned_seqs

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while executing MAFFT: {e}")
            print(f"MAFFT stderr: {e.stderr}")
            raise
        finally:
            if os.path.exists(temp_fasta):
                os.unlink(temp_fasta)


class VariantAnalyzer:
    """Class for performing mutation analysis"""

    def __init__(self, ref_fasta: str):
        self.references = self._load_references(ref_fasta)
        self.aligner = SequenceAligner()
        self.genetic_code = {
            'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
            'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
            'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Ter', 'TAG': 'Ter',
            'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Ter', 'TGG': 'Trp',
            'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
            'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
            'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
            'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
            'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
            'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
            'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
            'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
            'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
            'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
            'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
            'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
        }

    def _load_references(self, ref_fasta: str) -> Dict[str, str]:
        """Load reference FASTA file"""
        references = {}
        for record in SeqIO.parse(ref_fasta, 'fasta'):
            references[record.id] = str(record.seq)
        return references

    def get_unaligned_position(self, aligned_seq: str, aligned_pos: int) -> int:
        """Calculate the position in the original array from the alignment position"""
        if aligned_pos < 0 or aligned_pos >= len(aligned_seq):
            return 0

        unaligned_pos = 0
        for i in range(aligned_pos):
            if aligned_seq[i] != '-':
                unaligned_pos += 1
        return unaligned_pos

    def is_valid_codon(self, codon: str) -> bool:
        """Check whether the codon can be translated into an effective amino acid."""
        if '-' in codon or len(codon) != 3:
            return False
        return codon.upper() in self.genetic_code

    def get_amino_acid(self, codon: str) -> str:
        """Conversion from codon to amino acid"""
        if not self.is_valid_codon(codon):
            return '---'
        return self.genetic_code.get(codon.upper(), 'Unk')

    def determine_mutation_type(self, ref_codon: str, sample_codon: str) -> str:
        """Determine the type of mutation (synonymous/non-synonymous)."""
        ref_aa = self.get_amino_acid(ref_codon)
        sample_aa = self.get_amino_acid(sample_codon)

        if ref_aa == '---' or sample_aa == '---' or ref_aa == 'Unk' or sample_aa == 'Unk':
            return 'Unknown'

        return 'Synonymous' if ref_aa == sample_aa else 'Nonsynonymous'

    def determine_variant_type(self, ref_codon: str, sample_codon: str) -> str:
        """Determine the mutation type (SNP/MNP/INS/DEL)"""
        if '-' in ref_codon:
            return 'INS'
        elif '-' in sample_codon:
            return 'DEL'

        diff_count = sum(1 for r, s in zip(ref_codon, sample_codon) if r != s)
        return 'SNP' if diff_count == 1 else 'MNP'

    def analyze_variants(self, sample_file: str) -> List[Variant]:
        """Mutations in the sample file"""
        variants = []
        sample_sequences = {}

        try:
            for record in SeqIO.parse(sample_file, 'fasta'):
                sample_sequences[record.id] = str(record.seq)
        except Exception as e:
            print(f"Error reading FASTA file: {sample_file}")
            print(f"Error details: {e}")
            return variants

        sample_name = os.path.splitext(os.path.basename(sample_file))[0]
        if sample_name.endswith('_transcripts'):
            sample_name = sample_name[:-12]

        for sample_id, sample_seq in sample_sequences.items():
            try:
                transcript_id = sample_id.split('_')[-1]
                if transcript_id not in self.references:
                    print(f"Transcripts that do not exist in the reference: {transcript_id}")
                    continue

                ref_seq = self.references[transcript_id]
                aligned_ref, aligned_sample = self.aligner.align_sequences(ref_seq, sample_seq)

                if not self.aligner.validate_alignment(aligned_ref, aligned_sample):
                    print(f"The alignment results are incorrect.: {sample_id}")
                    continue

                for i in range(0, len(aligned_ref) - 2, 3):
                    if i + 3 > len(aligned_ref) or i + 3 > len(aligned_sample):
                        break

                    ref_codon = aligned_ref[i:i + 3]
                    sample_codon = aligned_sample[i:i + 3]

                    if len(ref_codon) != 3 or len(sample_codon) != 3:
                        continue

                    if ref_codon != sample_codon:
                        unaligned_pos1 = self.get_unaligned_position(aligned_sample, i) + 1
                        unaligned_pos2 = self.get_unaligned_position(aligned_sample, i + 1) + 1
                        unaligned_pos3 = self.get_unaligned_position(aligned_sample, i + 2) + 1

                        ref_pos1 = self.get_unaligned_position(aligned_ref, i) + 1
                        ref_pos2 = self.get_unaligned_position(aligned_ref, i + 1) + 1
                        ref_pos3 = self.get_unaligned_position(aligned_ref, i + 2) + 1

                        clean_ref_codon = ref_codon.replace('-', '')
                        clean_sample_codon = sample_codon.replace('-', '')

                        variants.append(Variant(
                            sample=sample_name,
                            transcript=transcript_id,
                            valid_translation=self.is_valid_codon(clean_sample_codon),
                            cds_pos1=unaligned_pos1,
                            cds_pos2=unaligned_pos2,
                            cds_pos3=unaligned_pos3,
                            ref_codon=clean_ref_codon,
                            sample_codon=clean_sample_codon,
                            ref_pos1=ref_pos1,
                            ref_pos2=ref_pos2,
                            ref_pos3=ref_pos3,
                            ref_aa=self.get_amino_acid(clean_ref_codon),
                            sample_aa=self.get_amino_acid(clean_sample_codon),
                            mutation_type=self.determine_mutation_type(clean_ref_codon, clean_sample_codon),
                            variant_type=self.determine_variant_type(ref_codon, sample_codon)
                        ))

            except Exception as e:
                print(f"Error occurred: {sample_id}")
                print(f"Error details: {e}")
                continue

        return variants


def process_sample_file(args) -> List[Variant]:
    """Function to process a single sample file (for parallel processing)"""
    sample_file, ref_fasta = args
    analyzer = VariantAnalyzer(ref_fasta)
    return analyzer.analyze_variants(sample_file)


def clean_output_directory(output_dir):
    """Clean up the output directory"""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)


def main():
    # Parameter Settings
    output_dir = "./0610_list_variants"
    input_dir = "./0430_transcripts"
    ref_fasta = "./ref/ref_trimmed_cds.fasta"

    try:
        # Creating the output directory
        clean_output_directory(output_dir)

        # Create a list of input files
        sample_files = [(f, ref_fasta) for f in Path(input_dir).glob("*.fasta")]
        total_files = len(sample_files)

        if total_files == 0:
            print("The file to be processed cannot be found.")
            return

        print(f"Number of files to be processed: {total_files}")

        # プロセス数の設定（CPU数の90%を使用）
        num_processes = max(1, int(multiprocessing.cpu_count() * 0.90))
        print(f"Number of processes used: {num_processes}")

        # ProcessPoolExecutorで並列処理
        variants = []
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            futures = [
                executor.submit(process_sample_file, args)
                for args in sample_files
            ]

            # 進捗表示用のカウンター
            completed = 0

            for future in as_completed(futures):
                completed += 1
                print(f"Progress: {completed}/{total_files} files Finished")

                try:
                    sample_variants = future.result()
                    variants.extend(sample_variants)
                except Exception as e:
                    print(f"An error occurred while processing the file.: {str(e)}")
                    continue

        # Output file path
        output_file = os.path.join(output_dir, "variants.tsv")

        # 結果の保存
        with open(output_file, 'w') as f:
            # ヘッダーの書き込み
            f.write("Sample\tTranscript\tValid_translation\t"
                    "Ref_CDS_pos1\tRef_CDS_pos2\tRef_CDS_pos3\t"
                    "CDS_pos1\tCDS_pos2\tCDS_pos3\t"
                    "Ref_codon\tSample_codon\t"
                    "Ref_AA\tSample_AA\t"
                    "Mutation_type\tVariant_type\n")

            # Writing mutation data
            for variant in variants:
                f.write(f"{variant.sample}\t{variant.transcript}\t{variant.valid_translation}\t"
                        f"{variant.ref_pos1}\t{variant.ref_pos2}\t{variant.ref_pos3}\t"
                        f"{variant.cds_pos1}\t{variant.cds_pos2}\t{variant.cds_pos3}\t"
                        f"{variant.ref_codon}\t{variant.sample_codon}\t"
                        f"{variant.ref_aa}\t{variant.sample_aa}\t"
                        f"{variant.mutation_type}\t{variant.variant_type}\n")

            print("\nThe processing is complete.")
            print(f"Total {len(variants)} mutations detected")
            print(f"The results were saved to {output_file}.")

    except Exception as e:
        print(f"An error occurred during execution: {str(e)}")
        raise


if __name__ == "__main__":
    main()
