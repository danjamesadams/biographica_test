from Bio import SeqIO
import pandas as pd


# Column labels for the BlastP alignment output
HEADER = ["qseqid", "sseqid", "nident", "pident", "length", "mismatch", "gapopen", "qlen", "slen", "evalue", "bitscore"]


class ProteinBLASTMapper:
    """Given two protein FASTA files and the output of their BlastP alignment, select the optimum matches for each
    unique sequence, join up their header lines from the original FASTA files and ultimately create a single file
    containing ALL unique sequences only once.

    It should be noted that this class is EXTREMELY tied to the input files and could do with some tweaks to make it
    generalisable for different input FASTA files.
    """

    def __init__(self, assembly1: str, assembly2: str, alignment: str) -> None:
        """Instantiate the class by reading in the sequence files and filtered BlastP output.

        :param assembly1: Path to the first of the two protein FASTA files
        :param assembly2: Path to the second of the two protein FASTA files
        :param alignment: Path to the BlastP output of the two assemblies
        """
        self.assembly1 = assembly1
        self.assembly2 = assembly2
        self.alignment = alignment

    @staticmethod
    def select_group_by_evalue_then_bitscore(df: pd.DataFrame, group: str) -> pd.DataFrame:
        """Sort the protein BLAST output by a group (either query or subject sequence ID), then select the alignment
        with the highest E-value (and bit-score if the E-values are tied).

        :param df: A pandas Dataframe containing BlastP matches
        :param group: The name of the column on which to call df.groupby()
        :return: The Dataframe row with the highest E-value and bit-score for a given group
        """
        return df.sort_values(by=[group, "evalue", "bitscore"], ascending=[True, True, False]).groupby(group).head(1)

    def parse_protein_blast_results_and_select_best_matches(self) -> pd.DataFrame:
        """Read in the BlastP output file and select the match with the highest bit score.

        This function works based on the fact that we have already selected matches with a coverage greater than or
        equal to 90% and sequence identity of greater than or equal to 95%. Where there are multiple matches, the one
        with the lowest E-value is our best match. If the E-values are equal (i.e. both 0 due to rounding), then we
        should take the alignment with the highest bit-score - this ensures that our final dataframe contains entirely
        unique query and subject sequences.

        :return: A pandas Dataframe containing the 'best' matches from the protein BLAST
        """
        # Read the BlastP results into a pandas dataframe (more efficient at handling larger datasets)
        blastp_output = pd.read_csv(self.alignment, sep="\t")
        blastp_output.columns = HEADER

        # For each unique query sequence, take the protein alignment with the lowest E-value (and bit-score if tied)
        matching_proteins = self.select_group_by_evalue_then_bitscore(blastp_output, "qseqid")

        # If we have duplicate matches for any subject sequences, we need to take their 'best' match
        duplicate_matches = matching_proteins[matching_proteins.duplicated("sseqid", keep=False)]
        deduplicated_matches = self.select_group_by_evalue_then_bitscore(duplicate_matches, "sseqid")

        # Pull out all unique matches for rejoining with the subset of deduplicated matches
        unique_matches = matching_proteins.drop_duplicates(subset="sseqid", keep=False)
        best_matches = pd.concat([unique_matches, deduplicated_matches])

        return best_matches

    def write_output_fasta_file(self, aligned_sequences: pd.DataFrame) -> None:
        """Create a protein FASTA file with each unique sequence (counting those that successfully mapped according to
        BlastP as a single sequence per pair) recorded once. Those mapped pairs are written with concatenated headers.

        :param aligned_sequences: A pandas Dataframe containing the best alignments from a BlastP output
        """
        # Read the sequences from the original FASTA files
        assembly1_dict = SeqIO.to_dict(SeqIO.parse(self.assembly1, "fasta"))
        assembly2_dict = SeqIO.to_dict(SeqIO.parse(self.assembly2, "fasta"))

        # Get the identifiers of sequences that are matches in the 'best_matches' DataFrame
        assembly1_mapped_sequences = aligned_sequences["sseqid"].tolist()
        assembly2_mapped_sequences = aligned_sequences["qseqid"].tolist()

        # Open the output FASTA file for writing the sequences
        with open("output/final.fa", "wt") as final_fasta:
            # Start with the larger, first assembly (our database/subject proteins)
            for sequence in assembly1_dict:
                # Store the ID in a variable
                sequence_id = assembly1_dict[sequence].id

                # Check to see if the identifier was one of the mapped proteins
                if sequence_id in assembly1_mapped_sequences:
                    # Pull out the corresponding sequence ID from the alternative assembly and concatenate the two
                    alt_sequence_id = aligned_sequences[aligned_sequences["sseqid"] == sequence_id].iloc[0]["qseqid"]
                    sequence_id = f"assembly_1.prot.fa={sequence_id}, assembly_2.prot.fa={alt_sequence_id}"
                else:
                    sequence_id = f"assembly_1.prot.fa={sequence_id}"

                # Write the record to the FASTA file
                record = SeqIO.SeqRecord(assembly1_dict[sequence].seq, id=sequence_id, description="")
                SeqIO.write(record, final_fasta, "fasta")

            # Move onto the smaller, second assembly (our query proteins)
            for sequence in assembly2_dict:
                # Store the ID in a variable
                sequence_id = assembly2_dict[sequence].id

                # As we have already written out the mapped proteins, skip any matches in the query assembly
                if sequence_id in assembly2_mapped_sequences:
                    continue

                # If we did not map the protein, add it to the final FASTA file
                sequence_id = f"assembly_2.prot.fa={sequence_id}"
                record = SeqIO.SeqRecord(assembly2_dict[sequence].seq, id=sequence_id, description="")
                SeqIO.write(record, final_fasta, "fasta")


if __name__ == "__main__":
    # TODO: Add argparse functionality to call script from the command line
    pass
