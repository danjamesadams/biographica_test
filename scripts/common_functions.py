from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_fasta_file(fasta_file: str) -> dict[SeqRecord]:
    """Read in a protein FASTA file and parse the sequences into Bio.SeqRecord.SeqRecord objects.

    :param fasta_file: Path to a protein FASTA file
    :return: A dictionary of Fasta SeqRecord objects
    """
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
