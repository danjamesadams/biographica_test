import argparse

import requests
from Bio.SeqRecord import SeqRecord

from common_functions import parse_fasta_file


class UniProtMapper:
    """Given a FASTA file where at least some sequence headers contain UniProt protein sequence IDs, map those
    identifiers to entries in UniProt. Afterwards, those mapped proteins can be interrogated for passed features to
    select proteins with a specific function.
    """

    UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self, fasta_file: str, go_id: str) -> None:
        """Instantiate the class with a protein FASTA file.

        :param fasta_file: Path to a protein FASTA file where at least some headers contain UniProt protein sequence IDs
        :param go_id: Gene Ontology ID denoting the function that we want to find in our protein
        """
        self.input = parse_fasta_file(fasta_file)
        self.go_id = go_id

    # TODO: Update task 1 to store the header lines in a more parser-friendly manner
    def select_records_that_were_in_both_assemblies(self) -> dict[SeqRecord]:
        """Loop through an input FASTA file and keep only those records that were present in both input assembly FASTA
        files (meaning that they should all have UniProt accession IDs).

        :return: A subset dictionary of FASTA SeqRecords, where each entry has a putative UniProt accession ID
        """
        subset_dict = {}
        # If we have references to both assemblies in the ID, we know we have a mapped protein with a UniProt ID
        for sequence in self.input:
            record = self.input[sequence]
            # Due to the way the headers were saved in task 1, the multi-assembly lines only appear in the description
            if "assembly_1" in record.description and "assembly_2" in record.description:
                # Extract the UniProt ID for easier access later (annotations are empty here so this is safe)
                # FIXME: Naively assumes that the second record is the second assembly (which it is in this case but
                #  might not always be)
                record.annotations = record.description.split(",")[1].split("=")[1]
                subset_dict[sequence] = record

        return subset_dict

    def retrieve_uniprot_entries_and_pull_out_by_function(self, sequence_records: dict[SeqRecord]) -> str:
        """Query UniProt with the provided protein IDs and return a list of proteins that match the desired GO term.

        :param sequence_records: A dictionary of Biopython SeqRecords denoting protein sequences
        :return: A list of proteins if are found with the sought GO term, else a statement saying that none were found
        """
        proteins_with_desired_function = []
        for sequence in sequence_records:
            # We stored our protein IDs in the annotations attribute in the previous function
            protein_id = sequence_records[sequence].annotations
            # Attempt to retrieve a UniProt accession by searching with a protein ID
            response = requests.get(self.UNIPROT_BASE_URL + "search?fields=id&format=json&query={}".format(protein_id))
            response.raise_for_status()

            # There are quite a few things that could go wrong here, so wrap them all in a try-block
            try:
                # Interrogate the JSON output to extract the UniProt accession
                uniprot_accession = response.json()["results"][0]["primaryAccession"]
                uniprot_entry = requests.get(self.UNIPROT_BASE_URL + uniprot_accession)
                uniprot_entry.raise_for_status()

                for cross_reference in uniprot_entry.json()["uniProtKBCrossReferences"]:
                    # Check to see if we have any Gene Ontology entries for our UniProt entry
                    if cross_reference["database"] == "GO":
                        if cross_reference["id"] == self.go_id:
                            print(sequence_records[sequence].description)
                            return

            # If anything went amiss, raise it as an exception
            except requests.RequestException as e:
                return f"The following error occurred whilst trying to search UniProt: {str(e)}"

        # Print a message to the console to say that nothing matching was found
        return f"No proteins were found with the GO term {self.go_id}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify proteins with a specific GO term from a provided FASTA file")

    parser.add_argument("sequence_file", help="Path to the first, larger protein FASTA file")
    parser.add_argument("go_term", help="Path to the second, smaller protein FASTA file")
    args = parser.parse_args()

    mapper = UniProtMapper(args.sequence_file, args.go_term)
    sequences_in_both_assemblies = mapper.select_records_that_were_in_both_assemblies()
    mapper.retrieve_uniprot_entries_and_pull_out_by_function(sequences_in_both_assemblies)
