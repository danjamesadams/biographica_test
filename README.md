# Biographica interview task

This repo contains the thought processes, code and output for the Bioinformatician role at Biographica.

An initial inspection of the data shows that *assembly_1.prot.fa* contains 108190 sequences and *assembly_2.prot.fa* 
contains 61153 sequences, meaning that not every sequence in the first file will have a direct map in the second. Also
worth nothing is the fact that there may well be sequences in the larger, first file which do not map to the second, as
there is nothing to indicate that *assembly_2.prot.fa* is a subset of *assembly_1.prot.fa*.


## Task 1

In order to disambiguate the proteins, we need to perform a complete set of pairwise alignments between the two FASTA 
files. This is most easily done with the NCBI BLAST tool, which can be run locally (ideally in a production system this
would be done hosted in the cloud).

### Creating a local database to search against

Using the larger of the two assemblies, we can create a local database to search against with the smaller assembly. The
following command creates a database under the *db* directory:

`makeblastdb -in data/assembly_1.prot.fa -dbtype prot -title "Brassica napus db" -out db/brassica_db`

### Searching against the database

Using the parameters specified in the task (coverage of 90%), we can run a protein BLAST against our database:

`blastp -db db/brassica_db -query data/assembly_2.prot.fa -out output/blastp_raw.txt -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid nident pident length mismatch gapopen qlen slen evalue" -max_target_seqs 25 -num_threads 8`

Then apply the 95% identity (identical positions, rather than positive-scoring matches) filter:

`awk '$4>=95' output/blastp_raw.txt > output/blastp_results.txt`
