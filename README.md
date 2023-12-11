# Biographica interview task

This repo contains the thought processes, code and output for the Bioinformatician interview task at Biographica.

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

`blastp -db db/brassica_db -query data/assembly_2.prot.fa -out output/blastp_raw.txt -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid nident pident length mismatch gapopen qlen slen evalue bitscore" -max_target_seqs 25 -num_threads 8`

This command took about 53 minutes on my local machine. Next apply the 95% identity (identical positions, rather than 
positive-scoring matches) filter:

`awk '$4>=95' output/blastp_raw.txt > output/blastp_results.txt`

*Up to this part of the task has been incorporated into the NextFlow directory.

### Parsing the output and generating the final file

To complete the first task, we just need to call the task1.py script on our input files and BlastP alignment (note that 
the files uploaded to the *output* directory have been zipped for upload):

`python scripts/task1.py data/assembly_1.prot.fa data/assembly_2.prot.fa output/blastp.txt output/final.fa`

This creates the final merged FASTA file, which has been included in the *output* directory.


## Task 2

The second task asks to find proteins with 'hydrolase activity, acting on ester bonds', which is a Gene Ontology term
the maps to GO:0016788. Running task2.py with this GO term and the sequence file generated in the first task looks like
this:

`python scripts/task2.py output/final.fa GO:0016788`

The output of the script indicates that the protein _A01p09130.1_ is the first protein mapped from _assembly_1.prot.fa_
that has hydrolase activity, acting on ester bonds.


## Further considerations and future work

There are a number of improvements that could be made to this process that would result in a better product. 
Documented here are some thoughts about what the next steps would be for writing out this pipeline.

Firstly, ~~the work here can't really be considered a true pipeline as such, as the scripts must be called manually. A 
good immediate next step would be to implement some proper pipeline functionality by wrapping the code with an 
orchestrator such as NextFlow~~ the pipeline should be extended to incorporate the two task scripts; a basic NextFlow 
script has been written around the BLAST portion of the task to demonstrate functionality, but due to the constraints of
time I did not implement the wrapper for the python segment. Moreover, there is a lot of work that could be done to both
scripts (but especially _task2.py_) to generalise them, as right now they are not at all robust and would not cope with 
using any other format than the one present in the two input files. A potential solution would be to have a modular 
system for reading in different types of FASTA file headers.

Secondly, another easy win would be to make the script for the second task store and output EVERY protein that matched
a given GO term - this is easy to code, but my internet access has been poor these last few days and the API calls were
taking a while, so I decided it was best to print the first matching protein. The _task2.py_ script could also be
simplified by removing the input-subsetting function, as this is super specific and it would be cleaner to alter
_task1.py_ to output a file JUST the sequences that were present in both files, rather than all of the sequences from 
both files.

### Pipeline versioning, data governance, metadata logging and scaling

Early on in the README I mention that it would be ideal if this were hosted in the cloud, as that gives easy access to 
scalability. For instance if these scripts were to be run on the cloud, it would be very easy to horizontally scale out 
the BLAST search using compute-specific resources (easily accessed by setting the number of threads in the BLAST call) 
as well as parallelising the API calls.

Similarly, cloud-native resources such as resource monitoring would be a good first step for metadata logging, 
collecting usage statistics, such as the time taken for API responses or the amount of compute time and power used. 
This could also be done at the simplest level using the python _logger_ module and writing usage statistics to a file,
but there are plenty of third-party solutions that would fit the bill (for example CircleCI or Splunk).

Cloud-hosting would also allow for easy handling of data governance, as one could set up security groups specific to 
certain customers (i.e. with AWS's IAM or an implementation of Active Directory), which would prevent inappropriate 
access to commercially sensitive materials. 

A naive pipeline versioning has been implemented here with a git tag, which could be incremented appropriately with 
features, fixes and release work. This would work with a NextFlow-based pipeline as one can pull based on a git tag.