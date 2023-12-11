// Pipeline module that takes a pair of protein FASTA files, generates a database from the first and runs BlastP
//  against it with the second. The final output of the module can be used by downstream processes.


// Using a separate config file allows us to generalise the BlastP portion of the pipeline
params.configFile = file("protein_blast.config")

// Parse the configuration file
blastConfig = readConfig(params.configFile)


// Generate a database using the given input FASTA file
process makeBlastDB {
    output:
    file "${blastConfig.makeblastdb.output_dir}/${blastConfig.makeblastdb.output}"

    script:
    """
    makeblastdb -in ${blastConfig.makeblastdb.input} -dbtype ${blastConfig.makeblastdb.db_type} -title "${blastConfig.makeblastdb.db_title}" -out ${blastConfig.makeblastdb.output_dir}/${blastConfig.makeblastdb.output}
    """
}

// Run the protein blast of our second FASTA file against the database generated in the first step
process blastP {
    input:
    path "${blastConfig.makeblastdb.output_dir}/${blastConfig.makeblastdb.output}"

    output:
    file "${blastConfig.blastp.output}"

    script:
    """
    blastp -db ${blastDB} -query ${blastConfig.blastp.input} -out ${blastConfig.blastp.output} -qcov_hsp_perc ${blastConfig.blastp.qcov_hsp_perc} -outfmt "${blastConfig.blastp.outfmt}" -max_target_seqs ${blastConfig.blastp.max_target_seqs} -num_threads ${blastConfig.blastp.num_threads}
    """
}

// Filter the output with our specified percentage sequence identity
process filterBlastResults {
    input:
    file "${blastConfig.blastp.output}"

    output:
    file "${blastConfig.awk_filter.output}"

    script:
    """
    awk '$4>=${blastConfig.awk_filter.pident}' ${blastp_raw} > ${blastConfig.awk_filter.output}
    """
}


// Define the Nextflow workflow
workflow blastWorkflow {
    // Execute the steps
    makeBlastDB
    blastP(makeBlastDB.out)
    filterBlastResults(blastP.out)
}
