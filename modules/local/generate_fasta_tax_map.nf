process GENERATE_FASTA_TAX_MAP {
    label 'process_high'

    container "${ 'https://depot.galaxyproject.org/singularity/pysam%3A0.24.1--py39hbe927e6_0' }"

    output:
	path('*tax_map.tsv'), emit: tax_map

    script:
    
    """

    python ${projectDir}/bin/generate_fasta_tax_map.py ${params.ref_fasta} ${projectDir}/assets/taxonomy/names_scientific.dmp ${projectDir}/assets/taxonomy/names_other.dmp ${projectDir}/assets/taxonomy/nucl_gb_seqid2taxid.map

    """
}