process GENERATE_TAXONOMY_TSV {
	
    tag "$meta.id"
    label 'process_low'
    errorStrategy 'ignore'

    container "${ 'https://depot.galaxyproject.org/singularity/pysam%3A0.24.1--py39hbe927e6_0' }"

    input:
    tuple val(meta), path(tax_summary), path(coverage), path(read_count)
	path(tax_map)

    // publishDir "${params.outdir}/tax_id_filtered"

    output:
    tuple val(meta), path('*report.tsv'), optional:true, emit: tax_report_tsv

    script:
    def prefix = "${meta.id}"
    
    """

    python ${projectDir}/bin/generate_tax_report.py $prefix $tax_summary $coverage $tax_map $read_count

    """
}