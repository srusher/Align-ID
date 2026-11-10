process FILTER_ALIGNMENTS_BY_ID {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    container "${ 'https://depot.galaxyproject.org/singularity/pysam%3A0.24.1--py39hbe927e6_0' }"

    input:
    tuple val(meta), path(bam), path(db)
    path(nodes)
    path(my_tax_ids)
    val(include_children)

    // publishDir "${params.outdir}/tax_id_filtered"

    output:
    tuple val(meta), path('*primary-all-filtered_sorted.bam'), optional:true, emit: tax_filtered_bam  
    tuple val(meta), path('*read_pairs_count.tsv'), optional:true, emit: count_tsv 

    script:
    def prefix = "${meta.id}"
    
    """

    python ${projectDir}/bin/filter_alignments_by_id.py $prefix $bam ${params.my_tax_ids} $include_children ${params.ncbi_taxonomy_nodes} ${projectDir}/bin/collect_child_nodes.py $db

    """
}