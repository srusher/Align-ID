process FILTER_ALIGNMENTS_BY_ID {
    tag "$meta.id"
    label 'process_high'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam), path(db)
    path(seq2tax_map)
    path(my_tax_ids)
    val(include_children)

    publishDir "${params.outdir}/tax_id_filtered"

    output:
    tuple val(meta), path('*primary-all-filtered_sorted.bam') , optional:true, emit: tax_filtered_bam  

    script:
    def prefix = "${meta.id}"
    
    """

    bash "${projectDir}/bin/filter_alignments_by_id.sh" $prefix $bam $seq2tax_map $my_tax_ids $include_children ${params.local_nodes_db} $db

    """
}