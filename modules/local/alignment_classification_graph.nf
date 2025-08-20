process ALIGNMENT_CLASSIFICATION_GRAPH {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(alignment_report)
    val(data_type)

    output:
    tuple val(meta), path('*.png') , optional:true, emit: plot

    script:
    def prefix = "${meta.id}"

    """

    bash "${projectDir}/bin/alignment_classification_graph.sh" $prefix $alignment_report "${projectDir}/bin/alignment_plot.R" ${params.mapping_quality} ${data_type}

    """

}