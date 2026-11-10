process GENERATE_TAXONOMY_GRAPH {

    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(alignment_report), path(tax_report_tsv), path(read_count_tsv)
    val(data_type)

    output:
    tuple val(meta), path('*.png') , optional:true, emit: plot

    script:
    def prefix = "${meta.id}"

    """
    if [[ "${data_type}" == "Reads" && "${meta.single_end}" == 'false' ]]; then

        data_type='Read Pairs'

    else

        data_type=${data_type}

    fi

    if [[ "${params.skip_filter_alignment_by_id}" == 'false' ]]; then

        total_reads=\$(head -1 $read_count_tsv)

        bash "${projectDir}/bin/generate_taxonomy_graph.sh" $prefix $alignment_report "${projectDir}/bin/alignment_plot_plus_coverage.R" "\$data_type" "${params.skip_filter_alignment_by_id}" \$total_reads $tax_report_tsv
    
    else

        bash "${projectDir}/bin/generate_taxonomy_graph.sh" $prefix $alignment_report "${projectDir}/bin/alignment_plot_plus_coverage.R" "\$data_type" "${params.skip_filter_alignment_by_id}" 1 $tax_report_tsv

    fi

    """

}