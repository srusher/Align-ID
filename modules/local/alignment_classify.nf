process ALIGNMENT_CLASSIFY {
    tag "$meta.id"
    label 'process_high'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(sam_db)


    output:
    tuple val(meta), path('*unmapped-sorted.bam') , optional:true, emit: unmapped_bam
    tuple val(meta), path('*summary.tsv') , optional:true, emit: summary_tsv
    tuple val(meta), path('*primary-all-sorted.bam') , optional:true, emit: primary_all
    tuple val(meta), path('*primary_unambiguous-sorted.bam') , optional:true, emit: primary_unambiguous
    tuple val(meta), path('*primary_ambiguous_single_genome-sorted.bam') , optional:true, emit: primary_ambiguous_single_genome
    tuple val(meta), path('*primary_ambiguous_multi_genome-sorted.bam') , optional:true, emit: primary_ambiguous_multi_genome

    script:
    def prefix = "${meta.id}"

    if (params.workflow == "long-read") {
    
        """

        bash "${projectDir}/bin/alignment_classify_SE.sh" $prefix $bam $sam_db ${params.ncbi_taxonomy_names} ${params.mapping_quality}

        """

    } else if (params.workflow == "short-read") {

        """
        if [ ${meta.single_end} == 'true' ]; then

            bash "${projectDir}/bin/alignment_classify_SE.sh" $prefix $bam $sam_db ${params.ncbi_taxonomy_names} ${params.mapping_quality}

        else

            bash "${projectDir}/bin/alignment_classify_PE.sh" $prefix $bam $sam_db ${params.ncbi_taxonomy_names} ${params.mapping_quality} 

        fi

        """

    }
}