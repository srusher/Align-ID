process ALIGNMENT_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    container "${ 'https://depot.galaxyproject.org/singularity/pysam%3A0.24.1--py39hbe927e6_0' }"
        
    input:
    tuple val(meta), path(bam), path(sam_db)


    output:
    tuple val(meta), path('*unmapped-sorted.bam') , optional:true, emit: unmapped_bam
    tuple val(meta), path('*summary.tsv') , optional:true, emit: summary_tsv
    tuple val(meta), path('*primary_all_sorted.bam') , optional:true, emit: primary_all
    tuple val(meta), path('*primary_unambiguous_sorted.bam') , optional:true, emit: primary_unambiguous
    tuple val(meta), path('*primary_ambiguous_single_genome_sorted.bam') , optional:true, emit: primary_ambiguous_single_genome
    tuple val(meta), path('*primary_ambiguous_multi_genome_sorted.bam') , optional:true, emit: primary_ambiguous_multi_genome

    script:
    def prefix = "${meta.id}"

    if (params.workflow == "long-read") {
    
        """

        python "${projectDir}/bin/alignment_classify_SE.py" $prefix $bam $sam_db ${projectDir}/${params.ncbi_taxonomy_names} ${params.mapping_quality}

        """

    } else if (params.workflow == "short-read") {

        """
        if [ ${meta.single_end} == 'true' ]; then

            python "${projectDir}/bin/alignment_classify_SE.py" $prefix $bam $sam_db ${projectDir}/${params.ncbi_taxonomy_names} ${params.mapping_quality}

        else

            python "${projectDir}/bin/alignment_classify_PE.py" $prefix $bam $sam_db ${projectDir}/${params.ncbi_taxonomy_names} ${params.mapping_quality} 

        fi

        """

    }
}