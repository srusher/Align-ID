process IMPORT_SAM_TO_DB {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam)
    path(seq2tax_map)


    output:
    tuple val(meta), path('*sam.db') , optional:true, emit: sam_db 

    script:
    def prefix = "${meta.id}"

    """

    bash "${projectDir}/bin/import_sam_to_db.sh" $prefix $bam $seq2tax_map "${meta.single_end}"

    """

    
}