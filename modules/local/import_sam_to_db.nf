process IMPORT_SAM_TO_DB {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    container "${ 'https://depot.galaxyproject.org/singularity/pysam%3A0.24.1--py39hbe927e6_0' }"

    input:
    tuple val(meta), path(bam)
	path(tax_map)


    output:
    tuple val(meta), path('*sam.db') , optional:true, emit: sam_db 

    script:
    def prefix = "${meta.id}"
    // def seqid2taxid_map = "${params.seqid2taxid_map}"

    // if (!params.seqid2taxid_map.startsWith('/')) {

    //     seqid2taxid_map = "${projectDir}/${params.seqid2taxid_map}"

    // }

    """

    python "${projectDir}/bin/import_sam_to_db.py" $prefix $bam $tax_map "${meta.single_end}" "${params.include_primary_seq}"

    """

    
}