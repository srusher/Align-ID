process BOWTIE2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    errorStrategy 'ignore'

    //Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2%3A2.5.4--h7071971_4' :
        'https://depot.galaxyproject.org/singularity/bowtie2%3A2.5.4--h7071971_4' }"

    input:
    tuple val(meta), path(reads)
    val(index)

    output:
    tuple val(meta), path("*.sam"), optional: true, emit: sam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if ${params.memory_saver}; then
    
        start="false"

        while [[ \$start == "false" ]]; do

            if [[ \$(ls ${projectDir}/queue/bowtie2 | wc -l) -gt 1 ]]; then

                sleep 5
            
            else

                start="true"
                touch ${projectDir}/queue/bowtie2/$prefix-bowtie2

            fi

        done
    fi

    bowtie2 \\
        -x $index \\
        -U $reads \\
        -S ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1)
    END_VERSIONS

    if ${params.memory_saver}; then

        rm -f ${projectDir}/queue/bowtie2/$prefix-bowtie2
    
    fi

    """


}