process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_low'
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/gatk4/gatk4%3A4.6.1.0--py310hdfd78af_0' :
        'https://depot.galaxyproject.org/singularity/gatk4%3A4.6.1.0--py310hdfd78af_0' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val(meta), path("*.tbi")          , optional:true, emit: tbi
    tuple val(meta), path("*.realigned.bam"), optional:true, emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals $intervals" : ""

    """
    gatk --java-options "-Xmx16g" \\
        CreateSequenceDictionary \\
        -R $fasta

    gatk --java-options "-Xmx16g" \\
        HaplotypeCaller \\
        -R $fasta \\
        -I $bam \\
        -O ${prefix}.vcf.gz \\
        --native-pair-hmm-threads 8 \\
        $interval_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}