process MANUAL_TRIM {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_manual_trim*.fastq.gz")   , emit: trimmed_reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
	def seqs = "${params.man_trim_seqs}" ?: "placeholder"
	def headcrop = "${params.man_trim_headcrop}" ?: "0"
	def tailcrop = "${params.man_trim_tailcrop}" ?: "0"
	def which_reads = "${params.man_trim_reads}" ?: "r1"

	"""

	bash "${projectDir}/bin/manual_trim.sh" $prefix "${meta.single_end}" "${projectDir}/bin/manual_trim.py" ${params.man_trim_mode} $headcrop $tailcrop $which_reads $seqs $reads

	"""

}