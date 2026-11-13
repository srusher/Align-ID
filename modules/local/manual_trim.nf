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
	def prefix 	  = "${meta.id}"
	def seqs      = params.man_trim_seqs      ?: 'placeholder'
	def headcrop  = params.man_trim_headcrop  ?: 0
	def tailcrop  = params.man_trim_tailcrop  ?: 0
	def which     = meta.single_end ? 'r1' : (params.man_trim_reads ?: 'r1')
	def r1 = reads[0]
	def r2 = meta.single_end ? 'NA' : reads[1]

	"""

	bash "${projectDir}/bin/manual_trim.sh" $prefix ${meta.single_end} "${projectDir}/bin/manual_trim.py" ${params.man_trim_mode} $headcrop $tailcrop $which $seqs $r1 $r2

	"""

}