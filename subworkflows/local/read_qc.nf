// short read QC modules
include { FASTQC as FASTQC_RAW                            } from '../../modules/nf-core/fastqc/main'
include { TRIMMOMATIC                                     } from '../../modules/local/trimmomatic'
include { SEQTK_LEN_TRIM                                  } from '../../modules/local/seqtk_len_trim'
include { BBMAP_BBMERGE                                   } from '../../modules/local/bbmap_merge'
include { FASTP                                           } from '../../modules/local/fastp'
include { FASTQC as FASTQC_TRIMMED                        } from '../../modules/nf-core/fastqc/main'

// long read QC modules
include { NANOPLOT as NANOPLOT_RAW                        } from '../../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_TRIMMED                    } from '../../modules/local/nanoplot'
include { PORECHOP_PORECHOP                               } from '../../modules/nf-core/porechop/porechop/main'
include { CHOPPER                                         } from '../../modules/nf-core/chopper/main'
include { CUTADAPT                                        } from '../../modules/nf-core/cutadapt/main'

workflow READ_QC {
    take:
    raw_reads // channel: [ val(meta), [ reads ] ]
	ch_multiqc_files

    main:

	if (params.workflow == 'short-read') {

		FASTQC_RAW (
			raw_reads
		)

		ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))

		if (!params.skip_trimming) {

			if (params.trim_tool == "fastp") {

				if (!params.adapter_auto_detect) {

					FASTP (
						raw_reads,
						["${params.adapt_ref}"],
						[],
						[],
						[]
					)

				} else {

					FASTP (
						raw_reads,
						[],
						[],
						[],
						[]
					)
					
				}

				qc_reads = FASTP.out.reads

				ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] }.ifEmpty([]))
			
			} else if (params.trim_tool == "trimmomatic") {

				TRIMMOMATIC (

					raw_reads

				)

				qc_reads = TRIMMOMATIC.out.trimmed_reads

				ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.trim_log.collect{it[1]}.ifEmpty([]))

			}

			if (!params.skip_merging) {

				BBMAP_BBMERGE (
					qc_reads,
					false
				)

				qc_reads = BBMAP_BBMERGE.out.merged
			}

			// filtering out reads below specified length
			// SEQTK_LEN_TRIM (

			//     qc_reads,
			//     params.seqtk_min_length

			// )
			
			// qc_reads = SEQTK_LEN_TRIM.out.fastq

			FASTQC_TRIMMED (

				qc_reads

			)

			ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))

			// ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

		} else {

			qc_reads = raw_reads

		}
	
	} else {

		NANOPLOT_RAW (

			raw_reads

		)

		ch_multiqc_files.mix(NANOPLOT_RAW.out.txt.collect{it[1]}.ifEmpty([]))

		if (!params.skip_trimming) {

			if (params.skip_porechop) {

				qc_reads = raw_reads
			
			} else {

				PORECHOP_PORECHOP (

					raw_reads    

				)

				qc_reads = PORECHOP_PORECHOP.out.reads

				ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_PORECHOP.out.log.collect().ifEmpty([]))

			}

			CHOPPER (

				qc_reads

			)

			if (!params.skip_cutadapt) {

				CUTADAPT (

					CHOPPER.out.fastq

				)

				qc_reads = CUTADAPT.out.reads

			} else {

				qc_reads = CHOPPER.out.fastq

			}

			NANOPLOT_TRIMMED (

				qc_reads

			)

			ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_TRIMMED.out.txt.collect{it[1]}.ifEmpty([]))  
		
		}

	}

	emit:

	qc_reads = qc_reads
	ch_multiqc_files = ch_multiqc_files

}