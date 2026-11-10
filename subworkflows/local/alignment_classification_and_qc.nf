// global modules
include { IMPORT_SAM_TO_DB														 } from '../../modules/local/import_sam_to_db'
include { ALIGNMENT_CLASSIFY                              						 } from '../../modules/local/alignment_classify'
include { BLAST_UNMAPPED_READS                            					     } from '../../modules/local/blast_unmapped_reads'
include { SAMTOOLS_DEPTH                                 						 } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_COVERAGE                              						 } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_STATS                                  						 } from '../../modules/local/samtools_stats'
include { SAMTOOLS_SORT_INDEX                             						 } from '../../modules/local/samtools_sort_index'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED         						 } from '../../modules/local/samtools_fastq'
include { FILTER_ALIGNMENTS_BY_ID                         						 } from '../../modules/local/filter_alignments_by_id'
include { GENERATE_TAXONOMY_TSV                        						 	 } from '../../modules/local/generate_taxonomy_tsv'
include { GENERATE_TAXONOMY_GRAPH as GENERATE_TAXONOMY_GRAPH_READS 				 } from '../../modules/local/generate_taxonomy_graph'

// short read QC modules
include { ALIGN_SHORT_READS_MINIMAP2        } from '../../modules/local/minimap2_short'
include { FASTQC as FASTQC_ALIGNED_FILTERED } from '../../modules/nf-core/fastqc/main'

// long read QC modules
include { ALIGN_LONG_READS_MINIMAP2        				  } from '../../modules/local/minimap2_long'
include { NANOPLOT as NANOPLOT_ALIGNMENT_TAXON_FILTERED   } from '../../modules/local/nanoplot'

workflow ALIGNMENT_CLASSIFICATION_AND_QC {
    take:
    qc_reads // channel: [ val(meta), [ reads ] ]
	tax_map
	ch_multiqc_files

    main:

	if (params.workflow == 'short-read') {

		ALIGN_SHORT_READS_MINIMAP2 (

			qc_reads,
			[[params.minimap2_meta],[params.minimap2_genome_reference]]
		)

		aligner_bam_ch = ALIGN_SHORT_READS_MINIMAP2.out.bam

	} else if (params.workflow == 'long-read') {


		ALIGN_LONG_READS_MINIMAP2 (

			qc_reads,
			[[params.minimap2_meta],[params.minimap2_genome_reference]]
		)

		aligner_bam_ch = ALIGN_LONG_READS_MINIMAP2.out.bam

	}

	IMPORT_SAM_TO_DB (

		aligner_bam_ch,
		tax_map

	)

	if (!params.skip_filter_alignment_by_id) {

		filter_alignment_ch = aligner_bam_ch.join(IMPORT_SAM_TO_DB.out.sam_db)

		FILTER_ALIGNMENTS_BY_ID (

			filter_alignment_ch,
			params.ncbi_taxonomy_nodes,
			params.my_tax_ids,
			params.include_children

		)

		align_class_ch = FILTER_ALIGNMENTS_BY_ID.out.tax_filtered_bam.join(IMPORT_SAM_TO_DB.out.sam_db)
		aligned_bam = FILTER_ALIGNMENTS_BY_ID.out.tax_filtered_bam

	} else {

		align_class_ch = aligner_bam_ch.join(IMPORT_SAM_TO_DB.out.sam_db)
		aligned_bam = aligner_bam_ch
	}


	ALIGNMENT_CLASSIFY (

		align_class_ch

	)

	// ALIGNMENT_CLASSIFICATION_GRAPH_READS (

	// 	align_class_graph_ch,
	// 	"Reads"

	// )

	// ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_CLASSIFICATION_GRAPH_READS.out.plot.collect{it[1]}.ifEmpty([]))

	if (!params.skip_blast_unmapped) {

		BLAST_UNMAPPED_READS (
			ALIGNMENT_CLASSIFY.out.unmapped_bam,
			params.blast_db
		)

	}

	SAMTOOLS_STATS (

		aligned_bam

	)

	SAMTOOLS_SORT_INDEX (

		aligned_bam

	)

	SAMTOOLS_DEPTH (

		aligned_bam

	)

	SAMTOOLS_COVERAGE (

		aligned_bam.join(SAMTOOLS_SORT_INDEX.out.bai)

	)

	// capturing aligned reads and converting to fastq
	SAMTOOLS_FASTQ_MAPPED (

		aligned_bam,
		false

	)


	if (!params.skip_filter_alignment_by_id) {

		tax_report_input_ch = ALIGNMENT_CLASSIFY.out.summary_tsv.join( SAMTOOLS_COVERAGE.out.coverage).join(FILTER_ALIGNMENTS_BY_ID.out.count_tsv)
	
	} else {

		tax_report_input_ch = ALIGNMENT_CLASSIFY.out.summary_tsv.join( SAMTOOLS_COVERAGE.out.coverage)
    		.map { meta, summary_tsv, coverage -> tuple(meta, summary_tsv, coverage, []) }

	}

	GENERATE_TAXONOMY_TSV (

		tax_report_input_ch,
		tax_map

	)

	if (!params.skip_filter_alignment_by_id) {

		tax_graph_input_ch = ALIGNMENT_CLASSIFY.out.summary_tsv.join(GENERATE_TAXONOMY_TSV.out.tax_report_tsv).join(FILTER_ALIGNMENTS_BY_ID.out.count_tsv)
	
	} else {

		tax_graph_input_ch = ALIGNMENT_CLASSIFY.out.summary_tsv.join(GENERATE_TAXONOMY_TSV.out.tax_report_tsv)
    		.map { meta, summary_tsv, tax_report_tsv -> tuple(meta, summary_tsv, tax_report_tsv, []) }

	}

	GENERATE_TAXONOMY_GRAPH_READS (

		tax_graph_input_ch,
		'Reads'

	)

	filtered_reads = SAMTOOLS_FASTQ_MAPPED.out.fastq

	if (params.workflow == 'short-read') {

		FASTQC_ALIGNED_FILTERED (

			filtered_reads

		)

		//ch_multiqc_files = ch_multiqc_files.mix(FASTQC_ALIGNED_FILTERED.out.zip.collect {it[1]}.ifEmpty([]))

	} else if (params.workflow == 'long-read') {

		// running nanoplot again to compare read stats pre and post filter
        NANOPLOT_ALIGNMENT_TAXON_FILTERED (

            filtered_reads

        )

		//ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_ALIGNMENT_TAXON_FILTERED.out.txt.collect{it[1]}.ifEmpty([]))

	}


	emit:

	sorted_bam = SAMTOOLS_SORT_INDEX.out.bam
	sam_db = IMPORT_SAM_TO_DB.out.sam_db
	filtered_reads = filtered_reads
	ch_multiqc_files = ch_multiqc_files

}