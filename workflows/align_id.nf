/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.workflow == 'short-read') {

    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/short-read/multiqc_config.yml", checkIfExists: true)

	ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/short-read/methods_description_template.yml", checkIfExists: true)

} else {

    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/long-read/multiqc_config.yml", checkIfExists: true)

	ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/long-read/methods_description_template.yml", checkIfExists: true)

}

ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK_SHORT 						} from '../subworkflows/local/input_check_short_read'
include { INPUT_CHECK_LONG							} from '../subworkflows/local/input_check_long_read'
include { READ_QC									} from '../subworkflows/local/read_qc'
include { ALIGNMENT_CLASSIFICATION_AND_QC			} from '../subworkflows/local/alignment_classification_and_qc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'


//
// MODULE: custom, local modules
//
include { UPDATE_NODES_DB                                 } from '../modules/local/update_nodes_db'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_SUBSAMPLE      } from '../modules/local/bbmap_reformat_subsample'

//clearing out minimap2 queues if memory_saver mode is enabled (only required for local compute; memory allocation should generally be handled by the job scheduler when submitting to the cluster)
if (params.memory_saver) {

    def minimap2_queue = new File("${projectDir}/queue/minimap2")

    if (minimap2_queue.exists() && minimap2_queue.isDirectory()) {
        minimap2_queue.eachFile { file ->
            file.delete()
        }
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

ch_versions = Channel.empty()
ch_multiqc_files = Channel.empty()

workflow ALIGN_ID {

    if (!params.skip_filter_alignment_by_id && params.include_children) {

        UPDATE_NODES_DB (

            params.local_nodes_db,
            params.ncbi_taxonomy_nodes,
            params.my_tax_ids

        )

        placeholder = UPDATE_NODES_DB.out.complete

    } else {

        placeholder = []

    }

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

	if (params.workflow == 'short-read') {

		INPUT_CHECK_SHORT (

			file(params.input),
			placeholder

		)

		ch_versions = ch_versions.mix(INPUT_CHECK_SHORT.out.versions)

		input_check_ch = INPUT_CHECK_SHORT.out.reads

	} else if (params.workflow == 'long-read') {

		INPUT_CHECK_LONG (

			file(params.input),
			placeholder

		)

		ch_versions = ch_versions.mix(INPUT_CHECK_LONG.out.versions)

		input_check_ch = INPUT_CHECK_LONG.out.reads

	}


    if (!params.skip_subsample) {

        BBMAP_REFORMAT_SUBSAMPLE (
            input_check_ch,
            params.num_subsamples
        )

        raw_reads = BBMAP_REFORMAT_SUBSAMPLE.out.fastq

    } else {

        raw_reads = input_check_ch

    }

    // Subworkflow: QC Reads
	READ_QC (

		raw_reads,
		ch_multiqc_files

	)

	qc_reads = READ_QC.out.qc_reads
	ch_multiqc_files = READ_QC.out.ch_multiqc_files

	// Subworkflow: Align reads to ref genome and assign taxonomy
	ALIGNMENT_CLASSIFICATION_AND_QC (

		qc_reads,
		ch_multiqc_files

	)

	sorted_bam = ALIGNMENT_CLASSIFICATION_AND_QC.out.sorted_bam
	sam_db = ALIGNMENT_CLASSIFICATION_AND_QC.out.sam_db
	filtered_reads = ALIGNMENT_CLASSIFICATION_AND_QC.out.filtered_reads
	ch_multiqc_files = ALIGNMENT_CLASSIFICATION_AND_QC.out.ch_multiqc_files

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowAlignID.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowAlignID.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
