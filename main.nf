#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnasplice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnasplice
    Website: https://nf-co.re/rnasplice
    Slack  : https://nfcore.slack.com/channels/rnasplice
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASPLICE               } from './workflows/rnasplice'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Explicitly populate structural genome references from institutional profiles using strict syntax attributes
params.fasta        = getGenomeAttribute('fasta')
params.gff          = getGenomeAttribute('gff')
params.gtf          = getGenomeAttribute('gtf')
params.salmon_index = getGenomeAttribute('salmon')
params.star_index   = getGenomeAttribute('star')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_RNASPLICE {

    main:

    def is_aws_igenome = params.fasta && params.gtf && (file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')

    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME(
        params.fasta,
        params.gtf,
        params.gff,
        params.transcript_fasta,
        params.star_index,
        params.salmon_index,
        params.gff_dexseq,
        params.suppa_tpm,
        params.gencode,
        is_aws_igenome,
    )

    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))
    ch_contrastsheet = channel.value(file(params.contrasts, checkIfExists: true))

    //
    // WORKFLOW: Run pipeline
    //
    RNASPLICE(
        ch_samplesheet,
        ch_contrastsheet,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.dexseq_gff,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.suppa_tpm,
        PREPARE_GENOME.out.chrom_sizes,
        is_aws_igenome,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    emit:
    multiqc_report = RNASPLICE.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNASPLICE()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_RNASPLICE.out.multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}
