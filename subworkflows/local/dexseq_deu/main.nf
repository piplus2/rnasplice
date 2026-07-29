//
// DEXSeq DEU subworkflow
//

include { DEXSEQ_ANNOTATION   } from '../../../modules/local/dexseq/annotation'
include { DEXSEQ_COUNT        } from '../../../modules/local/dexseq/count'
include { DEXSEQ_EXON         } from '../../../modules/local/dexseq/exon'

workflow DEXSEQ_DEU {

    take:

    gtf                // path gtf
    ch_genome_bam      // bam channel
    ch_dexseq_gff      // path dexseq gff
    ch_samplesheet     // channel: path(samplesheet)
    ch_contrastsheet   // channel: path(contrastsheet)

    main:

    def n_dexseq_plot     = params.n_dexseq_plot     // val: number of genes to plot
    def aggregation       = params.aggregation       // val: aggregate overlapping genes
    def alignment_quality = params.alignment_quality // val: minimum alignment quality

    if (!params.gff_dexseq) {

        //
        // MODULE: DEXSeq Annotation
        //

        DEXSEQ_ANNOTATION (
            gtf,
            aggregation
        )

        ch_dexseq_gff = DEXSEQ_ANNOTATION.out.gff

    }

    ch_genome_bam_gff = ch_genome_bam.combine(ch_dexseq_gff)

    //
    // MODULE: DEXSeq Count
    //

    DEXSEQ_COUNT (
        ch_genome_bam_gff,
        alignment_quality
    )

    //
    // MODULE: DEXSeq DEU
    //

    DEXSEQ_EXON (
        DEXSEQ_COUNT.out.dexseq_clean_txt.map{ it[1] }.collect(),
        ch_dexseq_gff,
        ch_samplesheet,
        ch_contrastsheet,
        n_dexseq_plot
    )

    emit:

    dexseq_clean_txt        = DEXSEQ_COUNT.out.dexseq_clean_txt.map{ it[1] }.collect()

    dexseq_exon_dataset_rds = DEXSEQ_EXON.out.dexseq_exon_dataset_rds
    dexseq_exon_results_rds = DEXSEQ_EXON.out.dexseq_exon_results_rds
    dexseq_gene_results_rds = DEXSEQ_EXON.out.dexseq_gene_results_rds
    dexseq_exon_results_csv = DEXSEQ_EXON.out.dexseq_exon_results_csv
    dexseq_gene_results_csv = DEXSEQ_EXON.out.dexseq_gene_results_csv
    dexseq_plot_results_pdf = DEXSEQ_EXON.out.dexseq_plot_results_pdf
}
