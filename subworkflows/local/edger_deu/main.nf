//
// edgeR DEU subworkflow
//

include { SUBREAD_FLATTENGTF    } from '../../../modules/local/subread/flattengtf'
include { SUBREAD_FEATURECOUNTS } from '../../../modules/nf-core/subread/featurecounts'
include { EDGER_EXON            } from '../../../modules/local/edger_exon'

workflow EDGER_DEU {
    take:
    gtf // path: gtf
    ch_genome_bam // channel: [ val(meta), path(bams) ]
    ch_samplesheet // channel.fromPath(params.input)
    ch_contrastsheet // channel.fromPath(params.contrasts)
    n_edger_plot // val: integer to plot

    main:

    // MODULE: SUBREAD_FLATTENGTF

    SUBREAD_FLATTENGTF(gtf.map { file -> [ [id: file.baseName], file ] })

    //
    // MODULE: SUBREAD_FEATURECOUNTS
    //

    ch_feature_counts = ch_genome_bam.combine(SUBREAD_FLATTENGTF.out.saf.map { _meta, saf -> saf })

    SUBREAD_FEATURECOUNTS(ch_feature_counts)

    //
    // MODULE: EDGER_COUNTS AND PLOT
    //
    EDGER_EXON(
        SUBREAD_FEATURECOUNTS.out.counts.collect { it -> it[1] },
        ch_samplesheet,
        ch_contrastsheet,
        n_edger_plot,
    )

    emit:
    featureCounts_summary = SUBREAD_FEATURECOUNTS.out.summary // path featureCounts.txt.summary
}
