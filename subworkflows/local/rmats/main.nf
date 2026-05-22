//
// Uncompress and prepare reference genome files
//

include { CREATE_BAMLIST        } from '../../../modules/local/create_bamlist'
include { RMATS_PREP            } from '../../../modules/local/rmats_prep'
include { RMATS_POST            } from '../../../modules/local/rmats_post'

include { CREATE_BAMLIST_SINGLE } from '../../../modules/local/create_bamlist_single'
include { RMATS_POST_SINGLE     } from '../../../modules/local/rmats_post_single'

workflow RMATS {
    take:
    ch_contrastsheet // channel: contrastsheet
    ch_genome_bam_conditions // channel: genome_bam_conditions
    gtf // channel: /path/to/genome.gtf
    single_condition_bool // channel: true/false
    rmats_read_len // params.rmats_read_len
    rmats_splice_diff_cutoff // params.rmats_splice_diff_cutoff
    rmats_novel_splice_site // params.rmats_novel_splice_site
    rmats_min_intron_len // params.rmats_min_intron_len
    rmats_max_exon_len // params.rmats_max_exon_len
    rmats_paired_stats // params.rmats_paired_stats

    main:

    ch_versions = channel.empty()

    //
    // Create contrasts channel for RMATS input
    //

    ch_contrasts = ch_contrastsheet.splitCsv(header: true)

    ch_contrasts_treatment = ch_contrasts
        .map { it -> [it['treatment'], it] }
        .combine(ch_genome_bam_conditions, by: 0)
        .map { it -> it[1] + ['meta1': it[2], 'bam1': it[3]] }

    ch_contrasts_full = ch_contrasts_treatment
        .map { it -> [it['control'], it] }
        .combine(ch_genome_bam_conditions, by: 0)
        .map { it -> it[1] + ['meta2': it[2], 'bam2': it[3]] }

    //
    // Create input channels
    //
    if (single_condition_bool) {
        ch_contrasts_bam = ch_contrasts_full.map { it -> [it.contrast, it.treatment, it.meta1, it.bam1, null] }
    }
    else {
        ch_contrasts_bam = ch_contrasts_full.map { it -> [it.contrast, it.treatment, it.control, it.meta1, it.meta2, it.bam1, it.bam2] }
    }

    //
    // Run RMATS in single or contrast mode
    //

    CREATE_BAMLIST(
        ch_contrasts_bam
    )

    if (single_condition_bool) {
        ch_contrasts_with_bamlist = ch_contrasts_full
            .map { it -> [it['contrast'], it] }
            .combine(CREATE_BAMLIST.out.bamlist1, by: 0)
            .map { it -> it[1] + ['bam1_text': it[2]] }

        ch_prep_input = ch_contrasts_with_bamlist.map { it -> [it.contrast, it.treatment, it.meta1, it.bam1, it.bam1_text, null, null, null, null] }
    }
    else {
        ch_contrasts_with_bamlist = ch_contrasts_full
            .map { it -> [it['contrast'], it] }
            .combine(CREATE_BAMLIST.out.bamlist1, by: 0)
            .combine(CREATE_BAMLIST.out.bamlist2, by: 0) { c1, c2 -> c1 + [bam2_text: c2[2]] }
            .map { it -> it[1] + ['bam1_text': it[2]] }
        ch_prep_input = ch_contrasts_with_bamlist.map { it -> [it.contrast, it.treatment, it.meta1, it.bam1, it.bam1_text, it.control, it.meta2, it.bam2, it.bam2_text] }
    }

    //
    // Run rMATS prep step
    //
    RMATS_PREP(
        gtf,
        ch_prep_input,
        rmats_read_len,
        rmats_splice_diff_cutoff,
        rmats_novel_splice_site,
        rmats_min_intron_len,
        rmats_max_exon_len,
    )

    //
    // Prepare input for rMATS post step
    //
    if (single_condition_bool) {
        ch_post_input = ch_contrasts_with_bamlist
            .map { it -> [it['contrast'], it] }
            .join(RMATS_PREP.out.rmats_temp, by: 0)
            .map { contrasts, meta, temp ->
                [contrasts, meta.treatment, meta.meta1, meta.bam1, meta.bam1_text, null, null, null, null, temp]
            }
    }
    else {
        ch_contrasts = ch_contrasts
            .map { it -> [it['contrast'], it] }
            .combine(RMATS_PREP.out.rmats_temp, by: 0)
            .map { it -> it[1] + ['rmats_temp': it[2]] }
    }

    RMATS_POST(
        gtf,
        ch_post_input,
        rmats_read_len,
        rmats_splice_diff_cutoff,
        rmats_novel_splice_site,
        rmats_min_intron_len,
        rmats_max_exon_len,
        rmats_paired_stats,
    )

    ch_versions = ch_versions.mix(RMATS_POST.out.versions)

    ch_rmats_prep = RMATS_PREP.out.rmats_temp
    ch_rmats_prep_log = RMATS_PREP.out.log
    ch_rmats_post = RMATS_POST.out.rmats_post
    ch_rmats_post_log = RMATS_POST.out.log

    emit:
    rmats_prep     = ch_rmats_prep
    rmats_prep_log = ch_rmats_prep_log
    rmats_post     = ch_rmats_post
    rmats_post_log = ch_rmats_post_log
    versions       = ch_versions.ifEmpty(null)
}
