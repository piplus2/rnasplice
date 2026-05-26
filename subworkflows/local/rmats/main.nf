//
// Uncompress and prepare reference genome files
//

include { CREATE_BAMLIST } from '../../../modules/local/create_bamlist'
include { RMATS_PREP     } from '../../../modules/local/rmats_prep'
include { RMATS_POST     } from '../../../modules/local/rmats_post'

workflow RMATS {
    take:
    ch_contrasts // channel: contrasts = [ contrast, treatment, control ]
    ch_genome_bam_conditions // channel: genome_bam_conditions
    gtf // channel: /path/to/genome.gtf
    is_single_condition // channel: true/false
    rmats_read_len // params.rmats_read_len
    rmats_splice_diff_cutoff // params.rmats_splice_diff_cutoff
    rmats_novel_splice_site // params.rmats_novel_splice_site
    rmats_min_intron_len // params.rmats_min_intron_len
    rmats_max_exon_len // params.rmats_max_exon_len
    rmats_paired_stats // params.rmats_paired_stats

    main:

    //
    // Create contrasts channel for RMATS input
    //

    ch_contrasts_treatment = ch_contrasts
        .map { it -> [it['treatment'], it] }
        .combine(ch_genome_bam_conditions, by: 0)
        .map { _key, contrast, metas, bams -> contrast + ['meta1': metas, 'bam1': bams] }

    ch_contrasts_full = ch_contrasts_treatment
        .map { it -> [it['control'], it] }
        .combine(ch_genome_bam_conditions, by: 0)
        .map { _key, contrast, metas, bams -> contrast + ['meta2': metas, 'bam2': bams] }

    //
    // Create input channels
    //
    if (is_single_condition) {
        ch_contrasts_bam = ch_contrasts_full.map { it -> [it.contrast, it.treatment, it.treatment, it.bam1, it.bam1] }
    }
    else {
        ch_contrasts_bam = ch_contrasts_full.map { it -> [it.contrast, it.treatment, it.control, it.bam1, it.bam2] }
    }

    //
    // Run RMATS in single or contrast mode
    //

    CREATE_BAMLIST(
        ch_contrasts_bam
    )

    if (is_single_condition) {
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
            .map { key, contrast, bamlist1 ->
                [key, contrast + ['bam1_text': bamlist1]]
            }
            .combine(CREATE_BAMLIST.out.bamlist2, by: 0)
            .map { _key, contrast, bamlist2 ->
                contrast + ['bam2_text': bamlist2]
            }

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
    if (is_single_condition) {
        ch_post_input = ch_contrasts_with_bamlist
            .map { it -> [it['contrast'], it] }
            .join(RMATS_PREP.out.rmats_temp, by: 0)
            .map { contrast, meta, temp ->
                [contrast, meta.treatment, meta.meta1, meta.bam1, meta.bam1_text, null, null, null, null, temp]
            }
    }
    else {
        ch_post_input = ch_contrasts_with_bamlist
            .map { it -> [it['contrast'], it] }
            .join(RMATS_PREP.out.rmats_temp, by: 0)
            .map { contrast, meta, temp ->
                [contrast, meta.treatment, meta.meta1, meta.bam1, meta.bam1_text, meta.control, meta.meta2, meta.bam2, meta.bam2_text, temp]
            }
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

    ch_rmats_prep = RMATS_PREP.out.rmats_temp
    ch_rmats_prep_log = RMATS_PREP.out.log
    ch_rmats_post = RMATS_POST.out.rmats_post
    ch_rmats_post_log = RMATS_POST.out.log

    emit:
    rmats_prep     = ch_rmats_prep
    rmats_prep_log = ch_rmats_prep_log
    rmats_post     = ch_rmats_post
    rmats_post_log = ch_rmats_post_log
}
