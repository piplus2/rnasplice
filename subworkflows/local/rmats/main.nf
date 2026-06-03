//
// Uncompress and prepare reference genome files
//

include { CREATE_BAMLIST         } from '../../../modules/local/create_bamlist'
include { RMATS_PREP             } from '../../../modules/local/rmats_prep'
include { RMATS_POST             } from '../../../modules/local/rmats_post'

include { CREATE_BAMLIST_SINGLE  } from '../../../modules/local/create_bamlist_single'
include { RMATS_PREP_SINGLE      } from '../../../modules/local/rmats_prep_single'
include { RMATS_POST_SINGLE      } from '../../../modules/local/rmats_post_single'

workflow RMATS {

    take:
    ch_contrastsheet          // channel: contrastsheet
    ch_genome_bam_conditions  // channel: genome_bam_conditions
    gtf                       // channel: /path/to/genome.gtf
    single_condition_bool     // channel: true/false
    rmats_read_len            // params.rmats_read_len
    rmats_splice_diff_cutoff  // params.rmats_splice_diff_cutoff
    rmats_novel_splice_site   // params.rmats_novel_splice_site
    rmats_min_intron_len      // params.rmats_min_intron_len
    rmats_max_exon_len        // params.rmats_max_exon_len
    rmats_paired_stats        // params.rmats_paired_stats

    main:

    //
    // Create contrasts channel for RMATS input
    //

    ch_contrasts = ch_contrastsheet.splitCsv(header: true)

    ch_contrasts_treatment = ch_contrasts
        .map { it -> [it['treatment'], it] }
        .combine(ch_genome_bam_conditions, by: 0)
        .map { _key, contrast, metas, bams -> contrast + ['meta1': metas, 'bam1': bams] }

    ch_contrasts_full = ch_contrasts_treatment
        .map { it -> [it['control'], it] }
        .combine(ch_genome_bam_conditions, by: 0)
        .map { _key, contrast, metas, bams -> contrast + ['meta2': metas, 'bam2': bams] }
        .tap { ch_contrasts_full_for_bamlist }

    //
    // Run RMATS in single or contrast mode
    //

    if (single_condition_bool) {

        //
        // Create input channels
        //

        ch_bam = ch_contrasts_full.map { it -> [ it.contrast, it.treatment, it.bam1 ] }

        //
        // Create input bam list file
        //

        CREATE_BAMLIST_SINGLE (
            ch_bam
        )

        //
        // Join bamlist with contrasts
        //

        ch_contrasts_full = ch_contrasts_full_for_bamlist
            .map { it -> [it['contrast'], it] }
            .combine ( CREATE_BAMLIST_SINGLE.out.bamlist, by: 0 )
            .map { it -> it[1] + ['bam1_text': it[2]] }

        //
        // Create input channels
        //

        ch_contrasts_bamlist = ch_contrasts_full.map { it -> [ it.contrast, it.treatment, it.meta1, it.bam1, it.bam1_text ] }

        //
        // Run rMATS prep step
        //

        RMATS_PREP_SINGLE (
            gtf,
            ch_contrasts_bamlist,
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len
        )

        ch_rmats_temp = RMATS_PREP_SINGLE.out.rmats_temp

        //
        // Join rmats temp with contrasts
        //

        ch_contrasts_full = ch_contrasts_full
            .map { it -> [it['contrast'], it] }
            .combine ( ch_rmats_temp, by: 0 )
            .map { it -> it[1] + ['rmats_temp': it[2]] }

        //
        // Create input channels
        //

        ch_contrasts_bamlist = ch_contrasts_full.map { it -> [ it.contrast, it.treatment, it.meta1, it.bam1, it.bam1_text, it.rmats_temp ] }

        //
        // Run rMATs post step
        //

        RMATS_POST_SINGLE (
            gtf,
            ch_contrasts_bamlist,
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len,
            rmats_paired_stats
        )

        ch_rmats_prep       = RMATS_PREP_SINGLE.out.rmats_temp  //    path: rmats_temp/*
        ch_rmats_prep_log   = RMATS_PREP_SINGLE.out.log         //    path: rmats_prep.log
        ch_rmats_post       = RMATS_POST_SINGLE.out.rmats_post  //    path: rmats_post/*
        ch_rmats_post_log   = RMATS_POST_SINGLE.out.log         //    path: rmats_post.log

    } else {

        //
        // Create input channel
        //

        ch_bam = ch_contrasts_full.map { it -> [ it.contrast, it.treatment, it.control, it.bam1, it.bam2 ] }

        //
        // Create input bam list file
        //

        CREATE_BAMLIST (
            ch_bam
        )

        //
        // Join bamlist with contrasts
        //

        ch_contrasts_full = ch_contrasts_full_for_bamlist
            .map { it -> [it['contrast'], it] }
            .combine ( CREATE_BAMLIST.out.bamlist, by: 0 )
            .map { _key, contrast, bam1_text, bam2_text -> contrast + ['bam1_text': bam1_text, 'bam2_text': bam2_text] }

        //
        // Create input channels
        //

        ch_contrasts_bamlist = ch_contrasts_full.map { it -> [ it.contrast, it.treatment, it.control, it.meta1, it.meta2, it.bam1, it.bam2, it.bam1_text, it.bam2_text ] }

        //
        // Run rMATS prep step
        //

        RMATS_PREP (
            gtf,
            ch_contrasts_bamlist,
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len
        )

        ch_rmats_temp = RMATS_PREP.out.rmats_temp

        //
        // Join rmats temp with contrasts
        //

        ch_contrasts_full = ch_contrasts_full
            .map { it -> [it['contrast'], it] }
            .combine ( ch_rmats_temp, by: 0 )
            .map { it -> it[1] + ['rmats_temp': it[2]] }

        //
        // Create input channels
        //

        ch_contrasts_bamlist = ch_contrasts_full.map { it -> [ it.contrast, it.treatment, it.control, it.meta1, it.meta2, it.bam1, it.bam2, it.bam1_text, it.bam2_text, it.rmats_temp ] }

        //
        // Run rMATs post step
        //

        RMATS_POST (
            gtf,
            ch_contrasts_bamlist,
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len,
            rmats_paired_stats
        )

        ch_rmats_prep         = RMATS_PREP.out.rmats_temp     //    path: rmats_prep/*
        ch_rmats_prep_log     = RMATS_PREP.out.log            //    path: rmats_prep.log
        ch_rmats_post         = RMATS_POST.out.rmats_post     //    path: rmats_post/*
        ch_rmats_post_log     = RMATS_POST.out.log            //    path: rmats_post.log

    }

    // Define output

    emit:

    rmats_prep         = ch_rmats_prep        //    path: rmats_temp/*
    rmats_prep_log     = ch_rmats_prep_log    //    path: rmats_prep.log
    rmats_post         = ch_rmats_post        //    path: rmats_post/*
    rmats_post_log     = ch_rmats_post_log    //    path: rmats_post.log
}
