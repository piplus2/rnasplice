//
// Uncompress and prepare reference genome files
//

include { CREATE_BAMLIST } from '../../../modules/local/create_bamlist'
include { RMATS_PREP     } from '../../../modules/local/rmats_prep'
include { RMATS_POST     } from '../../../modules/local/rmats_post'

workflow RMATS {
    take:
    ch_contrastsheet // channel: contrastsheet
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

    ch_contrasts = ch_contrastsheet.splitCsv(header: true)

    if (is_single_condition) {

        //
        // SINGLE CONDITION MODE: samples have all the same condition
        //

        ch_contrasts_bamlist = ch_genome_bam_conditions
            .groupTuple(by: 0) // group by condition
            .map { condition, metas, bams ->
                def contrast_name = "${condition}_profiling"
                return [ contrast_name, condition, metas, bams ]
            }

        CREATE_BAMLIST(
            ch_contrasts_bamlist
                .map {
                    contrast, cond1, _meta1, bam1 -> [contrast, cond1, bam1, '', [] ]
                }
        )

        ch_prep_ready = ch_contrasts_bamlist
            .join( CREATE_BAMLIST.out.bamlists, by: 0 )
            .map { contrast, remainder, list_txts ->
                def (cond1, meta1, bam1) = remainder
                def bam1_txt = list_txts[0]
                return [ contrast, cond1, meta1, bam1, bam1_txt ]
            }

        RMATS_PREP(
            gtf,
            ch_prep_ready
                .map { contrast, cond1, meta1, bam1, bam1_txt -> [ contrast, cond1, meta1, bam1, bam1_txt, '', [], [], '' ] },
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len,
        )

        ch_post_ready = ch_prep_ready
            .join(RMATS_PREP.out.rmats_temp, by: 0)
            .map { contrast, remainder, rmats_temp ->
                def (cond1, meta1, bam1, bam1_txt) = remainder
                return [ contrast, cond1, meta1, bam1, bam1_txt, rmats_temp ]
            }

        RMATS_POST(
            gtf,
            ch_post_ready
                .map { contrast, cond1, meta1, bam1, bam1_txt, rmats_temp -> [ contrast, cond1, meta1, bam1, bam1_txt, '', [], [], '', rmats_temp ] },
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len,
            rmats_paired_stats,
        )

    } else {

        //
        // TWO CONDITIONS MODE
        //

        if (rmats_paired_stats) {

            //
            // PAIRED SAMPLES
            //
            ch_indexed_bams = ch_genome_bam_conditions
                .groupTuple(by: 0) // group by condition
                .flatMap { condition, metas, bams ->
                    metas.withIndex().collect {
                        meta, idx -> [ condition, idx, meta, bams[idx] ]
                    }
                }

            ch_tx = ch_contrasts
                .map { row -> [ row.treatment, row ] }
                .combine( ch_indexed_bams, by: 0 )
                .map { _condition, row, idx, meta, bam ->
                    [ "${row.contrast}_pair${idx}", row + [tx_meta: meta, tx_bam: bam] ]
                }

            ch_ctrl = ch_contrasts
                .map { row -> [ row.control, row ] }
                .combine( ch_indexed_bams, by: 0 )
                .map { _condition, row, idx, meta, bam ->
                    [ "${row.contrast}_pair${idx}", row + [ctrl_meta: meta, ctrl_bam: bam] ]
                }

            ch_fully_paired = ch_tx
                .join ( ch_ctrl, by: 0 )
                .map { _pair_key, tx, ctrl -> tx + ctrl }

            ch_contrasts_bamlist = ch_fully_paired
                .map { it -> [ it.contrast, it ] }
                .groupTuple(by: 0)
                .map { contrast, pairs ->
                    def cond1 = pairs[0].treatment
                    def cond2 = pairs[0].control
                    def meta1 = pairs.collect { it -> it.tx_meta }
                    def bam1  = pairs.collect { it -> it.tx_bam }
                    def meta2 = pairs.collect { it -> it.ctrl_meta }
                    def bam2  = pairs.collect { it -> it.ctrl_bam }
                    return [ contrast, cond1, meta1, bam1, cond2, meta2, bam2 ]
                }

        } else {

            //
            // UNPAIRED SAMPLES
            //

            ch_grouped_bams = ch_genome_bam_conditions
                .groupTuple(by: 0) // Output: [ condition, [metas], [bams] ]

            ch_contrasts_bamlist = ch_contrasts
                .map { row -> [ row.treatment, row ] }
                .join( ch_grouped_bams, by: 0 )
                .map { _tx, row, tx_metas, tx_bams ->
                    [ row.control, row + [ tx_metas: tx_metas, tx_bams: tx_bams ] ]
                }
                .join( ch_grouped_bams, by: 0 )
                .map { _ctrl, row, ctrl_metas, ctrl_bams ->
                    [ row.contrast, row.treatment, row.tx_metas, row.tx_bams, row.ctrl, ctrl_metas, ctrl_bams ]
                }
        }

        CREATE_BAMLIST(
            ch_contrasts_bamlist
                .map {
                    contrast, cond1, _meta1, bam1, cond2, _meta2, bam2
                    -> [ contrast, cond1, bam1, cond2, bam2 ]
                }
        )

        ch_prep_ready = ch_contrasts_bamlist
            .join( CREATE_BAMLIST.out.bamlists, by: 0 )
            .map { contrast, remainder, list_txts ->
                def (cond1, meta1, bam1, cond2, meta2, bam2) = remainder
                def bam1_txt = list_txts[0]
                def bam2_txt = list_txts[1]
                return [ contrast, cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt ]
            }

        RMATS_PREP(
            gtf,
            ch_prep_ready,
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len,
        )

        ch_post_ready = ch_prep_ready
            .join( RMATS_PREP.out.rmats_temp, by: 0 )
            .map { contrast, remainder, rmats_temp ->
                def (cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt) = remainder
                return [ contrast, cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt, rmats_temp ]
            }

        RMATS_POST(
            gtf,
            ch_post_ready,
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len,
            rmats_paired_stats,
        )

    }

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
