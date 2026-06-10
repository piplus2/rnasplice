//
// rMATS differential splicing analysis
//

include { CREATE_BAMLIST         } from '../../../modules/local/create_bamlist'
include { RMATS_PREP             } from '../../../modules/local/rmats_prep'
include { RMATS_POST             } from '../../../modules/local/rmats_post'


workflow RMATS {
    take:
    ch_samplesheet           // channel: path(samplesheet.csv)
    ch_contrastsheet         // channel: path(contrastsheet.csv)
    ch_genome_bam_conditions // channel: [ condition, meta, bam ]
    gtf                      // channel: path(genome.gtf)
    is_single_condition      // val: true/false
    rmats_read_len
    rmats_splice_diff_cutoff
    rmats_novel_splice_site
    rmats_min_intron_len
    rmats_max_exon_len
    rmats_paired_stats

    main:

    if (is_single_condition) {

        //
        // SINGLE CONDITION MODE
        //

        ch_contrasts_bamlist = ch_genome_bam_conditions
            .groupTuple(by: 0)
            .map { condition, metas, bams ->
                def contrast_name = "${condition}_profiling"
                return [ contrast_name, condition, metas, bams ]
            }

        CREATE_BAMLIST(
            ch_contrasts_bamlist
                .map { contrast, cond1, _meta1, bam1 ->
                    [ contrast, cond1, bam1, '', [] ]
                }
        )

        ch_prep_ready = ch_contrasts_bamlist
            .join( CREATE_BAMLIST.out.bam_list1, by: 0 )
            .map { contrast, cond1, meta1, bam1, bam1_txt ->
                return [ contrast, cond1, meta1, bam1, bam1_txt ]
            }

        RMATS_PREP(
            gtf,
            ch_prep_ready
                .map { contrast, cond1, meta1, bam1, bam1_txt ->
                    [ contrast, cond1, meta1, bam1, bam1_txt, '', [], [], '' ]
                },
            rmats_read_len,
            rmats_splice_diff_cutoff,
            rmats_novel_splice_site,
            rmats_min_intron_len,
            rmats_max_exon_len,
        )

        ch_post_ready = ch_prep_ready
            .join( RMATS_PREP.out.rmats_temp, by: 0 )
            .map { contrast, cond1, meta1, bam1, bam1_txt, rmats_temp ->
                [ contrast, cond1, meta1, bam1, bam1_txt, '', [], [], '', rmats_temp ]
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

    } else {

        //
        // TWO CONDITIONS MODE
        //

        if (rmats_paired_stats) {

            //
            // PAIRED SAMPLES
            //

            // Re-parse samplesheet to get stable per-sample index
            // Deduplicating first handles multi-lane samples
            ch_sample_index = ch_samplesheet
                .splitCsv(header: true)
                .collect()
                .flatMap { rows ->
                    def by_condition = [:]
                    rows.each { row ->
                        by_condition.computeIfAbsent(row.condition) { [] } << row.sample
                    }
                    by_condition.collectMany { _condition, samples ->
                        samples.withIndex().collect { sample, idx -> [ sample, idx ] }
                    }
                }

            ch_genome_bam_conditions
                .multiMap { condition, meta, bam ->
                    tx:   [ condition, meta, bam ]
                    ctrl: [ condition, meta, bam ]
                }
                .set { ch_bams_fork }

            ch_indexed_bams_tx = ch_bams_fork.tx
                .map { condition, meta, bam -> [ meta.id, condition, meta, bam ] }
                .join( ch_sample_index, by: 0 )
                .map { _id, condition, meta, bam, idx -> [ condition, idx, meta, bam ] }

            ch_indexed_bams_ctrl = ch_bams_fork.ctrl
                .map { condition, meta, bam -> [ meta.id, condition, meta, bam ] }
                .join( ch_sample_index, by: 0 )
                .map { _id, condition, meta, bam, idx -> [ condition, idx, meta, bam ] }


            // Build per-contrast pairs keyed by contrast__sample_id
            ch_tx = ch_contrastsheet.splitCsv(header: true)
                .map { row -> [ row.treatment, row ] }           // key: treatment condition
                .combine( ch_indexed_bams_tx, by: 0 )            // match BAMs by condition
                .map { _condition, row, idx, meta, bam ->
                    [ row.contrast, meta.id, idx, row, meta, bam ]  // key: contrast + sample id (separate)
                }

            ch_ctrl = ch_contrastsheet.splitCsv(header: true)
                .map { row -> [ row.control, row ] }             // key: control condition
                .combine( ch_indexed_bams_ctrl, by: 0 )          // match BAMs by condition
                .map { _condition, row, idx, meta, bam ->
                    [ row.contrast, meta.id, idx, row, meta, bam ]  // same contrast + sample id
                }

            // Join treatment and control on contrast__sample_id
            ch_fully_paired = ch_tx
                .join( ch_ctrl, by: [0, 2] )  // join on both contrast AND sample id
                .map { _contrast, _sample_id, tx_idx, tx_row, tx_meta, tx_bam, ctrl_idx, _ctrl_row, ctrl_meta, ctrl_bam ->
                    tx_row + [
                        tx_idx   : tx_idx,
                        tx_meta  : tx_meta,
                        tx_bam   : tx_bam,
                        ctrl_idx : ctrl_idx,
                        ctrl_meta: ctrl_meta,
                        ctrl_bam : ctrl_bam
                    ]
                }

            ch_contrasts_bamlist = ch_fully_paired
                .map { it -> [ it.contrast, it ] }
                .groupTuple(by: 0)
                .map { contrast, pairs ->
                    def cond1  = pairs[0].treatment
                    def cond2  = pairs[0].control
                    def sorted = pairs.sort { it -> it.tx_idx }
                    def meta1  = sorted.collect { it -> it.tx_meta }
                    def bam1   = sorted.collect { it -> it.tx_bam }
                    def meta2  = sorted.collect { it -> it.ctrl_meta }
                    def bam2   = sorted.collect { it -> it.ctrl_bam }

                    if (meta1.size() != meta2.size()) {
                        error(
                            "Paired rMATS contrast '${contrast}': unequal sample counts — " +
                            "${cond1}: ${meta1.collect { it -> it.id }}, " +
                            "${cond2}: ${meta2.collect { it -> it.id }}. " +
                            "Each sample ID must appear in both conditions."
                        )
                    }

                    return [ contrast, cond1, meta1, bam1, cond2, meta2, bam2 ]
                }

        } else {

            //
            // UNPAIRED SAMPLES
            //

            ch_contrasts = ch_contrastsheet.splitCsv(header: true)

            ch_grouped_bams = ch_genome_bam_conditions
                .groupTuple(by: 0)

            ch_contrasts_bamlist = ch_contrasts
                .map { row -> [ row.treatment, row ] }
                .join( ch_grouped_bams, by: 0 )
                .map { _tx, row, tx_metas, tx_bams ->
                    [ row.control, row + [ tx_metas: tx_metas, tx_bams: tx_bams ] ]
                }
                .join( ch_grouped_bams, by: 0 )
                .map { _ctrl, row, ctrl_metas, ctrl_bams ->
                    [ row.contrast, row.treatment, row.tx_metas, row.tx_bams, row.control, ctrl_metas, ctrl_bams ]
                }
        }

        CREATE_BAMLIST(
            ch_contrasts_bamlist
                .map { contrast, cond1, _meta1, bam1, cond2, _meta2, bam2 ->
                    [ contrast, cond1, bam1, cond2, bam2 ]
                }
        )

        ch_prep_ready = ch_contrasts_bamlist
            .join( CREATE_BAMLIST.out.bam_list1, by: 0 )
            .join( CREATE_BAMLIST.out.bam_list2, by: 0 )
            .map { contrast, cond1, meta1, bam1, cond2, meta2, bam2, bam1_txt, bam2_txt ->
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
            .map { contrast, cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt, rmats_temp ->
                [ contrast, cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt, rmats_temp ]
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

    ch_rmats_prep     = RMATS_PREP.out.rmats_temp
    ch_rmats_prep_log = RMATS_PREP.out.log
    ch_rmats_post     = RMATS_POST.out.rmats_post
    ch_rmats_post_log = RMATS_POST.out.log

    emit:
    rmats_prep     = ch_rmats_prep
    rmats_prep_log = ch_rmats_prep_log
    rmats_post     = ch_rmats_post
    rmats_post_log = ch_rmats_post_log
}
