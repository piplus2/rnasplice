//
// Alignment with STAR
//

include { STAR_ALIGN              } from '../../../modules/nf-core/star/align/main'
include { STAR_ALIGN_IGENOMES     } from '../../../modules/local/star_align_igenomes'
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'

include { isAwsIgenome            } from '../utils_nfcore_rnasplice_pipeline'

workflow ALIGN_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/star/index/
    gtf // channel: /path/to/genome.gtf
    fasta // channel: /path/to/fasta

    main:

    // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    def star_ignore_sjdbgtf = params.star_ignore_sjdbgtf
    def seq_platform = '' // string : sequencing platform
    def seq_center = params.seq_center ?: '' // string : sequencing center
    def is_aws_igenome = isAwsIgenome() // boolean: whether the genome files are from AWS iGenomes

    //
    // Map reads with STAR
    //

    ch_orig_bam = channel.empty()
    ch_log_final = channel.empty()
    ch_log_out = channel.empty()
    ch_log_progress = channel.empty()
    ch_bam_sorted = channel.empty()
    ch_bam_transcript = channel.empty()
    ch_fastq = channel.empty()
    ch_tab = channel.empty()

    if (is_aws_igenome) {
        STAR_ALIGN_IGENOMES(reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        ch_orig_bam = STAR_ALIGN_IGENOMES.out.bam
        ch_log_final = STAR_ALIGN_IGENOMES.out.log_final
        ch_log_out = STAR_ALIGN_IGENOMES.out.log_out
        ch_log_progress = STAR_ALIGN_IGENOMES.out.log_progress
        ch_bam_sorted = STAR_ALIGN_IGENOMES.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN_IGENOMES.out.bam_transcript
        ch_fastq = STAR_ALIGN_IGENOMES.out.fastq
        ch_tab = STAR_ALIGN_IGENOMES.out.tab
    }
    else {
        STAR_ALIGN(reads, index, gtf, star_ignore_sjdbgtf)
        ch_orig_bam = STAR_ALIGN.out.bam
        ch_log_final = STAR_ALIGN.out.log_final
        ch_log_out = STAR_ALIGN.out.log_out
        ch_log_progress = STAR_ALIGN.out.log_progress
        ch_bam_sorted = STAR_ALIGN.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN.out.bam_transcript
        ch_fastq = STAR_ALIGN.out.fastq
        ch_tab = STAR_ALIGN.out.tab
    }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //

    BAM_SORT_STATS_SAMTOOLS(
        ch_orig_bam,
        fasta,
    )

    emit:
    orig_bam       = ch_orig_bam // channel: [ val(meta), bam            ]
    log_final      = ch_log_final // channel: [ val(meta), log_final      ]
    log_out        = ch_log_out // channel: [ val(meta), log_out        ]
    log_progress   = ch_log_progress // channel: [ val(meta), log_progress   ]
    bam_sorted     = ch_bam_sorted // channel: [ val(meta), bam_sorted     ]
    bam_transcript = ch_bam_transcript // channel: [ val(meta), bam_transcript ]
    fastq          = ch_fastq // channel: [ val(meta), fastq          ]
    tab            = ch_tab // channel: [ val(meta), tab            ]
    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam // channel: [ val(meta), [ bam ] ]
    index          = BAM_SORT_STATS_SAMTOOLS.out.index // channel: [ val(meta), [ bai/csi ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
