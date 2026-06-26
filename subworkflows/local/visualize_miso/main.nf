//
// Visualise miso subworkflow
//

include { GTF_2_GFF3    } from '../../../modules/local/gtf_2_gff3'
include { MISOPY_INDEX  } from '../../../modules/local/misopy/index'
include { MISOPY_RUN      } from '../../../modules/local/misopy/run'
include { MISO_SETTINGS } from '../../../modules/local/miso_settings'
include { MISO_SASHIMI  } from '../../../modules/local/miso_sashimi'


workflow VISUALISE_MISO {

    take:

    gtf                // path gtf
    ch_genome_bam      // channel: [ val(meta), path(bams) ]
    ch_genome_bai      // channel: [ val(meta), path(bais) ]
    fig_width          // 7
    fig_height         // 5
    miso_genes         // params.miso_genes
    miso_genes_file    // params.miso_genes_file

    main:

   //
   // MODULE: gtf_2_gff3
   //

    GTF_2_GFF3 (gtf)

   //
   // MODULE: DEXSeq Annotation
   //

    MISOPY_INDEX(GTF_2_GFF3.out.gff3.map{ gff -> [ [:], gff ] })

    //
    // MODULE: MISOPY_RUN
    //

    ch_bam_join = ch_genome_bam.join(ch_genome_bai)
    ch_miso_index = MISOPY_INDEX.out.miso_index.collect().map{ it -> it[1] }

    MISOPY_RUN (
        ch_bam_join,
        ch_miso_index
    )

    //
    // MODULE: MISO_SETTINGS
    //

    ch_bams = ch_genome_bam.collect({it -> it[1]})
    ch_miso_run = MISOPY_RUN.out.miso.map{it -> it[1]}.collect()

    MISO_SETTINGS (
        ch_miso_run,
        ch_bams,
        fig_width,
        fig_height
    )

    //
    // MODULE: MISO_SASHIMI
    //

    def miso_genes_list = miso_genes ? miso_genes.split(',').collect{ it -> it.trim() } : [""]
    ch_miso_genes_list = channel.fromList( miso_genes_list )

    if (miso_genes_file && miso_genes) {
        ch_miso_genes_file = channel.fromPath(miso_genes_file)
            .splitCsv()
        ch_miso_genes_list.concat( ch_miso_genes_file )
        .set{ ch_miso_genes }
    } else if (miso_genes_file) {
        ch_miso_genes_file = channel.fromPath(miso_genes_file)
            .splitCsv()
            .set{ ch_miso_genes }
        } else {
        ch_miso_genes = ch_miso_genes_list
    }
    ch_miso_input = MISO_SETTINGS.out.miso_settings.combine(ch_miso_genes)

    ch_bam_bai = ch_bam_join
        .map { it -> [it[1], it[2] ] }
        .collect()

    MISO_SASHIMI (
        ch_miso_index,
        ch_miso_input,
        ch_bam_bai,
        ch_miso_run
    )

    emit:

    gff3                        = GTF_2_GFF3.out.gff3               // path *.gff3
    miso_index                  = ch_miso_index                     // channel: [ ch_miso_index ]
    miso_data                   = ch_miso_run                       // channel: [ ch_miso_run ]
    miso_settings               = MISO_SETTINGS.out.miso_settings   // path miso_setting.txt
    miso_sashimi                = MISO_SASHIMI.out.sashimi          // path sashimi/*.pdf
}
