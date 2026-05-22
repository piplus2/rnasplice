//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA                          } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF                            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF                            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA               } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF_DEXSEQ                     } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_SUPPA_TPM                      } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_STAR_INDEX                       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX                     } from '../../../modules/nf-core/untar'

include { CUSTOM_GETCHROMSIZES                            } from '../../../modules/nf-core/custom/getchromsizes'
include { GFFREAD                                         } from '../../../modules/nf-core/gffread'
include { STAR_GENOMEGENERATE                             } from '../../../modules/nf-core/star/genomegenerate'
include { STAR_GENOMEGENERATE_IGENOMES                    } from '../../../modules/local/star_genomegenerate_igenomes'
include { SALMON_INDEX                                    } from '../../../modules/nf-core/salmon/index'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../../modules/nf-core/rsem/preparereference'

include { GTF_GENE_FILTER                                 } from '../../../modules/local/gtf_gene_filter'
include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE            } from '../../../modules/local/preprocess_transcripts_fasta_gencode'

workflow PREPARE_GENOME {
    take:
    fasta            //      file: /path/to/genome.fasta
    gtf              //      file: /path/to/genome.gtf
    gff              //      file: /path/to/genome.gff
    transcript_fasta //      file: /path/to/transcript.fasta
    star_index       // directory: /path/to/star/index/
    salmon_index     // directory: /path/to/salmon/index/
    gff_dexseq       //      file: /path/to/dexseq/genome.gff
    suppa_tpm        //      file: /path/to/suppa/quant.tpm
    gencode          //   boolean: whether gene annotation is from gencode

    main:

    ch_versions = channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA([[:], fasta])
        ch_fasta    = GUNZIP_FASTA.gunzip.map { it -> it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            GUNZIP_GTF([[:], gtf])
            ch_gtf      = GUNZIP_GTF.out.gunzip.map { it -> it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = channel.value(file(gtf))
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            GUNZIP_GFF([[:], gff])
            ch_gff      = GUNZIP_GFF.out.gunzip.map { it -> it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = channel.value(file(gff))
        }
        ch_gtf      = GFFREAD(ch_gff, null).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            GUNZIP_TRANSCRIPT_FASTA([[:], transcript_fasta])
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA.out.gunzip.map { it -> it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = channel.value(file(transcript_fasta))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(ch_transcript_fasta)
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        ch_filter_gtf       = GTF_GENE_FILTER(ch_fasta, ch_gtf).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA(ch_fasta, ch_filter_gtf).transcript_fasta
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES(ch_fasta.map { it -> [[:], it] })
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map   { it -> it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it -> it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress STAR index or generate from scratch if required
    //
    def is_aws_igenome = params.fasta && params.gtf &&
        (file(params.fasta).getName() - '.gz' == 'genome.fa') &&
        (file(params.gtf).getName() - '.gz' == 'genes.gtf')

    ch_star_index = channel.empty()
    if (params.source == 'fastq' && !params.skip_alignment && (params.aligner == 'star' || params.aligner == 'star_salmon')) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX([[:], star_index]).untar.map { it -> it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = channel.value(file(star_index))
            }
        } else {
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES(ch_fasta, ch_gtf).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE(ch_fasta, ch_gtf).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = channel.empty()
    if (params.source == 'fastq' && (params.pseudo_aligner == 'salmon' || params.aligner == 'star_salmon')) {
        if (salmon_index) {
            if (salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR_SALMON_INDEX([[:], salmon_index]).untar.map { it -> it[1] }
                ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
            } else {
                ch_salmon_index = channel.value(file(salmon_index))
            }
        } else {
            if (params.pseudo_aligner == 'salmon') {
                SALMON_INDEX(ch_fasta, ch_transcript_fasta)
                ch_salmon_index = SALMON_INDEX.out.index
                ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
            }
        }
    }

    //
    // Uncompress DEXSeq GFF annotation file if required
    //
    ch_dexseq_gff = channel.empty()
    if (gff_dexseq) {
        if (gff_dexseq.endsWith('.gz')) {
            GUNZIP_GFF_DEXSEQ([[:], gff_dexseq])
            ch_dexseq_gff = GUNZIP_GFF_DEXSEQ.out.gunzip.map { it -> it[1] }
            ch_versions   = ch_versions.mix(GUNZIP_GFF_DEXSEQ.out.versions)
        } else {
            ch_dexseq_gff = channel.value(file(gff_dexseq))
        }
    }

    //
    // Uncompress SUPPA TPM file if required
    //
    ch_suppa_tpm = channel.empty()
    if (suppa_tpm) {
        if (suppa_tpm.endsWith('.gz')) {
            GUNZIP_SUPPA_TPM([[:], suppa_tpm])
            ch_suppa_tpm = GUNZIP_SUPPA_TPM.out.gunzip.map { it -> it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_SUPPA_TPM.out.versions)
        } else {
            ch_suppa_tpm = channel.value(file(suppa_tpm))
        }
    }

    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    fai              = ch_fai              //    path: genome.fai
    chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    gtf              = ch_gtf              //    path: genome.gtf
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    star_index       = ch_star_index       //    path: star/index/
    salmon_index     = ch_salmon_index     //    path: salmon/index/
    dexseq_gff       = ch_dexseq_gff       //    path: dexseq.gff
    suppa_tpm        = ch_suppa_tpm        //    path: suppa.tpm
    is_aws_igenome   = is_aws_igenome     //    boolean: whether the genome files are from AWS iGenomes
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
