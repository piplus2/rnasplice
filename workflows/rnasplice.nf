/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FORWARD                } from '../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REVERSE                } from '../modules/nf-core/bedtools/genomecov'
include { ISOFORMSWITCHANALYZER                                           } from '../modules/local/isoformswitchanalyzer'
include { ALIGN_STAR                                                      } from '../subworkflows/local/align_star'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_SALMON                     } from '../subworkflows/local/tx2gene_tximport'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_STAR_SALMON                } from '../subworkflows/local/tx2gene_tximport'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_SALMON                 } from '../subworkflows/local/drimseq_dexseq_dtu'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_STAR_SALMON            } from '../subworkflows/local/drimseq_dexseq_dtu'
include { RMATS                                                           } from '../subworkflows/local/rmats'
include { DEXSEQ_DEU                                                      } from '../subworkflows/local/dexseq_deu'
include { EDGER_DEU                                                       } from '../subworkflows/local/edger_deu'
include { SUPPA as SUPPA_SALMON                                           } from '../subworkflows/local/suppa'
include { SUPPA as SUPPA_STAR_SALMON                                      } from '../subworkflows/local/suppa'
include { VISUALISE_MISO                                                  } from '../subworkflows/local/visualize_miso'
include { LEAFCUTTER                                                      } from '../subworkflows/local/leafcutter'

include { validateInputSamplesheet                                        } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { validateInputContrastsheet                                      } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { rmatsReadError                                                  } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { rmatsStrandednessError                                          } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { multiqcTsvFromList                                              } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { paramsSummaryMultiqc                                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                          } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList                                               } from 'plugin/nf-schema'
include { FASTQC                                                          } from '../modules/nf-core/fastqc/main'
include { SALMON_QUANT as SALMON_QUANT_SALMON                             } from '../modules/nf-core/salmon/quant'
include { SALMON_QUANT as SALMON_QUANT_STAR                               } from '../modules/nf-core/salmon/quant'
include { MULTIQC                                                         } from '../modules/nf-core/multiqc/'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                     } from '../modules/nf-core/custom/dumpsoftwareversions'
include { CAT_FASTQ                                                       } from '../modules/nf-core/cat/fastq'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE                                } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BAM_SORT_STATS_SAMTOOLS                                         } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { softwareVersionsToYAML                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                                                } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT SOURCES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// One entry per accepted value of `--source`: the schema its samplesheet is
// validated against, and how a validated row is normalised to [ meta, [ files ] ].
// Supporting a new input type means adding one entry here and nothing else.
//
def inputSources() {
    return [
        fastq: [
            schema: "${projectDir}/assets/schema_input.json",
            reads: { row ->
                def (meta, fastq_1, fastq_2) = row
                [meta + [single_end: !fastq_2], fastq_2 ? [fastq_1, fastq_2] : [fastq_1]]
            },
        ],
        genome_bam: [
            schema: "${projectDir}/assets/schema_input_genome_bam.json",
            reads: { row ->
                def (meta, genome_bam) = row
                [[id: meta.id, condition: meta.condition], [genome_bam]]
            },
        ],
        transcriptome_bam: [
            // The samplesheet carries both a genome and a transcriptome BAM column;
            // only the transcriptome BAM is consumed downstream.
            schema: "${projectDir}/assets/schema_input_transcriptome_bam.json",
            reads: { row ->
                def (meta, _genome_bam, transcriptome_bam) = row
                [[id: meta.id, condition: meta.condition], [transcriptome_bam]]
            },
        ],
        salmon_results: [
            schema: "${projectDir}/assets/schema_input_salmon_results.json",
            reads: { row ->
                def (meta, salmon_results) = row
                [[id: meta.id, condition: meta.condition], [salmon_results]]
            },
        ],
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASPLICE {
    take:
    ch_samplesheet // channel: file(samplesheet)
    ch_contrastsheet // channel: file(contrastsheet)
    ch_fasta // channel: path of genome fasta
    ch_gtf // channel: path of genome gtf
    ch_transcript_fasta // channel: path of transcript fasta
    ch_dexseq_gff // channel: path of dexseq gff
    ch_salmon_index // channel: path of salmon index
    ch_star_index // channel: path of star index
    ch_suppa_tpm // channel: path of user-supplied suppa tpm
    ch_chrom_sizes // channel: path of chromosome sizes

    main:

    def ch_multiqc_files = channel.empty()
    def pass_trimmed_reads = [:]

    //
    // Resolve the routes through the pipeline that the input options imply, once.
    // Every step below tests one of these flags instead of re-deriving the same
    // compound `--source`/`--aligner` conditions.
    //
    def from_fastq = params.source == 'fastq'
    def from_genome_bam = params.source == 'genome_bam'
    def from_transcriptome_bam = params.source == 'transcriptome_bam'
    def from_salmon_results = params.source == 'salmon_results'

    // Genome BAMs are either read from the samplesheet or produced by STAR
    def align_with_star = from_fastq && !params.skip_alignment && params.aligner in ['star', 'star_salmon']
    def has_genome_bam = from_genome_bam || align_with_star

    // Salmon is run either against a transcriptome BAM or directly against reads
    def quant_from_bam = from_transcriptome_bam || (from_fastq && !params.skip_alignment && params.aligner == 'star_salmon')
    def quant_from_reads = from_fastq && params.pseudo_aligner == 'salmon'
    def has_salmon_results = quant_from_reads || from_salmon_results

    def input_sources = inputSources()

    if (!input_sources.containsKey(params.source)) {
        error("Invalid --source parameter: '${params.source}'. Must be one of: ${input_sources.keySet().join(', ')}")
    }

    //
    // Read the samplesheet with the schema for the declared `--source` and stage
    // its files as [ meta, [ files ] ], whatever the input type
    //
    def input_source = input_sources[params.source]

    ch_reads = channel
        .fromList(samplesheetToList(params.input, input_source.schema))
        .map(input_source.reads)

    //
    // Validate the contrastsheet: each contrast must compare two distinct conditions
    //
    channel
        .fromList(samplesheetToList(params.contrasts, "${projectDir}/assets/schema_contrasts.json"))
        .map { contrast -> validateInputContrastsheet(contrast) }

    // Check rMATS parameter configuration mapping checks
    if (params.rmats && from_fastq) {
        rmatsReadError(ch_reads)
        rmatsStrandednessError(ch_reads)
    }

    //
    // MODULE: Merge the runs of each sample, then concatenate FastQ technical
    // replicates (ids suffixed `_T1`, `_T2`, ...) if required
    //
    if (from_fastq) {
        ch_reads
            .map { meta, fastqs -> [meta.id, meta, fastqs] }
            .groupTuple()
            .map { samplesheet -> validateInputSamplesheet(samplesheet) }
            .map { meta, fastqs -> [meta + [id: meta.id - ~/_T\d+/], fastqs.flatten()] }
            .groupTuple()
            .branch { meta, fastq ->
                single: fastq.size() == 1
                return [meta, fastq.flatten()]
                multiple: fastq.size() > 1
                return [meta, fastq.flatten()]
            }
            .set { ch_fastq }

        CAT_FASTQ(ch_fastq.multiple)
        CAT_FASTQ.out.reads.mix(ch_fastq.single).set { ch_cat_fastq }

        //
        // SUBWORKFLOW: Read QC and trim adapters with TrimGalore!
        //
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
            ch_cat_fastq,
            params.skip_fastqc,
            false,
            true,
            params.skip_trimming,
            0,
            0,
        )
        ch_trim_reads = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count

        // Parse custom evaluation warnings if reading maps fall below specific metric criteria
        ch_trim_read_count
            .map { meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return ["${meta.id}\t${num_reads}"]
                }
            }
            .collect()
            .map { tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_trimming_multiqc }
    }

    //
    // SUBWORKFLOW: Sort, index and QC BAMs read straight from the samplesheet
    //
    if (from_genome_bam || from_transcriptome_bam) {
        BAM_SORT_STATS_SAMTOOLS(ch_reads, ch_fasta)
        ch_samtools_stats = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_samtools_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_samtools_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats

        if (from_genome_bam) {
            ch_genome_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
            ch_genome_bam_index = BAM_SORT_STATS_SAMTOOLS.out.index
        }
        else {
            // Salmon quantifies against the transcriptome BAM as supplied
            ch_transcriptome_bam_for_salmon = ch_reads
        }
    }

    //
    // SUBWORKFLOW: Alignment with STAR
    //
    if (align_with_star) {
        ALIGN_STAR(
            ch_trim_reads,
            ch_star_index.map { index -> [[:], index] },
            ch_gtf.map { gtf -> [[:], gtf] },
            ch_fasta.map { fasta -> [[:], fasta, []] },
        )
        ch_genome_bam = ALIGN_STAR.out.bam
        ch_genome_bam_index = ALIGN_STAR.out.index
        ch_transcriptome_bam_for_salmon = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc = ALIGN_STAR.out.log_final
    }

    // ====================================================================
    // EXON SPLICING ENGINES
    // ====================================================================

    if (has_genome_bam) {

        if (params.dexseq_exon) {
            DEXSEQ_DEU(
                ch_gtf,
                ch_genome_bam,
                ch_dexseq_gff,
                ch_samplesheet,
                ch_contrastsheet,
            )
        }

        if (params.edger_exon) {
            EDGER_DEU(
                ch_gtf,
                ch_genome_bam,
                ch_samplesheet,
                ch_contrastsheet,
            )
        }

        if (params.rmats) {
            RMATS(
                ch_samplesheet,
                ch_contrastsheet,
                ch_genome_bam.map { meta, bam -> [meta.condition, meta, bam] },
                ch_gtf,
            )
        }

        if (params.sashimi_plot) {
            VISUALISE_MISO(
                ch_gtf,
                ch_genome_bam,
                ch_genome_bam_index,
            )
        }

        if (params.leafcutter) {
            LEAFCUTTER(ch_genome_bam, ch_genome_bam_index, ch_gtf)
        }
    }

    // ====================================================================
    // TRANSCRIPT ALIGNMENT & QUANTIFICATION ENGINES (SALMON)
    // ====================================================================

    if (quant_from_bam) {
        ch_dummy_salmon_index = channel.value(file("${projectDir}/assets/dummy_file.txt", checkIfExists: true))

        SALMON_QUANT_STAR(
            ch_transcriptome_bam_for_salmon,
            ch_dummy_salmon_index,
            ch_gtf,
            ch_transcript_fasta,
            true,
            params.salmon_quant_libtype ?: '',
        )
        ch_salmon_results = SALMON_QUANT_STAR.out.results

        TX2GENE_TXIMPORT_STAR_SALMON(ch_salmon_results, ch_gtf)

        if (params.dexseq_dtu) {
            DRIMSEQ_DEXSEQ_DTU_STAR_SALMON(
                params.dtu_txi == "dtuScaledTPM" ? TX2GENE_TXIMPORT_STAR_SALMON.out.txi_dtu : TX2GENE_TXIMPORT_STAR_SALMON.out.txi_s,
                TX2GENE_TXIMPORT_STAR_SALMON.out.tximport_tx2gene,
                ch_samplesheet,
                ch_contrastsheet,
            )
        }

        if (params.suppa) {
            // Get Suppa tpm either from tximport or user supplied
            def ch_star_salmon_suppa_tpm = params.suppa_tpm ? ch_suppa_tpm : TX2GENE_TXIMPORT_STAR_SALMON.out.suppa_tpm

            SUPPA_STAR_SALMON(
                ch_gtf.map { gtf -> [[id: gtf.baseName], gtf] },
                ch_star_salmon_suppa_tpm.map { tpm_psi -> [[id: tpm_psi.baseName], tpm_psi] },
                ch_samplesheet,
                ch_contrastsheet,
            )
        }
    }

    if (quant_from_reads) {
        ch_dummy_transcript_fasta = channel.value(file("${projectDir}/assets/dummy_file2.txt", checkIfExists: true))

        SALMON_QUANT_SALMON(
            ch_trim_reads,
            ch_salmon_index,
            ch_gtf,
            ch_dummy_transcript_fasta,
            false,
            params.salmon_quant_libtype ?: '',
        )
        ch_salmon_results = SALMON_QUANT_SALMON.out.results
    }
    else if (from_salmon_results) {
        ch_salmon_results = ch_reads
    }

    if (has_salmon_results) {
        TX2GENE_TXIMPORT_SALMON(ch_salmon_results, ch_gtf)

        if (params.dexseq_dtu) {
            DRIMSEQ_DEXSEQ_DTU_SALMON(
                params.dtu_txi == "dtuScaledTPM" ? TX2GENE_TXIMPORT_SALMON.out.txi_dtu : TX2GENE_TXIMPORT_SALMON.out.txi_s,
                TX2GENE_TXIMPORT_SALMON.out.tximport_tx2gene,
                ch_samplesheet,
                ch_contrastsheet,
            )
        }

        if (params.suppa) {
            // Get Suppa tpm either from tximport or user supplied
            def ch_salmon_suppa_tpm = params.suppa_tpm ? ch_suppa_tpm : TX2GENE_TXIMPORT_SALMON.out.suppa_tpm

            SUPPA_SALMON(
                ch_gtf.map { gtf -> [[id: gtf.baseName], gtf] },
                ch_salmon_suppa_tpm.map { tpm_psi -> [[id: tpm_psi.baseName], tpm_psi] },
                ch_samplesheet,
                ch_contrastsheet,
            )
        }
    }

    if (params.isoformswitchanalyzer && (from_fastq || from_salmon_results)) {
        ISOFORMSWITCHANALYZER(
            ch_salmon_results.collect { it -> it[1] },
            ch_gtf,
            ch_transcript_fasta,
            ch_samplesheet,
            ch_contrastsheet,
            params.isoformswitchanalyzer_alpha,
            params.isoformswitchanalyzer_dIF,
        )
    }

    //
    // Module: Generate stranded coverage bigWig files for visualisation in genome browsers
    //
    if (!params.skip_bigwig && (from_genome_bam || (from_fastq && !params.skip_alignment))) {
        ch_sizes = ch_chrom_sizes.map { _meta, sizes -> sizes }
        ch_genome_bam_scaled = ch_genome_bam.map { meta, bam -> [meta, bam, 0] }

        BEDTOOLS_GENOMECOV_FORWARD(ch_genome_bam_scaled, ch_sizes, 'bedGraph', true)
        BEDTOOLS_GENOMECOV_REVERSE(ch_genome_bam_scaled, ch_sizes, 'bedGraph', true)

        BEDGRAPH_TO_BIGWIG_FORWARD(BEDTOOLS_GENOMECOV_FORWARD.out.genomecov, ch_sizes)
        BEDGRAPH_TO_BIGWIG_REVERSE(BEDTOOLS_GENOMECOV_REVERSE.out.genomecov, ch_sizes)
    }

    // ====================================================================
    // REPORTING LAYER & VERSION COLLATOR SYNC
    // ====================================================================

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(topic_versions_string).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'nf_core_' + 'rnasplice_software_' + 'mqc_' + 'versions.yml',
        sort: true,
        newLine: true,
    )

    // Build static summary models for MultiQC ingestion space
    def summary_params = paramsSummaryMap(workflow)
    def workflow_summary = paramsSummaryMultiqc(summary_params)
    ch_workflow_summary = channel.value(workflow_summary)

    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def methods_description = methodsDescriptionText(ch_multiqc_custom_methods_description)
    ch_methods_description = channel.value(methods_description)

    // Blend final channel maps
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    if (from_fastq) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]))
    }

    if (quant_from_reads) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_salmon_results.collect { it -> it[1] }.ifEmpty([]))
    }

    if (align_with_star) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_multiqc.collect { it -> it[1] }.ifEmpty([]))
    }

    if (from_genome_bam || from_transcriptome_bam || align_with_star) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_idxstats.collect { it -> it[1] }.ifEmpty([]))

        if (params.edger_exon && !from_transcriptome_bam) {
            ch_multiqc_files = ch_multiqc_files.mix(EDGER_DEU.out.featureCounts_summary.collect { it -> it[1] }.ifEmpty([]))
        }
    }

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [[id: 'rnasplice'], files, params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true), params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [], [], []]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
}
