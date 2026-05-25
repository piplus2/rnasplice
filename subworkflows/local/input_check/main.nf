//
// Subworkflow to check the input samplesheet configuration and emit structured streams
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    take:
    samplesheet // string: path to input samplesheet matching assets/schema_input.json
    source      //    val: execution input format type profile [fastq, genome_bam, transcriptome_bam, salmon_results]

    main:

    // Select schema matching the declared input source type
    def source_schemas = [
        fastq:             "${projectDir}/assets/schema_input_fastq.json",
        genome_bam:        "${projectDir}/assets/schema_input_genome_bam.json",
        transcriptome_bam: "${projectDir}/assets/schema_input_transcriptome_bam.json",
        salmon_results:    "${projectDir}/assets/schema_input_salmon_results.json"
    ]

    def raw_samplesheet_list = samplesheetToList(samplesheet, source_schemas[source])

    channel
        .fromList(raw_samplesheet_list)
        .branch {
            fastq:             source == 'fastq'
            genome_bam:        source == 'genome_bam'
            transcriptome_bam: source == 'transcriptome_bam'
            salmon_results:    source == 'salmon_results'
            none:              true
        }
        .set { ch_branched_inputs }

    // ====================================================================
    // BRANCH TRANSFORMATIONS
    // ====================================================================

    // 1. Standard FASTQ structural lane generation
    ch_reads_fastq = ch_branched_inputs.fastq
        .map { meta, fastq_1, fastq_2 ->
            def meta_map = [
                id:           meta.id,
                single_end:   !fastq_2,
                strandedness: meta.strandedness,
                condition:    meta.condition
            ]
            def files = meta_map.single_end ? [ file(fastq_1) ] : [ file(fastq_1), file(fastq_2) ]
            return [ meta_map, files ]
        }

    // 2. Alignment mapping profile (Genome BAM)
    ch_reads_genome = ch_branched_inputs.genome_bam
        .map { meta, genome_bam ->
            def meta_map = [
                id:        meta.id,
                condition: meta.condition
            ]
            return [ meta_map, [ file(genome_bam) ] ]
        }

    // 3. Transcript quantification profile (Transcriptome BAM)
    // Samplesheet has both genome_bam and transcriptome_bam columns; only
    // transcriptome_bam is used here
    ch_reads_transcriptome = ch_branched_inputs.transcriptome_bam
        .map { meta, _genome_bam, transcriptome_bam ->
            def meta_map = [
                id:        meta.id,
                condition: meta.condition
            ]
            return [ meta_map, [ file(transcriptome_bam) ] ]
        }

    // 4. Pre-calculated Salmon quant ingestion target structures
    ch_reads_salmon = ch_branched_inputs.salmon_results
        .map { meta, salmon_dir ->
            def meta_map = [
                id:        meta.id,
                condition: meta.condition
            ]
            return [ meta_map, [ file(salmon_dir) ] ]
        }

    // Synthesize structured output indices based on selected mode
    ch_final_reads = channel.empty()
        .mix(ch_reads_fastq)
        .mix(ch_reads_genome)
        .mix(ch_reads_transcriptome)
        .mix(ch_reads_salmon)

    emit:
    reads    = ch_final_reads   // channel: [ val(meta), [ files ] ]
}
