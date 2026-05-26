# nf-core/rnasplice: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Dev v1.0.6dev - TBD

### Added

- Add first steps of leafcutter splicing quantification
- Ignore transcript version in tximport by default to improve compatibility with tx2gene files from GTFs (e.g., from GENCODE) that may not include version numbers in transcript IDs. This can be overridden with `--ignore_tx_version false` if desired.

### Changed

- Synced pipeline with nf-core template 4.0.2
- Migrated samplesheet and contrastsheet validation from Python scripts to nf-schema
- Replaced deprecated `CUSTOM_GETCHROMSIZES` module with `SAMTOOLS_FAIDX`
- Migrated all modules to use `topic: versions` pattern for software version tracking
- Migrated local modules to `TOOL/SUBTOOL/main.nf` directory structure
- Updated `PIPELINE_INITIALISATION` to handle all 4 input source types (fastq, genome_bam, transcriptome_bam, salmon_results)

### Fixed

- Fixed Nextflow 26+ compatibility: removed deprecated `if/else` blocks from `nextflow.config`
- Fixed Nextflow 26+ compatibility: migrated `conf/modules.config` from `if` blocks to `ext.when` closures
- Fixed Nextflow 26+ compatibility: replaced `switch/case` statements with `if/else if` chains
- Fixed channel reuse bugs in SUPPA subworkflow
- Fixed `CREATE_BAMLIST` null path handling for single-condition rMATS runs
- Fixed `SUBREAD_FLATTENGTF` multi-line version output breaking YAML parsing
- Fixed `STAR_ALIGN` deprecated `--quantTranscriptomeBan` parameter renamed to `--quantTranscriptomeSAMoutput`
- Fixed ignored arguments `--clip_r1` and `--clip_r2` in `TRIMGALORE` module for NextSeq trimming and read clipping

## v1.0.5 - 2024-11-03

- Added IsoformSwitchAnalyzeR to pipeline.
- Updated to migrate from nf-validation to nf-schema.
- Updated to remove max_memory, max_cpus, max_time and replace with process resourceLimits.
- Updated for nf-core template version 3.0.2.
- Updated to remove lib folder and add utils subworkflows.

## Dev v1.0.5dev - TBD - TBD

- Add first steps of leafcutter splicing quantification

## v1.0.4 - 2024-04-21

- Fixed incorrect assignment of cluster groups (Issue #131).

## v1.0.3 - 2024-02-23

- Improved TPM file splitting performance (Issue #120).
- Fixed an issue where R scripts altered sample names upon loading (Issue #122).

## v1.0.2 - 2024-01-08

Patch for run_stager.R (#108) and template update v2.11.1 (#109).

## v1.0.1 - 2023-11-15

Patch for run_drimseq_filter.R to cast command line arguments to numeric. See issue #98 on nf-core/rnasplice.

## v1.0.0 - 2023-05-22

First release of nf-core/rnasplice, created with the [nf-core](https://nf-co.re/) template.

### `Added`

Implemented pipeline:

- Merge re-sequenced FastQ files (cat)
- Read QC (FastQC)
- Adapter and quality trimming (TrimGalore)
- Alignment with STAR:
  - STAR -> Salmon
  - STAR -> featureCounts
  - STAR -> HTSeq (DEXSeq count)
- Sort and index alignments (SAMtools)
- Create bigWig coverage files (BEDTools, bedGraphToBigWig)
- Pseudo-alignment and quantification (Salmon; optional)
- Summarize QC (MultiQC)
- Differential Exon Usage (DEU):
  - HTSeq -> DEXSeq
  - featureCounts -> edgeR
  - Quantification with featureCounts or HTSeq
- Differential exon usage with DEXSeq or edgeR
  - Differential Transcript Usage (DTU):
  - Salmon -> DRIMSeq -> DEXSeq
  - Filtering with DRIMSeq
- Differential transcript usage with DEXSeq
- Event-based splicing analysis:
  - STAR -> rMATS
  - Salmon -> SUPPA2

Updated pipeline:

- Visualization of differential results with edgeR, DEXSeq, and MISO
- Contrasts specified using contrastsheet.csv
- Allow users to specify input data type and start point (e.g., fastq, genome_bam, transcript_bam, salmon_results)
- Pipeline schematic updated
