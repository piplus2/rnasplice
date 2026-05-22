process DEXSEQ_EXON {
    label 'process_high'

    conda "bioconda::bioconductor-dexseq=1.36.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dexseq:1.36.0--r40_0' :
        'biocontainers/bioconductor-dexseq:1.36.0--r40_0' }"

    input:
    path ("dexseq_clean_counts/*")    // path: dexseq_clean_counts
    path gff                          // path: dexseq_gff
    path samplesheet                  // path: samplesheet
    path contrastsheet                // path: contrastsheet
    val ntop                          // val: n_dexseq_plot

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "DEXSeqResults.*.csv"  , emit: dexseq_exon_results_csv
    path "perGeneQValue.*.csv"  , emit: dexseq_gene_results_csv
    path "plotDEXSeq.*.pdf"     , emit: dexseq_plot_results_pdf
    tuple val("${task.process}"), val('r-base'), eval('R --version 2>&1 | head -n 1 | sed "s/^.*version //; s/ .*$//"'), topic: versions, emit: versions_R
    tuple val("${task.process}"), val('bioconductor-dexseq'), eval('Rscript -e "library(DEXSeq); cat(as.character(packageVersion(\'DEXSeq\')))"'), topic: versions, emit: versions_dexseq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_dexseq_exon.R dexseq_clean_counts $gff $samplesheet $contrastsheet $ntop
    """
}
