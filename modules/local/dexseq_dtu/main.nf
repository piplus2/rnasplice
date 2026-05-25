process DEXSEQ_DTU {
    label 'process_high'

    conda "bioconda::bioconductor-dexseq=1.36.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dexseq:1.36.0--r40_0' :
        'biocontainers/bioconductor-dexseq:1.36.0--r40_0' }"

    input:
    path drimseq_sample_data
    path drimseq_d_counts
    path drimseq_contrast_data
    val ntop

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "DEXSeqResults.*.tsv"  , emit: dexseq_exon_results_tsv
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "perGeneQValue.*.tsv"  , emit: dexseq_gene_results_tsv
    tuple val("${task.process}"), val('r-base'), eval('R --version 2>&1 | head -n 1 | sed "s/^.*version //; s/ .*$//"'), topic: versions, emit: versions_R
    tuple val("${task.process}"), val('bioconductor-dexseq'), eval('Rscript -e "library(DEXSeq); cat(as.character(packageVersion(\'DEXSeq\')))"'), topic: versions, emit: versions_dexseq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_dexseq_dtu.R $drimseq_sample_data $drimseq_contrast_data $drimseq_d_counts $ntop
    """

}
