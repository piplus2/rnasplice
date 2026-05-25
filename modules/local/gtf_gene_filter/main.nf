process GTF_GENE_FILTER {
    tag "$fasta"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta2), path("*.gtf")       , emit: gtf
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script: // filter_gtf_for_genes_in_genome.py is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    filter_gtf_for_genes_in_genome.py \\
        --gtf $gtf \\
        --fasta $fasta \\
        -o ${fasta.baseName}_genes.gtf
    """
}
