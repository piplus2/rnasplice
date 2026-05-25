process GTF_2_GFF3 {
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_3' :
        'biocontainers/gffread:0.12.7--hdcf5f25_3' }"

    input:
    path gtf

    output:
    path "*.gff3"           , emit: gff3
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gffread $gtf -L --keep-genes | awk -F'\\t' -vOFS='\\t' '{ gsub("transcript", "mRNA", \$3); print}' > ${gtf.baseName}_genes.gff3
    """
}
