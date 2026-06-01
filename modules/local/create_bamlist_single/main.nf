process CREATE_BAMLIST_SINGLE {
    label 'process_single'

    conda "conda-forge::sed=4.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'biocontainers/sed:4.7.0' }"

    input:
    tuple val(contrast), val(cond1), path(bam1)

    output:
    tuple val(contrast), path("${cond1}_bamlist.txt"), emit: bamlist
    tuple val("${task.process}"), val('sed'), eval('sed --version 2>&1 | head -1 | sed "s/sed (GNU sed) //g"'), topic: versions, emit: versions_sed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo ${bam1} | sed 's: :,:g' > ${cond1}_bamlist.txt
    """
}
