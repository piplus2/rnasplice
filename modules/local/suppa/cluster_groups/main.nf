process CLUSTERGROUPS {
    tag "${cond1}-${cond2}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    tuple val(cond1), val(cond2), path(psivec)

    output:
    tuple val(cond1), val(cond2), path(groups), emit: groups
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script:
    groups = "${cond1}-${cond2}_groups.txt"
    """
    suppa_groups.py $psivec > ${groups}
    """

}
