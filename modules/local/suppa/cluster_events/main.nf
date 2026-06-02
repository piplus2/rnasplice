process CLUSTEREVENTS {
    tag "${cond1}-${cond2}"
    label 'process_high'
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    tuple val(cond1), val(cond2), path(dpsi), path(psivec), val(group_ranges) // e.g. 1-3,4-6
    val prefix
    val clusterevents_dpsithreshold    // val params.clusterevents_dpsithreshold
    val clusterevents_eps              // val params.clusterevents_eps
    val clusterevents_metric           // val params.clusterevents_metric
    val clusterevents_min_pts          // val params.clusterevents_min_pts
    val clusterevents_method           // val params.clusterevents_method
    val clusterevents_sigthreshold     // val params.clusterevents_sigthreshold
    val clusterevents_separation       // val params.clusterevents_separation

    output:
    path "*.clustvec"   , emit: clustvec
    path "*.log"        , emit: cluster_log
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script: //  Cluster events between conditions

    def thresh_arg  = clusterevents_sigthreshold ? "-st ${clusterevents_sigthreshold}" : ''
    def separation_arg = clusterevents_separation ? "-s ${clusterevents_separation}" : ''

    """
    touch ${cond1}-${cond2}_${prefix}_cluster.clustvec

    suppa.py \\
        clusterEvents \\
        --dpsi ${dpsi} \\
        --psivec ${psivec} \\
        --dpsi-threshold ${clusterevents_dpsithreshold} \\
        --eps ${clusterevents_eps} \\
        --metric ${clusterevents_metric} \\
        --min-pts ${clusterevents_min_pts} \\
        --groups \$(cat ${group_ranges}) \\
        --clustering ${clusterevents_method} \\
        ${thresh_arg} ${separation_arg} -o ${cond1}-${cond2}_${prefix}_cluster
    """

}
