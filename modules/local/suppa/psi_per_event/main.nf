process PSIPEREVENT {
    tag "$tpm"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    path ioe
    path tpm
    val psiperevent_total_filter   // val params.psiperevent_total_filter

    output:
    path "suppa_local.psi"    , emit: psi
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script: // Calculate the psi values of local events

    """
    suppa.py \\
        psiPerEvent \\
        -i $ioe \\
        -e $tpm \\
        -f $psiperevent_total_filter \\
        -o suppa_local
    """
}
