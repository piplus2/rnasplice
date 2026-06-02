process MQC_RMATS {
    tag "${contrast}"

    input:
    tuple val(contrast), path(summary)

    output:
    path "${contrast}_summary.rmats.txt"

    script:
    """
    awk -v contrast=${contrast} \\
        'NR==1 {print "Contrast_Event\\tContrast\\t" \$0} \\
        NR>1 {print contrast"_"\$1 "\\t" contrast "\\t" \$0}' \\
        ${summary} > ${contrast}_summary.rmats.txt
    """
}