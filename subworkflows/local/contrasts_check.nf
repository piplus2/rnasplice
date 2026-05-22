//
// Validate and parse input contrastsheet
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow CONTRASTS_CHECK {

    take:
    contrastsheet // file: /path/to/contrastsheet.csv

    main:

    channel
        .fromList(samplesheetToList(contrastsheet, "${projectDir}/assets/schema_contrasts.json"))
        .map { meta ->
            return [
                contrast:  meta.contrast,
                treatment: meta.treatment,
                control:   meta.control
            ]
        }
        .set { ch_contrasts }

    emit:
    contrasts = ch_contrasts                        // channel: [ val(contrast_map) ]
    versions = channel.empty()       // channel: [ versions.yml ]
}
