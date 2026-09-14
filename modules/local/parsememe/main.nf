// Parses the XML output from MEME into user-friendly CSV and PWM files.

process PARSE_MEME {
    tag "Parse MEME [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/8_meme_analysis_parsed", mode: 'copy'

    conda "${projectDir}/conf/conda.yml"

    input:
    tuple val(fam_id), path(meme_dir)

    output:
    tuple val(fam_id), path("parse_meme"), emit: parsed_meme_families

    script:
    """
    source ${projectDir}/bin/helpers.sh
    log_message "Preparing to parse MEME results for ${fam_id}..."

    rm -rf parse_meme
    mkdir -p parse_meme
    
    cp -r ${meme_dir}/* parse_meme/
    
    log_message "Parsing MEME files"

    python -u ${projectDir}/bin/parse_meme.py parse_meme

    log_message "MEME parsing complete."
    """
}
