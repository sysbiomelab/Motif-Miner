// Extracts COG sequences from the ProteinOrtho results.

process EXTRACT_COGS {
    tag "Extract COGs [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/3_proteinortho_cogs", mode: 'copy'

    conda "${projectDir}/conf/conda.yml"

    input:
    tuple val(fam_id), path(faa_dir), path(proteinortho_tsv)

    output:
    tuple val(fam_id), path("COG_files"), emit: cog_files

    path "*.tsv", emit: cog_stats, optional: true    

    script:
    """
    source ${projectDir}/bin/helpers.sh
    log_message "Extracting COG sequences for ${fam_id}..."

    python -u ${projectDir}/bin/extract_cog_proteins.py extract \\
        "${proteinortho_tsv}" \\
        "${faa_dir}" \\
        --output_dir COG_files \\
        --stats_output .

    log_message "COG extraction complete."
    """
}
