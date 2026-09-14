// Copies parsed MEME results into gene-kind category folders (regulatory, biosynthetic, biosynthetic_additional, unknown).

process ORGANIZE_RESULTS {
    tag "Organize[${fam_id}]"

    publishDir "${params.outdir}/${fam_id}/", mode: 'copy'

    conda "${projectDir}/conf/conda.yml"

    input:
    tuple val(fam_id), path(parsed_meme_dir), path(gene_info_tsv)
    val threshold

    output:
    tuple val(fam_id), path("organised_results"), emit: final_results

    script:
    """
    source ${projectDir}/bin/helpers.sh
    
    rm -rf organised_results
    mkdir -p organised_results
    log_message "Starting final organization..." 
    python -u ${projectDir}/bin/organize_genekind.py \\
        "${parsed_meme_dir}" \\
        organised_results \\
        "${gene_info_tsv}" \\
        --threshold ${threshold}
    """
}
