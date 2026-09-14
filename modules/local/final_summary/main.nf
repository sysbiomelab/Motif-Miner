
process FINAL_SUMMARY {
    tag "${fam_id}"
    publishDir "${params.outdir}/${fam_id}", mode: 'copy'
    conda "${projectDir}/conf/conda.yml"
    input:
    tuple val(fam_id), path(final_results_dir), path(cds_summary)

    output:
    tuple val(fam_id), path("${fam_id}_complete_summary.tsv"), emit: final_summary

    script:
    """
    source ${projectDir}/bin/helpers.sh
    
    log_message "Creating final summary table for ${fam_id}..."

    python -u ${projectDir}/bin/final_summary.py \\
        --fam_id ${fam_id} \\
        --results_dir ${final_results_dir} \\
        --cds_summary ${cds_summary} \\
        --output ${fam_id}_complete_summary.tsv \\
        --verbose
    log_message "Final summary complete."
    log_message "Done"
    """
}
