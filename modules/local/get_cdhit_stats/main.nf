// Generates a statistics summary of the CD-HIT clustering results.

process GET_CDHIT_STATS {
    tag "CD-HIT Stats [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/5_cdhit", mode: 'copy', pattern: '*.tsv'
    
    conda "${projectDir}/conf/conda.yml"
    
    input:
    tuple val(fam_id), path(cdhit_fasta_dir)

    output:
    tuple val(fam_id), path(cdhit_fasta_dir), emit: cdhit_fastas_with_stats
    path "*_cdhit_stats.tsv", emit: cdhit_stats, optional: true

    script:
    def stats_filename = "${fam_id}_cdhit_stats.tsv"
    """
    source ${projectDir}/bin/helpers.sh

    log_message "Generating CD-HIT statistics..."
    python -u ${projectDir}/bin/cd_stats.py \\
        --family_dir "${cdhit_fasta_dir}" \\
        --output_name "${stats_filename}"
    log_message "Statistics generation complete."

    """
}
