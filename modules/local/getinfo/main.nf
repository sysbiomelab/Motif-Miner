// Classifies genes in each clustered COG using the master CDS summary file,
// producing a table of gene kind percentages for the family.

process GET_GENE_INFO {
    tag "Gene Info [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/6_gene_info", mode: 'copy'

    conda "${projectDir}/conf/conda.yml"

    input:
    tuple val(fam_id), path(cdhit_fasta_dir), path(cds_summary)

    output:
    tuple val(fam_id), path("${fam_id}_gene_analysis_results.tsv"), emit: gene_info_tables

    script:
    def output_filename = "${fam_id}_gene_analysis_results.tsv"
    """
    source ${projectDir}/bin/helpers.sh
    
    log_message "Analyzing gene composition for ${fam_id}..."

    python -u ${projectDir}/bin/get_info.py \\
        --family_dir "${cdhit_fasta_dir}" \\
        --familyname "${fam_id}" \
        --reference_file "${cds_summary}" \\
        --output "${output_filename}"
    
    log_message "Gene info analysis complete."

    """
}
