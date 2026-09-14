// For a single BGC family, this module:
// 1. Creates the CDS summary TSV file.
// 2. Converts all GBK files to protein FAA files.

process PREPARE_INPUTS {
    tag "Pre-processing [${fam_id}]"

    publishDir "${params.tempdir}/${fam_id}/1_cds_summary", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.tempdir}/${fam_id}/2_faa_files", mode: 'copy', pattern: "faa_files"

    conda "${projectDir}/conf/conda.yml"

    input:
    tuple val(fam_id), path(fam_dir)

    output:
    tuple val(fam_id), path("faa_files") , emit: faa_files
    tuple val(fam_id), path("${fam_id}_cds_summary.tsv"), emit: cds_summary

    script:
    def summary_filename = "${fam_id}_cds_summary.tsv"
    """
    source ${projectDir}/bin/helpers.sh
    python -u ${projectDir}/bin/gbk_to_faa.py ${fam_dir} --output-dir faa_files --include-partial
    log_message "-----------------------------------"
    log_message "gbk_to_faa DONE"
    log_message "-----------------------------------"
    python -u ${projectDir}/bin/gbk_pre.py ${fam_dir} . ${summary_filename}
    log_message "-----------------------------------"
    log_message "DONE"
    log_message "-----------------------------------"
    """
}
