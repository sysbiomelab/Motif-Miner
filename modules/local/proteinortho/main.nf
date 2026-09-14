// Runs ProteinOrtho on the family's protein FASTA files.

process PROTEINORTHO {
    tag "ProteinOrtho [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/3_proteinortho_raw", mode: 'copy'

    container 'docker://quay.io/biocontainers/proteinortho:6.3.5--h2b77389_1'

    input:
    tuple val(fam_id), path(faa_dir)

    output:
    tuple val(fam_id), path("safe_faa_copy"), path("*.proteinortho.tsv"), emit: proteinortho_results
    
    path "*", emit: all_reports_for_publishing, optional: true
    
    script:
    """
    source ${projectDir}/bin/helpers.sh

    log_message "Starting main ProteinOrtho analysis for ${fam_id}..."

    rm -rf safe_faa_copy
    log_message "Making the directory safe_faa_copy..."
    mkdir -p safe_faa_copy

    log_message "Copying the faa files to safe_faa_copy..."
    cp ${faa_dir}/*.faa safe_faa_copy/
    log_message "Done copying."

    log_message "Running Proteinortho"
    proteinortho \\
        -project="${fam_id}" \\
        safe_faa_copy/*.faa \\
        -cpus=${task.cpus} \\
        -verbose=1
    log_message "Main ProteinOrtho analysis complete."

    """
}
