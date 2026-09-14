// Clusters upstream DNA sequences using CD-HIT-EST and generates a stats summary.

process CDHIT {
    tag "CD-HIT [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/5_cdhit", mode: 'copy'

    container 'docker://nanozoo/cdhit:4.8.1--c697693'

    input:
    tuple val(fam_id), path(upstream_fasta_dir)

    output:
    tuple val(fam_id), path("filtered_cdhit"), emit: cdhit_fastas

    script:
    """
    source ${projectDir}/bin/helpers.sh
    rm -rf filtered_cdhit
    mkdir -p filtered_cdhit

    log_message "Clustering upstream sequences for ${fam_id}..."

    for fasta_file in ${upstream_fasta_dir}/*.fasta; do
        base=\$(basename "\$fasta_file" .fasta)

        cd-hit-est \\
            -i "\$fasta_file" \\
            -o "filtered_cdhit/\${base}.fasta" \\
            -c 0.90 \\
            -n 5 \\
            -d 0 \\
            -T ${task.cpus}

    done

    log_message "CD-HIT clustering complete."
    """
}
