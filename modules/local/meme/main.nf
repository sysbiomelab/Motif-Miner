// Runs MEME to perform motif discovery on upstream sequences that meet the count threshold.

process MEME {
    tag "MEME [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/7_meme_analysis", mode: 'copy'

    container 'docker://memesuite/memesuite:5.5.8'

    input:
    tuple val(fam_id), path(cdhit_fasta_dir), path(gene_info_tsv)

    output:
    tuple val(fam_id), path("meme_output"), emit: meme_families

    script:
    """
    source ${projectDir}/bin/helpers.sh
    rm -rf meme_output
    mkdir -p meme_output

    tail -n +2 ${gene_info_tsv} | while IFS=\$'\\t' read -r family cog total rest; do
        if [[ "\$total" -gt 4 ]]; then
            fasta_file="${cdhit_fasta_dir}/\${cog}.fasta"

            if [[ -f "\$fasta_file" ]]; then
                log_message "Running MEME on \${fasta_file} with \$total sequences..."
                
                meme "\$fasta_file" -oc "meme_output/\${cog}_anr_maxw30_\${total}seqs"  -dna -mod anr   -maxw 30 -nmotifs 10 -revcomp -p ${task.cpus} < /dev/null
                meme "\$fasta_file" -oc "meme_output/\${cog}_zoops_maxw30_\${total}seqs" -dna -mod zoops -maxw 30 -nmotifs 10 -revcomp -p ${task.cpus} < /dev/null
                meme "\$fasta_file" -oc "meme_output/\${cog}_anr_maxw20_\${total}seqs"  -dna -mod anr   -maxw 20 -nmotifs 10 -revcomp -p ${task.cpus} < /dev/null
                meme "\$fasta_file" -oc "meme_output/\${cog}_zoops_maxw20_\${total}seqs"  -dna -mod zoops   -maxw 20 -nmotifs 10 -revcomp -p ${task.cpus} < /dev/null
                meme "\$fasta_file" -oc "meme_output/\${cog}_anr_maxw10_\${total}seqs"  -dna -mod anr   -maxw 10 -nmotifs 10 -revcomp -p ${task.cpus} < /dev/null
                meme "\$fasta_file" -oc "meme_output/\${cog}_zoops_maxw10_\${total}seqs" -dna -mod zoops -maxw 10 -nmotifs 10 -revcomp -p ${task.cpus} < /dev/null

            else
                log_message "WARNING: FASTA file not found for \${cog}, skipping."
            fi
        fi
    done

    log_message "-----------------------"
    log_message "-----------------------"
    log_message "-----------------------"
    log_message "All Done"
    log_message "-----------------------"
    log_message "-----------------------"
    log_message "-----------------------"

    """
}
