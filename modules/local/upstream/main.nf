// Extracts upstream intergenic regions based on COG/GBK files and filters them by length.

process EXTRACT_AND_FILTER_UPSTREAM {
    tag "Upstream Extractor [${fam_id}]"
    publishDir "${params.tempdir}/${fam_id}/4_upstream_filtered", mode: 'copy'

    input:
    tuple val(fam_id), path(cog_dir), path(gbk_dir)

    val min_len
    val max_len
    val inside_len

    output:
    tuple val(fam_id), path("upstream_fastas"), emit: filtered_fastas
    path "*_summary.tsv", emit: summary_tsv, optional: true

    script:
    """
    source ${projectDir}/bin/helpers.sh
    
    log_message "Extracting and filtering upstream sequences for ${fam_id}..."
    rm -rf upstream_fastas
    mkdir -p upstream_fastas
    
    python -u ${projectDir}/bin/upstream_extractor.py \
        --gbk-path "${gbk_dir}" \
        --cog-path "${cog_dir}" \
        --output-path "upstream_fastas" \
        --stats-output . \
        --min-length ${min_len} \
        --upstream ${max_len} \
        --inside ${inside_len}

    log_message "Upstream extraction complete."    
    """
}
