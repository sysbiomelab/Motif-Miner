#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
BGC Family Analysis Pipeline
Author: Idris Matine
Contact: idris.matine@kcl.ac.uk
Repository: https://github.com/sysbiomelab/Motif-Miner
Version: 1.0.0
This pipeline analyzes a pre-grouped BGC family directory to discover upstream
regulatory motifs, focusing on COGs containing regulatory genes.

*/

include { PREPARE_INPUTS }              from './modules/local/prepare_inputs/main'
include { PROTEINORTHO }                from './modules/local/proteinortho/main'
include { EXTRACT_COGS }                from './modules/local/extract_cogs/main'
include { EXTRACT_AND_FILTER_UPSTREAM } from './modules/local/upstream/main'
include { CDHIT }                       from './modules/local/cdhit/main'
include { GET_CDHIT_STATS }             from './modules/local/get_cdhit_stats/main'
include { GET_GENE_INFO }               from './modules/local/getinfo/main'
include { MEME }                        from './modules/local/meme/main'
include { PARSE_MEME }                  from './modules/local/parsememe/main'
include { ORGANIZE_RESULTS }            from './modules/local/organize_results/main'
include { FINAL_SUMMARY }               from './modules/local/final_summary/main'

def helpMessage() {
    log.info"""
    BGC Family Analysis Pipeline - Version 1.0.0
    ==========================================
    Analyzes a pre-grouped BGC family directory to discover upstream motifs.

    Usage:
        nextflow run ${workflow.scriptFile} --input_dir <path> [options]

    Required arguments:
        --input_dir               Path to the directory containing all BGC family subdirectories.
                                  Each family must be a subdirectory named FAM* containing GenBank files.

    Output options:
        --outdir                  Path to final output directory. (Default: './results')
        --tempdir                 Path for intermediate pipeline files. (Default: './temp')

    Upstream Extraction options:
        --min_upstream_len        Minimum length of extracted sequence to keep. (Default: 100)
        --max_upstream_len        Max length of upstream region to extract from start of gene. (Default: 350)
        --inside_gene_len         Length of sequence to include from within the gene. (Default: 20)

    Classification options:
        --org_threshold           Percentage of a gene kind above which a COG is assigned to that
                                  category (regulatory / biosynthetic / biosynthetic_additional). (Default: 1.0)

    Other options:
        --help                    Display this help message and exit.
        -profile                  Configuration profile to use (standard, hpc, docker; combinable, e.g. standard,docker).
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.input_dir) {
    log.error "ERROR: Missing required argument --input_dir. Use --help for usage."
    exit 1
}

def filtered_families = []

workflow.onComplete {
    def duration_hours = workflow.duration.toHours()
    def duration_mins = workflow.duration.toMinutes() % 60
    
    log.info "\n" + "="*70
    if (workflow.success) {
        log.info "[--- Pipeline Complete ---]"
        log.info "="*70
        log.info "Status          : ✓ SUCCESS"
        log.info "Completed at    : ${workflow.complete}"
        log.info "Duration        : ${duration_hours}h ${duration_mins}m"
        log.info "CPU hours       : ${workflow.stats.computeTimeFmt ?: 'N/A'}"
        log.info "Tasks completed : ${workflow.stats.succeededCount}"
        
        if (filtered_families.size() > 0) {
            log.info "\n" + "-"*70
            log.info "⚠️  IMPORTANT NOTICE"
            log.info "-"*70
            log.info "Number of families excluded from final results: ${filtered_families.size()}"
            log.info "Reason: No motifs discovered by MEME"
            log.info "\nExcluded families:"
            filtered_families.sort().each { fam ->
                log.info "  • ${fam}"
            }
            log.info "\nA complete list has been saved to:"
            log.info "  ${params.outdir}/families_excluded_no_meme_results.txt"
            log.info "-"*70
        } else {
            log.info "\nAll families successfully processed with MEME results."
        }
        
        log.info "\nOutput directory: ${params.outdir}"
        log.info "="*70 + "\n"
        
    } else {
        log.info "[--- Pipeline FAILED ---]"
        log.info "="*70
        log.info "Status       : ✗ FAILED"
        log.info "Completed at : ${workflow.complete}"
        log.info "Duration     : ${duration_hours}h ${duration_mins}m"
        log.info "Error report : ${workflow.errorReport ?: 'See log for details'}"
        log.info "="*70 + "\n"
        log.info "Please check the execution logs for more information.\n"
    }
}

/*
========================================================================================
WORKFLOW DEFINITION
========================================================================================
*/
workflow {

    // --- STEP 0: CREATE AND SPLIT THE INITIAL INPUT CHANNEL ---
    Channel
        .fromPath( "${params.input_dir}/FAM*", type: 'dir' )
        .ifEmpty { exit 1, "No 'FAM*' directories found in input: ${params.input_dir}" }
        .map { dir -> tuple(dir.name, dir) }
        .set { families_ch }

    families_for_prepare = families_ch
    families_for_upstream = families_ch

    // --- STEP 1: PREPARE THE INPUTS ---
    PREPARE_INPUTS(families_for_prepare)

    // --- STEP 2: RUN PROTEINORTHO & EXTRACT COGS ---
    PROTEINORTHO(PREPARE_INPUTS.out.faa_files)
    EXTRACT_COGS(PROTEINORTHO.out.proteinortho_results)

    // --- STEP 3: EXTRACT UPSTREAM SEQUENCES ---
    upstream_input_ch = EXTRACT_COGS.out.cog_files.join(families_for_upstream)
    EXTRACT_AND_FILTER_UPSTREAM(
        upstream_input_ch,
        params.min_upstream_len,
        params.max_upstream_len,
        params.inside_gene_len
    )

    // --- STEP 4: CLUSTER SEQUENCES WITH CD-HIT ---
    CDHIT(EXTRACT_AND_FILTER_UPSTREAM.out.filtered_fastas)
    GET_CDHIT_STATS(CDHIT.out.cdhit_fastas)

    // --- STEP 5: ANALYZE GENE KIND COMPOSITION ---
    gene_info_input_ch = GET_CDHIT_STATS.out.cdhit_fastas_with_stats.join( PREPARE_INPUTS.out.cds_summary )
    GET_GENE_INFO(gene_info_input_ch)

    // --- STEP 6: RUN MEME FOR MOTIF DISCOVERY ---
    meme_input_ch = CDHIT.out.cdhit_fastas.join(GET_GENE_INFO.out.gene_info_tables)
    MEME(meme_input_ch)

    // --- STEP 7: FILTER & PARSE MEME RESULTS ---
    
    MEME.out.meme_families
        .filter { fam_id, meme_dir ->
            meme_dir.list().size() > 0
        }
        .set { ch_meme_results_to_parse }

    MEME.out.meme_families
        .filter { fam_id, meme_dir -> meme_dir.list().size() == 0 }
        .subscribe { fam_id, meme_dir -> 
            filtered_families << fam_id
            log.info "⚠️  Family ${fam_id}: No MEME motifs discovered - excluded from final results"
        }

    MEME.out.meme_families
        .filter { fam_id, meme_dir -> meme_dir.list().size() == 0 }
        .map { fam_id, meme_dir -> fam_id }
        .collectFile(
            name: 'families_excluded_no_meme_results.txt', 
            newLine: true,
            storeDir: params.outdir,
            seed: "# BGC Families Excluded from Final Results\n# Reason: No motifs discovered by MEME\n# Generated: ${new Date()}\n#\n"
        )

    PARSE_MEME(ch_meme_results_to_parse)

    // --- STEP 8: ORGANIZE FINAL RESULTS ---
    final_org_input_ch = PARSE_MEME.out.parsed_meme_families.join(GET_GENE_INFO.out.gene_info_tables)
    ORGANIZE_RESULTS(final_org_input_ch, params.org_threshold)

    // --- STEP 9: CREATE FINAL COMPREHENSIVE SUMMARY ---
    summary_input_ch = ORGANIZE_RESULTS.out.final_results.join( PREPARE_INPUTS.out.cds_summary )
    FINAL_SUMMARY(summary_input_ch)
}
