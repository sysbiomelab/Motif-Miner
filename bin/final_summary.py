#!/usr/bin/env python3
"""
Final Summary Generator for BGC Family Analysis Pipeline
=========================================================
Author: Idris Matine
Contact: idris.matine@kcl.ac.uk

This script consolidates motif discovery results with gene annotations
to create a comprehensive summary table for each BGC family.
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
import re


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate comprehensive summary table combining motif data with gene annotations'
    )
    parser.add_argument(
        '--fam_id',
        required=True,
        help='Family identifier (e.g., FAM9209)'
    )
    parser.add_argument(
        '--results_dir',
        required=True,
        type=Path,
        help='Path to the organized results directory for this family'
    )
    parser.add_argument(
        '--cds_summary',
        required=True,
        type=Path,
        help='Path to the CDS summary TSV file'
    )
    parser.add_argument(
        '--output',
        required=True,
        type=Path,
        help='Output file path for the summary table'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    return parser.parse_args()


def log_message(message, verbose=False):
    """Print log message if verbose is enabled."""
    if verbose:
        print(f"[INFO] {message}", file=sys.stderr)


def read_cds_summary(cds_file, verbose=False):
    """Read the CDS summary file."""
    log_message(f"Reading CDS summary from {cds_file}", verbose)
    
    try:
        cds_df = pd.read_csv(cds_file, sep='\t')
        log_message(f"Loaded {len(cds_df)} genes from CDS summary", verbose)
        return cds_df
    except Exception as e:
        print(f"ERROR: Failed to read CDS summary: {e}", file=sys.stderr)
        sys.exit(1)


def find_motif_files(results_dir, verbose=False):
    """Find all combined_sequence_motif_binding.csv files in the results directory."""
    log_message(f"Searching for motif files in {results_dir}", verbose)
    
    motif_files = []
    
    if not results_dir.exists():
        print(f"ERROR: Results directory does not exist: {results_dir}", file=sys.stderr)
        sys.exit(1)
    
    for category_dir in results_dir.iterdir():
        if not category_dir.is_dir():
            continue
        
        cog_type = category_dir.name
        
        for cog_dir in category_dir.iterdir():
            if not cog_dir.is_dir():
                continue
            
            cog_match = re.match(r'(COG_\d+)_', cog_dir.name)
            if not cog_match:
                log_message(f"Skipping directory (no COG pattern): {cog_dir.name}", verbose)
                continue
            
            cog_name = cog_match.group(1)
            
            csv_file = cog_dir / "combined_sequence_motif_binding.csv"
            if csv_file.exists():
                motif_files.append((cog_type, cog_name, csv_file, cog_dir))
                log_message(f"Found motif file: {cog_type}/{cog_name}", verbose)
            else:
                log_message(f"No motif file in: {cog_dir.name}", verbose)
    
    log_message(f"Found {len(motif_files)} motif files total", verbose)
    return motif_files


def read_motif_file(csv_file, verbose=False):
    """Read a motif binding CSV file."""
    try:
        # The files are COMMA-separated, not tab-separated
        df = pd.read_csv(csv_file, sep=',')
        if verbose and len(df) > 0:
            log_message(f"  Loaded {len(df)} rows with columns: {', '.join(df.columns[:5])}...", verbose)
        return df
    except Exception as e:
        log_message(f"Warning: Failed to read {csv_file}: {e}", verbose)
        return pd.DataFrame()


def get_gene_info(gene_name, cds_df):
    """Get gene information from CDS summary for a given gene name."""
    gene_info = cds_df[cds_df['Locus_Tag'] == gene_name]
    
    if len(gene_info) == 0:
        return {
            'File': '',
            'Start': '',
            'End': '',
            'Strand': '',
            'Gene_Functions': '',
            'Gene_Kind': '',
            'Protein_ID': '',
            'Product': '',
            'Translation': ''
        }
    else:
        return gene_info.iloc[0].to_dict()


def process_motif_data(fam_id, motif_files, cds_df, results_dir, verbose=False):
    """Process all motif files and combine with gene information."""
    all_data = []
    
    for cog_type, cog_name, csv_file, cog_dir in motif_files:
        log_message(f"Processing {cog_type}/{cog_name}", verbose)
        
        motif_df = read_motif_file(csv_file, verbose)
        
        if motif_df.empty:
            log_message(f"Skipping empty motif file: {csv_file}", verbose)
            continue
        
        try:
            relative_dir = cog_dir.relative_to(results_dir)
        except ValueError:
            relative_dir = cog_dir.name
        
        for _, motif_row in motif_df.iterrows():
            gene_name = motif_row.get('gene_name', '')
            
            gene_info_dict = get_gene_info(gene_name, cds_df)
            
            row_data = {
                'Family_ID': fam_id,
                'COG_Name': cog_name,
                'COG_Type': cog_type,
                'Gene_Name': gene_name,
                'Source_File': gene_info_dict.get('File', ''),
                'Gene_Start': gene_info_dict.get('Start', ''),
                'Gene_End': gene_info_dict.get('End', ''),
                'Gene_Strand': gene_info_dict.get('Strand', ''),
                'Gene_Functions': gene_info_dict.get('Gene_Functions', ''),
                'Gene_Kind': gene_info_dict.get('Gene_Kind', ''),
                'Protein_ID': gene_info_dict.get('Protein_ID', ''),
                'Product': gene_info_dict.get('Product', ''),
                'Translation': gene_info_dict.get('Translation', ''),
                'Motif_Directory': str(relative_dir),
                'Seq_ID': motif_row.get('seq_id', ''),
                'Sequence_Name': motif_row.get('sequence_name', ''),
                'Seq_Length': motif_row.get('seq_len', ''),
                'Seq_Strand': motif_row.get('seq_strand', ''),
                'Seq_Weight': motif_row.get('seq_weight', ''),
                'Seq_Pvalue': motif_row.get('seq_pvalue', ''),
                'Num_Sites_In_Seq': motif_row.get('num_sites_in_seq', ''),
                'Motif_ID': motif_row.get('motif_id', ''),
                'Motif_Sequence': motif_row.get('motif_seq', ''),
                'Motif_Width': motif_row.get('motif_width', ''),
                'Motif_Evalue': motif_row.get('motif_evalue', ''),
                'Motif_Start_Position': motif_row.get('start_position', ''),
                'Motif_End_Position': motif_row.get('end_position', ''),
                'Motif_Relative_Position_Percent': motif_row.get('relative_position_percent', ''),
                'Distance_To_Gene': motif_row.get('distance_to_gene', ''),
                'Seq_Motif_Pvalue': motif_row.get('seq_motif_pvalue', '')
            }
            
            all_data.append(row_data)
    
    log_message(f"Processed {len(all_data)} total motif occurrences", verbose)
    return pd.DataFrame(all_data)


def write_summary_table(df, output_file, verbose=False):
    """Write the final summary table to file."""
    if df.empty:
        log_message("No data to write, creating empty file with headers", verbose)
        
        headers = [
            'Family_ID', 'COG_Name', 'COG_Type', 'Gene_Name', 'Source_File',
            'Gene_Start', 'Gene_End', 'Gene_Strand', 'Gene_Functions', 'Gene_Kind',
            'Protein_ID', 'Product', 'Translation', 'Motif_Directory',
            'Seq_ID', 'Sequence_Name', 'Seq_Length', 'Seq_Strand', 'Seq_Weight', 'Seq_Pvalue',
            'Num_Sites_In_Seq', 'Motif_ID', 'Motif_Sequence', 'Motif_Width',
            'Motif_Evalue', 'Motif_Start_Position', 'Motif_End_Position',
            'Motif_Relative_Position_Percent', 'Distance_To_Gene', 'Seq_Motif_Pvalue'
        ]
        df = pd.DataFrame(columns=headers)
    else:
        df = df.sort_values(['COG_Name', 'Gene_Name'])
        log_message(f"Writing {len(df)} rows to {output_file}", verbose)
    
    try:
        df.to_csv(output_file, sep='\t', index=False)
        log_message(f"Successfully wrote summary to {output_file}", verbose)
    except Exception as e:
        print(f"ERROR: Failed to write output file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main execution function."""
    args = parse_args()
    
    log_message(f"Starting final summary generation for {args.fam_id}", args.verbose)
    log_message(f"Results directory: {args.results_dir}", args.verbose)
    log_message(f"CDS summary: {args.cds_summary}", args.verbose)
    log_message(f"Output file: {args.output}", args.verbose)
    
    cds_df = read_cds_summary(args.cds_summary, args.verbose)
    
    motif_files = find_motif_files(args.results_dir, args.verbose)
    
    if not motif_files:
        log_message("No motif files found, creating empty summary", args.verbose)
        write_summary_table(pd.DataFrame(), args.output, args.verbose)
        return
    
    final_df = process_motif_data(
        args.fam_id,
        motif_files,
        cds_df,
        args.results_dir,
        args.verbose
    )
    
    write_summary_table(final_df, args.output, args.verbose)
    
    log_message("Final summary generation complete!", args.verbose)


if __name__ == '__main__':
    main()
