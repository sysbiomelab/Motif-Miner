#!/usr/bin/env python3
"""
Get stats after the CD-HIT is done for a single family
"""

import os
import sys
import logging
import argparse
import csv
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

def setup_logging(level=logging.INFO) -> logging.Logger:
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def validate_args(family_path: str, output_dir: str):
    """Validate command line arguments."""
    logger = logging.getLogger(__name__)
    
    if not os.path.exists(family_path):
        logger.error(f"Family directory does not exist: {family_path}")
        sys.exit(1)
    
    if not os.path.isdir(family_path):
        logger.error(f"Family path is not a directory: {family_path}")
        sys.exit(1)
    
    if not os.path.exists(output_dir):
        logger.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.isdir(output_dir):
        logger.error(f"Output path is not a directory: {output_dir}")
        sys.exit(1)

def extract_cog_name(filename: str) -> str:
    """Extract COG name from filename"""
    base_name = Path(filename).stem
    match = re.match(r'^(COG_\d+)', base_name)
    if match:
        return match.group(1)
    else:
        parts = base_name.split('_')
        return '_'.join(parts[:2]) if len(parts) >= 2 else base_name

def get_stat(fasta_file: str) -> Dict[str, Any]:
    """Analyze a FASTA file and return statistics."""
    logger = logging.getLogger(__name__)
    
    input_filename = Path(fasta_file).name
    input_filename_base = input_filename.replace("_cdhit.fasta", "")
    
    cog_name = extract_cog_name(input_filename_base)
    
    stats = {
        'cog_name': cog_name,
        'input_file': input_filename,
        'total_sequences': 0,
    }
    
    try:
        with open(fasta_file, 'r') as f:
            sequence_count = 0
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    sequence_count += 1
        
        stats['total_sequences'] = sequence_count
        logger.info(f"Processed {input_filename}: {sequence_count} sequences")
        
    except Exception as e:
        logger.error(f"Error processing file {fasta_file}: {str(e)}")
        stats['total_sequences'] = 0
    
    return stats

def main():
    logger = setup_logging()
    
    parser = argparse.ArgumentParser(description='Analyze CD-HIT FASTA files for a single gene family')
    parser.add_argument('--family_dir', '-f', required=True, 
                        help='Path to the family directory containing FASTA files')
    parser.add_argument('--output_dir', '-d', default='.',
                        help='Directory to save the output TSV file (default: current directory)')
    parser.add_argument('--output_name', '-o', default='cdhit_stats.tsv',
                        help='Name of the output TSV file (default: cdhit_stats.tsv)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging (DEBUG level)')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")
    
    logger.info("Starting single family CD-HIT statistics analysis")
    logger.info(f"Family directory: {args.family_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Output filename: {args.output_name}")
    
    validate_args(args.family_dir, args.output_dir)
    
    family_name = os.path.basename(args.family_dir.rstrip('/'))
    logger.info(f"Processing family: {family_name}")
    
    family_results = []
    
    column_order = [
        'cog_name', 'input_file', 'total_sequences'
    ]
    
    fasta_files = []
    for file_path in Path(args.family_dir).iterdir():
        if file_path.is_file() and file_path.suffix.lower() in ['.fasta', '.fa']:
            fasta_files.append(file_path)
    
    if not fasta_files:
        logger.error(f"No FASTA files found in {family_name}")
        return 1
    
    logger.info(f"Found {len(fasta_files)} FASTA files in {family_name}")
    
    for i, file_path in enumerate(fasta_files, 1):
        logger.info(f"Processing file {i}/{len(fasta_files)}: {file_path.name}")
        
        stats = get_stat(str(file_path))
        family_results.append(stats)
    
    if family_results:
        family_df = pd.DataFrame(family_results, columns=column_order)
        
        family_df = family_df.sort_values('cog_name')
        
        output_path = os.path.join(args.output_dir, args.output_name)
        
        family_df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved results to: {output_path}")
        
        logger.info("=" * 50)
        logger.info("CD-HIT STATISTICS SUMMARY")
        logger.info("=" * 50)
        logger.info(f"Family: {family_name}")
        logger.info(f"Total files processed: {len(family_results)}")
        
        total_sequences = family_df['total_sequences'].sum()
        avg_sequences = family_df['total_sequences'].mean()
        min_sequences = family_df['total_sequences'].min()
        max_sequences = family_df['total_sequences'].max()
        
        logger.info(f"Total sequences across all files: {total_sequences:,}")
        logger.info(f"Average sequences per file: {avg_sequences:.1f}")
        logger.info(f"Minimum sequences in a file: {min_sequences}")
        logger.info(f"Maximum sequences in a file: {max_sequences}")
        
        logger.info(f"\nAll {len(family_results)} file results:")
        for _, row in family_df.iterrows():
            logger.info(f"  {row['cog_name']}: {row['total_sequences']} sequences")
        
        logger.info("CD-HIT statistics analysis completed successfully!")
        return 0
    else:
        logger.error("No results generated. Please check your input files.")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
