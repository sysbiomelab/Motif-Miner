"""
Combined Upstream Sequence Extractor and Filter

This script extracts upstream sequences from genes in GenBank files based on 
COG (Clusters of Orthologous Groups) annotations and filters them based on 
minimum length requirements. It extracts intergenic regions between CDS features
and allows for flexible upstream and inside sequence lengths.

Usage:
python upstream_extractor.py \
  --gbk-path "/path/to/genbank/files/" \
  --cog-path "/path/to/cog/files/" \
  --output-path "/path/to/output/" \
  --stats-output "/path/to/stats/" \
  --min-length 100 \
  --upstream 350 \
  --inside 20
"""

import os
import sys
import logging
import argparse
import csv
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

def validate_paths(gbk_path: str, cog_path: str, output_path: str, stats_output: str) -> bool:
    """Validate that required paths exist."""
    logger = logging.getLogger(__name__)
    
    if not os.path.exists(gbk_path):
        logger.error(f"GenBank path does not exist: {gbk_path}")
        return False
    
    if not os.path.exists(cog_path):
        logger.error(f"COG path does not exist: {cog_path}")
        return False
    
    os.makedirs(output_path, exist_ok=True)
    logger.info(f"Output directory ready: {output_path}")
    
    os.makedirs(stats_output, exist_ok=True)
    logger.info(f"Stats output directory ready: {stats_output}")
    
    return True

def parse_header(header: str) -> Optional[Tuple[str, str, int, int, int]]:
    """Parse FASTA header to extract gene information."""
    logger = logging.getLogger(__name__)
    
    try:
        parts = header.split()
        if len(parts) < 3:
            logger.warning(f"Invalid header format: {header}")
            return None
        
        file_id = parts[1]
        info = parts[2]
        parts_info = info.split('|')
        
        if len(parts_info) < 3:
            logger.warning(f"Invalid info format in header: {header}")
            return None
        
        locus_tag = parts_info[0]
        start_str, end_str = parts_info[1].split('-')
        start, end = int(start_str), int(end_str)
        strand = int(parts_info[2])
        
        return file_id, locus_tag, start, end, strand
    
    except (ValueError, IndexError) as e:
        logger.error(f"Error parsing header '{header}': {e}")
        return None

def get_cds_features(record: SeqRecord) -> Tuple[List[Tuple[int, int, int, str]], List[Tuple[int, int, int, str]]]:
    """Extract all CDS features from a GenBank record, separated by strand."""
    forward_cds = []
    reverse_cds = []
    
    for feature in record.features:
        if feature.type == "CDS":
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            
            locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
            
            cds_tuple = (start, end, strand, locus_tag)
            
            if strand == 1:
                forward_cds.append(cds_tuple)
            else:
                reverse_cds.append(cds_tuple)
    
    forward_cds.sort(key=lambda x: x[0])
    reverse_cds.sort(key=lambda x: x[0])
    
    return forward_cds, reverse_cds

def find_previous_cds_same_strand(target_start: int, target_end: int, target_strand: int, 
                                 forward_cds: List[Tuple[int, int, int, str]], 
                                 reverse_cds: List[Tuple[int, int, int, str]]) -> Optional[Tuple[int, int, int, str]]:
    """Find the previous CDS on the same strand."""
    if target_strand == 1:
        # Forward strand: find previous CDS that ends before target starts
        previous_cds = None
        for start, end, strand, locus_tag in forward_cds:
            if end <= target_start:
                previous_cds = (start, end, strand, locus_tag)
            else:
                break
        return previous_cds
    else:
        # For reverse strand, "previous" means downstream in genomic coordinates
        for start, end, strand, locus_tag in reverse_cds:
            if start >= target_end:
                return (start, end, strand, locus_tag)
        return None

def extract_upstream_sequence(record: SeqRecord, target_start: int, target_end: int, 
                            target_strand: int, prev_cds: Optional[Tuple[int, int, int, str]],
                            upstream_length: int, inside_length: int) -> Optional[Seq]:
    """Extract the upstream sequence for a target CDS with specified lengths."""
    logger = logging.getLogger(__name__)
    seq_length = len(record.seq)
    
    if target_strand == 1:
        if prev_cds is None:
            # First CDS on forward strand - extract from position 0 to target start + inside_length
            upstream_start = max(0, target_start + inside_length - upstream_length)
            upstream_end = target_start + inside_length
            logger.debug(f"First forward CDS: extracting from {upstream_start} to {upstream_end}")
        else:
            prev_start, prev_end, prev_strand, prev_locus = prev_cds
            upstream_start = max(prev_end, target_start + inside_length - upstream_length)
            upstream_end = target_start + inside_length
            logger.debug(f"Forward CDS: extracting from {upstream_start} to {upstream_end}")
        
            
        upstream_seq = record.seq[upstream_start:upstream_end]
        
    else:
        if prev_cds is None:
            # Last CDS on reverse strand - extract from target_end - inside_length backwards
            upstream_start = target_end - inside_length
            upstream_end = min(seq_length, target_end - inside_length + upstream_length)
            logger.debug(f"Last reverse CDS: extracting from {upstream_start} to {upstream_end}")
        else:
            prev_start, prev_end, prev_strand, prev_locus = prev_cds
            upstream_start = target_end - inside_length
            upstream_end = min(prev_start, target_end - inside_length + upstream_length)
            logger.debug(f"Reverse CDS: extracting from {upstream_start} to {upstream_end}")
        

        upstream_seq = record.seq[upstream_start:upstream_end].reverse_complement()
    
    # Return the extracted sequence, which might be shorter than upstream_length
    return upstream_seq

def extract_cog_name(filename: str) -> str:
    """Extract COG name from filename."""
    base_name = filename.replace(".fasta", "").replace(".fa", "").replace(".fas", "")
    return base_name

def process_cog_file(cog_file_path: str, gbk_path: str, min_length: int, 
                    upstream_length: int, inside_length: int) -> Tuple[List[SeqRecord], Dict[str, Any]]:
    """Process a single COG FASTA file and extract upstream sequences."""
    logger = logging.getLogger(__name__)
    upstream_records = []
    
    cog_name = extract_cog_name(os.path.basename(cog_file_path))
    
    stats = {
        'cog_name': cog_name,
        'total_genes': 0,
        'missing_gbk': 0,
        'no_previous_cds': 0,
        'insufficient_length': 0,
        'sequences_kept': 0,
        'parse_errors': 0,
        'processing_errors': 0
    }
    
    try:
        with open(cog_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if not line.startswith('>'):
                    continue
                
                stats['total_genes'] += 1
                header = line.strip()[1:]
                parsed = parse_header(header)
                
                if parsed is None:
                    stats['parse_errors'] += 1
                    logger.warning(f"Skipping line {line_num} due to parsing error")
                    continue
                
                file_id, locus_tag, start, end, strand = parsed
                
                gbk_file = file_id.replace('.faa', '.gbk')
                gbk_full_path = os.path.join(gbk_path, gbk_file)
                
                if not os.path.isfile(gbk_full_path):
                    stats['missing_gbk'] += 1
                    logger.warning(f"Missing GenBank file: {gbk_file}")
                    continue
                
                try:
                    logger.debug(f"Reading {gbk_file}")
                    record = next(SeqIO.parse(gbk_full_path, "genbank"))
                    
                    forward_cds, reverse_cds = get_cds_features(record)
                    logger.debug(f"Found {len(forward_cds)} forward and {len(reverse_cds)} reverse CDS features")
                    
                    prev_cds = find_previous_cds_same_strand(start, end, strand, forward_cds, reverse_cds)
                    
                    upstream_seq = extract_upstream_sequence(record, start, end, strand, prev_cds,
                                                           upstream_length, inside_length)
                    
                    if upstream_seq is None:
                        stats['insufficient_length'] += 1
                        logger.debug(f"Insufficient upstream sequence for {locus_tag}")
                        continue
                    
                    if len(upstream_seq) < min_length:
                        stats['insufficient_length'] += 1
                        logger.debug(f"Sequence too short for {locus_tag}: {len(upstream_seq)} < {min_length}")
                        continue
                    
                    if prev_cds is None:
                        if strand == 1:
                            seq_type = "sequence_start_to_gene"
                        else:
                            seq_type = "gene_to_sequence_end"
                    else:
                        seq_type = "intergenic"
                    
                    clean_id = f"{locus_tag}_upstream"
                    if prev_cds is None:
                        description = f"{file_id} {locus_tag}|{start}-{end}|{strand} {seq_type} length={len(upstream_seq)}bp upstream={upstream_length}bp inside={inside_length}bp"
                    else:
                        prev_locus = prev_cds[3]
                        description = f"{file_id} {locus_tag}|{start}-{end}|{strand} {seq_type} length={len(upstream_seq)}bp upstream={upstream_length}bp inside={inside_length}bp prev_CDS={prev_locus}"
                    
                    new_record = SeqRecord(
                        upstream_seq,
                        id=clean_id,
                        description=description
                    )
                    
                    upstream_records.append(new_record)
                    stats['sequences_kept'] += 1
                    logger.debug(f"Added upstream sequence for {locus_tag} ({len(upstream_seq)} bp, {seq_type})")
                
                except Exception as e:
                    stats['processing_errors'] += 1
                    logger.error(f"Error processing {gbk_file}: {e}")
                    continue
        
        print(f"\nStatistics for {cog_name}:")
        print(f"  Total genes processed: {stats['total_genes']}")
        print(f"  Sequences kept: {stats['sequences_kept']}")
        print(f"  Missing GenBank files: {stats['missing_gbk']}")
        print(f"  Insufficient length: {stats['insufficient_length']}")
        print(f"  Header parsing errors: {stats['parse_errors']}")
        print(f"  Processing errors: {stats['processing_errors']}")
        
        if stats['sequences_kept'] == 0:
            if stats['missing_gbk'] > 0:
                logger.warning(f"No sequences extracted - {stats['missing_gbk']} missing GenBank files")
            elif stats['insufficient_length'] > 0:
                logger.warning(f"No sequences extracted - {stats['insufficient_length']} sequences too short")
            else:
                logger.warning("No sequences extracted - check file formats and paths")
        
    except FileNotFoundError:
        logger.error(f"COG file not found: {cog_file_path}")
    except Exception as e:
        logger.error(f"Error reading COG file {cog_file_path}: {e}")
    
    return upstream_records, stats

def write_cog_summary_tsv(tsv_path: str, all_stats: List[Dict[str, Any]], min_length: int, 
                         upstream_length: int, inside_length: int):
    """Write COG summary to TSV file."""
    logger = logging.getLogger(__name__)
    
    try:
        os.makedirs(os.path.dirname(tsv_path), exist_ok=True)
        
        with open(tsv_path, 'w', newline='', encoding='utf-8') as tsvfile:
            fieldnames = [
                'cog_name', 'total_genes', 'sequences_kept', 'missing_gbk', 
                'insufficient_length', 'parse_errors', 'processing_errors', 'retention_rate_percent'
            ]
            
            writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
            
            writer.writerow({
                'cog_name': f'# Parameters: min_length={min_length}bp, upstream={upstream_length}bp, inside={inside_length}bp',
                'total_genes': '', 'sequences_kept': '', 'missing_gbk': '', 
                'insufficient_length': '', 'parse_errors': '', 'processing_errors': '', 'retention_rate_percent': ''
            })
            
            writer.writeheader()
            
            total_genes_all_cogs = 0
            total_kept_all_cogs = 0
            
            for stats in all_stats:
                retention_rate = (stats['sequences_kept'] / stats['total_genes'] * 100) if stats['total_genes'] > 0 else 0.0
                
                writer.writerow({
                    'cog_name': stats['cog_name'],
                    'total_genes': stats['total_genes'],
                    'sequences_kept': stats['sequences_kept'],
                    'missing_gbk': stats['missing_gbk'],
                    'insufficient_length': stats['insufficient_length'],
                    'parse_errors': stats['parse_errors'],
                    'processing_errors': stats['processing_errors'],
                    'retention_rate_percent': round(retention_rate, 2)
                })
                
                total_genes_all_cogs += stats['total_genes']
                total_kept_all_cogs += stats['sequences_kept']
            
            overall_retention = (total_kept_all_cogs / total_genes_all_cogs * 100) if total_genes_all_cogs > 0 else 0.0
            writer.writerow({
                'cog_name': 'TOTAL',
                'total_genes': total_genes_all_cogs,
                'sequences_kept': total_kept_all_cogs,
                'missing_gbk': sum(s['missing_gbk'] for s in all_stats),
                'insufficient_length': sum(s['insufficient_length'] for s in all_stats),
                'parse_errors': sum(s['parse_errors'] for s in all_stats),
                'processing_errors': sum(s['processing_errors'] for s in all_stats),
                'retention_rate_percent': round(overall_retention, 2)
            })
                
        logger.info(f"📊 Wrote COG summary to {tsv_path}")
        logger.info(f"Summary contains {len(all_stats)} COGs with min_length≥{min_length}bp, upstream={upstream_length}bp, inside={inside_length}bp")
        
    except Exception as e:
        logger.error(f"Error writing COG summary TSV file {tsv_path}: {e}")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Extract and filter upstream sequences from genes based on COG annotations. '
                   'Extracts intergenic regions between CDS features with flexible length parameters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--gbk-path', 
        required=True,
        help='Path to directory containing GenBank (.gbk) files'
    )
    
    parser.add_argument(
        '--cog-path', 
        required=True,
        help='Path to directory containing COG FASTA files'
    )
    
    parser.add_argument(
        '--output-path', 
        required=True,
        help='Path to output directory for upstream sequences'
    )
    
    parser.add_argument(
        '--stats-output',
        help='Path to directory for saving statistics summary file (default: same as output-path)'
    )
    
    parser.add_argument(
        '--min-length',
        type=int,
        default=100,
        help='Minimum sequence length to keep (sequences shorter than this will be discarded)'
    )
    
    parser.add_argument(
        '--upstream',
        type=int,
        default=350,
        help='Length of upstream sequence to extract (bp)'
    )
    
    parser.add_argument(
        '--inside',
        type=int,
        default=20,
        help='Length from inside the gene to include (bp)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    return parser.parse_args()

def main():
    """Main function to process all COG files."""
    args = parse_arguments()
    
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = setup_logging(log_level)
    logger.info("Starting upstream sequence extraction and filtering")
    
    stats_output = args.stats_output if args.stats_output else args.output_path
    logger.info(f"Statistics will be saved to: {stats_output}")
    
    if not validate_paths(args.gbk_path, args.cog_path, args.output_path, stats_output):
        logger.error("Path validation failed. Exiting.")
        sys.exit(1)
    
    if args.min_length <= 0:
        logger.error(f"Minimum length must be positive, got: {args.min_length}")
        sys.exit(1)
    
    if args.upstream <= 0:
        logger.error(f"Upstream length must be positive, got: {args.upstream}")
        sys.exit(1)
    
    if args.inside < 0:
        logger.error(f"Inside length must be non-negative, got: {args.inside}")
        sys.exit(1)
    
    logger.info(f"Parameters: min_length={args.min_length}bp, upstream={args.upstream}bp, inside={args.inside}bp")
    
    cog_files = [f for f in os.listdir(args.cog_path) if f.endswith('.fasta')]
    
    if not cog_files:
        logger.warning("No FASTA files found in COG directory")
        return
    
    logger.info(f"Found {len(cog_files)} COG files to process")
    
    total_sequences = 0
    all_stats = []
    
    for cog_file in cog_files:
        logger.info(f"Processing: {cog_file}")
        
        cog_file_path = os.path.join(args.cog_path, cog_file)
        upstream_records, stats = process_cog_file(cog_file_path, args.gbk_path, 
                                                  args.min_length, args.upstream, args.inside)
        
        all_stats.append(stats)
        
        if upstream_records:
            base_name = os.path.splitext(cog_file)[0]
            output_file = os.path.join(args.output_path, f'{base_name}_upstream_{args.upstream}bp_inside_{args.inside}bp_min_{args.min_length}bp.fasta')
            try: 
                SeqIO.write(upstream_records, output_file, "fasta")
                logger.info(f"✅ Saved {len(upstream_records)} sequences to {output_file}")
                total_sequences += len(upstream_records)
            except Exception as e:
                logger.error(f"Error writing output file {output_file}: {e}")
        else:
            logger.warning(f"⚠️ No valid sequences found for {cog_file}")
    
    summary_tsv_path = os.path.join(stats_output, f"cog_summary_upstream_{args.upstream}bp_inside_{args.inside}bp_min_{args.min_length}bp.tsv")
    write_cog_summary_tsv(summary_tsv_path, all_stats, args.min_length, args.upstream, args.inside)
    
    logger.info("="*60)
    logger.info("PROCESSING SUMMARY")
    logger.info("="*60)
    logger.info(f"COG files processed: {len(cog_files)}")
    logger.info(f"Total sequences extracted: {total_sequences}")
    
    total_genes = sum(s['total_genes'] for s in all_stats)
    total_kept = sum(s['sequences_kept'] for s in all_stats)
    
    logger.info(f"Total genes across all COGs: {total_genes}")
    logger.info(f"Total sequences kept: {total_kept}")
    if total_genes > 0:
        overall_retention = (total_kept / total_genes * 100)
        logger.info(f"Overall retention rate: {overall_retention:.2f}%")
    
    logger.info(f"✅ Processing complete. Check {summary_tsv_path} for detailed statistics.")

if __name__ == "__main__":
    main()
