import logging
import sys
import os
import re
import argparse
import pandas as pd
from typing import Optional, Tuple, List, Dict

def setup_logging(level=logging.INFO) -> logging.Logger:
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

logger = setup_logging()

def parse_header(header: str) -> Optional[Tuple[str, str, int, int, int, int]]:
    """Parse FASTA header to extract gene information"""
    try:
        parts = header.split()
        if len(parts) < 3:
            logger.warning(f"Header has insufficient parts (expected ≥3, got {len(parts)}): {header[:100]}...")
            return None
        
        file_id = parts[1]
        info = parts[2]
        parts_info = info.split('|')
        
        if len(parts_info) < 3:
            logger.warning(f"Info section has insufficient parts (expected ≥3, got {len(parts_info)}): {info}")
            return None
        
        locus_tag = parts_info[0]
        
        try:
            start_str, end_str = parts_info[1].split('-')
            start, end = int(start_str), int(end_str)
        except ValueError:
            logger.error(f"Invalid coordinate format in '{parts_info[1]}' from header: {header[:100]}...")
            return None
            
        try:
            strand = int(parts_info[2])
        except ValueError:
            logger.error(f"Invalid strand value '{parts_info[2]}' from header: {header[:100]}...")
            return None

        length_match = re.search(r'length=(\d+)bp', header)
        length = int(length_match.group(1)) if length_match else None
        
        if length is None:
            logger.debug(f"No length information found in header: {header[:100]}...")

        return file_id, locus_tag, start, end, strand, length
    
    except (ValueError, IndexError) as e:
        logger.error(f"Unexpected error parsing header '{header[:100]}...': {type(e).__name__}: {e}")
        return None

def debug_fasta_headers(file_path: str, num_headers: int = 3) -> None:
    """Debug function to show sample FASTA headers"""
    logger.info(f"Debugging FASTA headers from {os.path.basename(file_path)}")
    try:
        with open(file_path, 'r') as f:
            count = 0
            for line in f:
                if line.startswith('>') and count < num_headers:
                    logger.debug(f"Header {count + 1}: {line.strip()}")
                    parsed = parse_header(line[1:].strip())
                    if parsed:
                        logger.debug(f"  Successfully parsed locus_tag: {parsed[1]}")
                    else:
                        logger.warning(f"  Failed to parse header {count + 1}")
                    count += 1
                    if count >= num_headers:
                        break
        logger.info(f"Completed debugging {count} headers from {os.path.basename(file_path)}")
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {type(e).__name__}: {e}")


def parse_fasta_file(file_path: str, debug: bool = False) -> List[Dict]:
    """Parse a single FASTA file and return gene information."""
    genes = []
    logger.info(f"Starting to parse FASTA file: {os.path.basename(file_path)}")
    
    try:
        if debug:
            debug_fasta_headers(file_path)

        parsed_count = 0
        count = 0
        
        with open(file_path, "r") as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    count += 1
                    parsed_header = parse_header(line[1:])

                    if parsed_header:
                        file_id, locus_tag, start, end, strand, length = parsed_header
                        genes.append({
                            'file_id': file_id,
                            'locus_tag': locus_tag,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'length': length
                        })
                        parsed_count += 1
                    else:
                        logger.warning(f"Failed to parse header on line {line_num}")

        for gene in genes:
            gene['total'] = count

        logger.info(
            f"Completed parsing {os.path.basename(file_path)}: "
            f"{parsed_count}/{count} headers successfully parsed "
            f"({count - parsed_count} failed)"
        )

    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
    except PermissionError:
        logger.error(f"Permission denied reading file: {file_path}")
    except Exception as e:
        logger.error(f"Unexpected error reading file {file_path}: {type(e).__name__}: {e}")
    
    return genes

def analyze_gene_types(genes: List[Dict], df: pd.DataFrame, debug=False) -> Dict[str, int]:
    """Analyze gene types based on locus tags"""
    gene_counts = {
        'regulatory': 0,
        'biosynthetic': 0,
        'biosynthetic_additional': 0,
        'total': len(genes)
    }
    
    logger.info(f"Starting gene type analysis for {len(genes)} genes")
    
    if debug:
        logger.debug(f"DataFrame shape: {df.shape}")
        logger.debug(f"DataFrame columns: {list(df.columns)}")
        if 'Locus_Tag' in df.columns:
            logger.debug(f"Sample DataFrame locus tags: {df['Locus_Tag'].head().tolist()}")
        else:
            logger.warning("No 'Locus_Tag' column found in DataFrame")
        logger.debug(f"Sample FASTA locus tags: {[g['locus_tag'] for g in genes[:5]]}")
    
    locus_col = "Locus_Tag"
    gene_kind_col = "Gene_Kind"
    
    if locus_col not in df.columns:
        logger.error(f"Required column '{locus_col}' not found in DataFrame. Available columns: {list(df.columns)}")
        return gene_counts
        
    if gene_kind_col not in df.columns:
        logger.error(f"Required column '{gene_kind_col}' not found in DataFrame. Available columns: {list(df.columns)}")
        return gene_counts
    
    logger.debug(f"Using locus column: '{locus_col}', gene kind column: '{gene_kind_col}'")
    
    matches_found = 0
    for i, gene in enumerate(genes):
        locus_tag = gene['locus_tag']
        
        matching_row = df[df[locus_col] == locus_tag]
        
        if not matching_row.empty:
            gene_kind = matching_row[gene_kind_col].iloc[0]
            matches_found += 1
            
            if debug and matches_found <= 5:
                logger.debug(f"Match {matches_found}: {locus_tag} -> {gene_kind}")
            
            gene_kind_lower = str(gene_kind).lower()
            if 'biosynthetic-additional' in gene_kind_lower:
                gene_counts['biosynthetic_additional'] += 1
            elif 'biosynthetic' in gene_kind_lower:
                gene_counts['biosynthetic'] += 1
            elif 'regulatory' in gene_kind_lower:
                gene_counts['regulatory'] += 1
        else:
            if debug and i < 5:
                logger.debug(f"No match found for locus_tag: {locus_tag}")
    
    match_percentage = (matches_found / len(genes) * 100) if genes else 0
    logger.info(f"Gene type analysis complete: {matches_found}/{len(genes)} genes matched ({match_percentage:.1f}%)")
    
    if matches_found == 0:
        logger.warning("No gene matches found. Check that locus tags in FASTA match those in reference file.")
    elif match_percentage < 50:
        logger.warning(f"Low match percentage ({match_percentage:.1f}%). Verify data compatibility.")
    
    return gene_counts

def calculate_percentages(gene_counts: Dict[str, int]) -> Dict[str, float]:
    """Calculate percentages for each gene type"""
    total = gene_counts['total']
    percentages = {}
    
    if total == 0:
        logger.warning("Cannot calculate percentages: total gene count is 0")
        for gene_type in gene_counts:
            if gene_type != 'total':
                percentages[gene_type] = 0.0
        return percentages
    
    for gene_type, count in gene_counts.items():
        if gene_type != 'total':
            percentages[gene_type] = (count / total * 100)
    
    logger.debug(f"Calculated percentages: {percentages}")
    return percentages

def main():
    parser = argparse.ArgumentParser(description='Analyze FASTA files for a single gene family')
    parser.add_argument('--family_dir', '-f', required=True, 
                        help='Path to the family directory containing FASTA files')
    parser.add_argument('--familyname', '-n', required=True,
                        help='Name of the family to use in the output (first column)')
    parser.add_argument('--reference_file', '-r', required=True,
                        help='Path to the reference TSV file (bigslice_cds_summary.tsv)')
    parser.add_argument('--output', '-o', default='gene_analysis_results.tsv',
                        help='Output file path. If just filename, saves in current directory. If full path, saves there. (default: gene_analysis_results.tsv)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging (DEBUG level)')
    
    args = parser.parse_args()
    
    output_path = args.output
    if not os.path.isabs(output_path):
        output_path = os.path.abspath(output_path)
    
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        logger.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")
    
    logger.info("Starting single family gene analysis pipeline")
    logger.info(f"Family name: {args.familyname}")
    logger.info(f"Family directory: {args.family_dir}")
    logger.info(f"Reference file: {args.reference_file}")
    logger.info(f"Output file: {output_path}")
    
    if not os.path.exists(args.family_dir):
        logger.error(f"Family directory does not exist: {args.family_dir}")
        return 1
    
    if not os.path.isdir(args.family_dir):
        logger.error(f"Family path is not a directory: {args.family_dir}")
        return 1
    
    if not os.path.exists(args.reference_file):
        logger.error(f"Reference file does not exist: {args.reference_file}")
        return 1
    
    try:
        logger.info("Loading reference file...")
        df = pd.read_csv(args.reference_file, sep='\t')
        logger.info(f"Successfully loaded reference file: {len(df)} rows, {len(df.columns)} columns")
        logger.debug(f"Reference file columns: {list(df.columns)}")
    except FileNotFoundError:
        logger.error(f"Reference file not found: {args.reference_file}")
        return 1
    except pd.errors.EmptyDataError:
        logger.error(f"Reference file is empty: {args.reference_file}")
        return 1
    except Exception as e:
        logger.error(f"Error loading reference file: {type(e).__name__}: {e}")
        return 1
    
    column_order = [
        'family', 'cog', 'total', 'regulatory', 'biosynthetic', 'biosynthetic_additional'
    ]
    df_info = pd.DataFrame(columns=column_order)
    
    family_name = args.familyname
    logger.info(f"Processing family: {family_name}")
    
    fasta_files = [f for f in os.listdir(args.family_dir) if f.endswith('.fasta')]
    if not fasta_files:
        logger.error(f"No FASTA files found in {args.family_dir}")
        return 1
        
    logger.info(f"Found {len(fasta_files)} FASTA files in {args.family_dir}")
    
    for j, file_name in enumerate(fasta_files, 1):
        file_path = os.path.join(args.family_dir, file_name)
        logger.info(f"Processing file {j}/{len(fasta_files)}: {file_name}")
        
        debug_mode = j == 1 and args.verbose
        genes = parse_fasta_file(file_path, debug=debug_mode)
        
        if not genes:
            logger.warning(f"No genes extracted from {file_name}")
            continue
        
        gene_counts = analyze_gene_types(genes, df, debug=debug_mode)
        
        percentages = calculate_percentages(gene_counts)
        
        result_row = {
            'family': family_name,
            'cog': file_name.replace('.fasta', ''),
            'total': gene_counts['total'],
            'regulatory': f"{percentages['regulatory']:.1f}",
            'biosynthetic': f"{percentages['biosynthetic']:.1f}",
            'biosynthetic_additional': f"{percentages['biosynthetic_additional']:.1f}"
        }
        
        df_info = pd.concat([df_info, pd.DataFrame([result_row])], ignore_index=True)
        
        logger.info(f"Completed {file_name}: {gene_counts['total']} total genes, "
                   f"{gene_counts['regulatory']} regulatory, "
                   f"{gene_counts['biosynthetic']} biosynthetic, "
                   f"{gene_counts['biosynthetic_additional']} biosynthetic-additional")
    
    if not df_info.empty:
        try:
            df_info.to_csv(output_path, sep='\t', index=False)
            logger.info(f"Saved {len(df_info)} results to: {output_path}")
        except Exception as e:
            logger.error(f"Failed to save results: {type(e).__name__}: {e}")
            return 1
        
        logger.info("=" * 50)
        logger.info("GENE ANALYSIS SUMMARY")
        logger.info("=" * 50)
        logger.info(f"Family: {family_name}")
        logger.info(f"Total COG files processed: {len(df_info)}")
        logger.info(f"Results saved to: {output_path}")
        
        total_genes = df_info['total'].astype(int).sum()
        logger.info(f"Total genes analyzed: {total_genes:,}")
        
        avg_regulatory = df_info['regulatory'].astype(float).mean()
        avg_biosynthetic = df_info['biosynthetic'].astype(float).mean()
        avg_biosynthetic_additional = df_info['biosynthetic_additional'].astype(float).mean()
        
        logger.info(f"\nFamily-wide averages:")
        logger.info(f"  Regulatory: {avg_regulatory:.1f}%")
        logger.info(f"  Biosynthetic: {avg_biosynthetic:.1f}%")
        logger.info(f"  Biosynthetic-additional: {avg_biosynthetic_additional:.1f}%")
        
        logger.info(f"\nAll {len(df_info)} results:")
        for _, row in df_info.iterrows():
            logger.info(f"  {row['cog']}: {row['total']} genes "
                       f"(R:{row['regulatory']}%, B:{row['biosynthetic']}%, BA:{row['biosynthetic_additional']}%)")
        
        logger.info("Analysis completed successfully!")
        return 0
    else:
        logger.error("No results generated. Please check your input files and reference data.")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
