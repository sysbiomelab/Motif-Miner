"""
GenBank to Protein FASTA Converter

Converts GenBank (.gbk) files to protein FASTA (.faa) files by extracting
CDS features.
"""

import pandas as pd
import os
import argparse
import shutil
import logging
import re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path
from typing import Optional, List, Tuple
import sys


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def extract_locus_tag_location_strand(feature) -> str:
    """Extract locus_tag, location (start-end), and strand from a CDS feature."""
    qualifiers = feature.qualifiers

    if 'locus_tag' in qualifiers:
        locus_tag = qualifiers['locus_tag'][0]
    else:
        start = int(feature.location.start)
        end = int(feature.location.end)
        locus_tag = f"unknown_locus_{start}_{end}"

    start = int(feature.location.start)
    end = int(feature.location.end)

    strand_num = feature.location.strand
    if strand_num == 1:
        strand = "1"
    elif strand_num == -1:
        strand = "-1"
    else:
        strand = "."

    return f"{locus_tag}|{start}-{end}|{strand}"


def convert_gbk_to_faa(gbk_file: Path, output_file: Path, 
                      include_partial: bool = False) -> Tuple[int, int, int]:
    """Convert a single GenBank file to protein FASTA format."""
    # Define valid amino acids pattern (inverse match - find invalid chars)
    valid_aa_pattern = re.compile(r'[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]')
    
    total_cds = 0
    extracted_proteins = 0
    invalid_chars_removed = 0
    
    try:
        with open(output_file, "w") as faa_out:
            for record in SeqIO.parse(gbk_file, "genbank"):
                logging.debug(f"Processing record: {record.id}")
                
                for feature in record.features:
                    if feature.type == "CDS":
                        total_cds += 1
                        
                        if not include_partial:
                            start_pos = feature.location.start
                            end_pos = feature.location.end
                            
                            # Skip if positions are fuzzy (BeforePosition, AfterPosition, etc.)
                            if (hasattr(start_pos, '__class__') and 
                                start_pos.__class__.__name__ in ['BeforePosition', 'AfterPosition', 'WithinPosition'] or
                                hasattr(end_pos, '__class__') and 
                                end_pos.__class__.__name__ in ['BeforePosition', 'AfterPosition', 'WithinPosition']):
                                continue
                        
                        if "translation" in feature.qualifiers:
                            protein_seq = feature.qualifiers["translation"][0]
                            
                            invalid_matches = valid_aa_pattern.findall(protein_seq)
                            if invalid_matches:
                                num_invalid = len(invalid_matches)
                                invalid_chars_removed += num_invalid
                                unique_invalid = set(invalid_matches)
                                header = extract_locus_tag_location_strand(feature)
                                logging.debug(f"Found {num_invalid} invalid character(s) in {header}: {unique_invalid}")
                            
                            clean_protein_seq = valid_aa_pattern.sub('', protein_seq)
                            
                            header = extract_locus_tag_location_strand(feature)
                            
                            faa_out.write(f">{header}\n{clean_protein_seq}\n")
                            extracted_proteins += 1
                        else:
                            logging.warning(f"CDS feature without translation in {gbk_file.name}: {extract_locus_tag_location_strand(feature)}")
    
    except Exception as e:
        logging.error(f"Error processing {gbk_file}: {e}")
        raise
    
    return total_cds, extracted_proteins, invalid_chars_removed


def process_directory(input_dir: Path, output_dir: Optional[Path] = None,
                     include_partial: bool = False, overwrite: bool = False) -> None:
    """Process all GenBank files in a directory."""
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    if output_dir is None:
        output_dir = input_dir / "faa_files"
    
    output_dir.mkdir(exist_ok=True)
    logging.info(f"Output directory: {output_dir}")
    
    gbk_files = list(input_dir.glob("*.gbk")) + list(input_dir.glob("*.gb"))
    
    if not gbk_files:
        logging.warning(f"No GenBank files found in {input_dir}")
        return
    
    logging.info(f"Found {len(gbk_files)} GenBank files to process")
    
    total_files_processed = 0
    total_proteins_extracted = 0
    total_invalid_chars = 0
    files_with_invalid_chars = []
    
    for gbk_file in gbk_files:
        output_file = output_dir / (gbk_file.stem + ".faa")
        
        if output_file.exists() and not overwrite:
            logging.info(f"Skipping {gbk_file.name} (output exists, use --overwrite to replace)")
            continue
        
        logging.info(f"Processing: {gbk_file.name}")
        
        try:
            total_cds, extracted_proteins, invalid_chars = convert_gbk_to_faa(
                gbk_file, output_file, include_partial
            )
            
            if invalid_chars > 0:
                logging.info(f"  {gbk_file.name}: {extracted_proteins}/{total_cds} proteins extracted, {invalid_chars} invalid character(s) removed")
                files_with_invalid_chars.append((gbk_file.name, invalid_chars))
            else:
                logging.info(f"  {gbk_file.name}: {extracted_proteins}/{total_cds} proteins extracted")
            
            total_files_processed += 1
            total_proteins_extracted += extracted_proteins
            total_invalid_chars += invalid_chars
            
        except Exception as e:
            logging.error(f"Failed to process {gbk_file.name}: {e}")
            continue
    
    logging.info(f"\nSummary:")
    logging.info(f"  Files processed: {total_files_processed}/{len(gbk_files)}")
    logging.info(f"  Total proteins extracted: {total_proteins_extracted}")
    
    if total_invalid_chars > 0:
        logging.info(f"  Total invalid characters removed: {total_invalid_chars}")
        logging.info(f"  Files with invalid characters: {len(files_with_invalid_chars)}")
        if files_with_invalid_chars:
            logging.info(f"  Files cleaned:")
            for filename, count in files_with_invalid_chars:
                logging.info(f"    - {filename}: {count} character(s)")


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Convert GenBank files to protein FASTA format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s /path/to/gbk/files
  %(prog)s /path/to/gbk/files -o /path/to/output
  %(prog)s /path/to/gbk/files --include-partial --overwrite
  %(prog)s /path/to/gbk/files -v  # verbose mode to see which sequences had invalid chars
        """
    )
    
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing GenBank (.gbk) files"
    )
    
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        help="Output directory (default: input_dir/faa_files)"
    )
    
    parser.add_argument(
        "--include-partial",
        action="store_true",
        help="Include partial CDS features"
    )
    
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging (shows which specific sequences had invalid characters)"
    )
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        process_directory(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            include_partial=args.include_partial,
            overwrite=args.overwrite
        )
        
    except Exception as e:
        logging.error(f"Script failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
