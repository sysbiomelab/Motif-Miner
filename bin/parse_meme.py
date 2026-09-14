"""
MEME XML Batch Processor with PWM Output

Processes MEME XML files in subdirectories and outputs results in the same directory.
For each subdirectory containing a meme.xml file:
1. Parse the meme.xml file
2. Extract sequence, motif, and binding site data
3. Save results in the same subdirectory
4. Save PWMs in proper MEME format for TOMTOM compatibility

Expected XML structure:
- Root element with training_set, motifs, and scanned_sites_summary
- training_set contains sequence elements with id, name, length, weight
- motifs contains motif elements with PWM data in probabilities/alphabet_matrix
- scanned_sites_summary contains binding site information
"""

import os
import sys
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import logging
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import warnings

warnings.filterwarnings('ignore')


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def find_meme_directories(input_dir: Path) -> List[Path]:
    """Find all subdirectories containing meme.xml files."""
    meme_dirs = []
    
    logging.info(f"Searching for subdirectories with MEME XML files in: {input_dir}")
    
    for subdir in input_dir.iterdir():
        if subdir.is_dir():
            meme_file = subdir / 'meme.xml'
            if meme_file.exists():
                meme_dirs.append(subdir)
                logging.debug(f"Found MEME directory: {subdir} (contains {meme_file.name})")
    
    logging.info(f"Found {len(meme_dirs)} directories with MEME XML files")
    return meme_dirs


def extract_sequences(root: ET.Element, dir_name: str) -> Tuple[List[Dict], Dict]:
    """Extract sequence information from XML."""
    sequences_data = []
    sequence_info = {}
    
    training_set = root.find('training_set')
    if training_set is None:
        logging.warning(f"No training_set found in {dir_name}")
        return sequences_data, sequence_info
    
    for sequence in training_set.findall('sequence'):
        seq_id = sequence.get('id')
        seq_name = sequence.get('name')
        seq_length = int(sequence.get('length')) if sequence.get('length') else None
        seq_weight = float(sequence.get('weight')) if sequence.get('weight') else None
        gene_name = seq_name.replace("_upstream", "") if seq_name else None

        seq_data = {
            'seq_id': seq_id,
            'sequence_name': seq_name,
            'seq_length': seq_length,
            'seq_weight': seq_weight,
            'gene_name': gene_name,
        }
        
        sequences_data.append(seq_data)
        sequence_info[seq_id] = seq_data
    
    scanned_sites_summary = root.find('scanned_sites_summary')
    if scanned_sites_summary is not None:
        for scanned_sites in scanned_sites_summary.findall('scanned_sites'):
            seq_id = scanned_sites.get('sequence_id')
            seq_pvalue = float(scanned_sites.get('pvalue')) if scanned_sites.get('pvalue') else None
            num_sites = int(scanned_sites.get('num_sites')) if scanned_sites.get('num_sites') else 0
            
            if seq_id in sequence_info:
                sequence_info[seq_id]['seq_pvalue'] = seq_pvalue
                sequence_info[seq_id]['num_sites'] = num_sites
                
                for seq_data in sequences_data:
                    if seq_data['seq_id'] == seq_id:
                        seq_data['seq_pvalue'] = seq_pvalue
                        seq_data['num_sites'] = num_sites
                        break
    
    logging.debug(f"Extracted {len(sequences_data)} sequences from {dir_name}")
    return sequences_data, sequence_info


def validate_pwm(pwm_df: pd.DataFrame, motif_id: str, tolerance: float = 0.01) -> bool:
    """Validate that PWM probabilities sum to approximately 1.0 for each position."""
    row_sums = pwm_df.sum(axis=1)
    invalid_rows = row_sums[(row_sums < 1.0 - tolerance) | (row_sums > 1.0 + tolerance)]
    
    if len(invalid_rows) > 0:
        logging.warning(f"Motif {motif_id}: {len(invalid_rows)} positions have probabilities "
                       f"not summing to 1.0 (tolerance={tolerance})")
        return False
    
    return True


def save_pwm_meme_format(motif_id: str, pwm_df: pd.DataFrame, motif_info: Dict, 
                         output_dir: Path, dir_name: str) -> None:
    """Save PWM in MEME format for TOMTOM compatibility with updated naming convention."""
    new_motif_name = f"{dir_name}_{motif_id}"
    pwm_file = output_dir / f"{new_motif_name}.pwm"
    
    validate_pwm(pwm_df, motif_id)
    
    with open(pwm_file, 'w') as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        
        motif_width = motif_info.get('motif_width', len(pwm_df))
        motif_sites = motif_info.get('motif_sites', 10)
        motif_evalue = motif_info.get('motif_evalue', 1e-5)
        
        f.write(f"MOTIF {new_motif_name}\n")
        f.write(f"letter-probability matrix: alength= 4 w= {motif_width} nsites= {motif_sites} E= {motif_evalue:.2e}\n")
        
        # Write PWM matrix (ensure order is A C G T)
        for index, row in pwm_df.iterrows():
            f.write(f"{row['A']:.6f} {row['C']:.6f} {row['G']:.6f} {row['T']:.6f}\n")
        
        f.write("\n")
    
    logging.debug(f"Saved PWM in MEME format: {pwm_file}")


def extract_motifs(root: ET.Element, dir_name: str) -> Tuple[List[Dict], Dict]:
    """Extract motif information from XML including PWM data."""
    motifs_data = []
    motif_pwms = {}
    
    motifs = root.find('motifs')
    if motifs is None:
        logging.warning(f"No motifs found in {dir_name}")
        return motifs_data, motif_pwms
    
    for motif in motifs.findall('motif'):
        motif_id = motif.get('id')
        motif_consensus = motif.get('name')
        motif_width = int(motif.get('width')) if motif.get('width') else None
        motif_sites = int(motif.get('sites')) if motif.get('sites') else None
        motif_ic = float(motif.get('ic')) if motif.get('ic') else None
        motif_re = float(motif.get('re')) if motif.get('re') else None
        motif_llr = float(motif.get('llr')) if motif.get('llr') else None
        motif_evalue = float(motif.get('e_value')) if motif.get('e_value') else None
        
        motif_alt = motif.get('alt')
        bayes_threshold = float(motif.get('bayes_threshold')) if motif.get('bayes_threshold') else None
        elapsed_time = float(motif.get('elapsed_time')) if motif.get('elapsed_time') else None
        
        motif_data = {
            'motif_id': motif_id,
            'motif_seq': motif_consensus,
            'motif_alt': motif_alt,
            'motif_width': motif_width,
            'motif_sites': motif_sites,
            'motif_ic': motif_ic,
            'motif_re': motif_re,
            'motif_llr': motif_llr,
            'motif_evalue': motif_evalue,
            'bayes_threshold': bayes_threshold,
            'elapsed_time': elapsed_time,
            'motif_regex': None
        }
        
        for child in motif:
            if child.tag == 'regular_expression':
                motif_data['motif_regex'] = child.text
            
            elif child.tag == 'probabilities':
                prob_matrix = []
                
                for alphabet_matrix in child.findall('alphabet_matrix'):
                    for alphabet_array in alphabet_matrix.findall('alphabet_array'):
                        row = []
                        for value in alphabet_array.findall('value'):
                            row.append(float(value.text))
                        if row:
                            prob_matrix.append(row)
                
                if prob_matrix:
                    pwm_df = pd.DataFrame(
                        prob_matrix, 
                        columns=['A', 'C', 'G', 'T'],
                        index=[f"Position_{i+1}" for i in range(len(prob_matrix))]
                    )
                    motif_pwms[motif_id] = {
                        'pwm_df': pwm_df,
                        'motif_info': motif_data
                    }
        
        motifs_data.append(motif_data)
    
    logging.debug(f"Extracted {len(motifs_data)} motifs from {dir_name}")
    return motifs_data, motif_pwms


def extract_binding_sites(root: ET.Element, sequence_info: Dict, 
                         motifs_data: List[Dict], dir_name: str) -> List[Dict]:
    """Extract binding site information from XML."""
    combined_data = []
    
    scanned_sites_summary = root.find('scanned_sites_summary')
    if scanned_sites_summary is None:
        logging.warning(f"No scanned_sites_summary found in {dir_name}")
        return combined_data
    
    motifs_df = pd.DataFrame(motifs_data)
    
    for scanned_sites in scanned_sites_summary.findall('scanned_sites'):
        seq_id = scanned_sites.get('sequence_id')
        seq_pvalue = float(scanned_sites.get('pvalue')) if scanned_sites.get('pvalue') else None
        num_sites = int(scanned_sites.get('num_sites')) if scanned_sites.get('num_sites') else 0
        
        seq_info = sequence_info.get(seq_id, {})
        seq_length = seq_info.get('seq_length', None)
        
        for site in scanned_sites.findall('scanned_site'):
            motif_id = site.get('motif_id')
            seq_strand = site.get('strand')
            position = int(site.get('position')) if site.get('position') else None
            site_pvalue = float(site.get('pvalue')) if site.get('pvalue') else None
            
            if position is None:
                continue
            
            motif_matches = motifs_df[motifs_df['motif_id'] == motif_id]
            if len(motif_matches) > 0:
                motif_info = motif_matches.iloc[0]
                motif_width = motif_info.get('motif_width', 0)
                motif_seq = motif_info.get('motif_seq', None)
                motif_evalue = motif_info.get('motif_evalue', None)
            else:
                motif_width = 0
                motif_seq = None
                motif_evalue = None
            
            end_position = position + motif_width - 1 if motif_width else position
            relative_position = round((position / seq_length) * 100, 1) if seq_length else None
            
            # Distance to gene (for upstream regions)
            distance_to_gene_start = seq_length - position if seq_length else None
            distance_to_gene_end = seq_length - end_position if seq_length else None
            closest_distance_to_gene = distance_to_gene_end
            
            combined_row = {            
                'seq_id': seq_id,
                'sequence_name': seq_info.get('sequence_name', seq_id),
                'seq_len': seq_length,
                'seq_strand': seq_strand,
                'seq_weight': seq_info.get('seq_weight', None),
                'seq_pvalue': seq_pvalue,
                'num_sites_in_seq': num_sites,
                'gene_name': seq_info.get('gene_name', None),
                'motif_id': motif_id,
                'motif_seq': motif_seq,
                'motif_width': motif_width,
                'motif_evalue': motif_evalue,
                'start_position': position,
                'end_position': end_position,
                'relative_position_percent': relative_position,
                'distance_to_gene': closest_distance_to_gene,
                'seq_motif_pvalue': site_pvalue,
            }
            
            combined_data.append(combined_row)
    
    logging.debug(f"Extracted {len(combined_data)} binding sites from {dir_name}")
    return combined_data


def save_results_to_directory(sequences_data: List[Dict], motifs_data: List[Dict], 
                             combined_data: List[Dict], motif_pwms: Dict, 
                             output_dir: Path, dir_name: str) -> None:
    """Save all results to the same directory as the meme.xml file."""
    try:
        if sequences_data:
            sequences_df = pd.DataFrame(sequences_data)
            sequences_file = output_dir / "sequences_info.csv"
            sequences_df.to_csv(sequences_file, index=False)
            logging.debug(f"Saved sequences data: {sequences_file}")
        
        if motifs_data:
            motifs_df = pd.DataFrame(motifs_data)
            motifs_file = output_dir / "motifs_info.csv"
            motifs_df.to_csv(motifs_file, index=False)
            logging.debug(f"Saved motifs data: {motifs_file}")
        
        if combined_data:
            combined_df = pd.DataFrame(combined_data)
            combined_file = output_dir / "combined_sequence_motif_binding.csv"
            combined_df.to_csv(combined_file, index=False)
            logging.debug(f"Saved combined data: {combined_file}")
        
        if motif_pwms:
            for motif_id, pwm_data in motif_pwms.items():
                pwm_df = pwm_data['pwm_df']
                motif_info = pwm_data['motif_info']
                
                new_motif_name = f"{dir_name}_{motif_id}"
                
                pwm_csv_file = output_dir / f"{new_motif_name}.csv"
                pwm_df.to_csv(pwm_csv_file, index=True)
                
                save_pwm_meme_format(motif_id, pwm_df, motif_info, output_dir, dir_name)
            
            logging.debug(f"Saved {len(motif_pwms)} PWM files in both CSV and MEME formats")
        
    except Exception as e:
        logging.error(f"Error saving results to {output_dir}: {e}")
        raise


def process_single_directory(meme_dir: Path) -> bool:
    """Process a single directory containing a meme.xml file."""
    try:
        meme_file = None
        for filename in ['meme.xml', 'meme.out']:
            potential_file = meme_dir / filename
            if potential_file.exists():
                meme_file = potential_file
                break
        
        if not meme_file:
            logging.warning(f"No meme.xml file found in {meme_dir}")
            return False
        
        logging.info(f"Processing: {meme_dir.name}/{meme_file.name}")
        
        tree = ET.parse(meme_file)
        root = tree.getroot()
        
        sequences_data, sequence_info = extract_sequences(root, meme_dir.name)
        motifs_data, motif_pwms = extract_motifs(root, meme_dir.name)
        combined_data = extract_binding_sites(root, sequence_info, motifs_data, meme_dir.name)
        
        save_results_to_directory(sequences_data, motifs_data, combined_data, 
                                 motif_pwms, meme_dir, meme_dir.name)
        
        logging.info(f"  ✓ {meme_dir.name}: {len(sequences_data)} sequences, "
                    f"{len(motifs_data)} motifs, {len(combined_data)} binding sites")
        
        return True
        
    except Exception as e:
        logging.error(f"Failed to process {meme_dir.name}: {e}", exc_info=True)
        return False


def generate_overall_summary(successful_dirs: List[str], failed_dirs: List[str]) -> None:
    """Generate an overall summary of processing results."""
    logging.info("\n" + "="*50)
    logging.info("PROCESSING SUMMARY")
    logging.info("="*50)
    
    total_dirs = len(successful_dirs) + len(failed_dirs)
    logging.info(f"Directories processed: {len(successful_dirs)}/{total_dirs}")
    
    if successful_dirs:
        logging.info("\nSuccessfully processed:")
        for dir_name in successful_dirs:
            logging.info(f"  ✓ {dir_name}")
    
    if failed_dirs:
        logging.info("\nFailed to process:")
        for dir_name in failed_dirs:
            logging.info(f"  ✗ {dir_name}")
    
    logging.info(f"\nResults saved in each respective subdirectory:")
    logging.info("  - sequences_info.csv")
    logging.info("  - motifs_info.csv")
    logging.info("  - combined_sequence_motif_binding.csv")
    logging.info("  - [subdirectory]_[motif_id].csv files (CSV format)")
    logging.info("  - [subdirectory]_[motif_id].pwm files (MEME format for TOMTOM)")


def process_directories(input_dir: Path) -> None:
    """Process all directories containing MEME XML files."""
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    meme_dirs = find_meme_directories(input_dir)
    
    if not meme_dirs:
        logging.warning(f"No subdirectories with MEME XML files found in {input_dir}")
        return
    
    successful_dirs = []
    failed_dirs = []
    
    for meme_dir in meme_dirs:
        success = process_single_directory(meme_dir)
        
        if success:
            successful_dirs.append(meme_dir.name)
        else:
            failed_dirs.append(meme_dir.name)
    
    generate_overall_summary(successful_dirs, failed_dirs)


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Process MEME XML files in subdirectories, outputting results in the same directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script processes each subdirectory containing a meme.xml file and saves 
the results in the same subdirectory.

Examples:
  %(prog)s /path/to/meme/analyses
  %(prog)s /path/to/meme/analyses --verbose

Directory structure example:
  input_dir/
  ├── analysis1/
  │   ├── meme.xml
  │   ├── sequences_info.csv          ← generated
  │   ├── motifs_info.csv             ← generated
  │   ├── combined_sequence_motif_binding.csv ← generated
  │   ├── analysis1_motif_1.csv       ← generated (CSV format)
  │   └── analysis1_motif_1.pwm       ← generated (MEME format for TOMTOM)
  ├── analysis2/
  │   ├── meme.xml
  │   └── ... (results saved here)
  └── analysis3/
      ├── meme.xml
      └── ... (results saved here)
        """
    )
    
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing subdirectories with MEME XML files"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        process_directories(args.input_dir)
        
    except Exception as e:
        logging.error(f"Script failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
