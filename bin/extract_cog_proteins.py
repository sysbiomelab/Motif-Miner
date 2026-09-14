#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import argparse
import logging
from datetime import datetime
import glob

def setup_logging(log_level='INFO'):
    """Setup logging configuration with timestamps and proper formatting."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def save_cog_statistics(cog_stats, stats_output_path):
    """Save COG statistics to TSV file at specified path."""
    if os.path.isdir(stats_output_path):
        stats_file = os.path.join(stats_output_path, 'cog_statistics.tsv')
    else:
        stats_file = stats_output_path
    
    os.makedirs(os.path.dirname(stats_file), exist_ok=True)
    
    stats_df = pd.DataFrame(list(cog_stats.items()), columns=['COG_ID', 'Gene_Count'])
    stats_df = stats_df.sort_values('COG_ID')
    
    stats_df.to_csv(stats_file, sep='\t', index=False)
    
    return stats_file, stats_df

def extract_sequences(proteinortho_file, fasta_path, output_dir, stats_output=None, log_level='INFO'):
    """Extract sequences from Proteinortho output."""
    logger = setup_logging(log_level)
    
    logger.info(f"Starting Proteinortho sequence extraction")
    logger.info(f"Input file: {proteinortho_file}")
    logger.info(f"FASTA directory: {fasta_path}")
    logger.info(f"Output directory: {output_dir}")
    
    if stats_output is None:
        stats_output = output_dir
    logger.info(f"Statistics output: {stats_output}")
    
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Created/verified output directory: {output_dir}")

    logger.info("Loading Proteinortho TSV file...")
    df = pd.read_csv(proteinortho_file, sep="\t")
    logger.info(f"Loaded {len(df)} COGs from proteinortho file")

    genome_cols = df.columns[3:]
    logger.info(f"Found {len(genome_cols)} genomes: {', '.join(genome_cols)}")
    
    genome_proteins = {}
    missing_files = []

    logger.info("Loading protein sequences from FASTA files...")
    for genome_file in genome_cols:
        full_path = os.path.join(fasta_path, genome_file)
        if not os.path.exists(full_path):
            logger.warning(f"Missing file: {full_path}")
            missing_files.append(genome_file)
            continue
            
        logger.info(f"Loading {genome_file}...")
        try:
            protein_dict = SeqIO.to_dict(SeqIO.parse(full_path, "fasta"))
            genome_proteins[genome_file] = protein_dict
            logger.debug(f"Loaded {len(protein_dict)} proteins from {genome_file}")
        except Exception as e:
            logger.error(f"Error loading {genome_file}: {str(e)}")
            missing_files.append(genome_file)

    if missing_files:
        logger.warning(f"Missing {len(missing_files)} genome files: {', '.join(missing_files)}")

    logger.info("Processing COGs and extracting sequences...")
    cog_statistics = {}
    successful_cogs = 0
    empty_cogs = 0
    total_genes_extracted = 0

    for i, row in df.iterrows():
        cog_id = f"COG_{i + 1:05d}"
        fasta_records = []
        missing_proteins = []
        
        logger.debug(f"Processing {cog_id}...")
        
        for genome_file in genome_cols:
            if genome_file not in genome_proteins:
                continue
                
            ids = row[genome_file]
            if isinstance(ids, str) and ids not in ['*', '_GENOME_HIT_']:
                for protein_id in ids.split(','):
                    protein_id = protein_id.strip()
                    logger.debug(f"Looking for protein: {protein_id} in {genome_file}")
                    
                    if protein_id in genome_proteins[genome_file]:
                        record = genome_proteins[genome_file][protein_id]
                        new_id = f"{genome_file}|{record.id}"
                        new_description = f"{genome_file} {record.description}"
                        fasta_records.append(SeqRecord(record.seq, id=new_id, description=new_description))
                        logger.debug(f"Added protein {protein_id} from {genome_file}")
                    else:
                        logger.warning(f"Missing protein: {protein_id} in {genome_file}")
                        missing_proteins.append(f"{protein_id} ({genome_file})")
        
        gene_count = len(fasta_records)
        cog_statistics[cog_id] = gene_count
        
        if fasta_records:
            out_file = os.path.join(output_dir, f"{cog_id}.fasta")
            write(fasta_records, out_file, "fasta")
            logger.info(f"✅ Wrote {gene_count} genes to {cog_id}.fasta")
            successful_cogs += 1
            total_genes_extracted += gene_count
        else:
            logger.warning(f"❌ No sequences found for {cog_id}")
            empty_cogs += 1

        if missing_proteins:
            logger.debug(f"Missing proteins in {cog_id}: {', '.join(missing_proteins)}")

    logger.info("Saving COG statistics...")
    stats_file, stats_df = save_cog_statistics(cog_statistics, stats_output)
    logger.info(f"Saved COG statistics to: {stats_file}")
    
    logger.info("=" * 50)
    logger.info("EXTRACTION SUMMARY")
    logger.info("=" * 50)
    logger.info(f"Total COGs processed: {len(df)}")
    logger.info(f"Successful COGs: {successful_cogs}")
    logger.info(f"Empty COGs: {empty_cogs}")
    logger.info(f"Total genes extracted: {total_genes_extracted}")
    logger.info(f"Average genes per COG: {total_genes_extracted/successful_cogs:.2f}" if successful_cogs > 0 else "Average genes per COG: 0")
    logger.info(f"Missing genome files: {len(missing_files)}")
    
    if not stats_df.empty:
        logger.info(f"COG size statistics:")
        logger.info(f"  - Min genes per COG: {stats_df['Gene_Count'].min()}")
        logger.info(f"  - Max genes per COG: {stats_df['Gene_Count'].max()}")
        logger.info(f"  - Mean genes per COG: {stats_df['Gene_Count'].mean():.2f}")
        logger.info(f"  - Median genes per COG: {stats_df['Gene_Count'].median():.2f}")
    
    logger.info("Extraction completed successfully!")
    return 0

def count_sequences_in_fasta(fasta_file):
    """Count sequences in a FASTA file."""
    try:
        count = 0
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                count += 1
        return count
    except Exception as e:
        logging.error(f"Error reading {fasta_file}: {str(e)}")
        return 0

def get_genome_breakdown(fasta_file):
    """Get breakdown of sequences by genome for a COG file."""
    genome_counts = {}
    try:
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Extract genome name from sequence ID (format: genome|protein_id)
                genome = record.id.split('|')[0] if '|' in record.id else 'unknown'
                genome_counts[genome] = genome_counts.get(genome, 0) + 1
    except Exception as e:
        logging.error(f"Error analyzing {fasta_file}: {str(e)}")
    
    return genome_counts

def count_sequences(cog_directory, output_file='cog_counts.tsv', pattern='*.fasta', detailed=False, log_level='INFO'):
    """Count sequences in COG FASTA files."""
    logger = setup_logging(log_level)
    
    logger.info(f"Starting COG sequence counting")
    logger.info(f"COG directory: {cog_directory}")
    logger.info(f"File pattern: {pattern}")
    logger.info(f"Output file: {output_file}")
    logger.info(f"Detailed analysis: {detailed}")
    
    if not os.path.exists(cog_directory):
        logger.error(f"Directory does not exist: {cog_directory}")
        return 1
    
    search_pattern = os.path.join(cog_directory, pattern)
    fasta_files = glob.glob(search_pattern)
    
    if not fasta_files:
        logger.error(f"No FASTA files found matching pattern: {search_pattern}")
        return 1
    
    logger.info(f"Found {len(fasta_files)} FASTA files to process")
    
    cog_data = []
    total_sequences = 0
    processed_files = 0
    empty_files = 0
    
    for fasta_file in sorted(fasta_files):
        filename = os.path.basename(fasta_file)
        cog_id = os.path.splitext(filename)[0]
        
        logger.debug(f"Processing {filename}...")
        
        seq_count = count_sequences_in_fasta(fasta_file)
        
        if seq_count == 0:
            logger.warning(f"Empty file: {filename}")
            empty_files += 1
        else:
            logger.debug(f"{filename}: {seq_count} sequences")
            processed_files += 1
        
        total_sequences += seq_count
        
        row_data = {
            'COG_ID': cog_id,
            'Filename': filename,
            'Sequence_Count': seq_count,
            'File_Size_KB': round(os.path.getsize(fasta_file) / 1024, 2)
        }
        
        if detailed and seq_count > 0:
            genome_counts = get_genome_breakdown(fasta_file)
            row_data['Genome_Count'] = len(genome_counts)
            row_data['Genomes'] = ';'.join([f"{genome}:{count}" for genome, count in sorted(genome_counts.items())])
            logger.debug(f"{filename}: {len(genome_counts)} genomes represented")
        
        cog_data.append(row_data)
    
    df = pd.DataFrame(cog_data)
    
    df = df.sort_values('COG_ID')
    
    output_path = os.path.join(cog_directory, output_file) if not os.path.dirname(output_file) else output_file
    df.to_csv(output_path, sep='\t', index=False)
    
    logger.info(f"Saved results to: {output_path}")
    
    logger.info("=" * 60)
    logger.info("COG COUNTING SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total files processed: {len(fasta_files)}")
    logger.info(f"Files with sequences: {processed_files}")
    logger.info(f"Empty files: {empty_files}")
    logger.info(f"Total sequences counted: {total_sequences}")
    
    if processed_files > 0:
        logger.info(f"Average sequences per COG: {total_sequences/processed_files:.2f}")
        logger.info(f"Min sequences per COG: {df[df['Sequence_Count'] > 0]['Sequence_Count'].min()}")
        logger.info(f"Max sequences per COG: {df['Sequence_Count'].max()}")
        logger.info(f"Median sequences per COG: {df[df['Sequence_Count'] > 0]['Sequence_Count'].median():.2f}")
    
    if not df.empty:
        logger.info("\nTop 10 largest COGs:")
        top_10 = df.nlargest(10, 'Sequence_Count')[['COG_ID', 'Sequence_Count']]
        for _, row in top_10.iterrows():
            logger.info(f"  {row['COG_ID']}: {row['Sequence_Count']} sequences")
    
    if processed_files > 0:
        logger.info("\nSequence count distribution:")
        bins = [0, 1, 5, 10, 20, 50, 100, float('inf')]
        labels = ['0', '1', '2-5', '6-10', '11-20', '21-50', '51-100', '>100']
        
        for i in range(len(bins)-1):
            count = len(df[(df['Sequence_Count'] > bins[i]) & (df['Sequence_Count'] <= bins[i+1])])
            if count > 0:
                logger.info(f"  {labels[i+1]} sequences: {count} COGs")
    
    logger.info("Counting completed successfully!")
    return 0

def main():
    parser = argparse.ArgumentParser(description='Extract sequences from Proteinortho output and count COG sequences.')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    extract_parser = subparsers.add_parser('extract', help='Extract sequences from Proteinortho output')
    extract_parser.add_argument('proteinortho_file', help='Path to Proteinortho TSV file')
    extract_parser.add_argument('fasta_path', help='Path to directory with .faa protein files')
    extract_parser.add_argument('--output_dir', default='output_fasta', help='Directory to save output FASTA files')
    extract_parser.add_argument('--stats_output', help='Path to save COG statistics file (default: same as output_dir)')
    extract_parser.add_argument('--log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                               help='Set logging level (default: INFO)')
    
    count_parser = subparsers.add_parser('count', help='Count sequences in COG FASTA files')
    count_parser.add_argument('cog_directory', help='Directory containing COG FASTA files')
    count_parser.add_argument('--output_file', default='cog_counts.tsv', help='Output TSV file name (default: cog_counts.tsv)')
    count_parser.add_argument('--pattern', default='*.fasta', help='File pattern to match (default: *.fasta)')
    count_parser.add_argument('--detailed', action='store_true', help='Include detailed genome breakdown')
    count_parser.add_argument('--log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                             help='Set logging level (default: INFO)')
    
    both_parser = subparsers.add_parser('both', help='Extract sequences then count them')
    both_parser.add_argument('proteinortho_file', help='Path to Proteinortho TSV file')
    both_parser.add_argument('fasta_path', help='Path to directory with .faa protein files')
    both_parser.add_argument('--output_dir', default='output_fasta', help='Directory to save output FASTA files')
    both_parser.add_argument('--stats_output', help='Path to save COG statistics file (default: same as output_dir)')
    both_parser.add_argument('--pattern', default='*.fasta', help='File pattern to match for counting (default: *.fasta)')
    both_parser.add_argument('--detailed', action='store_true', help='Include detailed genome breakdown in counting')
    both_parser.add_argument('--log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                            help='Set logging level (default: INFO)')

    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    if args.command == 'extract':
        return extract_sequences(args.proteinortho_file, args.fasta_path, args.output_dir, args.stats_output, args.log_level)
    
    elif args.command == 'count':
        return count_sequences(args.cog_directory, args.output_file, args.pattern, args.detailed, args.log_level)
    
    elif args.command == 'both':
        logger = setup_logging(args.log_level)
        logger.info("Running both extraction and counting...")
        
        extract_result = extract_sequences(args.proteinortho_file, args.fasta_path, args.output_dir, args.stats_output, args.log_level)
        if extract_result != 0:
            logger.error("Extraction failed, skipping counting step")
            return extract_result
        
        logger.info("Extraction completed, starting counting...")
        
        count_result = count_sequences(args.output_dir, 'cog_counts.tsv', args.pattern, args.detailed, args.log_level)
        
        return count_result

if __name__ == "__main__":
    exit(main())
