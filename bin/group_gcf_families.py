import pandas as pd
import os
import argparse
import shutil
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('bgc_grouping.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def count_cds_in_genbank(file_path):
    """Count the number of CDS features in a GenBank file."""
    try:
        record = SeqIO.read(file_path, "genbank")
        cds_count = sum(1 for feature in record.features if feature.type == "CDS")
        return cds_count
    except Exception as e:
        logger.warning(f"Error reading {file_path}: {e}")
        return 0

parser = argparse.ArgumentParser(description="Group BGCs by GCF family.")
parser.add_argument("tsv_file", help="Path to input TSV file")
parser.add_argument("file_dir", help="Path to directory containing GenBank files")
parser.add_argument("output_dir", help="Path to output directory")

args = parser.parse_args()
tsv_file = args.tsv_file
file_dir = args.file_dir
output_dir = args.output_dir

logger.info("Starting BGC grouping script")
logger.info(f"Input TSV: {tsv_file}")
logger.info(f"File directory: {file_dir}")
logger.info(f"Output directory: {output_dir}")

grouped_gcf_dir = f'{output_dir}/'
os.makedirs(grouped_gcf_dir, exist_ok=True)
logger.info(f"Created output directory: {grouped_gcf_dir}")

try:
    df = pd.read_csv(tsv_file, sep='\t')
    logger.info(f"Successfully loaded {len(df)} records from TSV file")
    logger.info(f"Found {df['gcf_id'].nunique()} unique GCF families")
except Exception as e:
    logger.error(f"Error reading TSV file: {e}")
    exit(1)

summary_data = []
total_files_processed = 0
total_files_missing = 0

for fam in df['gcf_id'].unique():
    fam_name = f'FAM{str(fam)}'
    fam_dir = os.path.join(grouped_gcf_dir, fam_name)
    
    logger.info(f"Processing family {fam_name}")
    
    df_tmp = df[df["gcf_id"] == fam]
    file_count = len(df_tmp)
    
    logger.info(f"  Family {fam_name} contains {file_count} files")
    
    family_summary_data = []
    files_copied = 0
    files_missing = 0
    
    files_exist = False
    for file in df_tmp['new_name']:
        src = os.path.join(file_dir, f"{file}")
        if os.path.exists(src):
            files_exist = True
            break
    
    if not files_exist:
        logger.warning(f"  Family {fam_name} skipped - no files found")
        total_files_missing += file_count
        continue
    
    os.makedirs(fam_dir, exist_ok=True)
    
    for file in df_tmp['new_name']:
        src = os.path.join(file_dir, f"{file}")
        dst = os.path.join(fam_dir, f"{file}")
        
        try:
            shutil.copyfile(src, dst)
            files_copied += 1
            total_files_processed += 1
            
            cds_count = count_cds_in_genbank(src)
            
            family_summary_data.append({
                'bgc_name': file,
                'cds_count': cds_count
            })
            
            logger.debug(f"    Copied {file} (CDS count: {cds_count})")
            
        except FileNotFoundError:
            logger.warning(f"    File not found: {src}")
            files_missing += 1
            total_files_missing += 1
            continue
        except Exception as e:
            logger.error(f"    Error processing {file}: {e}")
            files_missing += 1
            total_files_missing += 1
            continue
    
    logger.info(f"  Family {fam_name} summary: {files_copied} files copied, {files_missing} files missing")
    
    if family_summary_data:
        family_df = pd.DataFrame(family_summary_data)
        family_df = family_df.sort_values('cds_count', ascending=False)
        
        total_cds = family_df['cds_count'].sum()
        avg_cds = family_df['cds_count'].mean()
        median_cds = family_df['cds_count'].median()
        
        family_summary_path = os.path.join(fam_dir, f"{fam_name}_summary.tsv")
        family_df.to_csv(family_summary_path, sep="\t", index=False)
        
        logger.info(f"  Family {fam_name} CDS statistics: Total={total_cds}, Average={avg_cds:.1f}, Median={median_cds}")
        logger.info(f"  Family summary saved to: {family_summary_path}")
        
        summary_data.append({
            'gcf_id': fam,
            'family_name': fam_name,
            'total_files': file_count,
            'files_copied': files_copied,
            'files_missing': files_missing,
            'total_cds': sum(item['cds_count'] for item in family_summary_data),
            'avg_cds_per_file': sum(item['cds_count'] for item in family_summary_data) / files_copied if files_copied > 0 else 0
        })

summary_df = pd.DataFrame(summary_data)

summary_df = summary_df.sort_values('files_copied', ascending=False)
summary_path = os.path.join(output_dir, "gcf_summary.tsv")
summary_df.to_csv(summary_path, sep="\t", index=False)

logger.info("="*50)
logger.info("FINAL SUMMARY")
logger.info("="*50)
logger.info(f"Total families processed: {len(df['gcf_id'].unique())}")
logger.info(f"Total files processed successfully: {total_files_processed}")
logger.info(f"Total files missing: {total_files_missing}")
logger.info(f"Overall summary saved to: {summary_path}")

logger.info("\nTop 5 families by file count:")
for idx, row in summary_df.head().iterrows():
    logger.info(f"  {row['family_name']}: {row['files_copied']} files, {row['total_cds']} total CDS")

logger.info("BGC grouping script completed successfully!")
