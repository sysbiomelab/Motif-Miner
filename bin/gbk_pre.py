import os
import sys
import argparse
import logging
import pandas as pd
from Bio import SeqIO

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

def extract_cds_info(gbk_file, logger):
    cds_records = []
    try:
        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    qualifiers = feature.qualifiers
                    cds_data = {
                        "File": os.path.basename(gbk_file),
                        "Locus_Tag": qualifiers.get("locus_tag", ["NA"])[0],
                        "Protein_ID": qualifiers.get("protein_id", ["NA"])[0],
                        "Product": qualifiers.get("product", ["NA"])[0],
                        "Start": int(feature.location.start),
                        "End": int(feature.location.end),
                        "Strand": feature.location.strand,
                        "Translation": qualifiers.get("translation", ["NA"])[0],
                        "Gene_Functions": qualifiers.get("gene_functions", ["NA"])[0],
                        "Gene_Kind": qualifiers.get("gene_kind", ["NA"])[0],
                    }
                    cds_records.append(cds_data)
        logger.info(f"Extracted {len(cds_records)} CDS entries from {gbk_file}")
    except Exception as e:
        logger.error(f"Failed to parse {gbk_file}: {e}")
    return cds_records

def process_single_family(family_dir, output_dir, output_name, logger=None):
    all_cds = []
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    if not os.path.exists(family_dir):
        logger.error(f"Family directory does not exist: {family_dir}")
        return
    
    if not os.path.isdir(family_dir):
        logger.error(f"Path is not a directory: {family_dir}")
        return

    logger.info(f"Processing single family")
    
    gbk_files = [f for f in os.listdir(family_dir) if f.endswith((".gbk", ".gbff"))]
    
    if not gbk_files:
        logger.warning(f"No GenBank files found in {family_dir}")
        return
    
    logger.info(f"Found {len(gbk_files)} GenBank files to process")
    
    for file in gbk_files:
        file_path = os.path.join(family_dir, file)
        logger.debug(f"Parsing file: {file_path}")
        cds_info = extract_cds_info(file_path, logger)
        all_cds.extend(cds_info)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")
    
    output_file = os.path.join(output_dir, output_name)
    
    if all_cds:
        df = pd.DataFrame(all_cds)
        column_order = [
            "File", "Locus_Tag", "Start", "End", "Strand",
            "Gene_Functions", "Gene_Kind", "Protein_ID", "Product", "Translation"
        ]
        df = df[column_order]
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"CDS info saved to: {output_file}")
        logger.info(f"Total CDS records processed: {len(all_cds)}")
    else:
        logger.warning("No CDS data found to write.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract CDS info from a single BigSlice family GenBank files.")
    parser.add_argument("family_dir", help="Path to a single family directory (e.g., /path/to/FAM001).")
    parser.add_argument("output_dir", help="Path to output directory")
    parser.add_argument("output_name", nargs='?', default="cds_summary.tsv", help="name of output file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    
    args = parser.parse_args()
    
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = setup_logging(log_level)
    
    process_single_family(args.family_dir, args.output_dir, args.output_name ,logger)
