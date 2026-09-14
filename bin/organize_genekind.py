"""
MEME Analysis Organizer - Enhanced Version
Organizes MEME analysis directories based on percentage thresholds across all columns.
Creates folders for each column type and an 'unknown' folder for entries with all zeros.
"""

import os
import sys
import shutil
import argparse
from datetime import datetime
from pathlib import Path
import re


def log_message(message):
    """Log a message with timestamp."""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")


def validate_inputs(input_dir, output_dir, table_file):
    """Validate input parameters."""
    if not os.path.isdir(input_dir):
        print(f"❌ Input directory does not exist: {input_dir}")
        sys.exit(1)
    
    if not os.path.isfile(table_file):
        print(f"❌ Gene table not found: {table_file}")
        sys.exit(1)
    
    if not os.access(table_file, os.R_OK):
        print(f"❌ Gene table is not readable: {table_file}")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)


def load_table_data(table_file):
    """Load and parse the gene table file."""
    table_data = {}
    entries_count = 0
    
    log_message("🔍 Reading table file...")
    
    try:
        with open(table_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) < 6:
                    log_message(f"⚠️ Skipping malformed line {line_num}: insufficient columns")
                    continue
                
                family, cog, total, regulatory, biosynthetic, biosynthetic_additional = parts[:6]
                
                if family == "family":
                    continue
                
                if not cog:
                    continue
                
                table_data[cog] = {
                    'total': total,
                    'regulatory': regulatory,
                    'biosynthetic': biosynthetic,
                    'biosynthetic_additional': biosynthetic_additional
                }
                entries_count += 1
                
    except Exception as e:
        log_message(f"❌ Error reading table file: {e}")
        sys.exit(1)
    
    log_message(f"📊 Loaded {entries_count} entries from table")
    return table_data


def is_numeric(value):
    """Check if a value is numeric (including decimals)."""
    try:
        float(value)
        return True
    except ValueError:
        return False


def calculate_percentages(data):
    """Calculate percentages for regulatory, biosynthetic, and biosynthetic_additional columns."""
    try:
        total = float(data['total'])
        if total == 0:
            return {'regulatory': 0, 'biosynthetic': 0, 'biosynthetic_additional': 0}
        
        regulatory_pct = (float(data['regulatory']) / total) * 100
        biosynthetic_pct = (float(data['biosynthetic']) / total) * 100
        biosynthetic_additional_pct = (float(data['biosynthetic_additional']) / total) * 100
        
        return {
            'regulatory': regulatory_pct,
            'biosynthetic': biosynthetic_pct,
            'biosynthetic_additional': biosynthetic_additional_pct
        }
    except (ValueError, ZeroDivisionError):
        return {'regulatory': 0, 'biosynthetic': 0, 'biosynthetic_additional': 0}


def determine_categories(percentages, threshold=1.0):
    """Determine which categories a sample belongs to based on percentage threshold."""
    categories = []
    
    if percentages['regulatory'] > threshold:
        categories.append('regulatory')
    if percentages['biosynthetic'] > threshold:
        categories.append('biosynthetic')
    if percentages['biosynthetic_additional'] > threshold:
        categories.append('biosynthetic_additional')
    
    if not categories:
        if all(pct == 0 for pct in percentages.values()):
            categories.append('unknown')
        else:
            # Values exist but are below threshold - still categorize as unknown
            categories.append('unknown')
    
    return categories


def extract_base_filename(dir_name):
    """Extract base filename without trailing motif mode (_anr_maxwXX or _zoops_maxwXX)"""
    base = re.sub(r'_(anr|zoops)_maxw\d+.*$', '', dir_name)
    return base


def copy_directory_to_categories(dir_path, output_dir, categories, dir_name):
    """Copy directory to all applicable category folders."""
    copied_paths = []
    
    for category in categories:
        category_dir = Path(output_dir) / category
        category_dir.mkdir(parents=True, exist_ok=True)
        
        dest_path = category_dir / dir_name
        
        try:
            if dest_path.exists():
                shutil.rmtree(dest_path)
            
            shutil.copytree(dir_path, dest_path)
            copied_paths.append(str(dest_path))
            
        except Exception as e:
            log_message(f"❌ Failed to copy directory to {dest_path}: {e}")
            return []
    
    return copied_paths


def process_directories(input_dir, output_dir, table_data, threshold=1.0):
    """Process MEME directories based on percentage criteria."""
    total_processed = 0
    total_copied = 0
    skipped_count = 0
    dirs_found = 0
    category_counts = {'regulatory': 0, 'biosynthetic': 0, 'biosynthetic_additional': 0, 'unknown': 0}
    
    log_message(f"🔍 Scanning input directory for MEME directories (threshold: {threshold}%)...")
    
    input_path = Path(input_dir)
    
    for dir_path in input_path.iterdir():
        if not dir_path.is_dir():
            continue
        
        dir_name = dir_path.name
        dirs_found += 1
        
        file_key = extract_base_filename(dir_name)
        
        if dirs_found <= 5:
            log_message(f"🔍 Directory: {dir_name} → extracted file: {file_key}")
        
        if file_key in table_data:
            data = table_data[file_key]
            
            if not all(is_numeric(data[col]) for col in ['total', 'regulatory', 'biosynthetic', 'biosynthetic_additional']):
                log_message(f"⚠️ Skipping {file_key}: non-numeric values")
                skipped_count += 1
                continue
            
            percentages = calculate_percentages(data)
            
            categories = determine_categories(percentages, threshold)
            
            log_message(f"📋 {file_key}: regulatory={percentages['regulatory']:.1f}%, "
                       f"biosynthetic={percentages['biosynthetic']:.1f}%, "
                       f"biosynthetic_additional={percentages['biosynthetic_additional']:.1f}% "
                       f"→ Categories: {', '.join(categories)}")
            
            copied_paths = copy_directory_to_categories(dir_path, output_dir, categories, dir_name)
            
            if copied_paths:
                total_processed += 1
                total_copied += len(copied_paths)
                
                for category in categories:
                    category_counts[category] += 1
                
                log_message(f"✅ Copied {dir_name} to: {', '.join(copied_paths)}")
            else:
                log_message(f"❌ Failed to copy {dir_name}")
                skipped_count += 1
        else:
            log_message(f"⚠️ No table entry found for {file_key} (from directory {dir_name})")
            skipped_count += 1
    
    return dirs_found, total_processed, total_copied, skipped_count, category_counts


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Organize MEME analysis directories based on percentage thresholds across all columns',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python meme_organizer.py /path/to/input /path/to/output gene_table.tsv
    python meme_organizer.py /path/to/input /path/to/output gene_table.tsv --threshold 2.0
        """
    )
    
    parser.add_argument('input_dir', help='Input directory containing MEME directories')
    parser.add_argument('output_dir', help='Output directory for organized results')
    parser.add_argument('table_file', help='Gene table file (TSV format)')
    parser.add_argument('--threshold', type=float, default=1.0,
                       help='Percentage threshold for categorization (default: 1.0)')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    validate_inputs(args.input_dir, args.output_dir, args.table_file)
    
    log_message("📁 Starting MEME analysis organization")
    log_message(f"📥 Input directory: {args.input_dir}")
    log_message(f"📤 Output directory: {args.output_dir}")
    log_message(f"📄 Using gene table: {args.table_file}")
    log_message(f"📊 Percentage threshold: {args.threshold}%")
    
    table_data = load_table_data(args.table_file)
    
    dirs_found, total_processed, total_copied, skipped_count, category_counts = process_directories(
        args.input_dir, args.output_dir, table_data, args.threshold
    )
    
    log_message("✅ MEME organization complete.")
    log_message(f"📊 Summary: Found {dirs_found} directories, {total_processed} processed, "
               f"{total_copied} total copies created, {skipped_count} skipped")
    
    log_message("📁 Category distribution:")
    for category, count in category_counts.items():
        if count > 0:
            log_message(f"  - {category}: {count} directories")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        log_message("❌ Script interrupted by user")
        sys.exit(1)
    except Exception as e:
        log_message(f"❌ Script failed: {e}")
        sys.exit(1)
