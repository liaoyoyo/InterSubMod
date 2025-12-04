#!/usr/bin/env python3
"""
Cluster Heatmap Visualization for InterSubMod
==============================================

This script generates cluster heatmaps from methylation matrices and read metadata
produced by the InterSubMod C++ pipeline.

Features:
- Hierarchical clustering with UPGMA (default) or other methods
- Annotation bars for biological labels (HP, Tumor/Normal, Strand, Alt-support)
- Parallel processing for multiple regions
- Customizable output formats and parameters

Usage:
    # Single region
    python plot_cluster_heatmap.py --region-dir /path/to/region_dir
    
    # All regions in output directory (parallel)
    python plot_cluster_heatmap.py --output-dir /path/to/output --threads 64
    
    # With custom parameters
    python plot_cluster_heatmap.py --output-dir /path/to/output --linkage ward --metric euclidean

Author: InterSubMod Development Team
Date: 2025-12-02
"""

import os
import sys
import argparse
import warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional, Tuple, Dict, List
import time

import numpy as np
import pandas as pd

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Import visualization libraries with fallback
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for server
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import squareform
    HAS_VIZ = True
except ImportError as e:
    HAS_VIZ = False
    IMPORT_ERROR = str(e)


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_methylation_matrix(region_dir: str, strand: str = "all") -> Tuple[Optional[pd.DataFrame], Optional[np.ndarray]]:
    """
    Load methylation matrix from CSV file.
    
    Args:
        region_dir: Path to region directory
        strand: "all", "forward", or "reverse"
    
    Returns:
        Tuple of (DataFrame with read_id index, numpy array of CpG positions)
    """
    meth_dir = os.path.join(region_dir, "methylation")
    
    if strand == "forward":
        filepath = os.path.join(meth_dir, "methylation_forward.csv")
    elif strand == "reverse":
        filepath = os.path.join(meth_dir, "methylation_reverse.csv")
    else:
        filepath = os.path.join(meth_dir, "methylation.csv")
    
    if not os.path.exists(filepath):
        return None, None
    
    try:
        df = pd.read_csv(filepath, index_col=0)
        
        # Handle strand-specific files which have original_read_id column
        if strand in ["forward", "reverse"] and "original_read_id" in df.columns:
            df = df.drop(columns=["original_read_id"])
        
        # Get CpG positions from column names
        cpg_positions = np.array([int(col) for col in df.columns])
        
        # Replace "NA" strings with NaN
        df = df.replace("NA", np.nan)
        df = df.astype(float)
        
        return df, cpg_positions
    except Exception as e:
        print(f"Error loading methylation matrix: {e}")
        return None, None


def load_reads_metadata(region_dir: str) -> Optional[pd.DataFrame]:
    """
    Load read metadata from TSV file.
    
    Returns DataFrame with columns: read_id, hp, alt_support, is_tumor, strand
    """
    reads_file = os.path.join(region_dir, "reads", "reads.tsv")
    
    if not os.path.exists(reads_file):
        return None
    
    try:
        df = pd.read_csv(reads_file, sep="\t")
        df = df.set_index("read_id")
        return df
    except Exception as e:
        print(f"Error loading reads metadata: {e}")
        return None


def load_distance_matrix(region_dir: str, metric: str = "NHD", strand: str = "all") -> Optional[np.ndarray]:
    """
    Load precomputed distance matrix from CSV file.
    
    Args:
        region_dir: Path to region directory
        metric: Distance metric name (e.g., "NHD", "L1")
        strand: "all", "forward", or "reverse"
    
    Returns:
        Square distance matrix as numpy array
    """
    dist_dir = os.path.join(region_dir, "distance", metric)
    
    if strand == "forward":
        filepath = os.path.join(dist_dir, "matrix_forward.csv")
    elif strand == "reverse":
        filepath = os.path.join(dist_dir, "matrix_reverse.csv")
    else:
        filepath = os.path.join(dist_dir, "matrix.csv")
    
    if not os.path.exists(filepath):
        return None
    
    try:
        df = pd.read_csv(filepath, index_col=0)
        return df.values.astype(float)
    except Exception as e:
        print(f"Error loading distance matrix: {e}")
        return None


def load_region_info(region_dir: str) -> Dict:
    """Load region metadata for plot titles."""
    metadata_file = os.path.join(region_dir, "metadata.txt")
    info = {"region": "Unknown", "snv": "Unknown", "num_reads": 0, "num_cpgs": 0}
    
    if os.path.exists(metadata_file):
        try:
            with open(metadata_file, 'r') as f:
                for line in f:
                    if line.startswith("Region:"):
                        info["region"] = line.split(":")[1].strip()
                    elif line.startswith("SNV Position:"):
                        info["snv"] = line.split(":", 1)[1].strip()
                    elif line.startswith("Num Reads:"):
                        info["num_reads"] = int(line.split(":")[1].strip())
                    elif line.startswith("Num CpG Sites:"):
                        info["num_cpgs"] = int(line.split(":")[1].strip())
        except:
            pass
    
    return info


# ============================================================================
# Clustering Functions
# ============================================================================

def compute_linkage(distance_matrix: np.ndarray, method: str = "average") -> np.ndarray:
    """
    Compute hierarchical clustering linkage from distance matrix.
    
    Args:
        distance_matrix: Square distance matrix
        method: Linkage method ("average" for UPGMA, "ward", "single", "complete")
    
    Returns:
        Linkage matrix (n-1 x 4)
    """
    # Handle NaN values
    dist_clean = np.nan_to_num(distance_matrix, nan=1.0)
    
    # Ensure symmetry
    dist_clean = (dist_clean + dist_clean.T) / 2
    np.fill_diagonal(dist_clean, 0)
    
    # Convert to condensed form
    try:
        condensed = squareform(dist_clean, checks=False)
        Z = linkage(condensed, method=method)
        return Z
    except Exception as e:
        print(f"Linkage computation failed: {e}")
        return None


def get_cluster_order(Z: np.ndarray) -> np.ndarray:
    """Get optimal leaf ordering from linkage matrix."""
    return leaves_list(Z)


# ============================================================================
# Visualization Functions
# ============================================================================

def create_annotation_colors(reads_df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict]:
    """
    Create annotation DataFrame and color palette for seaborn clustermap.
    
    Returns:
        Tuple of (annotation DataFrame, color dictionary)
    """
    annotations = pd.DataFrame(index=reads_df.index)
    colors = {}
    
    # Haplotype (HP) annotation
    if "hp" in reads_df.columns:
        hp_values = reads_df["hp"].astype(str)
        annotations["HP"] = hp_values
        # Color palette for HP
        unique_hp = hp_values.unique()
        hp_palette = {
            "0": "#CCCCCC",      # Unphased - gray
            "1": "#E74C3C",      # HP1 - red
            "2": "#3498DB",      # HP2 - blue
            "1-1": "#E74C3C",    # HP1-1 - red
            "1-2": "#9B59B6",    # HP1-2 - purple
            "2-1": "#2ECC71",    # HP2-1 - green
            "2-2": "#3498DB",    # HP2-2 - blue
        }
        colors["HP"] = {hp: hp_palette.get(hp, "#AAAAAA") for hp in unique_hp}
    
    # Strand annotation
    if "strand" in reads_df.columns:
        annotations["Strand"] = reads_df["strand"]
        colors["Strand"] = {
            "+": "#FF6B6B",      # Forward - coral red
            "-": "#4ECDC4",      # Reverse - teal
            "?": "#95A5A6"       # Unknown - gray
        }
    
    # Tumor/Normal annotation
    if "is_tumor" in reads_df.columns:
        annotations["Source"] = reads_df["is_tumor"].map({1: "Tumor", 0: "Normal", "1": "Tumor", "0": "Normal"})
        colors["Source"] = {
            "Tumor": "#E74C3C",   # Red
            "Normal": "#27AE60"   # Green
        }
    
    # Alt/Ref support annotation
    if "alt_support" in reads_df.columns:
        annotations["Allele"] = reads_df["alt_support"]
        colors["Allele"] = {
            "ALT": "#F39C12",     # Orange
            "REF": "#8E44AD",     # Purple
            "UNKNOWN": "#BDC3C7" # Light gray
        }
    
    return annotations, colors


def plot_cluster_heatmap(
    meth_matrix: pd.DataFrame,
    reads_df: pd.DataFrame,
    linkage_matrix: np.ndarray,
    output_path: str,
    region_info: Dict,
    figsize: Tuple[int, int] = (14, 10),
    cmap: str = "RdYlBu_r",
    show_dendrogram: bool = True,
    dpi: int = 150
) -> bool:
    """
    Create and save cluster heatmap with annotations.
    
    Args:
        meth_matrix: Methylation matrix (reads x CpGs)
        reads_df: Read metadata DataFrame
        linkage_matrix: Hierarchical clustering linkage
        output_path: Path to save the figure
        region_info: Region metadata for title
        figsize: Figure size (width, height)
        cmap: Colormap for heatmap
        show_dendrogram: Whether to show row dendrogram
        dpi: Output resolution
    
    Returns:
        True if successful, False otherwise
    """
    try:
        # Get cluster order
        order = get_cluster_order(linkage_matrix)
        
        # Reorder data
        meth_ordered = meth_matrix.iloc[order]
        reads_ordered = reads_df.loc[meth_ordered.index]
        
        # Create annotation colors
        annotations, color_dict = create_annotation_colors(reads_ordered)
        
        # Create row colors for clustermap
        row_colors = None
        if not annotations.empty:
            row_colors_list = []
            for col in annotations.columns:
                row_colors_list.append(annotations[col].map(color_dict[col]))
            row_colors = pd.concat(row_colors_list, axis=1)
            row_colors.columns = annotations.columns
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        # Compute mask for NaN values
        mask = meth_ordered.isna()
        
        # Fill NaN with 0.5 for visualization (will be masked)
        data_filled = meth_ordered.fillna(0.5)
        
        # Create clustermap
        g = sns.clustermap(
            data_filled,
            row_cluster=False,  # Already ordered by our linkage
            col_cluster=False,  # Keep CpG order by position
            row_colors=row_colors,
            mask=mask,
            cmap=cmap,
            vmin=0, vmax=1,
            figsize=figsize,
            xticklabels=True,
            yticklabels=False,
            cbar_kws={'label': 'Methylation Level'},
            dendrogram_ratio=(0.1, 0) if show_dendrogram else (0, 0)
        )
        
        # Add title
        title = f"Cluster Heatmap: {region_info.get('snv', 'Unknown')}\n"
        title += f"Reads: {region_info.get('num_reads', 0)}, CpGs: {region_info.get('num_cpgs', 0)}"
        g.fig.suptitle(title, y=1.02, fontsize=12, fontweight='bold')
        
        # Rotate x-axis labels
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        
        # Add legend for annotations
        if row_colors is not None:
            handles = []
            labels = []
            for col in annotations.columns:
                for val, color in color_dict[col].items():
                    if val in annotations[col].values:
                        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color))
                        labels.append(f"{col}: {val}")
            
            if handles:
                g.ax_heatmap.legend(
                    handles, labels,
                    loc='upper left',
                    bbox_to_anchor=(1.02, 1),
                    fontsize=8,
                    title='Annotations'
                )
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        g.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        plt.close(g.fig)
        
        return True
        
    except Exception as e:
        print(f"Error creating heatmap: {e}")
        return False


# ============================================================================
# Main Processing Functions
# ============================================================================

def process_single_region(
    region_dir: str,
    distance_metric: str = "NHD",
    linkage_method: str = "average",
    strand: str = "all",
    output_format: str = "png",
    min_reads: int = 10,
    min_cpgs: int = 3,
    figsize: Tuple[int, int] = (14, 10),
    dpi: int = 150
) -> Tuple[str, bool, str]:
    """
    Process a single region and generate cluster heatmap.
    
    Args:
        region_dir: Path to region directory
        distance_metric: Distance metric to use
        linkage_method: Hierarchical clustering method
        strand: "all", "forward", or "reverse"
        output_format: Image format ("png", "pdf", "svg")
        min_reads: Minimum reads required
        min_cpgs: Minimum CpG sites required
        figsize: Figure size
        dpi: Output resolution
    
    Returns:
        Tuple of (region_dir, success, message)
    """
    region_name = os.path.basename(region_dir)
    
    # Load methylation matrix
    meth_df, cpg_positions = load_methylation_matrix(region_dir, strand)
    if meth_df is None:
        return region_dir, False, "No methylation matrix found"
    
    # Check minimum requirements
    if len(meth_df) < min_reads:
        return region_dir, False, f"Too few reads ({len(meth_df)} < {min_reads})"
    
    if len(cpg_positions) < min_cpgs:
        return region_dir, False, f"Too few CpGs ({len(cpg_positions)} < {min_cpgs})"
    
    # Load reads metadata
    reads_df = load_reads_metadata(region_dir)
    if reads_df is None:
        return region_dir, False, "No reads metadata found"
    
    # Load distance matrix
    dist_matrix = load_distance_matrix(region_dir, distance_metric, strand)
    if dist_matrix is None:
        return region_dir, False, f"No distance matrix found for {distance_metric}"
    
    # Handle size mismatch (strand-specific may have different sizes)
    if dist_matrix.shape[0] != len(meth_df):
        return region_dir, False, f"Size mismatch: distance({dist_matrix.shape[0]}) vs methylation({len(meth_df)})"
    
    # Compute linkage
    Z = compute_linkage(dist_matrix, method=linkage_method)
    if Z is None:
        return region_dir, False, "Linkage computation failed"
    
    # Load region info
    region_info = load_region_info(region_dir)
    
    # Align reads_df with meth_df
    common_idx = meth_df.index.intersection(reads_df.index)
    if len(common_idx) < min_reads:
        return region_dir, False, f"Index alignment failed"
    
    meth_df = meth_df.loc[common_idx]
    reads_df = reads_df.loc[common_idx]
    
    # Create output directory
    plot_dir = os.path.join(region_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # Generate output filename
    strand_suffix = f"_{strand}" if strand != "all" else ""
    output_filename = f"cluster_heatmap{strand_suffix}.{output_format}"
    output_path = os.path.join(plot_dir, output_filename)
    
    # Create heatmap
    success = plot_cluster_heatmap(
        meth_df, reads_df, Z, output_path,
        region_info, figsize=figsize, dpi=dpi
    )
    
    if success:
        return region_dir, True, output_path
    else:
        return region_dir, False, "Heatmap generation failed"


def find_region_dirs(output_dir: str) -> List[str]:
    """Find all region directories in output directory."""
    region_dirs = []
    output_path = Path(output_dir)
    
    # Walk through the directory tree
    for root, dirs, files in os.walk(output_path):
        # Check if this directory is a region directory
        # For cluster heatmap: needs methylation/methylation.csv
        if (Path(root) / "methylation" / "methylation.csv").exists():
            region_dirs.append(root)
            
    return sorted(region_dirs)


def process_all_regions(
    output_dir: str,
    num_threads: int = 64,
    **kwargs
) -> Tuple[int, int, float]:
    """
    Process all regions in parallel.
    
    Returns:
        Tuple of (success_count, fail_count, elapsed_time)
    """
    region_dirs = find_region_dirs(output_dir)
    total = len(region_dirs)
    
    if total == 0:
        print(f"No region directories found in {output_dir}")
        return 0, 0, 0.0
    
    print(f"Found {total} regions to process")
    print(f"Using {num_threads} threads")
    
    start_time = time.time()
    success_count = 0
    fail_count = 0
    
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = {
            executor.submit(process_single_region, region_dir, **kwargs): region_dir
            for region_dir in region_dirs
        }
        
        for i, future in enumerate(as_completed(futures), 1):
            region_dir, success, message = future.result()
            
            if success:
                success_count += 1
            else:
                fail_count += 1
            
            # Progress update every 100 regions
            if i % 100 == 0 or i == total:
                elapsed = time.time() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                print(f"Progress: {i}/{total} ({100*i/total:.1f}%) - "
                      f"Success: {success_count}, Failed: {fail_count} - "
                      f"Rate: {rate:.1f} regions/sec")
    
    elapsed_time = time.time() - start_time
    print(f"\nCompleted in {elapsed_time:.1f} seconds")
    print(f"Success: {success_count}/{total} ({100*success_count/total:.1f}%)")
    print(f"Failed: {fail_count}/{total}")
    
    return success_count, fail_count, elapsed_time


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate cluster heatmaps from InterSubMod output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single region
  %(prog)s --region-dir /path/to/chr1_12345/chr1_11345_13345
  
  # Process all regions in output directory
  %(prog)s --output-dir /path/to/output --threads 64
  
  # With custom parameters
  %(prog)s --output-dir /path/to/output --linkage ward --min-reads 20
        """
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--region-dir", type=str,
                            help="Path to single region directory")
    input_group.add_argument("--output-dir", type=str,
                            help="Path to output directory with multiple regions")
    
    # Processing options
    parser.add_argument("-t", "--threads", type=int, default=64,
                       help="Number of parallel threads (default: 64)")
    parser.add_argument("--metric", type=str, default="NHD",
                       help="Distance metric to use (default: NHD)")
    parser.add_argument("--linkage", type=str, default="average",
                       choices=["average", "single", "complete", "ward"],
                       help="Linkage method (default: average/UPGMA)")
    parser.add_argument("--strand", type=str, default="all",
                       choices=["all", "forward", "reverse"],
                       help="Which strand to plot (default: all)")
    
    # Filter options
    parser.add_argument("--min-reads", type=int, default=10,
                       help="Minimum reads required (default: 10)")
    parser.add_argument("--min-cpgs", type=int, default=3,
                       help="Minimum CpG sites required (default: 3)")
    
    # Output options
    parser.add_argument("--format", type=str, default="png",
                       choices=["png", "pdf", "svg"],
                       help="Output image format (default: png)")
    parser.add_argument("--dpi", type=int, default=150,
                       help="Output resolution (default: 150)")
    parser.add_argument("--figsize", type=str, default="14,10",
                       help="Figure size as 'width,height' (default: 14,10)")
    
    args = parser.parse_args()
    
    # Check dependencies
    if not HAS_VIZ:
        print(f"ERROR: Required visualization libraries not available: {IMPORT_ERROR}")
        print("Please install: pip install matplotlib seaborn scipy pandas numpy")
        sys.exit(1)
    
    # Parse figure size
    try:
        figsize = tuple(map(int, args.figsize.split(",")))
    except:
        figsize = (14, 10)
    
    # Common kwargs
    kwargs = {
        "distance_metric": args.metric,
        "linkage_method": args.linkage,
        "strand": args.strand,
        "output_format": args.format,
        "min_reads": args.min_reads,
        "min_cpgs": args.min_cpgs,
        "figsize": figsize,
        "dpi": args.dpi
    }
    
    if args.region_dir:
        # Single region mode
        print(f"Processing single region: {args.region_dir}")
        region_dir, success, message = process_single_region(args.region_dir, **kwargs)
        
        if success:
            print(f"✓ Success: {message}")
            sys.exit(0)
        else:
            print(f"✗ Failed: {message}")
            sys.exit(1)
    else:
        # Batch mode
        print(f"Processing all regions in: {args.output_dir}")
        success, failed, elapsed = process_all_regions(
            args.output_dir,
            num_threads=args.threads,
            **kwargs
        )
        
        if failed == 0:
            sys.exit(0)
        else:
            sys.exit(1 if success == 0 else 0)


if __name__ == "__main__":
    main()

