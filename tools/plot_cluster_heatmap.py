#!/usr/bin/env python3
"""
Cluster Heatmap Visualization for InterSubMod
==============================================

This script generates methylation cluster heatmaps from methylation matrices
with dendrogram visualization based on read-read distance clustering.

Features:
- Methylation heatmap (Reads × CpGs) with Y-axis dendrogram
- Reads are sorted by hierarchical clustering (using distance matrix)
- CpG sites are ordered by genomic position
- Annotation bars for biological labels (HP, Tumor/Normal, Strand, Alt-support)
- Parallel processing for multiple regions

Key difference from distance_heatmap:
- distance_heatmap: Read × Read distance matrix
- cluster_heatmap: Read × CpG methylation matrix with clustering-ordered Y axis

Usage:
    # Single region
    python plot_cluster_heatmap.py --region-dir /path/to/region_dir

    # All regions in output directory (parallel)
    python plot_cluster_heatmap.py --output-dir /path/to/output --threads 64

Author: InterSubMod Development Team
Date: 2025-12-18 (Updated with dendrogram support)
"""

import os
import sys
import argparse
import warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional, Tuple, Dict, List
import time
import logging


import numpy as np
import pandas as pd

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# Import visualization libraries with fallback
try:
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend for server
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import squareform

    HAS_VIZ = True
except ImportError as e:
    HAS_VIZ = False
    IMPORT_ERROR = str(e)


# ============================================================================
# Logging Setup
# ============================================================================


def setup_logging(log_file: Optional[str] = None, verbose: bool = False):
    """
    Configure logging to console and optional file.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    
    # Remove existing handlers to avoid duplicates
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console Handler
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # File Handler
    if log_file:
        try:
            file_handler = logging.FileHandler(log_file, mode='w')
            file_handler.setFormatter(formatter)
            root_logger.addHandler(file_handler)
            logging.info(f"Logging to file: {log_file}")
        except Exception as e:
            logging.error(f"Failed to setup file logging to {log_file}: {e}")


# ============================================================================
# Color Palettes
# ============================================================================

HP_COLORS = {
    "0": "#CCCCCC",  # Unphased - gray
    "1": "#E74C3C",  # HP1 - red
    "2": "#3498DB",  # HP2 - blue
    "1-1": "#E74C3C",  # HP1-1 - red
    "1-2": "#9B59B6",  # HP1-2 - purple
    "2-1": "#2ECC71",  # HP2-1 - green
    "2-2": "#3498DB",  # HP2-2 - blue
    "unphased": "#CCCCCC",
}

STRAND_COLORS = {
    "+": "#FF6B6B",  # Forward - coral red
    "-": "#4ECDC4",  # Reverse - teal
    "?": "#95A5A6",  # Unknown - gray
}

SOURCE_COLORS = {"Tumor": "#E74C3C", "Normal": "#27AE60"}  # Red  # Green

ALLELE_COLORS = {
    "ALT": "#F39C12",  # Orange
    "REF": "#8E44AD",  # Purple
    "UNKNOWN": "#BDC3C7",  # Light gray
}


# ============================================================================
# Data Loading Functions
# ============================================================================


def load_methylation_matrix(
    region_dir: str, strand: str = "all"
) -> Tuple[Optional[pd.DataFrame], Optional[np.ndarray]]:
    """
    Load methylation matrix from CSV file.

    Args:
        region_dir: Path to region directory
        strand: "all", "forward", or "reverse"

    Returns:
        Tuple of (DataFrame with read_id as integer index, numpy array of CpG positions)
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

        # Ensure index is integer (read_id)
        df.index = df.index.astype(int)

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
        # We rely on the caller to handle failures, but we can log debug info if needed
        # logging.debug(f"Error loading methylation matrix: {e}") 
        # (Commented out to avoid noise in subprocess if logging is not configured same way)
        return None, None


def load_reads_metadata(region_dir: str) -> Optional[pd.DataFrame]:
    """
    Load read metadata from TSV file.

    Returns DataFrame with read_id (integer) as index and columns: hp, alt_support, is_tumor, strand
    """
    reads_file = os.path.join(region_dir, "reads", "reads.tsv")

    if not os.path.exists(reads_file):
        return None

    try:
        df = pd.read_csv(reads_file, sep="\t")
        # Set read_id as index (integer ID matching methylation/distance matrices)
        if "read_id" in df.columns:
            df = df.set_index("read_id")
            df.index = df.index.astype(int)
        return df
    except Exception as e:
        return None


def load_distance_matrix(
    region_dir: str, metric: str = "NHD", strand: str = "all"
) -> Tuple[Optional[pd.DataFrame], Optional[List[int]]]:
    """
    Load precomputed distance matrix from CSV file.

    Args:
        region_dir: Path to region directory
        metric: Distance metric name (e.g., "NHD", "L1")
        strand: "all", "forward", or "reverse"

    Returns:
        Tuple of (DataFrame with distances, list of read IDs as integers)
    """
    dist_dir = os.path.join(region_dir, "distance", metric)

    if strand == "forward":
        filepath = os.path.join(dist_dir, "matrix_forward.csv")
    elif strand == "reverse":
        filepath = os.path.join(dist_dir, "matrix_reverse.csv")
    else:
        filepath = os.path.join(dist_dir, "matrix.csv")

    if not os.path.exists(filepath):
        return None, None

    try:
        df = pd.read_csv(filepath, index_col=0)
        # Ensure both index and columns are integers (read_id)
        df.index = df.index.astype(int)
        df.columns = df.columns.astype(int)
        read_ids = list(df.index)
        return df, read_ids
    except Exception as e:
        return None, None


def load_region_info(region_dir: str) -> Dict:
    """Load region metadata for plot titles."""
    metadata_file = os.path.join(region_dir, "metadata.txt")
    info = {"region": "Unknown", "snv": "Unknown", "num_reads": 0, "num_cpgs": 0}

    if os.path.exists(metadata_file):
        try:
            with open(metadata_file, "r") as f:
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


def compute_linkage_from_distance(
    dist_df: pd.DataFrame, method: str = "average"
) -> Optional[np.ndarray]:
    """
    Compute hierarchical clustering linkage from distance matrix DataFrame.

    Args:
        dist_df: Square distance matrix as DataFrame
        method: Linkage method ("average" for UPGMA, "ward", "single", "complete")

    Returns:
        Linkage matrix (n-1 x 4)
    """
    # Handle NaN values
    dist_clean = dist_df.values.copy()
    dist_clean = np.nan_to_num(dist_clean, nan=1.0)

    # Ensure symmetry
    dist_clean = (dist_clean + dist_clean.T) / 2
    np.fill_diagonal(dist_clean, 0)

    # Convert to condensed form
    try:
        condensed = squareform(dist_clean, checks=False)
        Z = linkage(condensed, method=method)
        return Z
    except Exception as e:
        # print(f"Linkage computation failed: {e}")
        return None


def get_cluster_order(Z: np.ndarray) -> np.ndarray:
    """Get optimal leaf ordering from linkage matrix."""
    return leaves_list(Z)


# ============================================================================
# Visualization Functions
# ============================================================================


def create_annotation_colors(
    reads_df: pd.DataFrame, read_ids: List[int]
) -> Tuple[pd.DataFrame, Dict]:
    """
    Create annotation DataFrame and color palette for seaborn clustermap.

    Returns:
        Tuple of (annotation DataFrame aligned with read_ids, color dictionary)
    """
    # Filter to only include reads in our list
    available_reads = [r for r in read_ids if r in reads_df.index]
    if not available_reads:
        return pd.DataFrame(), {}

    reads_filtered = reads_df.loc[available_reads]
    annotations = pd.DataFrame(index=available_reads)
    colors = {}

    # Haplotype (HP) annotation
    if "hp" in reads_filtered.columns:
        hp_values = reads_filtered["hp"].astype(str)
        annotations["HP"] = hp_values
        unique_hp = hp_values.unique()
        colors["HP"] = {hp: HP_COLORS.get(hp, "#AAAAAA") for hp in unique_hp}

    # Strand annotation
    if "strand" in reads_filtered.columns:
        annotations["Strand"] = reads_filtered["strand"]
        unique_strand = reads_filtered["strand"].unique()
        colors["Strand"] = {s: STRAND_COLORS.get(s, "#95A5A6") for s in unique_strand}

    # Tumor/Normal annotation
    if "is_tumor" in reads_filtered.columns:
        annotations["Source"] = reads_filtered["is_tumor"].map(
            {
                1: "Tumor",
                0: "Normal",
                "1": "Tumor",
                "0": "Normal",
                True: "Tumor",
                False: "Normal",
            }
        )
        unique_source = annotations["Source"].unique()
        colors["Source"] = {s: SOURCE_COLORS.get(s, "#95A5A6") for s in unique_source}

    # Alt/Ref support annotation
    if "alt_support" in reads_filtered.columns:
        annotations["Allele"] = reads_filtered["alt_support"]
        unique_allele = reads_filtered["alt_support"].unique()
        colors["Allele"] = {a: ALLELE_COLORS.get(a, "#BDC3C7") for a in unique_allele}

    return annotations, colors


def plot_cluster_heatmap(
    meth_df: pd.DataFrame,
    reads_df: pd.DataFrame,
    linkage_matrix: np.ndarray,
    output_path: str,
    region_info: Dict,
    cpg_positions: np.ndarray,
    metric_name: str = "NHD",
    linkage_method: str = "UPGMA",
    figsize: Tuple[int, int] = (14, 10),
    cmap: str = "RdYlBu_r",
    dpi: int = 150,
) -> bool:
    """
    Create and save methylation cluster heatmap with Y-axis dendrogram.

    This heatmap shows:
    - X axis: CpG sites (ordered by genomic position, no clustering)
    - Y axis: Reads (clustered by distance matrix, with dendrogram)
    - Color: Methylation level (0 = unmethylated, 1 = methylated)

    Args:
        meth_df: Methylation matrix (reads × CpGs)
        reads_df: Read metadata DataFrame
        linkage_matrix: Hierarchical clustering linkage from distance matrix
        output_path: Path to save the figure
        region_info: Region metadata for title
        cpg_positions: Array of CpG genomic positions
        metric_name: Name of distance metric used for clustering
        linkage_method: Name of linkage method
        figsize: Figure size (width, height)
        cmap: Colormap for heatmap
        dpi: Output resolution

    Returns:
        True if successful, False otherwise
    """
    try:
        n_reads = len(meth_df)
        n_cpgs = len(cpg_positions)
        read_ids = list(meth_df.index)

        # Create annotation colors
        annotations, color_dict = create_annotation_colors(reads_df, read_ids)

        # Create row colors for clustermap
        row_colors = None
        if not annotations.empty:
            row_colors_list = []
            for col in annotations.columns:
                color_series = annotations[col].map(color_dict[col])
                row_colors_list.append(color_series)
            row_colors = pd.concat(row_colors_list, axis=1)
            row_colors.columns = annotations.columns
            # Ensure index matches methylation matrix
            row_colors = row_colors.reindex(meth_df.index)

        # Prepare methylation data (handle NaN for visualization)
        # NaN values will be shown as gray using mask
        data_for_plot = meth_df.copy()

        # Create clustermap with row dendrogram only
        # row_linkage: uses our pre-computed linkage from distance matrix
        # col_cluster=False: keep CpG positions in genomic order
        g = sns.clustermap(
            data_for_plot,
            row_linkage=linkage_matrix,  # Use distance-based linkage for Y-axis
            col_cluster=False,  # Keep CpGs in position order
            row_colors=row_colors,
            cmap=cmap,
            vmin=0,
            vmax=1,
            figsize=figsize,
            mask=data_for_plot.isna(),  # Mask NaN values (show as gray)
            cbar_kws={"label": "Methylation Level"},
            dendrogram_ratio=(0.15, 0),  # Show row dendrogram, no column dendrogram
            colors_ratio=0.03 if row_colors is not None else 0,
            linewidths=0,
            xticklabels=True if n_cpgs <= 50 else False,
            yticklabels=False,  # Too many reads to show labels
        )

        # Set title
        title = f"Methylation Cluster Heatmap\n"
        title += f"SNV: {region_info.get('snv', 'Unknown')} | "
        title += f"Reads: {n_reads} | CpGs: {n_cpgs} | Cluster: {linkage_method}"
        g.fig.suptitle(title, y=1.02, fontsize=12, fontweight="bold")

        # Axis labels
        g.ax_heatmap.set_xlabel("CpG Sites (genomic position)", fontsize=10)
        g.ax_heatmap.set_ylabel("Reads (clustered)", fontsize=10)

        # Rotate x-axis labels if shown
        if n_cpgs <= 50:
            plt.setp(
                g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7
            )

        # Add legend for annotations
        if row_colors is not None and not annotations.empty:
            handles = []
            labels = []
            for col in annotations.columns:
                for val, color in color_dict[col].items():
                    if val in annotations[col].values:
                        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color))
                        labels.append(f"{col}: {val}")

            if handles:
                # Position legend outside the heatmap
                g.fig.legend(
                    handles,
                    labels,
                    loc="upper left",
                    bbox_to_anchor=(0.02, 0.98),
                    fontsize=8,
                    framealpha=0.9,
                    title="Annotations",
                )

        # Adjust layout
        plt.tight_layout()

        # Save figure
        g.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        plt.close(g.fig)

        return True

    except Exception as e:
        # Return False so the caller can log the error message
        return False
        import traceback

        traceback.print_exc()
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
    dpi: int = 150,
) -> Tuple[str, bool, str]:
    """
    Process a single region and generate methylation cluster heatmap.

    Args:
        region_dir: Path to region directory
        distance_metric: Distance metric to use for clustering
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

    n_reads = len(meth_df)
    n_cpgs = len(cpg_positions) if cpg_positions is not None else 0

    # Check minimum requirements
    if n_reads < min_reads:
        return region_dir, False, f"Too few reads ({n_reads} < {min_reads})"

    if n_cpgs < min_cpgs:
        return region_dir, False, f"Too few CpGs ({n_cpgs} < {min_cpgs})"

    # Load distance matrix (for clustering)
    dist_df, dist_read_ids = load_distance_matrix(region_dir, distance_metric, strand)
    if dist_df is None:
        return region_dir, False, f"No distance matrix found for {distance_metric}"

    # Load reads metadata
    reads_df = load_reads_metadata(region_dir)
    if reads_df is None:
        return region_dir, False, "No reads metadata found"

    # Align methylation matrix with distance matrix (use common reads)
    common_reads = list(set(meth_df.index) & set(dist_df.index))
    if len(common_reads) < min_reads:
        return (
            region_dir,
            False,
            f"Too few common reads ({len(common_reads)} < {min_reads})",
        )

    # Reorder to match distance matrix order (important for linkage alignment)
    ordered_reads = [r for r in dist_df.index if r in common_reads]
    meth_df = meth_df.loc[ordered_reads]
    dist_df = dist_df.loc[ordered_reads, ordered_reads]

    # Compute linkage from distance matrix
    Z = compute_linkage_from_distance(dist_df, method=linkage_method)
    if Z is None:
        return region_dir, False, "Linkage computation failed"

    # Load region info
    region_info = load_region_info(region_dir)

    # Map linkage method to display name
    linkage_display = {
        "average": "UPGMA",
        "ward": "Ward",
        "single": "Single",
        "complete": "Complete",
    }.get(linkage_method, linkage_method.upper())

    # Create output directory
    plot_dir = os.path.join(region_dir, "plots", distance_metric)
    os.makedirs(plot_dir, exist_ok=True)

    # Generate output filename
    strand_suffix = f"_{strand}" if strand != "all" else ""
    output_filename = f"cluster_heatmap{strand_suffix}.{output_format}"
    output_path = os.path.join(plot_dir, output_filename)

    # Create heatmap
    success = plot_cluster_heatmap(
        meth_df,
        reads_df,
        Z,
        output_path,
        region_info,
        cpg_positions,
        metric_name=distance_metric,
        linkage_method=linkage_display,
        figsize=figsize,
        dpi=dpi,
    )

    if success:
        return region_dir, True, output_path
    else:
        return region_dir, False, "Cluster heatmap generation failed"


def find_region_dirs(output_dir: str) -> List[str]:
    """Find all region directories in output directory."""
    region_dirs = []
    output_path = Path(output_dir)

    # Walk through the directory tree
    for root, dirs, files in os.walk(output_path):
        # Check if this directory is a region directory
        # Need both methylation and distance folder for cluster heatmap
        if (Path(root) / "methylation" / "methylation.csv").exists() and (
            Path(root) / "distance"
        ).exists():
            region_dirs.append(root)

    return sorted(region_dirs)


def process_all_regions(
    output_dir: str, num_threads: int = 64, log_file: Optional[str] = None, **kwargs
) -> Tuple[int, int, float]:
    """
    Process all regions in parallel.

    Returns:
        Tuple of (success_count, fail_count, elapsed_time)
    """
    region_dirs = find_region_dirs(output_dir)
    total = len(region_dirs)

    if total == 0:
        logging.warning(f"No region directories found in {output_dir}")
        return 0, 0, 0.0

    logging.info(f"Found {total} regions to process")
    logging.info(f"Using {num_threads} threads")

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
                logging.info(
                    f"Progress: {i}/{total} ({100*i/total:.1f}%) - "
                    f"Success: {success_count}, Failed: {fail_count} - "
                    f"Rate: {rate:.1f} regions/sec"
                )
            
            # Log failures individually
            if not success:
               logging.warning(f"FAILED {region_dir}: {message}")

    elapsed_time = time.time() - start_time
    logging.info(f"Completed in {elapsed_time:.1f} seconds")
    logging.info(f"Success: {success_count}/{total} ({100*success_count/total:.1f}%)")
    logging.info(f"Failed: {fail_count}/{total}")

    return success_count, fail_count, elapsed_time


# ============================================================================
# Main Entry Point
# ============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Generate methylation cluster heatmaps with dendrogram from InterSubMod output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script generates methylation cluster heatmaps:
- X axis: CpG sites (ordered by genomic position)
- Y axis: Reads (clustered by distance matrix, with dendrogram)
- Color: Methylation level (0=unmethylated, 1=methylated)

The Y-axis dendrogram shows read clustering based on the distance matrix,
allowing visualization of how methylation patterns relate to read similarity.

Examples:
  # Process single region
  %(prog)s --region-dir /path/to/chr1_12345/chr1_11345_13345
  
  # Process all regions in output directory
  %(prog)s --output-dir /path/to/output --threads 64
  
  # With custom parameters
  %(prog)s --output-dir /path/to/output --linkage ward --min-reads 20
        """,
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--region-dir", type=str, help="Path to single region directory"
    )
    input_group.add_argument(
        "--output-dir", type=str, help="Path to output directory with multiple regions"
    )

    # Processing options
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=64,
        help="Number of parallel threads (default: 64)",
    )
    parser.add_argument(
        "--metric",
        type=str,
        default="NHD",
        help="Distance metric to use for clustering (default: NHD)",
    )
    parser.add_argument(
        "--linkage",
        type=str,
        default="average",
        choices=["average", "single", "complete", "ward"],
        help="Linkage method (default: average/UPGMA)",
    )
    parser.add_argument(
        "--strand",
        type=str,
        default="all",
        choices=["all", "forward", "reverse"],
        help="Which strand to plot (default: all)",
    )

    # Filter options
    parser.add_argument(
        "--min-reads", type=int, default=10, help="Minimum reads required (default: 10)"
    )
    parser.add_argument(
        "--min-cpgs",
        type=int,
        default=3,
        help="Minimum CpG sites required (default: 3)",
    )

    # Output options
    parser.add_argument(
        "--format",
        type=str,
        default="png",
        choices=["png", "pdf", "svg"],
        help="Output image format (default: png)",
    )
    parser.add_argument(
        "--dpi", type=int, default=150, help="Output resolution (default: 150)"
    )
    parser.add_argument(
        "--figsize",
        type=str,
        default="14,10",
        help="Figure size as 'width,height' (default: 14,10)",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file (default: heatmap_generation.log in output dir)",
    )

    args = parser.parse_args()

    # Setup logging
    log_file = args.log_file
    if not log_file and args.output_dir:
        # Default log file in output directory
        log_file = os.path.join(args.output_dir, "heatmap_generation.log")
    elif not log_file:
         # Fallback default
         log_file = "heatmap_generation.log"
         
    setup_logging(log_file)

    # Check dependencies
    if not HAS_VIZ:
        logging.error(f"Required visualization libraries not available: {IMPORT_ERROR}")
        logging.error("Please install: pip install matplotlib seaborn scipy pandas numpy")
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
        "dpi": args.dpi,
    }

    if args.region_dir:
        # Single region mode
        logging.info(f"Processing single region: {args.region_dir}")
        region_dir, success, message = process_single_region(args.region_dir, **kwargs)

        if success:
            logging.info(f"✓ Success: {message}")
            sys.exit(0)
        else:
            logging.error(f"✗ Failed: {message}")
            sys.exit(1)
    else:
        # Batch mode
        logging.info(f"Processing all regions in: {args.output_dir}")
        success, failed, elapsed = process_all_regions(
            args.output_dir, num_threads=args.threads, log_file=log_file, **kwargs
        )

        if failed == 0:
            sys.exit(0)
        else:
            sys.exit(1 if success == 0 else 0)


if __name__ == "__main__":
    main()
