#!/usr/bin/env python3
"""
Distance-based Cluster Heatmap Visualization for InterSubMod
============================================================

This script generates **proper** cluster heatmaps from the distance matrices
produced by the InterSubMod C++ pipeline.

Key Features:
- Plots Read × Read distance matrix (NOT methylation × CpGs)
- Shows dendrogram on both axes (rows and columns)
- Annotation bars for biological labels (HP, Tumor/Normal, Strand, Allele)
- Parallel processing for multiple regions
- Supports using pre-computed linkage from C++ or computing from distance matrix

Output: Distance-based cluster heatmaps showing read similarity patterns

Author: InterSubMod Development Team
Date: 2025-12-03
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


def load_distance_matrix(
    region_dir: str, metric: str = "NHD", strand: str = "all"
) -> Tuple[Optional[pd.DataFrame], Optional[List[str]]]:
    """
    Load precomputed distance matrix from CSV file.

    Args:
        region_dir: Path to region directory
        metric: Distance metric name (e.g., "NHD", "L1")
        strand: "all", "forward", or "reverse"

    Returns:
        Tuple of (DataFrame with distances, list of read IDs)
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
        read_ids = list(df.index)
        return df, read_ids
    except Exception as e:
        print(f"Error loading distance matrix: {e}")
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
        # Set read_name as index (matches distance matrix)
        if "read_name" in df.columns:
            df = df.set_index("read_name")
        elif "read_id" in df.columns:
            df = df.set_index("read_id")
        return df
    except Exception as e:
        print(f"Error loading reads metadata: {e}")
        return None


def load_linkage_matrix(region_dir: str, strand: str = "all") -> Optional[np.ndarray]:
    """
    Load pre-computed linkage matrix from C++ output.

    Returns scipy-compatible linkage matrix if available.
    """
    clustering_dir = os.path.join(region_dir, "clustering")

    if strand == "forward":
        filepath = os.path.join(clustering_dir, "linkage_matrix_forward.csv")
    elif strand == "reverse":
        filepath = os.path.join(clustering_dir, "linkage_matrix_reverse.csv")
    else:
        filepath = os.path.join(clustering_dir, "linkage_matrix.csv")

    if not os.path.exists(filepath):
        return None

    try:
        df = pd.read_csv(filepath, sep="\t")
        # Format: cluster_i, cluster_j, distance, size
        # Convert to scipy linkage format
        Z = df[["cluster_i", "cluster_j", "distance", "size"]].values
        return Z
    except Exception as e:
        print(f"Error loading linkage matrix: {e}")
        return None


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
) -> np.ndarray:
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
        print(f"Linkage computation failed: {e}")
        return None


def get_cluster_order(Z: np.ndarray) -> np.ndarray:
    """Get optimal leaf ordering from linkage matrix."""
    return leaves_list(Z)


# ============================================================================
# Visualization Functions
# ============================================================================


def create_annotation_colors(
    reads_df: pd.DataFrame, read_ids: List[str]
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


def plot_distance_heatmap(
    dist_df: pd.DataFrame,
    read_ids: List[str],
    reads_df: pd.DataFrame,
    linkage_matrix: np.ndarray,
    output_path: str,
    region_info: Dict,
    metric_name: str = "NHD",
    linkage_method: str = "UPGMA",
    figsize: Tuple[int, int] = (12, 10),
    dpi: int = 150,
) -> bool:
    """
    Create and save distance-based cluster heatmap with annotations.

    This is the CORRECT cluster heatmap:
    - X axis: Reads (clustered)
    - Y axis: Reads (clustered)
    - Color: Distance values (0 = identical, 1 = completely different)
    - Dendrograms: Shown on both sides

    Args:
        dist_df: Square distance matrix DataFrame
        read_ids: List of read IDs
        reads_df: Read metadata DataFrame
        linkage_matrix: Hierarchical clustering linkage
        output_path: Path to save the figure
        region_info: Region metadata for title
        metric_name: Name of distance metric
        linkage_method: Name of linkage method
        figsize: Figure size (width, height)
        dpi: Output resolution

    Returns:
        True if successful, False otherwise
    """
    try:
        n_reads = len(read_ids)

        # Create annotation colors
        annotations, color_dict = create_annotation_colors(reads_df, read_ids)

        # Create row/col colors DataFrame for clustermap
        row_colors = None
        if not annotations.empty:
            row_colors_list = []
            for col in annotations.columns:
                color_series = annotations[col].map(color_dict[col])
                row_colors_list.append(color_series)
            row_colors = pd.concat(row_colors_list, axis=1)
            row_colors.columns = annotations.columns
            # Ensure index matches distance matrix
            row_colors = row_colors.reindex(dist_df.index)

        # Create clustermap with dendrograms
        g = sns.clustermap(
            dist_df,
            row_linkage=linkage_matrix,
            col_linkage=linkage_matrix,
            row_colors=row_colors,
            col_colors=row_colors,
            cmap="viridis_r",  # Reversed: dark = similar (low distance)
            vmin=0,
            vmax=1,
            figsize=figsize,
            cbar_kws={"label": f"Distance ({metric_name})"},
            dendrogram_ratio=(0.15, 0.15),
            colors_ratio=0.03 if row_colors is not None else 0,
            linewidths=0,
            xticklabels=False if n_reads > 50 else True,
            yticklabels=False if n_reads > 50 else True,
        )

        # Title
        title = f"Read-Read Distance Cluster Heatmap\n"
        title += f"SNV: {region_info.get('snv', 'Unknown')} | "
        title += f"Reads: {n_reads} | Method: {linkage_method}"
        g.fig.suptitle(title, y=1.02, fontsize=12, fontweight="bold")

        # Axis labels
        g.ax_heatmap.set_xlabel("Reads (clustered)", fontsize=10)
        g.ax_heatmap.set_ylabel("Reads (clustered)", fontsize=10)

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
        print(f"Error creating distance heatmap: {e}")
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
    figsize: Tuple[int, int] = (12, 10),
    dpi: int = 150,
) -> Tuple[str, bool, str]:
    """
    Process a single region and generate distance-based cluster heatmap.

    Args:
        region_dir: Path to region directory
        distance_metric: Distance metric to use
        linkage_method: Hierarchical clustering method
        strand: "all", "forward", or "reverse"
        output_format: Image format ("png", "pdf", "svg")
        min_reads: Minimum reads required
        figsize: Figure size
        dpi: Output resolution

    Returns:
        Tuple of (region_dir, success, message)
    """
    region_name = os.path.basename(region_dir)

    # Load distance matrix
    dist_df, read_ids = load_distance_matrix(region_dir, distance_metric, strand)
    if dist_df is None:
        return region_dir, False, f"No distance matrix found for {distance_metric}"

    n_reads = len(read_ids)

    # Check minimum requirements
    if n_reads < min_reads:
        return region_dir, False, f"Too few reads ({n_reads} < {min_reads})"

    # Load reads metadata
    reads_df = load_reads_metadata(region_dir)
    if reads_df is None:
        return region_dir, False, "No reads metadata found"

    # Try to load pre-computed linkage from C++
    Z = load_linkage_matrix(region_dir, strand)

    # If not available, compute from distance matrix
    if Z is None:
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
    output_filename = f"distance_heatmap{strand_suffix}.{output_format}"
    output_path = os.path.join(plot_dir, output_filename)

    # Create heatmap
    success = plot_distance_heatmap(
        dist_df,
        read_ids,
        reads_df,
        Z,
        output_path,
        region_info,
        metric_name=distance_metric,
        linkage_method=linkage_display,
        figsize=figsize,
        dpi=dpi,
    )

    if success:
        return region_dir, True, output_path
    else:
        return region_dir, False, "Distance heatmap generation failed"


def find_region_dirs(output_dir: str) -> List[str]:
    """Find all region directories in output directory."""
    region_dirs = []
    output_path = Path(output_dir)

    # Walk through the directory tree
    for root, dirs, files in os.walk(output_path):
        # Check if this directory is a region directory
        # For distance heatmap: needs distance folder
        if (Path(root) / "distance").exists():
            region_dirs.append(root)

    return sorted(region_dirs)


def process_all_regions(
    output_dir: str, num_threads: int = 64, **kwargs
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
                print(
                    f"Progress: {i}/{total} ({100*i/total:.1f}%) - "
                    f"Success: {success_count}, Failed: {fail_count} - "
                    f"Rate: {rate:.1f} regions/sec"
                )

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
        description="Generate distance-based cluster heatmaps from InterSubMod output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script generates the CORRECT cluster heatmap:
- Plots Read × Read distance matrix (not methylation matrix)
- Shows dendrograms on both axes
- Supports biological annotations (HP, Tumor/Normal, Strand, Allele)

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
        help="Distance metric to use (default: NHD)",
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
        default="12,10",
        help="Figure size as 'width,height' (default: 12,10)",
    )

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
        figsize = (12, 10)

    # Common kwargs
    kwargs = {
        "distance_metric": args.metric,
        "linkage_method": args.linkage,
        "strand": args.strand,
        "output_format": args.format,
        "min_reads": args.min_reads,
        "figsize": figsize,
        "dpi": args.dpi,
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
            args.output_dir, num_threads=args.threads, **kwargs
        )

        if failed == 0:
            sys.exit(0)
        else:
            sys.exit(1 if success == 0 else 0)


if __name__ == "__main__":
    main()
