#!/usr/bin/env python3
"""
Common utilities for InterSubMod visualization tools.

This module provides shared functions and constants used by:
- plot_distance_heatmap.py
- plot_cluster_heatmap.py

Author: InterSubMod Development Team
Date: 2026-01-13
"""

import os
import sys
import logging
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd

# Import visualization libraries with fallback
try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ============================================================================
# Color Palettes
# ============================================================================

HP_COLORS = {
    "0": "#CCCCCC",      # Unphased - gray
    "1": "#5773CC",      # HP1 - blue
    "2": "#F0E685",      # HP2 - yellow
    "1-1": "#0A47FF",    # HP1-1 - light blue
    "2-1": "#FFB900",    # HP2-1 - orange
    "unphased": "#A9A9A9",  # Unphased - gray
}

STRAND_COLORS = {
    "+": "#D1D1D1",      # Forward - light green
    "-": "#D1DCD1",      # Reverse - light blue
    "?": "#A9A9A9",      # Unknown - light gray
}

SOURCE_COLORS = {
    "Tumor": "#FF6B6B",  # Red
    "Normal": "#27AE60", # Green
}

ALLELE_COLORS = {
    "ALT": "#949494",    # dark gray
    "REF": "#CCCCCC",    # light gray
    "UNKNOWN": "#A9A9A9",  # Light gray
}


# ============================================================================
# Logging Setup
# ============================================================================


def setup_logging(log_file: Optional[str] = None, verbose: bool = False):
    """
    Configure logging to console and optional file.

    Args:
        log_file: Optional path to log file
        verbose: If True, set DEBUG level; otherwise INFO
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if verbose else logging.INFO)

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
# Data Loading Functions
# ============================================================================


def load_reads_metadata(region_dir: str) -> Optional[pd.DataFrame]:
    """
    Load read metadata from TSV file.

    Args:
        region_dir: Path to region directory

    Returns:
        DataFrame with read_id as index and columns: hp, alt_support, is_tumor, strand
        Returns None if file not found or error occurs.
    """
    reads_file = os.path.join(region_dir, "reads", "reads.tsv")

    if not os.path.exists(reads_file):
        return None

    try:
        df = pd.read_csv(reads_file, sep="\t")
        # Set read_id as index (matches distance matrix and methylation matrix)
        if "read_id" in df.columns:
            df = df.set_index("read_id")
            df.index = df.index.astype(int)
        elif "read_name" in df.columns:
            df = df.set_index("read_name")
        return df
    except Exception as e:
        logging.debug(f"Error loading reads metadata: {e}")
        return None


def load_region_info(region_dir: str) -> Dict:
    """
    Load region metadata for plot titles.

    Args:
        region_dir: Path to region directory

    Returns:
        Dictionary with keys: region, snv, num_reads, num_cpgs
    """
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
        except Exception:
            pass

    return info


def load_distance_matrix(
    region_dir: str, metric: str = "NHD", strand: str = "all"
) -> Tuple[Optional[pd.DataFrame], Optional[List[int]]]:
    """
    Load precomputed distance matrix from CSV file.

    Args:
        region_dir: Path to region directory
        metric: Distance metric name (e.g., "NHD", "L1", "BERNOULLI")
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
        logging.debug(f"Error loading distance matrix: {e}")
        return None, None


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
        logging.debug(f"Error loading methylation matrix: {e}")
        return None, None


def load_linkage_matrix(region_dir: str, strand: str = "all") -> Optional[np.ndarray]:
    """
    Load pre-computed linkage matrix from C++ output.

    Args:
        region_dir: Path to region directory
        strand: "all", "forward", or "reverse"

    Returns:
        Scipy-compatible linkage matrix if available, None otherwise.
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
        logging.debug(f"Error loading linkage matrix: {e}")
        return None


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
        Linkage matrix (n-1 x 4), or None if computation fails
    """
    if not HAS_SCIPY:
        logging.error("scipy is required for linkage computation")
        return None

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
        logging.debug(f"Linkage computation failed: {e}")
        return None


def get_cluster_order(Z: np.ndarray) -> np.ndarray:
    """
    Get optimal leaf ordering from linkage matrix.

    Args:
        Z: Linkage matrix from hierarchical clustering

    Returns:
        Array of indices representing the optimal leaf order
    """
    if not HAS_SCIPY:
        return np.arange(len(Z) + 1)
    return leaves_list(Z)


# ============================================================================
# Annotation Functions
# ============================================================================


def create_annotation_colors(
    reads_df: pd.DataFrame, read_ids: List[int]
) -> Tuple[pd.DataFrame, Dict]:
    """
    Create annotation DataFrame and color palette for seaborn clustermap.

    Args:
        reads_df: DataFrame with read metadata (indexed by read_id)
        read_ids: List of read IDs to include

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


# ============================================================================
# Region Discovery Functions
# ============================================================================


def find_region_dirs(
    output_dir: str,
    require_methylation: bool = True,
    require_distance: bool = False,
    distance_metric: str = "NHD"
) -> List[str]:
    """
    Find all region directories in output directory.

    Args:
        output_dir: Path to output directory
        require_methylation: Only include regions with methylation.csv
        require_distance: Only include regions with distance matrix
        distance_metric: Distance metric to check for

    Returns:
        List of paths to region directories
    """
    region_dirs = []

    for item in os.listdir(output_dir):
        region_path = os.path.join(output_dir, item)

        # Check if it's a directory and looks like a region
        if not os.path.isdir(region_path):
            continue
        if not item.startswith("region_"):
            continue

        # Check required files
        if require_methylation:
            meth_file = os.path.join(region_path, "methylation", "methylation.csv")
            if not os.path.exists(meth_file):
                continue

        if require_distance:
            dist_file = os.path.join(region_path, "distance", distance_metric, "matrix.csv")
            if not os.path.exists(dist_file):
                continue

        region_dirs.append(region_path)

    return sorted(region_dirs)
