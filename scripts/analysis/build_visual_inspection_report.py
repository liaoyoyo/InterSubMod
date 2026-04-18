#!/usr/bin/env python3
"""M6: Visual inspection report — ISM per-region visualization + IGV batch (Wave 3 Script G).

Samples 40 representative sites (4 categories × TP/FP × 5 each) from HCC1395 TO mode,
generates methylation heatmaps, distance matrix heatmaps, and read structure plots
from ISM per-region output. Also prepares IGV batch script for screenshots.

Excludes centromere and telomere regions.

Output: {OUTPUT_ROOT}/20260406_visual_inspection/
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")

# Configure CJK font support (same as tools/plot_*.py)
_tools_dir = str(Path(__file__).resolve().parent.parent / "tools")
sys.path.insert(0, _tools_dir)
try:
    from font_config import configure_matplotlib_fonts
    configure_matplotlib_fonts()
except ImportError:
    pass

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, leaves_list, fcluster
from scipy.cluster.hierarchy import linkage as compute_linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(Path(__file__).parent))
from observation_common import (
    OUTPUT_ROOT, SAMPLE_ORDER,
    load_master_dataset, setup_plot_style, save_figure, ensure_dir,
    encode_truth_binary,
)

# ---------------------------------------------------------------------------
# ISM Standard Color Palettes (from tools/plot_distance_heatmap.py)
# ---------------------------------------------------------------------------

HP_COLORS = {
    "0": "#CCCCCC",       # Unphased - gray
    "1": "#5773CC",       # HP1 - blue
    "2": "#F0E685",       # HP2 - yellow
    "1-1": "#0A47FF",     # HP1-1 - deep blue
    "2-1": "#FFB900",     # HP2-1 - orange
    "unphased": "#A9A9A9",
}

STRAND_COLORS = {
    "+": "#D1D1D1",
    "-": "#D1DCD1",
    "?": "#A9A9A9",
}

SOURCE_COLORS = {"Tumor": "#FF6B6B", "Normal": "#27AE60"}

ALLELE_COLORS = {
    "ALT": "#949494",     # dark gray
    "REF": "#CCCCCC",     # light gray
    "UNKNOWN": "#A9A9A9",
}

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_visual_inspection")

# ISM output paths
ISM_ROOT_TP = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260315_hcc1395_to_pilot/step05_intersubmod/intersubmod_tp/filtered_snv_tp"
)
ISM_ROOT_FP = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260315_hcc1395_to_pilot/step05_intersubmod/intersubmod_fp/filtered_snv_fp"  # 注意：FP 可能路徑不同
)

# IGV
IGV_EXECUTABLE = "/bip7_disk/IGV_Linux_2.11.1/igv.sh"
IGV_SESSION_TEMPLATE = "/big8_disk/liaoyoyo2001/IGV_session/template.xml"

# hg38 centromere coordinates (UCSC, approximate)
# Source: https://genome.ucsc.edu/cgi-bin/hgTables (group: Mapping, track: Centromeres)
CENTROMERES_HG38 = {
    "chr1": (121700000, 125100000),
    "chr2": (91800000, 96000000),
    "chr3": (87800000, 94000000),
    "chr4": (48200000, 51800000),
    "chr5": (46100000, 51400000),
    "chr6": (58500000, 62600000),
    "chr7": (58100000, 62100000),
    "chr8": (43200000, 47200000),
    "chr9": (42200000, 45500000),
    "chr10": (38000000, 41600000),
    "chr11": (51000000, 55800000),
    "chr12": (33200000, 37800000),
    "chr13": (16500000, 18900000),
    "chr14": (16100000, 18200000),
    "chr15": (17500000, 20500000),
    "chr16": (35300000, 38400000),
    "chr17": (22700000, 27400000),
    "chr18": (15400000, 21500000),
    "chr19": (24200000, 28100000),
    "chr20": (25700000, 30400000),
    "chr21": (10900000, 13000000),
    "chr22": (13700000, 17400000),
    "chrX": (58100000, 63800000),
    "chrY": (10300000, 10600000),
}

# hg38 telomere regions (first/last 500kb of each chromosome, approximate)
TELOMERE_MARGIN = 500_000

# Chromosome sizes for telomere end detection
CHROM_SIZES_HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def is_centromere_or_telomere(chrom: str, pos: int) -> bool:
    """Check if position is in centromere or telomere region."""
    # Centromere
    if chrom in CENTROMERES_HG38:
        cen_start, cen_end = CENTROMERES_HG38[chrom]
        if cen_start <= pos <= cen_end:
            return True

    # Telomere (first/last TELOMERE_MARGIN bp)
    if pos <= TELOMERE_MARGIN:
        return True
    if chrom in CHROM_SIZES_HG38:
        if pos >= CHROM_SIZES_HG38[chrom] - TELOMERE_MARGIN:
            return True

    return False


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

CATEGORIES = {
    "loh_high_allele": {
        "desc": "LOH + AlleleDelta high (>75th percentile)",
        "filter": lambda df, ad_q75: (
            df["is_loh"] &
            (pd.to_numeric(df["AlleleDelta"], errors="coerce") > ad_q75)
        ),
    },
    "loh_low_allele": {
        "desc": "LOH + AlleleDelta low (<25th percentile) — hopeless",
        "filter": lambda df, ad_q25: (
            df["is_loh"] &
            (pd.to_numeric(df["AlleleDelta"], errors="coerce") < ad_q25)
        ),
    },
    "nonloh_hp_sig": {
        "desc": "Non-LOH + HPMergedSig=True",
        "filter": lambda df, _: (
            ~df["is_loh"] &
            df["HPMergedSig"].astype(str).str.lower().isin(["true", "1", "yes"])
        ),
    },
    "nonloh_no_signal": {
        "desc": "Non-LOH + HPSig=F + AlleleSig=F — hardest",
        "filter": lambda df, _: (
            ~df["is_loh"] &
            ~df["HPMergedSig"].astype(str).str.lower().isin(["true", "1", "yes"]) &
            ~df["AlleleSig"].astype(str).str.lower().isin(["true", "1", "yes"])
        ),
    },
}

N_PER_CATEGORY = 5


def sample_sites(df: pd.DataFrame) -> pd.DataFrame:
    """Sample representative sites from HCC1395 TO mode."""
    print("\n" + "=" * 60)
    print("Sampling representative sites")
    print("=" * 60)

    # Filter: HCC1395, TO mode, valid truth label
    df = df.copy()
    df["is_loh"] = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1", "yes"])

    hcc = df[(df["sample"] == "HCC1395") & (df["mode"] == "to")].copy()
    hcc["truth_binary"] = encode_truth_binary(hcc)
    hcc = hcc[hcc["truth_binary"].notna()]

    # Read count and CpG filters
    nr = pd.to_numeric(hcc["NumReads"], errors="coerce")
    nc = pd.to_numeric(hcc["NumCpGs"], errors="coerce")
    hcc = hcc[(nr >= 30) & (nr <= 100) & (nc >= 10)]

    # Exclude centromere and telomere
    mask_ct = hcc.apply(
        lambda row: is_centromere_or_telomere(str(row["Chr"]), int(row["Pos"])),
        axis=1,
    )
    n_excluded = mask_ct.sum()
    hcc = hcc[~mask_ct]
    print(f"  Excluded {n_excluded} centromere/telomere sites")
    print(f"  Candidate pool: {len(hcc):,} sites (HCC1395 TO, 30≤NR≤100, NC≥10)")

    # AlleleDelta quantiles for LOH categories
    loh_sites = hcc[hcc["is_loh"]]
    ad = pd.to_numeric(loh_sites["AlleleDelta"], errors="coerce").dropna()
    ad_q25 = ad.quantile(0.25)
    ad_q75 = ad.quantile(0.75)
    print(f"  AlleleDelta quantiles (LOH): Q25={ad_q25:.4f}, Q75={ad_q75:.4f}")

    all_sampled = []

    for cat_name, cat_info in CATEGORIES.items():
        # Apply category filter
        if "high" in cat_name:
            mask = cat_info["filter"](hcc, ad_q75)
        elif "low" in cat_name:
            mask = cat_info["filter"](hcc, ad_q25)
        else:
            mask = cat_info["filter"](hcc, None)

        pool = hcc[mask]
        tp_pool = pool[pool["truth_label"] == "TP"]
        fp_pool = pool[pool["truth_label"] == "FP"]

        print(f"\n  {cat_name}: {cat_info['desc']}")
        print(f"    Pool: {len(pool):,} (TP={len(tp_pool):,}, FP={len(fp_pool):,})")

        # Sample
        n_tp = min(N_PER_CATEGORY, len(tp_pool))
        n_fp = min(N_PER_CATEGORY, len(fp_pool))

        tp_sampled = tp_pool.sample(n=n_tp, random_state=42) if n_tp > 0 else pd.DataFrame()
        fp_sampled = fp_pool.sample(n=n_fp, random_state=42) if n_fp > 0 else pd.DataFrame()

        for i, (_, row) in enumerate(tp_sampled.iterrows()):
            row_copy = row.copy()
            row_copy["category"] = cat_name
            row_copy["label"] = f"tp_{i+1:02d}"
            row_copy["site_id"] = f"tp_{i+1:02d}_{row['Chr']}_{int(row['Pos'])}"
            all_sampled.append(row_copy)

        for i, (_, row) in enumerate(fp_sampled.iterrows()):
            row_copy = row.copy()
            row_copy["category"] = cat_name
            row_copy["label"] = f"fp_{i+1:02d}"
            row_copy["site_id"] = f"fp_{i+1:02d}_{row['Chr']}_{int(row['Pos'])}"
            all_sampled.append(row_copy)

        print(f"    Sampled: TP={n_tp}, FP={n_fp}")

    sampled_df = pd.DataFrame(all_sampled)
    sampled_df.to_csv(OUT_DIR / "sampled_positions.tsv", sep="\t", index=False)
    print(f"\n  Total sampled: {len(sampled_df)} sites")
    return sampled_df


# ---------------------------------------------------------------------------
# ISM Region Data Reader
# ---------------------------------------------------------------------------

def find_region_dir(chrom: str, pos: int, is_tp: bool) -> Optional[Path]:
    """Find ISM per-region output directory for a given position."""
    root = ISM_ROOT_TP if is_tp else ISM_ROOT_FP

    # Check if FP root exists, if not try alternative path
    if not root.exists():
        # Try without filtered_snv prefix
        alt = root.parent / "filtered_snv_fp"
        if alt.exists():
            root = alt
        else:
            # Try intersubmod_fp directory
            fp_parent = root.parent.parent / "intersubmod_fp"
            if fp_parent.exists():
                for sub in fp_parent.iterdir():
                    if sub.is_dir() and "filtered" in sub.name:
                        root = sub
                        break

    chr_dir = root / chrom
    if not chr_dir.exists():
        return None

    pos_dir = chr_dir / f"{chrom}_{pos}"
    if not pos_dir.exists():
        return None

    # Find the region subdirectory (chr_start_end)
    for sub in pos_dir.iterdir():
        if sub.is_dir() and sub.name.startswith(chrom):
            return sub

    return None


def read_metadata(region_dir: Path) -> Dict:
    """Read metadata.txt from ISM region output."""
    meta_file = region_dir / "metadata.txt"
    if not meta_file.exists():
        return {}
    result = {}
    for line in meta_file.read_text().splitlines():
        if ":" in line:
            key, _, val = line.partition(":")
            result[key.strip()] = val.strip()
    return result


def read_reads_tsv(region_dir: Path) -> Optional[pd.DataFrame]:
    """Read reads.tsv from ISM region output (read_id as integer index)."""
    reads_file = region_dir / "reads" / "reads.tsv"
    if not reads_file.exists():
        return None
    df = pd.read_csv(reads_file, sep="\t")
    if "read_id" in df.columns:
        df = df.set_index("read_id")
        df.index = df.index.astype(int)
    return df


def read_methylation_matrix(region_dir: Path) -> Tuple[Optional[pd.DataFrame], Optional[np.ndarray]]:
    """Read methylation.csv from ISM region output (read_id as index, CpG positions as columns).

    Returns:
        Tuple of (DataFrame with read_id integer index, CpG positions array)
    """
    meth_file = region_dir / "methylation" / "methylation.csv"
    if not meth_file.exists():
        return None, None
    df = pd.read_csv(meth_file, index_col=0)
    df.index = df.index.astype(int)
    df = df.replace("NA", np.nan).astype(float)
    cpg_positions = np.array([int(col) for col in df.columns])
    return df, cpg_positions


def read_significance(region_dir: Path) -> Dict:
    """Read significance.json from ISM region output."""
    sig_file = region_dir / "clustering" / "significance.json"
    if not sig_file.exists():
        return {}
    return json.loads(sig_file.read_text())


def read_linkage_matrix(region_dir: Path) -> Optional[np.ndarray]:
    """Read linkage_matrix.csv from ISM region output.

    ISM outputs 5 columns with header:
        cluster_i, cluster_j, distance, new_cluster_id, size
    scipy.hierarchy.dendrogram expects 4 columns:
        [idx1, idx2, distance, count]
    We read with header, then select columns [cluster_i, cluster_j, distance, size].
    """
    link_file = region_dir / "clustering" / "linkage_matrix.csv"
    if not link_file.exists():
        return None
    df = pd.read_csv(link_file, sep="\t")
    if df.empty:
        return None
    # Map ISM columns to scipy linkage format: [idx1, idx2, dist, count]
    return df[["cluster_i", "cluster_j", "distance", "size"]].values.astype(float)


def read_distance_matrix(region_dir: Path) -> Optional[pd.DataFrame]:
    """Read BERNOULLI distance matrix from ISM region output (read_id as index/columns)."""
    dist_dir = region_dir / "distance" / "BERNOULLI"
    if not dist_dir.exists():
        return None
    matrix_file = dist_dir / "matrix.csv"
    if not matrix_file.exists():
        return None
    df = pd.read_csv(matrix_file, index_col=0)
    df.index = df.index.astype(int)
    df.columns = df.columns.astype(int)
    df = df.replace("NA", np.nan).astype(float)
    return df


# ---------------------------------------------------------------------------
# Annotation helpers (from tools/plot_distance_heatmap.py standard)
# ---------------------------------------------------------------------------

def create_annotation_colors(
    reads_df: pd.DataFrame, read_ids: List[int],
) -> Tuple[pd.DataFrame, Dict]:
    """Create annotation DataFrame and colour palette for sns.clustermap row_colors.

    Produces annotation bars for HP, Strand, Source (Tumor/Normal), Allele.
    Same logic as tools/plot_distance_heatmap.py.
    """
    available = [r for r in read_ids if r in reads_df.index]
    if not available:
        return pd.DataFrame(), {}

    sub = reads_df.loc[available]
    ann = pd.DataFrame(index=available)
    colors: Dict[str, Dict] = {}

    # HP
    if "hp" in sub.columns:
        hp_vals = sub["hp"].astype(str)
        ann["HP"] = hp_vals
        colors["HP"] = {h: HP_COLORS.get(h, "#AAAAAA") for h in hp_vals.unique()}

    # Strand
    if "strand" in sub.columns:
        ann["Strand"] = sub["strand"]
        colors["Strand"] = {s: STRAND_COLORS.get(s, "#95A5A6") for s in sub["strand"].unique()}

    # Source (Tumor/Normal)
    if "is_tumor" in sub.columns:
        ann["Source"] = sub["is_tumor"].map(
            {1: "Tumor", 0: "Normal", "1": "Tumor", "0": "Normal",
             True: "Tumor", False: "Normal"})
        colors["Source"] = {s: SOURCE_COLORS.get(s, "#95A5A6") for s in ann["Source"].unique()}

    # Allele
    if "alt_support" in sub.columns:
        ann["Allele"] = sub["alt_support"]
        colors["Allele"] = {a: ALLELE_COLORS.get(a, "#BDC3C7") for a in sub["alt_support"].unique()}

    return ann, colors


def _build_row_colors(annotations: pd.DataFrame, color_dict: Dict, idx) -> Optional[pd.DataFrame]:
    """Convert annotation + color_dict into a row_colors DataFrame aligned to *idx*."""
    if annotations.empty:
        return None
    rc_list = []
    for col in annotations.columns:
        rc_list.append(annotations[col].map(color_dict[col]))
    rc = pd.concat(rc_list, axis=1)
    rc.columns = annotations.columns
    return rc.reindex(idx)


def _compute_linkage(dist_df: pd.DataFrame) -> Optional[np.ndarray]:
    """Compute scipy linkage from distance DataFrame (NaN→1.0, ensure symmetry)."""
    mat = dist_df.values.copy()
    mat = np.nan_to_num(mat, nan=1.0)
    mat = (mat + mat.T) / 2
    np.fill_diagonal(mat, 0)
    try:
        condensed = squareform(mat, checks=False)
        return compute_linkage(condensed, method="average")
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Visualization (using ISM standard tools/plot_*.py approach)
# ---------------------------------------------------------------------------

def plot_methylation_heatmap(
    meth_df: pd.DataFrame,
    cpg_positions: np.ndarray,
    reads_df: Optional[pd.DataFrame],
    dist_df: Optional[pd.DataFrame],
    title: str,
    out_path: Path,
):
    """V1: Methylation cluster heatmap — same as tools/plot_cluster_heatmap.py.

    Uses sns.clustermap with:
      - Y-axis dendrogram from distance-matrix linkage (reads clustered)
      - X-axis: CpG positions in genomic order (no clustering)
      - Annotation bars: HP, Strand, Source, Allele
      - Colormap: RdYlBu_r (ISM standard)
    """
    setup_plot_style()

    n_reads, n_cpgs = meth_df.shape
    read_ids = list(meth_df.index)

    # Compute linkage from distance matrix for Y-axis ordering
    Z = None
    if dist_df is not None:
        common = sorted(set(meth_df.index) & set(dist_df.index))
        if len(common) >= 2:
            ordered = [r for r in dist_df.index if r in common]
            meth_df = meth_df.loc[ordered]
            dist_sub = dist_df.loc[ordered, ordered]
            read_ids = list(meth_df.index)
            Z = _compute_linkage(dist_sub)

    # Annotation bars
    row_colors = None
    ann_colors: Dict = {}
    if reads_df is not None:
        ann, ann_colors = create_annotation_colors(reads_df, read_ids)
        row_colors = _build_row_colors(ann, ann_colors, meth_df.index)

    data_for_plot = meth_df.copy()

    # Build clustermap kwargs
    # cbar_pos: right side of figure, vertically centered — avoids overlap with dendrogram
    cm_kw = dict(
        data=data_for_plot,
        col_cluster=False,            # CpGs in genomic order
        row_colors=row_colors,
        cmap="RdYlBu_r",             # ISM standard colormap
        vmin=0, vmax=1,
        figsize=(max(12, n_cpgs * 0.14), max(8, n_reads * 0.12)),
        mask=data_for_plot.isna(),
        cbar_kws={"label": "Methylation Level"},
        cbar_pos=(0.98, 0.3, 0.015, 0.4),
        dendrogram_ratio=(0.15, 0),
        colors_ratio=0.03 if row_colors is not None else 0,
        linewidths=0,
        xticklabels=n_cpgs <= 50,
        yticklabels=False,
    )
    if Z is not None:
        cm_kw["row_linkage"] = Z
    else:
        cm_kw["row_cluster"] = False  # No distance data → no clustering

    g = sns.clustermap(**cm_kw)

    g.fig.suptitle(title, y=1.02, fontsize=11, fontweight="bold")
    g.ax_heatmap.set_xlabel(f"CpG Sites (n={n_cpgs})", fontsize=10)
    g.ax_heatmap.set_ylabel(f"Reads (n={n_reads})", fontsize=10)

    if n_cpgs <= 50:
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7)

    # Legend for annotations — placed below the figure (horizontal, no overlap)
    if row_colors is not None and ann_colors:
        handles, labels = [], []
        for col_name, cmap_dict in ann_colors.items():
            for val, color in cmap_dict.items():
                handles.append(plt.Rectangle((0, 0), 1, 1, fc=color))
                labels.append(f"{col_name}: {val}")
        if handles:
            g.fig.legend(handles, labels, loc="upper center",
                         bbox_to_anchor=(0.5, -0.01), fontsize=7,
                         ncol=min(len(handles), 5),
                         framealpha=0.9, title="Annotations",
                         title_fontsize=8)

    g.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(g.fig)
    print(f"  [V1] {out_path.name} ({out_path.stat().st_size / 1024:.0f} KB)")


def plot_distance_heatmap(
    dist_df: Optional[pd.DataFrame],
    reads_df: Optional[pd.DataFrame],
    title: str,
    out_path: Path,
):
    """V2: Distance matrix heatmap — same as tools/plot_distance_heatmap.py.

    Uses sns.clustermap with:
      - Dendrograms on both axes (symmetric)
      - Annotation bars: HP, Strand, Source, Allele
      - Colormap: viridis_r (ISM standard — dark=similar)
    """
    setup_plot_style()

    if dist_df is None or dist_df.shape[0] < 2:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "Distance matrix not available", ha="center", va="center")
        ax.set_title(title)
        save_figure(fig, out_path)
        return

    read_ids = list(dist_df.index)
    n = len(read_ids)

    Z = _compute_linkage(dist_df)
    if Z is None:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "Linkage computation failed", ha="center", va="center")
        ax.set_title(title)
        save_figure(fig, out_path)
        return

    # Annotation bars
    row_colors = None
    ann_colors: Dict = {}
    if reads_df is not None:
        ann, ann_colors = create_annotation_colors(reads_df, read_ids)
        row_colors = _build_row_colors(ann, ann_colors, dist_df.index)

    # cbar_pos: right side of figure — avoids overlap with annotation bars
    g = sns.clustermap(
        dist_df,
        row_linkage=Z,
        col_linkage=Z,
        row_colors=row_colors,
        col_colors=None,
        cmap="viridis_r",             # ISM standard — dark=similar (low distance)
        vmin=0, vmax=1,
        figsize=(max(12, n * 0.18), max(10, n * 0.15)),
        cbar_kws={"label": "Bernoulli Distance"},
        dendrogram_ratio=(0, 0.15),
        cbar_pos=(0.98, 0.3, 0.015, 0.4),
        colors_ratio=0.03 if row_colors is not None else 0,
        linewidths=0,
        xticklabels=n <= 50,
        yticklabels=n <= 50,
    )

    g.fig.suptitle(title, y=1.02, fontsize=11, fontweight="bold")
    g.ax_heatmap.set_xlabel("Reads (clustered)", fontsize=10)
    g.ax_heatmap.set_ylabel("Reads (clustered)", fontsize=10)

    # Legend for annotations — placed below the figure (horizontal, no overlap)
    if row_colors is not None and ann_colors:
        handles, labels = [], []
        for col_name, cmap_dict in ann_colors.items():
            for val, color in cmap_dict.items():
                handles.append(plt.Rectangle((0, 0), 1, 1, fc=color))
                labels.append(f"{col_name}: {val}")
        if handles:
            g.fig.legend(handles, labels, loc="upper center",
                         bbox_to_anchor=(0.5, -0.01), fontsize=7,
                         ncol=min(len(handles), 5),
                         framealpha=0.9, title="Annotations",
                         title_fontsize=8)

    g.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(g.fig)
    print(f"  [V2] {out_path.name} ({out_path.stat().st_size / 1024:.0f} KB)")


def plot_read_structure(
    reads_df: pd.DataFrame,
    title: str,
    out_path: Path,
):
    """V3: Read structure overview — HP tag × ALT/REF stacked bar (ISM standard colors)."""
    setup_plot_style()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    data = reads_df.copy()
    hp_str = data["hp"].astype(str)
    data["hp_label"] = hp_str

    # HP tag distribution using ISM standard palette
    ax = axes[0]
    hp_counts = hp_str.value_counts()
    # Sort: 0/unphased first, then 1, 1-1, 2, 2-1
    order = [h for h in ["0", "1", "1-1", "2", "2-1"] if h in hp_counts.index]
    order += [h for h in hp_counts.index if h not in order]
    hp_counts = hp_counts.reindex(order).dropna()

    bars = ax.bar(
        [str(k) for k in hp_counts.index],
        hp_counts.values,
        color=[HP_COLORS.get(str(k), "#AAAAAA") for k in hp_counts.index],
    )
    for bar, count in zip(bars, hp_counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(int(count)), ha="center", fontsize=9)
    ax.set_xlabel("HP Tag")
    ax.set_ylabel("Count")
    ax.set_title("HP Tag Distribution")

    # ALT/REF by HP tag using ISM standard palette
    ax = axes[1]
    if "alt_support" in data.columns:
        ct = pd.crosstab(data["hp_label"], data["alt_support"])
        ct = ct.reindex(order).dropna(how="all")
        allele_cols = [c for c in ["ALT", "REF", "UNKNOWN"] if c in ct.columns]
        ct[allele_cols].plot(
            kind="bar", stacked=True, ax=ax,
            color=[ALLELE_COLORS.get(c, "#BDC3C7") for c in allele_cols],
        )
        ax.set_xlabel("HP Tag")
        ax.set_ylabel("Count")
        ax.set_title("ALT/REF by HP Tag")
        ax.legend(title="Allele")
    else:
        ax.text(0.5, 0.5, "alt_support not available", ha="center")

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    save_figure(fig, out_path)


def compute_clustering_metrics(
    dist_matrix: Optional[pd.DataFrame],
    reads_df: Optional[pd.DataFrame],
    cluster_k: int = 2,
) -> Dict:
    """Compute quantitative clustering structure metrics.

    Returns dict with:
    - silhouette: silhouette score (-1 to 1, >0.25 = structure)
    - n_clusters: optimal k (2 or 3)
    - cluster_sizes: list of cluster sizes
    - hp_cluster_purity: how well clusters separate HP tags (0-1)
    - allele_cluster_purity: how well clusters separate ALT/REF (0-1)
    - structure_verdict: 'CLEAR' / 'WEAK' / 'NONE'
    """
    result = {
        "silhouette": None, "n_clusters": None, "cluster_sizes": [],
        "hp_cluster_purity": None, "allele_cluster_purity": None,
        "structure_verdict": "N/A",
    }

    if dist_matrix is None:
        return result

    mat = dist_matrix.values.astype(float)
    n = mat.shape[0]
    if n < 4:
        result["structure_verdict"] = "TOO_FEW_READS"
        return result

    try:
        condensed = squareform(mat, checks=False)
        Z = compute_linkage(condensed, method="average")
    except Exception:
        return result

    # Try k=2 and k=3, pick best silhouette
    best_sil = -1.0
    best_k = 2
    best_labels = None
    for k in [2, 3]:
        labels = fcluster(Z, t=k, criterion="maxclust")
        if len(set(labels)) < 2:
            continue
        sil = silhouette_score(mat, labels, metric="precomputed")
        if sil > best_sil:
            best_sil = sil
            best_k = k
            best_labels = labels

    if best_labels is None:
        result["structure_verdict"] = "NONE"
        return result

    result["silhouette"] = round(best_sil, 4)
    result["n_clusters"] = best_k
    from collections import Counter
    result["cluster_sizes"] = sorted(Counter(best_labels).values(), reverse=True)

    # Align reads_df to distance matrix read order
    dist_read_ids = list(dist_matrix.index)
    if reads_df is not None:
        avail = [r for r in dist_read_ids if r in reads_df.index]
        reads_aligned = reads_df.loc[avail] if avail else None
    else:
        reads_aligned = None

    # HP purity: max fraction of dominant HP tag per cluster
    if reads_aligned is not None and "hp" in reads_aligned.columns:
        hp_parsed = reads_aligned["hp"].astype(str).str.split("-").str[0]
        hp_parsed = pd.to_numeric(hp_parsed, errors="coerce").fillna(0).astype(int)
        purities = []
        for c in set(best_labels):
            mask = best_labels == c
            hp_vals = hp_parsed.values[:n][mask]
            if len(hp_vals) > 0:
                counts = Counter(hp_vals)
                purities.append(max(counts.values()) / sum(counts.values()))
        result["hp_cluster_purity"] = round(np.mean(purities), 4) if purities else None

    # Allele purity
    if reads_aligned is not None and "alt_support" in reads_aligned.columns:
        purities = []
        for c in set(best_labels):
            mask = best_labels == c
            allele_vals = reads_aligned["alt_support"].values[:n][mask]
            if len(allele_vals) > 0:
                counts = Counter(allele_vals)
                purities.append(max(counts.values()) / sum(counts.values()))
        result["allele_cluster_purity"] = round(np.mean(purities), 4) if purities else None

    # Verdict
    if best_sil > 0.40:
        result["structure_verdict"] = "CLEAR"
    elif best_sil > 0.20:
        result["structure_verdict"] = "WEAK"
    else:
        result["structure_verdict"] = "NONE"

    return result


def plot_label_vs_cluster(
    dist_matrix: Optional[pd.DataFrame],
    reads_df: Optional[pd.DataFrame],
    metrics: Dict,
    title: str,
    out_path: Path,
    cluster_k: int = 2,
):
    """V4: Cluster vs HP + Cluster vs Allele Support (ISM standard colors).
    Also displays clustering metrics.
    """
    setup_plot_style()

    if dist_matrix is None or reads_df is None:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, "Data not available", ha="center", va="center")
        ax.set_title(title)
        save_figure(fig, out_path)
        return

    n = dist_matrix.shape[0]
    read_ids = list(dist_matrix.index)

    # Compute clusters from distance matrix
    Z = _compute_linkage(dist_matrix)
    if Z is not None:
        k = metrics.get("n_clusters", cluster_k) or cluster_k
        clusters = fcluster(Z, t=k, criterion="maxclust")
    else:
        clusters = np.ones(n, dtype=int)

    # Align reads_df to distance matrix read_ids
    available = [r for r in read_ids if r in reads_df.index]
    data = reads_df.loc[available].copy()
    data["cluster"] = clusters[:len(available)]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5),
                             gridspec_kw={"width_ratios": [1, 1, 0.8]})

    # Panel 1: Cluster vs HP (ISM standard palette)
    hp_str = data["hp"].astype(str)
    data["hp_label"] = hp_str
    hp_ct = pd.crosstab(data["cluster"], data["hp_label"])
    hp_order = [h for h in ["0", "1", "1-1", "2", "2-1"] if h in hp_ct.columns]
    hp_order += [h for h in hp_ct.columns if h not in hp_order]
    hp_ct = hp_ct[hp_order]
    hp_ct.plot(kind="bar", stacked=True, ax=axes[0],
               color=[HP_COLORS.get(c, "#AAAAAA") for c in hp_ct.columns])
    axes[0].set_title("Cluster vs HP Tag", fontsize=10)
    axes[0].set_xlabel("Cluster")
    axes[0].set_ylabel("Count")
    axes[0].legend(title="HP", fontsize=8)

    # Panel 2: Cluster vs Allele Support (ISM standard palette)
    if "alt_support" in data.columns:
        allele_ct = pd.crosstab(data["cluster"], data["alt_support"])
        allele_order = [a for a in ["ALT", "REF", "UNKNOWN"] if a in allele_ct.columns]
        allele_ct = allele_ct[allele_order]
        allele_ct.plot(kind="bar", stacked=True, ax=axes[1],
                       color=[ALLELE_COLORS.get(c, "#BDC3C7") for c in allele_ct.columns])
    axes[1].set_title("Cluster vs Allele Support", fontsize=10)
    axes[1].set_xlabel("Cluster")
    axes[1].set_ylabel("Count")
    axes[1].legend(title="Allele", fontsize=8)

    # Panel 3: Metrics text
    ax_m = axes[2]
    ax_m.axis("off")
    sil = metrics.get("silhouette", "N/A")
    verdict = metrics.get("structure_verdict", "N/A")
    hp_pur = metrics.get("hp_cluster_purity", "N/A")
    al_pur = metrics.get("allele_cluster_purity", "N/A")
    sizes = metrics.get("cluster_sizes", [])

    text_lines = [
        f"Clustering Metrics",
        f"",
        f"Silhouette: {sil}",
        f"Verdict: {verdict}",
        f"K: {metrics.get('n_clusters', 'N/A')}",
        f"Sizes: {sizes}",
        f"HP purity: {hp_pur}",
        f"Allele purity: {al_pur}",
    ]
    verdict_color = {"CLEAR": "green", "WEAK": "orange", "NONE": "red"}.get(verdict, "gray")
    ax_m.text(0.1, 0.95, "\n".join(text_lines), transform=ax_m.transAxes,
              fontsize=10, verticalalignment="top", fontfamily="monospace",
              bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.8))
    ax_m.text(0.5, 0.05, verdict, transform=ax_m.transAxes,
              fontsize=16, fontweight="bold", color=verdict_color,
              ha="center", va="bottom")

    fig.suptitle(title, fontsize=11, y=1.02)
    fig.tight_layout()
    save_figure(fig, out_path)


# ---------------------------------------------------------------------------
# Nearest TP lookup
# ---------------------------------------------------------------------------

def build_tp_position_index(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    """Build per-chromosome sorted arrays of TP positions for fast nearest-TP lookup.

    Args:
        df: Master dataset (must have Chr, Pos, truth_label, sample, mode columns).

    Returns:
        Dict mapping chromosome name to sorted numpy array of TP positions.
    """
    hcc_to = df[(df["sample"] == "HCC1395") & (df["mode"] == "to")].copy()
    tp_rows = hcc_to[hcc_to["truth_label"] == "TP"]
    idx: Dict[str, np.ndarray] = {}
    for chrom, grp in tp_rows.groupby("Chr"):
        positions = np.sort(grp["Pos"].astype(int).unique())
        idx[str(chrom)] = positions
    return idx


def find_nearest_tp(
    chrom: str, pos: int, is_tp: bool, tp_index: Dict[str, np.ndarray],
) -> Tuple[Optional[int], Optional[int]]:
    """Find the nearest TP position to a given site.

    For TP sites, returns the nearest *other* TP (excluding self).
    For FP sites, returns the nearest TP.

    Returns:
        (nearest_tp_pos, distance_bp) or (None, None) if no TP on this chromosome.
    """
    arr = tp_index.get(chrom)
    if arr is None or len(arr) == 0:
        return None, None

    idx = np.searchsorted(arr, pos)

    candidates = []
    for i in [idx - 1, idx, idx + 1]:
        if 0 <= i < len(arr):
            tp_pos = int(arr[i])
            if is_tp and tp_pos == pos:
                continue  # Skip self
            candidates.append((tp_pos, abs(tp_pos - pos)))

    if not candidates:
        return None, None

    candidates.sort(key=lambda x: x[1])
    return candidates[0]


# ---------------------------------------------------------------------------
# Summary generation
# ---------------------------------------------------------------------------

def write_site_summary(
    site_dir: Path,
    row: pd.Series,
    metadata: Dict,
    significance: Dict,
    reads_df: Optional[pd.DataFrame],
    clustering_metrics: Optional[Dict] = None,
    nearest_tp_info: Optional[Tuple[Optional[int], Optional[int]]] = None,
):
    """Write summary.md for a single site."""
    chrom = row["Chr"]
    pos = int(row["Pos"])
    cat = row["category"]
    label = row["label"]
    truth = row["truth_label"]

    lines = [
        f"# {truth} Site: {chrom}:{pos}",
        f"",
        f"**Category**: {cat}",
        f"**Label**: {label}",
        f"",
    ]

    # Nearest TP context
    if nearest_tp_info is not None:
        ntp_pos, ntp_dist = nearest_tp_info
        if ntp_pos is not None:
            if ntp_dist < 1000:
                dist_str = f"{ntp_dist:,} bp"
            elif ntp_dist < 1_000_000:
                dist_str = f"{ntp_dist / 1000:.1f} kb"
            else:
                dist_str = f"{ntp_dist / 1_000_000:.2f} Mb"
            lines.append(f"**Nearest TP**: {chrom}:{ntp_pos:,} ({dist_str} away)")
        else:
            lines.append(f"**Nearest TP**: N/A (no TP on {chrom})")
        lines.append(f"")

    lines.extend([
        f"## Basic Info",
        f"",
    ])

    # Metadata
    for k, v in metadata.items():
        lines.append(f"- **{k}**: {v}")

    # Key features from master dataset
    lines.extend([
        f"",
        f"## Key Features (from master dataset)",
        f"",
        f"| Feature | Value |",
        f"|---------|-------|",
    ])
    key_cols = [
        "NumReads", "NumCpGs", "CramersV", "PairwiseMedianDist",
        "HPMergedDelta", "HPMergedSig", "AlleleDelta", "AlleleSig",
        "caller_af", "Coverage_Multiple", "Potential_LOH",
        "HP_Ratio", "Quality_Score", "Quality_Tier",
        "LabelHPPermanovaValid", "LabelAllelePermanovaValid",
    ]
    for col in key_cols:
        val = row.get(col, "N/A")
        lines.append(f"| {col} | {val} |")

    # Significance from JSON
    if significance:
        lines.extend([
            f"",
            f"## Statistical Analysis (significance.json)",
            f"",
        ])
        sig = significance
        lines.append(f"- **Passed Gating**: {sig.get('passed_gating', 'N/A')}")
        lines.append(f"- **Heuristic Score**: {sig.get('heuristic_score', 'N/A')}")

        if "global_alt" in sig:
            ga = sig["global_alt"]
            lines.append(f"- **Global Alt χ²**: p={ga.get('p_value', 'N/A')}, "
                         f"CramersV={ga.get('cramers_v', 'N/A')}")
        if "global_hp" in sig:
            gh = sig["global_hp"]
            lines.append(f"- **Global HP χ²**: p={gh.get('p_value', 'N/A')}, "
                         f"CramersV={gh.get('cramers_v', 'N/A')}")

        if "label_structure" in sig:
            ls = sig["label_structure"]
            lines.append(f"- **HP PERMANOVA**: valid={ls.get('hp_permanova_valid', 'N/A')}, "
                         f"F={ls.get('hp_permanova_f', 'N/A')}, p={ls.get('hp_permanova_p', 'N/A')}")
            lines.append(f"- **Allele PERMANOVA**: valid={ls.get('allele_permanova_valid', 'N/A')}, "
                         f"F={ls.get('allele_permanova_f', 'N/A')}, p={ls.get('allele_permanova_p', 'N/A')}")

    # Read structure summary
    if reads_df is not None:
        lines.extend([
            f"",
            f"## Read Structure",
            f"",
            f"- Total reads: {len(reads_df)}",
        ])
        hp_parsed = reads_df["hp"].astype(str).str.split("-").str[0]
        hp_parsed = pd.to_numeric(hp_parsed, errors="coerce").fillna(0).astype(int)
        hp_counts = hp_parsed.value_counts().to_dict()
        for hp_val in sorted(hp_counts.keys()):
            lines.append(f"- HP{hp_val}: {hp_counts[hp_val]}")
        if "alt_support" in reads_df.columns:
            alt_counts = reads_df["alt_support"].value_counts().to_dict()
            for allele in sorted(alt_counts.keys()):
                lines.append(f"- {allele}: {alt_counts[allele]}")

    # Clustering metrics
    if clustering_metrics:
        lines.extend([
            f"",
            f"## Clustering Structure Metrics",
            f"",
            f"| Metric | Value |",
            f"|--------|-------|",
            f"| Silhouette Score | {clustering_metrics.get('silhouette', 'N/A')} |",
            f"| Structure Verdict | **{clustering_metrics.get('structure_verdict', 'N/A')}** |",
            f"| Optimal K | {clustering_metrics.get('n_clusters', 'N/A')} |",
            f"| Cluster Sizes | {clustering_metrics.get('cluster_sizes', [])} |",
            f"| HP Cluster Purity | {clustering_metrics.get('hp_cluster_purity', 'N/A')} |",
            f"| Allele Cluster Purity | {clustering_metrics.get('allele_cluster_purity', 'N/A')} |",
        ])

    # Figure references
    lines.extend([
        f"",
        f"## Figures",
        f"",
        f"- [Methylation Heatmap](v1_methylation_heatmap.png)",
        f"- [Distance Matrix Heatmap](v2_distance_heatmap.png)",
        f"- [Read Structure](v3_read_structure.png)",
        f"- [Cluster vs Label](v4_label_vs_cluster.png)",
    ])

    summary_path = site_dir / "summary.md"
    summary_path.write_text("\n".join(lines))


def write_comparison_template(cat_dir: Path, cat_name: str, cat_desc: str):
    """Write comparison.md template for a category."""
    lines = [
        f"# {cat_desc} — TP vs FP 肉眼觀察",
        f"",
        f"## 甲基化模式差異",
        f"- TP 甲基化特徵: [觀察]",
        f"- FP 甲基化特徵: [觀察]",
        f"- 是否有肉眼可辨差異: [是/否/不確定]",
        f"",
        f"## HP / Allele 分佈差異",
        f"- TP 的 HP tag 分佈: [觀察]",
        f"- FP 的 HP tag 分佈: [觀察]",
        f"- HP 分配是否與 TP/FP 類別相關: [觀察]",
        f"",
        f"## 距離矩陣特徵",
        f"- TP 的 clustering 結構: [觀察]",
        f"- FP 的 clustering 結構: [觀察]",
        f"",
        f"## IGV 讀段特徵",
        f"- TP 在 IGV 中的特徵: [觀察]",
        f"- FP 在 IGV 中的特徵: [觀察]",
        f"",
        f"## 總結判斷",
        f"- 這個類別的 TP/FP 是否有人眼可辨的系統性差異: [判定]",
        f"- 如果有差異，差異的主要來源是: [甲基化/HP分佈/clustering/IGV]",
    ]
    (cat_dir / "comparison.md").write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# IGV Batch Script
# ---------------------------------------------------------------------------

def generate_igv_batch(sampled_df: pd.DataFrame):
    """Generate IGV batch script for automated screenshots."""
    igv_dir = ensure_dir(OUT_DIR / "igv_screenshots")

    lines = [
        "new",
        f"load {IGV_SESSION_TEMPLATE}",
        f"snapshotDirectory {igv_dir}",
        "maxPanelHeight 2000",
        "",
    ]

    for _, row in sampled_df.iterrows():
        chrom = row["Chr"]
        pos = int(row["Pos"])
        cat = row["category"]
        label = row["label"]
        site_id = row["site_id"]

        # Navigate to position with ±500bp window
        locus = f"{chrom}:{max(1, pos-500)}-{pos+500}"
        lines.append(f"goto {locus}")
        lines.append(f"snapshot {site_id}.png")
        lines.append("")

    batch_path = OUT_DIR / "igv_batch.txt"
    batch_path.write_text("\n".join(lines))
    print(f"  IGV batch script: {batch_path}")

    # Also write a locus list for manual inspection
    locus_lines = ["# Locus list for manual IGV inspection", ""]
    for _, row in sampled_df.iterrows():
        chrom = row["Chr"]
        pos = int(row["Pos"])
        locus_lines.append(
            f"{row['category']}\t{row['label']}\t{row['truth_label']}\t"
            f"{chrom}:{pos}\t"
            f"AlleleDelta={row.get('AlleleDelta', 'N/A')}\t"
            f"HPMergedSig={row.get('HPMergedSig', 'N/A')}"
        )
    (OUT_DIR / "igv_locus_list.tsv").write_text("\n".join(locus_lines))


# ---------------------------------------------------------------------------
# Index
# ---------------------------------------------------------------------------

def write_index(sampled_df: pd.DataFrame):
    """Write 00_inspection_index.md."""
    lines = [
        "# 肉眼檢視文件索引",
        "",
        f"**生成時間**: Wave 3 M6",
        f"**樣本**: HCC1395 TO mode",
        f"**總位點數**: {len(sampled_df)}",
        "",
        "## 分類",
        "",
    ]

    for cat_name, cat_info in CATEGORIES.items():
        cat_sites = sampled_df[sampled_df["category"] == cat_name]
        tp_count = (cat_sites["truth_label"] == "TP").sum()
        fp_count = (cat_sites["truth_label"] == "FP").sum()
        lines.append(f"### {cat_name} ({cat_info['desc']})")
        lines.append(f"")
        lines.append(f"TP: {tp_count}, FP: {fp_count}")
        lines.append(f"")
        lines.append(f"| Site | Truth | Chr:Pos | AlleleDelta | HPSig | NumReads |")
        lines.append(f"|------|-------|---------|-------------|-------|----------|")

        for _, row in cat_sites.iterrows():
            site_dir = f"{cat_name}/{row['site_id']}"
            lines.append(
                f"| [{row['site_id']}]({site_dir}/summary.md) | {row['truth_label']} | "
                f"{row['Chr']}:{int(row['Pos'])} | "
                f"{row.get('AlleleDelta', 'N/A')} | "
                f"{row.get('HPMergedSig', 'N/A')} | "
                f"{row.get('NumReads', 'N/A')} |"
            )
        lines.append(f"")
        lines.append(f"[TP vs FP 肉眼觀察]({cat_name}/comparison.md)")
        lines.append(f"")

    lines.extend([
        "## IGV",
        "",
        "- [IGV Batch Script](igv_batch.txt) — 使用 `xvfb-run /bip7_disk/IGV_Linux_2.11.1/igv.sh --batch igv_batch.txt`",
        "- [Locus List](igv_locus_list.tsv) — 供手動 IGV 檢視",
        "",
        "## 抽樣條件",
        "",
        "- NumReads ∈ [30, 100]",
        "- NumCpGs ≥ 10",
        "- 排除著絲點（centromere）與端粒（telomere）區域",
        "- 每類 TP/FP 各 5 個（如可用）",
    ])

    (OUT_DIR / "00_inspection_index.md").write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Wave 3 M6: Visual Inspection Report")
    print("=" * 70)

    # Load and sample
    df = load_master_dataset()
    sampled_df = sample_sites(df)

    # Build TP position index for nearest-TP lookup
    tp_index = build_tp_position_index(df)
    print(f"\n  TP index: {sum(len(v) for v in tp_index.values()):,} TPs across {len(tp_index)} chromosomes")

    # Generate index and IGV batch
    write_index(sampled_df)
    generate_igv_batch(sampled_df)

    # Process each site
    n_found = 0
    n_missing = 0

    for _, row in sampled_df.iterrows():
        cat_name = row["category"]
        site_id = row["site_id"]
        chrom = str(row["Chr"])
        pos = int(row["Pos"])
        is_tp = row["truth_label"] == "TP"

        site_dir = ensure_dir(OUT_DIR / cat_name / site_id)

        # Nearest TP lookup
        nearest_tp_info = find_nearest_tp(chrom, pos, is_tp, tp_index)

        # Find ISM region directory
        region_dir = find_region_dir(chrom, pos, is_tp)

        if region_dir is None:
            print(f"  MISSING: {site_id} ({chrom}:{pos}, {'TP' if is_tp else 'FP'})")
            n_missing += 1
            # Write minimal summary
            write_site_summary(site_dir, row, {}, {}, None,
                               nearest_tp_info=nearest_tp_info)
            continue

        n_found += 1
        print(f"  Processing: {site_id} ({chrom}:{pos})")

        # Read ISM data
        metadata = read_metadata(region_dir)
        reads_df = read_reads_tsv(region_dir)
        meth_df, cpg_positions = read_methylation_matrix(region_dir)
        significance = read_significance(region_dir)
        dist_matrix = read_distance_matrix(region_dir)

        # Compute clustering metrics (quantitative structure detection)
        clustering_metrics = compute_clustering_metrics(dist_matrix, reads_df)

        # Write summary
        write_site_summary(site_dir, row, metadata, significance, reads_df,
                           clustering_metrics, nearest_tp_info=nearest_tp_info)

        # V1: Methylation heatmap (ISM standard: sns.clustermap + dendrogram)
        if meth_df is not None and cpg_positions is not None:
            plot_methylation_heatmap(
                meth_df, cpg_positions, reads_df, dist_matrix,
                f"{row['truth_label']} — {chrom}:{pos}\n"
                f"AlleleDelta={row.get('AlleleDelta', 'N/A')}, "
                f"HPSig={row.get('HPMergedSig', 'N/A')}",
                site_dir / "v1_methylation_heatmap.png",
            )

        # V2: Distance heatmap (ISM standard: sns.clustermap + dual dendrogram)
        plot_distance_heatmap(
            dist_matrix, reads_df,
            f"{row['truth_label']} — {chrom}:{pos} (Bernoulli Distance)",
            site_dir / "v2_distance_heatmap.png",
        )

        # V3: Read structure (ISM standard colors)
        if reads_df is not None:
            plot_read_structure(
                reads_df,
                f"{row['truth_label']} — {chrom}:{pos}",
                site_dir / "v3_read_structure.png",
            )

        # V4: Cluster vs Label (ISM standard colors)
        plot_label_vs_cluster(
            dist_matrix, reads_df, clustering_metrics,
            f"{row['truth_label']} — {chrom}:{pos}",
            site_dir / "v4_label_vs_cluster.png",
        )

    # Write comparison templates
    for cat_name, cat_info in CATEGORIES.items():
        cat_dir = ensure_dir(OUT_DIR / cat_name)
        write_comparison_template(cat_dir, cat_name, cat_info["desc"])

    print(f"\n--- Summary ---")
    print(f"  Sites found: {n_found}")
    print(f"  Sites missing: {n_missing}")
    print(f"  Output: {OUT_DIR}")
    print(f"\nM6 complete.")


if __name__ == "__main__":
    main()
