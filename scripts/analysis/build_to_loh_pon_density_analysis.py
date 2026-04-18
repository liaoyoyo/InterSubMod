#!/usr/bin/env python3
"""
PON Hit Count + Neighbor Variant Density Analysis for TO regions.

Part 1: Check each TO variant against 4 PON databases (gnomAD, dbSNP, 1000G, CoLoRSdb)
        and analyze PON hit count vs TP/FP status.
Part 2: Compute neighbor variant density at 3 window sizes (5kb, 50kb, 500kb)
        and analyze density vs TP/FP status.

Output: output/to_loh_dual_definition_report/supplementary/
"""

import gzip
import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from sklearn.metrics import roc_auc_score

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
MASTER_TSV = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260327_loh_round1_cross_sample_audit/"
    "all_region_rows.tsv.gz"
)
PON_DIR = "/big8_disk/data/PON/clairs-to_databases"
PON_DEFS = {
    "gnomAD": {
        "file": f"{PON_DIR}/gnomad.r2.1.af-ge-0.001.sites.vcf.gz",
        "match": "allele",   # position + ALT
    },
    "dbSNP": {
        "file": f"{PON_DIR}/dbsnp.b138.non-somatic.sites.vcf.gz",
        "match": "allele",
    },
    "1000G": {
        "file": f"{PON_DIR}/1000g-pon.sites.vcf.gz",
        "match": "position",
    },
    "CoLoRSdb": {
        "file": f"{PON_DIR}/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz",
        "match": "position",
    },
}
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/output/"
    "to_loh_dual_definition_report/supplementary"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Style
TP_COLOR = "#2196F3"
FP_COLOR = "#F44336"
plt.rcParams.update({
    "font.size": 11,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
})


# ---------------------------------------------------------------------------
# Load master dataset (TO only)
# ---------------------------------------------------------------------------
def load_master():
    print("Loading master dataset...")
    t0 = time.time()
    df = pd.read_csv(MASTER_TSV, sep="\t", compression="gzip")
    # Filter to tumor-only
    df = df[df["mode"] == "to"].copy()
    print(f"  Loaded {len(df):,} TO regions in {time.time()-t0:.1f}s")
    print(f"  TP={sum(df['truth_label']=='TP'):,}  FP={sum(df['truth_label']=='FP'):,}")
    return df


# ---------------------------------------------------------------------------
# Part 1: PON Hit Count
# ---------------------------------------------------------------------------
def query_pon_tabix(pon_name, pon_cfg, unique_variants):
    """
    Query a single PON VCF using tabix -R (batch region file).
    Creates a temporary BED file with all query positions, then
    runs tabix -R to fetch all matching records at once.

    Returns a set of variant_keys that hit.
    """
    import subprocess
    import tempfile

    vcf_path = pon_cfg["file"]
    match_mode = pon_cfg["match"]
    print(f"  Querying {pon_name} ({match_mode} match)...")
    t0 = time.time()

    hits = set()

    # Build lookup structures
    if match_mode == "position":
        # (chrom, pos) -> list of variant_keys
        pos_to_keys = defaultdict(list)
        for chrom, pos, ref, alt, vkey in unique_variants:
            pos_to_keys[(chrom, pos)].append(vkey)
    else:
        # (chrom, pos, alt) -> list of variant_keys
        allele_to_keys = defaultdict(list)
        pos_set = set()
        for chrom, pos, ref, alt, vkey in unique_variants:
            allele_to_keys[(chrom, pos, alt)].append(vkey)
            pos_set.add((chrom, pos))

    # Create temp BED file with all positions (0-based start, 1-based end)
    unique_positions = set()
    for chrom, pos, ref, alt, vkey in unique_variants:
        unique_positions.add((chrom, pos))

    # Sort positions for BED file
    def chrom_sort_key(cp):
        c = cp[0].replace("chr", "")
        try:
            return (int(c), cp[1])
        except ValueError:
            return (999, cp[1])

    sorted_positions = sorted(unique_positions, key=chrom_sort_key)

    print(f"    Creating region file for {len(sorted_positions):,} unique positions...")

    # Write as BED file (0-based start, 1-based end, tab-delimited)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
        region_file = f.name
        for chrom, pos in sorted_positions:
            f.write(f"{chrom}\t{pos - 1}\t{pos}\n")

    print(f"    Running tabix -R ...")
    try:
        result = subprocess.run(
            ["tabix", "-R", region_file, vcf_path],
            capture_output=True, text=True, timeout=600
        )
    except subprocess.TimeoutExpired:
        print(f"    WARNING: tabix timed out for {pon_name}")
        import os
        os.unlink(region_file)
        return hits
    finally:
        import os
        if os.path.exists(region_file):
            os.unlink(region_file)

    if result.returncode != 0:
        print(f"    WARNING: tabix returned error for {pon_name}: {result.stderr[:200]}")
        return hits

    # Parse tabix output
    n_records = 0
    n_found = 0
    for line in result.stdout.split("\n"):
        if not line:
            continue
        n_records += 1
        parts = line.split("\t", 5)
        if len(parts) < 5:
            continue

        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            continue

        if match_mode == "position":
            key = (chrom, pos)
            if key in pos_to_keys:
                for vkey in pos_to_keys[key]:
                    hits.add(vkey)
                    n_found += 1
        else:
            # Allele match: check each ALT
            alts = parts[4].split(",")
            for alt in alts:
                key = (chrom, pos, alt)
                if key in allele_to_keys:
                    for vkey in allele_to_keys[key]:
                        hits.add(vkey)
                        n_found += 1

    elapsed = time.time() - t0
    print(f"    {pon_name} DONE: {len(hits):,} unique variant hits "
          f"({n_found:,} total matches from {n_records:,} PON records) "
          f"in {elapsed:.1f}s")
    return hits


def run_pon_analysis(df):
    """
    Run PON hit count analysis.
    Uses two matching strategies:
      - 'strict': allele match for gnomAD/dbSNP, position-only for 1000G/CoLoRSdb
      - 'position_only': position-only match for ALL 4 PON databases
    """
    print("\n=== Part 1: PON Hit Count Analysis ===")

    # Deduplicate by variant_key
    unique_rows = df.drop_duplicates(subset="variant_key")[
        ["Chr", "Pos", "Ref", "Alt", "variant_key"]
    ].copy()
    unique_variants = list(zip(
        unique_rows["Chr"], unique_rows["Pos"],
        unique_rows["Ref"], unique_rows["Alt"],
        unique_rows["variant_key"]
    ))
    print(f"Unique variants to query: {len(unique_variants):,}")

    # ----- Strategy 1: standard matching (allele for gnomAD/dbSNP) -----
    print("\n--- Strategy: standard matching ---")
    pon_hits_strict = {}
    for pon_name, pon_cfg in PON_DEFS.items():
        hits = query_pon_tabix(pon_name, pon_cfg, unique_variants)
        pon_hits_strict[pon_name] = hits

    # ----- Strategy 2: position-only for ALL databases -----
    print("\n--- Strategy: position-only matching for ALL databases ---")
    pon_hits_posonly = {}
    for pon_name, pon_cfg in PON_DEFS.items():
        cfg_posonly = {"file": pon_cfg["file"], "match": "position"}
        hits = query_pon_tabix(f"{pon_name}_posonly", cfg_posonly, unique_variants)
        pon_hits_posonly[pon_name] = hits

    # Annotate with BOTH strategies
    pon_names = list(PON_DEFS.keys())

    # Strict matching columns
    for pn in pon_names:
        df[f"pon_{pn}"] = df["variant_key"].isin(pon_hits_strict[pn]).astype(int)
    df["pon_hit_count"] = sum(df[f"pon_{pn}"] for pn in pon_names)

    # Position-only columns
    for pn in pon_names:
        df[f"pon_{pn}_posonly"] = df["variant_key"].isin(pon_hits_posonly[pn]).astype(int)
    df["pon_hit_count_posonly"] = sum(df[f"pon_{pn}_posonly"] for pn in pon_names)

    return df, pon_names


def safe_auc(y_true, y_score, label_positive="FP"):
    """Compute AUC safely, returning NaN if not computable."""
    y_bin = (y_true == label_positive).astype(int)
    if y_bin.sum() == 0 or y_bin.sum() == len(y_bin):
        return np.nan
    try:
        return roc_auc_score(y_bin, y_score)
    except Exception:
        return np.nan


def plot_pon_results(df, pon_names):
    """Generate PON hit analysis figure with both matching strategies."""
    print("\nGenerating PON hit figure...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle("PON Hit Count Analysis — TO PASS Variants (Post-PON-Filter)\n"
                 "Note: ClairS-TO already applied these 4 PON databases as filters",
                 fontsize=13, fontweight="bold")

    tp_mask = df["truth_label"] == "TP"
    fp_mask = df["truth_label"] == "FP"

    # --- Panel [0,0]: Position-only hit count distribution ---
    ax = axes[0, 0]
    hit_col = "pon_hit_count_posonly"
    hit_counts = sorted(df[hit_col].unique())
    tp_counts = []
    fp_counts = []
    for hc in hit_counts:
        sub = df[df[hit_col] == hc]
        tp_counts.append(sum(sub["truth_label"] == "TP"))
        fp_counts.append(sum(sub["truth_label"] == "FP"))

    x_pos = np.arange(len(hit_counts))
    ax.bar(x_pos, tp_counts, color=TP_COLOR, label="TP", alpha=0.85)
    ax.bar(x_pos, fp_counts, bottom=tp_counts, color=FP_COLOR, label="FP", alpha=0.85)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([str(hc) for hc in hit_counts])
    ax.set_xlabel("PON Hit Count (Position-Only Matching)")
    ax.set_ylabel("Number of Variants")
    ax.set_title("Position-Only PON Hits: All TO")
    ax.legend()
    for i, hc in enumerate(hit_counts):
        total = tp_counts[i] + fp_counts[i]
        if total > 0:
            fp_rate = fp_counts[i] / total
            ax.text(i, tp_counts[i] + fp_counts[i], f"FP={fp_rate:.1%}\nn={total:,}",
                    ha="center", va="bottom", fontsize=8, fontweight="bold")

    # --- Panel [0,1]: Per-PON position-only hit rate ---
    ax = axes[0, 1]
    pon_labels = pon_names
    tp_rates = []
    fp_rates = []
    for pn in pon_labels:
        col = f"pon_{pn}_posonly"
        tp_rates.append(df.loc[tp_mask, col].mean())
        fp_rates.append(df.loc[fp_mask, col].mean())

    x_pon = np.arange(len(pon_labels))
    w = 0.35
    ax.bar(x_pon - w / 2, tp_rates, w, color=TP_COLOR, label="TP", alpha=0.85)
    ax.bar(x_pon + w / 2, fp_rates, w, color=FP_COLOR, label="FP", alpha=0.85)
    ax.set_xticks(x_pon)
    ax.set_xticklabels(pon_labels, rotation=15, ha="right")
    ax.set_ylabel("Position-Only Hit Rate")
    ax.set_title("Per-PON Position Hit Rate: TP vs FP")
    ax.legend()
    for i, pn in enumerate(pon_labels):
        ax.text(i - w / 2, tp_rates[i] + 0.001, f"{tp_rates[i]:.4f}",
                ha="center", va="bottom", fontsize=8)
        ax.text(i + w / 2, fp_rates[i] + 0.001, f"{fp_rates[i]:.4f}",
                ha="center", va="bottom", fontsize=8)

    # --- Panel [1,0]: Position-only hits: LOH vs Non-LOH ---
    ax = axes[1, 0]
    loh_col = "to_loh_bed_hit"
    categories = [
        ("All", df),
        ("LOH.bed=T", df[df[loh_col] == True]),
        ("LOH.bed=F", df[df[loh_col] == False]),
    ]
    cat_labels = []
    cat_fp_rates_hit = []
    cat_fp_rates_nohit = []
    cat_n_hit = []
    cat_n_nohit = []
    for label, sub in categories:
        has_hit = sub[sub["pon_hit_count_posonly"] > 0]
        no_hit = sub[sub["pon_hit_count_posonly"] == 0]
        fp_hit = sum(has_hit["truth_label"] == "FP") / max(len(has_hit), 1)
        fp_nohit = sum(no_hit["truth_label"] == "FP") / max(len(no_hit), 1)
        cat_labels.append(label)
        cat_fp_rates_hit.append(fp_hit)
        cat_fp_rates_nohit.append(fp_nohit)
        cat_n_hit.append(len(has_hit))
        cat_n_nohit.append(len(no_hit))

    x_cat = np.arange(len(cat_labels))
    w = 0.35
    ax.bar(x_cat - w / 2, cat_fp_rates_nohit, w, color="#607D8B", label="No PON hit", alpha=0.85)
    ax.bar(x_cat + w / 2, cat_fp_rates_hit, w, color="#FF5722", label="Has PON hit", alpha=0.85)
    ax.set_xticks(x_cat)
    ax.set_xticklabels(cat_labels)
    ax.set_ylabel("FP Rate")
    ax.set_title("FP Rate: PON Hit vs No Hit (Position-Only)")
    ax.legend()
    for i in range(len(cat_labels)):
        ax.text(i - w / 2, cat_fp_rates_nohit[i] + 0.005,
                f"n={cat_n_nohit[i]:,}", ha="center", va="bottom", fontsize=7)
        ax.text(i + w / 2, cat_fp_rates_hit[i] + 0.005,
                f"n={cat_n_hit[i]:,}", ha="center", va="bottom", fontsize=7)

    # --- Panel [1,1]: Summary text ---
    ax = axes[1, 1]
    ax.axis("off")

    # Compute summary stats
    n_total = len(df)
    n_strict_hit = sum(df["pon_hit_count"] > 0)
    n_posonly_hit = sum(df["pon_hit_count_posonly"] > 0)

    # Position-only AUC
    posonly_auc = safe_auc(df["truth_label"], df["pon_hit_count_posonly"])

    lines = [
        "Key Finding",
        "",
        f"Total TO PASS: {n_total:,}",
        f"Unique positions: {df['variant_key'].nunique():,}",
        "",
        f"Strict match (allele): {n_strict_hit:,} ({n_strict_hit/n_total*100:.3f}%)",
        f"Position-only hits:    {n_posonly_hit:,} ({n_posonly_hit/n_total*100:.2f}%)",
        f"Position-only AUC:     {posonly_auc:.4f}",
        "",
        "Explanation:",
        "ClairS-TO already applied these 4 PON",
        "databases as filters. PASS variants",
        "survived PON filtering => allele match ~0%.",
        "",
        "Position-only overlap: different variant",
        "exists at the same genomic position in PON.",
        f"  gnomAD: {sum(df['pon_gnomAD_posonly']):,} hits",
        f"  dbSNP:  {sum(df['pon_dbSNP_posonly']):,} hits",
        f"  1000G:  {sum(df['pon_1000G_posonly']):,} hits",
        f"  CoLoRSdb: {sum(df['pon_CoLoRSdb_posonly']):,} hits",
        "",
        "Position-only hit FP rate: "
        f"{sum((df['pon_hit_count_posonly']>0) & (df['truth_label']=='FP')) / max(sum(df['pon_hit_count_posonly']>0),1):.1%}",
        f"vs overall FP rate: {sum(df['truth_label']=='FP')/len(df):.1%}",
    ]
    summary_text = "\n".join(lines)
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=9, verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#f5f5f5", alpha=0.8))

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    out_path = OUT_DIR / "sup_pon_hit_distribution.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def save_pon_tsv(df):
    """Save PON hit results TSV with both matching strategies."""
    rows = []

    # Section 1: Position-only hit count distribution
    for hc in sorted(df["pon_hit_count_posonly"].unique()):
        sub = df[df["pon_hit_count_posonly"] == hc]
        n_tp = sum(sub["truth_label"] == "TP")
        n_fp = sum(sub["truth_label"] == "FP")
        total = n_tp + n_fp
        fp_rate = n_fp / total if total > 0 else np.nan
        ad_auc = safe_auc(sub["truth_label"], sub["allele_delta"])
        af_auc = safe_auc(sub["truth_label"], sub["caller_af"])
        rows.append({
            "match_type": "position_only",
            "PON_hits": hc,
            "N_TP": n_tp,
            "N_FP": n_fp,
            "FP_rate": round(fp_rate, 4),
            "AlleleDelta_AUC": round(ad_auc, 4) if not np.isnan(ad_auc) else "NA",
            "caller_af_AUC": round(af_auc, 4) if not np.isnan(af_auc) else "NA",
        })

    # Section 2: Strict (allele) hit count distribution
    for hc in sorted(df["pon_hit_count"].unique()):
        sub = df[df["pon_hit_count"] == hc]
        n_tp = sum(sub["truth_label"] == "TP")
        n_fp = sum(sub["truth_label"] == "FP")
        total = n_tp + n_fp
        fp_rate = n_fp / total if total > 0 else np.nan
        ad_auc = safe_auc(sub["truth_label"], sub["allele_delta"])
        af_auc = safe_auc(sub["truth_label"], sub["caller_af"])
        rows.append({
            "match_type": "strict_allele",
            "PON_hits": hc,
            "N_TP": n_tp,
            "N_FP": n_fp,
            "FP_rate": round(fp_rate, 4),
            "AlleleDelta_AUC": round(ad_auc, 4) if not np.isnan(ad_auc) else "NA",
            "caller_af_AUC": round(af_auc, 4) if not np.isnan(af_auc) else "NA",
        })

    out = pd.DataFrame(rows)
    path = OUT_DIR / "sup_pon_hit_results.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"  Saved: {path}")
    print(out.to_string(index=False))

    # Also save per-PON breakdown
    pon_names = [pn for pn in ["gnomAD", "dbSNP", "1000G", "CoLoRSdb"]]
    per_pon_rows = []
    tp_mask = df["truth_label"] == "TP"
    fp_mask = df["truth_label"] == "FP"
    for pn in pon_names:
        for mtype, suffix in [("strict", ""), ("position_only", "_posonly")]:
            col = f"pon_{pn}{suffix}"
            n_hits_tp = df.loc[tp_mask, col].sum()
            n_hits_fp = df.loc[fp_mask, col].sum()
            per_pon_rows.append({
                "PON": pn,
                "match_type": mtype,
                "TP_hits": int(n_hits_tp),
                "FP_hits": int(n_hits_fp),
                "TP_rate": round(n_hits_tp / tp_mask.sum(), 6),
                "FP_rate": round(n_hits_fp / fp_mask.sum(), 6),
            })
    per_pon = pd.DataFrame(per_pon_rows)
    path2 = OUT_DIR / "sup_pon_per_database_breakdown.tsv"
    per_pon.to_csv(path2, sep="\t", index=False)
    print(f"  Saved: {path2}")
    print(per_pon.to_string(index=False))


# ---------------------------------------------------------------------------
# Part 2: Neighbor Variant Density
# ---------------------------------------------------------------------------
def compute_neighbor_density(df):
    """
    For each variant, count neighbors within +-5kb, +-50kb, +-500kb
    from the same sample. Uses sorted array + bisect for efficiency.
    """
    print("\n=== Part 2: Neighbor Variant Density ===")
    import bisect

    windows = [5000, 50000, 500000]
    window_labels = ["5kb", "50kb", "500kb"]

    # Group by (sample, chrom) and sort by position
    grouped = df.groupby(["sample", "Chr"])

    # Pre-build sorted position arrays
    print("  Building position indices...")
    pos_arrays = {}
    for (sample, chrom), grp in grouped:
        positions = np.sort(grp["Pos"].values)
        pos_arrays[(sample, chrom)] = positions

    # For each variant, count neighbors
    print("  Computing neighbor counts...")
    for w, wl in zip(windows, window_labels):
        col_name = f"density_{wl}"
        df[col_name] = 0  # Initialize

    t0 = time.time()
    n_total = len(df)
    results = {wl: np.zeros(n_total, dtype=int) for wl in window_labels}

    for idx_block, ((sample, chrom), grp) in enumerate(grouped):
        positions = pos_arrays[(sample, chrom)]
        # Get the dataframe indices for this group
        grp_indices = grp.index.values
        grp_positions = grp["Pos"].values

        for w, wl in zip(windows, window_labels):
            for i, (df_idx, pos) in enumerate(zip(grp_indices, grp_positions)):
                lo = bisect.bisect_left(positions, pos - w)
                hi = bisect.bisect_right(positions, pos + w)
                # Subtract 1 for the variant itself
                count = hi - lo - 1
                results[wl][df.index.get_loc(df_idx)] = count

        if (idx_block + 1) % 50 == 0:
            elapsed = time.time() - t0
            print(f"    Processed {idx_block+1} (sample, chrom) groups ({elapsed:.0f}s)")

    for wl in window_labels:
        df[f"density_{wl}"] = results[wl]

    elapsed = time.time() - t0
    print(f"  Density computation completed in {elapsed:.1f}s")
    return df, window_labels


def plot_density_results(df, window_labels):
    """Generate density analysis figure."""
    print("\nGenerating density figure...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle("Neighbor Variant Density Analysis — TO Regions",
                 fontsize=14, fontweight="bold")

    tp_mask = df["truth_label"] == "TP"
    fp_mask = df["truth_label"] == "FP"

    # --- Panel [0,0]: 5kb density histogram ---
    ax = axes[0, 0]
    clip_5k = 20
    tp_vals = np.clip(df.loc[tp_mask, "density_5kb"].values, 0, clip_5k)
    fp_vals = np.clip(df.loc[fp_mask, "density_5kb"].values, 0, clip_5k)
    bins = np.arange(-0.5, clip_5k + 1.5, 1)
    ax.hist(tp_vals, bins=bins, alpha=0.6, color=TP_COLOR, label=f"TP (n={len(tp_vals):,})",
            density=True)
    ax.hist(fp_vals, bins=bins, alpha=0.6, color=FP_COLOR, label=f"FP (n={len(fp_vals):,})",
            density=True)
    ax.set_xlabel("Neighbors within ±5kb")
    ax.set_ylabel("Density")
    ax.set_title("5kb Neighbor Count: TP vs FP")
    ax.legend(fontsize=9)

    # --- Panel [0,1]: 50kb density histogram ---
    ax = axes[0, 1]
    clip_50k = 100
    tp_vals_50 = np.clip(df.loc[tp_mask, "density_50kb"].values, 0, clip_50k)
    fp_vals_50 = np.clip(df.loc[fp_mask, "density_50kb"].values, 0, clip_50k)
    bins_50 = np.linspace(-0.5, clip_50k + 0.5, 51)
    ax.hist(tp_vals_50, bins=bins_50, alpha=0.6, color=TP_COLOR,
            label=f"TP (median={np.median(tp_vals_50):.1f})", density=True)
    ax.hist(fp_vals_50, bins=bins_50, alpha=0.6, color=FP_COLOR,
            label=f"FP (median={np.median(fp_vals_50):.1f})", density=True)
    ax.set_xlabel("Neighbors within ±50kb")
    ax.set_ylabel("Density")
    ax.set_title("50kb Neighbor Count: TP vs FP")
    ax.legend(fontsize=9)

    # --- Panel [1,0]: 5kb density vs FP rate ---
    ax = axes[1, 0]
    density_col = "density_5kb"
    # Create bins (0, 1, 2, 3, 4, 5-9, 10-19, 20+)
    bin_edges = [0, 1, 2, 3, 4, 5, 10, 20, 1000]
    bin_labels = ["0", "1", "2", "3", "4", "5-9", "10-19", "20+"]
    df["_dbin"] = pd.cut(df[density_col], bins=[-0.5] + [e - 0.5 for e in bin_edges[1:]],
                         labels=bin_labels, right=True)

    fp_rates = []
    ns = []
    valid_labels = []
    for bl in bin_labels:
        sub = df[df["_dbin"] == bl]
        if len(sub) > 0:
            fp_rates.append(sum(sub["truth_label"] == "FP") / len(sub))
            ns.append(len(sub))
            valid_labels.append(bl)

    x_d = np.arange(len(valid_labels))
    bars = ax.bar(x_d, fp_rates, color="#607D8B", alpha=0.8)
    ax.set_xticks(x_d)
    ax.set_xticklabels(valid_labels, rotation=0)
    ax.set_xlabel("5kb Neighbor Count Bin")
    ax.set_ylabel("FP Rate")
    ax.set_title("FP Rate by 5kb Neighbor Density")
    ax.axhline(sum(df["truth_label"] == "FP") / len(df), ls="--",
               color="red", alpha=0.6, label="Overall FP rate")
    ax.legend(fontsize=9)
    for i, (fr, n) in enumerate(zip(fp_rates, ns)):
        ax.text(i, fr + 0.005, f"n={n:,}", ha="center", va="bottom", fontsize=7, color="gray")

    df.drop(columns=["_dbin"], inplace=True)

    # --- Panel [1,1]: AUC by window size ---
    ax = axes[1, 1]
    aucs = []
    for wl in window_labels:
        auc_val = safe_auc(df["truth_label"], df[f"density_{wl}"])
        aucs.append(auc_val)

    x_w = np.arange(len(window_labels))
    bars = ax.bar(x_w, aucs, color=["#4CAF50", "#FF9800", "#9C27B0"], alpha=0.8)
    ax.set_xticks(x_w)
    ax.set_xticklabels([f"±{wl}" for wl in window_labels])
    ax.set_xlabel("Window Size")
    ax.set_ylabel("AUC (FP as positive)")
    ax.set_title("Density as FP Discriminator")
    ax.set_ylim(0.3, 0.8)
    ax.axhline(0.5, ls="--", color="gray", alpha=0.5, label="Random")
    ax.legend()
    for i, auc_val in enumerate(aucs):
        ax.text(i, auc_val + 0.01, f"{auc_val:.3f}", ha="center", fontsize=10,
                fontweight="bold")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = OUT_DIR / "sup_neighbor_density.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def save_density_tsv(df, window_labels):
    """Save density results TSV."""
    rows = []
    for wl in window_labels:
        col = f"density_{wl}"
        tp_mask = df["truth_label"] == "TP"
        fp_mask = df["truth_label"] == "FP"

        tp_med = df.loc[tp_mask, col].median()
        fp_med = df.loc[fp_mask, col].median()
        tp_mean = df.loc[tp_mask, col].mean()
        fp_mean = df.loc[fp_mask, col].mean()
        auc_val = safe_auc(df["truth_label"], df[col])

        rows.append({
            "window": wl,
            "TP_median": round(tp_med, 2),
            "FP_median": round(fp_med, 2),
            "TP_mean": round(tp_mean, 2),
            "FP_mean": round(fp_mean, 2),
            "AUC_FP_positive": round(auc_val, 4) if not np.isnan(auc_val) else "NA",
        })
    out = pd.DataFrame(rows)
    path = OUT_DIR / "sup_neighbor_density_results.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"  Saved: {path}")
    print(out.to_string(index=False))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t_start = time.time()

    df = load_master()

    # Part 1: PON
    df, pon_names = run_pon_analysis(df)
    plot_pon_results(df, pon_names)
    save_pon_tsv(df)

    # Part 2: Density
    df, window_labels = compute_neighbor_density(df)
    plot_density_results(df, window_labels)
    save_density_tsv(df, window_labels)

    total_time = time.time() - t_start
    print(f"\n=== All done in {total_time:.0f}s ===")


if __name__ == "__main__":
    main()
