#!/usr/bin/env python3
"""V5 Audit Suite — read intersection concordance analysis.

For each of 15 audit sites, fetch reads from 5 BAMs (Baseline, V2b, V3-Fixed,
V5, Paired tumor), restrict to primary-mapped reads, and compute the read_name
intersection of {BL, V5, Paired}. Within that shared set, compute four
concordance metrics that compare BL vs Paired and V5 vs Paired:

  L1  HP family agreement      — {1, 11} vs {2, 21} family consistency
  L2  HP exact agreement       — exact HP value (after canonical mapping)
  L3  HP1:HP2 ratio distance   — |ratio_X - ratio_PA| with orientation flip
  L4  Orientation-corrected    — per-PS best of (native, flipped) family match

This intersection-first design eliminates the read-set bias that confounds the
aggregate HP composition metric (different BAMs tag different read pools, so a
naive "% HP1" comparison is apples-to-oranges).

Outputs
-------
- TSV: per_site_concordance.tsv  (15 rows × all 4 metrics, both BL & V5)
- TSV: hp_family_exact.tsv       (long-format read-level breakdown)
- PNG: fig02a_4metric_heatmap.png
- PNG: fig02b_winloss_summary.png

Usage
-----
    python3 v5_read_intersection_concordance.py
"""
from __future__ import annotations

import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

plt.rcParams["font.family"] = "DejaVu Sans"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BAMS = {
    "BL": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
    "V2b": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_tagged.bam",
    "V3F": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam",
    "V5": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam",
    "PA": "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam",
}

SITES = [
    ("TP01", "chr6", 145444893),
    ("TP02", "chr4", 70548355),
    ("TP03", "chr5", 153209947),
    ("TP04", "chr16", 35118902),
    ("TP05", "chr7", 109185781),
    ("FPA1", "chr8", 93565727),
    ("FPA2", "chr9", 137953060),
    ("FPB1", "chr7", 52087777),
    ("FPB2", "chr9", 75383880),
    ("V5max1", "chr19", 4639528),
    ("V5max2", "chr19", 2235521),
    ("V5max3", "chr19", 7405500),
    ("SP1", "chr19", 17565944),
    ("SP2", "chr19", 12452332),
    ("SP3", "chr19", 12467180),
]

OUT_BASE = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite"
)
OUT_FIG_DIR = OUT_BASE / "figures" / "02_concordance"
OUT_DATA_DIR = OUT_BASE / "data"
OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DATA_DIR.mkdir(parents=True, exist_ok=True)

PAD = 50  # bp padding around target position when fetching reads

# ---------------------------------------------------------------------------
# HP normalization
# ---------------------------------------------------------------------------
# longphase / V5 BAM HP tag (integer):
#   1, 2          -> single-allele family marker (PS-implicit)
#   11, 21        -> 1-1, 2-1 in two-PS notation (HP family 1 / 2 respectively)
#   33            -> "Both" (Self-phasing 33 marker — no family)
#   None          -> untagged
# Paired BAM HP tag (string):
#   '1', '2'      -> family 1, family 2
#   '1-1','2-1'   -> family 1 (PS_id-allele); '1-2'/'2-2' similarly
#   '3'           -> "Both" / unphasable
#   None          -> untagged
LP_INT_TO_FAMILY = {1: 1, 2: 2, 11: 1, 21: 2, 33: 0}
LP_INT_TO_EXACT = {1: "1", 2: "2", 11: "1-1", 21: "2-1", 33: "3"}


def normalize_hp(raw):
    """Return (family, exact_canonical, has_tag).

    family : int    1 / 2 / 0(both/3) / -1(none)
    exact  : str    canonical exact HP string ("1", "2", "1-1", ..., "3", "NA")
    """
    if raw is None:
        return -1, "NA", False
    if isinstance(raw, int):
        family = LP_INT_TO_FAMILY.get(raw, -1)
        exact = LP_INT_TO_EXACT.get(raw, "NA")
        return family, exact, family != -1
    # string (paired)
    s = str(raw).strip()
    if s in ("", ".", "NA"):
        return -1, "NA", False
    if s == "3":
        return 0, "3", True
    # parse leading digit as family
    head = s.split("-")[0]
    try:
        f = int(head)
    except ValueError:
        return -1, "NA", False
    if f == 1:
        return 1, s, True
    if f == 2:
        return 2, s, True
    return 0, s, True


def fetch_reads(bam_path, chrom, pos, pad=PAD):
    """Return dict[read_name] -> (raw_hp, ps_tag).

    Restrict to primary, mapped, non-secondary, non-supplementary reads.
    """
    out = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    start = max(0, pos - pad)
    end = pos + pad
    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate:
                continue
            try:
                hp = read.get_tag("HP")
            except KeyError:
                hp = None
            try:
                ps = read.get_tag("PS")
            except KeyError:
                ps = None
            # if a read appears multiple times (split alignments already filtered),
            # keep the first occurrence
            if read.query_name not in out:
                out[read.query_name] = (hp, ps)
    finally:
        bam.close()
    return out


# ---------------------------------------------------------------------------
# Metric calculators (operate on the shared read set)
# ---------------------------------------------------------------------------
def metric_family(reads_x, reads_pa, shared):
    """L1: % of shared reads where X family == PA family (both must be 1 or 2)."""
    n_eval = 0
    n_match = 0
    for rn in shared:
        fx, _, hx = normalize_hp(reads_x[rn][0])
        fp, _, hp_ = normalize_hp(reads_pa[rn][0])
        # family must be 1 or 2 on both sides for this metric
        if fx in (1, 2) and fp in (1, 2):
            n_eval += 1
            if fx == fp:
                n_match += 1
    return n_match, n_eval


def metric_exact(reads_x, reads_pa, shared):
    """L2: % of shared reads where canonical exact HP matches.

    Mapping convention: BL/V5 11 -> "1-1", 21 -> "2-1", 33 -> "3", 1 -> "1",
    2 -> "2"; PA strings used as-is. So BL HP=1 vs PA HP="1-1" -> exact mismatch
    even though family agrees.
    """
    n_eval = 0
    n_match = 0
    for rn in shared:
        _, ex, hx = normalize_hp(reads_x[rn][0])
        _, ep, hp_ = normalize_hp(reads_pa[rn][0])
        if hx and hp_:
            n_eval += 1
            if ex == ep:
                n_match += 1
    return n_match, n_eval


def hp1_hp2_ratio(reads, names):
    """ratio = HP1_count / (HP1_count + HP2_count) over the given read-name set."""
    f1 = 0
    f2 = 0
    for rn in names:
        if rn not in reads:
            continue
        fam, _, _ = normalize_hp(reads[rn][0])
        if fam == 1:
            f1 += 1
        elif fam == 2:
            f2 += 1
    tot = f1 + f2
    if tot == 0:
        return None, f1, f2
    return f1 / tot, f1, f2


def metric_ratio_distance(reads_x, reads_pa, shared):
    """L3: 1 - min(|rx-rp|, |(1-rx)-rp|).

    Returns (concordance_score in [0,1], ratio_x, ratio_pa, n_x_phased,
    n_pa_phased).
    """
    rx, fx1, fx2 = hp1_hp2_ratio(reads_x, shared)
    rp, fp1, fp2 = hp1_hp2_ratio(reads_pa, shared)
    if rx is None or rp is None:
        return None, rx, rp, fx1 + fx2, fp1 + fp2
    d_native = abs(rx - rp)
    d_flip = abs((1 - rx) - rp)
    score = 1.0 - min(d_native, d_flip)
    return score, rx, rp, fx1 + fx2, fp1 + fp2


def metric_orientation_corrected(reads_x, reads_pa, shared):
    """L4: per-PS best-of (native, flipped) family match rate.

    Group shared reads by their PA PS tag (treat None as a single bucket).
    Within each PS bucket, compute family-match under native and flipped
    (swap 1<->2) X labels. Take the max. Aggregate weighted by bucket size.
    """
    buckets = defaultdict(list)
    for rn in shared:
        ps_pa = reads_pa[rn][1]
        buckets[str(ps_pa)].append(rn)
    n_eval_total = 0
    n_match_total = 0
    for _, names in buckets.items():
        n_eval = 0
        n_native = 0
        n_flip = 0
        for rn in names:
            fx, _, _ = normalize_hp(reads_x[rn][0])
            fp, _, _ = normalize_hp(reads_pa[rn][0])
            if fx in (1, 2) and fp in (1, 2):
                n_eval += 1
                if fx == fp:
                    n_native += 1
                if (fx == 1 and fp == 2) or (fx == 2 and fp == 1):
                    n_flip += 1
        if n_eval == 0:
            continue
        n_match_total += max(n_native, n_flip)
        n_eval_total += n_eval
    return n_match_total, n_eval_total


# ---------------------------------------------------------------------------
# Per-site analysis
# ---------------------------------------------------------------------------
def analyze_site(label, chrom, pos):
    print(f"[{label}] {chrom}:{pos}", flush=True)
    reads_per_bam = {k: fetch_reads(v, chrom, pos) for k, v in BAMS.items()}
    counts = {k: len(v) for k, v in reads_per_bam.items()}

    # Core intersection: BL ∩ V5 ∩ PA
    shared = (
        set(reads_per_bam["BL"].keys())
        & set(reads_per_bam["V5"].keys())
        & set(reads_per_bam["PA"].keys())
    )
    n_shared = len(shared)

    # ------------------------------------------------------------------
    # L1 family
    # ------------------------------------------------------------------
    bl_l1 = metric_family(reads_per_bam["BL"], reads_per_bam["PA"], shared)
    v5_l1 = metric_family(reads_per_bam["V5"], reads_per_bam["PA"], shared)
    bl_l2 = metric_exact(reads_per_bam["BL"], reads_per_bam["PA"], shared)
    v5_l2 = metric_exact(reads_per_bam["V5"], reads_per_bam["PA"], shared)
    bl_l3 = metric_ratio_distance(reads_per_bam["BL"], reads_per_bam["PA"], shared)
    v5_l3 = metric_ratio_distance(reads_per_bam["V5"], reads_per_bam["PA"], shared)
    bl_l4 = metric_orientation_corrected(reads_per_bam["BL"], reads_per_bam["PA"], shared)
    v5_l4 = metric_orientation_corrected(reads_per_bam["V5"], reads_per_bam["PA"], shared)

    # Build per-read records for hp_family_exact.tsv
    per_read_rows = []
    for rn in sorted(shared):
        fb, eb, hb = normalize_hp(reads_per_bam["BL"][rn][0])
        f5, e5, h5 = normalize_hp(reads_per_bam["V5"][rn][0])
        fp, ep, hp = normalize_hp(reads_per_bam["PA"][rn][0])
        per_read_rows.append(
            dict(
                site=label,
                chrom=chrom,
                pos=pos,
                read_name=rn,
                BL_raw=str(reads_per_bam["BL"][rn][0]),
                V5_raw=str(reads_per_bam["V5"][rn][0]),
                PA_raw=str(reads_per_bam["PA"][rn][0]),
                BL_family=fb,
                V5_family=f5,
                PA_family=fp,
                BL_exact=eb,
                V5_exact=e5,
                PA_exact=ep,
                family_BL_match_PA=int(fb in (1, 2) and fp in (1, 2) and fb == fp),
                family_V5_match_PA=int(f5 in (1, 2) and fp in (1, 2) and f5 == fp),
                exact_BL_match_PA=int(hb and hp and eb == ep),
                exact_V5_match_PA=int(h5 and hp and e5 == ep),
            )
        )

    summary = dict(
        site=label,
        chrom=chrom,
        pos=pos,
        n_BL=counts["BL"],
        n_V2b=counts["V2b"],
        n_V3F=counts["V3F"],
        n_V5=counts["V5"],
        n_PA=counts["PA"],
        n_shared_BL_V5_PA=n_shared,
        # L1 family
        L1_BL_match=bl_l1[0],
        L1_BL_eval=bl_l1[1],
        L1_BL_rate=(bl_l1[0] / bl_l1[1]) if bl_l1[1] else np.nan,
        L1_V5_match=v5_l1[0],
        L1_V5_eval=v5_l1[1],
        L1_V5_rate=(v5_l1[0] / v5_l1[1]) if v5_l1[1] else np.nan,
        # L2 exact
        L2_BL_match=bl_l2[0],
        L2_BL_eval=bl_l2[1],
        L2_BL_rate=(bl_l2[0] / bl_l2[1]) if bl_l2[1] else np.nan,
        L2_V5_match=v5_l2[0],
        L2_V5_eval=v5_l2[1],
        L2_V5_rate=(v5_l2[0] / v5_l2[1]) if v5_l2[1] else np.nan,
        # L3 ratio distance (concordance score)
        L3_BL_score=bl_l3[0] if bl_l3[0] is not None else np.nan,
        L3_V5_score=v5_l3[0] if v5_l3[0] is not None else np.nan,
        L3_BL_ratio=bl_l3[1] if bl_l3[1] is not None else np.nan,
        L3_V5_ratio=v5_l3[1] if v5_l3[1] is not None else np.nan,
        L3_PA_ratio=bl_l3[2] if bl_l3[2] is not None else np.nan,
        # L4 orientation-corrected
        L4_BL_match=bl_l4[0],
        L4_BL_eval=bl_l4[1],
        L4_BL_rate=(bl_l4[0] / bl_l4[1]) if bl_l4[1] else np.nan,
        L4_V5_match=v5_l4[0],
        L4_V5_eval=v5_l4[1],
        L4_V5_rate=(v5_l4[0] / v5_l4[1]) if v5_l4[1] else np.nan,
    )
    return summary, per_read_rows


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_4metric_heatmap(df: pd.DataFrame, out_path: Path):
    metrics = [
        ("L1_BL_rate", "L1_V5_rate", "L1 family"),
        ("L2_BL_rate", "L2_V5_rate", "L2 exact"),
        ("L3_BL_score", "L3_V5_score", "L3 ratio (1-d)"),
        ("L4_BL_rate", "L4_V5_rate", "L4 orient-corr"),
    ]
    sites = df["site"].tolist()
    n_sites = len(sites)
    fig, axes = plt.subplots(1, 4, figsize=(20, max(6, n_sites * 0.42)), sharey=True)
    cmap = plt.get_cmap("RdYlGn")
    cmap.set_bad("lightgrey")
    for ax, (col_bl, col_v5, title) in zip(axes, metrics):
        mat = df[[col_bl, col_v5]].to_numpy(dtype=float)
        im = ax.imshow(mat, aspect="auto", vmin=0, vmax=1, cmap=cmap)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["BL", "V5"])
        ax.set_yticks(range(n_sites))
        ax.set_yticklabels(sites, fontsize=9)
        ax.set_title(title, fontsize=11)
        for i in range(n_sites):
            for j, val in enumerate(mat[i]):
                if np.isnan(val):
                    txt = "NA"
                    color = "black"
                else:
                    txt = f"{val:.2f}"
                    color = "black" if 0.3 < val < 0.85 else "white"
                ax.text(j, i, txt, ha="center", va="center", fontsize=7, color=color)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, label="concordance with PA")
    fig.suptitle(
        "V5 Audit fig 02a — read-intersection concordance (BL vs V5 vs Paired tumor)\n"
        "shared reads = BL ∩ V5 ∩ PA per site",
        fontsize=13,
    )
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_winloss_summary(df: pd.DataFrame, out_path: Path):
    metrics = [
        ("L1_BL_rate", "L1_V5_rate", "L1 family"),
        ("L2_BL_rate", "L2_V5_rate", "L2 exact"),
        ("L3_BL_score", "L3_V5_score", "L3 ratio"),
        ("L4_BL_rate", "L4_V5_rate", "L4 orient-corr"),
    ]
    rows = []
    means_bl = []
    means_v5 = []
    medians_bl = []
    medians_v5 = []
    for col_bl, col_v5, label in metrics:
        wins = losses = ties = na = 0
        for _, r in df.iterrows():
            a = r[col_bl]
            b = r[col_v5]
            if pd.isna(a) or pd.isna(b):
                na += 1
                continue
            if abs(a - b) < 1e-9:
                ties += 1
            elif b > a:
                wins += 1
            else:
                losses += 1
        rows.append((label, wins, ties, losses, na))
        means_bl.append(np.nanmean(df[col_bl].to_numpy(dtype=float)))
        means_v5.append(np.nanmean(df[col_v5].to_numpy(dtype=float)))
        medians_bl.append(np.nanmedian(df[col_bl].to_numpy(dtype=float)))
        medians_v5.append(np.nanmedian(df[col_v5].to_numpy(dtype=float)))

    labels = [r[0] for r in rows]
    wins = [r[1] for r in rows]
    ties = [r[2] for r in rows]
    losses = [r[3] for r in rows]
    nas = [r[4] for r in rows]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: stacked Win/Tie/Loss/NA bar
    ax = axes[0, 0]
    x = np.arange(len(labels))
    bot = np.zeros(len(labels))
    for vals, color, name in [
        (wins, "#2ca02c", "V5 wins (closer to PA)"),
        (ties, "#bdbdbd", "tie"),
        (losses, "#d62728", "BL wins"),
        (nas, "#7f7f7f", "NA / no eval"),
    ]:
        ax.bar(x, vals, bottom=bot, width=0.6, label=name, color=color, edgecolor="white")
        for xi, vi, bi in zip(x, vals, bot):
            if vi > 0:
                ax.text(xi, bi + vi / 2, str(vi), ha="center", va="center", fontsize=10, color="black")
        bot = bot + np.array(vals)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(f"# sites (of {len(df)})")
    ax.set_title("A. Win/Loss/Tie distribution")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_ylim(0, len(df) + 1)

    # Panel B: mean concordance bars (BL vs V5)
    ax = axes[0, 1]
    w = 0.38
    ax.bar(x - w / 2, means_bl, w, label="BL", color="#1f77b4")
    ax.bar(x + w / 2, means_v5, w, label="V5", color="#ff7f0e")
    for xi, b, v in zip(x, means_bl, means_v5):
        ax.text(xi - w / 2, b + 0.01, f"{b:.2f}", ha="center", fontsize=9)
        ax.text(xi + w / 2, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("mean concordance")
    ax.set_title("B. Mean concordance across 15 sites")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    # Panel C: per-site delta (V5 - BL)
    ax = axes[1, 0]
    site_labels = df["site"].tolist()
    deltas = {}
    for col_bl, col_v5, label in metrics:
        d = df[col_v5].to_numpy(dtype=float) - df[col_bl].to_numpy(dtype=float)
        deltas[label] = d
    delta_df = pd.DataFrame(deltas, index=site_labels)
    im = ax.imshow(delta_df.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m[2] for m in metrics], fontsize=10)
    ax.set_yticks(range(len(site_labels)))
    ax.set_yticklabels(site_labels, fontsize=8)
    ax.set_title("C. delta = V5 − BL  (red=V5 better, blue=BL better)")
    for i in range(len(site_labels)):
        for j in range(len(metrics)):
            v = delta_df.iat[i, j]
            if pd.isna(v):
                txt = "NA"
            else:
                txt = f"{v:+.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="black")
    plt.colorbar(im, ax=ax, shrink=0.85, label="V5 − BL")

    # Panel D: scatter L1 vs L4 to show orientation-correction effect
    ax = axes[1, 1]
    ax.scatter(df["L1_BL_rate"], df["L4_BL_rate"], c="#1f77b4", marker="o", s=70,
               edgecolors="white", label="BL")
    ax.scatter(df["L1_V5_rate"], df["L4_V5_rate"], c="#ff7f0e", marker="s", s=70,
               edgecolors="white", label="V5")
    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", alpha=0.5)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("L1 family rate")
    ax.set_ylabel("L4 orientation-corrected rate")
    ax.set_title("D. L1 vs L4 — orientation correction lifts both BL and V5")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    for _, r in df.iterrows():
        ax.annotate(r["site"], (r["L1_V5_rate"], r["L4_V5_rate"]),
                    fontsize=6, alpha=0.6, xytext=(3, 3), textcoords="offset points")

    fig.suptitle("V5 Audit fig 02b — Win/Loss/Tie + concordance summary (15 sites)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    summaries = []
    per_read_all = []
    for label, chrom, pos in SITES:
        s, rows = analyze_site(label, chrom, pos)
        summaries.append(s)
        per_read_all.extend(rows)

    df = pd.DataFrame(summaries)
    site_tsv = OUT_DATA_DIR / "per_site_concordance.tsv"
    df.to_csv(site_tsv, sep="\t", index=False, float_format="%.6f")

    df_reads = pd.DataFrame(per_read_all)
    reads_tsv = OUT_DATA_DIR / "hp_family_exact.tsv"
    df_reads.to_csv(reads_tsv, sep="\t", index=False)

    plot_4metric_heatmap(df, OUT_FIG_DIR / "fig02a_4metric_heatmap.png")
    plot_winloss_summary(df, OUT_FIG_DIR / "fig02b_winloss_summary.png")

    # Console summary
    print("\n=== Win/Loss/Tie (V5 vs BL, higher = closer to PA) ===")
    for col_bl, col_v5, label in [
        ("L1_BL_rate", "L1_V5_rate", "L1 family"),
        ("L2_BL_rate", "L2_V5_rate", "L2 exact"),
        ("L3_BL_score", "L3_V5_score", "L3 ratio"),
        ("L4_BL_rate", "L4_V5_rate", "L4 orient-corr"),
    ]:
        wins = losses = ties = na = 0
        for _, r in df.iterrows():
            a = r[col_bl]
            b = r[col_v5]
            if pd.isna(a) or pd.isna(b):
                na += 1
                continue
            if abs(a - b) < 1e-9:
                ties += 1
            elif b > a:
                wins += 1
            else:
                losses += 1
        print(f"  {label:18s} V5_wins={wins:2d}  ties={ties:2d}  BL_wins={losses:2d}  NA={na:2d}")

    print(f"\nWrote {site_tsv}")
    print(f"Wrote {reads_tsv}")
    print(f"Wrote {OUT_FIG_DIR / 'fig02a_4metric_heatmap.png'}")
    print(f"Wrote {OUT_FIG_DIR / 'fig02b_winloss_summary.png'}")


if __name__ == "__main__":
    main()
