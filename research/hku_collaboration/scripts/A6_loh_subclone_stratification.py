#!/usr/bin/env python3
"""A6 — LOH x CNV stratified HP / subclone analysis (HKU collaboration, 5/23).

Scope: 22 autosomes (BED has no chrX/chrY/chrM coverage; chrX included as
"unannotated baseline" only — all reads fall into LOH-outside x CN-neutral).

Zones (per-read assignment from BED):
  Z1 = LOH inside  x CN-gain
  Z2 = LOH inside  x CN-loss
  Z3 = LOH inside  x CN-neutral   (LOH only)
  Z4 = LOH outside x CN-gain
  Z5 = LOH outside x CN-loss
  Z6 = LOH outside x CN-neutral   (baseline)

Per-zone metrics:
  - HP tag counts: 1 / 2 / 1-1 / 2-1 / 3 / no-HP
  - subclone proxy ratio (HP1-1 side) = HP1-1 / (HP1 + HP1-1)
  - subclone proxy ratio (HP2-1 side) = HP2-1 / (HP2 + HP2-1)
  - mean depth = n_reads / (zone length, Mb)
  - per-chromosome breakdown

Stats:
  - Kruskal-Wallis on subclone proxy across 6 zones (per-chrom binned)
  - Pairwise Mann-Whitney U with Bonferroni
  - effect size (rank-biserial r)

Outputs:
  data/A6_loh_cnv_stratified_stats.tsv
  data/A6_loh_cnv_per_chr.tsv
  data/A6_loh_cnv_pairwise_stats.tsv
  figures/A6_loh_cnv_stratified_hp_distribution.png
  figures/A6_per_chr_loh_cnv_hp.png
  figures/A6_subclone_ratio_violin.png

Provenance: BAM via ClairS-TO ssrs v0.4.x -> LongPhase v1.7.3 haplotag --somaticMode
(see _common.py header + findings_5_23.md §1).
"""
from __future__ import annotations

import bisect
import multiprocessing as mp
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    BAM_PATH,
    DATA_DIR,
    FIG_DIR,
    LOH_BED,
    log_timing,
    setup_cjk_font,
)
import matplotlib.pyplot as plt

# ---- Constants ----
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
INCLUDE_CHRX = True  # included as unannotated baseline; results expected all -> Z6
TARGET_CHROMS = AUTOSOMES + (["chrX"] if INCLUDE_CHRX else [])

HP_LABELS_ORDER = ["1", "2", "1-1", "2-1", "3", "no-HP"]
HP_COLORS = {
    "1": "#2E7D32",     # germline H1
    "2": "#66BB6A",     # germline H2
    "1-1": "#1565C0",   # subclone H1-1
    "2-1": "#42A5F5",   # subclone H2-1
    "3": "#E65100",     # H3 / ambiguous
    "no-HP": "#9E9E9E", # untagged
}
ZONE_ORDER = [
    "Z1_LOHin_CNgain",
    "Z2_LOHin_CNloss",
    "Z3_LOHin_CNneutral",
    "Z4_LOHout_CNgain",
    "Z5_LOHout_CNloss",
    "Z6_LOHout_CNneutral",
]
ZONE_LABEL = {
    "Z1_LOHin_CNgain":    "LOH in × gain",
    "Z2_LOHin_CNloss":    "LOH in × loss",
    "Z3_LOHin_CNneutral": "LOH in × neutral",
    "Z4_LOHout_CNgain":   "LOH out × gain",
    "Z5_LOHout_CNloss":   "LOH out × loss",
    "Z6_LOHout_CNneutral": "LOH out × neutral",
}
ZONE_CN_WEIGHT = {
    "Z1_LOHin_CNgain":    2.0,
    "Z2_LOHin_CNloss":    0.5,
    "Z3_LOHin_CNneutral": 1.0,
    "Z4_LOHout_CNgain":   2.0,
    "Z5_LOHout_CNloss":   0.5,
    "Z6_LOHout_CNneutral": 1.0,
}

MIN_MAPQ = 20
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ----------------- BED loading -----------------
def load_bed_by_type(bed_path: str = LOH_BED) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """Returns {chrom: {type: sorted list of (start, end)}} for type in {gain, loss, loh}."""
    out: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip().split("\t")
            if len(p) < 4:
                continue
            chrom, s, e, t = p[0], int(p[1]), int(p[2]), p[3]
            if t not in {"gain", "loss", "loh"}:
                continue
            out[chrom][t].append((s, e))
    for c in out:
        for t in out[c]:
            out[c][t].sort()
    return out


def point_in_intervals(intervals: list[tuple[int, int]], pos: int) -> bool:
    if not intervals:
        return False
    starts = [iv[0] for iv in intervals]
    idx = bisect.bisect_right(starts, pos) - 1
    if idx < 0:
        return False
    s, e = intervals[idx]
    return s <= pos < e


def compute_zone_length_mb(chrom: str, bed_by_chrom: dict[str, dict[str, list[tuple[int, int]]]],
                            chrom_len: int) -> dict[str, float]:
    """Approximate per-zone genomic length (Mb) by 10 kb scan."""
    step = 10_000
    bed = bed_by_chrom.get(chrom, {})
    loh = bed.get("loh", [])
    gain = bed.get("gain", [])
    loss = bed.get("loss", [])
    zone_bp = {z: 0 for z in ZONE_ORDER}
    n_bins = (chrom_len + step - 1) // step
    for i in range(n_bins):
        mid = i * step + step // 2
        if mid >= chrom_len:
            break
        in_loh = point_in_intervals(loh, mid)
        in_gain = point_in_intervals(gain, mid)
        in_loss = point_in_intervals(loss, mid)
        z = classify_zone(in_loh, in_gain, in_loss)
        zone_bp[z] += step
    return {z: zone_bp[z] / 1e6 for z in zone_bp}


def classify_zone(in_loh: bool, in_gain: bool, in_loss: bool) -> str:
    # CNV gain and loss overlap is rare but resolve by gain dominance (loss is rare anyway)
    if in_gain and in_loss:
        in_loss = False
    if in_loh:
        if in_gain:
            return "Z1_LOHin_CNgain"
        if in_loss:
            return "Z2_LOHin_CNloss"
        return "Z3_LOHin_CNneutral"
    else:
        if in_gain:
            return "Z4_LOHout_CNgain"
        if in_loss:
            return "Z5_LOHout_CNloss"
        return "Z6_LOHout_CNneutral"


# ----------------- Per-chrom worker -----------------
def process_chrom(chrom: str):
    """Worker: scan all reads on `chrom`, return per-zone HP tag Counter and read positions list."""
    bam = pysam.AlignmentFile(BAM_PATH, "rb")
    if chrom not in bam.references:
        bam.close()
        return chrom, None, None, None

    chrom_len = dict(zip(bam.references, bam.lengths))[chrom]
    bed_by_chrom = load_bed_by_type()  # cheap re-load per worker (BED is tiny)
    bed = bed_by_chrom.get(chrom, {})
    loh = bed.get("loh", [])
    gain = bed.get("gain", [])
    loss = bed.get("loss", [])

    # Precompute starts arrays for fast bisect inside hot loop
    loh_starts = [iv[0] for iv in loh]
    gain_starts = [iv[0] for iv in gain]
    loss_starts = [iv[0] for iv in loss]

    def in_iv(intervals, starts, pos):
        if not intervals:
            return False
        idx = bisect.bisect_right(starts, pos) - 1
        if idx < 0:
            return False
        s, e = intervals[idx]
        return s <= pos < e

    # Per-zone HP counts + per-zone subclone-proxy bins per 1Mb (for per-chrom stats)
    zone_hp = {z: Counter() for z in ZONE_ORDER}
    # Bin (1Mb) -> zone -> (n_HP1, n_HP11, n_HP2, n_HP21)
    BIN_MB = 1_000_000
    n_bins = (chrom_len + BIN_MB - 1) // BIN_MB
    # accumulators: per bin, per zone
    bin_zone_acc = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0]))
    # [HP1, HP1-1, HP2, HP2-1]

    n_reads = 0
    t0 = time.time()
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate or read.mapping_quality < MIN_MAPQ:
            continue
        # Reference midpoint
        mid = (read.reference_start + (read.reference_end or read.reference_start)) // 2
        in_loh = in_iv(loh, loh_starts, mid)
        in_gain = in_iv(gain, gain_starts, mid)
        in_loss = in_iv(loss, loss_starts, mid)
        zone = classify_zone(in_loh, in_gain, in_loss)
        # HP tag
        if read.has_tag("HP"):
            hp_raw = read.get_tag("HP")
            hp = str(hp_raw)
        else:
            hp = "no-HP"
        if hp not in HP_LABELS_ORDER:
            # Unexpected value -> bucket into "3" if it looks like an integer 3, else no-HP
            hp = "3" if hp == "3" else "no-HP"
        zone_hp[zone][hp] += 1

        # per-bin sub-acc for subclone proxy stats
        bidx = mid // BIN_MB
        acc = bin_zone_acc[bidx][zone]
        if hp == "1":
            acc[0] += 1
        elif hp == "1-1":
            acc[1] += 1
        elif hp == "2":
            acc[2] += 1
        elif hp == "2-1":
            acc[3] += 1

        n_reads += 1
        if n_reads % 1_000_000 == 0:
            print(f"  [{chrom}] {n_reads:,} reads ({time.time()-t0:.0f}s)",
                  file=sys.stderr, flush=True)

    bam.close()
    # zone lengths (Mb)
    zone_lengths = compute_zone_length_mb(chrom, bed_by_chrom, chrom_len)

    # Flatten per-bin records to TSV rows
    bin_rows = []
    for bidx, zone_d in bin_zone_acc.items():
        for zone, (n_h1, n_h11, n_h2, n_h21) in zone_d.items():
            denom1 = n_h1 + n_h11
            denom2 = n_h2 + n_h21
            r1 = n_h11 / denom1 if denom1 >= 5 else np.nan
            r2 = n_h21 / denom2 if denom2 >= 5 else np.nan
            bin_rows.append({
                "chrom": chrom,
                "bin_start_mb": bidx,
                "zone": zone,
                "n_HP1": n_h1, "n_HP1_1": n_h11,
                "n_HP2": n_h2, "n_HP2_1": n_h21,
                "subclone_ratio_h1side": r1,
                "subclone_ratio_h2side": r2,
            })

    print(f"  [{chrom}] DONE: {n_reads:,} reads ({time.time()-t0:.0f}s)",
          file=sys.stderr, flush=True)
    return chrom, zone_hp, zone_lengths, bin_rows


# ----------------- Main -----------------
def main():
    setup_cjk_font()
    t_all = time.time()

    print(f"[A6] Targets: {len(TARGET_CHROMS)} chromosomes "
          f"({', '.join(TARGET_CHROMS)})", file=sys.stderr, flush=True)

    # Parallel per-chrom
    n_proc = min(8, len(TARGET_CHROMS))
    print(f"[A6] Launching {n_proc} workers...", file=sys.stderr, flush=True)
    with mp.Pool(n_proc) as pool:
        results = pool.map(process_chrom, TARGET_CHROMS)

    log_timing("All-chrom scan", t_all, time.time())

    # Aggregate
    global_zone_hp = {z: Counter() for z in ZONE_ORDER}
    global_zone_len = {z: 0.0 for z in ZONE_ORDER}
    per_chr_rows = []   # zone-level per chrom
    all_bin_rows = []   # 1Mb bin rows for stats

    for chrom, zh, zl, br in results:
        if zh is None:
            continue
        for z in ZONE_ORDER:
            global_zone_hp[z].update(zh[z])
            global_zone_len[z] += zl.get(z, 0.0)
            n_reads = sum(zh[z].values())
            n_h1 = zh[z].get("1", 0)
            n_h11 = zh[z].get("1-1", 0)
            n_h2 = zh[z].get("2", 0)
            n_h21 = zh[z].get("2-1", 0)
            denom1 = n_h1 + n_h11
            denom2 = n_h2 + n_h21
            per_chr_rows.append({
                "chrom": chrom,
                "zone": z,
                "n_reads": n_reads,
                "zone_length_mb": zl.get(z, 0.0),
                "depth_per_mb": n_reads / zl[z] if zl.get(z, 0.0) > 0 else np.nan,
                "n_HP1": n_h1, "n_HP2": n_h2,
                "n_HP1_1": n_h11, "n_HP2_1": n_h21,
                "n_HP3": zh[z].get("3", 0),
                "n_noHP": zh[z].get("no-HP", 0),
                "subclone_ratio_h1side": n_h11 / denom1 if denom1 > 0 else np.nan,
                "subclone_ratio_h2side": n_h21 / denom2 if denom2 > 0 else np.nan,
            })
        all_bin_rows.extend(br)

    # ---- Global stats TSV ----
    stats_rows = []
    for z in ZONE_ORDER:
        c = global_zone_hp[z]
        total = sum(c.values())
        n_h1 = c.get("1", 0)
        n_h11 = c.get("1-1", 0)
        n_h2 = c.get("2", 0)
        n_h21 = c.get("2-1", 0)
        denom1 = n_h1 + n_h11
        denom2 = n_h2 + n_h21
        zone_len = global_zone_len[z]
        for hp in HP_LABELS_ORDER:
            stats_rows.append({
                "zone": z,
                "zone_label": ZONE_LABEL[z],
                "hp_tag": hp,
                "count": c.get(hp, 0),
                "ratio": c.get(hp, 0) / total if total > 0 else 0.0,
                "zone_total_reads": total,
                "zone_length_mb": zone_len,
                "mean_depth_per_mb": total / zone_len if zone_len > 0 else np.nan,
                "cn_weight": ZONE_CN_WEIGHT[z],
                "subclone_ratio_h1side": n_h11 / denom1 if denom1 > 0 else np.nan,
                "subclone_ratio_h2side": n_h21 / denom2 if denom2 > 0 else np.nan,
            })

    df_stats = pd.DataFrame(stats_rows)
    df_stats.to_csv(DATA_DIR / "A6_loh_cnv_stratified_stats.tsv", sep="\t", index=False)
    print(f"[A6] wrote A6_loh_cnv_stratified_stats.tsv ({len(df_stats)} rows)",
          file=sys.stderr, flush=True)

    df_perchr = pd.DataFrame(per_chr_rows)
    df_perchr.to_csv(DATA_DIR / "A6_loh_cnv_per_chr.tsv", sep="\t", index=False)
    print(f"[A6] wrote A6_loh_cnv_per_chr.tsv ({len(df_perchr)} rows)",
          file=sys.stderr, flush=True)

    df_bins = pd.DataFrame(all_bin_rows)
    df_bins.to_csv(DATA_DIR / "A6_loh_cnv_1mb_bins.tsv", sep="\t", index=False)
    print(f"[A6] wrote A6_loh_cnv_1mb_bins.tsv ({len(df_bins)} rows)",
          file=sys.stderr, flush=True)

    # ---- Statistical tests (on 1Mb bin subclone ratios) ----
    pairwise_rows = []
    # Use h1-side ratio as primary metric (memo: H1/H1-1 is the established subclone direction
    # in HCC1395 per project_v5_v6_tradeoff_sp123 / hpfinengroups marker; H2-side reported too).
    for ratio_col in ["subclone_ratio_h1side", "subclone_ratio_h2side"]:
        groups = []
        labels = []
        for z in ZONE_ORDER:
            vals = df_bins[(df_bins["zone"] == z) & df_bins[ratio_col].notna()][ratio_col].values
            groups.append(vals)
            labels.append(z)
        non_empty_idx = [i for i, g in enumerate(groups) if len(g) >= 3]
        if len(non_empty_idx) >= 2:
            try:
                kw_stat, kw_p = stats.kruskal(*[groups[i] for i in non_empty_idx])
            except ValueError:
                kw_stat, kw_p = np.nan, np.nan
        else:
            kw_stat, kw_p = np.nan, np.nan
        print(f"[A6] Kruskal-Wallis ({ratio_col}): H={kw_stat:.3f} p={kw_p:.3e} "
              f"(across {len(non_empty_idx)} non-empty zones)",
              file=sys.stderr, flush=True)

        # Pairwise Mann-Whitney U
        from itertools import combinations
        n_pairs = len(list(combinations(non_empty_idx, 2)))
        for i, j in combinations(non_empty_idx, 2):
            a, b = groups[i], groups[j]
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            except ValueError:
                continue
            # rank-biserial r = 1 - 2U / (n1 * n2)
            r = 1.0 - (2.0 * u) / (len(a) * len(b))
            pairwise_rows.append({
                "metric": ratio_col,
                "zone_a": labels[i], "zone_b": labels[j],
                "n_a": len(a), "n_b": len(b),
                "median_a": float(np.median(a)),
                "median_b": float(np.median(b)),
                "mwu_U": float(u),
                "p_value": float(p),
                "p_bonferroni": float(min(1.0, p * max(1, n_pairs))),
                "rank_biserial_r": float(r),
                "kw_stat_zone_panel": float(kw_stat),
                "kw_p_zone_panel": float(kw_p),
            })

    df_pw = pd.DataFrame(pairwise_rows)
    df_pw.to_csv(DATA_DIR / "A6_loh_cnv_pairwise_stats.tsv", sep="\t", index=False)
    print(f"[A6] wrote A6_loh_cnv_pairwise_stats.tsv ({len(df_pw)} rows)",
          file=sys.stderr, flush=True)

    # ============ Figures ============
    plot_figure_1(global_zone_hp, global_zone_len)
    plot_figure_2(df_perchr)
    plot_figure_3(df_bins, df_pw)

    log_timing("A6 total", t_all, time.time())


# ----------------- Figures -----------------
def plot_figure_1(zone_hp: dict, zone_len: dict):
    """6-panel grid: per-zone stacked bar of HP tag fractions."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.flatten()
    for k, z in enumerate(ZONE_ORDER):
        ax = axes[k]
        c = zone_hp[z]
        total = sum(c.values())
        if total == 0:
            ax.text(0.5, 0.5, "No reads", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(ZONE_LABEL[z], fontsize=11)
            continue
        bottom = 0.0
        for hp in HP_LABELS_ORDER:
            ratio = c.get(hp, 0) / total
            ax.bar(0, ratio, bottom=bottom, color=HP_COLORS[hp],
                   edgecolor="white", linewidth=0.4,
                   label=f"{hp} ({c.get(hp,0):,}; {ratio*100:.1f}%)")
            bottom += ratio
        depth = total / zone_len[z] if zone_len[z] > 0 else np.nan
        ax.set_title(
            f"{ZONE_LABEL[z]}\n"
            f"n_reads={total:,}, length={zone_len[z]:.1f} Mb, "
            f"depth={depth:.0f}/Mb",
            fontsize=10,
        )
        ax.set_xticks([])
        ax.set_ylabel("Fraction" if k % 3 == 0 else "")
        ax.set_ylim(0, 1)
        ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  frameon=False)

    fig.suptitle(
        "A6 — HP tag distribution by LOH × CNV zone (HCC1395, 22 autosomes + chrX)\n"
        "BAM: ClairS-TO ssrs v0.4 → LongPhase v1.7.3 haplotag --somaticMode",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG_DIR / "A6_loh_cnv_stratified_hp_distribution.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[A6] wrote {out}", file=sys.stderr, flush=True)


def plot_figure_2(df_perchr: pd.DataFrame):
    """24-chr (incl chrX) panel: in-vs-out LOH subclone proxy comparison per chr."""
    chroms = TARGET_CHROMS
    n = len(chroms)
    ncol = 6
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(18, 2.6 * nrow), sharey=True)
    axes = axes.flatten()

    for k, chrom in enumerate(chroms):
        ax = axes[k]
        sub = df_perchr[df_perchr["chrom"] == chrom].set_index("zone")
        if sub.empty:
            ax.set_title(f"{chrom} (no data)", fontsize=9)
            ax.set_xticks([])
            continue
        labels = [ZONE_LABEL[z].replace(" × ", "×\n") for z in ZONE_ORDER]
        vals_h1 = [sub.loc[z, "subclone_ratio_h1side"] if z in sub.index else np.nan
                   for z in ZONE_ORDER]
        vals_h2 = [sub.loc[z, "subclone_ratio_h2side"] if z in sub.index else np.nan
                   for z in ZONE_ORDER]
        x = np.arange(len(ZONE_ORDER))
        ax.bar(x - 0.2, vals_h1, width=0.4, color="#1565C0", label="HP1-1 / (HP1+HP1-1)")
        ax.bar(x + 0.2, vals_h2, width=0.4, color="#42A5F5", label="HP2-1 / (HP2+HP2-1)")
        ax.set_title(chrom, fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6, rotation=45, ha="right")
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3, axis="y")

    # hide remaining axes
    for k in range(n, len(axes)):
        axes[k].axis("off")
    # legend on first panel
    axes[0].legend(fontsize=7, loc="upper right")

    fig.suptitle(
        "A6 — per-chromosome subclone proxy ratio by LOH × CNV zone\n"
        "y = HP{1-1 / 2-1} / (HP{1 / 2} + HP{1-1 / 2-1})",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "A6_per_chr_loh_cnv_hp.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[A6] wrote {out}", file=sys.stderr, flush=True)


def plot_figure_3(df_bins: pd.DataFrame, df_pw: pd.DataFrame):
    """6-violin of 1Mb-bin subclone proxy per zone, with KW + pairwise p annotation."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))

    for ax, ratio_col, title_tag in zip(
        axes, ["subclone_ratio_h1side", "subclone_ratio_h2side"],
        ["HP1 side", "HP2 side"],
    ):
        groups = []
        positions = []
        labels = []
        for k, z in enumerate(ZONE_ORDER):
            vals = df_bins[(df_bins["zone"] == z) & df_bins[ratio_col].notna()][ratio_col].values
            if len(vals) >= 3:
                groups.append(vals)
                positions.append(k)
                labels.append(f"{ZONE_LABEL[z]}\nn={len(vals)}")

        if groups:
            parts = ax.violinplot(groups, positions=positions, showmedians=True,
                                  showextrema=False, widths=0.75)
            for i, pc in enumerate(parts["bodies"]):
                z = ZONE_ORDER[positions[i]]
                pc.set_facecolor("#e76f51" if "LOHin" in z else "#2a9d8f")
                pc.set_alpha(0.75)
            # scatter overlay
            for pos, g in zip(positions, groups):
                xs = np.random.normal(pos, 0.04, len(g))
                ax.scatter(xs, g, s=6, color="black", alpha=0.25)

        # Kruskal-Wallis p
        sub_pw = df_pw[df_pw["metric"] == ratio_col]
        kw_p = sub_pw["kw_p_zone_panel"].iloc[0] if not sub_pw.empty else np.nan
        kw_h = sub_pw["kw_stat_zone_panel"].iloc[0] if not sub_pw.empty else np.nan

        ax.set_xticks(list(range(len(ZONE_ORDER))))
        ax.set_xticklabels([ZONE_LABEL[z] for z in ZONE_ORDER], rotation=30, ha="right", fontsize=9)
        ax.set_ylabel(f"Subclone proxy ratio ({title_tag})")
        ax.set_title(f"{title_tag}: Kruskal-Wallis H={kw_h:.2f}, p={kw_p:.2e}",
                     fontsize=10)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3, axis="y")

        # Significant pairs (Bonferroni < 0.05) annotate with bracket
        if not sub_pw.empty:
            sig = sub_pw[sub_pw["p_bonferroni"] < 0.05].copy()
            sig["a_idx"] = sig["zone_a"].apply(lambda x: ZONE_ORDER.index(x))
            sig["b_idx"] = sig["zone_b"].apply(lambda x: ZONE_ORDER.index(x))
            ymax = 1.0
            y_step = 0.045
            for k, (_, row) in enumerate(sig.iterrows()):
                y = ymax + (k + 1) * y_step
                ax.plot([row["a_idx"], row["b_idx"]], [y, y], color="k", linewidth=0.8)
                p = row["p_bonferroni"]
                star = "***" if p < 1e-3 else ("**" if p < 1e-2 else "*")
                ax.text((row["a_idx"] + row["b_idx"]) / 2, y + 0.005,
                        f"{star} p_bonf={p:.1e} r={row['rank_biserial_r']:+.2f}",
                        ha="center", fontsize=6)
            if len(sig) > 0:
                ax.set_ylim(-0.02, ymax + (len(sig) + 1) * y_step + 0.02)

    fig.suptitle(
        "A6 — Subclone proxy ratio per 1Mb bin by LOH × CNV zone\n"
        "Red = LOH inside (Z1/Z2/Z3) · Teal = LOH outside (Z4/Z5/Z6); "
        "* = Mann-Whitney U with Bonferroni",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIG_DIR / "A6_subclone_ratio_violin.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[A6] wrote {out}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
