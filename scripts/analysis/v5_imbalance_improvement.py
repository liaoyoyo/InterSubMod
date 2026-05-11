#!/usr/bin/env python3
"""V5 Audit Suite — Imbalance ratio + Per-site improvement analysis.

Computes HP1:HP2 imbalance ratio for Baseline (BL), V5, and Paired tumor BAM
across 15 audited sites (Phase4 TP / FP / V5_max reassign / Self-phasing extreme).

Outputs:
  - imbalance_ratio.tsv : per-site BL/V5/PA ratios + distance metrics
  - improvement_quantification.tsv : per-site improvement classification
  - 4 figures (fig04a, fig04b, fig05a, fig05b)
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pysam  # noqa: E402

# CJK font fallback per project rule
plt.rcParams["font.sans-serif"] = [
    "DejaVu Sans",
    "Droid Sans Fallback",
    "Noto Sans CJK TC",
    "Microsoft YaHei",
]
plt.rcParams["axes.unicode_minus"] = False

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BL_BAM = "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam"
V5_BAM = "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam"
PA_BAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"

OUT_ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite"
DATA_DIR = os.path.join(OUT_ROOT, "data")
FIG_04 = os.path.join(OUT_ROOT, "figures", "04_imbalance")
FIG_08 = os.path.join(OUT_ROOT, "figures", "08_synthesis")
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_04, exist_ok=True)
os.makedirs(FIG_08, exist_ok=True)

# 15 audited sites. Position is 1-based (matching VCF/IGV display).
# Window = +/- HALF_WINDOW around the site for HP counting.
HALF_WINDOW = 0  # use exact pileup at the site (overlapping reads only)

SITES: List[Tuple[str, str, int, str, str]] = [
    # (case, chrom, pos_1based, alt_or_label, comment_short)
    ("A_TP01", "chr6", 145444893, "G>A", "Phase4 TP_01 allele-only"),
    ("A_TP02", "chr4", 70548355, "G>A", "Phase4 TP_02 HP0 background"),
    ("A_TP03", "chr5", 153209947, "C>A", "Phase4 TP_03 low-methyl"),
    ("A_TP04", "chr16", 35118902, "G>A", "Phase4 TP_04 bimodal"),
    ("A_TP05", "chr7", 109185781, "G>T", "Phase4 TP_05 HP+Allele"),
    ("B_FPA1", "chr8", 93565727, "C>T", "FP_A1 low VAF"),
    ("B_FPA2", "chr9", 137953060, "T>C", "FP_A2 high CpG"),
    ("B_FPB1", "chr7", 52087777, "A>T", "FP_B1 HP-driven"),
    ("B_FPB2", "chr9", 75383880, "T>A", "FP_B2 MNP HP-driven"),
    ("C_V5max1", "chr19", 4639528, "(reassign)", "V5 39 reads HP33->HP11"),
    ("C_V5max2", "chr19", 2235521, "(reassign)", "V5 26 reads HP33->HP11"),
    ("C_V5max3", "chr19", 7405500, "(reassign)", "V5 18 reads HP33->HP21"),
    ("D_SP1", "chr19", 17565944, "(extreme)", "Self-phasing HP2=113 HP1=0"),
    ("D_SP2", "chr19", 12452332, "(extreme)", "Self-phasing 109:1"),
    ("D_SP3", "chr19", 12467180, "(extreme)", "Self-phasing 108:0"),
]


# ---------------------------------------------------------------------------
# HP-tag readers
# ---------------------------------------------------------------------------
@dataclass
class HPCounts:
    HP1: int = 0
    HP2: int = 0
    HP11: int = 0
    HP21: int = 0
    HP33: int = 0
    HP0: int = 0  # untagged

    @property
    def hp1_side(self) -> int:
        """HP1-side reads (HP1 + HP11)."""
        return self.HP1 + self.HP11

    @property
    def hp2_side(self) -> int:
        """HP2-side reads (HP2 + HP21)."""
        return self.HP2 + self.HP21

    @property
    def total_assigned(self) -> int:
        """HP1-side + HP2-side (excludes HP33 ambiguous and HP0 untagged)."""
        return self.hp1_side + self.hp2_side

    def ratio_hp1(self) -> float:
        """HP1 / (HP1 + HP2). Returns NaN when total_assigned == 0."""
        n = self.total_assigned
        return self.hp1_side / n if n > 0 else float("nan")


def _decode_hp_int(tag_value) -> str:
    """Decode HP:i:{1,2,11,21,33} into a HPCounts attribute name."""
    if tag_value == 1:
        return "HP1"
    if tag_value == 2:
        return "HP2"
    if tag_value == 11:
        return "HP11"
    if tag_value == 21:
        return "HP21"
    if tag_value == 33:
        return "HP33"
    return "HP0"


def _decode_hp_str(tag_value: str) -> str:
    """Decode HP:Z paired-tumor strings: '1','2','1-1','2-1' -> HP1-side/HP2-side."""
    s = str(tag_value)
    if s in ("1", "1-1"):
        return "HP1"  # treat both as HP1-side
    if s in ("2", "2-1"):
        return "HP2"
    return "HP0"


def count_hp_at_site(bam_path: str, chrom: str, pos_1based: int, hp_kind: str) -> HPCounts:
    """Count HP tags of reads overlapping pos_1based.

    hp_kind = 'int' for HP:i (BL/V5), 'str' for HP:Z (paired).
    """
    counts = HPCounts()
    pos0 = pos_1based - 1
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, pos0, pos0 + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate:
                continue
            try:
                hp = read.get_tag("HP")
                if hp_kind == "int":
                    field = _decode_hp_int(hp)
                else:
                    field = _decode_hp_str(hp)
            except KeyError:
                field = "HP0"
            setattr(counts, field, getattr(counts, field) + 1)
    return counts


# ---------------------------------------------------------------------------
# Distance utilities (orientation-flip aware)
# ---------------------------------------------------------------------------
def imbalance_distance(r_query: float, r_ref: float) -> float:
    """min(|r_query - r_ref|, |1 - r_query - r_ref|).

    The second form handles arbitrary HP1/HP2 label flip across PS blocks.
    Returns NaN when either input is NaN.
    """
    if np.isnan(r_query) or np.isnan(r_ref):
        return float("nan")
    return float(min(abs(r_query - r_ref), abs((1.0 - r_query) - r_ref)))


def imbalance_magnitude(r: float) -> float:
    """|ratio - 0.5|, NaN-safe."""
    if np.isnan(r):
        return float("nan")
    return float(abs(r - 0.5))


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def collect_data() -> List[Dict]:
    rows: List[Dict] = []
    for case, chrom, pos, alt, comment in SITES:
        bl = count_hp_at_site(BL_BAM, chrom, pos, "int")
        v5 = count_hp_at_site(V5_BAM, chrom, pos, "int")
        pa = count_hp_at_site(PA_BAM, chrom, pos, "str")

        bl_r = bl.ratio_hp1()
        v5_r = v5.ratio_hp1()
        pa_r = pa.ratio_hp1()

        bl_dist = imbalance_distance(bl_r, pa_r)
        v5_dist = imbalance_distance(v5_r, pa_r)
        delta_dist = (bl_dist - v5_dist) if not (np.isnan(bl_dist) or np.isnan(v5_dist)) else float("nan")

        rows.append(
            {
                "case": case,
                "chrom": chrom,
                "pos": pos,
                "alt": alt,
                "comment": comment,
                # Counts
                "BL_HP1": bl.HP1, "BL_HP11": bl.HP11, "BL_HP2": bl.HP2, "BL_HP21": bl.HP21,
                "BL_HP33": bl.HP33, "BL_HP0": bl.HP0,
                "V5_HP1": v5.HP1, "V5_HP11": v5.HP11, "V5_HP2": v5.HP2, "V5_HP21": v5.HP21,
                "V5_HP33": v5.HP33, "V5_HP0": v5.HP0,
                "PA_HP1side": pa.hp1_side, "PA_HP2side": pa.hp2_side, "PA_HP0": pa.HP0,
                # Ratios
                "BL_ratio": bl_r,
                "V5_ratio": v5_r,
                "PA_ratio": pa_r,
                # Imbalance magnitudes
                "BL_imbal": imbalance_magnitude(bl_r),
                "V5_imbal": imbalance_magnitude(v5_r),
                "PA_imbal": imbalance_magnitude(pa_r),
                # Distances vs paired
                "BL_dist": bl_dist,
                "V5_dist": v5_dist,
                "delta_dist": delta_dist,
            }
        )
    return rows


# ---------------------------------------------------------------------------
# TSV writers
# ---------------------------------------------------------------------------
def write_imbalance_tsv(rows: List[Dict]) -> str:
    path = os.path.join(DATA_DIR, "imbalance_ratio.tsv")
    cols = [
        "case", "chrom", "pos", "alt", "comment",
        "BL_HP1side", "BL_HP2side", "BL_HP33", "BL_HP0",
        "V5_HP1side", "V5_HP2side", "V5_HP33", "V5_HP0",
        "PA_HP1side", "PA_HP2side", "PA_HP0",
        "BL_ratio", "V5_ratio", "PA_ratio",
        "BL_imbal", "V5_imbal", "PA_imbal",
        "BL_dist_to_PA", "V5_dist_to_PA", "delta_dist_BL_minus_V5",
    ]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow(
                [
                    r["case"], r["chrom"], r["pos"], r["alt"], r["comment"],
                    r["BL_HP1"] + r["BL_HP11"], r["BL_HP2"] + r["BL_HP21"], r["BL_HP33"], r["BL_HP0"],
                    r["V5_HP1"] + r["V5_HP11"], r["V5_HP2"] + r["V5_HP21"], r["V5_HP33"], r["V5_HP0"],
                    r["PA_HP1side"], r["PA_HP2side"], r["PA_HP0"],
                    f"{r['BL_ratio']:.4f}", f"{r['V5_ratio']:.4f}", f"{r['PA_ratio']:.4f}",
                    f"{r['BL_imbal']:.4f}", f"{r['V5_imbal']:.4f}", f"{r['PA_imbal']:.4f}",
                    f"{r['BL_dist']:.4f}", f"{r['V5_dist']:.4f}", f"{r['delta_dist']:.4f}",
                ]
            )
    return path


def classify(delta: float) -> str:
    if np.isnan(delta):
        return "NA"
    if delta > 0.10:
        return "strong_improve"
    if delta > 0.02:
        return "moderate_improve"
    if delta < -0.02:
        return "regression"
    return "neutral"


def write_improvement_tsv(rows: List[Dict]) -> str:
    path = os.path.join(DATA_DIR, "improvement_quantification.tsv")
    cols = [
        "rank", "case", "chrom", "pos", "comment",
        "BL_dist", "V5_dist", "delta_dist",
        "BL_imbal", "V5_imbal", "PA_imbal",
        "delta_imbal_BL_minus_V5",
        "category",
    ]
    sorted_rows = sorted(
        rows,
        key=lambda r: (-(r["delta_dist"] if not np.isnan(r["delta_dist"]) else -1e9)),
    )
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for i, r in enumerate(sorted_rows, 1):
            d_imbal = r["BL_imbal"] - r["V5_imbal"] if not (
                np.isnan(r["BL_imbal"]) or np.isnan(r["V5_imbal"])
            ) else float("nan")
            w.writerow(
                [
                    i, r["case"], r["chrom"], r["pos"], r["comment"],
                    f"{r['BL_dist']:.4f}", f"{r['V5_dist']:.4f}", f"{r['delta_dist']:.4f}",
                    f"{r['BL_imbal']:.4f}", f"{r['V5_imbal']:.4f}", f"{r['PA_imbal']:.4f}",
                    f"{d_imbal:.4f}",
                    classify(r["delta_dist"]),
                ]
            )
    return path


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
CAT_COLOR = {
    "strong_improve": "#2ca02c",
    "moderate_improve": "#88c999",
    "neutral": "#bdbdbd",
    "regression": "#d62728",
    "NA": "#000000",
}


def fig04a_ratio_scatter(rows: List[Dict], out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 6.2))
    xs = np.arange(len(rows))
    bl = [r["BL_ratio"] for r in rows]
    v5 = [r["V5_ratio"] for r in rows]
    pa = [r["PA_ratio"] for r in rows]

    # Connect BL->V5->PA per site with thin grey line
    for i, (b, v, p) in enumerate(zip(bl, v5, pa)):
        ys = [y for y in [b, v, p] if not np.isnan(y)]
        if len(ys) >= 2:
            ax.plot([i, i, i], [b, v, p], color="#cccccc", lw=0.6, zorder=1)

    ax.scatter(xs, bl, marker="o", s=50, color="#1f77b4", label="Baseline", zorder=3, edgecolor="black", lw=0.4)
    ax.scatter(xs, v5, marker="s", s=50, color="#ff7f0e", label="V5", zorder=3, edgecolor="black", lw=0.4)
    ax.scatter(xs, pa, marker="^", s=60, color="#2ca02c", label="Paired", zorder=3, edgecolor="black", lw=0.4)

    ax.axhline(0.5, color="grey", ls="--", lw=0.8, alpha=0.6)
    ax.set_xticks(xs)
    ax.set_xticklabels([r["case"] for r in rows], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("HP1 ratio = HP1side / (HP1side + HP2side)")
    ax.set_title("Per-site HP1 ratio: Baseline vs V5 vs Paired (15 audited sites)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25)
    # Annotate NaN PA sites
    for i, p in enumerate(pa):
        if np.isnan(p):
            ax.annotate("PA=NA", (i, 0.5), fontsize=7, color="red", ha="center",
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="red"))
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig04b_distance_distribution(rows: List[Dict], out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    bl_d = np.array([r["BL_dist"] for r in rows], dtype=float)
    v5_d = np.array([r["V5_dist"] for r in rows], dtype=float)
    bl_d = bl_d[~np.isnan(bl_d)]
    v5_d = v5_d[~np.isnan(v5_d)]
    bins = np.linspace(0, 0.55, 12)
    ax.hist(bl_d, bins=bins, alpha=0.55, color="#1f77b4",
            label=f"Baseline (mean={bl_d.mean():.3f}, median={np.median(bl_d):.3f})", edgecolor="black")
    ax.hist(v5_d, bins=bins, alpha=0.55, color="#ff7f0e",
            label=f"V5 (mean={v5_d.mean():.3f}, median={np.median(v5_d):.3f})", edgecolor="black")
    ax.axvline(bl_d.mean(), color="#1f77b4", ls="--", lw=1.0)
    ax.axvline(v5_d.mean(), color="#ff7f0e", ls="--", lw=1.0)
    ax.set_xlabel("Imbalance distance to paired (orientation-corrected)")
    ax.set_ylabel("Number of sites")
    ax.set_title("Distance-to-Paired distribution: Baseline vs V5 (15 sites)")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)

    # Inset: paired distance per-site comparison
    bl_d_full = np.array([r["BL_dist"] for r in rows], dtype=float)
    v5_d_full = np.array([r["V5_dist"] for r in rows], dtype=float)
    mask = ~(np.isnan(bl_d_full) | np.isnan(v5_d_full))
    inset = fig.add_axes([0.6, 0.42, 0.34, 0.40])
    inset.scatter(bl_d_full[mask], v5_d_full[mask], color="#444444", s=24, edgecolor="black", lw=0.4)
    lim = max(bl_d_full[mask].max(), v5_d_full[mask].max()) * 1.1
    inset.plot([0, lim], [0, lim], color="grey", ls="--", lw=0.7)
    inset.set_xlabel("BL_dist", fontsize=8)
    inset.set_ylabel("V5_dist", fontsize=8)
    inset.set_title("paired BL vs V5 distance", fontsize=9)
    inset.tick_params(labelsize=7)
    inset.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig05a_per_site_bar(rows: List[Dict], out_path: str) -> None:
    sorted_rows = sorted(
        rows, key=lambda r: -(r["delta_dist"] if not np.isnan(r["delta_dist"]) else -1e9)
    )
    deltas = [r["delta_dist"] for r in sorted_rows]
    cats = [classify(d) for d in deltas]
    colors = [CAT_COLOR[c] for c in cats]
    cases = [r["case"] for r in sorted_rows]

    fig, ax = plt.subplots(figsize=(10.4, 6.0))
    xs = np.arange(len(sorted_rows))
    ax.bar(xs, deltas, color=colors, edgecolor="black", lw=0.5)
    ax.axhline(0.0, color="black", lw=0.6)
    ax.axhline(0.10, color="#2ca02c", ls=":", lw=0.7, label="strong (>0.10)")
    ax.axhline(0.02, color="#88c999", ls=":", lw=0.7, label="moderate (>0.02)")
    ax.axhline(-0.02, color="#d62728", ls=":", lw=0.7, label="regression (<-0.02)")
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("delta_dist = BL_dist - V5_dist (positive = V5 closer to paired)")
    ax.set_title("Per-site V5 improvement (sorted by magnitude, 15 sites)")

    counts = {k: cats.count(k) for k in ["strong_improve", "moderate_improve", "neutral", "regression"]}
    txt = (
        f"strong: {counts['strong_improve']}   "
        f"moderate: {counts['moderate_improve']}   "
        f"neutral: {counts['neutral']}   "
        f"regression: {counts['regression']}"
    )
    ax.text(0.01, 0.98, txt, transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="grey"))
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig05b_rank_heatmap(rows: List[Dict], out_path: str) -> None:
    """4 metric x 15 sites V5-BL improvement heatmap.

    Metrics (all V5-BL signed; positive = V5 better):
      M1 delta_dist        : (BL_dist - V5_dist) to paired
      M2 delta_imbal       : (BL_imbal - V5_imbal) — moves toward 50:50
      M3 hp33_reduction    : (BL_HP33 - V5_HP33) / max(BL_HP33,1) — fraction of ambiguous resolved
      M4 hp0_change        : -(V5_HP0 - BL_HP0) / max(BL_HP0,1) — drop in untagged (positive = better)
    """
    sorted_rows = sorted(
        rows, key=lambda r: -(r["delta_dist"] if not np.isnan(r["delta_dist"]) else -1e9)
    )
    cases = [r["case"] for r in sorted_rows]
    M1 = [r["delta_dist"] for r in sorted_rows]
    M2 = [r["BL_imbal"] - r["V5_imbal"] for r in sorted_rows]
    M3 = [
        ((r["BL_HP33"] - r["V5_HP33"]) / max(r["BL_HP33"], 1)) if r["BL_HP33"] > 0 else 0.0
        for r in sorted_rows
    ]
    M4 = [
        ((r["BL_HP0"] - r["V5_HP0"]) / max(r["BL_HP0"], 1)) if r["BL_HP0"] > 0 else 0.0
        for r in sorted_rows
    ]
    mat = np.array([M1, M2, M3, M4], dtype=float)

    fig, ax = plt.subplots(figsize=(11.5, 4.4))
    vmax = float(np.nanmax(np.abs(mat))) or 0.1
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(cases)))
    ax.set_xticklabels(cases, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(
        [
            "M1 dist_to_paired (V5 closer if +)",
            "M2 imbalance reduction (V5 more balanced if +)",
            "M3 HP33 ambiguous resolved fraction",
            "M4 HP0 untagged reduction fraction",
        ],
        fontsize=9,
    )
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7,
                    color="black" if abs(v) < vmax * 0.55 else "white")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("V5 - BL (positive = V5 improvement)")
    ax.set_title("4-metric x 15-site V5 improvement heatmap (sorted by M1)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    rows = collect_data()
    p1 = write_imbalance_tsv(rows)
    p2 = write_improvement_tsv(rows)
    print(f"Wrote {p1}")
    print(f"Wrote {p2}")

    fig04a = os.path.join(FIG_04, "fig04a_ratio_scatter.png")
    fig04b = os.path.join(FIG_04, "fig04b_distance_distribution.png")
    fig05a = os.path.join(FIG_08, "fig05a_per_site_improvement_bar.png")
    fig05b = os.path.join(FIG_08, "fig05b_improvement_rank_heatmap.png")

    fig04a_ratio_scatter(rows, fig04a)
    fig04b_distance_distribution(rows, fig04b)
    fig05a_per_site_bar(rows, fig05a)
    fig05b_rank_heatmap(rows, fig05b)

    for p in [fig04a, fig04b, fig05a, fig05b]:
        sz = os.path.getsize(p)
        print(f"Wrote {p} ({sz} bytes)")

    # Summary stats
    bl_d = np.array([r["BL_dist"] for r in rows])
    v5_d = np.array([r["V5_dist"] for r in rows])
    delta = np.array([r["delta_dist"] for r in rows])
    cats = [classify(d) for d in delta]
    print("\n=== Summary ===")
    print(f"BL_dist mean={np.nanmean(bl_d):.3f} median={np.nanmedian(bl_d):.3f}")
    print(f"V5_dist mean={np.nanmean(v5_d):.3f} median={np.nanmedian(v5_d):.3f}")
    print(f"delta_dist mean={np.nanmean(delta):.3f}  V5_better={int((delta>0).sum())} sites")
    print(f"strong_improve={cats.count('strong_improve')}  moderate={cats.count('moderate_improve')}  "
          f"neutral={cats.count('neutral')}  regression={cats.count('regression')}")


if __name__ == "__main__":
    main()
