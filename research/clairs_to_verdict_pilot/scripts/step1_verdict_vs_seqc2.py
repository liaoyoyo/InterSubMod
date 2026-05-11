#!/usr/bin/env python3
"""Step 1: Verdict × SEQC2 truth cross-tabulation.

Context:
    ClairS-TO v0.3.0 outputs per-variant Verdict tags (Germline/Somatic/
    SubclonalSomatic) based on ASCAT-derived purity/ploidy + binomial tests.
    The HCC1395 t20_n30 subsample is the only sample in this project with
    populated Verdict annotations (14,875 tags total). We cross-tabulate these
    tags against SEQC2 sSNV truth to judge whether the Verdict module gives a
    meaningful germline-vs-somatic signal for this project.

    Hypotheses (declared in plan):
        H-V1  Verdict_Germline FP rate >= 70%  → Verdict correctly separates
              germline residue from somatic signal.
        H-V2  Verdict_Somatic TP rate >= 85%   → Verdict's somatic branch is
              not dropping true positives.

Outputs:
    data/step1_verdict_vs_seqc2.tsv         — per-variant annotation
    data/step1_confusion_matrix.tsv         — Verdict class × {TP, FP} + Wilson CI
    figures/step1_verdict_fp_rate_per_class.png
"""
from __future__ import annotations

import argparse
import bisect
import math
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pysam


VERDICT_FLAGS = ("Verdict_Germline", "Verdict_Somatic", "Verdict_SubclonalSomatic")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="ClairS-TO VCF with Verdict tags")
    p.add_argument("--truth", required=True, help="SEQC2 sSNV truth VCF")
    p.add_argument("--bed", required=True, help="SEQC2 HighConf BED")
    p.add_argument("--out-dir", required=True, help="Output directory (contains data/ and figures/)")
    p.add_argument("--truth-filters", default="PASS,HighConf",
                   help="Comma-separated FILTER values to accept as truth (default: PASS,HighConf)")
    return p.parse_args()


def load_highconf_bed(path: str) -> dict[str, tuple[list[int], list[int]]]:
    """Return {chrom: (starts, ends)} where starts/ends are sorted for bisect."""
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            by_chrom[chrom].append((start, end))
    flat: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        starts = [s for s, _ in ivs]
        ends = [e for _, e in ivs]
        flat[chrom] = (starts, ends)
    return flat


def in_bed(chrom: str, pos0: int, bed: dict[str, tuple[list[int], list[int]]]) -> bool:
    """Check pos0 (0-based) in BED intervals [start, end)."""
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    idx = bisect.bisect_right(starts, pos0) - 1
    if idx < 0:
        return False
    return pos0 < ends[idx]


def load_truth_set(path: str, accepted_filters: set[str]) -> set[tuple[str, int, str, str]]:
    truth: set[tuple[str, int, str, str]] = set()
    vcf = pysam.VariantFile(path)
    for rec in vcf:
        # FILTER field: normalize
        if rec.filter.keys():
            filters = set(rec.filter.keys())
        else:
            filters = {"PASS"}
        if accepted_filters and not (filters & accepted_filters):
            continue
        if rec.ref is None or not rec.alts:
            continue
        ref = rec.ref.upper()
        for alt in rec.alts:
            # sSNV only (single base)
            if len(ref) != 1 or len(alt) != 1:
                continue
            truth.add((rec.chrom, rec.pos, ref, alt.upper()))
    return truth


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def classify_verdict(info) -> str:
    flags = [f for f in VERDICT_FLAGS if info.get(f, False)]
    if len(flags) == 0:
        return "None"
    if len(flags) == 1:
        return flags[0]
    return "Multi(" + "+".join(flags) + ")"


def main() -> None:
    args = parse_args()
    data_dir = os.path.join(args.out_dir, "data") if os.path.isdir(os.path.join(args.out_dir, "data")) else args.out_dir
    fig_dir_candidate = os.path.join(os.path.dirname(data_dir), "figures")
    fig_dir = fig_dir_candidate if os.path.isdir(fig_dir_candidate) else data_dir

    accepted = set(f.strip() for f in args.truth_filters.split(",") if f.strip())

    print(f"[step1] loading HighConf BED … {args.bed}")
    bed = load_highconf_bed(args.bed)
    n_regions = sum(len(v[0]) for v in bed.values())
    print(f"        {len(bed)} chroms, {n_regions:,} intervals")

    print(f"[step1] loading SEQC2 truth (FILTER ∈ {sorted(accepted)}) … {args.truth}")
    truth = load_truth_set(args.truth, accepted)
    print(f"        {len(truth):,} truth sSNVs")

    print(f"[step1] scanning ClairS-TO VCF … {args.vcf}")
    rows: list[dict] = []
    vcf = pysam.VariantFile(args.vcf)
    for rec in vcf:
        if rec.ref is None or not rec.alts:
            continue
        ref = rec.ref.upper()
        # take first alt only (ClairS-TO outputs single-alt records)
        alt = rec.alts[0].upper()
        if len(ref) != 1 or len(alt) != 1:
            continue

        filters = list(rec.filter.keys()) if rec.filter.keys() else ["PASS"]
        filter_str = ";".join(sorted(filters))
        is_pass = "PASS" in filters and len(filters) == 1

        # Verdict class
        verdict_cls = classify_verdict(rec.info)

        # in HighConf
        pos0 = rec.pos - 1  # BED is 0-based; VCF pos is 1-based
        in_hc = in_bed(rec.chrom, pos0, bed)

        # truth label (only meaningful when in_hc=True)
        key = (rec.chrom, rec.pos, ref, alt)
        seqc2_label = "TP" if key in truth else "FP"

        # AF / DP (FORMAT field; sample index 0)
        try:
            fmt = rec.samples[0]
            af = fmt.get("AF")
            dp = fmt.get("DP")
            if isinstance(af, tuple):
                af = af[0] if af else None
        except Exception:
            af, dp = None, None

        qual = rec.qual if rec.qual is not None else ""

        rows.append({
            "chrom": rec.chrom,
            "pos": rec.pos,
            "ref": ref,
            "alt": alt,
            "qual": qual,
            "filter": filter_str,
            "is_pass": int(is_pass),
            "verdict": verdict_cls,
            "in_highconf": int(in_hc),
            "seqc2_label": seqc2_label,
            "af": af,
            "dp": dp,
        })

    df = pd.DataFrame(rows)
    print(f"        {len(df):,} sSNV records parsed")

    out_tsv = os.path.join(data_dir, "step1_verdict_vs_seqc2.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[step1] wrote {out_tsv}")

    # Confusion matrix restricted to in_highconf=True (truth set complete only there).
    sub = df[df["in_highconf"] == 1].copy()
    print(f"[step1] in-HighConf subset: {len(sub):,}")

    cm_rows = []
    for subset_name, subset_df in [
        ("ALL", sub),
        ("PASS_only", sub[sub["is_pass"] == 1]),
        ("LowQual_only", sub[sub["filter"].str.contains("LowQual")]),
    ]:
        for v in ["Verdict_Germline", "Verdict_Somatic", "Verdict_SubclonalSomatic", "None"]:
            vdf = subset_df[subset_df["verdict"] == v]
            n = len(vdf)
            tp = (vdf["seqc2_label"] == "TP").sum()
            fp = (vdf["seqc2_label"] == "FP").sum()
            fp_rate = fp / n if n else float("nan")
            tp_rate = tp / n if n else float("nan")
            lo, hi = wilson_ci(fp, n)
            cm_rows.append({
                "subset": subset_name,
                "verdict": v,
                "n": n,
                "tp": int(tp),
                "fp": int(fp),
                "tp_rate": tp_rate,
                "fp_rate": fp_rate,
                "fp_rate_ci_lo": lo,
                "fp_rate_ci_hi": hi,
            })
    cm = pd.DataFrame(cm_rows)
    cm_tsv = os.path.join(data_dir, "step1_confusion_matrix.tsv")
    cm.to_csv(cm_tsv, sep="\t", index=False)
    print(f"[step1] wrote {cm_tsv}")

    # ---- Figure: FP rate per Verdict class within in-HighConf ALL subset ----
    plot_df = cm[cm["subset"] == "ALL"].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    xs = list(range(len(plot_df)))
    heights = plot_df["fp_rate"].fillna(0).tolist()
    err_lo = [h - lo if pd.notna(lo) else 0 for h, lo in zip(heights, plot_df["fp_rate_ci_lo"])]
    err_hi = [hi - h if pd.notna(hi) else 0 for h, hi in zip(heights, plot_df["fp_rate_ci_hi"])]
    colors = {"Verdict_Germline": "#d62728", "Verdict_Somatic": "#2ca02c",
              "Verdict_SubclonalSomatic": "#ff7f0e", "None": "#7f7f7f"}
    bar_colors = [colors.get(v, "#cccccc") for v in plot_df["verdict"]]
    ax.bar(xs, heights, yerr=[err_lo, err_hi], color=bar_colors, alpha=0.85,
           capsize=4, edgecolor="black")
    for x, h, n in zip(xs, heights, plot_df["n"]):
        ax.text(x, h + 0.02, f"n={n}\n{h:.1%}", ha="center", fontsize=9)
    ax.axhline(0.70, color="red", linestyle="--", alpha=0.6, label="H-V1 gate (70%)")
    ax.set_xticks(xs)
    ax.set_xticklabels(plot_df["verdict"], rotation=15, ha="right")
    ax.set_ylabel("FP rate (vs SEQC2 truth)")
    ax.set_ylim(0, 1.08)
    ax.set_title("Verdict class × SEQC2 FP rate (HCC1395 t20_n30, in HighConf BED, all FILTER)")
    ax.legend(loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "step1_verdict_fp_rate_per_class.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[step1] wrote {fig_path}")

    # ---- Stdout summary ----
    print("\n" + "=" * 70)
    print("Step 1 summary — Verdict × SEQC2 (in_highconf=1)")
    print("=" * 70)
    for _, r in cm.iterrows():
        if r["n"] == 0:
            continue
        print(f"  [{r['subset']:<13}] {r['verdict']:<27}  n={r['n']:>6}  "
              f"TP={r['tp']:>5}  FP={r['fp']:>5}  "
              f"FP_rate={r['fp_rate']:.3f} [CI {r['fp_rate_ci_lo']:.3f}, {r['fp_rate_ci_hi']:.3f}]")

    # ---- Hypothesis gate check ----
    print("\n--- Hypothesis gates ---")
    germ_all = cm[(cm["subset"] == "ALL") & (cm["verdict"] == "Verdict_Germline")].iloc[0]
    som_pass = cm[(cm["subset"] == "PASS_only") & (cm["verdict"] == "Verdict_Somatic")]
    germ_fp_rate = germ_all["fp_rate"] if germ_all["n"] else float("nan")
    h_v1 = "POSITIVE" if germ_fp_rate >= 0.70 else ("WEAK" if germ_fp_rate >= 0.55 else "NEGATIVE")
    print(f"  H-V1  Verdict_Germline FP rate = {germ_fp_rate:.3f}  → {h_v1} (gate ≥0.70)")
    if len(som_pass):
        r = som_pass.iloc[0]
        som_tp = r["tp_rate"] if r["n"] else float("nan")
        h_v2 = "POSITIVE" if som_tp >= 0.85 else ("WEAK" if r["fp_rate"] <= 0.20 else "NEGATIVE")
        print(f"  H-V2  Verdict_Somatic TP rate (PASS) = {som_tp:.3f}  → {h_v2} (gate ≥0.85)")
    print("=" * 70)


if __name__ == "__main__":
    main()
