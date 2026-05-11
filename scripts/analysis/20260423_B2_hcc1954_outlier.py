#!/usr/bin/env python3
"""
B2 analysis: HCC1954 Outer cross-het TP rate outlier root-cause.

Context (from obs18):
  HCC1954 Inner same-hap TP = 0.43 (other TO samples 0.76-0.96)
  HCC1954 Outer cross-het TP = 0.08 (other TO samples 0.50-0.88)
Both far below the cross-sample cohort.

Hypotheses:
  H1  Potential_LOH annotation reliability issue specific to HCC1954
  H2  AF / CovM distribution shape differs (more low-AF variants -> more germline FPs)
  H3  HER2 (chr17) / MYC (chr8) amplicons drive FP enrichment

Outputs:
  4-panel figure: AF dist / CovM dist / chr stacked / TP/FP ratio by chr
  Per-chr stats table with HCC1954 vs cross-sample comparison
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUT_FIG_DIR = ROOT / "docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier"
OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV_DIR = OUT_FIG_DIR

BUCKET_COLS = ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]
NEED_COLS = [
    "RegionID", "Chr", "Pos", "AlleleDelta", "HPFineNGroups",
    "Potential_LOH", "Coverage_Multiple", "NTumorReads", "NNormalReads",
    "NumReads",
] + BUCKET_COLS

SAMPLES_TO = {
    "HCC1395": {
        "tp": "/tmp/ism_hp_fix_phase1/tp_off/significance_summary.csv",
        "fp": "/tmp/ism_hp_fix_phase1/fp_off/significance_summary.csv",
    },
    "HCC1395_DORADO": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "HCC1937": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "HCC1954": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "H2009": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "H1437": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
}


def categorize_combo(has_hp1: int, has_hp1s: int, has_hp2: int, has_hp2s: int) -> str:
    if has_hp1 and has_hp1s and not has_hp2 and not has_hp2s:
        return "same_HP1 (HP1 + HP1-1)"
    if has_hp2 and has_hp2s and not has_hp1 and not has_hp1s:
        return "same_HP2 (HP2 + HP2-1)"
    if has_hp1 and has_hp2s and not has_hp1s and not has_hp2:
        return "cross_het (HP1 + HP2-1)"
    if has_hp1s and has_hp2 and not has_hp1 and not has_hp2s:
        return "cross_het_inv (HP1-1 + HP2)"
    return "other"


def load_sample(name: str, paths: dict) -> pd.DataFrame:
    frames = []
    for label, col in [("tp", 1), ("fp", 0)]:
        p = Path(paths[label])
        if not p.exists():
            print(f"  [!] {name} {label}: missing {p}")
            continue
        # Only keep NEED_COLS that actually exist
        head = pd.read_csv(p, nrows=0, low_memory=False)
        usecols = [c for c in NEED_COLS if c in head.columns]
        df = pd.read_csv(p, low_memory=False, usecols=usecols)
        df["tp_label"] = col
        df["sample"] = name
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def prep(df: pd.DataFrame) -> pd.DataFrame:
    """Add combo / loh_side / AF_class columns (following obs18 convention)."""
    if df.empty:
        return df
    ng2 = df[df["HPFineNGroups"] == 2].copy()
    if len(ng2) == 0:
        return ng2
    for c in BUCKET_COLS:
        ng2[c] = ng2[c].fillna(0).astype(int)
    ng2["combo"] = ng2.apply(
        lambda r: categorize_combo(
            int(r["HPFineN_HP1"] > 0), int(r["HPFineN_HP1S"] > 0),
            int(r["HPFineN_HP2"] > 0), int(r["HPFineN_HP2S"] > 0),
        ), axis=1,
    )
    ng2["loh_side"] = ng2["Potential_LOH"].apply(
        lambda x: "Inner" if str(x).lower() in ("true", "1") else "Outer"
    )
    ng2["AF_class"] = ng2["AlleleDelta"].clip(0, 1).apply(
        lambda af: "Extreme" if (af < 0.1 or af > 0.9) else
        ("Near-half" if 0.4 <= af <= 0.6 else "Intermediate")
    )
    return ng2


def plot_4panel(hcc: pd.DataFrame, others: pd.DataFrame, out_png: Path) -> dict:
    """4-panel: (A) AF dist; (B) CovM dist; (C) chr stacked TP/FP;
    (D) FP fraction by chr for HCC1954 vs cohort."""
    stats = {}

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        "B2 · HCC1954 Outer × NG=2 × cross_het TP outlier root cause\n"
        "(HCC1954 Outer cross-het TP = 0.08; cohort 0.24-0.88)",
        fontsize=14, fontweight="bold", y=0.98,
    )

    # (A) AF distribution
    ax = axes[0, 0]
    bins = np.linspace(0, 1, 41)
    ax.hist(hcc["AlleleDelta"].dropna(), bins=bins, density=True,
            alpha=0.6, color="#C1272D", label=f"HCC1954 (n={len(hcc)})")
    ax.hist(others["AlleleDelta"].dropna(), bins=bins, density=True,
            alpha=0.4, color="#4A6FA5", label=f"Other 5 TO (n={len(others)})")
    ax.set_xlabel("AlleleDelta (|AF|)")
    ax.set_ylabel("density")
    ax.set_title("(A) AF distribution — Outer × NG=2 × cross_het")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    hcc_af_low = float((hcc["AlleleDelta"] < 0.1).mean()) if len(hcc) else np.nan
    oth_af_low = float((others["AlleleDelta"] < 0.1).mean()) if len(others) else np.nan
    stats["hcc1954_af_low_frac"] = hcc_af_low
    stats["cohort_af_low_frac"] = oth_af_low
    ax.text(0.98, 0.95,
            f"AF<0.1 frac:\nHCC1954 = {hcc_af_low:.2%}\nCohort = {oth_af_low:.2%}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    # (B) CovM distribution
    ax = axes[0, 1]
    cov_bins = np.linspace(0, 5, 51)
    hcc_cov = hcc["Coverage_Multiple"].dropna().clip(0, 5)
    oth_cov = others["Coverage_Multiple"].dropna().clip(0, 5)
    ax.hist(hcc_cov, bins=cov_bins, density=True, alpha=0.6, color="#C1272D",
            label=f"HCC1954 (median={hcc_cov.median():.2f})")
    ax.hist(oth_cov, bins=cov_bins, density=True, alpha=0.4, color="#4A6FA5",
            label=f"Other 5 TO (median={oth_cov.median():.2f})")
    ax.set_xlabel("Coverage_Multiple (clipped to 5)")
    ax.set_ylabel("density")
    ax.set_title("(B) CovM distribution — Outer × NG=2 × cross_het")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    stats["hcc1954_covm_median"] = float(hcc_cov.median()) if len(hcc_cov) else np.nan
    stats["cohort_covm_median"] = float(oth_cov.median()) if len(oth_cov) else np.nan

    # (C) chr distribution: HCC1954 TP vs FP stacked, normalized
    ax = axes[1, 0]
    chrom_order = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    hcc_chr = hcc.groupby(["Chr", "tp_label"]).size().unstack(fill_value=0)
    hcc_chr = hcc_chr.reindex(chrom_order, fill_value=0)
    hcc_chr.columns = ["FP" if c == 0 else "TP" for c in hcc_chr.columns]
    for col in ["TP", "FP"]:
        if col not in hcc_chr.columns:
            hcc_chr[col] = 0
    x = np.arange(len(chrom_order))
    ax.bar(x, hcc_chr["TP"], color="#2E7D5B", label="TP", edgecolor="black", linewidth=0.3)
    ax.bar(x, hcc_chr["FP"], bottom=hcc_chr["TP"], color="#C26B3A", label="FP",
           edgecolor="black", linewidth=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order], fontsize=8, rotation=0)
    ax.set_ylabel("count")
    ax.set_title("(C) HCC1954 Outer × NG=2 × cross_het — TP/FP count by chr")
    # Highlight chr8 (MYC) and chr17 (HER2)
    for highlight_chr, _ in [("chr8", "MYC"), ("chr17", "HER2")]:
        if highlight_chr in chrom_order:
            idx = chrom_order.index(highlight_chr)
            ax.axvspan(idx - 0.5, idx + 0.5, color="gold", alpha=0.2)
    ax.legend(fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)

    # (D) FP fraction by chr: HCC1954 vs cohort
    ax = axes[1, 1]
    def fp_by_chr(df):
        g = df.groupby("Chr").agg(n=("tp_label", "size"),
                                  n_tp=("tp_label", "sum"))
        g["fp_frac"] = 1 - g["n_tp"] / g["n"]
        return g
    hcc_fpc = fp_by_chr(hcc).reindex(chrom_order)
    oth_fpc = fp_by_chr(others).reindex(chrom_order)
    width = 0.4
    ax.bar(x - width / 2, hcc_fpc["fp_frac"].fillna(0), width,
           color="#C1272D", alpha=0.85, label="HCC1954 FP frac",
           edgecolor="black", linewidth=0.3)
    ax.bar(x + width / 2, oth_fpc["fp_frac"].fillna(0), width,
           color="#4A6FA5", alpha=0.7, label="Cohort FP frac",
           edgecolor="black", linewidth=0.3)
    for highlight_chr in ["chr8", "chr17"]:
        if highlight_chr in chrom_order:
            idx = chrom_order.index(highlight_chr)
            ax.axvspan(idx - 0.5, idx + 0.5, color="gold", alpha=0.2)
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order], fontsize=8, rotation=0)
    ax.set_ylabel("FP fraction")
    ax.set_ylim(0, 1.05)
    ax.set_title("(D) FP fraction per chr — HCC1954 vs cohort (gold = MYC chr8, HER2 chr17)")
    ax.legend(fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)

    # Stats for (C)(D)
    total_fp = int(hcc_chr["FP"].sum())
    total_n = int(hcc_chr["TP"].sum() + hcc_chr["FP"].sum())
    stats["hcc1954_total_fp"] = total_fp
    stats["hcc1954_total_n"] = total_n
    stats["hcc1954_overall_fp_frac"] = total_fp / max(total_n, 1)
    top5_chr = hcc_chr["FP"].sort_values(ascending=False).head(5)
    stats["hcc1954_top5_fp_chrs"] = top5_chr.to_dict()
    # chr8 & chr17 burden
    for ch in ["chr8", "chr17"]:
        if ch in hcc_chr.index:
            stats[f"hcc1954_{ch}_fp_count"] = int(hcc_chr.loc[ch, "FP"])
            stats[f"hcc1954_{ch}_tp_count"] = int(hcc_chr.loc[ch, "TP"])
            stats[f"hcc1954_{ch}_fp_share_of_total_fp"] = int(hcc_chr.loc[ch, "FP"]) / max(total_fp, 1)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return stats


def main() -> None:
    print("[B2] loading 6 TO samples...")
    frames = []
    for name, paths in SAMPLES_TO.items():
        df = load_sample(name, paths)
        if not df.empty:
            df = prep(df)
            frames.append(df)
            print(f"  {name}: {len(df):,} NG=2 rows")
    if not frames:
        print("[B2] no data")
        sys.exit(1)

    all_ng2 = pd.concat(frames, ignore_index=True)

    # Filter Outer × cross-het (include both cross_het and cross_het_inv) × Extreme AF
    target_combos = {"cross_het (HP1 + HP2-1)", "cross_het_inv (HP1-1 + HP2)"}
    mask = (
        (all_ng2["loh_side"] == "Outer") &
        (all_ng2["combo"].isin(target_combos)) &
        (all_ng2["AF_class"] == "Extreme")
    )
    subset = all_ng2[mask].copy()
    print(f"[B2] Outer × cross_het × Extreme AF subset: {len(subset):,} rows")

    hcc = subset[subset["sample"] == "HCC1954"].copy()
    others = subset[subset["sample"] != "HCC1954"].copy()
    print(f"  HCC1954: {len(hcc):,}  | others: {len(others):,}")

    # TP rate confirmation
    def tp_rate(df):
        return df["tp_label"].mean() if len(df) else np.nan
    stats_per_sample = (subset.groupby("sample")
                        .agg(n=("tp_label", "size"),
                             n_tp=("tp_label", "sum"))
                        .assign(tp_rate=lambda d: d["n_tp"] / d["n"]))
    print("\n=== TP rate by sample (Outer × NG=2 × cross_het × Extreme AF) ===")
    print(stats_per_sample.to_string())
    stats_per_sample.to_csv(OUT_TSV_DIR / "B2_tp_rate_by_sample.tsv", sep="\t")

    # Main 4-panel figure
    fig_path = OUT_FIG_DIR / "B2_hcc1954_outlier_4panel.png"
    panel_stats = plot_4panel(hcc, others, fig_path)
    print(f"\n[B2] wrote {fig_path}")

    # Per-chr detailed table
    chrom_order = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    per_chr_rows = []
    for ch in chrom_order:
        hcc_ch = hcc[hcc["Chr"] == ch]
        oth_ch = others[others["Chr"] == ch]
        per_chr_rows.append({
            "chr": ch,
            "hcc1954_n": len(hcc_ch),
            "hcc1954_tp": int(hcc_ch["tp_label"].sum()),
            "hcc1954_fp": int((1 - hcc_ch["tp_label"]).sum()) if len(hcc_ch) else 0,
            "hcc1954_tp_rate": float(hcc_ch["tp_label"].mean()) if len(hcc_ch) else np.nan,
            "cohort_n": len(oth_ch),
            "cohort_tp": int(oth_ch["tp_label"].sum()) if len(oth_ch) else 0,
            "cohort_tp_rate": float(oth_ch["tp_label"].mean()) if len(oth_ch) else np.nan,
        })
    per_chr_df = pd.DataFrame(per_chr_rows)
    per_chr_df.to_csv(OUT_TSV_DIR / "B2_per_chr_breakdown.tsv", sep="\t", index=False)
    print(f"[B2] wrote {OUT_TSV_DIR / 'B2_per_chr_breakdown.tsv'}")

    # Summary stats
    print("\n=== Panel stats ===")
    for k, v in panel_stats.items():
        print(f"  {k}: {v}")

    # Additional H1 check: Potential_LOH annotation consistency
    # For NG=2 cross-het, obs18 says Outer should dominate. For HCC1954 specifically,
    # check fraction of cross_het with Potential_LOH=True vs False across samples.
    print("\n=== H1 Potential_LOH annotation diagnostic (cross_het × Extreme AF) ===")
    ch_subset = all_ng2[
        (all_ng2["combo"].isin(target_combos)) &
        (all_ng2["AF_class"] == "Extreme")
    ]
    h1_tab = (ch_subset.groupby(["sample", "loh_side"])
              .size().unstack(fill_value=0))
    if "Inner" not in h1_tab.columns:
        h1_tab["Inner"] = 0
    h1_tab["Inner_frac"] = h1_tab["Inner"] / (h1_tab["Inner"] + h1_tab.get("Outer", 0))
    print(h1_tab.to_string())
    h1_tab.to_csv(OUT_TSV_DIR / "B2_h1_loh_annotation.tsv", sep="\t")

    # Dump summary JSON-ish TSV
    import json
    summary = {
        "per_sample_tp_rate": stats_per_sample.to_dict(),
        "panel_stats": panel_stats,
    }
    with open(OUT_TSV_DIR / "B2_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n[B2] wrote {OUT_TSV_DIR / 'B2_summary.json'}")


if __name__ == "__main__":
    main()
