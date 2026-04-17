#!/usr/bin/env python3
"""B.2-3: NR≥80 高功效 bin 單獨驗證 + NR-matched sampling.

Questions:
1. 原報告 NR bin 10-30 / 30-50 / 50-80 effect 遞增 (ρ=0.483→0.709)。
   是否 NR≥80 仍有顯著 AF-extreme vs intermediate Δ_NG？
2. NR-matched sampling 消除 NR 分佈差異後效應是否保留？

Data: all_region_rows.tsv.gz (TO + Paired, LOH only)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = OUT_DIR / "data"
FIG_DIR = OUT_DIR / "figures"

SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
NR_BINS = [(40, 60), (60, 80), (80, 100), (100, 150), (150, 500)]
MODE_TO = "to"  # Master dataset uses lowercase
MODE_PAIRED = "paired"


def load_master() -> pd.DataFrame:
    print(f"Loading {MASTER.name}...")
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)
    # LOH annotation only in TO mode (Paired mode has no to_loh_bed_hit)
    df = df[df["to_loh_bed_hit"] == True].copy()
    print(f"  LOH (TO only) available modes: {df['mode'].value_counts().to_dict()}")
    df["AF_class"] = pd.cut(
        df["caller_af"],
        bins=[-0.01, 0.3, 0.4, 0.6, 0.7, 1.01],
        labels=["ext_low", "near_ext", "intermediate", "near_ext2", "ext_high"],
    )
    df["AF_group"] = df["AF_class"].map({
        "ext_low": "extreme",
        "ext_high": "extreme",
        "intermediate": "intermediate",
        "near_ext": "near",
        "near_ext2": "near",
    })
    return df


def nr_bin_analysis(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    rows = []
    sub = df[df["mode"] == mode].copy()
    for lo, hi in NR_BINS:
        nr_sub = sub[(sub["NumReads"] >= lo) & (sub["NumReads"] < hi)]
        for sample in SAMPLES:
            s_sub = nr_sub[nr_sub["sample"] == sample]
            ext = s_sub[s_sub["AF_group"] == "extreme"]["HPFineNGroups"].dropna()
            inter = s_sub[s_sub["AF_group"] == "intermediate"]["HPFineNGroups"].dropna()
            if len(ext) < 20 or len(inter) < 20:
                rows.append({
                    "mode": mode, "sample": sample, "nr_lo": lo, "nr_hi": hi,
                    "n_ext": len(ext), "n_inter": len(inter),
                    "ng_ext": np.nan, "ng_inter": np.nan, "delta_ng": np.nan,
                    "mw_p": np.nan, "r": np.nan, "sig": False,
                })
                continue
            u, p = stats.mannwhitneyu(ext, inter, alternative="less")
            r = 1 - 2 * u / (len(ext) * len(inter))
            rows.append({
                "mode": mode, "sample": sample, "nr_lo": lo, "nr_hi": hi,
                "n_ext": len(ext), "n_inter": len(inter),
                "ng_ext": ext.mean(), "ng_inter": inter.mean(),
                "delta_ng": inter.mean() - ext.mean(),
                "mw_p": p, "r": r, "sig": p < 0.05,
            })
    return pd.DataFrame(rows)


def nr_matched_sampling(df: pd.DataFrame, mode: str, seed: int = 42) -> pd.DataFrame:
    """Match NR distribution via quantile binning + subsample."""
    rows = []
    rng = np.random.default_rng(seed)
    sub = df[(df["mode"] == mode) & (df["NumReads"] >= 40)].copy()
    for sample in SAMPLES:
        s = sub[sub["sample"] == sample]
        ext = s[s["AF_group"] == "extreme"].copy()
        inter = s[s["AF_group"] == "intermediate"].copy()
        if len(ext) < 50 or len(inter) < 50:
            rows.append({
                "mode": mode, "sample": sample, "n_matched": 0,
                "ng_ext_matched": np.nan, "ng_inter_matched": np.nan,
                "delta_matched": np.nan, "mw_p_matched": np.nan,
            })
            continue
        # Decile-based NR matching
        bins = np.quantile(s["NumReads"], np.linspace(0, 1, 11))
        bins = np.unique(bins)
        if len(bins) < 3:
            continue
        matched_ext = []
        matched_inter = []
        for b_lo, b_hi in zip(bins[:-1], bins[1:]):
            e = ext[(ext["NumReads"] >= b_lo) & (ext["NumReads"] < b_hi)]
            i = inter[(inter["NumReads"] >= b_lo) & (inter["NumReads"] < b_hi)]
            n_match = min(len(e), len(i))
            if n_match == 0:
                continue
            matched_ext.append(e.sample(n=n_match, random_state=rng.integers(1e9)))
            matched_inter.append(i.sample(n=n_match, random_state=rng.integers(1e9)))
        if not matched_ext:
            continue
        e_df = pd.concat(matched_ext)
        i_df = pd.concat(matched_inter)
        e_ng = e_df["HPFineNGroups"].dropna()
        i_ng = i_df["HPFineNGroups"].dropna()
        if len(e_ng) < 20 or len(i_ng) < 20:
            continue
        u, p = stats.mannwhitneyu(e_ng, i_ng, alternative="less")
        rows.append({
            "mode": mode, "sample": sample, "n_matched": len(e_ng),
            "ng_ext_matched": e_ng.mean(), "ng_inter_matched": i_ng.mean(),
            "delta_matched": i_ng.mean() - e_ng.mean(),
            "mw_p_matched": p,
        })
    return pd.DataFrame(rows)


def plot_nr_bin_forest(to_df: pd.DataFrame, paired_df: pd.DataFrame, path: Path):
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    bin_labels = [f"[{lo},{hi})" for lo, hi in NR_BINS]
    for ax, df, mode in zip(axes, [to_df, paired_df], ["TO", "Paired"]):
        if df.empty or df["delta_ng"].notna().sum() == 0:
            ax.text(0.5, 0.5, f"{mode}: no data (LOH annotation unavailable)",
                    ha="center", va="center", transform=ax.transAxes, fontsize=14)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        pivot = df.pivot_table(index="sample", columns=["nr_lo", "nr_hi"], values="delta_ng")
        pivot.columns = bin_labels[:len(pivot.columns)]
        pivot = pivot.reindex(SAMPLES)
        im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(bin_labels)))
        ax.set_xticklabels(bin_labels, rotation=30, ha="right")
        ax.set_yticks(range(len(SAMPLES)))
        ax.set_yticklabels(SAMPLES)
        ax.set_title(f"{mode} mode: Δ_NG = intermediate - extreme (LOH)")
        for i, sample in enumerate(SAMPLES):
            for j, col in enumerate(bin_labels):
                v = pivot.loc[sample, col]
                if pd.notna(v):
                    ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                            color="white" if abs(v) > 0.5 else "black", fontsize=9)
        plt.colorbar(im, ax=ax, label="Δ_NG")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close()


def main():
    print("=== B.2-3: High NR bin validation + NR-matched sampling ===\n")
    df = load_master()
    print(f"  LOH rows: {len(df):,}\n")

    # NR bin analysis (LOH annotation only in TO mode; Paired LOH unavailable)
    to_nr = nr_bin_analysis(df, MODE_TO)
    paired_nr = nr_bin_analysis(df, MODE_PAIRED)

    print("[TO NR bin summary]")
    for (lo, hi) in NR_BINS:
        sub = to_nr[(to_nr["nr_lo"] == lo) & (to_nr["sig"] == True)]
        print(f"  NR [{lo},{hi}): {len(sub)}/7 sig pos, median Δ_NG={to_nr[(to_nr['nr_lo']==lo) & to_nr['delta_ng'].notna()]['delta_ng'].median():.3f}")

    print("\n[Paired NR bin summary]")
    for (lo, hi) in NR_BINS:
        sub = paired_nr[(paired_nr["nr_lo"] == lo) & (paired_nr["sig"] == True)]
        print(f"  NR [{lo},{hi}): {len(sub)}/7 sig pos, median Δ_NG={paired_nr[(paired_nr['nr_lo']==lo) & paired_nr['delta_ng'].notna()]['delta_ng'].median():.3f}")

    to_nr.to_csv(DATA_DIR / "b2_3_nr_bin_to.tsv", sep="\t", index=False)
    paired_nr.to_csv(DATA_DIR / "b2_3_nr_bin_paired.tsv", sep="\t", index=False)

    # NR-matched sampling
    print("\n[NR-matched sampling]")
    to_match = nr_matched_sampling(df, MODE_TO)
    paired_match = nr_matched_sampling(df, MODE_PAIRED)
    print("\nTO NR-matched:")
    print(to_match.to_string(index=False))
    print("\nPaired NR-matched:")
    print(paired_match.to_string(index=False))

    to_match.to_csv(DATA_DIR / "b2_3_nr_matched_to.tsv", sep="\t", index=False)
    paired_match.to_csv(DATA_DIR / "b2_3_nr_matched_paired.tsv", sep="\t", index=False)

    plot_nr_bin_forest(to_nr, paired_nr, FIG_DIR / "01_nr_bin_heatmap.png")

    # FINAL VERDICTS
    print("\n=== FINAL VERDICTS ===")

    to_high = to_nr[to_nr["nr_lo"] >= 80]
    paired_high = paired_nr[paired_nr["nr_lo"] >= 80]
    to_high_sig = to_high[to_high["sig"] == True]
    paired_high_sig = paired_high[paired_high["sig"] == True]
    to_high_pos = to_high[(to_high["delta_ng"] > 0) & to_high["sig"]]
    paired_high_pos = paired_high[(paired_high["delta_ng"] > 0) & paired_high["sig"]]

    print(f"\n[H_B2_3-a High NR]")
    print(f"  TO NR≥80 sig pos: {len(to_high_pos)}/{to_high['delta_ng'].notna().sum()} rows")
    print(f"  Paired NR≥80 sig pos: {len(paired_high_pos)}/{paired_high['delta_ng'].notna().sum()} rows")

    print(f"\n[H_B2_3-b NR-matched]")
    to_match_pos = to_match[(to_match["delta_matched"] > 0) & (to_match["mw_p_matched"] < 0.05)]
    paired_match_pos = paired_match[(paired_match["delta_matched"] > 0) & (paired_match["mw_p_matched"] < 0.05)]
    print(f"  TO matched sig pos: {len(to_match_pos)}/{to_match['delta_matched'].notna().sum()} samples")
    print(f"  Paired matched sig pos: {len(paired_match_pos)}/{paired_match['delta_matched'].notna().sum()} samples")

    if len(to_high_pos) >= 4 and len(paired_high_pos) >= 4 and len(to_match_pos) >= 4 and len(paired_match_pos) >= 4:
        print("\n  → ROBUST: 高 NR bin + NR-matched 皆保留效應 ≥4/7")
    else:
        print("\n  → ATTENUATED or CONFOUNDED")


if __name__ == "__main__":
    main()
