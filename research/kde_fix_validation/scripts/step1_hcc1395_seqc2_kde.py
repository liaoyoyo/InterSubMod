#!/usr/bin/env python3
"""Step 1 (KDE-fix rerun): HCC1395 Coverage_Multiple vs SEQC2 CNV truth.

Adapted from research/coverage_multiple_validation/scripts/step1_hcc1395_covm_vs_seqc2_truth.py
Changes:
  - Input: new master (all_region_rows_kde_B_tp.tsv.gz) — TP-only, per-sample KDE baseline
  - Modes: paired_pileup, paired_full (14 combos had no TO rerun; HCC1395 TO was covered in pilot)
  - Output: research/kde_fix_validation/outputs/step1_hcc1395_seqc2/{mode}/

Produces two mode-specific confusion matrices + scatter + per-CN-bin tables,
plus a summary comparing stale (14.6% Gain recall) vs fixed.
"""
import os
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/kde_rerun_B_14combos/all_region_rows_kde_B_tp.tsv.gz"
SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_LOSS_LOH_BED = f"{SEQC2_DIR}/ngs_benchmark_cnvs_gain_loss_loh.bed"
GAIN_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_gain_cn.bed"
LOSS_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_loss_cn.bed"
EXCLUSION_BED = f"{SEQC2_DIR}/exclusion.bed"

OUT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/kde_fix_validation/outputs/step1_hcc1395_seqc2")
STALE_GAIN_RECALL = 0.146  # Stale-binary baseline (75.0x) reference

COLS = ["Chr", "Pos", "NumReads", "Coverage_Multiple",
        "Diploid_Coverage_Used", "Coverage_Category",
        "sample", "mode", "truth_label"]


def load_bed(path, has_cn=False):
    return pd.read_csv(
        path, sep="\t", header=None,
        names=["chr", "start", "end", "val"] if has_cn else ["chr", "start", "end"],
        dtype={"chr": str},
    )


def build_idx(bed_df, extra_col=None):
    idx = {}
    for chrom, g in bed_df.groupby("chr"):
        g = g.sort_values("start").reset_index(drop=True)
        idx[chrom] = (g["start"].values, g["end"].values,
                      g[extra_col].values if extra_col else None)
    return idx


def lookup(idx, chrom, pos):
    if chrom not in idx:
        return None
    s, e, v = idx[chrom]
    j = np.searchsorted(s, pos, side="right") - 1
    if j < 0:
        return None
    if pos < e[j]:
        return v[j] if v is not None else True
    return None


def annotate(df, gain_cn_idx, loss_cn_idx, loh_idx, excl_idx):
    labels, cns = [], []
    for chrom, pos in zip(df["Chr"].values, df["Pos"].values):
        if pd.isna(pos):
            labels.append("Unknown"); cns.append(np.nan); continue
        pos = int(pos)
        if lookup(excl_idx, chrom, pos) is not None:
            labels.append("Excluded"); cns.append(np.nan); continue
        g = lookup(gain_cn_idx, chrom, pos)
        if g is not None:
            labels.append("Gain"); cns.append(float(g)); continue
        l = lookup(loss_cn_idx, chrom, pos)
        if l is not None:
            labels.append("Loss"); cns.append(float(l)); continue
        if lookup(loh_idx, chrom, pos) is not None:
            labels.append("LOH"); cns.append(2.0); continue
        labels.append("Neutral"); cns.append(2.0)
    return labels, cns


def analyze_mode(df, mode, gain_idx, loss_idx, loh_idx, excl_idx):
    sub = df[df["mode"] == mode].copy()
    if sub.empty:
        return None
    print(f"\n── Mode: {mode}   n_rows={len(sub):,}   baseline={sub['Diploid_Coverage_Used'].iloc[0]:.0f}x ──")

    labels, cns = annotate(sub, gain_idx, loss_idx, loh_idx, excl_idx)
    sub["truth_label_seqc2"] = labels
    sub["truth_cn"] = cns
    print(f"  Label dist: {sub['truth_label_seqc2'].value_counts().to_dict()}")

    out_dir = OUT_ROOT / mode
    out_dir.mkdir(parents=True, exist_ok=True)
    sub.to_csv(out_dir / "per_region_truth_aligned.tsv.gz", sep="\t",
               index=False, compression="gzip")

    # Confusion matrix
    keep = sub[sub["truth_label_seqc2"] != "Excluded"]
    pred_order = ["CNV_Loss", "Low", "Normal", "Elevated", "CNV_Gain", "High_Copy"]
    truth_order = ["Loss", "LOH", "Neutral", "Gain"]
    conf = pd.crosstab(keep["Coverage_Category"].fillna("Unknown"),
                       keep["truth_label_seqc2"],
                       margins=True, margins_name="Total")
    conf = conf.reindex(index=[c for c in pred_order if c in conf.index] + ["Total"],
                        columns=[c for c in truth_order if c in conf.columns] + ["Total"],
                        fill_value=0)
    conf.to_csv(out_dir / "confusion_matrix.tsv", sep="\t")

    # Per-CN-bin metrics
    bins = [(0, 1.0, "CN=0-1"), (1.0, 2.0, "CN=1-2"),
            (2.0, 2.01, "CN=2 (neutral+LOH)"),
            (2.01, 3.5, "CN=3"), (3.5, 5.0, "CN=4"),
            (5.0, 999, "CN≥5")]
    rows = []
    for lo, hi, lab in bins:
        if "2 (neutral" in lab:
            s = keep[keep["truth_cn"] == 2.0]
        else:
            s = keep[(keep["truth_cn"] >= lo) & (keep["truth_cn"] < hi)]
        if s.empty:
            continue
        rows.append({
            "cn_bin": lab, "n": len(s),
            "mean_CovM": s["Coverage_Multiple"].mean(),
            "median_CovM": s["Coverage_Multiple"].median(),
            "pred_gain_n": s["Coverage_Category"].isin(["CNV_Gain", "High_Copy"]).sum(),
            "pred_loss_n": (s["Coverage_Category"] == "CNV_Loss").sum(),
            "pred_norm_n": s["Coverage_Category"].isin(["Normal", "Elevated", "Low"]).sum(),
            "pred_gain_rate": s["Coverage_Category"].isin(["CNV_Gain", "High_Copy"]).sum() / len(s),
            "pred_loss_rate": (s["Coverage_Category"] == "CNV_Loss").sum() / len(s),
        })
    cn_metrics = pd.DataFrame(rows)
    cn_metrics.to_csv(out_dir / "per_cn_bin_metrics.tsv", sep="\t", index=False)

    # Overall recall
    gain_n = (keep["truth_label_seqc2"] == "Gain").sum()
    loss_n = (keep["truth_label_seqc2"] == "Loss").sum()
    neutral_n = (keep["truth_label_seqc2"] == "Neutral").sum()
    gain_tp = ((keep["truth_label_seqc2"] == "Gain") &
               keep["Coverage_Category"].isin(["CNV_Gain", "High_Copy"])).sum()
    loss_tp = ((keep["truth_label_seqc2"] == "Loss") &
               (keep["Coverage_Category"] == "CNV_Loss")).sum()
    neutral_tp = ((keep["truth_label_seqc2"] == "Neutral") &
                  keep["Coverage_Category"].isin(["Normal", "Elevated", "Low"])).sum()

    gain_rec = gain_tp / max(gain_n, 1)
    loss_rec = loss_tp / max(loss_n, 1)
    neu_prec = neutral_tp / max(keep["Coverage_Category"].isin(["Normal", "Elevated", "Low"]).sum(), 1)

    # Correlation
    full = keep.dropna(subset=["Coverage_Multiple", "truth_cn"])
    sr, sp = stats.spearmanr(full["Coverage_Multiple"], full["truth_cn"]) if len(full) > 10 else (np.nan, np.nan)
    pr, pp = stats.pearsonr(full["Coverage_Multiple"], full["truth_cn"]) if len(full) > 10 else (np.nan, np.nan)

    summary = {
        "mode": mode,
        "n_rows": len(sub),
        "baseline_x": float(sub["Diploid_Coverage_Used"].iloc[0]),
        "gain_n": int(gain_n),
        "loss_n": int(loss_n),
        "neutral_n": int(neutral_n),
        "gain_recall": float(gain_rec),
        "loss_recall": float(loss_rec),
        "neutral_precision": float(neu_prec),
        "spearman_r": float(sr),
        "pearson_r": float(pr),
        "stale_gain_recall": STALE_GAIN_RECALL,
        "gain_recall_ratio_vs_stale": float(gain_rec / STALE_GAIN_RECALL) if STALE_GAIN_RECALL > 0 else np.nan,
    }
    pd.DataFrame([summary]).to_csv(out_dir / "summary.tsv", sep="\t", index=False)

    print(f"  Gain recall: {gain_rec:.4f}  ({gain_tp}/{gain_n})"
          f"   Stale ref: {STALE_GAIN_RECALL:.4f}   ratio: {gain_rec/STALE_GAIN_RECALL:.2f}x")
    print(f"  Loss recall: {loss_rec:.4f}  ({loss_tp}/{loss_n})")
    print(f"  Neutral prec: {neu_prec:.4f}")
    print(f"  Spearman(CovM, CN): r={sr:.3f}  p={sp:.2e}")
    print(f"  → {out_dir}")

    # Figure: scatter + CN-bin bar
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    colors = {"Gain": "#e74c3c", "Loss": "#3498db",
              "LOH": "#f39c12", "Neutral": "#7f8c8d"}
    for lab, c in colors.items():
        s = keep[keep["truth_label_seqc2"] == lab].dropna(subset=["Coverage_Multiple", "truth_cn"])
        if s.empty:
            continue
        jitter = np.random.uniform(-0.08, 0.08, size=len(s))
        ax.scatter(s["truth_cn"].values + jitter, s["Coverage_Multiple"].values,
                   c=c, alpha=0.3, s=5, label=f"{lab} (n={len(s)})")
    ax.axhline(1.0, color="black", linestyle="--", alpha=0.4, label="CovM=1")
    ax.set_xlabel("SEQC2 truth CN"); ax.set_ylabel("Coverage_Multiple")
    ax.set_title(f"HCC1395 {mode} | baseline={summary['baseline_x']:.0f}x\n"
                 f"Gain recall={gain_rec:.3f} ({gain_rec/STALE_GAIN_RECALL:.1f}x stale)")
    ax.set_ylim(0, 4); ax.set_xlim(-0.3, 6)
    ax.legend(loc="upper left", fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    x = np.arange(len(cn_metrics))
    ax.bar(x, cn_metrics["mean_CovM"], color="#3498db", alpha=0.7, label="mean CovM")
    ax2 = ax.twinx()
    ax2.plot(x, cn_metrics["pred_gain_rate"], "ro-", label="pred-Gain rate")
    ax2.plot(x, cn_metrics["pred_loss_rate"], "bs-", label="pred-Loss rate")
    ax.set_xticks(x); ax.set_xticklabels(cn_metrics["cn_bin"], rotation=30, ha="right")
    ax.set_ylabel("mean Coverage_Multiple"); ax2.set_ylabel("pred rate")
    ax.set_title(f"Per-CN-bin | Spearman r={sr:.3f}")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "step1_summary.png", dpi=140, bbox_inches="tight")
    plt.close()

    return summary


def main():
    print("=" * 70)
    print("Step 1 KDE-FIX: HCC1395 Coverage_Multiple vs SEQC2 CNV truth")
    print("=" * 70)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    print(f"\nLoading {MASTER}...")
    chunks = []
    for chk in pd.read_csv(MASTER, sep="\t", usecols=COLS, chunksize=100_000, low_memory=False):
        m = chk["sample"].isin(["HCC1395", "HCC1395_DORADO"])
        if m.any():
            chunks.append(chk[m])
    df = pd.concat(chunks, ignore_index=True)
    for col in ["NumReads", "Coverage_Multiple", "Diploid_Coverage_Used", "Pos"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["Chr"] = df["Chr"].astype(str)
    print(f"  Loaded {len(df):,} HCC1395 rows")
    print(f"  By sample × mode:\n{df.groupby(['sample','mode']).size().to_string()}")

    # Filter to HCC1395 (not DORADO) for SEQC2 comparison — DORADO is same tumor, different basecaller
    df1395 = df[df["sample"] == "HCC1395"].copy()
    print(f"\n  HCC1395 only: {len(df1395):,} rows")

    print("\nLoading SEQC2 BEDs...")
    gll = load_bed(GAIN_LOSS_LOH_BED, has_cn=True)
    gcn = load_bed(GAIN_CN_BED, has_cn=True); gcn["val"] = gcn["val"].astype(float)
    lcn = load_bed(LOSS_CN_BED, has_cn=True); lcn["val"] = lcn["val"].astype(float)
    excl = load_bed(EXCLUSION_BED)

    gain_idx = build_idx(gcn, "val")
    loss_idx = build_idx(lcn, "val")
    loh_idx = build_idx(gll[gll["val"] == "loh"], "val")
    excl_idx = build_idx(excl)

    summaries = []
    for mode in ["paired_pileup", "paired_full"]:
        s = analyze_mode(df1395, mode, gain_idx, loss_idx, loh_idx, excl_idx)
        if s:
            summaries.append(s)

    if summaries:
        all_sum = pd.DataFrame(summaries)
        all_sum.to_csv(OUT_ROOT / "summary_all_modes.tsv", sep="\t", index=False)
        print("\n" + "=" * 70)
        print("ALL-MODE SUMMARY")
        print("=" * 70)
        print(all_sum.to_string(index=False))

        print("\n" + "=" * 70)
        print("H-CN1 VERDICT (based on new KDE baseline)")
        print("=" * 70)
        best_gain = max(s["gain_recall"] for s in summaries)
        best_spearman = max(s["spearman_r"] for s in summaries)
        if best_gain >= 0.6 and best_spearman >= 0.5:
            verdict = "FULL POSITIVE"
        elif best_gain >= 0.3 and best_spearman >= 0.3:
            verdict = "PARTIAL POSITIVE"
        else:
            verdict = "CONDITIONAL NEGATIVE"
        print(f"  Best Gain recall: {best_gain:.3f}  (stale {STALE_GAIN_RECALL:.3f}, ratio {best_gain/STALE_GAIN_RECALL:.1f}x)")
        print(f"  Best Spearman:    {best_spearman:.3f}")
        print(f"  VERDICT:          {verdict}")


if __name__ == "__main__":
    main()
