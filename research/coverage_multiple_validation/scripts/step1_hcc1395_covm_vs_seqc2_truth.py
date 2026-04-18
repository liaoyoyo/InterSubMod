#!/usr/bin/env python3
"""Step 1: HCC1395 Coverage_Multiple vs SEQC2 CNV truth alignment.

Context: Z3 amplicon blacklist pilot showed Coverage_Multiple (KDE-estimated
    diploid baseline) as CN proxy may be imprecise. Before concluding this is
    the root cause of HCC1954 cross-sample collateral damage, validate on
    HCC1395 where SEQC2 gold standard CNV truth is available.

Inputs:
  - Master dataset (HCC1395 + HCC1395_DORADO rows)
  - SEQC2 CNV truth BEDs: gain/loss/LOH + exclusion

Outputs:
  - data/step1_per_region_truth_aligned.tsv.gz  (per region × truth label/CN)
  - data/step1_confusion_matrix.tsv              (pred category × truth label)
  - data/step1_per_cn_bin_metrics.tsv            (per truth CN bin metrics)
  - data/step1_correlation_per_label.tsv          (Pearson/Spearman per subset)
  - figures/step1_confusion_matrix.png
  - figures/step1_covm_vs_truth_scatter.png
  - figures/step1_recall_per_cn.png
"""
import os
import sys
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')

MASTER = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
          "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_LOSS_LOH_BED = os.path.join(SEQC2_DIR, "ngs_benchmark_cnvs_gain_loss_loh.bed")
GAIN_CN_BED = os.path.join(SEQC2_DIR, "ngs_benchmark_cnv_gain_cn.bed")
LOSS_CN_BED = os.path.join(SEQC2_DIR, "ngs_benchmark_cnv_loss_cn.bed")
EXCLUSION_BED = os.path.join(SEQC2_DIR, "exclusion.bed")

OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/coverage_multiple_validation"
DATA_DIR = os.path.join(OUT_DIR, "data")
FIG_DIR = os.path.join(OUT_DIR, "figures")

COLS = ["Chr", "Pos", "NumReads", "Coverage_Multiple", "Coverage_Category",
        "sample", "mode", "truth_label"]


def load_bed(path, has_cn=False):
    df = pd.read_csv(path, sep='\t', header=None,
                     names=["chr", "start", "end", "val"] if has_cn else ["chr", "start", "end"],
                     dtype={"chr": str})
    return df


def build_interval_index(bed_df, extra_col=None):
    """Return dict chr → (starts_sorted, ends, vals_or_None)."""
    idx = {}
    for chrom, g in bed_df.groupby("chr"):
        g = g.sort_values("start").reset_index(drop=True)
        starts = g["start"].values
        ends = g["end"].values
        vals = g[extra_col].values if extra_col else None
        idx[chrom] = (starts, ends, vals)
    return idx


def lookup(idx, chrom, pos):
    """Return val (or True/False for membership) if pos in any interval."""
    if chrom not in idx:
        return None
    starts, ends, vals = idx[chrom]
    j = np.searchsorted(starts, pos, side='right') - 1
    if j < 0:
        return None
    if pos < ends[j]:
        return vals[j] if vals is not None else True
    return None


def main():
    print("=" * 70)
    print("Step 1: HCC1395 Coverage_Multiple vs SEQC2 CNV truth")
    print("=" * 70)

    # ── Load master dataset (HCC1395 + HCC1395_DORADO) ────────────────
    print("\nLoading master dataset...")
    chunks = []
    for chunk in pd.read_csv(MASTER, sep='\t', usecols=COLS,
                             chunksize=100_000, low_memory=False):
        mask = chunk["sample"].isin(["HCC1395", "HCC1395_DORADO"])
        if mask.any():
            chunks.append(chunk[mask].copy())
    df = pd.concat(chunks, ignore_index=True)
    print(f"  Loaded {len(df):,} HCC1395+DORADO rows")
    for col in ["NumReads", "Coverage_Multiple", "Pos"]:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df["Chr"] = df["Chr"].astype(str)
    print(f"  Samples × mode counts:")
    print(df.groupby(["sample", "mode"]).size().to_string())

    # ── Load SEQC2 truth BEDs ──────────────────────────────────────────
    print("\nLoading SEQC2 CNV truth BEDs...")
    gain_loss_loh = load_bed(GAIN_LOSS_LOH_BED, has_cn=True)
    gain_cn = load_bed(GAIN_CN_BED, has_cn=True)
    loss_cn = load_bed(LOSS_CN_BED, has_cn=True)
    excl = load_bed(EXCLUSION_BED)

    print(f"  Combined gain_loss_loh: {len(gain_loss_loh)} segments")
    print(f"    counts: {gain_loss_loh['val'].value_counts().to_dict()}")
    print(f"  Gain CN bed: {len(gain_cn)} segments  CN dist: "
          f"{gain_cn['val'].astype(float).value_counts().sort_index().to_dict()}")
    print(f"  Loss CN bed: {len(loss_cn)} segments  CN dist: "
          f"{loss_cn['val'].astype(float).value_counts().sort_index().to_dict()}")
    print(f"  Exclusion bed: {len(excl)} segments")

    # Cast CN cols to float
    gain_cn["val"] = gain_cn["val"].astype(float)
    loss_cn["val"] = loss_cn["val"].astype(float)

    # ── Build interval indexes ─────────────────────────────────────────
    excl_idx = build_interval_index(excl)
    gain_label_idx = build_interval_index(gain_loss_loh[gain_loss_loh["val"] == "gain"], extra_col="val")
    loss_label_idx = build_interval_index(gain_loss_loh[gain_loss_loh["val"] == "loss"], extra_col="val")
    loh_label_idx = build_interval_index(gain_loss_loh[gain_loss_loh["val"] == "loh"], extra_col="val")
    gain_cn_idx = build_interval_index(gain_cn, extra_col="val")
    loss_cn_idx = build_interval_index(loss_cn, extra_col="val")

    # ── Per-region annotation ──────────────────────────────────────────
    print("\nAnnotating per-region with SEQC2 truth labels...")
    labels = []
    truth_cns = []
    for chrom, pos in zip(df["Chr"].values, df["Pos"].values):
        if pd.isna(pos):
            labels.append("Unknown"); truth_cns.append(np.nan); continue
        pos = int(pos)
        if lookup(excl_idx, chrom, pos) is not None:
            labels.append("Excluded"); truth_cns.append(np.nan); continue
        g = lookup(gain_cn_idx, chrom, pos)
        if g is not None:
            labels.append("Gain"); truth_cns.append(float(g)); continue
        l = lookup(loss_cn_idx, chrom, pos)
        if l is not None:
            labels.append("Loss"); truth_cns.append(float(l)); continue
        if lookup(loh_label_idx, chrom, pos) is not None:
            labels.append("LOH"); truth_cns.append(2.0); continue
        labels.append("Neutral"); truth_cns.append(2.0)

    df["truth_label_seqc2"] = labels
    df["truth_cn"] = truth_cns
    print(f"  Label distribution:\n{df['truth_label_seqc2'].value_counts().to_string()}")

    # Save aligned table
    aligned_path = os.path.join(DATA_DIR, "step1_per_region_truth_aligned.tsv.gz")
    df.to_csv(aligned_path, sep='\t', index=False, compression='gzip')
    print(f"  → Saved {aligned_path}")

    # ── Confusion matrix (pred Coverage_Category × truth label) ───────
    print("\nBuilding confusion matrix (mode=to only)...")
    to_df = df[(df["mode"] == "to") & (df["truth_label_seqc2"] != "Excluded")].copy()
    pred_order = ["CNV_Loss", "Low", "Normal", "Elevated", "CNV_Gain", "High_Copy"]
    truth_order = ["Loss", "LOH", "Neutral", "Gain"]

    conf = pd.crosstab(
        to_df["Coverage_Category"].fillna("Unknown"),
        to_df["truth_label_seqc2"],
        margins=True, margins_name="Total",
    )
    # Reorder if present
    conf = conf.reindex(index=[c for c in pred_order if c in conf.index] + ["Total"],
                        columns=[c for c in truth_order if c in conf.columns] + ["Total"],
                        fill_value=0)
    conf.to_csv(os.path.join(DATA_DIR, "step1_confusion_matrix.tsv"), sep='\t')
    print(conf.to_string())

    # ── Per-CN-bin metrics ─────────────────────────────────────────────
    print("\nPer-CN-bin precision/recall for 'Gain' / 'Loss' calls...")
    # Define CN bins
    bins = [(0, 1.0, "CN=0-1"), (1.0, 2.0, "CN=1-2"),
            (2.0, 2.01, "CN=2 (neutral+LOH)"),
            (2.01, 3.5, "CN=3"), (3.5, 5.0, "CN=4"),
            (5.0, 999, "CN≥5")]
    rows = []
    for lo, hi, label in bins:
        if "2 (neutral" in label:
            sub = to_df[(to_df["truth_cn"] == 2.0)]
        else:
            sub = to_df[(to_df["truth_cn"] >= lo) & (to_df["truth_cn"] < hi)]
        if len(sub) == 0:
            continue
        pred_gain = sub["Coverage_Category"].isin(["CNV_Gain", "High_Copy"]).sum()
        pred_loss = sub["Coverage_Category"].isin(["CNV_Loss"]).sum()
        pred_norm = sub["Coverage_Category"].isin(["Normal", "Elevated", "Low"]).sum()
        mean_covm = sub["Coverage_Multiple"].mean()
        median_covm = sub["Coverage_Multiple"].median()
        rows.append({
            "cn_bin": label, "n": len(sub),
            "mean_CovM": mean_covm, "median_CovM": median_covm,
            "pred_gain_n": pred_gain, "pred_loss_n": pred_loss,
            "pred_norm_n": pred_norm,
            "pred_gain_rate": pred_gain / len(sub),
            "pred_loss_rate": pred_loss / len(sub),
        })
    cn_metrics = pd.DataFrame(rows)
    cn_metrics.to_csv(os.path.join(DATA_DIR, "step1_per_cn_bin_metrics.tsv"),
                      sep='\t', index=False)
    print(cn_metrics.to_string(index=False))

    # ── Recall per label ──────────────────────────────────────────────
    # Gain recall = pred_gain / total_gain
    gain_n = (to_df["truth_label_seqc2"] == "Gain").sum()
    loss_n = (to_df["truth_label_seqc2"] == "Loss").sum()
    neutral_n = (to_df["truth_label_seqc2"] == "Neutral").sum()
    loh_n = (to_df["truth_label_seqc2"] == "LOH").sum()

    gain_tp = ((to_df["truth_label_seqc2"] == "Gain")
               & to_df["Coverage_Category"].isin(["CNV_Gain", "High_Copy"])).sum()
    loss_tp = ((to_df["truth_label_seqc2"] == "Loss")
               & (to_df["Coverage_Category"] == "CNV_Loss")).sum()
    neutral_tp = ((to_df["truth_label_seqc2"] == "Neutral")
                  & to_df["Coverage_Category"].isin(["Normal", "Elevated", "Low"])).sum()

    gain_recall = gain_tp / max(gain_n, 1)
    loss_recall = loss_tp / max(loss_n, 1)
    neutral_precision = neutral_tp / max(
        to_df["Coverage_Category"].isin(["Normal", "Elevated", "Low"]).sum(), 1)

    print(f"\nHCC1395 (TO mode) metrics:")
    print(f"  Gain recall:       {gain_recall:.4f}  ({gain_tp}/{gain_n})")
    print(f"  Loss recall:       {loss_recall:.4f}  ({loss_tp}/{loss_n})")
    print(f"  Neutral precision: {neutral_precision:.4f}")
    print(f"  LOH total:         {loh_n}")

    # ── Correlation per label ─────────────────────────────────────────
    corr_rows = []
    for label in ["Gain", "Loss", "Neutral", "LOH"]:
        sub = to_df[(to_df["truth_label_seqc2"] == label)].dropna(
            subset=["Coverage_Multiple", "truth_cn"])
        if len(sub) < 10:
            continue
        pr, pp = stats.pearsonr(sub["Coverage_Multiple"], sub["truth_cn"])
        sr, sp = stats.spearmanr(sub["Coverage_Multiple"], sub["truth_cn"])
        corr_rows.append({
            "label": label, "n": len(sub),
            "pearson_r": pr, "pearson_p": pp,
            "spearman_r": sr, "spearman_p": sp,
            "mean_CovM": sub["Coverage_Multiple"].mean(),
            "mean_truth_cn": sub["truth_cn"].mean(),
        })
    # Global (mode=to, excluding Excluded/Unknown)
    sub_all = to_df.dropna(subset=["Coverage_Multiple", "truth_cn"])
    if len(sub_all) > 10:
        pr, pp = stats.pearsonr(sub_all["Coverage_Multiple"], sub_all["truth_cn"])
        sr, sp = stats.spearmanr(sub_all["Coverage_Multiple"], sub_all["truth_cn"])
        corr_rows.append({
            "label": "ALL (TO, excl. Excluded)", "n": len(sub_all),
            "pearson_r": pr, "pearson_p": pp,
            "spearman_r": sr, "spearman_p": sp,
            "mean_CovM": sub_all["Coverage_Multiple"].mean(),
            "mean_truth_cn": sub_all["truth_cn"].mean(),
        })
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(os.path.join(DATA_DIR, "step1_correlation_per_label.tsv"),
                   sep='\t', index=False)
    print("\nCorrelation per label:")
    print(corr_df.to_string(index=False))

    # ── Figure 1: Confusion matrix heatmap ────────────────────────────
    conf_plot = conf.drop(columns=["Total"]).drop(index=["Total"])
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(conf_plot.values, aspect='auto', cmap='Blues')
    ax.set_xticks(range(len(conf_plot.columns)))
    ax.set_xticklabels(conf_plot.columns, rotation=30, ha='right')
    ax.set_yticks(range(len(conf_plot.index)))
    ax.set_yticklabels(conf_plot.index)
    ax.set_xlabel("SEQC2 truth label")
    ax.set_ylabel("Coverage_Category (predicted)")
    ax.set_title("HCC1395 (TO) — Pred × Truth confusion matrix")
    for i in range(conf_plot.shape[0]):
        for j in range(conf_plot.shape[1]):
            v = conf_plot.values[i, j]
            txt_color = 'white' if v > conf_plot.values.max() * 0.5 else 'black'
            ax.text(j, i, f"{v}", ha='center', va='center',
                    color=txt_color, fontsize=9)
    plt.colorbar(im, ax=ax, label='count')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "step1_confusion_matrix.png"),
                dpi=150, bbox_inches='tight')
    plt.close()

    # ── Figure 2: CovM vs truth CN scatter ─────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = {"Gain": "#e74c3c", "Loss": "#3498db",
              "LOH": "#f39c12", "Neutral": "#7f8c8d"}
    for label, c in colors.items():
        sub = to_df[(to_df["truth_label_seqc2"] == label)].dropna(
            subset=["Coverage_Multiple", "truth_cn"])
        if len(sub) == 0:
            continue
        # Jitter truth_cn for visibility
        jitter = np.random.uniform(-0.08, 0.08, size=len(sub))
        ax.scatter(sub["truth_cn"].values + jitter,
                   sub["Coverage_Multiple"].values,
                   c=c, alpha=0.3, s=6, label=f"{label} (n={len(sub)})")
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.4, label='CovM=1 (pred CN=2)')
    ax.axvline(2.0, color='black', linestyle=':', alpha=0.4)
    ax.set_xlabel("SEQC2 truth CN")
    ax.set_ylabel("Coverage_Multiple (predicted CN proxy)")
    ax.set_title("HCC1395 (TO) — Coverage_Multiple vs SEQC2 truth CN")
    ax.set_ylim(0, 4)
    ax.set_xlim(-0.3, 6)
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "step1_covm_vs_truth_scatter.png"),
                dpi=150, bbox_inches='tight')
    plt.close()

    # ── Figure 3: Recall per CN bin ────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(cn_metrics))
    ax.bar(x, cn_metrics["mean_CovM"], alpha=0.7, color='#3498db',
           label='mean Coverage_Multiple')
    ax2 = ax.twinx()
    ax2.plot(x, cn_metrics["pred_gain_rate"], 'ro-', label='pred-Gain rate')
    ax2.plot(x, cn_metrics["pred_loss_rate"], 'bs-', label='pred-Loss rate')
    ax.set_xticks(x)
    ax.set_xticklabels(cn_metrics["cn_bin"], rotation=30, ha='right')
    ax.set_ylabel("Coverage_Multiple (mean)")
    ax2.set_ylabel("pred rate")
    ax.set_title(f"HCC1395 (TO) — CovM vs truth CN bin  "
                 f"| Gain recall={gain_recall:.2f}  Loss recall={loss_recall:.2f}")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "step1_recall_per_cn.png"),
                dpi=150, bbox_inches='tight')
    plt.close()

    # ── Summary stdout ─────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY (判定基準)")
    print("=" * 70)
    print(f"  Gain recall (H-CN2):  {gain_recall:.3f}  "
          f"→ {'POSITIVE (baseline 偏誤)' if gain_recall < 0.6 else 'NEGATIVE (proxy 足夠)' if gain_recall >= 0.8 else 'CONDITIONAL'}")
    print(f"  Loss recall:          {loss_recall:.3f}")
    print(f"  Gain-label CovM mean vs truth mean: "
          f"{cn_metrics[cn_metrics['cn_bin'].str.contains('3|4|≥5')]['mean_CovM'].mean():.3f} "
          f"vs expected {cn_metrics[cn_metrics['cn_bin'].str.contains('3|4|≥5')]['n'].sum()} Gain rows")
    print("=" * 70)


if __name__ == "__main__":
    main()
