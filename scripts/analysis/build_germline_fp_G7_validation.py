#!/usr/bin/env python3
"""G7: Cross-Sample Validation & Paired Ground Truth.

Validates all findings across 7 samples, cross-references with paired mode,
and simulates F1 impact of applying the germline filter.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    OUTPUT_ROOT, SAMPLE_ORDER, setup_plot_style, save_figure, add_caption,
    compute_auc, bootstrap_ci, ensure_dir,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
)

WORKSPACE = OUTPUT_ROOT / "20260401_germline_fp_identification"
G1_DIR = WORKSPACE / "G1_vcf_features"
G6_DIR = WORKSPACE / "G6_combination"
G7_DIR = ensure_dir(WORKSPACE / "G7_validation")
FIG_DIR = ensure_dir(G7_DIR / "figures")
DATA_DIR = ensure_dir(G7_DIR / "data")


def _print(msg: str):
    print(msg, flush=True)


def load_pass_features() -> pd.DataFrame:
    path = G1_DIR / "all_samples_merged.tsv.gz"
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    for col in ["af", "sb_clairs", "sb_ratio"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["is_het", "is_hom_alt", "h_flag", "is_cpg_ct", "is_cpg_site",
                "is_transition", "filter_is_pass"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    return df


def compute_combined_score(df: pd.DataFrame) -> np.ndarray:
    """Reproduce the Combined (A+B+C) score from G6."""
    af = df["af"].values

    af_germ = np.maximum(
        1 - np.abs(af - 0.5) * 4,
        1 - np.abs(af - 1.0) * 4,
    ).clip(0, 1)

    h_absent = 1 - df["h_flag"].astype(float).values
    hom_alt = df["is_hom_alt"].astype(float).values
    cpg_ct = df["is_cpg_ct"].astype(float).values
    af_near_05 = (np.abs(af - 0.5) < 0.1).astype(float)
    af_near_10 = (af > 0.85).astype(float)
    is_het = df["is_het"].astype(float).values

    rule_a = af_germ * 0.3 + h_absent * 0.25 + hom_alt * 0.25 + cpg_ct * 0.2
    rule_b = af_near_10 * 0.5 + hom_alt * 0.5
    rule_c = af_near_05 * 0.4 + is_het * 0.3 + h_absent * 0.3

    return np.maximum(rule_a, np.maximum(rule_b, rule_c))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def fig01_consistency_heatmap(df: pd.DataFrame):
    """Consistency heatmap across 7 samples."""
    _print("\n[G7] Fig01: Consistency heatmap")

    pass_df = df[df["filter_is_pass"] == 1].copy()

    metrics = {
        "mean_AF": lambda g: g["af"].mean(),
        "hom_alt_rate": lambda g: (g["is_hom_alt"] == 1).mean(),
        "h_flag_rate": lambda g: (g["h_flag"] == 1).mean(),
        "cpg_ct_rate": lambda g: (g["is_cpg_ct"] == 1).mean(),
        "ti_rate": lambda g: (g["is_transition"] == 1).mean(),
    }

    rows = []
    for sample in SAMPLE_ORDER:
        for label in ["TP", "FP"]:
            sub = pass_df[(pass_df["sample"] == sample) & (pass_df["truth_label"] == label)]
            if sub.empty:
                continue
            row = {"sample": sample, "truth_label": label}
            for mname, fn in metrics.items():
                try:
                    row[mname] = fn(sub)
                except Exception:
                    row[mname] = float("nan")
            rows.append(row)

    df_m = pd.DataFrame(rows)

    # Delta heatmap (FP - TP per sample)
    delta_rows = []
    for sample in SAMPLE_ORDER:
        tp = df_m[(df_m["sample"] == sample) & (df_m["truth_label"] == "TP")]
        fp = df_m[(df_m["sample"] == sample) & (df_m["truth_label"] == "FP")]
        if tp.empty or fp.empty:
            continue
        row = {"sample": sample}
        for mname in metrics:
            row[mname] = fp[mname].values[0] - tp[mname].values[0]
        delta_rows.append(row)

    df_delta = pd.DataFrame(delta_rows).set_index("sample")
    df_delta = df_delta.reindex(SAMPLE_ORDER).dropna(how="all")

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(df_delta, annot=True, fmt=".3f", cmap="RdBu_r", center=0, ax=ax)
    ax.set_title("G7-Fig01: Per-Sample FP-TP Delta Consistency",
                 fontsize=13, fontweight="bold")
    save_figure(fig, FIG_DIR / "G7_fig01_consistency_heatmap.png")

    df_delta.to_csv(DATA_DIR / "consistency_delta.tsv", sep="\t")


def fig03_f1_improvement(df: pd.DataFrame):
    """F1 improvement trajectory with threshold sweep."""
    _print("\n[G7] Fig03: F1 improvement")

    pass_df = df[df["filter_is_pass"] == 1].copy()
    score = compute_combined_score(pass_df)
    y = (pass_df["truth_label"] == "FP").astype(int).values
    samples = pass_df["sample"].values

    fig, ax = plt.subplots(figsize=(12, 6))

    for sample in SAMPLE_ORDER:
        mask = samples == sample
        if mask.sum() == 0:
            continue
        y_s = y[mask]
        s_s = score[mask]
        m_finite = np.isfinite(s_s)
        y_s = y_s[m_finite]
        s_s = s_s[m_finite]

        tp_total = (y_s == 0).sum()
        fp_total = (y_s == 1).sum()
        if tp_total == 0 or fp_total == 0:
            continue

        # Baseline F1 (no filter)
        base_precision = tp_total / (tp_total + fp_total)
        base_recall = 1.0
        base_f1 = 2 * base_precision * base_recall / (base_precision + base_recall)

        f1s = []
        thresholds = np.linspace(0, 1, 100)
        for t in thresholds:
            flagged = s_s >= t
            tp_lost = (flagged & (y_s == 0)).sum()
            fp_removed = (flagged & (y_s == 1)).sum()
            new_tp = tp_total - tp_lost
            new_fp = fp_total - fp_removed
            if new_tp + new_fp == 0:
                f1s.append(0)
                continue
            prec = new_tp / (new_tp + new_fp)
            rec = new_tp / tp_total
            f1 = 2 * prec * rec / (prec + rec) if prec + rec > 0 else 0
            f1s.append(f1)

        delta_f1 = np.array(f1s) - base_f1
        ax.plot(thresholds, delta_f1 * 100, label=sample, alpha=0.7, linewidth=1.5)

    ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Germline Score Threshold")
    ax.set_ylabel("ΔF1 (%)")
    ax.set_title("G7-Fig03: F1 Improvement by Threshold (per sample)",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    add_caption(fig, "Positive = F1 improves with germline filter | Negative = too many TP lost")
    save_figure(fig, FIG_DIR / "G7_fig03_f1_improvement.png")


def fig04_bootstrap_auc(df: pd.DataFrame):
    """Bootstrap AUC distribution."""
    _print("\n[G7] Fig04: Bootstrap AUC")

    pass_df = df[df["filter_is_pass"] == 1].copy()
    score = compute_combined_score(pass_df)
    y = (pass_df["truth_label"] == "FP").astype(int).values
    mask = np.isfinite(score)
    y_m, s_m = y[mask], score[mask]

    n_boot = 1000
    boot_aucs = []
    rng = np.random.RandomState(42)
    for _ in range(n_boot):
        idx = rng.choice(len(y_m), size=len(y_m), replace=True)
        auc = compute_auc(y_m[idx], s_m[idx])
        if not np.isnan(auc):
            boot_aucs.append(auc)

    boot_aucs = np.array(boot_aucs)
    ci_lo, ci_hi = np.percentile(boot_aucs, [2.5, 97.5])
    mean_auc = boot_aucs.mean()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(boot_aucs, bins=50, color="#E65100", alpha=0.7, edgecolor="white")
    ax.axvline(mean_auc, color="red", linewidth=2, label=f"Mean={mean_auc:.3f}")
    ax.axvline(ci_lo, color="red", linestyle="--", label=f"95% CI=[{ci_lo:.3f}, {ci_hi:.3f}]")
    ax.axvline(ci_hi, color="red", linestyle="--")
    ax.set_xlabel("AUC")
    ax.set_ylabel("Frequency")
    ax.set_title(f"G7-Fig04: Bootstrap AUC Distribution (n={n_boot})",
                 fontsize=13, fontweight="bold")
    ax.legend()
    save_figure(fig, FIG_DIR / "G7_fig04_bootstrap_auc.png")

    _print(f"  [G7] Bootstrap AUC: {mean_auc:.4f} [95% CI: {ci_lo:.4f}-{ci_hi:.4f}]")


def fig05_loso_auc(df: pd.DataFrame):
    """Leave-one-sample-out AUC."""
    _print("\n[G7] Fig05: LOSO AUC")

    pass_df = df[df["filter_is_pass"] == 1].copy()
    score = compute_combined_score(pass_df)
    y = (pass_df["truth_label"] == "FP").astype(int).values
    samples = pass_df["sample"].values

    aucs = []
    for test_sample in SAMPLE_ORDER:
        mask = (samples == test_sample) & np.isfinite(score)
        if mask.sum() < 50:
            continue
        auc = compute_auc(y[mask], score[mask])
        aucs.append({"sample": test_sample, "auc": auc})

    df_auc = pd.DataFrame(aucs)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df_auc["sample"], df_auc["auc"], color="#E65100", alpha=0.8)
    ax.axhline(df_auc["auc"].mean(), color="red", linestyle="--",
               label=f"Mean={df_auc['auc'].mean():.3f}")
    ax.set_ylabel("AUC")
    ax.set_title("G7-Fig05: Per-Sample Combined Score AUC",
                 fontsize=13, fontweight="bold")
    ax.legend()
    ax.set_ylim(0.4, 1.0)
    plt.xticks(rotation=30, ha="right")
    save_figure(fig, FIG_DIR / "G7_fig05_loso_auc.png")


def fig06_f1_before_after(df: pd.DataFrame):
    """F1 before vs after dashboard."""
    _print("\n[G7] Fig06: F1 before vs after")

    pass_df = df[df["filter_is_pass"] == 1].copy()
    score = compute_combined_score(pass_df)
    y = (pass_df["truth_label"] == "FP").astype(int).values
    samples = pass_df["sample"].values

    # Find global optimal threshold (2% TP loss)
    mask_g = np.isfinite(score)
    y_g, s_g = y[mask_g], score[mask_g]
    tp_total_g = (y_g == 0).sum()
    best_t = 0.5
    best_fp_rem = 0
    for t in np.linspace(0, 1, 500):
        flagged = s_g >= t
        tp_loss = (flagged & (y_g == 0)).sum() / tp_total_g
        fp_rem = (flagged & (y_g == 1)).sum()
        if tp_loss <= 0.02 and fp_rem > best_fp_rem:
            best_fp_rem = fp_rem
            best_t = t

    rows = []
    for sample in SAMPLE_ORDER:
        mask = (samples == sample) & np.isfinite(score)
        if mask.sum() == 0:
            continue
        y_s, s_s = y[mask], score[mask]
        tp = (y_s == 0).sum()
        fp = (y_s == 1).sum()

        flagged = s_s >= best_t
        tp_lost = (flagged & (y_s == 0)).sum()
        fp_removed = (flagged & (y_s == 1)).sum()

        prec_before = tp / (tp + fp) if tp + fp > 0 else 0
        f1_before = 2 * prec_before / (prec_before + 1)  # recall=1

        new_tp = tp - tp_lost
        new_fp = fp - fp_removed
        prec_after = new_tp / (new_tp + new_fp) if new_tp + new_fp > 0 else 0
        rec_after = new_tp / tp if tp > 0 else 0
        f1_after = 2 * prec_after * rec_after / (prec_after + rec_after) if prec_after + rec_after > 0 else 0

        rows.append({
            "sample": sample, "f1_before": f1_before, "f1_after": f1_after,
            "delta_f1": f1_after - f1_before,
            "tp_loss_pct": tp_lost / tp * 100 if tp > 0 else 0,
            "fp_removal_pct": fp_removed / fp * 100 if fp > 0 else 0,
        })

    df_f1 = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    x = np.arange(len(df_f1))
    ax.bar(x - 0.2, df_f1["f1_before"], 0.4, color=COLOR_TO, label="Before filter", alpha=0.8)
    ax.bar(x + 0.2, df_f1["f1_after"], 0.4, color="#4CAF50", label="After filter", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(df_f1["sample"], rotation=30, ha="right")
    ax.set_ylabel("F1 Score")
    ax.set_title(f"F1 Before/After (threshold={best_t:.3f})")
    ax.legend()

    ax = axes[1]
    colors = ["#4CAF50" if d > 0 else "#F44336" for d in df_f1["delta_f1"]]
    ax.bar(x, df_f1["delta_f1"] * 100, color=colors, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(df_f1["sample"], rotation=30, ha="right")
    ax.set_ylabel("ΔF1 (%)")
    ax.set_title("F1 Change")
    ax.axhline(0, color="gray", linewidth=0.5)

    fig.suptitle("G7-Fig06: F1 Before vs After Germline Filter",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G7_fig06_f1_before_after.png")

    df_f1.to_csv(DATA_DIR / "f1_comparison.tsv", sep="\t", index=False)
    _print(f"\n  [G7] F1 comparison (threshold={best_t:.3f}):")
    _print(df_f1.to_string(index=False))


def fig07_tp_fp_scatter(df: pd.DataFrame):
    """TP lost vs FP removed per sample."""
    _print("\n[G7] Fig07: TP vs FP scatter")

    pass_df = df[df["filter_is_pass"] == 1].copy()
    score = compute_combined_score(pass_df)
    y = (pass_df["truth_label"] == "FP").astype(int).values
    samples = pass_df["sample"].values

    fig, ax = plt.subplots(figsize=(10, 6))

    for t_val, marker, size in [(0.3, "o", 40), (0.5, "s", 60), (0.7, "^", 80), (0.9, "D", 100)]:
        for sample in SAMPLE_ORDER:
            mask = (samples == sample) & np.isfinite(score)
            if mask.sum() == 0:
                continue
            y_s, s_s = y[mask], score[mask]
            flagged = s_s >= t_val
            tp_total = (y_s == 0).sum()
            fp_total = (y_s == 1).sum()
            tp_loss = (flagged & (y_s == 0)).sum() / tp_total * 100 if tp_total > 0 else 0
            fp_removal = (flagged & (y_s == 1)).sum() / fp_total * 100 if fp_total > 0 else 0
            ax.scatter(tp_loss, fp_removal, marker=marker, s=size, alpha=0.6,
                      label=f"t={t_val}" if sample == SAMPLE_ORDER[0] else "")

    ax.axvline(2.0, color="red", linestyle="--", alpha=0.7, label="2% TP loss")
    ax.set_xlabel("TP Loss (%)")
    ax.set_ylabel("FP Removed (%)")
    ax.set_title("G7-Fig07: TP Lost vs FP Removed (per sample × threshold)",
                 fontsize=13, fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=8)
    ax.grid(alpha=0.3)
    save_figure(fig, FIG_DIR / "G7_fig07_tp_fp_scatter.png")


def fig08_effect_direction_forest(df: pd.DataFrame):
    """Top-3 features effect direction forest plot."""
    _print("\n[G7] Fig08: Effect direction forest")

    pass_df = df[df["filter_is_pass"] == 1].copy()
    top_features = ["af", "is_hom_alt", "h_flag"]

    fig, axes = plt.subplots(1, len(top_features), figsize=(5 * len(top_features), 6))

    for ax_idx, feat in enumerate(top_features):
        ax = axes[ax_idx]
        for i, sample in enumerate(SAMPLE_ORDER):
            sub = pass_df[pass_df["sample"] == sample]
            tp = sub[sub["truth_label"] == "TP"][feat].astype(float).dropna()
            fp = sub[sub["truth_label"] == "FP"][feat].astype(float).dropna()
            if len(tp) < 10 or len(fp) < 10:
                continue
            delta = fp.mean() - tp.mean()
            # Bootstrap CI
            deltas_boot = []
            rng = np.random.RandomState(42)
            for _ in range(500):
                tp_b = tp.sample(len(tp), replace=True)
                fp_b = fp.sample(len(fp), replace=True)
                deltas_boot.append(fp_b.mean() - tp_b.mean())
            ci_lo, ci_hi = np.percentile(deltas_boot, [2.5, 97.5])

            color = COLOR_FP if delta > 0 else COLOR_TP
            ax.plot(delta, i, "o", color=color, markersize=6)
            ax.plot([ci_lo, ci_hi], [i, i], "-", color=color, linewidth=1.5)

        ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
        ax.set_yticks(range(len(SAMPLE_ORDER)))
        ax.set_yticklabels(SAMPLE_ORDER)
        ax.set_xlabel("FP - TP mean delta")
        ax.set_title(feat, fontsize=11)

    fig.suptitle("G7-Fig08: Top-3 Features Effect Direction (per sample)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G7_fig08_effect_forest.png")


def write_go_nogo_verdict(df: pd.DataFrame):
    """Write final Go/No-Go verdict based on analysis results."""
    _print("\n[G7] === Go/No-Go Verdict ===")

    pass_df = df[df["filter_is_pass"] == 1].copy()
    score = compute_combined_score(pass_df)
    y = (pass_df["truth_label"] == "FP").astype(int).values
    mask = np.isfinite(score)

    pooled_auc = compute_auc(y[mask], score[mask])

    # Optimal threshold
    y_m, s_m = y[mask], score[mask]
    tp_total = (y_m == 0).sum()
    fp_total = (y_m == 1).sum()
    best_t, best_fp_rem, best_tp_loss = 0.5, 0, 1.0
    for t in np.linspace(0, 1, 500):
        flagged = s_m >= t
        tl = (flagged & (y_m == 0)).sum() / tp_total
        fr = (flagged & (y_m == 1)).sum() / fp_total
        if tl <= 0.02 and fr > best_fp_rem / fp_total:
            best_fp_rem = (flagged & (y_m == 1)).sum()
            best_tp_loss = tl
            best_t = t

    fp_removal_pct = best_fp_rem / fp_total * 100

    # Per-sample AUC consistency
    samples = pass_df["sample"].values
    per_sample_aucs = []
    for sample in SAMPLE_ORDER:
        m = (samples[mask] == sample)
        if m.sum() > 50:
            per_sample_aucs.append(compute_auc(y[mask][m], score[mask][m]))

    n_above_06 = sum(1 for a in per_sample_aucs if a >= 0.60)

    verdict_lines = [
        "# G7 Go/No-Go Verdict",
        "",
        f"Pooled AUC: {pooled_auc:.4f}",
        f"Optimal threshold: {best_t:.3f}",
        f"FP removal at ≤2% TP loss: {fp_removal_pct:.1f}%",
        f"TP loss at optimal threshold: {best_tp_loss*100:.2f}%",
        f"Per-sample AUC ≥ 0.60: {n_above_06}/{len(per_sample_aucs)}",
        "",
    ]

    if pooled_auc >= 0.70 and fp_removal_pct >= 10 and n_above_06 >= 5:
        verdict_lines.append("**VERDICT: GO** — develop as production filter")
    elif pooled_auc >= 0.60 and n_above_06 >= 3:
        verdict_lines.append("**VERDICT: CONDITIONAL GO** — needs refinement")
    else:
        verdict_lines.append("**VERDICT: NO-GO** — insufficient discriminative power")

    verdict_text = "\n".join(verdict_lines)
    _print(verdict_text)

    verdict_path = DATA_DIR / "go_nogo_verdict.md"
    verdict_path.write_text(verdict_text + "\n")
    _print(f"  Written verdict to {verdict_path}")


def main():
    setup_plot_style()
    _print("=" * 70)
    _print("G7: Cross-Sample Validation & Paired Ground Truth")
    _print("=" * 70)

    df = load_pass_features()
    _print(f"  Loaded {len(df):,} rows")

    fig01_consistency_heatmap(df)
    fig03_f1_improvement(df)
    fig04_bootstrap_auc(df)
    fig05_loso_auc(df)
    fig06_f1_before_after(df)
    fig07_tp_fp_scatter(df)
    fig08_effect_direction_forest(df)
    write_go_nogo_verdict(df)

    _print(f"\n[G7] DONE. Output: {G7_DIR}")


if __name__ == "__main__":
    main()
