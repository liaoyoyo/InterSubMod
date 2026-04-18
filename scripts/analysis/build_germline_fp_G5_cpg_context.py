#!/usr/bin/env python3
"""G5: CpG-SNP Trinucleotide Context — mutation spectrum and CpG enrichment analysis.

Uses trinucleotide context from G1 to compare mutation spectra of PASS TP vs FP.
Germline variants are enriched in CpG>TpG transitions; somatic variants follow
tumor-specific mutational signatures.
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
    compute_auc, ensure_dir,
    COLOR_TP, COLOR_FP,
)

WORKSPACE = OUTPUT_ROOT / "20260401_germline_fp_identification"
G1_DIR = WORKSPACE / "G1_vcf_features"
G5_DIR = ensure_dir(WORKSPACE / "G5_cpg_context")
FIG_DIR = ensure_dir(G5_DIR / "figures")
DATA_DIR = ensure_dir(G5_DIR / "data")

# SBS-96 canonical order (pyrimidine context)
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MUTATION_6 = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
MUTATION_6_COLORS = ["#1EBFF0", "#050708", "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5"]


def _print(msg: str):
    print(msg, flush=True)


def load_pass_features() -> pd.DataFrame:
    path = G1_DIR / "all_samples_merged.tsv.gz"
    _print(f"[G5] Loading from {path.name}")
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    df["af"] = pd.to_numeric(df["af"], errors="coerce")
    _print(f"  -> {len(df):,} rows")
    return df


def _canonicalize_trinuc(trinuc: str, ref: str, alt: str):
    """Convert to pyrimidine context (canonical SBS-96 form)."""
    if ref in ("C", "T"):
        return trinuc, ref, alt
    # Complement
    trinuc_rc = trinuc[::-1].translate(COMPLEMENT)
    ref_rc = ref.translate(COMPLEMENT)
    alt_rc = alt.translate(COMPLEMENT)
    return trinuc_rc, ref_rc, alt_rc


def build_sbs96_spectrum(df: pd.DataFrame) -> pd.DataFrame:
    """Build SBS-96 spectrum from trinucleotide context."""
    counts = {}
    for _, row in df.iterrows():
        trinuc = str(row.get("trinuc_context", ""))
        ref = str(row.get("ref", ""))
        alt = str(row.get("alt", ""))
        if len(trinuc) != 3 or "N" in trinuc or ref not in "ACGT" or alt not in "ACGT":
            continue
        ct, cr, ca = _canonicalize_trinuc(trinuc, ref, alt)
        mut_type = f"{cr}>{ca}"
        sbs96_key = f"{ct[0]}[{mut_type}]{ct[2]}"
        counts[sbs96_key] = counts.get(sbs96_key, 0) + 1

    total = sum(counts.values())
    rows = []
    for key, count in sorted(counts.items()):
        mut_type = key[2:5]  # e.g., "C>A"
        rows.append({
            "sbs96": key,
            "mutation_type": mut_type,
            "count": count,
            "fraction": count / total if total > 0 else 0,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def fig01_sbs96_spectrum(df: pd.DataFrame):
    """SBS-96 mutation spectrum for FP vs TP."""
    _print("\n[G5] Fig01: SBS-96 spectrum")

    pass_df = df[df["filter_is_pass"] == 1]

    fig, axes = plt.subplots(2, 1, figsize=(20, 8), sharex=True)

    for ax_idx, (label, title) in enumerate([("TP", "True Positives"), ("FP", "False Positives")]):
        ax = axes[ax_idx]
        sub = pass_df[pass_df["truth_label"] == label]
        spec = build_sbs96_spectrum(sub)
        if spec.empty:
            continue

        # Color by mutation type
        colors = []
        for _, row in spec.iterrows():
            mt = row["mutation_type"]
            idx = MUTATION_6.index(mt) if mt in MUTATION_6 else 0
            colors.append(MUTATION_6_COLORS[idx])

        ax.bar(range(len(spec)), spec["fraction"], color=colors, width=1.0, edgecolor="none")
        ax.set_ylabel("Fraction")
        ax.set_title(f"{title} (n={len(sub):,})", fontsize=11)
        ax.set_ylim(0, spec["fraction"].max() * 1.3)

    axes[-1].set_xlabel("96 trinucleotide contexts")
    fig.suptitle("G5-Fig01: SBS-96 Mutation Spectrum (PASS TP vs FP)",
                 fontsize=14, fontweight="bold")
    add_caption(fig, "Pyrimidine context | Colors: " + " ".join(MUTATION_6))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G5_fig01_sbs96_spectrum.png")


def fig02_cpg_fraction(df: pd.DataFrame):
    """CpG>TpG fraction per sample."""
    _print("\n[G5] Fig02: CpG>TpG fraction")

    pass_df = df[df["filter_is_pass"] == 1]

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(SAMPLE_ORDER))
    width = 0.35

    rows = []
    for offset, (label, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        fracs = []
        for sample in SAMPLE_ORDER:
            sub = pass_df[(pass_df["sample"] == sample) & (pass_df["truth_label"] == label)]
            n = len(sub)
            n_cpg_ct = (sub["is_cpg_ct"] == 1).sum()
            frac = n_cpg_ct / n * 100 if n > 0 else 0
            fracs.append(frac)
            rows.append({"sample": sample, "truth_label": label,
                         "n_cpg_ct": n_cpg_ct, "n_total": n, "cpg_ct_frac": frac})
        ax.bar(x + offset * width - width / 2, fracs, width=width,
               color=color, label=label, alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right")
    ax.set_ylabel("CpG>TpG Fraction (%)")
    ax.set_title("G5-Fig02: CpG>TpG Transition Fraction (PASS TP vs FP)",
                 fontsize=13, fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    add_caption(fig, "CpG>TpG = most common germline SNP type + somatic SBS1")
    save_figure(fig, FIG_DIR / "G5_fig02_cpg_fraction.png")

    pd.DataFrame(rows).to_csv(DATA_DIR / "cpg_ct_fraction.tsv", sep="\t", index=False)


def fig03_titv_ratio(df: pd.DataFrame):
    """Ti/Tv ratio comparison."""
    _print("\n[G5] Fig03: Ti/Tv ratio")

    pass_df = df[df["filter_is_pass"] == 1]

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(SAMPLE_ORDER))
    width = 0.35

    rows = []
    for offset, (label, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        ratios = []
        for sample in SAMPLE_ORDER:
            sub = pass_df[(pass_df["sample"] == sample) & (pass_df["truth_label"] == label)]
            n_ti = (sub["is_transition"] == 1).sum()
            n_tv = len(sub) - n_ti
            ratio = n_ti / n_tv if n_tv > 0 else 0
            ratios.append(ratio)
            rows.append({"sample": sample, "truth_label": label,
                         "n_ti": n_ti, "n_tv": n_tv, "titv_ratio": ratio})
        ax.bar(x + offset * width - width / 2, ratios, width=width,
               color=color, label=label, alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right")
    ax.set_ylabel("Ti/Tv Ratio")
    ax.set_title("G5-Fig03: Transition/Transversion Ratio (PASS TP vs FP)",
                 fontsize=13, fontweight="bold")
    ax.axhline(2.0, color="gray", linestyle="--", alpha=0.5, label="Germline typical (2.0)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    save_figure(fig, FIG_DIR / "G5_fig03_titv_ratio.png")

    pd.DataFrame(rows).to_csv(DATA_DIR / "titv_ratio.tsv", sep="\t", index=False)


def fig04_mutation_type_distribution(df: pd.DataFrame):
    """6-class mutation type distribution."""
    _print("\n[G5] Fig04: 6-class mutation type")

    pass_df = df[df["filter_is_pass"] == 1]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, (label, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        ax = axes[ax_idx]
        sub = pass_df[pass_df["truth_label"] == label]
        counts = sub["mutation_type_6class"].value_counts()
        total = counts.sum()

        for i, mt in enumerate(MUTATION_6):
            c = counts.get(mt, 0)
            ax.bar(i, c / total * 100, color=MUTATION_6_COLORS[i], edgecolor="white")

        ax.set_xticks(range(len(MUTATION_6)))
        ax.set_xticklabels(MUTATION_6)
        ax.set_ylabel("Fraction (%)")
        ax.set_title(f"{label} (n={total:,})")

    fig.suptitle("G5-Fig04: 6-Class Mutation Type Distribution",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G5_fig04_mutation_type_6class.png")


def fig05_cpg_enrichment_forest(df: pd.DataFrame):
    """CpG enrichment forest plot across 7 samples."""
    _print("\n[G5] Fig05: CpG enrichment forest plot")

    pass_df = df[df["filter_is_pass"] == 1]

    fig, ax = plt.subplots(figsize=(10, 6))

    rows = []
    for i, sample in enumerate(SAMPLE_ORDER):
        sub = pass_df[pass_df["sample"] == sample]
        tp = sub[sub["truth_label"] == "TP"]
        fp = sub[sub["truth_label"] == "FP"]

        tp_cpg = (tp["is_cpg_ct"] == 1).sum()
        tp_total = len(tp)
        fp_cpg = (fp["is_cpg_ct"] == 1).sum()
        fp_total = len(fp)

        if tp_total > 0 and fp_total > 0 and tp_cpg > 0:
            # Odds ratio
            a, b = fp_cpg, fp_total - fp_cpg
            c, d = tp_cpg, tp_total - tp_cpg
            odds_ratio = (a * d) / (b * c) if b * c > 0 else float("nan")
            # 95% CI (log scale)
            import math
            se_log_or = math.sqrt(1/max(a,1) + 1/max(b,1) + 1/max(c,1) + 1/max(d,1))
            ci_lo = math.exp(math.log(odds_ratio) - 1.96 * se_log_or)
            ci_hi = math.exp(math.log(odds_ratio) + 1.96 * se_log_or)
        else:
            odds_ratio = float("nan")
            ci_lo = ci_hi = float("nan")

        rows.append({
            "sample": sample, "odds_ratio": odds_ratio,
            "ci_lo": ci_lo, "ci_hi": ci_hi,
            "fp_cpg_frac": fp_cpg / fp_total * 100 if fp_total > 0 else 0,
            "tp_cpg_frac": tp_cpg / tp_total * 100 if tp_total > 0 else 0,
        })

        if not np.isnan(odds_ratio):
            ax.plot(odds_ratio, i, "o", color=COLOR_FP, markersize=8)
            ax.plot([ci_lo, ci_hi], [i, i], "-", color=COLOR_FP, linewidth=2)

    ax.axvline(1.0, color="gray", linestyle="--", alpha=0.5)
    ax.set_yticks(range(len(SAMPLE_ORDER)))
    ax.set_yticklabels(SAMPLE_ORDER)
    ax.set_xlabel("Odds Ratio (FP CpG>TpG enrichment vs TP)")
    ax.set_title("G5-Fig05: CpG>TpG Enrichment Forest Plot",
                 fontsize=13, fontweight="bold")
    add_caption(fig, "OR > 1 → FP enriched in CpG>TpG | 95% CI")
    save_figure(fig, FIG_DIR / "G5_fig05_cpg_enrichment_forest.png")

    pd.DataFrame(rows).to_csv(DATA_DIR / "cpg_enrichment.tsv", sep="\t", index=False)


def fig06_cpg_pr(df: pd.DataFrame):
    """CpG flag as classifier PR curve."""
    _print("\n[G5] Fig06: CpG PR curve")
    from sklearn.metrics import precision_recall_curve, roc_curve

    pass_df = df[df["filter_is_pass"] == 1].copy()
    y_true = (pass_df["truth_label"] == "FP").astype(int).values

    # Simple CpG-based score: is_cpg_ct + is_cpg_site + is_transition
    score = (pass_df["is_cpg_ct"].astype(float).fillna(0) * 0.5 +
             pass_df["is_cpg_site"].astype(float).fillna(0) * 0.3 +
             pass_df["is_transition"].astype(float).fillna(0) * 0.2).values

    mask = np.isfinite(score)
    auc = compute_auc(y_true[mask], score[mask])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    precision, recall, _ = precision_recall_curve(y_true[mask], score[mask])
    axes[0].plot(recall, precision, color=COLOR_FP, linewidth=2)
    axes[0].set_xlabel("Recall")
    axes[0].set_ylabel("Precision")
    axes[0].set_title("Precision-Recall")

    fpr, tpr, _ = roc_curve(y_true[mask], score[mask])
    axes[1].plot(fpr, tpr, color=COLOR_FP, linewidth=2, label=f"AUC={auc:.3f}")
    axes[1].plot([0, 1], [0, 1], "k--", alpha=0.3)
    axes[1].set_xlabel("FPR")
    axes[1].set_ylabel("TPR")
    axes[1].set_title("ROC")
    axes[1].legend()

    fig.suptitle("G5-Fig06: CpG Context as FP Classifier",
                 fontsize=14, fontweight="bold")
    add_caption(fig, f"Score = 0.5*is_cpg_ct + 0.3*is_cpg_site + 0.2*is_transition | AUC={auc:.3f}")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G5_fig06_cpg_pr.png")

    _print(f"  [G5] CpG context AUC = {auc:.4f}")


def fig07_sbs_heatmap(df: pd.DataFrame):
    """Per-sample SBS spectrum heatmap."""
    _print("\n[G5] Fig07: SBS heatmap")

    pass_df = df[df["filter_is_pass"] == 1]

    # Build spectrum per sample × truth_label
    rows = []
    for sample in SAMPLE_ORDER:
        for label in ["TP", "FP"]:
            sub = pass_df[(pass_df["sample"] == sample) & (pass_df["truth_label"] == label)]
            n = len(sub)
            if n == 0:
                continue
            for mt in MUTATION_6:
                count = (sub["mutation_type_6class"] == mt).sum()
                rows.append({
                    "sample": sample, "truth_label": label,
                    "mutation_type": mt, "fraction": count / n * 100,
                })

    if not rows:
        return

    df_spec = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax_idx, label in enumerate(["TP", "FP"]):
        ax = axes[ax_idx]
        sub = df_spec[df_spec["truth_label"] == label]
        pivot = sub.pivot_table(index="sample", columns="mutation_type", values="fraction", fill_value=0)
        pivot = pivot.reindex(index=SAMPLE_ORDER, columns=MUTATION_6).fillna(0)
        sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd", ax=ax, vmin=0)
        ax.set_title(f"{label}")
        ax.set_ylabel("")

    fig.suptitle("G5-Fig07: Per-Sample Mutation Type Heatmap (%)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G5_fig07_sbs_heatmap.png")


def fig08_cpg_titv_combined_roc(df: pd.DataFrame):
    """Combined CpG + Ti/Tv ROC."""
    _print("\n[G5] Fig08: Combined CpG + Ti/Tv ROC")
    from sklearn.metrics import roc_curve

    pass_df = df[df["filter_is_pass"] == 1].copy()
    y_true = (pass_df["truth_label"] == "FP").astype(int).values

    scores = {
        "is_cpg_ct": pass_df["is_cpg_ct"].astype(float).fillna(0).values,
        "is_cpg_site": pass_df["is_cpg_site"].astype(float).fillna(0).values,
        "is_transition": pass_df["is_transition"].astype(float).fillna(0).values,
        "Combined (CpG+Ti)": (pass_df["is_cpg_ct"].astype(float).fillna(0) * 0.5 +
                               pass_df["is_cpg_site"].astype(float).fillna(0) * 0.3 +
                               pass_df["is_transition"].astype(float).fillna(0) * 0.2).values,
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    for name, score in scores.items():
        mask = np.isfinite(score)
        auc = compute_auc(y_true[mask], score[mask])
        fpr, tpr, _ = roc_curve(y_true[mask], score[mask])
        lw = 2 if "Combined" in name else 1.5
        ax.plot(fpr, tpr, label=f"{name} (AUC={auc:.3f})", linewidth=lw)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.set_title("G5-Fig08: CpG + Ti/Tv Combined ROC", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    save_figure(fig, FIG_DIR / "G5_fig08_cpg_titv_roc.png")


def main():
    setup_plot_style()
    _print("=" * 70)
    _print("G5: CpG-SNP Trinucleotide Context")
    _print("=" * 70)

    df = load_pass_features()

    fig01_sbs96_spectrum(df)
    fig02_cpg_fraction(df)
    fig03_titv_ratio(df)
    fig04_mutation_type_distribution(df)
    fig05_cpg_enrichment_forest(df)
    fig06_cpg_pr(df)
    fig07_sbs_heatmap(df)
    fig08_cpg_titv_combined_roc(df)

    _print(f"\n[G5] DONE. Output: {G5_DIR}")


if __name__ == "__main__":
    main()
