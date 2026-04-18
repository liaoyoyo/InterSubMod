#!/usr/bin/env python3
"""PERMANOVA Applicability Audit in LOH Regions

Step 3.1: Assess whether PERMANOVA preconditions hold in LOH regions.

Core question: When HP balance disappears in LOH regions, is HP-based
PERMANOVA still valid? Is min_reads_per_group=3 commonly violated?

Data: Master dataset (748K rows), both modes, LongPhase-TO baseline.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure,
    compute_auc, compute_effect_size, compute_enrichment,
    bootstrap_ci, write_round_context, write_feature_statistics,
    ensure_dir, safe_div, encode_truth_binary,
    mannwhitney_test, chi2_test, format_p,
    SAMPLE_ORDER, MODE_ORDER, OUTPUT_ROOT,
    TRUTH_PALETTE, MODE_PALETTE,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from scipy import stats as sp_stats

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260408_permanova_loh_audit")


def classify_loh(df: pd.DataFrame) -> pd.DataFrame:
    """Add LOH classification columns."""
    df["hp_imbalance"] = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1", "yes"])
    df["loh_bed"] = df["to_loh_bed_hit"].astype(str).str.lower().isin(["true", "1", "yes"])
    # Use either definition as "LOH region"
    df["any_loh"] = df["hp_imbalance"] | df["loh_bed"]
    df["min_hp_group"] = df[["HP1FamilyN", "HP2FamilyN"]].min(axis=1)
    df["hp_valid_bool"] = df["LabelHPPermanovaValid"].astype(str).str.lower().isin(["true", "1", "yes"])
    df["hp_sig_bool"] = df["HPMergedSig"].astype(str).str.lower().isin(["true", "1", "yes"])
    df["allele_valid_bool"] = df["LabelAllelePermanovaValid"].astype(str).str.lower().isin(["true", "1", "yes"])
    df["allele_sig_bool"] = df["AlleleSig"].astype(str).str.lower().isin(["true", "1", "yes"])
    return df


def fig01_valid_rate_comparison(df: pd.DataFrame) -> None:
    """F1: HP PERMANOVA valid rate: LOH vs non-LOH, split by mode."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, mode in zip(axes, MODE_ORDER):
        dm = df[df["mode"] == mode]

        groups = ["HP Imbalance", "LOH.bed", "Any LOH", "Non-LOH"]
        masks = [dm["hp_imbalance"], dm["loh_bed"], dm["any_loh"], ~dm["any_loh"]]
        valid_rates = []
        counts = []

        for mask in masks:
            subset = dm[mask]
            rate = subset["hp_valid_bool"].mean() if len(subset) > 0 else 0
            valid_rates.append(rate * 100)
            counts.append(len(subset))

        colors = ["#D32F2F", "#1976D2", "#FF9800", "#4CAF50"]
        bars = ax.bar(range(len(groups)), valid_rates, color=colors, alpha=0.85)

        for bar, rate, n in zip(bars, valid_rates, counts):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f"{rate:.1f}%\n(n={n:,})", ha="center", va="bottom", fontsize=8)

        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, fontsize=9, rotation=15)
        ax.set_ylabel("HP PERMANOVA Valid Rate (%)")
        ax.set_title(f"{mode.upper()} mode", fontsize=12)
        ax.set_ylim(0, 110)
        ax.axhline(50, color="gray", linestyle=":", alpha=0.5)

    fig.suptitle("HP PERMANOVA Valid Rate: LOH vs Non-LOH Regions", fontsize=14, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "01_hp_permanova_valid_rate.png")


def fig02_min_hp_group_distribution(df: pd.DataFrame) -> None:
    """F2: min(HP1FamilyN, HP2FamilyN) distribution in LOH vs non-LOH."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, mode in zip(axes, MODE_ORDER):
        dm = df[df["mode"] == mode]

        for label, mask, color in [
            ("LOH (Any)", dm["any_loh"], "#D32F2F"),
            ("Non-LOH", ~dm["any_loh"], "#4CAF50"),
        ]:
            subset = dm[mask]["min_hp_group"].dropna()
            subset_clipped = subset.clip(upper=50)  # clip for visibility
            ax.hist(subset_clipped, bins=50, alpha=0.5, label=f"{label} (n={len(subset):,})",
                    color=color, density=True)

        ax.axvline(3, color="red", linestyle="--", linewidth=1.5,
                   label="min_reads_per_group=3")
        ax.set_xlabel("min(HP1FamilyN, HP2FamilyN)")
        ax.set_ylabel("Density")
        ax.set_title(f"{mode.upper()} mode", fontsize=12)
        ax.legend(fontsize=8)

        # Compute % below threshold
        for label, mask in [("LOH", dm["any_loh"]), ("Non-LOH", ~dm["any_loh"])]:
            subset = dm[mask]["min_hp_group"]
            below3 = (subset < 3).mean() * 100
            ax.text(0.98, 0.95 if label == "LOH" else 0.88,
                    f"{label} < 3: {below3:.1f}%",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    fig.suptitle("Minority HP Group Size Distribution\n"
                 "Red dashed: PERMANOVA min_reads_per_group threshold",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "02_min_hp_group_distribution.png")


def fig03_hp_merged_sig_rate(df: pd.DataFrame) -> None:
    """F3: HPMergedSig rate comparison (LOH vs non-LOH, by mode)."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, mode in zip(axes, MODE_ORDER):
        dm = df[df["mode"] == mode]

        groups = ["HP Imbalance", "LOH.bed", "Any LOH", "Non-LOH"]
        masks = [dm["hp_imbalance"], dm["loh_bed"], dm["any_loh"], ~dm["any_loh"]]

        for truth in ["TP", "FP"]:
            rates = []
            for mask in masks:
                subset = dm[mask & (dm["truth_label"] == truth)]
                rate = subset["hp_sig_bool"].mean() * 100 if len(subset) > 0 else 0
                rates.append(rate)

            x = np.arange(len(groups))
            offset = -0.2 if truth == "TP" else 0.2
            color = COLOR_TP if truth == "TP" else COLOR_FP
            bars = ax.bar(x + offset, rates, 0.35, label=truth, color=color, alpha=0.85)
            for bar, rate in zip(bars, rates):
                if rate > 0:
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                            f"{rate:.1f}", ha="center", va="bottom", fontsize=7)

        ax.set_xticks(x)
        ax.set_xticklabels(groups, fontsize=9, rotation=15)
        ax.set_ylabel("HPMergedSig Rate (%)")
        ax.set_title(f"{mode.upper()} mode", fontsize=12)
        ax.legend(fontsize=8)

    fig.suptitle("HP Merged Significance Rate: LOH vs Non-LOH\n(split by TP/FP)", fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "03_hp_merged_sig_rate.png")


def fig04_pseudo_f_distribution(df: pd.DataFrame) -> None:
    """F4: HP PERMANOVA pseudo-F distribution (LOH vs non-LOH, valid only)."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, mode in zip(axes, MODE_ORDER):
        dm = df[(df["mode"] == mode) & df["hp_valid_bool"]]

        data_loh = dm[dm["any_loh"]]["LabelHPPermanovaF"].dropna()
        data_non = dm[~dm["any_loh"]]["LabelHPPermanovaF"].dropna()

        parts = ax.violinplot(
            [data_loh.clip(upper=50), data_non.clip(upper=50)],
            positions=[0, 1], showmeans=True, showmedians=True,
        )

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["LOH\n(valid only)", "Non-LOH\n(valid only)"])
        ax.set_ylabel("HP PERMANOVA pseudo-F (clipped at 50)")
        ax.set_title(f"{mode.upper()} mode\n"
                     f"LOH median={data_loh.median():.2f} (n={len(data_loh):,}), "
                     f"Non-LOH median={data_non.median():.2f} (n={len(data_non):,})",
                     fontsize=10)

    fig.suptitle("HP PERMANOVA pseudo-F Distribution\n(Only valid PERMANOVA results)",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "04_pseudo_f_distribution.png")


def fig05_permanova_p_cdf(df: pd.DataFrame) -> None:
    """F5: HP PERMANOVA p-value CDF (LOH vs non-LOH)."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, mode in zip(axes, MODE_ORDER):
        dm = df[(df["mode"] == mode) & df["hp_valid_bool"]]

        for label, mask, color in [
            ("LOH", dm["any_loh"], "#D32F2F"),
            ("Non-LOH", ~dm["any_loh"], "#4CAF50"),
        ]:
            pvals = dm[mask]["LabelHPPermanovaP"].dropna()
            sorted_p = np.sort(pvals)
            cdf = np.arange(1, len(sorted_p) + 1) / len(sorted_p)
            ax.plot(sorted_p, cdf, color=color, label=f"{label} (n={len(pvals):,})", linewidth=1.5)

        ax.axvline(0.05, color="gray", linestyle="--", alpha=0.5, label="p=0.05")
        ax.set_xlabel("HP PERMANOVA p-value")
        ax.set_ylabel("CDF")
        ax.set_title(f"{mode.upper()} mode", fontsize=12)
        ax.legend(fontsize=8)

    fig.suptitle("HP PERMANOVA p-value CDF: LOH vs Non-LOH", fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "05_permanova_p_cdf.png")


def fig06_tp_fp_pseudo_f(df: pd.DataFrame) -> None:
    """F6: TP vs FP pseudo-F comparison within LOH regions."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, mode in zip(axes, MODE_ORDER):
        dm = df[(df["mode"] == mode) & df["any_loh"] & df["hp_valid_bool"]]

        tp_f = dm[dm["truth_label"] == "TP"]["LabelHPPermanovaF"].dropna().clip(upper=50)
        fp_f = dm[dm["truth_label"] == "FP"]["LabelHPPermanovaF"].dropna().clip(upper=50)

        bp = ax.boxplot([tp_f, fp_f], labels=["TP", "FP"],
                        patch_artist=True, showfliers=False)
        bp["boxes"][0].set_facecolor(COLOR_TP)
        bp["boxes"][1].set_facecolor(COLOR_FP)
        for b in bp["boxes"]:
            b.set_alpha(0.6)

        # AUC
        if len(tp_f) > 10 and len(fp_f) > 10:
            y_true = np.concatenate([np.ones(len(tp_f)), np.zeros(len(fp_f))])
            y_score = np.concatenate([tp_f, fp_f])
            auc_val = compute_auc(y_true, y_score)
            ax.text(0.98, 0.95, f"AUC={auc_val:.3f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow"))

        ax.set_ylabel("HP PERMANOVA pseudo-F (clipped at 50)")
        ax.set_title(f"{mode.upper()} mode (LOH regions only)\n"
                     f"TP n={len(tp_f):,}, FP n={len(fp_f):,}", fontsize=10)

    fig.suptitle("TP vs FP: HP PERMANOVA pseudo-F in LOH Regions", fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "06_tp_fp_pseudo_f_loh.png")


def fig07_per_sample_valid_heatmap(df: pd.DataFrame) -> None:
    """F7: Per-sample HP PERMANOVA valid rate heatmap."""
    setup_plot_style()

    rows = []
    for sample in SAMPLE_ORDER:
        for mode in MODE_ORDER:
            dm = df[(df["sample"] == sample) & (df["mode"] == mode)]
            for loh_label, mask in [("LOH", dm["any_loh"]), ("Non-LOH", ~dm["any_loh"])]:
                subset = dm[mask]
                rate = subset["hp_valid_bool"].mean() * 100 if len(subset) > 0 else float("nan")
                rows.append({
                    "sample": sample,
                    "mode": mode,
                    "loh_status": loh_label,
                    "valid_rate": round(rate, 1),
                    "label": f"{mode}_{loh_label}",
                })

    hm_df = pd.DataFrame(rows)
    pivot = hm_df.pivot(index="sample", columns="label", values="valid_rate")
    col_order = ["paired_LOH", "paired_Non-LOH", "to_LOH", "to_Non-LOH"]
    pivot = pivot.reindex(index=SAMPLE_ORDER, columns=col_order)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="RdYlGn", ax=ax,
                linewidths=0.5, vmin=0, vmax=100)
    ax.set_title("HP PERMANOVA Valid Rate (%) by Sample x Mode x LOH Status", fontsize=12)
    ax.set_ylabel("Sample")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "07_per_sample_valid_heatmap.png")


def fig08_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """F8: Summary statistics table (rendered as figure)."""
    setup_plot_style()

    rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for loh_label, mask in [("LOH", dm["any_loh"]), ("Non-LOH", ~dm["any_loh"])]:
            subset = dm[mask]
            n = len(subset)
            valid_rate = subset["hp_valid_bool"].mean() * 100
            sig_rate = subset["hp_sig_bool"].mean() * 100
            min_hp_below3 = (subset["min_hp_group"] < 3).mean() * 100
            median_f = subset.loc[subset["hp_valid_bool"], "LabelHPPermanovaF"].median()

            # Allele valid rate for comparison
            allele_valid = subset["allele_valid_bool"].mean() * 100
            allele_sig = subset["allele_sig_bool"].mean() * 100

            rows.append({
                "Mode": mode.upper(),
                "LOH Status": loh_label,
                "N": f"{n:,}",
                "HP Valid (%)": f"{valid_rate:.1f}",
                "HP Sig (%)": f"{sig_rate:.1f}",
                "min(HP) < 3 (%)": f"{min_hp_below3:.1f}",
                "median F": f"{median_f:.2f}" if np.isfinite(median_f) else "N/A",
                "Allele Valid (%)": f"{allele_valid:.1f}",
                "Allele Sig (%)": f"{allele_sig:.1f}",
            })

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(OUT_DIR / "summary_statistics.tsv", sep="\t", index=False)

    # Render as figure
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.axis("off")
    table = ax.table(
        cellText=summary_df.values,
        colLabels=summary_df.columns,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)

    # Color coding
    for i in range(len(summary_df)):
        for j in range(len(summary_df.columns)):
            cell = table[i + 1, j]
            if summary_df.columns[j] == "LOH Status":
                if summary_df.iloc[i]["LOH Status"] == "LOH":
                    cell.set_facecolor("#FFCDD2")
                else:
                    cell.set_facecolor("#C8E6C9")

    ax.set_title("PERMANOVA Applicability Summary", fontsize=13, pad=20)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "08_summary_table.png")
    return summary_df


# ===== Main =====

def main():
    print("=" * 70)
    print("PERMANOVA Applicability Audit in LOH Regions")
    print("Step 3.1: Validity, significance, and pseudo-F analysis")
    print("=" * 70)

    # Load data
    df = load_master_dataset()
    df = classify_loh(df)

    print(f"\nTotal rows: {len(df):,}")
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        loh_n = dm["any_loh"].sum()
        print(f"  {mode.upper()}: {len(dm):,} rows, "
              f"LOH={loh_n:,} ({loh_n/len(dm)*100:.1f}%), "
              f"Non-LOH={len(dm)-loh_n:,} ({(len(dm)-loh_n)/len(dm)*100:.1f}%)")

    # Key verification: min HP group < 3 rate in LOH
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        loh_below3 = (dm[dm["any_loh"]]["min_hp_group"] < 3).mean() * 100
        non_below3 = (dm[~dm["any_loh"]]["min_hp_group"] < 3).mean() * 100
        print(f"\n[{mode.upper()}] min(HP1N, HP2N) < 3:")
        print(f"  LOH:     {loh_below3:.1f}%")
        print(f"  Non-LOH: {non_below3:.1f}%")

    # Key verification: HP PERMANOVA valid rate
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        loh_valid = dm[dm["any_loh"]]["hp_valid_bool"].mean() * 100
        non_valid = dm[~dm["any_loh"]]["hp_valid_bool"].mean() * 100
        print(f"\n[{mode.upper()}] HP PERMANOVA valid rate:")
        print(f"  LOH:     {loh_valid:.1f}%")
        print(f"  Non-LOH: {non_valid:.1f}%")

    # AUC: LabelHPPermanovaP for TP vs FP discrimination
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for loh_label, mask in [("LOH", dm["any_loh"]), ("Non-LOH", ~dm["any_loh"])]:
            subset = dm[mask & dm["hp_valid_bool"]]
            y_true = encode_truth_binary(subset)
            y_score = -subset["LabelHPPermanovaP"]  # lower p -> more significant
            valid = y_true.notna() & np.isfinite(y_score)
            auc_val = compute_auc(y_true[valid].values, y_score[valid].values)
            print(f"  [{mode.upper()} {loh_label}] HP PERMANOVA P -> TP/FP AUC: {auc_val:.3f}")

    # Generate figures
    print("\nGenerating figures...")
    fig01_valid_rate_comparison(df)
    fig02_min_hp_group_distribution(df)
    fig03_hp_merged_sig_rate(df)
    fig04_pseudo_f_distribution(df)
    fig05_permanova_p_cdf(df)
    fig06_tp_fp_pseudo_f(df)
    fig07_per_sample_valid_heatmap(df)
    summary_df = fig08_summary_table(df)

    # Feature statistics
    write_feature_statistics(
        df,
        ["LabelHPPermanovaF", "LabelHPPermanovaP", "LabelHPPermanovaValid",
         "LabelAllelePermanovaF", "LabelAllelePermanovaP", "LabelAllelePermanovaValid",
         "HPMergedSig", "AlleleSig", "HPMergedDelta", "AlleleDelta",
         "HP1FamilyN", "HP2FamilyN", "min_hp_group"],
        OUT_DIR / "feature_statistics.tsv",
    )

    # Round context
    write_round_context(
        OUT_DIR,
        observation_id="permanova_loh_audit",
        title="PERMANOVA Applicability Audit in LOH Regions",
        description=(
            "Systematic audit of HP PERMANOVA validity in LOH regions. "
            "Checks min_reads_per_group precondition, valid rates, significance rates, "
            "and pseudo-F distributions. Compares HP vs Allele PERMANOVA applicability."
        ),
        script_path="scripts/analysis/build_permanova_loh_audit.py",
        data_sources=[{
            "name": "master_dataset",
            "path": str(load_master_dataset.__defaults__),
            "version": "2026-03-27 HP fix rebuild",
            "longphase_version": "baseline",
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "n_figures": 8,
            "n_tables": 2,
            "judgment_criteria": {
                "loh_hp_valid_rate_threshold": "< 50% = HP dimension mostly invalid in LOH",
                "min_hp_below3_threshold": "> 60% = group size precondition commonly violated",
                "hp_sig_auc_threshold": "< 0.55 = HP significance has no discriminating power in LOH",
            },
        },
    )

    print("\n" + "=" * 70)
    print(f"DONE. Output: {OUT_DIR}")
    print(f"Figures: 8, Tables: 2")
    print("=" * 70)


if __name__ == "__main__":
    main()
