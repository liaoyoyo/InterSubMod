#!/usr/bin/env python3
"""
S0 Self-Phasing Circular Dependency Verification + VCF Layer Analysis
=====================================================================
Verifies whether LongPhase-TO incorporates somatic variants into phasing scaffold
(self-phasing circular dependency) and quantifies its effects.

Sections:
  S0.1 - Somatic variants present in TO phased VCF
  S0.2 - Self-phasing effect quantification (HP imbalance)
  S0.3 - Paired vs TO imbalance comparison (control)
  S1.1 - VCF layer variant census
  S1.4 - Somatic-as-anchor direct count
"""

import gzip
import os
import sys
import warnings
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ──────────────────────────────────────────────────────────────────
OUT_DIR = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260402_phasing_causal_chain"

TO_PHASED_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased.vcf"
TRUTH_TP_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step02_benchmark_clairs_to/filtered_snv_tp.vcf.gz"
TRUTH_FP_VCF = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step02_benchmark_clairs_to/filtered_snv_fp.vcf.gz"
MASTER_DS = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]

os.makedirs(OUT_DIR, exist_ok=True)


# ── Helper: read VCF variant keys ──────────────────────────────────────────
def read_vcf_keys(path):
    """Read a VCF (plain or .gz) and return set of chr:pos:ref:alt keys."""
    keys = set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 6)
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            key = f"{chrom}:{pos}:{ref}:{alt}"
            keys.add(key)
    return keys


def cohens_d(a, b):
    """Cohen's d for two groups."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    ma, mb = np.mean(a), np.mean(b)
    sa, sb = np.std(a, ddof=1), np.std(b, ddof=1)
    pooled = np.sqrt(((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2))
    if pooled == 0:
        return 0.0
    return (ma - mb) / pooled


# ═══════════════════════════════════════════════════════════════════════════
# S0.1: Somatic variants in TO phased VCF
# ═══════════════════════════════════════════════════════════════════════════
def run_s0_1():
    print("=" * 70)
    print("S0.1: Somatic variants in TO phased VCF")
    print("=" * 70)

    # Load truth keys
    print("  Loading truth TP keys...", flush=True)
    tp_keys = read_vcf_keys(TRUTH_TP_VCF)
    print(f"  Truth TP variants: {len(tp_keys)}")

    print("  Loading truth FP keys...", flush=True)
    fp_keys = read_vcf_keys(TRUTH_FP_VCF)
    print(f"  Truth FP variants: {len(fp_keys)}")

    # Scan TO phased VCF line by line
    print("  Scanning TO phased VCF (655MB, streaming)...", flush=True)

    results = []
    total_lines = 0
    matched_tp = 0
    matched_fp = 0
    tp_phased = 0  # has PS and phased GT (|)
    fp_phased = 0
    tp_het_phased = 0  # GT is 0|1 or 1|0
    fp_het_phased = 0

    with open(TO_PHASED_VCF, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            total_lines += 1
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            filt = parts[6]
            fmt_keys = parts[8].split(":")
            fmt_vals = parts[9].split(":")
            fmt_dict = dict(zip(fmt_keys, fmt_vals))

            key = f"{chrom}:{pos}:{ref}:{alt}"
            gt = fmt_dict.get("GT", "./.")
            ps = fmt_dict.get("PS", ".")
            has_ps = ps != "." and ps != ""
            is_phased_gt = "|" in gt
            is_het_phased = gt in ("0|1", "1|0")

            if key in tp_keys:
                matched_tp += 1
                if is_phased_gt and has_ps:
                    tp_phased += 1
                if is_het_phased:
                    tp_het_phased += 1
                results.append({
                    "variant_key": key,
                    "truth_label": "TP",
                    "filter": filt,
                    "GT": gt,
                    "PS": ps,
                    "is_phased": is_phased_gt and has_ps,
                    "is_het_phased": is_het_phased,
                    "AF": fmt_dict.get("AF", "."),
                    "DP": fmt_dict.get("DP", "."),
                    "GQ": fmt_dict.get("GQ", "."),
                })
            elif key in fp_keys:
                matched_fp += 1
                if is_phased_gt and has_ps:
                    fp_phased += 1
                if is_het_phased:
                    fp_het_phased += 1
                results.append({
                    "variant_key": key,
                    "truth_label": "FP",
                    "filter": filt,
                    "GT": gt,
                    "PS": ps,
                    "is_phased": is_phased_gt and has_ps,
                    "is_het_phased": is_het_phased,
                    "AF": fmt_dict.get("AF", "."),
                    "DP": fmt_dict.get("DP", "."),
                    "GQ": fmt_dict.get("GQ", "."),
                })

    # Save results
    df = pd.DataFrame(results)
    out_path = os.path.join(OUT_DIR, "s0_1_somatic_in_phased_vcf.tsv")
    df.to_csv(out_path, sep="\t", index=False)

    print(f"\n  --- S0.1 Summary ---")
    print(f"  TO phased VCF total data lines: {total_lines:,}")
    print(f"  Truth TP found in phased VCF: {matched_tp} / {len(tp_keys)} ({100*matched_tp/max(len(tp_keys),1):.1f}%)")
    print(f"  Truth FP found in phased VCF: {matched_fp} / {len(fp_keys)} ({100*matched_fp/max(len(fp_keys),1):.1f}%)")
    print(f"  TP phased (GT has | + PS assigned): {tp_phased} / {matched_tp} ({100*tp_phased/max(matched_tp,1):.1f}%)")
    print(f"  FP phased (GT has | + PS assigned): {fp_phased} / {matched_fp} ({100*fp_phased/max(matched_fp,1):.1f}%)")
    print(f"  TP het-phased (0|1 or 1|0): {tp_het_phased} / {matched_tp} ({100*tp_het_phased/max(matched_tp,1):.1f}%)")
    print(f"  FP het-phased (0|1 or 1|0): {fp_het_phased} / {matched_fp} ({100*fp_het_phased/max(matched_fp,1):.1f}%)")
    print(f"  Output: {out_path}")

    return {
        "total_phased_vcf_lines": total_lines,
        "tp_total": len(tp_keys),
        "fp_total": len(fp_keys),
        "tp_in_phased": matched_tp,
        "fp_in_phased": matched_fp,
        "tp_phased": tp_phased,
        "fp_phased": fp_phased,
        "tp_het_phased": tp_het_phased,
        "fp_het_phased": fp_het_phased,
    }


# ═══════════════════════════════════════════════════════════════════════════
# S0.2: Self-phasing effect quantification
# ═══════════════════════════════════════════════════════════════════════════
def run_s0_2(master):
    print("\n" + "=" * 70)
    print("S0.2: Self-Phasing Effect Quantification (HP Imbalance)")
    print("=" * 70)

    # Compute HP imbalance
    hp1 = master["HP1FamilyN"].astype(float)
    hp2 = master["HP2FamilyN"].astype(float)
    master = master.copy()
    master["hp_imbalance"] = np.abs(hp1 - hp2) / (hp1 + hp2 + 0.001)

    # --- Per-sample analysis ---
    rows = []
    for samp in SAMPLE_ORDER:
        for mode in ["to", "paired"]:
            for truth in ["TP", "FP"]:
                mask = (master["sample"] == samp) & (master["mode"] == mode) & (master["truth_label"] == truth)
                subset = master.loc[mask, "hp_imbalance"]
                n = len(subset)
                if n == 0:
                    rows.append({
                        "sample": samp, "mode": mode, "truth_label": truth,
                        "n": 0, "mean": np.nan, "median": np.nan, "std": np.nan
                    })
                    continue
                rows.append({
                    "sample": samp, "mode": mode, "truth_label": truth,
                    "n": n,
                    "mean": subset.mean(),
                    "median": subset.median(),
                    "std": subset.std(),
                })

    df_sample = pd.DataFrame(rows)

    # Effect sizes (Cohen's d)
    effect_rows = []
    for samp in SAMPLE_ORDER:
        to_tp = master.loc[(master["sample"] == samp) & (master["mode"] == "to") & (master["truth_label"] == "TP"), "hp_imbalance"]
        to_fp = master.loc[(master["sample"] == samp) & (master["mode"] == "to") & (master["truth_label"] == "FP"), "hp_imbalance"]
        pa_tp = master.loc[(master["sample"] == samp) & (master["mode"] == "paired") & (master["truth_label"] == "TP"), "hp_imbalance"]
        pa_fp = master.loc[(master["sample"] == samp) & (master["mode"] == "paired") & (master["truth_label"] == "FP"), "hp_imbalance"]

        d_to_tp_vs_fp = cohens_d(to_tp.values, to_fp.values)
        d_to_tp_vs_pa_tp = cohens_d(to_tp.values, pa_tp.values)
        d_pa_tp_vs_fp = cohens_d(pa_tp.values, pa_fp.values)

        # Mann-Whitney U test for TO TP vs FP
        if len(to_tp) > 5 and len(to_fp) > 5:
            u_stat, u_p = stats.mannwhitneyu(to_tp, to_fp, alternative="two-sided")
        else:
            u_stat, u_p = np.nan, np.nan

        effect_rows.append({
            "sample": samp,
            "TO_TP_n": len(to_tp), "TO_FP_n": len(to_fp),
            "Paired_TP_n": len(pa_tp), "Paired_FP_n": len(pa_fp),
            "TO_TP_mean": to_tp.mean() if len(to_tp) > 0 else np.nan,
            "TO_FP_mean": to_fp.mean() if len(to_fp) > 0 else np.nan,
            "Paired_TP_mean": pa_tp.mean() if len(pa_tp) > 0 else np.nan,
            "Paired_FP_mean": pa_fp.mean() if len(pa_fp) > 0 else np.nan,
            "d_TO_TP_vs_FP": d_to_tp_vs_fp,
            "d_TO_TP_vs_Paired_TP": d_to_tp_vs_pa_tp,
            "d_Paired_TP_vs_FP": d_pa_tp_vs_fp,
            "U_TO_TP_vs_FP": u_stat,
            "p_TO_TP_vs_FP": u_p,
        })

    df_effect = pd.DataFrame(effect_rows)

    out1 = os.path.join(OUT_DIR, "s0_2_self_phasing_effect_by_sample.tsv")
    df_sample.to_csv(out1, sep="\t", index=False, float_format="%.4f")
    print(f"  Saved: {out1}")

    # --- AF-stratified analysis ---
    af_bins = [
        ("<0.1", 0.0, 0.1),
        ("0.1-0.3", 0.1, 0.3),
        ("0.3-0.5", 0.3, 0.5),
        (">0.5", 0.5, 1.01),
    ]

    af_rows = []
    for label, lo, hi in af_bins:
        for mode in ["to", "paired"]:
            for truth in ["TP", "FP"]:
                mask = (
                    (master["mode"] == mode)
                    & (master["truth_label"] == truth)
                    & (master["caller_af"] >= lo)
                    & (master["caller_af"] < hi)
                )
                subset = master.loc[mask, "hp_imbalance"]
                n = len(subset)
                af_rows.append({
                    "af_bin": label, "mode": mode, "truth_label": truth,
                    "n": n,
                    "mean": subset.mean() if n > 0 else np.nan,
                    "median": subset.median() if n > 0 else np.nan,
                    "std": subset.std() if n > 0 else np.nan,
                })

    # Effect sizes per AF bin
    for label, lo, hi in af_bins:
        to_tp = master.loc[
            (master["mode"] == "to") & (master["truth_label"] == "TP")
            & (master["caller_af"] >= lo) & (master["caller_af"] < hi),
            "hp_imbalance"
        ]
        to_fp = master.loc[
            (master["mode"] == "to") & (master["truth_label"] == "FP")
            & (master["caller_af"] >= lo) & (master["caller_af"] < hi),
            "hp_imbalance"
        ]
        pa_tp = master.loc[
            (master["mode"] == "paired") & (master["truth_label"] == "TP")
            & (master["caller_af"] >= lo) & (master["caller_af"] < hi),
            "hp_imbalance"
        ]
        d1 = cohens_d(to_tp.values, to_fp.values)
        d2 = cohens_d(to_tp.values, pa_tp.values)
        af_rows.append({
            "af_bin": label, "mode": "EFFECT_SIZE", "truth_label": "d_TO_TP_vs_FP",
            "n": len(to_tp), "mean": d1, "median": np.nan, "std": np.nan,
        })
        af_rows.append({
            "af_bin": label, "mode": "EFFECT_SIZE", "truth_label": "d_TO_TP_vs_Paired_TP",
            "n": len(to_tp), "mean": d2, "median": np.nan, "std": np.nan,
        })

    df_af = pd.DataFrame(af_rows)
    out2 = os.path.join(OUT_DIR, "s0_2_self_phasing_effect_by_af.tsv")
    df_af.to_csv(out2, sep="\t", index=False, float_format="%.4f")
    print(f"  Saved: {out2}")

    # --- Boxplot: HP imbalance by mode x truth, faceted by sample ---
    plot_data = master[master["sample"].isin(SAMPLE_ORDER)].copy()
    plot_data["group"] = plot_data["mode"].str.upper() + "_" + plot_data["truth_label"]

    fig, axes = plt.subplots(2, 4, figsize=(24, 10), sharey=True)
    axes_flat = axes.flatten()

    for i, samp in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        sdata = plot_data[plot_data["sample"] == samp]
        order = ["TO_TP", "TO_FP", "PAIRED_TP", "PAIRED_FP"]
        palette = {"TO_TP": "#e74c3c", "TO_FP": "#3498db", "PAIRED_TP": "#e67e22", "PAIRED_FP": "#2ecc71"}
        sns.boxplot(
            data=sdata, x="group", y="hp_imbalance", order=order,
            palette=palette, ax=ax, showfliers=False, width=0.6,
        )
        # Add n labels
        for j, grp in enumerate(order):
            n = len(sdata[sdata["group"] == grp])
            ax.text(j, ax.get_ylim()[1] * 0.95, f"n={n}", ha="center", va="top", fontsize=7)
        ax.set_title(samp, fontsize=11, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("HP Imbalance" if i % 4 == 0 else "")
        ax.tick_params(axis="x", rotation=30, labelsize=8)

    # Hide last subplot if not used
    if len(SAMPLE_ORDER) < 8:
        axes_flat[-1].set_visible(False)

    fig.suptitle("S0.2: HP Imbalance — TO vs Paired × TP vs FP\n(Self-Phasing Circular Dependency Test)", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig_path = os.path.join(OUT_DIR, "s0_2_hp_imbalance_to_vs_paired.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # Print summary
    print("\n  --- S0.2 Per-Sample Effect Summary ---")
    print(f"  {'Sample':<18} {'d(TO_TP-FP)':>12} {'d(TO_TP-Pa_TP)':>15} {'d(Pa_TP-FP)':>12} {'p(TO_TP-FP)':>12}")
    for _, r in df_effect.iterrows():
        print(f"  {r['sample']:<18} {r['d_TO_TP_vs_FP']:>12.3f} {r['d_TO_TP_vs_Paired_TP']:>15.3f} {r['d_Paired_TP_vs_FP']:>12.3f} {r['p_TO_TP_vs_FP']:>12.2e}")

    return df_effect, df_sample, master


# ═══════════════════════════════════════════════════════════════════════════
# S0.3: Paired vs TO imbalance comparison (control)
# ═══════════════════════════════════════════════════════════════════════════
def run_s0_3(master):
    print("\n" + "=" * 70)
    print("S0.3: Paired vs TO Imbalance Comparison (Control)")
    print("=" * 70)

    rows = []
    for mode in ["to", "paired"]:
        for truth in ["TP", "FP"]:
            subset = master.loc[(master["mode"] == mode) & (master["truth_label"] == truth), "hp_imbalance"]
            rows.append({
                "mode": mode, "truth_label": truth,
                "n": len(subset),
                "mean": subset.mean(),
                "median": subset.median(),
                "std": subset.std(),
            })

    # Effect sizes: global
    to_tp = master.loc[(master["mode"] == "to") & (master["truth_label"] == "TP"), "hp_imbalance"]
    to_fp = master.loc[(master["mode"] == "to") & (master["truth_label"] == "FP"), "hp_imbalance"]
    pa_tp = master.loc[(master["mode"] == "paired") & (master["truth_label"] == "TP"), "hp_imbalance"]
    pa_fp = master.loc[(master["mode"] == "paired") & (master["truth_label"] == "FP"), "hp_imbalance"]

    d_to = cohens_d(to_tp.values, to_fp.values)
    d_pa = cohens_d(pa_tp.values, pa_fp.values)
    d_to_pa_tp = cohens_d(to_tp.values, pa_tp.values)
    d_to_pa_fp = cohens_d(to_fp.values, pa_fp.values)

    # Statistical tests
    _, p_to = stats.mannwhitneyu(to_tp, to_fp, alternative="two-sided")
    _, p_pa = stats.mannwhitneyu(pa_tp, pa_fp, alternative="two-sided")
    _, p_tp_modes = stats.mannwhitneyu(to_tp, pa_tp, alternative="two-sided")

    rows.append({"mode": "EFFECT", "truth_label": "d_TO_TP_vs_FP", "n": len(to_tp) + len(to_fp), "mean": d_to, "median": p_to, "std": np.nan})
    rows.append({"mode": "EFFECT", "truth_label": "d_Paired_TP_vs_FP", "n": len(pa_tp) + len(pa_fp), "mean": d_pa, "median": p_pa, "std": np.nan})
    rows.append({"mode": "EFFECT", "truth_label": "d_TO_TP_vs_Paired_TP", "n": len(to_tp) + len(pa_tp), "mean": d_to_pa_tp, "median": p_tp_modes, "std": np.nan})
    rows.append({"mode": "EFFECT", "truth_label": "d_TO_FP_vs_Paired_FP", "n": len(to_fp) + len(pa_fp), "mean": d_to_pa_fp, "median": np.nan, "std": np.nan})

    df = pd.DataFrame(rows)
    out_path = os.path.join(OUT_DIR, "s0_3_paired_vs_to_imbalance_comparison.tsv")
    df.to_csv(out_path, sep="\t", index=False, float_format="%.4f")
    print(f"  Saved: {out_path}")

    print(f"\n  --- S0.3 Global Imbalance Summary ---")
    print(f"  TO_TP:     mean={to_tp.mean():.4f}  median={to_tp.median():.4f}  n={len(to_tp)}")
    print(f"  TO_FP:     mean={to_fp.mean():.4f}  median={to_fp.median():.4f}  n={len(to_fp)}")
    print(f"  Paired_TP: mean={pa_tp.mean():.4f}  median={pa_tp.median():.4f}  n={len(pa_tp)}")
    print(f"  Paired_FP: mean={pa_fp.mean():.4f}  median={pa_fp.median():.4f}  n={len(pa_fp)}")
    print(f"\n  Cohen's d (TO TP vs FP):     {d_to:.4f}  (p={p_to:.2e})")
    print(f"  Cohen's d (Paired TP vs FP): {d_pa:.4f}  (p={p_pa:.2e})")
    print(f"  Cohen's d (TO_TP vs Paired_TP): {d_to_pa_tp:.4f}  (p={p_tp_modes:.2e})")

    hypothesis_supported = (abs(d_to) > abs(d_pa)) and (abs(d_to) > 0.2)
    if hypothesis_supported:
        print(f"\n  >>> SELF-PHASING HYPOTHESIS SUPPORTED: |d_TO|={abs(d_to):.3f} > |d_Paired|={abs(d_pa):.3f}")
    else:
        print(f"\n  >>> SELF-PHASING HYPOTHESIS NOT CLEARLY SUPPORTED: |d_TO|={abs(d_to):.3f}, |d_Paired|={abs(d_pa):.3f}")

    return {
        "d_to_tp_vs_fp": d_to,
        "d_paired_tp_vs_fp": d_pa,
        "d_to_tp_vs_paired_tp": d_to_pa_tp,
        "p_to": p_to,
        "p_paired": p_pa,
        "hypothesis_supported": hypothesis_supported,
    }


# ═══════════════════════════════════════════════════════════════════════════
# S1.1: VCF layer variant census
# ═══════════════════════════════════════════════════════════════════════════
def run_s1_1():
    print("\n" + "=" * 70)
    print("S1.1: TO Phased VCF Variant Census")
    print("=" * 70)

    filter_counts = defaultdict(int)
    gt_counts = defaultdict(int)
    phased_count = 0
    unphased_count = 0
    total = 0
    het_phased = 0
    hom_phased = 0
    filter_x_phased = defaultdict(lambda: {"phased": 0, "unphased": 0})

    print("  Scanning VCF (streaming)...", flush=True)
    with open(TO_PHASED_VCF, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            parts = line.rstrip("\n").split("\t")
            filt = parts[6]
            fmt_keys = parts[8].split(":")
            fmt_vals = parts[9].split(":")
            fmt_dict = dict(zip(fmt_keys, fmt_vals))

            gt = fmt_dict.get("GT", "./.")
            ps = fmt_dict.get("PS", ".")
            has_ps = ps != "." and ps != ""
            is_phased = "|" in gt

            # Count filters (may be multi-value like "LowQual;NonSomatic")
            filter_counts[filt] += 1
            gt_counts[gt] += 1

            if is_phased and has_ps:
                phased_count += 1
                filter_x_phased[filt]["phased"] += 1
                if gt in ("0|1", "1|0"):
                    het_phased += 1
                else:
                    hom_phased += 1
            else:
                unphased_count += 1
                filter_x_phased[filt]["unphased"] += 1

    # Build output table
    rows = []
    for filt, cnt in sorted(filter_counts.items(), key=lambda x: -x[1]):
        rows.append({
            "FILTER": filt,
            "count": cnt,
            "pct": 100.0 * cnt / total,
            "phased": filter_x_phased[filt]["phased"],
            "unphased": filter_x_phased[filt]["unphased"],
            "phased_rate": 100.0 * filter_x_phased[filt]["phased"] / max(cnt, 1),
        })

    df = pd.DataFrame(rows)
    out_path = os.path.join(OUT_DIR, "s1_1_to_phased_vcf_variant_census.tsv")
    df.to_csv(out_path, sep="\t", index=False, float_format="%.2f")
    print(f"  Saved: {out_path}")

    print(f"\n  --- S1.1 Summary ---")
    print(f"  Total variants: {total:,}")
    print(f"  Phased (GT has | + PS): {phased_count:,} ({100*phased_count/total:.1f}%)")
    print(f"    Het phased (0|1/1|0): {het_phased:,}")
    print(f"    Hom phased (1|./etc): {hom_phased:,}")
    print(f"  Unphased: {unphased_count:,} ({100*unphased_count/total:.1f}%)")
    print(f"\n  Top 10 FILTER categories:")
    for _, r in df.head(10).iterrows():
        print(f"    {r['FILTER']:<40} {r['count']:>8,}  ({r['pct']:5.1f}%)  phased: {r['phased_rate']:5.1f}%")

    print(f"\n  Top 10 GT values:")
    for gt, cnt in sorted(gt_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"    {gt:<10} {cnt:>10,}  ({100*cnt/total:5.1f}%)")

    return {"total": total, "phased": phased_count, "het_phased": het_phased}


# ═══════════════════════════════════════════════════════════════════════════
# S1.4: Somatic-as-Anchor Direct Count
# ═══════════════════════════════════════════════════════════════════════════
def run_s1_4():
    print("\n" + "=" * 70)
    print("S1.4: Somatic-as-Anchor Direct Count")
    print("=" * 70)

    # Load truth keys
    tp_keys = read_vcf_keys(TRUTH_TP_VCF)
    fp_keys = read_vcf_keys(TRUTH_FP_VCF)
    somatic_keys = tp_keys | fp_keys

    print(f"  Total somatic variants (TP+FP): {len(somatic_keys)}")

    # Scan phased VCF: for each somatic variant, check if it has PS
    results = []
    tp_with_ps = 0
    fp_with_ps = 0
    tp_het_with_ps = 0
    fp_het_with_ps = 0
    tp_total_found = 0
    fp_total_found = 0

    # Also collect PS block sizes (how many variants share the same PS)
    ps_blocks = defaultdict(lambda: {"total": 0, "somatic_tp": 0, "somatic_fp": 0})

    print("  Scanning VCF for somatic anchor analysis...", flush=True)
    with open(TO_PHASED_VCF, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            fmt_keys = parts[8].split(":")
            fmt_vals = parts[9].split(":")
            fmt_dict = dict(zip(fmt_keys, fmt_vals))

            key = f"{chrom}:{pos}:{ref}:{alt}"
            gt = fmt_dict.get("GT", "./.")
            ps = fmt_dict.get("PS", ".")
            has_ps = ps != "." and ps != ""
            is_phased = "|" in gt
            is_het = gt in ("0|1", "1|0", "0/1", "1/0")
            is_het_phased = gt in ("0|1", "1|0")

            # Track PS blocks
            if has_ps and is_phased:
                ps_key = f"{chrom}:{ps}"
                ps_blocks[ps_key]["total"] += 1
                if key in tp_keys:
                    ps_blocks[ps_key]["somatic_tp"] += 1
                elif key in fp_keys:
                    ps_blocks[ps_key]["somatic_fp"] += 1

            if key in tp_keys:
                tp_total_found += 1
                if has_ps and is_phased:
                    tp_with_ps += 1
                if is_het_phased:
                    tp_het_with_ps += 1
            elif key in fp_keys:
                fp_total_found += 1
                if has_ps and is_phased:
                    fp_with_ps += 1
                if is_het_phased:
                    fp_het_with_ps += 1

    # Analyze PS blocks contamination
    blocks_with_somatic = 0
    blocks_total = len(ps_blocks)
    somatic_in_blocks = 0
    total_in_blocks = 0
    blocks_somatic_only = 0  # blocks where ALL variants are somatic

    for ps_key, info in ps_blocks.items():
        total_in_blocks += info["total"]
        n_somatic = info["somatic_tp"] + info["somatic_fp"]
        if n_somatic > 0:
            blocks_with_somatic += 1
            somatic_in_blocks += n_somatic
        if n_somatic == info["total"]:
            blocks_somatic_only += 1

    # Output
    summary_rows = [
        {"metric": "TP_found_in_phased_VCF", "value": tp_total_found},
        {"metric": "TP_with_PS_assignment", "value": tp_with_ps},
        {"metric": "TP_het_phased_0|1_or_1|0", "value": tp_het_with_ps},
        {"metric": "TP_phased_rate", "value": f"{100*tp_with_ps/max(tp_total_found,1):.1f}%"},
        {"metric": "FP_found_in_phased_VCF", "value": fp_total_found},
        {"metric": "FP_with_PS_assignment", "value": fp_with_ps},
        {"metric": "FP_het_phased_0|1_or_1|0", "value": fp_het_with_ps},
        {"metric": "FP_phased_rate", "value": f"{100*fp_with_ps/max(fp_total_found,1):.1f}%"},
        {"metric": "total_PS_blocks", "value": blocks_total},
        {"metric": "PS_blocks_with_somatic_variant", "value": blocks_with_somatic},
        {"metric": "PS_blocks_contamination_rate", "value": f"{100*blocks_with_somatic/max(blocks_total,1):.1f}%"},
        {"metric": "somatic_variants_in_blocks", "value": somatic_in_blocks},
        {"metric": "total_variants_in_blocks", "value": total_in_blocks},
        {"metric": "somatic_fraction_of_block_vars", "value": f"{100*somatic_in_blocks/max(total_in_blocks,1):.2f}%"},
    ]

    df = pd.DataFrame(summary_rows)
    out_path = os.path.join(OUT_DIR, "s1_4_somatic_as_anchor_count.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")

    print(f"\n  --- S1.4 Summary ---")
    print(f"  TP somatic phased as anchor: {tp_with_ps} / {tp_total_found} ({100*tp_with_ps/max(tp_total_found,1):.1f}%)")
    print(f"  FP somatic phased as anchor: {fp_with_ps} / {fp_total_found} ({100*fp_with_ps/max(fp_total_found,1):.1f}%)")
    print(f"  Total PS blocks: {blocks_total:,}")
    print(f"  PS blocks containing somatic: {blocks_with_somatic:,} ({100*blocks_with_somatic/max(blocks_total,1):.1f}%)")
    print(f"  Somatic variants in blocks: {somatic_in_blocks:,} / {total_in_blocks:,} ({100*somatic_in_blocks/max(total_in_blocks,1):.2f}%)")

    return {
        "tp_phased": tp_with_ps,
        "fp_phased": fp_with_ps,
        "blocks_total": blocks_total,
        "blocks_with_somatic": blocks_with_somatic,
        "somatic_in_blocks": somatic_in_blocks,
        "total_in_blocks": total_in_blocks,
    }


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("Self-Phasing Circular Dependency + VCF Layer Analysis")
    print(f"Output: {OUT_DIR}")
    print("=" * 70)

    # ── S0.1 ──
    s0_1_results = run_s0_1()

    # ── Load master dataset (for S0.2, S0.3) ──
    print("\n  Loading master dataset...", flush=True)
    master = pd.read_csv(MASTER_DS, sep="\t", compression="gzip")
    print(f"  Master dataset: {len(master):,} rows, {master.shape[1]} columns")
    print(f"  Modes: {master['mode'].value_counts().to_dict()}")
    print(f"  Truth: {master['truth_label'].value_counts().to_dict()}")

    # Compute hp_imbalance on master
    hp1 = master["HP1FamilyN"].astype(float)
    hp2 = master["HP2FamilyN"].astype(float)
    master["hp_imbalance"] = np.abs(hp1 - hp2) / (hp1 + hp2 + 0.001)

    # ── S0.2 ──
    s0_2_effect, _, master = run_s0_2(master)

    # ── S0.3 ──
    s0_3_results = run_s0_3(master)

    # ── S1.1 ──
    s1_1_results = run_s1_1()

    # ── S1.4 ──
    s1_4_results = run_s1_4()

    # ═══════════════════════════════════════════════════════════════════════
    # FINAL SUMMARY
    # ═══════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("SUMMARY: Self-Phasing Circular Dependency + VCF Analysis")
    print("=" * 70)

    print(f"\n[S0.1] Somatic variants in TO phased VCF:")
    print(f"  TP found: {s0_1_results['tp_in_phased']}/{s0_1_results['tp_total']} ({100*s0_1_results['tp_in_phased']/max(s0_1_results['tp_total'],1):.1f}%)")
    print(f"  FP found: {s0_1_results['fp_in_phased']}/{s0_1_results['fp_total']} ({100*s0_1_results['fp_in_phased']/max(s0_1_results['fp_total'],1):.1f}%)")
    print(f"  TP phased: {s0_1_results['tp_phased']} ({100*s0_1_results['tp_phased']/max(s0_1_results['tp_in_phased'],1):.1f}%)")
    print(f"  FP phased: {s0_1_results['fp_phased']} ({100*s0_1_results['fp_phased']/max(s0_1_results['fp_in_phased'],1):.1f}%)")
    print(f"  TP het-phased (0|1/1|0): {s0_1_results['tp_het_phased']}")
    print(f"  FP het-phased (0|1/1|0): {s0_1_results['fp_het_phased']}")

    print(f"\n[S0.2] Self-Phasing HP Imbalance Effect (pooled across samples):")
    for _, r in s0_2_effect.iterrows():
        print(f"  {r['sample']}: d(TO_TP-FP)={r['d_TO_TP_vs_FP']:.3f}, d(TO_TP-Pa_TP)={r['d_TO_TP_vs_Paired_TP']:.3f}, d(Pa_TP-FP)={r['d_Paired_TP_vs_FP']:.3f}")

    print(f"\n[S0.3] Self-Phasing Hypothesis Verdict:")
    print(f"  d(TO TP vs FP):     {s0_3_results['d_to_tp_vs_fp']:.4f}  (p={s0_3_results['p_to']:.2e})")
    print(f"  d(Paired TP vs FP): {s0_3_results['d_paired_tp_vs_fp']:.4f}  (p={s0_3_results['p_paired']:.2e})")
    print(f"  d(TO_TP vs Pa_TP):  {s0_3_results['d_to_tp_vs_paired_tp']:.4f}")
    print(f"  Hypothesis supported: {s0_3_results['hypothesis_supported']}")

    print(f"\n[S1.1] TO Phased VCF Census:")
    print(f"  Total variants: {s1_1_results['total']:,}")
    print(f"  Phased: {s1_1_results['phased']:,} ({100*s1_1_results['phased']/s1_1_results['total']:.1f}%)")
    print(f"  Het-phased: {s1_1_results['het_phased']:,}")

    print(f"\n[S1.4] Somatic-as-Anchor:")
    print(f"  TP phased as anchor: {s1_4_results['tp_phased']}")
    print(f"  FP phased as anchor: {s1_4_results['fp_phased']}")
    print(f"  PS blocks total: {s1_4_results['blocks_total']:,}")
    print(f"  PS blocks with somatic: {s1_4_results['blocks_with_somatic']:,} ({100*s1_4_results['blocks_with_somatic']/max(s1_4_results['blocks_total'],1):.1f}%)")
    print(f"  Somatic fraction of all block variants: {100*s1_4_results['somatic_in_blocks']/max(s1_4_results['total_in_blocks'],1):.2f}%")

    print(f"\n  All outputs saved to: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
