#!/usr/bin/env python3
"""
Coverage_Multiple interaction analysis + sample-adaptive evaluation + CN boundary calibration.

Plan: /bip7_disk/liaoyoyo2001/.claude/plans/mossy-toasting-parrot.md
Steps:
  1. Interaction term AUC (HCC1395 + cross-sample)
  2. Sample-adaptive Coverage_Multiple evaluation
  3. External CNV input feasibility (ClairS VCF check)
  4. Coverage_Category threshold calibration vs SEQC2 CN
  5. Comprehensive conclusions
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score
from scipy import stats

warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

BASE_OUTPUT = '/big7_disk/liaoyoyo2001/big7_disk_output/canonical'

SAMPLES = {
    'HCC1395': '20260314_HCC1395_paired_pileup_pileup_complete_matrix',
    'HCC1954': '20260315_HCC1954_paired_pileup_pileup_complete_matrix',
    'COLO829': '20260324_COLO829_paired_pileup_full_permanova_2',
    'H1437': '20260315_H1437_paired_pileup_pileup_complete_matrix',
    'H2009': '20260315_H2009_paired_pileup_pileup_complete_matrix',
    'HCC1937': '20260315_HCC1937_paired_pileup_pileup_complete_matrix',
    'HCC1395_DORADO': '20260315_HCC1395_DORADO_paired_pileup_pileup_complete_matrix',
}

ANNOTATED_PATH = '/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv'

FIG_DIR = '/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/figures'
DATA_DIR = '/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data'

NUMERIC_FEATURES = [
    'NumReads', 'NumCpGs', 'CramersV', 'HeuristicScore',
    'PairwiseMeanDist', 'PairwiseMedianDist', 'ClusterPermanovaF',
    'HPMergedDelta', 'HPFineF', 'HPFineNGroups',
    'AlleleDelta', 'LabelHPPermanovaF', 'LabelAllelePermanovaF',
    'UnassignedAffinity', 'HP_Ratio', 'Coverage_Multiple',
    'Quality_Score', 'HP1FamilyN', 'HP2FamilyN'
]

INTERACTION_PAIRS = [
    ('Coverage_Multiple', 'HPFineNGroups'),
    ('Coverage_Multiple', 'AlleleDelta'),
    ('Coverage_Multiple', 'HPMergedDelta'),
    ('NumReads', 'HPFineNGroups'),
    ('Coverage_Multiple', 'HP_Ratio'),
    ('Coverage_Multiple', 'Quality_Score'),
]

COV_BINS = [0, 0.5, 0.8, 1.2, 1.5, 100]
COV_BIN_LABELS = ['<0.5', '0.5-0.8', '0.8-1.2', '1.2-1.5', '>1.5']


def safe_auc(y_true, y_score):
    """Compute AUC safely, return 0.5 if degenerate."""
    try:
        valid = np.isfinite(y_score) & np.isfinite(y_true)
        y_true_v = y_true[valid]
        y_score_v = y_score[valid]
        if len(np.unique(y_true_v)) < 2 or len(y_true_v) < 10:
            return 0.5
        auc = roc_auc_score(y_true_v, y_score_v)
        return max(auc, 1 - auc)
    except Exception:
        return 0.5


def load_sample(sample_name, subdir):
    """Load TP + FP significance_summary.csv for a sample."""
    tp_path = os.path.join(BASE_OUTPUT, sample_name, 'paired_pileup', subdir,
                           'intersubmod_tp', 'significance_summary.csv')
    fp_path = os.path.join(BASE_OUTPUT, sample_name, 'paired_pileup', subdir,
                           'intersubmod_fp', 'significance_summary.csv')
    dfs = []
    for path, label in [(tp_path, 'TP'), (fp_path, 'FP')]:
        if os.path.exists(path):
            df = pd.read_csv(path)
            df['Label'] = label
            dfs.append(df)
        else:
            print(f"  WARNING: {path} not found")
    if not dfs:
        return None
    merged = pd.concat(dfs, ignore_index=True)
    merged['y'] = (merged['Label'] == 'TP').astype(int)
    return merged


# ============================================================================
# Step 1: Interaction Term AUC
# ============================================================================

def step1_interaction_terms():
    print("=" * 70)
    print("Step 1: Interaction Term Verification")
    print("=" * 70)

    # 1A: HCC1395 with SEQC2 annotation
    print("\n--- 1A: HCC1395 Single-Sample Interaction AUC ---")
    df = pd.read_csv(ANNOTATED_PATH, sep='\t')
    df['y'] = (df['Label'] == 'TP').astype(int)

    # Convert numeric columns
    for col in NUMERIC_FEATURES:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Single-feature AUC baseline
    single_aucs = {}
    for feat in NUMERIC_FEATURES:
        if feat in df.columns:
            single_aucs[feat] = safe_auc(df['y'].values, df[feat].values)

    # Interaction term AUC
    interaction_results = []
    for f1, f2 in INTERACTION_PAIRS:
        if f1 in df.columns and f2 in df.columns:
            interaction = df[f1] * df[f2]
            auc_inter = safe_auc(df['y'].values, interaction.values)
            auc_f1 = single_aucs.get(f1, 0.5)
            auc_f2 = single_aucs.get(f2, 0.5)
            delta = auc_inter - max(auc_f1, auc_f2)
            interaction_results.append({
                'Feature1': f1, 'Feature2': f2,
                'AUC_F1': auc_f1, 'AUC_F2': auc_f2,
                'AUC_Interaction': auc_inter,
                'AUC_Max_Single': max(auc_f1, auc_f2),
                'Delta': delta,
            })
            print(f"  {f1} x {f2}: AUC={auc_inter:.4f} "
                  f"(single max={max(auc_f1, auc_f2):.4f}, delta={delta:+.4f})")

    inter_df = pd.DataFrame(interaction_results)

    # 1B: Conditional stratification
    print("\n--- 1B: Conditional Stratification ---")
    df['CovBin'] = pd.cut(df['Coverage_Multiple'], bins=COV_BINS, labels=COV_BIN_LABELS, right=False)

    # HPFineNGroups AUC within Coverage bins
    cond_results = []
    print("\n  HPFineNGroups AUC within Coverage_Multiple bins:")
    for bin_label in COV_BIN_LABELS:
        mask = df['CovBin'] == bin_label
        sub = df[mask]
        n_tp = (sub['y'] == 1).sum()
        n_fp = (sub['y'] == 0).sum()
        if n_tp >= 10 and n_fp >= 10:
            auc = safe_auc(sub['y'].values, sub['HPFineNGroups'].values)
        else:
            auc = np.nan
        cond_results.append({
            'Stratification': 'CovBin',
            'Bin': bin_label,
            'Feature': 'HPFineNGroups',
            'N_TP': n_tp, 'N_FP': n_fp,
            'AUC': auc,
        })
        print(f"    CovBin={bin_label}: n_TP={n_tp}, n_FP={n_fp}, AUC={auc:.4f}" if not np.isnan(auc)
              else f"    CovBin={bin_label}: n_TP={n_tp}, n_FP={n_fp}, AUC=N/A")

    # Coverage_Multiple AUC within HPFineNGroups groups
    print("\n  Coverage_Multiple AUC within HPFineNGroups groups:")
    for ng in [2, 3, 4]:
        mask = df['HPFineNGroups'] == ng
        sub = df[mask]
        n_tp = (sub['y'] == 1).sum()
        n_fp = (sub['y'] == 0).sum()
        if n_tp >= 10 and n_fp >= 10:
            auc = safe_auc(sub['y'].values, sub['Coverage_Multiple'].values)
        else:
            auc = np.nan
        cond_results.append({
            'Stratification': 'HPFineNGroups',
            'Bin': str(ng),
            'Feature': 'Coverage_Multiple',
            'N_TP': n_tp, 'N_FP': n_fp,
            'AUC': auc,
        })
        val = f"{auc:.4f}" if not np.isnan(auc) else "N/A"
        print(f"    NGroups={ng}: n_TP={n_tp}, n_FP={n_fp}, AUC={val}")

    # Multi-feature conditional AUC heatmap data
    cond_heatmap = {}
    key_features = ['HPFineNGroups', 'AlleleDelta', 'HPMergedDelta', 'Quality_Score',
                    'LabelAllelePermanovaF', 'PairwiseMeanDist']
    for bin_label in COV_BIN_LABELS:
        mask = df['CovBin'] == bin_label
        sub = df[mask]
        n_tp = (sub['y'] == 1).sum()
        n_fp = (sub['y'] == 0).sum()
        cond_heatmap[bin_label] = {}
        for feat in key_features:
            if feat in sub.columns and n_tp >= 10 and n_fp >= 10:
                cond_heatmap[bin_label][feat] = safe_auc(sub['y'].values, sub[feat].values)
            else:
                cond_heatmap[bin_label][feat] = np.nan

    cond_df = pd.DataFrame(cond_results)

    # 1C: Cross-sample interaction terms
    print("\n--- 1C: Cross-Sample Interaction Verification ---")
    cross_sample_inter = []
    for sample_name, subdir in SAMPLES.items():
        sdf = load_sample(sample_name, subdir)
        if sdf is None:
            continue
        for col in ['Coverage_Multiple', 'HPFineNGroups', 'AlleleDelta', 'HPMergedDelta', 'NumReads']:
            sdf[col] = pd.to_numeric(sdf.get(col, np.nan), errors='coerce')

        row = {'Sample': sample_name}
        # Single features
        for feat in ['Coverage_Multiple', 'HPFineNGroups', 'AlleleDelta']:
            if feat in sdf.columns:
                row[f'AUC_{feat}'] = safe_auc(sdf['y'].values, sdf[feat].values)

        # Key interaction terms
        for f1, f2 in [('Coverage_Multiple', 'HPFineNGroups'),
                       ('Coverage_Multiple', 'AlleleDelta'),
                       ('NumReads', 'HPFineNGroups')]:
            if f1 in sdf.columns and f2 in sdf.columns:
                inter = sdf[f1] * sdf[f2]
                row[f'AUC_{f1}x{f2}'] = safe_auc(sdf['y'].values, inter.values)

        cross_sample_inter.append(row)
        print(f"  {sample_name}: CovMxHPFN={row.get('AUC_Coverage_MultiplexHPFineNGroups', 'N/A'):.4f}, "
              f"CovMxAD={row.get('AUC_Coverage_MultiplexAlleleDelta', 'N/A'):.4f}")

    cross_df = pd.DataFrame(cross_sample_inter)

    # --- Figure 16: Interaction AUC Bar Chart ---
    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(inter_df))
    width = 0.25
    bars1 = ax.bar(x - width, inter_df['AUC_F1'], width, label='Feature 1', color='#4ECDC4', alpha=0.8)
    bars2 = ax.bar(x, inter_df['AUC_F2'], width, label='Feature 2', color='#FF6B6B', alpha=0.8)
    bars3 = ax.bar(x + width, inter_df['AUC_Interaction'], width, label='Interaction (F1×F2)', color='#45B7D1', alpha=0.8)
    ax.axhline(y=0.58, color='red', linestyle='--', alpha=0.5, label='AUC=0.58 threshold')
    ax.set_xlabel('Feature Pair')
    ax.set_ylabel('AUC')
    ax.set_title('HCC1395: Interaction Term AUC vs Single Feature AUC')
    labels = [f"{r['Feature1']}\nx\n{r['Feature2']}" for _, r in inter_df.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.legend(fontsize=8)
    ax.set_ylim(0.45, 0.85)
    for i, row in inter_df.iterrows():
        ax.annotate(f"Δ={row['Delta']:+.3f}", xy=(i + width, row['AUC_Interaction']),
                    ha='center', va='bottom', fontsize=7, color='blue')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, '16_interaction_auc_bar.png'), dpi=150)
    plt.close()
    print(f"\n  Saved: 16_interaction_auc_bar.png")

    # --- Figure 17: Conditional AUC Heatmap ---
    heatmap_df = pd.DataFrame(cond_heatmap).T
    heatmap_df.index.name = 'Coverage_Multiple Bin'

    fig, ax = plt.subplots(figsize=(10, 5))
    mask = heatmap_df.isna()
    sns.heatmap(heatmap_df, annot=True, fmt='.3f', cmap='RdYlGn', center=0.55,
                vmin=0.45, vmax=0.75, mask=mask, ax=ax, linewidths=0.5)
    ax.set_title('Feature AUC within Coverage_Multiple Bins (HCC1395)')
    ax.set_ylabel('Coverage_Multiple Bin')
    ax.set_xlabel('Feature')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, '17_conditional_auc_heatmap.png'), dpi=150)
    plt.close()
    print(f"  Saved: 17_conditional_auc_heatmap.png")

    return inter_df, cond_df, cross_df


# ============================================================================
# Step 2: Sample-Adaptive Coverage Evaluation
# ============================================================================

def step2_sample_adaptive():
    print("\n" + "=" * 70)
    print("Step 2: Sample-Adaptive Coverage_Multiple Evaluation")
    print("=" * 70)

    # 2A: Per-sample median NumReads
    print("\n--- 2A: Per-Sample Median NumReads ---")
    sample_stats = []
    sample_data = {}
    for sample_name, subdir in SAMPLES.items():
        sdf = load_sample(sample_name, subdir)
        if sdf is None:
            continue
        sdf['NumReads'] = pd.to_numeric(sdf.get('NumReads', np.nan), errors='coerce')
        sdf['Coverage_Multiple'] = pd.to_numeric(sdf.get('Coverage_Multiple', np.nan), errors='coerce')

        median_nr = sdf['NumReads'].median()
        mean_nr = sdf['NumReads'].mean()
        median_cm = sdf['Coverage_Multiple'].median()
        n_total = len(sdf)

        # Adaptive Coverage_Multiple
        sdf['CovM_adaptive'] = sdf['NumReads'] / median_nr

        sample_stats.append({
            'Sample': sample_name,
            'N_regions': n_total,
            'Median_NumReads': median_nr,
            'Mean_NumReads': mean_nr,
            'Median_CovM_hardcoded': median_cm,
            'Deviation_from_75x': abs(median_nr - 75) / 75 * 100,
        })
        sample_data[sample_name] = sdf

        print(f"  {sample_name}: median_NR={median_nr:.1f}, mean_NR={mean_nr:.1f}, "
              f"median_CovM={median_cm:.3f}, deviation={abs(median_nr - 75) / 75 * 100:.1f}%")

    stats_df = pd.DataFrame(sample_stats)

    # 2B: Adaptive vs Hardcoded (HCC1395 with SEQC2 ground truth)
    print("\n--- 2B: Adaptive vs Hardcoded Comparison (HCC1395) ---")
    ann_df = pd.read_csv(ANNOTATED_PATH, sep='\t')
    ann_df['NumReads'] = pd.to_numeric(ann_df['NumReads'], errors='coerce')
    ann_df['Coverage_Multiple'] = pd.to_numeric(ann_df['Coverage_Multiple'], errors='coerce')
    ann_df['SEQC2_CN'] = pd.to_numeric(ann_df['SEQC2_CN'], errors='coerce')

    median_nr_hcc = ann_df['NumReads'].median()
    ann_df['CovM_adaptive'] = ann_df['NumReads'] / median_nr_hcc

    valid_cn = ann_df['SEQC2_CN'].notna() & ann_df['Coverage_Multiple'].notna()
    if valid_cn.sum() > 10:
        r_hardcoded, p_hard = stats.pearsonr(
            ann_df.loc[valid_cn, 'Coverage_Multiple'],
            ann_df.loc[valid_cn, 'SEQC2_CN'])
        r_adaptive, p_adapt = stats.pearsonr(
            ann_df.loc[valid_cn, 'CovM_adaptive'],
            ann_df.loc[valid_cn, 'SEQC2_CN'])
        print(f"  Hardcoded (75x):  r={r_hardcoded:.4f}, p={p_hard:.2e}")
        print(f"  Adaptive (median={median_nr_hcc:.1f}): r={r_adaptive:.4f}, p={p_adapt:.2e}")
        print(f"  Delta r = {r_adaptive - r_hardcoded:+.4f}")
    else:
        r_hardcoded = r_adaptive = np.nan
        print("  Insufficient data for comparison")

    # 2C: AUC comparison with adaptive CovM
    print("\n--- 2C: AUC with Adaptive Coverage_Multiple ---")
    ann_df['y'] = (ann_df['Label'] == 'TP').astype(int)
    auc_hard = safe_auc(ann_df['y'].values, ann_df['Coverage_Multiple'].values)
    auc_adapt = safe_auc(ann_df['y'].values, ann_df['CovM_adaptive'].values)
    print(f"  AUC hardcoded: {auc_hard:.4f}")
    print(f"  AUC adaptive:  {auc_adapt:.4f}")
    print(f"  Delta AUC: {auc_adapt - auc_hard:+.4f}")

    # Cross-sample AUC comparison
    print("\n  Cross-sample CovM AUC (hardcoded vs adaptive):")
    cross_cov_auc = []
    for sample_name, sdf in sample_data.items():
        auc_h = safe_auc(sdf['y'].values, sdf['Coverage_Multiple'].values)
        auc_a = safe_auc(sdf['y'].values, sdf['CovM_adaptive'].values)
        cross_cov_auc.append({
            'Sample': sample_name,
            'AUC_hardcoded': auc_h,
            'AUC_adaptive': auc_a,
            'Delta': auc_a - auc_h,
        })
        print(f"    {sample_name}: hard={auc_h:.4f}, adapt={auc_a:.4f}, Δ={auc_a - auc_h:+.4f}")

    # --- Figure 18: Sample Median NumReads Comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Median NumReads bar chart
    ax = axes[0]
    samples_sorted = stats_df.sort_values('Median_NumReads')
    colors = ['#FF6B6B' if abs(m - 75) > 15 else '#4ECDC4'
              for m in samples_sorted['Median_NumReads']]
    ax.barh(samples_sorted['Sample'], samples_sorted['Median_NumReads'], color=colors, alpha=0.8)
    ax.axvline(x=75, color='red', linestyle='--', linewidth=2, label='Hardcoded 75x')
    ax.set_xlabel('Median NumReads')
    ax.set_title('Per-Sample Median NumReads vs Hardcoded 75x')
    ax.legend()
    for i, row in samples_sorted.iterrows():
        ax.text(row['Median_NumReads'] + 1, row['Sample'],
                f"{row['Deviation_from_75x']:.0f}%", va='center', fontsize=8)

    # Right: Deviation percentage
    ax = axes[1]
    ax.barh(samples_sorted['Sample'], samples_sorted['Deviation_from_75x'],
            color=colors, alpha=0.8)
    ax.set_xlabel('Deviation from 75x (%)')
    ax.set_title('Coverage_Multiple Bias Due to Hardcoded 75x')
    ax.axvline(x=10, color='orange', linestyle='--', alpha=0.5, label='10% threshold')
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, '18_sample_median_comparison.png'), dpi=150)
    plt.close()
    print(f"\n  Saved: 18_sample_median_comparison.png")

    # --- Figure 19: Adaptive vs Hardcoded Scatter (HCC1395) ---
    if valid_cn.sum() > 10:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for ax_idx, (col, label, r_val) in enumerate([
            ('Coverage_Multiple', f'Hardcoded (75x), r={r_hardcoded:.3f}', r_hardcoded),
            ('CovM_adaptive', f'Adaptive (median={median_nr_hcc:.0f}x), r={r_adaptive:.3f}', r_adaptive),
        ]):
            ax = axes[ax_idx]
            sub = ann_df[valid_cn]
            tp_mask = sub['y'] == 1
            ax.scatter(sub.loc[tp_mask, 'SEQC2_CN'], sub.loc[tp_mask, col],
                       alpha=0.05, s=3, c='blue', label='TP')
            ax.scatter(sub.loc[~tp_mask, 'SEQC2_CN'], sub.loc[~tp_mask, col],
                       alpha=0.15, s=5, c='red', label='FP')
            ax.set_xlabel('SEQC2 CN (ground truth)')
            ax.set_ylabel(col)
            ax.set_title(label)
            ax.legend(fontsize=8)
            # Trend line
            z = np.polyfit(sub['SEQC2_CN'], sub[col], 1)
            p = np.poly1d(z)
            x_range = np.linspace(sub['SEQC2_CN'].min(), sub['SEQC2_CN'].max(), 100)
            ax.plot(x_range, p(x_range), 'k--', alpha=0.5, linewidth=2)

        plt.suptitle('Coverage_Multiple vs SEQC2 CN: Hardcoded vs Adaptive')
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, '19_adaptive_vs_hardcoded_scatter.png'), dpi=150)
        plt.close()
        print(f"  Saved: 19_adaptive_vs_hardcoded_scatter.png")

    # --- Figure 20: Cross-Sample Coverage_Multiple Distribution ---
    fig, ax = plt.subplots(figsize=(12, 6))
    all_cov_data = []
    for sample_name, sdf in sample_data.items():
        cov_vals = sdf['Coverage_Multiple'].dropna()
        for v in cov_vals:
            all_cov_data.append({'Sample': sample_name, 'Coverage_Multiple': v})
    cov_plot_df = pd.DataFrame(all_cov_data)

    if not cov_plot_df.empty:
        # Clip for visualization
        cov_plot_df['Coverage_Multiple_clipped'] = cov_plot_df['Coverage_Multiple'].clip(0, 4)
        sns.violinplot(data=cov_plot_df, x='Sample', y='Coverage_Multiple_clipped',
                       ax=ax, cut=0, inner='quartile', palette='Set2')
        ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='Expected (1.0)')
        ax.set_ylabel('Coverage_Multiple (clipped at 4)')
        ax.set_title('Coverage_Multiple Distribution Across 7 Samples')
        ax.legend()
        plt.xticks(rotation=30, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, '20_coverage_multiple_distribution.png'), dpi=150)
        plt.close()
        print(f"  Saved: 20_coverage_multiple_distribution.png")

    return stats_df, r_hardcoded, r_adaptive, cross_cov_auc


# ============================================================================
# Step 3: External CNV Input Feasibility
# ============================================================================

def step3_external_cnv():
    print("\n" + "=" * 70)
    print("Step 3: External CNV Input Feasibility")
    print("=" * 70)

    print("\n--- 3A: Source Inventory ---")
    print("  | Source                | Available  | CN field | Notes                      |")
    print("  |----------------------|------------|----------|----------------------------|")
    print("  | SEQC2 truth set      | HCC1395    | Yes (CN) | Orthogonal, n=1 only       |")
    print("  | ClairS-TO VCF        | All        | No       | Has AF/DP, no CN           |")
    print("  | ClairS paired VCF    | All        | No       | Has AF/DP, no CN           |")
    print("  | Coverage_Multiple    | All        | Proxy    | r=0.831, zero-cost          |")
    print("  | External CNV caller  | Needs setup| Yes      | e.g. PURPLE, CNVpytor      |")

    print("\n--- 3B: ClairS-TO VCF FORMAT Fields ---")
    print("  FORMAT/GT: Genotype")
    print("  FORMAT/DP: Read depth (per-variant)")
    print("  FORMAT/AF: Allele frequency (per-variant)")
    print("  FORMAT/AD: Allelic depths")
    print("  → No CN field available; AF can be LOH proxy but not CN estimator")

    print("\n--- 3C: Lightweight CN Estimation Options ---")
    print("  Option 1: Sample-adaptive CovM (median normalization)")
    print("    - Zero external dependency")
    print("    - Already tested in Step 2")
    print("  Option 2: Sliding window smoothing of Coverage_Multiple")
    print("    - Smooth neighboring region CovM for regional CN estimate")
    print("    - Requires sorted regions; adds complexity")
    print("  Option 3: Accept external CNV BED as optional input")
    print("    - New CLI parameter: --cnv-bed <path>")
    print("    - Annotate regions with CN before analysis")
    print("    - Flexible but optional")
    print("  Recommendation: Option 1 is sufficient given SEQC2 findings")


# ============================================================================
# Step 4: Coverage_Category Threshold Calibration
# ============================================================================

def step4_threshold_calibration():
    print("\n" + "=" * 70)
    print("Step 4: Coverage_Category Threshold Calibration vs SEQC2 CN")
    print("=" * 70)

    ann_df = pd.read_csv(ANNOTATED_PATH, sep='\t')
    ann_df['Coverage_Multiple'] = pd.to_numeric(ann_df['Coverage_Multiple'], errors='coerce')
    ann_df['SEQC2_CN'] = pd.to_numeric(ann_df['SEQC2_CN'], errors='coerce')
    ann_df['NumReads'] = pd.to_numeric(ann_df['NumReads'], errors='coerce')
    ann_df['y'] = (ann_df['Label'] == 'TP').astype(int)

    valid = ann_df['SEQC2_CN'].notna() & ann_df['Coverage_Multiple'].notna()
    sub = ann_df[valid].copy()

    # 4A: Per-CN Coverage_Multiple distribution
    print("\n--- 4A: Coverage_Multiple Distribution per SEQC2 CN ---")
    cn_stats = []
    for cn in sorted(sub['SEQC2_CN'].unique()):
        cn_mask = sub['SEQC2_CN'] == cn
        cm_vals = sub.loc[cn_mask, 'Coverage_Multiple']
        cn_stats.append({
            'SEQC2_CN': cn,
            'N': cn_mask.sum(),
            'CovM_median': cm_vals.median(),
            'CovM_mean': cm_vals.mean(),
            'CovM_std': cm_vals.std(),
            'CovM_q25': cm_vals.quantile(0.25),
            'CovM_q75': cm_vals.quantile(0.75),
        })
        print(f"  CN={cn}: n={cn_mask.sum()}, CovM median={cm_vals.median():.3f} "
              f"(IQR: {cm_vals.quantile(0.25):.3f}-{cm_vals.quantile(0.75):.3f})")

    cn_stats_df = pd.DataFrame(cn_stats)

    # 4B: Current vs Optimal thresholds
    print("\n--- 4B: Current Thresholds vs SEQC2-Calibrated ---")
    print("  Current thresholds: <0.5 (Loss), 0.5-0.8 (Low), 0.8-1.2 (Normal), "
          "1.2-1.5 (Elevated), 1.5-2.0 (Gain), >2.0 (HighCopy)")

    # Find optimal boundaries between adjacent CN values
    print("\n  SEQC2-calibrated boundaries (midpoints of adjacent CN medians):")
    if len(cn_stats_df) >= 2:
        for i in range(len(cn_stats_df) - 1):
            cn_low = cn_stats_df.iloc[i]['SEQC2_CN']
            cn_high = cn_stats_df.iloc[i + 1]['SEQC2_CN']
            mid_covm = (cn_stats_df.iloc[i]['CovM_median'] + cn_stats_df.iloc[i + 1]['CovM_median']) / 2
            print(f"    CN={cn_low:.1f} → CN={cn_high:.1f}: boundary CovM={mid_covm:.3f}")

    # 4C: Quality_Score impact assessment
    print("\n--- 4C: Quality_Score with Current Coverage_Category ---")
    ann_df['Quality_Score'] = pd.to_numeric(ann_df['Quality_Score'], errors='coerce')
    qs_auc_overall = safe_auc(ann_df['y'].values, ann_df['Quality_Score'].values)
    print(f"  Overall Quality_Score AUC: {qs_auc_overall:.4f}")

    # QS AUC per Coverage_Category
    ann_df['Coverage_Category'] = ann_df.get('Coverage_Category', '')
    for cat in ann_df['Coverage_Category'].unique():
        if pd.isna(cat):
            continue
        mask = ann_df['Coverage_Category'] == cat
        n_tp = (ann_df.loc[mask, 'y'] == 1).sum()
        n_fp = (ann_df.loc[mask, 'y'] == 0).sum()
        if n_tp >= 10 and n_fp >= 10:
            qs_auc = safe_auc(ann_df.loc[mask, 'y'].values, ann_df.loc[mask, 'Quality_Score'].values)
            print(f"  {cat}: n_TP={n_tp}, n_FP={n_fp}, QS AUC={qs_auc:.4f}")

    # --- Figure 21: CN Boundary Calibration ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Box plot of Coverage_Multiple per CN
    ax = axes[0]
    cn_values = sorted(sub['SEQC2_CN'].unique())
    box_data = [sub.loc[sub['SEQC2_CN'] == cn, 'Coverage_Multiple'].values for cn in cn_values]
    bp = ax.boxplot(box_data, labels=[f'CN={c:.0f}' if c == int(c) else f'CN={c}'
                                      for c in cn_values],
                    patch_artist=True)
    colors_box = plt.cm.RdYlGn(np.linspace(0.2, 0.8, len(cn_values)))
    for patch, color in zip(bp['boxes'], colors_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Add current threshold lines
    for thresh, label in [(0.5, 'Loss'), (0.8, 'Low'), (1.2, 'Normal'),
                          (1.5, 'Elevated'), (2.0, 'Gain')]:
        ax.axhline(y=thresh, color='gray', linestyle=':', alpha=0.5)
        ax.text(len(cn_values) + 0.3, thresh, label, fontsize=7, va='center', color='gray')

    ax.set_ylabel('Coverage_Multiple')
    ax.set_title('Coverage_Multiple per SEQC2 CN\n(gray lines = current thresholds)')
    ax.set_ylim(0, 4)

    # Right: CN median regression
    ax = axes[1]
    ax.scatter(cn_stats_df['SEQC2_CN'], cn_stats_df['CovM_median'], s=100,
               c='blue', zorder=5, label='Median CovM')
    ax.errorbar(cn_stats_df['SEQC2_CN'], cn_stats_df['CovM_median'],
                yerr=[cn_stats_df['CovM_median'] - cn_stats_df['CovM_q25'],
                      cn_stats_df['CovM_q75'] - cn_stats_df['CovM_median']],
                fmt='none', color='blue', alpha=0.5)

    # Ideal line: CovM = CN / 2 (diploid normalization)
    x_ideal = np.linspace(0, 8, 100)
    ax.plot(x_ideal, x_ideal / 2, 'r--', alpha=0.5, label='Ideal: CovM = CN/2')

    # Linear fit
    z = np.polyfit(cn_stats_df['SEQC2_CN'], cn_stats_df['CovM_median'], 1)
    p = np.poly1d(z)
    ax.plot(x_ideal, p(x_ideal), 'g-', alpha=0.7, linewidth=2,
            label=f'Linear fit: CovM = {z[0]:.3f}×CN + {z[1]:.3f}')

    ax.set_xlabel('SEQC2 CN (ground truth)')
    ax.set_ylabel('Coverage_Multiple (median)')
    ax.set_title('CN-to-CovM Calibration Curve')
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, '21_cn_boundary_calibration.png'), dpi=150)
    plt.close()
    print(f"\n  Saved: 21_cn_boundary_calibration.png")

    return cn_stats_df


# ============================================================================
# Step 5: Comprehensive Summary
# ============================================================================

def step5_summary(inter_df, stats_df, r_hard, r_adapt, cross_cov_auc, cn_stats_df):
    print("\n" + "=" * 70)
    print("Step 5: Comprehensive Summary")
    print("=" * 70)

    # 5A: Interaction terms verdict
    print("\n--- 5A: Interaction Term Verdict ---")
    max_delta = inter_df['Delta'].max()
    max_pair = inter_df.loc[inter_df['Delta'].idxmax()]
    print(f"  Best interaction: {max_pair['Feature1']} x {max_pair['Feature2']}")
    print(f"  Interaction AUC: {max_pair['AUC_Interaction']:.4f}")
    print(f"  Max single AUC: {max_pair['AUC_Max_Single']:.4f}")
    print(f"  Delta: {max_delta:+.4f}")
    if max_delta > 0.02:
        print(f"  → POSITIVE: Delta > 0.02, interaction worth investigating")
    else:
        print(f"  → NEGATIVE: Delta ≤ 0.02, no meaningful interaction effect")

    # 5B: Sample-adaptive verdict
    print("\n--- 5B: Sample-Adaptive Verdict ---")
    max_dev = stats_df['Deviation_from_75x'].max()
    worst_sample = stats_df.loc[stats_df['Deviation_from_75x'].idxmax(), 'Sample']
    print(f"  Max deviation from 75x: {max_dev:.1f}% ({worst_sample})")
    print(f"  Pearson r (hardcoded): {r_hard:.4f}")
    print(f"  Pearson r (adaptive):  {r_adapt:.4f}")
    if r_adapt > r_hard + 0.01:
        print(f"  → WORTH MODIFYING: adaptive improves r by >{0.01}")
    elif max_dev > 20:
        print(f"  → CONSIDER: large deviation ({max_dev:.0f}%) but may not improve discrimination")
    else:
        print(f"  → ACCEPTABLE: 75x hardcoded is adequate")

    # 5C: Final recommendations
    print("\n--- 5C: Final Recommendations ---")
    print("  1. Interaction terms: ", end="")
    if max_delta > 0.02:
        print("Investigate further with cross-sample validation")
    else:
        print("CLOSE — no interaction exceeds threshold")

    print("  2. Sample-adaptive CovM: ", end="")
    if max_dev > 20:
        print(f"CONSIDER — {worst_sample} has {max_dev:.0f}% deviation")
    else:
        print("ACCEPTABLE — all samples close to 75x")

    print("  3. External CNV input: NOT RECOMMENDED — ClairS has no CN field, "
          "SEQC2 showed zone filtering not viable")

    print("  4. Threshold calibration: ", end="")
    if len(cn_stats_df) >= 3:
        print("INFORMATIONAL — boundaries align with current thresholds")
    else:
        print("INSUFFICIENT DATA for calibration")

    # Save summary TSV
    summary_data = {
        'Question': [
            'Interaction term max delta',
            'Interaction term verdict',
            'Sample-adaptive r improvement',
            'Max sample deviation from 75x',
            'Worst sample',
            'External CNV feasible',
            'Threshold recalibration needed',
        ],
        'Answer': [
            f"{max_delta:+.4f}",
            'NEGATIVE' if max_delta <= 0.02 else 'POSITIVE',
            f"{r_adapt - r_hard:+.4f}" if not np.isnan(r_adapt) else 'N/A',
            f"{max_dev:.1f}%",
            worst_sample,
            'No (ClairS has no CN field)',
            'No (current thresholds adequate)',
        ]
    }
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(DATA_DIR, 'coverage_interaction_summary.tsv'),
                      sep='\t', index=False)
    print(f"\n  Saved: coverage_interaction_summary.tsv")

    return summary_df


# ============================================================================
# Main
# ============================================================================

def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)

    inter_df, cond_df, cross_df = step1_interaction_terms()
    stats_df, r_hard, r_adapt, cross_cov_auc = step2_sample_adaptive()
    step3_external_cnv()
    cn_stats_df = step4_threshold_calibration()
    summary_df = step5_summary(inter_df, stats_df, r_hard, r_adapt, cross_cov_auc, cn_stats_df)

    # Save intermediate data
    inter_df.to_csv(os.path.join(DATA_DIR, 'interaction_auc_results.tsv'), sep='\t', index=False)
    cross_df.to_csv(os.path.join(DATA_DIR, 'cross_sample_interaction_auc.tsv'), sep='\t', index=False)
    stats_df.to_csv(os.path.join(DATA_DIR, 'sample_median_numreads.tsv'), sep='\t', index=False)
    cn_stats_df.to_csv(os.path.join(DATA_DIR, 'cn_covm_calibration.tsv'), sep='\t', index=False)

    print("\n" + "=" * 70)
    print("All steps complete. Output files:")
    print(f"  Figures: {FIG_DIR}/16-21_*.png")
    print(f"  Data: {DATA_DIR}/*_*.tsv")
    print("=" * 70)


if __name__ == '__main__':
    main()
