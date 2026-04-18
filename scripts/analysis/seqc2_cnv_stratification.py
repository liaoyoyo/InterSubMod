#!/usr/bin/env python3
"""
SEQC2 CNV Stratification Analysis for HCC1395
===============================================
Annotate ISM regions with SEQC2 CNV zones and LongPhase-TO LOH,
then analyze feature distributions and TP/FP discrimination per zone.

Usage:
    python scripts/analysis/seqc2_cnv_stratification.py

Output:
    research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv
    research/seqc2_cnv_stratification/figures/*.png
"""

import os
import sys
import warnings
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / "research" / "seqc2_cnv_stratification"
FIG_DIR = OUTPUT_DIR / "figures"
DATA_DIR = OUTPUT_DIR / "data"

SEQC2_COMBINED = Path("/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed")
SEQC2_GAIN_CN = Path("/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_gain_cn.bed")
SEQC2_LOSS_CN = Path("/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_loss_cn.bed")

LONGPHASE_LOH = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
    "20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed"
)

CANONICAL_BASE = Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_pileup")

# Find the actual subdirectory
_subdirs = list(CANONICAL_BASE.glob("*complete_matrix"))
if _subdirs:
    ISM_TP = _subdirs[0] / "intersubmod_tp" / "significance_summary.csv"
    ISM_FP = _subdirs[0] / "intersubmod_fp" / "significance_summary.csv"
else:
    print("ERROR: Cannot find ISM data directory under", CANONICAL_BASE)
    sys.exit(1)

# ── Style ──────────────────────────────────────────────────────────────────────
ZONE_ORDER = ["Neutral", "Gain-only", "LOH-only", "Gain+LOH", "Loss-only", "Loss+LOH"]
ZONE_COLORS = {
    "Neutral": "#4CAF50",
    "Gain-only": "#FF9800",
    "LOH-only": "#9C27B0",
    "Gain+LOH": "#F44336",
    "Loss-only": "#2196F3",
    "Loss+LOH": "#795548",
}

KEY_FEATURES = [
    "NumReads", "Coverage_Multiple", "HP_Ratio", "HPFineNGroups",
    "Quality_Score", "PairwiseMeanDist", "AlleleDelta", "HPMergedDelta",
]

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.dpi': 150,
})


# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Data Preparation
# ══════════════════════════════════════════════════════════════════════════════

def load_seqc2_cnv():
    """Load SEQC2 CNV BED files and build non-overlapping zone BEDs."""
    combined = pd.read_csv(SEQC2_COMBINED, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'type'])
    gain_cn = pd.read_csv(SEQC2_GAIN_CN, sep='\t', header=None,
                          names=['chr', 'start', 'end', 'cn'])
    loss_cn = pd.read_csv(SEQC2_LOSS_CN, sep='\t', header=None,
                          names=['chr', 'start', 'end', 'cn'])
    print(f"  SEQC2 combined: {len(combined)} segments "
          f"(gain={sum(combined.type=='gain')}, loss={sum(combined.type=='loss')}, "
          f"loh={sum(combined.type=='loh')})")
    print(f"  SEQC2 gain+CN: {len(gain_cn)} segments (CN range: {gain_cn.cn.min()}-{gain_cn.cn.max()})")
    print(f"  SEQC2 loss+CN: {len(loss_cn)} segments (CN range: {loss_cn.cn.min()}-{loss_cn.cn.max()})")
    return combined, gain_cn, loss_cn


def load_longphase_loh():
    """Load LongPhase-TO LOH.bed."""
    loh = pd.read_csv(LONGPHASE_LOH, sep='\t', header=None, usecols=[0, 1, 2],
                      names=['chr', 'start', 'end'])
    print(f"  LongPhase LOH: {len(loh)} regions")
    return loh


def load_ism_data():
    """Load HCC1395 paired ISM significance_summary.csv for TP and FP."""
    tp = pd.read_csv(ISM_TP)
    tp['Label'] = 'TP'
    fp = pd.read_csv(ISM_FP)
    fp['Label'] = 'FP'
    df = pd.concat([tp, fp], ignore_index=True)
    print(f"  ISM data: {len(tp)} TP + {len(fp)} FP = {len(df)} total")
    return df


def annotate_with_bedtools(ism_df, seqc2_combined, longphase_loh, gain_cn, loss_cn):
    """
    Annotate each ISM region with SEQC2 CNV zone and LongPhase LOH status.
    Uses bedtools intersect via subprocess for speed.
    """
    # Create ISM BED (chr, pos, pos+1)
    ism_bed = ism_df[['Chr', 'Pos']].copy()
    ism_bed['end'] = ism_bed['Pos'] + 1
    ism_bed['idx'] = range(len(ism_bed))
    ism_bed = ism_bed[['Chr', 'Pos', 'end', 'idx']]

    with tempfile.TemporaryDirectory() as tmpdir:
        ism_path = os.path.join(tmpdir, 'ism.bed')
        ism_bed.to_csv(ism_path, sep='\t', header=False, index=False)

        # ── SEQC2 gain/loss/loh annotation ──
        # Write separate BEDs for gain, loss, loh from combined
        for cnv_type in ['gain', 'loss', 'loh']:
            subset = seqc2_combined[seqc2_combined.type == cnv_type]
            bed_path = os.path.join(tmpdir, f'{cnv_type}.bed')
            subset[['chr', 'start', 'end']].to_csv(bed_path, sep='\t', header=False, index=False)

        # Intersect each type
        annotations = {}
        for cnv_type in ['gain', 'loss', 'loh']:
            bed_path = os.path.join(tmpdir, f'{cnv_type}.bed')
            result = subprocess.run(
                ['bedtools', 'intersect', '-a', ism_path, '-b', bed_path, '-wa', '-u'],
                capture_output=True, text=True
            )
            hit_indices = set()
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    hit_indices.add(int(parts[3]))
            annotations[cnv_type] = hit_indices
            print(f"    SEQC2 {cnv_type}: {len(hit_indices)} ISM regions overlap")

        # ── LongPhase LOH annotation ──
        loh_path = os.path.join(tmpdir, 'longphase_loh.bed')
        longphase_loh[['chr', 'start', 'end']].to_csv(loh_path, sep='\t', header=False, index=False)
        result = subprocess.run(
            ['bedtools', 'intersect', '-a', ism_path, '-b', loh_path, '-wa', '-u'],
            capture_output=True, text=True
        )
        lp_loh_indices = set()
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                lp_loh_indices.add(int(parts[3]))
        print(f"    LongPhase LOH: {len(lp_loh_indices)} ISM regions overlap")

        # ── SEQC2 CN annotation (gain_cn + loss_cn) ──
        # Write gain_cn and loss_cn as BED with CN in col4
        gain_cn_path = os.path.join(tmpdir, 'gain_cn.bed')
        gain_cn[['chr', 'start', 'end', 'cn']].to_csv(gain_cn_path, sep='\t', header=False, index=False)
        loss_cn_path = os.path.join(tmpdir, 'loss_cn.bed')
        loss_cn[['chr', 'start', 'end', 'cn']].to_csv(loss_cn_path, sep='\t', header=False, index=False)

        # Get CN for each ISM region (use -wb to get the CN value)
        cn_map = {}
        for cn_path in [gain_cn_path, loss_cn_path]:
            result = subprocess.run(
                ['bedtools', 'intersect', '-a', ism_path, '-b', cn_path, '-wa', '-wb'],
                capture_output=True, text=True
            )
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    idx = int(parts[3])
                    cn = float(parts[7])
                    cn_map[idx] = cn

    # ── Assign zones ──
    zones = []
    for i in range(len(ism_df)):
        in_gain = i in annotations['gain']
        in_loss = i in annotations['loss']
        in_loh = i in annotations['loh']

        if in_gain and in_loh:
            zones.append("Gain+LOH")
        elif in_loss and in_loh:
            zones.append("Loss+LOH")
        elif in_gain:
            zones.append("Gain-only")
        elif in_loss:
            zones.append("Loss-only")
        elif in_loh:
            zones.append("LOH-only")
        else:
            zones.append("Neutral")

    ism_df = ism_df.copy()
    ism_df['SEQC2_Zone'] = zones
    ism_df['LongPhase_LOH'] = [i in lp_loh_indices for i in range(len(ism_df))]
    ism_df['SEQC2_CN'] = [cn_map.get(i, 2.0) for i in range(len(ism_df))]  # default CN=2 for neutral

    return ism_df


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Descriptive Statistics
# ══════════════════════════════════════════════════════════════════════════════

def plot_region_distribution(df):
    """Fig 1: Region distribution by SEQC2 zone and TP/FP."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # 1a: Count by zone
    counts = df.groupby(['SEQC2_Zone', 'Label']).size().unstack(fill_value=0)
    counts = counts.reindex(ZONE_ORDER)
    counts[['TP', 'FP']].plot(kind='bar', stacked=True, ax=axes[0],
                               color=['#4CAF50', '#F44336'])
    axes[0].set_title('ISM Region Count by SEQC2 Zone')
    axes[0].set_xlabel('SEQC2 Zone')
    axes[0].set_ylabel('Count')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].legend(title='Label')

    # 1b: FP rate by zone
    fp_rates = {}
    for zone in ZONE_ORDER:
        zdf = df[df.SEQC2_Zone == zone]
        n_tp = sum(zdf.Label == 'TP')
        n_fp = sum(zdf.Label == 'FP')
        total = n_tp + n_fp
        fp_rates[zone] = n_fp / total if total > 0 else 0
    overall_fp_rate = sum(df.Label == 'FP') / len(df)

    bars = axes[1].bar(range(len(ZONE_ORDER)),
                       [fp_rates[z] for z in ZONE_ORDER],
                       color=[ZONE_COLORS[z] for z in ZONE_ORDER])
    axes[1].axhline(y=overall_fp_rate, color='red', linestyle='--', alpha=0.7,
                    label=f'Overall: {overall_fp_rate:.3f}')
    axes[1].set_xticks(range(len(ZONE_ORDER)))
    axes[1].set_xticklabels(ZONE_ORDER, rotation=45, ha='right')
    axes[1].set_title('FP Rate by SEQC2 Zone')
    axes[1].set_ylabel('FP Rate')
    axes[1].legend()
    for i, (z, rate) in enumerate(zip(ZONE_ORDER, [fp_rates[z] for z in ZONE_ORDER])):
        axes[1].text(i, rate + 0.002, f'{rate:.3f}', ha='center', va='bottom', fontsize=8)

    # 1c: Coverage_Category vs SEQC2 Zone confusion matrix
    if 'Coverage_Category' in df.columns:
        ct = pd.crosstab(df['SEQC2_Zone'], df['Coverage_Category'], normalize='index')
        ct = ct.reindex(ZONE_ORDER)
        ct.plot(kind='bar', stacked=True, ax=axes[2], colormap='Set3')
        axes[2].set_title('ISM Coverage_Category within SEQC2 Zones')
        axes[2].set_xlabel('SEQC2 Zone')
        axes[2].set_ylabel('Proportion')
        axes[2].tick_params(axis='x', rotation=45)
        axes[2].legend(fontsize=7, loc='upper right')

    plt.tight_layout()
    fig.savefig(FIG_DIR / '01_region_distribution.png')
    plt.close()
    print("  Saved: 01_region_distribution.png")
    return fp_rates


def plot_feature_violins(df):
    """Fig 3: Key feature distributions by SEQC2 zone (violin plots)."""
    features_to_plot = [f for f in KEY_FEATURES if f in df.columns]
    n_features = len(features_to_plot)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    for i, feat in enumerate(features_to_plot):
        ax = axes[i]
        plot_df = df[['SEQC2_Zone', feat, 'Label']].dropna()
        if len(plot_df) == 0:
            ax.set_title(f'{feat} (no data)')
            continue

        # Clip extreme outliers for visualization
        q01, q99 = plot_df[feat].quantile([0.01, 0.99])
        plot_df = plot_df[(plot_df[feat] >= q01) & (plot_df[feat] <= q99)]

        palette = [ZONE_COLORS[z] for z in ZONE_ORDER if z in plot_df.SEQC2_Zone.unique()]
        order = [z for z in ZONE_ORDER if z in plot_df.SEQC2_Zone.unique()]
        sns.violinplot(data=plot_df, x='SEQC2_Zone', y=feat, order=order,
                       palette=ZONE_COLORS, ax=ax, inner='quartile', cut=0)
        ax.set_title(feat)
        ax.set_xlabel('')
        ax.tick_params(axis='x', rotation=45)

    for i in range(n_features, len(axes)):
        axes[i].set_visible(False)

    plt.suptitle('ISM Feature Distributions by SEQC2 CNV Zone (HCC1395)', fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / '03_feature_violin_by_zone.png')
    plt.close()
    print("  Saved: 03_feature_violin_by_zone.png")


def compute_zone_aucs(df):
    """Fig 4: AUC heatmap — features × zones."""
    numeric_features = [
        'NumReads', 'NumCpGs', 'CramersV', 'HeuristicScore',
        'PairwiseMeanDist', 'PairwiseMedianDist',
        'ClusterPermanovaF', 'HPMergedDelta', 'HPFineF', 'HPFineNGroups',
        'AlleleDelta', 'LabelHPPermanovaF', 'LabelAllelePermanovaF',
        'UnassignedAffinity', 'HP_Ratio', 'Coverage_Multiple',
        'Quality_Score', 'HP1FamilyN', 'HP2FamilyN',
    ]
    numeric_features = [f for f in numeric_features if f in df.columns]

    results = {}
    for zone in ZONE_ORDER + ['All']:
        if zone == 'All':
            zdf = df
        else:
            zdf = df[df.SEQC2_Zone == zone]

        y = (zdf.Label == 'FP').astype(int).values
        if len(set(y)) < 2 or len(y) < 20:
            results[zone] = {f: np.nan for f in numeric_features}
            continue

        zone_aucs = {}
        for feat in numeric_features:
            vals = zdf[feat].values
            mask = ~np.isnan(vals)
            if mask.sum() < 20 or len(set(y[mask])) < 2:
                zone_aucs[feat] = np.nan
            else:
                try:
                    auc = roc_auc_score(y[mask], vals[mask])
                    zone_aucs[feat] = max(auc, 1 - auc)  # direction-agnostic
                except:
                    zone_aucs[feat] = np.nan
        results[zone] = zone_aucs

    auc_df = pd.DataFrame(results).T
    auc_df = auc_df.reindex(ZONE_ORDER + ['All'])

    # Save as TSV
    auc_df.to_csv(DATA_DIR / 'zone_auc_matrix.tsv', sep='\t', float_format='%.4f')

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(16, 6))
    mask = auc_df.isna()
    sns.heatmap(auc_df.astype(float), annot=True, fmt='.3f', cmap='RdYlGn',
                vmin=0.45, vmax=0.75, mask=mask, ax=ax,
                linewidths=0.5, cbar_kws={'label': 'AUC (direction-agnostic)'})
    ax.set_title('TP/FP Discrimination AUC by SEQC2 Zone × ISM Feature (HCC1395)')
    ax.set_ylabel('SEQC2 Zone')
    ax.set_xlabel('ISM Feature')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(FIG_DIR / '04_auc_heatmap_by_zone.png')
    plt.close()
    print("  Saved: 04_auc_heatmap_by_zone.png")

    return auc_df


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: SEQC2 × LongPhase-TO Cross-Reference
# ══════════════════════════════════════════════════════════════════════════════

def cross_reference_analysis(df):
    """Analyze 4 cross-reference zones: SEQC2 LOH × LongPhase LOH."""
    # Define cross-reference zones
    seqc2_loh = df.SEQC2_Zone.isin(['LOH-only', 'Gain+LOH', 'Loss+LOH'])
    lp_loh = df.LongPhase_LOH

    conditions = [
        (seqc2_loh & lp_loh, 'Both-LOH'),
        (seqc2_loh & ~lp_loh, 'SEQC2-only'),
        (~seqc2_loh & lp_loh, 'LP-only'),
        (~seqc2_loh & ~lp_loh, 'Neither'),
    ]
    df = df.copy()
    df['CrossZone'] = 'Neither'
    for mask, label in conditions:
        df.loc[mask, 'CrossZone'] = label

    cross_order = ['Both-LOH', 'SEQC2-only', 'LP-only', 'Neither']
    cross_colors = {'Both-LOH': '#9C27B0', 'SEQC2-only': '#FF9800',
                    'LP-only': '#F44336', 'Neither': '#4CAF50'}

    # ── Fig 5: Cross-reference counts (pseudo Venn) ──
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # 5a: Count bar
    counts = df.groupby('CrossZone').size()
    counts = counts.reindex(cross_order, fill_value=0)
    bars = axes[0].bar(range(len(cross_order)), counts.values,
                       color=[cross_colors[z] for z in cross_order])
    axes[0].set_xticks(range(len(cross_order)))
    axes[0].set_xticklabels(cross_order, rotation=30)
    axes[0].set_title('ISM Regions by SEQC2×LongPhase LOH Cross-Zone')
    axes[0].set_ylabel('Count')
    for i, v in enumerate(counts.values):
        axes[0].text(i, v + 50, str(v), ha='center', fontsize=9)

    # 5b: HP_Ratio distribution by cross-zone
    if 'HP_Ratio' in df.columns:
        for zone in cross_order:
            zdf = df[df.CrossZone == zone]
            vals = zdf['HP_Ratio'].dropna()
            if len(vals) > 0:
                axes[1].hist(vals, bins=50, alpha=0.5, label=f'{zone} (n={len(vals)})',
                             color=cross_colors[zone], density=True)
        axes[1].set_title('HP_Ratio by Cross-Zone')
        axes[1].set_xlabel('HP_Ratio')
        axes[1].set_ylabel('Density')
        axes[1].legend(fontsize=8)

    # 5c: FP rate by cross-zone
    fp_rates = {}
    for zone in cross_order:
        zdf = df[df.CrossZone == zone]
        n_fp = sum(zdf.Label == 'FP')
        n_total = len(zdf)
        fp_rates[zone] = n_fp / n_total if n_total > 0 else 0

    axes[2].bar(range(len(cross_order)),
                [fp_rates[z] for z in cross_order],
                color=[cross_colors[z] for z in cross_order])
    axes[2].set_xticks(range(len(cross_order)))
    axes[2].set_xticklabels(cross_order, rotation=30)
    axes[2].set_title('FP Rate by Cross-Zone')
    axes[2].set_ylabel('FP Rate')
    for i, z in enumerate(cross_order):
        axes[2].text(i, fp_rates[z] + 0.002, f'{fp_rates[z]:.3f}', ha='center', fontsize=9)

    plt.tight_layout()
    fig.savefig(FIG_DIR / '05_cross_reference_analysis.png')
    plt.close()
    print("  Saved: 05_cross_reference_analysis.png")

    # ── Fig 6: Key features by cross-zone ──
    cross_features = ['HP_Ratio', 'HPFineNGroups', 'Coverage_Multiple', 'Quality_Score']
    cross_features = [f for f in cross_features if f in df.columns]

    fig, axes = plt.subplots(1, len(cross_features), figsize=(5 * len(cross_features), 5))
    if len(cross_features) == 1:
        axes = [axes]

    for i, feat in enumerate(cross_features):
        plot_df = df[['CrossZone', feat]].dropna()
        q01, q99 = plot_df[feat].quantile([0.01, 0.99])
        plot_df = plot_df[(plot_df[feat] >= q01) & (plot_df[feat] <= q99)]
        sns.violinplot(data=plot_df, x='CrossZone', y=feat, order=cross_order,
                       palette=cross_colors, ax=axes[i], inner='quartile', cut=0)
        axes[i].set_title(feat)
        axes[i].tick_params(axis='x', rotation=30)

    plt.suptitle('Feature Distributions by SEQC2×LongPhase Cross-Zone', fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / '06_cross_zone_features.png')
    plt.close()
    print("  Saved: 06_cross_zone_features.png")

    return df, fp_rates


# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Coverage_Multiple vs SEQC2 CN
# ══════════════════════════════════════════════════════════════════════════════

def coverage_cn_analysis(df):
    """Fig 2: Coverage_Multiple vs SEQC2 CN scatter plot."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 2a: Scatter — Coverage_Multiple vs SEQC2_CN
    for label, color, marker in [('TP', '#4CAF50', 'o'), ('FP', '#F44336', 'x')]:
        ldf = df[df.Label == label]
        axes[0].scatter(ldf.SEQC2_CN, ldf.Coverage_Multiple,
                        alpha=0.1, s=5, color=color, marker=marker, label=label)

    axes[0].set_xlabel('SEQC2 Copy Number')
    axes[0].set_ylabel('ISM Coverage_Multiple')
    axes[0].set_title('Coverage_Multiple vs SEQC2 CN')
    axes[0].legend()

    # Add correlation
    valid = df[['SEQC2_CN', 'Coverage_Multiple']].dropna()
    if len(valid) > 10:
        corr = valid.SEQC2_CN.corr(valid.Coverage_Multiple)
        axes[0].text(0.05, 0.95, f'Pearson r = {corr:.3f}',
                     transform=axes[0].transAxes, fontsize=10,
                     verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # 2b: Box plot — Coverage_Multiple by CN bin
    cn_bins = [0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.5]
    cn_labels = ['CN=0', 'CN=1', 'CN=2', 'CN=3', 'CN=4', 'CN=5', 'CN≥6']
    df_plot = df.copy()
    df_plot['CN_Bin'] = pd.cut(df_plot.SEQC2_CN, bins=cn_bins, labels=cn_labels, right=True)
    df_plot = df_plot.dropna(subset=['CN_Bin', 'Coverage_Multiple'])

    sns.boxplot(data=df_plot, x='CN_Bin', y='Coverage_Multiple', ax=axes[1],
                palette='YlOrRd', showfliers=False)
    axes[1].set_title('Coverage_Multiple by SEQC2 Copy Number')
    axes[1].set_xlabel('SEQC2 Copy Number')
    axes[1].set_ylabel('ISM Coverage_Multiple')

    # Add count labels
    for i, cn_label in enumerate(cn_labels):
        n = sum(df_plot.CN_Bin == cn_label)
        axes[1].text(i, axes[1].get_ylim()[1] * 0.95, f'n={n}', ha='center', fontsize=8)

    plt.tight_layout()
    fig.savefig(FIG_DIR / '02_coverage_vs_cn.png')
    plt.close()
    print("  Saved: 02_coverage_vs_cn.png")


# ══════════════════════════════════════════════════════════════════════════════
# Step 5: FP Rate Analysis & Summary
# ══════════════════════════════════════════════════════════════════════════════

def fp_rate_summary(df, fp_rates_zone, auc_df):
    """Fig 7: FP rate bar chart + summary statistics."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 7a: FP rate by zone with confidence intervals (Wilson)
    overall_fp = sum(df.Label == 'FP') / len(df)
    zone_stats = []
    for zone in ZONE_ORDER:
        zdf = df[df.SEQC2_Zone == zone]
        n_tp = sum(zdf.Label == 'TP')
        n_fp = sum(zdf.Label == 'FP')
        n = n_tp + n_fp
        p = n_fp / n if n > 0 else 0
        # Wilson CI
        if n > 0:
            z_val = 1.96
            denom = 1 + z_val**2 / n
            center = (p + z_val**2 / (2*n)) / denom
            margin = z_val * np.sqrt((p*(1-p) + z_val**2/(4*n)) / n) / denom
            ci_lo, ci_hi = center - margin, center + margin
        else:
            ci_lo, ci_hi = 0, 0
        zone_stats.append({
            'Zone': zone, 'TP': n_tp, 'FP': n_fp, 'Total': n,
            'FP_Rate': p, 'CI_Lo': ci_lo, 'CI_Hi': ci_hi
        })

    stats_df = pd.DataFrame(zone_stats)

    bars = axes[0].bar(range(len(ZONE_ORDER)), stats_df.FP_Rate,
                       yerr=[stats_df.FP_Rate - stats_df.CI_Lo,
                             stats_df.CI_Hi - stats_df.FP_Rate],
                       color=[ZONE_COLORS[z] for z in ZONE_ORDER],
                       capsize=5, error_kw={'linewidth': 1.5})
    axes[0].axhline(y=overall_fp, color='red', linestyle='--', alpha=0.7,
                    label=f'Overall: {overall_fp:.3f}')
    axes[0].set_xticks(range(len(ZONE_ORDER)))
    axes[0].set_xticklabels(ZONE_ORDER, rotation=45, ha='right')
    axes[0].set_title('FP Rate by SEQC2 Zone (with 95% Wilson CI)')
    axes[0].set_ylabel('FP Rate')
    axes[0].legend()

    # 7b: Top AUC per zone comparison
    if auc_df is not None and not auc_df.empty:
        top_aucs = auc_df.max(axis=1)
        top_feat = auc_df.idxmax(axis=1)
        x_pos = range(len(top_aucs))
        colors_list = [ZONE_COLORS.get(z, '#999999') for z in top_aucs.index]
        axes[1].bar(x_pos, top_aucs.values, color=colors_list)
        axes[1].axhline(y=0.58, color='red', linestyle='--', alpha=0.7, label='AUC=0.58 ceiling')
        axes[1].set_xticks(list(x_pos))
        axes[1].set_xticklabels(top_aucs.index, rotation=45, ha='right')
        axes[1].set_title('Best Single-Feature AUC per Zone')
        axes[1].set_ylabel('Max AUC')
        axes[1].legend()
        for i, (auc_val, feat_name) in enumerate(zip(top_aucs.values, top_feat.values)):
            if not np.isnan(auc_val):
                axes[1].text(i, auc_val + 0.005, f'{feat_name}\n{auc_val:.3f}',
                             ha='center', va='bottom', fontsize=7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / '07_fp_rate_summary.png')
    plt.close()
    print("  Saved: 07_fp_rate_summary.png")

    return stats_df


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("SEQC2 CNV Stratification Analysis — HCC1395")
    print("=" * 70)

    # ── Step 1: Load and annotate ──
    print("\n[Step 1] Loading data...")
    seqc2_combined, gain_cn, loss_cn = load_seqc2_cnv()
    longphase_loh = load_longphase_loh()
    ism_df = load_ism_data()

    print("\n[Step 1] Annotating ISM regions with SEQC2 zones...")
    df = annotate_with_bedtools(ism_df, seqc2_combined, longphase_loh, gain_cn, loss_cn)

    # Save annotated data
    df.to_csv(DATA_DIR / 'annotated_hcc1395_cnv.tsv', sep='\t', index=False)
    print(f"\n  Saved: annotated_hcc1395_cnv.tsv ({len(df)} rows)")

    # Print zone distribution
    print("\n  Zone distribution:")
    for zone in ZONE_ORDER:
        n = sum(df.SEQC2_Zone == zone)
        n_tp = sum((df.SEQC2_Zone == zone) & (df.Label == 'TP'))
        n_fp = sum((df.SEQC2_Zone == zone) & (df.Label == 'FP'))
        fp_rate = n_fp / n if n > 0 else 0
        print(f"    {zone:12s}: {n:6d} ({n_tp:5d} TP + {n_fp:4d} FP, FP rate={fp_rate:.4f})")

    # ── Step 2: Descriptive statistics ──
    print("\n[Step 2] Plotting region distribution...")
    fp_rates_zone = plot_region_distribution(df)

    print("\n[Step 2] Plotting feature violins...")
    plot_feature_violins(df)

    print("\n[Step 2] Computing zone-specific AUCs...")
    auc_df = compute_zone_aucs(df)

    # Print key AUC results
    print("\n  Key AUC results (direction-agnostic):")
    for zone in ZONE_ORDER + ['All']:
        if zone in auc_df.index:
            row = auc_df.loc[zone]
            top_feat = row.idxmax() if not row.isna().all() else 'N/A'
            top_val = row.max() if not row.isna().all() else np.nan
            print(f"    {zone:12s}: best={top_feat} AUC={top_val:.3f}" if not np.isnan(top_val)
                  else f"    {zone:12s}: insufficient data")

    # ── Step 3: Cross-reference analysis ──
    print("\n[Step 3] SEQC2 × LongPhase cross-reference...")
    df, fp_rates_cross = cross_reference_analysis(df)

    # ── Step 4: Coverage vs CN ──
    print("\n[Step 4] Coverage_Multiple vs SEQC2 CN analysis...")
    coverage_cn_analysis(df)

    # ── Step 5: Summary ──
    print("\n[Step 5] FP rate summary...")
    stats_df = fp_rate_summary(df, fp_rates_zone, auc_df)

    # Save summary stats
    stats_df.to_csv(DATA_DIR / 'zone_fp_rate_summary.tsv', sep='\t', index=False)

    # ── Final summary ──
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"Figures: {len(list(FIG_DIR.glob('*.png')))} PNG files")
    print(f"Data: {len(list(DATA_DIR.glob('*.tsv')))} TSV files")

    # Key findings summary
    print("\n── Key Observations ──")
    print(f"  Total regions: {len(df)} ({sum(df.Label=='TP')} TP + {sum(df.Label=='FP')} FP)")
    print(f"  Overall FP rate: {sum(df.Label=='FP')/len(df):.4f}")
    for _, row in stats_df.iterrows():
        zone = row['Zone']
        print(f"  {zone:12s}: FP={row['FP_Rate']:.4f} (n={int(row['Total'])})")

    # Check if Neutral zone has better AUC
    if 'Neutral' in auc_df.index and 'All' in auc_df.index:
        neutral_best = auc_df.loc['Neutral'].max()
        all_best = auc_df.loc['All'].max()
        print(f"\n  Neutral zone best AUC: {neutral_best:.3f}")
        print(f"  All-data best AUC:     {all_best:.3f}")
        if neutral_best > all_best + 0.02:
            print("  → CNV may be a confounding factor (Neutral zone has better discrimination)")
        else:
            print("  → No evidence that CNV confounds discrimination")


if __name__ == '__main__':
    main()
