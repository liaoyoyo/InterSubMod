#!/usr/bin/env python3
"""
TO LOH Dual-Definition Report — Complete Visualization
======================================================
Two LOH definitions:
  1. LOH.bed (longphase-to): to_loh_bed_hit column
  2. HP Imbalance (ISM Potential_LOH): Potential_LOH column

Generates: 12+ publication-quality figures for the report
"""

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import matplotlib.font_manager as fm
import seaborn as sns
from sklearn.metrics import roc_auc_score
from pathlib import Path

# ── CJK Font Setup ──────────────────────────────────────────
def _get_cjk_font():
    """Find available CJK font, return FontProperties or None."""
    candidates = ['Noto Sans CJK TC', 'Noto Sans CJK JP', 'Noto Sans CJK SC',
                  'WenQuanYi Micro Hei', 'SimHei']
    available = {f.name for f in fm.fontManager.ttflist}
    for name in candidates:
        if name in available:
            return fm.FontProperties(family=name)
    return None

CJK_FONT = _get_cjk_font()

# ── Config ──────────────────────────────────────────────────
MASTER = '/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz'
OUTDIR = Path('/big7_disk/liaoyoyo2001/InterSubMod/output/to_loh_dual_definition_report')
OUTDIR.mkdir(parents=True, exist_ok=True)

SAMPLES_ORDER = ['HCC1395', 'HCC1395_DORADO', 'HCC1937', 'HCC1954', 'COLO829', 'H1437', 'H2009']

# Key features to analyze
FEATURES = {
    'HP-dependent': ['HPFineSig', 'HPFineP', 'HPFineNGroups', 'HPMergedDelta', 'HPMergedSig',
                     'HP1FamilyN', 'HP2FamilyN', 'HP_Ratio'],
    'Allele': ['AlleleDelta', 'AlleleP', 'AlleleSig'],
    'Methylation': ['PairwiseMeanDist', 'PairwiseMedianDist', 'CramersV', 'HeuristicScore', 'GlobalP'],
    'Structure': ['NumReads', 'NumCpGs', 'Coverage_Multiple', 'NHP0', 'NHP3'],
    'Caller': ['caller_af', 'caller_gq', 'caller_dp'],
    'Score': ['Quality_Score'],
}
ALL_FEATURES = [f for group in FEATURES.values() for f in group]

plt.rcParams.update({
    'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 11,
    'figure.dpi': 150, 'savefig.dpi': 150, 'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans',
})

# ── Helpers ─────────────────────────────────────────────────
def safe_auc(y_true, y_score):
    """Compute AUC safely, return NaN if not computable."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    mask = np.isfinite(y_score)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 10:
        return np.nan
    try:
        return roc_auc_score(y_true, y_score)
    except:
        return np.nan

def load_data():
    print("Loading master dataset...")
    df = pd.read_csv(MASTER, sep='\t', low_memory=False)
    to = df[df['mode'] == 'to'].copy()
    to['is_tp'] = (to['truth_label'] == 'TP').astype(int)
    # Ensure boolean columns
    for col in ['Potential_LOH', 'to_loh_bed_hit', 'HPFineSig', 'AlleleSig', 'HPMergedSig', 'PassedGating']:
        if col in to.columns:
            to[col] = to[col].astype(bool)
    # Create 4-quadrant LOH classification
    to['LOH_quadrant'] = 'Neither'
    to.loc[to['Potential_LOH'] & to['to_loh_bed_hit'], 'LOH_quadrant'] = 'Both'
    to.loc[to['Potential_LOH'] & ~to['to_loh_bed_hit'], 'LOH_quadrant'] = 'HP_only'
    to.loc[~to['Potential_LOH'] & to['to_loh_bed_hit'], 'LOH_quadrant'] = 'Bed_only'
    print(f"  TO: {len(to):,} rows, {to['is_tp'].sum():,} TP, {(1-to['is_tp']).sum():,.0f} FP")
    return to

# ── Figure 1: LOH Definition Venn / Cross-Table ────────────
def fig01_loh_definitions_overview(to):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel A: Stacked bar by quadrant
    ax = axes[0]
    quads = ['Both', 'HP_only', 'Bed_only', 'Neither']
    quad_labels = ['Both\n(LOH.bed ∩ HP)', 'HP Imbalance\nonly', 'LOH.bed\nonly', 'Neither']
    tp_counts = [to[to['LOH_quadrant']==q]['is_tp'].sum() for q in quads]
    fp_counts = [(to['LOH_quadrant']==q).sum() - to[to['LOH_quadrant']==q]['is_tp'].sum() for q in quads]

    x = np.arange(len(quads))
    bars_tp = ax.bar(x, tp_counts, color='#2ecc71', alpha=0.8, label='TP')
    bars_fp = ax.bar(x, fp_counts, bottom=tp_counts, color='#e74c3c', alpha=0.8, label='FP')
    ax.set_xticks(x)
    ax.set_xticklabels(quad_labels, fontsize=9)
    ax.set_ylabel('Count')
    ax.set_title('A. LOH Quadrant: TP/FP Distribution')
    ax.legend()
    # Add FP% labels
    for i, (tp, fp) in enumerate(zip(tp_counts, fp_counts)):
        total = tp + fp
        if total > 0:
            pct = fp / total * 100
            ax.text(i, total + total*0.01, f'FP {pct:.1f}%', ha='center', fontsize=9, fontweight='bold',
                    color='#e74c3c' if pct > 30 else '#2ecc71')

    # Panel B: Per-sample breakdown
    ax = axes[1]
    sample_data = []
    for s in SAMPLES_ORDER:
        ss = to[to['sample'] == s]
        for q in quads:
            sq = ss[ss['LOH_quadrant'] == q]
            sample_data.append({'Sample': s, 'Quadrant': q, 'N': len(sq),
                                'FP%': (1-sq['is_tp']).mean()*100 if len(sq)>0 else 0})
    sdf = pd.DataFrame(sample_data)
    pivot = sdf.pivot(index='Sample', columns='Quadrant', values='N').reindex(columns=quads)
    pivot = pivot.reindex(SAMPLES_ORDER)
    pivot.plot(kind='barh', stacked=True, ax=ax, color=['#3498db', '#9b59b6', '#e67e22', '#bdc3c7'], alpha=0.8)
    ax.set_title('B. Per-Sample LOH Quadrant Size')
    ax.set_xlabel('Region Count')
    ax.legend(title='Quadrant', fontsize=8, title_fontsize=9, loc='lower right')

    # Panel C: FP rate by quadrant per sample
    ax = axes[2]
    for q, color, marker in zip(quads, ['#3498db', '#9b59b6', '#e67e22', '#bdc3c7'], ['o', 's', '^', 'D']):
        vals = []
        for s in SAMPLES_ORDER:
            sq = to[(to['sample']==s) & (to['LOH_quadrant']==q)]
            vals.append((1-sq['is_tp']).mean()*100 if len(sq)>10 else np.nan)
        ax.plot(range(len(SAMPLES_ORDER)), vals, marker=marker, label=q, color=color, linewidth=2, markersize=8)
    ax.set_xticks(range(len(SAMPLES_ORDER)))
    ax.set_xticklabels(SAMPLES_ORDER, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('FP Rate (%)')
    ax.set_title('C. FP Rate by LOH Quadrant × Sample')
    ax.legend(fontsize=8)
    ax.axhline(y=30, ls='--', color='gray', alpha=0.5, label='30% ref')
    ax.grid(axis='y', alpha=0.3)

    fig.suptitle('Figure 1: TO LOH Dual-Definition Overview', fontsize=15, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig01_loh_definitions_overview.png')
    plt.close()
    print("  fig01 done")

# ── Figure 2: AUC Heatmap — Features × LOH Definition ──────
def fig02_auc_heatmap_dual_loh(to):
    # Compute AUC for each feature × LOH definition × LOH/Non-LOH
    loh_defs = {
        'LOH.bed=T': to['to_loh_bed_hit'] == True,
        'LOH.bed=F': to['to_loh_bed_hit'] == False,
        'HP_Imb=T': to['Potential_LOH'] == True,
        'HP_Imb=F': to['Potential_LOH'] == False,
        'Both=T': to['LOH_quadrant'] == 'Both',
        'HP_only': to['LOH_quadrant'] == 'HP_only',
        'Neither': to['LOH_quadrant'] == 'Neither',
    }

    feature_list = ['AlleleDelta', 'HPFineP', 'HPFineSig', 'HPFineNGroups',
                    'PairwiseMeanDist', 'CramersV', 'NumReads', 'Quality_Score',
                    'caller_af', 'caller_gq', 'Coverage_Multiple', 'HPMergedDelta',
                    'HP1FamilyN', 'HP2FamilyN']

    results = {}
    for lname, lmask in loh_defs.items():
        subset = to[lmask]
        for feat in feature_list:
            if feat in subset.columns:
                vals = pd.to_numeric(subset[feat], errors='coerce')
                # For boolean features, convert
                if feat in ['HPFineSig', 'AlleleSig', 'HPMergedSig']:
                    vals = vals.astype(float)
                auc = safe_auc(subset['is_tp'], vals)
                results[(feat, lname)] = auc

    # Build matrix
    mat = pd.DataFrame(index=feature_list, columns=list(loh_defs.keys()), dtype=float)
    for (feat, lname), auc in results.items():
        mat.loc[feat, lname] = auc

    fig, ax = plt.subplots(figsize=(14, 10))
    # Custom colormap: red < 0.48, white ~ 0.50, green > 0.52
    cmap = sns.diverging_palette(10, 130, s=80, l=55, as_cmap=True)
    sns.heatmap(mat.astype(float), annot=True, fmt='.3f', cmap=cmap, center=0.50,
                vmin=0.35, vmax=0.70, ax=ax, linewidths=0.5,
                cbar_kws={'label': 'AUC (TP-favored > 0.5)'})
    ax.set_title('Figure 2: Feature AUC by LOH Definition (TO mode, 7 samples pooled)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Feature')
    ax.set_xlabel('LOH Definition Subset')

    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig02_auc_heatmap_dual_loh.png')
    plt.close()
    print("  fig02 done")
    return mat

# ── Figure 3: Per-Sample AUC — AlleleDelta × Two Definitions ──
def fig03_per_sample_allele_auc(to):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for idx, (loh_col, title) in enumerate([
        ('to_loh_bed_hit', 'LOH.bed (longphase-to)'),
        ('Potential_LOH', 'HP Imbalance (ISM)')
    ]):
        ax = axes[idx]
        data_loh = []
        data_nonloh = []
        for s in SAMPLES_ORDER:
            for loh_val, label, data_list in [(True, 'LOH', data_loh), (False, 'Non-LOH', data_nonloh)]:
                subset = to[(to['sample']==s) & (to[loh_col]==loh_val)]
                auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['AlleleDelta'], errors='coerce'))
                n = len(subset)
                data_list.append({'Sample': s, 'AUC': auc, 'N': n})

        x = np.arange(len(SAMPLES_ORDER))
        width = 0.35
        loh_aucs = [d['AUC'] for d in data_loh]
        nonloh_aucs = [d['AUC'] for d in data_nonloh]

        bars1 = ax.bar(x - width/2, loh_aucs, width, label='LOH', color='#e74c3c', alpha=0.7)
        bars2 = ax.bar(x + width/2, nonloh_aucs, width, label='Non-LOH', color='#3498db', alpha=0.7)

        ax.axhline(y=0.50, ls='--', color='gray', alpha=0.5)
        ax.axhline(y=0.58, ls='--', color='red', alpha=0.3, label='0.58 threshold')
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLES_ORDER, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('AlleleDelta AUC')
        ax.set_title(f'AlleleDelta AUC by {title}')
        ax.legend(fontsize=8)
        ax.set_ylim(0.35, 0.75)
        ax.grid(axis='y', alpha=0.3)

        # Add sample sizes
        for i, (ld, nld) in enumerate(zip(data_loh, data_nonloh)):
            ax.text(i - width/2, 0.36, f'n={ld["N"]//1000}K', ha='center', fontsize=6, color='#e74c3c')
            ax.text(i + width/2, 0.36, f'n={nld["N"]//1000}K', ha='center', fontsize=6, color='#3498db')

    fig.suptitle('Figure 3: AlleleDelta AUC — LOH.bed vs HP Imbalance × Per-Sample', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig03_per_sample_allele_auc.png')
    plt.close()
    print("  fig03 done")

# ── Figure 4: HPFineP AUC — dual definition per sample ─────
def fig04_per_sample_hpfinep_auc(to):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for idx, (loh_col, title) in enumerate([
        ('to_loh_bed_hit', 'LOH.bed (longphase-to)'),
        ('Potential_LOH', 'HP Imbalance (ISM)')
    ]):
        ax = axes[idx]
        for loh_val, label, color, marker in [
            (True, 'LOH', '#e74c3c', 'o'),
            (False, 'Non-LOH', '#3498db', 's')
        ]:
            aucs = []
            for s in SAMPLES_ORDER:
                subset = to[(to['sample']==s) & (to[loh_col]==loh_val)]
                vals = pd.to_numeric(subset['HPFineP'], errors='coerce')
                auc = safe_auc(subset['is_tp'], -vals)  # lower P = more significant = TP
                aucs.append(auc)
            ax.plot(range(len(SAMPLES_ORDER)), aucs, marker=marker, label=label, color=color, linewidth=2, markersize=8)

        ax.axhline(y=0.50, ls='--', color='gray', alpha=0.5)
        ax.axhline(y=0.58, ls='--', color='red', alpha=0.3)
        ax.set_xticks(range(len(SAMPLES_ORDER)))
        ax.set_xticklabels(SAMPLES_ORDER, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('-HPFineP AUC')
        ax.set_title(f'-HPFineP AUC by {title}')
        ax.legend(fontsize=9)
        ax.set_ylim(0.40, 0.75)
        ax.grid(axis='y', alpha=0.3)

    fig.suptitle('Figure 4: HPFineP AUC — LOH.bed vs HP Imbalance × Per-Sample', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig04_per_sample_hpfinep_auc.png')
    plt.close()
    print("  fig04 done")

# ── Figure 5: Feature Top-10 AUC Comparison ────────────────
def fig05_top10_auc_comparison(to):
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    conditions = [
        ('to_loh_bed_hit', True, 'LOH.bed = True'),
        ('to_loh_bed_hit', False, 'LOH.bed = False (Non-LOH)'),
        ('Potential_LOH', True, 'HP Imbalance = True'),
        ('Potential_LOH', False, 'HP Imbalance = False (Non-LOH)'),
    ]

    feature_list = ['AlleleDelta', 'HPFineP', 'HPFineNGroups', 'HPMergedDelta',
                    'PairwiseMeanDist', 'PairwiseMedianDist', 'CramersV', 'NumReads',
                    'Quality_Score', 'caller_af', 'caller_gq', 'Coverage_Multiple',
                    'HP1FamilyN', 'HP2FamilyN', 'HeuristicScore', 'GlobalP', 'NHP0', 'NHP3']

    for ax, (col, val, title) in zip(axes.flat, conditions):
        subset = to[to[col] == val]
        n_tp = subset['is_tp'].sum()
        n_fp = len(subset) - n_tp

        aucs = {}
        for feat in feature_list:
            if feat in subset.columns:
                vals = pd.to_numeric(subset[feat], errors='coerce')
                if feat == 'HPFineP':
                    vals = -vals
                auc = safe_auc(subset['is_tp'], vals)
                if not np.isnan(auc):
                    aucs[feat] = auc

        # Sort by |AUC - 0.5| (distance from random)
        sorted_feats = sorted(aucs.keys(), key=lambda f: abs(aucs[f]-0.5), reverse=True)[:12]
        sorted_aucs = [aucs[f] for f in sorted_feats]

        colors = ['#2ecc71' if a > 0.52 else '#e74c3c' if a < 0.48 else '#bdc3c7' for a in sorted_aucs]

        y = range(len(sorted_feats))
        ax.barh(y, [a - 0.5 for a in sorted_aucs], left=0.5, color=colors, alpha=0.8, height=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(sorted_feats, fontsize=8)
        ax.axvline(x=0.5, ls='-', color='gray', alpha=0.5)
        ax.axvline(x=0.58, ls='--', color='red', alpha=0.3)
        ax.axvline(x=0.42, ls='--', color='red', alpha=0.3)
        ax.set_xlabel('AUC')
        ax.set_title(f'{title}\n(TP={n_tp:,}, FP={n_fp:,})', fontsize=11)
        ax.set_xlim(0.35, 0.70)
        ax.grid(axis='x', alpha=0.3)

        # Annotate AUC values
        for i, (f, a) in enumerate(zip(sorted_feats, sorted_aucs)):
            ax.text(max(a, 0.5) + 0.005, i, f'{a:.3f}', va='center', fontsize=8, fontweight='bold')

    fig.suptitle('Figure 5: Top-12 Feature AUC — LOH.bed vs HP Imbalance (TO, 7 samples)', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig05_top10_auc_comparison.png')
    plt.close()
    print("  fig05 done")

# ── Figure 6: HP_only quadrant deep dive ───────────────────
def fig06_hp_only_quadrant_analysis(to):
    """The 63,610 regions that are HP Imbalance but NOT LOH.bed — what are they?"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    hp_only = to[to['LOH_quadrant'] == 'HP_only']
    both = to[to['LOH_quadrant'] == 'Both']
    neither = to[to['LOH_quadrant'] == 'Neither']

    # Panel A: HP_Ratio distribution by quadrant
    ax = axes[0, 0]
    for subset, label, color in [(both, 'Both', '#3498db'), (hp_only, 'HP_only', '#9b59b6'), (neither, 'Neither', '#bdc3c7')]:
        vals = pd.to_numeric(subset['HP_Ratio'], errors='coerce').dropna()
        if len(vals) > 100:
            ax.hist(vals, bins=50, alpha=0.5, label=f'{label} (n={len(vals):,})', color=color, density=True)
    ax.set_xlabel('HP_Ratio')
    ax.set_ylabel('Density')
    ax.set_title('A. HP_Ratio Distribution by Quadrant')
    ax.legend(fontsize=9)

    # Panel B: Coverage_Multiple by quadrant
    ax = axes[0, 1]
    for subset, label, color in [(both, 'Both', '#3498db'), (hp_only, 'HP_only', '#9b59b6'), (neither, 'Neither', '#bdc3c7')]:
        vals = pd.to_numeric(subset['Coverage_Multiple'], errors='coerce').dropna()
        vals = vals[vals < 3]  # clip
        if len(vals) > 100:
            ax.hist(vals, bins=50, alpha=0.5, label=f'{label} (n={len(vals):,})', color=color, density=True)
    ax.set_xlabel('Coverage_Multiple')
    ax.set_ylabel('Density')
    ax.set_title('B. Coverage Distribution by Quadrant')
    ax.legend(fontsize=9)
    ax.axvline(x=0.8, ls='--', color='red', alpha=0.3, label='deletion threshold')

    # Panel C: Per-sample HP_only size
    ax = axes[1, 0]
    sample_data = []
    for s in SAMPLES_ORDER:
        ss = to[to['sample'] == s]
        for q in ['Both', 'HP_only', 'Neither']:
            sq = ss[ss['LOH_quadrant'] == q]
            sample_data.append({'Sample': s, 'Quadrant': q, 'Fraction': len(sq) / len(ss) * 100})
    sdf = pd.DataFrame(sample_data)
    pivot = sdf.pivot(index='Sample', columns='Quadrant', values='Fraction').reindex(columns=['Both', 'HP_only', 'Neither'])
    pivot.plot(kind='bar', stacked=True, ax=ax, color=['#3498db', '#9b59b6', '#bdc3c7'], alpha=0.8)
    ax.set_ylabel('Fraction (%)')
    ax.set_title('C. LOH Quadrant Fraction by Sample')
    ax.legend(fontsize=8)
    ax.set_xticklabels(SAMPLES_ORDER, rotation=45, ha='right', fontsize=8)

    # Panel D: Key feature AUC in HP_only vs Both
    ax = axes[1, 1]
    features = ['AlleleDelta', 'HPFineNGroups', 'PairwiseMeanDist', 'CramersV', 'NumReads', 'Quality_Score', 'caller_af']
    both_aucs = []
    hp_only_aucs = []
    for feat in features:
        vals_b = pd.to_numeric(both[feat], errors='coerce')
        vals_h = pd.to_numeric(hp_only[feat], errors='coerce')
        auc_b = safe_auc(both['is_tp'], vals_b)
        auc_h = safe_auc(hp_only['is_tp'], vals_h)
        both_aucs.append(auc_b)
        hp_only_aucs.append(auc_h)

    x = np.arange(len(features))
    width = 0.35
    ax.bar(x - width/2, both_aucs, width, label='Both (LOH.bed ∩ HP)', color='#3498db', alpha=0.7)
    ax.bar(x + width/2, hp_only_aucs, width, label='HP_only', color='#9b59b6', alpha=0.7)
    ax.axhline(y=0.5, ls='--', color='gray', alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('AUC')
    ax.set_title('D. Feature AUC: Both vs HP_only')
    ax.legend(fontsize=9)
    ax.set_ylim(0.35, 0.70)
    ax.grid(axis='y', alpha=0.3)

    fig.suptitle('Figure 6: HP_only Quadrant Deep Dive (HP Imbalance but NOT LOH.bed)', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig06_hp_only_quadrant_analysis.png')
    plt.close()
    print("  fig06 done")

# ── Figure 7: Composition Effect for Both Definitions ──────
def fig07_composition_effect(to):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for idx, (loh_col, title) in enumerate([
        ('to_loh_bed_hit', 'LOH.bed'),
        ('Potential_LOH', 'HP Imbalance')
    ]):
        ax = axes[idx]
        features = ['AlleleDelta', 'HPFineNGroups', 'PairwiseMeanDist', 'CramersV', 'NumReads']

        loh = to[to[loh_col] == True]
        nonloh = to[to[loh_col] == False]

        loh_aucs = []
        nonloh_aucs = []
        all_aucs = []

        for feat in features:
            vals_l = pd.to_numeric(loh[feat], errors='coerce')
            vals_n = pd.to_numeric(nonloh[feat], errors='coerce')
            vals_a = pd.to_numeric(to[feat], errors='coerce')
            loh_aucs.append(safe_auc(loh['is_tp'], vals_l))
            nonloh_aucs.append(safe_auc(nonloh['is_tp'], vals_n))
            all_aucs.append(safe_auc(to['is_tp'], vals_a))

        x = np.arange(len(features))
        width = 0.25
        ax.bar(x - width, loh_aucs, width, label='LOH', color='#e74c3c', alpha=0.7)
        ax.bar(x, nonloh_aucs, width, label='Non-LOH', color='#3498db', alpha=0.7)
        ax.bar(x + width, all_aucs, width, label='All', color='#95a5a6', alpha=0.7)
        ax.axhline(y=0.5, ls='--', color='gray', alpha=0.5)
        ax.axhline(y=0.58, ls='--', color='red', alpha=0.3)
        ax.set_xticks(x)
        ax.set_xticklabels(features, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('AUC')
        ax.set_title(f'{title} Definition')
        ax.legend(fontsize=9)
        ax.set_ylim(0.40, 0.70)
        ax.grid(axis='y', alpha=0.3)

    fig.suptitle('Figure 7: LOH vs Non-LOH Feature AUC (Two Definitions)', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig07_composition_effect.png')
    plt.close()
    print("  fig07 done")

# ── Figure 8: HCC1395 Detailed Breakdown ──────────────────
def fig08_hcc1395_detailed(to):
    hcc = to[to['sample'] == 'HCC1395']

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    # Panel A: LOH quadrant TP/FP for HCC1395
    ax = axes[0, 0]
    quads = ['Both', 'HP_only', 'Bed_only', 'Neither']
    tp_counts = [hcc[hcc['LOH_quadrant']==q]['is_tp'].sum() for q in quads]
    fp_counts = [(hcc['LOH_quadrant']==q).sum() - hcc[hcc['LOH_quadrant']==q]['is_tp'].sum() for q in quads]
    x = np.arange(len(quads))
    ax.bar(x, tp_counts, color='#2ecc71', alpha=0.8, label='TP')
    ax.bar(x, fp_counts, bottom=tp_counts, color='#e74c3c', alpha=0.8, label='FP')
    ax.set_xticks(x)
    ax.set_xticklabels(['Both\n(LOH.bed∩HP)', 'HP\nonly', 'Bed\nonly', 'Neither'], fontsize=9)
    ax.set_title('A. HCC1395 LOH Quadrant')
    ax.legend()
    for i, (tp, fp) in enumerate(zip(tp_counts, fp_counts)):
        total = tp + fp
        if total > 0:
            ax.text(i, total+100, f'FP{fp/total*100:.0f}%', ha='center', fontsize=9, fontweight='bold')

    # Panel B: AlleleDelta distribution in LOH.bed
    ax = axes[0, 1]
    for lbl, mask, color in [('TP', hcc['is_tp']==1, '#2ecc71'), ('FP', hcc['is_tp']==0, '#e74c3c')]:
        subset = hcc[mask & (hcc['to_loh_bed_hit']==True)]
        vals = pd.to_numeric(subset['AlleleDelta'], errors='coerce').dropna()
        if len(vals) > 10:
            ax.hist(vals, bins=50, alpha=0.5, label=f'{lbl} LOH.bed (n={len(vals):,})', color=color, density=True)
    ax.set_xlabel('AlleleDelta')
    ax.set_title('B. AlleleDelta in LOH.bed (HCC1395)')
    ax.legend(fontsize=8)

    # Panel C: AlleleDelta in HP Imbalance
    ax = axes[0, 2]
    for lbl, mask, color in [('TP', hcc['is_tp']==1, '#2ecc71'), ('FP', hcc['is_tp']==0, '#e74c3c')]:
        subset = hcc[mask & (hcc['Potential_LOH']==True)]
        vals = pd.to_numeric(subset['AlleleDelta'], errors='coerce').dropna()
        if len(vals) > 10:
            ax.hist(vals, bins=50, alpha=0.5, label=f'{lbl} HP_Imb (n={len(vals):,})', color=color, density=True)
    ax.set_xlabel('AlleleDelta')
    ax.set_title('C. AlleleDelta in HP Imbalance (HCC1395)')
    ax.legend(fontsize=8)

    # Panel D: HPFineP distribution in LOH.bed
    ax = axes[1, 0]
    for lbl, mask, color in [('TP', hcc['is_tp']==1, '#2ecc71'), ('FP', hcc['is_tp']==0, '#e74c3c')]:
        subset = hcc[mask & (hcc['to_loh_bed_hit']==True)]
        vals = pd.to_numeric(subset['HPFineP'], errors='coerce').dropna()
        vals = vals[vals > 0]
        if len(vals) > 10:
            ax.hist(np.log10(vals), bins=50, alpha=0.5, label=f'{lbl} (n={len(vals):,})', color=color, density=True)
    ax.set_xlabel('log10(HPFineP)')
    ax.set_title('D. HPFineP in LOH.bed (HCC1395)')
    ax.legend(fontsize=8)
    ax.axvline(x=np.log10(0.05), ls='--', color='black', alpha=0.5, label='p=0.05')

    # Panel E: Feature AUC comparison for HCC1395 — two definitions
    ax = axes[1, 1]
    features = ['AlleleDelta', 'HPFineNGroups', 'PairwiseMeanDist', 'CramersV', 'Quality_Score', 'caller_af']
    bed_aucs, hp_aucs = [], []
    for feat in features:
        for loh_col, aucs_list in [('to_loh_bed_hit', bed_aucs), ('Potential_LOH', hp_aucs)]:
            subset = hcc[hcc[loh_col] == True]
            vals = pd.to_numeric(subset[feat], errors='coerce')
            aucs_list.append(safe_auc(subset['is_tp'], vals))

    x = np.arange(len(features))
    width = 0.35
    ax.bar(x - width/2, bed_aucs, width, label='LOH.bed', color='#e74c3c', alpha=0.7)
    ax.bar(x + width/2, hp_aucs, width, label='HP Imbalance', color='#9b59b6', alpha=0.7)
    ax.axhline(y=0.5, ls='--', color='gray', alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('AUC')
    ax.set_title('E. Feature AUC in LOH: LOH.bed vs HP Imb (HCC1395)')
    ax.legend(fontsize=9)
    ax.set_ylim(0.35, 0.70)
    ax.grid(axis='y', alpha=0.3)

    # Panel F: Feature AUC in Non-LOH
    ax = axes[1, 2]
    bed_nonloh, hp_nonloh = [], []
    for feat in features:
        for loh_col, aucs_list in [('to_loh_bed_hit', bed_nonloh), ('Potential_LOH', hp_nonloh)]:
            subset = hcc[hcc[loh_col] == False]
            vals = pd.to_numeric(subset[feat], errors='coerce')
            aucs_list.append(safe_auc(subset['is_tp'], vals))

    x = np.arange(len(features))
    ax.bar(x - width/2, bed_nonloh, width, label='LOH.bed=F', color='#e74c3c', alpha=0.7)
    ax.bar(x + width/2, hp_nonloh, width, label='HP_Imb=F', color='#9b59b6', alpha=0.7)
    ax.axhline(y=0.5, ls='--', color='gray', alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('AUC')
    ax.set_title('F. Feature AUC in Non-LOH (HCC1395)')
    ax.legend(fontsize=9)
    ax.set_ylim(0.35, 0.70)
    ax.grid(axis='y', alpha=0.3)

    fig.suptitle('Figure 8: HCC1395 Detailed — LOH.bed vs HP Imbalance', fontsize=15, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig08_hcc1395_detailed.png')
    plt.close()
    print("  fig08 done")

# ── Figure 9: Multi-Feature Combo AUC by Definition ───────
def fig09_combo_auc(to):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    combo_features = ['AlleleDelta', 'HPFineNGroups', 'PairwiseMeanDist', 'CramersV', 'NumReads']

    for idx, (loh_col, title) in enumerate([
        ('to_loh_bed_hit', 'LOH.bed'),
        ('Potential_LOH', 'HP Imbalance')
    ]):
        ax = axes[idx]
        results = {}

        for loh_val, loh_label in [(True, 'LOH'), (False, 'Non-LOH')]:
            subset = to[to[loh_col] == loh_val].copy()

            # Single features
            for feat in combo_features:
                vals = pd.to_numeric(subset[feat], errors='coerce')
                auc = safe_auc(subset['is_tp'], vals)
                results[(f'Single: {feat}', loh_label)] = auc

            # Combo: all 5 features LR
            try:
                X = subset[combo_features].apply(pd.to_numeric, errors='coerce').fillna(0)
                y = subset['is_tp'].values
                mask = np.isfinite(X.values).all(axis=1)
                X, y = X.values[mask], y[mask]
                if len(np.unique(y)) == 2 and len(y) > 100:
                    scaler = StandardScaler()
                    X_s = scaler.fit_transform(X)
                    lr = LogisticRegression(max_iter=1000, C=1.0)
                    lr.fit(X_s, y)
                    pred = lr.predict_proba(X_s)[:, 1]
                    results[('Combo: LR(5feat)', loh_label)] = roc_auc_score(y, pred)
            except:
                pass

        # Plot
        combo_keys = sorted(set(k[0] for k in results.keys()))
        loh_vals = [results.get((k, 'LOH'), np.nan) for k in combo_keys]
        nonloh_vals = [results.get((k, 'Non-LOH'), np.nan) for k in combo_keys]

        x = np.arange(len(combo_keys))
        width = 0.35
        ax.barh(x - width/2, loh_vals, width, label='LOH', color='#e74c3c', alpha=0.7)
        ax.barh(x + width/2, nonloh_vals, width, label='Non-LOH', color='#3498db', alpha=0.7)
        ax.axvline(x=0.5, ls='--', color='gray', alpha=0.5)
        ax.axvline(x=0.58, ls='--', color='red', alpha=0.3)
        ax.set_yticks(x)
        ax.set_yticklabels([k.replace('Single: ', '') for k in combo_keys], fontsize=9)
        ax.set_xlabel('AUC')
        ax.set_title(f'{title} Definition')
        ax.legend(fontsize=9)
        ax.set_xlim(0.40, 0.70)
        ax.grid(axis='x', alpha=0.3)

    fig.suptitle('Figure 9: Single vs Combo Feature AUC by LOH Definition', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig09_combo_auc.png')
    plt.close()
    print("  fig09 done")

# ── Figure 10: Quality Score by LOH Definition ─────────────
def fig10_quality_score_analysis(to):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel A: QS distribution by LOH.bed
    ax = axes[0, 0]
    for loh_val, label, ls in [(True, 'LOH.bed=T', '-'), (False, 'LOH.bed=F', '--')]:
        for tp_val, tp_label, color in [(1, 'TP', '#2ecc71'), (0, 'FP', '#e74c3c')]:
            subset = to[(to['to_loh_bed_hit']==loh_val) & (to['is_tp']==tp_val)]
            vals = pd.to_numeric(subset['Quality_Score'], errors='coerce').dropna()
            if len(vals) > 100:
                ax.hist(vals, bins=50, alpha=0.3, label=f'{label} {tp_label} (n={len(vals):,})',
                        color=color, linestyle=ls, density=True, histtype='step', linewidth=2)
    ax.set_xlabel('Quality_Score')
    ax.set_title('A. QS Distribution by LOH.bed')
    ax.legend(fontsize=7)

    # Panel B: QS AUC per sample × LOH definition
    ax = axes[0, 1]
    for loh_col, label, color, marker in [
        ('to_loh_bed_hit', 'LOH.bed=T', '#e74c3c', 'o'),
        ('to_loh_bed_hit', 'LOH.bed=F', '#3498db', 's'),
        ('Potential_LOH', 'HP_Imb=T', '#e74c3c', '^'),
        ('Potential_LOH', 'HP_Imb=F', '#9b59b6', 'D'),
    ]:
        loh_val = True if '=T' in label else False
        aucs = []
        for s in SAMPLES_ORDER:
            subset = to[(to['sample']==s) & (to[loh_col]==loh_val)]
            vals = pd.to_numeric(subset['Quality_Score'], errors='coerce')
            aucs.append(safe_auc(subset['is_tp'], vals))
        ax.plot(range(len(SAMPLES_ORDER)), aucs, marker=marker, label=label,
                color=color, linewidth=1.5, markersize=7,
                linestyle='-' if 'bed' in label.lower() else '--')
    ax.axhline(y=0.5, ls='--', color='gray', alpha=0.5)
    ax.set_xticks(range(len(SAMPLES_ORDER)))
    ax.set_xticklabels(SAMPLES_ORDER, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Quality_Score AUC')
    ax.set_title('B. QS AUC by LOH Definition × Sample')
    ax.legend(fontsize=7)
    ax.set_ylim(0.35, 0.65)
    ax.grid(axis='y', alpha=0.3)

    # Panel C: caller_af AUC per LOH definition
    ax = axes[1, 0]
    for loh_col, colors in [('to_loh_bed_hit', ['#e74c3c', '#3498db']), ('Potential_LOH', ['#9b59b6', '#f39c12'])]:
        for loh_val, label, color in [(True, f'{loh_col.split("_")[0]}=T', colors[0]),
                                       (False, f'{loh_col.split("_")[0]}=F', colors[1])]:
            aucs = []
            for s in SAMPLES_ORDER:
                subset = to[(to['sample']==s) & (to[loh_col]==loh_val)]
                vals = pd.to_numeric(subset['caller_af'], errors='coerce')
                aucs.append(safe_auc(subset['is_tp'], vals))
            style = '-' if 'bed' in loh_col else '--'
            ax.plot(range(len(SAMPLES_ORDER)), aucs, marker='o', label=f'{loh_col}: {label}',
                    color=color, linewidth=1.5, linestyle=style)
    ax.axhline(y=0.5, ls='--', color='gray', alpha=0.5)
    ax.set_xticks(range(len(SAMPLES_ORDER)))
    ax.set_xticklabels(SAMPLES_ORDER, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('caller_af AUC')
    ax.set_title('C. caller_af AUC (note: < 0.5 = FP-favored)')
    ax.legend(fontsize=7)
    ax.grid(axis='y', alpha=0.3)

    # Panel D: Summary table
    ax = axes[1, 1]
    ax.axis('off')
    summary_data = []
    for loh_col, def_name in [('to_loh_bed_hit', 'LOH.bed'), ('Potential_LOH', 'HP Imbalance')]:
        for loh_val in [True, False]:
            subset = to[to[loh_col] == loh_val]
            n = len(subset)
            fp_rate = (1 - subset['is_tp']).mean() * 100
            qs_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['Quality_Score'], errors='coerce'))
            ad_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['AlleleDelta'], errors='coerce'))
            af_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['caller_af'], errors='coerce'))
            label = f'{def_name}={"T" if loh_val else "F"}'
            summary_data.append([label, f'{n:,}', f'{fp_rate:.1f}%', f'{qs_auc:.3f}', f'{ad_auc:.3f}', f'{af_auc:.3f}'])

    table = ax.table(cellText=summary_data,
                     colLabels=['Subset', 'N', 'FP%', 'QS AUC', 'AD AUC', 'AF AUC'],
                     loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)
    ax.set_title('D. Summary: Key AUCs by LOH Definition', fontsize=12, fontweight='bold', pad=20)

    fig.suptitle('Figure 10: Quality Score & Caller Features by LOH Definition', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig10_quality_score_analysis.png')
    plt.close()
    print("  fig10 done")

# ── Figure 11: Complete Evidence Chain Visual ──────────────
def fig11_evidence_chain(to):
    fig, ax = plt.subplots(figsize=(20, 14))
    ax.axis('off')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)

    # CJK font kwargs for text with Chinese characters
    cjk_kw = {'fontproperties': CJK_FONT} if CJK_FONT else {}

    # Title
    ax.text(5, 13.5, 'TO ISM LOH 區分力 — 完整推論鏈', fontsize=18, fontweight='bold',
            ha='center', va='center', **cjk_kw)

    # Evidence blocks
    blocks = [
        (1, 12, 4, 1.2, '#3498db', 'Step 1: O1-O10 系統觀察\n82 圖表 × 9 主題\n→ TO 無單一 AUC > 0.58\n→ QS AUC=0.507 隨機'),
        (5.5, 12, 4, 1.2, '#9b59b6', 'Step 2: O11-O13 甲基化假說\nO11 heterogeneity: NEGATIVE\nO12 LOH scenario: NEGATIVE\nO13 cross-region: NEGATIVE'),
        (1, 10.2, 4, 1.2, '#e67e22', 'Step 3: G1-G7 VCF 特徵\n60+ features AUC < 0.64\nTP loss ≤2% → FP removal=0%\n根因: 殘餘FP是罕見germline'),
        (5.5, 10.2, 4, 1.2, '#e74c3c', 'Step 4: LOH 深入 (Q1-Q6 + W1-3)\nHP2FamilyN 0.72→0.54 (循環)\nAlleleDelta AUC=0.556 (弱)\nVoting combo 0.577 < 0.58'),
        (1, 8.4, 4, 1.2, '#1abc9c', 'Step 5: Self-Phasing 因果鏈\n62% LOH 消失 (d=-1.20)\nSomatic bias 17.3:1→消除\nHP features 被系統性汙染'),
        (5.5, 8.4, 4, 1.2, '#f39c12', 'Step 6: R1-R5 + Option C\nCramersV 93% 零=2×2缺陷\nHP-free AUC=0.512 (隨機)\n所有信號來自HP tags'),
        (3.25, 6.6, 4, 1.2, '#c0392b', 'Step 7: TO-pure LOSO + O9 FN\ncaller_af=0.654 > 全ISM\nISM增量+0.003~0.030\nFN≡TP in methylation'),
    ]

    for x, y, w, h, color, text in blocks:
        rect = FancyBboxPatch((x-w/2, y-h/2), w, h,
                              boxstyle="round,pad=0.1", facecolor=color, alpha=0.15, edgecolor=color, linewidth=2)
        ax.add_patch(rect)
        ax.text(x, y, text, fontsize=8, ha='center', va='center', **cjk_kw)

    # Arrows
    for (x1, y1), (x2, y2) in [
        ((3, 11.3), (3, 10.8)), ((7.5, 11.3), (7.5, 10.8)),
        ((3, 9.5), (3, 9.0)), ((7.5, 9.5), (7.5, 9.0)),
        ((3, 7.7), (4.25, 7.2)), ((7.5, 7.7), (6.25, 7.2)),
    ]:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))

    # Final conclusion box
    conclusion_text = (
        '╔═══════════════════════════════════════════════════════════════╗\n'
        '║  結論: TO ISM 在 LOH 內外均無法區分 TP/FP                   ║\n'
        '║  • LOH.bed 區域: max AUC=0.556 (AlleleDelta, confound-free) ║\n'
        '║  • HP Imb 區域:  max AUC=0.556 (同上)                       ║\n'
        '║  • Non-LOH:      max AUC=0.643 (read count proxy)           ║\n'
        '║  • 根因1: Self-phasing 汙染 HP features                     ║\n'
        '║  • 根因2: Germline ASM ≈ Somatic heterogeneity              ║\n'
        '║  • 唯一有效: caller_af (AUC=0.654)                          ║\n'
        '╚═══════════════════════════════════════════════════════════════╝'
    )
    ax.text(5, 3.5, conclusion_text, fontsize=9, ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='#fadbd8', edgecolor='#c0392b', linewidth=2),
            **cjk_kw)

    # Positive findings
    positive_text = (
        '正面發現:\n'
        '✓ HPFineNGroups: somatic heterogeneity marker (N≥4 TP 89.1%)\n'
        '✓ ASM 廣泛存在 (32-66%), ISM PERMANOVA 可偵測 entropy imbalance\n'
        '✓ AlleleDelta: confound-free true signal (7/7 consistent)\n'
        '✓ Phase 1A paired-pure: ΔF1=+0.0112'
    )
    ax.text(5, 1.2, positive_text, fontsize=9, ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='#d5f5e3', edgecolor='#27ae60', linewidth=1.5),
            **cjk_kw)

    fig.savefig(OUTDIR / 'fig11_evidence_chain.png')
    plt.close()
    print("  fig11 done")

# ── Figure 12: Per-Sample × Dual Definition Summary Table ──
def fig12_summary_table(to):
    fig, ax = plt.subplots(figsize=(22, 10))
    ax.axis('off')

    # Build comprehensive table
    rows = []
    for s in SAMPLES_ORDER:
        ss = to[to['sample'] == s]
        for loh_col, def_name in [('to_loh_bed_hit', 'LOH.bed'), ('Potential_LOH', 'HP_Imb')]:
            for loh_val, loh_label in [(True, 'LOH'), (False, 'Non-LOH')]:
                subset = ss[(ss[loh_col] == loh_val)]
                n = len(subset)
                if n < 20:
                    continue
                n_tp = subset['is_tp'].sum()
                n_fp = n - n_tp
                fp_pct = n_fp / n * 100

                ad_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['AlleleDelta'], errors='coerce'))
                hp_auc = safe_auc(subset['is_tp'], -pd.to_numeric(subset['HPFineP'], errors='coerce'))
                ng_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['HPFineNGroups'], errors='coerce'))
                pw_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['PairwiseMeanDist'], errors='coerce'))
                qs_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['Quality_Score'], errors='coerce'))
                af_auc = safe_auc(subset['is_tp'], pd.to_numeric(subset['caller_af'], errors='coerce'))

                rows.append([
                    s, f'{def_name}:{loh_label}', f'{n:,}', f'{fp_pct:.1f}%',
                    f'{ad_auc:.3f}' if not np.isnan(ad_auc) else '—',
                    f'{hp_auc:.3f}' if not np.isnan(hp_auc) else '—',
                    f'{ng_auc:.3f}' if not np.isnan(ng_auc) else '—',
                    f'{pw_auc:.3f}' if not np.isnan(pw_auc) else '—',
                    f'{qs_auc:.3f}' if not np.isnan(qs_auc) else '—',
                    f'{af_auc:.3f}' if not np.isnan(af_auc) else '—',
                ])

    cols = ['Sample', 'Definition', 'N', 'FP%', 'AlleleDelta', '-HPFineP', 'NGroups', 'PwMeanDist', 'QS', 'caller_af']

    table = ax.table(cellText=rows, colLabels=cols, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.3)

    # Color cells based on AUC value
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor('#2c3e50')
            cell.set_text_props(color='white', fontweight='bold')
        elif col >= 4:  # AUC columns
            try:
                val = float(cell.get_text().get_text())
                if val > 0.58:
                    cell.set_facecolor('#d5f5e3')
                elif val < 0.48:
                    cell.set_facecolor('#fadbd8')
                elif val > 0.52:
                    cell.set_facecolor('#fdebd0')
            except:
                pass

    ax.set_title('Figure 12: Complete Per-Sample × Dual-LOH-Definition AUC Table (TO mode)',
                 fontsize=14, fontweight='bold', pad=30)

    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig12_summary_table.png')
    plt.close()

    # Also save as TSV
    pd.DataFrame(rows, columns=cols).to_csv(OUTDIR / 'summary_table.tsv', sep='\t', index=False)
    print("  fig12 done")

# ── Main ────────────────────────────────────────────────────
def main():
    to = load_data()

    print("\nGenerating figures...")
    fig01_loh_definitions_overview(to)
    auc_mat = fig02_auc_heatmap_dual_loh(to)
    fig03_per_sample_allele_auc(to)
    fig04_per_sample_hpfinep_auc(to)
    fig05_top10_auc_comparison(to)
    fig06_hp_only_quadrant_analysis(to)
    fig07_composition_effect(to)
    fig08_hcc1395_detailed(to)
    fig09_combo_auc(to)
    fig10_quality_score_analysis(to)
    fig11_evidence_chain(to)
    fig12_summary_table(to)

    # Save AUC matrix
    auc_mat.to_csv(OUTDIR / 'auc_matrix_dual_loh.tsv', sep='\t')

    print(f"\nAll figures saved to {OUTDIR}/")
    print("Done!")

if __name__ == '__main__':
    main()
