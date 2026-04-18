#!/usr/bin/env python3
"""
TO LOH Individual Feature Plots
================================
Generates 19 individual feature distribution figures:
  - 15 single-feature 2x2 panels (LOH.bed T/F x HP_Imb T/F)
  - 4 per-sample comparison plots (7 samples x 2 LOH strata)

Uses the same master dataset and conventions as build_to_loh_dual_definition_report.py.
"""

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns
from sklearn.metrics import roc_auc_score
from pathlib import Path

# ── Config ──────────────────────────────────────────────────
MASTER = '/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz'
OUTDIR = Path('/big7_disk/liaoyoyo2001/InterSubMod/output/to_loh_dual_definition_report/individual')
OUTDIR.mkdir(parents=True, exist_ok=True)

SAMPLES_ORDER = ['HCC1395', 'HCC1395_DORADO', 'HCC1937', 'HCC1954', 'COLO829', 'H1437', 'H2009']

TP_COLOR = '#2196F3'
FP_COLOR = '#F44336'

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
    except Exception:
        return np.nan


def load_data():
    """Load master dataset, filter to TO mode, create labels."""
    print("Loading master dataset...")
    df = pd.read_csv(MASTER, sep='\t', low_memory=False)
    to = df[df['mode'] == 'to'].copy()
    to['is_tp'] = (to['truth_label'] == 'TP').astype(int)
    to['label'] = to['truth_label']
    # Ensure boolean columns
    for col in ['Potential_LOH', 'to_loh_bed_hit']:
        if col in to.columns:
            to[col] = to[col].astype(bool)
    print(f"  TO: {len(to):,} rows, {to['is_tp'].sum():,} TP, {(1-to['is_tp']).sum():,.0f} FP")
    return to


def _annotate_panel(ax, subset, feature, is_pvalue=False):
    """Add N_TP, N_FP, medians, and AUC annotation to a panel."""
    vals = pd.to_numeric(subset[feature], errors='coerce')
    tp_vals = vals[subset['is_tp'] == 1].dropna()
    fp_vals = vals[subset['is_tp'] == 0].dropna()

    score_vals = -vals if is_pvalue else vals
    auc = safe_auc(subset['is_tp'], score_vals)

    text = (
        f"N_TP={len(tp_vals):,}  N_FP={len(fp_vals):,}\n"
        f"med_TP={tp_vals.median():.3f}  med_FP={fp_vals.median():.3f}\n"
        f"AUC={auc:.3f}" if not np.isnan(auc) else
        f"N_TP={len(tp_vals):,}  N_FP={len(fp_vals):,}\n"
        f"med_TP={tp_vals.median():.3f}  med_FP={fp_vals.median():.3f}\n"
        f"AUC=N/A"
    )
    ax.text(0.98, 0.98, text, transform=ax.transAxes, fontsize=7,
            va='top', ha='right', bbox=dict(boxstyle='round,pad=0.3',
            facecolor='white', alpha=0.8, edgecolor='gray'))


def _get_subsets(to):
    """Return the 4 LOH-definition subsets with labels."""
    return [
        ('LOH.bed=T', to[to['to_loh_bed_hit'] == True]),
        ('LOH.bed=F', to[to['to_loh_bed_hit'] == False]),
        ('HP_Imb=T', to[to['Potential_LOH'] == True]),
        ('HP_Imb=F', to[to['Potential_LOH'] == False]),
    ]


# ── Plot types ──────────────────────────────────────────────
def _plot_violin(ax, subset, feature, clip_range=None):
    """Split violin plot: TP left half, FP right half, with box overlay."""
    vals = pd.to_numeric(subset[feature], errors='coerce')
    plot_df = pd.DataFrame({
        'value': vals,
        'label': subset['label'].values,
    }).dropna(subset=['value'])

    if clip_range is not None:
        plot_df['value'] = plot_df['value'].clip(clip_range[0], clip_range[1])

    if len(plot_df) < 10:
        ax.text(0.5, 0.5, 'N < 10', transform=ax.transAxes, ha='center')
        return

    # Use split violin
    try:
        sns.violinplot(data=plot_df, x=['']*len(plot_df), y='value', hue='label',
                       split=True, hue_order=['TP', 'FP'],
                       palette={'TP': TP_COLOR, 'FP': FP_COLOR},
                       inner='box', density_norm='width', ax=ax, linewidth=0.8)
    except Exception:
        # Fallback to overlaid histograms if violin fails
        tp_v = plot_df[plot_df['label'] == 'TP']['value']
        fp_v = plot_df[plot_df['label'] == 'FP']['value']
        ax.hist(tp_v, bins=50, alpha=0.5, color=TP_COLOR, label='TP', density=True)
        ax.hist(fp_v, bins=50, alpha=0.5, color=FP_COLOR, label='FP', density=True)

    ax.set_xlabel('')
    ax.set_ylabel(feature)
    if ax.get_legend():
        ax.get_legend().set_visible(False)


def _plot_histogram_log10(ax, subset, feature, clip_max=10):
    """-log10 transformed histogram for p-value features."""
    vals = pd.to_numeric(subset[feature], errors='coerce').dropna()
    labels = subset.loc[vals.index, 'label']

    # Apply -log10 with clipping
    vals_log = -np.log10(vals.clip(lower=1e-300))
    vals_log = vals_log.clip(upper=clip_max)

    tp_v = vals_log[labels == 'TP']
    fp_v = vals_log[labels == 'FP']

    if len(tp_v) < 5 or len(fp_v) < 5:
        ax.text(0.5, 0.5, 'N < 5', transform=ax.transAxes, ha='center')
        return

    bins = np.linspace(0, clip_max, 40)
    ax.hist(tp_v, bins=bins, alpha=0.6, color=TP_COLOR, label='TP', density=True)
    ax.hist(fp_v, bins=bins, alpha=0.6, color=FP_COLOR, label='FP', density=True)
    ax.set_xlabel(f'-log10({feature})')
    ax.set_ylabel('Density')


def _plot_grouped_bar(ax, subset, feature, max_val=4):
    """Grouped bar chart for discrete features: proportion at each value."""
    vals = pd.to_numeric(subset[feature], errors='coerce').dropna()
    labels = subset.loc[vals.index, 'label']

    tp_vals = vals[labels == 'TP']
    fp_vals = vals[labels == 'FP']

    all_categories = sorted(vals.unique())
    # Limit to reasonable range
    all_categories = [c for c in all_categories if c <= max_val + 2]

    tp_counts = tp_vals.value_counts().reindex(all_categories, fill_value=0)
    fp_counts = fp_vals.value_counts().reindex(all_categories, fill_value=0)

    tp_prop = tp_counts / tp_counts.sum() if tp_counts.sum() > 0 else tp_counts
    fp_prop = fp_counts / fp_counts.sum() if fp_counts.sum() > 0 else fp_counts

    x = np.arange(len(all_categories))
    width = 0.35
    ax.bar(x - width/2, tp_prop.values, width, color=TP_COLOR, alpha=0.8, label='TP')
    ax.bar(x + width/2, fp_prop.values, width, color=FP_COLOR, alpha=0.8, label='FP')
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(c)) for c in all_categories], fontsize=8)
    ax.set_xlabel(feature)
    ax.set_ylabel('Proportion')


# ── Individual feature figure (2x2 panels) ──────────────────
def generate_individual_feature(to, feature, plot_type, clip_range=None,
                                clip_max=10, discrete_max=4, filename=None):
    """Generate a 2x2 figure for one feature across 4 LOH subsets."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    subsets = _get_subsets(to)

    is_pvalue = plot_type == 'histogram'

    for idx, (label, subset) in enumerate(subsets):
        row, col = divmod(idx, 2)
        ax = axes[row, col]

        if plot_type == 'violin':
            _plot_violin(ax, subset, feature, clip_range=clip_range)
        elif plot_type == 'histogram':
            _plot_histogram_log10(ax, subset, feature, clip_max=clip_max)
        elif plot_type == 'grouped_bar':
            _plot_grouped_bar(ax, subset, feature, max_val=discrete_max)

        ax.set_title(f'{label}  (N={len(subset):,})', fontsize=12, fontweight='bold')
        _annotate_panel(ax, subset, feature, is_pvalue=is_pvalue)

    # Shared legend
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=TP_COLOR, alpha=0.8, label='TP'),
        plt.Rectangle((0, 0), 1, 1, fc=FP_COLOR, alpha=0.8, label='FP'),
    ]
    fig.legend(handles=handles, loc='upper right', fontsize=11, framealpha=0.9)

    fig.suptitle(f'{feature} — TP vs FP Distribution by LOH Definition',
                 fontsize=15, fontweight='bold', y=1.01)
    fig.tight_layout()

    out_path = OUTDIR / (filename or f'ind_{feature}.png')
    fig.savefig(out_path)
    plt.close()
    print(f"  {out_path.name} done")


# ── Per-sample comparison figure (2 rows x 7 cols) ──────────
def generate_per_sample_comparison(to, feature, filename, is_pvalue=False,
                                   clip_range=None):
    """Generate 2-row (LOH/Non-LOH) x 7-col (samples) comparison."""
    fig, axes = plt.subplots(2, 7, figsize=(20, 8), sharey='row')

    row_defs = [
        ('LOH.bed=T', to['to_loh_bed_hit'] == True),
        ('LOH.bed=F', to['to_loh_bed_hit'] == False),
    ]

    for row_idx, (row_label, row_mask) in enumerate(row_defs):
        for col_idx, sample in enumerate(SAMPLES_ORDER):
            ax = axes[row_idx, col_idx]
            subset = to[row_mask & (to['sample'] == sample)]

            vals = pd.to_numeric(subset[feature], errors='coerce')
            if is_pvalue:
                vals = -np.log10(vals.clip(lower=1e-300)).clip(upper=10)
                feat_label = f'-log10({feature})'
            else:
                if clip_range is not None:
                    vals = vals.clip(clip_range[0], clip_range[1])
                feat_label = feature

            plot_df = pd.DataFrame({
                'value': vals,
                'label': subset['label'].values,
            }).dropna(subset=['value'])

            n_tp = (plot_df['label'] == 'TP').sum()
            n_fp = (plot_df['label'] == 'FP').sum()

            if len(plot_df) >= 10 and plot_df['label'].nunique() == 2:
                try:
                    sns.violinplot(data=plot_df, x=['']*len(plot_df), y='value',
                                   hue='label', split=True, hue_order=['TP', 'FP'],
                                   palette={'TP': TP_COLOR, 'FP': FP_COLOR},
                                   inner='box', density_norm='width', ax=ax,
                                   linewidth=0.5)
                except Exception:
                    tp_v = plot_df[plot_df['label'] == 'TP']['value']
                    fp_v = plot_df[plot_df['label'] == 'FP']['value']
                    ax.hist(tp_v, bins=30, alpha=0.5, color=TP_COLOR, density=True)
                    ax.hist(fp_v, bins=30, alpha=0.5, color=FP_COLOR, density=True)
            else:
                ax.text(0.5, 0.5, f'N={len(plot_df)}', transform=ax.transAxes,
                        ha='center', fontsize=8)

            # AUC annotation
            score = -pd.to_numeric(subset[feature], errors='coerce') if is_pvalue \
                    else pd.to_numeric(subset[feature], errors='coerce')
            auc = safe_auc(subset['is_tp'], score)
            auc_str = f'{auc:.3f}' if not np.isnan(auc) else 'N/A'
            ax.set_title(f'{sample}\n{row_label}' if row_idx == 0 else row_label,
                         fontsize=7, fontweight='bold')
            ax.text(0.5, 0.02, f'AUC={auc_str}\nTP={n_tp} FP={n_fp}',
                    transform=ax.transAxes, fontsize=6, ha='center', va='bottom',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray', pad=1))

            if ax.get_legend():
                ax.get_legend().set_visible(False)
            ax.set_xlabel('')
            if col_idx > 0:
                ax.set_ylabel('')
            else:
                ax.set_ylabel(feat_label, fontsize=8)
            ax.tick_params(axis='both', labelsize=6)

    # Global legend
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=TP_COLOR, alpha=0.8, label='TP'),
        plt.Rectangle((0, 0), 1, 1, fc=FP_COLOR, alpha=0.8, label='FP'),
    ]
    fig.legend(handles=handles, loc='upper right', fontsize=10, framealpha=0.9)

    fig.suptitle(f'{feature} — Per-Sample TP vs FP (LOH.bed LOH / Non-LOH)',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(OUTDIR / filename)
    plt.close()
    print(f"  {filename} done")


# ── Main ────────────────────────────────────────────────────
def main():
    to = load_data()

    print("\nGenerating 15 individual feature figures...")

    # 1. AlleleDelta — violin, range [-1, 1]
    generate_individual_feature(to, 'AlleleDelta', 'violin', clip_range=(-1, 1),
                                filename='ind_AlleleDelta.png')
    # 2. HPFineP — histogram, -log10
    generate_individual_feature(to, 'HPFineP', 'histogram', clip_max=10,
                                filename='ind_HPFineP.png')
    # 3. HPFineNGroups — grouped bar (discrete 0-4)
    generate_individual_feature(to, 'HPFineNGroups', 'grouped_bar', discrete_max=4,
                                filename='ind_HPFineNGroups.png')
    # 4. PairwiseMeanDist — violin
    generate_individual_feature(to, 'PairwiseMeanDist', 'violin',
                                filename='ind_PairwiseMeanDist.png')
    # 5. CramersV — violin, range [0, 1]
    generate_individual_feature(to, 'CramersV', 'violin', clip_range=(0, 1),
                                filename='ind_CramersV.png')
    # 6. Quality_Score — violin, range [0, 100]
    generate_individual_feature(to, 'Quality_Score', 'violin', clip_range=(0, 100),
                                filename='ind_Quality_Score.png')
    # 7. caller_af — violin, range [0, 1]
    generate_individual_feature(to, 'caller_af', 'violin', clip_range=(0, 1),
                                filename='ind_caller_af.png')
    # 8. NumReads — violin, clip at 500
    generate_individual_feature(to, 'NumReads', 'violin', clip_range=(0, 500),
                                filename='ind_NumReads.png')
    # 9. HP_Ratio — violin, range [0, 1]
    generate_individual_feature(to, 'HP_Ratio', 'violin', clip_range=(0, 1),
                                filename='ind_HP_Ratio.png')
    # 10. Coverage_Multiple — violin, clip at 5
    generate_individual_feature(to, 'Coverage_Multiple', 'violin', clip_range=(0, 5),
                                filename='ind_Coverage_Multiple.png')
    # 11. HPMergedDelta — violin
    generate_individual_feature(to, 'HPMergedDelta', 'violin',
                                filename='ind_HPMergedDelta.png')
    # 12. HP1FamilyN — violin, clip at 300
    generate_individual_feature(to, 'HP1FamilyN', 'violin', clip_range=(0, 300),
                                filename='ind_HP1FamilyN.png')
    # 13. HP2FamilyN — violin, clip at 300
    generate_individual_feature(to, 'HP2FamilyN', 'violin', clip_range=(0, 300),
                                filename='ind_HP2FamilyN.png')
    # 14. GlobalP — histogram, -log10
    generate_individual_feature(to, 'GlobalP', 'histogram', clip_max=10,
                                filename='ind_GlobalP.png')
    # 15. HeuristicScore — violin
    generate_individual_feature(to, 'HeuristicScore', 'violin',
                                filename='ind_HeuristicScore.png')

    print("\nGenerating 4 per-sample comparison figures...")

    # 16. AlleleDelta by sample
    generate_per_sample_comparison(to, 'AlleleDelta', 'cmp_AlleleDelta_by_sample.png',
                                   clip_range=(-1, 1))
    # 17. HPFineP by sample (-log10)
    generate_per_sample_comparison(to, 'HPFineP', 'cmp_HPFineP_by_sample.png',
                                   is_pvalue=True)
    # 18. caller_af by sample
    generate_per_sample_comparison(to, 'caller_af', 'cmp_caller_af_by_sample.png',
                                   clip_range=(0, 1))
    # 19. HPFineNGroups by sample
    generate_per_sample_comparison(to, 'HPFineNGroups', 'cmp_HPFineNGroups_by_sample.png')

    print(f"\nAll 19 figures saved to {OUTDIR}/")
    print("Done!")


if __name__ == '__main__':
    main()
