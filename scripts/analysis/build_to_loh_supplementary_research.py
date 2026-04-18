#!/usr/bin/env python3
"""
TO LOH Supplementary Research — C1/C2/C3 Analyses
==================================================
C1: caller_af Reversal Analysis
C2: HP_only Quadrant Deep Dive
C3: Masking Strategy Analysis

Outputs go to: output/to_loh_dual_definition_report/supplementary/
"""

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
import seaborn as sns
from sklearn.metrics import roc_auc_score
from pathlib import Path

# ── CJK Font Setup ──────────────────────────────────────────
def _get_cjk_font():
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
OUTDIR = Path('/big7_disk/liaoyoyo2001/InterSubMod/output/to_loh_dual_definition_report/supplementary')
OUTDIR.mkdir(parents=True, exist_ok=True)

SAMPLES_ORDER = ['HCC1395', 'HCC1395_DORADO', 'HCC1937', 'HCC1954', 'COLO829', 'H1437', 'H2009']

# Colors
C_TP = '#2196F3'
C_FP = '#F44336'
C_BOTH = '#1565C0'
C_HP_ONLY = '#FF8F00'
C_NEITHER = '#757575'

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
    print("Loading master dataset...")
    df = pd.read_csv(MASTER, sep='\t', low_memory=False)
    to = df[df['mode'] == 'to'].copy()
    to['is_tp'] = (to['truth_label'] == 'TP').astype(int)
    # Ensure boolean columns
    for col in ['Potential_LOH', 'to_loh_bed_hit', 'PassedGating']:
        if col in to.columns:
            to[col] = to[col].astype(bool)
    # Create 4-quadrant LOH classification
    to['LOH_quadrant'] = 'Neither'
    to.loc[to['Potential_LOH'] & to['to_loh_bed_hit'], 'LOH_quadrant'] = 'Both'
    to.loc[to['Potential_LOH'] & ~to['to_loh_bed_hit'], 'LOH_quadrant'] = 'HP_only'
    to.loc[~to['Potential_LOH'] & to['to_loh_bed_hit'], 'LOH_quadrant'] = 'Bed_only'
    print(f"  TO: {len(to):,} rows, {to['is_tp'].sum():,} TP, {(1-to['is_tp']).sum():,.0f} FP")
    return to


# ── LOH definition masks ───────────────────────────────────
def get_loh_defs(to):
    return {
        'LOH.bed=T': to['to_loh_bed_hit'] == True,
        'LOH.bed=F': to['to_loh_bed_hit'] == False,
        'HP_Imb=T': to['Potential_LOH'] == True,
        'HP_Imb=F': to['Potential_LOH'] == False,
    }


# ═════════════════════════════════════════════════════════════
# C1: caller_af Reversal Analysis
# ═════════════════════════════════════════════════════════════

def c1_caller_af_distribution(to):
    """C1-1: caller_af distribution histograms across LOH definitions."""
    print("  C1-1: caller_af distribution...")
    loh_defs = get_loh_defs(to)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    panels = [('LOH.bed=T', axes[0, 0]), ('LOH.bed=F', axes[0, 1]),
              ('HP_Imb=T', axes[1, 0]), ('HP_Imb=F', axes[1, 1])]

    for def_name, ax in panels:
        mask = loh_defs[def_name]
        sub = to[mask]
        tp_af = sub.loc[sub['is_tp'] == 1, 'caller_af'].dropna()
        fp_af = sub.loc[sub['is_tp'] == 0, 'caller_af'].dropna()

        ax.hist(tp_af, bins=50, alpha=0.5, color=C_TP, label=f'TP (N={len(tp_af):,})', density=True)
        ax.hist(fp_af, bins=50, alpha=0.5, color=C_FP, label=f'FP (N={len(fp_af):,})', density=True)

        # Annotate stats
        tp_med, tp_mean = tp_af.median(), tp_af.mean()
        fp_med, fp_mean = fp_af.median(), fp_af.mean()
        stats_text = (f'TP: med={tp_med:.3f}, mean={tp_mean:.3f}\n'
                      f'FP: med={fp_med:.3f}, mean={fp_mean:.3f}')
        ax.text(0.02, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                va='top', ha='left', bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

        ax.set_title(f'{def_name}')
        ax.set_xlabel('caller_af')
        ax.set_ylabel('Density')
        ax.legend(fontsize=9, loc='upper left', bbox_to_anchor=(0.0, 0.75))

    fig.suptitle('C1-1: caller_af Distribution by TP/FP across LOH Definitions', fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUTDIR / 'sup_caller_af_distribution_loh.png')
    plt.close()
    print("    saved sup_caller_af_distribution_loh.png")


def c1_threshold_sweep(to):
    """C1-2: caller_af threshold sweep — removal analysis."""
    print("  C1-2: caller_af threshold sweep...")
    loh_defs = get_loh_defs(to)
    thresholds = [0.80, 0.85, 0.90, 0.95, 0.98, 0.99, 1.00]

    all_results = []

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    panels = [('LOH.bed=T', axes[0, 0]), ('LOH.bed=F', axes[0, 1]),
              ('HP_Imb=T', axes[1, 0]), ('HP_Imb=F', axes[1, 1])]

    for def_name, ax in panels:
        mask = loh_defs[def_name]
        sub = to[mask].copy()
        n_total = len(sub)
        n_fp_total = (sub['is_tp'] == 0).sum()
        n_tp_total = (sub['is_tp'] == 1).sum()

        fp_pcts = []
        tp_pcts = []
        remaining_fp_rates = []

        for thr in thresholds:
            removed_mask = sub['caller_af'] >= thr
            n_removed = removed_mask.sum()
            n_fp_removed = (removed_mask & (sub['is_tp'] == 0)).sum()
            n_tp_lost = (removed_mask & (sub['is_tp'] == 1)).sum()

            remaining = sub[~removed_mask]
            n_remaining = len(remaining)
            remaining_fp = (remaining['is_tp'] == 0).sum() if n_remaining > 0 else 0
            remaining_fp_rate = remaining_fp / n_remaining * 100 if n_remaining > 0 else 0

            fp_removal_pct = n_fp_removed / n_fp_total * 100 if n_fp_total > 0 else 0
            tp_loss_pct = n_tp_lost / n_tp_total * 100 if n_tp_total > 0 else 0

            fp_pcts.append(fp_removal_pct)
            tp_pcts.append(tp_loss_pct)
            remaining_fp_rates.append(remaining_fp_rate)

            all_results.append({
                'Definition': def_name,
                'Threshold': thr,
                'N_removed': int(n_removed),
                'N_FP_removed': int(n_fp_removed),
                'N_TP_lost': int(n_tp_lost),
                'FP_removal_rate': round(fp_removal_pct, 2),
                'TP_loss_rate': round(tp_loss_pct, 2),
                'Remaining_FP_rate': round(remaining_fp_rate, 2),
                'Remaining_N': int(n_remaining),
            })

        # Plot
        ax.plot(thresholds, fp_pcts, 'o-', color=C_FP, linewidth=2, markersize=7, label='% FP removed')
        ax.plot(thresholds, tp_pcts, 's-', color=C_TP, linewidth=2, markersize=7, label='% TP lost')
        ax.set_xlabel('caller_af threshold (remove if >= thr)')
        ax.set_ylabel('% removed')
        ax.legend(loc='upper left', fontsize=9)

        ax2 = ax.twinx()
        ax2.plot(thresholds, remaining_fp_rates, 'D--', color='#4CAF50', linewidth=1.5, markersize=5, label='Remaining FP rate')
        ax2.set_ylabel('Remaining FP rate (%)', color='#4CAF50')
        ax2.tick_params(axis='y', labelcolor='#4CAF50')
        ax2.legend(loc='upper right', fontsize=9)

        ax.set_title(f'{def_name} (N={n_total:,})')
        ax.grid(axis='y', alpha=0.3)

        # Add compact table below the plot area — use text annotation
        table_text = 'Thr   | FP_rm% | TP_loss% | FP_rate%\n'
        for i, thr in enumerate(thresholds):
            table_text += f'{thr:.2f}  | {fp_pcts[i]:5.1f}  | {tp_pcts[i]:6.1f}   | {remaining_fp_rates[i]:5.1f}\n'
        ax.text(0.02, -0.25, table_text, transform=ax.transAxes, fontsize=7,
                va='top', ha='left', family='monospace',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    fig.suptitle('C1-2: caller_af Threshold Sweep — FP Removal vs TP Loss', fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0.08, 1, 0.96])
    fig.savefig(OUTDIR / 'sup_caller_af_threshold_sweep.png')
    plt.close()
    print("    saved sup_caller_af_threshold_sweep.png")

    # Save results TSV
    res_df = pd.DataFrame(all_results)
    res_df.to_csv(OUTDIR / 'sup_caller_af_results.tsv', sep='\t', index=False)
    print(f"    saved sup_caller_af_results.tsv ({len(res_df)} rows)")
    return res_df


def c1_af1_by_sample(to):
    """C1-3: AF=1.0 prevalence by sample, split by LOH definitions."""
    print("  C1-3: AF=1.0 by sample...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    for ax, (def_col, def_label) in zip(axes, [('to_loh_bed_hit', 'LOH.bed'),
                                                 ('Potential_LOH', 'HP Imbalance')]):
        rows = []
        for s in SAMPLES_ORDER:
            for loh_val, loh_label in [(True, f'{def_label}=T'), (False, f'{def_label}=F')]:
                sub = to[(to['sample'] == s) & (to[def_col] == loh_val)]
                tp_sub = sub[sub['is_tp'] == 1]
                fp_sub = sub[sub['is_tp'] == 0]
                tp_af1 = (tp_sub['caller_af'] == 1.0).sum() / len(tp_sub) * 100 if len(tp_sub) > 0 else 0
                fp_af1 = (fp_sub['caller_af'] == 1.0).sum() / len(fp_sub) * 100 if len(fp_sub) > 0 else 0
                rows.append({'Sample': s, 'LOH': loh_label, 'TP_AF1_pct': tp_af1, 'FP_AF1_pct': fp_af1})

        rdf = pd.DataFrame(rows)

        # Split by LOH definition
        x = np.arange(len(SAMPLES_ORDER))
        width = 0.18

        for i, loh_label in enumerate([f'{def_label}=T', f'{def_label}=F']):
            sub_r = rdf[rdf['LOH'] == loh_label]
            offset = (i - 0.5) * 2 * width
            tp_vals = sub_r['TP_AF1_pct'].values
            fp_vals = sub_r['FP_AF1_pct'].values

            bars_tp = ax.bar(x + offset - width/2, tp_vals, width, label=f'{loh_label} TP',
                            color=C_TP, alpha=0.6 + i*0.2, edgecolor='white')
            bars_fp = ax.bar(x + offset + width/2, fp_vals, width, label=f'{loh_label} FP',
                            color=C_FP, alpha=0.6 + i*0.2, edgecolor='white')

        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLES_ORDER, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('% with caller_af == 1.0')
        ax.set_title(f'AF=1.0 Prevalence by {def_label}')
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(axis='y', alpha=0.3)

    fig.suptitle('C1-3: Exact AF=1.0 Prevalence by Sample and LOH Definition', fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUTDIR / 'sup_caller_af_af1_by_sample.png')
    plt.close()
    print("    saved sup_caller_af_af1_by_sample.png")


# ═════════════════════════════════════════════════════════════
# C2: HP_only Quadrant Deep Dive
# ═════════════════════════════════════════════════════════════

def c2_feature_comparison(to):
    """C2-1: Feature distributions across 3 quadrants (Both, HP_only, Neither)."""
    print("  C2-1: HP_only feature comparison...")
    features = ['AlleleDelta', 'HPFineP', 'HPFineNGroups', 'PairwiseMeanDist',
                'CramersV', 'Quality_Score', 'caller_af', 'NumReads',
                'HP_Ratio', 'Coverage_Multiple', 'HPMergedDelta', 'NumCpGs']

    fig, axes = plt.subplots(3, 4, figsize=(18, 14))
    quad_colors = {'Both': C_BOTH, 'HP_only': C_HP_ONLY, 'Neither': C_NEITHER}
    quads = ['Both', 'HP_only', 'Neither']

    for idx, feat in enumerate(features):
        row, col_idx = divmod(idx, 4)
        ax = axes[row, col_idx]

        if feat not in to.columns:
            ax.set_title(f'{feat} (N/A)')
            ax.text(0.5, 0.5, 'Column not found', ha='center', va='center', transform=ax.transAxes)
            continue

        for q in quads:
            vals = pd.to_numeric(to.loc[to['LOH_quadrant'] == q, feat], errors='coerce').dropna()
            if feat == 'HPFineP':
                vals = -np.log10(vals.clip(lower=1e-300))
                label_feat = '-log10(HPFineP)'
            else:
                label_feat = feat

            if len(vals) == 0:
                continue

            # Clip extreme outliers for visualization
            p1, p99 = vals.quantile(0.01), vals.quantile(0.99)
            vals_clipped = vals.clip(p1, p99)

            ax.hist(vals_clipped, bins=50, alpha=0.4, color=quad_colors[q],
                    label=f'{q} (med={vals.median():.3f})', density=True)

        if feat == 'HPFineP':
            ax.set_xlabel('-log10(HPFineP)')
        else:
            ax.set_xlabel(feat)
        ax.set_ylabel('Density')
        ax.set_title(feat)
        ax.legend(fontsize=7, loc='upper right')

    fig.suptitle('C2-1: Feature Distribution by LOH Quadrant (Both / HP_only / Neither)',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUTDIR / 'sup_hp_only_feature_comparison.png')
    plt.close()
    print("    saved sup_hp_only_feature_comparison.png")


def c2_auc_comparison(to):
    """C2-2: AUC comparison across 3 quadrants."""
    print("  C2-2: HP_only AUC comparison...")
    features = ['AlleleDelta', 'HPFineP', 'HPFineNGroups', 'PairwiseMeanDist',
                'CramersV', 'Quality_Score', 'caller_af', 'NumReads']

    quads = ['Both', 'HP_only', 'Neither']
    quad_colors = {'Both': C_BOTH, 'HP_only': C_HP_ONLY, 'Neither': C_NEITHER}

    results = []

    for feat in features:
        if feat not in to.columns:
            continue
        for q in quads:
            sub = to[to['LOH_quadrant'] == q]
            vals = pd.to_numeric(sub[feat], errors='coerce')
            if feat == 'HPFineP':
                vals = -vals  # negate for AUC (lower p = more significant = TP)
            auc = safe_auc(sub['is_tp'], vals)
            n_tp = sub['is_tp'].sum()
            n_fp = (sub['is_tp'] == 0).sum()
            results.append({
                'Feature': feat if feat != 'HPFineP' else '-HPFineP',
                'Quadrant': q,
                'AUC': auc,
                'N_TP': int(n_tp),
                'N_FP': int(n_fp),
                'N_total': int(len(sub)),
                'FP_rate': round(n_fp / len(sub) * 100, 2) if len(sub) > 0 else 0,
            })

    res_df = pd.DataFrame(results)

    # Plot grouped bar chart
    fig, ax = plt.subplots(figsize=(14, 8))

    feat_labels = res_df['Feature'].unique()
    x = np.arange(len(feat_labels))
    width = 0.25

    for i, q in enumerate(quads):
        q_data = res_df[res_df['Quadrant'] == q]
        auc_vals = [q_data.loc[q_data['Feature'] == f, 'AUC'].values[0]
                    if f in q_data['Feature'].values else np.nan for f in feat_labels]
        bars = ax.bar(x + (i - 1) * width, auc_vals, width, label=q,
                     color=quad_colors[q], alpha=0.85, edgecolor='white')
        # Add value labels
        for j, v in enumerate(auc_vals):
            if not np.isnan(v):
                ax.text(x[j] + (i - 1) * width, v + 0.005, f'{v:.3f}',
                       ha='center', va='bottom', fontsize=7, rotation=90)

    ax.axhline(y=0.50, ls='--', color='gray', alpha=0.6, linewidth=1, label='Random (0.50)')
    ax.axhline(y=0.58, ls='--', color='orange', alpha=0.6, linewidth=1, label='Threshold (0.58)')

    ax.set_xticks(x)
    ax.set_xticklabels(feat_labels, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('AUC')
    ax.set_title('C2-2: AUC Comparison — Both vs HP_only vs Neither')
    ax.legend(fontsize=9)
    ax.set_ylim(0.35, max(0.80, res_df['AUC'].max() + 0.05))
    ax.grid(axis='y', alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUTDIR / 'sup_hp_only_auc_comparison.png')
    plt.close()
    print("    saved sup_hp_only_auc_comparison.png")

    # Save results TSV
    res_df.to_csv(OUTDIR / 'sup_hp_only_results.tsv', sep='\t', index=False)
    print(f"    saved sup_hp_only_results.tsv ({len(res_df)} rows)")
    return res_df


# ═════════════════════════════════════════════════════════════
# C3: Masking Strategy Analysis
# ═════════════════════════════════════════════════════════════

def c3_masking_analysis(to):
    """C3: Masking strategy — exclude subsets and measure AUC improvement."""
    print("  C3: Masking strategy analysis...")

    # Check which column name exists for PassedGating
    gating_col = None
    for candidate in ['PassedGating', 'passed_gating']:
        if candidate in to.columns:
            gating_col = candidate
            break
    if gating_col is None:
        print("    WARNING: Neither PassedGating nor passed_gating found; mask 4 will be empty.")

    # Define masks (True = EXCLUDE this region)
    def get_masks(sub):
        masks = {}
        masks['M1: NumReads<20'] = sub['NumReads'] < 20
        masks['M2: CNV_anomaly'] = (sub['Coverage_Multiple'] < 0.5) | (sub['Coverage_Multiple'] > 2.0)
        masks['M3: HPFineNGroups==0'] = sub['HPFineNGroups'] == 0
        if gating_col is not None:
            masks['M4: !PassedGating'] = ~sub[gating_col].astype(bool)
        else:
            masks['M4: !PassedGating'] = pd.Series(False, index=sub.index)
        masks['M5: NumCpGs<5'] = sub['NumCpGs'] < 5
        masks['M6: AF>0.95'] = sub['caller_af'] > 0.95
        masks['M7: HP_balanced'] = (sub['HP_Ratio'] < 0.3) | (sub['HP_Ratio'] > 0.7)
        masks['M8: M1+M2+M3'] = masks['M1: NumReads<20'] | masks['M2: CNV_anomaly'] | masks['M3: HPFineNGroups==0']
        masks['M9: M1+M2+M3+M5'] = masks['M8: M1+M2+M3'] | masks['M5: NumCpGs<5']
        masks['M10: M1+M2+M3+M5+M6'] = masks['M9: M1+M2+M3+M5'] | masks['M6: AF>0.95']
        return masks

    loh_defs = get_loh_defs(to)
    auc_features = ['AlleleDelta', 'HPFineNGroups', 'HPFineP', 'Quality_Score', 'caller_af']
    auc_labels = ['AlleleDelta', 'HPFineNGroups', '-HPFineP', 'Quality_Score', 'caller_af']

    all_results = []

    for def_name, def_mask in loh_defs.items():
        sub = to[def_mask].copy()
        n_total = len(sub)
        fp_rate_original = (sub['is_tp'] == 0).sum() / n_total * 100 if n_total > 0 else 0

        # Baseline (no masking)
        row = {
            'Mask': 'Baseline',
            'Definition': def_name,
            'N_total': n_total,
            'N_excluded': 0,
            'Pct_excluded': 0.0,
            'FP_rate': round(fp_rate_original, 3),
            'FP_rate_change': 0.0,
        }
        for feat, label in zip(auc_features, auc_labels):
            vals = pd.to_numeric(sub[feat], errors='coerce')
            if label == '-HPFineP':
                vals = -vals
            row[f'AUC_{label}'] = round(safe_auc(sub['is_tp'], vals), 4)
        all_results.append(row)

        # Each mask
        masks = get_masks(sub)
        for mask_name, exclude_mask in masks.items():
            remaining = sub[~exclude_mask]
            n_excluded = exclude_mask.sum()
            n_remaining = len(remaining)
            pct_excluded = n_excluded / n_total * 100 if n_total > 0 else 0
            fp_rate_new = (remaining['is_tp'] == 0).sum() / n_remaining * 100 if n_remaining > 0 else 0

            row = {
                'Mask': mask_name,
                'Definition': def_name,
                'N_total': n_remaining,
                'N_excluded': int(n_excluded),
                'Pct_excluded': round(pct_excluded, 2),
                'FP_rate': round(fp_rate_new, 3),
                'FP_rate_change': round(fp_rate_new - fp_rate_original, 3),
            }
            for feat, label in zip(auc_features, auc_labels):
                vals = pd.to_numeric(remaining[feat], errors='coerce')
                if label == '-HPFineP':
                    vals = -vals
                row[f'AUC_{label}'] = round(safe_auc(remaining['is_tp'], vals), 4)
            all_results.append(row)

    res_df = pd.DataFrame(all_results)
    res_df.to_csv(OUTDIR / 'sup_masking_results.tsv', sep='\t', index=False)
    print(f"    saved sup_masking_results.tsv ({len(res_df)} rows)")

    # ── C3-2: AUC Heatmap ──────────────────────────────────
    print("  C3-2: Masking AUC heatmap...")
    # Build a wide matrix: rows=masks, cols=(feature × definition)
    mask_names = ['Baseline'] + list(get_masks(to.head(1)).keys())
    col_pairs = [(label, d) for d in loh_defs.keys() for label in auc_labels]
    col_names = [f'{label}\n{d}' for label, d in col_pairs]

    heatmap_data = np.full((len(mask_names), len(col_pairs)), np.nan)
    for i, mn in enumerate(mask_names):
        for j, (label, d) in enumerate(col_pairs):
            row_match = res_df[(res_df['Mask'] == mn) & (res_df['Definition'] == d)]
            if len(row_match) > 0:
                auc_col = f'AUC_{label}'
                if auc_col in row_match.columns:
                    heatmap_data[i, j] = row_match[auc_col].values[0]

    fig, ax = plt.subplots(figsize=(18, 12))
    cmap = sns.diverging_palette(10, 130, s=80, l=55, as_cmap=True)

    sns.heatmap(heatmap_data, ax=ax, cmap=cmap, center=0.50, vmin=0.40, vmax=0.70,
                annot=True, fmt='.3f', linewidths=0.5,
                xticklabels=col_names, yticklabels=mask_names, annot_kws={'fontsize': 8})
    ax.set_xlabel('Feature × LOH Definition')
    ax.set_ylabel('Masking Strategy')
    ax.set_title('C3-2: AUC by Masking Strategy — Feature × LOH Definition', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTDIR / 'sup_masking_auc_heatmap.png')
    plt.close()
    print("    saved sup_masking_auc_heatmap.png")

    # ── C3-3: Pareto frontier (LOH.bed=T) ──────────────────
    print("  C3-3: Masking Pareto frontier...")
    loh_bed_t = res_df[res_df['Definition'] == 'LOH.bed=T'].copy()
    loh_bed_t = loh_bed_t[loh_bed_t['Mask'] != 'Baseline']

    # For each mask: best AUC across all features
    auc_cols = [c for c in loh_bed_t.columns if c.startswith('AUC_')]
    loh_bed_t['Best_AUC'] = loh_bed_t[auc_cols].max(axis=1)

    fig, ax = plt.subplots(figsize=(16, 10))
    x_vals = loh_bed_t['Pct_excluded'].values
    y_vals = loh_bed_t['Best_AUC'].values
    labels = loh_bed_t['Mask'].values

    scatter = ax.scatter(x_vals, y_vals, s=120, c=range(len(x_vals)),
                        cmap='viridis', edgecolors='black', linewidth=0.5, zorder=5)

    for i, label in enumerate(labels):
        short_label = label.split(': ')[1] if ': ' in label else label
        ax.annotate(short_label, (x_vals[i], y_vals[i]),
                   xytext=(8, 8), textcoords='offset points', fontsize=8,
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', alpha=0.7))

    # Compute and highlight Pareto frontier
    # Pareto = maximize y (AUC) while minimizing x (exclusion %)
    # Sort by x ascending, then find points where y is max so far from right
    sorted_idx = np.argsort(x_vals)
    pareto_x = []
    pareto_y = []
    max_y = -np.inf
    # For Pareto: we want lowest exclusion for highest AUC
    # Go from right to left (most excluded), track best AUC
    for i in reversed(sorted_idx):
        if y_vals[i] >= max_y:
            max_y = y_vals[i]
            pareto_x.append(x_vals[i])
            pareto_y.append(y_vals[i])
    pareto_x.reverse()
    pareto_y.reverse()

    if len(pareto_x) > 1:
        ax.plot(pareto_x, pareto_y, 'r--', linewidth=2, alpha=0.7, label='Pareto frontier')

    ax.axhline(y=0.58, ls=':', color='orange', alpha=0.5, label='AUC=0.58 threshold')
    ax.axhline(y=0.50, ls=':', color='gray', alpha=0.5, label='AUC=0.50 random')

    ax.set_xlabel('% Regions Excluded', fontsize=12)
    ax.set_ylabel('Best AUC Achieved', fontsize=12)
    ax.set_title('C3-3: Masking Strategy Pareto — LOH.bed=T\n(x = cost in exclusion, y = best discriminability)',
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUTDIR / 'sup_masking_strategy.png')
    plt.close()
    print("    saved sup_masking_strategy.png")

    return res_df


# ═════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("TO LOH Supplementary Research — C1/C2/C3")
    print("=" * 70)

    to = load_data()

    # Quick overview
    print("\n--- Data Overview ---")
    for q in ['Both', 'HP_only', 'Bed_only', 'Neither']:
        sub = to[to['LOH_quadrant'] == q]
        n_tp = sub['is_tp'].sum()
        n_fp = (sub['is_tp'] == 0).sum()
        print(f"  {q:10s}: N={len(sub):>7,}  TP={n_tp:>6,}  FP={n_fp:>6,}  FP%={n_fp/len(sub)*100:.1f}%")

    loh_defs = get_loh_defs(to)
    for def_name, mask in loh_defs.items():
        sub = to[mask]
        print(f"  {def_name:12s}: N={len(sub):>7,}")

    # ── C1: caller_af Reversal ──
    print("\n=== C1: caller_af Reversal Analysis ===")
    c1_caller_af_distribution(to)
    af_results = c1_threshold_sweep(to)
    c1_af1_by_sample(to)

    # Key C1 findings
    print("\n  C1 Key Findings:")
    for def_name in ['LOH.bed=T', 'LOH.bed=F', 'HP_Imb=T', 'HP_Imb=F']:
        sub = to[get_loh_defs(to)[def_name]]
        tp_af = sub.loc[sub['is_tp'] == 1, 'caller_af']
        fp_af = sub.loc[sub['is_tp'] == 0, 'caller_af']
        print(f"    {def_name}: TP med={tp_af.median():.3f}, FP med={fp_af.median():.3f}, "
              f"TP AF=1.0: {(tp_af==1.0).mean()*100:.1f}%, FP AF=1.0: {(fp_af==1.0).mean()*100:.1f}%")

    # ── C2: HP_only Quadrant ──
    print("\n=== C2: HP_only Quadrant Deep Dive ===")
    c2_feature_comparison(to)
    c2_results = c2_auc_comparison(to)

    # Key C2 findings
    print("\n  C2 Key Findings (AUC by quadrant):")
    for feat in c2_results['Feature'].unique():
        vals = c2_results[c2_results['Feature'] == feat]
        line = f"    {feat:20s}: "
        for _, r in vals.iterrows():
            line += f"{r['Quadrant']}={r['AUC']:.3f}  "
        print(line)

    # ── C3: Masking Strategy ──
    print("\n=== C3: Masking Strategy Analysis ===")
    c3_results = c3_masking_analysis(to)

    # Key C3 findings
    print("\n  C3 Key Findings (LOH.bed=T, best AUC per mask):")
    auc_cols = [c for c in c3_results.columns if c.startswith('AUC_')]
    loh_bed_t = c3_results[c3_results['Definition'] == 'LOH.bed=T']
    for _, r in loh_bed_t.iterrows():
        best_auc = max([r[c] for c in auc_cols if not np.isnan(r[c])])
        best_feat = [c.replace('AUC_', '') for c in auc_cols if r[c] == best_auc][0]
        print(f"    {r['Mask']:25s}: excl={r['Pct_excluded']:5.1f}%, "
              f"FP_rate={r['FP_rate']:.1f}%, best_AUC={best_auc:.4f} ({best_feat})")

    # ── Summary ──
    print("\n" + "=" * 70)
    print("All outputs saved to:", OUTDIR)
    print("Files generated:")
    for f in sorted(OUTDIR.glob('sup_*')):
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name:45s} ({size_kb:.0f} KB)")
    print("=" * 70)


if __name__ == '__main__':
    main()
