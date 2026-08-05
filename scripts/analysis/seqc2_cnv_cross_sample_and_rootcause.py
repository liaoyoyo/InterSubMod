#!/usr/bin/env python3
"""
SEQC2 CNV Cross-Sample Analysis + Gain+LOH Root Cause Investigation
=====================================================================
Part A: Replicate CNV stratification across all 7 samples using Coverage_Multiple as CN proxy
Part B: Deep-dive into Gain+LOH FP root cause (HCC1395)
Part C: Evaluate whether Coverage_Multiple-based CNV zones can aid FP discrimination

Usage:
    python scripts/analysis/seqc2_cnv_cross_sample_and_rootcause.py
"""

import argparse
import json
import os
import sys
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    UNKNOWN_CURRENT_CLASS,
    V1_CURRENT_CLASSES,
    SchemaContractError,
    read_evidence,
    select_current_view,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
OUTPUT_DIR = PROJECT_ROOT / "research" / "seqc2_cnv_stratification"
FIG_DIR = OUTPUT_DIR / "figures"
DATA_DIR = OUTPUT_DIR / "data"
CANONICAL_BASE = Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical")

# Sample list — use the *_pileup_complete_matrix directories (standard runs)
SAMPLES = {
    'HCC1395': '20260314_HCC1395_paired_pileup_pileup_complete_matrix',
    'HCC1954': '20260315_HCC1954_paired_pileup_pileup_complete_matrix',
    'COLO829': '20260315_COLO829_paired_pileup_pileup_complete_matrix',
    'H1437':   '20260315_H1437_paired_pileup_pileup_complete_matrix',
    'H2009':   '20260315_H2009_paired_pileup_pileup_complete_matrix',
    'HCC1937': '20260315_HCC1937_paired_pileup_pileup_complete_matrix',
    'HCC1395_DORADO': '20260315_HCC1395_DORADO_paired_pileup_pileup_complete_matrix',
}

# Coverage_Multiple-based CN zone thresholds
# Based on HCC1395 validation: r=0.831 with SEQC2 CN
CN_THRESHOLDS = {
    'CN_Loss':    (0, 0.5),       # CN < 1 (deep deletion)
    'CN_Low':     (0.5, 0.75),    # CN ~ 1 (hemizygous deletion)
    'CN_Normal':  (0.75, 1.25),   # CN ~ 2 (diploid)
    'CN_Gain':    (1.25, 1.75),   # CN ~ 3 (single copy gain)
    'CN_HighGain':(1.75, 999),    # CN ≥ 4 (amplification)
}

CN_ZONE_ORDER = ['CN_Loss', 'CN_Low', 'CN_Normal', 'CN_Gain', 'CN_HighGain']
CN_ZONE_COLORS = {
    'CN_Loss':    '#2196F3',
    'CN_Low':     '#03A9F4',
    'CN_Normal':  '#4CAF50',
    'CN_Gain':    '#FF9800',
    'CN_HighGain':'#F44336',
}

KEY_FEATURES = [
    'NumReads', 'Coverage_Multiple', 'HP_Ratio', 'HPFineNGroups',
    'Quality_Score', 'PairwiseMeanDist', 'AlleleDelta', 'HPMergedDelta',
    'LabelAllelePermanovaF', 'LabelHPPermanovaF', 'CramersV', 'HeuristicScore',
]

plt.rcParams.update({
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 9,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.dpi': 150,
})

PROVENANCE_FIELDS = [
    'VerificationClass_V1_Deprecated', 'VerificationClass_Legacy',
    'LabelFirstSupport', 'ClusterFirstSupport', 'WithinHPSupport',
    'DispersionWarning', 'EvidencePath', 'EvidenceDerivation',
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--canonical-base', default=str(CANONICAL_BASE))
    parser.add_argument('--output-dir', default=str(OUTPUT_DIR))
    parser.add_argument(
        '--allow-unversioned-v1', action='store_true',
        help='Explicitly authorize historical unversioned v1 current taxonomy inputs.',
    )
    return parser.parse_args()


def normalize_current_taxonomy(df, allow_unversioned_v1=False, context='SEQC2 input'):
    """Return a current-class frame with explicit schema identity and direct evidence checks."""
    result = df.copy()
    if 'VerificationSchemaVersion' in result.columns:
        missing = [field for field in PROVENANCE_FIELDS if field not in result.columns]
        if missing:
            raise SchemaContractError(f"{context}: missing v2 provenance fields: {', '.join(missing)}")
        view = select_current_view(result)
        read_evidence(result)
        result['VerificationClass'] = view.values
        metadata = view.metadata()
        metadata['evidence_fields'] = list(PROVENANCE_FIELDS)
        return result, metadata

    if 'VerificationClass' not in result.columns:
        raise SchemaContractError(f'{context}: VerificationClass is missing')
    if not allow_unversioned_v1:
        raise SchemaContractError(f'{context}: unversioned current taxonomy requires --allow-unversioned-v1')

    raw = result['VerificationClass'].astype('string')
    valid = raw.isin(V1_CURRENT_CLASSES)
    unknown_counts = {
        str(key): int(value)
        for key, value in raw[~valid].fillna('<MISSING>').value_counts().sort_index().items()
    }
    result['VerificationClass'] = raw.where(valid, UNKNOWN_CURRENT_CLASS)
    message = (
        f'UNVERSIONED: {context} accepted as historical v1 current taxonomy under explicit '
        '--allow-unversioned-v1; v2 evidence provenance is unavailable'
    )
    warnings.warn(message, UserWarning, stacklevel=2)
    return result, {
        'selection_field': 'VerificationClass',
        'schema_status': 'UNVERSIONED_V1',
        'categories': list(V1_CURRENT_CLASSES) + [UNKNOWN_CURRENT_CLASS],
        'unknown_counts': unknown_counts,
        'warnings': [message],
        'evidence_fields': [],
    }


# ══════════════════════════════════════════════════════════════════════════════
# Part A: Cross-Sample CNV Stratification
# ══════════════════════════════════════════════════════════════════════════════

def load_all_samples(canonical_base=CANONICAL_BASE, allow_unversioned_v1=False):
    """Load significance_summary.csv for all samples."""
    all_data = []
    source_metadata = []
    canonical_base = Path(canonical_base)
    for sample_name, subdir in SAMPLES.items():
        tp_path = (
            canonical_base / sample_name / "paired_pileup" / subdir
            / "intersubmod_tp" / "significance_summary.csv"
        )
        fp_path = (
            canonical_base / sample_name / "paired_pileup" / subdir
            / "intersubmod_fp" / "significance_summary.csv"
        )

        if not tp_path.exists():
            print(f"  SKIP {sample_name}: TP file not found")
            continue

        tp = pd.read_csv(tp_path)
        tp['Label'] = 'TP'
        tp['Sample'] = sample_name

        if fp_path.exists():
            fp = pd.read_csv(fp_path)
            fp['Label'] = 'FP'
            fp['Sample'] = sample_name
            df = pd.concat([tp, fp], ignore_index=True)
        else:
            df = tp

        df, taxonomy = normalize_current_taxonomy(
            df,
            allow_unversioned_v1=allow_unversioned_v1,
            context=f'{sample_name} significance summaries',
        )
        taxonomy['sample'] = sample_name
        taxonomy['source_files'] = [str(tp_path)] + ([str(fp_path)] if fp_path.exists() else [])
        source_metadata.append(taxonomy)

        # Assign CN zone based on Coverage_Multiple
        if 'Coverage_Multiple' in df.columns:
            zones = []
            for cm in df['Coverage_Multiple']:
                assigned = 'Unknown'
                if pd.isna(cm):
                    assigned = 'Unknown'
                else:
                    for zone_name, (lo, hi) in CN_THRESHOLDS.items():
                        if lo <= cm < hi:
                            assigned = zone_name
                            break
                zones.append(assigned)
            df['CN_Zone'] = zones

        # Also add LOH indicator from ISM's own Potential_LOH
        if 'Potential_LOH' in df.columns:
            df['ISM_LOH'] = df['Potential_LOH'].fillna(0).astype(bool)

        all_data.append(df)
        n_tp = sum(df.Label == 'TP')
        n_fp = sum(df.Label == 'FP')
        print(f"  {sample_name}: {n_tp} TP + {n_fp} FP = {len(df)} total")

    if not all_data:
        raise SchemaContractError(f'SEQC2 inputs: no readable significance summaries below {canonical_base}')
    statuses = sorted({item['schema_status'] for item in source_metadata})
    if len(statuses) != 1:
        raise SchemaContractError(f'SEQC2 inputs mix taxonomy schema statuses: {statuses}')
    combined = pd.concat(all_data, ignore_index=True)
    combined.attrs['verification_taxonomy'] = {
        'schema_status': statuses[0],
        'selection_field': 'VerificationClass',
        'source_count': len(source_metadata),
        'sources': source_metadata,
    }
    print(f"\n  Total: {len(combined)} regions across {combined.Sample.nunique()} samples")
    return combined


def cross_sample_cn_zone_analysis(df):
    """Analyze FP rate and AUC by CN_Zone across all samples."""

    # ── Fig 8: FP rate by CN_Zone × Sample ──
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # 8a: FP rate by CN_Zone (all samples pooled)
    fp_rates = {}
    for zone in CN_ZONE_ORDER:
        zdf = df[df.CN_Zone == zone]
        n_fp = sum(zdf.Label == 'FP')
        n_total = len(zdf)
        fp_rates[zone] = n_fp / n_total if n_total > 0 else 0

    overall_fp = sum(df.Label == 'FP') / len(df)
    bars = axes[0].bar(range(len(CN_ZONE_ORDER)), [fp_rates[z] for z in CN_ZONE_ORDER],
                       color=[CN_ZONE_COLORS[z] for z in CN_ZONE_ORDER])
    axes[0].axhline(y=overall_fp, color='red', linestyle='--', alpha=0.7,
                    label=f'Overall: {overall_fp:.3f}')
    axes[0].set_xticks(range(len(CN_ZONE_ORDER)))
    axes[0].set_xticklabels(CN_ZONE_ORDER, rotation=30, ha='right')
    axes[0].set_title('FP Rate by Coverage_Multiple CN Zone\n(All Samples Pooled)')
    axes[0].set_ylabel('FP Rate')
    axes[0].legend()
    for i, z in enumerate(CN_ZONE_ORDER):
        n = sum(df.CN_Zone == z)
        axes[0].text(i, fp_rates[z] + 0.002, f'{fp_rates[z]:.3f}\nn={n}',
                     ha='center', fontsize=8)

    # 8b: FP rate by CN_Zone × Sample heatmap
    sample_order = sorted(df.Sample.unique())
    fp_matrix = pd.DataFrame(index=sample_order, columns=CN_ZONE_ORDER, dtype=float)
    count_matrix = pd.DataFrame(index=sample_order, columns=CN_ZONE_ORDER, dtype=int)

    for sample in sample_order:
        for zone in CN_ZONE_ORDER:
            zdf = df[(df.Sample == sample) & (df.CN_Zone == zone)]
            n_fp = sum(zdf.Label == 'FP')
            n_total = len(zdf)
            fp_matrix.loc[sample, zone] = n_fp / n_total if n_total > 0 else np.nan
            count_matrix.loc[sample, zone] = n_total

    # Create annotation with FP rate + count
    annot = fp_matrix.copy()
    for s in sample_order:
        for z in CN_ZONE_ORDER:
            if pd.notna(fp_matrix.loc[s, z]):
                annot.loc[s, z] = f'{fp_matrix.loc[s, z]:.3f}\n({count_matrix.loc[s, z]})'
            else:
                annot.loc[s, z] = 'N/A'

    sns.heatmap(fp_matrix.astype(float), annot=annot, fmt='', cmap='YlOrRd',
                vmin=0, vmax=0.15, ax=axes[1], linewidths=0.5,
                cbar_kws={'label': 'FP Rate'})
    axes[1].set_title('FP Rate by CN Zone × Sample\n(count in parentheses)')

    # 8c: CN_Zone distribution per sample
    zone_counts = df.groupby(['Sample', 'CN_Zone']).size().unstack(fill_value=0)
    zone_counts = zone_counts.reindex(columns=CN_ZONE_ORDER, fill_value=0)
    zone_props = zone_counts.div(zone_counts.sum(axis=1), axis=0)
    zone_props.plot(kind='barh', stacked=True, ax=axes[2],
                    color=[CN_ZONE_COLORS[z] for z in CN_ZONE_ORDER])
    axes[2].set_title('CN Zone Distribution per Sample')
    axes[2].set_xlabel('Proportion')
    axes[2].legend(fontsize=7, loc='lower right')

    plt.tight_layout()
    fig.savefig(FIG_DIR / '08_cross_sample_cn_zone_fp.png')
    plt.close()
    print("  Saved: 08_cross_sample_cn_zone_fp.png")

    return fp_matrix, count_matrix


def cross_sample_auc_analysis(df):
    """Compute AUC by CN_Zone × Sample for key features."""

    sample_order = sorted(df.Sample.unique())
    features = [f for f in KEY_FEATURES if f in df.columns]

    # ── Fig 9: AUC by CN_Zone (per-sample, then mean) ──
    # Compute per-sample AUC for each zone×feature
    results = defaultdict(lambda: defaultdict(list))

    for sample in sample_order:
        sdf = df[df.Sample == sample]
        for zone in CN_ZONE_ORDER + ['All']:
            if zone == 'All':
                zdf = sdf
            else:
                zdf = sdf[sdf.CN_Zone == zone]

            y = (zdf.Label == 'FP').astype(int).values
            if len(set(y)) < 2 or len(y) < 20:
                continue

            for feat in features:
                vals = zdf[feat].values
                mask = ~np.isnan(vals)
                if mask.sum() < 20 or len(set(y[mask])) < 2:
                    continue
                try:
                    auc = roc_auc_score(y[mask], vals[mask])
                    auc = max(auc, 1 - auc)  # direction-agnostic
                    results[(zone, feat)][sample] = auc
                except:
                    pass

    # Build mean AUC matrix
    zones_to_show = CN_ZONE_ORDER + ['All']
    auc_mean = pd.DataFrame(index=zones_to_show, columns=features, dtype=float)
    auc_std = pd.DataFrame(index=zones_to_show, columns=features, dtype=float)
    auc_n = pd.DataFrame(index=zones_to_show, columns=features, dtype=int)

    for zone in zones_to_show:
        for feat in features:
            vals = list(results[(zone, feat)].values())
            if len(vals) >= 2:
                auc_mean.loc[zone, feat] = np.mean(vals)
                auc_std.loc[zone, feat] = np.std(vals)
                auc_n.loc[zone, feat] = len(vals)
            elif len(vals) == 1:
                auc_mean.loc[zone, feat] = vals[0]
                auc_std.loc[zone, feat] = 0
                auc_n.loc[zone, feat] = 1

    # Save
    auc_mean.to_csv(DATA_DIR / 'cross_sample_cn_zone_auc_mean.tsv', sep='\t', float_format='%.4f')

    # Plot
    fig, ax = plt.subplots(figsize=(16, 6))
    mask = auc_mean.isna()

    # Annotation: mean ± std
    annot_text = auc_mean.copy()
    for z in zones_to_show:
        for f in features:
            m = auc_mean.loc[z, f]
            s = auc_std.loc[z, f]
            n = auc_n.loc[z, f]
            if pd.notna(m):
                if pd.notna(s) and s > 0:
                    annot_text.loc[z, f] = f'{m:.3f}\n±{s:.3f}'
                else:
                    annot_text.loc[z, f] = f'{m:.3f}'
            else:
                annot_text.loc[z, f] = ''

    sns.heatmap(auc_mean.astype(float), annot=annot_text, fmt='', cmap='RdYlGn',
                vmin=0.45, vmax=0.75, mask=mask, ax=ax,
                linewidths=0.5, cbar_kws={'label': 'Mean AUC (direction-agnostic)'})
    ax.set_title('Per-Sample Mean AUC by Coverage_Multiple CN Zone × Feature\n'
                 '(7 samples, direction-agnostic, ± std shown)')
    ax.set_ylabel('CN Zone')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(FIG_DIR / '09_cross_sample_auc_heatmap.png')
    plt.close()
    print("  Saved: 09_cross_sample_auc_heatmap.png")

    # ── Fig 10: Per-sample AUC dot plot for top features ──
    top_features = ['AlleleDelta', 'HPMergedDelta', 'Quality_Score',
                    'HPFineNGroups', 'LabelAllelePermanovaF']
    top_features = [f for f in top_features if f in features]

    fig, axes = plt.subplots(1, len(top_features), figsize=(4*len(top_features), 5))
    if len(top_features) == 1:
        axes = [axes]

    for i, feat in enumerate(top_features):
        ax = axes[i]
        for j, zone in enumerate(zones_to_show):
            per_sample = results[(zone, feat)]
            if per_sample:
                vals = list(per_sample.values())
                samples_list = list(per_sample.keys())
                ax.scatter([j]*len(vals), vals, alpha=0.6, s=30,
                           color=CN_ZONE_COLORS.get(zone, '#999999'))
                ax.plot([j-0.2, j+0.2], [np.mean(vals), np.mean(vals)],
                        color='black', linewidth=2)

        ax.axhline(y=0.58, color='red', linestyle='--', alpha=0.5, linewidth=0.8)
        ax.set_xticks(range(len(zones_to_show)))
        ax.set_xticklabels(zones_to_show, rotation=45, ha='right', fontsize=7)
        ax.set_title(feat)
        ax.set_ylabel('AUC' if i == 0 else '')
        ax.set_ylim(0.4, 0.9)

    plt.suptitle('Per-Sample AUC by CN Zone (dots=samples, line=mean)', fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / '10_per_sample_auc_dotplot.png')
    plt.close()
    print("  Saved: 10_per_sample_auc_dotplot.png")

    return auc_mean, results


def cross_sample_loh_interaction(df):
    """Fig 11: CN_Zone × LOH interaction analysis."""
    if 'ISM_LOH' not in df.columns:
        print("  SKIP: No ISM_LOH column")
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # 11a: FP rate by CN_Zone × LOH status
    combo_labels = []
    combo_fp_rates = []
    combo_counts = []
    combo_colors = []

    for zone in CN_ZONE_ORDER:
        for loh_status, loh_label in [(True, '+LOH'), (False, '-LOH')]:
            zdf = df[(df.CN_Zone == zone) & (df.ISM_LOH == loh_status)]
            n_fp = sum(zdf.Label == 'FP')
            n_total = len(zdf)
            fp_rate = n_fp / n_total if n_total > 0 else 0
            combo_labels.append(f'{zone}\n{loh_label}')
            combo_fp_rates.append(fp_rate)
            combo_counts.append(n_total)
            color = CN_ZONE_COLORS[zone]
            # Lighter for -LOH
            if not loh_status:
                from matplotlib.colors import to_rgba
                rgba = to_rgba(color)
                color = (rgba[0], rgba[1], rgba[2], 0.4)
            combo_colors.append(color)

    bars = axes[0].bar(range(len(combo_labels)), combo_fp_rates, color=combo_colors)
    axes[0].set_xticks(range(len(combo_labels)))
    axes[0].set_xticklabels(combo_labels, rotation=45, ha='right', fontsize=7)
    axes[0].set_title('FP Rate by CN Zone × LOH Status')
    axes[0].set_ylabel('FP Rate')
    for i, (rate, count) in enumerate(zip(combo_fp_rates, combo_counts)):
        if count > 50:
            axes[0].text(i, rate + 0.003, f'{rate:.3f}\n({count})',
                         ha='center', fontsize=6)

    # 11b: FP rate heatmap: CN_Zone × LOH × Sample
    sample_order = sorted(df.Sample.unique())
    for ax_idx, loh_val, loh_title in [(1, True, 'LOH Regions'), (2, False, 'Non-LOH Regions')]:
        matrix = pd.DataFrame(index=sample_order, columns=CN_ZONE_ORDER, dtype=float)
        for sample in sample_order:
            for zone in CN_ZONE_ORDER:
                zdf = df[(df.Sample == sample) & (df.CN_Zone == zone) & (df.ISM_LOH == loh_val)]
                n_fp = sum(zdf.Label == 'FP')
                n_total = len(zdf)
                matrix.loc[sample, zone] = n_fp / n_total if n_total > 20 else np.nan

        sns.heatmap(matrix.astype(float), annot=True, fmt='.3f', cmap='YlOrRd',
                    vmin=0, vmax=0.20, ax=axes[ax_idx], linewidths=0.5)
        axes[ax_idx].set_title(f'FP Rate: {loh_title}\n(blank if n≤20)')

    plt.tight_layout()
    fig.savefig(FIG_DIR / '11_cn_loh_interaction.png')
    plt.close()
    print("  Saved: 11_cn_loh_interaction.png")


# ══════════════════════════════════════════════════════════════════════════════
# Part B: Gain+LOH Root Cause Investigation (HCC1395)
# ══════════════════════════════════════════════════════════════════════════════

def load_hcc1395_annotated(allow_unversioned_v1=False, expected_schema_status=None):
    """Load the SEQC2-annotated HCC1395 data from Step 1."""
    path = DATA_DIR / 'annotated_hcc1395_cnv.tsv'
    if path.exists():
        df = pd.read_csv(path, sep='\t', low_memory=False)
        if 'VerificationClass' in df.columns or 'VerificationSchemaVersion' in df.columns:
            df, taxonomy = normalize_current_taxonomy(
                df,
                allow_unversioned_v1=allow_unversioned_v1,
                context='HCC1395 SEQC2 annotated table',
            )
            if expected_schema_status is not None and taxonomy['schema_status'] != expected_schema_status:
                raise SchemaContractError(
                    'SEQC2 inputs mix taxonomy schema statuses: '
                    f"canonical={expected_schema_status}, annotated={taxonomy['schema_status']}"
                )
            df.attrs['verification_taxonomy'] = taxonomy
        return df
    return None


def gain_loh_rootcause(df):
    """Deep investigation of why Gain+LOH has 10.2% FP rate."""
    gl = df[df.SEQC2_Zone == 'Gain+LOH'].copy()
    other = df[df.SEQC2_Zone != 'Gain+LOH'].copy()

    print(f"\n  Gain+LOH: {len(gl)} regions ({sum(gl.Label=='TP')} TP + {sum(gl.Label=='FP')} FP)")
    print(f"  Other zones: {len(other)} regions ({sum(other.Label=='TP')} TP + {sum(other.Label=='FP')} FP)")

    # ── Fig 12: Gain+LOH deep feature comparison ──
    fig = plt.figure(figsize=(20, 16))
    gs = gridspec.GridSpec(3, 4, hspace=0.4, wspace=0.3)

    features_deep = [
        ('NumReads', 'Read Count'),
        ('Coverage_Multiple', 'Coverage (CN proxy)'),
        ('HP_Ratio', 'HP Ratio'),
        ('AlleleDelta', 'Allele Delta'),
        ('HPMergedDelta', 'HP Merged Delta'),
        ('HPFineNGroups', 'HP Fine N Groups'),
        ('Quality_Score', 'Quality Score'),
        ('LabelAllelePermanovaF', 'Label Allele PERMANOVA F'),
        ('PairwiseMeanDist', 'Pairwise Mean Dist'),
        ('CramersV', "Cramer's V"),
        ('HeuristicScore', 'Heuristic Score'),
        ('SEQC2_CN', 'SEQC2 Copy Number'),
    ]

    for i, (feat, title) in enumerate(features_deep):
        if feat not in gl.columns:
            continue
        ax = fig.add_subplot(gs[i // 4, i % 4])

        # Compare TP vs FP within Gain+LOH
        tp_vals = gl[gl.Label == 'TP'][feat].dropna()
        fp_vals = gl[gl.Label == 'FP'][feat].dropna()

        if len(tp_vals) > 0 and len(fp_vals) > 0:
            # Clip outliers
            q01 = min(tp_vals.quantile(0.01), fp_vals.quantile(0.01))
            q99 = max(tp_vals.quantile(0.99), fp_vals.quantile(0.99))
            tp_clipped = tp_vals[(tp_vals >= q01) & (tp_vals <= q99)]
            fp_clipped = fp_vals[(fp_vals >= q01) & (fp_vals <= q99)]

            ax.hist(tp_clipped, bins=40, alpha=0.5, density=True, color='#4CAF50',
                    label=f'TP (n={len(tp_vals)})')
            ax.hist(fp_clipped, bins=40, alpha=0.5, density=True, color='#F44336',
                    label=f'FP (n={len(fp_vals)})')

            # Compute AUC
            try:
                y = np.concatenate([np.zeros(len(tp_vals)), np.ones(len(fp_vals))])
                x = np.concatenate([tp_vals.values, fp_vals.values])
                auc = roc_auc_score(y, x)
                auc = max(auc, 1 - auc)
                ax.set_title(f'{title}\nAUC={auc:.3f}', fontsize=9)
            except:
                ax.set_title(title, fontsize=9)

            ax.legend(fontsize=7)
        else:
            ax.set_title(f'{title} (insufficient data)', fontsize=9)

    fig.suptitle('Gain+LOH Zone: TP vs FP Feature Distributions (HCC1395)', fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(FIG_DIR / '12_gain_loh_rootcause_features.png')
    plt.close()
    print("  Saved: 12_gain_loh_rootcause_features.png")

    # ── Fig 13: CN-stratified analysis within Gain+LOH ──
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 13a: FP rate by CN level within Gain+LOH
    cn_bins = [2.5, 3.5, 4.5, 5.5, 6.5, 12]
    cn_labels_plot = ['CN=3', 'CN=4', 'CN=5', 'CN=6', 'CN≥7']
    gl['CN_Bin'] = pd.cut(gl.SEQC2_CN, bins=cn_bins, labels=cn_labels_plot)

    cn_fp_rates = {}
    for cn_label in cn_labels_plot:
        cdf = gl[gl.CN_Bin == cn_label]
        n_fp = sum(cdf.Label == 'FP')
        n_total = len(cdf)
        cn_fp_rates[cn_label] = n_fp / n_total if n_total > 0 else 0

    ax = axes[0, 0]
    bars = ax.bar(range(len(cn_labels_plot)), [cn_fp_rates[c] for c in cn_labels_plot],
                  color=plt.cm.Reds(np.linspace(0.3, 0.9, len(cn_labels_plot))))
    ax.set_xticks(range(len(cn_labels_plot)))
    ax.set_xticklabels(cn_labels_plot)
    ax.set_title('FP Rate by CN within Gain+LOH')
    ax.set_ylabel('FP Rate')
    for i, cn in enumerate(cn_labels_plot):
        n = sum(gl.CN_Bin == cn)
        ax.text(i, cn_fp_rates[cn] + 0.003, f'{cn_fp_rates[cn]:.3f}\n(n={n})',
                ha='center', fontsize=8)

    # 13b: HP_Ratio distribution TP vs FP in Gain+LOH
    ax = axes[0, 1]
    for label, color in [('TP', '#4CAF50'), ('FP', '#F44336')]:
        vals = gl[gl.Label == label]['HP_Ratio'].dropna()
        ax.hist(vals, bins=50, alpha=0.5, density=True, color=color,
                label=f'{label} (n={len(vals)})')
    ax.set_title('HP_Ratio: TP vs FP in Gain+LOH')
    ax.set_xlabel('HP_Ratio')
    ax.legend()

    # 13c: VerificationClass distribution
    ax = axes[0, 2]
    if 'VerificationClass' in gl.columns:
        vc_tp = gl[gl.Label == 'TP']['VerificationClass'].value_counts(normalize=True)
        vc_fp = gl[gl.Label == 'FP']['VerificationClass'].value_counts(normalize=True)
        vc_all = sorted(set(vc_tp.index) | set(vc_fp.index))
        x_pos = np.arange(len(vc_all))
        width = 0.35
        ax.bar(x_pos - width/2, [vc_tp.get(v, 0) for v in vc_all], width,
               label='TP', color='#4CAF50')
        ax.bar(x_pos + width/2, [vc_fp.get(v, 0) for v in vc_all], width,
               label='FP', color='#F44336')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(vc_all, rotation=30)
        ax.set_title('VerificationClass in Gain+LOH')
        ax.legend()

    # 13d: Quality_Score TP vs FP
    ax = axes[1, 0]
    for label, color in [('TP', '#4CAF50'), ('FP', '#F44336')]:
        vals = gl[gl.Label == label]['Quality_Score'].dropna()
        ax.hist(vals, bins=40, alpha=0.5, density=True, color=color,
                label=f'{label} (n={len(vals)})')
    ax.set_title('Quality_Score: TP vs FP in Gain+LOH')
    ax.set_xlabel('Quality_Score')
    ax.legend()

    # 13e: AlleleDelta TP vs FP
    ax = axes[1, 1]
    for label, color in [('TP', '#4CAF50'), ('FP', '#F44336')]:
        vals = gl[gl.Label == label]['AlleleDelta'].dropna()
        ax.hist(vals, bins=40, alpha=0.5, density=True, color=color,
                label=f'{label} (n={len(vals)})')
    ax.set_title('AlleleDelta: TP vs FP in Gain+LOH')
    ax.set_xlabel('AlleleDelta')
    ax.legend()

    # 13f: HPMergedDelta TP vs FP
    ax = axes[1, 2]
    for label, color in [('TP', '#4CAF50'), ('FP', '#F44336')]:
        vals = gl[gl.Label == label]['HPMergedDelta'].dropna()
        ax.hist(vals, bins=40, alpha=0.5, density=True, color=color,
                label=f'{label} (n={len(vals)})')
    ax.set_title('HPMergedDelta: TP vs FP in Gain+LOH')
    ax.set_xlabel('HPMergedDelta')
    ax.legend()

    plt.suptitle('Gain+LOH Root Cause: CN Stratification & Feature Analysis', fontsize=13, y=1.01)
    plt.tight_layout()
    fig.savefig(FIG_DIR / '13_gain_loh_cn_stratification.png')
    plt.close()
    print("  Saved: 13_gain_loh_cn_stratification.png")

    return cn_fp_rates


def gain_loh_vs_others_comparison(hcc_df):
    """Fig 14: Compare Gain+LOH FP characteristics vs other zones."""
    zones_compare = ['Neutral', 'Gain-only', 'LOH-only', 'Gain+LOH']

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # For each zone, compare FP characteristics
    features_compare = [
        ('Coverage_Multiple', 'Coverage_Multiple'),
        ('HP_Ratio', 'HP_Ratio'),
        ('AlleleDelta', 'AlleleDelta'),
        ('HPMergedDelta', 'HPMergedDelta'),
        ('Quality_Score', 'Quality_Score'),
        ('HPFineNGroups', 'HPFineNGroups'),
    ]

    for i, (feat, title) in enumerate(features_compare):
        ax = axes[i // 3, i % 3]
        if feat not in hcc_df.columns:
            continue

        # Box plot: FP only, across zones
        fp_data = hcc_df[hcc_df.Label == 'FP']
        plot_data = []
        for zone in zones_compare:
            vals = fp_data[fp_data.SEQC2_Zone == zone][feat].dropna()
            for v in vals:
                plot_data.append({'Zone': zone, feat: v})

        if plot_data:
            plot_df = pd.DataFrame(plot_data)
            # Clip extremes
            q01, q99 = plot_df[feat].quantile([0.01, 0.99])
            plot_df = plot_df[(plot_df[feat] >= q01) & (plot_df[feat] <= q99)]
            sns.violinplot(data=plot_df, x='Zone', y=feat, ax=ax,
                           palette={z: '#F44336' for z in zones_compare},
                           inner='quartile', cut=0, alpha=0.7)
            ax.set_title(f'FP {title} by Zone')
            ax.tick_params(axis='x', rotation=30)

    plt.suptitle('FP Characteristics Across SEQC2 Zones (HCC1395)', fontsize=13, y=1.01)
    plt.tight_layout()
    fig.savefig(FIG_DIR / '14_fp_characteristics_by_zone.png')
    plt.close()
    print("  Saved: 14_fp_characteristics_by_zone.png")


# ══════════════════════════════════════════════════════════════════════════════
# Part C: Feasibility of Coverage_Multiple-based CNV Zone for FP discrimination
# ══════════════════════════════════════════════════════════════════════════════

def feasibility_analysis(df, verification_metadata=None):
    """Fig 15: Can Coverage_Multiple-based zones improve FP discrimination?"""

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 15a: Zone-conditional AUC comparison across samples
    # For each sample, compute best AUC in CN_Normal vs All
    sample_order = sorted(df.Sample.unique())
    features = [f for f in KEY_FEATURES if f in df.columns]

    normal_best = {}
    all_best = {}
    for sample in sample_order:
        sdf = df[df.Sample == sample]
        # All
        y_all = (sdf.Label == 'FP').astype(int).values
        best_all = 0.5
        if len(set(y_all)) >= 2:
            for feat in features:
                vals = sdf[feat].values
                mask = ~np.isnan(vals)
                if mask.sum() >= 20 and len(set(y_all[mask])) >= 2:
                    try:
                        auc = roc_auc_score(y_all[mask], vals[mask])
                        best_all = max(best_all, max(auc, 1-auc))
                    except:
                        pass
        all_best[sample] = best_all

        # CN_Normal only
        ndf = sdf[sdf.CN_Zone == 'CN_Normal']
        y_n = (ndf.Label == 'FP').astype(int).values
        best_n = 0.5
        if len(set(y_n)) >= 2 and len(y_n) >= 30:
            for feat in features:
                vals = ndf[feat].values
                mask = ~np.isnan(vals)
                if mask.sum() >= 20 and len(set(y_n[mask])) >= 2:
                    try:
                        auc = roc_auc_score(y_n[mask], vals[mask])
                        best_n = max(best_n, max(auc, 1-auc))
                    except:
                        pass
        normal_best[sample] = best_n

    ax = axes[0, 0]
    x = np.arange(len(sample_order))
    width = 0.35
    ax.bar(x - width/2, [all_best[s] for s in sample_order], width,
           label='All Zones', color='#2196F3')
    ax.bar(x + width/2, [normal_best[s] for s in sample_order], width,
           label='CN_Normal Only', color='#4CAF50')
    ax.axhline(y=0.58, color='red', linestyle='--', alpha=0.5, label='0.58 ceiling')
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=45, ha='right', fontsize=7)
    ax.set_title('Best Feature AUC: All vs CN_Normal')
    ax.set_ylabel('Best AUC')
    ax.legend(fontsize=7)

    # 15b: FP proportion contributed by each CN zone
    ax = axes[0, 1]
    fp_by_zone = df[df.Label == 'FP'].groupby('CN_Zone').size()
    fp_by_zone = fp_by_zone.reindex(CN_ZONE_ORDER, fill_value=0)
    total_fp = fp_by_zone.sum()
    ax.pie(fp_by_zone.values, labels=[f'{z}\n{v} ({v/total_fp:.1%})'
           for z, v in zip(fp_by_zone.index, fp_by_zone.values)],
           colors=[CN_ZONE_COLORS[z] for z in fp_by_zone.index],
           autopct='', startangle=140)
    ax.set_title(f'FP Distribution by CN Zone\n(Total FP={total_fp})')

    # 15c: If we exclude CN_HighGain+LOH regions, how much FP do we remove vs TP loss?
    ax = axes[0, 2]
    # Simulate different exclusion strategies
    strategies = []
    # Strategy 1: Exclude CN_HighGain
    for zone_to_exclude in CN_ZONE_ORDER:
        excluded = df[df.CN_Zone == zone_to_exclude]
        remaining = df[df.CN_Zone != zone_to_exclude]
        fp_removed = sum(excluded.Label == 'FP')
        tp_lost = sum(excluded.Label == 'TP')
        strategies.append({
            'Strategy': f'Exclude {zone_to_exclude}',
            'FP_Removed': fp_removed,
            'TP_Lost': tp_lost,
            'FP_Removal%': fp_removed / sum(df.Label == 'FP') * 100,
            'TP_Loss%': tp_lost / sum(df.Label == 'TP') * 100,
        })

    # Combined: exclude CN_HighGain + CN_Loss
    excluded_combo = df[df.CN_Zone.isin(['CN_HighGain', 'CN_Loss'])]
    strategies.append({
        'Strategy': 'Exclude HighGain+Loss',
        'FP_Removed': sum(excluded_combo.Label == 'FP'),
        'TP_Lost': sum(excluded_combo.Label == 'TP'),
        'FP_Removal%': sum(excluded_combo.Label == 'FP') / sum(df.Label == 'FP') * 100,
        'TP_Loss%': sum(excluded_combo.Label == 'TP') / sum(df.Label == 'TP') * 100,
    })

    strat_df = pd.DataFrame(strategies)
    ax.scatter(strat_df['TP_Loss%'], strat_df['FP_Removal%'], s=100, c='#2196F3')
    for _, row in strat_df.iterrows():
        ax.annotate(row['Strategy'].replace('Exclude ', ''),
                    (row['TP_Loss%'], row['FP_Removal%']),
                    fontsize=7, ha='left', va='bottom',
                    xytext=(5, 5), textcoords='offset points')
    ax.set_xlabel('TP Loss (%)')
    ax.set_ylabel('FP Removal (%)')
    ax.set_title('Zone Exclusion Trade-off')
    ax.axline((0,0), slope=1, color='gray', linestyle='--', alpha=0.3, label='Break-even')
    ax.legend(fontsize=8)

    # 15d: Zone-aware threshold simulation
    ax = axes[1, 0]
    # For the best feature in each zone, compute threshold at various FP removal rates
    # Focus on AlleleDelta (strongest in Gain+LOH)
    if 'AlleleDelta' in df.columns:
        for zone in ['CN_Normal', 'CN_Gain', 'CN_HighGain', 'All']:
            if zone == 'All':
                zdf = df
            else:
                zdf = df[df.CN_Zone == zone]

            y = (zdf.Label == 'FP').astype(int).values
            vals = zdf['AlleleDelta'].values
            mask = ~np.isnan(vals)
            if mask.sum() < 50 or len(set(y[mask])) < 2:
                continue

            # ROC curve
            from sklearn.metrics import roc_curve
            fpr, tpr, _ = roc_curve(y[mask], vals[mask])
            ax.plot(fpr, tpr, label=f'{zone} (n={mask.sum()})', alpha=0.8)

        ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title('AlleleDelta ROC by CN Zone')
        ax.legend(fontsize=7)

    # 15e: Feature importance shift by zone
    ax = axes[1, 1]
    # Top-3 features per zone
    zone_top = {}
    for zone in CN_ZONE_ORDER:
        zdf = df[df.CN_Zone == zone]
        y = (zdf.Label == 'FP').astype(int).values
        if len(set(y)) < 2 or len(y) < 50:
            continue
        zone_aucs = {}
        for feat in features:
            vals = zdf[feat].values
            mask = ~np.isnan(vals)
            if mask.sum() < 30 and len(set(y[mask])) >= 2:
                continue
            try:
                auc = roc_auc_score(y[mask], vals[mask])
                zone_aucs[feat] = max(auc, 1-auc)
            except:
                pass
        if zone_aucs:
            top3 = sorted(zone_aucs.items(), key=lambda x: x[1], reverse=True)[:3]
            zone_top[zone] = top3

    # Plot as grouped bar
    if zone_top:
        text_lines = []
        for zone, top3 in zone_top.items():
            line = f"{zone}: " + ", ".join([f"{f}={a:.3f}" for f, a in top3])
            text_lines.append(line)
        ax.text(0.05, 0.95, '\n'.join(text_lines), transform=ax.transAxes,
                fontsize=8, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax.set_title('Top-3 Features per CN Zone')
        ax.axis('off')

    # 15f: Summary statistics table
    ax = axes[1, 2]
    summary = []
    for zone in CN_ZONE_ORDER:
        zdf = df[df.CN_Zone == zone]
        n_tp = sum(zdf.Label == 'TP')
        n_fp = sum(zdf.Label == 'FP')
        fp_rate = n_fp / (n_tp + n_fp) if (n_tp + n_fp) > 0 else 0
        n_samples = zdf.Sample.nunique()
        summary.append([zone, n_tp + n_fp, n_fp, f'{fp_rate:.4f}', n_samples])

    ax.axis('off')
    table = ax.table(cellText=summary,
                     colLabels=['CN Zone', 'Total', 'FP', 'FP Rate', 'Samples'],
                     loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    ax.set_title('Summary: CN Zone Statistics (7 Samples)')

    plt.suptitle('Feasibility: Coverage_Multiple-based CNV Zones for FP Discrimination',
                 fontsize=13, y=1.01)
    plt.tight_layout()
    fig.savefig(FIG_DIR / '15_feasibility_analysis.png')
    plt.close()
    print("  Saved: 15_feasibility_analysis.png")

    if verification_metadata is not None:
        strat_df['VerificationSchemaStatus'] = verification_metadata['schema_status']
        strat_df['VerificationSelectionField'] = verification_metadata['selection_field']
    strat_df.to_csv(DATA_DIR / 'zone_exclusion_tradeoff.tsv', sep='\t', index=False)
    return strat_df


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    global OUTPUT_DIR, FIG_DIR, DATA_DIR, CANONICAL_BASE
    args = parse_args()
    OUTPUT_DIR = Path(args.output_dir).resolve()
    FIG_DIR = OUTPUT_DIR / 'figures'
    DATA_DIR = OUTPUT_DIR / 'data'
    CANONICAL_BASE = Path(args.canonical_base).resolve()

    print("=" * 70)
    print("SEQC2 CNV Cross-Sample Analysis + Root Cause Investigation")
    print("=" * 70)

    # ══════════════════════════════════════════════════════════════════════════
    # Part A: Cross-Sample Analysis
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Part A] Loading all samples...")
    try:
        all_df = load_all_samples(
            CANONICAL_BASE,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
        verification_metadata = all_df.attrs['verification_taxonomy']
        hcc_df = load_hcc1395_annotated(
            allow_unversioned_v1=args.allow_unversioned_v1,
            expected_schema_status=verification_metadata['schema_status'],
        )
    except SchemaContractError as exc:
        raise SystemExit(f'[seqc2-cnv][schema-contract] {exc}') from exc

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(DATA_DIR / 'verification_taxonomy.json', 'w') as handle:
        json.dump(verification_metadata, handle, indent=2)

    print("\n[Part A] CN Zone × Sample FP rate analysis...")
    fp_matrix, count_matrix = cross_sample_cn_zone_analysis(all_df)

    print("\n  FP Rate Matrix (per-sample):")
    print(fp_matrix.to_string())
    print("\n  Count Matrix:")
    print(count_matrix.to_string())

    print("\n[Part A] Per-sample AUC by CN Zone...")
    auc_mean, auc_results = cross_sample_auc_analysis(all_df)

    print("\n  Mean AUC (key features):")
    key_cols = ['AlleleDelta', 'HPMergedDelta', 'Quality_Score', 'HPFineNGroups']
    key_cols = [c for c in key_cols if c in auc_mean.columns]
    print(auc_mean[key_cols].to_string())

    print("\n[Part A] CN Zone × LOH interaction...")
    cross_sample_loh_interaction(all_df)

    # ══════════════════════════════════════════════════════════════════════════
    # Part B: Gain+LOH Root Cause
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Part B] Loading HCC1395 SEQC2-annotated data...")
    if hcc_df is not None:
        print("\n[Part B] Gain+LOH root cause analysis...")
        cn_fp_rates = gain_loh_rootcause(hcc_df)

        print("\n  FP rate by CN within Gain+LOH:")
        for cn, rate in cn_fp_rates.items():
            print(f"    {cn}: {rate:.4f}")

        print("\n[Part B] FP characteristics comparison...")
        gain_loh_vs_others_comparison(hcc_df)

    # ══════════════════════════════════════════════════════════════════════════
    # Part C: Feasibility Analysis
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Part C] Feasibility of Coverage_Multiple-based zone discrimination...")
    strat_df = feasibility_analysis(all_df, verification_metadata)

    print("\n  Zone exclusion trade-off:")
    print(strat_df.to_string(index=False))

    # ══════════════════════════════════════════════════════════════════════════
    # Summary
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nFigures: {len(list(FIG_DIR.glob('*.png')))} PNG files")
    print(f"Data: {len(list(DATA_DIR.glob('*.tsv')))} TSV files")

    # Key conclusions
    print("\n── Key Conclusions ──")
    print("1. FP rate pattern (CN_HighGain > CN_Gain > CN_Normal > CN_Loss): "
          "consistent across samples?")
    fp_highgain = fp_matrix['CN_HighGain'].dropna()
    fp_normal = fp_matrix['CN_Normal'].dropna()
    if len(fp_highgain) > 0 and len(fp_normal) > 0:
        n_higher = sum(fp_highgain > fp_normal)
        print(f"   CN_HighGain > CN_Normal in {n_higher}/{len(fp_highgain)} samples")

    print("\n2. Best feature shifts by CN zone:")
    for zone in CN_ZONE_ORDER:
        if zone in auc_mean.index:
            best = auc_mean.loc[zone].idxmax()
            val = auc_mean.loc[zone].max()
            if pd.notna(val):
                print(f"   {zone}: {best} = {val:.3f}")

    # Simpson's Paradox check
    print("\n3. Simpson's Paradox check:")
    for feat in ['AlleleDelta', 'HPMergedDelta', 'Quality_Score']:
        if feat in auc_mean.columns:
            all_auc = auc_mean.loc['All', feat] if 'All' in auc_mean.index else np.nan
            zone_aucs = [auc_mean.loc[z, feat] for z in CN_ZONE_ORDER
                         if z in auc_mean.index and pd.notna(auc_mean.loc[z, feat])]
            if zone_aucs and pd.notna(all_auc):
                mean_zone = np.mean(zone_aucs)
                print(f"   {feat}: All={all_auc:.3f}, mean(zones)={mean_zone:.3f}, "
                      f"diff={mean_zone - all_auc:+.3f}")


if __name__ == '__main__':
    main()
