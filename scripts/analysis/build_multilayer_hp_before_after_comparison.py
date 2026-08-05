#!/usr/bin/env python3
"""
Multilayer HP Before/After Comparison
Compare old (before multi-layer HP) vs new (after) benchmark output.
Generates figures and summary data for visual inspection.
"""

import argparse
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    select_current_view,
    select_legacy_view,
)

# ============================================================================
# Configuration
# ============================================================================

BASE = Path("/big7_disk/liaoyoyo2001/InterSubMod/output/multilayer_hp_benchmark")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/output/multilayer_hp_benchmark/comparison_figures")

BEFORE_TP = BASE / "before_paired_tp" / "significance_summary.csv"
BEFORE_FP = BASE / "before_paired_fp" / "significance_summary.csv"
AFTER_TP  = BASE / "paired_tp" / "significance_summary.csv"
AFTER_FP  = BASE / "paired_fp" / "significance_summary.csv"

# ============================================================================
# Load and merge data
# ============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verification-view",
        required=True,
        choices=("current-v2", "legacy"),
        help=(
            "Select one semantic class view for both matrix axes. current-v2 rejects "
            "v1 input; legacy reads VerificationClass_Legacy or an explicitly allowed H1 field."
        ),
    )
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Authorize unversioned VerificationClass only when --verification-view=legacy",
    )
    return parser.parse_args()


def attach_verification_view(
    frame: pd.DataFrame,
    verification_view: str,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    """Attach one explicitly selected semantic view; never mix v1 raw with v2 current."""
    if verification_view == "current-v2":
        if allow_unversioned_v1:
            raise SchemaContractError(
                "--allow-unversioned-v1 is only valid with --verification-view=legacy"
            )
        view = select_current_view(frame)
        label = "Current VerificationClass v2"
    elif verification_view == "legacy":
        view = select_legacy_view(
            frame,
            allow_unversioned_v1=allow_unversioned_v1,
            unversioned_unknown_policy="fail",
        )
        label = "Legacy VerificationClass"
    else:
        raise SchemaContractError(f"unknown verification view: {verification_view!r}")

    prepared = frame.copy()
    prepared["_verification_class_selected"] = view.values
    prepared.attrs["verification_schema_contract"] = {
        "requested_view": verification_view,
        "view_label": label,
        "selection_field": view.field,
        "schema_status": view.schema_status,
        "categories": list(view.categories),
        "unknown_counts": dict(view.unknown_counts),
        "allow_unversioned_v1": allow_unversioned_v1,
        "warnings": list(view.warning_messages),
    }
    return prepared


def load_and_merge(
    before_path,
    after_path,
    label,
    verification_view,
    allow_unversioned_v1=False,
):
    """Load before/after CSVs and merge on variant key."""
    df_before = attach_verification_view(
        pd.read_csv(before_path),
        verification_view,
        allow_unversioned_v1=allow_unversioned_v1,
    )
    df_after = attach_verification_view(
        pd.read_csv(after_path),
        verification_view,
        allow_unversioned_v1=allow_unversioned_v1,
    )
    before_contract = dict(df_before.attrs["verification_schema_contract"])
    after_contract = dict(df_after.attrs["verification_schema_contract"])
    if before_contract["categories"] != after_contract["categories"]:
        raise SchemaContractError(
            "before/after verification views expose different category contracts"
        )

    # Create unique key
    for df in [df_before, df_after]:
        df['key'] = df['Chr'] + ':' + df['Pos'].astype(str) + ':' + df['Ref'] + ':' + df['Alt']

    merged = pd.merge(df_before, df_after, on='key', suffixes=('_before', '_after'))
    print(f"[{label}] Before: {len(df_before)}, After: {len(df_after)}, Merged: {len(merged)}")
    contract = {
        "label": label,
        "before_path": str(before_path),
        "after_path": str(after_path),
        "before": before_contract,
        "after": after_contract,
        "comparison_semantics": before_contract["view_label"],
        "direct_v1_v2_current_mix": False,
    }
    return merged, df_before, df_after, contract


# ============================================================================
# Figure 1: PassedGating transition matrix
# ============================================================================

def fig1_gating_transition(merged_tp, merged_fp):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, (df, label) in zip(axes, [(merged_tp, 'TP'), (merged_fp, 'FP')]):
        # Build transition matrix
        cats = [False, True]
        mat = np.zeros((2, 2), dtype=int)
        for i, bv in enumerate(cats):
            for j, av in enumerate(cats):
                mat[i, j] = ((df['PassedGating_before'] == bv) & (df['PassedGating_after'] == av)).sum()

        im = ax.imshow(mat, cmap='Blues', aspect='auto')
        ax.set_xticks([0, 1]); ax.set_xticklabels(['False', 'True'])
        ax.set_yticks([0, 1]); ax.set_yticklabels(['False', 'True'])
        ax.set_xlabel('After (New)')
        ax.set_ylabel('Before (Old)')
        ax.set_title(f'PassedGating Transition — {label} (n={len(df)})')

        # Annotate cells
        for i in range(2):
            for j in range(2):
                pct = mat[i, j] / len(df) * 100
                txt = f"{mat[i, j]}\n({pct:.1f}%)"
                color = 'white' if mat[i, j] > len(df) * 0.3 else 'black'
                ax.text(j, i, txt, ha='center', va='center', fontsize=12, fontweight='bold', color=color)

        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle('Figure 1: PassedGating Transition (Before → After)', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig1_gating_transition.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("[Fig1] Saved gating transition matrix")


# ============================================================================
# Figure 2: VerificationClass transition (Sankey-style table)
# ============================================================================

def fig2_verification_class_transition(merged_tp, merged_fp, classes, view_label):
    width = 18 if len(classes) > 6 else 14
    fig, axes = plt.subplots(1, 2, figsize=(width, 7))

    for ax, (df, label) in zip(axes, [(merged_tp, 'TP'), (merged_fp, 'FP')]):
        mat = np.zeros((len(classes), len(classes)), dtype=int)
        for i, bv in enumerate(classes):
            for j, av in enumerate(classes):
                mat[i, j] = (
                    (df['_verification_class_selected_before'] == bv)
                    & (df['_verification_class_selected_after'] == av)
                ).sum()

        im = ax.imshow(mat, cmap='YlOrRd', aspect='auto')
        ax.set_xticks(range(len(classes)))
        ax.set_xticklabels(classes, rotation=45, ha='right')
        ax.set_yticks(range(len(classes)))
        ax.set_yticklabels(classes)
        ax.set_xlabel('After (New)')
        ax.set_ylabel('Before (Old)')
        ax.set_title(f'{view_label} Transition — {label}')

        for i in range(len(classes)):
            for j in range(len(classes)):
                if mat[i, j] > 0:
                    pct = mat[i, j] / len(df) * 100
                    txt = f"{mat[i, j]}\n({pct:.1f}%)"
                    color = 'white' if mat[i, j] > len(df) * 0.15 else 'black'
                    ax.text(j, i, txt, ha='center', va='center', fontsize=9, fontweight='bold', color=color)

        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle(
        f'Figure 2: {view_label} Transition (Before → After)',
        fontsize=14,
        fontweight='bold',
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig2_verification_class_transition.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("[Fig2] Saved VerificationClass transition matrix")


# ============================================================================
# Figure 3: GlobalP scatter + delta histogram
# ============================================================================

def fig3_globalp_comparison(merged_tp, merged_fp):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for col_idx, (df, label) in enumerate([(merged_tp, 'TP'), (merged_fp, 'FP')]):
        # Scatter: -log10(p) before vs after
        ax = axes[0, col_idx]
        p_before = df['GlobalP_before'].clip(lower=1e-20)
        p_after  = df['GlobalP_after'].clip(lower=1e-20)
        logp_b = -np.log10(p_before)
        logp_a = -np.log10(p_after)

        ax.scatter(logp_b, logp_a, alpha=0.1, s=2, c='steelblue')
        max_val = max(logp_b.max(), logp_a.max())
        ax.plot([0, max_val], [0, max_val], 'r--', lw=1, label='y=x')
        ax.set_xlabel('-log10(GlobalP) Before')
        ax.set_ylabel('-log10(GlobalP) After')
        ax.set_title(f'GlobalP Scatter — {label}')
        ax.legend()

        # Count changes
        n_improved = (logp_a > logp_b + 0.01).sum()
        n_same = ((logp_a - logp_b).abs() <= 0.01).sum()
        n_worse = (logp_a < logp_b - 0.01).sum()
        ax.text(0.05, 0.95, f"Improved: {n_improved}\nSame: {n_same}\nWorse: {n_worse}",
                transform=ax.transAxes, va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        # Delta histogram
        ax = axes[1, col_idx]
        delta = logp_a - logp_b
        ax.hist(delta, bins=100, color='steelblue', alpha=0.7, edgecolor='none')
        ax.axvline(0, color='red', ls='--', lw=1)
        ax.set_xlabel('Delta -log10(GlobalP) (After - Before)')
        ax.set_ylabel('Count')
        ax.set_title(f'GlobalP Delta Distribution — {label}')
        ax.text(0.95, 0.95, f"mean={delta.mean():.4f}\nmedian={delta.median():.4f}",
                transform=ax.transAxes, va='top', ha='right', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    fig.suptitle('Figure 3: GlobalP Before vs After', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig3_globalp_comparison.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("[Fig3] Saved GlobalP comparison")


# ============================================================================
# Figure 4: HeuristicScore comparison
# ============================================================================

def fig4_heuristic_score(merged_tp, merged_fp):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for col_idx, (df, label) in enumerate([(merged_tp, 'TP'), (merged_fp, 'FP')]):
        # Scatter
        ax = axes[0, col_idx]
        hs_b = df['HeuristicScore_before']
        hs_a = df['HeuristicScore_after']
        ax.scatter(hs_b, hs_a, alpha=0.1, s=2, c='darkorange')
        max_val = max(hs_b.max(), hs_a.max())
        ax.plot([0, max_val], [0, max_val], 'r--', lw=1)
        ax.set_xlabel('HeuristicScore Before')
        ax.set_ylabel('HeuristicScore After')
        ax.set_title(f'HeuristicScore Scatter — {label}')

        n_improved = (hs_a > hs_b + 0.01).sum()
        n_same = ((hs_a - hs_b).abs() <= 0.01).sum()
        ax.text(0.05, 0.95, f"Improved: {n_improved}\nSame: {n_same}",
                transform=ax.transAxes, va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        # Delta histogram
        ax = axes[1, col_idx]
        delta = hs_a - hs_b
        ax.hist(delta, bins=100, color='darkorange', alpha=0.7, edgecolor='none')
        ax.axvline(0, color='red', ls='--', lw=1)
        ax.set_xlabel('Delta HeuristicScore (After - Before)')
        ax.set_ylabel('Count')
        ax.set_title(f'HeuristicScore Delta — {label}')
        ax.text(0.95, 0.95, f"mean={delta.mean():.4f}\nmedian={delta.median():.4f}",
                transform=ax.transAxes, va='top', ha='right', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    fig.suptitle('Figure 4: HeuristicScore Before vs After', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig4_heuristic_score.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("[Fig4] Saved HeuristicScore comparison")


# ============================================================================
# Figure 5: Quality Score + Tier comparison
# ============================================================================

def fig5_quality(merged_tp, merged_fp):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for col_idx, (df, label) in enumerate([(merged_tp, 'TP'), (merged_fp, 'FP')]):
        # Quality_Score scatter
        ax = axes[0, col_idx]
        qs_b = df['Quality_Score_before']
        qs_a = df['Quality_Score_after']
        ax.scatter(qs_b, qs_a, alpha=0.1, s=2, c='forestgreen')
        max_val = max(qs_b.max(), qs_a.max())
        ax.plot([0, max_val], [0, max_val], 'r--', lw=1)
        ax.set_xlabel('Quality_Score Before')
        ax.set_ylabel('Quality_Score After')
        ax.set_title(f'Quality_Score Scatter — {label}')

        n_improved = (qs_a > qs_b + 0.01).sum()
        n_same = ((qs_a - qs_b).abs() <= 0.01).sum()
        n_worse = (qs_a < qs_b - 0.01).sum()
        ax.text(0.05, 0.95, f"Improved: {n_improved}\nSame: {n_same}\nWorse: {n_worse}",
                transform=ax.transAxes, va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        # Quality_Tier bar comparison
        ax = axes[1, col_idx]
        tiers = ['Excellent', 'Good', 'Fair', 'Poor']
        before_counts = df['Quality_Tier_before'].value_counts()
        after_counts  = df['Quality_Tier_after'].value_counts()

        x = np.arange(len(tiers))
        w = 0.35
        b_vals = [before_counts.get(t, 0) for t in tiers]
        a_vals = [after_counts.get(t, 0) for t in tiers]

        bars1 = ax.bar(x - w/2, b_vals, w, label='Before', color='lightcoral', edgecolor='gray')
        bars2 = ax.bar(x + w/2, a_vals, w, label='After', color='lightseagreen', edgecolor='gray')
        ax.set_xticks(x); ax.set_xticklabels(tiers)
        ax.set_ylabel('Count')
        ax.set_title(f'Quality_Tier Distribution — {label}')
        ax.legend()

        # Annotate differences
        for i, (bv, av) in enumerate(zip(b_vals, a_vals)):
            diff = av - bv
            if diff != 0:
                sign = '+' if diff > 0 else ''
                ax.text(i, max(bv, av) + len(df)*0.01, f'{sign}{diff}',
                       ha='center', fontsize=9, fontweight='bold', color='red' if diff < 0 else 'green')

    fig.suptitle('Figure 5: Quality Score & Tier Before vs After', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig5_quality_comparison.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("[Fig5] Saved Quality comparison")


# ============================================================================
# Figure 6: New columns distribution (HP Family & HP Fine)
# ============================================================================

def fig6_new_columns(df_after_tp, df_after_fp):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for row, (df, label) in enumerate([(df_after_tp, 'TP'), (df_after_fp, 'FP')]):
        # GlobalP_HPFamily histogram
        ax = axes[row, 0]
        p_vals = df['GlobalP_HPFamily'].clip(lower=1e-20)
        logp = -np.log10(p_vals)
        ax.hist(logp, bins=80, color='mediumpurple', alpha=0.7, edgecolor='none')
        ax.axvline(-np.log10(0.05), color='red', ls='--', lw=1, label='p=0.05')
        ax.set_xlabel('-log10(GlobalP_HPFamily)')
        ax.set_ylabel('Count')
        ax.set_title(f'HP Family p-value — {label}')
        n_sig = (df['GlobalP_HPFamily'] < 0.05).sum()
        ax.text(0.95, 0.95, f"Sig: {n_sig} ({n_sig/len(df)*100:.1f}%)",
                transform=ax.transAxes, va='top', ha='right', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax.legend()

        # CramersV_HPFamily histogram
        ax = axes[row, 1]
        ax.hist(df['CramersV_HPFamily'], bins=80, color='mediumpurple', alpha=0.7, edgecolor='none')
        ax.set_xlabel("Cramer's V (HP Family)")
        ax.set_ylabel('Count')
        ax.set_title(f"HP Family Cramer's V — {label}")

        # HPFine_NGroups_CF distribution
        ax = axes[row, 2]
        ngroups = df['HPFine_NGroups_CF'].value_counts().sort_index()
        bars = ax.bar(ngroups.index, ngroups.values, color='teal', alpha=0.7, edgecolor='gray')
        ax.set_xlabel('HPFine_NGroups_CF')
        ax.set_ylabel('Count')
        ax.set_title(f'HP Fine N Groups (Cluster-First) — {label}')
        for bar, val in zip(bars, ngroups.values):
            pct = val / len(df) * 100
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                   f'{pct:.1f}%', ha='center', va='bottom', fontsize=10)

    fig.suptitle('Figure 6: New Multi-Layer HP Columns', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig6_new_columns.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("[Fig6] Saved new columns distribution")


# ============================================================================
# Figure 7: Sanity check — unchanged columns
# ============================================================================

def fig7_sanity_check(merged_tp):
    """Verify columns that should NOT change are identical."""
    check_cols_before = ['NumReads', 'NumCpGs', 'ClusterPermanovaF', 'ClusterPermanovaP',
                         'HPMergedP', 'AlleleP', 'HP1FamilyN', 'HP2FamilyN',
                         'LabelHPPermanovaF', 'LabelAllelePermanovaF']

    results = []
    for col in check_cols_before:
        cb = f'{col}_before'
        ca = f'{col}_after'
        if cb in merged_tp.columns and ca in merged_tp.columns:
            # Compare with tolerance for float columns
            b = merged_tp[cb].fillna(-999)
            a = merged_tp[ca].fillna(-999)
            if b.dtype in [np.float64, np.float32]:
                n_diff = (np.abs(b - a) > 1e-6).sum()
            else:
                n_diff = (b != a).sum()
            results.append({'Column': col, 'N_Different': n_diff, 'Pct': n_diff / len(merged_tp) * 100})

    df_check = pd.DataFrame(results)

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ['forestgreen' if x == 0 else 'red' for x in df_check['N_Different']]
    bars = ax.barh(df_check['Column'], df_check['N_Different'], color=colors, alpha=0.7)
    ax.set_xlabel('Number of Different Values')
    ax.set_title('Figure 7: Sanity Check — Columns That Should NOT Change (TP)')

    for bar, val, pct in zip(bars, df_check['N_Different'], df_check['Pct']):
        txt = f'{val} ({pct:.2f}%)' if val > 0 else '0 (identical)'
        ax.text(bar.get_width() + len(merged_tp)*0.001, bar.get_y() + bar.get_height()/2,
               txt, va='center', fontsize=10)

    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig7_sanity_check.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("[Fig7] Saved sanity check")
    return df_check


# ============================================================================
# Summary statistics
# ============================================================================

def compute_summary(
    merged_tp,
    merged_fp,
    df_after_tp,
    df_after_fp,
    classes,
    view_label,
):
    """Compute and return summary statistics as text."""
    lines = []
    lines.append("# Multilayer HP Before/After Comparison Summary\n")
    lines.append(f"- Verification comparison semantics: {view_label}\n")

    for df, label in [(merged_tp, 'TP'), (merged_fp, 'FP')]:
        lines.append(f"\n## {label} (n={len(df)} matched regions)\n")

        # PassedGating
        g_b = df['PassedGating_before'].sum()
        g_a = df['PassedGating_after'].sum()
        lines.append(f"### PassedGating")
        lines.append(f"- Before: {g_b} ({g_b/len(df)*100:.2f}%)")
        lines.append(f"- After:  {g_a} ({g_a/len(df)*100:.2f}%)")
        lines.append(f"- Delta:  +{g_a - g_b} ({(g_a-g_b)/len(df)*100:.2f}pp)")
        n_new_pass = ((~df['PassedGating_before']) & df['PassedGating_after']).sum()
        n_lost = (df['PassedGating_before'] & (~df['PassedGating_after'])).sum()
        lines.append(f"- New passes (False→True): {n_new_pass}")
        lines.append(f"- Lost passes (True→False): {n_lost}")

        # Explicitly selected verification view
        lines.append(f"\n### {view_label}")
        for cls in classes:
            nb = (df['_verification_class_selected_before'] == cls).sum()
            na = (df['_verification_class_selected_after'] == cls).sum()
            lines.append(f"- {cls}: {nb} → {na} (delta: {na-nb:+d})")

        # Transitions (non-diagonal)
        lines.append(f"\n### Key Transitions")
        for src in classes:
            for dst in classes:
                if src == dst:
                    continue
                n = (
                    (df['_verification_class_selected_before'] == src)
                    & (df['_verification_class_selected_after'] == dst)
                ).sum()
                if n > 0:
                    lines.append(f"- {src} → {dst}: {n}")

        # GlobalP
        p_b = -np.log10(df['GlobalP_before'].clip(lower=1e-20))
        p_a = -np.log10(df['GlobalP_after'].clip(lower=1e-20))
        delta = p_a - p_b
        lines.append(f"\n### GlobalP (-log10)")
        lines.append(f"- Mean delta: {delta.mean():.4f}")
        lines.append(f"- Median delta: {delta.median():.4f}")
        lines.append(f"- Improved: {(delta > 0.01).sum()}")
        lines.append(f"- Same: {(delta.abs() <= 0.01).sum()}")
        lines.append(f"- Worse: {(delta < -0.01).sum()}")

        # HeuristicScore
        hs_b = df['HeuristicScore_before']
        hs_a = df['HeuristicScore_after']
        hs_delta = hs_a - hs_b
        lines.append(f"\n### HeuristicScore")
        lines.append(f"- Mean delta: {hs_delta.mean():.4f}")
        lines.append(f"- Median delta: {hs_delta.median():.4f}")
        lines.append(f"- Improved: {(hs_delta > 0.01).sum()}")
        lines.append(f"- Same: {(hs_delta.abs() <= 0.01).sum()}")

        # Quality_Score
        qs_b = df['Quality_Score_before']
        qs_a = df['Quality_Score_after']
        qs_delta = qs_a - qs_b
        lines.append(f"\n### Quality_Score")
        lines.append(f"- Mean delta: {qs_delta.mean():.4f}")
        lines.append(f"- Median delta: {qs_delta.median():.4f}")
        lines.append(f"- Improved: {(qs_delta > 0.01).sum()}")
        lines.append(f"- Same: {(qs_delta.abs() <= 0.01).sum()}")
        lines.append(f"- Worse: {(qs_delta < -0.01).sum()}")

        # Quality_Tier
        lines.append(f"\n### Quality_Tier")
        tiers = ['Excellent', 'Good', 'Fair', 'Poor']
        for t in tiers:
            nb = (df['Quality_Tier_before'] == t).sum()
            na = (df['Quality_Tier_after'] == t).sum()
            lines.append(f"- {t}: {nb} → {na} (delta: {na-nb:+d})")

    # New columns summary
    lines.append(f"\n\n## New Columns Summary (After only)\n")
    for df, label in [(df_after_tp, 'TP'), (df_after_fp, 'FP')]:
        lines.append(f"\n### {label}")
        n = len(df)
        hp_fam_sig = (df['GlobalP_HPFamily'] < 0.05).sum()
        lines.append(f"- HP Family sig (p<0.05): {hp_fam_sig} ({hp_fam_sig/n*100:.1f}%)")
        hp_fine_valid = (df['HPFine_NGroups_CF'] >= 2).sum()
        lines.append(f"- HP Fine valid (NGroups>=2): {hp_fine_valid} ({hp_fine_valid/n*100:.1f}%)")
        for ng in sorted(df['HPFine_NGroups_CF'].unique()):
            cnt = (df['HPFine_NGroups_CF'] == ng).sum()
            lines.append(f"  - NGroups={int(ng)}: {cnt} ({cnt/n*100:.1f}%)")

    return '\n'.join(lines)


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Check if before data exists
    if not BEFORE_TP.exists() or not BEFORE_FP.exists():
        print(f"ERROR: Before data not found at {BEFORE_TP} or {BEFORE_FP}")
        print("Run old binary first to generate before data.")
        sys.exit(1)

    print("Loading data...")
    merged_tp, _, df_after_tp, tp_contract = load_and_merge(
        BEFORE_TP,
        AFTER_TP,
        'TP',
        args.verification_view,
        allow_unversioned_v1=args.allow_unversioned_v1,
    )
    merged_fp, _, df_after_fp, fp_contract = load_and_merge(
        BEFORE_FP,
        AFTER_FP,
        'FP',
        args.verification_view,
        allow_unversioned_v1=args.allow_unversioned_v1,
    )
    classes = tp_contract["before"]["categories"]
    if classes != fp_contract["before"]["categories"]:
        raise SchemaContractError("TP and FP comparisons use different category contracts")
    view_label = tp_contract["comparison_semantics"]
    if view_label != fp_contract["comparison_semantics"]:
        raise SchemaContractError("TP and FP comparisons use different semantic views")

    provenance = {
        "requested_view": args.verification_view,
        "allow_unversioned_v1": args.allow_unversioned_v1,
        "categories": classes,
        "TP": tp_contract,
        "FP": fp_contract,
    }
    with open(OUT_DIR / 'verification_schema_provenance.json', 'w') as f:
        json.dump(provenance, f, indent=2)

    print("\nGenerating figures...")
    fig1_gating_transition(merged_tp, merged_fp)
    fig2_verification_class_transition(merged_tp, merged_fp, classes, view_label)
    fig3_globalp_comparison(merged_tp, merged_fp)
    fig4_heuristic_score(merged_tp, merged_fp)
    fig5_quality(merged_tp, merged_fp)
    fig6_new_columns(df_after_tp, df_after_fp)
    df_sanity = fig7_sanity_check(merged_tp)

    print("\nSanity check results:")
    print(df_sanity.to_string(index=False))

    summary = compute_summary(
        merged_tp,
        merged_fp,
        df_after_tp,
        df_after_fp,
        classes,
        view_label,
    )
    summary_path = OUT_DIR / 'summary_stats.txt'
    with open(summary_path, 'w') as f:
        f.write(summary)
    print(f"\nSummary written to {summary_path}")
    print(summary)

    print(f"\nAll figures saved to {OUT_DIR}/")


if __name__ == '__main__':
    main()
