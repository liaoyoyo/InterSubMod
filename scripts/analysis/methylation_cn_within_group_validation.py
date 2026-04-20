#!/usr/bin/env python3
"""
P0-B C16 Within-Group Validation: HPFineNGroups × SEQC2_CN confound check

Rationale:
  methylation_cn_validation.py uses POOLED OLS (LinearRegression on all regions)
  to residualize HPFineNGroups against NumReads, then reports 68% confound.
  If NumReads distribution differs across CN strata (e.g., high-CN regions have
  more reads), pooled residualization can introduce collider bias / Simpson
  paradox, artificially shrinking the residual correlation.

Method:
  1. Stratify all HCC1395 TP regions by NumReads deciles (10 bins).
  2. Within each bin, compute raw Spearman(HPFineNGroups, SEQC2_CN).
  3. Compute bootstrap 95% CI per bin (1000 iterations).
  4. If cross-bin correlation is consistent and nonzero → C16 signal is REAL.
     If zero across all bins → C16 pooled residualized is spurious, confound dominates.

Output:
  research/seqc2_cnv_stratification/data/methylation_cn_within_group.tsv
  research/seqc2_cnv_stratification/figures/31_within_group_cn_correlation.png
  Console: per-bin r / CI + overall verdict (REAL / CONFOUND / INCONCLUSIVE)

Decision Framework (CHECKLIST P0-B, R-04 within-group OLS standard):
  - Median |r| >= 0.15 across ≥7/10 bins AND CI excludes 0 → C16 REAL (upgrade POSITIVE)
  - Median |r| < 0.05 across ≥7/10 bins OR all CI include 0 → C16 CONFOUND (downgrade to NEGATIVE/CONDITIONAL)
  - Otherwise → INCONCLUSIVE (need larger samples or alternative stratum)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_PATH = Path('/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv')
OUT_DATA = Path('/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data/methylation_cn_within_group.tsv')
OUT_FIG = Path('/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/figures/31_within_group_cn_correlation.png')

N_BINS = 10
N_BOOT = 1000
MIN_BIN_SIZE = 30
SEED = 42


def bootstrap_spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int, rng: np.random.Generator) -> tuple[float, float]:
    """Return (ci_low, ci_high) for Spearman r via case-resampling bootstrap."""
    rs = []
    n = len(x)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r_b, _ = spearmanr(x[idx], y[idx])
        if not np.isnan(r_b):
            rs.append(r_b)
    if not rs:
        return (np.nan, np.nan)
    return (float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5)))


def main() -> int:
    if not DATA_PATH.exists():
        print(f'[ERROR] Data not found: {DATA_PATH}', file=sys.stderr)
        return 1

    df = pd.read_csv(DATA_PATH, sep='\t')
    tp = df[df['Label'] == 'TP'].copy()
    tp = tp[tp['SEQC2_CN'].notna() & (tp['SEQC2_CN'] > 0)].copy()
    tp = tp.dropna(subset=['HPFineNGroups', 'NumReads']).copy()

    print(f'Loaded {len(tp)} HCC1395 TP regions with valid CN + HPFineNGroups + NumReads')

    overall_r, overall_p = spearmanr(tp['HPFineNGroups'], tp['SEQC2_CN'])
    print(f'\nOverall Spearman(HPFineNGroups, SEQC2_CN) = {overall_r:+.3f}  p={overall_p:.2e}')

    # Stratify by NumReads deciles
    try:
        tp['nr_bin'] = pd.qcut(tp['NumReads'], q=N_BINS, labels=False, duplicates='drop')
    except ValueError as exc:
        print(f'[ERROR] qcut failed: {exc}')
        return 2

    rng = np.random.default_rng(SEED)

    rows = []
    print(f'\n{"Bin":>4} {"NR_min":>8} {"NR_max":>8} {"n":>6} {"r":>8} {"p":>10} {"CI95_low":>9} {"CI95_high":>10} {"Excl_0":>8}')
    print('-' * 80)
    for bin_id in sorted(tp['nr_bin'].dropna().unique()):
        sub = tp[tp['nr_bin'] == bin_id]
        if len(sub) < MIN_BIN_SIZE:
            continue
        r_bin, p_bin = spearmanr(sub['HPFineNGroups'], sub['SEQC2_CN'])
        ci_lo, ci_hi = bootstrap_spearman_ci(
            sub['HPFineNGroups'].values.astype(float),
            sub['SEQC2_CN'].values.astype(float),
            N_BOOT, rng,
        )
        excludes_zero = (ci_lo > 0) or (ci_hi < 0)
        nr_min = float(sub['NumReads'].min())
        nr_max = float(sub['NumReads'].max())
        print(f'{int(bin_id):>4} {nr_min:>8.0f} {nr_max:>8.0f} {len(sub):>6} '
              f'{r_bin:>+8.3f} {p_bin:>10.2e} {ci_lo:>+9.3f} {ci_hi:>+10.3f} {"YES" if excludes_zero else "NO":>8}')
        rows.append({
            'bin': int(bin_id),
            'NR_min': nr_min, 'NR_max': nr_max,
            'n': int(len(sub)),
            'spearman_r': float(r_bin), 'p_value': float(p_bin),
            'ci95_low': ci_lo, 'ci95_high': ci_hi,
            'excludes_zero': bool(excludes_zero),
        })

    result = pd.DataFrame(rows)
    OUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT_DATA, sep='\t', index=False)
    print(f'\nSaved per-bin table: {OUT_DATA}')

    # Verdict
    median_abs_r = float(result['spearman_r'].abs().median())
    n_bins_exclude_zero = int(result['excludes_zero'].sum())
    n_bins_total = int(len(result))
    print(f'\n=== Verdict summary ===')
    print(f'  Bins total: {n_bins_total}')
    print(f'  Bins with CI excluding 0: {n_bins_exclude_zero}')
    print(f'  Median |r| across bins: {median_abs_r:.3f}')

    verdict = 'INCONCLUSIVE'
    reason = ''
    if median_abs_r >= 0.15 and n_bins_exclude_zero >= max(7, int(0.7 * n_bins_total)):
        verdict = 'REAL_SIGNAL'
        reason = f'median |r|={median_abs_r:.3f} >= 0.15 AND {n_bins_exclude_zero}/{n_bins_total} bins exclude 0'
    elif median_abs_r < 0.05:
        verdict = 'CONFOUND_DOMINATED'
        reason = f'median |r|={median_abs_r:.3f} < 0.05; pooled residualized claim is spurious'
    elif n_bins_exclude_zero == 0:
        verdict = 'CONFOUND_DOMINATED'
        reason = 'no bin CI excludes 0 → within-group correlation indistinguishable from noise'
    else:
        reason = f'median |r|={median_abs_r:.3f}, {n_bins_exclude_zero}/{n_bins_total} bins exclude 0 → mixed evidence'

    print(f'\n  C16 within-group verdict: {verdict}')
    print(f'  Reason: {reason}')

    # Figure
    fig, ax = plt.subplots(figsize=(9, 6))
    xs = result['bin'].values
    ys = result['spearman_r'].values
    ci_lo = result['ci95_low'].values
    ci_hi = result['ci95_high'].values
    ax.errorbar(xs, ys, yerr=[ys - ci_lo, ci_hi - ys], fmt='o-', color='steelblue', capsize=4, label='Within-bin Spearman r')
    ax.axhline(0.0, color='black', linestyle='-', alpha=0.4)
    ax.axhline(0.15, color='orange', linestyle=':', alpha=0.5, label='|r|=0.15 (signal threshold)')
    ax.axhline(-0.15, color='orange', linestyle=':', alpha=0.5)
    ax.axhline(overall_r, color='red', linestyle='--', alpha=0.5, label=f'Overall pooled r={overall_r:.3f}')
    ax.set_xlabel('NumReads decile bin (0=low coverage, 9=high)')
    ax.set_ylabel('Within-bin Spearman(HPFineNGroups, SEQC2_CN)')
    ax.set_title(f'P0-B C16 Within-Group Validation\n'
                 f'HCC1395 TP, N={len(tp):,} regions, Verdict: {verdict}',
                 fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_FIG, dpi=150, bbox_inches='tight')
    print(f'\nSaved figure: {OUT_FIG}')

    print(f'\n{"="*60}')
    print(f'Next step: update docs/reports/audit/cards/C16_*.md with verdict: {verdict}')
    print(f'{"="*60}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
