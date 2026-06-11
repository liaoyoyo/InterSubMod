#!/usr/bin/env python3
"""Step 4 stage 2 — Per-sample 3-axis grid (LOH × HP bucket × Coverage_Multiple) for
HCC1395 + 4 phaseD samples. 50 cells × 5 samples, V6-only.

Per cell:
  - n, n_TP, n_FP, TP_rate, FP_rate, Wilson 95% CI
  - log-odds vs rest, Fisher exact p
  - BH-FDR on fisher_p (within sample)
  - lr_beta0 (logit P(TP) ~ NG + AF + n_reads); intercept = adjusted log-odds
  - lr_dev_explained, LRT per covariate
  - Confound guard (simplified 5-axis): chr-stratified Mantel-Haenszel + permutation +
    within-group OLS residualize on n_reads + on caller_af
  - powered/marginal/underpowered flag

Falls back to V6_off HP bucket + n_reads quantile proxy when master.tsv covariate
(LOH/CN) is missing.

Outputs:
  step4_cross_sample_extension/per_sample_grid/{sample}_grid.tsv
  step4_cross_sample_extension/intermediate/per_sample_build_log.txt
  step4_cross_sample_extension/figures/{sample}_facets.png
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from _common import (
    COV_LEVELS,
    FIGURES_DIR,
    HP_LEVELS,
    INTERMEDIATE_DIR,
    LOH_LEVELS,
    MARGINAL_THRESHOLD,
    PER_SAMPLE_GRID,
    POWERED_THRESHOLD,
    SAMPLES,
    annotate_axes,
    cell_id,
    joined_subset,
    load_sample_master,
    matplotlib_setup,
    wilson_ci,
)

PER_SAMPLE_GRID.mkdir(parents=True, exist_ok=True)
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    mask = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    if mask.sum() == 0:
        return q
    p_valid = p[mask]
    n = len(p_valid)
    order = np.argsort(p_valid)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    adj = p_valid * n / ranks
    adj_sorted = adj[order]
    adj_cummin = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj[order] = np.minimum(adj_cummin, 1.0)
    q[mask] = adj
    return q


def fit_cell_lr(cell_df: pd.DataFrame) -> dict:
    y = (cell_df["label"] == "TP").astype(int).values
    base = {
        "lr_beta0": float("nan"),
        "lr_beta_NG": float("nan"),
        "lr_beta_AF": float("nan"),
        "lr_beta_NR": float("nan"),
        "lr_dev_null": float("nan"),
        "lr_dev_full": float("nan"),
        "lr_dev_explained": float("nan"),
        "lr_lrt_p_NG": float("nan"),
        "lr_lrt_p_AF": float("nan"),
        "lr_lrt_p_NR": float("nan"),
        "lr_converged": False,
        "lr_skip_reason": "",
    }
    if len(y) < 10 or len(np.unique(y)) < 2:
        base["lr_skip_reason"] = (
            "insufficient_y_variation" if len(y) >= 10 else "n<10"
        )
        return base
    feat = cell_df[["NG_v6off", "caller_af_num", "n_reads_v6off"]].copy()
    if feat["caller_af_num"].isna().all():
        # fall back to NG + n_reads only
        feat = feat.drop(columns=["caller_af_num"]).assign(caller_af_num=0.5)
        base["lr_skip_reason"] = "af_missing_imputed_0.5"
    feat = feat.fillna(feat.median(numeric_only=True))
    X_full = sm.add_constant(feat.astype(float).values, has_constant="add")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            m_full = sm.Logit(y, X_full).fit(disp=False, maxiter=200)
        except Exception as e:
            base["lr_skip_reason"] = f"logit_fit_failed:{type(e).__name__}"
            return base

    X_null = np.ones((len(y), 1))
    try:
        m_null = sm.Logit(y, X_null).fit(disp=False, maxiter=200)
        dev_null = -2 * m_null.llf
    except Exception:
        dev_null = float("nan")
    dev_full = -2 * m_full.llf
    dev_explained = (
        1 - dev_full / dev_null if not np.isnan(dev_null) and dev_null > 0 else float("nan")
    )

    def lrt_drop(idx_drop: int):
        cols = [0] + [i for i in [1, 2, 3] if i != idx_drop]
        X_red = X_full[:, cols]
        try:
            m_red = sm.Logit(y, X_red).fit(disp=False, maxiter=200)
            stat = -2 * (m_red.llf - m_full.llf)
            return float(stats.chi2.sf(stat, df=1))
        except Exception:
            return float("nan")

    base.update(
        {
            "lr_beta0": float(m_full.params[0]),
            "lr_beta_NG": float(m_full.params[1]),
            "lr_beta_AF": float(m_full.params[2]),
            "lr_beta_NR": float(m_full.params[3]),
            "lr_dev_null": float(dev_null),
            "lr_dev_full": float(dev_full),
            "lr_dev_explained": float(dev_explained),
            "lr_lrt_p_NG": lrt_drop(1),
            "lr_lrt_p_AF": lrt_drop(2),
            "lr_lrt_p_NR": lrt_drop(3),
            "lr_converged": True,
        }
    )
    return base


def confound_guards(cell_df: pd.DataFrame, df_universe: pd.DataFrame, rng: np.random.Generator) -> dict:
    """Light-weight confound suite (Step 4 simplified):
      - chr-stratified Mantel-Haenszel-like Fisher (chr by chr)
      - permutation (shuffle labels within sample), 500 reps
      - within-group OLS residualize on n_reads
      - within-group OLS residualize on caller_af
    """
    y = (cell_df["label"] == "TP").astype(int).values
    out = {
        "residual_TP_rate_NR": float("nan"),
        "residual_TP_rate_AF": float("nan"),
        "permutation_extreme_flag": False,
        "permutation_p": float("nan"),
        "mh_chr_stratified_min_p": float("nan"),
        "n_chr_with_signal": 0,
    }
    if len(y) < 10:
        return out
    raw_rate = y.mean()

    # Guard 1: residualize on n_reads
    try:
        x = cell_df["n_reads_v6off"].astype(float).values.reshape(-1, 1)
        X = sm.add_constant(x, has_constant="add")
        m = sm.OLS(y, X).fit()
        pred = m.predict(X)
        out["residual_TP_rate_NR"] = float(raw_rate - (pred.mean() - raw_rate))
    except Exception:
        pass

    # Guard 2: residualize on caller_af
    try:
        af = cell_df["caller_af_num"].astype(float).values
        mask = ~np.isnan(af)
        if mask.sum() >= 10 and len(np.unique(y[mask])) >= 2:
            X = sm.add_constant(af[mask].reshape(-1, 1), has_constant="add")
            m = sm.OLS(y[mask], X).fit()
            pred = m.predict(X)
            out["residual_TP_rate_AF"] = float(y[mask].mean() - (pred.mean() - y[mask].mean()))
    except Exception:
        pass

    # Guard 3: permutation test (shuffle within universe label set; compare cell rate
    # vs distribution of same-size random draws from universe).
    n_cell = len(y)
    universe_y = (df_universe["label"] == "TP").astype(int).values
    if len(universe_y) >= n_cell + 100:
        rates = np.empty(500)
        for i in range(500):
            idx = rng.choice(len(universe_y), size=n_cell, replace=False)
            rates[i] = universe_y[idx].mean()
        lower = np.percentile(rates, 2.5)
        upper = np.percentile(rates, 97.5)
        out["permutation_extreme_flag"] = bool(raw_rate < lower or raw_rate > upper)
        # two-sided p
        p_tail = (
            (np.sum(rates >= raw_rate) + 1) / (len(rates) + 1)
            if raw_rate >= rates.mean()
            else (np.sum(rates <= raw_rate) + 1) / (len(rates) + 1)
        )
        out["permutation_p"] = float(min(1.0, 2 * p_tail))

    # Guard 4: chr-stratified Fisher
    chr_ps = []
    rest_y = (df_universe["label"] == "TP").values
    rest_chr = df_universe["chr"].values
    for chrom, sub in cell_df.groupby("chr"):
        sub_y = (sub["label"] == "TP").astype(int).values
        if len(sub_y) < 5:
            continue
        rest_mask = (rest_chr == chrom)
        rest_tp = int(rest_y[rest_mask].sum() - sub_y.sum())
        rest_fp = int(rest_mask.sum() - rest_y[rest_mask].sum() - (len(sub_y) - sub_y.sum()))
        n_tp = int(sub_y.sum())
        n_fp = int(len(sub_y) - n_tp)
        if rest_tp + rest_fp < 10:
            continue
        try:
            _, p = stats.fisher_exact([[n_tp, n_fp], [rest_tp, rest_fp]])
            chr_ps.append(p)
        except Exception:
            continue
    if chr_ps:
        out["mh_chr_stratified_min_p"] = float(np.min(chr_ps))
        out["n_chr_with_signal"] = int(sum(p < 0.05 for p in chr_ps))
    return out


def build_grid_for_sample(sample: str) -> tuple[pd.DataFrame, dict]:
    print(f"\n[Step 4 grid] {sample} START")
    df_all = load_sample_master(sample)
    df = joined_subset(df_all)
    df = annotate_axes(df, prefix="V6_off")
    n_total = len(df)
    # When caller_af completely missing in joined subset, fall back to NumReads quantiles
    # to keep covariate logic running (LR will skip AF if all-NaN).
    if df["caller_af_num"].isna().all():
        print(f"  [warn] {sample}: all caller_af NA in master-joined subset → LR uses NG+n_reads only")

    global_tp = (df["label"] == "TP").mean() if n_total else float("nan")
    global_fp = (df["label"] == "FP").mean() if n_total else float("nan")
    print(
        f"  joined universe n={n_total:,}, global TP_rate={global_tp:.4f}, FP_rate={global_fp:.4f}"
    )

    rng = np.random.default_rng(seed=42)
    cells = []
    for loh in LOH_LEVELS:
        for hp in HP_LEVELS:
            for cov in COV_LEVELS:
                mask = (
                    (df["loh_side_norm"] == loh)
                    & (df["hp_bucket"] == hp)
                    & (df["cov_bin"] == cov)
                )
                cell_df = df[mask]
                n = len(cell_df)
                n_tp = int((cell_df["label"] == "TP").sum())
                n_fp = int((cell_df["label"] == "FP").sum())
                tp_rate = (n_tp / n) if n else float("nan")
                fp_rate = (n_fp / n) if n else float("nan")
                lo, hi = wilson_ci(n_tp, n) if n else (float("nan"), float("nan"))
                rest_tp = (df["label"] == "TP").sum() - n_tp
                rest_fp = (df["label"] == "FP").sum() - n_fp
                fisher_p = float("nan")
                log_odds = float("nan")
                if n > 0 and rest_tp > 0 and rest_fp > 0:
                    try:
                        _, fisher_p = stats.fisher_exact([[n_tp, n_fp], [rest_tp, rest_fp]])
                    except Exception:
                        pass
                    if n_tp > 0 and n_fp > 0:
                        log_odds = float(np.log((n_tp / n_fp) / (rest_tp / rest_fp)))
                lr = fit_cell_lr(cell_df) if n >= POWERED_THRESHOLD else {
                    "lr_beta0": float("nan"),
                    "lr_beta_NG": float("nan"),
                    "lr_beta_AF": float("nan"),
                    "lr_beta_NR": float("nan"),
                    "lr_dev_null": float("nan"),
                    "lr_dev_full": float("nan"),
                    "lr_dev_explained": float("nan"),
                    "lr_lrt_p_NG": float("nan"),
                    "lr_lrt_p_AF": float("nan"),
                    "lr_lrt_p_NR": float("nan"),
                    "lr_converged": False,
                    "lr_skip_reason": "not_powered",
                }
                guards = confound_guards(cell_df, df, rng) if n >= POWERED_THRESHOLD else {
                    "residual_TP_rate_NR": float("nan"),
                    "residual_TP_rate_AF": float("nan"),
                    "permutation_extreme_flag": False,
                    "permutation_p": float("nan"),
                    "mh_chr_stratified_min_p": float("nan"),
                    "n_chr_with_signal": 0,
                }
                rec = {
                    "sample": sample,
                    "cell_id": cell_id(loh, hp, cov),
                    "loh_side": loh,
                    "hp_bucket": hp,
                    "cov_bin": cov,
                    "n": n,
                    "n_TP": n_tp,
                    "n_FP": n_fp,
                    "TP_rate": tp_rate,
                    "FP_rate": fp_rate,
                    "TP_wilson_lo": lo,
                    "TP_wilson_hi": hi,
                    "TP_enrichment": (tp_rate / global_tp) if (n and global_tp) else float("nan"),
                    "FP_enrichment": (fp_rate / global_fp) if (n and global_fp) else float("nan"),
                    "log_odds": log_odds,
                    "fisher_p": fisher_p,
                    **lr,
                    **guards,
                    "powered_flag": n >= POWERED_THRESHOLD,
                    "marginal_flag": MARGINAL_THRESHOLD <= n < POWERED_THRESHOLD,
                    "underpowered_flag": 0 < n < MARGINAL_THRESHOLD,
                    "empty_flag": n == 0,
                }
                cells.append(rec)
    grid = pd.DataFrame(cells)
    grid["fisher_q_bh"] = bh_fdr(grid["fisher_p"].values)

    out_path = PER_SAMPLE_GRID / f"{sample}_grid.tsv"
    grid.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
    print(f"  wrote {out_path.name}: cells={len(grid)}, powered={int(grid['powered_flag'].sum())}")

    summary = {
        "sample": sample,
        "n_universe": int(n_total),
        "global_TP_rate": float(global_tp) if global_tp == global_tp else None,
        "global_FP_rate": float(global_fp) if global_fp == global_fp else None,
        "powered_cells": int(grid["powered_flag"].sum()),
        "marginal_cells": int(grid["marginal_flag"].sum()),
        "underpowered_cells": int(grid["underpowered_flag"].sum()),
        "empty_cells": int(grid["empty_flag"].sum()),
        "lr_converged": int(grid["lr_converged"].fillna(False).sum()),
    }
    return grid, summary


def make_facets(sample: str, grid: pd.DataFrame):
    plt = matplotlib_setup()
    fig, axes = plt.subplots(1, len(COV_LEVELS), figsize=(20, 4.2), sharey=True)
    for ax, cov in zip(axes, COV_LEVELS):
        sub = grid[grid["cov_bin"] == cov]
        mat = np.full((len(LOH_LEVELS), len(HP_LEVELS)), np.nan)
        ann = np.full((len(LOH_LEVELS), len(HP_LEVELS)), "", dtype=object)
        for i, loh in enumerate(LOH_LEVELS):
            for j, hp in enumerate(HP_LEVELS):
                r = sub[(sub["loh_side"] == loh) & (sub["hp_bucket"] == hp)]
                if len(r) and r.iloc[0]["n"] > 0:
                    mat[i, j] = r.iloc[0]["TP_rate"]
                    n = int(r.iloc[0]["n"])
                    flag = "*" if r.iloc[0]["powered_flag"] else (
                        "." if r.iloc[0]["marginal_flag"] else ""
                    )
                    ann[i, j] = f"{r.iloc[0]['TP_rate']:.2f}\n(n={n}){flag}"
        im = ax.imshow(mat, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
        for i in range(len(LOH_LEVELS)):
            for j in range(len(HP_LEVELS)):
                if ann[i, j]:
                    ax.text(j, i, ann[i, j], ha="center", va="center", fontsize=7, color="black")
        ax.set_xticks(range(len(HP_LEVELS)))
        ax.set_xticklabels(HP_LEVELS, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(LOH_LEVELS)))
        ax.set_yticklabels(LOH_LEVELS)
        ax.set_title(cov, fontsize=10)
    fig.suptitle(f"{sample} TP rate (LOH × HP × cov)  *=powered, .=marginal", fontsize=11)
    cbar = fig.colorbar(im, ax=axes, fraction=0.015, pad=0.02)
    cbar.set_label("TP rate")
    out = FIGURES_DIR / f"{sample}_facets.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.name}")


def main():
    summaries = []
    all_rows = []
    for s in SAMPLES:
        grid, summary = build_grid_for_sample(s)
        summaries.append(summary)
        all_rows.append(grid)
        make_facets(s, grid)
    summary_path = INTERMEDIATE_DIR / "step4_grid_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summaries, f, indent=2)
    # combined long
    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(
        PER_SAMPLE_GRID.parent / "step4_per_sample_grid.tsv",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    print(f"\n[Step 4 grid] wrote combined long: step4_per_sample_grid.tsv ({len(combined)} rows)")
    print(f"[Step 4 grid] wrote summary → {summary_path}")


if __name__ == "__main__":
    main()
