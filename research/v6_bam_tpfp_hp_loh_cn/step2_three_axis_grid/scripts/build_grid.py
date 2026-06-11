"""Step 2 — build 3-axis 50-cell grid with per-cell stats + LR + V3F/V5/V6 trajectory.

Outputs:
  - step2_grid_3d.tsv (per-cell full schema)
  - step2_trajectory_per_cell.tsv (V3F/V5/V6 region count per cell + ratios)
  - step2_marginal_vs_observed.tsv (multiplicative null vs observed)
  - step2_power_curve.png
  - step2_facets.png (5 CN bins × Inner/Outer heatmap)

The companion confound_guard.py adds Guard 1-7 columns into step2_grid_3d.tsv
and produces step2_confound_guard.tsv (per-cell × 7 guards) + step2_top_lists.tsv.
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
    HP_LEVELS,
    HP_PREFIX_PRIMARY,
    LOH_LEVELS,
    MARGINAL_THRESHOLD,
    POWERED_THRESHOLD,
    STEP2_DIR,
    annotate_axes,
    cell_id,
    joined_subset,
    load_master,
    log_msg,
    matplotlib_setup,
    wilson_ci,
)


def fit_cell_lr(cell_df: pd.DataFrame, label_col: str = "label"):
    """Fit logit(P(TP=1)) ~ NG + AF + NumReads. Return dict with β, LRT_p per covariate,
    deviance, dev_explained. If degenerate (no variation in y or <10 samples), returns NaNs.
    """
    y = (cell_df[label_col] == "TP").astype(int).values
    if len(y) < 10 or len(np.unique(y)) < 2:
        return {
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
            "lr_skip_reason": "insufficient_y_variation" if len(y) >= 10 else "n<10",
        }

    X_full = sm.add_constant(
        cell_df[["NG_v6off", "caller_af", "n_reads_v6off"]].astype(float).values,
        has_constant="add",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            m_full = sm.Logit(y, X_full).fit(disp=False, maxiter=200)
        except Exception as e:
            return {
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
                "lr_skip_reason": f"logit_fit_failed:{type(e).__name__}",
            }

    # Null model (intercept only)
    X_null = np.ones((len(y), 1))
    try:
        m_null = sm.Logit(y, X_null).fit(disp=False, maxiter=200)
        dev_null = -2 * m_null.llf
    except Exception:
        dev_null = float("nan")

    dev_full = -2 * m_full.llf
    if not np.isnan(dev_null) and dev_null > 0:
        dev_explained = 1 - dev_full / dev_null
    else:
        dev_explained = float("nan")

    # LRT for each covariate by removing it
    def lrt_drop(idx_drop: int):
        # X_full has columns [const, NG, AF, NR]; drop the idx_drop-th data column (1=NG, 2=AF, 3=NR)
        cols = [0] + [i for i in [1, 2, 3] if i != idx_drop]
        X_red = X_full[:, cols]
        try:
            m_red = sm.Logit(y, X_red).fit(disp=False, maxiter=200)
            stat = -2 * (m_red.llf - m_full.llf)
            return float(stats.chi2.sf(stat, df=1))
        except Exception:
            return float("nan")

    return {
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
        "lr_skip_reason": "",
    }


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values. NaN-safe."""
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
    # enforce monotonicity (descending in rank)
    adj_sorted = adj[order]
    adj_cummin = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj[order] = np.minimum(adj_cummin, 1.0)
    q[mask] = adj
    return q


def compute_trajectory(df_full: pd.DataFrame, cell_id_str: str) -> dict:
    """Compute V3F/V5/V6 region count + ratios for a given cell using each version's
    HP bucket + cov bin partition (re-derived per BAM version). The cell_id corresponds
    to (loh, hp, cov) in V6_off coords; we report 'how many regions fall in same coord
    under V3F_off / V5_off / V6_off respectively'.
    """
    loh, hp, cov = cell_id_str.split("|")
    out = {"cell_id": cell_id_str, "loh_side": loh, "hp_bucket_v6off": hp, "cov_bin": cov}
    for prefix, label in [("V3F_off", "V3F"), ("V5_off", "V5"), ("V6_off", "V6")]:
        ann = annotate_axes(df_full, hp_prefix=prefix)
        mask = (
            (ann["loh_side_norm"] == loh)
            & (ann["hp_bucket"] == hp)
            & (ann["cov_bin"] == cov)
        )
        n = int(mask.sum())
        n_tp = int(((ann["label"] == "TP") & mask).sum())
        out[f"n_{label}"] = n
        out[f"n_TP_{label}"] = n_tp
    out["v5_over_v3f"] = out["n_V5"] / out["n_V3F"] if out["n_V3F"] > 0 else float("nan")
    out["v6_over_v5"] = out["n_V6"] / out["n_V5"] if out["n_V5"] > 0 else float("nan")
    out["v6_over_v3f"] = out["n_V6"] / out["n_V3F"] if out["n_V3F"] > 0 else float("nan")
    out["v5_over_promote_flag"] = (
        out["v5_over_v3f"] > 1.3 if not np.isnan(out["v5_over_v3f"]) else False
    )
    return out


def main():
    log_msg("=== Step 2 build_grid START ===")

    df_all = load_master()
    df = joined_subset(df_all)
    df = annotate_axes(df, hp_prefix=HP_PREFIX_PRIMARY)

    global_tp = (df["label"] == "TP").mean()
    global_fp = (df["label"] == "FP").mean()
    n_total = len(df)
    log_msg(
        f"Universe: n={n_total:,}, global TP_rate={global_tp:.4f}, FP_rate={global_fp:.4f}"
    )

    # ---- Build per-cell records ----
    cells = []
    df_full_for_trajectory = joined_subset(df_all)  # raw subset (no annotation yet)
    trajectories = []

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
                cid = cell_id(loh, hp, cov)
                n_tp = int((cell_df["label"] == "TP").sum())
                n_fp = int((cell_df["label"] == "FP").sum())
                tp_rate = (n_tp / n) if n else float("nan")
                fp_rate = (n_fp / n) if n else float("nan")
                lo, hi = wilson_ci(n_tp, n) if n else (float("nan"), float("nan"))

                # Fisher vs global (2x2: cell TP/FP vs rest)
                rest_tp = (df["label"] == "TP").sum() - n_tp
                rest_fp = (df["label"] == "FP").sum() - n_fp
                if n > 0 and (n_tp + n_fp) > 0:
                    try:
                        _, fisher_p = stats.fisher_exact(
                            [[n_tp, n_fp], [rest_tp, rest_fp]]
                        )
                    except Exception:
                        fisher_p = float("nan")
                else:
                    fisher_p = float("nan")

                # Log-odds (cell TP/FP vs rest)
                if n_tp > 0 and n_fp > 0 and rest_tp > 0 and rest_fp > 0:
                    log_odds = np.log((n_tp / n_fp) / (rest_tp / rest_fp))
                else:
                    log_odds = float("nan")

                # LR per cell (powered only)
                if n >= POWERED_THRESHOLD:
                    lr = fit_cell_lr(cell_df)
                else:
                    lr = {
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

                # Trajectory (only for non-empty)
                if n >= 1:
                    trajectories.append(compute_trajectory(df_full_for_trajectory, cid))

                rec = {
                    "cell_id": cid,
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
                    "TP_enrichment": (tp_rate / global_tp) if n else float("nan"),
                    "FP_enrichment": (fp_rate / global_fp) if n else float("nan"),
                    "log_odds": log_odds,
                    "fisher_p": fisher_p,
                    **lr,
                    "powered_flag": n >= POWERED_THRESHOLD,
                    "marginal_flag": MARGINAL_THRESHOLD <= n < POWERED_THRESHOLD,
                    "underpowered_flag": 0 < n < MARGINAL_THRESHOLD,
                    "empty_flag": n == 0,
                }
                cells.append(rec)

    grid = pd.DataFrame(cells)

    # BH-FDR on fisher_p
    grid["fisher_q_bh"] = bh_fdr(grid["fisher_p"].values)

    # Marginal expectation (3-axis multiplicative null on TP rate)
    # P(TP | loh, hp, cov) ≈ TP_rate_loh × TP_rate_hp × TP_rate_cov / global^2
    p_loh = df.groupby("loh_side_norm")["label"].apply(lambda s: (s == "TP").mean()).to_dict()
    p_hp = df.groupby("hp_bucket")["label"].apply(lambda s: (s == "TP").mean()).to_dict()
    p_cov = df.groupby("cov_bin")["label"].apply(lambda s: (s == "TP").mean()).to_dict()

    def marg_expect(loh, hp, cov):
        if global_tp == 0:
            return float("nan")
        return p_loh.get(loh, np.nan) * p_hp.get(hp, np.nan) * p_cov.get(cov, np.nan) / (global_tp ** 2)

    grid["expected_TP_rate_multiplicative"] = grid.apply(
        lambda r: marg_expect(r["loh_side"], r["hp_bucket"], r["cov_bin"]), axis=1
    )
    grid["surprise"] = grid["TP_rate"] - grid["expected_TP_rate_multiplicative"]

    # Write main grid
    grid_out = STEP2_DIR / "step2_grid_3d.tsv"
    grid.to_csv(grid_out, sep="\t", index=False, float_format="%.6f")
    log_msg(f"Wrote {grid_out.name} ({grid_out.stat().st_size:,} bytes; {len(grid)} cells)")

    # Marginal vs observed table (separate file)
    marg_df = grid[["cell_id", "loh_side", "hp_bucket", "cov_bin", "n", "TP_rate",
                    "expected_TP_rate_multiplicative", "surprise"]].copy()
    marg_df.to_csv(STEP2_DIR / "step2_marginal_vs_observed.tsv", sep="\t", index=False, float_format="%.6f")

    # Trajectory table
    traj_df = pd.DataFrame(trajectories)
    traj_df.to_csv(STEP2_DIR / "step2_trajectory_per_cell.tsv", sep="\t", index=False, float_format="%.6f")
    log_msg(f"Wrote trajectory ({len(traj_df)} cells)")

    # ---- Power curve plot ----
    plt = matplotlib_setup()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Histogram of cell n
    nonempty = grid[grid["n"] > 0]
    axes[0].hist(np.log10(nonempty["n"].clip(lower=1)), bins=20, color="steelblue", edgecolor="black")
    axes[0].axvline(np.log10(POWERED_THRESHOLD), color="red", ls="--", label=f"powered ≥ {POWERED_THRESHOLD}")
    axes[0].axvline(np.log10(MARGINAL_THRESHOLD), color="orange", ls="--", label=f"marginal ≥ {MARGINAL_THRESHOLD}")
    axes[0].set_xlabel("log10(cell n)")
    axes[0].set_ylabel("# cells")
    axes[0].set_title("Cell size distribution (non-empty cells)")
    axes[0].legend()

    # Bar chart of categories
    main_grid = grid[grid["loh_side"] != "UNKNOWN"]
    cats = ["powered", "marginal", "underpowered", "empty"]
    counts_main = [int(main_grid["powered_flag"].sum()),
                   int(main_grid["marginal_flag"].sum()),
                   int(main_grid["underpowered_flag"].sum()),
                   int(main_grid["empty_flag"].sum())]
    counts_full = [int(grid["powered_flag"].sum()),
                   int(grid["marginal_flag"].sum()),
                   int(grid["underpowered_flag"].sum()),
                   int(grid["empty_flag"].sum())]
    x = np.arange(len(cats))
    w = 0.4
    axes[1].bar(x - w/2, counts_main, w, label="50-cell main", color="steelblue")
    axes[1].bar(x + w/2, counts_full, w, label="75-cell incl. UNKNOWN", color="lightcoral")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cats)
    axes[1].set_ylabel("# cells")
    axes[1].set_title("Power stratification")
    axes[1].legend()
    for i, v in enumerate(counts_main):
        axes[1].text(i - w/2, v + 0.3, str(v), ha="center", fontsize=9)
    for i, v in enumerate(counts_full):
        axes[1].text(i + w/2, v + 0.3, str(v), ha="center", fontsize=9)

    plt.tight_layout()
    fig.savefig(STEP2_DIR / "step2_power_curve.png", dpi=120)
    plt.close(fig)
    log_msg("Wrote step2_power_curve.png")

    # ---- Facet heatmap (CN bin × LOH × HP TP-rate) ----
    fig, axes = plt.subplots(1, len(COV_LEVELS), figsize=(20, 4.2), sharey=True)
    for ax, cov in zip(axes, COV_LEVELS):
        sub = grid[grid["cov_bin"] == cov]
        mat = np.full((2, len(HP_LEVELS)), np.nan)
        ann = np.full((2, len(HP_LEVELS)), "", dtype=object)
        for i, loh in enumerate(["Inner", "Outer"]):
            for j, hp in enumerate(HP_LEVELS):
                r = sub[(sub["loh_side"] == loh) & (sub["hp_bucket"] == hp)]
                if len(r) and r.iloc[0]["n"] > 0:
                    mat[i, j] = r.iloc[0]["TP_rate"]
                    n = int(r.iloc[0]["n"])
                    flag = ""
                    if r.iloc[0]["powered_flag"]:
                        flag = "★"
                    elif r.iloc[0]["marginal_flag"]:
                        flag = "·"
                    ann[i, j] = f"{r.iloc[0]['TP_rate']:.2f}\n(n={n}){flag}"
        im = ax.imshow(mat, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
        for i in range(2):
            for j in range(len(HP_LEVELS)):
                if ann[i, j]:
                    ax.text(j, i, ann[i, j], ha="center", va="center", fontsize=7, color="black")
        ax.set_xticks(range(len(HP_LEVELS)))
        ax.set_xticklabels(HP_LEVELS, rotation=45, ha="right", fontsize=8)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["Inner", "Outer"])
        ax.set_title(cov, fontsize=10)
    fig.suptitle("TP rate by (Inner/Outer × HP bucket), faceted by CN bin (★=powered, ·=marginal)", fontsize=11)
    cbar = fig.colorbar(im, ax=axes, fraction=0.015, pad=0.02)
    cbar.set_label("TP rate")
    fig.savefig(STEP2_DIR / "step2_facets.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    log_msg("Wrote step2_facets.png")

    # Summary JSON
    summary = {
        "n_universe": int(n_total),
        "global_TP_rate": float(global_tp),
        "global_FP_rate": float(global_fp),
        "grid_cells_total": int(len(grid)),
        "non_empty_cells": int((grid["n"] > 0).sum()),
        "powered_cells": int(grid["powered_flag"].sum()),
        "marginal_cells": int(grid["marginal_flag"].sum()),
        "underpowered_cells": int(grid["underpowered_flag"].sum()),
        "empty_cells": int(grid["empty_flag"].sum()),
        "main_50_powered": int(grid[grid["loh_side"] != "UNKNOWN"]["powered_flag"].sum()),
        "lr_converged_cells": int(grid["lr_converged"].fillna(False).sum()),
    }
    with open(STEP2_DIR / "intermediate" / "step2_build_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log_msg("=== Step 2 build_grid END ===")
    return summary


if __name__ == "__main__":
    main()
