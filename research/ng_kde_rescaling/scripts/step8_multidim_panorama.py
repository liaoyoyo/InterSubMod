#!/usr/bin/env python3
"""
Step 8: Multi-dimensional panorama — combine ALL parameters at once.

Rationale: §14 single-feature sub-schemes fail (max fold 1.17x on S4).
User asks: if we look at ALL parameter dimensions jointly (not one at a time),
what are the high-TP cells and their Pareto frontier?

Dimensions:
  LOH_Subtype (5) x AF_class (3) x cn_tier_F (5) x NG (1-4) x NR_band (3)
  Total 5 * 3 * 5 * 4 * 3 = 900 bins (before sparsity)

Outputs:
  tpfp_5d_cube_HCC1395_TO.tsv   — all cells with n>=1
  tpfp_5d_pareto_HCC1395_TO.tsv — Pareto-sorted high-TP cells
  tpfp_5d_cumulative_pareto.tsv — cumulative (sorted by TP rate desc) envelope
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import DATA_DIR, FIG_DIR, PALETTE, CN_STRATEGIES, assign_cn_tier


MASTER_TSV = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"

NR_BANDS = [("low", 0, 60), ("mid", 60, 120), ("high", 120, 10_000)]


def nr_band_label(nr: float) -> str:
    for name, lo, hi in NR_BANDS:
        if lo <= nr < hi:
            return name
    return "unknown"


def build_5d_cube(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("None")
    df["cn_tier_F"] = assign_cn_tier(df["CovM_used"], CN_STRATEGIES["F_SEQC2_grounded"])
    df["nr_band"] = df["NumReads"].apply(nr_band_label)

    keys = ["LOH_Subtype", "AF_class", "cn_tier_F", "HPFineNGroups", "nr_band"]
    g = df.groupby(keys, observed=True)
    agg = g.agg(
        n=("tp_label", "size"),
        n_tp=("tp_label", "sum"),
    ).reset_index()
    agg["n_fp"] = agg["n"] - agg["n_tp"]
    agg["tp_rate"] = agg["n_tp"] / agg["n"]
    agg["tp_fp_ratio"] = np.where(agg["n_fp"] == 0, np.inf, agg["n_tp"] / agg["n_fp"])
    return agg


def pareto_sort(cube: pd.DataFrame, min_n: int = 20) -> pd.DataFrame:
    sub = cube[cube["n"] >= min_n].copy()
    # Primary sort: tp_rate desc, tiebreak: n desc
    sub = sub.sort_values(["tp_rate", "n"], ascending=[False, False]).reset_index(drop=True)
    sub["rank"] = np.arange(1, len(sub) + 1)
    return sub


def cumulative_envelope(cube: pd.DataFrame, min_n: int = 20,
                        total_tp: int = 0, total_fp: int = 0) -> pd.DataFrame:
    """Sort cells by tp_rate desc, accumulate n_tp, n_fp, compute cum tp_rate + recall."""
    sub = cube[cube["n"] >= min_n].copy().sort_values("tp_rate", ascending=False).reset_index(drop=True)
    sub["cum_n"] = sub["n"].cumsum()
    sub["cum_n_tp"] = sub["n_tp"].cumsum()
    sub["cum_n_fp"] = sub["n_fp"].cumsum()
    sub["cum_tp_rate"] = sub["cum_n_tp"] / sub["cum_n"]
    sub["cum_tp_recall"] = sub["cum_n_tp"] / total_tp if total_tp > 0 else np.nan
    sub["cum_fp_recall"] = sub["cum_n_fp"] / total_fp if total_fp > 0 else np.nan
    sub["cum_tp_fp_ratio"] = np.where(sub["cum_n_fp"] == 0, np.inf, sub["cum_n_tp"] / sub["cum_n_fp"])
    return sub


def fig11_panorama(cube: pd.DataFrame, env: pd.DataFrame, title_sample: str,
                   baseline_ratio: float, baseline_prec: float, out_path: Path) -> None:
    plt.rcParams.update({
        "figure.facecolor": PALETTE["bg"],
        "axes.facecolor": PALETTE["bg"],
        "axes.edgecolor": PALETTE["dark"],
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
    })

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Panel A: scatter of all cells (tp_rate vs n), colored by fold-improvement
    ax = axes[0]
    valid = cube[cube["n"] >= 20].copy()
    ratios = np.where(valid["n_fp"] == 0, 100.0, valid["tp_fp_ratio"])
    folds = ratios / baseline_ratio if baseline_ratio > 0 else ratios
    sc = ax.scatter(valid["n"], valid["tp_rate"],
                    c=np.log10(folds.clip(0.1, 100)),
                    cmap="RdYlGn", s=valid["n"].clip(20, 500) / 4,
                    alpha=0.7, edgecolors=PALETTE["dark"], linewidths=0.5)
    ax.axhline(baseline_prec, color=PALETTE["accent"], linestyle="--", linewidth=2,
               label=f"baseline precision = {baseline_prec:.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("n (cell size, log scale)")
    ax.set_ylabel("TP rate (precision within cell)")
    ax.set_title(f"{title_sample} | 5D cube cells (n>=20)\n"
                 f"colored by log10(fold-improvement vs baseline)")
    ax.set_ylim(0, 1.05)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("log10(fold)")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel B: cumulative envelope (cum_tp_rate vs cum_tp_recall)
    ax = axes[1]
    ax.plot(env["cum_tp_recall"], env["cum_tp_rate"],
            color=PALETTE["dark"], linewidth=2.5, label="cum envelope (cells sorted by tp_rate desc)")
    ax.scatter(env["cum_tp_recall"].iloc[:5], env["cum_tp_rate"].iloc[:5],
               s=90, color=PALETTE["green"], edgecolor=PALETTE["dark"], zorder=5,
               label="top-5 cells")
    ax.axhline(baseline_prec, color=PALETTE["accent"], linestyle="--", linewidth=2,
               label=f"baseline precision = {baseline_prec:.3f}")
    # Mark 80%/90%/95% purity cutoffs
    for threshold in [0.95, 0.90, 0.80]:
        eligible = env[env["cum_tp_rate"] >= threshold]
        if not eligible.empty:
            last = eligible.iloc[-1]
            ax.scatter([last["cum_tp_recall"]], [last["cum_tp_rate"]],
                       s=120, marker="*", color=PALETTE["red"], zorder=6,
                       edgecolor=PALETTE["dark"])
            ax.annotate(f"cum_TP_rate≥{threshold:.2f}\nrecall={last['cum_tp_recall']:.2%}, "
                        f"n={int(last['cum_n'])}",
                        (last["cum_tp_recall"], last["cum_tp_rate"]),
                        textcoords="offset points", xytext=(8, -18),
                        fontsize=9, color=PALETTE["dark"])
    ax.set_xlabel("cumulative TP recall")
    ax.set_ylabel("cumulative TP rate (purity)")
    ax.set_title(f"{title_sample} | Pareto envelope\n"
                 f"(trading purity vs recall by taking top-k cells)")
    ax.set_xlim(-0.02, 1.0)
    ax.set_ylim(0.4, 1.05)
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Figure v2-11 | Multi-dimensional panorama (LOH x AF x CN x NG x NR)",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig11] wrote {out_path}")


def main() -> None:
    df = pd.read_csv(MASTER_TSV, sep="\t", compression="gzip", low_memory=False)
    print(f"[step8] loaded master: {df.shape}")

    # ---- HCC1395 TO (main testbed) ----
    to = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")].copy()
    total_tp = int(to["tp_label"].sum()); total_fp = int((to["tp_label"] == 0).sum())
    baseline_ratio = total_tp / total_fp if total_fp > 0 else np.nan
    baseline_prec = total_tp / len(to) if len(to) > 0 else np.nan
    print(f"[HCC1395 TO] N={len(to)}  TP={total_tp}  FP={total_fp}  baseline TP:FP={baseline_ratio:.2f}")

    cube_to = build_5d_cube(to)
    cube_to.to_csv(DATA_DIR / "tpfp_5d_cube_HCC1395_TO.tsv", sep="\t", index=False)
    print(f"[step8] 5D cube HCC1395 TO: {cube_to.shape}, max cell n={cube_to['n'].max()}")

    pareto = pareto_sort(cube_to, min_n=20)
    pareto.to_csv(DATA_DIR / "tpfp_5d_pareto_HCC1395_TO.tsv", sep="\t", index=False)

    print(f"\n[panorama] Top 15 high-TP cells (n>=20) HCC1395 TO:")
    cols_show = ["rank", "LOH_Subtype", "AF_class", "cn_tier_F", "HPFineNGroups", "nr_band",
                 "n", "n_tp", "n_fp", "tp_rate", "tp_fp_ratio"]
    print(pareto.head(15)[cols_show].to_string(index=False))

    env_to = cumulative_envelope(cube_to, min_n=20, total_tp=total_tp, total_fp=total_fp)
    env_to.to_csv(DATA_DIR / "tpfp_5d_cumulative_envelope_HCC1395_TO.tsv", sep="\t", index=False)

    print(f"\n[envelope] at cum_TP_rate>=0.95:")
    elig = env_to[env_to["cum_tp_rate"] >= 0.95]
    if not elig.empty:
        r = elig.iloc[-1]
        print(f"  k={int(r.name)+1} cells, cum_n={int(r['cum_n'])}, cum_n_tp={int(r['cum_n_tp'])}, "
              f"cum_n_fp={int(r['cum_n_fp'])}, cum_tp_rate={r['cum_tp_rate']:.4f}, "
              f"recall={r['cum_tp_recall']:.3%}, fold={r['cum_tp_fp_ratio']/baseline_ratio:.2f}x")

    for th in [0.90, 0.85, 0.80]:
        elig = env_to[env_to["cum_tp_rate"] >= th]
        if not elig.empty:
            r = elig.iloc[-1]
            print(f"  cum_TP>={th}: k={int(r.name)+1}, cum_n={int(r['cum_n'])}, "
                  f"recall={r['cum_tp_recall']:.3%}, fold={r['cum_tp_fp_ratio']/baseline_ratio:.2f}x")

    fig11_panorama(cube_to, env_to, "HCC1395 TO", baseline_ratio, baseline_prec,
                   FIG_DIR / "fig_v2_11_panorama_HCC1395_TO.png")

    # ---- COLO829 paired_full (secondary testbed) ----
    co = df[(df["sample"] == "COLO829") & (df["mode"] == "paired_full")].copy()
    c_tp = int(co["tp_label"].sum()); c_fp = int((co["tp_label"] == 0).sum())
    c_baseline_ratio = c_tp / c_fp if c_fp > 0 else np.nan
    c_baseline_prec = c_tp / len(co) if len(co) > 0 else np.nan
    print(f"\n[COLO829 pf] N={len(co)}  TP={c_tp}  FP={c_fp}  baseline TP:FP={c_baseline_ratio:.2f}")

    cube_co = build_5d_cube(co)
    cube_co.to_csv(DATA_DIR / "tpfp_5d_cube_COLO829_pf.tsv", sep="\t", index=False)
    pareto_co = pareto_sort(cube_co, min_n=20)
    pareto_co.to_csv(DATA_DIR / "tpfp_5d_pareto_COLO829_pf.tsv", sep="\t", index=False)

    print(f"\n[panorama] Top 10 high-TP cells (n>=20) COLO829 pf:")
    print(pareto_co.head(10)[cols_show].to_string(index=False))

    env_co = cumulative_envelope(cube_co, min_n=20, total_tp=c_tp, total_fp=c_fp)
    env_co.to_csv(DATA_DIR / "tpfp_5d_cumulative_envelope_COLO829_pf.tsv", sep="\t", index=False)

    for th in [0.98, 0.95, 0.90]:
        elig = env_co[env_co["cum_tp_rate"] >= th]
        if not elig.empty:
            r = elig.iloc[-1]
            print(f"  [COLO829] cum_TP>={th}: k={int(r.name)+1}, cum_n={int(r['cum_n'])}, "
                  f"recall={r['cum_tp_recall']:.3%}, fold={r['cum_tp_fp_ratio']/c_baseline_ratio:.2f}x")

    fig11_panorama(cube_co, env_co, "COLO829 paired_full", c_baseline_ratio, c_baseline_prec,
                   FIG_DIR / "fig_v2_11b_panorama_COLO829.png")

    print("\n[step8] DONE")


if __name__ == "__main__":
    main()
