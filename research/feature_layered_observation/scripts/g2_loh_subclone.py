"""Phase C G2 — LOH & Subclone Annotation observation.

Features in scope (registry group G2):
  Potential_LOH, LOH_Subtype, LOH_Bed_Overlap, LOH_Source,
  LOH_Bed_Annotation, Subclone_ID.

Only LOH_Subtype (5-class: None/Noise/Weak/Strong/Subclone), Potential_LOH
(binary), and LOH_Bed_Overlap (binary, almost all False) are populated
in the master TSV. LOH_Source / Subclone_ID / LOH_Bed_Annotation are not
present (require an extended master) — handled by graceful skip.

LOH_Subtype is a *stratification dimension*, not a continuous AUC target.
The main analysis is per-class TP rate + Wilson CI across 7 samples,
optional ordinal AUC treating the 5 classes as ordered.

Outputs (figures/G2/ + data/G2/):
  Step 1  : fig01_lohsubtype_tp_rate_bars.png  + csv
  Step 2  : fig02_loh_af_cn_heatmap.png        + csv
  Step 3  : fig03_per_sample_grid.png          + csv
  Step 4  : fig04_ordinal_auc.png              + csv
  Step 5  : fig05_residualized_auc.png         + csv
  Step 6  : fig06_per_chr_distribution.png     + csv

Usage:
  python research/feature_layered_observation/scripts/g2_loh_subclone.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

warnings.filterwarnings("ignore", category=RuntimeWarning)

# --- paths ------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[3]
MASTER = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
FIG_DIR = ROOT / "research/feature_layered_observation/figures/G2"
DATA_DIR = ROOT / "research/feature_layered_observation/data/G2"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

# --- constants --------------------------------------------------------------
SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
    "H2009", "H1437", "COLO829",
]
# Subtype order from "most benign/quiet" to "most subclonal"
LOH_ORDER = ["None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"]
LOH_ORD_MAP = {k: i for i, k in enumerate(LOH_ORDER)}
AF_CLASS_ORDER = ["Extreme", "Intermediate", "Near-half"]
COV_ORDER = ["Low", "Normal", "Elevated", "CNV_Loss", "CNV_Gain", "High_Copy"]
MIN_N = 20  # minimum effective n per cell

# --- helpers ----------------------------------------------------------------
def wilson_ci(k: int, n: int, z: float = 1.96):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def auc_binary(labels: np.ndarray, scores: np.ndarray) -> float:
    """Mann-Whitney U based AUC (binary labels 0/1)."""
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=float)
    mask = ~np.isnan(scores)
    labels = labels[mask]
    scores = scores[mask]
    if len(np.unique(labels)) < 2:
        return np.nan
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    # rank-based
    order = np.argsort(np.concatenate([pos, neg]))
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(order) + 1)
    # tie average via scipy-style rankdata
    from scipy.stats import rankdata
    combined = np.concatenate([pos, neg])
    ranks = rankdata(combined)
    r_pos = ranks[: len(pos)].sum()
    n1, n2 = len(pos), len(neg)
    u = r_pos - n1 * (n1 + 1) / 2
    return u / (n1 * n2)


def ols_residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return y - y_hat using ordinary least squares (intercept added)."""
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    mask = ~np.isnan(y) & np.all(~np.isnan(X), axis=1)
    y_out = np.full_like(y, np.nan)
    if mask.sum() < 10:
        return y_out
    X_m = np.c_[np.ones(mask.sum()), X[mask]]
    beta, *_ = np.linalg.lstsq(X_m, y[mask], rcond=None)
    y_pred = X_m @ beta
    y_out[mask] = y[mask] - y_pred
    return y_out


def save_df(df: pd.DataFrame, name: str):
    p = DATA_DIR / name
    df.to_csv(p, index=False)
    print(f"  wrote {p}")


# --- data load --------------------------------------------------------------
def load_master() -> pd.DataFrame:
    print(f"[load] {MASTER}")
    df = pd.read_csv(MASTER, sep="\t", low_memory=False)
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("None")
    df["LOH_Subtype"] = pd.Categorical(df["LOH_Subtype"], categories=LOH_ORDER, ordered=True)
    df["AF_class"] = pd.Categorical(df["AF_class"], categories=AF_CLASS_ORDER, ordered=True)
    df["Coverage_Category"] = pd.Categorical(df["Coverage_Category"], categories=COV_ORDER, ordered=True)
    df["sample"] = pd.Categorical(df["sample"], categories=SAMPLE_ORDER, ordered=True)
    df["LOH_ord"] = df["LOH_Subtype"].map(LOH_ORD_MAP).astype(float)
    print(f"       rows={len(df):,}  tp={int((df.tp_label==1).sum()):,}  fp={int((df.tp_label==0).sum()):,}")
    return df


# --- Step 1: per-class TP rate + Wilson CI ---------------------------------
def step1_tp_rate(df: pd.DataFrame) -> pd.DataFrame:
    print("[step1] per-class TP rate + Wilson CI")
    rows = []
    for (sample, mode, sub), g in df.groupby(["sample", "mode", "LOH_Subtype"], observed=True):
        n = len(g)
        k = int((g.tp_label == 1).sum())
        p, lo, hi = wilson_ci(k, n)
        rows.append({
            "sample": str(sample), "mode": mode, "LOH_Subtype": str(sub),
            "n": n, "n_tp": k, "tp_rate": p, "ci_lo": lo, "ci_hi": hi,
        })
    out = pd.DataFrame(rows)
    save_df(out, "step1_tp_rate_per_class.csv")

    # bar chart: 2 cols (paired_full, to_pileup) x 7 samples
    modes = ["paired_full", "to_pileup"]
    fig, axes = plt.subplots(len(SAMPLE_ORDER), 2, figsize=(10, 14), sharex=True, sharey=True)
    colors = {"None": "#888888", "LOH_Noise": "#b15928",
              "LOH_Weak": "#fdbf6f", "LOH_Strong": "#ff7f00", "LOH_Subclone": "#6a3d9a"}
    for i, sample in enumerate(SAMPLE_ORDER):
        for j, mode in enumerate(modes):
            ax = axes[i, j]
            sub = out[(out["sample"] == sample) & (out["mode"] == mode)]
            if sub.empty:
                ax.axis("off"); continue
            sub = sub.set_index("LOH_Subtype").reindex(LOH_ORDER).reset_index()
            x = np.arange(len(LOH_ORDER))
            heights = sub["tp_rate"].values
            err_lo = np.clip(heights - sub["ci_lo"].values, 0, None)
            err_hi = np.clip(sub["ci_hi"].values - heights, 0, None)
            ax.bar(x, heights, color=[colors[k] for k in LOH_ORDER],
                   yerr=[np.where(np.isnan(err_lo), 0, err_lo),
                         np.where(np.isnan(err_hi), 0, err_hi)],
                   capsize=2, ecolor="black")
            ax.axhline(0.5, color="gray", lw=0.6, ls="--")
            ax.set_ylim(0, 1.05)
            if i == 0:
                ax.set_title(mode, fontsize=11)
            if j == 0:
                ax.set_ylabel(sample, fontsize=9)
            # n annotations
            for xi, nval in zip(x, sub["n"].values):
                ax.text(xi, 1.02, f"{int(nval) if not np.isnan(nval) else 0}",
                        ha="center", va="bottom", fontsize=6, rotation=0)
    for ax in axes[-1]:
        ax.set_xticks(np.arange(len(LOH_ORDER)))
        ax.set_xticklabels(LOH_ORDER, rotation=40, ha="right", fontsize=8)
    fig.suptitle("Step 1: LOH_Subtype TP rate + Wilson 95% CI (per sample × mode)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_path = FIG_DIR / "fig01_lohsubtype_tp_rate_bars.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")
    return out


# --- Step 2: LOH × AF × CN heatmap (TP rate + n) ---------------------------
def step2_heatmap(df: pd.DataFrame) -> pd.DataFrame:
    print("[step2] LOH_Subtype × AF_class × Coverage_Category heatmap")
    rows = []
    for (sub, af, cov), g in df.groupby(["LOH_Subtype", "AF_class", "Coverage_Category"], observed=True):
        n = len(g)
        k = int((g.tp_label == 1).sum())
        p, _, _ = wilson_ci(k, n)
        rows.append({
            "LOH_Subtype": str(sub), "AF_class": str(af), "Coverage_Category": str(cov),
            "n": n, "n_tp": k, "tp_rate": p,
        })
    out = pd.DataFrame(rows)
    save_df(out, "step2_loh_af_cov_cells.csv")

    # Heatmap: rows = LOH x AF (5x3=15), cols = Coverage_Category (6)
    row_index = pd.MultiIndex.from_product([LOH_ORDER, AF_CLASS_ORDER],
                                           names=["LOH_Subtype", "AF_class"])
    tp_mat = np.full((len(row_index), len(COV_ORDER)), np.nan)
    n_mat = np.zeros_like(tp_mat)
    for (sub, af), g in out.groupby(["LOH_Subtype", "AF_class"], observed=True):
        if (sub, af) not in row_index:
            continue
        ri = row_index.get_loc((sub, af))
        for _, r in g.iterrows():
            if r["Coverage_Category"] in COV_ORDER:
                ci = COV_ORDER.index(r["Coverage_Category"])
                n_mat[ri, ci] = r["n"]
                if r["n"] >= MIN_N:
                    tp_mat[ri, ci] = r["tp_rate"]

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(tp_mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(COV_ORDER)))
    ax.set_xticklabels(COV_ORDER, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(row_index)))
    ax.set_yticklabels([f"{a} | {b}" for a, b in row_index], fontsize=8)
    ax.set_xlabel("Coverage_Category (CN proxy)")
    ax.set_ylabel("LOH_Subtype × AF_class")
    ax.set_title("Step 2: TP rate per (LOH × AF × CN) cell  (cells with n<20 blank)")
    # n annotations
    for i in range(tp_mat.shape[0]):
        for j in range(tp_mat.shape[1]):
            if n_mat[i, j] > 0:
                color = "black" if not np.isnan(tp_mat[i, j]) and 0.3 < tp_mat[i, j] < 0.7 else "white"
                if np.isnan(tp_mat[i, j]):
                    color = "#888"
                ax.text(j, i, f"{int(n_mat[i,j])}", ha="center", va="center",
                        fontsize=6, color=color)
    fig.colorbar(im, ax=ax, label="TP rate")
    fig.tight_layout()
    out_path = FIG_DIR / "fig02_loh_af_cn_heatmap.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")
    return out


# --- Step 3: per-sample LOH_Subtype grid ------------------------------------
def step3_per_sample(df: pd.DataFrame, step1_df: pd.DataFrame):
    print("[step3] per-sample LOH_Subtype consistency grid")
    # heatmap: sample (rows) x LOH_Subtype (cols), split paired/TO
    modes = ["paired_full", "to_pileup"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, mode in zip(axes, modes):
        mat = np.full((len(SAMPLE_ORDER), len(LOH_ORDER)), np.nan)
        n_mat = np.zeros_like(mat)
        sub = step1_df[step1_df["mode"] == mode]
        for _, r in sub.iterrows():
            if r["sample"] in SAMPLE_ORDER and r["LOH_Subtype"] in LOH_ORDER and r["n"] >= MIN_N:
                ri = SAMPLE_ORDER.index(r["sample"])
                ci = LOH_ORDER.index(r["LOH_Subtype"])
                mat[ri, ci] = r["tp_rate"]
                n_mat[ri, ci] = r["n"]
        im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(np.arange(len(LOH_ORDER)))
        ax.set_xticklabels(LOH_ORDER, rotation=35, ha="right", fontsize=8)
        ax.set_yticks(np.arange(len(SAMPLE_ORDER)))
        ax.set_yticklabels(SAMPLE_ORDER)
        ax.set_title(f"{mode}: TP rate")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f"{mat[i,j]:.2f}\n(n={int(n_mat[i,j])})",
                            ha="center", va="center", fontsize=6, color="black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Step 3: per-sample LOH_Subtype TP rate consistency")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = FIG_DIR / "fig03_per_sample_grid.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")


# --- Step 4: ordinal AUC + Spearman ----------------------------------------
def step4_ordinal_auc(df: pd.DataFrame) -> pd.DataFrame:
    print("[step4] ordinal AUC + Spearman (LOH_ord vs tp_label)")
    from scipy.stats import spearmanr
    rows = []
    for (sample, mode), g in df.groupby(["sample", "mode"], observed=True):
        a = auc_binary(g.tp_label.values, g.LOH_ord.values)
        rho, p = spearmanr(g.tp_label.values, g.LOH_ord.values, nan_policy="omit")
        rows.append({"sample": str(sample), "mode": mode, "n": len(g),
                     "auc_ordinal": a, "spearman_rho": rho, "spearman_p": p})
    # Also Potential_LOH AUC (binary)
    rows_pl = []
    for (sample, mode), g in df.groupby(["sample", "mode"], observed=True):
        a = auc_binary(g.tp_label.values, g.Potential_LOH.astype(float).values)
        rows_pl.append({"sample": str(sample), "mode": mode, "n": len(g),
                        "auc_potential_loh": a})
    out = pd.DataFrame(rows).merge(pd.DataFrame(rows_pl), on=["sample", "mode", "n"])
    save_df(out, "step4_ordinal_auc.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    modes = ["paired_full", "to_pileup"]
    for ax, mode in zip(axes, modes):
        sub = out[out["mode"] == mode].set_index("sample").reindex(SAMPLE_ORDER)
        x = np.arange(len(SAMPLE_ORDER))
        ax.bar(x - 0.2, sub["auc_ordinal"].values, width=0.4, label="LOH_ord AUC",
               color="#1f78b4")
        ax.bar(x + 0.2, sub["auc_potential_loh"].values, width=0.4,
               label="Potential_LOH AUC", color="#33a02c")
        ax.axhline(0.5, color="gray", lw=0.6, ls="--", label="random")
        ax.axhline(0.58, color="red", lw=0.6, ls=":", label="Beyond-AUC ceiling")
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLE_ORDER, rotation=40, ha="right", fontsize=8)
        ax.set_ylim(0.3, 0.9)
        ax.set_title(mode)
        ax.set_ylabel("AUC")
        ax.legend(fontsize=8)
    fig.suptitle("Step 4: LOH ordinal / Potential_LOH AUC per sample × mode")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = FIG_DIR / "fig04_ordinal_auc.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")
    return out


# --- Step 5: residualize LOH_ord on AF + Coverage_Multiple ------------------
def step5_residualized(df: pd.DataFrame) -> pd.DataFrame:
    print("[step5] residualized AUC (within-sample OLS on vcf_AF + CovM)")
    rows = []
    for (sample, mode), g in df.groupby(["sample", "mode"], observed=True):
        g = g.copy()
        X = g[["vcf_AF", "Coverage_Multiple"]].values
        loh_resid = ols_residualize(g["LOH_ord"].values, X)
        pl_resid = ols_residualize(g["Potential_LOH"].astype(float).values, X)
        a_loh = auc_binary(g.tp_label.values, loh_resid)
        a_pl = auc_binary(g.tp_label.values, pl_resid)
        a_raw_loh = auc_binary(g.tp_label.values, g["LOH_ord"].values)
        a_raw_pl = auc_binary(g.tp_label.values, g["Potential_LOH"].astype(float).values)
        rows.append({
            "sample": str(sample), "mode": mode, "n": len(g),
            "auc_loh_raw": a_raw_loh, "auc_loh_resid": a_loh,
            "auc_pl_raw": a_raw_pl, "auc_pl_resid": a_pl,
        })
    out = pd.DataFrame(rows)
    save_df(out, "step5_residualized_auc.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, mode in zip(axes, ["paired_full", "to_pileup"]):
        sub = out[out["mode"] == mode].set_index("sample").reindex(SAMPLE_ORDER)
        x = np.arange(len(SAMPLE_ORDER))
        ax.plot(x, sub["auc_loh_raw"].values, "o-", label="LOH_ord raw", color="#1f78b4")
        ax.plot(x, sub["auc_loh_resid"].values, "s--", label="LOH_ord resid", color="#a6cee3")
        ax.plot(x, sub["auc_pl_raw"].values, "o-", label="Potential_LOH raw", color="#33a02c")
        ax.plot(x, sub["auc_pl_resid"].values, "s--", label="Potential_LOH resid", color="#b2df8a")
        ax.axhline(0.5, color="gray", lw=0.6, ls="--")
        ax.axhline(0.58, color="red", lw=0.6, ls=":")
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLE_ORDER, rotation=40, ha="right", fontsize=8)
        ax.set_ylim(0.3, 0.9)
        ax.set_title(mode)
        ax.set_ylabel("AUC")
        ax.legend(fontsize=7)
    fig.suptitle("Step 5: AUC raw vs residualized on (vcf_AF, Coverage_Multiple)")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = FIG_DIR / "fig05_residualized_auc.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")
    return out


# --- Step 6: per-chromosome LOH_Subtype distribution ------------------------
def step6_per_chr(df: pd.DataFrame) -> pd.DataFrame:
    print("[step6] per-chromosome LOH_Subtype distribution")
    # Aggregate: chr × LOH_Subtype frequency per sample (focus on paired_full
    # where Subtype definitions are cleanest)
    chr_order = [f"chr{i}" for i in list(range(1, 23))] + ["chrX", "chrY"]
    rows = []
    for (sample, mode, ch), g in df.groupby(["sample", "mode", "Chr"], observed=True):
        cnts = g["LOH_Subtype"].value_counts(dropna=False)
        tot = len(g)
        tp_rate = (g["tp_label"] == 1).mean() if tot > 0 else np.nan
        row = {"sample": str(sample), "mode": mode, "Chr": str(ch), "n": tot,
               "tp_rate": tp_rate}
        for lc in LOH_ORDER:
            row[f"frac_{lc}"] = cnts.get(lc, 0) / tot if tot > 0 else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    save_df(out, "step6_per_chr_loh.csv")

    # Plot focus: fraction of LOH_Subclone + LOH_Strong per chr (paired_full)
    fig, axes = plt.subplots(len(SAMPLE_ORDER), 1, figsize=(12, 14), sharex=True)
    mode = "paired_full"
    for ax, sample in zip(axes, SAMPLE_ORDER):
        sub = out[(out["sample"] == sample) & (out["mode"] == mode)]
        sub = sub.set_index("Chr").reindex(chr_order).reset_index()
        x = np.arange(len(chr_order))
        bottom = np.zeros(len(x))
        colors = {"None": "#cccccc", "LOH_Noise": "#b15928",
                  "LOH_Weak": "#fdbf6f", "LOH_Strong": "#ff7f00", "LOH_Subclone": "#6a3d9a"}
        for lc in LOH_ORDER:
            vals = sub[f"frac_{lc}"].fillna(0).values
            ax.bar(x, vals, bottom=bottom, color=colors[lc], label=lc if sample == SAMPLE_ORDER[0] else "")
            bottom += vals
        ax.set_ylim(0, 1.02)
        ax.set_ylabel(sample, fontsize=9)
        ax.set_yticks([0, 0.5, 1.0])
    axes[-1].set_xticks(np.arange(len(chr_order)))
    axes[-1].set_xticklabels(chr_order, rotation=60, fontsize=7)
    axes[0].legend(loc="upper right", fontsize=7, ncol=5)
    fig.suptitle("Step 6: per-chromosome LOH_Subtype fraction (paired_full)")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_path = FIG_DIR / "fig06_per_chr_distribution.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")
    return out


# --- main -------------------------------------------------------------------
def main():
    df = load_master()
    s1 = step1_tp_rate(df)
    _ = step2_heatmap(df)
    step3_per_sample(df, s1)
    _ = step4_ordinal_auc(df)
    _ = step5_residualized(df)
    _ = step6_per_chr(df)
    print("[done] G2 analysis complete")


if __name__ == "__main__":
    main()
