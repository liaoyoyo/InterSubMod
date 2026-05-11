#!/usr/bin/env python3
"""
Phase C · G1 — Read Counting & Coverage observation.

Scope:
  7 features: NumReads, NumCpGs, NTumorReads, NNormalReads,
              Coverage_Multiple, Diploid_Coverage_Used, Coverage_Category
  Stratification: LOH_Subtype (5) x AF_class (3, from vcf_AF) x cn_tier_F (CovM-derived, 5)
  TP/FP contrast, 7 samples, 2 modes (paired_full, to_pileup).

Steps (from 02_methodology.md):
  1. Global distribution (violin + KDE) + AUC, Cohen's d, MW-U
  2. LOH x AF x CN 32-cell heatmaps (TP rate, feat_TP, feat_FP, delta)
  3. Per-sample consistency grid + Spearman concordance
  4. Stratified AUC bar (Global / per-LOH / per-AF / per-CN)
  5. Confound guard: within-cell OLS residualization on (NumReads, vcf_AF, Coverage_Multiple)
  6. Spatial autocorrelation (chr+pos 5Mb bin)

Output:
  - figures/G1_coverage/fig01..fig06 per-feature (+ group combined heatmap)
  - data/G1_auc_table.tsv, data/G1_cell_delta.tsv
"""
from __future__ import annotations

import gzip
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
from sklearn.linear_model import LinearRegression

warnings.filterwarnings("ignore")

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/feature_layered_observation")
DATA_TSV = ROOT / "data" / "merged_with_vcf.tsv.gz"
FIG_DIR = ROOT / "figures" / "G1_coverage"
DATA_DIR = ROOT / "data"
FIG_DIR.mkdir(parents=True, exist_ok=True)

FEATURES_NUM = ["NumReads", "NumCpGs", "NTumorReads", "NNormalReads",
                "Coverage_Multiple", "Diploid_Coverage_Used"]
FEATURES_CAT = ["Coverage_Category"]
ALL_FEATURES = FEATURES_NUM + FEATURES_CAT

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]

LOH_ORDER = ["LOH_None", "LOH_Weak", "LOH_Noise", "LOH_Strong", "LOH_Subclone"]
AF_ORDER = ["Extreme", "Near-half", "Intermediate"]
COV_CAT_ORDER = ["Low", "Normal", "Elevated", "CNV_Gain", "CNV_Loss", "High_Copy"]

# CN tier (CovM) thresholds: 0.65 / 0.99 / 1.33 / 1.82
CN_BREAKS = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
CN_LABELS = ["CN_Loss", "CN_Near1", "CN_Diploid", "CN_Gain", "CN_HighGain"]


# --------------------------------------------------------------------------- #
# I/O                                                                         #
# --------------------------------------------------------------------------- #
def load() -> pd.DataFrame:
    print(f"[load] reading {DATA_TSV}", flush=True)
    usecols = ["sample", "mode", "Chr", "Pos", "tp_label",
               "LOH_Subtype", "AF_class", "vcf_AF", "AF",
               "Coverage_Multiple", "Coverage_Category",
               "NumReads", "NumCpGs", "NTumorReads", "NNormalReads",
               "Diploid_Coverage_Used"]
    df = pd.read_csv(DATA_TSV, sep="\t", usecols=usecols, low_memory=False)
    print(f"[load] rows={len(df):,} cols={df.shape[1]}", flush=True)

    # cn_tier_F from Coverage_Multiple (master thresholds)
    df["cn_tier_F"] = pd.cut(df["Coverage_Multiple"].astype(float),
                             bins=CN_BREAKS, labels=CN_LABELS, right=True)

    # Ensure categorical ordering for strat keys
    df["LOH_Subtype"] = pd.Categorical(df["LOH_Subtype"].fillna("LOH_None"),
                                       categories=LOH_ORDER, ordered=True)
    df["AF_class"] = pd.Categorical(df["AF_class"].fillna("Unknown"),
                                    categories=AF_ORDER + ["Unknown"], ordered=False)
    df["tp_label"] = df["tp_label"].astype(int)
    return df


# --------------------------------------------------------------------------- #
# Stats helpers                                                               #
# --------------------------------------------------------------------------- #
def auc_wilson_ci(y_true: np.ndarray, scores: np.ndarray) -> tuple[float, float, float, int, int]:
    """Return (auc, lo, hi, n_pos, n_neg) with Wilson-style (Hanley-McNeil) CI."""
    mask = ~np.isnan(scores)
    y = y_true[mask]
    s = scores[mask]
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos < 5 or n_neg < 5:
        return np.nan, np.nan, np.nan, n_pos, n_neg
    # Mann-Whitney-based AUC (no sklearn dep)
    order = np.argsort(s, kind="mergesort")
    y_o = y[order]
    ranks = np.empty_like(s, dtype=float)
    # Use scipy rankdata on scores (average rank)
    rnks = stats.rankdata(s)
    sum_rank_pos = rnks[y == 1].sum()
    auc = (sum_rank_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    # Hanley-McNeil SE (conservative)
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    se = np.sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc**2) + (n_neg - 1) * (q2 - auc**2))
                 / (n_pos * n_neg))
    lo = max(0.0, auc - 1.96 * se)
    hi = min(1.0, auc + 1.96 * se)
    return auc, lo, hi, n_pos, n_neg


def cohen_d(x: np.ndarray, y: np.ndarray) -> float:
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) < 2 or len(y) < 2:
        return np.nan
    nx, ny = len(x), len(y)
    sx, sy = x.std(ddof=1), y.std(ddof=1)
    sp = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    if sp == 0:
        return 0.0
    return (x.mean() - y.mean()) / sp


def mwu_p(x: np.ndarray, y: np.ndarray) -> float:
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) < 5 or len(y) < 5:
        return np.nan
    try:
        return float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except Exception:
        return np.nan


def cat_auc(y_true: np.ndarray, cat: pd.Series) -> float:
    """One-hot aggregate: TP fraction per category used as score."""
    df = pd.DataFrame({"y": y_true, "c": cat})
    rate = df.groupby("c")["y"].mean()
    scores = df["c"].map(rate).to_numpy(dtype=float)
    auc, *_ = auc_wilson_ci(y_true, scores)
    return auc


# --------------------------------------------------------------------------- #
# Step 1 — Global distribution                                                #
# --------------------------------------------------------------------------- #
def step1_global_panel(df: pd.DataFrame, feat_list: list[str], out_png: Path):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    axes = axes.flatten()
    summary_rows = []
    for i, feat in enumerate(feat_list):
        ax = axes[i]
        if feat in FEATURES_CAT:
            # stacked bar of TP rate per category
            if feat == "Coverage_Category":
                order = [c for c in COV_CAT_ORDER if c in df[feat].dropna().unique()]
                tp_rate = df.groupby(feat)["tp_label"].mean().reindex(order)
                n_per = df.groupby(feat).size().reindex(order)
                ax.bar(range(len(order)), tp_rate.values, color="#2E86AB")
                ax.axhline(0.5, ls="--", c="gray", lw=0.7)
                ax.set_xticks(range(len(order)))
                ax.set_xticklabels(order, rotation=30, ha="right", fontsize=8)
                for j, (v, n) in enumerate(zip(tp_rate.values, n_per.values)):
                    ax.text(j, v + 0.01, f"n={n:,}", ha="center", fontsize=7)
                ax.set_ylim(0, 1.05)
                auc = cat_auc(df["tp_label"].to_numpy(), df[feat])
                ax.set_title(f"{feat}\nAUC={auc:.3f}  TP rate by category", fontsize=9)
                ax.set_ylabel("TP rate")
                summary_rows.append({"feature": feat, "auc": auc, "lo": np.nan, "hi": np.nan,
                                     "cohen_d": np.nan, "mwu_p": np.nan,
                                     "mean_tp": np.nan, "mean_fp": np.nan,
                                     "n_tp": int((df["tp_label"] == 1).sum()),
                                     "n_fp": int((df["tp_label"] == 0).sum())})
            continue
        vals = df[feat].astype(float).to_numpy()
        tp = vals[df["tp_label"] == 1]
        fp = vals[df["tp_label"] == 0]
        tp = tp[~np.isnan(tp)]
        fp = fp[~np.isnan(fp)]
        # log-scale for counts-like features
        log_scale = feat in {"NumReads", "NumCpGs", "NTumorReads", "NNormalReads"}
        plot_tp = np.log10(tp + 1) if log_scale else tp
        plot_fp = np.log10(fp + 1) if log_scale else fp
        parts = ax.violinplot([plot_fp, plot_tp], positions=[0, 1], widths=0.7,
                              showmeans=False, showmedians=True)
        for pc, color in zip(parts["bodies"], ["#E07A5F", "#2E86AB"]):
            pc.set_facecolor(color)
            pc.set_alpha(0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"], fontsize=8)
        auc, lo, hi, _, _ = auc_wilson_ci(df["tp_label"].to_numpy(), vals)
        d = cohen_d(tp, fp)
        p = mwu_p(tp, fp)
        ax.set_title(f"{feat}{' (log10)' if log_scale else ''}\nAUC={auc:.3f} "
                     f"[{lo:.3f},{hi:.3f}]  d={d:.2f}", fontsize=9)
        summary_rows.append({"feature": feat, "auc": auc, "lo": lo, "hi": hi,
                             "cohen_d": d, "mwu_p": p,
                             "mean_tp": float(np.nanmean(tp)) if len(tp) else np.nan,
                             "mean_fp": float(np.nanmean(fp)) if len(fp) else np.nan,
                             "n_tp": int(len(tp)), "n_fp": int(len(fp))})
    # hide any unused axes
    for i in range(len(feat_list), len(axes)):
        axes[i].axis("off")
    fig.suptitle("G1 Coverage · Step 1 Global TP vs FP distribution (all samples, both modes)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return pd.DataFrame(summary_rows)


# --------------------------------------------------------------------------- #
# Step 2 — LOH × AF × CN heatmap (feature-wise)                              #
# --------------------------------------------------------------------------- #
def step2_heatmap(df: pd.DataFrame, feat: str, out_stub: Path, min_n: int = 20):
    df = df.copy()
    # rows: LOH x AF ; cols: CN tier
    df["row_key"] = df["LOH_Subtype"].astype(str) + " | " + df["AF_class"].astype(str)
    row_order = [f"{l} | {a}" for l in LOH_ORDER for a in AF_ORDER]
    col_order = CN_LABELS

    piv_n = df.pivot_table(index="row_key", columns="cn_tier_F",
                           values="tp_label", aggfunc="size").reindex(row_order, axis=0) \
              .reindex(col_order, axis=1)
    piv_rate = df.pivot_table(index="row_key", columns="cn_tier_F",
                              values="tp_label", aggfunc="mean").reindex(row_order, axis=0) \
                 .reindex(col_order, axis=1)

    # Only use cells with n >= min_n
    mask_ok = piv_n >= min_n

    if feat in FEATURES_CAT:
        # use numeric encoding: TP rate in the category as a proxy
        return None
    vals = df[feat].astype(float)
    piv_tp = df[df["tp_label"] == 1].pivot_table(
        index="row_key", columns="cn_tier_F", values=feat, aggfunc="median"
    ).reindex(row_order, axis=0).reindex(col_order, axis=1)
    piv_fp = df[df["tp_label"] == 0].pivot_table(
        index="row_key", columns="cn_tier_F", values=feat, aggfunc="median"
    ).reindex(row_order, axis=0).reindex(col_order, axis=1)
    piv_delta = (piv_tp - piv_fp).where(mask_ok)

    fig, axes = plt.subplots(1, 4, figsize=(24, 11))
    # A: TP rate
    im0 = axes[0].imshow(piv_rate.where(mask_ok).values, aspect="auto",
                         cmap="RdBu_r", vmin=0, vmax=1)
    axes[0].set_title(f"A · TP rate (n≥{min_n})")
    fig.colorbar(im0, ax=axes[0], fraction=0.03)
    # B: feature median TP
    v_tp = piv_tp.where(mask_ok)
    im1 = axes[1].imshow(v_tp.values, aspect="auto", cmap="Blues")
    axes[1].set_title(f"B · median({feat}) | TP only")
    fig.colorbar(im1, ax=axes[1], fraction=0.03)
    # C: feature median FP
    v_fp = piv_fp.where(mask_ok)
    im2 = axes[2].imshow(v_fp.values, aspect="auto", cmap="Oranges")
    axes[2].set_title(f"C · median({feat}) | FP only")
    fig.colorbar(im2, ax=axes[2], fraction=0.03)
    # D: delta (TP - FP)
    dlt = piv_delta.values.astype(float)
    finite = dlt[np.isfinite(dlt)]
    if len(finite):
        vlim = max(abs(np.nanpercentile(finite, 5)), abs(np.nanpercentile(finite, 95)), 1e-6)
        norm = TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim)
        im3 = axes[3].imshow(dlt, aspect="auto", cmap="coolwarm", norm=norm)
    else:
        im3 = axes[3].imshow(dlt, aspect="auto", cmap="coolwarm")
    axes[3].set_title(f"D · Δ(TP−FP) median")
    fig.colorbar(im3, ax=axes[3], fraction=0.03)

    for ax in axes:
        ax.set_yticks(range(len(row_order)))
        ax.set_yticklabels(row_order, fontsize=7)
        ax.set_xticks(range(len(col_order)))
        ax.set_xticklabels(col_order, rotation=30, ha="right", fontsize=8)
    # annotate n on TP rate panel
    for i, r in enumerate(row_order):
        for j, c in enumerate(col_order):
            n = piv_n.loc[r, c] if (r in piv_n.index and c in piv_n.columns) else 0
            if pd.isna(n):
                n = 0
            if n >= min_n:
                axes[0].text(j, i, f"{int(n):,}", ha="center", va="center",
                             fontsize=6, color="black")
    fig.suptitle(f"{feat} · Step 2 LOH × AF × CN heatmap (all samples pooled)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_png = out_stub.with_name(f"{feat}_fig02_heatmap.png")
    fig.savefig(out_png, dpi=130)
    plt.close(fig)

    # emit flat delta table
    flat = piv_delta.reset_index().melt(id_vars="row_key", var_name="cn_tier",
                                        value_name="delta_tp_fp")
    flat["n"] = flat.apply(lambda r: piv_n.loc[r["row_key"], r["cn_tier"]]
                           if r["row_key"] in piv_n.index and r["cn_tier"] in piv_n.columns
                           else 0, axis=1)
    flat["feature"] = feat
    return flat


# --------------------------------------------------------------------------- #
# Step 3 — Per-sample consistency                                             #
# --------------------------------------------------------------------------- #
def step3_per_sample(df: pd.DataFrame, feat: str, out_stub: Path, min_n: int = 10):
    if feat in FEATURES_CAT:
        return None
    samples = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    n_samp = len(samples)
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()
    per_sample_rate = {}
    for i, s in enumerate(samples):
        sub = df[df["sample"] == s]
        piv = sub.pivot_table(index="LOH_Subtype", columns="AF_class",
                              values="tp_label", aggfunc="mean") \
                 .reindex(LOH_ORDER, axis=0).reindex(AF_ORDER, axis=1)
        piv_n = sub.pivot_table(index="LOH_Subtype", columns="AF_class",
                                values="tp_label", aggfunc="size") \
                   .reindex(LOH_ORDER, axis=0).reindex(AF_ORDER, axis=1)
        mask = piv_n >= min_n
        im = axes[i].imshow(piv.where(mask).values, aspect="auto", cmap="RdBu_r",
                            vmin=0, vmax=1)
        axes[i].set_title(f"{s}\n(n={len(sub):,})", fontsize=9)
        axes[i].set_xticks(range(len(AF_ORDER)))
        axes[i].set_xticklabels(AF_ORDER, rotation=30, ha="right", fontsize=7)
        axes[i].set_yticks(range(len(LOH_ORDER)))
        axes[i].set_yticklabels(LOH_ORDER, fontsize=7)
        for ii, r in enumerate(LOH_ORDER):
            for jj, c in enumerate(AF_ORDER):
                n = piv_n.loc[r, c] if r in piv_n.index and c in piv_n.columns else 0
                if pd.isna(n):
                    n = 0
                if n >= min_n:
                    v = piv.loc[r, c]
                    axes[i].text(jj, ii, f"{v:.2f}\nn={int(n)}", ha="center",
                                 va="center", fontsize=6)
        per_sample_rate[s] = piv.values.flatten()
    # concordance subplot
    ax_conc = axes[-1]
    mat = np.full((n_samp, n_samp), np.nan)
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            v1 = per_sample_rate[s1]
            v2 = per_sample_rate[s2]
            m = ~(np.isnan(v1) | np.isnan(v2))
            if m.sum() >= 3:
                mat[i, j] = stats.spearmanr(v1[m], v2[m]).correlation
    im2 = ax_conc.imshow(mat, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax_conc.set_title("Cross-sample Spearman concordance\n(TP rate per LOH×AF cell)",
                      fontsize=9)
    ax_conc.set_xticks(range(n_samp))
    ax_conc.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax_conc.set_yticks(range(n_samp))
    ax_conc.set_yticklabels(samples, fontsize=7)
    fig.colorbar(im2, ax=ax_conc, fraction=0.04)
    fig.suptitle(f"{feat} · Step 3 Per-sample TP rate consistency (LOH × AF)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_stub.with_name(f"{feat}_fig03_per_sample.png"), dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Step 4 — Stratified AUC panel                                               #
# --------------------------------------------------------------------------- #
def step4_stratified(df: pd.DataFrame, feat: str, out_stub: Path):
    if feat in FEATURES_CAT:
        return None
    y = df["tp_label"].to_numpy()
    s = df[feat].astype(float).to_numpy()
    rows = [("global", "all", *auc_wilson_ci(y, s))]
    for loh in LOH_ORDER:
        m = df["LOH_Subtype"] == loh
        if m.sum() >= 100:
            rows.append(("LOH", loh, *auc_wilson_ci(y[m], s[m])))
    for af in AF_ORDER:
        m = df["AF_class"] == af
        if m.sum() >= 100:
            rows.append(("AF", af, *auc_wilson_ci(y[m], s[m])))
    for cn in CN_LABELS:
        m = df["cn_tier_F"] == cn
        if m.sum() >= 100:
            rows.append(("CN", cn, *auc_wilson_ci(y[m], s[m])))
    dfres = pd.DataFrame(rows, columns=["layer", "group", "auc", "lo", "hi", "n_pos", "n_neg"])
    dfres["feature"] = feat

    fig, ax = plt.subplots(figsize=(11, 5.5))
    xpos = np.arange(len(dfres))
    colors = {"global": "#333333", "LOH": "#2E86AB", "AF": "#F4A261", "CN": "#6A994E"}
    for i, (_, row) in enumerate(dfres.iterrows()):
        yerr = [[max(0, row["auc"] - row["lo"])], [max(0, row["hi"] - row["auc"])]]
        ax.bar(i, row["auc"], color=colors.get(row["layer"], "gray"),
               yerr=yerr, capsize=3, edgecolor="black", linewidth=0.4)
        ax.text(i, row["auc"] + 0.01,
                f"{row['auc']:.3f}\nn={int(row['n_pos'])+int(row['n_neg']):,}",
                ha="center", fontsize=6)
    ax.axhline(0.5, ls="--", c="gray", lw=0.8, label="random")
    ax.axhline(0.58, ls=":", c="red", lw=0.8, label="Beyond-AUC ceiling 0.58")
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"{r.layer}:{r.group}" for _, r in dfres.iterrows()],
                       rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0.3, 1.0)
    ax.set_ylabel("AUC")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"{feat} · Step 4 Stratified AUC (Global / LOH / AF / CN)")
    fig.tight_layout()
    fig.savefig(out_stub.with_name(f"{feat}_fig04_stratified_auc.png"), dpi=130)
    plt.close(fig)
    return dfres


# --------------------------------------------------------------------------- #
# Step 5 — Confound guard (OLS residualization)                              #
# --------------------------------------------------------------------------- #
def step5_confound(df: pd.DataFrame, feat: str, out_stub: Path):
    if feat in FEATURES_CAT:
        return None
    # Within-cell OLS: residualize feat on (NumReads, vcf_AF, Coverage_Multiple)
    covars = [c for c in ["NumReads", "vcf_AF", "Coverage_Multiple"] if c != feat]
    sub = df.dropna(subset=[feat, "vcf_AF", "NumReads", "Coverage_Multiple", "tp_label"]).copy()
    if len(sub) < 200:
        return None
    X = sub[covars].astype(float).to_numpy()
    y_feat = sub[feat].astype(float).to_numpy()
    y_tp = sub["tp_label"].to_numpy()
    try:
        lr = LinearRegression().fit(X, y_feat)
        resid = y_feat - lr.predict(X)
    except Exception:
        return None
    raw_auc = auc_wilson_ci(y_tp, y_feat)[0]
    res_auc = auc_wilson_ci(y_tp, resid)[0]

    # AF-bin cross
    af_rows = []
    for af in AF_ORDER:
        m = sub["AF_class"] == af
        if m.sum() >= 200:
            r_auc = auc_wilson_ci(y_tp[m], y_feat[m])[0]
            af_rows.append((af, r_auc, int(m.sum())))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # Panel A: raw vs residualized histogram by tp
    ax = axes[0]
    bins = 50
    ax.hist(y_feat[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
            density=True, label=f"raw FP")
    ax.hist(y_feat[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
            density=True, label=f"raw TP")
    ax.set_title(f"raw {feat}  AUC={raw_auc:.3f}")
    ax.legend(fontsize=8)
    ax.set_xlabel(feat)
    # Panel B: residualized
    ax2 = axes[1]
    ax2.hist(resid[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
             density=True, label="resid FP")
    ax2.hist(resid[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
             density=True, label="resid TP")
    ax2.set_title(f"resid. on {covars}  AUC={res_auc:.3f}\n"
                  f"AF-bin AUCs: " + ", ".join([f"{a}={v:.3f}" for a, v, _ in af_rows]))
    ax2.legend(fontsize=8)
    fig.suptitle(f"{feat} · Step 5 Confound guard (within-cell OLS residualization)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_stub.with_name(f"{feat}_fig05_confound.png"), dpi=130)
    plt.close(fig)
    return {"feature": feat, "raw_auc": raw_auc, "resid_auc": res_auc,
            "af_bin_aucs": {a: v for a, v, _ in af_rows}}


# --------------------------------------------------------------------------- #
# Step 6 — Spatial autocorrelation                                            #
# --------------------------------------------------------------------------- #
def step6_spatial(df: pd.DataFrame, feat: str, out_stub: Path, bin_mb: int = 5):
    if feat in FEATURES_CAT:
        return None
    d = df.dropna(subset=["Chr", "Pos", feat, "tp_label"]).copy()
    d["pos_bin"] = (d["Pos"].astype(int) // (bin_mb * 1_000_000))
    d["bin_key"] = d["Chr"].astype(str) + ":" + d["pos_bin"].astype(str)
    grp = d.groupby("bin_key")
    rows = []
    for key, g in grp:
        if len(g) < 50:
            continue
        y = g["tp_label"].to_numpy()
        if y.sum() < 5 or (len(y) - y.sum()) < 5:
            continue
        a, *_ = auc_wilson_ci(y, g[feat].astype(float).to_numpy())
        rows.append({"bin": key, "n": len(g), "tp_rate": y.mean(), "auc": a})
    if not rows:
        return None
    s = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    ax.scatter(s["tp_rate"], s["auc"], alpha=0.3, s=8)
    ax.axhline(0.5, ls="--", c="gray")
    ax.axhline(0.58, ls=":", c="red")
    ax.axvline(0.8, ls="--", c="purple", alpha=0.5)
    ax.set_xlabel("bin TP rate")
    ax.set_ylabel(f"bin-level AUC for {feat}")
    ax.set_title(f"{feat} · Step 6 Spatial auto-correlation ({bin_mb} Mb bins)")
    ax2 = axes[1]
    ax2.hist(s["auc"], bins=30, color="#2E86AB", alpha=0.7)
    ax2.axvline(0.5, ls="--", c="gray")
    ax2.axvline(0.58, ls=":", c="red")
    ax2.set_xlabel("per-bin AUC")
    ax2.set_ylabel("# bins")
    # artifact flag
    high_tp = s[s["tp_rate"] > 0.8]
    low_tp = s[s["tp_rate"] < 0.5]
    flag = ""
    if len(high_tp) >= 10 and len(low_tp) >= 10:
        hi_mean = high_tp["auc"].mean()
        lo_mean = low_tp["auc"].mean()
        if hi_mean - lo_mean > 0.08:
            flag = f"⚠ artifact suspect (ΔAUC={hi_mean - lo_mean:+.3f} hi_tp vs lo_tp)"
        else:
            flag = f"ok (ΔAUC={hi_mean - lo_mean:+.3f})"
    fig.suptitle(f"{feat} · spatial auto-correlation  {flag}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_stub.with_name(f"{feat}_fig06_spatial.png"), dpi=130)
    plt.close(fig)
    return {"feature": feat, "n_bins": len(s), "flag": flag}


# --------------------------------------------------------------------------- #
# Group combined: Coverage_Multiple KDE baseline per sample-mode              #
# --------------------------------------------------------------------------- #
def group_covm_kde(df: pd.DataFrame, out_png: Path):
    fig, axes = plt.subplots(2, 7, figsize=(24, 8), sharex=True, sharey=True)
    samples = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    for col, s in enumerate(samples):
        for row, mode in enumerate(["paired_full", "to_pileup"]):
            ax = axes[row, col]
            sub = df[(df["sample"] == s) & (df["mode"] == mode)]
            if len(sub) < 50:
                ax.axis("off")
                continue
            dc = sub["Diploid_Coverage_Used"].dropna().unique()
            dc_str = f"dc={dc[0]:.0f}" if len(dc) > 0 else "dc=?"
            for tp, color in [(0, "#E07A5F"), (1, "#2E86AB")]:
                vals = sub.loc[sub["tp_label"] == tp, "Coverage_Multiple"].dropna()
                if len(vals) < 20:
                    continue
                vals = vals[(vals >= 0) & (vals <= 6)]
                ax.hist(vals, bins=60, density=True, alpha=0.45,
                        color=color, label="TP" if tp == 1 else "FP")
            ax.axvline(1.0, ls=":", c="gray", lw=0.6)
            ax.axvline(0.65, ls=":", c="red", lw=0.5)
            ax.axvline(1.33, ls=":", c="red", lw=0.5)
            ax.axvline(1.82, ls=":", c="red", lw=0.5)
            ax.set_title(f"{s} · {mode}\n{dc_str} n={len(sub):,}", fontsize=8)
            if col == 0:
                ax.set_ylabel("density")
            if row == 1:
                ax.set_xlabel("Coverage_Multiple")
            if col == 0 and row == 0:
                ax.legend(fontsize=7)
    fig.suptitle("G1 group · Coverage_Multiple distribution per sample × mode "
                 "(baseline comparability)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Group combined: AUC matrix heatmap                                          #
# --------------------------------------------------------------------------- #
def group_auc_matrix(df: pd.DataFrame, out_png: Path):
    samples = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    rows = []
    for s in samples:
        for mode in ["paired_full", "to_pileup"]:
            sub = df[(df["sample"] == s) & (df["mode"] == mode)]
            if len(sub) < 100:
                continue
            y = sub["tp_label"].to_numpy()
            rec = {"sample": s, "mode": mode, "n": len(sub), "tp_rate": y.mean()}
            for feat in FEATURES_NUM:
                if feat in sub.columns:
                    a, *_ = auc_wilson_ci(y, sub[feat].astype(float).to_numpy())
                    rec[feat] = a
            rec["Coverage_Category"] = cat_auc(y, sub["Coverage_Category"])
            rows.append(rec)
    m = pd.DataFrame(rows)
    if m.empty:
        return None
    m["key"] = m["sample"] + " · " + m["mode"]
    # Order rows by sample order + mode
    key_order = [f"{s} · {md}" for s in samples for md in ["paired_full", "to_pileup"]]
    m = m.set_index("key").reindex([k for k in key_order if k in m["key"].values])
    mat = m[ALL_FEATURES].astype(float)
    fig, ax = plt.subplots(figsize=(10, max(6, 0.35 * len(mat))))
    im = ax.imshow(mat.values, aspect="auto", cmap="coolwarm", vmin=0.35, vmax=0.75)
    ax.set_xticks(range(len(ALL_FEATURES)))
    ax.set_xticklabels(ALL_FEATURES, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7,
                        color="white" if (v > 0.63 or v < 0.45) else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, label="per-sample per-mode AUC")
    ax.set_title("G1 per-feature AUC · 7 samples × 2 modes\n"
                 "(dashed 0.58 Beyond-AUC ceiling)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return m


# --------------------------------------------------------------------------- #
# Main driver                                                                 #
# --------------------------------------------------------------------------- #
def main():
    df = load()
    print("[main] producing group Step1 panel", flush=True)
    step1 = step1_global_panel(df, ALL_FEATURES, FIG_DIR / "fig01_group_global.png")

    auc_rows = []
    confound_rows = []
    cell_deltas = []

    for feat in ALL_FEATURES:
        print(f"[feat] {feat}", flush=True)
        stub = FIG_DIR / feat  # used to derive filenames
        if feat in FEATURES_NUM:
            flat = step2_heatmap(df, feat, stub)
            if flat is not None:
                cell_deltas.append(flat)
            step3_per_sample(df, feat, stub)
            r = step4_stratified(df, feat, stub)
            if r is not None:
                auc_rows.append(r)
            c = step5_confound(df, feat, stub)
            if c is not None:
                confound_rows.append(c)
            step6_spatial(df, feat, stub)
        else:
            # Coverage_Category categorical: keep only Step1+per-sample TP rate panel
            pass

    # Group combined figures
    print("[main] group Coverage_Multiple KDE", flush=True)
    group_covm_kde(df, FIG_DIR / "fig07_group_CovM_kde.png")
    print("[main] group AUC matrix", flush=True)
    group_auc_matrix(df, FIG_DIR / "fig08_group_auc_matrix.png")

    # Save tables
    step1.to_csv(DATA_DIR / "G1_global_stats.tsv", sep="\t", index=False)
    if auc_rows:
        big = pd.concat(auc_rows, ignore_index=True)
        big.to_csv(DATA_DIR / "G1_auc_table.tsv", sep="\t", index=False)
    if cell_deltas:
        big = pd.concat(cell_deltas, ignore_index=True)
        big.to_csv(DATA_DIR / "G1_cell_delta.tsv", sep="\t", index=False)
    if confound_rows:
        pd.DataFrame(confound_rows).to_csv(DATA_DIR / "G1_confound.tsv",
                                           sep="\t", index=False)
    print("[main] done. Outputs under:")
    print(" ", FIG_DIR)
    print(" ", DATA_DIR)


if __name__ == "__main__":
    main()
