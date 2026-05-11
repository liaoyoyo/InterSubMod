#!/usr/bin/env python3
"""Phase C · G5 — HP Merged (2-bucket) observation.

Features in scope (Group G5, 7 features):
  HP_Ratio, HP1FamilyN, HP2FamilyN, HPMergedDelta, HPMergedP, HPMergedSig,
  and per-group combined (GlobalP_HPFamily, CramersV_HPFamily).

Special handling (per task spec):
  - HP_Ratio: known extreme in LOH (~0.09 reported in Phase 1); stratify by
    Potential_LOH / LOH_Subtype.  Test normalization min(HP_Ratio, 1-HP_Ratio).
  - HPMergedDelta: high correlation with AlleleDelta expected (master AF is
    ~|AlleleDelta|). Compute Pearson/Spearman vs AlleleDelta and vs vcf_AF.
    Residualize on (AlleleDelta, NumReads, vcf_AF) and re-AUC.
    [L2 collider bias warning — noted in 02_methodology §Step 5.]
  - HPMergedSig: boolean; per-class TP rate + Wilson CI.
  - HP1FamilyN + HP2FamilyN: check if sum == NumReads (i.e. effectively a
    NumReads proxy) by computing (HP1FamilyN+HP2FamilyN)/NumReads ratio.

Input:  data/G5/master_g5.tsv.gz  (produced by g5_build_extended_master.py)
Output: figures/G5_hp_merged/fig01..fig08
        data/G5/G5_*.tsv
"""
from __future__ import annotations

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
DATA_TSV = ROOT / "data/G5/master_g5.tsv.gz"
FIG_DIR = ROOT / "figures/G5_hp_merged"
DATA_DIR = ROOT / "data/G5"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]
LOH_ORDER = ["LOH_None", "LOH_Weak", "LOH_Noise", "LOH_Strong", "LOH_Subclone"]
AF_ORDER = ["Extreme", "Near-half", "Intermediate"]
CN_BREAKS = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
CN_LABELS = ["CN_Loss", "CN_Near1", "CN_Diploid", "CN_Gain", "CN_HighGain"]

# Continuous / discrete feature sets (used by Step 1 global panel etc.)
FEATURES_CONT = ["HP_Ratio", "HP_Ratio_norm", "HPMergedDelta",
                 "HPMergedP", "HP1FamilyN", "HP2FamilyN", "HP_FamilyN_sum",
                 "HP_FamilyN_frac", "GlobalP_HPFamily", "CramersV_HPFamily"]
FEATURES_BIN = ["HPMergedSig"]

# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #
def auc_wilson_ci(y_true: np.ndarray, scores: np.ndarray):
    mask = ~np.isnan(scores)
    y = y_true[mask]
    s = scores[mask]
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos < 5 or n_neg < 5:
        return np.nan, np.nan, np.nan, n_pos, n_neg
    rnks = stats.rankdata(s)
    sum_rank_pos = rnks[y == 1].sum()
    auc = (sum_rank_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    se = np.sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc**2) +
                  (n_neg - 1) * (q2 - auc**2)) / (n_pos * n_neg))
    lo = max(0.0, auc - 1.96 * se)
    hi = min(1.0, auc + 1.96 * se)
    return auc, lo, hi, n_pos, n_neg


def cohen_d(x, y):
    x = np.asarray(x, dtype=float); x = x[~np.isnan(x)]
    y = np.asarray(y, dtype=float); y = y[~np.isnan(y)]
    if len(x) < 2 or len(y) < 2:
        return np.nan
    sp = np.sqrt(((len(x)-1)*x.std(ddof=1)**2 + (len(y)-1)*y.std(ddof=1)**2)
                 / (len(x)+len(y)-2))
    if sp == 0:
        return 0.0
    return (x.mean() - y.mean()) / sp


def mwu_p(x, y):
    x = np.asarray(x, dtype=float); x = x[~np.isnan(x)]
    y = np.asarray(y, dtype=float); y = y[~np.isnan(y)]
    if len(x) < 5 or len(y) < 5:
        return np.nan
    try:
        return float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except Exception:
        return np.nan


def wilson_ci_k_n(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    denom = 1 + z**2/n
    centre = (p + z**2/(2*n))/denom
    half = z*np.sqrt(p*(1-p)/n + z**2/(4*n**2))/denom
    return p, max(0., centre-half), min(1., centre+half)


# --------------------------------------------------------------------------- #
# Load & enrich                                                               #
# --------------------------------------------------------------------------- #
def load() -> pd.DataFrame:
    print(f"[load] {DATA_TSV}", flush=True)
    df = pd.read_csv(DATA_TSV, sep="\t", low_memory=False)
    print(f"[load] rows={len(df):,}  cols={df.shape[1]}", flush=True)

    # Ensure tp_label int
    df["tp_label"] = pd.to_numeric(df["tp_label"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["tp_label"]).copy()
    df["tp_label"] = df["tp_label"].astype(int)

    # Derive cn_tier_F
    df["cn_tier_F"] = pd.cut(df["Coverage_Multiple"].astype(float),
                             bins=CN_BREAKS, labels=CN_LABELS, right=True)
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("LOH_None")
    df["LOH_Subtype"] = pd.Categorical(df["LOH_Subtype"], categories=LOH_ORDER, ordered=True)
    df["AF_class"] = pd.Categorical(df["AF_class"].fillna("Unknown"),
                                    categories=AF_ORDER + ["Unknown"])
    # Normalized HP ratio: min distance to 0.5 (0 == balanced, 0.5 == extreme)
    df["HP_Ratio_norm"] = np.where(df["HP_Ratio"].notna(),
                                   (df["HP_Ratio"] - 0.5).abs(),
                                   np.nan)
    # Family count sanity metrics
    df["HP_FamilyN_sum"] = df["HP1FamilyN"].fillna(0) + df["HP2FamilyN"].fillna(0)
    df["HP_FamilyN_frac"] = df["HP_FamilyN_sum"] / df["NumReads"].replace(0, np.nan)
    return df


# --------------------------------------------------------------------------- #
# Step 0 — confound / relation sanity                                         #
# --------------------------------------------------------------------------- #
def step0_relations(df: pd.DataFrame):
    """Build relation sanity table + figure.

    Focus on three questions:
      Q1  HPMergedDelta vs AlleleDelta / vcf_AF  (L2 collider risk)
      Q2  HP1FamilyN + HP2FamilyN vs NumReads   (is it NumReads proxy?)
      Q3  HP_Ratio vs LOH_Subtype                (extreme in LOH?)
    """
    rows = []
    pairs_scatter = [
        ("HPMergedDelta", "AlleleDelta"),
        ("HPMergedDelta", "vcf_AF"),
        ("HP1FamilyN", "NumReads"),
        ("HP_FamilyN_sum", "NumReads"),
        ("HP_Ratio", "Potential_LOH"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    axes = axes.flatten()

    for idx, (a, b) in enumerate(pairs_scatter):
        ax = axes[idx]
        sub = df[[a, b]].dropna()
        if len(sub) < 50 or sub[a].std() == 0 or sub[b].std() == 0:
            ax.axis("off"); continue
        try:
            rho_p = stats.pearsonr(sub[a], sub[b])[0]
        except Exception:
            rho_p = np.nan
        try:
            rho_s = stats.spearmanr(sub[a], sub[b]).correlation
        except Exception:
            rho_s = np.nan
        if a == "HP_Ratio" and b == "Potential_LOH":
            # boolean b -> boxplot
            data = [sub.loc[sub[b] == 0, a].to_numpy(),
                    sub.loc[sub[b] == 1, a].to_numpy()]
            ax.boxplot(data, labels=["LOH=F", "LOH=T"],
                       showfliers=False)
            ax.axhline(0.5, ls=":", c="gray", lw=0.7)
            ax.set_title(f"{a} by {b}\n"
                         f"median(LOH=F)={np.nanmedian(data[0]):.3f}  "
                         f"median(LOH=T)={np.nanmedian(data[1]):.3f}")
            ax.set_ylabel(a)
        else:
            sl = sub.sample(n=min(50_000, len(sub)), random_state=7)
            ax.scatter(sl[b], sl[a], s=1, alpha=0.05)
            # Fit line
            try:
                k, c = np.polyfit(sl[b], sl[a], 1)
                xr = np.linspace(sl[b].min(), sl[b].max(), 20)
                ax.plot(xr, k*xr + c, c="red", lw=1,
                        label=f"y={k:.2f}x+{c:.2f}")
            except Exception:
                pass
            ax.set_xlabel(b); ax.set_ylabel(a)
            ax.set_title(f"{a} vs {b}\n"
                         f"Pearson={rho_p:.3f}  Spearman={rho_s:.3f}  n={len(sub):,}",
                         fontsize=9)
        rows.append({"a": a, "b": b, "pearson": rho_p, "spearman": rho_s,
                     "n": int(len(sub))})
    # last axes: note box
    axes[-1].axis("off")
    txt = ("G5 relation sanity notes\n\n"
           "Q1: HPMergedDelta vs AlleleDelta — if |r|>0.7,\n"
           "    HPMergedDelta is essentially AF echo → L2 collider\n"
           "    with tp_label (already flagged for AlleleDelta).\n\n"
           "Q2: (HP1FamilyN + HP2FamilyN)/NumReads ≈ 1 if\n"
           "    HP tags cover all reads.  frac_median reported below.\n\n"
           f"  HP_FamilyN_frac median: {df['HP_FamilyN_frac'].median():.3f}\n"
           f"  fraction frac>0.95:    {(df['HP_FamilyN_frac']>0.95).mean():.3f}\n"
           f"  fraction frac<0.20:    {(df['HP_FamilyN_frac']<0.20).mean():.3f}\n"
           "\nQ3: HP_Ratio — median shown in boxplot.  LOH ~0.09\n"
           "    means most reads → one haplotype, as expected.\n")
    axes[-1].text(0.01, 0.99, txt, transform=axes[-1].transAxes,
                  va="top", ha="left", family="monospace", fontsize=8)
    fig.suptitle("G5 · Step 0 Relation sanity (7 samples, both modes pooled)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / "fig00_relation_sanity.png", dpi=140)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(DATA_DIR / "G5_relations.tsv", sep="\t", index=False)
    return out


# --------------------------------------------------------------------------- #
# Step 1 — Global distribution & AUC panel                                    #
# --------------------------------------------------------------------------- #
def step1_global_panel(df: pd.DataFrame) -> pd.DataFrame:
    feats = FEATURES_CONT + FEATURES_BIN
    n_rows = 3; n_cols = 4
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 13))
    axes = axes.flatten()
    rows = []
    for i, feat in enumerate(feats):
        ax = axes[i]
        if feat not in df.columns:
            ax.axis("off"); continue
        vals = pd.to_numeric(df[feat], errors="coerce").to_numpy(dtype=float)
        y = df["tp_label"].to_numpy()
        tp = vals[y == 1]; fp = vals[y == 0]
        tp = tp[~np.isnan(tp)]; fp = fp[~np.isnan(fp)]
        if len(tp) < 10 or len(fp) < 10:
            ax.axis("off"); continue
        log_scale = feat in {"HP1FamilyN", "HP2FamilyN", "HP_FamilyN_sum"}
        p_tp = np.log10(tp + 1) if log_scale else tp
        p_fp = np.log10(fp + 1) if log_scale else fp
        if feat == "HPMergedSig":
            rate_tp = tp.mean()
            rate_fp = fp.mean()
            ax.bar([0, 1], [rate_fp, rate_tp],
                   color=["#E07A5F", "#2E86AB"], edgecolor="black")
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"])
            ax.set_ylabel("HPMergedSig rate")
            ax.set_ylim(0, 1.0)
            auc, lo, hi, _, _ = auc_wilson_ci(y, vals)
            ax.set_title(f"{feat}\nrate_TP={rate_tp:.3f} rate_FP={rate_fp:.3f}\n"
                         f"AUC={auc:.3f}")
            rows.append({"feature": feat, "auc": auc, "lo": lo, "hi": hi,
                         "cohen_d": np.nan, "mwu_p": np.nan,
                         "mean_tp": rate_tp, "mean_fp": rate_fp,
                         "n_tp": int(len(tp)), "n_fp": int(len(fp))})
            continue
        parts = ax.violinplot([p_fp, p_tp], positions=[0, 1], widths=0.7,
                              showmeans=False, showmedians=True)
        for pc, color in zip(parts["bodies"], ["#E07A5F", "#2E86AB"]):
            pc.set_facecolor(color); pc.set_alpha(0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"],
                           fontsize=8)
        auc, lo, hi, _, _ = auc_wilson_ci(y, vals)
        d = cohen_d(tp, fp); p = mwu_p(tp, fp)
        ax.set_title(f"{feat}{' (log10)' if log_scale else ''}\n"
                     f"AUC={auc:.3f} [{lo:.3f},{hi:.3f}] d={d:.2f}\n"
                     f"mw-p={p:.1e}", fontsize=9)
        rows.append({"feature": feat, "auc": auc, "lo": lo, "hi": hi,
                     "cohen_d": d, "mwu_p": p,
                     "mean_tp": float(np.nanmean(tp)),
                     "mean_fp": float(np.nanmean(fp)),
                     "n_tp": int(len(tp)), "n_fp": int(len(fp))})
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle("G5 HP Merged · Step 1 Global TP vs FP distribution "
                 "(7 samples × 2 modes pooled)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / "fig01_global_distribution.png", dpi=140)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(DATA_DIR / "G5_global_stats.tsv", sep="\t", index=False)
    return out


# --------------------------------------------------------------------------- #
# Step 2 — LOH × AF × CN 32-cell heatmap (for flagship continuous feats)     #
# --------------------------------------------------------------------------- #
def step2_heatmap(df: pd.DataFrame, feat: str, min_n: int = 20):
    if feat not in df.columns:
        return None
    df = df.copy()
    df["row_key"] = df["LOH_Subtype"].astype(str) + " | " + df["AF_class"].astype(str)
    row_order = [f"{l} | {a}" for l in LOH_ORDER for a in AF_ORDER]
    col_order = CN_LABELS

    piv_n = df.pivot_table(index="row_key", columns="cn_tier_F",
                           values="tp_label", aggfunc="size") \
              .reindex(row_order, axis=0).reindex(col_order, axis=1)
    piv_rate = df.pivot_table(index="row_key", columns="cn_tier_F",
                              values="tp_label", aggfunc="mean") \
                 .reindex(row_order, axis=0).reindex(col_order, axis=1)
    mask_ok = piv_n >= min_n
    piv_tp = df[df["tp_label"] == 1].pivot_table(
        index="row_key", columns="cn_tier_F", values=feat, aggfunc="median"
    ).reindex(row_order, axis=0).reindex(col_order, axis=1)
    piv_fp = df[df["tp_label"] == 0].pivot_table(
        index="row_key", columns="cn_tier_F", values=feat, aggfunc="median"
    ).reindex(row_order, axis=0).reindex(col_order, axis=1)
    piv_delta = (piv_tp - piv_fp).where(mask_ok)

    fig, axes = plt.subplots(1, 4, figsize=(24, 11))
    im0 = axes[0].imshow(piv_rate.where(mask_ok).values, aspect="auto",
                         cmap="RdBu_r", vmin=0, vmax=1)
    axes[0].set_title(f"A · TP rate (n≥{min_n})")
    fig.colorbar(im0, ax=axes[0], fraction=0.03)
    im1 = axes[1].imshow(piv_tp.where(mask_ok).values, aspect="auto",
                         cmap="Blues")
    axes[1].set_title(f"B · median({feat}) | TP only")
    fig.colorbar(im1, ax=axes[1], fraction=0.03)
    im2 = axes[2].imshow(piv_fp.where(mask_ok).values, aspect="auto",
                         cmap="Oranges")
    axes[2].set_title(f"C · median({feat}) | FP only")
    fig.colorbar(im2, ax=axes[2], fraction=0.03)
    dlt = piv_delta.values.astype(float)
    finite = dlt[np.isfinite(dlt)]
    if len(finite):
        vlim = max(abs(np.nanpercentile(finite, 5)),
                   abs(np.nanpercentile(finite, 95)), 1e-6)
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
    for i, r in enumerate(row_order):
        for j, c in enumerate(col_order):
            n = piv_n.loc[r, c] if (r in piv_n.index and c in piv_n.columns) else 0
            if pd.isna(n): n = 0
            if n >= min_n:
                axes[0].text(j, i, f"{int(n):,}", ha="center", va="center",
                             fontsize=6)
    fig.suptitle(f"{feat} · G5 Step 2 LOH × AF × CN heatmap (pooled)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / f"fig02_{feat}_heatmap.png", dpi=130)
    plt.close(fig)

    flat = piv_delta.reset_index().melt(id_vars="row_key", var_name="cn_tier",
                                        value_name="delta_tp_fp")
    flat["feature"] = feat
    return flat


# --------------------------------------------------------------------------- #
# Step 3 — Per-sample consistency                                             #
# --------------------------------------------------------------------------- #
def step3_per_sample(df: pd.DataFrame, feat: str, min_n: int = 10):
    if feat not in df.columns:
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
        im = axes[i].imshow(piv.where(mask).values, aspect="auto",
                            cmap="RdBu_r", vmin=0, vmax=1)
        axes[i].set_title(f"{s}\n(n={len(sub):,})", fontsize=9)
        axes[i].set_xticks(range(len(AF_ORDER)))
        axes[i].set_xticklabels(AF_ORDER, rotation=30, ha="right", fontsize=7)
        axes[i].set_yticks(range(len(LOH_ORDER)))
        axes[i].set_yticklabels(LOH_ORDER, fontsize=7)
        for ii, r in enumerate(LOH_ORDER):
            for jj, c in enumerate(AF_ORDER):
                if r in piv_n.index and c in piv_n.columns:
                    n = piv_n.loc[r, c]
                    if pd.isna(n): n = 0
                    if n >= min_n:
                        axes[i].text(jj, ii, f"{piv.loc[r, c]:.2f}\nn={int(n)}",
                                     ha="center", va="center", fontsize=6)
        per_sample_rate[s] = piv.values.flatten()
    ax_conc = axes[-1]
    mat = np.full((n_samp, n_samp), np.nan)
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            v1 = per_sample_rate[s1]; v2 = per_sample_rate[s2]
            m = ~(np.isnan(v1) | np.isnan(v2))
            if m.sum() >= 3:
                mat[i, j] = stats.spearmanr(v1[m], v2[m]).correlation
    im2 = ax_conc.imshow(mat, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax_conc.set_title("Cross-sample Spearman concordance", fontsize=9)
    ax_conc.set_xticks(range(n_samp))
    ax_conc.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax_conc.set_yticks(range(n_samp))
    ax_conc.set_yticklabels(samples, fontsize=7)
    fig.colorbar(im2, ax=ax_conc, fraction=0.04)
    fig.suptitle(f"{feat} · G5 Step 3 Per-sample TP rate consistency (LOH × AF)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / f"fig03_{feat}_per_sample.png", dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Step 4 — Stratified AUC                                                     #
# --------------------------------------------------------------------------- #
def step4_stratified(df: pd.DataFrame, feat: str):
    if feat not in df.columns:
        return None
    y = df["tp_label"].to_numpy()
    s = pd.to_numeric(df[feat], errors="coerce").to_numpy(dtype=float)
    rows = [("global", "all", *auc_wilson_ci(y, s))]
    for loh in LOH_ORDER:
        m = (df["LOH_Subtype"] == loh).to_numpy()
        if m.sum() >= 100:
            rows.append(("LOH", loh, *auc_wilson_ci(y[m], s[m])))
    for af in AF_ORDER:
        m = (df["AF_class"] == af).to_numpy()
        if m.sum() >= 100:
            rows.append(("AF", af, *auc_wilson_ci(y[m], s[m])))
    for cn in CN_LABELS:
        m = (df["cn_tier_F"] == cn).to_numpy()
        if m.sum() >= 100:
            rows.append(("CN", cn, *auc_wilson_ci(y[m], s[m])))
    for mode in ["paired_full", "to_pileup"]:
        m = (df["mode"] == mode).to_numpy()
        if m.sum() >= 100:
            rows.append(("mode", mode, *auc_wilson_ci(y[m], s[m])))
    dfres = pd.DataFrame(rows, columns=["layer", "group", "auc", "lo", "hi",
                                        "n_pos", "n_neg"])
    dfres["feature"] = feat

    fig, ax = plt.subplots(figsize=(12, 5.5))
    colors = {"global": "#333", "LOH": "#2E86AB", "AF": "#F4A261",
              "CN": "#6A994E", "mode": "#8E44AD"}
    for i, (_, row) in enumerate(dfres.iterrows()):
        yerr = [[max(0, row["auc"]-row["lo"])], [max(0, row["hi"]-row["auc"])]]
        ax.bar(i, row["auc"], color=colors.get(row["layer"], "gray"),
               yerr=yerr, capsize=3, edgecolor="black", linewidth=0.4)
        ax.text(i, row["auc"] + 0.01,
                f"{row['auc']:.3f}\nn={int(row['n_pos'])+int(row['n_neg']):,}",
                ha="center", fontsize=6)
    ax.axhline(0.5, ls="--", c="gray", lw=0.8, label="random")
    ax.axhline(0.58, ls=":", c="red", lw=0.8, label="Beyond-AUC 0.58")
    ax.set_xticks(range(len(dfres)))
    ax.set_xticklabels([f"{r.layer}:{r.group}" for _, r in dfres.iterrows()],
                       rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0.3, 1.0)
    ax.set_ylabel("AUC")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"{feat} · G5 Step 4 Stratified AUC")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"fig04_{feat}_stratified_auc.png", dpi=130)
    plt.close(fig)
    return dfres


# --------------------------------------------------------------------------- #
# Step 5 — Confound guard (OLS residualization)                              #
# --------------------------------------------------------------------------- #
def step5_confound(df: pd.DataFrame, feat: str):
    if feat not in df.columns:
        return None
    # Residualize on (AlleleDelta, NumReads, vcf_AF, Coverage_Multiple)
    covars = [c for c in ["AlleleDelta", "NumReads", "vcf_AF", "Coverage_Multiple"]
              if c != feat and c in df.columns]
    sub = df.dropna(subset=[feat] + covars + ["tp_label"]).copy()
    if len(sub) < 200:
        return None
    X = sub[covars].astype(float).to_numpy()
    y_feat = pd.to_numeric(sub[feat], errors="coerce").to_numpy(dtype=float)
    y_tp = sub["tp_label"].to_numpy()
    try:
        lr = LinearRegression().fit(X, y_feat)
        resid = y_feat - lr.predict(X)
    except Exception:
        return None
    raw_auc = auc_wilson_ci(y_tp, y_feat)[0]
    res_auc = auc_wilson_ci(y_tp, resid)[0]
    af_rows = []
    for af in AF_ORDER:
        m = (sub["AF_class"] == af).to_numpy()
        if m.sum() >= 200:
            r_auc = auc_wilson_ci(y_tp[m], y_feat[m])[0]
            af_rows.append((af, r_auc, int(m.sum())))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    bins = 50
    ax.hist(y_feat[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
            density=True, label="raw FP")
    ax.hist(y_feat[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
            density=True, label="raw TP")
    ax.set_title(f"raw {feat}  AUC={raw_auc:.3f}")
    ax.legend(fontsize=8)
    ax.set_xlabel(feat)
    ax2 = axes[1]
    ax2.hist(resid[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
             density=True, label="resid FP")
    ax2.hist(resid[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
             density=True, label="resid TP")
    ax2.set_title(f"resid. on {covars}\nAUC={res_auc:.3f}\n"
                  f"AF-bin raw AUCs: " +
                  ", ".join([f"{a}={v:.3f}" for a, v, _ in af_rows]),
                  fontsize=9)
    ax2.legend(fontsize=8)
    fig.suptitle(f"{feat} · G5 Step 5 Confound guard "
                 f"(within-cell OLS residualization)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / f"fig05_{feat}_confound.png", dpi=130)
    plt.close(fig)
    return {"feature": feat, "raw_auc": raw_auc, "resid_auc": res_auc,
            "covars": ",".join(covars),
            "af_bin_aucs": {a: v for a, v, _ in af_rows}}


# --------------------------------------------------------------------------- #
# Step 6 — Spatial auto-correlation                                           #
# --------------------------------------------------------------------------- #
def step6_spatial(df: pd.DataFrame, feat: str, bin_mb: int = 5):
    if feat not in df.columns:
        return None
    d = df.dropna(subset=["Chr", "Pos", feat, "tp_label"]).copy()
    d["pos_bin"] = (d["Pos"].astype(int) // (bin_mb * 1_000_000))
    d["bin_key"] = d["Chr"].astype(str) + ":" + d["pos_bin"].astype(str)
    rows = []
    for key, g in d.groupby("bin_key"):
        if len(g) < 50:
            continue
        y = g["tp_label"].to_numpy()
        if y.sum() < 5 or (len(y) - y.sum()) < 5:
            continue
        a, *_ = auc_wilson_ci(y, pd.to_numeric(g[feat], errors="coerce")
                              .to_numpy(dtype=float))
        rows.append({"bin": key, "n": len(g), "tp_rate": y.mean(), "auc": a})
    if not rows:
        return None
    s = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    ax.scatter(s["tp_rate"], s["auc"], alpha=0.3, s=8)
    ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
    ax.axvline(0.8, ls="--", c="purple", alpha=0.5)
    ax.set_xlabel("bin TP rate"); ax.set_ylabel(f"bin AUC for {feat}")
    ax2 = axes[1]
    ax2.hist(s["auc"], bins=30, color="#2E86AB", alpha=0.7)
    ax2.axvline(0.5, ls="--", c="gray"); ax2.axvline(0.58, ls=":", c="red")
    ax2.set_xlabel("per-bin AUC"); ax2.set_ylabel("# bins")
    high_tp = s[s["tp_rate"] > 0.8]; low_tp = s[s["tp_rate"] < 0.5]
    flag = ""
    if len(high_tp) >= 10 and len(low_tp) >= 10:
        hi_mean = high_tp["auc"].mean(); lo_mean = low_tp["auc"].mean()
        if hi_mean - lo_mean > 0.08:
            flag = f"⚠ artifact suspect (ΔAUC={hi_mean-lo_mean:+.3f})"
        else:
            flag = f"ok (ΔAUC={hi_mean-lo_mean:+.3f})"
    fig.suptitle(f"{feat} · G5 Step 6 Spatial auto-corr ({bin_mb} Mb) {flag}",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / f"fig06_{feat}_spatial.png", dpi=130)
    plt.close(fig)
    return {"feature": feat, "n_bins": len(s), "flag": flag}


# --------------------------------------------------------------------------- #
# Step 7 — HP_Ratio normalization comparison                                  #
# --------------------------------------------------------------------------- #
def step7_hp_ratio_norm(df: pd.DataFrame):
    """Compare HP_Ratio vs HP_Ratio_norm = |HP_Ratio - 0.5| per LOH stratum."""
    sub = df.dropna(subset=["HP_Ratio", "tp_label"]).copy()
    strata = [("ALL", np.ones(len(sub), dtype=bool))]
    strata += [(l, (sub["LOH_Subtype"] == l).to_numpy()) for l in LOH_ORDER]
    strata += [("Potential_LOH=1", (sub["Potential_LOH"] == 1).to_numpy()),
               ("Potential_LOH=0", (sub["Potential_LOH"] == 0).to_numpy())]
    rows = []
    for name, m in strata:
        if m.sum() < 100:
            continue
        y = sub.loc[m, "tp_label"].to_numpy()
        r = sub.loc[m, "HP_Ratio"].to_numpy()
        rn = (sub.loc[m, "HP_Ratio"] - 0.5).abs().to_numpy()
        auc_r, *_ = auc_wilson_ci(y, r)
        auc_rn, *_ = auc_wilson_ci(y, rn)
        rows.append({"stratum": name, "n": int(m.sum()),
                     "tp_rate": float(y.mean()),
                     "auc_HP_Ratio": auc_r, "auc_HP_Ratio_norm": auc_rn,
                     "median_HP_Ratio_TP": float(np.nanmedian(r[y == 1])) if y.sum() > 0 else np.nan,
                     "median_HP_Ratio_FP": float(np.nanmedian(r[y == 0])) if (len(y) - y.sum()) > 0 else np.nan,
                     "median_norm_TP": float(np.nanmedian(rn[y == 1])) if y.sum() > 0 else np.nan,
                     "median_norm_FP": float(np.nanmedian(rn[y == 0])) if (len(y) - y.sum()) > 0 else np.nan})
    out = pd.DataFrame(rows)
    out.to_csv(DATA_DIR / "G5_hp_ratio_normalization.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    ax = axes[0]
    x = np.arange(len(out))
    ax.bar(x - 0.18, out["auc_HP_Ratio"], width=0.35,
           color="#F4A261", label="HP_Ratio (raw)", edgecolor="black")
    ax.bar(x + 0.18, out["auc_HP_Ratio_norm"], width=0.35,
           color="#2E86AB", label="|HP_Ratio − 0.5|", edgecolor="black")
    ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
    ax.set_xticks(x)
    ax.set_xticklabels(out["stratum"], rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0.3, 1.0)
    ax.set_ylabel("AUC (TP discrimination)")
    ax.set_title("HP_Ratio raw vs |HP_Ratio−0.5| across strata")
    ax.legend(fontsize=8)
    for xi, (_, row) in zip(x, out.iterrows()):
        ax.text(xi - 0.18, row["auc_HP_Ratio"] + 0.005, f"{row['auc_HP_Ratio']:.2f}",
                ha="center", fontsize=6)
        ax.text(xi + 0.18, row["auc_HP_Ratio_norm"] + 0.005, f"{row['auc_HP_Ratio_norm']:.2f}",
                ha="center", fontsize=6)
    # Panel B: distribution of HP_Ratio per LOH subtype (TP vs FP)
    ax = axes[1]
    loh_subs = LOH_ORDER
    positions = []; labels = []
    offset = 0
    for loh in loh_subs:
        for tp, color in [(0, "#E07A5F"), (1, "#2E86AB")]:
            v = sub.loc[(sub["LOH_Subtype"] == loh) & (sub["tp_label"] == tp),
                        "HP_Ratio"].dropna().to_numpy()
            if len(v) < 20:
                offset += 1
                continue
            pos = offset
            ax.boxplot([v], positions=[pos], widths=0.7, showfliers=False,
                       boxprops=dict(color=color, lw=1.5),
                       medianprops=dict(color="black"))
            offset += 1
        offset += 0.5
    ax.axhline(0.5, ls="--", c="gray", lw=0.6)
    ax.set_ylabel("HP_Ratio")
    ax.set_title("HP_Ratio distribution by LOH × TP/FP\n"
                 "blue=TP orange=FP  groups: LOH_None, Weak, Noise, Strong, Subclone")
    ax.set_xticks([])
    fig.suptitle("G5 Step 7 HP_Ratio normalization benchmark", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG_DIR / "fig07_hp_ratio_normalization.png", dpi=130)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Step 8 — HPMergedSig per-class TP rate                                      #
# --------------------------------------------------------------------------- #
def step8_hpmergedsig_per_class(df: pd.DataFrame):
    sub = df.dropna(subset=["HPMergedSig", "tp_label"]).copy()
    rows = []
    # global
    for cls in [0, 1]:
        m = (sub["HPMergedSig"] == cls).to_numpy()
        y = sub.loc[m, "tp_label"].to_numpy()
        k = int(y.sum()); n = int(len(y))
        p, lo, hi = wilson_ci_k_n(k, n)
        rows.append({"stratum": "ALL", "HPMergedSig": cls,
                     "tp_rate": p, "lo": lo, "hi": hi, "n": n})
    for loh in LOH_ORDER:
        for cls in [0, 1]:
            m = ((sub["LOH_Subtype"] == loh) &
                 (sub["HPMergedSig"] == cls)).to_numpy()
            y = sub.loc[m, "tp_label"].to_numpy()
            if len(y) < 20:
                continue
            k = int(y.sum()); n = int(len(y))
            p, lo, hi = wilson_ci_k_n(k, n)
            rows.append({"stratum": loh, "HPMergedSig": cls,
                         "tp_rate": p, "lo": lo, "hi": hi, "n": n})
    for mode in ["paired_full", "to_pileup"]:
        for cls in [0, 1]:
            m = ((sub["mode"] == mode) &
                 (sub["HPMergedSig"] == cls)).to_numpy()
            y = sub.loc[m, "tp_label"].to_numpy()
            if len(y) < 20:
                continue
            k = int(y.sum()); n = int(len(y))
            p, lo, hi = wilson_ci_k_n(k, n)
            rows.append({"stratum": mode, "HPMergedSig": cls,
                         "tp_rate": p, "lo": lo, "hi": hi, "n": n})
    out = pd.DataFrame(rows)
    out.to_csv(DATA_DIR / "G5_HPMergedSig_per_class.tsv", sep="\t", index=False)

    # figure
    fig, ax = plt.subplots(figsize=(14, 6))
    strata = out["stratum"].unique().tolist()
    x = np.arange(len(strata))
    widths = 0.35
    for cls, color, offset in [(0, "#E07A5F", -0.18), (1, "#2E86AB", 0.18)]:
        sub_c = out[out["HPMergedSig"] == cls].set_index("stratum") \
                     .reindex(strata)
        yerr = [sub_c["tp_rate"] - sub_c["lo"], sub_c["hi"] - sub_c["tp_rate"]]
        ax.bar(x + offset, sub_c["tp_rate"], width=widths,
               color=color, edgecolor="black",
               yerr=yerr, capsize=3,
               label=f"HPMergedSig={cls}")
        for xi, v, n in zip(x + offset, sub_c["tp_rate"], sub_c["n"]):
            if pd.notna(v):
                ax.text(xi, v + 0.01, f"{v:.2f}\nn={int(n)}",
                        ha="center", fontsize=6)
    ax.axhline(0.5, ls="--", c="gray", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(strata, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("TP rate")
    ax.legend(fontsize=9)
    ax.set_title("G5 Step 8 HPMergedSig TP rate by stratum (Wilson 95% CI)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig08_hpmergedsig_per_class.png", dpi=130)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #
def main():
    df = load()
    print("[main] Step 0 relations", flush=True)
    step0_relations(df)

    print("[main] Step 1 global", flush=True)
    step1_global_panel(df)

    focus_feats = ["HP_Ratio", "HP_Ratio_norm", "HPMergedDelta", "HPMergedSig"]
    cell_deltas = []
    auc_rows = []
    confound_rows = []
    spatial_rows = []
    for feat in focus_feats:
        print(f"[main] Step 2 heatmap {feat}", flush=True)
        f = step2_heatmap(df, feat)
        if f is not None:
            cell_deltas.append(f)
        print(f"[main] Step 3 per-sample {feat}", flush=True)
        step3_per_sample(df, feat)
        print(f"[main] Step 4 stratified {feat}", flush=True)
        r = step4_stratified(df, feat)
        if r is not None:
            auc_rows.append(r)
        if feat not in {"HPMergedSig"}:
            print(f"[main] Step 5 confound {feat}", flush=True)
            c = step5_confound(df, feat)
            if c is not None:
                confound_rows.append(c)
            print(f"[main] Step 6 spatial {feat}", flush=True)
            sp = step6_spatial(df, feat)
            if sp is not None:
                spatial_rows.append(sp)

    print("[main] Step 7 HP_Ratio normalization", flush=True)
    step7_hp_ratio_norm(df)

    print("[main] Step 8 HPMergedSig per-class", flush=True)
    step8_hpmergedsig_per_class(df)

    if auc_rows:
        pd.concat(auc_rows, ignore_index=True).to_csv(
            DATA_DIR / "G5_auc_table.tsv", sep="\t", index=False)
    if cell_deltas:
        pd.concat(cell_deltas, ignore_index=True).to_csv(
            DATA_DIR / "G5_cell_delta.tsv", sep="\t", index=False)
    if confound_rows:
        pd.DataFrame(confound_rows).to_csv(
            DATA_DIR / "G5_confound.tsv", sep="\t", index=False)
    if spatial_rows:
        pd.DataFrame(spatial_rows).to_csv(
            DATA_DIR / "G5_spatial.tsv", sep="\t", index=False)
    print("[main] done. Outputs:", flush=True)
    print(f"  figures: {FIG_DIR}")
    print(f"  tables : {DATA_DIR}")


if __name__ == "__main__":
    main()
