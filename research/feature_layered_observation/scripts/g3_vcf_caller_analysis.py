"""Phase C G3 — Caller VCF Signals feature group analysis.

Executes SOP Step 1-6 on caller-derived VCF features:
    vcf_QUAL, vcf_AF, vcf_NAF, vcf_DP, vcf_NDP, vcf_AD_{ref,alt}, vcf_NAD_{ref,alt},
    vcf_GQ, vcf_SB, vcf_PoN_1..4, vcf_Verdict_{Germline,Somatic,SubclonalSomatic},
    vcf_H, vcf_{F,R}{A,C,G,T}U (strand-allele counts)

Notes
-----
- master dataset is PRE-filtered to PASS only (vcf_FILTER == 'PASS' for 100% rows)
    -> Step 4 per-FILTER stratum reports "PASS-only (non-PASS rows absent by construction)"
- vcf_NAF / vcf_NDP / vcf_NAD_* / vcf_GQ exist only in paired_full
- vcf_SB / vcf_PoN_* / vcf_Verdict_* exist only in to_pileup
- vcf_Verdict_* are all 0 in PASS subset (Verdict_Germline lives in LowQual-only rows)
    -> re-validates the 0420 Pilot finding on all 7 samples
- Step 2 stratifier uses vcf_AF (caller VAF) NOT master AF (|AlleleDelta|)
"""
from __future__ import annotations

import gzip
import json
import warnings
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, norm
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "merged_with_vcf.tsv.gz"
FIG_DIR = ROOT / "figures" / "G3_vcf_caller"
FEAT_DIR = ROOT / "features"
FIG_DIR.mkdir(parents=True, exist_ok=True)
FEAT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
    "H2009",
    "H1437",
    "COLO829",
]

# ---------------------------------------------------------------------------
# Feature groups
# ---------------------------------------------------------------------------
CONTINUOUS_BOTH = [
    "vcf_QUAL",
    "vcf_AF",
    "vcf_DP",
    "vcf_AD_ref",
    "vcf_AD_alt",
    "vcf_H",  # 0/1 flag but also a quick binary AUC target
]
CONTINUOUS_PAIRED_ONLY = [
    "vcf_NAF",
    "vcf_NDP",
    "vcf_NAD_ref",
    "vcf_NAD_alt",
    "vcf_GQ",
]
CONTINUOUS_TO_ONLY = [
    "vcf_SB",
]
BINARY_TO_ONLY = [
    "vcf_PoN_1",
    "vcf_PoN_2",
    "vcf_PoN_3",
    "vcf_PoN_4",
    "vcf_Verdict_Germline",
    "vcf_Verdict_Somatic",
    "vcf_Verdict_SubclonalSomatic",
]
STRAND_ALLELE_COUNTS = [
    "vcf_FAU", "vcf_FCU", "vcf_FGU", "vcf_FTU",
    "vcf_RAU", "vcf_RCU", "vcf_RGU", "vcf_RTU",
]
ALL_FEATURES = (
    CONTINUOUS_BOTH
    + CONTINUOUS_PAIRED_ONLY
    + CONTINUOUS_TO_ONLY
    + BINARY_TO_ONLY
    + STRAND_ALLELE_COUNTS
)

STRATIFIERS = ["LOH_Subtype", "AF_class", "cn_tier_F"]
CONFOUNDERS = ["NumReads", "vcf_AF", "Coverage_Multiple"]  # Step 5

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def wilson_auc_ci(auc: float, n_pos: int, n_neg: int, alpha: float = 0.05) -> tuple[float, float]:
    """Hanley-McNeil asymptotic SE for AUC -> normal CI (no continuity correction)."""
    if n_pos == 0 or n_neg == 0 or np.isnan(auc):
        return (float("nan"), float("nan"))
    q1 = auc / (2 - auc)
    q2 = 2 * auc * auc / (1 + auc)
    var = (
        auc * (1 - auc)
        + (n_pos - 1) * (q1 - auc * auc)
        + (n_neg - 1) * (q2 - auc * auc)
    ) / (n_pos * n_neg)
    se = float(np.sqrt(max(var, 0)))
    z = norm.ppf(1 - alpha / 2)
    return (max(auc - z * se, 0.0), min(auc + z * se, 1.0))


def compute_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    mask = ~(np.isnan(y_score) | np.isnan(y_true))
    if mask.sum() < 10:
        return float("nan")
    y_t = y_true[mask].astype(int)
    y_s = y_score[mask]
    if y_t.sum() == 0 or y_t.sum() == len(y_t):
        return float("nan")
    try:
        return float(roc_auc_score(y_t, y_s))
    except ValueError:
        return float("nan")


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    if pooled == 0:
        return float("nan")
    return float((a.mean() - b.mean()) / pooled)


def mw_u_pval(tp: np.ndarray, fp: np.ndarray) -> float:
    tp = tp[~np.isnan(tp)]
    fp = fp[~np.isnan(fp)]
    if len(tp) < 5 or len(fp) < 5:
        return float("nan")
    try:
        _, p = mannwhitneyu(tp, fp, alternative="two-sided")
        return float(p)
    except ValueError:
        return float("nan")


def make_cn_tier_F(cov_mult: pd.Series) -> pd.Series:
    """5-tier CN bucket from Coverage_Multiple (per methodology)."""
    bins = [-np.inf, 0.6, 0.85, 1.2, 1.8, np.inf]
    labels = ["Loss", "Low", "Normal", "Elevated", "Gain"]
    return pd.cut(cov_mult, bins=bins, labels=labels).astype(str)


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """OLS residualize y on X; returns residuals (y - Xbeta)."""
    mask = ~(np.isnan(y) | np.any(np.isnan(X), axis=1))
    if mask.sum() < 20:
        return np.full_like(y, np.nan, dtype=float)
    reg = LinearRegression()
    reg.fit(X[mask], y[mask])
    pred = np.full_like(y, np.nan, dtype=float)
    pred[mask] = reg.predict(X[mask])
    return y - pred


# ---------------------------------------------------------------------------
# Core analyses
# ---------------------------------------------------------------------------

def step1_global_distribution(df: pd.DataFrame, feature: str, mode: str) -> dict:
    sub = df[df["mode"] == mode].copy()
    if feature not in sub.columns:
        return {}
    vals = pd.to_numeric(sub[feature], errors="coerce")
    y = sub["tp_label"].values.astype(int)
    tp_vals = vals[y == 1].values
    fp_vals = vals[y == 0].values
    if np.isnan(tp_vals).all() or np.isnan(fp_vals).all():
        return {
            "feature": feature,
            "mode": mode,
            "available": False,
            "note": "all-NaN in one class",
        }
    auc = compute_auc(y, vals.values)
    auc_inv = 1 - auc if not np.isnan(auc) else float("nan")
    # record both; report the one farther from 0.5 as "signed"
    final_auc = max(auc, auc_inv) if not np.isnan(auc) else float("nan")
    direction = "TP>FP" if (not np.isnan(auc) and auc >= 0.5) else "FP>TP"
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    ci_lo, ci_hi = wilson_auc_ci(final_auc, n_pos, n_neg)
    return {
        "feature": feature,
        "mode": mode,
        "available": True,
        "n_tp": int((~np.isnan(tp_vals)).sum()),
        "n_fp": int((~np.isnan(fp_vals)).sum()),
        "mean_tp": float(np.nanmean(tp_vals)),
        "mean_fp": float(np.nanmean(fp_vals)),
        "median_tp": float(np.nanmedian(tp_vals)),
        "median_fp": float(np.nanmedian(fp_vals)),
        "cohens_d": cohens_d(tp_vals, fp_vals),
        "mw_u_p": mw_u_pval(tp_vals, fp_vals),
        "auc_raw": float(auc) if not np.isnan(auc) else float("nan"),
        "auc_abs": float(final_auc) if not np.isnan(final_auc) else float("nan"),
        "direction": direction,
        "auc_ci_lo": ci_lo,
        "auc_ci_hi": ci_hi,
    }


def plot_global_distribution(df: pd.DataFrame, features: list[str], mode: str, png_path: Path, stats_rows: list[dict]):
    sub = df[df["mode"] == mode].copy()
    n_feat = len(features)
    if n_feat == 0:
        return
    ncols = 3
    nrows = int(np.ceil(n_feat / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.5, nrows * 3.5))
    axes = np.array(axes).reshape(-1)
    stats_map = {(r["feature"], r["mode"]): r for r in stats_rows if r.get("available")}
    for i, feat in enumerate(features):
        ax = axes[i]
        if feat not in sub.columns:
            ax.axis("off")
            continue
        vals = pd.to_numeric(sub[feat], errors="coerce")
        tp_mask = sub["tp_label"] == 1
        fp_mask = sub["tp_label"] == 0
        tp_vals = vals[tp_mask].dropna()
        fp_vals = vals[fp_mask].dropna()
        if tp_vals.empty or fp_vals.empty:
            ax.text(0.5, 0.5, f"{feat}: N/A in {mode}", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue
        # violin
        data = [fp_vals.values, tp_vals.values]
        parts = ax.violinplot(data, positions=[0, 1], showmeans=False, showmedians=True, widths=0.8)
        for pc, color in zip(parts["bodies"], ["#ff9f43", "#3498db"]):
            pc.set_facecolor(color)
            pc.set_alpha(0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"FP (n={len(fp_vals)})", f"TP (n={len(tp_vals)})"], fontsize=8)
        st = stats_map.get((feat, mode), {})
        auc = st.get("auc_abs", float("nan"))
        d = st.get("cohens_d", float("nan"))
        dirn = st.get("direction", "")
        ax.set_title(f"{feat} [{mode}]\nAUC={auc:.3f} d={d:+.2f} {dirn}", fontsize=9)
        ax.set_ylabel(feat, fontsize=8)
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"G3 · Step 1 Global TP vs FP — {mode}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(png_path, dpi=130)
    plt.close(fig)


def step2_strata_heatmaps(df: pd.DataFrame, feature: str, mode: str, png_dir: Path):
    sub = df[df["mode"] == mode].copy()
    if feature not in sub.columns:
        return
    # ensure stratifier columns
    for s in STRATIFIERS:
        if s not in sub.columns:
            return
    vals = pd.to_numeric(sub[feature], errors="coerce")
    if vals.notna().sum() < 100:
        return
    sub = sub.assign(_f=vals)
    # Cross stratify: (LOH_Subtype × AF_class × cn_tier_F) -> long form, then show as LOH×cn pivot for 3 AF_class panels
    af_levels = ["Extreme", "Near-half", "Intermediate"]
    loh_levels = ["None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"]
    cn_levels = ["Loss", "Low", "Normal", "Elevated", "Gain"]

    agg = (
        sub.groupby(["AF_class", "LOH_Subtype", "cn_tier_F"])
        .agg(n=("tp_label", "size"),
             tp=("tp_label", "sum"),
             tp_feat=("_f", lambda x: np.nanmean(x[sub.loc[x.index, "tp_label"] == 1])),
             fp_feat=("_f", lambda x: np.nanmean(x[sub.loc[x.index, "tp_label"] == 0])))
        .reset_index()
    )
    agg["tp_rate"] = agg["tp"] / agg["n"]
    agg["delta"] = agg["tp_feat"] - agg["fp_feat"]
    agg.loc[agg["n"] < 20, ["tp_rate", "tp_feat", "fp_feat", "delta"]] = np.nan

    # Plot 4 heatmaps (tp_rate / tp_feat / fp_feat / delta) with 3 AF_class columns × LOH rows × cn cols
    fig, axes = plt.subplots(4, 3, figsize=(14, 16))
    heat_types = [("tp_rate", "TP rate", "viridis", 0, 1),
                  ("tp_feat", f"{feature} mean (TP)", "Blues", None, None),
                  ("fp_feat", f"{feature} mean (FP)", "Oranges", None, None),
                  ("delta", f"Δ(TP-FP) {feature}", "RdBu_r", None, None)]
    for row_i, (col, title, cmap, vmin, vmax) in enumerate(heat_types):
        for col_j, af in enumerate(af_levels):
            ax = axes[row_i, col_j]
            piv = (
                agg[agg["AF_class"] == af]
                .pivot(index="LOH_Subtype", columns="cn_tier_F", values=col)
                .reindex(index=loh_levels, columns=cn_levels)
            )
            n_piv = (
                agg[agg["AF_class"] == af]
                .pivot(index="LOH_Subtype", columns="cn_tier_F", values="n")
                .reindex(index=loh_levels, columns=cn_levels)
            )
            if col == "delta":
                absmax = np.nanmax(np.abs(piv.values)) if np.isfinite(np.nanmax(np.abs(piv.values))) else 1
                sns.heatmap(piv, ax=ax, cmap=cmap, center=0, vmin=-absmax, vmax=absmax,
                            annot=True, fmt=".2f", cbar_kws={"shrink": 0.7}, linewidths=0.3)
            elif vmin is not None:
                sns.heatmap(piv, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
                            annot=True, fmt=".2f", cbar_kws={"shrink": 0.7}, linewidths=0.3)
            else:
                sns.heatmap(piv, ax=ax, cmap=cmap,
                            annot=True, fmt=".2f", cbar_kws={"shrink": 0.7}, linewidths=0.3)
            # annotate n in cell corner
            for yi in range(piv.shape[0]):
                for xi in range(piv.shape[1]):
                    n = n_piv.iat[yi, xi] if not np.isnan(n_piv.iat[yi, xi]) else 0
                    ax.text(xi + 0.5, yi + 0.85, f"n={int(n)}", ha="center", va="top",
                            fontsize=6, color="gray")
            ax.set_title(f"{title} | AF={af}", fontsize=9)
            ax.set_xlabel("CN tier (CovM)")
            ax.set_ylabel("LOH_Subtype" if col_j == 0 else "")
    fig.suptitle(f"G3 · Step 2 · LOH × AF × CN stratification — {feature} [{mode}]", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = png_dir / f"fig02_{feature}_{mode}_strata.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)


def step3_cross_sample(df: pd.DataFrame, feature: str, mode: str, png_dir: Path):
    sub = df[df["mode"] == mode].copy()
    if feature not in sub.columns:
        return
    if pd.to_numeric(sub[feature], errors="coerce").notna().sum() < 100:
        return
    # per sample TP rate on 3×5 LOH×AF
    samples = [s for s in SAMPLE_ORDER if s in sub["sample"].unique()]
    if not samples:
        return
    sub = sub.assign(_f=pd.to_numeric(sub[feature], errors="coerce"))

    # build per-sample × (LOH×AF×CN) TP rate signatures
    agg = (
        sub.groupby(["sample", "LOH_Subtype", "AF_class", "cn_tier_F"])
        .agg(n=("tp_label", "size"), tp=("tp_label", "sum"))
        .reset_index()
    )
    agg["tp_rate"] = agg["tp"] / agg["n"]
    agg.loc[agg["n"] < 10, "tp_rate"] = np.nan
    agg["cell"] = agg["LOH_Subtype"] + "|" + agg["AF_class"] + "|" + agg["cn_tier_F"]

    wide = agg.pivot(index="cell", columns="sample", values="tp_rate")
    # drop cells absent in >=4 samples
    wide = wide.dropna(thresh=max(1, int(len(samples) * 0.5)))
    if wide.shape[0] < 3:
        spear = None
    else:
        spear = wide.corr(method="spearman").reindex(index=samples, columns=samples)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    # feature mean (TP only) per sample × 5 LOH × 3 AF grid
    ax = axes[0]
    pvt = (
        sub[sub["tp_label"] == 1]
        .groupby(["sample", "LOH_Subtype", "AF_class"])["_f"].mean()
        .reset_index()
    )
    pvt["row"] = pvt["LOH_Subtype"] + "|" + pvt["AF_class"]
    ptable = pvt.pivot(index="row", columns="sample", values="_f").reindex(columns=samples)
    sns.heatmap(ptable, ax=ax, cmap="Blues", annot=True, fmt=".2f",
                cbar_kws={"shrink": 0.7}, linewidths=0.3)
    ax.set_title(f"Per-sample {feature} mean (TP only)\nrows = LOH|AF_class", fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")

    ax = axes[1]
    if spear is not None and spear.notna().sum().sum() > 0:
        sns.heatmap(spear, ax=ax, cmap="RdBu_r", vmin=-1, vmax=1, center=0,
                    annot=True, fmt=".2f", cbar_kws={"shrink": 0.7}, linewidths=0.3)
        ax.set_title(f"Spearman(cell TP rate) across samples", fontsize=10)
    else:
        ax.text(0.5, 0.5, "insufficient overlap", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
    fig.suptitle(f"G3 · Step 3 · Cross-sample consistency — {feature} [{mode}]", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = png_dir / f"fig03_{feature}_{mode}_per_sample.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)


def step4_stratified_auc(df: pd.DataFrame, feature: str, mode: str, png_dir: Path) -> list[dict]:
    sub = df[df["mode"] == mode].copy()
    if feature not in sub.columns:
        return []
    y_all = sub["tp_label"].values.astype(int)
    v_all = pd.to_numeric(sub[feature], errors="coerce").values
    if np.isnan(v_all).all():
        return []
    records: list[dict] = []
    # Global
    auc = compute_auc(y_all, v_all)
    n_pos = int((y_all == 1).sum())
    n_neg = int((y_all == 0).sum())
    auc_abs = max(auc, 1 - auc) if not np.isnan(auc) else float("nan")
    ci = wilson_auc_ci(auc_abs, n_pos, n_neg)
    records.append({"feature": feature, "mode": mode, "stratum": "Global", "label": "ALL",
                    "n": len(y_all), "n_tp": n_pos, "n_fp": n_neg,
                    "auc_raw": auc, "auc_abs": auc_abs, "ci_lo": ci[0], "ci_hi": ci[1]})
    # Per-LOH
    for lvl in sub["LOH_Subtype"].dropna().unique():
        m = sub["LOH_Subtype"] == lvl
        if m.sum() < 50:
            continue
        yy = y_all[m.values]
        vv = v_all[m.values]
        if len(set(yy.tolist())) < 2:
            continue
        a = compute_auc(yy, vv)
        aa = max(a, 1 - a) if not np.isnan(a) else float("nan")
        ci = wilson_auc_ci(aa, int((yy == 1).sum()), int((yy == 0).sum()))
        records.append({"feature": feature, "mode": mode, "stratum": "LOH", "label": lvl,
                        "n": int(m.sum()), "n_tp": int((yy == 1).sum()), "n_fp": int((yy == 0).sum()),
                        "auc_raw": a, "auc_abs": aa, "ci_lo": ci[0], "ci_hi": ci[1]})
    # Per-AF_class
    for lvl in sub["AF_class"].dropna().unique():
        m = sub["AF_class"] == lvl
        if m.sum() < 50:
            continue
        yy = y_all[m.values]
        vv = v_all[m.values]
        if len(set(yy.tolist())) < 2:
            continue
        a = compute_auc(yy, vv)
        aa = max(a, 1 - a) if not np.isnan(a) else float("nan")
        ci = wilson_auc_ci(aa, int((yy == 1).sum()), int((yy == 0).sum()))
        records.append({"feature": feature, "mode": mode, "stratum": "AF_class", "label": lvl,
                        "n": int(m.sum()), "n_tp": int((yy == 1).sum()), "n_fp": int((yy == 0).sum()),
                        "auc_raw": a, "auc_abs": aa, "ci_lo": ci[0], "ci_hi": ci[1]})
    # Per-CN
    for lvl in sub["cn_tier_F"].dropna().unique():
        m = sub["cn_tier_F"] == lvl
        if m.sum() < 50:
            continue
        yy = y_all[m.values]
        vv = v_all[m.values]
        if len(set(yy.tolist())) < 2:
            continue
        a = compute_auc(yy, vv)
        aa = max(a, 1 - a) if not np.isnan(a) else float("nan")
        ci = wilson_auc_ci(aa, int((yy == 1).sum()), int((yy == 0).sum()))
        records.append({"feature": feature, "mode": mode, "stratum": "CN_tier", "label": lvl,
                        "n": int(m.sum()), "n_tp": int((yy == 1).sum()), "n_fp": int((yy == 0).sum()),
                        "auc_raw": a, "auc_abs": aa, "ci_lo": ci[0], "ci_hi": ci[1]})
    # Per-FILTER (PASS only in master; report diagnostic)
    for lvl in sub["vcf_FILTER"].dropna().unique():
        m = sub["vcf_FILTER"] == lvl
        if m.sum() < 50:
            continue
        yy = y_all[m.values]
        vv = v_all[m.values]
        if len(set(yy.tolist())) < 2:
            continue
        a = compute_auc(yy, vv)
        aa = max(a, 1 - a) if not np.isnan(a) else float("nan")
        ci = wilson_auc_ci(aa, int((yy == 1).sum()), int((yy == 0).sum()))
        records.append({"feature": feature, "mode": mode, "stratum": "FILTER", "label": lvl,
                        "n": int(m.sum()), "n_tp": int((yy == 1).sum()), "n_fp": int((yy == 0).sum()),
                        "auc_raw": a, "auc_abs": aa, "ci_lo": ci[0], "ci_hi": ci[1]})
    return records


def plot_stratified_auc(records: list[dict], feature: str, mode: str, png_dir: Path):
    df_r = pd.DataFrame([r for r in records if r["feature"] == feature and r["mode"] == mode])
    if df_r.empty:
        return
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    strata_order = ["Global", "LOH", "AF_class", "CN_tier"]
    for i, strat in enumerate(strata_order):
        ax = axes[i]
        sub = df_r[df_r["stratum"] == strat]
        if sub.empty:
            ax.text(0.5, 0.5, "N/A", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(strat)
            ax.axis("off")
            continue
        labels = sub["label"].tolist()
        aucs = sub["auc_abs"].tolist()
        lo = sub["ci_lo"].tolist()
        hi = sub["ci_hi"].tolist()
        ypos = np.arange(len(labels))
        ax.barh(ypos, aucs, color="#3498db", alpha=0.7)
        for yp, a, l, h in zip(ypos, aucs, lo, hi):
            if not np.isnan(a):
                ax.plot([l, h], [yp, yp], color="black", lw=1)
        ax.axvline(0.5, color="gray", linestyle="--", lw=1)
        ax.axvline(0.58, color="red", linestyle=":", lw=1)
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(0.4, 1.0)
        ax.set_title(strat, fontsize=10)
        ax.set_xlabel("|AUC|")
    fig.suptitle(f"G3 · Step 4 · Stratified AUC — {feature} [{mode}]", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = png_dir / f"fig04_{feature}_{mode}_stratified_auc.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)


def step5_confound_guard(df: pd.DataFrame, feature: str, mode: str, png_dir: Path) -> dict:
    sub = df[df["mode"] == mode].copy()
    if feature not in sub.columns:
        return {}
    v = pd.to_numeric(sub[feature], errors="coerce").values
    if np.isnan(v).all():
        return {}
    # Build confounder matrix X (drop confounder == feature itself, e.g. vcf_AF)
    conf = [c for c in CONFOUNDERS if c != feature and c in sub.columns]
    if not conf:
        return {}
    X = sub[conf].apply(pd.to_numeric, errors="coerce").values
    y = sub["tp_label"].values.astype(int)
    raw_auc = compute_auc(y, v)
    raw_abs = max(raw_auc, 1 - raw_auc) if not np.isnan(raw_auc) else float("nan")
    resid = residualize(v, X)
    res_auc = compute_auc(y, resid)
    res_abs = max(res_auc, 1 - res_auc) if not np.isnan(res_auc) else float("nan")
    # AF-bin x-validation
    af_res = {}
    for lvl in ["Extreme", "Near-half", "Intermediate"]:
        m = sub["AF_class"].values == lvl
        if m.sum() < 50:
            continue
        a = compute_auc(y[m], v[m])
        aa = max(a, 1 - a) if not np.isnan(a) else float("nan")
        af_res[lvl] = aa
    # bar chart
    fig, ax = plt.subplots(figsize=(7, 4.5))
    cats = ["Raw (global)", "Residualized", *[f"AF:{k}" for k in af_res]]
    vals = [raw_abs, res_abs, *af_res.values()]
    colors = ["#3498db", "#2ecc71", *["#9b59b6"] * len(af_res)]
    ax.bar(cats, vals, color=colors, alpha=0.75)
    ax.axhline(0.5, color="gray", linestyle="--")
    ax.axhline(0.55, color="orange", linestyle=":")
    ax.axhline(0.58, color="red", linestyle=":")
    ax.set_ylim(0.4, 1.0)
    ax.set_ylabel("|AUC|")
    ax.set_title(f"G3 · Step 5 · Confound guard — {feature} [{mode}]\n(residualize on {','.join(conf)})", fontsize=10)
    for tick in ax.get_xticklabels():
        tick.set_rotation(20)
        tick.set_ha("right")
    fig.tight_layout()
    out = png_dir / f"fig05_{feature}_{mode}_confound.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return {
        "feature": feature,
        "mode": mode,
        "raw_auc_abs": raw_abs,
        "residualized_auc_abs": res_abs,
        **{f"af_{k}_auc": v for k, v in af_res.items()},
        "confounders": ",".join(conf),
        "confound_collapsed": (not np.isnan(res_abs)) and res_abs < 0.55 and raw_abs >= 0.58,
    }


def step6_spatial_autocorr(df: pd.DataFrame, feature: str, mode: str, png_dir: Path) -> dict:
    sub = df[df["mode"] == mode].copy()
    if feature not in sub.columns:
        return {}
    v = pd.to_numeric(sub[feature], errors="coerce").values
    if np.isnan(v).all():
        return {}
    sub = sub.assign(_f=v, _bin=((pd.to_numeric(sub["Pos"], errors="coerce") // 5_000_000).astype("Int64")))
    sub = sub.dropna(subset=["_f", "_bin"])
    if sub.empty:
        return {}
    agg = (
        sub.groupby(["Chr", "_bin"])
        .apply(lambda g: pd.Series({
            "n": len(g),
            "tp_rate": g["tp_label"].mean(),
            "auc": compute_auc(g["tp_label"].values.astype(int), g["_f"].values),
        }))
        .reset_index()
    )
    agg = agg.dropna(subset=["auc"])
    if agg.empty:
        return {}
    agg["auc_abs"] = agg["auc"].apply(lambda x: max(x, 1 - x))
    # scatter auc_abs vs tp_rate
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sc = ax.scatter(agg["tp_rate"], agg["auc_abs"], c=np.log10(agg["n"] + 1), s=14, cmap="viridis")
    ax.axhline(0.5, color="gray", linestyle="--")
    ax.axhline(0.58, color="red", linestyle=":")
    ax.axvline(0.8, color="orange", linestyle=":")
    ax.set_xlabel("baseline TP rate per 5Mb bin")
    ax.set_ylabel("|AUC| per bin")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
    cbar.set_label("log10(n+1)")
    frac_high_tp = ((agg["tp_rate"] > 0.8) & (agg["auc_abs"] >= 0.58)).sum()
    frac_all = (agg["auc_abs"] >= 0.58).sum()
    ratio = (frac_high_tp / frac_all) if frac_all > 0 else float("nan")
    ax.set_title(f"G3 · Step 6 · Spatial — {feature} [{mode}]\nhigh-TP artifact ratio = {ratio:.2f}", fontsize=10)
    fig.tight_layout()
    out = png_dir / f"fig06_{feature}_{mode}_spatial.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return {
        "feature": feature,
        "mode": mode,
        "n_bins": int(len(agg)),
        "bins_auc_ge_58": int(frac_all),
        "bins_auc_ge_58_high_tp": int(frac_high_tp),
        "high_tp_artifact_ratio": ratio,
    }


# ---------------------------------------------------------------------------
# Verdict PASS-subset diagnostic (paper question #2)
# ---------------------------------------------------------------------------

def verdict_pass_diagnostic(df: pd.DataFrame) -> pd.DataFrame:
    """Per-sample Verdict prevalence and TP rate in the PASS subset (TO-only).

    Master is PASS-only by construction, so Verdict != 0 is essentially absent.
    This function formalises that diagnostic and re-tests the 0420 Pilot conclusion
    across 7 samples.
    """
    to = df[df["mode"] == "to_pileup"].copy()
    rows = []
    for samp in SAMPLE_ORDER:
        s = to[to["sample"] == samp]
        if s.empty:
            continue
        n = len(s)
        tp = int((s["tp_label"] == 1).sum())
        for v in ["vcf_Verdict_Germline", "vcf_Verdict_Somatic", "vcf_Verdict_SubclonalSomatic"]:
            flagged = int((s[v].fillna(0) == 1).sum())
            flagged_tp = int(((s[v].fillna(0) == 1) & (s["tp_label"] == 1)).sum())
            flagged_fp = flagged - flagged_tp
            rows.append({
                "sample": samp,
                "verdict": v.replace("vcf_", ""),
                "n_PASS": n,
                "n_tp_PASS": tp,
                "n_tp_rate_PASS": tp / n if n else float("nan"),
                "n_verdict_flag": flagged,
                "pct_verdict_flag": flagged / n if n else float("nan"),
                "tp_rate_in_flagged": flagged_tp / flagged if flagged else float("nan"),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("[G3] Loading merged_with_vcf.tsv.gz ...")
    df = pd.read_csv(DATA, sep="\t", low_memory=False)
    # Build cn_tier_F
    df["cn_tier_F"] = make_cn_tier_F(pd.to_numeric(df["Coverage_Multiple"], errors="coerce"))
    # Force numeric on feature cols
    for c in ALL_FEATURES:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    logs_dir = ROOT / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 1 ----
    step1_rows: list[dict] = []
    for mode in ["paired_full", "to_pileup"]:
        feats = [f for f in ALL_FEATURES if f in df.columns
                 and df[df["mode"] == mode][f].notna().any()]
        for f in feats:
            s = step1_global_distribution(df, f, mode)
            if s:
                step1_rows.append(s)
        plot_global_distribution(df, feats, mode, FIG_DIR / f"fig01_global_{mode}.png", step1_rows)
    pd.DataFrame(step1_rows).to_csv(logs_dir / "g3_step1_global.tsv", sep="\t", index=False)
    print(f"[G3] Step 1 done: {len(step1_rows)} rows")

    # ---- Step 2 + 3 + 4 + 5 + 6 for "top" features ----
    # Select strongest per mode (by |AUC| absolute)
    s1 = pd.DataFrame(step1_rows)
    top_per_mode: dict[str, list[str]] = {}
    for mode in ["paired_full", "to_pileup"]:
        sub = s1[(s1["mode"] == mode) & s1["available"]]
        sub = sub.sort_values("auc_abs", ascending=False)
        # include top 4 + vcf_QUAL + vcf_AF + NAF(paired) + Verdict(to) for completeness
        pick = list(sub["feature"].head(6))
        for must in ["vcf_QUAL", "vcf_AF", "vcf_NAF", "vcf_Verdict_Germline", "vcf_SB"]:
            if must in sub["feature"].values and must not in pick:
                pick.append(must)
        top_per_mode[mode] = [f for f in pick if f in df.columns]

    step4_rows: list[dict] = []
    step5_rows: list[dict] = []
    step6_rows: list[dict] = []
    for mode, feats in top_per_mode.items():
        print(f"[G3] Deep steps for {mode}: {feats}")
        for f in feats:
            step2_strata_heatmaps(df, f, mode, FIG_DIR)
            step3_cross_sample(df, f, mode, FIG_DIR)
            recs = step4_stratified_auc(df, f, mode, FIG_DIR)
            step4_rows.extend(recs)
            plot_stratified_auc(recs, f, mode, FIG_DIR)
            r5 = step5_confound_guard(df, f, mode, FIG_DIR)
            if r5:
                step5_rows.append(r5)
            r6 = step6_spatial_autocorr(df, f, mode, FIG_DIR)
            if r6:
                step6_rows.append(r6)

    pd.DataFrame(step4_rows).to_csv(logs_dir / "g3_step4_stratified_auc.tsv", sep="\t", index=False)
    pd.DataFrame(step5_rows).to_csv(logs_dir / "g3_step5_confound.tsv", sep="\t", index=False)
    pd.DataFrame(step6_rows).to_csv(logs_dir / "g3_step6_spatial.tsv", sep="\t", index=False)

    # ---- Verdict PASS diagnostic ----
    verdict_df = verdict_pass_diagnostic(df)
    verdict_df.to_csv(logs_dir / "g3_verdict_pass_diagnostic.tsv", sep="\t", index=False)
    # barplot
    if not verdict_df.empty:
        fig, ax = plt.subplots(figsize=(10, 5))
        piv = verdict_df.pivot(index="sample", columns="verdict", values="n_verdict_flag").reindex(SAMPLE_ORDER)
        piv.plot(kind="bar", ax=ax, color=["#e74c3c", "#2ecc71", "#9b59b6"])
        ax.set_ylabel("# Verdict flag == 1 within PASS subset")
        ax.set_title("G3 · Verdict PASS-subset prevalence (TO)\n(0 across all 7 samples = replicates 0420 Pilot conclusion)")
        ax.legend(loc="upper right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig07_verdict_pass_prevalence.png", dpi=120)
        plt.close(fig)

    # ---- save consolidated summary ----
    summary = {
        "step1_feature_count": len(step1_rows),
        "step4_rows": len(step4_rows),
        "step5_rows": len(step5_rows),
        "step6_rows": len(step6_rows),
        "top_per_mode": top_per_mode,
    }
    (logs_dir / "g3_summary.json").write_text(json.dumps(summary, indent=2))
    print("[G3] Summary:", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
