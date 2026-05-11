#!/usr/bin/env python3
"""Phase C · G10 — Quality Summary & Verification feature group observation.

Scope (27 registered features, 24 available from canonical runs):
  Composite score:       HeuristicScore, Quality_Score, Quality_Tier, Stability
  Verification outcome:  VerificationClass, DominantLabel, Significant,
                         SuggestFilter, PassedGating
  PERMANOVA label:       LabelHPPermanovaF/P/Valid, LabelHPDispersionP/Warn,
                         LabelAllelePermanovaF/P/Valid, LabelAlleleDispersionP/Warn
  Unassigned affinity:   UnassignedAffinity, UnassignedAffinityP, UnassignedDir,
                         NHP0, NHP3
  UNAVAILABLE in canonical (post-fix addition):
    NHP_Somatic11, NHP_Somatic21, NHP_Somatic33 — will be reported as missing

Special handling per task spec:
  - Quality_Score: suspected TO AUC=0.497 (QS redesign evidence).  Report per-
    sample AUC across 7 samples × 2 modes.
  - LabelAllelePermanovaF: suspected AF proxy — residualize on vcf_AF and
    re-AUC.  If resid AUC falls ≤0.53, mark CONFOUND_COLLAPSED.
  - VerificationClass: 4 categorical (Strong/Subclone/Weak/Noise).  Compute
    per-class TP rate × 7 samples; check cross-sample stability.
  - SuggestFilter (binary): compute precision/recall at the filter decision.
  - NHP0/NHP3 boundary: fraction of reads unassigned/conflicting.
  - PassedGating: binary gate, per-class TP rate.

Input:
  data/G10/master_g10.tsv.gz  (built inline by _build_master() if missing)
Output:
  features/G10_quality_verification.md        (written separately by caller)
  figures/G10_quality/fig01..fig09            PNG panels
  data/G10/G10_*.tsv                          tables
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
from sklearn.linear_model import LinearRegression

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------
# Constants
ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
BIG7 = Path("/big7_disk/liaoyoyo2001/big7_disk_output")
MASTER_IN = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
OUT_DATA = ROOT / "research/feature_layered_observation/data/G10"
OUT_FIG = ROOT / "research/feature_layered_observation/figures/G10_quality"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)
MASTER_G10 = OUT_DATA / "master_g10.tsv.gz"

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]
LOH_ORDER = ["LOH_None", "LOH_Weak", "LOH_Noise", "LOH_Strong", "LOH_Subclone"]
AF_ORDER = ["Extreme", "Near-half", "Intermediate"]
CN_BREAKS = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
CN_LABELS = ["CN_Loss", "CN_Near1", "CN_Diploid", "CN_Gain", "CN_HighGain"]

# Feature groups (by dtype for panels)
FEATURES_CONT = [
    "HeuristicScore", "Quality_Score", "Stability",
    "LabelHPPermanovaF", "LabelHPPermanovaP", "LabelHPDispersionP",
    "LabelAllelePermanovaF", "LabelAllelePermanovaP", "LabelAlleleDispersionP",
    "UnassignedAffinity", "UnassignedAffinityP",
    "NHP0", "NHP3",
]
FEATURES_BIN = [
    "PassedGating", "Significant", "SuggestFilter",
    "LabelHPPermanovaValid", "LabelHPDispersionWarn",
    "LabelAllelePermanovaValid", "LabelAlleleDispersionWarn",
]
FEATURES_CAT = ["VerificationClass", "DominantLabel", "UnassignedDir", "Quality_Tier"]

# Run path maps
PAIRED_FULL = {
    "HCC1395":        BIG7 / "canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO": BIG7 / "canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full",
    "HCC1937":        BIG7 / "canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full",
    "HCC1954":        BIG7 / "canonical/HCC1954/paired_full/20260421_HCC1954_paired_full_full",
    "H2009":          BIG7 / "canonical/H2009/paired_full/20260421_H2009_paired_full_full",
    "H1437":          BIG7 / "canonical/H1437/paired_full/20260421_H1437_paired_full_full",
    "COLO829":        BIG7 / "canonical/COLO829/paired_full/20260421_COLO829_paired_full_full",
}
TO_SOURCES = {
    "HCC1395":        BIG7 / "bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/step05_intersubmod",
    "HCC1395_DORADO": BIG7 / "synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod",
    "HCC1937":        BIG7 / "synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod",
    "HCC1954":        BIG7 / "synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod",
    "H2009":          BIG7 / "synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod",
    "H1437":          BIG7 / "synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod",
}

# All G10 cols in significance_summary.csv (not including NHP_Somatic11/21/33)
G10_COLS_AVAILABLE = [
    "HeuristicScore", "PassedGating",
    "LabelHPPermanovaF", "LabelHPPermanovaP", "LabelHPPermanovaValid",
    "LabelHPDispersionP", "LabelHPDispersionWarn",
    "LabelAllelePermanovaF", "LabelAllelePermanovaP", "LabelAllelePermanovaValid",
    "LabelAlleleDispersionP", "LabelAlleleDispersionWarn",
    "UnassignedAffinity", "UnassignedAffinityP", "UnassignedDir",
    "NHP3", "NHP0",
    "DominantLabel", "Stability", "VerificationClass",
    "Quality_Score", "Quality_Tier",
    "Significant", "SuggestFilter",
]


# ----------------------------------------------------------------------------
# Helpers

def auc_wilson(y: np.ndarray, s: np.ndarray):
    mask = ~np.isnan(s) & ~pd.isna(y)
    y = y[mask].astype(int); s = s[mask]
    n1 = int((y == 1).sum()); n0 = int((y == 0).sum())
    if n1 < 5 or n0 < 5:
        return np.nan, np.nan, np.nan, n1, n0
    r = stats.rankdata(s)
    auc = (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    q1 = auc / (2 - auc); q2 = 2 * auc * auc / (1 + auc)
    var = (auc * (1 - auc) + (n1 - 1) * (q1 - auc**2) + (n0 - 1) * (q2 - auc**2)) / (n1 * n0)
    se = np.sqrt(max(var, 0))
    return auc, max(0., auc - 1.96 * se), min(1., auc + 1.96 * se), n1, n0


def cohen_d(pos, neg):
    pos = np.asarray(pos, float); pos = pos[~np.isnan(pos)]
    neg = np.asarray(neg, float); neg = neg[~np.isnan(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return np.nan
    s2 = ((len(pos)-1)*pos.var(ddof=1) + (len(neg)-1)*neg.var(ddof=1)) / (len(pos)+len(neg)-2)
    if s2 <= 0:
        return 0.0
    return (pos.mean() - neg.mean()) / np.sqrt(s2)


def mwu_p(pos, neg):
    pos = np.asarray(pos, float); pos = pos[~np.isnan(pos)]
    neg = np.asarray(neg, float); neg = neg[~np.isnan(neg)]
    if len(pos) < 5 or len(neg) < 5:
        return np.nan
    try:
        return float(stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue)
    except Exception:
        return np.nan


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    denom = 1 + z*z/n
    center = (p + z*z/(2*n)) / denom
    half = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return p, max(0., center - half), min(1., center + half)


# ----------------------------------------------------------------------------
# Master build

def _read_ss(sample: str, mode: str, run_dir: Path) -> pd.DataFrame:
    frames = []
    for tp, subdir in [(1, "intersubmod_tp"), (0, "intersubmod_fp")]:
        ss = run_dir / subdir / "significance_summary.csv"
        if not ss.exists():
            print(f"[warn] missing {ss}", flush=True)
            continue
        try:
            header = pd.read_csv(ss, nrows=0).columns.tolist()
            keep = ["Chr", "Pos"] + [c for c in G10_COLS_AVAILABLE if c in header]
            d = pd.read_csv(ss, usecols=keep, low_memory=False)
        except Exception as e:
            print(f"[err] {ss}: {e}", flush=True); continue
        d["sample"] = sample; d["mode"] = mode; d["_ss_tp"] = tp
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_master() -> pd.DataFrame:
    print(f"[load] master {MASTER_IN}", flush=True)
    master = pd.read_csv(MASTER_IN, sep="\t", low_memory=False)
    print(f"[load] rows={len(master):,} cols={master.shape[1]}", flush=True)
    frames = []
    for s, p in PAIRED_FULL.items():
        if not p.exists():
            print(f"[skip] paired_full {s}: {p}", flush=True); continue
        df = _read_ss(s, "paired_full", p)
        if not df.empty:
            frames.append(df)
            print(f"[ok]   paired_full {s}: rows={len(df):,}", flush=True)
    for s, p in TO_SOURCES.items():
        if not p.exists():
            print(f"[skip] to_pileup {s}: {p}", flush=True); continue
        df = _read_ss(s, "to_pileup", p)
        if not df.empty:
            frames.append(df)
            print(f"[ok]   to_pileup  {s}: rows={len(df):,}", flush=True)
    if not frames:
        print("[err] no G10 data collected; abort"); sys.exit(1)
    g10 = pd.concat(frames, ignore_index=True)
    g10["Chr"] = g10["Chr"].astype(str)
    master["Chr"] = master["Chr"].astype(str)
    g10["Pos"] = pd.to_numeric(g10["Pos"], errors="coerce").astype("Int64")
    master["Pos"] = pd.to_numeric(master["Pos"], errors="coerce").astype("Int64")
    g10 = g10.drop_duplicates(subset=["sample", "mode", "Chr", "Pos"], keep="first") \
             .drop(columns=["_ss_tp"], errors="ignore")
    print(f"[g10]  unique rows={len(g10):,}", flush=True)
    merged = master.merge(g10, on=["sample", "mode", "Chr", "Pos"], how="left",
                          suffixes=("", "_ss"))
    # Coverage stats
    for c in G10_COLS_AVAILABLE:
        if c in merged.columns:
            n_ok = merged[c].notna().sum()
            pct = n_ok / len(merged) * 100
            print(f"[cov]  {c}: {n_ok:,}/{len(merged):,} ({pct:.1f}%)", flush=True)
    merged.to_csv(MASTER_G10, sep="\t", index=False, compression="gzip")
    print(f"[save] {MASTER_G10}", flush=True)
    return merged


def load_master() -> pd.DataFrame:
    if MASTER_G10.exists():
        print(f"[load] {MASTER_G10}", flush=True)
        df = pd.read_csv(MASTER_G10, sep="\t", low_memory=False)
    else:
        df = build_master()
    df["tp_label"] = pd.to_numeric(df["tp_label"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["tp_label"]).copy()
    df["tp_label"] = df["tp_label"].astype(int)
    df["cn_tier_F"] = pd.cut(df["Coverage_Multiple"].astype(float),
                             bins=CN_BREAKS, labels=CN_LABELS, right=True)
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("LOH_None")
    df["LOH_Subtype"] = pd.Categorical(df["LOH_Subtype"], categories=LOH_ORDER, ordered=True)
    df["AF_class"] = pd.Categorical(df["AF_class"].fillna("Unknown"),
                                    categories=AF_ORDER + ["Unknown"])

    # Booleanize string/int flags
    for c in ["PassedGating", "Significant", "SuggestFilter",
              "LabelHPPermanovaValid", "LabelHPDispersionWarn",
              "LabelAllelePermanovaValid", "LabelAlleleDispersionWarn"]:
        if c in df.columns:
            df[c] = df[c].map({"true": 1, "True": 1, "1": 1, 1: 1, True: 1,
                               "false": 0, "False": 0, "0": 0, 0: 0, False: 0})
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


# ----------------------------------------------------------------------------
# Step 1: global distribution panel

def step1_global(df: pd.DataFrame) -> pd.DataFrame:
    feats = FEATURES_CONT + FEATURES_BIN
    rows = []
    n_cols = 5
    n_rows = int(np.ceil(len(feats) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.3*n_cols, 3.8*n_rows))
    axes = axes.flatten()
    for i, feat in enumerate(feats):
        ax = axes[i]
        if feat not in df.columns:
            ax.set_title(f"{feat}\n(MISSING)", fontsize=9)
            ax.axis("off"); continue
        vals = pd.to_numeric(df[feat], errors="coerce").to_numpy(float)
        y = df["tp_label"].to_numpy(int)
        tp = vals[y == 1]; fp = vals[y == 0]
        tp = tp[~np.isnan(tp)]; fp = fp[~np.isnan(fp)]
        if len(tp) < 10 or len(fp) < 10:
            ax.set_title(f"{feat}\n(insufficient n)", fontsize=9); ax.axis("off"); continue

        if feat in FEATURES_BIN:
            rate_tp = float(tp.mean()); rate_fp = float(fp.mean())
            ax.bar([0, 1], [rate_fp, rate_tp], color=["#E07A5F", "#2E86AB"], edgecolor="black")
            ax.set_xticks([0, 1]); ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"])
            ax.set_ylim(0, 1.0); ax.set_ylabel(f"{feat} rate")
            auc, lo, hi, _, _ = auc_wilson(y, vals)
            ax.set_title(f"{feat}\nrateTP={rate_tp:.3f} rateFP={rate_fp:.3f}\nAUC={auc:.3f}", fontsize=9)
            rows.append(dict(feature=feat, dtype="binary", auc=auc, auc_lo=lo, auc_hi=hi,
                             cohen_d=np.nan, mwu_p=np.nan,
                             mean_tp=rate_tp, mean_fp=rate_fp,
                             n_tp=len(tp), n_fp=len(fp)))
            continue
        log_scale = feat in {"NHP0", "NHP3"}
        p_tp = np.log10(tp + 1) if log_scale else tp
        p_fp = np.log10(fp + 1) if log_scale else fp
        try:
            parts = ax.violinplot([p_fp, p_tp], positions=[0, 1], widths=0.7,
                                  showmeans=False, showmedians=True)
            for pc, color in zip(parts["bodies"], ["#E07A5F", "#2E86AB"]):
                pc.set_facecolor(color); pc.set_alpha(0.6)
        except Exception:
            pass
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"], fontsize=8)
        auc, lo, hi, _, _ = auc_wilson(y, vals)
        d = cohen_d(tp, fp); p = mwu_p(tp, fp)
        ax.set_title(f"{feat}{' (log10)' if log_scale else ''}\nAUC={auc:.3f}[{lo:.2f},{hi:.2f}] d={d:.2f}\nmwu-p={p:.1e}",
                     fontsize=9)
        rows.append(dict(feature=feat, dtype="cont", auc=auc, auc_lo=lo, auc_hi=hi,
                         cohen_d=d, mwu_p=p,
                         mean_tp=float(np.nanmean(tp)), mean_fp=float(np.nanmean(fp)),
                         n_tp=len(tp), n_fp=len(fp)))
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle("G10 Quality/Verification · Step 1 Global TP vs FP distribution "
                 "(7 samples × 2 modes pooled)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG / "fig01_global_distribution.png", dpi=130)
    plt.close(fig)
    out = pd.DataFrame(rows).sort_values("auc", ascending=False)
    out.to_csv(OUT_DATA / "G10_global_stats.tsv", sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
# Step 2: LOH × AF × CN 4-panel heatmap for flagship features

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

    fig, axes = plt.subplots(1, 4, figsize=(24, 10))
    im0 = axes[0].imshow(piv_rate.where(mask_ok).values, aspect="auto",
                         cmap="RdBu_r", vmin=0, vmax=1)
    axes[0].set_title(f"A · TP rate (n≥{min_n})")
    fig.colorbar(im0, ax=axes[0], fraction=0.03)
    im1 = axes[1].imshow(piv_tp.where(mask_ok).values, aspect="auto", cmap="Blues")
    axes[1].set_title(f"B · median({feat}) | TP")
    fig.colorbar(im1, ax=axes[1], fraction=0.03)
    im2 = axes[2].imshow(piv_fp.where(mask_ok).values, aspect="auto", cmap="Oranges")
    axes[2].set_title(f"C · median({feat}) | FP")
    fig.colorbar(im2, ax=axes[2], fraction=0.03)
    dlt = piv_delta.values.astype(float)
    finite = dlt[np.isfinite(dlt)]
    if len(finite):
        vlim = max(abs(np.nanpercentile(finite, 5)), abs(np.nanpercentile(finite, 95)), 1e-6)
        norm = TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim)
        im3 = axes[3].imshow(dlt, aspect="auto", cmap="coolwarm", norm=norm)
    else:
        im3 = axes[3].imshow(dlt, aspect="auto", cmap="coolwarm")
    axes[3].set_title("D · Δ(TP−FP) median")
    fig.colorbar(im3, ax=axes[3], fraction=0.03)
    for ax in axes:
        ax.set_yticks(range(len(row_order)))
        ax.set_yticklabels(row_order, fontsize=7)
        ax.set_xticks(range(len(col_order)))
        ax.set_xticklabels(col_order, rotation=30, ha="right", fontsize=8)
    for i, r in enumerate(row_order):
        for j, c in enumerate(col_order):
            if r in piv_n.index and c in piv_n.columns:
                n = piv_n.loc[r, c]
                if pd.isna(n): n = 0
                if n >= min_n:
                    axes[0].text(j, i, f"{int(n):,}", ha="center", va="center", fontsize=6)
    fig.suptitle(f"{feat} · G10 Step 2 LOH × AF × CN heatmap", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_FIG / f"fig02_{feat}_heatmap.png", dpi=130)
    plt.close(fig)
    flat = piv_delta.reset_index().melt(id_vars="row_key", var_name="cn_tier",
                                        value_name="delta_tp_fp")
    flat["feature"] = feat
    return flat


# ----------------------------------------------------------------------------
# Step 3: per-sample consistency grid + VerificationClass × 7 samples

def step3_per_sample(df: pd.DataFrame) -> pd.DataFrame:
    """Quality_Score AUC per (sample, mode) + VerificationClass per-class TP rate."""
    rows = []
    samples = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    for s in samples:
        for mode in ["paired_full", "to_pileup"]:
            sub = df[(df["sample"] == s) & (df["mode"] == mode)]
            if len(sub) < 100:
                continue
            y = sub["tp_label"].to_numpy(int)
            for feat in ["Quality_Score", "HeuristicScore", "Stability",
                         "LabelAllelePermanovaF", "LabelHPPermanovaF",
                         "UnassignedAffinity", "NHP0", "NHP3"]:
                if feat not in sub.columns:
                    continue
                v = pd.to_numeric(sub[feat], errors="coerce").to_numpy(float)
                auc, lo, hi, n1, n0 = auc_wilson(y, v)
                rows.append(dict(sample=s, mode=mode, feature=feat,
                                 auc=auc, auc_lo=lo, auc_hi=hi,
                                 n_tp=n1, n_fp=n0, tp_rate=n1/(n1+n0) if (n1+n0) else np.nan))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G10_persample_auc.tsv", sep="\t", index=False)

    # Figure: heatmap AUC × sample × mode for each feature row
    feats_plot = ["Quality_Score", "HeuristicScore", "Stability",
                  "LabelAllelePermanovaF", "LabelHPPermanovaF",
                  "UnassignedAffinity", "NHP0", "NHP3"]
    mat = np.full((len(feats_plot), len(samples) * 2), np.nan)
    col_labels = []
    for mi, mode in enumerate(["paired_full", "to_pileup"]):
        for si, s in enumerate(samples):
            col_labels.append(f"{s}\n{mode.split('_')[0]}")
            for fi, feat in enumerate(feats_plot):
                sel = out[(out["sample"] == s) & (out["mode"] == mode) & (out["feature"] == feat)]
                if not sel.empty:
                    mat[fi, mi * len(samples) + si] = sel["auc"].iloc[0]
    fig, ax = plt.subplots(figsize=(16, 6))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0.40, vmax=0.80)
    ax.set_yticks(range(len(feats_plot)))
    ax.set_yticklabels(feats_plot, fontsize=9)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                        fontsize=7, color="black")
    ax.axvline(len(samples) - 0.5, ls="--", c="black", lw=1.5)
    ax.set_title("G10 Step 3 · Per-sample AUC (left: paired_full / right: to_pileup)\n"
                 "Colorbar 0.40-0.80, dashed line = paired|TO split", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.02)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig03_persample_auc_heatmap.png", dpi=130)
    plt.close(fig)
    return out


# ----------------------------------------------------------------------------
# Step 4: stratified AUC for flagship features

def step4_stratified(df: pd.DataFrame, feats: List[str]) -> pd.DataFrame:
    all_rows = []
    for feat in feats:
        if feat not in df.columns:
            continue
        y = df["tp_label"].to_numpy(int)
        s = pd.to_numeric(df[feat], errors="coerce").to_numpy(float)
        rows = []
        rows.append(("global", "all", *auc_wilson(y, s)))
        for loh in LOH_ORDER:
            m = (df["LOH_Subtype"] == loh).to_numpy()
            if m.sum() >= 100:
                rows.append(("LOH", loh, *auc_wilson(y[m], s[m])))
        for af in AF_ORDER:
            m = (df["AF_class"] == af).to_numpy()
            if m.sum() >= 100:
                rows.append(("AF", af, *auc_wilson(y[m], s[m])))
        for cn in CN_LABELS:
            m = (df["cn_tier_F"] == cn).to_numpy()
            if m.sum() >= 100:
                rows.append(("CN", cn, *auc_wilson(y[m], s[m])))
        for mode in ["paired_full", "to_pileup"]:
            m = (df["mode"] == mode).to_numpy()
            if m.sum() >= 100:
                rows.append(("mode", mode, *auc_wilson(y[m], s[m])))
        dfres = pd.DataFrame(rows, columns=["layer", "group", "auc", "lo", "hi", "n_pos", "n_neg"])
        dfres["feature"] = feat
        all_rows.append(dfres)
    out = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    out.to_csv(OUT_DATA / "G10_auc_table.tsv", sep="\t", index=False)

    # Figure: grid of bar charts
    n_feats = len(feats)
    n_cols = 2
    n_rows = int(np.ceil(n_feats / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4.5 * n_rows))
    axes = axes.flatten()
    colors = {"global": "#333", "LOH": "#2E86AB", "AF": "#F4A261",
              "CN": "#6A994E", "mode": "#8E44AD"}
    for i, feat in enumerate(feats):
        ax = axes[i]
        sub = out[out["feature"] == feat]
        if sub.empty:
            ax.axis("off"); continue
        for j, (_, row) in enumerate(sub.iterrows()):
            lo = max(0., row["auc"] - row["lo"]) if pd.notna(row["lo"]) else 0
            hi = max(0., row["hi"] - row["auc"]) if pd.notna(row["hi"]) else 0
            ax.bar(j, row["auc"], color=colors.get(row["layer"], "gray"),
                   yerr=[[lo], [hi]], capsize=3, edgecolor="black", linewidth=0.4)
            if pd.notna(row["auc"]):
                ax.text(j, row["auc"] + 0.01, f"{row['auc']:.2f}", ha="center", fontsize=6)
        ax.axhline(0.5, ls="--", c="gray", lw=0.7)
        ax.axhline(0.58, ls=":", c="red", lw=0.7)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels([f"{r.layer}:{r.group}" for _, r in sub.iterrows()],
                           rotation=45, ha="right", fontsize=7)
        ax.set_ylim(0.3, 1.0)
        ax.set_title(f"{feat} stratified AUC", fontsize=10)
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle("G10 Step 4 · Stratified AUC (global / LOH / AF / CN / mode)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG / "fig04_stratified_auc.png", dpi=130)
    plt.close(fig)
    return out


# ----------------------------------------------------------------------------
# Step 5: Confound guard — LabelAllelePermanovaF residualized on vcf_AF

def step5_confound(df: pd.DataFrame, feats: List[str]) -> pd.DataFrame:
    rows = []
    covars_default = ["NumReads", "vcf_AF", "Coverage_Multiple", "AlleleDelta"]

    n_feats = len(feats)
    fig, axes = plt.subplots(n_feats, 2, figsize=(14, 4.5 * n_feats))
    if n_feats == 1:
        axes = axes[np.newaxis, :]
    for fi, feat in enumerate(feats):
        if feat not in df.columns:
            axes[fi, 0].axis("off"); axes[fi, 1].axis("off"); continue
        covars = [c for c in covars_default if c in df.columns and c != feat]
        sub = df.dropna(subset=[feat] + covars + ["tp_label"]).copy()
        if len(sub) < 200:
            axes[fi, 0].axis("off"); axes[fi, 1].axis("off"); continue
        X = sub[covars].astype(float).to_numpy()
        y_feat = pd.to_numeric(sub[feat], errors="coerce").to_numpy(float)
        y_tp = sub["tp_label"].to_numpy(int)
        try:
            lr = LinearRegression().fit(X, y_feat)
            resid = y_feat - lr.predict(X)
        except Exception:
            axes[fi, 0].axis("off"); axes[fi, 1].axis("off"); continue
        raw_auc = auc_wilson(y_tp, y_feat)[0]
        res_auc = auc_wilson(y_tp, resid)[0]
        af_bin = {}
        for af in AF_ORDER:
            m = (sub["AF_class"] == af).to_numpy()
            if m.sum() >= 200:
                af_bin[af] = auc_wilson(y_tp[m], y_feat[m])[0]
        verdict = "COLLAPSED" if res_auc < 0.53 else ("ATTENUATED" if abs(raw_auc - res_auc) > 0.03 else "STABLE")
        rows.append(dict(feature=feat, raw_auc=raw_auc, resid_auc=res_auc, verdict=verdict,
                         covars=",".join(covars),
                         **{f"auc_{k}": v for k, v in af_bin.items()}))
        # Plot
        bins = 50
        ax = axes[fi, 0]
        ax.hist(y_feat[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F", density=True, label="FP")
        ax.hist(y_feat[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB", density=True, label="TP")
        ax.set_title(f"raw {feat}  AUC={raw_auc:.3f}")
        ax.legend(fontsize=8); ax.set_xlabel(feat)
        ax2 = axes[fi, 1]
        ax2.hist(resid[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F", density=True, label="FP")
        ax2.hist(resid[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB", density=True, label="TP")
        ax2.set_title(f"resid on {covars}\nAUC={res_auc:.3f} verdict={verdict}\n"
                      f"AF-bin AUC: " + ", ".join(f"{k}={v:.2f}" for k, v in af_bin.items()),
                      fontsize=8)
        ax2.legend(fontsize=8)
    fig.suptitle("G10 Step 5 · Confound guard (OLS residualization on NumReads/vcf_AF/CovM/AlleleDelta)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG / "fig05_confound_residualized.png", dpi=130)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G10_confound.tsv", sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
# Step 6: Spatial auto-correlation

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
        y = g["tp_label"].to_numpy(int)
        if y.sum() < 5 or (len(y) - y.sum()) < 5:
            continue
        a, *_ = auc_wilson(y, pd.to_numeric(g[feat], errors="coerce").to_numpy(float))
        rows.append(dict(bin=key, n=len(g), tp_rate=y.mean(), auc=a))
    if not rows:
        return None
    s = pd.DataFrame(rows)
    s.to_csv(OUT_DATA / f"G10_spatial_{feat}.tsv", sep="\t", index=False)
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
        diff = high_tp["auc"].mean() - low_tp["auc"].mean()
        flag = f"artifact-suspect (ΔAUC={diff:+.3f})" if diff > 0.08 else f"ok (ΔAUC={diff:+.3f})"
    fig.suptitle(f"{feat} · G10 Step 6 spatial ({bin_mb} Mb) {flag}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_FIG / f"fig06_{feat}_spatial.png", dpi=130)
    plt.close(fig)
    return dict(feature=feat, n_bins=len(s), flag=flag)


# ----------------------------------------------------------------------------
# Step 7: VerificationClass × sample TP rate (categorical consistency)

def step7_verification_class(df: pd.DataFrame) -> pd.DataFrame:
    if "VerificationClass" not in df.columns:
        return None
    sub = df.dropna(subset=["VerificationClass", "tp_label"]).copy()
    classes = sorted(sub["VerificationClass"].astype(str).unique().tolist())
    samples = [s for s in SAMPLE_ORDER if s in sub["sample"].unique()]
    rows = []
    for s in samples:
        for mode in ["paired_full", "to_pileup"]:
            g = sub[(sub["sample"] == s) & (sub["mode"] == mode)]
            if len(g) < 50:
                continue
            for cls in classes:
                m = g["VerificationClass"].astype(str) == cls
                y = g.loc[m, "tp_label"].to_numpy(int)
                n = len(y); k = int(y.sum())
                p, lo, hi = wilson_ci(k, n)
                rows.append(dict(sample=s, mode=mode, verification_class=cls,
                                 n=n, k_tp=k, tp_rate=p, lo=lo, hi=hi))
    # Global pooled
    for cls in classes:
        m = sub["VerificationClass"].astype(str) == cls
        y = sub.loc[m, "tp_label"].to_numpy(int)
        n = len(y); k = int(y.sum())
        p, lo, hi = wilson_ci(k, n)
        rows.append(dict(sample="POOLED", mode="all", verification_class=cls,
                         n=n, k_tp=k, tp_rate=p, lo=lo, hi=hi))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G10_verification_class.tsv", sep="\t", index=False)

    # Figure: grid of class × sample
    pv = out[out["sample"] != "POOLED"].pivot_table(
        index="verification_class", columns=["sample", "mode"],
        values="tp_rate", aggfunc="first"
    )
    # reorder columns
    pv = pv.reindex(columns=pd.MultiIndex.from_product([samples, ["paired_full", "to_pileup"]]),
                    fill_value=np.nan)
    fig, axes = plt.subplots(1, 2, figsize=(20, 6))
    ax = axes[0]
    im = ax.imshow(pv.values, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
    ax.set_yticks(range(len(pv.index))); ax.set_yticklabels(pv.index, fontsize=9)
    ax.set_xticks(range(len(pv.columns)))
    ax.set_xticklabels([f"{a}/{b.split('_')[0]}" for a, b in pv.columns],
                       rotation=45, ha="right", fontsize=8)
    for i in range(pv.shape[0]):
        for j in range(pv.shape[1]):
            v = pv.values[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("VerificationClass TP rate by (sample, mode)")
    fig.colorbar(im, ax=ax, fraction=0.03)
    ax2 = axes[1]
    pooled = out[out["sample"] == "POOLED"]
    ax2.barh(range(len(pooled)), pooled["tp_rate"], color="#6A994E",
             xerr=[(pooled["tp_rate"] - pooled["lo"]).values,
                   (pooled["hi"] - pooled["tp_rate"]).values],
             edgecolor="black")
    ax2.set_yticks(range(len(pooled)))
    ax2.set_yticklabels(pooled["verification_class"])
    for i, (_, row) in enumerate(pooled.iterrows()):
        ax2.text(row["tp_rate"] + 0.01, i,
                 f"{row['tp_rate']:.3f}\nn={int(row['n']):,}", va="center", fontsize=8)
    ax2.set_xlim(0, 1.0)
    ax2.axvline(0.5, ls="--", c="gray")
    ax2.set_xlabel("TP rate (Wilson 95% CI)")
    ax2.set_title("VerificationClass pooled TP rate")
    fig.suptitle("G10 Step 7 · VerificationClass consistency", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_FIG / "fig07_verification_class.png", dpi=130)
    plt.close(fig)
    return out


# ----------------------------------------------------------------------------
# Step 8: SuggestFilter precision/recall
#   SuggestFilter=1 means "flagged as likely FP (label_delta > 0.3)"
#   In F1-filter context: treat the filter as predicting FP (tp_label=0).
#   Compute precision = P(FP | SuggestFilter=1) and recall = P(SuggestFilter=1 | FP).

def step8_suggest_filter(df: pd.DataFrame) -> pd.DataFrame:
    if "SuggestFilter" not in df.columns:
        return None
    sub = df.dropna(subset=["SuggestFilter", "tp_label"]).copy()
    rows = []

    def prec_rec(g):
        n = len(g); k_flag = int((g["SuggestFilter"] == 1).sum())
        n_fp = int((g["tp_label"] == 0).sum()); n_tp = int((g["tp_label"] == 1).sum())
        if k_flag == 0 or n_fp == 0:
            return dict(n=n, n_flag=k_flag, n_fp=n_fp, n_tp=n_tp,
                        precision=np.nan, recall=np.nan,
                        tp_loss=np.nan, flag_rate=k_flag / n if n else np.nan)
        tp_in_flag = int(((g["SuggestFilter"] == 1) & (g["tp_label"] == 0)).sum())
        fp_in_flag = tp_in_flag  # same count
        tp_lost = int(((g["SuggestFilter"] == 1) & (g["tp_label"] == 1)).sum())
        precision_fp = fp_in_flag / k_flag
        recall_fp = fp_in_flag / n_fp
        return dict(n=n, n_flag=k_flag, n_fp=n_fp, n_tp=n_tp,
                    precision=precision_fp, recall=recall_fp,
                    tp_loss_rate=tp_lost / n_tp if n_tp else np.nan,
                    flag_rate=k_flag / n)
    # Global pooled
    rows.append(dict(stratum="GLOBAL", **prec_rec(sub)))
    # Per mode
    for mode in ["paired_full", "to_pileup"]:
        g = sub[sub["mode"] == mode]
        if len(g) < 50:
            continue
        rows.append(dict(stratum=f"mode={mode}", **prec_rec(g)))
    # Per sample
    for s in SAMPLE_ORDER:
        for mode in ["paired_full", "to_pileup"]:
            g = sub[(sub["sample"] == s) & (sub["mode"] == mode)]
            if len(g) < 50:
                continue
            rows.append(dict(stratum=f"{s}/{mode}", **prec_rec(g)))
    # Per LOH
    for loh in LOH_ORDER:
        g = sub[sub["LOH_Subtype"] == loh]
        if len(g) < 50:
            continue
        rows.append(dict(stratum=f"LOH={loh}", **prec_rec(g)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G10_suggestfilter_prec_rec.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    ax = axes[0]
    plotted = out.dropna(subset=["precision"]).copy()
    x = np.arange(len(plotted))
    ax.bar(x - 0.2, plotted["precision"], width=0.38,
           color="#2E86AB", label="precision (FP among flagged)", edgecolor="black")
    ax.bar(x + 0.2, plotted["recall"], width=0.38,
           color="#E07A5F", label="recall (FP captured)", edgecolor="black")
    ax.axhline(0.5, ls="--", c="gray", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(plotted["stratum"], rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0, 1.0)
    ax.set_title("SuggestFilter as FP predictor (precision & recall, FP=positive)")
    ax.legend(fontsize=8)
    for xi, (_, row) in zip(x, plotted.iterrows()):
        ax.text(xi - 0.2, row["precision"] + 0.01, f"{row['precision']:.2f}",
                ha="center", fontsize=6)
        ax.text(xi + 0.2, row["recall"] + 0.01, f"{row['recall']:.2f}",
                ha="center", fontsize=6)
    ax2 = axes[1]
    if "tp_loss_rate" in plotted.columns:
        ax2.bar(x, plotted["tp_loss_rate"], color="#8E44AD", edgecolor="black")
        ax2.axhline(0.05, ls="--", c="red", lw=0.7, label="5% TP loss threshold")
        ax2.set_xticks(x)
        ax2.set_xticklabels(plotted["stratum"], rotation=45, ha="right", fontsize=7)
        ax2.set_ylabel("TP loss rate (TP flagged / total TP)")
        ax2.set_title("SuggestFilter TP loss rate per stratum")
        for xi, (_, row) in zip(x, plotted.iterrows()):
            if pd.notna(row["tp_loss_rate"]):
                ax2.text(xi, row["tp_loss_rate"] + 0.01, f"{row['tp_loss_rate']:.2f}",
                         ha="center", fontsize=6)
        ax2.legend(fontsize=8)
    fig.suptitle("G10 Step 8 · SuggestFilter precision/recall", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_FIG / "fig08_suggestfilter_prec_rec.png", dpi=130)
    plt.close(fig)
    return out


# ----------------------------------------------------------------------------
# Step 9: Quality_Score vs NumReads baseline gain

def step9_qs_vs_numreads(df: pd.DataFrame) -> pd.DataFrame:
    if "Quality_Score" not in df.columns:
        return None
    rows = []
    for mode in ["paired_full", "to_pileup"]:
        for s in SAMPLE_ORDER:
            sub = df[(df["sample"] == s) & (df["mode"] == mode)] \
                  .dropna(subset=["Quality_Score", "NumReads", "tp_label"])
            if len(sub) < 100:
                continue
            y = sub["tp_label"].to_numpy(int)
            auc_qs = auc_wilson(y, sub["Quality_Score"].to_numpy(float))[0]
            auc_nr = auc_wilson(y, sub["NumReads"].to_numpy(float))[0]
            rows.append(dict(sample=s, mode=mode, n=len(sub),
                             auc_QS=auc_qs, auc_NumReads=auc_nr,
                             delta=auc_qs - auc_nr if pd.notna(auc_qs) and pd.notna(auc_nr) else np.nan))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G10_qs_vs_numreads.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(12, 6))
    if len(out):
        labels = [f"{r.sample}/{r['mode'].split('_')[0]}" for _, r in out.iterrows()]
        x = np.arange(len(out))
        ax.bar(x - 0.2, out["auc_QS"], width=0.38, color="#2E86AB",
               edgecolor="black", label="Quality_Score AUC")
        ax.bar(x + 0.2, out["auc_NumReads"], width=0.38, color="#F4A261",
               edgecolor="black", label="NumReads AUC")
        for xi, (_, row) in zip(x, out.iterrows()):
            if pd.notna(row["auc_QS"]):
                ax.text(xi - 0.2, row["auc_QS"] + 0.005, f"{row['auc_QS']:.2f}",
                        ha="center", fontsize=6)
            if pd.notna(row["auc_NumReads"]):
                ax.text(xi + 0.2, row["auc_NumReads"] + 0.005, f"{row['auc_NumReads']:.2f}",
                        ha="center", fontsize=6)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.axhline(0.5, ls="--", c="gray", lw=0.7)
        ax.axhline(0.58, ls=":", c="red", lw=0.7)
        ax.set_ylabel("AUC"); ax.set_ylim(0.3, 1.0)
        ax.set_title("Quality_Score vs NumReads baseline AUC per (sample, mode)")
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig09_qs_vs_numreads.png", dpi=130)
    plt.close(fig)
    return out


# ----------------------------------------------------------------------------
# Main

def main():
    df = load_master()
    print(f"[main] loaded master rows={len(df):,}", flush=True)

    print("[main] Step 1 global distribution", flush=True)
    step1_out = step1_global(df)
    print(step1_out.head(15).to_string(index=False))

    # Select flagship features (likely to matter)
    focus = [f for f in ["Quality_Score", "HeuristicScore", "Stability",
                         "LabelAllelePermanovaF", "LabelHPPermanovaF",
                         "UnassignedAffinity", "NHP3", "NHP0"]
             if f in df.columns]

    for feat in focus[:4]:
        print(f"[main] Step 2 heatmap {feat}", flush=True)
        step2_heatmap(df, feat)

    print("[main] Step 3 per-sample", flush=True)
    step3_per_sample(df)

    print("[main] Step 4 stratified AUC", flush=True)
    step4_stratified(df, focus)

    print("[main] Step 5 confound (residualize on AF/cov)", flush=True)
    step5_confound(df, ["Quality_Score", "HeuristicScore",
                        "LabelAllelePermanovaF", "LabelHPPermanovaF",
                        "UnassignedAffinity"])

    spatial_rows = []
    for feat in ["Quality_Score", "LabelAllelePermanovaF"]:
        print(f"[main] Step 6 spatial {feat}", flush=True)
        sp = step6_spatial(df, feat)
        if sp is not None:
            spatial_rows.append(sp)
    if spatial_rows:
        pd.DataFrame(spatial_rows).to_csv(OUT_DATA / "G10_spatial.tsv", sep="\t", index=False)

    print("[main] Step 7 VerificationClass", flush=True)
    step7_verification_class(df)

    print("[main] Step 8 SuggestFilter prec/rec", flush=True)
    step8_suggest_filter(df)

    print("[main] Step 9 QS vs NumReads", flush=True)
    step9_qs_vs_numreads(df)

    print("[main] done.")
    print(f"  figures: {OUT_FIG}")
    print(f"  tables : {OUT_DATA}")


if __name__ == "__main__":
    main()
