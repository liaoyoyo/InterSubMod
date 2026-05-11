#!/usr/bin/env python3
"""Phase C · G7 — Cluster & Global Methylation Structure observation.

Scope:
  11 features (G7 registry):
    GlobalP, CramersV,
    PairwiseMeanDist, PairwiseMedianDist,
    ClusterPermanovaF, ClusterPermanovaP, ClusterPermanovaValid,
    ClusterDispersionP, ClusterDispersionWarn,
    LocalBestCluster, LocalBestP

Data:
  - Master: research/feature_layered_observation/data/merged_with_vcf.tsv.gz
           (LOH/AF/CN axes + tp_label)
  - Raw G7 columns enriched from per-(sample, mode) significance_summary.csv
    (same provenance as G6 enrichment).

Key derivations:
  - neg_log10_GlobalP      = -log10(max(GlobalP, 1e-300))
  - neg_log10_ClusterPermP = -log10(max(ClusterPermanovaP, 1e-300))
  - neg_log10_LocalBestP   = -log10(max(LocalBestP, 1e-300))
  - ClusterPermanovaF is scored only on ClusterPermanovaValid==True subset.

Methodology: research/feature_layered_observation/02_methodology.md Step 1-6.
  Step 1 global dist (violin+KDE)        -> fig01*
  Step 2 LOH×AF×CN 32-cell heatmap        -> fig02*
  Step 3 per-sample consistency (7 grid)  -> fig03*
  Step 4 stratified AUC                   -> fig04*
  Step 5 confound guard (OLS on NumReads+vcf_AF+Coverage_Multiple) -> fig05*
  Step 6 spatial autocorrelation 5Mb bins -> fig06*

Outputs:
  figures/G7_cluster/*.png
  data/G7_auc_table.tsv, G7_cell_delta.tsv, G7_global_stats.tsv,
  G7_confound.tsv, G7_pairwise_collinearity.tsv
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
from sklearn.linear_model import LinearRegression

warnings.filterwarnings("ignore")

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER_TSV = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
FIG_DIR = ROOT / "research/feature_layered_observation/figures/G7_cluster"
DATA_DIR = ROOT / "research/feature_layered_observation/data"
OUT_SUB = DATA_DIR / "G7"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_SUB.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]

# Master uses "LOH_None/Weak/Noise/Strong/Subclone" prefixes (see G1 script).
LOH_ORDER = ["LOH_None", "LOH_Weak", "LOH_Noise", "LOH_Strong", "LOH_Subclone"]
AF_ORDER = ["Extreme", "Near-half", "Intermediate"]
CN_BREAKS = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
CN_LABELS = ["CN_Loss", "CN_Near1", "CN_Diploid", "CN_Gain", "CN_HighGain"]

# Raw enrichment source paths (mirrors G6_hp_fine.py) ---------------------------
NEW_KDE_PAIRED_FULL = {
    "HCC1395": ROOT / "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO": ROOT / "output/canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full",
    "H2009": ROOT / "output/canonical/H2009/paired_full/20260421_H2009_paired_full_full",
    "H1437": ROOT / "output/canonical/H1437/paired_full/20260421_H1437_paired_full_full",
    "COLO829": ROOT / "output/canonical/COLO829/paired_full/20260421_COLO829_paired_full_full",
    "HCC1937": ROOT / "output/canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full",
    "HCC1954": ROOT / "output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix",
}
HCC1395_TO = Path("/tmp/ism_hp_fix_phase1")
NEW_KDE_TO_ARCHIVE = {
    "HCC1395_DORADO": Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod"),
    "H1437": Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod"),
    "H2009": Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod"),
    "HCC1937": Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod"),
    "HCC1954": Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod"),
}
COLO829_TO = ROOT / "output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod"

G7_RAW_COLS = [
    "GlobalP", "CramersV",
    "PairwiseMeanDist", "PairwiseMedianDist",
    "ClusterPermanovaF", "ClusterPermanovaP", "ClusterPermanovaValid",
    "ClusterDispersionP", "ClusterDispersionWarn",
    "LocalBestCluster", "LocalBestP",
]

# Features used in Step 1-6 plots (after derivations)
G7_CONTINUOUS = [
    "GlobalP", "neg_log10_GlobalP", "CramersV",
    "PairwiseMeanDist", "PairwiseMedianDist",
    "ClusterPermanovaF",                # filtered to Valid==True
    "ClusterPermanovaP", "neg_log10_ClusterPermP",
    "ClusterDispersionP",
    "LocalBestP", "neg_log10_LocalBestP",
]
G7_BINARY = ["ClusterPermanovaValid", "ClusterDispersionWarn"]
G7_CATEGORICAL = ["LocalBestCluster"]


# ---------------------------------------------------------------------------- #
# I/O                                                                          #
# ---------------------------------------------------------------------------- #
def load_master() -> pd.DataFrame:
    print(f"[load] master {MASTER_TSV}", flush=True)
    usecols = ["sample", "mode", "Chr", "Pos", "RegionID", "tp_label",
               "LOH_Subtype", "AF_class", "vcf_AF", "AF",
               "Coverage_Multiple", "NumReads", "NumCpGs"]
    df = pd.read_csv(MASTER_TSV, sep="\t", usecols=usecols, low_memory=False)
    df["cn_tier_F"] = pd.cut(df["Coverage_Multiple"].astype(float),
                             bins=CN_BREAKS, labels=CN_LABELS, right=True)
    df["LOH_Subtype"] = pd.Categorical(df["LOH_Subtype"].fillna("LOH_None"),
                                       categories=LOH_ORDER, ordered=True)
    df["AF_class"] = pd.Categorical(df["AF_class"].fillna("Unknown"),
                                    categories=AF_ORDER + ["Unknown"], ordered=False)
    df["tp_label"] = df["tp_label"].astype(int)
    print(f"[load] master rows={len(df):,}", flush=True)
    return df


def load_raw_run(run_dir: Path) -> pd.DataFrame:
    tp = pd.read_csv(run_dir / "intersubmod_tp" / "significance_summary.csv", low_memory=False)
    fp = pd.read_csv(run_dir / "intersubmod_fp" / "significance_summary.csv", low_memory=False)
    keep = ["RegionID", "Chr", "Pos"] + [c for c in G7_RAW_COLS if c in tp.columns]
    tp = tp[keep].copy(); tp["_tp"] = 1
    fp = fp[keep].copy(); fp["_tp"] = 0
    return pd.concat([tp, fp], ignore_index=True)


def enrich(master: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for (sample, mode), _grp in master.groupby(["sample", "mode"], sort=False):
        run_dir = None
        if mode == "paired_full" and sample in NEW_KDE_PAIRED_FULL:
            run_dir = NEW_KDE_PAIRED_FULL[sample]
        elif mode == "to_pileup":
            if sample == "HCC1395":
                tp_path = HCC1395_TO / "tp_off" / "significance_summary.csv"
                fp_path = HCC1395_TO / "fp_off" / "significance_summary.csv"
                if tp_path.exists() and fp_path.exists():
                    tp = pd.read_csv(tp_path, low_memory=False); tp["_tp"] = 1
                    fp = pd.read_csv(fp_path, low_memory=False); fp["_tp"] = 0
                    keep = ["RegionID", "Chr", "Pos"] + [c for c in G7_RAW_COLS if c in tp.columns] + ["_tp"]
                    df_raw = pd.concat([tp[keep], fp[keep]], ignore_index=True)
                    df_raw["sample"] = sample; df_raw["mode"] = mode
                    frames.append(df_raw); continue
            elif sample == "COLO829":
                run_dir = COLO829_TO
            elif sample in NEW_KDE_TO_ARCHIVE:
                run_dir = NEW_KDE_TO_ARCHIVE[sample]
        if run_dir is None:
            print(f"[warn] no run_dir for ({sample},{mode})"); continue
        try:
            df_raw = load_raw_run(Path(run_dir))
        except FileNotFoundError as e:
            print(f"[warn] missing raw for ({sample},{mode}): {e}"); continue
        df_raw["sample"] = sample; df_raw["mode"] = mode
        frames.append(df_raw)
    raw = pd.concat(frames, ignore_index=True)
    print(f"[enrich] raw rows: {len(raw):,}", flush=True)
    merged = master.merge(
        raw, on=["sample", "mode", "RegionID", "Chr", "Pos"],
        how="left", suffixes=("", "_r"),
    )
    mismatch = (merged["tp_label"] != merged["_tp"]).sum()
    print(f"[enrich] merge done; tp_label raw vs master mismatch: {mismatch:,}", flush=True)
    merged.drop(columns=["_tp"], inplace=True)
    return merged


def derive(df: pd.DataFrame) -> pd.DataFrame:
    """Create -log10 transforms and validity-filtered versions."""
    eps = 1e-300
    for p_col, out in [
        ("GlobalP", "neg_log10_GlobalP"),
        ("ClusterPermanovaP", "neg_log10_ClusterPermP"),
        ("LocalBestP", "neg_log10_LocalBestP"),
    ]:
        if p_col in df.columns:
            df[out] = -np.log10(pd.to_numeric(df[p_col], errors="coerce").clip(lower=eps))
    # Boolean flags to {0,1}
    for b in ["ClusterPermanovaValid", "ClusterDispersionWarn"]:
        if b in df.columns:
            v = df[b]
            if v.dtype == object:
                df[b] = v.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(int)
            else:
                df[b] = pd.to_numeric(v, errors="coerce").fillna(0).astype(int)
    return df


# ---------------------------------------------------------------------------- #
# Stats helpers                                                                #
# ---------------------------------------------------------------------------- #
def auc_wilson_ci(y_true: np.ndarray, scores: np.ndarray):
    mask = ~(np.isnan(scores) | np.isnan(y_true.astype(float)))
    y = y_true[mask]; s = scores[mask]
    n_pos = int((y == 1).sum()); n_neg = int((y == 0).sum())
    if n_pos < 5 or n_neg < 5:
        return np.nan, np.nan, np.nan, n_pos, n_neg
    rnks = stats.rankdata(s)
    sum_rank_pos = rnks[y == 1].sum()
    auc = (sum_rank_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    q1 = auc / (2 - auc); q2 = 2 * auc**2 / (1 + auc)
    se = np.sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc**2)
                  + (n_neg - 1) * (q2 - auc**2)) / (n_pos * n_neg))
    return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se), n_pos, n_neg


def cohen_d(x, y):
    x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
    if len(x) < 2 or len(y) < 2: return np.nan
    nx, ny = len(x), len(y)
    sx, sy = x.std(ddof=1), y.std(ddof=1)
    sp = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    return 0.0 if sp == 0 else (x.mean() - y.mean()) / sp


def mwu_p(x, y):
    x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
    if len(x) < 5 or len(y) < 5: return np.nan
    try:
        return float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except Exception:
        return np.nan


def scores_for_feat(df: pd.DataFrame, feat: str) -> np.ndarray:
    """Return the AUC score vector; filters ClusterPermanovaF to Valid==True."""
    if feat == "ClusterPermanovaF":
        s = pd.to_numeric(df["ClusterPermanovaF"], errors="coerce").copy()
        if "ClusterPermanovaValid" in df.columns:
            s = s.where(df["ClusterPermanovaValid"] == 1, np.nan)
        return s.to_numpy()
    if feat in df.columns:
        return pd.to_numeric(df[feat], errors="coerce").to_numpy()
    return np.full(len(df), np.nan)


# ---------------------------------------------------------------------------- #
# Step 1 · global distribution                                                 #
# ---------------------------------------------------------------------------- #
def step1_global(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # Continuous
    for feat in G7_CONTINUOUS:
        s = scores_for_feat(df, feat)
        y = df["tp_label"].to_numpy()
        auc, lo, hi, n_pos, n_neg = auc_wilson_ci(y, s)
        tp = s[y == 1]; fp = s[y == 0]
        rows.append(dict(feature=feat, kind="cont",
                         auc=auc, lo=lo, hi=hi,
                         cohens_d=cohen_d(tp, fp),
                         mwu_p=mwu_p(tp, fp),
                         mean_tp=float(np.nanmean(tp)) if np.isfinite(tp).any() else np.nan,
                         mean_fp=float(np.nanmean(fp)) if np.isfinite(fp).any() else np.nan,
                         median_tp=float(np.nanmedian(tp)) if np.isfinite(tp).any() else np.nan,
                         median_fp=float(np.nanmedian(fp)) if np.isfinite(fp).any() else np.nan,
                         n_tp_valid=int(np.isfinite(tp).sum()),
                         n_fp_valid=int(np.isfinite(fp).sum())))
    # Binary — TP rate by flag level as score
    for feat in G7_BINARY:
        if feat not in df.columns: continue
        y = df["tp_label"].to_numpy()
        s = df[feat].astype(float).to_numpy()
        auc, lo, hi, n_pos, n_neg = auc_wilson_ci(y, s)
        rate_0 = df.loc[df[feat] == 0, "tp_label"].mean() if (df[feat] == 0).any() else np.nan
        rate_1 = df.loc[df[feat] == 1, "tp_label"].mean() if (df[feat] == 1).any() else np.nan
        rows.append(dict(feature=feat, kind="binary",
                         auc=auc, lo=lo, hi=hi,
                         cohens_d=np.nan, mwu_p=np.nan,
                         mean_tp=rate_1, mean_fp=rate_0,
                         median_tp=rate_1, median_fp=rate_0,
                         n_tp_valid=int((df[feat] == 1).sum()),
                         n_fp_valid=int((df[feat] == 0).sum())))
    # Categorical — map category -> TP rate, then AUC
    for feat in G7_CATEGORICAL:
        if feat not in df.columns: continue
        y = df["tp_label"].to_numpy()
        cat = df[feat].astype("category")
        rate = df.groupby(cat, observed=True)["tp_label"].mean()
        scores = df[feat].map(rate).to_numpy(dtype=float)
        auc, lo, hi, *_ = auc_wilson_ci(y, scores)
        rows.append(dict(feature=feat, kind="cat",
                         auc=auc, lo=lo, hi=hi,
                         cohens_d=np.nan, mwu_p=np.nan,
                         mean_tp=np.nan, mean_fp=np.nan,
                         median_tp=np.nan, median_fp=np.nan,
                         n_tp_valid=int((y == 1).sum()),
                         n_fp_valid=int((y == 0).sum())))
    out = pd.DataFrame(rows).sort_values("auc", ascending=False)
    out.to_csv(DATA_DIR / "G7_global_stats.tsv", sep="\t", index=False)
    return out


def fig01_auc_bar(stats_df: pd.DataFrame):
    d = stats_df.dropna(subset=["auc"]).copy().sort_values("auc")
    fig, ax = plt.subplots(figsize=(11, 7), dpi=130)
    yp = np.arange(len(d))
    xerr = np.vstack([(d["auc"] - d["lo"]).clip(lower=0).values,
                      (d["hi"] - d["auc"]).clip(lower=0).values])
    colors = ["#A85540" if v >= 0.58 else ("#5B8DB7" if v >= 0.53 else "#8A8A8A") for v in d["auc"]]
    ax.barh(yp, d["auc"], xerr=xerr, color=colors, ecolor="#333", capsize=3)
    for i, row in enumerate(d.itertuples()):
        ax.text(row.auc + 0.005, i, f"{row.auc:.3f}", va="center", fontsize=8)
    ax.set_yticks(yp); ax.set_yticklabels(d["feature"], fontsize=9)
    ax.axvline(0.5, ls="--", c="#888", lw=1)
    ax.axvline(0.58, ls="--", c="#c00", lw=1, label="Beyond-AUC 0.58")
    ax.set_xlabel("AUC (TP vs FP, Hanley-McNeil 95% CI)")
    ax.set_xlim(0.35, 0.85)
    ax.set_title("G7 · Cluster & Global Methyl · Global AUC ranking (7 samples × 2 modes pooled)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG_DIR / "fig01_global_auc.png"); plt.close(fig)


def fig01_violin(df: pd.DataFrame, feats: list[str]):
    """Grid of raw vs -log10 distributions for continuous features."""
    n = len(feats)
    ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.6 * nrow), dpi=120)
    axes = axes.flatten() if nrow * ncol > 1 else [axes]
    y = df["tp_label"].to_numpy()
    for i, feat in enumerate(feats):
        ax = axes[i]
        s = scores_for_feat(df, feat)
        tp = s[(y == 1) & np.isfinite(s)]
        fp = s[(y == 0) & np.isfinite(s)]
        if len(tp) < 10 or len(fp) < 10:
            ax.text(0.5, 0.5, f"{feat}\n(insufficient)", ha="center"); ax.axis("off"); continue
        # subsample for speed
        tp_p = np.random.default_rng(0).choice(tp, size=min(20000, len(tp)), replace=False)
        fp_p = np.random.default_rng(0).choice(fp, size=min(20000, len(fp)), replace=False)
        parts = ax.violinplot([fp_p, tp_p], positions=[0, 1], widths=0.7,
                              showmeans=False, showmedians=True)
        for pc, c in zip(parts["bodies"], ["#D55E00", "#5B8DB7"]):
            pc.set_facecolor(c); pc.set_alpha(0.6)
        auc, lo, hi, *_ = auc_wilson_ci(y, s)
        d = cohen_d(tp, fp)
        ax.set_xticks([0, 1]); ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"], fontsize=8)
        ax.set_title(f"{feat}\nAUC={auc:.3f} [{lo:.3f},{hi:.3f}] d={d:.2f}", fontsize=9)
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle("G7 · Step 1 Global TP vs FP distributions", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / "fig01_violin_grid.png"); plt.close(fig)


def fig01_pvalue_comparison(df: pd.DataFrame):
    """Explicit P-value vs -log10 comparison (answers question 2)."""
    pairs = [("GlobalP", "neg_log10_GlobalP"),
             ("ClusterPermanovaP", "neg_log10_ClusterPermP"),
             ("LocalBestP", "neg_log10_LocalBestP")]
    y = df["tp_label"].to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=120)
    for ax, (raw, log) in zip(axes, pairs):
        s_raw = scores_for_feat(df, raw)
        s_log = scores_for_feat(df, log)
        auc_raw, *_ = auc_wilson_ci(y, s_raw)
        # AUC is rank-invariant to monotone transform -> should match (up to NaN).
        auc_log, *_ = auc_wilson_ci(y, s_log)
        # Plot -log10 histogram by tp
        tp = s_log[(y == 1) & np.isfinite(s_log)]
        fp = s_log[(y == 0) & np.isfinite(s_log)]
        if len(tp) > 0 and len(fp) > 0:
            lim = float(np.nanpercentile(np.concatenate([tp, fp]), 99))
            ax.hist(fp.clip(max=lim), bins=60, density=True, alpha=0.5,
                    color="#D55E00", label=f"FP n={len(fp):,}")
            ax.hist(tp.clip(max=lim), bins=60, density=True, alpha=0.5,
                    color="#5B8DB7", label=f"TP n={len(tp):,}")
        ax.set_xlabel(f"-log10({raw})"); ax.set_ylabel("density")
        ax.set_title(f"raw AUC={auc_raw:.4f}   log AUC={auc_log:.4f}\n"
                     f"(rank-invariant ⇒ identical up to NaN)", fontsize=9)
        ax.legend(fontsize=8)
    fig.suptitle("G7 · Step 1b  P-value vs -log10(P) AUC comparison", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "fig01b_pvalue_vs_log.png"); plt.close(fig)


# ---------------------------------------------------------------------------- #
# Step 2 · LOH × AF × CN 32-cell heatmap                                       #
# ---------------------------------------------------------------------------- #
def step2_heatmap(df: pd.DataFrame, feat: str, out_name: str, min_n: int = 20):
    d = df.copy()
    d["row_key"] = d["LOH_Subtype"].astype(str) + " | " + d["AF_class"].astype(str)
    row_order = [f"{l} | {a}" for l in LOH_ORDER for a in AF_ORDER]
    col_order = CN_LABELS
    piv_n = d.pivot_table(index="row_key", columns="cn_tier_F",
                          values="tp_label", aggfunc="size").reindex(row_order, axis=0) \
             .reindex(col_order, axis=1)
    piv_rate = d.pivot_table(index="row_key", columns="cn_tier_F",
                             values="tp_label", aggfunc="mean").reindex(row_order, axis=0) \
                .reindex(col_order, axis=1)
    mask_ok = piv_n >= min_n

    s = scores_for_feat(d, feat)
    d["_feat"] = s
    piv_tp = d[d["tp_label"] == 1].pivot_table(
        index="row_key", columns="cn_tier_F", values="_feat", aggfunc="median"
    ).reindex(row_order, axis=0).reindex(col_order, axis=1)
    piv_fp = d[d["tp_label"] == 0].pivot_table(
        index="row_key", columns="cn_tier_F", values="_feat", aggfunc="median"
    ).reindex(row_order, axis=0).reindex(col_order, axis=1)
    piv_delta = (piv_tp - piv_fp).where(mask_ok)

    fig, axes = plt.subplots(1, 4, figsize=(24, 11), dpi=120)
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
        ax.set_yticks(range(len(row_order))); ax.set_yticklabels(row_order, fontsize=7)
        ax.set_xticks(range(len(col_order))); ax.set_xticklabels(col_order, rotation=30, ha="right", fontsize=8)
    for i, r in enumerate(row_order):
        for j, c in enumerate(col_order):
            n = piv_n.loc[r, c] if (r in piv_n.index and c in piv_n.columns) else 0
            if pd.isna(n): n = 0
            if n >= min_n:
                axes[0].text(j, i, f"{int(n):,}", ha="center", va="center", fontsize=6)
    fig.suptitle(f"{feat} · Step 2 LOH × AF × CN heatmap (7 samples × 2 modes pooled)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / f"fig02_{out_name}.png"); plt.close(fig)

    flat = piv_delta.reset_index().melt(id_vars="row_key", var_name="cn_tier",
                                        value_name="delta_tp_fp")
    flat["n"] = flat.apply(lambda r: piv_n.loc[r["row_key"], r["cn_tier"]]
                           if r["row_key"] in piv_n.index and r["cn_tier"] in piv_n.columns else 0,
                           axis=1)
    flat["feature"] = feat
    return flat


# ---------------------------------------------------------------------------- #
# Step 3 · per-sample consistency                                              #
# ---------------------------------------------------------------------------- #
def step3_per_sample(df: pd.DataFrame, feat: str, out_name: str, min_n: int = 10):
    samples = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    fig, axes = plt.subplots(2, 4, figsize=(22, 10), dpi=120)
    axes = axes.flatten()
    per_sample_rate = {}
    for i, sm in enumerate(samples):
        sub = df[df["sample"] == sm]
        piv = sub.pivot_table(index="LOH_Subtype", columns="AF_class",
                              values="tp_label", aggfunc="mean") \
                 .reindex(LOH_ORDER, axis=0).reindex(AF_ORDER, axis=1)
        piv_n = sub.pivot_table(index="LOH_Subtype", columns="AF_class",
                                values="tp_label", aggfunc="size") \
                   .reindex(LOH_ORDER, axis=0).reindex(AF_ORDER, axis=1)
        mask = piv_n >= min_n
        im = axes[i].imshow(piv.where(mask).values, aspect="auto",
                            cmap="RdBu_r", vmin=0, vmax=1)
        axes[i].set_title(f"{sm}\n(n={len(sub):,})", fontsize=9)
        axes[i].set_xticks(range(len(AF_ORDER))); axes[i].set_xticklabels(AF_ORDER, rotation=30, ha="right", fontsize=7)
        axes[i].set_yticks(range(len(LOH_ORDER))); axes[i].set_yticklabels(LOH_ORDER, fontsize=7)
        for ii, r in enumerate(LOH_ORDER):
            for jj, c in enumerate(AF_ORDER):
                n = piv_n.loc[r, c] if r in piv_n.index and c in piv_n.columns else 0
                if pd.isna(n): n = 0
                if n >= min_n:
                    v = piv.loc[r, c]
                    axes[i].text(jj, ii, f"{v:.2f}\nn={int(n)}", ha="center", va="center", fontsize=6)
        per_sample_rate[sm] = piv.values.flatten()
    # Spearman concordance on per-cell AUC for the feature itself
    ax_conc = axes[-1]
    auc_per_cell = {}
    for sm in samples:
        sub = df[df["sample"] == sm]
        s_feat = scores_for_feat(sub, feat)
        y = sub["tp_label"].to_numpy()
        # Per (LOH, AF) cell AUC
        vec = []
        for l in LOH_ORDER:
            for a in AF_ORDER:
                m = (sub["LOH_Subtype"] == l) & (sub["AF_class"] == a)
                if m.sum() >= 50:
                    au, *_ = auc_wilson_ci(y[m], s_feat[m.values])
                    vec.append(au)
                else:
                    vec.append(np.nan)
        auc_per_cell[sm] = np.array(vec)
    nn = len(samples)
    mat = np.full((nn, nn), np.nan)
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            v1 = auc_per_cell[s1]; v2 = auc_per_cell[s2]
            m = ~(np.isnan(v1) | np.isnan(v2))
            if m.sum() >= 3:
                mat[i, j] = stats.spearmanr(v1[m], v2[m]).correlation
    im2 = ax_conc.imshow(mat, cmap="coolwarm", vmin=-1, vmax=1)
    ax_conc.set_title(f"Cross-sample Spearman\n(per-cell AUC of {feat})", fontsize=9)
    ax_conc.set_xticks(range(nn)); ax_conc.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax_conc.set_yticks(range(nn)); ax_conc.set_yticklabels(samples, fontsize=7)
    fig.colorbar(im2, ax=ax_conc, fraction=0.04)
    fig.suptitle(f"{feat} · Step 3 Per-sample TP rate consistency (LOH × AF)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / f"fig03_{out_name}_per_sample.png"); plt.close(fig)


# ---------------------------------------------------------------------------- #
# Step 4 · stratified AUC                                                      #
# ---------------------------------------------------------------------------- #
def step4_stratified(df: pd.DataFrame, feat: str, out_name: str):
    y = df["tp_label"].to_numpy()
    s = scores_for_feat(df, feat)
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
    for mode in df["mode"].unique():
        m = (df["mode"] == mode).to_numpy()
        if m.sum() >= 100:
            rows.append(("mode", mode, *auc_wilson_ci(y[m], s[m])))
    dfres = pd.DataFrame(rows, columns=["layer", "group", "auc", "lo", "hi", "n_pos", "n_neg"])
    dfres["feature"] = feat

    fig, ax = plt.subplots(figsize=(12, 5.5), dpi=120)
    xpos = np.arange(len(dfres))
    colors_m = {"global": "#333", "LOH": "#2E86AB", "AF": "#F4A261",
                "CN": "#6A994E", "mode": "#9467BD"}
    for i, (_, row) in enumerate(dfres.iterrows()):
        if np.isnan(row["auc"]): continue
        yerr = [[max(0, row["auc"] - row["lo"])], [max(0, row["hi"] - row["auc"])]]
        ax.bar(i, row["auc"], color=colors_m.get(row["layer"], "gray"),
               yerr=yerr, capsize=3, edgecolor="black", linewidth=0.4)
        ax.text(i, row["auc"] + 0.01,
                f"{row['auc']:.3f}\nn={int(row['n_pos'])+int(row['n_neg']):,}",
                ha="center", fontsize=6)
    ax.axhline(0.5, ls="--", c="gray", lw=0.8, label="random")
    ax.axhline(0.58, ls=":", c="red", lw=0.8, label="Beyond-AUC 0.58")
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"{r.layer}:{r.group}" for _, r in dfres.iterrows()],
                       rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0.3, 1.0); ax.set_ylabel("AUC")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"{feat} · Step 4 Stratified AUC (Global / LOH / AF / CN / Mode)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"fig04_{out_name}_stratified_auc.png"); plt.close(fig)
    return dfres


# ---------------------------------------------------------------------------- #
# Step 5 · confound guard                                                      #
# ---------------------------------------------------------------------------- #
def step5_confound(df: pd.DataFrame, feat: str, out_name: str):
    covars = ["NumReads", "vcf_AF", "Coverage_Multiple"]
    sub = df.copy()
    sub["_feat"] = scores_for_feat(sub, feat)
    sub = sub.dropna(subset=["_feat"] + covars + ["tp_label"]).copy()
    if len(sub) < 500:
        return None
    X = sub[covars].astype(float).to_numpy()
    yfeat = sub["_feat"].astype(float).to_numpy()
    ytp = sub["tp_label"].to_numpy()
    try:
        lr = LinearRegression().fit(X, yfeat)
        resid = yfeat - lr.predict(X)
    except Exception:
        return None
    raw_auc, *_ = auc_wilson_ci(ytp, yfeat)
    res_auc, *_ = auc_wilson_ci(ytp, resid)

    # AF-bin cross
    af_rows = []
    for af in AF_ORDER:
        m = (sub["AF_class"] == af).to_numpy()
        if m.sum() >= 300:
            r_auc, *_ = auc_wilson_ci(ytp[m], yfeat[m])
            r_auc_res, *_ = auc_wilson_ci(ytp[m], resid[m])
            af_rows.append((af, r_auc, r_auc_res, int(m.sum())))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120)
    bins = 60
    for ax, vec, title, auc in [
        (axes[0], yfeat, f"raw {feat}  AUC={raw_auc:.3f}", raw_auc),
        (axes[1], resid, f"resid. on {covars}\nAUC={res_auc:.3f}", res_auc),
    ]:
        tp_v = vec[ytp == 1]; fp_v = vec[ytp == 0]
        if len(tp_v) > 0 and len(fp_v) > 0:
            lim_lo = float(np.nanpercentile(np.concatenate([tp_v, fp_v]), 0.5))
            lim_hi = float(np.nanpercentile(np.concatenate([tp_v, fp_v]), 99.5))
            ax.hist(np.clip(fp_v, lim_lo, lim_hi), bins=bins, alpha=0.45,
                    density=True, color="#D55E00", label="FP")
            ax.hist(np.clip(tp_v, lim_lo, lim_hi), bins=bins, alpha=0.45,
                    density=True, color="#5B8DB7", label="TP")
        ax.set_title(title, fontsize=9); ax.legend(fontsize=8)
    af_txt = "  ".join([f"{a}: raw={r:.3f} res={rr:.3f} (n={n:,})"
                        for a, r, rr, n in af_rows])
    fig.suptitle(f"{feat} · Step 5 Confound guard (within-cell OLS residualization)\n"
                 f"AF-bin AUC — {af_txt}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(FIG_DIR / f"fig05_{out_name}_confound.png"); plt.close(fig)
    return dict(feature=feat, raw_auc=raw_auc, resid_auc=res_auc,
                n=len(sub),
                af_bin=";".join([f"{a}:raw={r:.3f},res={rr:.3f},n={n}"
                                 for a, r, rr, n in af_rows]))


# ---------------------------------------------------------------------------- #
# Step 6 · spatial autocorrelation                                             #
# ---------------------------------------------------------------------------- #
def step6_spatial(df: pd.DataFrame, feat: str, out_name: str, bin_mb: int = 5):
    d = df.copy()
    d["_feat"] = scores_for_feat(d, feat)
    d = d.dropna(subset=["Chr", "Pos", "_feat", "tp_label"])
    d["pos_bin"] = (d["Pos"].astype(int) // (bin_mb * 1_000_000))
    d["bin_key"] = d["Chr"].astype(str) + ":" + d["pos_bin"].astype(str)
    rows = []
    for key, g in d.groupby("bin_key"):
        if len(g) < 50: continue
        y = g["tp_label"].to_numpy()
        if y.sum() < 5 or (len(y) - y.sum()) < 5: continue
        a, *_ = auc_wilson_ci(y, g["_feat"].astype(float).to_numpy())
        rows.append({"bin": key, "n": len(g), "tp_rate": y.mean(), "auc": a})
    if not rows: return None
    s = pd.DataFrame(rows).dropna(subset=["auc"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120)
    ax = axes[0]
    ax.scatter(s["tp_rate"], s["auc"], alpha=0.3, s=8)
    ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
    ax.axvline(0.8, ls="--", c="purple", alpha=0.5)
    ax.set_xlabel("bin TP rate"); ax.set_ylabel(f"bin AUC of {feat}")
    ax.set_title(f"{feat} · Step 6 spatial auto-correlation ({bin_mb} Mb bins)")
    axes[1].hist(s["auc"], bins=30, color="#5B8DB7", alpha=0.7)
    axes[1].axvline(0.5, ls="--", c="gray"); axes[1].axvline(0.58, ls=":", c="red")
    axes[1].set_xlabel("per-bin AUC"); axes[1].set_ylabel("# bins")
    hi = s[s["tp_rate"] > 0.8]; lo = s[s["tp_rate"] < 0.5]
    flag = ""
    if len(hi) >= 10 and len(lo) >= 10:
        delta = hi["auc"].mean() - lo["auc"].mean()
        flag = (f"⚠ artifact suspect (ΔAUC={delta:+.3f} hi vs lo)"
                if delta > 0.08 else f"ok (ΔAUC={delta:+.3f})")
    fig.suptitle(f"{feat} · spatial auto-correlation  {flag}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / f"fig06_{out_name}_spatial.png"); plt.close(fig)
    return dict(feature=feat, n_bins=len(s), flag=flag)


# ---------------------------------------------------------------------------- #
# Extra: CramersV zero-rate + PairwiseDist collinearity + HPFine overlap       #
# ---------------------------------------------------------------------------- #
def fig_cramers_zero_rate(df: pd.DataFrame, hpfine_df: pd.DataFrame | None):
    """Quantify 2x2 sparsity and overlap with HPFineNGroups / HPFine_NGroups_CF."""
    sub = df.dropna(subset=["CramersV"]).copy()
    n_total = len(sub)
    zero_mask = sub["CramersV"] <= 1e-9
    zero_rate = zero_mask.mean() if n_total else np.nan

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), dpi=120)
    # A: zero vs non-zero rates per sample
    ax = axes[0]
    grp = sub.groupby("sample").apply(lambda g: (g["CramersV"] <= 1e-9).mean())
    ax.bar(range(len(grp)), grp.values, color="#5B8DB7")
    ax.set_xticks(range(len(grp))); ax.set_xticklabels(grp.index, rotation=40, ha="right", fontsize=8)
    ax.axhline(0.93, ls="--", c="red", label="Prior: 93% zero ceiling")
    ax.set_ylabel("CramersV == 0 rate"); ax.set_ylim(0, 1.05)
    ax.set_title("A · CramersV zero-rate per sample"); ax.legend(fontsize=8)

    # B: TP rate zero vs nonzero
    ax = axes[1]
    tp_zero = sub.loc[zero_mask, "tp_label"].mean() if zero_mask.any() else np.nan
    tp_nz = sub.loc[~zero_mask, "tp_label"].mean() if (~zero_mask).any() else np.nan
    ax.bar([0, 1], [tp_zero, tp_nz], color=["#8A8A8A", "#A85540"])
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"CramersV=0\nn={int(zero_mask.sum()):,}",
                                                f"CramersV>0\nn={int((~zero_mask).sum()):,}"])
    ax.set_ylim(0, 1); ax.set_ylabel("TP rate")
    ax.axhline(sub["tp_label"].mean(), ls="--", c="black", label=f"global {sub['tp_label'].mean():.2f}")
    ax.set_title(f"B · TP rate by CramersV zero (zero_rate={zero_rate:.3f})")
    ax.legend(fontsize=8)

    # C: AUC comparison CramersV (raw) vs HPFine family if avail
    ax = axes[2]
    y = sub["tp_label"].to_numpy()
    auc_cv, *_ = auc_wilson_ci(y, sub["CramersV"].to_numpy())
    vals = [("CramersV", auc_cv)]
    for alt_col in ["CramersV_HPFine", "HPFineNGroups", "HPFine_NGroups_CF"]:
        if hpfine_df is not None and alt_col in hpfine_df.columns:
            merged = sub[["sample", "mode", "RegionID", "Chr", "Pos", "tp_label"]].merge(
                hpfine_df[["sample", "mode", "RegionID", "Chr", "Pos", alt_col]],
                on=["sample", "mode", "RegionID", "Chr", "Pos"], how="left"
            )
            auc_a, *_ = auc_wilson_ci(merged["tp_label"].to_numpy(),
                                      pd.to_numeric(merged[alt_col], errors="coerce").to_numpy())
            vals.append((alt_col, auc_a))
    xs = list(range(len(vals)))
    ax.bar(xs, [v[1] for v in vals],
           color=["#8A8A8A"] + ["#A85540"] * (len(vals) - 1))
    for i, (name, a) in enumerate(vals):
        ax.text(i, a + 0.01, f"{a:.3f}", ha="center", fontsize=8)
    ax.set_xticks(xs); ax.set_xticklabels([v[0] for v in vals], rotation=25, ha="right", fontsize=8)
    ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
    ax.set_ylim(0.35, 0.85); ax.set_ylabel("AUC")
    ax.set_title("C · CramersV vs HPFine family (replacement evidence)")
    fig.suptitle("G7 · CramersV zero-rate & HPFine replacement", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "fig07_cramersv_diagnostic.png"); plt.close(fig)
    return {"cramers_zero_rate": zero_rate,
            "tp_rate_at_zero": tp_zero,
            "tp_rate_at_nonzero": tp_nz,
            "auc_comparison": vals}


def fig_pairwise_collinearity(df: pd.DataFrame):
    """PairwiseMeanDist vs PairwiseMedianDist collinearity + vs NumReads."""
    sub = df.dropna(subset=["PairwiseMeanDist", "PairwiseMedianDist"]).copy()
    if len(sub) < 1000:
        return None
    # subsample for plotting
    smp = sub.sample(min(50_000, len(sub)), random_state=0)
    rho, _ = stats.spearmanr(sub["PairwiseMeanDist"], sub["PairwiseMedianDist"])
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=120)
    ax = axes[0]
    ax.scatter(smp["PairwiseMeanDist"], smp["PairwiseMedianDist"],
               c=smp["tp_label"].map({0: "#D55E00", 1: "#5B8DB7"}), alpha=0.3, s=4)
    ax.plot([0, 1], [0, 1], ls="--", c="black")
    ax.set_xlabel("PairwiseMeanDist"); ax.set_ylabel("PairwiseMedianDist")
    ax.set_title(f"A · Collinearity (Spearman ρ = {rho:.4f})")

    ax = axes[1]
    mean_num = stats.spearmanr(sub["PairwiseMeanDist"], sub["NumReads"])[0]
    med_num = stats.spearmanr(sub["PairwiseMedianDist"], sub["NumReads"])[0]
    ax.bar([0, 1], [abs(mean_num), abs(med_num)], color="#5B8DB7")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Mean vs NumReads", "Median vs NumReads"])
    ax.set_ylabel("|Spearman ρ|")
    ax.set_title(f"B · Pairwise dist vs NumReads\nmean={mean_num:+.3f}  median={med_num:+.3f}")

    ax = axes[2]
    # Per mode: AUC of each
    rows = []
    for mode in sub["mode"].unique():
        m = sub[sub["mode"] == mode]
        y = m["tp_label"].to_numpy()
        a_mn, *_ = auc_wilson_ci(y, m["PairwiseMeanDist"].to_numpy())
        a_md, *_ = auc_wilson_ci(y, m["PairwiseMedianDist"].to_numpy())
        rows.append((mode, a_mn, a_md))
    x = np.arange(len(rows)); w = 0.35
    ax.bar(x - w/2, [r[1] for r in rows], w, color="#2E86AB", label="Mean")
    ax.bar(x + w/2, [r[2] for r in rows], w, color="#F4A261", label="Median")
    ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
    for i, (_, a, b) in enumerate(rows):
        ax.text(i - w/2, a + 0.005, f"{a:.3f}", ha="center", fontsize=7)
        ax.text(i + w/2, b + 0.005, f"{b:.3f}", ha="center", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels([r[0] for r in rows]); ax.set_ylim(0.35, 0.8)
    ax.set_ylabel("AUC"); ax.set_title("C · AUC mean vs median per mode"); ax.legend(fontsize=8)
    fig.suptitle("G7 · Pairwise distance collinearity & AUC direction", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "fig08_pairwise_collinearity.png"); plt.close(fig)
    pd.DataFrame({"spearman_mean_median": [rho],
                  "mean_vs_numreads": [mean_num],
                  "median_vs_numreads": [med_num]}).to_csv(
                      DATA_DIR / "G7_pairwise_collinearity.tsv", sep="\t", index=False)
    return {"spearman": rho, "mean_vs_numreads": mean_num, "median_vs_numreads": med_num}


def fig_permanova_valid_gate(df: pd.DataFrame):
    """Effect of ClusterPermanovaValid gating on F-stat AUC."""
    sub = df.dropna(subset=["ClusterPermanovaF"]).copy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120)
    # A: valid rate per sample x mode
    ax = axes[0]
    piv = sub.groupby(["sample", "mode"])["ClusterPermanovaValid"].mean().unstack()
    piv = piv.reindex([s for s in SAMPLE_ORDER if s in piv.index])
    im = ax.imshow(piv.values, aspect="auto", cmap="Greens", vmin=0, vmax=1)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=25, ha="right")
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=9)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if v > 0.5 else "black", fontsize=8)
    plt.colorbar(im, ax=ax, fraction=0.04)
    ax.set_title("A · ClusterPermanovaValid rate per (sample × mode)")

    # B: AUC of F-stat with vs without Valid gate
    ax = axes[1]
    y_all = sub["tp_label"].to_numpy()
    auc_all, *_ = auc_wilson_ci(y_all, sub["ClusterPermanovaF"].to_numpy())
    mask_v = sub["ClusterPermanovaValid"] == 1
    auc_v, *_ = auc_wilson_ci(sub.loc[mask_v, "tp_label"].to_numpy(),
                              sub.loc[mask_v, "ClusterPermanovaF"].to_numpy())
    ax.bar([0, 1], [auc_all, auc_v], color=["#8A8A8A", "#A85540"])
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"all\nn={len(sub):,}",
                                                f"Valid=1\nn={int(mask_v.sum()):,}"])
    ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
    ax.set_ylabel("AUC of ClusterPermanovaF"); ax.set_ylim(0.35, 0.8)
    for i, v in enumerate([auc_all, auc_v]):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)
    ax.set_title("B · F-stat AUC raw vs Valid-gated")
    fig.suptitle("G7 · ClusterPermanovaValid gating impact", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "fig09_permanova_valid_gate.png"); plt.close(fig)
    return {"auc_all": auc_all, "auc_valid_only": auc_v,
            "delta": auc_v - auc_all if np.isfinite(auc_all) and np.isfinite(auc_v) else np.nan}


# ---------------------------------------------------------------------------- #
# Main driver                                                                  #
# ---------------------------------------------------------------------------- #
def _load_hpfine_for_overlap() -> pd.DataFrame | None:
    """Try to read G6 output for HPFine overlap, if it exists."""
    cand = ROOT / "research/feature_layered_observation/data/G6/step1_global_stats.tsv"
    if not cand.exists():
        return None
    # For AUC comparison we need per-variant; try enriched master cache
    return None  # heavier; skipped unless required


def main():
    master = load_master()
    df = enrich(master)
    df = derive(df)
    missing = [c for c in G7_RAW_COLS if c not in df.columns]
    if missing:
        print(f"[warn] missing after enrich: {missing}")

    # Step 1 global -----------------------------------------------------------
    print("[step1] global distributions", flush=True)
    stats_df = step1_global(df)
    fig01_auc_bar(stats_df)
    fig01_violin(df, [f for f in G7_CONTINUOUS if f in df.columns or f.startswith("neg_log10")])
    fig01_pvalue_comparison(df)

    # Step 2 heatmap for flagship features ------------------------------------
    flagship = [
        ("ClusterPermanovaF", "ClusterPermF"),
        ("neg_log10_GlobalP", "nlogGlobalP"),
        ("CramersV", "CramersV"),
        ("PairwiseMeanDist", "PairwiseMean"),
        ("PairwiseMedianDist", "PairwiseMedian"),
        ("neg_log10_LocalBestP", "nlogLocalBestP"),
    ]
    cell_deltas = []
    print("[step2] heatmaps", flush=True)
    for feat, out_name in flagship:
        if feat in df.columns or feat.startswith("neg_log10"):
            print(f"  fig02_{out_name}", flush=True)
            flat = step2_heatmap(df, feat, out_name)
            cell_deltas.append(flat)
    if cell_deltas:
        pd.concat(cell_deltas, ignore_index=True).to_csv(
            DATA_DIR / "G7_cell_delta.tsv", sep="\t", index=False)

    # Step 3 per-sample consistency -------------------------------------------
    print("[step3] per-sample consistency", flush=True)
    for feat, out_name in flagship[:3]:  # top-3 only to limit figure count
        step3_per_sample(df, feat, out_name)

    # Step 4 stratified AUC ---------------------------------------------------
    print("[step4] stratified AUC", flush=True)
    auc_rows = []
    for feat, out_name in flagship:
        r = step4_stratified(df, feat, out_name)
        if r is not None:
            auc_rows.append(r)
    if auc_rows:
        pd.concat(auc_rows, ignore_index=True).to_csv(
            DATA_DIR / "G7_auc_table.tsv", sep="\t", index=False)

    # Step 5 confound guard ---------------------------------------------------
    print("[step5] confound guard", flush=True)
    conf_rows = []
    for feat, out_name in flagship:
        c = step5_confound(df, feat, out_name)
        if c is not None:
            conf_rows.append(c)
    if conf_rows:
        pd.DataFrame(conf_rows).to_csv(DATA_DIR / "G7_confound.tsv", sep="\t", index=False)

    # Step 6 spatial ----------------------------------------------------------
    print("[step6] spatial", flush=True)
    for feat, out_name in flagship[:3]:
        step6_spatial(df, feat, out_name)

    # Extra diagnostics -------------------------------------------------------
    print("[extra] CramersV zero rate + HPFine replacement", flush=True)
    cramers_info = fig_cramers_zero_rate(df, hpfine_df=None)
    print("[extra] Pairwise dist collinearity", flush=True)
    coll_info = fig_pairwise_collinearity(df)
    print("[extra] PermanovaValid gate", flush=True)
    valid_info = fig_permanova_valid_gate(df)

    # Summary sidecar ---------------------------------------------------------
    summary = {
        "n_rows": len(df),
        "cramers_zero_rate": cramers_info["cramers_zero_rate"] if cramers_info else None,
        "cramers_auc_comparison": cramers_info["auc_comparison"] if cramers_info else None,
        "pairwise_spearman_mean_median": coll_info["spearman"] if coll_info else None,
        "permanova_valid_gate_delta_auc": valid_info["delta"] if valid_info else None,
    }
    pd.DataFrame([{"key": k, "value": str(v)} for k, v in summary.items()]).to_csv(
        DATA_DIR / "G7_summary.tsv", sep="\t", index=False)

    # Persist enriched subset for reproducibility (RegionID-level, small cols)
    keep_cols = ["sample", "mode", "RegionID", "Chr", "Pos", "tp_label",
                 "LOH_Subtype", "AF_class", "cn_tier_F", "vcf_AF", "NumReads",
                 "Coverage_Multiple"] + G7_RAW_COLS + [
                 "neg_log10_GlobalP", "neg_log10_ClusterPermP", "neg_log10_LocalBestP"]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df[keep_cols].to_csv(OUT_SUB / "G7_enriched.tsv.gz",
                         sep="\t", index=False, compression="gzip")

    print("[done] Outputs:")
    print(" ", FIG_DIR)
    print(" ", DATA_DIR)


if __name__ == "__main__":
    main()
