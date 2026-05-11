#!/usr/bin/env python3
"""Phase C · G8 · Entropy / Epipolymorphism / PerCpgASM feature group analysis.

Features (11):
  - NME_HP1, NME_HP2                         (Normalized Methylation Entropy per haplotype)
  - Entropy_Imbalance                        (|NME_HP1 - NME_HP2|)
  - Epipoly_HP1, Epipoly_HP2, Epipoly_Delta  (within-HP read-level pattern diversity)
  - PerCpgASM_Valid                          (boolean: eligible region for per-CpG Fisher ASM)
  - Fisher_N_Sig, Fisher_Frac_Sig            (per-site ASM counts / fraction significant)
  - Fisher_N_Tested                          (per-site ASM sites tested)
  - Fisher_MaxNegLogFDR                      (-log10 min FDR)

Source: raw significance_summary.csv per (sample, mode). Master TSV does not carry these cols.

SOP Step 1-6 (matches G5/G6 observe scripts):
  1. Global distribution (TP/FP + AUC with Wilson CI + Cohen's d)
  2. LOH × AF × CN stratification (32-cell heatmap; target: Fisher_Frac_Sig + Entropy_Imbalance)
  3. Per-sample consistency: 7-sample AUC + Spearman concordance
  4. Stratified AUC (LOH / AF / CN / mode / sample)
  5. Confound guard: within-(sample,mode,LOH) residualization on vcf_AF + Coverage_Multiple + NumReads
  6. Spatial autocorrelation (chr × 5 Mb bin)

Special asks from task spec:
  - Fisher_Frac_Sig paired AUC 0.726 cross-sample stability (7/7) + CI
  - Entropy_Imbalance vs AlleleDelta independence (residualize on AlleleDelta)
  - Epipoly_Delta vs NME_HP1/HP2 independence
  - Paired vs TO direction concordance

Outputs:
  - data/G8/*.tsv          per-step stats tables (including G8_auc_table.tsv headline)
  - figures/G8_entropy/*.png
  - features/G8_entropy_epipoly.md  (written separately)

Known sample gap:
  - HCC1954 paired_full raw has no G8 features (older pipeline). Gracefully degrades to
    6 samples paired_full + 7 samples TO = 13 sample-mode cells.

Author: InterSubMod Phase C
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ----------------------------------------------------------------------------
ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER_TSV = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
OUT_DATA = ROOT / "research/feature_layered_observation/data/G8"
OUT_FIG = ROOT / "research/feature_layered_observation/figures/G8_entropy"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
    "H2009", "H1437", "COLO829",
]

PALETTE = {
    "tp": "#5B8DB7",
    "fp": "#D55E00",
    "dark": "#1E2A44",
    "accent": "#A85540",
    "green": "#009E73",
    "grey": "#8A8A8A",
}

G8_FEATURES = [
    "NME_HP1", "NME_HP2", "Entropy_Imbalance",
    "Epipoly_HP1", "Epipoly_HP2", "Epipoly_Delta",
    "PerCpgASM_Valid",
    "Fisher_N_Sig", "Fisher_Frac_Sig", "Fisher_N_Tested", "Fisher_MaxNegLogFDR",
]

CONTINUOUS_G8 = [f for f in G8_FEATURES if f != "PerCpgASM_Valid"]

# Source paths (mirror G6_hp_fine.py)
PAIRED_FULL_PATHS = {
    "HCC1395":        ROOT / "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO": ROOT / "output/canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full",
    "H2009":          ROOT / "output/canonical/H2009/paired_full/20260421_H2009_paired_full_full",
    "H1437":          ROOT / "output/canonical/H1437/paired_full/20260421_H1437_paired_full_full",
    "COLO829":        ROOT / "output/canonical/COLO829/paired_full/20260421_COLO829_paired_full_full",
    "HCC1937":        ROOT / "output/canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full",
    # HCC1954 paired_full is older pipeline — no G8 features (gracefully skipped)
    "HCC1954":        ROOT / "output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix",
}

HCC1395_TO = Path("/tmp/ism_hp_fix_phase1")  # tp_off / fp_off
TO_ARCHIVE_PATHS = {
    "HCC1395_DORADO": ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod",
    "H1437":          ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod",
    "H2009":          ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod",
    "HCC1937":        ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod",
    "HCC1954":        ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod",
}
COLO829_TO = ROOT / "output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod"


# ----------------------------------------------------------------------------
# Loading helpers

def _load_pair(tp_csv: Path, fp_csv: Path) -> pd.DataFrame | None:
    if not (tp_csv.exists() and fp_csv.exists()):
        return None
    tp = pd.read_csv(tp_csv, low_memory=False)
    fp = pd.read_csv(fp_csv, low_memory=False)
    keep_base = ["RegionID", "Chr", "Pos"]
    present = [f for f in G8_FEATURES if f in tp.columns]
    if not present:
        return None  # sample lacks G8 features entirely (e.g. HCC1954 paired_full)
    keep = keep_base + present
    tp = tp[keep].copy(); tp["_tp"] = 1
    fp = fp[keep].copy(); fp["_tp"] = 0
    return pd.concat([tp, fp], ignore_index=True)


def enrich_g8(master: pd.DataFrame) -> pd.DataFrame:
    """Merge raw G8 features onto master by (sample, mode, RegionID, Chr, Pos)."""
    frames = []
    skipped = []
    for (sample, mode), _grp in master.groupby(["sample", "mode"], sort=False):
        df_raw = None
        if mode == "paired_full" and sample in PAIRED_FULL_PATHS:
            rd = PAIRED_FULL_PATHS[sample]
            df_raw = _load_pair(rd / "intersubmod_tp" / "significance_summary.csv",
                                rd / "intersubmod_fp" / "significance_summary.csv")
        elif mode == "to_pileup":
            if sample == "HCC1395":
                df_raw = _load_pair(HCC1395_TO / "tp_off" / "significance_summary.csv",
                                    HCC1395_TO / "fp_off" / "significance_summary.csv")
            elif sample == "COLO829":
                df_raw = _load_pair(COLO829_TO / "intersubmod_tp" / "significance_summary.csv",
                                    COLO829_TO / "intersubmod_fp" / "significance_summary.csv")
            elif sample in TO_ARCHIVE_PATHS:
                rd = TO_ARCHIVE_PATHS[sample].resolve()
                df_raw = _load_pair(rd / "intersubmod_tp" / "significance_summary.csv",
                                    rd / "intersubmod_fp" / "significance_summary.csv")
        if df_raw is None:
            skipped.append(f"({sample},{mode})")
            continue
        df_raw["sample"] = sample
        df_raw["mode"] = mode
        frames.append(df_raw)
    if skipped:
        print(f"[enrich] skipped (no G8 features): {skipped}")
    raw = pd.concat(frames, ignore_index=True)
    print(f"[enrich] raw G8 rows: {len(raw)}")
    merged = master.merge(
        raw,
        on=["sample", "mode", "RegionID", "Chr", "Pos"],
        how="left",
        suffixes=("", "_r"),
    )
    mismatch = ((merged["tp_label"] != merged["_tp"]) & merged["_tp"].notna()).sum()
    print(f"[enrich] tp_label mismatch vs _tp: {mismatch}")
    merged.drop(columns=["_tp"], inplace=True)
    # Coerce PerCpgASM_Valid to boolean int
    if "PerCpgASM_Valid" in merged.columns:
        merged["PerCpgASM_Valid_int"] = merged["PerCpgASM_Valid"].map(
            {True: 1, False: 0, "True": 1, "False": 0}
        )
    return merged


# ----------------------------------------------------------------------------
# Stats utilities

def _mw_auc(pos: np.ndarray, neg: np.ndarray) -> float:
    x = np.concatenate([pos, neg])
    r = stats.rankdata(x)
    r_pos = r[: len(pos)].sum()
    n1 = len(pos); n0 = len(neg)
    return (r_pos - n1 * (n1 + 1) / 2) / (n1 * n0)


def auc_wilson(y: np.ndarray, score: np.ndarray) -> Tuple[float, float, float]:
    """AUC + Hanley-McNeil 95% CI."""
    y = np.asarray(y); score = np.asarray(score, dtype=float)
    mask = ~np.isnan(score) & ~pd.isna(y)
    y = y[mask]; score = score[mask]
    if y.sum() == 0 or (1 - y).sum() == 0 or len(y) < 5:
        return (np.nan, np.nan, np.nan)
    pos = score[y == 1]; neg = score[y == 0]
    auc = _mw_auc(pos, neg)
    n1 = len(pos); n0 = len(neg)
    q1 = auc / (2 - auc); q2 = (2 * auc * auc) / (1 + auc)
    var = (auc * (1 - auc) + (n1 - 1) * (q1 - auc ** 2) + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0)
    se = np.sqrt(max(var, 0))
    return (auc, auc - 1.96 * se, auc + 1.96 * se)


def cohens_d(pos: np.ndarray, neg: np.ndarray) -> float:
    pos = pos[~np.isnan(pos)]; neg = neg[~np.isnan(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return np.nan
    s2 = ((len(pos) - 1) * np.var(pos, ddof=1) + (len(neg) - 1) * np.var(neg, ddof=1)) / (len(pos) + len(neg) - 2)
    if s2 <= 0:
        return np.nan
    return (np.mean(pos) - np.mean(neg)) / np.sqrt(s2)


def cn_tier_F(covm: pd.Series) -> pd.Series:
    edges = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
    labels = ["T0_sub", "T1_near", "T2_dup", "T3_quad", "T4_amp"]
    return pd.cut(covm.astype(float), bins=edges, labels=labels, include_lowest=True)


# ----------------------------------------------------------------------------
# Step 1: global distribution

def step1_global(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feat in G8_FEATURES:
        if feat not in df.columns:
            continue
        if feat == "PerCpgASM_Valid":
            s = pd.to_numeric(df.get("PerCpgASM_Valid_int"), errors="coerce")
        else:
            s = pd.to_numeric(df[feat], errors="coerce")
        y = df["tp_label"].astype(int).values
        pos = s[y == 1].values; neg = s[y == 0].values
        auc, lo, hi = auc_wilson(y, s.values)
        d = cohens_d(pos, neg)
        try:
            u, p = stats.mannwhitneyu(pos[~np.isnan(pos)], neg[~np.isnan(neg)], alternative="two-sided")
        except ValueError:
            u, p = np.nan, np.nan
        rows.append(dict(
            feature=feat,
            n_tp_valid=int((~np.isnan(pos)).sum()),
            n_fp_valid=int((~np.isnan(neg)).sum()),
            mean_tp=float(np.nanmean(pos)) if len(pos) else np.nan,
            mean_fp=float(np.nanmean(neg)) if len(neg) else np.nan,
            median_tp=float(np.nanmedian(pos)) if len(pos) else np.nan,
            median_fp=float(np.nanmedian(neg)) if len(neg) else np.nan,
            cohens_d=float(d) if not np.isnan(d) else np.nan,
            mwu_p=float(p) if not np.isnan(p) else np.nan,
            auc=float(auc) if not np.isnan(auc) else np.nan,
            auc_lo=float(lo) if not np.isnan(lo) else np.nan,
            auc_hi=float(hi) if not np.isnan(hi) else np.nan,
        ))
    out = pd.DataFrame(rows).sort_values("auc", ascending=False)
    out.to_csv(OUT_DATA / "G8_auc_table.tsv", sep="\t", index=False)
    out.to_csv(OUT_DATA / "step1_global_stats.tsv", sep="\t", index=False)
    return out


def fig_step1(stats_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    d = stats_df.dropna(subset=["auc"]).copy().sort_values("auc")
    y_pos = np.arange(len(d))
    xerr = np.vstack([(d["auc"] - d["auc_lo"]).clip(lower=0).values,
                      (d["auc_hi"] - d["auc"]).clip(lower=0).values])
    colors = [PALETTE["accent"] if v >= 0.58 else PALETTE["grey"] for v in d["auc"]]
    ax.barh(y_pos, d["auc"], xerr=xerr, color=colors, ecolor="#333", capsize=3)
    ax.set_yticks(y_pos); ax.set_yticklabels(d["feature"])
    ax.axvline(0.5, ls="--", c="#888", lw=1)
    ax.axvline(0.58, ls="--", c="#c00", lw=1, label="Beyond-AUC 0.58")
    ax.set_xlabel("AUC (TP vs FP, Hanley-McNeil 95% CI)")
    ax.set_xlim(0.35, 0.90)
    ax.set_title("G8 Entropy/Epipoly/ASM · Global AUC ranking")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig01_global_auc.png"); plt.close(fig)


def fig_step1_dist(df: pd.DataFrame, features: List[str]) -> None:
    n = len(features)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4.5), dpi=120)
    if n == 1:
        axes = [axes]
    for ax, feat in zip(axes, features):
        if feat not in df.columns:
            ax.text(0.5, 0.5, f"{feat} missing", ha="center"); continue
        sub = df[["tp_label", feat]].dropna()
        if len(sub) < 10:
            ax.text(0.5, 0.5, f"{feat}: insufficient", ha="center"); continue
        tp_vals = pd.to_numeric(sub.loc[sub["tp_label"] == 1, feat], errors="coerce").dropna()
        fp_vals = pd.to_numeric(sub.loc[sub["tp_label"] == 0, feat], errors="coerce").dropna()
        bins = 40
        ax.hist(fp_vals, bins=bins, alpha=0.55, color=PALETTE["fp"], label=f"FP (n={len(fp_vals)})", density=True)
        ax.hist(tp_vals, bins=bins, alpha=0.55, color=PALETTE["tp"], label=f"TP (n={len(tp_vals)})", density=True)
        ax.set_title(feat, fontsize=10)
        ax.set_xlabel("value"); ax.set_ylabel("density")
        ax.legend(fontsize=7)
    fig.suptitle("G8 · Step 1 · TP vs FP distributions", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig01b_dist_tp_fp.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Step 2: LOH × AF × CN heatmap

def step2_heatmap(df: pd.DataFrame, feature: str, out_prefix: str) -> pd.DataFrame:
    d = df.copy()
    d["CN_tier"] = cn_tier_F(d["Coverage_Multiple"])
    loh_order = ["None", "Noise", "Weak", "Strong", "Subclone"]
    d["LOH_Subtype"] = d.get("LOH_Subtype", pd.Series(["None"] * len(d))).fillna("None")
    af_order = ["Extreme", "Intermediate", "Near-half"]
    d["AF_class"] = d.get("AF_class", "Unknown").astype(str)

    # Drop rows with NaN feature
    d = d.dropna(subset=[feature])
    cells = d.groupby(["LOH_Subtype", "AF_class", "CN_tier"], observed=True).apply(
        lambda g: pd.Series({
            "n": len(g),
            "tp_rate": g["tp_label"].mean(),
            "feat_tp": g.loc[g["tp_label"] == 1, feature].mean(),
            "feat_fp": g.loc[g["tp_label"] == 0, feature].mean(),
        })
    ).reset_index()
    cells = cells[cells["n"] >= 20].copy()
    cells["delta"] = cells["feat_tp"] - cells["feat_fp"]
    cells.to_csv(OUT_DATA / f"step2_{out_prefix}_cells.tsv", sep="\t", index=False)

    loh_present = [x for x in loh_order if x in cells["LOH_Subtype"].unique()]
    cn_order = ["T0_sub", "T1_near", "T2_dup", "T3_quad", "T4_amp"]
    cn_present = [x for x in cn_order if x in cells["CN_tier"].astype(str).unique()]

    fig, axes = plt.subplots(1, 4, figsize=(22, 7), dpi=120)
    panels = [
        ("tp_rate", "TP rate", 0, 1, "viridis"),
        ("feat_tp", f"{feature} (TP only)", None, None, "Blues"),
        ("feat_fp", f"{feature} (FP only)", None, None, "Oranges"),
        ("delta",   f"Δ(TP−FP) {feature}", None, None, "RdBu_r"),
    ]
    for ax, (col, title, vmin, vmax, cmap) in zip(axes, panels):
        pivot = cells.pivot_table(
            index=["LOH_Subtype", "CN_tier"], columns="AF_class", values=col, aggfunc="first"
        )
        idx = [(l, c) for l in loh_present for c in cn_present if (l, c) in pivot.index]
        pivot = pivot.reindex(idx).reindex(columns=[c for c in af_order if c in pivot.columns])
        ncells = cells.pivot_table(
            index=["LOH_Subtype", "CN_tier"], columns="AF_class", values="n", aggfunc="first"
        ).reindex(idx).reindex(columns=[c for c in af_order if c in pivot.columns])
        arr = pivot.values.astype(float)
        if col == "delta":
            lim = max(abs(np.nanmin(arr)), abs(np.nanmax(arr)), 1e-6)
            im = ax.imshow(arr, cmap=cmap, vmin=-lim, vmax=lim, aspect="auto")
        else:
            im = ax.imshow(arr, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels([f"{l}/{c}" for l, c in pivot.index], fontsize=8)
        ax.set_title(title, fontsize=10)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = arr[i, j]; nv = ncells.values[i, j]
                if np.isnan(v):
                    continue
                ax.text(j, i, f"{v:.2f}\n({int(nv)})", ha="center", va="center", fontsize=6, color="#111")
        plt.colorbar(im, ax=ax, shrink=0.7)
    fig.suptitle(f"G8 Step 2 · LOH × AF × CN stratification · {feature}", fontsize=13)
    fig.tight_layout(); fig.savefig(OUT_FIG / f"fig02_{out_prefix}.png"); plt.close(fig)
    return cells


# ----------------------------------------------------------------------------
# Step 3: per-sample consistency

def step3_per_sample(df: pd.DataFrame, features: List[str]) -> pd.DataFrame:
    rows = []
    for feat in features:
        if feat not in df.columns:
            continue
        s_all = pd.to_numeric(df[feat], errors="coerce")
        for (sample, mode), g in df.groupby(["sample", "mode"]):
            s = pd.to_numeric(g[feat], errors="coerce")
            if s.notna().sum() < 30:
                continue
            y = g["tp_label"].astype(int).values
            auc, lo, hi = auc_wilson(y, s.values)
            pos = s[y == 1].values; neg = s[y == 0].values
            rows.append(dict(
                feature=feat, sample=sample, mode=mode,
                n=len(g), n_valid=int(s.notna().sum()),
                n_tp=int((y == 1).sum()), n_fp=int((y == 0).sum()),
                mean_tp=float(np.nanmean(pos)) if len(pos) else np.nan,
                mean_fp=float(np.nanmean(neg)) if len(neg) else np.nan,
                auc=float(auc) if not np.isnan(auc) else np.nan,
                auc_lo=float(lo) if not np.isnan(lo) else np.nan,
                auc_hi=float(hi) if not np.isnan(hi) else np.nan,
            ))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "step3_per_sample_auc.tsv", sep="\t", index=False)
    return out


def fig_step3(per_sample: pd.DataFrame, focus_feats: List[str]) -> None:
    n = len(focus_feats)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3.4 * n), dpi=120, sharex=True)
    if n == 1:
        axes = [axes]
    for ax, feat in zip(axes, focus_feats):
        sub = per_sample[per_sample["feature"] == feat]
        modes = {"paired_full": (PALETTE["tp"], -0.18), "to_pileup": (PALETTE["fp"], 0.18)}
        x = np.arange(len(SAMPLE_ORDER))
        for mode, (c, offs) in modes.items():
            dm = sub[sub["mode"] == mode].set_index("sample").reindex(SAMPLE_ORDER)
            ax.bar(x + offs, dm["auc"].values, 0.35, color=c, label=mode, alpha=0.9)
            err = np.vstack([(dm["auc"] - dm["auc_lo"]).clip(lower=0).fillna(0).values,
                             (dm["auc_hi"] - dm["auc"]).clip(lower=0).fillna(0).values])
            ax.errorbar(x + offs, dm["auc"].values, yerr=err, fmt="none", ecolor="#333", lw=0.8, capsize=2)
            for xi, yi, n_ in zip(x + offs, dm["auc"].values, dm["n_valid"].fillna(0).values):
                if not np.isnan(yi):
                    ax.text(xi, yi + 0.01, f"{int(n_)}", ha="center", fontsize=6, rotation=90)
        ax.axhline(0.5, ls="--", c="#888", lw=0.8)
        ax.axhline(0.58, ls="--", c="#c00", lw=0.8)
        ax.set_ylim(0.35, 0.92)
        ax.set_xticks(x); ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel(f"{feat}\nAUC")
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle("G8 Step 3 · Per-sample AUC (paired_full vs to_pileup)", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig03_per_sample_auc.png"); plt.close(fig)


def fig_step3_concordance(per_sample: pd.DataFrame, feature: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=120)
    for ax, mode in zip(axes, ["paired_full", "to_pileup"]):
        sub = per_sample[(per_sample["feature"] == feature) & (per_sample["mode"] == mode)]
        if len(sub) < 2:
            ax.set_visible(False); continue
        aucs = sub.set_index("sample")["auc"].reindex(SAMPLE_ORDER)
        n_above = int((aucs >= 0.58).sum())
        n_pos = int((aucs > 0.5).sum())
        total = int(aucs.notna().sum())
        ax.bar(range(len(aucs)), aucs.values,
               color=[PALETTE["accent"] if (v and v >= 0.58) else PALETTE["grey"] for v in aucs.values])
        ax.axhline(0.5, ls="--", c="#888", lw=0.8)
        ax.axhline(0.58, ls="--", c="#c00", lw=0.8)
        ax.set_ylim(0.35, 0.92)
        ax.set_xticks(range(len(aucs))); ax.set_xticklabels(aucs.index, rotation=30, ha="right", fontsize=8)
        ax.set_title(f"{mode}: {feature}\n{n_above}/{total} ≥ 0.58,  {n_pos}/{total} > 0.5")
        for i, v in enumerate(aucs.values):
            if not np.isnan(v):
                ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=7)
    fig.suptitle(f"G8 Step 3b · Cross-sample AUC · {feature}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_FIG / f"fig03b_concordance_{feature}.png")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Step 4: stratified AUC (LOH / AF / CN / mode / sample)

def step4_stratified_auc(df: pd.DataFrame, feature: str) -> pd.DataFrame:
    d = df.copy()
    d["CN_tier"] = cn_tier_F(d["Coverage_Multiple"])
    rows = []
    a, lo, hi = auc_wilson(d["tp_label"].values, pd.to_numeric(d[feature], errors="coerce").values)
    rows.append(dict(stratum="Global", value="all", n=len(d), auc=a, auc_lo=lo, auc_hi=hi))
    if "LOH_Subtype" in d.columns:
        for v, g in d.groupby("LOH_Subtype"):
            a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
            rows.append(dict(stratum="LOH_Subtype", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    for v, g in d.groupby("AF_class"):
        a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
        rows.append(dict(stratum="AF_class", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    for v, g in d.groupby("CN_tier", observed=True):
        a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
        rows.append(dict(stratum="CN_tier", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    for v, g in d.groupby("mode"):
        a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
        rows.append(dict(stratum="mode", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    for v, g in d.groupby("sample"):
        a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
        rows.append(dict(stratum="sample", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / f"step4_stratified_auc_{feature}.tsv", sep="\t", index=False)
    return out


def fig_step4(stratified: pd.DataFrame, feature: str) -> None:
    fig, ax = plt.subplots(figsize=(11, 7), dpi=120)
    d = stratified.dropna(subset=["auc"]).copy()
    labels = [f"{row['stratum']}:{row['value']} (n={int(row['n'])})" for _, row in d.iterrows()]
    y_pos = np.arange(len(d))
    xerr = np.vstack([(d["auc"] - d["auc_lo"]).clip(lower=0).values,
                      (d["auc_hi"] - d["auc"]).clip(lower=0).values])
    cmap = {"Global": "#1E2A44", "LOH_Subtype": "#A85540", "AF_class": "#009E73",
            "CN_tier": "#5B8DB7", "mode": "#D55E00", "sample": "#8A8A8A"}
    colors = [cmap.get(row["stratum"], "#444") for _, row in d.iterrows()]
    ax.barh(y_pos, d["auc"], xerr=xerr, color=colors, ecolor="#333", capsize=2)
    ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0.5, ls="--", c="#888"); ax.axvline(0.58, ls="--", c="#c00")
    ax.set_xlim(0.35, 0.92)
    ax.set_xlabel("AUC")
    ax.set_title(f"G8 Step 4 · Stratified AUC · {feature}")
    fig.tight_layout(); fig.savefig(OUT_FIG / f"fig04_stratified_auc_{feature}.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Step 5: confound guard
#
# Runs three residualization regimes per feature:
#   (a) raw AUC
#   (b) within-(sample, mode, LOH_Subtype) OLS residualization on [NumReads, vcf_AF, Coverage_Multiple]
#       (global pool after group de-meaning)
#   (c) AdditionalExtra: residualize Entropy_Imbalance on AlleleDelta (task spec ask)
#                        residualize Epipoly_Delta on NME_HP1, NME_HP2
# Plus AF-bin cross (Extreme / Intermediate / Near-half) raw AUC per feature.

def _within_group_residualize(d: pd.DataFrame, feature: str,
                              covariates: List[str],
                              group_keys: List[str]) -> np.ndarray:
    """Within-(group) OLS residualization — returns residual series aligned to d.index."""
    resid = pd.Series(np.nan, index=d.index, dtype=float)
    target = pd.to_numeric(d[feature], errors="coerce")
    for _, g in d.groupby(group_keys, observed=True):
        idx = g.index
        X_cols = covariates
        X = g[X_cols].apply(pd.to_numeric, errors="coerce")
        y = target.loc[idx]
        mask = X.notna().all(axis=1) & y.notna()
        if mask.sum() < 30 or X.loc[mask].std().min() == 0:
            continue
        Xm = X.loc[mask].values
        # Add intercept
        Xm = np.column_stack([np.ones(len(Xm)), Xm])
        ym = y.loc[mask].values
        try:
            beta, *_ = np.linalg.lstsq(Xm, ym, rcond=None)
            yhat = Xm @ beta
            resid.loc[mask[mask].index] = ym - yhat
        except np.linalg.LinAlgError:
            continue
    return resid.values


def step5_confound(df: pd.DataFrame, features: List[str]) -> pd.DataFrame:
    d = df.copy()
    for c in ["vcf_AF", "NumReads", "Coverage_Multiple", "AlleleDelta"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    if "LOH_Subtype" in d.columns:
        d["LOH_Subtype"] = d["LOH_Subtype"].fillna("None")

    rows = []
    for feat in features:
        if feat not in d.columns:
            continue
        s = pd.to_numeric(d[feat], errors="coerce")
        mask = s.notna()
        if mask.sum() < 200:
            continue
        # Raw
        y = d.loc[mask, "tp_label"].values
        a_raw, a_raw_lo, a_raw_hi = auc_wilson(y, s.loc[mask].values)

        # (b) within-(sample, mode, LOH_Subtype) OLS residualization on NumReads + vcf_AF + CovM
        covariates_b = [c for c in ["NumReads", "vcf_AF", "Coverage_Multiple"] if c in d.columns]
        resid_b = _within_group_residualize(d, feat, covariates_b, ["sample", "mode", "LOH_Subtype"])
        a_res_b, a_res_b_lo, a_res_b_hi = auc_wilson(
            d.loc[~np.isnan(resid_b), "tp_label"].values,
            resid_b[~np.isnan(resid_b)],
        )

        # (c) Feature-specific additional residualization
        a_res_extra = np.nan; a_res_extra_lo = np.nan; a_res_extra_hi = np.nan; extra_label = ""
        if feat == "Entropy_Imbalance" and "AlleleDelta" in d.columns:
            resid_c = _within_group_residualize(
                d, feat,
                covariates=["AlleleDelta"] + covariates_b,
                group_keys=["sample", "mode", "LOH_Subtype"],
            )
            a_res_extra, a_res_extra_lo, a_res_extra_hi = auc_wilson(
                d.loc[~np.isnan(resid_c), "tp_label"].values,
                resid_c[~np.isnan(resid_c)],
            )
            extra_label = "+AlleleDelta"
        elif feat == "Epipoly_Delta" and "NME_HP1" in d.columns and "NME_HP2" in d.columns:
            # residualize on NME_HP1, NME_HP2 + baseline covariates within groups
            d_tmp = d.copy()
            resid_c = _within_group_residualize(
                d_tmp, feat,
                covariates=["NME_HP1", "NME_HP2"] + covariates_b,
                group_keys=["sample", "mode", "LOH_Subtype"],
            )
            a_res_extra, a_res_extra_lo, a_res_extra_hi = auc_wilson(
                d_tmp.loc[~np.isnan(resid_c), "tp_label"].values,
                resid_c[~np.isnan(resid_c)],
            )
            extra_label = "+NME_HP1+NME_HP2"

        # AF-bin cross (raw AUC within each AF_class for reality check)
        af_results = {}
        for af_cls in ["Extreme", "Intermediate", "Near-half"]:
            m2 = mask & (d["AF_class"] == af_cls)
            if m2.sum() < 100:
                af_results[af_cls] = np.nan
                continue
            a_sub, _, _ = auc_wilson(d.loc[m2, "tp_label"].values, s.loc[m2].values)
            af_results[af_cls] = a_sub

        rows.append(dict(
            feature=feat,
            n_used=int(mask.sum()),
            auc_raw=a_raw, auc_raw_lo=a_raw_lo, auc_raw_hi=a_raw_hi,
            auc_resid_b=a_res_b, auc_resid_b_lo=a_res_b_lo, auc_resid_b_hi=a_res_b_hi,
            extra_label=extra_label,
            auc_resid_extra=a_res_extra, auc_resid_extra_lo=a_res_extra_lo, auc_resid_extra_hi=a_res_extra_hi,
            af_Extreme=af_results.get("Extreme", np.nan),
            af_Intermediate=af_results.get("Intermediate", np.nan),
            af_NearHalf=af_results.get("Near-half", np.nan),
        ))

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "step5_confound.tsv", sep="\t", index=False)
    return out


def fig_step5(conf_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11, 6.5), dpi=120)
    d = conf_df.dropna(subset=["auc_raw"]).sort_values("auc_raw")
    y_pos = np.arange(len(d))
    bar_h = 0.26
    ax.barh(y_pos - bar_h, d["auc_raw"], bar_h, color=PALETTE["accent"], label="raw AUC")
    ax.barh(y_pos, d["auc_resid_b"], bar_h, color=PALETTE["grey"],
            label="within-(sample,mode,LOH) resid on NumReads+vcf_AF+CovM")
    has_extra = d["auc_resid_extra"].notna()
    if has_extra.any():
        ax.barh(y_pos[has_extra] + bar_h, d.loc[has_extra, "auc_resid_extra"], bar_h,
                color=PALETTE["green"], label="+extra covariate (see label)")
        for yi, lbl in zip(y_pos[has_extra], d.loc[has_extra, "extra_label"].values):
            ax.text(0.36, yi + bar_h, f"{lbl}", fontsize=6, va="center", color=PALETTE["green"])
    ax.axvline(0.5, ls="--", c="#888"); ax.axvline(0.58, ls="--", c="#c00")
    ax.set_yticks(y_pos); ax.set_yticklabels(d["feature"], fontsize=9)
    ax.set_xlim(0.30, 0.90); ax.set_xlabel("AUC")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("G8 Step 5 · Confound guard (raw vs within-group residualized)")
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig05_confound_residualized.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Step 6: spatial autocorrelation (chr × 5 Mb bin)

def step6_spatial(df: pd.DataFrame, feature: str) -> pd.DataFrame:
    d = df[["Chr", "Pos", "tp_label", feature]].dropna().copy()
    if len(d) == 0:
        return pd.DataFrame()
    d["bin"] = (d["Pos"] // 5_000_000).astype(int)
    rows = []
    for (chrom, b), g in d.groupby(["Chr", "bin"]):
        if len(g) < 80 or g["tp_label"].sum() == 0 or (1 - g["tp_label"]).sum() == 0:
            continue
        auc, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
        rows.append(dict(chrom=chrom, bin=b, n=len(g), tp_rate=float(g["tp_label"].mean()),
                         auc=auc, auc_lo=lo, auc_hi=hi))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / f"step6_spatial_{feature}.tsv", sep="\t", index=False)
    return out


def fig_step6(sp: pd.DataFrame, feature: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 6), dpi=120)
    if sp.empty:
        ax.text(0.5, 0.5, "insufficient per-bin data", ha="center")
    else:
        sc = ax.scatter(sp["tp_rate"], sp["auc"], c=sp["n"], cmap="viridis", alpha=0.7, s=30)
        ax.axhline(0.5, ls="--", c="#888"); ax.axhline(0.58, ls="--", c="#c00")
        ax.set_xlim(0, 1.02); ax.set_ylim(0.3, 1.0)
        ax.set_xlabel("per-bin TP rate (baseline)")
        ax.set_ylabel(f"per-bin AUC · {feature}")
        plt.colorbar(sc, ax=ax, label="n variants per 5 Mb bin")
        mid = sp[(sp["tp_rate"] >= 0.4) & (sp["tp_rate"] <= 0.7)]
        if len(mid) > 0:
            ax.annotate(
                f"mid-TP-rate bins (0.4-0.7): n={len(mid)}, median AUC {mid['auc'].median():.3f}",
                xy=(0.05, 0.93), xycoords="axes fraction", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", edgecolor="#888"),
            )
    ax.set_title(f"G8 Step 6 · Spatial autocorrelation AUC vs TP rate · {feature}")
    fig.tight_layout(); fig.savefig(OUT_FIG / f"fig06_spatial_{feature}.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Extra analyses for task asks

def fisher_frac_sig_cross_sample(per_sample: pd.DataFrame) -> pd.DataFrame:
    """7/7 cross-sample stability table for Fisher_Frac_Sig."""
    sub = per_sample[per_sample["feature"] == "Fisher_Frac_Sig"].copy()
    out = sub.sort_values(["mode", "sample"])[
        ["mode", "sample", "n", "n_valid", "mean_tp", "mean_fp", "auc", "auc_lo", "auc_hi"]
    ]
    out.to_csv(OUT_DATA / "fisher_frac_sig_crosssample.tsv", sep="\t", index=False)
    return out


def epipoly_vs_nme_independence(df: pd.DataFrame) -> pd.DataFrame:
    """Does Epipoly_Delta carry signal independent of NME_HP1/NME_HP2?

    Approach: compute Pearson + Spearman correlation Epipoly_Delta vs
      (a) |NME_HP1-NME_HP2| = Entropy_Imbalance
      (b) NME_HP1
      (c) NME_HP2
    """
    d = df.dropna(subset=["Epipoly_Delta", "NME_HP1", "NME_HP2", "Entropy_Imbalance"])
    rows = []
    rows.append(dict(pair="Epipoly_Delta vs Entropy_Imbalance",
                     pearson=d["Epipoly_Delta"].corr(d["Entropy_Imbalance"]),
                     spearman=d["Epipoly_Delta"].corr(d["Entropy_Imbalance"], method="spearman"),
                     n=len(d)))
    rows.append(dict(pair="Epipoly_Delta vs NME_HP1",
                     pearson=d["Epipoly_Delta"].corr(d["NME_HP1"]),
                     spearman=d["Epipoly_Delta"].corr(d["NME_HP1"], method="spearman"),
                     n=len(d)))
    rows.append(dict(pair="Epipoly_Delta vs NME_HP2",
                     pearson=d["Epipoly_Delta"].corr(d["NME_HP2"]),
                     spearman=d["Epipoly_Delta"].corr(d["NME_HP2"], method="spearman"),
                     n=len(d)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "epipoly_vs_nme_correlation.tsv", sep="\t", index=False)
    return out


def entropy_vs_allele_delta(df: pd.DataFrame) -> pd.DataFrame:
    """Entropy_Imbalance vs AlleleDelta (concerns: Entropy_Imbalance may ≈ |m1-m2| after entropy xform)."""
    if "AlleleDelta" not in df.columns:
        return pd.DataFrame()
    d = df.dropna(subset=["Entropy_Imbalance", "AlleleDelta"]).copy()
    d["AlleleDelta_abs"] = d["AlleleDelta"].abs()
    rows = []
    for (scope, keyvals) in [
        ("Global", [("all", d)]),
        ("mode", list(d.groupby("mode"))),
        ("LOH_Subtype", list(d.groupby(d["LOH_Subtype"].fillna("None")))),
    ]:
        for k, g in keyvals:
            if len(g) < 30:
                continue
            rows.append(dict(
                scope=scope, value=str(k), n=len(g),
                pearson_absAD=g["Entropy_Imbalance"].corr(g["AlleleDelta_abs"]),
                spearman_absAD=g["Entropy_Imbalance"].corr(g["AlleleDelta_abs"], method="spearman"),
                pearson_AD=g["Entropy_Imbalance"].corr(g["AlleleDelta"]),
            ))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "entropy_vs_alleledelta.tsv", sep="\t", index=False)
    return out


def fig_paired_vs_to_scatter(per_sample: pd.DataFrame, features: List[str]) -> None:
    n = len(features)
    cols = min(n, 3); rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), dpi=120, squeeze=False)
    for ax, feat in zip(axes.flatten(), features):
        sub = per_sample[per_sample["feature"] == feat]
        p = sub[sub["mode"] == "paired_full"].set_index("sample")["auc"]
        t = sub[sub["mode"] == "to_pileup"].set_index("sample")["auc"]
        common = [s for s in SAMPLE_ORDER if s in p.index and s in t.index and not np.isnan(p[s]) and not np.isnan(t[s])]
        if not common:
            ax.set_visible(False); continue
        px = p.loc[common].values; tx = t.loc[common].values
        ax.scatter(px, tx, s=60, c=PALETTE["dark"], alpha=0.8)
        for s, xv, yv in zip(common, px, tx):
            ax.annotate(s, (xv, yv), fontsize=7, xytext=(4, 3), textcoords="offset points")
        lo = min(px.min(), tx.min(), 0.4) - 0.02
        hi = max(px.max(), tx.max(), 0.6) + 0.02
        ax.plot([lo, hi], [lo, hi], ls="--", c="#888", lw=0.8)
        ax.axhline(0.5, ls=":", c="#999", lw=0.6); ax.axvline(0.5, ls=":", c="#999", lw=0.6)
        ax.axhline(0.58, ls=":", c="#c00", lw=0.6); ax.axvline(0.58, ls=":", c="#c00", lw=0.6)
        r = np.corrcoef(px, tx)[0, 1] if len(common) >= 3 else np.nan
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("paired_full AUC"); ax.set_ylabel("to_pileup AUC")
        ax.set_title(f"{feat}\nPearson r={r:.2f} (n={len(common)})", fontsize=9)
    for ax in axes.flatten()[len(features):]:
        ax.set_visible(False)
    fig.suptitle("G8 · Paired vs TO per-sample AUC concordance", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig07_paired_vs_to_scatter.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# PerCpgASM_Valid as gate

def step_percpgasm_gate(df: pd.DataFrame) -> pd.DataFrame:
    """How does TP rate change conditional on PerCpgASM_Valid = True vs False?"""
    rows = []
    if "PerCpgASM_Valid" not in df.columns:
        return pd.DataFrame()
    d = df.copy()
    d["PerCpgASM_Valid_norm"] = d["PerCpgASM_Valid"].map({True: True, False: False,
                                                          "True": True, "False": False})
    for (sample, mode), g in d.groupby(["sample", "mode"]):
        for v in [True, False]:
            sub = g[g["PerCpgASM_Valid_norm"] == v]
            if len(sub) == 0:
                continue
            rows.append(dict(
                sample=sample, mode=mode, PerCpgASM_Valid=v,
                n=len(sub), n_tp=int(sub["tp_label"].sum()),
                n_fp=int((1 - sub["tp_label"]).sum()),
                tp_rate=float(sub["tp_label"].mean()),
            ))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "step_percpgasm_gate.tsv", sep="\t", index=False)
    return out


def fig_percpgasm(gate: pd.DataFrame) -> None:
    if gate.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120, sharey=True)
    for ax, mode in zip(axes, ["paired_full", "to_pileup"]):
        sub = gate[gate["mode"] == mode]
        x = np.arange(len(SAMPLE_ORDER))
        width = 0.38
        for i, v in enumerate([True, False]):
            dm = sub[sub["PerCpgASM_Valid"] == v].set_index("sample").reindex(SAMPLE_ORDER)
            ax.bar(x + (i - 0.5) * width, dm["tp_rate"].values, width,
                   label=f"Valid={v}", color=PALETTE["tp"] if v else PALETTE["grey"], alpha=0.9)
            for xi, yi, ni in zip(x + (i - 0.5) * width, dm["tp_rate"].values, dm["n"].fillna(0).values):
                if not np.isnan(yi):
                    ax.text(xi, yi + 0.02, f"n={int(ni)}", ha="center", fontsize=6, rotation=90)
        ax.set_xticks(x); ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right", fontsize=8)
        ax.set_ylim(0, 1.05); ax.axhline(0.5, ls=":", c="#888")
        ax.set_title(f"{mode} · TP rate | PerCpgASM_Valid"); ax.set_ylabel("TP rate")
        ax.legend(loc="lower right", fontsize=8)
    fig.suptitle("G8 · PerCpgASM_Valid as eligibility gate", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig08_percpgasm_gate.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Main

def main():
    print("[G8] loading master ...")
    master = pd.read_csv(MASTER_TSV, sep="\t", low_memory=False)
    print(f"[G8] master shape: {master.shape}")

    print("[G8] enriching with raw G8 features ...")
    df = enrich_g8(master)
    print(f"[G8] enriched shape: {df.shape}")
    present = [c for c in G8_FEATURES if c in df.columns]
    print(f"[G8] G8 features present: {len(present)}/{len(G8_FEATURES)}: {present}")

    # Provenance dump
    keep_cols = ["sample", "mode", "tp_label", "Chr", "Pos", "RegionID",
                 "LOH_Subtype", "AF_class", "vcf_AF", "AF", "NumReads",
                 "Coverage_Multiple", "Potential_LOH", "AlleleDelta"] + present
    keep_cols = [c for c in keep_cols if c in df.columns]
    df[keep_cols].to_csv(OUT_DATA / "g8_enriched.tsv.gz", sep="\t", index=False, compression="gzip")

    # Step 1
    print("[G8] Step 1 · global distribution ...")
    s1 = step1_global(df)
    fig_step1(s1)
    fig_step1_dist(df, ["Fisher_Frac_Sig", "Entropy_Imbalance", "Epipoly_Delta", "NME_HP1"])

    # Step 2: heatmaps for top-interest features
    print("[G8] Step 2 · LOH × AF × CN heatmaps ...")
    for feat, prefix in [
        ("Fisher_Frac_Sig", "fisher_frac_sig"),
        ("Entropy_Imbalance", "entropy_imb"),
        ("Epipoly_Delta", "epipoly_delta"),
    ]:
        if feat in df.columns:
            step2_heatmap(df, feat, prefix)

    # Step 3: per-sample consistency
    print("[G8] Step 3 · per-sample consistency ...")
    ps = step3_per_sample(df, CONTINUOUS_G8)
    fig_step3(ps, ["Fisher_Frac_Sig", "Entropy_Imbalance", "Epipoly_Delta"])
    fig_step3_concordance(ps, "Fisher_Frac_Sig")
    fig_step3_concordance(ps, "Entropy_Imbalance")
    fig_step3_concordance(ps, "Epipoly_Delta")

    # Step 4
    print("[G8] Step 4 · stratified AUC ...")
    for feat in ["Fisher_Frac_Sig", "Entropy_Imbalance", "Epipoly_Delta", "NME_HP1", "NME_HP2", "Fisher_MaxNegLogFDR"]:
        if feat in df.columns:
            s4 = step4_stratified_auc(df, feat)
            fig_step4(s4, feat)

    # Step 5: confound guard
    print("[G8] Step 5 · confound guard ...")
    conf_feats = [f for f in CONTINUOUS_G8 if f in df.columns]
    s5 = step5_confound(df, conf_feats)
    fig_step5(s5)

    # Step 6: spatial
    print("[G8] Step 6 · spatial autocorrelation ...")
    sp = step6_spatial(df, "Fisher_Frac_Sig")
    fig_step6(sp, "Fisher_Frac_Sig")
    sp2 = step6_spatial(df, "Entropy_Imbalance")
    fig_step6(sp2, "Entropy_Imbalance")

    # Extra analyses for task asks
    print("[G8] Extra · Fisher_Frac_Sig cross-sample table ...")
    ff = fisher_frac_sig_cross_sample(ps)
    print(ff.to_string(index=False))

    print("[G8] Extra · Epipoly vs NME correlation ...")
    corr_en = epipoly_vs_nme_independence(df)
    print(corr_en.to_string(index=False))

    print("[G8] Extra · Entropy_Imbalance vs AlleleDelta ...")
    ea = entropy_vs_allele_delta(df)
    if not ea.empty:
        print(ea.head(12).to_string(index=False))

    print("[G8] Extra · paired vs TO scatter ...")
    fig_paired_vs_to_scatter(ps, ["Fisher_Frac_Sig", "Entropy_Imbalance", "Epipoly_Delta",
                                  "NME_HP1", "NME_HP2", "Fisher_MaxNegLogFDR"])

    print("[G8] Extra · PerCpgASM_Valid gate ...")
    gate = step_percpgasm_gate(df)
    fig_percpgasm(gate)

    n_fig = len(list(OUT_FIG.glob("*.png")))
    print(f"[G8] done — {n_fig} figures in {OUT_FIG}")
    print(f"[G8] data tsv in {OUT_DATA}")


if __name__ == "__main__":
    main()
