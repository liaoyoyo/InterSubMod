#!/usr/bin/env python3
"""Phase C · G6 · HP Fine (4-bucket LOH-phasing) feature group analysis.

Inputs:
  - Master TSV: research/feature_layered_observation/data/merged_with_vcf.tsv.gz
    (contains HPFineNGroups, HPFine_NGroups_CF + LOH/AF/CN axes + tp_label)
  - Raw significance_summary.csv per (sample, mode) run: full G6 columns
    (HPFineF/P/Sig, HPFineN_*, HPFineD_*, GlobalP_HPFine, CramersV_HPFine)

Outputs (research/feature_layered_observation/):
  - data/G6/*.tsv               statistics per step
  - figures/G6_hp_fine/*.png    15+ figures across Step 1-6 + HPFineD panels + fingerprint
  - features/G6_hp_fine.md      10-section report (written by caller)

Methodology: research/feature_layered_observation/02_methodology.md Step 1-6.
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
from matplotlib.colors import TwoSlopeNorm
from scipy import stats

# ----------------------------------------------------------------------------
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER_TSV = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
OUT_DATA = ROOT / "research/feature_layered_observation/data/G6"
OUT_FIG = ROOT / "research/feature_layered_observation/figures/G6_hp_fine"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
    "H2009", "H1437", "COLO829",
]

PALETTE = {
    "tp": "#5B8DB7",   # blue
    "fp": "#D55E00",   # orange-red
    "dark": "#1E2A44",
    "accent": "#A85540",
    "green": "#009E73",
    "grey": "#8A8A8A",
}

# Full G6 feature list
G6_FEATURES = [
    "GlobalP_HPFine", "CramersV_HPFine", "HPFine_NGroups_CF",
    "HPFineF", "HPFineP", "HPFineSig", "HPFineNGroups",
    "HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S",
    "HPFineD_HP1_HP1S", "HPFineD_HP1_HP2", "HPFineD_HP1_HP2S",
    "HPFineD_HP1S_HP2", "HPFineD_HP1S_HP2S", "HPFineD_HP2_HP2S",
]
HPFINE_D_PAIRS = [
    "HPFineD_HP1_HP1S", "HPFineD_HP1_HP2", "HPFineD_HP1_HP2S",
    "HPFineD_HP1S_HP2", "HPFineD_HP1S_HP2S", "HPFineD_HP2_HP2S",
]
HPFINE_N_COLS = ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]

# Source paths for enrichment (from utils_io.py)
NEW_KDE_PAIRED_FULL = {
    "HCC1395": ROOT / "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO": ROOT / "output/canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full",
    "H2009": ROOT / "output/canonical/H2009/paired_full/20260421_H2009_paired_full_full",
    "H1437": ROOT / "output/canonical/H1437/paired_full/20260421_H1437_paired_full_full",
    "COLO829": ROOT / "output/canonical/COLO829/paired_full/20260421_COLO829_paired_full_full",
    "HCC1937": ROOT / "output/canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full",
    "HCC1954": ROOT / "output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix",
}
HCC1395_TO = Path("/tmp/ism_hp_fix_phase1")  # tp_off / fp_off
NEW_KDE_TO_ARCHIVE = {
    "HCC1395_DORADO": ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod",
    "H1437":          ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod",
    "H2009":          ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod",
    "HCC1937":        ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod",
    "HCC1954":        ROOT / "../big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod",
}
COLO829_TO = ROOT / "output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod"


# ----------------------------------------------------------------------------
# Loading

def load_g6_raw(run_dir: Path) -> pd.DataFrame:
    tp = pd.read_csv(run_dir / "intersubmod_tp" / "significance_summary.csv", low_memory=False)
    fp = pd.read_csv(run_dir / "intersubmod_fp" / "significance_summary.csv", low_memory=False)
    keep = ["RegionID", "Chr", "Pos"] + [c for c in G6_FEATURES if c in tp.columns]
    tp = tp[keep].copy(); tp["_tp"] = 1
    fp = fp[keep].copy(); fp["_tp"] = 0
    df = pd.concat([tp, fp], ignore_index=True)
    return df


def enrich_g6(master: pd.DataFrame) -> pd.DataFrame:
    """Merge G6 raw features onto master by (sample, mode, RegionID, Chr, Pos)."""
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
                    keep = ["RegionID", "Chr", "Pos"] + [c for c in G6_FEATURES if c in tp.columns] + ["_tp"]
                    df_raw = pd.concat([tp[keep], fp[keep]], ignore_index=True)
                    df_raw["sample"] = sample; df_raw["mode"] = mode
                    frames.append(df_raw); continue
            elif sample == "COLO829":
                run_dir = COLO829_TO
            elif sample in NEW_KDE_TO_ARCHIVE:
                run_dir = NEW_KDE_TO_ARCHIVE[sample].resolve()
        if run_dir is None:
            print(f"[warn] no run_dir for ({sample},{mode})")
            continue
        try:
            df_raw = load_g6_raw(Path(run_dir))
        except FileNotFoundError as e:
            print(f"[warn] missing raw for ({sample},{mode}): {e}")
            continue
        df_raw["sample"] = sample; df_raw["mode"] = mode
        frames.append(df_raw)
    raw = pd.concat(frames, ignore_index=True)
    print(f"[enrich] raw G6 rows: {len(raw)}")
    # Avoid col collisions; drop the 2 already in master then merge
    dupe_cols = [c for c in ["HPFineNGroups", "HPFine_NGroups_CF"] if c in master.columns and c in raw.columns]
    merged = master.merge(
        raw.rename(columns={c: f"__raw_{c}" for c in dupe_cols}),
        on=["sample", "mode", "RegionID", "Chr", "Pos"],
        how="left", suffixes=("", "_r"),
    )
    # Replace master columns with raw values where raw is present (ensures consistency)
    for c in dupe_cols:
        r = f"__raw_{c}"
        if r in merged.columns:
            merged[c] = merged[r].combine_first(merged[c])
            merged.drop(columns=[r], inplace=True)
    # Sanity vs tp_label
    merged["tp_label_raw"] = merged["_tp"]
    mismatch = (merged["tp_label"] != merged["_tp"]).sum()
    print(f"[enrich] merge done; tp_label raw vs master mismatch: {mismatch}")
    merged.drop(columns=["_tp"], inplace=True)
    return merged


# ----------------------------------------------------------------------------
# Utility stats

def auc_wilson(y: np.ndarray, score: np.ndarray) -> Tuple[float, float, float]:
    """AUC + Wilson-style 95% CI via Hanley-McNeil approx."""
    y = np.asarray(y); score = np.asarray(score)
    mask = ~np.isnan(score)
    y = y[mask]; score = score[mask]
    if y.sum() == 0 or (1 - y).sum() == 0 or len(y) < 5:
        return (np.nan, np.nan, np.nan)
    pos = score[y == 1]; neg = score[y == 0]
    auc = _mw_auc(pos, neg)
    n1 = len(pos); n0 = len(neg)
    # Hanley-McNeil SE
    q1 = auc / (2 - auc); q2 = (2 * auc * auc) / (1 + auc)
    var = (auc * (1 - auc) + (n1 - 1) * (q1 - auc ** 2) + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0)
    se = np.sqrt(max(var, 0))
    return (auc, auc - 1.96 * se, auc + 1.96 * se)


def _mw_auc(pos: np.ndarray, neg: np.ndarray) -> float:
    x = np.concatenate([pos, neg])
    r = stats.rankdata(x)
    r_pos = r[: len(pos)].sum()
    n1 = len(pos); n0 = len(neg)
    return (r_pos - n1 * (n1 + 1) / 2) / (n1 * n0)


def cohens_d(pos: np.ndarray, neg: np.ndarray) -> float:
    pos = pos[~np.isnan(pos)]; neg = neg[~np.isnan(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return np.nan
    s2 = ((len(pos) - 1) * np.var(pos, ddof=1) + (len(neg) - 1) * np.var(neg, ddof=1)) / (len(pos) + len(neg) - 2)
    if s2 <= 0:
        return np.nan
    return (np.mean(pos) - np.mean(neg)) / np.sqrt(s2)


def cn_tier_F(covm: pd.Series) -> pd.Series:
    # F_SEQC2_grounded: [0.65, 0.99, 1.33, 1.82]
    edges = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
    labels = ["T0_sub", "T1_near", "T2_dup", "T3_quad", "T4_amp"]
    return pd.cut(covm.astype(float), bins=edges, labels=labels, include_lowest=True)


# ----------------------------------------------------------------------------
# Step 1: global distribution

def step1_global(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # Restrict to continuous/ordinal numeric features
    for feat in G6_FEATURES:
        if feat not in df.columns:
            continue
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
            n_tp=int((y == 1).sum()),
            n_fp=int((y == 0).sum()),
            mean_tp=float(np.nanmean(pos)), mean_fp=float(np.nanmean(neg)),
            median_tp=float(np.nanmedian(pos)), median_fp=float(np.nanmedian(neg)),
            cohens_d=float(d) if not np.isnan(d) else np.nan,
            mwu_p=float(p) if not np.isnan(p) else np.nan,
            auc=float(auc) if not np.isnan(auc) else np.nan,
            auc_lo=float(lo) if not np.isnan(lo) else np.nan,
            auc_hi=float(hi) if not np.isnan(hi) else np.nan,
        ))
    out = pd.DataFrame(rows).sort_values("auc", ascending=False)
    out.to_csv(OUT_DATA / "step1_global_stats.tsv", sep="\t", index=False)
    return out


def fig_step1(stats_df: pd.DataFrame) -> None:
    # Figure 01: AUC bar with CI
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    d = stats_df.dropna(subset=["auc"]).copy().sort_values("auc")
    y_pos = np.arange(len(d))
    xerr = np.vstack([(d["auc"] - d["auc_lo"]).clip(lower=0).values, (d["auc_hi"] - d["auc"]).clip(lower=0).values])
    colors = [PALETTE["accent"] if v >= 0.58 else PALETTE["grey"] for v in d["auc"]]
    ax.barh(y_pos, d["auc"], xerr=xerr, color=colors, ecolor="#333", capsize=3)
    ax.set_yticks(y_pos); ax.set_yticklabels(d["feature"])
    ax.axvline(0.5, ls="--", c="#888", lw=1); ax.axvline(0.58, ls="--", c="#c00", lw=1, label="Beyond-AUC 0.58")
    ax.set_xlabel("AUC (TP vs FP, 95% CI Hanley-McNeil)"); ax.set_xlim(0.40, 0.85)
    ax.set_title("G6 HP Fine · Global AUC ranking (all 7 samples × 2 modes pooled)")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig01_global_auc.png"); plt.close(fig)

    # Figure 01b: HPFineNGroups distribution TP vs FP
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), dpi=120, sharey=True)
    for ax, feat in zip(axes, ["HPFineNGroups", "HPFine_NGroups_CF"]):
        if feat not in stats_df["feature"].values:
            continue
        # placeholder: plot from master
    # We need real data; handled below in main with df
    fig.tight_layout(); plt.close(fig)


def fig_step1_ng_dist(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), dpi=120)
    for ax, feat in zip(axes, ["HPFineNGroups", "HPFine_NGroups_CF"]):
        if feat not in df.columns:
            ax.text(0.5, 0.5, f"{feat} missing", ha="center"); continue
        sub = df[["tp_label", feat]].dropna()
        sub[feat] = sub[feat].astype(int)
        order = sorted(sub[feat].unique())
        tp_rate = sub.groupby(feat)["tp_label"].mean()
        n_per = sub.groupby(feat).size()
        ax.bar(tp_rate.index - 0.15, tp_rate.values, 0.3, color=PALETTE["tp"], label="TP rate")
        ax2 = ax.twinx()
        ax2.bar(n_per.index + 0.15, n_per.values, 0.3, color=PALETTE["grey"], alpha=0.6, label="n")
        ax.set_xlabel(feat); ax.set_ylabel("TP rate"); ax2.set_ylabel("n (variants)")
        ax.set_xticks(list(order))
        ax.set_ylim(0, 1); ax.axhline(sub["tp_label"].mean(), ls="--", c="#888", lw=1)
        ax.set_title(f"{feat}: TP rate by value (global pool)")
        for v, r, n in zip(order, tp_rate.values, n_per.values):
            ax.text(v, r + 0.02, f"{r:.2f}\nn={n}", ha="center", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig01b_NG_distribution.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Step 2: LOH × AF × CN heatmaps (primary: HPFineNGroups as feature of interest)

def step2_heatmap(df: pd.DataFrame, feature: str, out_prefix: str) -> pd.DataFrame:
    """32-cell heatmap (LOH × AF × CN) produced for TP rate + feature mean TP/FP + delta."""
    d = df.copy()
    d["CN_tier"] = cn_tier_F(d["Coverage_Multiple"])
    loh_order = ["None", "Noise", "Weak", "Strong", "Subclone"]
    if "LOH_Subtype" in d.columns:
        d["LOH_Subtype"] = d["LOH_Subtype"].fillna("None")
    else:
        d["LOH_Subtype"] = "None"
    af_order = ["Extreme", "Intermediate", "Near-half"]
    # Build long form
    d["AF_class"] = d.get("AF_class", "Unknown").astype(str)
    cells = d.groupby(["LOH_Subtype", "AF_class", "CN_tier"], observed=True).agg(
        n=("tp_label", "size"),
        tp_rate=("tp_label", "mean"),
        feat_tp=(feature, lambda s: s[d.loc[s.index, "tp_label"] == 1].mean()),
        feat_fp=(feature, lambda s: s[d.loc[s.index, "tp_label"] == 0].mean()),
    ).reset_index()
    cells = cells[cells["n"] >= 20].copy()
    cells["delta"] = cells["feat_tp"] - cells["feat_fp"]
    cells.to_csv(OUT_DATA / f"step2_{out_prefix}_cells.tsv", sep="\t", index=False)
    # Plot heatmaps as 2D: rows = LOH_Subtype (5) × CN_tier (5) = 25 rows, cols = AF_class (3)
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
    fig.suptitle(f"G6 Step 2 · LOH × AF × CN stratification · {feature}", fontsize=13)
    fig.tight_layout(); fig.savefig(OUT_FIG / f"fig02_{out_prefix}.png"); plt.close(fig)
    return cells


# ----------------------------------------------------------------------------
# Step 3: per-sample consistency (TP rate by NG value) + concordance matrix

def step3_per_sample(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["HPFineNGroups"] = pd.to_numeric(d["HPFineNGroups"], errors="coerce").astype("Int64")
    grp = d.groupby(["sample", "mode", "HPFineNGroups"], observed=True).agg(
        n=("tp_label", "size"),
        tp_rate=("tp_label", "mean"),
    ).reset_index()
    grp.to_csv(OUT_DATA / "step3_per_sample_ng.tsv", sep="\t", index=False)
    return grp


def fig_step3(per_sample: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), dpi=120, sharex=True)
    for ax, mode, title in zip(axes, ["paired_full", "to_pileup"], ["Paired Full", "TO Pileup"]):
        sub = per_sample[per_sample["mode"] == mode]
        x = sorted(sub["HPFineNGroups"].dropna().unique().tolist())
        width = 0.11
        for i, sample in enumerate(SAMPLE_ORDER):
            row = sub[sub["sample"] == sample].set_index("HPFineNGroups")
            yv = [row.loc[v, "tp_rate"] if v in row.index else np.nan for v in x]
            nv = [row.loc[v, "n"] if v in row.index else 0 for v in x]
            xp = np.array(x) + (i - 3) * width
            bars = ax.bar(xp, yv, width, label=sample)
            for xi, yi, ni in zip(xp, yv, nv):
                if not np.isnan(yi) and ni >= 10:
                    ax.text(xi, yi + 0.01, f"{ni}", ha="center", fontsize=6, rotation=90)
        ax.set_xticks(x); ax.set_xticklabels([str(v) for v in x])
        ax.set_ylabel("TP rate"); ax.set_ylim(0, 1.05)
        ax.set_title(f"{title} · TP rate by HPFineNGroups per sample")
        ax.axhline(0.5, ls=":", c="#888", lw=0.8)
    axes[0].legend(ncol=7, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.12))
    axes[-1].set_xlabel("HPFineNGroups (n non-empty buckets among HP1/HP1-1/HP2/HP2-1)")
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig03_per_sample_ng_tprate.png"); plt.close(fig)


def fig_concordance(per_sample: pd.DataFrame, mode: str) -> None:
    sub = per_sample[per_sample["mode"] == mode]
    pivot = sub.pivot_table(index="HPFineNGroups", columns="sample", values="tp_rate")
    pivot = pivot[[s for s in SAMPLE_ORDER if s in pivot.columns]]
    corr = pivot.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=120)
    im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr.columns))); ax.set_xticklabels(corr.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(corr.index))); ax.set_yticklabels(corr.index, fontsize=9)
    for i in range(len(corr)):
        for j in range(len(corr)):
            v = corr.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 0.5 else "black")
    ax.set_title(f"Step 3 · Spearman concordance (TP-rate vs HPFineNGroups) · {mode}")
    plt.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout(); fig.savefig(OUT_FIG / f"fig03b_concordance_{mode}.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Step 4: stratified AUC

def step4_stratified_auc(df: pd.DataFrame, feature: str) -> pd.DataFrame:
    rows = []
    d = df.copy()
    d["CN_tier"] = cn_tier_F(d["Coverage_Multiple"])
    # Global
    a, lo, hi = auc_wilson(d["tp_label"].values, pd.to_numeric(d[feature], errors="coerce").values)
    rows.append(dict(stratum="Global", value="all", n=len(d), auc=a, auc_lo=lo, auc_hi=hi))
    # Per LOH
    if "LOH_Subtype" in d.columns:
        for v, g in d.groupby("LOH_Subtype"):
            a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
            rows.append(dict(stratum="LOH_Subtype", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    # Per AF class
    for v, g in d.groupby("AF_class"):
        a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
        rows.append(dict(stratum="AF_class", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    # Per CN tier
    for v, g in d.groupby("CN_tier", observed=True):
        a, lo, hi = auc_wilson(g["tp_label"].values, pd.to_numeric(g[feature], errors="coerce").values)
        rows.append(dict(stratum="CN_tier", value=str(v), n=len(g), auc=a, auc_lo=lo, auc_hi=hi))
    # Per mode / sample
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
    xerr = np.vstack([(d["auc"] - d["auc_lo"]).clip(lower=0).values, (d["auc_hi"] - d["auc"]).clip(lower=0).values])
    cmap = {"Global": "#1E2A44", "LOH_Subtype": "#A85540", "AF_class": "#009E73",
            "CN_tier": "#5B8DB7", "mode": "#D55E00", "sample": "#8A8A8A"}
    colors = [cmap.get(row["stratum"], "#444") for _, row in d.iterrows()]
    ax.barh(y_pos, d["auc"], xerr=xerr, color=colors, ecolor="#333", capsize=2)
    ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0.5, ls="--", c="#888"); ax.axvline(0.58, ls="--", c="#c00")
    ax.set_xlim(0.40, 0.85)
    ax.set_xlabel("AUC")
    ax.set_title(f"G6 Step 4 · Stratified AUC for {feature}")
    fig.tight_layout(); fig.savefig(OUT_FIG / f"fig04_stratified_auc_{feature}.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Step 5: confound guard (within-cell OLS on NumReads + vcf_AF + Coverage_Multiple)

def step5_confound(df: pd.DataFrame, features: List[str]) -> pd.DataFrame:
    from sklearn.linear_model import LinearRegression
    d = df.copy()
    d["vcf_AF"] = pd.to_numeric(d["vcf_AF"], errors="coerce")
    d["NumReads"] = pd.to_numeric(d["NumReads"], errors="coerce")
    d["Coverage_Multiple"] = pd.to_numeric(d["Coverage_Multiple"], errors="coerce")
    rows = []
    for feat in features:
        if feat not in d.columns:
            continue
        s = pd.to_numeric(d[feat], errors="coerce")
        mask = s.notna() & d["vcf_AF"].notna() & d["NumReads"].notna() & d["Coverage_Multiple"].notna()
        if mask.sum() < 200:
            continue
        X = d.loc[mask, ["NumReads", "vcf_AF", "Coverage_Multiple"]].values
        y = s.loc[mask].values
        reg = LinearRegression().fit(X, y)
        resid = y - reg.predict(X)
        a_raw, _, _ = auc_wilson(d.loc[mask, "tp_label"].values, y)
        a_res, _, _ = auc_wilson(d.loc[mask, "tp_label"].values, resid)
        # AF-bin cross
        af_results = []
        for af_cls in ["Extreme", "Intermediate", "Near-half"]:
            m2 = mask & (d["AF_class"] == af_cls)
            if m2.sum() < 100:
                continue
            a_sub, _, _ = auc_wilson(d.loc[m2, "tp_label"].values, s.loc[m2].values)
            af_results.append((af_cls, a_sub, int(m2.sum())))
        rows.append(dict(feature=feat, auc_raw=a_raw, auc_resid=a_res,
                         af_extreme=next((x[1] for x in af_results if x[0] == "Extreme"), np.nan),
                         af_inter=next((x[1] for x in af_results if x[0] == "Intermediate"), np.nan),
                         af_nearhalf=next((x[1] for x in af_results if x[0] == "Near-half"), np.nan)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "step5_confound.tsv", sep="\t", index=False)
    return out


def fig_step5(conf_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    d = conf_df.dropna(subset=["auc_raw"]).sort_values("auc_raw")
    y_pos = np.arange(len(d))
    ax.barh(y_pos - 0.2, d["auc_raw"], 0.4, color=PALETTE["accent"], label="raw AUC")
    ax.barh(y_pos + 0.2, d["auc_resid"], 0.4, color=PALETTE["grey"], label="residualized AUC (NumReads+vcf_AF+CovM)")
    ax.axvline(0.5, ls="--", c="#888"); ax.axvline(0.58, ls="--", c="#c00")
    ax.set_yticks(y_pos); ax.set_yticklabels(d["feature"], fontsize=9)
    ax.set_xlim(0.35, 0.85); ax.set_xlabel("AUC")
    ax.legend(loc="lower right")
    ax.set_title("G6 Step 5 · Confound guard (raw vs within-cell residualized)")
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig05_confound_residualized.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Step 6: spatial autocorrelation (chr + 5Mb bin)

def step6_spatial(df: pd.DataFrame, feature: str) -> pd.DataFrame:
    d = df[["Chr", "Pos", "tp_label", feature]].dropna().copy()
    d["bin"] = (d["Pos"] // 5_000_000).astype(int)
    grp = d.groupby(["Chr", "bin"])
    rows = []
    for (chrom, b), g in grp:
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
        ax.set_ylabel(f"per-bin AUC {feature}")
        plt.colorbar(sc, ax=ax, label="n variants per 5Mb bin")
        high_rate = sp[sp["tp_rate"] >= 0.8]
        if len(high_rate) > 0:
            ax.annotate(f"high-TP-rate bins: {len(high_rate)}, median AUC {high_rate['auc'].median():.3f}",
                        xy=(0.82, 0.72), xycoords="axes fraction", fontsize=8,
                        bbox=dict(boxstyle="round", facecolor="white", edgecolor="#888"))
    ax.set_title(f"G6 Step 6 · Spatial autocorrelation AUC vs TP rate · {feature}")
    fig.tight_layout(); fig.savefig(OUT_FIG / f"fig06_spatial_{feature}.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# HPFineD 6-pair analysis (G6-specific)

def hpfineD_analysis(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feat in HPFINE_D_PAIRS:
        if feat not in df.columns:
            continue
        # Global per-mode
        for mode in ["paired_full", "to_pileup"]:
            g = df[df["mode"] == mode]
            s = pd.to_numeric(g[feat], errors="coerce")
            a, lo, hi = auc_wilson(g["tp_label"].values, s.values)
            rows.append(dict(mode=mode, scope="global", feature=feat, n=len(g),
                             mean_tp=float(s[g["tp_label"] == 1].mean()),
                             mean_fp=float(s[g["tp_label"] == 0].mean()),
                             auc=a, auc_lo=lo, auc_hi=hi))
            # Per sample
            for smp, g2 in g.groupby("sample"):
                s2 = pd.to_numeric(g2[feat], errors="coerce")
                a2, lo2, hi2 = auc_wilson(g2["tp_label"].values, s2.values)
                rows.append(dict(mode=mode, scope=f"sample:{smp}", feature=feat, n=len(g2),
                                 mean_tp=float(s2[g2["tp_label"] == 1].mean()),
                                 mean_fp=float(s2[g2["tp_label"] == 0].mean()),
                                 auc=a2, auc_lo=lo2, auc_hi=hi2))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "hpfineD_auc.tsv", sep="\t", index=False)
    return out


def fig_hpfineD(hpd: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17, 9), dpi=120, sharey=True)
    feats = HPFINE_D_PAIRS
    for ax, feat in zip(axes.flatten(), feats):
        d = hpd[(hpd["feature"] == feat) & (hpd["scope"].str.startswith("sample:"))].copy()
        if d.empty:
            ax.set_visible(False); continue
        d["sample"] = d["scope"].str.replace("sample:", "", regex=False)
        modes = {"paired_full": ("#5B8DB7", -0.18), "to_pileup": ("#D55E00", 0.18)}
        x = np.arange(len(SAMPLE_ORDER))
        for mode, (c, offs) in modes.items():
            dm = d[d["mode"] == mode].set_index("sample").reindex(SAMPLE_ORDER)
            ax.bar(x + offs, dm["auc"].values, 0.35, color=c, label=mode, alpha=0.9)
            err = np.vstack([(dm["auc"] - dm["auc_lo"]).clip(lower=0).values,
                             (dm["auc_hi"] - dm["auc"]).clip(lower=0).values])
            ax.errorbar(x + offs, dm["auc"].values, yerr=err, fmt="none", ecolor="#333", lw=0.8, capsize=2)
        ax.axhline(0.5, ls="--", c="#888", lw=0.8); ax.axhline(0.58, ls="--", c="#c00", lw=0.8)
        ax.set_title(feat, fontsize=10)
        ax.set_xticks(x); ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right", fontsize=8)
        ax.set_ylim(0.35, 0.85)
    axes[0, 0].legend(loc="upper left", fontsize=8)
    axes[0, 0].set_ylabel("AUC"); axes[1, 0].set_ylabel("AUC")
    fig.suptitle("G6 · HPFineD pairwise-distance AUC across 7 samples × 2 modes", fontsize=13)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig07_hpfineD_per_sample_auc.png"); plt.close(fig)


def fig_hpfineD_heatmap(df: pd.DataFrame) -> None:
    # TP vs FP heatmap: rows = 6 pairs, cols = samples, value = mean
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=120)
    for ax, label, tp_val in zip(axes, ["TP only", "FP only"], [1, 0]):
        mat = np.full((len(HPFINE_D_PAIRS), len(SAMPLE_ORDER)), np.nan)
        for i, feat in enumerate(HPFINE_D_PAIRS):
            if feat not in df.columns:
                continue
            for j, smp in enumerate(SAMPLE_ORDER):
                g = df[(df["sample"] == smp) & (df["tp_label"] == tp_val)]
                v = pd.to_numeric(g[feat], errors="coerce").mean()
                mat[i, j] = v
        im = ax.imshow(mat, cmap="magma", aspect="auto")
        ax.set_yticks(range(len(HPFINE_D_PAIRS))); ax.set_yticklabels(HPFINE_D_PAIRS, fontsize=8)
        ax.set_xticks(range(len(SAMPLE_ORDER))); ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right", fontsize=8)
        ax.set_title(f"HPFineD mean · {label}")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isnan(v): continue
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="white")
        plt.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("G6 · HPFineD mean distance heatmap (per sample × 6 pairs)", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig07b_hpfineD_heatmap.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Fingerprint analysis: same_hap vs cross_het markers at NG=2

def fingerprint(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = df.copy()
    for c in HPFINE_N_COLS:
        d[c] = pd.to_numeric(d.get(c, 0), errors="coerce").fillna(0).astype(int)
    d["HPFineNGroups"] = pd.to_numeric(d["HPFineNGroups"], errors="coerce")
    ng2 = d[d["HPFineNGroups"] == 2].copy()
    # fingerprint: populated buckets (count > 0)
    ng2["pop"] = ng2[HPFINE_N_COLS].gt(0).sum(axis=1)
    ng2["same_hap_marker"] = (
        ((ng2["HPFineN_HP1"] > 0) & (ng2["HPFineN_HP1S"] > 0) & (ng2["HPFineN_HP2"] == 0) & (ng2["HPFineN_HP2S"] == 0))
        | ((ng2["HPFineN_HP2"] > 0) & (ng2["HPFineN_HP2S"] > 0) & (ng2["HPFineN_HP1"] == 0) & (ng2["HPFineN_HP1S"] == 0))
    )
    ng2["cross_het_marker"] = (
        ((ng2["HPFineN_HP1"] > 0) & (ng2["HPFineN_HP2S"] > 0) & (ng2["HPFineN_HP1S"] == 0) & (ng2["HPFineN_HP2"] == 0))
        | ((ng2["HPFineN_HP1S"] > 0) & (ng2["HPFineN_HP2"] > 0) & (ng2["HPFineN_HP1"] == 0) & (ng2["HPFineN_HP2S"] == 0))
    )
    ng2["other"] = ~(ng2["same_hap_marker"] | ng2["cross_het_marker"])
    # Per-sample × mode stats
    rows = []
    for (smp, mode), g in ng2.groupby(["sample", "mode"]):
        for flag in ["same_hap_marker", "cross_het_marker", "other"]:
            sub = g[g[flag]]
            if len(sub) == 0:
                continue
            rows.append(dict(
                sample=smp, mode=mode, fingerprint=flag,
                n=len(sub), n_tp=int(sub["tp_label"].sum()),
                n_fp=int((1 - sub["tp_label"]).sum()),
                tp_rate=float(sub["tp_label"].mean()),
            ))
    fp_df = pd.DataFrame(rows)
    fp_df.to_csv(OUT_DATA / "fingerprint_ng2.tsv", sep="\t", index=False)

    # Precision/recall of same_hap_marker as positive predictor
    pr_rows = []
    for (smp, mode), g in ng2.groupby(["sample", "mode"]):
        tp_in_marker = int(g.loc[g["same_hap_marker"], "tp_label"].sum())
        marker_total = int(g["same_hap_marker"].sum())
        total_tp_ng2 = int(g["tp_label"].sum())
        precision = tp_in_marker / marker_total if marker_total else np.nan
        recall = tp_in_marker / total_tp_ng2 if total_tp_ng2 else np.nan
        pr_rows.append(dict(
            sample=smp, mode=mode,
            n_ng2=len(g), n_marker=marker_total, n_tp_ng2=total_tp_ng2,
            tp_in_marker=tp_in_marker,
            precision=precision, recall=recall,
        ))
    pr_df = pd.DataFrame(pr_rows)
    pr_df.to_csv(OUT_DATA / "fingerprint_precision_recall.tsv", sep="\t", index=False)
    return fp_df, pr_df


def fig_fingerprint(fp_df: pd.DataFrame, pr_df: pd.DataFrame) -> None:
    # Plot A: TP rate by fingerprint × sample × mode
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=120, sharey=True)
    markers = ["same_hap_marker", "cross_het_marker", "other"]
    colors = {"same_hap_marker": "#009E73", "cross_het_marker": "#D55E00", "other": "#8A8A8A"}
    width = 0.25
    for ax, mode in zip(axes, ["paired_full", "to_pileup"]):
        sub = fp_df[fp_df["mode"] == mode]
        x = np.arange(len(SAMPLE_ORDER))
        for i, m in enumerate(markers):
            ss = sub[sub["fingerprint"] == m].set_index("sample").reindex(SAMPLE_ORDER)
            ax.bar(x + (i - 1) * width, ss["tp_rate"].values, width, color=colors[m], label=m, alpha=0.9)
            for xi, yi, ni in zip(x + (i - 1) * width, ss["tp_rate"].values, ss["n"].values):
                if not np.isnan(yi):
                    ax.text(xi, yi + 0.02, f"n={int(ni)}", ha="center", fontsize=6, rotation=90)
        ax.set_xticks(x); ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right", fontsize=8)
        ax.set_ylim(0, 1.05); ax.axhline(0.5, ls=":", c="#888")
        ax.set_title(f"NG=2 fingerprint TP rate · {mode}")
    axes[0].set_ylabel("TP rate"); axes[0].legend(loc="lower left", fontsize=8)
    fig.suptitle("G6 · NG=2 same_hap vs cross_het marker TP rate (Thread D replication)", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig08_fingerprint_tprate.png"); plt.close(fig)

    # Plot B: same_hap precision / recall scatter
    fig, ax = plt.subplots(figsize=(8, 6), dpi=120)
    for mode, marker in [("paired_full", "o"), ("to_pileup", "s")]:
        sub = pr_df[pr_df["mode"] == mode]
        ax.scatter(sub["recall"], sub["precision"], s=80, marker=marker, label=mode,
                   c=PALETTE["tp"] if mode == "paired_full" else PALETTE["fp"])
        for _, row in sub.iterrows():
            ax.annotate(row["sample"], (row["recall"], row["precision"]), fontsize=7,
                        xytext=(5, 3), textcoords="offset points")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.axhline(0.9, ls="--", c="#009E73", lw=1, label="precision ≥ 0.9")
    ax.set_xlabel("recall (same_hap TP / all TP in NG=2)")
    ax.set_ylabel("precision (same_hap TP / same_hap total)")
    ax.set_title("G6 · same_hap_marker precision-recall (7 samples × 2 modes)")
    ax.legend(loc="lower left")
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig08b_same_hap_precision_recall.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# CF vs NGroups collinearity

def cf_ngroups_collinearity(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    d = df.dropna(subset=["HPFineNGroups", "HPFine_NGroups_CF"]).copy()
    d["HPFineNGroups"] = d["HPFineNGroups"].astype(int)
    d["HPFine_NGroups_CF"] = d["HPFine_NGroups_CF"].astype(int)
    for (smp, mode), g in d.groupby(["sample", "mode"]):
        total = len(g)
        eq = int((g["HPFineNGroups"] == g["HPFine_NGroups_CF"]).sum())
        spearman = g[["HPFineNGroups", "HPFine_NGroups_CF"]].corr(method="spearman").iloc[0, 1]
        rows.append(dict(sample=smp, mode=mode, n=total, n_equal=eq, frac_equal=eq / total,
                         spearman=float(spearman)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "cf_ngroups_collinearity.tsv", sep="\t", index=False)
    return out


def fig_cf_vs_ngroups(df: pd.DataFrame) -> None:
    d = df.dropna(subset=["HPFineNGroups", "HPFine_NGroups_CF"]).copy()
    d["HPFineNGroups"] = d["HPFineNGroups"].astype(int)
    d["HPFine_NGroups_CF"] = d["HPFine_NGroups_CF"].astype(int)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=120)
    for ax, mode in zip(axes, ["paired_full", "to_pileup"]):
        sub = d[d["mode"] == mode]
        tab = pd.crosstab(sub["HPFine_NGroups_CF"], sub["HPFineNGroups"])
        im = ax.imshow(tab.values, cmap="Blues", aspect="auto")
        ax.set_xticks(range(tab.shape[1])); ax.set_xticklabels(tab.columns)
        ax.set_yticks(range(tab.shape[0])); ax.set_yticklabels(tab.index)
        ax.set_xlabel("HPFineNGroups (LabelTest path)")
        ax.set_ylabel("HPFine_NGroups_CF (Cluster-First)")
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                v = tab.values[i, j]
                if v > 0:
                    ax.text(j, i, f"{v}", ha="center", va="center", fontsize=8,
                            color="white" if v > tab.values.max() / 2 else "black")
        ax.set_title(f"{mode} · n={len(sub):,}")
        plt.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("G6 · CF (Cluster-First) vs NGroups (LabelTest) crosstab", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_FIG / "fig09_cf_vs_ngroups.png"); plt.close(fig)


# ----------------------------------------------------------------------------
# Main

def main():
    print("[G6] loading master ...")
    master = pd.read_csv(MASTER_TSV, sep="\t", low_memory=False)
    print(f"[G6] master shape: {master.shape}")
    print("[G6] enriching with raw G6 features ...")
    df = enrich_g6(master)
    print(f"[G6] enriched shape: {df.shape}")
    present = [c for c in G6_FEATURES if c in df.columns]
    print(f"[G6] G6 features present: {len(present)}/{len(G6_FEATURES)}: {present}")

    # Save enriched slice for provenance
    keep_cols = ["sample", "mode", "tp_label", "Chr", "Pos", "RegionID",
                 "LOH_Subtype", "AF_class", "vcf_AF", "AF", "NumReads",
                 "Coverage_Multiple", "Potential_LOH"] + present
    df[keep_cols].to_csv(OUT_DATA / "g6_enriched.tsv.gz", sep="\t", index=False, compression="gzip")

    # Step 1
    print("[G6] step 1: global distribution ...")
    s1 = step1_global(df)
    fig_step1(s1)
    fig_step1_ng_dist(df)

    # Step 2 for HPFineNGroups and HPFineF and HPFineD_HP1_HP2
    print("[G6] step 2: heatmaps ...")
    step2_heatmap(df, "HPFineNGroups", "ngroups")
    if "HPFineF" in df.columns:
        step2_heatmap(df, "HPFineF", "hpfinef")

    # Step 3
    print("[G6] step 3: per-sample consistency ...")
    ps = step3_per_sample(df)
    fig_step3(ps)
    for mode in ["paired_full", "to_pileup"]:
        fig_concordance(ps, mode)

    # Step 4
    print("[G6] step 4: stratified AUC ...")
    for feat in ["HPFineNGroups", "HPFineF", "HPFineD_HP1_HP2", "HPFineD_HP1_HP1S"]:
        if feat in df.columns:
            s4 = step4_stratified_auc(df, feat)
            fig_step4(s4, feat)

    # Step 5
    print("[G6] step 5: confound guard ...")
    conf_feats = [f for f in G6_FEATURES if f in df.columns and f not in ("HPFineSig",)]
    s5 = step5_confound(df, conf_feats)
    fig_step5(s5)

    # Step 6
    print("[G6] step 6: spatial autocorrelation ...")
    sp = step6_spatial(df, "HPFineNGroups")
    fig_step6(sp, "HPFineNGroups")

    # HPFineD analysis
    print("[G6] HPFineD per-pair analysis ...")
    hpd = hpfineD_analysis(df)
    fig_hpfineD(hpd)
    fig_hpfineD_heatmap(df)

    # Fingerprint (NG=2)
    print("[G6] NG=2 fingerprint ...")
    fp_df, pr_df = fingerprint(df)
    fig_fingerprint(fp_df, pr_df)

    # CF vs NGroups
    print("[G6] CF vs NGroups collinearity ...")
    coll = cf_ngroups_collinearity(df)
    print(coll.to_string())
    fig_cf_vs_ngroups(df)

    print(f"[G6] done — {len(list(OUT_FIG.glob('*.png')))} figures in {OUT_FIG}")
    print(f"[G6] data tsv in {OUT_DATA}")


if __name__ == "__main__":
    main()
