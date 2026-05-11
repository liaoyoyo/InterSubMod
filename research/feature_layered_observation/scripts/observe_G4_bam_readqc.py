#!/usr/bin/env python3
"""Phase C · G4 — BAM Read-Level QC feature group observation.

Scope (9 features from bam_readqc TSVs, extracted Phase B2):
  Mapping quality:  MapQ_mean, MapQ_median, MapQ_p10, MapQ_p90, LowMQ20_Frac
  Alignment error:  NM_mean (mismatch count per read)
  Soft clipping:    Softclip_Frac (fraction of reads with any softclip)
  Strand bias:      Strand_Bias (|frac_forward - 0.5| * 2)
  Read length:      Read_Length_mean

Data sources:
  14 TSVs at data/bam_readqc/bam_readqc_{sample}_{mode}.tsv
    (7 paired_full + 7 to_pileup; COLO829 to_pileup empty-header-only)
  Master: data/merged_with_vcf.tsv.gz (748,676 rows x 60 cols)
  Key:    sample + mode + Chr + Pos

Output:
  data/G4/master_g4.tsv.gz           extended master with 9 BAM QC cols
  data/G4/G4_global_stats.tsv        Step 1 global AUC table
  data/G4/G4_cell_delta.tsv          Step 2 LOH x AF x CN 32-cell delta
  data/G4/G4_auc_table.tsv           Step 4 stratified AUC
  data/G4/G4_confound.tsv            Step 5 residualized AUC
  data/G4/G4_spatial.tsv             Step 6 spatial summary
  data/G4/G4_persample_auc.tsv       Step 3 per-sample AUC
  data/G4/G4_corr_matrix.tsv         9x9 Spearman correlation
  figures/G4_bam_readqc/fig01..fig10 PNG panels
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
BAMQC_DIR = ROOT / "research/feature_layered_observation/data/bam_readqc"
MASTER_IN = ROOT / "research/feature_layered_observation/data/merged_with_vcf.tsv.gz"
OUT_DATA = ROOT / "research/feature_layered_observation/data/G4"
OUT_FIG = ROOT / "research/feature_layered_observation/figures/G4_bam_readqc"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)
MASTER_G4 = OUT_DATA / "master_g4.tsv.gz"

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]
MODES = ["paired_full", "to_pileup"]
LOH_ORDER = ["LOH_None", "LOH_Weak", "LOH_Noise", "LOH_Strong", "LOH_Subclone"]
AF_ORDER = ["Extreme", "Near-half", "Intermediate"]
CN_BREAKS = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
CN_LABELS = ["CN_Loss", "CN_Near1", "CN_Diploid", "CN_Gain", "CN_HighGain"]

G4_FEATURES = [
    "MapQ_mean", "MapQ_median", "MapQ_p10", "MapQ_p90",
    "LowMQ20_Frac", "NM_mean", "Softclip_Frac", "Strand_Bias",
    "Read_Length_mean",
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


# ----------------------------------------------------------------------------
# Master build

def build_master() -> pd.DataFrame:
    print(f"[load] master {MASTER_IN}", flush=True)
    master = pd.read_csv(MASTER_IN, sep="\t", low_memory=False)
    print(f"[load] master rows={len(master):,} cols={master.shape[1]}", flush=True)

    frames = []
    data_gaps = []
    for sample in SAMPLE_ORDER:
        for mode in MODES:
            p = BAMQC_DIR / f"bam_readqc_{sample}_{mode}.tsv"
            if not p.exists():
                data_gaps.append(f"{sample}/{mode} (missing)")
                print(f"[skip] {sample}/{mode}: file missing", flush=True); continue
            try:
                d = pd.read_csv(p, sep="\t", low_memory=False)
            except Exception as e:
                print(f"[err] {p}: {e}", flush=True); continue
            if len(d) == 0:
                data_gaps.append(f"{sample}/{mode} (empty)")
                print(f"[skip] {sample}/{mode}: empty", flush=True); continue
            # Normalise sample/mode
            d["sample"] = sample; d["mode"] = mode
            frames.append(d)
            print(f"[ok]   {sample}/{mode}: rows={len(d):,}", flush=True)

    if data_gaps:
        print(f"[warn] data gaps: {data_gaps}", flush=True)

    bam = pd.concat(frames, ignore_index=True)
    bam["Chr"] = bam["Chr"].astype(str)
    bam["Pos"] = pd.to_numeric(bam["Pos"], errors="coerce").astype("Int64")
    keep = ["sample", "mode", "Chr", "Pos"] + [c for c in G4_FEATURES if c in bam.columns]
    missing_feats = [c for c in G4_FEATURES if c not in bam.columns]
    if missing_feats:
        print(f"[warn] G4 features missing from TSV header: {missing_feats}", flush=True)
    bam = bam[keep].drop_duplicates(subset=["sample", "mode", "Chr", "Pos"], keep="first")
    print(f"[bam]  unique rows={len(bam):,}", flush=True)

    master["Chr"] = master["Chr"].astype(str)
    master["Pos"] = pd.to_numeric(master["Pos"], errors="coerce").astype("Int64")
    merged = master.merge(bam, on=["sample", "mode", "Chr", "Pos"], how="left",
                          suffixes=("", "_bam"))
    print(f"[merge] rows={len(merged):,}", flush=True)
    for c in G4_FEATURES:
        if c in merged.columns:
            n_ok = merged[c].notna().sum()
            pct = n_ok / len(merged) * 100
            print(f"[cov]  {c}: {n_ok:,}/{len(merged):,} ({pct:.1f}%)", flush=True)

    # Save compact extended master (master + G4 features only)
    merged.to_csv(MASTER_G4, sep="\t", index=False, compression="gzip")
    print(f"[save] {MASTER_G4}", flush=True)
    return merged


def load_master() -> pd.DataFrame:
    if MASTER_G4.exists():
        print(f"[load] {MASTER_G4}", flush=True)
        df = pd.read_csv(MASTER_G4, sep="\t", low_memory=False)
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
    return df


# ----------------------------------------------------------------------------
# Step 1: global distribution

def step1_global(df: pd.DataFrame) -> pd.DataFrame:
    feats = [f for f in G4_FEATURES if f in df.columns]
    rows = []
    n_cols = 3
    n_rows = int(np.ceil(len(feats) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.0*n_cols, 4.2*n_rows))
    axes = axes.flatten()
    for i, feat in enumerate(feats):
        ax = axes[i]
        vals = pd.to_numeric(df[feat], errors="coerce").to_numpy(float)
        y = df["tp_label"].to_numpy(int)
        tp = vals[y == 1]; fp = vals[y == 0]
        tp = tp[~np.isnan(tp)]; fp = fp[~np.isnan(fp)]
        if len(tp) < 10 or len(fp) < 10:
            ax.set_title(f"{feat}\n(insufficient n)", fontsize=9); ax.axis("off"); continue
        try:
            parts = ax.violinplot([fp, tp], positions=[0, 1], widths=0.7,
                                  showmeans=False, showmedians=True)
            for pc, color in zip(parts["bodies"], ["#E07A5F", "#2E86AB"]):
                pc.set_facecolor(color); pc.set_alpha(0.6)
        except Exception:
            pass
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"], fontsize=8)
        auc, lo, hi, _, _ = auc_wilson(y, vals)
        d = cohen_d(tp, fp); p = mwu_p(tp, fp)
        ax.set_title(f"{feat}\nAUC={auc:.3f}[{lo:.2f},{hi:.2f}] d={d:.2f}\n"
                     f"mean TP={np.nanmean(tp):.3g} FP={np.nanmean(fp):.3g}\n"
                     f"mwu-p={p:.1e}", fontsize=9)
        rows.append(dict(feature=feat, auc=auc, auc_lo=lo, auc_hi=hi,
                         cohen_d=d, mwu_p=p,
                         mean_tp=float(np.nanmean(tp)), mean_fp=float(np.nanmean(fp)),
                         median_tp=float(np.nanmedian(tp)),
                         median_fp=float(np.nanmedian(fp)),
                         std_tp=float(np.nanstd(tp)), std_fp=float(np.nanstd(fp)),
                         n_tp=len(tp), n_fp=len(fp)))
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle("G4 BAM Read-QC · Step 1 Global TP vs FP distribution "
                 "(7 samples x 2 modes pooled; COLO829 TO missing)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG / "fig01_global_distribution.png", dpi=130)
    plt.close(fig)
    out = pd.DataFrame(rows).sort_values("auc", ascending=False)
    out.to_csv(OUT_DATA / "G4_global_stats.tsv", sep="\t", index=False)

    # Global AUC bar panel
    fig, ax = plt.subplots(figsize=(10, 5))
    o = out.sort_values("auc", ascending=True)
    bars = ax.barh(range(len(o)), o["auc"], color="#2E86AB", edgecolor="black")
    ax.errorbar(o["auc"], range(len(o)),
                xerr=[o["auc"] - o["auc_lo"], o["auc_hi"] - o["auc"]],
                fmt="none", ecolor="black", capsize=3)
    ax.axvline(0.5, ls="--", c="gray", lw=0.7)
    ax.axvline(0.58, ls=":", c="red", lw=0.7, label="Beyond-AUC ceiling 0.58")
    ax.set_yticks(range(len(o))); ax.set_yticklabels(o["feature"])
    ax.set_xlim(0.3, 1.0)
    ax.set_xlabel("Global AUC (with 95% Wilson CI)")
    for i, (_, r) in enumerate(o.iterrows()):
        ax.text(r["auc"] + 0.005, i, f"{r['auc']:.3f}", va="center", fontsize=8)
    ax.set_title("G4 · 9 BAM QC features global AUC")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig02_global_auc_bar.png", dpi=130)
    plt.close(fig)
    return out


# ----------------------------------------------------------------------------
# Step 2: LOH x AF x CN 32-cell heatmaps per feature (4 panels)

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
    axes[0].set_title(f"A · TP rate (n>={min_n})")
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
    axes[3].set_title("D · delta(TP - FP) median")
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
    fig.suptitle(f"{feat} · G4 Step 2 LOH x AF x CN heatmap", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_FIG / f"fig03_{feat}_heatmap.png", dpi=130)
    plt.close(fig)
    flat = piv_delta.reset_index().melt(id_vars="row_key", var_name="cn_tier",
                                        value_name="delta_tp_fp")
    flat["feature"] = feat
    return flat


# ----------------------------------------------------------------------------
# Step 3: per-sample consistency (7 samples x 2 modes grid)

def step3_per_sample(df: pd.DataFrame, feats: List[str]) -> pd.DataFrame:
    rows = []
    samples = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    for s in samples:
        for mode in MODES:
            sub = df[(df["sample"] == s) & (df["mode"] == mode)]
            if len(sub) < 100:
                continue
            y = sub["tp_label"].to_numpy(int)
            for feat in feats:
                if feat not in sub.columns:
                    continue
                v = pd.to_numeric(sub[feat], errors="coerce").to_numpy(float)
                if np.sum(~np.isnan(v)) < 100:
                    continue
                auc, lo, hi, n1, n0 = auc_wilson(y, v)
                rows.append(dict(sample=s, mode=mode, feature=feat,
                                 auc=auc, auc_lo=lo, auc_hi=hi,
                                 n_tp=n1, n_fp=n0,
                                 tp_rate=n1/(n1+n0) if (n1+n0) else np.nan,
                                 n_valid=int(np.sum(~np.isnan(v)))))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G4_persample_auc.tsv", sep="\t", index=False)

    # Heatmap AUC x sample x mode for each feature
    mat = np.full((len(feats), len(samples) * 2), np.nan)
    col_labels = []
    for mi, mode in enumerate(MODES):
        for si, s in enumerate(samples):
            col_labels.append(f"{s}\n{mode.split('_')[0]}")
            for fi, feat in enumerate(feats):
                sel = out[(out["sample"] == s) & (out["mode"] == mode) & (out["feature"] == feat)]
                if not sel.empty:
                    mat[fi, mi * len(samples) + si] = sel["auc"].iloc[0]
    fig, ax = plt.subplots(figsize=(16, 6))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0.30, vmax=0.85)
    ax.set_yticks(range(len(feats))); ax.set_yticklabels(feats, fontsize=9)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not np.isnan(mat[i, j]):
                color = "black" if 0.40 <= mat[i, j] <= 0.70 else "white"
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                        fontsize=7, color=color)
    ax.axvline(len(samples) - 0.5, ls="--", c="black", lw=1.5)
    ax.set_title("G4 Step 3 · Per-sample AUC for 9 BAM QC features "
                 "(left: paired_full / right: to_pileup)", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.02)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig04_persample_auc_heatmap.png", dpi=130)
    plt.close(fig)

    # Cross-sample consistency stats
    cons_rows = []
    for feat in feats:
        sub = out[out["feature"] == feat]
        for mode in MODES:
            sm = sub[sub["mode"] == mode]
            if sm.empty:
                continue
            n_ok = (sm["auc"] >= 0.58).sum()
            n_samples = len(sm)
            cons_rows.append(dict(feature=feat, mode=mode,
                                  n_samples=n_samples, n_auc_gt_058=int(n_ok),
                                  median_auc=float(sm["auc"].median()),
                                  min_auc=float(sm["auc"].min()),
                                  max_auc=float(sm["auc"].max()),
                                  range_auc=float(sm["auc"].max() - sm["auc"].min())))
    pd.DataFrame(cons_rows).to_csv(OUT_DATA / "G4_cross_sample_consistency.tsv",
                                    sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
# Step 4: stratified AUC

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
        for mode in MODES:
            m = (df["mode"] == mode).to_numpy()
            if m.sum() >= 100:
                rows.append(("mode", mode, *auc_wilson(y[m], s[m])))
        dfres = pd.DataFrame(rows, columns=["layer", "group", "auc", "lo", "hi",
                                              "n_pos", "n_neg"])
        dfres["feature"] = feat
        all_rows.append(dfres)
    out = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    out.to_csv(OUT_DATA / "G4_auc_table.tsv", sep="\t", index=False)

    # Grid bar charts
    n_feats = len(feats)
    n_cols = 2
    n_rows = int(np.ceil(n_feats / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4.5 * n_rows))
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
                ax.text(j, row["auc"] + 0.01, f"{row['auc']:.2f}",
                        ha="center", fontsize=6)
        ax.axhline(0.5, ls="--", c="gray", lw=0.7)
        ax.axhline(0.58, ls=":", c="red", lw=0.7)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels([f"{r.layer}:{r.group}" for _, r in sub.iterrows()],
                           rotation=45, ha="right", fontsize=7)
        ax.set_ylim(0.3, 1.0)
        ax.set_title(f"{feat} stratified AUC", fontsize=10)
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle("G4 Step 4 · Stratified AUC (global / LOH / AF / CN / mode)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG / "fig05_stratified_auc.png", dpi=130)
    plt.close(fig)
    return out


# ----------------------------------------------------------------------------
# Step 5: confound guard (residualize on NumReads + vcf_AF + Coverage_Multiple + vcf_QUAL)

def step5_confound(df: pd.DataFrame, feats: List[str]) -> pd.DataFrame:
    rows = []
    covars_default = ["NumReads", "vcf_AF", "Coverage_Multiple", "vcf_QUAL"]
    n_feats = len(feats)
    fig, axes = plt.subplots(n_feats, 2, figsize=(14, 3.8 * n_feats))
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
        if pd.isna(res_auc):
            verdict = "NA"
        elif res_auc < 0.53:
            verdict = "COLLAPSED"
        elif abs((raw_auc or 0) - (res_auc or 0)) > 0.03:
            verdict = "ATTENUATED"
        else:
            verdict = "STABLE"
        rows.append(dict(feature=feat, raw_auc=raw_auc, resid_auc=res_auc,
                         delta=(res_auc - raw_auc) if (pd.notna(raw_auc) and pd.notna(res_auc)) else np.nan,
                         verdict=verdict, covars=",".join(covars), n_used=len(sub),
                         **{f"auc_{k}": v for k, v in af_bin.items()}))
        bins = 50
        ax = axes[fi, 0]
        ax.hist(y_feat[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
                density=True, label="FP")
        ax.hist(y_feat[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
                density=True, label="TP")
        ax.set_title(f"raw {feat} AUC={raw_auc:.3f}", fontsize=9)
        ax.legend(fontsize=8); ax.set_xlabel(feat)
        ax2 = axes[fi, 1]
        ax2.hist(resid[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
                 density=True, label="FP")
        ax2.hist(resid[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
                 density=True, label="TP")
        ax2.set_title(f"resid on {covars}\nAUC={res_auc:.3f} verdict={verdict}\n"
                      f"AF-bin AUC: " + ", ".join(f"{k}={v:.2f}" for k, v in af_bin.items()),
                      fontsize=8)
        ax2.legend(fontsize=8)
    fig.suptitle("G4 Step 5 · Confound guard (OLS residualize on "
                 "NumReads+vcf_AF+CovM+vcf_QUAL)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG / "fig06_confound_residualized.png", dpi=130)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G4_confound.tsv", sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
# Step 6: spatial autocorrelation

def step6_spatial(df: pd.DataFrame, feats: List[str], bin_mb: int = 5) -> pd.DataFrame:
    rows = []
    fig, axes = plt.subplots(int(np.ceil(len(feats) / 3)), 3,
                              figsize=(15, 4.2 * int(np.ceil(len(feats) / 3))))
    axes = axes.flatten()
    for fi, feat in enumerate(feats):
        ax = axes[fi]
        if feat not in df.columns:
            ax.axis("off"); continue
        d = df.dropna(subset=["Chr", "Pos", feat, "tp_label"]).copy()
        d["pos_bin"] = (d["Pos"].astype(int) // (bin_mb * 1_000_000))
        d["bin_key"] = d["Chr"].astype(str) + ":" + d["pos_bin"].astype(str)
        bin_rows = []
        for key, g in d.groupby("bin_key"):
            if len(g) < 50:
                continue
            y = g["tp_label"].to_numpy(int)
            if y.sum() < 5 or (len(y) - y.sum()) < 5:
                continue
            a, *_ = auc_wilson(y, pd.to_numeric(g[feat], errors="coerce").to_numpy(float))
            bin_rows.append(dict(bin=key, n=len(g), tp_rate=y.mean(), auc=a))
        if not bin_rows:
            ax.axis("off"); continue
        s = pd.DataFrame(bin_rows)
        ax.scatter(s["tp_rate"], s["auc"], alpha=0.3, s=8)
        ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
        ax.axvline(0.8, ls="--", c="purple", alpha=0.4)
        ax.set_xlabel("bin TP rate"); ax.set_ylabel(f"bin AUC")
        high_tp = s[s["tp_rate"] > 0.8]; low_tp = s[s["tp_rate"] < 0.5]
        flag = "ok"
        diff = np.nan
        if len(high_tp) >= 10 and len(low_tp) >= 10:
            diff = high_tp["auc"].mean() - low_tp["auc"].mean()
            flag = f"artifact-suspect (dAUC={diff:+.3f})" if diff > 0.08 \
                    else f"ok (dAUC={diff:+.3f})"
        ax.set_title(f"{feat}\n{flag}", fontsize=9)
        rows.append(dict(feature=feat, n_bins=len(s), flag=flag,
                         median_auc=float(s["auc"].median()),
                         dAUC_hi_lo=diff))
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"G4 Step 6 · Spatial ({bin_mb} Mb) per-bin AUC vs TP rate",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_FIG / "fig07_spatial.png", dpi=130)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA / "G4_spatial.tsv", sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
# Group-level visualisations

def fig_corr_matrix(df: pd.DataFrame, feats: List[str]):
    sub = df[feats].apply(pd.to_numeric, errors="coerce").dropna()
    if len(sub) < 100:
        return None
    corr = sub.corr(method="spearman")
    corr.to_csv(OUT_DATA / "G4_corr_matrix.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(feats))); ax.set_yticks(range(len(feats)))
    ax.set_xticklabels(feats, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(feats, fontsize=9)
    for i in range(len(feats)):
        for j in range(len(feats)):
            v = corr.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=7, color="white" if abs(v) > 0.5 else "black")
    fig.colorbar(im, ax=ax, fraction=0.04)
    ax.set_title("G4 · 9 BAM QC features · Spearman correlation")
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig08_corr_matrix.png", dpi=130)
    plt.close(fig)
    return corr


def fig_top_feature_persample(df: pd.DataFrame, top_feat: str):
    """Best feature: 7x2 sample/mode TP vs FP median bars."""
    if top_feat not in df.columns:
        return None
    rows = []
    samples = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    for s in samples:
        for mode in MODES:
            g = df[(df["sample"] == s) & (df["mode"] == mode)]
            if len(g) < 50:
                continue
            tp = pd.to_numeric(g.loc[g["tp_label"] == 1, top_feat], errors="coerce").dropna()
            fp = pd.to_numeric(g.loc[g["tp_label"] == 0, top_feat], errors="coerce").dropna()
            if len(tp) < 10 or len(fp) < 10:
                continue
            rows.append(dict(sample=s, mode=mode,
                             median_tp=float(tp.median()),
                             median_fp=float(fp.median()),
                             mean_tp=float(tp.mean()),
                             mean_fp=float(fp.mean()),
                             auc=auc_wilson(g["tp_label"].to_numpy(int),
                                             pd.to_numeric(g[top_feat], errors="coerce").to_numpy(float))[0],
                             n_tp=len(tp), n_fp=len(fp)))
    tbl = pd.DataFrame(rows)
    tbl.to_csv(OUT_DATA / f"G4_{top_feat}_persample.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    for mi, mode in enumerate(MODES):
        ax = axes[mi]
        sub = tbl[tbl["mode"] == mode]
        x = np.arange(len(sub))
        ax.bar(x - 0.2, sub["median_fp"], width=0.38, color="#E07A5F",
               edgecolor="black", label="median FP")
        ax.bar(x + 0.2, sub["median_tp"], width=0.38, color="#2E86AB",
               edgecolor="black", label="median TP")
        for xi, (_, r) in zip(x, sub.iterrows()):
            ax.text(xi, max(r["median_tp"], r["median_fp"]) * 1.02,
                    f"AUC={r['auc']:.2f}", ha="center", fontsize=7)
        ax.set_xticks(x); ax.set_xticklabels(sub["sample"], rotation=30, ha="right")
        ax.set_ylabel(top_feat); ax.legend()
        ax.set_title(f"{top_feat} per-sample ({mode})")
    fig.suptitle(f"G4 · {top_feat} per-sample TP vs FP medians", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG / f"fig09_{top_feat}_persample.png", dpi=130)
    plt.close(fig)


def fig_summary_vs_g3(df: pd.DataFrame, g4_feats: List[str]):
    """Bar chart comparing top G4 features AUC vs G3 vcf_QUAL paired baseline."""
    # Quick paired-mode AUC for each G4 feat + vcf_QUAL baseline
    sub = df[df["mode"] == "paired_full"]
    rows = []
    y = sub["tp_label"].to_numpy(int)
    for feat in g4_feats + ["vcf_QUAL"]:
        if feat not in sub.columns:
            continue
        v = pd.to_numeric(sub[feat], errors="coerce").to_numpy(float)
        if np.sum(~np.isnan(v)) < 100:
            continue
        a, lo, hi, _, _ = auc_wilson(y, v)
        rows.append(dict(feature=feat, auc=a, auc_lo=lo, auc_hi=hi))
    tbl = pd.DataFrame(rows).sort_values("auc", ascending=True)
    tbl.to_csv(OUT_DATA / "G4_vs_g3_paired.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#2E86AB" if f != "vcf_QUAL" else "#F4A261" for f in tbl["feature"]]
    ax.barh(range(len(tbl)), tbl["auc"], color=colors, edgecolor="black")
    ax.errorbar(tbl["auc"], range(len(tbl)),
                xerr=[tbl["auc"] - tbl["auc_lo"], tbl["auc_hi"] - tbl["auc"]],
                fmt="none", ecolor="black", capsize=3)
    ax.axvline(0.5, ls="--", c="gray", lw=0.7)
    ax.axvline(0.58, ls=":", c="red", lw=0.7)
    ax.axvline(0.813, ls="-.", c="purple", lw=0.7, label="G3 vcf_QUAL paired 0.813")
    ax.set_yticks(range(len(tbl))); ax.set_yticklabels(tbl["feature"])
    ax.set_xlabel("Paired-mode AUC")
    for i, (_, r) in enumerate(tbl.iterrows()):
        ax.text(r["auc"] + 0.005, i, f"{r['auc']:.3f}", va="center", fontsize=8)
    ax.set_title("G4 vs G3 vcf_QUAL paired_full AUC (reference line at 0.813)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig10_vs_g3_paired.png", dpi=130)
    plt.close(fig)
    return tbl


# ----------------------------------------------------------------------------
# Main

def main():
    df = load_master()
    print(f"[main] loaded master rows={len(df):,}", flush=True)
    feats = [f for f in G4_FEATURES if f in df.columns]
    print(f"[main] G4 features present: {feats}", flush=True)

    print("[main] Step 1 global distribution", flush=True)
    step1_out = step1_global(df)
    print(step1_out.to_string(index=False))

    # Focus features for heatmaps: top-4 by global AUC (or |AUC-0.5|)
    step1_rank = step1_out.assign(abs_dev=(step1_out["auc"] - 0.5).abs()) \
                          .sort_values("abs_dev", ascending=False)
    focus = step1_rank["feature"].head(5).tolist()
    print(f"[main] Step 2 heatmaps for {focus}", flush=True)
    cell_frames = []
    for feat in focus:
        r = step2_heatmap(df, feat)
        if r is not None:
            cell_frames.append(r)
    if cell_frames:
        pd.concat(cell_frames, ignore_index=True).to_csv(
            OUT_DATA / "G4_cell_delta.tsv", sep="\t", index=False)

    print("[main] Step 3 per-sample", flush=True)
    step3_per_sample(df, feats)

    print("[main] Step 4 stratified AUC", flush=True)
    step4_stratified(df, feats)

    print("[main] Step 5 confound guard", flush=True)
    step5_out = step5_confound(df, feats)
    print(step5_out.to_string(index=False))

    print("[main] Step 6 spatial", flush=True)
    step6_spatial(df, feats)

    print("[main] group viz · correlation matrix", flush=True)
    fig_corr_matrix(df, feats)

    top_feat = step1_out.iloc[0]["feature"]
    print(f"[main] group viz · per-sample for top feature {top_feat}", flush=True)
    fig_top_feature_persample(df, top_feat)

    print("[main] group viz · vs G3 paired", flush=True)
    fig_summary_vs_g3(df, feats)

    print("[main] done.")
    print(f"  figures: {OUT_FIG}")
    print(f"  tables : {OUT_DATA}")


if __name__ == "__main__":
    main()
