#!/usr/bin/env python3
"""Phase C · G9 — Allele-Specific Methylation (ASM) observation.

Group G9 covers 25 features centered on allele- and sample-level ASM plus HP
residual/direction-aware deltas.  This script executes Step 0-6 of the
feature_layered_observation SOP (02_methodology.md) and produces:

  figures/G9_asm/fig00..fig08*.png     (15+ panels)
  data/G9/G9_*.tsv                     per-step tables
  features/G9_asm.md                   companion 10-section narrative (written
                                       separately alongside this run)

Features in scope (25):
  Allele trio : AlleleDelta, AlleleP, AlleleSig
  Sample ASM  : SampleASM_Delta, SampleASM_P, SampleASM_Sig,
                SampleASM_NTumor, SampleASM_NNormal
  Normal base : NormalBaseline_Mean, NormalBaseline_Coverage
  HP residual : HP_Residual_Delta, HP_Residual_P, HP_Residual_Sig
  Tumor HP    : Tumor_HP_Delta, Tumor_HP_Valid, Tumor_HP1, Tumor_HP2,
                Tumor_HP_Signed_Delta
  Normal HP   : Normal_HP_Delta, Normal_HP_Valid, Normal_HP1, Normal_HP2,
                Normal_HP_Signed_Delta
  Combined    : HP_Signed_Residual, Combined_HP_Signed_Delta

Critical traps handled
~~~~~~~~~~~~~~~~~~~~~~
1. `AF` column in the master == `|AlleleDelta|` (L2 collider bias with
   tp_label — CLAUDE memory feedback_L2_collider_bias). We therefore
   compute raw AlleleDelta AUC and residualise on `vcf_AF` to confirm
   the collapse.
2. SampleASM_* / NormalBaseline_* / HP_Residual_* / Normal_HP_* are
   paired-only: archive TO significance_summary predates these columns.
   Step 4 stratified AUC reports per-mode coverage to make that explicit.
3. Normal_HP_Delta / Normal_HP_Signed_Delta / HP_Signed_Residual are
   currently 0% populated (even in paired_full canonical), because the
   ReadParser path that emits them requires a normal-HP-tagged stream
   that is not yet standard. We therefore treat Combined_HP_Signed_Delta
   == Tumor_HP_Signed_Delta (verified corr=1.0, diff=0 on HCC1395 TP).

Input : data/G9/master_g9.tsv.gz  (produced by g9_build_extended_master.py)
"""
from __future__ import annotations

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
DATA_TSV = ROOT / "data/G9/master_g9.tsv.gz"
FIG_DIR = ROOT / "figures/G9_asm"
DATA_DIR = ROOT / "data/G9"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]
LOH_ORDER = ["LOH_None", "LOH_Weak", "LOH_Noise", "LOH_Strong", "LOH_Subclone"]
AF_ORDER = ["Extreme", "Near-half", "Intermediate"]
CN_BREAKS = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
CN_LABELS = ["CN_Loss", "CN_Near1", "CN_Diploid", "CN_Gain", "CN_HighGain"]

CONT_FEATURES = [
    "AlleleDelta", "AlleleP",
    "SampleASM_Delta", "SampleASM_P",
    "SampleASM_NTumor", "SampleASM_NNormal",
    "NormalBaseline_Mean", "NormalBaseline_Coverage",
    "HP_Residual_Delta", "HP_Residual_P",
    "Tumor_HP_Delta", "Tumor_HP1", "Tumor_HP2", "Tumor_HP_Signed_Delta",
    "Normal_HP_Delta", "Normal_HP1", "Normal_HP2", "Normal_HP_Signed_Delta",
    "HP_Signed_Residual", "Combined_HP_Signed_Delta",
]
BIN_FEATURES = [
    "AlleleSig", "SampleASM_Sig", "HP_Residual_Sig",
    "Tumor_HP_Valid", "Normal_HP_Valid",
]

# Features for deep-dive (Steps 2,3,5,6)
FOCUS_CONT = [
    "AlleleDelta",
    "SampleASM_Delta",
    "NormalBaseline_Mean",
    "HP_Residual_Delta",
    "Tumor_HP_Signed_Delta",
    "Combined_HP_Signed_Delta",
]

# -------------------------------------------------------------------------- #
def auc_wilson_ci(y_true, scores):
    scores = np.asarray(scores, dtype=float)
    y_true = np.asarray(y_true)
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
    return (x.mean() - y.mean()) / sp if sp else 0.0


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


# -------------------------------------------------------------------------- #
def load() -> pd.DataFrame:
    print(f"[load] {DATA_TSV}", flush=True)
    df = pd.read_csv(DATA_TSV, sep="\t", low_memory=False)
    print(f"[load] rows={len(df):,}  cols={df.shape[1]}", flush=True)

    df["tp_label"] = pd.to_numeric(df["tp_label"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["tp_label"]).copy()
    df["tp_label"] = df["tp_label"].astype(int)

    df["cn_tier_F"] = pd.cut(df["Coverage_Multiple"].astype(float),
                             bins=CN_BREAKS, labels=CN_LABELS, right=True)
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("LOH_None")
    df["LOH_Subtype"] = pd.Categorical(df["LOH_Subtype"], categories=LOH_ORDER, ordered=True)
    df["AF_class"] = pd.Categorical(df["AF_class"].fillna("Unknown"),
                                    categories=AF_ORDER + ["Unknown"])

    # AlleleDelta is |mean_ALT - mean_REF|; the raw SS column is signed in
    # some schemas — coerce to numeric and expose its abs version for AUC use.
    for c in ["AlleleDelta"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
            df[c + "_abs"] = df[c].abs()
    # Tumor_HP_Delta and Normal_HP_Delta are unsigned in docs but the
    # canonical SS stores signed values — expose absolute helpers.
    for c in ["Tumor_HP_Delta", "Normal_HP_Delta", "HP_Residual_Delta",
              "SampleASM_Delta"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
            df[c + "_abs"] = df[c].abs()
    return df


# -------------------------------------------------------------------------- #
def step0_relations(df: pd.DataFrame):
    """Relation sanity: AlleleDelta vs vcf_AF + Combined vs Tumor_HP_Signed."""
    pairs_scatter = [
        ("AlleleDelta_abs", "vcf_AF"),
        ("AlleleDelta_abs", "AF"),
        ("SampleASM_Delta_abs", "AlleleDelta_abs"),
        ("HP_Residual_Delta_abs", "Tumor_HP_Delta_abs"),
        ("Combined_HP_Signed_Delta", "Tumor_HP_Signed_Delta"),
        ("NormalBaseline_Mean", "Tumor_HP_Delta_abs"),
    ]
    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    axes = axes.flatten()
    for idx, (a, b) in enumerate(pairs_scatter):
        ax = axes[idx]
        if a not in df.columns or b not in df.columns:
            ax.axis("off"); continue
        sub = df[[a, b]].dropna()
        if len(sub) < 50 or sub[a].std() == 0 or sub[b].std() == 0:
            ax.axis("off"); continue
        rho_p = stats.pearsonr(sub[a], sub[b])[0]
        rho_s = stats.spearmanr(sub[a], sub[b]).correlation
        sl = sub.sample(n=min(50_000, len(sub)), random_state=7)
        ax.scatter(sl[b], sl[a], s=1, alpha=0.05)
        try:
            k, c = np.polyfit(sl[b], sl[a], 1)
            xr = np.linspace(sl[b].min(), sl[b].max(), 20)
            ax.plot(xr, k*xr + c, c="red", lw=1, label=f"y={k:.2f}x+{c:.2f}")
        except Exception:
            pass
        ax.set_xlabel(b); ax.set_ylabel(a)
        ax.set_title(f"{a} vs {b}\nPearson={rho_p:.3f}  Spearman={rho_s:.3f}  n={len(sub):,}",
                     fontsize=9)
        ax.legend(fontsize=7, loc="best")
        rows.append({"a": a, "b": b, "pearson": rho_p, "spearman": rho_s,
                     "n": int(len(sub))})

    fig.suptitle("G9 · Step 0 Relation sanity — AF collider + Combined collapse check",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / "fig00_relation_sanity.png", dpi=140)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(DATA_DIR / "G9_relations.tsv", sep="\t", index=False)

    # Auxiliary: confirm Combined == Tumor_HP_Signed on non-NaN rows
    sub = df[["Combined_HP_Signed_Delta", "Tumor_HP_Signed_Delta",
              "Normal_HP_Signed_Delta"]].dropna(how="all")
    has_comb = df["Combined_HP_Signed_Delta"].notna()
    has_normsig = df["Normal_HP_Signed_Delta"].notna()
    has_hpresid = df["HP_Signed_Residual"].notna() if "HP_Signed_Residual" in df.columns else pd.Series(False, index=df.index)
    with open(DATA_DIR / "G9_combined_collapse_note.txt", "w") as f:
        f.write(f"Combined_HP_Signed_Delta non-NaN rows : {int(has_comb.sum()):,}\n")
        f.write(f"Normal_HP_Signed_Delta  non-NaN rows : {int(has_normsig.sum()):,}\n")
        f.write(f"HP_Signed_Residual      non-NaN rows : {int(has_hpresid.sum()):,}\n")
        common = df[["Combined_HP_Signed_Delta", "Tumor_HP_Signed_Delta"]].dropna()
        if len(common) >= 100:
            equal_frac = (common["Combined_HP_Signed_Delta"] ==
                          common["Tumor_HP_Signed_Delta"]).mean()
            diff = (common["Combined_HP_Signed_Delta"] -
                    common["Tumor_HP_Signed_Delta"]).abs().max()
            f.write(f"Identity test (Combined == Tumor_HP_Signed): "
                    f"equal_frac={equal_frac:.4f}  max_abs_diff={diff:.6e}\n")
            f.write("Interpretation: With Normal_HP_Signed_Delta = NaN, the "
                    "aggregate collapses to tumor-only. Combined provides no "
                    "independent signal in current canonical runs.\n")


# -------------------------------------------------------------------------- #
def step1_global_panel(df: pd.DataFrame):
    """Step 1 · global TP vs FP distribution panel for every continuous +
    binary feature. AUC + Wilson CI + Cohen d + Mann-Whitney."""
    feats = CONT_FEATURES + BIN_FEATURES
    # Replace absolute versions where appropriate for AUC (|delta| richer)
    rename_map = {"AlleleDelta": "AlleleDelta_abs",
                  "SampleASM_Delta": "SampleASM_Delta_abs",
                  "HP_Residual_Delta": "HP_Residual_Delta_abs",
                  "Tumor_HP_Delta": "Tumor_HP_Delta_abs",
                  "Normal_HP_Delta": "Normal_HP_Delta_abs"}

    n_rows = 5; n_cols = 5
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(26, 20))
    axes = axes.flatten()
    rows = []
    for i, feat in enumerate(feats):
        ax = axes[i]
        col = rename_map.get(feat, feat)
        if col not in df.columns:
            ax.text(0.5, 0.5, f"{feat}\nN/A", ha="center", va="center",
                    transform=ax.transAxes); ax.axis("off"); continue
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        y = df["tp_label"].to_numpy()
        tp = vals[y == 1]; fp = vals[y == 0]
        tp = tp[~np.isnan(tp)]; fp = fp[~np.isnan(fp)]
        if len(tp) < 10 or len(fp) < 10:
            ax.text(0.5, 0.5, f"{feat}\ninsufficient data\n(n_tp={len(tp)} n_fp={len(fp)})",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=8)
            ax.axis("off"); continue

        if feat in BIN_FEATURES:
            rate_tp = tp.mean(); rate_fp = fp.mean()
            ax.bar([0, 1], [rate_fp, rate_tp],
                   color=["#E07A5F", "#2E86AB"], edgecolor="black")
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"])
            ax.set_ylabel("rate of flag")
            ax.set_ylim(0, 1.0)
            auc, lo, hi, _, _ = auc_wilson_ci(y, vals)
            ax.set_title(f"{feat}\nrate_TP={rate_tp:.3f} rate_FP={rate_fp:.3f}\nAUC={auc:.3f}",
                         fontsize=9)
            rows.append({"feature": feat, "auc": auc, "lo": lo, "hi": hi,
                         "cohen_d": np.nan, "mwu_p": np.nan,
                         "mean_tp": rate_tp, "mean_fp": rate_fp,
                         "n_tp": int(len(tp)), "n_fp": int(len(fp))})
            continue

        log_scale = feat in {"AlleleP", "SampleASM_P", "HP_Residual_P"}
        p_tp = -np.log10(np.clip(tp, 1e-6, 1.0)) if log_scale else tp
        p_fp = -np.log10(np.clip(fp, 1e-6, 1.0)) if log_scale else fp
        parts = ax.violinplot([p_fp, p_tp], positions=[0, 1], widths=0.7,
                              showmeans=False, showmedians=True)
        for pc, color in zip(parts["bodies"], ["#E07A5F", "#2E86AB"]):
            pc.set_facecolor(color); pc.set_alpha(0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"FP\nn={len(fp):,}", f"TP\nn={len(tp):,}"], fontsize=8)
        auc, lo, hi, _, _ = auc_wilson_ci(y, vals)
        d = cohen_d(tp, fp); pv = mwu_p(tp, fp)
        ax.set_title(f"{feat}{' (-log10 p)' if log_scale else ''}\n"
                     f"AUC={auc:.3f} [{lo:.3f},{hi:.3f}] d={d:.2f}\n"
                     f"mw-p={pv:.1e}", fontsize=9)
        rows.append({"feature": feat, "auc": auc, "lo": lo, "hi": hi,
                     "cohen_d": d, "mwu_p": pv,
                     "mean_tp": float(np.nanmean(tp)),
                     "mean_fp": float(np.nanmean(fp)),
                     "n_tp": int(len(tp)), "n_fp": int(len(fp))})
    for j in range(len(feats), len(axes)):
        axes[j].axis("off")
    fig.suptitle("G9 ASM · Step 1 Global TP vs FP distribution "
                 "(7 samples × 2 modes pooled; SampleASM/NormalBaseline/HP_Residual are paired-only)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIG_DIR / "fig01_global_distribution.png", dpi=130)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(DATA_DIR / "G9_global_stats.tsv", sep="\t", index=False)
    return out


# -------------------------------------------------------------------------- #
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
    im1 = axes[1].imshow(piv_tp.where(mask_ok).values, aspect="auto", cmap="Blues")
    axes[1].set_title(f"B · median({feat}) | TP only")
    fig.colorbar(im1, ax=axes[1], fraction=0.03)
    im2 = axes[2].imshow(piv_fp.where(mask_ok).values, aspect="auto", cmap="Oranges")
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
    axes[3].set_title("D · Δ(TP−FP) median")
    fig.colorbar(im3, ax=axes[3], fraction=0.03)

    for ax in axes:
        ax.set_yticks(range(len(row_order)))
        ax.set_yticklabels(row_order, fontsize=7)
        ax.set_xticks(range(len(col_order)))
        ax.set_xticklabels(col_order, rotation=30, ha="right", fontsize=8)
    for i, r in enumerate(row_order):
        for j, c in enumerate(col_order):
            n = piv_n.loc[r, c] if (r in piv_n.index and c in piv_n.columns) else 0
            if pd.isna(n):
                n = 0
            if n >= min_n:
                axes[0].text(j, i, f"{int(n):,}", ha="center", va="center", fontsize=6)
    fig.suptitle(f"{feat} · G9 Step 2 LOH × AF × CN heatmap (pooled)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / f"fig02_{feat}_heatmap.png", dpi=130)
    plt.close(fig)

    flat = piv_delta.reset_index().melt(id_vars="row_key", var_name="cn_tier",
                                        value_name="delta_tp_fp")
    flat["feature"] = feat
    return flat


# -------------------------------------------------------------------------- #
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
        axes[i].imshow(piv.where(mask).values, aspect="auto",
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
    fig.suptitle(f"{feat} · G9 Step 3 Per-sample TP rate consistency (LOH × AF)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / f"fig03_{feat}_per_sample.png", dpi=130)
    plt.close(fig)


# -------------------------------------------------------------------------- #
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
    for sample in SAMPLE_ORDER:
        m = (df["sample"] == sample).to_numpy()
        if m.sum() >= 100:
            rows.append(("sample", sample, *auc_wilson_ci(y[m], s[m])))

    dfres = pd.DataFrame(rows, columns=["layer", "group", "auc", "lo", "hi",
                                        "n_pos", "n_neg"])
    dfres["feature"] = feat

    fig, ax = plt.subplots(figsize=(14, 6))
    colors = {"global": "#333", "LOH": "#2E86AB", "AF": "#F4A261",
              "CN": "#6A994E", "mode": "#8E44AD", "sample": "#444"}
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
    ax.set_title(f"{feat} · G9 Step 4 Stratified AUC")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"fig04_{feat}_stratified_auc.png", dpi=130)
    plt.close(fig)
    return dfres


# -------------------------------------------------------------------------- #
def step5_confound(df: pd.DataFrame, feat: str):
    """Residualise feat on vcf_AF + NumReads + Coverage_Multiple.

    For AlleleDelta the residualisation on vcf_AF is the critical L2 collider
    test: if raw AUC collapses to ~0.5 after residualisation, the signal is
    pure AF echo.
    """
    if feat not in df.columns:
        return None
    covars = [c for c in ["vcf_AF", "NumReads", "Coverage_Multiple"]
              if c in df.columns and c != feat]
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
            raw_af = auc_wilson_ci(y_tp[m], y_feat[m])[0]
            res_af = auc_wilson_ci(y_tp[m], resid[m])[0]
            af_rows.append((af, raw_af, res_af, int(m.sum())))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    bins = 50
    ax.hist(y_feat[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
            density=True, label="raw FP")
    ax.hist(y_feat[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
            density=True, label="raw TP")
    ax.set_title(f"raw {feat}  AUC={raw_auc:.3f}")
    ax.legend(fontsize=8)
    ax2 = axes[1]
    ax2.hist(resid[y_tp == 0], bins=bins, alpha=0.4, color="#E07A5F",
             density=True, label="resid FP")
    ax2.hist(resid[y_tp == 1], bins=bins, alpha=0.4, color="#2E86AB",
             density=True, label="resid TP")
    ax2.set_title(f"resid. on {covars}\nAUC={res_auc:.3f}\n"
                  "AF-bin (raw -> resid): " +
                  ", ".join([f"{a} {r:.3f}->{s:.3f}" for a, r, s, _ in af_rows]),
                  fontsize=9)
    ax2.legend(fontsize=8)
    fig.suptitle(f"{feat} · G9 Step 5 Confound guard (OLS residualisation on vcf_AF, NumReads, CovM)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / f"fig05_{feat}_confound.png", dpi=130)
    plt.close(fig)
    return {"feature": feat, "raw_auc": raw_auc, "resid_auc": res_auc,
            "collapse_delta": raw_auc - res_auc,
            "covars": ",".join(covars),
            "af_bin_raw": {a: r for a, r, _, _ in af_rows},
            "af_bin_resid": {a: s for a, _, s, _ in af_rows},
            "n": int(len(sub))}


# -------------------------------------------------------------------------- #
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
            flag = f"WARN artifact suspect (ΔAUC={hi_mean-lo_mean:+.3f})"
        else:
            flag = f"ok (ΔAUC={hi_mean-lo_mean:+.3f})"
    fig.suptitle(f"{feat} · G9 Step 6 Spatial auto-corr ({bin_mb} Mb) {flag}",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / f"fig06_{feat}_spatial.png", dpi=130)
    plt.close(fig)
    return {"feature": feat, "n_bins": len(s), "flag": flag}


# -------------------------------------------------------------------------- #
def step7_signed_vs_unsigned(df: pd.DataFrame):
    """Compare signed vs unsigned delta AUC per sample.

    If signed delta has noticeably different AUC from |delta|, then direction
    carries extra information (e.g. tumor hypomethylation vs hyper).
    """
    pairs = [
        ("Tumor_HP_Delta_abs", "Tumor_HP_Signed_Delta"),
        ("HP_Residual_Delta_abs", "HP_Signed_Residual"),
        ("SampleASM_Delta_abs", "SampleASM_Delta"),  # delta is signed natively
    ]
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    for ax, (u, s) in zip(axes, pairs):
        per_sample_u = []; per_sample_s = []; labels = []
        for sample in SAMPLE_ORDER:
            for mode in ["paired_full", "to_pileup"]:
                sub = df[(df["sample"] == sample) & (df["mode"] == mode)]
                if u not in sub.columns or s not in sub.columns:
                    continue
                y = sub["tp_label"].to_numpy()
                if y.sum() < 50 or (len(y) - y.sum()) < 50:
                    continue
                au = auc_wilson_ci(y, sub[u].to_numpy())[0]
                asg = auc_wilson_ci(y, sub[s].to_numpy())[0]
                if np.isnan(au) or np.isnan(asg):
                    continue
                per_sample_u.append(au); per_sample_s.append(asg)
                labels.append(f"{sample}/{mode[:3]}")
                rows.append({"feature_unsigned": u, "feature_signed": s,
                             "sample": sample, "mode": mode,
                             "auc_unsigned": au, "auc_signed": asg})
        if not labels:
            ax.axis("off"); continue
        x = np.arange(len(labels))
        ax.bar(x - 0.18, per_sample_u, width=0.35, color="#F4A261",
               edgecolor="black", label=u)
        ax.bar(x + 0.18, per_sample_s, width=0.35, color="#2E86AB",
               edgecolor="black", label=s)
        ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=75, fontsize=6)
        ax.set_ylim(0.3, 1.0)
        ax.set_ylabel("AUC")
        ax.set_title(f"{s} (signed) vs |unsigned|", fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle("G9 Step 7 Signed vs unsigned delta AUC per (sample,mode)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG_DIR / "fig07_signed_vs_unsigned.png", dpi=130)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(DATA_DIR / "G9_signed_vs_unsigned.tsv",
                              sep="\t", index=False)


# -------------------------------------------------------------------------- #
def step8_binary_per_class(df: pd.DataFrame):
    """Per-class TP rate for AlleleSig, SampleASM_Sig, HP_Residual_Sig."""
    rows = []
    for feat in ["AlleleSig", "SampleASM_Sig", "HP_Residual_Sig",
                 "Tumor_HP_Valid", "Normal_HP_Valid"]:
        if feat not in df.columns:
            continue
        sub = df.dropna(subset=[feat, "tp_label"]).copy()
        for stratum, mask in [("ALL", np.ones(len(sub), dtype=bool))] + \
                             [(l, (sub["LOH_Subtype"] == l).to_numpy())
                              for l in LOH_ORDER] + \
                             [(f"mode:{m}", (sub["mode"] == m).to_numpy())
                              for m in ["paired_full", "to_pileup"]]:
            if mask.sum() < 100:
                continue
            for cls in [0, 1]:
                m2 = mask & (sub[feat] == cls).to_numpy()
                if m2.sum() < 20:
                    continue
                y = sub.loc[m2, "tp_label"].to_numpy()
                k = int(y.sum()); n = int(len(y))
                p, lo, hi = wilson_ci_k_n(k, n)
                rows.append({"feature": feat, "stratum": stratum, "value": cls,
                             "tp_rate": p, "lo": lo, "hi": hi, "n": n})
    out = pd.DataFrame(rows)
    out.to_csv(DATA_DIR / "G9_binary_per_class.tsv", sep="\t", index=False)

    # Figure: ALL-stratum bars per feature + paired vs TO split
    fig, axes = plt.subplots(1, 5, figsize=(22, 5))
    feats = ["AlleleSig", "SampleASM_Sig", "HP_Residual_Sig",
             "Tumor_HP_Valid", "Normal_HP_Valid"]
    for ax, feat in zip(axes, feats):
        if feat not in out["feature"].unique():
            ax.axis("off"); ax.set_title(f"{feat}\nN/A"); continue
        sub = out[(out["feature"] == feat) &
                  (out["stratum"].isin(["ALL", "mode:paired_full", "mode:to_pileup"]))]
        strata = sub["stratum"].unique().tolist()
        x = np.arange(len(strata))
        for cls, color, offset in [(0, "#E07A5F", -0.18), (1, "#2E86AB", 0.18)]:
            sub_c = sub[sub["value"] == cls].set_index("stratum").reindex(strata)
            yerr = [sub_c["tp_rate"] - sub_c["lo"], sub_c["hi"] - sub_c["tp_rate"]]
            ax.bar(x + offset, sub_c["tp_rate"], width=0.35, color=color,
                   edgecolor="black", yerr=yerr, capsize=3, label=f"flag={cls}")
            for xi, v, n in zip(x + offset, sub_c["tp_rate"], sub_c["n"]):
                if pd.notna(v):
                    ax.text(xi, v + 0.01, f"{v:.2f}\nn={int(n)}",
                            ha="center", fontsize=6)
        ax.axhline(0.5, ls="--", c="gray", lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(strata, rotation=30, ha="right", fontsize=7)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("TP rate")
        ax.set_title(feat, fontsize=9)
        ax.legend(fontsize=6)
    fig.suptitle("G9 Step 8 Binary flag TP rate (Wilson 95% CI) — flag=1 means evidence-positive",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "fig08_binary_per_class.png", dpi=130)
    plt.close(fig)


# -------------------------------------------------------------------------- #
def step9_normbaseline_paired(df: pd.DataFrame):
    """Deep dive on NormalBaseline_Mean restricted to paired_full.

    Hypothesis: high NormalBaseline_Mean (normal-like methylation) = TP,
    low = FP (germline).
    """
    sub = df[(df["mode"] == "paired_full") &
             df["NormalBaseline_Mean"].notna() &
             df["tp_label"].notna()].copy()
    if len(sub) < 500:
        return
    y = sub["tp_label"].to_numpy()
    vals = sub["NormalBaseline_Mean"].to_numpy(dtype=float)
    auc_g, lo, hi, _, _ = auc_wilson_ci(y, vals)
    rows = [("global", auc_g, lo, hi, int(len(sub)))]
    for sample in SAMPLE_ORDER:
        s = sub[sub["sample"] == sample]
        if len(s) < 100:
            continue
        a, l, h, _, _ = auc_wilson_ci(s["tp_label"].to_numpy(),
                                      s["NormalBaseline_Mean"].to_numpy())
        rows.append((sample, a, l, h, int(len(s))))

    dfres = pd.DataFrame(rows, columns=["stratum", "auc", "lo", "hi", "n"])
    dfres.to_csv(DATA_DIR / "G9_normbaseline_paired.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    ax.hist(vals[y == 0], bins=50, alpha=0.5, color="#E07A5F",
            density=True, label=f"FP n={int((y==0).sum()):,}")
    ax.hist(vals[y == 1], bins=50, alpha=0.5, color="#2E86AB",
            density=True, label=f"TP n={int((y==1).sum()):,}")
    ax.set_xlabel("NormalBaseline_Mean (methylation 0-1)")
    ax.set_ylabel("density")
    ax.set_title(f"Global paired_full AUC={auc_g:.3f} [{lo:.3f},{hi:.3f}]")
    ax.legend(fontsize=9)

    ax = axes[1]
    x = np.arange(len(dfres))
    ax.bar(x, dfres["auc"], color="#6A994E", edgecolor="black",
           yerr=[dfres["auc"] - dfres["lo"], dfres["hi"] - dfres["auc"]],
           capsize=3)
    ax.axhline(0.5, ls="--", c="gray"); ax.axhline(0.58, ls=":", c="red")
    for xi, (_, row) in zip(x, dfres.iterrows()):
        ax.text(xi, row["auc"] + 0.005, f"{row['auc']:.3f}\nn={int(row['n']):,}",
                ha="center", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(dfres["stratum"], rotation=35, ha="right", fontsize=8)
    ax.set_ylim(0.3, 1.0)
    ax.set_title("NormalBaseline_Mean AUC per sample (paired_full only)")

    fig.suptitle("G9 Step 9 NormalBaseline_Mean paired-full deep dive",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "fig09_normbaseline_paired.png", dpi=130)
    plt.close(fig)


# -------------------------------------------------------------------------- #
def step10_collinearity(df: pd.DataFrame):
    """Spearman matrix across the G9 continuous features."""
    cols = [c for c in [
        "AlleleDelta_abs", "vcf_AF", "AF",
        "SampleASM_Delta_abs", "SampleASM_NTumor", "SampleASM_NNormal",
        "NormalBaseline_Mean", "NormalBaseline_Coverage",
        "HP_Residual_Delta_abs",
        "Tumor_HP_Delta_abs", "Tumor_HP_Signed_Delta",
        "Tumor_HP1", "Tumor_HP2",
        "Combined_HP_Signed_Delta",
        "NumReads", "Coverage_Multiple",
    ] if c in df.columns]
    sub = df[cols].dropna(how="all")
    if len(sub) < 1000:
        return
    sub = sub.sample(n=min(200_000, len(sub)), random_state=7)
    mat = sub.corr(method="spearman").astype(float)
    mat.to_csv(DATA_DIR / "G9_collinearity_spearman.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(13, 11))
    im = ax.imshow(mat.values, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(cols, fontsize=8)
    for i in range(len(cols)):
        for j in range(len(cols)):
            v = mat.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6, color="white" if abs(v) > 0.5 else "black")
    ax.set_title("G9 Step 10 Spearman collinearity matrix "
                 "(downsample 200k rows, all modes pooled)", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.03)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig10_collinearity.png", dpi=130)
    plt.close(fig)


# -------------------------------------------------------------------------- #
def main():
    df = load()

    print("[main] Step 0 relations", flush=True)
    step0_relations(df)

    print("[main] Step 1 global", flush=True)
    step1_global_panel(df)

    cell_deltas = []; auc_rows = []; confound_rows = []; spatial_rows = []
    for feat in FOCUS_CONT:
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
        print(f"[main] Step 5 confound {feat}", flush=True)
        c = step5_confound(df, feat)
        if c is not None:
            confound_rows.append(c)
        print(f"[main] Step 6 spatial {feat}", flush=True)
        sp = step6_spatial(df, feat)
        if sp is not None:
            spatial_rows.append(sp)

    print("[main] Step 7 signed vs unsigned", flush=True)
    step7_signed_vs_unsigned(df)

    print("[main] Step 8 binary per class", flush=True)
    step8_binary_per_class(df)

    print("[main] Step 9 NormalBaseline paired-full", flush=True)
    step9_normbaseline_paired(df)

    print("[main] Step 10 collinearity", flush=True)
    step10_collinearity(df)

    if auc_rows:
        pd.concat(auc_rows, ignore_index=True).to_csv(
            DATA_DIR / "G9_auc_table.tsv", sep="\t", index=False)
    if cell_deltas:
        pd.concat(cell_deltas, ignore_index=True).to_csv(
            DATA_DIR / "G9_cell_delta.tsv", sep="\t", index=False)
    if confound_rows:
        df_conf = pd.DataFrame(confound_rows)
        df_conf.to_csv(DATA_DIR / "G9_confound.tsv", sep="\t", index=False)
    if spatial_rows:
        pd.DataFrame(spatial_rows).to_csv(
            DATA_DIR / "G9_spatial.tsv", sep="\t", index=False)
    print("[main] done. Outputs:", flush=True)
    print(f"  figures: {FIG_DIR}")
    print(f"  tables : {DATA_DIR}")


if __name__ == "__main__":
    main()
