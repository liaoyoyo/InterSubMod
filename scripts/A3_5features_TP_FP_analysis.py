"""
A3_5features_TP_FP_analysis.py
================================
Phase A3 — PI 4-goal G3：5 features 分布 + 組合圖 + cross-sample direction consistency.

Outputs:
  phase2_completeness_audit/
    A3_per_feature_AUC_cohend_5sample.tsv
    A3_9cell_heatmap_data_5sample.tsv
    A3_direction_consistency_spearman.tsv
  figures/
    A3_5features_boxplot_per_sample.png
    A3_methyl_5sub_boxplot.png
    A3_AUC_cohend_bar_5sample.png
    A3_9cell_heatmap_5sample.png
    A3_direction_consistency_heatmap.png

Source TSVs (5 samples):
  HCC1395  -> research/methyl_augmented_filter_phase2/cycle2/data/per_bam_master_V6.tsv (bam_off_* prefix)
  HCC1937  -> research/methyl_augmented_filter_phase2/cycle2/data/HCC1937_master_augmented.tsv (V6_off_* prefix)
  HCC1954  -> .../HCC1954_master_augmented.tsv
  H1437    -> .../H1437_master_augmented.tsv
  H2009    -> .../H2009_master_augmented.tsv

Features (8 individual + composite 9-cell):
  caller_af, loh_inner_flag, Coverage_Multiple, NG,
  HPMergedDelta, HPFineF, NME_imbalance, Epipoly_Delta, ClusterPermanovaF

Stats:
  AUC = sklearn.metrics.roc_auc_score (TP=1, FP=0; nan-safe drop)
  Cohen's d = (mean_TP - mean_FP) / pooled_sd  (Cohen 1988)
  Mann-Whitney U two-sided p-value (scipy.stats.mannwhitneyu)
  Direction consistency = Spearman rho on per-sample Cohen's d sign-agreement.

9-cell binning:
  LOH:    inner (loh_inner_flag==1)  vs outer (==0)         -> 2 bin
  NG:     NG==2 / NG==3 / NG>=4                              -> 3 bin
  AF:     tertile on per-sample caller_af                   -> 3 bin (L/M/H)
  -> 2 x 3 x 3 = 18 conceptual cells; "9-cell" = LOH-averaged NG x AF (3x3 per LOH side).
     We report both: full 18 (in TSV) + 9-cell aggregate (heatmap = NG x AF, LOH split as
     facet).

CJK font: scripts.lib.plot_setup
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import roc_auc_score

# Bootstrap repo root + plot_setup
REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
sys.path.insert(0, str(REPO_ROOT))
from scripts.lib.plot_setup import setup_plot_style  # noqa: E402

setup_plot_style()
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

OUT_TSV = REPO_ROOT / "phase2_completeness_audit"
OUT_FIG = REPO_ROOT / "figures"
OUT_TSV.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Sample registry  (sample -> path, prefix)
# ---------------------------------------------------------------------------
DATA_DIR = REPO_ROOT / "research/methyl_augmented_filter_phase2/cycle2/data"

SAMPLES = [
    ("HCC1395", DATA_DIR / "per_bam_master_V6.tsv",        "bam_off_"),
    ("HCC1937", DATA_DIR / "HCC1937_master_augmented.tsv", "V6_off_"),
    ("HCC1954", DATA_DIR / "HCC1954_master_augmented.tsv", "V6_off_"),
    ("H1437",   DATA_DIR / "H1437_master_augmented.tsv",   "V6_off_"),
    ("H2009",   DATA_DIR / "H2009_master_augmented.tsv",   "V6_off_"),
]

# Canonical -> per-prefix column suffix
CANONICAL = {
    "caller_af":         ("caller_af", None),           # no prefix
    "loh_inner_flag":    ("loh_inner_flag", None),
    "Coverage_Multiple": ("Coverage_Multiple", None),
    "NG":                (None, "NG"),                  # bam_off_NG / V6_off_NG
    "HPMergedDelta":     (None, "meth_HPMergedDelta"),
    "HPFineF":           (None, "meth_HPFineF"),
    "NME_imbalance":     (None, "meth_NME_imbalance"),
    "Epipoly_Delta":     (None, "meth_Epipoly_Delta"),
    "ClusterPermanovaF": (None, "meth_ClusterPermanovaF"),
}

METHYL_FEATURES = ["HPMergedDelta", "HPFineF", "NME_imbalance", "Epipoly_Delta", "ClusterPermanovaF"]
CORE_FEATURES   = ["caller_af", "loh_inner_flag", "Coverage_Multiple", "NG"] + METHYL_FEATURES


def resolve_col(prefix: str, canonical: str) -> str:
    direct, suffix = CANONICAL[canonical]
    if direct is not None:
        return direct
    return f"{prefix}{suffix}"


def load_sample(sample: str, path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    sub_cols = {}
    for feat in CORE_FEATURES:
        col = resolve_col(prefix, feat)
        if col in df.columns:
            sub_cols[feat] = df[col]
        else:
            sub_cols[feat] = pd.Series([np.nan] * len(df))
    sub = pd.DataFrame(sub_cols)
    sub["label"] = df["label"].values if "label" in df.columns else (
        df["y"].map({1: "TP", 0: "FP"}).values if "y" in df.columns else None
    )
    sub["sample"] = sample
    # Ensure NG is numeric
    sub["NG"] = pd.to_numeric(sub["NG"], errors="coerce")
    sub["loh_inner_flag"] = pd.to_numeric(sub["loh_inner_flag"], errors="coerce").fillna(0).astype(int)
    return sub


print("=" * 70)
print("Phase A3 — loading 5 samples")
print("=" * 70)

frames = []
for sample, path, prefix in SAMPLES:
    print(f"  loading {sample} from {path.name} (prefix={prefix!r}) ...")
    sub = load_sample(sample, path, prefix)
    n_tp = int((sub["label"] == "TP").sum())
    n_fp = int((sub["label"] == "FP").sum())
    print(f"    -> rows={len(sub)}  TP={n_tp}  FP={n_fp}")
    frames.append(sub)

merged = pd.concat(frames, ignore_index=True)
print(f"\nMerged dataframe shape: {merged.shape}\n")


# ---------------------------------------------------------------------------
# Task 1 — Individual feature AUC / Cohen's d / Mann-Whitney p per sample
# ---------------------------------------------------------------------------
def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    n1, n2 = len(a), len(b)
    s1, s2 = a.var(ddof=1), b.var(ddof=1)
    pooled = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    if pooled == 0:
        return np.nan
    return (a.mean() - b.mean()) / pooled


def feature_stats(sub: pd.DataFrame, feat: str):
    tp = sub.loc[sub["label"] == "TP", feat].to_numpy(dtype=float)
    fp = sub.loc[sub["label"] == "FP", feat].to_numpy(dtype=float)
    tp = tp[~np.isnan(tp)]
    fp = fp[~np.isnan(fp)]
    if len(tp) < 5 or len(fp) < 5:
        return dict(AUC=np.nan, Cohens_d=np.nan, pvalue=np.nan,
                    n_TP=len(tp), n_FP=len(fp),
                    TP_median=np.nan, FP_median=np.nan)
    y = np.r_[np.ones_like(tp), np.zeros_like(fp)]
    x = np.r_[tp, fp]
    try:
        auc = roc_auc_score(y, x)
    except ValueError:
        auc = np.nan
    try:
        _, p = stats.mannwhitneyu(tp, fp, alternative="two-sided")
    except ValueError:
        p = np.nan
    return dict(
        AUC=auc,
        Cohens_d=cohen_d(tp, fp),
        pvalue=p,
        n_TP=len(tp),
        n_FP=len(fp),
        TP_median=float(np.median(tp)),
        FP_median=float(np.median(fp)),
    )


print("=" * 70)
print("Task 1 — per-feature AUC / Cohen's d / MWU per sample")
print("=" * 70)

rows = []
for sample, _, _ in SAMPLES:
    sub = merged[merged["sample"] == sample]
    for feat in CORE_FEATURES:
        st = feature_stats(sub, feat)
        st.update({"sample": sample, "feature": feat})
        rows.append(st)

per_feat_df = pd.DataFrame(rows)[
    ["sample", "feature", "AUC", "Cohens_d", "pvalue", "n_TP", "n_FP", "TP_median", "FP_median"]
]
out_path = OUT_TSV / "A3_per_feature_AUC_cohend_5sample.tsv"
per_feat_df.to_csv(out_path, sep="\t", index=False, float_format="%.4f")
print(f"  -> wrote {out_path.name}  rows={len(per_feat_df)}")

print("\n  AUC pivot (rows=feature, cols=sample):")
print(per_feat_df.pivot(index="feature", columns="sample", values="AUC").round(3))


# ---------------------------------------------------------------------------
# Task 2 — Cross-sample direction consistency (Spearman rho)
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Task 2 — direction consistency (Spearman) across 5 samples")
print("=" * 70)

# Per-feature: Spearman across samples on Cohen's d magnitude+sign.
# For pairwise concordance, compute rho of d-vector vs reference (HCC1395).
d_pivot = per_feat_df.pivot(index="feature", columns="sample", values="Cohens_d")

consistency_rows = []
samples_order = [s for s, _, _ in SAMPLES]
ref_sample = "HCC1395"

# (a) per-feature: sign-agreement count + range
for feat in CORE_FEATURES:
    ds = d_pivot.loc[feat, samples_order].astype(float)
    signs = np.sign(ds.fillna(0).values)
    n_pos = int((signs > 0).sum())
    n_neg = int((signs < 0).sum())
    n_zero = int((signs == 0).sum())
    sign_majority = max(n_pos, n_neg)
    consistency_rows.append({
        "feature": feat,
        "Cohens_d_HCC1395": ds.get(ref_sample, np.nan),
        "Cohens_d_HCC1937": ds.get("HCC1937", np.nan),
        "Cohens_d_HCC1954": ds.get("HCC1954", np.nan),
        "Cohens_d_H1437":   ds.get("H1437", np.nan),
        "Cohens_d_H2009":   ds.get("H2009", np.nan),
        "n_positive_d":     n_pos,
        "n_negative_d":     n_neg,
        "n_zero_d":         n_zero,
        "sign_majority":    sign_majority,
        "sign_consistent":  int(sign_majority >= 4),
        "d_range":          float(ds.max() - ds.min()) if ds.notna().any() else np.nan,
    })

# (b) pairwise Spearman rho on feature ranking across samples (samples = obs, feat = vars)
rank_df = d_pivot[samples_order].abs()  # use |d| for "predictive strength" ranking
ranked = rank_df.rank(axis=0, method="average")
pair_rho_rows = []
for i, s1 in enumerate(samples_order):
    for s2 in samples_order[i + 1:]:
        a = ranked[s1].dropna()
        b = ranked[s2].dropna()
        common = a.index.intersection(b.index)
        if len(common) < 3:
            rho, p = np.nan, np.nan
        else:
            rho, p = stats.spearmanr(ranked.loc[common, s1], ranked.loc[common, s2])
        pair_rho_rows.append({"sample_A": s1, "sample_B": s2, "spearman_rho_abs_d": rho, "pvalue": p})

cons_df = pd.DataFrame(consistency_rows)
pair_rho_df = pd.DataFrame(pair_rho_rows)

out_path = OUT_TSV / "A3_direction_consistency_spearman.tsv"
with open(out_path, "w") as fh:
    fh.write("# Section 1: per-feature Cohen's d across 5 samples + sign-agreement\n")
    cons_df.to_csv(fh, sep="\t", index=False, float_format="%.4f")
    fh.write("\n# Section 2: pairwise Spearman rho on |Cohen's d| ranking across features\n")
    pair_rho_df.to_csv(fh, sep="\t", index=False, float_format="%.4f")
print(f"  -> wrote {out_path.name}")
print("\n  Sign-consistency (>=4/5 same direction):")
print(cons_df[["feature", "n_positive_d", "n_negative_d", "sign_consistent"]].to_string(index=False))
print("\n  Pairwise Spearman rho (|d| ranking):")
print(pair_rho_df.to_string(index=False))


# ---------------------------------------------------------------------------
# Task 3 — 9-cell heatmap (LOH x NG x AF) per sample + aggregate
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Task 3 — 9-cell (LOH x NG x AF) heatmap")
print("=" * 70)


def assign_bins(sub: pd.DataFrame) -> pd.DataFrame:
    sub = sub.copy()
    # LOH: inner=1 / outer=0
    sub["LOH_bin"] = sub["loh_inner_flag"].map({1: "inner", 0: "outer"})
    # NG: 2 / 3 / >=4
    def ng_bin(v):
        if pd.isna(v):
            return np.nan
        v = int(v)
        if v <= 2:
            return "NG2"
        if v == 3:
            return "NG3"
        return "NGge4"
    sub["NG_bin"] = sub["NG"].map(ng_bin)
    # AF: tertile per sample
    af = sub["caller_af"].astype(float)
    try:
        af_q = pd.qcut(af, q=3, labels=["AF_L", "AF_M", "AF_H"], duplicates="drop")
    except ValueError:
        af_q = pd.Series(np.nan, index=sub.index)
    sub["AF_bin"] = af_q.astype("object")
    return sub


cell_rows = []
for sample, _, _ in SAMPLES:
    sub = assign_bins(merged[merged["sample"] == sample])
    grp = sub.groupby(["LOH_bin", "NG_bin", "AF_bin"], dropna=True)["label"]
    for (loh, ng, af), labels in grp:
        if pd.isna(loh) or pd.isna(ng) or pd.isna(af):
            continue
        n_tp = int((labels == "TP").sum())
        n_fp = int((labels == "FP").sum())
        n = n_tp + n_fp
        cell_rows.append({
            "sample": sample,
            "LOH_bin": loh,
            "NG_bin": ng,
            "AF_bin": af,
            "n_TP": n_tp,
            "n_FP": n_fp,
            "n_total": n,
            "TP_rate": n_tp / n if n > 0 else np.nan,
        })

cell_df = pd.DataFrame(cell_rows)
out_path = OUT_TSV / "A3_9cell_heatmap_data_5sample.tsv"
cell_df.to_csv(out_path, sep="\t", index=False, float_format="%.4f")
print(f"  -> wrote {out_path.name}  rows={len(cell_df)}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Generating figures")
print("=" * 70)

NG_ORDER  = ["NG2", "NG3", "NGge4"]
AF_ORDER  = ["AF_L", "AF_M", "AF_H"]
LOH_ORDER = ["outer", "inner"]
SAMPLE_ORDER = [s for s, _, _ in SAMPLES]

# ---- Fig 1: 4 core features boxplot — 4 rows feature × 5 cols sample
def _safe_log(x, eps=1e-3):
    x = np.asarray(x, dtype=float)
    x[x <= 0] = eps
    return np.log10(x)


def boxplot_grid(features, fname, log_y_for=None):
    log_y_for = log_y_for or set()
    n_feat = len(features)
    n_sample = len(SAMPLE_ORDER)
    fig, axes = plt.subplots(n_feat, n_sample, figsize=(3.0 * n_sample, 2.6 * n_feat),
                             squeeze=False)
    for i, feat in enumerate(features):
        for j, sample in enumerate(SAMPLE_ORDER):
            ax = axes[i][j]
            sub = merged[merged["sample"] == sample][[feat, "label"]].dropna()
            if len(sub) == 0:
                ax.set_axis_off()
                continue
            if feat in log_y_for:
                sub = sub[sub[feat] > 0]
                sub = sub.assign(**{feat: np.log10(sub[feat])})
            sns.boxplot(data=sub, x="label", y=feat, order=["FP", "TP"],
                        palette={"FP": "#d62728", "TP": "#1f77b4"},
                        showfliers=False, ax=ax)
            # AUC annotation
            st = feature_stats(merged[merged["sample"] == sample], feat)
            ax.set_title(f"{sample}\nAUC={st['AUC']:.3f}  d={st['Cohens_d']:.2f}",
                         fontsize=9)
            ax.set_xlabel("")
            if j == 0:
                ax.set_ylabel(feat + (" (log10)" if feat in log_y_for else ""), fontsize=9)
            else:
                ax.set_ylabel("")
    fig.suptitle(f"Phase A3 — TP vs FP boxplot per sample × feature", y=1.0, fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_FIG / fname, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> wrote {fname}")


boxplot_grid(
    ["caller_af", "loh_inner_flag", "Coverage_Multiple", "NG"],
    "A3_5features_boxplot_per_sample.png",
)
boxplot_grid(
    METHYL_FEATURES,
    "A3_methyl_5sub_boxplot.png",
)


# ---- Fig 3: AUC + Cohen's d bar chart per sample
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

auc_long = per_feat_df[["sample", "feature", "AUC"]].copy()
sns.barplot(data=auc_long, x="feature", y="AUC", hue="sample",
            order=CORE_FEATURES, hue_order=SAMPLE_ORDER, ax=axes[0])
axes[0].axhline(0.5, color="grey", linestyle="--", linewidth=0.8)
axes[0].axhline(0.58, color="#d62728", linestyle=":", linewidth=0.8,
                label="AUC=0.58 ceiling")
axes[0].set_title("AUC per feature × sample (TP=1 vs FP=0)")
axes[0].set_ylim(0.3, 1.0)
axes[0].set_xlabel("")
axes[0].legend(loc="lower right", fontsize=8, ncol=3)

d_long = per_feat_df[["sample", "feature", "Cohens_d"]].copy()
sns.barplot(data=d_long, x="feature", y="Cohens_d", hue="sample",
            order=CORE_FEATURES, hue_order=SAMPLE_ORDER, ax=axes[1], legend=False)
axes[1].axhline(0, color="grey", linestyle="--", linewidth=0.8)
axes[1].set_title("Cohen's d per feature × sample  (positive → TP>FP)")
axes[1].set_xlabel("Feature")
plt.setp(axes[1].get_xticklabels(), rotation=20, ha="right")

fig.tight_layout()
fig.savefig(OUT_FIG / "A3_AUC_cohend_bar_5sample.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("  -> wrote A3_AUC_cohend_bar_5sample.png")


# ---- Fig 4: 9-cell heatmap — 5 samples × 2 LOH facets
fig, axes = plt.subplots(len(SAMPLE_ORDER), 2,
                         figsize=(7.5, 2.6 * len(SAMPLE_ORDER)),
                         squeeze=False)
for i, sample in enumerate(SAMPLE_ORDER):
    for j, loh in enumerate(LOH_ORDER):
        ax = axes[i][j]
        sub = cell_df[(cell_df["sample"] == sample) & (cell_df["LOH_bin"] == loh)]
        if len(sub) == 0:
            ax.set_axis_off()
            continue
        pivot = sub.pivot(index="NG_bin", columns="AF_bin", values="TP_rate")
        pivot = pivot.reindex(index=NG_ORDER, columns=AF_ORDER)
        cnt   = sub.pivot(index="NG_bin", columns="AF_bin", values="n_total").reindex(
            index=NG_ORDER, columns=AF_ORDER)
        annot = pivot.copy().astype(object)
        for r in NG_ORDER:
            for c in AF_ORDER:
                v = pivot.at[r, c] if (r in pivot.index and c in pivot.columns) else np.nan
                n = cnt.at[r, c]   if (r in cnt.index   and c in cnt.columns)   else np.nan
                annot.at[r, c] = f"{v:.2f}\n(n={int(n)})" if not pd.isna(v) and not pd.isna(n) else "—"
        sns.heatmap(pivot, vmin=0.0, vmax=1.0, cmap="RdYlBu",
                    annot=annot.values, fmt="", cbar=(j == 1),
                    ax=ax, linewidths=0.4, linecolor="white",
                    annot_kws={"fontsize": 7})
        ax.set_title(f"{sample}  LOH={loh}", fontsize=9)
        ax.set_xlabel("AF tertile" if i == len(SAMPLE_ORDER) - 1 else "")
        ax.set_ylabel("NG bin" if j == 0 else "")
fig.suptitle("Phase A3 — 9-cell TP-rate heatmap (NG × AF × LOH) per sample", y=1.0, fontsize=12)
fig.tight_layout()
fig.savefig(OUT_FIG / "A3_9cell_heatmap_5sample.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("  -> wrote A3_9cell_heatmap_5sample.png")


# ---- Fig 5: direction consistency heatmap (Cohen's d sign, 9 feature × 5 sample)
fig, ax = plt.subplots(figsize=(7, 5))
d_mat = d_pivot.reindex(index=CORE_FEATURES, columns=SAMPLE_ORDER).astype(float)
sns.heatmap(d_mat, annot=True, fmt=".2f", center=0, cmap="RdBu_r",
            vmin=-1.5, vmax=1.5, cbar_kws={"label": "Cohen's d"},
            linewidths=0.5, linecolor="white", ax=ax)
ax.set_title("Phase A3 — Cohen's d (TP - FP) heatmap : 9 features × 5 samples")
ax.set_xlabel("Sample")
ax.set_ylabel("Feature")
fig.tight_layout()
fig.savefig(OUT_FIG / "A3_direction_consistency_heatmap.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("  -> wrote A3_direction_consistency_heatmap.png")


# ---------------------------------------------------------------------------
# Status summary (< 300 chars)
# ---------------------------------------------------------------------------
auc_min = float(np.nanmin(per_feat_df["AUC"]))
auc_max = float(np.nanmax(per_feat_df["AUC"]))
auc_med = float(np.nanmedian(per_feat_df["AUC"]))
consistent_feats = cons_df.loc[cons_df["sign_consistent"] == 1, "feature"].tolist()
rho_med = float(np.nanmedian(pair_rho_df["spearman_rho_abs_d"]))
rho_min = float(np.nanmin(pair_rho_df["spearman_rho_abs_d"]))
rho_max = float(np.nanmax(pair_rho_df["spearman_rho_abs_d"]))
tpr_min = float(np.nanmin(cell_df["TP_rate"]))
tpr_max = float(np.nanmax(cell_df["TP_rate"]))

print("\n" + "=" * 70)
print("STATUS SUMMARY")
print("=" * 70)
print(
    f"5 features × 5 samples: AUC range {auc_min:.3f}-{auc_max:.3f} (median {auc_med:.3f}). "
    f"Sign-consistent (>=4/5 samples agree): {len(consistent_feats)}/9 -> {consistent_feats}. "
    f"Pairwise Spearman rho on |Cohen's d| ranking: {rho_min:.2f}-{rho_max:.2f} (median {rho_med:.2f}). "
    f"9-cell TP-rate range across LOH×NG×AF×5 sample = {tpr_min:.2f}-{tpr_max:.2f}. "
    f"caveat: caller_af direction differs across samples (HCC1395 vs HCC1937/1954 tumour purity)."
)
