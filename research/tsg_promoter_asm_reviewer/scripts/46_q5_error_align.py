#!/usr/bin/env python3
"""
Q5 - Error / alignment confound + H-E clustering-vs-CN/error.

H-D: ASM signal NOT explained by error/alignment.
     Falsifier: residual(|db| | n_paired_cpg + median_cn) correlates with
                low-MAPQ / high-NM, rho>0.2 p<0.05.
     Threshold: rho<=0.1 -> technical noise not main driver.

H-E: read-level clustering blind_ari does NOT track CN/error.
     Falsifier: ari~median_cn OR ari~error significant (p<0.05, |rho|>0.2).
     null -> clustering reflects real methyl structure.

Method:
  H-D: rank-OLS residualize abs_delta on [n_paired_cpg, median_cn];
       Spearman residual vs {median_mapq, frac_mapq_lt20, median_nm,
       mean_nm_per_kb, frac_supp}.
  H-E: Spearman blind_ari vs {median_cn, median_mapq, median_nm};
       placebo_ari distribution + blind_ari vs placebo_ari.

Single sample (HCC1395). A pilot. Read-only on source; writes json + png.
"""
import json
import os
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from scipy import stats

# ----------------------------------------------------------------------
IN_TSV = (
    "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
    "genome_survey_v2/cn_confound/master_o2_error.tsv"
)
OUT_JSON = (
    "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
    "genome_survey_v2/cn_confound/q5_error_align.json"
)
OUT_PNG = (
    "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
    "figures/cn_confound/q5_error_align.png"
)
CJK_FONT = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"

# Decision thresholds (pre-registered in task)
RHO_FALSIFY = 0.2   # |rho|>this AND p<0.05 -> falsifier triggered
P_SIG = 0.05
RHO_NOISE_OK = 0.1  # rho<=this -> technical noise not main driver (H-D pass)

os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)


def rank(x):
    """Average-rank transform, NaN preserved."""
    x = np.asarray(x, float)
    out = np.full_like(x, np.nan, dtype=float)
    m = ~np.isnan(x)
    out[m] = stats.rankdata(x[m])
    return out


def spearman_safe(a, b):
    """Pairwise-complete Spearman. Returns dict; handles zero-variance / n<3."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = ~(np.isnan(a) | np.isnan(b))
    n = int(m.sum())
    if n < 3:
        return {"rho": None, "p": None, "n": n, "note": "n<3"}
    av, bv = a[m], b[m]
    if np.nanstd(av) == 0 or np.nanstd(bv) == 0:
        return {
            "rho": None,
            "p": None,
            "n": n,
            "note": "zero_variance_predictor",
        }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rho, p = stats.spearmanr(av, bv)
    return {"rho": float(rho), "p": float(p), "n": n, "note": "ok"}


def fmt(d):
    if d["rho"] is None:
        return f"rho=NA ({d['note']}, n={d['n']})"
    return f"rho={d['rho']:+.4f} p={d['p']:.3g} n={d['n']}"


# ----------------------------------------------------------------------
df = pd.read_csv(IN_TSV, sep="\t")
N = len(df)

# axis_type (file is all HP per inspection; derive robustly anyway)
df["axis_type"] = df["axis"].apply(lambda s: "HP" if "HP" in str(s) else "ALLELE")
axis_counts = df["axis_type"].value_counts().to_dict()

error_cols = ["median_mapq", "frac_mapq_lt20", "median_nm", "mean_nm_per_kb", "frac_supp"]
# Flag zero-variance predictors up front (e.g. median_mapq constant 60, frac_mapq_lt20 constant 0)
zero_var = {c: bool(df[c].nunique(dropna=True) <= 1) for c in error_cols}

# =====================================================================
# H-D: rank-OLS residualize abs_delta on [n_paired_cpg, median_cn]
# =====================================================================
y = rank(df["abs_delta"].values)
X = np.column_stack(
    [np.ones(N), rank(df["n_paired_cpg"].values), rank(df["median_cn"].values)]
)
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
yhat = X @ beta
resid = y - yhat  # residual abs_delta after removing coverage + CN (rank space)

# correlation of the rank-OLS fit (how much coverage+CN explains)
ss_tot = np.sum((y - y.mean()) ** 2)
ss_res = np.sum(resid**2)
r2_cov_cn = float(1 - ss_res / ss_tot)

hd = {}
for c in error_cols:
    # NOTE: expect NEGATIVE rho vs MAPQ if confound (low MAPQ -> inflated delta);
    # POSITIVE rho vs NM / frac_supp if alignment error drives the signal.
    hd[c] = spearman_safe(resid, df[c].values)

# H-D verdict: falsified if ANY live (non-zero-variance) error proxy has
# |rho|>0.2 & p<0.05. Pass (noise not main driver) if ALL live proxies |rho|<=0.1.
hd_live = {c: v for c, v in hd.items() if not zero_var[c] and v["rho"] is not None}
hd_falsifiers = [
    c
    for c, v in hd_live.items()
    if abs(v["rho"]) > RHO_FALSIFY and v["p"] < P_SIG
]
hd_all_below_noise = all(abs(v["rho"]) <= RHO_NOISE_OK for v in hd_live.values())
# also flag the "grey zone" 0.1<|rho|<=0.2 or p<0.05 with small effect
hd_grey = [
    c
    for c, v in hd_live.items()
    if (abs(v["rho"]) > RHO_NOISE_OK) and (c not in hd_falsifiers)
]

if hd_falsifiers:
    hd_verdict = "FALSIFIED_confound"  # error/alignment IS a driver
    hd_call = "REFUTE"  # refutes H-D (signal IS partly error)
elif hd_all_below_noise:
    hd_verdict = "PASS_noise_not_main_driver"
    hd_call = "CONFIRM"  # H-D holds: ASM not explained by error/alignment
else:
    hd_verdict = "INCONCLUSIVE_grey_zone"
    hd_call = "INCONCLUSIVE"

# =====================================================================
# H-E: blind_ari vs CN / error  + placebo distribution
# =====================================================================
he = {}
he["ari_vs_median_cn"] = spearman_safe(df["blind_ari"].values, df["median_cn"].values)
he["ari_vs_median_mapq"] = spearman_safe(df["blind_ari"].values, df["median_mapq"].values)
he["ari_vs_median_nm"] = spearman_safe(df["blind_ari"].values, df["median_nm"].values)
# extra live error proxies for completeness
he["ari_vs_mean_nm_per_kb"] = spearman_safe(
    df["blind_ari"].values, df["mean_nm_per_kb"].values
)
he["ari_vs_frac_supp"] = spearman_safe(df["blind_ari"].values, df["frac_supp"].values)

# placebo distribution + collider check
plc = df["placebo_ari"].dropna().values
n_plc = int(plc.size)
n_plc_nan = int(df["placebo_ari"].isna().sum())
placebo_stats = {
    "n_valid": n_plc,
    "n_nan": n_plc_nan,
    "min": float(np.min(plc)) if n_plc else None,
    "median": float(np.median(plc)) if n_plc else None,
    "mean": float(np.mean(plc)) if n_plc else None,
    "p90": float(np.percentile(plc, 90)) if n_plc else None,
    "max": float(np.max(plc)) if n_plc else None,
    "frac_gt_0.2": float(np.mean(plc > 0.2)) if n_plc else None,
}
# blind vs placebo: large placebo + correlated blind => collider/label-leak artifact
he["blind_vs_placebo_ari"] = spearman_safe(
    df["blind_ari"].values, df["placebo_ari"].values
)
# paired comparison blind vs placebo where both present (is real clustering > shuffled?)
both = df[["blind_ari", "placebo_ari"]].dropna()
if len(both) >= 3:
    w_stat, w_p = stats.wilcoxon(both["blind_ari"], both["placebo_ari"])
    blind_gt_placebo = {
        "n": int(len(both)),
        "median_blind": float(both["blind_ari"].median()),
        "median_placebo": float(both["placebo_ari"].median()),
        "wilcoxon_stat": float(w_stat),
        "wilcoxon_p": float(w_p),
        "median_diff_blind_minus_placebo": float(
            (both["blind_ari"] - both["placebo_ari"]).median()
        ),
    }
else:
    blind_gt_placebo = {"note": "n<3"}

# H-E verdict: tracks CN/error if ANY of {median_cn, the live error proxies} has
# |rho|>0.2 & p<0.05. null (clustering = real methyl structure) otherwise.
he_targets = [
    "ari_vs_median_cn",
    "ari_vs_median_nm",
    "ari_vs_mean_nm_per_kb",
    "ari_vs_frac_supp",
]  # ari_vs_median_mapq excluded (zero-variance predictor -> NA)
he_falsifiers = [
    k
    for k in he_targets
    if he[k]["rho"] is not None and abs(he[k]["rho"]) > RHO_FALSIFY and he[k]["p"] < P_SIG
]
# collider flag: high placebo (median>0.1 or frac>0.2 nontrivial) AND blind~placebo correlated
collider_flag = (
    placebo_stats["median"] is not None
    and placebo_stats["median"] > 0.1
    and he["blind_vs_placebo_ari"]["rho"] is not None
    and he["blind_vs_placebo_ari"]["rho"] > RHO_FALSIFY
    and he["blind_vs_placebo_ari"]["p"] < P_SIG
)

if he_falsifiers:
    he_verdict = "FALSIFIED_tracks_cn_or_error"
    he_call = "REFUTE"  # clustering DOES track CN/error -> not pure methyl structure
else:
    he_verdict = "NULL_clustering_independent_of_cn_error"
    he_call = "CONFIRM"  # H-E holds

# ----------------------------------------------------------------------
# JSON
# ----------------------------------------------------------------------
result = {
    "question": "Q5 - error/alignment confound (H-D) + clustering-vs-CN/error (H-E)",
    "sample": "HCC1395",
    "task_type": "A_pilot_single_sample",
    "n_rows": N,
    "axis_type_counts": axis_counts,
    "note_all_HP_axis": axis_counts.get("ALLELE", 0) == 0,
    "zero_variance_predictors": zero_var,
    "thresholds": {
        "rho_falsify": RHO_FALSIFY,
        "p_sig": P_SIG,
        "rho_noise_ok_HD": RHO_NOISE_OK,
    },
    "H_D": {
        "description": "residual(|db| | n_paired_cpg+median_cn) vs error/alignment proxies",
        "rank_ols_r2_coverage_plus_cn": r2_cov_cn,
        "spearman_resid_vs_error": hd,
        "live_proxies": list(hd_live.keys()),
        "falsifiers_rho_gt_0.2_p_lt_0.05": hd_falsifiers,
        "grey_zone_0.1_to_0.2": hd_grey,
        "all_live_below_noise_0.1": bool(hd_all_below_noise),
        "verdict": hd_verdict,
        "call": hd_call,
    },
    "H_E": {
        "description": "blind_ari vs CN/error + placebo distribution (collider check)",
        "spearman": he,
        "placebo_ari_distribution": placebo_stats,
        "blind_vs_placebo_paired_wilcoxon": blind_gt_placebo,
        "targets_tested": he_targets,
        "falsifiers_rho_gt_0.2_p_lt_0.05": he_falsifiers,
        "collider_artifact_flag": bool(collider_flag),
        "verdict": he_verdict,
        "call": he_call,
    },
    "outputs": {"json": OUT_JSON, "png": OUT_PNG},
    "caveats": [
        "Single sample HCC1395 -> A pilot; not cross-sample validated.",
        "median_mapq constant 60 & frac_mapq_lt20 constant 0 (zero variance): "
        "low-MAPQ confound impossible by construction, not testable here.",
        "placebo_ari has 258/596 NaN; placebo stats & blind-vs-placebo on complete pairs only.",
        "All rows are HP-axis (somatic-controlled); no ALLELE rows in this file.",
    ],
}

with open(OUT_JSON, "w") as f:
    json.dump(result, f, indent=2)

# ----------------------------------------------------------------------
# Figure: 5 panels
#   (1) |db| residual vs median_nm
#   (2) |db| residual vs mean_nm_per_kb (also annotate MAPQ degeneracy)
#   (3) blind_ari vs median_cn
#   (4) blind_ari vs median_mapq (will be a vertical line -> degenerate)
#   (5) placebo_ari hist
# ----------------------------------------------------------------------
# Labels are all ENGLISH (task: "prefer ENGLISH axis labels") -> use default
# DejaVu font (has Latin glyphs). CJK_FONT only registered as fallback, not forced.
fp = None
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.family"] = "DejaVu Sans"

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
ax = axes.ravel()


def scatter_panel(a, x, ydata, xlabel, ylabel, title, d, color="#2c6fbb"):
    xv = np.asarray(x, float)
    yv = np.asarray(ydata, float)
    m = ~(np.isnan(xv) | np.isnan(yv))
    a.scatter(xv[m], yv[m], s=14, alpha=0.5, color=color, edgecolors="none")
    a.set_xlabel(xlabel, fontsize=11)
    a.set_ylabel(ylabel, fontsize=11)
    if d["rho"] is not None:
        sub = f"Spearman rho={d['rho']:+.3f}  p={d['p']:.2g}  n={d['n']}"
        # fit line for visual
        if np.nanstd(xv[m]) > 0:
            try:
                z = np.polyfit(stats.rankdata(xv[m]), stats.rankdata(yv[m]), 1)
                xr = np.linspace(np.nanmin(xv[m]), np.nanmax(xv[m]), 50)
                # show trend in raw x via rank-linear on sorted order (illustrative)
            except Exception:
                pass
    else:
        sub = f"NA ({d['note']}, n={d['n']})"
    a.set_title(f"{title}\n{sub}", fontsize=10.5)
    a.grid(alpha=0.25)


# Panel 1: residual vs median_nm
scatter_panel(
    ax[0],
    df["median_nm"].values,
    resid,
    "median_nm (alignment edit distance)",
    "|db| residual (rank, | cov+CN)",
    "H-D: |db| residual vs median_nm",
    hd["median_nm"],
)
# Panel 2: residual vs mean_nm_per_kb
scatter_panel(
    ax[1],
    df["mean_nm_per_kb"].values,
    resid,
    "mean_nm_per_kb",
    "|db| residual (rank, | cov+CN)",
    "H-D: |db| residual vs mean_nm_per_kb",
    hd["mean_nm_per_kb"],
    color="#b5651d",
)
# Panel 3: residual vs frac_supp
scatter_panel(
    ax[2],
    df["frac_supp"].values,
    resid,
    "frac_supp (supplementary aln frac)",
    "|db| residual (rank, | cov+CN)",
    "H-D: |db| residual vs frac_supp",
    hd["frac_supp"],
    color="#7a4fa0",
)
# Panel 4: blind_ari vs median_cn (jitter cn for visibility)
cn = df["median_cn"].values.astype(float)
jit = cn + np.random.default_rng(0).normal(0, 0.06, size=cn.size)
ax[3].scatter(jit, df["blind_ari"].values, s=14, alpha=0.5, color="#2a9d8f", edgecolors="none")
ax[3].set_xlabel("median_cn (jittered)", fontsize=11)
ax[3].set_ylabel("blind_ari", fontsize=11)
d = he["ari_vs_median_cn"]
ax[3].set_title(
    f"H-E: blind_ari vs median_cn\nSpearman rho={d['rho']:+.3f} p={d['p']:.2g} n={d['n']}",
    fontsize=10.5,
)
ax[3].grid(alpha=0.25)
# Panel 5: blind_ari vs median_mapq (degenerate - annotate)
mapq = df["median_mapq"].values
ax[4].scatter(mapq, df["blind_ari"].values, s=14, alpha=0.5, color="#999999", edgecolors="none")
ax[4].set_xlabel("median_mapq", fontsize=11)
ax[4].set_ylabel("blind_ari", fontsize=11)
ax[4].set_title(
    "H-E: blind_ari vs median_mapq\n(median_mapq constant=60 -> rho NA, no MAPQ confound)",
    fontsize=10.5,
)
ax[4].grid(alpha=0.25)
ax[4].set_xlim(55, 65)
# Panel 6: placebo_ari histogram + blind_ari overlay
ax[5].hist(plc, bins=30, color="#e76f51", alpha=0.7, label=f"placebo_ari (n={n_plc})")
ax[5].axvline(np.median(plc), color="#b03020", ls="--", lw=1.5, label=f"placebo median={np.median(plc):.3f}")
ax[5].axvline(
    df["blind_ari"].median(),
    color="#264653",
    ls="-",
    lw=1.8,
    label=f"blind median={df['blind_ari'].median():.3f}",
)
ax[5].set_xlabel("ARI", fontsize=11)
ax[5].set_ylabel("count", fontsize=11)
ax[5].set_title("Placebo (label-shuffled) ARI distribution\n(low placebo -> blind ARI not collider artifact)", fontsize=10.5)
ax[5].legend(fontsize=8.5)
ax[5].grid(alpha=0.25)

suptitle = (
    f"Q5 Error/Alignment Confound - HCC1395 (n={N}, all HP-axis)  |  "
    f"H-D={hd_call} ({hd_verdict})   H-E={he_call} ({he_verdict})"
)
if fp:
    fig.suptitle(suptitle, fontsize=12.5, fontproperties=fp)
else:
    fig.suptitle(suptitle, fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_PNG, dpi=130)
plt.close(fig)

# ----------------------------------------------------------------------
print("=== H-D (error/alignment confound) ===")
print(f"rank-OLS R2 (coverage+CN explains |db|): {r2_cov_cn:.4f}")
for c in error_cols:
    zv = " [ZERO-VAR]" if zero_var[c] else ""
    print(f"  resid vs {c:16s}: {fmt(hd[c])}{zv}")
print(f"  falsifiers (|rho|>0.2,p<.05): {hd_falsifiers}")
print(f"  grey-zone (0.1<|rho|<=0.2): {hd_grey}")
print(f"  VERDICT: {hd_verdict} -> call={hd_call}")
print()
print("=== H-E (clustering vs CN/error) ===")
for k in ["ari_vs_median_cn", "ari_vs_median_mapq", "ari_vs_median_nm",
          "ari_vs_mean_nm_per_kb", "ari_vs_frac_supp", "blind_vs_placebo_ari"]:
    print(f"  {k:24s}: {fmt(he[k])}")
print(f"  placebo: median={placebo_stats['median']} p90={placebo_stats['p90']} "
      f"max={placebo_stats['max']} frac>0.2={placebo_stats['frac_gt_0.2']} "
      f"(n_valid={placebo_stats['n_valid']}, n_nan={placebo_stats['n_nan']})")
print(f"  blind vs placebo paired: {blind_gt_placebo}")
print(f"  falsifiers: {he_falsifiers}  collider_flag={collider_flag}")
print(f"  VERDICT: {he_verdict} -> call={he_call}")
print()
print(f"JSON -> {OUT_JSON}")
print(f"PNG  -> {OUT_PNG}")
