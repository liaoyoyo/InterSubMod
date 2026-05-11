#!/usr/bin/env python3
"""R1 residualization + visualization.

For HCC1395 Normal BAM pilot (TO + paired), remove known confounds and
re-compute AUC with bootstrap CI. Produce 6 figures.

Residualizations:
  - NormalBaseline_Coverage ~ NumReads + Coverage_Multiple (OLS residual)
  - Fisher_Frac_Sig ~ Fisher_N_Tested + NumReads + CovM (OLS residual)
  - SampleASM_Delta ~ NumReads + Coverage_Multiple

Output figures:
  1. feature_auc_bar_TO.png, feature_auc_bar_paired.png
  2. feature_residualized_auc_bar.png (two-panel)
  3. top_features_tp_fp_distribution_TO.png
  4. top_features_tp_fp_distribution_paired.png
  5. auc_bootstrap_ci.png (all features + residualized comparison)

Also saves a summary TSV with point-estimate + CI for each feature.
"""
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LinearRegression

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PILOT = ROOT / "output/hcc1395_normal_pilot"
DATA = ROOT / "research/F_hpfinengroups_deepening/data"
FIG = ROOT / "research/F_hpfinengroups_deepening/figures"
FIG.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")

RNG = np.random.default_rng(2026)


def load_arm(mode):
    sig = pd.read_csv(PILOT / mode / "significance_summary.csv", low_memory=False)
    truth = pd.read_csv(DATA / f"r1_hcc1395_filter_snvs_{'to' if mode == 'TO' else mode}.tsv", sep="\t")
    sig["key"] = sig["Chr"].astype(str) + ":" + sig["Pos"].astype(str)
    truth["key"] = truth["Chr"].astype(str) + ":" + truth["Pos"].astype(str)
    merged = sig.merge(truth[["key", "truth_label"]], on="key", how="left")
    merged["is_tp"] = (merged["truth_label"].astype(str) == "TP").astype(int)
    return merged


def safe_auc(y, x):
    mask = np.isfinite(x) & np.isfinite(y) & (~pd.isna(x)) & (~pd.isna(y))
    if mask.sum() < 20 or len(np.unique(y[mask])) < 2:
        return np.nan
    try:
        auc = roc_auc_score(y[mask], x[mask])
        return max(auc, 1 - auc)
    except Exception:
        return np.nan


def bootstrap_auc(y, x, n=500, seed=2026):
    mask = np.isfinite(x) & np.isfinite(y) & (~pd.isna(x)) & (~pd.isna(y))
    if mask.sum() < 20:
        return np.nan, np.nan, np.nan
    y2, x2 = y[mask].values, x[mask].values
    idx = np.arange(len(y2))
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n):
        sample = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y2[sample])) < 2:
            continue
        auc = roc_auc_score(y2[sample], x2[sample])
        aucs.append(max(auc, 1 - auc))
    if not aucs:
        return np.nan, np.nan, np.nan
    arr = np.array(aucs)
    return float(np.mean(arr)), float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))


def residualize(df, target, covariates):
    sub = df.dropna(subset=[target] + covariates).copy()
    if len(sub) < 20:
        return None
    X = sub[covariates].astype(float).values
    y = sub[target].astype(float).values
    reg = LinearRegression().fit(X, y)
    resid = y - reg.predict(X)
    out = df.copy()
    out[f"{target}__resid"] = np.nan
    out.loc[sub.index, f"{target}__resid"] = resid
    return out


def fig_feature_auc_bar(df_to, df_pa, out):
    features = [
        "NormalBaseline_Coverage", "SampleASM_Delta", "Fisher_Frac_Sig",
        "Epipoly_Delta", "Tumor_HP_Delta", "Entropy_Imbalance",
        "HP_Residual_Delta", "Normal_HP_Delta", "HPFineF",
        "HPFineNGroups", "PairwiseMedianDist", "CramersV",
    ]
    rows = []
    for mode_name, df in [("TO", df_to), ("paired", df_pa)]:
        flt = df[(df["NumReads"] >= 80) & (df["HPFineNGroups"] >= 4)]
        for f in features:
            if f not in flt.columns:
                continue
            mean, lo, hi = bootstrap_auc(flt["is_tp"], flt[f])
            rows.append({"feature": f, "mode": mode_name, "auc": mean, "lo": lo, "hi": hi})
    rdf = pd.DataFrame(rows).dropna(subset=["auc"])
    rdf = rdf.sort_values(["mode", "auc"], ascending=[True, False])

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)
    for ax, mode_name in zip(axes, ["TO", "paired"]):
        sub = rdf[rdf["mode"] == mode_name]
        yerr = np.array([sub["auc"] - sub["lo"], sub["hi"] - sub["auc"]])
        colors = ["#2c7fb8" if a >= 0.65 else ("#7fcdbb" if a >= 0.60 else "#edf8b1") for a in sub["auc"]]
        ax.barh(sub["feature"], sub["auc"], xerr=yerr, color=colors, edgecolor="k")
        ax.axvline(0.5, color="red", ls="--", alpha=0.5, label="random")
        ax.axvline(0.60, color="gray", ls=":", alpha=0.5, label="practical threshold")
        ax.axvline(0.70, color="green", ls=":", alpha=0.5, label="strong signal")
        ax.set_xlabel("AUC (bootstrap mean, 95% CI)")
        ax.set_title(f"{mode_name} arm — F pilot subset (NR≥80, NG≥4)")
        ax.set_xlim(0.45, 0.85)
        ax.invert_yaxis()
        if mode_name == "TO":
            ax.legend(loc="lower right", fontsize=10)
    plt.suptitle("HCC1395 Normal BAM pilot · per-feature AUC with bootstrap CI", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    return rdf


def fig_residualized_comparison(df_to, df_pa, out):
    features_confound = {
        "NormalBaseline_Coverage": ["NumReads", "Coverage_Multiple"],
        "SampleASM_Delta": ["NumReads", "Coverage_Multiple"],
        "Fisher_Frac_Sig": ["Fisher_N_Tested", "NumReads"],
        "Epipoly_Delta": ["NumReads", "Coverage_Multiple"],
    }
    rows = []
    for mode_name, df in [("TO", df_to), ("paired", df_pa)]:
        flt = df[(df["NumReads"] >= 80) & (df["HPFineNGroups"] >= 4)].copy()
        for feat, covs in features_confound.items():
            if feat not in flt.columns:
                continue
            valid_covs = [c for c in covs if c in flt.columns]
            if not valid_covs:
                continue
            m0, l0, h0 = bootstrap_auc(flt["is_tp"], flt[feat])
            rflt = residualize(flt, feat, valid_covs)
            if rflt is None:
                continue
            m1, l1, h1 = bootstrap_auc(rflt["is_tp"], rflt[f"{feat}__resid"])
            rows.append({"feature": feat, "mode": mode_name, "type": "raw", "auc": m0, "lo": l0, "hi": h0})
            rows.append({"feature": feat, "mode": mode_name, "type": "residualized", "auc": m1, "lo": l1, "hi": h1})
    rdf = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)
    for ax, mode_name in zip(axes, ["TO", "paired"]):
        sub = rdf[rdf["mode"] == mode_name]
        if sub.empty:
            continue
        feats = sub["feature"].unique()
        width = 0.35
        x = np.arange(len(feats))
        raw = sub[sub["type"] == "raw"].set_index("feature").reindex(feats)
        res = sub[sub["type"] == "residualized"].set_index("feature").reindex(feats)
        ax.bar(x - width/2, raw["auc"], width, yerr=[raw["auc"]-raw["lo"], raw["hi"]-raw["auc"]],
               label="raw", color="#fdb462", edgecolor="k")
        ax.bar(x + width/2, res["auc"], width, yerr=[res["auc"]-res["lo"], res["hi"]-res["auc"]],
               label="residualized", color="#80b1d3", edgecolor="k")
        ax.axhline(0.5, color="red", ls="--", alpha=0.5)
        ax.axhline(0.60, color="gray", ls=":", alpha=0.5)
        ax.axhline(0.70, color="green", ls=":", alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(feats, rotation=30, ha="right")
        ax.set_ylabel("AUC (bootstrap mean, 95% CI)")
        ax.set_title(f"{mode_name} arm")
        ax.set_ylim(0.45, 0.85)
        if mode_name == "TO":
            ax.legend(loc="upper right")
    plt.suptitle("Residualized AUC vs Raw AUC (排除 coverage / Fisher_N_Tested 等 confound)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    return rdf


def fig_tp_fp_distribution(df, mode_name, out):
    flt = df[(df["NumReads"] >= 80) & (df["HPFineNGroups"] >= 4)].copy()
    feats = ["NormalBaseline_Coverage", "SampleASM_Delta", "Fisher_Frac_Sig", "Epipoly_Delta"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, f in zip(axes.flat, feats):
        if f not in flt.columns:
            ax.axis("off")
            continue
        tp = flt[flt["is_tp"] == 1][f].dropna()
        fp = flt[flt["is_tp"] == 0][f].dropna()
        if len(tp) < 5 or len(fp) < 5:
            ax.axis("off")
            continue
        bins = np.linspace(min(tp.min(), fp.min()), max(tp.max(), fp.max()), 30)
        ax.hist(tp, bins=bins, alpha=0.55, label=f"TP (n={len(tp)})", color="#2c7fb8", density=True)
        ax.hist(fp, bins=bins, alpha=0.55, label=f"FP (n={len(fp)})", color="#d7191c", density=True)
        ax.axvline(tp.median(), color="#2c7fb8", ls="--", lw=2, alpha=0.8)
        ax.axvline(fp.median(), color="#d7191c", ls="--", lw=2, alpha=0.8)
        auc = safe_auc(flt["is_tp"].values, flt[f].values)
        ax.set_title(f"{f}\nAUC={auc:.3f}")
        ax.set_xlabel(f)
        ax.set_ylabel("density")
        ax.legend(loc="upper right", fontsize=10)
    plt.suptitle(f"HCC1395 Normal BAM pilot · {mode_name} arm · TP vs FP distribution\n(F pilot subset: NR≥80, NG≥4)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()


def fig_5axis_heatmap(out):
    matrix_path = ROOT / "docs/reports/research_landscape/data/10_Registry_5axis_matrix.tsv"
    m = pd.read_csv(matrix_path, sep="\t")
    m = m[m["n_rows"] >= 20].copy()
    # Aggregate to mode × zone × af_band (collapse loh, cn)
    agg = m.groupby(["mode", "zone", "af_band"]).agg(
        n_rows=("n_rows", "sum"), n_TP=("n_TP", "sum"), n_FP=("n_FP", "sum")
    ).reset_index()
    agg["TP_rate"] = agg["n_TP"] / (agg["n_TP"] + agg["n_FP"]).replace(0, 1)

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    for ax, mode_name in zip(axes, ["to", "paired"]):
        sub = agg[agg["mode"] == mode_name].pivot(index="zone", columns="af_band", values="TP_rate")
        sub = sub.reindex(["Z1", "Z2", "Z3", "Z4", "Z5"])
        sub = sub[["extreme", "half", "intermediate"]]
        sns.heatmap(sub, annot=True, fmt=".3f", cmap="RdYlGn", vmin=0.5, vmax=1.0, ax=ax,
                    cbar_kws={"label": "TP rate"}, linewidths=0.5, linecolor="white")
        ax.set_title(f"{mode_name.upper()} mode")
        ax.set_xlabel("AF band")
        ax.set_ylabel("Zone")
    plt.suptitle("5-axis Registry matrix · mode × zone × AF band (collapsed LOH/CN)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    df_to = load_arm("TO")
    df_pa = load_arm("paired")
    print(f"TO rows: {len(df_to)}, paired rows: {len(df_pa)}")

    # Figure 1: Full feature AUC bar
    print("Figure 1: feature AUC bar ...")
    rdf1 = fig_feature_auc_bar(df_to, df_pa, FIG / "fig1_feature_auc_bar.png")
    rdf1.to_csv(FIG / "fig1_feature_auc_bar.tsv", sep="\t", index=False)

    # Figure 2: Residualized comparison
    print("Figure 2: residualized comparison ...")
    rdf2 = fig_residualized_comparison(df_to, df_pa, FIG / "fig2_residualized_auc_comparison.png")
    rdf2.to_csv(FIG / "fig2_residualized_auc_comparison.tsv", sep="\t", index=False)

    # Figure 3-4: TP vs FP distribution
    print("Figure 3-4: TP/FP distribution ...")
    fig_tp_fp_distribution(df_to, "TO", FIG / "fig3_tp_fp_distribution_TO.png")
    fig_tp_fp_distribution(df_pa, "paired", FIG / "fig4_tp_fp_distribution_paired.png")

    # Figure 5: 5-axis heatmap
    print("Figure 5: 5-axis heatmap ...")
    fig_5axis_heatmap(FIG / "fig5_5axis_tp_rate_heatmap.png")

    # Summary TSV
    print("\n=== AUC summary (raw) ===")
    print(rdf1.to_string(index=False))
    print("\n=== AUC summary (residualized) ===")
    print(rdf2.to_string(index=False))

    summary = FIG / "summary.txt"
    with open(summary, "w") as f:
        f.write("Raw AUC (bootstrap, F pilot subset)\n")
        f.write(rdf1.to_string(index=False) + "\n\n")
        f.write("Residualized comparison\n")
        f.write(rdf2.to_string(index=False) + "\n")
    print(f"\n[+] Figures: {FIG}")
    print(f"[+] Summary: {summary}")


if __name__ == "__main__":
    main()
