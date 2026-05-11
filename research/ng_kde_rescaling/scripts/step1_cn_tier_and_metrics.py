#!/usr/bin/env python3
"""
Step 1+2: Apply multiple CN cut strategies, stratify by LOH inner/outer,
          then compute NG distribution + TP rate + Cohen's h per sample.

Outputs:
  - data/cn_tier_coverage_matrix.tsv    (n, %, median CovM per sample × strategy × tier)
  - data/ng_dist_stratified.tsv          (NG distribution per sample × strategy × tier × loh)
  - data/ng_tprate_stratified.tsv        (TP rate + CI)
  - data/af_bin10_ng_crosstab.tsv        (AF 10-bin × NG value TP rate)
  - data/ng_cohen_h_inner_vs_outer.tsv   (effect size with 20260414 comparison)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import (
    AF_BIN_EDGES_FINE,
    AF_BIN_LABELS_FINE,
    CN_STRATEGIES,
    DATA_DIR,
    SAMPLE_ORDER,
    assign_cn_tier,
    quantile_boundaries,
)


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z * np.sqrt((p * (1 - p) + z * z / (4 * n)) / n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def cohens_h(p1: float, p2: float) -> float:
    p1 = np.clip(p1, 1e-9, 1 - 1e-9)
    p2 = np.clip(p2, 1e-9, 1 - 1e-9)
    return 2 * (np.arcsin(np.sqrt(p1)) - np.arcsin(np.sqrt(p2)))


def main() -> None:
    master_path = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
    master = pd.read_csv(master_path, sep="\t", low_memory=False)
    print(f"[step1] loaded master n={len(master)}")

    # Restrict main analysis to paired_full (7 samples)
    # HCC1395 TO handled separately as sidebar
    paired = master[master["mode"] == "paired_full"].copy()
    to_side = master[master["mode"] == "to_pileup"].copy()
    print(f"[step1] paired_full n={len(paired)}, to_pileup sidebar n={len(to_side)}")

    # ------------------------------------------------------------------
    # Assign CN tier for each strategy
    # ------------------------------------------------------------------
    strat_records = []
    for strat_name, edges in CN_STRATEGIES.items():
        col = f"CN_tier_{strat_name}"
        paired[col] = assign_cn_tier(paired["CovM_used"], edges).astype(str)
        to_side[col] = assign_cn_tier(to_side["CovM_used"], edges).astype(str)
        # coverage matrix
        for (sample, tier), grp in paired.groupby(["sample", col], observed=True):
            strat_records.append({
                "sample": sample,
                "strategy": strat_name,
                "tier": tier,
                "n": len(grp),
                "pct_of_sample": len(grp) / (paired["sample"] == sample).sum() * 100,
                "covm_median": float(grp["CovM_used"].median()),
                "covm_mean": float(grp["CovM_used"].mean()),
                "covm_p05": float(grp["CovM_used"].quantile(0.05)),
                "covm_p95": float(grp["CovM_used"].quantile(0.95)),
                "loh_inner_frac": float((grp["loh_side"] == "Inner").mean()),
            })

    # Strategy D: per-sample quantile
    for sample, grp in paired.groupby("sample"):
        qs = quantile_boundaries(grp["CovM_used"])
        col = f"CN_tier_D_quantile"
        paired.loc[paired["sample"] == sample, col] = assign_cn_tier(
            paired.loc[paired["sample"] == sample, "CovM_used"], qs
        ).astype(str)
        for (s, tier), sub in paired[paired["sample"] == sample].groupby(["sample", col], observed=True):
            strat_records.append({
                "sample": sample,
                "strategy": "D_PerSampleQuantile",
                "tier": tier,
                "n": len(sub),
                "pct_of_sample": len(sub) / len(grp) * 100,
                "covm_median": float(sub["CovM_used"].median()),
                "covm_mean": float(sub["CovM_used"].mean()),
                "covm_p05": float(sub["CovM_used"].quantile(0.05)),
                "covm_p95": float(sub["CovM_used"].quantile(0.95)),
                "loh_inner_frac": float((sub["loh_side"] == "Inner").mean()),
            })

    # Strategy E: LOH-aware hybrid
    inner_edges = [0.6, 0.9]
    outer_edges = [0.85, 1.15]
    paired["CN_tier_E_LOHaware"] = np.where(
        paired["loh_side"] == "Inner",
        assign_cn_tier(paired["CovM_used"], inner_edges).astype(str),
        assign_cn_tier(paired["CovM_used"], outer_edges).astype(str),
    )
    # coverage record for E
    for (sample, tier), grp in paired.groupby(["sample", "CN_tier_E_LOHaware"]):
        strat_records.append({
            "sample": sample,
            "strategy": "E_LOHaware",
            "tier": tier,
            "n": len(grp),
            "pct_of_sample": len(grp) / (paired["sample"] == sample).sum() * 100,
            "covm_median": float(grp["CovM_used"].median()),
            "covm_mean": float(grp["CovM_used"].mean()),
            "covm_p05": float(grp["CovM_used"].quantile(0.05)),
            "covm_p95": float(grp["CovM_used"].quantile(0.95)),
            "loh_inner_frac": float((grp["loh_side"] == "Inner").mean()),
        })

    cov_matrix = pd.DataFrame(strat_records)
    cov_matrix.to_csv(DATA_DIR / "cn_tier_coverage_matrix.tsv", sep="\t", index=False)
    print(f"[step1] wrote cn_tier_coverage_matrix ({len(cov_matrix)} rows)")

    # ------------------------------------------------------------------
    # NG distribution + TP rate stratified
    # ------------------------------------------------------------------
    strategies_for_stratify = ["A_Legacy", "B_ClinicalStrict", "C_FineGrained", "F_SEQC2_grounded", "D_PerSampleQuantile", "E_LOHaware"]

    ng_dist_rows = []
    ng_tp_rows = []
    for strat in strategies_for_stratify:
        col = f"CN_tier_{strat}" if strat in CN_STRATEGIES else (
            "CN_tier_D_quantile" if strat == "D_PerSampleQuantile" else "CN_tier_E_LOHaware"
        )
        for (sample, tier, loh), grp in paired.groupby(["sample", col, "loh_side"], observed=True):
            if len(grp) == 0:
                continue
            ng_vc = grp["HPFineNGroups"].value_counts().sort_index()
            for ng_val, cnt in ng_vc.items():
                ng_dist_rows.append({
                    "sample": sample,
                    "strategy": strat,
                    "tier": tier,
                    "loh_side": loh,
                    "NG": int(ng_val),
                    "n": int(cnt),
                    "pct_within_cell": cnt / len(grp) * 100,
                })
            n_all = len(grp)
            tp_all = int(grp["tp_label"].sum())
            tp_rate_all = tp_all / n_all
            lo, hi = wilson_ci(tp_all, n_all)
            row_common = {
                "sample": sample, "strategy": strat, "tier": tier, "loh_side": loh,
                "n": n_all, "n_tp": tp_all, "tp_rate": tp_rate_all,
                "tp_ci_low": lo, "tp_ci_high": hi,
            }
            ng_tp_rows.append({**row_common, "ng_filter": "all"})
            for ng_th in [3, 4]:
                sub = grp[grp["HPFineNGroups"] >= ng_th]
                if len(sub) == 0:
                    continue
                k = int(sub["tp_label"].sum())
                n = len(sub)
                lo2, hi2 = wilson_ci(k, n)
                ng_tp_rows.append({
                    **row_common,
                    "ng_filter": f"NG>={ng_th}",
                    "n": n, "n_tp": k, "tp_rate": k / n,
                    "tp_ci_low": lo2, "tp_ci_high": hi2,
                })
            # canonical F pilot condition: NG=4 + NR>=80
            canonical = grp[(grp["HPFineNGroups"] == 4) & (grp["NumReads"] >= 80)]
            if len(canonical) > 0:
                k = int(canonical["tp_label"].sum())
                n = len(canonical)
                lo2, hi2 = wilson_ci(k, n)
                ng_tp_rows.append({
                    **row_common,
                    "ng_filter": "NG=4_NR>=80",
                    "n": n, "n_tp": k, "tp_rate": k / n,
                    "tp_ci_low": lo2, "tp_ci_high": hi2,
                })

    pd.DataFrame(ng_dist_rows).to_csv(DATA_DIR / "ng_dist_stratified.tsv", sep="\t", index=False)
    tp_df = pd.DataFrame(ng_tp_rows)
    tp_df.to_csv(DATA_DIR / "ng_tprate_stratified.tsv", sep="\t", index=False)
    print(f"[step1] wrote ng_dist_stratified ({len(ng_dist_rows)}) and ng_tprate_stratified ({len(tp_df)})")

    # ------------------------------------------------------------------
    # AF 10-bin × NG crosstab per sample × CN tier A (legacy reference)
    # ------------------------------------------------------------------
    paired["AF_bin10_str"] = paired["AF_bin10"].astype(str)
    af_rows = []
    for (sample, cn_tier, af_bin), grp in paired.groupby(
        ["sample", "CN_tier_A_Legacy", "AF_bin10_str"], observed=True
    ):
        if len(grp) == 0:
            continue
        for ng_val, sub in grp.groupby("HPFineNGroups"):
            if len(sub) == 0:
                continue
            af_rows.append({
                "sample": sample,
                "cn_tier_A": cn_tier,
                "AF_bin10": af_bin,
                "NG": int(ng_val),
                "n": len(sub),
                "n_tp": int(sub["tp_label"].sum()),
                "tp_rate": float(sub["tp_label"].mean()),
            })
    af_df = pd.DataFrame(af_rows)
    af_df.to_csv(DATA_DIR / "af_bin10_ng_crosstab.tsv", sep="\t", index=False)
    print(f"[step1] wrote af_bin10_ng_crosstab ({len(af_df)})")

    # ------------------------------------------------------------------
    # Replication Test A: 20260414 canonical — Intermediate AF vs Extreme AF
    # within CN1 LOH regions → ΔNG_mean + Cohen's h(NG≥2)
    # ------------------------------------------------------------------
    from scipy import stats
    replication_rows = []
    for strat in ["A_Legacy", "B_ClinicalStrict", "F_SEQC2_grounded"]:
        col = f"CN_tier_{strat}"
        for sample, grp in paired.groupby("sample"):
            cn1 = grp[grp[col] == "T0"]
            # 20260414 restricts to LOH inner (Potential_LOH==True) for clean subclonal LOH signal
            cn1_loh = cn1[cn1["loh_side"] == "Inner"]
            if len(cn1_loh) == 0:
                continue
            inter = cn1_loh[cn1_loh["AF_class"] == "Intermediate"]
            extr = cn1_loh[cn1_loh["AF_class"] == "Extreme"]
            nearh = cn1_loh[cn1_loh["AF_class"] == "Near-half"]
            row = {
                "sample": sample, "strategy": strat,
                "n_cn1_loh": len(cn1_loh),
                "n_intermediate": len(inter),
                "n_extreme": len(extr),
                "n_nearhalf": len(nearh),
                "ng_mean_intermediate": float(inter["HPFineNGroups"].mean()) if len(inter) else np.nan,
                "ng_mean_extreme": float(extr["HPFineNGroups"].mean()) if len(extr) else np.nan,
                "ng_mean_nearhalf": float(nearh["HPFineNGroups"].mean()) if len(nearh) else np.nan,
                "delta_ng_inter_vs_extreme": np.nan,
                "cohen_h_ng2_inter_vs_extreme": np.nan,
                "mw_pvalue": np.nan,
            }
            if len(inter) >= 10 and len(extr) >= 10:
                row["delta_ng_inter_vs_extreme"] = row["ng_mean_intermediate"] - row["ng_mean_extreme"]
                p_i_ng2 = (inter["HPFineNGroups"] >= 2).mean()
                p_e_ng2 = (extr["HPFineNGroups"] >= 2).mean()
                row["cohen_h_ng2_inter_vs_extreme"] = float(cohens_h(p_i_ng2, p_e_ng2))
                try:
                    mw = stats.mannwhitneyu(
                        inter["HPFineNGroups"],
                        extr["HPFineNGroups"],
                        alternative="greater",
                    )
                    row["mw_pvalue"] = float(mw.pvalue)
                except ValueError:
                    pass
            replication_rows.append(row)
    rep_df = pd.DataFrame(replication_rows)
    rep_df.to_csv(DATA_DIR / "ng_20260414_replication_inter_vs_extreme.tsv", sep="\t", index=False)
    print(f"[step1] wrote 20260414 replication ({len(rep_df)})")

    # ------------------------------------------------------------------
    # Extension Test B: LOH Inner vs Outer within each CN tier × AF class
    # ------------------------------------------------------------------
    extension_rows = []
    for strat in ["A_Legacy", "B_ClinicalStrict", "F_SEQC2_grounded"]:
        col = f"CN_tier_{strat}"
        for (sample, tier, af_cls), grp in paired.groupby(
            ["sample", col, "AF_class"], observed=True
        ):
            inner = grp[grp["loh_side"] == "Inner"]
            outer = grp[grp["loh_side"] == "Outer"]
            if len(inner) < 20 or len(outer) < 20:
                continue
            row = {
                "sample": sample, "strategy": strat, "tier": tier, "af_class": af_cls,
                "n_inner": len(inner), "n_outer": len(outer),
                "ng_mean_inner": float(inner["HPFineNGroups"].mean()),
                "ng_mean_outer": float(outer["HPFineNGroups"].mean()),
                "delta_ng_inner_vs_outer": float(inner["HPFineNGroups"].mean() - outer["HPFineNGroups"].mean()),
                "tprate_inner": float(inner["tp_label"].mean()),
                "tprate_outer": float(outer["tp_label"].mean()),
            }
            extension_rows.append(row)
    ext_df = pd.DataFrame(extension_rows)
    ext_df.to_csv(DATA_DIR / "ng_loh_inner_vs_outer_extension.tsv", sep="\t", index=False)
    print(f"[step1] wrote LOH inner vs outer extension ({len(ext_df)})")

    # ------------------------------------------------------------------
    # Summary print: 20260414 replication
    # ------------------------------------------------------------------
    print("\n[step1] 20260414 REPLICATION: ΔNG mean (Intermediate AF vs Extreme AF) within CN1 LOH:")
    for strat in ["A_Legacy", "B_ClinicalStrict", "F_SEQC2_grounded"]:
        sub = rep_df[rep_df["strategy"] == strat].set_index("sample").reindex(SAMPLE_ORDER)
        print(f"  === {strat} (CN1 = lowest bin) ===")
        for s, r in sub.iterrows():
            if pd.isna(r.get("delta_ng_inter_vs_extreme")):
                print(f"    {s}: SKIP (n_inter={r.get('n_intermediate', 0)}, n_extreme={r.get('n_extreme', 0)})")
            else:
                pval = r.get("mw_pvalue", np.nan)
                pstr = f"p={pval:.2e}" if pd.notna(pval) else "p=NA"
                print(f"    {s}: ΔNG={r['delta_ng_inter_vs_extreme']:+.3f} (Inter={r['ng_mean_intermediate']:.3f} vs Extreme={r['ng_mean_extreme']:.3f}), h_NG2={r['cohen_h_ng2_inter_vs_extreme']:+.3f}, n_i={int(r['n_intermediate'])}/n_e={int(r['n_extreme'])}, {pstr}")

    print("\n[step1] EXTENSION: LOH Inner vs Outer ΔNG per sample × tier × AF class (summary):")
    for strat in ["A_Legacy"]:
        sub = ext_df[ext_df["strategy"] == strat]
        print(f"  === {strat} ===")
        print(sub.pivot_table(index=["sample", "tier"], columns="af_class",
                              values="delta_ng_inner_vs_outer", aggfunc="first").round(3).to_string())


if __name__ == "__main__":
    main()
