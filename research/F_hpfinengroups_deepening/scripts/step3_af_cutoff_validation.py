#!/usr/bin/env python3
"""F Pilot Step 3: NG=4 + AF cutoff 新 filter 驗證

目的：驗證 Step 2 提出的 "NG=4 + AF<0.4 + NR>=80 NonLOH" 作為新 canonical filter

測試：
  1. 全樣本 per-sample TP rate × AF cutoff grid
  2. 新 filter vs 舊 filter 對比（TP rate + coverage trade-off）
  3. AF cutoff sensitivity (0.2 / 0.3 / 0.35 / 0.4 / 0.45 / 0.5)
  4. Coverage_Multiple 分層驗證（CN confound 檢查）
  5. 跨 chr shuffle sanity check（確認不是 spatial confound）

成功標準：
  - HCC1954 TP rate 從 0.497 提升至 > 0.75（NG=4+AF<0.4+NR>=80）
  - 6/7 樣本 TP rate ≥ 0.85
  - Chr shuffle null AUC < 0.55（避免 spatial artifact）
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = OUT_DIR / "data"


def wilson_ci(n_tp, n):
    if n == 0:
        return np.nan, np.nan, np.nan
    rate = n_tp / n
    lo, hi = stats.binomtest(n_tp, n).proportion_ci(confidence_level=0.95)
    return rate, lo, hi


def load_data():
    return pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)


def af_cutoff_grid_per_sample(df: pd.DataFrame, nr_min=80, ng_min=4):
    """Per-sample TP rate across AF upper-bound cutoffs."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= nr_min) & (sub["HPFineNGroups"] >= ng_min)]
    af_cutoffs = [0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1.01]  # 1.01 = all
    rows = []
    for sample in sorted(sub["sample"].unique()):
        s = sub[sub["sample"] == sample]
        for af_cut in af_cutoffs:
            s_cut = s[s["caller_af"] < af_cut]
            n = len(s_cut)
            if n < 5:
                continue
            n_tp = (s_cut["truth_label"] == "TP").sum()
            rate, lo, hi = wilson_ci(n_tp, n)
            rows.append({"sample": sample, "af_cut": af_cut, "n": n, "n_tp": int(n_tp),
                         "tp_rate": rate, "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def compare_old_new(df: pd.DataFrame, af_cut=0.4, nr_min=80, ng_min=4):
    """Compare old vs new canonical filters."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[sub["to_loh_bed_hit"] == False]

    old = sub[(sub["NumReads"] >= nr_min) & (sub["HPFineNGroups"] >= ng_min)]
    new = old[old["caller_af"] < af_cut]

    overall = sub
    rows = []
    for name, s in [("Overall NonLOH", overall),
                    ("Old: NG>=4 + NR>=80", old),
                    (f"New: NG>=4 + NR>=80 + AF<{af_cut}", new)]:
        n = len(s)
        n_tp = (s["truth_label"] == "TP").sum()
        rate, lo, hi = wilson_ci(n_tp, n)
        rows.append({"condition": name, "n": n, "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": lo, "ci_high": hi})

    # Per-sample
    per_sample_rows = []
    for sample in sorted(sub["sample"].unique()):
        s_all = sub[sub["sample"] == sample]
        s_old = old[old["sample"] == sample]
        s_new = new[new["sample"] == sample]
        per_sample_rows.append({
            "sample": sample,
            "overall_n": len(s_all), "overall_tp_rate": (s_all["truth_label"]=="TP").mean(),
            "old_n": len(s_old), "old_tp_rate": (s_old["truth_label"]=="TP").mean() if len(s_old) else np.nan,
            "new_n": len(s_new), "new_tp_rate": (s_new["truth_label"]=="TP").mean() if len(s_new) else np.nan,
        })
    return pd.DataFrame(rows), pd.DataFrame(per_sample_rows)


def coverage_multiple_stratification(df: pd.DataFrame, af_cut=0.4, nr_min=80, ng_min=4):
    """CN confound check via Coverage_Multiple bins."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= nr_min) &
              (sub["HPFineNGroups"] >= ng_min) & (sub["caller_af"] < af_cut)]
    cov_bins = [(0, 1.1), (1.1, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 10.0)]
    rows = []
    for lo, hi in cov_bins:
        s = sub[(sub["Coverage_Multiple"] >= lo) & (sub["Coverage_Multiple"] < hi)]
        if len(s) < 20:
            continue
        n_tp = (s["truth_label"] == "TP").sum()
        rate, cl, ch = wilson_ci(n_tp, len(s))
        rows.append({"cov_bin": f"{lo:.1f}-{hi:.1f}", "n": len(s), "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": cl, "ci_high": ch})
    return pd.DataFrame(rows)


def chr_shuffle_null(df: pd.DataFrame, af_cut=0.4, nr_min=80, ng_min=4, n_shuffle=20):
    """Chr shuffle null: if we randomly relabel truth within each chr, does the filter
    still appear to work? If yes → spatial confound."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[sub["to_loh_bed_hit"] == False].copy()
    sub["y"] = (sub["truth_label"] == "TP").astype(int)
    observed = sub[(sub["NumReads"] >= nr_min) & (sub["HPFineNGroups"] >= ng_min) &
                   (sub["caller_af"] < af_cut)]
    observed_rate = observed["y"].mean()

    null_rates = []
    rng = np.random.default_rng(42)
    for i in range(n_shuffle):
        shuffled = sub.copy()
        shuffled["y_shuffle"] = shuffled.groupby("Chr")["y"].transform(
            lambda x: rng.permutation(x.values))
        filt = shuffled[(shuffled["NumReads"] >= nr_min) & (shuffled["HPFineNGroups"] >= ng_min) &
                        (shuffled["caller_af"] < af_cut)]
        null_rates.append(filt["y_shuffle"].mean())
    null_mean = np.mean(null_rates)
    null_std = np.std(null_rates)
    z = (observed_rate - null_mean) / null_std if null_std > 0 else np.nan
    return {"observed_rate": observed_rate, "null_mean": null_mean,
            "null_std": null_std, "z": z, "n_shuffle": n_shuffle}


def main():
    print("=" * 70)
    print("F Pilot Step 3: NG=4 + AF cutoff 新 filter 驗證")
    print("=" * 70)
    df = load_data()

    print("\n[A] Per-sample AF-cutoff grid (NG>=4, NR>=80, NonLOH)")
    af_grid = af_cutoff_grid_per_sample(df)
    pv_rate = af_grid.pivot(index="sample", columns="af_cut", values="tp_rate").round(4)
    pv_n = af_grid.pivot(index="sample", columns="af_cut", values="n")
    print("  TP rate:")
    print(pv_rate.to_string())
    print("\n  n:")
    print(pv_n.to_string())
    af_grid.to_csv(DATA_DIR / "step3_af_cutoff_grid.tsv", sep="\t", index=False)

    print("\n[B] Old vs New canonical filter (AF<0.4)")
    summary, per_sample = compare_old_new(df, af_cut=0.4)
    print(summary.round(4).to_string(index=False))
    print("\n[B-per-sample] Per-sample old vs new:")
    per_sample["improve_pp"] = (per_sample["new_tp_rate"] - per_sample["old_tp_rate"]) * 100
    per_sample["coverage_loss_pp"] = (1 - per_sample["new_n"] / per_sample["old_n"]) * 100
    print(per_sample.round(4).to_string(index=False))
    summary.to_csv(DATA_DIR / "step3_old_vs_new.tsv", sep="\t", index=False)
    per_sample.to_csv(DATA_DIR / "step3_per_sample_compare.tsv", sep="\t", index=False)

    print("\n[C] Coverage_Multiple stratification (NG>=4 + NR>=80 + AF<0.4 NonLOH)")
    cov = coverage_multiple_stratification(df)
    print(cov.round(4).to_string(index=False))
    cov.to_csv(DATA_DIR / "step3_coverage_strat.tsv", sep="\t", index=False)

    print("\n[D] Chr-shuffle null check (within-chr label permutation)")
    null = chr_shuffle_null(df, n_shuffle=20)
    print(f"  Observed TP rate: {null['observed_rate']:.4f}")
    print(f"  Null mean (shuffled within chr): {null['null_mean']:.4f} ± {null['null_std']:.4f}")
    print(f"  Z-score: {null['z']:.2f}")
    print(f"  n_shuffle: {null['n_shuffle']}")
    if null['z'] > 5:
        print("  → PASS: filter gain not spatial artifact")
    else:
        print("  → WARN: filter gain similar to null distribution")

    # Verdict
    print("\n=== STEP 3 VERDICT ===")
    hc1954 = per_sample[per_sample["sample"] == "HCC1954"]
    colo = per_sample[per_sample["sample"] == "COLO829"]
    if not hc1954.empty:
        print(f"  HCC1954: old {hc1954['old_tp_rate'].values[0]:.3f} → new {hc1954['new_tp_rate'].values[0]:.3f}"
              f" ({hc1954['improve_pp'].values[0]:+.1f}pp)")
    if not colo.empty:
        print(f"  COLO829: old {colo['old_tp_rate'].values[0]:.3f} → new {colo['new_tp_rate'].values[0]:.3f}"
              f" ({colo['improve_pp'].values[0]:+.1f}pp)")
    n_pass = (per_sample["new_tp_rate"] >= 0.85).sum()
    print(f"  Samples with new TP rate >= 0.85: {n_pass}/{len(per_sample)}")


if __name__ == "__main__":
    main()
