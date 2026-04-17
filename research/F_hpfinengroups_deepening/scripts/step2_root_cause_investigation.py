#!/usr/bin/env python3
"""F Pilot Step 2: 根因調查

針對 Step 1 發現的 3 個關鍵問題：

Q1: NGroups 非單調（N=2 TP rate 0.64 < N=1 的 0.76）— 為什麼？
  假說 H1: NGroups=2 富集 germline ASM（biallelic simple pattern）
  假說 H2: NGroups=2 富集特定 AF bin（germline ≈ 0.5 AF）
  方法: 比較 NGroups=1,2,3,4 的 caller_af 分佈、AlleleDelta、LOH flag

Q2: HCC1954 NG≥4+NR≥80 TP rate 僅 0.497（n=1,622） — 為什麼？
  假說 H1: HCC1954 有樣本特異的 FP pattern
  假說 H2: HCC1954 FP 在特定 AF/LOH/Coverage_Multiple 區段
  方法: HCC1954 NG=4 FP 的特徵 profile vs TP profile

Q3: Paired NG≥4+NR≥80 TP rate 99.85% — 是 filter gain 還是 baseline 就高？
  方法: 比較 Paired overall TP rate、Paired + NG<4、Paired + NR<80
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = OUT_DIR / "data"


def load_data():
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)
    return df


# ============================================================================
# Q1: NGroups 非單調調查
# ============================================================================
def q1_ngroups_profile(df: pd.DataFrame):
    print("\n" + "=" * 70)
    print("Q1: NGroups 非單調根因調查")
    print("=" * 70)

    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80)]

    # Per-NGroups profile
    rows = []
    for ng in [1, 2, 3, 4]:
        s = sub[sub["HPFineNGroups"] == ng]
        if len(s) < 20:
            continue
        rows.append({
            "ng": ng, "n": len(s),
            "tp_rate": (s["truth_label"] == "TP").mean(),
            "af_mean": s["caller_af"].mean(),
            "af_median": s["caller_af"].median(),
            "af_near_half_frac": ((s["caller_af"] >= 0.4) & (s["caller_af"] <= 0.6)).mean(),
            "ad_mean": s["AlleleDelta"].mean(),
            "ad_abs_mean": s["AlleleDelta"].abs().mean(),
            "allele_sig_frac": (s["AlleleSig"] == "SIG").mean() if "AlleleSig" in s.columns else np.nan,
            "hp_fine_sig_frac": (s["HPFineSig"] == "SIG").mean() if "HPFineSig" in s.columns else np.nan,
            "numcpgs_median": s["NumCpGs"].median(),
        })
    summary = pd.DataFrame(rows)
    print("\n[Q1-A] Per-NGroups summary (TO NonLOH NR>=80)")
    print(summary.round(4).to_string(index=False))

    # AF bin 分層
    print("\n[Q1-B] TP rate per (NGroups, AF bin)")
    af_bins = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]
    rows2 = []
    for ng in [1, 2, 3, 4]:
        for af_lo, af_hi in af_bins:
            s = sub[(sub["HPFineNGroups"] == ng) & (sub["caller_af"] >= af_lo) & (sub["caller_af"] < af_hi)]
            if len(s) < 20:
                continue
            rows2.append({"ng": ng, "af_bin": f"{af_lo:.1f}-{af_hi:.1f}",
                          "n": len(s), "tp_rate": (s["truth_label"] == "TP").mean()})
    af_tbl = pd.DataFrame(rows2)
    print(af_tbl.pivot(index="ng", columns="af_bin", values="tp_rate").round(3).to_string())
    print("\n  n:")
    print(af_tbl.pivot(index="ng", columns="af_bin", values="n").to_string())

    # 假說：NGroups=2 是否為 germline ASM proxy?
    print("\n[Q1-C] NGroups=2 是否 germline ASM 富集？")
    n2 = sub[sub["HPFineNGroups"] == 2]
    n1 = sub[sub["HPFineNGroups"] == 1]
    n4 = sub[sub["HPFineNGroups"] == 4]
    for name, s in [("NG=1", n1), ("NG=2", n2), ("NG=4", n4)]:
        af_near_half = ((s["caller_af"] >= 0.45) & (s["caller_af"] <= 0.55)).mean()
        ad_abs_mean = s["AlleleDelta"].abs().mean()
        print(f"  {name}: n={len(s):,}  AF∈[0.45,0.55] frac={af_near_half:.3f}  |AlleleDelta| mean={ad_abs_mean:.3f}")

    # 儲存
    summary.to_csv(DATA_DIR / "step2_q1_ngroups_profile.tsv", sep="\t", index=False)
    af_tbl.to_csv(DATA_DIR / "step2_q1_ng_af_grid.tsv", sep="\t", index=False)


# ============================================================================
# Q2: HCC1954 失效根因
# ============================================================================
def q2_hcc1954_investigation(df: pd.DataFrame):
    print("\n" + "=" * 70)
    print("Q2: HCC1954 NG>=4+NR>=80 TP rate 0.497 根因")
    print("=" * 70)

    # 全域 HCC1954 filter
    hc = df[(df["sample"] == "HCC1954") & (df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    hc = hc[(hc["to_loh_bed_hit"] == False) & (hc["NumReads"] >= 80) & (hc["HPFineNGroups"] >= 4)]
    hc_tp = hc[hc["truth_label"] == "TP"]
    hc_fp = hc[hc["truth_label"] == "FP"]

    print(f"\n  TP: n={len(hc_tp):,}")
    print(f"  FP: n={len(hc_fp):,}")

    # Compare numerical features
    features = ["caller_af", "AlleleDelta", "Coverage_Multiple", "NumReads", "NumCpGs",
                "HPFineF", "HPMergedDelta", "Quality_Score", "HP_Ratio", "HPFineNGroups",
                "PairwiseMedianDist", "caller_dp", "caller_naf"]
    print("\n[Q2-A] Feature profile TP vs FP (HCC1954, NG>=4 NR>=80 NonLOH)")
    rows = []
    for f in features:
        if f not in hc.columns:
            continue
        tp_v = hc_tp[f].dropna()
        fp_v = hc_fp[f].dropna()
        if len(tp_v) < 10 or len(fp_v) < 10:
            continue
        rows.append({
            "feature": f,
            "tp_mean": tp_v.mean(), "fp_mean": fp_v.mean(),
            "tp_median": tp_v.median(), "fp_median": fp_v.median(),
            "delta_mean": tp_v.mean() - fp_v.mean(),
        })
    prof = pd.DataFrame(rows)
    print(prof.round(4).to_string(index=False))

    # AF bin 分層
    print("\n[Q2-B] HCC1954 TP rate per AF bin (NG>=4 NR>=80 NonLOH)")
    rows2 = []
    for af_lo, af_hi in [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]:
        s = hc[(hc["caller_af"] >= af_lo) & (hc["caller_af"] < af_hi)]
        if len(s) == 0:
            continue
        n_tp = (s["truth_label"] == "TP").sum()
        rows2.append({"af_bin": f"{af_lo:.1f}-{af_hi:.1f}", "n": len(s), "n_tp": n_tp,
                      "tp_rate": n_tp / len(s) if len(s) > 0 else np.nan})
    af_tbl = pd.DataFrame(rows2)
    print(af_tbl.round(3).to_string(index=False))

    # HCC1954 overall TP rate comparison
    hc_all = df[(df["sample"] == "HCC1954") & (df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    print(f"\n  HCC1954 overall TO TP rate:     {(hc_all['truth_label']=='TP').mean():.4f} (n={len(hc_all):,})")
    hc_nonloh_all = hc_all[hc_all["to_loh_bed_hit"] == False]
    print(f"  HCC1954 overall NonLOH TP rate: {(hc_nonloh_all['truth_label']=='TP').mean():.4f} (n={len(hc_nonloh_all):,})")
    print(f"  HCC1954 NG>=4+NR>=80 TP rate:   {(hc['truth_label']=='TP').mean():.4f} (n={len(hc):,})")

    prof.to_csv(DATA_DIR / "step2_q2_hcc1954_profile.tsv", sep="\t", index=False)
    af_tbl.to_csv(DATA_DIR / "step2_q2_hcc1954_af.tsv", sep="\t", index=False)


# ============================================================================
# Q3: Paired 99.85% 是否真 gain
# ============================================================================
def q3_paired_baseline(df: pd.DataFrame):
    print("\n" + "=" * 70)
    print("Q3: Paired + NG>=4 + NR>=80 TP rate 99.85% 是真 gain 嗎？")
    print("=" * 70)

    pr = df[(df["mode"] == "paired") & df["truth_label"].isin(["TP", "FP"])]

    conditions = [
        ("All paired (baseline)", pr),
        ("Paired + NR>=80",        pr[pr["NumReads"] >= 80]),
        ("Paired + NG>=4",         pr[pr["HPFineNGroups"] >= 4]),
        ("Paired + NG>=4 + NR>=80", pr[(pr["HPFineNGroups"] >= 4) & (pr["NumReads"] >= 80)]),
        ("Paired + NG<4 + NR<80",  pr[(pr["HPFineNGroups"] < 4) & (pr["NumReads"] < 80)]),
    ]
    rows = []
    for name, s in conditions:
        n = len(s)
        n_tp = (s["truth_label"] == "TP").sum()
        rate = n_tp / n if n > 0 else np.nan
        rows.append({"condition": name, "n": n, "n_tp": int(n_tp), "tp_rate": rate})
    tbl = pd.DataFrame(rows)
    print(tbl.round(4).to_string(index=False))

    baseline = tbl.iloc[0]["tp_rate"]
    filter_rate = tbl.iloc[3]["tp_rate"]
    print(f"\n  Paired overall baseline: {baseline:.4f}")
    print(f"  Paired NG>=4 + NR>=80:   {filter_rate:.4f}")
    print(f"  Gain: {(filter_rate - baseline)*100:+.2f}pp")
    if filter_rate - baseline < 0.02:
        print("  → 結論: 幾乎無 filter gain；Paired baseline 本身就高")
    elif filter_rate - baseline < 0.05:
        print("  → 結論: 微小 gain；主要來自 baseline")
    else:
        print("  → 結論: 真實 gain，filter 有效")

    # Per-sample Paired 查驗
    print("\n[Q3-B] Per-sample Paired baseline vs NG>=4+NR>=80")
    rows2 = []
    for sample in sorted(pr["sample"].unique()):
        s_all = pr[pr["sample"] == sample]
        s_f = s_all[(s_all["HPFineNGroups"] >= 4) & (s_all["NumReads"] >= 80)]
        if len(s_all) < 10:
            continue
        rows2.append({
            "sample": sample,
            "paired_n": len(s_all), "paired_tp_rate": (s_all["truth_label"]=="TP").mean(),
            "filter_n": len(s_f), "filter_tp_rate": (s_f["truth_label"]=="TP").mean() if len(s_f) else np.nan,
        })
    per_s = pd.DataFrame(rows2)
    per_s["gain_pp"] = (per_s["filter_tp_rate"] - per_s["paired_tp_rate"]) * 100
    print(per_s.round(4).to_string(index=False))

    tbl.to_csv(DATA_DIR / "step2_q3_paired_baseline.tsv", sep="\t", index=False)
    per_s.to_csv(DATA_DIR / "step2_q3_paired_per_sample.tsv", sep="\t", index=False)


def main():
    df = load_data()
    q1_ngroups_profile(df)
    q2_hcc1954_investigation(df)
    q3_paired_baseline(df)


if __name__ == "__main__":
    main()
