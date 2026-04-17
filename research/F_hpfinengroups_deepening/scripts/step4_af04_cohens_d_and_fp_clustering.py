#!/usr/bin/env python3
"""F Pilot Step 4 — AF<0.4 stratified per-sample Cohen's d + FP clustering 分析

三個目的：
  A. 重算 B.1-3 per-sample Cohen's d 但加 AF<0.4 stratification
     - 比較 NG=4 vs NG<4（都在 AF<0.4+NR≥80+NonLOH 內）
     - 確認 HCC1954/COLO829 是否從「特殊」→「POS」
  B. 量化 AF<0.4 對整體與 per-sample 數量影響
     - 若 n(AF<0.4)/n(total NG≥4+NR≥80+NonLOH) 過低 → power 不足
  C. HCC1954 AF≥0.4 FP 是否聚集
     - 染色體分布（2×7 表）
     - 位置聚集度（inter-FP distance）
     - feature profile（NumCpGs、Coverage_Multiple、NR、AlleleDelta）
     - 與 TP（AF<0.4）的 feature 差異

背景：
  Step 3 新 canonical filter NG=4+AF<0.4+NR≥80 NonLOH TP rate=0.9281（舊 0.8912 +3.7pp）
  HCC1954 舊 filter 0.497 → 新 0.707 挽救 +21.0pp
  Step 2 HCC1954 AF binning: AF<0.2 TP=0.874 正常；AF 0.8-1.0 TP=0.022 → FP 極端富集在 AF≥0.4

  本 Step 解答：(1) AF<0.4 重算 Cohen's d 升級 (2) power 是否足夠 (3) FP 聚集根因
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
DATA_DIR.mkdir(exist_ok=True)


def cohens_h(p1, p2):
    """Cohen's h for two proportions (TP rate)."""
    if p1 < 0 or p1 > 1 or p2 < 0 or p2 > 1:
        return np.nan
    # Safe arcsine
    phi1 = 2 * np.arcsin(np.sqrt(p1))
    phi2 = 2 * np.arcsin(np.sqrt(p2))
    return phi1 - phi2


def wilson_ci(n_tp, n):
    if n == 0:
        return np.nan, np.nan, np.nan
    rate = n_tp / n
    lo, hi = stats.binomtest(n_tp, n).proportion_ci(confidence_level=0.95)
    return rate, lo, hi


def load_data():
    return pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)


# ============================================================
# Part A: AF<0.4 stratified per-sample Cohen's d
# ============================================================
def part_a_cohens_d_af04(df):
    """Per-sample Cohen's h (for TP rate) under AF<0.4+NR≥80+NonLOH, NG=4 vs NG<4.

    Compare against the original B.1-3 unstratified baseline.
    """
    print("\n" + "=" * 70)
    print("Part A: AF<0.4 stratified per-sample Cohen's d (TP rate effect size)")
    print("=" * 70)

    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80)]

    rows = []
    for sample in sorted(sub["sample"].unique()):
        s = sub[sub["sample"] == sample]

        # Unstratified (B.1-3 original methodology replica, NR≥80)
        n4_all = s[s["HPFineNGroups"] >= 4]
        nlt_all = s[s["HPFineNGroups"] < 4]
        p4_all = (n4_all["truth_label"] == "TP").mean() if len(n4_all) >= 10 else np.nan
        plt_all = (nlt_all["truth_label"] == "TP").mean() if len(nlt_all) >= 10 else np.nan
        h_all = cohens_h(p4_all, plt_all) if not (np.isnan(p4_all) or np.isnan(plt_all)) else np.nan

        # AF<0.4 stratified
        s_af = s[s["caller_af"] < 0.4]
        n4_af = s_af[s_af["HPFineNGroups"] >= 4]
        nlt_af = s_af[s_af["HPFineNGroups"] < 4]
        p4_af = (n4_af["truth_label"] == "TP").mean() if len(n4_af) >= 10 else np.nan
        plt_af = (nlt_af["truth_label"] == "TP").mean() if len(nlt_af) >= 10 else np.nan
        h_af = cohens_h(p4_af, plt_af) if not (np.isnan(p4_af) or np.isnan(plt_af)) else np.nan

        # AF<0.2 stratified
        s_af20 = s[s["caller_af"] < 0.2]
        n4_af20 = s_af20[s_af20["HPFineNGroups"] >= 4]
        nlt_af20 = s_af20[s_af20["HPFineNGroups"] < 4]
        p4_af20 = (n4_af20["truth_label"] == "TP").mean() if len(n4_af20) >= 10 else np.nan
        plt_af20 = (nlt_af20["truth_label"] == "TP").mean() if len(nlt_af20) >= 10 else np.nan
        h_af20 = cohens_h(p4_af20, plt_af20) if not (np.isnan(p4_af20) or np.isnan(plt_af20)) else np.nan

        rows.append({
            "sample": sample,
            "n4_all": len(n4_all), "plt_all": plt_all, "p4_all": p4_all, "h_all": h_all,
            "n4_af04": len(n4_af), "plt_af04": plt_af, "p4_af04": p4_af, "h_af04": h_af,
            "n4_af02": len(n4_af20), "plt_af02": plt_af20, "p4_af02": p4_af20, "h_af02": h_af20,
        })

    out = pd.DataFrame(rows)

    print("\nPer-sample Cohen's h (for proportion difference):")
    print("  h ~ 0.2 small, 0.5 medium, 0.8 large")
    print("  Positive h means NG=4 has higher TP rate than NG<4.")
    print()
    display_cols = ["sample", "n4_all", "p4_all", "plt_all", "h_all",
                    "n4_af04", "p4_af04", "plt_af04", "h_af04",
                    "n4_af02", "p4_af02", "plt_af02", "h_af02"]
    print(out[display_cols].round(3).to_string(index=False))

    # Classification
    def classify(h):
        if pd.isna(h):
            return "NA"
        if h >= 0.5:
            return "medium+"
        elif h >= 0.2:
            return "small"
        elif h > -0.2:
            return "negligible"
        else:
            return "NEGATIVE"

    out["class_all"] = out["h_all"].apply(classify)
    out["class_af04"] = out["h_af04"].apply(classify)
    out["class_af02"] = out["h_af02"].apply(classify)

    print("\nClassification transitions (all → AF<0.4):")
    for _, row in out.iterrows():
        trans = "" if row["class_all"] == row["class_af04"] else f"  → CHANGED"
        print(f"  {row['sample']:20s}: {row['class_all']:10s} (h={row['h_all']:+.2f}) → {row['class_af04']:10s} (h={row['h_af04']:+.2f}){trans}")

    out.to_csv(DATA_DIR / "step4a_cohens_d_af04.tsv", sep="\t", index=False)
    return out


# ============================================================
# Part B: Coverage impact quantification
# ============================================================
def part_b_coverage_impact(df):
    print("\n" + "=" * 70)
    print("Part B: AF<0.4 coverage impact on filter n")
    print("=" * 70)

    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80)]

    rows = []
    for sample in sorted(sub["sample"].unique()):
        s = sub[sub["sample"] == sample]
        n_old = len(s[s["HPFineNGroups"] >= 4])
        n_new = len(s[(s["HPFineNGroups"] == 4) & (s["caller_af"] < 0.4)])
        n_af02 = len(s[(s["HPFineNGroups"] == 4) & (s["caller_af"] < 0.2)])
        n_af05 = len(s[(s["HPFineNGroups"] == 4) & (s["caller_af"] < 0.5)])

        # TP/FP breakdown under new filter
        new_sub = s[(s["HPFineNGroups"] == 4) & (s["caller_af"] < 0.4)]
        n_tp_new = (new_sub["truth_label"] == "TP").sum()
        n_fp_new = (new_sub["truth_label"] == "FP").sum()

        rows.append({
            "sample": sample,
            "n_old_NGge4_NR80": n_old,
            "n_NG4_AF04_NR80": n_new, "coverage_loss_pct": 100 * (1 - n_new / n_old) if n_old > 0 else np.nan,
            "n_NG4_AF02_NR80": n_af02, "pct_af02_vs_all": 100 * n_af02 / n_old if n_old > 0 else np.nan,
            "n_NG4_AF05_NR80": n_af05,
            "n_TP_new": int(n_tp_new), "n_FP_new": int(n_fp_new),
            "power_adequate": "YES" if n_new >= 100 else ("MARGINAL" if n_new >= 30 else "LOW"),
        })

    out = pd.DataFrame(rows)
    print()
    print(out.round(1).to_string(index=False))
    print("\n説明：")
    print("  power_adequate: YES(n≥100) / MARGINAL(30≤n<100) / LOW(n<30)")
    print("  coverage_loss_pct: 新 filter vs 舊 filter 相對丟失率")

    out.to_csv(DATA_DIR / "step4b_coverage_impact.tsv", sep="\t", index=False)
    return out


# ============================================================
# Part C: HCC1954 AF≥0.4 FP clustering analysis
# ============================================================
def part_c_hcc1954_fp_clustering(df):
    print("\n" + "=" * 70)
    print("Part C: HCC1954 AF≥0.4 FP 聚集與 feature 分析")
    print("=" * 70)

    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80) & (sub["HPFineNGroups"] >= 4)]
    s = sub[sub["sample"] == "HCC1954"].copy()
    s["af_group"] = np.where(s["caller_af"] < 0.4, "AF<0.4", "AF≥0.4")

    print(f"\nHCC1954 NG≥4+NR≥80+NonLOH total: n={len(s)}")

    # (1) Chr distribution of FP across AF groups
    print("\n[C1] 染色體 × AF group 的 FP 計數")
    fp = s[s["truth_label"] == "FP"]
    chr_af = fp.groupby(["Chr", "af_group"]).size().unstack(fill_value=0)
    chr_af["total_fp"] = chr_af.sum(axis=1)
    chr_af["af_ge04_frac"] = chr_af.get("AF≥0.4", 0) / chr_af["total_fp"].replace(0, np.nan)
    chr_af = chr_af.sort_values("total_fp", ascending=False)
    print(chr_af.head(15).round(3).to_string())
    chr_af.to_csv(DATA_DIR / "step4c_hcc1954_chr_af_fp.tsv", sep="\t")

    # (2) AF≥0.4 FP spatial clustering via pairwise distance within-chr
    print("\n[C2] AF≥0.4 FP 空間聚集度（per-chr，inter-FP distance）")
    fp_hi = fp[fp["af_group"] == "AF≥0.4"].copy()
    if "Pos" in fp_hi.columns:
        pos_col = "Pos"
    else:
        pos_col = "Start" if "Start" in fp_hi.columns else None

    if pos_col:
        fp_hi = fp_hi.sort_values(["Chr", pos_col]).copy()
        fp_hi["next_diff"] = fp_hi.groupby("Chr")[pos_col].diff(-1).abs()
        # Compare with TP AF<0.4 inter-distance (baseline dispersion)
        tp_lo = s[(s["truth_label"] == "TP") & (s["af_group"] == "AF<0.4")].sort_values(["Chr", pos_col]).copy()
        tp_lo["next_diff"] = tp_lo.groupby("Chr")[pos_col].diff(-1).abs()

        print(f"\n  AF≥0.4 FP inter-neighbor distance (bp):")
        print(f"    median: {fp_hi['next_diff'].median():,.0f}")
        print(f"    Q1/Q3:  {fp_hi['next_diff'].quantile(0.25):,.0f} / {fp_hi['next_diff'].quantile(0.75):,.0f}")
        print(f"    n pairs: {fp_hi['next_diff'].notna().sum()}")
        print(f"\n  AF<0.4 TP inter-neighbor distance (baseline):")
        print(f"    median: {tp_lo['next_diff'].median():,.0f}")
        print(f"    Q1/Q3:  {tp_lo['next_diff'].quantile(0.25):,.0f} / {tp_lo['next_diff'].quantile(0.75):,.0f}")

        # <100kb "clustered" definition
        clust_frac_fp = (fp_hi["next_diff"] < 1e5).mean()
        clust_frac_tp = (tp_lo["next_diff"] < 1e5).mean()
        print(f"\n  inter-distance < 100kb fraction:")
        print(f"    AF≥0.4 FP: {clust_frac_fp:.3f}")
        print(f"    AF<0.4 TP: {clust_frac_tp:.3f}")
        print(f"    ratio: {clust_frac_fp / max(clust_frac_tp, 1e-6):.2f}× (>1 means FP more clustered)")

    # (3) Feature profile AF<0.4 TP vs AF≥0.4 FP
    print("\n[C3] Feature profile: AF<0.4 TP vs AF≥0.4 FP")
    features = ["NumReads", "NumCpGs", "Coverage_Multiple", "AlleleDelta", "caller_af",
                "PairwiseMedianDist", "HPMergedDelta"]
    features = [f for f in features if f in s.columns]

    tp_lo_sub = s[(s["truth_label"] == "TP") & (s["af_group"] == "AF<0.4")]
    fp_hi_sub = s[(s["truth_label"] == "FP") & (s["af_group"] == "AF≥0.4")]
    tp_hi_sub = s[(s["truth_label"] == "TP") & (s["af_group"] == "AF≥0.4")]
    fp_lo_sub = s[(s["truth_label"] == "FP") & (s["af_group"] == "AF<0.4")]

    feat_rows = []
    for f in features:
        def stats_of(df_sub):
            vals = pd.to_numeric(df_sub[f], errors="coerce").dropna()
            if len(vals) == 0:
                return (np.nan, np.nan, np.nan)
            return (vals.mean(), vals.median(), len(vals))

        m_tp_lo, md_tp_lo, n_tp_lo = stats_of(tp_lo_sub)
        m_fp_hi, md_fp_hi, n_fp_hi = stats_of(fp_hi_sub)
        m_tp_hi, md_tp_hi, n_tp_hi = stats_of(tp_hi_sub)
        m_fp_lo, md_fp_lo, n_fp_lo = stats_of(fp_lo_sub)

        feat_rows.append({
            "feature": f,
            "TP_AF<0.4_mean": m_tp_lo, "TP_AF<0.4_median": md_tp_lo, "n_TP_lo": n_tp_lo,
            "FP_AF≥0.4_mean": m_fp_hi, "FP_AF≥0.4_median": md_fp_hi, "n_FP_hi": n_fp_hi,
            "TP_AF≥0.4_mean": m_tp_hi, "TP_AF≥0.4_median": md_tp_hi, "n_TP_hi": n_tp_hi,
            "FP_AF<0.4_mean": m_fp_lo, "FP_AF<0.4_median": md_fp_lo, "n_FP_lo": n_fp_lo,
        })

    feat_out = pd.DataFrame(feat_rows)
    print()
    print(feat_out.round(3).to_string(index=False))
    feat_out.to_csv(DATA_DIR / "step4c_hcc1954_feature_profile.tsv", sep="\t", index=False)

    # (4) AF≥0.4 FP AlleleDelta 分布 vs TP：
    # 如果 FP AlleleDelta > TP → FP 更像 phased/biallelic（germline-like）
    if "AlleleDelta" in s.columns:
        print("\n[C4] AlleleDelta distribution |AD|（germline-like magnitude）:")
        ad_fp_hi = pd.to_numeric(fp_hi_sub["AlleleDelta"], errors="coerce").dropna().abs()
        ad_tp_lo = pd.to_numeric(tp_lo_sub["AlleleDelta"], errors="coerce").dropna().abs()
        print(f"  |AlleleDelta| AF≥0.4 FP:  mean={ad_fp_hi.mean():.3f}, median={ad_fp_hi.median():.3f}, n={len(ad_fp_hi)}")
        print(f"  |AlleleDelta| AF<0.4 TP:  mean={ad_tp_lo.mean():.3f}, median={ad_tp_lo.median():.3f}, n={len(ad_tp_lo)}")
        # Mann-Whitney U
        if len(ad_fp_hi) > 10 and len(ad_tp_lo) > 10:
            u, p = stats.mannwhitneyu(ad_fp_hi, ad_tp_lo, alternative="two-sided")
            print(f"  Mann-Whitney U p-value: {p:.2e}")

    # (5) Quick cross-sample check: does clustering repeat in HCC1937?
    print("\n[C5] Cross-sample FP AF≥0.4 聚集檢查 (5 samples)")
    cross_rows = []
    for sample in ["HCC1954", "HCC1937", "HCC1395", "H1437", "H2009"]:
        ss = sub[sub["sample"] == sample].copy()
        if len(ss) < 50:
            continue
        ss["af_group"] = np.where(ss["caller_af"] < 0.4, "AF<0.4", "AF≥0.4")
        fp_hi_s = ss[(ss["truth_label"] == "FP") & (ss["af_group"] == "AF≥0.4")]
        tp_lo_s = ss[(ss["truth_label"] == "TP") & (ss["af_group"] == "AF<0.4")]

        if pos_col and len(fp_hi_s) >= 10:
            fp_hi_s = fp_hi_s.sort_values(["Chr", pos_col])
            fp_hi_s_nd = fp_hi_s.groupby("Chr")[pos_col].diff(-1).abs()
            tp_lo_s = tp_lo_s.sort_values(["Chr", pos_col])
            tp_lo_s_nd = tp_lo_s.groupby("Chr")[pos_col].diff(-1).abs()
            clust_fp = (fp_hi_s_nd < 1e5).mean() if len(fp_hi_s_nd) > 0 else np.nan
            clust_tp = (tp_lo_s_nd < 1e5).mean() if len(tp_lo_s_nd) > 0 else np.nan
            ratio = clust_fp / max(clust_tp, 1e-6)
            cross_rows.append({
                "sample": sample, "n_FP_hi": len(fp_hi_s), "n_TP_lo": len(tp_lo_s),
                "clust_FP_hi": clust_fp, "clust_TP_lo": clust_tp, "ratio": ratio,
            })
    if cross_rows:
        print(pd.DataFrame(cross_rows).round(3).to_string(index=False))
        pd.DataFrame(cross_rows).to_csv(DATA_DIR / "step4c_cross_sample_clustering.tsv", sep="\t", index=False)


def main():
    print("=" * 70)
    print("F Pilot Step 4 — AF<0.4 stratified Cohen's d + FP clustering")
    print("=" * 70)
    df = load_data()

    out_a = part_a_cohens_d_af04(df)
    out_b = part_b_coverage_impact(df)
    part_c_hcc1954_fp_clustering(df)

    print("\n" + "=" * 70)
    print("STEP 4 VERDICT")
    print("=" * 70)
    print("(A) AF<0.4 stratified Cohen's h 分類更新 → 看 class_af04 欄")
    print("(B) Per-sample coverage: 見 step4b_coverage_impact.tsv")
    print("(C) HCC1954 FP 聚集分析: 見 step4c_*.tsv")


if __name__ == "__main__":
    main()
