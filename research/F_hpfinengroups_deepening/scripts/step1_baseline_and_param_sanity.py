#!/usr/bin/env python3
"""F Pilot Step 1: Baseline reproduction + parameter sanity check

目的：
  1. 重現 HPFineNGroups N≥4 + NR≥80 NonLOH TP rate 89.1% baseline
  2. 檢查參數選擇的合理性：
     - NGroups threshold: 1, 2, 3, 4 在固定 NR≥80 的 TP rate
     - NR threshold: 20, 40, 60, 80, 100, 150 在固定 NGroups=4 的 TP rate
     - 兩者聯合 2D grid
  3. 檢查「範圍內外」：
     - 範圍內 (N≥4 + NR≥80 NonLOH)：預期 TP rate ~89%
     - 範圍外：TP rate 多少？coverage 多少？是否符合預期？

假設：
  H_baseline: TO NonLOH + HPFineNGroups=4 + NR≥80 TP rate ≈ 89.1% (memory 記載)
  H_monotone_N: TP rate 隨 NGroups 增加而上升（N=1 < N=2 < N=3 < N=4）
  H_monotone_NR: TP rate 在 N=4 條件下隨 NR 增加而上升（till plateau）

數據源：all_region_rows.tsv.gz (748K rows, 2026-03-27 audit)
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

NR_THRESHOLDS = [20, 40, 60, 80, 100, 150]
NG_THRESHOLDS = [1, 2, 3, 4]


def tp_rate_with_ci(n_tp: int, n_total: int) -> tuple[float, float, float]:
    """Wilson 95% CI."""
    if n_total == 0:
        return np.nan, np.nan, np.nan
    rate = n_tp / n_total
    ci_low, ci_high = stats.binomtest(n_tp, n_total).proportion_ci(confidence_level=0.95)
    return rate, ci_low, ci_high


def load_data():
    print(f"Loading {MASTER.name}...")
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)
    print(f"  rows: {len(df):,}, cols: {df.shape[1]}\n")
    return df


def baseline_reproduction(df: pd.DataFrame) -> dict:
    """Reproduce 89.1% TP rate claim."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80) & (sub["HPFineNGroups"] >= 4)]
    n_tp = (sub["truth_label"] == "TP").sum()
    n_total = len(sub)
    rate, lo, hi = tp_rate_with_ci(n_tp, n_total)
    return {
        "condition": "TO NonLOH + NGroups>=4 + NR>=80",
        "n_tp": int(n_tp), "n_fp": n_total - int(n_tp), "n_total": n_total,
        "tp_rate": rate, "ci_low": lo, "ci_high": hi,
    }


def overall_baseline(df: pd.DataFrame) -> dict:
    """Overall TO NonLOH TP rate (no NG/NR filter) — baseline for gain computation."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[sub["to_loh_bed_hit"] == False]
    n_tp = (sub["truth_label"] == "TP").sum()
    n_total = len(sub)
    rate, lo, hi = tp_rate_with_ci(n_tp, n_total)
    return {
        "condition": "TO NonLOH (no filter)",
        "n_tp": int(n_tp), "n_fp": n_total - int(n_tp), "n_total": n_total,
        "tp_rate": rate, "ci_low": lo, "ci_high": hi,
    }


def scan_ng_threshold(df: pd.DataFrame, nr_min: int) -> pd.DataFrame:
    """At fixed NR threshold, scan NGroups threshold."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= nr_min)]
    rows = []
    # Each exact NGroups value
    for ng_exact in [1, 2, 3, 4]:
        s = sub[sub["HPFineNGroups"] == ng_exact]
        if len(s) < 5:
            continue
        n_tp = (s["truth_label"] == "TP").sum()
        rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
        rows.append({"ng_cut": f"=={ng_exact}", "n": len(s), "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": lo, "ci_high": hi})
    # Threshold-based (>=)
    for ng_min in NG_THRESHOLDS:
        s = sub[sub["HPFineNGroups"] >= ng_min]
        if len(s) < 5:
            continue
        n_tp = (s["truth_label"] == "TP").sum()
        rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
        rows.append({"ng_cut": f">={ng_min}", "n": len(s), "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def scan_nr_threshold(df: pd.DataFrame, ng_min: int) -> pd.DataFrame:
    """At fixed NGroups threshold, scan NR threshold."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["HPFineNGroups"] >= ng_min)]
    rows = []
    for nr_cut in NR_THRESHOLDS:
        s = sub[sub["NumReads"] >= nr_cut]
        if len(s) < 5:
            continue
        n_tp = (s["truth_label"] == "TP").sum()
        rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
        rows.append({"nr_cut": f">={nr_cut}", "n": len(s), "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def grid_2d(df: pd.DataFrame) -> pd.DataFrame:
    """2D NG × NR grid (both threshold-based)."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[sub["to_loh_bed_hit"] == False]
    rows = []
    for ng_min in NG_THRESHOLDS:
        for nr_min in NR_THRESHOLDS:
            s = sub[(sub["HPFineNGroups"] >= ng_min) & (sub["NumReads"] >= nr_min)]
            if len(s) < 20:
                continue
            n_tp = (s["truth_label"] == "TP").sum()
            rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
            rows.append({
                "ng_min": ng_min, "nr_min": nr_min, "n": len(s), "n_tp": int(n_tp),
                "tp_rate": rate, "ci_low": lo, "ci_high": hi,
            })
    return pd.DataFrame(rows)


def out_of_scope_check(df: pd.DataFrame) -> pd.DataFrame:
    """Check TP rate 在範圍外：
      - LOH regions (same NG/NR cuts)
      - Paired mode (same cuts)
      - NGroups<4 regions at NR>=80
      - NGroups=4 at NR<80
    """
    rows = []
    # 1. TO LOH + NG>=4 + NR>=80
    loh = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    loh = loh[(loh["to_loh_bed_hit"] == True) & (loh["NumReads"] >= 80) & (loh["HPFineNGroups"] >= 4)]
    if len(loh) >= 5:
        n_tp = (loh["truth_label"] == "TP").sum()
        rate, lo, hi = tp_rate_with_ci(n_tp, len(loh))
        rows.append({"scope": "TO LOH + NG>=4 + NR>=80", "n": len(loh), "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": lo, "ci_high": hi})

    # 2. Paired + NG>=4 + NR>=80 (no LOH annotation; check all paired)
    pr = df[(df["mode"] == "paired") & df["truth_label"].isin(["TP", "FP"])]
    pr = pr[(pr["NumReads"] >= 80) & (pr["HPFineNGroups"] >= 4)]
    if len(pr) >= 5:
        n_tp = (pr["truth_label"] == "TP").sum()
        rate, lo, hi = tp_rate_with_ci(n_tp, len(pr))
        rows.append({"scope": "Paired (all) + NG>=4 + NR>=80", "n": len(pr), "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": lo, "ci_high": hi})

    # 3. NG<4 at NR>=80 (NonLOH)
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    nonloh = sub[sub["to_loh_bed_hit"] == False]
    s = nonloh[(nonloh["NumReads"] >= 80) & (nonloh["HPFineNGroups"] < 4)]
    n_tp = (s["truth_label"] == "TP").sum()
    rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
    rows.append({"scope": "TO NonLOH + NG<4 + NR>=80", "n": len(s), "n_tp": int(n_tp),
                 "tp_rate": rate, "ci_low": lo, "ci_high": hi})

    # 4. NG>=4 at NR<80
    s = nonloh[(nonloh["NumReads"] < 80) & (nonloh["HPFineNGroups"] >= 4)]
    n_tp = (s["truth_label"] == "TP").sum()
    rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
    rows.append({"scope": "TO NonLOH + NG>=4 + NR<80", "n": len(s), "n_tp": int(n_tp),
                 "tp_rate": rate, "ci_low": lo, "ci_high": hi})

    # 5. NG<4 + NR<80
    s = nonloh[(nonloh["NumReads"] < 80) & (nonloh["HPFineNGroups"] < 4)]
    n_tp = (s["truth_label"] == "TP").sum()
    rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
    rows.append({"scope": "TO NonLOH + NG<4 + NR<80", "n": len(s), "n_tp": int(n_tp),
                 "tp_rate": rate, "ci_low": lo, "ci_high": hi})

    return pd.DataFrame(rows)


def per_sample_baseline(df: pd.DataFrame) -> pd.DataFrame:
    """Per-sample TP rate at the baseline condition."""
    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80) & (sub["HPFineNGroups"] >= 4)]
    rows = []
    for sample in sorted(sub["sample"].unique()):
        s = sub[sub["sample"] == sample]
        n_tp = (s["truth_label"] == "TP").sum()
        rate, lo, hi = tp_rate_with_ci(n_tp, len(s))
        rows.append({"sample": sample, "n": len(s), "n_tp": int(n_tp),
                     "tp_rate": rate, "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def main():
    print("=" * 70)
    print("F Pilot Step 1: Baseline Reproduction + Parameter Sanity")
    print("=" * 70)
    df = load_data()

    print("\n[A] Overall baseline (TO NonLOH, no filter)")
    ob = overall_baseline(df)
    print(f"  n={ob['n_total']:,}  TP rate = {ob['tp_rate']:.4f}  CI=[{ob['ci_low']:.4f}, {ob['ci_high']:.4f}]")

    print("\n[B] Claimed POSITIVE: TO NonLOH + NGroups>=4 + NR>=80")
    br = baseline_reproduction(df)
    print(f"  n={br['n_total']:,}  n_TP={br['n_tp']:,}  n_FP={br['n_fp']:,}")
    print(f"  TP rate = {br['tp_rate']:.4f}  CI=[{br['ci_low']:.4f}, {br['ci_high']:.4f}]")
    print(f"  Gain over overall: +{(br['tp_rate']-ob['tp_rate'])*100:.2f}pp")
    print(f"  claim 記載 = 0.891 (memory)  | 實測 = {br['tp_rate']:.4f}  |  match={abs(br['tp_rate']-0.891)<0.02}")

    print("\n[C] NGroups threshold scan at NR>=80 (NonLOH)")
    ng_scan = scan_ng_threshold(df, nr_min=80)
    print(ng_scan.round(4).to_string(index=False))

    print("\n[D] NR threshold scan at NGroups>=4 (NonLOH)")
    nr_scan = scan_nr_threshold(df, ng_min=4)
    print(nr_scan.round(4).to_string(index=False))

    print("\n[E] 2D Grid (NG_min × NR_min)")
    grid = grid_2d(df)
    pv = grid.pivot(index="ng_min", columns="nr_min", values="tp_rate").round(4)
    print("  TP rate:")
    print(pv.to_string())
    pv_n = grid.pivot(index="ng_min", columns="nr_min", values="n")
    print("\n  n regions:")
    print(pv_n.to_string())

    print("\n[F] Out-of-scope sanity checks")
    oos = out_of_scope_check(df)
    print(oos.round(4).to_string(index=False))

    print("\n[G] Per-sample baseline condition (TO NonLOH + NG>=4 + NR>=80)")
    ps = per_sample_baseline(df)
    print(ps.round(4).to_string(index=False))

    # Save all
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    ng_scan.to_csv(DATA_DIR / "step1_ng_scan.tsv", sep="\t", index=False)
    nr_scan.to_csv(DATA_DIR / "step1_nr_scan.tsv", sep="\t", index=False)
    grid.to_csv(DATA_DIR / "step1_grid_2d.tsv", sep="\t", index=False)
    oos.to_csv(DATA_DIR / "step1_out_of_scope.tsv", sep="\t", index=False)
    ps.to_csv(DATA_DIR / "step1_per_sample.tsv", sep="\t", index=False)

    print("\n=== Summary for Step 1 ===")
    print(f"  Baseline 89.1% reproduction: {br['tp_rate']:.4f} ({'MATCH' if abs(br['tp_rate']-0.891)<0.02 else 'MISMATCH'})")
    print(f"  Overall TO NonLOH TP rate:    {ob['tp_rate']:.4f}")
    print(f"  Gain at (NG>=4, NR>=80):      +{(br['tp_rate']-ob['tp_rate'])*100:.2f}pp")
    print(f"  N samples passing:            {(ps['n']>=20).sum()}/7")


if __name__ == "__main__":
    main()
