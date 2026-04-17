#!/usr/bin/env python3
"""P3 confound check: window-level AUC 提升是否由局部 hot spot 驅動？

Hypothesis:
  H0 (confound): Window-level AUC 提升是因為 TP/FP 在基因組上的 auto-correlation
                 （例如 HCC1395 chr8 熱點），而非 aggregation 本身有訊息增益。
  H1 (true gain): 即使排除 top 10% 最極端 windows (TP rate ≥0.9 or ≤0.1)，
                 window AUC 仍 > region AUC。

Method:
  - 重做 H2009 / HCC1937 / HCC1954 (P3 正向最強 3 樣本)
  - 排除 TP rate ≥0.9 or ≤0.1 的 window
  - 比較 window AUC 重算
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT = REPO_ROOT / "research/P3_window_aggregation_pilot/data/hotspot_confound.tsv"

FOCUS = ["H2009", "HCC1937", "HCC1954", "H1437"]
WINDOW_SIZE = 1_000_000


def safe_auc(y, s):
    try:
        return roc_auc_score(y, s)
    except Exception:
        return np.nan


def main():
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)
    rows = []
    for mode in ["to", "paired"]:
        for loh in ["NonLOH"]:
            sub = df[df["mode"] == mode].copy()
            sub = sub[sub["to_loh_bed_hit"] == False] if loh == "NonLOH" else sub
            sub = sub[(sub["NumReads"] >= 40) & sub["truth_label"].isin(["TP", "FP"])].copy()
            sub["y"] = (sub["truth_label"] == "TP").astype(int)
            sub["window_id"] = sub["Chr"].astype(str) + "_" + (sub["Pos"] // WINDOW_SIZE).astype(str)

            for sample in FOCUS:
                s = sub[sub["sample"] == sample]
                if len(s) < 50:
                    continue
                for metric in ["HPFineNGroups", "AlleleDelta"]:
                    s_m = s[[metric, "y", "window_id"]].dropna()
                    if len(s_m) < 50:
                        continue
                    reg_auc = safe_auc(s_m["y"], s_m[metric])
                    agg = s_m.groupby("window_id").agg(
                        n=("y", "size"), n_tp=("y", "sum"),
                        mean=(metric, "mean"), max=(metric, "max"),
                    )
                    agg["tp_rate"] = agg["n_tp"] / agg["n"]
                    agg["y_maj"] = (agg["tp_rate"] >= 0.5).astype(int)
                    clean = agg[agg["n"] >= 3]
                    # All windows
                    full_auc_mean = safe_auc(clean["y_maj"], clean["mean"]) if len(clean) >= 20 else np.nan
                    # Exclude extremes
                    filt = clean[(clean["tp_rate"] > 0.1) & (clean["tp_rate"] < 0.9)]
                    filt_auc_mean = safe_auc(filt["y_maj"], filt["mean"]) if len(filt) >= 20 else np.nan
                    # TP rate matched bins - equalize TP/FP ratio within bins
                    mid_bin = clean[(clean["tp_rate"] >= 0.3) & (clean["tp_rate"] <= 0.7)]
                    mid_auc_mean = safe_auc(mid_bin["y_maj"], mid_bin["mean"]) if len(mid_bin) >= 20 else np.nan

                    rows.append({
                        "sample": sample, "mode": mode, "loh": loh, "metric": metric,
                        "region_auc": reg_auc,
                        "window_auc_all": full_auc_mean,
                        "window_auc_excl_extremes": filt_auc_mean,
                        "window_auc_mid_tp_rate": mid_auc_mean,
                        "n_windows_all": len(clean),
                        "n_windows_excl": len(filt),
                        "n_windows_mid": len(mid_bin),
                    })

    result = pd.DataFrame(rows)
    result.to_csv(OUT, sep="\t", index=False)
    print(result.round(3).to_string(index=False))

    print("\n=== CONFOUND CHECK VERDICT ===")
    # 排除 extremes 後 window_AUC vs region_AUC
    # 計算 Δ = window_auc_excl_extremes - region_auc
    result["delta_excl"] = result["window_auc_excl_extremes"] - result["region_auc"]
    result["delta_mid"] = result["window_auc_mid_tp_rate"] - result["region_auc"]

    print(f"\n  Δ (排除 TP rate 0.1-0.9 之外):")
    print(result.groupby("sample")[["delta_excl", "delta_mid"]].mean().round(3).to_string())

    n_still_pos = ((result["delta_excl"] > 0.03) & result["delta_excl"].notna()).sum()
    print(f"\n  樣本×metric with Δ >+0.03 after excluding extremes: {n_still_pos}/{len(result[result['delta_excl'].notna()])}")
    if n_still_pos >= len(result) * 0.5:
        print("  → 結論: Gain 有潛力不完全是 hot spot confound")
    else:
        print("  → 結論: 大部分 gain 由 hot spot 驅動 — confound 風險高")


if __name__ == "__main__":
    main()
