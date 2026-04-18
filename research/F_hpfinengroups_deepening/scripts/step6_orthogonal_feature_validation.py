#!/usr/bin/env python3
"""F Pilot Step 6B — Orthogonal feature validation (HPFineNGroups 是否獨佔訊號)

目的：在 NG=4+AF<0.4+NR≥80+NonLOH 子集內檢查其他 heterogeneity-related 特徵
     （PairwiseMedianDist、HPMergedDelta、AlleleDelta）是否提供**獨立於 HPFineNGroups** 的訊號，
     或全部都塌縮到同一潛在 axis。

輸出：
  (a) 特徵 per-feature AUC (TP vs FP, pooled)
  (b) per-sample AUC
  (c) 特徵 vs HPFineNGroups 的 Spearman ρ（檢查冗餘性）
  (d) 特徵在 TP / FP 群的分布差異 (Mann-Whitney U)

背景：
  Step 4 已發現 AF≥0.4 FP 有低 PairwiseMedianDist（0.168 vs TP 0.189，-11%），
  暗示 FP 來自 germline-like reads（單 haplotype 主導）。
  本 Step 探索在 AF<0.4 內 PairwiseMedianDist 是否仍區分 TP/FP，
  以判斷是否值得加入 filter 或作為 orthogonal validation layer。

注意：本 Step 僅做 exploratory analysis，不做 filter 決策。
      filter 升級需等 Step 6A actual purity mixture 驗證 5A simulation。
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = OUT_DIR / "data"
DATA_DIR.mkdir(exist_ok=True)

FEATURES = ["PairwiseMedianDist", "PairwiseMeanDist", "HPMergedDelta", "AlleleDelta", "HPFineF"]


def load_data():
    return pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)


def safe_auc(y_true, scores):
    if len(np.unique(y_true)) < 2 or np.isnan(scores).all():
        return np.nan
    mask = ~np.isnan(scores)
    if mask.sum() < 10:
        return np.nan
    try:
        return roc_auc_score(y_true[mask], scores[mask])
    except ValueError:
        return np.nan


def part_a_pooled_auc(sub):
    print("\n" + "=" * 70)
    print("Part A — Pooled AUC (TP vs FP) within NG=4+AF<0.4+NR≥80+NonLOH")
    print("=" * 70)

    y = (sub["truth_label"] == "TP").astype(int).values
    rows = []
    for feat in FEATURES:
        x = sub[feat].values
        auc = safe_auc(y, x)
        auc_neg = safe_auc(y, -x)  # 取反向 (可能 FP 方向更大)
        best_auc = max(auc, auc_neg) if not (np.isnan(auc) or np.isnan(auc_neg)) else np.nan
        direction = "TP>FP" if (not np.isnan(auc) and auc >= 0.5) else "FP>TP"

        tp_vals = sub.loc[sub["truth_label"] == "TP", feat].dropna()
        fp_vals = sub.loc[sub["truth_label"] == "FP", feat].dropna()
        if len(tp_vals) >= 10 and len(fp_vals) >= 10:
            try:
                u, p = stats.mannwhitneyu(tp_vals, fp_vals, alternative="two-sided")
            except ValueError:
                u, p = np.nan, np.nan
            tp_mean = tp_vals.mean()
            fp_mean = fp_vals.mean()
        else:
            u, p, tp_mean, fp_mean = np.nan, np.nan, np.nan, np.nan

        rows.append({
            "feature": feat,
            "n_TP": int((sub["truth_label"] == "TP").sum()),
            "n_FP": int((sub["truth_label"] == "FP").sum()),
            "auc_raw": auc, "auc_reversed": auc_neg, "best_auc": best_auc,
            "direction": direction,
            "TP_mean": tp_mean, "FP_mean": fp_mean,
            "delta": (tp_mean - fp_mean) if not np.isnan(tp_mean) else np.nan,
            "mannU_p": p,
        })

    out = pd.DataFrame(rows)
    print(out.round(4).to_string(index=False))
    out.to_csv(DATA_DIR / "step6b_pooled_auc.tsv", sep="\t", index=False)
    return out


def part_b_per_sample_auc(sub):
    print("\n" + "=" * 70)
    print("Part B — Per-sample AUC within NG=4+AF<0.4+NR≥80+NonLOH")
    print("=" * 70)

    rows = []
    for sample in sorted(sub["sample"].unique()):
        s = sub[sub["sample"] == sample]
        if len(s) < 50:
            continue
        y = (s["truth_label"] == "TP").astype(int).values
        if len(np.unique(y)) < 2:
            rows.append({"sample": sample, "n": len(s), **{f: np.nan for f in FEATURES}})
            continue
        row = {"sample": sample, "n": len(s),
               "n_TP": int(y.sum()), "n_FP": int((1 - y).sum())}
        for feat in FEATURES:
            x = s[feat].values
            auc = safe_auc(y, x)
            if not np.isnan(auc) and auc < 0.5:
                auc = 1 - auc  # 取較強方向
            row[feat] = auc
        rows.append(row)

    out = pd.DataFrame(rows)
    print(out.round(3).to_string(index=False))
    out.to_csv(DATA_DIR / "step6b_per_sample_auc.tsv", sep="\t", index=False)
    return out


def part_c_correlation_with_ng(sub):
    print("\n" + "=" * 70)
    print("Part C — Spearman ρ(feature, HPFineNGroups) in NG=4+AF<0.4+NR≥80+NonLOH")
    print("=" * 70)
    print("注意：子集已鎖 NG=4 → ρ 應接近 0（無變異）；改比 NG≥2 更廣條件")

    # 擴大到 NG≥2 + AF<0.4 + NR≥80 + NonLOH 看 orthogonality
    wide = sub.copy()
    # 加回 NG≥2 to have NG variation
    print("\n(改用 NG≥2+AF<0.4+NR≥80+NonLOH 以確保 NG variation 存在)")

    rows = []
    for feat in FEATURES:
        mask = wide[feat].notna() & wide["HPFineNGroups"].notna()
        if mask.sum() < 50:
            rho, p = np.nan, np.nan
        else:
            rho, p = stats.spearmanr(wide.loc[mask, feat], wide.loc[mask, "HPFineNGroups"])
        rows.append({"feature": feat, "n": int(mask.sum()), "spearman_rho": rho, "p_value": p})

    out = pd.DataFrame(rows)
    print(out.round(4).to_string(index=False))
    out.to_csv(DATA_DIR / "step6b_correlation_with_ng.tsv", sep="\t", index=False)
    return out


def part_d_combined_score(sub):
    """嘗試 HPFineNGroups + top orthogonal feature 的線性組合是否 AUC 超過單一特徵."""
    print("\n" + "=" * 70)
    print("Part D — Combined score (HPFineNGroups + orthogonal feature)")
    print("=" * 70)

    y = (sub["truth_label"] == "TP").astype(int).values
    base = sub["HPFineNGroups"].values.astype(float)
    auc_base = safe_auc(y, base)
    print(f"Baseline HPFineNGroups AUC: {auc_base:.4f}")

    rows = [{"combo": "HPFineNGroups_only", "auc": auc_base}]
    for feat in FEATURES:
        if feat == "HPFineNGroups":
            continue
        x = sub[feat].values.astype(float)
        mask = ~np.isnan(x) & ~np.isnan(base)
        if mask.sum() < 50:
            continue
        xs = x[mask]
        bs = base[mask]
        ys = y[mask]
        # normalize to [0,1] z-like
        def std_norm(arr):
            m, s = np.nanmean(arr), np.nanstd(arr)
            return (arr - m) / (s if s > 1e-9 else 1)

        bn = std_norm(bs)
        xn = std_norm(xs)
        for alpha in [0.0, 0.3, 0.5, 0.7, 1.0]:
            # alpha=0 → NG only, 1 → feat only
            # test both directions of feat
            for sign in [+1, -1]:
                combo = (1 - alpha) * bn + alpha * sign * xn
                auc_c = safe_auc(ys, combo)
                rows.append({
                    "combo": f"NG*(1-{alpha})+{sign}*{feat}*{alpha}",
                    "auc": auc_c,
                })

    out = pd.DataFrame(rows)
    out_sorted = out.sort_values("auc", ascending=False).head(15)
    print("\nTop 15 combinations:")
    print(out_sorted.round(4).to_string(index=False))
    out.to_csv(DATA_DIR / "step6b_combined_score.tsv", sep="\t", index=False)
    return out


def main():
    print("=" * 70)
    print("F Pilot Step 6B — Orthogonal feature validation")
    print("=" * 70)
    df = load_data()

    # 核心子集：NG=4+AF<0.4+NR≥80+NonLOH
    sub_core = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"]) &
                  (df["to_loh_bed_hit"] == False) & (df["NumReads"] >= 80) &
                  (df["HPFineNGroups"] == 4) & (df["caller_af"] < 0.4)].copy()
    print(f"\nCore subset (NG=4+AF<0.4+NR≥80+NonLOH): n={len(sub_core)}")
    print(f"  TP: {(sub_core['truth_label']=='TP').sum()}, FP: {(sub_core['truth_label']=='FP').sum()}")

    _ = part_a_pooled_auc(sub_core)
    _ = part_b_per_sample_auc(sub_core)

    # Part C: 需要 NG variation，用 NG≥2+AF<0.4+NR≥80+NonLOH
    sub_wide = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"]) &
                  (df["to_loh_bed_hit"] == False) & (df["NumReads"] >= 80) &
                  (df["HPFineNGroups"] >= 2) & (df["caller_af"] < 0.4)].copy()
    print(f"\nWide subset (NG≥2+AF<0.4+NR≥80+NonLOH): n={len(sub_wide)}")
    _ = part_c_correlation_with_ng(sub_wide)

    _ = part_d_combined_score(sub_core)

    print("\n" + "=" * 70)
    print("STEP 6B VERDICT")
    print("=" * 70)
    print("(A) Pooled AUC: feature TP/FP 區分力")
    print("(B) Per-sample AUC: 跨樣本一致性")
    print("(C) Spearman ρ: 與 HPFineNGroups 獨立性")
    print("(D) Combined: HPFineNGroups + feature 組合是否 AUC 升")


if __name__ == "__main__":
    main()
