#!/usr/bin/env python3
"""F Pilot Step 5 — Purity<1 dilution simulation + AF 低端 strata 穩健性

目的：回應 Opus 4.7 plan B.2-4「Cell line purity=1 如何外推臨床 purity 0.3-0.8」+
  B.2-? 「AF<0.4 filter 在低 AF 區是否有 artifact」。

兩部分（Step 5C CpG island annotation 因 reference 不存在，暫緩）：

Part 5A — Purity dilution simulation：
  細胞系 purity=1 → observed_AF = tumor_cell_VAF。
  對臨床 purity=p：obs_AF(somatic_TP)  = p × tumor_cell_VAF
                    obs_AF(germline_FP) = p × tumor_cell_VAF + (1-p) × 0.5
  用 HCC1954 / HCC1937 / HCC1395 的 NG=4+NR≥80+NonLOH rows 做模擬。
  輸出：purity ∈ {0.3, 0.5, 0.7, 0.9, 1.0} 下 AF<0.4 filter 的 TP retention、FP rejection。

Part 5B — AF 低端 strata 檢查：
  Per-sample TP rate in AF bins [0, 0.1] / [0.1, 0.2] / [0.2, 0.3] / [0.3, 0.4]。
  under NG=4 + NR≥80 + NonLOH。
  目的：確認 AF<0.1 是否有 sequencing error / strand-specific artifact
        （若 AF<0.1 TP rate 突降 → 需考慮再加 AF>0.05 下限）。

背景：
  Step 4 確認 HCC1954 AF<0.4 後 h=+0.654 medium+ POS；但全部訊號來自 purity=1 cell line。
  臨床應用必須考量 normal contamination 稀釋 + caller AF 下限（ClairS-TO 通常 AF≥0.05-0.1）。
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

PURITY_GRID = [0.3, 0.5, 0.7, 0.9, 1.0]
CALLER_AF_LOWER = 0.05  # ClairS-TO 典型下限


def load_data():
    return pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)


# ============================================================
# Part 5A — Purity dilution simulation
# ============================================================
def simulate_observed_af(cell_line_af: np.ndarray, is_tp: np.ndarray, purity: float) -> np.ndarray:
    """Simulate observed_AF under purity p.

    Model:
      - Cell line (purity=1): observed_AF == tumor_cell_VAF.
      - Clinical sample (purity=p):
          TP (somatic):   obs_AF = p × tumor_cell_VAF       (normal tumor-cell free)
          FP (germline):  obs_AF = p × tumor_cell_VAF + (1-p) × 0.5
                          (germline het in normal cells stays at 0.5, diluting CNV drift)

    注意：此為一階近似。真實情況 FP 可能混合 germline het、sequencing error、artifact，
         此模擬假設 FP 主要來自 germline het + CNV drift（Step 4 HCC1954 結論支撐）。
    """
    obs = np.where(is_tp,
                   purity * cell_line_af,
                   purity * cell_line_af + (1 - purity) * 0.5)
    # Apply caller AF lower bound: AF < 0.05 → 視為未被 caller 偵測到
    obs_clipped = np.where(obs >= CALLER_AF_LOWER, obs, np.nan)
    return obs_clipped


def part_5a_purity_simulation(df):
    print("\n" + "=" * 70)
    print("Part 5A — Purity<1 dilution simulation")
    print("=" * 70)
    print("Model: obs_AF_TP = p × cell_line_AF")
    print("       obs_AF_FP = p × cell_line_AF + (1-p) × 0.5  (germline het normal dilution)")
    print(f"Caller AF lower bound: {CALLER_AF_LOWER} (below → dropped as uncalled)")
    print()

    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80) &
              (sub["HPFineNGroups"] == 4)]

    focus_samples = ["HCC1954", "HCC1937", "HCC1395", "H1437", "H2009"]

    rows = []
    for sample in focus_samples:
        s = sub[sub["sample"] == sample].copy()
        if len(s) < 50:
            continue
        is_tp = (s["truth_label"] == "TP").values
        cell_af = s["caller_af"].values

        n_total = len(s)
        n_tp_original = int(is_tp.sum())
        n_fp_original = int((~is_tp).sum())

        for p in PURITY_GRID:
            obs_af = simulate_observed_af(cell_af, is_tp, p)
            valid = ~np.isnan(obs_af)

            # AF<0.4 retention after purity dilution + caller lower bound
            pass_filter = valid & (obs_af < 0.4)
            tp_kept = int((pass_filter & is_tp).sum())
            fp_kept = int((pass_filter & ~is_tp).sum())
            n_kept = tp_kept + fp_kept
            tp_rate = tp_kept / n_kept if n_kept > 0 else np.nan

            # 同時追蹤 caller dropout
            tp_dropped = int((~valid & is_tp).sum())  # 低 AF 被 caller 濾掉（lose TP）
            fp_dropped = int((~valid & ~is_tp).sum())

            rows.append({
                "sample": sample, "purity": p,
                "n_total": n_total,
                "n_tp_in": n_tp_original, "n_fp_in": n_fp_original,
                "tp_kept_af04": tp_kept, "fp_kept_af04": fp_kept,
                "tp_dropped_caller": tp_dropped, "fp_dropped_caller": fp_dropped,
                "tp_rate_af04": tp_rate,
                "tp_recovery": tp_kept / max(n_tp_original, 1),
                "fp_rejection": 1 - fp_kept / max(n_fp_original, 1),
            })

    out = pd.DataFrame(rows)
    print("Per sample × purity × AF<0.4 filter result:")
    display_cols = ["sample", "purity", "n_total", "n_tp_in", "n_fp_in",
                    "tp_kept_af04", "fp_kept_af04",
                    "tp_dropped_caller", "tp_rate_af04",
                    "tp_recovery", "fp_rejection"]
    print(out[display_cols].round(3).to_string(index=False))
    out.to_csv(DATA_DIR / "step5a_purity_simulation.tsv", sep="\t", index=False)

    # Summary: tp_rate_af04 trend across purity
    print("\nTP rate under AF<0.4 filter across purity (summary):")
    pivot = out.pivot(index="sample", columns="purity", values="tp_rate_af04")
    print(pivot.round(3).to_string())

    print("\nTP recovery (TP kept / total TP input):")
    pivot_rec = out.pivot(index="sample", columns="purity", values="tp_recovery")
    print(pivot_rec.round(3).to_string())

    print("\nFP rejection (1 - FP kept / total FP input):")
    pivot_fpr = out.pivot(index="sample", columns="purity", values="fp_rejection")
    print(pivot_fpr.round(3).to_string())

    return out


# ============================================================
# Part 5B — AF 低端 strata 檢查
# ============================================================
def part_5b_af_lowend_strata(df):
    print("\n" + "=" * 70)
    print("Part 5B — AF 低端 strata TP rate 檢查 (NG=4+NR≥80+NonLOH)")
    print("=" * 70)
    print("目的：確認 AF<0.1 是否有 sequencing error / strand artifact 壓低 TP rate")
    print()

    sub = df[(df["mode"] == "to") & df["truth_label"].isin(["TP", "FP"])]
    sub = sub[(sub["to_loh_bed_hit"] == False) & (sub["NumReads"] >= 80) &
              (sub["HPFineNGroups"] == 4)]

    bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 1.0)]

    rows = []
    for sample in sorted(sub["sample"].unique()):
        s = sub[sub["sample"] == sample]
        for lo, hi in bins:
            bin_s = s[(s["caller_af"] >= lo) & (s["caller_af"] < hi)]
            n = len(bin_s)
            if n == 0:
                rows.append({
                    "sample": sample, "af_bin": f"[{lo:.1f},{hi:.1f})",
                    "n": 0, "n_tp": 0, "n_fp": 0, "tp_rate": np.nan,
                })
                continue
            n_tp = int((bin_s["truth_label"] == "TP").sum())
            n_fp = int((bin_s["truth_label"] == "FP").sum())
            tp_rate = n_tp / n
            rows.append({
                "sample": sample, "af_bin": f"[{lo:.1f},{hi:.1f})",
                "n": n, "n_tp": n_tp, "n_fp": n_fp, "tp_rate": tp_rate,
            })

    out = pd.DataFrame(rows)
    print("Per-sample × AF-bin TP rate:")
    pivot = out.pivot(index="sample", columns="af_bin", values="tp_rate")
    print(pivot.round(3).to_string())
    print("\nPer-sample × AF-bin n (count):")
    pivot_n = out.pivot(index="sample", columns="af_bin", values="n")
    print(pivot_n.to_string())

    out.to_csv(DATA_DIR / "step5b_af_lowend_strata.tsv", sep="\t", index=False)

    # Flag any sample where AF[0,0.1] TP rate << AF[0.1,0.2] TP rate
    print("\nAF<0.1 artifact flag (AF[0,0.1] TP rate vs [0.1,0.2]):")
    for sample in sorted(sub["sample"].unique()):
        r1 = out[(out["sample"] == sample) & (out["af_bin"] == "[0.0,0.1)")]["tp_rate"].values
        r2 = out[(out["sample"] == sample) & (out["af_bin"] == "[0.1,0.2)")]["tp_rate"].values
        if len(r1) == 0 or len(r2) == 0:
            continue
        r1v, r2v = r1[0], r2[0]
        if pd.isna(r1v) or pd.isna(r2v):
            flag = "NA"
        elif r1v < r2v - 0.1:
            flag = "⚠ SUSPECT (Δ ≥ 0.1 drop)"
        else:
            flag = "OK"
        print(f"  {sample:20s}: AF[0,0.1]={r1v:.3f} vs [0.1,0.2]={r2v:.3f} → {flag}")
    return out


def main():
    print("=" * 70)
    print("F Pilot Step 5 — Purity dilution simulation + AF 低端 strata")
    print("=" * 70)
    df = load_data()

    _ = part_5a_purity_simulation(df)
    _ = part_5b_af_lowend_strata(df)

    print("\n" + "=" * 70)
    print("STEP 5 VERDICT")
    print("=" * 70)
    print("(5A) Purity dilution: 見 step5a_purity_simulation.tsv → tp_rate_af04 across purity")
    print("(5B) AF 低端 strata: 見 step5b_af_lowend_strata.tsv → AF[0,0.1] artifact flag")
    print("(5C) CpG island annotation: DEFERRED (no reference BED available)")


if __name__ == "__main__":
    main()
