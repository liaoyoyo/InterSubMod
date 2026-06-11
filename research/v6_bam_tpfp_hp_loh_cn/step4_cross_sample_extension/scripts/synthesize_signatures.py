#!/usr/bin/env python3
"""Step 4 stage 5 — Signature candidates synthesis.

Produce step4_signature_candidates.md listing all cells that pass:
  - n=5 valid samples
  - direction concordance ≥ 4/5
  - Wilcoxon p ≤ 0.0625 (n=5 exact min)

Plus effect-size + per-sample rates + HCC1937 outlier sensitivity.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from _common import INTERMEDIATE_DIR, SAMPLES, STEP4_DIR


def main():
    cons = pd.read_csv(STEP4_DIR / "step4_consistency.tsv", sep="\t", low_memory=False)
    sens = pd.read_csv(
        INTERMEDIATE_DIR / "HCC1937_signature_sensitivity.tsv", sep="\t", low_memory=False
    )

    summary_path = INTERMEDIATE_DIR / "step4_grid_summary.json"
    with open(summary_path) as f:
        summaries = json.load(f)
    rates_summary = {s["sample"]: s["global_TP_rate"] for s in summaries}

    candidates = cons[cons["signature_candidate_flag"] == True].copy()
    candidates = candidates.merge(
        sens[["cell_id", "n4_signature_flag_relaxed", "n4_wilcoxon_p", "n4_direction_match"]],
        on="cell_id",
        how="left",
    )

    md = []
    md.append("<!--")
    md.append("build_date: 2026-05-15")
    md.append("agent: Step 4 cross-sample extension — signature candidates synthesis")
    md.append("scope: 5 samples (HCC1395 + H1437 + H2009 + HCC1954 + HCC1937), V6-only")
    md.append("decision_rule: n=5 + direction ≥4/5 + Wilcoxon p ≤ 0.0625 (n=5 exact min)")
    md.append("-->")
    md.append("")
    md.append("# Step 4 — Cross-Sample Signature Candidates")
    md.append("")
    md.append("## 0. TL;DR")
    md.append("")
    md.append(f"- 跨樣本 n=5 signature candidates: **{len(candidates)} cells**")
    md.append(
        f"- 4/5+ direction concordant cells (任意 Wilcoxon): "
        f"**{int((cons['direction_match'] >= 4).sum())}**"
    )
    md.append(
        f"- HCC1937 sensitivity (n=4 排除 HCC1937, p≤0.125 relaxed) candidates: "
        f"**{int(sens['n4_signature_flag_relaxed'].sum())}**"
    )
    md.append("")
    md.append("## 1. Global TP rates (per sample baseline)")
    md.append("")
    md.append("| sample | global_TP_rate |")
    md.append("|---|---:|")
    for s in SAMPLES:
        md.append(
            f"| {s} | "
            f"{rates_summary.get(s, 'NA') if rates_summary.get(s) is None else f'{rates_summary[s]:.4f}'} |"
        )
    md.append("")
    md.append("## 2. Signature candidates (n=5, direction ≥4/5, Wilcoxon p ≤ 0.0625)")
    md.append("")
    if len(candidates):
        md.append(
            "| cell_id | majority_sign | direction (n_above/n_below) | Wilcoxon p | mean Δ vs global | mean n | HCC1395 | H1437 | H2009 | HCC1954 | HCC1937 | n4-sens flag |"
        )
        md.append(
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"
        )
        for _, r in candidates.iterrows():
            md.append(
                f"| {r['cell_id']} | {r['majority_sign']} | "
                f"{int(r['n_above'])}/{int(r['n_below'])} | "
                f"{r['wilcoxon_p']:.4f} | {r['mean_delta_vs_global']:+.4f} | "
                f"{r['mean_n']:.0f} | "
                + " | ".join(
                    f"{r.get(f'TP_rate_{s}', float('nan')):.3f} (n={int(r.get(f'n_{s}', 0))})"
                    for s in SAMPLES
                )
                + f" | {r.get('n4_signature_flag_relaxed', '')} |"
            )
    else:
        md.append("**無 cell 達標。**")
        md.append("")
        md.append("退一步看：direction ≥4/5（無 Wilcoxon 門檻）cells：")
        near = cons[cons["direction_match"] >= 4].sort_values(
            ["mean_delta_vs_global"], ascending=False
        )
        md.append("")
        if len(near):
            md.append(
                "| cell_id | direction (above/below) | Wilcoxon p | mean Δ vs global | mean n | majority |"
            )
            md.append("|---|---|---:|---:|---:|---|")
            for _, r in near.head(20).iterrows():
                md.append(
                    f"| {r['cell_id']} | {int(r['n_above'])}/{int(r['n_below'])} | "
                    f"{r['wilcoxon_p']:.4f} | {r['mean_delta_vs_global']:+.4f} | "
                    f"{r['mean_n']:.0f} | {r['majority_sign']} |"
                )
        else:
            md.append("(none)")
    md.append("")
    md.append("## 3. n=5 exact Wilcoxon 門檻說明")
    md.append("")
    md.append(
        "- n=5 Wilcoxon signed-rank exact 兩側最小 p = 0.0625（5 個全部同號），"
        "故本 step 採 p ≤ 0.0625 作門檻。若 cell n_samples_valid < 5（缺樣本資料），"
        "改用 direction concordance 加效應量。"
    )
    md.append(
        "- HCC1937 移除後 n=4 exact min p = 0.125（4 個全部同號），故 sensitivity 用 p≤0.125 relaxed。"
    )
    md.append("")
    md.append("## 4. 解釋與後續")
    md.append("")
    md.append("- **若 candidates ≥ 1**：升級 candidate cells 為 cross-sample characterization signature，")
    md.append("  並依 plan §H7 conf-guard pass cells 進入 Step 5（H7 needs FDR + permutation pass）。")
    md.append("- **若 candidates = 0**：")
    md.append("  - 仍可看 direction ≥4/5 cells 作 sample-specific characterization 集合")
    md.append("  - HCC1395 + HCC1937 雙 outlier 可能在 same_HP1/cov_normal 與 cross_het 異常")
    md.append("  - 跨樣本 ceiling effect（V6 marker rate 全 ≥0.85）使 TP rate 已接近 1.0 → small absolute delta 難達 Wilcoxon 顯著")
    md.append("")
    md.append("## 5. 完整檔案")
    md.append("")
    md.append("- `step4_consistency.tsv` — 50 cells × {direction, Wilcoxon, per-sample rates}")
    md.append("- `intermediate/HCC1937_outlier_per_cell.tsv` — HCC1937 vs others deviation per cell")
    md.append("- `intermediate/HCC1937_fp_per_chr.tsv` — per-chr FP rate breakdown")
    md.append("- `intermediate/HCC1937_signature_sensitivity.tsv` — n=4 排除 HCC1937 重算")
    md.append("- `step4_per_sample_grid.tsv` — 5 樣本 × 50 cell 完整 grid")
    md.append("- `step4_HCC1937_outlier_analysis.md` — HCC1937 outlier 詳細分析")

    out = STEP4_DIR / "step4_signature_candidates.md"
    with open(out, "w") as f:
        f.write("\n".join(md))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
