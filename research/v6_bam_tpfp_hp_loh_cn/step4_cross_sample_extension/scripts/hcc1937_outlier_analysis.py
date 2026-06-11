#!/usr/bin/env python3
"""Step 4 stage 4 — HCC1937 outlier analysis.

HCC1937 has BRCA1 mutant + high ploidy CNV-driven FP outlier (marker rate 0.8165
slightly below 0.85 gate). Analyse:
  - Per-cell HCC1937 TP rate vs mean of other 4 samples (HCC1395+H1437+H2009+HCC1954)
  - Identify cells where HCC1937 deviates by > 0.1 absolute
  - Per-chromosome decomposition for HCC1937 FP regions (which chr drives FP)
  - Effect on signature candidates (with vs without HCC1937 sensitivity analysis)

Outputs:
  step4_cross_sample_extension/intermediate/HCC1937_outlier_per_cell.tsv
  step4_cross_sample_extension/intermediate/HCC1937_fp_per_chr.tsv
  step4_cross_sample_extension/step4_HCC1937_outlier_analysis.md
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from _common import (
    INTERMEDIATE_DIR,
    PER_SAMPLE_GRID,
    PER_SAMPLE_MASTER,
    SAMPLES,
    STEP4_DIR,
    annotate_axes,
    joined_subset,
    load_sample_master,
)


def per_cell_outlier():
    grid = pd.read_csv(STEP4_DIR / "step4_per_sample_grid.tsv", sep="\t", low_memory=False)
    others = ["HCC1395", "H1437", "H2009", "HCC1954"]
    rows = []
    for cid, sub in grid.groupby("cell_id"):
        sub_idx = sub.set_index("sample")
        hcc1937_row = sub_idx.loc["HCC1937"] if "HCC1937" in sub_idx.index else None
        other_rates = []
        other_ns = []
        for o in others:
            if o in sub_idx.index and sub_idx.loc[o, "n"] >= 10:
                other_rates.append(sub_idx.loc[o, "TP_rate"])
                other_ns.append(sub_idx.loc[o, "n"])
        if hcc1937_row is None or hcc1937_row["n"] < 10 or not other_rates:
            continue
        other_mean = float(np.nanmean(other_rates))
        other_std = float(np.nanstd(other_rates))
        h_rate = float(hcc1937_row["TP_rate"])
        deviation = h_rate - other_mean
        rows.append({
            "cell_id": cid,
            "loh_side": hcc1937_row["loh_side"],
            "hp_bucket": hcc1937_row["hp_bucket"],
            "cov_bin": hcc1937_row["cov_bin"],
            "HCC1937_TP_rate": h_rate,
            "HCC1937_n": int(hcc1937_row["n"]),
            "HCC1937_n_FP": int(hcc1937_row["n_FP"]),
            "other_mean_TP_rate": other_mean,
            "other_std": other_std,
            "deviation": deviation,
            "abs_deviation": abs(deviation),
            "z_score": deviation / other_std if other_std > 0 else float("nan"),
            "n_others_valid": len(other_rates),
        })
    df = pd.DataFrame(rows).sort_values("abs_deviation", ascending=False)
    out = INTERMEDIATE_DIR / "HCC1937_outlier_per_cell.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {out.name} ({len(df)} cells)")
    return df


def per_chr_fp_analysis():
    """Per-chromosome FP region count + FP rate for HCC1937 vs other 4 samples."""
    rows = []
    for s in SAMPLES:
        try:
            df_all = load_sample_master(s)
            df = annotate_axes(joined_subset(df_all))
        except FileNotFoundError:
            continue
        for chrom, sub in df.groupby("chr"):
            n = len(sub)
            n_fp = int((sub["label"] == "FP").sum())
            n_tp = int((sub["label"] == "TP").sum())
            rows.append({
                "sample": s,
                "chr": chrom,
                "n": n,
                "n_TP": n_tp,
                "n_FP": n_fp,
                "FP_rate": (n_fp / n) if n else float("nan"),
            })
    df = pd.DataFrame(rows)
    pivot = df.pivot_table(index="chr", columns="sample", values="FP_rate")
    out = INTERMEDIATE_DIR / "HCC1937_fp_per_chr.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {out.name} ({len(df)} rows)")
    return df, pivot


def signature_sensitivity():
    """Recompute Wilcoxon excluding HCC1937 to see if signature candidates change."""
    cons = pd.read_csv(STEP4_DIR / "step4_consistency.tsv", sep="\t", low_memory=False)
    n_with = int(cons["signature_candidate_flag"].sum())
    # Hand-rebuild via n=4 (HCC1395+H1437+H2009+HCC1954)
    others = ["HCC1395", "H1437", "H2009", "HCC1954"]
    rows = []
    from scipy import stats as _stats
    summary_path = INTERMEDIATE_DIR / "step4_grid_summary.json"
    with open(summary_path) as f:
        summaries = json.load(f)
    global_rates = {s["sample"]: s["global_TP_rate"] for s in summaries}
    for _, row in cons.iterrows():
        rates = []
        for s in others:
            tp = row.get(f"TP_rate_{s}", float("nan"))
            n = row.get(f"n_{s}", 0)
            if n >= 10 and not np.isnan(tp):
                rates.append(tp - global_rates.get(s, np.nan))
        if len(rates) < 3:
            continue
        deltas = np.array(rates)
        if np.any(deltas != 0):
            try:
                w = _stats.wilcoxon(
                    deltas[deltas != 0],
                    zero_method="wilcox",
                    alternative="two-sided",
                    method="exact",
                )
                wp = float(w.pvalue)
            except Exception:
                wp = float("nan")
        else:
            wp = float("nan")
        n_above = int(np.sum(deltas > 0))
        n_below = int(np.sum(deltas < 0))
        direction = max(n_above, n_below)
        # n=4 exact min p = 0.125 (not achievable strict <0.0625), use direction ≥3 as proxy
        flag = len(deltas) == 4 and direction >= 4 and not np.isnan(wp) and wp <= 0.125
        rows.append({
            "cell_id": row["cell_id"],
            "n4_direction_match": direction,
            "n4_wilcoxon_p": wp,
            "n4_signature_flag_relaxed": flag,
            "original_flag": bool(row["signature_candidate_flag"]),
        })
    df = pd.DataFrame(rows)
    n_without = int(df["n4_signature_flag_relaxed"].sum())
    return df, n_with, n_without


def main():
    print("== HCC1937 outlier — per-cell deviation ==")
    outlier_df = per_cell_outlier()
    print("\nTop 10 cells where HCC1937 deviates most from other 4:")
    print(outlier_df.head(10).to_string(index=False))

    print("\n== HCC1937 per-chromosome FP rate ==")
    per_chr_df, pivot = per_chr_fp_analysis()
    # HCC1937-specific top chr by FP_rate
    h1937 = per_chr_df[per_chr_df["sample"] == "HCC1937"].sort_values("FP_rate", ascending=False)
    print(h1937.head(10).to_string(index=False))

    print("\n== Signature sensitivity (n=5 with HCC1937 vs n=4 without) ==")
    sens_df, n_with, n_without = signature_sensitivity()
    sens_df.to_csv(
        INTERMEDIATE_DIR / "HCC1937_signature_sensitivity.tsv",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    print(f"n=5 signature candidates: {n_with}")
    print(f"n=4 (without HCC1937) signature candidates (≤0.125 relaxed): {n_without}")

    # Write findings .md
    md = []
    md.append("<!--")
    md.append("build_date: 2026-05-15")
    md.append("agent: Step 4 cross-sample extension — HCC1937 outlier analysis")
    md.append("scope: HCC1937 (BRCA1 mutant, high-ploidy CNV-driven FP) vs other 4 samples")
    md.append("-->")
    md.append("")
    md.append("# Step 4 — HCC1937 Outlier Analysis")
    md.append("")
    md.append("## 0. TL;DR")
    md.append("")
    md.append(
        f"- HCC1937 marker rate 0.8165 邊緣低於 0.85 gate（phase D 唯一未過驗收的樣本）。"
    )
    md.append(
        "- 本檔比較 HCC1937 vs 其他 4 樣本（HCC1395/H1437/H2009/HCC1954）的 per-cell TP rate 差異，"
        "並用 n=4 sensitivity 評估訊號是否依賴 HCC1937。"
    )
    md.append("")
    md.append("## 1. Per-cell deviation (top 10 cells)")
    md.append("")
    md.append("| cell_id | HCC1937_TP_rate | other_mean | deviation | HCC1937_n | n_FP | z_score |")
    md.append("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in outlier_df.head(10).iterrows():
        md.append(
            f"| {r['cell_id']} | {r['HCC1937_TP_rate']:.4f} | {r['other_mean_TP_rate']:.4f} | "
            f"{r['deviation']:+.4f} | {int(r['HCC1937_n'])} | {int(r['HCC1937_n_FP'])} | "
            f"{r['z_score']:.2f} |"
        )
    md.append("")
    md.append("## 2. HCC1937 per-chromosome FP rate (top 10 by FP_rate)")
    md.append("")
    md.append("| chr | n | n_FP | FP_rate |")
    md.append("|---|---:|---:|---:|")
    for _, r in h1937.head(10).iterrows():
        md.append(
            f"| {r['chr']} | {int(r['n'])} | {int(r['n_FP'])} | {r['FP_rate']:.4f} |"
        )
    md.append("")
    md.append("## 3. Signature sensitivity (with vs without HCC1937)")
    md.append("")
    md.append(f"- n=5 (含 HCC1937) signature candidates: **{n_with}**")
    md.append(
        f"- n=4 (排除 HCC1937, Wilcoxon p≤0.125 relaxed) signature candidates: **{n_without}**"
    )
    if n_without > n_with:
        md.append(
            "- → 排除 HCC1937 後 candidates 增加，提示 HCC1937 是部分 cell 的訊號擾動者；"
            "n=5 結論需註明 HCC1937 影響。"
        )
    elif n_without < n_with:
        md.append(
            "- → HCC1937 對訊號有貢獻，但因 ploidy/BRCA1 outlier 特性，需以 sensitivity-pass 通過 cells 為主結論。"
        )
    else:
        md.append("- → HCC1937 不影響 signature candidates 集合。")
    md.append("")
    md.append("## 4. BRCA1 / chr17 / 高 ploidy 染色體解釋")
    md.append("")
    # Look for chr17 FP rate in HCC1937 vs others
    chr17 = per_chr_df[per_chr_df["chr"] == "chr17"].sort_values("sample")
    if not chr17.empty:
        md.append("| sample | chr17 n | chr17 n_FP | chr17 FP_rate |")
        md.append("|---|---:|---:|---:|")
        for _, r in chr17.iterrows():
            md.append(
                f"| {r['sample']} | {int(r['n'])} | {int(r['n_FP'])} | {r['FP_rate']:.4f} |"
            )
    md.append("")
    md.append("> BRCA1 在 chr17q21。HCC1937 chr17 FP_rate 與其他樣本對比可揭露 mutation hotspot 對 FP 富集是否有貢獻。")
    md.append("")
    md.append("## 5. 結論")
    md.append("")
    md.append("- 完整數值見 `intermediate/HCC1937_outlier_per_cell.tsv` + `intermediate/HCC1937_fp_per_chr.tsv`")
    md.append(
        "- HCC1937 signature 影響 sensitivity 結論寫入 `step4_signature_candidates.md`"
    )
    with open(STEP4_DIR / "step4_HCC1937_outlier_analysis.md", "w") as f:
        f.write("\n".join(md))
    print(f"\nwrote step4_HCC1937_outlier_analysis.md")


if __name__ == "__main__":
    main()
