#!/usr/bin/env python3
"""
跨樣本 TO AF 分佈與 LOH-like FP 富集分析
目的：驗證 HCC1395 5kHz 的 LOH-like FP 富集現象（AF≥0.9 FP/TP 2.14倍）是否普遍存在

輸出：
  research/fp_provenance/20260322_cross_sample_to_loh_analysis/
  - cross_sample_af_summary.tsv
  - cross_sample_loh_enrichment.tsv
  - [sample]_af_bins.tsv（每個樣本的 AF 分佈）
"""
import os
import sys
import gzip
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/fp_provenance/20260322_cross_sample_to_loh_analysis")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    "HCC1395_5kHz": "20260315_hcc1395_to_pilot",
    "HCC1395_DORADO": "20260315_hcc1395_dorado_to_pilot",
    "COLO829": "20260317_colo829_to_pilot_1",
    "H1437": "20260318_h1437_to_pilot_fastresume",
    "H2009": "20260318_h2009_to_pilot_fastresume",
    "HCC1937": "20260318_hcc1937_to_pilot_fastresume",
    "HCC1954": "20260318_hcc1954_to_pilot_fastresume",
}

AF_BINS = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
AF_LABELS = ["0.0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",
             "0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9",">=0.9"]

def get_af_from_vcfgz(vcf_path):
    """從 gzipped VCF 提取 AF 值"""
    afs = []
    result = subprocess.run(
        ["bcftools", "query", "-f", "%AF\n", str(vcf_path)],
        capture_output=True, text=True
    )
    for line in result.stdout.strip().split("\n"):
        line = line.strip()
        if not line or line == ".":
            continue
        try:
            # AF 可能在 FORMAT 欄位（ClairS-TO 輸出 GT:GQ:DP:AF:...）
            afs.append(float(line))
        except ValueError:
            pass
    return afs

def get_af_from_format(vcf_path):
    """從 FORMAT/AF 欄位提取（ClairS-TO 格式）"""
    afs = []
    result = subprocess.run(
        ["bcftools", "query", "-f", "[%AF]\n", str(vcf_path)],
        capture_output=True, text=True
    )
    for line in result.stdout.strip().split("\n"):
        line = line.strip()
        if not line or line == ".":
            continue
        try:
            afs.append(float(line))
        except ValueError:
            pass
    return afs

def analyze_sample(sample_name, run_dir):
    """分析單個樣本的 TP/FP AF 分佈"""
    step_dir = BASE / run_dir / "step04_benchmark_longphase_to"
    tp_vcf = step_dir / "filtered_snv_tp.vcf.gz"
    fp_vcf = step_dir / "filtered_snv_fp.vcf.gz"

    if not tp_vcf.exists():
        print(f"  [SKIP] {sample_name}: {tp_vcf} not found")
        return None

    print(f"  Processing {sample_name}...")

    tp_afs = get_af_from_format(tp_vcf)
    fp_afs = get_af_from_format(fp_vcf)

    if not tp_afs:
        tp_afs = get_af_from_vcfgz(tp_vcf)
    if not fp_afs:
        fp_afs = get_af_from_vcfgz(fp_vcf)

    if not tp_afs or not fp_afs:
        print(f"  [WARN] {sample_name}: No AF values extracted")
        return None

    tp_arr = np.array(tp_afs)
    fp_arr = np.array(fp_afs)

    # AF 分佈（10個 bins）
    tp_hist, _ = np.histogram(tp_arr, bins=AF_BINS)
    fp_hist, _ = np.histogram(fp_arr, bins=AF_BINS)

    # 關鍵 cutoff 統計
    loh_cutoffs = [0.7, 0.8, 0.9]
    stats = {
        "sample": sample_name,
        "tp_total": len(tp_arr),
        "fp_total": len(fp_arr),
        "tp_af_median": float(np.median(tp_arr)),
        "fp_af_median": float(np.median(fp_arr)),
        "tp_af_mean": float(np.mean(tp_arr)),
        "fp_af_mean": float(np.mean(fp_arr)),
        "fp_tp_ratio": len(fp_arr) / len(tp_arr) if len(tp_arr) > 0 else 0,
    }

    for cut in loh_cutoffs:
        tp_loh_n = int((tp_arr >= cut).sum())
        fp_loh_n = int((fp_arr >= cut).sum())
        tp_loh_pct = tp_loh_n / len(tp_arr) * 100 if len(tp_arr) > 0 else 0
        fp_loh_pct = fp_loh_n / len(fp_arr) * 100 if len(fp_arr) > 0 else 0
        enrichment = fp_loh_pct / tp_loh_pct if tp_loh_pct > 0 else 0
        cut_str = str(int(cut*10))
        stats[f"tp_af_ge{cut_str}_n"] = tp_loh_n
        stats[f"fp_af_ge{cut_str}_n"] = fp_loh_n
        stats[f"tp_af_ge{cut_str}_pct"] = round(tp_loh_pct, 2)
        stats[f"fp_af_ge{cut_str}_pct"] = round(fp_loh_pct, 2)
        stats[f"loh{cut_str}_enrichment"] = round(enrichment, 3)

    # 儲存每個樣本的 AF bins TSV
    bins_rows = []
    for i, lab in enumerate(AF_LABELS):
        tp_n = int(tp_hist[i])
        fp_n = int(fp_hist[i])
        bins_rows.append({
            "af_bin": lab,
            "tp_count": tp_n,
            "fp_count": fp_n,
            "tp_pct": round(tp_n / len(tp_arr) * 100, 2) if len(tp_arr) > 0 else 0,
            "fp_pct": round(fp_n / len(fp_arr) * 100, 2) if len(fp_arr) > 0 else 0,
            "fp_tp_pct_ratio": round((fp_n / len(fp_arr)) / (tp_n / len(tp_arr)), 3) if (tp_n > 0 and len(fp_arr) > 0) else 0,
        })
    bins_df = pd.DataFrame(bins_rows)
    bins_df.to_csv(OUT_DIR / f"{sample_name}_af_bins.tsv", sep="\t", index=False)

    return stats

def main():
    print("=== 跨樣本 TO AF LOH 富集分析 ===\n")

    all_stats = []
    for sample_name, run_dir in SAMPLES.items():
        result = analyze_sample(sample_name, run_dir)
        if result:
            all_stats.append(result)

    if not all_stats:
        print("ERROR: No samples processed")
        sys.exit(1)

    # 儲存總覽表
    summary_df = pd.DataFrame(all_stats)
    summary_path = OUT_DIR / "cross_sample_af_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    # 儲存 LOH 富集比較表
    enrichment_cols = ["sample","tp_total","fp_total","fp_tp_ratio",
                       "tp_af_median","fp_af_median",
                       "tp_af_ge9_pct","fp_af_ge9_pct","loh9_enrichment",
                       "tp_af_ge8_pct","fp_af_ge8_pct","loh8_enrichment",
                       "tp_af_ge7_pct","fp_af_ge7_pct","loh7_enrichment"]
    enrich_df = summary_df[[c for c in enrichment_cols if c in summary_df.columns]]
    enrich_path = OUT_DIR / "cross_sample_loh_enrichment.tsv"
    enrich_df.to_csv(enrich_path, sep="\t", index=False)

    # 列印結果
    print("\n" + "="*90)
    print("跨樣本 LOH-like FP 富集摘要（AF ≥ 0.9）")
    print("="*90)
    print(f"{'樣本':20s} {'TP':>7s} {'FP':>7s} {'FP/TP':>6s} {'TP_med':>7s} {'FP_med':>7s} {'TP≥0.9%':>8s} {'FP≥0.9%':>8s} {'富集倍數':>8s}")
    print("-"*90)
    for row in all_stats:
        print(f"{row['sample']:20s} {row['tp_total']:>7d} {row['fp_total']:>7d} "
              f"{row['fp_tp_ratio']:>6.2f} "
              f"{row['tp_af_median']:>7.3f} {row['fp_af_median']:>7.3f} "
              f"{row['tp_af_ge9_pct']:>8.1f} {row['fp_af_ge9_pct']:>8.1f} "
              f"{row['loh9_enrichment']:>8.2f}×")
    print("="*90)

    print(f"\n[HCC1395 5kHz 參考值: TP≥0.9=17.5%, FP≥0.9=37.5%, 富集=2.14×]")
    print(f"\n輸出目錄: {OUT_DIR}")
    print(f"  - cross_sample_af_summary.tsv")
    print(f"  - cross_sample_loh_enrichment.tsv")
    print(f"  - [sample]_af_bins.tsv")

if __name__ == "__main__":
    main()
