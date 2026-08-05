#!/usr/bin/env python3
"""
跨樣本 ISM 甲基化特徵 AF 梯度分析
目的：驗證 HCC1395 5kHz 的 ISM 判別力 AF 梯度（低 AF 有效 → LOH 失效 + PassedGating 反效果）
     是否在 HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954 中一致

方法：
  1. 讀取各樣本 step05_intersubmod/intersubmod_tp|fp/significance_summary.csv（含 Quality_Score, CramersV, PMD, PassedGating）
  2. 從 step04_benchmark_longphase_to/filtered_snv_tp|fp.vcf.gz 提取 AF
  3. Join on chr:pos:ref:alt → 按 AF bins 分層統計

輸出：research/fp_provenance/20260322_cross_sample_to_loh_analysis/
  - ism_gradient_[sample].tsv        （各 AF bin 的 TP/FP 特徵統計）
  - ism_gradient_summary.tsv         （跨樣本 AF bin 比較表）
  - ism_passgating_enrichment.tsv    （PassedGating 富集比例）
"""
import argparse
import os
import sys
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import SchemaContractError, select_legacy_view

SYNTH_ROOT = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/fp_provenance/20260322_cross_sample_to_loh_analysis")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# HCC1395_5kHz 使用特殊的 rescue_joined_features.tsv（ISM 候選池）
HCC1395_5KHZ_RJF = Path("/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv")

SAMPLES = {
    "HCC1395_DORADO":  "20260315_hcc1395_dorado_to_pilot",
    "COLO829":         "20260317_colo829_to_pilot_1",
    "H1437":           "20260318_h1437_to_pilot_fastresume",
    "H2009":           "20260318_h2009_to_pilot_fastresume",
    "HCC1937":         "20260318_hcc1937_to_pilot_fastresume",
    "HCC1954":         "20260318_hcc1954_to_pilot_fastresume",
}

AF_BINS    = [0, 0.2, 0.4, 0.6, 0.9, 1.01]
AF_LABELS  = ["<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.9", ">=0.9"]
FEATURES   = ["Quality_Score", "PairwiseMedianDist", "CramersV", "GlobalP"]


def attach_historical_legacy_view(df, allow_unversioned_v1=False):
    """Select the historical four-state metric without folding v2 classes."""
    view = select_legacy_view(
        df,
        allow_unversioned_v1=allow_unversioned_v1,
        unversioned_unknown_policy="fail",
    )
    selected = df.copy()
    selected["_verification_legacy_class"] = view.values
    selected["_verification_selection_field"] = view.field
    selected["_verification_schema_status"] = view.schema_status
    return selected

def af_from_vcf(vcf_path):
    """提取 chr:pos:ref:alt → AF 映射，使用 FORMAT/AF 欄位"""
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n", str(vcf_path)],
        capture_output=True, text=True
    )
    af_map = {}
    for line in result.stdout.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, ref, alt, af_str = parts
        try:
            key = f"{chrom}:{pos}:{ref}:{alt}"
            af_map[key] = float(af_str)
        except ValueError:
            pass
    return af_map

def load_significance_summary(csv_path, truth_status, allow_unversioned_v1=False):
    """讀 significance_summary.csv 並加上 truth_status 和 variant_key"""
    df = pd.read_csv(csv_path)
    df["truth_status"] = truth_status
    df["variant_key"] = df["Chr"].astype(str) + ":" + df["Pos"].astype(str) + ":" + df["Ref"].astype(str) + ":" + df["Alt"].astype(str)
    # 標準化 PassedGating
    if "PassedGating" in df.columns:
        df["PassedGating"] = df["PassedGating"].map(lambda x: str(x).lower() in ["true","1","yes"])
    return attach_historical_legacy_view(df, allow_unversioned_v1)

def analyze_sample_full(sample_name, run_dir, allow_unversioned_v1=False):
    """從 significance_summary + VCF 的 AF 做完整 AF 梯度分析"""
    step04 = SYNTH_ROOT / run_dir / "step04_benchmark_longphase_to"
    step05 = SYNTH_ROOT / run_dir / "step05_intersubmod"

    tp_csv  = step05 / "intersubmod_tp" / "significance_summary.csv"
    fp_csv  = step05 / "intersubmod_fp" / "significance_summary.csv"
    tp_vcf  = step04 / "filtered_snv_tp.vcf.gz"
    fp_vcf  = step04 / "filtered_snv_fp.vcf.gz"

    if not tp_csv.exists() or not fp_csv.exists():
        print(f"  [SKIP] {sample_name}: significance_summary.csv not found")
        return None

    print(f"  Loading {sample_name}...")

    # 讀 ISM 特徵
    tp_df = load_significance_summary(tp_csv, "TP", allow_unversioned_v1)
    fp_df = load_significance_summary(fp_csv, "FP", allow_unversioned_v1)
    df = pd.concat([tp_df, fp_df], ignore_index=True)
    print(f"    ISM: TP={len(tp_df)}, FP={len(fp_df)}")

    # 提取 AF
    tp_af = af_from_vcf(tp_vcf)
    fp_af = af_from_vcf(fp_vcf)
    af_map = {**tp_af, **fp_af}
    print(f"    AF map: TP={len(tp_af)}, FP={len(fp_af)}")

    df["af"] = df["variant_key"].map(af_map)
    missing = df["af"].isna().sum()
    if missing > 0:
        print(f"    WARN: {missing}/{len(df)} rows missing AF (key format mismatch)")
    df = df.dropna(subset=["af"])

    return df

def analyze_hcc1395_5khz(allow_unversioned_v1=False):
    """HCC1395 5kHz 使用已有的 rescue_joined_features.tsv"""
    if not HCC1395_5KHZ_RJF.exists():
        return None
    df = pd.read_csv(HCC1395_5KHZ_RJF, sep="\t")
    df["truth_status"] = df["downstream_status"].map({
        "caller_lost_tp": "TP",
        "caller_removed_fp": "FP",
    })
    df = df.dropna(subset=["truth_status", "af"])
    df = attach_historical_legacy_view(df, allow_unversioned_v1)
    print(f"  HCC1395_5kHz: TP={len(df[df.truth_status=='TP'])}, FP={len(df[df.truth_status=='FP'])} (rescue candidate pool)")
    return df

def compute_af_gradient(df, sample_name):
    """計算 AF bins 的 TP/FP 甲基化特徵比較"""
    required = {
        "_verification_legacy_class",
        "_verification_selection_field",
        "_verification_schema_status",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise SchemaContractError(
            "AF gradient requires an explicit historical legacy view; missing "
            + ", ".join(missing)
        )
    selection_fields = df["_verification_selection_field"].dropna().unique().tolist()
    schema_statuses = df["_verification_schema_status"].dropna().unique().tolist()
    if len(selection_fields) != 1 or len(schema_statuses) != 1:
        raise SchemaContractError(
            "AF gradient cannot combine different verification selections: "
            f"fields={selection_fields}, statuses={schema_statuses}"
        )

    rows = []
    for i, label in enumerate(AF_LABELS):
        mask = (df["af"] >= AF_BINS[i]) & (df["af"] < AF_BINS[i+1])
        tp_g = df[mask & (df["truth_status"] == "TP")]
        fp_g = df[mask & (df["truth_status"] == "FP")]

        row = {
            "sample": sample_name,
            "af_bin": label,
            "tp_n": len(tp_g),
            "fp_n": len(fp_g),
            "verification_view": "legacy4_historical",
            "verification_selection_field": selection_fields[0],
            "verification_schema_status": schema_statuses[0],
        }

        # 特徵統計
        for feat in FEATURES:
            if feat in df.columns:
                tv = tp_g[feat].dropna()
                fv = fp_g[feat].dropna()
                row[f"tp_{feat}_med"] = round(float(np.median(tv)), 4) if len(tv) > 0 else None
                row[f"fp_{feat}_med"] = round(float(np.median(fv)), 4) if len(fv) > 0 else None
                if len(tv) >= 3 and len(fv) >= 3:
                    try:
                        _, p = stats.mannwhitneyu(tv, fv, alternative="two-sided")
                        row[f"{feat}_pval"] = round(float(p), 4)
                    except:
                        row[f"{feat}_pval"] = None
                else:
                    row[f"{feat}_pval"] = None

        # PassedGating 富集
        if "PassedGating" in df.columns:
            tp_total = len(tp_g)
            fp_total = len(fp_g)
            tp_pass = tp_g["PassedGating"].sum() if tp_total > 0 else 0
            fp_pass = fp_g["PassedGating"].sum() if fp_total > 0 else 0
            tp_pass_rate = tp_pass / tp_total if tp_total > 0 else 0
            fp_pass_rate = fp_pass / fp_total if fp_total > 0 else 0
            row["tp_pass_rate"] = round(tp_pass_rate * 100, 1)
            row["fp_pass_rate"] = round(fp_pass_rate * 100, 1)
            # TP:FP ratio before/after
            ratio_before = tp_total / fp_total if fp_total > 0 else None
            ratio_after  = tp_pass / fp_pass if fp_pass > 0 else None
            row["ratio_before"] = round(ratio_before, 3) if ratio_before else None
            row["ratio_after"]  = round(ratio_after, 3)  if ratio_after  else None
            if ratio_before and ratio_after:
                row["gate_enrichment"] = round(ratio_after / ratio_before, 3)
            else:
                row["gate_enrichment"] = None

        # Historical four-state VerificationClass metrics. Empty cohorts are N/A,
        # never silently reported as zero percent.
        legacy_field = "_verification_legacy_class"
        tp_noise = (tp_g[legacy_field] == "Noise").mean() * 100 if len(tp_g) > 0 else None
        fp_noise = (fp_g[legacy_field] == "Noise").mean() * 100 if len(fp_g) > 0 else None
        tp_strong = (tp_g[legacy_field] == "Strong").mean() * 100 if len(tp_g) > 0 else None
        fp_strong = (fp_g[legacy_field] == "Strong").mean() * 100 if len(fp_g) > 0 else None
        row["tp_noise_pct"] = round(tp_noise, 1) if tp_noise is not None else None
        row["fp_noise_pct"] = round(fp_noise, 1) if fp_noise is not None else None
        row["tp_strong_pct"] = round(tp_strong, 1) if tp_strong is not None else None
        row["fp_strong_pct"] = round(fp_strong, 1) if fp_strong is not None else None

        rows.append(row)
    return pd.DataFrame(rows)

def print_gradient_table(results):
    """列印跨樣本 PassedGating 富集對比表"""
    print("\n" + "="*100)
    print("PassedGating 富集倍數 vs AF 區間（各樣本）")
    print("格式: 倍數（TP:FP before → after）  * = 反效果（<0.9×）  !! = 強反效果（<0.5×）")
    print("="*100)

    all_samples = sorted(results.keys())
    print(f"  {'AF':10s}", end="")
    for s in all_samples:
        print(f"  {s[:14]:>14s}", end="")
    print()
    print("-"*100)

    for lab in AF_LABELS:
        print(f"  {lab:10s}", end="")
        for s in all_samples:
            df = results.get(s)
            if df is None:
                print(f"  {'N/A':>14s}", end="")
                continue
            row = df[df["af_bin"] == lab]
            if row.empty or row["gate_enrichment"].isna().all():
                print(f"  {'N/A':>14s}", end="")
                continue
            ge = row["gate_enrichment"].values[0]
            if ge is None:
                print(f"  {'N/A':>14s}", end="")
                continue
            if ge < 0.5:
                marker = "!!"
            elif ge < 0.9:
                marker = " *"
            elif ge >= 1.5:
                marker = " ✓"
            else:
                marker = "  "
            print(f"  {ge:6.2f}×{marker}    ", end="")
        print()
    print("="*100)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical raw VerificationClass as the legacy-four-state field; "
            "canonical v2 input does not need this flag"
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    print("=== 跨樣本 ISM 甲基化特徵 AF 梯度分析 ===\n")
    results = {}

    # HCC1395_5kHz 特殊處理
    df_5khz = analyze_hcc1395_5khz(args.allow_unversioned_v1)
    if df_5khz is not None:
        grad = compute_af_gradient(df_5khz, "HCC1395_5kHz")
        grad.to_csv(OUT_DIR / "ism_gradient_HCC1395_5kHz.tsv", sep="\t", index=False)
        results["HCC1395_5kHz"] = grad

    # 其他樣本
    for sample_name, run_dir in SAMPLES.items():
        df = analyze_sample_full(sample_name, run_dir, args.allow_unversioned_v1)
        if df is None:
            continue
        grad = compute_af_gradient(df, sample_name)
        grad.to_csv(OUT_DIR / f"ism_gradient_{sample_name}.tsv", sep="\t", index=False)
        results[sample_name] = grad

    if not results:
        print("ERROR: No results produced")
        sys.exit(1)

    # 合併彙總表
    all_grad = pd.concat(list(results.values()), ignore_index=True)
    all_grad.to_csv(OUT_DIR / "ism_gradient_summary.tsv", sep="\t", index=False)

    # 列印 PassedGating 富集表
    print_gradient_table(results)

    # 列印 QS 中位數比較
    print("\n=== Quality_Score TP vs FP 中位數（各樣本各 AF 區間）===")
    print(f"  {'AF':10s} {'樣本':18s} {'TP_QS':8s} {'FP_QS':8s} {'diff':8s} {'p-val':8s}")
    print("-"*70)
    for s, df in results.items():
        for _, row in df.iterrows():
            tp_qs = row.get("tp_Quality_Score_med")
            fp_qs = row.get("fp_Quality_Score_med")
            pv = row.get("Quality_Score_pval")
            if tp_qs is not None and fp_qs is not None:
                diff = fp_qs - tp_qs
                sig = "**" if pv and pv < 0.01 else ("*" if pv and pv < 0.05 else "  ")
                pv_str = f"{pv:.4f}" if pv is not None else "N/A"
                print(f"  {row['af_bin']:10s} {s:18s} {tp_qs:8.1f} {fp_qs:8.1f} {diff:+8.1f} {pv_str:8s}{sig}")

    print(f"\n輸出目錄: {OUT_DIR}")
    print(f"  - ism_gradient_[sample].tsv")
    print(f"  - ism_gradient_summary.tsv")

if __name__ == "__main__":
    main()
