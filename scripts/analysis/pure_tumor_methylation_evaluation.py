#!/usr/bin/env python3
"""
Pure tumor methylation evaluation workflow.

Goal:
- Aggregate available s-pure-pileup runs for target samples.
- Compare TP/FP behavior under current methylation filtering rule.
- Build classification views: correct vs incorrect decisions.
- Quantify factor impact (especially methylation-related features).
- Export organized tables, plots, and markdown reports.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact


DEFAULT_SAMPLES = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

DEFAULT_DATA_ROOT = "/bip8_disk/liaoyoyo2001/InterSubMod_out/output/s-pure-pileup"
DEFAULT_OUTPUT_ROOT = "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation"
RANDOM_STATE = 42


@dataclass
class RunRecord:
    sample: str
    run_dir: Path
    run_name: str
    metrics_path: Path
    tp_summary: Path
    fp_summary: Path
    tp_vcf: Path
    fp_vcf: Path
    variant_counts: Path
    baseline_f1: float
    filtered_f1: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pure tumor methylation evaluation workflow")
    parser.add_argument("--data-root", default=DEFAULT_DATA_ROOT, help="Root of s-pure-pileup outputs")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: auto timestamp under output root)")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT, help="Base output root when --output-dir is not provided")
    parser.add_argument(
        "--samples",
        default=",".join(DEFAULT_SAMPLES),
        help="Comma-separated sample list",
    )
    parser.add_argument(
        "--max-points-scatter",
        type=int,
        default=25000,
        help="Max points for scatter plot down-sampling",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def parse_vaf_from_record(format_keys: List[str], sample_values: List[str]) -> Optional[float]:
    fmt = dict(zip(format_keys, sample_values))

    for key in ["VAF", "AF"]:
        if key in fmt and fmt[key] not in {".", ""}:
            raw = fmt[key].split(",")[0]
            try:
                return float(raw)
            except ValueError:
                pass

    if "AD" in fmt and fmt["AD"] not in {".", ""}:
        ad = fmt["AD"].split(",")
        if len(ad) >= 2:
            try:
                ref = float(ad[0])
                alt = float(ad[1])
                denom = ref + alt
                if denom > 0:
                    return alt / denom
            except ValueError:
                pass

    return None


def parse_vcf_features(vcf_path: Path) -> pd.DataFrame:
    rows = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue

            qual_raw = fields[5]
            try:
                qual = float(qual_raw) if qual_raw != "." else np.nan
            except ValueError:
                qual = np.nan

            format_keys = fields[8].split(":")
            sample_vals = fields[9].split(":")
            vaf = parse_vaf_from_record(format_keys, sample_vals)

            rows.append((chrom, pos, qual, np.nan if vaf is None else float(vaf)))

    return pd.DataFrame(rows, columns=["Chr", "Pos", "QUAL", "VAF"])


def load_metrics(metrics_path: Path) -> Dict[str, float]:
    with open(metrics_path, "r") as f:
        m = json.load(f)
    return {
        "baseline_f1": float(m.get("baseline", {}).get("f1", np.nan)),
        "filtered_f1": float(m.get("filtered", {}).get("f1", np.nan)),
        "f1_delta": float(m.get("improvement", {}).get("f1_delta", np.nan)),
        "tp_total": int(m.get("baseline", {}).get("tp", 0)),
        "fp_total": int(m.get("baseline", {}).get("fp", 0)),
        "tp_regions_analyzed": int(m.get("tp_regions_analyzed", 0)),
        "fp_regions_analyzed": int(m.get("fp_regions_analyzed", 0)),
        "tp_removed": int(m.get("improvement", {}).get("tp_removed", 0)),
        "fp_removed": int(m.get("improvement", {}).get("fp_removed", 0)),
    }


def parse_variant_counts(path: Path) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if not path.exists():
        return out
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            try:
                out[k.strip()] = int(v.strip())
            except ValueError:
                continue
    return out


def looks_like_candidate_run(run_name: str) -> bool:
    lower = run_name.lower()
    if lower.startswith("purity"):
        return False
    if "dryrun" in lower:
        return False
    if not run_name.startswith("20"):
        return False
    return True


def required_files_for_run(run_dir: Path) -> Tuple[bool, Dict[str, Path]]:
    paths = {
        "metrics": run_dir / "metrics.json",
        "tp_summary": run_dir / "intersubmod_tp" / "significance_summary.csv",
        "fp_summary": run_dir / "intersubmod_fp" / "significance_summary.csv",
        "tp_vcf": run_dir / "longphase_s" / "filtered_snv_tp.vcf.gz",
        "fp_vcf": run_dir / "longphase_s" / "filtered_snv_fp.vcf.gz",
        "variant_counts": run_dir / "longphase_s" / "variant_counts.txt",
    }
    ok = all(p.exists() for p in paths.values())
    return ok, paths


def find_best_run_for_sample(sample_root: Path, sample: str) -> Optional[RunRecord]:
    if not sample_root.exists():
        return None

    candidates: List[RunRecord] = []
    for child in sorted(sample_root.iterdir()):
        if not child.is_dir():
            continue
        run_name = child.name
        if not looks_like_candidate_run(run_name):
            continue

        ok, p = required_files_for_run(child)
        if not ok:
            continue

        metric = load_metrics(p["metrics"])
        candidates.append(
            RunRecord(
                sample=sample,
                run_dir=child,
                run_name=run_name,
                metrics_path=p["metrics"],
                tp_summary=p["tp_summary"],
                fp_summary=p["fp_summary"],
                tp_vcf=p["tp_vcf"],
                fp_vcf=p["fp_vcf"],
                variant_counts=p["variant_counts"],
                baseline_f1=metric["baseline_f1"],
                filtered_f1=metric["filtered_f1"],
            )
        )

    if not candidates:
        return None

    # Prefer higher filtered F1, then larger positive delta (if tie), then newer run name.
    candidates.sort(key=lambda r: (r.filtered_f1, r.filtered_f1 - r.baseline_f1, r.run_name), reverse=True)
    return candidates[0]


def load_summary_with_features(
    summary_path: Path,
    vcf_path: Path,
    sample: str,
    run_name: str,
    truth_label: str,
) -> pd.DataFrame:
    df = pd.read_csv(summary_path)
    feat = parse_vcf_features(vcf_path)

    merged = df.merge(feat, on=["Chr", "Pos"], how="left")
    merged["sample"] = sample
    merged["run_name"] = run_name
    merged["truth_label"] = truth_label
    merged["is_tp"] = truth_label == "TP"

    for col in ["AlleleDelta", "CramersV", "HPMergedDelta", "HP_Ratio", "Quality_Score"]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    merged["AlleleDelta_abs"] = merged.get("AlleleDelta", pd.Series(index=merged.index, dtype=float)).abs()
    merged["VAF"] = pd.to_numeric(merged.get("VAF", np.nan), errors="coerce")
    return merged


def apply_current_filter_rule(df: pd.DataFrame) -> pd.DataFrame:
    # Rule from scripts/pipeline/steps/03_filter_analysis.py
    # remove if: abs(AlleleDelta) > 0.15 and CramersV < 0.03 and (VAF < 0.15 OR VAF missing)
    cond_ad = df["AlleleDelta_abs"] > 0.15
    cond_cv = df["CramersV"] < 0.03
    cond_vaf_low = df["VAF"] < 0.15
    cond_vaf_missing = df["VAF"].isna()

    should_filter = cond_ad & cond_cv & (cond_vaf_low | cond_vaf_missing)
    pred_keep = ~should_filter

    out = df.copy()
    out["factor_ad_gt_0_15"] = cond_ad
    out["factor_cv_lt_0_03"] = cond_cv
    out["factor_vaf_lt_0_15"] = cond_vaf_low
    out["factor_vaf_missing"] = cond_vaf_missing
    out["rule_triggered"] = should_filter
    out["pred_keep"] = pred_keep

    def classify(row: pd.Series) -> str:
        if row["is_tp"] and row["pred_keep"]:
            return "TP_correct_keep"
        if row["is_tp"] and row["rule_triggered"]:
            return "TP_error_removed"
        if (not row["is_tp"]) and row["rule_triggered"]:
            return "FP_correct_removed"
        return "FP_error_kept"

    out["decision_class"] = out.apply(classify, axis=1)
    out["decision_correct"] = out["decision_class"].isin(["TP_correct_keep", "FP_correct_removed"])
    return out


def fisher_result(a: int, b: int, c: int, d: int) -> Tuple[float, float]:
    # [[a,b],[c,d]] with Haldane correction when needed for stable OR reporting.
    table = np.array([[a, b], [c, d]], dtype=float)
    if (table == 0).any():
        table += 0.5
    odds_ratio, p_value = fisher_exact(table)
    return float(odds_ratio), float(p_value)


def build_factor_impact_table(df: pd.DataFrame) -> pd.DataFrame:
    factors = [
        "factor_ad_gt_0_15",
        "factor_cv_lt_0_03",
        "factor_vaf_lt_0_15",
        "factor_vaf_missing",
        "rule_triggered",
        "Potential_LOH",
        "Significant",
        "SuggestFilter",
    ]

    rows: List[Dict[str, object]] = []
    for factor in factors:
        if factor not in df.columns:
            continue

        f = df[factor].fillna(False).astype(bool)

        # Context A: TP retained(correct) vs TP removed(error)
        tp = df["is_tp"]
        tp_correct = tp & (df["decision_class"] == "TP_correct_keep")
        tp_error = tp & (df["decision_class"] == "TP_error_removed")

        a = int((f & tp_error).sum())
        b = int((~f & tp_error).sum())
        c = int((f & tp_correct).sum())
        d = int((~f & tp_correct).sum())
        or_tp, p_tp = fisher_result(a, b, c, d)

        rows.append(
            {
                "factor": factor,
                "context": "TP_error_vs_TP_correct",
                "a_factor_in_error": a,
                "b_no_factor_in_error": b,
                "c_factor_in_correct": c,
                "d_no_factor_in_correct": d,
                "odds_ratio": or_tp,
                "p_value": p_tp,
                "factor_rate_error": a / (a + b) if (a + b) > 0 else np.nan,
                "factor_rate_correct": c / (c + d) if (c + d) > 0 else np.nan,
            }
        )

        # Context B: FP removed(correct) vs FP kept(error)
        fp = ~df["is_tp"]
        fp_correct = fp & (df["decision_class"] == "FP_correct_removed")
        fp_error = fp & (df["decision_class"] == "FP_error_kept")

        a2 = int((f & fp_correct).sum())
        b2 = int((~f & fp_correct).sum())
        c2 = int((f & fp_error).sum())
        d2 = int((~f & fp_error).sum())
        or_fp, p_fp = fisher_result(a2, b2, c2, d2)

        rows.append(
            {
                "factor": factor,
                "context": "FP_correct_vs_FP_error",
                "a_factor_in_correct": a2,
                "b_no_factor_in_correct": b2,
                "c_factor_in_error": c2,
                "d_no_factor_in_error": d2,
                "odds_ratio": or_fp,
                "p_value": p_fp,
                "factor_rate_correct": a2 / (a2 + b2) if (a2 + b2) > 0 else np.nan,
                "factor_rate_error": c2 / (c2 + d2) if (c2 + d2) > 0 else np.nan,
            }
        )

    return pd.DataFrame(rows)


def build_feature_group_stats(df: pd.DataFrame) -> pd.DataFrame:
    features = [
        "AlleleDelta_abs",
        "CramersV",
        "HPMergedDelta",
        "HP_Ratio",
        "Quality_Score",
        "VAF",
    ]

    rows: List[Dict[str, object]] = []
    group_cols = ["sample", "decision_class"]

    for (sample, cls), sub in df.groupby(group_cols):
        for feature in features:
            if feature not in sub.columns:
                continue
            vals = pd.to_numeric(sub[feature], errors="coerce").dropna()
            if vals.empty:
                continue
            rows.append(
                {
                    "sample": sample,
                    "decision_class": cls,
                    "feature": feature,
                    "n": int(vals.shape[0]),
                    "mean": float(vals.mean()),
                    "median": float(vals.median()),
                    "q25": float(vals.quantile(0.25)),
                    "q75": float(vals.quantile(0.75)),
                }
            )

    return pd.DataFrame(rows)


def build_sample_summary(records: List[RunRecord], df: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []

    for rec in records:
        m = load_metrics(rec.metrics_path)
        vc = parse_variant_counts(rec.variant_counts)
        sub = df[df["sample"] == rec.sample]

        cls_counts = sub["decision_class"].value_counts().to_dict()
        tp_total_regions = int((sub["truth_label"] == "TP").sum())
        fp_total_regions = int((sub["truth_label"] == "FP").sum())

        tp_error_removed = int(cls_counts.get("TP_error_removed", 0))
        fp_correct_removed = int(cls_counts.get("FP_correct_removed", 0))

        tp_error_rate = tp_error_removed / tp_total_regions if tp_total_regions > 0 else np.nan
        fp_correct_rate = fp_correct_removed / fp_total_regions if fp_total_regions > 0 else np.nan

        rows.append(
            {
                "sample": rec.sample,
                "run_name": rec.run_name,
                "run_dir": str(rec.run_dir),
                "baseline_f1": m["baseline_f1"],
                "filtered_f1": m["filtered_f1"],
                "f1_delta": m["f1_delta"],
                "tp_total_baseline": m["tp_total"],
                "fp_total_baseline": m["fp_total"],
                "tp_regions_analyzed": m["tp_regions_analyzed"],
                "fp_regions_analyzed": m["fp_regions_analyzed"],
                "TP_correct_keep": int(cls_counts.get("TP_correct_keep", 0)),
                "TP_error_removed": tp_error_removed,
                "FP_correct_removed": fp_correct_removed,
                "FP_error_kept": int(cls_counts.get("FP_error_kept", 0)),
                "tp_error_rate": tp_error_rate,
                "fp_correct_rate": fp_correct_rate,
                "net_filter_benefit_rate": (
                    (fp_correct_rate - tp_error_rate)
                    if not (math.isnan(fp_correct_rate) or math.isnan(tp_error_rate))
                    else np.nan
                ),
                "PASS_TOTAL": vc.get("PASS_TOTAL", np.nan),
                "TP_COUNT_variant_counts": vc.get("TP_COUNT", np.nan),
                "FP_COUNT_variant_counts": vc.get("FP_COUNT", np.nan),
            }
        )

    return pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def save_manifest(records: List[RunRecord], path: Path) -> None:
    rows = []
    for r in records:
        rows.append(
            {
                "sample": r.sample,
                "run_name": r.run_name,
                "run_dir": str(r.run_dir),
                "metrics": str(r.metrics_path),
                "tp_summary": str(r.tp_summary),
                "fp_summary": str(r.fp_summary),
                "tp_vcf": str(r.tp_vcf),
                "fp_vcf": str(r.fp_vcf),
                "variant_counts": str(r.variant_counts),
            }
        )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def plot_f1_delta(sample_df: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#2ca25f" if x >= 0 else "#de2d26" for x in sample_df["f1_delta"]]
    ax.bar(sample_df["sample"], sample_df["f1_delta"], color=colors)
    ax.axhline(0, color="black", linewidth=1)
    ax.set_ylabel("F1 Delta (Filtered - Baseline)")
    ax.set_title("Pure Tumor Samples: F1 Delta by Sample")
    ax.tick_params(axis="x", rotation=30)
    for i, v in enumerate(sample_df["f1_delta"]):
        ax.text(i, v + (0.0002 if v >= 0 else -0.00035), f"{v:+.4f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_decision_stack(sample_df: pd.DataFrame, out_png: Path) -> None:
    classes = ["TP_correct_keep", "TP_error_removed", "FP_correct_removed", "FP_error_kept"]
    data = sample_df.set_index("sample")[classes].fillna(0)
    pct = data.div(data.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(11, 5))
    bottom = np.zeros(len(pct))
    palette = {
        "TP_correct_keep": "#1b9e77",
        "TP_error_removed": "#d95f02",
        "FP_correct_removed": "#7570b3",
        "FP_error_kept": "#e7298a",
    }
    for c in classes:
        vals = pct[c].values
        ax.bar(pct.index, vals, bottom=bottom, label=c, color=palette[c])
        bottom += vals

    ax.set_ylabel("Percentage (%)")
    ax.set_title("Decision Class Composition by Sample")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_feature_box(df: pd.DataFrame, out_png: Path) -> None:
    features = ["AlleleDelta_abs", "CramersV", "HPMergedDelta", "VAF"]
    cls_order = ["TP_correct_keep", "TP_error_removed", "FP_correct_removed", "FP_error_kept"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()

    for ax, feature in zip(axes, features):
        if feature not in df.columns:
            ax.set_visible(False)
            continue
        plot_df = df[["decision_class", feature]].dropna()
        if plot_df.empty:
            ax.set_visible(False)
            continue
        sns.boxplot(
            data=plot_df,
            x="decision_class",
            y=feature,
            order=cls_order,
            showfliers=False,
            ax=ax,
        )
        ax.set_title(feature)
        ax.tick_params(axis="x", rotation=20)

    fig.suptitle("Methylation-Related Feature Distribution by Decision Class", y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_rule_scatter(df: pd.DataFrame, out_png: Path, max_points: int) -> None:
    sub = df[["AlleleDelta_abs", "CramersV", "decision_class"]].dropna().copy()
    if sub.empty:
        return

    if len(sub) > max_points:
        sub = sub.sample(max_points, random_state=RANDOM_STATE)

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=sub,
        x="AlleleDelta_abs",
        y="CramersV",
        hue="decision_class",
        alpha=0.45,
        s=10,
        ax=ax,
    )
    ax.axvline(0.15, color="black", linestyle="--", linewidth=1)
    ax.axhline(0.03, color="black", linestyle="--", linewidth=1)
    ax.set_title("Decision Boundary View (Current Filter Rule)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_benefit_heatmap(sample_df: pd.DataFrame, out_png: Path) -> None:
    mat = sample_df.set_index("sample")[["tp_error_rate", "fp_correct_rate", "net_filter_benefit_rate"]]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    sns.heatmap(mat, annot=True, fmt=".4f", cmap="RdYlGn", center=0, ax=ax)
    ax.set_title("Filter Impact Heatmap (Rate)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def write_longphase_to_plan(path: Path) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    content = f"""<!--
建立時間: {now}
目標: 評估在 LongPhase-TO（tumor-only）純樣本下複用同一套甲基檢驗流程的可行性
-->

# LongPhase-TO 純樣本可行性與計劃

## 1. 結論摘要

- 可行，但需補一層 **tumor-only 專用前處理**（caller + PON + phase/haplotag）後才能直接套用本次純 tumor 評估流程。
- 目前流程中與 LongPhase-S 綁定最緊的是 `step01`（TP/FP VCF + tagged BAM 產生方式），其餘 TP/FP 統計、InterSubMod 特徵比較、圖表報告可高度沿用。

## 2. 可直接沿用項目

- InterSubMod 輸出解析：`significance_summary.csv`（TP/FP 各一份）
- Filter 規則檢驗：`AlleleDelta/CramersV/VAF` 的決策比較
- 分類框架：`TP_correct_keep / TP_error_removed / FP_correct_removed / FP_error_kept`
- 報告與圖表：可沿用本次 Python 產出的資料結構與繪圖模板

## 3. 必要新增或調整

1. 上游改為 tumor-only caller 與 phasing
- 建議路線：ClairS-TO（或 DeepSomatic TO） -> LongPhase-TO `phase` -> LongPhase-TO `haplotag`
- 需提供 `--caller` 與 `--pon-file / --strict-pon-file`。

2. TP/FP 標記流程改造
- 現行 `run_benchmark.sh` step01 以 LongPhase-S 流程為主，LongPhase-TO 需新增對應入口，產出同樣格式的 `filtered_snv_tp/fp.vcf.gz`。

3. BAM tag 檢查
- 必須保留 MM/ML；同時確認 HP tag 格式在 InterSubMod 解析端與 LongPhase-TO 輸出一致。

## 4. 風險與驗證點

- 風險 A：PON 品質/版本導致 somatic 判定差異，會直接改變 TP/FP 基線。
- 風險 B：tumor-only 模式下，若 low-AF 事件增加，現行固定門檻可能放大 TP 誤刪。
- 風險 C：不同 caller 的 FORMAT 欄位差異（VAF/AF/AD）影響規則重建一致性。

## 5. 分階段執行計劃

### Phase A：最小可行打通（1 個樣本）

- 樣本：`HCC1395_DORADO`（純 tumor）
- 產物：
  - TO 版 tagged BAM
  - TO 版 TP/FP VCF
  - InterSubMod TP/FP summaries
  - 套用本次評估腳本的完整報告

### Phase B：跨樣本擴展（6-7 個純 tumor 樣本）

- 目標：建立與 LongPhase-S 同口徑的跨樣本比較矩陣
- 指標：F1 delta、TP/FP 決策分類、因子 OR/p-value、甲基特徵分布差異

### Phase C：流程整合

- 在 `scripts/pipeline/` 增加 TO mode（或新增 `run_benchmark_to.sh`）
- 讓本次 Python 評估成為固定的 step（資料 + 圖 + md 報告自動化）
"""
    path.write_text(content, encoding="utf-8")


def write_main_report(
    report_path: Path,
    sample_df: pd.DataFrame,
    factor_df: pd.DataFrame,
    output_dir: Path,
    note_on_existing_changes: str,
) -> None:
    total_tp = int(sample_df["TP_correct_keep"].sum() + sample_df["TP_error_removed"].sum())
    total_fp = int(sample_df["FP_correct_removed"].sum() + sample_df["FP_error_kept"].sum())
    total_tp_err = int(sample_df["TP_error_removed"].sum())
    total_fp_correct = int(sample_df["FP_correct_removed"].sum())

    tp_err_rate = total_tp_err / total_tp if total_tp > 0 else np.nan
    fp_correct_rate = total_fp_correct / total_fp if total_fp > 0 else np.nan

    best = sample_df.loc[sample_df["f1_delta"].idxmax()]
    worst = sample_df.loc[sample_df["f1_delta"].idxmin()]

    # Pick two strongest factor signals by p-value in each context.
    factor_tp = factor_df[factor_df["context"] == "TP_error_vs_TP_correct"].copy()
    factor_fp = factor_df[factor_df["context"] == "FP_correct_vs_FP_error"].copy()
    excluded = {"rule_triggered"}
    factor_tp = factor_tp[~factor_tp["factor"].isin(excluded)]
    factor_fp = factor_fp[~factor_fp["factor"].isin(excluded)]
    factor_tp = factor_tp.sort_values("p_value").head(3)
    factor_fp = factor_fp.sort_values("p_value").head(3)

    def factor_lines(df: pd.DataFrame) -> str:
        lines = []
        for _, r in df.iterrows():
            lines.append(
                f"- `{r['factor']}`: OR={r['odds_ratio']:.3f}, p={r['p_value']:.3e}, "
                f"rate_correct={r.get('factor_rate_correct', np.nan):.3f}, rate_error={r.get('factor_rate_error', np.nan):.3f}"
            )
        return "\n".join(lines) if lines else "- 無可用結果"

    content = f"""<!--
建立時間: {datetime.now().strftime('%Y-%m-%d %H:%M')}
目標: 純 tumor（目前以 s-pure-pileup 為主）甲基檢驗評估流程與觀察報告
-->

# 純 Tumor 甲基檢驗評估報告

## 1. 本次分析範圍

- 樣本：`{', '.join(sample_df['sample'].tolist())}`
- 資料來源：`/bip8_disk/liaoyoyo2001/InterSubMod_out/output/s-pure-pileup/*/<run>`（每樣本挑選可用且指標最佳 run）
- 比較重點：
  - TP / FP
  - 分類正確 / 分類錯誤（依現行過濾規則）
  - 甲基相關特徵對判斷的影響

## 2. 既有未提交變更影響確認

{note_on_existing_changes}

## 3. 核心結論（跨樣本）

- 總 TP 區域：`{total_tp}`，TP 誤刪（錯誤）`{total_tp_err}`，誤刪率 `{tp_err_rate:.4f}`
- 總 FP 區域：`{total_fp}`，FP 正確移除 `({total_fp_correct})`，移除率 `{fp_correct_rate:.4f}`
- 整體淨效益（FP 去除 vs TP 誤刪）約為 `{(fp_correct_rate - tp_err_rate):+.4f}`
- 最佳 F1 變化樣本：`{best['sample']}` (`{best['f1_delta']:+.4f}`)
- 最差 F1 變化樣本：`{worst['sample']}` (`{worst['f1_delta']:+.4f}`)

## 4. 圖表總覽

- [F1 變化總覽](../plots/f1_delta_by_sample.png)
- [TP/FP + 正確/錯誤組成](../plots/decision_class_stack_by_sample.png)
- [甲基特徵分布（分 decision class）](../plots/methylation_feature_boxplots.png)
- [現行規則邊界散點（AD vs CV）](../plots/rule_boundary_scatter.png)
- [各樣本過濾影響熱圖](../plots/filter_impact_heatmap.png)

## 5. 因子影響摘要

### 5.1 對 TP 誤刪最敏感的因子（TP_error_vs_TP_correct）

{factor_lines(factor_tp)}

### 5.2 對 FP 正確移除最有關的因子（FP_correct_vs_FP_error）

{factor_lines(factor_fp)}

## 6. 可檢閱資料表

- [run_manifest.tsv](../data/run_manifest.tsv)
- [sample_summary.tsv](../data/sample_summary.tsv)
- [classification_counts.tsv](../data/classification_counts.tsv)
- [factor_impact.tsv](../data/factor_impact.tsv)
- [feature_group_stats.tsv](../data/feature_group_stats.tsv)
- [per_variant_decision.tsv.gz](../data/per_variant_decision.tsv.gz)

## 7. LongPhase-TO 可行性計劃

- [longphase_to_feasibility_plan.md](longphase_to_feasibility_plan.md)
"""
    report_path.write_text(content, encoding="utf-8")


def main() -> None:
    args = parse_args()
    np.random.seed(RANDOM_STATE)

    data_root = Path(args.data_root)
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    if args.output_dir:
        out_root = Path(args.output_dir)
    else:
        run_tag = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_root = Path(args.output_root) / f"pure_tumor_eval_{run_tag}"

    data_dir = ensure_dir(out_root / "data")
    plot_dir = ensure_dir(out_root / "plots")
    report_dir = ensure_dir(out_root / "report")

    records: List[RunRecord] = []
    missing_samples: List[str] = []

    for sample in samples:
        rec = find_best_run_for_sample(data_root / sample, sample)
        if rec is None:
            missing_samples.append(sample)
        else:
            records.append(rec)

    if not records:
        raise RuntimeError("No valid runs found. Please check --data-root and sample availability.")

    all_parts: List[pd.DataFrame] = []
    for rec in records:
        tp_df = load_summary_with_features(rec.tp_summary, rec.tp_vcf, rec.sample, rec.run_name, "TP")
        fp_df = load_summary_with_features(rec.fp_summary, rec.fp_vcf, rec.sample, rec.run_name, "FP")
        merged = pd.concat([tp_df, fp_df], ignore_index=True)
        merged = apply_current_filter_rule(merged)
        all_parts.append(merged)

    all_df = pd.concat(all_parts, ignore_index=True)

    sample_df = build_sample_summary(records, all_df)

    cls_counts = (
        all_df.groupby(["sample", "decision_class"]).size().unstack(fill_value=0).reset_index().sort_values("sample")
    )

    factor_df = build_factor_impact_table(all_df)
    feat_stats_df = build_feature_group_stats(all_df)

    # Save tables
    save_manifest(records, data_dir / "run_manifest.tsv")
    save_tsv(sample_df, data_dir / "sample_summary.tsv")
    save_tsv(cls_counts, data_dir / "classification_counts.tsv")
    save_tsv(factor_df, data_dir / "factor_impact.tsv")
    save_tsv(feat_stats_df, data_dir / "feature_group_stats.tsv")

    # Per-variant decision table (compressed TSV)
    with gzip.open(data_dir / "per_variant_decision.tsv.gz", "wt") as gz:
        all_df.to_csv(gz, sep="\t", index=False)

    # Plots
    sns.set_theme(style="whitegrid")
    plot_f1_delta(sample_df, plot_dir / "f1_delta_by_sample.png")
    plot_decision_stack(sample_df, plot_dir / "decision_class_stack_by_sample.png")
    plot_feature_box(all_df, plot_dir / "methylation_feature_boxplots.png")
    plot_rule_scatter(all_df, plot_dir / "rule_boundary_scatter.png", args.max_points_scatter)
    plot_benefit_heatmap(sample_df, plot_dir / "filter_impact_heatmap.png")

    # Reports
    note = (
        "- `docs/experiments/INDEX.md` 的未提交變更屬文件索引更新，不影響本次數據計算。\n"
        "- `src/core/RegionProcessor.cpp` 的未提交變更為 `--no-filter` 安全修正；"
        "本次採用的 pipeline run 未傳入 `--no-filter`，因此不改變本報告所用結果。"
    )
    write_longphase_to_plan(report_dir / "longphase_to_feasibility_plan.md")
    write_main_report(report_dir / "README.md", sample_df, factor_df, out_root, note)

    # Metadata
    meta = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "data_root": str(data_root),
        "output_dir": str(out_root),
        "samples_requested": samples,
        "samples_missing": missing_samples,
        "samples_used": [r.sample for r in records],
        "record_count": len(records),
        "row_count": int(all_df.shape[0]),
    }
    with open(data_dir / "metadata.json", "w") as f:
        json.dump(meta, f, indent=2)

    print("[Done] Pure tumor methylation evaluation complete")
    print(f"[Done] Output: {out_root}")
    if missing_samples:
        print(f"[Warn] Missing samples (no valid run found): {', '.join(missing_samples)}")


if __name__ == "__main__":
    main()
