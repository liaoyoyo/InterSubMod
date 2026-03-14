#!/usr/bin/env python3
"""
診斷純 tumor 規則提升幅度偏小的根因：
1. 各樣本理論上限（若可移除全部 FP、且不誤刪 TP）
2. 現行規則命中品質（trigger precision / F1 delta）
3. 主要特徵對 FP/TP 可分性的 AUC
4. 單一特徵最佳門檻可達到的 F1 delta（僅作上界參考）
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


DEFAULT_FEATURES = [
    "VAF",
    "AlleleDelta_abs",
    "CramersV",
    "HPMergedDelta",
    "Quality_Score",
]


def f1_score(tp: int, fp: int, truth_total: int) -> float:
    fn = max(truth_total - tp, 0)
    denom = 2 * tp + fp + fn
    if denom == 0:
        return 0.0
    return (2 * tp) / denom


def load_truth_and_baseline(manifest_path: Path) -> Tuple[Dict[str, int], Dict[str, Tuple[int, int]]]:
    manifest = pd.read_csv(manifest_path, sep="\t")
    truth_total: Dict[str, int] = {}
    baseline_counts: Dict[str, Tuple[int, int]] = {}

    for _, row in manifest.iterrows():
        sample = row["sample"]
        metrics_path = Path(row["metrics"])
        metrics = json.loads(metrics_path.read_text())
        truth_total[sample] = int(metrics["truth_total"])
        baseline_counts[sample] = (int(metrics["baseline"]["tp"]), int(metrics["baseline"]["fp"]))

    return truth_total, baseline_counts


def normalize_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().eq("true")


def evaluate_features(
    sample_df: pd.DataFrame,
    sample: str,
    truth_total: int,
    baseline_tp: int,
    baseline_fp: int,
    features: List[str],
) -> List[dict]:
    y_fp = (~sample_df["is_tp"]).astype(int).values
    baseline_f1 = f1_score(baseline_tp, baseline_fp, truth_total)
    rows: List[dict] = []

    for feature in features:
        feature_values = pd.to_numeric(sample_df[feature], errors="coerce")
        valid = ~feature_values.isna()
        x = feature_values[valid].values
        y = y_fp[valid.values]

        result = {
            "sample": sample,
            "feature": feature,
            "auc_fp": np.nan,
            "best_oriented_auc": np.nan,
            "direction": "na",
            "best_single_feature_f1_delta": np.nan,
            "best_single_feature_fp_removed": 0,
            "best_single_feature_tp_removed": 0,
            "best_single_feature_rule": "na",
        }

        if len(np.unique(y)) < 2 or len(np.unique(x)) < 2:
            rows.append(result)
            continue

        auc_raw = roc_auc_score(y, x)
        if auc_raw >= 0.5:
            direction = "higher->FP"
            oriented_auc = auc_raw
        else:
            direction = "lower->FP"
            oriented_auc = 1 - auc_raw

        values = np.unique(x)
        candidates = np.quantile(values, np.linspace(0.005, 0.995, 400))

        best_delta = -1e18
        best_rule = ""
        best_fp_removed = 0
        best_tp_removed = 0

        for threshold in candidates:
            if direction == "lower->FP":
                remove = x <= threshold
                rule = f"{feature}<={threshold:.6g}"
            else:
                remove = x >= threshold
                rule = f"{feature}>={threshold:.6g}"

            tp_removed = int(((y == 0) & remove).sum())
            fp_removed = int(((y == 1) & remove).sum())
            new_f1 = f1_score(baseline_tp - tp_removed, baseline_fp - fp_removed, truth_total)
            delta = new_f1 - baseline_f1

            if delta > best_delta:
                best_delta = float(delta)
                best_rule = rule
                best_fp_removed = fp_removed
                best_tp_removed = tp_removed

        result.update(
            {
                "auc_fp": float(auc_raw),
                "best_oriented_auc": float(oriented_auc),
                "direction": direction,
                "best_single_feature_f1_delta": best_delta,
                "best_single_feature_fp_removed": best_fp_removed,
                "best_single_feature_tp_removed": best_tp_removed,
                "best_single_feature_rule": best_rule,
            }
        )
        rows.append(result)

    return rows


def evaluate_global_auc(df: pd.DataFrame, features: List[str]) -> pd.DataFrame:
    y_fp = (~df["is_tp"]).astype(int)
    rows = []
    for feature in features:
        x = pd.to_numeric(df[feature], errors="coerce")
        valid = ~x.isna()
        xv = x[valid].values
        yv = y_fp[valid].values
        if len(np.unique(yv)) < 2 or len(np.unique(xv)) < 2:
            rows.append(
                {
                    "feature": feature,
                    "global_auc_fp": np.nan,
                    "global_best_oriented_auc": np.nan,
                    "global_direction": "na",
                }
            )
            continue

        auc_raw = roc_auc_score(yv, xv)
        if auc_raw >= 0.5:
            best = auc_raw
            direction = "higher->FP"
        else:
            best = 1 - auc_raw
            direction = "lower->FP"

        rows.append(
            {
                "feature": feature,
                "global_auc_fp": float(auc_raw),
                "global_best_oriented_auc": float(best),
                "global_direction": direction,
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Small-gain root cause diagnosis.")
    parser.add_argument(
        "--input",
        required=True,
        help="Path to per_variant_decision.tsv(.gz)",
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to run_manifest.tsv",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for diagnosis tables",
    )
    parser.add_argument(
        "--features",
        nargs="+",
        default=DEFAULT_FEATURES,
        help="Feature columns to evaluate",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    manifest_path = Path(args.manifest)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    df["is_tp"] = normalize_bool(df["is_tp"])
    df["rule_triggered"] = normalize_bool(df["rule_triggered"])

    truth_total_by_sample, baseline_counts_by_sample = load_truth_and_baseline(manifest_path)

    headroom_rows: List[dict] = []
    current_rule_rows: List[dict] = []
    feature_rows: List[dict] = []

    for sample, group in df.groupby("sample"):
        truth_total = truth_total_by_sample[sample]
        baseline_tp, baseline_fp = baseline_counts_by_sample[sample]
        baseline_f1 = f1_score(baseline_tp, baseline_fp, truth_total)
        upper_f1 = f1_score(baseline_tp, 0, truth_total)

        headroom_rows.append(
            {
                "sample": sample,
                "baseline_f1_exact": baseline_f1,
                "f1_upper_remove_all_fp": upper_f1,
                "max_delta_if_zero_fp": upper_f1 - baseline_f1,
                "baseline_tp": baseline_tp,
                "baseline_fp": baseline_fp,
                "fp_to_tp_ratio": baseline_fp / max(baseline_tp, 1),
            }
        )

        rule = group["rule_triggered"]
        tp_removed = int((group["is_tp"] & rule).sum())
        fp_removed = int(((~group["is_tp"]) & rule).sum())
        triggered_total = tp_removed + fp_removed
        triggered_f1 = f1_score(baseline_tp - tp_removed, baseline_fp - fp_removed, truth_total)

        current_rule_rows.append(
            {
                "sample": sample,
                "triggered_total": triggered_total,
                "trigger_precision_fp": (fp_removed / triggered_total) if triggered_total > 0 else np.nan,
                "trigger_tp_removed": tp_removed,
                "trigger_fp_removed": fp_removed,
                "trigger_f1_delta": triggered_f1 - baseline_f1,
            }
        )

        feature_rows.extend(
            evaluate_features(
                sample_df=group,
                sample=sample,
                truth_total=truth_total,
                baseline_tp=baseline_tp,
                baseline_fp=baseline_fp,
                features=args.features,
            )
        )

    feature_auc_df = pd.DataFrame(feature_rows)
    headroom_df = pd.DataFrame(headroom_rows)
    current_rule_df = pd.DataFrame(current_rule_rows)
    global_auc_df = evaluate_global_auc(df, args.features)

    feature_auc_df.to_csv(out_dir / "small_gain_diagnosis_feature_auc.tsv", sep="\t", index=False)
    headroom_df.to_csv(out_dir / "small_gain_diagnosis_headroom.tsv", sep="\t", index=False)
    global_auc_df.to_csv(out_dir / "small_gain_diagnosis_global_auc.tsv", sep="\t", index=False)
    current_rule_df.to_csv(out_dir / "small_gain_diagnosis_current_rule_precision.tsv", sep="\t", index=False)

    print(f"[DONE] Wrote {out_dir / 'small_gain_diagnosis_feature_auc.tsv'}")
    print(f"[DONE] Wrote {out_dir / 'small_gain_diagnosis_headroom.tsv'}")
    print(f"[DONE] Wrote {out_dir / 'small_gain_diagnosis_global_auc.tsv'}")
    print(f"[DONE] Wrote {out_dir / 'small_gain_diagnosis_current_rule_precision.tsv'}")


if __name__ == "__main__":
    main()

