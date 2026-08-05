#!/usr/bin/env python3
"""
Compare TP/FP enrichment patterns using explicit legacy verification support.

Outputs:
1) Metrics table (TSV)
2) Rule performance table (TSV)
3) Enriched sample table with derived features (TSV)
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    select_legacy_view,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare TP/FP enrichment using the explicit legacy four-state verification field."
    )
    parser.add_argument(
        "--advanced-samples",
        default="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/advanced_analysis_20260119/advanced_samples.csv",
        help="Path to advanced_samples.csv",
    )
    parser.add_argument(
        "--tp-vcf",
        default="data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz",
        help="TP VCF path",
    )
    parser.add_argument(
        "--fp-vcf",
        default="data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz",
        help="FP VCF path",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=25000,
        help="Half window size used for cluster counting (bp)",
    )
    parser.add_argument(
        "--high-cluster-threshold",
        type=int,
        default=20,
        help="Threshold for HighCluster feature",
    )
    parser.add_argument(
        "--low-af-threshold",
        type=float,
        default=0.3,
        help="Threshold for LowAF feature",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix for generated TSV files (required; no historical output default)",
    )
    return parser.parse_args()


def derive_legacy_support_columns(df: pd.DataFrame):
    """Attach the two locked historical support definitions from the legacy field."""
    view = select_legacy_view(df, allow_unversioned_v1=False)
    result = df.copy()
    result["VerificationClass_Legacy_Selected"] = view.values
    result["LegacyVerificationSupport"] = view.values.isin(["Strong", "Subclone"])
    result["LegacyClusterFirstOnly"] = view.values.eq("Subclone")
    return result, view


def query_vcf_af_h(vcf_path: str) -> pd.DataFrame:
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%INFO/H\t[%AF]\n", vcf_path]
    out = subprocess.check_output(cmd, text=True)
    rows: List[Tuple[str, int, int, float]] = []
    for line in out.splitlines():
        chrom, pos, hflag, af = line.split("\t")
        rows.append((chrom, int(pos), 1 if hflag == "1" else 0, float(af)))
    return pd.DataFrame(rows, columns=["Chr", "Pos", "HFlag", "AF"])


def build_position_index(df: pd.DataFrame) -> Dict[str, List[int]]:
    pos_map: Dict[str, List[int]] = defaultdict(list)
    for row in df.itertuples(index=False):
        pos_map[row.Chr].append(int(row.Pos))
    for chrom in pos_map:
        pos_map[chrom].sort()
    return dict(pos_map)


def count_neighbors(chrom: str, pos: int, pos_map: Dict[str, List[int]], window: int) -> int:
    arr = pos_map.get(chrom, [])
    left = bisect_left(arr, pos - window)
    right = bisect_right(arr, pos + window)
    return right - left


def fisher_stats(df: pd.DataFrame, feature_col: str, label_col: str = "IsTP") -> Dict[str, float]:
    cond = df[feature_col].astype(bool)
    y = df[label_col].astype(bool)

    a = int(((cond) & (y)).sum())
    b = int(((cond) & (~y)).sum())
    c = int(((~cond) & (y)).sum())
    d = int(((~cond) & (~y)).sum())

    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")

    tp_rate_cond = float(a / (a + b)) if (a + b) > 0 else float("nan")
    tp_rate_non = float(c / (c + d)) if (c + d) > 0 else float("nan")

    return {
        "a_cond_tp": a,
        "b_cond_fp": b,
        "c_noncond_tp": c,
        "d_noncond_fp": d,
        "odds_ratio": float(odds_ratio),
        "p_value": float(p_value),
        "tp_rate_cond": tp_rate_cond,
        "tp_rate_non": tp_rate_non,
    }


def f1_stats(y_true: Iterable[bool], y_pred: Iterable[bool]) -> Dict[str, float]:
    truth = np.asarray(list(y_true), dtype=bool)
    pred = np.asarray(list(y_pred), dtype=bool)

    tp = int((truth & pred).sum())
    fp = int((~truth & pred).sum())
    fn = int((truth & ~pred).sum())
    tn = int((~truth & ~pred).sum())

    precision = float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0
    recall = float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0
    f1 = float((2 * precision * recall) / (precision + recall)) if (precision + recall) > 0 else 0.0
    specificity = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
    accuracy = float((tp + tn) / (tp + tn + fp + fn)) if (tp + tn + fp + fn) > 0 else 0.0

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "accuracy": accuracy,
        "f1": f1,
    }


def main() -> None:
    args = parse_args()
    out_prefix = Path(args.out_prefix)

    advanced_df = pd.read_csv(args.advanced_samples)
    try:
        advanced_df, legacy_view = derive_legacy_support_columns(advanced_df)
    except SchemaContractError as exc:
        print(f"ERROR: legacy verification schema contract failed: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc

    advanced_df["IsTP"] = advanced_df["Category"].astype(str).str.startswith("TP")

    tp_df = query_vcf_af_h(args.tp_vcf)
    fp_df = query_vcf_af_h(args.fp_vcf)

    tp_map = {(row.Chr, int(row.Pos)): (int(row.HFlag), float(row.AF)) for row in tp_df.itertuples(index=False)}
    fp_map = {(row.Chr, int(row.Pos)): (int(row.HFlag), float(row.AF)) for row in fp_df.itertuples(index=False)}

    h_values: List[int] = []
    af_values: List[float] = []
    missing = 0
    for row in advanced_df.itertuples(index=False):
        key = (row.Chr, int(row.Pos))
        source_map = tp_map if str(row.Category).startswith("TP") else fp_map
        if key in source_map:
            hflag, af = source_map[key]
        else:
            hflag, af = 0, np.nan
            missing += 1
        h_values.append(hflag)
        af_values.append(af)

    advanced_df["HFlag"] = h_values
    advanced_df["AF"] = af_values

    tp_pos_map = build_position_index(tp_df)
    fp_pos_map = build_position_index(fp_df)

    cluster_values: List[int] = []
    for row in advanced_df.itertuples(index=False):
        pos_map = tp_pos_map if str(row.Category).startswith("TP") else fp_pos_map
        cluster_values.append(count_neighbors(str(row.Chr), int(row.Pos), pos_map, args.window))
    advanced_df["ClusterCount50kb"] = cluster_values

    advanced_df["LowAF"] = advanced_df["AF"] < args.low_af_threshold
    advanced_df["HighCluster"] = advanced_df["ClusterCount50kb"] > args.high_cluster_threshold
    advanced_df["LOH_like"] = advanced_df["Potential_LOH"].fillna(False) | (advanced_df["HFlag"] == 1)

    lowaf_df = advanced_df[advanced_df["LowAF"]].copy()

    metrics_rows: List[Dict[str, object]] = []

    metrics_rows.append({"group": "meta", "metric": "sample_size_total", "value": int(len(advanced_df))})
    metrics_rows.append({"group": "meta", "metric": "sample_size_tp", "value": int(advanced_df["IsTP"].sum())})
    metrics_rows.append({"group": "meta", "metric": "sample_size_fp", "value": int((~advanced_df["IsTP"]).sum())})
    metrics_rows.append({"group": "meta", "metric": "sample_size_lowAF", "value": int(len(lowaf_df))})
    metrics_rows.append({"group": "meta", "metric": "missing_af_h_lookup", "value": int(missing)})
    metrics_rows.append(
        {"group": "meta", "metric": "verification_selection_field", "value": legacy_view.field}
    )
    metrics_rows.append(
        {"group": "meta", "metric": "verification_schema_status", "value": legacy_view.schema_status}
    )

    for feat in ["LOH_like", "HighCluster", "LowAF"]:
        stats = fisher_stats(advanced_df, feat)
        for key, value in stats.items():
            metrics_rows.append({"group": feat, "metric": key, "value": value})

    lowaf_support_stats = fisher_stats(lowaf_df, "LegacyVerificationSupport")
    for key, value in lowaf_support_stats.items():
        metrics_rows.append({"group": "lowAF_LegacyVerificationSupport", "metric": key, "value": value})

    lowaf_cluster_first_stats = fisher_stats(lowaf_df, "LegacyClusterFirstOnly")
    for key, value in lowaf_cluster_first_stats.items():
        metrics_rows.append({"group": "lowAF_LegacyClusterFirstOnly", "metric": key, "value": value})

    lowaf_baseline_tp_rate = float(lowaf_df["IsTP"].mean()) if len(lowaf_df) > 0 else float("nan")
    lowaf_support_tp_rate = (
        float(lowaf_df.loc[lowaf_df["LegacyVerificationSupport"], "IsTP"].mean())
        if int(lowaf_df["LegacyVerificationSupport"].sum()) > 0
        else float("nan")
    )
    lowaf_nonsupport_tp_rate = (
        float(lowaf_df.loc[~lowaf_df["LegacyVerificationSupport"], "IsTP"].mean())
        if int((~lowaf_df["LegacyVerificationSupport"]).sum()) > 0
        else float("nan")
    )
    metrics_rows.append({"group": "lowAF_rates", "metric": "baseline_tp_rate", "value": lowaf_baseline_tp_rate})
    metrics_rows.append(
        {"group": "lowAF_rates", "metric": "legacy_verification_support_tp_rate", "value": lowaf_support_tp_rate}
    )
    metrics_rows.append(
        {
            "group": "lowAF_rates",
            "metric": "non_legacy_verification_support_tp_rate",
            "value": lowaf_nonsupport_tp_rate,
        }
    )

    metrics_df = pd.DataFrame(metrics_rows)

    rule_rows: List[Dict[str, object]] = []
    y = advanced_df["IsTP"].astype(bool)

    rule_1 = advanced_df["AF"] >= args.low_af_threshold
    rule_2 = (advanced_df["AF"] >= args.low_af_threshold) | advanced_df["LegacyClusterFirstOnly"]
    rule_3 = advanced_df["Significant"].fillna(False).astype(bool)
    rule_4 = (advanced_df["AF"] >= args.low_af_threshold) | advanced_df["LegacyVerificationSupport"]

    for name, pred in [
        (f"AF>={args.low_af_threshold}", rule_1),
        (f"AF>={args.low_af_threshold}_OR_LegacyClusterFirstOnly", rule_2),
        ("SignificantOnly", rule_3),
        (f"AF>={args.low_af_threshold}_OR_LegacyVerificationSupport", rule_4),
    ]:
        stats = f1_stats(y, pred)
        stats["rule"] = name
        rule_rows.append(stats)

    rules_df = pd.DataFrame(rule_rows)[
        ["rule", "tp", "fp", "fn", "tn", "precision", "recall", "specificity", "accuracy", "f1"]
    ]

    metrics_out = out_prefix.with_suffix(".tsv")
    rules_out = out_prefix.parent / f"{out_prefix.name}_rules.tsv"
    samples_out = out_prefix.parent / f"{out_prefix.name}_samples.tsv"

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(metrics_out, sep="\t", index=False)
    rules_df.to_csv(rules_out, sep="\t", index=False)
    advanced_df.to_csv(samples_out, sep="\t", index=False)

    print(f"[OK] metrics: {metrics_out}")
    print(f"[OK] rules:   {rules_out}")
    print(f"[OK] samples: {samples_out}")


if __name__ == "__main__":
    main()
