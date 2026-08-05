#!/usr/bin/env python3
"""Compare candidate methylation filtering rules across one or more benchmark runs."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence

import pandas as pd

from research_common import infer_platform, read_json, to_float, write_tsv_rows

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import select_current_view


RuleFunc = Callable[[Dict[str, object]], bool]


COMPARISON_FIELDS = [
    "sample",
    "platform",
    "run_dir",
    "rule_id",
    "trigger_count",
    "tp_removed",
    "fp_removed",
    "fp_tp_ratio",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "baseline_f1",
    "truth_total",
]

SUBSET_FIELDS = [
    "sample",
    "platform",
    "rule_id",
    "subset",
    "count",
    "median_qual",
    "median_vaf",
    "median_allele_delta",
    "median_cramers_v",
    "median_quality_score",
    "median_pairwise_median_dist",
    "dominant_label_top",
    "verification_class_top",
    "verification_selection_field",
    "verification_schema_status",
    "unknown_current_class_count",
]

SCHEMA_FIELDS = [
    "run_dir",
    "scope",
    "selection_field",
    "schema_status",
    "row_count",
    "unknown_current_class_count",
    "unknown_current_class_values",
]

BEST_FIELDS = [
    "sample",
    "platform",
    "run_dir",
    "baseline_f1",
    "best_rule_id",
    "best_f1",
    "best_delta_f1_vs_baseline",
]


def median(values: Sequence[float]) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", action="append", required=True, help="Benchmark run directory")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "--focus-rule",
        action="append",
        default=["old_combo_only", "current_core_ad015_vaf015"],
        help="Rule IDs to emit subset summaries for",
    )
    return parser.parse_args()


def open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def parse_vcf_features(path: Path, scope: str) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with open_vcf(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            chrom, pos, _vid, ref, alt, qual, filt, info, *rest = raw.rstrip("\n").split("\t")
            info_map: Dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info_map[key] = value
            fmt_fields: List[str] = []
            sample_fields: List[str] = []
            if len(rest) >= 2:
                fmt_fields = rest[0].split(":")
                sample_fields = rest[1].split(":")
            fmt_map = dict(zip(fmt_fields, sample_fields))
            vaf = math.nan
            for key in ("VAF", "AF"):
                if key in fmt_map:
                    vaf = to_float(fmt_map[key].split(",")[0])
                    break
                if key in info_map:
                    vaf = to_float(info_map[key].split(",")[0])
                    break

            rows[f"{chrom}:{pos}:{ref}:{alt}"] = {
                "label": scope,
                "QUAL": to_float(qual),
                "VAF": vaf,
                "FILTER": filt,
            }
    return rows


def parse_summary(
    path: Path,
    scope: str,
) -> tuple[Dict[str, Dict[str, object]], Dict[str, object]]:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    view = select_current_view(frame)
    frame = frame.copy()
    frame["_verification_class_current"] = view.values
    frame["_verification_selection_field"] = view.field
    frame["_verification_schema_status"] = view.schema_status
    rows: Dict[str, Dict[str, object]] = {}
    for row in frame.to_dict(orient="records"):
        key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
        rows[key] = {
            **row,
            "label": scope,
            "AlleleDelta": to_float(row.get("AlleleDelta")),
            "CramersV": to_float(row.get("CramersV")),
            "Quality_Score": to_float(row.get("Quality_Score")),
            "PairwiseMedianDist": to_float(row.get("PairwiseMedianDist")),
        }
    metadata = {
        "scope": scope,
        "selection_field": view.field,
        "schema_status": view.schema_status,
        "row_count": len(frame),
        "unknown_current_class_count": sum(view.unknown_counts.values()),
        "unknown_current_class_values": "; ".join(
            f"{key}:{value}" for key, value in sorted(view.unknown_counts.items())
        ),
    }
    return rows, metadata


def resolve_longphase_vcf(run_dir: Path, scope: str) -> Path:
    direct = run_dir / "longphase_s" / f"filtered_snv_{scope.lower()}.vcf.gz"
    if direct.exists():
        return direct

    parent = run_dir.parent
    candidates: List[Path] = []
    for sibling in sorted(parent.iterdir(), reverse=True):
        if not sibling.is_dir() or sibling == run_dir:
            continue
        path = sibling / "longphase_s" / f"filtered_snv_{scope.lower()}.vcf.gz"
        if path.exists():
            candidates.append(path)
    if candidates:
        return candidates[0]

    raise FileNotFoundError(f"filtered_snv_{scope.lower()}.vcf.gz not found for {run_dir}")


def build_rules() -> Dict[str, RuleFunc]:
    return {
        "old_written_rule_raw_qual": lambda row: (
            (not math.isnan(row["QUAL"]) and row["QUAL"] < 0.75)
            or (
                (not math.isnan(row["AlleleDelta"]) and row["AlleleDelta"] > 0.25)
                and (not math.isnan(row["CramersV"]) and row["CramersV"] < 0.05)
                and (not math.isnan(row["VAF"]) and row["VAF"] < 0.24)
            )
        ),
        "old_combo_only": lambda row: (
            (not math.isnan(row["AlleleDelta"]) and row["AlleleDelta"] > 0.25)
            and (not math.isnan(row["CramersV"]) and row["CramersV"] < 0.05)
            and (not math.isnan(row["VAF"]) and row["VAF"] < 0.24)
        ),
        "current_core_ad015_vaf015": lambda row: (
            (not math.isnan(row["AlleleDelta"]) and row["AlleleDelta"] > 0.15)
            and (not math.isnan(row["VAF"]) and row["VAF"] < 0.15)
        ),
        "ad025_vaf024": lambda row: (
            (not math.isnan(row["AlleleDelta"]) and row["AlleleDelta"] > 0.25)
            and (not math.isnan(row["VAF"]) and row["VAF"] < 0.24)
        ),
        "ad015_cv005_vaf015": lambda row: (
            (not math.isnan(row["AlleleDelta"]) and row["AlleleDelta"] > 0.15)
            and (not math.isnan(row["CramersV"]) and row["CramersV"] < 0.05)
            and (not math.isnan(row["VAF"]) and row["VAF"] < 0.15)
        ),
    }


def compute_metrics(tp: int, fp: int, truth_total: int) -> Dict[str, float]:
    fn = max(truth_total - tp, 0)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / truth_total if truth_total else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "fn": fn,
    }


def summarize_subset(sample: str, platform: str, rule_id: str, subset: str, rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    quals = [row["QUAL"] for row in rows if not math.isnan(row["QUAL"])]
    vafs = [row["VAF"] for row in rows if not math.isnan(row["VAF"])]
    allele_delta = [row["AlleleDelta"] for row in rows if not math.isnan(row["AlleleDelta"])]
    cramers_v = [row["CramersV"] for row in rows if not math.isnan(row["CramersV"])]
    quality_score = [row["Quality_Score"] for row in rows if not math.isnan(row["Quality_Score"])]
    pairwise = [row["PairwiseMedianDist"] for row in rows if not math.isnan(row["PairwiseMedianDist"])]
    dominant = Counter(str(row.get("DominantLabel", "")) for row in rows if row.get("DominantLabel"))
    classes = Counter(
        str(row.get("_verification_class_current", ""))
        for row in rows
        if row.get("_verification_class_current")
    )
    selection_fields = sorted(
        {str(row.get("_verification_selection_field")) for row in rows if row.get("_verification_selection_field")}
    )
    schema_statuses = sorted(
        {str(row.get("_verification_schema_status")) for row in rows if row.get("_verification_schema_status")}
    )
    return {
        "sample": sample,
        "platform": platform,
        "rule_id": rule_id,
        "subset": subset,
        "count": len(rows),
        "median_qual": f"{median(quals):.4f}" if quals else "",
        "median_vaf": f"{median(vafs):.4f}" if vafs else "",
        "median_allele_delta": f"{median(allele_delta):.4f}" if allele_delta else "",
        "median_cramers_v": f"{median(cramers_v):.4f}" if cramers_v else "",
        "median_quality_score": f"{median(quality_score):.4f}" if quality_score else "",
        "median_pairwise_median_dist": f"{median(pairwise):.4f}" if pairwise else "",
        "dominant_label_top": "; ".join(f"{key}:{value}" for key, value in dominant.most_common(3)),
        "verification_class_top": "; ".join(f"{key}:{value}" for key, value in classes.most_common(3)),
        "verification_selection_field": ";".join(selection_fields),
        "verification_schema_status": ";".join(schema_statuses),
        "unknown_current_class_count": classes.get("UnknownCurrentClass", 0),
    }


def load_joined_rows(
    run_dir: Path,
) -> tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    joined: Dict[str, Dict[str, object]] = {}
    schema_metadata: List[Dict[str, object]] = []
    for scope in ("TP", "FP"):
        summary_path = run_dir / f"intersubmod_{scope.lower()}" / "significance_summary.csv"
        vcf_path = resolve_longphase_vcf(run_dir, scope)
        summary_rows, metadata = parse_summary(summary_path, scope)
        schema_metadata.append({"run_dir": str(run_dir), **metadata})
        vcf_rows = parse_vcf_features(vcf_path, scope)
        for key, vcf_row in vcf_rows.items():
            if key not in summary_rows:
                continue
            joined[key] = {
                **summary_rows[key],
                **vcf_row,
                "key": key,
            }
    return list(joined.values()), schema_metadata


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rules = build_rules()
    focus_rules = list(dict.fromkeys(args.focus_rule))
    comparison_rows: List[Dict[str, object]] = []
    subset_rows: List[Dict[str, object]] = []
    best_rows: List[Dict[str, object]] = []
    schema_rows: List[Dict[str, object]] = []

    for run_dir_str in args.run_dir:
        run_dir = Path(run_dir_str).resolve()
        metrics = read_json(run_dir / "metrics.json")
        sample = str(metrics.get("sample") or run_dir.parent.name)
        platform = infer_platform(sample)
        truth_total = int(metrics.get("truth_total", 0))
        baseline = metrics.get("baseline", {})
        baseline_tp = int(baseline.get("tp", 0))
        baseline_fp = int(baseline.get("fp", 0))
        baseline_f1 = float(baseline.get("f1", 0.0))
        rows, run_schema_rows = load_joined_rows(run_dir)
        schema_rows.extend(run_schema_rows)

        best_rule_id = ""
        best_f1 = baseline_f1
        best_delta = 0.0

        for rule_id, rule_func in rules.items():
            triggered = [row for row in rows if rule_func(row)]
            tp_removed = sum(1 for row in triggered if row["label"] == "TP")
            fp_removed = sum(1 for row in triggered if row["label"] == "FP")
            current = compute_metrics(baseline_tp - tp_removed, baseline_fp - fp_removed, truth_total)
            delta = current["f1"] - baseline_f1
            if delta > best_delta:
                best_rule_id = rule_id
                best_f1 = current["f1"]
                best_delta = delta

            comparison_rows.append(
                {
                    "sample": sample,
                    "platform": platform,
                    "run_dir": str(run_dir),
                    "rule_id": rule_id,
                    "trigger_count": len(triggered),
                    "tp_removed": tp_removed,
                    "fp_removed": fp_removed,
                    "fp_tp_ratio": f"{(fp_removed / tp_removed):.4f}" if tp_removed else "inf",
                    "precision": f"{current['precision']:.6f}",
                    "recall": f"{current['recall']:.6f}",
                    "f1": f"{current['f1']:.6f}",
                    "delta_f1_vs_baseline": f"{delta:.6f}",
                    "baseline_f1": f"{baseline_f1:.6f}",
                    "truth_total": truth_total,
                }
            )

            if rule_id in focus_rules:
                tp_rows = [row for row in triggered if row["label"] == "TP"]
                fp_rows = [row for row in triggered if row["label"] == "FP"]
                subset_rows.append(summarize_subset(sample, platform, rule_id, "TP", tp_rows))
                subset_rows.append(summarize_subset(sample, platform, rule_id, "FP", fp_rows))

        best_rows.append(
            {
                "sample": sample,
                "platform": platform,
                "run_dir": str(run_dir),
                "baseline_f1": f"{baseline_f1:.6f}",
                "best_rule_id": best_rule_id or "baseline",
                "best_f1": f"{best_f1:.6f}",
                "best_delta_f1_vs_baseline": f"{best_delta:.6f}",
            }
        )

    write_tsv_rows(output_dir / "candidate_rule_comparison.tsv", COMPARISON_FIELDS, comparison_rows)
    write_tsv_rows(output_dir / "focus_rule_subset_summary.tsv", SUBSET_FIELDS, subset_rows)
    write_tsv_rows(output_dir / "best_rule_by_sample.tsv", BEST_FIELDS, best_rows)
    write_tsv_rows(output_dir / "verification_schema_provenance.tsv", SCHEMA_FIELDS, schema_rows)

    md_lines = [
        "# Candidate Rule Summary",
        "",
        "- 規則來源：舊 combo、目前 `AD+VAF` 核心、以及其簡化/加回 `CramersV` 版本",
        "- 用途：快速比較各樣本的 TP/FP 移除模式與 F1 delta",
        "",
        "## Best Rule By Sample",
        "",
        "| sample | platform | baseline_f1 | best_rule_id | best_f1 | best_delta_f1_vs_baseline |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in best_rows:
        md_lines.append(
            f"| {row['sample']} | {row['platform']} | {row['baseline_f1']} | {row['best_rule_id']} | {row['best_f1']} | {row['best_delta_f1_vs_baseline']} |"
        )
    md_lines.extend(
        [
            "",
            "- 詳細表格：`candidate_rule_comparison.tsv`",
            "- Focus subset：`focus_rule_subset_summary.tsv`",
        ]
    )
    (output_dir / "candidate_rule_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[analyze_candidate_rules] Wrote {output_dir / 'candidate_rule_comparison.tsv'}")
    print(f"[analyze_candidate_rules] Wrote {output_dir / 'focus_rule_subset_summary.tsv'}")
    print(f"[analyze_candidate_rules] Wrote {output_dir / 'best_rule_by_sample.tsv'}")


if __name__ == "__main__":
    main()
