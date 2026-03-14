#!/usr/bin/env python3
"""Refine label-first Strong calls into suspect/trusted subgroups and evaluate filter impact."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, List

from research_common import compute_metrics, infer_platform, load_tsv_rows, read_json, to_float, write_tsv_rows


PER_REGION_FIELDS = [
    "sample",
    "platform",
    "sample_dir",
    "source_scope",
    "region_key",
    "cluster_class",
    "label_class",
    "agreement_type",
    "class_shift",
    "refined_label",
    "suspect_flags",
    "reads",
    "vaf",
    "allele_delta",
    "cramers_v",
    "pairwise_median_dist",
    "quality_score",
]

RULE_FIELDS = [
    "sample",
    "platform",
    "sample_dir",
    "rule_id",
    "trigger_count",
    "tp_removed",
    "fp_removed",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_intersubmod",
    "baseline_f1",
    "truth_total",
    "notes",
]

BEST_FIELDS = [
    "sample",
    "platform",
    "sample_dir",
    "baseline_f1",
    "best_rule_id",
    "best_f1",
    "best_delta_f1_vs_intersubmod",
]

SUBSET_FIELDS = [
    "sample",
    "platform",
    "rule_id",
    "source_scope",
    "count",
    "median_vaf",
    "median_allele_delta",
    "median_cramers_v",
    "median_pairwise_median_dist",
    "top_class_shift",
    "top_agreement_type",
]


RuleFunc = Callable[[Dict[str, object]], bool]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", action="append", required=True, help="Round bundle sample dir")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--ad-threshold", type=float, default=0.15, help="AlleleDelta suspect threshold")
    parser.add_argument("--vaf-threshold", type=float, default=0.15, help="VAF suspect threshold")
    parser.add_argument("--cv-threshold", type=float, default=0.05, help="CramersV suspect threshold")
    return parser.parse_args()


def load_summary_rows(path: Path, scope: str) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            region_key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
            rows[region_key] = {
                "source_scope": scope,
                "AlleleDelta": to_float(row.get("AlleleDelta")),
                "CramersV": to_float(row.get("CramersV")),
                "PairwiseMedianDist": to_float(row.get("PairwiseMedianDist")),
                "Quality_Score": to_float(row.get("Quality_Score")),
            }
    return rows


def open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def resolve_longphase_vcf(run_dir: Path, scope: str) -> Path:
    direct = run_dir / "longphase_s" / f"filtered_snv_{scope}.vcf.gz"
    if direct.exists():
        return direct
    plain = run_dir / "longphase_s" / f"filtered_snv_{scope}.vcf"
    if plain.exists():
        return plain
    parent = run_dir.parent
    for sibling in sorted(parent.iterdir(), reverse=True):
        if not sibling.is_dir() or sibling == run_dir:
            continue
        direct = sibling / "longphase_s" / f"filtered_snv_{scope}.vcf.gz"
        if direct.exists():
            return direct
        plain = sibling / "longphase_s" / f"filtered_snv_{scope}.vcf"
        if plain.exists():
            return plain
    raise FileNotFoundError(f"filtered_snv_{scope}.vcf(.gz) not found for {run_dir}")


def load_vcf_vaf(path: Path) -> Dict[str, float]:
    rows: Dict[str, float] = {}
    with open_vcf(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            chrom, pos, _vid, ref, alt, _qual, _filt, info, *rest = raw.rstrip("\n").split("\t")
            info_map: Dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info_map[key] = value
            fmt_map: Dict[str, str] = {}
            if len(rest) >= 2:
                fmt_map = dict(zip(rest[0].split(":"), rest[1].split(":")))
            vaf = math.nan
            for key in ("VAF", "AF"):
                if key in fmt_map:
                    vaf = to_float(fmt_map[key].split(",")[0])
                    break
                if key in info_map:
                    vaf = to_float(info_map[key].split(",")[0])
                    break
            rows[f"{chrom}:{pos}:{ref}:{alt}"] = vaf
    return rows


def median(values: List[float]) -> str:
    clean = sorted(value for value in values if not math.isnan(value))
    if not clean:
        return ""
    n = len(clean)
    mid = n // 2
    if n % 2:
        return f"{clean[mid]:.6f}"
    return f"{((clean[mid - 1] + clean[mid]) / 2.0):.6f}"


def build_rules() -> Dict[str, RuleFunc]:
    return {
        "strong_noise_upgrade": lambda row: row["label_class"] == "Strong" and row["class_shift"] == "Noise->Strong",
        "strong_any_upgrade": lambda row: row["label_class"] == "Strong"
        and row["class_shift"] in {"Weak->Strong", "Noise->Strong"},
        "strong_low_vaf_high_ad": lambda row: row["label_class"] == "Strong"
        and row["low_vaf_high_ad"],
        "strong_upgrade_low_vaf_high_ad": lambda row: row["label_class"] == "Strong"
        and row["class_shift"] in {"Weak->Strong", "Noise->Strong"}
        and row["low_vaf_high_ad"],
        "strong_noise_or_low_vaf_high_ad": lambda row: row["label_class"] == "Strong"
        and (row["class_shift"] == "Noise->Strong" or row["low_vaf_high_ad"]),
        "strong_noise_or_low_vaf_high_ad_cv005": lambda row: row["label_class"] == "Strong"
        and (row["class_shift"] == "Noise->Strong" or row["low_vaf_high_ad_cv"]),
    }


def load_baseline_metrics(sample_dir: Path) -> Dict[str, float]:
    rows = load_tsv_rows(sample_dir / "benchmark_comparison.tsv")
    for row in rows:
        if row.get("method") == "InterSubMod":
            return {
                "truth_total": int(float(row.get("truth_total", 0))),
                "tp": int(float(row.get("tp", 0))),
                "fp": int(float(row.get("fp", 0))),
                "f1": float(row.get("f1", 0.0)),
            }
    raise FileNotFoundError(f"InterSubMod row not found in {sample_dir / 'benchmark_comparison.tsv'}")


def classify_refined_label(row: Dict[str, object]) -> str:
    if row["label_class"] != "Strong":
        return str(row["label_class"])
    if row["class_shift"] == "Noise->Strong":
        if row["low_vaf_high_ad_cv"]:
            return "Strong-NoiseUpgradeLowVAF"
        return "Strong-NoiseUpgrade"
    if row["class_shift"] == "Weak->Strong":
        if row["low_vaf_high_ad"]:
            return "Strong-UpgradeSuspectLowVAF"
        return "Strong-Upgraded"
    if row["low_vaf_high_ad"]:
        return "Strong-SuspectLowVAF"
    return "Strong-Trusted"


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rules = build_rules()
    per_region_rows: List[Dict[str, object]] = []
    rule_rows: List[Dict[str, object]] = []
    best_rows: List[Dict[str, object]] = []
    subset_rows: List[Dict[str, object]] = []

    for sample_dir_str in args.sample_dir:
        sample_dir = Path(sample_dir_str).resolve()
        context = read_json(sample_dir / "round_context.json")
        sample = str(context.get("sample") or sample_dir.name)
        platform = str(context.get("platform") or infer_platform(sample))
        source_run_dir = Path(str(context["source_run_dir"]))
        agreement_rows = load_tsv_rows(sample_dir / "label_cluster_agreement.tsv")
        baseline = load_baseline_metrics(sample_dir)
        tp_summary = load_summary_rows(source_run_dir / "intersubmod_tp" / "significance_summary.csv", "tp")
        fp_summary = load_summary_rows(source_run_dir / "intersubmod_fp" / "significance_summary.csv", "fp")
        summary_index = {**tp_summary, **fp_summary}
        vaf_index = {
            **load_vcf_vaf(resolve_longphase_vcf(source_run_dir, "tp")),
            **load_vcf_vaf(resolve_longphase_vcf(source_run_dir, "fp")),
        }

        enriched_rows: List[Dict[str, object]] = []
        for row in agreement_rows:
            region_key = str(row.get("region_key", ""))
            summary = summary_index.get(region_key, {})
            vaf = vaf_index.get(region_key, math.nan)
            allele_delta = to_float(summary.get("AlleleDelta"))
            cramers_v = to_float(summary.get("CramersV"))
            pairwise_median_dist = to_float(summary.get("PairwiseMedianDist"))
            quality_score = to_float(summary.get("Quality_Score"))

            low_vaf_high_ad = (
                not math.isnan(vaf)
                and vaf < args.vaf_threshold
                and not math.isnan(allele_delta)
                and allele_delta > args.ad_threshold
            )
            low_vaf_high_ad_cv = low_vaf_high_ad and not math.isnan(cramers_v) and cramers_v < args.cv_threshold

            enriched = {
                **row,
                "sample": sample,
                "platform": platform,
                "sample_dir": str(sample_dir),
                "vaf": vaf,
                "allele_delta": allele_delta,
                "cramers_v": cramers_v,
                "pairwise_median_dist": pairwise_median_dist,
                "quality_score": quality_score,
                "low_vaf_high_ad": low_vaf_high_ad,
                "low_vaf_high_ad_cv": low_vaf_high_ad_cv,
            }
            suspect_flags: List[str] = []
            if row.get("label_class") == "Strong":
                if row.get("class_shift") == "Noise->Strong":
                    suspect_flags.append("noise_upgrade")
                elif row.get("class_shift") == "Weak->Strong":
                    suspect_flags.append("weak_upgrade")
                if low_vaf_high_ad:
                    suspect_flags.append("low_vaf_high_ad")
                if low_vaf_high_ad_cv:
                    suspect_flags.append("low_vaf_high_ad_cv")

            enriched["refined_label"] = classify_refined_label(enriched)
            enriched["suspect_flags"] = ",".join(suspect_flags)
            enriched_rows.append(enriched)
            per_region_rows.append(
                {
                    "sample": sample,
                    "platform": platform,
                    "sample_dir": str(sample_dir),
                    "source_scope": row.get("source_scope", ""),
                    "region_key": region_key,
                    "cluster_class": row.get("cluster_class", ""),
                    "label_class": row.get("label_class", ""),
                    "agreement_type": row.get("agreement_type", ""),
                    "class_shift": row.get("class_shift", ""),
                    "refined_label": enriched["refined_label"],
                    "suspect_flags": enriched["suspect_flags"],
                    "reads": row.get("reads", ""),
                    "vaf": f"{vaf:.6f}" if not math.isnan(vaf) else "",
                    "allele_delta": f"{allele_delta:.6f}" if not math.isnan(allele_delta) else "",
                    "cramers_v": f"{cramers_v:.6f}" if not math.isnan(cramers_v) else "",
                    "pairwise_median_dist": f"{pairwise_median_dist:.6f}" if not math.isnan(pairwise_median_dist) else "",
                    "quality_score": f"{quality_score:.6f}" if not math.isnan(quality_score) else "",
                }
            )

        best_rule_id = "baseline"
        best_f1 = baseline["f1"]
        best_delta = 0.0
        for rule_id, rule in rules.items():
            triggered = [row for row in enriched_rows if rule(row)]
            tp_removed = sum(1 for row in triggered if row.get("source_scope") == "tp")
            fp_removed = sum(1 for row in triggered if row.get("source_scope") == "fp")
            metrics = compute_metrics(
                baseline["tp"] - tp_removed,
                baseline["fp"] - fp_removed,
                baseline["truth_total"],
            )
            delta = metrics["f1"] - baseline["f1"]
            if delta > best_delta:
                best_delta = delta
                best_rule_id = rule_id
                best_f1 = metrics["f1"]

            notes = []
            if rule_id.endswith("cv005"):
                notes.append("requires low CramersV")
            if "noise" in rule_id:
                notes.append("targets Noise->Strong risk")
            if "low_vaf_high_ad" in rule_id:
                notes.append("targets low VAF + high AlleleDelta")

            rule_rows.append(
                {
                    "sample": sample,
                    "platform": platform,
                    "sample_dir": str(sample_dir),
                    "rule_id": rule_id,
                    "trigger_count": len(triggered),
                    "tp_removed": tp_removed,
                    "fp_removed": fp_removed,
                    "precision": f"{metrics['precision']:.6f}",
                    "recall": f"{metrics['recall']:.6f}",
                    "f1": f"{metrics['f1']:.6f}",
                    "delta_f1_vs_intersubmod": f"{delta:.6f}",
                    "baseline_f1": f"{baseline['f1']:.6f}",
                    "truth_total": baseline["truth_total"],
                    "notes": "; ".join(notes),
                }
            )

            for scope in ("tp", "fp"):
                subset = [row for row in triggered if row.get("source_scope") == scope]
                class_shifts = Counter(str(row.get("class_shift", "")) for row in subset if row.get("class_shift"))
                agreement_types = Counter(str(row.get("agreement_type", "")) for row in subset if row.get("agreement_type"))
                subset_rows.append(
                    {
                        "sample": sample,
                        "platform": platform,
                        "rule_id": rule_id,
                        "source_scope": scope,
                        "count": len(subset),
                        "median_vaf": median([float(row["vaf"]) for row in subset if not math.isnan(float(row["vaf"]))]),
                        "median_allele_delta": median(
                            [float(row["allele_delta"]) for row in subset if not math.isnan(float(row["allele_delta"]))]
                        ),
                        "median_cramers_v": median(
                            [float(row["cramers_v"]) for row in subset if not math.isnan(float(row["cramers_v"]))]
                        ),
                        "median_pairwise_median_dist": median(
                            [float(row["pairwise_median_dist"]) for row in subset if not math.isnan(float(row["pairwise_median_dist"]))]
                        ),
                        "top_class_shift": "; ".join(f"{key}:{value}" for key, value in class_shifts.most_common(3)),
                        "top_agreement_type": "; ".join(f"{key}:{value}" for key, value in agreement_types.most_common(3)),
                    }
                )

        best_rows.append(
            {
                "sample": sample,
                "platform": platform,
                "sample_dir": str(sample_dir),
                "baseline_f1": f"{baseline['f1']:.6f}",
                "best_rule_id": best_rule_id,
                "best_f1": f"{best_f1:.6f}",
                "best_delta_f1_vs_intersubmod": f"{best_delta:.6f}",
            }
        )

    write_tsv_rows(output_dir / "refined_label_per_region.tsv", PER_REGION_FIELDS, per_region_rows)
    write_tsv_rows(output_dir / "refined_label_rule_comparison.tsv", RULE_FIELDS, rule_rows)
    write_tsv_rows(output_dir / "refined_label_best_by_sample.tsv", BEST_FIELDS, best_rows)
    write_tsv_rows(output_dir / "refined_label_subset_summary.tsv", SUBSET_FIELDS, subset_rows)

    md_lines = [
        "# Refined Strong Label Summary",
        "",
        "- 目的：將 `label-first Strong` 細分為較可信與較可疑子群，避免過早把 5kHz 特化規則寫入全域流程。",
        "- 主要檔案：`refined_label_rule_comparison.tsv`、`refined_label_per_region.tsv`、`refined_label_best_by_sample.tsv`",
        "",
        "## 每個樣本最佳 refined rule",
        "",
        "| sample | platform | baseline_f1 | best_rule_id | best_f1 | best_delta_f1_vs_intersubmod |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in best_rows:
        md_lines.append(
            f"| {row['sample']} | {row['platform']} | {row['baseline_f1']} | {row['best_rule_id']} | {row['best_f1']} | {row['best_delta_f1_vs_intersubmod']} |"
        )
    (output_dir / "refined_label_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[refine_strong_labels] Wrote {output_dir / 'refined_label_per_region.tsv'}")
    print(f"[refine_strong_labels] Wrote {output_dir / 'refined_label_rule_comparison.tsv'}")
    print(f"[refine_strong_labels] Wrote {output_dir / 'refined_label_best_by_sample.tsv'}")


if __name__ == "__main__":
    main()
