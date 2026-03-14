#!/usr/bin/env python3
"""Analyze caller TP lost by LongPhase and evaluate simple rescue rules."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence

from research_common import compute_metrics, infer_caller_name, infer_platform, to_float, to_int, write_tsv_rows


RuleFunc = Callable[[Dict[str, object]], bool]


VARIANT_FIELDS = [
    "sample",
    "platform",
    "caller",
    "variant_group",
    "region_key",
    "qual",
    "af",
    "dp",
    "gq",
    "alt_count",
    "filter",
    "has_h_flag",
]

SUMMARY_FIELDS = [
    "sample",
    "platform",
    "caller",
    "variant_group",
    "count",
    "median_qual",
    "median_af",
    "median_dp",
    "median_gq",
    "median_alt_count",
    "h_flag_fraction",
]

RULE_FIELDS = [
    "sample",
    "platform",
    "caller",
    "rule_id",
    "trigger_count",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_longphase",
    "baseline_f1",
    "truth_total",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--caller-vcf", required=True, help="Caller VCF or benchmark TP/FP source root")
    parser.add_argument("--caller-tp-vcf", required=True, help="Caller benchmark TP VCF")
    parser.add_argument("--caller-fp-vcf", required=True, help="Caller benchmark FP VCF")
    parser.add_argument("--longphase-tp-vcf", required=True, help="LongPhase kept TP VCF")
    parser.add_argument("--longphase-fp-vcf", required=True, help="LongPhase kept FP VCF")
    parser.add_argument("--truth-total", type=int, required=True, help="Truth total for F1 computation")
    parser.add_argument("--baseline-tp", type=int, required=True, help="LongPhase baseline TP count")
    parser.add_argument("--baseline-fp", type=int, required=True, help="LongPhase baseline FP count")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def parse_vcf(path: Path) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with open_vcf(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            chrom, pos, _vid, ref, alt, qual, filt, info, *rest = raw.rstrip("\n").split("\t")
            info_set = set()
            info_map: Dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info_map[key] = value
                elif item:
                    info_set.add(item)
            fmt_map: Dict[str, str] = {}
            if len(rest) >= 2:
                fmt_map = dict(zip(rest[0].split(":"), rest[1].split(":")))

            af = math.nan
            for key in ("AF", "VAF"):
                if key in fmt_map:
                    af = to_float(fmt_map[key].split(",")[0])
                    break
                if key in info_map:
                    af = to_float(info_map[key].split(",")[0])
                    break

            dp = to_int(fmt_map.get("DP"))
            gq = to_int(fmt_map.get("GQ"))
            alt_count = 0
            ad = fmt_map.get("AD", "")
            if ad:
                parts = ad.split(",")
                if len(parts) >= 2:
                    alt_count = to_int(parts[1])

            rows[f"{chrom}:{pos}:{ref}:{alt}"] = {
                "region_key": f"{chrom}:{pos}:{ref}:{alt}",
                "qual": to_float(qual),
                "af": af,
                "dp": dp,
                "gq": gq,
                "alt_count": alt_count,
                "filter": filt,
                "has_h_flag": "H" in info_set,
            }
    return rows


def filter_vcf_to_keys(src: Path, keys: set[str], dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open_vcf(src) as handle, dst.open("w", encoding="utf-8") as out:
        for raw in handle:
            if raw.startswith("#"):
                out.write(raw)
                continue
            chrom, pos, _vid, ref, alt, *_rest = raw.rstrip("\n").split("\t")
            key = f"{chrom}:{pos}:{ref}:{alt}"
            if key in keys:
                out.write(raw)


def median(values: Sequence[float]) -> str:
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
        "af_ge_0_10": lambda row: row["af"] >= 0.10,
        "af_ge_0_15": lambda row: row["af"] >= 0.15,
        "af_ge_0_20": lambda row: row["af"] >= 0.20,
        "gq_ge_20": lambda row: row["gq"] >= 20,
        "dp_ge_30": lambda row: row["dp"] >= 30,
        "alt_ge_10": lambda row: row["alt_count"] >= 10,
        "h_flag": lambda row: row["has_h_flag"],
        "af_ge_0_15_and_alt_ge_10": lambda row: row["af"] >= 0.15 and row["alt_count"] >= 10,
        "af_ge_0_20_and_gq_ge_20": lambda row: row["af"] >= 0.20 and row["gq"] >= 20,
        "af_ge_0_20_and_h_flag": lambda row: row["af"] >= 0.20 and row["has_h_flag"],
    }


def group_summary(
    sample: str,
    platform: str,
    caller: str,
    variant_group: str,
    rows: Iterable[Dict[str, object]],
) -> Dict[str, object]:
    rows = list(rows)
    quals = [float(row["qual"]) for row in rows if not math.isnan(float(row["qual"]))]
    afs = [float(row["af"]) for row in rows if not math.isnan(float(row["af"]))]
    dps = [float(row["dp"]) for row in rows]
    gqs = [float(row["gq"]) for row in rows]
    alt_counts = [float(row["alt_count"]) for row in rows]
    h_fraction = sum(1 for row in rows if row["has_h_flag"]) / len(rows) if rows else 0.0
    return {
        "sample": sample,
        "platform": platform,
        "caller": caller,
        "variant_group": variant_group,
        "count": len(rows),
        "median_qual": median(quals),
        "median_af": median(afs),
        "median_dp": median(dps),
        "median_gq": median(gqs),
        "median_alt_count": median(alt_counts),
        "h_flag_fraction": f"{h_fraction:.6f}",
    }


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sample = args.sample
    platform = infer_platform(sample)
    caller = infer_caller_name(args.caller_vcf)

    caller_tp = parse_vcf(Path(args.caller_tp_vcf))
    caller_fp = parse_vcf(Path(args.caller_fp_vcf))
    longphase_tp = parse_vcf(Path(args.longphase_tp_vcf))
    longphase_fp = parse_vcf(Path(args.longphase_fp_vcf))

    caller_kept_tp_keys = set(caller_tp) & set(longphase_tp)
    caller_lost_tp_keys = set(caller_tp) - set(longphase_tp)
    caller_kept_fp_keys = set(caller_fp) & set(longphase_fp)
    caller_removed_fp_keys = set(caller_fp) - set(longphase_fp)

    group_index = {
        "caller_kept_tp": [caller_tp[key] for key in sorted(caller_kept_tp_keys)],
        "caller_lost_tp": [caller_tp[key] for key in sorted(caller_lost_tp_keys)],
        "caller_kept_fp": [caller_fp[key] for key in sorted(caller_kept_fp_keys)],
        "caller_removed_fp": [caller_fp[key] for key in sorted(caller_removed_fp_keys)],
    }

    variant_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    for group_name, rows in group_index.items():
        summary_rows.append(group_summary(sample, platform, caller, group_name, rows))
        for row in rows:
            variant_rows.append(
                {
                    "sample": sample,
                    "platform": platform,
                    "caller": caller,
                    "variant_group": group_name,
                    **row,
                }
            )

    baseline_metrics = compute_metrics(args.baseline_tp, args.baseline_fp, args.truth_total)
    rule_rows: List[Dict[str, object]] = []
    rules = build_rules()
    best_rule_id = "baseline"
    best_f1 = baseline_metrics["f1"]
    best_delta = 0.0

    for rule_id, rule in rules.items():
        lost_tp = group_index["caller_lost_tp"]
        removed_fp = group_index["caller_removed_fp"]
        tp_rescued = sum(1 for row in lost_tp if rule(row))
        fp_reintroduced = sum(1 for row in removed_fp if rule(row))
        metrics = compute_metrics(
            args.baseline_tp + tp_rescued,
            args.baseline_fp + fp_reintroduced,
            args.truth_total,
        )
        delta = metrics["f1"] - baseline_metrics["f1"]
        if delta > best_delta:
            best_delta = delta
            best_rule_id = rule_id
            best_f1 = metrics["f1"]
        notes = []
        if "af" in rule_id:
            notes.append("caller AF rescue")
        if "gq" in rule_id:
            notes.append("caller GQ rescue")
        if "alt" in rule_id:
            notes.append("caller alt-depth rescue")
        if "h_flag" in rule_id:
            notes.append("haplotype-unique rescue")
        rule_rows.append(
            {
                "sample": sample,
                "platform": platform,
                "caller": caller,
                "rule_id": rule_id,
                "trigger_count": tp_rescued + fp_reintroduced,
                "tp_rescued": tp_rescued,
                "fp_reintroduced": fp_reintroduced,
                "precision": f"{metrics['precision']:.6f}",
                "recall": f"{metrics['recall']:.6f}",
                "f1": f"{metrics['f1']:.6f}",
                "delta_f1_vs_longphase": f"{delta:.6f}",
                "baseline_f1": f"{baseline_metrics['f1']:.6f}",
                "truth_total": args.truth_total,
                "notes": "; ".join(notes),
            }
        )

    write_tsv_rows(output_dir / "longphase_rescue_variants.tsv", VARIANT_FIELDS, variant_rows)
    write_tsv_rows(output_dir / "longphase_rescue_group_summary.tsv", SUMMARY_FIELDS, summary_rows)
    write_tsv_rows(output_dir / "longphase_rescue_rule_comparison.tsv", RULE_FIELDS, rule_rows)
    filter_vcf_to_keys(Path(args.caller_vcf), caller_lost_tp_keys, output_dir / "caller_lost_tp.vcf")
    filter_vcf_to_keys(Path(args.caller_vcf), caller_removed_fp_keys, output_dir / "caller_removed_fp.vcf")

    md_lines = [
        "# LongPhase Rescue Summary",
        "",
        f"- sample: `{sample}`",
        f"- platform: `{platform}`",
        f"- caller: `{caller}`",
        "",
        "## 群組摘要",
        "",
        "| variant_group | count | median_af | median_dp | median_gq | median_alt_count | h_flag_fraction |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in summary_rows:
        md_lines.append(
            f"| {row['variant_group']} | {row['count']} | {row['median_af']} | {row['median_dp']} | {row['median_gq']} | {row['median_alt_count']} | {row['h_flag_fraction']} |"
        )
    md_lines.extend(
        [
            "",
            "## 最佳 rescue rule",
            "",
            f"- baseline LongPhase-S F1: `{baseline_metrics['f1']:.6f}`",
            f"- best rule: `{best_rule_id}`",
            f"- best F1: `{best_f1:.6f}`",
            f"- delta: `{best_delta:.6f}`",
            "",
            "詳細結果請見 `longphase_rescue_rule_comparison.tsv`。",
        ]
    )
    (output_dir / "longphase_rescue_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[analyze_longphase_rescue] Wrote {output_dir / 'longphase_rescue_group_summary.tsv'}")
    print(f"[analyze_longphase_rescue] Wrote {output_dir / 'longphase_rescue_rule_comparison.tsv'}")


if __name__ == "__main__":
    main()
