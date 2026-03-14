#!/usr/bin/env python3
"""Extract caller borderline rescue candidates from final VCF and downstream benchmark splits."""

from __future__ import annotations

import argparse
import gzip
import math
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Set, Tuple

from research_common import (
    infer_caller_name,
    infer_platform,
    to_float,
    to_int,
    write_tsv_rows,
)


CANDIDATE_FIELDS = [
    "sample",
    "platform",
    "caller",
    "mode",
    "region_key",
    "source_stage",
    "truth_status",
    "downstream_status",
    "filter",
    "qual",
    "gq",
    "dp",
    "af",
    "ad_ref",
    "ad_alt",
    "has_h_flag",
    "verdict_somatic",
    "verdict_subclonal",
    "verdict_germline",
    "pon_hit",
    "pon_hit_count",
    "multihap_flag",
    "noancestry_flag",
    "rescue_filter_eligible",
    "meets_numeric_thresholds",
    "candidate_eligible",
    "rescue_exclusion_reason",
]

SUMMARY_FIELDS = [
    "sample",
    "platform",
    "caller",
    "mode",
    "truth_status",
    "downstream_status",
    "count",
    "candidate_eligible_count",
    "rescue_filter_eligible_count",
    "median_qual",
    "median_gq",
    "median_dp",
    "median_af",
    "median_ad_alt",
    "h_flag_fraction",
    "verdict_somatic_fraction",
    "verdict_subclonal_fraction",
    "pon_hit_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--mode", required=True, help="paired | to")
    parser.add_argument("--caller-vcf", required=True, help="Caller final VCF (ClairS / ClairS-TO)")
    parser.add_argument("--truth-vcf", required=True, help="Truth VCF")
    parser.add_argument("--truth-bed", default="", help="Optional truth BED for scope restriction")
    parser.add_argument("--baseline-tp-vcf", required=True, help="Downstream kept TP VCF")
    parser.add_argument("--baseline-fp-vcf", required=True, help="Downstream kept FP VCF")
    parser.add_argument("--caller-name", default="", help="Optional caller label override")
    parser.add_argument("--dp-min", type=int, default=20, help="Minimum DP for candidate eligibility")
    parser.add_argument("--af-min", type=float, default=0.05, help="Minimum AF for candidate eligibility")
    parser.add_argument("--ad-alt-min", type=int, default=5, help="Minimum ALT count for candidate eligibility")
    parser.add_argument(
        "--include-filter",
        action="append",
        default=None,
        help="Rescue-eligible FILTER term (repeatable). Default: PASS,LowQual,MultiHap,NoAncestry",
    )
    parser.add_argument(
        "--exclude-filter",
        action="append",
        default=None,
        help="Excluded FILTER term (repeatable). Default: NonSomatic,ReadStartEnd,Realignment,LowAltMQ,LowAltBQ,VariantCluster,StrandBias,LowSeqEntropy,RefCall",
    )
    parser.add_argument("--bcftools", default=shutil.which("bcftools") or "bcftools", help="bcftools binary")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def run_cmd(cmd: List[str]) -> None:
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def vcf_key(fields: List[str]) -> str:
    return f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"


def parse_flag_info(info: str) -> Tuple[Set[str], Dict[str, str]]:
    flag_set: Set[str] = set()
    key_value: Dict[str, str] = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            key_value[key] = value
        elif item:
            flag_set.add(item)
    return flag_set, key_value


def parse_vcf_records(path: Path) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with open_vcf(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, qual, filt, info = fields[:8]
            fmt_map: Dict[str, str] = {}
            if len(fields) >= 10:
                fmt_map = dict(zip(fields[8].split(":"), fields[9].split(":")))
            info_flags, info_map = parse_flag_info(info)
            key = f"{chrom}:{pos}:{ref}:{alt}"

            filter_terms = []
            if filt and filt != ".":
                filter_terms = filt.split(";")
            filter_terms = ["PASS"] if not filter_terms else filter_terms

            ad_ref = 0
            ad_alt = 0
            ad = fmt_map.get("AD", "")
            if ad:
                parts = ad.split(",")
                if len(parts) >= 2:
                    ad_ref = to_int(parts[0])
                    ad_alt = to_int(parts[1])

            af = math.nan
            for field_name in ("AF", "VAF"):
                if field_name in fmt_map:
                    af = to_float(fmt_map[field_name].split(",")[0])
                    break
                if field_name in info_map:
                    af = to_float(info_map[field_name].split(",")[0])
                    break

            pon_hit_count = sum(1 for idx in range(1, 5) if f"PoN_{idx}" in info_flags)
            rows[key] = {
                "region_key": key,
                "filter_terms": filter_terms,
                "filter": ";".join(filter_terms),
                "qual": to_float(qual, 0.0),
                "gq": to_int(fmt_map.get("GQ")),
                "dp": to_int(fmt_map.get("DP")),
                "af": af,
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
                "has_h_flag": "H" in info_flags,
                "verdict_somatic": "Verdict_Somatic" in info_flags,
                "verdict_subclonal": "Verdict_SubclonalSomatic" in info_flags,
                "verdict_germline": "Verdict_Germline" in info_flags,
                "pon_hit": pon_hit_count > 0,
                "pon_hit_count": pon_hit_count,
                "multihap_flag": "MultiHap" in filter_terms,
                "noancestry_flag": "NoAncestry" in filter_terms,
            }
    return rows


def parse_vcf_keys(path: Path) -> Set[str]:
    keys: Set[str] = set()
    with open_vcf(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            keys.add(vcf_key(raw.rstrip("\n").split("\t")))
    return keys


def compute_summary_value(values: Iterable[float]) -> str:
    clean = [value for value in values if not math.isnan(value)]
    if not clean:
        return ""
    return f"{median(clean):.6f}"


def build_group_summary(
    sample: str,
    platform: str,
    caller: str,
    mode: str,
    truth_status: str,
    downstream_status: str,
    rows: List[Dict[str, object]],
) -> Dict[str, object]:
    if not rows:
        return {}
    total = len(rows)
    return {
        "sample": sample,
        "platform": platform,
        "caller": caller,
        "mode": mode,
        "truth_status": truth_status,
        "downstream_status": downstream_status,
        "count": total,
        "candidate_eligible_count": sum(1 for row in rows if row["candidate_eligible"]),
        "rescue_filter_eligible_count": sum(1 for row in rows if row["rescue_filter_eligible"]),
        "median_qual": compute_summary_value(float(row["qual"]) for row in rows),
        "median_gq": compute_summary_value(float(row["gq"]) for row in rows),
        "median_dp": compute_summary_value(float(row["dp"]) for row in rows),
        "median_af": compute_summary_value(float(row["af"]) for row in rows),
        "median_ad_alt": compute_summary_value(float(row["ad_alt"]) for row in rows),
        "h_flag_fraction": f"{sum(1 for row in rows if row['has_h_flag']) / total:.6f}",
        "verdict_somatic_fraction": f"{sum(1 for row in rows if row['verdict_somatic']) / total:.6f}",
        "verdict_subclonal_fraction": f"{sum(1 for row in rows if row['verdict_subclonal']) / total:.6f}",
        "pon_hit_fraction": f"{sum(1 for row in rows if row['pon_hit']) / total:.6f}",
    }


def main() -> None:
    args = parse_args()
    sample = args.sample
    platform = infer_platform(sample)
    caller = args.caller_name or infer_caller_name(args.caller_vcf)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    include_filters = set(args.include_filter or ["PASS", "LowQual", "MultiHap", "NoAncestry"])
    exclude_filters = set(
        args.exclude_filter
        or [
            "NonSomatic",
            "ReadStartEnd",
            "Realignment",
            "LowAltMQ",
            "LowAltBQ",
            "VariantCluster",
            "StrandBias",
            "LowSeqEntropy",
            "RefCall",
        ]
    )

    with tempfile.TemporaryDirectory(prefix="borderline_rescue_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        caller_scoped = tmpdir_path / "caller.scoped.vcf.gz"
        truth_scoped = tmpdir_path / "truth.scoped.vcf.gz"
        caller_view_cmd = [args.bcftools, "view", "-v", "snps"]
        truth_view_cmd = [args.bcftools, "view", "-v", "snps"]
        if args.truth_bed:
            caller_view_cmd.extend(["-R", args.truth_bed])
            truth_view_cmd.extend(["-R", args.truth_bed])
        caller_view_cmd.extend([args.caller_vcf, "-Oz", "-o", str(caller_scoped)])
        truth_view_cmd.extend([args.truth_vcf, "-Oz", "-o", str(truth_scoped)])
        run_cmd(caller_view_cmd)
        run_cmd([args.bcftools, "index", "-f", str(caller_scoped)])
        run_cmd(truth_view_cmd)
        run_cmd([args.bcftools, "index", "-f", str(truth_scoped)])

        isec_dir = tmpdir_path / "isec"
        run_cmd([args.bcftools, "isec", "-p", str(isec_dir), str(caller_scoped), str(truth_scoped)])

        caller_rows = parse_vcf_records(caller_scoped)
        caller_tp_keys = parse_vcf_keys(isec_dir / "0002.vcf")
        caller_fp_keys = parse_vcf_keys(isec_dir / "0000.vcf")

    baseline_tp_keys = parse_vcf_keys(Path(args.baseline_tp_vcf))
    baseline_fp_keys = parse_vcf_keys(Path(args.baseline_fp_vcf))

    candidate_rows: List[Dict[str, object]] = []
    grouped_rows: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)

    for region_key in sorted(caller_rows):
        row = caller_rows[region_key]
        if region_key in caller_tp_keys:
            truth_status = "TP"
            downstream_status = "caller_kept_tp" if region_key in baseline_tp_keys else "caller_lost_tp"
        elif region_key in caller_fp_keys:
            truth_status = "FP"
            downstream_status = "caller_kept_fp" if region_key in baseline_fp_keys else "caller_removed_fp"
        else:
            continue

        filter_terms = set(row["filter_terms"])
        has_excluded_filter = any(term in exclude_filters for term in filter_terms)
        rescue_filter_eligible = not has_excluded_filter and all(term in include_filters for term in filter_terms)

        reasons: List[str] = []
        if has_excluded_filter:
            reasons.extend(sorted(term for term in filter_terms if term in exclude_filters))
        if row["dp"] < args.dp_min:
            reasons.append("low_dp")
        if row["ad_alt"] < args.ad_alt_min:
            reasons.append("low_alt_count")
        if math.isnan(float(row["af"])) or float(row["af"]) < args.af_min:
            reasons.append("low_af")

        meets_numeric = row["dp"] >= args.dp_min and row["ad_alt"] >= args.ad_alt_min and not math.isnan(
            float(row["af"])
        ) and float(row["af"]) >= args.af_min

        candidate = {
            "sample": sample,
            "platform": platform,
            "caller": caller,
            "mode": args.mode,
            "region_key": region_key,
            "source_stage": "caller_final_vcf",
            "truth_status": truth_status,
            "downstream_status": downstream_status,
            "filter": row["filter"],
            "qual": f"{float(row['qual']):.6f}",
            "gq": int(row["gq"]),
            "dp": int(row["dp"]),
            "af": f"{float(row['af']):.6f}" if not math.isnan(float(row["af"])) else "",
            "ad_ref": int(row["ad_ref"]),
            "ad_alt": int(row["ad_alt"]),
            "has_h_flag": row["has_h_flag"],
            "verdict_somatic": row["verdict_somatic"],
            "verdict_subclonal": row["verdict_subclonal"],
            "verdict_germline": row["verdict_germline"],
            "pon_hit": row["pon_hit"],
            "pon_hit_count": int(row["pon_hit_count"]),
            "multihap_flag": row["multihap_flag"],
            "noancestry_flag": row["noancestry_flag"],
            "rescue_filter_eligible": rescue_filter_eligible,
            "meets_numeric_thresholds": meets_numeric,
            "candidate_eligible": rescue_filter_eligible and meets_numeric,
            "rescue_exclusion_reason": ";".join(reasons),
        }
        candidate_rows.append(candidate)
        grouped_rows[(truth_status, downstream_status)].append(candidate)

    summary_rows = []
    for (truth_status, downstream_status), rows in sorted(grouped_rows.items()):
        summary = build_group_summary(sample, platform, caller, args.mode, truth_status, downstream_status, rows)
        if summary:
            summary_rows.append(summary)

    downstream_pool_rows = [
        row for row in candidate_rows if row["downstream_status"] in {"caller_lost_tp", "caller_removed_fp"}
    ]

    write_tsv_rows(output_dir / "borderline_candidate_pool.tsv", CANDIDATE_FIELDS, candidate_rows)
    write_tsv_rows(output_dir / "candidate_group_summary.tsv", SUMMARY_FIELDS, summary_rows)
    write_tsv_rows(output_dir / "downstream_lost_candidate_pool.tsv", CANDIDATE_FIELDS, downstream_pool_rows)

    md_lines = [
        f"# Borderline Rescue Candidate Summary: {sample} ({args.mode})",
        "",
        f"- Caller: `{caller}`",
        f"- Caller VCF: `{Path(args.caller_vcf).resolve()}`",
        f"- Baseline TP VCF: `{Path(args.baseline_tp_vcf).resolve()}`",
        f"- Baseline FP VCF: `{Path(args.baseline_fp_vcf).resolve()}`",
        f"- Candidate thresholds: `DP>={args.dp_min}`, `AD_alt>={args.ad_alt_min}`, `AF>={args.af_min}`",
        f"- Rescue-eligible FILTER: `{','.join(sorted(include_filters))}`",
        f"- Excluded FILTER: `{','.join(sorted(exclude_filters))}`",
        "",
        "## Downstream pool summary",
        "",
        "| truth_status | downstream_status | count | candidate_eligible_count | median_gq | median_af | median_ad_alt |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in summary_rows:
        if row["downstream_status"] not in {"caller_lost_tp", "caller_removed_fp"}:
            continue
        md_lines.append(
            f"| {row['truth_status']} | {row['downstream_status']} | {row['count']} | {row['candidate_eligible_count']} | {row['median_gq']} | {row['median_af']} | {row['median_ad_alt']} |"
        )
    (output_dir / "candidate_group_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[extract_borderline_rescue_candidates] Wrote {output_dir / 'borderline_candidate_pool.tsv'}")
    print(f"[extract_borderline_rescue_candidates] Wrote {output_dir / 'candidate_group_summary.tsv'}")


if __name__ == "__main__":
    main()
