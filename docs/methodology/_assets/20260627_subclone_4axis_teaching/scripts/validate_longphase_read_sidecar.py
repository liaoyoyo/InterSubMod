#!/usr/bin/env python3
"""Validate a LongPhase-S read log and recalibrated VCF without rescanning BAM."""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import pysam


EXPECTED_TAGS = {".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"}
CHROM_RANK = {f"chr{value}": value for value in range(1, 23)} | {"chrX": 23, "chrY": 24}


def normalize_tag(value):
    return value[1:] if value.startswith("H") else value


def parse_execution_log(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "tumor_snp_count": r"Tumor SNP count:\s*(\d+)",
        "total_alignment": r"total alignment\s*:\s*(\d+)",
        "lower_mapping_quality": r"lower mapping quality\s*:\s*(\d+)",
        "start_after_last_variant": r"start pos > last variant pos\s*:\s*(\d+)",
    }
    values = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        values[name] = int(match.group(1)) if match else None
    values["truth_vcf_empty"] = bool(re.search(r"truth VCF file\s*:\s*$", text, re.MULTILINE))
    values["truth_bed_empty"] = bool(re.search(r"truth BED file\s*:\s*$", text, re.MULTILINE))
    values["benchmark_removal_absent"] = "removing tumor & truth somatic variants outside bed regions" not in text
    return values


def vcf_keys(path, region=None):
    keys = Counter()
    filters = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            if region and record.contig != region:
                continue
            keys[(record.contig, int(record.pos), record.ref, tuple(record.alts or ()))] += 1
            for label in list(record.filter.keys()) or ["."]:
                filters[label] += 1
    return keys, filters


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sidecar", required=True, type=Path)
    parser.add_argument("--execution-log", required=True, type=Path)
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--recalibrated-vcf", required=True, type=Path)
    parser.add_argument("--region", help="Optional chromosome for a bounded smoke validation")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    counts = Counter()
    tag_counts = Counter()
    tag_with_ps = Counter()
    unknown = Counter()
    sorted_ok = True
    duplicate_rows = duplicate_conflicts = 0
    previous_order = None
    locus = None
    locus_tags = {}
    with args.sidecar.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                counts["malformed_rows"] += 1
                continue
            read_id, chrom = fields[0], fields[1]
            try:
                start = int(fields[2])
            except ValueError:
                counts["malformed_rows"] += 1
                continue
            tag = normalize_tag(fields[5])
            ps = fields[6]
            counts["rows"] += 1
            tag_counts[tag] += 1
            if ps != ".":
                tag_with_ps[tag] += 1
            if tag not in EXPECTED_TAGS:
                unknown[tag] += 1
            order = (CHROM_RANK.get(chrom, 1000), start)
            if previous_order is not None and order < previous_order:
                sorted_ok = False
            previous_order = order
            current_locus = (chrom, start)
            if current_locus != locus:
                locus = current_locus
                locus_tags = {}
            key = read_id
            value = (tag, ps)
            if key in locus_tags:
                duplicate_rows += 1
                duplicate_conflicts += int(locus_tags[key] != value)
            else:
                locus_tags[key] = value
    execution = parse_execution_log(args.execution_log)
    expected_rows = None
    if all(execution.get(name) is not None for name in
           ("total_alignment", "lower_mapping_quality", "start_after_last_variant")):
        expected_rows = (execution["total_alignment"] - execution["lower_mapping_quality"]
                         - execution["start_after_last_variant"])
    input_keys, input_filters = vcf_keys(args.input_vcf, args.region)
    output_keys, output_filters = vcf_keys(args.recalibrated_vcf, args.region)
    missing = input_keys - output_keys
    extra = output_keys - input_keys
    checks = {
        "truth_flags_absent": execution["truth_vcf_empty"] and execution["truth_bed_empty"]
                              and execution["benchmark_removal_absent"],
        "parser_count_matches_input": execution["tumor_snp_count"] == sum(input_keys.values()),
        "sidecar_row_count_exact": expected_rows == counts["rows"],
        "sidecar_coordinate_sorted": sorted_ok,
        "sidecar_no_malformed_rows": counts["malformed_rows"] == 0,
        "sidecar_no_unknown_HP": not unknown,
        "sidecar_no_duplicate_key_conflicts": duplicate_conflicts == 0,
        "recalibrated_preserves_all_input_keys": not missing and not extra,
    }
    output = {
        "schema_version": "1.0",
        "sidecar": str(args.sidecar),
        "region": args.region or "all",
        "execution": execution,
        "expected_rows": expected_rows,
        "counts": dict(counts),
        "HP_counts": dict(tag_counts),
        "HP_with_PS_counts": dict(tag_with_ps),
        "unknown_HP_counts": dict(unknown),
        "duplicate_alignment_key_rows": duplicate_rows,
        "duplicate_alignment_key_conflicts": duplicate_conflicts,
        "input_filters": dict(input_filters),
        "recalibrated_filters": dict(output_filters),
        "record_key_missing": sum(missing.values()),
        "record_key_extra": sum(extra.values()),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"SIDECAR VALIDATION: rows={counts['rows']} HP={dict(tag_counts)} duplicate_conflicts={duplicate_conflicts}")
    print(f"  checks={checks} -> {'PASS' if output['pass'] else 'FAIL'}")
    raise SystemExit(0 if output["pass"] else 1)


if __name__ == "__main__":
    main()
