#!/usr/bin/env python3
"""Validate exact alignment tags captured from a LongPhase-S BAM stream."""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import pysam


EXPECTED_HP = {".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"}


def parse_log(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "tumor_snp_count": r"Tumor SNP count:\s*(\d+)",
        "total_alignment": r"total alignment\s*:\s*(\d+)",
        "total_unmapped": r"total unmapped\s*:\s*(\d+)",
        "total_tagged": r"total tagged alignments\s*:\s*(\d+)",
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


def scan_sidecar(path):
    counts = Counter()
    hp = Counter()
    hp_ps = Counter()
    unknown = Counter()
    duplicate_rows = duplicate_conflicts = 0
    locus = None
    identities = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                counts["malformed"] += 1
                continue
            chrom, start, end, qname, flag, _mapq, cigar_digest, tag, ps = fields
            counts["rows"] += 1
            hp[tag] += 1
            if ps != ".":
                hp_ps[tag] += 1
            if tag not in EXPECTED_HP:
                unknown[tag] += 1
            current_locus = (chrom, start)
            if current_locus != locus:
                locus = current_locus
                identities = {}
            key = (qname, end, flag, cigar_digest)
            value = (tag, ps)
            if key in identities:
                duplicate_rows += 1
                duplicate_conflicts += int(identities[key] != value)
            else:
                identities[key] = value
    return counts, hp, hp_ps, unknown, duplicate_rows, duplicate_conflicts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sidecar", required=True, type=Path)
    parser.add_argument("--capture-summary", required=True, type=Path)
    parser.add_argument("--execution-log", required=True, type=Path)
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--recalibrated-vcf", required=True, type=Path)
    parser.add_argument("--region", help="Optional chromosome for a bounded smoke validation")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    capture = json.loads(args.capture_summary.read_text(encoding="utf-8"))
    execution = parse_log(args.execution_log)
    counts, hp, hp_ps, unknown, duplicate_rows, duplicate_conflicts = scan_sidecar(args.sidecar)
    input_all_keys, _ = vcf_keys(args.input_vcf)
    input_keys, input_filters = vcf_keys(args.input_vcf, args.region)
    output_keys, output_filters = vcf_keys(args.recalibrated_vcf, args.region)
    missing, extra = input_keys - output_keys, output_keys - input_keys
    tagged = sum(value for tag, value in hp.items() if tag != ".")
    checks = {
        "truth_flags_absent": execution["truth_vcf_empty"] and execution["truth_bed_empty"]
                              and execution["benchmark_removal_absent"],
        "parser_count_matches_input": execution["tumor_snp_count"] == sum(input_all_keys.values()),
        "capture_pass": capture.get("pass") is True,
        "execution_alignment_count_matches_capture": execution["total_alignment"] == capture["alignment_classes"]["total"],
        "sidecar_row_count_matches_capture": counts["rows"] == capture["rows_mapped"],
        "tagged_count_matches_execution": execution["total_tagged"] == tagged,
        "sidecar_no_malformed_rows": counts["malformed"] == 0,
        "sidecar_no_unknown_HP": not unknown,
        "sidecar_no_exact_identity_conflicts": duplicate_conflicts == 0,
        "recalibrated_preserves_all_input_keys": not missing and not extra,
    }
    output = {
        "schema_version": "1.0",
        "sidecar": str(args.sidecar),
        "region": args.region or "all",
        "capture": capture,
        "execution": execution,
        "HP_counts": dict(hp),
        "HP_with_PS_counts": dict(hp_ps),
        "unknown_HP_counts": dict(unknown),
        "duplicate_exact_alignment_rows": duplicate_rows,
        "duplicate_exact_alignment_conflicts": duplicate_conflicts,
        "input_filters": dict(input_filters),
        "recalibrated_filters": dict(output_filters),
        "record_key_missing": sum(missing.values()),
        "record_key_extra": sum(extra.values()),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"STREAM SIDECAR: rows={counts['rows']} HP={dict(hp)} exact_conflicts={duplicate_conflicts}")
    print(f"  checks={checks} -> {'PASS' if output['pass'] else 'FAIL'}")
    raise SystemExit(0 if output["pass"] else 1)


if __name__ == "__main__":
    main()
