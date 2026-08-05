#!/usr/bin/env python3
"""Audit record conservation and FILTER transitions across LongPhase-S."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pysam


def json_value(value: Any) -> Any:
    if isinstance(value, tuple):
        return [json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    return value


def filter_label(record: pysam.VariantRecord) -> str:
    values = list(record.filter.keys())
    return ";".join(values) if values else "."


def is_pass(record: pysam.VariantRecord) -> bool:
    return "PASS" in set(record.filter.keys())


def record_key(record: pysam.VariantRecord) -> tuple[str, int, str, tuple[str, ...]]:
    return record.contig, int(record.pos), record.ref, tuple(record.alts or ())


def payload_without_filter(record: pysam.VariantRecord) -> dict[str, Any]:
    return {
        "id": record.id,
        "qual": record.qual,
        "info": {key: json_value(value) for key, value in record.info.items()},
        "format_keys": list(record.format.keys()),
        "samples": {
            sample: {key: json_value(value) for key, value in values.items()}
            for sample, values in record.samples.items()
        },
    }


def scan(path: Path) -> dict[str, Any]:
    records: dict[tuple[str, int, str, tuple[str, ...]], dict[str, Any]] = {}
    key_counts: Counter[tuple[str, int, str, tuple[str, ...]]] = Counter()
    filters: Counter[str] = Counter()
    biallelic_snvs = 0
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            key = record_key(record)
            key_counts[key] += 1
            filters[filter_label(record)] += 1
            if len(record.ref) == 1 and len(record.alts or ()) == 1 and len(record.alts[0]) == 1:
                biallelic_snvs += 1
            records[key] = {
                "filter": filter_label(record),
                "pass": is_pass(record),
                "payload": payload_without_filter(record),
            }
    return {
        "path": str(path),
        "records": records,
        "record_count": sum(key_counts.values()),
        "unique_key_count": len(key_counts),
        "duplicate_key_excess": sum(count - 1 for count in key_counts.values() if count > 1),
        "biallelic_snv_count": biallelic_snvs,
        "filter_counts": dict(filters),
    }


def key_json(key: tuple[str, int, str, tuple[str, ...]]) -> list[Any]:
    return [key[0], key[1], key[2], list(key[3])]


def audit(input_vcf: Path, output_vcf: Path) -> dict[str, Any]:
    before = scan(input_vcf)
    after = scan(output_vcf)
    before_keys = set(before["records"])
    after_keys = set(after["records"])
    missing = sorted(before_keys - after_keys)
    extra = sorted(after_keys - before_keys)
    transitions: Counter[str] = Counter()
    payload_mismatches = []
    for key in sorted(before_keys & after_keys):
        left = before["records"][key]
        right = after["records"][key]
        transitions[f"{left['filter']}->{right['filter']}"] += 1
        if left["payload"] != right["payload"]:
            payload_mismatches.append(key)
    rescue = sum(
        1 for key in before_keys & after_keys
        if not before["records"][key]["pass"] and after["records"][key]["pass"]
    )
    removed = sum(
        1 for key in before_keys & after_keys
        if before["records"][key]["pass"] and not after["records"][key]["pass"]
    )
    checks = {
        "input_record_keys_unique": before["duplicate_key_excess"] == 0,
        "output_record_keys_unique": after["duplicate_key_excess"] == 0,
        "record_key_multiset_equal": not missing and not extra and before["record_count"] == after["record_count"],
        "input_all_biallelic_snv": before["record_count"] == before["biallelic_snv_count"],
        "output_all_biallelic_snv": after["record_count"] == after["biallelic_snv_count"],
        "non_filter_payload_preserved": not payload_mismatches,
        "transition_conservation": sum(transitions.values()) == before["record_count"],
    }
    return {
        "schema_version": "1.0",
        "input": {key: value for key, value in before.items() if key != "records"},
        "output": {key: value for key, value in after.items() if key != "records"},
        "filter_transition_counts": dict(sorted(transitions.items())),
        "rescued_nonpass_to_pass": rescue,
        "removed_pass_to_nonpass": removed,
        "missing_key_count": len(missing),
        "extra_key_count": len(extra),
        "payload_mismatch_count": len(payload_mismatches),
        "missing_key_examples": [key_json(key) for key in missing[:5]],
        "extra_key_examples": [key_json(key) for key in extra[:5]],
        "payload_mismatch_examples": [key_json(key) for key in payload_mismatches[:5]],
        "checks": checks,
        "pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()
    if args.output_json.exists():
        raise SystemExit(f"immutable output exists: {args.output_json}")
    result = audit(args.input_vcf, args.output_vcf)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        f"FILTER TRANSITIONS: input={result['input']['record_count']} output={result['output']['record_count']} "
        f"rescued={result['rescued_nonpass_to_pass']} removed={result['removed_pass_to_nonpass']} "
        f"pass={result['pass']} -> {args.output_json}"
    )
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
