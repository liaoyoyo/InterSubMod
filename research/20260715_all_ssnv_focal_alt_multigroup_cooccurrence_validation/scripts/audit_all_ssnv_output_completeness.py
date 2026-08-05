#!/usr/bin/env python3
"""Reconcile every all-sSNV output artifact against the frozen site universe."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping


SITE_PATTERN = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]))_(\d+)$")
EXPECTED_TOTAL = 469_849
VariantKey = tuple[str, str, int, str, str]
PositionKey = tuple[str, str, int]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def variant_key(sample: str, chrom: str, pos: int | str, ref: str, alt: str) -> VariantKey:
    return sample, chrom, int(pos), ref.upper(), alt.upper()


def expected_keys(path: Path, expected_sample: str | None = None) -> set[VariantKey]:
    keys: set[VariantKey] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = variant_key(row["sample"], row["chrom"], row["pos"], row["ref"], row["alt"])
            if expected_sample is not None and key[0] != expected_sample:
                raise RuntimeError(
                    f"Frozen site manifest sample mismatch in {path}: {key[0]!r} != {expected_sample!r}"
                )
            if key in keys:
                raise RuntimeError(f"Duplicate frozen site key in {path}: {key}")
            keys.add(key)
    return keys


def index_expected_positions(keys: set[VariantKey]) -> dict[PositionKey, VariantKey]:
    """Map path-encoded positions to alleles, rejecting identities the output cannot prove."""
    indexed: dict[PositionKey, VariantKey] = {}
    for key in keys:
        position = key[:3]
        previous = indexed.get(position)
        if previous is not None and previous != key:
            raise RuntimeError(
                "Frozen manifest has multiple alleles at one sample/chrom/pos, but InterSubMod "
                f"artifact paths do not encode REF/ALT: {previous} versus {key}"
            )
        indexed[position] = key
    return indexed


def path_site(path: Path) -> tuple[str, int]:
    for parent in path.parents:
        match = SITE_PATTERN.fullmatch(parent.name)
        if match:
            return match.group(1), int(match.group(2))
    raise ValueError(f"Cannot derive site key from {path}")


def observed_keys(
    paths: Iterable[Path],
    sample: str,
    expected_by_position: Mapping[PositionKey, VariantKey],
) -> tuple[set[VariantKey], list[VariantKey], int, list[PositionKey]]:
    keys: set[VariantKey] = set()
    duplicates: list[VariantKey] = []
    unresolved: list[PositionKey] = []
    empty = 0
    for path in paths:
        if path.stat().st_size == 0:
            empty += 1
        chrom, pos = path_site(path)
        position = (sample, chrom, pos)
        key = expected_by_position.get(position)
        if key is None:
            unresolved.append(position)
            continue
        if key in keys:
            duplicates.append(key)
        keys.add(key)
    return keys, duplicates, empty, unresolved


def compact_examples(
    values: set[VariantKey] | list[VariantKey] | list[PositionKey], limit: int = 20
) -> list[str]:
    examples: list[str] = []
    for value in sorted(values)[:limit]:
        if len(value) == 5:
            sample, chrom, pos, ref, alt = value
            examples.append(f"{sample}:{chrom}:{pos}:{ref}>{alt}")
        else:
            sample, chrom, pos = value
            examples.append(f"{sample}:{chrom}:{pos}:ALLELE_UNRESOLVED")
    return examples


def load_status(path: Path) -> dict[str, str] | None:
    if not path.exists():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0] if len(rows) == 1 else None


def receipt_path(value: Any) -> Path | None:
    if not isinstance(value, str) or not value.strip():
        return None
    return Path(value)


def same_path(value: Any, expected: Path) -> bool:
    observed = receipt_path(value)
    return observed is not None and observed.resolve() == expected.resolve()


def nonempty_absolute_file(value: Any) -> bool:
    path = receipt_path(value)
    return bool(path is not None and path.is_absolute() and path.is_file() and path.stat().st_size > 0)


def receipt_audit(
    receipt: dict[str, Any] | None,
    entry: dict[str, Any],
    sample_root: Path,
    expected_sites: int,
) -> dict[str, Any]:
    sample = entry["sample"]
    if receipt is None:
        return {"gates": {"receipt_present_and_json_object": False}, "limitations": [], "pass": False}

    validation = receipt.get("validation")
    validation = validation if isinstance(validation, dict) else {}
    input_vcf = Path(entry["all_ssnv_vcf"]["path"])
    binary = receipt_path(receipt.get("binary"))
    reference = receipt_path(receipt.get("reference"))
    binary_hash_matches = bool(
        binary is not None
        and nonempty_absolute_file(str(binary))
        and isinstance(receipt.get("binary_sha256"), str)
        and sha256(binary) == receipt["binary_sha256"]
    )
    gates = {
        "receipt_present_and_json_object": True,
        "schema": receipt.get("schema_name") == "intersubmod.all_ssnv_site_run",
        "sample": receipt.get("sample") == sample,
        "receipt_pass": receipt.get("pass") is True,
        "exit_code": receipt.get("exit_code") == 0,
        "input_vcf_path": same_path(receipt.get("input_vcf"), input_vcf),
        "input_vcf_sha256": receipt.get("input_vcf_sha256") == entry["all_ssnv_vcf"]["sha256"],
        "input_vcf_nonempty_absolute": nonempty_absolute_file(str(input_vcf)),
        "input_vcf_current_sha256": bool(
            nonempty_absolute_file(str(input_vcf))
            and sha256(input_vcf) == entry["all_ssnv_vcf"]["sha256"]
        ),
        "binary_nonempty_absolute": nonempty_absolute_file(receipt.get("binary")),
        "binary_sha256": binary_hash_matches,
        "reference_nonempty_absolute": nonempty_absolute_file(receipt.get("reference")),
        "output_dir": same_path(receipt.get("output_dir"), sample_root),
        "validation_pass": validation.get("pass") is True,
        "validation_expected_sites": validation.get("expected_vcf_sites") == expected_sites,
    }
    limitations: list[str] = []
    if "reference_sha256" in receipt:
        gates["reference_sha256"] = bool(
            reference is not None
            and nonempty_absolute_file(str(reference))
            and isinstance(receipt.get("reference_sha256"), str)
            and sha256(reference) == receipt["reference_sha256"]
        )
    else:
        limitations.append(
            "Receipt schema has no reference_sha256; reference provenance is limited to its absolute nonempty path."
        )
    return {
        "recorded": {
            "schema_name": receipt.get("schema_name"),
            "sample": receipt.get("sample"),
            "input_vcf": receipt.get("input_vcf"),
            "input_vcf_sha256": receipt.get("input_vcf_sha256"),
            "binary": receipt.get("binary"),
            "binary_sha256": receipt.get("binary_sha256"),
            "reference": receipt.get("reference"),
            "reference_sha256": receipt.get("reference_sha256"),
            "output_dir": receipt.get("output_dir"),
        },
        "gates": gates,
        "limitations": limitations,
        "pass": all(gates.values()),
    }


def audit_sample(entry: dict[str, Any], output_root: Path) -> dict[str, Any]:
    sample = entry["sample"]
    sample_root = output_root / sample
    expected = expected_keys(Path(entry["site_manifest"]["path"]), sample)
    expected_by_position = index_expected_positions(expected)
    patterns = {
        "reads": "reads.tsv",
        "methylation": "methylation.csv",
        "bernoulli": "matrix.csv",
    }
    candidates_by_role: dict[str, list[Path]] = {role: [] for role in patterns}
    for path in sample_root.rglob("*"):
        if not path.is_file():
            continue
        if path.name == "reads.tsv" and path.parent.name == "reads":
            candidates_by_role["reads"].append(path)
        elif path.name == "methylation.csv" and path.parent.name == "methylation":
            candidates_by_role["methylation"].append(path)
        elif path.name == "matrix.csv" and path.parent.name == "BERNOULLI":
            candidates_by_role["bernoulli"].append(path)
    observations: dict[str, Any] = {}
    for role in patterns:
        observed, duplicates, empty, unresolved = observed_keys(
            candidates_by_role[role], sample, expected_by_position
        )
        missing = expected - observed
        observations[role] = {
            "key_fields": ["dataset", "chrom", "pos", "ref", "alt"],
            "legacy_key_alias_fields": ["sample", "chrom", "pos", "ref", "alt"],
            "expected_keys": len(expected),
            "observed_unique_keys": len(observed),
            "duplicate_keys": len(duplicates),
            "empty_files": empty,
            "missing_keys": len(missing),
            "extra_keys": len(unresolved),
            "missing_examples": compact_examples(missing),
            "extra_examples": compact_examples(unresolved),
            "duplicate_examples": compact_examples(duplicates),
            "pass": not duplicates and not empty and not missing and not unresolved,
        }
    log_path = sample_root / "inter_sub_mod.log"
    log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
    failed_lines = [
        line for line in log_text.splitlines() if "[ERROR] Region " in line and " failed:" in line
    ]
    status_path = sample_root / "region_stratification_status.tsv"
    status = load_status(status_path)
    status_value = status.get("status") if status else None
    receipt_path = sample_root / "run_receipt.json"
    receipt: dict[str, Any] | None = None
    receipt_load_error: str | None = None
    if receipt_path.exists():
        try:
            loaded = json.loads(receipt_path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                receipt = loaded
            else:
                receipt_load_error = "Receipt JSON root is not an object"
        except (OSError, json.JSONDecodeError) as error:
            receipt_load_error = repr(error)
    provenance = receipt_audit(receipt, entry, sample_root, len(expected))
    passed = (
        all(result["pass"] for result in observations.values())
        and not failed_lines
        and status_value in {"SUCCESS", "NOT_APPLICABLE_TUMOR_ONLY"}
        and provenance["pass"]
    )
    return {
        "dataset": sample,
        "sample": sample,
        "input_site_manifest": entry["site_manifest"],
        "output_root": str(sample_root.resolve()),
        "expected_sites": len(expected),
        "artifacts": observations,
        "region_failure_count": len(failed_lines),
        "region_failure_examples": failed_lines[:20],
        "status_path": str(status_path),
        "status": status,
        "receipt_path": str(receipt_path),
        "receipt_pass": receipt.get("pass") if receipt else None,
        "receipt_load_error": receipt_load_error,
        "receipt_audit": provenance,
        "pass": passed,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    samples = [audit_sample(entry, args.output_root) for entry in manifest["samples"]]
    totals: Counter[str] = Counter()
    totals["expected_sites"] = sum(item["expected_sites"] for item in samples)
    totals["region_failures"] = sum(item["region_failure_count"] for item in samples)
    for role in ("reads", "methylation", "bernoulli"):
        totals[f"{role}_keys"] = sum(item["artifacts"][role]["observed_unique_keys"] for item in samples)
        totals[f"{role}_missing"] = sum(item["artifacts"][role]["missing_keys"] for item in samples)
        totals[f"{role}_extra"] = sum(item["artifacts"][role]["extra_keys"] for item in samples)
        totals[f"{role}_duplicates"] = sum(item["artifacts"][role]["duplicate_keys"] for item in samples)
        totals[f"{role}_empty"] = sum(item["artifacts"][role]["empty_files"] for item in samples)
    passed = totals["expected_sites"] == EXPECTED_TOTAL and all(item["pass"] for item in samples)
    payload = {
        "schema_name": "intersubmod.all_ssnv_output_reconciliation",
        "schema_version": "1.1.0",
        "created_at_utc": now_utc(),
        "status": "EXECUTION_PASS" if passed else "EXECUTION_FAIL",
        "pass_semantics": "output_completeness_and_provenance_only_not_scientific_evidence",
        "input_manifest": str(args.manifest.resolve()),
        "input_output_root": str(args.output_root.resolve()),
        "site_key_contract": ["dataset", "chrom", "pos", "ref", "alt"],
        "legacy_site_key_alias": ["sample", "chrom", "pos", "ref", "alt"],
        "totals": dict(totals),
        "datasets": samples,
        "samples": samples,
        "pass": passed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit: {args.output}")
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output.resolve()), "totals": dict(totals), "pass": passed}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
