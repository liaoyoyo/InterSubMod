#!/usr/bin/env python3
"""Aggregate chromosome-level Python/C++ strict graph parity receipts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Sequence


SCHEMA_NAME = "intersubmod.strict_endpoint_python_cpp_genome_parity"
SCHEMA_VERSION = "1.0.0"
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256_path(path)}


def aggregate(
    *, input_root: Path, output: Path, dataset: str, chromosomes: Sequence[str]
) -> dict[str, Any]:
    input_root = input_root.resolve(strict=True)
    totals = {
        "observed_edge_rows": 0,
        "components_all": 0,
        "edge_mismatch_preview_count": 0,
        "component_mismatch_preview_count": 0,
    }
    inputs = {}
    checks_by_chrom = {}
    for chrom in chromosomes:
        path = input_root / chrom / "parity.json"
        receipt = json.loads(path.read_text(encoding="utf-8"))
        if (
            receipt.get("schema_name") != "intersubmod.strict_endpoint_python_cpp_parity"
            or receipt.get("schema_version") != "1.0.0"
        ):
            raise ValueError(f"{path}: schema mismatch")
        counts = receipt.get("counts", {})
        if counts.get("python_edges") != counts.get("cpp_edges"):
            raise ValueError(f"{path}: edge count mismatch")
        if counts.get("python_components") != counts.get("cpp_components"):
            raise ValueError(f"{path}: component count mismatch")
        totals["observed_edge_rows"] += int(counts["python_edges"])
        totals["components_all"] += int(counts["python_components"])
        totals["edge_mismatch_preview_count"] += int(counts["edge_mismatch_preview_count"])
        totals["component_mismatch_preview_count"] += int(
            counts["component_mismatch_preview_count"]
        )
        checks_by_chrom[chrom] = receipt.get("all_pass") is True and all(
            receipt.get("checks", {}).values()
        )
        inputs[chrom] = identity(path)
    checks = {
        "all_chromosomes_present": len(inputs) == len(chromosomes),
        "all_chromosome_parity_receipts_pass": all(checks_by_chrom.values()),
        "edge_mismatch_zero": totals["edge_mismatch_preview_count"] == 0,
        "component_and_role_mismatch_zero": totals["component_mismatch_preview_count"] == 0,
    }
    document = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": all(checks.values()),
        "scope": {"dataset": dataset, "chromosomes": list(chromosomes), "threshold": 3},
        "counts": totals,
        "checks": checks,
        "checks_by_chromosome": checks_by_chrom,
        "inputs": inputs,
    }
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(document, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    return document


def parse_chromosomes(value: str) -> tuple[str, ...]:
    result = tuple(token for token in value.split(",") if token)
    if not result or len(result) != len(set(result)):
        raise argparse.ArgumentTypeError("chromosomes must be unique and nonempty")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--chromosomes", type=parse_chromosomes, default=AUTOSOMES)
    args = parser.parse_args()
    result = aggregate(
        input_root=args.input_root,
        output=args.output,
        dataset=args.dataset,
        chromosomes=args.chromosomes,
    )
    print(json.dumps({"all_pass": result["all_pass"], "counts": result["counts"]}, sort_keys=True))
    return 0 if result["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
