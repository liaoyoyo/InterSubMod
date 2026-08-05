#!/usr/bin/env python3
"""Finalize an exact-PS strict-region partition root with fail-closed receipts."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_NAME = "intersubmod.strict_exact_ps_partition_run"
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


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return value


def read_tsv_gz(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(required) - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{path}: missing columns {','.join(missing)}")
        return [dict(row) for row in reader]


def write_json(path: Path, document: Mapping[str, Any]) -> None:
    if path.exists():
        raise FileExistsError(path)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(document, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    sidecar = path.with_name(path.name + ".sha256")
    with sidecar.open("x", encoding="ascii") as handle:
        handle.write(f"{sha256_path(path)}  {path.name}\n")
        handle.flush()
        os.fsync(handle.fileno())


def require_same_identity(left: Mapping[str, Any], right: Mapping[str, Any], label: str) -> None:
    if any(left.get(key) != right.get(key) for key in ("path", "size_bytes", "sha256")):
        raise ValueError(f"{label}: input identity mismatch")


def finalize(
    *,
    sample: str,
    root: Path,
    chromosomes: Sequence[str],
    strict_parity_root: Path | None,
) -> dict[str, Any]:
    root = root.resolve(strict=True)
    totals: Counter[str] = Counter()
    chrom_records = []
    for chrom in chromosomes:
        chrom_dir = root / "chromosomes" / chrom
        strict_path = chrom_dir / "strict_regions" / "receipt.json"
        partition_path = chrom_dir / "python_partition" / "receipt.json"
        comparison_path = chrom_dir / "comparison.json"
        strict = read_json(strict_path)
        partition = read_json(partition_path)
        comparison = read_json(comparison_path)
        if (
            strict.get("schema_name") != "intersubmod.strict_exact_ps_hp_regions"
            or strict.get("schema_version") != "1.1.0"
            or strict.get("all_pass") is not True
            or strict.get("scope", {}).get("dataset") != sample
            or strict.get("scope", {}).get("chrom") != chrom
        ):
            raise ValueError(f"{strict_path}: strict receipt mismatch")
        if (
            partition.get("all_pass") is not True
            or partition.get("scope", {}).get("dataset") != sample
            or partition.get("scope", {}).get("chrom") != chrom
        ):
            raise ValueError(f"{partition_path}: partition receipt mismatch")
        if comparison.get("all_pass") is not True or int(comparison.get("mismatch_count", -1)) != 0:
            raise ValueError(f"{comparison_path}: Python/C++ partition parity failed")
        strict_membership = strict.get("outputs", {}).get("membership")
        partition_membership = partition.get("inputs", {}).get("site_component_membership")
        if not isinstance(strict_membership, dict) or not isinstance(partition_membership, dict):
            raise ValueError(f"{chrom}: membership identities are missing")
        require_same_identity(strict_membership, partition_membership, f"{chrom} strict membership")

        strict_summary = strict.get("summary_by_threshold", {}).get("3", {})
        units_path = chrom_dir / "python_partition" / "units.tsv.gz"
        units = read_tsv_gz(
            units_path,
            ("dataset", "chrom", "linkage_basis", "phase_set", "threshold", "k"),
        )
        if any(
            row["dataset"] != sample
            or row["chrom"] != chrom
            or row["linkage_basis"] not in {"PS_HP1", "PS_HP2"}
            or not row["phase_set"]
            or int(row["threshold"]) != 3
            or int(row["k"]) <= 1
            for row in units
        ):
            raise ValueError(f"{units_path}: non-strict or singleton unit entered partition")
        membership_rows = read_tsv_gz(
            Path(str(strict_membership["path"])),
            ("threshold", "chrom", "pos1", "inference_role", "tree_eligible"),
        )
        selected_memberships = [row for row in membership_rows if int(row["threshold"]) == 3]
        active_unique = len({(row["chrom"], int(row["pos1"])) for row in selected_memberships})
        counts = partition.get("counts", {})
        metrics = {
            "S": int(strict["scope"]["candidate_loci_S"]),
            "unique_sites": active_unique,
            "strict_active_node_memberships": int(strict_summary["active_node_memberships"]),
            "strict_components_all": int(strict_summary["components_all"]),
            "strict_singleton_unlinked_components": int(
                strict_summary["singleton_unlinked_components"]
            ),
            "strict_tree_eligible_regions": int(strict_summary["tree_eligible_regions"]),
            "units": int(counts["eligible_units"]),
            "unit_memberships": int(counts["eligible_unit_sites"]),
            "blocks": int(counts["blocks"]),
            "python_cpp_mismatch_count": int(comparison["mismatch_count"]),
        }
        checks = {
            "strict_component_conservation": metrics["strict_components_all"]
            == metrics["strict_singleton_unlinked_components"]
            + metrics["strict_tree_eligible_regions"],
            "singletons_excluded_before_partition": metrics["units"]
            == metrics["strict_tree_eligible_regions"],
            "strict_membership_conservation": metrics["strict_active_node_memberships"]
            == metrics["unit_memberships"]
            + metrics["strict_singleton_unlinked_components"],
            "partition_units_all_k_gt1": len(units) == metrics["units"],
            "partition_cpp_parity_zero": metrics["python_cpp_mismatch_count"] == 0,
        }
        graph_parity_identity = None
        if strict_parity_root is not None:
            graph_parity_path = strict_parity_root / chrom / "parity.json"
            graph_parity = read_json(graph_parity_path)
            checks["strict_graph_python_cpp_parity"] = graph_parity.get("all_pass") is True
            graph_parity_identity = identity(graph_parity_path)
        if not all(checks.values()):
            raise ValueError(f"{chrom}: finalization check failed: {checks}")
        stage_receipt = {
            "schema_name": f"{SCHEMA_NAME}.chromosome_stage",
            "schema_version": SCHEMA_VERSION,
            "all_pass": True,
            "sample": sample,
            "chrom": chrom,
            "metrics": metrics,
            "checks": checks,
            "inputs": {
                "strict_regions": identity(strict_path),
                "partition": identity(partition_path),
                "partition_comparison": identity(comparison_path),
                "strict_graph_parity": graph_parity_identity,
            },
        }
        stage_path = chrom_dir / "stage_receipt.json"
        write_json(stage_path, stage_receipt)
        for key, value in metrics.items():
            totals[key] += value
        chrom_records.append({"chrom": chrom, "all_pass": True, "metrics": metrics, "receipt": identity(stage_path)})

    checks = {
        "all_requested_chromosomes_finalized": len(chrom_records) == len(chromosomes),
        "singleton_exclusion_conserved_genomewide": totals["strict_components_all"]
        == totals["strict_singleton_unlinked_components"]
        + totals["strict_tree_eligible_regions"],
        "partition_unit_denominator_is_tree_eligible_regions": totals["units"]
        == totals["strict_tree_eligible_regions"],
        "partition_python_cpp_mismatch_zero": totals["python_cpp_mismatch_count"] == 0,
    }
    document = {
        "schema_name": f"{SCHEMA_NAME}.run_receipt",
        "schema_version": SCHEMA_VERSION,
        "task_type": "production_deployment",
        "claim_status": "VALIDATED_PARTITION",
        "validation_evidence_eligible": True,
        "all_pass": all(checks.values()),
        "sample": sample,
        "scope": {"chromosomes": list(chromosomes), "autosomes_only": True},
        "parameters": {
            "container": "chromosome x exact PS x HP1/HP2",
            "edge_threshold": 3,
            "singleton_policy": "ABSTAIN_SINGLETON_UNLINKED",
            "max_block_size": 12,
        },
        "aggregate": dict(sorted(totals.items())),
        "checks": checks,
        "chromosomes": chrom_records,
        "claim_ceiling": (
            "Validated through strict read-linked region construction and k<=12 partition. "
            "This receipt does not claim topology, VAF ranking, or cellular clone truth."
        ),
    }
    write_json(root / "run_receipt.json", document)
    return document


def parse_chromosomes(value: str) -> tuple[str, ...]:
    result = tuple(token for token in value.split(",") if token)
    if not result or len(result) != len(set(result)):
        raise argparse.ArgumentTypeError("chromosomes must be unique and nonempty")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--chromosomes", type=parse_chromosomes, default=AUTOSOMES)
    parser.add_argument("--strict-parity-root", type=Path)
    args = parser.parse_args()
    receipt = finalize(
        sample=args.sample,
        root=args.root,
        chromosomes=args.chromosomes,
        strict_parity_root=(
            args.strict_parity_root.resolve(strict=True)
            if args.strict_parity_root is not None
            else None
        ),
    )
    print(json.dumps({"all_pass": receipt["all_pass"], "aggregate": receipt["aggregate"]}, sort_keys=True))
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
