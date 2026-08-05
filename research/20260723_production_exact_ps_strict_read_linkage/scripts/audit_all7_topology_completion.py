#!/usr/bin/env python3
"""Audit the latest all-dataset completion boundary for strict topology.

This audit intentionally separates:

* L1: exact-PS × HP strict read-linkage regions;
* L2: older or exploratory candidate mutation-state-tree observations;
* L3: the current production v4 strict directed-topology stage; and
* L4: cellular clone count, exact parent-child lineage, and fused trees.

The script does not launch topology reconstruction.  It inventories bounded
search roots and fails closed on the receipts used to describe current status.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence


CANONICAL_DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
SCHEMA_NAME = "intersubmod.all7_topology_completion_audit"
SCHEMA_VERSION = "1.0.0"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return value


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def write_tsv_exclusive(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError("topology completion TSV cannot be empty")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def parse_dataset_root(value: str) -> tuple[str, Path]:
    dataset, separator, raw_path = value.partition("=")
    if separator != "=" or dataset not in CANONICAL_DATASETS:
        raise argparse.ArgumentTypeError(
            "--dataset-root must be CANONICAL_DATASET=/absolute/path"
        )
    path = Path(raw_path)
    if not path.is_absolute():
        raise argparse.ArgumentTypeError("--dataset-root path must be absolute")
    return dataset, path


def verify_l1_roots(dataset_roots: Mapping[str, Path]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if set(dataset_roots) != set(CANONICAL_DATASETS):
        raise ValueError("strict L1 roots must contain exactly the seven canonical datasets")
    rows: list[dict[str, Any]] = []
    inputs: dict[str, Any] = {}
    for dataset in CANONICAL_DATASETS:
        root = dataset_roots[dataset].resolve(strict=True)
        summary_path = root / "summary" / "summary.json"
        summary = read_json(summary_path)
        checks = summary.get("checks")
        if (
            summary.get("schema_name")
            != "intersubmod.strict_exact_ps_hp_genome_summary"
            or summary.get("all_pass") is not True
            or not isinstance(checks, Mapping)
            or not checks
            or any(value is not True for value in checks.values())
        ):
            raise ValueError(f"{summary_path}: strict L1 summary contract failed")
        summary_inputs = summary.get("inputs")
        if (
            not isinstance(summary_inputs, Mapping)
            or set(summary_inputs) != set(AUTOSOMES)
        ):
            raise ValueError(f"{summary_path}: strict L1 chromosome grid is not chr1-22")
        aggregate = summary.get("aggregate")
        if not isinstance(aggregate, Mapping):
            raise ValueError(f"{summary_path}: strict L1 aggregate is absent")
        w = aggregate.get("tree_eligible_read_linked_regions")
        if not isinstance(w, int) or isinstance(w, bool) or w < 0:
            raise ValueError(f"{summary_path}: strict L1 W is invalid")
        rows.append(
            {
                "dataset": dataset,
                "strict_linkage_status": "COMPLETE_22_OF_22",
                "strict_W": w,
            }
        )
        inputs[dataset] = {
            "root": str(root),
            "summary": identity(summary_path),
        }
    return rows, inputs


def discover_v4_receipts(search_roots: Sequence[Path]) -> list[tuple[Path, dict[str, Any]]]:
    receipts: list[tuple[Path, dict[str, Any]]] = []
    seen: set[Path] = set()
    for raw_root in search_roots:
        root = raw_root.resolve(strict=True)
        for pattern in (
            "*layered_v4_strict*/run_receipt.json",
            "*layered_v4_strict*/**/run_receipt.json",
        ):
            for path in sorted(root.glob(pattern)):
                resolved = path.resolve()
                if resolved in seen or not resolved.is_file():
                    continue
                seen.add(resolved)
                value = read_json(resolved)
                if (
                    value.get("schema_name")
                    == "intersubmod.layered_v4_strict_production_run"
                ):
                    receipts.append((resolved, value))
    return sorted(receipts, key=lambda item: (item[0].stat().st_mtime_ns, str(item[0])))


def topology_identity_for_dataset(receipt: Mapping[str, Any], dataset: str) -> Mapping[str, Any] | None:
    sample_receipts = receipt.get("sample_receipts")
    if not isinstance(sample_receipts, Mapping):
        return None
    sample = sample_receipts.get(dataset)
    if not isinstance(sample, Mapping):
        return None
    topology = sample.get("topology")
    return topology if isinstance(topology, Mapping) else None


def is_full_production_topology(receipt: Mapping[str, Any]) -> bool:
    scope = receipt.get("scope")
    return (
        receipt.get("all_pass") is True
        and receipt.get("stage_through") == "topology"
        and receipt.get("task_type") == "B_comprehensive_validation"
        and receipt.get("validation_evidence_eligible") is True
        and isinstance(scope, Mapping)
        and scope.get("datasets") == list(CANONICAL_DATASETS)
        and scope.get("chromosomes") == list(AUTOSOMES)
        and all(
            topology_identity_for_dataset(receipt, dataset) is not None
            for dataset in CANONICAL_DATASETS
        )
    )


def verify_hcc_pilot_receipt(path: Path) -> dict[str, Any]:
    receipt = read_json(path)
    evidence = receipt.get("evidence_tier")
    if (
        receipt.get("schema_name")
        != "intersubmod.exact_ps_topology_observation.receipt"
        or receipt.get("all_technical_checks_pass") is not True
        or receipt.get("sample") != "HCC1395"
        or not isinstance(evidence, Mapping)
        or evidence.get("publication_or_cohort_final_eligible") is not False
        or evidence.get("upstream_partition_receipt", {}).get(
            "validation_evidence_eligible"
        )
        is not False
    ):
        raise ValueError(f"{path}: HCC1395 exploratory topology receipt mismatch")
    return receipt


def verify_legacy_summary(path: Path) -> dict[str, Any]:
    summary = read_json(path)
    canonical = summary.get("canonical")
    samples = canonical.get("samples") if isinstance(canonical, Mapping) else None
    if (
        summary.get("schema_name") != "intersubmod.current_layered_topology_summary"
        or summary.get("all_pass") is not True
        or summary.get("scope") != "7 datasets / 6 biological samples / chr1-22"
        or summary.get("claim_scope")
        != "regional mutation-state candidate trees; not confirmed cell clones"
        or not isinstance(samples, list)
        or {row.get("sample") for row in samples} != set(CANONICAL_DATASETS)
        or any(row.get("pass") is not True for row in samples)
    ):
        raise ValueError(f"{path}: legacy all-dataset topology summary mismatch")
    return summary


def build_audit(
    *,
    dataset_roots: Mapping[str, Path],
    search_roots: Sequence[Path],
    hcc_pilot_receipt_path: Path,
    legacy_summary_path: Path,
) -> dict[str, Any]:
    l1_rows, l1_inputs = verify_l1_roots(dataset_roots)
    v4_receipts = discover_v4_receipts(search_roots)
    if not v4_receipts:
        raise ValueError("no v4 strict run receipt was found in the bounded search roots")
    hcc_pilot = verify_hcc_pilot_receipt(hcc_pilot_receipt_path)
    legacy = verify_legacy_summary(legacy_summary_path)

    full_receipts = [
        (path, receipt)
        for path, receipt in v4_receipts
        if is_full_production_topology(receipt)
    ]
    latest_path, latest_receipt = v4_receipts[-1]
    rows: list[dict[str, Any]] = []
    for base in l1_rows:
        dataset = base["dataset"]
        partial_topology = [
            (path, receipt)
            for path, receipt in v4_receipts
            if topology_identity_for_dataset(receipt, dataset) is not None
        ]
        complete_topology = [
            (path, receipt)
            for path, receipt in full_receipts
            if topology_identity_for_dataset(receipt, dataset) is not None
        ]
        if complete_topology:
            strict_topology_status = "COMPLETE_PRODUCTION"
            topology_reason = "full-scope v4 strict topology receipt present"
        elif partial_topology:
            strict_topology_status = "PARTIAL_NOT_PRODUCTION_VALIDATED"
            topology_reason = "v4 topology output exists only outside full-scope eligible production"
        elif (
            dataset in latest_receipt.get("scope", {}).get("datasets", [])
            and topology_identity_for_dataset(latest_receipt, dataset) is None
        ):
            strict_topology_status = "NOT_RUN_LATEST_V4_PARTITION_ONLY"
            topology_reason = (
                f"latest v4 receipt ends at {latest_receipt.get('stage_through')}; "
                "sample topology identity is null"
            )
        else:
            strict_topology_status = "NOT_RUN"
            topology_reason = "no v4 strict topology receipt found in bounded search roots"
        rows.append(
            {
                **base,
                "legacy_candidate_tree_status": "LEGACY_REFERENCE_ONLY",
                "exploratory_exact_ps_topology_status": (
                    "HCC1395_TECHNICAL_PILOT_NOT_ELIGIBLE"
                    if dataset == "HCC1395"
                    else "NOT_AVAILABLE"
                ),
                "strict_directed_topology_status": strict_topology_status,
                "strict_topology_reason": topology_reason,
                "vaf_or_read_likelihood_ranking_status": "NOT_PRODUCTION_VALIDATED",
                "cellular_clone_count_status": "NOT_VALIDATED",
                "exact_cellular_parent_child_status": "NOT_VALIDATED",
                "cross_hp_or_technical_fusion_status": "NOT_VALIDATED",
            }
        )

    strict_topology_complete = sum(
        row["strict_directed_topology_status"] == "COMPLETE_PRODUCTION"
        for row in rows
    )
    checks = {
        "strict_l1_seven_dataset_chr1_22_grid_verified": len(rows) == 7,
        "bounded_v4_receipt_inventory_nonempty": bool(v4_receipts),
        "no_full_scope_production_v4_topology_receipt_found": not full_receipts,
        "latest_v4_receipt_is_partition_only": (
            latest_receipt.get("stage_through") == "partition"
            and all(
                topology_identity_for_dataset(latest_receipt, dataset) is None
                for dataset in latest_receipt.get("scope", {}).get("datasets", [])
            )
        ),
        "hcc1395_exact_ps_topology_is_exploratory_not_final": (
            hcc_pilot["evidence_tier"]["publication_or_cohort_final_eligible"]
            is False
        ),
        "legacy_all7_claim_is_candidate_mutation_state_not_clone": (
            legacy.get("claim_scope")
            == "regional mutation-state candidate trees; not confirmed cell clones"
        ),
        "strict_topology_production_complete_count_is_zero": (
            strict_topology_complete == 0
        ),
        "cellular_clone_parent_child_and_fusion_not_validated_for_all7": all(
            row["cellular_clone_count_status"] == "NOT_VALIDATED"
            and row["exact_cellular_parent_child_status"] == "NOT_VALIDATED"
            and row["cross_hp_or_technical_fusion_status"] == "NOT_VALIDATED"
            for row in rows
        ),
    }
    if not all(checks.values()):
        raise ValueError(f"topology completion audit checks failed: {checks}")

    return {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "scope": {
            "task_type": "B_comprehensive_validation",
            "datasets": list(CANONICAL_DATASETS),
            "chromosomes": list(AUTOSOMES),
            "search_roots": [str(path.resolve()) for path in search_roots],
            "search_pattern": "*layered_v4_strict*/**/run_receipt.json",
            "absence_claim_boundary": (
                "NOT_RUN means no qualifying v4 strict topology receipt was found "
                "inside the declared bounded search roots at audit time."
            ),
        },
        "summary": {
            "strict_linkage_complete_datasets": 7,
            "strict_directed_topology_production_complete_datasets": strict_topology_complete,
            "cellular_clone_count_validated_datasets": 0,
            "exact_cellular_parent_child_validated_datasets": 0,
            "fused_tree_validated_datasets": 0,
            "v4_receipts_in_inventory": len(v4_receipts),
            "full_scope_v4_topology_receipts": len(full_receipts),
        },
        "rows": rows,
        "checks": checks,
        "inputs": {
            "strict_l1": l1_inputs,
            "v4_receipts": [identity(path) for path, _ in v4_receipts],
            "latest_v4_receipt": identity(latest_path),
            "hcc1395_exploratory_topology_receipt": identity(
                hcc_pilot_receipt_path
            ),
            "legacy_all7_topology_summary": identity(legacy_summary_path),
        },
        "verdict": (
            "L1_STRICT_READ_LINKAGE_COMPLETE_7_OF_7__"
            "L3_PRODUCTION_STRICT_DIRECTED_TOPOLOGY_COMPLETE_0_OF_7__"
            "L4_CLONE_PARENT_CHILD_FUSION_VALIDATED_0_OF_7"
        ),
        "claim_ceiling": (
            "All seven datasets have verified exact-PS × HP strict read-linkage "
            "regions. No dataset currently has a full-scope eligible v4 strict "
            "directed-topology receipt. Clone counts, exact cellular parent-child "
            "edges, and cross-HP/technical fused trees remain unvalidated."
        ),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", action="append", required=True, type=parse_dataset_root)
    parser.add_argument("--search-root", action="append", required=True, type=Path)
    parser.add_argument("--hcc-pilot-receipt", required=True, type=Path)
    parser.add_argument("--legacy-summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tsv", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    dataset_roots = dict(args.dataset_root)
    if len(dataset_roots) != len(args.dataset_root):
        raise ValueError("duplicate --dataset-root dataset")
    output_path = args.output.resolve()
    tsv_path = args.tsv.resolve()
    receipt_path = args.receipt.resolve()
    sidecar_path = receipt_path.with_name(f"{receipt_path.name}.sha256")
    for path in (output_path, tsv_path, receipt_path, sidecar_path):
        if path.exists():
            raise ValueError(f"output must not exist: {path}")

    audit = build_audit(
        dataset_roots=dataset_roots,
        search_roots=args.search_root,
        hcc_pilot_receipt_path=args.hcc_pilot_receipt.resolve(strict=True),
        legacy_summary_path=args.legacy_summary.resolve(strict=True),
    )
    write_json_exclusive(output_path, audit)
    write_tsv_exclusive(tsv_path, audit["rows"])
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "checks": audit["checks"],
        "inputs": audit["inputs"],
        "outputs": {
            output_path.name: identity(output_path),
            tsv_path.name: identity(tsv_path),
        },
    }
    write_json_exclusive(receipt_path, receipt)
    with sidecar_path.open("x", encoding="ascii") as handle:
        handle.write(f"{sha256_path(receipt_path)}  {receipt_path.name}\n")
        handle.flush()
        os.fsync(handle.fileno())
    print(
        json.dumps(
            {
                "all_pass": audit["all_pass"],
                "strict_linkage_complete": audit["summary"][
                    "strict_linkage_complete_datasets"
                ],
                "strict_topology_complete": audit["summary"][
                    "strict_directed_topology_production_complete_datasets"
                ],
                "clone_count_validated": audit["summary"][
                    "cellular_clone_count_validated_datasets"
                ],
                "output": str(output_path),
                "receipt": str(receipt_path),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
