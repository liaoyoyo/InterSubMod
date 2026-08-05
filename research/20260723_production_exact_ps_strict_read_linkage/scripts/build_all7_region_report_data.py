#!/usr/bin/env python3
"""Build an auditable all-dataset data package for the strict region HTML report.

The input unit is one completed strict-region root per canonical dataset.  Each
root must contain ``summary/summary.json`` and 22 chromosome receipts.  The
script independently re-reads the primary-threshold component, membership and
edge files so that the report does not merely repeat pre-computed totals.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
from typing import Any, Iterable, Mapping, Sequence


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
PRIMARY_THRESHOLD = 3
VALID_LINKAGE_BASIS = {"PS_HP1": "1", "PS_HP2": "2"}
MISSING_PHASE_SET_TOKENS = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}
SCHEMA_NAME = "intersubmod.strict_region_all7_report_data"
SCHEMA_VERSION = "1.1.0"
THRESHOLD_SUMMARY_FIELDS = (
    "active_node_memberships",
    "components_all",
    "containers",
    "k1_regions",
    "k_gt1_regions",
    "regions",
    "retained_endpoint_edges",
    "singleton_unlinked_components",
    "tree_eligible_regions",
)
REQUIRED_SUMMARY_CHECKS = frozenset(
    {
        "all_components_independently_connected",
        "all_requested_chromosomes_present",
        "all_retained_edges_have_direct_support_at_threshold",
        "component_membership_mass_conserved",
        "component_partition_conserved",
        "cross_hp_violations_zero",
        "cross_ps_violations_zero",
        "singletons_excluded_from_tree_regions",
    }
)
REQUIRED_CHROM_RECEIPT_CHECKS = frozenset(
    {
        "all_reported_edges_have_direct_fixed_endpoint_support",
        "component_membership_mass_conserved",
        "distance_not_used_for_connectivity",
        "edge_state_mass_conserved",
        "every_active_node_has_one_component_per_threshold",
        "every_multisite_component_connected_by_retained_edges",
        "only_exact_known_ps_hp_primary",
        "primary_components_are_multisite_and_read_linked",
        "primary_threshold_present",
        "singleton_components_never_primary",
        "unique_molecule_ids",
    }
)
REQUIRED_EXTRACTION_RECEIPT_CHECKS = frozenset(
    {
        "all_linkage_units_and_thresholds_conserve_their_active_sites",
        "alt_mass_conserved",
        "canonical_molecule_unique",
        "duplicate_identity_collapse_conserved",
        "duplicate_identity_policy_matches_manifest",
        "eligible_alignment_rows_written",
        "every_eligible_alignment_exact_sidecar_joined",
        "fixed_ra_mass_conserved",
        "missing_ps_never_primary",
        "multi_region_alignment_unique",
        "primary_components_only_use_known_ps",
        "primary_membership_is_active_fixed_ra_only",
        "raw_alignment_class_conserved",
        "raw_filter_funnel_conserved",
        "samtools_exit_zero",
        "site_call_reason_mass_conserved",
        "site_catalog_cardinality",
    }
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


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


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return value


def read_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def verify_identity(spec: Mapping[str, Any], path: Path) -> None:
    observed = identity(path)
    keys = ("path", "size_bytes", "sha256")
    if any(spec.get(key) != observed[key] for key in keys):
        raise ValueError(f"{path}: identity mismatch")


def verify_path_sha_provenance(spec: Any, *, label: str) -> None:
    """Require a current, content-addressed provenance file."""
    if not isinstance(spec, Mapping):
        raise ValueError(f"{label}: provenance identity is absent")
    raw_path = spec.get("path")
    declared_sha = spec.get("sha256")
    if (
        not isinstance(raw_path, str)
        or not raw_path
        or not Path(raw_path).is_absolute()
        or not isinstance(declared_sha, str)
        or len(declared_sha) != 64
        or any(char not in "0123456789abcdef" for char in declared_sha)
    ):
        raise ValueError(f"{label}: provenance must contain absolute path + sha256")
    path = Path(raw_path)
    if not path.is_file() or sha256_path(path) != declared_sha:
        raise ValueError(f"{label}: provenance path/sha256 mismatch")


def validate_extraction_receipt(
    *,
    dataset: str,
    chrom: str,
    strict_inputs: Any,
) -> dict[str, Any]:
    """Validate the extraction receipt that produced both strict graph inputs."""
    if not isinstance(strict_inputs, Mapping):
        raise ValueError(f"{dataset}/{chrom}: strict receipt inputs are absent")
    bound_inputs: dict[str, Mapping[str, Any]] = {}
    input_paths: dict[str, Path] = {}
    for name in ("molecule_calls", "site_catalog"):
        spec = strict_inputs.get(name)
        if (
            not isinstance(spec, Mapping)
            or not isinstance(spec.get("path"), str)
            or not Path(spec["path"]).is_absolute()
            or not isinstance(spec.get("size_bytes"), int)
            or isinstance(spec.get("size_bytes"), bool)
            or spec["size_bytes"] < 0
            or not isinstance(spec.get("sha256"), str)
            or len(spec["sha256"]) != 64
        ):
            raise ValueError(
                f"{dataset}/{chrom}: strict {name} identity is incomplete"
            )
        bound_inputs[name] = spec
        input_paths[name] = Path(spec["path"])
    extraction_dirs = {path.parent for path in input_paths.values()}
    if len(extraction_dirs) != 1:
        raise ValueError(
            f"{dataset}/{chrom}: strict graph inputs do not share one extraction dir"
        )
    extraction_dir = next(iter(extraction_dirs))
    receipt_path = extraction_dir / "receipt.json"
    sidecar_path = extraction_dir / "receipt.json.sha256"
    if not receipt_path.is_file() or not sidecar_path.is_file():
        raise ValueError(
            f"{dataset}/{chrom}: extraction receipt or checksum sidecar is missing"
        )
    expected_sidecar = f"{sha256_path(receipt_path)}  receipt.json"
    if sidecar_path.read_text(encoding="ascii").strip() != expected_sidecar:
        raise ValueError(f"{sidecar_path}: extraction receipt checksum mismatch")

    receipt = read_json(receipt_path)
    if (
        receipt.get("schema_name")
        != "intersubmod.lossless_read_linkage_chromosome_receipt"
        or receipt.get("schema_version") != "1.3.0"
        or receipt.get("all_pass") is not True
        or receipt.get("scope", {}).get("dataset") != dataset
        or receipt.get("scope", {}).get("chrom") != chrom
    ):
        raise ValueError(f"{receipt_path}: extraction receipt contract mismatch")
    checks = receipt.get("checks")
    if (
        not isinstance(checks, Mapping)
        or set(checks) != REQUIRED_EXTRACTION_RECEIPT_CHECKS
        or any(value is not True for value in checks.values())
    ):
        raise ValueError(
            f"{receipt_path}: canonical extraction checks are missing, extra or failed"
        )

    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping):
        raise ValueError(f"{receipt_path}: extraction outputs are absent")
    for name, strict_spec in bound_inputs.items():
        output_name = input_paths[name].name
        if outputs.get(output_name) != strict_spec:
            raise ValueError(
                f"{receipt_path}: extraction {name} identity differs from strict input"
            )

    provenance = receipt.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError(f"{receipt_path}: extraction provenance is absent")
    for name in ("extractor", "manifest"):
        verify_path_sha_provenance(
            provenance.get(name),
            label=f"{receipt_path}: provenance.{name}",
        )
    return identity(receipt_path)


def pct(numerator: int | float, denominator: int | float) -> float | None:
    return round(100.0 * numerator / denominator, 4) if denominator else None


def percentile(values: Sequence[int], quantile: float) -> int | None:
    """Return a deterministic nearest-rank percentile."""
    if not values:
        return None
    ordered = sorted(values)
    rank = max(1, math.ceil(quantile * len(ordered)))
    return ordered[rank - 1]


def require_exact_primary_container(
    *,
    dataset: str,
    chrom: str,
    row: Mapping[str, str],
    path: Path,
) -> tuple[str, str]:
    """Validate and return the exact-PS primary linkage container."""
    if row.get("dataset") != dataset or row.get("chrom") != chrom:
        raise ValueError(f"{path}: dataset/chrom row scope mismatch")
    basis = row.get("linkage_basis", "")
    if basis not in VALID_LINKAGE_BASIS:
        raise ValueError(f"{path}: invalid primary linkage_basis={basis!r}")
    phase_set = row.get("phase_set", "")
    if phase_set.strip().upper() in MISSING_PHASE_SET_TOKENS:
        raise ValueError(f"{path}: missing/non-exact phase_set in primary row")
    return basis, phase_set


def k_band(k: int) -> str:
    if k <= 5:
        return str(k)
    if k <= 8:
        return "6–8"
    if k <= 12:
        return "9–12"
    return ">12"


def span_band(span_bp: int) -> str:
    if span_bp <= 10_000:
        return "≤10 kb"
    if span_bp <= 25_000:
        return "10–25 kb"
    if span_bp <= 50_000:
        return "25–50 kb"
    if span_bp <= 100_000:
        return "50–100 kb"
    return ">100 kb"


def make_output_dir(path: Path) -> Path:
    if path.exists():
        if not path.is_dir() or next(path.iterdir(), None) is not None:
            raise ValueError(f"output directory must be new or empty: {path}")
    else:
        path.mkdir(parents=True, exist_ok=False)
    return path.resolve()


def write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def load_topology_completion_audit(
    path: Path,
    *,
    loaded_datasets: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Verify and load the separately inventoried L1-L4 completion boundary."""
    audit_path = path.resolve(strict=True)
    audit = read_json(audit_path)
    if (
        audit.get("schema_name")
        != "intersubmod.all7_topology_completion_audit"
        or audit.get("schema_version") != "1.0.0"
        or audit.get("all_pass") is not True
    ):
        raise ValueError(f"{audit_path}: topology completion audit contract mismatch")
    scope = audit.get("scope")
    rows = audit.get("rows")
    summary = audit.get("summary")
    checks = audit.get("checks")
    if (
        not isinstance(scope, Mapping)
        or scope.get("task_type") != "B_comprehensive_validation"
        or scope.get("datasets") != list(CANONICAL_DATASETS)
        or scope.get("chromosomes") != list(AUTOSOMES)
        or not isinstance(rows, list)
        or [row.get("dataset") for row in rows] != list(CANONICAL_DATASETS)
        or not isinstance(summary, Mapping)
        or not isinstance(checks, Mapping)
        or not checks
        or any(value is not True for value in checks.values())
    ):
        raise ValueError(f"{audit_path}: topology completion audit scope/check mismatch")

    dataset_by_name = {
        item["dataset_row"]["dataset"]: item for item in loaded_datasets
    }
    audit_l1 = audit.get("inputs", {}).get("strict_l1")
    if not isinstance(audit_l1, Mapping) or set(audit_l1) != set(CANONICAL_DATASETS):
        raise ValueError(f"{audit_path}: topology audit strict-L1 provenance is incomplete")
    for row in rows:
        dataset = row["dataset"]
        current = dataset_by_name[dataset]
        if (
            row.get("strict_linkage_status") != "COMPLETE_22_OF_22"
            or row.get("strict_W")
            != current["dataset_row"]["tree_eligible_W"]
            or row.get("strict_directed_topology_status")
            not in {
                "NOT_RUN",
                "NOT_RUN_LATEST_V4_PARTITION_ONLY",
                "PARTIAL_NOT_PRODUCTION_VALIDATED",
                "COMPLETE_PRODUCTION",
            }
            or row.get("cellular_clone_count_status") != "NOT_VALIDATED"
            or row.get("exact_cellular_parent_child_status") != "NOT_VALIDATED"
            or row.get("cross_hp_or_technical_fusion_status") != "NOT_VALIDATED"
            or audit_l1[dataset].get("summary") != current["summary_identity"]
        ):
            raise ValueError(
                f"{audit_path}: {dataset} topology status does not match strict L1"
            )
    if (
        summary.get("strict_linkage_complete_datasets") != 7
        or summary.get(
            "strict_directed_topology_production_complete_datasets"
        )
        != 0
        or summary.get("cellular_clone_count_validated_datasets") != 0
        or summary.get("exact_cellular_parent_child_validated_datasets") != 0
        or summary.get("fused_tree_validated_datasets") != 0
    ):
        raise ValueError(
            f"{audit_path}: expected current completion boundary is not preserved"
        )

    receipt_path = audit_path.parent / "receipt.json"
    sidecar_path = receipt_path.with_name("receipt.json.sha256")
    receipt = read_json(receipt_path)
    expected_sidecar = f"{sha256_path(receipt_path)}  receipt.json"
    if sidecar_path.read_text(encoding="ascii").strip() != expected_sidecar:
        raise ValueError(f"{sidecar_path}: topology audit checksum mismatch")
    output_identity = receipt.get("outputs", {}).get(audit_path.name)
    if (
        receipt.get("schema_name")
        != "intersubmod.all7_topology_completion_audit.receipt"
        or receipt.get("schema_version") != "1.0.0"
        or receipt.get("all_pass") is not True
        or not isinstance(output_identity, Mapping)
    ):
        raise ValueError(f"{receipt_path}: topology audit receipt contract mismatch")
    verify_identity(output_identity, audit_path)
    return {
        "audit": audit,
        "audit_identity": identity(audit_path),
        "receipt_identity": identity(receipt_path),
        "receipt_sidecar_identity": identity(sidecar_path),
    }


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fields = list(rows[0])
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def parse_dataset_root(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("dataset root must be DATASET=/absolute/path")
    dataset, raw_path = value.split("=", 1)
    if dataset not in CANONICAL_DATASETS:
        raise argparse.ArgumentTypeError(f"unsupported dataset: {dataset}")
    path = Path(raw_path)
    if not path.is_absolute():
        raise argparse.ArgumentTypeError("dataset root must be absolute")
    return dataset, path


def load_dataset(dataset: str, root: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    summary_path = root / "summary" / "summary.json"
    summary = read_json(summary_path)
    if (
        summary.get("schema_name") != "intersubmod.strict_exact_ps_hp_genome_summary"
        or summary.get("schema_version") != "1.0.0"
        or summary.get("all_pass") is not True
        or summary.get("scope", {}).get("dataset") != dataset
        or summary.get("scope", {}).get("chromosomes") != list(AUTOSOMES)
        or summary.get("scope", {}).get("primary_threshold") != PRIMARY_THRESHOLD
    ):
        raise ValueError(f"{summary_path}: full-scope strict summary contract mismatch")
    summary_checks = summary.get("checks")
    if (
        not isinstance(summary_checks, Mapping)
        or set(summary_checks) != REQUIRED_SUMMARY_CHECKS
        or any(value is not True for value in summary_checks.values())
    ):
        raise ValueError(
            f"{summary_path}: canonical summary checks are missing, extra or failed"
        )
    summary_inputs = summary.get("inputs")
    if (
        not isinstance(summary_inputs, Mapping)
        or set(summary_inputs) != set(AUTOSOMES)
    ):
        raise ValueError(f"{summary_path}: chromosome input identity grid is incomplete")

    active_loci: set[tuple[str, int]] = set()
    loci_in_w: set[tuple[str, int]] = set()
    loci_in_singleton: set[tuple[str, int]] = set()
    exact_ps: set[tuple[str, str]] = set()
    exact_ps_hp: set[tuple[str, str, str]] = set()
    dataset_spans: list[int] = []
    dataset_ks: list[int] = []
    dataset_direct_edge_distances: list[int] = []
    dataset_counts: Counter[str] = Counter()
    hp_counts: dict[str, Counter[str]] = defaultdict(Counter)
    k_counts: Counter[int] = Counter()
    k_band_counts: Counter[str] = Counter()
    span_band_counts: Counter[str] = Counter()
    chromosome_rows: list[dict[str, Any]] = []
    input_receipts: dict[str, Any] = {}
    extraction_receipts: dict[str, Any] = {}
    receipt_threshold_totals: dict[int, Counter[str]] = {
        threshold: Counter() for threshold in (1, 2, 3, 5)
    }

    for chrom in AUTOSOMES:
        receipt_path = root / "chromosomes" / chrom / "receipt.json"
        sidecar_path = receipt_path.with_name("receipt.json.sha256")
        receipt = read_json(receipt_path)
        expected_sidecar = f"{sha256_path(receipt_path)}  receipt.json"
        if sidecar_path.read_text(encoding="ascii").strip() != expected_sidecar:
            raise ValueError(f"{sidecar_path}: checksum mismatch")
        if (
            receipt.get("schema_name") != "intersubmod.strict_exact_ps_hp_regions"
            or receipt.get("schema_version") != "1.1.0"
            or receipt.get("all_pass") is not True
            or receipt.get("scope", {}).get("dataset") != dataset
            or receipt.get("scope", {}).get("chrom") != chrom
        ):
            raise ValueError(f"{receipt_path}: chromosome receipt contract mismatch")
        receipt_checks = receipt.get("checks")
        if (
            not isinstance(receipt_checks, Mapping)
            or set(receipt_checks) != REQUIRED_CHROM_RECEIPT_CHECKS
            or any(value is not True for value in receipt_checks.values())
        ):
            raise ValueError(
                f"{receipt_path}: canonical chromosome checks are missing, extra or failed"
            )
        extraction_receipts[chrom] = validate_extraction_receipt(
            dataset=dataset,
            chrom=chrom,
            strict_inputs=receipt.get("inputs"),
        )
        current_receipt_identity = identity(receipt_path)
        if summary_inputs.get(chrom) != current_receipt_identity:
            raise ValueError(
                f"{summary_path}: {chrom} input identity differs from current receipt"
            )
        input_receipts[chrom] = current_receipt_identity
        receipt_thresholds = receipt.get("summary_by_threshold")
        if (
            not isinstance(receipt_thresholds, Mapping)
            or set(receipt_thresholds) != {"1", "2", "3", "5"}
        ):
            raise ValueError(f"{receipt_path}: threshold summary grid mismatch")
        for threshold in (1, 2, 3, 5):
            threshold_summary = receipt_thresholds[str(threshold)]
            if (
                not isinstance(threshold_summary, Mapping)
                or not set(THRESHOLD_SUMMARY_FIELDS).issubset(threshold_summary)
            ):
                raise ValueError(
                    f"{receipt_path}: threshold={threshold} summary fields missing"
                )
            for field in THRESHOLD_SUMMARY_FIELDS:
                value = threshold_summary[field]
                if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                    raise ValueError(
                        f"{receipt_path}: threshold={threshold} {field} is invalid"
                    )
                receipt_threshold_totals[threshold][field] += value

        outputs = receipt["outputs"]
        component_path = Path(outputs["components"]["path"])
        membership_path = Path(outputs["membership"]["path"])
        edge_path = Path(outputs["edges"]["path"])
        for key, path in (
            ("components", component_path),
            ("membership", membership_path),
            ("edges", edge_path),
        ):
            verify_identity(outputs[key], path)

        component_k: dict[str, int] = {}
        component_meta: dict[str, tuple[str, str, bool]] = {}
        component_row_by_id: dict[str, dict[str, str]] = {}
        chrom_spans: list[int] = []
        chrom_ks: list[int] = []
        chrom_counts: Counter[str] = Counter()
        for row in read_tsv_gz(component_path):
            if int(row["threshold"]) != PRIMARY_THRESHOLD:
                continue
            component_id = row["component_id"]
            if not component_id or component_id in component_k:
                raise ValueError(f"{component_path}: duplicate/empty primary component id")
            basis, ps = require_exact_primary_container(
                dataset=dataset,
                chrom=chrom,
                row=row,
                path=component_path,
            )
            k = int(row["k"])
            if k < 1:
                raise ValueError(f"{component_path}: component k must be positive")
            span_bp = int(row["span_bp"])
            max_gap_bp = int(row["max_adjacent_gap_bp"])
            if span_bp < 0 or max_gap_bp < 0 or max_gap_bp > span_bp:
                raise ValueError(f"{component_path}: invalid span/gap geometry")
            if row["tree_eligible"] not in {"true", "false"}:
                raise ValueError(f"{component_path}: invalid tree_eligible value")
            tree_eligible = row["tree_eligible"] == "true"
            if tree_eligible != (k > 1):
                raise ValueError(f"{component_path}: k/tree eligibility mismatch")
            component_k[component_id] = k
            component_meta[component_id] = (basis, ps, tree_eligible)
            component_row_by_id[component_id] = row
            chrom_counts["components_all"] += 1
            if k == 1:
                chrom_counts["singleton_components"] += 1
            else:
                chrom_counts["tree_regions"] += 1
                chrom_counts["w_memberships"] += k
                chrom_counts["w_span_gt_50kb"] += span_bp > 50_000
                chrom_counts["w_adjacent_gap_gt_50kb"] += max_gap_bp > 50_000
                chrom_counts["w_k_gt12"] += k > 12
                chrom_spans.append(span_bp)
                chrom_ks.append(k)
                dataset_spans.append(span_bp)
                dataset_ks.append(k)
                k_counts[k] += 1
                k_band_counts[k_band(k)] += 1
                span_band_counts[span_band(span_bp)] += 1
            hp = "HP1" if basis == "PS_HP1" else "HP2"
            hp_counts[hp]["components_all"] += 1
            hp_counts[hp]["singleton_components"] += k == 1
            hp_counts[hp]["tree_regions"] += k > 1
            hp_counts[hp]["memberships"] += k
            hp_counts[hp]["w_memberships"] += k if k > 1 else 0

        site_position: dict[tuple[str, str, int], int] = {}
        chrom_active_loci: set[int] = set()
        chrom_w_loci: set[int] = set()
        chrom_singleton_loci: set[int] = set()
        membership_count = 0
        component_member_positions: dict[str, list[int]] = defaultdict(list)
        component_member_site_keys: dict[str, set[tuple[str, str, int]]] = defaultdict(set)
        primary_membership_site_keys: set[tuple[str, str, int]] = set()
        site_component: dict[tuple[str, str, int], str] = {}
        for row in read_tsv_gz(membership_path):
            if int(row["threshold"]) != PRIMARY_THRESHOLD:
                continue
            component_id = row["component_id"]
            if component_id not in component_k:
                raise ValueError(f"{membership_path}: unknown component id")
            k = component_k[component_id]
            pos = int(row["pos1"])
            site_index = int(row["site_index"])
            basis, ps = require_exact_primary_container(
                dataset=dataset,
                chrom=chrom,
                row=row,
                path=membership_path,
            )
            expected_basis, expected_ps, expected_tree_eligible = component_meta[component_id]
            if (basis, ps) != (expected_basis, expected_ps):
                raise ValueError(
                    f"{membership_path}: membership/container differs from component"
                )
            if row["tree_eligible"] not in {"true", "false"}:
                raise ValueError(f"{membership_path}: invalid tree_eligible value")
            if (row["tree_eligible"] == "true") != expected_tree_eligible:
                raise ValueError(
                    f"{membership_path}: membership tree eligibility differs from component"
                )
            key = (basis, ps, site_index)
            if key in primary_membership_site_keys:
                raise ValueError(
                    f"{membership_path}: primary site membership occurs more than once"
                )
            if pos < 1:
                raise ValueError(f"{membership_path}: pos1 must be positive")
            primary_membership_site_keys.add(key)
            site_position[key] = pos
            site_component[key] = component_id
            component_member_positions[component_id].append(pos)
            component_member_site_keys[component_id].add(key)
            membership_count += 1
            chrom_active_loci.add(pos)
            active_loci.add((chrom, pos))
            exact_ps.add((chrom, ps))
            exact_ps_hp.add((chrom, basis, ps))
            if k > 1:
                chrom_w_loci.add(pos)
                loci_in_w.add((chrom, pos))
            else:
                chrom_singleton_loci.add(pos)
                loci_in_singleton.add((chrom, pos))

        for component_id, expected_k in component_k.items():
            positions = sorted(component_member_positions[component_id])
            site_keys = component_member_site_keys[component_id]
            if (
                len(positions) != expected_k
                or len(set(positions)) != expected_k
                or len(site_keys) != expected_k
            ):
                raise ValueError(
                    f"{membership_path}: component {component_id} membership count "
                    f"{len(positions)} != declared k={expected_k}"
                )
            component_row = component_row_by_id[component_id]
            observed_span = positions[-1] - positions[0]
            observed_max_gap = max(
                (right - left for left, right in zip(positions, positions[1:])),
                default=0,
            )
            if (
                int(component_row["start1"]) != positions[0]
                or int(component_row["end1"]) != positions[-1]
                or int(component_row["span_bp"]) != observed_span
                or int(component_row["max_adjacent_gap_bp"]) != observed_max_gap
            ):
                raise ValueError(
                    f"{component_path}: component {component_id} geometry mismatch"
                )

        retained_edges = 0
        chrom_direct_edge_distances: list[int] = []
        retained_adjacency: dict[str, dict[tuple[str, str, int], set[tuple[str, str, int]]]] = (
            defaultdict(lambda: defaultdict(set))
        )
        short_edge_adjacency: dict[
            str, dict[tuple[str, str, int], set[tuple[str, str, int]]]
        ] = defaultdict(lambda: defaultdict(set))
        components_with_long_edge: set[str] = set()
        edge_pairs: set[tuple[str, str, int, int]] = set()
        component_edge_supports: dict[str, list[int]] = defaultdict(list)
        for row in read_tsv_gz(edge_path):
            basis, ps = require_exact_primary_container(
                dataset=dataset,
                chrom=chrom,
                row=row,
                path=edge_path,
            )
            support = int(row["support_total"])
            state_support = {
                state: int(row[f"support_{state}"])
                for state in ("RR", "RA", "AR", "AA")
            }
            if min(support, *state_support.values()) < 0:
                raise ValueError(f"{edge_path}: edge support must be nonnegative")
            if support != sum(state_support.values()):
                raise ValueError(f"{edge_path}: edge state support mass mismatch")
            if row.get("hp_family") != VALID_LINKAGE_BASIS[basis]:
                raise ValueError(f"{edge_path}: HP family differs from linkage basis")
            if row.get("passes_primary_threshold") not in {"true", "false"}:
                raise ValueError(f"{edge_path}: invalid primary-threshold flag")
            passes_primary = row["passes_primary_threshold"] == "true"
            if passes_primary != (support >= PRIMARY_THRESHOLD):
                raise ValueError(f"{edge_path}: primary-threshold flag mismatch")
            if support < PRIMARY_THRESHOLD:
                continue
            site_i = int(row["site_i_index"])
            site_j = int(row["site_j_index"])
            pair_key = (basis, ps, site_i, site_j)
            if site_i >= site_j or pair_key in edge_pairs:
                raise ValueError(f"{edge_path}: duplicate or noncanonical retained edge pair")
            edge_pairs.add(pair_key)
            left_key = (basis, ps, site_i)
            right_key = (basis, ps, site_j)
            left = site_position.get(left_key)
            right = site_position.get(right_key)
            if left is None or right is None:
                raise ValueError(f"{edge_path}: retained edge endpoint missing from membership")
            if site_component[left_key] != site_component[right_key]:
                raise ValueError(f"{edge_path}: retained edge crosses graph components")
            declared_left = int(row["pos_i1"])
            declared_right = int(row["pos_j1"])
            declared_gap = int(row["gap_bp"])
            distance = right - left
            if (
                distance < 0
                or (declared_left, declared_right, declared_gap)
                != (left, right, distance)
            ):
                raise ValueError(f"{edge_path}: retained edge coordinate/gap mismatch")
            component_id = site_component[left_key]
            component_edge_supports[component_id].append(support)
            retained_adjacency[component_id][left_key].add(right_key)
            retained_adjacency[component_id][right_key].add(left_key)
            if distance <= 50_000:
                short_edge_adjacency[component_id][left_key].add(right_key)
                short_edge_adjacency[component_id][right_key].add(left_key)
            else:
                components_with_long_edge.add(component_id)
            retained_edges += 1
            chrom_direct_edge_distances.append(distance)
            dataset_direct_edge_distances.append(distance)
            chrom_counts["direct_edges_gt_50kb"] += distance > 50_000
            if distance > 50_000:
                if support == 3:
                    dataset_counts["direct_edges_gt_50kb_support_3"] += 1
                elif support == 4:
                    dataset_counts["direct_edges_gt_50kb_support_4"] += 1
                else:
                    dataset_counts["direct_edges_gt_50kb_support_ge_5"] += 1
                alt_informative = any(
                    state_support[state] > 0 for state in ("RA", "AR", "AA")
                )
                dataset_counts[
                    "direct_edges_gt_50kb_alt_informative"
                    if alt_informative
                    else "direct_edges_gt_50kb_RR_only"
                ] += 1

        for component_id, expected_k in component_k.items():
            declared_component = component_row_by_id[component_id]
            observed_supports = component_edge_supports.get(component_id, [])
            declared_edge_count = int(declared_component["retained_edge_count"])
            declared_min_support = declared_component[
                "min_retained_endpoint_support"
            ]
            declared_max_support = declared_component[
                "max_retained_endpoint_support"
            ]
            if (
                declared_edge_count != len(observed_supports)
                or (
                    observed_supports
                    and (
                        int(declared_min_support) != min(observed_supports)
                        or int(declared_max_support) != max(observed_supports)
                    )
                )
                or (
                    not observed_supports
                    and (declared_min_support != "" or declared_max_support != "")
                )
            ):
                raise ValueError(
                    f"{component_path}: component {component_id} retained-edge "
                    "count/support summary mismatch"
                )
            if expected_k == 1:
                if retained_adjacency.get(component_id):
                    raise ValueError(
                        f"{edge_path}: singleton component has a retained edge"
                    )
                continue
            nodes = component_member_site_keys[component_id]
            adjacency = retained_adjacency[component_id]
            if sum(len(neighbours) for neighbours in adjacency.values()) // 2 < expected_k - 1:
                raise ValueError(
                    f"{edge_path}: component {component_id} has too few retained edges"
                )
            start = next(iter(nodes))
            visited = {start}
            stack = [start]
            while stack:
                current = stack.pop()
                for neighbour in adjacency.get(current, set()):
                    if neighbour not in visited:
                        visited.add(neighbour)
                        stack.append(neighbour)
            if visited != nodes:
                raise ValueError(
                    f"{edge_path}: component {component_id} is not connected by retained edges"
                )
            remaining = set(nodes)
            part_sizes: list[int] = []
            while remaining:
                seed = next(iter(remaining))
                part = {seed}
                stack = [seed]
                remaining.remove(seed)
                while stack:
                    current = stack.pop()
                    for neighbour in short_edge_adjacency[component_id].get(
                        current, set()
                    ):
                        if neighbour in remaining:
                            remaining.remove(neighbour)
                            part.add(neighbour)
                            stack.append(neighbour)
                part_sizes.append(len(part))
            counterfactual_w = sum(size >= 2 for size in part_sizes)
            singleton_fragments = sum(size == 1 for size in part_sizes)
            chrom_counts["counterfactual_W_after_50kb_cut"] += counterfactual_w
            chrom_counts[
                "W_memberships_lost_to_singletons_if_50kb_cut"
            ] += singleton_fragments
            if component_id in components_with_long_edge:
                chrom_counts["W_with_direct_edge_gt_50kb"] += 1
            if len(part_sizes) > 1:
                chrom_counts["W_partition_changed_if_50kb_cut"] += 1
                if component_id not in components_with_long_edge:
                    raise ValueError(
                        f"{edge_path}: 50 kb counterfactual changed a W with no long edge"
                    )
                if counterfactual_w == 0:
                    chrom_counts["W_fully_lost_if_50kb_cut"] += 1
                elif counterfactual_w == 1:
                    chrom_counts["W_reduced_to_one_W_if_50kb_cut"] += 1
                else:
                    chrom_counts["W_split_to_multiple_W_if_50kb_cut"] += 1

        chrom_counts["counterfactual_W_delta_after_50kb_cut"] = (
            chrom_counts["counterfactual_W_after_50kb_cut"]
            - chrom_counts["tree_regions"]
        )
        if (
            chrom_counts["W_partition_changed_if_50kb_cut"]
            > chrom_counts["W_with_direct_edge_gt_50kb"]
            or chrom_counts["W_with_direct_edge_gt_50kb"]
            > chrom_counts["tree_regions"]
            or chrom_counts["W_partition_changed_if_50kb_cut"]
            != chrom_counts["W_fully_lost_if_50kb_cut"]
            + chrom_counts["W_reduced_to_one_W_if_50kb_cut"]
            + chrom_counts["W_split_to_multiple_W_if_50kb_cut"]
        ):
            raise ValueError(f"{chrom}: invalid 50 kb counterfactual partition")

        candidate_s = int(receipt["scope"]["candidate_loci_S"])
        chrom_counts["active_memberships"] = membership_count
        chrom_counts["active_unique_loci"] = len(chrom_active_loci)
        chrom_counts["unique_loci_in_w"] = len(chrom_w_loci)
        chrom_counts["retained_edges"] = retained_edges
        if chrom_counts["components_all"] != (
            chrom_counts["singleton_components"] + chrom_counts["tree_regions"]
        ):
            raise ValueError(f"{chrom}: component partition is not conserved")
        if membership_count != (
            chrom_counts["singleton_components"] + chrom_counts["w_memberships"]
        ):
            raise ValueError(f"{chrom}: membership mass is not conserved")

        chromosome_rows.append(
            {
                "dataset": dataset,
                "chrom": chrom,
                "candidate_loci_S": candidate_s,
                "active_unique_loci": len(chrom_active_loci),
                "unique_loci_in_any_W": len(chrom_w_loci),
                "active_node_memberships": membership_count,
                "all_components": chrom_counts["components_all"],
                "singleton_k1_components": chrom_counts["singleton_components"],
                "tree_eligible_W": chrom_counts["tree_regions"],
                "W_memberships": chrom_counts["w_memberships"],
                "retained_direct_edges": retained_edges,
                "direct_edges_gt_50kb": chrom_counts["direct_edges_gt_50kb"],
                "W_span_gt_50kb": chrom_counts["w_span_gt_50kb"],
                "W_adjacent_gap_gt_50kb": chrom_counts["w_adjacent_gap_gt_50kb"],
                "W_with_direct_edge_gt_50kb": chrom_counts[
                    "W_with_direct_edge_gt_50kb"
                ],
                "W_partition_changed_if_50kb_cut": chrom_counts[
                    "W_partition_changed_if_50kb_cut"
                ],
                "W_fully_lost_if_50kb_cut": chrom_counts[
                    "W_fully_lost_if_50kb_cut"
                ],
                "W_reduced_to_one_W_if_50kb_cut": chrom_counts[
                    "W_reduced_to_one_W_if_50kb_cut"
                ],
                "W_split_to_multiple_W_if_50kb_cut": chrom_counts[
                    "W_split_to_multiple_W_if_50kb_cut"
                ],
                "counterfactual_W_after_50kb_cut": chrom_counts[
                    "counterfactual_W_after_50kb_cut"
                ],
                "counterfactual_W_delta_after_50kb_cut": chrom_counts[
                    "counterfactual_W_delta_after_50kb_cut"
                ],
                "W_memberships_lost_to_singletons_if_50kb_cut": chrom_counts[
                    "W_memberships_lost_to_singletons_if_50kb_cut"
                ],
                "W_k_gt12": chrom_counts["w_k_gt12"],
                "median_W_k": statistics.median(chrom_ks) if chrom_ks else None,
                "median_W_span_bp": statistics.median(chrom_spans) if chrom_spans else None,
                "max_W_span_bp": max(chrom_spans, default=None),
                "max_direct_edge_distance_bp": max(chrom_direct_edge_distances, default=None),
                "all_pass": True,
            }
        )
        dataset_counts.update(chrom_counts)

    aggregate = summary["aggregate"]
    summary_thresholds = summary.get("threshold_sensitivity")
    if (
        not isinstance(summary_thresholds, Mapping)
        or set(summary_thresholds) != {"1", "2", "3", "5"}
    ):
        raise ValueError(f"{summary_path}: threshold sensitivity grid mismatch")
    for threshold in (1, 2, 3, 5):
        summary_threshold = summary_thresholds[str(threshold)]
        recomputed_threshold = {
            field: receipt_threshold_totals[threshold][field]
            for field in THRESHOLD_SUMMARY_FIELDS
        }
        if (
            not isinstance(summary_threshold, Mapping)
            or {
                field: summary_threshold.get(field)
                for field in THRESHOLD_SUMMARY_FIELDS
            }
            != recomputed_threshold
        ):
            raise ValueError(
                f"{summary_path}: threshold={threshold} differs from current "
                "chromosome receipt reaggregation"
            )
    expected = {
        "candidate_loci_S": sum(row["candidate_loci_S"] for row in chromosome_rows),
        "active_unique_loci": len(active_loci),
        "active_node_memberships": dataset_counts["active_memberships"],
        "components_all_including_singletons": dataset_counts["components_all"],
        "singleton_unlinked_components_excluded_from_tree": dataset_counts[
            "singleton_components"
        ],
        "tree_eligible_read_linked_regions": dataset_counts["tree_regions"],
        "retained_endpoint_edges": dataset_counts["retained_edges"],
        "k_gt12_regions": dataset_counts["w_k_gt12"],
        "regions_spanning_adjacent_gap_gt_50kb": dataset_counts[
            "w_adjacent_gap_gt_50kb"
        ],
    }
    for key, value in expected.items():
        if aggregate.get(key) != value:
            raise ValueError(
                f"{summary_path}: independent {key}={value} != summary {aggregate.get(key)}"
            )
    if len(exact_ps_hp) != aggregate["exact_ps_hp_containers_with_active_nodes"]:
        raise ValueError(f"{dataset}: independent exact PS×HP container count mismatch")

    tree_only = len(loci_in_w - loci_in_singleton)
    singleton_only = len(loci_in_singleton - loci_in_w)
    both = len(loci_in_w & loci_in_singleton)
    all_components = dataset_counts["components_all"]
    tree_regions = dataset_counts["tree_regions"]
    active_memberships = dataset_counts["active_memberships"]
    retained_edges = dataset_counts["retained_edges"]
    dataset_row = {
        "dataset": dataset,
        "candidate_loci_S": aggregate["candidate_loci_S"],
        "active_unique_loci": len(active_loci),
        "active_unique_pct_of_S": pct(len(active_loci), aggregate["candidate_loci_S"]),
        "unique_loci_in_any_W": len(loci_in_w),
        "unique_loci_in_W_pct_of_S": pct(len(loci_in_w), aggregate["candidate_loci_S"]),
        "unique_loci_in_W_pct_of_active": pct(len(loci_in_w), len(active_loci)),
        "tree_only_unique_loci": tree_only,
        "singleton_only_unique_loci": singleton_only,
        "both_tree_and_singleton_unique_loci": both,
        "exact_chromosome_PS": len(exact_ps),
        "exact_PS_HP_containers": len(exact_ps_hp),
        "active_node_memberships": active_memberships,
        "singleton_memberships": dataset_counts["singleton_components"],
        "W_memberships": dataset_counts["w_memberships"],
        "W_membership_pct_of_active_memberships": pct(
            dataset_counts["w_memberships"], active_memberships
        ),
        "all_components": all_components,
        "singleton_k1_components": dataset_counts["singleton_components"],
        "singleton_pct_of_all_components": pct(
            dataset_counts["singleton_components"], all_components
        ),
        "tree_eligible_W": tree_regions,
        "W_pct_of_all_components": pct(tree_regions, all_components),
        "HP1_W": hp_counts["HP1"]["tree_regions"],
        "HP2_W": hp_counts["HP2"]["tree_regions"],
        "mean_W_k": round(statistics.fmean(dataset_ks), 4) if dataset_ks else None,
        "median_W_k": statistics.median(dataset_ks) if dataset_ks else None,
        "p90_W_k": percentile(dataset_ks, 0.90),
        "p95_W_k": percentile(dataset_ks, 0.95),
        "p99_W_k": percentile(dataset_ks, 0.99),
        "max_W_k": max(dataset_ks, default=None),
        "W_k_gt12": dataset_counts["w_k_gt12"],
        "W_k_gt12_pct": pct(dataset_counts["w_k_gt12"], tree_regions),
        "median_W_span_bp": statistics.median(dataset_spans) if dataset_spans else None,
        "p90_W_span_bp": percentile(dataset_spans, 0.90),
        "p95_W_span_bp": percentile(dataset_spans, 0.95),
        "p99_W_span_bp": percentile(dataset_spans, 0.99),
        "max_W_span_bp": max(dataset_spans, default=None),
        "W_span_gt_50kb": dataset_counts["w_span_gt_50kb"],
        "W_span_gt_50kb_pct": pct(dataset_counts["w_span_gt_50kb"], tree_regions),
        "W_adjacent_gap_gt_50kb": dataset_counts["w_adjacent_gap_gt_50kb"],
        "W_adjacent_gap_gt_50kb_pct": pct(
            dataset_counts["w_adjacent_gap_gt_50kb"], tree_regions
        ),
        "W_span_gt_50kb_without_adjacent_gap_gt_50kb": (
            dataset_counts["w_span_gt_50kb"]
            - dataset_counts["w_adjacent_gap_gt_50kb"]
        ),
        "W_with_direct_edge_gt_50kb": dataset_counts[
            "W_with_direct_edge_gt_50kb"
        ],
        "W_with_direct_edge_gt_50kb_pct": pct(
            dataset_counts["W_with_direct_edge_gt_50kb"], tree_regions
        ),
        "W_partition_changed_if_50kb_cut": dataset_counts[
            "W_partition_changed_if_50kb_cut"
        ],
        "W_partition_changed_if_50kb_cut_pct": pct(
            dataset_counts["W_partition_changed_if_50kb_cut"], tree_regions
        ),
        "W_fully_lost_if_50kb_cut": dataset_counts[
            "W_fully_lost_if_50kb_cut"
        ],
        "W_reduced_to_one_W_if_50kb_cut": dataset_counts[
            "W_reduced_to_one_W_if_50kb_cut"
        ],
        "W_split_to_multiple_W_if_50kb_cut": dataset_counts[
            "W_split_to_multiple_W_if_50kb_cut"
        ],
        "counterfactual_W_after_50kb_cut": dataset_counts[
            "counterfactual_W_after_50kb_cut"
        ],
        "counterfactual_W_delta_after_50kb_cut": dataset_counts[
            "counterfactual_W_delta_after_50kb_cut"
        ],
        "W_memberships_lost_to_singletons_if_50kb_cut": dataset_counts[
            "W_memberships_lost_to_singletons_if_50kb_cut"
        ],
        "W_memberships_lost_to_singletons_if_50kb_cut_pct": pct(
            dataset_counts["W_memberships_lost_to_singletons_if_50kb_cut"],
            dataset_counts["w_memberships"],
        ),
        "retained_direct_edges": retained_edges,
        "direct_edges_gt_50kb": dataset_counts["direct_edges_gt_50kb"],
        "direct_edges_gt_50kb_pct": pct(
            dataset_counts["direct_edges_gt_50kb"], retained_edges
        ),
        "direct_edges_gt_50kb_support_3": dataset_counts[
            "direct_edges_gt_50kb_support_3"
        ],
        "direct_edges_gt_50kb_support_4": dataset_counts[
            "direct_edges_gt_50kb_support_4"
        ],
        "direct_edges_gt_50kb_support_ge_5": dataset_counts[
            "direct_edges_gt_50kb_support_ge_5"
        ],
        "direct_edges_gt_50kb_RR_only": dataset_counts[
            "direct_edges_gt_50kb_RR_only"
        ],
        "direct_edges_gt_50kb_alt_informative": dataset_counts[
            "direct_edges_gt_50kb_alt_informative"
        ],
        "median_direct_edge_distance_bp": statistics.median(
            dataset_direct_edge_distances
        )
        if dataset_direct_edge_distances
        else None,
        "p99_direct_edge_distance_bp": percentile(dataset_direct_edge_distances, 0.99),
        "max_direct_edge_distance_bp": max(dataset_direct_edge_distances, default=None),
        "canonical_molecule_rows": aggregate["canonical_molecule_rows"],
        "primary_known_PS_molecule_rows": aggregate["primary_known_ps_molecule_rows"],
        "excluded_missing_PS_rows": aggregate["excluded_missing_ps_rows"],
        "excluded_nonprimary_HP_rows": aggregate["excluded_nonprimary_hp_rows"],
        "all_pass": True,
    }

    k_band_rows = [
        {
            "dataset": dataset,
            "k_band": band,
            "W": k_band_counts[band],
            "pct_of_W": pct(k_band_counts[band], tree_regions),
        }
        for band in ("2", "3", "4", "5", "6–8", "9–12", ">12")
    ]
    span_band_rows = [
        {
            "dataset": dataset,
            "span_band": band,
            "W": span_band_counts[band],
            "pct_of_W": pct(span_band_counts[band], tree_regions),
        }
        for band in ("≤10 kb", "10–25 kb", "25–50 kb", "50–100 kb", ">100 kb")
    ]
    k_rows = [
        {
            "dataset": dataset,
            "k": k,
            "W": count,
            "pct_of_W": pct(count, tree_regions),
        }
        for k, count in sorted(k_counts.items())
    ]
    threshold_rows = [
        {
            "dataset": dataset,
            "minimum_distinct_molecules": int(threshold),
            "all_components": values["regions"],
            "singleton_k1_components": values["k1_regions"],
            "multisite_components": values["k_gt1_regions"],
            "retained_direct_edges": values["retained_endpoint_edges"],
            "validation_basis": (
                "independent_row_recomputation_and_receipt_reconciliation"
                if int(threshold) == PRIMARY_THRESHOLD
                else "chromosome_receipt_reaggregation_and_summary_reconciliation"
            ),
        }
        for threshold, values in sorted(receipt_threshold_totals.items())
    ]
    return {
        "dataset_row": dataset_row,
        "chromosome_rows": chromosome_rows,
        "k_band_rows": k_band_rows,
        "span_band_rows": span_band_rows,
        "k_rows": k_rows,
        "threshold_rows": threshold_rows,
        "summary_identity": identity(summary_path),
        "input_receipts": input_receipts,
        "extraction_receipts": extraction_receipts,
        "summary_checks": summary["checks"],
        "primary_threshold": summary["scope"]["primary_threshold"],
        "claim_ceiling": summary["claim_ceiling"],
    }


def build(
    dataset_roots: Mapping[str, Path],
    output_dir: Path,
    topology_audit_path: Path,
) -> dict[str, Any]:
    if set(dataset_roots) != set(CANONICAL_DATASETS):
        raise ValueError(
            f"full canonical dataset set required: {sorted(dataset_roots)} "
            f"!= {sorted(CANONICAL_DATASETS)}"
        )
    output_dir = make_output_dir(output_dir)
    loaded = [
        load_dataset(dataset, dataset_roots[dataset]) for dataset in CANONICAL_DATASETS
    ]
    topology = load_topology_completion_audit(
        topology_audit_path,
        loaded_datasets=loaded,
    )
    topology_rows = topology["audit"]["rows"]
    dataset_rows = [item["dataset_row"] for item in loaded]
    chromosome_rows = [
        row for item in loaded for row in item["chromosome_rows"]
    ]
    k_band_rows = [row for item in loaded for row in item["k_band_rows"]]
    span_band_rows = [row for item in loaded for row in item["span_band_rows"]]
    k_rows = [row for item in loaded for row in item["k_rows"]]
    threshold_rows = [row for item in loaded for row in item["threshold_rows"]]

    checks = {
        "canonical_7_datasets_present": len(dataset_rows) == 7,
        "all_154_autosome_rows_present": len(chromosome_rows) == 7 * 22,
        "all_154_chromosome_receipt_identities_verified": all(
            len(item["input_receipts"]) == 22 for item in loaded
        ),
        "all_154_extraction_receipt_identities_verified": (
            sum(len(item["extraction_receipts"]) for item in loaded) == 7 * 22
            and all(
                set(item["extraction_receipts"]) == set(AUTOSOMES)
                for item in loaded
            )
        ),
        "row_level_exact_ps_hp_graph_contracts_validated": len(loaded) == 7,
        "all_dataset_summaries_pass": all(row["all_pass"] for row in dataset_rows),
        "all_chromosome_rows_pass": all(row["all_pass"] for row in chromosome_rows),
        "latest_topology_completion_audit_verified": (
            topology["audit"]["summary"][
                "strict_linkage_complete_datasets"
            ]
            == 7
            and topology["audit"]["summary"][
                "strict_directed_topology_production_complete_datasets"
            ]
            == 0
        ),
        "component_partition_conserved": all(
            row["all_components"]
            == row["singleton_k1_components"] + row["tree_eligible_W"]
            for row in dataset_rows
        ),
        "membership_mass_conserved": all(
            row["active_node_memberships"]
            == row["singleton_memberships"] + row["W_memberships"]
            for row in dataset_rows
        ),
        "physical_locus_role_union_conserved": all(
            row["active_unique_loci"]
            == row["tree_only_unique_loci"]
            + row["singleton_only_unique_loci"]
            + row["both_tree_and_singleton_unique_loci"]
            for row in dataset_rows
        ),
        "physical_locus_funnel_ordered": all(
            row["unique_loci_in_any_W"]
            <= row["active_unique_loci"]
            <= row["candidate_loci_S"]
            for row in dataset_rows
        ),
        "canonical_molecule_row_routing_conserved": all(
            row["canonical_molecule_rows"]
            == row["primary_known_PS_molecule_rows"]
            + row["excluded_missing_PS_rows"]
            + row["excluded_nonprimary_HP_rows"]
            for row in dataset_rows
        ),
        "distance_qc_nesting_conserved": all(
            row["W_adjacent_gap_gt_50kb"]
            <= row["W_span_gt_50kb"]
            <= row["tree_eligible_W"]
            and row["W_span_gt_50kb_without_adjacent_gap_gt_50kb"]
            == row["W_span_gt_50kb"] - row["W_adjacent_gap_gt_50kb"]
            for row in dataset_rows
        ),
        "long_edge_support_classes_conserved": all(
            row["direct_edges_gt_50kb"]
            == row["direct_edges_gt_50kb_support_3"]
            + row["direct_edges_gt_50kb_support_4"]
            + row["direct_edges_gt_50kb_support_ge_5"]
            for row in dataset_rows
        ),
        "long_edge_state_classes_conserved": all(
            row["direct_edges_gt_50kb"]
            == row["direct_edges_gt_50kb_RR_only"]
            + row["direct_edges_gt_50kb_alt_informative"]
            for row in dataset_rows
        ),
        "long_edges_are_subset_of_retained_edges": all(
            row["direct_edges_gt_50kb"] <= row["retained_direct_edges"]
            for row in dataset_rows
        ),
        "distance_cut_counterfactual_partition_conserved": all(
            row["W_partition_changed_if_50kb_cut"]
            <= row["W_with_direct_edge_gt_50kb"]
            <= row["tree_eligible_W"]
            and row["W_partition_changed_if_50kb_cut"]
            == row["W_fully_lost_if_50kb_cut"]
            + row["W_reduced_to_one_W_if_50kb_cut"]
            + row["W_split_to_multiple_W_if_50kb_cut"]
            and row["counterfactual_W_after_50kb_cut"]
            == row["tree_eligible_W"]
            + row["counterfactual_W_delta_after_50kb_cut"]
            and 0
            <= row["W_memberships_lost_to_singletons_if_50kb_cut"]
            <= row["W_memberships"]
            for row in dataset_rows
        ),
        "k_bands_sum_to_W": all(
            sum(
                row["W"]
                for row in k_band_rows
                if row["dataset"] == dataset
            )
            == next(
                row["tree_eligible_W"]
                for row in dataset_rows
                if row["dataset"] == dataset
            )
            for dataset in CANONICAL_DATASETS
        ),
        "span_bands_sum_to_W": all(
            sum(
                row["W"]
                for row in span_band_rows
                if row["dataset"] == dataset
            )
            == next(
                row["tree_eligible_W"]
                for row in dataset_rows
                if row["dataset"] == dataset
            )
            for dataset in CANONICAL_DATASETS
        ),
        "exact_k_distribution_sums_to_W": all(
            sum(
                row["W"]
                for row in k_rows
                if row["dataset"] == dataset
            )
            == next(
                row["tree_eligible_W"]
                for row in dataset_rows
                if row["dataset"] == dataset
            )
            for dataset in CANONICAL_DATASETS
        ),
        "chromosome_totals_reconcile_to_dataset": all(
            all(
                sum(
                    chrom_row[chrom_field]
                    for chrom_row in chromosome_rows
                    if chrom_row["dataset"] == dataset_row["dataset"]
                )
                == dataset_row[dataset_field]
                for chrom_field, dataset_field in (
                    ("candidate_loci_S", "candidate_loci_S"),
                    ("active_node_memberships", "active_node_memberships"),
                    ("all_components", "all_components"),
                    ("singleton_k1_components", "singleton_k1_components"),
                    ("tree_eligible_W", "tree_eligible_W"),
                    ("W_memberships", "W_memberships"),
                    ("retained_direct_edges", "retained_direct_edges"),
                    ("direct_edges_gt_50kb", "direct_edges_gt_50kb"),
                    ("W_span_gt_50kb", "W_span_gt_50kb"),
                    ("W_adjacent_gap_gt_50kb", "W_adjacent_gap_gt_50kb"),
                    (
                        "W_with_direct_edge_gt_50kb",
                        "W_with_direct_edge_gt_50kb",
                    ),
                    (
                        "W_partition_changed_if_50kb_cut",
                        "W_partition_changed_if_50kb_cut",
                    ),
                    (
                        "W_fully_lost_if_50kb_cut",
                        "W_fully_lost_if_50kb_cut",
                    ),
                    (
                        "W_reduced_to_one_W_if_50kb_cut",
                        "W_reduced_to_one_W_if_50kb_cut",
                    ),
                    (
                        "W_split_to_multiple_W_if_50kb_cut",
                        "W_split_to_multiple_W_if_50kb_cut",
                    ),
                    (
                        "counterfactual_W_after_50kb_cut",
                        "counterfactual_W_after_50kb_cut",
                    ),
                    (
                        "counterfactual_W_delta_after_50kb_cut",
                        "counterfactual_W_delta_after_50kb_cut",
                    ),
                    (
                        "W_memberships_lost_to_singletons_if_50kb_cut",
                        "W_memberships_lost_to_singletons_if_50kb_cut",
                    ),
                    ("W_k_gt12", "W_k_gt12"),
                )
            )
            for dataset_row in dataset_rows
        ),
        "primary_threshold_is_three": all(
            item["primary_threshold"] == PRIMARY_THRESHOLD == 3 for item in loaded
        ),
        "expected_threshold_sensitivity_grid_present": all(
            {
                row["minimum_distinct_molecules"]
                for row in threshold_rows
                if row["dataset"] == dataset
            }
            == {1, 2, 3, 5}
            for dataset in CANONICAL_DATASETS
        ),
        "retained_edges_monotone_with_threshold": all(
            all(
                left["retained_direct_edges"] >= right["retained_direct_edges"]
                for left, right in zip(rows, rows[1:])
            )
            for dataset in CANONICAL_DATASETS
            for rows in (
                [
                    row
                    for row in threshold_rows
                    if row["dataset"] == dataset
                ],
            )
        ),
    }
    if not all(checks.values()):
        raise ValueError(f"all-7 report data checks failed: {checks}")

    tables = {
        "dataset_summary.tsv": dataset_rows,
        "chromosome_summary.tsv": chromosome_rows,
        "k_band_distribution.tsv": k_band_rows,
        "span_band_distribution.tsv": span_band_rows,
        "exact_k_distribution.tsv": k_rows,
        "threshold_sensitivity.tsv": threshold_rows,
        "topology_completion_status.tsv": topology_rows,
    }
    for name, rows in tables.items():
        write_tsv(output_dir / name, rows)

    document = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": all(checks.values()),
        "scope": {
            "task_type": "B_comprehensive_validation",
            "datasets": list(CANONICAL_DATASETS),
            "chromosomes": list(AUTOSOMES),
            "primary_minimum_distinct_molecules": PRIMARY_THRESHOLD,
        },
        "definitions": {
            "container": "dataset × chromosome × exact nonmissing PS × HP1/HP2",
            "node": "sSNV locus with a fixed R/A observation in the container",
            "edge": "two endpoints fixed-observed by at least 3 distinct canonical molecules",
            "W": "maximal connected component of qualifying edges with k≥2",
            "singleton": "k=1 audit component; excluded from tree/topology inference",
            "span": "max(position) − min(position); a descriptive envelope, not a membership rule",
            "distance_policy": "no hard distance cutoff; distance is retained as QC",
            "threshold_sensitivity_provenance": (
                "threshold=3 is independently recomputed from component, membership, "
                "and edge rows; thresholds 1/2/5 are reaggregated from the current 22 "
                "chromosome receipts and reconciled to the dataset summary"
            ),
        },
        "datasets": dataset_rows,
        "chromosomes": chromosome_rows,
        "k_bands": k_band_rows,
        "span_bands": span_band_rows,
        "exact_k": k_rows,
        "threshold_sensitivity": threshold_rows,
        "topology_completion": topology_rows,
        "topology_completion_summary": topology["audit"]["summary"],
        "checks": checks,
        "inputs": {
            item["dataset_row"]["dataset"]: {
                "summary": item["summary_identity"],
                "chromosome_receipts": item["input_receipts"],
                "extraction_receipts": item["extraction_receipts"],
            }
            for item in loaded
        },
        "supporting_inputs": {
            "topology_completion_audit": topology["audit_identity"],
            "topology_completion_audit_receipt": topology["receipt_identity"],
            "topology_completion_audit_receipt_sidecar": topology[
                "receipt_sidecar_identity"
            ],
        },
        "claim_ceiling": (
            "This package validates exact-PS/HP strict read-linked regions. "
            "The current production v4 strict directed-topology stage is complete "
            "for 0/7 datasets. It does not establish a unique mutation-state tree, "
            "cellular clone count, exact cellular parent-child edge, cross-HP/"
            "technical fused tree, or subclone truth."
        ),
    }
    data_path = output_dir / "all7_report_data.json"
    write_json(data_path, document)
    receipt = {
        "schema_name": "intersubmod.strict_region_all7_report_data_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": document["all_pass"],
        "checks": checks,
        "outputs": {
            path.name: identity(path)
            for path in sorted(output_dir.iterdir())
            if path.is_file()
        },
    }
    receipt_path = output_dir / "receipt.json"
    write_json(receipt_path, receipt)
    sidecar_path = receipt_path.with_name("receipt.json.sha256")
    with sidecar_path.open("x", encoding="ascii") as handle:
        handle.write(f"{sha256_path(receipt_path)}  receipt.json\n")
        handle.flush()
        os.fsync(handle.fileno())
    return document


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-root",
        action="append",
        required=True,
        type=parse_dataset_root,
        help="Repeat exactly seven times as DATASET=/absolute/strict/root",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--topology-audit",
        required=True,
        type=Path,
        help="Verified all-7 L1-L4 topology completion audit JSON",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    dataset_roots = dict(args.dataset_root)
    if len(dataset_roots) != len(args.dataset_root):
        raise ValueError("duplicate --dataset-root dataset")
    result = build(dataset_roots, args.output_dir, args.topology_audit)
    print(
        json.dumps(
            {
                "all_pass": result["all_pass"],
                "datasets": len(result["datasets"]),
                "chromosome_rows": len(result["chromosomes"]),
                "tree_eligible_W_total": sum(
                    row["tree_eligible_W"] for row in result["datasets"]
                ),
                "strict_topology_complete": result[
                    "topology_completion_summary"
                ]["strict_directed_topology_production_complete_datasets"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
