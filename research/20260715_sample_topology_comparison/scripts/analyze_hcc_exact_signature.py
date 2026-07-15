#!/usr/bin/env python3
"""Validate exact-topology concordance for HCC1395 and HCC1395_DORADO.

The current-v5 read-AF sidecar intentionally omits payloads for fixed
single-candidate primary units.  This analysis therefore joins each sidecar to
its SHA-bound layered reconstruction and layered region view.  It compares
primary-unit signature *multisets* rather than HP1/HP2 names so a harmless HP
label flip cannot create a false discordance.

The two reported topology levels are deliberately different:

* ``shape`` is an unlabeled rooted-tree shape and does not identify mutations.
* ``exact_labeled_edges`` preserves mutation-state labels and is only compared
  when the exact region, internal sSNV set, and selected candidate are all
  comparable.

Neither level confirms a biological clone, true ancestry, CCF, or sample
identity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shlex
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
EXPECTED_SAMPLES = ("HCC1395", "HCC1395_DORADO")


class ContractError(RuntimeError):
    """Raised when an input or comparison contract is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def load_json(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot read JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"JSON root is not an object: {path}")
    return payload


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def resolve_existing(path: Path, kind: str) -> Path:
    resolved = path.expanduser().resolve()
    if kind == "file":
        require(resolved.is_file(), f"missing input file: {resolved}")
    elif kind == "directory":
        require(resolved.is_dir(), f"missing input directory: {resolved}")
    else:
        raise AssertionError(f"unsupported path kind: {kind}")
    return resolved


def canonical_shape(edges: Sequence[Sequence[str]]) -> str:
    """Return the builder-compatible unlabeled rooted-forest signature."""

    children: Dict[str, List[str]] = defaultdict(list)
    child_nodes = set()
    all_nodes = set()
    for edge in edges:
        require(len(edge) == 2, f"topology edge must have two nodes: {edge}")
        parent, child = str(edge[0]), str(edge[1])
        children[parent].append(child)
        child_nodes.add(child)
        all_nodes.update((parent, child))
    roots = sorted(all_nodes - child_nodes)

    def visit(node: str, active: Tuple[str, ...]) -> str:
        require(node not in active, f"cycle detected while canonicalizing tree at {node}")
        return "(" + "".join(
            sorted(visit(child, active + (node,)) for child in children[node])
        ) + ")"

    return "|".join(sorted(visit(root, ()) for root in roots)) if roots else "()"


def canonical_edges(edges: Sequence[Sequence[str]]) -> Tuple[Tuple[str, str], ...]:
    """Return a candidate's exact labeled-edge signature."""

    normalized = []
    for edge in edges:
        require(len(edge) == 2, f"topology edge must have two nodes: {edge}")
        normalized.append((str(edge[0]), str(edge[1])))
    return tuple(sorted(normalized))


def region_multiset(signatures: Iterable[Any]) -> Tuple[Any, ...]:
    """Ignore HP family labels while preserving unit multiplicity."""

    return tuple(sorted(signatures))


def unit_complete(unit: Mapping[str, Any]) -> bool:
    """Mirror the current-v5 builder's fail-closed completeness contract."""

    return (
        unit.get("analysis_candidate_set_complete") is True
        and unit.get("capped") is False
        and unit.get("verification_status") == "full_pass"
        and unit.get("verification_complete") is True
        and unit.get("verify_pass") is True
        and int(unit.get("analysis_trees_generated") or -1)
        == int(unit.get("n_trees") or 0)
        and unit.get("n_distinct_shapes_exact") is not None
    )


def derive_unit_signatures(
    unit: Mapping[str, Any], ranking: Optional[Mapping[str, Any]]
) -> Dict[str, Any]:
    """Resolve shape and exact-edge signatures for one primary unit.

    Fixed ``n_trees=1`` units are supplemented from the layered reconstruction
    because current-v5 sidecars intentionally omit their unit payload.  For an
    ambiguous unit, the source tree supplies a shape only when its structural
    shape is already unique; otherwise the SHA-bound read-AF top-shape
    representative is used.  Exact labeled edges require one exact candidate.
    """

    result: Dict[str, Any] = {
        "shape": None,
        "shape_source": None,
        "exact_labeled_edges": None,
        "edge_source": None,
    }
    if not unit_complete(unit):
        return result

    trees = unit.get("trees") or []
    n_trees = int(unit.get("n_trees") or 0)
    n_shapes = int(unit.get("n_distinct_shapes_exact") or 0)

    if n_shapes == 1 and trees:
        result["shape"] = canonical_shape(trees[0].get("edges") or [])
        result["shape_source"] = (
            "layered_reconstruction_fixed_unit"
            if n_trees == 1
            else "layered_reconstruction_structural_shape_unique"
        )
    elif (
        ranking
        and ranking.get("status") == "ranked"
        and int(ranking.get("n_top_shapes_exact") or 0) == 1
    ):
        representatives = ranking.get("top_shape_representatives") or []
        signatures = {str(item.get("shape_signature")) for item in representatives}
        require(len(signatures) == 1, "ranked unique top shape has non-unique representatives")
        result["shape"] = next(iter(signatures))
        result["shape_source"] = "sidecar_read_af_top_shape"

    if n_trees == 1 and trees:
        result["exact_labeled_edges"] = canonical_edges(trees[0].get("edges") or [])
        result["edge_source"] = "layered_reconstruction_fixed_unit"
    elif (
        ranking
        and ranking.get("status") == "ranked"
        and int(ranking.get("n_top_exact") or 0) == 1
    ):
        representatives = ranking.get("top_shape_representatives") or []
        require(
            len(representatives) == 1,
            "ranked unique exact candidate must have one top-shape representative",
        )
        result["exact_labeled_edges"] = canonical_edges(
            representatives[0].get("edges") or []
        )
        result["edge_source"] = "sidecar_read_af_unique_candidate"
    return result


def region_key(region: Mapping[str, Any]) -> str:
    return f"{region['chrom']}:{int(region['start'])}-{int(region['end'])}"


def ratio_payload(numerator: int, denominator: int) -> Dict[str, Any]:
    require(denominator > 0, "ratio denominator must be positive")
    require(0 <= numerator <= denominator, "ratio numerator is outside denominator")
    return {
        "numerator": numerator,
        "denominator": denominator,
        "rate": numerator / denominator,
    }


def set_overlap(left: set, right: set) -> Dict[str, Any]:
    intersection = len(left & right)
    union = len(left | right)
    require(union > 0, "set overlap union is empty")
    return {
        "left": len(left),
        "right": len(right),
        "intersection": intersection,
        "union": union,
        "jaccard": intersection / union,
        "left_coverage": intersection / len(left),
        "right_coverage": intersection / len(right),
    }


def path_is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def validate_shared_provenance(
    *,
    sidecars: Mapping[str, Dict[str, Any]],
    sidecar_paths: Mapping[str, Path],
    sidecar_index: Dict[str, Any],
    sidecar_index_path: Path,
    current_summary: Dict[str, Any],
    current_summary_path: Path,
    run_root: Path,
    manifest: Dict[str, Any],
    manifest_path: Path,
) -> Dict[str, Any]:
    require(
        sidecar_index.get("schema_name")
        == "intersubmod.current_v5_read_af_topology_index",
        "unexpected sidecar index schema",
    )
    require(sidecar_index.get("all_checks_pass") is True, "sidecar index checks failed")
    require(
        current_summary.get("schema_name") == "intersubmod.current_layered_topology_summary",
        "unexpected current summary schema",
    )
    require(current_summary.get("all_pass") is True, "current summary all_pass is not true")
    canonical = current_summary.get("canonical") or {}
    require(
        canonical.get("label") == "longphase_s_recalibrated_FILTER_PASS",
        "current summary canonical label is not FILTER_PASS",
    )
    require(Path(canonical.get("run_root", "")).resolve() == run_root, "run-root mismatch")

    current_sha = sha256_file(current_summary_path)
    manifest_sha = sha256_file(manifest_path)
    index_records = {row["sample"]: row for row in sidecar_index.get("samples", [])}
    current_records = {
        row["sample"]: row for row in (canonical.get("samples") or [])
    }
    require(set(EXPECTED_SAMPLES).issubset(index_records), "HCC pair absent from sidecar index")
    require(
        set(EXPECTED_SAMPLES).issubset(current_records),
        "HCC pair absent from current summary",
    )
    index_provenance = sidecar_index.get("provenance") or {}
    require(
        Path(index_provenance.get("run_root", "")).resolve() == run_root,
        "sidecar index run-root mismatch",
    )
    require(
        Path(index_provenance.get("current_summary", "")).resolve() == current_summary_path,
        "sidecar index current-summary path mismatch",
    )
    require(
        index_provenance.get("current_summary_sha256") == current_sha,
        "sidecar index current-summary SHA mismatch",
    )
    require(
        Path(index_provenance.get("input_manifest", "")).resolve() == manifest_path,
        "sidecar index manifest path mismatch",
    )
    require(
        index_provenance.get("input_manifest_sha256") == manifest_sha,
        "sidecar index manifest SHA mismatch",
    )

    for sample in EXPECTED_SAMPLES:
        sidecar = sidecars[sample]
        path = sidecar_paths[sample]
        require(
            sidecar.get("schema_name") == "intersubmod.current_v5_read_af_topology_sample",
            f"{sample}: unexpected sidecar schema",
        )
        require(sidecar.get("sample") == sample, f"{sample}: sidecar sample mismatch")
        require(
            (sidecar.get("summary") or {}).get("all_checks_pass") is True,
            f"{sample}: sidecar checks failed",
        )
        provenance = sidecar.get("provenance") or {}
        require(
            Path(provenance.get("run_root", "")).resolve() == run_root,
            f"{sample}: sidecar run-root mismatch",
        )
        require(
            Path(provenance.get("current_summary", "")).resolve() == current_summary_path,
            f"{sample}: current-summary path mismatch",
        )
        require(
            provenance.get("current_summary_sha256") == current_sha,
            f"{sample}: current-summary SHA mismatch",
        )
        require(
            Path(provenance.get("input_manifest", "")).resolve() == manifest_path,
            f"{sample}: manifest path mismatch",
        )
        require(
            provenance.get("input_manifest_sha256") == manifest_sha,
            f"{sample}: manifest SHA mismatch",
        )
        record = index_records[sample]
        require(Path(record.get("output", "")).resolve() == path, f"{sample}: index path mismatch")
        require(
            record.get("output_sha256") == sha256_file(path),
            f"{sample}: sidecar SHA mismatch",
        )
        require(record.get("all_checks_pass") is True, f"{sample}: index sample checks failed")
        current_record = current_records[sample]
        require(current_record.get("pass") is True, f"{sample}: current summary sample failed")
        require(
            int(record.get("W_tree")) == int(sidecar["summary"]["W_tree"])
            == int(current_record.get("W_tree")),
            f"{sample}: W_tree differs across index, sidecar, and current summary",
        )
        require(
            int(record.get("W_primary")) == int(sidecar["summary"]["W_primary"])
            == int(current_record.get("W_primary")),
            f"{sample}: W_primary differs across index, sidecar, and current summary",
        )
        current_paths = current_record.get("paths") or {}
        current_hashes = current_record.get("sha256") or {}
        require(
            Path(current_paths.get("layered_reconstruction", "")).resolve()
            == Path(provenance.get("layered_reconstruction", "")).resolve(),
            f"{sample}: layered reconstruction path differs from current summary",
        )
        require(
            current_hashes.get("layered_reconstruction")
            == provenance.get("layered_reconstruction_sha256"),
            f"{sample}: layered reconstruction SHA differs from current summary",
        )
        require(
            Path(current_paths.get("layered_region_view", "")).resolve()
            == Path(provenance.get("layered_region_view", "")).resolve(),
            f"{sample}: region-view path differs from current summary",
        )
        require(
            current_hashes.get("layered_region_view")
            == provenance.get("layered_region_view_sha256"),
            f"{sample}: region-view SHA differs from current summary",
        )

    manifest_records = {row["sample"]: row for row in manifest.get("samples", [])}
    require(set(EXPECTED_SAMPLES).issubset(manifest_records), "HCC pair absent from manifest")
    require(
        manifest_records["HCC1395"].get("biological_id") == "HCC1395"
        and manifest_records["HCC1395_DORADO"].get("biological_id") == "HCC1395",
        "HCC pair is not bound to one biological_id",
    )
    require(
        manifest_records["HCC1395"].get("replicate_role") == "canonical"
        and manifest_records["HCC1395_DORADO"].get("replicate_role")
        == "platform_replica",
        "unexpected HCC replicate roles",
    )
    require(
        manifest_records["HCC1395"].get("alignment_payload", {}).get("storage_identity_v1", {}).get("identity_sha256")
        != manifest_records["HCC1395_DORADO"].get("alignment_payload", {}).get("storage_identity_v1", {}).get("identity_sha256"),
        "HCC pair unexpectedly uses the same alignment payload identity",
    )
    return {
        "current_summary_sha256": current_sha,
        "input_manifest_sha256": manifest_sha,
        "sidecar_index_sha256": sha256_file(sidecar_index_path),
        "manifest_records": manifest_records,
    }


def load_sample_contract(
    *, sample: str, sidecar: Dict[str, Any], sidecar_path: Path, run_root: Path
) -> Dict[str, Any]:
    provenance = sidecar.get("provenance") or {}
    reconstruction_path = resolve_existing(Path(provenance["layered_reconstruction"]), "file")
    region_view_path = resolve_existing(Path(provenance["layered_region_view"]), "file")
    expected_sample_root = (run_root / "samples" / sample).resolve()
    require(
        path_is_within(reconstruction_path, expected_sample_root),
        f"{sample}: layered reconstruction is outside sample run-root",
    )
    require(
        path_is_within(region_view_path, expected_sample_root),
        f"{sample}: region view is outside sample run-root",
    )
    reconstruction_sha = sha256_file(reconstruction_path)
    region_view_sha = sha256_file(region_view_path)
    require(
        reconstruction_sha == provenance.get("layered_reconstruction_sha256"),
        f"{sample}: layered reconstruction SHA mismatch",
    )
    require(
        region_view_sha == provenance.get("layered_region_view_sha256"),
        f"{sample}: region view SHA mismatch",
    )

    reconstruction = load_json(reconstruction_path)
    region_view = load_json(region_view_path)
    require(reconstruction.get("sample") == sample, f"{sample}: reconstruction sample mismatch")
    require(region_view.get("sample") == sample, f"{sample}: region-view sample mismatch")

    sidecar_regions_list = sidecar.get("regions") or []
    require(
        len(sidecar_regions_list) == int(sidecar["summary"]["W_tree"]),
        f"{sample}: sidecar region count differs from W_tree",
    )
    sidecar_regions: Dict[str, Dict[str, Any]] = {}
    malformed_region_keys = 0
    for region in sidecar_regions_list:
        key = str(region.get("region"))
        require(key not in sidecar_regions, f"{sample}: duplicate sidecar region {key}")
        malformed_region_keys += key != region_key(region)
        require(region.get("chrom") in AUTOSOMES, f"{sample}: non-autosomal region {key}")
        sidecar_regions[key] = region
    require(malformed_region_keys == 0, f"{sample}: malformed sidecar region keys")

    view_regions: Dict[str, Dict[str, Any]] = {}
    internal_sites: Dict[str, frozenset] = {}
    for region in region_view.get("regions") or []:
        key = str(region.get("region"))
        require(key == region_key(region), f"{sample}: malformed region-view key {key}")
        require(key not in view_regions, f"{sample}: duplicate region-view region {key}")
        positions = {
            int(position)
            for lineage in region.get("lineages") or []
            for position in (lineage.get("obs_col_coverage") or {}).keys()
        }
        require(positions, f"{sample}: empty internal sSNV set at {key}")
        require(
            min(positions) == int(region["start"]) and max(positions) == int(region["end"]),
            f"{sample}: region endpoints do not equal internal sSNV extrema at {key}",
        )
        view_regions[key] = region
        internal_sites[key] = frozenset(positions)
    require(
        set(view_regions) == set(sidecar_regions),
        f"{sample}: sidecar and region-view universes differ",
    )

    primary_units: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
    for unit in reconstruction.get("detail") or []:
        if unit.get("is_primary_lineage") is not True:
            continue
        key = str(unit.get("region"))
        family = str(unit.get("family"))
        require(key in sidecar_regions, f"{sample}: primary unit region absent from sidecar: {key}")
        require(
            family not in primary_units[key],
            f"{sample}: duplicate primary unit {key} HP{family}",
        )
        primary_units[key][family] = unit

    fixed_supplemented = 0
    fixed_payload_omission_violations = 0
    internal_ssnv_unit_mismatches = 0
    stored_shape_signatures_checked = 0
    stored_shape_signature_mismatches = 0
    shape_source_counts: Counter = Counter()
    edge_source_counts: Counter = Counter()
    resolved: Dict[str, Dict[str, Any]] = {}
    for key, region in sidecar_regions.items():
        source_units = primary_units.get(key, {})
        expected_families = sorted(str(value) for value in region.get("primary_families") or [])
        require(
            expected_families == sorted(source_units),
            f"{sample}: primary-family join mismatch at {key}",
        )
        sidecar_units = region.get("units") or {}
        shape_signatures = []
        edge_signatures = []
        shape_resolved = bool(expected_families)
        edge_resolved = bool(expected_families)
        for family in expected_families:
            unit = source_units[family]
            ranking = sidecar_units.get(family)
            n_trees = int(unit.get("n_trees") or 0)
            if n_trees == 1:
                fixed_payload_omission_violations += ranking is not None
                if unit_complete(unit):
                    fixed_supplemented += 1
            else:
                require(ranking is not None, f"{sample}: missing ambiguous payload at {key} HP{family}")
            internal_ssnv_unit_mismatches += int(unit.get("n_sSNV") or -1) != len(
                internal_sites[key]
            )
            if ranking:
                for representative in ranking.get("top_shape_representatives") or []:
                    stored_shape_signatures_checked += 1
                    stored_shape_signature_mismatches += (
                        canonical_shape(representative.get("edges") or [])
                        != representative.get("shape_signature")
                    )
            signatures = derive_unit_signatures(unit, ranking)
            if signatures["shape"] is None:
                shape_resolved = False
            else:
                shape_signatures.append(signatures["shape"])
                shape_source_counts[signatures["shape_source"]] += 1
            if signatures["exact_labeled_edges"] is None:
                edge_resolved = False
            else:
                edge_signatures.append(signatures["exact_labeled_edges"])
                edge_source_counts[signatures["edge_source"]] += 1
        resolved[key] = {
            "chrom": region["chrom"],
            "internal_sites": internal_sites[key],
            "primary_unit_count": len(expected_families),
            "shape_multiset": region_multiset(shape_signatures) if shape_resolved else None,
            "exact_labeled_edge_multiset": (
                region_multiset(edge_signatures) if edge_resolved else None
            ),
        }

    summary = sidecar["summary"]
    expected_shape_resolved = int(summary["W_primary"]) - int(
        summary["morphology_classes"]["unresolved"]
    )
    expected_edge_resolved = int(summary["selection_classes"]["structural_exact_unique"]) + int(
        summary["selection_classes"]["read_af_unique_first"]
    )
    actual_shape_resolved = sum(
        item["shape_multiset"] is not None for item in resolved.values()
    )
    actual_edge_resolved = sum(
        item["exact_labeled_edge_multiset"] is not None for item in resolved.values()
    )
    checks = {
        "region_key_unique_and_well_formed": malformed_region_keys == 0,
        "sidecar_region_count_matches_W_tree": len(sidecar_regions) == int(summary["W_tree"]),
        "sidecar_and_region_view_universe_equal": set(sidecar_regions) == set(view_regions),
        "primary_family_join_exact": all(
            sorted(str(value) for value in sidecar_regions[key].get("primary_families") or [])
            == sorted(primary_units.get(key, {}))
            for key in sidecar_regions
        ),
        "internal_ssnv_matches_unit_n_sSNV": internal_ssnv_unit_mismatches == 0,
        "fixed_unit_payload_is_omitted": fixed_payload_omission_violations == 0,
        "fixed_unit_supplement_count_matches_summary": fixed_supplemented
        == int(summary["unit_status"]["fixed_exact_unique"]),
        "stored_shape_signatures_recompute": stored_shape_signature_mismatches == 0,
        "shape_resolved_count_matches_morphology_partition": actual_shape_resolved
        == expected_shape_resolved,
        "edge_resolved_count_matches_selection_partition": actual_edge_resolved
        == expected_edge_resolved,
    }
    require(all(checks.values()), f"{sample}: sample contract checks failed: {checks}")

    site_universe = {
        (sidecar_regions[key]["chrom"], position)
        for key, positions in internal_sites.items()
        for position in positions
    }
    return {
        "sample": sample,
        "sidecar_path": sidecar_path,
        "sidecar_sha256": sha256_file(sidecar_path),
        "reconstruction_path": reconstruction_path,
        "reconstruction_sha256": reconstruction_sha,
        "region_view_path": region_view_path,
        "region_view_sha256": region_view_sha,
        "regions": resolved,
        "site_universe": site_universe,
        "counts": {
            "W_tree": int(summary["W_tree"]),
            "W_primary": int(summary["W_primary"]),
            "source_primary_units": sum(len(units) for units in primary_units.values()),
            "fixed_units_supplemented_from_reconstruction": fixed_supplemented,
            "shape_resolved_regions": actual_shape_resolved,
            "exact_candidate_resolved_regions": actual_edge_resolved,
            "stored_shape_signatures_checked": stored_shape_signatures_checked,
            "stored_shape_signature_mismatches": stored_shape_signature_mismatches,
            "internal_ssnv_unit_mismatches": internal_ssnv_unit_mismatches,
        },
        "signature_source_counts": {
            "shape": dict(sorted(shape_source_counts.items())),
            "exact_labeled_edges": dict(sorted(edge_source_counts.items())),
        },
        "checks": checks,
    }


def signature_metric(
    keys: Iterable[str],
    left: Mapping[str, Dict[str, Any]],
    right: Mapping[str, Dict[str, Any]],
    signature_key: str,
    *,
    require_same_sites: bool,
    require_same_primary_unit_count: bool,
) -> Dict[str, Any]:
    eligible = []
    for key in keys:
        left_row, right_row = left[key], right[key]
        if left_row[signature_key] is None or right_row[signature_key] is None:
            continue
        if require_same_sites and left_row["internal_sites"] != right_row["internal_sites"]:
            continue
        if (
            require_same_primary_unit_count
            and left_row["primary_unit_count"] != right_row["primary_unit_count"]
        ):
            continue
        eligible.append(key)
    agreement_count = sum(left[key][signature_key] == right[key][signature_key] for key in eligible)
    payload = ratio_payload(agreement_count, len(eligible))
    payload.update(
        {
            "eligible_regions": len(eligible),
            "agreement_count": agreement_count,
            "require_same_internal_ssnv_set": require_same_sites,
            "require_same_primary_unit_count": require_same_primary_unit_count,
        }
    )
    return payload


def prequalification_count(
    keys: Iterable[str],
    left: Mapping[str, Dict[str, Any]],
    right: Mapping[str, Dict[str, Any]],
    signature_key: str,
) -> int:
    return sum(
        left[key][signature_key] is not None and right[key][signature_key] is not None
        for key in keys
    )


def build_by_chromosome(
    *,
    left_regions: Mapping[str, Dict[str, Any]],
    right_regions: Mapping[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    rows = []
    for chrom in AUTOSOMES:
        left_keys = {key for key, value in left_regions.items() if value["chrom"] == chrom}
        right_keys = {key for key, value in right_regions.items() if value["chrom"] == chrom}
        common = sorted(left_keys & right_keys)
        same_sites = [
            key
            for key in common
            if left_regions[key]["internal_sites"] == right_regions[key]["internal_sites"]
        ]
        shape = signature_metric(
            same_sites,
            left_regions,
            right_regions,
            "shape_multiset",
            require_same_sites=False,
            require_same_primary_unit_count=False,
        )
        edge = signature_metric(
            same_sites,
            left_regions,
            right_regions,
            "exact_labeled_edge_multiset",
            require_same_sites=False,
            require_same_primary_unit_count=False,
        )
        overlap = set_overlap(left_keys, right_keys)
        rows.append(
            {
                "chrom": chrom,
                "left_regions": len(left_keys),
                "right_regions": len(right_keys),
                "exact_region_intersection": len(common),
                "exact_region_union": overlap["union"],
                "exact_region_jaccard": overlap["jaccard"],
                "same_internal_ssnv_regions": len(same_sites),
                "shape_eligible": shape["eligible_regions"],
                "shape_agreement": shape["agreement_count"],
                "exact_edge_eligible": edge["eligible_regions"],
                "exact_edge_agreement": edge["agreement_count"],
            }
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate HCC1395 current-v5 exact topology signatures."
    )
    parser.add_argument("--hcc-sidecar", required=True, type=Path)
    parser.add_argument("--dorado-sidecar", required=True, type=Path)
    parser.add_argument("--sidecar-index", required=True, type=Path)
    parser.add_argument("--current-summary", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def analyze(args: argparse.Namespace) -> Dict[str, Any]:
    script_path = Path(__file__).resolve()
    run_root = resolve_existing(args.run_root, "directory")
    sidecar_paths = {
        "HCC1395": resolve_existing(args.hcc_sidecar, "file"),
        "HCC1395_DORADO": resolve_existing(args.dorado_sidecar, "file"),
    }
    sidecar_index_path = resolve_existing(args.sidecar_index, "file")
    current_summary_path = resolve_existing(args.current_summary, "file")
    manifest_path = resolve_existing(args.input_manifest, "file")
    output_path = args.output.expanduser().resolve()

    sidecars = {sample: load_json(path) for sample, path in sidecar_paths.items()}
    sidecar_index = load_json(sidecar_index_path)
    current_summary = load_json(current_summary_path)
    manifest = load_json(manifest_path)
    shared = validate_shared_provenance(
        sidecars=sidecars,
        sidecar_paths=sidecar_paths,
        sidecar_index=sidecar_index,
        sidecar_index_path=sidecar_index_path,
        current_summary=current_summary,
        current_summary_path=current_summary_path,
        run_root=run_root,
        manifest=manifest,
        manifest_path=manifest_path,
    )
    samples = {
        sample: load_sample_contract(
            sample=sample,
            sidecar=sidecars[sample],
            sidecar_path=sidecar_paths[sample],
            run_root=run_root,
        )
        for sample in EXPECTED_SAMPLES
    }
    left, right = samples["HCC1395"], samples["HCC1395_DORADO"]
    left_regions, right_regions = left["regions"], right["regions"]
    common = sorted(set(left_regions) & set(right_regions))
    same_sites = [
        key
        for key in common
        if left_regions[key]["internal_sites"] == right_regions[key]["internal_sites"]
    ]

    shape_all = signature_metric(
        common,
        left_regions,
        right_regions,
        "shape_multiset",
        require_same_sites=False,
        require_same_primary_unit_count=False,
    )
    shape_same_sites = signature_metric(
        common,
        left_regions,
        right_regions,
        "shape_multiset",
        require_same_sites=True,
        require_same_primary_unit_count=False,
    )
    shape_same_sites_units = signature_metric(
        common,
        left_regions,
        right_regions,
        "shape_multiset",
        require_same_sites=True,
        require_same_primary_unit_count=True,
    )
    edge_same_sites = signature_metric(
        common,
        left_regions,
        right_regions,
        "exact_labeled_edge_multiset",
        require_same_sites=True,
        require_same_primary_unit_count=False,
    )
    edge_same_sites_units = signature_metric(
        common,
        left_regions,
        right_regions,
        "exact_labeled_edge_multiset",
        require_same_sites=True,
        require_same_primary_unit_count=True,
    )
    exact_candidate_prequalified = prequalification_count(
        common,
        left_regions,
        right_regions,
        "exact_labeled_edge_multiset",
    )

    same_site_payload = ratio_payload(len(same_sites), len(common))
    same_site_payload.update(
        {
            "exact_common_regions": len(common),
            "same_internal_ssnv_regions": len(same_sites),
        }
    )
    region_overlap = set_overlap(set(left_regions), set(right_regions))
    site_overlap = set_overlap(left["site_universe"], right["site_universe"])
    by_chromosome = build_by_chromosome(
        left_regions=left_regions, right_regions=right_regions
    )

    checks = {
        "shared_provenance_hashes_match": True,
        "sample_contracts_pass": all(all(item["checks"].values()) for item in samples.values()),
        "exact_region_intersection_conserved_by_chromosome": sum(
            row["exact_region_intersection"] for row in by_chromosome
        )
        == region_overlap["intersection"],
        "same_internal_ssnv_conserved_by_chromosome": sum(
            row["same_internal_ssnv_regions"] for row in by_chromosome
        )
        == len(same_sites),
        "shape_eligible_conserved_by_chromosome": sum(
            row["shape_eligible"] for row in by_chromosome
        )
        == shape_same_sites["eligible_regions"],
        "shape_agreement_conserved_by_chromosome": sum(
            row["shape_agreement"] for row in by_chromosome
        )
        == shape_same_sites["agreement_count"],
        "edge_eligible_conserved_by_chromosome": sum(
            row["exact_edge_eligible"] for row in by_chromosome
        )
        == edge_same_sites["eligible_regions"],
        "edge_agreement_conserved_by_chromosome": sum(
            row["exact_edge_agreement"] for row in by_chromosome
        )
        == edge_same_sites["agreement_count"],
        "conditional_shape_denominator_is_nested": shape_same_sites_units["denominator"]
        <= shape_same_sites["denominator"]
        <= shape_all["denominator"],
        "conditional_edge_denominator_is_nested": edge_same_sites_units["denominator"]
        <= edge_same_sites["denominator"]
        <= exact_candidate_prequalified,
        "hp_labels_not_used_as_pairing_key": True,
    }
    require(all(checks.values()), f"pair checks failed: {checks}")

    manifest_records = shared.pop("manifest_records")
    command_argv = [
        sys.executable,
        str(script_path),
        "--hcc-sidecar",
        str(sidecar_paths["HCC1395"]),
        "--dorado-sidecar",
        str(sidecar_paths["HCC1395_DORADO"]),
        "--sidecar-index",
        str(sidecar_index_path),
        "--current-summary",
        str(current_summary_path),
        "--run-root",
        str(run_root),
        "--input-manifest",
        str(manifest_path),
        "--output",
        str(output_path),
    ]
    stdout_lines = [
        f"INPUT HCC1395={sidecar_paths['HCC1395']}",
        f"INPUT HCC1395_DORADO={sidecar_paths['HCC1395_DORADO']}",
        f"INPUT current_summary={current_summary_path}",
        f"INPUT run_root={run_root}",
        f"OUTPUT artifact={output_path}",
        (
            "RESULT region_jaccard={:.6f} same_internal_ssnv={}/{} "
            "shape={}/{} exact_edges={}/{}"
        ).format(
            region_overlap["jaccard"],
            len(same_sites),
            len(common),
            shape_same_sites["agreement_count"],
            shape_same_sites["denominator"],
            edge_same_sites["agreement_count"],
            edge_same_sites["denominator"],
        ),
        f"CHECKS all_checks_pass={str(all(checks.values())).lower()}",
    ]

    artifact: Dict[str, Any] = {
        "schema_name": "intersubmod.hcc1395_exact_signature_validation",
        "schema_version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": "HCC1395 canonical vs HCC1395_DORADO / GRCh38 chr1-22 / current canonical v5",
        "claim_ceiling": (
            "cross-platform pipeline-output concordance for regional mutation-state topology; "
            "not biological sample-identity proof, confirmed clone, true ancestry, CCF, or probability"
        ),
        "comparison_contract": {
            "region_pairing": "exact stored region key chrom:start-end",
            "internal_ssnv_pairing": (
                "exact equality of the per-region genomic-position set reconstructed from "
                "region-view obs_col_coverage"
            ),
            "primary_unit_pairing": (
                "compare a sorted multiset of primary-unit signatures; preserve multiplicity "
                "but ignore HP1/HP2 family names to tolerate label flips"
            ),
            "fixed_unit_supplement": (
                "n_trees=1 unit payloads are intentionally absent from the sidecar and are "
                "supplemented from the SHA-bound layered reconstruction"
            ),
            "shape_signature": (
                "unlabeled rooted-forest shape; comparable as geometry but does not encode mutation identity"
            ),
            "exact_labeled_edges": (
                "sorted mutation-state parent-child edges; recurrence metadata is outside this edge-only metric"
            ),
            "edge_comparison_gate": (
                "exact region + identical internal sSNV set + one selected exact candidate in each dataset"
            ),
            "candidate_completeness": (
                "same fail-closed unit_complete predicate as current-v5 builder"
            ),
        },
        "provenance": {
            "analysis_script": str(script_path),
            "analysis_script_sha256": sha256_file(script_path),
            "sidecar_index": str(sidecar_index_path),
            "sidecar_index_sha256": shared["sidecar_index_sha256"],
            "sidecar_builder": (sidecar_index.get("provenance") or {}).get("builder"),
            "sidecar_builder_sha256": (sidecar_index.get("provenance") or {}).get(
                "builder_sha256"
            ),
            "topology_solver": (sidecar_index.get("provenance") or {}).get("solver"),
            "topology_solver_sha256": (sidecar_index.get("provenance") or {}).get(
                "solver_sha256"
            ),
            "current_summary": str(current_summary_path),
            "current_summary_sha256": shared["current_summary_sha256"],
            "run_root": str(run_root),
            "input_manifest": str(manifest_path),
            "input_manifest_sha256": shared["input_manifest_sha256"],
            "samples": {
                sample: {
                    "sidecar": str(samples[sample]["sidecar_path"]),
                    "sidecar_sha256": samples[sample]["sidecar_sha256"],
                    "layered_reconstruction": str(samples[sample]["reconstruction_path"]),
                    "layered_reconstruction_sha256": samples[sample]["reconstruction_sha256"],
                    "layered_region_view": str(samples[sample]["region_view_path"]),
                    "layered_region_view_sha256": samples[sample]["region_view_sha256"],
                    "biological_id": manifest_records[sample].get("biological_id"),
                    "platform": manifest_records[sample].get("platform"),
                    "replicate_role": manifest_records[sample].get("replicate_role"),
                    "alignment_payload_identity_sha256": manifest_records[sample]
                    .get("alignment_payload", {})
                    .get("storage_identity_v1", {})
                    .get("identity_sha256"),
                    "tree_vcf_sha256": manifest_records[sample]
                    .get("somatic", {})
                    .get("tree_vcf", {})
                    .get("identity", {})
                    .get("sha256"),
                    "read_tag_sidecar_sha256": manifest_records[sample]
                    .get("read_tags", {})
                    .get("sidecar", {})
                    .get("identity", {})
                    .get("sha256"),
                }
                for sample in EXPECTED_SAMPLES
            },
        },
        "sample_contracts": {
            sample: {
                "counts": samples[sample]["counts"],
                "signature_source_counts": samples[sample]["signature_source_counts"],
                "checks": samples[sample]["checks"],
                "all_checks_pass": all(samples[sample]["checks"].values()),
            }
            for sample in EXPECTED_SAMPLES
        },
        "region_universe": region_overlap,
        "retained_ssnv_universe": site_overlap,
        "internal_ssnv_pairing": {
            "same_set_within_exact_common_region": same_site_payload,
        },
        "shape_agreement": {
            "exact_region_both_shape_resolved": shape_all,
            "same_internal_ssnv_both_shape_resolved": shape_same_sites,
            "same_internal_ssnv_same_primary_unit_count": shape_same_sites_units,
        },
        "exact_labeled_edge_agreement": {
            "exact_region_both_candidate_unique_prequalification": {
                "regions": exact_candidate_prequalified,
                "comparison_valid": False,
                "reason": "internal sSNV equality is required before labeled edges are comparable",
            },
            "same_internal_ssnv_both_candidate_unique": edge_same_sites,
            "same_internal_ssnv_same_primary_unit_count": edge_same_sites_units,
        },
        "by_chromosome": by_chromosome,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "execution": {
            "command_argv": command_argv,
            "command_shell": " ".join(shlex.quote(value) for value in command_argv),
            "output": str(output_path),
            "stdout_lines": stdout_lines,
        },
    }
    return artifact


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        artifact = analyze(args)
        output_path = args.output.expanduser().resolve()
        write_json(output_path, artifact)
    except ContractError as exc:
        print(f"ERROR {exc}", file=sys.stderr)
        return 2
    for line in artifact["execution"]["stdout_lines"]:
        print(line)
    print(f"SHA256 artifact={sha256_file(output_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
