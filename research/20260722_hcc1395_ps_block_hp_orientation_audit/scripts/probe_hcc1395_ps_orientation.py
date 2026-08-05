#!/usr/bin/env python3
"""Probe HCC1395 topology sensitivity to independent PS-block HP flips.

The first observed PS is fixed as the orientation anchor.  Every remaining PS
is independently tested in native and HP1/HP2-swapped orientation.  This removes
the biologically irrelevant whole-region global swap while exposing whether the
legacy ``region x HP-family`` pooling depends on arbitrary between-block labels.

This is a sensitivity analysis, not an estimator of the true cross-block phase.
Reads without PS remain in their native family and are reported separately.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import importlib.util
import itertools
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from types import ModuleType
from typing import Any

import pysam


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def hp_family(raw_hp: Any) -> str:
    if raw_hp is None:
        return "none"
    value = str(raw_hp)
    if value.startswith("1"):
        return "1"
    if value.startswith("2"):
        return "2"
    if value in {"3", "4"}:
        return value
    return "none"


def canonical_shape(edges: list[tuple[str, str]]) -> str:
    children: defaultdict[str, list[str]] = defaultdict(list)
    child_nodes: set[str] = set()
    all_nodes: set[str] = set()
    for parent, child in edges:
        children[parent].append(child)
        child_nodes.add(child)
        all_nodes.update((parent, child))
    roots = sorted(all_nodes - child_nodes)

    def visit(node: str) -> str:
        return "(" + "".join(sorted(visit(child) for child in children[node])) + ")"

    return "|".join(sorted(visit(root) for root in roots)) if roots else "()"


def shape_features(edges: list[tuple[str, str]]) -> tuple[bool, bool]:
    nodes: set[str] = {"ROOT"}
    children: defaultdict[str, list[str]] = defaultdict(list)
    incoming = Counter()
    for parent, child in edges:
        nodes.update((parent, child))
        children[parent].append(child)
        incoming[child] += 1
    roots = sorted(node for node in nodes if incoming[node] == 0)
    depth = {root: 0 for root in roots}
    for _ in range(len(nodes) + 1):
        changed = False
        for parent, child in edges:
            if parent not in depth:
                continue
            candidate = depth[parent] + 1
            if candidate > depth.get(child, -1):
                depth[child] = candidate
                changed = True
        if not changed:
            break
    return max(depth.values(), default=0) >= 2, max((len(children[node]) for node in nodes), default=0) >= 2


def morphology(features: tuple[bool, bool]) -> str:
    direct, sister = features
    if direct and sister:
        return "direct_and_sister"
    if direct:
        return "direct_chain"
    if sister:
        return "sister_branch"
    return "single_no_within_hp_relation"


def structural_class(candidate_count: int, topology_count: int) -> str:
    if candidate_count == 1 and topology_count == 1:
        return "exact_and_topology_unique"
    if candidate_count > 1 and topology_count == 1:
        return "topology_unique_exact_multiple"
    if candidate_count > 1 and topology_count > 1:
        return "topology_multiple_exact_multiple"
    return "invalid_or_no_primary"


def semantic_digest(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def normalized_nested_counts(value: dict[str, dict[str, int]] | None) -> dict[str, dict[str, int]]:
    return {
        str(family): {str(pattern): int(count) for pattern, count in sorted(counts.items())}
        for family, counts in sorted((value or {}).items())
    }


def find_group(mlhp_glob: str, region: str) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    paths = [Path(path) for path in sorted(glob.glob(mlhp_glob))]
    if not paths:
        raise FileNotFoundError(mlhp_glob)
    found: list[tuple[dict[str, Any], dict[str, Any]]] = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            document = json.load(handle)
        for group in document.get("groups", []):
            key = f"{group['chrom']}:{group['start']}-{group['end']}"
            if key == region:
                found.append((document, group))
    if len(found) != 1:
        raise ValueError(f"expected one MLHP group for {region}; observed {len(found)}")
    return found[0][0], found[0][1], paths


def manifest_entry(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    entries = [entry for entry in manifest.get("samples", []) if entry.get("sample") == "HCC1395"]
    if len(entries) != 1:
        raise ValueError(f"HCC1395 manifest cardinality: {len(entries)}")
    return entries[0]


def load_variants(vcf_path: Path, chrom: str, positions: list[int]) -> list[tuple[int, str, str, str]]:
    wanted = set(positions)
    observed: dict[int, tuple[str, str]] = {}
    with pysam.VariantFile(str(vcf_path)) as source:
        for record in source.fetch(chrom, min(positions) - 1, max(positions)):
            if int(record.pos) not in wanted:
                continue
            alts = tuple(record.alts or ())
            if len(record.ref) != 1 or len(alts) != 1 or len(alts[0]) != 1:
                raise ValueError(f"non-biallelic SNV at {chrom}:{record.pos}")
            observed[int(record.pos)] = (record.ref.upper(), alts[0].upper())
    missing = sorted(wanted - set(observed))
    if missing:
        raise ValueError(f"positions absent from tree VCF: {missing}")
    return [(position, observed[position][0], observed[position][1], "PASS") for position in positions]


def aggregate_patterns(
    reads: list[dict[str, Any]], flipped_ps: set[str], minread: int
) -> tuple[dict[str, Any], dict[str, dict[str, int]], dict[str, dict[str, int]]]:
    subread: defaultdict[str, Counter[str]] = defaultdict(Counter)
    full: defaultdict[str, Counter[str]] = defaultdict(Counter)
    reads_by_family = Counter()
    for read in reads:
        family = read["family"]
        if read["phase_set"] in flipped_ps and family in {"1", "2"}:
            family = "2" if family == "1" else "1"
        reads_by_family[family] += 1
        vector = read["vector"]
        subread[family][vector] += 1
        if "X" not in vector:
            full[family][vector] += 1
    retained_subread = {
        family: {pattern: count for pattern, count in counts.items() if count >= minread}
        for family, counts in subread.items()
    }
    retained_full = {
        family: {pattern: count for pattern, count in counts.items() if count >= minread}
        for family, counts in full.items()
    }
    retained_subread = {family: values for family, values in retained_subread.items() if values}
    retained_full = {family: values for family, values in retained_full.items() if values}
    return (
        {"reads_by_family": dict(sorted(reads_by_family.items()))},
        normalized_nested_counts(retained_full),
        normalized_nested_counts(retained_subread),
    )


def solve_configuration(
    solver: ModuleType,
    reads: list[dict[str, Any]],
    flipped_ps: set[str],
    minread: int,
    k: int,
) -> dict[str, Any]:
    counts, full_by_family, partial_by_family = aggregate_patterns(reads, flipped_ps, minread)
    units: dict[str, Any] = {}
    morphology_features_by_family: list[list[tuple[bool, bool]]] = []
    for family in ("1", "2"):
        full = full_by_family.get(family, {})
        partial_counts = partial_by_family.get(family, {})
        mutation_bearing = any("A" in pattern for pattern in full) or any(
            "A" in pattern for pattern in partial_counts
        )
        if not mutation_bearing:
            continue
        result = solver.enumerate_min_trees(full, list(partial_counts), k, tree_cap=0)
        representatives: dict[str, dict[str, Any]] = {}
        for tree in result.get("trees", []):
            signature = canonical_shape(tree["edges"])
            representatives.setdefault(
                signature,
                {
                    "shape_signature": signature,
                    "edges": [list(edge) for edge in tree["edges"]],
                    "morphology": morphology(shape_features(tree["edges"])),
                },
            )
        signatures = set(representatives)
        feature_options = sorted({shape_features(tree["edges"]) for tree in result.get("trees", [])})
        units[family] = {
            "n_reads": int(counts["reads_by_family"].get(family, 0)),
            "full_patterns": full,
            "partial_patterns": partial_counts,
            "n_full_patterns": len(full),
            "n_partial_patterns": len(partial_counts),
            "n_trees": int(result.get("n_trees", 0)),
            "n_topologies": len(signatures),
            "n_hidden": int(result.get("n_hidden", 0)),
            "capped": bool(result.get("capped")),
            "trees_complete": bool(result.get("trees_complete")),
            "morphology_options": [morphology(option) for option in feature_options],
            "topology_representatives": list(representatives.values()),
        }
        morphology_features_by_family.append(feature_options)

    if not units:
        return {
            "flipped_ps": sorted(flipped_ps),
            **counts,
            "units": {},
            "region": {"n_primary_families": 0, "T": None, "Topo": None, "structural_class": "no_primary"},
        }
    if any(unit["capped"] or not unit["trees_complete"] for unit in units.values()):
        complete = False
        candidate_count = None
        topology_count = None
        region_morphologies: list[str] = []
    else:
        complete = True
        candidate_count = math.prod(unit["n_trees"] for unit in units.values())
        topology_count = math.prod(unit["n_topologies"] for unit in units.values())
        region_morphologies = sorted(
            {
                morphology((any(option[0] for option in combination), any(option[1] for option in combination)))
                for combination in itertools.product(*morphology_features_by_family)
            }
        )
    return {
        "flipped_ps": sorted(flipped_ps),
        **counts,
        "units": units,
        "region": {
            "n_primary_families": len(units),
            "candidate_set_complete": complete,
            "T": candidate_count,
            "Topo": topology_count,
            "structural_class": (
                structural_class(candidate_count, topology_count) if complete else "incomplete"
            ),
            "morphology_options": region_morphologies,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--region", required=True)
    parser.add_argument("--mlhp-glob", required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--frozen-linkage", type=Path, required=True)
    parser.add_argument("--frozen-solver", type=Path, required=True)
    parser.add_argument("--current-topology", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    document, group, mlhp_paths = find_group(args.mlhp_glob, args.region)
    entry = manifest_entry(args.manifest)
    linkage = load_module("frozen_sm_linkage_genomewide", args.frozen_linkage)
    solver = load_module("frozen_tree_enumeration_solver", args.frozen_solver)
    params = document.get("params", {})
    minread = int(params.get("MINREAD", 3))
    linkage.MAPQ_MIN = int(params.get("MAPQ_MIN", 20))
    linkage.BASEQ_MIN = int(params.get("BASEQ_MIN", 0))
    linkage.READ_TAG_SIDECAR = entry["read_tags"]["sidecar"]["path"]
    linkage._TAG_TABIX = None

    positions = [int(value) for value in group["positions"]]
    variants = load_variants(
        Path(entry["somatic"]["tree_vcf"]["path"]), group["chrom"], positions
    )
    with pysam.AlignmentFile(entry["alignment_payload"]["path"], "rb") as bam:
        calls, hp, ps = linkage.per_read_calls(bam, group["chrom"], variants)
    diagnostics = linkage.last_tag_diagnostics()

    reads: list[dict[str, Any]] = []
    raw_phase_set_counts = Counter()
    raw_phase_set_hp_counts = Counter()
    for read_key, read_calls in calls.items():
        vector = "".join(
            "A" if read_calls.get(position) == "ALT" else "R" if read_calls.get(position) == "REF" else "X"
            for position in positions
        )
        if set(vector) == {"X"}:
            continue
        raw_hp = str(hp.get(read_key, "."))
        phase_set = None if read_key not in ps else str(ps[read_key])
        if phase_set is not None:
            raw_phase_set_counts[phase_set] += 1
            raw_phase_set_hp_counts[f"{phase_set}|{raw_hp}"] += 1
        reads.append(
            {
                "family": hp_family(hp.get(read_key)),
                "raw_hp": raw_hp,
                "phase_set": phase_set,
                "vector": vector,
            }
        )

    native_counts, native_full, native_partial = aggregate_patterns(reads, set(), minread)
    reproduction_checks = {
        "alignment_exposures_match": int(diagnostics.get("alignment_group_exposures", -1))
        == int(group.get("read_tag_diagnostics", {}).get("alignment_group_exposures", -2)),
        "exact_matches_match_exposures": int(diagnostics.get("sidecar_exact_matches", -1))
        == int(diagnostics.get("alignment_group_exposures", -2)),
        "sidecar_missing_zero": int(diagnostics.get("sidecar_missing", -1)) == 0,
        "sidecar_conflicts_zero": int(diagnostics.get("sidecar_conflicts", -1)) == 0,
        "phase_set_counts_match": dict(raw_phase_set_counts) == (group.get("phase_set_counts") or {}),
        "phase_set_hp_counts_match": dict(raw_phase_set_hp_counts)
        == (group.get("phase_set_HP_counts") or {}),
        "reads_by_hp_match": native_counts["reads_by_family"] == (group.get("reads_by_hp") or {}),
        "populations_by_hp_match": native_full
        == normalized_nested_counts(group.get("populations_by_hp") or {}),
        "subread_groups_by_hp_match": native_partial
        == normalized_nested_counts(group.get("subread_groups_by_hp") or {}),
    }

    phase_sets = sorted(raw_phase_set_counts)
    if len(phase_sets) < 2:
        raise ValueError(f"region is not mixed-PS: {args.region}")
    anchor = phase_sets[0]
    configurations: list[dict[str, Any]] = []
    for bits in itertools.product((0, 1), repeat=len(phase_sets) - 1):
        flipped = {phase_sets[offset + 1] for offset, bit in enumerate(bits) if bit}
        configuration = solve_configuration(solver, reads, flipped, minread, len(positions))
        configuration["configuration_id"] = f"F{len(configurations)}"
        configuration["family_unit_signature_sha256"] = semantic_digest(configuration["units"])
        configurations.append(configuration)

    with args.current_topology.open(encoding="utf-8") as handle:
        topology_document = json.load(handle)
    current_matches = [row for row in topology_document.get("regions", []) if row.get("region") == args.region]
    if len(current_matches) != 1:
        raise ValueError(f"current-v5 region cardinality: {len(current_matches)}")
    current = current_matches[0]
    native = configurations[0]["region"]
    current_reproduction_checks = {
        "native_T_matches_current": native["T"] == current.get("C"),
        "native_Topo_matches_current": native["Topo"] == current.get("Topo"),
        "native_structural_class_matches_current": native["structural_class"]
        == current.get("structural_class"),
    }

    comparison_fields = ("n_primary_families", "T", "Topo", "structural_class", "morphology_options")
    changes = []
    for configuration in configurations[1:]:
        changed = {
            field: {
                "native": native.get(field),
                "alternative": configuration["region"].get(field),
            }
            for field in comparison_fields
            if native.get(field) != configuration["region"].get(field)
        }
        if configurations[0]["family_unit_signature_sha256"] != configuration["family_unit_signature_sha256"]:
            changed["family_unit_signature"] = {
                "native": configurations[0]["family_unit_signature_sha256"],
                "alternative": configuration["family_unit_signature_sha256"],
            }
        changes.append(
            {
                "configuration_id": configuration["configuration_id"],
                "flipped_ps": configuration["flipped_ps"],
                "changed_fields": changed,
                "orientation_sensitive": bool(changed),
            }
        )

    output = {
        "schema_name": "intersubmod.hcc1395_ps_orientation_probe",
        "schema_version": "1.0.0",
        "scope": {
            "sample": "HCC1395",
            "region": args.region,
            "positions": positions,
            "k": len(positions),
            "minread": minread,
            "orientation_anchor_ps": anchor,
            "orientation_space": f"2^({len(phase_sets)}-1)={len(configurations)} relative orientations",
            "claim_boundary": (
                "sensitivity to arbitrary relative PS labels; not identification of the true orientation"
            ),
        },
        "inputs": {
            "mlhp_parts": [str(path.resolve()) for path in mlhp_paths],
            "manifest": str(args.manifest.resolve()),
            "bam": entry["alignment_payload"]["path"],
            "tree_vcf": entry["somatic"]["tree_vcf"]["path"],
            "read_tag_sidecar": entry["read_tags"]["sidecar"]["path"],
            "frozen_linkage": str(args.frozen_linkage.resolve()),
            "frozen_solver": str(args.frozen_solver.resolve()),
            "current_topology": str(args.current_topology.resolve()),
        },
        "read_evidence": {
            "alignment_exposures": int(diagnostics.get("alignment_group_exposures", 0)),
            "exact_sidecar_matches": int(diagnostics.get("sidecar_exact_matches", 0)),
            "missing_sidecar": int(diagnostics.get("sidecar_missing", 0)),
            "conflicting_sidecar": int(diagnostics.get("sidecar_conflicts", 0)),
            "reads_used": len(reads),
            "reads_with_known_ps": sum(read["phase_set"] is not None for read in reads),
            "reads_without_ps": sum(read["phase_set"] is None for read in reads),
            "phase_set_counts": dict(sorted(raw_phase_set_counts.items())),
            "phase_set_HP_counts": dict(sorted(raw_phase_set_hp_counts.items())),
        },
        "canonical_reproduction": {
            "mlhp_checks": reproduction_checks,
            "current_v5_checks": current_reproduction_checks,
            "all_pass": all(reproduction_checks.values()) and all(current_reproduction_checks.values()),
        },
        "configurations": configurations,
        "comparisons_to_native": changes,
        "orientation_sensitive": any(change["orientation_sensitive"] for change in changes),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(
        json.dumps(
            {
                "output_json": str(args.output_json),
                "region": args.region,
                "phase_sets": phase_sets,
                "configuration_count": len(configurations),
                "canonical_reproduction_all_pass": output["canonical_reproduction"]["all_pass"],
                "orientation_sensitive": output["orientation_sensitive"],
                "comparisons_to_native": changes,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
