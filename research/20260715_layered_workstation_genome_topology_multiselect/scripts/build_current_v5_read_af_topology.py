#!/usr/bin/env python3
"""Build current-canonical-v5 read-AF ranking and topology morphology sidecars.

The output is descriptive engineering evidence.  Exact candidate ties use
Fraction-valued family-specific ALT/(REF+ALT) and the ordering score is not a
calibrated likelihood, posterior, CCF, or biological clone confirmation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime
from fractions import Fraction
from pathlib import Path


SAMPLE_ORDER = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
STRUCTURAL_KEYS = (
    "exact_and_topology_unique",
    "topology_unique_exact_multiple",
    "topology_multiple_exact_multiple",
    "incomplete",
)
SELECTION_KEYS = (
    "structural_exact_unique",
    "read_af_unique_first",
    "read_af_tied_same_topology",
    "read_af_tied_different_topology",
    "read_af_unavailable",
    "incomplete",
)
MORPHOLOGY_KEYS = (
    "single_no_within_hp_relation",
    "direct_chain",
    "sister_branch",
    "direct_and_sister",
    "unresolved",
)


def fail(message: str) -> None:
    raise SystemExit(f"CURRENT-V5 READ-AF BUILD FAILED: {message}")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object, *, compact: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {"ensure_ascii": False}
    if compact:
        kwargs["separators"] = (",", ":")
    else:
        kwargs["indent"] = 2
    path.write_text(json.dumps(value, **kwargs) + "\n", encoding="utf-8")


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "fraction": f"{value.numerator}/{value.denominator}",
        "value": float(value),
    }


def canonical_shape(edges: list[list[str]] | list[tuple[str, str]]) -> str:
    children: dict[str, list[str]] = defaultdict(list)
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


def shape_id(signature: str) -> str:
    return "TS-" + hashlib.sha1(signature.encode("utf-8")).hexdigest()[:10]


def tree_id(edges: list[list[str]] | list[tuple[str, str]]) -> str:
    encoded = json.dumps(sorted([list(edge) for edge in edges]), separators=(",", ":")).encode()
    return "TR-" + hashlib.sha1(encoded).hexdigest()[:12]


def unlabel(node: str) -> frozenset[int]:
    if node == "ROOT":
        return frozenset()
    value = node[2:] if node.startswith("H_") else node
    return frozenset(index for index, allele in enumerate(value) if allele == "A")


def exact_ordering_score(
    edges: list[list[str]] | list[tuple[str, str]], read_af: dict[int, Fraction]
) -> tuple[Fraction, list[Fraction]]:
    score = Fraction(0, 1)
    comparisons: list[Fraction] = []
    for parent, child in edges:
        parent_set, child_set = unlabel(parent), unlabel(child)
        acquired = child_set - parent_set
        if len(acquired) != 1:
            continue
        child_mutation = next(iter(acquired))
        for ancestor_mutation in parent_set:
            delta = read_af[ancestor_mutation] - read_af[child_mutation]
            score += delta
            comparisons.append(delta)
    return score, comparisons


def exact_read_af(
    col_coverage: dict, positions: list[int]
) -> tuple[dict[int, Fraction] | None, list[dict], str | None]:
    values: dict[int, Fraction] = {}
    rows: list[dict] = []
    for index, position in enumerate(positions):
        counts = col_coverage.get(str(position)) or col_coverage.get(position)
        if not counts:
            return None, rows, f"missing coverage at {position}"
        n_ref, n_alt = int(counts[0]), int(counts[1])
        total = n_ref + n_alt
        if total <= 0:
            return None, rows, f"zero denominator at {position}"
        fraction = Fraction(n_alt, total)
        values[index] = fraction
        rows.append(
            {
                "site_index": index + 1,
                "position": int(position),
                "ref_reads": n_ref,
                "alt_reads": n_alt,
                "read_af": fraction_payload(fraction),
            }
        )
    return values, rows, None


def shape_features(edges: list[list[str]] | list[tuple[str, str]]) -> dict:
    nodes: set[str] = {"ROOT"}
    children: dict[str, list[str]] = defaultdict(list)
    incoming: Counter = Counter()
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
    max_depth = max(depth.values(), default=0)
    max_outdegree = max((len(children[node]) for node in nodes), default=0)
    return {
        "max_depth": max_depth,
        "max_outdegree": max_outdegree,
        "has_direct_chain": max_depth >= 2,
        "has_sister_branch": max_outdegree >= 2,
    }


def morphology_class(has_direct: bool, has_sister: bool) -> str:
    if has_direct and has_sister:
        return "direct_and_sister"
    if has_direct:
        return "direct_chain"
    if has_sister:
        return "sister_branch"
    return "single_no_within_hp_relation"


def unit_complete(unit: dict) -> bool:
    return (
        unit.get("analysis_candidate_set_complete") is True
        and unit.get("capped") is False
        and unit.get("verification_status") == "full_pass"
        and unit.get("verification_complete") is True
        and unit.get("verify_pass") is True
        and int(unit.get("analysis_trees_generated") or -1) == int(unit.get("n_trees") or 0)
        and unit.get("n_distinct_shapes_exact") is not None
    )


def structural_class(units: list[dict]) -> tuple[str, int | None, int | None]:
    if not units:
        return "no_primary_lineage", None, None
    if not all(unit_complete(unit) for unit in units):
        return "incomplete", None, None
    candidate_count = math.prod(int(unit["n_trees"]) for unit in units)
    topology_count = math.prod(int(unit["n_distinct_shapes_exact"]) for unit in units)
    if candidate_count == 1 and topology_count == 1:
        return "exact_and_topology_unique", candidate_count, topology_count
    if candidate_count > 1 and topology_count == 1:
        return "topology_unique_exact_multiple", candidate_count, topology_count
    if candidate_count > 1 and topology_count > 1:
        return "topology_multiple_exact_multiple", candidate_count, topology_count
    fail(f"impossible structural state C={candidate_count}, Topo={topology_count}")


def load_groups(sample_dir: Path) -> dict[str, dict]:
    groups: dict[str, dict] = {}
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        for group in read_json(path).get("groups", []):
            key = f"{group['chrom']}:{group['start']}-{group['end']}"
            if key in groups:
                fail(f"duplicate multilocus region: {key}")
            groups[key] = group
    return groups


def analyze_ambiguous_unit(unit: dict, group: dict | None, solver) -> dict:
    base = {
        "family": str(unit["family"]),
        "n_trees_expected": int(unit["n_trees"]),
        "n_shapes_expected": int(unit["n_distinct_shapes_exact"]),
    }
    if group is None:
        return {**base, "status": "missing_group", "reason": "region absent from mlhp parts"}
    family = str(unit["family"])
    full = (group.get("populations_by_hp", {}) or {}).get(family, {}) or {}
    partial = list(((group.get("subread_groups_by_hp", {}) or {}).get(family, {}) or {}).keys())
    positions = [int(value) for value in group.get("positions", [])]
    result = solver.enumerate_min_trees(full, partial, len(positions), tree_cap=0)
    signatures = [canonical_shape(tree.get("edges", [])) for tree in result.get("trees", [])]
    enumeration_ok = (
        result.get("capped") is False
        and result.get("trees_complete") is True
        and int(result.get("n_trees") or 0) == base["n_trees_expected"]
        and len(set(signatures)) == base["n_shapes_expected"]
    )
    if not enumeration_ok:
        return {
            **base,
            "status": "candidate_or_shape_mismatch",
            "n_trees_enumerated": int(result.get("n_trees") or 0),
            "n_trees_stored": len(result.get("trees", [])),
            "n_shapes_enumerated": len(set(signatures)),
            "solver_capped": bool(result.get("capped")),
        }
    col_coverage = ((group.get("col_coverage_by_hp", {}) or {}).get(family, {}) or {})
    read_af, read_af_rows, missing_reason = exact_read_af(col_coverage, positions)
    if read_af is None:
        return {
            **base,
            "status": "missing_read_af",
            "reason": missing_reason,
            "n_trees_enumerated": int(result["n_trees"]),
            "n_shapes_enumerated": len(set(signatures)),
            "read_af_site_count_before_failure": len(read_af_rows),
        }

    candidates: list[dict] = []
    for candidate_index, tree in enumerate(result["trees"], 1):
        edges = [list(edge) for edge in tree.get("edges", [])]
        score, comparisons = exact_ordering_score(edges, read_af)
        signature = canonical_shape(edges)
        candidates.append(
            {
                "candidate_index": candidate_index,
                "candidate_id": tree_id(edges),
                "shape_id": shape_id(signature),
                "shape_signature": signature,
                "score": score,
                "n_ancestry_comparisons": len(comparisons),
                "edges": edges,
                "recurrence": list(tree.get("recurrence", [])),
                "n_hidden": int(tree.get("n_hidden") or 0),
            }
        )
    ranked = sorted(candidates, key=lambda candidate: (-candidate["score"], candidate["candidate_id"]))
    max_score = ranked[0]["score"]
    exact_top = [candidate for candidate in ranked if candidate["score"] == max_score]
    exact_top_shapes = sorted({candidate["shape_id"] for candidate in exact_top})
    dense_rank = 0
    previous: Fraction | None = None
    preview: list[dict] = []
    for order, candidate in enumerate(ranked, 1):
        if previous is None or candidate["score"] != previous:
            dense_rank += 1
            previous = candidate["score"]
        if order <= 5:
            preview.append(
                {
                    "candidate_order": order,
                    "score_rank": dense_rank,
                    "candidate_index": candidate["candidate_index"],
                    "candidate_id": candidate["candidate_id"],
                    "shape_id": candidate["shape_id"],
                    "score_fraction": (
                        f"{candidate['score'].numerator}/{candidate['score'].denominator}"
                    ),
                    "score_value": float(candidate["score"]),
                    "is_exact_top": candidate["score"] == max_score,
                    "n_ancestry_comparisons": candidate["n_ancestry_comparisons"],
                }
            )
    representatives: dict[str, dict] = {}
    for candidate in exact_top:
        representatives.setdefault(candidate["shape_id"], candidate)
    top_candidates = []
    for candidate in sorted(representatives.values(), key=lambda item: item["candidate_id"]):
        top_candidates.append(
            {
                "candidate_index": candidate["candidate_index"],
                "candidate_id": candidate["candidate_id"],
                "shape_id": candidate["shape_id"],
                "shape_signature": candidate["shape_signature"],
                "edges": candidate["edges"],
                "recurrence": candidate["recurrence"],
                "n_hidden": candidate["n_hidden"],
                "morphology_features": shape_features(candidate["edges"]),
            }
        )
    top_class = (
        "unique_first"
        if len(exact_top) == 1
        else "tied_first_same_topology"
        if len(exact_top_shapes) == 1
        else "tied_first_different_topology"
    )
    return {
        **base,
        "status": "ranked",
        "n_trees_enumerated": int(result["n_trees"]),
        "n_shapes_enumerated": len(set(signatures)),
        "read_af_site_count": len(read_af_rows),
        "top_score_fraction": f"{max_score.numerator}/{max_score.denominator}",
        "top_score_value": float(max_score),
        "n_top_exact": len(exact_top),
        "n_top_shapes_exact": len(exact_top_shapes),
        "exact_top_class": top_class,
        "candidate_comparison_count_varies": len(
            {candidate["n_ancestry_comparisons"] for candidate in candidates}
        )
        > 1,
        "ranked_preview": preview,
        "top_shape_representatives": top_candidates,
    }


def representative_for_morphology(unit: dict, ranking: dict | None) -> tuple[dict | None, str]:
    if not unit_complete(unit):
        return None, "candidate_incomplete"
    if int(unit["n_distinct_shapes_exact"]) == 1:
        if ranking and ranking.get("status") == "ranked":
            representatives = ranking.get("top_shape_representatives") or []
            if representatives:
                return representatives[0], "structural_topology_unique"
        trees = unit.get("trees") or []
        if trees:
            edges = [list(edge) for edge in trees[0].get("edges", [])]
            signature = canonical_shape(edges)
            return {
                "candidate_index": 1,
                "candidate_id": tree_id(edges),
                "shape_id": shape_id(signature),
                "shape_signature": signature,
                "edges": edges,
                "recurrence": list(trees[0].get("recurrence", [])),
                "n_hidden": int(trees[0].get("n_hidden") or 0),
                "morphology_features": shape_features(edges),
            }, "structural_topology_unique"
        return None, "topology_unique_without_stored_representative"
    if ranking and ranking.get("status") == "ranked" and int(
        ranking.get("n_top_shapes_exact") or 0
    ) == 1:
        representatives = ranking.get("top_shape_representatives") or []
        if representatives:
            return representatives[0], "read_af_topology_unique"
    return None, "multiple_topologies_remain"


def analyze_sample(
    sample: str,
    run_root: Path,
    sample_record: dict,
    solver,
    manifest_path: Path,
    summary_path: Path,
    output_dir: Path,
) -> dict:
    sample_dir = run_root / "samples" / sample
    layered_path = Path(sample_record["paths"]["layered_reconstruction"]).resolve()
    region_view_path = Path(sample_record["paths"]["layered_region_view"]).resolve()
    if layered_path != (sample_dir / f"layered_reconstruction_{sample}.json").resolve():
        fail(f"{sample}: current layered path is outside requested run root")
    if not layered_path.is_file() or not region_view_path.is_file():
        fail(f"{sample}: canonical layered inputs missing")
    if file_sha256(layered_path) != sample_record["sha256"]["layered_reconstruction"]:
        fail(f"{sample}: layered reconstruction SHA-256 mismatch")
    if file_sha256(region_view_path) != sample_record["sha256"]["layered_region_view"]:
        fail(f"{sample}: layered region-view SHA-256 mismatch")

    layered = read_json(layered_path)
    region_view = read_json(region_view_path)
    groups = load_groups(sample_dir)
    primary_by_region: dict[str, list[dict]] = defaultdict(list)
    for unit in layered.get("detail", []):
        if unit.get("is_primary_lineage") is True:
            primary_by_region[unit["region"]].append(unit)
    for units in primary_by_region.values():
        units.sort(key=lambda unit: str(unit["family"]))

    ranking_by_unit: dict[tuple[str, str], dict] = {}
    unit_status = Counter()
    for region, units in sorted(primary_by_region.items()):
        for unit in units:
            family = str(unit["family"])
            if not unit_complete(unit):
                ranking = {"family": family, "status": "incomplete"}
            elif int(unit["n_trees"]) == 1:
                ranking = {"family": family, "status": "fixed_exact_unique"}
            elif str(unit.get("L1_base_class", "")).startswith("recurrence"):
                ranking = {"family": family, "status": "recurrence_not_evaluable"}
            elif unit.get("L1_base_class") == "ambiguous":
                ranking = analyze_ambiguous_unit(unit, groups.get(region), solver)
            else:
                fail(
                    f"{sample} {region} HP{family}: unsupported complete T>1 class "
                    f"{unit.get('L1_base_class')}"
                )
            ranking_by_unit[(region, family)] = ranking
            unit_status[ranking["status"]] += 1

    regions_output: list[dict] = []
    structural_counts = Counter()
    selection_counts = Counter()
    morphology_counts = Counter()
    no_primary = 0
    top_edge_outside_display_cap = 0
    for region in region_view.get("regions", []):
        region_id = region["region"]
        primary_units = primary_by_region.get(region_id, [])
        category, candidate_count, topology_count = structural_class(primary_units)
        structural_counts[category] += 1
        base = {
            "region": region_id,
            "chrom": region["chrom"],
            "start": int(region["start"]),
            "end": int(region["end"]),
            "structural_class": category,
            "C": candidate_count,
            "Topo": topology_count,
            "primary_families": [str(unit["family"]) for unit in primary_units],
        }
        if category == "no_primary_lineage":
            no_primary += 1
            regions_output.append(
                {
                    **base,
                    "selection_class": "no_primary",
                    "selection_note": "no mutation-bearing HP1/HP2 primary unit",
                    "morphology_class": "not_applicable",
                    "morphology_note": "no primary HP component",
                    "units": {},
                }
            )
            continue

        unit_payload = {
            str(unit["family"]): ranking_by_unit[(region_id, str(unit["family"]))]
            for unit in primary_units
            if int(unit.get("n_trees") or 0) > 1
        }
        if category == "incomplete":
            selection_class = "incomplete"
            selection_note = "candidate set incomplete; read-AF selection not evaluated"
        elif category == "exact_and_topology_unique":
            selection_class = "structural_exact_unique"
            selection_note = "C=1 and Topo=1; no read-AF ranking required"
        else:
            ambiguous = [
                ranking_by_unit[(region_id, str(unit["family"]))]
                for unit in primary_units
                if int(unit["n_trees"]) > 1
            ]
            if not ambiguous or not all(item.get("status") == "ranked" for item in ambiguous):
                selection_class = "read_af_unavailable"
                selection_note = "; ".join(
                    sorted({item.get("status", "unknown") for item in ambiguous})
                )
            else:
                n_top = math.prod(int(item["n_top_exact"]) for item in ambiguous)
                n_top_shapes = math.prod(int(item["n_top_shapes_exact"]) for item in ambiguous)
                base["n_top_joint_exact"] = n_top
                base["n_top_shapes_joint_exact"] = n_top_shapes
                if n_top == 1:
                    selection_class = "read_af_unique_first"
                    selection_note = "one exact first-ranked candidate combination under read-AF ordering"
                elif n_top_shapes == 1:
                    selection_class = "read_af_tied_same_topology"
                    selection_note = "multiple exact first-ranked candidates share one topology shape"
                else:
                    selection_class = "read_af_tied_different_topology"
                    selection_note = "exact first-rank tie spans multiple topology shapes"

        has_direct = False
        has_sister = False
        unresolved_reasons = []
        for unit in primary_units:
            family = str(unit["family"])
            representative, source = representative_for_morphology(
                unit, ranking_by_unit.get((region_id, family))
            )
            if representative is None:
                unresolved_reasons.append(f"HP{family}:{source}")
                continue
            features = representative["morphology_features"]
            has_direct = has_direct or bool(features["has_direct_chain"])
            has_sister = has_sister or bool(features["has_sister_branch"])
        if unresolved_reasons:
            morphology = "unresolved"
            morphology_note = "; ".join(unresolved_reasons)
        else:
            morphology = morphology_class(has_direct, has_sister)
            morphology_note = (
                "OR across independent primary HP components; no cross-HP edge is inferred"
            )
        for ranking in unit_payload.values():
            if ranking.get("status") != "ranked":
                continue
            if any(
                int(candidate.get("candidate_index") or 0) > 32
                for candidate in ranking.get("top_shape_representatives", [])
            ):
                top_edge_outside_display_cap += 1

        selection_counts[selection_class] += 1
        morphology_counts[morphology] += 1
        regions_output.append(
            {
                **base,
                "selection_class": selection_class,
                "selection_note": selection_note,
                "morphology_class": morphology,
                "morphology_note": morphology_note,
                "units": unit_payload,
            }
        )

    expected_structural = sample_record["topology_classes"]
    checks = {
        "region_count_matches_W_tree": len(regions_output) == int(sample_record["W_tree"]),
        "no_primary_matches_current": no_primary == int(sample_record["no_primary_lineage"]),
        "primary_count_matches_W_primary": len(regions_output) - no_primary
        == int(sample_record["W_primary"]),
        "structural_classes_match_current": all(
            structural_counts[key] == int(expected_structural[key]) for key in STRUCTURAL_KEYS
        ),
        "selection_partition_conserves_W_primary": sum(selection_counts.values())
        == int(sample_record["W_primary"]),
        "morphology_partition_conserves_W_primary": sum(morphology_counts.values())
        == int(sample_record["W_primary"]),
        "no_candidate_or_shape_mismatch": unit_status["candidate_or_shape_mismatch"] == 0,
        "ranked_units_have_top_edges": all(
            ranking.get("status") != "ranked"
            or bool(ranking.get("top_shape_representatives"))
            for ranking in ranking_by_unit.values()
        ),
    }
    if not all(checks.values()):
        fail(f"{sample}: validation checks failed: {checks}")

    output = {
        "schema_name": "intersubmod.current_v5_read_af_topology_sample",
        "schema_version": "1.0.0",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "sample": sample,
        "scope": "GRCh38 chr1-22 current canonical v5",
        "claim_ceiling": (
            "descriptive family-specific read-AF first-rank and regional mutation-state "
            "topology morphology; not confirmed clone, true ancestry, CCF, or calibrated likelihood"
        ),
        "method": {
            "read_af": "family-specific ALT/(REF+ALT) from col_coverage_by_hp",
            "candidate_enumeration": "all candidates with tree_cap=0; fail closed on count/shape mismatch",
            "tie": "exact equality of Fraction-valued ordering scores",
            "morphology": (
                "Topo=1 or read-AF co-top shape-set size=1; direct=max_depth>=2, "
                "sister=max_outdegree>=2; region is OR across independent primary HP components"
            ),
            "prohibited_interpretation": "score is not probability/posterior and morphology is not a clone census",
        },
        "summary": {
            "W_tree": int(sample_record["W_tree"]),
            "W_primary": int(sample_record["W_primary"]),
            "no_primary": no_primary,
            "structural_classes": {key: structural_counts[key] for key in STRUCTURAL_KEYS},
            "selection_classes": {key: selection_counts[key] for key in SELECTION_KEYS},
            "morphology_classes": {key: morphology_counts[key] for key in MORPHOLOGY_KEYS},
            "unit_status": dict(sorted(unit_status.items())),
            "top_shape_representative_outside_original_display_cap_32": top_edge_outside_display_cap,
            "checks": checks,
            "all_checks_pass": all(checks.values()),
        },
        "regions": regions_output,
        "provenance": {
            "run_root": str(run_root),
            "input_manifest": str(manifest_path),
            "input_manifest_sha256": file_sha256(manifest_path),
            "current_summary": str(summary_path),
            "current_summary_sha256": file_sha256(summary_path),
            "layered_reconstruction": str(layered_path),
            "layered_reconstruction_sha256": file_sha256(layered_path),
            "layered_region_view": str(region_view_path),
            "layered_region_view_sha256": file_sha256(region_view_path),
        },
    }
    output_path = output_dir / f"{sample}.current_v5_read_af_topology.json"
    write_json(output_path, output, compact=True)
    return {
        "sample": sample,
        "output": str(output_path.resolve()),
        "output_sha256": file_sha256(output_path),
        **output["summary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--current-summary", required=True, type=Path)
    parser.add_argument("--method-script-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--samples", default=",".join(SAMPLE_ORDER))
    args = parser.parse_args()

    run_root = args.run_root.resolve()
    manifest_path = args.input_manifest.resolve()
    summary_path = args.current_summary.resolve()
    method_script_dir = args.method_script_dir.resolve()
    output_dir = args.output_dir.resolve()
    for path in (run_root, manifest_path, summary_path, method_script_dir):
        if not path.exists():
            fail(f"missing input: {path}")
    if str(method_script_dir) not in sys.path:
        sys.path.insert(0, str(method_script_dir))
    solver = importlib.import_module("tree_enumeration_solver")

    manifest = read_json(manifest_path)
    current = read_json(summary_path)
    if current.get("all_pass") is not True:
        fail("current summary all_pass is not true")
    canonical = current.get("canonical") or {}
    if canonical.get("label") != "longphase_s_recalibrated_FILTER_PASS":
        fail(f"unexpected canonical label: {canonical.get('label')}")
    if Path(canonical.get("run_root", "")).resolve() != run_root:
        fail("current summary canonical run_root differs from requested run root")
    manifest_samples = {row["sample"]: row for row in manifest.get("samples", [])}
    current_samples = {row["sample"]: row for row in canonical.get("samples", [])}
    requested = [value.strip() for value in args.samples.split(",") if value.strip()]
    unknown = [sample for sample in requested if sample not in manifest_samples or sample not in current_samples]
    if unknown:
        fail(f"unknown requested samples: {unknown}")

    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for sample in requested:
        summary = analyze_sample(
            sample,
            run_root,
            current_samples[sample],
            solver,
            manifest_path,
            summary_path,
            output_dir,
        )
        summaries.append(summary)
        print(
            f"{sample}: W_primary={summary['W_primary']} "
            f"ranked_units={summary['unit_status'].get('ranked', 0)} "
            f"missing_read_af={summary['unit_status'].get('missing_read_af', 0)} "
            f"outside_cap32={summary['top_shape_representative_outside_original_display_cap_32']}"
        )

    aggregate = {
        "schema_name": "intersubmod.current_v5_read_af_topology_index",
        "schema_version": "1.0.0",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "scope": "GRCh38 chr1-22 current canonical v5",
        "samples": summaries,
        "dataset_count": len(summaries),
        "all_checks_pass": all(sample["all_checks_pass"] for sample in summaries),
        "provenance": {
            "run_root": str(run_root),
            "input_manifest": str(manifest_path),
            "input_manifest_sha256": file_sha256(manifest_path),
            "current_summary": str(summary_path),
            "current_summary_sha256": file_sha256(summary_path),
            "solver": str((method_script_dir / "tree_enumeration_solver.py").resolve()),
            "solver_sha256": file_sha256(method_script_dir / "tree_enumeration_solver.py"),
            "builder": str(Path(__file__).resolve()),
            "builder_sha256": file_sha256(Path(__file__).resolve()),
        },
    }
    index_path = output_dir / "current_v5_read_af_topology.index.json"
    write_json(index_path, aggregate)
    print(f"CURRENT-V5 READ-AF TOPOLOGY -> {index_path} ({len(summaries)} datasets)")


if __name__ == "__main__":
    main()
