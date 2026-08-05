#!/usr/bin/env python3
"""Postprocess the HCC1395 site-topology comparison into report tables.

This script is intentionally independent from the tree-enumeration builder. It
only consumes the frozen region table, the machine-readable component
signature sets, and the frozen input manifest/VCFs.  It produces:

* a per-shared-site-pair evidence ledger under one read-selected whole-region
  HP mapping;
* one mutually exclusive, fixed-denominator outcome per region and layer;
* complexity strata and deterministic representative examples; and
* an independent CHROM/POS/REF/ALT identity audit from the two source VCFs.

The output describes technical reproducibility of read-compatible candidate
spaces.  It is not biological clone-tree truth.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import math
from itertools import combinations
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pysam


SAMPLE_A = "HCC1395"
SAMPLE_B = "HCC1395_DORADO"
LAYERS = ("read_full", "vaf_official", "vaf_normalized")
OUTCOMES = (
    "strict_full_exact",
    "A_induced_substructure",
    "B_induced_substructure",
    "resolution_A_more_specific",
    "resolution_B_more_specific",
    "candidate_overlap",
    "conflict",
    "shared_core_only",
    "not_evaluable",
)
PAIR_OUTCOMES = (
    "determined_same",
    "ambiguous_equal",
    "one_sided_contained",
    "ambiguous_overlap",
    "conflict",
    "not_evaluable",
)
SITE_RELATIONS = (
    "equal",
    "a_proper_subset_b",
    "b_proper_subset_a",
    "partial_overlap",
    "disjoint",
)
FIXED_DENOMINATOR = 5720


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--signatures", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_jsonl_gz(path: Path) -> list[dict]:
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def write_tsv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def deterministic_gzip_text(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="")

    class Writer:
        def __enter__(self):
            return text

        def __exit__(self, exc_type, exc, traceback):
            try:
                text.close()
            finally:
                if not raw.closed:
                    raw.close()
            return False

    return Writer()


def write_tsv_gz(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    with deterministic_gzip_text(path) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start_text, end_text = span.split("-", 1)
    return chrom, int(start_text), int(end_text)


def chromosome_key(chrom: str) -> tuple[int, str]:
    value = chrom.removeprefix("chr")
    return (int(value), "") if value.isdigit() else (10_000, value)


def region_key(row: Mapping[str, object]) -> tuple:
    chrom, start, end = parse_region(str(row["region"]))
    return chromosome_key(chrom), start, end, str(row.get("match_id", ""))


def bool_text(value: object) -> bool:
    return value is True or value == "True"


def set_relation(a_values: Iterable[object], b_values: Iterable[object]) -> str:
    a, b = set(a_values), set(b_values)
    if a == b:
        return "equal"
    if a < b:
        return "a_proper_subset_b"
    if b < a:
        return "b_proper_subset_a"
    if a & b:
        return "partial_overlap"
    return "disjoint"


def load_vcf_paths(manifest_path: Path) -> dict[str, Path]:
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    result = {
        row["sample"]: Path(row["somatic_vcf"])
        for row in document["samples"]
        if row["sample"] in {SAMPLE_A, SAMPLE_B}
    }
    if set(result) != {SAMPLE_A, SAMPLE_B}:
        raise RuntimeError(f"manifest does not contain both technical datasets: {sorted(result)}")
    for sample, path in result.items():
        if not path.is_file():
            raise RuntimeError(f"VCF missing for {sample}: {path}")
    return result


def variants_in_region(vcf: pysam.VariantFile, chrom: str, start: int, end: int) -> dict[int, set[tuple]]:
    result: dict[int, set[tuple]] = collections.defaultdict(set)
    for record in vcf.fetch(chrom, start - 1, end):
        if not start <= record.pos <= end:
            continue
        for alt in record.alts or ():
            result[int(record.pos)].add((str(record.chrom), int(record.pos), str(record.ref), str(alt)))
    return dict(result)


def audit_vcf_sites(
    region_rows: Sequence[dict[str, str]], vcf_paths: Mapping[str, Path]
) -> tuple[dict[str, dict], list[dict], dict]:
    vcfs = {sample: pysam.VariantFile(str(path)) for sample, path in vcf_paths.items()}
    by_match: dict[str, dict] = {}
    locus_rows: list[dict] = []
    counters: collections.Counter = collections.Counter()
    try:
        for row in sorted(region_rows, key=region_key):
            chrom, start, end = parse_region(row["region"])
            variants_a = variants_in_region(vcfs[SAMPLE_A], chrom, start, end)
            variants_b = variants_in_region(vcfs[SAMPLE_B], chrom, start, end)
            positions_a, positions_b = set(variants_a), set(variants_b)
            shared = sorted(positions_a & positions_b)
            relation = set_relation(positions_a, positions_b)
            counters[f"region_site_relation|{relation}"] += 1
            counters["caller_sites_a"] += len(positions_a)
            counters["caller_sites_b"] += len(positions_b)
            counters["caller_shared_sites"] += len(shared)
            mismatch_positions = []
            for position in shared:
                same = variants_a[position] == variants_b[position]
                counters["caller_shared_allele_identity"] += same
                counters["caller_shared_allele_mismatch"] += not same
                if not same:
                    mismatch_positions.append(position)
                locus_rows.append(
                    {
                        "match_id": row["match_id"],
                        "chrom": chrom,
                        "region": row["region"],
                        "position": position,
                        "alleles_a": canonical_json(sorted(variants_a[position])),
                        "alleles_b": canonical_json(sorted(variants_b[position])),
                        "allele_identity": same,
                    }
                )
            by_match[row["match_id"]] = {
                "caller_k_a": len(positions_a),
                "caller_k_b": len(positions_b),
                "caller_shared_k": len(shared),
                "caller_site_set_relation": relation,
                "shared_positions": shared,
                "allele_identity_count": len(shared) - len(mismatch_positions),
                "allele_mismatch_count": len(mismatch_positions),
                "allele_mismatch_positions": mismatch_positions,
                "alleles_a": variants_a,
                "alleles_b": variants_b,
            }
    finally:
        for vcf in vcfs.values():
            vcf.close()
    return by_match, locus_rows, dict(counters)


def q_relation_set(projection: Mapping[str, object], k: int, pair_index: int) -> tuple[str, ...] | None:
    if projection.get("status") not in {"ok", "noninformative"}:
        return None
    expected_length = math.comb(k, 2)
    values = set()
    for signature in projection.get("q", []):
        if len(signature) != expected_length:
            raise RuntimeError(f"signature length {len(signature)} != C({k},2)={expected_length}")
        values.add(signature[pair_index])
    if not values:
        raise RuntimeError("evaluable projection has an empty relation set")
    return tuple(sorted(values))


def relation_labels(codes: Sequence[str], left: int, right: int) -> list[str]:
    labels = {"F": f"{left}>{right}", "R": f"{right}>{left}", "P": "parallel"}
    return [labels[code] for code in codes]


def pair_set_outcome(a: tuple[str, ...] | None, b: tuple[str, ...] | None) -> tuple[str, str]:
    if a is None or b is None:
        return "not_evaluable", "endpoint_projection_not_evaluable"
    relation = set_relation(a, b)
    if relation == "equal":
        if len(a) == 1:
            return "determined_same", "both_relation_sets_are_the_same_singleton"
        return "ambiguous_equal", "equal_but_multi_valued_candidate_relation_sets"
    if relation in {"a_proper_subset_b", "b_proper_subset_a"}:
        return "one_sided_contained", relation
    if relation == "partial_overlap":
        return "ambiguous_overlap", "candidate_relation_sets_overlap"
    return "conflict", "candidate_relation_sets_disjoint"


def evidence_rows(signature_rows: Sequence[dict]) -> list[dict]:
    rows: list[dict] = []
    for record in sorted(signature_rows, key=region_key):
        read_selected = record["swap_tolerant"]["read_full"].get("selected_mapping", "")
        vaf_selected = record["swap_tolerant"]["vaf_official"].get("selected_mapping", "")
        if not read_selected:
            continue
        alignment = record["alignments"][read_selected]
        for component_index, component in enumerate(alignment["components"], start=1):
            shared = [int(value) for value in component["shared_positions"]]
            pair_order = list(combinations(range(len(shared)), 2))
            for pair_index, (left_index, right_index) in enumerate(pair_order):
                left, right = shared[left_index], shared[right_index]
                output = {
                    "match_id": record["match_id"],
                    "chrom": record["chrom"],
                    "region": record["region"],
                    "read_selected_mapping": read_selected,
                    "vaf_selected_mapping": vaf_selected,
                    "read_vaf_mapping_same": read_selected == vaf_selected,
                    "component_index": component_index,
                    "family_a": component["a_family"],
                    "family_b": component["b_family"],
                    "component_site_set_relation": component["site_set_relation"],
                    "component_shared_k": len(shared),
                    "pair_index": pair_index,
                    "left_position": left,
                    "right_position": right,
                }
                for layer in LAYERS:
                    projection_a = component[layer]["a"]
                    projection_b = component[layer]["b"]
                    values_a = q_relation_set(projection_a, len(shared), pair_index)
                    values_b = q_relation_set(projection_b, len(shared), pair_index)
                    outcome, reason = pair_set_outcome(values_a, values_b)
                    output[f"{layer}_relation_codes_a"] = "" if values_a is None else canonical_json(values_a)
                    output[f"{layer}_relation_codes_b"] = "" if values_b is None else canonical_json(values_b)
                    output[f"{layer}_relations_a"] = (
                        "" if values_a is None else canonical_json(relation_labels(values_a, left, right))
                    )
                    output[f"{layer}_relations_b"] = (
                        "" if values_b is None else canonical_json(relation_labels(values_b, left, right))
                    )
                    output[f"{layer}_outcome"] = outcome
                    output[f"{layer}_reason"] = reason
                    output[f"{layer}_determined_same"] = (
                        values_a is not None and values_b is not None and values_a == values_b and len(values_a) == 1
                    )
                rows.append(output)
    return rows


def selected_alignment(record: Mapping[str, object], layer: str) -> tuple[dict | None, dict]:
    summary = record["swap_tolerant"][layer]
    mapping = summary.get("selected_mapping", "")
    if not mapping:
        return None, summary
    return record["alignments"][mapping], summary


def fixed_outcome(record: Mapping[str, object], layer: str, caller_relation: str) -> tuple[str, str]:
    alignment, summary = selected_alignment(record, layer)
    if alignment is None or summary.get("status") != "evaluable":
        return "not_evaluable", str(summary.get("reason", "mapping_or_endpoint_not_evaluable"))
    components = alignment["components"]
    candidate_relation = summary["category"]
    # A disjoint forest product must come from at least one informative
    # component: k<2 components project to the same neutral empty relation
    # vector and cannot create disjointness.  Preserve that hard negative even
    # when a different component is noninformative.
    if candidate_relation == "disjoint":
        return "conflict", "projected_candidate_topology_sets_disjoint"
    if not components or any(int(component["shared_k"]) < 2 for component in components):
        return "shared_core_only", "at_least_one_component_has_fewer_than_two_shared_sites"

    site_relations = [caller_relation] + [component["site_set_relation"] for component in components]
    if any(value in {"partial_overlap", "disjoint"} for value in site_relations):
        return "shared_core_only", "non_nested_site_inventory"
    has_a_direction = "a_proper_subset_b" in site_relations
    has_b_direction = "b_proper_subset_a" in site_relations
    if has_a_direction and has_b_direction:
        return "shared_core_only", "site_containment_direction_mixed_across_caller_or_HP_components"

    if not has_a_direction and not has_b_direction:
        if candidate_relation == "exact":
            return "strict_full_exact", "equal_site_inventories_and_equal_projected_candidate_sets"
        if candidate_relation == "a_proper_subset_b":
            return "resolution_A_more_specific", "equal_sites_but_A_candidate_space_is_narrower"
        if candidate_relation == "b_proper_subset_a":
            return "resolution_B_more_specific", "equal_sites_but_B_candidate_space_is_narrower"
        if candidate_relation == "overlap":
            return "candidate_overlap", "equal_sites_with_partially_overlapping_candidate_spaces"
        raise RuntimeError(f"unknown candidate relation: {candidate_relation}")

    if has_a_direction:
        if candidate_relation in {"exact", "a_proper_subset_b"}:
            return "A_induced_substructure", "A_sites_nested_in_B_and_candidate_relation_is_equal_or_same_direction"
        if candidate_relation == "overlap":
            return "candidate_overlap", "A_sites_nested_in_B_but_candidate_spaces_only_overlap"
        return "shared_core_only", "A_site_nesting_and_candidate_set_direction_disagree"

    if candidate_relation in {"exact", "b_proper_subset_a"}:
        return "B_induced_substructure", "B_sites_nested_in_A_and_candidate_relation_is_equal_or_same_direction"
    if candidate_relation == "overlap":
        return "candidate_overlap", "B_sites_nested_in_A_but_candidate_spaces_only_overlap"
    return "shared_core_only", "B_site_nesting_and_candidate_set_direction_disagree"


def k_bin(value: int) -> str:
    if value < 2:
        return "k<2"
    if value == 2:
        return "k=2"
    if value == 3:
        return "k=3"
    if value == 4:
        return "k=4"
    return "k>=5"


def pair_outcome_rows(
    region_rows: Sequence[dict[str, str]], signature_by_match: Mapping[str, dict], vcf_audit: Mapping[str, dict]
) -> list[dict]:
    rows = []
    for source in sorted(region_rows, key=region_key):
        match_id = source["match_id"]
        signature = signature_by_match[match_id]
        vcf = vcf_audit[match_id]
        output = {
            "match_id": match_id,
            "chrom": source["chrom"],
            "region": source["region"],
            "fixed_denominator": FIXED_DENOMINATOR,
            "caller_k_a": vcf["caller_k_a"],
            "caller_k_b": vcf["caller_k_b"],
            "caller_shared_k": vcf["caller_shared_k"],
            "k_bin": k_bin(vcf["caller_shared_k"]),
            "caller_site_set_relation": vcf["caller_site_set_relation"],
            "shared_allele_identity_count": vcf["allele_identity_count"],
            "shared_allele_mismatch_count": vcf["allele_mismatch_count"],
            "allele_identity_pass": vcf["allele_mismatch_count"] == 0,
            "hp_count_a": source["hp_count_a"],
            "hp_count_b": source["hp_count_b"],
            "hp_families_a": source["hp_families_a"],
            "hp_families_b": source["hp_families_b"],
            "hp_count_mismatch": source["hp_count_mismatch"],
            "hp_specific_site_set_relation": source["site_set_relation"],
            "hp_specific_shared_k": source["shared_k"],
        }
        selected_mappings = {}
        for layer in LAYERS:
            alignment, summary = selected_alignment(signature, layer)
            outcome, reason = fixed_outcome(signature, layer, vcf["caller_site_set_relation"])
            selected_mapping = summary.get("selected_mapping", "")
            selected_mappings[layer] = selected_mapping
            output[f"{layer}_selected_mapping"] = selected_mapping
            output[f"{layer}_evaluable"] = summary.get("status") == "evaluable"
            output[f"{layer}_reason"] = summary.get("reason", "") or reason
            output[f"{layer}_candidate_relation"] = summary.get("category", "")
            output[f"{layer}_candidate_tree_product_a"] = summary.get("candidate_tree_product_a", "")
            output[f"{layer}_candidate_tree_product_b"] = summary.get("candidate_tree_product_b", "")
            output[f"{layer}_topology_set_size_a"] = summary.get("topology_set_size_a", "")
            output[f"{layer}_topology_set_size_b"] = summary.get("topology_set_size_b", "")
            output[f"{layer}_jaccard"] = summary.get("jaccard", "")
            output[f"{layer}_outcome"] = outcome
            output[f"{layer}_outcome_reason"] = reason
        output["read_vaf_selected_mapping_same"] = (
            bool(selected_mappings["read_full"])
            and selected_mappings["read_full"] == selected_mappings["vaf_official"]
        )
        rows.append(output)
    return rows


def fixed_metrics(rows: Sequence[dict]) -> list[dict]:
    output = []
    for layer in LAYERS:
        counts = collections.Counter(row[f"{layer}_outcome"] for row in rows)
        for outcome in OUTCOMES:
            n = counts[outcome]
            output.append(
                {
                    "layer": layer,
                    "mapping": "whole_region_HP_swap_tolerant",
                    "outcome": outcome,
                    "n": n,
                    "denominator": FIXED_DENOMINATOR,
                    "share": f"{100 * n / FIXED_DENOMINATOR:.6f}",
                }
            )
    return output


def complexity_strata(rows: Sequence[dict]) -> list[dict]:
    output = []
    bins = ("k<2", "k=2", "k=3", "k=4", "k>=5")
    for layer in LAYERS:
        for label in bins:
            selected = [row for row in rows if row["k_bin"] == label]
            counts = collections.Counter(row[f"{layer}_outcome"] for row in selected)
            for outcome in OUTCOMES:
                n = counts[outcome]
                output.append(
                    {
                        "layer": layer,
                        "mapping": "whole_region_HP_swap_tolerant",
                        "k_bin": label,
                        "outcome": outcome,
                        "n": n,
                        "denominator": len(selected),
                        "share": f"{100 * n / len(selected):.6f}" if selected else "",
                    }
                )
    return output


def deterministic_examples(rows: Sequence[dict], per_outcome: int = 2) -> list[dict]:
    output = []
    ordered = sorted(rows, key=region_key)
    for layer in LAYERS:
        for outcome in OUTCOMES:
            selected = [row for row in ordered if row[f"{layer}_outcome"] == outcome][:per_outcome]
            for rank, row in enumerate(selected, start=1):
                output.append(
                    {
                        "layer": layer,
                        "outcome": outcome,
                        "example_rank": rank,
                        "match_id": row["match_id"],
                        "chrom": row["chrom"],
                        "region": row["region"],
                        "caller_site_set_relation": row["caller_site_set_relation"],
                        "caller_shared_k": row["caller_shared_k"],
                        "hp_count_a": row["hp_count_a"],
                        "hp_count_b": row["hp_count_b"],
                        "selected_mapping": row[f"{layer}_selected_mapping"],
                        "candidate_relation": row[f"{layer}_candidate_relation"],
                        "candidate_tree_product_a": row[f"{layer}_candidate_tree_product_a"],
                        "candidate_tree_product_b": row[f"{layer}_candidate_tree_product_b"],
                        "topology_set_size_a": row[f"{layer}_topology_set_size_a"],
                        "topology_set_size_b": row[f"{layer}_topology_set_size_b"],
                        "outcome_reason": row[f"{layer}_outcome_reason"],
                    }
                )
    return output


def main() -> None:
    args = parse_args()
    region_rows = read_tsv(args.regions)
    signature_rows = read_jsonl_gz(args.signatures)
    if len(region_rows) != FIXED_DENOMINATOR or len(signature_rows) != FIXED_DENOMINATOR:
        raise RuntimeError(
            f"fixed denominator mismatch: regions={len(region_rows)}, signatures={len(signature_rows)}"
        )
    signature_by_match = {row["match_id"]: row for row in signature_rows}
    if len(signature_by_match) != FIXED_DENOMINATOR:
        raise RuntimeError("duplicate match_id in signature JSONL")
    if {row["match_id"] for row in region_rows} != set(signature_by_match):
        raise RuntimeError("region TSV and signature JSONL match_id sets differ")

    vcf_paths = load_vcf_paths(args.input_manifest)
    vcf_audit, allele_rows, vcf_counts = audit_vcf_sites(region_rows, vcf_paths)
    outcomes = pair_outcome_rows(region_rows, signature_by_match, vcf_audit)
    evidence = evidence_rows(signature_rows)
    metrics = fixed_metrics(outcomes)
    strata = complexity_strata(outcomes)
    examples = deterministic_examples(outcomes)

    checks: list[dict] = []

    def check(name: str, observed, expected) -> None:
        passed = observed == expected
        checks.append({"check": name, "pass": passed, "observed": observed, "expected": expected})
        if not passed:
            raise RuntimeError(f"QA failed {name}: {observed!r} != {expected!r}")

    check("region_input_fixed_denominator", len(region_rows), FIXED_DENOMINATOR)
    check("signature_input_fixed_denominator", len(signature_rows), FIXED_DENOMINATOR)
    check("pair_outcome_fixed_denominator", len(outcomes), FIXED_DENOMINATOR)
    check("match_id_join_conservation", len(signature_by_match), FIXED_DENOMINATOR)
    check("caller_site_relation_conservation", sum(vcf_counts.get(f"region_site_relation|{v}", 0) for v in SITE_RELATIONS), FIXED_DENOMINATOR)
    check("caller_equal_site_regions", vcf_counts.get("region_site_relation|equal", 0), 5534)
    check(
        "caller_strict_subset_site_regions",
        vcf_counts.get("region_site_relation|a_proper_subset_b", 0)
        + vcf_counts.get("region_site_relation|b_proper_subset_a", 0),
        185,
    )
    check("caller_partial_overlap_site_regions", vcf_counts.get("region_site_relation|partial_overlap", 0), 1)
    check("caller_disjoint_site_regions", vcf_counts.get("region_site_relation|disjoint", 0), 0)
    check("caller_shared_site_loci", vcf_counts.get("caller_shared_sites", 0), 15713)
    check("caller_shared_allele_identity", vcf_counts.get("caller_shared_allele_identity", 0), 15713)
    check("caller_shared_allele_mismatch", vcf_counts.get("caller_shared_allele_mismatch", 0), 0)
    check("allele_ledger_row_conservation", len(allele_rows), vcf_counts.get("caller_shared_sites", 0))
    check("evidence_pair_outcome_domain", sum(row[f"{layer}_outcome"] in PAIR_OUTCOMES for row in evidence for layer in LAYERS), len(evidence) * len(LAYERS))
    for layer in LAYERS:
        check(
            f"{layer}_fixed_outcome_conservation",
            sum(row[f"{layer}_outcome"] in OUTCOMES for row in outcomes),
            FIXED_DENOMINATOR,
        )
        check(
            f"{layer}_fixed_metric_conservation",
            sum(int(row["n"]) for row in metrics if row["layer"] == layer),
            FIXED_DENOMINATOR,
        )
        check(
            f"{layer}_complexity_stratum_conservation",
            sum(int(row["n"]) for row in strata if row["layer"] == layer),
            FIXED_DENOMINATOR,
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    outcomes_path = args.output_dir / "hcc1395_site_topology_pair_outcomes.tsv"
    evidence_path = args.output_dir / "hcc1395_site_relation_evidence.tsv.gz"
    allele_path = args.output_dir / "hcc1395_site_allele_identity.tsv.gz"
    metrics_path = args.output_dir / "hcc1395_topology_compatibility_metrics.tsv"
    strata_path = args.output_dir / "hcc1395_topology_outcome_complexity_strata.tsv"
    examples_path = args.output_dir / "hcc1395_topology_examples.tsv"
    summary_path = args.output_dir / "hcc1395_topology_containment_summary.json"
    checks_path = args.output_dir / "hcc1395_topology_containment_checks.tsv"

    write_tsv(outcomes_path, outcomes, list(outcomes[0]))
    write_tsv_gz(evidence_path, evidence, list(evidence[0]))
    write_tsv_gz(
        allele_path,
        allele_rows,
        ["match_id", "chrom", "region", "position", "alleles_a", "alleles_b", "allele_identity"],
    )
    write_tsv(metrics_path, metrics, ["layer", "mapping", "outcome", "n", "denominator", "share"])
    write_tsv(
        strata_path,
        strata,
        ["layer", "mapping", "k_bin", "outcome", "n", "denominator", "share"],
    )
    write_tsv(
        examples_path,
        examples,
        [
            "layer",
            "outcome",
            "example_rank",
            "match_id",
            "chrom",
            "region",
            "caller_site_set_relation",
            "caller_shared_k",
            "hp_count_a",
            "hp_count_b",
            "selected_mapping",
            "candidate_relation",
            "candidate_tree_product_a",
            "candidate_tree_product_b",
            "topology_set_size_a",
            "topology_set_size_b",
            "outcome_reason",
        ],
    )
    write_tsv(checks_path, checks, ["check", "pass", "observed", "expected"])

    outcome_summary = {
        layer: {
            outcome: sum(row[f"{layer}_outcome"] == outcome for row in outcomes)
            for outcome in OUTCOMES
        }
        for layer in LAYERS
    }
    site_pair_summary = {
        layer: {
            outcome: sum(row[f"{layer}_outcome"] == outcome for row in evidence)
            for outcome in PAIR_OUTCOMES
        }
        for layer in LAYERS
    }
    summary = {
        "schema_version": "1.0",
        "analysis_date": "2026-07-12",
        "fixed_denominator": FIXED_DENOMINATOR,
        "scope": "HCC1395 vs HCC1395_DORADO; chr1-22; exact-coordinate complete-both regions",
        "claim_ceiling": (
            "Same-cell-line cross-source technical reproducibility of shared-site mutation constraints only; "
            "not biological clone truth, accuracy, or independent validation."
        ),
        "outcome_precedence": list(OUTCOMES),
        "outcome_definitions": {
            "strict_full_exact": "Caller and every mapped HP site inventory are equal; projected candidate sets are equal.",
            "A_induced_substructure": "A site inventories are consistently nested in B; every component has shared_k>=2; candidate relation is equal or A-subset-B.",
            "B_induced_substructure": "B site inventories are consistently nested in A; every component has shared_k>=2; candidate relation is equal or B-subset-A.",
            "resolution_A_more_specific": "Equal site inventories; A has the narrower projected candidate space.",
            "resolution_B_more_specific": "Equal site inventories; B has the narrower projected candidate space.",
            "candidate_overlap": "Projected candidate spaces overlap but neither contains the other.",
            "conflict": "Projected candidate spaces are disjoint on the shared sites.",
            "shared_core_only": "Some shared relation exists, but non-nested/mixed/private or k<2 component scope blocks a full-region claim.",
            "not_evaluable": "HP count/mapping, fewer than two shared relations, recurrence, or VAF availability blocks evaluation.",
        },
        "mapping_contract": {
            "fixed_outcomes": "Each layer uses its own precomputed identity-or-one-global-swap selected mapping.",
            "site_pair_evidence": "Read-full selected mapping is frozen first; VAF is displayed under that same mapping to avoid per-site cherry-picking.",
        },
        "site_inventory": {
            "caller_region_relations": {
                relation: vcf_counts.get(f"region_site_relation|{relation}", 0) for relation in SITE_RELATIONS
            },
            "caller_shared_sites": vcf_counts.get("caller_shared_sites", 0),
            "caller_shared_allele_identity": vcf_counts.get("caller_shared_allele_identity", 0),
            "caller_shared_allele_mismatch": vcf_counts.get("caller_shared_allele_mismatch", 0),
        },
        "layers": outcome_summary,
        "site_pair_evidence": site_pair_summary,
        "evidence_rows": len(evidence),
        "examples_per_nonempty_layer_outcome": 2,
        "sources": [
            {"path": str(args.regions.resolve()), "sha256": sha256(args.regions)},
            {"path": str(args.signatures.resolve()), "sha256": sha256(args.signatures)},
            {"path": str(args.input_manifest.resolve()), "sha256": sha256(args.input_manifest)},
            *[
                {"sample": sample, "path": str(path.resolve()), "sha256": sha256(path)}
                for sample, path in sorted(vcf_paths.items())
            ],
        ],
        "outputs": {
            "pair_outcomes": str(outcomes_path.resolve()),
            "site_relation_evidence": str(evidence_path.resolve()),
            "allele_identity": str(allele_path.resolve()),
            "compatibility_metrics": str(metrics_path.resolve()),
            "outcome_complexity_strata": str(strata_path.resolve()),
            "examples": str(examples_path.resolve()),
            "checks": str(checks_path.resolve()),
        },
        "checks": checks,
    }
    write_json(summary_path, summary)

    print(f"INPUT regions: {args.regions.resolve()}")
    print(f"INPUT signatures: {args.signatures.resolve()}")
    print(f"INPUT manifest: {args.input_manifest.resolve()}")
    print(f"OUTPUT pair outcomes: {outcomes_path.resolve()}")
    print(f"OUTPUT site-pair evidence: {evidence_path.resolve()}")
    print(f"OUTPUT allele identity: {allele_path.resolve()}")
    print(f"OUTPUT metrics: {metrics_path.resolve()}")
    print(f"OUTPUT complexity strata: {strata_path.resolve()}")
    print(f"OUTPUT examples: {examples_path.resolve()}")
    print(f"OUTPUT checks: {checks_path.resolve()}")
    print(f"OUTPUT summary: {summary_path.resolve()}")
    print(
        f"RESULT regions={len(outcomes)} evidence_pairs={len(evidence)} "
        f"shared_alleles={vcf_counts.get('caller_shared_allele_identity', 0)}/"
        f"{vcf_counts.get('caller_shared_sites', 0)} checks={len(checks)}/{len(checks)} PASS"
    )


if __name__ == "__main__":
    main()
