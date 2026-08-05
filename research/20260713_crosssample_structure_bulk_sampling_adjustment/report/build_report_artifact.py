#!/usr/bin/env python3
"""Build the canonical portable-report artifact for cross-sample structure validation.

The report is deliberately fail-closed.  It verifies the current 7-sample
historical layered-v2 census, the bulk-sampling adjustment snapshot, and the
regional PyClone-VI annotation snapshot before authoring ``artifact.json``.
It does not render HTML; the shared Data Analytics portable renderer remains
the sole renderer.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCRIPT_PATH = Path(__file__).resolve()
TOPIC_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]

RERUN_ROOT = TOPIC_ROOT / "results" / "rerun_a"
REGION_ROOT = TOPIC_ROOT / "results" / "region_possible_clone_annotations_v1"
CENSUS_ROOT = REPO_ROOT / "research" / "20260712_vaf_selected_shape_four_class_census" / "data"
BRIDGE_ROOT = (
    REPO_ROOT
    / "research"
    / "20260712_vaf_ccf_subclone_crosssoftware_validation"
    / "results"
    / "clone_region_bridge_v1"
)
TOPOLOGY_PAIR_ROOT = REPO_ROOT / "research" / "20260712_hcc1395_pair_coarse_topology_gene_drug_validation"
TOPOLOGY_PAIR_JSON = TOPOLOGY_PAIR_ROOT / "data" / "topology_pair_analysis.json"
TOPOLOGY_PAIR_SCRIPT = TOPOLOGY_PAIR_ROOT / "scripts" / "topology_pair_analysis.py"

RERUN_REQUIRED = {
    "summary": "summary.json",
    "provenance": "provenance.json",
    "checks": "checks.tsv",
    "ranks": "hcc_pair_relative_rank.tsv",
    "composition": "biological_id_compositions.tsv",
    "raw_pairwise": "raw_pairwise_distances.tsv",
}
REGION_REQUIRED = {
    "summary": "summary.json",
    "provenance": "provenance.json",
    "checks": "checks.tsv",
    "annotations": "region_possible_clone_annotations.tsv.gz",
    "pair": "hcc1395_pair_region_possible_clone_comparison.tsv.gz",
    "sample_summary": "sample_possible_clone_summary.tsv",
}
CENSUS_REQUIRED = {
    "summary": "20260712_vaf_final_single_shape_four_class_summary.tsv",
    "census": "20260712_vaf_final_single_shape_four_class_census.json",
    "checks": "20260712_vaf_final_single_shape_checks.tsv",
}

DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

CATEGORY_SPECS = [
    ("single_only", "Single", "single_only"),
    ("sister_only", "Sister", "sister_only"),
    ("direct_only", "Direct", "direct_only"),
    ("sister_and_direct", "Sister＋Direct", "sister_and_direct"),
    ("unresolved", "Unresolved", "unresolved_regions"),
]

SUBCLONE_STATES = {
    "single_modeled_subclonal_like_cluster_represented",
    "multiple_subclonal_cluster_candidate",
    "clonal_plus_single_subclonal_cluster_candidate",
    "clonal_plus_multiple_subclonal_clusters_candidate",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rerun-dir", type=Path, default=RERUN_ROOT)
    parser.add_argument("--region-dir", type=Path, default=REGION_ROOT)
    parser.add_argument("--census-dir", type=Path, default=CENSUS_ROOT)
    parser.add_argument("--bridge-dir", type=Path, default=BRIDGE_ROOT)
    parser.add_argument("--output", type=Path, default=SCRIPT_PATH.parent / "artifact.json")
    parser.add_argument(
        "--representative-max-rows",
        type=int,
        default=120,
        help="Hard ceiling for embedded region-level rows; all informative rows must fit.",
    )
    return parser.parse_args()


def require_files(root: Path, required: Mapping[str, str], label: str) -> dict[str, Path]:
    paths = {key: root / relative for key, relative in required.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"{label} is incomplete: {missing}")
    return paths


def read_json(path: Path) -> Any:
    with path.open() as handle:
        value = json.load(handle)
    if value is None:
        raise ValueError(f"{path} is empty JSON")
    return value


def typed(value: str) -> Any:
    value = value.strip()
    if value == "":
        return None
    if value in {"True", "False"}:
        return value == "True"
    try:
        number = int(value)
        return number
    except ValueError:
        pass
    try:
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"Non-finite value: {value}")
        return number
    except ValueError:
        return value


def read_tsv(path: Path, required_columns: Iterable[str]) -> list[dict[str, Any]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = set(required_columns).difference(columns)
        if missing:
            raise ValueError(f"{path} lacks required columns: {sorted(missing)}")
        rows = [{key: typed(value) for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path} has no data rows")
    return rows


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError as exc:
        raise ValueError(f"Source paths must remain repo-relative: {path}") from exc


def pct(value: float, digits: int = 2) -> str:
    return f"{100 * value:.{digits}f}%"


def jensen_shannon_distance_base2(counts_a: Mapping[str, int], counts_b: Mapping[str, int]) -> float:
    keys = sorted(set(counts_a).union(counts_b))
    total_a = sum(counts_a.values())
    total_b = sum(counts_b.values())
    p = [counts_a.get(key, 0) / total_a for key in keys]
    q = [counts_b.get(key, 0) / total_b for key in keys]
    midpoint = [(a + b) / 2 for a, b in zip(p, q)]
    divergence = 0.0
    for values in (p, q):
        divergence += 0.5 * sum(value * math.log2(value / middle) for value, middle in zip(values, midpoint) if value > 0)
    return math.sqrt(divergence)


def count_checks(path: Path, pass_field: str) -> tuple[int, int]:
    rows = read_tsv(path, {pass_field})
    passed = sum(bool(row[pass_field]) for row in rows)
    return passed, len(rows)


def table_spec(
    table_id: str,
    title: str,
    subtitle: str,
    dataset: str,
    source_id: str,
    columns: Sequence[tuple[str, str, str]],
    *,
    default_sort: tuple[str, str] | None = None,
) -> dict[str, Any]:
    table: dict[str, Any] = {
        "id": table_id,
        "title": title,
        "subtitle": subtitle,
        "dataset": dataset,
        "sourceId": source_id,
        "density": "compact",
        "columns": [
            {"field": field, "label": label, **({"format": fmt} if fmt != "text" else {"type": "text"})}
            for field, label, fmt in columns
        ],
    }
    sort = default_sort or (columns[0][0], "asc")
    table["defaultSort"] = {"field": sort[0], "direction": sort[1]}
    return table


def normalize_inputs(args: argparse.Namespace) -> dict[str, Any]:
    rerun_paths = require_files(args.rerun_dir, RERUN_REQUIRED, "sampling-adjustment snapshot")
    region_paths = require_files(args.region_dir, REGION_REQUIRED, "regional annotation snapshot")
    census_paths = require_files(args.census_dir, CENSUS_REQUIRED, "historical topology census")
    bridge_summary_path = args.bridge_dir / "summary.json"
    for required_path in (bridge_summary_path, TOPOLOGY_PAIR_JSON, TOPOLOGY_PAIR_SCRIPT):
        if not required_path.is_file():
            raise FileNotFoundError(required_path)

    rerun_summary = read_json(rerun_paths["summary"])
    region_summary = read_json(region_paths["summary"])
    census_json = read_json(census_paths["census"])
    bridge_summary = read_json(bridge_summary_path)
    topology_pair = read_json(TOPOLOGY_PAIR_JSON)

    if (
        rerun_summary.get("schema_name") != "intersubmod.crosssample_structure_bulk_sampling_adjustment"
        or rerun_summary.get("schema_version") != "1.0.0"
        or rerun_summary.get("status") != "PASS_WITH_CAVEATS"
        or rerun_summary.get("task_type") != "B_comprehensive_validation"
    ):
        raise RuntimeError("Sampling-adjustment snapshot is not the reviewed comprehensive v1 snapshot")
    if (
        region_summary.get("status") != "PASS_WITH_CLAIM_CEILING"
        or bridge_summary.get("status") != "PASS"
        or census_json.get("status") != "PASS"
    ):
        raise RuntimeError("One or more evidence snapshots are not PASS")
    if topology_pair.get("status") != "PASS" or topology_pair.get("task_type") != "B_comprehensive_validation":
        raise RuntimeError("Matched topology-pair analysis is not the reviewed comprehensive PASS snapshot")

    rerun_pass, rerun_total = count_checks(rerun_paths["checks"], "pass")
    region_pass, region_total = count_checks(region_paths["checks"], "pass")
    if (rerun_pass, rerun_total) != (81, 81):
        raise RuntimeError(f"Sampling QA changed: {rerun_pass}/{rerun_total}")
    if (region_pass, region_total) != (32, 32):
        raise RuntimeError(f"Regional annotation QA changed: {region_pass}/{region_total}")

    for output_key, path_key in (
        ("annotations", "annotations"),
        ("checks", "checks"),
        ("pair_comparison", "pair"),
        ("sample_summary", "sample_summary"),
    ):
        recorded = region_summary["outputs"][output_key]
        path = region_paths[path_key]
        if recorded["sha256"] != sha256_file(path) or int(recorded["size_bytes"]) != path.stat().st_size:
            raise RuntimeError(f"Regional output receipt mismatch: {output_key}")

    census_rows = read_tsv(
        census_paths["summary"],
        {
            "sample", "primary_regions", "complete_regions", "incomplete_regions",
            "final_single_shape_regions", "unresolved_regions", "single_only", "sister_only",
            "direct_only", "sister_and_direct",
        },
    )
    dataset_rows = [row for row in census_rows if row["sample"] != "aggregate"]
    if [row["sample"] for row in dataset_rows] != DATASETS:
        raise RuntimeError("Seven-dataset census order changed")
    for row in dataset_rows:
        resolved_sum = sum(int(row[field]) for field in ("single_only", "sister_only", "direct_only", "sister_and_direct"))
        if resolved_sum != int(row["final_single_shape_regions"]):
            raise RuntimeError(f"Resolved category conservation failed: {row['sample']}")
        if resolved_sum + int(row["unresolved_regions"]) != int(row["complete_regions"]):
            raise RuntimeError(f"Complete category conservation failed: {row['sample']}")

    annotation_rows = read_tsv(
        region_paths["annotations"],
        {
            "sample", "region", "complete", "C_read_groups_total", "T_exact_candidate_forest_count",
            "Topo_shape_count", "tree_selection_source", "vaf_most_likely_topology_status",
            "vaf_inferred_final_shape", "claim_ceiling",
        },
    )
    if len(annotation_rows) != 47_377:
        raise RuntimeError(f"Annotation population changed: {len(annotation_rows)}")
    if len({(row["sample"], row["region"]) for row in annotation_rows}) != len(annotation_rows):
        raise RuntimeError("Annotation sample+region keys are not unique")

    complete_rows = [row for row in annotation_rows if row["complete"] is True]
    if len(complete_rows) != 39_885:
        raise RuntimeError(f"Complete-region population changed: {len(complete_rows)}")
    c0 = sum(int(row["C_read_groups_total"]) == 0 for row in annotation_rows)
    c1 = sum(int(row["C_read_groups_total"]) == 1 for row in annotation_rows)
    cgt1 = len(annotation_rows) - c0 - c1
    t1 = sum(int(row["T_exact_candidate_forest_count"]) == 1 for row in complete_rows)
    topo1 = sum(int(row["Topo_shape_count"]) == 1 for row in complete_rows)
    tgt1_topo1 = sum(
        int(row["T_exact_candidate_forest_count"]) > 1 and int(row["Topo_shape_count"]) == 1
        for row in complete_rows
    )
    vaf_resolved = sum(row["tree_selection_source"] == "vaf_resolved_topogt1" for row in annotation_rows)
    vaf_unresolved = sum(row["tree_selection_source"] == "vaf_unresolved_topogt1" for row in annotation_rows)
    observed_ladder = (c0, c1, cgt1, t1, len(complete_rows) - t1, topo1, len(complete_rows) - topo1, vaf_resolved, vaf_unresolved)
    expected_ladder = (26_593, 10_322, 10_462, 10_832, 29_053, 21_976, 17_909, 15_063, 2_846)
    if observed_ladder != expected_ladder:
        raise RuntimeError(f"C/T/Topo/VAF ladder changed: {observed_ladder}")
    if tgt1_topo1 != 11_144:
        raise RuntimeError(f"T>1 but Topo=1 count changed: {tgt1_topo1}")

    sample_summary_rows = read_tsv(
        region_paths["sample_summary"], {"sample", "endpoint", "status", "regions", "sample_regions", "share"}
    )
    if len(sample_summary_rows) != 83:
        raise RuntimeError(f"Sample state summary row count changed: {len(sample_summary_rows)}")
    for sample in DATASETS:
        for endpoint in ("all_joined_diagnostic", "full_exact_join", "full_exact_assignment_ge_0.8"):
            rows = [row for row in sample_summary_rows if row["sample"] == sample and row["endpoint"] == endpoint]
            if not rows or sum(int(row["regions"]) for row in rows) != int(rows[0]["sample_regions"]):
                raise RuntimeError(f"Sample state conservation failed: {sample}/{endpoint}")

    pair_rows = read_tsv(
        region_paths["pair"],
        {
            "match_id", "chrom", "region_a", "region_b", "final_shape_category_a",
            "final_shape_category_b", "C_read_groups_total_a", "C_read_groups_total_b",
            "regional_possible_state_a", "regional_possible_state_b", "possible_subclone_signal_a",
            "possible_subclone_signal_b", "cross_source_possible_state", "external_partition_exact",
            "external_partition_informative", "external_partition_pattern",
            "external_partition_vacuity_reason", "external_common_joined_mutations",
            "external_subclonal_intersection_n", "external_subclonal_union_n",
            "external_subclonal_jaccard_defined", "external_subclonal_jaccard",
            "external_subclonal_jaccard_vacuity_reason",
            "body_gene_symbols", "promoter_gene_symbols", "cgc_body_symbols", "cgc_promoter_symbols",
            "dgidb_approved_antineoplastic_body_symbols",
            "dgidb_approved_antineoplastic_promoter_symbols",
            "approved_antineoplastic_body_drug_claim_names", "clp_all_variant_count",
            "clp_confirmed_somatic_variant_count",
        },
    )
    if len(pair_rows) != 5_720 or len({row["match_id"] for row in pair_rows}) != 5_720:
        raise RuntimeError("HCC fixed pair population or key uniqueness changed")
    evaluable = [row for row in pair_rows if row["cross_source_possible_state"] != "not_evaluable_missing_external_cluster"]
    same = [row for row in evaluable if row["cross_source_possible_state"] == "same_possible_state"]
    both_single_clonal = [
        row for row in same
        if row["regional_possible_state_a"] == "single_modeled_clonal_like_cluster_represented"
        and row["regional_possible_state_b"] == "single_modeled_clonal_like_cluster_represented"
    ]
    either_subclone = [row for row in pair_rows if row["possible_subclone_signal_a"] or row["possible_subclone_signal_b"]]
    both_subclone = [row for row in pair_rows if row["possible_subclone_signal_a"] and row["possible_subclone_signal_b"]]
    informative = [row for row in pair_rows if row["external_partition_informative"] is True]
    informative_exact = [row for row in informative if row["external_partition_exact"] is True]
    jaccard_undefined = [row for row in pair_rows if row["external_subclonal_jaccard_defined"] is False]
    undefined_union_zero = [row for row in jaccard_undefined if int(row["external_subclonal_union_n"]) == 0]
    gene_drug_informative = [
        row for row in pair_rows
        if (row["possible_subclone_signal_a"] or row["possible_subclone_signal_b"])
        and any(
            row[field]
            for field in (
                "cgc_body_symbols", "cgc_promoter_symbols",
                "dgidb_approved_antineoplastic_body_symbols",
                "dgidb_approved_antineoplastic_promoter_symbols",
            )
        )
    ]
    observed_pair = (len(evaluable), len(same), len(both_single_clonal), len(either_subclone), len(both_subclone), len(informative), len(informative_exact))
    expected_pair = (4_438, 4_296, 4_267, 172, 40, 34, 21)
    if observed_pair != expected_pair:
        raise RuntimeError(f"HCC pair decomposition changed: {observed_pair}")
    if len(jaccard_undefined) != 5_520 or len(undefined_union_zero) != 5_520:
        raise RuntimeError(f"Jaccard vacuity semantics changed: {len(jaccard_undefined)}/{len(undefined_union_zero)}")
    if len(informative) > args.representative_max_rows:
        raise RuntimeError(f"Informative region table exceeds bound: {len(informative)}")
    if len(gene_drug_informative) != 15 or len(informative) + len(gene_drug_informative) > args.representative_max_rows:
        raise RuntimeError(f"Gene/drug shortlist or combined row bound changed: {len(gene_drug_informative)}")

    exact_pair_metrics = next(
        row for row in topology_pair["hcc1395_pair_metrics"] if row["scenario"] == "exact_coordinate"
    )
    if (
        int(exact_pair_metrics["complete_both"]) != 5_720
        or round(float(exact_pair_metrics["raw_agreement"]) * 5_720) != 3_969
        or abs(float(exact_pair_metrics["permutation_null"]["agreement_mean"]) - 0.3951152797202797) > 1e-12
    ):
        raise RuntimeError("Matched pre-VAF topology endpoint changed")

    matched_final_agreement_n = sum(
        row["final_shape_category_a"] == row["final_shape_category_b"] for row in pair_rows
    )
    final_counts_a: dict[str, int] = {}
    final_counts_b: dict[str, int] = {}
    for row in pair_rows:
        final_counts_a[str(row["final_shape_category_a"])] = final_counts_a.get(str(row["final_shape_category_a"]), 0) + 1
        final_counts_b[str(row["final_shape_category_b"])] = final_counts_b.get(str(row["final_shape_category_b"]), 0) + 1
    matched_final_jsd = jensen_shannon_distance_base2(final_counts_a, final_counts_b)
    if matched_final_agreement_n != 4_243 or abs(matched_final_jsd - 0.19666903152894216) > 1e-12:
        raise RuntimeError(f"Matched final-shape descriptive endpoint changed: {matched_final_agreement_n}/{matched_final_jsd}")

    raw_pairwise = read_tsv(
        rerun_paths["raw_pairwise"],
        {"comparison_level", "estimand", "entity_a", "entity_b", "biological_id_a", "biological_id_b", "pair_type", "jensen_shannon_distance_base2"},
    )
    hcc_jsd = float(rerun_summary["hcc_raw"]["complete_five_class_primary"]["jensen_shannon_distance_base2"])
    cross_pairs = [
        row for row in raw_pairwise
        if row["comparison_level"] == "dataset_row"
        and row["estimand"] == "complete_five_class_primary"
        and row["pair_type"] == "cross_biological_id"
    ]
    closer_pairs = [row for row in cross_pairs if float(row["jensen_shannon_distance_base2"]) < hcc_jsd]
    closer_biological_pairs = {
        tuple(sorted((str(row["biological_id_a"]), str(row["biological_id_b"])))) for row in closer_pairs
    }
    if (len(cross_pairs), len(closer_pairs), len(closer_biological_pairs)) != (20, 8, 7):
        raise RuntimeError("HCC rank decomposition changed")

    return {
        "paths": {
            "rerun": rerun_paths,
            "region": region_paths,
            "census": census_paths,
            "bridge_summary": bridge_summary_path,
            "topology_pair_json": TOPOLOGY_PAIR_JSON,
            "topology_pair_script": TOPOLOGY_PAIR_SCRIPT,
        },
        "rerun_summary": rerun_summary,
        "region_summary": region_summary,
        "bridge_summary": bridge_summary,
        "topology_pair": topology_pair,
        "matched": {
            "pre_vaf": exact_pair_metrics,
            "final_agreement_n": matched_final_agreement_n,
            "final_jsd": matched_final_jsd,
            "cross_pairs": len(cross_pairs),
            "closer_pairs": len(closer_pairs),
            "closer_biological_pairs": len(closer_biological_pairs),
        },
        "census_rows": dataset_rows,
        "annotation_rows": annotation_rows,
        "complete_rows": complete_rows,
        "sample_summary_rows": sample_summary_rows,
        "pair_rows": pair_rows,
        "informative_rows": informative,
        "gene_drug_informative_rows": gene_drug_informative,
        "qa": {"sampling": (rerun_pass, rerun_total), "regional": (region_pass, region_total)},
        "ladder": {
            "c0": c0, "c1": c1, "cgt1": cgt1, "t1": t1, "tgt1": len(complete_rows) - t1,
            "topo1": topo1, "topogt1": len(complete_rows) - topo1,
            "tgt1_topo1": tgt1_topo1,
            "vaf_resolved": vaf_resolved, "vaf_unresolved": vaf_unresolved,
        },
        "pair": {
            "evaluable": len(evaluable), "same": len(same), "both_single_clonal": len(both_single_clonal),
            "either_subclone": len(either_subclone), "both_subclone": len(both_subclone),
            "informative": len(informative), "informative_exact": len(informative_exact),
            "jaccard_undefined_union_zero": len(undefined_union_zero),
        },
    }


def endpoint_counts(rows: Sequence[Mapping[str, Any]], sample: str, endpoint: str) -> dict[str, int]:
    selected = [row for row in rows if row["sample"] == sample and row["endpoint"] == endpoint]
    counts = {str(row["status"]): int(row["regions"]) for row in selected}
    unavailable = sum(value for key, value in counts.items() if key.startswith("unavailable"))
    clonal_like = counts.get("single_modeled_clonal_like_cluster_represented", 0)
    subclone = sum(counts.get(key, 0) for key in SUBCLONE_STATES)
    total = int(selected[0]["sample_regions"])
    if unavailable + clonal_like + subclone != total:
        raise RuntimeError(f"Unexpected state classification: {sample}/{endpoint}/{counts}")
    return {"total": total, "unavailable": unavailable, "clonal_like": clonal_like, "subclone": subclone}


def build_artifact(data: Mapping[str, Any]) -> dict[str, Any]:
    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    paths = data["paths"]
    rerun = data["rerun_summary"]
    region = data["region_summary"]
    bridge = data["bridge_summary"]
    matched = data["matched"]
    pre_vaf = matched["pre_vaf"]
    ladder = data["ladder"]
    pair = data["pair"]

    composition_rows: list[dict[str, Any]] = []
    topology_table: list[dict[str, Any]] = []
    direct_complete_shares: list[float] = []
    for row in data["census_rows"]:
        complete = int(row["complete_regions"])
        sample = str(row["sample"])
        for category, label, count_field in CATEGORY_SPECS:
            count = int(row[count_field])
            composition_rows.append({"sample": sample, "category": label, "category_id": category, "count": count})
        direct_share = int(row["direct_only"]) / complete
        direct_complete_shares.append(direct_share)
        topology_table.append(
            {
                "sample": sample,
                "primary_regions": int(row["primary_regions"]),
                "complete_regions": complete,
                "incomplete_regions": int(row["incomplete_regions"]),
                "single": f"{int(row['single_only']):,} ({pct(int(row['single_only']) / complete)})",
                "sister": f"{int(row['sister_only']):,} ({pct(int(row['sister_only']) / complete)})",
                "direct": f"{int(row['direct_only']):,} ({pct(direct_share)})",
                "mixed": f"{int(row['sister_and_direct']):,} ({pct(int(row['sister_and_direct']) / complete)})",
                "unresolved": f"{int(row['unresolved_regions']):,} ({pct(int(row['unresolved_regions']) / complete)})",
            }
        )
    if not all(max(row["count"] for row in composition_rows if row["sample"] == sample) == next(
        row["count"] for row in composition_rows if row["sample"] == sample and row["category_id"] == "direct_only"
    ) for sample in DATASETS):
        raise RuntimeError("Direct is no longer the dominant complete five-class category in all seven datasets")

    rank_labels = {
        "raw_unadjusted": "Raw",
        "dirichlet_posterior": "Dirichlet",
        "rarefaction_without_replacement_equal_n_2557": "Rarefaction",
        "paired_common_genomic_block_bootstrap_5mb": "5 Mb block",
        "paired_common_genomic_block_bootstrap_10mb": "10 Mb block",
        "empirical_bayes:leave_hcc1395_out_5_biological_ids:x2": "EB leave-HCC-out ×2",
        "source_standardized:global_complete_three_source_mix_sensitivity": "Source-standardized",
    }
    rank_rows = []
    methods = rerun["primary_jsd_adjustment_stability"]["methods"]
    for method_id, label in rank_labels.items():
        method = methods[method_id]
        rank_rows.append(
            {
                "method": label,
                "method_id": method_id,
                "jsd_median": method["hcc_distance_median"],
                "rank": method["hcc_rank_smaller_is_more_similar_median"],
                "rank_q025": method["hcc_rank_q025"],
                "rank_q975": method["hcc_rank_q975"],
                "rank_1_probability": method["probability_hcc_rank_1"],
                "replicates": method["replicates"],
            }
        )

    coverage_chart = []
    coverage_table = []
    assignment_table = []
    for sample in DATASETS:
        primary = endpoint_counts(data["sample_summary_rows"], sample, "full_exact_join")
        high = endpoint_counts(data["sample_summary_rows"], sample, "full_exact_assignment_ge_0.8")
        sample_meta = region["samples"][sample]
        sample_annotations = [row for row in data["annotation_rows"] if row["sample"] == sample]
        cn_available = bool(sample_annotations[0].get("external_cn_fit_available"))
        cn_independent = sample_annotations[0].get("external_cn_source_independent")
        cluster_count = sample_annotations[0].get("external_global_fit_cluster_count")
        for label, count in (
            ("Clonal-like", primary["clonal_like"]),
            ("Subclone candidate", primary["subclone"]),
            ("Unavailable", primary["unavailable"]),
        ):
            coverage_chart.append({"sample": sample, "status": label, "regions": count})
        coverage_table.append(
            {
                "sample": sample,
                "primary_regions": primary["total"],
                "region_linked_sites": int(sample_meta["sites"]),
                "mutation_join_fraction": float(sample_meta["join_fraction"]),
                "full_exact_regions": int(sample_meta["regions_full_join"]),
                "full_exact_share": int(sample_meta["regions_full_join"]) / primary["total"],
                "clonal_like_regions": primary["clonal_like"],
                "possible_subclone_regions": primary["subclone"],
                "unavailable_regions": primary["unavailable"],
                "modeled_clusters": cluster_count,
                "cn_gate": (
                    "PASS, source-specific CN" if cn_available and cn_independent is True
                    else "PASS, shared HCC CN sensitivity" if cn_available
                    else "BLOCKED: no reviewed allele-specific CN"
                ),
            }
        )
        assignment_table.append(
            {
                "sample": sample,
                "primary_evaluable": primary["total"] - primary["unavailable"],
                "primary_evaluable_share": (primary["total"] - primary["unavailable"]) / primary["total"],
                "primary_subclone": primary["subclone"],
                "assignment_ge_0_8_evaluable": high["total"] - high["unavailable"],
                "assignment_ge_0_8_evaluable_share": (high["total"] - high["unavailable"]) / high["total"],
                "assignment_ge_0_8_subclone": high["subclone"],
                "fail_closed": high["unavailable"],
            }
        )

    aggregate_ladder = [
        {
            "layer": "Region population",
            "state": "primary / complete / incomplete",
            "counts": "47,377 / 39,885 / 7,492",
            "meaning": "七 dataset rows 的 historical primary regions；complete 才進入 T/Topo 判讀。",
        },
        {
            "layer": "C_read_groups",
            "state": "C=0 / C=1 / C>1",
            "counts": f"{ladder['c0']:,} / {ladder['c1']:,} / {ladder['cgt1']:,}",
            "meaning": "read-supported R/A mutation-state groups；不是 clone 數。",
        },
        {
            "layer": "T candidate structure",
            "state": "T=1 / T>1 among complete",
            "counts": f"{ladder['t1']:,} / {ladder['tgt1']:,}",
            "meaning": "mutation-labeled exact candidate forest count；多個 T 可共享一個無標籤 shape。",
        },
        {
            "layer": "Topo shape",
            "state": "Topo=1 / Topo>1 among complete",
            "counts": f"{ladder['topo1']:,} / {ladder['topogt1']:,}",
            "meaning": f"rooted-unlabeled shape multiplicity；Topo=1 中只有 {ladder['t1']:,} 個 T=1，另有 {ladder['tgt1_topo1']:,} 個 T>1。",
        },
        {
            "layer": "VAF selection within Topo>1",
            "state": "resolved most-likely / unresolved",
            "counts": f"{ladder['vaf_resolved']:,} / {ladder['vaf_unresolved']:,}",
            "meaning": "用各位點 VAF 排序候選 T/shape 的推測；不是獨立 truth validation。",
        },
    ]

    definitions = [
        {"term": "S / sSNV entries", "definition": "每個 dataset 的 region-linked sSNV entries；跨 dataset 不合併成唯一生物位點分母。", "not_equal": "clone count"},
        {"term": "C_read_groups", "definition": region["definitions"]["C_read_groups"], "not_equal": "external cluster / clone count"},
        {"term": "T", "definition": region["definitions"]["T_exact_candidate_forest_count"], "not_equal": "Topo shape count"},
        {"term": "Topo", "definition": region["definitions"]["Topo_shape_count"], "not_equal": "mutation-labeled unique tree truth"},
        {"term": "VAF-selected result", "definition": region["definitions"]["tree_selection_source"], "not_equal": "experimentally confirmed tree"},
        {"term": "External PyClone cluster", "definition": region["definitions"]["external_cluster_count"], "not_equal": "cross-fit stable clone identity or ancestry"},
        {"term": "Clonal / subclonal", "definition": "PyClone conditional CP≥0.90 / CP<0.90 in this analysis.", "not_equal": "single-cell observed state"},
        {"term": "Subclonal Jaccard", "definition": "|A∩B|/|A∪B| only when the subclonal union is >0; union=0 is NA.", "not_equal": "perfect agreement when neither side has a subclonal member"},
    ]

    hcc_decomposition = [
        {"endpoint": "Full exact external state evaluable", "basis": "full-exact region state", "numerator": pair["evaluable"], "denominator": 5_720, "rate": pair["evaluable"] / 5_720, "interpretation": "1,282 regions fail closed before cross-source state comparison."},
        {"endpoint": "Same possible state", "basis": "full-exact region state", "numerator": pair["same"], "denominator": pair["evaluable"], "rate": pair["same"] / pair["evaluable"], "interpretation": "表面很高，但 4,267 rows 是雙側皆 single modeled clonal-like。"},
        {"endpoint": "Both-single-clonal share of same-state", "basis": "full-exact region state", "numerator": pair["both_single_clonal"], "denominator": pair["same"], "rate": pair["both_single_clonal"] / pair["same"], "interpretation": "主群 dominance 解釋 99.33% 的 same-state rows。"},
        {"endpoint": "Possible subclone regional Jaccard", "basis": "full-exact region signal", "numerator": pair["both_subclone"], "denominator": pair["either_subclone"], "rate": pair["both_subclone"] / pair["either_subclone"], "interpretation": "40/172；minor-signal-focused cross-source overlap。"},
        {"endpoint": "Informative both-multicluster partition exact", "basis": "all-joined common-mutation partition", "numerator": pair["informative_exact"], "denominator": pair["informative"], "rate": pair["informative_exact"] / pair["informative"], "interpretation": "21/34；可含 incomplete full-exact join，與 4,438 endpoint 不是同一分母。"},
        {"endpoint": "Assignment≥0.8 sensitivity evaluable", "basis": "full exact + assignment gate", "numerator": 2_370, "denominator": 5_720, "rate": 2_370 / 5_720, "interpretation": "全部 same-state 但 minor groups 幾乎被移除；selection-induced degeneracy。"},
    ]

    informative_rows = []
    for row in sorted(data["informative_rows"], key=lambda item: (str(item["chrom"]), str(item["region_a"]))):
        informative_rows.append(
            {
                "match_id": row["match_id"],
                "chrom": row["chrom"],
                "region": row["region_a"],
                "shape_hcc": row["final_shape_category_a"],
                "shape_dorado": row["final_shape_category_b"],
                "C_hcc": row["C_read_groups_total_a"],
                "C_dorado": row["C_read_groups_total_b"],
                "state_hcc": row["regional_possible_state_a"],
                "state_dorado": row["regional_possible_state_b"],
                "cross_source_state": row["cross_source_possible_state"],
                "partition": "exact" if row["external_partition_exact"] else "different",
                "common_joined_mutations": row["external_common_joined_mutations"],
                "subclonal_intersection": row["external_subclonal_intersection_n"],
                "subclonal_union": row["external_subclonal_union_n"],
                "jaccard_defined": row["external_subclonal_jaccard_defined"],
                "subclone_jaccard": row["external_subclonal_jaccard"],
                "vacuity_reason": row["external_subclonal_jaccard_vacuity_reason"] or "—",
            }
        )

    annotation_index = {
        (str(row["sample"]), str(row["region"])): row for row in data["annotation_rows"]
    }
    gene_drug_examples = []
    for row in sorted(data["gene_drug_informative_rows"], key=lambda item: (str(item["chrom"]), str(item["region_a"]))):
        annotation_a = annotation_index[("HCC1395", str(row["region_a"]))]
        annotation_b = annotation_index[("HCC1395_DORADO", str(row["region_b"]))]
        cgc = "|".join(
            value for value in (row["cgc_body_symbols"], row["cgc_promoter_symbols"]) if value
        ) or "—"
        approved_targets = "|".join(
            value
            for value in (
                row["dgidb_approved_antineoplastic_body_symbols"],
                row["dgidb_approved_antineoplastic_promoter_symbols"],
            )
            if value
        ) or "—"
        gene_drug_examples.append(
            {
                "chrom": row["chrom"],
                "region": row["region_a"],
                "body_genes": row["body_gene_symbols"] or "—",
                "cgc_symbols": cgc,
                "dgidb_targets": approved_targets,
                "drug_claim_names": row["approved_antineoplastic_body_drug_claim_names"] or "—",
                "clp_all": int(row["clp_all_variant_count"]),
                "clp_confirmed": int(row["clp_confirmed_somatic_variant_count"]),
                "C_hcc": int(annotation_a["C_read_groups_total"]),
                "C_dorado": int(annotation_b["C_read_groups_total"]),
                "assignment_min_hcc": annotation_a.get("external_assignment_min"),
                "assignment_min_dorado": annotation_b.get("external_assignment_min"),
                "state_hcc": row["regional_possible_state_a"],
                "state_dorado": row["regional_possible_state_b"],
                "partition": row["external_partition_pattern"],
            }
        )
    highlighted_regions = {
        "chr11:128802981-128853809", "chr2:141063678-141096354", "chr19:44784709-44850535"
    }
    highlighted = [row for row in gene_drug_examples if row["region"] in highlighted_regions]
    if (
        len(highlighted) != 3
        or any(row["C_hcc"] != 0 or row["C_dorado"] != 0 for row in highlighted)
        or any(min(float(row["assignment_min_hcc"]), float(row["assignment_min_dorado"])) >= 0.8 for row in highlighted)
    ):
        raise RuntimeError("Highlighted gene examples no longer satisfy C=0 and low-assignment caveat")

    gene_drug_summary = [
        {"annotation": "Body gene", "flagged_regions": 4_459, "coarse_agree": 3_104, "either_subclone": 128, "both_subclone": 29},
        {"annotation": "CGC body/promoter (any)", "flagged_regions": 275, "coarse_agree": 199, "either_subclone": 6, "both_subclone": 1},
        {"annotation": "DGIdb approved-antineoplastic (any)", "flagged_regions": 464, "coarse_agree": 334, "either_subclone": 10, "both_subclone": 0},
        {"annotation": "Cancer Licence Project (any)", "flagged_regions": 333, "coarse_agree": 233, "either_subclone": 12, "both_subclone": 3},
    ]

    global_concordance = bridge["global_concordance"]
    regional_concordance = bridge["regional_concordance"]
    if (
        (global_concordance["clusters_a"], global_concordance["clusters_b"],
         global_concordance["subclonal_intersection"], global_concordance["subclonal_union"])
        != (4, 3, 277, 727)
        or (regional_concordance["subclonal_intersection"], regional_concordance["subclonal_union"])
        != (62, 259)
    ):
        raise RuntimeError("Independent PyClone global/regional context changed")
    pyclone_context = [
        {
            "population": "Global separate-fit shared universe",
            "mutations": int(global_concordance["n"]),
            "fit_local_clusters": f"{global_concordance['clusters_a']} vs {global_concordance['clusters_b']}",
            "clonal_state_agreement": float(global_concordance["clonal_state_agreement"]),
            "kappa": float(global_concordance["clonal_state_kappa"]),
            "subclonal_intersection": int(global_concordance["subclonal_intersection"]),
            "subclonal_union": int(global_concordance["subclonal_union"]),
            "subclonal_jaccard": float(global_concordance["subclonal_jaccard"]),
            "interpretation": "主群 agreement 很高；κ 與 minor-set Jaccard 只屬中度再現。",
        },
        {
            "population": "Fixed-region all-joined mutation subset",
            "mutations": int(regional_concordance["n"]),
            "fit_local_clusters": f"{regional_concordance['clusters_a']} vs {regional_concordance['clusters_b']}",
            "clonal_state_agreement": float(regional_concordance["clonal_state_agreement"]),
            "kappa": float(regional_concordance["clonal_state_kappa"]),
            "subclonal_intersection": int(regional_concordance["subclonal_intersection"]),
            "subclonal_union": int(regional_concordance["subclonal_union"]),
            "subclonal_jaccard": float(regional_concordance["subclonal_jaccard"]),
            "interpretation": "區域 mutation subset 的 minor overlap 再下降；不是 full-exact region-state endpoint。",
        },
    ]

    claims = [
        {"claim": "七樣本皆 Direct-dominant", "verdict": "描述性共同輸出模式", "evidence": f"7/7；complete-region share {pct(min(direct_complete_shares))}–{pct(max(direct_complete_shares))}。", "ceiling": "尚無 solver/class-space/eligibility null；不可單獨當方法有效性證據。"},
        {"claim": "HCC same-locus pre-VAF coarse structure 跨來源再現", "verdict": "部分支持", "evidence": "3,969/5,720=69.39%（95% CI 68.18%–70.57%），高於 chromosome-preserving null 39.51%，p=1/5001。", "ceiling": "same-cell-line cross-source reproducibility；不是 final-shape 或 exact-tree accuracy。"},
        {"claim": "HCC pair 完整 composition 高度相似", "verdict": "不支持", "evidence": "JSD 在 raw/Dirichlet/rarefaction/5/10 Mb/EB rank=9/21；source-standardized rank=10/21；P(rank 1)=0。", "ceiling": "sampling/source 校正不會讓它成為最接近 pair。"},
        {"claim": "HCC pair external regional state 高一致", "verdict": "表面支持、需拆解", "evidence": "4,296/4,438=96.80%，其中 4,267 為雙 single-clonal-like。", "ceiling": "主要驗證 clonal-majority，不是 minor subclone identity。"},
        {"claim": "HCC minor/subclone 結構再現", "verdict": "部分支持", "evidence": "full-exact possible-subclone Jaccard=40/172=23.26%；另一路徑 all-joined common-mutation both-multicluster partition exact=21/34=61.76%。", "ceiling": "兩個不同 endpoint basis 的區域 candidate signal；不是 exact ancestry。"},
        {"claim": "每區域唯一準確真實樹", "verdict": "不支持", "evidence": "complete regions 中 T>1=29,053、Topo>1=17,909；15,063 是 VAF-ranked 推測而非 truth。", "ceiling": "需要 single-cell / synthetic truth / orthogonal lineage evidence。"},
    ]

    estimand_comparison = [
        {
            "estimand": "Aggregate complete five-class final-shape composition",
            "population": "HCC1395 n=6,940 vs DORADO n=6,750; unmatched aggregate",
            "result": f"JSD={rerun['hcc_raw']['complete_five_class_primary']['jensen_shannon_distance_base2']:.4f}; rank 9/21",
            "null_status": "No solver/class-space null",
            "interpretation": f"{matched['closer_pairs']}/20 dataset-row comparisons are closer, spanning {matched['closer_biological_pairs']} unique biological-ID pairs.",
        },
        {
            "estimand": "Matched pre-VAF coarse topology",
            "population": "5,720 exact-coordinate complete-both regions",
            "result": "3,969/5,720=69.39% agreement",
            "null_status": "Chromosome-preserving null mean=39.51%; p=1/5001",
            "interpretation": "Strongest positive reproducibility endpoint; coarse category, not exact tree.",
        },
        {
            "estimand": "Matched VAF-selected final shape",
            "population": "same 5,720 matched regions",
            "result": f"{matched['final_agreement_n']:,}/5,720={pct(matched['final_agreement_n']/5720)}; composition JSD={matched['final_jsd']:.4f}",
            "null_status": "No matched final-shape permutation null",
            "interpretation": "Descriptive same-population sensitivity; cannot inherit the pre-VAF null.",
        },
        {
            "estimand": "External PyClone regional state / partition",
            "population": "full-exact state 4,438; all-joined informative partition 34",
            "result": "same-state 4,296/4,438; informative exact 21/34",
            "null_status": "No biological clone-tree truth",
            "interpretation": "Different model and eligibility population; clonal-majority dominated.",
        },
    ]

    qa_rows = [
        {"artifact": "bulk-sampling adjustment", "status": rerun["status"], "checks": "81/81 PASS", "population": "7 dataset rows; 39,885 complete regions", "hash_anchor": sha256_file(paths["rerun"]["summary"])},
        {"artifact": "regional possible-clone annotations", "status": region["status"], "checks": "32/32 PASS", "population": "47,377 annotations; 5,720 HCC pair matches", "hash_anchor": sha256_file(paths["region"]["summary"])},
        {"artifact": "matched pre-VAF coarse topology", "status": "PASS", "checks": "5,000 chromosome-preserving permutations", "population": "5,720 exact-coordinate complete-both regions", "hash_anchor": sha256_file(paths["topology_pair_json"])},
        {"artifact": "historical final-shape census", "status": "PASS", "checks": "source census PASS", "population": "47,377 primary; 39,885 complete", "hash_anchor": sha256_file(paths["census"]["census"])},
    ]

    sources = [
        {
            "id": "historical_topology_census",
            "label": "Seven-dataset historical layered-v2 final-shape census",
            "path": repo_relative(paths["census"]["census"]),
            "query": {
                "engine": "Python validated census snapshot",
                "sql": "SELECT * FROM snapshot.historical_topology_census",
                "language": "python",
                "description": "Seven dataset rows, complete/unresolved denominators, VAF-selected final single-shape categories.",
                "tables_used": [repo_relative(path) for path in paths["census"].values()],
                "filters": ["GRCh38 chr1-22", "historical layered-v2", "primary regions"],
                "metric_definitions": [
                    "complete five-class denominator includes unresolved",
                    "Direct/Sister/Single/Mixed are rooted-unlabeled graph-pattern categories, not clone fractions",
                ],
            },
        },
        {
            "id": "bulk_sampling_adjustment",
            "label": "Cross-sample bulk-sampling and source adjustment rerun A",
            "path": repo_relative(paths["rerun"]["summary"]),
            "query": {
                "engine": "Python Dirichlet + rarefaction + genomic block bootstrap + empirical Bayes",
                "sql": "SELECT * FROM snapshot.bulk_sampling_adjustment",
                "language": "python",
                "description": "Primary complete-five-class composition distances and HCC pair rank against 20 cross-biological pairs.",
                "tables_used": [repo_relative(path) for path in paths["rerun"].values()],
                "filters": ["7 dataset rows", "6 biological IDs", "5/10 Mb common genomic blocks", "81/81 checks PASS"],
                "metric_definitions": [
                    "Jensen-Shannon distance uses base-2 logs; lower is more similar",
                    "rank inserts the HCC same-cell-line cross-source pair beside 20 cross-biological pairs; rank 1 is closest",
                    "source-standardization is sensitivity-only, not an observationally identified correction",
                ],
            },
        },
        {
            "id": "regional_clone_annotations",
            "label": "Region-level C/T/Topo/VAF-selection and PyClone-VI possible-state annotations",
            "path": repo_relative(paths["region"]["summary"]),
            "query": {
                "engine": "Python exact-position bridge to PyClone-VI 0.2.0",
                "sql": "SELECT * FROM snapshot.region_possible_clone_annotations",
                "language": "python",
                "description": "Full 47,377-region annotation, seven-sample coverage, fail-closed CN gates, and fixed 5,720 HCC comparison.",
                "tables_used": [repo_relative(path) for path in paths["region"].values()],
                "filters": ["GRCh38 chr1-22", "CP clonal threshold 0.90", "assignment>=0.80 sensitivity", "32/32 checks PASS"],
                "metric_definitions": [
                    "possible subclone signal is conditional on selected PyClone fit, CN/purity, joined mutations, and region definition",
                    "subclonal Jaccard is defined only when the subclonal union is positive; union=0 is NA rather than 1",
                    "COLO829/HCC1937 fail closed because reviewed allele-specific CN is unavailable",
                    "DORADO reuses HCC1395 CN as a sensitivity assumption and is not source-independent",
                ],
            },
        },
        {
            "id": "matched_topology_pair_analysis",
            "label": "HCC1395 exact-coordinate matched pre-VAF coarse topology analysis",
            "path": repo_relative(paths["topology_pair_json"]),
            "query": {
                "engine": "Python exact-coordinate matching + chromosome-preserving permutation",
                "sql": "SELECT * FROM snapshot.topology_pair_analysis WHERE scenario='exact_coordinate'",
                "language": "python",
                "description": "Fixed 5,720 complete-both regions; five-class pre-VAF coarse agreement, CI, kappa, exact-tree digest, and 5,000-permutation chromosome-preserving null.",
                "tables_used": [repo_relative(paths["topology_pair_json"]), repo_relative(paths["topology_pair_script"])],
                "filters": ["HCC1395 same-cell-line cross-source pair", "GRCh38 chr1-22", "exact coordinate", "both complete"],
                "metric_definitions": [
                    "pre-VAF coarse agreement = 3,969/5,720",
                    "chromosome-preserving null mean = 39.51%; add-one p = 1/5001",
                    "matched final-shape agreement is a separate descriptive endpoint and has no corresponding permutation null",
                ],
            },
        },
        {
            "id": "clone_bridge_redteam",
            "label": "Independent PyClone mutation-cluster to region bridge red-team summary",
            "path": repo_relative(paths["bridge_summary"]),
            "query": {
                "engine": "Python label-invariant partition audit",
                "sql": "SELECT * FROM snapshot.clone_region_bridge WHERE endpoint='both_multicluster'",
                "language": "python",
                "description": "Uses informative both-multicluster partitions; excludes duplicate-key global-vs-regional cluster concordance table.",
                "filters": ["fixed 5,720 HCC regions", "independent separate fits", "label-invariant partitions"],
                "metric_definitions": ["all-joined common-mutation both-multicluster exact partition = 21/34; this is distinct from the full-exact region-state denominator"],
            },
        },
        {
            "id": "pyclone_vi_paper",
            "label": "PyClone-VI: scalable inference of clonal population structures using whole genome data (Gillis & Roth)",
            "href": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7730797/",
        },
        {
            "id": "pyclone_original_paper",
            "label": "PyClone: statistical inference of clonal population structure in cancer (Roth et al.)",
            "href": "https://www.nature.com/articles/nmeth.2883",
        },
    ]

    headline_metrics = [{
        "direct_dominant": 7,
        "hcc_jsd_rank": 9,
        "matched_pre_vaf_coarse": 3_969 / 5_720,
        "same_state": pair["same"] / pair["evaluable"],
        "subclone_jaccard": pair["both_subclone"] / pair["either_subclone"],
        "informative_partition": pair["informative_exact"] / pair["informative"],
    }]
    cards = [
        {"id": "dominant_card", "description": "Shared output mode only; no solver/class-space null.", "dataset": "headline_metrics", "sourceId": "historical_topology_census", "metrics": [{"label": "Direct-dominant datasets", "field": "direct_dominant", "format": "number"}]},
        {"id": "rank_card", "description": "Lower rank is more similar; denominator is 21 pairwise comparisons.", "dataset": "headline_metrics", "sourceId": "bulk_sampling_adjustment", "metrics": [{"label": "HCC JSD rank / 21", "field": "hcc_jsd_rank", "format": "number"}]},
        {"id": "matched_card", "description": "Same-locus pre-VAF coarse endpoint; null mean 39.51%, p=1/5001.", "dataset": "headline_metrics", "sourceId": "matched_topology_pair_analysis", "metrics": [{"label": "Matched coarse agreement", "field": "matched_pre_vaf_coarse", "format": "percent"}]},
        {"id": "same_state_card", "description": "Majority-dominated full-exact external state endpoint.", "dataset": "headline_metrics", "sourceId": "regional_clone_annotations", "metrics": [{"label": "Same possible state", "field": "same_state", "format": "percent"}]},
        {"id": "minor_card", "description": "Minor-signal-focused overlap and informative partition endpoint.", "dataset": "headline_metrics", "sourceId": "clone_bridge_redteam", "metrics": [{"label": "Subclone Jaccard", "field": "subclone_jaccard", "format": "percent"}, {"label": "Informative partition exact", "field": "informative_partition", "format": "percent"}]},
    ]

    charts = [
        {
            "id": "complete_five_class_chart",
            "title": "七個 dataset 的 complete-region 五分類組成",
            "subtitle": "每個 dataset 內正規化為 100%；unresolved 保留在 primary denominator，Direct 仍為 7/7 最大類別。",
            "type": "stackedBar100",
            "dataset": "complete_composition",
            "sourceId": "historical_topology_census",
            "intent": "composition",
            "question": "全樣本是否共享同一個 broad output mode？",
            "rationale": "100% stacked bars reveal both Direct dominance and unresolved/source-dependent shifts without dropping hard regions; no solver/class-space null is implied.",
            "comparisonContext": {"denominator": "complete regions within each dataset", "grain": "dataset × five-class category", "unit": "region share"},
            "encodings": {
                "x": {"field": "sample", "type": "nominal", "label": "Dataset"},
                "y": {"field": "count", "type": "quantitative", "label": "Regions"},
                "color": {"field": "category", "type": "ordinal", "label": "Final shape / unresolved"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
        {
            "id": "hcc_jsd_rank_chart",
            "title": "HCC complete-five-class JSD：調整後仍不是最近 pair",
            "subtitle": "HCC same-cell-line cross-source pair 與 20 個 cross-biological dataset-row pairs 並列；rank 越小越相似。Raw 到 EB 均 rank 9，source-standardized rank 10。",
            "type": "bar",
            "dataset": "hcc_jsd_ranks",
            "sourceId": "bulk_sampling_adjustment",
            "intent": "comparison",
            "question": "sampling/source correction 是否能把 HCC pair 推成最接近的組合？",
            "rationale": "Showing every reviewed adjustment makes the stable non-rank-1 result visually explicit.",
            "comparisonContext": {"baseline": "rank 1 is closest", "denominator": "21 comparisons", "grain": "adjustment method", "unit": "rank"},
            "encodings": {
                "x": {"field": "method", "type": "ordinal", "label": "Adjustment"},
                "y": {"field": "rank", "type": "quantitative", "label": "JSD similarity rank"},
                "tooltip": [
                    {"field": "jsd_median", "type": "quantitative", "label": "Median JSD"},
                    {"field": "rank_q025", "type": "quantitative", "label": "Rank q2.5%"},
                    {"field": "rank_q975", "type": "quantitative", "label": "Rank q97.5%"},
                ],
            },
            "valueFormat": "number",
            "palette": {"kind": "categorical"},
        },
        {
            "id": "regional_state_coverage_chart",
            "title": "七樣本 full-exact-join PyClone regional state coverage",
            "subtitle": "每個樣本以全部 primary regions 為 100%；possible subclone 是區域候選訊號，不是細胞 clone fraction。",
            "type": "stackedBar100",
            "dataset": "regional_state_coverage",
            "sourceId": "regional_clone_annotations",
            "intent": "composition",
            "question": "外部工具在每個樣本能判讀多少 regions，minor signal 佔多少，哪些因 CN/join fail closed？",
            "rationale": "A three-state 100% composition exposes both modeled states and unavailable mass, preventing coverage from being mistaken for biology.",
            "comparisonContext": {"denominator": "all primary regions within each dataset", "grain": "dataset × external-state group", "unit": "region share"},
            "encodings": {
                "x": {"field": "sample", "type": "nominal", "label": "Dataset"},
                "y": {"field": "regions", "type": "quantitative", "label": "Regions"},
                "color": {"field": "status", "type": "nominal", "label": "Conditional regional status"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
    ]

    tables = [
        table_spec("topology_summary_table", "七樣本完整五分類精確數字", "括號為各 dataset complete-region 分母內百分比。", "topology_summary", "historical_topology_census", [
            ("sample", "Dataset", "text"), ("primary_regions", "Primary W", "number"), ("complete_regions", "Complete W", "number"), ("incomplete_regions", "Incomplete", "number"),
            ("single", "Single", "text"), ("sister", "Sister", "text"), ("direct", "Direct", "text"), ("mixed", "Sister＋direct", "text"), ("unresolved", "Unresolved", "text"),
        ]),
        table_spec("rank_table", "Bulk-sampling / source sensitivity 的 HCC JSD rank", "Median rank 與 2.5–97.5% interval；rank 1 才是 closest。", "hcc_jsd_ranks", "bulk_sampling_adjustment", [
            ("method", "Method", "text"), ("jsd_median", "Median JSD", "number"), ("rank", "Median rank", "number"), ("rank_q025", "Rank q2.5%", "number"), ("rank_q975", "Rank q97.5%", "number"), ("rank_1_probability", "P(rank 1)", "percent"), ("replicates", "Replicates", "number"),
        ]),
        table_spec("definitions_table", "C、T、Topo、VAF 與 external cluster 的語意邊界", "這些量不能互相替代；尤其 C 與 PyClone cluster 都不是真實 clone 數。", "definitions", "regional_clone_annotations", [
            ("term", "Term", "text"), ("definition", "This report means", "text"), ("not_equal", "Does not equal", "text"),
        ]),
        table_spec("ladder_table", "47,377 regions 的 C → T → Topo → VAF selection ladder", "T/Topo 分母為 39,885 complete regions；VAF selection 只在 Topo>1 中判斷 most-likely。", "aggregate_ladder", "regional_clone_annotations", [
            ("layer", "Layer", "text"), ("state", "States", "text"), ("counts", "Counts", "text"), ("meaning", "Interpretation", "text"),
        ]),
        table_spec("coverage_table", "七樣本 external regional state 與 coverage", "Full-exact join 為 primary endpoint；modeled cluster count 是 fit-local component count。", "regional_coverage_table", "regional_clone_annotations", [
            ("sample", "Dataset", "text"), ("primary_regions", "Regions W", "number"), ("region_linked_sites", "Region-linked S", "number"), ("mutation_join_fraction", "Mutation join", "percent"),
            ("full_exact_regions", "Full exact W", "number"), ("full_exact_share", "Full exact %", "percent"), ("clonal_like_regions", "Single clonal-like", "number"), ("possible_subclone_regions", "Possible subclone", "number"),
            ("unavailable_regions", "Unavailable", "number"), ("modeled_clusters", "Fit-local clusters", "number"), ("cn_gate", "CN gate", "text"),
        ]),
        table_spec("assignment_table", "Assignment≥0.80 sensitivity：coverage collapse 必須可見", "High-assignment endpoint 非 primary；region 中任一 mutation 未達門檻即 fail closed。", "assignment_sensitivity", "regional_clone_annotations", [
            ("sample", "Dataset", "text"), ("primary_evaluable", "Primary evaluable", "number"), ("primary_evaluable_share", "Primary eval %", "percent"), ("primary_subclone", "Primary subclone W", "number"),
            ("assignment_ge_0_8_evaluable", "≥0.80 evaluable", "number"), ("assignment_ge_0_8_evaluable_share", "≥0.80 eval %", "percent"), ("assignment_ge_0_8_subclone", "≥0.80 subclone W", "number"), ("fail_closed", "Fail closed", "number"),
        ]),
        table_spec("hcc_decomposition_table", "HCC1395 pair：從表面一致拆到 informative minor signal", "每列保留自己的 denominator，避免 majority class 與 selection 造成 inflated agreement。", "hcc_decomposition", "clone_bridge_redteam", [
            ("endpoint", "Endpoint", "text"), ("basis", "Endpoint basis", "text"), ("numerator", "Numerator", "number"), ("denominator", "Denominator", "number"), ("rate", "Rate", "percent"), ("interpretation", "Interpretation", "text"),
        ]),
        table_spec("pyclone_context_table", "Independent PyClone separate fits：global 與 regional mutation context", "Cluster count 是 fit-local component count；agreement 必須與 κ / minor Jaccard 一起讀。", "pyclone_context", "clone_bridge_redteam", [
            ("population", "Population", "text"), ("mutations", "Mutations", "number"), ("fit_local_clusters", "Fit-local clusters A vs B", "text"), ("clonal_state_agreement", "Clonal-state agreement", "percent"),
            ("kappa", "κ", "number"), ("subclonal_intersection", "Subclonal ∩", "number"), ("subclonal_union", "Subclonal ∪", "number"), ("subclonal_jaccard", "Subclonal Jaccard", "percent"), ("interpretation", "Interpretation", "text"),
        ]),
        table_spec("claims_table", "可以說與不能說的結論", "Verdict 是對目前 evidence scope 的 claim audit，不是未來方法上限。", "claims", "bulk_sampling_adjustment", [
            ("claim", "Claim", "text"), ("verdict", "Verdict", "text"), ("evidence", "Evidence", "text"), ("ceiling", "Claim ceiling", "text"),
        ]),
        table_spec("estimand_table", "四個不可混用的 HCC estimands", "Aggregate、matched pre-VAF、matched final-shape、external PyClone 的母體與 null 均不同。", "estimand_comparison", "matched_topology_pair_analysis", [
            ("estimand", "Estimand", "text"), ("population", "Population / grain", "text"), ("result", "Result", "text"), ("null_status", "Null / truth status", "text"), ("interpretation", "Interpretation", "text"),
        ]),
        table_spec("informative_regions_table", "HCC pair 的 34 個 all-joined informative both-multicluster regions", "Common-joined-mutation partition endpoint；21 exact、13 different，可含 incomplete full-exact join。", "informative_regions", "regional_clone_annotations", [
            ("chrom", "Chr", "text"), ("region", "Region", "text"), ("shape_hcc", "Shape HCC", "text"), ("shape_dorado", "Shape DORADO", "text"), ("C_hcc", "C HCC", "number"), ("C_dorado", "C DORADO", "number"),
            ("state_hcc", "External state HCC", "text"), ("state_dorado", "External state DORADO", "text"), ("common_joined_mutations", "Common joined", "number"),
            ("subclonal_intersection", "Subclonal ∩", "number"), ("subclonal_union", "Subclonal ∪", "number"), ("partition", "Partition", "text"), ("jaccard_defined", "Jaccard defined", "text"), ("subclone_jaccard", "Subclone Jaccard", "number"),
        ], default_sort=("chrom", "asc")),
        table_spec("gene_drug_summary_table", "HCC 5,720 regions 的 gene / cancer / drug annotation context", "Annotation coverage 只提供 biological prioritization；不是 topology accuracy 或 treatment evidence。", "gene_drug_summary", "regional_clone_annotations", [
            ("annotation", "Annotation", "text"), ("flagged_regions", "Flagged W", "number"), ("coarse_agree", "Coarse agree", "number"), ("either_subclone", "Either subclone", "number"), ("both_subclone", "Both subclone", "number"),
        ]),
        table_spec("gene_drug_examples_table", "15 個 CGC / approved-antineoplastic-context 且任一側有 possible-subclone signal 的 regions", "DGIdb 是 gene–drug claim annotation；不代表 region 內 variant 可用藥，也不驗證 topology。", "gene_drug_examples", "regional_clone_annotations", [
            ("region", "Region", "text"), ("body_genes", "Body genes", "text"), ("cgc_symbols", "CGC", "text"), ("dgidb_targets", "DGIdb targets", "text"), ("drug_claim_names", "Drug claim names", "text"),
            ("C_hcc", "C HCC", "number"), ("C_dorado", "C DORADO", "number"), ("assignment_min_hcc", "Assign min HCC", "number"), ("assignment_min_dorado", "Assign min DORADO", "number"),
            ("clp_all", "CLP any", "number"), ("clp_confirmed", "CLP confirmed", "number"), ("state_hcc", "HCC state", "text"), ("state_dorado", "DORADO state", "text"), ("partition", "Partition context", "text"),
        ]),
        table_spec("qa_table", "可追溯 QA 與 hash anchors", "產生器會在 schema、population、conservation 或 receipt hash 改變時停止。", "qa_receipts", "regional_clone_annotations", [
            ("artifact", "Artifact", "text"), ("status", "Status", "text"), ("checks", "Checks", "text"), ("population", "Population", "text"), ("hash_anchor", "Summary SHA-256", "text"),
        ]),
    ]

    blocks = [
        {"id": "title", "type": "markdown", "body": "# 全樣本結構一致性：Bulk Sampling 與區域 Clone 候選標記"},
        {"id": "framework", "type": "markdown", "body": "**敘述框架：SCQA + evidence ladder。** 先回答 broad structure 是否重現，再拆解 sampling-adjusted composition、VAF-ranked tree candidate 與 external subclone model。"},
        {
            "id": "partial_historical_warning",
            "type": "markdown",
            "sourceId": "historical_topology_census",
            "body": "**PARTIAL / HISTORICAL LAYERED-v2 — Internal engineering: PASS_WITH_CAVEATS；Scientific / external release: NO-GO。** chr1–22 與 7 dataset rows 已全量納入，但結構層是 historical layered-v2 engineering snapshot。它只描述 VAF-selected rooted-unlabeled mutation-state graph-pattern composition；不等於 clean layered-v3、真實 clone/subclone、biological equivalence 或已驗證的演化樹。External PyClone 亦是 bulk sequencing 下、CN/purity/model-conditional 的 putative clusters 與 cellular prevalence，不是 tree/ancestry truth。",
        },
        {
            "id": "technical_summary",
            "type": "markdown",
            "body": (
                "## 結論先講：matched coarse 訊號部分跨來源再現，但方法的科學有效性仍是 NO-GO\n\n"
                f"- **全七 dataset：** Direct 是 7/7 complete five-class 的最大類別，佔 {pct(min(direct_complete_shares))}–{pct(max(direct_complete_shares))}；這只是共同輸出模式。沒有 solver/class-space/eligibility null，不能單獨當有效性證據。\n"
                "- **最強正向 endpoint：** HCC same-locus pre-VAF coarse agreement=3,969/5,720=69.39%（95% CI 68.18%–70.57%），高於 chromosome-preserving null 39.51%，p=1/5001；支持 partial same-cell-line cross-source reproducibility，不是 final-shape 或 exact-tree accuracy。\n"
                f"- **Matched final shape：** 4,243/5,720=74.18%，composition JSD={matched['final_jsd']:.4f}；尚無對應 permutation null，不能借用 pre-VAF null。\n"
                f"- **HCC aggregate composition：** raw、Dirichlet、rarefaction、5/10 Mb block bootstrap、EB 的 JSD rank 都是 9/21；共有 {matched['closer_pairs']}/20 dataset-row comparisons 更近，涵蓋 {matched['closer_biological_pairs']} 個 unique biological-ID pairs，所以不能稱完整 fingerprint 高度相似。\n"
                "- **HCC external possible state：** 4,296/4,438=96.80% same-state，但 4,267/4,296=99.33% 是雙側皆 single modeled clonal-like；高值主要是 majority class。\n"
                "- **Minor/subclone focus：** full-exact possible-subclone region Jaccard=40/172=23.26%；另一路徑 all-joined common-mutation both-multicluster partition exact=21/34=61.76%。兩個 denominator 不可混用；都只支持部分再現。\n"
                "- **VAF 的角色：** complete regions 中 Topo>1=17,909；其中 15,063 由各位點 VAF 排序得到 most-likely shape，2,846 unresolved。這是推測選擇，不能再拿同一 VAF 當獨立 accuracy truth。\n"
                "- **外部工具 coverage：** PyClone-VI 能對 HCC pair、H1437、H2009、HCC1954 建模；COLO829/HCC1937 缺 reviewed allele-specific CN，全部 fail closed。DORADO 沿用 HCC CN，只是 sensitivity。\n"
                "- **最終主張：** internal diagnostic 可用於定位結構訊號與跨來源部分再現；fresh layered-v3 7/7 root `_SUCCESS` 與 orthogonal truth 前，scientific/external release 維持 NO-GO。"
            ),
        },
        {"id": "headline", "type": "metric-strip", "cardIds": ["dominant_card", "rank_card", "matched_card", "same_state_card", "minor_card"]},
        {"id": "claims_heading", "type": "markdown", "body": "## 一頁判讀：matched coarse 有 null 支持；其他 estimands 不能借用這個證據"},
        {"id": "claims_block", "type": "table", "tableId": "claims_table"},
        {"id": "estimand_intro", "type": "markdown", "sourceId": "matched_topology_pair_analysis", "body": "### 四個 endpoint 不是同一件事\n\nAggregate final-shape composition、matched pre-VAF coarse、matched VAF-selected final shape、external PyClone state/partition 使用不同母體、分類與 null。數字方向可並列，但不能互相當作 validation truth。"},
        {"id": "estimand_block", "type": "table", "tableId": "estimand_table"},
        {
            "id": "topology_heading",
            "type": "markdown",
            "sourceId": "historical_topology_census",
            "body": "## 1. 全樣本粗結構：Direct-dominant 是共同輸出模式，不是 standalone validity evidence\n\nDirect 在七個 dataset 都最大，但目前沒有 solver/class-space/eligibility null，無法排除演算法或可評估區域的 base-rate bias。HCC1395/HCC1395_DORADO 若只看 resolved rows，Direct 是 73.20%/72.26%；放回完整分母後是 71.70%/65.11%，且 unresolved 是 2.05%/9.90%。這張 aggregate 圖描述輸出分佈，不能證明同位點或同 clone evolution。",
        },
        {"id": "topology_chart", "type": "chart", "chartId": "complete_five_class_chart"},
        {"id": "topology_note", "type": "markdown", "body": "圖的每根柱都包含 unresolved，因此不會因排除難解 regions 而人工拉高已解析類別的一致。類別比例是 region graph-pattern composition，不是 clone 細胞比例。"},
        {"id": "topology_table_block", "type": "table", "tableId": "topology_summary_table"},
        {
            "id": "rank_heading",
            "type": "markdown",
            "sourceId": "bulk_sampling_adjustment",
            "body": "## 2. Bulk sampling 校正：HCC same-cell-line cross-source pair 的相對位置穩定，但不是最近鄰\n\n如果差異主要只是區域抽樣、有限樣本量或來源 mix，校正後 HCC 應該向 rank 1 大幅移動。實際上六種 primary/sensitivity adjustment 都停在 rank 9，source-standardization 甚至是 rank 10；這排除了『只因抽樣所以看起來不同』這個過度簡化解釋。",
        },
        {"id": "rank_chart", "type": "chart", "chartId": "hcc_jsd_rank_chart"},
        {"id": "rank_note", "type": "markdown", "body": f"5/10 Mb block bootstrap 的 rank interval 可到 10，但 median 仍為 9；所有 stochastic endpoints 的 P(rank 1)=0。Rank 9/21 具體表示 {matched['closer_pairs']}/20 cross-biological dataset-row comparisons 更近，涵蓋 {matched['closer_biological_pairs']} 個 unique biological-ID pair combinations。Source-standardized 是敏感度分析，不能視為可識別的 causal correction。"},
        {"id": "rank_table_block", "type": "table", "tableId": "rank_table"},
        {
            "id": "semantics_heading",
            "type": "markdown",
            "sourceId": "regional_clone_annotations",
            "body": "## 3. 從 reads 到 external subclone：五個層次不能混稱\n\n`C_read_groups` 是 reads 支持的 mutation-state groups；`T` 是 mutation-labeled candidate forest；`Topo` 是無標籤 shape；VAF 只在候選中選 most-likely；PyClone cluster 則是另一套 bulk VAF/allele-count＋CN/purity 模型的 putative components。它們彼此可以比較訊號，但數值不能直接視為相同的 clone count。尤其 **Topo=1 的 21,976 regions 並不等於 T=1**：其中只有 10,832 個 T=1，另有 11,144 個是 T>1 但候選樹共享同一無標籤 shape；所以拓撲唯一不代表 mutation-labeled 樹唯一。",
        },
        {"id": "definitions_block", "type": "table", "tableId": "definitions_table"},
        {"id": "ladder_block", "type": "table", "tableId": "ladder_table"},
        {
            "id": "vaf_explanation",
            "type": "markdown",
            "body": "**為何可以在多顆候選樹中說『most likely』？** 每個位點的 VAF 提供相對豐度排序訊號，候選 T/shape 可依與這些排序的一致性評分；因此 15,063 個 Topo>1 regions 得到 VAF-ranked final shape。這是模型內推測，不是觀察真實祖先關係：CN/LOH、purity、subclone mixture、allelic imbalance 與 sampling noise 都可能讓 VAF order 偏離真實 lineage order。",
        },
        {
            "id": "external_heading",
            "type": "markdown",
            "sourceId": "regional_clone_annotations",
            "body": "## 4. 其他軟體的 subclone 結果：條件式模型顯示 clonal-majority，minor signal 與 coverage 高度依樣本而變\n\nPyClone-VI 論文將輸出定位為 bulk sequencing 下、校正 copy number 與 normal contamination 後的 putative mutation clusters / cellular prevalence inference；它不是 phylogeny builder，也不是本方法的獨立真值。下圖將 unavailable 一起納入，避免把只有少量可評估 regions 的樣本誤稱為高一致。",
        },
        {"id": "external_chart", "type": "chart", "chartId": "regional_state_coverage_chart"},
        {"id": "external_note", "type": "markdown", "body": "HCC1954 的 possible-subclone region share 明顯較高；H1437/H2009 居中；HCC pair 較低。但這是『有 subclonal-like PyClone cluster 進入該 region』的區域比例，不是 subclone 細胞比例，也不能跨 fit 直接比 cluster ID。COLO829/HCC1937 的 100% unavailable 是正確的 fail-closed 結果。"},
        {"id": "coverage_block", "type": "table", "tableId": "coverage_table"},
        {"id": "assignment_block", "type": "table", "tableId": "assignment_table"},
        {"id": "pyclone_context_intro", "type": "markdown", "body": "### Global 與 regional PyClone mutation endpoints 說同一件事：clonal majority 很像，minor set 只部分再現\n\nSeparate fits 得到 4 vs 3 個 fit-local clusters；數字本身不是 4 vs 3 個真實 clones。Global clonal-state agreement 98.51% 但 κ=0.544、subclonal mutation Jaccard=277/727=38.10%；固定區域 all-joined subset agreement 98.63% 但 κ=0.379、Jaccard=62/259=23.94%。"},
        {"id": "pyclone_context_block", "type": "table", "tableId": "pyclone_context_table"},
        {
            "id": "hcc_heading",
            "type": "markdown",
            "sourceId": "clone_bridge_redteam",
            "body": f"## 5. HCC same-cell-line cross-source pair：96.80% same-state 不能單獨證明 subclone 一致\n\nHCC1395 與 HCC1395_DORADO 是同 cell line 的不同來源資料，但尚未證明 aliquot、passage、library 與 preprocessing 全部相同，因此不是嚴格 technical replicate，也不是 biological replicate。同狀態 4,296/4,438 看似很高，但其中 4,267 是兩側都只有 single modeled clonal-like state。把 clonal majority 拿掉後，172 個任一側有 possible-subclone signal 的 full-exact regions 只有 40 個兩側都有。另一路徑的 21/34 則是較寬鬆 all-joined common-mutation partition endpoint，可包含 incomplete full-exact join；它不是 4,438 primary state endpoint 的子分母。語意 hotfix 後，{pair['jaccard_undefined_union_zero']:,} 個 subclonal union=0 rows 的 Jaccard 明確為 NA，而不是假設成 1。兩個可判讀 endpoint 共同顯示 minor structure 只有部分跨來源再現。",
        },
        {"id": "hcc_decomposition_block", "type": "table", "tableId": "hcc_decomposition_table"},
        {"id": "representative_heading", "type": "markdown", "body": "### 34 個 all-joined informative partition regions（完整列出，未抽樣）\n\n此表是 both-multicluster all-joined common-mutation denominator，可含 incomplete full-exact join，不能和 4,438 full-exact state denominator混用。完整資料請見 [47,377-row region annotation TSV.GZ](research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1/region_possible_clone_annotations.tsv.gz)。"},
        {"id": "representative_block", "type": "table", "tableId": "informative_regions_table"},
        {"id": "gene_drug_heading", "type": "markdown", "body": "## 6. 癌症基因與藥物資料：可做 region prioritization，不能做 accuracy validation\n\n5,720 pair rows 已有 GENCODE、CGC、DGIdb approved-antineoplastic 與 Cancer Licence Project flags。CGC any=275 regions、approved-antineoplastic gene context=464、CLP any=333，但 **CLP confirmed somatic=0**；因此這些只是座標／基因 annotation。特別是表中醒目的 FLI1、LRP1B、CBLC regions，兩來源的 `C_read_groups` 都是 0，且至少一側最小 assignment probability <0.80，沒有 InterSubMod read-level grouping 支持。DGIdb 的藥名代表 database claim context，不代表該 region variant 可用藥，更不能因落在癌症基因就證明 tree 正確。"},
        {"id": "gene_drug_summary_block", "type": "table", "tableId": "gene_drug_summary_table"},
        {"id": "gene_drug_examples_block", "type": "table", "tableId": "gene_drug_examples_table"},
        {
            "id": "methods",
            "type": "markdown",
            "body": (
                "## 7. 方法與判讀條件\n\n"
                "1. **結構 population：** GRCh38 chr1–22、7 dataset rows、47,377 primary regions；39,885 complete regions進 five-class/T/Topo分析。\n"
                "2. **Composition estimand：** complete five-class（Single、Sister、Direct、Sister＋Direct、Unresolved）；HCC same-cell-line cross-source pair 與 20 cross-biological dataset-row pairs 比較 JSD。\n"
                "3. **Sampling sensitivity：** Dirichlet posterior 5,000 replicates、without-replacement rarefaction to n=2,557、共同 5/10 Mb genomic-block bootstrap各5,000 replicates、leave-HCC-out empirical Bayes。\n"
                "4. **Matched topology estimand：** 5,720 個 exact-coordinate complete-both regions 的 pre-VAF coarse agreement 使用 chromosome-preserving 5,000-permutation null；VAF-selected final-shape agreement 另列為 descriptive sensitivity，沒有對應 null，不能借用前者 p-value。\n"
                "5. **External state：** PyClone CP≥0.90 定為 modeled clonal-like，CP<0.90 定為 modeled subclonal-like；必須 full exact mutation join。\n"
                "6. **High assignment sensitivity：** cluster_assignment_prob≥0.80；任一 mutation 未達門檻則整個 region unavailable，避免 cherry-pick。\n"
                "7. **CN gate：** COLO829/HCC1937 無 reviewed allele-specific CN，禁止以 CN=2 代替；DORADO 沿用 HCC1395 CN 只作 sensitivity。\n"
                "8. **Cross-source comparison：** 分開 fits 的 cluster ID 不直接對齊；使用 clonal/subclonal state、subclone-set Jaccard 與 label-invariant region partition。"
            ),
        },
        {
            "id": "limitations",
            "type": "markdown",
            "sourceId": "pyclone_vi_paper",
            "body": (
                "## 8. 限制：為何不能把現有數字升級成真實演化樹證明\n\n"
                "- **Bulk non-identifiability：** 不同 clone mixture、CN、purity 與 mutation multiplicity 可產生相似 VAF/CP。\n"
                "- **同 cell line 不等於嚴格 replicate：** HCC1395/HCC1395_DORADO 是不同來源資料；aliquot、passage、library 與 preprocessing 同一性未被證明，因此只能稱 same-cell-line cross-source pair。\n"
                "- **Direct base-rate 未校正：** 7/7 Direct-dominant 可能受 solver、候選 class space 或 eligibility 分佈影響；沒有相應 null 前只是共同輸出模式。\n"
                "- **Estimand 不可互借：** aggregate composition、matched pre-VAF、matched VAF-selected final shape、PyClone state/partition 的 population、模型與 null 不同。\n"
                "- **VAF circularity：** VAF 用來選最可能 T 後，不能用同一 VAF 當獨立 validation。\n"
                "- **Majority inflation：** same-state 96.80% 由 4,267 both-single-clonal-like rows 主導。\n"
                "- **Selection collapse：** assignment≥0.80 只剩 2,370/5,720 HCC regions 可評估且全同狀態，這是退化 selection，不是完美再現。\n"
                "- **CN gate/source dependence：** 兩樣本 blocked；DORADO CN 非 source-independent。\n"
                "- **Historical release：** clean layered-v3 尚未取代本 snapshot；7/7 sample roots 尚無 fresh `_SUCCESS` release gate。\n"
                "- **No orthogonal tree truth：** 沒有 single-cell lineage、spike-in synthetic tree 或多區域時間序列做 per-region ancestry ground truth。"
            ),
        },
        {
            "id": "next_steps",
            "type": "markdown",
            "body": (
                "## 9. 最短補強路徑\n\n"
                "1. 建立 solver/class-space/eligibility null，先檢查 7/7 Direct-dominant 是否超出方法本身的 output base rate。\n"
                "2. 取得 HCC1395_DORADO source-specific allele-specific CN，重跑 separate PyClone fits；先驗證 minor-cluster membership，不先驗證 tree。\n"
                "3. 對 34 個 informative regions 做 BAM/read-level人工 audit，確認 site set、HP mapping、depth 與 VAF order；將 21 exact/13 different逐區域解釋。\n"
                "4. 建立 within-source split-half ceiling：同一來源隨機分 reads，估計方法在沒有技術切換時可達的 subclone Jaccard/partition exact。\n"
                "5. 取得 fresh clean layered-v3 **7/7 sample-root `_SUCCESS`** 後，以相同固定 denominators 重跑，禁止只比較 resolved subset。\n"
                "6. 若要聲稱 accuracy，加入 synthetic mixtures 或 single-cell/orthogonal lineage truth，預先定義 exact/induced subtree scoring 與 coverage gate。"
            ),
        },
        {"id": "qa_heading", "type": "markdown", "body": "## 10. QA、來源與可重現性\n\n本報告的 7-sample、47,377-region 與 5,720-pair population 都有 conservation checks；bulk-sampling rerun 為 81/81 PASS，regional semantic-hotfix snapshot 為 32/32 PASS，matched pre-VAF endpoint 另有 5,000 次 chromosome-preserving permutation。產生器也核對 region output receipt 的 SHA-256 與 size。"},
        {"id": "qa_block", "type": "table", "tableId": "qa_table"},
        {"id": "footer", "type": "markdown", "body": "**決策 gate：Internal engineering = PASS_WITH_CAVEATS；Scientific / external release = NO-GO。** 本報告只支持 historical layered-v2 結構訊號診斷、same-locus coarse topology 的跨來源部分再現，以及 conditional PyClone possible regional state；不可宣稱真實 clone identity、clone count、ancestry、唯一 tree 或 method accuracy 已被證明。升級至少需要 fresh clean layered-v3 7/7 sample-root `_SUCCESS`、solver/class-space null 與 orthogonal truth。"},
    ]

    snapshot = {
        "version": 1,
        "generatedAt": generated_at,
        "status": "partial",
        "datasets": {
            "headline_metrics": headline_metrics,
            "complete_composition": composition_rows,
            "topology_summary": topology_table,
            "hcc_jsd_ranks": rank_rows,
            "regional_state_coverage": coverage_chart,
            "regional_coverage_table": coverage_table,
            "assignment_sensitivity": assignment_table,
            "aggregate_ladder": aggregate_ladder,
            "definitions": definitions,
            "hcc_decomposition": hcc_decomposition,
            "pyclone_context": pyclone_context,
            "claims": claims,
            "estimand_comparison": estimand_comparison,
            "informative_regions": informative_rows,
            "gene_drug_summary": gene_drug_summary,
            "gene_drug_examples": gene_drug_examples,
            "qa_receipts": qa_rows,
        },
    }
    manifest = {
        "version": 1,
        "surface": "report",
        "title": "全樣本結構一致性：Bulk Sampling 與區域 Clone 候選標記",
        "description": "七 dataset historical layered-v2 complete five-class composition、HCC same-locus pre-VAF permutation endpoint、bulk-sampling/source adjustment、C/T/Topo/VAF ladder 與 PyClone-VI regional possible-state red-team report；scientific release NO-GO。",
        "generatedAt": generated_at,
        "cards": cards,
        "charts": charts,
        "tables": tables,
        "sources": sources,
        "blocks": blocks,
    }
    return {"surface": "report", "manifest": manifest, "snapshot": snapshot, "sources": sources}


def main() -> int:
    args = parse_args()
    data = normalize_inputs(args)
    artifact = build_artifact(data)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "status": "PASS",
        "output": str(args.output.resolve()),
        "snapshot_status": artifact["snapshot"]["status"],
        "datasets": len(artifact["snapshot"]["datasets"]),
        "blocks": len(artifact["manifest"]["blocks"]),
        "charts": len(artifact["manifest"]["charts"]),
        "tables": len(artifact["manifest"]["tables"]),
        "sources": len(artifact["sources"]),
        "representative_regions": len(artifact["snapshot"]["datasets"]["informative_regions"]),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
