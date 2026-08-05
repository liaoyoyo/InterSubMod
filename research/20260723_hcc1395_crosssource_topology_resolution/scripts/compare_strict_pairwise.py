#!/usr/bin/env python3
"""Compare phase-invariant exact-PS strict read-linkage sets across seven datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
CHROMS = tuple(f"chr{i}" for i in range(1, 23))
TARGET_PAIR = frozenset(("HCC1395", "HCC1395_DORADO"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all7-root", type=Path, required=True)
    parser.add_argument("--hcc-extraction-root", type=Path, required=True)
    parser.add_argument("--hcc-strict-root", type=Path, required=True)
    parser.add_argument("--summary-db", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    return parser.parse_args()


def read_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def allele_key(row: dict[str, str]) -> str:
    return f"{row['chrom']}:{row['pos1']}:{row['ref']}:{row['alt']}"


def set_metrics(left: set[Any], right: set[Any]) -> dict[str, float | int]:
    intersection = len(left & right)
    union = len(left | right)
    return {
        "left_n": len(left),
        "right_n": len(right),
        "intersection_n": intersection,
        "union_n": union,
        "jaccard": intersection / union if union else 1.0,
        "left_recall": intersection / len(left) if left else 1.0,
        "right_recall": intersection / len(right) if right else 1.0,
        "overlap_coefficient": intersection / min(len(left), len(right))
        if left and right
        else 1.0,
    }


def co_membership_pairs(
    components: set[tuple[str, ...]], allowed_sites: set[str] | None = None
) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for component in components:
        members = (
            sorted(set(component) & allowed_sites)
            if allowed_sites is not None
            else list(component)
        )
        pairs.update(itertools.combinations(members, 2))
    return pairs


def project_edges(
    edges: set[tuple[str, str]], allowed_sites: set[str]
) -> set[tuple[str, str]]:
    return {
        edge for edge in edges if edge[0] in allowed_sites and edge[1] in allowed_sites
    }


def source_paths(
    dataset: str,
    chrom: str,
    all7_root: Path,
    hcc_extraction_root: Path,
    hcc_strict_root: Path,
) -> tuple[Path, Path, Path, Path]:
    if dataset == "HCC1395":
        extraction = hcc_extraction_root / "chromosomes" / chrom / "extraction"
        strict = hcc_strict_root / "chromosomes" / chrom
    else:
        sample_root = all7_root / "samples" / dataset
        extraction = sample_root / "chromosomes" / chrom / "extraction"
        strict = sample_root / "strict_regions_v1" / "chromosomes" / chrom
    return (
        extraction / f"{dataset}.{chrom}.site_catalog.tsv.gz",
        strict / f"{dataset}.{chrom}.site_component_membership.tsv.gz",
        strict / f"{dataset}.{chrom}.endpoint_edges.tsv.gz",
        strict / f"{dataset}.{chrom}.components.tsv.gz",
    )


def load_expected(summary_db: Path) -> dict[str, dict[str, int]]:
    connection = sqlite3.connect(summary_db)
    connection.row_factory = sqlite3.Row
    rows = connection.execute(
        """
        SELECT dataset, candidate_loci_S, active_unique_loci,
               unique_loci_in_any_W, tree_eligible_W, retained_direct_edges
        FROM dataset_summary
        """
    ).fetchall()
    connection.close()
    return {
        row["dataset"]: {
            key: int(row[key])
            for key in (
                "candidate_loci_S",
                "active_unique_loci",
                "unique_loci_in_any_W",
                "tree_eligible_W",
                "retained_direct_edges",
            )
        }
        for row in rows
    }


def load_dataset(
    dataset: str,
    all7_root: Path,
    hcc_extraction_root: Path,
    hcc_strict_root: Path,
) -> dict[str, Any]:
    candidate_sites: set[str] = set()
    active_sites: set[str] = set()
    w_sites: set[str] = set()
    primary_edges: set[tuple[str, str]] = set()
    alt_informative_edges: set[tuple[str, str]] = set()
    aa_edges: set[tuple[str, str]] = set()
    component_members: dict[str, set[str]] = defaultdict(set)
    site_maps: dict[str, dict[str, str]] = {}
    primary_edge_rows = 0
    component_rows = 0
    files_read = 0

    for chrom in CHROMS:
        site_path, membership_path, edge_path, component_path = source_paths(
            dataset, chrom, all7_root, hcc_extraction_root, hcc_strict_root
        )
        missing = [
            str(path)
            for path in (site_path, membership_path, edge_path, component_path)
            if not path.is_file()
        ]
        if missing:
            raise FileNotFoundError(f"{dataset}/{chrom} missing: {missing}")

        index_to_allele: dict[str, str] = {}
        for row in read_tsv_gz(site_path):
            key = allele_key(row)
            index = row["site_index"]
            if index in index_to_allele and index_to_allele[index] != key:
                raise ValueError(f"{dataset}/{chrom}: site_index collision {index}")
            index_to_allele[index] = key
            candidate_sites.add(key)
        site_maps[chrom] = index_to_allele
        files_read += 1

        for row in read_tsv_gz(membership_path):
            if row["threshold"] != "3":
                continue
            try:
                key = index_to_allele[row["site_index"]]
            except KeyError as exc:
                raise ValueError(
                    f"{dataset}/{chrom}: membership site missing from catalog: {exc}"
                ) from exc
            active_sites.add(key)
            if row["tree_eligible"].lower() == "true":
                w_sites.add(key)
                component_members[row["component_id"]].add(key)
        files_read += 1

        for row in read_tsv_gz(edge_path):
            if row["passes_primary_threshold"].lower() != "true":
                continue
            left = index_to_allele[row["site_i_index"]]
            right = index_to_allele[row["site_j_index"]]
            edge = tuple(sorted((left, right)))
            primary_edges.add(edge)
            primary_edge_rows += 1
            if sum(int(row[key]) for key in ("support_RA", "support_AR", "support_AA")) > 0:
                alt_informative_edges.add(edge)
            if int(row["support_AA"]) > 0:
                aa_edges.add(edge)
        files_read += 1

        for row in read_tsv_gz(component_path):
            if row["threshold"] == "3" and row["tree_eligible"].lower() == "true":
                component_rows += 1
        files_read += 1

    exact_components = {tuple(sorted(members)) for members in component_members.values()}
    if len(component_members) != component_rows:
        raise ValueError(
            f"{dataset}: membership components={len(component_members)} "
            f"but component rows={component_rows}"
        )

    return {
        "sets": {
            "candidate_sites": candidate_sites,
            "active_sites": active_sites,
            "w_sites": w_sites,
            "primary_edges": primary_edges,
            "alt_informative_edges": alt_informative_edges,
            "aa_edges": aa_edges,
            "exact_components": exact_components,
        },
        "counts": {
            "candidate_sites": len(candidate_sites),
            "active_sites": len(active_sites),
            "w_sites": len(w_sites),
            "primary_edge_rows": primary_edge_rows,
            "primary_edge_unique": len(primary_edges),
            "alt_informative_edge_unique": len(alt_informative_edges),
            "aa_edge_unique": len(aa_edges),
            "tree_eligible_component_rows": component_rows,
            "exact_component_unique": len(exact_components),
            "input_files_read": files_read,
        },
    }


def main() -> None:
    args = parse_args()
    expected = load_expected(args.summary_db)
    loaded: dict[str, dict[str, Any]] = {}
    validations: list[dict[str, Any]] = []

    for dataset in DATASETS:
        loaded[dataset] = load_dataset(
            dataset,
            args.all7_root,
            args.hcc_extraction_root,
            args.hcc_strict_root,
        )
        counts = loaded[dataset]["counts"]
        checks = {
            "candidate_loci_S": counts["candidate_sites"]
            == expected[dataset]["candidate_loci_S"],
            "active_unique_loci": counts["active_sites"]
            == expected[dataset]["active_unique_loci"],
            "unique_loci_in_any_W": counts["w_sites"]
            == expected[dataset]["unique_loci_in_any_W"],
            "tree_eligible_W": counts["tree_eligible_component_rows"]
            == expected[dataset]["tree_eligible_W"],
            "retained_direct_edges": counts["primary_edge_rows"]
            == expected[dataset]["retained_direct_edges"],
            "files_88": counts["input_files_read"] == 88,
        }
        validations.append(
            {
                "dataset": dataset,
                "checks": checks,
                "all_pass": all(checks.values()),
            }
        )

    metric_names = tuple(next(iter(loaded.values()))["sets"].keys())
    pairwise: list[dict[str, Any]] = []
    for left, right in itertools.combinations(DATASETS, 2):
        record: dict[str, Any] = {
            "left": left,
            "right": right,
            "is_target_pair": frozenset((left, right)) == TARGET_PAIR,
        }
        for metric in metric_names:
            record[metric] = set_metrics(
                loaded[left]["sets"][metric], loaded[right]["sets"][metric]
            )
        candidate_intersection = record["candidate_sites"]["intersection_n"]
        record["joint_w_coverage_of_shared_candidates"] = (
            record["w_sites"]["intersection_n"] / candidate_intersection
            if candidate_intersection
            else 0.0
        )
        pairwise.append(record)

    ranks: dict[str, dict[str, Any]] = {}
    for metric in metric_names:
        ordered = sorted(
            pairwise,
            key=lambda row: (
                row[metric]["jaccard"],
                row[metric]["intersection_n"],
            ),
            reverse=True,
        )
        target_index = next(i for i, row in enumerate(ordered) if row["is_target_pair"])
        target = ordered[target_index]
        next_unrelated = next(row for row in ordered if not row["is_target_pair"])
        ranks[metric] = {
            "target_rank_of_21": target_index + 1,
            "target_jaccard": target[metric]["jaccard"],
            "next_unrelated_pair": [next_unrelated["left"], next_unrelated["right"]],
            "next_unrelated_jaccard": next_unrelated[metric]["jaccard"],
            "target_to_next_ratio": (
                target[metric]["jaccard"] / next_unrelated[metric]["jaccard"]
                if next_unrelated[metric]["jaccard"]
                else None
            ),
        }

    target = next(row for row in pairwise if row["is_target_pair"])
    left_sets = loaded[target["left"]]["sets"]
    right_sets = loaded[target["right"]]["sets"]
    shared_candidates = left_sets["candidate_sites"] & right_sets["candidate_sites"]
    target_projection: dict[str, Any] = {
        "shared_candidate_sites_n": len(shared_candidates)
    }
    for metric in ("active_sites", "w_sites"):
        target_projection[metric] = set_metrics(
            left_sets[metric] & shared_candidates,
            right_sets[metric] & shared_candidates,
        )
    for metric in ("primary_edges", "alt_informative_edges", "aa_edges"):
        target_projection[metric] = set_metrics(
            project_edges(left_sets[metric], shared_candidates),
            project_edges(right_sets[metric], shared_candidates),
        )

    left_co_membership = co_membership_pairs(left_sets["exact_components"])
    right_co_membership = co_membership_pairs(right_sets["exact_components"])
    target["co_membership_pairs"] = set_metrics(
        left_co_membership, right_co_membership
    )
    target_projection["co_membership_pairs"] = set_metrics(
        co_membership_pairs(left_sets["exact_components"], shared_candidates),
        co_membership_pairs(right_sets["exact_components"], shared_candidates),
    )
    target["shared_candidate_projection"] = target_projection

    payload = {
        "schema_version": "strict_pairwise_phase_invariant_v1",
        "scope": {
            "datasets": list(DATASETS),
            "chromosomes": list(CHROMS),
            "pair_count": len(pairwise),
            "strict_threshold": 3,
            "target_pair": sorted(TARGET_PAIR),
        },
        "definitions": {
            "allele_key": "chrom:pos1:ref:alt",
            "phase_invariance": "PS and HP labels are intentionally excluded across datasets",
            "candidate_sites": "all site_catalog rows",
            "active_sites": "unique allele keys in threshold=3 membership",
            "w_sites": "unique allele keys in threshold=3 tree_eligible components",
            "primary_edges": "undirected allele endpoint pairs with passes_primary_threshold=true",
            "alt_informative_edges": "primary edges with RA+AR+AA support > 0",
            "aa_edges": "primary edges with AA support > 0",
            "exact_components": "exact allele-member sets of threshold=3 tree_eligible components",
            "claim_ceiling": (
                "read-linkage observability only; not mutation-state topology, ancestry, "
                "clone count, or a global clone tree"
            ),
        },
        "inputs": {
            "all7_root": str(args.all7_root.resolve()),
            "hcc_extraction_root": str(args.hcc_extraction_root.resolve()),
            "hcc_strict_root": str(args.hcc_strict_root.resolve()),
            "summary_db": str(args.summary_db.resolve()),
        },
        "dataset_counts": {
            dataset: loaded[dataset]["counts"] for dataset in DATASETS
        },
        "validations": validations,
        "all_validations_pass": all(row["all_pass"] for row in validations),
        "target_pair": target,
        "target_ranks": ranks,
        "pairwise": pairwise,
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    columns = ["left", "right", "is_target_pair"]
    for metric in metric_names:
        columns.extend(
            (
                f"{metric}_intersection_n",
                f"{metric}_union_n",
                f"{metric}_jaccard",
                f"{metric}_left_recall",
                f"{metric}_right_recall",
                f"{metric}_overlap_coefficient",
            )
        )
    columns.append("joint_w_coverage_of_shared_candidates")
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in pairwise:
            flat: dict[str, Any] = {
                "left": row["left"],
                "right": row["right"],
                "is_target_pair": row["is_target_pair"],
                "joint_w_coverage_of_shared_candidates": row[
                    "joint_w_coverage_of_shared_candidates"
                ],
            }
            for metric in metric_names:
                for key in (
                    "intersection_n",
                    "union_n",
                    "jaccard",
                    "left_recall",
                    "right_recall",
                    "overlap_coefficient",
                ):
                    flat[f"{metric}_{key}"] = row[metric][key]
            writer.writerow(flat)

    print(
        json.dumps(
            {
                "all_validations_pass": payload["all_validations_pass"],
                "target_pair": target,
                "target_ranks": ranks,
                "output_json": str(args.output_json),
                "output_tsv": str(args.output_tsv),
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
