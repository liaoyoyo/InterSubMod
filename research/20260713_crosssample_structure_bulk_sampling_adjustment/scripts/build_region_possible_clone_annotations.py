#!/usr/bin/env python3
"""Annotate historical regional topology with conditional PyClone-VI states.

The output deliberately keeps read-supported group count C separate from
external PyClone-VI cluster counts.  All clone/subclone labels are conditional
model annotations, not observed biological truth.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


CLONAL_CP_THRESHOLD = 0.90
HIGH_ASSIGNMENT_THRESHOLD = 0.80
EXPECTED_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
CN_READY_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "H1437",
    "H2009",
    "HCC1954",
)
CN_BLOCKED_SAMPLES = ("COLO829", "HCC1937")


def parse_args() -> argparse.Namespace:
    script = Path(__file__).resolve()
    repo = script.parents[3]
    topic = script.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=repo)
    parser.add_argument(
        "--run-root",
        type=Path,
        default=repo
        / "output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2",
    )
    parser.add_argument(
        "--coarse-regions",
        type=Path,
        default=repo
        / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv",
    )
    parser.add_argument(
        "--selected-regions",
        type=Path,
        default=repo
        / "research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_regions.tsv",
    )
    parser.add_argument(
        "--ct-report",
        type=Path,
        default=repo
        / "research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json",
    )
    parser.add_argument(
        "--pair-matches",
        type=Path,
        default=repo
        / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv",
    )
    parser.add_argument(
        "--bridge-patterns",
        type=Path,
        default=repo
        / "research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/region_cluster_patterns.tsv.gz",
    )
    parser.add_argument(
        "--gene-drug-flags",
        type=Path,
        default=repo
        / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/hcc1395_exact_complete_pair_gene_drug_flags.tsv",
    )
    parser.add_argument(
        "--raw-vaf-records",
        type=Path,
        default=repo
        / "research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/raw_vaf_validation_v1/data/raw_vaf_records.tsv.gz",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=topic / "results/region_possible_clone_annotations_v1",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_gzip_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="")
    return raw, zipped, text


def write_tsv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
    gzip_output: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if gzip_output:
        raw, zipped, text = stable_gzip_writer(path)
        try:
            writer = csv.DictWriter(text, delimiter="\t", fieldnames=list(columns), extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
            text.flush()
        finally:
            text.close()
            zipped.close()
            raw.close()
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_groups(sample_dir: Path) -> dict[str, dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        for group in document.get("groups", []):
            key = f"{group['chrom']}:{group['start']}-{group['end']}"
            if key in groups:
                raise RuntimeError(f"Duplicate group: {sample_dir.name} {key}")
            groups[key] = group
    return groups


def c_alt_for_family(group: Mapping[str, Any], family: str) -> int:
    full = ((group.get("populations_by_hp") or {}).get(family) or {})
    return sum("A" in genotype for genotype in full)


def count_bins(values: Iterable[int], cap: int = 6) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for value in values:
        counts[str(value) if value <= cap else f">{cap}"] += 1
    return {key: counts.get(key, 0) for key in [*(str(i) for i in range(cap + 1)), f">{cap}"]}


def pyclone_paths(repo: Path) -> dict[str, Path]:
    root = repo / "research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi"
    return {
        "HCC1395": root / "hcc1395_pair_primary_separate_HCC1395_main/results.tsv.gz",
        "HCC1395_DORADO": root
        / "hcc1395_pair_primary_separate_HCC1395_DORADO_main/results.tsv.gz",
        "H1437": root / "H1437_individual_main/results.tsv.gz",
        "H2009": root / "H2009_individual_main/results.tsv.gz",
        "HCC1954": root / "HCC1954_individual_main/results.tsv.gz",
    }


def load_pyclone(path: Path, expected_sample: str) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    frame = pd.read_csv(path, sep="\t")
    if set(frame["sample_id"].astype(str)) != {expected_sample}:
        raise RuntimeError(f"Unexpected sample_id in {path}")
    parsed = frame["mutation_id"].str.extract(r"^(chr[^:]+):(\d+):([^>]+)>(.+)$")
    if parsed.isna().any().any():
        raise RuntimeError(f"Unparseable mutation_id in {path}")
    frame["chrom"] = parsed[0]
    frame["position"] = parsed[1].astype(int)
    duplicated = int(frame.duplicated(["chrom", "position"], keep=False).sum())
    if duplicated:
        raise RuntimeError(f"Position-only join is ambiguous in {path}: {duplicated} rows")
    records: dict[str, dict[str, Any]] = {}
    for row in frame.to_dict("records"):
        key = str(row["mutation_id"])
        if key in records:
            raise RuntimeError(f"Duplicate exact mutation_id in {path}: {key}")
        records[key] = row
    return records, {
        "rows": int(len(frame)),
        "unique_positions": int(len(records)),
        "position_duplicates": duplicated,
        "clusters": int(frame["cluster_id"].nunique()),
    }


def unavailable_state(reason: str) -> dict[str, Any]:
    state = classify_records([])
    state["status"] = reason
    return state


def classify_records(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if not records:
        return {
            "status": "unavailable_no_joined_mutation",
            "possible_clone_signal": False,
            "possible_subclone_signal": False,
            "joined": 0,
            "cluster_count": 0,
            "clonal_cluster_count": 0,
            "subclonal_cluster_count": 0,
            "clonal_mutations": 0,
            "subclonal_mutations": 0,
            "cp_boundary_uncertain_mutations": 0,
            "cluster_ids": "",
            "subclonal_cluster_ids": "",
            "cp_min": "",
            "cp_max": "",
            "assignment_min": "",
            "assignment_median": "",
        }
    clonal = [row for row in records if float(row["cellular_prevalence"]) >= CLONAL_CP_THRESHOLD]
    subclonal = [row for row in records if float(row["cellular_prevalence"]) < CLONAL_CP_THRESHOLD]
    all_clusters = sorted({int(row["cluster_id"]) for row in records})
    clonal_clusters = sorted({int(row["cluster_id"]) for row in clonal})
    subclonal_clusters = sorted({int(row["cluster_id"]) for row in subclonal})
    if len(records) < 2:
        status = "not_evaluable_lt2_joined_mutations"
    elif len(all_clusters) == 1 and not subclonal:
        status = "single_modeled_clonal_like_cluster_represented"
    elif len(all_clusters) == 1:
        status = "single_modeled_subclonal_like_cluster_represented"
    elif not subclonal:
        status = "multiple_clonal_like_clusters_candidate"
    elif not clonal:
        status = "multiple_subclonal_cluster_candidate"
    elif len(subclonal_clusters) == 1:
        status = "clonal_plus_single_subclonal_cluster_candidate"
    else:
        status = "clonal_plus_multiple_subclonal_clusters_candidate"
    cp = np.array([float(row["cellular_prevalence"]) for row in records], dtype=float)
    assignment = np.array([float(row["cluster_assignment_prob"]) for row in records], dtype=float)
    boundary_uncertain = sum(
        float(row["cellular_prevalence"]) - 1.96 * float(row["cellular_prevalence_std"])
        <= CLONAL_CP_THRESHOLD
        <= float(row["cellular_prevalence"]) + 1.96 * float(row["cellular_prevalence_std"])
        for row in records
    )
    return {
        "status": status,
        "possible_clone_signal": bool(clonal),
        "possible_subclone_signal": bool(subclonal),
        "joined": len(records),
        "cluster_count": len(all_clusters),
        "clonal_cluster_count": len(clonal_clusters),
        "subclonal_cluster_count": len(subclonal_clusters),
        "clonal_mutations": len(clonal),
        "subclonal_mutations": len(subclonal),
        "cp_boundary_uncertain_mutations": boundary_uncertain,
        "cluster_ids": ",".join(map(str, all_clusters)),
        "subclonal_cluster_ids": ",".join(map(str, subclonal_clusters)),
        "cp_min": f"{cp.min():.6f}",
        "cp_max": f"{cp.max():.6f}",
        "assignment_min": f"{assignment.min():.6f}",
        "assignment_median": f"{np.median(assignment):.6f}",
    }


def pair_state(status_a: str, status_b: str, sub_a: bool, sub_b: bool) -> str:
    if status_a.startswith(("unavailable", "not_evaluable")) or status_b.startswith(
        ("unavailable", "not_evaluable")
    ):
        return "not_evaluable_missing_external_cluster"
    if status_a == status_b:
        return "same_possible_state"
    if sub_a and sub_b:
        return "both_subclone_signal_different_state"
    if sub_a:
        return "HCC1395_only_subclone_signal"
    if sub_b:
        return "DORADO_only_subclone_signal"
    return "different_non_subclone_state"


def main() -> int:
    args = parse_args()
    repo = args.repo_root.resolve()
    run_root = args.run_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    input_paths = {
        "coarse_regions": args.coarse_regions.resolve(),
        "selected_regions": args.selected_regions.resolve(),
        "ct_report": args.ct_report.resolve(),
        "pair_matches": args.pair_matches.resolve(),
        "bridge_patterns": args.bridge_patterns.resolve(),
        "gene_drug_flags": args.gene_drug_flags.resolve(),
        "raw_vaf_records": args.raw_vaf_records.resolve(),
    }
    for name, path in input_paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name}: {path}")

    coarse = pd.read_csv(input_paths["coarse_regions"], sep="\t")
    selected = pd.read_csv(input_paths["selected_regions"], sep="\t")
    if coarse.duplicated(["sample", "region"]).any():
        raise RuntimeError("coarse sample+region is not unique")
    if selected.duplicated(["sample", "region"]).any():
        raise RuntimeError("selected sample+region is not unique")
    if tuple(coarse["sample"].drop_duplicates()) != EXPECTED_SAMPLES:
        raise RuntimeError("Unexpected sample order or membership")
    selected_index = {
        (str(row.sample), str(row.region)): row
        for row in selected.itertuples(index=False)
    }

    ct_report = json.loads(input_paths["ct_report"].read_text(encoding="utf-8"))
    ct_by_sample = {row["sample"]: row for row in ct_report["samples"]}
    raw_vaf = pd.read_csv(
        input_paths["raw_vaf_records"],
        sep="\t",
        compression="gzip",
        usecols=["sample", "source", "chrom", "pos", "ref", "alt"],
    )
    raw_vaf = raw_vaf[raw_vaf["source"] == "baseline_caller_pass"].copy()
    if raw_vaf.duplicated(["sample", "chrom", "pos"], keep=False).any():
        duplicated = raw_vaf[raw_vaf.duplicated(["sample", "chrom", "pos"], keep=False)]
        raise RuntimeError(f"Raw exact allele map has ambiguous positions: {len(duplicated)} rows")
    raw_mutation_id = {
        (str(row.sample), str(row.chrom), int(row.pos)): f"{row.chrom}:{int(row.pos)}:{row.ref}>{row.alt}"
        for row in raw_vaf.itertuples(index=False)
    }
    py_paths = pyclone_paths(repo)
    py_indices: dict[str, dict[str, dict[str, Any]]] = {}
    py_qa: dict[str, dict[str, Any]] = {}
    for sample, path in py_paths.items():
        index, qa = load_pyclone(path.resolve(), sample)
        py_indices[sample] = index
        py_qa[sample] = {**qa, "path": str(path.resolve()), "sha256": sha256_file(path.resolve())}

    annotation_rows: list[dict[str, Any]] = []
    c_observed: dict[str, list[int]] = defaultdict(list)
    construction_qa: dict[str, dict[str, int]] = {}
    for sample in EXPECTED_SAMPLES:
        sample_dir = run_root / "samples" / sample
        groups = load_groups(sample_dir)
        layered_path = sample_dir / f"layered_reconstruction_{sample}.json"
        layered = json.loads(layered_path.read_text(encoding="utf-8"))
        primary_by_region: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for unit in layered["detail"]:
            if unit.get("is_primary_lineage"):
                primary_by_region[str(unit["region"])].append(unit)
        sample_coarse = coarse[coarse["sample"] == sample]
        missing_groups = 0
        missing_units = 0
        missing_raw_alleles = 0
        external_index = py_indices.get(sample, {})
        for coarse_row in sample_coarse.itertuples(index=False):
            region = str(coarse_row.region)
            group = groups.get(region)
            units = sorted(primary_by_region.get(region, []), key=lambda row: str(row["family"]))
            if group is None:
                missing_groups += 1
                raise RuntimeError(f"Missing mlhp group {sample} {region}")
            if not units:
                missing_units += 1
                raise RuntimeError(f"Missing primary units {sample} {region}")
            c_by_hp = {str(unit["family"]): c_alt_for_family(group, str(unit["family"])) for unit in units}
            c_total = sum(c_by_hp.values())
            c_observed[sample].append(c_total)
            complete = all(unit.get("analysis_candidate_set_complete") is True for unit in units)
            t_joint = math.prod(int(unit["n_trees"]) for unit in units) if complete else ""
            topo_joint = math.prod(int(unit["n_distinct_shapes_exact"]) for unit in units) if complete else ""
            hidden_total = sum(int(unit["n_hidden"]) for unit in units) if complete else ""
            positions = [int(value) for value in group.get("positions", [])]
            mutation_ids = []
            for position in positions:
                mutation_id = raw_mutation_id.get((sample, str(coarse_row.chrom), position))
                if mutation_id is None:
                    missing_raw_alleles += 1
                    continue
                mutation_ids.append(mutation_id)
            joined_records = [external_index[mutation_id] for mutation_id in mutation_ids if mutation_id in external_index]
            high_records = [
                row
                for row in joined_records
                if float(row["cluster_assignment_prob"]) >= HIGH_ASSIGNMENT_THRESHOLD
            ]
            diagnostic_state = classify_records(joined_records)
            primary_state = (
                diagnostic_state
                if len(mutation_ids) == len(positions) and len(joined_records) == len(positions)
                else unavailable_state("unavailable_incomplete_exact_join")
            )
            high_state = (
                classify_records(high_records)
                if len(mutation_ids) == len(positions)
                and len(joined_records) == len(positions)
                and len(high_records) == len(positions)
                else unavailable_state("unavailable_incomplete_join_or_low_assignment")
            )
            if sample in CN_BLOCKED_SAMPLES:
                diagnostic_state = unavailable_state("unavailable_cn_gate")
                primary_state = unavailable_state("unavailable_cn_gate")
                high_state = unavailable_state("unavailable_cn_gate")
            selected_row = selected_index.get((sample, region))
            if selected_row is not None:
                final_category = str(selected_row.final_shape_category)
                selection_source = str(selected_row.selection_source)
                vaf_exact_top_class = (
                    "" if pd.isna(selected_row.vaf_exact_top_class) else str(selected_row.vaf_exact_top_class)
                )
                vaf_top_direction_state = (
                    "" if pd.isna(selected_row.vaf_top_direction_state) else str(selected_row.vaf_top_direction_state)
                )
                selected_region_shape_id = str(selected_row.selected_region_shape_id)
            elif bool(coarse_row.complete):
                final_category = "unresolved"
                selection_source = "vaf_unresolved_topogt1"
                vaf_exact_top_class = "not_selected_unresolved"
                vaf_top_direction_state = ""
                selected_region_shape_id = ""
            else:
                final_category = "incomplete"
                selection_source = "not_evaluable_incomplete"
                vaf_exact_top_class = ""
                vaf_top_direction_state = ""
                selected_region_shape_id = ""
            row = {
                "sample": sample,
                "biological_id": "HCC1395" if sample == "HCC1395_DORADO" else sample,
                "region": region,
                "chrom": str(coarse_row.chrom),
                "start": int(coarse_row.start),
                "end": int(coarse_row.end),
                "complete": bool(coarse_row.complete),
                "coarse_category": str(coarse_row.coarse_category),
                "final_shape_category": final_category,
                "tree_selection_source": selection_source,
                "vaf_most_likely_topology_status": vaf_exact_top_class,
                "vaf_direction_consistency": vaf_top_direction_state,
                "selected_region_shape_id": selected_region_shape_id,
                "vaf_inferred_final_shape": selection_source == "vaf_resolved_topogt1",
                "primary_hp_units": len(units),
                "C_read_groups_total": c_total,
                "C_read_groups_by_HP": ";".join(f"HP{key}={value}" for key, value in sorted(c_by_hp.items())),
                "T_exact_candidate_forest_count": t_joint,
                "Topo_shape_count": topo_joint,
                "hidden_node_count": hidden_total,
                "mutation_state_nodes_C_plus_hidden": c_total + hidden_total if complete else "",
                "region_site_count": len(positions),
                "region_mutation_ids": ",".join(mutation_ids),
                "external_cn_fit_available": sample in CN_READY_SAMPLES,
                "external_cn_source_independent": "" if sample in CN_BLOCKED_SAMPLES else sample != "HCC1395_DORADO",
                "external_global_fit_cluster_count": py_qa[sample]["clusters"] if sample in py_qa else "",
                "external_joined_mutations": diagnostic_state["joined"],
                "external_joined_mutation_ids": ",".join(str(row["mutation_id"]) for row in joined_records),
                "external_join_fraction": f"{diagnostic_state['joined'] / len(positions):.6f}" if positions else "",
                "diagnostic_possible_state": diagnostic_state["status"],
                "diagnostic_possible_clone_signal": diagnostic_state["possible_clone_signal"],
                "diagnostic_possible_subclone_signal": diagnostic_state["possible_subclone_signal"],
                "external_cluster_count": diagnostic_state["cluster_count"],
                "external_clonal_cluster_count": diagnostic_state["clonal_cluster_count"],
                "external_subclonal_cluster_count": diagnostic_state["subclonal_cluster_count"],
                "external_clonal_mutations": diagnostic_state["clonal_mutations"],
                "external_subclonal_mutations": diagnostic_state["subclonal_mutations"],
                "external_cp_boundary_uncertain_mutations": diagnostic_state["cp_boundary_uncertain_mutations"],
                "external_cluster_ids": diagnostic_state["cluster_ids"],
                "external_subclonal_cluster_ids": diagnostic_state["subclonal_cluster_ids"],
                "possible_clone_signal": primary_state["possible_clone_signal"],
                "possible_subclone_signal": primary_state["possible_subclone_signal"],
                "regional_possible_state": primary_state["status"],
                "external_cp_min": diagnostic_state["cp_min"],
                "external_cp_max": diagnostic_state["cp_max"],
                "external_assignment_min": diagnostic_state["assignment_min"],
                "external_assignment_median": diagnostic_state["assignment_median"],
                "high_assignment_joined_mutations": high_state["joined"],
                "high_assignment_possible_clone_signal": high_state["possible_clone_signal"],
                "high_assignment_possible_subclone_signal": high_state["possible_subclone_signal"],
                "high_assignment_possible_state": high_state["status"],
                "claim_ceiling": "conditional_PyClone_possible_state_not_observed_clone_truth",
            }
            annotation_rows.append(row)
        construction_qa[sample] = {
            "coarse_regions": int(len(sample_coarse)),
            "mlhp_groups": int(len(groups)),
            "primary_regions": int(len(primary_by_region)),
            "missing_groups": missing_groups,
            "missing_primary_units": missing_units,
            "missing_raw_alleles": missing_raw_alleles,
        }

    annotation_rows.sort(key=lambda row: (EXPECTED_SAMPLES.index(str(row["sample"])), str(row["chrom"]), int(row["start"]), str(row["region"])))
    annotation_columns = list(annotation_rows[0].keys())
    annotation_path = output_dir / "region_possible_clone_annotations.tsv.gz"
    write_tsv(annotation_path, annotation_rows, annotation_columns, gzip_output=True)

    summary_rows: list[dict[str, Any]] = []
    for sample in EXPECTED_SAMPLES:
        rows = [row for row in annotation_rows if row["sample"] == sample]
        for endpoint, field in (
            ("all_joined_diagnostic", "diagnostic_possible_state"),
            ("full_exact_join", "regional_possible_state"),
            ("full_exact_assignment_ge_0.8", "high_assignment_possible_state"),
        ):
            counts = Counter(str(row[field]) for row in rows)
            for status, count in sorted(counts.items()):
                summary_rows.append(
                    {
                        "sample": sample,
                        "endpoint": endpoint,
                        "status": status,
                        "regions": count,
                        "sample_regions": len(rows),
                        "share": f"{count / len(rows):.8f}",
                    }
                )
    summary_columns = ["sample", "endpoint", "status", "regions", "sample_regions", "share"]
    summary_path = output_dir / "sample_possible_clone_summary.tsv"
    write_tsv(summary_path, summary_rows, summary_columns)

    pair_matches = pd.read_csv(input_paths["pair_matches"], sep="\t")
    pair_matches = pair_matches[
        (pair_matches["scenario"] == "exact_coordinate")
        & pair_matches["complete_a"].astype(bool)
        & pair_matches["complete_b"].astype(bool)
    ].copy()
    ann_index = {(str(row["sample"]), str(row["region"])): row for row in annotation_rows}
    bridge = pd.read_csv(input_paths["bridge_patterns"], sep="\t", compression="gzip")
    bridge = bridge[bridge["endpoint"] == "all_joined"].copy()
    if bridge.duplicated(["match_id"]).any():
        raise RuntimeError("all_joined bridge pattern match_id is not unique")
    bridge_index = {str(row.match_id): row for row in bridge.itertuples(index=False)}
    gene_drug = pd.read_csv(input_paths["gene_drug_flags"], sep="\t")
    gene_drug = gene_drug[gene_drug["scenario"] == "exact_coordinate"].copy()
    if gene_drug.duplicated(["match_id"]).any():
        raise RuntimeError("exact gene/drug flag match_id is not unique")
    gene_drug_index = {str(row.match_id): row for row in gene_drug.itertuples(index=False)}

    def text_or_blank(value: Any) -> str:
        return "" if pd.isna(value) else str(value)

    pair_rows: list[dict[str, Any]] = []
    pair_common_join_mismatches = 0
    for match in pair_matches.itertuples(index=False):
        a = ann_index[("HCC1395", str(match.region_a))]
        b = ann_index[("HCC1395_DORADO", str(match.region_b))]
        pattern = bridge_index[str(match.match_id)]
        gene = gene_drug_index[str(match.match_id)]
        joined_ids_a = {value for value in str(a["external_joined_mutation_ids"]).split(",") if value}
        joined_ids_b = {value for value in str(b["external_joined_mutation_ids"]).split(",") if value}
        common_joined_ids = joined_ids_a & joined_ids_b
        if len(common_joined_ids) != int(pattern.joined_mutations):
            pair_common_join_mismatches += 1
        subclonal_ids_a = {
            mutation_id
            for mutation_id in common_joined_ids
            if float(py_indices["HCC1395"][mutation_id]["cellular_prevalence"])
            < CLONAL_CP_THRESHOLD
        }
        subclonal_ids_b = {
            mutation_id
            for mutation_id in common_joined_ids
            if float(py_indices["HCC1395_DORADO"][mutation_id]["cellular_prevalence"])
            < CLONAL_CP_THRESHOLD
        }
        subclonal_intersection = len(subclonal_ids_a & subclonal_ids_b)
        subclonal_union = len(subclonal_ids_a | subclonal_ids_b)
        partition_informative = str(pattern.pattern_category) in {
            "both_multicluster_exact_partition",
            "both_multicluster_different_partition",
        }
        partition_vacuity_reason = {
            "both_single_cluster": "both_single_cluster_no_partition_variation",
            "one_single_one_multicluster": "resolution_asymmetric",
            "not_evaluable_lt2_joined_mutations": "lt2_common_joined_mutations",
        }.get(str(pattern.pattern_category), "")
        pair_rows.append(
            {
                "match_id": str(match.match_id),
                "chrom": str(match.chrom),
                "region_a": str(match.region_a),
                "region_b": str(match.region_b),
                "coarse_category_a": str(match.category_a),
                "coarse_category_b": str(match.category_b),
                "coarse_category_agree": bool(match.category_agree),
                "final_shape_category_a": a["final_shape_category"],
                "final_shape_category_b": b["final_shape_category"],
                "C_read_groups_total_a": a["C_read_groups_total"],
                "C_read_groups_total_b": b["C_read_groups_total"],
                "external_joined_mutations_a": a["external_joined_mutations"],
                "external_joined_mutations_b": b["external_joined_mutations"],
                "external_join_fraction_a": a["external_join_fraction"],
                "external_join_fraction_b": b["external_join_fraction"],
                "regional_possible_state_a": a["regional_possible_state"],
                "regional_possible_state_b": b["regional_possible_state"],
                "possible_subclone_signal_a": a["possible_subclone_signal"],
                "possible_subclone_signal_b": b["possible_subclone_signal"],
                "cross_source_possible_state": pair_state(
                    str(a["regional_possible_state"]),
                    str(b["regional_possible_state"]),
                    bool(a["possible_subclone_signal"]),
                    bool(b["possible_subclone_signal"]),
                ),
                "high_assignment_possible_state_a": a["high_assignment_possible_state"],
                "high_assignment_possible_state_b": b["high_assignment_possible_state"],
                "high_assignment_cross_source_state": pair_state(
                    str(a["high_assignment_possible_state"]),
                    str(b["high_assignment_possible_state"]),
                    bool(a["high_assignment_possible_subclone_signal"]),
                    bool(b["high_assignment_possible_subclone_signal"]),
                ),
                "external_partition_evaluable": bool(pattern.evaluable_multilocus),
                "external_partition_informative": partition_informative,
                "external_partition_exact": bool(pattern.partition_exact) if partition_informative else "",
                "external_partition_pattern": str(pattern.pattern_category),
                "external_partition_vacuity_reason": partition_vacuity_reason,
                "external_common_joined_mutations": len(common_joined_ids),
                "external_subclonal_intersection_n": subclonal_intersection,
                "external_subclonal_union_n": subclonal_union,
                "external_subclonal_jaccard_defined": subclonal_union > 0,
                "external_subclonal_jaccard": (
                    f"{subclonal_intersection / subclonal_union:.6f}" if subclonal_union > 0 else ""
                ),
                "external_subclonal_jaccard_vacuity_reason": "" if subclonal_union > 0 else "both_absent",
                "body_gene_symbols": text_or_blank(gene.body_gene_symbols),
                "promoter_gene_symbols": text_or_blank(gene.promoter_gene_symbols),
                "cgc_body_symbols": text_or_blank(gene.cgc_body_symbols),
                "cgc_promoter_symbols": text_or_blank(gene.cgc_promoter_symbols),
                "dgidb_approved_antineoplastic_body_symbols": text_or_blank(
                    gene.dgidb_approved_antineoplastic_body_symbols
                ),
                "dgidb_approved_antineoplastic_promoter_symbols": text_or_blank(
                    gene.dgidb_approved_antineoplastic_promoter_symbols
                ),
                "approved_antineoplastic_body_drug_claim_names": text_or_blank(
                    gene.approved_antineoplastic_body_drug_claim_names
                ),
                "clp_all_variant_count": int(gene.clp_all_variant_count),
                "clp_confirmed_somatic_variant_count": int(gene.clp_confirmed_somatic_variant_count),
                "claim_ceiling": "conditional_cross_source_possible_state_not_tree_truth",
            }
        )
    pair_rows.sort(key=lambda row: (str(row["chrom"]), str(row["region_a"])))
    pair_columns = list(pair_rows[0].keys())
    pair_path = output_dir / "hcc1395_pair_region_possible_clone_comparison.tsv.gz"
    write_tsv(pair_path, pair_rows, pair_columns, gzip_output=True)

    checks: list[dict[str, Any]] = []

    def add_check(name: str, observed: Any, expected: Any, passed: bool) -> None:
        checks.append({"check": name, "observed": observed, "expected": expected, "pass": bool(passed)})

    add_check("sample_membership", sorted(coarse["sample"].unique()), sorted(EXPECTED_SAMPLES), set(coarse["sample"]) == set(EXPECTED_SAMPLES))
    add_check("annotation_row_conservation", len(annotation_rows), len(coarse), len(annotation_rows) == len(coarse) == 47377)
    add_check("annotation_key_uniqueness", len({(row["sample"], row["region"]) for row in annotation_rows}), len(annotation_rows), len({(row["sample"], row["region"]) for row in annotation_rows}) == len(annotation_rows))
    add_check("selected_shape_join_conservation", sum(row["final_shape_category"] not in {"unresolved", "incomplete"} for row in annotation_rows), len(selected), sum(row["final_shape_category"] not in {"unresolved", "incomplete"} for row in annotation_rows) == len(selected) == 37039)
    source_counts = Counter(str(row["tree_selection_source"]) for row in annotation_rows)
    add_check(
        "tree_selection_source_distribution",
        dict(sorted(source_counts.items())),
        {
            "not_evaluable_incomplete": 7492,
            "structural_topo1": 21976,
            "vaf_resolved_topogt1": 15063,
            "vaf_unresolved_topogt1": 2846,
        },
        source_counts
        == Counter(
            {
                "not_evaluable_incomplete": 7492,
                "structural_topo1": 21976,
                "vaf_resolved_topogt1": 15063,
                "vaf_unresolved_topogt1": 2846,
            }
        ),
    )
    add_check(
        "vaf_inferred_flag_conservation",
        sum(bool(row["vaf_inferred_final_shape"]) for row in annotation_rows),
        15063,
        sum(bool(row["vaf_inferred_final_shape"]) for row in annotation_rows) == 15063,
    )
    add_check("pair_fixed_region_count", len(pair_rows), 5720, len(pair_rows) == 5720)
    add_check("pair_bridge_pattern_coverage", len(bridge_index), 5720, len(bridge_index) == 5720)
    add_check("pair_gene_drug_flag_coverage", len(gene_drug_index), 5720, len(gene_drug_index) == 5720)
    add_check("pair_common_join_matches_bridge", pair_common_join_mismatches, 0, pair_common_join_mismatches == 0)
    informative_rows = [row for row in pair_rows if bool(row["external_partition_informative"])]
    add_check("pair_informative_multicluster_count", len(informative_rows), 34, len(informative_rows) == 34)
    add_check(
        "pair_informative_partition_exact_count",
        sum(row["external_partition_exact"] is True for row in informative_rows),
        21,
        sum(row["external_partition_exact"] is True for row in informative_rows) == 21,
    )
    add_check(
        "pair_vacuous_subclonal_jaccard_is_undefined",
        sum(
            bool(row["external_partition_evaluable"])
            and int(row["external_subclonal_union_n"]) == 0
            and row["external_subclonal_jaccard"] == ""
            for row in pair_rows
        ),
        4996,
        sum(
            bool(row["external_partition_evaluable"])
            and int(row["external_subclonal_union_n"]) == 0
            and row["external_subclonal_jaccard"] == ""
            for row in pair_rows
        )
        == 4996,
    )
    add_check("blocked_samples_fail_closed", sum(row["regional_possible_state"] != "unavailable_cn_gate" for row in annotation_rows if row["sample"] in CN_BLOCKED_SAMPLES), 0, all(row["regional_possible_state"] == "unavailable_cn_gate" for row in annotation_rows if row["sample"] in CN_BLOCKED_SAMPLES))
    add_check("cn_ready_position_join_unambiguous", sum(value["position_duplicates"] for value in py_qa.values()), 0, all(value["position_duplicates"] == 0 for value in py_qa.values()))
    add_check("construction_missing_groups", sum(value["missing_groups"] for value in construction_qa.values()), 0, all(value["missing_groups"] == 0 for value in construction_qa.values()))
    add_check("construction_missing_primary_units", sum(value["missing_primary_units"] for value in construction_qa.values()), 0, all(value["missing_primary_units"] == 0 for value in construction_qa.values()))
    add_check("historical_region_sites_have_exact_raw_alleles", sum(value["missing_raw_alleles"] for value in construction_qa.values()), 0, all(value["missing_raw_alleles"] == 0 for value in construction_qa.values()))
    for sample in EXPECTED_SAMPLES:
        expected_bins = {str(key): int(value) for key, value in ct_by_sample[sample]["C"]["primary_region_sum"].items()}
        observed_bins = count_bins(c_observed[sample])
        add_check(f"{sample}_C_read_group_distribution", observed_bins, expected_bins, observed_bins == expected_bins)
        sample_rows = [row for row in annotation_rows if row["sample"] == sample]
        add_check(f"{sample}_status_conservation", len(sample_rows), int((coarse["sample"] == sample).sum()), len(sample_rows) == int((coarse["sample"] == sample).sum()))

    checks_path = output_dir / "checks.tsv"
    write_tsv(checks_path, checks, ["check", "observed", "expected", "pass"])
    checks_passed = sum(bool(row["pass"]) for row in checks)
    if checks_passed != len(checks):
        raise RuntimeError(f"Validation failed: {checks_passed}/{len(checks)}")

    sample_summary: dict[str, Any] = {}
    for sample in EXPECTED_SAMPLES:
        rows = [row for row in annotation_rows if row["sample"] == sample]
        status_counts = Counter(str(row["regional_possible_state"]) for row in rows)
        high_counts = Counter(str(row["high_assignment_possible_state"]) for row in rows)
        site_count = sum(int(row["region_site_count"]) for row in rows)
        joined_count = sum(int(row["external_joined_mutations"]) for row in rows)
        sample_summary[sample] = {
            "regions": len(rows),
            "sites": site_count,
            "joined_mutations": joined_count,
            "join_fraction": joined_count / site_count if site_count else None,
            "regions_any_join": sum(int(row["external_joined_mutations"]) > 0 for row in rows),
            "regions_full_join": sum(int(row["external_joined_mutations"]) == int(row["region_site_count"]) for row in rows),
            "possible_state_counts": dict(sorted(status_counts.items())),
            "high_assignment_state_counts": dict(sorted(high_counts.items())),
        }
    pair_state_counts = Counter(str(row["cross_source_possible_state"]) for row in pair_rows)
    pair_high_state_counts = Counter(str(row["high_assignment_cross_source_state"]) for row in pair_rows)
    summary = {
        "schema_name": "intersubmod.region_possible_clone_annotations",
        "schema_version": "1.0.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS_WITH_CLAIM_CEILING",
        "task_type": "B_comprehensive_validation",
        "goals": ["G4", "G5"],
        "scope": "chr1-22 historical layered-v2 primary regions; five CN-ready PyClone main fits; two CN-blocked samples fail closed",
        "definitions": {
            "C_read_groups": "ALT-containing MINREAD-supported full R/A groups summed across primary HP1/HP2; not a clone count",
            "T_exact_candidate_forest_count": "Product of exact candidate-tree counts across primary HP units; not a topology count",
            "Topo_shape_count": "Product of rooted-unlabeled exact shape counts across primary HP units",
            "tree_selection_source": "structural_topo1 means the rooted-unlabeled topology shape is unique, although exact labeled tree structures may remain multiple; vaf_resolved_topogt1 is a VAF-ranked most-likely shape inference; unresolved/incomplete remain fail-closed",
            "external_cluster_count": "PyClone-VI fit-local cluster labels joined by exact CHROM:POS:REF>ALT mutation_id after enforcing unique CHROM+POS within each fit",
            "external_partition_informative": "Partition exactness is defined only when both sources represent at least two fit-local clusters on the common-joined mutation endpoint",
            "external_subclonal_jaccard": "Defined only when the cross-source subclonal mutation union is non-empty; otherwise NA with both_absent reason",
            "clonal": f"PyClone cellular_prevalence >= {CLONAL_CP_THRESHOLD:.2f}",
            "subclonal": f"PyClone cellular_prevalence < {CLONAL_CP_THRESHOLD:.2f}",
            "high_assignment": f"cluster_assignment_prob >= {HIGH_ASSIGNMENT_THRESHOLD:.2f}; sensitivity only",
        },
        "claim_ceiling": "Possible regional clonal/subclonal state conditional on selected PyClone-VI model, shared/source-specific CN assumptions, joined mutations, and historical region definition. Not observed clone identity, clone count, ancestry, or tree truth.",
        "samples": sample_summary,
        "hcc1395_pair": {
            "fixed_regions": len(pair_rows),
            "cross_source_possible_state_counts": dict(sorted(pair_state_counts.items())),
            "high_assignment_cross_source_state_counts": dict(sorted(pair_high_state_counts.items())),
            "both_subclone_signal": sum(bool(row["possible_subclone_signal_a"]) and bool(row["possible_subclone_signal_b"]) for row in pair_rows),
            "either_subclone_signal": sum(bool(row["possible_subclone_signal_a"]) or bool(row["possible_subclone_signal_b"]) for row in pair_rows),
        },
        "validation": {"passed": checks_passed, "total": len(checks), "failed": len(checks) - checks_passed},
        "outputs": {},
    }

    provenance_inputs = {
        name: {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        for name, path in input_paths.items()
    }
    provenance_inputs["historical_run_validation_scope"] = {
        "path": str(run_root / "VALIDATION_SCOPE.md"),
        "sha256": sha256_file(run_root / "VALIDATION_SCOPE.md"),
        "size_bytes": (run_root / "VALIDATION_SCOPE.md").stat().st_size,
    }
    for sample, path in py_paths.items():
        provenance_inputs[f"pyclone_{sample}"] = {
            "path": str(path.resolve()),
            "sha256": sha256_file(path.resolve()),
            "size_bytes": path.resolve().stat().st_size,
        }

    output_paths = {
        "annotations": annotation_path,
        "sample_summary": summary_path,
        "pair_comparison": pair_path,
        "checks": checks_path,
    }
    for key, path in output_paths.items():
        summary["outputs"][key] = {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
    summary_path_json = output_dir / "summary.json"
    summary_path_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_paths["summary"] = summary_path_json

    provenance = {
        "schema_name": "intersubmod.region_possible_clone_annotations.provenance",
        "schema_version": "1.0.0",
        "command": [str(Path(__file__).resolve()), *(__import__("sys").argv[1:])],
        "script": {"path": str(Path(__file__).resolve()), "sha256": sha256_file(Path(__file__).resolve())},
        "inputs": provenance_inputs,
        "pyclone_qa": py_qa,
        "construction_qa": construction_qa,
        "outputs": {
            key: {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for key, path in output_paths.items()
        },
    }
    provenance_path = output_dir / "provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    execution = output_dir / "execution_record.md"
    execution.write_text(
        "\n".join(
            [
                "<!--",
                f"建立時間: {datetime.now().astimezone().isoformat()}",
                "目標: 產生 region-level possible clone/subclone conditional annotations",
                "處理範圍: chr1-22；historical layered-v2；5 CN-ready samples + 2 fail-closed samples",
                "關聯檔案: InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/scripts/build_region_possible_clone_annotations.py",
                "-->",
                "",
                "# Region possible clone/subclone annotation execution record",
                "",
                "## 結論",
                "",
                f"PASS WITH CLAIM CEILING — {checks_passed}/{len(checks)} checks PASS。`C_read_groups` 與 `external_cluster_count` 分離；possible state 不是 clone truth。",
                "",
                "## 輸入路徑",
                "",
                *[f"- `{path}`" for path in provenance_inputs.values() for path in [path["path"]]],
                "",
                "## 執行命令",
                "",
                "```bash",
                f"python3 {Path(__file__).resolve()} --output-dir {output_dir}",
                "```",
                "",
                "## 輸出路徑",
                "",
                *[f"- `{path}`" for path in [*output_paths.values(), provenance_path]],
                "",
                "## 實際輸出片段",
                "",
                "```text",
                f"regions={len(annotation_rows)}",
                f"hcc_fixed_pair_regions={len(pair_rows)}",
                f"checks={checks_passed}/{len(checks)} PASS",
                f"pair_state_counts={dict(sorted(pair_state_counts.items()))}",
                "```",
                "",
                "## Claim ceiling",
                "",
                "PyClone-VI CP／cluster只提供條件式 possible-state annotation；不能確認 clone 數、subclone identity、祖先方向或真實演化樹。High-assignment subset 可能因 subclonal union 退化而產生 selection-induced perfect agreement。",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "status": "PASS",
                "regions": len(annotation_rows),
                "pair_regions": len(pair_rows),
                "checks_passed": checks_passed,
                "checks_total": len(checks),
                "output_dir": str(output_dir),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
