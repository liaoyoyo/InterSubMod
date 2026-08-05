#!/usr/bin/env python3
"""Tests for the fail-closed professor-facing HTML report builder."""

from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "scripts" / "build_validated_html_report.py"
FIXTURE = HERE / "fixtures" / "validated_html_report_bundle.json"
SPEC = importlib.util.spec_from_file_location("validated_report", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = REPORT
SPEC.loader.exec_module(REPORT)


class HtmlProbe(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.tags: dict[str, int] = {}
        self.lang = None
        self.title_depth = 0
        self.title_text: list[str] = []
        self.ids: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        self.tags[tag] = self.tags.get(tag, 0) + 1
        values = dict(attrs)
        if values.get("id"):
            self.ids.append(values["id"])
        if tag == "html":
            self.lang = values.get("lang")
        if tag == "title":
            self.title_depth += 1

    def handle_endtag(self, tag: str) -> None:
        if tag == "title":
            self.title_depth -= 1

    def handle_data(self, data: str) -> None:
        if self.title_depth:
            self.title_text.append(data)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def split_count_tree(value, parts: int):
    """Split production-shaped cells while preserving invariant metadata."""
    if isinstance(value, dict):
        children = [{} for _ in range(parts)]
        for key, child_value in value.items():
            split_values = split_count_tree(child_value, parts)
            for index, split_value in enumerate(split_values):
                children[index][key] = split_value
        return children
    if isinstance(value, bool):
        return [value for _ in range(parts)]
    if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
        quotient, remainder = divmod(value, parts)
        return [quotient + (index < remainder) for index in range(parts)]
    if isinstance(value, str):
        return [value for _ in range(parts)]
    raise TypeError(f"fixture count tree contains unsupported value: {value!r}")


def ratio_payload(numerator: int | float, denominator: int | float, label: str):
    if denominator == 0:
        return {
            "numerator": numerator,
            "denominator": denominator,
            "value": None,
            "percent": None,
            "denominator_label": label,
            "reason": "denominator_is_zero",
        }
    return {
        "numerator": numerator,
        "denominator": denominator,
        "value": numerator / denominator,
        "percent": 100.0 * numerator / denominator,
        "denominator_label": label,
        "reason": None,
    }


def numeric_component_cell(distribution: dict[str, int]):
    normalized = {str(key): int(value) for key, value in distribution.items() if value}
    n_components = sum(normalized.values())
    k1 = normalized.get("1", 0)
    k_gt1 = sum(value for key, value in normalized.items() if int(key) > 1)
    return {
        "n_components": n_components,
        "k_equals_1": k1,
        "k_greater_than_1": k_gt1,
        "k_distribution": normalized,
        "max_k_component_sites": max((int(key) for key in normalized), default=0),
        "active_site_membership_mass": sum(
            int(key) * value for key, value in normalized.items()
        ),
        "k_equals_1_share": ratio_payload(
            k1, n_components, "components_in_same_dataset_basis_threshold"
        ),
        "k_greater_than_1_share": ratio_payload(
            k_gt1, n_components, "components_in_same_dataset_basis_threshold"
        ),
    }


def numeric_candidate_cell(
    raw: dict, *, h_zero: int | None = None,
    tied_consistent: int | None = None,
    tree_partition_counts: dict[str, int] | None = None,
):
    complete = raw["solver_complete_units"]
    unique = raw["quality_primary_unique_vertex_units"]
    tied = raw["quality_primary_tied_vertex_units"]
    incomplete = raw["solver_incomplete_or_not_run_units"]
    optimizer_abstain = raw["rank_abstain_units"] - incomplete
    if h_zero is None:
        h_zero = complete // 2
    h_one = complete - h_zero
    if tied_consistent is None:
        tied_consistent = tied // 2
    tied_inconsistent = tied - tied_consistent
    evaluated = raw["topology_evaluated_units"]
    coarse_unique = raw["coarse_topology_class_unique_units"]
    coarse_multiple = raw["coarse_topology_multiple_class_units"]
    parent_unique = raw["parent_edge_assignment_unique_units"]
    if tree_partition_counts is None:
        # Used only by legacy helper callers.  Materialized FINAL fixtures pass
        # the exact per-unit partition reconstructed while generating rows.
        tree_partition_counts = {
            "T_EQ_1_V_EQ_1": 0,
            "T_GT_1_V_EQ_1": 0,
            "T_GT_1_V_GT_1": complete,
        }
    return {
        "n_solver_complete_candidate_units": complete,
        "n_candidate_vertex_sets_V": raw["distinct_vertex_sets_V_complete_units"],
        "n_parent_edge_trees_T": raw["raw_tree_candidates_T_complete_units"],
        "tree_vertex_partition": {
            "counts": copy.deepcopy(tree_partition_counts),
            "denominator": complete,
            "shares": {
                key: ratio_payload(value, complete, "solver_complete_candidate_units")
                for key, value in tree_partition_counts.items()
            },
            "definition": REPORT.TREE_VERTEX_PARTITION_DEFINITION,
            "bucket_definitions": copy.deepcopy(REPORT.TREE_VERTEX_BUCKET_DEFINITIONS),
        },
        "unique_first": unique,
        "tied_first": tied,
        "solver_complete_optimizer_abstain": optimizer_abstain,
        "unique_first_share": ratio_payload(unique, complete, "solver_complete_candidate_units"),
        "tied_first_share": ratio_payload(tied, complete, "solver_complete_candidate_units"),
        "optimizer_abstain_share": ratio_payload(
            optimizer_abstain, complete, "solver_complete_candidate_units"
        ),
        "tied_by_coarse_topology": {
            "consistent": tied_consistent,
            "inconsistent": tied_inconsistent,
            "denominator": tied,
            "consistent_share": ratio_payload(
                tied_consistent, tied, "tied_first_solver_complete_units"
            ),
            "inconsistent_share": ratio_payload(
                tied_inconsistent, tied, "tied_first_solver_complete_units"
            ),
            "definition": (
                "consistent iff the union of coarse topology classes across all tied winning "
                "vertex sets has cardinality one; this does not claim exact parent-edge uniqueness"
            ),
        },
        "topology": {
            "evaluated_units": evaluated,
            "coarse_class_unique_units": coarse_unique,
            "coarse_class_multiple_units": coarse_multiple,
            "coarse_unique_class_counts": copy.deepcopy(
                raw["coarse_topology_unique_class_counts"]
            ),
            "coarse_ambiguous_class_set_counts": copy.deepcopy(
                raw["coarse_topology_ambiguous_class_set_counts"]
            ),
            "parent_edge_assignment_unique_units": parent_unique,
            "exact_topology_proven_unique_units": raw[
                "exact_topology_proven_unique_units"
            ],
        },
        "h_star_distribution": {"0": h_zero, "1": h_one},
        "coarse_topology_class_inclusion_counts": copy.deepcopy(
            raw.get("topology_class_inclusion_counts", {})
        ),
    }


def numeric_rank_cell(raw: dict, *, candidate: dict | None = None):
    units = raw["n_component_hp_units"]
    complete = raw["solver_complete_units"]
    projections = raw["molecule_component_projections"]
    informative = raw["informative_scoring_molecules"]
    structural = raw["structural_retained_molecules"]
    below = raw["below_minread_scoring_molecules"]
    candidate = candidate or numeric_candidate_cell(raw)
    unique = raw["quality_primary_unique_vertex_units"]
    tied = raw["quality_primary_tied_vertex_units"]
    abstain = raw["rank_abstain_units"]
    component_mass = raw["k_component_sites_total"]
    active_mass = raw["k_observed_alt_active_total"]
    nonstructural_mass = raw["not_structural_alt_active_sites_total"]
    evaluated = raw["topology_evaluated_units"]
    coarse_unique = raw["coarse_topology_class_unique_units"]
    coarse_multiple = raw["coarse_topology_multiple_class_units"]
    return {
        "units": {
            "n_component_hp_unit_evaluations": units,
            "solver_complete": complete,
            "solver_incomplete_or_not_run": raw["solver_incomplete_or_not_run_units"],
            "solver_complete_share": ratio_payload(
                complete, units, "component_hp_unit_evaluations"
            ),
        },
        "molecule_funnel": {
            "component_projections": projections,
            "informative_scoring_molecules": informative,
            "all_X_excluded_molecules": raw["all_x_excluded_molecules"],
            "structural_retained_molecules": structural,
            "below_structural_minread_but_scored_molecules": below,
            "informative_share_of_projections": ratio_payload(
                informative, projections, "molecule_component_projections"
            ),
            "structural_share_of_informative": ratio_payload(
                structural, informative, "informative_scoring_molecules"
            ),
            "below_minread_share_of_informative": ratio_payload(
                below, informative, "informative_scoring_molecules"
            ),
        },
        "partial_read_funnel": {
            "structural_partial_pattern_groups": raw["structural_partial_pattern_groups"],
            "coverage_denominator": raw["partial_group_coverage_denominator"],
            "covered": raw["partial_groups_covered"],
            "unsatisfied": raw["partial_groups_unsatisfied"],
            "covered_share": ratio_payload(
                raw["partial_groups_covered"],
                raw["partial_group_coverage_denominator"],
                "structural_partial_group_constraints",
            ),
            "full_detail": copy.deepcopy(raw["partial_pattern_funnel"]),
        },
        "candidate_structure": {
            "raw_parent_edge_trees_T_complete_units": raw[
                "raw_tree_candidates_T_complete_units"
            ],
            "distinct_optimal_vertex_sets_V_complete_units": raw[
                "distinct_vertex_sets_V_complete_units"
            ],
            "mean_T_per_solver_complete_unit": ratio_payload(
                raw["raw_tree_candidates_T_complete_units"],
                complete,
                "solver_complete_units",
            ),
            "mean_V_per_solver_complete_unit": ratio_payload(
                raw["distinct_vertex_sets_V_complete_units"],
                complete,
                "solver_complete_units",
            ),
            "candidate_table": candidate,
        },
        "ranking_outcome": {
            "unique_first": unique,
            "tied_first": tied,
            "abstain_all_causes": abstain,
            "unique_share": ratio_payload(unique, units, "component_hp_unit_evaluations"),
            "tied_share": ratio_payload(tied, units, "component_hp_unit_evaluations"),
            "abstain_share": ratio_payload(abstain, units, "component_hp_unit_evaluations"),
            "selection_status_counts": copy.deepcopy(raw["selection_status_counts"]),
        },
        "topology": {
            "evaluated_units": evaluated,
            "coarse_class_unique_units": coarse_unique,
            "coarse_class_multiple_units": coarse_multiple,
            "coarse_unique_class_counts": copy.deepcopy(
                raw["coarse_topology_unique_class_counts"]
            ),
            "coarse_ambiguous_class_set_counts": copy.deepcopy(
                raw["coarse_topology_ambiguous_class_set_counts"]
            ),
            "parent_edge_assignment_unique_units": raw[
                "parent_edge_assignment_unique_units"
            ],
            "exact_topology_proven_unique_units": raw[
                "exact_topology_proven_unique_units"
            ],
            "coarse_class_unique_share": ratio_payload(
                coarse_unique, evaluated, "topology_evaluated_units"
            ),
            "coarse_class_multiple_share": ratio_payload(
                coarse_multiple, evaluated, "topology_evaluated_units"
            ),
            "exact_topology_proven_unique_share": ratio_payload(
                raw["exact_topology_proven_unique_units"],
                evaluated,
                "topology_evaluated_units",
            ),
        },
        "effective_k": {
            "component_site_mass": component_mass,
            "observed_ALT_active_mass": active_mass,
            "not_structural_ALT_active_mass": nonstructural_mass,
            "k_route_counts": copy.deepcopy(raw["k_route_counts"]),
            "observed_ALT_active_share_of_component_site_mass": ratio_payload(
                active_mass, component_mass, "component_site_mass"
            ),
            "not_structural_ALT_active_share_of_component_site_mass": ratio_payload(
                nonstructural_mass, component_mass, "component_site_mass"
            ),
            "route_shares_of_unit_evaluations": {
                key: ratio_payload(value, units, "component_hp_unit_evaluations")
                for key, value in raw["k_route_counts"].items()
            },
        },
    }


def _expand_count_mapping(mapping: dict[str, int]) -> list[tuple[str, ...]]:
    expanded: list[tuple[str, ...]] = []
    for key, count in mapping.items():
        values = tuple(key.split("|"))
        expanded.extend([values] * count)
    return expanded


def generate_candidate_rows_for_stratum(
    dataset: str, basis: str, threshold: str, raw: dict,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Generate a complete candidate stream matching one production-shaped cell."""

    complete = raw["solver_complete_units"]
    n_v = raw["distinct_vertex_sets_V_complete_units"]
    n_t = raw["raw_tree_candidates_T_complete_units"]
    unique = raw["quality_primary_unique_vertex_units"]
    tied = raw["quality_primary_tied_vertex_units"]
    optimizer_abstain = complete - unique - tied
    assert complete >= 0 and n_t >= n_v >= complete
    assert optimizer_abstain >= 0

    outcomes = ["unique"] * unique + ["tied"] * tied + ["abstain"] * optimizer_abstain
    base_v = [2 if outcome == "tied" else 1 for outcome in outcomes]
    assert sum(base_v) <= n_v
    for index in range(n_v - sum(base_v)):
        base_v[index % max(1, complete)] += 1

    topology_specs = _expand_count_mapping(raw["coarse_topology_unique_class_counts"])
    topology_specs.extend(
        _expand_count_mapping(raw["coarse_topology_ambiguous_class_set_counts"])
    )
    assert len(topology_specs) == raw["topology_evaluated_units"] == unique + tied
    parent_unique_remaining = raw["parent_edge_assignment_unique_units"]
    assert parent_unique_remaining <= unique

    parent_counts_by_unit = [[1] * count for count in base_v]
    # A unique winner with one parent assignment is the only way a unit can be
    # parent-edge/exact-topology proven unique.  Force every other unique unit
    # to have >1 winning parent choices before allocating arbitrary extra T.
    required_extra = unique - parent_unique_remaining
    assert n_t - n_v >= required_extra
    for unit_index in range(parent_unique_remaining, unique):
        parent_counts_by_unit[unit_index][0] += 1
    remaining_extra = n_t - n_v - required_extra
    for index in range(remaining_extra):
        unit_index = index % max(1, complete)
        row_index = -1
        parent_counts_by_unit[unit_index][row_index] += 1

    rows: list[dict[str, str]] = []
    partition_counts = {key: 0 for key in REPORT.TREE_VERTEX_BUCKETS}
    h_distribution: dict[str, int] = {"0": 0, "1": 0}
    tied_consistent = 0
    topology_index = 0
    hp_family = basis.removeprefix("PS_HP")
    for unit_index, (outcome, row_count, parent_counts) in enumerate(
        zip(outcomes, base_v, parent_counts_by_unit)
    ):
        unit_t = sum(parent_counts)
        if unit_t == 1 and row_count == 1:
            partition_counts["T_EQ_1_V_EQ_1"] += 1
        elif unit_t > 1 and row_count == 1:
            partition_counts["T_GT_1_V_EQ_1"] += 1
        else:
            partition_counts["T_GT_1_V_GT_1"] += 1
        h_star = unit_index % 2
        h_distribution[str(h_star)] += 1
        classes = () if outcome == "abstain" else topology_specs[topology_index]
        if outcome != "abstain":
            topology_index += 1
        if outcome == "tied" and len(classes) == 1:
            tied_consistent += 1
        unit_key = (
            f"{dataset}|chr1|c{threshold}_{basis}_{unit_index:05d}|{threshold}|"
            f"HP{hp_family}|PS{1000 + unit_index}"
        )
        for row_index in range(row_count):
            if outcome == "unique":
                status = "UNIQUE_WINNER" if row_index == 0 else "NON_WINNER"
            elif outcome == "tied":
                status = "TIED_WINNER" if row_index < 2 else "NON_WINNER"
            else:
                status = "ABSTAIN_UNIT_OPTIMIZER"
            winner_classes = classes
            if outcome == "tied" and len(classes) > 1:
                winner_classes = (classes[row_index % len(classes)],) if row_index < 2 else classes
            elif status == "NON_WINNER":
                winner_classes = ("single-only",)
            topology_json = json.dumps(list(winner_classes), separators=(",", ":"))
            full_mask = 1 + (row_index % 6)
            states = {"0": "RRRR", str(full_mask): "ARRR"}
            roles = {"0": ["root"], str(full_mask): ["full-observed"]}
            if h_star:
                connector = str(8 + row_index)
                states[connector] = "AAAR"
                roles[connector] = ["connector"]
            rows.append(
                {
                    "unit_key": unit_key,
                    "dataset": dataset,
                    "chrom": "chr1",
                    "component_id": f"c{threshold}_{basis}_{unit_index:05d}",
                    "threshold": threshold,
                    "hp_family": hp_family,
                    "ps": str(1000 + unit_index),
                    "candidate_id": f"c{row_index:05d}",
                    "vertex_set_id": f"{unit_key}|v{row_index:05d}",
                    "vertex_states": json.dumps(states, separators=(",", ":")),
                    "vertex_roles": json.dumps(roles, separators=(",", ":")),
                    "parent_choice_count": str(parent_counts[row_index]),
                    "profile_log_likelihood": str(-10.0 - row_index),
                    "relative_log_likelihood": "0.0" if status in {"UNIQUE_WINNER", "TIED_WINNER"} else "-1.0",
                    "mixture_weights_pi": '{"0":1.0}',
                    "winner_status": status,
                    "tie_group": "TOP_TIE" if status == "TIED_WINNER" else "",
                    "coarse_topology_class": topology_json,
                    "candidate_set_complete": "true",
                }
            )
    assert topology_index == len(topology_specs)
    assert sum(partition_counts.values()) == complete
    return rows, {
        "tree_partition_counts": partition_counts,
        "h_distribution": h_distribution,
        "tied_consistent": tied_consistent,
    }


class ReportBuilderTest(unittest.TestCase):
    def setUp(self) -> None:
        self.bundle = json.loads(FIXTURE.read_text(encoding="utf-8"))
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self) -> None:
        # The S13 fixture deliberately models an immutable 0555/0444 release
        # contract.  Restore temporary directories so tempfile can remove its
        # own sandbox after each test.
        for path in (self.root, *self.root.rglob("*")):
            if path.is_dir() and not path.is_symlink():
                path.chmod(0o700)
        self.temp.cleanup()

    def materialize(self, bundle=None):
        bundle = copy.deepcopy(bundle if bundle is not None else self.bundle)
        paths = {}

        canonical = bundle["canonical"]
        canonical["canonical"]["aggregate"]["topology_classes"].setdefault(
            "impossible_exact_unique_topology_multiple", 0
        )
        for row in canonical["canonical"]["samples"]:
            row.setdefault("no_primary_lineage", 1)
            row.setdefault("primary_units", 11)
            row.setdefault(
                "topology_classes",
                {
                    "exact_and_topology_unique": 2,
                    "topology_unique_exact_multiple": 2,
                    "topology_multiple_exact_multiple": 4,
                    "incomplete": 1,
                    "impossible_exact_unique_topology_multiple": 0,
                },
            )

        m2_pilot_counts = bundle["m2_pilot"]["counts"]
        m2_pilot_counts["canonical_eligible_alignments"] = sum(m2_pilot_counts["raw_HP_counts"].values())
        pilot_runner = self.root / "run_exact_symbolic_pilot.py"
        pilot_solver = self.root / "hypercube_exact.py"
        pilot_runner.write_text("# fixture pilot runner\n", encoding="utf-8")
        pilot_solver.write_text("# fixture research solver\n", encoding="utf-8")
        bundle["pilot"]["implementation"] = {
            "runner": str(pilot_runner.resolve()),
            "runner_sha256": sha256(pilot_runner),
            "research_solver": str(pilot_solver.resolve()),
            "research_solver_sha256": sha256(pilot_solver),
        }

        extraction = bundle["m2_extraction"]
        extraction["n_results"] = len(REPORT.DATASETS) * len(REPORT.AUTOSOMES)
        extraction_component_summary = extraction["aggregate"][
            "component_summary_by_linkage_basis"
        ]
        child_component_cells: dict[tuple[str, str, str, str], dict] = {}
        for basis in ("PS_HP1", "PS_HP2"):
            for threshold in ("1", "2", "3", "5"):
                global_distribution = extraction_component_summary[basis][threshold][
                    "k_component_sites_distribution"
                ]
                dataset_distributions = split_count_tree(
                    global_distribution, len(REPORT.DATASETS)
                )
                for dataset_index, dataset in enumerate(REPORT.DATASETS):
                    chrom_distributions = split_count_tree(
                        dataset_distributions[dataset_index], len(REPORT.AUTOSOMES)
                    )
                    for chrom_index, chrom in enumerate(REPORT.AUTOSOMES):
                        distribution = {
                            key: value
                            for key, value in chrom_distributions[chrom_index].items()
                            if value
                        }
                        n_components = sum(distribution.values())
                        child_component_cells[(dataset, chrom, basis, threshold)] = {
                            "n_components": n_components,
                            "n_singletons_k1": distribution.get("1", 0),
                            "n_multisite_k_gt1": sum(
                                value for key, value in distribution.items() if int(key) > 1
                            ),
                            "max_k_component_sites": max(
                                (int(key) for key, value in distribution.items() if value),
                                default=0,
                            ),
                            "n_active_site_memberships": sum(
                                int(key) * value for key, value in distribution.items()
                            ),
                            "k_component_sites_distribution": distribution,
                        }
                aggregate_cell = extraction_component_summary[basis][threshold]
                aggregate_cell["n_active_site_memberships"] = sum(
                    int(key) * value for key, value in global_distribution.items()
                )

        extraction_children: dict[tuple[str, str], tuple[Path, dict]] = {}
        extraction_results = []
        for dataset in REPORT.DATASETS:
            for chrom in REPORT.AUTOSOMES:
                child_path = self.root / "full_extraction" / "samples" / dataset / chrom / "receipt.json"
                child_path.parent.mkdir(parents=True, exist_ok=True)
                child = {
                    "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
                    "schema_version": "1.2.0",
                    "all_pass": True,
                    "scope": {"dataset": dataset, "chrom": chrom},
                    "component_summary_by_linkage_basis": {
                        basis: {
                            threshold: copy.deepcopy(
                                child_component_cells[(dataset, chrom, basis, threshold)]
                            )
                            for threshold in ("1", "2", "3", "5")
                        }
                        for basis in ("PS_HP1", "PS_HP2")
                    },
                    "checks": {"fixture_component_summary_conserved": True},
                }
                child_path.write_text(
                    json.dumps(child, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8",
                )
                child_sidecar = child_path.with_name(f"{child_path.name}.sha256")
                child_sidecar.write_text(
                    f"{sha256(child_path)}  {child_path.name}\n", encoding="ascii"
                )
                extraction_children[(dataset, chrom)] = (child_path, child)
                extraction_results.append(
                    {
                        "dataset": dataset,
                        "chrom": chrom,
                        "status": "PASS",
                        "receipt": child,
                    }
                )
        extraction["results"] = extraction_results
        extraction["task_status_counts"] = {"PASS": 154, "REUSED_PASS": 0}
        split_extraction_counts = {
            key: split_count_tree(value, len(REPORT.DATASETS))
            for key, value in extraction["aggregate"]["counts"].items()
        }
        extraction["aggregate"]["by_dataset"] = {
            dataset: {
                "task_status_counts": {"PASS": 22},
                "counts": {key: values[index] for key, values in split_extraction_counts.items()},
            }
            for index, dataset in enumerate(REPORT.DATASETS)
        }

        ranking = bundle["m2_ranking"]
        m2_ranker = self.root / "build_m2_patterns_and_rank.py"
        m2_ranker.write_text("# fixture M2 ranker\n", encoding="utf-8")
        ranking["run_contract"] = {
            "ranker": {
                "path": str(m2_ranker.resolve()),
                "sha256": sha256(m2_ranker),
            },
            "method_contract": copy.deepcopy(REPORT.EXPECTED_METHOD_CONTRACT),
        }
        combined_template = ranking["aggregate"]["by_linkage_basis_and_threshold"].pop("primary_hp")["3"]
        combined_template["coarse_topology_class_unique_units"] = combined_template[
            "coarse_topology_unique_units"
        ]
        combined_template["coarse_topology_multiple_class_units"] = combined_template[
            "coarse_topology_ambiguous_units"
        ]
        topology_key_map = {
            "single_only": "single-only",
            "sister_only": "sister-only",
            "direct_only": "direct-only",
            "sister_and_direct": "sister+direct",
        }
        combined_template["coarse_topology_unique_class_counts"] = {
            topology_key_map[key]: value
            for key, value in combined_template["coarse_topology_unique_class_counts"].items()
        }
        combined_template["topology_class_inclusion_counts"] = {
            topology_key_map[key]: value
            for key, value in combined_template["topology_class_inclusion_counts"].items()
        }
        ambiguous_key_map = {
            "single_or_direct": "single-only|direct-only",
            "sister_or_sister_and_direct": "sister-only|sister+direct",
        }
        combined_template["ambiguous_topology_class_set_counts"] = {
            ambiguous_key_map[key]: value
            for key, value in combined_template["ambiguous_topology_class_set_counts"].items()
        }
        combined_template["coarse_topology_ambiguous_class_set_counts"] = copy.deepcopy(
            combined_template["ambiguous_topology_class_set_counts"]
        )
        inclusion_counts = copy.deepcopy(
            combined_template["coarse_topology_unique_class_counts"]
        )
        for class_set, count in combined_template[
            "coarse_topology_ambiguous_class_set_counts"
        ].items():
            for topology_class in class_set.split("|"):
                inclusion_counts[topology_class] = (
                    inclusion_counts.get(topology_class, 0) + count
                )
        combined_template["topology_class_inclusion_counts"] = inclusion_counts
        basis_cells = split_count_tree(combined_template, 2)
        ranking["aggregate"]["by_linkage_basis_and_threshold"] = {
            basis: {threshold: copy.deepcopy(basis_cells[index]) for threshold in ("1", "2", "3", "5")}
            for index, basis in enumerate(("PS_HP1", "PS_HP2"))
        }
        for dataset_index, dataset in enumerate(REPORT.DATASETS):
            dataset_nested = {}
            for basis in ("PS_HP1", "PS_HP2"):
                dataset_nested[basis] = {}
                for threshold in ("1", "2", "3", "5"):
                    parts = split_count_tree(
                        ranking["aggregate"]["by_linkage_basis_and_threshold"][basis][threshold],
                        len(REPORT.DATASETS),
                    )
                    dataset_nested[basis][threshold] = parts[dataset_index]
            ranking["by_dataset"][dataset]["by_linkage_basis_and_threshold"] = dataset_nested

        generated_candidate_rows: list[dict[str, str]] = []
        candidate_cell_metadata: dict[tuple[str, str, str], dict[str, object]] = {}
        for dataset in REPORT.DATASETS:
            for basis in ("PS_HP1", "PS_HP2"):
                for threshold in ("1", "2", "3", "5"):
                    raw = ranking["by_dataset"][dataset]["by_linkage_basis_and_threshold"][basis][threshold]
                    rows, metadata = generate_candidate_rows_for_stratum(
                        dataset, basis, threshold, raw
                    )
                    generated_candidate_rows.extend(rows)
                    candidate_cell_metadata[(dataset, basis, threshold)] = metadata
        generated_candidate_rows.sort(
            key=lambda row: (row["unit_key"], row["candidate_id"])
        )
        bundle["candidate_rows"] = generated_candidate_rows

        for key in ("canonical", "pilot", "m2_extraction", "m2_pilot"):
            path = self.root / f"{key}.json"
            path.write_text(json.dumps(bundle[key], ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            paths[key] = path

        funnel_rows = []
        funnel_sources = []
        branch_totals = {
            "excluded_by_longphase_filter": 0,
            "out_of_scope_non_autosomal": 0,
            "max_snv_excluded": 0,
            "positional_singleton": 0,
            "retained": 0,
        }
        for row in canonical["canonical"]["samples"]:
            dataset = row["sample"]
            raw = row["tree_input_records"] + 10
            branches = {
                "excluded_by_longphase_filter": raw - row["tree_input_records"],
                "out_of_scope_non_autosomal": row["tree_input_records"] - row["autosomal_biallelic_sSNV"],
                "max_snv_excluded": 0,
                "positional_singleton": row["autosomal_biallelic_sSNV"] - row["retained_sSNV"],
                "retained": row["retained_sSNV"],
            }
            for key, value in branches.items():
                branch_totals[key] += value
            ledger_path = self.root / "ledgers" / dataset / f"ssnv_site_ledger_{dataset}.summary.json"
            ledger_path.parent.mkdir(parents=True, exist_ok=True)
            ledger_path.write_text(json.dumps({"dataset": dataset, "fixture": True}) + "\n", encoding="utf-8")
            funnel_sources.append({"dataset": dataset, "path": str(ledger_path), "sha256": sha256(ledger_path)})
            funnel_rows.append(
                {
                    "dataset": dataset,
                    "raw_clairs_records": raw,
                    "longphase_s_recalibrated_pass": row["tree_input_records"],
                    "autosomal_biallelic_sSNV": row["autosomal_biallelic_sSNV"],
                    "retained_sSNV": row["retained_sSNV"],
                    "branch_counts": branches,
                }
            )
        funnel_aggregate = {
            key: sum(row[key] for row in funnel_rows)
            for key in ("raw_clairs_records", "longphase_s_recalibrated_pass", "autosomal_biallelic_sSNV", "retained_sSNV")
        }
        raw = funnel_aggregate["raw_clairs_records"]
        tree = funnel_aggregate["longphase_s_recalibrated_pass"]
        autosomal = funnel_aggregate["autosomal_biallelic_sSNV"]
        retained = funnel_aggregate["retained_sSNV"]
        funnel_aggregate.update(
            {
                "branch_counts": branch_totals,
                "relative_ratios": {
                    "longphase_pass_over_raw": tree / raw,
                    "autosomal_over_longphase_pass": autosomal / tree,
                    "retained_over_autosomal": retained / autosomal,
                },
                "total_ratios_over_raw": {"autosomal": autosomal / raw, "retained": retained / raw},
            }
        )
        funnel = {
            "schema_name": "intersubmod_current_ssnv_funnel_receipt",
            "schema_version": "1.0.0",
            "task_type": "B_comprehensive_validation",
            "all_pass": True,
            "scope": {
                "datasets": list(REPORT.DATASETS),
                "dataset_count": 7,
                "biological_sample_count": 6,
            },
            "inputs": {
                "canonical_json": {"path": str(paths["canonical"]), "sha256": sha256(paths["canonical"])},
                "site_ledger_summaries": funnel_sources,
            },
            "aggregate": funnel_aggregate,
            "datasets": funnel_rows,
            "checks": {
                "all_7_dataset_ledgers_present": True,
                "all_ledger_checks_pass": True,
                "all_record_keys_unique": True,
                "all_branch_conservation_checks_pass": True,
                "all_dataset_counts_match_canonical": True,
                "aggregate_counts_match_canonical": True,
            },
        }
        paths["funnel"] = self.root / "funnel.json"
        paths["funnel"].write_text(json.dumps(funnel, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        funnel_verification_checks = {f"check_{index:03d}": True for index in range(1, 323)}
        funnel_verification = {
            "schema_name": "intersubmod_current_ssnv_funnel_independent_verification",
            "schema_version": "1.0.0",
            "task_type": "B_comprehensive_validation",
            "all_pass": True,
            "scope": {
                "datasets": list(REPORT.DATASETS),
                "dataset_count": 7,
                "biological_sample_count": 6,
                "chromosomes": "chr1-chr22 for autosomal_biallelic_sSNV and downstream",
            },
            "inputs": {
                "canonical_json": {
                    "path": str(paths["canonical"].resolve()),
                    "sha256": sha256(paths["canonical"]),
                },
                "receipt_under_test": {
                    "path": str(paths["funnel"].resolve()),
                    "sha256": sha256(paths["funnel"]),
                },
            },
            "recomputed": {"aggregate": copy.deepcopy(funnel_aggregate)},
            "checks": funnel_verification_checks,
            "n_checks": 322,
            "n_pass": 322,
            "n_fail": 0,
            "failures": [],
        }
        paths["funnel_verification"] = self.root / "funnel_verification.json"
        paths["funnel_verification"].write_text(
            json.dumps(funnel_verification, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

        m0_rows = self.root / "m0_rows.tsv.gz"
        with gzip.open(m0_rows, "wt", encoding="utf-8") as handle:
            handle.write("unit_key\n")
            for index in range(bundle["m0"]["aggregate"]["n_hp_lineage_units"]):
                handle.write(f"u{index:04d}\n")
        bundle["m0"]["outputs"] = {
            "rows": str(m0_rows),
            "rows_size_bytes": m0_rows.stat().st_size,
            "rows_sha256": sha256(m0_rows),
        }
        paths["m0"] = self.root / "m0.json"
        paths["m0"].write_text(json.dumps(bundle["m0"], ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        m0_verification = {
            "schema_name": "intersubmod.hypercube_m0_independent_verification",
            "schema_version": "1.1.0",
            "verdict": "PASS",
            "candidate_mode": "deep",
            "scope": {
                "selected_datasets": list(REPORT.DATASETS),
                "full_task_b_scope": True,
                "row_schema_matches": True,
            },
            "inputs": {
                "receipt_sha256": sha256(paths["m0"]),
                "receipt_size_bytes": paths["m0"].stat().st_size,
                "rows_sha256": sha256(m0_rows),
                "rows_size_bytes": m0_rows.stat().st_size,
            },
            "independently_recomputed_aggregate": {
                "n_hp_lineage_units": bundle["m0"]["aggregate"]["n_hp_lineage_units"],
                "selection_status_counts": bundle["m0"]["aggregate"]["selection_status_counts"],
                "by_dataset": bundle["m0"]["aggregate"]["by_dataset"],
            },
            "categorical_conservation": {"partition_conserves": True, "T_V_partition_conserves": True},
            "canonical_reconciliation": {
                "eligible_plus_excluded_equals_primary_hp_units": True,
                "missing_tsv_units": 0,
                "extra_tsv_units": 0,
                "n_candidate_units_deep_checked": bundle["m0"]["aggregate"]["n_hp_lineage_units"],
            },
            "checks": {
                "n_errors": 0,
                "receipt_tsv_hash_matches": True,
                "receipt_aggregate_matches_tsv": True,
                "strong_eligible_excluded_conservation": True,
            },
        }
        paths["m0_verification"] = self.root / "m0_verification.json"
        paths["m0_verification"].write_text(json.dumps(m0_verification, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

        ranking["upstream_extraction_receipt"]["sha256"] = sha256(paths["m2_extraction"])
        candidate_path = self.root / "m2_ps_aware_candidate_table.tsv.gz"
        with gzip.open(candidate_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, REPORT.CANDIDATE_TABLE_COLUMNS, delimiter="\t", extrasaction="raise")
            writer.writeheader()
            writer.writerows(bundle["candidate_rows"])
        ranking["candidate_table"] = {
            "schema_version": "2.0.0",
            "format": "tsv.gz",
            "path": candidate_path.name,
            "size_bytes": candidate_path.stat().st_size,
            "sha256": sha256(candidate_path),
            "semantic_sha256": hashlib.sha256(json.dumps(bundle["candidate_rows"], sort_keys=True).encode("utf-8")).hexdigest(),
            "n_rows": len(bundle["candidate_rows"]),
            "n_units": len({row["unit_key"] for row in bundle["candidate_rows"]}),
            "columns": list(REPORT.CANDIDATE_TABLE_COLUMNS),
            "sort_order": "unit_key,candidate_id",
        }
        ranking["task_index"] = [
            {
                "dataset": dataset,
                "chrom": chrom,
                "all_pass": True,
                "schema_version": "2.0.0",
                "parameters_match_extraction": True,
                "input_hashes_match_extraction": True,
                "upstream_outputs_verified": True,
                "no_cross_ps_pattern_pooling": True,
                "known_ps_never_mixed": True,
                "missing_ps_separate_diagnostic": True,
                "runtime_diagnostics_contract_valid": True,
                "method_contract_matches": True,
                "ranker_source_bound": True,
            }
            for dataset in REPORT.DATASETS
            for chrom in REPORT.AUTOSOMES
        ]
        paths["m2_ranking"] = self.root / "m2_ranking.json"
        paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

        release_root = self.root / "release_contract"
        if release_root.exists():
            # Some tests rematerialize into one TemporaryDirectory.  Thaw only
            # this disposable fixture before recreating the same contract.
            release_root.chmod(0o755)
            for existing in release_root.rglob("*"):
                if existing.is_dir() and not existing.is_symlink():
                    existing.chmod(0o755)
                elif existing.is_file() and not existing.is_symlink():
                    existing.chmod(0o644)
        snapshot_root = release_root / "code_snapshot"
        snapshot_root.mkdir(parents=True, exist_ok=True)
        release_roles = (
            "extractor", "lossless_contract", "full_extraction_runner",
            "ranker", "hypercube_solver", "full_ranking_runner",
            "single_task_pilot_verifier", "full_verifier",
            "input_identity_verifier", "release_contract_freezer",
            "canonical_manifest_schema",
        )
        release_entries = []
        release_bound_sources = {}
        for index, role in enumerate(release_roles):
            relative = Path("code_snapshot") / f"{index:02d}_{role}.py"
            physical = release_root / relative
            physical.write_text(f"# immutable fixture snapshot: {role}\n", encoding="utf-8")
            physical.chmod(0o444)
            observed = physical.lstat()
            digest = sha256(physical)
            release_entries.append(
                {
                    "role": role,
                    "snapshot": {
                        "path": relative.as_posix(),
                        "size_bytes": physical.stat().st_size,
                        "sha256": digest,
                        "mode_octal": "0444",
                        "st_nlink": observed.st_nlink,
                    },
                }
            )
            release_bound_sources[role] = {
                "path": str(physical.resolve()),
                "sha256": digest,
            }
        canonical_copy = release_root / "input_contract" / "canonical_manifest.json"
        canonical_copy.parent.mkdir(parents=True, exist_ok=True)
        canonical_copy.write_text('{"fixture":"canonical-v5"}\n', encoding="utf-8")
        canonical_copy.chmod(0o444)
        canonical_observed = canonical_copy.lstat()
        canonical_digest = sha256(canonical_copy)
        release_manifest_path = release_root / "m2_run_manifest.json"
        release_manifest = {
            "schema_name": "intersubmod.m2_release_run_manifest",
            "schema_version": "1.0.0",
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "authority_mode": "CANONICAL_V5_FROZEN",
            "validation_evidence_eligible": True,
            "scope": {
                "technical_datasets": 7,
                "biological_samples": 6,
                "chromosomes": 22,
                "chromosome_names": list(REPORT.AUTOSOMES),
                "tasks": 154,
                "datasets": list(REPORT.DATASETS),
            },
            "canonical_manifest": {
                "expected_sha256": canonical_digest,
                "immutable_copy": {
                    "path": "input_contract/canonical_manifest.json",
                    "size_bytes": canonical_copy.stat().st_size,
                    "sha256": canonical_digest,
                    "mode_octal": "0444",
                    "st_nlink": canonical_observed.st_nlink,
                },
            },
            "source_snapshot": {
                "n_files": 11,
                "entries": release_entries,
            },
            "checks": {
                "fixture_canonical_authority": True,
                "fixture_exact_eleven_snapshot": True,
            },
            "all_pass": True,
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": "m2_run_manifest.json.sha256",
                "covers": "m2_run_manifest.json",
            },
        }
        release_manifest_path.write_text(
            json.dumps(release_manifest, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        release_manifest_path.chmod(0o444)
        release_sidecar_path = release_manifest_path.with_name(
            f"{release_manifest_path.name}.sha256"
        )
        release_sidecar_path.write_text(
            f"{sha256(release_manifest_path)}  {release_manifest_path.name}\n",
            encoding="ascii",
        )
        release_sidecar_path.chmod(0o444)
        freezer_identity = release_bound_sources["release_contract_freezer"]
        release_binding = {
            "schema_name": "intersubmod.m2_release_binding",
            "schema_version": "1.0.0",
            "release_manifest": {
                "path": str(release_manifest_path.resolve()),
                "sha256": sha256(release_manifest_path),
                "semantic_sha256": REPORT._semantic_json_sha256(release_manifest),
                "sidecar": {
                    "path": str(release_sidecar_path.resolve()),
                    "sha256": sha256(release_sidecar_path),
                },
            },
            "authority_mode": "CANONICAL_V5_FROZEN",
            "validation_evidence_eligible": True,
            "canonical_input_manifest": {
                "path": str(canonical_copy.resolve()),
                "sha256": canonical_digest,
            },
            "snapshot_sources": release_bound_sources,
            "deep_release_verification": {
                "mode": "FROZEN_FREEZER_VERIFY_RELEASE_CONTRACT",
                "freezer_path": freezer_identity["path"],
                "freezer_sha256": freezer_identity["sha256"],
                "release_manifest_sha256": sha256(release_manifest_path),
                "verified_snapshot_files": 11,
                "all_pass": True,
            },
        }
        for directory in sorted(
            (path for path in release_root.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts), reverse=True,
        ):
            directory.chmod(0o555)
        release_root.chmod(0o555)

        extraction_session_path = self.root / "extraction_orchestration" / "session.json"
        extraction_session_path.parent.mkdir(parents=True, exist_ok=True)
        zero_conflict_snapshot = {
            "process_count": 0,
            "root_count": 0,
            "representatives": [],
        }
        extraction_session = {
            "schema_name": "intersubmod.m2_orchestration_session",
            "schema_version": "1.0.0",
            "stage": "extraction",
            "session_id": "6" * 64,
            "session_nonce": "fixture-session-nonce",
            "created_at_utc": "2026-07-16T00:00:01Z",
            "created_monotonic_ns": 20,
            "host_boot_id": "fixture-host-boot-id",
            "release_manifest": copy.deepcopy(release_binding["release_manifest"]),
            "release_binding_semantic_sha256": REPORT._semantic_json_sha256(release_binding),
            "run_contract_semantic_sha256": "e" * 64,
            "scope": {
                "datasets": list(REPORT.DATASETS),
                "chromosomes": list(REPORT.AUTOSOMES),
                "expected_tasks": 154,
            },
            "output_root": "/fixture/output",
            "producer_sources": [],
            "scheduler_policy": {"batch_size": 16},
            "parent_extraction": None,
            "resource_gate": {
                "ignore_resource_gate": False,
                **zero_conflict_snapshot,
                "observed_at_utc": "2026-07-16T00:00:00Z",
                "observed_monotonic_ns": 10,
                "host_boot_id": "fixture-host-boot-id",
                "process_snapshot_sha256": REPORT._semantic_json_sha256(
                    zero_conflict_snapshot
                ),
            },
            "receipt_integrity": {"mode": "fixture"},
        }
        extraction_session_path.write_text(
            json.dumps(extraction_session, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        paths["m2_extraction_session"] = extraction_session_path
        m2_verification = {
            "schema_name": "intersubmod.m2_full_independent_verification",
            "schema_version": "1.0.0",
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "release_binding": copy.deepcopy(release_binding),
            "post_input_identity": {
                "post_receipt": {
                    "path": str((self.root / "post.json").resolve()),
                    "sha256": "1" * 64,
                    "semantic_sha256": "2" * 64,
                },
                "immutable_pre_receipt": {
                    "path": str((self.root / "pre.json").resolve()),
                    "sha256": "3" * 64,
                    "semantic_sha256": "4" * 64,
                },
                "input_identity_snapshot_sha256": "5" * 64,
                "n_artifact_roles": 42,
                "exact_snapshot_equal": True,
                "claim_boundary": "Persistent PRE-to-POST equality; no transient-restore proof.",
            },
            "scope": {
                "datasets": list(REPORT.DATASETS),
                "chromosomes": list(REPORT.AUTOSOMES),
                "expected_tasks": 154,
                "n_technical_datasets": 7,
            },
            "verification_independence": {
                "imports_production_aggregator": False,
                "imports_production_ranker": False,
                "reads_bam": False,
                "recomputes_from_child_receipts_and_tables": True,
                "candidate_table_compared_rowwise_to_independent_child_reconstruction": True,
                "method_contract_verification_mode": (
                    "STATIC_EXACT_CONTRACT_AND_ACTUAL_SOURCE_SHA_BINDING"
                ),
                "profile_likelihood_recomputed_from_patterns_states_pi": True,
                "profile_likelihood_claim_boundary": (
                    "Fixture independently recomputes profile likelihood from R/A/X+BQ patterns, states, and persisted pi."
                ),
            },
            "extraction": {
                "receipt_sha256": sha256(paths["m2_extraction"]),
                "n_tasks": 154,
                "orchestration": {
                    "stage": "extraction",
                    "session_id": "6" * 64,
                    "session_receipt": {
                        "path": str(extraction_session_path.resolve()),
                        "sha256": sha256(extraction_session_path),
                    },
                    "terminal_receipt": {"path": "/fixture/extraction/full.json", "sha256": "8" * 64},
                    "legal_cumulative_counts": [8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 154],
                    "n_batches": 11,
                    "n_attested_children": 154,
                    "all_child_receipts_and_outputs_rehashed": True,
                    "all_grants_and_completions_session_bound": True,
                    "all_checkpoints_and_terminal_chain_bound": True,
                    "no_open_orphan_or_preseeded_child_accepted": True,
                    "claim_boundary": "Parent attestation under non-hostile same-UID assumption.",
                },
            },
            "ranking": {
                "receipt_sha256": sha256(paths["m2_ranking"]),
                "n_tasks": 154,
                "all_aggregate_cells_conserved": True,
                "orchestration": {
                    "stage": "ranking",
                    "session_id": "9" * 64,
                    "session_receipt": {"path": "/fixture/ranking/session.json", "sha256": "a" * 64},
                    "terminal_receipt": {"path": "/fixture/ranking/full.json", "sha256": "b" * 64},
                    "legal_cumulative_counts": [8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 154],
                    "n_batches": 11,
                    "n_attested_children": 154,
                    "all_child_receipts_and_outputs_rehashed": True,
                    "all_grants_and_completions_session_bound": True,
                    "all_checkpoints_and_terminal_chain_bound": True,
                    "no_open_orphan_or_preseeded_child_accepted": True,
                    "claim_boundary": "Parent attestation under non-hostile same-UID assumption.",
                },
                "candidate_table": {
                    "sha256": sha256(candidate_path),
                    "n_rows": ranking["candidate_table"]["n_rows"],
                    "n_units": ranking["candidate_table"]["n_units"],
                    "all_rows_match_independent_child_reconstruction": True,
                    "winner_partitions_conserved": True,
                },
                "method_contract_verification": {
                    "contract": copy.deepcopy(REPORT.EXPECTED_METHOD_CONTRACT),
                    "ranker_source_path": str(m2_ranker.resolve()),
                    "ranker_source_sha256": sha256(m2_ranker),
                    "n_children_exactly_matched_and_source_bound": 154,
                    "all_children_exactly_matched_and_source_bound": True,
                    "verification_mode": "STATIC_EXACT_CONTRACT_AND_ACTUAL_SOURCE_SHA_BINDING",
                },
                "profile_likelihood_independent_recomputation": {
                    "n_children": 154,
                    "n_units": ranking["candidate_table"]["n_units"],
                    "n_candidates": ranking["candidate_table"]["n_rows"],
                    "n_pattern_rows": 7,
                    "n_scoring_molecules": 19,
                    "max_abs_ll_delta": 1e-12,
                    "max_abs_relative_weight_delta": 2e-13,
                    "max_abs_gap_delta": 3e-14,
                    "max_abs_simplex_residual_delta": 4e-15,
                    "peak_pattern_rows_per_unit": 4,
                    "peak_candidates_per_unit": 3,
                    "peak_states_per_candidate": 5,
                    "peak_emission_cells_per_candidate": 20,
                    "all_profile_likelihoods_and_certificates_match": True,
                    "all_relative_weights_match": True,
                    "all_winner_tie_partitions_match": True,
                    "count_unit_definitions": {
                        "n_scoring_molecules": "molecule projection sum, not global unique molecules"
                    },
                    "numeric_contract": {"persisted_scalar_format": ".12g"},
                    "streaming_memory_contract": "one child and one primary unit at a time",
                    "child_summaries": [
                        {"dataset": dataset, "chrom": chrom}
                        for dataset in REPORT.DATASETS
                        for chrom in REPORT.AUTOSOMES
                    ],
                },
                "runtime_diagnostics": {
                    "n_child_runtime_files": ranking["runtime_diagnostics"]["n_child_runtime_files"],
                    "n_unit_evaluations": ranking["runtime_diagnostics"]["n_unit_evaluations"],
                    "metrics": copy.deepcopy(ranking["runtime_diagnostics"]["metrics"]),
                    "metrics_when_invoked": copy.deepcopy(
                        ranking["runtime_diagnostics"]["metrics_when_invoked"]
                    ),
                    "all_child_and_full_runtime_summaries_independently_recomputed": True,
                },
            },
            "checks": {
                "all_expected_extraction_tasks_verified": True,
                "all_expected_ranking_tasks_verified": True,
                "extraction_full_aggregate_exactly_recomputed": True,
                "ranking_core_aggregates_exactly_recomputed": True,
                "ranking_all_recomputed_cells_conserve": True,
                "candidate_table_physical_and_semantic_hashes_verified": True,
                "candidate_table_rows_match_child_sources": True,
                "candidate_table_winner_partitions_conserve": True,
                "all_child_method_contracts_exactly_compared_and_source_bound": True,
                "runtime_diagnostics_independently_recomputed": True,
                "profile_likelihood_recomputed_independently": True,
                "profile_likelihood_certificates_match": True,
                "profile_relative_weights_match": True,
                "profile_winner_tie_partitions_match": True,
                "release_contract_authenticated_and_eligible": True,
                "release_contract_all_snapshot_sources_rehashed": True,
                "extraction_and_ranking_bound_to_same_release": True,
                "full_runner_dependency_paths_and_shas_match_release": True,
                "frozen_scientific_and_scheduler_parameters_match": True,
                "post_input_identity_authenticated_and_exactly_equals_frozen_pre": True,
                "extraction_orchestration_session_batch_grant_completion_chain_verified": True,
                "ranking_orchestration_session_batch_grant_completion_chain_verified": True,
                "ranking_orchestration_bound_to_verified_extraction_session": True,
            },
            "all_pass": True,
        }
        paths["m2_verification"] = self.root / "m2_verification.json"
        paths["m2_verification"].write_text(json.dumps(m2_verification, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

        numeric_summary = bundle["m2_numeric_summary"]
        extraction_components = extraction["aggregate"][
            "component_summary_by_linkage_basis"
        ]
        numeric_extraction_by_dataset = {
            dataset: {"counts": {}, "component_by_linkage_basis_threshold": {}}
            for dataset in REPORT.DATASETS
        }
        overall_components = {}
        for basis in ("PS_HP1", "PS_HP2"):
            overall_components[basis] = {}
            for threshold in ("1", "2", "3", "5"):
                distribution = extraction_components[basis][threshold][
                    "k_component_sites_distribution"
                ]
                overall_components[basis][threshold] = numeric_component_cell(
                    distribution
                )
                distribution_parts = split_count_tree(
                    distribution, len(REPORT.DATASETS)
                )
                for dataset_index, dataset in enumerate(REPORT.DATASETS):
                    numeric_extraction_by_dataset[dataset][
                        "component_by_linkage_basis_threshold"
                    ].setdefault(basis, {})[threshold] = numeric_component_cell(
                        distribution_parts[dataset_index]
                    )
        numeric_summary["extraction"] = {
            "overall_counts": copy.deepcopy(extraction["aggregate"]["counts"]),
            "overall_component_by_linkage_basis_threshold": overall_components,
            "by_dataset": numeric_extraction_by_dataset,
        }

        numeric_ranking_by_dataset = {}
        dataset_candidate_cells = {}
        for dataset in REPORT.DATASETS:
            nested = ranking["by_dataset"][dataset][
                "by_linkage_basis_and_threshold"
            ]
            display_nested = {}
            for basis in ("PS_HP1", "PS_HP2"):
                display_nested[basis] = {}
                for threshold in ("1", "2", "3", "5"):
                    raw = nested[basis][threshold]
                    metadata = candidate_cell_metadata[(dataset, basis, threshold)]
                    candidate = numeric_candidate_cell(
                        raw,
                        h_zero=metadata["h_distribution"].get("0", 0),
                        tied_consistent=metadata["tied_consistent"],
                        tree_partition_counts=metadata["tree_partition_counts"],
                    )
                    dataset_candidate_cells[(dataset, basis, threshold)] = candidate
                    display_nested[basis][threshold] = numeric_rank_cell(
                        raw, candidate=candidate
                    )
            numeric_ranking_by_dataset[dataset] = {
                "input_call_funnel": {},
                "input_sparse_call_code_counts": {},
                "input_hp_family_rows": {},
                "by_HP_basis_and_bridge_threshold": display_nested,
                "runtime_all_primary_unit_evaluations": {},
                "runtime_by_HP_basis_and_bridge_threshold": {},
            }

        overall_rank = {}
        aggregate_cells = ranking["aggregate"]["by_linkage_basis_and_threshold"]
        for basis in ("PS_HP1", "PS_HP2"):
            overall_rank[basis] = {}
            for threshold in ("1", "2", "3", "5"):
                raw = aggregate_cells[basis][threshold]
                parts = [
                    dataset_candidate_cells[(dataset, basis, threshold)]
                    for dataset in REPORT.DATASETS
                ]
                h_zero = sum(part["h_star_distribution"]["0"] for part in parts)
                tied_consistent = sum(
                    part["tied_by_coarse_topology"]["consistent"] for part in parts
                )
                tree_partition_counts = {
                    key: sum(part["tree_vertex_partition"]["counts"][key] for part in parts)
                    for key in REPORT.TREE_VERTEX_BUCKETS
                }
                candidate = numeric_candidate_cell(
                    raw, h_zero=h_zero, tied_consistent=tied_consistent,
                    tree_partition_counts=tree_partition_counts,
                )
                overall_rank[basis][threshold] = numeric_rank_cell(
                    raw, candidate=candidate
                )
        all_raw_cells = [
            aggregate_cells[basis][threshold]
            for basis in ("PS_HP1", "PS_HP2")
            for threshold in ("1", "2", "3", "5")
        ]
        numeric_summary["ranking"] = {
            "overall_input_call_funnel": {},
            "overall_input_sparse_call_code_counts": {},
            "overall_by_HP_basis_and_bridge_threshold": overall_rank,
            "overall_candidate_table": numeric_candidate_cell(
                REPORT._sum_count_trees(all_raw_cells)
            ),
            "overall_runtime_all_primary_unit_evaluations": copy.deepcopy(
                ranking["runtime_diagnostics"]
            ),
            "by_dataset": numeric_ranking_by_dataset,
        }
        frozen_numeric_producer = SCRIPT.with_name(
            "build_final_numeric_summary.py"
        ).resolve()
        numeric_producer = self.root / "live_build_final_numeric_summary.py"
        numeric_producer.write_bytes(frozen_numeric_producer.read_bytes())
        self.assertNotEqual(numeric_producer.resolve(), frozen_numeric_producer)
        numeric_summary["source_ledger"] = {
            "producer": {
                "path": str(numeric_producer),
                "size_bytes": numeric_producer.stat().st_size,
                "sha256": sha256(numeric_producer),
            },
            "final_independent_verification": {
                "path": str(paths["m2_verification"].resolve()),
                "size_bytes": paths["m2_verification"].stat().st_size,
                "sha256": sha256(paths["m2_verification"]),
            },
            "terminal_extraction": {
                "path": str(paths["m2_extraction"].resolve()),
                "size_bytes": paths["m2_extraction"].stat().st_size,
                "sha256": sha256(paths["m2_extraction"]),
            },
            "terminal_ranking": {
                "path": str(paths["m2_ranking"].resolve()),
                "size_bytes": paths["m2_ranking"].stat().st_size,
                "sha256": sha256(paths["m2_ranking"]),
            },
            "candidate_table": {
                "path": str(candidate_path.resolve()),
                "size_bytes": candidate_path.stat().st_size,
                "sha256": sha256(candidate_path),
            },
            "extraction_children": [
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "path": str(extraction_children[(dataset, chrom)][0].resolve()),
                    "size_bytes": extraction_children[(dataset, chrom)][0].stat().st_size,
                    "sha256": sha256(extraction_children[(dataset, chrom)][0]),
                }
                for dataset in REPORT.DATASETS
                for chrom in REPORT.AUTOSOMES
            ],
            "ranking_children": [],
            "runtime_files": [],
        }
        paths["m2_numeric_summary"] = self.root / "m2_numeric_summary.json"
        numeric_summary["receipt_integrity"] = {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": "m2_numeric_summary.json.sha256",
        }
        paths["m2_numeric_summary"].write_text(
            json.dumps(numeric_summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        paths["m2_numeric_summary_sidecar"] = self.root / "m2_numeric_summary.json.sha256"
        paths["m2_numeric_summary_sidecar"].write_text(
            f'{sha256(paths["m2_numeric_summary"])}  m2_numeric_summary.json\n',
            encoding="ascii",
        )
        paths["method"] = self.root / "partial-read-method-audit.md"
        paths["method"].write_text("# Fixture method audit\nSymbolic group; no first-success.\n", encoding="utf-8")
        paths["literature"] = self.root / "primary-literature-audit.md"
        paths["literature"].write_text("# Fixture literature audit\nClaim ceiling only.\n", encoding="utf-8")
        return paths, bundle

    @staticmethod
    def kwargs(paths, output):
        return {
            "canonical_json": paths["canonical"],
            "funnel_receipt": paths["funnel"],
            "funnel_verification_receipt": paths["funnel_verification"],
            "m0_receipt": paths["m0"],
            "m0_verification_receipt": paths["m0_verification"],
            "pilot_receipt": paths["pilot"],
            "method_audit": paths["method"],
            "literature_audit": paths["literature"],
            "m2_extraction_receipt": paths["m2_extraction"],
            "m2_ranking_receipt": paths["m2_ranking"],
            "m2_verification_receipt": paths["m2_verification"],
            "m2_numeric_summary": paths["m2_numeric_summary"],
            "m2_pilot_extraction_receipt": paths["m2_pilot"],
            "output_path": output,
        }

    @staticmethod
    def rewrite_numeric_summary(paths, payload) -> None:
        paths["m2_numeric_summary"].write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        paths["m2_numeric_summary_sidecar"].write_text(
            f'{sha256(paths["m2_numeric_summary"])}  {paths["m2_numeric_summary"].name}\n',
            encoding="ascii",
        )

    def test_fixture_values_flow_to_final_html_without_hand_copied_kpis(self) -> None:
        # Mutate fixture values while preserving invariants.  Assertions are
        # derived from the mutated receipts, so a hard-coded report fails.
        bundle = copy.deepcopy(self.bundle)
        bundle["canonical"]["canonical"]["aggregate"]["autosomal_biallelic_sSNV"] += 3
        bundle["canonical"]["canonical"]["samples"][0]["autosomal_biallelic_sSNV"] += 3
        bundle["m2_extraction"]["aggregate"]["counts"]["raw_overlapping_alignments"] += 7
        cell = bundle["m2_ranking"]["aggregate"]["by_linkage_basis_and_threshold"]["primary_hp"]["3"]
        cell["raw_tree_candidates_T_complete_units"] += 11
        paths, materialized = self.materialize(bundle)
        output = self.root / REPORT.FINAL_REPORT_NAME
        result = REPORT.build_report(**self.kwargs(paths, output))
        self.assertTrue(result.final_ready)
        self.assertEqual(result.output_path, output)
        document = output.read_text(encoding="utf-8")
        expected_ssnv = materialized["canonical"]["canonical"]["aggregate"]["autosomal_biallelic_sSNV"]
        expected_raw = materialized["m2_extraction"]["aggregate"]["counts"]["raw_overlapping_alignments"]
        expected_t = sum(
            materialized["m2_ranking"]["aggregate"]["by_linkage_basis_and_threshold"][basis]["3"]["raw_tree_candidates_T_complete_units"]
            for basis in ("PS_HP1", "PS_HP2")
        )
        self.assertIn(f"{expected_ssnv:,}", document)
        self.assertIn(f"{expected_raw:,}", document)
        self.assertIn(f"{expected_t:,}", document)
        self.assertIn("總分母：", document)
        self.assertIn("相對分母：", document)
        self.assertIn("7 datasets ≠ 7 biological samples", document)
        self.assertIn("1-1/2-1", document)
        self.assertIn("somatic sub-haplotype", document)
        self.assertIn("HP_family×known_PS×read_linked_component×threshold", document)
        self.assertIn("教授版 display convention：bridge ≥3 unique molecules", document)
        self.assertIn("這個 ≥3 不是 frozen primary bridge threshold", document)
        self.assertIn("PS_HP1", document)
        self.assertIn("PS_HP2", document)
        self.assertIn("Component-site k&gt;12", document)
        self.assertIn("不等於 effective solver k", document)
        self.assertIn("symmetric-substitution conditional R/A model", document)
        self.assertIn("conditional-on-fixed-candidate-set ranking stability", document)
        self.assertIn("獨立 verifier 已從 pattern、BQ、state 與 π 重算 likelihood", document)
        self.assertIn("candidate-bearing primary unit", document)
        self.assertIn("molecule projection", document)
        self.assertIn("max |Δ profile LL|", document)
        self.assertIn("topology_class_inclusion_counts", document)
        self.assertIn("coarse_topology_unique_class_counts", document)
        self.assertIn("Unique patterns", document)
        self.assertIn("BQ-quality groups", document)
        self.assertIn("Molecule projections", document)
        self.assertIn("Structural unique R/A/X patterns", document)
        self.assertIn("2ᵘ=4", document)
        self.assertIn("Structural evidence 與 scoring evidence 分開守恆", document)
        self.assertIn("Solver complete/incomplete 與 effective-k 路由", document)
        self.assertIn("Partial-read conceptual 與 effective completions", document)
        self.assertIn("Exact topology 與 parent-choice uniqueness", document)
        self.assertIn("NOT_EVALUATED_CANONICAL_SHAPE_ISOMORPHISM", document)
        self.assertIn("Ranking runtime 與 resource-conflict attestation", document)
        self.assertIn("12 <small>zero-conflict launch attestation", document)
        self.assertIn("不是 peak RSS", document)
        self.assertIn("Authenticated final numeric summary：每 dataset × HP × bridge", document)
        self.assertIn("Extraction component 與 k=1/k&gt;1（8 overall + 56 dataset strata）", document)
        self.assertIn("完整 T/V、effective route、h*、ranking 與 topology ledger", document)
        self.assertIn("T_EQ_1_V_EQ_1", document)
        self.assertIn("EXACT_COMPLETE unit", document)
        self.assertIn("OVERALL · 7 datasets", document)
        structure_ledger = document.split(
            "完整 T/V、effective route、h*、ranking 與 topology ledger", 1
        )[1].split("</details>", 1)[0]
        self.assertEqual(structure_ledger.count("OVERALL · 7 datasets"), 8)
        zero_bucket = (
            '<span class="metric-main">0 <small>T_EQ_1_V_EQ_1 unit</small></span>'
            '<span class="denom">總分母：7-dataset PS_HP1 bridge≥1 '
            'solver-complete units 49（0.00%）</span>'
        )
        self.assertIn(zero_bucket, structure_ledger)  # 0/49 = 0%, not N/A.
        self.assertNotIn(
            '<span class="metric-main">0 <small>T_EQ_1_V_EQ_1 unit</small></span>'
            '<span class="denom">總分母：7-dataset PS_HP1 bridge≥1 '
            'solver-complete units 0（不適用）</span>',
            structure_ledger,
        )
        self.assertIn(
            '<span class="metric-main">7 <small>single-only unit</small></span>'
            '<span class="denom">總分母：同 HP×threshold 的 7-dataset '
            'topology denominator 35（20.00%）</span>',
            structure_ledger,
        )
        exact_unique = (
            '<span class="metric-main">21 <small>exact-topology proven unique unit</small></span>'
            '<span class="denom">總分母：7-dataset PS_HP1 bridge≥1 '
            'topology-evaluated units 49（42.86%）</span>'
        )
        self.assertIn(exact_unique, structure_ledger)
        self.assertNotIn(
            '<span class="metric-main">21 <small>exact-topology proven unique unit</small></span>'
            '<span class="denom">總分母：7-dataset PS_HP1 bridge≥1 '
            'exact-topology unique units 21（100.00%）</span>',
            structure_ledger,
        )
        self.assertNotIn("該 dataset×HP×threshold", structure_ledger)
        for route in (
            "EXACT_COMPLETE", "EXACT_INCOMPLETE", "NOT_RUN_NO_STRUCTURE",
            "GT_EXACT_LIMIT",
        ):
            self.assertEqual(structure_ledger.count(f"{route} unit"), 64)
        for topology_class in REPORT.TOPOLOGY_CLASS_ORDER:
            self.assertGreaterEqual(structure_ledger.count(f"{topology_class} unit"), 64)
        self.assertIn("目前不可由正式證據辨識的量：null + reason", document)
        self.assertIn("未估計（null）", document)
        self.assertIn("threshold strata reuse evidence", document)
        self.assertNotIn(
            "各 dataset × PS_HP1/PS_HP2 的 solver、likelihood、T/V 與 topology",
            document,
        )
        self.assertIn(sha256(paths["m2_numeric_summary"]), document)
        self.assertIn("latent expected proportions；不是 cell clone assignment", document)
        self.assertIn(materialized["candidate_rows"][0]["unit_key"], document)
        self.assertIn("322 / 322 checks PASS", document)
        self.assertIn(sha256(paths["funnel_verification"]), document)
        candidate_examples = document.split(
            "實際區域候選：hashed candidate table 的 bounded examples", 1
        )[1].split("</details>", 1)[0]
        n_candidate_rows = len(materialized["candidate_rows"])
        n_candidate_units = len(
            {row["unit_key"] for row in materialized["candidate_rows"]}
        )
        self.assertIn(f"{n_candidate_rows:,} <small>candidate rows", candidate_examples)
        self.assertIn(f"{n_candidate_units:,} <small>PS-aware units", candidate_examples)
        self.assertIn(
            f"{REPORT._n(n_candidate_rows / n_candidate_units)} <small>rows / unit",
            candidate_examples,
        )
        self.assertIn("parent choices", candidate_examples)
        self.assertNotIn("%", candidate_examples.split("<h3>Minimum-extra-state", 1)[0])
        self.assertIn("Minimum-extra-state count h*", candidate_examples)
        self.assertIn("partial-compatible representative", candidate_examples)
        self.assertIn("不是 hidden clone 數", candidate_examples)
        self.assertIn("h*+1、h*+2", candidate_examples)
        self.assertIn("parent_choice_count", candidate_examples)
        self.assertIn("並未在 canonical candidate table 展開", candidate_examples)
        self.assertIn(
            "7 datasets / 6 biological samples：raw molecule、PS、primary-unit inclusion",
            document,
        )
        self.assertIn("HCC1395 的 Dorado pipeline dataset；非獨立生物樣本", document)
        self.assertIn("把 T、V 與 likelihood 排名拆成互斥層級", document)
        self.assertIn("各 dataset 的 T/V、likelihood unique/tie 與 abstain", document)
        self.assertIn("T&gt;1、V&gt;1：likelihood 並列第一 vertex sets", document)
        self.assertIn("名詞與符號完整定義", document)
        self.assertIn("C / structural pattern groups", document)
        self.assertIn("N / hidden-extra vertex", document)
        self.assertIn("C 不是 clone 數", document)
        self.assertIn("分開處理不是因為ONT不能讀chrX", document)
        self.assertIn("50 kb 是 positional heuristic，不是 ONT 物理讀長上限", document)
        self.assertIn("舊 MAX_SNV=8 是計算 cap", document)
        self.assertIn("需另做 genomic matched-unit 分析", document)
        self.assertIn("eligible alignments by raw HP value（含 HP=. missing）", document)
        self.assertIn("明確包含 missing <code>HP=.</code>", document)
        self.assertNotIn("39,885", document)
        self.assertNotIn("48,963", document)
        self.assertNotIn("PARTIAL preview 不可作論文", document)
        self.assertIn("FINAL artifact 已通過本報告宣告的 Task-B receipt gates", document)
        self.assertNotIn("https://", document)
        self.assertNotIn("<script", document.lower())
        self.assertGreaterEqual(document.count("<svg"), 6)
        self.assertIn(sha256(paths["canonical"]), document)
        probe = HtmlProbe()
        probe.feed(document)
        self.assertEqual(probe.lang, "zh-Hant")
        self.assertGreaterEqual(probe.tags.get("section", 0), 8)
        self.assertGreaterEqual(probe.tags.get("details", 0), 5)
        self.assertEqual(len(probe.ids), len(set(probe.ids)), "HTML/SVG id 必須唯一")
        self.assertIn("Read-linked Hypercube", "".join(probe.title_text))

    def test_production_shape_ranking_cells_merge_counts_and_contracts_semantically(self) -> None:
        paths, _ = self.materialize()
        ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
        aggregate_cell = ranking["aggregate"]["by_linkage_basis_and_threshold"]["PS_HP1"]["3"]
        dataset_cells = [
            ranking["by_dataset"][dataset]["by_linkage_basis_and_threshold"]["PS_HP1"]["3"]
            for dataset in REPORT.DATASETS
        ]
        merged = REPORT._sum_count_trees(dataset_cells)
        self.assertEqual(merged, aggregate_cell)
        self.assertEqual(
            merged["partial_pattern_funnel"]["definitions"]["conceptual"],
            aggregate_cell["partial_pattern_funnel"]["definitions"]["conceptual"],
        )
        self.assertTrue(merged["all_conserved"])
        self.assertTrue(all(merged["conservation_checks"].values()))
        self.assertEqual(
            set(merged["partial_pattern_funnel"]),
            {
                "definitions",
                *REPORT.PARTIAL_FUNNEL_GRAINS,
                "units_denominator",
                "units_with_partial_structural_groups",
                "partial_group_coverage_denominator",
                "partial_groups_covered",
                "partial_groups_unsatisfied",
            },
        )
        self.assertEqual(
            ranking["runtime_diagnostics"]["n_unit_evaluations"],
            sum(
                ranking["aggregate"]["by_linkage_basis_and_threshold"][basis][threshold][
                    "n_component_hp_units"
                ]
                for basis in ("PS_HP1", "PS_HP2")
                for threshold in ("1", "2", "3", "5")
            ),
        )
        self.assertEqual(
            merged["parent_edge_assignment_unique_units"],
            merged["exact_topology_proven_unique_units"],
        )

        self.assertEqual(
            REPORT._sum_count_trees(
                [
                    {"count": 2, "definition": "fixed", "checks": {"ok": True}},
                    {"count": 3, "definition": "fixed", "checks": {"ok": False}},
                ]
            ),
            {"count": 5, "definition": "fixed", "checks": {"ok": False}},
        )
        with self.assertRaisesRegex(REPORT.ReportGateError, "invariant definition 不一致"):
            REPORT._sum_count_trees(
                [
                    {"count": 1, "definition": "contract-a"},
                    {"count": 1, "definition": "contract-b"},
                ]
            )

    def test_numeric_summary_sidecar_and_cli_argument_are_required(self) -> None:
        paths, _ = self.materialize()
        paths["m2_numeric_summary_sidecar"].write_text(
            f'{"0" * 64}  {paths["m2_numeric_summary"].name}\n', encoding="ascii"
        )
        with self.assertRaisesRegex(REPORT.ReportGateError, "sidecar 未精確綁定"):
            REPORT.build_report(
                **self.kwargs(paths, self.root / "bad_numeric_sidecar.html"),
                allow_partial=True,
            )

        paths, _ = self.materialize()
        command = [
            sys.executable,
            str(SCRIPT),
            "--canonical-json", str(paths["canonical"]),
            "--pilot-receipt", str(paths["pilot"]),
            "--method-audit", str(paths["method"]),
            "--literature-audit", str(paths["literature"]),
            "--output", str(self.root / "missing_numeric_cli.html"),
            "--allow-partial",
        ]
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("--m2-numeric-summary", completed.stderr)

    def test_numeric_summary_schema_scope_and_source_bindings_fail_closed(self) -> None:
        cases = (
            (
                "schema",
                lambda payload: payload.__setitem__("schema_version", "0.9.0"),
                "schema_version 不符",
            ),
            (
                "all_pass",
                lambda payload: payload.__setitem__("all_pass", False),
                "all_pass 不是 true",
            ),
            (
                "scope",
                lambda payload: payload["scope"].__setitem__("n_technical_datasets", 6),
                "canonical 7-dataset×chr1-22 scope",
            ),
            (
                "check",
                lambda payload: payload["checks"].__setitem__(
                    "all_component_k_distributions_conserve", False
                ),
                "checks 不完整或未全數 true",
            ),
            (
                "producer_size",
                lambda payload: payload["source_ledger"]["producer"].__setitem__(
                    "size_bytes", 1
                ),
                "producer size_bytes 與 recorded live producer 不符",
            ),
            (
                "producer_sha",
                lambda payload: payload["source_ledger"]["producer"].__setitem__(
                    "sha256", "0" * 64
                ),
                "producer SHA-256 與 recorded live producer 不符",
            ),
            (
                "extraction_binding",
                lambda payload: payload["source_ledger"]["terminal_extraction"].__setitem__(
                    "sha256", "0" * 64
                ),
                "terminal_extraction SHA-256 未精確綁定",
            ),
            (
                "ranking_binding",
                lambda payload: payload["source_ledger"]["terminal_ranking"].__setitem__(
                    "path", "/fixture/wrong-ranking.json"
                ),
                "terminal_ranking path 未綁定",
            ),
            (
                "verifier_binding",
                lambda payload: payload["source_ledger"][
                    "final_independent_verification"
                ].__setitem__("size_bytes", 0),
                "final_independent_verification size_bytes 不符",
            ),
        )
        for label, mutate, expected_error in cases:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                payload = json.loads(
                    paths["m2_numeric_summary"].read_text(encoding="utf-8")
                )
                mutate(payload)
                self.rewrite_numeric_summary(paths, payload)
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_numeric_{label}.html"),
                        allow_partial=True,
                    )

    def test_numeric_summary_live_producer_must_match_frozen_sibling_bytes(self) -> None:
        paths, _ = self.materialize()
        payload = json.loads(
            paths["m2_numeric_summary"].read_text(encoding="utf-8")
        )
        identity = payload["source_ledger"]["producer"]
        live_producer = Path(identity["path"])
        frozen_sibling = SCRIPT.with_name("build_final_numeric_summary.py").resolve()
        self.assertNotEqual(live_producer.resolve(), frozen_sibling)
        live_producer.write_bytes(live_producer.read_bytes() + b"# producer drift\n")
        # Preserve ledger↔live self-consistency: only the required
        # live↔presentation-snapshot byte equivalence should reject this.
        identity["size_bytes"] = live_producer.stat().st_size
        identity["sha256"] = sha256(live_producer)
        self.rewrite_numeric_summary(paths, payload)
        with self.assertRaisesRegex(
            REPORT.ReportGateError,
            "recorded live producer bytes .* frozen sibling 不符",
        ):
            REPORT.build_report(
                **self.kwargs(paths, self.root / "bad_numeric_producer_drift.html"),
                allow_partial=True,
            )

    def test_numeric_summary_component_hstar_tied_topology_and_null_fail_closed(self) -> None:
        cases = (
            (
                "component_k",
                lambda payload: payload["extraction"]["by_dataset"]["COLO829"][
                    "component_by_linkage_basis_threshold"
                ]["PS_HP1"]["3"].__setitem__("k_equals_1", 99),
                "k=1 不符 distribution|k=1\+k>1",
            ),
            (
                "effective_route",
                lambda payload: payload["ranking"]["by_dataset"]["H1437"][
                    "by_HP_basis_and_bridge_threshold"
                ]["PS_HP2"]["2"]["effective_k"]["k_route_counts"].__setitem__(
                    "EXACT_COMPLETE", 99
                ),
                "effective-k route partition 不守恆",
            ),
            (
                "h_star",
                lambda payload: payload["ranking"]["by_dataset"]["H2009"][
                    "by_HP_basis_and_bridge_threshold"
                ]["PS_HP1"]["5"]["candidate_structure"]["candidate_table"][
                    "h_star_distribution"
                ].__setitem__("1", 99),
                "h\* distribution 不守恆",
            ),
            (
                "tied_topology",
                lambda payload: payload["ranking"]["by_dataset"]["HCC1395"][
                    "by_HP_basis_and_bridge_threshold"
                ]["PS_HP2"]["1"]["candidate_structure"]["candidate_table"][
                    "tied_by_coarse_topology"
                ].__setitem__("consistent", 99),
                "tied×coarse-Topo partition 不守恆",
            ),
            (
                "unsupported_null",
                lambda payload: payload["unsupported_or_nonidentifiable"][
                    "cellular_HP1_HP2_clone_pairing"
                ].__setitem__("value", 0),
                "unsupported.cellular_HP1_HP2_clone_pairing 未以 null\+reason",
            ),
        )
        for label, mutate, expected_error in cases:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                payload = json.loads(
                    paths["m2_numeric_summary"].read_text(encoding="utf-8")
                )
                mutate(payload)
                self.rewrite_numeric_summary(paths, payload)
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_numeric_data_{label}.html"),
                        allow_partial=True,
                    )

    def test_s16_t_v_partition_mean_and_missing_fields_fail_closed(self) -> None:
        def mutate_partition(payload):
            candidate = payload["ranking"]["by_dataset"]["COLO829"][
                "by_HP_basis_and_bridge_threshold"
            ]["PS_HP1"]["3"]["candidate_structure"]["candidate_table"]
            partition = candidate["tree_vertex_partition"]
            source = next(
                key for key, value in partition["counts"].items() if value > 0
            )
            target = next(key for key in REPORT.TREE_VERTEX_BUCKETS if key != source)
            partition["counts"][source] -= 1
            partition["counts"][target] += 1
            for key, value in partition["counts"].items():
                partition["shares"][key] = ratio_payload(
                    value, partition["denominator"],
                    "solver_complete_candidate_units",
                )

        cases = (
            (
                "missing_T",
                lambda payload: payload["ranking"]["by_dataset"]["COLO829"]
                ["by_HP_basis_and_bridge_threshold"]["PS_HP1"]["3"]
                ["candidate_structure"]["candidate_table"].pop(
                    "n_parent_edge_trees_T"
                ),
                "candidate table keys 不完整|candidate count 無效",
            ),
            (
                "partition",
                mutate_partition,
                "tree_vertex_partition 未精確 reconcile S9 stream",
            ),
            (
                "mean_T",
                lambda payload: payload["ranking"]["by_dataset"]["H1437"]
                ["by_HP_basis_and_bridge_threshold"]["PS_HP2"]["5"]
                ["candidate_structure"]["mean_T_per_solver_complete_unit"]
                .__setitem__("value", 999.0),
                "mean T ratio value 不可重算",
            ),
        )
        for label, mutate, expected in cases:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                payload = json.loads(
                    paths["m2_numeric_summary"].read_text(encoding="utf-8")
                )
                mutate(payload)
                self.rewrite_numeric_summary(paths, payload)
                with self.assertRaisesRegex(REPORT.ReportGateError, expected):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_s16_{label}.html"),
                        allow_partial=True,
                    )

    def test_s16_dataset_route_t_and_topology_exchange_fail_closed(self) -> None:
        # Route exchange: preserve both dataset unit totals and the 7-dataset
        # total, while exchanging one complete/incomplete route assignment.
        paths, _ = self.materialize()
        payload = json.loads(paths["m2_numeric_summary"].read_text(encoding="utf-8"))
        first = payload["ranking"]["by_dataset"]["COLO829"][
            "by_HP_basis_and_bridge_threshold"
        ]["PS_HP1"]["3"]
        second = payload["ranking"]["by_dataset"]["H1437"][
            "by_HP_basis_and_bridge_threshold"
        ]["PS_HP1"]["3"]
        first["effective_k"]["k_route_counts"]["EXACT_COMPLETE"] -= 1
        first["effective_k"]["k_route_counts"]["EXACT_INCOMPLETE"] += 1
        second["effective_k"]["k_route_counts"]["EXACT_COMPLETE"] += 1
        second["effective_k"]["k_route_counts"]["EXACT_INCOMPLETE"] -= 1
        for cell in (first, second):
            denominator = cell["units"]["n_component_hp_unit_evaluations"]
            cell["effective_k"]["route_shares_of_unit_evaluations"] = {
                key: ratio_payload(value, denominator, "component_hp_unit_evaluations")
                for key, value in cell["effective_k"]["k_route_counts"].items()
            }
        self.rewrite_numeric_summary(paths, payload)
        with self.assertRaisesRegex(
            REPORT.ReportGateError,
            "EXACT_COMPLETE != solver complete|未精確 reconcile S7 raw cell",
        ):
            REPORT.build_report(
                **self.kwargs(paths, self.root / "bad_route_exchange.html"),
                allow_partial=True,
            )

        # T exchange: construct unequal COLO829/H1437 T cells, then swap every
        # duplicated S16 T/mean field.  Aggregate arithmetic still conserves;
        # only the independent S9 per-stratum stream catches the route swap.
        bundle = copy.deepcopy(self.bundle)
        raw = bundle["m2_ranking"]["aggregate"][
            "by_linkage_basis_and_threshold"
        ]["primary_hp"]["3"]
        raw["raw_tree_candidates_T_complete_units"] += 2
        paths, _ = self.materialize(bundle)
        payload = json.loads(paths["m2_numeric_summary"].read_text(encoding="utf-8"))
        first = payload["ranking"]["by_dataset"]["COLO829"][
            "by_HP_basis_and_bridge_threshold"
        ]["PS_HP1"]["3"]["candidate_structure"]
        second = payload["ranking"]["by_dataset"]["H1437"][
            "by_HP_basis_and_bridge_threshold"
        ]["PS_HP1"]["3"]["candidate_structure"]
        self.assertNotEqual(
            first["raw_parent_edge_trees_T_complete_units"],
            second["raw_parent_edge_trees_T_complete_units"],
        )
        for key in (
            "raw_parent_edge_trees_T_complete_units",
            "mean_T_per_solver_complete_unit",
        ):
            first[key], second[key] = second[key], first[key]
        first["candidate_table"]["n_parent_edge_trees_T"], second[
            "candidate_table"
        ]["n_parent_edge_trees_T"] = (
            second["candidate_table"]["n_parent_edge_trees_T"],
            first["candidate_table"]["n_parent_edge_trees_T"],
        )
        self.rewrite_numeric_summary(paths, payload)
        with self.assertRaisesRegex(
            REPORT.ReportGateError,
            "candidate.n_parent_edge_trees_T 未精確 reconcile S9 stream",
        ):
            REPORT.build_report(
                **self.kwargs(paths, self.root / "bad_T_exchange.html"),
                allow_partial=True,
            )

        # Topology exchange preserves both per-cell topology denominators and
        # the HP×threshold overall class totals.
        paths, _ = self.materialize()
        payload = json.loads(paths["m2_numeric_summary"].read_text(encoding="utf-8"))
        cells = [
            payload["ranking"]["by_dataset"][dataset]
            ["by_HP_basis_and_bridge_threshold"]["PS_HP1"]["3"]
            for dataset in ("COLO829", "H1437")
        ]
        for location in (
            lambda cell: cell["topology"]["coarse_unique_class_counts"],
            lambda cell: cell["candidate_structure"]["candidate_table"]
            ["topology"]["coarse_unique_class_counts"],
        ):
            first_map, second_map = location(cells[0]), location(cells[1])
            first_map["single-only"] += 1
            first_map["direct-only"] -= 1
            second_map["single-only"] -= 1
            second_map["direct-only"] += 1
        self.rewrite_numeric_summary(paths, payload)
        with self.assertRaisesRegex(
            REPORT.ReportGateError,
            "candidate topology.coarse_unique_class_counts 未精確 reconcile S9 stream|未精確 reconcile S7 raw cell",
        ):
            REPORT.build_report(
                **self.kwargs(paths, self.root / "bad_topology_exchange.html"),
                allow_partial=True,
            )

    def test_s13_physical_release_manifest_sidecar_and_snapshot_tamper_fail_closed(self) -> None:
        for label in ("manifest", "sidecar", "snapshot"):
            with self.subTest(label=label):
                paths, _ = self.materialize()
                verification = json.loads(
                    paths["m2_verification"].read_text(encoding="utf-8")
                )
                binding = verification["release_binding"]
                if label == "manifest":
                    target = Path(binding["release_manifest"]["path"])
                elif label == "sidecar":
                    target = Path(binding["release_manifest"]["sidecar"]["path"])
                else:
                    target = Path(
                        binding["snapshot_sources"]["extractor"]["path"]
                    )
                target.chmod(0o644)
                target.write_bytes(target.read_bytes() + b"tamper")
                with self.assertRaisesRegex(
                    REPORT.ReportGateError,
                    "release manifest (?:physical mode|physical SHA-256|JSON)|release sidecar|snapshot/extractor physical",
                ):
                    REPORT.build_report(
                        **self.kwargs(
                            paths, self.root / f"bad_release_{label}.html"
                        ),
                        allow_partial=True,
                    )

    def test_production_shape_partial_solver_and_topology_drift_fail_closed(self) -> None:
        cases = (
            (
                "conceptual_mass",
                lambda cell: cell["partial_pattern_funnel"]["molecule_projections"].__setitem__(
                    "conceptual_completions_weighted_total",
                    cell["partial_pattern_funnel"]["molecule_projections"][
                        "conceptual_completions_weighted_total"
                    ]
                    + 1,
                ),
                "conceptual weighted total",
            ),
            (
                "effective_distribution",
                lambda cell: cell["partial_pattern_funnel"]["unique_rax_pattern_groups"][
                    "observed_alt_effective_completions_distribution"
                ].__setitem__(
                    "1",
                    cell["partial_pattern_funnel"]["unique_rax_pattern_groups"][
                        "observed_alt_effective_completions_distribution"
                    ]["1"]
                    + 1,
                ),
                "effective distribution 不守恆",
            ),
            (
                "structural_scoring",
                lambda cell: cell.__setitem__(
                    "structural_retained_molecules",
                    cell["structural_retained_molecules"] + 1,
                ),
                "structural/scoring funnel 不守恆",
            ),
            (
                "effective_k_mass",
                lambda cell: cell.__setitem__(
                    "k_observed_alt_active_total",
                    cell["k_observed_alt_active_total"] + 1,
                ),
                "effective-k site mass 不守恆",
            ),
            (
                "exact_topology_parent_choice",
                lambda cell: cell.__setitem__(
                    "parent_edge_assignment_unique_units",
                    cell["parent_edge_assignment_unique_units"] + 1,
                ),
                "parent-edge unique != exact-topology proven unique",
            ),
        )
        for label, mutate, expected_error in cases:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
                cell = ranking["aggregate"]["by_linkage_basis_and_threshold"]["PS_HP1"]["3"]
                mutate(cell)
                paths["m2_ranking"].write_text(
                    json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8"
                )
                verification = json.loads(
                    paths["m2_verification"].read_text(encoding="utf-8")
                )
                verification["ranking"]["receipt_sha256"] = sha256(paths["m2_ranking"])
                paths["m2_verification"].write_text(
                    json.dumps(verification, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_{label}.html"),
                        allow_partial=True,
                    )

    def test_production_shape_runtime_contract_drift_fails_closed(self) -> None:
        cases = (
            (
                "all_unit_n",
                lambda ranking: ranking["runtime_diagnostics"]["metrics"][
                    "unit_total_elapsed_seconds"
                ].__setitem__("n", 447),
                "runtime n 不符 expected units",
            ),
            (
                "invoked_n",
                lambda ranking: ranking["runtime_diagnostics"]["metrics_when_invoked"][
                    "candidate_generation_elapsed_seconds"
                ].__setitem__("n", 449),
                "n > all units",
            ),
            (
                "required_check",
                lambda ranking: ranking["checks"].pop(
                    "full_runtime_aggregate_recomputed_from_streamed_unit_rows"
                ),
                "checks.full_runtime_aggregate_recomputed_from_streamed_unit_rows",
            ),
        )
        for label, mutate, expected_error in cases:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
                mutate(ranking)
                paths["m2_ranking"].write_text(
                    json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8"
                )
                verification = json.loads(
                    paths["m2_verification"].read_text(encoding="utf-8")
                )
                verification["ranking"]["receipt_sha256"] = sha256(paths["m2_ranking"])
                paths["m2_verification"].write_text(
                    json.dumps(verification, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_runtime_{label}.html"),
                        allow_partial=True,
                    )

    def test_resource_conflict_session_attestation_fails_closed(self) -> None:
        cases = (
            (
                "nonzero_process",
                lambda gate: gate.__setitem__("process_count", 1),
                "zero-conflict attestation",
            ),
            (
                "bypass",
                lambda gate: gate.__setitem__("ignore_resource_gate", True),
                "resource gate 被 bypass",
            ),
        )
        for label, mutate, expected_error in cases:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                session = json.loads(
                    paths["m2_extraction_session"].read_text(encoding="utf-8")
                )
                gate = session["resource_gate"]
                mutate(gate)
                gate["process_snapshot_sha256"] = REPORT._semantic_json_sha256(
                    {
                        "process_count": gate["process_count"],
                        "root_count": gate["root_count"],
                        "representatives": gate["representatives"],
                    }
                )
                paths["m2_extraction_session"].write_text(
                    json.dumps(session, ensure_ascii=False) + "\n", encoding="utf-8"
                )
                verification = json.loads(
                    paths["m2_verification"].read_text(encoding="utf-8")
                )
                verification["extraction"]["orchestration"]["session_receipt"][
                    "sha256"
                ] = sha256(paths["m2_extraction_session"])
                paths["m2_verification"].write_text(
                    json.dumps(verification, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_resource_{label}.html"),
                        allow_partial=True,
                    )

    def test_cli_missing_full_m2_exits_nonzero_and_creates_no_final(self) -> None:
        paths, _ = self.materialize()
        final = self.root / REPORT.FINAL_REPORT_NAME
        command = [
            sys.executable,
            str(SCRIPT),
            "--canonical-json", str(paths["canonical"]),
            "--funnel-receipt", str(paths["funnel"]),
            "--funnel-verification-receipt", str(paths["funnel_verification"]),
            "--m0-receipt", str(paths["m0"]),
            "--m0-verification-receipt", str(paths["m0_verification"]),
            "--pilot-receipt", str(paths["pilot"]),
            "--method-audit", str(paths["method"]),
            "--literature-audit", str(paths["literature"]),
            "--m2-numeric-summary", str(paths["m2_numeric_summary"]),
            "--m2-pilot-extraction-receipt", str(paths["m2_pilot"]),
            "--output", str(final),
        ]
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        self.assertEqual(completed.returncode, 2, completed.stdout + completed.stderr)
        self.assertFalse(final.exists())
        self.assertIn("缺少 M2 full extraction receipt", completed.stdout)
        self.assertIn("缺少 M2 full ranking receipt", completed.stdout)
        self.assertNotIn("funnel independent", completed.stdout.lower())

    def test_allow_partial_redirects_final_name_and_adds_ribbon(self) -> None:
        paths, _ = self.materialize()
        final = self.root / REPORT.FINAL_REPORT_NAME
        kwargs = self.kwargs(paths, final)
        kwargs["m2_extraction_receipt"] = None
        kwargs["m2_ranking_receipt"] = None
        result = REPORT.build_report(**kwargs, allow_partial=True)
        self.assertFalse(result.final_ready)
        self.assertFalse(final.exists())
        self.assertIn("partial-preview", result.output_path.name)
        document = result.output_path.read_text(encoding="utf-8")
        self.assertIn("PARTIAL PREVIEW · NOT VALIDATION EVIDENCE", document)
        self.assertIn("PARTIAL preview 不可作論文、口試或 validation evidence", document)
        self.assertIn("M2 154-task extraction 尚未完成或未提供", document)
        self.assertIn("缺值保持「未評估」", document)
        self.assertIn("M0 legacy baseline 已完成", document)
        self.assertIn("P(p | N, π) = Σᵥ∈N πᵥ P(p | v, error)", document)
        self.assertIn("只有 count 達 <code>structural_exact_pattern_minread</code>", document)
        self.assertIn("每次 MILP construction", document)
        self.assertIn("相同 N、不同 parent edges", document)
        self.assertNotIn("P(p | V)", document)
        self.assertNotIn("列出全部 optimal V", document)
        self.assertNotIn("solver 只存一個 group", document)
        self.assertNotIn("PS-aware full M0", document)

    def test_report_publication_race_never_replaces_competing_file(self) -> None:
        paths, _ = self.materialize()
        output = self.root / "publication_race.html"
        original = REPORT._rename_noreplace

        def inject_competing_file(source: Path, destination: Path) -> None:
            destination.write_text("competing writer\n", encoding="utf-8")
            original(source, destination)

        REPORT._rename_noreplace = inject_competing_file
        try:
            with self.assertRaisesRegex(FileExistsError, "refusing to overwrite"):
                REPORT.build_report(**self.kwargs(paths, output))
        finally:
            REPORT._rename_noreplace = original
        self.assertEqual(output.read_text(encoding="utf-8"), "competing writer\n")
        staged = list(self.root.glob(f".{output.name}.*.tmp"))
        self.assertEqual(len(staged), 1)
        self.assertEqual(staged[0].stat().st_mode & 0o777, 0o444)

    def test_explicit_overwrite_preserves_superseded_report(self) -> None:
        paths, _ = self.materialize()
        output = self.root / "superseded.html"
        output.write_text("previous report\n", encoding="utf-8")
        result = REPORT.build_report(
            **self.kwargs(paths, output),
            overwrite=True,
        )
        self.assertTrue(result.final_ready)
        self.assertNotEqual(output.read_text(encoding="utf-8"), "previous report\n")
        self.assertEqual(output.stat().st_mode & 0o777, 0o444)
        archived = list(self.root.glob(f"{output.name}.superseded.*"))
        self.assertEqual(len(archived), 1)
        self.assertEqual(archived[0].read_text(encoding="utf-8"), "previous report\n")

    def test_cli_partial_generation_never_reports_all_pass(self) -> None:
        paths, _ = self.materialize()
        final = self.root / REPORT.FINAL_REPORT_NAME
        command = [
            sys.executable,
            str(SCRIPT),
            "--canonical-json", str(paths["canonical"]),
            "--m0-receipt", str(paths["m0"]),
            "--pilot-receipt", str(paths["pilot"]),
            "--method-audit", str(paths["method"]),
            "--literature-audit", str(paths["literature"]),
            "--m2-numeric-summary", str(paths["m2_numeric_summary"]),
            "--output", str(final),
            "--allow-partial",
        ]
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        receipt = json.loads(completed.stdout)
        self.assertTrue(receipt["generation_pass"])
        self.assertFalse(receipt["all_pass"])
        self.assertFalse(receipt["final_ready"])
        self.assertEqual(receipt["artifact_status"], "PARTIAL_PREVIEW_NOT_VALIDATION_EVIDENCE")
        self.assertFalse(final.exists())

    def test_pilot_receipt_binds_live_runner_and_research_solver_hashes(self) -> None:
        for key, expected_error in (
            ("runner_sha256", "runner live SHA-256 不符"),
            ("research_solver_sha256", "research solver live SHA-256 不符"),
        ):
            with self.subTest(key=key):
                paths, _ = self.materialize()
                pilot = json.loads(paths["pilot"].read_text(encoding="utf-8"))
                pilot["implementation"][key] = "0" * 64
                paths["pilot"].write_text(
                    json.dumps(pilot, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_pilot_{key}.html"),
                        allow_partial=True,
                    )

    def test_ranking_must_bind_exact_extraction_receipt_hash(self) -> None:
        paths, _ = self.materialize()
        ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
        ranking["upstream_extraction_receipt"]["sha256"] = "0" * 64
        paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8")
        final = self.root / REPORT.FINAL_REPORT_NAME
        with self.assertRaisesRegex(REPORT.ReportGateError, "extraction receipt SHA-256"):
            REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)
        self.assertFalse(final.exists())

    def test_ranking_status_and_solver_funnels_fail_closed(self) -> None:
        paths, _ = self.materialize()
        ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
        cell = ranking["aggregate"]["by_linkage_basis_and_threshold"]["PS_HP1"]["3"]
        cell["selection_status_counts"]["RANK_ABSTAIN"] += 1
        paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8")
        final = self.root / REPORT.FINAL_REPORT_NAME
        with self.assertRaisesRegex(REPORT.ReportGateError, "selection_status_counts != HP units"):
            REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)
        self.assertFalse(final.exists())

    def test_ps_isolation_check_failure_fails_closed(self) -> None:
        paths, _ = self.materialize()
        ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
        ranking["task_index"][0]["known_ps_never_mixed"] = False
        paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8")
        final = self.root / REPORT.FINAL_REPORT_NAME
        with self.assertRaisesRegex(REPORT.ReportGateError, "known PS"):
            REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)
        self.assertFalse(final.exists())

    def test_old_partial_key_only_cannot_satisfy_final_method_gate(self) -> None:
        paths, _ = self.materialize()
        ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
        ranking["run_contract"].pop("method_contract")
        ranking["checks"].pop(
            "all_154_child_method_contracts_identical_and_source_bound"
        )
        ranking["checks"]["no_" + "partial_completions_materialized"] = True
        paths["m2_ranking"].write_text(
            json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(REPORT.ReportGateError, "method contract|all_154_child"):
            REPORT.build_report(
                **self.kwargs(paths, self.root / REPORT.FINAL_REPORT_NAME),
                allow_partial=True,
            )

    def test_method_contract_drift_vaf_and_edge_scoring_fail_closed(self) -> None:
        for field in (
            "structural_group_scope",
            "same_molecule_vaf_added_as_separate_term",
            "parent_edges_or_transitions_scored",
        ):
            with self.subTest(field=field):
                paths, _ = self.materialize()
                ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
                ranking["run_contract"]["method_contract"][field] = (
                    "DRIFT" if field == "structural_group_scope" else True
                )
                paths["m2_ranking"].write_text(
                    json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8"
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, "method contract"):
                    REPORT.build_report(
                        **self.kwargs(
                            paths, self.root / f"bad_contract_{field}.html"
                        ),
                        allow_partial=True,
                    )

    def test_verifier_must_independently_recompute_profile_likelihood(self) -> None:
        paths, _ = self.materialize()
        verification = json.loads(paths["m2_verification"].read_text(encoding="utf-8"))
        verification["verification_independence"][
            "profile_likelihood_recomputed_from_patterns_states_pi"
        ] = False
        paths["m2_verification"].write_text(
            json.dumps(verification, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(REPORT.ReportGateError, "未獨立由 patterns/states/pi 重算"):
            REPORT.build_report(
                **self.kwargs(paths, self.root / "false_likelihood_claim.html"),
                allow_partial=True,
            )

    def test_final_requires_release_post_and_orchestration_evidence(self) -> None:
        mutations = (
            (
                "release",
                lambda payload: payload["release_binding"].__setitem__(
                    "validation_evidence_eligible", False
                ),
                "canonical release contract",
            ),
            (
                "post",
                lambda payload: payload["post_input_identity"].__setitem__(
                    "exact_snapshot_equal", False
                ),
                "POST.*frozen PRE",
            ),
            (
                "extraction_orphan",
                lambda payload: payload["extraction"]["orchestration"].__setitem__(
                    "no_open_orphan_or_preseeded_child_accepted", False
                ),
                "extraction orchestration.*no_open_orphan",
            ),
            (
                "ranking_chain",
                lambda payload: payload["ranking"]["orchestration"].__setitem__(
                    "n_batches", 10
                ),
                "ranking orchestration 8→16 batch/154-child chain",
            ),
        )
        for label, mutate, expected_error in mutations:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                verification = json.loads(
                    paths["m2_verification"].read_text(encoding="utf-8")
                )
                mutate(verification)
                paths["m2_verification"].write_text(
                    json.dumps(verification, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(
                            paths, self.root / f"bad_{label}_evidence.html"
                        ),
                        allow_partial=True,
                    )

    def test_candidate_table_mutation_after_receipt_fails_closed(self) -> None:
        paths, _ = self.materialize()
        candidate_path = self.root / "m2_ps_aware_candidate_table.tsv.gz"
        with candidate_path.open("ab") as handle:
            handle.write(b"post-receipt mutation")
        final = self.root / REPORT.FINAL_REPORT_NAME
        with self.assertRaisesRegex(REPORT.ReportGateError, "candidate table (?:size_bytes|SHA-256)"):
            REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)
        self.assertFalse(final.exists())

    def test_final_requires_funnel_m0_and_m2_independent_verifiers(self) -> None:
        paths, _ = self.materialize()
        missing_cases = (
            ("funnel_receipt", "current sSNV funnel"),
            ("funnel_verification_receipt", "funnel independent verification"),
            ("m0_verification_receipt", "M0 independent"),
            ("m2_verification_receipt", "M2 full independent"),
            ("m2_numeric_summary", "authenticated final numeric summary"),
        )
        for key, expected_issue in missing_cases:
            with self.subTest(key=key):
                output = self.root / f"missing_{key}.html"
                kwargs = self.kwargs(paths, output)
                kwargs[key] = None
                result = REPORT.build_report(**kwargs, allow_partial=True)
                self.assertFalse(result.final_ready)
                self.assertTrue(any(expected_issue in issue for issue in result.final_issues))
                self.assertIn("partial-preview", result.output_path.name)

    def test_partial_with_real_funnel_verifier_has_only_four_m2_issues(self) -> None:
        paths, _ = self.materialize()
        kwargs = self.kwargs(paths, self.root / REPORT.FINAL_REPORT_NAME)
        kwargs["m2_extraction_receipt"] = None
        kwargs["m2_ranking_receipt"] = None
        result = REPORT.build_report(**kwargs, allow_partial=True)
        self.assertEqual(
            result.final_issues,
            (
                "缺少 M2 full extraction receipt",
                "缺少 M2 full ranking receipt",
                "M2 independent verifier 無法套用：full extraction/ranking/candidate sources 尚未全部通過",
                "M2 final numeric summary 無法套用：terminal extraction/ranking/final verifier 尚未全部通過",
            ),
        )
        self.assertNotIn("funnel", " ".join(result.final_issues).lower())

    def test_funnel_verifier_contract_scope_and_bindings_fail_closed(self) -> None:
        mutations = (
            ("schema", "schema_version", lambda payload, paths: payload.__setitem__("schema_version", "0.9.0")),
            ("all_pass", "all_pass", lambda payload, paths: payload.__setitem__("all_pass", False)),
            ("322_checks", "checks 不是 322", lambda payload, paths: payload["checks"].pop("check_322")),
            ("zero_fail", "n_fail 不是 0", lambda payload, paths: payload.__setitem__("n_fail", 1)),
            (
                "task_scope",
                "dataset_count 不是 7",
                lambda payload, paths: payload["scope"].__setitem__("dataset_count", 6),
            ),
            (
                "canonical_path",
                "canonical path binding",
                lambda payload, paths: payload["inputs"]["canonical_json"].__setitem__(
                    "path", str(paths["funnel"].resolve())
                ),
            ),
            (
                "canonical_sha",
                "canonical SHA-256 binding",
                lambda payload, paths: payload["inputs"]["canonical_json"].__setitem__("sha256", "0" * 64),
            ),
            (
                "producer_path",
                "producer funnel path binding",
                lambda payload, paths: payload["inputs"]["receipt_under_test"].__setitem__(
                    "path", str(paths["canonical"].resolve())
                ),
            ),
            (
                "producer_sha",
                "producer funnel SHA-256 binding",
                lambda payload, paths: payload["inputs"]["receipt_under_test"].__setitem__("sha256", "0" * 64),
            ),
        )
        for label, expected_error, mutate in mutations:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                verification = json.loads(paths["funnel_verification"].read_text(encoding="utf-8"))
                mutate(verification, paths)
                paths["funnel_verification"].write_text(
                    json.dumps(verification, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_funnel_verifier_{label}.html"),
                        allow_partial=True,
                    )

    def test_per_dataset_canonical_invariants_catch_aggregate_cancellation(self) -> None:
        topology_template = {
            "exact_and_topology_unique": 2,
            "topology_unique_exact_multiple": 2,
            "topology_multiple_exact_multiple": 4,
            "incomplete": 1,
            "impossible_exact_unique_topology_multiple": 0,
        }
        mutations = (
            (
                "complete",
                r"complete \+ incomplete 不等於 W_primary",
                lambda first, second: (
                    first.__setitem__("complete_regions", first["complete_regions"] + 1),
                    second.__setitem__("complete_regions", second["complete_regions"] - 1),
                ),
            ),
            (
                "no_primary",
                r"W_primary \+ no_primary_lineage 不等於 W_tree",
                lambda first, second: (
                    first.__setitem__("no_primary_lineage", first["no_primary_lineage"] + 1),
                    second.__setitem__("no_primary_lineage", second["no_primary_lineage"] - 1),
                ),
            ),
            (
                "topology",
                "topology 互斥分類總和不等於 W_primary",
                lambda first, second: (
                    first["topology_classes"].__setitem__(
                        "exact_and_topology_unique",
                        first["topology_classes"]["exact_and_topology_unique"] + 1,
                    ),
                    second["topology_classes"].__setitem__(
                        "exact_and_topology_unique",
                        second["topology_classes"]["exact_and_topology_unique"] - 1,
                    ),
                ),
            ),
        )
        for label, expected_error, mutate in mutations:
            with self.subTest(label=label):
                bundle = copy.deepcopy(self.bundle)
                sample_rows = bundle["canonical"]["canonical"]["samples"]
                for row in sample_rows:
                    row["no_primary_lineage"] = 1
                    row["primary_units"] = 11
                    row["topology_classes"] = copy.deepcopy(topology_template)
                mutate(sample_rows[0], sample_rows[1])
                self.assertEqual(
                    sum(row["complete_regions"] for row in sample_rows),
                    bundle["canonical"]["canonical"]["aggregate"]["complete_regions"],
                )
                self.assertEqual(
                    sum(row["no_primary_lineage"] for row in sample_rows),
                    bundle["canonical"]["canonical"]["aggregate"]["no_primary_lineage"],
                )
                self.assertEqual(
                    sum(row["topology_classes"]["exact_and_topology_unique"] for row in sample_rows),
                    bundle["canonical"]["canonical"]["aggregate"]["topology_classes"][
                        "exact_and_topology_unique"
                    ],
                )
                paths, _ = self.materialize(bundle)
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_canonical_{label}.html"),
                        allow_partial=True,
                    )

    def test_m0_rows_post_receipt_mutation_fails_closed(self) -> None:
        paths, _ = self.materialize()
        with (self.root / "m0_rows.tsv.gz").open("ab") as handle:
            handle.write(b"tampered")
        final = self.root / REPORT.FINAL_REPORT_NAME
        with self.assertRaisesRegex(REPORT.ReportGateError, "M0 rows TSV (?:size_bytes|SHA-256)"):
            REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)

    def test_legacy_s8_raw_hp_denominator_includes_missing_category(self) -> None:
        mutations = (
            ("missing_category", lambda counts: counts["raw_HP_counts"].pop("."), "缺 missing HP='.'"),
            (
                "denominator_mismatch",
                lambda counts: counts.__setitem__(
                    "canonical_eligible_alignments", counts["canonical_eligible_alignments"] + 1
                ),
                "總和不等於 canonical eligible alignments",
            ),
        )
        for label, mutate, expected_error in mutations:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                pilot = json.loads(paths["m2_pilot"].read_text(encoding="utf-8"))
                mutate(pilot["counts"])
                paths["m2_pilot"].write_text(
                    json.dumps(pilot, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(
                        **self.kwargs(paths, self.root / f"bad_s8_{label}.html"),
                        allow_partial=True,
                    )

    def test_extraction_rejects_pass_154_plus_extra_skip(self) -> None:
        paths, _ = self.materialize()
        extraction = json.loads(paths["m2_extraction"].read_text(encoding="utf-8"))
        extraction["task_status_counts"]["SKIP"] = 1
        extraction["n_results"] = 155
        extraction["results"].append({"dataset": "COLO829", "chrom": "chr1", "status": "SKIP", "receipt": {}})
        paths["m2_extraction"].write_text(json.dumps(extraction, ensure_ascii=False) + "\n", encoding="utf-8")
        ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
        ranking["upstream_extraction_receipt"]["sha256"] = sha256(paths["m2_extraction"])
        paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8")
        final = self.root / REPORT.FINAL_REPORT_NAME
        with self.assertRaisesRegex(REPORT.ReportGateError, "PASS/REUSED_PASS|恰好 154"):
            REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)

    def test_ranking_requires_both_ps_strata_and_per_dataset_conservation(self) -> None:
        mutations = (
            (
                "missing_ps_hp2",
                lambda ranking: ranking["aggregate"]["by_linkage_basis_and_threshold"].pop("PS_HP2"),
                "缺 PS_HP1/PS_HP2",
            ),
            (
                "dataset_sum",
                lambda ranking: ranking["by_dataset"]["COLO829"]["by_linkage_basis_and_threshold"]["PS_HP1"]["3"].__setitem__(
                    "raw_tree_candidates_T_complete_units",
                    ranking["by_dataset"]["COLO829"]["by_linkage_basis_and_threshold"]["PS_HP1"]["3"]["raw_tree_candidates_T_complete_units"] + 1,
                ),
                "per-dataset PS_HP1/3 加總不符",
            ),
        )
        for label, mutate, expected in mutations:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
                mutate(ranking)
                paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8")
                final = self.root / REPORT.FINAL_REPORT_NAME
                with self.assertRaisesRegex(REPORT.ReportGateError, expected):
                    REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)

    def test_funnel_branch_and_m2_verifier_hashes_fail_closed(self) -> None:
        paths, _ = self.materialize()
        funnel = json.loads(paths["funnel"].read_text(encoding="utf-8"))
        funnel["aggregate"]["branch_counts"]["positional_singleton"] += 1
        paths["funnel"].write_text(json.dumps(funnel, ensure_ascii=False) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(REPORT.ReportGateError, "branches 不守恆"):
            REPORT.build_report(**self.kwargs(paths, self.root / "bad_funnel.html"), allow_partial=True)

        paths, _ = self.materialize()
        verification = json.loads(paths["m2_verification"].read_text(encoding="utf-8"))
        verification["ranking"]["candidate_table"]["sha256"] = "0" * 64
        paths["m2_verification"].write_text(json.dumps(verification, ensure_ascii=False) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(REPORT.ReportGateError, "candidate table SHA-256"):
            REPORT.build_report(**self.kwargs(paths, self.root / "bad_m2_verify.html"), allow_partial=True)

    def test_topology_and_partial_denominator_nonconservation_fail_closed(self) -> None:
        mutations = (
            ("topology", lambda cell: cell["coarse_topology_unique_class_counts"].__setitem__("single-only", cell["coarse_topology_unique_class_counts"]["single-only"] + 1), "coarse unique classes"),
            ("partial", lambda cell: cell["partial_u_distribution"]["molecule_projections"].__setitem__("2", cell["partial_u_distribution"]["molecule_projections"]["2"] + 1), "partial u distribution molecule_projections"),
        )
        for label, mutate, expected_error in mutations:
            with self.subTest(label=label):
                paths, _ = self.materialize()
                ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
                cell = ranking["aggregate"]["by_linkage_basis_and_threshold"]["PS_HP1"]["3"]
                mutate(cell)
                paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8")
                final = self.root / REPORT.FINAL_REPORT_NAME
                with self.assertRaisesRegex(REPORT.ReportGateError, expected_error):
                    REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)
                self.assertFalse(final.exists())

    def test_legacy_m2_schemas_are_diagnostic_only_and_cannot_create_final(self) -> None:
        paths, _ = self.materialize()
        extraction = json.loads(paths["m2_extraction"].read_text(encoding="utf-8"))
        extraction["schema_version"] = "1.1.0"
        extraction["scope"].pop("child_schema_version", None)
        extraction["checks"].pop("all_child_receipts_schema_1_2", None)
        paths["m2_extraction"].write_text(json.dumps(extraction, ensure_ascii=False) + "\n", encoding="utf-8")
        ranking = json.loads(paths["m2_ranking"].read_text(encoding="utf-8"))
        ranking["schema_version"] = "1.0.0"
        ranking["upstream_extraction_receipt"]["sha256"] = sha256(paths["m2_extraction"])
        paths["m2_ranking"].write_text(json.dumps(ranking, ensure_ascii=False) + "\n", encoding="utf-8")
        final = self.root / REPORT.FINAL_REPORT_NAME
        result = REPORT.build_report(**self.kwargs(paths, final), allow_partial=True)
        self.assertFalse(result.final_ready)
        self.assertFalse(final.exists())
        self.assertIn("partial-preview", result.output_path.name)
        document = result.output_path.read_text(encoding="utf-8")
        self.assertIn("schema 1.1.0；尚未證明", document)
        self.assertIn("舊 ranking schema 僅作 diagnostic", document)
        self.assertNotIn("FINAL · TASK-B EVIDENCE GATE PASS", document)


if __name__ == "__main__":
    unittest.main()
