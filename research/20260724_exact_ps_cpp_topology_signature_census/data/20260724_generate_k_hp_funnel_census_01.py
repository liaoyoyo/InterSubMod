#!/usr/bin/env python3
"""Recompute the seven-dataset exact-PS k/HP/funnel census from immutable inputs."""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
STRICT_REPORT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260723_production_exact_ps_strict_read_linkage/all7_production_v1/"
    "all7_report_data_v1/all7_report_data.json"
)
INPUT_MANIFEST = (
    REPO_ROOT
    / "research/20260724_exact_ps_cpp_topology_af_all_samples/input_contract/"
    "all7_exact_ps_inputs.local.json"
)
TOPOLOGY_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples"
)
OUTPUT_DIR = (
    REPO_ROOT / "research/20260724_exact_ps_cpp_topology_signature_census/data"
)

OVERVIEW_TSV = OUTPUT_DIR / "20260724_exactPS_source_overview_by_sample_01.tsv"
SOURCE_K_TSV = OUTPUT_DIR / "20260724_exactPS_source_component_k_by_sample_01.tsv"
FINAL_K_TSV = OUTPUT_DIR / "20260724_exactPS_final_group_k_by_sample_01.tsv"
SOURCE_K_MEMBERSHIP_TSV = (
    OUTPUT_DIR / "20260724_exactPS_source_component_k_memberships_by_sample_01.tsv"
)
FINAL_K_MEMBERSHIP_TSV = (
    OUTPUT_DIR / "20260724_exactPS_final_group_k_memberships_by_sample_01.tsv"
)
HP_TSV = OUTPUT_DIR / "20260724_exactPS_hp_split_by_sample_01.tsv"
FUNNEL_TSV = OUTPUT_DIR / "20260724_exactPS_funnel_by_sample_01.tsv"
OUTPUT_JSON = OUTPUT_DIR / "20260724_exactPS_k_hp_funnel_census_01.json"


def pct(count: int, denominator: int) -> float:
    return round(100.0 * count / denominator, 6) if denominator else 0.0


def hp_from_linkage_basis(linkage_basis: str) -> str:
    if linkage_basis == "PS_HP1":
        return "HP1"
    if linkage_basis == "PS_HP2":
        return "HP2"
    raise ValueError(f"Unexpected linkage_basis: {linkage_basis}")


def dump_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def counter_to_bands(
    counter: Counter[int],
    region_denominator: int,
    membership_denominator: int,
) -> dict[str, dict]:
    result: dict[str, dict] = {}
    for k in range(1, 13):
        count = int(counter[k])
        membership_count = k * count
        result[str(k)] = {
            "count": count,
            "denominator": region_denominator,
            "percent": pct(count, region_denominator),
            "membership_count": membership_count,
            "membership_denominator": membership_denominator,
            "membership_percent": pct(membership_count, membership_denominator),
        }
    gt12 = sum(count for k, count in counter.items() if k > 12)
    gt12_memberships = sum(k * count for k, count in counter.items() if k > 12)
    result[">12"] = {
        "count": int(gt12),
        "denominator": region_denominator,
        "percent": pct(int(gt12), region_denominator),
        "membership_count": int(gt12_memberships),
        "membership_denominator": membership_denominator,
        "membership_percent": pct(int(gt12_memberships), membership_denominator),
    }
    return result


def sum_dicts(rows: list[dict], keys: list[str]) -> dict[str, int]:
    return {key: sum(int(row[key]) for row in rows) for key in keys}


def main() -> None:
    strict_report = json.loads(STRICT_REPORT.read_text(encoding="utf-8"))
    manifest = json.loads(INPUT_MANIFEST.read_text(encoding="utf-8"))
    dataset_order = manifest["dataset_order"]
    report_by_dataset = {row["dataset"]: row for row in strict_report["datasets"]}

    records: dict[str, dict] = {}
    all_checks: list[dict] = []
    total_source_k: Counter[int] = Counter()
    total_final_k: Counter[int] = Counter()

    for dataset in dataset_order:
        report_row = report_by_dataset[dataset]
        strict_root = Path(manifest["samples"][dataset]["strict_root"])
        component_files = sorted(strict_root.glob("chromosomes/chr*/*.components.tsv.gz"))
        container_files = sorted(
            strict_root.glob("chromosomes/chr*/*.container_summary.tsv.gz")
        )
        if len(component_files) != 22 or len(container_files) != 22:
            raise RuntimeError(
                f"{dataset}: expected 22 component/container files, found "
                f"{len(component_files)}/{len(container_files)}"
            )

        source_hp = {
            hp: {
                "containers": 0,
                "components": 0,
                "memberships": 0,
                "singleton_k1_components": 0,
                "tree_eligible_W": 0,
                "W_memberships": 0,
                "k_counter": Counter(),
            }
            for hp in ("HP1", "HP2")
        }

        for component_file in component_files:
            with gzip.open(component_file, "rt", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if int(row["threshold"]) != 3:
                        continue
                    hp = hp_from_linkage_basis(row["linkage_basis"])
                    k = int(row["k"])
                    source_hp[hp]["components"] += 1
                    source_hp[hp]["memberships"] += k
                    source_hp[hp]["k_counter"][k] += 1
                    if k == 1:
                        source_hp[hp]["singleton_k1_components"] += 1
                    else:
                        source_hp[hp]["tree_eligible_W"] += 1
                        source_hp[hp]["W_memberships"] += k

        for container_file in container_files:
            with gzip.open(container_file, "rt", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if int(row["threshold"]) != 3:
                        continue
                    hp = f"HP{row['hp_family']}"
                    if hp not in source_hp:
                        raise ValueError(f"{dataset}: unexpected hp_family={row['hp_family']}")
                    source_hp[hp]["containers"] += 1

        mlhp_path = TOPOLOGY_ROOT / dataset / f"{dataset}.exact_ps_mlhp.json"
        mlhp = json.loads(mlhp_path.read_text(encoding="utf-8"))
        final_hp = {
            hp: {"groups": 0, "memberships": 0, "k_counter": Counter()}
            for hp in ("HP1", "HP2")
        }
        for group in mlhp["groups"]:
            hp = f"HP{group['hp_family']}"
            k = int(group["n_sSNV"])
            final_hp[hp]["groups"] += 1
            final_hp[hp]["memberships"] += k
            final_hp[hp]["k_counter"][k] += 1

        source_k = source_hp["HP1"]["k_counter"] + source_hp["HP2"]["k_counter"]
        final_k = final_hp["HP1"]["k_counter"] + final_hp["HP2"]["k_counter"]
        total_source_k.update(source_k)
        total_final_k.update(final_k)

        source_totals = {
            key: sum(int(source_hp[hp][key]) for hp in ("HP1", "HP2"))
            for key in (
                "containers",
                "components",
                "memberships",
                "singleton_k1_components",
                "tree_eligible_W",
                "W_memberships",
            )
        }
        final_totals = {
            key: sum(int(final_hp[hp][key]) for hp in ("HP1", "HP2"))
            for key in ("groups", "memberships")
        }
        funnel = mlhp["input_funnel"]

        checks = {
            "22_autosomes_present": len(component_files) == 22
            and len(container_files) == 22,
            "containers_match_strict_report": (
                source_totals["containers"] == report_row["exact_PS_HP_containers"]
            ),
            "components_match_strict_report": (
                source_totals["components"] == report_row["all_components"]
            ),
            "active_memberships_match_strict_report": (
                source_totals["memberships"] == report_row["active_node_memberships"]
            ),
            "singleton_match_strict_report": (
                source_totals["singleton_k1_components"]
                == report_row["singleton_k1_components"]
            ),
            "W_match_strict_report": (
                source_totals["tree_eligible_W"] == report_row["tree_eligible_W"]
            ),
            "W_memberships_match_strict_report": (
                source_totals["W_memberships"] == report_row["W_memberships"]
            ),
            "component_partition_conserved": (
                source_totals["components"]
                == source_totals["singleton_k1_components"]
                + source_totals["tree_eligible_W"]
            ),
            "membership_partition_conserved": (
                source_totals["memberships"]
                == source_totals["singleton_k1_components"]
                + source_totals["W_memberships"]
            ),
            "final_group_count_matches_mlhp": (
                final_totals["groups"] == mlhp["n_groups_analyzed"]
            ),
            "source_k_band_memberships_conserved": (
                sum(k * count for k, count in source_k.items())
                == source_totals["memberships"]
            ),
            "final_k_band_memberships_conserved": (
                sum(k * count for k, count in final_k.items())
                == final_totals["memberships"]
            ),
            "final_original_k_in_2_12": (
                final_k[1] == 0
                and sum(count for k, count in final_k.items() if k < 2 or k > 12) == 0
            ),
            "bounded_block_funnel_conserved": (
                funnel["bounded_blocks"]
                - funnel["k1_blocks_not_tree_eligible"]
                - funnel["pattern_unsupported_blocks"]
                == funnel["tree_input_groups"]
                == final_totals["groups"]
            ),
            "membership_funnel_conserved": (
                funnel["unit_memberships"]
                - funnel["k1_blocks_not_tree_eligible"]
                - funnel["pattern_unsupported_memberships"]
                == funnel["tree_input_memberships"]
                == final_totals["memberships"]
            ),
        }
        checks["all_pass"] = all(checks.values())
        all_checks.append({"dataset": dataset, **checks})
        if not checks["all_pass"]:
            raise RuntimeError(f"{dataset}: reconciliation failure: {checks}")

        source_overview = {
            "candidate_loci_S": int(report_row["candidate_loci_S"]),
            "active_unique_loci": int(report_row["active_unique_loci"]),
            "active_unique_pct_of_S": pct(
                int(report_row["active_unique_loci"]),
                int(report_row["candidate_loci_S"]),
            ),
            "active_node_memberships": source_totals["memberships"],
            "exact_PS_HP_containers": source_totals["containers"],
            "all_components": source_totals["components"],
            "singleton_k1_components": source_totals["singleton_k1_components"],
            "singleton_pct_of_all_components": pct(
                source_totals["singleton_k1_components"],
                source_totals["components"],
            ),
            "tree_eligible_W": source_totals["tree_eligible_W"],
            "W_pct_of_all_components": pct(
                source_totals["tree_eligible_W"], source_totals["components"]
            ),
            "W_memberships": source_totals["W_memberships"],
            "W_membership_pct_of_active_memberships": pct(
                source_totals["W_memberships"], source_totals["memberships"]
            ),
        }

        hp_output: dict[str, dict] = {}
        for hp in ("HP1", "HP2"):
            hp_output[hp] = {}
            for key in (
                "containers",
                "components",
                "memberships",
                "singleton_k1_components",
                "tree_eligible_W",
                "W_memberships",
            ):
                count = int(source_hp[hp][key])
                hp_output[hp][key] = {
                    "count": count,
                    "denominator": source_totals[key],
                    "percent": pct(count, source_totals[key]),
                }
            for key in ("groups", "memberships"):
                count = int(final_hp[hp][key])
                output_key = f"final_{key}"
                hp_output[hp][output_key] = {
                    "count": count,
                    "denominator": final_totals[key],
                    "percent": pct(count, final_totals[key]),
                }

        records[dataset] = {
            "inputs": {
                "strict_root": str(strict_root),
                "mlhp": str(mlhp_path),
            },
            "source_overview": source_overview,
            "source_component_k": counter_to_bands(
                source_k,
                source_totals["components"],
                source_totals["memberships"],
            ),
            "final_group_k": counter_to_bands(
                final_k,
                final_totals["groups"],
                final_totals["memberships"],
            ),
            "hp_split": hp_output,
            "funnel": {
                "source_W": int(funnel["exact_ps_hp_units"]),
                "source_W_memberships": int(funnel["unit_memberships"]),
                "bounded_blocks": int(funnel["bounded_blocks"]),
                "post_cut_k1_blocks_excluded": int(
                    funnel["k1_blocks_not_tree_eligible"]
                ),
                "pattern_unsupported_blocks_excluded": int(
                    funnel["pattern_unsupported_blocks"]
                ),
                "pattern_unsupported_memberships_excluded": int(
                    funnel["pattern_unsupported_memberships"]
                ),
                "final_groups": int(funnel["tree_input_groups"]),
                "final_memberships": int(funnel["tree_input_memberships"]),
                "final_membership_retention_pct": pct(
                    int(funnel["tree_input_memberships"]),
                    int(funnel["unit_memberships"]),
                ),
            },
            "checks": checks,
        }

    overview_rows = [
        {"dataset": dataset, **records[dataset]["source_overview"]}
        for dataset in dataset_order
    ]
    overview_total_keys = [
        "candidate_loci_S",
        "active_unique_loci",
        "active_node_memberships",
        "exact_PS_HP_containers",
        "all_components",
        "singleton_k1_components",
        "tree_eligible_W",
        "W_memberships",
    ]
    overview_total = sum_dicts(overview_rows, overview_total_keys)
    overview_total.update(
        {
            "dataset": "TOTAL",
            "active_unique_pct_of_S": pct(
                overview_total["active_unique_loci"],
                overview_total["candidate_loci_S"],
            ),
            "singleton_pct_of_all_components": pct(
                overview_total["singleton_k1_components"],
                overview_total["all_components"],
            ),
            "W_pct_of_all_components": pct(
                overview_total["tree_eligible_W"],
                overview_total["all_components"],
            ),
            "W_membership_pct_of_active_memberships": pct(
                overview_total["W_memberships"],
                overview_total["active_node_memberships"],
            ),
        }
    )
    overview_rows.append(overview_total)

    source_total_denominator = overview_total["all_components"]
    final_total_denominator = sum(
        records[dataset]["funnel"]["final_groups"] for dataset in dataset_order
    )
    final_membership_denominator = sum(
        records[dataset]["funnel"]["final_memberships"] for dataset in dataset_order
    )
    source_k_rows: list[dict] = []
    final_k_rows: list[dict] = []
    source_k_membership_rows: list[dict] = []
    final_k_membership_rows: list[dict] = []
    for dataset in dataset_order + ["TOTAL"]:
        if dataset == "TOTAL":
            source_bands = counter_to_bands(
                total_source_k,
                source_total_denominator,
                overview_total["active_node_memberships"],
            )
            final_bands = counter_to_bands(
                total_final_k,
                final_total_denominator,
                final_membership_denominator,
            )
        else:
            source_bands = records[dataset]["source_component_k"]
            final_bands = records[dataset]["final_group_k"]
        for band in [str(k) for k in range(1, 13)] + [">12"]:
            source_band = source_bands[band]
            final_band = final_bands[band]
            source_k_rows.append(
                {
                    "dataset": dataset,
                    "k_band": band,
                    "count": source_band["count"],
                    "denominator": source_band["denominator"],
                    "percent": source_band["percent"],
                }
            )
            final_k_rows.append(
                {
                    "dataset": dataset,
                    "k_band": band,
                    "count": final_band["count"],
                    "denominator": final_band["denominator"],
                    "percent": final_band["percent"],
                }
            )
            source_k_membership_rows.append(
                {
                    "dataset": dataset,
                    "k_band": band,
                    "membership_count": source_band["membership_count"],
                    "membership_denominator": source_band[
                        "membership_denominator"
                    ],
                    "membership_percent": source_band["membership_percent"],
                }
            )
            final_k_membership_rows.append(
                {
                    "dataset": dataset,
                    "k_band": band,
                    "membership_count": final_band["membership_count"],
                    "membership_denominator": final_band[
                        "membership_denominator"
                    ],
                    "membership_percent": final_band["membership_percent"],
                }
            )

    hp_rows: list[dict] = []
    hp_metric_names = [
        "containers",
        "components",
        "memberships",
        "singleton_k1_components",
        "tree_eligible_W",
        "W_memberships",
        "final_groups",
        "final_memberships",
    ]
    for dataset in dataset_order:
        for hp in ("HP1", "HP2"):
            row = {"dataset": dataset, "hp_family": hp}
            for metric in hp_metric_names:
                value = records[dataset]["hp_split"][hp][metric]
                row[f"{metric}_count"] = value["count"]
                row[f"{metric}_denominator"] = value["denominator"]
                row[f"{metric}_percent"] = value["percent"]
            hp_rows.append(row)
    for hp in ("HP1", "HP2"):
        total_row = {"dataset": "TOTAL", "hp_family": hp}
        for metric in hp_metric_names:
            count = sum(
                records[dataset]["hp_split"][hp][metric]["count"]
                for dataset in dataset_order
            )
            denominator = sum(
                records[dataset]["hp_split"][hp][metric]["denominator"]
                for dataset in dataset_order
            )
            total_row[f"{metric}_count"] = count
            total_row[f"{metric}_denominator"] = denominator
            total_row[f"{metric}_percent"] = pct(count, denominator)
        hp_rows.append(total_row)

    funnel_rows = [
        {"dataset": dataset, **records[dataset]["funnel"]}
        for dataset in dataset_order
    ]
    funnel_total_keys = [
        "source_W",
        "source_W_memberships",
        "bounded_blocks",
        "post_cut_k1_blocks_excluded",
        "pattern_unsupported_blocks_excluded",
        "pattern_unsupported_memberships_excluded",
        "final_groups",
        "final_memberships",
    ]
    funnel_total = sum_dicts(funnel_rows, funnel_total_keys)
    funnel_total.update(
        {
            "dataset": "TOTAL",
            "final_membership_retention_pct": pct(
                funnel_total["final_memberships"],
                funnel_total["source_W_memberships"],
            ),
        }
    )
    funnel_rows.append(funnel_total)

    total_hp_split = {
        hp: {
            metric: {
                "count": next(
                    row[f"{metric}_count"]
                    for row in hp_rows
                    if row["dataset"] == "TOTAL" and row["hp_family"] == hp
                ),
                "denominator": next(
                    row[f"{metric}_denominator"]
                    for row in hp_rows
                    if row["dataset"] == "TOTAL" and row["hp_family"] == hp
                ),
                "percent": next(
                    row[f"{metric}_percent"]
                    for row in hp_rows
                    if row["dataset"] == "TOTAL" and row["hp_family"] == hp
                ),
            }
            for metric in hp_metric_names
        }
        for hp in ("HP1", "HP2")
    }

    cohort_checks = {
        "all_sample_checks_pass": all(row["all_pass"] for row in all_checks),
        "source_component_total_is_255752": (
            overview_total["all_components"] == 255_752
        ),
        "source_partition_is_170131_plus_85621": (
            overview_total["all_components"]
            == overview_total["singleton_k1_components"]
            + overview_total["tree_eligible_W"]
            == 170_131 + 85_621
        ),
        "source_membership_partition_is_170131_plus_443349": (
            overview_total["active_node_memberships"]
            == overview_total["singleton_k1_components"]
            + overview_total["W_memberships"]
            == 170_131 + 443_349
        ),
        "final_group_total_is_98955": final_total_denominator == 98_955,
        "final_k1_is_zero": total_final_k[1] == 0,
        "source_k_band_memberships_sum_to_613480": (
            sum(k * count for k, count in total_source_k.items()) == 613_480
        ),
        "final_k_band_memberships_sum_to_439685": (
            sum(k * count for k, count in total_final_k.items()) == 439_685
        ),
        "all_final_groups_have_original_k_2_12": (
            sum(total_final_k[k] for k in range(2, 13)) == 98_955
            and sum(count for k, count in total_final_k.items() if k > 12) == 0
        ),
        "cohort_block_funnel_conserved": (
            funnel_total["bounded_blocks"]
            - funnel_total["post_cut_k1_blocks_excluded"]
            - funnel_total["pattern_unsupported_blocks_excluded"]
            == funnel_total["final_groups"]
            == 98_955
        ),
        "cohort_membership_funnel_conserved": (
            funnel_total["source_W_memberships"]
            - funnel_total["post_cut_k1_blocks_excluded"]
            - funnel_total["pattern_unsupported_memberships_excluded"]
            == funnel_total["final_memberships"]
            == 439_685
        ),
    }
    cohort_checks["all_pass"] = all(cohort_checks.values())
    if not cohort_checks["all_pass"]:
        raise RuntimeError(f"Cohort reconciliation failure: {cohort_checks}")

    output = {
        "schema_name": "intersubmod.exact_ps_k_hp_funnel_census",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "task_type": "comprehensive_validation",
            "datasets": dataset_order,
            "chromosomes": "chr1-22",
            "strict_endpoint_support_threshold": 3,
            "source_analysis_unit": "exact_PS_x_primary_HP_x_read_linked_component",
            "final_analysis_unit": "exact_PS_x_primary_HP_x_component_x_bounded_block",
            "final_k_definition": "original_bit_count/n_sSNV before active-column reduction",
            "mutation_model": "not evaluated in this census",
        },
        "definitions": {
            "HP1": "merged raw family 1/1-1/1-2",
            "HP2": "merged raw family 2/2-1/2-2",
            "source_W": "threshold-3 retained-endpoint connected component with k>=2",
            "singleton": "source component with k=1; excluded before tree partition",
            "percent_denominators": {
                "source_component_k": "all source components within dataset",
                "source_component_k_memberships": (
                    "all source component sSNV memberships within dataset"
                ),
                "final_group_k": "all final topology groups within dataset",
                "final_group_k_memberships": (
                    "all final group sSNV memberships within dataset"
                ),
                "hp_split": "HP1+HP2 total for the same metric within dataset",
            },
        },
        "inputs": {
            "strict_report": str(STRICT_REPORT),
            "input_manifest": str(INPUT_MANIFEST),
            "topology_root": str(TOPOLOGY_ROOT),
        },
        "samples": records,
        "totals": {
            "source_overview": overview_total,
            "source_component_k": counter_to_bands(
                total_source_k,
                source_total_denominator,
                overview_total["active_node_memberships"],
            ),
            "final_group_k": counter_to_bands(
                total_final_k,
                final_total_denominator,
                funnel_total["final_memberships"],
            ),
            "hp_split": total_hp_split,
            "funnel": funnel_total,
        },
        "checks": {
            "per_sample": all_checks,
            "cohort": cohort_checks,
        },
    }

    dump_tsv(OVERVIEW_TSV, list(overview_rows[0].keys()), overview_rows)
    dump_tsv(
        SOURCE_K_TSV,
        ["dataset", "k_band", "count", "denominator", "percent"],
        source_k_rows,
    )
    dump_tsv(
        FINAL_K_TSV,
        ["dataset", "k_band", "count", "denominator", "percent"],
        final_k_rows,
    )
    dump_tsv(
        SOURCE_K_MEMBERSHIP_TSV,
        [
            "dataset",
            "k_band",
            "membership_count",
            "membership_denominator",
            "membership_percent",
        ],
        source_k_membership_rows,
    )
    dump_tsv(
        FINAL_K_MEMBERSHIP_TSV,
        [
            "dataset",
            "k_band",
            "membership_count",
            "membership_denominator",
            "membership_percent",
        ],
        final_k_membership_rows,
    )
    dump_tsv(HP_TSV, list(hp_rows[0].keys()), hp_rows)
    dump_tsv(FUNNEL_TSV, list(funnel_rows[0].keys()), funnel_rows)
    OUTPUT_JSON.write_text(
        json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "outputs": [
                    str(OVERVIEW_TSV),
                    str(SOURCE_K_TSV),
                    str(FINAL_K_TSV),
                    str(SOURCE_K_MEMBERSHIP_TSV),
                    str(FINAL_K_MEMBERSHIP_TSV),
                    str(HP_TSV),
                    str(FUNNEL_TSV),
                    str(OUTPUT_JSON),
                ],
                "cohort_checks": cohort_checks,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
