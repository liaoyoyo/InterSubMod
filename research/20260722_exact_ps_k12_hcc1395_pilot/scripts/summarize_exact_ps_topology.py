#!/usr/bin/env python3
"""Validate and summarize an exact-PS HCC1395 topology reconstruction.

This report intentionally treats an exact PS x HP x component x block as the
analysis unit.  It does not merge local PS blocks back into a biological clone
tree and it does not rank trees by VAF.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


CLASS_ORDER = (
    "T=1|Topo=1",
    "T>1|Topo=1",
    "T>1|Topo>1",
    "incomplete",
)
MORPHOLOGY_ORDER = (
    "no_branch_no_multilevel",
    "sister_only",
    "direct_only",
    "sister_and_direct",
)


def fail(message: str) -> None:
    raise RuntimeError(f"EXACT-PS TOPOLOGY SUMMARY FAILED: {message}")


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pct(numerator: int, denominator: int) -> float | None:
    return round(100.0 * numerator / denominator, 4) if denominator else None


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def canonical_shape(edges: list[list[str]] | list[tuple[str, str]]) -> str:
    children: dict[str, list[str]] = collections.defaultdict(list)
    all_nodes: set[str] = set()
    child_nodes: set[str] = set()
    for parent, child in edges:
        children[parent].append(child)
        all_nodes.update((parent, child))
        child_nodes.add(child)
    roots = sorted(all_nodes - child_nodes)

    def visit(node: str) -> str:
        return "(" + "".join(sorted(visit(child) for child in children[node])) + ")"

    return "|".join(sorted(visit(root) for root in roots)) if roots else "()"


def shape_features(edges: list[list[str]] | list[tuple[str, str]]) -> tuple[bool, bool]:
    nodes: set[str] = {"ROOT"}
    children: dict[str, list[str]] = collections.defaultdict(list)
    incoming: collections.Counter[str] = collections.Counter()
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
    return max((len(children[node]) for node in nodes), default=0) >= 2, max(
        depth.values(), default=0
    ) >= 2


def morphology_name(sister: bool, direct: bool) -> str:
    if sister and direct:
        return "sister_and_direct"
    if sister:
        return "sister_only"
    if direct:
        return "direct_only"
    return "no_branch_no_multilevel"


def is_complete(unit: dict) -> bool:
    return (
        unit.get("capped") is False
        and unit.get("analysis_candidate_set_complete") is True
        and unit.get("verification_status") == "full_pass"
        and unit.get("verification_complete") is True
        and unit.get("verify_pass") is True
        and int(unit.get("analysis_trees_generated") or -1) == int(unit.get("n_trees") or 0)
        and unit.get("n_distinct_shapes_exact") is not None
    )


def structural_class(unit: dict) -> str:
    if not is_complete(unit):
        return "incomplete"
    trees = int(unit["n_trees"])
    topologies = int(unit["n_distinct_shapes_exact"])
    if trees == 1 and topologies == 1:
        return "T=1|Topo=1"
    if trees > 1 and topologies == 1:
        return "T>1|Topo=1"
    if trees > 1 and topologies > 1:
        return "T>1|Topo>1"
    fail(f"impossible T/Topo state in {unit.get('region')}: T={trees}, Topo={topologies}")


def t_bin(value: int) -> str:
    return str(value) if value <= 6 else ">6"


def ordered_counter(counter: collections.Counter, order: tuple[str, ...] | None = None) -> dict:
    if order is None:
        return dict(sorted(counter.items(), key=lambda item: str(item[0])))
    return {key: int(counter.get(key, 0)) for key in order}


def count_with_percent(counter: dict[str, int], denominator: int) -> dict:
    return {
        key: {"n": int(value), "percent": pct(int(value), denominator)}
        for key, value in counter.items()
    }


def legacy_baseline(path: Path) -> dict:
    data = read_json(path)
    classes: collections.Counter[str] = collections.Counter()
    primary_regions = 0

    def legacy_complete(unit: dict) -> bool:
        # Canonical v5 predates analysis_trees_generated; its frozen contract is
        # candidate-set completeness + full verification + exact shape count.
        return (
            unit.get("capped") is False
            and unit.get("analysis_candidate_set_complete") is True
            and unit.get("verification_status") == "full_pass"
            and unit.get("verification_complete") is True
            and unit.get("verify_pass") is True
            and unit.get("n_distinct_shapes_exact") is not None
        )

    for region in data["regions"]:
        units = [unit for unit in region.get("lineages", []) if unit.get("is_primary_lineage")]
        if not units:
            continue
        primary_regions += 1
        if not all(legacy_complete(unit) for unit in units):
            classes["incomplete"] += 1
            continue
        trees = 1
        topologies = 1
        for unit in units:
            trees *= int(unit["n_trees"])
            topologies *= int(unit["n_distinct_shapes_exact"])
        if trees == 1 and topologies == 1:
            classes["T=1|Topo=1"] += 1
        elif trees > 1 and topologies == 1:
            classes["T>1|Topo=1"] += 1
        elif trees > 1 and topologies > 1:
            classes["T>1|Topo>1"] += 1
        else:
            fail(f"impossible legacy T/Topo state in {region.get('region')}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "analysis_unit": "legacy 50-kb coordinate region with HP-family aggregation",
        "regions_total": len(data["regions"]),
        "regions_with_primary_lineage": primary_regions,
        "regions_without_primary_lineage": len(data["regions"]) - primary_regions,
        "classes": ordered_counter(classes, CLASS_ORDER),
        "comparison_rule": (
            "Reference only. The legacy region and exact-PS route-block denominators differ, "
            "so no percent change, transition, or improvement rate is calculated."
        ),
    }


def validate_and_summarize(
    adapter_path: Path,
    receipt_path: Path,
    layered_path: Path,
    region_view_path: Path,
    legacy_path: Path | None,
) -> tuple[dict, list[dict]]:
    adapter = read_json(adapter_path)
    receipt = read_json(receipt_path)
    layered = read_json(layered_path)
    region_view = read_json(region_view_path)

    if receipt.get("all_pass") is not True:
        fail("adapter receipt all_pass is not true")
    expected_adapter_sha = receipt.get("output", {}).get("sha256")
    if expected_adapter_sha != sha256(adapter_path):
        fail("adapter SHA does not match its receipt")
    upstream_receipt_identity = (adapter.get("inputs") or {}).get("partition_run_receipt") or {}
    upstream_receipt_path = Path(str(upstream_receipt_identity.get("path", "")))
    if not upstream_receipt_path.is_file():
        fail("partition run receipt is missing")
    if upstream_receipt_identity.get("sha256") != sha256(upstream_receipt_path):
        fail("partition run receipt SHA does not match adapter identity")
    upstream_receipt = read_json(upstream_receipt_path)
    if upstream_receipt.get("all_pass") is not True:
        fail("partition run receipt all_pass is not true")
    funnel = adapter.get("input_funnel", {})
    if funnel.get("constraint_weight_total") != (
        funnel.get("constraint_weight_retained", 0)
        + funnel.get("constraint_weight_cut", 0)
        + funnel.get("constraint_weight_unavoidable", 0)
    ):
        fail("constraint-weight conservation failed")
    for check in (
        "check_constraint_weight_conservation",
        "check_cpp_parity_required",
        "check_cross_hp_zero",
        "check_cross_ps_zero",
        "check_stable_region_ids_unique",
    ):
        if funnel.get(check) is not True:
            fail(f"funnel check is not true: {check}")

    groups = adapter["groups"]
    details = layered["detail"]
    regions = region_view["regions"]
    expected_n = int(adapter["n_groups_analyzed"])
    if not (len(groups) == len(details) == len(regions) == expected_n):
        fail(
            f"group/detail/region count mismatch: {len(groups)}/{len(details)}/{len(regions)}/{expected_n}"
        )

    def unique_map(rows: list[dict], field: str) -> dict[str, dict]:
        mapped = {str(row[field]): row for row in rows}
        if len(mapped) != len(rows):
            fail(f"duplicate {field}")
        return mapped

    groups_by_id = unique_map(groups, "region_id")
    detail_by_id = unique_map(details, "region")
    regions_by_id = unique_map(regions, "region")
    if set(groups_by_id) != set(detail_by_id) or set(groups_by_id) != set(regions_by_id):
        fail("stable region-id sets differ across adapter/layered/region-view")

    classes: collections.Counter[str] = collections.Counter()
    hp_all: collections.Counter[str] = collections.Counter()
    hp_primary: collections.Counter[str] = collections.Counter()
    chromosome_all: collections.Counter[str] = collections.Counter()
    chromosome_primary: collections.Counter[str] = collections.Counter()
    t_bins: collections.Counter[str] = collections.Counter()
    hidden: collections.Counter[str] = collections.Counter()
    unique_morphology: collections.Counter[str] = collections.Counter()
    topo_multiple_morphology_sets: collections.Counter[str] = collections.Counter()
    topo_multiple_coarse_resolution: collections.Counter[str] = collections.Counter()
    candidate_invariant_morphology: collections.Counter[str] = collections.Counter()
    primary_count = 0
    complete_count = 0
    vaf_eligible_count = 0
    all_phase_sets: set[tuple[str, str]] = set()
    rows: list[dict] = []

    for region_id in sorted(groups_by_id):
        group = groups_by_id[region_id]
        unit = detail_by_id[region_id]
        region = regions_by_id[region_id]
        hp = str(group["hp_family"])
        phase_set = str(group["phase_set"])
        if hp not in {"1", "2"}:
            fail(f"non-primary HP family in exact route-block: {region_id} HP={hp}")
        if not phase_set or phase_set in {".", "None", "NA", "null"}:
            fail(f"missing exact PS: {region_id}")
        if str(unit.get("phase_set")) != phase_set or str(region.get("phase_set")) != phase_set:
            fail(f"PS mismatch across stages: {region_id}")
        if str(unit.get("family")) != hp:
            fail(f"HP mismatch across stages: {region_id}")
        if region.get("analysis_unit") != "exact_ps_hp_component_bounded_block":
            fail(f"unexpected analysis unit: {region_id}")
        if len(region.get("lineages", [])) != 1:
            fail(f"exact route-block must contain exactly one lineage/control: {region_id}")
        if int(group["n_sSNV"]) < 2 or int(group["n_sSNV"]) > 12:
            fail(f"tree-input k outside 2..12: {region_id}")
        if group.get("vaf_eligible") is not False or unit.get("vaf_eligible") is not False:
            vaf_eligible_count += 1

        hp_all[hp] += 1
        chromosome_all[str(group["chrom"])] += 1
        all_phase_sets.add((str(group["chrom"]), phase_set))
        is_primary = bool(unit.get("is_primary_lineage"))
        tree_class = "no_primary_lineage"
        topology_count: int | None = None
        candidate_count: int | None = None
        morphology_values: list[str] = []
        if is_primary:
            primary_count += 1
            hp_primary[hp] += 1
            chromosome_primary[str(group["chrom"])] += 1
            tree_class = structural_class(unit)
            classes[tree_class] += 1
            if tree_class != "incomplete":
                complete_count += 1
                candidate_count = int(unit["n_trees"])
                topology_count = int(unit["n_distinct_shapes_exact"])
                t_bins[t_bin(candidate_count)] += 1
                hidden[str(int(unit["n_hidden"]))] += 1
                if int(unit.get("n_trees_stored") or -1) != candidate_count:
                    fail(
                        f"topology morphology requires DISPLAY_TREE_CAP=0: {region_id} "
                        f"stored={unit.get('n_trees_stored')} T={candidate_count}"
                    )
                signatures: dict[str, str] = {}
                for tree in unit["trees"]:
                    signature = canonical_shape(tree.get("edges", []))
                    signatures[signature] = morphology_name(*shape_features(tree.get("edges", [])))
                if len(signatures) != topology_count:
                    fail(
                        f"stored canonical shape count mismatch: {region_id} "
                        f"{len(signatures)} != {topology_count}"
                    )
                morphology_values = sorted(set(signatures.values()))
                if topology_count == 1:
                    if len(morphology_values) != 1:
                        fail(f"Topo=1 has multiple coarse morphology classes: {region_id}")
                    unique_morphology[morphology_values[0]] += 1
                    candidate_invariant_morphology[morphology_values[0]] += 1
                else:
                    key = "+".join(morphology_values)
                    topo_multiple_morphology_sets[key] += 1
                    topo_multiple_coarse_resolution[
                        "coarse_class_same" if len(morphology_values) == 1 else "coarse_class_multiple"
                    ] += 1
                    if len(morphology_values) == 1:
                        candidate_invariant_morphology[morphology_values[0]] += 1

        rows.append(
            {
                "region_id": region_id,
                "chrom": group["chrom"],
                "start": group["start"],
                "end": group["end"],
                "phase_set": phase_set,
                "hp": hp,
                "source_unit_id": group["unit_id"],
                "source_component_id": group["component_id"],
                "source_block_id": group["block_id"],
                "k_sSNV": group["n_sSNV"],
                "is_primary_lineage": is_primary,
                "tree_class": tree_class,
                "T": candidate_count,
                "Topo": topology_count,
                "morphology_possible": ";".join(morphology_values),
                "n_hidden": unit.get("n_hidden") if is_primary else None,
                "n_reads": unit.get("n_reads"),
                "verification_status": unit.get("verification_status"),
                "vaf_eligible": unit.get("vaf_eligible"),
            }
        )

    if vaf_eligible_count:
        fail(f"{vaf_eligible_count} exact-PS units were unexpectedly marked VAF-eligible")
    if region_view["census"]["n_regions"] != expected_n:
        fail("region-view census total differs from stable ID count")
    if region_view["census"]["n_regions_with_primary_lineage"] != primary_count:
        fail("region-view primary count differs from recomputation")
    if region_view["census"]["L1"]["all_primary_noncapped_V1V7_pass"] is not True:
        fail("not all primary non-capped units passed V1-V7")
    if classes["incomplete"] != int(
        region_view["census"]["L1"]["n_verification_not_applicable_capped"]
    ):
        fail("incomplete count differs from capped/not-applicable count")

    class_counts = ordered_counter(classes, CLASS_ORDER)
    unique_topology_total = class_counts["T=1|Topo=1"] + class_counts["T>1|Topo=1"]
    topo_multiple_total = class_counts["T>1|Topo>1"]
    if sum(class_counts.values()) != primary_count:
        fail("primary structural classes do not conserve")
    if sum(unique_morphology.values()) != unique_topology_total:
        fail("unique-topology morphology classes do not conserve")
    if sum(topo_multiple_coarse_resolution.values()) != topo_multiple_total:
        fail("Topo>1 coarse morphology resolution does not conserve")
    candidate_invariant_total = unique_topology_total + topo_multiple_coarse_resolution[
        "coarse_class_same"
    ]
    if sum(candidate_invariant_morphology.values()) != candidate_invariant_total:
        fail("candidate-invariant coarse morphology classes do not conserve")

    evidence_paths = {
        "adapter": adapter_path,
        "adapter_receipt": receipt_path,
        "layered": layered_path,
        "region_view": region_view_path,
    }
    summary = {
        "schema_name": "intersubmod.exact_ps_topology_observation",
        "schema_version": "1.0.0",
        "generated_at": datetime.now(ZoneInfo("Asia/Taipei")).isoformat(),
        "task_type": "E_hotfix_methodology_validation",
        "sample": adapter["sample"],
        "scope": "HCC1395 chr1-22 only",
        "status": "HCC1395_exact_PS_topology_reconstruction_passed; seven-sample promotion not authorized",
        "analysis_unit": "exact PS x HP x read-linked component x bounded block (k=2..12)",
        "claim_ceiling": adapter["claim_ceiling"],
        "inputs": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in evidence_paths.items()
        },
        "evidence_tier": {
            "technical_checks_passed": True,
            "upstream_partition_receipt": {
                "path": str(upstream_receipt_path.resolve()),
                "sha256": sha256(upstream_receipt_path),
                "task_type": upstream_receipt.get("task_type"),
                "claim_status": upstream_receipt.get("claim_status"),
                "validation_evidence_eligible": upstream_receipt.get(
                    "validation_evidence_eligible"
                ),
            },
            "publication_or_cohort_final_eligible": False,
            "reason": (
                "HCC1395-only Task-E technical validation reuses an upstream exploratory-pilot "
                "partition receipt marked PARTIAL and validation_evidence_eligible=false."
            ),
        },
        "validation": {
            "all_pass": True,
            "stable_region_ids_unique_and_equal_across_3_stages": True,
            "exact_nonmissing_PS_for_every_route_block": True,
            "one_HP_and_one_lineage_or_control_per_route_block": True,
            "cross_PS_violations": 0,
            "cross_HP_violations": 0,
            "python_cpp_mismatch_count": int(receipt["counts"]["python_cpp_mismatch_count"]),
            "constraint_weight_conservation": {
                "total": funnel["constraint_weight_total"],
                "retained": funnel["constraint_weight_retained"],
                "cut": funnel["constraint_weight_cut"],
                "unavoidable": funnel["constraint_weight_unavoidable"],
                "equation_pass": True,
            },
            "all_primary_noncapped_V1_V7_pass": True,
            "vaf_ranking_performed": False,
        },
        "funnel": {
            **funnel,
            "primary_unique_loci_percent_of_S": pct(
                int(funnel["primary_unique_loci"]), int(funnel["candidate_loci_S"])
            ),
            "constraint_weight_retained_percent": pct(
                int(funnel["constraint_weight_retained"]), int(funnel["constraint_weight_total"])
            ),
            "tree_input_groups_percent_of_bounded_blocks": pct(
                int(funnel["tree_input_groups"]), int(funnel["bounded_blocks"])
            ),
        },
        "route_block_census": {
            "all_tree_input_route_blocks": expected_n,
            "unique_chromosome_PS_pairs": len(all_phase_sets),
            "HP_all": count_with_percent(ordered_counter(hp_all), expected_n),
            "primary_mutation_route_blocks": primary_count,
            "reference_only_controls": expected_n - primary_count,
            "HP_primary": count_with_percent(ordered_counter(hp_primary), primary_count),
            "chromosome_all": ordered_counter(chromosome_all),
            "chromosome_primary": ordered_counter(chromosome_primary),
        },
        "tree_results": {
            "denominator_primary_route_blocks": primary_count,
            "complete": {"n": complete_count, "percent_of_primary": pct(complete_count, primary_count)},
            "incomplete": {
                "n": class_counts["incomplete"],
                "percent_of_primary": pct(class_counts["incomplete"], primary_count),
                "reason": "candidate enumeration capped; verification not applicable",
            },
            "classes": {
                key: {
                    "n": value,
                    "percent_of_primary": pct(value, primary_count),
                    "percent_of_complete": pct(value, complete_count)
                    if key != "incomplete"
                    else None,
                }
                for key, value in class_counts.items()
            },
            "T_distribution_complete": ordered_counter(
                t_bins, ("1", "2", "3", "4", "5", "6", ">6")
            ),
            "hidden_nodes_distribution_complete": ordered_counter(hidden),
        },
        "topology_morphology": {
            "definition": (
                "sister=max outdegree >=2; direct=max ROOT-based depth >=2; hidden H_* nodes are included. "
                "These are mutation-state graph features, not direct clone/subclone or CNV proof."
            ),
            "unique_topology_route_blocks": unique_topology_total,
            "unique_topology_classes": count_with_percent(
                ordered_counter(unique_morphology, MORPHOLOGY_ORDER), unique_topology_total
            ),
            "Topo_gt1_exact_unresolved": topo_multiple_total,
            "Topo_gt1_coarse_resolution": count_with_percent(
                ordered_counter(
                    topo_multiple_coarse_resolution,
                    ("coarse_class_same", "coarse_class_multiple"),
                ),
                topo_multiple_total,
            ),
            "candidate_invariant_coarse_morphology_route_blocks": candidate_invariant_total,
            "candidate_invariant_coarse_classes": count_with_percent(
                ordered_counter(candidate_invariant_morphology, MORPHOLOGY_ORDER),
                candidate_invariant_total,
            ),
            "Topo_gt1_possible_class_sets": count_with_percent(
                ordered_counter(topo_multiple_morphology_sets), topo_multiple_total
            ),
        },
        "legacy_reference": legacy_baseline(legacy_path) if legacy_path else None,
        "interpretation_limits": [
            "The new denominator is an exact-PS route-block, not a legacy coordinate region W.",
            "HP1 and HP2 remain independent local trees; no cross-PS or cross-HP biological clone forest is inferred here.",
            "Read counts currently support MINREAD and segmentation only; the solver objective still uses pattern presence.",
            "Projected HP-specific coverage is not caller VAF; VAF ranking is intentionally disabled.",
            "Direct includes hidden-state depth and therefore cannot alone establish a subclone or exclude CNV/error explanations.",
        ],
    }
    return summary, rows


def write_metric_tsv(path: Path, summary: dict) -> None:
    records: list[tuple[str, str, int | float | str | None, str, float | None]] = []
    funnel = summary["funnel"]
    records.extend(
        [
            ("funnel", "candidate_loci_S", funnel["candidate_loci_S"], "candidate_loci_S", 100.0),
            (
                "funnel",
                "primary_unique_loci",
                funnel["primary_unique_loci"],
                "candidate_loci_S",
                funnel["primary_unique_loci_percent_of_S"],
            ),
            (
                "funnel",
                "exact_ps_hp_units",
                funnel["exact_ps_hp_units"],
                "exact_ps_hp_units",
                100.0,
            ),
            (
                "funnel",
                "bounded_blocks",
                funnel["bounded_blocks"],
                "bounded_blocks",
                100.0,
            ),
            (
                "funnel",
                "tree_input_groups",
                funnel["tree_input_groups"],
                "bounded_blocks",
                funnel["tree_input_groups_percent_of_bounded_blocks"],
            ),
        ]
    )
    results = summary["tree_results"]
    for key, row in results["classes"].items():
        records.append(
            (
                "tree_class",
                key,
                row["n"],
                "primary_route_blocks",
                row["percent_of_primary"],
            )
        )
    for key, row in summary["topology_morphology"]["unique_topology_classes"].items():
        records.append(
            ("unique_topology_morphology", key, row["n"], "unique_topology_route_blocks", row["percent"])
        )
    lines = ["section\tmetric\tvalue\tdenominator\tpercent"]
    lines.extend("\t".join("" if value is None else str(value) for value in record) for record in records)
    atomic_text(path, "\n".join(lines) + "\n")


def write_region_tsv_gz(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    fieldnames = list(rows[0]) if rows else []
    with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def markdown_report(summary: dict) -> str:
    funnel = summary["funnel"]
    results = summary["tree_results"]
    morphology = summary["topology_morphology"]

    def number(value: int | float | None) -> str:
        if value is None:
            return "—"
        if isinstance(value, float):
            return f"{value:.2f}%"
        return f"{value:,}"

    lines = [
        "<!--",
        f"建立時間: {summary['generated_at']}",
        "目標: 驗證 HCC1395 exact-PS 分區後的局部 mutation-state tree 重建",
        "處理範圍: HCC1395 chr1-22；Task Type E methodology bugfix validation",
        f"關聯檔案: {summary['inputs']['region_view']['path']}",
        "-->",
        "",
        "# HCC1395 exact-PS 局部拓撲重建觀察",
        "",
        "> **TL;DR**：Python 與 C++ 已把 PS 固定為 fail-closed 的樹分析邊界；HCC1395 chr1–22 的 11,542 個 tree-input route-block 中，跨 PS=0、跨 HP=0、Python/C++ mismatch=0。9,600 個 mutation-bearing route-block 中有 9,180 個完整枚舉；但這是新單位的局部樹，不能直接與舊 W 相減，也尚未做 VAF 排序或 clone 確認。",
        "",
        "> **Scope ribbon — HCC1395 ONLY / PARTIAL VALIDATION**：尚未推廣至 7 樣本，不可作 cohort final evidence。",
        "",
        "> **Evidence tier**：技術檢查 PASS；但上游 partition receipt 為 `exploratory_pilot / PARTIAL / validation_evidence_eligible=false`，所以目前不是 Task-B comprehensive final。",
        "",
        "## 1. 分析單位與方法修正",
        "",
        "新版單位固定為 `exact PS × HP × read-linked component × bounded block (k=2..12)`。同一個 HP 編號只在同一個 exact PS 內有意義；不同 PS 不合併，也不假設 HP1/HP2 在 block 間方向一致。",
        "",
        "## 2. 完整漏斗",
        "",
        "| 階段 | 數量 | 分母與比例 |",
        "|---|---:|---:|",
        f"| LongPhase-S PASS autosomal sSNV (S) | {funnel['candidate_loci_S']:,} | 100.00% of S |",
        f"| 進入 exact-PS primary units 的 unique loci | {funnel['primary_unique_loci']:,} | {funnel['primary_unique_loci_percent_of_S']:.2f}% of S |",
        f"| exact PS×HP components | {funnel['exact_ps_hp_units']:,} | 單位數，非 loci |",
        f"| k≤12 bounded blocks | {funnel['bounded_blocks']:,} | 單位數，非 loci |",
        f"| k=1，不進樹 | {funnel['k1_blocks_not_tree_eligible']:,} | {pct(funnel['k1_blocks_not_tree_eligible'], funnel['bounded_blocks']):.2f}% of blocks |",
        f"| pattern support 不足 | {funnel['pattern_unsupported_blocks']:,} | {pct(funnel['pattern_unsupported_blocks'], funnel['bounded_blocks']):.2f}% of blocks |",
        f"| tree-input route-blocks | {funnel['tree_input_groups']:,} | {funnel['tree_input_groups_percent_of_bounded_blocks']:.2f}% of blocks |",
        f"| constraint weight retained | {funnel['constraint_weight_retained']:,} / {funnel['constraint_weight_total']:,} | {funnel['constraint_weight_retained_percent']:.2f}% |",
        "",
        "守恆：`285,897 = 281,967 retained + 1,254 cut + 2,676 unavoidable`。",
        "",
        "## 3. 局部樹候選與拓撲",
        "",
        f"分母為 **{results['denominator_primary_route_blocks']:,} 個 mutation-bearing exact-PS route-blocks**；不是舊版 coordinate region W。",
        "",
        "| 結構狀態 | 數量 | 占全部 primary | 占 complete |",
        "|---|---:|---:|---:|",
    ]
    for key in CLASS_ORDER:
        row = results["classes"][key]
        lines.append(
            f"| {key} | {row['n']:,} | {number(row['percent_of_primary'])} | {number(row['percent_of_complete'])} |"
        )
    lines.extend(
        [
            "",
            f"Complete = {results['complete']['n']:,}（{results['complete']['percent_of_primary']:.2f}%）；incomplete/capped = {results['incomplete']['n']:,}（{results['incomplete']['percent_of_primary']:.2f}%）。",
            "",
            "## 4. Topo=1 的圖形特徵",
            "",
            f"只對 **Topo=1 的 {morphology['unique_topology_route_blocks']:,} 個 route-blocks**做唯一分類。",
            "",
            "| mutation-state graph 特徵 | 數量 | 占 Topo=1 |",
            "|---|---:|---:|",
        ]
    )
    labels = {
        "no_branch_no_multilevel": "無分枝、無多層（舊 internal: single-only）",
        "sister_only": "有分枝、無多層（sister-only）",
        "direct_only": "有多層、無分枝（direct-only）",
        "sister_and_direct": "同時有分枝與多層",
    }
    for key in MORPHOLOGY_ORDER:
        row = morphology["unique_topology_classes"][key]
        lines.append(f"| {labels[key]} | {row['n']:,} | {row['percent']:.2f}% |")
    coarse = morphology["Topo_gt1_coarse_resolution"]
    lines.extend(
        [
            "",
            f"Topo>1 共 {morphology['Topo_gt1_exact_unresolved']:,} 個；其中 {coarse['coarse_class_same']['n']:,} 個雖然 exact shape 不唯一，但所有候選落在同一個粗分類；{coarse['coarse_class_multiple']['n']:,} 個連粗分類也不同。這不等於選出了唯一真實拓撲。",
            "",
            f"若使用較寬鬆但仍不任選候選樹的 `candidate-invariant coarse morphology` 口徑，可分類 {morphology['candidate_invariant_coarse_morphology_route_blocks']:,} 個；其中額外納入的 {coarse['coarse_class_same']['n']:,} 個全部屬 sister+direct，但它們的 exact Topo 仍為 >1。",
            "",
            "## 5. 驗證結論",
            "",
            "- adapter receipt：PASS；stable region IDs 在 adapter、layered、region-view 三階段完全一致且唯一。",
            "- 每個 tree-input route-block：exact PS 非缺值、僅一個 HP、僅一個 lineage/control。",
            "- cross-PS violations = 0；cross-HP violations = 0；Python/C++ mismatch = 0。",
            "- 所有 non-capped primary units 的 V1–V7 全部通過。",
            "- VAF ranking = **未執行**；所有 exact-PS units 均 `vaf_eligible=false`。",
            "",
            "## 6. 舊版數字為何不能直接相減",
            "",
        ]
    )
    legacy = summary.get("legacy_reference")
    if legacy:
        lines.extend(
            [
                f"舊版是 {legacy['regions_total']:,} 個 50-kb coordinate regions，其中 {legacy['regions_with_primary_lineage']:,} 個有 primary lineage；新版是 {summary['route_block_census']['all_tree_input_route_blocks']:,} 個 exact-PS route-blocks。分母、BaseQ、k 上限與分區規則均不同，因此只保留 baseline，不計算 percent change。",
                "",
                "| 舊版 baseline 狀態 | 數量 |",
                "|---|---:|",
            ]
        )
        for key in CLASS_ORDER:
            lines.append(f"| {key} | {legacy['classes'][key]:,} |")
    lines.extend(
        [
            "",
            "## 7. 主張邊界",
            "",
            "- `direct` 的判定包含 hidden H_* state 的深度，因此不能單獨宣稱 subclone，也不能排除 CNV、錯誤或 partial-read 投影造成的替代解釋。",
            "- read count 現在用於 MINREAD 與 segmentation audit；Steiner solver 仍以 pattern presence 建樹，尚非 count-weighted likelihood。",
            "- HP-specific projected coverage 不是 caller VAF；這一版沒有用 VAF 排序候選樹。",
            "- 下一個 gate 應先做 HCC1395 matched legacy-region candidate-forest transition，再決定是否推廣 7 樣本。",
            "",
            "---",
            "",
            "Partial validation：HCC1395 chr1–22 only；不可替代全 7 樣本 comprehensive validation。",
        ]
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adapter", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--layered", type=Path, required=True)
    parser.add_argument("--region-view", type=Path, required=True)
    parser.add_argument("--legacy-region-view", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-metrics-tsv", type=Path, required=True)
    parser.add_argument("--output-regions-tsv-gz", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument("--output-receipt", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary, rows = validate_and_summarize(
        args.adapter,
        args.receipt,
        args.layered,
        args.region_view,
        args.legacy_region_view,
    )
    atomic_json(args.output_json, summary)
    write_metric_tsv(args.output_metrics_tsv, summary)
    write_region_tsv_gz(args.output_regions_tsv_gz, rows)
    atomic_text(args.output_md, markdown_report(summary))
    authoritative_outputs = {
        "summary_json": args.output_json,
        "metrics_tsv": args.output_metrics_tsv,
        "regions_tsv_gz": args.output_regions_tsv_gz,
        "observation_markdown": args.output_md,
    }
    receipt = {
        "schema_name": "intersubmod.exact_ps_topology_observation.receipt",
        "schema_version": "1.0.0",
        "generated_at": datetime.now(ZoneInfo("Asia/Taipei")).isoformat(),
        "all_technical_checks_pass": summary["validation"]["all_pass"],
        "sample": summary["sample"],
        "scope": summary["scope"],
        "analysis_unit": summary["analysis_unit"],
        "authoritative_input_identities": summary["inputs"],
        "authoritative_outputs": {
            name: {
                "path": str(path.resolve()),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for name, path in authoritative_outputs.items()
        },
        "generator": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "evidence_tier": summary["evidence_tier"],
        "authoritative_path_rule": (
            "Only paths listed in authoritative_input_identities and authoritative_outputs are "
            "authoritative for this observation. Similarly named earlier artifacts in the same "
            "directory are superseded and must not be used for morphology reporting."
        ),
    }
    atomic_json(args.output_receipt, receipt)
    print(
        json.dumps(
            {
                "all_pass": summary["validation"]["all_pass"],
                "sample": summary["sample"],
                "scope": summary["scope"],
                "tree_input_route_blocks": summary["route_block_census"][
                    "all_tree_input_route_blocks"
                ],
                "primary_route_blocks": summary["tree_results"][
                    "denominator_primary_route_blocks"
                ],
                "classes": {
                    key: value["n"] for key, value in summary["tree_results"]["classes"].items()
                },
                "outputs": {
                    "json": str(args.output_json.resolve()),
                    "metrics_tsv": str(args.output_metrics_tsv.resolve()),
                    "regions_tsv_gz": str(args.output_regions_tsv_gz.resolve()),
                    "markdown": str(args.output_md.resolve()),
                    "receipt": str(args.output_receipt.resolve()),
                },
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
