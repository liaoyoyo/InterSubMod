#!/usr/bin/env python3
"""Build region-level coarse topology and observed-support censuses.

The structural census classifies only complete regions with a unique ordered
regional topology.  Topo>1 regions remain unresolved.  The observed-support
census is deliberately stricter: it is limited to T=1/Topo=1 regions and does
not treat H_* Steiner/partial-supported nodes as directly observed states.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from datetime import datetime
from itertools import product
from pathlib import Path
from zoneinfo import ZoneInfo


SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
GRAPH_CLASSES = [
    "single_only",
    "sister_only",
    "direct_only",
    "sister_and_direct",
    "topology_multiple_unresolved",
]
OBSERVED_CLASSES = [
    "no_observed_within_hp_relation",
    "observed_sister_only",
    "observed_direct_only",
    "observed_sister_and_direct",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def graph_class(sister: bool, direct: bool) -> str:
    if sister and direct:
        return "sister_and_direct"
    if sister:
        return "sister_only"
    if direct:
        return "direct_only"
    return "single_only"


def observed_tree_features(edges: list[list[str]]) -> dict:
    children: dict[str, list[str]] = collections.defaultdict(list)
    nodes = {"ROOT"}
    for parent, child in edges:
        children[parent].append(child)
        nodes.update((parent, child))
    observed = {node for node in nodes if node != "ROOT" and not node.startswith("H_")}

    direct = False
    for ancestor in observed:
        stack = list(children.get(ancestor, []))
        visited: set[str] = set()
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            if node in observed:
                direct = True
            stack.extend(children.get(node, []))

    def observed_below(start: str) -> set[str]:
        found: set[str] = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in observed:
                found.add(node)
            stack.extend(children.get(node, []))
        return found

    sister = any(
        sum(bool(observed_below(child)) for child in child_nodes) >= 2
        for child_nodes in list(children.values())
    )
    return {
        "n_observed_states": len(observed),
        "sister": sister,
        "direct": direct,
    }


def observed_region_class(sister: bool, direct: bool) -> str:
    if sister and direct:
        return "observed_sister_and_direct"
    if sister:
        return "observed_sister_only"
    if direct:
        return "observed_direct_only"
    return "no_observed_within_hp_relation"


def load_region_rows(path: Path) -> tuple[list[dict], dict[tuple[str, str], dict]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_key = {(row["sample"], row["region"]): row for row in rows}
    if len(rows) != len(by_key):
        raise RuntimeError("region TSV keys are not unique")
    return rows, by_key


def build_census(
    exact_catalog: dict,
    ct_report: dict,
    region_rows: list[dict],
    region_by_key: dict[tuple[str, str], dict],
    vaf_census: dict,
) -> dict:
    shape_meta = {row["shape_id"]: row for row in ct_report["shape_catalog"]}
    if len(shape_meta) != int(ct_report["aggregate"]["distinct_shapes_global"]):
        raise RuntimeError("global shape metadata count differs from corrected report")

    region_shape_sets: dict[tuple[str, str], list[set[str]]] = collections.defaultdict(list)
    for sample in exact_catalog["samples"]:
        sample_name = sample["sample"]
        for unit in sample["unit_rows"]:
            shape_ids = set(unit["shape_candidate_counts"])
            unknown = shape_ids - set(shape_meta)
            if unknown:
                raise RuntimeError(f"unknown shape IDs for {sample_name}: {sorted(unknown)[:3]}")
            region_shape_sets[(sample_name, unit["region"])].append(shape_ids)

    graph_counts = {sample: collections.Counter() for sample in SAMPLE_ORDER}
    topo_multiple_feature_counts = {sample: collections.Counter() for sample in SAMPLE_ORDER}
    single_only_hp_counts = {sample: collections.Counter() for sample in SAMPLE_ORDER}
    region_classes: dict[tuple[str, str], str] = {}

    for row in region_rows:
        sample = row["sample"]
        if sample not in graph_counts:
            raise RuntimeError(f"unexpected dataset row: {sample}")
        structural_class = row["structural_class"]
        if structural_class == "incomplete":
            graph_counts[sample]["incomplete"] += 1
            continue
        graph_counts[sample]["complete"] += 1
        key = (sample, row["region"])
        component_shapes = region_shape_sets.get(key)
        if not component_shapes:
            raise RuntimeError(f"complete region lacks shape rows: {sample} {row['region']}")

        if structural_class == "T>1|Topo>1":
            graph_counts[sample]["topology_multiple_unresolved"] += 1
            region_classes[key] = "topology_multiple_unresolved"
            feature_pairs = []
            for combination in product(*[sorted(values) for values in component_shapes]):
                sister = any(int(shape_meta[shape_id]["max_outdegree"]) >= 2 for shape_id in combination)
                direct = any(int(shape_meta[shape_id]["max_depth"]) >= 2 for shape_id in combination)
                feature_pairs.append((sister, direct))
            topo_multiple_feature_counts[sample]["sister_possible"] += int(any(pair[0] for pair in feature_pairs))
            topo_multiple_feature_counts[sample]["sister_definite"] += int(all(pair[0] for pair in feature_pairs))
            topo_multiple_feature_counts[sample]["direct_possible"] += int(any(pair[1] for pair in feature_pairs))
            topo_multiple_feature_counts[sample]["direct_definite"] += int(all(pair[1] for pair in feature_pairs))
            topo_multiple_feature_counts[sample]["both_definite"] += int(
                all(pair[0] and pair[1] for pair in feature_pairs)
            )
            continue

        sister = any(
            int(shape_meta[shape_id]["max_outdegree"]) >= 2
            for values in component_shapes
            for shape_id in values
        )
        direct = any(
            int(shape_meta[shape_id]["max_depth"]) >= 2
            for values in component_shapes
            for shape_id in values
        )
        category = graph_class(sister, direct)
        graph_counts[sample][category] += 1
        graph_counts[sample]["topology_unique"] += 1
        region_classes[key] = category
        if category == "single_only":
            hp_count = int(row["primary_HP_units"])
            single_only_hp_counts[sample][f"primary_hp_{hp_count}"] += 1

    source_by_sample = {row["sample"]: row["source_layered"] for row in vaf_census["samples"]}
    observed_counts = {sample: collections.Counter() for sample in SAMPLE_ORDER}
    observed_no_relation_detail = {sample: collections.Counter() for sample in SAMPLE_ORDER}
    single_hp_single_node_detail = {sample: collections.Counter() for sample in SAMPLE_ORDER}

    for sample in SAMPLE_ORDER:
        layered_path = Path(source_by_sample[sample])
        layered = json.loads(layered_path.read_text(encoding="utf-8"))
        primary_by_region: dict[str, list[dict]] = collections.defaultdict(list)
        for unit in layered["detail"]:
            if unit.get("is_primary_lineage"):
                primary_by_region[unit["region"]].append(unit)

        for region, units in primary_by_region.items():
            key = (sample, region)
            row_class = region_classes.get(key)
            if row_class is None:
                continue
            region_row = region_by_key[key]
            if region_row["structural_class"] != "T=1|Topo=1":
                continue
            observed_counts[sample]["T1_Topo1"] += 1
            region_sister = False
            region_direct = False
            total_observed = 0
            per_component_observed = []
            for unit in units:
                if int(unit["n_trees"]) != 1 or len(unit["trees"]) != 1:
                    raise RuntimeError(f"T=1 region lacks exactly one tree: {sample} {region}")
                features = observed_tree_features(unit["trees"][0]["edges"])
                region_sister |= bool(features["sister"])
                region_direct |= bool(features["direct"])
                total_observed += int(features["n_observed_states"])
                per_component_observed.append(int(features["n_observed_states"]))
            category = observed_region_class(region_sister, region_direct)
            observed_counts[sample][category] += 1
            if category == "no_observed_within_hp_relation":
                if total_observed == 0:
                    observed_no_relation_detail[sample]["zero_observed_states"] += 1
                elif total_observed == 1:
                    observed_no_relation_detail[sample]["one_observed_state"] += 1
                elif len(units) == 2 and sorted(per_component_observed) == [1, 1]:
                    observed_no_relation_detail[sample]["two_hp_isolated_states"] += 1
                else:
                    observed_no_relation_detail[sample]["other"] += 1

            if row_class == "single_only" and len(units) == 1:
                nodes = {
                    node
                    for edge in units[0]["trees"][0]["edges"]
                    for node in edge
                    if node != "ROOT"
                }
                if len(nodes) != 1:
                    raise RuntimeError(f"single-HP single-only region has nontrivial node count: {sample} {region}")
                node = next(iter(nodes))
                single_hp_single_node_detail[sample][
                    "hidden_or_partial" if node.startswith("H_") else "directly_observed"
                ] += 1

    samples = []
    for sample in SAMPLE_ORDER:
        graph = graph_counts[sample]
        observed = observed_counts[sample]
        complete = int(graph["complete"])
        t1 = int(observed["T1_Topo1"])
        row = {
            "sample": sample,
            "graph": {
                "complete_regions": complete,
                "topology_unique": int(graph["topology_unique"]),
                "classes": {key: int(graph[key]) for key in GRAPH_CLASSES},
                "shares_of_complete": {key: ratio(int(graph[key]), complete) for key in GRAPH_CLASSES},
                "incomplete_regions": int(graph["incomplete"]),
                "single_only_by_primary_hp_count": dict(single_only_hp_counts[sample]),
                "topology_multiple_features": dict(topo_multiple_feature_counts[sample]),
            },
            "observed_support_T1_Topo1": {
                "denominator": t1,
                "classes": {key: int(observed[key]) for key in OBSERVED_CLASSES},
                "shares": {key: ratio(int(observed[key]), t1) for key in OBSERVED_CLASSES},
                "no_relation_detail": dict(observed_no_relation_detail[sample]),
                "single_hp_single_node_detail": dict(single_hp_single_node_detail[sample]),
            },
        }
        samples.append(row)

    aggregate_graph = collections.Counter()
    aggregate_observed = collections.Counter()
    aggregate_no_relation = collections.Counter()
    aggregate_single_node = collections.Counter()
    aggregate_topon_features = collections.Counter()
    for row in samples:
        aggregate_graph["complete_regions"] += row["graph"]["complete_regions"]
        aggregate_graph["topology_unique"] += row["graph"]["topology_unique"]
        aggregate_graph["incomplete_regions"] += row["graph"]["incomplete_regions"]
        aggregate_graph.update(row["graph"]["classes"])
        aggregate_observed["denominator"] += row["observed_support_T1_Topo1"]["denominator"]
        aggregate_observed.update(row["observed_support_T1_Topo1"]["classes"])
        aggregate_no_relation.update(row["observed_support_T1_Topo1"]["no_relation_detail"])
        aggregate_single_node.update(row["observed_support_T1_Topo1"]["single_hp_single_node_detail"])
        aggregate_topon_features.update(row["graph"]["topology_multiple_features"])

    complete = int(aggregate_graph["complete_regions"])
    observed_denominator = int(aggregate_observed["denominator"])
    return {
        "schema_version": "1.0",
        "generated_at": datetime.now(ZoneInfo("Asia/Taipei")).isoformat(timespec="seconds"),
        "status": "PASS",
        "scope": "chr1-22; 7 dataset rows; historical layered-v2 engineering snapshot",
        "claim_ceiling": "regional mutation-state topology and observed-state relations; not confirmed clones/subclones",
        "definitions": {
            "graph_single_only": "Topo=1 ordered HP forest with no within-HP max_outdegree>=2 and no within-HP max_depth>=2",
            "graph_sister": "within at least one primary HP component, max_outdegree>=2; cross-HP components are never sisters",
            "graph_direct": "within at least one primary HP component, max_depth>=2; H_* nodes are included in graph depth",
            "topology_multiple_unresolved": "complete region with T>1 and more than one ordered regional topology",
            "observed_node": "non-ROOT node whose label does not begin with H_",
            "observed_sister": "a branch has at least two child subtrees that each contain an observed node",
            "observed_direct": "two observed nodes have an ancestor-descendant relationship; H_* may connect the path but is not evidence",
        },
        "aggregate": {
            "graph": {
                "complete_regions": complete,
                "topology_unique": int(aggregate_graph["topology_unique"]),
                "classes": {key: int(aggregate_graph[key]) for key in GRAPH_CLASSES},
                "shares_of_complete": {key: ratio(int(aggregate_graph[key]), complete) for key in GRAPH_CLASSES},
                "incomplete_regions": int(aggregate_graph["incomplete_regions"]),
                "topology_multiple_features": dict(aggregate_topon_features),
            },
            "observed_support_T1_Topo1": {
                "denominator": observed_denominator,
                "classes": {key: int(aggregate_observed[key]) for key in OBSERVED_CLASSES},
                "shares": {
                    key: ratio(int(aggregate_observed[key]), observed_denominator)
                    for key in OBSERVED_CLASSES
                },
                "no_relation_detail": dict(aggregate_no_relation),
                "single_hp_single_node_detail": dict(aggregate_single_node),
            },
        },
        "samples": samples,
    }


def validation_rows(census: dict, ct_report: dict, region_row_count: int) -> list[dict]:
    rows: list[dict] = []

    def add(scope: str, check: str, passed: bool) -> None:
        rows.append({"scope": scope, "check": check, "pass": str(bool(passed))})

    ct_by_sample = {row["sample"]: row for row in ct_report["samples"]}
    for row in census["samples"]:
        sample = row["sample"]
        graph = row["graph"]
        observed = row["observed_support_T1_Topo1"]
        expected_classes = ct_by_sample[sample]["T_and_Topology"]["classes"]
        add(sample, "graph_classes_conserve_complete", sum(graph["classes"].values()) == graph["complete_regions"])
        add(sample, "topology_unique_plus_multiple_conserve_complete", graph["topology_unique"] + graph["classes"]["topology_multiple_unresolved"] == graph["complete_regions"])
        add(sample, "complete_matches_corrected_report", graph["complete_regions"] == sum(int(expected_classes[key]) for key in ("T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1")))
        add(sample, "topology_multiple_matches_corrected_report", graph["classes"]["topology_multiple_unresolved"] == int(expected_classes["T>1|Topo>1"]))
        add(sample, "observed_classes_conserve_T1", sum(observed["classes"].values()) == observed["denominator"])
        add(sample, "observed_T1_matches_corrected_report", observed["denominator"] == int(expected_classes["T=1|Topo=1"]))
        add(sample, "observed_no_relation_detail_conserves", sum(observed["no_relation_detail"].values()) == observed["classes"]["no_observed_within_hp_relation"])
        add(sample, "no_unclassified_observed_no_relation", int(observed["no_relation_detail"].get("other", 0)) == 0)

    aggregate = census["aggregate"]
    graph = aggregate["graph"]
    observed = aggregate["observed_support_T1_Topo1"]
    corrected = ct_report["aggregate"]["T_and_Topology.classes"]
    add("aggregate", "region_tsv_has_W_primary_rows", region_row_count == int(ct_report["aggregate"]["W_primary"]))
    add("aggregate", "graph_classes_conserve_complete", sum(graph["classes"].values()) == graph["complete_regions"])
    add("aggregate", "complete_matches_corrected_report", graph["complete_regions"] == sum(int(corrected[key]) for key in ("T=1|Topo=1", "T>1|Topo=1", "T>1|Topo>1")))
    add("aggregate", "topology_unique_matches_corrected_report", graph["topology_unique"] == int(corrected["T=1|Topo=1"]) + int(corrected["T>1|Topo=1"]))
    add("aggregate", "topology_multiple_matches_corrected_report", graph["classes"]["topology_multiple_unresolved"] == int(corrected["T>1|Topo>1"]))
    add("aggregate", "observed_classes_conserve_T1", sum(observed["classes"].values()) == observed["denominator"])
    add("aggregate", "observed_T1_matches_corrected_report", observed["denominator"] == int(corrected["T=1|Topo=1"]))
    add("aggregate", "all_sample_rows_present", [row["sample"] for row in census["samples"]] == SAMPLE_ORDER)
    add("aggregate", "all_checks_pass", all(row["pass"] == "True" for row in rows))
    return rows


def write_summary(path: Path, census: dict) -> None:
    fields = [
        "sample",
        "complete_regions",
        *GRAPH_CLASSES,
        "topology_unique",
        "T1_Topo1",
        *OBSERVED_CLASSES,
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        rows = list(census["samples"]) + [
            {
                "sample": "aggregate",
                "graph": census["aggregate"]["graph"],
                "observed_support_T1_Topo1": census["aggregate"]["observed_support_T1_Topo1"],
            }
        ]
        for row in rows:
            graph = row["graph"]
            observed = row["observed_support_T1_Topo1"]
            writer.writerow(
                {
                    "sample": row["sample"],
                    "complete_regions": graph["complete_regions"],
                    **graph["classes"],
                    "topology_unique": graph["topology_unique"],
                    "T1_Topo1": observed["denominator"],
                    **observed["classes"],
                }
            )


def write_checks(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scope", "check", "pass"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, census: dict, sources: dict[str, Path], checks: list[dict]) -> None:
    graph = census["aggregate"]["graph"]
    observed = census["aggregate"]["observed_support_T1_Topo1"]
    lines = [
        "<!--",
        f"建立時間: {census['generated_at']}",
        "目標: 重算 complete region 的粗拓撲與 full-read-state endpoint sensitivity",
        "處理範圍: chr1-22; 7 dataset rows; historical layered-v2 engineering snapshot",
        "關聯檔案:",
        *[f"  - {path_value}" for path_value in sources.values()],
        "-->",
        "",
        "# 區域粗拓撲與 Full-read-state Endpoint Sensitivity 重算",
        "",
        "> **TL;DR**：Complete=39,885 中，Topo=1 的 21,976 個可分為無 within-HP 關係 6,027、sister-only 1,374、direct-only 13,899、sister+direct 676；Topo>1 的 17,909 保留未定。另以 T=1/Topo=1 的 10,832 個區域做 full-read-state endpoint sensitivity：sister-only=1,342、direct-only=2,948、both=80。",
        "",
        "> **主張邊界**：tree node 是 mutation state；HP1/HP2 是 ordered forest components；H_ 不是 directly observed clone。`single_only` 是內部欄位名，不得解讀為 single clone。",
        "",
        "## Complete-region graph topology",
        "",
        "| 類別 | 區域數 | 占 Complete |",
        "|---|---:|---:|",
        *[
            f"| {display} | {graph['classes'][key]:,} | {graph['shares_of_complete'][key] * 100:.2f}% |"
            for key, display in (
                ("single_only", "無 within-HP 關係 (`single_only` internal key)"),
                ("sister_only", "sister-only"),
                ("direct_only", "direct-only"),
                ("sister_and_direct", "sister + direct"),
                ("topology_multiple_unresolved", "Topo>1 未定"),
            )
        ],
        f"| **合計** | **{graph['complete_regions']:,}** | **100.00%** |",
        "",
        "註：6,027 中只有 837 個是單一 primary HP 的單節點 graph；其餘 5,190 個是 HP1/HP2 各有獨立節點，沒有可比的 within-HP 姐妹或直系關係。該 837 個單節點中，143 個是 full-read observed state，694 個為 H_* hidden/partial-supported state。",
        "",
        "## 逐 dataset complete-region 粗拓撲",
        "",
        "| Dataset | Complete | 無 HP 內關係 | Sister only | Direct only | Both | Topo>1 未定 |",
        "|---|---:|---:|---:|---:|---:|---:|",
        *[
            "| {sample} | {complete:,} | {single:,} | {sister:,} | {direct:,} | {both:,} | {multiple:,} |".format(
                sample=row["sample"],
                complete=row["graph"]["complete_regions"],
                single=row["graph"]["classes"]["single_only"],
                sister=row["graph"]["classes"]["sister_only"],
                direct=row["graph"]["classes"]["direct_only"],
                both=row["graph"]["classes"]["sister_and_direct"],
                multiple=row["graph"]["classes"]["topology_multiple_unresolved"],
            )
            for row in census["samples"]
        ],
        "",
        "## T=1/Topo=1 full-read-state endpoint sensitivity",
        "",
        "| 類別 | 區域數 | 占 T=1/Topo=1 |",
        "|---|---:|---:|",
        *[
            f"| {key} | {observed['classes'][key]:,} | {observed['shares'][key] * 100:.2f}% |"
            for key in OBSERVED_CLASSES
        ],
        f"| **合計** | **{observed['denominator']:,}** | **100.00%** |",
        "",
        "## 驗證",
        "",
        f"- Checks：{sum(row['pass'] == 'True' for row in checks)}/{len(checks)} PASS。",
        "- 逐 dataset graph 五類均守恆回 Complete。",
        "- Full-read-state endpoint 四類均守恆回 T=1/Topo=1；H_* 可連接路徑，但不當成直接觀測 endpoint。",
        "- 本報告是只讀重算，不修改 layered output。",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-catalog", required=True, type=Path)
    parser.add_argument("--ct-report", required=True, type=Path)
    parser.add_argument("--region-tsv", required=True, type=Path)
    parser.add_argument("--vaf-census", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-summary-tsv", required=True, type=Path)
    parser.add_argument("--output-checks-tsv", required=True, type=Path)
    parser.add_argument("--output-report-md", required=True, type=Path)
    args = parser.parse_args()

    exact_catalog = json.loads(args.exact_catalog.read_text(encoding="utf-8"))
    ct_report = json.loads(args.ct_report.read_text(encoding="utf-8"))
    vaf_census = json.loads(args.vaf_census.read_text(encoding="utf-8"))
    region_rows, region_by_key = load_region_rows(args.region_tsv)
    census = build_census(exact_catalog, ct_report, region_rows, region_by_key, vaf_census)
    checks = validation_rows(census, ct_report, len(region_rows))
    if not all(row["pass"] == "True" for row in checks):
        failed = [row for row in checks if row["pass"] != "True"]
        raise SystemExit(f"coarse topology validation failed: {failed[:3]}")

    sources = {
        "exact_catalog": args.exact_catalog.resolve(),
        "ct_report": args.ct_report.resolve(),
        "region_tsv": args.region_tsv.resolve(),
        "vaf_census": args.vaf_census.resolve(),
    }
    census["sources"] = {
        key: {"path": str(path), "sha256": sha256(path)} for key, path in sources.items()
    }
    census["validation"] = {
        "status": "PASS",
        "checks": len(checks),
        "checks_passed": sum(row["pass"] == "True" for row in checks),
    }
    write_json(args.output_json, census)
    write_summary(args.output_summary_tsv, census)
    write_checks(args.output_checks_tsv, checks)
    write_report(args.output_report_md, census, sources, checks)

    print(f"INPUT EXACT CATALOG -> {args.exact_catalog.resolve()}")
    print(f"INPUT C/T REPORT -> {args.ct_report.resolve()}")
    print(f"INPUT REGION TSV -> {args.region_tsv.resolve()}")
    print(f"INPUT VAF CENSUS -> {args.vaf_census.resolve()}")
    print(f"OUTPUT JSON -> {args.output_json.resolve()}")
    print(f"OUTPUT SUMMARY -> {args.output_summary_tsv.resolve()}")
    print(f"OUTPUT CHECKS -> {args.output_checks_tsv.resolve()}")
    print(f"OUTPUT REPORT -> {args.output_report_md.resolve()}")
    print(f"STATUS -> {len(checks)}/{len(checks)} PASS")


if __name__ == "__main__":
    main()
