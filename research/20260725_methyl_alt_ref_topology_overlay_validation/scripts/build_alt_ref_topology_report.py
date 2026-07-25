#!/usr/bin/env python3
"""Build the ALT/REF methylation × latest exact-PS topology audit report.

This deterministic, read-only analysis joins four evidence layers:

1. v8 formal methyl-group × linked-partner REF/ALT associations;
2. focal tumor-ALT versus tumor-REF joint methyl controls;
3. threshold-3 exact-PS × HP strict W components and direct edges;
4. the 2026-07-24 C++ candidate-topology/read-AF output.
5. the 2026-07-25 exact factorization of every AF-global-best tree.

The generated report deliberately keeps the read populations separate.  It does
not promote a model-ranked candidate tree to a true cellular lineage.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUTPUT = Path("/big7_disk/liaoyoyo2001/big7_disk_output")

TOPIC = REPO / "research/20260725_methyl_alt_ref_topology_overlay_validation"
DATA_DIR = TOPIC / "data"
RESULTS_DIR = TOPIC / "results"

PAIR_ROOT = (
    OUTPUT
    / "synthesis/observation_workspaces"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
)
PAIR_TSV = PAIR_ROOT / "methyl_ssnv_pair_results.tsv.gz"

M1_ROOT = (
    OUTPUT
    / "synthesis/observation_workspaces"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
)
M1_TSV = M1_ROOT / "all_ssnv_site_results.tsv.gz"

REF_ROOT = (
    OUTPUT
    / "synthesis/observation_workspaces"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
)
REF_TSV = REF_ROOT / "all_ssnv_tumor_ref_control_site_results.tsv.gz"
REF_SUMMARY = REF_ROOT / "all_ssnv_tumor_ref_control_summary.json"

STRICT_ROOT = (
    OUTPUT
    / "synthesis/research_rounds"
    / "20260723_production_exact_ps_strict_read_linkage"
)

CPP_ROOT = (
    OUTPUT
    / "synthesis/research_rounds"
    / "20260724_exact_ps_cpp_topology_af_all_samples"
    / "all7_strict_guard1000_v1"
)
CPP_COHORT_RECEIPT = CPP_ROOT / "cohort_receipt.json"

FACTORIZATION_ROOT = (
    OUTPUT
    / "synthesis/observation_workspaces"
    / "20260725_exact_ps_candidate_factorization"
    / "all7_v2"
)
FACTORIZATION_RECEIPT = FACTORIZATION_ROOT / "receipt.json"
FACTORIZATION_SOURCE_CURRENT = (
    REPO
    / "research/20260724_exact_ps_cpp_topology_signature_census"
    / "cpp/exact_topology_candidate_factorization.cpp"
)

TOPOLOGY_AUDIT_0723 = (
    REPO
    / "research/20260723_production_exact_ps_strict_read_linkage"
    / "20260723_exactPS嚴格ReadLinkage全資料集報告_01"
    / "topology_status_audit/all7_topology_completion_audit.json"
)

DORADO_SAMPLE = "HCC1395_DORADO"
HCC_FOCAL = 127986757
HCC_PARTNER = 127981978

PAIR_OUT = DATA_DIR / "formal_pair_alt_ref_topology_join.tsv"
LANE_OUT = DATA_DIR / "strict_hp_lane_cpp_topology_join.tsv"
FOCAL_OUT = DATA_DIR / "focal_alt_ref_control_join.tsv"
CROSS_OUT = DATA_DIR / "hcc1395_crosssource_pair_substructure.tsv"
ARTIFACT_OUT = TOPIC / "artifact.json"
RECEIPT_OUT = RESULTS_DIR / "validation_receipt.json"
SQLITE_OUT = DATA_DIR / "report.sqlite"
QUERY_DIR = TOPIC / "queries"


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def as_int(value: Any, default: int | None = None) -> int | None:
    if value in (None, ""):
        return default
    return int(value)


def as_float(value: Any, default: float | None = None) -> float | None:
    if value in (None, ""):
        return default
    result = float(value)
    return result if math.isfinite(result) else default


def json_value(value: str, fallback: Any) -> Any:
    if not value:
        return fallback
    return json.loads(value)


def read_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def sqlite_value(value: Any) -> Any:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return value


def write_sqlite_snapshot(datasets: Mapping[str, list[dict[str, Any]]]) -> dict[str, str]:
    """Materialize every report dataset and execute its reader-facing SELECT."""

    SQLITE_OUT.parent.mkdir(parents=True, exist_ok=True)
    QUERY_DIR.mkdir(parents=True, exist_ok=True)
    queries: dict[str, str] = {}
    with sqlite3.connect(SQLITE_OUT) as connection:
        for dataset, rows in datasets.items():
            if not rows:
                raise ValueError(f"dataset {dataset} must not be empty")
            fields = list(rows[0])
            for row in rows:
                if set(row) != set(fields):
                    raise ValueError(f"dataset {dataset} has inconsistent row fields")
            definitions = []
            for field in fields:
                values = [row[field] for row in rows if row[field] is not None]
                if values and all(isinstance(value, (bool, int)) for value in values):
                    sql_type = "INTEGER"
                elif values and all(isinstance(value, (bool, int, float)) for value in values):
                    sql_type = "REAL"
                else:
                    sql_type = "TEXT"
                definitions.append(f'"{field}" {sql_type}')
            connection.execute(f'DROP TABLE IF EXISTS "{dataset}"')
            connection.execute(f'CREATE TABLE "{dataset}" ({", ".join(definitions)})')
            placeholders = ", ".join("?" for _ in fields)
            quoted_fields = ", ".join(f'"{field}"' for field in fields)
            connection.executemany(
                f'INSERT INTO "{dataset}" ({quoted_fields}) VALUES ({placeholders})',
                [
                    tuple(sqlite_value(row[field]) for field in fields)
                    for row in rows
                ],
            )
            query = f'SELECT * FROM "{dataset}" ORDER BY rowid;'
            observed = connection.execute(query).fetchall()
            if len(observed) != len(rows):
                raise AssertionError(f"SQLite snapshot row count mismatch for {dataset}")
            query_path = QUERY_DIR / f"{dataset}.sql"
            query_path.write_text(query + "\n", encoding="utf-8")
            queries[dataset] = query
        connection.commit()
    return queries


def wilson(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return (math.nan, math.nan)
    p = success / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denominator
    return center - half, center + half


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Exact two-sided Fisher p using probability <= observed probability."""

    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2
    denominator = math.comb(total, row1)

    def probability(x: int) -> Fraction:
        return Fraction(math.comb(col1, x) * math.comb(total - col1, row1 - x), denominator)

    lower = max(0, row1 - (total - col1))
    upper = min(row1, col1)
    observed = probability(a)
    return float(sum((probability(x) for x in range(lower, upper + 1) if probability(x) <= observed), Fraction()))


def strict_chrom_dir(sample: str, chrom: str) -> Path:
    if sample == "HCC1395":
        return STRICT_ROOT / "hcc1395_strict_regions_v2/chromosomes" / chrom
    return (
        STRICT_ROOT
        / "all7_production_v1/samples"
        / sample
        / "strict_regions_v1/chromosomes"
        / chrom
    )


def strict_files(sample: str, chrom: str) -> tuple[Path, Path, Path]:
    root = strict_chrom_dir(sample, chrom)
    prefix = f"{sample}.{chrom}"
    return (
        root / f"{prefix}.site_component_membership.tsv.gz",
        root / f"{prefix}.components.tsv.gz",
        root / f"{prefix}.endpoint_edges.tsv.gz",
    )


def load_formal_pairs() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in read_tsv_gz(PAIR_TSV):
        if not parse_bool(row["endpoint_a_formal_pair_by_confirmed"]):
            continue
        groups = json_value(row["endpoint_a_groups"], [])
        table = json_value(row["endpoint_a_table"], [])
        if len(groups) != len(table):
            raise ValueError("endpoint_a_groups/table length mismatch")
        alt_fractions = [
            (counts[1] / sum(counts)) if sum(counts) else math.nan for counts in table
        ]
        enriched_index = max(range(len(groups)), key=lambda index: alt_fractions[index])
        rows.append(
            {
                **row,
                "focal_pos_i": int(row["focal_pos"]),
                "partner_pos_i": int(row["partner_pos"]),
                "groups_json": groups,
                "table_json": table,
                "enriched_group": groups[enriched_index],
                "enriched_alt_fraction": alt_fractions[enriched_index],
                "locus": (
                    f"{row['sample']} {row['focal_chrom']}:{int(row['focal_pos']):,}"
                    f"→{int(row['partner_pos']):,}"
                ),
                "locus_short": f"{row['sample']} {row['focal_chrom']}:{row['focal_pos']}",
            }
        )
    rows.sort(key=lambda row: (row["sample"], row["focal_chrom"], row["focal_pos_i"]))
    if len(rows) != 7:
        raise AssertionError(f"expected 7 formal pairs, observed {len(rows)}")
    return rows


def load_focal_controls(pairs: list[dict[str, Any]]) -> dict[tuple[str, str, int], dict[str, str]]:
    keys = {(row["sample"], row["focal_chrom"], row["focal_pos_i"]) for row in pairs}
    controls: dict[tuple[str, str, int], dict[str, str]] = {}
    for row in read_tsv_gz(REF_TSV):
        key = (row["dataset"], row["chrom"], int(row["pos"]))
        if key in keys:
            controls[key] = row
    if set(controls) != keys:
        raise AssertionError(f"missing focal controls: {sorted(keys - set(controls))}")
    return controls


def load_m1_rows(keys: set[tuple[str, str, int]]) -> dict[tuple[str, str, int], dict[str, str]]:
    rows: dict[tuple[str, str, int], dict[str, str]] = {}
    for row in read_tsv_gz(M1_TSV):
        key = (row["dataset"], row["chrom"], int(row["pos"]))
        if key in keys:
            rows[key] = row
    return rows


def focal_control_status(row: Mapping[str, str]) -> str:
    if parse_bool(row["joint_allele_axis_aligned"]):
        return "ALIGNED"
    if parse_bool(row["joint_allele_testable"]):
        v = as_float(row["joint_allele_v"], 0.0) or 0.0
        p = as_float(row["joint_allele_p_perm"], 1.0) or 1.0
        return "BELOW_EFFECT" if v < 0.30 else ("NOT_SIGNIFICANT" if p >= 0.05 else "NOT_ALIGNED")
    reason = row["joint_allele_not_testable_reason"]
    if "fewer_than_MIN_SIZE_tumor_REF" in reason:
        return "INSUFFICIENT_REF"
    if "joint_clustering_not_stable" in reason:
        return "UNSTABLE"
    return "NOT_TESTABLE"


def convert_edge_states(edge: Mapping[str, str], focal_pos: int) -> dict[str, int]:
    raw = {
        "RR": int(edge["support_RR"]),
        "RA": int(edge["support_RA"]),
        "AR": int(edge["support_AR"]),
        "AA": int(edge["support_AA"]),
    }
    if int(edge["pos_i1"]) == focal_pos:
        return raw
    if int(edge["pos_j1"]) == focal_pos:
        return {"RR": raw["RR"], "RA": raw["AR"], "AR": raw["RA"], "AA": raw["AA"]}
    raise ValueError("focal position is not an edge endpoint")


def load_cpp_rows(sample: str) -> list[dict[str, Any]]:
    path = CPP_ROOT / "samples" / sample / f"{sample}.topology.jsonl"
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            rows.append(json.loads(line))
    return rows


def load_factorization_rows(
    target_regions: Mapping[str, set[str]],
) -> dict[tuple[str, str], dict[str, Any]]:
    """Load only factorization rows needed by the formal-pair topology lanes."""

    rows: dict[tuple[str, str], dict[str, Any]] = {}
    for sample, region_ids in target_regions.items():
        if not region_ids:
            continue
        path = FACTORIZATION_ROOT / f"{sample}.candidate_factorization.jsonl"
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                region_id = row["region_id"]
                if region_id in region_ids:
                    rows[(sample, region_id)] = row
                    if len(
                        {
                            current_region
                            for current_sample, current_region in rows
                            if current_sample == sample
                        }
                    ) == len(region_ids):
                        break
    return rows


def topology_pair_status(row: Mapping[str, Any], focal: int, partner: int) -> str:
    if row.get("family_status") == "ABSTAIN_RESOURCE_LIMIT":
        return "ABSTAIN_RESOURCE_LIMIT"
    active = set(row.get("active_positions") or [])
    if not ({focal, partner} & active):
        return "PAIR_RR_BACKGROUND"
    if row.get("read_af_status") != "ranked_complete":
        return str(row.get("read_af_status") or row.get("unit_status"))
    return "UNIQUE_AF_CANDIDATE" if row.get("best_tree_unique") is True else "AF_TIED"


def topology_shape(row: Mapping[str, Any]) -> str:
    morphology = row.get("representative_best_morphology") or {}
    if not morphology:
        return "NA"
    if row.get("best_tree_unique") is not True:
        return "UNRESOLVED_REPRESENTATIVE_ONLY"
    sister = morphology.get("has_sister") is True
    direct = morphology.get("has_direct") is True
    if sister and direct:
        return "SISTER+DIRECT"
    if sister:
        return "SISTER_ONLY"
    if direct:
        return "DIRECT_ONLY"
    return "SINGLE_ONLY"


def candidate_pair_order(row: Mapping[str, Any], focal: int, partner: int) -> str:
    if row.get("best_tree_unique") is not True:
        return "UNRESOLVED"
    positions = {int(edge["acquired_position"]): index for index, edge in enumerate(row.get("representative_best_edges") or [])}
    if focal not in positions or partner not in positions:
        return "PAIR_NOT_ACTIVE_IN_CANDIDATE"
    if positions[focal] < positions[partner]:
        return "FOCAL_BEFORE_PARTNER"
    if positions[partner] < positions[focal]:
        return "PARTNER_BEFORE_FOCAL"
    return "UNRESOLVED"


def acquired_active_bit(parent: int, child: int) -> int:
    delta = parent ^ child
    if delta <= 0 or delta & (delta - 1):
        raise ValueError(f"edge {parent}->{child} is not a one-bit acquisition")
    return delta.bit_length() - 1


def factorized_pair_relation(
    row: Mapping[str, Any],
    focal: int,
    partner: int,
) -> dict[str, Any]:
    """Classify focal/partner relation across every AF-global-best tree.

    Each global-best tree acquires every active mutation exactly once.  An edge
    acquiring the partner from a parent vertex that already contains the focal
    is focal-before-partner.  The reverse test gives partner-before-focal; if
    neither mutation is present in the other's acquisition parent, they are
    incomparable branches.
    """

    positions = [int(position) for position in row["active_positions"]]
    if focal not in positions or partner not in positions:
        return {
            "status": "PAIR_NOT_ACTIVE_IN_GLOBAL_BEST",
            "order": "PAIR_NOT_ACTIVE_IN_CANDIDATE",
            "best_trees": int(row["best_tree_tie_count"]),
            "focal_before": 0,
            "partner_before": 0,
            "incomparable": 0,
        }
    focal_bit = positions.index(focal)
    partner_bit = positions.index(partner)
    best_trees = int(row["best_tree_tie_count"])
    focal_before = 0
    partner_before = 0
    partner_acquisition_total = 0
    focal_acquisition_total = 0
    for parent_raw, child_raw, count_raw in row["global_best_edge_incidence"]:
        parent = int(parent_raw)
        child = int(child_raw)
        count = int(count_raw)
        acquired_bit = acquired_active_bit(parent, child)
        if acquired_bit == partner_bit:
            partner_acquisition_total += count
            if parent & (1 << focal_bit):
                focal_before += count
        if acquired_bit == focal_bit:
            focal_acquisition_total += count
            if parent & (1 << partner_bit):
                partner_before += count
    if partner_acquisition_total != best_trees or focal_acquisition_total != best_trees:
        raise AssertionError(
            "global-best edge incidence does not acquire each active pair endpoint once"
        )
    incomparable = best_trees - focal_before - partner_before
    if incomparable < 0:
        raise AssertionError("pair relation counts exceed global-best tree count")
    if focal_before == best_trees:
        status = "RESOLVED_COMMON_ACROSS_BEST_TREES"
        order = "FOCAL_BEFORE_PARTNER"
    elif partner_before == best_trees:
        status = "RESOLVED_COMMON_ACROSS_BEST_TREES"
        order = "PARTNER_BEFORE_FOCAL"
    elif incomparable == best_trees:
        status = "RESOLVED_COMMON_ACROSS_BEST_TREES"
        order = "INCOMPARABLE_BRANCHES"
    else:
        status = "UNRESOLVED_ACROSS_BEST_TREES"
        order = "UNRESOLVED"
    return {
        "status": status,
        "order": order,
        "best_trees": best_trees,
        "focal_before": focal_before,
        "partner_before": partner_before,
        "incomparable": incomparable,
    }


def attach_factorized_relations(lanes: list[dict[str, Any]]) -> None:
    target_regions: dict[str, set[str]] = {}
    for lane in lanes:
        if lane["cpp_read_af_status"] != "ranked_complete":
            continue
        target_regions.setdefault(lane["sample"], set()).add(lane["cpp_region_id"])
    factor_rows = load_factorization_rows(target_regions)
    for lane in lanes:
        key = (lane["sample"], lane["cpp_region_id"])
        factor = factor_rows.get(key)
        if factor is None:
            lane.update(
                {
                    "factorization_resolution_class": "NOT_AVAILABLE_OR_NOT_RANKED",
                    "cpp_pair_relation_status": "NOT_RESOLVED",
                    "cpp_pair_order_across_best": "UNRESOLVED",
                    "cpp_pair_order_focal_before_count": 0,
                    "cpp_pair_order_partner_before_count": 0,
                    "cpp_pair_order_incomparable_count": 0,
                    "factorization_global_best_signature_count": None,
                    "factorization_coarse_classes": "NA",
                }
            )
            continue
        relation = factorized_pair_relation(
            factor,
            int(lane["focal_pos"]),
            int(lane["partner_pos"]),
        )
        lane.update(
            {
                "factorization_resolution_class": factor["resolution_class"],
                "cpp_pair_relation_status": relation["status"],
                "cpp_pair_order_across_best": relation["order"],
                "cpp_pair_order_focal_before_count": relation["focal_before"],
                "cpp_pair_order_partner_before_count": relation["partner_before"],
                "cpp_pair_order_incomparable_count": relation["incomparable"],
                "factorization_global_best_signature_count": int(
                    factor["global_best_signature_count"]
                ),
                "factorization_coarse_classes": ", ".join(
                    f"{name}:{count}"
                    for name, count in sorted(
                        factor["global_best_coarse_class_tree_counts"].items()
                    )
                ),
            }
        )


def build_strict_and_cpp(
    pairs: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[tuple[str, str, int, int], list[dict[str, Any]]]]:
    lanes: list[dict[str, Any]] = []
    by_pair: dict[tuple[str, str, int, int], list[dict[str, Any]]] = {}
    cpp_cache = {sample: load_cpp_rows(sample) for sample in {row["sample"] for row in pairs}}

    for pair in pairs:
        sample = pair["sample"]
        chrom = pair["focal_chrom"]
        focal = pair["focal_pos_i"]
        partner = pair["partner_pos_i"]
        membership_path, components_path, edges_path = strict_files(sample, chrom)

        memberships = [
            row
            for row in read_tsv_gz(membership_path)
            if row["threshold"] == "3" and int(row["pos1"]) in {focal, partner}
        ]
        components = {
            (row["linkage_basis"], row["phase_set"], row["component_id"]): row
            for row in read_tsv_gz(components_path)
            if row["threshold"] == "3"
        }
        edges = [
            row
            for row in read_tsv_gz(edges_path)
            if parse_bool(row["passes_primary_threshold"])
            and {int(row["pos_i1"]), int(row["pos_j1"])} == {focal, partner}
        ]

        focal_memberships = [row for row in memberships if int(row["pos1"]) == focal]
        partner_memberships = [row for row in memberships if int(row["pos1"]) == partner]
        matches = []
        for left in focal_memberships:
            for right in partner_memberships:
                identity = (left["linkage_basis"], left["phase_set"], left["component_id"])
                if identity == (right["linkage_basis"], right["phase_set"], right["component_id"]):
                    matches.append((left, components[identity]))
        if not matches:
            raise AssertionError(f"no SAME_W match for {pair['locus']}")

        pair_key = (sample, chrom, focal, partner)
        pair_lanes: list[dict[str, Any]] = []
        for membership, component in matches:
            edge_candidates = [
                edge
                for edge in edges
                if edge["linkage_basis"] == membership["linkage_basis"]
                and edge["phase_set"] == membership["phase_set"]
            ]
            if len(edge_candidates) != 1:
                raise AssertionError(f"expected one direct edge for {pair['locus']} lane")
            edge = edge_candidates[0]
            states = convert_edge_states(edge, focal)

            topology_candidates = []
            for topology in cpp_cache[sample]:
                coverage_positions = {int(item["position"]) for item in topology.get("af_coverage") or []}
                if (
                    topology["chrom"] == chrom
                    and str(topology["phase_set"]) == membership["phase_set"]
                    and str(topology["hp_family"]) == edge["hp_family"]
                    and {focal, partner}.issubset(coverage_positions)
                ):
                    topology_candidates.append(topology)
            if len(topology_candidates) != 1:
                raise AssertionError(
                    f"expected one C++ topology block for {pair['locus']} {membership['linkage_basis']}, "
                    f"observed {len(topology_candidates)}"
                )
            topology = topology_candidates[0]
            coverage = {int(item["position"]): item for item in topology["af_coverage"]}
            lane = {
                "sample": sample,
                "chrom": chrom,
                "focal_pos": focal,
                "partner_pos": partner,
                "pair_label": pair["locus_short"],
                "linkage_basis": membership["linkage_basis"],
                "hp_family": edge["hp_family"],
                "phase_set": membership["phase_set"],
                "component_id": membership["component_id"],
                "W_k": int(component["k"]),
                "W_start": int(component["start1"]),
                "W_end": int(component["end1"]),
                "direct_support": int(edge["support_total"]),
                "state_RR": states["RR"],
                "state_RA": states["RA"],
                "state_AR": states["AR"],
                "state_AA": states["AA"],
                "lane_role": (
                    "FORMAL_SIGNAL"
                    if states["RA"] + states["AR"] + states["AA"] > 0
                    else "PAIRED_BACKGROUND"
                ),
                "cpp_region_id": topology["region_id"],
                "cpp_block_k": int(topology["original_bit_count"]),
                "cpp_active_k": int(topology["active_bit_count"]),
                "cpp_family_status": topology["family_status"],
                "cpp_unit_status": topology["unit_status"],
                "cpp_read_af_status": topology["read_af_status"],
                "cpp_pair_status": topology_pair_status(topology, focal, partner),
                "cpp_minimum_vertex_sets": as_int(topology.get("minimum_vertex_set_count")),
                "cpp_total_trees": as_int(topology.get("total_tree_count")),
                "cpp_best_tree_ties": as_int(topology.get("best_tree_tie_count")),
                "cpp_best_tree_unique": topology.get("best_tree_unique"),
                "cpp_best_score": topology.get("best_score_fraction"),
                "cpp_shape": topology_shape(topology),
                "cpp_pair_order": candidate_pair_order(topology, focal, partner),
                "cpp_reason": topology.get("reason"),
                "cpp_focal_ref": int(coverage[focal]["ref_reads"]),
                "cpp_focal_alt": int(coverage[focal]["alt_reads"]),
                "cpp_partner_ref": int(coverage[partner]["ref_reads"]),
                "cpp_partner_alt": int(coverage[partner]["alt_reads"]),
                "cpp_input_vaf_eligible": topology.get("input_vaf_eligible"),
                "cpp_edges_json": json.dumps(topology.get("representative_best_edges") or [], separators=(",", ":")),
                "cpp_vertices_json": json.dumps(topology.get("representative_best_vertices") or [], separators=(",", ":")),
            }
            lanes.append(lane)
            pair_lanes.append(lane)
        by_pair[pair_key] = pair_lanes

    if len(lanes) != 10:
        raise AssertionError(f"expected 10 shared W lanes, observed {len(lanes)}")
    if sum(row["direct_support"] for row in lanes) != 790:
        raise AssertionError("direct support total drifted from 790")
    attach_factorized_relations(lanes)
    return lanes, by_pair


def build_dorado_crosssource(
    hcc_pair: dict[str, Any],
    hcc_lanes: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    chrom = hcc_pair["focal_chrom"]
    focal = hcc_pair["focal_pos_i"]
    partner = hcc_pair["partner_pos_i"]
    membership_path, components_path, edges_path = strict_files(DORADO_SAMPLE, chrom)

    memberships = [
        row
        for row in read_tsv_gz(membership_path)
        if row["threshold"] == "3" and int(row["pos1"]) in {focal, partner}
    ]
    shared = {
        (row["linkage_basis"], row["phase_set"], row["component_id"])
        for row in memberships
        if int(row["pos1"]) == focal
    } & {
        (row["linkage_basis"], row["phase_set"], row["component_id"])
        for row in memberships
        if int(row["pos1"]) == partner
    }
    if len(shared) != 1:
        raise AssertionError("DORADO HCC pair must map to exactly one shared W")
    identity = next(iter(shared))
    component = next(
        row
        for row in read_tsv_gz(components_path)
        if row["threshold"] == "3"
        and (row["linkage_basis"], row["phase_set"], row["component_id"]) == identity
    )
    edge = next(
        row
        for row in read_tsv_gz(edges_path)
        if parse_bool(row["passes_primary_threshold"])
        and row["linkage_basis"] == identity[0]
        and row["phase_set"] == identity[1]
        and {int(row["pos_i1"]), int(row["pos_j1"])} == {focal, partner}
    )
    states = convert_edge_states(edge, focal)
    topology = next(
        row
        for row in load_cpp_rows(DORADO_SAMPLE)
        if row["chrom"] == chrom
        and str(row["phase_set"]) == identity[1]
        and str(row["hp_family"]) == edge["hp_family"]
        and {focal, partner}.issubset(
            {int(item["position"]) for item in row.get("af_coverage") or []}
        )
    )
    coverage = {int(item["position"]): item for item in topology["af_coverage"]}

    hcc_lane = next(row for row in hcc_lanes if row["hp_family"] == "2")
    cross_rows = []
    for sample, lane, topology_row, coverage_row, W_k in [
        (
            "HCC1395",
            hcc_lane,
            None,
            {
                focal: {"ref_reads": hcc_lane["cpp_focal_ref"], "alt_reads": hcc_lane["cpp_focal_alt"]},
                partner: {"ref_reads": hcc_lane["cpp_partner_ref"], "alt_reads": hcc_lane["cpp_partner_alt"]},
            },
            hcc_lane["W_k"],
        ),
        (
            DORADO_SAMPLE,
            {
                "phase_set": identity[1],
                "hp_family": edge["hp_family"],
                "direct_support": int(edge["support_total"]),
                "state_RR": states["RR"],
                "state_RA": states["RA"],
                "state_AR": states["AR"],
                "state_AA": states["AA"],
                "cpp_pair_status": topology_pair_status(topology, focal, partner),
                "cpp_pair_order": candidate_pair_order(topology, focal, partner),
                "cpp_shape": topology_shape(topology),
                "cpp_minimum_vertex_sets": as_int(topology.get("minimum_vertex_set_count")),
                "cpp_best_tree_ties": as_int(topology.get("best_tree_tie_count")),
            },
            topology,
            coverage,
            int(component["k"]),
        ),
    ]:
        focal_alt_pair = int(lane["state_AR"]) + int(lane["state_AA"])
        partner_alt_in_focal_alt = int(lane["state_AA"])
        low, high = wilson(partner_alt_in_focal_alt, focal_alt_pair)
        cross_rows.append(
            {
                "sample": sample,
                "biological_id": "HCC1395",
                "technology_role": "canonical" if sample == "HCC1395" else "DORADO technical source",
                "phase_set": lane["phase_set"],
                "hp_family": lane["hp_family"],
                "W_k": W_k,
                "direct_support": lane["direct_support"],
                "state_RR": lane["state_RR"],
                "state_RA": lane["state_RA"],
                "state_AR": lane["state_AR"],
                "state_AA": lane["state_AA"],
                "partner_alt_given_focal_alt": partner_alt_in_focal_alt / focal_alt_pair,
                "partner_alt_ci_low": low,
                "partner_alt_ci_high": high,
                "cpp_focal_af": coverage_row[focal]["alt_reads"]
                / (coverage_row[focal]["ref_reads"] + coverage_row[focal]["alt_reads"]),
                "cpp_partner_af": coverage_row[partner]["alt_reads"]
                / (coverage_row[partner]["ref_reads"] + coverage_row[partner]["alt_reads"]),
                "cpp_pair_status": lane["cpp_pair_status"],
                "cpp_pair_order": lane["cpp_pair_order"],
                "cpp_shape": lane["cpp_shape"],
                "cpp_minimum_vertex_sets": lane["cpp_minimum_vertex_sets"],
                "cpp_best_tree_ties": lane["cpp_best_tree_ties"],
            }
        )

    hcc = cross_rows[0]
    dorado = cross_rows[1]
    p = fisher_two_sided(hcc["state_AR"], hcc["state_AA"], dorado["state_AR"], dorado["state_AA"])
    for row in cross_rows:
        row["conditional_fisher_p"] = p
        row["absolute_fraction_difference"] = abs(
            hcc["partner_alt_given_focal_alt"] - dorado["partner_alt_given_focal_alt"]
        )
    return cross_rows


def pair_topology_summary(lanes: list[dict[str, Any]]) -> tuple[str, dict[str, Any] | None]:
    signal_lanes = [
        row
        for row in lanes
        if row["state_AR"] + row["state_RA"] + row["state_AA"] > 0
    ]
    ranked = [row for row in signal_lanes if row["cpp_pair_status"] in {"UNIQUE_AF_CANDIDATE", "AF_TIED"}]
    if any(row["cpp_pair_status"] == "UNIQUE_AF_CANDIDATE" for row in ranked):
        best = next(row for row in ranked if row["cpp_pair_status"] == "UNIQUE_AF_CANDIDATE")
        return "UNIQUE_AF_CANDIDATE", best
    if any(row["cpp_pair_status"] == "AF_TIED" for row in ranked):
        best = next(row for row in ranked if row["cpp_pair_status"] == "AF_TIED")
        return "AF_TIED", best
    if any(row["cpp_pair_status"] == "ABSTAIN_RESOURCE_LIMIT" for row in signal_lanes):
        return "ABSTAIN_RESOURCE_LIMIT", next(
            row for row in signal_lanes if row["cpp_pair_status"] == "ABSTAIN_RESOURCE_LIMIT"
        )
    return "NOT_EVALUABLE", None


def pair_relation_summary(lanes: list[dict[str, Any]]) -> tuple[str, str, dict[str, Any] | None]:
    resolved = [
        row
        for row in lanes
        if row["cpp_pair_relation_status"] == "RESOLVED_COMMON_ACROSS_BEST_TREES"
    ]
    if not resolved:
        return "UNRESOLVED", "UNRESOLVED", None
    orders = {row["cpp_pair_order_across_best"] for row in resolved}
    if len(orders) != 1:
        raise AssertionError("resolved HP lanes disagree on focal/partner relation")
    representative = resolved[0]
    return "RESOLVED_COMMON_ACROSS_BEST_TREES", next(iter(orders)), representative


def build_pair_rows(
    pairs: list[dict[str, Any]],
    controls: Mapping[tuple[str, str, int], dict[str, str]],
    lanes_by_pair: Mapping[tuple[str, str, int, int], list[dict[str, Any]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    pair_rows: list[dict[str, Any]] = []
    focal_rows: list[dict[str, Any]] = []
    for pair in pairs:
        key = (pair["sample"], pair["focal_chrom"], pair["focal_pos_i"])
        lane_key = (*key, pair["partner_pos_i"])
        control = controls[key]
        control_status = focal_control_status(control)
        topology_status, representative = pair_topology_summary(lanes_by_pair[lane_key])
        relation_status, relation_order, relation_representative = pair_relation_summary(
            lanes_by_pair[lane_key]
        )
        pair_rows.append(
            {
                "sample": pair["sample"],
                "focal": f"{pair['focal_chrom']}:{pair['focal_pos_i']} {pair['focal_ref']}>{pair['focal_alt']}",
                "partner": f"{pair['partner_chrom']}:{pair['partner_pos_i']} {pair['partner_ref']}>{pair['partner_alt']}",
                "locus_short": pair["locus_short"],
                "g1_informative_reads": int(pair["endpoint_a_n_informative"]),
                "g1_groups": ", ".join(pair["groups_json"]),
                "g1_table": pair["endpoint_a_table"],
                "g1_cramers_v": float(pair["endpoint_a_cramers_v"]),
                "g1_delta_alt_fraction": float(pair["endpoint_a_delta_alt_fraction"]),
                "g1_q_by": float(pair["endpoint_a_q_global_by"]),
                "g1_conditional_p": float(pair["endpoint_a_p_conditional_perm"]),
                "g1_enriched_group": pair["enriched_group"],
                "same_W": True,
                "W_containers": len(lanes_by_pair[lane_key]),
                "direct_support_total": sum(row["direct_support"] for row in lanes_by_pair[lane_key]),
                "focal_ref_alt_status": control_status,
                "focal_alt_reads": int(control["n_tumor_alt"]),
                "focal_ref_reads": int(control["n_tumor_ref"]),
                "focal_joint_v": as_float(control["joint_allele_v"]),
                "focal_joint_p": as_float(control["joint_allele_p_perm"]),
                "candidate_T_status": topology_status,
                "candidate_T_minimum_trees": representative["cpp_minimum_vertex_sets"] if representative else None,
                "candidate_T_best_ties": representative["cpp_best_tree_ties"] if representative else None,
                "candidate_T_shape": (
                    representative["factorization_coarse_classes"]
                    if representative
                    and representative["factorization_coarse_classes"] != "NA"
                    else (representative["cpp_shape"] if representative else "NA")
                ),
                "candidate_topology_resolution": (
                    representative["factorization_resolution_class"]
                    if representative
                    else "NOT_AVAILABLE_OR_NOT_RANKED"
                ),
                "candidate_pair_relation_status": relation_status,
                "candidate_pair_order": relation_order,
                "candidate_pair_order_focal_before_count": (
                    relation_representative["cpp_pair_order_focal_before_count"]
                    if relation_representative
                    else 0
                ),
                "candidate_pair_order_partner_before_count": (
                    relation_representative["cpp_pair_order_partner_before_count"]
                    if relation_representative
                    else 0
                ),
                "candidate_pair_order_incomparable_count": (
                    relation_representative["cpp_pair_order_incomparable_count"]
                    if relation_representative
                    else 0
                ),
                "endpoint_b_status": pair["endpoint_b_status"],
            }
        )
        focal_rows.append(
            {
                "sample": pair["sample"],
                "locus": pair["locus_short"],
                "n_ALT": int(control["n_tumor_alt"]),
                "n_REF": int(control["n_tumor_ref"]),
                "joint_status": control["joint_status"],
                "joint_testable": parse_bool(control["joint_allele_testable"]),
                "joint_V": as_float(control["joint_allele_v"]),
                "joint_p_perm": as_float(control["joint_allele_p_perm"]),
                "axis_aligned": parse_bool(control["joint_allele_axis_aligned"]),
                "classification": control_status,
                "REF_only_evaluable": parse_bool(control["ref_evaluable"]),
                "REF_only_stable_multigroup": parse_bool(control["ref_stable_null_multigroup"]),
                "background_interpretation": control["background_control_interpretation"],
                "not_testable_reason": control["joint_allele_not_testable_reason"],
            }
        )
    return pair_rows, focal_rows


def source(
    source_id: str,
    label: str,
    path: str,
    description: str,
    tables: list[str],
    filters: list[str],
    definitions: list[str],
) -> dict[str, Any]:
    return {
        "id": source_id,
        "label": label,
        "path": path,
        "query": {
            "engine": "python-stdlib",
            "language": "python",
            "description": description,
            "tables_used": tables,
            "filters": filters,
            "metric_definitions": definitions,
            "executed_at": now_iso(),
        },
    }


def report_dataset_source(
    dataset: str,
    label: str,
    description: str,
    definitions: list[str],
) -> dict[str, Any]:
    fields = {
        "id": f"source_{dataset}",
        "label": label,
        "path": f"research/20260725_methyl_alt_ref_topology_overlay_validation/queries/{dataset}.sql",
        "query": {
            "engine": "sqlite",
            "language": "sql",
            "description": description,
            "sql": f'SELECT * FROM "{dataset}" ORDER BY rowid;',
            "tables_used": [dataset],
            "filters": ["bounded reviewed snapshot generated by the deterministic join script"],
            "metric_definitions": definitions,
            "executed_at": now_iso(),
        },
    }
    return fields


def make_artifact(
    pairs: list[dict[str, Any]],
    pair_rows: list[dict[str, Any]],
    focal_rows: list[dict[str, Any]],
    lanes: list[dict[str, Any]],
    cross_rows: list[dict[str, Any]],
    hcc_m1: Mapping[tuple[str, str, int], dict[str, str]],
) -> dict[str, Any]:
    generated_at = now_iso()
    topology_counts = Counter(row["candidate_T_status"] for row in pair_rows)
    focal_counts = Counter(row["classification"] for row in focal_rows)
    unique_count = topology_counts["UNIQUE_AF_CANDIDATE"]
    resolved_pair_count = sum(
        row["candidate_pair_relation_status"]
        == "RESOLVED_COMMON_ACROSS_BEST_TREES"
        for row in pair_rows
    )
    aligned_count = focal_counts["ALIGNED"]
    signal_support = sum(
        row["direct_support"]
        for row in lanes
        if row["lane_role"] == "FORMAL_SIGNAL"
    )
    background_support = sum(
        row["direct_support"]
        for row in lanes
        if row["lane_role"] == "PAIRED_BACKGROUND"
    )
    if (signal_support, background_support) != (638, 152):
        raise AssertionError(
            "signal/background direct-support split drifted: "
            f"{signal_support}/{background_support}"
        )

    hcc_pair = next(row for row in pairs if row["sample"] == "HCC1395")
    hcc_groups = []
    for group, counts in zip(hcc_pair["groups_json"], hcc_pair["table_json"]):
        hcc_groups.extend(
            [
                {"methyl_group": group, "partner_allele": "REF", "reads": counts[0]},
                {"methyl_group": group, "partner_allele": "ALT", "reads": counts[1]},
            ]
        )

    evidence_ladder = [
        {"stage": "Formal partner MG×REF/ALT", "pair_count": 7, "share": 1.0},
        {"stage": "Same strict W", "pair_count": 7, "share": 1.0},
        {"stage": "Threshold-qualified direct edge", "pair_count": 7, "share": 1.0},
        {"stage": "Focal REF/ALT joint testable", "pair_count": 3, "share": 3 / 7},
        {"stage": "Focal allele-axis aligned", "pair_count": aligned_count, "share": aligned_count / 7},
        {"stage": "Unique latest AF candidate T", "pair_count": unique_count, "share": unique_count / 7},
        {
            "stage": "Pair relation resolved across best T",
            "pair_count": resolved_pair_count,
            "share": resolved_pair_count / 7,
        },
    ]
    topology_status = [
        {"status": "Unique AF candidate", "pairs": topology_counts["UNIQUE_AF_CANDIDATE"]},
        {"status": "AF tied", "pairs": topology_counts["AF_TIED"]},
        {"status": "Resource abstain", "pairs": topology_counts["ABSTAIN_RESOURCE_LIMIT"]},
    ]
    focal_status = [
        {"status": "Aligned", "sites": focal_counts["ALIGNED"]},
        {
            "status": "Below effect / non-significant",
            "sites": focal_counts["BELOW_EFFECT"] + focal_counts["NOT_SIGNIFICANT"],
        },
        {"status": "Unstable", "sites": focal_counts["UNSTABLE"]},
        {"status": "Insufficient REF", "sites": focal_counts["INSUFFICIENT_REF"]},
    ]
    partner_effect = [
        {
            "locus": row["locus_short"].replace("chr", "c"),
            "cramers_v": row["g1_cramers_v"],
            "delta_alt_fraction": row["g1_delta_alt_fraction"],
            "informative_reads": row["g1_informative_reads"],
        }
        for row in pair_rows
    ]
    definitions = [
        {
            "layer": "Focal M1",
            "read_population": "focal-ALT methyl-matrix core reads",
            "question": "ALT-bearing reads內是否有穩定甲基多群？",
            "tree_role": "候選分群；不是 clone node",
        },
        {
            "layer": "G1 partner association",
            "read_population": "focal-ALT core且partner allele可判定",
            "question": "甲基群與partner REF/ALT是否共分離？",
            "tree_role": "pairwise state projection overlay",
        },
        {
            "layer": "Tumor focal REF/ALT control",
            "read_population": "ALT+REF重新joint clustering",
            "question": "focal allele本身是否與甲基partition對齊？",
            "tree_role": "獨立徽章；群標不可與M1一對一",
        },
        {
            "layer": "Strict W / endpoint edge",
            "read_population": "exact-PS×HP fixed R/A canonical molecules",
            "question": "兩位點是否同W且有直接read linkage？",
            "tree_role": "無向結構骨架",
        },
        {
            "layer": "C++ candidate T / read-AF",
            "read_population": "whole-block patterns且pattern weight≥3",
            "question": "模型內最低樹族可否完整枚舉並由read-AF排序？",
            "tree_role": "model-conditional candidate；非獨立VAF驗證",
        },
    ]

    hcc_canonical = cross_rows[0]
    hcc_dorado = cross_rows[1]
    hcc_tree_text = (
        "ROOT  RRRRR  （模型根；不等於已觀察到的REF clone）\n"
        "└─ chr3:127,986,757 G>A  focal\n"
        "   └─ RRRRA\n"
        "      ├─ chr3:127,962,164 → RRARA\n"
        "      └─ chr3:127,981,978 C>G  partner → RRRAA\n"
        "         └─ chr3:127,948,631 → RARAA\n"
        "            └─ chr3:127,912,458 → AARAA\n\n"
        "Pairwise methyl overlay（只投影 focal / partner，不是完整5-bit node assignment）\n"
        "partner REF：MG 1-1 = 21，MG 1-2 = 0\n"
        "partner ALT：MG 1-1 = 11，MG 1-2 = 19"
    )

    dataset_source_specs = {
        "headline": (
            "Executed SQLite snapshot query: headline",
            "Reviewed headline metrics from the formal-pair, focal-control, strict-W, C++ topology, and factorization left join.",
        ),
        "evidence_ladder": (
            "Executed SQLite snapshot query: evidence_ladder",
            "Pair counts retained at each evidence layer; denominator is seven selected formal G1 positives.",
        ),
        "topology_status": (
            "Executed SQLite snapshot query: topology_status",
            "Full candidate-tree status from the 2026-07-24 bounded C++ run.",
        ),
        "partner_effect": (
            "Executed SQLite snapshot query: partner_effect",
            "Formal v8 methyl-group × partner-allele effect sizes selected from 407,738 evaluated rows.",
        ),
        "pair_master": (
            "Executed SQLite snapshot query: pair_master",
            "Seven-row coordinate-preserving join across G1, focal ALT/REF control, strict W, topology, and factorization.",
        ),
        "hcc_methyl_counts": (
            "Executed SQLite snapshot query: hcc_methyl_counts",
            "HCC1395 focal-ALT methyl-core reads tabulated by methyl group and linked-partner allele.",
        ),
        "hcc_layers": (
            "Executed SQLite snapshot query: hcc_layers",
            "Reader-facing denominator audit for the four HCC1395 read populations.",
        ),
        "hcc_crosssource": (
            "Executed SQLite snapshot query: hcc_crosssource",
            "HCC1395 canonical versus DORADO strict pair states and model-conditional candidate relation.",
        ),
        "focal_status": (
            "Executed SQLite snapshot query: focal_status",
            "Classification counts for seven focal ALT+REF methyl joint controls.",
        ),
        "focal_controls": (
            "Executed SQLite snapshot query: focal_controls",
            "Seven focal tumor ALT+REF joint-clustering control rows; the current selected subset has no new global multiplicity correction.",
        ),
        "definitions": (
            "Executed SQLite snapshot query: definitions",
            "Formal read-population and claim-boundary definitions used in this report.",
        ),
        "strict_lanes": (
            "Executed SQLite snapshot query: strict_lanes",
            "Ten exact-PS×HP lanes with focal-first RR/RA/AR/AA states and candidate factorization fields.",
        ),
    }
    sources = [
        report_dataset_source(
            dataset,
            label,
            description,
            [
                "Every query reads the deterministic bounded SQLite snapshot generated by this report.",
                "Upstream source identities and SHA256 values are fixed in validation_receipt.json.",
                "C++ read-AF is supported-pattern ALT/(REF+ALT), not caller VAF, CCF, or an independent likelihood.",
            ],
        )
        for dataset, (label, description) in dataset_source_specs.items()
    ]

    pair_columns = [
        {"field": "sample", "label": "Dataset", "type": "text"},
        {"field": "focal", "label": "Focal", "type": "text"},
        {"field": "partner", "label": "Partner", "type": "text"},
        {"field": "g1_cramers_v", "label": "G1 V", "format": "number"},
        {"field": "g1_delta_alt_fraction", "label": "G1 Δ ALT fraction", "format": "number"},
        {"field": "W_containers", "label": "Shared W lanes", "format": "number"},
        {"field": "direct_support_total", "label": "Direct support", "format": "number"},
        {"field": "focal_ref_alt_status", "label": "Focal REF/ALT", "type": "text"},
        {"field": "candidate_T_status", "label": "Latest candidate T", "type": "text"},
        {"field": "candidate_T_best_ties", "label": "Best T ties", "format": "number"},
        {
            "field": "candidate_pair_relation_status",
            "label": "Pair relation",
            "type": "text",
        },
        {"field": "candidate_pair_order", "label": "Pair order", "type": "text"},
    ]

    artifact = {
        "surface": "report",
        "manifest": {
            "title": "ALT/REF 甲基差異 × 最新 exact-PS 候選拓撲對應驗證",
            "version": 1,
            "surface": "report",
            "description": (
                "Seven formal methyl–partner associations joined to focal ALT/REF controls, "
                "strict W/edges, the 2026-07-24 C++ candidate topology, and an HCC1395 "
                "cross-technical pairwise substructure comparison."
            ),
            "generatedAt": generated_at,
            "blocks": [
                {
                    "id": "title",
                    "type": "markdown",
                    "body": "# ALT/REF 甲基差異 × 最新 exact-PS 候選拓撲對應驗證",
                },
                {
                    "id": "scope_ribbon",
                    "type": "markdown",
                    "body": (
                        "**TASK B / COMPLETE JOIN, PARTIAL BIOLOGICAL VALIDATION — "
                        "7 個 formal G1 pair 已完整對接；7/24 C++ cohort 技術 PASS，"
                        "但 `validation_evidence_eligible=false`，不可解讀為真實 cellular lineage。**"
                    ),
                },
                {
                    "id": "technical_summary",
                    "type": "markdown",
                    "body": (
                        "## 技術結論：能標到候選樹，但必須把三種 ALT/REF 問題分開\n\n"
                        "- **Partner REF/ALT：7/7 formal G1 pairs** 都有甲基群共分離，且 7/7 同 strict W、"
                        "10/10 HP lanes 有 direct edge（signal 638 reads；配對 RR-only "
                        "background 152；合計 790）。\n"
                        "- **Focal REF/ALT：只有 1/7** 通過 joint allele-axis gate；HCC1395 為 ALT=98、REF=2，"
                        "完全不可評估 focal REF/ALT methyl difference。\n"
                        "- **最新 candidate T：2/7 unique、1/7 AF tied、4/7 resource abstain。** "
                        "**局部 pair relation 則有 3/7 已解**：兩個 unique T，加上 H2009 chr18 "
                        "雖有 6 棵並列 best trees，但 6/6 都是 focal→partner。Unique 只表示目前"
                        "模型與 supported-pattern read-AF 下有單一 winner。\n"
                        "- HCC1395 可把 partner-ALT 的 MG enrichment 疊在 candidate focal→partner branch 上；"
                        "這是 **pairwise projection association**，不是突變造成甲基改變的證明。"
                    ),
                },
                {
                    "id": "headline_metrics",
                    "type": "metric-strip",
                    "cardIds": [
                        "formal_pairs_card",
                        "strict_map_card",
                        "focal_control_card",
                        "candidate_tree_card",
                        "release_ceiling_card",
                    ],
                },
                {
                    "id": "evidence_heading",
                    "type": "markdown",
                    "body": (
                        "## Read-linkage 全部對上，但 focal REF/ALT 與唯一候選樹只涵蓋少數\n\n"
                        "下圖的分母固定為 7 個已篩出的 formal G1 positives；它不是 407,738 個候選 pair 的"
                        "方法準確率。Direct edge 是無向結構證據，candidate T 則再受完整枚舉與 read-AF "
                        "tie 影響。"
                    ),
                    "sourceId": "source_evidence_ladder",
                },
                {"id": "evidence_chart_block", "type": "chart", "chartId": "evidence_ladder_chart"},
                {
                    "id": "topology_status_heading",
                    "type": "markdown",
                    "body": (
                        "## 7/24 最新 C++ 對位：2 個 unique candidate、1 個 tied、4 個 fail-closed abstain\n\n"
                        "HCC1395 與 HCC1954 各有一個 unique read-AF winner；H2009 chr18 有 6 個並列 best "
                        "trees，但 7/25 exact factorization 顯示 6/6 都包含 focal-bearing state→partner "
                        "acquisition，故局部邊已解、整樹仍不唯一。其餘四個 pair 的 mutation-bearing lane "
                        "觸發 1,000 search-node guard。"
                    ),
                    "sourceId": "source_topology_status",
                },
                {"id": "topology_status_chart_block", "type": "chart", "chartId": "topology_status_chart"},
                {
                    "id": "partner_effect_heading",
                    "type": "markdown",
                    "body": (
                        "## 7 個 partner REF/ALT association 的效應都大，但屬於被選中的陽性集合\n\n"
                        "Cramér's V 為 0.614–0.964；所有 pair 已通過 global BY 與 conditional permutation "
                        "gate。這證明 methyl groups 與 partner allele 在同一 focal-ALT read subset 中共分離，"
                        "不等於 focal mutation 導致甲基變化。"
                    ),
                    "sourceId": "source_partner_effect",
                },
                {"id": "partner_effect_chart_block", "type": "chart", "chartId": "partner_effect_chart"},
                {
                    "id": "pair_table_block",
                    "type": "table",
                    "tableId": "pair_master_table",
                    "layout": "full",
                },
                {
                    "id": "hcc_heading",
                    "type": "markdown",
                    "body": (
                        "## HCC1395：可把 partner-ALT enrichment 疊在唯一候選 branch，但不是 node truth\n\n"
                        "HCC1395 HP2 的 k=5 unit 有 44 棵最低樹；supported-pattern read-AF 選出 1 棵 unique "
                        "winner，整體形狀為 Sister＋direct。Focal 先出現，partner 再由 focal-bearing state "
                        "分支。下方 MG×partner allele 圖只用 51 個 methyl-core、partner-callable reads；"
                        "不可當成完整 5-bit node membership。"
                    ),
                    "sourceId": "source_hcc_methyl_counts",
                },
                {"id": "hcc_methyl_chart_block", "type": "chart", "chartId": "hcc_methyl_chart"},
                {
                    "id": "hcc_tree",
                    "type": "markdown",
                    "body": "### 最新 HCC1395 candidate T 與 pairwise methyl overlay\n\n```\n"
                    + hcc_tree_text
                    + "\n```",
                },
                {
                    "id": "hcc_layers_table_block",
                    "type": "table",
                    "tableId": "hcc_layers_table",
                    "layout": "full",
                },
                {
                    "id": "crosssource_heading",
                    "type": "markdown",
                    "body": (
                        "## HCC1395_DORADO 重現 focal→partner candidate 子結構，但沒有重現甲基多群\n\n"
                        f"兩來源 strict pairwise state 在 focal-ALT reads 中的 partner-ALT 比例為 "
                        f"{hcc_canonical['partner_alt_given_focal_alt']:.1%} 與 "
                        f"{hcc_dorado['partner_alt_given_focal_alt']:.1%}（差 "
                        f"{hcc_canonical['absolute_fraction_difference']:.1%}；Fisher p="
                        f"{hcc_canonical['conditional_fisher_p']:.3f}）。兩者 candidate order 都是 "
                        "`focal before partner`，DORADO k=2 是 HCC k=5 tree 的 pairwise contraction。"
                        "但 DORADO focal M1 是穩定單群 K=1，formal methyl pair 狀態為 ONE_PLATFORM_ONLY；"
                        "因此只支持 genetic candidate-substructure 的技術再現，不支持 methyl association 再現。"
                    ),
                    "sourceId": "source_hcc_crosssource",
                },
                {"id": "crosssource_chart_block", "type": "chart", "chartId": "crosssource_chart"},
                {
                    "id": "crosssource_table_block",
                    "type": "table",
                    "tableId": "crosssource_table",
                    "layout": "full",
                },
                {
                    "id": "focal_control_heading",
                    "type": "markdown",
                    "body": (
                        "## Focal REF/ALT control：只有 H2009 chr4 通過，HCC1395 因 REF=2 不可評估\n\n"
                        "Joint control 會把 ALT+REF reads 重新 clustering，沒有沿用 ALT-only M1 centroid。"
                        "因此其群標不能直接接到原 MG 1-1/1-2。H2009 chr4 的 V=0.618、p=0.002、"
                        "axis-aligned=true 可加一個 focal allele-associated methyl badge；但該 pair 的最新"
                        "候選樹剛好 resource abstain，只能標在 W，不能加方向箭頭。"
                    ),
                    "sourceId": "source_focal_controls",
                },
                {"id": "focal_status_chart_block", "type": "chart", "chartId": "focal_status_chart"},
                {
                    "id": "focal_table_block",
                    "type": "table",
                    "tableId": "focal_control_table",
                    "layout": "full",
                },
                {
                    "id": "definitions_heading",
                    "type": "markdown",
                    "body": (
                        "## 五個分母不能混用：同一座標不代表同一批 reads\n\n"
                        "M1、G1、tumor-REF、strict edge 與 C++ read-AF 都有不同 eligibility。"
                        "特別是 C++ AF 只累計整個 block pattern weight≥3 的 patterns；它不是 caller VAF。"
                        "H2009 chr13 的 strict pair edge 有 11 個 AA molecules，但多位點 pattern 分散後，"
                        "C++ supported-pattern marginal 對 partner ALT 為 0，正好顯示 selection effect。"
                    ),
                    "sourceId": "source_definitions",
                },
                {
                    "id": "definitions_table_block",
                    "type": "table",
                    "tableId": "definitions_table",
                    "layout": "full",
                },
                {
                    "id": "methods",
                    "type": "markdown",
                    "body": (
                        "## 對位方法：先鎖 coordinate，再保留 PS、HP、W、block 與 read population\n\n"
                        "1. 從 v8 取 `endpoint_a_formal_pair_by_confirmed=true` 的 7 rows。\n"
                        "2. 以 `(dataset, chrom, focal_pos)` 接 tumor-REF control。\n"
                        "3. 以兩 endpoint 在 threshold=3 membership 中共享的 "
                        "`(linkage_basis, exact PS, component_id)` 認定 SAME_W。\n"
                        "4. 逐 HP lane 查 direct endpoint edge，將 raw genomic-order states 轉成 focal-first。\n"
                        "5. 以 sample/chrom/PS/HP 且 C++ `af_coverage` 同時含兩座標，接到 bounded block。"
                    ),
                },
                {
                    "id": "lane_table_block",
                    "type": "table",
                    "tableId": "lane_table",
                    "layout": "full",
                },
                {
                    "id": "limitations",
                    "type": "markdown",
                    "body": (
                        "## 限制與 robustness 判定\n\n"
                        "- **可發布為資料對應：** 7/7 same W、10/10 direct edge、逐 pair G1 表與 focal control 狀態。\n"
                        "- **只能稱 candidate：** 7/24 C++ cohort 雖 technical all-pass，但 "
                        "`all_mutation_bearing_families_complete=false`、`validation_evidence_eligible=false`。\n"
                        "- **局部已解不等於整樹唯一：** H2009 chr18 的 6/6 best trees 共享 focal→partner "
                        "relation，但上游 parent mapping 仍有 6 種。\n"
                        "- **非獨立 VAF 驗證：** read-AF 與 tree constraints 來自同批 supported patterns；"
                        "`input_vaf_eligible=false`，也未校正 CN/LOH、purity、multiplicity 或 CCF。\n"
                        "- **非因果：** methyl group enrichment 可跟 state 共分離，不能證明 partner mutation "
                        "造成 methyl state，亦不能把 MG 直接命名為 clone/subclone。\n"
                        "- **跨技術只到結構：** HCC/DORADO 是同一 cell line 的 technical sources；"
                        "pairwise candidate order 一致，但 methyl multigroup 在 DORADO 不重現。\n"
                        "- **matched-normal 未完成：** 這 7 pair 沒有 formal matched-normal methyl control，"
                        "不能排除既存 ASM 或背景 epiallele。\n"
                        "- **Factorization source identity 已封閉：** current `all7_v2` receipt 綁定"
                        "目前 workspace source，current SHA 與 receipt-bound SHA byte-identical；"
                        "renderer 仍須逐次 fail-closed 驗證 receipt 與 source identity。"
                    ),
                },
                {
                    "id": "next_steps",
                    "type": "markdown",
                    "body": (
                        "## 建議下一步\n\n"
                        "1. 對四個 `ABSTAIN_RESOURCE_LIMIT` pair 做 targeted higher-guard／exact pruning rerun，"
                        "但保留 fail-closed completion receipt。\n"
                        "2. 將 raw marginal、strict pairwise、supported-pattern AF 三種 read denominator "
                        "做 per-locus retention audit，再決定 AF ranking 是否要改成 error-aware likelihood。\n"
                        "3. 對 HCC1395 與 DORADO 同 pair 使用相同 centroid／classifier 做預註冊 methyl replication；"
                        "目前 DORADO K=1 只能作 negative technical result。\n"
                        "4. 補 matched-normal、allele-specific CN/LOH、purity/CCF；最後才用 single-cell、colony "
                        "或 multi-region truth 驗證 cellular lineage。"
                    ),
                },
                {
                    "id": "further_questions",
                    "type": "markdown",
                    "body": (
                        "## 尚待回答\n\n"
                        "- HCC candidate tree 的 focal→partner order 在 error-aware likelihood、CN/CCF校正後是否仍唯一？\n"
                        "- DORADO 的 K=1 是真實甲基單群，還是 coverage、CpG callability 或 clustering power 差異？\n"
                        "- H2009 chr13 的 whole-pattern support filtering 是否會系統性移除低頻 derived state？"
                    ),
                },
            ],
            "cards": [
                {
                    "id": "formal_pairs_card",
                    "description": "從 407,738 個 evaluated rows 中選出的 formal positives；不是準確率。",
                    "dataset": "headline",
                    "sourceId": "source_headline",
                    "metrics": [
                        {"label": "Formal G1 pairs", "field": "formal_pairs", "format": "number"},
                        {"label": "Evaluated pair rows", "field": "evaluated_pairs", "format": "number"},
                    ],
                },
                {
                    "id": "strict_map_card",
                    "description": "全部 formal pairs 在 latest strict graph 中同 W 且有 direct edge。",
                    "dataset": "headline",
                    "sourceId": "source_headline",
                    "metrics": [
                        {"label": "Same W pairs", "field": "same_W_pairs", "format": "number"},
                        {"label": "Direct HP lanes", "field": "direct_lanes", "format": "number"},
                        {
                            "label": "Signal direct reads",
                            "field": "signal_direct_reads",
                            "format": "number",
                        },
                        {
                            "label": "RR-only background reads",
                            "field": "background_direct_reads",
                            "format": "number",
                        },
                    ],
                },
                {
                    "id": "focal_control_card",
                    "description": "Focal ALT+REF 重新 clustering 後通過 allele-axis gate。",
                    "dataset": "headline",
                    "sourceId": "source_headline",
                    "metrics": [
                        {"label": "Focal REF/ALT aligned", "field": "focal_aligned", "format": "number"},
                        {"label": "Joint testable", "field": "focal_testable", "format": "number"},
                    ],
                },
                {
                    "id": "candidate_tree_card",
                    "description": "7/24 C++ model內 unique read-AF winner；不是 true-tree count。",
                    "dataset": "headline",
                    "sourceId": "source_headline",
                    "metrics": [
                        {"label": "Unique candidate T pairs", "field": "unique_candidate_pairs", "format": "number"},
                        {
                            "label": "Resolved pair relations",
                            "field": "resolved_pair_relations",
                            "format": "number",
                        },
                        {"label": "Resource abstain pairs", "field": "resource_abstain_pairs", "format": "number"},
                    ],
                },
                {
                    "id": "release_ceiling_card",
                    "description": "Cohort有 guard-triggered abstain，不能宣稱所有family皆解出。",
                    "dataset": "headline",
                    "sourceId": "source_headline",
                    "metrics": [
                        {"label": "Release-grade complete datasets", "field": "release_grade_complete", "format": "number"},
                    ],
                },
            ],
            "charts": [
                {
                    "id": "evidence_ladder_chart",
                    "title": "七個 formal pair 在各證據層的可確認數",
                    "subtitle": "分母=7 formal G1 positives；不是全體候選 pair 的準確率。",
                    "type": "bar",
                    "dataset": "evidence_ladder",
                    "sourceId": "source_evidence_ladder",
                    "intent": "comparison",
                    "encodings": {
                        "x": {"field": "stage", "type": "nominal", "label": "Evidence stage"},
                        "y": {"field": "pair_count", "type": "quantitative", "label": "Pairs"},
                    },
                    "palette": {"kind": "categorical"},
                },
                {
                    "id": "topology_status_chart",
                    "title": "Formal pairs 的最新 C++ candidate-topology 狀態",
                    "subtitle": "2 unique、1 tied、4 resource abstain；unique 為 model-conditional。",
                    "type": "bar",
                    "dataset": "topology_status",
                    "sourceId": "source_topology_status",
                    "intent": "comparison",
                    "encodings": {
                        "x": {"field": "status", "type": "nominal", "label": "Status"},
                        "y": {"field": "pairs", "type": "quantitative", "label": "Pairs"},
                    },
                    "palette": {"kind": "categorical"},
                },
                {
                    "id": "partner_effect_chart",
                    "title": "Partner REF/ALT × methyl-group 關聯效應",
                    "subtitle": "Cramér's V；7 個已通過 formal G1 gate 的陽性 pair。",
                    "type": "bar",
                    "dataset": "partner_effect",
                    "sourceId": "source_partner_effect",
                    "intent": "comparison",
                    "encodings": {
                        "x": {"field": "locus", "type": "nominal", "label": "Focal locus"},
                        "y": {"field": "cramers_v", "type": "quantitative", "label": "Cramér's V"},
                    },
                    "palette": {"kind": "categorical"},
                },
                {
                    "id": "hcc_methyl_chart",
                    "title": "HCC1395 methyl group × partner allele read counts",
                    "subtitle": "n=51 focal-ALT methyl-core且partner-callable reads；MG 1-2 全為 partner ALT。",
                    "type": "bar",
                    "dataset": "hcc_methyl_counts",
                    "sourceId": "source_hcc_methyl_counts",
                    "intent": "composition",
                    "encodings": {
                        "x": {"field": "methyl_group", "type": "nominal", "label": "Methyl group"},
                        "y": {"field": "reads", "type": "quantitative", "label": "Reads"},
                        "color": {"field": "partner_allele", "type": "nominal", "label": "Partner allele"},
                    },
                    "options": {"grouping": "grouped"},
                    "palette": {"kind": "categorical"},
                },
                {
                    "id": "crosssource_chart",
                    "title": "HCC1395 兩技術來源的 focal-ALT 條件下 partner-ALT 比例",
                    "subtitle": "Strict direct-edge pair states；60.4% vs 53.1%，Fisher p=0.549。",
                    "type": "bar",
                    "dataset": "hcc_crosssource",
                    "sourceId": "source_hcc_crosssource",
                    "intent": "comparison",
                    "encodings": {
                        "x": {"field": "sample", "type": "nominal", "label": "Technical dataset"},
                        "y": {
                            "field": "partner_alt_given_focal_alt",
                            "type": "quantitative",
                            "label": "Partner ALT | focal ALT",
                        },
                    },
                    "valueFormat": "percent",
                    "palette": {"kind": "categorical"},
                },
                {
                    "id": "focal_status_chart",
                    "title": "七個 focal 位點的 ALT+REF methyl control 狀態",
                    "subtitle": "只有一個 axis-aligned；三個因 focal REF不足而不可評估。",
                    "type": "bar",
                    "dataset": "focal_status",
                    "sourceId": "source_focal_status",
                    "intent": "comparison",
                    "encodings": {
                        "x": {"field": "status", "type": "nominal", "label": "Control status"},
                        "y": {"field": "sites", "type": "quantitative", "label": "Sites"},
                    },
                    "palette": {"kind": "categorical"},
                },
            ],
            "tables": [
                {
                    "id": "pair_master_table",
                    "title": "七個 formal G1 pair 的 ALT/REF、W 與最新 candidate T 對應",
                    "subtitle": "Pair order 只在所有 AF-global-best trees 共享同一 relation 時顯示；整樹可仍不唯一。",
                    "dataset": "pair_master",
                    "sourceId": "source_pair_master",
                    "layout": "full",
                    "density": "dense",
                    "defaultSort": {"field": "g1_cramers_v", "direction": "desc"},
                    "columns": pair_columns,
                },
                {
                    "id": "hcc_layers_table",
                    "title": "HCC1395 同一座標的四個 read populations",
                    "subtitle": "每層 denominator 不同，不可直接把數字相加或互換。",
                    "dataset": "hcc_layers",
                    "sourceId": "source_hcc_layers",
                    "layout": "full",
                    "density": "spacious",
                    "defaultSort": {"field": "order", "direction": "asc"},
                    "columns": [
                        {"field": "order", "label": "#", "format": "number"},
                        {"field": "layer", "label": "Layer", "type": "text"},
                        {"field": "population", "label": "Read population", "type": "text"},
                        {"field": "counts", "label": "Observed counts", "type": "text"},
                        {"field": "valid_claim", "label": "Valid annotation", "type": "text"},
                    ],
                },
                {
                    "id": "crosssource_table",
                    "title": "HCC1395 與 DORADO 同 pair 的 strict state 與 candidate T",
                    "subtitle": "同 biological cell line、不同 technical dataset；PS ID 不要求相同。",
                    "dataset": "hcc_crosssource",
                    "sourceId": "source_hcc_crosssource",
                    "layout": "full",
                    "density": "dense",
                    "defaultSort": {"field": "sample", "direction": "asc"},
                    "columns": [
                        {"field": "sample", "label": "Dataset", "type": "text"},
                        {"field": "phase_set", "label": "Exact PS", "type": "text"},
                        {"field": "hp_family", "label": "HP", "type": "text"},
                        {"field": "W_k", "label": "W k", "format": "number"},
                        {"field": "direct_support", "label": "Direct support", "format": "number"},
                        {"field": "state_RR", "label": "RR", "format": "number"},
                        {"field": "state_RA", "label": "RA", "format": "number"},
                        {"field": "state_AR", "label": "AR", "format": "number"},
                        {"field": "state_AA", "label": "AA", "format": "number"},
                        {
                            "field": "partner_alt_given_focal_alt",
                            "label": "P(partner ALT | focal ALT)",
                            "format": "percent",
                        },
                        {"field": "cpp_partner_af", "label": "Supported-pattern partner AF", "format": "percent"},
                        {"field": "cpp_pair_order", "label": "Candidate order", "type": "text"},
                        {"field": "cpp_shape", "label": "Full shape", "type": "text"},
                    ],
                },
                {
                    "id": "focal_control_table",
                    "title": "七個 focal 位點的 tumor ALT+REF methyl joint control",
                    "subtitle": "Permutation p 未做本報告新增的跨位點 multiple correction；axis gate 同時要求 V≥0.30。",
                    "dataset": "focal_controls",
                    "sourceId": "source_focal_controls",
                    "layout": "full",
                    "density": "dense",
                    "defaultSort": {"field": "n_REF", "direction": "desc"},
                    "columns": [
                        {"field": "sample", "label": "Dataset", "type": "text"},
                        {"field": "locus", "label": "Focal", "type": "text"},
                        {"field": "n_ALT", "label": "ALT reads", "format": "number"},
                        {"field": "n_REF", "label": "REF reads", "format": "number"},
                        {"field": "joint_testable", "label": "Joint testable", "type": "boolean"},
                        {"field": "joint_V", "label": "V", "format": "number"},
                        {"field": "joint_p_perm", "label": "p_perm", "format": "number"},
                        {"field": "classification", "label": "Classification", "type": "text"},
                        {"field": "REF_only_stable_multigroup", "label": "REF-only stable multigroup", "type": "boolean"},
                    ],
                },
                {
                    "id": "definitions_table",
                    "title": "ALT/REF 與 topology 各層的正式定義",
                    "subtitle": "同一座標可被不同 read eligibility 與 filtering 觀察。",
                    "dataset": "definitions",
                    "sourceId": "source_definitions",
                    "layout": "full",
                    "density": "spacious",
                    "defaultSort": {"field": "layer", "direction": "asc"},
                    "columns": [
                        {"field": "layer", "label": "Layer", "type": "text"},
                        {"field": "read_population", "label": "Read population", "type": "text"},
                        {"field": "question", "label": "Question", "type": "text"},
                        {"field": "tree_role", "label": "Allowed tree role", "type": "text"},
                    ],
                },
                {
                    "id": "lane_table",
                    "title": "十個 exact-PS×HP lanes 的 strict edge 與 C++ 對位",
                    "subtitle": "RR/RA/AR/AA 均為 focal-first；C++ AF 為 whole-pattern support-filtered marginal。",
                    "dataset": "strict_lanes",
                    "sourceId": "source_strict_lanes",
                    "layout": "full",
                    "density": "dense",
                    "defaultSort": {"field": "direct_support", "direction": "desc"},
                    "columns": [
                        {"field": "pair_label", "label": "Pair", "type": "text"},
                        {"field": "hp_family", "label": "HP", "type": "text"},
                        {"field": "phase_set", "label": "Exact PS", "type": "text"},
                        {"field": "W_k", "label": "W k", "format": "number"},
                        {"field": "cpp_block_k", "label": "Block k", "format": "number"},
                        {"field": "direct_support", "label": "Support", "format": "number"},
                        {"field": "state_RR", "label": "RR", "format": "number"},
                        {"field": "state_RA", "label": "RA", "format": "number"},
                        {"field": "state_AR", "label": "AR", "format": "number"},
                        {"field": "state_AA", "label": "AA", "format": "number"},
                        {"field": "cpp_pair_status", "label": "C++ pair status", "type": "text"},
                        {"field": "cpp_focal_alt", "label": "C++ focal ALT", "format": "number"},
                        {"field": "cpp_partner_alt", "label": "C++ partner ALT", "format": "number"},
                        {
                            "field": "factorization_resolution_class",
                            "label": "Factorization",
                            "type": "text",
                        },
                        {
                            "field": "cpp_pair_order_across_best",
                            "label": "Pair order across best T",
                            "type": "text",
                        },
                    ],
                },
            ],
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "status": "ready",
            "generatedAt": generated_at,
            "datasets": {
                "headline": [
                    {
                        "formal_pairs": 7,
                        "evaluated_pairs": 407738,
                        "same_W_pairs": 7,
                        "direct_lanes": 10,
                        "signal_direct_reads": signal_support,
                        "background_direct_reads": background_support,
                        "focal_aligned": aligned_count,
                        "focal_testable": sum(row["joint_testable"] for row in focal_rows),
                        "unique_candidate_pairs": unique_count,
                        "resolved_pair_relations": resolved_pair_count,
                        "resource_abstain_pairs": topology_counts["ABSTAIN_RESOURCE_LIMIT"],
                        "release_grade_complete": 0,
                    }
                ],
                "evidence_ladder": evidence_ladder,
                "topology_status": topology_status,
                "partner_effect": partner_effect,
                "pair_master": pair_rows,
                "hcc_methyl_counts": hcc_groups,
                "hcc_layers": [
                    {
                        "order": 1,
                        "layer": "Tumor focal REF/ALT control",
                        "population": "focal methyl matrix",
                        "counts": "ALT=98, REF=2",
                        "valid_claim": "NOT_EVALUABLE：REF<3 joint、REF<6 REF-only",
                    },
                    {
                        "order": 2,
                        "layer": "G1 methyl × partner allele",
                        "population": "focal-ALT methyl core + partner callable",
                        "counts": "n=51；[[21,11],[0,19]]",
                        "valid_claim": "MG 1-2 與 partner-ALT enrichment 共分離",
                    },
                    {
                        "order": 3,
                        "layer": "Strict HP2 direct edge",
                        "population": "exact PS fixed endpoint molecules",
                        "counts": "support=55；RR/RA/AR/AA=2/0/21/32",
                        "valid_claim": "同 W 且 direct read-link；無向",
                    },
                    {
                        "order": 4,
                        "layer": "C++ supported-pattern AF",
                        "population": "whole-block pattern weight≥3",
                        "counts": "focal 95/95 ALT；partner 43/79 ALT",
                        "valid_claim": "44 minimum trees中 unique AF candidate",
                    },
                ],
                "hcc_crosssource": cross_rows,
                "focal_status": focal_status,
                "focal_controls": focal_rows,
                "definitions": definitions,
                "strict_lanes": lanes,
            },
        },
        "sources": sources,
    }
    return artifact


def main() -> None:
    pairs = load_formal_pairs()
    controls = load_focal_controls(pairs)
    lanes, lanes_by_pair = build_strict_and_cpp(pairs)
    pair_rows, focal_rows = build_pair_rows(pairs, controls, lanes_by_pair)

    hcc_pair = next(row for row in pairs if row["sample"] == "HCC1395")
    hcc_lane_key = (
        "HCC1395",
        hcc_pair["focal_chrom"],
        hcc_pair["focal_pos_i"],
        hcc_pair["partner_pos_i"],
    )
    cross_rows = build_dorado_crosssource(hcc_pair, lanes_by_pair[hcc_lane_key])

    m1_keys = {
        ("HCC1395", "chr3", HCC_FOCAL),
        ("HCC1395_DORADO", "chr3", HCC_FOCAL),
    }
    hcc_m1 = load_m1_rows(m1_keys)
    if set(hcc_m1) != m1_keys:
        raise AssertionError("missing HCC/DORADO focal M1 rows")
    if not parse_bool(hcc_m1[("HCC1395", "chr3", HCC_FOCAL)]["stable_null_multigroup"]):
        raise AssertionError("HCC focal M1 is expected to be stable multigroup")
    if parse_bool(hcc_m1[("HCC1395_DORADO", "chr3", HCC_FOCAL)]["stable_null_multigroup"]):
        raise AssertionError("DORADO focal M1 is expected to be stable single-group")

    topology_counts = Counter(row["candidate_T_status"] for row in pair_rows)
    focal_counts = Counter(row["classification"] for row in focal_rows)
    cpp_receipt = json.loads(CPP_COHORT_RECEIPT.read_text(encoding="utf-8"))
    factorization_receipt = json.loads(
        FACTORIZATION_RECEIPT.read_text(encoding="utf-8")
    )
    factorization_recorded_source_sha = factorization_receipt["implementation"][
        "candidate_factorization_source"
    ]["sha256"]
    factorization_current_source_sha = sha256(FACTORIZATION_SOURCE_CURRENT)
    topology_0723 = json.loads(TOPOLOGY_AUDIT_0723.read_text(encoding="utf-8"))
    resolved_pair_rows = [
        row
        for row in pair_rows
        if row["candidate_pair_relation_status"]
        == "RESOLVED_COMMON_ACROSS_BEST_TREES"
    ]
    h2009_chr18_pair = next(
        row
        for row in pair_rows
        if row["sample"] == "H2009" and row["focal"].startswith("chr18:")
    )
    signal_direct_support = sum(
        row["direct_support"]
        for row in lanes
        if row["lane_role"] == "FORMAL_SIGNAL"
    )
    background_direct_support = sum(
        row["direct_support"]
        for row in lanes
        if row["lane_role"] == "PAIRED_BACKGROUND"
    )

    checks = {
        "formal_pairs_eq_7": len(pairs) == 7,
        "formal_pair_rows_eq_7": len(pair_rows) == 7,
        "same_W_pairs_eq_7": sum(row["same_W"] for row in pair_rows) == 7,
        "strict_lanes_eq_10": len(lanes) == 10,
        "formal_signal_lanes_eq_7": sum(
            row["lane_role"] == "FORMAL_SIGNAL" for row in lanes
        )
        == 7,
        "paired_background_lanes_eq_3": sum(
            row["lane_role"] == "PAIRED_BACKGROUND" for row in lanes
        )
        == 3,
        "each_pair_has_exactly_one_signal_lane": all(
            sum(row["lane_role"] == "FORMAL_SIGNAL" for row in pair_lanes) == 1
            for pair_lanes in lanes_by_pair.values()
        ),
        "direct_support_sum_eq_790": sum(row["direct_support"] for row in lanes) == 790,
        "signal_direct_support_sum_eq_638": signal_direct_support == 638,
        "background_direct_support_sum_eq_152": background_direct_support == 152,
        "unique_candidate_pairs_eq_2": topology_counts["UNIQUE_AF_CANDIDATE"] == 2,
        "tied_candidate_pairs_eq_1": topology_counts["AF_TIED"] == 1,
        "resource_abstain_pairs_eq_4": topology_counts["ABSTAIN_RESOURCE_LIMIT"] == 4,
        "resolved_pair_relations_eq_3": len(resolved_pair_rows) == 3,
        "all_resolved_pairs_focal_before_partner": {
            row["candidate_pair_order"] for row in resolved_pair_rows
        }
        == {"FOCAL_BEFORE_PARTNER"},
        "h2009_chr18_tied_but_pair_relation_6_of_6": (
            h2009_chr18_pair["candidate_T_status"] == "AF_TIED"
            and h2009_chr18_pair["candidate_T_best_ties"] == 6
            and h2009_chr18_pair["candidate_pair_order_focal_before_count"] == 6
            and h2009_chr18_pair["candidate_pair_order"] == "FOCAL_BEFORE_PARTNER"
        ),
        "focal_joint_testable_eq_3": sum(row["joint_testable"] for row in focal_rows) == 3,
        "focal_axis_aligned_eq_1": focal_counts["ALIGNED"] == 1,
        "hcc_focal_ref_eq_2": next(
            row for row in focal_rows if row["sample"] == "HCC1395"
        )["n_REF"]
        == 2,
        "hcc_candidate_unique": next(
            row for row in pair_rows if row["sample"] == "HCC1395"
        )["candidate_T_status"]
        == "UNIQUE_AF_CANDIDATE",
        "hcc_dorado_candidate_orders_match": len(
            {row["cpp_pair_order"] for row in cross_rows}
        )
        == 1,
        "hcc_dorado_candidate_order_focal_first": cross_rows[0]["cpp_pair_order"]
        == "FOCAL_BEFORE_PARTNER",
        "hcc_m1_multigroup_only_canonical": parse_bool(
            hcc_m1[("HCC1395", "chr3", HCC_FOCAL)]["stable_null_multigroup"]
        )
        and not parse_bool(
            hcc_m1[("HCC1395_DORADO", "chr3", HCC_FOCAL)][
                "stable_null_multigroup"
            ]
        ),
        "cpp_technical_all_pass": cpp_receipt["technical_all_pass"] is True,
        "cpp_validation_evidence_ineligible": cpp_receipt["validation_evidence_eligible"]
        is False,
        "cpp_family_completeness_false": cpp_receipt[
            "all_mutation_bearing_families_complete"
        ]
        is False,
        "factorization_all_pass": factorization_receipt["all_pass"] is True,
        "factorization_ranked_units_eq_71955": factorization_receipt["totals"][
            "ranked_units"
        ]
        == 71955,
        "factorization_source_identity_matches_current": (
            factorization_current_source_sha == factorization_recorded_source_sha
        ),
        "july23_directed_topology_was_zero": topology_0723["summary"][
            "strict_directed_topology_production_complete_datasets"
        ]
        == 0,
    }
    if not all(checks.values()):
        raise AssertionError(f"validation checks failed: {[k for k, v in checks.items() if not v]}")

    pair_fields = [
        "sample",
        "focal",
        "partner",
        "locus_short",
        "g1_informative_reads",
        "g1_groups",
        "g1_table",
        "g1_cramers_v",
        "g1_delta_alt_fraction",
        "g1_q_by",
        "g1_conditional_p",
        "g1_enriched_group",
        "same_W",
        "W_containers",
        "direct_support_total",
        "focal_ref_alt_status",
        "focal_alt_reads",
        "focal_ref_reads",
        "focal_joint_v",
        "focal_joint_p",
        "candidate_T_status",
        "candidate_T_minimum_trees",
        "candidate_T_best_ties",
        "candidate_T_shape",
        "candidate_topology_resolution",
        "candidate_pair_relation_status",
        "candidate_pair_order",
        "candidate_pair_order_focal_before_count",
        "candidate_pair_order_partner_before_count",
        "candidate_pair_order_incomparable_count",
        "endpoint_b_status",
    ]
    lane_fields = [
        "sample",
        "chrom",
        "focal_pos",
        "partner_pos",
        "pair_label",
        "linkage_basis",
        "hp_family",
        "phase_set",
        "component_id",
        "W_k",
        "W_start",
        "W_end",
        "direct_support",
        "lane_role",
        "state_RR",
        "state_RA",
        "state_AR",
        "state_AA",
        "cpp_region_id",
        "cpp_block_k",
        "cpp_active_k",
        "cpp_family_status",
        "cpp_unit_status",
        "cpp_read_af_status",
        "cpp_pair_status",
        "cpp_minimum_vertex_sets",
        "cpp_total_trees",
        "cpp_best_tree_ties",
        "cpp_best_tree_unique",
        "cpp_best_score",
        "cpp_shape",
        "cpp_pair_order",
        "factorization_resolution_class",
        "factorization_global_best_signature_count",
        "factorization_coarse_classes",
        "cpp_pair_relation_status",
        "cpp_pair_order_across_best",
        "cpp_pair_order_focal_before_count",
        "cpp_pair_order_partner_before_count",
        "cpp_pair_order_incomparable_count",
        "cpp_reason",
        "cpp_focal_ref",
        "cpp_focal_alt",
        "cpp_partner_ref",
        "cpp_partner_alt",
        "cpp_input_vaf_eligible",
    ]
    focal_fields = list(focal_rows[0])
    cross_fields = list(cross_rows[0])
    write_tsv(PAIR_OUT, pair_rows, pair_fields)
    write_tsv(LANE_OUT, lanes, lane_fields)
    write_tsv(FOCAL_OUT, focal_rows, focal_fields)
    write_tsv(CROSS_OUT, cross_rows, cross_fields)

    artifact = make_artifact(pairs, pair_rows, focal_rows, lanes, cross_rows, hcc_m1)
    executed_queries = write_sqlite_snapshot(artifact["snapshot"]["datasets"])
    if set(executed_queries) != set(artifact["snapshot"]["datasets"]):
        raise AssertionError("SQLite query coverage does not match artifact datasets")
    checks["sqlite_dataset_query_coverage"] = True
    checks["sqlite_snapshot_row_counts_verified"] = True
    write_json(ARTIFACT_OUT, artifact)

    receipt = {
        "schema_name": "intersubmod.methyl_alt_ref_topology_overlay.validation_receipt",
        "schema_version": "1.1.0",
        "created_at": now_iso(),
        "task_type": "B_comprehensive_validation",
        "goals": ["G3", "G4"],
        "scope": {
            "formal_G1_pairs": 7,
            "strict_HP_lanes": 10,
            "focal_REF_ALT_controls": 7,
            "HCC_crosssource_pair": True,
        },
        "inputs": {
            "formal_pairs": {"path": str(PAIR_TSV), "sha256": sha256(PAIR_TSV)},
            "focal_controls": {"path": str(REF_TSV), "sha256": sha256(REF_TSV)},
            "focal_control_summary": {
                "path": str(REF_SUMMARY),
                "sha256": sha256(REF_SUMMARY),
            },
            "cpp_cohort_receipt": {
                "path": str(CPP_COHORT_RECEIPT),
                "sha256": sha256(CPP_COHORT_RECEIPT),
            },
            "candidate_factorization_receipt": {
                "path": str(FACTORIZATION_RECEIPT),
                "sha256": sha256(FACTORIZATION_RECEIPT),
            },
            "candidate_factorization_current_source": {
                "path": str(FACTORIZATION_SOURCE_CURRENT),
                "current_sha256": factorization_current_source_sha,
                "receipt_bound_sha256": factorization_recorded_source_sha,
                "byte_identical_to_receipt": (
                    factorization_current_source_sha
                    == factorization_recorded_source_sha
                ),
            },
            "topology_audit_20260723": {
                "path": str(TOPOLOGY_AUDIT_0723),
                "sha256": sha256(TOPOLOGY_AUDIT_0723),
            },
        },
        "checks": checks,
        "all_pass": all(checks.values()),
        "headline": {
            "formal_source_rows_evaluated": 407738,
            "formal_partner_G1_pairs": 7,
            "same_W_pairs": 7,
            "direct_HP_lanes": 10,
            "formal_signal_lanes": 7,
            "paired_background_lanes": 3,
            "direct_support_sum": 790,
            "signal_direct_support_sum": signal_direct_support,
            "background_direct_support_sum": background_direct_support,
            "focal_REF_ALT_joint_testable": 3,
            "focal_REF_ALT_axis_aligned": 1,
            "latest_unique_AF_candidate_pairs": 2,
            "latest_AF_tied_pairs": 1,
            "latest_resource_abstain_pairs": 4,
            "latest_pair_relations_resolved_across_all_best_trees": 3,
            "H2009_chr18_focal_before_partner_best_trees": "6/6",
            "HCC_DORADO_pairwise_candidate_order_match": True,
            "HCC_DORADO_methyl_multigroup_replication": False,
        },
        "claim_ceiling": (
            "May report exact data-to-structure overlay, same-W/direct-edge concordance, "
            "model-conditional candidate topology, and HCC/DORADO pairwise candidate-order "
            "agreement. A pair relation shared by every best tree may be reported even when "
            "the full tree is tied. Must not claim causal methylation, a true unique cellular tree, "
            "independent VAF validation, clone identity, or methyl replication in DORADO."
        ),
        "reproducibility": {
            "factorization_output_and_binary_identity": "receipt_verified",
            "factorization_current_source_byte_identical_to_receipt": True,
            "required_action_for_byte_identical_replay": (
                "None; all7_v2 receipt binds the current source at "
                f"{factorization_recorded_source_sha}."
            ),
        },
        "outputs": {},
    }
    write_json(RECEIPT_OUT, receipt)
    receipt["outputs"] = {
        "pair_join": {"path": str(PAIR_OUT), "sha256": sha256(PAIR_OUT)},
        "lane_join": {"path": str(LANE_OUT), "sha256": sha256(LANE_OUT)},
        "focal_control_join": {"path": str(FOCAL_OUT), "sha256": sha256(FOCAL_OUT)},
        "crosssource_join": {"path": str(CROSS_OUT), "sha256": sha256(CROSS_OUT)},
        "artifact": {"path": str(ARTIFACT_OUT), "sha256": sha256(ARTIFACT_OUT)},
        "sqlite_snapshot": {"path": str(SQLITE_OUT), "sha256": sha256(SQLITE_OUT)},
        "query_directory": {
            "path": str(QUERY_DIR),
            "query_count": len(executed_queries),
            "query_sha256": {
                dataset: sha256(QUERY_DIR / f"{dataset}.sql")
                for dataset in sorted(executed_queries)
            },
        },
    }
    write_json(RECEIPT_OUT, receipt)

    print(
        json.dumps(
            {
                "all_pass": receipt["all_pass"],
                "checks": len(checks),
                "formal_pairs": len(pair_rows),
                "strict_lanes": len(lanes),
                "unique_candidate_pairs": topology_counts["UNIQUE_AF_CANDIDATE"],
                "tied_candidate_pairs": topology_counts["AF_TIED"],
                "resource_abstain_pairs": topology_counts["ABSTAIN_RESOURCE_LIMIT"],
                "resolved_pair_relations": len(resolved_pair_rows),
                "focal_axis_aligned": focal_counts["ALIGNED"],
                "signal_direct_support": signal_direct_support,
                "background_direct_support": background_direct_support,
                "hcc_dorado_order": cross_rows[0]["cpp_pair_order"],
                "artifact": str(ARTIFACT_OUT),
                "receipt": str(RECEIPT_OUT),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
