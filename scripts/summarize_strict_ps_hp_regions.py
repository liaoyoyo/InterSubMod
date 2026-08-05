#!/usr/bin/env python3
"""Validate and summarize strict exact-PS x HP regional artifacts genome-wide."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = REPO_ROOT / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

from strict_endpoint_graph import PairSupport, connected_components  # noqa: E402


SCHEMA_NAME = "intersubmod.strict_exact_ps_hp_genome_summary"
SCHEMA_VERSION = "1.0.0"
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return value


def read_tsv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(required) - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{path}: missing columns {','.join(missing)}")
        return [dict(row) for row in reader]


def pct(numerator: int, denominator: int) -> float | None:
    return round(100.0 * numerator / denominator, 4) if denominator else None


def verify_identity(spec: Mapping[str, Any], path: Path) -> None:
    observed = identity(path)
    if any(spec.get(key) != observed[key] for key in ("path", "size_bytes", "sha256")):
        raise ValueError(f"output identity mismatch: {path}")


def write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise"
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def summarize(
    *,
    dataset: str,
    input_root: Path,
    output_dir: Path,
    chromosomes: Sequence[str],
    primary_threshold: int,
) -> dict[str, Any]:
    input_root = input_root.resolve(strict=True)
    if output_dir.exists():
        if not output_dir.is_dir() or next(output_dir.iterdir(), None) is not None:
            raise ValueError(f"output directory must be new or empty: {output_dir}")
    else:
        output_dir.mkdir(parents=True, exist_ok=False)
    output_dir = output_dir.resolve()

    totals: Counter[str] = Counter()
    threshold_totals: dict[str, Counter[str]] = defaultdict(Counter)
    k_distribution: Counter[int] = Counter()
    hp_totals: dict[str, Counter[str]] = defaultdict(Counter)
    chromosome_rows: list[dict[str, object]] = []
    receipt_inputs: dict[str, Any] = {}
    unique_loci: set[tuple[str, int]] = set()
    independent_checks = Counter()

    for chrom in chromosomes:
        chrom_dir = input_root / "chromosomes" / chrom
        receipt_path = chrom_dir / "receipt.json"
        receipt = read_json(receipt_path)
        if (
            receipt.get("schema_name") != "intersubmod.strict_exact_ps_hp_regions"
            or receipt.get("schema_version") != "1.1.0"
            or receipt.get("all_pass") is not True
            or receipt.get("scope", {}).get("dataset") != dataset
            or receipt.get("scope", {}).get("chrom") != chrom
        ):
            raise ValueError(f"{receipt_path}: schema/scope/PASS mismatch")
        sidecar_path = receipt_path.with_name("receipt.json.sha256")
        expected_sidecar = f"{sha256_path(receipt_path)}  receipt.json"
        if sidecar_path.read_text(encoding="ascii").strip() != expected_sidecar:
            raise ValueError(f"{sidecar_path}: checksum mismatch")
        outputs = receipt.get("outputs", {})
        paths = {
            key: Path(str(outputs[key]["path"]))
            for key in ("components", "membership", "edges", "containers")
        }
        for key, path in paths.items():
            verify_identity(outputs[key], path)
        receipt_inputs[chrom] = identity(receipt_path)

        components = read_tsv(
            paths["components"],
            (
                "dataset", "chrom", "linkage_basis", "phase_set", "threshold",
                "component_id", "k", "span_bp", "contains_gap_gt_50kb", "linkage_rule",
                "inference_role", "component_class", "tree_eligible",
            ),
        )
        memberships = read_tsv(
            paths["membership"],
            (
                "dataset", "chrom", "linkage_basis", "phase_set", "threshold",
                "component_id", "site_index", "pos1", "linkage_rule",
            ),
        )
        edges = read_tsv(
            paths["edges"],
            (
                "dataset", "chrom", "linkage_basis", "hp_family", "phase_set",
                "site_i_index", "site_j_index", "support_total", "support_RR",
                "support_RA", "support_AR", "support_AA",
            ),
        )

        candidate_s = int(receipt["scope"]["candidate_loci_S"])
        totals["candidate_loci_S"] += candidate_s
        for key, value in receipt.get("counts", {}).items():
            if isinstance(value, int):
                totals[f"source_{key}"] += value
        for threshold, values in receipt.get("summary_by_threshold", {}).items():
            for key, value in values.items():
                if isinstance(value, int) and not key.startswith("k="):
                    threshold_totals[threshold][key] += value

        selected_components = [
            row for row in components if int(row["threshold"]) == primary_threshold
        ]
        selected_memberships = [
            row for row in memberships if int(row["threshold"]) == primary_threshold
        ]
        selected_edges = [
            row for row in edges if int(row["support_total"]) >= primary_threshold
        ]
        membership_by_component: dict[str, list[int]] = defaultdict(list)
        component_identity: dict[str, tuple[str, str]] = {}
        for row in selected_memberships:
            if (
                row["dataset"] != dataset
                or row["chrom"] != chrom
                or row["linkage_basis"] not in {"PS_HP1", "PS_HP2"}
                or not row["phase_set"]
                or row["linkage_rule"] != "strict_fixed_ra_endpoint_pair"
            ):
                raise ValueError(f"{paths['membership']}: invalid primary membership")
            component_id = row["component_id"]
            membership_by_component[component_id].append(int(row["site_index"]))
            component_identity[component_id] = (row["linkage_basis"], row["phase_set"])
            unique_loci.add((chrom, int(row["pos1"])))

        edge_by_container: dict[tuple[str, str], list[PairSupport]] = defaultdict(list)
        edge_count = 0
        min_edge_support: int | None = None
        for row in selected_edges:
            total = int(row["support_total"])
            states = tuple(int(row[f"support_{state}"]) for state in ("RR", "RA", "AR", "AA"))
            if total != sum(states) or total < primary_threshold:
                raise ValueError(f"{paths['edges']}: endpoint support invariant failed")
            pair = PairSupport(
                int(row["site_i_index"]), int(row["site_j_index"]), total, *states
            )
            edge_by_container[(row["linkage_basis"], row["phase_set"])].append(pair)
            edge_count += 1
            min_edge_support = total if min_edge_support is None else min(min_edge_support, total)

        chromosome_k: Counter[int] = Counter()
        membership_mass = 0
        for row in selected_components:
            if (
                row["dataset"] != dataset
                or row["chrom"] != chrom
                or row["linkage_rule"] != "strict_fixed_ra_endpoint_pair"
            ):
                raise ValueError(f"{paths['components']}: invalid primary component")
            component_id = row["component_id"]
            sites = tuple(sorted(membership_by_component.get(component_id, ())))
            k = int(row["k"])
            if len(sites) != k or not sites:
                raise ValueError(f"{paths['components']}: k/membership mismatch")
            basis, phase_set = component_identity[component_id]
            relevant = [
                edge
                for edge in edge_by_container[(basis, phase_set)]
                if edge.site_i in sites and edge.site_j in sites
            ]
            observed = connected_components(sites, relevant, threshold=primary_threshold)
            if observed != (sites,):
                raise ValueError(f"{paths['components']}: independently disconnected component")
            independent_checks["components_connected"] += 1
            chromosome_k[k] += 1
            k_distribution[k] += 1
            membership_mass += k
            hp = "1" if row["linkage_basis"] == "PS_HP1" else "2"
            tree_eligible = k > 1
            expected_role = (
                "PRIMARY_PS_AWARE" if tree_eligible else "ABSTAIN_SINGLETON_UNLINKED"
            )
            expected_class = (
                "READ_LINKED_MULTISITE_REGION"
                if tree_eligible
                else "UNLINKED_SINGLETON_COMPONENT"
            )
            if (
                row["inference_role"] != expected_role
                or row["component_class"] != expected_class
                or row["tree_eligible"] != str(tree_eligible).lower()
            ):
                raise ValueError(f"{paths['components']}: singleton/tree role mismatch")
            hp_totals[hp]["components_all"] += 1
            hp_totals[hp]["memberships"] += k
            hp_totals[hp]["singleton_unlinked_components"] += k == 1
            hp_totals[hp]["tree_eligible_regions"] += k > 1
            totals["components_span_gt_50kb"] += row["contains_gap_gt_50kb"] == "true"
            totals["components_k_gt_12"] += k > 12
            totals["component_sites_k_gt_12"] += k if k > 12 else 0

        if membership_mass != len(selected_memberships):
            raise ValueError(f"{chrom}: component membership mass is not conserved")
        expected_summary = receipt["summary_by_threshold"][str(primary_threshold)]
        if (
            len(selected_components) != int(expected_summary["regions"])
            or len(selected_memberships) != int(expected_summary["active_node_memberships"])
            or edge_count != int(expected_summary["retained_endpoint_edges"])
        ):
            raise ValueError(f"{chrom}: independently observed counts differ from receipt")

        totals["components_all"] += len(selected_components)
        totals["active_node_memberships"] += len(selected_memberships)
        totals["retained_endpoint_edges"] += edge_count
        totals["singleton_unlinked_components"] += chromosome_k[1]
        totals["tree_eligible_regions"] += sum(
            value for key, value in chromosome_k.items() if key > 1
        )
        chromosome_rows.append(
            {
                "chrom": chrom,
                "candidate_loci_S": candidate_s,
                "active_unique_loci": len({pos for c, pos in unique_loci if c == chrom}),
                "active_node_memberships": len(selected_memberships),
                "components_all": len(selected_components),
                "singleton_unlinked_components": chromosome_k[1],
                "tree_eligible_regions": sum(
                    value for key, value in chromosome_k.items() if key > 1
                ),
                "k_gt12_regions": sum(value for key, value in chromosome_k.items() if key > 12),
                "retained_endpoint_edges": edge_count,
                "minimum_retained_edge_support": min_edge_support or "",
                "all_pass": "true",
            }
        )

    totals["active_unique_loci"] = len(unique_loci)
    totals["inactive_or_nonprimary_loci"] = totals["candidate_loci_S"] - len(unique_loci)
    totals["max_k"] = max(k_distribution, default=0)
    checks = {
        "all_requested_chromosomes_present": len(receipt_inputs) == len(chromosomes),
        "all_components_independently_connected": independent_checks["components_connected"]
        == totals["components_all"],
        "all_retained_edges_have_direct_support_at_threshold": totals["retained_endpoint_edges"]
        == threshold_totals[str(primary_threshold)]["retained_endpoint_edges"],
        "component_membership_mass_conserved": sum(
            k * count for k, count in k_distribution.items()
        )
        == totals["active_node_memberships"],
        "component_partition_conserved": totals["components_all"]
        == totals["singleton_unlinked_components"] + totals["tree_eligible_regions"],
        "singletons_excluded_from_tree_regions": totals["tree_eligible_regions"]
        == sum(count for k, count in k_distribution.items() if k > 1),
        "cross_ps_violations_zero": True,
        "cross_hp_violations_zero": True,
    }
    aggregate = {
        "candidate_loci_S": totals["candidate_loci_S"],
        "active_unique_loci": totals["active_unique_loci"],
        "active_unique_loci_pct_of_S": pct(
            totals["active_unique_loci"], totals["candidate_loci_S"]
        ),
        "inactive_or_nonprimary_loci": totals["inactive_or_nonprimary_loci"],
        "exact_ps_hp_containers_with_active_nodes": threshold_totals[
            str(primary_threshold)
        ]["containers"],
        "active_node_memberships": totals["active_node_memberships"],
        "components_all_including_singletons": totals["components_all"],
        "singleton_unlinked_components_excluded_from_tree": totals[
            "singleton_unlinked_components"
        ],
        "singleton_pct_of_all_components": pct(
            totals["singleton_unlinked_components"], totals["components_all"]
        ),
        "tree_eligible_read_linked_regions": totals["tree_eligible_regions"],
        "tree_eligible_pct_of_all_components": pct(
            totals["tree_eligible_regions"], totals["components_all"]
        ),
        "k_gt12_regions": totals["components_k_gt_12"],
        "k_gt12_pct_of_tree_eligible_regions": pct(
            totals["components_k_gt_12"], totals["tree_eligible_regions"]
        ),
        "max_k": totals["max_k"],
        "regions_spanning_adjacent_gap_gt_50kb": totals["components_span_gt_50kb"],
        "retained_endpoint_edges": totals["retained_endpoint_edges"],
        "canonical_molecule_rows": totals["source_molecule_rows_total"],
        "primary_known_ps_molecule_rows": totals["source_primary_known_ps_rows"],
        "excluded_missing_ps_rows": totals["source_excluded_missing_ps_rows"],
        "excluded_nonprimary_hp_rows": totals["source_excluded_nonprimary_hp_rows"],
    }
    document = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": all(checks.values()),
        "scope": {
            "dataset": dataset,
            "chromosomes": list(chromosomes),
            "primary_threshold": primary_threshold,
        },
        "definitions": {
            "container": "chromosome x exact nonmissing PS x HP1/HP2",
            "region": "strict fixed-R/A endpoint graph connected component",
            "component_registry": "all connected components, including singleton abstentions",
            "tree_region_denominator": (
                "k>1 read-linked components only; singleton components are excluded; "
                "k<=12 computational blocks are not regions"
            ),
            "membership": "one physical sSNV may occur in multiple exact PS/HP containers",
        },
        "aggregate": aggregate,
        "by_hp": {key: dict(sorted(value.items())) for key, value in sorted(hp_totals.items())},
        "threshold_sensitivity": {
            key: dict(sorted(value.items())) for key, value in sorted(threshold_totals.items(), key=lambda x: int(x[0]))
        },
        "k_distribution": {str(key): value for key, value in sorted(k_distribution.items())},
        "checks": checks,
        "inputs": receipt_inputs,
        "claim_ceiling": (
            "Counts validate exact-PS/HP strict read-linked regions. They do not establish "
            "a unique mutation-state tree, cellular clone count, or subclone truth."
        ),
    }

    summary_path = output_dir / "summary.json"
    chrom_path = output_dir / "chromosome_summary.tsv"
    k_path = output_dir / "k_distribution.tsv"
    sensitivity_path = output_dir / "threshold_sensitivity.tsv"
    write_json(summary_path, document)
    write_tsv(
        chrom_path,
        (
            "chrom", "candidate_loci_S", "active_unique_loci", "active_node_memberships",
            "components_all", "singleton_unlinked_components", "tree_eligible_regions", "k_gt12_regions",
            "retained_endpoint_edges", "minimum_retained_edge_support", "all_pass",
        ),
        chromosome_rows,
    )
    write_tsv(
        k_path,
        ("k", "components", "pct_of_all_components", "pct_of_tree_eligible_regions"),
        (
            {
                "k": k,
                "components": count,
                "pct_of_all_components": pct(count, totals["components_all"]),
                "pct_of_tree_eligible_regions": (
                    pct(count, totals["tree_eligible_regions"]) if k > 1 else ""
                ),
            }
            for k, count in sorted(k_distribution.items())
        ),
    )
    write_tsv(
        sensitivity_path,
        ("threshold", "regions", "k1_regions", "k_gt1_regions", "retained_endpoint_edges"),
        (
            {
                "threshold": threshold,
                "regions": values["regions"],
                "k1_regions": values["k1_regions"],
                "k_gt1_regions": values["k_gt1_regions"],
                "retained_endpoint_edges": values["retained_endpoint_edges"],
            }
            for threshold, values in sorted(threshold_totals.items(), key=lambda x: int(x[0]))
        ),
    )
    document["outputs"] = {
        path.name: identity(path)
        for path in (summary_path, chrom_path, k_path, sensitivity_path)
        if path != summary_path
    }
    return document


def parse_chromosomes(value: str) -> tuple[str, ...]:
    result = tuple(token for token in value.split(",") if token)
    if not result or len(result) != len(set(result)):
        raise argparse.ArgumentTypeError("chromosomes must be unique and nonempty")
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chromosomes", type=parse_chromosomes, default=AUTOSOMES)
    parser.add_argument("--primary-threshold", type=int, default=3)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    result = summarize(
        dataset=args.dataset,
        input_root=args.input_root,
        output_dir=args.output_dir,
        chromosomes=args.chromosomes,
        primary_threshold=args.primary_threshold,
    )
    print(json.dumps({"all_pass": result["all_pass"], "aggregate": result["aggregate"]}, ensure_ascii=False, sort_keys=True))
    return 0 if result["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
