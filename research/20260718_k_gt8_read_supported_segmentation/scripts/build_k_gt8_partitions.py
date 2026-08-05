#!/usr/bin/env python3
"""Build deterministic read-supported k<=8 partitions for legacy k>8 regions.

The input is the lossless sparse-call output from the frozen M2 extractor.
Legacy regions are reconstructed exactly from the ordered site catalog using
the historical adjacent-gap rule (new region when gap > 50,000 bp).

Primary constraints are unique-molecule, fixed R/A observations from HP1/HP2
with an exact known phase set. Identical sparse patterns are aggregated, but
HP families and phase sets are never pooled into one pattern. The common
engineering partition may sum independent pattern weights across units; this
does not authorize downstream pooling of their haplotype orientations.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from read_support_partition import (
    CUT,
    RETAINED,
    UNAVOIDABLE,
    Hyperedge,
    PartitionResult,
    partition_ordered_hypergraph,
)


SCHEMA_NAME = "intersubmod.k_gt8_read_supported_segmentation"
SCHEMA_VERSION = "0.1.0"
PRIMARY_HP_FAMILIES = frozenset({"1", "2"})
FIXED_CALLS = frozenset({"R", "A"})


@dataclass(frozen=True)
class Site:
    site_index: int
    chrom: str
    pos1: int
    ref: str
    alt: str


@dataclass(frozen=True)
class LegacyComponent:
    component_id: str
    group_ordinal: int
    sites: tuple[Site, ...]

    @property
    def positions(self) -> tuple[int, ...]:
        return tuple(site.pos1 for site in self.sites)

    @property
    def k(self) -> int:
        return len(self.sites)


@dataclass(frozen=True)
class PatternKey:
    hp_family: str
    phase_set: str
    positions: tuple[int, ...]
    call_codes: str


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def semantic_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _split_csv_ints(value: str) -> tuple[int, ...]:
    if not value:
        return ()
    return tuple(int(token) for token in value.split(","))


def _format_decimal(value: Decimal) -> str:
    return format(value, "f")


def _ratio_or_blank(numerator: Decimal, denominator: Decimal) -> str:
    if denominator == 0:
        return ""
    return _format_decimal(numerator / denominator)


def _deterministic_gzip_writer(path: Path):
    raw = path.open("xb")
    compressed = gzip.GzipFile(
        filename="",
        mode="wb",
        compresslevel=6,
        fileobj=raw,
        mtime=0,
    )
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")

    class _Context:
        def __enter__(self):
            return text

        def __exit__(self, exc_type, exc, traceback):
            try:
                text.close()
            finally:
                if not raw.closed:
                    raw.close()
            return False

    return _Context()


def _write_tsv_gz(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with _deterministic_gzip_writer(path) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_sites(path: Path, expected_chrom: str) -> tuple[Site, ...]:
    sites = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"site_index", "chrom", "pos1", "ref", "alt"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise RuntimeError(f"site catalog missing columns: {sorted(required)}")
        for row in reader:
            site = Site(
                site_index=int(row["site_index"]),
                chrom=row["chrom"],
                pos1=int(row["pos1"]),
                ref=row["ref"],
                alt=row["alt"],
            )
            if site.chrom != expected_chrom:
                raise RuntimeError(f"site catalog chromosome mismatch: {site.chrom}")
            sites.append(site)
    if tuple(site.site_index for site in sites) != tuple(range(len(sites))):
        raise RuntimeError("site_index must be contiguous, zero-based, and ordered")
    positions = tuple(site.pos1 for site in sites)
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise RuntimeError("site catalog positions must be strictly increasing")
    return tuple(sites)


def build_legacy_components(
    dataset: str,
    chrom: str,
    sites: Sequence[Site],
    *,
    legacy_gap_bp: int,
    target_k_min: int,
) -> tuple[tuple[LegacyComponent, ...], int, int]:
    groups: list[list[Site]] = []
    current: list[Site] = []
    for site in sites:
        if current and site.pos1 - current[-1].pos1 > legacy_gap_bp:
            groups.append(current)
            current = []
        current.append(site)
    if current:
        groups.append(current)

    targets = []
    for ordinal, group in enumerate(groups, start=1):
        if len(group) < target_k_min:
            continue
        component_id = (
            f"{dataset}:{chrom}:legacy_gap_{legacy_gap_bp}:"
            f"{ordinal}:{group[0].pos1}-{group[-1].pos1}"
        )
        targets.append(
            LegacyComponent(
                component_id=component_id,
                group_ordinal=ordinal,
                sites=tuple(group),
            )
        )
    n_multilocus = sum(len(group) >= 2 for group in groups)
    n_singleton = sum(len(group) == 1 for group in groups)
    return tuple(targets), n_multilocus, n_singleton


def densest_window(positions: Sequence[int], max_block_size: int) -> tuple[int, ...]:
    if len(positions) <= max_block_size:
        return tuple(positions)
    best_start = min(
        range(len(positions) - max_block_size + 1),
        key=lambda start: positions[start + max_block_size - 1] - positions[start],
    )
    return tuple(positions[best_start : best_start + max_block_size])


def load_primary_patterns(
    calls_path: Path,
    *,
    expected_dataset: str,
    expected_chrom: str,
    target_site_to_component: Mapping[int, tuple[str, int]],
    component_positions: Mapping[str, tuple[int, ...]],
) -> tuple[dict[str, Counter[PatternKey]], dict[str, Any]]:
    patterns: dict[str, Counter[PatternKey]] = defaultdict(Counter)
    counters: Counter[str] = Counter()
    hp_counts: Counter[str] = Counter()
    phase_set_counts: Counter[str] = Counter()
    active_sites: dict[str, set[int]] = defaultdict(set)

    with gzip.open(calls_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset",
            "chrom",
            "molecule_id",
            "hp_family",
            "phase_set",
            "site_indices",
            "call_codes",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise RuntimeError(f"molecule calls missing columns: {sorted(required)}")
        for row in reader:
            counters["molecule_rows_seen"] += 1
            if row["dataset"] != expected_dataset or row["chrom"] != expected_chrom:
                raise RuntimeError("molecule calls dataset/chromosome mismatch")
            hp_family = row["hp_family"]
            phase_set = row["phase_set"]
            if hp_family not in PRIMARY_HP_FAMILIES:
                counters["excluded_nonprimary_hp_rows"] += 1
                continue
            if not phase_set:
                counters["excluded_missing_ps_rows"] += 1
                continue
            counters["primary_known_ps_rows"] += 1
            hp_counts[hp_family] += 1
            phase_set_counts[f"HP{hp_family}|{phase_set}"] += 1

            indices = _split_csv_ints(row["site_indices"])
            codes = tuple(row["call_codes"])
            if len(indices) != len(codes):
                raise RuntimeError("site_indices and call_codes cardinality mismatch")
            by_component: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
            for site_index, code in zip(indices, codes):
                if code not in FIXED_CALLS:
                    continue
                membership = target_site_to_component.get(site_index)
                if membership is None:
                    continue
                component_id, local_index = membership
                by_component[component_id].append((local_index, site_index, code))

            for component_id, observations in by_component.items():
                observations.sort()
                if len(observations) < 2:
                    counters["primary_single_target_site_observations"] += 1
                    continue
                positions = component_positions[component_id]
                pattern_positions = tuple(positions[item[0]] for item in observations)
                pattern_codes = "".join(item[2] for item in observations)
                key = PatternKey(
                    hp_family=hp_family,
                    phase_set=phase_set,
                    positions=pattern_positions,
                    call_codes=pattern_codes,
                )
                patterns[component_id][key] += 1
                active_sites[component_id].update(pattern_positions)
                counters["primary_component_constraints_molecule_weight"] += 1

    diagnostics = {
        "counts": dict(sorted(counters.items())),
        "primary_rows_by_hp_family": dict(sorted(hp_counts.items())),
        "primary_rows_by_hp_phase_set_count": len(phase_set_counts),
        "primary_active_sites_by_component": {
            component_id: tuple(sorted(values))
            for component_id, values in sorted(active_sites.items())
        },
    }
    return dict(patterns), diagnostics


def _pattern_constraint_id(component_id: str, key: PatternKey) -> str:
    payload = {
        "component_id": component_id,
        "hp_family": key.hp_family,
        "phase_set": key.phase_set,
        "positions": key.positions,
        "call_codes": key.call_codes,
    }
    return "P" + semantic_sha256(payload)[:24]


def build_hyperedges(
    component_id: str,
    pattern_counts: Mapping[PatternKey, int],
    weighting: str,
) -> tuple[Hyperedge, ...]:
    edges = []
    for key, count in sorted(
        pattern_counts.items(),
        key=lambda item: (
            item[0].hp_family,
            item[0].phase_set,
            item[0].positions,
            item[0].call_codes,
        ),
    ):
        if weighting == "raw_molecule":
            weight = Decimal(count)
        elif weighting == "equal_pattern":
            weight = Decimal(1)
        elif weighting == "log1p_molecule":
            with localcontext() as context:
                context.prec = 28
                weight = Decimal(count + 1).ln()
        else:
            raise ValueError(f"unknown weighting: {weighting}")
        edges.append(
            Hyperedge(
                constraint_id=_pattern_constraint_id(component_id, key),
                sites=key.positions,
                weight=weight,
                pattern_count=1,
            )
        )
    return tuple(edges)


def partition_status(
    raw_result: PartitionResult,
    *,
    active_site_count: int,
    component_k: int,
    weight_stable: bool,
) -> str:
    if raw_result.total_pattern_count == 0:
        return "ABSTAIN_NO_PRIMARY_LINKAGE"
    if raw_result.lost_pattern_count:
        return (
            "PARTIAL_LOCAL_ONLY_WEIGHT_STABLE"
            if weight_stable
            else "PARTIAL_LOCAL_ONLY_WEIGHT_SENSITIVE"
        )
    coverage = "FULL_SITE_COVERAGE" if active_site_count == component_k else "PARTIAL_SITE_COVERAGE"
    stability = "WEIGHT_STABLE" if weight_stable else "WEIGHT_SENSITIVE"
    return f"OBSERVED_CONSTRAINT_NO_LOSS_{coverage}_{stability}"


def baseline_retention(
    result_edges: Sequence[Hyperedge],
    selected_positions: Sequence[int],
) -> tuple[Decimal, int]:
    selected = set(selected_positions)
    retained = tuple(edge for edge in result_edges if set(edge.sites).issubset(selected))
    return (
        sum((edge.weight for edge in retained), Decimal(0)),
        sum(edge.pattern_count for edge in retained),
    )


def output_identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if any(args.output_dir.iterdir()):
        raise RuntimeError(f"output directory must be empty: {args.output_dir}")

    stage_started = time.perf_counter()
    sites = load_sites(args.site_catalog, args.chrom)
    components, n_multilocus, n_singleton = build_legacy_components(
        args.dataset,
        args.chrom,
        sites,
        legacy_gap_bp=args.legacy_gap_bp,
        target_k_min=args.max_block_size + 1,
    )
    site_load_seconds = time.perf_counter() - stage_started

    component_positions = {
        component.component_id: component.positions for component in components
    }
    site_to_component = {}
    for component in components:
        for local_index, site in enumerate(component.sites):
            if site.site_index in site_to_component:
                raise AssertionError("target site belongs to multiple legacy components")
            site_to_component[site.site_index] = (component.component_id, local_index)

    stage_started = time.perf_counter()
    patterns_by_component, pattern_diagnostics = load_primary_patterns(
        args.molecule_calls,
        expected_dataset=args.dataset,
        expected_chrom=args.chrom,
        target_site_to_component=site_to_component,
        component_positions=component_positions,
    )
    pattern_load_seconds = time.perf_counter() - stage_started

    stage_started = time.perf_counter()
    component_rows = []
    block_rows = []
    membership_rows = []
    constraint_rows = []
    status_counts: Counter[str] = Counter()
    global_counts: Counter[str] = Counter()
    semantic_components = []
    all_sites_assigned: list[int] = []

    for component in components:
        pattern_counts = patterns_by_component.get(component.component_id, Counter())
        raw_edges = build_hyperedges(component.component_id, pattern_counts, "raw_molecule")
        equal_edges = build_hyperedges(component.component_id, pattern_counts, "equal_pattern")
        log_edges = build_hyperedges(component.component_id, pattern_counts, "log1p_molecule")
        raw_result = partition_ordered_hypergraph(
            component.positions, raw_edges, max_block_size=args.max_block_size
        )
        equal_result = partition_ordered_hypergraph(
            component.positions, equal_edges, max_block_size=args.max_block_size
        )
        log_result = partition_ordered_hypergraph(
            component.positions, log_edges, max_block_size=args.max_block_size
        )
        weight_stable = (
            raw_result.cut_indices == equal_result.cut_indices == log_result.cut_indices
        )
        active_sites = set(
            pattern_diagnostics["primary_active_sites_by_component"].get(
                component.component_id, ()
            )
        )
        status = partition_status(
            raw_result,
            active_site_count=len(active_sites),
            component_k=component.k,
            weight_stable=weight_stable,
        )
        status_counts[status] += 1

        selected = densest_window(component.positions, args.max_block_size)
        baseline_weight, baseline_patterns = baseline_retention(raw_edges, selected)
        result_by_constraint = {
            disposition.constraint_id: disposition
            for disposition in raw_result.dispositions
        }
        pattern_by_id = {
            _pattern_constraint_id(component.component_id, key): (key, count)
            for key, count in pattern_counts.items()
        }
        block_by_position = {}
        for block in raw_result.blocks:
            block_id = (
                f"{component.component_id}:B{block.block_index + 1}:"
                f"{block.positions[0]}-{block.positions[-1]}"
            )
            for position in block.positions:
                if position in block_by_position:
                    raise AssertionError("position assigned to multiple blocks")
                block_by_position[position] = block_id
            block_rows.append(
                {
                    "dataset": args.dataset,
                    "chrom": args.chrom,
                    "legacy_component_id": component.component_id,
                    "block_index": block.block_index + 1,
                    "block_id": block_id,
                    "start1": block.positions[0],
                    "end1": block.positions[-1],
                    "k": block.k,
                    "positions": ",".join(map(str, block.positions)),
                    "component_status": status,
                    "raw_retained_molecule_weight": _format_decimal(
                        block.retained_weight
                    ),
                    "retained_exact_pattern_count": block.retained_pattern_count,
                    "primary_active_site_count": sum(
                        position in active_sites for position in block.positions
                    ),
                }
            )
        for local_index, site in enumerate(component.sites):
            block_id = block_by_position[site.pos1]
            all_sites_assigned.append(site.site_index)
            membership_rows.append(
                {
                    "dataset": args.dataset,
                    "chrom": args.chrom,
                    "legacy_component_id": component.component_id,
                    "site_index": site.site_index,
                    "component_local_index": local_index,
                    "pos1": site.pos1,
                    "ref": site.ref,
                    "alt": site.alt,
                    "block_id": block_id,
                    "primary_linkage_observed": str(site.pos1 in active_sites).lower(),
                    "old_densest8_disposition": (
                        "selected" if site.pos1 in set(selected) else "cap_excluded"
                    ),
                }
            )

        for constraint_id in sorted(pattern_by_id):
            key, molecule_weight = pattern_by_id[constraint_id]
            disposition = result_by_constraint[constraint_id]
            unit_id = f"HP{key.hp_family}|PS{key.phase_set}"
            constraint_rows.append(
                {
                    "dataset": args.dataset,
                    "chrom": args.chrom,
                    "legacy_component_id": component.component_id,
                    "constraint_id": constraint_id,
                    "unit_id": unit_id,
                    "hp_family": key.hp_family,
                    "phase_set": key.phase_set,
                    "positions": ",".join(map(str, key.positions)),
                    "call_codes": key.call_codes,
                    "n_fixed_ra": len(key.positions),
                    "span_sites": disposition.span_sites,
                    "molecule_weight": molecule_weight,
                    "disposition": disposition.status,
                    "crossed_cut_count": disposition.crossed_cut_count,
                    "retained_block_index": (
                        ""
                        if disposition.block_index is None
                        else disposition.block_index + 1
                    ),
                }
            )

        component_row = {
            "dataset": args.dataset,
            "chrom": args.chrom,
            "legacy_component_id": component.component_id,
            "legacy_group_ordinal": component.group_ordinal,
            "start1": component.positions[0],
            "end1": component.positions[-1],
            "span_bp": component.positions[-1] - component.positions[0],
            "pre_cap_k": component.k,
            "old_densest8_selected": len(selected),
            "old_cap_excluded": component.k - len(selected),
            "new_block_count": len(raw_result.blocks),
            "new_site_retained": component.k,
            "primary_active_site_count": len(active_sites),
            "primary_active_site_fraction": (
                f"{len(active_sites) / component.k:.12f}"
            ),
            "exact_pattern_count": raw_result.total_pattern_count,
            "raw_total_molecule_weight": _format_decimal(raw_result.total_weight),
            "raw_retained_molecule_weight": _format_decimal(
                raw_result.retained_weight
            ),
            "raw_lost_molecule_weight": _format_decimal(raw_result.lost_weight),
            "raw_retention_ratio": _ratio_or_blank(
                raw_result.retained_weight, raw_result.total_weight
            ),
            "old_densest8_retained_molecule_weight": _format_decimal(
                baseline_weight
            ),
            "old_densest8_retention_ratio": _ratio_or_blank(
                baseline_weight, raw_result.total_weight
            ),
            "retained_exact_pattern_count": raw_result.retained_pattern_count,
            "lost_exact_pattern_count": raw_result.lost_pattern_count,
            "unavoidable_pattern_count": len(
                raw_result.unavoidable_constraint_ids
            ),
            "old_densest8_retained_pattern_count": baseline_patterns,
            "weight_stable": str(weight_stable).lower(),
            "raw_cut_indices": ",".join(map(str, raw_result.cut_indices)),
            "equal_pattern_cut_indices": ",".join(
                map(str, equal_result.cut_indices)
            ),
            "log1p_cut_indices": ",".join(map(str, log_result.cut_indices)),
            "cut_gap_sum_bp": sum(raw_result.cut_gaps),
            "status": status,
            "positions_sha256": semantic_sha256(component.positions),
        }
        component_rows.append(component_row)
        global_counts["target_components"] += 1
        global_counts["target_sites"] += component.k
        global_counts["old_selected_sites"] += len(selected)
        global_counts["old_cap_excluded_sites"] += component.k - len(selected)
        global_counts["new_blocks"] += len(raw_result.blocks)
        global_counts["primary_active_sites_component_sum"] += len(active_sites)
        global_counts["exact_patterns"] += raw_result.total_pattern_count
        global_counts["raw_total_molecule_weight"] += int(raw_result.total_weight)
        global_counts["raw_retained_molecule_weight"] += int(
            raw_result.retained_weight
        )
        global_counts["raw_lost_molecule_weight"] += int(raw_result.lost_weight)
        global_counts["unavoidable_patterns"] += len(
            raw_result.unavoidable_constraint_ids
        )
        global_counts["weight_stable_components"] += int(weight_stable)
        semantic_components.append(
            {
                "component_id": component.component_id,
                "positions": component.positions,
                "cuts": raw_result.cut_indices,
                "status": status,
                "retained": raw_result.retained_constraint_ids,
                "lost": raw_result.lost_constraint_ids,
            }
        )
    partition_seconds = time.perf_counter() - stage_started

    expected_sites = tuple(
        sorted(site.site_index for component in components for site in component.sites)
    )
    checks = {
        "target_component_count_matches_expected": (
            args.expected_target_components is None
            or len(components) == args.expected_target_components
        ),
        "target_site_count_matches_expected": (
            args.expected_target_sites is None
            or len(expected_sites) == args.expected_target_sites
        ),
        "target_sites_assigned_once": tuple(sorted(all_sites_assigned))
        == expected_sites,
        "every_block_k_lte_max": all(
            int(row["k"]) <= args.max_block_size for row in block_rows
        ),
        "old_selected_plus_excluded_equals_target": (
            global_counts["old_selected_sites"]
            + global_counts["old_cap_excluded_sites"]
            == global_counts["target_sites"]
        ),
        "new_retained_sites_equals_target": (
            global_counts["target_sites"] == len(membership_rows)
        ),
        "constraint_molecule_mass_conserved": (
            global_counts["raw_total_molecule_weight"]
            == global_counts["raw_retained_molecule_weight"]
            + global_counts["raw_lost_molecule_weight"]
        ),
        "constraint_rows_equal_exact_patterns": (
            len(constraint_rows) == global_counts["exact_patterns"]
        ),
        "hp_ps_columns_nonempty_and_isolated": all(
            row["hp_family"] in PRIMARY_HP_FAMILIES and bool(row["phase_set"])
            for row in constraint_rows
        ),
    }
    all_pass = all(checks.values())
    if not all_pass:
        raise RuntimeError(
            "partition verification failed: "
            + ", ".join(key for key, value in checks.items() if not value)
        )

    stage_started = time.perf_counter()
    prefix = f"{args.dataset}.{args.chrom}"
    component_path = args.output_dir / f"{prefix}.legacy_components.tsv.gz"
    blocks_path = args.output_dir / f"{prefix}.blocks.tsv.gz"
    membership_path = args.output_dir / f"{prefix}.site_membership.tsv.gz"
    constraints_path = args.output_dir / f"{prefix}.cut_constraints.tsv.gz"

    component_fields = list(component_rows[0]) if component_rows else [
        "dataset",
        "chrom",
        "legacy_component_id",
        "legacy_group_ordinal",
        "start1",
        "end1",
        "span_bp",
        "pre_cap_k",
        "old_densest8_selected",
        "old_cap_excluded",
        "new_block_count",
        "new_site_retained",
        "primary_active_site_count",
        "primary_active_site_fraction",
        "exact_pattern_count",
        "raw_total_molecule_weight",
        "raw_retained_molecule_weight",
        "raw_lost_molecule_weight",
        "raw_retention_ratio",
        "old_densest8_retained_molecule_weight",
        "old_densest8_retention_ratio",
        "retained_exact_pattern_count",
        "lost_exact_pattern_count",
        "unavoidable_pattern_count",
        "old_densest8_retained_pattern_count",
        "weight_stable",
        "raw_cut_indices",
        "equal_pattern_cut_indices",
        "log1p_cut_indices",
        "cut_gap_sum_bp",
        "status",
        "positions_sha256",
    ]
    block_fields = list(block_rows[0]) if block_rows else [
        "dataset",
        "chrom",
        "legacy_component_id",
        "block_index",
        "block_id",
        "start1",
        "end1",
        "k",
        "positions",
        "component_status",
        "raw_retained_molecule_weight",
        "retained_exact_pattern_count",
        "primary_active_site_count",
    ]
    membership_fields = list(membership_rows[0]) if membership_rows else [
        "dataset",
        "chrom",
        "legacy_component_id",
        "site_index",
        "component_local_index",
        "pos1",
        "ref",
        "alt",
        "block_id",
        "primary_linkage_observed",
        "old_densest8_disposition",
    ]
    constraint_fields = list(constraint_rows[0]) if constraint_rows else [
        "dataset",
        "chrom",
        "legacy_component_id",
        "constraint_id",
        "unit_id",
        "hp_family",
        "phase_set",
        "positions",
        "call_codes",
        "n_fixed_ra",
        "span_sites",
        "molecule_weight",
        "disposition",
        "crossed_cut_count",
        "retained_block_index",
    ]
    _write_tsv_gz(component_path, component_fields, component_rows)
    _write_tsv_gz(blocks_path, block_fields, block_rows)
    _write_tsv_gz(membership_path, membership_fields, membership_rows)
    _write_tsv_gz(constraints_path, constraint_fields, constraint_rows)
    write_seconds = time.perf_counter() - stage_started

    receipt_path = args.output_dir / "receipt.json"
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": all_pass,
        "scope": {
            "dataset": args.dataset,
            "chrom": args.chrom,
            "site_catalog_sites": len(sites),
            "legacy_multilocus_components": n_multilocus,
            "legacy_singletons": n_singleton,
        },
        "parameters": {
            "legacy_gap_bp": args.legacy_gap_bp,
            "max_block_size": args.max_block_size,
            "primary_hp_families": sorted(PRIMARY_HP_FAMILIES),
            "require_exact_known_phase_set": True,
            "fixed_call_codes": sorted(FIXED_CALLS),
            "primary_weighting": "raw_molecule",
            "sensitivity_weightings": [
                "equal_pattern",
                "log1p_molecule",
            ],
            "tie_break": [
                "retained_weight_desc",
                "retained_pattern_count_desc",
                "block_count_asc",
                "cut_gap_sum_desc",
                "cut_indices_lexicographic_asc",
            ],
        },
        "counts": dict(sorted(global_counts.items())),
        "status_counts": dict(sorted(status_counts.items())),
        "pattern_diagnostics": pattern_diagnostics["counts"],
        "primary_rows_by_hp_family": pattern_diagnostics[
            "primary_rows_by_hp_family"
        ],
        "primary_hp_phase_set_unit_count": pattern_diagnostics[
            "primary_rows_by_hp_phase_set_count"
        ],
        "checks": checks,
        "timings_seconds": {
            "load_site_catalog": site_load_seconds,
            "load_and_aggregate_primary_patterns": pattern_load_seconds,
            "ordered_hypergraph_dp": partition_seconds,
            "write_outputs": write_seconds,
            "total_before_receipt": time.perf_counter() - started,
        },
        "inputs": {
            "site_catalog": output_identity(args.site_catalog),
            "molecule_calls": output_identity(args.molecule_calls),
        },
        "outputs": {
            "legacy_components": output_identity(component_path),
            "blocks": output_identity(blocks_path),
            "site_membership": output_identity(membership_path),
            "cut_constraints": output_identity(constraints_path),
        },
        "semantic_result_sha256": semantic_sha256(semantic_components),
        "command": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "claim_ceiling": (
            "optimal local disjoint k<=8 partition under observed read constraints; "
            "not a unique true evolutionary tree"
        ),
    }
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    digest_path = args.output_dir / "receipt.json.sha256"
    digest_path.write_text(
        f"{sha256_path(receipt_path)}  {receipt_path.name}\n",
        encoding="utf-8",
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--site-catalog", required=True, type=Path)
    parser.add_argument("--molecule-calls", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--legacy-gap-bp", type=int, default=50_000)
    parser.add_argument("--max-block-size", type=int, default=8)
    parser.add_argument("--expected-target-components", type=int)
    parser.add_argument("--expected-target-sites", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    receipt = run(args)
    summary = {
        "all_pass": receipt["all_pass"],
        "scope": receipt["scope"],
        "counts": receipt["counts"],
        "status_counts": receipt["status_counts"],
        "timings_seconds": receipt["timings_seconds"],
        "semantic_result_sha256": receipt["semantic_result_sha256"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
