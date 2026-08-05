#!/usr/bin/env python3
"""Build exact PS x HP read constraints and deterministic k<=12 partitions.

This pilot consumes the lossless extractor's ``site_catalog.tsv.gz``,
``site_component_membership.tsv.gz``, and ``molecule_sparse_calls.tsv.gz``.
It deliberately accepts only threshold 3, ``PRIMARY_PS_AWARE`` memberships,
``PS_HP1``/``PS_HP2`` linkage, and non-empty phase sets.  A schema 1.2 input
may also contain its legitimate threshold 1/2/5 rows; those are audited under
``excluded_membership_other_threshold`` and cannot enter threshold-3 units.
For the selected threshold, a malformed primary row fails closed while a true
non-primary row is counted and excluded.  Thus every inference unit is one
source component with exactly one phase-set and haplotype identity.

Cross-language normalized schemas (column order is contractual):

``units.tsv.gz``
    dataset, chrom, unit_id, component_id, linkage_basis, hp_family,
    phase_set, threshold, k, positions1

``constraints.tsv.gz``
    dataset, chrom, unit_id, constraint_id, hp_family, phase_set, positions1,
    call_codes, molecule_weight, pattern_count

``molecule_weight`` is a base-10 non-negative integer string representing raw
unique-molecule support.  ``pattern_count`` is always ``1`` because every row
is one aggregated fixed R/A pattern.  X/L/O/D/S observations never establish
linkage.  A read pattern becomes a constraint only when at least two fixed R/A
calls route to the exact same unit.

All eligible units are emitted, including k<=12 units, which are represented
by exactly one output block.  Larger units are partitioned with the existing
ordered-hypergraph dynamic program and ``max_block_size=12``.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
from typing import Iterable, Mapping, Sequence


RESEARCH_ROOT = Path(__file__).resolve().parents[2]
PARTITION_SCRIPT_DIR = (
    RESEARCH_ROOT
    / "20260718_k_gt8_read_supported_segmentation"
    / "scripts"
)
if str(PARTITION_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(PARTITION_SCRIPT_DIR))

from read_support_partition import (  # noqa: E402
    CUT,
    RETAINED,
    UNAVOIDABLE,
    Hyperedge,
    partition_ordered_hypergraph,
)


SCHEMA_NAME = "intersubmod.exact_ps_k12_partition"
SCHEMA_VERSION = "0.1.0"
REQUIRED_THRESHOLD = 3
REQUIRED_MAX_BLOCK_SIZE = 12
FIXED_CALLS = frozenset({"R", "A"})
NON_LINKING_CALLS = frozenset({"X", "L", "O", "D", "S"})
VALID_CALLS = FIXED_CALLS | NON_LINKING_CALLS
BASIS_TO_HP = {"PS_HP1": "1", "PS_HP2": "2"}

UNITS_FIELDS = (
    "dataset",
    "chrom",
    "unit_id",
    "component_id",
    "linkage_basis",
    "hp_family",
    "phase_set",
    "threshold",
    "k",
    "positions1",
)
CONSTRAINT_FIELDS = (
    "dataset",
    "chrom",
    "unit_id",
    "constraint_id",
    "hp_family",
    "phase_set",
    "positions1",
    "call_codes",
    "molecule_weight",
    "pattern_count",
)
BLOCK_FIELDS = (
    "dataset",
    "chrom",
    "unit_id",
    "component_id",
    "linkage_basis",
    "hp_family",
    "phase_set",
    "threshold",
    "block_id",
    "block_index",
    "start1",
    "end1",
    "k",
    "positions1",
    "retained_molecule_weight",
    "retained_pattern_count",
)
MEMBERSHIP_FIELDS = (
    "dataset",
    "chrom",
    "unit_id",
    "component_id",
    "linkage_basis",
    "hp_family",
    "phase_set",
    "threshold",
    "site_index",
    "pos1",
    "unit_local_index",
    "block_id",
    "block_index",
)
DISPOSITION_FIELDS = (
    "dataset",
    "chrom",
    "unit_id",
    "constraint_id",
    "hp_family",
    "phase_set",
    "positions1",
    "call_codes",
    "molecule_weight",
    "pattern_count",
    "disposition",
    "span_sites",
    "crossed_cut_count",
    "retained_block_index",
)


@dataclass(frozen=True)
class Site:
    site_index: int
    chrom: str
    pos1: int
    ref: str
    alt: str


@dataclass(frozen=True)
class Unit:
    unit_id: str
    component_id: str
    linkage_basis: str
    hp_family: str
    phase_set: str
    threshold: int
    sites: tuple[Site, ...]

    @property
    def positions(self) -> tuple[int, ...]:
        return tuple(site.pos1 for site in self.sites)


@dataclass(frozen=True)
class Constraint:
    unit_id: str
    constraint_id: str
    hp_family: str
    phase_set: str
    positions: tuple[int, ...]
    call_codes: str
    molecule_weight: int


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--site-catalog", required=True, type=Path)
    parser.add_argument("--site-component-membership", required=True, type=Path)
    parser.add_argument("--molecule-calls", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threshold", required=True, type=int)
    parser.add_argument("--max-block-size", required=True, type=int)
    return parser.parse_args()


def _require_fields(
    path: Path, fieldnames: Sequence[str] | None, required: Iterable[str]
) -> None:
    missing = sorted(set(required) - set(fieldnames or ()))
    if missing:
        raise ValueError(f"{path}: missing required columns: {','.join(missing)}")


def _read_tsv_gz(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_fields(path, reader.fieldnames, required)
        return [dict(row) for row in reader]


def _write_tsv_gz(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    payload = buffer.getvalue().encode("utf-8")
    with path.open("xb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            gzip_handle.write(payload)


def _canonical_json_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _semantic_id(prefix: str, value: object, length: int = 24) -> str:
    digest = hashlib.sha256(_canonical_json_bytes(value)).hexdigest()[:length]
    return f"{prefix}{digest}"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_identity(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": _sha256_file(path),
    }


def _parse_nonnegative_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be an integer: {value!r}") from exc
    if parsed < 0:
        raise ValueError(f"{label} must be non-negative: {parsed}")
    return parsed


def _parse_csv_ints(value: str, label: str) -> tuple[int, ...]:
    if value == "":
        return ()
    return tuple(
        _parse_nonnegative_int(token, label)
        for token in value.split(",")
    )


def _load_sites(path: Path, chrom: str) -> tuple[Site, ...]:
    rows = _read_tsv_gz(path, ("site_index", "chrom", "pos1", "ref", "alt"))
    sites = []
    for row in rows:
        if row["chrom"] != chrom:
            raise ValueError(
                f"site catalog chrom mismatch: expected {chrom!r}, got {row['chrom']!r}"
            )
        sites.append(
            Site(
                site_index=_parse_nonnegative_int(row["site_index"], "site_index"),
                chrom=row["chrom"],
                pos1=_parse_nonnegative_int(row["pos1"], "pos1"),
                ref=row["ref"],
                alt=row["alt"],
            )
        )
    sites.sort(key=lambda site: site.site_index)
    indices = tuple(site.site_index for site in sites)
    if len(indices) != len(set(indices)):
        raise ValueError("site catalog contains duplicate site_index values")
    if indices != tuple(range(len(indices))):
        raise ValueError("site catalog site_index values must be contiguous from zero")
    positions = tuple(site.pos1 for site in sites)
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise ValueError("site catalog positions must be strictly increasing by site_index")
    return tuple(sites)


def _load_units(
    path: Path,
    *,
    dataset: str,
    chrom: str,
    threshold: int,
    site_by_index: Mapping[int, Site],
    counts: Counter[str],
) -> tuple[Unit, ...]:
    required = (
        "dataset",
        "chrom",
        "linkage_basis",
        "phase_set",
        "phase_set_status",
        "inference_role",
        "threshold",
        "site_index",
        "pos1",
        "component_id",
    )
    rows = _read_tsv_gz(path, required)
    accepted: list[dict[str, object]] = []
    for row in rows:
        if row["dataset"] != dataset or row["chrom"] != chrom:
            raise ValueError("membership dataset/chrom does not match requested scope")
        row_threshold = _parse_nonnegative_int(
            row["threshold"], "membership threshold"
        )
        if row_threshold != threshold:
            counts["excluded_membership_other_threshold"] += 1
            continue
        if row["inference_role"] != "PRIMARY_PS_AWARE":
            counts["excluded_membership_nonprimary_or_missing_ps"] += 1
            continue
        basis = row["linkage_basis"]
        phase_set = row["phase_set"]
        primary_contract_errors = []
        if row["phase_set_status"] != "KNOWN_PS_PRIMARY":
            primary_contract_errors.append(
                f"phase_set_status={row['phase_set_status']!r}"
            )
        if basis not in BASIS_TO_HP:
            primary_contract_errors.append(f"linkage_basis={basis!r}")
        if phase_set == "":
            primary_contract_errors.append("phase_set=missing")
        if primary_contract_errors:
            raise ValueError(
                "PRIMARY_PS_AWARE membership violates exact primary contract: "
                + ";".join(primary_contract_errors)
            )
        component_id = row["component_id"]
        if component_id == "":
            raise ValueError("accepted membership component_id must not be empty")
        site_index = _parse_nonnegative_int(row["site_index"], "membership site_index")
        if site_index not in site_by_index:
            raise ValueError(f"membership references unknown site_index {site_index}")
        site = site_by_index[site_index]
        pos1 = _parse_nonnegative_int(row["pos1"], "membership pos1")
        if site.pos1 != pos1:
            raise ValueError(
                f"membership pos1 mismatch for site_index {site_index}: {pos1} != {site.pos1}"
            )
        accepted.append(
            {
                "component_id": component_id,
                "linkage_basis": basis,
                "hp_family": BASIS_TO_HP[basis],
                "phase_set": phase_set,
                "threshold": row_threshold,
                "site": site,
            }
        )

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in accepted:
        grouped[str(row["component_id"])].append(row)

    units = []
    seen_route_sites: dict[tuple[str, str, int], str] = {}
    for component_id, component_rows in sorted(grouped.items()):
        identities = {
            (
                str(row["linkage_basis"]),
                str(row["hp_family"]),
                str(row["phase_set"]),
                int(row["threshold"]),
            )
            for row in component_rows
        }
        if len(identities) != 1:
            raise ValueError(
                f"source component is not exact PS x HP: component_id={component_id!r}"
            )
        linkage_basis, hp_family, phase_set, unit_threshold = next(iter(identities))
        component_sites = tuple(
            sorted(
                (row["site"] for row in component_rows),
                key=lambda site: site.site_index,
            )
        )
        indices = tuple(site.site_index for site in component_sites)
        if len(indices) != len(set(indices)):
            raise ValueError(
                f"source component contains duplicate sites: component_id={component_id!r}"
            )
        unit_semantics = {
            "dataset": dataset,
            "chrom": chrom,
            "component_id": component_id,
            "linkage_basis": linkage_basis,
            "hp_family": hp_family,
            "phase_set": phase_set,
            "threshold": unit_threshold,
            "site_indices": indices,
        }
        unit_id = _semantic_id("U", unit_semantics)
        for site in component_sites:
            route_key = (hp_family, phase_set, site.site_index)
            previous = seen_route_sites.get(route_key)
            if previous is not None:
                raise ValueError(
                    "ambiguous exact PS x HP site membership: "
                    f"{route_key!r} belongs to {previous!r} and {unit_id!r}"
                )
            seen_route_sites[route_key] = unit_id
        units.append(
            Unit(
                unit_id=unit_id,
                component_id=component_id,
                linkage_basis=linkage_basis,
                hp_family=hp_family,
                phase_set=phase_set,
                threshold=unit_threshold,
                sites=component_sites,
            )
        )
    return tuple(sorted(units, key=lambda unit: unit.unit_id))


def _load_constraints(
    path: Path,
    *,
    dataset: str,
    chrom: str,
    units: Sequence[Unit],
    site_by_index: Mapping[int, Site],
    counts: Counter[str],
) -> tuple[Constraint, ...]:
    required = (
        "dataset",
        "chrom",
        "molecule_id",
        "hp_family",
        "phase_set",
        "site_indices",
        "positions1",
        "call_codes",
    )
    rows = _read_tsv_gz(path, required)
    rows.sort(
        key=lambda row: (
            row["molecule_id"],
            row["hp_family"],
            row["phase_set"],
            row["site_indices"],
            row["call_codes"],
        )
    )
    molecule_ids = tuple(row["molecule_id"] for row in rows)
    if any(molecule_id == "" for molecule_id in molecule_ids):
        raise ValueError("molecule_id must not be empty")
    if len(molecule_ids) != len(set(molecule_ids)):
        raise ValueError("molecule_sparse_calls contains duplicate molecule_id values")

    unit_by_id = {unit.unit_id: unit for unit in units}
    route: dict[tuple[str, str, int], str] = {}
    local_order: dict[str, dict[int, int]] = {}
    for unit in units:
        local_order[unit.unit_id] = {
            site.site_index: index for index, site in enumerate(unit.sites)
        }
        for site in unit.sites:
            route[(unit.hp_family, unit.phase_set, site.site_index)] = unit.unit_id

    patterns: Counter[tuple[str, tuple[int, ...], str]] = Counter()
    for row in rows:
        counts["molecule_rows_total"] += 1
        if row["dataset"] != dataset or row["chrom"] != chrom:
            raise ValueError("molecule dataset/chrom does not match requested scope")
        hp_family = row["hp_family"]
        phase_set = row["phase_set"]
        if hp_family not in {"1", "2"}:
            counts["excluded_molecule_nonprimary_hp"] += 1
            continue
        if phase_set == "":
            counts["excluded_molecule_missing_phase_set"] += 1
            continue
        indices = _parse_csv_ints(row["site_indices"], "molecule site_index")
        positions = _parse_csv_ints(row["positions1"], "molecule positions1")
        codes = row["call_codes"]
        if not (len(indices) == len(positions) == len(codes)):
            raise ValueError(
                f"molecule {row['molecule_id']!r}: "
                "site_indices/positions1/call_codes length mismatch"
            )
        if len(indices) != len(set(indices)):
            raise ValueError(f"molecule {row['molecule_id']!r}: duplicate site_index")
        unknown_codes = sorted(set(codes) - VALID_CALLS)
        if unknown_codes:
            raise ValueError(
                f"molecule {row['molecule_id']!r}: unknown call codes {unknown_codes}"
            )
        for site_index, position in zip(indices, positions):
            if site_index not in site_by_index:
                raise ValueError(
                    f"molecule {row['molecule_id']!r}: unknown site_index {site_index}"
                )
            expected_position = site_by_index[site_index].pos1
            if position != expected_position:
                raise ValueError(
                    f"molecule {row['molecule_id']!r}: positions1 mismatch for "
                    f"site_index {site_index}: {position} != {expected_position}"
                )
        fixed_seen = 0
        routed_by_unit: dict[str, list[tuple[int, str]]] = defaultdict(list)
        for site_index, code in zip(indices, codes):
            if code not in FIXED_CALLS:
                counts["non_linking_calls_ignored"] += 1
                continue
            fixed_seen += 1
            unit_id = route.get((hp_family, phase_set, site_index))
            if unit_id is None:
                counts["fixed_calls_not_routed"] += 1
                continue
            routed_by_unit[unit_id].append((site_index, code))
        if fixed_seen and not routed_by_unit:
            counts["molecule_rows_no_matching_unit"] += 1
        for unit_id, observations in sorted(routed_by_unit.items()):
            observations.sort(key=lambda item: local_order[unit_id][item[0]])
            if len(observations) < 2:
                counts["routed_single_fixed_call_no_constraint"] += 1
                continue
            positions = tuple(site_by_index[index].pos1 for index, _ in observations)
            pattern_codes = "".join(code for _, code in observations)
            patterns[(unit_id, positions, pattern_codes)] += 1
            counts["molecule_unit_constraints_contributed"] += 1

    constraints = []
    for (unit_id, positions, codes), weight in sorted(patterns.items()):
        unit = unit_by_id[unit_id]
        semantics = {
            "unit_id": unit_id,
            "positions1": positions,
            "call_codes": codes,
        }
        constraints.append(
            Constraint(
                unit_id=unit_id,
                constraint_id=_semantic_id("C", semantics),
                hp_family=unit.hp_family,
                phase_set=unit.phase_set,
                positions=positions,
                call_codes=codes,
                molecule_weight=weight,
            )
        )
    constraint_ids = tuple(constraint.constraint_id for constraint in constraints)
    if len(constraint_ids) != len(set(constraint_ids)):
        raise AssertionError("semantic constraint_id collision")
    return tuple(sorted(constraints, key=lambda item: (item.unit_id, item.constraint_id)))


def _decimal_string(value: Decimal) -> str:
    if value == value.to_integral_value():
        return str(value.to_integral_value())
    return format(value.normalize(), "f")


def _unit_rows(dataset: str, chrom: str, units: Sequence[Unit]) -> list[dict[str, object]]:
    return [
        {
            "dataset": dataset,
            "chrom": chrom,
            "unit_id": unit.unit_id,
            "component_id": unit.component_id,
            "linkage_basis": unit.linkage_basis,
            "hp_family": unit.hp_family,
            "phase_set": unit.phase_set,
            "threshold": str(unit.threshold),
            "k": str(len(unit.sites)),
            "positions1": ",".join(map(str, unit.positions)),
        }
        for unit in units
    ]


def _constraint_rows(
    dataset: str, chrom: str, constraints: Sequence[Constraint]
) -> list[dict[str, object]]:
    return [
        {
            "dataset": dataset,
            "chrom": chrom,
            "unit_id": constraint.unit_id,
            "constraint_id": constraint.constraint_id,
            "hp_family": constraint.hp_family,
            "phase_set": constraint.phase_set,
            "positions1": ",".join(map(str, constraint.positions)),
            "call_codes": constraint.call_codes,
            "molecule_weight": str(constraint.molecule_weight),
            "pattern_count": "1",
        }
        for constraint in constraints
    ]


def _partition_units(
    dataset: str,
    chrom: str,
    units: Sequence[Unit],
    constraints: Sequence[Constraint],
    max_block_size: int,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    constraints_by_unit: dict[str, list[Constraint]] = defaultdict(list)
    for constraint in constraints:
        constraints_by_unit[constraint.unit_id].append(constraint)
    block_rows = []
    membership_rows = []
    disposition_rows = []
    for unit in units:
        unit_constraints = tuple(
            sorted(constraints_by_unit.get(unit.unit_id, ()), key=lambda item: item.constraint_id)
        )
        constraint_by_id = {
            constraint.constraint_id: constraint for constraint in unit_constraints
        }
        result = partition_ordered_hypergraph(
            unit.positions,
            (
                Hyperedge(
                    constraint_id=constraint.constraint_id,
                    sites=constraint.positions,
                    weight=Decimal(constraint.molecule_weight),
                    pattern_count=1,
                )
                for constraint in unit_constraints
            ),
            max_block_size=max_block_size,
        )
        block_id_by_index = {}
        block_for_position = {}
        for block in result.blocks:
            block_number = block.block_index + 1
            block_id = f"{unit.unit_id}:B{block_number:04d}"
            block_id_by_index[block.block_index] = block_id
            for position in block.positions:
                block_for_position[position] = (block_id, block_number)
            block_rows.append(
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "unit_id": unit.unit_id,
                    "component_id": unit.component_id,
                    "linkage_basis": unit.linkage_basis,
                    "hp_family": unit.hp_family,
                    "phase_set": unit.phase_set,
                    "threshold": str(unit.threshold),
                    "block_id": block_id,
                    "block_index": str(block_number),
                    "start1": str(block.positions[0]),
                    "end1": str(block.positions[-1]),
                    "k": str(block.k),
                    "positions1": ",".join(map(str, block.positions)),
                    "retained_molecule_weight": _decimal_string(block.retained_weight),
                    "retained_pattern_count": str(block.retained_pattern_count),
                }
            )
        for local_index, site in enumerate(unit.sites, start=1):
            block_id, block_number = block_for_position[site.pos1]
            membership_rows.append(
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "unit_id": unit.unit_id,
                    "component_id": unit.component_id,
                    "linkage_basis": unit.linkage_basis,
                    "hp_family": unit.hp_family,
                    "phase_set": unit.phase_set,
                    "threshold": str(unit.threshold),
                    "site_index": str(site.site_index),
                    "pos1": str(site.pos1),
                    "unit_local_index": str(local_index),
                    "block_id": block_id,
                    "block_index": str(block_number),
                }
            )
        for disposition in result.dispositions:
            constraint = constraint_by_id[disposition.constraint_id]
            retained_block = (
                ""
                if disposition.block_index is None
                else str(disposition.block_index + 1)
            )
            disposition_rows.append(
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "unit_id": unit.unit_id,
                    "constraint_id": constraint.constraint_id,
                    "hp_family": constraint.hp_family,
                    "phase_set": constraint.phase_set,
                    "positions1": ",".join(map(str, constraint.positions)),
                    "call_codes": constraint.call_codes,
                    "molecule_weight": str(constraint.molecule_weight),
                    "pattern_count": "1",
                    "disposition": disposition.status,
                    "span_sites": str(disposition.span_sites),
                    "crossed_cut_count": str(disposition.crossed_cut_count),
                    "retained_block_index": retained_block,
                }
            )
    block_rows.sort(key=lambda row: (str(row["unit_id"]), int(row["block_index"])))
    membership_rows.sort(
        key=lambda row: (str(row["unit_id"]), int(row["unit_local_index"]))
    )
    disposition_rows.sort(
        key=lambda row: (str(row["unit_id"]), str(row["constraint_id"]))
    )
    return block_rows, membership_rows, disposition_rows


def _validate_outputs(
    units: Sequence[Unit],
    constraints: Sequence[Constraint],
    block_rows: Sequence[Mapping[str, object]],
    membership_rows: Sequence[Mapping[str, object]],
    disposition_rows: Sequence[Mapping[str, object]],
    max_block_size: int,
) -> tuple[dict[str, bool], dict[str, str]]:
    unit_by_id = {unit.unit_id: unit for unit in units}
    expected_sites = Counter(
        (unit.unit_id, site.site_index) for unit in units for site in unit.sites
    )
    observed_sites = Counter(
        (str(row["unit_id"]), int(row["site_index"])) for row in membership_rows
    )
    unit_sites_once = expected_sites == observed_sites and all(
        multiplicity == 1 for multiplicity in observed_sites.values()
    )
    blocks_bounded = all(int(row["k"]) <= max_block_size for row in block_rows)
    blocks_per_unit = Counter(str(row["unit_id"]) for row in block_rows)
    small_units_single = all(
        len(unit.sites) > max_block_size or blocks_per_unit[unit.unit_id] == 1
        for unit in units
    )

    cross_ps = 0
    cross_hp = 0
    for row in (*block_rows, *membership_rows, *disposition_rows):
        unit = unit_by_id[str(row["unit_id"])]
        cross_ps += str(row["phase_set"]) != unit.phase_set
        cross_hp += str(row["hp_family"]) != unit.hp_family
    for constraint in constraints:
        unit = unit_by_id[constraint.unit_id]
        cross_ps += constraint.phase_set != unit.phase_set
        cross_hp += constraint.hp_family != unit.hp_family

    constraint_ids = {constraint.constraint_id for constraint in constraints}
    disposition_ids = [str(row["constraint_id"]) for row in disposition_rows]
    disposition_once = (
        set(disposition_ids) == constraint_ids
        and len(disposition_ids) == len(set(disposition_ids))
    )
    mass = {RETAINED: 0, CUT: 0, UNAVOIDABLE: 0}
    for row in disposition_rows:
        status = str(row["disposition"])
        if status not in mass:
            raise AssertionError(f"unknown constraint disposition {status!r}")
        mass[status] += int(row["molecule_weight"])
    total = sum(constraint.molecule_weight for constraint in constraints)
    mass_conserved = total == sum(mass.values())
    checks = {
        "unit_sites_assigned_exactly_once": unit_sites_once,
        "block_k_at_most_12": blocks_bounded,
        "k_at_most_12_has_one_block": small_units_single,
        "cross_ps_zero": cross_ps == 0,
        "cross_hp_zero": cross_hp == 0,
        "constraint_ids_disposed_exactly_once": disposition_once,
        "constraint_mass_conserved": mass_conserved,
    }
    mass_receipt = {
        "total": str(total),
        "retained": str(mass[RETAINED]),
        "cut": str(mass[CUT]),
        "unavoidable": str(mass[UNAVOIDABLE]),
    }
    return checks, mass_receipt


def main() -> int:
    args = _parse_args()
    if args.threshold != REQUIRED_THRESHOLD:
        raise ValueError("threshold must be exactly 3")
    if args.max_block_size != REQUIRED_MAX_BLOCK_SIZE:
        raise ValueError("max-block-size must be exactly 12 for this pilot")
    if not args.dataset or not args.chrom:
        raise ValueError("dataset and chrom must not be empty")
    if args.output_dir.exists():
        if any(args.output_dir.iterdir()):
            raise FileExistsError(f"output directory is not empty: {args.output_dir}")
    else:
        args.output_dir.mkdir(parents=True)

    counts: Counter[str] = Counter()
    sites = _load_sites(args.site_catalog, args.chrom)
    site_by_index = {site.site_index: site for site in sites}
    units = _load_units(
        args.site_component_membership,
        dataset=args.dataset,
        chrom=args.chrom,
        threshold=args.threshold,
        site_by_index=site_by_index,
        counts=counts,
    )
    constraints = _load_constraints(
        args.molecule_calls,
        dataset=args.dataset,
        chrom=args.chrom,
        units=units,
        site_by_index=site_by_index,
        counts=counts,
    )

    unit_rows = _unit_rows(args.dataset, args.chrom, units)
    constraint_rows = _constraint_rows(args.dataset, args.chrom, constraints)
    units_path = args.output_dir / "units.tsv.gz"
    constraints_path = args.output_dir / "constraints.tsv.gz"
    _write_tsv_gz(units_path, UNITS_FIELDS, unit_rows)
    _write_tsv_gz(constraints_path, CONSTRAINT_FIELDS, constraint_rows)

    block_rows, membership_rows, disposition_rows = _partition_units(
        args.dataset,
        args.chrom,
        units,
        constraints,
        args.max_block_size,
    )
    blocks_path = args.output_dir / "blocks.tsv.gz"
    membership_path = args.output_dir / "membership.tsv.gz"
    dispositions_path = args.output_dir / "dispositions.tsv.gz"
    _write_tsv_gz(blocks_path, BLOCK_FIELDS, block_rows)
    _write_tsv_gz(membership_path, MEMBERSHIP_FIELDS, membership_rows)
    _write_tsv_gz(dispositions_path, DISPOSITION_FIELDS, disposition_rows)

    checks, constraint_mass = _validate_outputs(
        units,
        constraints,
        block_rows,
        membership_rows,
        disposition_rows,
        args.max_block_size,
    )
    all_pass = all(checks.values())
    if not all_pass:
        failed = sorted(name for name, passed in checks.items() if not passed)
        raise AssertionError(f"output validation failed: {','.join(failed)}")

    status_counts = Counter(str(row["disposition"]) for row in disposition_rows)
    counts.update(
        {
            "site_catalog_rows": len(sites),
            "eligible_units": len(units),
            "eligible_unit_sites": sum(len(unit.sites) for unit in units),
            "constraints": len(constraints),
            "blocks": len(block_rows),
            "membership_rows": len(membership_rows),
            "disposition_rows": len(disposition_rows),
            "retained_constraints": status_counts[RETAINED],
            "cut_constraints": status_counts[CUT],
            "unavoidable_constraints": status_counts[UNAVOIDABLE],
        }
    )
    for expected_counter in (
        "excluded_membership_nonprimary_or_missing_ps",
        "excluded_membership_other_threshold",
        "excluded_molecule_missing_phase_set",
        "excluded_molecule_nonprimary_hp",
        "molecule_rows_no_matching_unit",
        "molecule_rows_total",
        "non_linking_calls_ignored",
        "fixed_calls_not_routed",
        "routed_single_fixed_call_no_constraint",
        "molecule_unit_constraints_contributed",
    ):
        counts.setdefault(expected_counter, 0)

    semantic_payload = {
        "units": unit_rows,
        "constraints": constraint_rows,
        "blocks": block_rows,
        "membership": membership_rows,
        "dispositions": disposition_rows,
    }
    output_paths = (
        units_path,
        constraints_path,
        blocks_path,
        membership_path,
        dispositions_path,
    )
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": all_pass,
        "scope": {"dataset": args.dataset, "chrom": args.chrom},
        "parameters": {
            "threshold": args.threshold,
            "max_block_size": args.max_block_size,
            "accepted_inference_role": "PRIMARY_PS_AWARE",
            "accepted_linkage_basis": ["PS_HP1", "PS_HP2"],
            "phase_set_required": True,
            "fixed_linkage_calls": ["A", "R"],
            "non_linking_calls": ["D", "L", "O", "S", "X"],
        },
        "counts": dict(sorted(counts.items())),
        "checks": checks,
        "constraint_mass": constraint_mass,
        "inputs": {
            "site_catalog": _file_identity(args.site_catalog),
            "site_component_membership": _file_identity(
                args.site_component_membership
            ),
            "molecule_calls": _file_identity(args.molecule_calls),
        },
        "outputs": {path.name: _file_identity(path) for path in output_paths},
        "semantic_result_sha256": hashlib.sha256(
            _canonical_json_bytes(semantic_payload)
        ).hexdigest(),
        "claim_ceiling": (
            "Synthetic/pilot partition evidence only; no HCC1395 biological "
            "claim is authorized until the full chromosome and sample run."
        ),
    }
    receipt_path = args.output_dir / "receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir.resolve()),
                "all_pass": all_pass,
                "counts": receipt["counts"],
                "constraint_mass": constraint_mass,
                "semantic_result_sha256": receipt["semantic_result_sha256"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
