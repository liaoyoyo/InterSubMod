#!/usr/bin/env python3
"""Audit HP/PS-specific retention for observed k>8 partition constraints.

This is an independent, read-only red-team audit.  It does not reconstruct
unobserved HP/PS opportunities from sparse reads.  Its smallest statistical
unit is therefore:

    dataset x chromosome x legacy component x HP family x exact phase set

Every molecule weight in this audit is a molecule-by-component incidence, not
a unique molecule count.  A molecule can contribute to more than one legacy
component and is then counted once in each component-specific exact pattern.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping, Sequence


SCHEMA_NAME = "intersubmod.hp_ps_observed_constraint_retention_audit"
SCHEMA_VERSION = "1.0.0"
SOURCE_PARTITION_SCHEMA = "intersubmod.k_gt8_read_supported_segmentation"
SOURCE_RUN_SCHEMA = "intersubmod.hcc1395_full_k_gt8_segmentation.run_receipt"
SOURCE_SUCCESS_SCHEMA = (
    "intersubmod.hcc1395_full_k_gt8_segmentation.success_marker"
)
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))
PRIMARY_HP_FAMILIES = ("1", "2")
RETAINED = "retained"
CUT_LOST = "cut"
UNAVOIDABLE_PREFIX = "unavoidable_"

UNIT_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "hp_family",
    "phase_set",
    "unit_id",
    "component_k",
    "component_primary_active_site_count",
    "unit_active_site_count",
    "unit_active_site_fraction",
    "total_pattern_rows",
    "retained_pattern_rows",
    "cut_lost_pattern_rows",
    "unavoidable_pattern_rows",
    "nonretained_pattern_rows",
    "total_molecule_component_incidence_weight",
    "retained_molecule_component_incidence_weight",
    "cut_lost_molecule_component_incidence_weight",
    "unavoidable_molecule_component_incidence_weight",
    "nonretained_molecule_component_incidence_weight",
    "retention_ratio",
    "ratio_status",
    "support_stratum",
    "eligible_headline",
)

PAIR_FIELDS = (
    "dataset",
    "chrom",
    "legacy_component_id",
    "component_k",
    "hp1_phase_set_unit_count",
    "hp1_total_pattern_rows",
    "hp1_total_molecule_component_incidence_weight",
    "hp1_retained_molecule_component_incidence_weight",
    "hp1_retention_ratio",
    "hp2_phase_set_unit_count",
    "hp2_total_pattern_rows",
    "hp2_total_molecule_component_incidence_weight",
    "hp2_retained_molecule_component_incidence_weight",
    "hp2_retention_ratio",
    "hp1_minus_hp2_retention_delta",
    "absolute_retention_delta",
    "both_hp_headline_eligible",
)

SUMMARY_FIELDS = (
    "scope_type",
    "scope_id",
    "metric",
    "value",
    "unit",
    "denominator_note",
)


class AuditError(RuntimeError):
    """Raised when any audit contract is violated."""


@dataclass(frozen=True)
class SourcePartition:
    dataset: str
    chrom: str
    partition_dir: Path
    receipt: Mapping[str, Any]
    receipt_identity: Mapping[str, Any]
    legacy_components_path: Path
    membership_path: Path
    constraints_path: Path
    source_context: Mapping[str, Any]


@dataclass(frozen=True)
class ComponentInfo:
    dataset: str
    chrom: str
    component_id: str
    positions: tuple[int, ...]
    primary_active_positions: frozenset[int]

    @property
    def k(self) -> int:
        return len(self.positions)


@dataclass
class UnitAccumulator:
    dataset: str
    chrom: str
    component_id: str
    hp_family: str
    phase_set: str
    component_k: int
    component_primary_active_site_count: int
    active_positions: set[int] = field(default_factory=set)
    pattern_rows: Counter[str] = field(default_factory=Counter)
    molecule_weights: Counter[str] = field(default_factory=Counter)


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
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise AuditError(f"expected regular non-symlink file: {path}")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def strict_json_load(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot read JSON {path}: {exc}") from exc


def verify_sha256_sidecar(path: Path) -> str:
    sidecar = path.with_name(path.name + ".sha256")
    if not sidecar.is_file() or sidecar.is_symlink():
        raise AuditError(f"missing SHA-256 sidecar: {sidecar}")
    tokens = sidecar.read_text(encoding="utf-8").strip().split()
    if len(tokens) != 2 or tokens[1] != path.name:
        raise AuditError(f"malformed SHA-256 sidecar: {sidecar}")
    observed = sha256_path(path)
    if tokens[0] != observed:
        raise AuditError(f"SHA-256 mismatch: {path}")
    return observed


def verify_declared_identity(
    declared: Mapping[str, Any], *, expected_parent: Path | None = None
) -> Path:
    required = {"path", "size_bytes", "sha256"}
    if not required.issubset(declared):
        raise AuditError(f"declared identity lacks {sorted(required)}")
    path = Path(str(declared["path"]))
    if expected_parent is not None and path.resolve().parent != expected_parent.resolve():
        raise AuditError(f"declared output escapes expected directory: {path}")
    observed = file_identity(path)
    if observed != {
        "path": str(path.resolve()),
        "size_bytes": declared["size_bytes"],
        "sha256": declared["sha256"],
    }:
        raise AuditError(f"declared file identity drift: {path}")
    return path


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


def write_tsv_gz(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with _deterministic_gzip_writer(path) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def write_sha256_sidecar(path: Path) -> None:
    sidecar = path.with_name(path.name + ".sha256")
    sidecar.write_text(
        f"{sha256_path(path)}  {path.name}\n", encoding="utf-8"
    )


def parse_positive_integer(value: str, label: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise AuditError(f"{label} is not an integer: {value!r}") from exc
    if number <= 0:
        raise AuditError(f"{label} must be positive: {number}")
    return number


def parse_positions(value: str, label: str) -> tuple[int, ...]:
    if not value:
        raise AuditError(f"{label} positions are empty")
    try:
        positions = tuple(int(token) for token in value.split(","))
    except ValueError as exc:
        raise AuditError(f"{label} positions are malformed") from exc
    if len(positions) < 2:
        raise AuditError(f"{label} must contain at least two positions")
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise AuditError(f"{label} positions are not strictly increasing")
    return positions


def decimal_text(value: Decimal | None, places: int = 12) -> str:
    if value is None:
        return ""
    quantum = Decimal(1).scaleb(-places)
    with localcontext() as context:
        context.prec = max(40, len(value.as_tuple().digits) + places + 8)
        return format(value.quantize(quantum), "f")


def ratio(numerator: int, denominator: int) -> Decimal | None:
    if denominator == 0:
        return None
    with localcontext() as context:
        context.prec = 50
        return Decimal(numerator) / Decimal(denominator)


def quantile_type7(values: Sequence[Decimal], probability: Decimal) -> Decimal | None:
    """Return the R-7/NumPy-linear quantile with exact Decimal arithmetic."""
    if not values:
        return None
    if probability < 0 or probability > 1:
        raise ValueError("probability must be in [0, 1]")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    with localcontext() as context:
        context.prec = 50
        index = Decimal(len(ordered) - 1) * probability
        lower = int(index)
        fraction = index - Decimal(lower)
        if fraction == 0:
            return ordered[lower]
        return ordered[lower] + fraction * (ordered[lower + 1] - ordered[lower])


def support_stratum(total_weight: int) -> str:
    if total_weight <= 0:
        raise AuditError("observed unit has non-positive total weight")
    if total_weight <= 4:
        return "1-4"
    if total_weight <= 19:
        return "5-19"
    if total_weight <= 49:
        return "20-49"
    return ">=50"


def disposition_class(value: str) -> str:
    if value == RETAINED:
        return RETAINED
    if value == CUT_LOST:
        return CUT_LOST
    if value.startswith(UNAVOIDABLE_PREFIX):
        return "unavoidable"
    raise AuditError(f"unknown constraint disposition: {value!r}")


def _required_columns(
    fieldnames: Sequence[str] | None, required: set[str], path: Path
) -> None:
    if not fieldnames or not required.issubset(fieldnames):
        missing = sorted(required.difference(fieldnames or ()))
        raise AuditError(f"{path} missing required columns: {missing}")


def load_membership(
    source: SourcePartition,
) -> dict[str, ComponentInfo]:
    required = {
        "dataset",
        "chrom",
        "legacy_component_id",
        "component_local_index",
        "pos1",
        "primary_linkage_observed",
    }
    rows: MutableMapping[str, list[tuple[int, int, bool]]] = defaultdict(list)
    seen_sites: set[tuple[str, int]] = set()
    with gzip.open(
        source.membership_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _required_columns(reader.fieldnames, required, source.membership_path)
        for line_number, row in enumerate(reader, start=2):
            if row["dataset"] != source.dataset or row["chrom"] != source.chrom:
                raise AuditError(
                    f"membership scope mismatch at {source.membership_path}:{line_number}"
                )
            component_id = row["legacy_component_id"]
            local_index = parse_positive_integer(
                str(int(row["component_local_index"]) + 1),
                "component_local_index+1",
            ) - 1
            pos1 = parse_positive_integer(row["pos1"], "pos1")
            observed_text = row["primary_linkage_observed"].lower()
            if observed_text not in {"true", "false"}:
                raise AuditError(
                    f"invalid primary_linkage_observed at line {line_number}"
                )
            key = (component_id, local_index)
            if key in seen_sites:
                raise AuditError(f"duplicate membership key: {key}")
            seen_sites.add(key)
            rows[component_id].append(
                (local_index, pos1, observed_text == "true")
            )

    components: dict[str, ComponentInfo] = {}
    for component_id, values in rows.items():
        values.sort()
        indices = tuple(value[0] for value in values)
        if indices != tuple(range(len(values))):
            raise AuditError(
                f"component_local_index is not contiguous: {component_id}"
            )
        positions = tuple(value[1] for value in values)
        if any(right <= left for left, right in zip(positions, positions[1:])):
            raise AuditError(f"component positions not increasing: {component_id}")
        components[component_id] = ComponentInfo(
            dataset=source.dataset,
            chrom=source.chrom,
            component_id=component_id,
            positions=positions,
            primary_active_positions=frozenset(
                pos for _, pos, observed in values if observed
            ),
        )
    return components


def load_constraints(
    source: SourcePartition,
    components: Mapping[str, ComponentInfo],
) -> tuple[list[dict[str, Any]], dict[str, set[int]]]:
    required = {
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
        "molecule_weight",
        "disposition",
        "crossed_cut_count",
        "retained_block_index",
    }
    records: list[dict[str, Any]] = []
    active_by_component: dict[str, set[int]] = defaultdict(set)
    seen_constraints: set[str] = set()
    with gzip.open(
        source.constraints_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _required_columns(reader.fieldnames, required, source.constraints_path)
        for line_number, row in enumerate(reader, start=2):
            label = f"{source.constraints_path}:{line_number}"
            if row["dataset"] != source.dataset or row["chrom"] != source.chrom:
                raise AuditError(f"constraint scope mismatch at {label}")
            component_id = row["legacy_component_id"]
            component = components.get(component_id)
            if component is None:
                raise AuditError(f"constraint references unknown component at {label}")
            constraint_id = row["constraint_id"]
            if not constraint_id or constraint_id in seen_constraints:
                raise AuditError(f"duplicate/empty constraint_id at {label}")
            seen_constraints.add(constraint_id)
            hp_family = row["hp_family"]
            phase_set = row["phase_set"]
            if hp_family not in PRIMARY_HP_FAMILIES:
                raise AuditError(f"non-primary HP family at {label}")
            if not phase_set or phase_set.upper() in {"MISSING", "__MISSING__"}:
                raise AuditError(f"missing phase set in primary unit at {label}")
            expected_unit = f"HP{hp_family}|PS{phase_set}"
            if row["unit_id"] != expected_unit:
                raise AuditError(f"unit_id drift at {label}")
            positions = parse_positions(row["positions"], label)
            if not set(positions).issubset(component.positions):
                raise AuditError(f"constraint positions escape component at {label}")
            codes = row["call_codes"]
            if len(codes) != len(positions) or set(codes).difference({"R", "A"}):
                raise AuditError(f"constraint call_codes are not fixed R/A at {label}")
            if parse_positive_integer(row["n_fixed_ra"], "n_fixed_ra") != len(
                positions
            ):
                raise AuditError(f"n_fixed_ra mismatch at {label}")
            weight = parse_positive_integer(row["molecule_weight"], "molecule_weight")
            category = disposition_class(row["disposition"])
            crossed = int(row["crossed_cut_count"])
            if crossed < 0:
                raise AuditError(f"negative crossed_cut_count at {label}")
            block_index = row["retained_block_index"]
            if category == RETAINED:
                if crossed != 0 or not block_index:
                    raise AuditError(f"invalid retained disposition at {label}")
                parse_positive_integer(block_index, "retained_block_index")
            elif block_index:
                raise AuditError(f"non-retained constraint has block index at {label}")
            active_by_component[component_id].update(positions)
            records.append(
                {
                    "dataset": source.dataset,
                    "chrom": source.chrom,
                    "legacy_component_id": component_id,
                    "constraint_id": constraint_id,
                    "hp_family": hp_family,
                    "phase_set": phase_set,
                    "positions": positions,
                    "weight": weight,
                    "category": category,
                }
            )
    return records, active_by_component


def load_legacy_component_metrics(
    source: SourcePartition,
) -> dict[str, dict[str, int]]:
    required = {
        "dataset",
        "chrom",
        "legacy_component_id",
        "pre_cap_k",
        "primary_active_site_count",
        "exact_pattern_count",
        "raw_total_molecule_weight",
        "raw_retained_molecule_weight",
        "raw_lost_molecule_weight",
        "retained_exact_pattern_count",
        "lost_exact_pattern_count",
        "unavoidable_pattern_count",
    }
    metrics: dict[str, dict[str, int]] = {}
    with gzip.open(
        source.legacy_components_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _required_columns(
            reader.fieldnames, required, source.legacy_components_path
        )
        for line_number, row in enumerate(reader, start=2):
            label = f"{source.legacy_components_path}:{line_number}"
            if row["dataset"] != source.dataset or row["chrom"] != source.chrom:
                raise AuditError(f"legacy component scope mismatch at {label}")
            component_id = row["legacy_component_id"]
            if not component_id or component_id in metrics:
                raise AuditError(f"duplicate/empty legacy component at {label}")
            parsed = {}
            for name in required.difference(
                {"dataset", "chrom", "legacy_component_id"}
            ):
                try:
                    value = int(row[name])
                except ValueError as exc:
                    raise AuditError(
                        f"non-integer legacy metric {name} at {label}"
                    ) from exc
                if value < 0:
                    raise AuditError(
                        f"negative legacy metric {name} at {label}"
                    )
                parsed[name] = value
            metrics[component_id] = parsed
    return metrics


def aggregate_units(
    components: Mapping[str, ComponentInfo],
    constraints: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    accumulators: dict[tuple[str, str, str, str, str], UnitAccumulator] = {}
    for record in constraints:
        component = components[str(record["legacy_component_id"])]
        key = (
            str(record["dataset"]),
            str(record["chrom"]),
            component.component_id,
            str(record["hp_family"]),
            str(record["phase_set"]),
        )
        accumulator = accumulators.get(key)
        if accumulator is None:
            accumulator = UnitAccumulator(
                dataset=key[0],
                chrom=key[1],
                component_id=key[2],
                hp_family=key[3],
                phase_set=key[4],
                component_k=component.k,
                component_primary_active_site_count=len(
                    component.primary_active_positions
                ),
            )
            accumulators[key] = accumulator
        category = str(record["category"])
        accumulator.pattern_rows[category] += 1
        accumulator.molecule_weights[category] += int(record["weight"])
        accumulator.active_positions.update(record["positions"])

    rows: list[dict[str, Any]] = []
    for key in sorted(accumulators):
        item = accumulators[key]
        retained_patterns = item.pattern_rows[RETAINED]
        cut_patterns = item.pattern_rows[CUT_LOST]
        unavoidable_patterns = item.pattern_rows["unavoidable"]
        total_patterns = retained_patterns + cut_patterns + unavoidable_patterns
        retained_weight = item.molecule_weights[RETAINED]
        cut_weight = item.molecule_weights[CUT_LOST]
        unavoidable_weight = item.molecule_weights["unavoidable"]
        total_weight = retained_weight + cut_weight + unavoidable_weight
        if total_patterns <= 0 or total_weight <= 0:
            raise AuditError(f"empty observed unit: {key}")
        retention = ratio(retained_weight, total_weight)
        rows.append(
            {
                "dataset": item.dataset,
                "chrom": item.chrom,
                "legacy_component_id": item.component_id,
                "hp_family": item.hp_family,
                "phase_set": item.phase_set,
                "unit_id": f"HP{item.hp_family}|PS{item.phase_set}",
                "component_k": item.component_k,
                "component_primary_active_site_count": (
                    item.component_primary_active_site_count
                ),
                "unit_active_site_count": len(item.active_positions),
                "unit_active_site_fraction": decimal_text(
                    ratio(len(item.active_positions), item.component_k)
                ),
                "total_pattern_rows": total_patterns,
                "retained_pattern_rows": retained_patterns,
                "cut_lost_pattern_rows": cut_patterns,
                "unavoidable_pattern_rows": unavoidable_patterns,
                "nonretained_pattern_rows": cut_patterns + unavoidable_patterns,
                "total_molecule_component_incidence_weight": total_weight,
                "retained_molecule_component_incidence_weight": retained_weight,
                "cut_lost_molecule_component_incidence_weight": cut_weight,
                "unavoidable_molecule_component_incidence_weight": (
                    unavoidable_weight
                ),
                "nonretained_molecule_component_incidence_weight": (
                    cut_weight + unavoidable_weight
                ),
                "retention_ratio": decimal_text(retention),
                "ratio_status": "OBSERVED_CONSTRAINT_DENOMINATOR",
                "support_stratum": support_stratum(total_weight),
                "eligible_headline": str(
                    total_weight >= 20 and total_patterns >= 5
                ).lower(),
            }
        )
    return rows


def aggregate_component_pairs(
    unit_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], dict[str, Counter[str]]] = defaultdict(
        lambda: {"1": Counter(), "2": Counter()}
    )
    component_k: dict[tuple[str, str, str], int] = {}
    for row in unit_rows:
        key = (
            str(row["dataset"]),
            str(row["chrom"]),
            str(row["legacy_component_id"]),
        )
        hp = str(row["hp_family"])
        component_k[key] = int(row["component_k"])
        counts = grouped[key][hp]
        counts["phase_sets"] += 1
        counts["patterns"] += int(row["total_pattern_rows"])
        counts["total_weight"] += int(
            row["total_molecule_component_incidence_weight"]
        )
        counts["retained_weight"] += int(
            row["retained_molecule_component_incidence_weight"]
        )

    rows: list[dict[str, Any]] = []
    for key in sorted(grouped):
        hp1 = grouped[key]["1"]
        hp2 = grouped[key]["2"]
        if hp1["phase_sets"] == 0 or hp2["phase_sets"] == 0:
            continue
        hp1_ratio = ratio(hp1["retained_weight"], hp1["total_weight"])
        hp2_ratio = ratio(hp2["retained_weight"], hp2["total_weight"])
        assert hp1_ratio is not None and hp2_ratio is not None
        delta = hp1_ratio - hp2_ratio
        rows.append(
            {
                "dataset": key[0],
                "chrom": key[1],
                "legacy_component_id": key[2],
                "component_k": component_k[key],
                "hp1_phase_set_unit_count": hp1["phase_sets"],
                "hp1_total_pattern_rows": hp1["patterns"],
                "hp1_total_molecule_component_incidence_weight": hp1[
                    "total_weight"
                ],
                "hp1_retained_molecule_component_incidence_weight": hp1[
                    "retained_weight"
                ],
                "hp1_retention_ratio": decimal_text(hp1_ratio),
                "hp2_phase_set_unit_count": hp2["phase_sets"],
                "hp2_total_pattern_rows": hp2["patterns"],
                "hp2_total_molecule_component_incidence_weight": hp2[
                    "total_weight"
                ],
                "hp2_retained_molecule_component_incidence_weight": hp2[
                    "retained_weight"
                ],
                "hp2_retention_ratio": decimal_text(hp2_ratio),
                "hp1_minus_hp2_retention_delta": decimal_text(delta),
                "absolute_retention_delta": decimal_text(abs(delta)),
                "both_hp_headline_eligible": str(
                    hp1["total_weight"] >= 20
                    and hp1["patterns"] >= 5
                    and hp2["total_weight"] >= 20
                    and hp2["patterns"] >= 5
                ).lower(),
            }
        )
    return rows


def exact_phase_set_pair_deltas(
    unit_rows: Sequence[Mapping[str, Any]],
) -> list[Decimal]:
    units: dict[tuple[str, str, str, str], dict[str, Decimal]] = defaultdict(dict)
    for row in unit_rows:
        key = (
            str(row["dataset"]),
            str(row["chrom"]),
            str(row["legacy_component_id"]),
            str(row["phase_set"]),
        )
        units[key][str(row["hp_family"])] = Decimal(str(row["retention_ratio"]))
    return [
        values["1"] - values["2"]
        for _, values in sorted(units.items())
        if "1" in values and "2" in values
    ]


def ratio_distribution(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    ratios = [Decimal(str(row["retention_ratio"])) for row in rows]
    quantiles = {
        label: decimal_text(quantile_type7(ratios, probability))
        for label, probability in (
            ("min", Decimal("0")),
            ("p05", Decimal("0.05")),
            ("p25", Decimal("0.25")),
            ("median", Decimal("0.5")),
            ("p75", Decimal("0.75")),
            ("p95", Decimal("0.95")),
            ("max", Decimal("1")),
        )
    }
    cumulative = {
        "lt_0_5": sum(value < Decimal("0.5") for value in ratios),
        "lt_0_8": sum(value < Decimal("0.8") for value in ratios),
        "lt_0_9": sum(value < Decimal("0.9") for value in ratios),
        "lt_0_95": sum(value < Decimal("0.95") for value in ratios),
        "eq_1": sum(value == Decimal("1") for value in ratios),
    }
    exclusive = {
        "lt_0_5": sum(value < Decimal("0.5") for value in ratios),
        "ge_0_5_lt_0_8": sum(
            Decimal("0.5") <= value < Decimal("0.8") for value in ratios
        ),
        "ge_0_8_lt_0_9": sum(
            Decimal("0.8") <= value < Decimal("0.9") for value in ratios
        ),
        "ge_0_9_lt_0_95": sum(
            Decimal("0.9") <= value < Decimal("0.95") for value in ratios
        ),
        "ge_0_95_lt_1": sum(
            Decimal("0.95") <= value < Decimal("1") for value in ratios
        ),
        "eq_1": sum(value == Decimal("1") for value in ratios),
    }
    return {
        "n_units": len(rows),
        "quantile_method": "type7_linear_interpolation_unweighted_units",
        "quantiles": quantiles,
        "cumulative_threshold_counts": cumulative,
        "exclusive_bucket_counts": exclusive,
    }


def subset_metrics(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    eligible = [row for row in rows if row["eligible_headline"] == "true"]
    totals = Counter()
    for row in rows:
        for name in (
            "total_pattern_rows",
            "retained_pattern_rows",
            "cut_lost_pattern_rows",
            "unavoidable_pattern_rows",
            "total_molecule_component_incidence_weight",
            "retained_molecule_component_incidence_weight",
            "cut_lost_molecule_component_incidence_weight",
            "unavoidable_molecule_component_incidence_weight",
        ):
            totals[name] += int(row[name])
    total_weight = totals["total_molecule_component_incidence_weight"]
    return {
        "observed_constraint_units": len(rows),
        "eligible_headline_units": len(eligible),
        "pattern_rows": {
            "total": totals["total_pattern_rows"],
            "retained": totals["retained_pattern_rows"],
            "cut_lost": totals["cut_lost_pattern_rows"],
            "unavoidable": totals["unavoidable_pattern_rows"],
        },
        "molecule_component_incidences": {
            "total": total_weight,
            "retained": totals[
                "retained_molecule_component_incidence_weight"
            ],
            "cut_lost": totals[
                "cut_lost_molecule_component_incidence_weight"
            ],
            "unavoidable": totals[
                "unavoidable_molecule_component_incidence_weight"
            ],
            "weighted_retention_ratio": decimal_text(
                ratio(
                    totals["retained_molecule_component_incidence_weight"],
                    total_weight,
                )
            ),
        },
        "retention_distribution": ratio_distribution(rows),
    }


def summarize_rows(
    unit_rows: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
    *,
    source_components: Mapping[str, set[str]],
    skipped_chromosomes: Sequence[str],
) -> dict[str, Any]:
    component_with_units = {
        f"{row['chrom']}|{row['legacy_component_id']}" for row in unit_rows
    }
    component_with_hp1 = {
        f"{row['chrom']}|{row['legacy_component_id']}"
        for row in unit_rows
        if str(row["hp_family"]) == "1"
    }
    component_with_hp2 = {
        f"{row['chrom']}|{row['legacy_component_id']}"
        for row in unit_rows
        if str(row["hp_family"]) == "2"
    }
    all_components = {
        f"{chrom}|{component}"
        for chrom, components in source_components.items()
        for component in components
    }
    component_coverage_by_chromosome = {}
    for chrom in AUTOSOMES:
        chrom_components = {
            f"{chrom}|{component}"
            for component in source_components.get(chrom, set())
        }
        if not chrom_components:
            continue
        hp1 = chrom_components.intersection(component_with_hp1)
        hp2 = chrom_components.intersection(component_with_hp2)
        component_coverage_by_chromosome[chrom] = {
            "components_in_partition_scope": len(chrom_components),
            "components_with_any_observed_unit": len(hp1 | hp2),
            "components_hp1_only": len(hp1 - hp2),
            "components_hp2_only": len(hp2 - hp1),
            "components_hp1_and_hp2": len(hp1 & hp2),
            "components_without_observed_unit": len(
                chrom_components.difference(hp1 | hp2)
            ),
        }
    component_hp_coverage = {
        "components_in_partition_scope": len(all_components),
        "components_with_any_observed_unit": len(
            component_with_hp1 | component_with_hp2
        ),
        "components_hp1_only": len(component_with_hp1 - component_with_hp2),
        "components_hp2_only": len(component_with_hp2 - component_with_hp1),
        "components_hp1_and_hp2": len(component_with_hp1 & component_with_hp2),
        "components_without_observed_unit": len(
            all_components.difference(component_with_hp1 | component_with_hp2)
        ),
        "by_chromosome": component_coverage_by_chromosome,
    }
    totals = Counter()
    for row in unit_rows:
        for name in (
            "total_pattern_rows",
            "retained_pattern_rows",
            "cut_lost_pattern_rows",
            "unavoidable_pattern_rows",
            "nonretained_pattern_rows",
            "total_molecule_component_incidence_weight",
            "retained_molecule_component_incidence_weight",
            "cut_lost_molecule_component_incidence_weight",
            "unavoidable_molecule_component_incidence_weight",
            "nonretained_molecule_component_incidence_weight",
        ):
            totals[name] += int(row[name])
    total_weight = totals["total_molecule_component_incidence_weight"]
    retained_weight = totals["retained_molecule_component_incidence_weight"]
    eligible = [row for row in unit_rows if row["eligible_headline"] == "true"]
    hp_counts = Counter(str(row["hp_family"]) for row in unit_rows)
    support_counts = Counter(str(row["support_stratum"]) for row in unit_rows)
    exact_ps_deltas = exact_phase_set_pair_deltas(unit_rows)
    pair_abs_deltas = [
        Decimal(str(row["absolute_retention_delta"])) for row in pair_rows
    ]
    pair_eligible = [
        row for row in pair_rows if row["both_hp_headline_eligible"] == "true"
    ]
    by_chromosome = {
        chrom: subset_metrics(
            [row for row in unit_rows if str(row["chrom"]) == chrom]
        )
        for chrom in AUTOSOMES
        if any(str(row["chrom"]) == chrom for row in unit_rows)
    }
    by_hp_family = {
        f"HP{hp}": subset_metrics(
            [row for row in unit_rows if str(row["hp_family"]) == hp]
        )
        for hp in PRIMARY_HP_FAMILIES
    }
    by_support_stratum = {
        stratum: subset_metrics(
            [
                row
                for row in unit_rows
                if str(row["support_stratum"]) == stratum
            ]
        )
        for stratum in ("1-4", "5-19", "20-49", ">=50")
    }

    def delta_quantiles(values: Sequence[Decimal]) -> dict[str, str]:
        return {
            label: decimal_text(quantile_type7(values, probability))
            for label, probability in (
                ("min", Decimal("0")),
                ("p25", Decimal("0.25")),
                ("median", Decimal("0.5")),
                ("p75", Decimal("0.75")),
                ("p95", Decimal("0.95")),
                ("max", Decimal("1")),
            )
        }

    return {
        "scope_contract": {
            "unit_grain": (
                "dataset x chromosome x legacy_component_id x hp_family x "
                "exact known phase_set"
            ),
            "scope_ceiling": "observed_constraint_units_only",
            "aggregation_weight": "molecule_x_component_incidence",
            "unobserved_opportunity_policy": (
                "not reconstructed; no synthetic zero/one retention ratios"
            ),
            "skipped_zero_target_chromosomes": sorted(skipped_chromosomes),
        },
        "counts": {
            "observed_constraint_units": len(unit_rows),
            "components_in_partition_scope": len(all_components),
            "components_with_observed_constraint_units": len(component_with_units),
            "components_without_observed_constraint_units": len(
                all_components.difference(component_with_units)
            ),
            "eligible_headline_units": len(eligible),
            "support_ge_20_units": sum(
                int(row["total_molecule_component_incidence_weight"]) >= 20
                for row in unit_rows
            ),
            "support_ge_50_units": sum(
                int(row["total_molecule_component_incidence_weight"]) >= 50
                for row in unit_rows
            ),
            "hp1_units": hp_counts["1"],
            "hp2_units": hp_counts["2"],
            "hp1_hp2_paired_components": len(pair_rows),
            "hp1_hp2_paired_components_both_eligible": len(pair_eligible),
            "hp1_hp2_exact_phase_set_pairs": len(exact_ps_deltas),
        },
        "component_hp_unit_coverage": component_hp_coverage,
        "support_strata_disjoint": {
            name: support_counts[name]
            for name in ("1-4", "5-19", "20-49", ">=50")
        },
        "by_chromosome": by_chromosome,
        "by_hp_family": by_hp_family,
        "by_support_stratum": by_support_stratum,
        "molecule_component_incidence_totals": {
            "total": total_weight,
            "retained": retained_weight,
            "cut_lost": totals[
                "cut_lost_molecule_component_incidence_weight"
            ],
            "unavoidable": totals[
                "unavoidable_molecule_component_incidence_weight"
            ],
            "nonretained": totals[
                "nonretained_molecule_component_incidence_weight"
            ],
            "weighted_retention_ratio": decimal_text(
                ratio(retained_weight, total_weight)
            ),
        },
        "pattern_row_totals": {
            "total": totals["total_pattern_rows"],
            "retained": totals["retained_pattern_rows"],
            "cut_lost": totals["cut_lost_pattern_rows"],
            "unavoidable": totals["unavoidable_pattern_rows"],
            "nonretained": totals["nonretained_pattern_rows"],
        },
        "retention_distribution_all_observed_units": ratio_distribution(
            unit_rows
        ),
        "retention_distribution_eligible_headline_units": ratio_distribution(
            eligible
        ),
        "hp1_hp2_component_paired": {
            "pair_grain": (
                "dataset x chromosome x legacy_component_id; exact phase-set "
                "units are summed separately within HP1 and HP2"
            ),
            "n_pairs": len(pair_rows),
            "n_both_headline_eligible": len(pair_eligible),
            "absolute_delta_quantiles": delta_quantiles(pair_abs_deltas),
            "absolute_delta_lt_0_05": sum(
                value < Decimal("0.05") for value in pair_abs_deltas
            ),
            "absolute_delta_lt_0_10": sum(
                value < Decimal("0.10") for value in pair_abs_deltas
            ),
            "absolute_delta_ge_0_25": sum(
                value >= Decimal("0.25") for value in pair_abs_deltas
            ),
        },
        "hp1_hp2_exact_phase_set_paired": {
            "pair_grain": (
                "dataset x chromosome x legacy_component_id x exact phase_set"
            ),
            "n_pairs": len(exact_ps_deltas),
            "signed_delta_quantiles": delta_quantiles(exact_ps_deltas),
        },
    }


def build_summary_tsv(summary: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def add(
        scope_type: str,
        scope_id: str,
        metric: str,
        value: Any,
        unit: str,
        note: str,
    ) -> None:
        rows.append(
            {
                "scope_type": scope_type,
                "scope_id": scope_id,
                "metric": metric,
                "value": value,
                "unit": unit,
                "denominator_note": note,
            }
        )

    for metric, value in summary["counts"].items():
        add("overall", "observed_constraint_units", metric, value, "count", "")
    for metric, value in summary["component_hp_unit_coverage"].items():
        if metric == "by_chromosome":
            continue
        add(
            "component_coverage",
            "overall",
            metric,
            value,
            "count",
            "k>8 legacy components; observed constraint units only",
        )
    for chrom, metrics in summary[
        "component_hp_unit_coverage"
    ]["by_chromosome"].items():
        for metric, value in metrics.items():
            add(
                "component_coverage",
                chrom,
                metric,
                value,
                "count",
                "k>8 legacy components; observed constraint units only",
            )
    for metric, value in summary["support_strata_disjoint"].items():
        add(
            "support_stratum",
            metric,
            "observed_constraint_units",
            value,
            "count",
            "disjoint total molecule-by-component incidence stratum",
        )
    for scope_type, sections in (
        ("chromosome", summary["by_chromosome"]),
        ("hp_family", summary["by_hp_family"]),
        ("support_stratum_metrics", summary["by_support_stratum"]),
    ):
        for scope_id, section in sections.items():
            add(
                scope_type,
                scope_id,
                "observed_constraint_units",
                section["observed_constraint_units"],
                "count",
                "observed HP/PS units only",
            )
            add(
                scope_type,
                scope_id,
                "eligible_headline_units",
                section["eligible_headline_units"],
                "count",
                "weight >=20 and patterns >=5",
            )
            for metric, value in section[
                "molecule_component_incidences"
            ].items():
                add(
                    scope_type,
                    scope_id,
                    metric,
                    value,
                    "ratio"
                    if metric == "weighted_retention_ratio"
                    else "incidence_weight",
                    "molecule-by-component incidences",
                )
            for metric, value in section[
                "retention_distribution"
            ]["quantiles"].items():
                add(
                    scope_type,
                    scope_id,
                    f"retention_{metric}",
                    value,
                    "ratio",
                    "type7 unweighted observed HP/PS units",
                )
    for metric, value in summary["molecule_component_incidence_totals"].items():
        add(
            "overall",
            "molecule_x_component_incidence",
            metric,
            value,
            "ratio" if metric.endswith("ratio") else "incidence_weight",
            "not a unique-molecule count",
        )
    for metric, value in summary["pattern_row_totals"].items():
        add(
            "overall",
            "exact_pattern_rows",
            metric,
            value,
            "pattern_rows",
            "",
        )
    for scope_name in (
        "retention_distribution_all_observed_units",
        "retention_distribution_eligible_headline_units",
    ):
        distribution = summary[scope_name]
        add(
            "unit_distribution",
            scope_name,
            "n_units",
            distribution["n_units"],
            "count",
            "unweighted observed HP/PS units",
        )
        for metric, value in distribution["quantiles"].items():
            add(
                "unit_distribution",
                scope_name,
                f"retention_{metric}",
                value,
                "ratio",
                distribution["quantile_method"],
            )
        for metric, value in distribution["cumulative_threshold_counts"].items():
            add(
                "unit_distribution",
                scope_name,
                f"cumulative_{metric}",
                value,
                "count",
                "overlapping cumulative threshold count",
            )
    for metric, value in summary["hp1_hp2_component_paired"].items():
        if isinstance(value, (str, dict)):
            continue
        add(
            "hp_pair",
            "component_aggregated",
            metric,
            value,
            "count",
            "HP1/HP2 paired within the same legacy component",
        )
    for metric, value in summary[
        "hp1_hp2_component_paired"
    ]["absolute_delta_quantiles"].items():
        add(
            "hp_pair",
            "component_aggregated",
            f"absolute_delta_{metric}",
            value,
            "ratio_delta",
            "type7 unweighted paired components",
        )
    return rows


def verify_partition_source(
    partition_dir: Path, source_context: Mapping[str, Any]
) -> SourcePartition:
    if not partition_dir.is_dir() or partition_dir.is_symlink():
        raise AuditError(f"partition directory is invalid: {partition_dir}")
    receipt_path = partition_dir / "receipt.json"
    receipt_sha = verify_sha256_sidecar(receipt_path)
    receipt = strict_json_load(receipt_path)
    if not isinstance(receipt, Mapping):
        raise AuditError(f"partition receipt is not an object: {receipt_path}")
    if receipt.get("schema_name") != SOURCE_PARTITION_SCHEMA:
        raise AuditError(f"unexpected partition schema: {receipt_path}")
    if receipt.get("all_pass") is not True:
        raise AuditError(f"partition receipt is not all_pass: {receipt_path}")
    required_source_checks = {
        "constraint_molecule_mass_conserved",
        "constraint_rows_equal_exact_patterns",
        "hp_ps_columns_nonempty_and_isolated",
        "target_component_count_matches_expected",
        "target_site_count_matches_expected",
        "target_sites_assigned_once",
    }
    source_check_values = receipt.get("checks", {})
    if not isinstance(source_check_values, Mapping) or not required_source_checks.issubset(
        source_check_values
    ):
        raise AuditError(f"partition source checks are incomplete: {receipt_path}")
    if not all(source_check_values.values()):
        raise AuditError(f"partition source check is false: {receipt_path}")
    scope = receipt.get("scope")
    parameters = receipt.get("parameters")
    counts = receipt.get("counts")
    outputs = receipt.get("outputs")
    inputs = receipt.get("inputs")
    if not all(
        isinstance(value, Mapping)
        for value in (scope, parameters, counts, outputs, inputs)
    ):
        raise AuditError(f"partition receipt sections are malformed: {receipt_path}")
    dataset = str(scope.get("dataset", ""))
    chrom = str(scope.get("chrom", ""))
    if not dataset or chrom not in AUTOSOMES:
        raise AuditError(f"partition scope is invalid: {receipt_path}")
    if tuple(map(str, parameters.get("primary_hp_families", ()))) != (
        "1",
        "2",
    ):
        raise AuditError(f"unexpected primary HP families: {receipt_path}")
    if parameters.get("require_exact_known_phase_set") is not True:
        raise AuditError(f"phase-set contract is not exact-known: {receipt_path}")
    if int(parameters.get("max_block_size", 0)) != 8:
        raise AuditError(f"unexpected max block size: {receipt_path}")
    legacy_components_path = verify_declared_identity(
        outputs.get("legacy_components", {}), expected_parent=partition_dir
    )
    verify_declared_identity(
        outputs.get("blocks", {}), expected_parent=partition_dir
    )
    membership_path = verify_declared_identity(
        outputs.get("site_membership", {}), expected_parent=partition_dir
    )
    constraints_path = verify_declared_identity(
        outputs.get("cut_constraints", {}), expected_parent=partition_dir
    )
    for declared in inputs.values():
        if not isinstance(declared, Mapping):
            raise AuditError(f"malformed partition input identity: {receipt_path}")
        verify_declared_identity(declared)
    return SourcePartition(
        dataset=dataset,
        chrom=chrom,
        partition_dir=partition_dir.resolve(),
        receipt=receipt,
        receipt_identity={
            "path": str(receipt_path.resolve()),
            "size_bytes": receipt_path.stat().st_size,
            "sha256": receipt_sha,
        },
        legacy_components_path=legacy_components_path,
        membership_path=membership_path,
        constraints_path=constraints_path,
        source_context=dict(source_context),
    )


def verify_full_root(source_root: Path) -> tuple[list[Path], list[str], dict[str, Any]]:
    if not source_root.is_dir() or source_root.is_symlink():
        raise AuditError(f"full source root is invalid: {source_root}")
    marker_path = source_root / "_SUCCESS"
    marker = strict_json_load(marker_path)
    if not isinstance(marker, Mapping):
        raise AuditError("_SUCCESS is not an object")
    if marker.get("schema_name") != SOURCE_SUCCESS_SCHEMA:
        raise AuditError("unexpected _SUCCESS schema")
    expected_marker = {
        "all_pass": True,
        "comprehensive": True,
        "sample": "HCC1395",
        "scope": list(AUTOSOMES),
    }
    for key, value in expected_marker.items():
        if marker.get(key) != value:
            raise AuditError(f"_SUCCESS mismatch for {key}")
    run_declared = marker.get("run_receipt")
    if not isinstance(run_declared, Mapping):
        raise AuditError("_SUCCESS lacks run_receipt")
    receipt_path = Path(str(run_declared.get("path", "")))
    if receipt_path.resolve() != (source_root / "receipt.json").resolve():
        raise AuditError("_SUCCESS run receipt escapes source root")
    receipt_sha = verify_sha256_sidecar(receipt_path)
    if run_declared.get("sha256") != receipt_sha:
        raise AuditError("_SUCCESS run receipt SHA-256 mismatch")
    receipt = strict_json_load(receipt_path)
    if (
        not isinstance(receipt, Mapping)
        or receipt.get("schema_name") != SOURCE_RUN_SCHEMA
        or receipt.get("all_pass") is not True
        or receipt.get("comprehensive_all_pass") is not True
    ):
        raise AuditError("canonical run receipt is not comprehensive all-pass")
    if receipt.get("sample") != "HCC1395" or receipt.get("scope") != {
        "chromosomes": list(AUTOSOMES),
        "test_mode": False,
    }:
        raise AuditError("canonical run scope mismatch")
    expected_totals = {
        "chromosomes": 22,
        "ssnv_sites": 79687,
        "k_gt8_components": 408,
        "k_gt8_sites": 47570,
        "k_gt8_max_k": 3574,
        "partitioned_chromosomes": 21,
        "zero_target_skipped_chromosomes": 1,
    }
    if receipt.get("totals") != expected_totals:
        raise AuditError("canonical run totals mismatch")
    entries = receipt.get("chromosomes")
    if not isinstance(entries, list) or len(entries) != 22:
        raise AuditError("canonical chromosome receipts are incomplete")
    by_chrom = {str(entry.get("chrom")): entry for entry in entries}
    if tuple(by_chrom) != AUTOSOMES:
        raise AuditError("canonical chromosome receipt order/scope drift")
    partition_dirs: list[Path] = []
    skipped: list[str] = []
    for chrom in AUTOSOMES:
        entry = by_chrom[chrom]
        status = entry.get("partition_status")
        declared = entry.get("partition")
        if not isinstance(declared, Mapping):
            raise AuditError(f"missing partition stage identity for {chrom}")
        stage_path = verify_declared_identity(declared)
        stage = strict_json_load(stage_path)
        if (
            not isinstance(stage, Mapping)
            or stage.get("all_pass") is not True
            or stage.get("chrom") != chrom
            or stage.get("status") != status
        ):
            raise AuditError(f"partition stage receipt mismatch for {chrom}")
        if status == "COMPLETED":
            partition_dirs.append(source_root / "chromosomes" / chrom / "partition")
        elif status == "SKIP_NO_K_GT8_TARGET":
            skipped.append(chrom)
        else:
            raise AuditError(f"unexpected partition status for {chrom}: {status}")
    context = {
        "mode": "full",
        "source_root": str(source_root.resolve()),
        "success_marker": file_identity(marker_path),
        "run_receipt": {
            "path": str(receipt_path.resolve()),
            "size_bytes": receipt_path.stat().st_size,
            "sha256": receipt_sha,
        },
    }
    return partition_dirs, skipped, context


def validate_source_counts(
    source: SourcePartition,
    components: Mapping[str, ComponentInfo],
    legacy_metrics: Mapping[str, Mapping[str, int]],
    constraints: Sequence[Mapping[str, Any]],
    active_by_component: Mapping[str, set[int]],
) -> dict[str, bool]:
    counts = source.receipt["counts"]
    total_weight = sum(int(row["weight"]) for row in constraints)
    retained_weight = sum(
        int(row["weight"]) for row in constraints if row["category"] == RETAINED
    )
    lost_weight = sum(
        int(row["weight"]) for row in constraints if row["category"] != RETAINED
    )
    unavoidable_patterns = sum(
        row["category"] == "unavoidable" for row in constraints
    )
    constraints_by_component: dict[str, Counter[str]] = defaultdict(Counter)
    for row in constraints:
        component_counts = constraints_by_component[
            str(row["legacy_component_id"])
        ]
        component_counts["patterns"] += 1
        component_counts["total_weight"] += int(row["weight"])
        if row["category"] == RETAINED:
            component_counts["retained_patterns"] += 1
            component_counts["retained_weight"] += int(row["weight"])
        else:
            component_counts["lost_patterns"] += 1
            component_counts["lost_weight"] += int(row["weight"])
        if row["category"] == "unavoidable":
            component_counts["unavoidable_patterns"] += 1

    per_component_legacy_match = set(legacy_metrics) == set(components)
    if per_component_legacy_match:
        for component_id, component in components.items():
            legacy = legacy_metrics[component_id]
            observed = constraints_by_component[component_id]
            expected_values = {
                "pre_cap_k": component.k,
                "primary_active_site_count": len(
                    active_by_component.get(component_id, set())
                ),
                "exact_pattern_count": observed["patterns"],
                "raw_total_molecule_weight": observed["total_weight"],
                "raw_retained_molecule_weight": observed["retained_weight"],
                "raw_lost_molecule_weight": observed["lost_weight"],
                "retained_exact_pattern_count": observed[
                    "retained_patterns"
                ],
                "lost_exact_pattern_count": observed["lost_patterns"],
                "unavoidable_pattern_count": observed[
                    "unavoidable_patterns"
                ],
            }
            if any(legacy[name] != value for name, value in expected_values.items()):
                per_component_legacy_match = False
                break
    checks = {
        "component_count_matches_receipt": len(components)
        == int(counts["target_components"]),
        "site_count_matches_receipt": sum(component.k for component in components.values())
        == int(counts["target_sites"]),
        "every_component_k_gt_max_block_size": all(
            component.k > int(source.receipt["parameters"]["max_block_size"])
            for component in components.values()
        ),
        "constraint_rows_match_receipt": len(constraints)
        == int(counts["exact_patterns"]),
        "total_weight_matches_receipt": total_weight
        == int(counts["raw_total_molecule_weight"]),
        "retained_weight_matches_receipt": retained_weight
        == int(counts["raw_retained_molecule_weight"]),
        "lost_weight_matches_receipt": lost_weight
        == int(counts["raw_lost_molecule_weight"]),
        "unavoidable_patterns_match_receipt": unavoidable_patterns
        == int(counts["unavoidable_patterns"]),
        "component_active_site_sum_matches_receipt": sum(
            len(values) for values in active_by_component.values()
        )
        == int(counts["primary_active_sites_component_sum"]),
        "constraint_active_sites_match_membership_flags": all(
            active_by_component.get(component_id, set())
            == set(component.primary_active_positions)
            for component_id, component in components.items()
        ),
        "per_component_legacy_metrics_match_constraints": (
            per_component_legacy_match
        ),
    }
    if not all(checks.values()):
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise AuditError(f"source conservation failed for {source.chrom}: {failed}")
    return checks


def run(args: argparse.Namespace) -> dict[str, Any]:
    if args.output_dir.exists() or args.output_dir.is_symlink():
        raise AuditError(f"output directory must not exist: {args.output_dir}")
    if args.mode == "full":
        if args.source_root is None or args.partition_dir:
            raise AuditError("full mode requires only --source-root")
        partition_dirs, skipped_chromosomes, source_context = verify_full_root(
            args.source_root
        )
    else:
        if args.source_root is not None or not args.partition_dir:
            raise AuditError("probe mode requires one or more --partition-dir")
        partition_dirs = [path.resolve() for path in args.partition_dir]
        skipped_chromosomes = []
        source_context = {
            "mode": "probe",
            "partition_directories": [str(path) for path in partition_dirs],
        }

    sources: list[SourcePartition] = []
    seen_scopes: set[tuple[str, str]] = set()
    for partition_dir in partition_dirs:
        source = verify_partition_source(partition_dir, source_context)
        scope_key = (source.dataset, source.chrom)
        if scope_key in seen_scopes:
            raise AuditError(f"duplicate partition scope: {scope_key}")
        seen_scopes.add(scope_key)
        sources.append(source)
    sources.sort(key=lambda item: AUTOSOMES.index(item.chrom))
    datasets = {source.dataset for source in sources}
    if len(datasets) != 1:
        raise AuditError(f"multiple datasets in one audit: {sorted(datasets)}")

    unit_rows: list[dict[str, Any]] = []
    source_checks: dict[str, dict[str, bool]] = {}
    source_components: dict[str, set[str]] = {}
    input_records: list[dict[str, Any]] = []
    for source in sources:
        components = load_membership(source)
        legacy_metrics = load_legacy_component_metrics(source)
        constraints, active_by_component = load_constraints(source, components)
        checks = validate_source_counts(
            source,
            components,
            legacy_metrics,
            constraints,
            active_by_component,
        )
        source_checks[source.chrom] = checks
        source_components[source.chrom] = set(components)
        unit_rows.extend(aggregate_units(components, constraints))
        input_records.append(
            {
                "dataset": source.dataset,
                "chrom": source.chrom,
                "partition_receipt": dict(source.receipt_identity),
                "legacy_components": file_identity(
                    source.legacy_components_path
                ),
                "site_membership": file_identity(source.membership_path),
                "cut_constraints": file_identity(source.constraints_path),
            }
        )

    unit_rows.sort(
        key=lambda row: (
            AUTOSOMES.index(str(row["chrom"])),
            str(row["legacy_component_id"]),
            int(row["hp_family"]),
            str(row["phase_set"]),
        )
    )
    pair_rows = aggregate_component_pairs(unit_rows)
    pair_rows.sort(
        key=lambda row: (
            AUTOSOMES.index(str(row["chrom"])),
            str(row["legacy_component_id"]),
        )
    )
    summary = summarize_rows(
        unit_rows,
        pair_rows,
        source_components=source_components,
        skipped_chromosomes=skipped_chromosomes,
    )
    summary_rows = build_summary_tsv(summary)
    semantic_hash = semantic_sha256(
        {"units": unit_rows, "pairs": pair_rows, "summary": summary}
    )

    args.output_dir.mkdir(parents=True, exist_ok=False)
    dataset = next(iter(datasets))
    unit_path = args.output_dir / f"{dataset}.hp_ps_observed_units.tsv.gz"
    pair_path = args.output_dir / f"{dataset}.hp1_hp2_paired_components.tsv.gz"
    summary_json_path = args.output_dir / "summary.json"
    summary_tsv_path = args.output_dir / "summary.tsv"
    write_tsv_gz(unit_path, UNIT_FIELDS, unit_rows)
    write_tsv_gz(pair_path, PAIR_FIELDS, pair_rows)
    write_json(summary_json_path, summary)
    write_tsv(summary_tsv_path, SUMMARY_FIELDS, summary_rows)

    output_identities = {
        "unit_table": file_identity(unit_path),
        "paired_component_table": file_identity(pair_path),
        "summary_json": file_identity(summary_json_path),
        "summary_tsv": file_identity(summary_tsv_path),
    }
    checks = {
        "source_checks_all_pass": all(
            all(chrom_checks.values())
            for chrom_checks in source_checks.values()
        ),
        "unit_rows_have_observed_denominator": all(
            row["ratio_status"] == "OBSERVED_CONSTRAINT_DENOMINATOR"
            and row["retention_ratio"] != ""
            for row in unit_rows
        ),
        "unit_pattern_mass_conserved": all(
            int(row["total_pattern_rows"])
            == int(row["retained_pattern_rows"])
            + int(row["cut_lost_pattern_rows"])
            + int(row["unavoidable_pattern_rows"])
            for row in unit_rows
        ),
        "unit_molecule_incidence_mass_conserved": all(
            int(row["total_molecule_component_incidence_weight"])
            == int(row["retained_molecule_component_incidence_weight"])
            + int(row["cut_lost_molecule_component_incidence_weight"])
            + int(row["unavoidable_molecule_component_incidence_weight"])
            for row in unit_rows
        ),
        "no_unobserved_unit_rows_synthesized": len(unit_rows)
        == summary["counts"]["observed_constraint_units"],
        "full_scope_is_exact_autosomes": (
            args.mode != "full"
            or (
                set(source.chrom for source in sources)
                | set(skipped_chromosomes)
            )
            == set(AUTOSOMES)
            and len(sources) + len(skipped_chromosomes) == len(AUTOSOMES)
        ),
    }
    if not all(checks.values()):
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise AuditError(f"final audit checks failed: {failed}")

    receipt_path = args.output_dir / "receipt.json"
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "scope": {
            "mode": args.mode,
            "dataset": dataset,
            "chromosomes_with_partition_targets": [
                source.chrom for source in sources
            ],
            "skipped_zero_target_chromosomes": sorted(skipped_chromosomes),
            "unit_grain": (
                "dataset x chromosome x legacy_component_id x hp_family x "
                "exact known phase_set"
            ),
            "scope_ceiling": "observed_constraint_units_only",
        },
        "definitions": {
            "aggregation_weight": "molecule_x_component_incidence",
            "aggregation_warning": (
                "sums are molecule-by-component incidences, not unique reads; "
                "one molecule may contribute once in each component it spans"
            ),
            "total": "retained + cut_lost + unavoidable",
            "cut_lost": "disposition == cut",
            "unavoidable": "disposition starts with unavoidable_",
            "nonretained": "cut_lost + unavoidable",
            "retention_ratio": "retained / total for observed units only",
            "support_strata": ["1-4", "5-19", "20-49", ">=50"],
            "headline_eligible": "total_weight >= 20 and total_pattern_rows >= 5",
            "unobserved_opportunities": (
                "not reconstructed from sparse calls because an independent "
                "zero-opportunity HP/PS denominator is not semantically "
                "identified by the partition outputs; no synthetic ratio is emitted"
            ),
            "vaf_role": "not used by this engineering retention audit",
        },
        "parameters": {
            "eligible_min_total_weight": 20,
            "eligible_min_pattern_rows": 5,
            "quantile_method": "type7_linear_interpolation_unweighted_units",
            "upstream_input_hash_verification": True,
        },
        "summary": summary,
        "source_checks": source_checks,
        "checks": checks,
        "tool": file_identity(Path(__file__).resolve()),
        "source_context": source_context,
        "inputs": input_records,
        "outputs": output_identities,
        "semantic_result_sha256": semantic_hash,
        "command": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "claim_ceiling": (
            "descriptive retention of observed HP/PS exact-pattern constraints; "
            "not unique-molecule coverage, not unobserved opportunity recovery, "
            "and not proof of a unique true evolutionary tree"
        ),
    }
    write_json(receipt_path, receipt)
    write_sha256_sidecar(receipt_path)
    print(
        json.dumps(
            {
                "all_pass": True,
                "mode": args.mode,
                "dataset": dataset,
                "chromosomes": [source.chrom for source in sources],
                "observed_constraint_units": len(unit_rows),
                "eligible_headline_units": summary["counts"][
                    "eligible_headline_units"
                ],
                "weighted_retention_ratio": summary[
                    "molecule_component_incidence_totals"
                ]["weighted_retention_ratio"],
                "hp1_hp2_paired_components": len(pair_rows),
                "receipt": str(receipt_path.resolve()),
                "receipt_sha256": sha256_path(receipt_path),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read-only HP/PS observed-constraint retention audit for k>8 "
            "partitions"
        )
    )
    parser.add_argument("--mode", choices=("full", "probe"), default="full")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument(
        "--partition-dir",
        type=Path,
        action="append",
        default=[],
        help="Direct partition directory; repeat in probe mode",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    try:
        run(parse_args())
    except (AuditError, OSError, ValueError, KeyError, csv.Error) as exc:
        print(f"FAIL-CLOSED: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
