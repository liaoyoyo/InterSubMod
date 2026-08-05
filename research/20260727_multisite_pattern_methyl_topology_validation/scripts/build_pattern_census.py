#!/usr/bin/env python3
"""Build an exact-raw-HP read-pattern census over exact-PS topology regions.

Topology rows are located by exact ``phase_set × hp_family``.  Every sparse
read call is then projected onto the ordered ``active_positions`` of every
matching topology row that contains at least one observed position.  The
analysis grain remains exact ``hp_raw``; family labels are never used to merge
raw HP strata.

Candidate-read output is sharded by dataset and chromosome.  Each shard is
externally sorted by ``dataset, chrom, region_id, hp_raw, qname_sha256`` so the
all-seven-dataset run does not require holding all candidate rows in memory.
The manifest defines the deterministic global shard order.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import shutil
import subprocess
import threading
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TextIO


AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))
DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
DEFAULT_TOPOLOGY_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1"
)
DEFAULT_HCC1395_SPARSE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2"
)
DEFAULT_PRODUCTION_SPARSE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260723_production_exact_ps_strict_read_linkage/all7_production_v1"
)

EXPECTED_HP_FAMILY = {
    ".": "none",
    "1": "1",
    "1-1": "1",
    "1-2": "1",
    "2": "2",
    "2-1": "2",
    "2-2": "2",
    "3": "3",
    "4": "4",
}
SPARSE_REQUIRED_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "qname_sha256",
    "read_group",
    "start0",
    "end0",
    "mapq",
    "strand",
    "hp_raw",
    "hp_family",
    "phase_set",
    "positions1",
    "call_codes",
)
CANDIDATE_FIELDS = (
    "dataset",
    "chrom",
    "region_id",
    "unit_id",
    "phase_set",
    "hp_family",
    "hp_raw",
    "qname_sha256",
    "molecule_id",
    "read_group",
    "strand",
    "mapq",
    "start0",
    "end0",
    "active_positions",
    "pattern",
    "complete_pattern",
    "n_active_bits",
)
PATTERN_COUNT_FIELDS = (
    "dataset",
    "chrom",
    "region_id",
    "unit_id",
    "phase_set",
    "hp_family",
    "hp_raw",
    "active_positions",
    "n_active_bits",
    "k_ge_2",
    "pair_full4",
    "k_ge_3",
    "n_total",
    "n_complete",
    "n_partial",
    "complete_state_count",
    "partial_state_count",
    "state_count_json",
    "partial_state_count_json",
    "formal_n5",
    "formal_n8",
    "formal_n10",
)
MARKER_FIELDS = (
    "dataset",
    "chrom",
    "region_id",
    "unit_id",
    "phase_set",
    "hp_family",
    "position1",
    "marker_index",
    "n_active_bits",
)
MANIFEST_FIELDS = (
    "sort_ordinal",
    "dataset",
    "chrom",
    "relative_path",
    "source_sparse_path",
    "n_source_rows",
    "n_projected_rows",
    "sha256",
    "size_bytes",
)


class CensusContractError(RuntimeError):
    """Raised when an input or output violates the census contract."""


@dataclass(frozen=True)
class TopologyRegion:
    dataset: str
    chrom: str
    region_id: str
    unit_id: str
    phase_set: str
    hp_family: str
    active_positions: tuple[int, ...]

    @property
    def n_active_bits(self) -> int:
        return len(self.active_positions)

    @property
    def active_positions_text(self) -> str:
        return ",".join(str(position) for position in self.active_positions)


@dataclass
class PatternAccumulator:
    region: TopologyRegion
    hp_raw: str
    n_total: int = 0
    complete: Counter[str] = field(default_factory=Counter)
    partial: Counter[str] = field(default_factory=Counter)

    def add(self, pattern: str, complete_pattern: bool) -> None:
        self.n_total += 1
        if complete_pattern:
            self.complete[pattern] += 1
        else:
            self.partial[pattern] += 1

    def formal_gate(self, minimum_state_n: int) -> bool:
        return (
            self.region.n_active_bits >= 2
            and sum(self.complete.values()) >= 40
            and sum(count >= minimum_state_n for count in self.complete.values()) >= 2
        )

    def as_row(self) -> dict[str, object]:
        pair_full4 = (
            self.region.n_active_bits == 2
            and sum(self.complete.values()) >= 40
            and all(self.complete[state] >= 5 for state in ("RR", "RA", "AR", "AA"))
        )
        return {
            "dataset": self.region.dataset,
            "chrom": self.region.chrom,
            "region_id": self.region.region_id,
            "unit_id": self.region.unit_id,
            "phase_set": self.region.phase_set,
            "hp_family": self.region.hp_family,
            "hp_raw": self.hp_raw,
            "active_positions": self.region.active_positions_text,
            "n_active_bits": self.region.n_active_bits,
            "k_ge_2": bool_text(self.region.n_active_bits >= 2),
            "pair_full4": bool_text(pair_full4),
            "k_ge_3": bool_text(self.region.n_active_bits >= 3),
            "n_total": self.n_total,
            "n_complete": sum(self.complete.values()),
            "n_partial": sum(self.partial.values()),
            "complete_state_count": len(self.complete),
            "partial_state_count": len(self.partial),
            "state_count_json": counter_json(self.complete),
            "partial_state_count_json": counter_json(self.partial),
            "formal_n5": bool_text(self.formal_gate(5)),
            "formal_n8": bool_text(self.formal_gate(8)),
            "formal_n10": bool_text(self.formal_gate(10)),
        }


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def counter_json(counter: Mapping[str, int]) -> str:
    payload = {key: int(counter[key]) for key in sorted(counter)}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def chromosome_key(chrom: str) -> int:
    if chrom not in AUTOSOMES:
        raise CensusContractError(f"unsupported chromosome: {chrom!r}")
    return int(chrom.removeprefix("chr"))


def parse_csv_values(values: Sequence[str] | None) -> tuple[str, ...]:
    parsed: list[str] = []
    for raw in values or ():
        parsed.extend(value.strip() for value in raw.split(",") if value.strip())
    return tuple(parsed)


def resolve_samples(*, all7: bool, sample_values: Sequence[str] | None) -> tuple[str, ...]:
    requested = parse_csv_values(sample_values)
    if all7:
        if requested:
            raise CensusContractError("--all7 cannot be combined with --sample")
        return DATASETS
    if not requested:
        raise CensusContractError("choose --all7 or provide at least one --sample")
    unknown = sorted(set(requested) - set(DATASETS))
    if unknown:
        raise CensusContractError(f"unsupported datasets: {unknown}")
    if len(requested) != len(set(requested)):
        raise CensusContractError("duplicate --sample selection")
    return tuple(sorted(requested))


def resolve_chromosomes(chrom_values: Sequence[str] | None) -> tuple[str, ...]:
    requested = parse_csv_values(chrom_values)
    if not requested:
        return AUTOSOMES
    unknown = sorted(set(requested) - set(AUTOSOMES))
    if unknown:
        raise CensusContractError(f"unsupported chromosomes: {unknown}")
    if len(requested) != len(set(requested)):
        raise CensusContractError("duplicate --chrom selection")
    return tuple(sorted(requested, key=chromosome_key))


def require_fields(
    path: Path, fieldnames: Sequence[str] | None, required: Sequence[str]
) -> None:
    present = set(fieldnames or ())
    missing = sorted(set(required) - present)
    if missing:
        raise CensusContractError(f"{path}: missing required columns {missing}")


def require_scalar_text(
    value: object, *, path: Path, line_number: int, field_name: str
) -> str:
    if not isinstance(value, str) or value == "":
        raise CensusContractError(
            f"{path}:{line_number}: {field_name} must be a non-empty string"
        )
    if "\t" in value or "\r" in value or "\n" in value:
        raise CensusContractError(
            f"{path}:{line_number}: {field_name} contains a TSV line delimiter"
        )
    return value


def parse_integer(
    value: str,
    *,
    path: Path,
    line_number: int,
    field_name: str,
    minimum: int,
    maximum: int | None = None,
) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise CensusContractError(
            f"{path}:{line_number}: {field_name} is not an integer: {value!r}"
        ) from exc
    if parsed < minimum or (maximum is not None and parsed > maximum):
        bounds = f">={minimum}" if maximum is None else f"in [{minimum}, {maximum}]"
        raise CensusContractError(
            f"{path}:{line_number}: {field_name} must be {bounds}: {parsed}"
        )
    return parsed


def parse_positions(
    value: str, *, path: Path, line_number: int
) -> tuple[int, ...]:
    if value == "":
        raise CensusContractError(f"{path}:{line_number}: positions1 must not be empty")
    positions: list[int] = []
    for token in value.split(","):
        positions.append(
            parse_integer(
                token,
                path=path,
                line_number=line_number,
                field_name="positions1",
                minimum=1,
            )
        )
    if positions != sorted(positions) or len(positions) != len(set(positions)):
        raise CensusContractError(
            f"{path}:{line_number}: positions1 must be strictly increasing and unique"
        )
    return tuple(positions)


def require_sha256(
    value: str, *, path: Path, line_number: int, field_name: str
) -> str:
    require_scalar_text(
        value, path=path, line_number=line_number, field_name=field_name
    )
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise CensusContractError(
            f"{path}:{line_number}: {field_name} must be a lowercase SHA-256 digest"
        )
    return value


def topology_path(topology_root: Path, sample: str) -> Path:
    path = topology_root / "samples" / sample / f"{sample}.topology.jsonl"
    if not path.is_file():
        raise CensusContractError(f"missing topology JSONL: {path}")
    return path


def sparse_path(
    *,
    sample: str,
    chrom: str,
    hcc1395_sparse_root: Path,
    production_sparse_root: Path,
) -> Path:
    if sample == "HCC1395":
        path = (
            hcc1395_sparse_root
            / "chromosomes"
            / chrom
            / "extraction"
            / f"{sample}.{chrom}.molecule_sparse_calls.tsv.gz"
        )
    else:
        path = (
            production_sparse_root
            / "samples"
            / sample
            / "chromosomes"
            / chrom
            / "extraction"
            / f"{sample}.{chrom}.molecule_sparse_calls.tsv.gz"
        )
    if not path.is_file():
        raise CensusContractError(f"missing sparse read-call TSV: {path}")
    return path


def load_topology_regions(
    path: Path, *, sample: str, selected_chromosomes: set[str]
) -> dict[str, tuple[TopologyRegion, ...]]:
    by_chrom: dict[str, list[TopologyRegion]] = defaultdict(list)
    seen_region_ids: set[str] = set()
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                raise CensusContractError(f"{path}:{line_number}: blank JSONL row")
            try:
                payload = json.loads(raw_line)
            except json.JSONDecodeError as exc:
                raise CensusContractError(
                    f"{path}:{line_number}: malformed topology JSON"
                ) from exc
            if not isinstance(payload, dict):
                raise CensusContractError(
                    f"{path}:{line_number}: topology row must be a JSON object"
                )
            if payload.get("schema_name") != "intersubmod.exact_ps_cpp_topology_af.unit":
                raise CensusContractError(
                    f"{path}:{line_number}: unexpected topology schema_name"
                )
            if payload.get("schema_version") != "1.0.0":
                raise CensusContractError(
                    f"{path}:{line_number}: unexpected topology schema_version"
                )
            if payload.get("sample") != sample:
                raise CensusContractError(
                    f"{path}:{line_number}: sample mismatch "
                    f"{payload.get('sample')!r} != {sample!r}"
                )
            chrom = require_scalar_text(
                payload.get("chrom"),
                path=path,
                line_number=line_number,
                field_name="chrom",
            )
            chromosome_key(chrom)
            region_id = require_scalar_text(
                payload.get("region_id"),
                path=path,
                line_number=line_number,
                field_name="region_id",
            )
            if region_id in seen_region_ids:
                raise CensusContractError(
                    f"{path}:{line_number}: duplicate region_id {region_id!r}"
                )
            seen_region_ids.add(region_id)
            unit_id = require_scalar_text(
                payload.get("unit_id"),
                path=path,
                line_number=line_number,
                field_name="unit_id",
            )
            phase_set = require_scalar_text(
                payload.get("phase_set"),
                path=path,
                line_number=line_number,
                field_name="phase_set",
            )
            hp_family = require_scalar_text(
                payload.get("hp_family"),
                path=path,
                line_number=line_number,
                field_name="hp_family",
            )
            if hp_family not in {"1", "2"}:
                raise CensusContractError(
                    f"{path}:{line_number}: topology hp_family must be 1 or 2"
                )
            active_raw = payload.get("active_positions")
            if not isinstance(active_raw, list) or any(
                isinstance(position, bool) or not isinstance(position, int)
                for position in active_raw
            ):
                raise CensusContractError(
                    f"{path}:{line_number}: active_positions must be an integer list"
                )
            active_positions = tuple(active_raw)
            if (
                any(position <= 0 for position in active_positions)
                or active_positions != tuple(sorted(active_positions))
                or len(active_positions) != len(set(active_positions))
            ):
                raise CensusContractError(
                    f"{path}:{line_number}: active_positions must be positive, "
                    "strictly increasing, and unique"
                )
            active_bit_count = payload.get("active_bit_count")
            if (
                isinstance(active_bit_count, bool)
                or not isinstance(active_bit_count, int)
                or active_bit_count != len(active_positions)
            ):
                raise CensusContractError(
                    f"{path}:{line_number}: active_bit_count does not match "
                    "active_positions"
                )
            if chrom in selected_chromosomes:
                by_chrom[chrom].append(
                    TopologyRegion(
                        dataset=sample,
                        chrom=chrom,
                        region_id=region_id,
                        unit_id=unit_id,
                        phase_set=phase_set,
                        hp_family=hp_family,
                        active_positions=active_positions,
                    )
                )
    return {
        chrom: tuple(sorted(regions, key=lambda region: region.region_id))
        for chrom, regions in by_chrom.items()
    }


def build_position_index(
    regions: Sequence[TopologyRegion],
) -> dict[tuple[str, str, int], tuple[TopologyRegion, ...]]:
    mutable: dict[tuple[str, str, int], list[TopologyRegion]] = defaultdict(list)
    for region in regions:
        for position in region.active_positions:
            mutable[(region.phase_set, region.hp_family, position)].append(region)
    return {
        key: tuple(sorted(value, key=lambda region: region.region_id))
        for key, value in mutable.items()
    }


@contextmanager
def deterministic_gzip_writer(path: Path) -> Iterator[TextIO]:
    with path.open("xb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            with io.TextIOWrapper(
                gzip_handle, encoding="utf-8", newline="", write_through=True
            ) as text_handle:
                yield text_handle


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_tsv_row(handle: TextIO, fields: Sequence[str], row: Mapping[str, object]) -> None:
    values: list[str] = []
    for field_name in fields:
        value = str(row[field_name])
        if "\t" in value or "\r" in value or "\n" in value:
            raise CensusContractError(
                f"output field {field_name!r} contains a TSV line delimiter"
            )
        values.append(value)
    handle.write("\t".join(values) + "\n")


def candidate_sort_command(sort_binary: str, temporary_directory: Path) -> list[str]:
    return [
        sort_binary,
        "--stable",
        "--field-separator=\t",
        "--key=3,3",
        "--key=7,7",
        "--key=8,8",
        "--key=9,9",
        "--key=4,4",
        "--key=16,16",
        f"--temporary-directory={temporary_directory}",
    ]


def normalize_call(code: str) -> str:
    return code if code in {"R", "A"} else "X"


def candidate_row(
    *,
    region: TopologyRegion,
    sparse_row: Mapping[str, str],
    pattern: str,
    complete_pattern: bool,
) -> dict[str, object]:
    return {
        "dataset": region.dataset,
        "chrom": region.chrom,
        "region_id": region.region_id,
        "unit_id": region.unit_id,
        "phase_set": region.phase_set,
        "hp_family": region.hp_family,
        "hp_raw": sparse_row["hp_raw"],
        "qname_sha256": sparse_row["qname_sha256"],
        "molecule_id": sparse_row["molecule_id"],
        "read_group": sparse_row["read_group"],
        "strand": sparse_row["strand"],
        "mapq": sparse_row["mapq"],
        "start0": sparse_row["start0"],
        "end0": sparse_row["end0"],
        "active_positions": region.active_positions_text,
        "pattern": pattern,
        "complete_pattern": bool_text(complete_pattern),
        "n_active_bits": region.n_active_bits,
    }


def stream_sparse_candidates(
    *,
    path: Path,
    sample: str,
    chrom: str,
    regions: Sequence[TopologyRegion],
    output_handle: TextIO,
) -> tuple[int, int, dict[tuple[str, str], PatternAccumulator]]:
    position_index = build_position_index(regions)
    accumulators: dict[tuple[str, str], PatternAccumulator] = {}
    source_rows = 0
    projected_rows = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, SPARSE_REQUIRED_FIELDS)
        for line_number, row in enumerate(reader, start=2):
            source_rows += 1
            if row["dataset"] != sample or row["chrom"] != chrom:
                raise CensusContractError(
                    f"{path}:{line_number}: dataset/chrom mismatch"
                )
            molecule_id = require_sha256(
                row["molecule_id"],
                path=path,
                line_number=line_number,
                field_name="molecule_id",
            )
            qname_sha256 = require_sha256(
                row["qname_sha256"],
                path=path,
                line_number=line_number,
                field_name="qname_sha256",
            )
            hp_raw = require_scalar_text(
                row["hp_raw"],
                path=path,
                line_number=line_number,
                field_name="hp_raw",
            )
            hp_family = require_scalar_text(
                row["hp_family"],
                path=path,
                line_number=line_number,
                field_name="hp_family",
            )
            expected_family = EXPECTED_HP_FAMILY.get(hp_raw)
            if expected_family is None or expected_family != hp_family:
                raise CensusContractError(
                    f"{path}:{line_number}: inconsistent hp_raw/hp_family "
                    f"{hp_raw!r}/{hp_family!r}"
                )
            phase_set = row["phase_set"]
            if "\t" in phase_set or "\r" in phase_set or "\n" in phase_set:
                raise CensusContractError(
                    f"{path}:{line_number}: phase_set contains a TSV line delimiter"
                )
            positions = parse_positions(
                row["positions1"], path=path, line_number=line_number
            )
            codes = row["call_codes"]
            if len(codes) != len(positions):
                raise CensusContractError(
                    f"{path}:{line_number}: positions1/call_codes length mismatch"
                )
            start0 = parse_integer(
                row["start0"],
                path=path,
                line_number=line_number,
                field_name="start0",
                minimum=0,
            )
            end0 = parse_integer(
                row["end0"],
                path=path,
                line_number=line_number,
                field_name="end0",
                minimum=start0 + 1,
            )
            mapq = parse_integer(
                row["mapq"],
                path=path,
                line_number=line_number,
                field_name="mapq",
                minimum=0,
                maximum=255,
            )
            strand = row["strand"]
            if strand not in {"+", "-"}:
                raise CensusContractError(
                    f"{path}:{line_number}: strand must be '+' or '-'"
                )
            read_group = row["read_group"]
            if "\t" in read_group or "\r" in read_group or "\n" in read_group:
                raise CensusContractError(
                    f"{path}:{line_number}: read_group contains a TSV line delimiter"
                )

            normalized_by_position = {
                position: normalize_call(code)
                for position, code in zip(positions, codes)
            }
            matched_regions: dict[str, TopologyRegion] = {}
            for position in positions:
                for region in position_index.get(
                    (phase_set, hp_family, position), ()
                ):
                    matched_regions[region.region_id] = region
            normalized_row = dict(row)
            normalized_row.update(
                {
                    "molecule_id": molecule_id,
                    "qname_sha256": qname_sha256,
                    "hp_raw": hp_raw,
                    "hp_family": hp_family,
                    "phase_set": phase_set,
                    "start0": str(start0),
                    "end0": str(end0),
                    "mapq": str(mapq),
                    "strand": strand,
                    "read_group": read_group,
                }
            )
            for region in sorted(
                matched_regions.values(), key=lambda value: value.region_id
            ):
                pattern = "".join(
                    normalized_by_position.get(position, "X")
                    for position in region.active_positions
                )
                complete_pattern = all(code in {"R", "A"} for code in pattern)
                output_row = candidate_row(
                    region=region,
                    sparse_row=normalized_row,
                    pattern=pattern,
                    complete_pattern=complete_pattern,
                )
                write_tsv_row(output_handle, CANDIDATE_FIELDS, output_row)
                accumulator_key = (region.region_id, hp_raw)
                accumulator = accumulators.setdefault(
                    accumulator_key,
                    PatternAccumulator(region=region, hp_raw=hp_raw),
                )
                accumulator.add(pattern, complete_pattern)
                projected_rows += 1
    return source_rows, projected_rows, accumulators


def write_sorted_candidate_shard(
    *,
    sparse_calls: Path,
    sample: str,
    chrom: str,
    regions: Sequence[TopologyRegion],
    output_path: Path,
    sort_binary: str,
    temporary_directory: Path,
) -> tuple[int, int, dict[tuple[str, str], PatternAccumulator]]:
    partial_path = output_path.with_suffix(output_path.suffix + ".partial")
    if output_path.exists() or partial_path.exists():
        raise CensusContractError(f"refusing to overwrite candidate shard: {output_path}")
    command = candidate_sort_command(sort_binary, temporary_directory)
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    process: subprocess.Popen[str] | None = None
    drain_errors: list[BaseException] = []
    result: tuple[int, int, dict[tuple[str, str], PatternAccumulator]] | None = None
    with deterministic_gzip_writer(partial_path) as output_handle:
        output_handle.write("\t".join(CANDIDATE_FIELDS) + "\n")
        process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            env=environment,
        )
        if process.stdin is None or process.stdout is None or process.stderr is None:
            raise CensusContractError("failed to create external-sort pipes")

        def drain_sorted_output() -> None:
            try:
                shutil.copyfileobj(process.stdout, output_handle)
            except BaseException as exc:  # pragma: no cover - exceptional I/O path
                drain_errors.append(exc)

        drain_thread = threading.Thread(target=drain_sorted_output, daemon=True)
        drain_thread.start()
        producer_error: BaseException | None = None
        try:
            result = stream_sparse_candidates(
                path=sparse_calls,
                sample=sample,
                chrom=chrom,
                regions=regions,
                output_handle=process.stdin,
            )
        except BaseException as exc:
            producer_error = exc
        finally:
            process.stdin.close()

        return_code = process.wait()
        drain_thread.join()
        stderr = process.stderr.read().strip()
        process.stdout.close()
        process.stderr.close()
        if producer_error is not None:
            raise producer_error
        if return_code != 0:
            raise CensusContractError(
                f"candidate external sort failed ({return_code}): {stderr}"
            )
        if drain_errors:
            raise CensusContractError(
                f"failed while writing sorted candidate shard: {drain_errors[0]}"
            )
    if result is None:
        raise CensusContractError("candidate producer did not return a result")
    partial_path.rename(output_path)
    return result


def write_marker_rows(
    handle: TextIO, regions: Sequence[TopologyRegion]
) -> int:
    row_count = 0
    for region in sorted(regions, key=lambda value: value.region_id):
        for marker_index, position in enumerate(region.active_positions):
            write_tsv_row(
                handle,
                MARKER_FIELDS,
                {
                    "dataset": region.dataset,
                    "chrom": region.chrom,
                    "region_id": region.region_id,
                    "unit_id": region.unit_id,
                    "phase_set": region.phase_set,
                    "hp_family": region.hp_family,
                    "position1": position,
                    "marker_index": marker_index,
                    "n_active_bits": region.n_active_bits,
                },
            )
            row_count += 1
    return row_count


def write_pattern_count_rows(
    handle: TextIO,
    accumulators: Mapping[tuple[str, str], PatternAccumulator],
) -> tuple[int, Counter[str]]:
    row_count = 0
    gates: Counter[str] = Counter()
    for key in sorted(accumulators):
        row = accumulators[key].as_row()
        write_tsv_row(handle, PATTERN_COUNT_FIELDS, row)
        row_count += 1
        for flag in ("formal_n5", "formal_n8", "formal_n10"):
            gates[flag] += row[flag] == "true"
        gates["pair_full4"] += row["pair_full4"] == "true"
        gates["k_ge_3"] += row["k_ge_3"] == "true"
    return row_count, gates


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(
            payload,
            handle,
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
        )
        handle.write("\n")


def build_pattern_census(
    *,
    topology_root: Path,
    hcc1395_sparse_root: Path,
    production_sparse_root: Path,
    output_dir: Path,
    samples: Sequence[str],
    chromosomes: Sequence[str],
    sort_binary: str,
) -> dict[str, object]:
    if output_dir.exists():
        raise CensusContractError(f"output directory already exists: {output_dir}")
    staging_dir = output_dir.with_name(output_dir.name + ".partial")
    if staging_dir.exists():
        raise CensusContractError(f"staging directory already exists: {staging_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir.mkdir()
    shard_dir = staging_dir / "candidate_read_join"
    shard_dir.mkdir()

    selected_samples = tuple(sorted(samples))
    if not selected_samples or set(selected_samples) - set(DATASETS):
        raise CensusContractError("samples must be a non-empty subset of DATASETS")
    selected_chromosomes = tuple(sorted(chromosomes, key=chromosome_key))
    if not selected_chromosomes:
        raise CensusContractError("chromosomes must not be empty")
    sort_path = shutil.which(sort_binary)
    if sort_path is None:
        raise CensusContractError(f"external sort binary not found: {sort_binary}")

    topology_by_sample: dict[str, dict[str, tuple[TopologyRegion, ...]]] = {}
    topology_sources: dict[str, str] = {}
    selected_chromosome_set = set(selected_chromosomes)
    for sample in selected_samples:
        source_path = topology_path(topology_root, sample)
        topology_sources[sample] = str(source_path.resolve())
        topology_by_sample[sample] = load_topology_regions(
            source_path,
            sample=sample,
            selected_chromosomes=selected_chromosome_set,
        )

    marker_path = staging_dir / "marker_universe.tsv"
    counts_path = staging_dir / "pattern_counts.tsv"
    manifest_path = staging_dir / "candidate_read_join.manifest.tsv"
    marker_rows = 0
    count_rows = 0
    source_rows_total = 0
    projected_rows_total = 0
    gate_totals: Counter[str] = Counter()
    manifest_rows: list[dict[str, object]] = []

    with marker_path.open("x", encoding="utf-8", newline="") as marker_handle:
        marker_handle.write("\t".join(MARKER_FIELDS) + "\n")
        with counts_path.open("x", encoding="utf-8", newline="") as counts_handle:
            counts_handle.write("\t".join(PATTERN_COUNT_FIELDS) + "\n")
            sort_ordinal = 0
            for sample in selected_samples:
                for chrom in selected_chromosomes:
                    regions = topology_by_sample[sample].get(chrom, ())
                    marker_rows += write_marker_rows(marker_handle, regions)
                    source_path = sparse_path(
                        sample=sample,
                        chrom=chrom,
                        hcc1395_sparse_root=hcc1395_sparse_root,
                        production_sparse_root=production_sparse_root,
                    )
                    shard_path = (
                        shard_dir / f"{sample}.{chrom}.candidate_read_join.tsv.gz"
                    )
                    source_rows, projected_rows, accumulators = (
                        write_sorted_candidate_shard(
                            sparse_calls=source_path,
                            sample=sample,
                            chrom=chrom,
                            regions=regions,
                            output_path=shard_path,
                            sort_binary=sort_path,
                            temporary_directory=staging_dir,
                        )
                    )
                    shard_count_rows, shard_gates = write_pattern_count_rows(
                        counts_handle, accumulators
                    )
                    count_rows += shard_count_rows
                    gate_totals.update(shard_gates)
                    source_rows_total += source_rows
                    projected_rows_total += projected_rows
                    relative_path = shard_path.relative_to(staging_dir).as_posix()
                    manifest_rows.append(
                        {
                            "sort_ordinal": sort_ordinal,
                            "dataset": sample,
                            "chrom": chrom,
                            "relative_path": relative_path,
                            "source_sparse_path": str(source_path.resolve()),
                            "n_source_rows": source_rows,
                            "n_projected_rows": projected_rows,
                            "sha256": sha256_path(shard_path),
                            "size_bytes": shard_path.stat().st_size,
                        }
                    )
                    sort_ordinal += 1

    with manifest_path.open("x", encoding="utf-8", newline="") as manifest_handle:
        manifest_handle.write("\t".join(MANIFEST_FIELDS) + "\n")
        for row in manifest_rows:
            write_tsv_row(manifest_handle, MANIFEST_FIELDS, row)

    receipt_path = staging_dir / "pattern_census.receipt.json"
    receipt: dict[str, object] = {
        "schema_name": "intersubmod.exact_raw_hp_pattern_census",
        "schema_version": "1.0.0",
        "scope": {
            "datasets": list(selected_samples),
            "chromosomes": list(selected_chromosomes),
            "task_type": "B_comprehensive_validation",
        },
        "contracts": {
            "topology_join_key": ["phase_set", "hp_family"],
            "analysis_stratum": "exact_hp_raw",
            "candidate_primary_sort": [
                "dataset",
                "chrom",
                "region_id",
                "hp_raw",
                "qname_sha256",
            ],
            "non_ra_call_normalization": "X",
            "formal_gate": {
                "minimum_active_bits": 2,
                "minimum_total_n": 40,
                "minimum_complete_n": 40,
                "minimum_complete_states": 2,
                "state_minimums": [5, 8, 10],
            },
        },
        "inputs": {
            "topology_jsonl": topology_sources,
            "sparse_tsv_count": len(manifest_rows),
        },
        "counts": {
            "topology_marker_rows": marker_rows,
            "sparse_source_rows": source_rows_total,
            "candidate_projected_rows": projected_rows_total,
            "pattern_count_rows": count_rows,
            **{key: int(gate_totals[key]) for key in sorted(gate_totals)},
        },
        "outputs": {
            "candidate_manifest": {
                "path": manifest_path.relative_to(staging_dir).as_posix(),
                "sha256": sha256_path(manifest_path),
                "size_bytes": manifest_path.stat().st_size,
            },
            "marker_universe": {
                "path": marker_path.relative_to(staging_dir).as_posix(),
                "sha256": sha256_path(marker_path),
                "size_bytes": marker_path.stat().st_size,
            },
            "pattern_counts": {
                "path": counts_path.relative_to(staging_dir).as_posix(),
                "sha256": sha256_path(counts_path),
                "size_bytes": counts_path.stat().st_size,
            },
        },
        "all_pass": True,
    }
    write_json(receipt_path, receipt)
    staging_dir.rename(output_dir)
    return receipt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build exact-raw-HP read-pattern census shards, counts, and marker universe."
        )
    )
    parser.add_argument(
        "--topology-root",
        type=Path,
        default=DEFAULT_TOPOLOGY_ROOT,
        help="Root containing samples/<sample>/<sample>.topology.jsonl.",
    )
    parser.add_argument(
        "--hcc1395-sparse-root",
        type=Path,
        default=DEFAULT_HCC1395_SPARSE_ROOT,
        help="HCC1395 exact-PS pilot root containing chromosomes/.",
    )
    parser.add_argument(
        "--production-sparse-root",
        type=Path,
        default=DEFAULT_PRODUCTION_SPARSE_ROOT,
        help="Six-dataset production root containing samples/.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    scope = parser.add_mutually_exclusive_group(required=True)
    scope.add_argument(
        "--all7",
        action="store_true",
        help="Process all seven technical datasets.",
    )
    scope.add_argument(
        "--sample",
        action="append",
        help="Dataset name; repeat or provide a comma-separated list.",
    )
    parser.add_argument(
        "--chrom",
        action="append",
        help="Autosome; repeat or provide a comma-separated list (default: chr1-22).",
    )
    parser.add_argument(
        "--sort-binary",
        default="sort",
        help="GNU-compatible external sort binary (default: sort).",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    samples = resolve_samples(all7=args.all7, sample_values=args.sample)
    chromosomes = resolve_chromosomes(args.chrom)
    receipt = build_pattern_census(
        topology_root=args.topology_root.resolve(),
        hcc1395_sparse_root=args.hcc1395_sparse_root.resolve(),
        production_sparse_root=args.production_sparse_root.resolve(),
        output_dir=args.output_dir.resolve(),
        samples=samples,
        chromosomes=chromosomes,
        sort_binary=args.sort_binary,
    )
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
