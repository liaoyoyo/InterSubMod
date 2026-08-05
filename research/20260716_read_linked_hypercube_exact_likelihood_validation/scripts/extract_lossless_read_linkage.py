#!/usr/bin/env python3
"""Extract sparse molecule×sSNV calls and read-linked boundary evidence.

The extractor uses ``samtools view -u -M -L`` so one coordinate-sorted primary
alignment overlapping any target sSNV is emitted once.  Every selected alignment
is then walked through CIGAR exactly once and merge-joined to the authoritative
LongPhase-S HP/PS sidecar by full alignment identity.

This is the M2 producer.  It preserves R/A/O/D/S/L/X call reasons and per-call
base quality, but deliberately does not choose evolutionary edges.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Sequence, TextIO

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lossless_read_contract import component_digest, linkage_components, molecule_id  # noqa: E402


AUTOSOMES = tuple(f"chr{value}" for value in range(1, 23))
DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
EXCLUDE_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
HP_VOCABULARY = frozenset({".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"})
DIAGNOSTIC_LINKAGE_BASES = (
    ("pooled", "pooled"),
    ("HP1", "1"),
    ("HP2", "2"),
    ("HP3", "3"),
    ("HP4", "4"),
    ("unphased", "none"),
)
PRIMARY_LINKAGE_BASES = ("PS_HP1", "PS_HP2")
MISSING_PS_LINKAGE_BASES = ("MISSING_PS_HP1", "MISSING_PS_HP2")
DECLARED_LINKAGE_BASES = tuple(basis for basis, _ in DIAGNOSTIC_LINKAGE_BASES) + PRIMARY_LINKAGE_BASES + MISSING_PS_LINKAGE_BASES
CALL_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "qname_sha256",
    "read_group",
    "alignment_id",
    "start0",
    "end0",
    "flag",
    "mapq",
    "strand",
    "hp_raw",
    "hp_family",
    "phase_set",
    "site_indices",
    "positions1",
    "call_codes",
    "base_qualities",
    "n_sites_in_span",
    "n_fixed_ra",
    "n_alt",
)


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos1: int
    ref: str
    alt: str


@dataclass(frozen=True)
class SidecarTags:
    hp: str
    ps: str | None


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def cigar_digest(cigar: str | None) -> str:
    return hashlib.blake2b((cigar or "*").encode(), digest_size=8).hexdigest()


def full_alignment_key(alignment: pysam.AlignedSegment) -> tuple[str, str, int, int, int, str]:
    start = int(alignment.reference_start)
    return (
        alignment.query_name,
        alignment.reference_name,
        start,
        int(alignment.reference_end or start + 1),
        int(alignment.flag),
        cigar_digest(alignment.cigarstring),
    )


def alignment_id(key: tuple[str, str, int, int, int, str]) -> str:
    return hashlib.sha256("\0".join(str(value) for value in key).encode()).hexdigest()


def hp_family(hp: str | None) -> str:
    value = "." if hp in {None, "", "."} else str(hp)
    if value not in HP_VOCABULARY:
        raise ValueError(f"unexpected LongPhase-S HP tag: {value!r}")
    if value in {"1", "1-1", "1-2"}:
        return "1"
    if value in {"2", "2-1", "2-2"}:
        return "2"
    if value in {"3", "4"}:
        return value
    return "none"


def parse_sidecar_line(line: str) -> tuple[tuple[str, str, int, int, int, str], SidecarTags]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 9:
        raise ValueError(f"malformed sidecar row with {len(fields)} columns")
    chrom, start, end, qname, flag, _mapq, cigar_b2, hp, ps = fields
    key = (qname, chrom, int(start), int(end), int(flag), cigar_b2)
    return key, SidecarTags(hp="." if hp in {"", "."} else hp, ps=None if ps in {"", "."} else ps)


class CoordinateSidecarJoiner:
    """Low-memory merge join for coordinate-sorted alignments and sidecar rows."""

    def __init__(self, lines: Iterator[str], chrom: str):
        self._lines = iter(lines)
        self.chrom = chrom
        self._pending: tuple[tuple[str, str, int, int, int, str], SidecarTags] | None = None
        self._group_start: int | None = None
        self._group: dict[tuple[str, str, int, int, int, str], SidecarTags] = {}
        self._last_query_start = -1
        self.rows_consumed = 0
        self.duplicate_identical = 0
        self.matched = 0
        self.missing = 0
        self._advance()

    def _advance(self) -> None:
        try:
            self._pending = parse_sidecar_line(next(self._lines))
            self.rows_consumed += 1
        except StopIteration:
            self._pending = None

    def _load_start(self, target: int) -> None:
        if target < self._last_query_start:
            raise RuntimeError("alignment stream is not coordinate sorted")
        self._last_query_start = target
        if self._group_start == target:
            return
        self._group_start = target
        self._group = {}
        while self._pending is not None and self._pending[0][2] < target:
            self._advance()
        while self._pending is not None and self._pending[0][2] == target:
            key, tags = self._pending
            if key[1] != self.chrom:
                raise RuntimeError(f"sidecar chromosome drift: {key[1]} != {self.chrom}")
            if key in self._group:
                if self._group[key] != tags:
                    raise RuntimeError(f"conflicting sidecar tags for exact identity: {key}")
                self.duplicate_identical += 1
            self._group[key] = tags
            self._advance()

    def lookup(self, alignment: pysam.AlignedSegment) -> SidecarTags:
        key = full_alignment_key(alignment)
        self._load_start(key[2])
        tags = self._group.get(key)
        if tags is None:
            self.missing += 1
            raise KeyError(f"sidecar missing exact alignment identity: {key}")
        self.matched += 1
        return tags


def load_variants(vcf_path: Path, chrom: str) -> list[Variant]:
    variants = []
    with pysam.VariantFile(str(vcf_path)) as source:
        try:
            records = source.fetch(chrom)
        except (ValueError, OSError):
            records = ()
        for record in records:
            alts = tuple(record.alts or ())
            if len(record.ref) != 1 or len(alts) != 1 or len(alts[0]) != 1:
                continue
            if set(record.filter.keys()) not in (set(), {"PASS"}):
                raise RuntimeError(f"non-PASS record in tree VCF: {chrom}:{record.pos}")
            variants.append(Variant(chrom, int(record.pos), record.ref.upper(), alts[0].upper()))
    variants.sort(key=lambda item: item.pos1)
    positions = [item.pos1 for item in variants]
    if len(positions) != len(set(positions)):
        raise RuntimeError(f"duplicate positional sSNV in {vcf_path} {chrom}")
    return variants


def call_alignment(
    alignment: pysam.AlignedSegment,
    variants: Sequence[Variant],
    positions0: Sequence[int],
    min_base_quality: int,
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[int | None, ...]]:
    """Call all target sSNVs inside one alignment span with one CIGAR walk."""
    start0 = int(alignment.reference_start)
    end0 = int(alignment.reference_end or start0)
    left = bisect.bisect_left(positions0, start0)
    right = bisect.bisect_left(positions0, end0)
    target_indices = tuple(range(left, right))
    if not target_indices:
        return (), (), ()
    calls = {index: "X" for index in target_indices}
    qualities: dict[int, int | None] = {index: None for index in target_indices}
    sequence = alignment.query_sequence or ""
    query_qualities = alignment.query_qualities
    ref_cursor = start0
    query_cursor = 0
    for operation, length in alignment.cigartuples or ():
        ref_consumes = operation in {0, 2, 3, 7, 8}
        query_consumes = operation in {0, 1, 4, 7, 8}
        if ref_consumes:
            block_end = ref_cursor + length
            block_left = bisect.bisect_left(positions0, ref_cursor, left, right)
            block_right = bisect.bisect_left(positions0, block_end, block_left, right)
            for index in range(block_left, block_right):
                if operation in {2, 3}:
                    calls[index] = "D" if operation == 2 else "S"
                    continue
                query_position = query_cursor + (positions0[index] - ref_cursor)
                if query_position < 0 or query_position >= len(sequence):
                    calls[index] = "X"
                    continue
                quality = None if query_qualities is None else int(query_qualities[query_position])
                qualities[index] = quality
                if quality is None or quality < min_base_quality:
                    calls[index] = "L"
                    continue
                base = sequence[query_position].upper()
                variant = variants[index]
                calls[index] = "R" if base == variant.ref else "A" if base == variant.alt else "O"
            ref_cursor = block_end
        if query_consumes:
            query_cursor += length
        if operation not in {0, 1, 2, 3, 4, 5, 6, 7, 8}:
            raise RuntimeError(f"unsupported CIGAR operation: {operation}")
    return (
        target_indices,
        tuple(calls[index] for index in target_indices),
        tuple(qualities[index] for index in target_indices),
    )


def apply_bridge(difference: list[int], fixed_indices: Sequence[int]) -> None:
    if len(fixed_indices) < 2:
        return
    left, right = min(fixed_indices), max(fixed_indices)
    difference[left] += 1
    difference[right] -= 1


def finish_difference(difference: Sequence[int], n_sites: int) -> tuple[int, ...]:
    if n_sites <= 1:
        return ()
    running = 0
    support = []
    for cut in range(n_sites - 1):
        running += int(difference[cut])
        support.append(running)
    return tuple(support)


def sparse_bridge_events(events: Counter[int], fixed_indices: Sequence[int]) -> None:
    """Accumulate one molecule span without allocating one array per phase set."""
    if len(fixed_indices) < 2:
        return
    left, right = min(fixed_indices), max(fixed_indices)
    events[left] += 1
    events[right] -= 1


def sparse_support_at_active_boundaries(
    active_indices: Sequence[int], events: Counter[int], n_sites: int
) -> tuple[int, ...]:
    """Return support between adjacent active sites in one HP-family x PS unit."""
    if tuple(active_indices) != tuple(sorted(set(active_indices))):
        raise ValueError("active site indices must be sorted and unique")
    if any(index < 0 or index >= n_sites for index in active_indices):
        raise ValueError("active site index outside chromosome catalog")
    if len(active_indices) <= 1:
        return ()
    event_items = sorted(
        (int(index), int(delta)) for index, delta in events.items() if int(delta) != 0
    )
    if any(index < 0 or index >= n_sites for index, _ in event_items):
        raise ValueError("bridge event index outside chromosome catalog")
    if sum(delta for _, delta in event_items) != 0:
        raise ValueError("bridge events must conserve interval mass")
    validation_running = 0
    for _, delta in event_items:
        validation_running += delta
        if validation_running < 0:
            raise ValueError("bridge support became negative")
    requested_cuts = tuple(active_indices[:-1])
    values = []
    running = 0
    event_offset = 0
    for cut in requested_cuts:
        while event_offset < len(event_items) and event_items[event_offset][0] <= cut:
            running += event_items[event_offset][1]
            if running < 0:
                raise ValueError("bridge support became negative")
            event_offset += 1
        values.append(running)
    return tuple(values)


def sparse_any_phase_support_at_active_boundaries(
    active_indices: Sequence[int],
    phase_events: Sequence[Counter[int]],
    thresholds: Sequence[int],
    n_sites: int,
) -> dict[int, tuple[bool, ...]]:
    """Return whether any phase set reaches each threshold at requested boundaries.

    Each phase set is swept only across its non-zero interval events.  Thresholded
    positive intervals are merged with another sparse difference sweep, avoiding
    a ``number_of_PS x chromosome_sites`` matrix for the legacy cross-PS audit.
    """
    requested_thresholds = tuple(sorted(set(int(value) for value in thresholds)))
    if any(value < 1 for value in requested_thresholds):
        raise ValueError("bridge thresholds must be positive")
    # Reuse the same validation contract even when no phase set has an event.
    if tuple(active_indices) != tuple(sorted(set(active_indices))):
        raise ValueError("active site indices must be sorted and unique")
    if any(index < 0 or index >= n_sites for index in active_indices):
        raise ValueError("active site index outside chromosome catalog")
    union_events = {threshold: Counter() for threshold in requested_thresholds}
    for events in phase_events:
        event_items = sorted(
            (int(index), int(delta)) for index, delta in events.items() if int(delta) != 0
        )
        if any(index < 0 or index >= n_sites for index, _ in event_items):
            raise ValueError("bridge event index outside chromosome catalog")
        if sum(delta for _, delta in event_items) != 0:
            raise ValueError("bridge events must conserve interval mass")
        running = 0
        for offset, (start, delta) in enumerate(event_items):
            running += delta
            if running < 0:
                raise ValueError("bridge support became negative")
            stop = event_items[offset + 1][0] if offset + 1 < len(event_items) else n_sites - 1
            if stop <= start or running == 0:
                continue
            for threshold in requested_thresholds:
                if running >= threshold:
                    union_events[threshold][start] += 1
                    union_events[threshold][stop] -= 1
    return {
        threshold: tuple(
            value > 0
            for value in sparse_support_at_active_boundaries(
                active_indices, union_events[threshold], n_sites
            )
        )
        for threshold in requested_thresholds
    }


def phase_set_token(phase_set: str | None) -> str:
    """Stable non-identifying token for component IDs; raw PS remains a column."""
    value = "MISSING" if phase_set is None else phase_set
    return hashlib.sha256(value.encode()).hexdigest()[:16]


def write_bed(path: Path, variants: Sequence[Variant]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for variant in variants:
            handle.write(f"{variant.chrom}\t{variant.pos1 - 1}\t{variant.pos1}\n")


def csv_value(value: Any) -> Any:
    if value is None:
        return ""
    return value


def input_entry(manifest: dict[str, Any], dataset: str) -> dict[str, Any]:
    matches = [entry for entry in manifest.get("samples", []) if entry.get("sample") == dataset]
    if len(matches) != 1:
        raise RuntimeError(f"manifest dataset cardinality mismatch: {dataset}")
    return matches[0]


def output_identity(path: Path) -> dict[str, Any]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256_path(path)}


def write_sha256_sidecar(path: Path) -> Path:
    """Write an exact-byte checksum outside the JSON it authenticates.

    A JSON document cannot contain the SHA-256 of its own final bytes without a
    self-reference problem.  The adjacent sidecar is therefore the canonical
    identity of ``receipt.json``.
    """
    checksum_path = path.with_name(f"{path.name}.sha256")
    checksum_path.write_text(f"{sha256_path(path)}  {path.name}\n", encoding="ascii")
    return checksum_path


def prepare_output_directory(
    output_dir: Path, *, require_existing_empty: bool = False
) -> None:
    """Create a fresh output directory or authenticate a preflight-created one.

    The explicit ``require_existing_empty`` mode is reserved for the direct
    release pilots: their launch-time resource receipt binds the inode before
    this producer starts.  Keeping the same empty directory preserves that
    binding while still failing closed on a symlink, file, missing path, or
    prior/partial output.
    """
    if require_existing_empty:
        if output_dir.is_symlink() or not output_dir.exists() or not output_dir.is_dir():
            raise FileExistsError(
                f"required preflight output directory is unavailable or not a real directory: {output_dir}"
            )
        if next(output_dir.iterdir(), None) is not None:
            raise FileExistsError(
                f"required preflight output directory is not empty: {output_dir}"
            )
        return
    if output_dir.is_symlink() or output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)


def extract_one(
    entry: dict[str, Any],
    chrom: str,
    output_dir: Path,
    *,
    mapq_min: int,
    baseq_min: int,
    thresholds: Sequence[int],
    samtools_threads: int,
    provenance: dict[str, Any] | None = None,
    require_existing_empty_output_dir: bool = False,
) -> dict[str, Any]:
    dataset = entry["sample"]
    if provenance is None:
        script_path = Path(__file__).resolve()
        provenance = {
            "extractor": {"path": str(script_path), "sha256": sha256_path(script_path)},
            "manifest": None,
        }
    bam_path = Path(entry["alignment_payload"]["path"])
    vcf_path = Path(entry["somatic"]["tree_vcf"]["path"])
    sidecar_path = Path(entry["read_tags"]["sidecar"]["path"])
    sidecar_index = Path(entry["read_tags"]["index"]["path"])
    for path in (bam_path, vcf_path, sidecar_path, sidecar_index):
        if not path.is_file():
            raise FileNotFoundError(path)
    variants = load_variants(vcf_path, chrom)
    if not variants:
        raise RuntimeError(f"no target variants: {dataset} {chrom}")
    positions1 = tuple(variant.pos1 for variant in variants)
    positions0 = tuple(position - 1 for position in positions1)
    prepare_output_directory(
        output_dir,
        require_existing_empty=require_existing_empty_output_dir,
    )
    bed_path = output_dir / f"{dataset}.{chrom}.targets.bed"
    calls_path = output_dir / f"{dataset}.{chrom}.molecule_sparse_calls.tsv.gz"
    sites_path = output_dir / f"{dataset}.{chrom}.site_catalog.tsv.gz"
    cuts_path = output_dir / f"{dataset}.{chrom}.cut_support.tsv.gz"
    components_path = output_dir / f"{dataset}.{chrom}.components.tsv.gz"
    membership_path = output_dir / f"{dataset}.{chrom}.site_component_membership.tsv.gz"
    stderr_path = output_dir / "samtools_view.stderr.log"
    write_bed(bed_path, variants)

    command = [
        "samtools",
        "view",
        "-u",
        "-M",
        "-L",
        str(bed_path),
        "-@",
        str(samtools_threads),
        "--no-PG",
        str(bam_path),
    ]
    difference = {family: [0] * len(variants) for family in ("pooled", "1", "2", "3", "4", "none")}
    phase_fixed_sites: dict[tuple[str, str | None], set[int]] = {}
    phase_bridge_events: dict[tuple[str, str | None], Counter[int]] = {}
    phase_molecule_counts: Counter[tuple[str, str | None]] = Counter()
    coverage = {code: [0] * len(variants) for code in ("R", "A", "O", "D", "S", "L", "X")}
    counts = Counter()
    raw_hp_counts = Counter()
    seen_alignment_ids: set[bytes] = set()
    seen_molecule_ids: set[bytes] = set()

    with pysam.TabixFile(str(sidecar_path), index=str(sidecar_index)) as sidecar, stderr_path.open("wb") as stderr_handle:
        joiner = CoordinateSidecarJoiner(iter(sidecar.fetch(chrom)), chrom)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr_handle)
        if process.stdout is None:
            raise RuntimeError("samtools stdout pipe unavailable")
        try:
            with pysam.AlignmentFile(process.stdout, "rb", check_sq=False) as stream, gzip.open(
                calls_path, "wt", encoding="utf-8", newline=""
            ) as calls_handle:
                writer = csv.DictWriter(calls_handle, CALL_FIELDS, delimiter="\t", extrasaction="raise")
                writer.writeheader()
                for alignment in stream.fetch(until_eof=True):
                    counts["raw_overlapping_alignments"] += 1
                    if alignment.is_secondary:
                        counts["alignment_class_secondary"] += 1
                    elif alignment.is_supplementary:
                        counts["alignment_class_supplementary"] += 1
                    elif alignment.is_unmapped:
                        counts["alignment_class_unmapped"] += 1
                    else:
                        counts["alignment_class_primary"] += 1
                    counts["flag_duplicate"] += bool(alignment.is_duplicate)
                    counts["flag_qcfail"] += bool(alignment.is_qcfail)
                    if alignment.flag & EXCLUDE_FLAGS:
                        counts["excluded_by_flag"] += 1
                        continue
                    if alignment.mapping_quality < mapq_min:
                        counts["mapq_rejected_after_flag"] += 1
                        continue
                    counts["canonical_eligible_alignments"] += 1
                    tags = joiner.lookup(alignment)
                    indices, codes, qualities = call_alignment(alignment, variants, positions0, baseq_min)
                    if not indices:
                        raise RuntimeError("samtools -L emitted alignment with no target site in reference span")
                    key = full_alignment_key(alignment)
                    aid = alignment_id(key)
                    aid_raw = bytes.fromhex(aid)
                    if aid_raw in seen_alignment_ids:
                        raise RuntimeError(f"duplicate alignment emitted by multi-region iterator: {key}")
                    seen_alignment_ids.add(aid_raw)
                    read_group = str(alignment.get_tag("RG")) if alignment.has_tag("RG") else "."
                    mid = molecule_id(dataset, read_group, alignment.query_name)
                    mid_raw = bytes.fromhex(mid)
                    if mid_raw in seen_molecule_ids:
                        raise RuntimeError(f"multiple canonical primary alignments for molecule: {dataset} {alignment.query_name}")
                    seen_molecule_ids.add(mid_raw)
                    family = hp_family(tags.hp)
                    raw_hp_counts[tags.hp] += 1
                    fixed = [index for index, code in zip(indices, codes) if code in {"R", "A"}]
                    apply_bridge(difference["pooled"], fixed)
                    apply_bridge(difference[family], fixed)
                    if family in {"1", "2"}:
                        phase_key = (family, tags.ps)
                        phase_molecule_counts[phase_key] += 1
                        phase_fixed_sites.setdefault(phase_key, set()).update(fixed)
                        sparse_bridge_events(
                            phase_bridge_events.setdefault(phase_key, Counter()), fixed
                        )
                        if tags.ps is None:
                            counts["missing_ps_hp12_molecule_rows"] += 1
                            counts[f"missing_ps_hp{family}_molecule_rows"] += 1
                            counts["missing_ps_hp12_fixed_ra_calls"] += len(fixed)
                        else:
                            counts["known_ps_hp12_molecule_rows"] += 1
                            counts[f"known_ps_hp{family}_molecule_rows"] += 1
                            counts["known_ps_hp12_fixed_ra_calls"] += len(fixed)
                    for index, code in zip(indices, codes):
                        coverage[code][index] += 1
                    counts["site_call_rows_sparse"] += len(indices)
                    counts["fixed_ra_calls"] += len(fixed)
                    counts["alt_calls"] += sum(code == "A" for code in codes)
                    counts["molecule_sparse_rows_written"] += 1
                    writer.writerow(
                        {
                            "dataset": dataset,
                            "chrom": chrom,
                            "molecule_id": mid,
                            "qname_sha256": hashlib.sha256(alignment.query_name.encode()).hexdigest(),
                            "read_group": read_group,
                            "alignment_id": aid,
                            "start0": key[2],
                            "end0": key[3],
                            "flag": key[4],
                            "mapq": int(alignment.mapping_quality),
                            "strand": "-" if alignment.is_reverse else "+",
                            "hp_raw": tags.hp,
                            "hp_family": family,
                            "phase_set": csv_value(tags.ps),
                            "site_indices": ",".join(str(index) for index in indices),
                            "positions1": ",".join(str(positions1[index]) for index in indices),
                            "call_codes": "".join(codes),
                            "base_qualities": ",".join("" if value is None else str(value) for value in qualities),
                            "n_sites_in_span": len(indices),
                            "n_fixed_ra": len(fixed),
                            "n_alt": sum(code == "A" for code in codes),
                        }
                    )
        finally:
            if process.stdout is not None:
                process.stdout.close()
            return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(f"samtools view failed with exit={return_code}; see {stderr_path}")

    support = {family: finish_difference(values, len(variants)) for family, values in difference.items()}
    with gzip.open(sites_path, "wt", encoding="utf-8", newline="") as handle:
        fields = ("site_index", "chrom", "pos1", "ref", "alt", "R", "A", "O", "D", "S", "L", "X")
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for index, variant in enumerate(variants):
            writer.writerow({
                "site_index": index,
                "chrom": chrom,
                "pos1": variant.pos1,
                "ref": variant.ref,
                "alt": variant.alt,
                **{code: coverage[code][index] for code in coverage},
            })

    cut_fields = (
        "cut_index",
        "left_pos1",
        "right_pos1",
        "gap_bp",
        "pooled_support",
        "hp1_support",
        "hp2_support",
        "hp3_support",
        "hp4_support",
        "unphased_support",
        "legacy_gap_le_50kb",
        *(f"boundary_at_threshold_{threshold}" for threshold in thresholds),
    )
    with gzip.open(cuts_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, cut_fields, delimiter="\t")
        writer.writeheader()
        for index, (left, right) in enumerate(zip(positions1, positions1[1:])):
            pooled = support["pooled"][index]
            writer.writerow(
                {
                    "cut_index": index,
                    "left_pos1": left,
                    "right_pos1": right,
                    "gap_bp": right - left,
                    "pooled_support": pooled,
                    "hp1_support": support["1"][index],
                    "hp2_support": support["2"][index],
                    "hp3_support": support["3"][index],
                    "hp4_support": support["4"][index],
                    "unphased_support": support["none"][index],
                    "legacy_gap_le_50kb": str(right - left <= 50_000).lower(),
                    **{
                        f"boundary_at_threshold_{threshold}": str(pooled < threshold).lower()
                        for threshold in thresholds
                    },
                }
            )

    component_fields = (
        "dataset",
        "chrom",
        "linkage_basis",
        "phase_set",
        "phase_set_status",
        "inference_role",
        "threshold",
        "component_index",
        "component_id",
        "start1",
        "end1",
        "k",
        "span_bp",
        "max_adjacent_gap_bp",
        "min_internal_bridge_support",
        "max_internal_bridge_support",
        "contains_gap_gt_50kb",
        "solver_route",
    )
    membership_fields = (
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
    component_rows: list[dict[str, Any]] = []
    membership_rows: list[dict[str, Any]] = []
    unit_conservation: dict[str, dict[str, dict[str, bool]]] = {}
    unit_summary: dict[str, dict[str, dict[str, dict[str, Any]]]] = {}
    unit_digests: dict[str, dict[str, dict[str, str]]] = {}
    aggregate_k_counts: dict[str, dict[str, Counter[int]]] = {
        basis: {str(threshold): Counter() for threshold in thresholds}
        for basis in DECLARED_LINKAGE_BASES
    }
    aggregate_component_metrics: dict[str, dict[str, Counter[str]]] = {
        basis: {str(threshold): Counter() for threshold in thresholds}
        for basis in DECLARED_LINKAGE_BASES
    }

    def add_component_unit(
        *,
        basis: str,
        phase_set: str | None,
        phase_status: str,
        inference_role: str,
        active_indices: Sequence[int],
        boundary_support: Sequence[int],
    ) -> None:
        active_indices = tuple(active_indices)
        active_positions = tuple(positions1[index] for index in active_indices)
        if len(boundary_support) != max(0, len(active_indices) - 1):
            raise AssertionError("active boundary support length mismatch")
        ps_key = "__MISSING__" if phase_set is None else str(phase_set)
        unit_conservation.setdefault(basis, {}).setdefault(ps_key, {})
        unit_summary.setdefault(basis, {}).setdefault(ps_key, {})
        unit_digests.setdefault(basis, {}).setdefault(ps_key, {})
        active_position_index = {position: index for index, position in enumerate(active_positions)}
        for threshold in thresholds:
            threshold_key = str(threshold)
            components = linkage_components(active_positions, boundary_support, threshold)
            conserved = tuple(position for component in components for position in component) == active_positions
            unit_conservation[basis][ps_key][threshold_key] = conserved
            k_counts = Counter(len(component) for component in components)
            summary = {
                "n_active_sites": len(active_positions),
                "n_components": len(components),
                "n_singletons_k1": k_counts.get(1, 0),
                "n_multisite_k_gt1": sum(value for key, value in k_counts.items() if key > 1),
                "n_effective_k_route_deferred_to_ranker": len(components),
                "max_k_component_sites": max(k_counts) if k_counts else 0,
                "k_component_sites_distribution": {
                    str(key): value for key, value in sorted(k_counts.items())
                },
            }
            unit_summary[basis][ps_key][threshold_key] = summary
            unit_digests[basis][ps_key][threshold_key] = component_digest(
                active_positions, boundary_support, threshold
            )
            aggregate_k_counts[basis][threshold_key].update(k_counts)
            aggregate_component_metrics[basis][threshold_key].update(
                {
                    "n_phase_set_units": 1,
                    "n_active_site_memberships": len(active_positions),
                    "n_components": len(components),
                    "n_singletons_k1": k_counts.get(1, 0),
                    "n_multisite_k_gt1": sum(value for key, value in k_counts.items() if key > 1),
                    "n_effective_k_route_deferred_to_ranker": len(components),
                }
            )
            for component_index, component in enumerate(components, start=1):
                selected_offsets = [active_position_index[position] for position in component]
                first_offset, last_offset = selected_offsets[0], selected_offsets[-1]
                internal_support = boundary_support[first_offset:last_offset]
                gaps = [right - left for left, right in zip(component, component[1:])]
                cid = (
                    f"{dataset}:{chrom}:{basis}:PS{phase_set_token(phase_set)}:B{threshold}:"
                    f"{component_index}:{component[0]}-{component[-1]}"
                )
                base = {
                    "dataset": dataset,
                    "chrom": chrom,
                    "linkage_basis": basis,
                    "phase_set": csv_value(phase_set),
                    "phase_set_status": phase_status,
                    "inference_role": inference_role,
                    "threshold": threshold,
                    "component_id": cid,
                }
                component_rows.append(
                    {
                        **base,
                        "component_index": component_index,
                        "start1": component[0],
                        "end1": component[-1],
                        "k": len(component),
                        "span_bp": component[-1] - component[0],
                        "max_adjacent_gap_bp": max(gaps) if gaps else 0,
                        "min_internal_bridge_support": min(internal_support) if internal_support else "",
                        "max_internal_bridge_support": max(internal_support) if internal_support else "",
                        "contains_gap_gt_50kb": str(any(gap > 50_000 for gap in gaps)).lower(),
                        "solver_route": "DEFER_TO_RANKER_EFFECTIVE_OBSERVED_ALT_K",
                    }
                )
                for position in component:
                    membership_rows.append(
                        {
                            **base,
                            "site_index": active_indices[active_position_index[position]],
                            "pos1": position,
                        }
                    )

    # Old pooled/HP-family bases are retained strictly as non-primary diagnostics.
    for basis, support_key in DIAGNOSTIC_LINKAGE_BASES:
        add_component_unit(
            basis=basis,
            phase_set=None,
            phase_status="MIXED_OR_NOT_APPLICABLE_DIAGNOSTIC",
            inference_role="DIAGNOSTIC_NONPRIMARY",
            active_indices=tuple(range(len(positions1))),
            boundary_support=support[support_key],
        )

    for (family, phase_set), active in sorted(
        phase_fixed_sites.items(), key=lambda item: (item[0][0], "" if item[0][1] is None else item[0][1])
    ):
        active_indices = tuple(sorted(active))
        if not active_indices:
            continue
        events = phase_bridge_events.get((family, phase_set), Counter())
        boundary_support = sparse_support_at_active_boundaries(
            active_indices, events, len(positions1)
        )
        if phase_set is None:
            basis = f"MISSING_PS_HP{family}"
            status = "MISSING_PS_SENSITIVITY_ABSTAIN"
            role = "SENSITIVITY_ABSTAIN_NONPRIMARY"
        else:
            basis = f"PS_HP{family}"
            status = "KNOWN_PS_PRIMARY"
            role = "PRIMARY_PS_AWARE"
        add_component_unit(
            basis=basis,
            phase_set=phase_set,
            phase_status=status,
            inference_role=role,
            active_indices=active_indices,
            boundary_support=boundary_support,
        )

    with gzip.open(components_path, "wt", encoding="utf-8", newline="") as component_handle, gzip.open(
        membership_path, "wt", encoding="utf-8", newline=""
    ) as membership_handle:
        component_writer = csv.DictWriter(component_handle, component_fields, delimiter="\t")
        membership_writer = csv.DictWriter(membership_handle, membership_fields, delimiter="\t")
        component_writer.writeheader()
        membership_writer.writeheader()
        component_writer.writerows(component_rows)
        membership_writer.writerows(membership_rows)

    component_summary: dict[str, dict[str, dict[str, Any]]] = {}
    for basis in DECLARED_LINKAGE_BASES:
        component_summary[basis] = {}
        for threshold in thresholds:
            key = str(threshold)
            metrics = dict(aggregate_component_metrics[basis][key])
            k_counts = aggregate_k_counts[basis][key]
            component_summary[basis][key] = {
                **metrics,
                "max_k_component_sites": max(k_counts) if k_counts else 0,
                "k_component_sites_distribution": {
                    str(k): value for k, value in sorted(k_counts.items())
                },
                # Compatibility aliases are component-site k, never effective solver k.
                "max_k": max(k_counts) if k_counts else 0,
                "k_distribution": {str(k): value for k, value in sorted(k_counts.items())},
            }

    legacy_cross_ps: dict[str, dict[str, dict[str, int]]] = {}
    for family in ("1", "2"):
        active = tuple(sorted(
            set().union(*(sites for (fam, _), sites in phase_fixed_sites.items() if fam == family))
            if any(fam == family for fam, _ in phase_fixed_sites)
            else set()
        ))
        legacy_cross_ps[family] = {}
        known_ps_linked = sparse_any_phase_support_at_active_boundaries(
            active,
            [
                events for (fam, ps), events in phase_bridge_events.items()
                if fam == family and ps is not None
            ],
            thresholds,
            len(positions1),
        )
        active_cut_indices = active[:-1]
        for threshold in thresholds:
            threshold_known = known_ps_linked[threshold]
            if len(threshold_known) != len(active_cut_indices):
                raise AssertionError("known-PS legacy boundary audit length mismatch")
            legacy_linked = 0
            any_known_linked = 0
            aggregate_only = 0
            for cut, known in zip(active_cut_indices, threshold_known):
                old = support[family][cut]
                legacy_linked += old >= threshold
                any_known_linked += known
                aggregate_only += old >= threshold and not known
            active_positions = tuple(positions1[index] for index in active)
            legacy_active_support = tuple(support[family][index] for index in active[:-1])
            legacy_components = linkage_components(active_positions, legacy_active_support, threshold)
            known_ps_component_sum = sum(
                int(values[str(threshold)]["n_components"])
                for ps, values in unit_summary.get(f"PS_HP{family}", {}).items()
            )
            active_by_known_ps = {
                ps: {positions1[index] for index in sites}
                for (fam, ps), sites in phase_fixed_sites.items()
                if fam == family and ps is not None
            }
            mixed_components = sum(
                sum(bool(set(component) & sites) for sites in active_by_known_ps.values()) > 1
                for component in legacy_components
            )
            legacy_cross_ps[family][str(threshold)] = {
                "active_family_boundary_denominator": max(0, len(active) - 1),
                "legacy_family_linked_boundaries": legacy_linked,
                "boundaries_linked_within_at_least_one_known_ps": any_known_linked,
                "legacy_boundaries_supported_only_after_cross_ps_or_missing_ps_aggregation": aggregate_only,
                "legacy_active_site_components": len(legacy_components),
                "legacy_components_spanning_multiple_known_phase_sets": mixed_components,
                "known_ps_components_sum": known_ps_component_sum,
                "known_ps_minus_legacy_component_count": known_ps_component_sum - len(legacy_components),
                "unique_family_active_sites": len(active),
                "known_ps_active_site_memberships": sum(len(sites) for sites in active_by_known_ps.values()),
            }
    component_digests_by_basis: dict[str, dict[str, str]] = {}
    for basis, phase_values in unit_digests.items():
        component_digests_by_basis[basis] = {}
        for threshold in thresholds:
            key = str(threshold)
            if set(phase_values) == {"__MISSING__"}:
                component_digests_by_basis[basis][key] = phase_values["__MISSING__"][key]
            else:
                component_digests_by_basis[basis][key] = hashlib.sha256(
                    json.dumps(
                        {ps: values[key] for ps, values in sorted(phase_values.items())},
                        sort_keys=True,
                        separators=(",", ":"),
                    ).encode()
                ).hexdigest()

    outputs = [bed_path, calls_path, sites_path, cuts_path, components_path, membership_path, stderr_path]
    call_reason_mass = sum(sum(values) for values in coverage.values())
    fixed_reason_mass = sum(coverage[code][index] for code in ("R", "A") for index in range(len(variants)))
    alt_reason_mass = sum(coverage["A"])
    all_basis_thresholds_conserve = all(
        conserved
        for basis_values in unit_conservation.values()
        for phase_values in basis_values.values()
        for conserved in phase_values.values()
    )
    receipt_path = output_dir / "receipt.json"
    checksum_path = receipt_path.with_name(f"{receipt_path.name}.sha256")
    receipt = {
        "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
        "schema_version": "1.2.0",
        "scope": {"dataset": dataset, "chrom": chrom, "n_sSNV": len(variants)},
        "provenance": provenance,
        "parameters": {
            "mapq_min": mapq_min,
            "baseq_min": baseq_min,
            "samtools_threads": samtools_threads,
            "excluded_flags_decimal": EXCLUDE_FLAGS,
            "excluded_flags_hex": hex(EXCLUDE_FLAGS),
            "bridge_thresholds": list(thresholds),
            "bridge_definition": "unique canonical molecule with >=1 fixed R/A site on each side of cut",
            "legacy_50kb_role": "sensitivity_only_not_hard_boundary",
            "component_linkage_bases": list(DECLARED_LINKAGE_BASES),
            "primary_component_linkage_bases": list(PRIMARY_LINKAGE_BASES),
            "phase_set_contract": (
                "primary unit = HP family x exact known PS x read-linked component x threshold; "
                "missing PS never enters primary; pooled/legacy HP bases are diagnostic only"
            ),
            "component_active_site_policy": (
                "PS-aware and missing-PS units contain only sites with >=1 fixed R/A call in that exact unit"
            ),
            "effective_k_policy": (
                "extractor reports component-site k only; observed-ALT effective k and exact-solver route "
                "are determined downstream after structural exact-pattern minread"
            ),
        },
        "command": command,
        "inputs": {
            "bam": entry["alignment_payload"],
            "tree_vcf": entry["somatic"]["tree_vcf"],
            "read_tag_sidecar": entry["read_tags"]["sidecar"],
            "read_tag_sidecar_index": entry["read_tags"]["index"],
        },
        "counts": {
            **dict(counts),
            "unique_alignment_ids": len(seen_alignment_ids),
            "unique_molecule_ids": len(seen_molecule_ids),
            "sidecar_rows_consumed_through_last_target_start": joiner.rows_consumed,
            "sidecar_exact_matches": joiner.matched,
            "sidecar_missing": joiner.missing,
            "sidecar_duplicate_identical": joiner.duplicate_identical,
            "raw_HP_counts": dict(raw_hp_counts),
        },
        "phase_set_contract_counts": {
            "known_phase_sets_by_hp_family": {
                family: len({ps for fam, ps in phase_fixed_sites if fam == family and ps is not None})
                for family in ("1", "2")
            },
            "known_ps_active_site_memberships_by_hp_family": {
                family: sum(
                    len(sites) for (fam, ps), sites in phase_fixed_sites.items()
                    if fam == family and ps is not None
                )
                for family in ("1", "2")
            },
            "missing_ps_active_sites_by_hp_family": {
                family: len(phase_fixed_sites.get((family, None), set())) for family in ("1", "2")
            },
            "molecule_rows_by_hp_family_phase_set": {
                f"HP{family}|{'MISSING' if ps is None else ps}": value
                for (family, ps), value in sorted(
                    phase_molecule_counts.items(),
                    key=lambda item: (item[0][0], "" if item[0][1] is None else item[0][1]),
                )
            },
        },
        "legacy_cross_phase_set_aggregation_audit": legacy_cross_ps,
        "component_conservation_by_linkage_basis_phase_set": unit_conservation,
        "component_summary_by_linkage_basis": component_summary,
        "component_summary_by_linkage_basis_phase_set": unit_summary,
        "component_digests_by_linkage_basis_phase_set": unit_digests,
        "component_digests_by_linkage_basis": component_digests_by_basis,
        # Backward-compatible aliases are explicitly pooled-only.  Tree solvers
        # must use the linkage-basis table above for per-HP inference.
        "component_conservation_by_threshold": unit_conservation["pooled"]["__MISSING__"],
        "component_summary_by_threshold": component_summary["pooled"],
        "component_digests_by_threshold": unit_digests["pooled"]["__MISSING__"],
        "outputs": {path.name: output_identity(path) for path in outputs},
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": checksum_path.name,
            "covers": receipt_path.name,
        },
        "checks": {
            "samtools_exit_zero": return_code == 0,
            "raw_alignment_class_conserved": counts["raw_overlapping_alignments"]
            == counts["alignment_class_primary"]
            + counts["alignment_class_secondary"]
            + counts["alignment_class_supplementary"]
            + counts["alignment_class_unmapped"],
            "raw_filter_funnel_conserved": counts["raw_overlapping_alignments"]
            == counts["excluded_by_flag"]
            + counts["mapq_rejected_after_flag"]
            + counts["canonical_eligible_alignments"],
            "every_eligible_alignment_exact_sidecar_joined": joiner.matched
            == counts["canonical_eligible_alignments"]
            and joiner.missing == 0,
            "multi_region_alignment_unique": len(seen_alignment_ids) == counts["canonical_eligible_alignments"],
            "canonical_molecule_unique": len(seen_molecule_ids) == counts["canonical_eligible_alignments"],
            "eligible_alignment_rows_written": counts["molecule_sparse_rows_written"]
            == counts["canonical_eligible_alignments"],
            "site_call_reason_mass_conserved": call_reason_mass == counts["site_call_rows_sparse"],
            "fixed_ra_mass_conserved": fixed_reason_mass == counts["fixed_ra_calls"],
            "alt_mass_conserved": alt_reason_mass == counts["alt_calls"],
            "all_linkage_units_and_thresholds_conserve_their_active_sites": all_basis_thresholds_conserve,
            "primary_components_only_use_known_ps": all(
                row["phase_set_status"] == "KNOWN_PS_PRIMARY" and row["phase_set"] not in {None, ""}
                for row in component_rows if row["inference_role"] == "PRIMARY_PS_AWARE"
            ),
            "missing_ps_never_primary": all(
                row["inference_role"] != "PRIMARY_PS_AWARE"
                for row in component_rows if row["phase_set"] in {None, ""}
            ),
            "primary_membership_is_active_fixed_ra_only": all(
                int(row["site_index"]) in phase_fixed_sites[(
                    "1" if row["linkage_basis"] == "PS_HP1" else "2",
                    str(row["phase_set"]),
                )]
                for row in membership_rows if row["inference_role"] == "PRIMARY_PS_AWARE"
            ),
            "site_catalog_cardinality": len(variants) == len(positions1),
        },
    }
    receipt["all_pass"] = all(receipt["checks"].values())
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_sha256_sidecar(receipt_path)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample", required=True, choices=DATASETS)
    parser.add_argument("--chrom", required=True, choices=AUTOSOMES)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mapq-min", type=int, default=20)
    parser.add_argument("--baseq-min", type=int, default=20)
    parser.add_argument("--bridge-thresholds", default="1,2,3,5")
    parser.add_argument("--samtools-threads", type=int, default=1)
    parser.add_argument(
        "--require-existing-empty-output-dir",
        action="store_true",
        help=(
            "use a real empty output directory already created and inode-bound "
            "by the release-pilot resource gate"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    thresholds = tuple(sorted({int(value) for value in args.bridge_thresholds.split(",")}))
    if args.mapq_min < 0 or args.baseq_min < 0 or args.samtools_threads < 0 or not thresholds or min(thresholds) < 1:
        raise ValueError("invalid quality, thread, or threshold parameter")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    script_path = Path(__file__).resolve()
    provenance = {
        "extractor": {"path": str(script_path), "sha256": sha256_path(script_path)},
        "manifest": {"path": str(args.manifest.resolve()), "sha256": sha256_path(args.manifest)},
    }
    receipt = extract_one(
        input_entry(manifest, args.sample),
        args.chrom,
        args.output_dir,
        mapq_min=args.mapq_min,
        baseq_min=args.baseq_min,
        thresholds=thresholds,
        samtools_threads=args.samtools_threads,
        provenance=provenance,
        require_existing_empty_output_dir=args.require_existing_empty_output_dir,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    raise SystemExit(0 if receipt["all_pass"] else 1)


if __name__ == "__main__":
    main()
