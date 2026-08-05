#!/usr/bin/env python3
"""Join InterSubMod read rows to the authoritative LongPhase-S HP/PS sidecar."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pysam


SIDECAR_COLUMNS = ("#CHROM", "START0", "END0", "QNAME", "FLAG", "MAPQ", "CIGAR_B2", "HP", "PS")
READ_COLUMNS = ("read_name", "chr", "start", "end", "mapq", "hp", "strand")
EXCLUDED_BY_INTERSUBMOD_FLAGS = 0x100 | 0x800 | 0x400 | 0x4
ProjectionKey = tuple[str, str, int, int, int, str]
FullKey = tuple[str, str, int, int, int, str]


@dataclass(frozen=True)
class LatestTags:
    hp: str
    ps: int | None


@dataclass(frozen=True)
class SidecarLookup:
    tags_by_projection: dict[ProjectionKey, LatestTags]
    full_identity_count_by_projection: dict[ProjectionKey, int]
    rows_fetched: int
    rows_eligible: int


def normalize_hp(value: str | None) -> str:
    return "0" if value in {None, "", "."} else str(value)


def normalize_ps(value: str | None) -> int | None:
    if value in {None, "", "."}:
        return None
    return int(value)


def strand_from_flag(flag: int) -> str:
    return "-" if flag & 16 else "+"


def cxx_read_name(qname: str, flag: int) -> str:
    if flag & 0x1:
        if flag & 0x40:
            return f"{qname}/1"
        if flag & 0x80:
            return f"{qname}/2"
    return qname


def projection_key(
    qname: str,
    chrom: str,
    start0: int,
    end0: int,
    mapq: int,
    strand: str,
) -> ProjectionKey:
    if strand not in {"+", "-"}:
        raise ValueError(f"Unexpected strand {strand!r}")
    return qname, chrom, int(start0), int(end0), int(mapq), strand


def read_projection(row: dict[str, str]) -> ProjectionKey:
    missing = [column for column in READ_COLUMNS if column not in row]
    if missing:
        raise ValueError(f"reads.tsv missing columns: {missing}")
    return projection_key(
        row["read_name"],
        row["chr"],
        int(row["start"]),
        int(row["end"]),
        int(row["mapq"]),
        row["strand"],
    )


def parse_sidecar_line(line: str) -> tuple[ProjectionKey, FullKey, LatestTags]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != len(SIDECAR_COLUMNS):
        raise ValueError(f"Malformed sidecar row with {len(fields)} fields")
    chrom, start, end, qname, flag, mapq, cigar_b2, hp, ps = fields
    flag_value = int(flag)
    projected = projection_key(
        cxx_read_name(qname, flag_value),
        chrom,
        int(start),
        int(end),
        int(mapq),
        strand_from_flag(flag_value),
    )
    full = (qname, chrom, int(start), int(end), flag_value, cigar_b2)
    return projected, full, LatestTags(normalize_hp(hp), normalize_ps(ps))


def fetch_site_lookup(tabix: pysam.TabixFile, chrom: str, pos1: int) -> SidecarLookup:
    """Fetch alignments overlapping one 1-based focal site and enforce tag uniqueness."""
    if pos1 < 1:
        raise ValueError(f"Position must be one-based and positive: {pos1}")
    by_projection: dict[ProjectionKey, LatestTags] = {}
    full_seen: dict[FullKey, LatestTags] = {}
    full_keys_by_projection: dict[ProjectionKey, set[FullKey]] = {}
    rows = 0
    eligible_rows = 0
    try:
        lines: Iterator[str] = tabix.fetch(chrom, pos1 - 1, pos1)
    except (ValueError, OSError) as error:
        raise RuntimeError(f"Cannot fetch sidecar at {chrom}:{pos1}: {error}") from error
    for line in lines:
        projected, full, tags = parse_sidecar_line(line)
        rows += 1
        if full[4] & EXCLUDED_BY_INTERSUBMOD_FLAGS:
            continue
        eligible_rows += 1
        if full in full_seen and full_seen[full] != tags:
            raise RuntimeError(f"Conflicting tags for exact sidecar identity: {full}")
        full_seen[full] = tags
        if projected in by_projection and by_projection[projected] != tags:
            raise RuntimeError(f"Ambiguous HP/PS after reads.tsv identity projection: {projected}")
        by_projection[projected] = tags
        full_keys_by_projection.setdefault(projected, set()).add(full)
    return SidecarLookup(
        tags_by_projection=by_projection,
        full_identity_count_by_projection={key: len(values) for key, values in full_keys_by_projection.items()},
        rows_fetched=rows,
        rows_eligible=eligible_rows,
    )


def iter_read_rows(path: Path) -> Iterator[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in READ_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} missing columns: {missing}")
        yield from reader


def join_read_rows(
    path: Path,
    lookup: SidecarLookup,
) -> list[tuple[dict[str, str], LatestTags, int]]:
    joined: list[tuple[dict[str, str], LatestTags, int]] = []
    for row in iter_read_rows(path):
        key = read_projection(row)
        if key not in lookup.tags_by_projection:
            raise KeyError(f"Latest LongPhase-S sidecar missing reads.tsv identity: {key}")
        joined.append((row, lookup.tags_by_projection[key], lookup.full_identity_count_by_projection[key]))
    return joined
