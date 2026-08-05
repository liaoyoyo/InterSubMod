#!/usr/bin/env python3
"""Lossless molecule-call and read-linked boundary primitives.

The module is storage/backend agnostic.  BAM, mpileup, or existing read tables
must first produce one canonical row per ``molecule × site``.  These helpers then
encode partial patterns without completion expansion and calculate molecule-level
cut support without constructing all site pairs.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from typing import Iterable, Sequence


CALL_CODES = frozenset({"R", "A", "O", "D", "S", "L", "X"})
FIXED_CODES = frozenset({"R", "A"})


def molecule_id(dataset_uuid: str, read_group: str | None, qname: str) -> str:
    """Stable molecule identity scoped to one dataset and read group."""
    if not dataset_uuid or not qname:
        raise ValueError("dataset_uuid and qname are required")
    payload = "\0".join((dataset_uuid, read_group or ".", qname))
    return hashlib.sha256(payload.encode()).hexdigest()


@dataclass(frozen=True)
class SiteCall:
    dataset: str
    chrom: str
    pos1: int
    molecule_id: str
    call: str
    hp_family: str
    phase_set: str | None = None
    base_quality: int | None = None
    mapping_quality: int | None = None

    def __post_init__(self) -> None:
        if self.pos1 < 1:
            raise ValueError("pos1 must be positive and one-based")
        if self.call not in CALL_CODES:
            raise ValueError(f"invalid call code: {self.call!r}")
        if self.hp_family not in {"1", "2", "3", "4", "none"}:
            raise ValueError(f"invalid HP family: {self.hp_family!r}")
        if self.base_quality is not None and self.base_quality < 0:
            raise ValueError("base quality cannot be negative")
        if self.mapping_quality is not None and self.mapping_quality < 0:
            raise ValueError("mapping quality cannot be negative")


@dataclass(frozen=True)
class EncodedPattern:
    k: int
    pattern: str
    ref_mask: int
    alt_mask: int
    other_mask: int
    deletion_mask: int
    refskip_mask: int
    low_baseq_mask: int
    uncovered_mask: int
    fixed_mask: int
    free_mask: int

    @classmethod
    def from_codes(cls, codes: Sequence[str]) -> "EncodedPattern":
        if not codes:
            raise ValueError("codes must not be empty")
        invalid = sorted(set(codes) - CALL_CODES)
        if invalid:
            raise ValueError(f"invalid call codes: {invalid}")
        masks = {code: 0 for code in CALL_CODES}
        for bit, code in enumerate(codes):
            masks[code] |= 1 << bit
        k = len(codes)
        all_mask = (1 << k) - 1
        fixed = masks["R"] | masks["A"]
        free = all_mask ^ fixed
        return cls(
            k=k,
            pattern="".join(code if code in FIXED_CODES else "X" for code in codes),
            ref_mask=masks["R"],
            alt_mask=masks["A"],
            other_mask=masks["O"],
            deletion_mask=masks["D"],
            refskip_mask=masks["S"],
            low_baseq_mask=masks["L"],
            uncovered_mask=masks["X"],
            fixed_mask=fixed,
            free_mask=free,
        )

    @property
    def n_free(self) -> int:
        return bin(self.free_mask).count("1")

    @property
    def n_completions(self) -> int:
        return 1 << self.n_free

    def compatible(self, vertex: int) -> bool:
        if vertex < 0 or vertex >= (1 << self.k):
            return False
        return ((vertex ^ self.alt_mask) & self.fixed_mask) == 0

    def to_dict(self) -> dict:
        return asdict(self)


def validate_unique_site_calls(calls: Iterable[SiteCall]) -> list[SiteCall]:
    rows = list(calls)
    seen = {}
    for row in rows:
        key = (row.dataset, row.chrom, row.pos1, row.molecule_id)
        if key in seen:
            raise ValueError(f"duplicate molecule-site row: {key}")
        seen[key] = row
    return rows


def molecule_patterns(
    calls: Iterable[SiteCall],
    positions: Sequence[int],
) -> dict[tuple[str, str, str | None], EncodedPattern]:
    """Create one mutually exclusive pattern per molecule×HP×PS.

    Missing rows are encoded as ``X``.  Explicit O/D/S/L reasons remain in the
    detailed masks while the solver-compatible display uses X.
    """
    if not positions or list(positions) != sorted(set(positions)):
        raise ValueError("positions must be nonempty, sorted, and unique")
    rows = validate_unique_site_calls(calls)
    position_index = {position: index for index, position in enumerate(positions)}
    grouped: dict[tuple[str, str, str | None], list[str]] = {}
    for row in rows:
        if row.pos1 not in position_index:
            raise ValueError(f"site call outside requested positions: {row.pos1}")
        key = (row.molecule_id, row.hp_family, row.phase_set)
        codes = grouped.setdefault(key, ["X"] * len(positions))
        codes[position_index[row.pos1]] = row.call
    return {key: EncodedPattern.from_codes(codes) for key, codes in grouped.items()}


def aggregate_symbolic_patterns(
    calls: Iterable[SiteCall],
    positions: Sequence[int],
) -> list[dict]:
    patterns = molecule_patterns(calls, positions)
    counts = Counter(
        (
            hp,
            ps,
            encoded.pattern,
            encoded.ref_mask,
            encoded.alt_mask,
            encoded.other_mask,
            encoded.deletion_mask,
            encoded.refskip_mask,
            encoded.low_baseq_mask,
            encoded.uncovered_mask,
        )
        for (_, hp, ps), encoded in patterns.items()
    )
    rows = []
    for key, count in sorted(counts.items(), key=lambda item: tuple(str(value) for value in item[0])):
        hp, ps, pattern, ref, alt, other, deletion, refskip, low, uncovered = key
        encoded = EncodedPattern.from_codes(list(pattern))
        rows.append(
            {
                "hp_family": hp,
                "phase_set": ps,
                "pattern_rax": pattern,
                "ref_mask": ref,
                "alt_mask": alt,
                "other_mask": other,
                "deletion_mask": deletion,
                "refskip_mask": refskip,
                "low_baseq_mask": low,
                "uncovered_mask": uncovered,
                "fixed_mask": ref | alt,
                "free_mask": ((1 << len(positions)) - 1) ^ (ref | alt),
                "n_free": encoded.n_free,
                "n_completions": 1 << encoded.n_free,
                "n_molecules": count,
            }
        )
    if sum(row["n_molecules"] for row in rows) != len(patterns):
        raise AssertionError("molecule pattern mass conservation failed")
    return rows


def cut_bridge_support(
    calls: Iterable[SiteCall],
    positions: Sequence[int],
    *,
    hp_family: str | None = None,
) -> tuple[int, ...]:
    """Count unique molecules with a fixed R/A observation on both sides of each cut.

    A molecule whose leftmost and rightmost fixed sites are indices ``lo`` and
    ``hi`` contributes exactly once to every cut ``lo .. hi-1``.  Difference-array
    accumulation makes this O(number of calls + number of sites), rather than a
    pairwise O(k²) site graph.
    """
    if not positions or list(positions) != sorted(set(positions)):
        raise ValueError("positions must be nonempty, sorted, and unique")
    index = {position: offset for offset, position in enumerate(positions)}
    fixed_by_molecule: dict[str, set[int]] = defaultdict(set)
    for row in validate_unique_site_calls(calls):
        if row.pos1 not in index:
            raise ValueError(f"site call outside requested positions: {row.pos1}")
        if hp_family is not None and row.hp_family != hp_family:
            continue
        if row.call in FIXED_CODES:
            fixed_by_molecule[row.molecule_id].add(index[row.pos1])
    if len(positions) == 1:
        return ()
    difference = [0] * len(positions)
    for indices in fixed_by_molecule.values():
        if len(indices) < 2:
            continue
        lo, hi = min(indices), max(indices)
        difference[lo] += 1
        difference[hi] -= 1
    support = []
    running = 0
    for cut in range(len(positions) - 1):
        running += difference[cut]
        support.append(running)
    return tuple(support)


def linkage_components(
    positions: Sequence[int],
    cut_support: Sequence[int],
    threshold: int,
) -> tuple[tuple[int, ...], ...]:
    """Split only where molecule bridge support is below the declared threshold."""
    if threshold < 1:
        raise ValueError("threshold must be positive")
    if positions and list(positions) != sorted(set(positions)):
        raise ValueError("positions must be sorted and unique")
    if len(cut_support) != max(0, len(positions) - 1):
        raise ValueError("cut support length mismatch")
    if not positions:
        return ()
    components = []
    current = [positions[0]]
    for cut, next_position in zip(cut_support, positions[1:]):
        if int(cut) < threshold:
            components.append(tuple(current))
            current = []
        current.append(next_position)
    components.append(tuple(current))
    if tuple(position for component in components for position in component) != tuple(positions):
        raise AssertionError("linkage component conservation failed")
    return tuple(components)


def component_digest(positions: Sequence[int], cut_support: Sequence[int], threshold: int) -> str:
    payload = {
        "positions": list(positions),
        "cut_support": list(cut_support),
        "threshold": int(threshold),
        "components": [list(component) for component in linkage_components(positions, cut_support, threshold)],
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
