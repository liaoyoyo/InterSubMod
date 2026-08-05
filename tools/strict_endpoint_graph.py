#!/usr/bin/env python3
"""Build deterministic strict endpoint co-call graphs.

This module deliberately does *not* decide chromosome, phase-set, or
haplotype membership.  The caller must invoke it independently for exactly
one ``chromosome x exact PS x HP`` container.  It also does not use genomic
distance: ``site_index`` is an opaque, sortable node identifier.

An undirected edge receives one vote from a molecule only when that same
canonical molecule has a fixed ``R`` or ``A`` call at both endpoints.  Calls
must already be identity-deduplicated by the caller.  Duplicate molecule IDs
fail closed here so an alignment cannot silently contribute twice.

The connected components may be non-contiguous in site order.  This is
intentional: strict endpoint evidence, rather than interval/cut-span
coverage, defines connectivity.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
from itertools import combinations
from typing import Iterable, Mapping, Sequence


SCHEMA_NAME = "intersubmod.strict_endpoint_graph"
SCHEMA_VERSION = "1.0.0"
FIXED_CODES = frozenset({"R", "A"})
STATE_ORDER = ("RR", "RA", "AR", "AA")


def _require_plain_int(value: object, *, label: str) -> int:
    """Return an integer while rejecting bool and implicit coercion."""

    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an int, got {type(value).__name__}")
    return value


@dataclass(frozen=True, slots=True)
class MoleculeCalls:
    """One canonical molecule and its strictly sorted fixed endpoint calls."""

    molecule_id: str
    calls: Sequence[tuple[int, str]]

    def __post_init__(self) -> None:
        if not isinstance(self.molecule_id, str) or not self.molecule_id:
            raise ValueError("molecule_id must be a non-empty string")

        normalized = tuple(self.calls)
        previous_site: int | None = None
        for ordinal, call in enumerate(normalized):
            if not isinstance(call, tuple) or len(call) != 2:
                raise TypeError(
                    f"molecule {self.molecule_id!r}: call {ordinal} must be a (site_index, code) tuple"
                )
            site_index = _require_plain_int(
                call[0],
                label=f"molecule {self.molecule_id!r} call {ordinal} site_index",
            )
            code = call[1]
            if code not in FIXED_CODES:
                raise ValueError(
                    f"molecule {self.molecule_id!r} site {site_index}: code must be R or A, got {code!r}"
                )
            if previous_site is not None and site_index <= previous_site:
                raise ValueError(
                    f"molecule {self.molecule_id!r}: calls must have unique, strictly increasing site_index values"
                )
            previous_site = site_index

        object.__setattr__(self, "calls", normalized)


@dataclass(frozen=True, order=True, slots=True)
class PairSupport:
    """Unique-molecule support for one ordered endpoint pair."""

    site_i: int
    site_j: int
    total: int
    rr: int
    ra: int
    ar: int
    aa: int

    def __post_init__(self) -> None:
        site_i = _require_plain_int(self.site_i, label="site_i")
        site_j = _require_plain_int(self.site_j, label="site_j")
        if site_i >= site_j:
            raise ValueError("PairSupport endpoints must satisfy site_i < site_j")

        counts = {
            "total": self.total,
            "RR": self.rr,
            "RA": self.ra,
            "AR": self.ar,
            "AA": self.aa,
        }
        for label, value in counts.items():
            _require_plain_int(value, label=label)
            if value < 0:
                raise ValueError(f"{label} must be non-negative")
        if self.total <= 0:
            raise ValueError("total must be positive")
        if self.total != self.rr + self.ra + self.ar + self.aa:
            raise ValueError("state conservation failed: total != RR + RA + AR + AA")

    @property
    def state_counts(self) -> Mapping[str, int]:
        """Return the four-state support in contractual display order."""

        return {
            "RR": self.rr,
            "RA": self.ra,
            "AR": self.ar,
            "AA": self.aa,
        }

    def as_row(self) -> dict[str, int]:
        """Return the canonical serialization row."""

        return {
            "site_i": self.site_i,
            "site_j": self.site_j,
            "total": self.total,
            "RR": self.rr,
            "RA": self.ra,
            "AR": self.ar,
            "AA": self.aa,
        }


@dataclass(frozen=True, slots=True)
class ThresholdPartition:
    """A deterministic connected-component partition at one threshold."""

    threshold: int
    components: tuple[tuple[int, ...], ...]


class EndpointPairAccumulator:
    """Incrementally accumulate one exact chromosome x PS x HP container.

    The streaming form avoids retaining every molecule object until the end of
    a chromosome.  It has the same semantics as :func:`accumulate_pair_support`:
    one canonical molecule contributes at most one vote to each endpoint pair.
    """

    def __init__(self) -> None:
        self._seen_molecule_ids: set[str] = set()
        self._active_sites: set[int] = set()
        self._state_counts: dict[tuple[int, int], dict[str, int]] = defaultdict(
            lambda: {state: 0 for state in STATE_ORDER}
        )
        self._fixed_call_count = 0

    def add(self, molecule: MoleculeCalls) -> None:
        if not isinstance(molecule, MoleculeCalls):
            raise TypeError("molecule must be a MoleculeCalls instance")
        if molecule.molecule_id in self._seen_molecule_ids:
            raise ValueError(f"duplicate molecule_id: {molecule.molecule_id!r}")
        self._seen_molecule_ids.add(molecule.molecule_id)
        self._fixed_call_count += len(molecule.calls)
        self._active_sites.update(site_index for site_index, _ in molecule.calls)
        for (site_i, code_i), (site_j, code_j) in combinations(molecule.calls, 2):
            self._state_counts[(site_i, site_j)][code_i + code_j] += 1

    @property
    def molecule_count(self) -> int:
        return len(self._seen_molecule_ids)

    @property
    def fixed_call_count(self) -> int:
        return self._fixed_call_count

    @property
    def site_indices(self) -> tuple[int, ...]:
        return tuple(sorted(self._active_sites))

    def pair_support(self) -> tuple[PairSupport, ...]:
        rows: list[PairSupport] = []
        for (site_i, site_j), counts in sorted(self._state_counts.items()):
            rows.append(
                PairSupport(
                    site_i=site_i,
                    site_j=site_j,
                    total=sum(counts.values()),
                    rr=counts["RR"],
                    ra=counts["RA"],
                    ar=counts["AR"],
                    aa=counts["AA"],
                )
            )
        return tuple(rows)


def accumulate_pair_support(molecules: Iterable[MoleculeCalls]) -> tuple[PairSupport, ...]:
    """Count strict endpoint support once per unique canonical molecule.

    Every unordered pair among a molecule's fixed calls receives exactly one
    state vote.  The returned rows are sorted by ``(site_i, site_j)`` and are
    independent of input molecule order.
    """

    accumulator = EndpointPairAccumulator()
    for molecule in molecules:
        accumulator.add(molecule)
    return accumulator.pair_support()


def _normalized_sites(site_indices: Iterable[int]) -> tuple[int, ...]:
    sites = tuple(site_indices)
    for ordinal, site_index in enumerate(sites):
        _require_plain_int(site_index, label=f"site_indices[{ordinal}]")
    if len(set(sites)) != len(sites):
        raise ValueError("site_indices must be unique")
    return tuple(sorted(sites))


def _normalized_support(
    pair_support: Iterable[PairSupport],
    *,
    sites: tuple[int, ...],
) -> tuple[PairSupport, ...]:
    rows = tuple(pair_support)
    site_set = set(sites)
    seen_pairs: set[tuple[int, int]] = set()
    for row in rows:
        if not isinstance(row, PairSupport):
            raise TypeError("pair_support must contain PairSupport instances")
        key = (row.site_i, row.site_j)
        if key in seen_pairs:
            raise ValueError(f"duplicate pair support row: {key}")
        seen_pairs.add(key)
        if row.site_i not in site_set or row.site_j not in site_set:
            raise ValueError(f"pair support endpoint is absent from site_indices: {key}")
    return tuple(sorted(rows))


class _DisjointSet:
    """Small deterministic union-find implementation."""

    def __init__(self, sites: tuple[int, ...]) -> None:
        self._parent = {site: site for site in sites}
        self._rank = {site: 0 for site in sites}

    def find(self, site: int) -> int:
        parent = self._parent[site]
        if parent != site:
            self._parent[site] = self.find(parent)
        return self._parent[site]

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return

        left_rank = self._rank[left_root]
        right_rank = self._rank[right_root]
        if left_rank < right_rank or (left_rank == right_rank and left_root > right_root):
            left_root, right_root = right_root, left_root
            left_rank, right_rank = right_rank, left_rank
        self._parent[right_root] = left_root
        if left_rank == right_rank:
            self._rank[left_root] += 1


def connected_components(
    site_indices: Iterable[int],
    pair_support: Iterable[PairSupport],
    *,
    threshold: int,
) -> tuple[tuple[int, ...], ...]:
    """Return strict endpoint components using ``total >= threshold`` edges."""

    threshold = _require_plain_int(threshold, label="threshold")
    if threshold <= 0:
        raise ValueError("threshold must be positive")
    sites = _normalized_sites(site_indices)
    support = _normalized_support(pair_support, sites=sites)
    dsu = _DisjointSet(sites)
    for row in support:
        if row.total >= threshold:
            dsu.union(row.site_i, row.site_j)

    groups: dict[int, list[int]] = defaultdict(list)
    for site in sites:
        groups[dsu.find(site)].append(site)
    return tuple(sorted(tuple(group) for group in groups.values()))


def components_for_thresholds(
    site_indices: Iterable[int],
    pair_support: Iterable[PairSupport],
    *,
    thresholds: Iterable[int],
) -> tuple[ThresholdPartition, ...]:
    """Return nested component partitions in increasing threshold order."""

    sites = _normalized_sites(site_indices)
    support = _normalized_support(pair_support, sites=sites)
    normalized_thresholds: list[int] = []
    for ordinal, threshold in enumerate(thresholds):
        threshold = _require_plain_int(threshold, label=f"thresholds[{ordinal}]")
        if threshold <= 0:
            raise ValueError("thresholds must be positive")
        normalized_thresholds.append(threshold)
    if not normalized_thresholds:
        raise ValueError("at least one threshold is required")
    if len(set(normalized_thresholds)) != len(normalized_thresholds):
        raise ValueError("thresholds must be unique")

    partitions = tuple(
        ThresholdPartition(
            threshold=threshold,
            components=connected_components(sites, support, threshold=threshold),
        )
        for threshold in sorted(normalized_thresholds)
    )

    for lower, higher in zip(partitions, partitions[1:]):
        lower_owner = {
            site: component_index
            for component_index, component in enumerate(lower.components)
            for site in component
        }
        for component in higher.components:
            if len({lower_owner[site] for site in component}) != 1:
                raise RuntimeError(
                    "threshold nesting invariant failed: a higher-threshold component spans lower-threshold components"
                )
    return partitions


def strict_graph_digest(
    site_indices: Iterable[int],
    pair_support: Iterable[PairSupport],
) -> str:
    """Return a SHA-256 digest of the canonical graph evidence table."""

    sites = _normalized_sites(site_indices)
    support = _normalized_support(pair_support, sites=sites)
    payload = {
        "schema": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "site_indices": list(sites),
        "pair_support": [row.as_row() for row in support],
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()
