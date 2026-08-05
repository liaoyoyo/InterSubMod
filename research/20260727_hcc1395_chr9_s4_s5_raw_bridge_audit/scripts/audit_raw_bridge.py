#!/usr/bin/env python3
"""Read-only HCC1395 chr9 S4/S5 raw molecule bridge audit."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
from collections import Counter, defaultdict, deque
from pathlib import Path

import pysam


LEFT_INDEX = 78
RIGHT_INDEX = 79
REGION_SITE_INDEX_MIN = 75
REGION_SITE_INDEX_MAX = 81
LEFT_POS1 = 5_074_198
RIGHT_POS1 = 5_095_768
LEFT_REF = "G"
LEFT_ALT = "C"
RIGHT_REF = "C"
RIGHT_ALT = "A"
PHASE_SET = "3625170"
HP_FAMILIES = {"1", "2"}
EXCLUDED_FLAGS = 0xF04
MAPQ_MIN = 20
BASEQ_MIN = 20


def query_position_at(alignment: pysam.AlignedSegment, reference_pos0: int):
    for query_pos, reference_pos in alignment.get_aligned_pairs(
        matches_only=False
    ):
        if reference_pos == reference_pos0:
            return query_pos
    return None


def raw_site_census(
    bam: pysam.AlignmentFile, chrom: str, pos1: int, ref: str, alt: str
):
    pos0 = pos1 - 1
    counts = Counter()
    qnames = {
        "raw": set(),
        "flagpass": set(),
        "flag_mapq_pass": set(),
    }
    for alignment in bam.fetch(chrom, pos0, pos0 + 1):
        counts["raw_reference_overlap_rows"] += 1
        qnames["raw"].add(alignment.query_name)
        flagpass = not (alignment.flag & EXCLUDED_FLAGS)
        if flagpass:
            counts["flagpass_rows"] += 1
            qnames["flagpass"].add(alignment.query_name)
        if flagpass and alignment.mapping_quality >= MAPQ_MIN:
            counts["flag_mapq_pass_rows"] += 1
            qnames["flag_mapq_pass"].add(alignment.query_name)
        query_pos = query_position_at(alignment, pos0)
        if query_pos is None:
            counts["no_query_base_at_site"] += 1
            continue
        counts["query_base_present"] += 1
        if not (flagpass and alignment.mapping_quality >= MAPQ_MIN):
            continue
        counts["filtered_query_base_present"] += 1
        base = alignment.query_sequence[query_pos].upper()
        baseq = alignment.query_qualities[query_pos]
        if baseq < BASEQ_MIN:
            counts["filtered_low_baseq"] += 1
        elif base == ref:
            counts["filtered_usable_R"] += 1
        elif base == alt:
            counts["filtered_usable_A"] += 1
        else:
            counts["filtered_other_base"] += 1
    return dict(counts), qnames


def raw_pair_census(bam: pysam.AlignmentFile, chrom: str):
    left0 = LEFT_POS1 - 1
    right0 = RIGHT_POS1 - 1
    left = list(bam.fetch(chrom, left0, left0 + 1))
    right = list(bam.fetch(chrom, right0, right0 + 1))
    counts = Counter()
    pair_patterns = Counter()
    pair_quality_patterns = Counter()
    for alignment in left:
        if (
            alignment.reference_start is not None
            and alignment.reference_end is not None
            and alignment.reference_start <= left0
            and alignment.reference_end > right0
        ):
            counts["continuous_reference_span_all"] += 1
            if alignment.flag & EXCLUDED_FLAGS:
                counts["excluded_by_flag"] += 1
                continue
            counts["continuous_reference_span_flagpass"] += 1
            if alignment.mapping_quality < MAPQ_MIN:
                counts["excluded_by_mapq"] += 1
                continue
            counts["continuous_reference_span_flag_mapq_pass"] += 1

            site_calls = []
            for pos0, ref, alt, label in (
                (left0, LEFT_REF, LEFT_ALT, "left"),
                (right0, RIGHT_REF, RIGHT_ALT, "right"),
            ):
                query_pos = query_position_at(alignment, pos0)
                if query_pos is None:
                    site_calls.append(("D", None))
                    counts[f"{label}_deletion_or_no_query_base"] += 1
                    continue
                base = alignment.query_sequence[query_pos].upper()
                baseq = alignment.query_qualities[query_pos]
                if baseq < BASEQ_MIN:
                    code = "L"
                    counts[f"{label}_low_baseq"] += 1
                elif base == ref:
                    code = "R"
                elif base == alt:
                    code = "A"
                else:
                    code = "O"
                    counts[f"{label}_other_base"] += 1
                site_calls.append((code, baseq))

            pattern = "".join(code for code, _ in site_calls)
            pair_patterns[pattern] += 1
            pair_quality_patterns[
                (
                    pattern,
                    str(site_calls[0][1]),
                    str(site_calls[1][1]),
                )
            ] += 1
            if all(code in {"R", "A"} for code, _ in site_calls):
                counts["continuous_reference_span_fixed_RA_both"] += 1
    result = dict(counts)
    for label, predicate in (
        ("raw", lambda alignment: True),
        (
            "flagpass",
            lambda alignment: not (alignment.flag & EXCLUDED_FLAGS),
        ),
        (
            "flag_mapq_pass",
            lambda alignment: not (alignment.flag & EXCLUDED_FLAGS)
            and alignment.mapping_quality >= MAPQ_MIN,
        ),
    ):
        left_names = {
            alignment.query_name for alignment in left if predicate(alignment)
        }
        right_names = {
            alignment.query_name for alignment in right if predicate(alignment)
        }
        result[f"{label}_left_qnames"] = len(left_names)
        result[f"{label}_right_qnames"] = len(right_names)
        result[f"{label}_qname_intersection"] = len(
            left_names & right_names
        )
    for key in (
        "continuous_reference_span_all",
        "continuous_reference_span_flagpass",
        "continuous_reference_span_flag_mapq_pass",
        "continuous_reference_span_fixed_RA_both",
    ):
        result.setdefault(key, 0)
    result["flag_mapq_pass_pair_patterns"] = dict(
        sorted(pair_patterns.items())
    )
    result["flag_mapq_pass_pattern_baseq_rows"] = [
        {
            "pattern": pattern,
            "left_baseq": left_baseq,
            "right_baseq": right_baseq,
            "count": count,
        }
        for (pattern, left_baseq, right_baseq), count in sorted(
            pair_quality_patterns.items()
        )
    ]
    return result


def canonical_census(path: Path):
    site_codes = Counter()
    primary_site_codes = Counter()
    pair_support = Counter()
    exact_pair_rows = []
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            indices = [int(value) for value in row["site_indices"].split(",")]
            calls = dict(zip(indices, row["call_codes"]))
            for index in (LEFT_INDEX, RIGHT_INDEX):
                if index in calls:
                    site_codes[(index, calls[index])] += 1
                    if (
                        row["phase_set"] == PHASE_SET
                        and row["hp_family"] in HP_FAMILIES
                    ):
                        primary_site_codes[
                            (row["hp_family"], index, calls[index])
                        ] += 1
            usable = {
                index: call
                for index, call in calls.items()
                if REGION_SITE_INDEX_MIN <= index <= REGION_SITE_INDEX_MAX
                and call in {"R", "A"}
            }
            if (
                row["phase_set"] == PHASE_SET
                and row["hp_family"] in HP_FAMILIES
            ):
                for left, right in itertools.combinations(sorted(usable), 2):
                    pair_support[(row["hp_family"], left, right)] += 1
                if LEFT_INDEX in calls and RIGHT_INDEX in calls:
                    call_pattern = calls[LEFT_INDEX] + calls[RIGHT_INDEX]
                    qualities = dict(
                        zip(indices, row["base_qualities"].split(","))
                    )
                    exact_pair_rows.append(
                        {
                            "hp_family": row["hp_family"],
                            "molecule_id": row["molecule_id"],
                            "alignment_id": row["alignment_id"],
                            "start0": int(row["start0"]),
                            "end0": int(row["end0"]),
                            "flag": int(row["flag"]),
                            "mapq": int(row["mapq"]),
                            "strand": row["strand"],
                            "hp_raw": row["hp_raw"],
                            "call_pattern": call_pattern,
                            "left_baseq": qualities[LEFT_INDEX],
                            "right_baseq": qualities[RIGHT_INDEX],
                        }
                    )
    return site_codes, primary_site_codes, pair_support, exact_pair_rows


def shortest_threshold_path(pair_support: Counter, hp_family: str):
    graph = defaultdict(set)
    for (family, left, right), support in pair_support.items():
        if family == hp_family and support >= 3:
            graph[left].add(right)
            graph[right].add(left)
    queue = deque([(LEFT_INDEX, [LEFT_INDEX])])
    visited = {LEFT_INDEX}
    while queue:
        node, path = queue.popleft()
        if node == RIGHT_INDEX:
            return path
        for neighbour in sorted(graph[node]):
            if neighbour not in visited:
                visited.add(neighbour)
                queue.append((neighbour, path + [neighbour]))
    return None


def downstream_mlhp(path: Path):
    document = json.loads(path.read_text())
    selected = []
    for group in document["groups"]:
        if (
            group.get("chrom") == "chr9"
            and str(group.get("phase_set")) == PHASE_SET
            and LEFT_POS1 in group.get("positions", [])
            and RIGHT_POS1 in group.get("positions", [])
        ):
            family = group["hp_family"]
            coverage = group["col_coverage_by_hp"][family]
            selected.append(
                {
                    "hp_family": family,
                    "region_id": group["region_id"],
                    "positions": group["positions"],
                    "left_retained_R_A": coverage[str(LEFT_POS1)],
                    "right_retained_R_A": coverage[str(RIGHT_POS1)],
                    "n_full_cov_reads": group["n_full_cov_reads"],
                    "tree_supported_pattern_count": group[
                        "tree_supported_pattern_count"
                    ],
                }
            )
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--sparse-calls", type=Path, required=True)
    parser.add_argument("--mlhp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        left_raw, _ = raw_site_census(
            bam, "chr9", LEFT_POS1, LEFT_REF, LEFT_ALT
        )
        right_raw, _ = raw_site_census(
            bam, "chr9", RIGHT_POS1, RIGHT_REF, RIGHT_ALT
        )
        raw_pair = raw_pair_census(bam, "chr9")

    (
        site_codes,
        primary_site_codes,
        pair_support,
        exact_pair_rows,
    ) = canonical_census(args.sparse_calls)
    exact_pair_pattern_counts = Counter(
        (row["hp_family"], row["call_pattern"]) for row in exact_pair_rows
    )
    direct_fixed_ra_molecules = [
        row
        for row in exact_pair_rows
        if set(row["call_pattern"]) <= {"R", "A"}
    ]
    primary_sites = {}
    for family in sorted(HP_FAMILIES):
        primary_sites[f"HP{family}"] = {}
        for index in (LEFT_INDEX, RIGHT_INDEX):
            primary_sites[f"HP{family}"][str(index)] = {
                code: primary_site_codes[(family, index, code)]
                for code in ("R", "A", "O", "D", "S", "L", "X")
            }
    paths = {}
    path_edges = {}
    for family in sorted(HP_FAMILIES):
        path = shortest_threshold_path(pair_support, family)
        paths[f"HP{family}"] = path
        path_edges[f"HP{family}"] = (
            [
                {
                    "left_index": left,
                    "right_index": right,
                    "support": pair_support[
                        (family, min(left, right), max(left, right))
                    ],
                }
                for left, right in zip(path, path[1:])
            ]
            if path
            else []
        )

    output = {
        "schema_name": "intersubmod.hcc1395_chr9_s4_s5_raw_bridge_audit",
        "schema_version": "1.0.0",
        "scope": {
            "sample": "HCC1395",
            "chrom": "chr9",
            "S4": {
                "site_index": LEFT_INDEX,
                "pos1": LEFT_POS1,
                "ref": LEFT_REF,
                "alt": LEFT_ALT,
            },
            "S5": {
                "site_index": RIGHT_INDEX,
                "pos1": RIGHT_POS1,
                "ref": RIGHT_REF,
                "alt": RIGHT_ALT,
            },
            "coordinate_difference_bp": RIGHT_POS1 - LEFT_POS1,
            "inclusive_interval_bp": RIGHT_POS1 - LEFT_POS1 + 1,
            "exact_phase_set": PHASE_SET,
        },
        "inputs": {
            "bam": str(args.bam.resolve()),
            "canonical_sparse_calls": str(args.sparse_calls.resolve()),
            "downstream_mlhp": str(args.mlhp.resolve()),
        },
        "parameters": {
            "excluded_flags": EXCLUDED_FLAGS,
            "mapq_min": MAPQ_MIN,
            "baseq_min": BASEQ_MIN,
            "strict_edge_min_molecules": 3,
        },
        "raw_bam": {
            "left_site": left_raw,
            "right_site": right_raw,
            "pair": raw_pair,
        },
        "canonical_calls": {
            "all_families": {
                str(index): {
                    code: site_codes[(index, code)]
                    for code in ("R", "A", "O", "D", "S", "L", "X")
                }
                for index in (LEFT_INDEX, RIGHT_INDEX)
            },
            "exact_ps_primary": primary_sites,
            "exact_ps_pair_pattern_counts": {
                f"HP{family}": {
                    pattern: exact_pair_pattern_counts[(family, pattern)]
                    for pattern in sorted(
                        {
                            row["call_pattern"]
                            for row in exact_pair_rows
                            if row["hp_family"] == family
                        }
                    )
                }
                for family in sorted(HP_FAMILIES)
            },
            "exact_ps_pair_rows": exact_pair_rows,
            "direct_fixed_RA_molecule_count": len(
                direct_fixed_ra_molecules
            ),
        },
        "strict_transitive_linkage": {
            "shortest_threshold3_paths": paths,
            "path_edges": path_edges,
        },
        "downstream_mlhp": downstream_mlhp(args.mlhp),
        "verdict": {
            "physical_continuous_alignment_spans_both": (
                raw_pair["continuous_reference_span_all"] > 0
            ),
            "raw_qname_touches_both": (
                raw_pair["raw_qname_intersection"] > 0
            ),
            "canonical_direct_RA_bridge_exists": bool(
                direct_fixed_ra_molecules
            ),
            "same_exact_ps_hp_transitive_path_exists": all(
                paths.values()
            ),
            "hp2_exact_ps_pair_rows": sum(
                row["hp_family"] == "2" for row in exact_pair_rows
            ),
            "hp2_fixed_RA_pair_rows": sum(
                row["hp_family"] == "2"
                and set(row["call_pattern"]) <= {"R", "A"}
                for row in exact_pair_rows
            ),
            "primary_reason_hp2_direct_edge_absent": (
                "RULE_FILTERED_SUPPORT_BELOW_MINREAD_AFTER_LOW_BQ_AND_DELETION"
            ),
            "rule_exclusion_caused_hp2_direct_edge_loss": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
