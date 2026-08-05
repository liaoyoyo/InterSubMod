#!/usr/bin/env python3
"""Instrument the current ranker to verify one primary fit per vertex set."""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
from typing import Any


def load_module(path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def make_unit(module, pattern_counts: list[tuple[str, int]], k: int):
    component = module.Component(
        "PROBE",
        "chr1",
        3,
        "hp1",
        "COMP",
        10,
        10 + k - 1,
        k,
        tuple(range(k)),
    )
    unit = module.UnitEvidence(component, "1")
    for pattern, count in pattern_counts:
        codes = {
            index: symbol
            for index, symbol in enumerate(pattern)
            if symbol in {"R", "A"}
        }
        qualities = {index: 30 for index in codes}
        for _ in range(count):
            unit.add_projection(codes, qualities)
    return unit


def instrument_case(module, case_id: str, patterns: list[tuple[str, int]], k: int) -> dict[str, Any]:
    original = module.fit_quality_aware_mixture
    calls: list[tuple[int, ...]] = []

    def counted(pattern_counts, vertices, **kwargs):
        calls.append(tuple(sorted(int(value) for value in vertices)))
        return original(pattern_counts, vertices, **kwargs)

    module.fit_quality_aware_mixture = counted
    try:
        result = module.rank_unit(
            make_unit(module, patterns, k),
            fixed_error_grid=(),
            conditional_candidate_ranking_bootstrap_replicates=0,
            solver_time_limit_seconds=5,
        )
    finally:
        module.fit_quality_aware_mixture = original
    distinct_calls = sorted(set(calls))
    return {
        "case_id": case_id,
        "raw_tree_candidates_T": result["raw_tree_candidates_T"],
        "distinct_vertex_sets_V": result["distinct_vertex_sets_V"],
        "top_edge_variants": result["top_edge_variants"],
        "selection_status": result["selection_status"],
        "primary_fit_calls": len(calls),
        "distinct_fit_vertex_sets": len(distinct_calls),
        "fit_vertices": [list(values) for values in calls],
        "check_one_primary_fit_per_vertex_set": (
            len(calls) == result["distinct_vertex_sets_V"] == len(distinct_calls)
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranker-path", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    ranker_path = pathlib.Path(args.ranker_path).resolve()
    current_scripts = str(ranker_path.parent)
    if current_scripts not in sys.path:
        sys.path.insert(0, current_scripts)
    module = load_module(str(ranker_path), "current_build_m2_patterns_and_rank_probe")
    cases = [
        instrument_case(
            module,
            "same_V_two_parent_assignments",
            [("AR", 3), ("RA", 3), ("AA", 3)],
            2,
        ),
        instrument_case(
            module,
            "two_distinct_V",
            [("AA", 3)],
            2,
        ),
    ]
    receipt = {
        "schema": "intersubmod.solver_methyl_edge_probe.read_af_cache.v1",
        "scope": "PRIMARY_QUALITY_AWARE_FIT_ONLY",
        "ranker_path": str(ranker_path),
        "cases": cases,
        "all_pass": all(
            case["check_one_primary_fit_per_vertex_set"] for case in cases
        ),
        "interpretation": (
            "Fixed-error sensitivity grids and bootstrap replicates intentionally "
            "refit each vertex set under a different model or resampled data; they "
            "are not duplicate parent-assignment fits."
        ),
    }
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_pass": receipt["all_pass"], "cases": cases}, ensure_ascii=False))
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
