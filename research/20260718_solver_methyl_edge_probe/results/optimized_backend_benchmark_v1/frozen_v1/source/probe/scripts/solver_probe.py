#!/usr/bin/env python3
"""Isolated exact-solver prototypes for the 2026-07-18 bounded probe.

This module intentionally does not modify or import into the canonical pipeline.
It implements the same minimum-extra-vertex contract in two independent ways:

1. persistent HiGHS with incremental cardinality/no-good rows;
2. exact obligation-driven branch-and-bound.

The branch-and-bound keeps both selected and excluded sets in its memo key.
Every partial group and every selected orphan vertex is an obligation that must
be hit by at least one available vertex.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import pathlib
import sys
import time
from dataclasses import asdict, dataclass
from itertools import combinations
from typing import Any, Iterable, Sequence


def popcount(value: int) -> int:
    return bin(int(value)).count("1")


def parse_pattern(pattern: str, *, allow_x: bool) -> tuple[int, int, int]:
    if not pattern:
        raise ValueError("pattern must not be empty")
    allowed = {"R", "A", "X"} if allow_x else {"R", "A"}
    invalid = set(pattern) - allowed
    if invalid:
        raise ValueError(f"invalid symbols in {pattern!r}: {sorted(invalid)}")
    fixed_mask = 0
    alt_mask = 0
    for bit, symbol in enumerate(pattern):
        if symbol in {"R", "A"}:
            fixed_mask |= 1 << bit
        if symbol == "A":
            alt_mask |= 1 << bit
    return len(pattern), fixed_mask, alt_mask


def submasks(mask: int) -> tuple[int, ...]:
    values: list[int] = []
    current = mask
    while True:
        values.append(current)
        if current == 0:
            break
        current = (current - 1) & mask
    return tuple(sorted(values))


def vertex_set_digest(k: int, vertices: Iterable[int]) -> str:
    payload = {"k": int(k), "vertices": sorted(int(value) for value in vertices)}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class ExactInstance:
    k: int
    universe_mask: int
    vertices: tuple[int, ...]
    mandatory: frozenset[int]
    groups: tuple[frozenset[int], ...]
    full_patterns: tuple[str, ...]
    partial_patterns: tuple[str, ...]
    universe_mode: str

    @property
    def extras(self) -> tuple[int, ...]:
        return tuple(vertex for vertex in self.vertices if vertex not in self.mandatory)

    def predecessors(self, vertex: int) -> tuple[int, ...]:
        return tuple(
            vertex ^ (1 << bit)
            for bit in range(self.k)
            if vertex & (1 << bit) and self.universe_mask & (1 << bit)
        )


def build_instance(
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
    k: int,
    *,
    universe_mode: str = "observed_alt",
) -> ExactInstance:
    if k < 1:
        raise ValueError("k must be positive")
    if universe_mode not in {"observed_alt", "all_loci"}:
        raise ValueError("universe_mode must be observed_alt or all_loci")
    full = tuple(full_patterns)
    partial = tuple(partial_patterns)
    for pattern in full:
        if len(pattern) != k:
            raise ValueError("full pattern length differs from k")
        parse_pattern(pattern, allow_x=False)
    for pattern in partial:
        if len(pattern) != k:
            raise ValueError("partial pattern length differs from k")
        parse_pattern(pattern, allow_x=True)

    universe_mask = (1 << k) - 1
    if universe_mode == "observed_alt":
        universe_mask = 0
        for pattern in full + partial:
            _, _, alt_mask = parse_pattern(pattern, allow_x=("X" in pattern))
            universe_mask |= alt_mask
    vertices = submasks(universe_mask)
    mandatory = {0}
    for pattern in full:
        _, _, alt_mask = parse_pattern(pattern, allow_x=False)
        mandatory.add(alt_mask)
    groups: list[frozenset[int]] = []
    for pattern in partial:
        _, fixed_mask, alt_mask = parse_pattern(pattern, allow_x=True)
        if alt_mask & ~universe_mask:
            raise ValueError("partial ALT lies outside solver universe")
        compatible = frozenset(
            vertex
            for vertex in vertices
            if ((vertex ^ alt_mask) & fixed_mask & universe_mask) == 0
        )
        if not compatible:
            raise ValueError(f"partial pattern has no compatible state: {pattern}")
        groups.append(compatible)
    return ExactInstance(
        k=int(k),
        universe_mask=int(universe_mask),
        vertices=vertices,
        mandatory=frozenset(mandatory),
        groups=tuple(groups),
        full_patterns=full,
        partial_patterns=partial,
        universe_mode=universe_mode,
    )


def check_selected(instance: ExactInstance, selected: Iterable[int]) -> dict[str, Any]:
    selected_set = frozenset(int(value) for value in selected)
    missing_mandatory = sorted(instance.mandatory - selected_set)
    outside = sorted(selected_set - set(instance.vertices))
    uncovered = [
        index for index, group in enumerate(instance.groups)
        if selected_set.isdisjoint(group)
    ]
    orphans = [
        vertex for vertex in sorted(selected_set - {0})
        if selected_set.isdisjoint(instance.predecessors(vertex))
    ]
    return {
        "pass": not missing_mandatory and not outside and not uncovered and not orphans,
        "missing_mandatory": missing_mandatory,
        "outside_universe": outside,
        "uncovered_group_indices": uncovered,
        "orphan_vertices": orphans,
    }


def brute_force_optimal(
    instance: ExactInstance,
    *,
    combination_cap: int = 2_000_000,
) -> dict[str, Any]:
    checked = 0
    for size in range(len(instance.extras) + 1):
        solutions: list[tuple[int, ...]] = []
        for extra in combinations(instance.extras, size):
            checked += 1
            if checked > combination_cap:
                return {
                    "complete": False,
                    "status": "COMBINATION_CAP",
                    "checked": checked,
                    "objective": None,
                    "vertex_sets": [],
                }
            selected = tuple(sorted(instance.mandatory | set(extra)))
            if check_selected(instance, selected)["pass"]:
                solutions.append(selected)
        if solutions:
            return {
                "complete": True,
                "status": "OPTIMAL_COMPLETE",
                "checked": checked,
                "objective": size,
                "vertex_sets": sorted(set(solutions)),
            }
    raise AssertionError("rooted instance unexpectedly infeasible")


@dataclass
class BnbStats:
    visited_states: int = 0
    memo_hits: int = 0
    propagated_singletons: int = 0
    pruned_infeasible: int = 0
    pruned_bound: int = 0
    feasible_leaves: int = 0
    max_depth: int = 0
    deadline_checks: int = 0


@dataclass
class BnbResult:
    status: str
    complete: bool
    objective: int | None
    vertex_sets: list[tuple[int, ...]]
    elapsed_seconds: float
    stats: BnbStats
    lower_bound_name: str

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["vertex_set_ids"] = []
        return result


class ExactObligationBnb:
    """Exact selected/excluded obligation search with safe lower bounds."""

    def __init__(
        self,
        instance: ExactInstance,
        *,
        time_limit_seconds: float = 30.0,
        max_sets: int | None = None,
    ) -> None:
        self.instance = instance
        self.deadline = time.perf_counter() + float(time_limit_seconds)
        self.max_sets = max_sets
        self.best = math.inf
        self.solutions: set[tuple[int, ...]] = set()
        self.memo: set[tuple[frozenset[int], frozenset[int]]] = set()
        self.stats = BnbStats()
        self.timed_out = False
        self.set_cap_reached = False

    def _obligations(
        self,
        selected: frozenset[int],
        excluded: frozenset[int],
    ) -> tuple[list[tuple[str, int, frozenset[int]]], bool]:
        obligations: list[tuple[str, int, frozenset[int]]] = []
        for index, group in enumerate(self.instance.groups):
            if selected.isdisjoint(group):
                domain = group - excluded
                if not domain:
                    return [], False
                obligations.append(("group", index, frozenset(domain)))
        for vertex in sorted(selected - {0}):
            predecessors = frozenset(self.instance.predecessors(vertex))
            if selected.isdisjoint(predecessors):
                domain = predecessors - excluded
                if not domain:
                    return [], False
                obligations.append(("parent", vertex, frozenset(domain)))
        return obligations, True

    def _root_reachable(
        self,
        vertex: int,
        excluded: frozenset[int],
        memo: dict[int, bool],
    ) -> bool:
        if vertex in memo:
            return memo[vertex]
        if vertex in excluded:
            memo[vertex] = False
        elif vertex == 0:
            memo[vertex] = True
        else:
            memo[vertex] = any(
                self._root_reachable(parent, excluded, memo)
                for parent in self.instance.predecessors(vertex)
            )
        return memo[vertex]

    def _propagate(
        self,
        selected: frozenset[int],
        excluded: frozenset[int],
    ) -> tuple[frozenset[int], frozenset[int], list[tuple[str, int, frozenset[int]]]] | None:
        if selected & excluded:
            return None
        current = set(selected)
        while True:
            reachable_memo: dict[int, bool] = {}
            if any(
                not self._root_reachable(vertex, excluded, reachable_memo)
                for vertex in current
            ):
                return None
            obligations, possible = self._obligations(frozenset(current), excluded)
            if not possible:
                return None
            singleton_values = {
                next(iter(domain))
                for _, _, domain in obligations
                if len(domain) == 1
            }
            if not singleton_values:
                return frozenset(current), excluded, obligations
            if singleton_values & excluded:
                return None
            before = len(current)
            current.update(singleton_values)
            self.stats.propagated_singletons += len(current) - before

    def _connection_cost(
        self,
        vertex: int,
        selected: frozenset[int],
        excluded: frozenset[int],
        memo: dict[int, float],
    ) -> float:
        if vertex in memo:
            return memo[vertex]
        if vertex in excluded:
            value = math.inf
        elif vertex == 0:
            value = 0.0 if vertex in selected else 1.0
        else:
            predecessor_cost = min(
                (
                    self._connection_cost(parent, selected, excluded, memo)
                    for parent in self.instance.predecessors(vertex)
                ),
                default=math.inf,
            )
            value = predecessor_cost + (0.0 if vertex in selected else 1.0)
        memo[vertex] = value
        return value

    def _lower_bound(
        self,
        selected: frozenset[int],
        excluded: frozenset[int],
        obligations: Sequence[tuple[str, int, frozenset[int]]],
    ) -> int:
        current_cost = len(selected - self.instance.mandatory)
        if not obligations:
            return current_cost

        domains = [domain for _, _, domain in obligations]
        used: set[int] = set()
        packing = 0
        for domain in sorted(domains, key=lambda values: (len(values), tuple(sorted(values)))):
            if used.isdisjoint(domain):
                used.update(domain)
                packing += 1

        coverage: dict[int, int] = {}
        for domain in domains:
            for vertex in domain:
                coverage[vertex] = coverage.get(vertex, 0) + 1
        maximum_coverage = max(coverage.values(), default=1)
        cover_bound = math.ceil(len(domains) / maximum_coverage)

        connection_memo: dict[int, float] = {}
        selected_connection = max(
            (
                self._connection_cost(vertex, selected, excluded, connection_memo)
                for vertex in selected
            ),
            default=0.0,
        )
        group_connection = max(
            (
                min(
                    self._connection_cost(vertex, selected, excluded, connection_memo)
                    for vertex in domain
                )
                for kind, _, domain in obligations
                if kind == "group"
            ),
            default=0.0,
        )
        additional = max(packing, cover_bound, selected_connection, group_connection)
        if not math.isfinite(additional):
            return math.inf
        return current_cost + int(math.ceil(additional - 1e-12))

    def _candidate_order(
        self,
        domain: frozenset[int],
        selected: frozenset[int],
        excluded: frozenset[int],
        obligations: Sequence[tuple[str, int, frozenset[int]]],
    ) -> list[int]:
        coverage = {
            vertex: sum(vertex in other for _, _, other in obligations)
            for vertex in domain
        }
        connection_memo: dict[int, float] = {}
        return sorted(
            domain,
            key=lambda vertex: (
                self._connection_cost(vertex, selected, excluded, connection_memo),
                -coverage[vertex],
                popcount(vertex),
                vertex,
            ),
        )

    def _search(
        self,
        selected: frozenset[int],
        excluded: frozenset[int],
        depth: int,
    ) -> None:
        self.stats.deadline_checks += 1
        if time.perf_counter() >= self.deadline:
            self.timed_out = True
            return
        if self.max_sets is not None and len(self.solutions) >= self.max_sets:
            self.set_cap_reached = True
            return
        propagated = self._propagate(selected, excluded)
        if propagated is None:
            self.stats.pruned_infeasible += 1
            return
        selected, excluded, obligations = propagated
        key = (selected, excluded)
        if key in self.memo:
            self.stats.memo_hits += 1
            return
        self.memo.add(key)
        self.stats.visited_states += 1
        self.stats.max_depth = max(self.stats.max_depth, depth)

        lower_bound = self._lower_bound(selected, excluded, obligations)
        if lower_bound > self.best:
            self.stats.pruned_bound += 1
            return
        if not obligations:
            self.stats.feasible_leaves += 1
            cost = len(selected - self.instance.mandatory)
            values = tuple(sorted(selected))
            if not check_selected(self.instance, values)["pass"]:
                raise AssertionError("B&B emitted an invalid vertex set")
            if cost < self.best:
                self.best = cost
                self.solutions = {values}
            elif cost == self.best:
                self.solutions.add(values)
            return

        chosen = min(
            obligations,
            key=lambda item: (
                len(item[2]),
                0 if item[0] == "group" else 1,
                item[1],
            ),
        )
        ordered = self._candidate_order(chosen[2], selected, excluded, obligations)
        earlier: set[int] = set()
        for candidate in ordered:
            if self.timed_out or self.set_cap_reached:
                return
            self._search(
                frozenset(set(selected) | {candidate}),
                frozenset(set(excluded) | earlier),
                depth + 1,
            )
            earlier.add(candidate)

    def run(self) -> BnbResult:
        started = time.perf_counter()
        self._search(self.instance.mandatory, frozenset(), 0)
        elapsed = time.perf_counter() - started
        complete = not self.timed_out and not self.set_cap_reached
        if self.timed_out:
            status = "TIME_LIMIT_INCOMPLETE"
        elif self.set_cap_reached:
            status = "MAX_SETS_INCOMPLETE"
        elif math.isfinite(self.best):
            status = "OPTIMAL_COMPLETE"
        else:
            status = "INFEASIBLE"
        return BnbResult(
            status=status,
            complete=complete,
            objective=(None if not math.isfinite(self.best) else int(self.best)),
            vertex_sets=sorted(self.solutions),
            elapsed_seconds=elapsed,
            stats=self.stats,
            lower_bound_name="max(disjoint_packing,coverage,connection)",
        )


def load_highspy(highspy_path: str | None):
    if highspy_path:
        path = str(pathlib.Path(highspy_path).resolve())
        if path not in sys.path:
            sys.path.insert(0, path)
    try:
        import highspy  # type: ignore
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "highspy is required; pass --highspy-path pointing at an isolated install"
        ) from error
    return highspy


def persistent_highs_enumerate(
    instance: ExactInstance,
    *,
    highspy_path: str | None,
    time_limit_seconds_per_run: float = 30.0,
    max_sets: int = 256,
) -> dict[str, Any]:
    highspy = load_highspy(highspy_path)
    started = time.perf_counter()
    highs = highspy.Highs()
    for name, value in (
        ("output_flag", False),
        ("threads", 1),
        ("random_seed", 20260718),
        ("mip_rel_gap", 0.0),
        ("time_limit", float(time_limit_seconds_per_run)),
    ):
        status = highs.setOptionValue(name, value)
        if status != highspy.HighsStatus.kOk:
            raise RuntimeError(f"failed to set HiGHS option {name}={value}")

    index = {vertex: idx for idx, vertex in enumerate(instance.vertices)}
    lower = [1.0 if vertex in instance.mandatory else 0.0 for vertex in instance.vertices]
    upper = [1.0 for _ in instance.vertices]
    objective = [0.0 if vertex in instance.mandatory else 1.0 for vertex in instance.vertices]
    variables = highs.addVariables(
        len(instance.vertices),
        lb=lower,
        ub=upper,
        obj=objective,
        type=highspy.HighsVarType.kInteger,
        out_array=True,
    )
    structural_rows = 0
    for vertex in instance.vertices:
        if vertex == 0:
            continue
        expression = variables[index[vertex]]
        for parent in instance.predecessors(vertex):
            expression = expression - variables[index[parent]]
        highs.addConstr(expression <= 0.0)
        structural_rows += 1
    for group in instance.groups:
        expression = sum((variables[index[vertex]] for vertex in group), 0)
        highs.addConstr(expression >= 1.0)
        structural_rows += 1
    model_builds = 1
    solve_calls = 0
    no_good_rows = 0
    cardinality_rows = 0
    sets: list[tuple[int, ...]] = []
    run_statuses: list[str] = []

    highs.run()
    solve_calls += 1
    model_status = highs.getModelStatus()
    run_statuses.append(highs.modelStatusToString(model_status))
    if model_status != highspy.HighsModelStatus.kOptimal:
        return {
            "status": f"PHASE_A_{highs.modelStatusToString(model_status).upper()}",
            "complete": False,
            "objective": None,
            "vertex_sets": [],
            "model_builds": model_builds,
            "solve_calls": solve_calls,
            "structural_rows": structural_rows,
            "cardinality_rows": cardinality_rows,
            "no_good_rows": no_good_rows,
            "elapsed_seconds": time.perf_counter() - started,
            "highspy_version": highs.version(),
            "run_statuses": run_statuses,
        }
    optimum = int(round(highs.getObjectiveValue()))
    if optimum == 0:
        values = highs.allVariableValues()
        selected = tuple(
            vertex
            for vertex, value in zip(instance.vertices, values)
            if float(value) > 0.5
        )
        checker = check_selected(instance, selected)
        if not checker["pass"] or set(selected) != set(instance.mandatory):
            raise AssertionError(f"invalid zero-extra optimum: {checker}")
        return {
            "status": "OPTIMAL_COMPLETE",
            "complete": True,
            "objective": optimum,
            "vertex_sets": [selected],
            "model_builds": model_builds,
            "solve_calls": solve_calls,
            "structural_rows": structural_rows,
            "cardinality_rows": cardinality_rows,
            "no_good_rows": no_good_rows,
            "elapsed_seconds": time.perf_counter() - started,
            "highspy_version": highs.version(),
            "run_statuses": run_statuses,
        }
    extra_expression = sum(
        (variables[index[vertex]] for vertex in instance.extras),
        0,
    )
    highs.addConstr(extra_expression == float(optimum))
    cardinality_rows = 1

    while True:
        values = highs.allVariableValues()
        selected = tuple(
            vertex
            for vertex, value in zip(instance.vertices, values)
            if float(value) > 0.5
        )
        checker = check_selected(instance, selected)
        if not checker["pass"]:
            raise AssertionError(f"persistent HiGHS returned invalid set: {checker}")
        if len(selected) - len(instance.mandatory) != optimum:
            raise AssertionError("persistent HiGHS returned wrong fixed cardinality")
        if selected in sets:
            raise AssertionError("incremental no-good failed to exclude a prior set")
        sets.append(selected)
        if len(sets) >= max_sets:
            return {
                "status": "MAX_SETS_INCOMPLETE",
                "complete": False,
                "objective": optimum,
                "vertex_sets": sorted(sets),
                "model_builds": model_builds,
                "solve_calls": solve_calls,
                "structural_rows": structural_rows,
                "cardinality_rows": cardinality_rows,
                "no_good_rows": no_good_rows,
                "elapsed_seconds": time.perf_counter() - started,
                "highspy_version": highs.version(),
                "run_statuses": run_statuses,
            }
        selected_extra = [vertex for vertex in selected if vertex not in instance.mandatory]
        no_good_expression = sum(
            (variables[index[vertex]] for vertex in selected_extra),
            0,
        )
        highs.addConstr(no_good_expression <= float(optimum - 1))
        no_good_rows += 1
        highs.run()
        solve_calls += 1
        model_status = highs.getModelStatus()
        run_statuses.append(highs.modelStatusToString(model_status))
        if model_status == highspy.HighsModelStatus.kInfeasible:
            return {
                "status": "OPTIMAL_COMPLETE",
                "complete": True,
                "objective": optimum,
                "vertex_sets": sorted(sets),
                "model_builds": model_builds,
                "solve_calls": solve_calls,
                "structural_rows": structural_rows,
                "cardinality_rows": cardinality_rows,
                "no_good_rows": no_good_rows,
                "elapsed_seconds": time.perf_counter() - started,
                "highspy_version": highs.version(),
                "run_statuses": run_statuses,
            }
        if model_status != highspy.HighsModelStatus.kOptimal:
            return {
                "status": f"ENUMERATION_{highs.modelStatusToString(model_status).upper()}",
                "complete": False,
                "objective": optimum,
                "vertex_sets": sorted(sets),
                "model_builds": model_builds,
                "solve_calls": solve_calls,
                "structural_rows": structural_rows,
                "cardinality_rows": cardinality_rows,
                "no_good_rows": no_good_rows,
                "elapsed_seconds": time.perf_counter() - started,
                "highspy_version": highs.version(),
                "run_statuses": run_statuses,
            }


def load_module(path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def current_scipy_enumerate(
    solver_path: str,
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
    k: int,
    *,
    max_sets: int,
    time_limit_seconds: float,
) -> dict[str, Any]:
    module = load_module(solver_path, "current_hypercube_exact_probe")
    original = module._build_problem
    counter = {"calls": 0}

    def counted(*args, **kwargs):
        counter["calls"] += 1
        return original(*args, **kwargs)

    module._build_problem = counted
    started = time.perf_counter()
    result = module.enumerate_optimal_vertex_sets(
        list(full_patterns),
        list(partial_patterns),
        int(k),
        universe_mode="observed_alt",
        max_sets=int(max_sets),
        time_limit_seconds=float(time_limit_seconds),
    )
    return {
        "status": (
            "OPTIMAL_COMPLETE"
            if result.get("complete")
            else str(result.get("stop_status") or result.get("first", {}).get("status"))
        ),
        "complete": bool(result.get("complete")),
        "objective": result.get("objective", result.get("first", {}).get("objective")),
        "vertex_sets": [tuple(values) for values in result.get("vertex_sets", [])],
        "model_builds": counter["calls"],
        "solve_calls": counter["calls"],
        "elapsed_seconds": time.perf_counter() - started,
    }


def summarize_result(instance: ExactInstance, result: dict[str, Any]) -> dict[str, Any]:
    vertex_sets = [tuple(sorted(values)) for values in result.get("vertex_sets", [])]
    summarized = dict(result)
    summarized["vertex_sets"] = [list(values) for values in vertex_sets]
    summarized["vertex_set_ids"] = [
        vertex_set_digest(instance.k, values) for values in vertex_sets
    ]
    summarized["n_vertex_sets"] = len(vertex_sets)
    summarized["all_sets_independently_valid"] = all(
        check_selected(instance, values)["pass"] for values in vertex_sets
    )
    return summarized


def read_real_fixtures(path: str) -> list[dict[str, Any]]:
    with open(path, encoding="utf-8") as handle:
        document = json.load(handle)
    return list(document["units"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver-path", required=True)
    parser.add_argument("--highspy-path", required=True)
    parser.add_argument("--real-fixtures", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--time-limit", type=float, default=30.0)
    args = parser.parse_args()

    cases = [
        {
            "case_id": "root_only_k2",
            "k": 2,
            "full": [],
            "partial": [],
            "expected_h": 0,
            "expected_sets": 1,
        },
        {
            "case_id": "full_AA",
            "k": 2,
            "full": ["AA"],
            "partial": [],
            "expected_h": 1,
            "expected_sets": 2,
        },
        {
            "case_id": "partial_AX_XA",
            "k": 2,
            "full": [],
            "partial": ["AX", "XA"],
            "expected_h": 2,
            "expected_sets": 3,
        },
        {
            "case_id": "full_AAAA_24",
            "k": 4,
            "full": ["AAAA"],
            "partial": [],
            "expected_h": 3,
            "expected_sets": 24,
        },
        {
            "case_id": "full_RRA_AAA_partial_RAX",
            "k": 3,
            "full": ["RRA", "AAA"],
            "partial": ["RAX"],
            "expected_h": 1,
            "expected_sets": 1,
        },
        {
            "case_id": "recurrence_allowed_AR_RA_AA",
            "k": 2,
            "full": ["AR", "RA", "AA"],
            "partial": [],
            "expected_h": 0,
            "expected_sets": 1,
        },
    ]
    synthetic: list[dict[str, Any]] = []
    for case in cases:
        instance = build_instance(case["full"], case["partial"], case["k"])
        brute = brute_force_optimal(instance)
        bnb_obj = ExactObligationBnb(instance, time_limit_seconds=args.time_limit).run()
        bnb = summarize_result(instance, bnb_obj.to_dict())
        scipy_result = current_scipy_enumerate(
            args.solver_path,
            case["full"],
            case["partial"],
            case["k"],
            max_sets=512,
            time_limit_seconds=args.time_limit,
        )
        scipy_summary = summarize_result(instance, scipy_result)
        persistent = summarize_result(
            instance,
            persistent_highs_enumerate(
                instance,
                highspy_path=args.highspy_path,
                time_limit_seconds_per_run=args.time_limit,
                max_sets=512,
            ),
        )
        expected_digest = sorted(
            vertex_set_digest(instance.k, values) for values in brute["vertex_sets"]
        )
        row = {
            "case_id": case["case_id"],
            "k": case["k"],
            "effective_m": popcount(instance.universe_mask),
            "n_vertices": len(instance.vertices),
            "expected_h": case["expected_h"],
            "expected_sets": case["expected_sets"],
            "brute_force": summarize_result(instance, brute),
            "bnb": bnb,
            "current_scipy": scipy_summary,
            "persistent_highs": persistent,
        }
        row["checks"] = {
            "brute_matches_expected": (
                brute["objective"] == case["expected_h"]
                and len(brute["vertex_sets"]) == case["expected_sets"]
            ),
            "bnb_complete_digest_match": (
                bnb["complete"] and sorted(bnb["vertex_set_ids"]) == expected_digest
            ),
            "scipy_complete_digest_match": (
                scipy_summary["complete"]
                and sorted(scipy_summary["vertex_set_ids"]) == expected_digest
            ),
            "persistent_complete_digest_match": (
                persistent["complete"]
                and sorted(persistent["vertex_set_ids"]) == expected_digest
            ),
        }
        synthetic.append(row)

    real_rows: list[dict[str, Any]] = []
    for fixture in read_real_fixtures(args.real_fixtures):
        instance = build_instance(
            list(fixture["full_counts"]),
            list(fixture["partial_counts"]),
            int(fixture["k"]),
        )
        bnb_obj = ExactObligationBnb(
            instance,
            time_limit_seconds=min(args.time_limit, 15.0),
            max_sets=int(fixture.get("max_sets", 256)),
        ).run()
        bnb = summarize_result(instance, bnb_obj.to_dict())
        persistent = summarize_result(
            instance,
            persistent_highs_enumerate(
                instance,
                highspy_path=args.highspy_path,
                time_limit_seconds_per_run=min(args.time_limit, 5.0),
                max_sets=int(fixture.get("max_sets", 256)),
            ),
        )
        real_rows.append(
            {
                "case_id": fixture["case_id"],
                "sample": fixture["sample"],
                "region": fixture["region"],
                "family": fixture["family"],
                "k": fixture["k"],
                "effective_m": popcount(instance.universe_mask),
                "n_vertices": len(instance.vertices),
                "expected_h": fixture["expected_h"],
                "bnb": bnb,
                "persistent_highs": persistent,
                "checks": {
                    "bnb_objective_if_available": (
                        bnb["objective"] is None
                        or bnb["objective"] == fixture["expected_h"]
                    ),
                    "persistent_objective_if_available": (
                        persistent["objective"] is None
                        or persistent["objective"] == fixture["expected_h"]
                    ),
                    "complete_digest_match_if_both_complete": (
                        not (bnb["complete"] and persistent["complete"])
                        or sorted(bnb["vertex_set_ids"]) == sorted(persistent["vertex_set_ids"])
                    ),
                },
            }
        )

    all_synthetic_pass = all(
        all(row["checks"].values()) for row in synthetic
    )
    receipt = {
        "schema": "intersubmod.solver_methyl_edge_probe.solver.v1",
        "scope": "PARTIAL_EXPLORATORY",
        "inputs": {
            "solver_path": str(pathlib.Path(args.solver_path).resolve()),
            "highspy_path": str(pathlib.Path(args.highspy_path).resolve()),
            "real_fixtures": str(pathlib.Path(args.real_fixtures).resolve()),
        },
        "synthetic": synthetic,
        "real": real_rows,
        "summary": {
            "n_synthetic": len(synthetic),
            "all_synthetic_pass": all_synthetic_pass,
            "persistent_highs_version": synthetic[0]["persistent_highs"]["highspy_version"],
            "aaaa_current_model_builds": next(
                row["current_scipy"]["model_builds"]
                for row in synthetic if row["case_id"] == "full_AAAA_24"
            ),
            "aaaa_persistent_model_builds": next(
                row["persistent_highs"]["model_builds"]
                for row in synthetic if row["case_id"] == "full_AAAA_24"
            ),
        },
    }
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt["summary"], ensure_ascii=False, sort_keys=True))
    return 0 if all_synthetic_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
