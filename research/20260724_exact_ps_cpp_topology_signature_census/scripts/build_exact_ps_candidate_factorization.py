#!/usr/bin/env python3
"""Build and receipt the all-seven exact-PS candidate-factorization sidecar."""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


HERE = Path(__file__).resolve().parent
TOPIC_ROOT = HERE.parent
REPO_ROOT = HERE.parents[2]
CPP_SOURCE = TOPIC_ROOT / "cpp" / "exact_topology_candidate_factorization.cpp"
FROZEN_SOURCE = (
    REPO_ROOT
    / "research"
    / "20260724_exact_ps_cpp_topology_af_all_samples"
    / "cpp"
    / "exact_ps_topology_af.cpp"
)
LONG_LINEAGE = Path("/big7_disk/liaoyoyo2001/LongLineage")
OBLIGATION_SOURCE = LONG_LINEAGE / "src/solver/obligation_bnb.cpp"
PARENT_MAPPING_SOURCE = LONG_LINEAGE / "src/solver/parent_mapping.cpp"
OBLIGATION_HEADER = LONG_LINEAGE / "include/longlineage/solver/obligation_bnb.hpp"
PARENT_MAPPING_HEADER = LONG_LINEAGE / "include/longlineage/solver/parent_mapping.hpp"
TOPOLOGY_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1"
)
CENSUS_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1"
)
TOPOLOGY_RECEIPT = TOPOLOGY_ROOT / "cohort_receipt.json"
TOPOLOGY_SUMMARY = TOPOLOGY_ROOT / "summary" / "all7_summary.json"
CENSUS_RECEIPT = CENSUS_ROOT / "receipt.v2.json"
EXPECTED_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
EXPECTED_RANKED = 71_955
EXPECTED_MINIMUM_TREES = 972_592
EXPECTED_GLOBAL_BEST_TREES = 680_527
EXPECTED_CHROMOSOMES = [f"chr{value}" for value in range(1, 23)]
EXPECTED_SAMPLE_COUNTS = {
    "HCC1395": (9_130, 48_337, 48_337, 15_350),
    "HCC1395_DORADO": (5_308, 14_919, 14_919, 8_292),
    "COLO829": (10_757, 53_345, 53_345, 47_025),
    "H1437": (13_740, 190_549, 190_549, 123_734),
    "H2009": (23_128, 630_580, 630_580, 471_910),
    "HCC1937": (4_245, 15_527, 15_527, 6_674),
    "HCC1954": (5_647, 19_335, 19_335, 7_542),
}
IMPLEMENTATION_PATHS = {
    "candidate_factorization_source": CPP_SOURCE,
    "frozen_topology_source": FROZEN_SOURCE,
    "longlineage_obligation_bnb_source": OBLIGATION_SOURCE,
    "longlineage_parent_mapping_source": PARENT_MAPPING_SOURCE,
    "longlineage_obligation_bnb_header": OBLIGATION_HEADER,
    "longlineage_parent_mapping_header": PARENT_MAPPING_HEADER,
}
UNIT_SCHEMA = (
    "intersubmod.exact_ps_topology_candidate_factorization.unit",
    "1.0.0",
)
RECEIPT_SCHEMA = (
    "intersubmod.exact_ps_topology_candidate_factorization.receipt",
    "1.1.0",
)


class FactorizationBuildError(RuntimeError):
    """Raised when a sidecar cannot be produced or validated exactly."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FactorizationBuildError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "sha256": sha256_path(resolved),
    }


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise FactorizationBuildError(f"invalid JSON {path}: {exc}") from exc
    require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def iter_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            require(line.endswith("\n"), f"{path}:{line_number} lacks newline")
            require(line.strip() != "", f"{path}:{line_number} is blank")
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise FactorizationBuildError(
                    f"{path}:{line_number} invalid JSON: {exc}"
                ) from exc
            require(isinstance(value, dict), f"{path}:{line_number} is not object")
            yield value


def write_new_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def rename_no_replace(source: Path, target: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    require(renameat2 is not None, "libc.renameat2 is unavailable")
    renameat2.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    renameat2.restype = ctypes.c_int
    result = renameat2(
        -100,
        os.fsencode(source),
        -100,
        os.fsencode(target),
        1,
    )
    if result != 0:
        error = ctypes.get_errno()
        raise FactorizationBuildError(
            f"atomic no-replace rename failed: {source} -> {target}: "
            f"{os.strerror(error)}"
        )


def verify_identity_spec(spec: dict[str, Any]) -> None:
    path = Path(str(spec["path"]))
    observed = identity(path)
    require(
        observed["size_bytes"] == int(spec["size_bytes"])
        and observed["sha256"] == str(spec["sha256"]),
        f"published identity mismatch: {path}",
    )


def validate_upstream_authority() -> dict[str, dict[str, Any]]:
    cohort = load_json(TOPOLOGY_RECEIPT)
    summary = load_json(TOPOLOGY_SUMMARY)
    census = load_json(CENSUS_RECEIPT)
    require(
        (
            cohort.get("schema_name"),
            cohort.get("schema_version"),
        )
        == (
            "intersubmod.exact_ps_cpp_topology_cohort.cohort_receipt",
            "1.0.0",
        ),
        "unsupported exact-PS cohort receipt",
    )
    require(
        (
            summary.get("schema_name"),
            summary.get("schema_version"),
        )
        == (
            "intersubmod.exact_ps_cpp_topology_af.cohort_summary",
            "1.0.0",
        ),
        "unsupported exact-PS cohort summary",
    )
    require(
        (
            census.get("schema_name"),
            census.get("schema_version"),
        )
        == (
            "intersubmod.exact_ps_cpp_topology_signature_census.receipt",
            "1.0.0",
        ),
        "unsupported signature-census receipt",
    )
    cohort_scope = cohort.get("scope") or {}
    runner_scope = (summary.get("scope") or {}).get("runner_scope") or {}
    require(
        cohort.get("all_pass") is True
        and cohort.get("technical_all_pass") is True
        and cohort.get("validation_evidence_eligible") is False
        and cohort.get("all_mutation_bearing_families_complete") is False,
        "upstream exact-PS claim ceiling drifted",
    )
    require(
        summary.get("all_pass") is True
        and summary.get("technical_all_pass") is True
        and summary.get("validation_evidence_eligible") is False
        and summary.get("all_mutation_bearing_families_complete") is False,
        "upstream summary claim ceiling drifted",
    )
    require(
        cohort_scope.get("samples") == list(EXPECTED_SAMPLES)
        and runner_scope.get("samples") == list(EXPECTED_SAMPLES)
        and cohort_scope.get("chromosomes") == EXPECTED_CHROMOSOMES
        and runner_scope.get("chromosomes") == EXPECTED_CHROMOSOMES
        and cohort_scope.get("autosomes_complete") is True
        and cohort_scope.get("canonical_all7") is True,
        "upstream exact-PS all7 chr1-22 scope drifted",
    )
    parameters = cohort.get("parameters") or {}
    require(
        int(parameters.get("threshold", -1)) == 3
        and int(parameters.get("max_family_size", -1)) == 100_000
        and int(parameters.get("max_search_nodes", -1)) == 1_000,
        "upstream solver/linkage parameters drifted",
    )
    totals = summary.get("totals") or {}
    require(
        int(totals.get("groups_total", -1)) == 98_955
        and int(totals.get("mutation_bearing_units", -1)) == 85_941
        and int(totals.get("ranked_units", -1)) == EXPECTED_RANKED
        and int(totals.get("mutation_family_abstain_units", -1)) == 10_717
        and int(totals.get("objective_abstain_units", -1)) == 10_717,
        "upstream exact-PS denominators drifted",
    )
    census_scope = census.get("scope") or {}
    census_checks = census.get("checks") or {}
    require(
        census_scope.get("datasets") == list(EXPECTED_SAMPLES)
        and int(census_scope.get("ranked_units", -1)) == EXPECTED_RANKED
        and int(census_scope.get("global_best_trees_enumerated", -1))
        == EXPECTED_GLOBAL_BEST_TREES
        and census_checks.get("all_pass") is True,
        "signature-census all7 authority drifted",
    )
    return {
        "topology_cohort_receipt": identity(TOPOLOGY_RECEIPT),
        "topology_all7_summary": identity(TOPOLOGY_SUMMARY),
        "signature_census_receipt": identity(CENSUS_RECEIPT),
    }


def compile_binary(implementation_dir: Path) -> tuple[Path, list[str], float]:
    require(CPP_SOURCE.is_file(), f"missing C++ source: {CPP_SOURCE}")
    require(FROZEN_SOURCE.is_file(), f"missing frozen source: {FROZEN_SOURCE}")
    implementation_dir.mkdir(parents=True, exist_ok=True)
    binary = implementation_dir / "exact_topology_candidate_factorization"
    command = [
        "g++",
        "-std=c++17",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Werror",
        f"-I{LONG_LINEAGE / 'include'}",
        str(CPP_SOURCE),
        str(OBLIGATION_SOURCE),
        str(PARENT_MAPPING_SOURCE),
        "-ljansson",
        "-lcrypto",
        "-o",
        str(binary),
    ]
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    elapsed = time.perf_counter() - started
    require(
        completed.returncode == 0,
        "compile failed:\n" + completed.stdout + completed.stderr,
    )
    require(binary.is_file(), f"compile did not produce {binary}")
    return binary, command, elapsed


def validate_sample_output(sample: str, path: Path) -> dict[str, Any]:
    ranked = 0
    minimum_sets = 0
    minimum_trees = 0
    global_best_trees = 0
    target: dict[str, Any] | None = None
    seen_indices: set[int] = set()
    for row in iter_jsonl(path):
        require(
            (row.get("schema_name"), row.get("schema_version")) == UNIT_SCHEMA,
            f"{sample} unexpected unit schema",
        )
        require(row.get("sample") == sample, f"{sample} row sample mismatch")
        index = int(row.get("group_index", -1))
        require(index >= 0 and index not in seen_indices, f"{sample} duplicate index")
        seen_indices.add(index)
        sets = row.get("minimum_vertex_sets")
        require(isinstance(sets, list) and sets, f"{sample}:{index} empty family")
        require(
            len(sets) == int(row.get("minimum_vertex_set_count", -1)),
            f"{sample}:{index} minimum-set count mismatch",
        )
        require(
            row.get("canonical_reproduction_pass") is True
            and row.get("census_reproduction_pass") is True,
            f"{sample}:{index} reproduction did not pass",
        )
        all_edge_incidence: Counter[tuple[int, int]] = Counter()
        best_edge_incidence: Counter[tuple[int, int]] = Counter()
        global_best_score = Fraction(str(row["best_score_fraction"]))
        best_sets = sum(item.get("is_global_best") is True for item in sets)
        require(
            best_sets >= 1
            and best_sets == int(row.get("best_vertex_set_count", -1)),
            f"{sample}:{index} global-best set count mismatch",
        )
        active_count = int(row.get("active_bit_count", -1))
        active_positions = row.get("active_positions") or []
        active_indices = row.get("active_original_indices") or []
        observed_vertices = [int(value) for value in row.get("observed_vertices") or []]
        require(
            active_count >= 1
            and len(active_positions) == active_count
            and len(active_indices) == active_count
            and observed_vertices == sorted(set(observed_vertices))
            and all(value >= 0 for value in observed_vertices),
            f"{sample}:{index} active/observed vertex metadata mismatch",
        )
        row_minimum_trees = 0
        row_best_trees = 0
        expected_all_edge_total = 0
        expected_best_edge_total = 0
        for set_index, item in enumerate(sets):
            vertices = [int(value) for value in item.get("vertices") or []]
            require(
                vertices == sorted(set(vertices))
                and vertices
                and vertices[0] == 0,
                f"{sample}:{index}:{set_index} invalid vertex set",
            )
            vertex_values = set(vertices)
            require(
                set(observed_vertices) <= vertex_values
                and max(vertices) < (1 << active_count),
                f"{sample}:{index}:{set_index} observed/bit-range mismatch",
            )
            parent_rows = item.get("parent_factorization") or []
            require(
                len(parent_rows) == len(vertices) - 1,
                f"{sample}:{index}:{set_index} child factorization incomplete",
            )
            child_values = set()
            tree_product = 1
            best_product = 1
            legal_edges = set()
            for parent_row in parent_rows:
                require(
                    isinstance(parent_row, list) and len(parent_row) == 3,
                    f"{sample}:{index}:{set_index} malformed parent row",
                )
                child = int(parent_row[0])
                legal = [int(value) for value in parent_row[1]]
                best = [int(value) for value in parent_row[2]]
                require(
                    child in vertex_values
                    and child != 0
                    and child not in child_values
                    and legal == sorted(set(legal))
                    and best == sorted(set(best))
                    and legal
                    and best
                    and set(best) <= set(legal),
                    f"{sample}:{index}:{set_index} invalid parent choices",
                )
                child_values.add(child)
                for parent in legal:
                    require(
                        parent in vertex_values
                        and parent < child
                        and (parent & child) == parent
                        and bin(parent ^ child).count("1") == 1,
                        f"{sample}:{index}:{set_index} illegal hypercube edge",
                    )
                    legal_edges.add((parent, child))
                tree_product *= len(legal)
                best_product *= len(best)
            require(
                child_values == vertex_values - {0},
                f"{sample}:{index}:{set_index} missing/extra child",
            )
            score_rows = item.get("legal_edge_score_fractions") or []
            require(
                all(
                    isinstance(value, list) and len(value) == 3
                    for value in score_rows
                )
                and
                {(int(value[0]), int(value[1])) for value in score_rows}
                == legal_edges
                and all(
                    Fraction(str(value[2])).denominator > 0
                    for value in score_rows
                ),
                f"{sample}:{index}:{set_index} legal edge score coverage mismatch",
            )
            score_map = {
                (int(value[0]), int(value[1])): Fraction(str(value[2]))
                for value in score_rows
            }
            reproduced_score = Fraction(0, 1)
            for child, legal, best in parent_rows:
                maximum = max(
                    score_map[(int(parent), int(child))] for parent in legal
                )
                reproduced_best = sorted(
                    int(parent)
                    for parent in legal
                    if score_map[(int(parent), int(child))] == maximum
                )
                require(
                    reproduced_best == [int(parent) for parent in best],
                    f"{sample}:{index}:{set_index} best-parent score mismatch",
                )
                reproduced_score += maximum
            set_trees = int(item["total_tree_count"])
            set_best_trees = int(item["best_tree_count"])
            require(
                set_trees == tree_product
                and set_best_trees == best_product,
                f"{sample}:{index}:{set_index} parent-product mismatch",
            )
            for child, legal, best in parent_rows:
                for parent in legal:
                    all_edge_incidence[(int(parent), int(child))] += (
                        set_trees // len(legal)
                    )
                if item.get("is_global_best") is True:
                    for parent in best:
                        best_edge_incidence[(int(parent), int(child))] += (
                            set_best_trees // len(best)
                        )
            set_score = Fraction(str(item["best_score_fraction"]))
            require(
                set_score == reproduced_score
                and
                (item.get("is_global_best") is True)
                == (set_score == global_best_score),
                f"{sample}:{index}:{set_index} global-best flag/score mismatch",
            )
            row_minimum_trees += set_trees
            expected_all_edge_total += set_trees * (len(vertices) - 1)
            if item.get("is_global_best") is True:
                row_best_trees += set_best_trees
                expected_best_edge_total += set_best_trees * (len(vertices) - 1)
        require(
            row_minimum_trees == int(row["total_tree_count"]),
            f"{sample}:{index} minimum-tree conservation failed",
        )
        require(
            row_best_trees == int(row["best_tree_tie_count"]),
            f"{sample}:{index} best-tree conservation failed",
        )
        require(
            sum(int(item["tree_count"]) for item in row["global_best_signatures"])
            == row_best_trees,
            f"{sample}:{index} signature conservation failed",
        )
        observed_all_rows = row["all_minimum_tree_edge_incidence"]
        observed_best_rows = row["global_best_edge_incidence"]
        require(
            all(isinstance(value, list) and len(value) == 3 for value in observed_all_rows)
            and all(
                isinstance(value, list) and len(value) == 3
                for value in observed_best_rows
            ),
            f"{sample}:{index} malformed aggregate edge incidence",
        )
        all_edge_keys = [
            (int(parent), int(child)) for parent, child, _ in observed_all_rows
        ]
        best_edge_keys = [
            (int(parent), int(child)) for parent, child, _ in observed_best_rows
        ]
        require(
            len(all_edge_keys) == len(set(all_edge_keys))
            and len(best_edge_keys) == len(set(best_edge_keys)),
            f"{sample}:{index} duplicate aggregate edge incidence",
        )
        observed_all_edges = Counter(
            {
                (int(parent), int(child)): int(count)
                for parent, child, count in observed_all_rows
            }
        )
        observed_best_edges = Counter(
            {
                (int(parent), int(child)): int(count)
                for parent, child, count in observed_best_rows
            }
        )
        require(
            observed_all_edges == all_edge_incidence
            and observed_best_edges == best_edge_incidence
            and set(observed_best_edges) <= set(observed_all_edges),
            f"{sample}:{index} generic edge-incidence conservation failed",
        )
        require(
            sum(observed_all_edges.values()) == expected_all_edge_total
            and sum(observed_best_edges.values()) == expected_best_edge_total,
            f"{sample}:{index} aggregate edge-incidence denominator failed",
        )
        coarse_counts = Counter(
            {
                str(key): int(value)
                for key, value in (
                    row.get("global_best_coarse_class_tree_counts") or {}
                ).items()
            }
        )
        signature_coarse = Counter()
        for item in row["global_best_signatures"]:
            signature_coarse[str(item["coarse_class"])] += int(
                item["tree_count"]
            )
        require(
            signature_coarse == coarse_counts,
            f"{sample}:{index} coarse/signature conservation failed",
        )
        ranked += 1
        minimum_sets += len(sets)
        minimum_trees += row_minimum_trees
        global_best_trees += row_best_trees
        if sample == "HCC1395" and index == 7469:
            target = row

    if sample == "HCC1395":
        require(target is not None, "HCC1395 target group 7469 missing")
        require(
            target["active_positions"] == [87818272, 87840023],
            "HCC1395 target active positions drifted",
        )
        require(
            [item["vertices"] for item in target["minimum_vertex_sets"]]
            == [[0, 1, 2], [0, 1, 3], [0, 2, 3]],
            "HCC1395 target minimum family drifted",
        )
        require(
            target["all_minimum_tree_edge_incidence"]
            == [[0, 1, "2"], [0, 2, "2"], [1, 3, "1"], [2, 3, "1"]],
            "HCC1395 target candidate edge union drifted",
        )
        require(
            target["global_best_edge_incidence"]
            == [[0, 1, "1"], [1, 3, "1"]],
            "HCC1395 target selected edge overlay drifted",
        )
    result = {
        "ranked_units": ranked,
        "minimum_vertex_sets": minimum_sets,
        "minimum_trees": minimum_trees,
        "global_best_trees": global_best_trees,
    }
    require(
        tuple(result.values()) == EXPECTED_SAMPLE_COUNTS[sample],
        f"{sample} exact fixture drifted: {result}",
    )
    return result


def run_sample(
    binary: Path, sample: str, output: Path
) -> tuple[dict[str, Any], list[str], str, float]:
    sample_root = TOPOLOGY_ROOT / "samples" / sample
    mlhp = sample_root / f"{sample}.exact_ps_mlhp.json"
    canonical = sample_root / f"{sample}.topology.jsonl"
    census = CENSUS_ROOT / f"{sample}.census.jsonl"
    for path in (mlhp, canonical, census):
        require(path.is_file(), f"missing {sample} input: {path}")
    mlhp_document = load_json(mlhp)
    require(
        (
            mlhp_document.get("schema_name"),
            mlhp_document.get("schema_version"),
        )
        == ("intersubmod.exact_ps_layered_tree_input", "1.0.0")
        and mlhp_document.get("sample") == sample,
        f"{sample} MLHP schema/sample drifted",
    )
    require(
        (mlhp_document.get("scope") or {}).get("chromosomes")
        == EXPECTED_CHROMOSOMES
        and int((mlhp_document.get("params") or {}).get("MINREAD", -1)) == 3
        and (mlhp_document.get("params") or {}).get("region_linkage_rule")
        == "strict fixed-R/A endpoint-pair support connected component",
        f"{sample} MLHP chr1-22/MINREAD/linkage semantics drifted",
    )
    input_paths = {
        "mlhp": mlhp,
        "canonical_topology": canonical,
        "signature_census": census,
    }
    pre_identities = {
        key: identity(path) for key, path in input_paths.items()
    }
    command = [
        str(binary),
        "--input",
        str(mlhp),
        "--canonical",
        str(canonical),
        "--census",
        str(census),
        "--output",
        str(output),
        "--max-family-size",
        "100000",
        "--max-search-nodes",
        "1000",
    ]
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    elapsed = time.perf_counter() - started
    require(
        completed.returncode == 0,
        f"{sample} factorization failed:\n{completed.stdout}{completed.stderr}",
    )
    require(output.is_file(), f"{sample} output missing: {output}")
    counts = validate_sample_output(sample, output)
    post_identities = {
        key: identity(path) for key, path in input_paths.items()
    }
    require(
        pre_identities == post_identities,
        f"{sample} input identity changed during factorization",
    )
    counts["elapsed_seconds"] = round(elapsed, 6)
    counts["inputs"] = pre_identities
    counts["output"] = identity(output)
    return counts, command, completed.stderr.strip(), elapsed


def build(output_root: Path) -> dict[str, Any]:
    output_root = output_root.resolve()
    require(not output_root.exists(), f"output root already exists: {output_root}")
    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = output_root.parent / (
        f".{output_root.name}.pending.{os.getpid()}.{time.time_ns()}"
    )
    require(not staging.exists(), f"staging root already exists: {staging}")
    staging.mkdir()
    try:
        upstream_pre = validate_upstream_authority()
        sources_pre = {
            key: identity(path) for key, path in IMPLEMENTATION_PATHS.items()
        }
        compiler = subprocess.run(
            ["g++", "--version"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        require(
            compiler.returncode == 0 and compiler.stdout.strip() != "",
            "cannot capture compiler version",
        )
        binary, compile_command, compile_elapsed = compile_binary(
            staging / "implementation"
        )
        samples: dict[str, Any] = {}
        run_commands: dict[str, dict[str, Any]] = {}
        stderr: dict[str, str] = {}
        for sample in EXPECTED_SAMPLES:
            sample_output = staging / f"{sample}.candidate_factorization.jsonl"
            counts, command, sample_stderr, _ = run_sample(
                binary, sample, sample_output
            )
            counts["output"]["path"] = str(
                output_root / sample_output.name
            )
            samples[sample] = counts
            replay_command = list(command)
            replay_command[0] = str(
                output_root
                / "implementation"
                / "exact_topology_candidate_factorization"
            )
            replay_command[replay_command.index("--output") + 1] = (
                "<NEW_OUTPUT_JSONL>"
            )
            run_commands[sample] = {
                "executed_command": command,
                "published_output": str(output_root / sample_output.name),
                "replay_command_template": replay_command,
            }
            stderr[sample] = sample_stderr
        upstream_post = validate_upstream_authority()
        sources_post = {
            key: identity(path) for key, path in IMPLEMENTATION_PATHS.items()
        }
        require(
            upstream_pre == upstream_post,
            "upstream authority changed during all7 factorization",
        )
        require(
            sources_pre == sources_post,
            "implementation source/header identity changed during build",
        )
        totals = {
            key: sum(int(row[key]) for row in samples.values())
            for key in (
                "ranked_units",
                "minimum_vertex_sets",
                "minimum_trees",
                "global_best_trees",
            )
        }
        checks = {
            "sample_scope_exact": tuple(samples) == EXPECTED_SAMPLES,
            "ranked_units_exact": totals["ranked_units"] == EXPECTED_RANKED,
            "minimum_trees_exact": (
                totals["minimum_trees"] == EXPECTED_MINIMUM_TREES
            ),
            "global_best_trees_exact": (
                totals["global_best_trees"] == EXPECTED_GLOBAL_BEST_TREES
            ),
            "current_all_minimum_sets_have_one_parent_mapping": (
                totals["minimum_vertex_sets"] == totals["minimum_trees"]
            ),
            "all_outputs_nonempty": all(
                row["output"]["size_bytes"] > 0 for row in samples.values()
            ),
            "per_sample_exact_fixtures": all(
                tuple(
                    samples[sample][key]
                    for key in (
                        "ranked_units",
                        "minimum_vertex_sets",
                        "minimum_trees",
                        "global_best_trees",
                    )
                )
                == EXPECTED_SAMPLE_COUNTS[sample]
                for sample in EXPECTED_SAMPLES
            ),
            "upstream_claim_ceiling_inherited": True,
            "source_and_input_pre_post_identity_exact": True,
        }
        require(all(checks.values()), f"cohort checks failed: {checks}")
        binary_spec = identity(binary)
        binary_spec["path"] = str(
            output_root
            / "implementation"
            / "exact_topology_candidate_factorization"
        )
        receipt = {
            "schema_name": RECEIPT_SCHEMA[0],
            "schema_version": RECEIPT_SCHEMA[1],
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "scope": {
                "datasets": list(EXPECTED_SAMPLES),
                "chromosomes": EXPECTED_CHROMOSOMES,
                "ranked_only": True,
                "exact_ps_primary_hp": True,
                "strict_endpoint_read_linkage_threshold": 3,
            },
            "denominators": {
                "all_exact_ps_hp_units": 98_955,
                "mutation_bearing_units": 85_941,
                "ranked_complete_units_included": EXPECTED_RANKED,
                "resource_abstain_units_excluded": 10_717,
                "no_active_alt_units_excluded": 13_014,
            },
            "claim_boundary": (
                "This bundle exhaustively factorizes minimum vertex sets and "
                "legal parent mappings only for the 71,955 ranked-complete "
                "exact-PS x primary-HP mathematical units. It does not prove "
                "cellular clones, biological ancestry, CN/LOH correctness, "
                "calibrated posterior probability, or completeness of the "
                "10,717 resource-abstain mutation-bearing units."
            ),
            "upstream_authority": upstream_pre,
            "implementation": {
                **sources_pre,
                "binary": binary_spec,
                "compile_command": compile_command,
                "compile_elapsed_seconds": round(compile_elapsed, 6),
                "compiler_version": compiler.stdout.splitlines()[0],
            },
            "samples": samples,
            "totals": totals,
            "checks": checks,
            "commands": run_commands,
            "stderr": stderr,
            "technical_all_pass": True,
            "validation_evidence_eligible": False,
            "all_mutation_bearing_families_complete": False,
            "all_pass": True,
        }
        write_new_json(staging / "receipt.json", receipt)
        rename_no_replace(staging, output_root)
        for row in receipt["samples"].values():
            verify_identity_spec(row["output"])
        verify_identity_spec(receipt["implementation"]["binary"])
        require(
            (output_root / "receipt.json").is_file(),
            "published receipt is missing",
        )
        return receipt
    except Exception as exc:
        if staging.exists():
            failed = output_root.parent / (
                f"{output_root.name}.failed."
                f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}."
                f"{os.getpid()}.{time.time_ns()}"
            )
            try:
                rename_no_replace(staging, failed)
            except Exception as archive_exc:
                raise FactorizationBuildError(
                    f"{exc}; staging preserved at {staging}; "
                    f"failure archive rename also failed: {archive_exc}"
                ) from exc
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        receipt = build(args.output_root)
    except (
        FactorizationBuildError,
        IndexError,
        KeyError,
        OSError,
        TypeError,
        ValueError,
        subprocess.TimeoutExpired,
    ) as exc:
        print(f"FAIL CLOSED: {exc}", file=sys.stderr)
        return 2
    print(f"INPUT_TOPOLOGY_ROOT={TOPOLOGY_ROOT}")
    print(f"INPUT_CENSUS_ROOT={CENSUS_ROOT}")
    print(f"COMMAND_SOURCE={CPP_SOURCE}")
    print(f"OUTPUT={args.output_root.resolve()}")
    print(json.dumps(receipt["totals"], sort_keys=True))
    print(f"ALL_PASS={str(receipt['all_pass']).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
