"""Synthetic correctness tests for the bounded exact-PS C++ topology/AF CLI."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[3]
SOURCE = (
    REPO_ROOT
    / "research/20260724_exact_ps_cpp_topology_af_all_samples/cpp"
    / "exact_ps_topology_af.cpp"
)
LONG_LINEAGE_ROOT = Path(
    os.environ.get("LONG_LINEAGE_ROOT", "/big7_disk/liaoyoyo2001/LongLineage")
)


_BUILD_DIRECTORY: tempfile.TemporaryDirectory[str] | None = None
BINARY: Path | None = None


def setUpModule() -> None:
    global _BUILD_DIRECTORY, BINARY
    _BUILD_DIRECTORY = tempfile.TemporaryDirectory(prefix="exact_ps_cpp_build.")
    output = Path(_BUILD_DIRECTORY.name) / "exact_ps_topology_af"
    command = [
        "g++",
        "-std=c++17",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        f"-I{LONG_LINEAGE_ROOT / 'include'}",
        str(SOURCE),
        str(LONG_LINEAGE_ROOT / "src/solver/obligation_bnb.cpp"),
        str(LONG_LINEAGE_ROOT / "src/solver/parent_mapping.cpp"),
        "-ljansson",
        "-lcrypto",
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True, cwd=REPO_ROOT)
    BINARY = output


def tearDownModule() -> None:
    global _BUILD_DIRECTORY
    if _BUILD_DIRECTORY is not None:
        _BUILD_DIRECTORY.cleanup()


def group(
    *,
    patterns: dict[str, int] | None = None,
    partials: dict[str, int] | None = None,
    coverage: list[tuple[int, int]] | None = None,
    positions: list[int] | None = None,
    region: str = "chr1|PS=100|HP=1|U1:B0001",
) -> dict:
    patterns = patterns or {}
    partials = partials or {}
    inferred_k = len(next(iter(patterns or partials)))
    positions = positions or [101 + index for index in range(inferred_k)]
    coverage = coverage or [(5, 5) for _ in positions]
    assert len(positions) == len(coverage)
    return {
        "analysis_unit": "exact_ps_hp_component_bounded_block",
        "block_id": "U1:B0001",
        "block_index": 1,
        "chrom": "chr1",
        "cn": "unavailable",
        "cn_data_available": False,
        "col_coverage_by_hp": {
            "1": {
                str(position): [counts[0], counts[1]]
                for position, counts in zip(positions, coverage)
            }
        },
        "component_id": "S:chr1:PS_HP1:PS1:B3:1:101-200",
        "coverage_interpretation": "synthetic test",
        "dataset": "S",
        "end": positions[-1],
        "hp_family": "1",
        "linkage_basis": "PS_HP1",
        "min_read": 3,
        "n_full_cov_reads": sum(patterns.values()),
        "n_sSNV": len(positions),
        "phase_set": "100",
        "phase_set_status": "KNOWN_PS_PRIMARY",
        "populations_by_hp": {"1": patterns} if patterns else {},
        "positions": positions,
        "projected_molecule_block_incidences": 10,
        "projected_molecule_rows": 10,
        "reads_by_hp": {"1": 10},
        "region_id": region,
        "retained_segmentation_constraint_weight": 10,
        "span": positions[-1] - positions[0],
        "start": positions[0],
        "subread_groups_by_hp": {"1": partials} if partials else {},
        "tree_supported_molecule_block_incidences": 10,
        "tree_supported_pattern_count": len(patterns) + len(partials),
        "unit_id": "U1",
        "vaf_eligible": False,
    }


def document(groups: list[dict]) -> dict:
    return {
        "claim_ceiling": "synthetic",
        "groups": groups,
        "n_groups_analyzed": len(groups),
        "sample": "S",
        "schema_name": "intersubmod.exact_ps_layered_tree_input",
        "schema_version": "1.0.0",
    }


def run_cli(
    binary: Path,
    tmp_path: Path,
    payload: dict,
    *,
    suffix: str = "",
    extra: list[str] | None = None,
) -> tuple[subprocess.CompletedProcess[str], list[dict], dict | None]:
    input_path = tmp_path / f"input{suffix}.json"
    output_path = tmp_path / f"output{suffix}.jsonl"
    receipt_path = tmp_path / f"receipt{suffix}.json"
    input_path.write_text(json.dumps(payload), encoding="utf-8")
    command = [
        str(binary),
        "--input",
        str(input_path),
        "--output",
        str(output_path),
        "--receipt",
        str(receipt_path),
        *(extra or []),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, cwd=REPO_ROOT)
    rows = (
        [json.loads(line) for line in output_path.read_text().splitlines()]
        if output_path.exists()
        else []
    )
    receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else None
    return completed, rows, receipt


def active_problem(
    patterns: list[str], partials: list[str]
) -> tuple[int, set[int], list[set[int]]]:
    k = len(next(iter(patterns or partials)))
    active = [
        index
        for index in range(k)
        if any(value[index] == "A" for value in patterns + partials)
    ]

    def mask(value: str) -> int:
        return sum(
            1 << compact
            for compact, original in enumerate(active)
            if value[original] == "A"
        )

    mandatory = {0, *(mask(value) for value in patterns)}
    groups: list[set[int]] = []
    for value in partials:
        fixed = mask(value)
        unknown = [
            compact
            for compact, original in enumerate(active)
            if value[original] == "X"
        ]
        groups.append(
            {
                fixed
                | sum(
                    1 << unknown[offset]
                    for offset in range(len(unknown))
                    if subset & (1 << offset)
                )
                for subset in range(1 << len(unknown))
            }
        )
    return len(active), mandatory, groups


def exhaustive_oracle(
    patterns: list[str], partials: list[str], af: list[Fraction] | None = None
) -> dict:
    k, mandatory, groups = active_problem(patterns, partials)
    vertices = range(1 << k)
    feasible: list[set[int]] = []
    for selected_bits in range(1 << (1 << k)):
        selected = {
            vertex for vertex in vertices if selected_bits & (1 << vertex)
        }
        if not mandatory <= selected or any(not (selected & group) for group in groups):
            continue
        if any(
            vertex
            and not any(
                (vertex ^ (1 << bit)) in selected
                for bit in range(k)
                if vertex & (1 << bit)
            )
            for vertex in selected
        ):
            continue
        feasible.append(selected)
    minimum_size = min(map(len, feasible))
    family = sorted(
        (tuple(sorted(value)) for value in feasible if len(value) == minimum_size)
    )
    tree_count = 0
    scored: list[Fraction] = []
    for selected_tuple in family:
        selected = set(selected_tuple)
        parent_lists = [
            [
                vertex ^ (1 << bit)
                for bit in range(k)
                if vertex & (1 << bit)
                and (vertex ^ (1 << bit)) in selected
            ]
            for vertex in selected_tuple
            if vertex
        ]
        tree_count += __import__("math").prod(map(len, parent_lists))
        if af is not None:
            for parents in product(*parent_lists):
                score = Fraction()
                children = [vertex for vertex in selected_tuple if vertex]
                for parent, child in zip(parents, children):
                    acquired = (child ^ parent).bit_length() - 1
                    score += sum(
                        af[bit] - af[acquired]
                        for bit in range(k)
                        if parent & (1 << bit)
                    )
                scored.append(score)
    canonical = (
        f"active_bit_count={k}\n"
        + "".join(",".join(map(str, selected)) + "\n" for selected in family)
    )
    return {
        "k": k,
        "family": family,
        "family_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "tree_count": tree_count,
        "best_score": max(scored) if scored else None,
        "best_ties": scored.count(max(scored)) if scored else None,
    }


class ExactPsTopologyAfTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory(prefix="exact_ps_cpp_test.")
        self.tmp_path = Path(self._temporary.name)
        assert BINARY is not None
        self.binary = BINARY

    def tearDown(self) -> None:
        self._temporary.cleanup()

    def test_diamond_is_counted_exactly_but_recurrence_is_not_ranked(self) -> None:
        value = group(patterns={"AR": 3, "RA": 3, "AA": 3})
        completed, rows, receipt = run_cli(
            self.binary, self.tmp_path, document([value])
        )
        assert completed.returncode == 0, completed.stderr
        assert receipt and receipt["all_pass"] is True
        assert receipt["all_units_well_formed"] is True
        row = rows[0]
        oracle = exhaustive_oracle(["AR", "RA", "AA"], [])
        assert row["active_bit_count"] == oracle["k"] == 2
        assert row["minimum_vertex_set_count"] == len(oracle["family"]) == 1
        assert row["minimum_family_sha256"] == oracle["family_sha256"]
        assert row["total_tree_count"] == str(oracle["tree_count"]) == "2"
        assert row["recurrence_required"] is True
        assert row["read_af_status"] == "recurrence_not_evaluable"
        assert row["best_tree_tie_count"] is None

    def test_cross_vertex_set_af_tie_is_exact_and_repeatable(self) -> None:
        value = group(partials={"AX": 3, "XA": 3})
        completed, rows, receipt = run_cli(
            self.binary, self.tmp_path, document([value])
        )
        assert completed.returncode == 0, completed.stderr
        row = rows[0]
        oracle = exhaustive_oracle(
            [], ["AX", "XA"], [Fraction(1, 2), Fraction(1, 2)]
        )
        assert row["minimum_vertex_set_count"] == len(oracle["family"]) == 3
        assert row["minimum_family_sha256"] == oracle["family_sha256"]
        assert row["total_tree_count"] == str(oracle["tree_count"]) == "3"
        assert row["best_score_fraction"] == "0/1"
        assert row["best_tree_tie_count"] == str(oracle["best_ties"]) == "3"
        assert row["best_tree_unique"] is False
        assert row["best_vertex_set_count"] == 3
        assert row["solver_elapsed_microseconds"] >= 0
        assert receipt and receipt["ranking_complete_count"] == 1

        completed2, rows2, _ = run_cli(
            self.binary, self.tmp_path, document([value]), suffix="_repeat"
        )
        assert completed2.returncode == 0, completed2.stderr
        stable_fields = [
            "minimum_family_sha256",
            "total_tree_count",
            "best_score_fraction",
            "best_tree_tie_count",
            "representative_best_edges",
            "representative_best_morphology",
        ]
        assert {key: rows[0][key] for key in stable_fields} == {
            key: rows2[0][key] for key in stable_fields
        }

    def test_af_direction_selects_unique_chain(self) -> None:
        value = group(
            partials={"AX": 3, "XA": 3},
            coverage=[(2, 8), (4, 6)],
        )
        completed, rows, _ = run_cli(
            self.binary, self.tmp_path, document([value])
        )
        assert completed.returncode == 0, completed.stderr
        row = rows[0]
        oracle = exhaustive_oracle(
            [], ["AX", "XA"], [Fraction(8, 10), Fraction(6, 10)]
        )
        assert row["best_score_fraction"] == (
            f"{oracle['best_score'].numerator}/{oracle['best_score'].denominator}"
        )
        assert row["best_tree_tie_count"] == str(oracle["best_ties"]) == "1"
        assert row["best_tree_unique"] is True
        assert [
            (edge["parent_label"], edge["child_label"])
            for edge in row["representative_best_edges"]
        ] == [("ROOT", "H_AR"), ("H_AR", "H_AA")]

    def test_zero_denominator_fails_closed_after_complete_structure(self) -> None:
        value = group(
            partials={"AX": 3, "XA": 3}, coverage=[(0, 0), (5, 5)]
        )
        completed, rows, receipt = run_cli(
            self.binary, self.tmp_path, document([value])
        )
        assert completed.returncode == 0, completed.stderr
        row = rows[0]
        assert row["family_status"] == "FAMILY_COMPLETE"
        assert row["read_af_status"] == "zero_denominator"
        assert row["best_score_fraction"] is None
        assert row["best_tree_unique"] is None
        assert receipt and receipt["all_pass"] is True
        assert receipt["all_units_well_formed"] is True
        assert receipt["ranking_complete_count"] == 0

    def test_family_cap_certifies_objective_but_withholds_family(self) -> None:
        value = group(patterns={"AAAA": 3})
        completed, rows, receipt = run_cli(
            self.binary,
            self.tmp_path,
            document([value]),
            extra=["--max-family-size", "2"],
        )
        assert completed.returncode == 0, completed.stderr
        row = rows[0]
        assert row["objective_status"] == "OBJECTIVE_CERTIFIED"
        assert row["objective_h"] == 3
        assert row["family_status"] == "FAMILY_INCOMPLETE_CAP"
        assert row["minimum_vertex_set_count"] is None
        assert row["total_tree_count"] is None
        assert row["read_af_status"] == "not_evaluable_family_incomplete"
        assert receipt and receipt["all_pass"] is True
        assert receipt["all_units_objective_certified"] is True
        assert receipt["all_mutation_bearing_families_complete"] is False

    def test_malformed_unit_is_explicit_abstention_not_run_failure(self) -> None:
        value = group(patterns={"AR": 3})
        value["populations_by_hp"] = {"1": {"AZ": 3}}
        completed, rows, receipt = run_cli(
            self.binary, self.tmp_path, document([value])
        )
        assert completed.returncode == 0, completed.stderr
        assert rows[0]["unit_status"] == "malformed_unit"
        assert rows[0]["objective_status"] == "ABSTAIN_NOT_IDENTIFIABLE"
        assert rows[0]["best_score_fraction"] is None
        assert receipt and receipt["all_pass"] is False
        assert receipt["all_units_well_formed"] is False
        assert receipt["all_units_objective_certified"] is False
        assert receipt["status_census"]["unit_status"] == {"malformed_unit": 1}

    def test_malformed_top_level_publishes_no_receipt(self) -> None:
        payload = document([])
        payload["schema_version"] = "broken"
        completed, rows, receipt = run_cli(self.binary, self.tmp_path, payload)
        assert completed.returncode == 1
        assert rows == []
        assert receipt is None
        assert "unsupported input schema" in completed.stderr


if __name__ == "__main__":
    unittest.main(verbosity=2)
