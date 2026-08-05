#!/usr/bin/env python3
"""Contract tests for the independent C++ exact-PS bounded partition kernel."""

from __future__ import annotations

import csv
from decimal import Decimal
import importlib.util
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE = ROOT / "cpp" / "exact_ps_k12_partition.cpp"
REFERENCE_SCRIPTS = ROOT.parents[1] / "research" / "20260718_k_gt8_read_supported_segmentation" / "scripts"
PRODUCER_PATH = ROOT / "scripts" / "exact_ps_k12_partition.py"
sys.path.insert(0, str(REFERENCE_SCRIPTS))

from read_support_partition import Hyperedge, partition_ordered_hypergraph  # noqa: E402


def load_python_producer():
    spec = importlib.util.spec_from_file_location("exact_ps_k12_python_producer", PRODUCER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import Python producer: {PRODUCER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module

UNIT_HEADER = (
    "dataset",
    "chrom",
    "unit_id",
    "component_id",
    "linkage_basis",
    "hp_family",
    "phase_set",
    "threshold",
    "k",
    "positions1",
)
CONSTRAINT_HEADER = (
    "dataset",
    "chrom",
    "unit_id",
    "constraint_id",
    "hp_family",
    "phase_set",
    "positions1",
    "call_codes",
    "molecule_weight",
    "pattern_count",
)
OUTPUT_NAMES = (
    "blocks.tsv",
    "membership.tsv",
    "dispositions.tsv",
    "summary.json",
)


def unit_row(
    unit_id: str,
    positions: tuple[int, ...],
    *,
    dataset: str = "HCC1395",
    chrom: str = "chr1",
    hp: str = "1",
    ps: str = "100",
) -> dict[str, str]:
    return {
        "dataset": dataset,
        "chrom": chrom,
        "unit_id": unit_id,
        "component_id": f"component:{unit_id}",
        "linkage_basis": f"PS_HP{hp}",
        "hp_family": hp,
        "phase_set": ps,
        "threshold": "3",
        "k": str(len(positions)),
        "positions1": ",".join(map(str, positions)),
    }


def constraint_row(
    unit_id: str,
    constraint_id: str,
    positions: tuple[int, ...],
    *,
    dataset: str = "HCC1395",
    chrom: str = "chr1",
    hp: str = "1",
    ps: str = "100",
    weight: str = "1",
) -> dict[str, str]:
    return {
        "dataset": dataset,
        "chrom": chrom,
        "unit_id": unit_id,
        "constraint_id": constraint_id,
        "hp_family": hp,
        "phase_set": ps,
        "positions1": ",".join(map(str, positions)),
        "call_codes": "A" * len(positions),
        "molecule_weight": weight,
        "pattern_count": "1",
    }


def write_tsv(path: pathlib.Path, header: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ExactPsK12CppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._build_tmp = tempfile.TemporaryDirectory()
        cls.binary = pathlib.Path(cls._build_tmp.name) / "exact_ps_k12_partition"
        command = [
            "g++",
            "-std=c++17",
            "-O2",
            "-Wall",
            "-Wextra",
            "-Werror",
            str(SOURCE),
            "-o",
            str(cls.binary),
        ]
        subprocess.run(command, check=True, text=True, capture_output=True)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._build_tmp.cleanup()

    def run_kernel(
        self,
        root: pathlib.Path,
        units: list[dict[str, str]],
        constraints: list[dict[str, str]],
        *,
        max_block_size: int = 12,
        label: str = "run",
    ) -> tuple[subprocess.CompletedProcess[str], pathlib.Path]:
        units_path = root / f"{label}.units.tsv"
        constraints_path = root / f"{label}.constraints.tsv"
        output_dir = root / f"{label}.output"
        write_tsv(units_path, UNIT_HEADER, units)
        write_tsv(constraints_path, CONSTRAINT_HEADER, constraints)
        process = subprocess.run(
            [
                str(self.binary),
                "--units",
                str(units_path),
                "--constraints",
                str(constraints_path),
                "--output-dir",
                str(output_dir),
                "--max-block-size",
                str(max_block_size),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        return process, output_dir

    def assert_contract_failure(
        self,
        units: list[dict[str, str]],
        constraints: list[dict[str, str]],
        expected: str,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            process, _ = self.run_kernel(pathlib.Path(temporary), units, constraints)
            self.assertNotEqual(process.returncode, 0)
            self.assertIn(expected, process.stderr)

    def test_rejects_cross_ps_constraint(self) -> None:
        positions = (10, 20, 30)
        self.assert_contract_failure(
            [unit_row("U1", positions)],
            [constraint_row("U1", "X1", (10, 20), ps="200")],
            "phase_set mismatch",
        )

    def test_rejects_duplicate_unit_site(self) -> None:
        self.assert_contract_failure(
            [unit_row("U1", (10, 20, 20))],
            [],
            "positions1 must be strictly increasing",
        )

    def test_rejects_unknown_constraint_site(self) -> None:
        self.assert_contract_failure(
            [unit_row("U1", (10, 20, 30))],
            [constraint_row("U1", "X1", (10, 40))],
            "references unknown unit position",
        )

    def test_rejects_duplicate_ids_and_nonunit_pattern_count(self) -> None:
        positions = (10, 20)
        duplicate_units = [unit_row("U1", positions), unit_row("U1", positions, ps="200")]
        self.assert_contract_failure(duplicate_units, [], "duplicate unit_id")

        duplicate_constraints = [
            constraint_row("U1", "C1", (10,)),
            constraint_row("U1", "C1", (20,)),
        ]
        self.assert_contract_failure(
            [unit_row("U1", positions)], duplicate_constraints, "duplicate constraint_id"
        )

        invalid_pattern_count = constraint_row("U1", "C1", (10,))
        invalid_pattern_count["pattern_count"] = "2"
        self.assert_contract_failure(
            [unit_row("U1", positions)],
            [invalid_pattern_count],
            "pattern_count must equal 1",
        )

    def test_rejects_duplicate_component_and_route_site_membership(self) -> None:
        positions = (10, 20)
        first = unit_row("U1", positions)
        duplicate_component = unit_row("U2", (30, 40), ps="200")
        duplicate_component["component_id"] = first["component_id"]
        self.assert_contract_failure(
            [first, duplicate_component], [], "duplicate component_id"
        )

        ambiguous_route = unit_row("U2", positions)
        self.assert_contract_failure(
            [first, ambiguous_route], [], "duplicates an exact PS x HP route-site membership"
        )

    def test_k13_split_uses_retained_weight_before_other_ties(self) -> None:
        positions = tuple(range(100, 1400, 100))
        constraints = [
            constraint_row("U1", "LEFT", positions[:6], weight="10"),
            constraint_row("U1", "RIGHT", positions[6:], weight="10"),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            process, output_dir = self.run_kernel(
                pathlib.Path(temporary), [unit_row("U1", positions)], constraints
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            blocks = read_tsv(output_dir / "blocks.tsv")
            self.assertEqual([row["positions1"] for row in blocks], [
                ",".join(map(str, positions[:6])),
                ",".join(map(str, positions[6:])),
            ])
            self.assertEqual({row["unit_cut_indices_zero_based"] for row in blocks}, {"6"})
            dispositions = read_tsv(output_dir / "dispositions.tsv")
            self.assertEqual(
                [row["disposition"] for row in dispositions], ["retained", "retained"]
            )

    def test_multicut_constraint_is_classified_once(self) -> None:
        positions = tuple(range(1, 26))
        huge = "100000000000000000000000000000000000000000000000001"
        with tempfile.TemporaryDirectory() as temporary:
            process, output_dir = self.run_kernel(
                pathlib.Path(temporary),
                [unit_row("U1", positions)],
                [constraint_row("U1", "LONG", positions, weight=huge)],
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            dispositions = read_tsv(output_dir / "dispositions.tsv")
            self.assertEqual(len(dispositions), 1)
            self.assertEqual(
                dispositions[0]["disposition"], "unavoidable_span_gt_max_block_size"
            )
            self.assertEqual(dispositions[0]["crossed_cut_count"], "2")
            self.assertEqual(dispositions[0]["molecule_weight"], huge)
            summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["counts"]["constraints_total"], 1)
            self.assertEqual(summary["counts"]["constraints_unavoidable"], 1)
            self.assertEqual(summary["weights"]["total"], huge)
            self.assertEqual(summary["weights"]["lost"], huge)

    def test_weight_ten_and_public_indices_match_python_adapter_contract(self) -> None:
        positions = (10, 20)
        with tempfile.TemporaryDirectory() as temporary:
            process, output_dir = self.run_kernel(
                pathlib.Path(temporary),
                [unit_row("U1", positions)],
                [constraint_row("U1", "C1", positions, weight="10")],
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            blocks = read_tsv(output_dir / "blocks.tsv")
            self.assertEqual(blocks[0]["block_id"], "U1:B0001")
            self.assertEqual(blocks[0]["block_index"], "1")
            self.assertEqual(blocks[0]["retained_molecule_weight"], "10")
            membership = read_tsv(output_dir / "membership.tsv")
            self.assertEqual([row["unit_local_index"] for row in membership], ["1", "2"])
            self.assertEqual([row["block_index"] for row in membership], ["1", "1"])
            dispositions = read_tsv(output_dir / "dispositions.tsv")
            self.assertEqual(dispositions[0]["molecule_weight"], "10")
            self.assertEqual(dispositions[0]["retained_block_index"], "1")
            self.assertEqual(dispositions[0]["retained_block_id"], "U1:B0001")
            summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["weights"], {"total": "10", "retained": "10", "lost": "0"})

    def test_header_only_inputs_produce_valid_zero_result(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            process, output_dir = self.run_kernel(pathlib.Path(temporary), [], [])
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertEqual(read_tsv(output_dir / "blocks.tsv"), [])
            self.assertEqual(read_tsv(output_dir / "membership.tsv"), [])
            self.assertEqual(read_tsv(output_dir / "dispositions.tsv"), [])
            summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertTrue(summary["all_pass"])
            self.assertEqual(
                summary["counts"],
                {
                    "units": 0,
                    "sites": 0,
                    "constraints_total": 0,
                    "blocks": 0,
                    "cuts": 0,
                    "constraints_retained": 0,
                    "constraints_cut": 0,
                    "constraints_unavoidable": 0,
                },
            )

    def test_rejects_nonprimary_threshold(self) -> None:
        row = unit_row("U1", (10, 20))
        row["threshold"] = "2"
        self.assert_contract_failure([row], [], "threshold must equal 3")

    def test_row_order_does_not_change_any_output_byte(self) -> None:
        first_positions = tuple(range(10, 150, 10))
        second_positions = (1000, 1100, 1200)
        units = [
            unit_row("U2", second_positions, chrom="chr2", hp="2", ps="900"),
            unit_row("U1", first_positions),
        ]
        constraints = [
            constraint_row(
                "U2", "Z", second_positions[:2], chrom="chr2", hp="2", ps="900", weight="7"
            ),
            constraint_row("U1", "B", first_positions[7:], weight="3"),
            constraint_row("U1", "A", first_positions[:7], weight="3"),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            first, first_output = self.run_kernel(root, units, constraints, label="first")
            second, second_output = self.run_kernel(
                root, list(reversed(units)), list(reversed(constraints)), label="second"
            )
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(second.returncode, 0, second.stderr)
            for name in OUTPUT_NAMES:
                self.assertEqual(
                    (first_output / name).read_bytes(),
                    (second_output / name).read_bytes(),
                    name,
                )

    def test_seeded_units_match_python_reference_exactly(self) -> None:
        rng = random.Random(20260722)
        for max_block_size in (1, 2, 3, 5, 12):
            units = []
            constraints = []
            expected = {}
            for unit_offset in range(5):
                n_sites = 1 + ((max_block_size * 3 + unit_offset * 5) % 19)
                positions = []
                cursor = 1000 * (unit_offset + 1)
                for _ in range(n_sites):
                    cursor += rng.randint(1, 100)
                    positions.append(cursor)
                positions_tuple = tuple(positions)
                unit_id = f"K{max_block_size}_U{unit_offset}"
                hp = "1" if unit_offset % 2 == 0 else "2"
                ps = str(100 + unit_offset)
                units.append(unit_row(unit_id, positions_tuple, hp=hp, ps=ps))
                python_edges = []
                for edge_offset in range(1 + unit_offset * 2):
                    edge_size = rng.randint(1, n_sites)
                    selected = tuple(sorted(rng.sample(positions, edge_size)))
                    constraint_id = f"{unit_id}_C{edge_offset:02d}"
                    weight = str(rng.randint(0, 10**6))
                    constraints.append(
                        constraint_row(
                            unit_id,
                            constraint_id,
                            selected,
                            hp=hp,
                            ps=ps,
                            weight=weight,
                        )
                    )
                    python_edges.append(
                        Hyperedge(
                            constraint_id,
                            selected,
                            Decimal(weight),
                            1,
                        )
                    )
                expected[unit_id] = partition_ordered_hypergraph(
                    positions_tuple,
                    python_edges,
                    max_block_size=max_block_size,
                )

            with tempfile.TemporaryDirectory() as temporary:
                process, output_dir = self.run_kernel(
                    pathlib.Path(temporary),
                    list(reversed(units)),
                    list(reversed(constraints)),
                    max_block_size=max_block_size,
                )
                self.assertEqual(process.returncode, 0, process.stderr)
                block_rows = read_tsv(output_dir / "blocks.tsv")
                disposition_rows = read_tsv(output_dir / "dispositions.tsv")
                for unit_id, reference in expected.items():
                    observed_blocks = [row for row in block_rows if row["unit_id"] == unit_id]
                    self.assertEqual(
                        [tuple(map(int, row["positions1"].split(","))) for row in observed_blocks],
                        [block.positions for block in reference.blocks],
                    )
                    self.assertEqual(
                        tuple(map(int, observed_blocks[0]["unit_cut_indices_zero_based"].split(",")))
                        if observed_blocks[0]["unit_cut_indices_zero_based"]
                        else (),
                        reference.cut_indices,
                    )
                    self.assertEqual(
                        int(observed_blocks[0]["unit_retained_molecule_weight"]),
                        int(reference.retained_weight),
                    )
                    self.assertEqual(
                        int(observed_blocks[0]["unit_retained_pattern_count"]),
                        reference.retained_pattern_count,
                    )
                    observed_dispositions = {
                        row["constraint_id"]: (
                            row["disposition"],
                            int(row["span_sites"]),
                            int(row["crossed_cut_count"]),
                            None
                            if row["retained_block_index"] == ""
                            else int(row["retained_block_index"]) - 1,
                        )
                        for row in disposition_rows
                        if row["unit_id"] == unit_id
                    }
                    reference_dispositions = {
                        row.constraint_id: (
                            row.status,
                            row.span_sites,
                            row.crossed_cut_count,
                            row.block_index,
                        )
                        for row in reference.dispositions
                    }
                    self.assertEqual(observed_dispositions, reference_dispositions)

    def test_output_adapter_matches_python_producer_shared_fields(self) -> None:
        producer = load_python_producer()
        positions = tuple(range(100, 1400, 100))
        sites = tuple(
            producer.Site(index, "chr1", position, "C", "T")
            for index, position in enumerate(positions)
        )
        python_unit = producer.Unit("U1", "component:U1", "PS_HP1", "1", "100", 3, sites)
        python_constraints = (
            producer.Constraint("U1", "LEFT", "1", "100", positions[:6], "A" * 6, 10),
            producer.Constraint("U1", "RIGHT", "1", "100", positions[6:], "R" * 7, 10),
        )
        python_blocks, python_membership, python_dispositions = producer._partition_units(
            "HCC1395", "chr1", (python_unit,), python_constraints, 12
        )

        normalized_constraints = [
            constraint_row("U1", "LEFT", positions[:6], weight="10"),
            constraint_row("U1", "RIGHT", positions[6:], weight="10"),
        ]
        normalized_constraints[1]["call_codes"] = "R" * 7
        with tempfile.TemporaryDirectory() as temporary:
            process, output_dir = self.run_kernel(
                pathlib.Path(temporary),
                [unit_row("U1", positions)],
                normalized_constraints,
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            cpp_blocks = read_tsv(output_dir / "blocks.tsv")
            cpp_membership = read_tsv(output_dir / "membership.tsv")
            cpp_dispositions = read_tsv(output_dir / "dispositions.tsv")

        self.assertEqual(
            [{field: row[field] for field in producer.BLOCK_FIELDS} for row in cpp_blocks],
            [{field: str(row[field]) for field in producer.BLOCK_FIELDS} for row in python_blocks],
        )
        shared_membership = tuple(
            field for field in producer.MEMBERSHIP_FIELDS if field != "site_index"
        )
        self.assertEqual(
            [{field: row[field] for field in shared_membership} for row in cpp_membership],
            [
                {field: str(row[field]) for field in shared_membership}
                for row in python_membership
            ],
        )
        self.assertEqual(
            [
                {field: row[field] for field in producer.DISPOSITION_FIELDS}
                for row in cpp_dispositions
            ],
            [
                {field: str(row[field]) for field in producer.DISPOSITION_FIELDS}
                for row in python_dispositions
            ],
        )

    def test_tie_break_uses_largest_gap_then_lexicographically_smallest_cut(self) -> None:
        positions = (10, 20, 30, 1000, 1010, 1020, 1030)
        with tempfile.TemporaryDirectory() as temporary:
            process, output_dir = self.run_kernel(
                pathlib.Path(temporary), [unit_row("U1", positions)], [], max_block_size=6
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            blocks = read_tsv(output_dir / "blocks.tsv")
            self.assertEqual({row["unit_cut_indices_zero_based"] for row in blocks}, {"3"})

        uniform = tuple(range(1, 14))
        with tempfile.TemporaryDirectory() as temporary:
            process, output_dir = self.run_kernel(
                pathlib.Path(temporary), [unit_row("U1", uniform)], []
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            blocks = read_tsv(output_dir / "blocks.tsv")
            self.assertEqual({row["unit_cut_indices_zero_based"] for row in blocks}, {"1"})


if __name__ == "__main__":
    unittest.main()
