#!/usr/bin/env python3
"""Tests for the Python/C++ exact-PS partition comparison receipt."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = TOPIC_ROOT / "scripts" / "compare_python_cpp_partitions.py"
PYTHON_PRODUCER = TOPIC_ROOT / "scripts" / "exact_ps_k12_partition.py"
CPP_SOURCE = TOPIC_ROOT / "cpp" / "exact_ps_k12_partition.cpp"

UNITS_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "k", "positions1",
)
CONSTRAINT_FIELDS = (
    "dataset", "chrom", "unit_id", "constraint_id", "hp_family",
    "phase_set", "positions1", "call_codes", "molecule_weight", "pattern_count",
)
PYTHON_BLOCK_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "block_id", "block_index",
    "start1", "end1", "k", "positions1", "retained_molecule_weight",
    "retained_pattern_count",
)
PYTHON_MEMBERSHIP_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "site_index", "pos1",
    "unit_local_index", "block_id", "block_index",
)
PYTHON_DISPOSITION_FIELDS = (
    "dataset", "chrom", "unit_id", "constraint_id", "hp_family", "phase_set",
    "positions1", "call_codes", "molecule_weight", "pattern_count", "disposition",
    "span_sites", "crossed_cut_count", "retained_block_index",
)
CPP_BLOCK_FIELDS = PYTHON_BLOCK_FIELDS + (
    "start_index_zero_based", "end_index_exclusive_zero_based",
    "retained_constraint_count", "unit_cut_indices_zero_based", "unit_cut_gaps",
    "unit_cut_gap_sum", "unit_retained_molecule_weight",
    "unit_retained_pattern_count", "unit_total_molecule_weight",
    "unit_total_pattern_count",
)
CPP_MEMBERSHIP_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "unit_local_index", "pos1",
    "block_id", "block_index", "site_index_in_unit_zero_based",
    "site_index_in_block_zero_based",
)
CPP_DISPOSITION_FIELDS = PYTHON_DISPOSITION_FIELDS + ("retained_block_id",)
SOURCE_SITE_FIELDS = ("site_index", "chrom", "pos1", "ref", "alt")
SOURCE_MEMBERSHIP_FIELDS = (
    "dataset", "chrom", "linkage_basis", "phase_set", "phase_set_status",
    "inference_role", "threshold", "site_index", "pos1", "component_id",
)
SOURCE_CALL_FIELDS = (
    "dataset", "chrom", "molecule_id", "hp_family", "phase_set",
    "site_indices", "positions1", "call_codes",
)


def tsv_bytes(fields: tuple[str, ...], rows: list[dict[str, str]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def write_gzip(path: Path, payload: bytes) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(payload)


def identity(path: Path) -> dict[str, object]:
    payload = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def write_plain_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> bytes:
    payload = tsv_bytes(fields, rows)
    path.write_bytes(payload)
    return payload


def rewrite_tsv_cell(path: Path, key_field: str, key: str, field: str, value: str) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        rows = list(reader)
    for row in rows:
        if row[key_field] == key:
            row[field] = value
            break
    else:
        raise AssertionError(f"missing fixture key {key!r}")
    write_plain_tsv(path, fields, rows)


def duplicate_first_tsv_row(path: Path) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        rows = list(reader)
    if not rows:
        raise AssertionError("fixture must contain at least one data row")
    write_plain_tsv(path, fields, rows + [dict(rows[0])])


def build_fixture(root: Path) -> tuple[Path, Path, Path, Path, Path]:
    python_dir = root / "python"
    cpp_dir = root / "cpp"
    python_dir.mkdir()
    cpp_dir.mkdir()

    units = [{
        "dataset": "HCC1395", "chrom": "chr1", "unit_id": "U1",
        "component_id": "component:U1", "linkage_basis": "PS_HP1",
        "hp_family": "1", "phase_set": "100", "threshold": "3",
        "k": "2", "positions1": "10,20",
    }]
    constraints = [{
        "dataset": "HCC1395", "chrom": "chr1", "unit_id": "U1",
        "constraint_id": "C1", "hp_family": "1", "phase_set": "100",
        "positions1": "10,20", "call_codes": "AA", "molecule_weight": "10",
        "pattern_count": "1",
    }]
    blocks = [{
        "dataset": "HCC1395", "chrom": "chr1", "unit_id": "U1",
        "component_id": "component:U1", "linkage_basis": "PS_HP1",
        "hp_family": "1", "phase_set": "100", "threshold": "3",
        "block_id": "U1:B0001", "block_index": "1", "start1": "10",
        "end1": "20", "k": "2", "positions1": "10,20",
        "retained_molecule_weight": "10", "retained_pattern_count": "1",
    }]
    membership = [
        {
            "dataset": "HCC1395", "chrom": "chr1", "unit_id": "U1",
            "component_id": "component:U1", "linkage_basis": "PS_HP1",
            "hp_family": "1", "phase_set": "100", "threshold": "3",
            "site_index": str(40 + offset), "pos1": str(position),
            "unit_local_index": str(offset + 1), "block_id": "U1:B0001",
            "block_index": "1",
        }
        for offset, position in enumerate((10, 20))
    ]
    dispositions = [{
        "dataset": "HCC1395", "chrom": "chr1", "unit_id": "U1",
        "constraint_id": "C1", "hp_family": "1", "phase_set": "100",
        "positions1": "10,20", "call_codes": "AA", "molecule_weight": "10",
        "pattern_count": "1", "disposition": "retained", "span_sites": "2",
        "crossed_cut_count": "0", "retained_block_index": "1",
    }]

    normalized = {
        "units.tsv.gz": (UNITS_FIELDS, units),
        "constraints.tsv.gz": (CONSTRAINT_FIELDS, constraints),
        "blocks.tsv.gz": (PYTHON_BLOCK_FIELDS, blocks),
        "membership.tsv.gz": (PYTHON_MEMBERSHIP_FIELDS, membership),
        "dispositions.tsv.gz": (PYTHON_DISPOSITION_FIELDS, dispositions),
    }
    for name, (fields, rows) in normalized.items():
        write_gzip(python_dir / name, tsv_bytes(fields, rows))

    python_receipt = {
        "schema_name": "intersubmod.exact_ps_k12_partition",
        "schema_version": "0.1.0",
        "all_pass": True,
        "parameters": {
            "threshold": 3,
            "max_block_size": 12,
            "accepted_inference_role": "PRIMARY_PS_AWARE",
            "accepted_linkage_basis": ["PS_HP1", "PS_HP2"],
            "phase_set_required": True,
            "fixed_linkage_calls": ["A", "R"],
            "non_linking_calls": ["D", "L", "O", "S", "X"],
        },
        "checks": {
            "unit_sites_assigned_exactly_once": True,
            "block_k_at_most_12": True,
            "k_at_most_12_has_one_block": True,
            "cross_ps_zero": True,
            "cross_hp_zero": True,
            "constraint_ids_disposed_exactly_once": True,
            "constraint_mass_conserved": True,
        },
        "counts": {
            "eligible_units": 1,
            "eligible_unit_sites": 2,
            "constraints": 1,
            "blocks": 1,
            "retained_constraints": 1,
            "cut_constraints": 0,
            "unavoidable_constraints": 0,
        },
        "constraint_mass": {
            "total": "10", "retained": "10", "cut": "0", "unavoidable": "0",
        },
        "outputs": {
            name: identity(python_dir / name)
            for name in normalized
        },
    }
    (python_dir / "receipt.json").write_text(
        json.dumps(python_receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )

    units_plain = root / "units.tsv"
    constraints_plain = root / "constraints.tsv"
    units_plain.write_bytes(gzip.decompress((python_dir / "units.tsv.gz").read_bytes()))
    constraints_plain.write_bytes(
        gzip.decompress((python_dir / "constraints.tsv.gz").read_bytes())
    )

    cpp_blocks = [{
        **blocks[0],
        "start_index_zero_based": "0", "end_index_exclusive_zero_based": "2",
        "retained_constraint_count": "1", "unit_cut_indices_zero_based": "",
        "unit_cut_gaps": "", "unit_cut_gap_sum": "0",
        "unit_retained_molecule_weight": "10", "unit_retained_pattern_count": "1",
        "unit_total_molecule_weight": "10", "unit_total_pattern_count": "1",
    }]
    cpp_membership = [
        {
            field: row[field]
            for field in CPP_MEMBERSHIP_FIELDS
            if field not in {"site_index_in_unit_zero_based", "site_index_in_block_zero_based"}
        }
        | {
            "site_index_in_unit_zero_based": str(offset),
            "site_index_in_block_zero_based": str(offset),
        }
        for offset, row in enumerate(membership)
    ]
    cpp_dispositions = [{**dispositions[0], "retained_block_id": "U1:B0001"}]
    write_plain_tsv(cpp_dir / "blocks.tsv", CPP_BLOCK_FIELDS, cpp_blocks)
    write_plain_tsv(cpp_dir / "membership.tsv", CPP_MEMBERSHIP_FIELDS, cpp_membership)
    write_plain_tsv(
        cpp_dir / "dispositions.tsv", CPP_DISPOSITION_FIELDS, cpp_dispositions
    )
    cpp_summary = {
        "schema_name": "intersubmod.exact_ps_k12_cpp_partition_summary",
        "schema_version": "1.0.0",
        "all_pass": True,
        "parameters": {
            "max_block_size": 12,
            "partition_type": "contiguous_nonoverlapping",
            "molecule_weight_type": "arbitrary_precision_nonnegative_integer",
            "objective_order": [
                "max_retained_molecule_weight", "max_retained_pattern_count",
                "min_blocks", "max_cut_gap_sum",
                "lexicographically_smaller_cut_tuple",
            ],
        },
        "counts": {
            "units": 1, "sites": 2, "constraints_total": 1, "blocks": 1,
            "cuts": 0, "constraints_retained": 1, "constraints_cut": 0,
            "constraints_unavoidable": 0,
        },
        "weights": {"total": "10", "retained": "10", "lost": "0"},
        "pattern_counts": {"total": "1", "retained": "1", "lost": "0"},
        "checks": {
            "exact_nonmissing_ps_hp_units": True,
            "strictly_increasing_unique_positions": True,
            "constraint_positions_within_unit": True,
            "unique_unit_and_constraint_ids": True,
            "max_block_size_respected": True,
            "site_membership_conserved": True,
            "constraint_dispositions_mutually_exclusive_and_conserved": True,
            "objective_matches_dispositions": True,
        },
    }
    (cpp_dir / "summary.json").write_text(
        json.dumps(cpp_summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return python_dir, units_plain, constraints_plain, cpp_dir, root / "comparison.json"


def run_comparator(
    python_dir: Path,
    units_plain: Path,
    constraints_plain: Path,
    cpp_dir: Path,
    output: Path,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--python-dir",
            str(python_dir),
            "--cpp-input-units",
            str(units_plain),
            "--cpp-input-constraints",
            str(constraints_plain),
            "--cpp-dir",
            str(cpp_dir),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


class ComparePythonCppPartitionsTests(unittest.TestCase):
    def test_tiny_real_producer_cpp_comparator_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_dir = root / "input"
            input_dir.mkdir()
            source_tables = {
                "site_catalog.tsv.gz": (
                    SOURCE_SITE_FIELDS,
                    [
                        {"site_index": "0", "chrom": "chr1", "pos1": "10", "ref": "C", "alt": "T"},
                        {"site_index": "1", "chrom": "chr1", "pos1": "20", "ref": "G", "alt": "A"},
                    ],
                ),
                "site_component_membership.tsv.gz": (
                    SOURCE_MEMBERSHIP_FIELDS,
                    [
                        {
                            "dataset": "SYNTHETIC", "chrom": "chr1",
                            "linkage_basis": "PS_HP1", "phase_set": "100",
                            "phase_set_status": "KNOWN_PS_PRIMARY",
                            "inference_role": "PRIMARY_PS_AWARE", "threshold": "3",
                            "site_index": str(index), "pos1": str(position),
                            "component_id": "component_ps100",
                        }
                        for index, position in enumerate((10, 20))
                    ],
                ),
                "molecule_sparse_calls.tsv.gz": (
                    SOURCE_CALL_FIELDS,
                    [{
                        "dataset": "SYNTHETIC", "chrom": "chr1",
                        "molecule_id": "read-1", "hp_family": "1",
                        "phase_set": "100", "site_indices": "0,1",
                        "positions1": "10,20", "call_codes": "AA",
                    }],
                ),
            }
            for name, (fields, rows) in source_tables.items():
                write_gzip(input_dir / name, tsv_bytes(fields, rows))

            python_dir = root / "python-output"
            producer = subprocess.run(
                [
                    sys.executable, str(PYTHON_PRODUCER),
                    "--dataset", "SYNTHETIC", "--chrom", "chr1",
                    "--site-catalog", str(input_dir / "site_catalog.tsv.gz"),
                    "--site-component-membership",
                    str(input_dir / "site_component_membership.tsv.gz"),
                    "--molecule-calls", str(input_dir / "molecule_sparse_calls.tsv.gz"),
                    "--output-dir", str(python_dir), "--threshold", "3",
                    "--max-block-size", "12",
                ],
                text=True, capture_output=True, check=False, timeout=60,
            )
            self.assertEqual(producer.returncode, 0, producer.stderr)

            units_plain = root / "units.tsv"
            constraints_plain = root / "constraints.tsv"
            units_plain.write_bytes(
                gzip.decompress((python_dir / "units.tsv.gz").read_bytes())
            )
            constraints_plain.write_bytes(
                gzip.decompress((python_dir / "constraints.tsv.gz").read_bytes())
            )

            binary = root / "exact_ps_k12_partition"
            compile_process = subprocess.run(
                [
                    "g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
                    str(CPP_SOURCE), "-o", str(binary),
                ],
                text=True, capture_output=True, check=False, timeout=120,
            )
            self.assertEqual(compile_process.returncode, 0, compile_process.stderr)
            cpp_dir = root / "cpp-output"
            cpp_process = subprocess.run(
                [
                    str(binary), "--units", str(units_plain), "--constraints",
                    str(constraints_plain), "--output-dir", str(cpp_dir),
                    "--max-block-size", "12",
                ],
                text=True, capture_output=True, check=False, timeout=60,
            )
            self.assertEqual(cpp_process.returncode, 0, cpp_process.stderr)

            output = root / "comparison.json"
            completed = run_comparator(
                python_dir, units_plain, constraints_plain, cpp_dir, output
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            receipt = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(receipt["counts"]["aggregate"]["units"], 1)
            self.assertEqual(receipt["counts"]["aggregate"]["sites"], 2)

    def test_pass_and_no_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            completed = run_comparator(*args)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(receipt["mismatch_count"], 0)
            self.assertEqual(receipt["counts"]["aggregate"]["weights"]["total"], "10")
            before = args[-1].read_bytes()
            repeated = run_comparator(*args)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual(args[-1].read_bytes(), before)

    def test_tampered_cpp_block_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            rewrite_tsv_cell(args[3] / "blocks.tsv", "block_id", "U1:B0001", "start1", "11")
            completed = run_comparator(*args)
            self.assertNotEqual(completed.returncode, 0)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertFalse(receipt["all_pass"])
            self.assertTrue(any(row["category"] == "blocks" for row in receipt["mismatch_samples"]))

    def test_tampered_cpp_disposition_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            rewrite_tsv_cell(
                args[3] / "dispositions.tsv",
                "constraint_id",
                "C1",
                "disposition",
                "cut",
            )
            completed = run_comparator(*args)
            self.assertNotEqual(completed.returncode, 0)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertFalse(receipt["all_pass"])
            self.assertTrue(
                any(row["category"] == "dispositions" for row in receipt["mismatch_samples"])
            )

    def test_tampered_plain_input_bytes_fail_before_semantic_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            args[1].write_bytes(args[1].read_bytes() + b"\n")
            completed = run_comparator(*args)
            self.assertNotEqual(completed.returncode, 0)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertFalse(receipt["all_pass"])
            self.assertFalse(receipt["checks"]["normalized_input_payload_bytes_equal"])
            self.assertTrue(
                any(row["category"] == "input_payload" for row in receipt["mismatch_samples"])
            )

    def test_tampered_header_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            blocks_path = args[3] / "blocks.tsv"
            blocks_path.write_bytes(
                blocks_path.read_bytes().replace(b"\tstart1\t", b"\twrong_start1\t", 1)
            )
            completed = run_comparator(*args)
            self.assertNotEqual(completed.returncode, 0)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertFalse(receipt["all_pass"])
            self.assertFalse(receipt["checks"]["headers_exact"])
            self.assertTrue(
                any(row["category"] == "header" for row in receipt["mismatch_samples"])
            )

    def test_duplicate_key_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            duplicate_first_tsv_row(args[3] / "dispositions.tsv")
            completed = run_comparator(*args)
            self.assertNotEqual(completed.returncode, 0)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertFalse(receipt["all_pass"])
            self.assertFalse(receipt["checks"]["keys_unique"])
            self.assertTrue(
                any(row["category"] == "dispositions" for row in receipt["mismatch_samples"])
            )

    def test_tampered_cpp_parameter_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            summary_path = args[3] / "summary.json"
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["parameters"]["max_block_size"] = 11
            summary_path.write_text(
                json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
            )
            completed = run_comparator(*args)
            self.assertNotEqual(completed.returncode, 0)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertFalse(receipt["checks"]["cpp_summary_valid"])
            self.assertTrue(
                any(row["category"] == "cpp_summary" for row in receipt["mismatch_samples"])
            )

    def test_tampered_cpp_pattern_count_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = build_fixture(Path(temporary))
            summary_path = args[3] / "summary.json"
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["pattern_counts"]["retained"] = "9"
            summary_path.write_text(
                json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
            )
            completed = run_comparator(*args)
            self.assertNotEqual(completed.returncode, 0)
            receipt = json.loads(args[-1].read_text(encoding="utf-8"))
            self.assertFalse(receipt["checks"]["declared_aggregates_exact"])
            self.assertTrue(
                any(
                    row["category"] == "declared_aggregates"
                    for row in receipt["mismatch_samples"]
                )
            )


if __name__ == "__main__":
    unittest.main()
