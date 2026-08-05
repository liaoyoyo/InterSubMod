from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = TOPIC_ROOT / "scripts" / "exact_ps_k12_partition.py"
DATASET = "SYNTHETIC"
CHROM = "chr1"

SITE_FIELDS = ("site_index", "chrom", "pos1", "ref", "alt")
MEMBERSHIP_FIELDS = (
    "dataset",
    "chrom",
    "linkage_basis",
    "phase_set",
    "phase_set_status",
    "inference_role",
    "threshold",
    "site_index",
    "pos1",
    "component_id",
)
CALL_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "hp_family",
    "phase_set",
    "site_indices",
    "positions1",
    "call_codes",
)
UNITS_FIELDS = (
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
CONSTRAINT_FIELDS = (
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


def _write_tsv_gz(path: Path, fieldnames: tuple[str, ...], rows: list[dict]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv_gz(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


def _site_rows(count: int) -> list[dict]:
    return [
        {
            "site_index": index,
            "chrom": CHROM,
            "pos1": 100 + 100 * index,
            "ref": "C",
            "alt": "T",
        }
        for index in range(count)
    ]


def _membership_rows(
    component_id: str,
    phase_set: str,
    site_indices: range | tuple[int, ...],
    *,
    hp_family: str = "1",
) -> list[dict]:
    return [
        {
            "dataset": DATASET,
            "chrom": CHROM,
            "linkage_basis": f"PS_HP{hp_family}",
            "phase_set": phase_set,
            "phase_set_status": "KNOWN_PS_PRIMARY" if phase_set else "MISSING_PS_SENSITIVITY_ABSTAIN",
            "inference_role": "PRIMARY_PS_AWARE" if phase_set else "SENSITIVITY_ABSTAIN_NONPRIMARY",
            "threshold": "3",
            "site_index": site_index,
            "pos1": 100 + 100 * site_index,
            "component_id": component_id,
        }
        for site_index in site_indices
    ]


def _call_row(
    molecule_id: str,
    phase_set: str,
    site_indices: tuple[int, ...],
    call_codes: str,
    *,
    hp_family: str = "1",
) -> dict:
    return {
        "dataset": DATASET,
        "chrom": CHROM,
        "molecule_id": molecule_id,
        "hp_family": hp_family,
        "phase_set": phase_set,
        "site_indices": ",".join(map(str, site_indices)),
        "positions1": ",".join(str(100 + 100 * index) for index in site_indices),
        "call_codes": call_codes,
    }


def _run(
    tmp_path: Path,
    label: str,
    sites: list[dict],
    memberships: list[dict],
    calls: list[dict],
    *,
    threshold: int = 3,
    expect_success: bool = True,
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    input_dir = tmp_path / f"{label}.input"
    output_dir = tmp_path / f"{label}.output"
    input_dir.mkdir()
    site_path = input_dir / "site_catalog.tsv.gz"
    membership_path = input_dir / "site_component_membership.tsv.gz"
    calls_path = input_dir / "molecule_sparse_calls.tsv.gz"
    _write_tsv_gz(site_path, SITE_FIELDS, sites)
    _write_tsv_gz(membership_path, MEMBERSHIP_FIELDS, memberships)
    _write_tsv_gz(calls_path, CALL_FIELDS, calls)

    command = [
        sys.executable,
        str(SCRIPT),
        "--dataset",
        DATASET,
        "--chrom",
        CHROM,
        "--site-catalog",
        str(site_path),
        "--site-component-membership",
        str(membership_path),
        "--molecule-calls",
        str(calls_path),
        "--output-dir",
        str(output_dir),
        "--threshold",
        str(threshold),
        "--max-block-size",
        "12",
    ]
    completed = subprocess.run(
        command,
        cwd=TOPIC_ROOT.parents[1],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    if expect_success:
        assert completed.returncode == 0, (
            f"command failed: {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    else:
        assert completed.returncode != 0
    return output_dir, completed


def test_same_hp_different_ps_units_are_never_merged(tmp_path: Path) -> None:
    sites = _site_rows(2)
    memberships = _membership_rows("component_ps100", "100", (0, 1))
    memberships += _membership_rows("component_ps200", "200", (0, 1))
    calls = [
        _call_row("m100", "100", (0, 1), "AA"),
        _call_row("m200", "200", (0, 1), "RR"),
    ]

    output_dir, _ = _run(tmp_path, "two_ps", sites, memberships, calls)
    unit_fields, units = _read_tsv_gz(output_dir / "units.tsv.gz")
    constraint_fields, constraints = _read_tsv_gz(output_dir / "constraints.tsv.gz")

    assert tuple(unit_fields) == UNITS_FIELDS
    assert tuple(constraint_fields) == CONSTRAINT_FIELDS
    assert len(units) == 2
    assert {row["phase_set"] for row in units} == {"100", "200"}
    assert len({row["unit_id"] for row in units}) == 2
    assert {(row["phase_set"], row["molecule_weight"]) for row in constraints} == {
        ("100", "1"),
        ("200", "1"),
    }


def test_source_component_with_multiple_ps_fails_closed(tmp_path: Path) -> None:
    sites = _site_rows(2)
    memberships = _membership_rows("mixed_component", "100", (0,))
    memberships += _membership_rows("mixed_component", "200", (1,))

    _, completed = _run(
        tmp_path,
        "mixed_source",
        sites,
        memberships,
        [],
        expect_success=False,
    )

    assert "source component is not exact PS x HP" in completed.stderr


@pytest.mark.parametrize(
    ("changes", "case_label"),
    (
        ({"phase_set_status": "KNOWN_PS_SECONDARY"}, "wrong_status"),
        ({"linkage_basis": "PS_HP3"}, "wrong_basis"),
        ({"phase_set": ""}, "missing_phase_set"),
    ),
)
def test_inconsistent_primary_membership_fails_closed(
    tmp_path: Path, changes: dict[str, str], case_label: str
) -> None:
    sites = _site_rows(2)
    memberships = _membership_rows("component_ps100", "100", (0, 1))
    for row in memberships:
        row.update(changes)
        row["inference_role"] = "PRIMARY_PS_AWARE"

    _, completed = _run(
        tmp_path,
        f"primary_contract_{case_label}",
        sites,
        memberships,
        [],
        expect_success=False,
    )

    assert "PRIMARY_PS_AWARE membership violates exact primary contract" in completed.stderr


def test_primary_membership_at_other_threshold_is_audited_and_skipped(
    tmp_path: Path,
) -> None:
    sites = _site_rows(2)
    memberships = _membership_rows("component_threshold2", "100", (0, 1))
    for row in memberships:
        row["threshold"] = "2"

    output_dir, _ = _run(tmp_path, "other_threshold", sites, memberships, [])
    _, units = _read_tsv_gz(output_dir / "units.tsv.gz")
    receipt = json.loads((output_dir / "receipt.json").read_text())

    assert units == []
    assert receipt["counts"]["excluded_membership_other_threshold"] == 2
    assert receipt["counts"]["excluded_membership_nonprimary_or_missing_ps"] == 0


def test_molecule_ps_mismatch_is_not_routed(tmp_path: Path) -> None:
    sites = _site_rows(2)
    memberships = _membership_rows("component_ps100", "100", (0, 1))
    calls = [_call_row("wrong_ps", "999", (0, 1), "AA")]

    output_dir, _ = _run(tmp_path, "ps_mismatch", sites, memberships, calls)
    _, constraints = _read_tsv_gz(output_dir / "constraints.tsv.gz")
    receipt = json.loads((output_dir / "receipt.json").read_text())

    assert constraints == []
    assert receipt["counts"]["molecule_rows_no_matching_unit"] == 1
    assert receipt["checks"]["cross_ps_zero"] is True


def test_missing_ps_membership_and_molecule_are_excluded(tmp_path: Path) -> None:
    sites = _site_rows(2)
    memberships = _membership_rows("missing_ps", "", (0, 1))
    calls = [_call_row("missing_ps_read", "", (0, 1), "AA")]

    output_dir, _ = _run(tmp_path, "missing_ps", sites, memberships, calls)
    _, units = _read_tsv_gz(output_dir / "units.tsv.gz")
    _, constraints = _read_tsv_gz(output_dir / "constraints.tsv.gz")
    receipt = json.loads((output_dir / "receipt.json").read_text())

    assert units == []
    assert constraints == []
    assert receipt["counts"]["excluded_membership_nonprimary_or_missing_ps"] == 2
    assert receipt["counts"]["excluded_molecule_missing_phase_set"] == 1


@pytest.mark.parametrize(
    ("positions1", "expected_error", "case_label"),
    (
        ("100", "site_indices/positions1/call_codes length mismatch", "length"),
        ("100,999", "positions1 mismatch", "coordinate"),
    ),
)
def test_molecule_positions_must_match_site_catalog(
    tmp_path: Path,
    positions1: str,
    expected_error: str,
    case_label: str,
) -> None:
    sites = _site_rows(2)
    memberships = _membership_rows("component_ps100", "100", (0, 1))
    call = _call_row("bad_positions", "100", (0, 1), "AA")
    call["positions1"] = positions1

    _, completed = _run(
        tmp_path,
        f"molecule_positions_{case_label}",
        sites,
        memberships,
        [call],
        expect_success=False,
    )

    assert expected_error in completed.stderr


def test_k13_unit_is_split_and_all_sites_are_conserved(tmp_path: Path) -> None:
    sites = _site_rows(13)
    memberships = _membership_rows("component_k13", "100", range(13))

    output_dir, _ = _run(tmp_path, "k13", sites, memberships, [])
    _, units = _read_tsv_gz(output_dir / "units.tsv.gz")
    _, blocks = _read_tsv_gz(output_dir / "blocks.tsv.gz")
    _, membership = _read_tsv_gz(output_dir / "membership.tsv.gz")
    receipt = json.loads((output_dir / "receipt.json").read_text())

    assert len(units) == 1
    assert units[0]["k"] == "13"
    assert len(blocks) == 2
    assert all(int(row["k"]) <= 12 for row in blocks)
    assert len(membership) == 13
    assert len({(row["unit_id"], row["site_index"]) for row in membership}) == 13
    assert receipt["checks"]["unit_sites_assigned_exactly_once"] is True


def test_constraint_crossing_multiple_cuts_is_counted_once(tmp_path: Path) -> None:
    sites = _site_rows(25)
    memberships = _membership_rows("component_k25", "100", range(25))
    calls = [_call_row("long", "100", (0, 24), "AA")]

    output_dir, _ = _run(tmp_path, "multi_cut", sites, memberships, calls)
    _, dispositions = _read_tsv_gz(output_dir / "dispositions.tsv.gz")
    receipt = json.loads((output_dir / "receipt.json").read_text())

    assert len(dispositions) == 1
    assert dispositions[0]["disposition"] == "unavoidable_span_gt_max_block_size"
    assert int(dispositions[0]["crossed_cut_count"]) == 2
    assert dispositions[0]["molecule_weight"] == "1"
    assert receipt["constraint_mass"] == {
        "total": "1",
        "retained": "0",
        "cut": "0",
        "unavoidable": "1",
    }
    assert receipt["checks"]["constraint_mass_conserved"] is True


def test_input_row_shuffle_preserves_all_normalized_outputs(tmp_path: Path) -> None:
    sites = _site_rows(14)
    memberships = _membership_rows("component_a", "100", range(13))
    memberships += _membership_rows("component_b", "200", (12, 13), hp_family="2")
    calls = [
        _call_row("left", "100", (0, 1, 2), "AAX"),
        _call_row("right", "100", (10, 11, 12), "RAA"),
        _call_row("hp2", "200", (12, 13), "RR", hp_family="2"),
    ]

    forward_dir, _ = _run(tmp_path, "forward", sites, memberships, calls)
    reverse_dir, _ = _run(
        tmp_path,
        "reverse",
        list(reversed(sites)),
        list(reversed(memberships)),
        list(reversed(calls)),
    )

    for filename in (
        "units.tsv.gz",
        "constraints.tsv.gz",
        "blocks.tsv.gz",
        "membership.tsv.gz",
        "dispositions.tsv.gz",
    ):
        assert (forward_dir / filename).read_bytes() == (reverse_dir / filename).read_bytes()
    forward_receipt = json.loads((forward_dir / "receipt.json").read_text())
    reverse_receipt = json.loads((reverse_dir / "receipt.json").read_text())
    assert forward_receipt["semantic_result_sha256"] == reverse_receipt["semantic_result_sha256"]


def test_threshold_other_than_three_is_rejected(tmp_path: Path) -> None:
    _, completed = _run(
        tmp_path,
        "threshold_two",
        _site_rows(2),
        _membership_rows("component", "100", (0, 1)),
        [],
        threshold=2,
        expect_success=False,
    )

    assert "threshold must be exactly 3" in completed.stderr
