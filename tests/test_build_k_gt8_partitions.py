import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPO_ROOT
    / "research/20260718_k_gt8_read_supported_segmentation"
    / "scripts/build_k_gt8_partitions.py"
)
DATASET = "SYNTHETIC"
CHROM = "chr1"

SITE_FIELDS = ("site_index", "chrom", "pos1", "ref", "alt")
CALL_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "hp_family",
    "phase_set",
    "site_indices",
    "call_codes",
)


def _write_tsv_gz(path, fieldnames, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _site_rows():
    # Component 1: the legacy 12-site case. Component 2: 24 sites so one
    # endpoint-to-endpoint constraint must cross two k=8 cut boundaries.
    positions = (
        tuple(range(100, 1_201, 100))
        + tuple(range(100_000, 102_301, 100))
        + (200_000,)
    )
    return [
        {
            "site_index": site_index,
            "chrom": CHROM,
            "pos1": pos1,
            "ref": "C",
            "alt": "T",
        }
        for site_index, pos1 in enumerate(positions)
    ]


def _call_row(molecule_id, hp_family, phase_set, indices, codes):
    return {
        "dataset": DATASET,
        "chrom": CHROM,
        "molecule_id": molecule_id,
        "hp_family": hp_family,
        "phase_set": phase_set,
        "site_indices": ",".join(map(str, indices)),
        "call_codes": codes,
    }


def _repeated_calls(prefix, count, hp_family, phase_set, indices, codes):
    return [
        _call_row(
            f"{prefix}_{ordinal:03d}",
            hp_family,
            phase_set,
            indices,
            codes,
        )
        for ordinal in range(count)
    ]


def _call_rows():
    rows = []
    rows.extend(
        _repeated_calls(
            "left",
            20,
            "1",
            "100",
            range(0, 6),
            "AAAAAA",
        )
    )
    rows.extend(
        _repeated_calls(
            "weak_bridge",
            3,
            "1",
            "100",
            (5, 6),
            "AR",
        )
    )
    rows.extend(
        _repeated_calls(
            "right",
            18,
            "1",
            "100",
            range(6, 12),
            "RRRAAA",
        )
    )
    rows.extend(
        _repeated_calls(
            "hp2_local",
            2,
            "2",
            "200",
            (0, 1),
            "RA",
        )
    )
    rows.extend(
        _repeated_calls(
            "long",
            5,
            "1",
            "300",
            (12, 35),
            "AA",
        )
    )

    # These rows exercise the production filters. HP3 and missing-PS rows
    # must not become constraints. The final two HP1 rows pass HP/PS filtering
    # but leave fewer than two fixed R/A calls and must also be excluded.
    rows.append(_call_row("hp3", "3", "100", (0, 1), "AA"))
    rows.append(_call_row("missing_ps", "1", "", (0, 1), "AA"))
    rows.append(_call_row("single_site", "1", "100", (2,), "A"))
    rows.append(_call_row("one_fixed_one_x", "1", "100", (2, 3), "AX"))
    return rows


def _run_cli(tmp_path, label, call_rows):
    site_catalog = tmp_path / f"{label}.sites.tsv.gz"
    molecule_calls = tmp_path / f"{label}.calls.tsv.gz"
    output_dir = tmp_path / f"{label}.output"
    _write_tsv_gz(site_catalog, SITE_FIELDS, _site_rows())
    _write_tsv_gz(molecule_calls, CALL_FIELDS, call_rows)

    command = [
        sys.executable,
        str(SCRIPT),
        "--dataset",
        DATASET,
        "--chrom",
        CHROM,
        "--site-catalog",
        str(site_catalog),
        "--molecule-calls",
        str(molecule_calls),
        "--output-dir",
        str(output_dir),
        "--expected-target-components",
        "2",
        "--expected-target-sites",
        "36",
    ]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, (
        f"command failed: {' '.join(command)}\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    receipt = json.loads((output_dir / "receipt.json").read_text())
    return output_dir, receipt


def _output_rows(output_dir, suffix):
    return _read_tsv_gz(output_dir / f"{DATASET}.{CHROM}.{suffix}.tsv.gz")


def test_cli_legacy_components_filters_mass_baseline_and_cut_accounting(
    tmp_path,
):
    output_dir, receipt = _run_cli(tmp_path, "forward", _call_rows())

    assert receipt["all_pass"] is True
    assert receipt["scope"] == {
        "dataset": DATASET,
        "chrom": CHROM,
        "site_catalog_sites": 37,
        "legacy_multilocus_components": 2,
        "legacy_singletons": 1,
    }
    assert all(receipt["checks"].values())

    counts = receipt["counts"]
    assert counts["target_components"] == 2
    assert counts["target_sites"] == 36
    assert counts["old_selected_sites"] == 16
    assert counts["old_cap_excluded_sites"] == 20
    assert counts["new_blocks"] == 5
    assert counts["exact_patterns"] == 5
    assert counts["raw_total_molecule_weight"] == 48
    assert counts["raw_retained_molecule_weight"] == 40
    assert counts["raw_lost_molecule_weight"] == 8
    assert counts["unavoidable_patterns"] == 1

    assert receipt["pattern_diagnostics"] == {
        "excluded_missing_ps_rows": 1,
        "excluded_nonprimary_hp_rows": 1,
        "molecule_rows_seen": 52,
        "primary_component_constraints_molecule_weight": 48,
        "primary_known_ps_rows": 50,
        "primary_single_target_site_observations": 2,
    }
    assert receipt["primary_rows_by_hp_family"] == {"1": 48, "2": 2}
    assert receipt["primary_hp_phase_set_unit_count"] == 3

    components = _output_rows(output_dir, "legacy_components")
    component_by_k = {int(row["pre_cap_k"]): row for row in components}
    legacy_12 = component_by_k[12]
    assert legacy_12["old_densest8_selected"] == "8"
    assert legacy_12["old_cap_excluded"] == "4"
    assert legacy_12["new_site_retained"] == "12"
    assert legacy_12["new_block_count"] == "2"
    assert legacy_12["raw_cut_indices"] == "6"
    assert legacy_12["raw_total_molecule_weight"] == "43"
    assert legacy_12["raw_retained_molecule_weight"] == "40"
    assert legacy_12["raw_lost_molecule_weight"] == "3"
    assert legacy_12["old_densest8_retained_molecule_weight"] == "25"
    assert legacy_12["old_densest8_retained_pattern_count"] == "3"

    legacy_24 = component_by_k[24]
    assert legacy_24["old_densest8_selected"] == "8"
    assert legacy_24["old_cap_excluded"] == "16"
    assert legacy_24["new_site_retained"] == "24"
    assert legacy_24["new_block_count"] == "3"
    assert legacy_24["raw_cut_indices"] == "8,16"
    assert legacy_24["raw_total_molecule_weight"] == "5"
    assert legacy_24["raw_retained_molecule_weight"] == "0"
    assert legacy_24["raw_lost_molecule_weight"] == "5"

    blocks = _output_rows(output_dir, "blocks")
    assert len(blocks) == 5
    assert all(int(row["k"]) <= 8 for row in blocks)
    blocks_by_component = {}
    for row in blocks:
        blocks_by_component.setdefault(row["legacy_component_id"], []).append(
            int(row["k"])
        )
    assert blocks_by_component[legacy_12["legacy_component_id"]] == [6, 6]
    assert blocks_by_component[legacy_24["legacy_component_id"]] == [8, 8, 8]

    membership = _output_rows(output_dir, "site_membership")
    assert len(membership) == 36
    assert sorted(int(row["site_index"]) for row in membership) == list(
        range(36)
    )
    disposition_counts = {}
    for row in membership:
        disposition_counts[row["old_densest8_disposition"]] = (
            disposition_counts.get(row["old_densest8_disposition"], 0) + 1
        )
    assert disposition_counts == {"selected": 16, "cap_excluded": 20}
    selected_legacy_12 = [
        int(row["pos1"])
        for row in membership
        if row["legacy_component_id"] == legacy_12["legacy_component_id"]
        and row["old_densest8_disposition"] == "selected"
    ]
    assert selected_legacy_12 == list(range(100, 801, 100))

    constraints = _output_rows(output_dir, "cut_constraints")
    assert len(constraints) == 5
    assert {
        (row["hp_family"], row["phase_set"]) for row in constraints
    } == {("1", "100"), ("1", "300"), ("2", "200")}
    assert all(row["hp_family"] in {"1", "2"} for row in constraints)
    assert all(row["phase_set"] for row in constraints)

    long_constraint = next(
        row for row in constraints if row["span_sites"] == "24"
    )
    assert long_constraint["molecule_weight"] == "5"
    assert (
        long_constraint["disposition"]
        == "unavoidable_span_gt_max_block_size"
    )
    assert long_constraint["crossed_cut_count"] == "2"
    assert long_constraint["retained_block_index"] == ""

    # The long constraint crosses two cuts but contributes exactly five once,
    # not five per cut. Together with the weak bridge it yields total loss 8.
    assert counts["raw_lost_molecule_weight"] == 5 + 3
    weak_bridge = next(
        row
        for row in constraints
        if row["positions"] == "600,700" and row["call_codes"] == "AR"
    )
    assert weak_bridge["molecule_weight"] == "3"
    assert weak_bridge["disposition"] == "cut"
    assert weak_bridge["crossed_cut_count"] == "1"


def test_cli_semantic_digest_is_deterministic_under_call_row_order(tmp_path):
    call_rows = _call_rows()
    forward_dir, forward_receipt = _run_cli(
        tmp_path,
        "forward",
        call_rows,
    )
    reverse_dir, reverse_receipt = _run_cli(
        tmp_path,
        "reverse",
        list(reversed(call_rows)),
    )

    assert (
        reverse_receipt["semantic_result_sha256"]
        == forward_receipt["semantic_result_sha256"]
    )
    assert reverse_receipt["counts"] == forward_receipt["counts"]
    assert reverse_receipt["status_counts"] == forward_receipt["status_counts"]

    for suffix in (
        "legacy_components",
        "blocks",
        "site_membership",
        "cut_constraints",
    ):
        assert _output_rows(reverse_dir, suffix) == _output_rows(
            forward_dir, suffix
        )
