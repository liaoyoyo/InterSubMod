from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "recover_all_ssnv_tumor_ref_controls.py"
)
SPEC = importlib.util.spec_from_file_location("recover_all_ssnv_tumor_ref_controls", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


SOURCE_FIELDS = ["sample", "chrom", "pos", "ref", "alt", "region_dir"]


def recovery_row(pos: int = 100) -> dict[str, str]:
    row = {
        "sample": "S1",
        "chrom": "chr1",
        "pos": str(pos),
        "ref": "A",
        "alt": "G",
        "region_dir": "/tmp/region",
    }
    row.update({field: "" for field in MODULE.B.CONTROL_FIELDS})
    row.update(
        {
            "primary_assignment_n_core_reads": "6",
            "primary_assignment_labels_sha256": "abc",
            "n_tumor_alt": "6",
            "n_tumor_ref": "6",
            "n_tumor_alt_matrix": "6",
            "n_tumor_ref_matrix": "6",
            "ref_status": "evaluable",
            "ref_evaluable": "true",
            "ref_n_after_peel": "6",
            "ref_n_peeled": "0",
            "ref_coarse_ng": "2",
            "ref_modal_fraction": "1.0",
            "ref_unstable": "false",
            "ref_stable_null_multigroup": "true",
            "ref_cluster_sizes": '{"1":3,"2":3}',
            "joint_status": "evaluable",
            "joint_evaluable": "true",
            "joint_n_after_peel": "12",
            "joint_n_peeled": "0",
            "joint_coarse_ng": "2",
            "joint_modal_fraction": "1.0",
            "joint_unstable": "false",
            "joint_stable_null_multigroup": "true",
            "joint_cluster_sizes": '{"1":6,"2":6}',
            "joint_allele_testable": "true",
            "joint_allele_v": "0.4",
            "joint_allele_p_perm": "0.01",
            "joint_allele_n": "12",
            "joint_allele_axis_aligned": "true",
            "background_control_interpretation": "REF_REPLICATION_WEAKENS_ALT_SPECIFICITY",
            "component_dedup_key": "old",
            "component_representative": "true",
            "readset_dedup_key": "old",
            "readset_representative": "true",
        }
    )
    return row


def task_for(row: dict[str, str]) -> dict:
    return {
        "site": {field: row[field] for field in SOURCE_FIELDS},
        "primary_assignment_n_core_reads": 6,
        "primary_assignment_labels_sha256": "abc",
    }


def write_prefix(path: Path, rows: list[dict[str, str]], *, invalid_tail: bool = False) -> None:
    fields = [*SOURCE_FIELDS, *MODULE.B.CONTROL_FIELDS]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if invalid_tail:
        with path.open("ab") as handle:
            handle.write(b"not-a-second-gzip-member")


def test_salvage_accepts_only_complete_canonical_rows_before_bad_tail(tmp_path: Path) -> None:
    row = recovery_row()
    prefix = tmp_path / "prefix.tsv.gz"
    write_prefix(prefix, [row], invalid_tail=True)
    accepted, audit = MODULE.salvage_canonical_prefix(
        prefix,
        [*SOURCE_FIELDS, *MODULE.B.CONTROL_FIELDS],
        [task_for(row), task_for(recovery_row(200))],
    )
    assert len(accepted) == 1
    assert audit["accepted_complete_rows"] == 1
    assert audit["truncated_gzip_eof_observed"] is True
    assert accepted[0]["ref_evaluable"] is True
    assert accepted[0]["joint_allele_axis_aligned"] is True


def test_salvage_rejects_noncanonical_prefix_order(tmp_path: Path) -> None:
    observed = recovery_row(200)
    expected = recovery_row(100)
    prefix = tmp_path / "wrong-order.tsv.gz"
    write_prefix(prefix, [observed])
    with pytest.raises(RuntimeError, match="not canonical"):
        MODULE.salvage_canonical_prefix(
            prefix,
            [*SOURCE_FIELDS, *MODULE.B.CONTROL_FIELDS],
            [task_for(expected), task_for(observed)],
        )


def test_recovery_dependencies_match_frozen_exact_equivalence_sources() -> None:
    observed = MODULE.verify_dependencies()
    assert observed["serial_base_analyzer"]["sha256"] == MODULE.EXPECTED_BASE_ANALYZER_SHA256
    assert observed["pinned_seed_parallel_phylo"]["sha256"] == MODULE.EXPECTED_PINNED_PHYLO_SHA256
    assert observed["focal_alt_cluster_lib"]["sha256"] == MODULE.EXPECTED_FOCAL_LIB_SHA256
