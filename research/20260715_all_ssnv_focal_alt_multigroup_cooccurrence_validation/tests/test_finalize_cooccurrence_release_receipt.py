from __future__ import annotations

import importlib.util
import csv
import gzip
import json
from pathlib import Path
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
SCRIPT = SCRIPT_DIR / "finalize_cooccurrence_release_receipt.py"
SPEC = importlib.util.spec_from_file_location("finalize_cooccurrence_release_receipt", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_current_release_paths_bind_v4_v8_and_v7_summary_hotfix() -> None:
    assert (
        MODULE.CANONICAL_PRIMARY_AUDIT.name
        == "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
    )
    assert (
        MODULE.CANONICAL_PREFLIGHT.name
        == "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
    )
    assert MODULE.CANONICAL_PRODUCER_DIR.name == (
        "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
        "source_locked_command_parity"
    )


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def make_release_inputs(tmp_path: Path, *, digest_match: bool) -> dict[str, Path]:
    output_dir = tmp_path / "cooccurrence"
    output_dir.mkdir()
    paths = {
        "preflight": tmp_path / "preflight.json",
        "producer_receipt": output_dir / "run_receipt.json",
        "summary": output_dir / "summary.json",
        "sites": output_dir / "methyl_ssnv_site_results.tsv.gz",
        "pairs": output_dir / "methyl_ssnv_pair_results.tsv.gz",
        "duplicates": output_dir / "raw_identity_duplicate_audit.tsv.gz",
        "oracle": output_dir / "oracle_cases.json",
        "runner_source": MODULE.RUNNER_SOURCE,
        "output": output_dir / "release_receipt.json",
    }
    for name in ("sites", "pairs", "duplicates"):
        paths[name].write_bytes(f"{name}\n".encode())
    write_json(paths["oracle"], {"fixture": True})

    runtime = {"contract": "fixture-runtime"}
    preflight_inputs: dict[str, object] = {}
    preflight_sources = {
        "preflight": MODULE.artifact(MODULE.PREFLIGHT_SOURCE),
        **{
            role: MODULE.artifact(path)
            for role, path in MODULE.CANONICAL_CODE_PATHS.items()
        },
    }
    preflight_digest = "a" * 64
    write_json(
        paths["preflight"],
        {
            "schema_name": "intersubmod.cooccurrence_task_contract_preflight",
            "schema_version": "2.0.0",
            "started_at_utc": "2026-07-15T00:00:00+00:00",
            "finished_at_utc": "2026-07-15T00:01:00+00:00",
            "task_type": "B_comprehensive_validation",
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
            "inputs": preflight_inputs,
            "input_lock": {
                "identity_before": preflight_inputs,
                "identity_after": preflight_inputs,
                "all_primary_inputs_unchanged": True,
            },
            "runtime_environment_lock": {
                "identity_before": runtime,
                "identity_after": runtime,
                "direct_runtime_unchanged": True,
            },
            "code": {"source_identity_after": preflight_sources},
            "checks": {key: True for key in MODULE.PREFLIGHT_CHECK_KEYS},
            "observed": {
                "raw_bam_identity_recovery": {
                    "analysis_scope_identity_contract": (
                        MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT
                    ),
                    "site_weighted_audit_sha256": preflight_digest,
                    "aggregate": {"site_tasks": 102_842},
                }
            },
            "pass": True,
        },
    )
    write_json(
        paths["summary"],
        {
            "schema_name": (
                "intersubmod.methyl_ssnv_cooccurrence_validation.summary"
            ),
            "schema_version": "2.0.0",
            "task_type": "B_comprehensive_validation",
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
            "raw_bam_identity_recovery_audit": {
                "analysis_scope_identity_contract": (
                    MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT
                ),
                "site_weighted_audit_sha256": (
                    preflight_digest if digest_match else "b" * 64
                ),
                "n_site_tasks": 102_842,
            },
            "pass": True,
        },
    )

    producer_inputs = {
        "preflight_receipt": MODULE.artifact(paths["preflight"]),
    }
    producer_code = {
        role: MODULE.artifact(path)
        for role, path in MODULE.CANONICAL_CODE_PATHS.items()
    }
    source_modes = {role: "0o444" for role in producer_code}
    producer_outputs = {
        "summary_json": MODULE.artifact(paths["summary"]),
        "site_tsv": MODULE.artifact(paths["sites"]),
        "pair_tsv": MODULE.artifact(paths["pairs"]),
        "raw_identity_duplicate_audit_tsv": MODULE.artifact(paths["duplicates"]),
        "case_json": MODULE.artifact(paths["oracle"]),
    }
    write_json(
        paths["producer_receipt"],
        {
            "schema_name": (
                "intersubmod.methyl_ssnv_cooccurrence_validation.run_receipt"
            ),
            "schema_version": "2.0.0",
            "started_at_utc": "2026-07-15T00:01:10+00:00",
            "finished_at_utc": "2026-07-15T00:02:00+00:00",
            "task_type": "B_comprehensive_validation",
            "scope": "all_manifest_samples",
            "full_scope": True,
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
            "release_status": "PRODUCER_PASS_PENDING_RUNNER_RECONCILIATION",
            "inputs": producer_inputs,
            "input_lock": {
                "identity_before": producer_inputs,
                "identity_after_compute": producer_inputs,
                "identity_after_output": producer_inputs,
                "all_primary_inputs_unchanged": True,
            },
            "runtime_environment_lock": {
                "identity_before": runtime,
                "identity_after_compute": runtime,
                "identity_after_output": runtime,
                "direct_runtime_unchanged": True,
            },
            "outputs": producer_outputs,
            "code": producer_code,
            "source_lock": {
                "source_identity_before": producer_code,
                "source_identity_after_compute": producer_code,
                "source_identity_after_output": producer_code,
                "source_modes_before": source_modes,
                "source_modes_after_compute": source_modes,
                "source_modes_after_output": source_modes,
                "all_sources_read_only_and_unchanged": True,
            },
            "counts": {"stable_sites": 102_842},
            "pass": True,
        },
    )
    return paths


def finalizer_argv(paths: dict[str, Path]) -> list[str]:
    return [
        "--preflight",
        str(paths["preflight"].resolve()),
        "--producer-receipt",
        str(paths["producer_receipt"].resolve()),
        "--summary",
        str(paths["summary"].resolve()),
        "--sites",
        str(paths["sites"].resolve()),
        "--pairs",
        str(paths["pairs"].resolve()),
        "--duplicates",
        str(paths["duplicates"].resolve()),
        "--oracle",
        str(paths["oracle"].resolve()),
        "--runner-source",
        str(paths["runner_source"].resolve()),
        "--output",
        str(paths["output"].resolve()),
    ]


def test_finalizer_canonical_command_uses_same_python_cache_contract(
    tmp_path: Path,
) -> None:
    paths = make_release_inputs(tmp_path, digest_match=True)
    args = MODULE.parse_args(finalizer_argv(paths))
    expected_cache = tmp_path / MODULE.ANALYZER.CANONICAL_PYTHON_CACHE_DIRNAME
    assert MODULE.canonical_command(args)[:5] == [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={expected_cache}",
        str(SCRIPT.resolve()),
    ]


def test_release_finalizer_rejects_old_synthetic_non_tsv_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = make_release_inputs(tmp_path, digest_match=True)
    monkeypatch.setattr(MODULE, "require_mode_0444", lambda *_: None)
    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY,
        "validate_release_source_authority",
        lambda *_: {"pass": True},
    )
    monkeypatch.setattr(MODULE, "CANONICAL_PREFLIGHT", paths["preflight"])
    monkeypatch.setattr(MODULE, "CANONICAL_PRODUCER_DIR", paths["producer_receipt"].parent)
    with pytest.raises((RuntimeError, gzip.BadGzipFile, EOFError)):
        MODULE.main(finalizer_argv(paths))
    assert not paths["output"].exists()


def test_release_finalizer_rejects_preflight_summary_digest_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = make_release_inputs(tmp_path, digest_match=False)
    monkeypatch.setattr(MODULE, "require_mode_0444", lambda *_: None)
    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY,
        "validate_release_source_authority",
        lambda *_: {"pass": True},
    )
    monkeypatch.setattr(MODULE, "CANONICAL_PREFLIGHT", paths["preflight"])
    monkeypatch.setattr(MODULE, "CANONICAL_PRODUCER_DIR", paths["producer_receipt"].parent)
    with pytest.raises((RuntimeError, gzip.BadGzipFile, EOFError)):
        MODULE.main(finalizer_argv(paths))
    assert not paths["output"].exists()


def make_site_row(pos: int, *, readset: str, representative: bool) -> dict[str, str]:
    row = {column: "" for column in MODULE.ANALYZER.SITE_COLUMNS}
    row.update(
        {
            "sample": "S",
            "biological_id": "S",
            "chrom": "chr1",
            "pos": str(pos),
            "ref": "A",
            "alt": "C",
            "m2_screen_eligible": "false",
            "n_all_focal_ref_alt_reads": "1",
            "raw_identity_expected_projections": "1",
            "raw_identity_matched_records": "1",
            "raw_identity_duplicate_projections_collapsed": "0",
            "raw_identity_duplicate_extra_records_collapsed": "0",
            "raw_identity_exact_duplicate_projections_collapsed": "0",
            "raw_identity_rg_only_duplicate_projections_collapsed": "0",
            "raw_identity_duplicate_projection_sha256": (
                MODULE.ANALYZER.projection_digest([])
            ),
            "raw_identity_alignment_payload_sha256": "1" * 64,
            "raw_identity_recovered_payload_sha256": "2" * 64,
            "raw_identity_analysis_site_payload_sha256": f"{pos:064x}"[-64:],
            "n_partner_markers": "0",
            "n_pair_rows_reconciled": "0",
            "pair_row_count_reconciliation_pass": "true",
            "n_endpoint_a_testable_markers": "0",
            "n_endpoint_a_exact_identifiable_markers": "0",
            "n_endpoint_a_exact_not_identifiable_markers": "0",
            "n_endpoint_a_conditional_permutable_markers": "0",
            "n_pair_bh_discoveries": "0",
            "pair_bh_discovery_positions": "[]",
            "n_pair_by_discoveries": "0",
            "pair_by_discovery_positions": "[]",
            "n_pair_by_confirmed": "0",
            "pair_by_confirmed_positions": "[]",
            "n_spatially_separated_pair_by_20bp": "0",
            "spatially_separated_pair_by_positions_20bp": "[]",
            "top_marker_positions": "[]",
            "n_top_marker_pair_by_confirmed": "0",
            "top_marker_pair_by_confirmed_positions": "[]",
            "joint_signature_complete_marker_effect_supported_positions": "[]",
            "n_same_complete_read_effect_supported_top_pair_by": "0",
            "same_complete_read_effect_supported_top_pair_by_positions": "[]",
            "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp": "0",
            "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp": "[]",
            "joint_signature_global_by_discovery": "false",
            "multi_marker_molecular_haplotype_base_candidate": "false",
            "n_endpoint_a_pre_candidates": "0",
            "n_confirmed_markers": "0",
            "confirmed_marker_positions": "[]",
            "n_independent_confirmed_markers_20bp": "0",
            "genetically_anchored_multi_marker_candidate": "false",
            "genetically_anchored_multi_marker_candidate_by_sensitivity": "false",
            "alt_readset_sha256": readset,
            "alt_readset_representative": "true" if representative else "false",
        }
    )
    return row


def write_tsv_gz(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_recompute_fixture(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    partial: bool = False,
) -> dict[str, Path]:
    rows = [
        make_site_row(100, readset="a" * 64, representative=True),
        make_site_row(200, readset="b" * 64, representative=True),
        make_site_row(201, readset="b" * 64, representative=False),
    ]
    paths = {
        "sites": tmp_path / "methyl_ssnv_site_results.tsv.gz",
        "pairs": tmp_path / "methyl_ssnv_pair_results.tsv.gz",
        "duplicates": tmp_path / "raw_identity_duplicate_audit.tsv.gz",
        "oracle": tmp_path / "oracle_cases.json",
    }
    write_tsv_gz(
        paths["sites"],
        MODULE.ANALYZER.SITE_COLUMNS,
        rows[:2] if partial else rows,
    )
    write_tsv_gz(paths["pairs"], MODULE.ANALYZER.PAIR_COLUMNS, [])
    write_tsv_gz(
        paths["duplicates"], MODULE.ANALYZER.RAW_IDENTITY_DUPLICATE_COLUMNS, []
    )
    monkeypatch.setattr(MODULE, "EXPECTED_STABLE_SITE_TASKS", 3)
    monkeypatch.setattr(
        MODULE.ANALYZER,
        "ORACLE_FOCALS",
        (
            {
                "case_id": "S_chr1_100",
                "sample": "S",
                "chrom": "chr1",
                "pos": 100,
                "expected_partner_positions": [],
            },
        ),
    )
    monkeypatch.setattr(
        MODULE.ANALYZER,
        "SHARED_READSET_ORACLE",
        {
            "case_id": "S_chr1_200_201_shared",
            "sample": "S",
            "chrom": "chr1",
            "positions": [200, 201],
        },
    )
    oracle = {
        "schema_name": f"{MODULE.ANALYZER.SCHEMA_NAME}.oracle_cases",
        "schema_version": MODULE.ANALYZER.SCHEMA_VERSION,
        "partner_window_bp": MODULE.ANALYZER.PAIR_WINDOW_BP,
        "focal_cases": [
            {
                "case_id": "S_chr1_100",
                "sample": "S",
                "focal": {"chrom": "chr1", "pos": 100},
                "present": True,
                "expected_partner_positions": [],
                "observed_partner_positions": [],
                "partner_window_oracle_pass": True,
                "site": {
                    "sample": "S",
                    "chrom": "chr1",
                    "pos": 100,
                    "raw_identity_analysis_site_payload_sha256": f"{100:064x}"[-64:],
                },
                "pairs": [],
            }
        ],
        "shared_readset_case": {
            "case_id": "S_chr1_200_201_shared",
            "positions": [200, 201],
            "present_positions": [200, 201],
            "same_alt_readset": True,
            "one_alt_readset_representative": True,
        },
    }
    write_json(paths["oracle"], oracle)
    return paths


def test_recompute_output_contract_accepts_complete_structured_fixture(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = make_recompute_fixture(tmp_path, monkeypatch)
    observed = MODULE.recompute_output_contract(
        sites_path=paths["sites"],
        pairs_path=paths["pairs"],
        duplicates_path=paths["duplicates"],
        oracle_path=paths["oracle"],
        expected_stable_sites=3,
    )
    assert observed["stable_sites"] == 3
    assert observed["focal_partner_pairs"] == 0
    assert observed["multi_marker_molecular_haplotype_base_candidates"] == 0
    assert observed["raw_identity_expected_projection_occurrences"] == 3
    assert len(observed["site_weighted_audit_sha256"]) == 64


def test_recompute_output_contract_rejects_partial_site_table(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = make_recompute_fixture(tmp_path, monkeypatch, partial=True)
    with pytest.raises(RuntimeError, match="Site row count is not expected scope"):
        MODULE.recompute_output_contract(
            sites_path=paths["sites"],
            pairs_path=paths["pairs"],
            duplicates_path=paths["duplicates"],
            oracle_path=paths["oracle"],
            expected_stable_sites=3,
        )


def test_recompute_output_contract_rejects_arbitrary_non_gzip_payload(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = make_recompute_fixture(tmp_path, monkeypatch)
    paths["sites"].write_bytes(b"arbitrary-nonempty-payload\n")
    with pytest.raises((gzip.BadGzipFile, EOFError, RuntimeError)):
        MODULE.recompute_output_contract(
            sites_path=paths["sites"],
            pairs_path=paths["pairs"],
            duplicates_path=paths["duplicates"],
            oracle_path=paths["oracle"],
        )
