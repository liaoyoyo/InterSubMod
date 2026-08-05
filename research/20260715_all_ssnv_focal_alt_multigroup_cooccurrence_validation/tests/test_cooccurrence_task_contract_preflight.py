from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "audit_cooccurrence_task_contract_preflight.py"
)
SPEC = importlib.util.spec_from_file_location("cooccurrence_task_contract_preflight", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_preflight_canonical_command_is_isolated_and_formal_import_is_rejected(
    tmp_path: Path,
) -> None:
    command = MODULE.ANALYZER.canonical_preflight_command(
        manifest=tmp_path / "manifest.json",
        assignments=tmp_path / "assignments.jsonl.gz",
        sites=tmp_path / "sites.tsv.gz",
        intersubmod_root=tmp_path / "runs",
        independent_m2_audit=tmp_path / "independent_m2.json",
        primary_artifact_audit_pre=tmp_path / "primary_pre.json",
        output=tmp_path / "preflight.json",
    )
    assert command[:3] == [sys.executable, "-I", str(SCRIPT.resolve())]
    with pytest.raises(ValueError, match="direct-CLI only"):
        MODULE.ANALYZER.resolve_release_command(
            argv=command[3:],
            script_path=SCRIPT,
            expected_command=command,
            source_authority={
                "authority_id": MODULE.SOURCE_AUTHORITY.AUTHORITY_ID,
            },
            role="preflight",
        )


def test_raw_constant_counts_include_group_limit_site_but_gate_counts_do_not(
    monkeypatch,
) -> None:
    tasks = [
        {
            "assignment": {"sample": "S", "chrom": "chr1", "pos": 1},
            "screen_row": {"ref": "A", "alt": "C", "kind": "evaluable"},
        },
        {
            "assignment": {"sample": "S", "chrom": "chr1", "pos": 2},
            "screen_row": {"ref": "G", "alt": "T", "kind": "group-limit"},
        },
    ]
    monkeypatch.setattr(
        MODULE.ANALYZER,
        "m2_categorical_level_counts",
        lambda *_: {"hp_exact": 1, "hp_family": 1, "strand": 2},
    )

    def gate(row, _levels):
        if row["kind"] == "group-limit":
            return {
                "eligible": False,
                "evaluable": False,
                "status": (
                    "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM:"
                    "observed=11:maximum=10"
                ),
                "constant_axes": [],
            }
        return {
            "eligible": True,
            "evaluable": True,
            "status": "ELIGIBLE_M2_RESIDUAL_UNEXPLAINED_AND_AXES_DETERMINATE",
            "constant_axes": ["hp_exact", "hp_family"],
        }

    monkeypatch.setattr(MODULE.ANALYZER, "m2_screen_eligibility", gate)
    observed = MODULE.summarize_tasks(tasks)

    assert observed["raw_assignment_constant_axis_site_counts"] == {
        "hp_exact": 2,
        "hp_family": 2,
    }
    assert observed["gate_evaluated_constant_axis_site_counts"] == {
        "hp_exact": 1,
        "hp_family": 1,
    }
    assert observed["group_limit_examples"] == [
        {
            "sample": "S",
            "chrom": "chr1",
            "pos": 2,
            "ref": "G",
            "alt": "T",
            "raw_categorical_level_counts": {
                "hp_exact": 1,
                "hp_family": 1,
                "strand": 2,
            },
            "gate_constant_axes": [],
        }
    ]


def test_runtime_statistical_api_probe_exercises_pearson_and_fisher() -> None:
    observed = MODULE.runtime_statistical_api_probe()

    assert observed["scipy_version"]
    assert observed["pearson_3x2"]["testable"] is True
    assert observed["pearson_3x2"]["analytic_test"] == "pearson_chi_square"
    assert observed["pearson_3x2"]["p_analytic"] == observed["chi2_tuple_index_1_p"]
    assert observed["library_marker_pearson_3x2"]["testable"] is True
    assert observed["library_marker_pearson_3x2"]["analytic_test"] == (
        "pearson_chi_square_descriptive_only"
    )
    assert observed["library_marker_pearson_3x2"]["p_analytic"] == (
        observed["marker_chi2_tuple_index_1_p"]
    )
    assert observed["fisher_2x2"]["testable"] is True
    assert observed["fisher_2x2"]["analytic_test"] == "fisher_exact_2x2"


def test_upstream_audit_binding_requires_current_assignments_and_sites(tmp_path) -> None:
    assignments = tmp_path / "assignments.jsonl.gz"
    sites = tmp_path / "sites.tsv.gz"
    assignments.write_bytes(b"assignments")
    sites.write_bytes(b"sites")
    payload = {
        "schema_name": "audit.schema",
        "schema_version": "1.0.0",
        "pass": True,
        "inputs": {
            "stable_assignments": MODULE.artifact(assignments),
            "site_results": MODULE.artifact(sites),
        },
    }

    MODULE.validate_upstream_audit_input_binding(
        payload,
        schema_name="audit.schema",
        schema_version="1.0.0",
        expected_input_roles={"stable_assignments", "site_results"},
        assignments=assignments,
        sites=sites,
    )

    assignments.write_bytes(b"drifted")
    with pytest.raises(RuntimeError, match="stable assignments identity drifted"):
        MODULE.validate_upstream_audit_input_binding(
            payload,
            schema_name="audit.schema",
            schema_version="1.0.0",
            expected_input_roles={"stable_assignments", "site_results"},
            assignments=assignments,
            sites=sites,
        )


def test_upstream_audit_binding_rejects_unexpected_input_roles(tmp_path) -> None:
    assignments = tmp_path / "assignments.jsonl.gz"
    sites = tmp_path / "sites.tsv.gz"
    assignments.write_bytes(b"assignments")
    sites.write_bytes(b"sites")
    payload = {
        "schema_name": "audit.schema",
        "schema_version": "1.0.0",
        "pass": True,
        "inputs": {
            "stable_assignments": MODULE.artifact(assignments),
            "site_results": MODULE.artifact(sites),
            "unexpected": {},
        },
    }
    with pytest.raises(RuntimeError, match="input role set drifted"):
        MODULE.validate_upstream_audit_input_binding(
            payload,
            schema_name="audit.schema",
            schema_version="1.0.0",
            expected_input_roles={"stable_assignments", "site_results"},
            assignments=assignments,
            sites=sites,
        )


def test_independent_m2_authority_rejects_same_basename_claim_contract(
    tmp_path: Path,
) -> None:
    fake_claim = tmp_path / "claim-contract-v5.md"
    fake_claim.write_text("same basename, wrong authority\n", encoding="utf-8")
    payload = {
        "inputs": {
            "claim_contract": MODULE.artifact(fake_claim),
        },
        "code": {
            "independent_recount": MODULE.artifact(
                MODULE.SOURCE_AUTHORITY.EXPECTED_SOURCE_PATHS[
                    "independent_m2_auditor"
                ]
            )
        },
    }
    with pytest.raises(RuntimeError, match="claim-contract-v5 identity drifted"):
        MODULE.validate_independent_m2_authority(payload)


def test_raw_identity_task_replays_same_recovery_and_alt_readset_contract(monkeypatch) -> None:
    focal = MODULE.ANALYZER.LIB.Variant("chr1", 101, "A", "C")
    recovered = [
        MODULE.ANALYZER.RecoveredRead(
            read_id="0",
            read_name="read-0",
            expected_focal_call="A",
            focal_call="A",
            latest_hp="1-1",
            latest_ps=10,
            strand="+",
            label="g1",
            partner_calls={},
            readset_identity="read-0|chr1|0|1200|+",
        )
    ]

    def recover(*args):
        audit = args[-1]
        projection = MODULE.ANALYZER.TAGS.projection_key(
            "read-0", "chr1", 0, 1200, 60, "+"
        )
        audit.update(
            {
                "equivalence_contract": MODULE.ANALYZER.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
                "allowed_differing_auxiliary_tags": ["RG"],
                "missing_projection_policy": MODULE.ANALYZER.RAW_IDENTITY_MISSING_POLICY,
                "conflicting_analysis_payload_policy": MODULE.ANALYZER.RAW_IDENTITY_CONFLICT_POLICY,
                "expected_projections": 1,
                "matched_records": 2,
                "duplicate_projections_collapsed": 1,
                "duplicate_extra_records_collapsed": 1,
                "exact_duplicate_projections_collapsed": 0,
                "rg_only_duplicate_projections_collapsed": 1,
                "duplicate_projection_sha256": MODULE.ANALYZER.projection_digest(
                    [projection]
                ),
                "alignment_payload_sha256": "1" * 64,
                "analysis_scope_identity_contract": (
                    MODULE.ANALYZER.ANALYSIS_SCOPE_IDENTITY_CONTRACT
                ),
                "recovered_payload_sha256": "2" * 64,
                "analysis_site_payload_sha256": "3" * 64,
                "duplicate_projection_examples": [list(projection)],
                "duplicate_projection_records": [
                    {
                        "projection": list(projection),
                        "record_count": 2,
                        "classification": "RG_ONLY_DUPLICATE",
                        "differing_auxiliary_tags": ["RG"],
                    }
                ],
            }
        )
        return recovered

    monkeypatch.setattr(MODULE.ANALYZER, "recover_site_reads", recover)
    monkeypatch.setattr(MODULE.ANALYZER, "get_vcf", lambda *_: object())
    monkeypatch.setattr(
        MODULE.ANALYZER,
        "discover_focal_and_partners",
        lambda *_: (focal, []),
    )
    result = MODULE.audit_raw_identity_task(
        {
            "assignment": {"sample": "S", "chrom": "chr1", "pos": 101},
            "screen_row": {
                "ref": "A",
                "alt": "C",
                "alt_readset_sha256": MODULE.ANALYZER.focal_alt_readset_digest(recovered),
            },
            "entry": {
                "sample": "S",
                "all_ssnv_vcf": {"path": "/unused.vcf.gz"},
                "all_ssnv_vcf_index": {"path": "/unused.vcf.gz.csi"},
            },
        }
    )

    assert result["raw_identity_expected_projections"] == 1
    assert result["raw_identity_matched_records"] == 2
    assert result["raw_identity_rg_only_duplicate_projections_collapsed"] == 1


def test_raw_identity_summary_reconciles_site_weighted_counts() -> None:
    duplicate_projection = MODULE.ANALYZER.TAGS.projection_key(
        "read", "chr1", 0, 1200, 60, "+"
    )
    rows = [
        {
            "sample": "S",
            "chrom": "chr1",
            "pos": 1,
            "n_all_focal_ref_alt_reads": 2,
            "raw_identity_expected_projections": 2,
            "raw_identity_matched_records": 3,
            "raw_identity_duplicate_projections_collapsed": 1,
            "raw_identity_duplicate_extra_records_collapsed": 1,
            "raw_identity_exact_duplicate_projections_collapsed": 0,
            "raw_identity_rg_only_duplicate_projections_collapsed": 1,
            "raw_identity_duplicate_projection_sha256": (
                MODULE.ANALYZER.projection_digest([duplicate_projection])
            ),
            "raw_identity_alignment_payload_sha256": "1" * 64,
            "raw_identity_recovered_payload_sha256": "2" * 64,
            "raw_identity_analysis_site_payload_sha256": "3" * 64,
            "duplicate_projection_examples": [list(duplicate_projection)],
            "duplicate_projection_records": [
                {
                    "projection": list(duplicate_projection),
                    "record_count": 2,
                    "classification": "RG_ONLY_DUPLICATE",
                    "differing_auxiliary_tags": ["RG"],
                }
            ],
        },
        {
            "sample": "S",
            "chrom": "chr1",
            "pos": 2,
            "n_all_focal_ref_alt_reads": 4,
            "raw_identity_expected_projections": 4,
            "raw_identity_matched_records": 4,
            "raw_identity_duplicate_projections_collapsed": 0,
            "raw_identity_duplicate_extra_records_collapsed": 0,
            "raw_identity_exact_duplicate_projections_collapsed": 0,
            "raw_identity_rg_only_duplicate_projections_collapsed": 0,
            "raw_identity_duplicate_projection_sha256": (
                MODULE.ANALYZER.projection_digest([])
            ),
            "raw_identity_alignment_payload_sha256": "4" * 64,
            "raw_identity_recovered_payload_sha256": "5" * 64,
            "raw_identity_analysis_site_payload_sha256": "6" * 64,
            "duplicate_projection_examples": [],
            "duplicate_projection_records": [],
        },
    ]
    observed = MODULE.summarize_raw_identity(rows)
    reversed_observed = MODULE.summarize_raw_identity(list(reversed(rows)))

    assert observed["aggregate"]["site_tasks"] == 2
    assert observed["aggregate"]["expected_projection_occurrences"] == 6
    assert observed["aggregate"]["matched_record_occurrences"] == 7
    assert observed["aggregate"]["sites_with_collapsed_duplicates"] == 1
    assert len(observed["site_weighted_audit_sha256"]) == 64
    assert reversed_observed["site_weighted_audit_sha256"] == observed[
        "site_weighted_audit_sha256"
    ]
    assert observed["all_result_rows_passed_invariant_validation"] is True
    assert observed["failure_counts_materialized"] is False


def test_raw_identity_summary_rejects_contradictory_worker_row() -> None:
    row = {
        "sample": "S",
        "chrom": "chr1",
        "pos": 1,
        "n_all_focal_ref_alt_reads": 2,
        "raw_identity_expected_projections": 2,
        "raw_identity_matched_records": 3,
        "raw_identity_duplicate_projections_collapsed": 0,
        "raw_identity_duplicate_extra_records_collapsed": 0,
        "raw_identity_exact_duplicate_projections_collapsed": 0,
        "raw_identity_rg_only_duplicate_projections_collapsed": 0,
        "raw_identity_duplicate_projection_sha256": (
            MODULE.ANALYZER.projection_digest([])
        ),
        "raw_identity_alignment_payload_sha256": "1" * 64,
        "raw_identity_recovered_payload_sha256": "2" * 64,
        "raw_identity_analysis_site_payload_sha256": "3" * 64,
        "duplicate_projection_examples": [],
        "duplicate_projection_records": [],
    }
    with pytest.raises(RuntimeError, match="matched-record"):
        MODULE.summarize_raw_identity([row])
