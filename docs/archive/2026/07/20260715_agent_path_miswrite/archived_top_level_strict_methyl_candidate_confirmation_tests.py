from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "scripts"
    / "run_strict_methyl_candidate_confirmation.py"
)
SPEC = importlib.util.spec_from_file_location("strict_methyl_partition_robustness", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_candidate_table(path: Path, selection_column: str) -> None:
    fields = [
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "region_dir",
        "alt_readset_sha256",
        "screen_contract",
        selection_column,
    ]
    rows = [
        {
            "sample": "S1",
            "chrom": "chr1",
            "pos": 100,
            "ref": "A",
            "alt": "C",
            "region_dir": "/unused/one",
            "alt_readset_sha256": "digest-one",
            "screen_contract": MODULE.SCREEN_CONTRACT,
            selection_column: "true",
        },
        {
            "sample": "S1",
            "chrom": "chr1",
            "pos": 200,
            "ref": "G",
            "alt": "T",
            "region_dir": "/unused/two",
            "alt_readset_sha256": "digest-two",
            "screen_contract": MODULE.SCREEN_CONTRACT,
            selection_column: "false",
        },
    ]
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_receipt_fixture(
    tmp_path: Path,
    *,
    selection_column: str = MODULE.FORMAL_SELECTION_COLUMN,
    family_size: int = 1,
    full_scope: bool = True,
) -> tuple[Path, Path, Path, MODULE.SelectionColumnContract, list[dict[str, str]]]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    candidate_table = tmp_path / "cooccurrence-sites.tsv"
    assignments = tmp_path / "screen-assignments.jsonl"
    summary_path = tmp_path / "cooccurrence-summary.json"
    receipt_path = tmp_path / "cooccurrence-receipt.json"
    write_candidate_table(candidate_table, selection_column)
    assignments.write_text(
        json.dumps(
            {
                "schema_name": MODULE.ASSIGNMENT_SCHEMA,
                "schema_version": MODULE.ASSIGNMENT_SCHEMA_VERSION,
                "screen_contract": MODULE.SCREEN_CONTRACT,
                "artifact_identity_contract": MODULE.ARTIFACT_IDENTITY_CONTRACT,
                "sample": "S1",
                "chrom": "chr1",
                "pos": 100,
            }
        )
        + "\n",
        encoding="ascii",
    )
    selection = MODULE.resolve_selection_column(MODULE.candidate_table_fields(candidate_table))
    _, selected, _ = MODULE.load_candidates(candidate_table, selection.resolved_column)
    family_field = (
        "n_genetically_anchored_multi_marker_candidates_by_sensitivity"
        if selection.legacy_fallback_used
        else "n_multi_marker_molecular_haplotype_base_candidates"
    )
    receipt_family_field = (
        "n_genetically_anchored_multi_marker_candidates_by_sensitivity"
        if selection.legacy_fallback_used
        else "multi_marker_molecular_haplotype_base_candidates"
    )
    summary = {
        "schema_name": MODULE.COOCCURRENCE_SUMMARY_SCHEMA,
        "schema_version": MODULE.COOCCURRENCE_FORMAL_SCHEMA_VERSION,
        "task_type": "B_comprehensive_validation" if full_scope else "A_partial_explicit_scope",
        "scope": "all_manifest_samples" if full_scope else "explicit_sample_subset",
        "partner_geometry_audit": {"full_scope": full_scope},
        "n_stable_sites": 2,
        "candidate_unique_key_count": 2,
        family_field: family_size,
    }
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    receipt = {
        "schema_name": MODULE.COOCCURRENCE_RECEIPT_SCHEMA,
        "schema_version": MODULE.COOCCURRENCE_FORMAL_SCHEMA_VERSION,
        "pass": True,
        "inputs": {"assignments": MODULE.artifact(assignments)},
        "outputs": {
            "site_tsv": MODULE.artifact(candidate_table),
            "summary_json": MODULE.artifact(summary_path),
        },
        "counts": {
            "stable_sites": 2,
            "candidate_unique_key_count": 2,
            receipt_family_field: family_size,
        },
    }
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    return candidate_table, assignments, receipt_path, selection, selected


def validate_fixture(
    candidate_table: Path,
    assignments: Path,
    receipt_path: Path,
    selection: MODULE.SelectionColumnContract,
    selected: list[dict[str, str]],
) -> dict[str, object]:
    keys = MODULE.candidate_table_keys(candidate_table)
    selected_keys = [MODULE.site_key(row, source="test selected") for row in selected]
    return MODULE.validate_cooccurrence_receipt(
        receipt_path,
        candidate_table=candidate_table,
        assignments_path=assignments,
        n_candidate_rows=len(keys),
        candidate_keys=keys,
        selected_keys=selected_keys,
        selection=selection,
    )


def make_flatten_result() -> dict[str, object]:
    site = MODULE.PreparedSite(
        input_row={
            "sample": "S1",
            "chrom": "chr1",
            "pos": "100",
            "ref": "A",
            "alt": "C",
        },
        source_row_number=2,
        sample="S1",
        chrom="chr1",
        pos=100,
        region_dir=Path("/unused"),
        distance=np.zeros((8, 8)),
        methylation=np.zeros((8, 4)),
        kept_ids=[str(index) for index in range(8)],
        screen_labels=["a"] * 4 + ["b"] * 4,
        n_reads_total=8,
        n_alt_raw=8,
        n_alt_after_peel=8,
    )
    metrics = {suffix: True for suffix in MODULE.MODE_OUTPUT_SUFFIXES}
    return {
        "site": site,
        "modes": {mode: dict(metrics) for mode in MODULE.NULL_MODES},
        "strict_combined_empirical_p_postselection_descriptive": 0.01,
        "strict_postselection_bh_q_descriptive": 0.02,
        "strict_postselection_by_q_descriptive": 0.03,
        "strict_cross_null_assignment_ari": 1.0,
        "strict_cross_null_exact_partition_gate": True,
        "strict_assignment_concordance_ari_min": 1.0,
    }


def test_formal_parameters_and_canonical_selection_are_both_required_for_pass() -> None:
    result = make_flatten_result()
    cooccurrence = {"full_scope_gate": True, "candidate_table_unique_key_count": 2}
    formal_selection = MODULE.SelectionColumnContract(
        MODULE.FORMAL_SELECTION_COLUMN, MODULE.FORMAL_SELECTION_COLUMN, False
    )
    formal = MODULE.flatten_result(
        result,
        formal_selection,
        cooccurrence,
        1,
        MODULE.ConfirmationConfig(permutations=999, seeds=10),
    )
    assert formal["strict_analysis_class"] == "FORMAL"
    assert formal["strict_null_robustness_pass"] is True
    assert formal["strict_methyl_partition_robustness_status"] == (
        "R1_STRICT_METHYL_PARTITION_ROBUSTNESS_PASS"
    )

    low_parameter = MODULE.flatten_result(
        result,
        formal_selection,
        cooccurrence,
        1,
        MODULE.ConfirmationConfig(permutations=7, seeds=3),
    )
    assert low_parameter["strict_analysis_class"] == "EXPLORATORY_ONLY"
    assert low_parameter["strict_null_robustness_pass"] is False
    assert low_parameter["strict_scientific_status"] == "EXPLORATORY_ONLY"
    assert "strict_formal_parameter_gate" in low_parameter["strict_failure_reason"]

    legacy_selection = MODULE.SelectionColumnContract(
        MODULE.FORMAL_SELECTION_COLUMN, MODULE.LEGACY_SELECTION_COLUMNS[0], True
    )
    legacy = MODULE.flatten_result(
        result,
        legacy_selection,
        cooccurrence,
        1,
        MODULE.ConfirmationConfig(permutations=999, seeds=10),
    )
    assert legacy["strict_selection_column"] == MODULE.FORMAL_SELECTION_COLUMN
    assert legacy["strict_selection_source_column"] == MODULE.LEGACY_SELECTION_COLUMNS[0]
    assert legacy["strict_analysis_class"] == "EXPLORATORY_ONLY"
    assert legacy["strict_null_robustness_pass"] is False


def test_selection_resolution_prefers_formal_and_marks_legacy_only_as_alias() -> None:
    both = [MODULE.LEGACY_SELECTION_COLUMNS[0], MODULE.FORMAL_SELECTION_COLUMN]
    resolved = MODULE.resolve_selection_column(both)
    assert resolved.resolved_column == MODULE.FORMAL_SELECTION_COLUMN
    assert resolved.legacy_fallback_used is False
    legacy = MODULE.resolve_selection_column([MODULE.LEGACY_SELECTION_COLUMNS[0]])
    assert legacy.formal_column == MODULE.FORMAL_SELECTION_COLUMN
    assert legacy.legacy_fallback_used is True


def test_formal_not_evaluable_is_not_counted_as_robustness_failure() -> None:
    result = make_flatten_result()
    result["modes"]["column"]["null_valid_gate"] = False
    selection = MODULE.SelectionColumnContract(
        MODULE.FORMAL_SELECTION_COLUMN, MODULE.FORMAL_SELECTION_COLUMN, False
    )
    row = MODULE.flatten_result(
        result,
        selection,
        {"full_scope_gate": True, "candidate_table_unique_key_count": 2},
        1,
        MODULE.ConfirmationConfig(),
    )
    assert row["strict_methyl_partition_robustness_evaluable"] is False
    assert row["strict_null_robustness_pass"] is False
    assert row["strict_methyl_partition_robustness_status"] == "NOT_EVALUABLE"
    assert row["strict_not_evaluable_reason"] == "strict_column_null_valid_gate"
    summary = MODULE.summarize_rows(
        [row],
        n_input_rows=1,
        n_assignment_records=1,
        selection=selection,
        cooccurrence_contract={"full_scope_gate": True},
        config=MODULE.ConfirmationConfig(),
    )
    assert summary["counts"]["n_methyl_partition_robustness_not_evaluable"] == 1
    assert summary["counts"]["n_null_robustness_fail_retained"] == 0


def test_cooccurrence_receipt_locks_full_scope_hashes_counts_and_family(tmp_path: Path) -> None:
    fixture = build_receipt_fixture(tmp_path)
    contract = validate_fixture(*fixture)
    assert contract["full_scope_gate"] is True
    assert contract["candidate_table_row_count"] == 2
    assert contract["candidate_table_unique_key_count"] == 2
    assert contract["selected_candidate_row_count"] == 1
    assert contract["selected_candidate_unique_key_count"] == 1
    assert contract["postselection_family_size"] == 1
    assert contract["screen_assignment_artifact_hash_lock"] is True


def test_cooccurrence_receipt_rejects_partial_scope(tmp_path: Path) -> None:
    fixture = build_receipt_fixture(tmp_path, full_scope=False)
    with pytest.raises(MODULE.ConfirmationContractError, match="not Task-B full scope"):
        validate_fixture(*fixture)


def test_cooccurrence_receipt_rejects_site_and_assignment_hash_drift(tmp_path: Path) -> None:
    site_fixture = build_receipt_fixture(tmp_path / "site")
    site_fixture[0].write_text(site_fixture[0].read_text() + "\n", encoding="ascii")
    with pytest.raises(MODULE.ConfirmationContractError, match="site_tsv artifact (size|sha256)"):
        validate_fixture(*site_fixture)

    assignment_fixture = build_receipt_fixture(tmp_path / "assignment")
    assignment_fixture[1].write_text("drift\n", encoding="ascii")
    with pytest.raises(
        MODULE.ConfirmationContractError, match="assignments artifact (size|sha256)"
    ):
        validate_fixture(*assignment_fixture)


def test_cooccurrence_receipt_rejects_candidate_family_mismatch(tmp_path: Path) -> None:
    fixture = build_receipt_fixture(tmp_path, family_size=2)
    with pytest.raises(MODULE.ConfirmationContractError, match="candidate family mismatch"):
        validate_fixture(*fixture)


def test_fixed_k_loo_output_names_and_narrative_do_not_claim_discovery() -> None:
    assert not any(suffix.startswith("loo_") for suffix in MODULE.MODE_OUTPUT_SUFFIXES)
    assert any(
        suffix.startswith("fixed_k_loo_partition_stability_")
        for suffix in MODULE.MODE_OUTPUT_SUFFIXES
    )
    selection = MODULE.SelectionColumnContract(
        MODULE.FORMAL_SELECTION_COLUMN, MODULE.FORMAL_SELECTION_COLUMN, False
    )
    row = MODULE.flatten_result(
        make_flatten_result(),
        selection,
        {"full_scope_gate": True, "candidate_table_unique_key_count": 2},
        1,
        MODULE.ConfirmationConfig(),
    )
    summary = MODULE.summarize_rows(
        [row],
        n_input_rows=2,
        n_assignment_records=1,
        selection=selection,
        cooccurrence_contract={"full_scope_gate": True},
        config=MODULE.ConfirmationConfig(),
    )
    narrative = summary["parameters"]["fixed_k_loo_partition_stability"]
    assert "keeping the already observed K fixed" in narrative
    assert "does not rerun discovery" in narrative


def test_main_retains_small_parameter_runs_as_exploratory_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    candidate_table, assignments, receipt_path, _, _ = build_receipt_fixture(tmp_path)

    def fake_prepare(candidate: dict[str, str], _assignment: dict[str, object]):
        template = make_flatten_result()["site"]
        template.input_row = {
            key: value for key, value in candidate.items() if not key.startswith("__")
        }
        template.source_row_number = int(candidate["__source_row_number"])
        return template

    def fake_analyze(site, _config):
        result = make_flatten_result()
        result["site"] = site
        return result

    monkeypatch.setattr(MODULE, "prepare_site", fake_prepare)
    monkeypatch.setattr(MODULE, "analyze_selected_site", fake_analyze)
    output_dir = tmp_path / "strict-output"
    assert (
        MODULE.main(
            [
                "--candidate-table",
                str(candidate_table),
                "--assignments",
                str(assignments),
                "--cooccurrence-receipt",
                str(receipt_path),
                "--output-dir",
                str(output_dir),
                "--permutations",
                "7",
                "--seeds",
                "3",
            ]
        )
        == 0
    )
    with gzip.open(
        output_dir / "strict_methyl_candidate_confirmation_sites.tsv.gz",
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["strict_analysis_class"] == "EXPLORATORY_ONLY"
    assert row["strict_null_robustness_pass"] == "False"
    assert row["strict_scientific_status"] == "EXPLORATORY_ONLY"
    summary = json.loads(
        (output_dir / "strict_methyl_candidate_confirmation_summary.json").read_text()
    )
    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert summary["scientific_output_class"] == "EXPLORATORY_ONLY"
    assert summary["counts"]["n_null_robustness_pass"] == 0
    assert manifest["contracts"]["formal_parameter_gate"] is False
    assert manifest["cooccurrence_receipt_contract"]["postselection_family_size"] == 1


def test_cli_requires_cooccurrence_receipt() -> None:
    with pytest.raises(SystemExit):
        MODULE.build_parser().parse_args(
            [
                "--candidate-table",
                "sites.tsv",
                "--assignments",
                "assignments.jsonl",
                "--output-dir",
                "out",
            ]
        )
