from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "run_strict_methyl_candidate_confirmation.py"
)
SPEC = importlib.util.spec_from_file_location("run_strict_methyl_candidate_confirmation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


@pytest.fixture(autouse=True)
def test_only_release_source_authority(monkeypatch: pytest.MonkeyPatch) -> None:
    authority = {"authority_id": "TEST_ONLY_UNSIGNED_AUTHORITY", "pass": True}
    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY,
        "validate_release_source_authority",
        lambda *_args, **_kwargs: authority,
    )


SELECTION_COLUMN = MODULE.DEFAULT_SELECTION_COLUMN


def test_canonical_command_includes_completion_runner_python_cache_prefix(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "strict-output"
    command = MODULE.canonical_command(
        candidate_table=tmp_path / "candidates.tsv.gz",
        assignments=tmp_path / "assignments.jsonl.gz",
        cooccurrence_receipt=tmp_path / "run_receipt.json",
        cooccurrence_release_receipt=tmp_path / "release_receipt.json",
        output_dir=output_dir,
        config=MODULE.ConfirmationConfig(),
    )
    expected_cache = tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME
    assert command[:5] == [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={expected_cache}",
        str(SCRIPT),
    ]


def test_canonical_runtime_rejects_python_cache_prefix_drift(tmp_path: Path) -> None:
    output_dir = tmp_path / "strict-output"
    expected_cache = str(tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME)
    identity = {
        "interpreter_flags": {
            "isolated": 1,
            "no_user_site": 1,
            "safe_path": True,
        },
        "python_cache": {
            "configured_prefix": expected_cache,
            "xoption_prefix": expected_cache,
        },
        "environment": {
            "PYTHONPATH": None,
            "PYTHONNOUSERSITE": "1",
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "BLIS_NUM_THREADS": "1",
        },
    }
    MODULE.require_canonical_runtime(identity, output_dir=output_dir)

    identity["python_cache"]["xoption_prefix"] = str(tmp_path / "wrong-cache")
    with pytest.raises(
        MODULE.ConfirmationContractError,
        match="canonical python -I -X pycache_prefix",
    ):
        MODULE.require_canonical_runtime(identity, output_dir=output_dir)


def write_matrix(path: Path, row_ids: list[str], columns: list[str], matrix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["read_id", *columns])
        for read_id, values in zip(row_ids, matrix):
            writer.writerow(
                [
                    read_id,
                    *[
                        "NA" if not np.isfinite(value) else f"{float(value):.12f}"
                        for value in values
                    ],
                ]
            )


def make_region(root: Path) -> tuple[Path, list[str], list[str], list[dict[str, object]]]:
    region = root / "region"
    read_ids = [str(index) for index in range(8)]
    labels = ["screen-A"] * 4 + ["screen-B"] * 4
    low = np.asarray([0.03, 0.06, 0.10, 0.14, 0.18, 0.23, 0.28, 0.34, 0.39, 0.44])
    offsets = np.asarray([-0.015, -0.005, 0.008, 0.017])
    first = np.vstack([np.clip(low + offset, 0.001, 0.999) for offset in offsets])
    second = np.vstack(
        [np.clip(1.0 - low + offset, 0.001, 0.999) for offset in offsets]
    )
    methylation = np.vstack([first, second])
    distance = MODULE.F.bernoulli_distance(methylation)

    reads_path = region / "reads/reads.tsv"
    reads_path.parent.mkdir(parents=True, exist_ok=True)
    identities: list[dict[str, object]] = []
    with reads_path.open("w", encoding="ascii", newline="") as handle:
        fields = [
            "read_id",
            "read_name",
            "chr",
            "start",
            "end",
            "mapq",
            "hp",
            "alt_support",
            "is_tumor",
            "strand",
        ]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, (read_id, label) in enumerate(zip(read_ids, labels)):
            row = {
                "read_id": read_id,
                "read_name": f"read-{index}",
                "chr": "chr1",
                "start": 90 + index,
                "end": 160 + index,
                "mapq": 60,
                "hp": "1-1" if index < 4 else "1-2",
                "alt_support": "ALT",
                "is_tumor": "1",
                "strand": "+" if index % 2 == 0 else "-",
            }
            writer.writerow(row)
            identities.append(
                {
                    "read_id": read_id,
                    "read_name": row["read_name"],
                    "chrom": row["chr"],
                    "start": row["start"],
                    "end": row["end"],
                    "mapq": row["mapq"],
                    "strand": row["strand"],
                    "label": label,
                }
            )

    write_matrix(region / "distance/BERNOULLI/matrix.csv", read_ids, read_ids, distance)
    write_matrix(
        region / "methylation/methylation.csv",
        read_ids,
        [str(1000 + index) for index in range(methylation.shape[1])],
        methylation,
    )
    return region, read_ids, labels, identities


def write_candidate_table(path: Path, region: Path, include_unselected: bool = True) -> None:
    _, primary_reads = MODULE.load_reads_strict(region / "reads/reads.tsv")
    digest = MODULE.readset_digest(primary_reads.values())
    fields = [
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "region_dir",
        "alt_readset_sha256",
        "screen_contract",
        "n_alt_after_peel",
        "cooccurrence_fdr_gate",
        "cooccurrence_effect_gate",
        "independent_marker_gate",
        SELECTION_COLUMN,
    ]
    rows = [
        {
            "sample": "SYNTH",
            "chrom": "chr1",
            "pos": "100",
            "ref": "A",
            "alt": "C",
            "region_dir": str(region),
            "alt_readset_sha256": digest,
            "screen_contract": MODULE.SCREEN_CONTRACT,
            "n_alt_after_peel": "8",
            "cooccurrence_fdr_gate": "true",
            "cooccurrence_effect_gate": "true",
            "independent_marker_gate": "true",
            SELECTION_COLUMN: "true",
        }
    ]
    if include_unselected:
        rows.append(
            {
                "sample": "SYNTH",
                "chrom": "chr1",
                "pos": "200",
                "ref": "G",
                "alt": "T",
                "region_dir": str(region.parent / "deliberately-missing-region"),
                "alt_readset_sha256": "unselected-digest",
                "screen_contract": MODULE.SCREEN_CONTRACT,
                "n_alt_after_peel": "",
                "cooccurrence_fdr_gate": "false",
                "cooccurrence_effect_gate": "true",
                "independent_marker_gate": "true",
                SELECTION_COLUMN: "false",
            }
        )
    with gzip.open(path, "wt", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_assignment(
    region: Path,
    read_ids: list[str],
    labels: list[str],
    identities: list[dict[str, object]],
) -> dict[str, object]:
    return {
        "schema_name": "intersubmod.all_ssnv_stable_multigroup_read_assignment",
        "schema_version": "1.0.0",
        "screen_contract": MODULE.SCREEN_CONTRACT,
        "artifact_identity_contract": MODULE.ARTIFACT_IDENTITY_CONTRACT,
        "sample": "SYNTH",
        "chrom": "chr1",
        "pos": 100,
        "region_dir": str(region),
        "all_after_peel_read_ids": read_ids,
        "coarse_labels_all_after_peel": labels,
        "read_ids": read_ids,
        "read_names": [f"read-{index}" for index in range(len(read_ids))],
        "labels": labels,
        "core_reads": identities,
        "strict_confirm_candidate": True,
        "primary_artifacts": {
            "reads": MODULE.artifact(region / "reads/reads.tsv"),
            "distance_matrix": MODULE.artifact(region / "distance/BERNOULLI/matrix.csv"),
            "methylation_matrix": MODULE.artifact(region / "methylation/methylation.csv"),
        },
        "posthoc": {"ref": "A", "alt": "C"},
    }


def write_assignments(path: Path, records: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")


def write_cooccurrence_receipt(
    candidate_table: Path,
    assignments: Path,
    *,
    full_scope: bool = True,
    family_size: int | None = None,
) -> Path:
    source_authority = MODULE.SOURCE_AUTHORITY.validate_release_source_authority()
    fields = MODULE.candidate_table_fields(candidate_table)
    selection = MODULE.resolve_selection_column(fields)
    keys = MODULE.candidate_table_keys(candidate_table)
    _, selected, _ = MODULE.load_candidates(candidate_table, selection.resolved_column)
    selected_count = len(selected) if family_size is None else family_size
    summary_path = candidate_table.parent / "cooccurrence_summary.json"
    summary = {
        "schema_name": MODULE.COOCCURRENCE_SUMMARY_SCHEMA,
        "schema_version": MODULE.COOCCURRENCE_FORMAL_SCHEMA_VERSION,
        "task_type": (
            "B_comprehensive_validation" if full_scope else "A_partial_explicit_scope"
        ),
        "scope": "all_manifest_samples" if full_scope else "explicit_sample_subset",
        "partner_geometry_audit": {"full_scope": full_scope},
        "n_stable_sites": len(keys),
        "candidate_unique_key_count": len(set(keys)),
        "n_multi_marker_molecular_haplotype_base_candidates": selected_count,
        "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
        "raw_bam_identity_recovery_audit": {
            "equivalence_contract": MODULE.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
            "allowed_differing_auxiliary_tags": ["RG"],
            "analysis_scope_identity_contract": (
                MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT
            ),
            "missing_projection_policy": MODULE.RAW_IDENTITY_MISSING_POLICY,
            "conflicting_analysis_payload_policy": MODULE.RAW_IDENTITY_CONFLICT_POLICY,
            "failure_counts_materialized": False,
            "all_site_results_passed_invariant_validation": True,
        },
        "pass": True,
    }
    summary_path.write_text(
        json.dumps(summary, sort_keys=True) + "\n", encoding="utf-8"
    )
    raw_duplicate_path = candidate_table.parent / "raw_identity_duplicate_audit.tsv.gz"
    with gzip.open(raw_duplicate_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("sample\tchrom\tpos\n")
    oracle_path = candidate_table.parent / "oracle_cases.json"
    oracle_path.write_text("{}\n", encoding="utf-8")
    code = {
        role: MODULE.artifact(path)
        for role, path in MODULE.COOCCURRENCE_CODE_PATHS.items()
    }
    source_modes = {role: "0o444" for role in code}
    preflight_path = candidate_table.parent / "cooccurrence_preflight.json"
    preflight_sources = {
        "preflight": MODULE.artifact(MODULE.COOCCURRENCE_PREFLIGHT_SOURCE),
        **code,
    }
    preflight_modes = {role: "0o444" for role in preflight_sources}
    preflight = {
        "schema_name": "intersubmod.cooccurrence_task_contract_preflight",
        "schema_version": "2.0.0",
        "task_type": "B_comprehensive_validation",
        "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
        "source_authority": source_authority,
        "observed": {
            "task_count": 102_842,
            "raw_bam_identity_recovery": {
                "aggregate": {"site_tasks": 102_842},
                "equivalence_contract": MODULE.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
                "allowed_differing_auxiliary_tags": ["RG"],
                "analysis_scope_identity_contract": (
                    MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT
                ),
                "missing_projection_policy": MODULE.RAW_IDENTITY_MISSING_POLICY,
                "conflicting_analysis_payload_policy": (
                    MODULE.RAW_IDENTITY_CONFLICT_POLICY
                ),
                "failure_counts_materialized": False,
                "all_result_rows_passed_invariant_validation": True,
            },
        },
        "code": {
            "source_identity_before": preflight_sources,
            "source_identity_after": preflight_sources,
            "source_modes_before": preflight_modes,
            "source_modes_after": preflight_modes,
        },
        "pass": True,
    }
    preflight_path.write_text(
        json.dumps(preflight, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt_path = candidate_table.parent / "cooccurrence_receipt.json"
    receipt = {
        "schema_name": MODULE.COOCCURRENCE_RECEIPT_SCHEMA,
        "schema_version": MODULE.COOCCURRENCE_FORMAL_SCHEMA_VERSION,
        "pass": True,
        "task_type": (
            "B_comprehensive_validation" if full_scope else "A_partial_explicit_scope"
        ),
        "scope": "all_manifest_samples" if full_scope else "explicit_sample_subset",
        "full_scope": full_scope,
        "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
        "source_authority": source_authority,
        "release_status": "PRODUCER_PASS_PENDING_RUNNER_RECONCILIATION",
        "preflight_gate": {
            "schema_name": "intersubmod.cooccurrence_task_contract_preflight",
            "schema_version": "2.0.0",
            "task_count": 102_842,
            "pass": True,
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
        },
        "inputs": {
            "assignments": MODULE.artifact(assignments),
            "preflight_receipt": MODULE.artifact(preflight_path),
        },
        "outputs": {
            "site_tsv": MODULE.artifact(candidate_table),
            "summary_json": MODULE.artifact(summary_path),
            "raw_identity_duplicate_audit_tsv": MODULE.artifact(raw_duplicate_path),
            "case_json": MODULE.artifact(oracle_path),
        },
        "code": code,
        "source_lock": {
            "source_identity_before": code,
            "source_identity_after_compute": code,
            "source_identity_after_output": code,
            "source_modes_before": source_modes,
            "source_modes_after_compute": source_modes,
            "source_modes_after_output": source_modes,
            "all_sources_read_only_and_unchanged": True,
        },
        "counts": {
            "stable_sites": len(keys),
            "candidate_unique_key_count": len(set(keys)),
            "multi_marker_molecular_haplotype_base_candidates": selected_count,
            "raw_identity_missing_projection_policy": (
                MODULE.RAW_IDENTITY_MISSING_POLICY
            ),
            "raw_identity_conflicting_analysis_payload_policy": (
                MODULE.RAW_IDENTITY_CONFLICT_POLICY
            ),
            "raw_identity_failure_counts_materialized": False,
            "raw_identity_all_site_results_passed_invariant_validation": True,
        },
    }
    receipt_path.write_text(
        json.dumps(receipt, sort_keys=True) + "\n", encoding="utf-8"
    )
    release_path = candidate_table.parent / "cooccurrence_release_receipt.json"
    release = {
        "schema_name": MODULE.COOCCURRENCE_RELEASE_SCHEMA,
        "schema_version": MODULE.COOCCURRENCE_RELEASE_SCHEMA_VERSION,
        "task_type": "B_comprehensive_validation",
        "scope": "all_manifest_samples",
        "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
        "release_status": "RELEASE_RECONCILIATION_PASS",
        "source_authority": source_authority,
        "inputs": {"producer_receipt": MODULE.artifact(receipt_path)},
        "code": {
            role: MODULE.artifact(path)
            for role, path in MODULE.COOCCURRENCE_RELEASE_CODE_PATHS.items()
        },
        "source_modes": {
            role: "0o444" for role in MODULE.COOCCURRENCE_RELEASE_CODE_PATHS
        },
        "checks": {"fixture_release_reconciled": True},
        "pass": True,
    }
    if release_path.exists():
        release_path.chmod(0o644)
    release_path.write_text(json.dumps(release, sort_keys=True) + "\n", encoding="utf-8")
    release_path.chmod(0o444)
    return receipt_path


def validate_cooccurrence_fixture(
    candidate_table: Path, assignments: Path, receipt_path: Path
) -> dict[str, object]:
    selection = MODULE.resolve_selection_column(MODULE.candidate_table_fields(candidate_table))
    _, selected, n_input_rows = MODULE.load_candidates(
        candidate_table, selection.resolved_column
    )
    candidate_keys = MODULE.candidate_table_keys(candidate_table)
    selected_keys = [MODULE.site_key(row, source="test selected") for row in selected]
    return MODULE.validate_cooccurrence_receipt(
        receipt_path,
        release_receipt_path=(
            receipt_path.parent / "cooccurrence_release_receipt.json"
        ),
        candidate_table=candidate_table,
        assignments_path=assignments,
        n_candidate_rows=n_input_rows,
        candidate_keys=candidate_keys,
        selected_keys=selected_keys,
        selection=selection,
    )


def make_inputs(
    tmp_path: Path, *, include_unselected: bool = True
) -> tuple[Path, Path, dict[str, object]]:
    region, read_ids, labels, identities = make_region(tmp_path)
    candidate_table = tmp_path / "candidates.tsv.gz"
    assignments = tmp_path / "assignments.jsonl.gz"
    write_candidate_table(candidate_table, region, include_unselected=include_unselected)
    assignment = make_assignment(region, read_ids, labels, identities)
    write_assignments(assignments, [assignment])
    write_cooccurrence_receipt(candidate_table, assignments)
    return candidate_table, assignments, assignment


def cli_args(candidate_table: Path, assignments: Path, output_dir: Path) -> list[str]:
    return [
        "--candidate-table",
        str(candidate_table),
        "--assignments",
        str(assignments),
        "--cooccurrence-receipt",
        str(candidate_table.parent / "cooccurrence_receipt.json"),
        "--cooccurrence-release-receipt",
        str(candidate_table.parent / "cooccurrence_release_receipt.json"),
        "--output-dir",
        str(output_dir),
        "--permutations",
        "7",
        "--seeds",
        "3",
        "--seed",
        "77",
        "--assignment-ari-threshold",
        "0.8",
        "--min-valid-null-fraction",
        "0.8",
    ]


def read_site_rows(output_dir: Path) -> list[dict[str, str]]:
    path = output_dir / "strict_methyl_candidate_confirmation_sites.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def test_formal_parameters_and_canonical_selection_are_required_for_pass() -> None:
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

    noncanonical_configs = (
        MODULE.ConfirmationConfig(permutations=1_000, seeds=10),
        MODULE.ConfirmationConfig(permutations=999, seeds=11),
        MODULE.ConfirmationConfig(permutations=999, seeds=10, seed=7),
        MODULE.ConfirmationConfig(
            permutations=999, seeds=10, assignment_ari_threshold=0.0
        ),
        MODULE.ConfirmationConfig(
            permutations=999, seeds=10, min_valid_null_fraction=0.01
        ),
    )
    for config in noncanonical_configs:
        noncanonical = MODULE.flatten_result(
            result,
            formal_selection,
            cooccurrence,
            1,
            config,
        )
        assert noncanonical["strict_analysis_class"] == "EXPLORATORY_ONLY"
        assert noncanonical["strict_null_robustness_pass"] is False
        assert "strict_formal_parameter_gate" in noncanonical["strict_failure_reason"]

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
    assert "strict_formal_selection_contract_gate" in legacy["strict_failure_reason"]


def test_selection_resolution_prefers_formal_and_marks_legacy_as_alias() -> None:
    resolved = MODULE.resolve_selection_column(
        [MODULE.LEGACY_SELECTION_COLUMNS[0], MODULE.FORMAL_SELECTION_COLUMN]
    )
    assert resolved.resolved_column == MODULE.FORMAL_SELECTION_COLUMN
    assert resolved.legacy_fallback_used is False
    legacy = MODULE.resolve_selection_column([MODULE.LEGACY_SELECTION_COLUMNS[0]])
    assert legacy.formal_column == MODULE.FORMAL_SELECTION_COLUMN
    assert legacy.legacy_fallback_used is True


def test_formal_not_evaluable_is_not_counted_as_robustness_failure() -> None:
    result = make_flatten_result()
    result["modes"]["column"].update(
        {
            "observed_between_within": MODULE.F.SEP_MIN + 1.0,
            "effect_gate": True,
            "null_valid_gate": False,
        }
    )
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


def test_finite_below_effect_threshold_is_evaluable_robustness_failure() -> None:
    result = make_flatten_result()
    result["modes"]["column"].update(
        {
            "observed_between_within": MODULE.F.SEP_MIN - 0.01,
            "effect_gate": False,
            "null_valid_gate": False,
            "null_threshold_gate": False,
            "failure": "below_sep_min_or_undefined",
        }
    )
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

    assert row["strict_methyl_partition_robustness_evaluable"] is True
    assert row["strict_null_robustness_pass"] is False
    assert row["strict_methyl_partition_robustness_status"] == (
        "R1_STRICT_METHYL_PARTITION_ROBUSTNESS_FAIL"
    )
    assert row["strict_scientific_status"] == "ROBUSTNESS_FAIL"
    assert row["strict_not_evaluable_reason"] == ""
    assert "strict_column_effect_gate" in row["strict_failure_reason"]


def test_cooccurrence_receipt_locks_scope_hashes_counts_and_family(tmp_path: Path) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path)
    receipt_path = candidate_table.parent / "cooccurrence_receipt.json"
    contract = validate_cooccurrence_fixture(candidate_table, assignments, receipt_path)
    assert contract["full_scope_gate"] is True
    assert contract["candidate_table_row_count"] == 2
    assert contract["candidate_table_unique_key_count"] == 2
    assert contract["selected_candidate_row_count"] == 1
    assert contract["selected_candidate_unique_key_count"] == 1
    assert contract["postselection_family_size"] == 1
    assert contract["screen_assignment_artifact_hash_lock"] is True


def test_cooccurrence_receipt_rejects_partial_scope(tmp_path: Path) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path)
    receipt_path = write_cooccurrence_receipt(
        candidate_table, assignments, full_scope=False
    )
    with pytest.raises(
        MODULE.ConfirmationContractError,
        match="(task_type is not full scope|not Task-B full scope)",
    ):
        validate_cooccurrence_fixture(candidate_table, assignments, receipt_path)


def test_cooccurrence_receipt_rejects_site_and_assignment_hash_drift(
    tmp_path: Path,
) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path / "site")
    receipt_path = candidate_table.parent / "cooccurrence_receipt.json"
    selection = MODULE.resolve_selection_column(MODULE.candidate_table_fields(candidate_table))
    _, selected, n_input_rows = MODULE.load_candidates(
        candidate_table, selection.resolved_column
    )
    candidate_keys = MODULE.candidate_table_keys(candidate_table)
    selected_keys = [MODULE.site_key(row, source="test selected") for row in selected]
    with candidate_table.open("ab") as handle:
        handle.write(b"drift")
    with pytest.raises(
        MODULE.ConfirmationContractError, match="site_tsv artifact (size|sha256)"
    ):
        MODULE.validate_cooccurrence_receipt(
            receipt_path,
            candidate_table=candidate_table,
            assignments_path=assignments,
            n_candidate_rows=n_input_rows,
            candidate_keys=candidate_keys,
            selected_keys=selected_keys,
            selection=selection,
        )

    candidate_table, assignments, _ = make_inputs(tmp_path / "assignment")
    receipt_path = candidate_table.parent / "cooccurrence_receipt.json"
    with assignments.open("ab") as handle:
        handle.write(b"drift")
    with pytest.raises(
        MODULE.ConfirmationContractError, match="assignments artifact (size|sha256)"
    ):
        validate_cooccurrence_fixture(candidate_table, assignments, receipt_path)


def test_cooccurrence_receipt_rejects_candidate_family_mismatch(tmp_path: Path) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path)
    receipt_path = write_cooccurrence_receipt(
        candidate_table, assignments, family_size=2
    )
    with pytest.raises(MODULE.ConfirmationContractError, match="candidate family mismatch"):
        validate_cooccurrence_fixture(candidate_table, assignments, receipt_path)


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


def test_empirical_null_plus_one_and_postselection_diagnostic_calculations() -> None:
    assert MODULE.empirical_p_value(2.0, [0.5, 2.0, 3.0]) == pytest.approx(0.75)
    assert MODULE.empirical_p_value(2.0, []) is None
    assert MODULE.benjamini_hochberg([0.01, 0.04, 0.03, 0.20]) == pytest.approx(
        [0.04, 0.05333333333333334, 0.05333333333333334, 0.20]
    )
    assert MODULE.benjamini_yekutieli([0.01, 0.04, 0.03, 0.20]) == pytest.approx(
        [0.08333333333333333, 0.1111111111111111, 0.1111111111111111, 0.41666666666666663]
    )


def test_partition_comparison_is_label_invariant_but_exact() -> None:
    assert MODULE.same_partition(["a", "a", "b", "b"], [2, 2, 1, 1])
    assert not MODULE.same_partition(["a", "a", "b", "b"], [2, 1, 1, 1])


def test_fixed_k_leave_one_out_rejects_original_group_below_minimum() -> None:
    labels = ["a"] * 3 + ["b"] * 4
    methylation = np.vstack(
        [np.full((3, 6), 0.1), np.full((4, 6), 0.9)]
    )
    distance = MODULE.F.bernoulli_distance(methylation)
    result = MODULE.leave_one_out_stability(distance, labels, 2, 0.8)
    assert result["exact_partition_gate"] is False
    assert result["invalid"] >= 3


def test_inference_uses_prefixed_seed_not_modal_multiseed_representative(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    labels = ["a"] * 4 + ["b"] * 4
    calls: list[int] = []

    def fake_analyze(*_args: object, **kwargs: object) -> dict[str, object]:
        seeds = int(kwargs["seeds"])
        calls.append(seeds)
        return {
            "coarse_labels": labels if seeds == 1 else ["x"] * 4 + ["y"] * 4,
            "coarse_split_trace": [
                {
                    "empirical_p": 0.125 if seeds == 1 else 0.875,
                    "n_valid_null": 7,
                    "observed_between_within": 2.0,
                    "null_threshold": 1.0,
                    "passed": True,
                    "failure": "",
                }
            ],
            "coarse_ng": 2,
            "unstable": False,
            "modal_fraction": 1.0,
            "modal_assignment_pair_count": seeds * (seeds - 1) // 2,
            "modal_assignment_ari_median": 1.0,
            "modal_assignment_ari_min": 1.0,
        }

    monkeypatch.setattr(MODULE.F, "analyze_phylo", fake_analyze)
    monkeypatch.setattr(
        MODULE,
        "leave_one_out_stability",
        lambda *_args, **_kwargs: {
            "feasible": True,
            "evaluated": 8,
            "invalid": 0,
            "ari_median": 1.0,
            "ari_min": 1.0,
            "all_group_size_gate": True,
            "stability_gate": True,
            "exact_partition_count": 8,
            "exact_partition_gate": True,
        },
    )
    site = MODULE.PreparedSite(
        input_row={},
        source_row_number=2,
        sample="S",
        chrom="chr1",
        pos=100,
        region_dir=Path("/unused"),
        distance=np.zeros((8, 8)),
        methylation=np.zeros((8, 4)),
        kept_ids=[str(index) for index in range(8)],
        screen_labels=labels,
        n_reads_total=8,
        n_alt_raw=8,
        n_alt_after_peel=8,
    )
    result = MODULE.run_null_mode(
        site, "column", MODULE.ConfirmationConfig(permutations=7, seeds=4, seed=77)
    )
    assert calls == [1, 4]
    assert result["empirical_p_postselection_descriptive"] == pytest.approx(0.125)
    assert result["labels"] == labels


def test_both_null_modes_are_seeded_deterministic_and_quantized(tmp_path: Path) -> None:
    candidate_table, assignments_path, _ = make_inputs(tmp_path)
    _, selected, _ = MODULE.load_candidates(candidate_table, SELECTION_COLUMN)
    assignments = MODULE.load_assignments(assignments_path)
    site = MODULE.prepare_site(selected[0], assignments[("SYNTH", "chr1", 100)])
    config = MODULE.ConfirmationConfig(permutations=7, seeds=3, seed=77)

    for mode in MODULE.NULL_MODES:
        first = MODULE.run_null_mode(site, mode, config)
        second = MODULE.run_null_mode(site, mode, config)
        comparable = [
            "base_seed",
            "observed_between_within",
            "null_threshold",
            "empirical_p_postselection_descriptive",
            "n_valid_null",
            "modal_groups",
            "multiseed_ari_min",
            "screen_assignment_ari",
            "fixed_k_loo_partition_stability_ari_min",
            "labels",
        ]
        assert {key: first[key] for key in comparable} == {
            key: second[key] for key in comparable
        }
        empirical = first["empirical_p_postselection_descriptive"]
        if empirical is not None:
            scaled = empirical * (first["n_valid_null"] + 1)
            assert math.isclose(scaled, round(scaled), abs_tol=1e-8)


def test_strict_replay_digest_matches_typed_and_tsv_string_rows() -> None:
    row = {
        field: None
        for field in ["sample", "chrom", "pos", "ref", "alt", *MODULE.STRICT_OUTPUT_FIELDS]
    }
    row.update(
        sample="SYNTH",
        chrom="chr1",
        pos=100,
        ref="A",
        alt="C",
        strict_candidate_selection_gate=True,
        strict_column_group_sizes=[3, 5],
        strict_column_observed_between_within=0.125,
    )
    string_row = {
        field: "" if value is None else str(value) for field, value in row.items()
    }
    assert MODULE.strict_analysis_replay_record([row]) == (
        MODULE.strict_analysis_replay_record([string_row])
    )


def test_selection_scope_and_exploratory_row_are_preserved(tmp_path: Path) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path)
    output_dir = tmp_path / "output"
    assert MODULE.main(cli_args(candidate_table, assignments, output_dir)) == 0

    rows = read_site_rows(output_dir)
    assert len(rows) == 1
    row = rows[0]
    assert row["pos"] == "100"
    assert row[SELECTION_COLUMN] == "true"
    assert row["cooccurrence_fdr_gate"] == "true"
    assert row["strict_candidate_selection_gate"] == "True"
    assert row["strict_artifact_identity_gate"] == "True"
    assert row["strict_postselection_family_size"] == "1"
    assert row["strict_postselection_scope"] == MODULE.POSTSELECTION_SCOPE
    assert row["strict_postselection_fdr_calibrated"] == "False"
    assert "strict_candidate_set_fdr_gate" not in row
    assert row["strict_analysis_class"] == "EXPLORATORY_ONLY"
    assert row["strict_formal_parameter_gate"] == "False"
    assert row["strict_null_robustness_pass"] == "False"
    assert row["strict_scientific_status"] == "EXPLORATORY_ONLY"
    assert "strict_formal_parameter_gate" in row["strict_failure_reason"]

    summary = json.loads(
        (output_dir / "strict_methyl_candidate_confirmation_summary.json").read_text()
    )
    assert summary["selection_contract"]["n_input_rows"] == 2
    assert summary["selection_contract"]["n_selected_candidates"] == 1
    assert summary["selection_contract"]["n_not_selected"] == 1
    assert summary["counts"]["n_null_robustness_pass"] == 0
    assert summary["counts"]["n_null_robustness_fail_retained"] == 0
    assert summary["counts"]["n_exploratory_only"] == 1
    assert summary["postselection_diagnostic_contract"]["family_size"] == 1
    assert summary["postselection_diagnostic_contract"]["fdr_calibrated"] is False
    assert summary["execution_status"] == "EXECUTION_PASS"
    assert summary["pass_semantics"] == "execution_integrity_only_not_scientific_confirmation"
    assert "methyl-partition robustness only" in summary["guardrail"]

    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["contracts"]["screen_positive_strict_fail_rows_retained"] is True
    assert manifest["contracts"]["candidate_selection_relaxed"] is False
    assert manifest["contracts"]["postselection_fdr_claim_allowed"] is False
    assert manifest["contracts"]["formal_parameter_gate"] is False
    assert manifest["execution_status"] == "EXECUTION_PASS"
    assert set(manifest["code"]) == {
        "strict_producer",
        "focal_alt_cluster_lib",
        "source_authority_validator",
    }
    assert all(record["sha256"] for record in manifest["code"].values())
    runtime = manifest["runtime_environment_lock"]
    assert runtime["identity_before"] == runtime["identity_after_compute"]
    assert runtime["identity_before"]["interpreter_flags"] == {
        "isolated": 1,
        "no_user_site": 1,
        "safe_path": True,
    }
    assert manifest["analysis_replay"]["n_rows"] == 1
    expected_cache = tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME
    assert manifest["command"][:5] == [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={expected_cache}",
        str(SCRIPT.resolve()),
    ]


def test_site_tsv_gzip_is_byte_deterministic(tmp_path: Path) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path)
    first = tmp_path / "first"
    second = tmp_path / "second"
    MODULE.main(cli_args(candidate_table, assignments, first))
    MODULE.main(cli_args(candidate_table, assignments, second))

    filename = "strict_methyl_candidate_confirmation_sites.tsv.gz"
    assert (first / filename).read_bytes() == (second / filename).read_bytes()


def test_zero_selected_candidates_is_a_valid_audited_result(tmp_path: Path) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path, include_unselected=False)
    with gzip.open(candidate_table, "rt", encoding="ascii", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    rows[0][SELECTION_COLUMN] = "false"
    with gzip.open(candidate_table, "wt", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_cooccurrence_receipt(candidate_table, assignments)

    output_dir = tmp_path / "zero-output"
    assert MODULE.main(cli_args(candidate_table, assignments, output_dir)) == 0
    assert not (output_dir / "strict_methyl_candidate_confirmation_sites.tsv.gz").exists()
    assert not (output_dir / "strict_methyl_candidate_confirmation_summary.json").exists()
    receipt = json.loads((output_dir / "not_applicable_receipt.json").read_text())
    assert receipt["status"] == "NOT_APPLICABLE"
    assert receipt["execution_status"] == "NOT_APPLICABLE"
    assert receipt["reason"] == "ZERO_SELECTED_CANDIDATES"
    assert receipt["selection_column"] == MODULE.DEFAULT_SELECTION_COLUMN
    assert receipt["counts"]["n_selected_candidates"] == 0
    assert receipt["counts"]["n_methyl_partition_robustness_evaluable"] == 0
    assert receipt["counts"]["n_methyl_partition_robustness_not_evaluable"] == 0
    assert receipt["counts"]["n_null_robustness_pass"] == 0
    assert receipt["counts"]["n_null_robustness_fail_retained"] == 0
    assert receipt["counts"]["n_exploratory_only"] == 0
    assert receipt["analysis_replay"] == MODULE.strict_analysis_replay_record([])


def test_assignment_identity_mismatch_hard_fails_without_output(tmp_path: Path) -> None:
    candidate_table, assignments, assignment = make_inputs(tmp_path)
    assignment["read_names"][0] = "wrong-read-name"
    write_assignments(assignments, [assignment])
    write_cooccurrence_receipt(candidate_table, assignments)
    output_dir = tmp_path / "must-not-exist"

    with pytest.raises(MODULE.ConfirmationContractError, match="read_names.*mismatch"):
        MODULE.main(cli_args(candidate_table, assignments, output_dir))
    assert not output_dir.exists()


def test_primary_artifact_hash_mismatch_hard_fails_without_output(tmp_path: Path) -> None:
    candidate_table, assignments, _ = make_inputs(tmp_path)
    region = tmp_path / "region"
    with (region / "methylation/methylation.csv").open("a", encoding="ascii") as handle:
        handle.write("\n")
    output_dir = tmp_path / "must-not-exist"

    with pytest.raises(MODULE.ConfirmationContractError, match="methylation_matrix artifact (size|sha256) mismatch"):
        MODULE.main(cli_args(candidate_table, assignments, output_dir))
    assert not output_dir.exists()


def test_candidate_assignment_allele_mismatch_hard_fails(tmp_path: Path) -> None:
    candidate_table, assignments, assignment = make_inputs(tmp_path)
    assignment["posthoc"]["alt"] = "G"
    write_assignments(assignments, [assignment])
    write_cooccurrence_receipt(candidate_table, assignments)
    output_dir = tmp_path / "must-not-exist"

    with pytest.raises(MODULE.ConfirmationContractError, match="Assignment/candidate alt mismatch"):
        MODULE.main(cli_args(candidate_table, assignments, output_dir))
    assert not output_dir.exists()


def test_existing_output_directory_is_rejected_before_reading_inputs(tmp_path: Path) -> None:
    output_dir = tmp_path / "existing"
    output_dir.mkdir()
    args = cli_args(tmp_path / "missing.tsv.gz", tmp_path / "missing.jsonl.gz", output_dir)

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.main(args)


def test_selection_column_cannot_be_overridden_from_cli(tmp_path: Path) -> None:
    args = cli_args(tmp_path / "candidates.tsv.gz", tmp_path / "assignments.jsonl.gz", tmp_path / "out")
    with pytest.raises(SystemExit):
        MODULE.build_parser().parse_args([*args, "--selection-column", "weaker_gate"])


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
