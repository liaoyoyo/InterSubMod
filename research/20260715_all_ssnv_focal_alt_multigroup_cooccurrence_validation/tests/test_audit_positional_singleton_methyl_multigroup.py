from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "audit_positional_singleton_methyl_multigroup.py"
)
SPEC = importlib.util.spec_from_file_location("singleton_auditor", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_audit_schema_uses_stable_main_authority_key() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    assert MODULE.SCHEMA_VERSION == "2.0.0"
    assert '"main_task_b_source_authority": source_authority_validation' in source
    assert '"v5_source_authority"' not in source


def item(dataset: str, chrom: str, pos: int, alt: str = "T") -> dict[str, object]:
    key = (dataset, chrom, pos, "C", alt)
    return {
        "key": key,
        "pos": pos,
        "component_id": "",
        "component_size": 0,
        "branch": "",
        "selected_for_read_census": False,
    }


def bind_component(
    rows: list[dict[str, object]], *, component_id: str, branch: str
) -> list[dict[str, object]]:
    for row in rows:
        row["component_id"] = component_id
        row["component_size"] = len(rows)
        row["branch"] = branch
    return rows


def test_component_boundary_50000_is_connected_and_50001_is_singleton() -> None:
    first = item("D", "chr1", 100)
    second = item("D", "chr1", 50_100)
    third = item("D", "chr1", 100_101)
    bind_component(
        [first, second],
        component_id="chr1:100-50100",
        branch="retained",
    )
    bind_component(
        [third],
        component_id="chr1:100101-100101",
        branch="positional_singleton",
    )
    singletons, counts = MODULE.recompute_components(
        {("D", "chr1"): [third, first, second]}
    )
    assert set(singletons) == {third["key"]}
    assert singletons[third["key"]]["left_gap_bp"] == 50_001
    assert singletons[third["key"]]["right_gap_bp"] is None
    assert counts["component_identity_mismatch"] == 0
    assert counts["singleton_branch_mismatch"] == 0


def test_transitive_component_can_span_more_than_50000() -> None:
    rows = [
        item("D", "chr2", 100),
        item("D", "chr2", 40_100),
        item("D", "chr2", 80_100),
    ]
    bind_component(
        rows,
        component_id="chr2:100-80100",
        branch="retained",
    )
    singletons, counts = MODULE.recompute_components({("D", "chr2"): rows})
    assert singletons == {}
    assert counts["component_identity_mismatch"] == 0
    assert counts["recomputed_singletons"] == 0


def test_same_coordinate_different_alt_is_not_singleton() -> None:
    rows = [item("D", "chr3", 500, "A"), item("D", "chr3", 500, "G")]
    bind_component(rows, component_id="chr3:500-500", branch="retained")
    singletons, counts = MODULE.recompute_components({("D", "chr3"): rows})
    assert singletons == {}
    assert counts["component_identity_mismatch"] == 0


def test_m1_recompute_requires_evaluable_stable_membership() -> None:
    base = {
        "dataset": "D",
        "chrom": "chr1",
        "pos": "1",
        "ref": "C",
        "alt": "T",
        "m1_evaluable": "true",
        "coarse_ng": "2",
        "unstable": "false",
        "modal_assignment_ari_min": "0.8",
    }
    assert MODULE.m1_recomputed(base)
    assert not MODULE.m1_recomputed({**base, "modal_assignment_ari_min": "0.799"})
    assert not MODULE.m1_recomputed({**base, "unstable": "true"})
    assert not MODULE.m1_recomputed(
        {
            **base,
            "m1_evaluable": "false",
            "coarse_ng": "",
            "unstable": "",
            "modal_assignment_ari_min": "",
        }
    )


def test_hcc1395_dorado_does_not_add_biological_n() -> None:
    assert MODULE.dataset_role("HCC1395") == ("technical_pair_primary", 1)
    assert MODULE.dataset_role("HCC1395_DORADO") == ("technical_replicate", 0)


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
def test_optional_float_rejects_non_finite_values(value: str) -> None:
    with pytest.raises(MODULE.SingletonAuditError, match="Non-finite"):
        MODULE.parse_optional_float(value)


def file_identity(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def write_json(path: Path, value: object, *, read_only: bool = False) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if read_only:
        path.chmod(0o444)


def configure_synthetic_expectations(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        MODULE,
        "EXPECTED_COUNTS",
        {
            "all_dataset_sites": 7,
            "m1_stable_all": 1,
            "positional_singleton": 7,
            "singleton_m1_evaluable": 7,
            "singleton_m1_not_evaluable": 0,
            "singleton_m1_flagged": 1,
            "singleton_m2_pass": 1,
            "singleton_m2_fail": 0,
            "singleton_m2_axis_indeterminate": 0,
            "singleton_m2_group_count_gt10": 0,
            "singleton_m2_not_run": 6,
            "singleton_min_finite_nearest_gap_bp": -1,
        },
    )
    monkeypatch.setattr(
        MODULE,
        "EXPECTED_BRANCH_COUNTS",
        {
            "max_snv_excluded": 0,
            "positional_singleton": 7,
            "retained": 0,
        },
    )
    monkeypatch.setattr(MODULE, "EXPECTED_M1_NOT_EVALUABLE_REASONS", {})
    monkeypatch.setattr(
        MODULE,
        "EXPECTED_M2_FAIL_SUBTYPES",
        {
            "HP_AXIS_CONFOUND": 0,
            "NOT_PHASE_ANCHORED_ROBUST": 0,
            "TECHNICAL_AXIS_CONFOUND": 0,
        },
    )
    monkeypatch.setattr(
        MODULE,
        "EXPECTED_M2_GLOBAL_COUNTS",
        {
            "all_rows": 7,
            "eligible": 1,
            "evaluable_ineligible": 0,
            "m1_stable_rows": 1,
            "not_evaluable_axis_indeterminate": 0,
            "not_evaluable_group_count_gt10": 0,
        },
    )
    monkeypatch.setattr(
        MODULE,
        "EXPECTED_GLOBAL_TRUTH_COUNTS",
        {"truth_tp": 7, "truth_fp": 0, "truth_unassessed": 0},
    )


def synthetic_release_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    first_truth_label: str = "TP",
    first_ari_min: object = 1.0,
    first_strict_status: str = "NOT_EVALUATED_AT_M1_SCREEN",
    mock_source_authority: bool = True,
) -> dict[str, Path]:
    configure_synthetic_expectations(monkeypatch)
    if mock_source_authority:
        monkeypatch.setattr(
            MODULE,
            "verify_source_authority",
            lambda **_kwargs: {
                "authority_id": "TEST_ONLY_SYNTHETIC_AUTHORITY",
                "source_set_sha256": "0" * 64,
                "authorized_git_head": "0" * 40,
                "authorized_source_roles": [],
                "external_approvals": 2,
                "pass": True,
            },
        )
    topic_root = SCRIPT.parents[1]
    repository_root = topic_root.parents[1]
    claim_contract = topic_root / "claim-contract-v5.md"
    m2_source = topic_root / "scripts" / "audit_independent_m2_gate_recount.py"
    source_authority = (
        repository_root
        / "docs"
        / "provenance"
        / "source_authorities"
        / "20260718_all_ssnv_focal_alt_release_source_authority.v4.json"
    )
    source_authority_approval = source_authority.with_suffix(
        ".approval.json"
    )
    source_authority_signature = Path(
        str(source_authority_approval) + ".ed25519.sig"
    )
    source_authority_public_key = Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
        "20260718_all_ssnv_v4/ed25519_public.pem"
    )

    datasets = MODULE.EXPECTED_DATASETS
    biological_ids = {
        dataset: ("HCC1395" if dataset == "HCC1395_DORADO" else dataset)
        for dataset in datasets
    }
    manifest_samples = []
    for dataset in datasets:
        manifest_samples.append(
            {
                "sample": dataset,
                "biological_id": biological_ids[dataset],
                "counts": {
                    "all_ssnv": 1,
                    "truth_tp": 1,
                    "truth_fp": 0,
                    "truth_unassessed": 0,
                    "ledger_branches": {
                        "max_snv_excluded": 0,
                        "positional_singleton": 1,
                        "retained": 0,
                    },
                },
                "checks": {"synthetic_input_pass": True},
            }
        )
    manifest = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_input_manifest",
        "schema_version": "1.0.0",
        "pass": True,
        "task_type": MODULE.EXPECTED_TASK_TYPE,
        "scope": MODULE.EXPECTED_SCOPE,
        "totals": {
            "all_ssnv": 7,
            "truth_tp": 7,
            "truth_fp": 0,
            "truth_unassessed": 0,
            "ledger_branches": {
                "max_snv_excluded": 0,
                "positional_singleton": 7,
                "retained": 0,
            },
        },
        "samples": manifest_samples,
    }
    input_manifest = tmp_path / "input_manifest.json"
    write_json(input_manifest, manifest)

    strict_reason = MODULE.STRICT_NOT_EVALUATED_REASON
    site_rows: list[dict[str, object]] = []
    stable_dataset = datasets[0]
    for index, dataset in enumerate(datasets):
        stable = dataset == stable_dataset
        site_rows.append(
            {
                "dataset": dataset,
                "sample": dataset,
                "biological_id": biological_ids[dataset],
                "truth_label": (
                    first_truth_label if index == 0 else "TP"
                ),
                "chrom": "chr1",
                "pos": 100,
                "ref": "C",
                "alt": "T",
                "filter": "PASS",
                "ssnv_branch": "positional_singleton",
                "component_id": "chr1:100-100",
                "component_size": 1,
                "ledger_selected_for_read_census": False,
                "latest_tag_join_status": "PASS",
                "m1_stability_gate_contract": (
                    "coarse_ng>=2 AND not unstable AND "
                    "modal_assignment_ari_min>=0.8"
                ),
                "m1_evaluable": True,
                "m1_not_evaluable_reason": "",
                "analysis_status": "evaluable",
                "n_alt_matrix": 1000,
                "n_alt_after_peel": 1000,
                "n_alt_peeled": 0,
                "coarse_ng": 2 if stable else 1,
                "unstable": False,
                "modal_assignment_ari_min": (
                    first_ari_min if index == 0 else 1.0
                ),
                "stable_null_multigroup": stable,
                "cluster_sizes": (
                    '{"1":500,"2":500}' if stable else '{"1":1000}'
                ),
                "hp_axis_confound": False,
                "technical_axis_confound": False,
                "residual_unexplained_multigroup": stable,
                "phase_anchored_robust_epigenetic_candidate": stable,
                "hp_exact_v": "",
                "hp_exact_p_perm": "",
                "hp_exact_n": 1000,
                "hp_exact_aligned": False,
                "hp_family_v": "",
                "hp_family_p_perm": "",
                "hp_family_n": 1000,
                "hp_family_aligned": False,
                "strand_v": "",
                "strand_p_perm": "",
                "strand_n": 1000,
                "strand_aligned": False,
                "start_eta2": 0.0,
                "start_p_perm": 1.0,
                "start_n": 1000,
                "start_aligned": False,
                "end_eta2": 0.0,
                "end_p_perm": 1.0,
                "end_n": 1000,
                "end_aligned": False,
                "length_eta2": 0.0,
                "length_p_perm": 1.0,
                "length_n": 1000,
                "length_aligned": False,
                "mapq_eta2": 0.0,
                "mapq_p_perm": 1.0,
                "mapq_n": 1000,
                "mapq_aligned": False,
                "cpg_called_eta2": 0.0,
                "cpg_called_p_perm": 1.0,
                "cpg_called_n": 1000,
                "cpg_called_aligned": False,
                "strict_methyl_partition_robustness_status": (
                    first_strict_status
                    if index == 0
                    else "NOT_EVALUATED_AT_M1_SCREEN"
                ),
                "strict_methyl_partition_robustness_not_evaluable_reason": (
                    strict_reason
                ),
                "strict_confirm_status": "NOT_RUN",
                "strict_confirm_candidate": stable,
                "strict_confirm_candidate_is_formal_r1_claim": False,
                "strict_confirm_reason": strict_reason,
            }
        )
    site_results = tmp_path / "site_results.tsv.gz"
    with gzip.open(site_results, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(site_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(site_rows)

    reads = [
        {
            "read_id": f"read-{index}",
            "read_name": f"name-{index}",
            "label": "1" if index < 500 else "2",
            "latest_hp": "1-1",
            "strand": "+",
        }
        for index in range(1000)
    ]
    assignment = {
        "schema_name": MODULE.ASSIGNMENT_SCHEMA,
        "schema_version": "1.0.0",
        "dataset": stable_dataset,
        "chrom": "chr1",
        "pos": 100,
        "ref": "C",
        "alt": "T",
        "posthoc": {"ref": "C", "alt": "T"},
        "core_reads": reads,
        "read_ids": [read["read_id"] for read in reads],
        "read_names": [read["read_name"] for read in reads],
    }
    stable_assignments = tmp_path / "stable_assignments.jsonl.gz"
    with gzip.open(stable_assignments, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(assignment, sort_keys=True) + "\n")

    completed = {
        dataset: {
            "validation": {
                "expected_vcf_sites": 1,
                "reads_files": 1,
                "bernoulli_matrix_files": 1,
                "methylation_files": 1,
            }
        }
        for dataset in datasets
    }
    screen_run_manifest = tmp_path / "screen_run_manifest.json"
    write_json(
        screen_run_manifest,
        {
            "schema_name": (
                "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest"
            ),
            "schema_version": "1.2.0",
            "pass": True,
            "status": "EXECUTION_PASS",
            "input_manifest": file_identity(input_manifest),
            "outputs": {
                "site_results": file_identity(site_results),
                "stable_assignments": file_identity(stable_assignments),
            },
            "counts": {
                "expected_sites": 7,
                "processed_sites": 7,
                "stable_assignment_records": 1,
            },
            "contracts": {
                "m1_stability_gate_contract": (
                    "coarse_ng>=2 AND not unstable AND "
                    "modal_assignment_ari_min>=0.8"
                ),
                "strict_methyl_partition_robustness_status": (
                    "NOT_EVALUATED_AT_M1_SCREEN"
                ),
                "strict_confirm_status_legacy_alias": "NOT_RUN",
                "strict_confirm_candidate_is_formal_r1_claim": False,
                "latest_hp_ps_terminal_join": {
                    "pass": True,
                    "all_sites_pass": True,
                },
            },
            "completed_dataset_runs": completed,
        },
    )

    tree_contract_audit = tmp_path / "tree_contract_audit.json"
    write_json(
        tree_contract_audit,
        {
            "schema_name": "intersubmod.latest_tree_input_contract_audit",
            "schema_version": "1.0.0",
            "pass": True,
            "scope": MODULE.EXPECTED_SCOPE,
            "totals": manifest["totals"],
            "input_manifest": file_identity(input_manifest),
            "top_level_checks": {
                key: True for key in MODULE.EXPECTED_TREE_CHECKS
            },
        },
    )

    primary_artifact_audit = tmp_path / "primary_artifact_audit.json"
    write_json(
        primary_artifact_audit,
        {
            "schema_name": "intersubmod.stable_primary_artifact_audit",
            "schema_version": "2.0.0",
            "pass": True,
            "task_type": MODULE.EXPECTED_TASK_TYPE,
            "formal_task_b_release_eligible": True,
            "inputs": {
                "site_results": file_identity(site_results),
                "stable_assignments": file_identity(stable_assignments),
            },
            "counts": {
                "stable_sites": 1,
                "assignment_records": 1,
                "primary_artifacts_expected": 308_526,
                "primary_artifacts_verified": 308_526,
            },
            "verification": {
                "stable_site_assignment_key_sets_exact": True,
                "path_size_sha256_verified": True,
                "errors": 0,
                "artifact_roles_exact": [
                    "distance_matrix",
                    "methylation_matrix",
                    "reads",
                ],
                "artifact_set_sha256": "synthetic-artifact-set",
            },
            "source_authority": {
                "authority_id": MODULE.EXPECTED_AUTHORITY_ID,
                "authority_manifest": file_identity(source_authority),
                "detached_approval": file_identity(source_authority_approval),
                "detached_approval_signature": file_identity(
                    source_authority_signature
                ),
                "external_public_key": file_identity(
                    source_authority_public_key
                ),
                    "source_set_sha256": "0" * 64,
                "claim_contract": file_identity(claim_contract),
                "pass": True,
            },
        },
        read_only=True,
    )

    independent_m2_receipt = tmp_path / "independent_m2_receipt.json"
    write_json(
        independent_m2_receipt,
        {
            "schema_name": "intersubmod.independent_m2_gate_recount",
            "schema_version": "2.0.0",
            "pass": True,
            "task_type": MODULE.EXPECTED_TASK_TYPE,
            "contract": MODULE.EXPECTED_M2_CONTRACT,
            "checks": {key: True for key in MODULE.EXPECTED_M2_CHECKS},
            "counts": MODULE.EXPECTED_M2_GLOBAL_COUNTS,
            "logic_independence": MODULE.EXPECTED_M2_LOGIC_INDEPENDENCE,
            "inputs": {
                "site_results": file_identity(site_results),
                "stable_assignments": file_identity(stable_assignments),
                "claim_contract": file_identity(claim_contract),
            },
            "code": {"independent_recount": file_identity(m2_source)},
        },
        read_only=True,
    )
    return {
        "site_results": site_results,
        "stable_assignments": stable_assignments,
        "input_manifest": input_manifest,
        "screen_run_manifest": screen_run_manifest,
        "tree_contract_audit": tree_contract_audit,
        "primary_artifact_audit": primary_artifact_audit,
        "source_authority": source_authority,
        "source_authority_approval": source_authority_approval,
        "source_authority_signature": source_authority_signature,
        "source_authority_public_key": source_authority_public_key,
        "independent_m2_auditor": m2_source,
        "independent_m2_receipt": independent_m2_receipt,
        "claim_contract": claim_contract,
        "output_dir": tmp_path / "output",
    }


def test_build_audit_end_to_end_uses_atomic_release(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = synthetic_release_inputs(tmp_path, monkeypatch)
    result = MODULE.build_audit(**inputs)
    assert result["pass"] is True
    assert result["atomic_release"] is True
    output_dir = Path(result["output_dir"])
    assert output_dir.stat().st_mode & 0o777 == 0o555
    assert (output_dir / "_SUCCESS.json").is_file()
    summary = json.loads(
        (output_dir / "positional_singleton_audit_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["checks"]["dataset_counts_match_bound_manifest"] is True
    assert summary["status_census"]["m2"] == {
        "PASS": 1,
        "FAIL": 0,
        "NOT_EVALUABLE": 0,
        "NOT_RUN": 6,
    }
    assert summary["rates"]["m2_pass_fraction_of_all_singletons"][
        "numerator"
    ] == 1
    assert summary["rates"]["m2_pass_fraction_of_all_singletons"][
        "denominator"
    ] == 7


def test_build_audit_rejects_unknown_truth_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = synthetic_release_inputs(
        tmp_path, monkeypatch, first_truth_label="UNKNOWN"
    )
    with pytest.raises(MODULE.SingletonAuditError, match="Unknown truth label"):
        MODULE.build_audit(**inputs)
    assert not inputs["output_dir"].exists()


def test_build_audit_rejects_non_finite_m1_value_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = synthetic_release_inputs(
        tmp_path, monkeypatch, first_ari_min="inf"
    )
    with pytest.raises(MODULE.SingletonAuditError, match="Non-finite"):
        MODULE.build_audit(**inputs)
    assert not inputs["output_dir"].exists()


def test_build_audit_rejects_strict_r1_drift_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = synthetic_release_inputs(
        tmp_path, monkeypatch, first_strict_status="PASS"
    )
    with pytest.raises(MODULE.SingletonAuditError, match="strict"):
        MODULE.build_audit(**inputs)
    assert not inputs["output_dir"].exists()


def test_build_audit_rejects_legacy_v4_authority_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = synthetic_release_inputs(
        tmp_path, monkeypatch, mock_source_authority=False
    )
    original = inputs["source_authority_signature"]
    tampered = tmp_path / "tampered_authority.ed25519.sig"
    tampered.write_bytes(original.read_bytes() + b"x")
    tampered.chmod(0o444)
    inputs["source_authority_signature"] = tampered
    with pytest.raises(
        MODULE.SingletonAuditError, match="Legacy or noncanonical"
    ):
        MODULE.build_audit(**inputs)
    assert not inputs["output_dir"].exists()


def test_build_audit_refuses_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = synthetic_release_inputs(tmp_path, monkeypatch)
    inputs["output_dir"].mkdir()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.build_audit(**inputs)


def test_rename_no_replace_preserves_concurrent_target(tmp_path: Path) -> None:
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.mkdir()
    target.mkdir()
    marker = target / "owner.txt"
    marker.write_text("existing publisher\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="No-replace"):
        MODULE.rename_no_replace(source, target)
    assert marker.read_text(encoding="utf-8") == "existing publisher\n"
    assert source.is_dir()
