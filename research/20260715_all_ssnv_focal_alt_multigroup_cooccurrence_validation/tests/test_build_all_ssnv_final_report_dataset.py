from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import importlib.util
import json
import sys
from collections import Counter
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "build_all_ssnv_final_report_dataset.py"
)
SPEC = importlib.util.spec_from_file_location("build_all_ssnv_final_report_dataset", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


@pytest.fixture(autouse=True)
def test_only_release_source_authority(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    authority = {"authority_id": "TEST_ONLY_UNSIGNED_AUTHORITY", "pass": True}
    authority_path = tmp_path / "TEST_ONLY_UNSIGNED_AUTHORITY.json"
    claim_contract_path = tmp_path / "TEST_ONLY_CLAIM_CONTRACT.md"
    authority_path.write_text(json.dumps(authority) + "\n", encoding="utf-8")
    claim_contract_path.write_text("TEST ONLY - NON-RELEASE CLAIM CONTRACT\n", encoding="utf-8")
    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY,
        "validate_release_source_authority",
        lambda *_args, **_kwargs: authority,
    )
    monkeypatch.setattr(MODULE.SOURCE_AUTHORITY, "AUTHORITY_PATH", authority_path)
    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY,
        "CLAIM_CONTRACT_PATH",
        claim_contract_path,
    )


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    restore_read_only = path.exists() and oct(path.stat().st_mode & 0o777) == "0o444"
    if restore_read_only:
        path.chmod(0o644)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if restore_read_only:
        path.chmod(0o444)


def test_downstream_primary_consumer_receipts_include_independent_m2_in_order(
    tmp_path: Path,
) -> None:
    bundle = MODULE.InputBundle(
        manifest=tmp_path / "manifest.json",
        screen_dir=tmp_path / "screen",
        cooccurrence_dir=tmp_path / "cooccurrence",
        strict_dir=tmp_path / "strict",
        tumor_ref_dir=tmp_path / "tumor-ref",
        primary_artifact_audit_pre=tmp_path / "pre.json",
        primary_artifact_audit_post=tmp_path / "post.json",
        cooccurrence_preflight=tmp_path / "cooccurrence-preflight.json",
        independent_m2_audit=tmp_path / "independent-m2.json",
    )
    strict_receipt = tmp_path / "strict" / "run_manifest.json"
    matched_paired = tmp_path / "matched" / "run_receipt.json"
    matched_analysis = tmp_path / "matched-analysis" / "run_receipt.json"

    assert MODULE.downstream_primary_consumer_receipts(
        bundle,
        {"receipt_path": strict_receipt},
        {
            "paired_receipt_path": matched_paired,
            "receipt_path": matched_analysis,
        },
    ) == [
            bundle.cooccurrence_preflight,
            bundle.cooccurrence_receipt,
            bundle.cooccurrence_release_receipt,
            bundle.tumor_ref_receipt,
        bundle.independent_m2_audit,
        strict_receipt,
        matched_paired,
        matched_analysis,
    ]


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    restore_read_only = path.exists() and oct(path.stat().st_mode & 0o777) == "0o444"
    if restore_read_only:
        path.chmod(0o644)
    fields = fields or list(rows[0])
    opener = gzip.open if path.name.endswith(".gz") else path.open
    kwargs = {"mode": "wt", "encoding": "utf-8", "newline": ""}
    with opener(path, **kwargs) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if restore_read_only:
        path.chmod(0o444)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_header(path: Path) -> tuple[str, ...]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return tuple(reader.fieldnames or ())


def artifact(path: Path) -> dict:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": MODULE.sha256(path),
    }


def producer_provenance(module: object, role: str) -> dict:
    code = module.capture_source_identity()
    modes = module.capture_source_modes()
    assert set(code) == {role}
    assert modes == {role: "0o444"}
    return {
        "source_authority": MODULE.current_release_source_authority(),
        "code": code,
        "source_lock": {
            "source_identity_before": code,
            "source_identity_after_compute": code,
            "source_modes_before": modes,
            "source_modes_after_compute": modes,
            "all_sources_read_only_and_unchanged": True,
        },
    }


def ast_assignment(path: Path, name: str):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if not any(isinstance(target, ast.Name) and target.id == name for target in targets):
            continue
        value = node.value
        if isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and value.func.id in {
            "set",
            "frozenset",
            "tuple",
        }:
            payload = ast.literal_eval(value.args[0])
            return {
                "set": set,
                "frozenset": frozenset,
                "tuple": tuple,
            }[value.func.id](payload)
        return ast.literal_eval(value)
    raise AssertionError(f"Missing AST assignment {name} in {path}")


def producer_string_literals(path: Path) -> set[str]:
    return {
        node.value
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }


def test_canonical_task_b_paths_bind_fresh_versions_and_both_matched_branches() -> None:
    expected_names = {
        "cooccurrence_dir": (
            "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
            "source_locked_command_parity"
        ),
        "strict_dir": (
            "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5"
        ),
        "primary_artifact_audit_pre": (
            "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
        ),
        "primary_artifact_audit_post": (
            "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json"
        ),
        "cooccurrence_preflight": (
            "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
        ),
        "matched_normal_dir": (
            "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5"
        ),
        "cn_ccf_annotations": (
            "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
        ),
    }
    for field, expected_name in expected_names.items():
        assert MODULE.CANONICAL_TASK_B_PATHS[field].name == expected_name

    matched_names = {path.name for path in MODULE.CANONICAL_MATCHED_NORMAL_DIRS}
    assert matched_names == {
        "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
        "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
    }
    for matched_normal_dir in MODULE.CANONICAL_MATCHED_NORMAL_DIRS:
        paths = dict(MODULE.CANONICAL_TASK_B_PATHS)
        paths["matched_normal_dir"] = matched_normal_dir
        MODULE.validate_canonical_task_b_paths(
            MODULE.InputBundle(**paths),
            MODULE.CANONICAL_FINAL_DATASET_DIR,
        )

    commands = MODULE.canonical_task_b_final_builder_commands(
        MODULE.CANONICAL_FINAL_DATASET_DIR
    )
    assert len(commands) == 2
    expected_prefix = MODULE.STRICT_PRODUCER.canonical_python_prefix(
        MODULE.CANONICAL_TASK_B_PATHS["strict_dir"]
    )
    assert all(command[:4] == expected_prefix for command in commands)
    assert all(command[4] == str(SCRIPT.resolve()) for command in commands)
    matched_argument_names = {
        Path(command[command.index("--matched-normal-dir") + 1]).name
        for command in commands
    }
    assert matched_argument_names == matched_names


def assigned_string_values(path: Path, function_name: str, variable_name: str) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == function_name
    )
    values = set()
    for node in ast.walk(function):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.Constant):
            continue
        if not isinstance(node.value.value, str):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == variable_name
            for target in node.targets
        ):
            values.add(node.value.value)
    return values


def screen_rows(*, eleven_group_site: bool = False) -> list[dict]:
    specs = [
        ("HCC1395", "chr1", 100, "A", "C", "TP", True, True, True),
        ("HCC1395", "chr2", 200, "G", "T", "FP", True, False, False),
        ("HCC1395_DORADO", "chr1", 100, "A", "C", "TP", True, True, True),
        ("COLO829", "chr3", 300, "C", "A", "TP", False, False, False),
        ("H1437", "chr4", 400, "T", "G", "FP", True, False, False),
        ("H2009", "chr5", 500, "A", "G", "UNASSESSED", True, True, False),
        ("HCC1937", "chr6", 600, "C", "T", "TP", False, False, False),
        ("HCC1954", "chr7", 700, "G", "A", "UNASSESSED", True, False, False),
    ]
    rows = []
    for index, (sample, chrom, pos, ref, alt, truth, evaluable, stable, robust) in enumerate(
        specs, 1
    ):
        hp_confound = stable and not robust
        n_reads = index + 2
        row = {
                "sample": sample,
                "biological_id": MODULE.BIOLOGICAL_IDS[sample],
                "truth_label": truth,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "analysis_status": "evaluable" if evaluable else "insufficient_alt_reads",
                "stable_null_multigroup": stable,
                "modal_assignment_ari_min": 1.0 if stable else "",
                "cluster_sizes": json.dumps({"g1": 100, "g2": 100}) if stable else "{}",
                "m2_categorical_level_counts": json.dumps(
                    {"hp_exact": 2, "hp_family": 2, "strand": 2}
                ),
                "hp_axis_confound": hp_confound,
                "technical_axis_confound": False,
                "residual_unexplained_multigroup": stable and not hp_confound,
                "phase_anchored_robust_epigenetic_candidate": robust,
                "ssnv_branch": "retained",
                "component_id": f"component_{sample}_{chrom}_{pos}",
                "selected_group_id": "group_1",
                "evidence_tier": "E4" if robust else "E2" if stable else "E0",
                "n_alt_raw": 12 if stable else 3,
                "n_alt_after_peel": 10 if stable else 3,
                "latest_tag_join_status": "PASS",
                "latest_tag_rows_fetched": n_reads + 2,
                "latest_tag_rows_eligible": n_reads + 1,
                "latest_tag_reads_joined": n_reads,
                "latest_tag_ps_present": n_reads - 1,
                "latest_tag_projection_multimatch_reads": 0,
                "source_hp_replaced_reads": 1,
                "n_reads_total": n_reads,
        }
        for prefix in ("hp_exact", "hp_family", "strand"):
            row.update(
                {
                    f"{prefix}_v": "",
                    f"{prefix}_p_perm": "",
                    f"{prefix}_n": "",
                    f"{prefix}_aligned": False,
                }
            )
        for prefix in ("start", "end", "length", "mapq", "cpg_called"):
            row.update(
                {
                    f"{prefix}_eta2": "",
                    f"{prefix}_p_perm": "",
                    f"{prefix}_n": "",
                    f"{prefix}_aligned": False,
                }
            )
        if stable:
            for prefix in ("hp_exact", "hp_family", "strand"):
                row.update(
                    {
                        f"{prefix}_v": 0.1,
                        f"{prefix}_p_perm": 0.5,
                        f"{prefix}_n": 200,
                    }
                )
            for prefix in ("start", "end", "length", "mapq", "cpg_called"):
                row.update(
                    {
                        f"{prefix}_eta2": 0.1,
                        f"{prefix}_p_perm": 0.5,
                        f"{prefix}_n": 200,
                    }
                )
            if hp_confound:
                row.update(
                    {
                        "hp_exact_v": 0.4,
                        "hp_exact_p_perm": 0.01,
                        "hp_exact_aligned": True,
                    }
                )
            if eleven_group_site and sample == "H2009":
                row["cluster_sizes"] = json.dumps(
                    {f"g{group}": 10 for group in range(1, 12)}
                )
                for prefix in (
                    "hp_exact",
                    "hp_family",
                    "strand",
                    "start",
                    "end",
                    "length",
                    "mapq",
                    "cpg_called",
                ):
                    row[f"{prefix}_n"] = 110
        rows.append(row)
    return rows


def latest_tag_audit(rows: list[dict]) -> dict:
    status_counts = Counter(str(row["latest_tag_join_status"]) for row in rows)
    n_sites = len(rows)
    reads = sum(int(row["n_reads_total"]) for row in rows)
    joined = sum(int(row["latest_tag_reads_joined"]) for row in rows)
    multimatch = sum(int(row["latest_tag_projection_multimatch_reads"]) for row in rows)
    all_sites_pass = status_counts == Counter({"PASS": n_sites})
    every_joined = reads == joined
    return {
        "authoritative_tag_source": "same_run_LongPhase_S_external_HP_PS_sidecar",
        "embedded_reads_tsv_hp_used_for_analysis": False,
        "join_occurs_before_focal_ALT_selection": True,
        "counting_unit": "site_read_observation_not_globally_unique_read",
        "n_sites": n_sites,
        "join_status_counts": dict(status_counts),
        "all_sites_pass": all_sites_pass,
        "n_reads_tsv_site_rows": reads,
        "n_exact_hp_ps_site_read_joins": joined,
        "every_reads_tsv_row_joined": every_joined,
        "n_ps_present_site_read_joins": sum(
            int(row["latest_tag_ps_present"]) for row in rows
        ),
        "n_source_hp_replaced_site_read_joins": sum(
            int(row["source_hp_replaced_reads"]) for row in rows
        ),
        "n_sidecar_rows_fetched_site_observations": sum(
            int(row["latest_tag_rows_fetched"]) for row in rows
        ),
        "n_sidecar_rows_eligible_site_observations": sum(
            int(row["latest_tag_rows_eligible"]) for row in rows
        ),
        "n_projection_multimatch_site_reads": multimatch,
        "all_projection_identities_unique": multimatch == 0,
        "n_sites_with_zero_reads_tsv_rows": sum(
            int(row["n_reads_total"]) == 0 for row in rows
        ),
        "pass": all_sites_pass and every_joined and multimatch == 0,
    }


def screen_summary_stratum(rows: list[dict]) -> dict:
    return {
        "n_sites": len(rows),
        "n_stable_null_multigroup": sum(bool(row["stable_null_multigroup"]) for row in rows),
        "latest_hp_ps_terminal_join_audit": latest_tag_audit(rows),
    }


def manifest_from_screen(rows: list[dict]) -> dict:
    entries = []
    totals = Counter()
    for sample in MODULE.DATASETS:
        selected = [row for row in rows if row["sample"] == sample]
        truth = Counter(row["truth_label"] for row in selected)
        counts = {
            "all_ssnv": len(selected),
            "truth_tp": truth["TP"],
            "truth_fp": truth["FP"],
            "truth_unassessed": truth["UNASSESSED"],
        }
        totals.update(counts)
        entries.append(
            {
                "sample": sample,
                "biological_id": MODULE.BIOLOGICAL_IDS[sample],
                "counts": counts,
            }
        )
    return {
        "schema_name": "intersubmod.all_ssnv_focal_alt_input_manifest",
        "schema_version": "1.0.0",
        "pass": True,
        "totals": dict(totals),
        "samples": entries,
    }


def assignment_row(row: dict) -> dict:
    return {
        "schema_name": MODULE.ASSIGNMENT_SCHEMA,
        "schema_version": "1.0.0",
        "strict_confirm_candidate": True,
        "sample": row["sample"],
        "chrom": row["chrom"],
        "pos": row["pos"],
        "posthoc": {
            "ref": row["ref"],
            "alt": row["alt"],
            "biological_id": row["biological_id"],
            "truth_label": row["truth_label"],
        },
    }


def cooccurrence_rows(stable: list[dict], with_candidates: bool) -> list[dict]:
    rows = []
    for source in stable:
        selected = with_candidates and source["sample"] in {"HCC1395", "HCC1395_DORADO"}
        confirmed = 2 if selected else 0
        positions = [120, 150] if selected else []
        m2_gate = MODULE.M2_GATE.evaluate_m2_screen(source)
        expected_projections = int(source["n_reads_total"])
        payload_seed = f"{source['sample']}:{source['chrom']}:{source['pos']}".encode()
        rows.append(
            {
                "sample": source["sample"],
                "biological_id": source["biological_id"],
                "chrom": source["chrom"],
                "pos": source["pos"],
                "ref": source["ref"],
                "alt": source["alt"],
                "truth_label": source["truth_label"],
                "ssnv_branch": source["ssnv_branch"],
                "component_id": source["component_id"],
                "selected_group_id": source["selected_group_id"],
                "stable_null_multigroup": source["stable_null_multigroup"],
                "modal_assignment_ari_min": source["modal_assignment_ari_min"],
                "hp_axis_confound": source["hp_axis_confound"],
                "technical_axis_confound": source["technical_axis_confound"],
                "residual_unexplained_multigroup": source[
                    "residual_unexplained_multigroup"
                ],
                "phase_anchored_robust_epigenetic_candidate": source[
                    "phase_anchored_robust_epigenetic_candidate"
                ],
                "n_all_focal_ref_alt_reads": expected_projections,
                "raw_identity_expected_projections": expected_projections,
                "raw_identity_matched_records": expected_projections,
                "raw_identity_duplicate_projections_collapsed": 0,
                "raw_identity_duplicate_extra_records_collapsed": 0,
                "raw_identity_exact_duplicate_projections_collapsed": 0,
                "raw_identity_rg_only_duplicate_projections_collapsed": 0,
                "raw_identity_duplicate_projection_sha256": MODULE.projection_digest([]),
                "raw_identity_alignment_payload_sha256": hashlib.sha256(
                    b"alignment:" + payload_seed
                ).hexdigest(),
                "raw_identity_recovered_payload_sha256": hashlib.sha256(
                    b"recovered:" + payload_seed
                ).hexdigest(),
                "raw_identity_analysis_site_payload_sha256": hashlib.sha256(
                    b"site:" + payload_seed
                ).hexdigest(),
                "m2_screen_gate_contract": m2_gate["contract"],
                "m2_screen_evaluable": m2_gate["evaluable"],
                "m2_screen_eligible": m2_gate["eligible"],
                "m2_screen_eligibility_status": m2_gate["status"],
                "m2_categorical_level_counts": source[
                    "m2_categorical_level_counts"
                ],
                "m2_axis_statuses": json.dumps(m2_gate["axis_statuses"]),
                "m2_indeterminate_axes": json.dumps(m2_gate["indeterminate_axes"]),
                "m2_low_power_axes": json.dumps(m2_gate["low_power_axes"]),
                "m2_aligned_axes": json.dumps(m2_gate["aligned_axes"]),
                "m2_constant_axes": json.dumps(m2_gate["constant_axes"]),
                "m2_aligned_below_negative_evaluability_power_axes": json.dumps(
                    m2_gate["aligned_below_negative_evaluability_power_axes"]
                ),
                "n_partner_markers": confirmed,
                "n_pair_rows_reconciled": confirmed,
                "pair_row_count_reconciliation_pass": True,
                "n_endpoint_a_testable_markers": confirmed,
                "n_endpoint_a_exact_identifiable_markers": confirmed,
                "n_endpoint_a_exact_not_identifiable_markers": 0,
                "n_endpoint_a_conditional_permutable_markers": confirmed,
                "n_pair_bh_discoveries": confirmed,
                "pair_bh_discovery_positions": json.dumps(positions),
                "n_pair_by_discoveries": confirmed,
                "pair_by_discovery_positions": json.dumps(positions),
                "n_pair_by_confirmed": confirmed,
                "pair_by_confirmed_positions": json.dumps(positions),
                "n_spatially_separated_pair_by_20bp": confirmed,
                "spatially_separated_pair_by_positions_20bp": json.dumps(positions),
                "top_marker_positions": json.dumps(positions),
                "n_top_marker_pair_by_confirmed": confirmed,
                "top_marker_pair_by_confirmed_positions": json.dumps(positions),
                "joint_signature_complete_marker_support": json.dumps(
                    [
                        {
                            "position": position,
                            "testable": True,
                            "reason": "TESTABLE",
                            "n_informative": 12,
                            "cramers_v": 0.7,
                            "delta_alt_fraction": 0.8,
                            "effect_gate_pass": True,
                        }
                        for position in positions
                    ]
                ),
                "joint_signature_n_complete_reads": 12 if selected else 0,
                "joint_signature_n_complete_marker_effect_supported": confirmed,
                "joint_signature_complete_marker_effect_supported_positions": json.dumps(
                    positions
                ),
                "joint_signature_testable": selected,
                "joint_signature_p_conditional_perm": 0.01 if selected else "",
                "joint_signature_permutations": 999 if selected else 0,
                "joint_signature_permutable": selected,
                "joint_signature_conditional_status": (
                    "PERMUTABLE"
                    if selected
                    else "NOT_IDENTIFIABLE_JOINT_SIGNATURE_NOT_TESTABLE"
                ),
                "joint_signature_sensitivity_pass": selected,
                "joint_signature_global_fdr_family_status": (
                    "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
                    if selected
                    else "ELIGIBLE_M2_JOINT_SIGNATURE_NOT_TESTABLE"
                    if m2_gate["eligible"]
                    else "INELIGIBLE_M2_SCREEN"
                ),
                "joint_signature_q_global_bh": 0.01 if selected else "",
                "joint_signature_q_global_by": 0.015 if selected else "",
                "joint_signature_global_bh_discovery": selected,
                "joint_signature_global_by_discovery": selected,
                "n_same_complete_read_effect_supported_top_pair_by": confirmed,
                "same_complete_read_effect_supported_top_pair_by_positions": json.dumps(
                    positions
                ),
                "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp": confirmed,
                "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp": json.dumps(
                    positions
                ),
                MODULE.FORMAL_SELECTION_COLUMN: selected,
                "n_endpoint_a_pre_candidates": confirmed,
                "n_confirmed_markers": confirmed,
                "confirmed_marker_positions": json.dumps(positions),
                "n_independent_confirmed_markers_20bp": confirmed,
                MODULE.DEFAULT_SELECTION_COLUMN: selected,
                "n_confirmed_markers_by_sensitivity": confirmed,
                "confirmed_marker_positions_by_sensitivity": json.dumps(positions),
                "n_independent_confirmed_markers_20bp_by_sensitivity": confirmed,
                MODULE.BY_SELECTION_COLUMN: selected,
                "alt_readset_sha256": hashlib.sha256(
                    b"alt-readset:" + payload_seed
                ).hexdigest(),
                "alt_readset_representative": True,
            }
        )
    return [
        {
            column: row.get(column, "")
            for column in MODULE.RELEASE_FINALIZER.ANALYZER.SITE_COLUMNS
        }
        for row in sorted(
            rows,
            key=lambda value: (
                value["sample"],
                value["chrom"],
                int(value["pos"]),
            ),
        )
    ]


def raw_identity_fixture_audits(
    rows: list[dict],
) -> tuple[dict, dict]:
    aggregate = Counter()
    site_weighted_digest = hashlib.sha256()
    for row in sorted(rows, key=lambda value: (value["sample"], value["chrom"], int(value["pos"]))):
        expected = int(row["raw_identity_expected_projections"])
        matched = int(row["raw_identity_matched_records"])
        duplicates = int(row["raw_identity_duplicate_projections_collapsed"])
        extras = int(row["raw_identity_duplicate_extra_records_collapsed"])
        exact = int(row["raw_identity_exact_duplicate_projections_collapsed"])
        rg_only = int(row["raw_identity_rg_only_duplicate_projections_collapsed"])
        payload = [
            row["sample"],
            row["chrom"],
            int(row["pos"]),
            expected,
            matched,
            duplicates,
            extras,
            exact,
            rg_only,
            row["raw_identity_duplicate_projection_sha256"],
            row["raw_identity_alignment_payload_sha256"],
            row["raw_identity_recovered_payload_sha256"],
            row["raw_identity_analysis_site_payload_sha256"],
        ]
        site_weighted_digest.update(MODULE.compact_json(payload).encode("utf-8"))
        site_weighted_digest.update(b"\n")
        aggregate.update(
            {
                "site_tasks": 1,
                "expected_projection_occurrences": expected,
                "matched_record_occurrences": matched,
                "sites_with_collapsed_duplicates": int(duplicates > 0),
                "duplicate_projection_occurrences_collapsed": duplicates,
                "duplicate_extra_record_occurrences_collapsed": extras,
                "exact_duplicate_projection_occurrences_collapsed": exact,
                "rg_only_duplicate_projection_occurrences_collapsed": rg_only,
            }
        )
    digest = site_weighted_digest.hexdigest()
    summary_audit = {
        **{f"n_{key}": int(value) for key, value in aggregate.items()},
        "equivalence_contract": MODULE.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
        "analysis_scope_identity_contract": MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT,
        "allowed_differing_auxiliary_tags": ["RG"],
        "site_weighted_audit_sha256": digest,
        "all_site_results_passed_invariant_validation": True,
        "missing_projection_policy": MODULE.RAW_IDENTITY_MISSING_POLICY,
        "conflicting_analysis_payload_policy": MODULE.RAW_IDENTITY_CONFLICT_POLICY,
        "failure_counts_materialized": False,
    }
    preflight_audit = {
        "aggregate": {key: int(value) for key, value in aggregate.items()},
        "equivalence_contract": MODULE.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
        "analysis_scope_identity_contract": MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT,
        "allowed_differing_auxiliary_tags": ["RG"],
        "site_weighted_audit_sha256": digest,
        "all_result_rows_passed_invariant_validation": True,
        "missing_projection_policy": MODULE.RAW_IDENTITY_MISSING_POLICY,
        "conflicting_analysis_payload_policy": MODULE.RAW_IDENTITY_CONFLICT_POLICY,
        "failure_counts_materialized": False,
    }
    return summary_audit, preflight_audit


def pair_rows(with_candidates: bool) -> list[dict]:
    if not with_candidates:
        return []
    state_counts = {"RR": 4, "AR": 4, "RA": 0, "AA": 203, "O": 1, "X": 0}
    four_state = MODULE.recompute_four_state_from_counts(state_counts)
    by_q = MODULE.benjamini_yekutieli([0.001] * 4)[0]
    rows = []
    for sample in ("HCC1395", "HCC1395_DORADO"):
        for partner_pos in (120, 150):
            replicated = partner_pos == 120
            rows.append(
                {
                    "sample": sample,
                    "focal_chrom": "chr1",
                    "focal_pos": 100,
                    "focal_ref": "A",
                    "focal_alt": "C",
                    "partner_chrom": "chr1",
                    "partner_pos": partner_pos,
                    "partner_ref": "G",
                    "partner_alt": "T",
                    "distance_bp": partner_pos - 100,
                    "partner_universe_contract": (
                        "all_frozen_latest_LongPhase-S_PASS_biallelic_sSNVs_"
                        "within_focal_plus_or_minus_5000bp"
                    ),
                    "n_all_focal_ref_alt_reads": (
                        3 if sample == "HCC1395" else 5
                    ),
                    "core_partner_call_counts": json.dumps(
                        {"R": 8, "A": 203, "O": 0, "X": 0}
                    ),
                    "all_partner_call_counts": json.dumps(
                        {"R": 8, "A": 203, "O": 1, "X": 0}
                    ),
                    "endpoint_a_testable": True,
                    "endpoint_a_n_informative": 211,
                    "endpoint_a_p_fixed_margins_exact": 0.001,
                    "endpoint_a_exact_state_space_status": "EXACT_ENUMERATED",
                    "endpoint_a_global_fdr_family_status": "ELIGIBLE_M2_EXACT_FAMILY",
                    "endpoint_a_q_global_bh": 0.001,
                    "endpoint_a_q_global_by": by_q,
                    "endpoint_a_cramers_v": 0.7,
                    "endpoint_a_delta_alt_fraction": 0.8,
                    "endpoint_a_effect_gate_pass": True,
                    "endpoint_a_exact_bh_discovery": True,
                    "endpoint_a_exact_by_discovery": True,
                    "endpoint_a_p_conditional_perm": 0.01,
                    "endpoint_a_permutations": 999,
                    "endpoint_a_permutable": True,
                    "endpoint_a_conditional_status": "PERMUTABLE",
                    "endpoint_a_conditional_sensitivity_pass": True,
                    "endpoint_a_formal_pair_by_confirmed": True,
                    "endpoint_a_confirmed_association": True,
                    "endpoint_a_confirmed_by_sensitivity": True,
                    "callability_testable": False,
                    "callability_p_analytic": "",
                    "callability_q_global_bh": "",
                    "callability_q_global_by": "",
                    "callability_cramers_v": "",
                    "callability_noncallable_core_reads": 0,
                    "callability_gate_status": "PASS_ALL_CORE_READS_CALLABLE",
                    "callability_gate_pass": True,
                    "endpoint_b_status": four_state["status"],
                    "endpoint_b_n_joint": four_state["n_joint"],
                    "endpoint_b_n_called_depth": four_state["n_called_depth"],
                    "endpoint_b_state_counts": json.dumps(state_counts),
                    "endpoint_b_n_focal_ref": four_state["n_focal_ref"],
                    "endpoint_b_n_focal_alt": four_state["n_focal_alt"],
                    "endpoint_b_n_partner_ref": four_state["n_partner_ref"],
                    "endpoint_b_n_partner_alt": four_state["n_partner_alt"],
                    "endpoint_b_focal_ancestor_violation_rate": four_state[
                        "focal_ancestor_violation_rate"
                    ],
                    "endpoint_b_partner_ancestor_violation_rate": four_state[
                        "partner_ancestor_violation_rate"
                    ],
                    "endpoint_b_error_ceiling": four_state["error_ceiling"],
                    "endpoint_b_error_model_confidence": four_state["confidence"],
                    "endpoint_b_familywise_confidence": four_state[
                        "familywise_confidence"
                    ],
                    "endpoint_b_relation_family_size": four_state["relation_family_size"],
                    "endpoint_b_multiplicity_method": four_state["multiplicity_method"],
                    "endpoint_b_minimum_zero_violation_depth": four_state[
                        "minimum_zero_violation_depth"
                    ],
                    "endpoint_b_focal_ancestor_violation_p_exact": four_state[
                        "focal_ancestor_violation"
                    ]["p_exact_greater"],
                    "endpoint_b_focal_ancestor_violation_upper_bound": four_state[
                        "focal_ancestor_violation"
                    ]["upper_bound"],
                    "endpoint_b_focal_ancestor_violation_threshold": four_state[
                        "focal_ancestor_violation"
                    ]["threshold"],
                    "endpoint_b_focal_ancestor_violation_status": four_state[
                        "focal_ancestor_violation"
                    ]["status"],
                    "endpoint_b_partner_ancestor_violation_p_exact": four_state[
                        "partner_ancestor_violation"
                    ]["p_exact_greater"],
                    "endpoint_b_partner_ancestor_violation_upper_bound": four_state[
                        "partner_ancestor_violation"
                    ]["upper_bound"],
                    "endpoint_b_partner_ancestor_violation_threshold": four_state[
                        "partner_ancestor_violation"
                    ]["threshold"],
                    "endpoint_b_partner_ancestor_violation_status": four_state[
                        "partner_ancestor_violation"
                    ]["status"],
                    "endpoint_b_branching_violation_p_exact": four_state[
                        "branching_violation"
                    ]["p_exact_greater"],
                    "endpoint_b_branching_violation_upper_bound": four_state[
                        "branching_violation"
                    ]["upper_bound"],
                    "endpoint_b_branching_violation_threshold": four_state[
                        "branching_violation"
                    ]["threshold"],
                    "endpoint_b_branching_violation_status": four_state[
                        "branching_violation"
                    ]["status"],
                    "endpoint_b_complete_four_state_testable": four_state[
                        "complete_four_state_testable"
                    ],
                    "endpoint_b_relation_compatibility": four_state[
                        "relation_compatibility"
                    ],
                    "endpoint_b_compatible_relation_models": json.dumps(
                        four_state["compatible_relation_models"]
                    ),
                    "endpoint_b_n_compatible_relation_models": four_state[
                        "n_compatible_relation_models"
                    ],
                    "top_coverage_marker": True,
                    "focal_truth_label": "TP",
                    "partner_truth_label": "TP",
                    "focal_ssnv_branch": "RECOVERED_PASS",
                    "partner_ssnv_branch": "RECOVERED_PASS",
                    "focal_component_id": "component-1",
                    "partner_component_id": "component-1",
                    "topology_scope": "DIRECT_RETAINED_SHARED_REGION",
                    "topology_region": "chr1:1-1000",
                    "topology_order_status": "FOCAL_BEFORE_PARTNER",
                    "cross_platform_exact_pair_present": True,
                    "cross_platform_biological_n": 1,
                    "cross_platform_conditional_by_effect_direction_replication_pass": replicated,
                    "cross_platform_replication_status": (
                        "DISTINCT_MOLECULE_SET_TECHNICAL_CONDITIONAL_BY_CONCORDANCE"
                        if replicated
                        else "EXACT_PAIR_EFFECT_NOT_COMPATIBLE"
                    ),
                    "cross_platform_core_read_name_overlap_present": False,
                    "cross_platform_molecule_independence_status": (
                        "DISTINCT_CORE_READ_NAME_SETS_TECHNICAL_ONLY_BIOLOGICAL_N1"
                    ),
                }
            )
    return [
        {
            column: row.get(column, "")
            for column in MODULE.RELEASE_FINALIZER.ANALYZER.PAIR_COLUMNS
        }
        for row in rows
    ]


def set_four_state_fields(row: dict, state_counts: dict[str, int]) -> None:
    summary = MODULE.recompute_four_state_from_counts(state_counts)
    row.update(
        {
            "endpoint_b_status": summary["status"],
            "endpoint_b_n_joint": summary["n_joint"],
            "endpoint_b_n_called_depth": summary["n_called_depth"],
            "endpoint_b_state_counts": json.dumps(state_counts),
            "endpoint_b_n_focal_ref": summary["n_focal_ref"],
            "endpoint_b_n_focal_alt": summary["n_focal_alt"],
            "endpoint_b_n_partner_ref": summary["n_partner_ref"],
            "endpoint_b_n_partner_alt": summary["n_partner_alt"],
            "endpoint_b_focal_ancestor_violation_rate": summary[
                "focal_ancestor_violation_rate"
            ],
            "endpoint_b_partner_ancestor_violation_rate": summary[
                "partner_ancestor_violation_rate"
            ],
            "endpoint_b_error_ceiling": summary["error_ceiling"],
            "endpoint_b_error_model_confidence": summary["confidence"],
            "endpoint_b_minimum_zero_violation_depth": summary[
                "minimum_zero_violation_depth"
            ],
            "endpoint_b_focal_ancestor_violation_p_exact": summary[
                "focal_ancestor_violation"
            ]["p_exact_greater"],
            "endpoint_b_focal_ancestor_violation_upper_bound": summary[
                "focal_ancestor_violation"
            ]["upper_bound"],
            "endpoint_b_focal_ancestor_violation_threshold": summary[
                "focal_ancestor_violation"
            ]["threshold"],
            "endpoint_b_focal_ancestor_violation_status": summary[
                "focal_ancestor_violation"
            ]["status"],
            "endpoint_b_partner_ancestor_violation_p_exact": summary[
                "partner_ancestor_violation"
            ]["p_exact_greater"],
            "endpoint_b_partner_ancestor_violation_upper_bound": summary[
                "partner_ancestor_violation"
            ]["upper_bound"],
            "endpoint_b_partner_ancestor_violation_threshold": summary[
                "partner_ancestor_violation"
            ]["threshold"],
            "endpoint_b_partner_ancestor_violation_status": summary[
                "partner_ancestor_violation"
            ]["status"],
            "endpoint_b_branching_violation_p_exact": summary["branching_violation"][
                "p_exact_greater"
            ],
            "endpoint_b_branching_violation_upper_bound": summary["branching_violation"][
                "upper_bound"
            ],
            "endpoint_b_branching_violation_threshold": summary["branching_violation"][
                "threshold"
            ],
            "endpoint_b_branching_violation_status": summary["branching_violation"][
                "status"
            ],
            "endpoint_b_complete_four_state_testable": summary[
                "complete_four_state_testable"
            ],
            "endpoint_b_relation_compatibility": summary["relation_compatibility"],
            "endpoint_b_compatible_relation_models": json.dumps(
                summary["compatible_relation_models"]
            ),
            "endpoint_b_n_compatible_relation_models": summary[
                "n_compatible_relation_models"
            ],
        }
    )


def strict_rows(selected: list[dict]) -> list[dict]:
    rows = []
    for source in selected:
        robustness = source["sample"] == "HCC1395"
        # Deliberately reverse descriptive q and robustness to prove q is not a gate.
        descriptive_q = 0.99 if robustness else 0.001
        rows.append(
            {
                **source,
                "strict_selection_column": MODULE.STRICT_SELECTION_COLUMN,
                "strict_analysis_class": "FORMAL",
                "strict_formal_parameter_gate": True,
                "strict_formal_selection_contract_gate": True,
                "strict_formal_selection_column": MODULE.STRICT_SELECTION_COLUMN,
                "strict_candidate_selection_gate": True,
                "strict_cooccurrence_receipt_gate": True,
                "strict_artifact_identity_gate": True,
                "strict_screen_artifact_hash_lock_gate": True,
                "strict_postselection_fdr_calibrated": False,
                "strict_null_robustness_pass": robustness,
                "strict_methyl_partition_robustness_status": (
                    "R1_STRICT_METHYL_PARTITION_ROBUSTNESS_PASS"
                    if robustness
                    else "NOT_EVALUABLE"
                ),
                "strict_scientific_status": (
                    "ROBUSTNESS_PASS_NOT_FDR_CALIBRATED"
                    if robustness
                    else "NOT_EVALUABLE"
                ),
                "strict_not_evaluable_reason": "" if robustness else "insufficient_null",
                "strict_combined_empirical_p_postselection_descriptive": descriptive_q,
                "strict_postselection_bh_q_descriptive": descriptive_q,
                "strict_postselection_by_q_descriptive": descriptive_q,
                "strict_assignment_concordance_ari_min": 1.0 if robustness else 0.4,
                "strict_failure_reason": "PASS" if robustness else "insufficient_null",
            }
        )
    return rows


def make_fixture(
    tmp_path: Path,
    *,
    with_candidates: bool = True,
    eleven_group_site: bool = False,
) -> MODULE.InputBundle:
    source_authority = MODULE.current_release_source_authority()
    root = tmp_path / "fixture"
    screen_dir = root / "screen"
    cooccurrence_dir = root / "cooccurrence"
    strict_dir = root / "strict"
    tumor_ref_dir = root / "tumor_ref"
    normal_dir = root / "normal" if with_candidates else None
    for directory in (screen_dir, cooccurrence_dir, strict_dir, tumor_ref_dir):
        directory.mkdir(parents=True, exist_ok=False)
    if normal_dir is not None:
        normal_dir.mkdir(parents=True, exist_ok=False)

    all_screen = screen_rows(eleven_group_site=eleven_group_site)
    manifest_path = root / "manifest.json"
    write_json(manifest_path, manifest_from_screen(all_screen))
    screen_sites = screen_dir / "all_ssnv_site_results.tsv.gz"
    write_tsv(screen_sites, all_screen)
    stable = [row for row in all_screen if row["stable_null_multigroup"]]
    assignment_path = screen_dir / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    with gzip.open(assignment_path, "wt", encoding="utf-8") as handle:
        for row in stable:
            handle.write(json.dumps(assignment_row(row)) + "\n")
    screen_summary = screen_dir / "all_ssnv_summary.json"
    pooled_tag_audit = latest_tag_audit(all_screen)
    clustering_contract = {
        "m1_stability_gate_contract": "coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8",
        "prior_screen_thresholds": {"modal_assignment_ari_min": 0.8},
        "stable_null_multigroup_basis": (
            "coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8"
        ),
    }
    write_json(
        screen_summary,
        {
            "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_screen",
            "schema_version": "1.2.0",
            "status": "EXECUTION_PASS",
            "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
            "scope": {
                "full_469849": True,
                "selected_datasets": list(MODULE.DATASETS),
                "selected_samples": list(MODULE.DATASETS),
                "expected_sites": len(all_screen),
                "processed_sites": len(all_screen),
            },
            "clustering_contract": clustering_contract,
            "pooled_site_weighted": screen_summary_stratum(all_screen),
            "per_dataset": {
                sample: screen_summary_stratum(
                    [row for row in all_screen if row["sample"] == sample]
                )
                for sample in MODULE.DATASETS
            },
            "per_sample": {
                sample: screen_summary_stratum(
                    [row for row in all_screen if row["sample"] == sample]
                )
                for sample in MODULE.DATASETS
            },
            "posthoc_truth_strata": {
                truth: screen_summary_stratum(
                    [row for row in all_screen if row["truth_label"] == truth]
                )
                for truth in MODULE.TRUTH_LABELS
            },
            "posthoc_biological_id_strata": {
                biological_id: screen_summary_stratum(
                    [
                        row
                        for row in all_screen
                        if row["biological_id"] == biological_id
                    ]
                )
                for biological_id in sorted(set(MODULE.BIOLOGICAL_IDS.values()))
            },
            "posthoc_ledger_branch_strata": {
                branch: screen_summary_stratum(
                    [row for row in all_screen if row["ssnv_branch"] == branch]
                )
                for branch in sorted({row["ssnv_branch"] for row in all_screen})
            },
            "n_stable_assignment_records": len(stable),
            "latest_hp_ps_terminal_join_audit": pooled_tag_audit,
            "pass": True,
        },
    )
    write_json(
        screen_dir / "run_manifest.json",
        {
            "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
            "schema_version": "1.2.0",
            "status": "EXECUTION_PASS",
            "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
            "input_manifest": artifact(manifest_path),
            "execution": {
                "selected_datasets": list(MODULE.DATASETS),
                "selected_samples": list(MODULE.DATASETS),
            },
            "contracts": {
                "truth_and_cooccurrence_enter_clustering": False,
                "screen_global_fdr_calibrated": False,
                "strict_methyl_partition_robustness_status": (
                    "NOT_EVALUATED_AT_M1_SCREEN"
                ),
                "strict_confirm_status_legacy_alias": "NOT_RUN",
                "strict_confirm_candidate_is_formal_r1_claim": False,
                "existing_results_overwritten": False,
                **clustering_contract,
                "latest_hp_ps_terminal_join": pooled_tag_audit,
            },
            "outputs": {
                "site_results": artifact(screen_sites),
                "stable_assignments": artifact(assignment_path),
                "summary": artifact(screen_summary),
            },
            "counts": {
                "expected_sites": len(all_screen),
                "processed_sites": len(all_screen),
                "stable_assignment_records": len(stable),
                "reads_tsv_site_rows": pooled_tag_audit["n_reads_tsv_site_rows"],
                "exact_hp_ps_site_read_joins": pooled_tag_audit[
                    "n_exact_hp_ps_site_read_joins"
                ],
                "ps_present_site_read_joins": pooled_tag_audit[
                    "n_ps_present_site_read_joins"
                ],
                "source_hp_replaced_site_read_joins": pooled_tag_audit[
                    "n_source_hp_replaced_site_read_joins"
                ],
            },
            "pass": True,
        },
    )

    coocc_sites_rows = cooccurrence_rows(stable, with_candidates)
    raw_summary_audit, raw_preflight_audit = raw_identity_fixture_audits(
        coocc_sites_rows
    )
    coocc_sites = cooccurrence_dir / "methyl_ssnv_site_results.tsv.gz"
    write_tsv(coocc_sites, coocc_sites_rows)
    pairs = pair_rows(with_candidates)
    pair_fields = list(pair_rows(True)[0])
    pair_path = cooccurrence_dir / "methyl_ssnv_pair_results.tsv.gz"
    write_tsv(pair_path, pairs, pair_fields)
    raw_duplicate_path = cooccurrence_dir / "raw_identity_duplicate_audit.tsv.gz"
    write_tsv(
        raw_duplicate_path,
        [],
        fields=list(MODULE.RAW_IDENTITY_DUPLICATE_COLUMNS),
    )
    oracle_path = cooccurrence_dir / "oracle_cases.json"
    write_json(
        oracle_path,
        {
            "schema_name": (
                "intersubmod.methyl_ssnv_cooccurrence_validation.oracle_cases"
            ),
            "schema_version": MODULE.COOCCURRENCE_SCHEMA_VERSION,
            "partner_window_bp": MODULE.RELEASE_FINALIZER.ANALYZER.PAIR_WINDOW_BP,
            "focal_cases": [
                {
                    "case_id": case_id,
                    "case_type": "focal_partner_window",
                    "present": True,
                    "expected_partner_positions": [],
                    "observed_partner_positions": [],
                    "partner_window_oracle_pass": True,
                }
                for case_id in sorted(MODULE.COOCCURRENCE_ORACLE_FOCAL_IDS)
            ],
            "shared_readset_case": {
                "case_id": MODULE.COOCCURRENCE_SHARED_READSET_ORACLE_ID,
                "case_type": "shared_alt_readset_dedup",
                "positions": [1, 2],
                "present_positions": [1, 2],
                "same_alt_readset": True,
                "one_alt_readset_representative": True,
            },
        },
    )
    frozen_identity_artifacts = []
    for sample in MODULE.DATASETS:
        for field_name in sorted(MODULE.MANIFEST_HASH_REQUIRED_FIELDS):
            frozen_identity_artifacts.append(
                {
                    "sample": sample,
                    "field": field_name,
                    "path": f"/fixture/{sample}/{field_name}",
                    "identity_mode": "sha256_size_path_v1",
                    "size_bytes": 1,
                    "mtime_ns": 1,
                    "sha256": "a" * 64,
                }
            )
        for field_name in sorted(MODULE.MANIFEST_LARGE_METADATA_IDENTITY_FIELDS):
            frozen_identity_artifacts.append(
                {
                    "sample": sample,
                    "field": field_name,
                    "path": f"/fixture/{sample}/{field_name}",
                    "identity_mode": "explicit_large_artifact_size_mtime_path_v1",
                    "size_bytes": 1,
                    "mtime_ns": 1,
                    "sha256": None,
                }
            )
    frozen_identity_policy = {
        "hash_required_fields": sorted(MODULE.MANIFEST_HASH_REQUIRED_FIELDS),
        "explicit_large_metadata_identity_fields": sorted(
            MODULE.MANIFEST_LARGE_METADATA_IDENTITY_FIELDS
        ),
        "n_sha256_identity_artifacts": 7
        * len(MODULE.MANIFEST_HASH_REQUIRED_FIELDS),
        "n_explicit_large_size_mtime_identity_artifacts": 7
        * len(MODULE.MANIFEST_LARGE_METADATA_IDENTITY_FIELDS),
        "artifacts": frozen_identity_artifacts,
        "claim_guardrail": (
            "Only explicitly listed large roles use size+mtime; this is weaker than full "
            "SHA-256 and is not bytewise identity."
        ),
    }
    m2_provenance = {
        "effect_and_permutation_p_source": (
            "frozen source-locked screen-producer site rows"
        ),
        "downstream_raw_read_axis_statistic_recomputed": False,
        "downstream_axis_classification_recomputed": True,
        "downstream_validation": [
            "axis sample-size reconciliation",
            "499-permutation add-one p-value grid",
            "effect threshold classification",
            "80-percent planning-power evaluability",
        ],
        "claim_guardrail": (
            "The cooccurrence stage independently reclassifies producer-derived axis "
            "statistics; it is not an independent raw-read remeasurement."
        ),
    }
    coocc_summary = cooccurrence_dir / "summary.json"
    write_json(
        coocc_summary,
        {
            "schema_name": "intersubmod.methyl_ssnv_cooccurrence_validation.summary",
            "schema_version": "2.0.0",
            "task_type": "B_comprehensive_validation",
            "samples": list(MODULE.DATASETS),
            "n_stable_sites": len(stable),
            "n_focal_partner_pairs": len(pairs),
            "n_endpoint_a_testable_pairs": len(pairs),
            "n_pair_by_confirmed": len(pairs),
            "n_multi_marker_molecular_haplotype_base_candidates": sum(
                row[MODULE.FORMAL_SELECTION_COLUMN] for row in coocc_sites_rows
            ),
            "site_pair_count_reconciliation": {
                "pass": True,
                "n_sites_reconciled": len(stable),
                "n_pair_rows_reconciled": len(pairs),
            },
            "cross_platform_replication": {
                "n_exact_pairs_present_both": len(pairs) // 2,
                "n_exact_pairs_confirmed_both": len(pairs) // 2,
                "n_exact_pairs_conditional_by_effect_direction_replicated": (
                    1 if with_candidates else 0
                ),
            },
            "frozen_input_identity_policy": frozen_identity_policy,
            "m2_axis_statistic_provenance": m2_provenance,
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
            "source_authority": source_authority,
            "raw_bam_identity_recovery_audit": raw_summary_audit,
            "pass": True,
        },
    )
    cooccurrence_code = {
        role: artifact(path)
        for role, path in MODULE.COOCCURRENCE_CODE_PATHS.items()
    }
    cooccurrence_source_modes = {role: "0o444" for role in cooccurrence_code}
    write_json(
        cooccurrence_dir / "run_receipt.json",
        {
            "schema_name": "intersubmod.methyl_ssnv_cooccurrence_validation.run_receipt",
            "schema_version": "2.0.0",
            "started_at_utc": "2026-07-15T00:01:00+00:00",
            "finished_at_utc": "2026-07-15T00:02:00+00:00",
            "task_type": "B_comprehensive_validation",
            "scope": "all_manifest_samples",
            "full_scope": True,
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
            "inputs": {
                "manifest": artifact(manifest_path),
                "assignments": artifact(assignment_path),
                "sites": artifact(screen_sites),
            },
            "outputs": {
                "site_tsv": artifact(coocc_sites),
                "pair_tsv": artifact(pair_path),
                "summary_json": artifact(coocc_summary),
                "raw_identity_duplicate_audit_tsv": artifact(raw_duplicate_path),
                "case_json": artifact(oracle_path),
            },
            "code": cooccurrence_code,
            "source_lock": {
                "source_identity_before": cooccurrence_code,
                "source_identity_after_compute": cooccurrence_code,
                "source_identity_after_output": cooccurrence_code,
                "source_modes_before": cooccurrence_source_modes,
                "source_modes_after_compute": cooccurrence_source_modes,
                "source_modes_after_output": cooccurrence_source_modes,
                "all_sources_read_only_and_unchanged": True,
            },
            "counts": {
                "stable_sites": len(stable),
                "focal_partner_pairs": len(pairs),
                "pair_by_confirmed": len(pairs),
                "multi_marker_molecular_haplotype_base_candidates": sum(
                    row[MODULE.FORMAL_SELECTION_COLUMN] for row in coocc_sites_rows
                ),
                "raw_identity_expected_projection_occurrences": raw_summary_audit[
                    "n_expected_projection_occurrences"
                ],
                "raw_identity_duplicate_projection_occurrences_collapsed": 0,
                "raw_identity_duplicate_extra_record_occurrences_collapsed": 0,
                "raw_identity_sparse_duplicate_rows": 0,
                "raw_identity_all_site_results_passed_invariant_validation": True,
                "raw_identity_missing_projection_policy": (
                    MODULE.RAW_IDENTITY_MISSING_POLICY
                ),
                "raw_identity_conflicting_analysis_payload_policy": (
                    MODULE.RAW_IDENTITY_CONFLICT_POLICY
                ),
                "raw_identity_failure_counts_materialized": False,
            },
            "frozen_manifest_input_identity_policy": frozen_identity_policy,
            "m2_axis_statistic_provenance": m2_provenance,
            "pass": True,
        },
    )

    selected = [
        row for row in coocc_sites_rows if row[MODULE.FORMAL_SELECTION_COLUMN]
    ]
    if with_candidates:
        strict_site_path = strict_dir / "strict_methyl_candidate_confirmation_sites.tsv.gz"
        strict_site_rows = strict_rows(selected)
        write_tsv(strict_site_path, strict_site_rows)
        robustness_n = sum(row["strict_null_robustness_pass"] for row in strict_site_rows)
        strict_summary = strict_dir / "strict_methyl_candidate_confirmation_summary.json"
        strict_counts = {
            "n_selected_candidates": len(selected),
            "n_methyl_partition_robustness_evaluable": robustness_n,
            "n_methyl_partition_robustness_not_evaluable": len(selected) - robustness_n,
            "n_null_robustness_pass": robustness_n,
            "n_null_robustness_fail_retained": 0,
            "n_exploratory_only": 0,
        }
        write_json(
            strict_summary,
            {
                "schema_name": "intersubmod.strict_methyl_candidate_confirmation_summary",
                "schema_version": "3.1.0",
                "execution_status": "EXECUTION_PASS",
                "status": "EXECUTION_PASS",
                "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
                "selection_contract": {
                    "selection_column": MODULE.STRICT_SELECTION_COLUMN,
                    "formal_selection_column": MODULE.STRICT_SELECTION_COLUMN,
                    "n_selected_candidates": len(selected),
                },
                "postselection_diagnostic_contract": {
                    "fdr_calibrated": False,
                    "bh_by_values_are_descriptive_only": True,
                },
                "counts": strict_counts,
                "pass": True,
            },
        )
        write_json(
            strict_dir / "run_manifest.json",
            {
                "schema_name": "intersubmod.strict_methyl_candidate_confirmation_run_manifest",
                "schema_version": "3.1.0",
                "started_at_utc": "2026-07-15T00:01:00+00:00",
                "finished_at_utc": "2026-07-15T00:02:00+00:00",
                "execution_status": "EXECUTION_PASS",
                "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
                "inputs": {
                    "candidate_table": artifact(coocc_sites),
                    "assignments": artifact(assignment_path),
                    "cooccurrence_receipt": artifact(
                        cooccurrence_dir / "run_receipt.json"
                    ),
                },
                "counts": strict_counts,
                "outputs": {
                    "site_results": artifact(strict_site_path),
                    "summary": artifact(strict_summary),
                },
                "pass": True,
            },
        )
    else:
        write_json(
            strict_dir / "run_receipt.json",
            {
                "schema_name": (
                    "intersubmod.strict_methyl_candidate_confirmation.not_applicable"
                ),
                "schema_version": "3.1.0",
                "started_at_utc": "2026-07-15T00:01:00+00:00",
                "finished_at_utc": "2026-07-15T00:02:00+00:00",
                "execution_status": "NOT_APPLICABLE",
                "status": "NOT_APPLICABLE",
                "reason": "ZERO_SELECTED_CANDIDATES",
                "selection_contract": {
                    "selection_column": MODULE.STRICT_SELECTION_COLUMN,
                    "formal_selection_column": MODULE.STRICT_SELECTION_COLUMN,
                    "n_selected_candidates": 0,
                },
                "counts": {"n_selected_candidates": 0},
                "inputs": {
                    "candidate_table": artifact(coocc_sites),
                    "assignments": artifact(assignment_path),
                    "cooccurrence_receipt": artifact(
                        cooccurrence_dir / "run_receipt.json"
                    ),
                },
                "pass_semantics": (
                    "execution_integrity_only_not_scientific_confirmation"
                ),
                "is_negative_result": False,
                "scientific_interpretation": {
                    "is_negative_result": False,
                    "statement": "Zero selected candidates are not a biological negative.",
                },
                "pass": True,
            },
        )

    tumor_rows = [
        {
            "sample": row["sample"],
            "chrom": row["chrom"],
            "pos": row["pos"],
            "ref": row["ref"],
            "alt": row["alt"],
            "n_tumor_ref": 12,
            "ref_evaluable": True,
            "ref_stable_null_multigroup": False,
            "ref_not_testable_reason": "",
        }
        for row in stable
    ]
    tumor_sites = tumor_ref_dir / "all_ssnv_tumor_ref_control_site_results.tsv.gz"
    write_tsv(tumor_sites, tumor_rows)
    tumor_summary = tumor_ref_dir / "all_ssnv_tumor_ref_control_summary.json"
    write_json(
        tumor_summary,
        {
            "schema_name": "intersubmod.all_ssnv_tumor_ref_controls.summary",
            "schema_version": "1.0.0",
            "task_type": "B_comprehensive_validation",
            "pass": True,
        },
    )
    write_json(
        tumor_ref_dir / "run_manifest.json",
        {
            "schema_name": "intersubmod.all_ssnv_tumor_ref_controls.run_manifest",
            "schema_version": "1.0.0",
            "started_at_utc": "2026-07-15T00:01:00+00:00",
            "finished_at_utc": "2026-07-15T00:02:00+00:00",
            "inputs": {
                "site_results": artifact(screen_sites),
                "stable_assignments": artifact(assignment_path),
            },
            "outputs": {
                "site_results": artifact(tumor_sites),
                "summary": artifact(tumor_summary),
            },
            "counts": {
                "primary_stable_sites": len(stable),
                "processed_sites": len(stable),
            },
            "pass": True,
        },
    )

    if normal_dir is not None:
        paired_receipt_path = root / "paired_run_receipt.json"
        paired_provenance = producer_provenance(
            MODULE.MATCHED_NORMAL_RUNNER, "matched_normal_runner"
        )
        selected_samples = sorted({row["sample"] for row in selected})
        paired_sample_receipts = [
            {
                "sample": sample,
                "validation": {"artifact_set_sha256": format(index + 1, "064x")},
            }
            for index, sample in enumerate(selected_samples)
        ]
        write_json(
            paired_receipt_path,
            {
                "schema_name": "intersubmod.matched_normal_candidate_control_run",
                "schema_version": "1.0.0",
                "started_at_utc": "2026-07-15T00:01:00+00:00",
                "finished_at_utc": "2026-07-15T00:02:00+00:00",
                "command": MODULE.MATCHED_NORMAL_RUNNER.canonical_task_b_command(),
                **paired_provenance,
                "receipts": paired_sample_receipts,
                "pass": True,
            },
        )
        paired_receipt_path.chmod(0o444)
        normal_rows = [
            {
                "sample": row["sample"],
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "normal_called_reads": 12,
                "normal_alt_reads": 0,
                "normal_ref_reads": 12,
                "normal_unknown_reads": 1,
                "normal_focal_callability_gate": True,
                "normal_control_evaluable": True,
                "normal_control_not_evaluable_reason": "",
                "normal_ref_methyl_stable_multigroup": False,
                "normal_ref_methyl_nonreplication_gate": True,
                "normal_control_status": "REF_METHYL_PARTITION_NOT_REPRODUCED",
                "normal_evaluable": True,
                "normal_stable_multigroup": False,
                "normal_genetic_alt_support_present": False,
                "normal_genetic_alt_fraction_called": "",
                "tumor_ref_evaluable": True,
                "tumor_ref_stable_multigroup": False,
                "primary_group_assignment_coverage": 1.0,
                "primary_identity_collision_count": 0,
                "primary_identity_missing_count": 0,
            }
            for row in selected
        ]
        normal_sites = normal_dir / "matched_normal_candidate_controls.tsv.gz"
        write_tsv(normal_sites, normal_rows)
        normal_summary = normal_dir / "matched_normal_candidate_controls_summary.json"
        analyzer_provenance = producer_provenance(
            MODULE.MATCHED_NORMAL_ANALYZER, "matched_normal_analyzer"
        )
        normal_counts = {
            "n_candidates": len(normal_rows),
            "all_primary_group_assignments_exact": True,
        }
        paired_artifact_validation = {
            "samples": {
                receipt["sample"]: {
                    "output_dir": str((root / receipt["sample"]).resolve()),
                    "candidate_sites": 1,
                    "artifacts_verified": 6,
                    "artifact_set_sha256": receipt["validation"]["artifact_set_sha256"],
                    "pass": True,
                }
                for receipt in paired_sample_receipts
            },
            "all_runner_artifact_set_sha256_recomputed": True,
        }
        write_json(
            normal_summary,
            {
                "schema_name": "intersubmod.matched_normal_candidate_control_analysis",
                "schema_version": "2.0.0",
                **analyzer_provenance,
                "inputs": {
                    "primary_assignments": artifact(assignment_path),
                    "paired_run_receipt": artifact(paired_receipt_path),
                },
                "pooled": normal_counts,
                "paired_artifact_identity_validation": paired_artifact_validation,
                "not_evaluable_is_negative_result": False,
                "pass_semantics": (
                    "execution_and_identity_integrity_only_not_background_negativity"
                ),
                "pass": True,
            },
        )
        normal_sites.chmod(0o444)
        normal_summary.chmod(0o444)
        normal_receipt_path = normal_dir / "run_receipt.json"
        write_json(
            normal_receipt_path,
            {
                "schema_name": "intersubmod.matched_normal_candidate_control_analysis_run",
                "schema_version": "2.0.0",
                "started_at_utc": "2026-07-15T00:01:00+00:00",
                "finished_at_utc": "2026-07-15T00:02:00+00:00",
                "command": MODULE.MATCHED_NORMAL_ANALYZER.canonical_task_b_command(),
                **analyzer_provenance,
                "inputs": {
                    "primary_assignments": artifact(assignment_path),
                    "paired_run_receipt": artifact(paired_receipt_path),
                },
                "outputs": {
                    "site_table": artifact(normal_sites),
                    "summary": artifact(normal_summary),
                },
                "counts": normal_counts,
                "paired_artifact_identity_validation": paired_artifact_validation,
                "not_evaluable_is_negative_result": False,
                "pass_semantics": (
                    "execution_and_identity_integrity_only_not_background_negativity"
                ),
                "pass": True,
            },
        )
        normal_receipt_path.chmod(0o444)

    audit_payload = {
        "schema_name": "intersubmod.stable_primary_artifact_audit",
        "schema_version": "2.0.0",
        "created_at_utc": "2026-07-15T00:00:00+00:00",
        "finished_at_utc": "2026-07-15T00:00:00+00:00",
        "started_at_utc": "2026-07-14T23:59:00+00:00",
        "task_type": "B_comprehensive_validation",
        "scope": "complete_primary_stable_null_multigroup_set",
        "source_authority": source_authority,
        "code": MODULE.PRIMARY_AUDITOR.capture_source_identity(),
        "source_lock": {
            "source_identity_before": MODULE.PRIMARY_AUDITOR.capture_source_identity(),
            "source_identity_after_compute": MODULE.PRIMARY_AUDITOR.capture_source_identity(),
            "source_modes_before": MODULE.PRIMARY_AUDITOR.capture_source_modes(),
            "source_modes_after_compute": MODULE.PRIMARY_AUDITOR.capture_source_modes(),
            "all_sources_read_only_and_unchanged": True,
        },
        "execution": dict(MODULE.PRIMARY_AUDITOR.CANONICAL_PARAMETERS),
        "inputs": {
            "site_results": artifact(screen_sites),
            "stable_assignments": artifact(assignment_path),
            "consumer_receipts": [],
        },
        "verification": {
            "stable_site_assignment_key_sets_exact": True,
            "artifact_roles_exact": ["distance_matrix", "methylation_matrix", "reads"],
            "path_size_sha256_verified": True,
            "artifact_set_sha256": "a" * 64,
            "errors": 0,
        },
        "counts": {
            "stable_sites": len(stable),
            "assignment_records": len(stable),
            "primary_artifacts_expected": len(stable) * 3,
            "primary_artifacts_verified": len(stable) * 3,
        },
        "pass_semantics": "execution_integrity_and_primary_artifact_identity_only",
        "scientific_interpretation": {"is_negative_result": False},
        "pass": True,
    }
    audit_pre = root / "primary_artifact_audit.pre.json"
    audit_post = root / "primary_artifact_audit.post.json"
    audit_payload["command"] = MODULE.PRIMARY_AUDITOR.canonical_command(
        site_results=screen_sites,
        assignments=assignment_path,
        consumer_receipts=(),
        output=audit_pre,
        **MODULE.PRIMARY_AUDITOR.CANONICAL_PARAMETERS,
    )
    write_json(audit_pre, audit_payload)
    audit_pre.chmod(0o444)

    cooccurrence_preflight = root / "cooccurrence_task_contract_preflight.json"
    preflight_code = {
        "preflight": artifact(MODULE.COOCCURRENCE_PREFLIGHT_SOURCE),
        **cooccurrence_code,
    }
    preflight_modes = {role: "0o444" for role in preflight_code}
    write_json(
        cooccurrence_preflight,
        {
            "schema_name": "intersubmod.cooccurrence_task_contract_preflight",
            "schema_version": MODULE.COOCCURRENCE_PREFLIGHT_SCHEMA_VERSION,
            "started_at_utc": "2026-07-15T00:00:10+00:00",
            "finished_at_utc": "2026-07-15T00:00:50+00:00",
            "task_type": "B_comprehensive_validation",
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
            "source_authority": source_authority,
            "inputs": {
                "manifest": artifact(manifest_path),
                "assignments": artifact(assignment_path),
                "sites": artifact(screen_sites),
                "primary_artifact_audit_pre": artifact(audit_pre),
            },
            "code": {
                "source_identity_before": preflight_code,
                "source_identity_after": preflight_code,
                "source_modes_before": preflight_modes,
                "source_modes_after": preflight_modes,
            },
            "observed": {
                "task_count": 102_842,
                "raw_bam_identity_recovery": raw_preflight_audit,
            },
            "checks": {"fixture_contract_complete": True},
            "pass": True,
        },
    )
    preflight_payload = json.loads(cooccurrence_preflight.read_text())
    preflight_payload["input_lock"] = {
        "identity_before": preflight_payload["inputs"],
        "identity_after": preflight_payload["inputs"],
        "all_primary_inputs_unchanged": True,
    }
    write_json(cooccurrence_preflight, preflight_payload)

    cooccurrence_receipt_path = cooccurrence_dir / "run_receipt.json"
    cooccurrence_receipt = json.loads(cooccurrence_receipt_path.read_text())
    cooccurrence_receipt["inputs"]["primary_artifact_audit_pre"] = artifact(audit_pre)
    cooccurrence_receipt["inputs"]["preflight_receipt"] = artifact(
        cooccurrence_preflight
    )
    cooccurrence_receipt["preflight_gate"] = {
        "schema_name": "intersubmod.cooccurrence_task_contract_preflight",
        "schema_version": MODULE.COOCCURRENCE_PREFLIGHT_SCHEMA_VERSION,
        "task_count": 102_842,
        "pass": True,
        "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
    }
    cooccurrence_receipt["release_status"] = (
        "PRODUCER_PASS_PENDING_RUNNER_RECONCILIATION"
    )
    cooccurrence_receipt["pass_semantics"] = (
        "execution_integrity_only_pending_release_reconciliation"
    )
    cooccurrence_receipt["source_authority"] = source_authority
    cooccurrence_receipt["input_lock"] = {
        "identity_before": cooccurrence_receipt["inputs"],
        "identity_after_compute": cooccurrence_receipt["inputs"],
        "identity_after_output": cooccurrence_receipt["inputs"],
        "all_primary_inputs_unchanged": True,
    }
    write_json(cooccurrence_receipt_path, cooccurrence_receipt)
    release_receipt_path = cooccurrence_dir / "release_receipt.json"
    release_inputs = {
        "preflight": artifact(cooccurrence_preflight),
        "producer_receipt": artifact(cooccurrence_receipt_path),
        "summary": artifact(coocc_summary),
        "sites": artifact(coocc_sites),
        "pairs": artifact(pair_path),
        "duplicates": artifact(raw_duplicate_path),
        "oracle": artifact(oracle_path),
    }
    release_args = SimpleNamespace(
        preflight=cooccurrence_preflight,
        producer_receipt=cooccurrence_receipt_path,
        summary=coocc_summary,
        sites=coocc_sites,
        pairs=pair_path,
        duplicates=raw_duplicate_path,
        oracle=oracle_path,
        runner_source=MODULE.COOCCURRENCE_RELEASE_CODE_PATHS["release_runner"],
        output=release_receipt_path,
    )
    release_recomputed = MODULE.RELEASE_FINALIZER.recompute_output_contract(
        sites_path=coocc_sites,
        pairs_path=pair_path,
        duplicates_path=raw_duplicate_path,
        oracle_path=oracle_path,
        expected_stable_sites=len(stable),
    )
    write_json(
        release_receipt_path,
        {
            "schema_name": "intersubmod.cooccurrence_release_receipt",
            "schema_version": MODULE.COOCCURRENCE_RELEASE_SCHEMA_VERSION,
            "started_at_utc": "2026-07-15T00:02:10+00:00",
            "finished_at_utc": "2026-07-15T00:02:20+00:00",
            "task_type": "B_comprehensive_validation",
            "scope": "all_manifest_samples",
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
            "producer_status": "EXECUTION_PASS",
            "release_status": "RELEASE_RECONCILIATION_PASS",
            "pass_semantics": (
                "runner_reconciled_release_integrity_only_not_scientific_claim"
            ),
            "command": MODULE.RELEASE_FINALIZER.canonical_command(release_args),
            "source_authority": source_authority,
            "inputs": release_inputs,
            "code": {
                role: artifact(path)
                for role, path in MODULE.COOCCURRENCE_RELEASE_CODE_PATHS.items()
            },
            "source_modes": {
                role: "0o444" for role in MODULE.COOCCURRENCE_RELEASE_CODE_PATHS
            },
            "recomputed": release_recomputed,
            "checks": {
                key: True for key in MODULE.RELEASE_FINALIZER.RELEASE_CHECK_KEYS
            },
            "pass": True,
        },
    )
    release_receipt_path.chmod(0o444)
    strict_receipt_path = (
        strict_dir / "run_manifest.json"
        if with_candidates
        else strict_dir / "run_receipt.json"
    )
    strict_receipt = json.loads(strict_receipt_path.read_text())
    strict_receipt["inputs"]["cooccurrence_receipt"] = artifact(
        cooccurrence_receipt_path
    )
    strict_receipt["inputs"]["cooccurrence_release_receipt"] = artifact(
        release_receipt_path
    )
    strict_config = MODULE.STRICT_PRODUCER.ConfirmationConfig()
    strict_code = MODULE.STRICT_PRODUCER.capture_source_identity()
    strict_modes = MODULE.STRICT_PRODUCER.capture_source_modes()
    strict_runtime = MODULE.STRICT_PRODUCER.capture_runtime_environment_identity()
    MODULE.STRICT_PRODUCER.require_canonical_runtime(
        strict_runtime,
        output_dir=strict_dir,
        formal_cache_required=False,
    )
    strict_receipt["source_authority"] = source_authority
    strict_receipt["code"] = strict_code
    strict_receipt["source_lock"] = {
        "source_identity_before": strict_code,
        "source_identity_after_compute": strict_code,
        "source_modes_before": strict_modes,
        "source_modes_after_compute": strict_modes,
        "all_sources_read_only_and_unchanged": True,
    }
    strict_receipt["runtime_environment_lock"] = {
        "identity_before": strict_runtime,
        "identity_after_compute": strict_runtime,
        "isolated_runtime_unchanged": True,
    }
    strict_receipt["analysis_replay"] = (
        MODULE.STRICT_PRODUCER.strict_analysis_replay_record(
            strict_site_rows if with_candidates else []
        )
    )
    strict_receipt["command"] = MODULE.STRICT_PRODUCER.canonical_command(
        candidate_table=coocc_sites,
        assignments=assignment_path,
        cooccurrence_receipt=cooccurrence_receipt_path,
        cooccurrence_release_receipt=release_receipt_path,
        output_dir=strict_dir,
        config=strict_config,
    )
    strict_receipt["parameters"] = {
        "permutations_per_seed_per_null": MODULE.STRICT_PRODUCER.DEFAULT_PERMUTATIONS,
        "seeds": MODULE.STRICT_PRODUCER.DEFAULT_SEEDS,
        "formal_minimum_permutations": MODULE.STRICT_PRODUCER.FORMAL_MIN_PERMUTATIONS,
        "formal_minimum_seeds": MODULE.STRICT_PRODUCER.FORMAL_MIN_SEEDS,
        "formal_parameter_gate": True,
    }
    write_json(strict_receipt_path, strict_receipt)
    strict_receipt_path.chmod(0o444)

    tumor_receipt_path = tumor_ref_dir / "run_manifest.json"
    tumor_receipt = json.loads(tumor_receipt_path.read_text())
    tumor_receipt["inputs"]["primary_artifact_audit_pre"] = artifact(audit_pre)
    write_json(tumor_receipt_path, tumor_receipt)

    post_audit_payload = deepcopy(audit_payload)
    post_audit_payload["started_at_utc"] = "2026-07-15T00:03:00+00:00"
    post_audit_payload["finished_at_utc"] = "2026-07-15T00:04:00+00:00"
    post_audit_payload["created_at_utc"] = "2026-07-15T00:04:00+00:00"
    consumer_receipt_paths = [
        cooccurrence_preflight,
        cooccurrence_receipt_path,
        release_receipt_path,
        tumor_receipt_path,
        strict_receipt_path,
    ]
    if normal_dir is not None:
        consumer_receipt_paths.extend(
            [paired_receipt_path, normal_dir / "run_receipt.json"]
        )
    post_audit_payload["inputs"]["consumer_receipts"] = [
        artifact(path) for path in consumer_receipt_paths
    ]
    post_audit_payload["command"] = MODULE.PRIMARY_AUDITOR.canonical_command(
        site_results=screen_sites,
        assignments=assignment_path,
        consumer_receipts=consumer_receipt_paths,
        output=audit_post,
        **MODULE.PRIMARY_AUDITOR.CANONICAL_PARAMETERS,
    )
    write_json(audit_post, post_audit_payload)
    audit_post.chmod(0o444)

    return MODULE.InputBundle(
        manifest=manifest_path,
        screen_dir=screen_dir,
        cooccurrence_dir=cooccurrence_dir,
        strict_dir=strict_dir,
        tumor_ref_dir=tumor_ref_dir,
        primary_artifact_audit_pre=audit_pre,
        primary_artifact_audit_post=audit_post,
        cooccurrence_preflight=cooccurrence_preflight,
        matched_normal_dir=normal_dir,
    )


def add_matched_normal_na_fixture(
    bundle: MODULE.InputBundle, root: Path
) -> MODULE.InputBundle:
    normal_dir = root / "matched_normal_na"
    normal_dir.mkdir(parents=True)
    normal_audit = root / "normal_audit.json"
    binary = root / "inter_sub_mod"
    reference = root / "reference.fa"
    runner = Path(MODULE.MATCHED_NORMAL_RUNNER.__file__).resolve()
    write_json(normal_audit, {"pass": True})
    for path, content in (
        (binary, "binary\n"),
        (reference, ">chr1\nA\n"),
    ):
        path.write_text(content, encoding="utf-8")
    receipt_path = normal_dir / "not_applicable_receipt.json"
    runner_provenance = producer_provenance(
        MODULE.MATCHED_NORMAL_RUNNER, "matched_normal_runner"
    )
    write_json(
        receipt_path,
        {
            "schema_name": "intersubmod.matched_normal_candidate_control.not_applicable",
            "schema_version": "1.0.0",
            "started_at_utc": "2026-07-15T00:01:00+00:00",
            "finished_at_utc": "2026-07-15T00:02:00+00:00",
            "task_type": "B_comprehensive_validation",
            "status": "NOT_APPLICABLE",
            "execution_status": "NOT_APPLICABLE",
            "reason": "ZERO_SELECTED_CANDIDATES",
            "n_selected_candidates": 0,
            "selection_column": MODULE.STRICT_SELECTION_COLUMN,
            "selection_value": "true",
            "selection_contract": {
                "selection_column": MODULE.STRICT_SELECTION_COLUMN,
                "selection_value": "true",
                "n_selected_candidates": 0,
            },
            "counts": {"n_selected_candidates": 0},
            "command": MODULE.MATCHED_NORMAL_RUNNER.canonical_task_b_command(),
            **runner_provenance,
            "inputs": {
                "candidate_table": artifact(bundle.cooccurrence_sites),
                "all_ssnv_manifest": artifact(bundle.manifest),
                "normal_audit": artifact(normal_audit),
                "binary": artifact(binary),
                "reference": artifact(reference),
                "runner_script": artifact(runner),
            },
            "outputs": {
                "output_root": str(normal_dir.resolve()),
                "not_applicable_receipt": str(receipt_path.resolve()),
                "sample_outputs": [],
            },
            "sample_outputs_created": False,
            "cpp_executed": False,
            "normal_control_executed": False,
            "not_evaluable_is_negative": False,
            "interpretation": "NOT_APPLICABLE is not a negative matched-normal result.",
            "pass_semantics": (
                "receipt_integrity_only_not_normal_control_execution_or_negative_evidence"
            ),
            "pass": True,
        },
    )
    receipt_path.chmod(0o444)
    post_audit = json.loads(bundle.primary_artifact_audit_post.read_text())
    post_audit["inputs"]["consumer_receipts"].append(artifact(receipt_path))
    save_primary_post_audit(bundle, post_audit)
    return replace(bundle, matched_normal_dir=normal_dir)


def add_zero_row_matched_normal_analysis_fixture(
    bundle: MODULE.InputBundle, root: Path
) -> MODULE.InputBundle:
    normal_dir = root / "matched_normal_analysis_zero"
    normal_dir.mkdir(parents=True)
    fields = [
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "normal_called_reads",
        "normal_alt_reads",
        "normal_ref_reads",
        "normal_unknown_reads",
        "normal_focal_callability_gate",
        "normal_control_evaluable",
        "normal_control_not_evaluable_reason",
        "normal_ref_methyl_stable_multigroup",
        "normal_ref_methyl_nonreplication_gate",
        "normal_control_status",
        "normal_evaluable",
        "normal_stable_multigroup",
        "normal_genetic_alt_support_present",
        "tumor_ref_evaluable",
        "tumor_ref_stable_multigroup",
        "primary_group_assignment_coverage",
        "primary_identity_collision_count",
        "primary_identity_missing_count",
    ]
    site_path = normal_dir / "matched_normal_candidate_controls.tsv.gz"
    write_tsv(site_path, [], fields)
    paired_bundle = add_matched_normal_na_fixture(bundle, root / "paired_provenance")
    assert paired_bundle.matched_normal_dir is not None
    paired_receipt = paired_bundle.matched_normal_dir / "not_applicable_receipt.json"
    counts = {"n_candidates": 0, "all_primary_group_assignments_exact": True}
    summary_path = normal_dir / "matched_normal_candidate_controls_summary.json"
    analyzer_provenance = producer_provenance(
        MODULE.MATCHED_NORMAL_ANALYZER, "matched_normal_analyzer"
    )
    write_json(
        summary_path,
        {
            "schema_name": "intersubmod.matched_normal_candidate_control_analysis",
            "schema_version": "2.0.0",
            **analyzer_provenance,
            "pooled": counts,
            "not_evaluable_is_negative_result": False,
            "pass_semantics": (
                "execution_and_identity_integrity_only_not_background_negativity"
            ),
            "pass": True,
        },
    )
    site_path.chmod(0o444)
    summary_path.chmod(0o444)
    receipt_path = normal_dir / "run_receipt.json"
    write_json(
        receipt_path,
        {
            "schema_name": "intersubmod.matched_normal_candidate_control_analysis_run",
            "schema_version": "2.0.0",
            "started_at_utc": "2026-07-15T00:01:00+00:00",
            "finished_at_utc": "2026-07-15T00:02:00+00:00",
            "command": MODULE.MATCHED_NORMAL_ANALYZER.canonical_task_b_command(),
            **analyzer_provenance,
            "inputs": {
                "primary_assignments": artifact(bundle.screen_assignments),
                "paired_run_receipt": artifact(paired_receipt),
            },
            "outputs": {
                "site_table": artifact(site_path),
                "summary": artifact(summary_path),
            },
            "counts": counts,
            "not_evaluable_is_negative_result": False,
            "pass_semantics": (
                "execution_and_identity_integrity_only_not_background_negativity"
            ),
            "pass": True,
        },
    )
    receipt_path.chmod(0o444)
    post_audit = json.loads(bundle.primary_artifact_audit_post.read_text())
    post_audit["inputs"]["consumer_receipts"].append(
        artifact(receipt_path)
    )
    save_primary_post_audit(bundle, post_audit)
    return replace(bundle, matched_normal_dir=normal_dir)


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_native_cn_fixture(
    bundle: MODULE.InputBundle, root: Path, *, zero_rows: bool = False
) -> MODULE.InputBundle:
    cn_dir = root / "native_cn"
    cn_dir.mkdir(parents=True)
    authority_dir = root / "cn_authority"
    authority_dir.mkdir(parents=True)
    config = authority_dir / "config.json"
    provenance = authority_dir / "provenance.json"
    analysis_summary = authority_dir / "analysis_summary.json"
    source = authority_dir / "source.tsv"
    write_json(config, {"authority": "config"})
    write_json(provenance, {"authority": "provenance"})
    write_json(analysis_summary, {"authority": "summary"})
    source.write_text("source\n", encoding="utf-8")
    source_sha = MODULE.sha256(source)
    selected = [] if zero_rows else [
        row
        for row in read_gzip_tsv(bundle.cooccurrence_sites)
        if row[MODULE.FORMAL_SELECTION_COLUMN].lower() == "true"
    ]
    rows = []
    for selected_row in selected:
        sample = selected_row["sample"]
        row = {field_name: "" for field_name in MODULE.CN_CCF_OUTPUT_COLUMNS}
        row.update(
            {
                "sample": sample,
                "chrom": selected_row["chrom"],
                "pos": selected_row["pos"],
                "ref": selected_row["ref"],
                "alt": selected_row["alt"],
                "mutation_id": (
                    f"{selected_row['chrom']}:{selected_row['pos']}:"
                    f"{selected_row['ref']}>{selected_row['alt']}"
                ),
                "callset_status": "INPUT_CANDIDATE_KEY_VALIDATED_NOT_RECHECKED",
                "cn_status": (
                    "SHARED_CN_SENSITIVITY"
                    if sample == "HCC1395_DORADO"
                    else "AVAILABLE_EXACT_SEGMENT"
                ),
                "cn_authority_sample": "HCC1395",
                "cn_transfer_policy": (
                    "SAME_CELL_LINE_SHARED_CN_SENSITIVITY"
                    if sample == "HCC1395_DORADO"
                    else "SAMPLE_SPECIFIC_MEASURED_CN"
                ),
                "independent_cn": "false" if sample == "HCC1395_DORADO" else "true",
                "cn_segment_id": "segment_1",
                "cn_segment_start": "1",
                "cn_segment_end": "1000",
                "savana_total_cn_raw": "2.2",
                "savana_minor_cn_raw": "1.0",
                "savana_total_cn_discrete": "2",
                "savana_major_cn_discrete": "1",
                "savana_minor_cn_discrete": "1",
                "purity_model_value": "0.85",
                "pyclone_status": "MATCHED_PRIMARY",
                "pyclone_primary_bundle_id": f"{sample}_main",
                "pyclone_primary_model_role": "primary",
                "pyclone_fit_local_cluster_id": "1",
                "pyclone_vi_cellular_prevalence": "0.42",
                "pyclone_vi_cellular_prevalence_std": "0.03",
                "pyclone_vi_assignment_probability": "0.99",
                "pyclone_sensitivity_status": "MATCHED_SENSITIVITY_ONLY",
                "pyclone_sensitivity_bundle_id": f"{sample}_sensitivity",
                "pyclone_sensitivity_fit_local_cluster_id": "1",
                "pyclone_sensitivity_cellular_prevalence": "0.40",
                "pyclone_sensitivity_cellular_prevalence_std": "0.04",
                "pyclone_sensitivity_assignment_probability": "0.98",
                "authority_config_sha256": MODULE.sha256(config),
                "input_provenance_sha256": MODULE.sha256(provenance),
                "analysis_summary_sha256": MODULE.sha256(analysis_summary),
                "cn_source_sha256": source_sha,
                "pyclone_primary_metadata_sha256": source_sha,
                "pyclone_primary_results_sha256": source_sha,
                "pyclone_primary_status_sha256": source_sha,
                "pyclone_sensitivity_metadata_sha256": source_sha,
                "pyclone_sensitivity_results_sha256": source_sha,
                "pyclone_sensitivity_status_sha256": source_sha,
                "claim_ceiling": MODULE.CN_CCF_CLAIM_CEILING,
            }
        )
        rows.append(row)
    output_path = cn_dir / "candidate_cn_ccf_annotations.tsv.gz"
    write_tsv(output_path, rows, list(MODULE.CN_CCF_OUTPUT_COLUMNS))
    receipt_path = cn_dir / "receipt.json"
    zero = not rows
    status_counts = {
        field_name: dict(sorted(Counter(row[field_name] for row in rows).items()))
        for field_name in ("cn_status", "pyclone_status", "pyclone_sensitivity_status")
    }
    input_rows = read_gzip_tsv(bundle.cooccurrence_sites)
    cn_provenance = producer_provenance(
        MODULE.CN_CCF_ANNOTATOR, "cn_ccf_annotator"
    )
    write_json(
        receipt_path,
        {
            "schema_name": "intersubmod.candidate_cn_ccf_annotation_receipt",
            "schema_version": "1.0.0",
            "pass": True,
            "status": "NOT_APPLICABLE" if zero else "PASS",
            "execution_status": "NOT_APPLICABLE" if zero else "EXECUTION_PASS",
            "reason": "ZERO_SELECTED_CANDIDATES" if zero else None,
            "n_selected_candidates": len(rows),
            "command": MODULE.CN_CCF_ANNOTATOR.canonical_task_b_command(),
            **cn_provenance,
            "input": {
                **artifact(bundle.cooccurrence_sites),
                "rows_read_total": len(input_rows),
                "rows_in": len(rows),
                "selection_column": MODULE.FORMAL_SELECTION_COLUMN,
                "selection_value": "true",
            },
            "output": {
                **artifact(output_path),
                "rows_out": len(rows),
                "columns": list(MODULE.CN_CCF_OUTPUT_COLUMNS),
            },
            "conservation": {
                "rows_in": len(rows),
                "rows_out": len(rows),
                "rows_in_equals_rows_out": True,
            },
            "authority": {
                "config": artifact(config),
                "input_provenance": artifact(provenance),
                "analysis_summary": artifact(analysis_summary),
                "all_source_hashes": {str(source.resolve()): source_sha},
            },
            "status_counts": status_counts,
            "status_enums": {
                "cn_status": sorted(MODULE.CN_STATUS_ENUM),
                "pyclone_status": sorted(MODULE.PYCLONE_STATUS_ENUM),
            },
            "claim_ceiling": MODULE.CN_CCF_CLAIM_CEILING,
            "pass_semantics": MODULE.CN_CCF_PASS_SEMANTICS,
            "scientific_interpretation": {
                "is_negative_result": False,
                "c1_formed": False,
                "statement": (
                    "ZERO_SELECTED_CANDIDATES is not a biological negative; "
                    "C1 was not evaluated or formed."
                    if zero
                    else "Conditional CN/PyClone annotation does not itself form C1."
                ),
            },
        },
    )
    output_path.chmod(0o444)
    receipt_path.chmod(0o444)
    return replace(bundle, cn_ccf_annotations=cn_dir)


def metric(dataset: dict, name: str) -> dict:
    return dataset["funnel_metrics"]["pooled"][name]


def refresh_cooccurrence_chain(
    bundle: MODULE.InputBundle, *, recompute_release: bool = False
) -> None:
    receipt = json.loads(bundle.cooccurrence_receipt.read_text())
    for key, path in (
        ("site_tsv", bundle.cooccurrence_sites),
        ("pair_tsv", bundle.cooccurrence_pairs),
        ("summary_json", bundle.cooccurrence_summary),
        ("raw_identity_duplicate_audit_tsv", bundle.cooccurrence_raw_identity_duplicates),
        ("case_json", bundle.cooccurrence_oracle_cases),
    ):
        receipt["outputs"][key] = artifact(path)
    receipt["inputs"]["preflight_receipt"] = artifact(
        bundle.cooccurrence_preflight
    )
    receipt["input_lock"] = {
        "identity_before": receipt["inputs"],
        "identity_after_compute": receipt["inputs"],
        "identity_after_output": receipt["inputs"],
        "all_primary_inputs_unchanged": True,
    }
    write_json(bundle.cooccurrence_receipt, receipt)
    release = json.loads(bundle.cooccurrence_release_receipt.read_text())
    release["inputs"].update(
        {
            "preflight": artifact(bundle.cooccurrence_preflight),
            "producer_receipt": artifact(bundle.cooccurrence_receipt),
            "summary": artifact(bundle.cooccurrence_summary),
            "sites": artifact(bundle.cooccurrence_sites),
            "pairs": artifact(bundle.cooccurrence_pairs),
            "duplicates": artifact(bundle.cooccurrence_raw_identity_duplicates),
            "oracle": artifact(bundle.cooccurrence_oracle_cases),
        }
    )
    if recompute_release:
        release["recomputed"] = MODULE.RELEASE_FINALIZER.recompute_output_contract(
            sites_path=bundle.cooccurrence_sites,
            pairs_path=bundle.cooccurrence_pairs,
            duplicates_path=bundle.cooccurrence_raw_identity_duplicates,
            oracle_path=bundle.cooccurrence_oracle_cases,
            expected_stable_sites=len(read_gzip_tsv(bundle.cooccurrence_sites)),
        )
    bundle.cooccurrence_release_receipt.chmod(0o644)
    write_json(bundle.cooccurrence_release_receipt, release)
    bundle.cooccurrence_release_receipt.chmod(0o444)
    strict_receipt = json.loads(bundle.strict_receipt.read_text())
    strict_receipt["inputs"]["candidate_table"] = artifact(
        bundle.cooccurrence_sites
    )
    strict_receipt["inputs"]["cooccurrence_receipt"] = artifact(
        bundle.cooccurrence_receipt
    )
    strict_receipt["inputs"]["cooccurrence_release_receipt"] = artifact(
        bundle.cooccurrence_release_receipt
    )
    write_json(bundle.strict_receipt, strict_receipt)
    bundle.strict_receipt.chmod(0o444)
    refresh_post_audit_consumer_identities(bundle)


def save_primary_post_audit(
    bundle: MODULE.InputBundle, post_audit: dict
) -> None:
    consumer_paths = [
        Path(str(reference["path"]))
        for reference in post_audit["inputs"]["consumer_receipts"]
    ]
    post_audit["command"] = MODULE.PRIMARY_AUDITOR.canonical_command(
        site_results=bundle.screen_sites,
        assignments=bundle.screen_assignments,
        consumer_receipts=consumer_paths,
        output=bundle.primary_artifact_audit_post,
        **MODULE.PRIMARY_AUDITOR.CANONICAL_PARAMETERS,
    )
    write_json(bundle.primary_artifact_audit_post, post_audit)
    bundle.primary_artifact_audit_post.chmod(0o444)


def refresh_post_audit_consumer_identities(bundle: MODULE.InputBundle) -> None:
    post_audit = json.loads(bundle.primary_artifact_audit_post.read_text())
    post_audit["inputs"]["consumer_receipts"] = [
        artifact(Path(str(reference["path"])))
        for reference in post_audit["inputs"]["consumer_receipts"]
    ]
    save_primary_post_audit(bundle, post_audit)


def install_single_raw_identity_duplicate(
    bundle: MODULE.InputBundle,
) -> tuple[str, str, int, int, int, str]:
    rows = read_gzip_tsv(bundle.cooccurrence_sites)
    fields = list(rows[0])
    target = next(row for row in rows if row["sample"] == "H2009")
    projection = (
        "duplicate-read",
        target["chrom"],
        int(target["pos"]) - 10,
        int(target["pos"]) + 60,
        60,
        "+",
    )
    expected = int(target["raw_identity_expected_projections"])
    target.update(
        {
            "raw_identity_matched_records": str(expected + 1),
            "raw_identity_duplicate_projections_collapsed": "1",
            "raw_identity_duplicate_extra_records_collapsed": "1",
            "raw_identity_exact_duplicate_projections_collapsed": "0",
            "raw_identity_rg_only_duplicate_projections_collapsed": "1",
            "raw_identity_duplicate_projection_sha256": MODULE.projection_digest(
                [projection]
            ),
        }
    )
    write_tsv(bundle.cooccurrence_sites, rows, fields)
    write_tsv(
        bundle.cooccurrence_raw_identity_duplicates,
        [
            {
                "sample": target["sample"],
                "chrom": target["chrom"],
                "pos": target["pos"],
                "projection_read_name": projection[0],
                "projection_chrom": projection[1],
                "projection_start": projection[2],
                "projection_end": projection[3],
                "projection_mapq": projection[4],
                "projection_strand": projection[5],
                "record_count": 2,
                "classification": "RG_ONLY_DUPLICATE",
                "differing_auxiliary_tags": json.dumps(["RG"]),
            }
        ],
        fields=list(MODULE.RAW_IDENTITY_DUPLICATE_COLUMNS),
    )
    summary_audit, preflight_audit = raw_identity_fixture_audits(rows)
    summary = json.loads(bundle.cooccurrence_summary.read_text())
    summary["raw_bam_identity_recovery_audit"] = summary_audit
    write_json(bundle.cooccurrence_summary, summary)
    preflight = json.loads(bundle.cooccurrence_preflight.read_text())
    preflight["observed"]["raw_bam_identity_recovery"] = preflight_audit
    write_json(bundle.cooccurrence_preflight, preflight)
    receipt = json.loads(bundle.cooccurrence_receipt.read_text())
    receipt["counts"].update(
        {
            "raw_identity_expected_projection_occurrences": summary_audit[
                "n_expected_projection_occurrences"
            ],
            "raw_identity_duplicate_projection_occurrences_collapsed": 1,
            "raw_identity_duplicate_extra_record_occurrences_collapsed": 1,
            "raw_identity_sparse_duplicate_rows": 1,
        }
    )
    write_json(bundle.cooccurrence_receipt, receipt)
    refresh_cooccurrence_chain(bundle, recompute_release=True)
    return projection


def refresh_tumor_ref_chain(bundle: MODULE.InputBundle) -> None:
    receipt = json.loads(bundle.tumor_ref_receipt.read_text())
    receipt["outputs"]["site_results"] = artifact(bundle.tumor_ref_sites)
    write_json(bundle.tumor_ref_receipt, receipt)
    refresh_post_audit_consumer_identities(bundle)


def refresh_matched_normal_chain(bundle: MODULE.InputBundle) -> None:
    assert bundle.matched_normal_receipt is not None
    assert bundle.matched_normal_sites is not None
    receipt = json.loads(bundle.matched_normal_receipt.read_text())
    receipt["outputs"]["site_table"] = artifact(bundle.matched_normal_sites)
    write_json(bundle.matched_normal_receipt, receipt)
    refresh_post_audit_consumer_identities(bundle)


def test_exact_joins_denominators_and_technical_replicate(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    output_dir = tmp_path / "final"
    outputs = MODULE.build_outputs(
        bundle, output_dir, expected_screen_sites=8, command=["synthetic-test"]
    )
    assert set(outputs) == {
        "final_report_dataset.json",
        "per_sample_metrics.tsv",
        "candidate_catalog.tsv",
        "candidate_witness_pairs.tsv",
        "claim_ladder.tsv",
        "run_receipt.json",
    }
    dataset = json.loads(outputs["final_report_dataset.json"].read_text())
    receipt = json.loads(outputs["run_receipt.json"].read_text())
    assert dataset["scope"]["dataset_count"] == 7
    assert dataset["scope"]["biological_sample_count"] == 6
    assert dataset["scope"]["technical_replicate"]["counts_as_independent_biological_n"] is False
    assert dataset["schema_version"] == "2.0.0"
    assert dataset["pass_semantics"] == MODULE.PASS_SEMANTICS
    assert dataset["task_type"] == "NON_RELEASE_TEST_FIXTURE_VALIDATION"
    assert dataset["formal_task_b_release_eligible"] is False
    assert receipt["schema_version"] == "2.0.0"
    assert receipt["pass_semantics"] == MODULE.PASS_SEMANTICS
    assert receipt["task_type"] == "NON_RELEASE_TEST_FIXTURE_VALIDATION"
    assert receipt["formal_task_b_release_eligible"] is False
    assert receipt["validations"]["task_type_b_full_scope"] is False
    producer = receipt["code"]["final_report_dataset_builder"]
    assert Path(producer["path"]).resolve() == SCRIPT.resolve()
    assert producer["sha256"] == MODULE.sha256(SCRIPT)
    assert metric(dataset, "M1") == {
        "numerator": 3,
        "denominator": 8,
        "ratio": 0.375,
        "not_evaluable": 0,
        "not_run": 0,
        "population": 8,
        "denominator_definition": MODULE.CLAIM_RULES["M1"]["denominator"],
    }
    assert dataset["m1_operational_screen"] == {
        "status_semantics": "FLAGGED_VS_NOT_FLAGGED_OPERATIONAL_SCREEN_ONLY",
        "denominator_definition": MODULE.CLAIM_RULES["M1"]["denominator"],
        "n_all_dataset_sites": 8,
        "n_screen_evaluable": 6,
        "n_screen_not_evaluable": 2,
        "n_flagged_stable_multigroup": 3,
        "n_not_flagged_all": 5,
        "flag_yield": 0.375,
        "flag_yield_among_screen_evaluable": 0.5,
        "n_evaluable_not_flagged": 3,
        "global_null_validity_exported_for_nonstable_sites": False,
        "nonflagged_scientific_interpretation": (
            "NOT_IDENTIFIABLE_AS_TRUE_NEGATIVE_VS_NULL_INVALID_FROM_SITE_TSV"
        ),
        "biological_prevalence_estimate": None,
    }
    assert dataset["background_control_replication_gate"] == {
        "contract": MODULE.BACKGROUND_CONTROL_REPLICATION_GATE_CONTRACT,
        "applies_to": ["tumor_REF", "matched_normal_REF"],
        "required_conditions": [
            "coarse_ng>=2",
            "modal_fraction>=0.7_via_unstable_false",
        ],
        "membership_ari_minimum_required": False,
        "relation_to_primary_m1_replication_flags": (
            MODULE.BACKGROUND_CONTROL_RELATION_TO_PRIMARY_M1
        ),
        "b1_pass_direction": "requires_no_lenient_background_replication",
        "false_positive_direction": (
            "cannot_increase_B1_passes_vs_ARI_qualified_predicate_on_same_background_payload"
        ),
        "false_negative_direction": (
            "may_conservatively_reduce_B1_passes_when_K_is_stable_but_membership_is_not"
        ),
        "scientific_interpretation": (
            "background nonreplication guardrail, not an exact primary-M1 replay"
        ),
    }
    assert metric(dataset, "M2")["numerator"] == 2
    assert metric(dataset, "M2")["denominator"] == 3
    assert dataset["m2_evaluability_contract"] == {
        "gate_contract": MODULE.M2_GATE.GATE_CONTRACT,
        "denominator_definition": MODULE.CLAIM_RULES["M2"]["denominator"],
        "minimum_supported_methyl_groups": 2,
        "maximum_supported_methyl_groups": 10,
        "categorical_planning_level_ceilings": {
            "hp_exact": 7,
            "hp_family": 5,
            "strand": 2,
        },
        "assignment_observed_levels_role": (
            "constant-axis proof only; observed levels do not replace the "
            "source-locked planning level ceilings"
        ),
        "n_m1_pass": 3,
        "n_m2_evaluable": 3,
        "n_m2_not_evaluable": 0,
        "not_evaluable_reason_counts": {},
        "n_group_count_exceeds_planning_model_maximum": 0,
        "group_count_exceeds_planning_model_examples": [],
        "group_count_exceeds_examples_complete": True,
        "group_count_exceeds_claim_behavior": {
            "M1": "PASS retained",
            "M2": "NOT_EVALUABLE excluded from PASS/FAIL denominator",
            "G1": "NOT_RUN",
            "G2": "NOT_RUN",
            "B1": "NOT_RUN",
        },
        "independent_logic_audit": {"status": "NOT_INCLUDED", "pass": None},
    }
    assert metric(dataset, "G2")["numerator"] == 2
    assert metric(dataset, "R1")["numerator"] == 1
    assert metric(dataset, "R1")["denominator"] == 1
    assert metric(dataset, "R1")["not_evaluable"] == 1

    candidates = read_tsv(outputs["candidate_catalog.tsv"])
    assert len(candidates) == 2
    by_sample = {row["sample"]: row for row in candidates}
    assert by_sample["HCC1395"]["r1_status"] == "PASS"
    assert by_sample["HCC1395"]["b1_status"] == "PASS"
    assert by_sample["HCC1395"]["c1_status"] == "NOT_EVALUABLE"
    assert by_sample["HCC1395"]["cn_total"] == ""
    assert by_sample["HCC1395"]["c1_reason"] == (
        "NO_NATIVE_AUTHORITY_LOCKED_CN_PYCLONE_ANNOTATION"
    )
    assert by_sample["HCC1395_DORADO"]["r1_status"] == "NOT_EVALUABLE"
    assert by_sample["HCC1395_DORADO"]["b1_status"] == "NOT_EVALUABLE"
    assert int(by_sample["HCC1395"]["n_same_pair_four_state_witnesses"]) == 1
    assert int(
        by_sample["HCC1395"]["n_four_state_compatible_formal_pair_opportunities"]
    ) == 2
    assert by_sample["HCC1395"]["b1_uses_posthoc_compatible_pair_search"] == "false"
    assert json.loads(by_sample["HCC1395"]["m2_low_power_axes"]) == []
    assert by_sample["HCC1395"]["joint_signature_global_fdr_family_status"] == (
        "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
    )
    assert float(by_sample["HCC1395"]["joint_signature_q_global_by"]) == 0.015
    assert by_sample["HCC1395"]["joint_signature_global_by_discovery"] == "true"

    witness_rows = read_tsv(outputs["candidate_witness_pairs.tsv"])
    assert witness_rows
    assert {float(row["four_state_familywise_confidence"]) for row in witness_rows} == {
        0.95
    }


    assert {int(row["four_state_relation_family_size"]) for row in witness_rows} == {3}
    assert {row["four_state_multiplicity_method"] for row in witness_rows} == {
        "bonferroni_three_relation_models"
    }

    assert dataset["technical_replication"] == {
        "endpoint": "HCC1395_vs_HCC1395_DORADO_exact_pair_technical_concordance",
        "status": "ANY_CONCORDANT_EXACT_PAIR_OBSERVED",
        "numerator": 1,
        "denominator": 2,
        "ratio": 0.5,
        "not_evaluable_one_platform_only": 0,
        "denominator_definition": (
            "exact shared focal+partner CHROM/POS/REF/ALT opportunities in both HCC1395 pipelines"
        ),
        "biological_n": 1,
        "independent_biological_replication_n": 0,
        "replication_claim_status": "NOT_EVALUABLE_BIOLOGICAL_N1",
        "inferential_confidence_interval": None,
        "pair_independence_assumption_met": False,
        "required_for_b1": False,
        "same_exact_pair_evidence_required": True,
    }

    metric_rows = read_tsv(outputs["per_sample_metrics.tsv"])
    replicate = next(
        row
        for row in metric_rows
        if row["sample"] == "HCC1395_DORADO" and row["claim_id"] == "M1"
    )
    assert replicate["dataset_role"] == "technical_replicate"
    assert replicate["biological_n_contribution"] == "0"
    primary = next(
        row
        for row in metric_rows
        if row["sample"] == "HCC1395" and row["claim_id"] == "M1"
    )
    assert primary["biological_n_contribution"] == "1"

    claims = {row["claim_id"]: row for row in read_tsv(outputs["claim_ladder.tsv"])}
    assert tuple(claims) == MODULE.CLAIM_IDS
    assert claims["M1"]["dataset_numerator"] == "3"
    assert claims["M1"]["dataset_denominator"] == "8"
    assert claims["M1"]["dataset_not_evaluable"] == "0"
    assert claims["G2"]["claim_name"] == "multi-marker molecular-haplotype base candidate"
    assert claims["B1"]["dataset_numerator"] == "1"
    assert claims["B1"]["dataset_denominator"] == "1"
    assert claims["C1"]["status"] == "NOT_EVALUABLE"
    assert all(row["automatic_upgrade_prohibited"] == "true" for row in claims.values())
    assert len(dataset["funnel_metrics"]["sample_by_truth"]) == 21

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)


def test_independent_primary_recount_rehashes_all_live_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assignments_path = tmp_path / "independent_primary_assignments.jsonl.gz"
    assignments = []
    expected_keys = set()
    digest_lines = []
    for index in range(2):
        key = ("HCC1395", "chr1", 100 + index, "A", "C")
        region_dir = tmp_path / f"region_{index}"
        contracts = {}
        for role, relative_path in MODULE.INDEPENDENT_PRIMARY_ARTIFACT_PATHS.items():
            path = region_dir / relative_path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"{index}:{role}\n", encoding="utf-8")
            contracts[role] = artifact(path)
            digest_lines.append(
                "|".join(
                    (
                        *(str(value) for value in key),
                        role,
                        str(path.resolve()),
                        str(path.stat().st_size),
                        MODULE.sha256(path),
                    )
                )
            )
        assignments.append(
            {
                "schema_name": MODULE.ASSIGNMENT_SCHEMA,
                "schema_version": MODULE.ASSIGNMENT_SCHEMA_VERSION,
                "screen_contract": (
                    "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"
                ),
                "artifact_identity_contract": "sha256_size_path_v1",
                "strict_confirm_candidate": True,
                "sample": key[0],
                "chrom": key[1],
                "pos": key[2],
                "posthoc": {"ref": key[3], "alt": key[4]},
                "region_dir": str(region_dir.resolve()),
                "primary_artifacts": contracts,
            }
        )
        expected_keys.add(key)
    with gzip.open(assignments_path, "wt", encoding="utf-8") as handle:
        for assignment in assignments:
            handle.write(json.dumps(assignment, sort_keys=True) + "\n")
    expected_digest = hashlib.sha256(
        "\n".join(sorted(digest_lines)).encode("utf-8")
    ).hexdigest()
    monkeypatch.setattr(MODULE, "INDEPENDENT_PRIMARY_RECOUNT_WORKERS", 2)
    record = MODULE.independently_recount_primary_artifacts(
        assignments_path,
        expected_keys,
        expected_digest,
    )
    assert record["implementation_independence"] is True
    assert record["primary_auditor_functions_called"] is False
    assert record["stable_sites"] == len(expected_keys)
    assert record["primary_artifacts_verified"] == len(expected_keys) * 3
    assert record["artifact_set_sha256"] == expected_digest


def test_eleven_group_m1_site_is_m2_not_evaluable_and_downstream_not_run(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(
        tmp_path,
        with_candidates=True,
        eleven_group_site=True,
    )
    outputs = MODULE.build_outputs(
        bundle,
        tmp_path / "final-eleven-groups",
        expected_screen_sites=8,
        command=["synthetic-eleven-group-test"],
    )
    dataset = json.loads(outputs["final_report_dataset.json"].read_text())
    assert metric(dataset, "M1")["numerator"] == 3
    assert metric(dataset, "M2") == {
        "numerator": 2,
        "denominator": 2,
        "ratio": 1.0,
        "not_evaluable": 1,
        "not_run": 5,
        "population": 8,
        "denominator_definition": MODULE.CLAIM_RULES["M2"]["denominator"],
    }
    h2009 = dataset["funnel_metrics"]["per_sample"]["H2009"]
    assert h2009["M1"]["numerator"] == 1
    assert h2009["M2"]["not_evaluable"] == 1
    assert h2009["M2"]["denominator"] == 0
    for claim_id in ("G1", "G2", "B1"):
        assert h2009[claim_id]["not_run"] == 1
        assert h2009[claim_id]["denominator"] == 0

    contract = dataset["m2_evaluability_contract"]
    reason = "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM"
    assert contract["n_group_count_exceeds_planning_model_maximum"] == 1
    assert contract["not_evaluable_reason_counts"] == {reason: 1}
    assert contract["group_count_exceeds_examples_complete"] is True
    assert contract["group_count_exceeds_planning_model_examples"] == [
        {
            "dataset": "H2009",
            "chrom": "chr5",
            "pos": 500,
            "ref": "A",
            "alt": "G",
            "observed_methyl_groups": 11,
            "maximum_supported_methyl_groups": 10,
            "m1_status": "PASS",
            "m2_status": "NOT_EVALUABLE",
            "g1_status": "NOT_RUN",
            "g2_status": "NOT_RUN",
            "b1_status": "NOT_RUN",
            "reason": f"{reason}:observed=11:maximum=10",
        }
    ]


@pytest.mark.parametrize(
    ("control", "expected_reason"),
    [
        (
            "tumor_ref",
            "TUMOR_REF_LENIENT_COARSE_MODAL_PARTITION_REPRODUCED",
        ),
        (
            "matched_normal",
            "NORMAL_REF_LENIENT_COARSE_MODAL_PARTITION_REPRODUCED",
        ),
    ],
)
def test_lenient_background_replication_monotonically_vetoes_b1(
    tmp_path: Path,
    control: str,
    expected_reason: str,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    if control == "tumor_ref":
        with gzip.open(
            bundle.tumor_ref_sites, "rt", encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = list(reader)
        target = next(row for row in rows if row["sample"] == "HCC1395")
        target["ref_stable_null_multigroup"] = "True"
        write_tsv(bundle.tumor_ref_sites, rows, fields)
        refresh_tumor_ref_chain(bundle)
    else:
        assert bundle.matched_normal_sites is not None
        with gzip.open(
            bundle.matched_normal_sites, "rt", encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = list(reader)
        target = next(row for row in rows if row["sample"] == "HCC1395")
        target["normal_ref_methyl_stable_multigroup"] = "True"
        target["normal_stable_multigroup"] = "True"
        target["normal_ref_methyl_nonreplication_gate"] = "False"
        target["normal_control_status"] = "REF_METHYL_PARTITION_REPRODUCED"
        write_tsv(bundle.matched_normal_sites, rows, fields)
        refresh_matched_normal_chain(bundle)

    outputs = MODULE.build_outputs(
        bundle,
        tmp_path / f"{control}_replication_veto",
        expected_screen_sites=8,
    )
    candidates = {
        row["sample"]: row for row in read_tsv(outputs["candidate_catalog.tsv"])
    }
    assert candidates["HCC1395"]["b1_status"] == "FAIL"
    assert candidates["HCC1395"]["b1_reason"] == expected_reason


def test_tumor_ref_source_identity_attestation_is_a_fail_closed_release_gate(
    tmp_path: Path,
) -> None:
    analyzer = tmp_path / "analyzer.py"
    library = tmp_path / "focal_alt_cluster_lib.py"
    analyzer.write_text("print('fixed')\n", encoding="utf-8")
    library.write_text("MIN_SIZE = 3\n", encoding="utf-8")

    def source_identity(path: Path) -> dict[str, object]:
        stat = path.stat()
        return {
            **artifact(path),
            "device": stat.st_dev,
            "inode": stat.st_ino,
            "mode": oct(stat.st_mode & 0o777),
            "mtime_ns": stat.st_mtime_ns,
            "ctime_ns": stat.st_ctime_ns,
        }

    during = {
        "analyzer": source_identity(analyzer),
        "focal_alt_cluster_lib": source_identity(library),
    }
    manifest_path = tmp_path / "tumor_ref_run_manifest.json"
    write_json(
        manifest_path,
        {
            "schema_name": "intersubmod.all_ssnv_tumor_ref_controls.run_manifest",
            "schema_version": "1.0.0",
            "finished_at_utc": "2026-07-17T00:00:00+00:00",
            "command": [str(analyzer.resolve()), "--fixed"],
            "source_code": {
                role: {
                    key: source[key] for key in ("path", "size_bytes", "sha256")
                }
                for role, source in during.items()
            },
            "pass": True,
        },
    )
    snapshot_path = tmp_path / "snapshot.json"
    write_json(
        snapshot_path,
        {
            "process": {
                "cmdline": [
                    "/usr/bin/python3",
                    str(analyzer.resolve()),
                    "--fixed",
                ]
            },
            "expected_command_fragment": str(analyzer.resolve()),
            "source_identity_during_execution": during,
            "snapshot_creator_source_identity": source_identity(analyzer),
            "pass": True,
        },
    )
    creator = source_identity(analyzer)
    verifier = source_identity(MODULE.TUMOR_REF_SOURCE_IDENTITY_VERIFIER)
    checks = {key: True for key in MODULE.TUMOR_REF_SOURCE_IDENTITY_TRUE_CHECKS}
    checks["process_start_clock_tolerance_seconds"] = 2.0
    receipt_path = tmp_path / "source_identity_receipt.json"
    write_json(
        receipt_path,
        {
            "schema_name": "intersubmod.retrospective_running_source_identity.receipt",
            "schema_version": MODULE.TUMOR_REF_SOURCE_IDENTITY_SCHEMA_VERSION,
            "created_at_utc": "2026-07-17T00:00:01+00:00",
            "task_type": "B_comprehensive_validation",
            "audit_class": "bounded_retrospective_source_file_identity",
            "snapshot": artifact(snapshot_path),
            "tumor_ref_run_manifest": artifact(manifest_path),
            "source_identity_during_execution": during,
            "source_identity_after_execution": during,
            "snapshot_creator_source_identity": creator,
            "snapshot_creator_source_identity_after_execution": creator,
            "post_run_verifier_source_identity": verifier,
            "auditor_source_identity_after_execution": verifier,
            "command_binding": {
                "live_python_launcher_token": "/usr/bin/python3",
                "manifest_script_token": str(analyzer.resolve()),
                "manifest_script_token_mode": "absolute",
                "attested_analyzer_path": str(analyzer.resolve()),
                "live_after_launcher_exactly_equals_manifest": True,
                "relative_token_rejects_dot_and_parent_segments": True,
                "repo_relative_token_must_equal_attested_source_relative_to_repo_root": True,
            },
            "checks": checks,
            "pass_semantics": (
                "Named producer source files were observed during execution, predated the run, "
                "remained identity-equal afterward, and match the passing producer manifest."
            ),
            "limitation": (
                "Retrospective bounded source-file attestation only; not a prelaunch lock or a "
                "complete interpreter, package, kernel, hardware, or environment attestation."
            ),
            "pass": True,
        },
    )
    attestation, paths = MODULE.load_tumor_ref_source_identity_attestation(
        receipt_path,
        manifest_path,
        json.loads(manifest_path.read_text(encoding="utf-8")),
    )
    assert attestation["release_gate_pass"] is True
    assert attestation["publishable_task_b_release"] is True
    assert attestation["source_roles"] == ["analyzer", "focal_alt_cluster_lib"]
    assert paths == [receipt_path, snapshot_path]

    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["post_run_verifier_source_identity"] = source_identity(library)
    receipt["auditor_source_identity_after_execution"] = source_identity(library)
    write_json(receipt_path, receipt)
    with pytest.raises(MODULE.ContractError, match="trusted v2 verifier"):
        MODULE.load_tumor_ref_source_identity_attestation(
            receipt_path,
            manifest_path,
            json.loads(manifest_path.read_text(encoding="utf-8")),
        )

    receipt["post_run_verifier_source_identity"] = verifier
    receipt["auditor_source_identity_after_execution"] = verifier
    receipt["command_binding"]["manifest_script_token"] = "analyzer.py"
    write_json(receipt_path, receipt)
    with pytest.raises(MODULE.ContractError, match="command-binding"):
        MODULE.load_tumor_ref_source_identity_attestation(
            receipt_path,
            manifest_path,
            json.loads(manifest_path.read_text(encoding="utf-8")),
        )

    receipt["command_binding"]["manifest_script_token"] = str(analyzer.resolve())
    write_json(receipt_path, receipt)

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source_code"]["analyzer"]["sha256"] = "0" * 64
    write_json(manifest_path, manifest)
    with pytest.raises(MODULE.ContractError, match="run manifest"):
        MODULE.load_tumor_ref_source_identity_attestation(
            receipt_path,
            manifest_path,
            manifest,
        )


def test_duplicate_candidate_and_exact_key_drift_fail_before_output(tmp_path: Path) -> None:
    duplicate_bundle = make_fixture(tmp_path / "duplicate", with_candidates=True)
    duplicate_path = duplicate_bundle.cooccurrence_sites
    with gzip.open(duplicate_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    rows.append(dict(rows[0]))
    write_tsv(duplicate_path, rows, fields)
    duplicate_output = tmp_path / "duplicate_output"
    with pytest.raises(MODULE.ContractError, match="Duplicate cooccurrence site table"):
        MODULE.build_outputs(
            duplicate_bundle, duplicate_output, expected_screen_sites=8
        )
    assert not duplicate_output.exists()

    drift_bundle = make_fixture(tmp_path / "drift", with_candidates=True)
    strict_path = drift_bundle.strict_sites
    with gzip.open(strict_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        strict = list(reader)
    strict[0]["alt"] = "G"
    write_tsv(strict_path, strict, fields)
    drift_output = tmp_path / "drift_output"
    with pytest.raises(MODULE.ContractError, match="key-set mismatch"):
        MODULE.build_outputs(drift_bundle, drift_output, expected_screen_sites=8)
    assert not drift_output.exists()


def test_g2_requires_two_effect_supported_markers_in_same_complete_read_set(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    with gzip.open(bundle.cooccurrence_sites, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    target = next(row for row in rows if row["sample"] == "HCC1395")
    target["joint_signature_n_complete_marker_effect_supported"] = "1"
    target["joint_signature_complete_marker_effect_supported_positions"] = json.dumps([120])
    support = json.loads(target["joint_signature_complete_marker_support"])
    support[1]["cramers_v"] = 0.1
    support[1]["effect_gate_pass"] = False
    target["joint_signature_complete_marker_support"] = json.dumps(support)
    target["n_same_complete_read_effect_supported_top_pair_by"] = "1"
    target["same_complete_read_effect_supported_top_pair_by_positions"] = json.dumps([120])
    target[
        "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp"
    ] = "1"
    target[
        "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp"
    ] = json.dumps([120])
    target[MODULE.FORMAL_SELECTION_COLUMN] = "false"
    target[MODULE.DEFAULT_SELECTION_COLUMN] = "false"
    target[MODULE.BY_SELECTION_COLUMN] = "false"
    write_tsv(bundle.cooccurrence_sites, rows, fields)
    manifest = MODULE.load_manifest(bundle.manifest, 8)
    screen = MODULE.load_screen(bundle.screen_sites, manifest, 8)
    _, _, g2_sites = MODULE.load_cooccurrence_sites(bundle.cooccurrence_sites, screen)
    assert ("HCC1395", "chr1", 100, "A", "C") not in g2_sites
    assert ("HCC1395_DORADO", "chr1", 100, "A", "C") in g2_sites


def test_g2_joint_signature_pass_is_recomputed_from_conditional_p(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    rows = read_gzip_tsv(bundle.cooccurrence_sites)
    fields = list(rows[0])
    target = next(row for row in rows if row["sample"] == "HCC1395")
    target["joint_signature_p_conditional_perm"] = "0.90"
    target["joint_signature_sensitivity_pass"] = "true"
    write_tsv(bundle.cooccurrence_sites, rows, fields)
    manifest = MODULE.load_manifest(bundle.manifest, 8)
    screen = MODULE.load_screen(bundle.screen_sites, manifest, 8)
    with pytest.raises(MODULE.ContractError, match="sensitivity-pass drift"):
        MODULE.load_cooccurrence_sites(bundle.cooccurrence_sites, screen)


def test_g2_complete_marker_effect_gate_is_recomputed(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    rows = read_gzip_tsv(bundle.cooccurrence_sites)
    fields = list(rows[0])
    target = next(row for row in rows if row["sample"] == "HCC1395")
    support = json.loads(target["joint_signature_complete_marker_support"])
    support[0]["cramers_v"] = 0.10
    support[0]["effect_gate_pass"] = True
    target["joint_signature_complete_marker_support"] = json.dumps(support)
    write_tsv(bundle.cooccurrence_sites, rows, fields)
    manifest = MODULE.load_manifest(bundle.manifest, 8)
    screen = MODULE.load_screen(bundle.screen_sites, manifest, 8)
    with pytest.raises(MODULE.ContractError, match="marker effect gate drift"):
        MODULE.load_cooccurrence_sites(bundle.cooccurrence_sites, screen)


def test_g2_top_marker_set_is_reconciled_to_pair_rows(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    rows = read_gzip_tsv(bundle.cooccurrence_pairs)
    fields = list(rows[0])
    target = next(
        row
        for row in rows
        if row["sample"] == "HCC1395" and int(row["partner_pos"]) == 150
    )
    target["top_coverage_marker"] = "false"
    write_tsv(bundle.cooccurrence_pairs, rows, fields)
    refresh_cooccurrence_chain(bundle)
    output_dir = tmp_path / "top_marker_drift"
    with pytest.raises(MODULE.ContractError, match="top-marker reconciliation"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_primary_artifact_audit_internal_time_order_is_fail_closed() -> None:
    pre = {
        "started_at_utc": "2026-07-15T00:01:00+00:00",
        "finished_at_utc": "2026-07-15T00:00:00+00:00",
        "created_at_utc": "2026-07-15T00:00:00+00:00",
    }
    post = {
        "started_at_utc": "2026-07-15T00:03:00+00:00",
        "finished_at_utc": "2026-07-15T00:04:00+00:00",
        "created_at_utc": "2026-07-15T00:04:00+00:00",
    }
    with pytest.raises(MODULE.ContractError, match="time window is reversed"):
        MODULE.validate_primary_artifact_audit_window(pre, post, {})


def test_matched_normal_per_sample_artifact_digest_mismatch_fails(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    assert bundle.matched_normal_summary is not None
    assert bundle.matched_normal_receipt is not None
    summary = json.loads(bundle.matched_normal_summary.read_text())
    receipt = json.loads(bundle.matched_normal_receipt.read_text())
    sample = next(iter(summary["paired_artifact_identity_validation"]["samples"]))
    summary["paired_artifact_identity_validation"]["samples"][sample][
        "artifact_set_sha256"
    ] = "f" * 64
    receipt["paired_artifact_identity_validation"] = deepcopy(
        summary["paired_artifact_identity_validation"]
    )
    write_json(bundle.matched_normal_summary, summary)
    receipt["outputs"]["summary"] = artifact(bundle.matched_normal_summary)
    write_json(bundle.matched_normal_receipt, receipt)
    output_dir = tmp_path / "matched_digest_drift"
    with pytest.raises(MODULE.ContractError, match="validation failed for"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_zero_candidate_requires_explicit_na_receipt(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=False)
    output_dir = tmp_path / "zero_final"
    outputs = MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    dataset = json.loads(outputs["final_report_dataset.json"].read_text())
    assert dataset["counts"]["g2_candidates"] == 0
    assert dataset["counts"]["candidate_catalog_rows"] == 0
    assert metric(dataset, "R1") == {
        "numerator": 0,
        "denominator": 0,
        "ratio": None,
        "not_evaluable": 0,
        "not_run": 8,
        "population": 8,
        "denominator_definition": MODULE.CLAIM_RULES["R1"]["denominator"],
    }
    assert read_tsv(outputs["candidate_catalog.tsv"]) == []
    assert read_tsv(outputs["candidate_witness_pairs.tsv"]) == []
    assert read_tsv_header(outputs["candidate_catalog.tsv"]) == MODULE.CANDIDATE_FIELDS
    assert read_tsv_header(outputs["candidate_witness_pairs.tsv"]) == MODULE.PAIR_DETAIL_FIELDS
    assert tuple(
        row["claim_id"] for row in read_tsv(outputs["claim_ladder.tsv"])
    ) == MODULE.CLAIM_IDS

    receipt_path = bundle.strict_dir / "run_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["reason"] = "NOT_RUN"
    write_json(receipt_path, receipt)
    bad_output = tmp_path / "bad_zero_final"
    with pytest.raises(MODULE.ContractError, match="reason must be"):
        MODULE.build_outputs(bundle, bad_output, expected_screen_sites=8)
    assert not bad_output.exists()


@pytest.mark.parametrize(
    ("field_path", "replacement", "error"),
    [
        (
            ("pass_semantics",),
            "scientific_negative_evidence",
            "Strict N/A pass semantics drift",
        ),
        (
            ("is_negative_result",),
            True,
            "is_negative_result=false",
        ),
        (
            ("scientific_interpretation", "is_negative_result"),
            True,
            "is_negative_result=false",
        ),
    ],
)
def test_strict_zero_candidate_cannot_be_recast_as_negative_evidence(
    tmp_path: Path,
    field_path: tuple[str, ...],
    replacement: object,
    error: str,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=False)
    receipt_path = bundle.strict_dir / "run_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    target = receipt
    for field_name in field_path[:-1]:
        target = target[field_name]
    target[field_path[-1]] = replacement
    write_json(receipt_path, receipt)
    receipt_path.chmod(0o444)

    with pytest.raises(MODULE.ContractError, match=error):
        MODULE.build_outputs(
            bundle,
            tmp_path / "negative_recast_output",
            expected_screen_sites=8,
        )


def test_strict_execution_status_is_required(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    summary_path = bundle.strict_summary
    summary = json.loads(summary_path.read_text())
    del summary["execution_status"]
    write_json(summary_path, summary)
    output_dir = tmp_path / "missing_execution_status"
    with pytest.raises(MODULE.ContractError, match="execution_status"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_strict_descriptive_schema_is_required_but_not_a_gate(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    strict_path = bundle.strict_sites
    with gzip.open(strict_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    fields.remove("strict_postselection_by_q_descriptive")
    for row in rows:
        row.pop("strict_postselection_by_q_descriptive")
    write_tsv(strict_path, rows, fields)

    output_dir = tmp_path / "missing_descriptive_schema"
    with pytest.raises(
        MODULE.ContractError, match="strict_postselection_by_q_descriptive"
    ):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_same_pair_witness_cannot_be_stitched_across_pairs(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    with gzip.open(bundle.cooccurrence_pairs, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    for row in rows:
        set_four_state_fields(
            row,
            {"RR": 4, "AR": 4, "RA": 4, "AA": 4, "O": 1, "X": 0},
        )
    write_tsv(bundle.cooccurrence_pairs, rows, fields)
    refresh_cooccurrence_chain(bundle)

    output_dir = tmp_path / "same_pair_final"
    outputs = MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    candidates = {
        row["sample"]: row for row in read_tsv(outputs["candidate_catalog.tsv"])
    }
    assert candidates["HCC1395"]["n_same_pair_four_state_witnesses"] == "0"
    assert candidates["HCC1395"]["b1_status"] == "FAIL"
    assert candidates["HCC1395"]["b1_reason"] == (
        "PRESPECIFIED_G1_PAIR_INCOMPATIBLE_OR_COMPLEX_UNDER_FIXED_ERROR_MODEL"
    )
    assert json.loads(outputs["final_report_dataset.json"].read_text())[
        "technical_replication"
    ]["numerator"] == 1


def test_b1_uses_one_prespecified_pair_and_preserves_fixed_error_not_evaluable(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    with gzip.open(bundle.cooccurrence_pairs, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    for row in rows:
        if row["partner_pos"] == "120":
            set_four_state_fields(
                row,
                {"RR": 8, "AR": 2, "RA": 1, "AA": 9, "O": 0, "X": 0},
            )
    write_tsv(bundle.cooccurrence_pairs, rows, fields)
    refresh_cooccurrence_chain(bundle)

    outputs = MODULE.build_outputs(
        bundle, tmp_path / "prespecified_pair_final", expected_screen_sites=8
    )
    candidates = {
        row["sample"]: row for row in read_tsv(outputs["candidate_catalog.tsv"])
    }
    primary = candidates["HCC1395"]
    assert primary["b1_status"] == "NOT_EVALUABLE"
    assert primary["b1_reason"] == (
        "PRESPECIFIED_G1_PAIR_NOT_IDENTIFIABLE_FIXED_ERROR_CEILING"
    )
    assert primary["n_same_pair_four_state_witnesses"] == "0"
    assert primary["n_four_state_compatible_formal_pair_opportunities"] == "1"
    assert json.loads(primary["b1_prespecified_pair_key"])[6] == 120
    assert primary["b1_uses_posthoc_compatible_pair_search"] == "false"


@pytest.mark.parametrize(
    ("mutation", "error_match"),
    [
        ("global_by_q", "global adjustment"),
        ("family_membership", "Endpoint-A global FDR family membership drift"),
        ("callability_gate", "Callability gate drift"),
        ("formal_false", "Formal G1 gate drift"),
        ("conditional_permutations", "Conditional permutation count drift"),
        ("four_state_upper_bound", "four-state.*upper_bound drift"),
    ],
)
def test_pair_scientific_gates_are_recomputed_fail_closed(
    tmp_path: Path, mutation: str, error_match: str
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    with gzip.open(bundle.cooccurrence_pairs, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    target = rows[0]
    if mutation == "global_by_q":
        target["endpoint_a_q_global_by"] = 0.9
    elif mutation == "family_membership":
        target["endpoint_a_global_fdr_family_status"] = (
            "ELIGIBLE_M2_EXACT_NOT_IDENTIFIABLE"
        )
    elif mutation == "callability_gate":
        target["callability_gate_pass"] = False
    elif mutation == "formal_false":
        for field_name in (
            "endpoint_a_formal_pair_by_confirmed",
            "endpoint_a_confirmed_association",
            "endpoint_a_confirmed_by_sensitivity",
        ):
            target[field_name] = False
    elif mutation == "conditional_permutations":
        target["endpoint_a_permutations"] = 998
    elif mutation == "four_state_upper_bound":
        target["endpoint_b_focal_ancestor_violation_upper_bound"] = 0.5
    else:
        raise AssertionError(mutation)
    write_tsv(bundle.cooccurrence_pairs, rows, fields)
    refresh_cooccurrence_chain(bundle)
    output_dir = tmp_path / f"tampered_{mutation}"
    with pytest.raises(MODULE.ContractError, match=error_match):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_joint_signature_global_by_adjustment_is_recomputed_fail_closed(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    rows = read_gzip_tsv(bundle.cooccurrence_sites)
    fields = list(rows[0])
    target = next(row for row in rows if row["sample"] == "HCC1395")
    target["joint_signature_q_global_by"] = "0.90"
    write_tsv(bundle.cooccurrence_sites, rows, fields)
    refresh_cooccurrence_chain(bundle)

    output_dir = tmp_path / "tampered_joint_global_by"
    with pytest.raises(MODULE.ContractError, match="joint-signature global adjustment"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_raw_identity_sparse_duplicate_is_independently_reconciled(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    install_single_raw_identity_duplicate(bundle)
    outputs = MODULE.build_outputs(
        bundle,
        tmp_path / "raw_duplicate_positive",
        expected_screen_sites=8,
    )
    assert outputs["final_report_dataset.json"].is_file()


def test_raw_identity_sparse_duplicate_class_tag_tamper_fails(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    install_single_raw_identity_duplicate(bundle)
    rows = read_gzip_tsv(bundle.cooccurrence_raw_identity_duplicates)
    rows[0]["classification"] = "EXACT_DUPLICATE"
    write_tsv(
        bundle.cooccurrence_raw_identity_duplicates,
        rows,
        fields=list(MODULE.RAW_IDENTITY_DUPLICATE_COLUMNS),
    )
    refresh_cooccurrence_chain(bundle)
    with pytest.raises(MODULE.ContractError, match="class/tag drift"):
        MODULE.build_outputs(
            bundle,
            tmp_path / "raw_duplicate_tampered",
            expected_screen_sites=8,
        )


def test_cooccurrence_preflight_source_role_superset_fails(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    preflight = json.loads(bundle.cooccurrence_preflight.read_text())
    unexpected = artifact(bundle.manifest)
    for field_name in ("source_identity_before", "source_identity_after"):
        preflight["code"][field_name]["unexpected"] = unexpected
    for field_name in ("source_modes_before", "source_modes_after"):
        preflight["code"][field_name]["unexpected"] = "0o444"
    write_json(bundle.cooccurrence_preflight, preflight)
    refresh_cooccurrence_chain(bundle)
    with pytest.raises(MODULE.ContractError, match="source identity map drift"):
        MODULE.build_outputs(
            bundle,
            tmp_path / "preflight_role_superset",
            expected_screen_sites=8,
        )


def test_cooccurrence_oracle_failure_is_not_accepted(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    oracle = json.loads(bundle.cooccurrence_oracle_cases.read_text())
    oracle["focal_cases"][0]["partner_window_oracle_pass"] = False
    write_json(bundle.cooccurrence_oracle_cases, oracle)
    refresh_cooccurrence_chain(bundle)
    with pytest.raises(MODULE.ContractError, match="focal oracle failed"):
        MODULE.build_outputs(
            bundle,
            tmp_path / "oracle_failure",
            expected_screen_sites=8,
        )


def test_g2_evaluable_requires_exact_markers_and_executable_joint_permutation() -> None:
    key = ("HCC1395", "chr1", 100, "A", "C")
    pair_template = {
        "top_coverage_marker": True,
        "endpoint_a_testable": True,
        "endpoint_a_global_fdr_family_status": "ELIGIBLE_M2_EXACT_FAMILY",
        "endpoint_a_p_fixed_margins_exact": 0.01,
    }
    pairs = [
        {**pair_template, "partner_pos": 120},
        {**pair_template, "partner_pos": 150},
    ]
    site = {
        "joint_signature_testable": True,
        "joint_signature_permutable": True,
        "joint_signature_p_conditional_perm": 0.01,
        "joint_signature_permutations": 999,
    }
    data = SimpleNamespace(
        cooccurrence_rows={key: site},
        pair_data=SimpleNamespace(by_focal={key: pairs}),
    )
    assert MODULE.g2_evaluable(data, key) is True
    pairs[1]["endpoint_a_global_fdr_family_status"] = "ELIGIBLE_M2_EXACT_NOT_IDENTIFIABLE"
    assert MODULE.g2_evaluable(data, key) is False
    pairs[1]["endpoint_a_global_fdr_family_status"] = "ELIGIBLE_M2_EXACT_FAMILY"
    site["joint_signature_permutations"] = 998
    assert MODULE.g2_evaluable(data, key) is False


def test_biological_claim_ladder_uses_hcc1395_primary_not_optimistic_union() -> None:
    primary_key = ("HCC1395", "chr1", 100, "A", "C")
    dorado_key = ("HCC1395_DORADO", "chr1", 100, "A", "C")
    statuses = {
        primary_key: {claim_id: "FAIL" for claim_id in MODULE.CLAIM_IDS},
        dorado_key: {claim_id: "PASS" for claim_id in MODULE.CLAIM_IDS},
    }
    rows = MODULE.build_claim_ladder_v2(None, statuses)
    for row in rows:
        assert row["biological_numerator"] == 0
        assert row["biological_denominator"] == 1
        assert row["hcc1395_primary_numerator"] == 0
        assert row["dorado_sensitivity_numerator"] == 1
        assert row["hcc1395_dorado_exact_site_overlap"] == 1
        assert row["hcc1395_dorado_pass_fail_discordant"] == 1


def test_receipt_input_hash_mismatch_fails_before_output(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    receipt = json.loads(bundle.cooccurrence_receipt.read_text())
    receipt["inputs"]["manifest"]["sha256"] = "0" * 64
    write_json(bundle.cooccurrence_receipt, receipt)
    output_dir = tmp_path / "hash_mismatch"
    with pytest.raises(MODULE.ContractError, match="SHA-256 mismatch"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_cooccurrence_receipt_requires_m2_gate_code_identity(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    receipt = json.loads(bundle.cooccurrence_receipt.read_text())
    receipt["code"].pop("m2_screen_gate")
    write_json(bundle.cooccurrence_receipt, receipt)
    output_dir = tmp_path / "missing_m2_gate_code"
    with pytest.raises(MODULE.ContractError, match="M2 screen gate code receipt lacks"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_cooccurrence_frozen_identity_policy_drift_fails_closed(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    receipt = json.loads(bundle.cooccurrence_receipt.read_text())
    receipt["frozen_manifest_input_identity_policy"]["artifacts"][0][
        "identity_mode"
    ] = "explicit_large_artifact_size_mtime_path_v1"
    write_json(bundle.cooccurrence_receipt, receipt)
    refresh_cooccurrence_chain(bundle)

    output_dir = tmp_path / "cooccurrence_identity_policy_drift"
    with pytest.raises(
        MODULE.ContractError,
        match="frozen-input identity policy is missing or drifted",
    ):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_pair_site_exact_reconciliation_fails_closed(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    with gzip.open(bundle.cooccurrence_sites, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    rows[0]["n_pair_rows_reconciled"] = "1"
    write_tsv(bundle.cooccurrence_sites, rows, fields)
    output_dir = tmp_path / "pair_reconciliation"
    with pytest.raises(MODULE.ContractError, match="site pair-count drift"):
        MODULE.build_outputs(bundle, output_dir, expected_screen_sites=8)
    assert not output_dir.exists()


def test_screen_latest_tag_row_and_aggregate_tampering_fail_closed(tmp_path: Path) -> None:
    row_bundle = make_fixture(tmp_path / "row", with_candidates=True)
    rows = read_gzip_tsv(row_bundle.screen_sites)
    fields = list(rows[0])
    rows[0]["latest_tag_reads_joined"] = str(int(rows[0]["n_reads_total"]) - 1)
    write_tsv(row_bundle.screen_sites, rows, fields)
    with pytest.raises(MODULE.ContractError, match="joined/read-row mismatch"):
        MODULE.build_outputs(
            row_bundle, tmp_path / "row_output", expected_screen_sites=8
        )

    aggregate_bundle = make_fixture(tmp_path / "aggregate", with_candidates=True)
    summary = json.loads(aggregate_bundle.screen_summary.read_text())
    summary["latest_hp_ps_terminal_join_audit"]["n_reads_tsv_site_rows"] += 1
    write_json(aggregate_bundle.screen_summary, summary)
    receipt = json.loads(aggregate_bundle.screen_receipt.read_text())
    receipt["outputs"]["summary"] = artifact(aggregate_bundle.screen_summary)
    write_json(aggregate_bundle.screen_receipt, receipt)
    with pytest.raises(MODULE.ContractError, match="terminal join aggregate mismatch"):
        MODULE.build_outputs(
            aggregate_bundle,
            tmp_path / "aggregate_output",
            expected_screen_sites=8,
        )

    receipt_bundle = make_fixture(tmp_path / "receipt", with_candidates=True)
    receipt = json.loads(receipt_bundle.screen_receipt.read_text())
    receipt["counts"]["exact_hp_ps_site_read_joins"] += 1
    write_json(receipt_bundle.screen_receipt, receipt)
    with pytest.raises(MODULE.ContractError, match="receipt count mismatch"):
        MODULE.build_outputs(
            receipt_bundle, tmp_path / "receipt_output", expected_screen_sites=8
        )


def test_g2_uses_top_spaced_positions_and_accepts_legal_false_candidate(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    rows = read_gzip_tsv(bundle.cooccurrence_sites)
    fields = list(rows[0])
    target = next(row for row in rows if row["sample"] == "HCC1395")
    target.update(
        {
            "pair_by_confirmed_positions": json.dumps([120, 130]),
            "n_spatially_separated_pair_by_20bp": "1",
            "spatially_separated_pair_by_positions_20bp": json.dumps([120]),
            "top_marker_positions": json.dumps([120, 130]),
            "n_top_marker_pair_by_confirmed": "2",
            "top_marker_pair_by_confirmed_positions": json.dumps([120, 130]),
            "joint_signature_complete_marker_support": json.dumps(
                [
                    {
                        "position": 120,
                        "testable": True,
                        "reason": "TESTABLE",
                        "n_informative": 12,
                        "cramers_v": 0.7,
                        "delta_alt_fraction": 0.8,
                        "effect_gate_pass": True,
                    },
                    {
                        "position": 130,
                        "testable": True,
                        "reason": "TESTABLE",
                        "n_informative": 12,
                        "cramers_v": 0.1,
                        "delta_alt_fraction": 0.8,
                        "effect_gate_pass": False,
                    },
                ]
            ),
            "joint_signature_n_complete_marker_effect_supported": "1",
            "joint_signature_complete_marker_effect_supported_positions": json.dumps([120]),
            "n_same_complete_read_effect_supported_top_pair_by": "1",
            "same_complete_read_effect_supported_top_pair_by_positions": json.dumps([120]),
            "n_spatially_separated_same_complete_read_effect_supported_top_pair_by_20bp": "1",
            "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp": json.dumps(
                [120]
            ),
            MODULE.FORMAL_SELECTION_COLUMN: "false",
            MODULE.DEFAULT_SELECTION_COLUMN: "false",
            MODULE.BY_SELECTION_COLUMN: "false",
        }
    )
    write_tsv(bundle.cooccurrence_sites, rows, fields)
    manifest = MODULE.load_manifest(bundle.manifest, 8)
    screen = MODULE.load_screen(bundle.screen_sites, manifest, 8)
    _, _, g2_sites = MODULE.load_cooccurrence_sites(bundle.cooccurrence_sites, screen)
    assert ("HCC1395", "chr1", 100, "A", "C") not in g2_sites


def test_four_state_full_producer_enum_is_required(tmp_path: Path) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    rows = read_gzip_tsv(bundle.cooccurrence_pairs)
    fields = list(rows[0])
    rows[0]["endpoint_b_relation_compatibility"] = "FOCAL_ANCESTOR_COMPATIBLE"
    write_tsv(bundle.cooccurrence_pairs, rows, fields)
    with pytest.raises(MODULE.ContractError, match="outside the producer enum"):
        MODULE.build_outputs(bundle, tmp_path / "output", expected_screen_sites=8)


def test_zero_candidate_matched_normal_explicit_na_and_native_analysis_are_validated(
    tmp_path: Path,
) -> None:
    na_bundle = add_matched_normal_na_fixture(
        make_fixture(tmp_path / "na", with_candidates=False), tmp_path / "na"
    )
    na_outputs = MODULE.build_outputs(
        na_bundle, tmp_path / "na_output", expected_screen_sites=8
    )
    na_dataset = json.loads(na_outputs["final_report_dataset.json"].read_text())
    assert na_dataset["validations"]["matched_normal_status"] == (
        "NOT_APPLICABLE_VALIDATED_RECEIPT"
    )
    assert na_dataset["validations"]["matched_normal_candidate_keys_exact"] is True

    native_bundle = add_zero_row_matched_normal_analysis_fixture(
        make_fixture(tmp_path / "native", with_candidates=False), tmp_path / "native"
    )
    native_outputs = MODULE.build_outputs(
        native_bundle, tmp_path / "native_output", expected_screen_sites=8
    )
    native_dataset = json.loads(native_outputs["final_report_dataset.json"].read_text())
    assert native_dataset["validations"]["matched_normal_status"] == (
        "NOT_APPLICABLE_VALIDATED_NATIVE_ZERO_ROWS"
    )


@pytest.mark.parametrize("target_name", ["summary", "receipt"])
def test_matched_normal_not_evaluable_cannot_be_recast_as_negative_evidence(
    tmp_path: Path, target_name: str
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    assert bundle.matched_normal_dir is not None
    summary_path = bundle.matched_normal_summary
    receipt_path = bundle.matched_normal_receipt
    assert summary_path is not None and receipt_path is not None

    target_path = summary_path if target_name == "summary" else receipt_path
    target = json.loads(target_path.read_text())
    target["not_evaluable_is_negative_result"] = True
    write_json(target_path, target)
    if target_name == "summary":
        receipt = json.loads(receipt_path.read_text())
        receipt["outputs"]["summary"] = artifact(summary_path)
        write_json(receipt_path, receipt)

    with pytest.raises(MODULE.ContractError, match="not_evaluable_is_negative_result=false"):
        MODULE.build_outputs(
            bundle,
            tmp_path / "negative_recast_output",
            expected_screen_sites=8,
        )


def test_matched_normal_runner_wrong_release_authority_fails_closed(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    assert bundle.matched_normal_receipt is not None
    analysis_receipt = json.loads(bundle.matched_normal_receipt.read_text())
    paired_path = Path(analysis_receipt["inputs"]["paired_run_receipt"]["path"])
    paired = json.loads(paired_path.read_text())
    paired["source_authority"] = {"authority_id": "WRONG_AUTHORITY", "pass": True}
    write_json(paired_path, paired)
    analysis_receipt["inputs"]["paired_run_receipt"] = artifact(paired_path)
    write_json(bundle.matched_normal_receipt, analysis_receipt)
    refresh_post_audit_consumer_identities(bundle)

    with pytest.raises(MODULE.ContractError, match="paired runner source authority drift"):
        MODULE.build_outputs(bundle, tmp_path / "output", expected_screen_sites=8)


def test_matched_normal_analyzer_wrong_producer_code_fails_closed(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    assert bundle.matched_normal_receipt is not None
    receipt = json.loads(bundle.matched_normal_receipt.read_text())
    receipt["code"]["matched_normal_analyzer"]["sha256"] = "0" * 64
    write_json(bundle.matched_normal_receipt, receipt)

    with pytest.raises(MODULE.ContractError, match="analyzer producer source identity drift"):
        MODULE.build_outputs(bundle, tmp_path / "output", expected_screen_sites=8)


def test_native_cn_receipt_annotations_and_c1_ceiling(tmp_path: Path) -> None:
    bundle = add_native_cn_fixture(
        make_fixture(tmp_path, with_candidates=True), tmp_path
    )
    outputs = MODULE.build_outputs(bundle, tmp_path / "output", expected_screen_sites=8)
    candidates = {
        row["sample"]: row for row in read_tsv(outputs["candidate_catalog.tsv"])
    }
    hcc = candidates["HCC1395"]
    assert hcc["cn_status"] == "AVAILABLE_EXACT_SEGMENT"
    assert hcc["savana_total_cn_raw"] == "2.2"
    assert hcc["pyclone_vi_cellular_prevalence"] == "0.42"
    assert hcc["c1_status"] == "NOT_EVALUABLE"
    assert hcc["c1_reason"] == (
        "FOCAL_ONLY_CN_PYCLONE_ANNOTATION_NO_PRESPECIFIED_JOINT_WITNESS_PAIR_MODEL"
    )
    assert MODULE.cn_ccf_claim_status(
        {
            "sample": "HCC1395_DORADO",
            "cn_status": "SHARED_CN_SENSITIVITY",
        }
    )[0] == "NOT_EVALUABLE"
    assert MODULE.cn_ccf_claim_status(
        {"sample": "COLO829", "cn_status": "BLOCKED_CN_MISFIT"}
    ) == ("NOT_EVALUABLE", "BLOCKED_CN_MISFIT_NO_CN2_IMPUTATION")


def test_native_cn_command_missing_python_isolation_fails_closed(tmp_path: Path) -> None:
    bundle = add_native_cn_fixture(
        make_fixture(tmp_path, with_candidates=True), tmp_path
    )
    assert bundle.cn_ccf_annotations is not None
    receipt_path = bundle.cn_ccf_annotations / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["command"].remove("-I")
    write_json(receipt_path, receipt)

    with pytest.raises(MODULE.ContractError, match="annotator producer command is not canonical"):
        MODULE.build_outputs(bundle, tmp_path / "output", expected_screen_sites=8)


def test_native_cn_zero_row_receipt_and_authority_hash_tamper(tmp_path: Path) -> None:
    zero_bundle = add_native_cn_fixture(
        make_fixture(tmp_path / "zero", with_candidates=False),
        tmp_path / "zero",
        zero_rows=True,
    )
    outputs = MODULE.build_outputs(
        zero_bundle, tmp_path / "zero_output", expected_screen_sites=8
    )
    dataset = json.loads(outputs["final_report_dataset.json"].read_text())
    assert dataset["cn_ccf_annotation"]["status"] == (
        "NOT_APPLICABLE_VALIDATED_NATIVE_ZERO_ROWS"
    )
    assert dataset["cn_ccf_annotation"]["c1_formed"] is False

    tampered_bundle = add_native_cn_fixture(
        make_fixture(tmp_path / "tampered", with_candidates=True),
        tmp_path / "tampered",
    )
    receipt = json.loads(
        (tampered_bundle.cn_ccf_annotations / "receipt.json").read_text()
    )
    source_path = Path(next(iter(receipt["authority"]["all_source_hashes"])))
    source_path.write_text("tampered\n", encoding="utf-8")
    with pytest.raises(MODULE.ContractError, match="authority source SHA-256 mismatch"):
        MODULE.build_outputs(
            tampered_bundle, tmp_path / "tampered_output", expected_screen_sites=8
        )


@pytest.mark.parametrize("drift", ["code", "command", "parameters"])
def test_strict_receipt_must_bind_authorized_source_command_and_parameters(
    tmp_path: Path, drift: str
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    receipt = json.loads(bundle.strict_receipt.read_text())
    if drift == "code":
        receipt["code"]["strict_producer"]["sha256"] = "0" * 64
    elif drift == "command":
        receipt["command"][-1] = "0.7"
    else:
        receipt["parameters"]["permutations_per_seed_per_null"] = 998
    write_json(bundle.strict_receipt, receipt)
    bundle.strict_receipt.chmod(0o444)

    with pytest.raises(MODULE.ContractError, match="Strict receipt"):
        MODULE.build_outputs(
            bundle,
            tmp_path / f"strict_{drift}_output",
            expected_screen_sites=8,
        )


def test_primary_audit_self_consistent_source_substitution_is_rejected(
    tmp_path: Path,
) -> None:
    bundle = make_fixture(tmp_path, with_candidates=True)
    audit = json.loads(bundle.primary_artifact_audit_pre.read_text())
    audit["code"]["primary_artifact_auditor"]["sha256"] = "0" * 64
    audit["source_lock"]["source_identity_before"] = audit["code"]
    audit["source_lock"]["source_identity_after_compute"] = audit["code"]
    write_json(bundle.primary_artifact_audit_pre, audit)
    bundle.primary_artifact_audit_pre.chmod(0o444)

    with pytest.raises(MODULE.ContractError, match="producer source identity drift"):
        MODULE.build_outputs(
            bundle,
            tmp_path / "primary_source_substitution_output",
            expected_screen_sites=8,
        )


def test_cross_producer_contract_constants_and_filenames_do_not_drift() -> None:
    scripts = SCRIPT.parent
    screen = scripts / "analyze_all_ssnv_focal_alt_multigroup.py"
    manifest = scripts / "prepare_all_ssnv_manifest.py"
    tumor_ref = scripts / "analyze_all_ssnv_tumor_ref_controls.py"
    matched_analysis = scripts / "analyze_matched_normal_candidate_controls.py"
    matched_runner = scripts / "run_matched_normal_candidate_controls.py"
    cooccurrence = scripts / "analyze_methyl_ssnv_cooccurrence.py"
    strict = scripts / "run_strict_methyl_candidate_confirmation.py"
    cn = scripts / "annotate_candidate_cn_ccf.py"
    cooccurrence_lib = SCRIPT.parents[3] / "scripts" / "ssnv_cooccurrence_lib.py"
    assert MODULE.INPUT_MANIFEST_SCHEMA_VERSION == "1.0.0"
    assert MODULE.SCREEN_OUTPUT_SCHEMA_VERSION == ast_assignment(
        screen, "OUTPUT_SCHEMA_VERSION"
    )
    assert MODULE.ASSIGNMENT_SCHEMA_VERSION == "1.0.0"
    assert MODULE.TUMOR_REF_SCHEMA_VERSION == ast_assignment(tumor_ref, "SCHEMA_VERSION")
    assert MODULE.MATCHED_NORMAL_ANALYSIS_SCHEMA_VERSION == ast_assignment(
        matched_analysis, "SCHEMA_VERSION"
    )
    assert MODULE.MATCHED_NORMAL_PAIRED_RUNNER_SCHEMA_VERSION == "1.0.0"
    assert MODULE.COOCCURRENCE_SCHEMA_VERSION == ast_assignment(
        cooccurrence, "SCHEMA_VERSION"
    )
    assert MODULE.STRICT_FORMAL_SCHEMA_VERSION == ast_assignment(strict, "SCHEMA_VERSION")
    assert MODULE.CN_CCF_ANNOTATION_SCHEMA_VERSION == ast_assignment(
        cn, "RECEIPT_SCHEMA_VERSION"
    )
    assert MODULE.CN_CCF_OUTPUT_COLUMNS == ast_assignment(cn, "OUTPUT_COLUMNS")
    assert MODULE.COMPATIBLE_RELATIONS == ast_assignment(
        cooccurrence, "CROSS_PLATFORM_FOUR_STATE_COMPATIBLE_RELATIONS"
    )
    assert MODULE.FOUR_STATE_RELATIONS == assigned_string_values(
        cooccurrence_lib, "four_state_summary", "relation"
    )
    required_screen_fields = {
        "latest_tag_join_status",
        "latest_tag_rows_fetched",
        "latest_tag_rows_eligible",
        "latest_tag_reads_joined",
        "latest_tag_ps_present",
        "latest_tag_projection_multimatch_reads",
        "source_hp_replaced_reads",
        "n_reads_total",
        "latest_hp_ps_terminal_join_audit",
    }
    assert required_screen_fields <= producer_string_literals(screen)
    expected_filenames = {
        screen: {
            "all_ssnv_site_results.tsv.gz",
            "all_ssnv_stable_multigroup_read_assignments.jsonl.gz",
            "all_ssnv_summary.json",
            "run_manifest.json",
        },
        tumor_ref: {
            "all_ssnv_tumor_ref_control_site_results.tsv.gz",
            "all_ssnv_tumor_ref_control_summary.json",
            "run_manifest.json",
        },
        matched_analysis: {
            "matched_normal_candidate_controls.tsv.gz",
            "matched_normal_candidate_controls_summary.json",
            "run_receipt.json",
        },
        matched_runner: {"not_applicable_receipt.json", "run_receipt.json"},
        cooccurrence: {
            "methyl_ssnv_pair_results.tsv.gz",
            "methyl_ssnv_site_results.tsv.gz",
            "summary.json",
            "run_receipt.json",
        },
        strict: {
            "strict_methyl_candidate_confirmation_sites.tsv.gz",
            "strict_methyl_candidate_confirmation_summary.json",
            "not_applicable_receipt.json",
            "run_manifest.json",
        },
        cn: {"candidate_cn_ccf_annotations.tsv.gz", "receipt.json"},
    }
    for producer, filenames in expected_filenames.items():
        assert filenames <= producer_string_literals(producer)
    assert "schema_version\": \"1.0.0" in manifest.read_text(encoding="utf-8")
    cooccurrence_source = cooccurrence.read_text(encoding="utf-8")
    assert "and len(top_spaced_positions) >= 2" in cooccurrence_source
    assert '"pass": True' in cooccurrence_source


def test_cli_legacy_directory_aliases_remain_readable(tmp_path: Path) -> None:
    args = MODULE.build_parser().parse_args(
        [
            "--input-manifest",
            str(tmp_path / "manifest.json"),
            "--screen-output-dir",
            str(tmp_path / "screen"),
            "--cooccurrence-output-dir",
            str(tmp_path / "cooccurrence"),
            "--strict-output-dir",
            str(tmp_path / "strict"),
            "--tumor-ref-output-dir",
            str(tmp_path / "tumor_ref"),
            "--primary-artifact-audit-pre",
            str(tmp_path / "primary_pre.json"),
            "--primary-artifact-audit-post",
            str(tmp_path / "primary_post.json"),
            "--cooccurrence-preflight",
            str(tmp_path / "cooccurrence_preflight.json"),
            "--normal-dir",
            str(tmp_path / "normal"),
            "--output-dir",
            str(tmp_path / "final"),
        ]
    )
    assert args.manifest == tmp_path / "manifest.json"
    assert args.primary_artifact_audit_pre == tmp_path / "primary_pre.json"
    assert args.primary_artifact_audit_post == tmp_path / "primary_post.json"
    assert args.cooccurrence_preflight == tmp_path / "cooccurrence_preflight.json"
    assert args.matched_normal_dir == tmp_path / "normal"
