from __future__ import annotations

import re
import subprocess
from pathlib import Path


TOPIC_ROOT = Path(__file__).resolve().parents[1]
RUNNER = TOPIC_ROOT / "scripts" / "run_terminal_validation_chain.sh"
RELEASE_RUNNER = TOPIC_ROOT / "scripts" / "run_source_attested_release_chain.sh"
M2V5_RUNNER = TOPIC_ROOT / "scripts" / "run_m2v5_recovered_completion_chain.sh"
COOCCURRENCE_V6_RUNNER = TOPIC_ROOT / "scripts" / "run_cooccurrence_v6_source_locked.sh"
FINAL_RESULT_FINALIZER = TOPIC_ROOT / "scripts" / "finalize_task_b_result_release.py"
FINAL_DATASET_BUILDER = (
    TOPIC_ROOT / "scripts" / "build_all_ssnv_final_report_dataset.py"
)
COOCCURRENCE_RELEASE_FINALIZER = (
    TOPIC_ROOT / "scripts" / "finalize_cooccurrence_release_receipt.py"
)
COOCCURRENCE_ANALYZER = (
    TOPIC_ROOT / "scripts" / "analyze_methyl_ssnv_cooccurrence.py"
)
STRICT_PRODUCER = TOPIC_ROOT / "scripts" / "run_strict_methyl_candidate_confirmation.py"
MATCHED_NORMAL_RUNNER = TOPIC_ROOT / "scripts" / "run_matched_normal_candidate_controls.py"
MATCHED_NORMAL_ANALYZER = (
    TOPIC_ROOT / "scripts" / "analyze_matched_normal_candidate_controls.py"
)
CN_CCF_ANNOTATOR = TOPIC_ROOT / "scripts" / "annotate_candidate_cn_ccf.py"
PRIMARY_AUDITOR = TOPIC_ROOT / "scripts" / "audit_stable_primary_artifacts.py"
REPORT_BUILDER = TOPIC_ROOT / "scripts" / "build_all_ssnv_report_artifact.py"
SOURCE_AUTHORITY_VALIDATOR = TOPIC_ROOT / "scripts" / "release_source_authority.py"
ONE_TIME_SIGNER = TOPIC_ROOT / "audit_notes" / "run_one_time_ed25519_signer_v2.sh"


def runner_text() -> str:
    return RUNNER.read_text(encoding="utf-8")


def release_runner_text() -> str:
    return RELEASE_RUNNER.read_text(encoding="utf-8")


def m2v5_runner_text() -> str:
    return M2V5_RUNNER.read_text(encoding="utf-8")


def cooccurrence_v6_runner_text() -> str:
    return COOCCURRENCE_V6_RUNNER.read_text(encoding="utf-8")


def cooccurrence_analyzer_text() -> str:
    return COOCCURRENCE_ANALYZER.read_text(encoding="utf-8")


def strict_producer_text() -> str:
    return STRICT_PRODUCER.read_text(encoding="utf-8")


def test_terminal_runner_has_valid_bash_syntax() -> None:
    subprocess.run(["bash", "-n", str(RUNNER)], check=True)


def test_terminal_runner_pins_all_numeric_thread_backends() -> None:
    text = runner_text()
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
    ):
        assert f"export {variable}=1" in text


def test_terminal_runner_locks_full_scope_and_formal_g2_selection() -> None:
    text = runner_text()
    assert ".counts.processed_sites == 469849" in text
    assert 'FORMAL_SELECTION_COLUMN="multi_marker_molecular_haplotype_base_candidate"' in text
    assert text.count('--selection-column "${FORMAL_SELECTION_COLUMN}"') == 2
    assert text.count("--selection-value true") == 2
    assert "--allow-partial-scope" not in text


def test_terminal_runner_preserves_primary_consumer_receipt_order() -> None:
    text = runner_text()
    expected = (
        '--consumer-receipt "${COOCCURRENCE}/run_receipt.json"',
        '--consumer-receipt "${TUMOR_REF}/run_manifest.json"',
        '--consumer-receipt "${STRICT_RECEIPT}"',
        '--consumer-receipt "${MATCHED_RUN_RECEIPT}"',
        '--consumer-receipt "${MATCHED_ANALYSIS}/run_receipt.json"',
    )
    positions = [text.index(token) for token in expected]
    assert positions == sorted(positions)


def test_terminal_runner_uses_post_downstream_frozen_audit_and_layout_gates() -> None:
    text = runner_text()
    assert '--post-immutability-audit "${FROZEN_POST_TERMINAL}"' in text
    assert text.count(".pass == true and .overlapCount == 0") == 2
    assert 'require_json_pass "${PORTABLE_RECEIPT}"' in text


def test_terminal_runner_uses_terminal_v3_claim_contract_without_mutating_v2() -> None:
    text = runner_text()
    assert 'TERMINAL_CLAIM_CONTRACT="${TOPIC_ROOT}/claim-contract-v3.md"' in text
    assert '--claim-contract "${TERMINAL_CLAIM_CONTRACT}"' in text
    assert '--claim-contract "${TOPIC_ROOT}/claim-contract-v2.md"' not in text


def test_terminal_runner_portable_runtime_and_qa_python_exist() -> None:
    text = runner_text()
    plugin_match = re.search(r'^readonly PLUGIN_ROOT="([^"]+)"$', text, re.MULTILINE)
    qa_python_match = re.search(
        r'^readonly QA_PYTHON="([^"]+)"$', text, re.MULTILINE
    )
    assert plugin_match is not None
    assert qa_python_match is not None
    plugin_root = Path(plugin_match.group(1))
    script_root = plugin_root / "skills" / "build-report" / "scripts"
    for name in (
        "build_portable_artifact.mjs",
        "extract_portable_chart_svgs.mjs",
        "verify_portable_artifact.mjs",
    ):
        assert (script_root / name).is_file()
    subprocess.run(
        [qa_python_match.group(1), "-c", "import playwright"],
        check=True,
    )


def test_source_attested_release_runner_is_v4_fail_closed() -> None:
    subprocess.run(["bash", "-n", str(RELEASE_RUNNER)], check=True)
    text = release_runner_text()
    assert 'CLAIM_CONTRACT="${TOPIC_ROOT}/claim-contract-v5.md"' in text
    assert 'SOURCE_SNAPSHOT="${SOURCE_AUDIT_ROOT}/observed_during_execution.snapshot.json"' in text
    assert 'SOURCE_RECEIPT="${SOURCE_AUDIT_ROOT}/post_run_source_identity.receipt.json"' in text
    assert '--tumor-ref-source-identity-receipt "${SOURCE_RECEIPT}"' in text
    assert "grep -Fxq 'Terminal validation chain PASS'" in text
    assert 'require_absent "${path}"' in text
    assert text.count(".pass == true and .overlapCount == 0") == 2
    assert "claim-contract-v3.md" not in text


def test_m2v5_recovered_completion_runner_is_full_scope_and_fail_closed() -> None:
    subprocess.run(["bash", "-n", str(M2V5_RUNNER)], check=True)
    text = m2v5_runner_text()
    assert (
        "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
        "source_locked_command_parity" in text
    )
    assert "all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel" in text
    assert 'CLAIM_CONTRACT="${TOPIC_ROOT}/claim-contract-v5.md"' in text
    assert (
        'PRIMARY_PRE="${RESULT_ROOT}/stable_primary_artifact_audit.'
        'v6_strict_command_parity_pre_downstream.json"' in text
    )
    assert 'INDEPENDENT_M2_AUDIT="${RESULT_ROOT}/independent_m2_gate_recount.v3.json"' in text
    assert (
        'COOCCURRENCE_PREFLIGHT="${RESULT_ROOT}/cooccurrence_task_contract_preflight.'
        'v9_command_parity_full_runtime.json"' in text
    )
    assert "all_ssnv_final_report_dataset_v5_m2v5_source_attested" in text
    assert "all_ssnv_final_report_v5_m2v5_source_attested" in text
    assert 'require_json_pass "${COOCCURRENCE_PREFLIGHT}"' in text
    assert 'require_json_pass "${COOCCURRENCE}/release_receipt.json"' in text
    assert '--cooccurrence-release-receipt "${COOCCURRENCE}/release_receipt.json"' in text
    assert "require_cooccurrence_preflight_reconciliation" in text
    assert "raw_identity_conflicting_analysis_payload_policy" in text
    assert "raw_identity_duplicate_audit.tsv.gz" in text
    assert "verify_retrospective_running_source_identity_v2.py" in text
    assert '--independent-m2-audit "${INDEPENDENT_M2_AUDIT}"' in text
    assert '--tumor-ref-source-identity-receipt "${SOURCE_RECEIPT}"' in text
    assert "--allow-partial-scope" not in text
    assert text.count('--selection-column "${FORMAL_SELECTION_COLUMN}"') == 2
    assert text.count("--selection-value true") == 2
    assert text.count(".pass == true and .overlapCount == 0") == 2
    assert "PASS_LOGIC_INDEPENDENT_RECOUNT" in text
    assert "verify_source_authority" in text
    assert "release_source_authority.v7.approval.json" in text
    assert "release_source_authority.v1.approval.json" not in text
    assert "release_source_authority.v2.approval.json" not in text
    assert "release_source_authority.v3.approval.json" not in text
    assert "release_source_authority.v6.approval.json" not in text
    assert "20260719_all_ssnv_result_v5_post_reboot_bootstrap/ed25519_public.pem" in text
    assert "20260719_all_ssnv_report_v5_post_reboot_bootstrap/ed25519_public.pem" in text
    assert "20260717_all_ssnv_result_v1" not in text
    assert "20260717_all_ssnv_report_v1" not in text
    assert "finalize_task_b_result_release.py" in text
    assert "task_b_final_dataset_release_receipt.v1.json" in text
    assert "--final-release-signature \"${FINAL_RELEASE_SIGNATURE}\"" in text
    assert "Awaiting one-time detached result signature" in text
    assert "task_b_final_report_release_receipt.v1.json" in text
    assert "Awaiting one-time detached report signature" in text
    assert "--create-report" in text
    assert "--verify-report" in text
    assert 'chmod 0444 "${DESKTOP_QA}" "${MOBILE_QA}"' in text
    assert 'export PATH="/usr/bin:/bin"' in text
    assert re.search(
        r'"\$\{(?:PYTHON|QA_PYTHON)\}" -I(?! -X "pycache_prefix=\$\{PYTHON_CACHE_ROOT\}")',
        text,
    ) is None
    assert 'require_absent "${path}"' in text
    assert '"${PYTHON_CACHE_ROOT}"' in text
    assert '/usr/bin/mkdir --mode=0700 -- "${PYTHON_CACHE_ROOT}"' in text


def test_m2v5_runner_and_strict_producer_share_python_cache_contract() -> None:
    runner = m2v5_runner_text()
    strict = strict_producer_text()
    runner_match = re.search(
        r'^readonly PYTHON_CACHE_ROOT="\$\{WORKSPACE_ROOT\}/([^"]+)"$',
        runner,
        re.MULTILINE,
    )
    strict_match = re.search(
        r'^CANONICAL_PYTHON_CACHE_DIRNAME = "([^"]+)"$',
        strict,
        re.MULTILINE,
    )
    assert runner_match is not None
    assert strict_match is not None
    assert runner_match.group(1) == strict_match.group(1)
    assert "*canonical_python_prefix(output_dir)" in strict
    assert "observed_process_command() != expected_command" in strict
    for path in (
        MATCHED_NORMAL_RUNNER,
        MATCHED_NORMAL_ANALYZER,
        CN_CCF_ANNOTATOR,
        PRIMARY_AUDITOR,
        FINAL_RESULT_FINALIZER,
        REPORT_BUILDER,
    ):
        source = path.read_text(encoding="utf-8")
        cache_match = re.search(
            r'^\s*"(\.python_cache_m2v5_completion_v2_bound_bootstrap)"$',
            source,
            re.MULTILINE,
        )
        assert cache_match is not None, path
        assert cache_match.group(1) == runner_match.group(1)
        assert "pycache_prefix={CANONICAL_PYTHON_CACHE_ROOT}" in source


def test_m2v5_dataset_and_report_paths_align_across_producers() -> None:
    runner = m2v5_runner_text()
    builder = FINAL_DATASET_BUILDER.read_text(encoding="utf-8")
    finalizer = FINAL_RESULT_FINALIZER.read_text(encoding="utf-8")
    dataset_name = "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
    report_name = "all_ssnv_final_report_v5_m2v5_source_attested"

    assert dataset_name in runner
    assert dataset_name in builder
    assert report_name in runner
    assert report_name in finalizer
    assert "all_ssnv_final_report_dataset_v4_m2v5_source_attested" not in builder
    assert "all_ssnv_final_report_v4_m2v5_source_attested" not in finalizer
    fresh_input_names = (
        "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
        "source_locked_command_parity",
        "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json",
        "strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5",
        "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json",
        "stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json",
        "matched_normal_candidate_controls_v3_m2v5_source_authority_v5",
        "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5",
        "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    )
    for name in fresh_input_names:
        assert name in runner
        assert name in builder
    stale_builder_input_names = (
        "methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked",
        "cooccurrence_task_contract_preflight.v6_analysis_payload_release_full_runtime.json",
        "strict_methyl_candidate_confirmation_v2_m2v5",
        "stable_primary_artifact_audit.v2_source_authorized_pre_downstream.json",
        "stable_primary_artifact_audit.v2_post_m2v5.json",
        "matched_normal_candidate_control_analysis_v2_m2v5",
        "candidate_cn_ccf_annotations_v2_m2v5",
    )
    assert all(name not in builder for name in stale_builder_input_names)
    assert "canonical_task_b_final_builder_commands" in finalizer
    assert "*STRICT_PRODUCER.canonical_python_prefix(bundle.strict_dir)" in builder


def test_current_release_leaf_producers_share_runner_path_versions() -> None:
    runner = m2v5_runner_text()
    cooccurrence_runner = cooccurrence_v6_runner_text()
    release_finalizer = COOCCURRENCE_RELEASE_FINALIZER.read_text(encoding="utf-8")
    matched_runner = MATCHED_NORMAL_RUNNER.read_text(encoding="utf-8")
    matched_analyzer = MATCHED_NORMAL_ANALYZER.read_text(encoding="utf-8")
    cn_ccf = CN_CCF_ANNOTATOR.read_text(encoding="utf-8")
    current_names = {
        "cooccurrence": (
            "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
            "source_locked_command_parity"
        ),
        "preflight": (
            "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
        ),
        "primary_pre": (
            "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
        ),
        "matched_run": (
            "matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
        ),
        "matched_analysis": (
            "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5"
        ),
        "cn_ccf": "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5",
    }
    assert current_names["cooccurrence"] in cooccurrence_runner
    assert current_names["cooccurrence"] in release_finalizer
    assert current_names["preflight"] in cooccurrence_runner
    assert current_names["preflight"] in release_finalizer
    assert current_names["primary_pre"] in cooccurrence_runner
    assert current_names["primary_pre"] in release_finalizer
    assert current_names["cooccurrence"] in matched_runner
    assert current_names["matched_run"] in runner
    assert current_names["matched_run"] in matched_runner
    assert current_names["matched_run"] in matched_analyzer
    assert current_names["matched_analysis"] in runner
    assert current_names["matched_analysis"] in matched_analyzer
    assert current_names["cooccurrence"] in cn_ccf
    assert current_names["cn_ccf"] in runner
    assert current_names["cn_ccf"] in cn_ccf
    for command_fragment in (
        '--candidate-table "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz"',
        '--output-root "${MATCHED_RUN}"',
        '--paired-output-root "${MATCHED_RUN}"',
        '--primary-assignments "${SCREEN_ASSIGNMENTS}"',
        '--output-dir "${MATCHED_ANALYSIS}"',
        '--input "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz"',
        '--output-dir "${CN_CCF}"',
    ):
        assert command_fragment in runner

    stale_names = (
        "methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked",
        "cooccurrence_task_contract_preflight.v6_analysis_payload_release_full_runtime.json",
        "stable_primary_artifact_audit.v2_source_authorized_pre_downstream.json",
        "matched_normal_candidate_controls_v2_m2v5",
        "matched_normal_candidate_control_analysis_v2_m2v5",
        "candidate_cn_ccf_annotations_v2_m2v5",
    )
    current_leaf_sources = (
        release_finalizer,
        matched_runner,
        matched_analyzer,
        cn_ccf,
    )
    assert all(
        stale_name not in source
        for stale_name in stale_names
        for source in current_leaf_sources
    )


def test_cooccurrence_v6_runner_binds_all_dependency_sources_and_sparse_audit() -> None:
    subprocess.run(["bash", "-n", str(COOCCURRENCE_V6_RUNNER)], check=True)
    text = cooccurrence_v6_runner_text()
    assert "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime" in text
    assert (
        "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
        "source_locked_command_parity" in text
    )
    assert "stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream" in text
    for source_name in (
        "preflight",
        "analyzer",
        "ssnv_cooccurrence_lib",
        "latest_tag_join",
        "m2_screen_gate",
    ):
        assert source_name in text
    assert "source_identity_before == .code.source_identity_after" in text
    assert "raw_identity_duplicate_audit.tsv.gz" in text
    assert "finalize_cooccurrence_release_receipt.py" in text
    assert "RELEASE_RECONCILIATION_PASS" in text
    assert "hard_fail_before_site_result" in text
    assert "raw_identity_missing_projections" not in text
    assert "raw_identity_conflicting_analysis_payload_projections" not in text
    assert "verify_source_authority" in text
    assert "source_authority_validator" in text
    assert 'export PATH="/usr/bin:/bin"' in text
    assert re.search(
        r'"\$\{PYTHON\}" -I(?! -X "pycache_prefix=\$\{PYTHON_CACHE_ROOT\}")',
        text,
    ) is None
    assert 'Refusing to reuse Python cache: %s\\n' in text
    assert '/usr/bin/mkdir --mode=0700 -- "${PYTHON_CACHE_ROOT}"' in text


def test_cooccurrence_runner_analyzer_and_finalizer_share_python_cache_contract() -> None:
    runner = cooccurrence_v6_runner_text()
    analyzer = cooccurrence_analyzer_text()
    finalizer = COOCCURRENCE_RELEASE_FINALIZER.read_text(encoding="utf-8")
    cache_name = ".python_cache_cooccurrence_v8_command_parity_bound_bootstrap"
    assert f'PYTHON_CACHE_ROOT="${{WORKSPACE_ROOT}}/{cache_name}"' in runner
    assert "CANONICAL_PYTHON_CACHE_DIRNAME = (" in analyzer
    assert f'"{cache_name}"' in analyzer
    assert "*canonical_python_prefix(output_dir)" in analyzer
    assert "*ANALYZER.canonical_python_prefix(args.output.resolve().parent)" in finalizer
    assert runner.count('-I -X "pycache_prefix=${PYTHON_CACHE_ROOT}"') == 4


def test_source_set_digest_uses_real_newline_delimiters() -> None:
    authority_id = (
        "20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
    )
    public_key_path = (
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
        "20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
    )
    public_key_sha256 = (
        "cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
    )
    for text in (cooccurrence_v6_runner_text(), m2v5_runner_text()):
        assert r'printf "\n"' in text
        assert r'printf "\\n"' not in text
        assert f'readonly SOURCE_AUTHORITY_ID="{authority_id}"' in text
        assert f'readonly SOURCE_AUTHORITY_PUBLIC_KEY="{public_key_path}"' in text
        assert f'readonly SOURCE_AUTHORITY_PUBLIC_KEY_SHA256="{public_key_sha256}"' in text
        assert ".authority_id == $authority_id" in text
        assert ".public_key.path == $public_key" in text
        assert ".public_key.sha256 == $public_key_sha256" in text
        assert "20260719_all_ssnv_focal_alt_task_b_release_v6_command_parity" not in text
        assert "20260719_all_ssnv_v9_post_reboot_bootstrap" not in text
        assert "cf6b8b5eef26e3a53dc6b40d7cd98b6a809257c324322cdbe9dee4e81fdd81b5" not in text


def test_m2v5_runner_preserves_primary_consumer_receipt_order() -> None:
    text = m2v5_runner_text()
    expected = (
        '--consumer-receipt "${COOCCURRENCE}/run_receipt.json"',
        '--consumer-receipt "${COOCCURRENCE}/release_receipt.json"',
        '--consumer-receipt "${TUMOR_REF}/run_manifest.json"',
        '--consumer-receipt "${INDEPENDENT_M2_AUDIT}"',
        '--consumer-receipt "${STRICT_RECEIPT}"',
        '--consumer-receipt "${MATCHED_RUN_RECEIPT}"',
        '--consumer-receipt "${MATCHED_ANALYSIS}/run_receipt.json"',
    )
    positions = [text.index(token) for token in expected]
    assert positions == sorted(positions)


def test_final_result_release_uses_direct_cli_and_detached_ed25519() -> None:
    text = FINAL_RESULT_FINALIZER.read_text(encoding="utf-8")
    assert 'Path("/proc/self/cmdline")' in text
    assert "BoundArtifactReader" in text
    assert "verify_ed25519_signature_fds" in text
    assert "RESULT_PUBLIC_KEY_SHA256" in text
    assert "REPORT_PUBLIC_KEY_SHA256" in text
    assert "REPORT_CREATE_CHECKS" in text
    assert "signed_final_dataset_release_reverified" in text
    assert "all_report_outputs_reverified" in text
    assert "retired one-time result private key" in text
    assert "formal_task_b_release_eligible" in text
    assert "20260719_all_ssnv_result_v5_post_reboot_bootstrap/ed25519_public.pem" in text
    assert "20260719_all_ssnv_report_v5_post_reboot_bootstrap/ed25519_public.pem" in text
    assert "20260717_all_ssnv_result_v1" not in text
    assert "20260717_all_ssnv_report_v1" not in text
    assert '"/dev/stdin"' not in text
    assert "verify_ed25519_signature_fds" in text


def test_source_authority_uses_seekable_payload_fd_for_ed25519() -> None:
    text = SOURCE_AUTHORITY_VALIDATOR.read_text(encoding="utf-8")
    signature_block = text[
        text.index("def verify_ed25519_signature_fds") :
        text.index("def _reject_duplicate_review_keys")
    ]
    git_start = text.index("current_head = subprocess.run(")
    git_block = text[git_start : text.index(").stdout.strip()", git_start)]
    assert 'f"/proc/self/fd/{data_fd}"' in text
    assert 'f"/proc/self/fd/{openssl_fd}"' in text
    assert "env=CLEAN_SUBPROCESS_ENV" in signature_block
    assert "GIT_CONFIG_NOSYSTEM" not in signature_block
    assert "GIT_CONFIG_GLOBAL" not in signature_block
    assert "retain_until_process_exit" in text
    assert '"GIT_CONFIG_NOSYSTEM": "1"' in git_block
    assert '"GIT_CONFIG_GLOBAL": "/dev/null"' in git_block
    assert "env=CLEAN_SUBPROCESS_ENV" not in git_block
    assert '"/dev/stdin"' not in text
    assert (
        "pass_fds=(openssl_fd, data_fd, public_key_fd, signature_fd)" in text
    )


def test_shell_release_runners_clear_loader_and_openssl_injection() -> None:
    for path in (COOCCURRENCE_V6_RUNNER, M2V5_RUNNER):
        text = path.read_text(encoding="utf-8")
        assert "LD_LIBRARY_PATH LD_PRELOAD" in text
        assert "OPENSSL_CONF OPENSSL_MODULES" in text
        assert "/usr/bin/env -i PATH=/usr/bin:/bin LC_ALL=C LANG=C" in text


def test_v2_one_time_signer_pins_system_openssl_and_retires_only_after_verify() -> None:
    subprocess.run(["bash", "-n", str(ONE_TIME_SIGNER)], check=True)
    text = ONE_TIME_SIGNER.read_text(encoding="utf-8")
    assert 'readonly OPENSSL="/usr/bin/openssl"' in text
    assert '"${OPENSSL}" pkeyutl -sign -rawin' in text
    assert '"${OPENSSL}" pkeyutl -verify -rawin' in text
    assert text.index('"${OPENSSL}" pkeyutl -verify -rawin') < text.index(
        '"${CHMOD}" 000 "${PRIVATE_KEY}"'
    )
    assert "SIGN_OPERATION_FAILED_KEY_RETAINED" in text
    assert "SIGN_VERIFY_FAILED_KEY_RETAINED" in text
