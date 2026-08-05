#!/usr/bin/env bash

set -Eeuo pipefail

export PATH="/usr/bin:/bin"
unset BASH_ENV ENV CDPATH GLOBIGNORE LD_AUDIT LD_DEBUG LD_LIBRARY_PATH LD_PRELOAD
unset OPENSSL_CONF OPENSSL_MODULES PYTHONHOME PYTHONINSPECT PYTHONPATH PYTHONSTARTUP
for inherited_function in jq stat sha256sum realpath sort awk git openssl date tee; do
    unset -f "${inherited_function}" 2>/dev/null || true
done

readonly REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
readonly TOPIC_ROOT="${REPO_ROOT}/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly SCRIPT_ROOT="${TOPIC_ROOT}/scripts"
readonly RESULT_ROOT="${TOPIC_ROOT}/results"
readonly WORKSPACE_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
readonly QA_PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/bin/python"
readonly NODE="/bip7_disk/liaoyoyo2001/.nvm/versions/node/v22.22.1/bin/node"
readonly NODE_SHA256="243fd8938011479f41b3de101842150fa990f33fbbb3f7aabd330857f2d79e1d"
readonly QA_CHROMIUM="/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium_headless_shell-1217/chrome-headless-shell-linux64/chrome-headless-shell"
readonly QA_CHROMIUM_SHA256="22d72b589f908111fa29ff4e61ead8a0601e925177f59c70e2e851c22e9b1886"
readonly INTERSUBMOD_BIN="${REPO_ROOT}/build/bin/inter_sub_mod"
readonly REFERENCE="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
readonly PLUGIN_ROOT="/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599"
readonly SOURCE_AUTHORITY="${REPO_ROOT}/docs/provenance/source_authorities/20260722_all_ssnv_focal_alt_release_source_authority.v7.json"
readonly SOURCE_AUTHORITY_APPROVAL="${REPO_ROOT}/docs/provenance/source_authorities/20260722_all_ssnv_focal_alt_release_source_authority.v7.approval.json"
readonly SOURCE_AUTHORITY_SIGNATURE="${SOURCE_AUTHORITY_APPROVAL}.ed25519.sig"
readonly SOURCE_AUTHORITY_ID="20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
readonly SOURCE_AUTHORITY_PUBLIC_KEY="/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
readonly SOURCE_AUTHORITY_PRIVATE_KEY="/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_private_encrypted_unrecoverable_after_signing.pem"
readonly SOURCE_AUTHORITY_PUBLIC_KEY_SHA256="cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
readonly OPENSSL="/usr/bin/openssl"
readonly OPENSSL_SHA256="38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"

readonly MANIFEST="${RESULT_ROOT}/all_ssnv_input_manifest.json"
readonly SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
readonly SCREEN_SITES="${SCREEN}/all_ssnv_site_results.tsv.gz"
readonly SCREEN_ASSIGNMENTS="${SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
readonly SCREEN_SUMMARY="${SCREEN}/all_ssnv_summary.json"
readonly EXTRACT_ROOT="${WORKSPACE_ROOT}/intersubmod_all_ssnv_v2_verification_fix"
readonly PRIMARY_PRE="${RESULT_ROOT}/stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"

readonly COOCCURRENCE="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
readonly TUMOR_REF="${WORKSPACE_ROOT}/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
readonly STRICT="${WORKSPACE_ROOT}/strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5"
readonly MATCHED_RUN="${WORKSPACE_ROOT}/matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
readonly MATCHED_ANALYSIS="${WORKSPACE_ROOT}/matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5"
readonly CN_CCF="${WORKSPACE_ROOT}/candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"

readonly SOURCE_SNAPSHOT="${WORKSPACE_ROOT}/tumor_ref_recovery_source_identity_v1/observed_during_execution.snapshot.json"
readonly SOURCE_RECEIPT="${WORKSPACE_ROOT}/tumor_ref_recovery_source_identity_v1/post_run_source_identity.receipt.json"
readonly INDEPENDENT_M2_AUDIT="${RESULT_ROOT}/independent_m2_gate_recount.v3.json"
readonly COOCCURRENCE_PREFLIGHT="${RESULT_ROOT}/cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
readonly PRIMARY_POST="${RESULT_ROOT}/stable_primary_artifact_audit.v6_strict_command_parity_post_downstream.json"
readonly FROZEN_POST="${RESULT_ROOT}/frozen_input_immutability.post_m2v5_downstream_v3_source_authority_v5.json"

readonly FINAL_DATASET="${WORKSPACE_ROOT}/all_ssnv_final_report_dataset_v5_m2v5_source_attested"
readonly FINAL_RELEASE_FINALIZER="${SCRIPT_ROOT}/finalize_task_b_result_release.py"
readonly FINAL_RELEASE_RECEIPT="${RESULT_ROOT}/task_b_final_dataset_release_receipt.v1.json"
readonly FINAL_RELEASE_SIGNATURE="${FINAL_RELEASE_RECEIPT}.ed25519.sig"
readonly FINAL_RELEASE_PUBLIC_KEY="/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/20260719_all_ssnv_result_v5_post_reboot_bootstrap/ed25519_public.pem"
readonly FINAL_RELEASE_PUBLIC_KEY_SHA256="5b7d5d026835ec6ec677bcd886bc16ac71444117dabc44a084ce3ede3a4db5a9"
readonly REPORT="${WORKSPACE_ROOT}/all_ssnv_final_report_v5_m2v5_source_attested"
readonly REPORT_RELEASE_RECEIPT="${RESULT_ROOT}/task_b_final_report_release_receipt.v1.json"
readonly REPORT_RELEASE_SIGNATURE="${REPORT_RELEASE_RECEIPT}.ed25519.sig"
readonly REPORT_RELEASE_PUBLIC_KEY="/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/20260719_all_ssnv_report_v5_post_reboot_bootstrap/ed25519_public.pem"
readonly REPORT_RELEASE_PRIVATE_KEY="/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/20260719_all_ssnv_report_v5_post_reboot_bootstrap/ed25519_private_encrypted_unrecoverable_after_signing.pem"
readonly REPORT_RELEASE_PUBLIC_KEY_SHA256="572e27167e1eea4c39ca53546ba33868b2874ec7fa3ed1682db821b2e50fa439"
readonly CLAIM_CONTRACT="${TOPIC_ROOT}/claim-contract-v5.md"
readonly EARLIER_FP_REPORT="${REPO_ROOT}/research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.md"
readonly LOG_PATH="${WORKSPACE_ROOT}/m2v5_recovered_completion_chain_v2_source_authority_v5.log"
readonly PYTHON_CACHE_ROOT="${WORKSPACE_ROOT}/.python_cache_m2v5_completion_v2_bound_bootstrap"

readonly FORMAL_SELECTION_COLUMN="multi_marker_molecular_haplotype_base_candidate"
readonly PORTABLE_HTML="${REPORT}/all_ssnv_focal_alt_multigroup_cooccurrence_report.standalone.html"
readonly PORTABLE_RECEIPT="${REPORT}/portable_report_delivery_receipt.json"
readonly PORTABLE_VERIFY_SCREENSHOT="${REPORT}/portable_report_official_verification.png"
readonly DESKTOP_SCREENSHOT="${REPORT}/portable_report_desktop_1440x1000.png"
readonly MOBILE_SCREENSHOT="${REPORT}/portable_report_mobile_390x844.png"
readonly DESKTOP_QA="${REPORT}/portable_report_desktop_qa.json"
readonly MOBILE_QA="${REPORT}/portable_report_mobile_qa.json"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
unset PYTHONPATH

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf 'Required file is missing or empty: %s\n' "${path}" >&2
        exit 1
    }
}

verify_source_authority() {
    [[ "$(stat -c '%a' "${SOURCE_AUTHORITY}")" == "444" ]]
    [[ "$(stat -c '%a' "${SOURCE_AUTHORITY_APPROVAL}")" == "444" ]]
    [[ "$(stat -c '%a' "${SOURCE_AUTHORITY_SIGNATURE}")" == "444" ]]
    [[ "$(stat -c '%a' "${SOURCE_AUTHORITY_PUBLIC_KEY}")" == "444" ]]
    [[ "$(stat -c '%a' "${SOURCE_AUTHORITY_PRIVATE_KEY}")" == "0" ]]
    [[ "$(sha256sum "${SOURCE_AUTHORITY_PUBLIC_KEY}" | awk '{print $1}')" == "${SOURCE_AUTHORITY_PUBLIC_KEY_SHA256}" ]]
    [[ "$(sha256sum "${OPENSSL}" | awk '{print $1}')" == "${OPENSSL_SHA256}" ]]
    /usr/bin/env -i PATH=/usr/bin:/bin LC_ALL=C LANG=C \
        "${OPENSSL}" pkeyutl -verify -rawin -pubin \
        -inkey "${SOURCE_AUTHORITY_PUBLIC_KEY}" \
        -sigfile "${SOURCE_AUTHORITY_SIGNATURE}" \
        -in "${SOURCE_AUTHORITY_APPROVAL}" >/dev/null
    jq -e \
        --arg authority_id "${SOURCE_AUTHORITY_ID}" \
        '
        .schema_name == "intersubmod.release_source_authority" and
        .schema_version == "1.0.0" and
        .authority_id == $authority_id and
        .approval_status == "APPROVED_FOR_FULL_TASK_B_RUN" and
        (.sources | length) == 23 and
        ([.sources[].mode] | all(. == "0o444"))
    ' "${SOURCE_AUTHORITY}" >/dev/null
    jq -e \
        --arg authority_id "${SOURCE_AUTHORITY_ID}" \
        --arg authority "$(realpath "${SOURCE_AUTHORITY}")" \
        --arg size "$(stat -c '%s' "${SOURCE_AUTHORITY}")" \
        --arg sha "$(sha256sum "${SOURCE_AUTHORITY}" | awk '{print $1}')" \
        --arg public_key "$(realpath "${SOURCE_AUTHORITY_PUBLIC_KEY}")" \
        --arg public_key_sha256 "${SOURCE_AUTHORITY_PUBLIC_KEY_SHA256}" \
        '
        .schema_name == "intersubmod.release_source_authority.approval" and
        .schema_version == "1.0.0" and
        .authority_id == $authority_id and
        .approval_status == "APPROVED_FOR_FULL_TASK_B_RUN" and
        (.review_approvals | length) >= 2 and
        ([.review_approvals[].verdict] | all(. == "APPROVE")) and
        .authority_manifest.path == $authority and
        (.authority_manifest.size_bytes | tostring) == $size and
        .authority_manifest.sha256 == $sha and
        .authority_manifest.mode == "0o444" and
        .public_key.path == $public_key and
        .public_key.sha256 == $public_key_sha256 and
        .public_key.mode == "0o444"
        ' "${SOURCE_AUTHORITY_APPROVAL}" >/dev/null
    [[ "$(jq -er '.repository.git_head_at_authorization' "${SOURCE_AUTHORITY}")" == "$(/usr/bin/git rev-parse HEAD)" ]]
    while IFS=$'\t' read -r role path size sha mode; do
        require_file "${path}"
        [[ "$(stat -c '%s' "${path}")" == "${size}" ]]
        [[ "$(sha256sum "${path}" | awk '{print $1}')" == "${sha}" ]]
        [[ "$(stat -c '%a' "${path}")" == "444" ]]
        [[ "${mode}" == "0o444" ]]
    done < <(jq -r '.sources | to_entries[] | [.key, .value.path, .value.size_bytes, .value.sha256, .value.mode] | @tsv' "${SOURCE_AUTHORITY}")
    observed_source_set_sha256="$({
        jq -r '.sources | to_entries[] | [.key, .value.path, (.value.size_bytes|tostring), .value.sha256, .value.mode] | join("|")' "${SOURCE_AUTHORITY}"
    } | sort | awk 'NR > 1 { printf "\n" } { printf "%s", $0 }' | sha256sum | awk '{print $1}')"
    [[ "$(jq -er '.source_set_sha256' "${SOURCE_AUTHORITY}")" == "${observed_source_set_sha256}" ]]
}

require_absent() {
    local path="$1"
    [[ ! -e "${path}" ]] || {
        printf 'Refusing to overwrite existing path: %s\n' "${path}" >&2
        exit 1
    }
}

require_json_pass() {
    local path="$1"
    require_file "${path}"
    jq -e '.pass == true' "${path}" >/dev/null || {
        printf 'JSON receipt does not declare pass=true: %s\n' "${path}" >&2
        exit 1
    }
}

verify_cooccurrence_declared_output() {
    local key="$1"
    local expected_path="$2"
    local declared_path declared_size declared_sha256
    declared_path="$(jq -er --arg key "${key}" '.outputs[$key].path' "${COOCCURRENCE}/run_receipt.json")"
    declared_size="$(jq -er --arg key "${key}" '.outputs[$key].size_bytes' "${COOCCURRENCE}/run_receipt.json")"
    declared_sha256="$(jq -er --arg key "${key}" '.outputs[$key].sha256' "${COOCCURRENCE}/run_receipt.json")"
    [[ "${declared_path}" == "$(realpath "${expected_path}")" ]]
    [[ "${declared_size}" == "$(stat -c '%s' "${expected_path}")" ]]
    [[ "${declared_sha256}" == "$(sha256sum "${expected_path}" | awk '{print $1}')" ]]
}

verify_cooccurrence_preflight_input() {
    local receipt="${COOCCURRENCE}/run_receipt.json"
    [[ "$(jq -er '.inputs.preflight_receipt.path' "${receipt}")" == "$(realpath "${COOCCURRENCE_PREFLIGHT}")" ]]
    [[ "$(jq -er '.inputs.preflight_receipt.size_bytes' "${receipt}")" == "$(stat -c '%s' "${COOCCURRENCE_PREFLIGHT}")" ]]
    [[ "$(jq -er '.inputs.preflight_receipt.sha256' "${receipt}")" == "$(sha256sum "${COOCCURRENCE_PREFLIGHT}" | awk '{print $1}')" ]]
}

require_cooccurrence_preflight_reconciliation() {
    jq -e \
        --slurpfile receipt "${COOCCURRENCE}/run_receipt.json" \
        --slurpfile summary "${COOCCURRENCE}/summary.json" \
        '
        .code.source_identity_after.analyzer == $receipt[0].code.analyzer and
        .code.source_identity_after.ssnv_cooccurrence_lib == $receipt[0].code.ssnv_cooccurrence_lib and
        .code.source_identity_after.latest_tag_join == $receipt[0].code.latest_tag_join and
        .code.source_identity_after.m2_screen_gate == $receipt[0].code.m2_screen_gate and
        .code.source_identity_after.release_finalizer == $receipt[0].code.release_finalizer and
        .code.source_identity_after.release_runner == $receipt[0].code.release_runner and
        .code.source_identity_after.source_authority_validator == $receipt[0].code.source_authority_validator and
        .source_authority == $receipt[0].source_authority and
        .schema_name == "intersubmod.cooccurrence_task_contract_preflight" and
        .schema_version == "2.0.0" and
        .task_type == "B_comprehensive_validation" and
        .observed.task_count == 102842 and
        (.code.source_identity_after | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "preflight", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
        (.code.source_identity_before | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "preflight", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
        (.code.source_modes_before | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "preflight", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
        (.code.source_modes_after | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "preflight", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
        .raw_identity_release_contract == "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3" and
        $summary[0].raw_identity_release_contract == "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3" and
        $receipt[0].raw_identity_release_contract == "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3" and
        $receipt[0].task_type == "B_comprehensive_validation" and
        $receipt[0].scope == "all_manifest_samples" and
        $receipt[0].full_scope == true and
        ($receipt[0].code | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
        $receipt[0].preflight_gate.schema_name == "intersubmod.cooccurrence_task_contract_preflight" and
        $receipt[0].preflight_gate.schema_version == "2.0.0" and
        $receipt[0].preflight_gate.task_count == 102842 and
        $receipt[0].preflight_gate.pass == true and
        $receipt[0].preflight_gate.raw_identity_release_contract == "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3" and
        .observed.raw_bam_identity_recovery.equivalence_contract ==
            $summary[0].raw_bam_identity_recovery_audit.equivalence_contract and
        .observed.raw_bam_identity_recovery.equivalence_contract ==
            "sam_core_and_all_aux_tags_except_RG_exact_v1" and
        .observed.raw_bam_identity_recovery.allowed_differing_auxiliary_tags == ["RG"] and
        $summary[0].raw_bam_identity_recovery_audit.allowed_differing_auxiliary_tags == ["RG"] and
        .observed.raw_bam_identity_recovery.site_weighted_audit_sha256 ==
            $summary[0].raw_bam_identity_recovery_audit.site_weighted_audit_sha256 and
        .observed.raw_bam_identity_recovery.aggregate.site_tasks ==
            $summary[0].raw_bam_identity_recovery_audit.n_site_tasks and
        .observed.raw_bam_identity_recovery.aggregate.expected_projection_occurrences ==
            $summary[0].raw_bam_identity_recovery_audit.n_expected_projection_occurrences and
        .observed.raw_bam_identity_recovery.aggregate.matched_record_occurrences ==
            $summary[0].raw_bam_identity_recovery_audit.n_matched_record_occurrences and
        .observed.raw_bam_identity_recovery.aggregate.sites_with_collapsed_duplicates ==
            $summary[0].raw_bam_identity_recovery_audit.n_sites_with_collapsed_duplicates and
        .observed.raw_bam_identity_recovery.aggregate.duplicate_projection_occurrences_collapsed ==
            $summary[0].raw_bam_identity_recovery_audit.n_duplicate_projection_occurrences_collapsed and
        .observed.raw_bam_identity_recovery.aggregate.duplicate_extra_record_occurrences_collapsed ==
            $summary[0].raw_bam_identity_recovery_audit.n_duplicate_extra_record_occurrences_collapsed and
        .observed.raw_bam_identity_recovery.aggregate.exact_duplicate_projection_occurrences_collapsed ==
            $summary[0].raw_bam_identity_recovery_audit.n_exact_duplicate_projection_occurrences_collapsed and
        .observed.raw_bam_identity_recovery.aggregate.rg_only_duplicate_projection_occurrences_collapsed ==
            $summary[0].raw_bam_identity_recovery_audit.n_rg_only_duplicate_projection_occurrences_collapsed and
        .observed.raw_bam_identity_recovery.aggregate.duplicate_projection_occurrences_collapsed ==
            $receipt[0].counts.raw_identity_sparse_duplicate_rows and
        .observed.raw_bam_identity_recovery.all_result_rows_passed_invariant_validation == true and
        $receipt[0].counts.raw_identity_all_site_results_passed_invariant_validation == true and
        .observed.raw_bam_identity_recovery.missing_projection_policy ==
            $receipt[0].counts.raw_identity_missing_projection_policy and
        .observed.raw_bam_identity_recovery.missing_projection_policy ==
            "hard_fail_before_site_result" and
        .observed.raw_bam_identity_recovery.conflicting_analysis_payload_policy ==
            $receipt[0].counts.raw_identity_conflicting_analysis_payload_policy and
        .observed.raw_bam_identity_recovery.conflicting_analysis_payload_policy ==
            "hard_fail_before_site_result" and
        .observed.raw_bam_identity_recovery.failure_counts_materialized == false and
        $receipt[0].counts.raw_identity_failure_counts_materialized == false and
        $receipt[0].source_lock.all_sources_read_only_and_unchanged == true and
        $receipt[0].source_lock.source_identity_before == $receipt[0].code and
        $receipt[0].source_lock.source_identity_after_compute == $receipt[0].code and
        $receipt[0].source_lock.source_identity_after_output == $receipt[0].code and
        ([$receipt[0].source_lock.source_modes_before[], $receipt[0].source_lock.source_modes_after_compute[], $receipt[0].source_lock.source_modes_after_output[]] | all(. == "0o444"))
        ' "${COOCCURRENCE_PREFLIGHT}" >/dev/null || {
        printf 'Cooccurrence output does not reconcile with the full raw-identity preflight\n' >&2
        exit 1
    }
}

run_step() {
    local label="$1"
    shift
    printf '\n[%s] %s\n' "$(date --iso-8601=seconds)" "${label}"
    printf 'Command:'
    printf ' %q' "$@"
    printf '\n'
    verify_source_authority
    "$@"
    verify_source_authority
}

on_error() {
    local exit_code=$?
    printf '\nFAILED exit=%d line=%d command=%q\n' \
        "${exit_code}" "${BASH_LINENO[0]}" "${BASH_COMMAND}" >&2
    exit "${exit_code}"
}
trap on_error ERR

cd "${REPO_ROOT}"
umask 0022

for path in \
    "${MANIFEST}" \
    "${SCREEN_SITES}" \
    "${SCREEN_ASSIGNMENTS}" \
    "${SCREEN_SUMMARY}" \
    "${PRIMARY_PRE}" \
    "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    "${COOCCURRENCE}/methyl_ssnv_pair_results.tsv.gz" \
    "${COOCCURRENCE}/raw_identity_duplicate_audit.tsv.gz" \
    "${COOCCURRENCE}/oracle_cases.json" \
    "${COOCCURRENCE}/summary.json" \
    "${COOCCURRENCE}/run_receipt.json" \
    "${COOCCURRENCE}/release_receipt.json" \
    "${TUMOR_REF}/all_ssnv_tumor_ref_control_site_results.tsv.gz" \
    "${TUMOR_REF}/summary.json" \
    "${TUMOR_REF}/run_manifest.json" \
    "${SOURCE_SNAPSHOT}" \
    "${INDEPENDENT_M2_AUDIT}" \
    "${COOCCURRENCE_PREFLIGHT}" \
    "${CLAIM_CONTRACT}" \
    "${SOURCE_AUTHORITY}" \
    "${SOURCE_AUTHORITY_APPROVAL}" \
    "${SOURCE_AUTHORITY_SIGNATURE}" \
    "${SOURCE_AUTHORITY_PUBLIC_KEY}" \
    "${SOURCE_AUTHORITY_PRIVATE_KEY}" \
    "${FINAL_RELEASE_FINALIZER}" \
    "${FINAL_RELEASE_PUBLIC_KEY}" \
    "${REPORT_RELEASE_PUBLIC_KEY}" \
    "${REPORT_RELEASE_PRIVATE_KEY}" \
    "${EARLIER_FP_REPORT}" \
    "${INTERSUBMOD_BIN}" \
    "${REFERENCE}" \
    "${REFERENCE}.fai"; do
    require_file "${path}"
done
verify_source_authority
[[ "$(sha256sum "${NODE}" | awk '{print $1}')" == "${NODE_SHA256}" ]]
[[ "$(sha256sum "${QA_CHROMIUM}" | awk '{print $1}')" == "${QA_CHROMIUM_SHA256}" ]]
[[ "$(sha256sum "${FINAL_RELEASE_PUBLIC_KEY}" | awk '{print $1}')" == "${FINAL_RELEASE_PUBLIC_KEY_SHA256}" ]]
[[ "$(stat -c '%a' "${REPORT_RELEASE_PUBLIC_KEY}")" == "444" ]]
[[ "$(stat -c '%a' "${REPORT_RELEASE_PRIVATE_KEY}")" == "400" ]]
[[ "$(sha256sum "${REPORT_RELEASE_PUBLIC_KEY}" | awk '{print $1}')" == "${REPORT_RELEASE_PUBLIC_KEY_SHA256}" ]]
require_json_pass "${PRIMARY_PRE}"
require_json_pass "${COOCCURRENCE}/run_receipt.json"
require_json_pass "${COOCCURRENCE}/release_receipt.json"
jq -e '.release_status == "RELEASE_RECONCILIATION_PASS"' \
    "${COOCCURRENCE}/release_receipt.json" >/dev/null
require_json_pass "${TUMOR_REF}/run_manifest.json"
require_json_pass "${SOURCE_SNAPSHOT}"
require_json_pass "${INDEPENDENT_M2_AUDIT}"
require_json_pass "${COOCCURRENCE_PREFLIGHT}"
require_cooccurrence_preflight_reconciliation
verify_cooccurrence_preflight_input
verify_cooccurrence_declared_output site_tsv "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz"
verify_cooccurrence_declared_output pair_tsv "${COOCCURRENCE}/methyl_ssnv_pair_results.tsv.gz"
verify_cooccurrence_declared_output summary_json "${COOCCURRENCE}/summary.json"
verify_cooccurrence_declared_output raw_identity_duplicate_audit_tsv "${COOCCURRENCE}/raw_identity_duplicate_audit.tsv.gz"
verify_cooccurrence_declared_output case_json "${COOCCURRENCE}/oracle_cases.json"

for path in \
    "${STRICT}" \
    "${MATCHED_RUN}" \
    "${MATCHED_ANALYSIS}" \
    "${CN_CCF}" \
    "${SOURCE_RECEIPT}" \
    "${PRIMARY_POST}" \
    "${FROZEN_POST}" \
    "${FINAL_DATASET}" \
    "${REPORT}" \
    "${FINAL_RELEASE_RECEIPT}" \
    "${FINAL_RELEASE_SIGNATURE}" \
    "${REPORT_RELEASE_RECEIPT}" \
    "${REPORT_RELEASE_SIGNATURE}" \
    "${LOG_PATH}" \
    "${PYTHON_CACHE_ROOT}"; do
    require_absent "${path}"
done

umask 0077
/usr/bin/mkdir --mode=0700 -- "${PYTHON_CACHE_ROOT}"
[[ "$(stat -c '%a' "${PYTHON_CACHE_ROOT}")" == "700" ]]
umask 0022

exec > >(tee "${LOG_PATH}") 2>&1

printf 'Input manifest: %s\n' "${MANIFEST}"
printf 'Input screen: %s\n' "${SCREEN}"
printf 'Input cooccurrence: %s\n' "${COOCCURRENCE}"
printf 'Input recovered tumor-REF: %s\n' "${TUMOR_REF}"
printf 'Input independent M2 audit: %s\n' "${INDEPENDENT_M2_AUDIT}"
printf 'Input cooccurrence runtime preflight: %s\n' "${COOCCURRENCE_PREFLIGHT}"
printf 'Output final dataset: %s\n' "${FINAL_DATASET}"
printf 'Output final report: %s\n' "${REPORT}"

run_step "Verify bounded tumor-REF recovery source identity" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/verify_retrospective_running_source_identity_v2.py" verify \
    --snapshot "${SOURCE_SNAPSHOT}" \
    --run-manifest "${TUMOR_REF}/run_manifest.json" \
    --output "${SOURCE_RECEIPT}"
require_json_pass "${SOURCE_RECEIPT}"

run_step "Run strict multi-seed and column-null robustness" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/run_strict_methyl_candidate_confirmation.py" \
    --candidate-table "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    --assignments "${SCREEN_ASSIGNMENTS}" \
    --cooccurrence-receipt "${COOCCURRENCE}/run_receipt.json" \
    --cooccurrence-release-receipt "${COOCCURRENCE}/release_receipt.json" \
    --output-dir "${STRICT}" \
    --permutations 999 \
    --seeds 10 \
    --seed 20260715 \
    --assignment-ari-threshold 0.8 \
    --min-valid-null-fraction 0.8

if [[ -s "${STRICT}/run_manifest.json" ]]; then
    STRICT_RECEIPT="${STRICT}/run_manifest.json"
else
    STRICT_RECEIPT="${STRICT}/not_applicable_receipt.json"
fi
readonly STRICT_RECEIPT
require_json_pass "${STRICT_RECEIPT}"

run_step "Run matched-normal controls for formal G2 selection" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/run_matched_normal_candidate_controls.py" \
    --candidate-table "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    --selection-column "${FORMAL_SELECTION_COLUMN}" \
    --selection-value true \
    --manifest "${MANIFEST}" \
    --normal-audit "${RESULT_ROOT}/matched_normal_methyl_tag_audit.v3_pre_candidate.json" \
    --binary "${INTERSUBMOD_BIN}" \
    --reference "${REFERENCE}" \
    --output-root "${MATCHED_RUN}" \
    --threads-per-sample 40

if [[ -s "${MATCHED_RUN}/run_receipt.json" ]]; then
    MATCHED_RUN_RECEIPT="${MATCHED_RUN}/run_receipt.json"
    require_json_pass "${MATCHED_RUN_RECEIPT}"
    run_step "Analyze matched-normal exact read-identity controls" \
        "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
        "${SCRIPT_ROOT}/analyze_matched_normal_candidate_controls.py" \
        --paired-output-root "${MATCHED_RUN}" \
        --primary-assignments "${SCREEN_ASSIGNMENTS}" \
        --output-dir "${MATCHED_ANALYSIS}"
    require_json_pass "${MATCHED_ANALYSIS}/run_receipt.json"
    MATCHED_FINAL="${MATCHED_ANALYSIS}"
else
    MATCHED_RUN_RECEIPT="${MATCHED_RUN}/not_applicable_receipt.json"
    require_json_pass "${MATCHED_RUN_RECEIPT}"
    MATCHED_FINAL="${MATCHED_RUN}"
fi
readonly MATCHED_RUN_RECEIPT MATCHED_FINAL

run_step "Attach authority-locked CN and fit-local CCF annotations" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/annotate_candidate_cn_ccf.py" \
    --input "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    --selection-column "${FORMAL_SELECTION_COLUMN}" \
    --selection-value true \
    --output-dir "${CN_CCF}"
require_json_pass "${CN_CCF}/receipt.json"

POST_AUDIT_COMMAND=(
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/audit_stable_primary_artifacts.py"
    --site-results "${SCREEN_SITES}"
    --stable-assignments "${SCREEN_ASSIGNMENTS}"
    --consumer-receipt "${COOCCURRENCE_PREFLIGHT}"
    --consumer-receipt "${COOCCURRENCE}/run_receipt.json"
    --consumer-receipt "${COOCCURRENCE}/release_receipt.json"
    --consumer-receipt "${TUMOR_REF}/run_manifest.json"
    --consumer-receipt "${INDEPENDENT_M2_AUDIT}"
    --consumer-receipt "${STRICT_RECEIPT}"
    --consumer-receipt "${MATCHED_RUN_RECEIPT}"
)
if [[ "${MATCHED_FINAL}" == "${MATCHED_ANALYSIS}" ]]; then
    POST_AUDIT_COMMAND+=(--consumer-receipt "${MATCHED_ANALYSIS}/run_receipt.json")
fi
POST_AUDIT_COMMAND+=(
    --output "${PRIMARY_POST}"
    --workers 40
    --chunk-size 64
    --max-pending-chunks 80
    --progress-every 1000
)
run_step "Re-audit primary artifacts and consumer chronology" "${POST_AUDIT_COMMAND[@]}"
require_json_pass "${PRIMARY_POST}"

run_step "Re-audit frozen inputs after downstream analyses" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/audit_frozen_input_immutability.py" \
    --manifest "${MANIFEST}" \
    --output "${FROZEN_POST}"
require_json_pass "${FROZEN_POST}"

run_step "Build source-attested final all-sSNV dataset" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/build_all_ssnv_final_report_dataset.py" \
    --manifest "${MANIFEST}" \
    --screen-dir "${SCREEN}" \
    --cooccurrence-dir "${COOCCURRENCE}" \
    --strict-dir "${STRICT}" \
    --tumor-ref-dir "${TUMOR_REF}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --primary-artifact-audit-post "${PRIMARY_POST}" \
    --cooccurrence-preflight "${COOCCURRENCE_PREFLIGHT}" \
    --tumor-ref-source-identity-receipt "${SOURCE_RECEIPT}" \
    --independent-m2-audit "${INDEPENDENT_M2_AUDIT}" \
    --matched-normal-dir "${MATCHED_FINAL}" \
    --cn-ccf-annotations "${CN_CCF}" \
    --output-dir "${FINAL_DATASET}"
require_json_pass "${FINAL_DATASET}/run_receipt.json"
jq -e \
    '.tumor_ref_source_identity_attestation.release_gate_pass == true and .m2_evaluability_contract.independent_logic_audit.status == "PASS_LOGIC_INDEPENDENT_RECOUNT"' \
    "${FINAL_DATASET}/final_report_dataset.json" >/dev/null

run_step "Create terminal Task-B final dataset release receipt" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${FINAL_RELEASE_FINALIZER}" \
    --create \
    --final-dataset-dir "${FINAL_DATASET}" \
    --builder-receipt "${FINAL_DATASET}/run_receipt.json" \
    --output "${FINAL_RELEASE_RECEIPT}"
require_json_pass "${FINAL_RELEASE_RECEIPT}"

printf '\n[%s] Awaiting one-time detached result signature: %s\n' \
    "$(date --iso-8601=seconds)" "${FINAL_RELEASE_SIGNATURE}"
for ((signature_wait=0; signature_wait<17280; signature_wait++)); do
    [[ -s "${FINAL_RELEASE_SIGNATURE}" ]] && break
    sleep 5
done
require_file "${FINAL_RELEASE_SIGNATURE}"
[[ "$(stat -c '%a' "${FINAL_RELEASE_SIGNATURE}")" == "444" ]]

run_step "Verify terminal Task-B detached result signature" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${FINAL_RELEASE_FINALIZER}" \
    --verify \
    --receipt "${FINAL_RELEASE_RECEIPT}" \
    --signature "${FINAL_RELEASE_SIGNATURE}"

run_step "Build complete Traditional Chinese report and canonical artifact" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/build_all_ssnv_report_artifact.py" \
    --final-dataset "${FINAL_DATASET}/final_report_dataset.json" \
    --final-receipt "${FINAL_DATASET}/run_receipt.json" \
    --final-release-receipt "${FINAL_RELEASE_RECEIPT}" \
    --final-release-signature "${FINAL_RELEASE_SIGNATURE}" \
    --final-release-public-key "${FINAL_RELEASE_PUBLIC_KEY}" \
    --candidate-catalog "${FINAL_DATASET}/candidate_catalog.tsv" \
    --candidate-witness-pairs "${FINAL_DATASET}/candidate_witness_pairs.tsv" \
    --claim-contract "${CLAIM_CONTRACT}" \
    --manifest "${MANIFEST}" \
    --screen-summary "${SCREEN_SUMMARY}" \
    --cooccurrence-sites "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    --cooccurrence-pairs "${COOCCURRENCE}/methyl_ssnv_pair_results.tsv.gz" \
    --cooccurrence-summary "${COOCCURRENCE}/summary.json" \
    --output-reconciliation "${RESULT_ROOT}/all_ssnv_output_reconciliation.v2_verification_fix.json" \
    --post-immutability-audit "${FROZEN_POST}" \
    --tree-input-audit "${RESULT_ROOT}/latest_tree_input_contract_audit.json" \
    --reference-identity-audit "${RESULT_ROOT}/extraction_reference_identity_audit.v1.json" \
    --earlier-fp-report "${EARLIER_FP_REPORT}" \
    --repo-root "${REPO_ROOT}" \
    --output-dir "${REPORT}"
require_json_pass "${REPORT}/report_build_receipt.json"

run_step "Package and officially verify portable HTML" \
    "${NODE}" "${SCRIPT_ROOT}/deliver_portable_report_scrollbar_safe.mjs" \
    --input "${REPORT}/artifact.json" \
    --output "${PORTABLE_HTML}" \
    --receipt "${PORTABLE_RECEIPT}" \
    --screenshot "${PORTABLE_VERIFY_SCREENSHOT}" \
    --plugin-root "${PLUGIN_ROOT}" \
    --ready-timeout-ms 120000 \
    --action-timeout-ms 60000 \
    --timeout-ms 300000
require_json_pass "${PORTABLE_RECEIPT}"

printf '\n[%s] Desktop portable HTML layout QA\n' "$(date --iso-8601=seconds)"
verify_source_authority
"${QA_PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 1440 \
    --height 1000 \
    --executable-path "${QA_CHROMIUM}" \
    --show-scrollbars \
    --screenshot "${DESKTOP_SCREENSHOT}" \
    --wait-ms 3000 | tee "${DESKTOP_QA}"
verify_source_authority

printf '\n[%s] Mobile portable HTML layout QA\n' "$(date --iso-8601=seconds)"
verify_source_authority
"${QA_PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" \
    "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 390 \
    --height 844 \
    --executable-path "${QA_CHROMIUM}" \
    --show-scrollbars \
    --screenshot "${MOBILE_SCREENSHOT}" \
    --wait-ms 3000 | tee "${MOBILE_QA}"
verify_source_authority

jq -e '.pass == true and .overlapCount == 0' "${DESKTOP_QA}" >/dev/null
jq -e '.pass == true and .overlapCount == 0' "${MOBILE_QA}" >/dev/null
chmod 0444 "${DESKTOP_QA}" "${MOBILE_QA}"

run_step "Create terminal Task-B final report release receipt" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${FINAL_RELEASE_FINALIZER}" \
    --create-report \
    --report-dir "${REPORT}" \
    --dataset-release-receipt "${FINAL_RELEASE_RECEIPT}" \
    --dataset-release-signature "${FINAL_RELEASE_SIGNATURE}" \
    --output "${REPORT_RELEASE_RECEIPT}"
require_json_pass "${REPORT_RELEASE_RECEIPT}"

printf '\n[%s] Awaiting one-time detached report signature: %s\n' \
    "$(date --iso-8601=seconds)" "${REPORT_RELEASE_SIGNATURE}"
for ((report_signature_wait=0; report_signature_wait<17280; report_signature_wait++)); do
    [[ -s "${REPORT_RELEASE_SIGNATURE}" ]] && break
    sleep 5
done
require_file "${REPORT_RELEASE_SIGNATURE}"
[[ "$(stat -c '%a' "${REPORT_RELEASE_SIGNATURE}")" == "444" ]]

run_step "Verify terminal Task-B detached report signature" \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${FINAL_RELEASE_FINALIZER}" \
    --verify-report \
    --receipt "${REPORT_RELEASE_RECEIPT}" \
    --signature "${REPORT_RELEASE_SIGNATURE}"

printf '\nM2v5 recovered completion chain PASS\n'
printf 'Source identity receipt: %s\n' "${SOURCE_RECEIPT}"
printf 'Final dataset: %s\n' "${FINAL_DATASET}/final_report_dataset.json"
printf 'Markdown report: %s\n' "${REPORT}/report.md"
printf 'Portable HTML: %s\n' "${PORTABLE_HTML}"
printf 'Signed report release receipt: %s\n' "${REPORT_RELEASE_RECEIPT}"
