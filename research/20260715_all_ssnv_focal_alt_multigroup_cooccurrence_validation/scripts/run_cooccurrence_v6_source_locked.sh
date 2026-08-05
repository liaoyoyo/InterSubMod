#!/usr/bin/env bash

set -Eeuo pipefail

export PATH="/usr/bin:/bin"
unset BASH_ENV ENV CDPATH GLOBIGNORE LD_AUDIT LD_DEBUG LD_LIBRARY_PATH LD_PRELOAD
unset OPENSSL_CONF OPENSSL_MODULES PYTHONHOME PYTHONINSPECT PYTHONPATH PYTHONSTARTUP
for inherited_function in jq stat sha256sum realpath sort awk git openssl; do
    unset -f "${inherited_function}" 2>/dev/null || true
done

readonly REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
readonly TOPIC_ROOT="${REPO_ROOT}/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly SCRIPT_ROOT="${TOPIC_ROOT}/scripts"
readonly RESULT_ROOT="${TOPIC_ROOT}/results"
readonly WORKSPACE_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"

readonly ANALYZER="${SCRIPT_ROOT}/analyze_methyl_ssnv_cooccurrence.py"
readonly PREFLIGHT_SOURCE="${SCRIPT_ROOT}/audit_cooccurrence_task_contract_preflight.py"
readonly RELEASE_FINALIZER="${SCRIPT_ROOT}/finalize_cooccurrence_release_receipt.py"
readonly RUNNER_SOURCE="${SCRIPT_ROOT}/run_cooccurrence_v6_source_locked.sh"
readonly SOURCE_AUTHORITY_VALIDATOR="${SCRIPT_ROOT}/release_source_authority.py"
readonly SOURCE_AUTHORITY="${REPO_ROOT}/docs/provenance/source_authorities/20260722_all_ssnv_focal_alt_release_source_authority.v7.json"
readonly SOURCE_AUTHORITY_APPROVAL="${REPO_ROOT}/docs/provenance/source_authorities/20260722_all_ssnv_focal_alt_release_source_authority.v7.approval.json"
readonly SOURCE_AUTHORITY_SIGNATURE="${SOURCE_AUTHORITY_APPROVAL}.ed25519.sig"
readonly SOURCE_AUTHORITY_ID="20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
readonly SOURCE_AUTHORITY_PUBLIC_KEY="/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
readonly SOURCE_AUTHORITY_PRIVATE_KEY="/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_private_encrypted_unrecoverable_after_signing.pem"
readonly SOURCE_AUTHORITY_PUBLIC_KEY_SHA256="cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
readonly OPENSSL="/usr/bin/openssl"
readonly OPENSSL_SHA256="38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
readonly COOCCURRENCE_LIB="${REPO_ROOT}/scripts/ssnv_cooccurrence_lib.py"
readonly LATEST_TAG_JOIN="${SCRIPT_ROOT}/latest_tag_join.py"
readonly M2_SCREEN_GATE="${SCRIPT_ROOT}/m2_screen_gate.py"
readonly MANIFEST="${RESULT_ROOT}/all_ssnv_input_manifest.json"
readonly SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
readonly ASSIGNMENTS="${SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
readonly SITES="${SCREEN}/all_ssnv_site_results.tsv.gz"
readonly PRIMARY_PRE="${RESULT_ROOT}/stable_primary_artifact_audit.v6_strict_command_parity_pre_downstream.json"
readonly INDEPENDENT_M2_AUDIT="${RESULT_ROOT}/independent_m2_gate_recount.v3.json"
readonly INTERSUBMOD_ROOT="${WORKSPACE_ROOT}/intersubmod_all_ssnv_v2_verification_fix"
readonly PREFLIGHT="${RESULT_ROOT}/cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
readonly OUTPUT="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
readonly LOG="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity.log"
readonly RELEASE_RECEIPT="${OUTPUT}/release_receipt.json"
readonly PYTHON_CACHE_ROOT="${WORKSPACE_ROOT}/.python_cache_cooccurrence_v8_command_parity_bound_bootstrap"

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
        .task_type == "B_comprehensive_validation" and
        .scope == "all_7_datasets_469849_sites_and_release_chain" and
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

verify_receipt_artifact() {
    local receipt="$1"
    local key="$2"
    local expected_path="$3"
    local declared_path declared_size declared_sha256
    declared_path="$(jq -er --arg key "${key}" '.outputs[$key].path' "${receipt}")"
    declared_size="$(jq -er --arg key "${key}" '.outputs[$key].size_bytes' "${receipt}")"
    declared_sha256="$(jq -er --arg key "${key}" '.outputs[$key].sha256' "${receipt}")"
    [[ "${declared_path}" == "$(realpath "${expected_path}")" ]]
    [[ "${declared_size}" == "$(stat -c '%s' "${expected_path}")" ]]
    [[ "${declared_sha256}" == "$(sha256sum "${expected_path}" | awk '{print $1}')" ]]
}

on_error() {
    local exit_code=$?
    printf 'FAILED exit=%d line=%d command=%q\n' \
        "${exit_code}" "${BASH_LINENO[0]}" "${BASH_COMMAND}" >&2
    exit "${exit_code}"
}
trap on_error ERR

cd "${REPO_ROOT}"
umask 0022

for path in \
    "${ANALYZER}" \
    "${RELEASE_FINALIZER}" \
    "${RUNNER_SOURCE}" \
    "${SOURCE_AUTHORITY_VALIDATOR}" \
    "${SOURCE_AUTHORITY}" \
    "${SOURCE_AUTHORITY_APPROVAL}" \
    "${SOURCE_AUTHORITY_SIGNATURE}" \
    "${SOURCE_AUTHORITY_PUBLIC_KEY}" \
    "${SOURCE_AUTHORITY_PRIVATE_KEY}" \
    "${MANIFEST}" \
    "${ASSIGNMENTS}" \
    "${SITES}" \
    "${PRIMARY_PRE}" \
    "${INDEPENDENT_M2_AUDIT}" \
    "${PREFLIGHT}"; do
    require_file "${path}"
done
verify_source_authority
[[ ! -e "${OUTPUT}" ]] || {
    printf 'Refusing to overwrite output: %s\n' "${OUTPUT}" >&2
    exit 1
}
[[ ! -e "${LOG}" ]] || {
    printf 'Refusing to overwrite log: %s\n' "${LOG}" >&2
    exit 1
}
[[ ! -e "${PYTHON_CACHE_ROOT}" ]] || {
    printf 'Refusing to reuse Python cache: %s\n' "${PYTHON_CACHE_ROOT}" >&2
    exit 1
}
umask 0077
/usr/bin/mkdir --mode=0700 -- "${PYTHON_CACHE_ROOT}"
[[ "$(stat -c '%a' "${PYTHON_CACHE_ROOT}")" == "700" ]]
umask 0022
jq -e '.pass == true' "${PRIMARY_PRE}" >/dev/null
jq -e \
    '
    .schema_name == "intersubmod.cooccurrence_task_contract_preflight" and
    .schema_version == "2.0.0" and
    .task_type == "B_comprehensive_validation" and
    .pass == true and
    .raw_identity_release_contract == "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3" and
    .observed.task_count == 102842 and
    .observed.raw_bam_identity_recovery.aggregate.site_tasks == 102842 and
    .observed.raw_bam_identity_recovery.all_result_rows_passed_invariant_validation == true and
    .observed.raw_bam_identity_recovery.equivalence_contract == "sam_core_and_all_aux_tags_except_RG_exact_v1" and
    .observed.raw_bam_identity_recovery.allowed_differing_auxiliary_tags == ["RG"] and
    .observed.raw_bam_identity_recovery.analysis_scope_identity_contract == "selected_focal_ref_alt_alignments_full_sam_core_all_typed_aux_except_RG_latest_HP_PS_assignment_label_and_partner_calls_v1" and
    .observed.raw_bam_identity_recovery.missing_projection_policy == "hard_fail_before_site_result" and
    .observed.raw_bam_identity_recovery.conflicting_analysis_payload_policy == "hard_fail_before_site_result" and
    .observed.raw_bam_identity_recovery.failure_counts_materialized == false and
    .code.source_identity_before == .code.source_identity_after and
    (.code.source_identity_before | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "preflight", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
    (.code.source_modes_before | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "preflight", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
    (.code.source_modes_after | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "preflight", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
    ([.code.source_modes_before[], .code.source_modes_after[]] | all(. == "0o444"))
    ' "${PREFLIGHT}" >/dev/null

for source_name in preflight analyzer ssnv_cooccurrence_lib latest_tag_join m2_screen_gate release_finalizer release_runner source_authority_validator; do
    source_path="$(jq -er --arg name "${source_name}" '.code.source_identity_after[$name].path' "${PREFLIGHT}")"
    source_size="$(jq -er --arg name "${source_name}" '.code.source_identity_after[$name].size_bytes' "${PREFLIGHT}")"
    source_sha256="$(jq -er --arg name "${source_name}" '.code.source_identity_after[$name].sha256' "${PREFLIGHT}")"
    require_file "${source_path}"
    [[ "$(stat -c '%s' "${source_path}")" == "${source_size}" ]] || {
        printf 'Source size differs from preflight: %s\n' "${source_path}" >&2
        exit 1
    }
    [[ "$(sha256sum "${source_path}" | awk '{print $1}')" == "${source_sha256}" ]] || {
        printf 'Source SHA-256 differs from preflight: %s\n' "${source_path}" >&2
        exit 1
    }
    [[ "$(stat -c '%a' "${source_path}")" == "444" ]] || {
        printf 'Source is not locked mode 0444: %s\n' "${source_path}" >&2
        exit 1
    }
done
for source_binding in \
    "preflight:${PREFLIGHT_SOURCE}" \
    "analyzer:${ANALYZER}" \
    "ssnv_cooccurrence_lib:${COOCCURRENCE_LIB}" \
    "latest_tag_join:${LATEST_TAG_JOIN}" \
    "m2_screen_gate:${M2_SCREEN_GATE}" \
    "release_finalizer:${RELEASE_FINALIZER}" \
    "release_runner:${RUNNER_SOURCE}" \
    "source_authority_validator:${SOURCE_AUTHORITY_VALIDATOR}"; do
    source_name="${source_binding%%:*}"
    trusted_path="${source_binding#*:}"
    [[ "$(jq -er --arg name "${source_name}" '.code.source_identity_after[$name].path' "${PREFLIGHT}")" == "$(realpath "${trusted_path}")" ]] || {
        printf 'Preflight source path differs from trusted role %s\n' "${source_name}" >&2
        exit 1
    }
done

exec > >(tee "${LOG}") 2>&1

printf 'Input manifest: %s\n' "${MANIFEST}"
printf 'Input assignments: %s\n' "${ASSIGNMENTS}"
printf 'Input sites: %s\n' "${SITES}"
printf 'Input full raw-identity preflight: %s\n' "${PREFLIGHT}"
printf 'Output cooccurrence: %s\n' "${OUTPUT}"
printf 'Command:'
printf ' %q' \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${ANALYZER}" \
    --manifest "${MANIFEST}" \
    --assignments "${ASSIGNMENTS}" \
    --sites "${SITES}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --independent-m2-audit "${INDEPENDENT_M2_AUDIT}" \
    --preflight-receipt "${PREFLIGHT}" \
    --intersubmod-root "${INTERSUBMOD_ROOT}" \
    --output-dir "${OUTPUT}" \
    --workers 40 \
    --chunk-size 8 \
    --max-pending-chunks 80 \
    --permutations 999 \
    --top-markers 3 \
    --exact-state-space-ceiling 250000
printf '\n'

"${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${ANALYZER}" \
    --manifest "${MANIFEST}" \
    --assignments "${ASSIGNMENTS}" \
    --sites "${SITES}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --independent-m2-audit "${INDEPENDENT_M2_AUDIT}" \
    --preflight-receipt "${PREFLIGHT}" \
    --intersubmod-root "${INTERSUBMOD_ROOT}" \
    --output-dir "${OUTPUT}" \
    --workers 40 \
    --chunk-size 8 \
    --max-pending-chunks 80 \
    --permutations 999 \
    --top-markers 3 \
    --exact-state-space-ceiling 250000

for path in \
    "${OUTPUT}/methyl_ssnv_site_results.tsv.gz" \
    "${OUTPUT}/methyl_ssnv_pair_results.tsv.gz" \
    "${OUTPUT}/raw_identity_duplicate_audit.tsv.gz" \
    "${OUTPUT}/oracle_cases.json" \
    "${OUTPUT}/summary.json" \
    "${OUTPUT}/run_receipt.json"; do
    require_file "${path}"
done
jq -e \
    --arg preflight "$(realpath "${PREFLIGHT}")" \
    '
    .pass == true and
    .task_type == "B_comprehensive_validation" and
    .scope == "all_manifest_samples" and
    .full_scope == true and
    .release_status == "PRODUCER_PASS_PENDING_RUNNER_RECONCILIATION" and
    .pass_semantics == "execution_integrity_only_pending_release_reconciliation" and
    .counts.stable_sites == 102842 and
    .counts.raw_identity_all_site_results_passed_invariant_validation == true and
    .counts.raw_identity_missing_projection_policy == "hard_fail_before_site_result" and
    .counts.raw_identity_conflicting_analysis_payload_policy == "hard_fail_before_site_result" and
    .counts.raw_identity_failure_counts_materialized == false and
    .raw_identity_release_contract == "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3" and
    .preflight_gate.schema_name == "intersubmod.cooccurrence_task_contract_preflight" and
    .preflight_gate.schema_version == "2.0.0" and
    .preflight_gate.task_count == 102842 and
    .preflight_gate.pass == true and
    .preflight_gate.raw_identity_release_contract == "typed_aux_dependency_attested_sparse_audit_analysis_payload_v3" and
    .inputs.preflight_receipt.path == $preflight and
    (.code | keys | sort) == ["analyzer", "latest_tag_join", "m2_screen_gate", "release_finalizer", "release_runner", "source_authority_validator", "ssnv_cooccurrence_lib"] and
    .source_lock.all_sources_read_only_and_unchanged == true and
    .source_lock.source_identity_before == .source_lock.source_identity_after_compute and
    .source_lock.source_identity_before == .source_lock.source_identity_after_output and
    .source_lock.source_identity_before == .code and
    ([.source_lock.source_modes_before[], .source_lock.source_modes_after_compute[], .source_lock.source_modes_after_output[]] | all(. == "0o444"))
    ' "${OUTPUT}/run_receipt.json" >/dev/null

verify_receipt_artifact "${OUTPUT}/run_receipt.json" pair_tsv "${OUTPUT}/methyl_ssnv_pair_results.tsv.gz"
verify_receipt_artifact "${OUTPUT}/run_receipt.json" site_tsv "${OUTPUT}/methyl_ssnv_site_results.tsv.gz"
verify_receipt_artifact "${OUTPUT}/run_receipt.json" summary_json "${OUTPUT}/summary.json"
verify_receipt_artifact "${OUTPUT}/run_receipt.json" case_json "${OUTPUT}/oracle_cases.json"
verify_receipt_artifact "${OUTPUT}/run_receipt.json" raw_identity_duplicate_audit_tsv "${OUTPUT}/raw_identity_duplicate_audit.tsv.gz"
jq -e \
    --slurpfile receipt "${OUTPUT}/run_receipt.json" \
    --slurpfile summary "${OUTPUT}/summary.json" \
    '
    .code.source_identity_after.analyzer == $receipt[0].code.analyzer and
    .code.source_identity_after.ssnv_cooccurrence_lib == $receipt[0].code.ssnv_cooccurrence_lib and
    .code.source_identity_after.latest_tag_join == $receipt[0].code.latest_tag_join and
    .code.source_identity_after.m2_screen_gate == $receipt[0].code.m2_screen_gate and
    .code.source_identity_after.release_finalizer == $receipt[0].code.release_finalizer and
    .code.source_identity_after.release_runner == $receipt[0].code.release_runner and
    .code.source_identity_after.source_authority_validator == $receipt[0].code.source_authority_validator and
    .source_authority == $receipt[0].source_authority and
    .observed.raw_bam_identity_recovery.site_weighted_audit_sha256 ==
        $summary[0].raw_bam_identity_recovery_audit.site_weighted_audit_sha256 and
    .observed.raw_bam_identity_recovery.aggregate.expected_projection_occurrences ==
        $summary[0].raw_bam_identity_recovery_audit.n_expected_projection_occurrences and
    .observed.raw_bam_identity_recovery.aggregate.duplicate_projection_occurrences_collapsed ==
        $summary[0].raw_bam_identity_recovery_audit.n_duplicate_projection_occurrences_collapsed and
    .observed.raw_bam_identity_recovery.aggregate.duplicate_extra_record_occurrences_collapsed ==
        $summary[0].raw_bam_identity_recovery_audit.n_duplicate_extra_record_occurrences_collapsed and
    .observed.raw_bam_identity_recovery.aggregate.duplicate_projection_occurrences_collapsed ==
        $receipt[0].counts.raw_identity_sparse_duplicate_rows and
    $summary[0].raw_bam_identity_recovery_audit.missing_projection_policy ==
        .observed.raw_bam_identity_recovery.missing_projection_policy and
    $summary[0].raw_bam_identity_recovery_audit.conflicting_analysis_payload_policy ==
        .observed.raw_bam_identity_recovery.conflicting_analysis_payload_policy
    ' "${PREFLIGHT}" >/dev/null

printf 'Command:'
printf ' %q' \
    "${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${RELEASE_FINALIZER}" \
    --preflight "${PREFLIGHT}" \
    --producer-receipt "${OUTPUT}/run_receipt.json" \
    --summary "${OUTPUT}/summary.json" \
    --sites "${OUTPUT}/methyl_ssnv_site_results.tsv.gz" \
    --pairs "${OUTPUT}/methyl_ssnv_pair_results.tsv.gz" \
    --duplicates "${OUTPUT}/raw_identity_duplicate_audit.tsv.gz" \
    --oracle "${OUTPUT}/oracle_cases.json" \
    --runner-source "${RUNNER_SOURCE}" \
    --output "${RELEASE_RECEIPT}"
printf '\n'
"${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" "${RELEASE_FINALIZER}" \
    --preflight "${PREFLIGHT}" \
    --producer-receipt "${OUTPUT}/run_receipt.json" \
    --summary "${OUTPUT}/summary.json" \
    --sites "${OUTPUT}/methyl_ssnv_site_results.tsv.gz" \
    --pairs "${OUTPUT}/methyl_ssnv_pair_results.tsv.gz" \
    --duplicates "${OUTPUT}/raw_identity_duplicate_audit.tsv.gz" \
    --oracle "${OUTPUT}/oracle_cases.json" \
    --runner-source "${RUNNER_SOURCE}" \
    --output "${RELEASE_RECEIPT}"
require_file "${RELEASE_RECEIPT}"
jq -e \
    --slurpfile preflight "${PREFLIGHT}" \
    '.pass == true and
     .schema_name == "intersubmod.cooccurrence_release_receipt" and
     .schema_version == "1.0.0" and
     .task_type == "B_comprehensive_validation" and
     .scope == "all_manifest_samples" and
     .release_status == "RELEASE_RECONCILIATION_PASS" and
     .code.release_finalizer == $preflight[0].code.source_identity_after.release_finalizer and
     .code.release_runner == $preflight[0].code.source_identity_after.release_runner and
     ([.checks[]] | all(. == true))' \
    "${RELEASE_RECEIPT}" >/dev/null

printf 'Cooccurrence v7 summary-hotfix PASS: %s\n' "${OUTPUT}"
