#!/usr/bin/env bash

set -Eeuo pipefail

readonly REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
readonly TOPIC_ROOT="${REPO_ROOT}/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly SCRIPT_ROOT="${TOPIC_ROOT}/scripts"
readonly RESULT_ROOT="${TOPIC_ROOT}/results"
readonly WORKSPACE_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
readonly QA_PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/bin/python"
readonly PLUGIN_ROOT="/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599"

readonly MANIFEST="${RESULT_ROOT}/all_ssnv_input_manifest.json"
readonly SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
readonly COOCCURRENCE="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v1_terminal_tag"
readonly STRICT="${WORKSPACE_ROOT}/strict_methyl_candidate_confirmation_v1_terminal_tag"
readonly TUMOR_REF="${WORKSPACE_ROOT}/all_ssnv_tumor_ref_controls_v1_terminal_tag"
readonly MATCHED_RUN="${WORKSPACE_ROOT}/matched_normal_candidate_controls_v1_terminal_tag"
readonly MATCHED_ANALYSIS="${WORKSPACE_ROOT}/matched_normal_candidate_control_analysis_v1_terminal_tag"
readonly CN_CCF="${WORKSPACE_ROOT}/candidate_cn_ccf_annotations_v1_terminal_tag"
readonly PRIMARY_PRE="${RESULT_ROOT}/stable_primary_artifact_audit.v1_pre_downstream.json"
readonly PRIMARY_POST="${RESULT_ROOT}/stable_primary_artifact_audit.v1_post_downstream.json"
readonly FROZEN_POST_TERMINAL="${RESULT_ROOT}/frozen_input_immutability.post_terminal_downstream_v1.json"
readonly TERMINAL_LOG="${WORKSPACE_ROOT}/terminal_validation_chain_v10.log"

readonly SOURCE_AUDIT_ROOT="${WORKSPACE_ROOT}/tumor_ref_source_identity_audit_v3"
readonly SOURCE_SNAPSHOT="${SOURCE_AUDIT_ROOT}/observed_during_execution.snapshot.json"
readonly SOURCE_RECEIPT="${SOURCE_AUDIT_ROOT}/post_run_source_identity.receipt.json"
readonly FINAL_DATASET="${WORKSPACE_ROOT}/all_ssnv_final_report_dataset_v3_source_attested_terminal_tag"
readonly REPORT="${WORKSPACE_ROOT}/all_ssnv_final_report_v3_source_attested_terminal_tag"
readonly RELEASE_LOG="${WORKSPACE_ROOT}/source_attested_release_chain_v3.log"

readonly CLAIM_CONTRACT="${TOPIC_ROOT}/claim-contract-v5.md"
readonly EARLIER_FP_REPORT="${REPO_ROOT}/research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.md"
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

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf 'Required file is missing or empty: %s\n' "${path}" >&2
        exit 1
    }
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

run_step() {
    local label="$1"
    shift
    printf '\n[%s] %s\n' "$(date --iso-8601=seconds)" "${label}"
    printf 'Command:'
    printf ' %q' "$@"
    printf '\n'
    "$@"
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
    "${SCREEN}/all_ssnv_site_results.tsv.gz" \
    "${SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz" \
    "${SCREEN}/all_ssnv_summary.json" \
    "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    "${COOCCURRENCE}/methyl_ssnv_pair_results.tsv.gz" \
    "${COOCCURRENCE}/summary.json" \
    "${COOCCURRENCE}/run_receipt.json" \
    "${TUMOR_REF}/run_manifest.json" \
    "${PRIMARY_PRE}" \
    "${PRIMARY_POST}" \
    "${FROZEN_POST_TERMINAL}" \
    "${CN_CCF}/receipt.json" \
    "${SOURCE_SNAPSHOT}" \
    "${CLAIM_CONTRACT}" \
    "${EARLIER_FP_REPORT}" \
    "${TERMINAL_LOG}"; do
    require_file "${path}"
done

grep -Fxq 'Terminal validation chain PASS' "${TERMINAL_LOG}" || {
    printf 'Terminal validation chain has not declared PASS: %s\n' "${TERMINAL_LOG}" >&2
    exit 1
}
require_json_pass "${COOCCURRENCE}/run_receipt.json"
require_json_pass "${TUMOR_REF}/run_manifest.json"
require_json_pass "${PRIMARY_PRE}"
require_json_pass "${PRIMARY_POST}"
require_json_pass "${FROZEN_POST_TERMINAL}"
require_json_pass "${CN_CCF}/receipt.json"

if [[ -s "${MATCHED_ANALYSIS}/run_receipt.json" ]]; then
    MATCHED_FINAL="${MATCHED_ANALYSIS}"
    require_json_pass "${MATCHED_ANALYSIS}/run_receipt.json"
else
    MATCHED_FINAL="${MATCHED_RUN}"
    require_json_pass "${MATCHED_RUN}/not_applicable_receipt.json"
fi
readonly MATCHED_FINAL

for path in "${SOURCE_RECEIPT}" "${FINAL_DATASET}" "${REPORT}" "${RELEASE_LOG}"; do
    require_absent "${path}"
done

exec > >(tee "${RELEASE_LOG}") 2>&1

printf 'Input terminal screen: %s\n' "${SCREEN}"
printf 'Input terminal cooccurrence: %s\n' "${COOCCURRENCE}"
printf 'Input source snapshot: %s\n' "${SOURCE_SNAPSHOT}"
printf 'Output source receipt: %s\n' "${SOURCE_RECEIPT}"
printf 'Output release dataset: %s\n' "${FINAL_DATASET}"
printf 'Output release report: %s\n' "${REPORT}"

run_step "Verify bounded retrospective tumor-REF source identity after producer completion" \
    "${PYTHON}" "${SCRIPT_ROOT}/audit_retrospective_running_source_identity.py" verify \
    --snapshot "${SOURCE_SNAPSHOT}" \
    --run-manifest "${TUMOR_REF}/run_manifest.json" \
    --output "${SOURCE_RECEIPT}"
require_json_pass "${SOURCE_RECEIPT}"

run_step "Build source-attested final all-sSNV machine dataset" \
    "${PYTHON}" "${SCRIPT_ROOT}/build_all_ssnv_final_report_dataset.py" \
    --manifest "${MANIFEST}" \
    --screen-dir "${SCREEN}" \
    --cooccurrence-dir "${COOCCURRENCE}" \
    --strict-dir "${STRICT}" \
    --tumor-ref-dir "${TUMOR_REF}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --primary-artifact-audit-post "${PRIMARY_POST}" \
    --tumor-ref-source-identity-receipt "${SOURCE_RECEIPT}" \
    --matched-normal-dir "${MATCHED_FINAL}" \
    --cn-ccf-annotations "${CN_CCF}" \
    --output-dir "${FINAL_DATASET}"
require_json_pass "${FINAL_DATASET}/run_receipt.json"
jq -e \
    '.tumor_ref_source_identity_attestation.release_gate_pass == true and .tumor_ref_source_identity_attestation.publishable_task_b_release == true' \
    "${FINAL_DATASET}/final_report_dataset.json" >/dev/null

run_step "Build source-attested Traditional Chinese Markdown report and canonical artifact" \
    "${PYTHON}" "${SCRIPT_ROOT}/build_all_ssnv_report_artifact.py" \
    --final-dataset "${FINAL_DATASET}/final_report_dataset.json" \
    --final-receipt "${FINAL_DATASET}/run_receipt.json" \
    --candidate-catalog "${FINAL_DATASET}/candidate_catalog.tsv" \
    --candidate-witness-pairs "${FINAL_DATASET}/candidate_witness_pairs.tsv" \
    --claim-contract "${CLAIM_CONTRACT}" \
    --manifest "${MANIFEST}" \
    --screen-summary "${SCREEN}/all_ssnv_summary.json" \
    --cooccurrence-sites "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    --cooccurrence-pairs "${COOCCURRENCE}/methyl_ssnv_pair_results.tsv.gz" \
    --cooccurrence-summary "${COOCCURRENCE}/summary.json" \
    --output-reconciliation "${RESULT_ROOT}/all_ssnv_output_reconciliation.v2_verification_fix.json" \
    --post-immutability-audit "${FROZEN_POST_TERMINAL}" \
    --tree-input-audit "${RESULT_ROOT}/latest_tree_input_contract_audit.json" \
    --reference-identity-audit "${RESULT_ROOT}/extraction_reference_identity_audit.v1.json" \
    --earlier-fp-report "${EARLIER_FP_REPORT}" \
    --repo-root "${REPO_ROOT}" \
    --output-dir "${REPORT}"
require_json_pass "${REPORT}/report_build_receipt.json"

run_step "Package and officially verify source-attested portable HTML" \
    node "${SCRIPT_ROOT}/deliver_portable_report_scrollbar_safe.mjs" \
    --input "${REPORT}/artifact.json" \
    --output "${PORTABLE_HTML}" \
    --receipt "${PORTABLE_RECEIPT}" \
    --screenshot "${PORTABLE_VERIFY_SCREENSHOT}" \
    --plugin-root "${PLUGIN_ROOT}" \
    --ready-timeout-ms 120000 \
    --action-timeout-ms 60000 \
    --timeout-ms 300000
require_json_pass "${PORTABLE_RECEIPT}"

printf '\n[%s] Desktop source-attested portable HTML layout QA\n' "$(date --iso-8601=seconds)"
"${QA_PYTHON}" "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 1440 \
    --height 1000 \
    --show-scrollbars \
    --screenshot "${DESKTOP_SCREENSHOT}" \
    --wait-ms 3000 | tee "${DESKTOP_QA}"

printf '\n[%s] Mobile source-attested portable HTML layout QA\n' "$(date --iso-8601=seconds)"
"${QA_PYTHON}" "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 390 \
    --height 844 \
    --show-scrollbars \
    --screenshot "${MOBILE_SCREENSHOT}" \
    --wait-ms 3000 | tee "${MOBILE_QA}"

jq -e '.pass == true and .overlapCount == 0' "${DESKTOP_QA}" >/dev/null
jq -e '.pass == true and .overlapCount == 0' "${MOBILE_QA}" >/dev/null

printf '\nSource-attested release chain PASS\n'
printf 'Source identity receipt: %s\n' "${SOURCE_RECEIPT}"
printf 'Final dataset: %s\n' "${FINAL_DATASET}/final_report_dataset.json"
printf 'Markdown report: %s\n' "${REPORT}/report.md"
printf 'Portable HTML: %s\n' "${PORTABLE_HTML}"
