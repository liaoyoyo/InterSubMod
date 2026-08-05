#!/usr/bin/env bash

set -Eeuo pipefail

readonly REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
readonly TOPIC_ROOT="${REPO_ROOT}/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly SCRIPT_ROOT="${TOPIC_ROOT}/scripts"
readonly RESULT_ROOT="${TOPIC_ROOT}/results"
readonly WORKSPACE_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
readonly QA_PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/bin/python"
readonly INTERSUBMOD_BIN="${REPO_ROOT}/build/bin/inter_sub_mod"
readonly REFERENCE="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
readonly PLUGIN_ROOT="/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599"

readonly MANIFEST="${RESULT_ROOT}/all_ssnv_input_manifest.json"
readonly SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
readonly SCREEN_SITES="${SCREEN}/all_ssnv_site_results.tsv.gz"
readonly SCREEN_ASSIGNMENTS="${SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
readonly SCREEN_SUMMARY="${SCREEN}/all_ssnv_summary.json"
readonly EXTRACT_ROOT="${WORKSPACE_ROOT}/intersubmod_all_ssnv_v2_verification_fix"
readonly PRIMARY_PRE="${RESULT_ROOT}/stable_primary_artifact_audit.v1_pre_downstream.json"

readonly COOCCURRENCE="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked"
readonly TUMOR_REF="${WORKSPACE_ROOT}/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel"
readonly STRICT="${WORKSPACE_ROOT}/strict_methyl_candidate_confirmation_v2_m2v5"
readonly MATCHED_RUN="${WORKSPACE_ROOT}/matched_normal_candidate_controls_v2_m2v5"
readonly MATCHED_ANALYSIS="${WORKSPACE_ROOT}/matched_normal_candidate_control_analysis_v2_m2v5"
readonly CN_CCF="${WORKSPACE_ROOT}/candidate_cn_ccf_annotations_v2_m2v5"

readonly SOURCE_SNAPSHOT="${WORKSPACE_ROOT}/tumor_ref_recovery_source_identity_v1/observed_during_execution.snapshot.json"
readonly SOURCE_RECEIPT="${WORKSPACE_ROOT}/tumor_ref_recovery_source_identity_v1/post_run_source_identity.receipt.json"
readonly INDEPENDENT_M2_AUDIT="${RESULT_ROOT}/independent_m2_gate_recount.v3.json"
readonly COOCCURRENCE_PREFLIGHT="${RESULT_ROOT}/cooccurrence_task_contract_preflight.v5_dependency_attested_raw_identity_full_runtime.json"
readonly PRIMARY_POST="${RESULT_ROOT}/stable_primary_artifact_audit.v2_post_m2v5.json"
readonly FROZEN_POST="${RESULT_ROOT}/frozen_input_immutability.post_m2v5_downstream_v2.json"

readonly FINAL_DATASET="${WORKSPACE_ROOT}/all_ssnv_final_report_dataset_v4_m2v5_source_attested"
readonly REPORT="${WORKSPACE_ROOT}/all_ssnv_final_report_v4_m2v5_source_attested"
readonly CLAIM_CONTRACT="${TOPIC_ROOT}/claim-contract-v5.md"
readonly EARLIER_FP_REPORT="${REPO_ROOT}/research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.md"
readonly LOG_PATH="${WORKSPACE_ROOT}/m2v5_recovered_completion_chain_v1.log"

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

require_cooccurrence_preflight_reconciliation() {
    jq -e \
        --slurpfile receipt "${COOCCURRENCE}/run_receipt.json" \
        --slurpfile summary "${COOCCURRENCE}/summary.json" \
        '
        .code.source_identity_after.analyzer == $receipt[0].code.analyzer and
        .code.source_identity_after.ssnv_cooccurrence_lib == $receipt[0].code.ssnv_cooccurrence_lib and
        .code.source_identity_after.latest_tag_join == $receipt[0].code.latest_tag_join and
        .code.source_identity_after.m2_screen_gate == $receipt[0].code.m2_screen_gate and
        .observed.raw_bam_identity_recovery.equivalence_contract ==
            $summary[0].raw_bam_identity_recovery_audit.equivalence_contract and
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
        .observed.raw_bam_identity_recovery.conflicting_analysis_payload_policy ==
            $receipt[0].counts.raw_identity_conflicting_analysis_payload_policy and
        .observed.raw_bam_identity_recovery.failure_counts_materialized == false and
        $receipt[0].counts.raw_identity_failure_counts_materialized == false and
        $receipt[0].source_lock.all_sources_read_only_and_unchanged == true
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
    "${SCREEN_SITES}" \
    "${SCREEN_ASSIGNMENTS}" \
    "${SCREEN_SUMMARY}" \
    "${PRIMARY_PRE}" \
    "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    "${COOCCURRENCE}/methyl_ssnv_pair_results.tsv.gz" \
    "${COOCCURRENCE}/raw_identity_duplicate_audit.tsv.gz" \
    "${COOCCURRENCE}/summary.json" \
    "${COOCCURRENCE}/run_receipt.json" \
    "${TUMOR_REF}/all_ssnv_tumor_ref_control_site_results.tsv.gz" \
    "${TUMOR_REF}/summary.json" \
    "${TUMOR_REF}/run_manifest.json" \
    "${SOURCE_SNAPSHOT}" \
    "${INDEPENDENT_M2_AUDIT}" \
    "${COOCCURRENCE_PREFLIGHT}" \
    "${CLAIM_CONTRACT}" \
    "${EARLIER_FP_REPORT}" \
    "${INTERSUBMOD_BIN}" \
    "${REFERENCE}" \
    "${REFERENCE}.fai"; do
    require_file "${path}"
done
require_json_pass "${PRIMARY_PRE}"
require_json_pass "${COOCCURRENCE}/run_receipt.json"
require_json_pass "${TUMOR_REF}/run_manifest.json"
require_json_pass "${SOURCE_SNAPSHOT}"
require_json_pass "${INDEPENDENT_M2_AUDIT}"
require_json_pass "${COOCCURRENCE_PREFLIGHT}"
require_cooccurrence_preflight_reconciliation

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
    "${LOG_PATH}"; do
    require_absent "${path}"
done

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
    "${PYTHON}" "${SCRIPT_ROOT}/verify_retrospective_running_source_identity_v2.py" verify \
    --snapshot "${SOURCE_SNAPSHOT}" \
    --run-manifest "${TUMOR_REF}/run_manifest.json" \
    --output "${SOURCE_RECEIPT}"
require_json_pass "${SOURCE_RECEIPT}"

run_step "Run strict multi-seed and column-null robustness" \
    "${PYTHON}" "${SCRIPT_ROOT}/run_strict_methyl_candidate_confirmation.py" \
    --candidate-table "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    --assignments "${SCREEN_ASSIGNMENTS}" \
    --cooccurrence-receipt "${COOCCURRENCE}/run_receipt.json" \
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
    "${PYTHON}" "${SCRIPT_ROOT}/run_matched_normal_candidate_controls.py" \
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
        "${PYTHON}" "${SCRIPT_ROOT}/analyze_matched_normal_candidate_controls.py" \
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
    "${PYTHON}" "${SCRIPT_ROOT}/annotate_candidate_cn_ccf.py" \
    --input "${COOCCURRENCE}/methyl_ssnv_site_results.tsv.gz" \
    --selection-column "${FORMAL_SELECTION_COLUMN}" \
    --selection-value true \
    --output-dir "${CN_CCF}"
require_json_pass "${CN_CCF}/receipt.json"

POST_AUDIT_COMMAND=(
    "${PYTHON}" "${SCRIPT_ROOT}/audit_stable_primary_artifacts.py"
    --site-results "${SCREEN_SITES}"
    --stable-assignments "${SCREEN_ASSIGNMENTS}"
    --consumer-receipt "${COOCCURRENCE}/run_receipt.json"
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
    "${PYTHON}" "${SCRIPT_ROOT}/audit_frozen_input_immutability.py" \
    --manifest "${MANIFEST}" \
    --output "${FROZEN_POST}"
require_json_pass "${FROZEN_POST}"

run_step "Build source-attested final all-sSNV dataset" \
    "${PYTHON}" "${SCRIPT_ROOT}/build_all_ssnv_final_report_dataset.py" \
    --manifest "${MANIFEST}" \
    --screen-dir "${SCREEN}" \
    --cooccurrence-dir "${COOCCURRENCE}" \
    --strict-dir "${STRICT}" \
    --tumor-ref-dir "${TUMOR_REF}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --primary-artifact-audit-post "${PRIMARY_POST}" \
    --tumor-ref-source-identity-receipt "${SOURCE_RECEIPT}" \
    --independent-m2-audit "${INDEPENDENT_M2_AUDIT}" \
    --matched-normal-dir "${MATCHED_FINAL}" \
    --cn-ccf-annotations "${CN_CCF}" \
    --output-dir "${FINAL_DATASET}"
require_json_pass "${FINAL_DATASET}/run_receipt.json"
jq -e \
    '.tumor_ref_source_identity_attestation.release_gate_pass == true and .m2_evaluability_contract.independent_logic_audit.status == "PASS_LOGIC_INDEPENDENT_RECOUNT"' \
    "${FINAL_DATASET}/final_report_dataset.json" >/dev/null

run_step "Build complete Traditional Chinese report and canonical artifact" \
    "${PYTHON}" "${SCRIPT_ROOT}/build_all_ssnv_report_artifact.py" \
    --final-dataset "${FINAL_DATASET}/final_report_dataset.json" \
    --final-receipt "${FINAL_DATASET}/run_receipt.json" \
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

printf '\n[%s] Desktop portable HTML layout QA\n' "$(date --iso-8601=seconds)"
"${QA_PYTHON}" "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 1440 \
    --height 1000 \
    --show-scrollbars \
    --screenshot "${DESKTOP_SCREENSHOT}" \
    --wait-ms 3000 | tee "${DESKTOP_QA}"

printf '\n[%s] Mobile portable HTML layout QA\n' "$(date --iso-8601=seconds)"
"${QA_PYTHON}" "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 390 \
    --height 844 \
    --show-scrollbars \
    --screenshot "${MOBILE_SCREENSHOT}" \
    --wait-ms 3000 | tee "${MOBILE_QA}"

jq -e '.pass == true and .overlapCount == 0' "${DESKTOP_QA}" >/dev/null
jq -e '.pass == true and .overlapCount == 0' "${MOBILE_QA}" >/dev/null

printf '\nM2v5 recovered completion chain PASS\n'
printf 'Source identity receipt: %s\n' "${SOURCE_RECEIPT}"
printf 'Final dataset: %s\n' "${FINAL_DATASET}/final_report_dataset.json"
printf 'Markdown report: %s\n' "${REPORT}/report.md"
printf 'Portable HTML: %s\n' "${PORTABLE_HTML}"
