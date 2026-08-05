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
readonly EXTRACT_ROOT="${WORKSPACE_ROOT}/intersubmod_all_ssnv_v2_verification_fix"
readonly PREFIX_SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v9_first_six_thread_pinned_source_locked"
readonly REPLACEMENT_SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v6_hcc1954_seed_parallel_200"
readonly EQUIVALENCE_RECEIPT="${RESULT_ROOT}/phylo_parallel_exact_equivalence.pinned_390228ce.v1.json"
readonly SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"

readonly SCREEN_SITES="${SCREEN}/all_ssnv_site_results.tsv.gz"
readonly SCREEN_ASSIGNMENTS="${SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
readonly SCREEN_SUMMARY="${SCREEN}/all_ssnv_summary.json"
readonly FROZEN_POST_SCREEN="${RESULT_ROOT}/frozen_input_immutability.post_screen_recovery_v1.json"
readonly PRIMARY_PRE="${RESULT_ROOT}/stable_primary_artifact_audit.v1_pre_downstream.json"

readonly TUMOR_REF="${WORKSPACE_ROOT}/all_ssnv_tumor_ref_controls_v1_terminal_tag"
readonly COOCCURRENCE="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v1_terminal_tag"
readonly STRICT="${WORKSPACE_ROOT}/strict_methyl_candidate_confirmation_v1_terminal_tag"
readonly MATCHED_RUN="${WORKSPACE_ROOT}/matched_normal_candidate_controls_v1_terminal_tag"
readonly MATCHED_ANALYSIS="${WORKSPACE_ROOT}/matched_normal_candidate_control_analysis_v1_terminal_tag"
readonly CN_CCF="${WORKSPACE_ROOT}/candidate_cn_ccf_annotations_v1_terminal_tag"
readonly PRIMARY_POST="${RESULT_ROOT}/stable_primary_artifact_audit.v1_post_downstream.json"
readonly FROZEN_POST_TERMINAL="${RESULT_ROOT}/frozen_input_immutability.post_terminal_downstream_v1.json"
readonly FINAL_DATASET="${WORKSPACE_ROOT}/all_ssnv_final_report_dataset_v2_terminal_tag"
readonly REPORT="${WORKSPACE_ROOT}/all_ssnv_final_report_v2_terminal_tag"

readonly FORMAL_SELECTION_COLUMN="multi_marker_molecular_haplotype_base_candidate"
readonly TERMINAL_CLAIM_CONTRACT="${TOPIC_ROOT}/claim-contract-v3.md"
readonly EARLIER_FP_REPORT="${REPO_ROOT}/research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.md"
readonly LOG_PATH="${WORKSPACE_ROOT}/terminal_validation_chain_v10.log"

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
    printf '\nFAILED exit=%d line=%d command=%q\n' "${exit_code}" "${BASH_LINENO[0]}" "${BASH_COMMAND}" >&2
    exit "${exit_code}"
}
trap on_error ERR

cd "${REPO_ROOT}"
umask 0022

for path in \
    "${MANIFEST}" \
    "${PREFIX_SCREEN}/all_ssnv_site_results.tsv.gz" \
    "${PREFIX_SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz" \
    "${PREFIX_SCREEN}/all_ssnv_summary.json" \
    "${PREFIX_SCREEN}/run_manifest.json" \
    "${PREFIX_SCREEN}/source_lock_receipt.json" \
    "${REPLACEMENT_SCREEN}/all_ssnv_site_results.tsv.gz" \
    "${REPLACEMENT_SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz" \
    "${REPLACEMENT_SCREEN}/all_ssnv_summary.json" \
    "${REPLACEMENT_SCREEN}/run_manifest.json" \
    "${EQUIVALENCE_RECEIPT}" \
    "${INTERSUBMOD_BIN}" \
    "${REFERENCE}" \
    "${REFERENCE}.fai" \
    "${TERMINAL_CLAIM_CONTRACT}" \
    "${EARLIER_FP_REPORT}"; do
    require_file "${path}"
done
require_json_pass "${PREFIX_SCREEN}/run_manifest.json"
require_json_pass "${PREFIX_SCREEN}/source_lock_receipt.json"
require_json_pass "${REPLACEMENT_SCREEN}/run_manifest.json"
require_json_pass "${EQUIVALENCE_RECEIPT}"

for path in \
    "${SCREEN}" \
    "${FROZEN_POST_SCREEN}" \
    "${PRIMARY_PRE}" \
    "${TUMOR_REF}" \
    "${COOCCURRENCE}" \
    "${STRICT}" \
    "${MATCHED_RUN}" \
    "${MATCHED_ANALYSIS}" \
    "${CN_CCF}" \
    "${PRIMARY_POST}" \
    "${FROZEN_POST_TERMINAL}" \
    "${FINAL_DATASET}" \
    "${REPORT}" \
    "${LOG_PATH}"; do
    require_absent "${path}"
done

exec > >(tee "${LOG_PATH}") 2>&1

printf 'Input prefix screen: %s\n' "${PREFIX_SCREEN}"
printf 'Input replacement screen: %s\n' "${REPLACEMENT_SCREEN}"
printf 'Input manifest: %s\n' "${MANIFEST}"
printf 'Output terminal screen: %s\n' "${SCREEN}"
printf 'Output report: %s\n' "${REPORT}"
printf 'Computational thread environment: OMP=%s OPENBLAS=%s MKL=%s NUMEXPR=%s BLIS=%s\n' \
    "${OMP_NUM_THREADS}" "${OPENBLAS_NUM_THREADS}" "${MKL_NUM_THREADS}" \
    "${NUMEXPR_NUM_THREADS}" "${BLIS_NUM_THREADS}"

run_step "Merge exact-equivalent HCC1954 replacement into source-locked first-six screen" \
    "${PYTHON}" "${SCRIPT_ROOT}/merge_all_ssnv_screen_recovery.py" \
    --manifest "${MANIFEST}" \
    --intersubmod-root "${EXTRACT_ROOT}" \
    --prefix-dir "${PREFIX_SCREEN}" \
    --replacement-dir "${REPLACEMENT_SCREEN}" \
    --replacement-sample HCC1954 \
    --equivalence-receipt "${EQUIVALENCE_RECEIPT}" \
    --output-dir "${SCREEN}"
require_json_pass "${SCREEN}/run_manifest.json"
jq -e '.counts.processed_sites == 469849' "${SCREEN}/run_manifest.json" >/dev/null

run_step "Audit frozen inputs after the terminal screen" \
    "${PYTHON}" "${SCRIPT_ROOT}/audit_frozen_input_immutability.py" \
    --manifest "${MANIFEST}" \
    --output "${FROZEN_POST_SCREEN}"
require_json_pass "${FROZEN_POST_SCREEN}"

run_step "Audit every stable primary artifact before downstream consumers" \
    "${PYTHON}" "${SCRIPT_ROOT}/audit_stable_primary_artifacts.py" \
    --site-results "${SCREEN_SITES}" \
    --stable-assignments "${SCREEN_ASSIGNMENTS}" \
    --output "${PRIMARY_PRE}" \
    --workers 40 \
    --chunk-size 64 \
    --max-pending-chunks 80 \
    --progress-every 1000
require_json_pass "${PRIMARY_PRE}"

run_step "Analyze tumor-REF methyl background controls" \
    "${PYTHON}" "${SCRIPT_ROOT}/analyze_all_ssnv_tumor_ref_controls.py" \
    --site-results "${SCREEN_SITES}" \
    --stable-assignments "${SCREEN_ASSIGNMENTS}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --output-dir "${TUMOR_REF}" \
    --workers 40 \
    --chunk-size 16 \
    --max-pending-chunks 80 \
    --progress-every 1000
require_json_pass "${TUMOR_REF}/run_manifest.json"

run_step "Analyze focal methyl-group associations with frozen LongPhase-S PASS sSNVs" \
    "${PYTHON}" "${SCRIPT_ROOT}/analyze_methyl_ssnv_cooccurrence.py" \
    --manifest "${MANIFEST}" \
    --assignments "${SCREEN_ASSIGNMENTS}" \
    --sites "${SCREEN_SITES}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --intersubmod-root "${EXTRACT_ROOT}" \
    --output-dir "${COOCCURRENCE}" \
    --workers 40 \
    --chunk-size 8 \
    --max-pending-chunks 80 \
    --permutations 999 \
    --top-markers 3 \
    --exact-state-space-ceiling 250000
require_json_pass "${COOCCURRENCE}/run_receipt.json"

run_step "Run strict multi-seed and column-null robustness for formal G2 candidates" \
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
require_json_pass "${STRICT_RECEIPT}"

run_step "Run matched-normal controls only for the formal G2 selection" \
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
    run_step "Analyze matched-normal controls and exact read-identity joins" \
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
run_step "Re-audit primary artifacts and downstream consumer chronology" "${POST_AUDIT_COMMAND[@]}"
require_json_pass "${PRIMARY_POST}"

run_step "Re-audit frozen inputs after every downstream analysis" \
    "${PYTHON}" "${SCRIPT_ROOT}/audit_frozen_input_immutability.py" \
    --manifest "${MANIFEST}" \
    --output "${FROZEN_POST_TERMINAL}"
require_json_pass "${FROZEN_POST_TERMINAL}"

run_step "Build fail-closed final all-sSNV machine dataset" \
    "${PYTHON}" "${SCRIPT_ROOT}/build_all_ssnv_final_report_dataset.py" \
    --manifest "${MANIFEST}" \
    --screen-dir "${SCREEN}" \
    --cooccurrence-dir "${COOCCURRENCE}" \
    --strict-dir "${STRICT}" \
    --tumor-ref-dir "${TUMOR_REF}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --primary-artifact-audit-post "${PRIMARY_POST}" \
    --matched-normal-dir "${MATCHED_FINAL}" \
    --cn-ccf-annotations "${CN_CCF}" \
    --output-dir "${FINAL_DATASET}"
require_json_pass "${FINAL_DATASET}/run_receipt.json"

run_step "Build complete Traditional Chinese Markdown report and canonical artifact" \
    "${PYTHON}" "${SCRIPT_ROOT}/build_all_ssnv_report_artifact.py" \
    --final-dataset "${FINAL_DATASET}/final_report_dataset.json" \
    --final-receipt "${FINAL_DATASET}/run_receipt.json" \
    --candidate-catalog "${FINAL_DATASET}/candidate_catalog.tsv" \
    --candidate-witness-pairs "${FINAL_DATASET}/candidate_witness_pairs.tsv" \
    --claim-contract "${TERMINAL_CLAIM_CONTRACT}" \
    --manifest "${MANIFEST}" \
    --screen-summary "${SCREEN_SUMMARY}" \
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
printf 'Command: %q %q --html %q --width 1440 --height 1000 --show-scrollbars --screenshot %q --wait-ms 3000\n' \
    "${QA_PYTHON}" \
    "${SCRIPT_ROOT}/qa_portable_report_layout.py" "${PORTABLE_HTML}" "${DESKTOP_SCREENSHOT}"
"${QA_PYTHON}" "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 1440 \
    --height 1000 \
    --show-scrollbars \
    --screenshot "${DESKTOP_SCREENSHOT}" \
    --wait-ms 3000 | tee "${DESKTOP_QA}"

printf '\n[%s] Mobile portable HTML layout QA\n' "$(date --iso-8601=seconds)"
printf 'Command: %q %q --html %q --width 390 --height 844 --show-scrollbars --screenshot %q --wait-ms 3000\n' \
    "${QA_PYTHON}" \
    "${SCRIPT_ROOT}/qa_portable_report_layout.py" "${PORTABLE_HTML}" "${MOBILE_SCREENSHOT}"
"${QA_PYTHON}" "${SCRIPT_ROOT}/qa_portable_report_layout.py" \
    --html "${PORTABLE_HTML}" \
    --width 390 \
    --height 844 \
    --show-scrollbars \
    --screenshot "${MOBILE_SCREENSHOT}" \
    --wait-ms 3000 | tee "${MOBILE_QA}"

jq -e '.pass == true and .overlapCount == 0' "${DESKTOP_QA}" >/dev/null
jq -e '.pass == true and .overlapCount == 0' "${MOBILE_QA}" >/dev/null

printf '\nTerminal validation chain PASS\n'
printf 'Final dataset: %s\n' "${FINAL_DATASET}/final_report_dataset.json"
printf 'Markdown report: %s\n' "${REPORT}/report.md"
printf 'Portable HTML: %s\n' "${PORTABLE_HTML}"
printf 'Delivery receipt: %s\n' "${PORTABLE_RECEIPT}"
