#!/usr/bin/env bash
# One-shot fail-closed continuation from raw-all LongPhase-S production to both layered-v3 tree backbones.
set -euo pipefail
export LC_ALL=C TZ=UTC PYTHONHASHSEED=0

REPO=/big7_disk/liaoyoyo2001/InterSubMod
METHOD="$REPO/docs/methodology/_assets/20260627_subclone_4axis_teaching"
PY=/bip7_disk/liaoyoyo2001/miniconda3/bin/python3
PRODUCTION=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2
RAW_RECEIPT_CLOSEOUT="$PRODUCTION/raw_all_receipt_closeout.json"
RAW_RECEIPT_MARKER="$PRODUCTION/_RAW_ALL_RECEIPTS_SUCCESS"
RUN_PARENT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds
BASE_MANIFEST="$METHOD/data/layered_v2_input_manifest.json"
REAL_RECEIPT="$REPO/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/equivalence_probe.json"
SYNTHETIC_RECEIPT_V1="$REPO/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/synthetic_contract_receipt.json"
SYNTHETIC_RECEIPT_V2="$REPO/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/synthetic_contract_receipt_v2.json"
SYNTHETIC_RECEIPT="$REPO/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/synthetic_contract_receipt_v3.json"
PREPARER="$METHOD/scripts/prepare_clean_layered_manifest_v3.py"
RECEIPT_FINALIZER="$METHOD/scripts/finalize_longphase_raw_all_capture_receipts.py"
RECEIPT_BUILDER="$METHOD/scripts/build_longphase_raw_all_capture_receipt_v2.py"
RUNNER="$REPO/scripts/run_layered_v3.py"
SUMMARIZER="$METHOD/scripts/summarize_backbone_sensitivity.py"
BAM_AUDITOR="$REPO/scripts/audit_canonical_longphase_bam_immutability.py"
SUPERVISOR="$REPO/scripts/complete_raw_all_layered_v3_validation.sh"

CANONICAL_MANIFEST="$REPO/research/20260710_layered_reconstruction_v2/data/layered_input_manifest_v3_raw_all_lps_pass_v2.json"
CANONICAL_FAILURE="$REPO/research/20260710_layered_reconstruction_v2/data/layered_input_manifest_v3_raw_all_lps_pass_v2.failure.json"
SENSITIVITY_MANIFEST="$REPO/research/20260710_layered_reconstruction_v2/data/layered_input_manifest_v3_raw_all_clairs_pass_sensitivity_v2.json"
SENSITIVITY_FAILURE="$REPO/research/20260710_layered_reconstruction_v2/data/layered_input_manifest_v3_raw_all_clairs_pass_sensitivity_v2.failure.json"
CANONICAL_RUN_ID="${INTERSUBMOD_CANONICAL_RUN_ID:-20260713_layered_reconstruction_v3_raw_all_lps_pass_v5}"
SENSITIVITY_RUN_ID="${INTERSUBMOD_SENSITIVITY_RUN_ID:-20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v5}"
CANONICAL_RUN="$RUN_PARENT/$CANONICAL_RUN_ID"
SENSITIVITY_RUN="$RUN_PARENT/$SENSITIVITY_RUN_ID"
COMPARISON="${INTERSUBMOD_COMPARISON:-$REPO/research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v5.json}"
BASELINE="$REPO/research/20260710_layered_reconstruction_v2/audit_receipts/20260711_canonical_longphase_tagged_bam_baseline_v1.json"
POST_BAM_RECEIPT="${INTERSUBMOD_POST_BAM_RECEIPT:-$REPO/research/20260710_layered_reconstruction_v2/audit_receipts/20260713_canonical_longphase_tagged_bam_post_full_validation_v5.json}"
RETRY_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_receipt_filter_audit_path_role_v1"
RETRY_RESOLUTION="$RETRY_ROOT/receipt_retry_resolution.json"
ARCHIVED_FAILURE="$RETRY_ROOT/raw_all_receipt_closeout.failure.attempt1.json"
FRESH_AUDIT="$RETRY_ROOT/HCC1395.fresh_packaged_filter_transition_audit.json"
PATH_COMPARISON="$RETRY_ROOT/HCC1395.audit_path_role_comparison.json"
CONTRACT_TEST_LOG="$RETRY_ROOT/raw_all_contract_tests.log"
REDUNDANCY_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_redundant_identity_contract_v1"
REDUNDANCY_RESOLUTION="$REDUNDANCY_ROOT/redundant_identity_resolution.json"
REDUNDANCY_FAILURE="$REDUNDANCY_ROOT/raw_all_receipt_closeout.failure.attempt2.json"
REDUNDANCY_TEST_LOG="$REDUNDANCY_ROOT/raw_all_contract_tests.89.log"
REDUNDANCY_TEST_SOURCES="$REDUNDANCY_ROOT/contract_test_sources.sha256"
REFERENCE_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_reference_fai_contract_v1"
REFERENCE_RESOLUTION="$REFERENCE_ROOT/reference_fai_resolution.json"
REFERENCE_FAILURE="$REFERENCE_ROOT/raw_all_receipt_closeout.failure.attempt3.json"
REFERENCE_TEST_LOG="$REFERENCE_ROOT/raw_all_contract_tests.90.log"
REFERENCE_TEST_SOURCES="$REFERENCE_ROOT/contract_test_sources.sha256"
SINGLE_SCAN_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_manifest_single_scan_v1"
SINGLE_SCAN_RESOLUTION="$SINGLE_SCAN_ROOT/manifest_single_scan_resolution.json"
SINGLE_SCAN_TEST_LOG="$SINGLE_SCAN_ROOT/raw_all_contract_tests.91.log"
SINGLE_SCAN_TEST_SOURCES="$SINGLE_SCAN_ROOT/contract_test_sources.sha256"
SINGLE_SCAN_FAILED_MARKER="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v4/_FAILED"
SINGLE_SCAN_EXECUTION_LOG="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v4/execution.log"
FAST_RECEIPT_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_receipt_restart_fast_verify_v1"
FAST_RECEIPT_RESOLUTION="$FAST_RECEIPT_ROOT/receipt_restart_fast_verify_resolution.json"
FAST_RECEIPT_TEST_LOG="$FAST_RECEIPT_ROOT/raw_all_contract_tests.93.log"
FAST_RECEIPT_TEST_SOURCES="$FAST_RECEIPT_ROOT/contract_test_sources.sha256"
SYNTHETIC_V2_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_synthetic_join_receipt_v2"
SYNTHETIC_V2_RESOLUTION="$SYNTHETIC_V2_ROOT/synthetic_join_receipt_v2_resolution.json"
SYNTHETIC_V2_FAILURE="$SYNTHETIC_V2_ROOT/layered_input_manifest_v3_raw_all_lps_pass_v1.failure.attempt1.json"
SYNTHETIC_V2_JOIN_LOG="$SYNTHETIC_V2_ROOT/test_read_tag_sidecar_contract.log"
SYNTHETIC_V2_TEST_LOG="$SYNTHETIC_V2_ROOT/raw_all_contract_tests.93.log"
SYNTHETIC_V2_TEST_SOURCES="$SYNTHETIC_V2_ROOT/contract_test_sources.sha256"
SYNTHETIC_V2_FAILED_MARKER="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v5/_FAILED"
SYNTHETIC_V2_EXECUTION_LOG="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v5/execution.log"
DETERMINISTIC_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_deterministic_environment_v1"
DETERMINISTIC_RESOLUTION="$DETERMINISTIC_ROOT/deterministic_environment_resolution.json"
DETERMINISTIC_TEST_LOG="$DETERMINISTIC_ROOT/raw_all_contract_tests.94.log"
DETERMINISTIC_TEST_SOURCES="$DETERMINISTIC_ROOT/contract_test_sources.sha256"
DETERMINISTIC_FAILED_MARKER="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v6/_FAILED"
DETERMINISTIC_EXECUTION_LOG="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v6/execution.log"
BUNDLE_RESOLUTION_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_bundle_relative_receipt_resolution_v1"
BUNDLE_RESOLUTION="$BUNDLE_RESOLUTION_ROOT/bundle_relative_receipt_resolution.json"
BUNDLE_FAILURE="$BUNDLE_RESOLUTION_ROOT/validation_failure.attempt1.json"
BUNDLE_MANIFEST="$BUNDLE_RESOLUTION_ROOT/source_bundle_manifest.attempt1.json"
BUNDLE_VALIDATOR_LOG="$BUNDLE_RESOLUTION_ROOT/validator.attempt1.log"
BUNDLE_TEST_LOG="$BUNDLE_RESOLUTION_ROOT/raw_all_contract_tests.96.log"
BUNDLE_TEST_SOURCES="$BUNDLE_RESOLUTION_ROOT/contract_test_sources.sha256"
BUNDLE_FAILED_MARKER="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v7/_FAILED"
BUNDLE_EXECUTION_LOG="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v7/execution.log"
SYSTEMEXIT_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_bundled_systemexit_contract_v1"
SYSTEMEXIT_RESOLUTION="$SYSTEMEXIT_ROOT/bundled_systemexit_resolution.json"
SYSTEMEXIT_TEST_LOG="$SYSTEMEXIT_ROOT/raw_all_contract_tests.97.log"
SYSTEMEXIT_TEST_SOURCES="$SYSTEMEXIT_ROOT/contract_test_sources.sha256"
SYSTEMEXIT_REPLAY_LEDGER="$SYSTEMEXIT_ROOT/real_data_replay/COLO829/ssnv_site_ledger_COLO829.tsv.gz"
SYSTEMEXIT_REPLAY_INDEX="$SYSTEMEXIT_REPLAY_LEDGER.tbi"
SYSTEMEXIT_REPLAY_SUMMARY="$SYSTEMEXIT_ROOT/real_data_replay/COLO829/ssnv_site_ledger_COLO829.summary.json"
SYSTEMEXIT_V8_FAILED="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v8/_FAILED"
SYSTEMEXIT_V8_EXECUTION_LOG="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v8/execution.log"
SYSTEMEXIT_V8_RUN_STATE="$RUN_PARENT/20260711_layered_reconstruction_v3_raw_all_lps_pass_v1/run_state.json"
SYSTEMEXIT_V8_WORKER_LOG="$RUN_PARENT/20260711_layered_reconstruction_v3_raw_all_lps_pass_v1/worker_COLO829.log"
SYSTEMEXIT_V8_SITE_LEDGER_LOG="$RUN_PARENT/20260711_layered_reconstruction_v3_raw_all_lps_pass_v1/samples/COLO829/site_ledger.log"
VERIFIER_FIX_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_verifier_logical_path_and_zero_funnel_v1"
VERIFIER_FIX_RESOLUTION="$VERIFIER_FIX_ROOT/verifier_failure_resolution.json"
VERIFIER_FIX_TEST_LOG="$VERIFIER_FIX_ROOT/raw_all_contract_tests.99.log"
VERIFIER_FIX_TEST_SOURCES="$VERIFIER_FIX_ROOT/contract_test_sources.sha256"
VERIFIER_V9_FAILED="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v9/_FAILED"
VERIFIER_V9_EXECUTION_LOG="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v9/execution.log"
VERIFIER_V2_RUN="$RUN_PARENT/20260712_layered_reconstruction_v3_raw_all_lps_pass_v2"
VERIFIER_V2_RUN_STATE="$VERIFIER_V2_RUN/run_state.json"
VERIFIER_V2_LOG="$VERIFIER_V2_RUN/verifier.log"
VERIFIER_V2_SUMMARY="$VERIFIER_V2_RUN/verification_summary.json"
VERIFIER_V2_TSV="$VERIFIER_V2_RUN/verification_summary.tsv"
VERIFIER_V2_HASHES="$VERIFIER_V2_RUN/verification_summary.sha256"
V11_EXECUTION_ROOT="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v11"
V11_FAILED_MARKER="$V11_EXECUTION_ROOT/_FAILED"
V11_EXECUTION_LOG="$V11_EXECUTION_ROOT/execution.log"
V11_VALIDATION_FAILURE="$RUN_PARENT/.20260712_layered_reconstruction_v3_raw_all_lps_pass_v3.staging.2138181.362b6425258b/validator_runtime/validation_failure.json"
V12_EXECUTION_ROOT="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v12"
V12_FAILED_MARKER="$V12_EXECUTION_ROOT/_FAILED"
V12_EXECUTION_LOG="$V12_EXECUTION_ROOT/execution.log"
V12_CANONICAL_RUN="$RUN_PARENT/20260712_layered_reconstruction_v3_raw_all_lps_pass_v3"
V12_RUN_STATE="$V12_CANONICAL_RUN/run_state.json"
V12_VERIFIER_LOG="$V12_CANONICAL_RUN/verifier.log"
V12_VERIFICATION_SUMMARY="$V12_CANONICAL_RUN/verification_summary.json"
V12_VERIFICATION_TSV="$V12_CANONICAL_RUN/verification_summary.tsv"
V12_VERIFICATION_HASHES="$V12_CANONICAL_RUN/verification_summary.sha256"
V12_OUTPUT_MANIFESTS="$V12_CANONICAL_RUN/output_manifests.json"
V12_HCC1937_OUTPUT_MANIFEST="$V12_CANONICAL_RUN/samples/HCC1937/output_manifest.json"
V13_EXECUTION_ROOT="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v13"
V13_FAILED_MARKER="$V13_EXECUTION_ROOT/_FAILED"
V13_EXECUTION_LOG="$V13_EXECUTION_ROOT/execution.log"
V13_SOURCE_AUTHORITIES="$V13_EXECUTION_ROOT/source_authorities.sha256"
V13_INPUT_AUTHORITIES="$V13_EXECUTION_ROOT/input_authorities.sha256"
V14_EXECUTION_ROOT="$REPO/research/20260710_layered_reconstruction_v2/execution/20260712_raw_all_to_layered_v3_completion_v14"
V14_FAILED_MARKER="$V14_EXECUTION_ROOT/_FAILED"
V14_EXECUTION_LOG="$V14_EXECUTION_ROOT/execution.log"
V14_SOURCE_AUTHORITIES="$V14_EXECUTION_ROOT/source_authorities.sha256"
V14_INPUT_AUTHORITIES="$V14_EXECUTION_ROOT/input_authorities.sha256"
V14_CANONICAL_RUN="$RUN_PARENT/20260712_layered_reconstruction_v3_raw_all_lps_pass_v4"
V14_RUN_STATE="$V14_CANONICAL_RUN/run_state.json"
V14_VERIFIER_LOG="$V14_CANONICAL_RUN/verifier.log"
V14_VERIFICATION_SUMMARY="$V14_CANONICAL_RUN/verification_summary.json"
V14_VERIFICATION_TSV="$V14_CANONICAL_RUN/verification_summary.tsv"
V14_VERIFICATION_HASHES="$V14_CANONICAL_RUN/verification_summary.sha256"
V14_OUTPUT_MANIFESTS="$V14_CANONICAL_RUN/output_manifests.json"
V14_HCC1937_OUTPUT_MANIFEST="$V14_CANONICAL_RUN/samples/HCC1937/output_manifest.json"
REGION_FIX_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260712_zero_population_region_conservation_v1"
REGION_FIX_RESOLUTION="$REGION_FIX_ROOT/zero_population_region_resolution.json"
REGION_FIX_PARALLEL_RESOLUTION="$REGION_FIX_ROOT/parallel4_runner_contract_resolution.json"
REGION_FIX_VERIFIER_RESOLUTION="$REGION_FIX_ROOT/parallel4_verifier_contract_resolution.json"
REGION_FIX_TEST_LOG="$REGION_FIX_ROOT/raw_all_contract_tests.102.v15.authority_final.log"
REGION_FIX_TEST_SOURCES="$REGION_FIX_ROOT/test_sources.sha256"
EXEC_ROOT="${INTERSUBMOD_EXEC_ROOT:-$REPO/research/20260710_layered_reconstruction_v2/execution/20260713_raw_all_to_layered_v3_completion_v15}"

for path in "$CANONICAL_FAILURE" "$SENSITIVITY_FAILURE"; do
    if [[ -e "$path" || -L "$path" ]]; then
        printf 'E_PRIOR_MANIFEST_FAILURE: refusing to continue with %s\n' "$path" >&2
        exit 7
    fi
done

USE_EXISTING_MANIFESTS=0
if [[ -f "$CANONICAL_MANIFEST" && ! -L "$CANONICAL_MANIFEST" \
   && -f "$SENSITIVITY_MANIFEST" && ! -L "$SENSITIVITY_MANIFEST" ]]; then
    USE_EXISTING_MANIFESTS=1
elif [[ -e "$CANONICAL_MANIFEST" || -L "$CANONICAL_MANIFEST" \
     || -e "$SENSITIVITY_MANIFEST" || -L "$SENSITIVITY_MANIFEST" ]]; then
    printf 'E_MANIFEST_PAIR: canonical and sensitivity manifests must both be new or both be regular files\n' >&2
    exit 7
fi

for path in \
    "$CANONICAL_RUN" "$SENSITIVITY_RUN" \
    "$COMPARISON" "${COMPARISON%.json}.tsv" "$POST_BAM_RECEIPT" "$EXEC_ROOT"; do
    if [[ -e "$path" || -L "$path" ]]; then
        printf 'E_OUTPUT_EXISTS: refusing to overwrite %s\n' "$path" >&2
        exit 7
    fi
done

mkdir -p "$EXEC_ROOT"
LOG="$EXEC_ROOT/execution.log"
exec > >(tee "$LOG") 2>&1

COMPLETED=0
on_exit() {
    local rc=$?
    trap - EXIT
    if [[ "$rc" -ne 0 && "$COMPLETED" -eq 0 && ! -e "$EXEC_ROOT/_FAILED" ]]; then
        printf 'status=FAILED\nexit_code=%s\ntimestamp=%s\n' \
            "$rc" "$(date --iso-8601=seconds)" > "$EXEC_ROOT/_FAILED"
    fi
    exit "$rc"
}
trap on_exit EXIT
trap 'exit 130' INT TERM

printf 'START %s\n' "$(date --iso-8601=seconds)"
printf 'INPUT producer=%s\n' "$PRODUCTION"
printf 'INPUT base_manifest=%s\n' "$BASE_MANIFEST"
printf 'OUTPUT canonical_run=%s\n' "$CANONICAL_RUN"
printf 'OUTPUT sensitivity_run=%s\n' "$SENSITIVITY_RUN"
SOURCE_PINS=(
    "$PY"
    "$SUPERVISOR"
    "$RECEIPT_FINALIZER"
    "$RECEIPT_BUILDER"
    "$PREPARER"
    "$METHOD/scripts/validate_layered_v3_inputs.py"
    "$METHOD/schemas/layered_input_manifest_v3.schema.json"
    "$METHOD/schemas/layered_input_lock_v1.schema.json"
    "$METHOD/schemas/longphase_raw_all_capture_receipt_v2.schema.json"
    "$RUNNER"
    "$REPO/scripts/layered_v3_lifecycle.py"
    "$REPO/scripts/verify_layered_v3.py"
    "$REPO/scripts/test_verify_layered_v3.py"
    "$METHOD/scripts/sm_linkage_genomewide.py"
    "$METHOD/scripts/sm_multilocus_combinations.py"
    "$METHOD/scripts/tree_enumeration_solver.py"
    "$METHOD/scripts/layered_tree_reconstruction.py"
    "$METHOD/scripts/build_region_view.py"
    "$METHOD/scripts/test_build_region_view_contract.py"
    "$METHOD/scripts/build_ssnv_site_ledger.py"
    "$REPO/scripts/run_layered_v3_raw_all_contract_tests.sh"
    "$REPO/scripts/test_run_layered_v3.py"
    "$SUMMARIZER"
    "$BAM_AUDITOR"
)
INPUT_PINS=(
    "$BASE_MANIFEST"
    "$METHOD/data/longphase_raw_all_production_manifest_v2.json"
    "$REAL_RECEIPT"
    "$SYNTHETIC_RECEIPT_V1"
    "$SYNTHETIC_RECEIPT_V2"
    "$SYNTHETIC_RECEIPT"
    "$BASELINE"
    "$RETRY_RESOLUTION"
    "$ARCHIVED_FAILURE"
    "$FRESH_AUDIT"
    "$PATH_COMPARISON"
    "$CONTRACT_TEST_LOG"
    "$REDUNDANCY_RESOLUTION"
    "$REDUNDANCY_FAILURE"
    "$REDUNDANCY_TEST_LOG"
    "$REDUNDANCY_TEST_SOURCES"
    "$REFERENCE_RESOLUTION"
    "$REFERENCE_FAILURE"
    "$REFERENCE_TEST_LOG"
    "$REFERENCE_TEST_SOURCES"
    "$SINGLE_SCAN_RESOLUTION"
    "$SINGLE_SCAN_TEST_LOG"
    "$SINGLE_SCAN_TEST_SOURCES"
    "$SINGLE_SCAN_FAILED_MARKER"
    "$SINGLE_SCAN_EXECUTION_LOG"
    "$FAST_RECEIPT_RESOLUTION"
    "$FAST_RECEIPT_TEST_LOG"
    "$FAST_RECEIPT_TEST_SOURCES"
    "$SYNTHETIC_V2_RESOLUTION"
    "$SYNTHETIC_V2_FAILURE"
    "$SYNTHETIC_V2_JOIN_LOG"
    "$SYNTHETIC_V2_TEST_LOG"
    "$SYNTHETIC_V2_TEST_SOURCES"
    "$SYNTHETIC_V2_FAILED_MARKER"
    "$SYNTHETIC_V2_EXECUTION_LOG"
    "$DETERMINISTIC_RESOLUTION"
    "$DETERMINISTIC_TEST_LOG"
    "$DETERMINISTIC_TEST_SOURCES"
    "$DETERMINISTIC_FAILED_MARKER"
    "$DETERMINISTIC_EXECUTION_LOG"
    "$BUNDLE_RESOLUTION"
    "$BUNDLE_FAILURE"
    "$BUNDLE_MANIFEST"
    "$BUNDLE_VALIDATOR_LOG"
    "$BUNDLE_TEST_LOG"
    "$BUNDLE_TEST_SOURCES"
    "$BUNDLE_FAILED_MARKER"
    "$BUNDLE_EXECUTION_LOG"
    "$SYSTEMEXIT_RESOLUTION"
    "$SYSTEMEXIT_TEST_LOG"
    "$SYSTEMEXIT_TEST_SOURCES"
    "$SYSTEMEXIT_REPLAY_LEDGER"
    "$SYSTEMEXIT_REPLAY_INDEX"
    "$SYSTEMEXIT_REPLAY_SUMMARY"
    "$SYSTEMEXIT_V8_FAILED"
    "$SYSTEMEXIT_V8_EXECUTION_LOG"
    "$SYSTEMEXIT_V8_RUN_STATE"
    "$SYSTEMEXIT_V8_WORKER_LOG"
    "$SYSTEMEXIT_V8_SITE_LEDGER_LOG"
    "$VERIFIER_FIX_RESOLUTION"
    "$VERIFIER_FIX_TEST_LOG"
    "$VERIFIER_FIX_TEST_SOURCES"
    "$VERIFIER_V9_FAILED"
    "$VERIFIER_V9_EXECUTION_LOG"
    "$VERIFIER_V2_RUN_STATE"
    "$VERIFIER_V2_LOG"
    "$VERIFIER_V2_SUMMARY"
    "$VERIFIER_V2_TSV"
    "$VERIFIER_V2_HASHES"
    "$V11_FAILED_MARKER"
    "$V11_EXECUTION_LOG"
    "$V11_VALIDATION_FAILURE"
    "$V12_FAILED_MARKER"
    "$V12_EXECUTION_LOG"
    "$V12_RUN_STATE"
    "$V12_VERIFIER_LOG"
    "$V12_VERIFICATION_SUMMARY"
    "$V12_VERIFICATION_TSV"
    "$V12_VERIFICATION_HASHES"
    "$V12_OUTPUT_MANIFESTS"
    "$V12_HCC1937_OUTPUT_MANIFEST"
    "$V13_FAILED_MARKER"
    "$V13_EXECUTION_LOG"
    "$V13_SOURCE_AUTHORITIES"
    "$V13_INPUT_AUTHORITIES"
    "$V14_FAILED_MARKER"
    "$V14_EXECUTION_LOG"
    "$V14_SOURCE_AUTHORITIES"
    "$V14_INPUT_AUTHORITIES"
    "$V14_RUN_STATE"
    "$V14_VERIFIER_LOG"
    "$V14_VERIFICATION_SUMMARY"
    "$V14_VERIFICATION_TSV"
    "$V14_VERIFICATION_HASHES"
    "$V14_OUTPUT_MANIFESTS"
    "$V14_HCC1937_OUTPUT_MANIFEST"
    "$REGION_FIX_RESOLUTION"
    "$REGION_FIX_PARALLEL_RESOLUTION"
    "$REGION_FIX_VERIFIER_RESOLUTION"
    "$REGION_FIX_TEST_LOG"
    "$REGION_FIX_TEST_SOURCES"
    "$RAW_RECEIPT_CLOSEOUT"
    "$RAW_RECEIPT_MARKER"
)
if [[ "$USE_EXISTING_MANIFESTS" -eq 1 ]]; then
    INPUT_PINS+=("$CANONICAL_MANIFEST" "$SENSITIVITY_MANIFEST")
fi
sha256sum "${SOURCE_PINS[@]}" > "$EXEC_ROOT/source_authorities.sha256"
sha256sum "${INPUT_PINS[@]}" > "$EXEC_ROOT/input_authorities.sha256"

verify_wait_pins() {
    sha256sum -c "$EXEC_ROOT/source_authorities.sha256"
    sha256sum -c "$EXEC_ROOT/input_authorities.sha256"
}

while true; do
    if [[ -e "$PRODUCTION/_FAILED" ]]; then
        printf 'E_PRODUCER_FAILED: %s\n' "$PRODUCTION/_FAILED" >&2
        exit 7
    fi
    if [[ -e "$PRODUCTION/_SUCCESS" ]]; then
        break
    fi
    printf 'WAIT producer terminal state %s\n' "$(date --iso-8601=seconds)"
    sleep 300
done

verify_wait_pins
printf 'STEP finalize raw-all receipts %s\n' "$(date --iso-8601=seconds)"
"$PY" "$RECEIPT_FINALIZER" \
    --run-root "$PRODUCTION"
[[ -f "$PRODUCTION/_RAW_ALL_RECEIPTS_SUCCESS" ]]

REGULAR_BAM="$(find "$PRODUCTION" -type f -name '*.bam' -print -quit)"
if [[ -n "$REGULAR_BAM" ]]; then
    printf 'E_PERSISTED_BAM: regular BAM found below %s\n' "$PRODUCTION" >&2
    exit 7
fi

printf 'STEP build canonical manifest %s\n' "$(date --iso-8601=seconds)"
verify_wait_pins
if [[ "$USE_EXISTING_MANIFESTS" -eq 1 ]]; then
    printf 'EXISTING canonical manifest pinned -> %s\n' "$CANONICAL_MANIFEST"
else
    "$PY" "$PREPARER" \
        --base-manifest "$BASE_MANIFEST" \
        --longphase-manifest "$PRODUCTION/input_manifest.json" \
        --production-root "$PRODUCTION" \
        --real-data-receipt "$REAL_RECEIPT" \
        --synthetic-receipt "$SYNTHETIC_RECEIPT" \
        --manifest-id 20260712_layered_v3_raw_all_lps_pass_v2 \
        --output "$CANONICAL_MANIFEST" \
        --failure-report "$CANONICAL_FAILURE"
fi

printf 'STEP build ClairS sensitivity manifest %s\n' "$(date --iso-8601=seconds)"
verify_wait_pins
if [[ "$USE_EXISTING_MANIFESTS" -eq 1 ]]; then
    printf 'EXISTING sensitivity manifest pinned -> %s\n' "$SENSITIVITY_MANIFEST"
else
    "$PY" "$PREPARER" \
        --base-manifest "$BASE_MANIFEST" \
        --longphase-manifest "$PRODUCTION/input_manifest.json" \
        --production-root "$PRODUCTION" \
        --real-data-receipt "$REAL_RECEIPT" \
        --synthetic-receipt "$SYNTHETIC_RECEIPT" \
        --manifest-id 20260712_layered_v3_raw_all_clairs_pass_sensitivity_v2 \
        --tree-input-contract clairs_FILTER_PASS_sensitivity \
        --output "$SENSITIVITY_MANIFEST" \
        --failure-report "$SENSITIVITY_FAILURE"
fi

RUNNER_OPTIONS=(
    --parallel-samples 4
    --verify-every 1
    --analysis-tree-cap 0
    --display-tree-cap 32
    --minread 3
    --max-snv 8
    --tier-r 50000
    --mapq-min 20
    --baseq-min 0
    --heartbeat-interval 60
    --min-logical-cpus 8
    --min-available-memory-gib 128
    --min-free-disk-gib 500
    --max-load-per-cpu 1.25
    --max-iowait-percent 20
    --resource-sample-seconds 300
    --max-nfs-read-mbps 80
    --nfs-mountpoint /big8_disk
)

printf 'STEP canonical LongPhase-S PASS tree %s\n' "$(date --iso-8601=seconds)"
verify_wait_pins
"$PY" "$RUNNER" \
    --manifest "$CANONICAL_MANIFEST" \
    --run-parent "$RUN_PARENT" \
    --run-id "$CANONICAL_RUN_ID" \
    "${RUNNER_OPTIONS[@]}"
[[ -f "$CANONICAL_RUN/_SUCCESS" ]]

printf 'STEP ClairS PASS sensitivity tree %s\n' "$(date --iso-8601=seconds)"
verify_wait_pins
"$PY" "$RUNNER" \
    --manifest "$SENSITIVITY_MANIFEST" \
    --run-parent "$RUN_PARENT" \
    --run-id "$SENSITIVITY_RUN_ID" \
    "${RUNNER_OPTIONS[@]}"
[[ -f "$SENSITIVITY_RUN/_SUCCESS" ]]

printf 'STEP compare tree backbones %s\n' "$(date --iso-8601=seconds)"
verify_wait_pins
"$PY" "$SUMMARIZER" \
    --base-run "$CANONICAL_RUN" \
    --clairs-run "$SENSITIVITY_RUN" \
    --input-manifest "$CANONICAL_MANIFEST" \
    --output "$COMPARISON"

printf 'STEP verify canonical tagged BAMs unchanged %s\n' "$(date --iso-8601=seconds)"
verify_wait_pins
"$PY" "$BAM_AUDITOR" \
    --manifest "$BASE_MANIFEST" \
    --baseline "$BASELINE" \
    --output "$POST_BAM_RECEIPT"

sha256sum \
    "$PRODUCTION/_SUCCESS" "$PRODUCTION/_RAW_ALL_RECEIPTS_SUCCESS" \
    "$CANONICAL_MANIFEST" "$SENSITIVITY_MANIFEST" \
    "$CANONICAL_RUN/_SUCCESS" "$CANONICAL_RUN/verification_summary.json" \
    "$SENSITIVITY_RUN/_SUCCESS" "$SENSITIVITY_RUN/verification_summary.json" \
    "$COMPARISON" "${COMPARISON%.json}.tsv" "$POST_BAM_RECEIPT" \
    > "$EXEC_ROOT/artifacts.sha256"
printf 'status=SUCCESS\ntimestamp=%s\nartifacts_sha256=%s\n' \
    "$(date --iso-8601=seconds)" "$(sha256sum "$EXEC_ROOT/artifacts.sha256" | cut -d' ' -f1)" \
    > "$EXEC_ROOT/_SUCCESS"
COMPLETED=1
printf 'COMPLETE %s\n' "$(date --iso-8601=seconds)"
