#!/usr/bin/env bash

set -Eeuo pipefail

readonly REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
readonly TOPIC_ROOT="${REPO_ROOT}/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly SCRIPT_ROOT="${TOPIC_ROOT}/scripts"
readonly RESULT_ROOT="${TOPIC_ROOT}/results"
readonly WORKSPACE_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
readonly PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"

readonly ANALYZER="${SCRIPT_ROOT}/analyze_methyl_ssnv_cooccurrence.py"
readonly MANIFEST="${RESULT_ROOT}/all_ssnv_input_manifest.json"
readonly SCREEN="${WORKSPACE_ROOT}/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
readonly ASSIGNMENTS="${SCREEN}/all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
readonly SITES="${SCREEN}/all_ssnv_site_results.tsv.gz"
readonly PRIMARY_PRE="${RESULT_ROOT}/stable_primary_artifact_audit.v1_pre_downstream.json"
readonly INTERSUBMOD_ROOT="${WORKSPACE_ROOT}/intersubmod_all_ssnv_v2_verification_fix"
readonly PREFLIGHT="${RESULT_ROOT}/cooccurrence_task_contract_preflight.v4_raw_identity_full_runtime.json"
readonly OUTPUT="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked"
readonly LOG="${WORKSPACE_ROOT}/methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked.log"

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
    "${MANIFEST}" \
    "${ASSIGNMENTS}" \
    "${SITES}" \
    "${PRIMARY_PRE}" \
    "${PREFLIGHT}"; do
    require_file "${path}"
done
[[ ! -e "${OUTPUT}" ]] || {
    printf 'Refusing to overwrite output: %s\n' "${OUTPUT}" >&2
    exit 1
}
[[ ! -e "${LOG}" ]] || {
    printf 'Refusing to overwrite log: %s\n' "${LOG}" >&2
    exit 1
}
jq -e '.pass == true' "${PRIMARY_PRE}" >/dev/null
jq -e \
    --arg analyzer "$(realpath "${ANALYZER}")" \
    --arg sha256 "$(sha256sum "${ANALYZER}" | awk '{print $1}')" \
    '
    .pass == true and
    .observed.task_count == 102842 and
    .observed.raw_bam_identity_recovery.aggregate.site_tasks == 102842 and
    .observed.raw_bam_identity_recovery.missing_projections == 0 and
    .observed.raw_bam_identity_recovery.conflicting_analysis_payload_projections == 0 and
    .code.source_identity_after.analyzer.path == $analyzer and
    .code.source_identity_after.analyzer.sha256 == $sha256
    ' "${PREFLIGHT}" >/dev/null
[[ "$(stat -c '%a' "${ANALYZER}")" == "444" ]] || {
    printf 'Analyzer is not source-locked mode 0444: %s\n' "${ANALYZER}" >&2
    exit 1
}

exec > >(tee "${LOG}") 2>&1

printf 'Input manifest: %s\n' "${MANIFEST}"
printf 'Input assignments: %s\n' "${ASSIGNMENTS}"
printf 'Input sites: %s\n' "${SITES}"
printf 'Input full raw-identity preflight: %s\n' "${PREFLIGHT}"
printf 'Output cooccurrence: %s\n' "${OUTPUT}"
printf 'Command:'
printf ' %q' \
    "${PYTHON}" "${ANALYZER}" \
    --manifest "${MANIFEST}" \
    --assignments "${ASSIGNMENTS}" \
    --sites "${SITES}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
    --intersubmod-root "${INTERSUBMOD_ROOT}" \
    --output-dir "${OUTPUT}" \
    --workers 40 \
    --chunk-size 8 \
    --max-pending-chunks 80 \
    --permutations 999 \
    --top-markers 3 \
    --exact-state-space-ceiling 250000
printf '\n'

"${PYTHON}" "${ANALYZER}" \
    --manifest "${MANIFEST}" \
    --assignments "${ASSIGNMENTS}" \
    --sites "${SITES}" \
    --primary-artifact-audit-pre "${PRIMARY_PRE}" \
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
    "${OUTPUT}/summary.json" \
    "${OUTPUT}/run_receipt.json"; do
    require_file "${path}"
done
jq -e '.pass == true and .counts.stable_sites == 102842' "${OUTPUT}/run_receipt.json" >/dev/null
jq -e \
    --slurpfile receipt "${OUTPUT}/run_receipt.json" \
    --slurpfile summary "${OUTPUT}/summary.json" \
    '
    .code.source_identity_after.analyzer == $receipt[0].code.analyzer and
    .observed.raw_bam_identity_recovery.site_weighted_audit_sha256 ==
        $summary[0].raw_bam_identity_recovery_audit.site_weighted_audit_sha256 and
    .observed.raw_bam_identity_recovery.aggregate.expected_projection_occurrences ==
        $summary[0].raw_bam_identity_recovery_audit.n_expected_projection_occurrences and
    .observed.raw_bam_identity_recovery.aggregate.duplicate_projection_occurrences_collapsed ==
        $summary[0].raw_bam_identity_recovery_audit.n_duplicate_projection_occurrences_collapsed and
    .observed.raw_bam_identity_recovery.aggregate.duplicate_extra_record_occurrences_collapsed ==
        $summary[0].raw_bam_identity_recovery_audit.n_duplicate_extra_record_occurrences_collapsed
    ' "${PREFLIGHT}" >/dev/null

printf 'Cooccurrence v6 PASS: %s\n' "${OUTPUT}"
