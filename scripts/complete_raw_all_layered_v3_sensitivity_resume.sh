#!/usr/bin/env bash
# Resume the fail-closed raw-all validation after a canonical run has already passed.
set -euo pipefail
export LC_ALL=C TZ=UTC PYTHONHASHSEED=0

REPO=/big7_disk/liaoyoyo2001/InterSubMod
METHOD="$REPO/docs/methodology/_assets/20260627_subclone_4axis_teaching"
PY=/bip7_disk/liaoyoyo2001/miniconda3/bin/python3
RUN_PARENT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds
RUNNER="$REPO/scripts/run_layered_v3.py"
SUMMARIZER="$METHOD/scripts/summarize_backbone_sensitivity.py"
BAM_AUDITOR="$REPO/scripts/audit_canonical_longphase_bam_immutability.py"
SUPERVISOR="$REPO/scripts/complete_raw_all_layered_v3_sensitivity_resume.sh"

CANONICAL_MANIFEST="$REPO/research/20260710_layered_reconstruction_v2/data/layered_input_manifest_v3_raw_all_lps_pass_v2.json"
SENSITIVITY_MANIFEST="$REPO/research/20260710_layered_reconstruction_v2/data/layered_input_manifest_v3_raw_all_clairs_pass_sensitivity_v2.json"
CANONICAL_RUN="$RUN_PARENT/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5"
SENSITIVITY_RUN_ID="${INTERSUBMOD_SENSITIVITY_RUN_ID:-20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6}"
SENSITIVITY_RUN="$RUN_PARENT/$SENSITIVITY_RUN_ID"
COMPARISON="${INTERSUBMOD_COMPARISON:-$REPO/research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.json}"
BASE_MANIFEST="$METHOD/data/layered_v2_input_manifest.json"
BASELINE="$REPO/research/20260710_layered_reconstruction_v2/audit_receipts/20260711_canonical_longphase_tagged_bam_baseline_v1.json"
POST_BAM_RECEIPT="${INTERSUBMOD_POST_BAM_RECEIPT:-$REPO/research/20260710_layered_reconstruction_v2/audit_receipts/20260713_canonical_longphase_tagged_bam_post_full_validation_v6.json}"
EXEC_ROOT="${INTERSUBMOD_EXEC_ROOT:-$REPO/research/20260710_layered_reconstruction_v2/execution/20260713_raw_all_to_layered_v3_sensitivity_resume_v16}"

FAILED_SENSITIVITY="$RUN_PARENT/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v5"
V15_EXEC="$REPO/research/20260710_layered_reconstruction_v2/execution/20260713_raw_all_to_layered_v3_completion_v15"
FIX_ROOT="$REPO/research/20260710_layered_reconstruction_v2/diagnostics/20260713_clairs_pass_raw_all_ledger_contract_v1"
FIX_RESOLUTION="$FIX_ROOT/clairs_pass_raw_all_ledger_resolution.json"
FIX_TEST_LOG="$FIX_ROOT/raw_all_contract_tests.103.log"
FIX_REPLAY_SUMMARY="$FIX_ROOT/real_data_replay/COLO829/ssnv_site_ledger_COLO829.summary.json"
FIX_REPLAY_LEDGER="$FIX_ROOT/real_data_replay/COLO829/ssnv_site_ledger_COLO829.tsv.gz"
FIX_REPLAY_INDEX="$FIX_REPLAY_LEDGER.tbi"
FAILED_STAGE_SUMMARY="$FAILED_SENSITIVITY/samples/COLO829/.site_ledger.stage.2640400/ssnv_site_ledger_COLO829.summary.json"

for path in "$SENSITIVITY_RUN" "$COMPARISON" "${COMPARISON%.json}.tsv" "$POST_BAM_RECEIPT" "$EXEC_ROOT"; do
    if [[ -e "$path" || -L "$path" ]]; then
        printf 'E_OUTPUT_EXISTS: refusing to overwrite %s\n' "$path" >&2
        exit 7
    fi
done

REQUIRED_INPUTS=(
    "$CANONICAL_MANIFEST"
    "$SENSITIVITY_MANIFEST"
    "$CANONICAL_RUN/_SUCCESS"
    "$CANONICAL_RUN/verification_summary.json"
    "$CANONICAL_RUN/verification_summary.tsv"
    "$CANONICAL_RUN/verification_summary.sha256"
    "$BASE_MANIFEST"
    "$BASELINE"
    "$FAILED_SENSITIVITY/run_state.json"
    "$FAILED_SENSITIVITY/worker_COLO829.log"
    "$FAILED_STAGE_SUMMARY"
    "$V15_EXEC/_FAILED"
    "$V15_EXEC/execution.log"
    "$FIX_RESOLUTION"
    "$FIX_TEST_LOG"
    "$FIX_REPLAY_SUMMARY"
    "$FIX_REPLAY_LEDGER"
    "$FIX_REPLAY_INDEX"
)
for path in "${REQUIRED_INPUTS[@]}"; do
    if [[ ! -f "$path" || -L "$path" ]]; then
        printf 'E_REQUIRED_INPUT: expected regular file %s\n' "$path" >&2
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
printf 'INPUT canonical_run=%s\n' "$CANONICAL_RUN"
printf 'INPUT failed_sensitivity=%s\n' "$FAILED_SENSITIVITY"
printf 'INPUT sensitivity_manifest=%s\n' "$SENSITIVITY_MANIFEST"
printf 'OUTPUT sensitivity_run=%s\n' "$SENSITIVITY_RUN"
printf 'OUTPUT comparison=%s\n' "$COMPARISON"
printf 'OUTPUT post_bam_receipt=%s\n' "$POST_BAM_RECEIPT"

SOURCE_PINS=(
    "$PY"
    "$SUPERVISOR"
    "$RUNNER"
    "$REPO/scripts/layered_v3_lifecycle.py"
    "$REPO/scripts/verify_layered_v3.py"
    "$REPO/scripts/test_verify_layered_v3.py"
    "$METHOD/scripts/validate_layered_v3_inputs.py"
    "$METHOD/scripts/sm_linkage_genomewide.py"
    "$METHOD/scripts/sm_multilocus_combinations.py"
    "$METHOD/scripts/tree_enumeration_solver.py"
    "$METHOD/scripts/layered_tree_reconstruction.py"
    "$METHOD/scripts/build_region_view.py"
    "$METHOD/scripts/build_ssnv_site_ledger.py"
    "$METHOD/scripts/test_ssnv_site_ledger.py"
    "$REPO/scripts/run_layered_v3_raw_all_contract_tests.sh"
    "$SUMMARIZER"
    "$BAM_AUDITOR"
)
sha256sum "${SOURCE_PINS[@]}" > "$EXEC_ROOT/source_authorities.sha256"
sha256sum "${REQUIRED_INPUTS[@]}" > "$EXEC_ROOT/input_authorities.sha256"

verify_pins() {
    sha256sum -c "$EXEC_ROOT/source_authorities.sha256"
    sha256sum -c "$EXEC_ROOT/input_authorities.sha256"
}

verify_pins
(cd "$CANONICAL_RUN" && sha256sum -c verification_summary.sha256)
"$PY" -c '
import json, pathlib, sys
summary = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
ok = summary.get("all_pass") is True and summary.get("n_pass") == 7 and summary.get("n_fail") == 0
raise SystemExit(0 if ok else 7)
' "$CANONICAL_RUN/verification_summary.json"
grep -Fx 'LAYERED V3 RAW-ALL CONTRACT TESTS: 103/103 PASS' "$FIX_TEST_LOG"
"$PY" -c '
import json, pathlib, sys
summary = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
checks = summary.get("checks", {})
ok = (summary.get("pass") is True
      and summary.get("tree_contract") == "clairs_PASS_input"
      and summary.get("longphase_input_contract") == "clairs_raw_all"
      and checks.get("raw_all_equals_longphase_input") is True
      and checks.get("caller_raw_PASS_equals_tree_input") is True
      and all(value is True for value in checks.values()))
raise SystemExit(0 if ok else 7)
' "$FIX_REPLAY_SUMMARY"

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

printf 'STEP ClairS PASS sensitivity tree %s\n' "$(date --iso-8601=seconds)"
verify_pins
"$PY" "$RUNNER" \
    --manifest "$SENSITIVITY_MANIFEST" \
    --run-parent "$RUN_PARENT" \
    --run-id "$SENSITIVITY_RUN_ID" \
    "${RUNNER_OPTIONS[@]}"
[[ -f "$SENSITIVITY_RUN/_SUCCESS" ]]

printf 'STEP compare tree backbones %s\n' "$(date --iso-8601=seconds)"
verify_pins
"$PY" "$SUMMARIZER" \
    --base-run "$CANONICAL_RUN" \
    --clairs-run "$SENSITIVITY_RUN" \
    --input-manifest "$CANONICAL_MANIFEST" \
    --output "$COMPARISON"

printf 'STEP verify canonical tagged BAMs unchanged %s\n' "$(date --iso-8601=seconds)"
verify_pins
"$PY" "$BAM_AUDITOR" \
    --manifest "$BASE_MANIFEST" \
    --baseline "$BASELINE" \
    --output "$POST_BAM_RECEIPT"

sha256sum \
    "$CANONICAL_RUN/_SUCCESS" "$CANONICAL_RUN/verification_summary.json" \
    "$SENSITIVITY_RUN/_SUCCESS" "$SENSITIVITY_RUN/verification_summary.json" \
    "$COMPARISON" "${COMPARISON%.json}.tsv" "$POST_BAM_RECEIPT" \
    "$FIX_RESOLUTION" "$FIX_TEST_LOG" "$FIX_REPLAY_SUMMARY" \
    > "$EXEC_ROOT/artifacts.sha256"
printf 'status=SUCCESS\ntimestamp=%s\nartifacts_sha256=%s\n' \
    "$(date --iso-8601=seconds)" "$(sha256sum "$EXEC_ROOT/artifacts.sha256" | cut -d' ' -f1)" \
    > "$EXEC_ROOT/_SUCCESS"
COMPLETED=1
printf 'COMPLETE %s\n' "$(date --iso-8601=seconds)"
