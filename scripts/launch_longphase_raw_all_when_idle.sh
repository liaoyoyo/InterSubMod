#!/usr/bin/env bash
# Resource-gated, hash-pinned launcher for the seven-dataset raw-all LongPhase-S producer.
set -euo pipefail

REPO=/big7_disk/liaoyoyo2001/InterSubMod
METHOD="$REPO/docs/methodology/_assets/20260627_subclone_4axis_teaching"
SCRIPTS="$METHOD/scripts"
PY=/bip7_disk/liaoyoyo2001/miniconda3/bin/python3
MANIFEST="$METHOD/data/longphase_raw_all_production_manifest_v2.json"
LPS=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_zero_read_patch_probe_v1/source_snapshot/longphase-s
PATCH_RECEIPT="$REPO/research/20260710_layered_reconstruction_v2/longphase_s_zero_read_patch_build_receipt.json"
PATCH_FILE="$REPO/research/20260710_layered_reconstruction_v2/longphase_s_zero_read_lowqual_v1.patch"
PRODUCER="$SCRIPTS/run_longphase_raw_all_production_sidecars.sh"
RUN_ROOT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2
GATE_ROOT="$REPO/research/20260710_layered_reconstruction_v2/execution/20260711_raw_all_producer_resource_gate_v1"

if [[ -e "$RUN_ROOT" || -L "$RUN_ROOT" ]]; then
    printf 'E_OUTPUT_EXISTS: %s\n' "$RUN_ROOT" >&2
    exit 7
fi
if [[ -e "$GATE_ROOT" || -L "$GATE_ROOT" ]]; then
    printf 'E_GATE_EXISTS: %s\n' "$GATE_ROOT" >&2
    exit 7
fi
mkdir -p "$GATE_ROOT"
exec > >(tee "$GATE_ROOT/execution.log") 2>&1

COMPLETED=0
on_exit() {
    local rc=$?
    trap - EXIT
    if [[ "$rc" -ne 0 && "$COMPLETED" -eq 0 && ! -e "$GATE_ROOT/_FAILED" ]]; then
        printf 'status=FAILED\nexit_code=%s\ntimestamp=%s\n' \
            "$rc" "$(date --iso-8601=seconds)" > "$GATE_ROOT/_FAILED"
    fi
    exit "$rc"
}
trap on_exit EXIT
trap 'exit 130' INT TERM

PINS=(
    "$PY"
    "$REPO/scripts/launch_longphase_raw_all_when_idle.sh"
    "$PRODUCER"
    "$SCRIPTS/capture_longphase_tagged_bam_sidecar.py"
    "$SCRIPTS/validate_streamed_longphase_sidecar.py"
    "$SCRIPTS/audit_longphase_filter_transitions.py"
    "$SCRIPTS/build_longphase_raw_all_capture_receipt_v2.py"
    "$SCRIPTS/finalize_longphase_raw_all_capture_receipts.py"
    "$SCRIPTS/build_longphase_production_capture_receipt_v3.py"
    "$SCRIPTS/validate_layered_v3_inputs.py"
    "$METHOD/schemas/longphase_raw_all_capture_receipt_v2.schema.json"
    "$MANIFEST"
    "$LPS"
    "$PATCH_RECEIPT"
    "$PATCH_FILE"
    "$(command -v jq)"
    "$(command -v samtools)"
    "$(command -v bcftools)"
    "$(command -v bgzip)"
    "$(command -v tabix)"
)
sha256sum "${PINS[@]}" > "$GATE_ROOT/launch_authorities.sha256"

printf 'START %s\n' "$(date --iso-8601=seconds)"
printf 'INPUT manifest=%s\n' "$MANIFEST"
printf 'INPUT longphase_s=%s\n' "$LPS"
printf 'OUTPUT run_root=%s\n' "$RUN_ROOT"

while true; do
    sha256sum -c "$GATE_ROOT/launch_authorities.sha256" >/dev/null
    if [[ -e "$RUN_ROOT" || -L "$RUN_ROOT" ]]; then
        printf 'E_OUTPUT_APPEARED: %s\n' "$RUN_ROOT" >&2
        exit 7
    fi
    read -r LOAD1 _ < /proc/loadavg
    IOWAIT="$(vmstat 1 2 | awk 'END {print $16}')"
    printf '%s RESOURCE_GATE load1=%s iowait=%s%%\n' \
        "$(date --iso-8601=seconds)" "$LOAD1" "$IOWAIT"
    if awk -v l="$LOAD1" -v w="$IOWAIT" 'BEGIN {exit !(l < 60 && w < 25)}'; then
        break
    fi
    sleep 300
done

sha256sum -c "$GATE_ROOT/launch_authorities.sha256"
printf 'LAUNCH %s\n' "$(date --iso-8601=seconds)"
set +e
env \
    SM_RUN_ROOT="$RUN_ROOT" \
    SM_PARALLEL_SAMPLES=1 \
    SM_LPS_THREADS=12 \
    SM_LPS_MANIFEST="$MANIFEST" \
    SM_LONGPHASE_S="$LPS" \
    bash "$PRODUCER"
RC=$?
set -e
if [[ "$RC" -ne 0 ]]; then
    printf 'E_PRODUCER_EXIT: %s\n' "$RC" >&2
    exit "$RC"
fi
[[ -f "$RUN_ROOT/_SUCCESS" && ! -e "$RUN_ROOT/_FAILED" ]]
printf 'status=SUCCESS\ntimestamp=%s\nproducer_success_sha256=%s\n' \
    "$(date --iso-8601=seconds)" "$(sha256sum "$RUN_ROOT/_SUCCESS" | cut -d' ' -f1)" \
    > "$GATE_ROOT/_SUCCESS"
COMPLETED=1
printf 'COMPLETE %s\n' "$(date --iso-8601=seconds)"
