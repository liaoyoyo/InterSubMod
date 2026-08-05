#!/usr/bin/env bash
# Bounded raw-all LongPhase-S probe; preserves every partial artifact on failure.
set -euo pipefail

REPO="/big7_disk/liaoyoyo2001/InterSubMod"
ROOT="$REPO/research/20260710_layered_reconstruction_v2/probes/20260711_HCC1954_chr22_raw_all_longphase_contract_v1"
RUN="$ROOT/raw_all_patched_regression"
SCRIPTS="$REPO/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"
PY="/bip7_disk/liaoyoyo2001/miniconda3/bin/python3"
LPS="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_zero_read_patch_probe_v1/source_snapshot/longphase-s"
RAW="$REPO/research/20260710_layered_reconstruction_v2/probes/20260711_HCC1954_chr22_longphase_pass_contract_v1/raw_clairs.chr22.vcf.gz"
GERMLINE="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
NORMAL="/big8_disk/data/HCC1954/ONT/HCC1954BL.bam"
TUMOR="/big8_disk/data/HCC1954/ONT/HCC1954.bam"
REFERENCE="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
THREADS=12

[[ ! -e "$RUN" ]] || { echo "immutable probe output exists: $RUN" >&2; exit 2; }
mkdir -p "$RUN"

record_unhandled_failure() {
    local rc="$?"
    trap - EXIT
    if [[ "$rc" -ne 0 && ! -e "$RUN/_FAILED" && ! -e "$RUN/_SUCCESS" ]]; then
        jq -n --argjson exit_code "$rc" --arg reason "Unhandled raw-all probe stage failure" \
            '{status:"FAILED",stage:"unhandled",exit_code:$exit_code,reason:$reason}' > "$RUN/_FAILED"
    fi
    exit "$rc"
}
trap record_unhandled_failure EXIT

NORM="$RUN/HCC1954.raw_all.normalized.chr22.vcf.gz"
PREFIX="$RUN/HCC1954_raw_all_chr22"
SIDECAR="$RUN/HCC1954.raw_all.read_tags.tsv"
LOG="$RUN/longphase_s.raw_all.log"

zcat "$RAW" \
    | sed 's/ID=GQ,Number=1,Type=Integer/ID=GQ,Number=1,Type=Float/' \
    | bcftools view -Oz -o "$NORM"
bcftools index -c "$NORM"

printf '%q ' "$LPS" somatic_haplotag -s "$GERMLINE" -b "$NORMAL" \
    --tumor-snv-file "$NORM" --tumor-bam-file "$TUMOR" -r "$REFERENCE" \
    -t "$THREADS" --tagSupplementary -q 20 --region chr22 --output-somatic-vcf -o "$PREFIX" \
    > "$RUN/command.sh.txt"
printf '\n' >> "$RUN/command.sh.txt"
sha256sum "$RAW" "${RAW}.csi" "$NORM" "${NORM}.csi" "$GERMLINE" "${GERMLINE}.csi" \
    "${NORMAL}.bai" "${TUMOR}.bai" "$SCRIPTS/capture_longphase_tagged_bam_sidecar.py" \
    "$SCRIPTS/validate_streamed_longphase_sidecar.py" "$SCRIPTS/audit_longphase_filter_transitions.py" \
    "$ROOT/run_probe.sh" "$LPS" > "$RUN/input_and_code.sha256"

mkfifo "${PREFIX}.bam"
"$PY" "$SCRIPTS/capture_longphase_tagged_bam_sidecar.py" \
    --input-bam "${PREFIX}.bam" --output "$SIDECAR" --summary "$RUN/stream_capture_summary.json" \
    > "$RUN/stream_capture.log" 2>&1 &
capture_pid="$!"

set +e
/usr/bin/time -v "$LPS" somatic_haplotag -s "$GERMLINE" -b "$NORMAL" \
    --tumor-snv-file "$NORM" --tumor-bam-file "$TUMOR" -r "$REFERENCE" \
    -t "$THREADS" --tagSupplementary -q 20 --region chr22 --output-somatic-vcf -o "$PREFIX" \
    > "$LOG" 2>&1
lps_rc="$?"
set -e
if [[ "$lps_rc" -ne 0 ]]; then
    kill "$capture_pid" 2>/dev/null || true
    wait "$capture_pid" 2>/dev/null || true
    jq -n --argjson exit_code "$lps_rc" --arg reason "LongPhase-S raw-all probe failed" \
        '{status:"FAILED",stage:"longphase",exit_code:$exit_code,reason:$reason}' > "$RUN/_FAILED"
    exit "$lps_rc"
fi

set +e
wait "$capture_pid"
capture_rc="$?"
set -e
if [[ "$capture_rc" -ne 0 ]]; then
    jq -n --argjson exit_code "$capture_rc" --arg reason "FIFO sidecar capture failed" \
        '{status:"FAILED",stage:"capture",exit_code:$exit_code,reason:$reason}' > "$RUN/_FAILED"
    exit "$capture_rc"
fi
mv "${PREFIX}.bam" "$RUN/consumed_tagged_bam.fifo"

"$PY" "$SCRIPTS/validate_streamed_longphase_sidecar.py" \
    --sidecar "$SIDECAR" --capture-summary "$RUN/stream_capture_summary.json" \
    --execution-log "$LOG" --input-vcf "$NORM" --recalibrated-vcf "${PREFIX}_sc.vcf" \
    --region chr22 --output "$RUN/sidecar_validation.json"
"$PY" "$SCRIPTS/audit_longphase_filter_transitions.py" \
    --input-vcf "$NORM" --output-vcf "${PREFIX}_sc.vcf" \
    --output-json "$RUN/filter_transition_audit.json"

bgzip -@ "$THREADS" -c "$SIDECAR" > "$RUN/HCC1954.raw_all.read_tags.tsv.gz"
tabix -0 -s 1 -b 2 -e 3 "$RUN/HCC1954.raw_all.read_tags.tsv.gz"
bgzip -@ "$THREADS" -c "${PREFIX}_sc.vcf" > "$RUN/HCC1954.raw_all.recalibrated.all.vcf.gz"
bcftools index -c "$RUN/HCC1954.raw_all.recalibrated.all.vcf.gz"
bcftools view -f PASS -Oz -o "$RUN/HCC1954.raw_all.recalibrated.pass.vcf.gz" "${PREFIX}_sc.vcf"
bcftools index -c "$RUN/HCC1954.raw_all.recalibrated.pass.vcf.gz"

sha256sum "$RUN"/*.json "$RUN"/*.tsv "$RUN"/*.tsv.gz "$RUN"/*.tbi \
    "$RUN"/*.vcf "$RUN"/*.vcf.gz "$RUN"/*.csi "$RUN/command.sh.txt" "$RUN/input_and_code.sha256" \
    > "$RUN/artifacts.sha256"
jq -n --arg input "$NORM" --arg output "${PREFIX}_sc.vcf" \
    --argjson input_records "$(bcftools view -H "$NORM" | wc -l)" \
    --argjson output_records "$(bcftools view -H "${PREFIX}_sc.vcf" | wc -l)" \
    --argjson rescued "$(jq '.rescued_nonpass_to_pass' "$RUN/filter_transition_audit.json")" \
    --argjson removed "$(jq '.removed_pass_to_nonpass' "$RUN/filter_transition_audit.json")" \
    '{status:"SUCCESS",scope:"HCC1954_chr22_raw_all_probe",input:$input,output:$output,input_records:$input_records,output_records:$output_records,rescued:$rescued,removed:$removed}' \
    > "$RUN/_SUCCESS"
echo "RAW-ALL PROBE PASS -> $RUN"
