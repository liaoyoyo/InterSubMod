#!/usr/bin/env bash
# Reproduce the historical raw-ClairS LowQual crash locus with the current binary.
set -euo pipefail

REPO="/big7_disk/liaoyoyo2001/InterSubMod"
ROOT="$REPO/research/20260710_layered_reconstruction_v2/probes/20260711_HCC1395_chrX_72880028_raw_all_crash_probe_v1"
RUN="$ROOT/run_patched_v3"
SCRIPTS="$REPO/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"
PY="/bip7_disk/liaoyoyo2001/miniconda3/bin/python3"
LPS="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_zero_read_patch_probe_v1/source_snapshot/longphase-s"
RAW="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz"
GERMLINE="/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
NORMAL="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam"
TUMOR="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam"
REFERENCE="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
REGION="chrX:72880028-72880028"
THREADS=12

[[ ! -e "$RUN" ]] || { echo "immutable probe output exists: $RUN" >&2; exit 2; }
mkdir -p "$RUN"
record_failure() {
    local rc="$?"
    trap - EXIT
    if [[ "$rc" -ne 0 && ! -e "$RUN/_FAILED" && ! -e "$RUN/_SUCCESS" ]]; then
        jq -n --argjson exit_code "$rc" '{status:"FAILED",exit_code:$exit_code}' > "$RUN/_FAILED"
    fi
    exit "$rc"
}
trap record_failure EXIT

INPUT="$RUN/HCC1395.chrX_72880028.normalized.vcf.gz"
PREFIX="$RUN/HCC1395_chrX_72880028"
SIDECAR="$RUN/HCC1395.chrX_72880028.read_tags.tsv"
LOG="$RUN/longphase_s.crash_probe.log"

bcftools view -r "$REGION" "$RAW" \
    | sed 's/ID=GQ,Number=1,Type=Integer/ID=GQ,Number=1,Type=Float/' \
    | bcftools view -Oz -o "$INPUT"
bcftools index -c "$INPUT"
[[ "$(bcftools view -H "$INPUT" | wc -l)" -eq 1 ]]

printf '%q ' "$LPS" somatic_haplotag -s "$GERMLINE" -b "$NORMAL" \
    --tumor-snv-file "$INPUT" --tumor-bam-file "$TUMOR" -r "$REFERENCE" \
    -t "$THREADS" --tagSupplementary -q 20 --region "$REGION" --output-somatic-vcf -o "$PREFIX" \
    > "$RUN/command.sh.txt"
printf '\n' >> "$RUN/command.sh.txt"
sha256sum "$RAW" "${RAW}.tbi" "$INPUT" "${INPUT}.csi" "$GERMLINE" "${GERMLINE}.csi" \
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
    --tumor-snv-file "$INPUT" --tumor-bam-file "$TUMOR" -r "$REFERENCE" \
    -t "$THREADS" --tagSupplementary -q 20 --region "$REGION" --output-somatic-vcf -o "$PREFIX" \
    > "$LOG" 2>&1
lps_rc="$?"
set -e
if [[ "$lps_rc" -ne 0 ]]; then
    kill "$capture_pid" 2>/dev/null || true
    wait "$capture_pid" 2>/dev/null || true
    jq -n --argjson exit_code "$lps_rc" '{status:"FAILED",stage:"longphase",exit_code:$exit_code}' > "$RUN/_FAILED"
    exit "$lps_rc"
fi
wait "$capture_pid"
mv "${PREFIX}.bam" "$RUN/consumed_tagged_bam.fifo"

"$PY" "$SCRIPTS/validate_streamed_longphase_sidecar.py" \
    --sidecar "$SIDECAR" --capture-summary "$RUN/stream_capture_summary.json" \
    --execution-log "$LOG" --input-vcf "$INPUT" --recalibrated-vcf "${PREFIX}_sc.vcf" \
    --region chrX --output "$RUN/sidecar_validation.json"
"$PY" "$SCRIPTS/audit_longphase_filter_transitions.py" \
    --input-vcf "$INPUT" --output-vcf "${PREFIX}_sc.vcf" \
    --output-json "$RUN/filter_transition_audit.json"

bgzip -@ "$THREADS" -c "$SIDECAR" > "$RUN/HCC1395.chrX_72880028.read_tags.tsv.gz"
tabix -0 -s 1 -b 2 -e 3 "$RUN/HCC1395.chrX_72880028.read_tags.tsv.gz"
sha256sum "$RUN"/*.json "$RUN"/*.tsv "$RUN"/*.tsv.gz "$RUN"/*.tbi \
    "$RUN"/*.vcf "$RUN"/*.vcf.gz "$RUN"/*.csi "$RUN/command.sh.txt" "$RUN/input_and_code.sha256" \
    > "$RUN/artifacts.sha256"
jq -n --arg region "$REGION" --argjson rescued "$(jq '.rescued_nonpass_to_pass' "$RUN/filter_transition_audit.json")" \
    '{status:"SUCCESS",scope:$region,input_records:1,output_records:1,rescued:$rescued}' > "$RUN/_SUCCESS"
echo "KNOWN-CRASH LOCUS PROBE PASS -> $RUN"
