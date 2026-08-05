#!/usr/bin/env bash
# Production normalized-raw-all LongPhase-S sidecars with bidirectional FILTER audit.
set -euo pipefail

SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SM_PYTHON:-/bip7_disk/liaoyoyo2001/miniconda3/bin/python3}"
LPS_SOURCE="${SM_LONGPHASE_S:-/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_zero_read_patch_probe_v1/source_snapshot/longphase-s}"
MANIFEST="${SM_LPS_MANIFEST:-$SCD/../data/longphase_raw_all_production_manifest_v2.json}"
RUN_ID="${SM_RUN_ID:-$(date '+%Y%m%d_%H%M%S')_longphase_s_raw_all_production_sidecars}"
RUN_ROOT="${SM_RUN_ROOT:-/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/$RUN_ID}"
THREADS="${SM_LPS_THREADS:-12}"
PARALLEL="${SM_PARALLEL_SAMPLES:-2}"
STATUS="$RUN_ROOT/run_status.tsv"
EXPECTED_BINARY_SHA256="5ceba723d31c52f01202478de19952219371f0b7f9c136b8882cdd93026789b8"

[[ -x "$PY" ]] || { echo "ERROR Python: $PY" >&2; exit 2; }
[[ -x "$LPS_SOURCE" ]] || { echo "ERROR LongPhase-S: $LPS_SOURCE" >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "ERROR manifest: $MANIFEST" >&2; exit 2; }
[[ ! -e "$RUN_ROOT" ]] || { echo "ERROR immutable output exists: $RUN_ROOT" >&2; exit 2; }
PATCH_RECEIPT_SOURCE="$(jq -r '.longphase_binary.patch_receipt' "$MANIFEST")"
PATCH_RECEIPT_EXPECTED_SHA256="$(jq -r '.longphase_binary.patch_receipt_sha256' "$MANIFEST")"
[[ "$(sha256sum "$LPS_SOURCE" | awk '{print $1}')" == "$EXPECTED_BINARY_SHA256" ]] || {
    echo "ERROR patched binary hash mismatch: $LPS_SOURCE" >&2; exit 2;
}
[[ -f "$PATCH_RECEIPT_SOURCE" && "$(sha256sum "$PATCH_RECEIPT_SOURCE" | awk '{print $1}')" == "$PATCH_RECEIPT_EXPECTED_SHA256" ]] || {
    echo "ERROR patch receipt hash mismatch" >&2; exit 2;
}
PATCH_FILE_SOURCE="$(jq -r '.patch' "$PATCH_RECEIPT_SOURCE")"
[[ -f "$PATCH_FILE_SOURCE" ]] || { echo "ERROR source patch missing: $PATCH_FILE_SOURCE" >&2; exit 2; }
jq -e '.status == "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION"
       and .approval.scope == "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY"' \
    "$PATCH_RECEIPT_SOURCE" >/dev/null || { echo "ERROR patch receipt is not approved" >&2; exit 2; }
mapfile -t SAMPLES < <(jq -r '.samples[].sample' "$MANIFEST")
[[ "${SAMPLES[*]}" == "HCC1395 HCC1395_DORADO COLO829 H1437 H2009 HCC1937 HCC1954" ]] || {
    echo "ERROR manifest dataset order/set differs" >&2; exit 2;
}
jq -e --arg hash "$EXPECTED_BINARY_SHA256" --arg path "$(realpath "$LPS_SOURCE")" \
    '.schema_version == "2.0" and .dataset_count == 7 and .biological_sample_count == 6
     and .contract.somatic_input_role == "normalized raw ClairS all records; bidirectional FILTER recalibration"
     and .longphase_binary.sha256 == $hash and .longphase_binary.path == $path
     and all(.samples[]; .longphase_input_contract == "normalized_ClairS_raw_all")' \
    "$MANIFEST" >/dev/null || { echo "ERROR raw-all manifest contract invalid" >&2; exit 2; }
mkdir -p "$RUN_ROOT/samples" "$RUN_ROOT/source_bundle"
cp "$MANIFEST" "$RUN_ROOT/input_manifest.json"
MANIFEST="$RUN_ROOT/input_manifest.json"
cp "$LPS_SOURCE" "$RUN_ROOT/source_bundle/longphase-s"
cp "$SCD/capture_longphase_tagged_bam_sidecar.py" "$RUN_ROOT/source_bundle/"
cp "$SCD/validate_streamed_longphase_sidecar.py" "$RUN_ROOT/source_bundle/"
cp "$SCD/audit_longphase_filter_transitions.py" "$RUN_ROOT/source_bundle/"
cp "$SCD/build_longphase_raw_all_capture_receipt_v2.py" "$RUN_ROOT/source_bundle/"
cp "$SCD/finalize_longphase_raw_all_capture_receipts.py" "$RUN_ROOT/source_bundle/"
cp "$SCD/build_longphase_production_capture_receipt_v3.py" "$RUN_ROOT/source_bundle/"
cp "$SCD/validate_layered_v3_inputs.py" "$RUN_ROOT/source_bundle/"
cp "$SCD/../schemas/longphase_raw_all_capture_receipt_v2.schema.json" "$RUN_ROOT/source_bundle/"
cp "$SCD/$(basename "$0")" "$RUN_ROOT/source_bundle/"
cp "$PATCH_RECEIPT_SOURCE" "$RUN_ROOT/source_bundle/"
cp "$PATCH_FILE_SOURCE" "$RUN_ROOT/source_bundle/"
LPS="$RUN_ROOT/source_bundle/longphase-s"
CAPTURE="$RUN_ROOT/source_bundle/capture_longphase_tagged_bam_sidecar.py"
VALIDATE="$RUN_ROOT/source_bundle/validate_streamed_longphase_sidecar.py"
TRANSITIONS="$RUN_ROOT/source_bundle/audit_longphase_filter_transitions.py"
printf 'timestamp\tsample\tstage\tstatus\tdetail\n' > "$STATUS"

write_root_failure() {
    local rc="$?"
    trap - EXIT
    if [[ "$rc" -ne 0 && ! -e "$RUN_ROOT/_FAILED" && ! -e "$RUN_ROOT/_SUCCESS" ]]; then
        jq -n --argjson exit_code "$rc" --arg stage "producer" \
            '{status:"FAILED",stage:$stage,exit_code:$exit_code}' > "$RUN_ROOT/_FAILED"
    fi
    exit "$rc"
}
trap write_root_failure EXIT

record_status() {
    {
        flock 9
        printf '%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "$1" "$2" "$3" "${4:-}" >> "$STATUS"
    } 9>"$RUN_ROOT/.status.lock"
}

value() {
    jq -r --arg sample "$1" --arg field "$2" '.samples[] | select(.sample==$sample) | .[$field] // ""' "$MANIFEST"
}

run_sample() {
    local sample="$1" wd="$RUN_ROOT/samples/$1"
    local germline normal tumor reference caller_raw caller_raw_index caller_pass caller_pass_index
    local expected prefix log normalized warning_count
    germline="$(value "$sample" germline_phased_vcf)"
    normal="$(value "$sample" normal_bam)"
    tumor="$(value "$sample" tumor_bam)"
    reference="$(value "$sample" reference)"
    caller_raw="$(value "$sample" caller_raw_vcf)"
    caller_raw_index="$(value "$sample" caller_raw_vcf_index)"
    caller_pass="$(value "$sample" caller_pass_vcf)"
    caller_pass_index="$(value "$sample" caller_pass_vcf_index)"
    expected="$(jq -r --arg sample "$sample" '.samples[] | select(.sample==$sample) | .caller_raw_scan.records' "$MANIFEST")"
    mkdir -p "$wd"
    for path in "$germline" "$normal" "$tumor" "$reference" "$caller_raw" "$caller_raw_index" \
                "$caller_pass" "$caller_pass_index"; do
        [[ -f "$path" ]] || { record_status "$sample" preflight FAIL "missing=$path"; return 1; }
    done
    normalized="$wd/${sample}.clairs.raw_all.normalized.vcf.gz"
    zcat "$caller_raw" \
        | sed 's/ID=GQ,Number=1,Type=Integer/ID=GQ,Number=1,Type=Float/' \
        | bcftools view -Oz -o "$normalized"
    bcftools index -c "$normalized"
    "$PY" "$TRANSITIONS" --input-vcf "$caller_raw" --output-vcf "$normalized" \
        --output-json "$wd/normalization_audit.json"
    jq -e --argjson expected "$expected" \
        '.pass == true and .input.record_count == $expected and .output.record_count == $expected
         and .rescued_nonpass_to_pass == 0 and .removed_pass_to_nonpass == 0' \
        "$wd/normalization_audit.json" >/dev/null
    prefix="$wd/${sample}_production"
    log="$wd/longphase_s.production.log"
    mkfifo "${prefix}.bam"
    printf '%q ' "$LPS" somatic_haplotag -s "$germline" -b "$normal" --tumor-snv-file "$normalized" \
        --tumor-bam-file "$tumor" -r "$reference" -t "$THREADS" --tagSupplementary -q 20 \
        --output-somatic-vcf -o "$prefix" > "$wd/command.sh.txt"
    printf '\n' >> "$wd/command.sh.txt"
    {
        printf 'role\tpath\tsize_bytes\tmtime\n'
        for spec in "germline:$germline" "normal_bam:$normal" "tumor_bam:$tumor" \
                    "reference:$reference" "caller_raw:$caller_raw" "longphase_input:$normalized" \
                    "caller_pass_baseline:$caller_pass"; do
            local role="${spec%%:*}" path="${spec#*:}"
            printf '%s\t%s\t%s\t%s\n' "$role" "$path" "$(stat -Lc %s "$path")" "$(stat -Lc %y "$path")"
        done
    } > "$wd/input_files.tsv"
    sha256sum "$caller_raw" "$caller_raw_index" "$normalized" "${normalized}.csi" \
        "$caller_pass" "$caller_pass_index" "$germline" "${germline}.csi" \
        "${normal}.bai" "${tumor}.bai" "$LPS" > "$wd/input.sha256"
    record_status "$sample" production_tagging START \
        "expected_raw_records=$expected;input_contract=normalized_raw_all;truth_flags=absent"
    "$PY" "$CAPTURE" --input-bam "${prefix}.bam" \
        --output "$wd/${sample}.read_tags.tsv" --summary "$wd/stream_capture_summary.json" \
        > "$wd/stream_capture.log" 2>&1 &
    local capture_pid="$!" rc capture_rc
    set -o pipefail
    set +e
    /usr/bin/time -v "$LPS" somatic_haplotag -s "$germline" -b "$normal" \
        --tumor-snv-file "$normalized" --tumor-bam-file "$tumor" -r "$reference" -t "$THREADS" \
        --tagSupplementary -q 20 --output-somatic-vcf -o "$prefix" 2>&1 | tee "$log"
    rc="${PIPESTATUS[0]}"
    set -e
    if [[ "$rc" -ne 0 ]]; then
        kill "$capture_pid" 2>/dev/null || true
        wait "$capture_pid" 2>/dev/null || true
        record_status "$sample" production_tagging FAIL "longphase_exit=$rc"
        return "$rc"
    fi
    set +e
    wait "$capture_pid"
    capture_rc="$?"
    set -e
    if [[ "$capture_rc" -ne 0 ]]; then
        record_status "$sample" production_tagging FAIL "capture_exit=$capture_rc"
        return "$capture_rc"
    fi
    mv "${prefix}.bam" "$wd/consumed_tagged_bam.fifo"
    [[ -s "${prefix}_sc.vcf" && -s "${prefix}_purity.out" && -s "$wd/${sample}.read_tags.tsv" ]]
    "$PY" "$VALIDATE" --sidecar "$wd/${sample}.read_tags.tsv" \
        --capture-summary "$wd/stream_capture_summary.json" --execution-log "$log" \
        --input-vcf "$normalized" --recalibrated-vcf "${prefix}_sc.vcf" \
        --output "$wd/sidecar_validation.json"
    "$PY" "$TRANSITIONS" --input-vcf "$normalized" --output-vcf "${prefix}_sc.vcf" \
        --output-json "$wd/filter_transition_audit.json"
    bgzip -@ "$THREADS" -c "$wd/${sample}.read_tags.tsv" > "$wd/${sample}.read_tags.tsv.gz"
    tabix -0 -s 1 -b 2 -e 3 "$wd/${sample}.read_tags.tsv.gz"
    bgzip -@ "$THREADS" -c "${prefix}_sc.vcf" > "$wd/${sample}.longphase_s.recalibrated.all.vcf.gz"
    bcftools index -c "$wd/${sample}.longphase_s.recalibrated.all.vcf.gz"
    bcftools view -f PASS -Oz -o "$wd/${sample}.longphase_s.recalibrated.pass.vcf.gz" "${prefix}_sc.vcf"
    bcftools index -c "$wd/${sample}.longphase_s.recalibrated.pass.vcf.gz"
    warning_count="$(grep -c 'no eligible tumor read at low-confidence site' "$log" || true)"
    jq -n --arg sample "$sample" --argjson expected "$expected" --argjson warnings "$warning_count" \
        --slurpfile sidecar "$wd/sidecar_validation.json" \
        --slurpfile transitions "$wd/filter_transition_audit.json" \
        --slurpfile normalization "$wd/normalization_audit.json" \
        '{sample:$sample,expected_raw_records:$expected,zero_read_warnings:$warnings,
          normalization:$normalization[0],filter_transitions:$transitions[0],sidecar:$sidecar[0],
          pass:($normalization[0].pass and $transitions[0].pass and $sidecar[0].pass
                and $transitions[0].input.record_count==$expected
                and $transitions[0].output.record_count==$expected)}' \
        > "$wd/sample_verification.json"
    jq -e '.pass == true' "$wd/sample_verification.json" >/dev/null
    [[ -p "$wd/consumed_tagged_bam.fifo" ]]
    [[ -z "$(find "$wd" -maxdepth 1 -type f -name '*.bam' -print -quit)" ]]
    sha256sum "$wd/${sample}.read_tags.tsv" "$wd/${sample}.read_tags.tsv.gz" \
        "$wd/${sample}.read_tags.tsv.gz.tbi" \
        "$wd/${sample}.longphase_s.recalibrated.all.vcf.gz" \
        "$wd/${sample}.longphase_s.recalibrated.all.vcf.gz.csi" \
        "$wd/${sample}.longphase_s.recalibrated.pass.vcf.gz" \
        "$wd/${sample}.longphase_s.recalibrated.pass.vcf.gz.csi" \
        "$wd/normalization_audit.json" "$wd/filter_transition_audit.json" \
        "$wd/sidecar_validation.json" "$wd/stream_capture_summary.json" \
        "$wd/sample_verification.json" > "$wd/output.sha256"
    jq -n --arg sample "$sample" --argjson rows "$(jq -r '.capture.rows_mapped' "$wd/sidecar_validation.json")" \
        --argjson rescued "$(jq -r '.rescued_nonpass_to_pass' "$wd/filter_transition_audit.json")" \
        --argjson removed "$(jq -r '.removed_pass_to_nonpass' "$wd/filter_transition_audit.json")" \
        --argjson warnings "$warning_count" \
        '{status:"SUCCESS",sample:$sample,sidecar_rows:$rows,rescued:$rescued,removed:$removed,
          zero_read_warnings:$warnings}' > "$wd/_SUCCESS"
    record_status "$sample" production_tagging PASS \
        "rows=$(jq -r '.sidecar_rows' "$wd/_SUCCESS");rescued=$(jq -r '.rescued' "$wd/_SUCCESS");removed=$(jq -r '.removed' "$wd/_SUCCESS");warnings=$warning_count"
    echo "[$sample] RAW-ALL PRODUCTION SIDECAR PASS -> $wd"
}

export -f record_status value run_sample
export SCD PY LPS CAPTURE VALIDATE TRANSITIONS MANIFEST RUN_ROOT THREADS STATUS
sha256sum "$RUN_ROOT/input_manifest.json" "$RUN_ROOT/source_bundle/"* "$PY" \
    > "$RUN_ROOT/code.sha256"
jq -n --arg run_root "$RUN_ROOT" --arg binary_sha256 "$EXPECTED_BINARY_SHA256" \
    --argjson threads "$THREADS" --argjson parallel "$PARALLEL" \
    '{run_root:$run_root,threads:$threads,parallel_samples:$parallel,truth_flags:false,mapq:20,
      tag_supplementary:true,output_mode:"read_tag_sidecar",longphase_input_contract:"normalized_ClairS_raw_all",
      filter_contract:"bidirectional_recalibration",patched_binary_sha256:$binary_sha256}' \
    > "$RUN_ROOT/params.json"
echo "=== LONGPHASE-S RAW-ALL PRODUCTION SIDECARS START $(date --iso-8601=seconds) ==="
echo "INPUT MANIFEST: $MANIFEST"
echo "OUTPUT ROOT: $RUN_ROOT"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$PARALLEL" -I {} bash -c 'set -euo pipefail; run_sample "$1"' _ {}
jq -s '{schema_version:"2.0",input_contract:"normalized_ClairS_raw_all",dataset_count:length,
        n_pass:(map(select(.pass==true))|length),all_pass:all(.pass==true),
        total_raw_records:(map(.filter_transitions.input.record_count)|add),
        total_rescued:(map(.filter_transitions.rescued_nonpass_to_pass)|add),
        total_removed:(map(.filter_transitions.removed_pass_to_nonpass)|add),
        total_zero_read_warnings:(map(.zero_read_warnings)|add),samples:.}' \
    "$RUN_ROOT"/samples/*/sample_verification.json > "$RUN_ROOT/verification_summary.json"
jq -e '.dataset_count == 7 and .n_pass == 7 and .all_pass == true' "$RUN_ROOT/verification_summary.json" >/dev/null
[[ "$(find "$RUN_ROOT/samples" -mindepth 2 -maxdepth 2 -name _SUCCESS -type f | wc -l)" -eq 7 ]]
record_status ALL verify PASS "7/7 production sidecars"
sha256sum "$RUN_ROOT/input_manifest.json" "$RUN_ROOT/params.json" "$RUN_ROOT/code.sha256" \
    "$RUN_ROOT/run_status.tsv" "$RUN_ROOT/verification_summary.json" \
    "$RUN_ROOT"/samples/*/input.sha256 "$RUN_ROOT"/samples/*/output.sha256 \
    "$RUN_ROOT"/samples/*/_SUCCESS > "$RUN_ROOT/closeout.sha256"
jq -n --arg status "SUCCESS" --arg input_contract "normalized_ClairS_raw_all" \
    --arg binary_sha256 "$EXPECTED_BINARY_SHA256" \
    --argjson summary "$(cat "$RUN_ROOT/verification_summary.json")" \
    '{schema_version:"2.0",status:$status,input_contract:$input_contract,
      patched_binary_sha256:$binary_sha256,verification:$summary}' \
    > "$RUN_ROOT/production_closeout.json"
sha256sum "$RUN_ROOT/production_closeout.json" "$RUN_ROOT/closeout.sha256" > "$RUN_ROOT/root_receipt.sha256"
jq -n --arg closeout "$RUN_ROOT/production_closeout.json" \
    --arg closeout_sha256 "$(sha256sum "$RUN_ROOT/production_closeout.json" | awk '{print $1}')" \
    --arg receipt "$RUN_ROOT/root_receipt.sha256" \
    '{status:"SUCCESS",closeout:$closeout,closeout_sha256:$closeout_sha256,root_receipt:$receipt}' \
    > "$RUN_ROOT/_SUCCESS.pending"
mv "$RUN_ROOT/_SUCCESS.pending" "$RUN_ROOT/_SUCCESS"
echo "=== LONGPHASE-S RAW-ALL PRODUCTION SIDECARS COMPLETE $(date --iso-8601=seconds) ==="
echo "OUTPUT ROOT: $RUN_ROOT"
