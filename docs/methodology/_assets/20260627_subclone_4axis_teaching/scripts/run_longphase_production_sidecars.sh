#!/usr/bin/env bash
# Production LongPhase-S tag sidecars: no truth flags and no persisted BAM payload.
set -euo pipefail

SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SM_PYTHON:-/bip7_disk/liaoyoyo2001/miniconda3/bin/python3}"
LPS="${SM_LONGPHASE_S:-/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s}"
MANIFEST="${SM_LPS_MANIFEST:-$SCD/../data/longphase_production_sidecar_manifest.json}"
RUN_ID="${SM_RUN_ID:-$(date '+%Y%m%d_%H%M%S')_longphase_s_production_sidecars}"
RUN_ROOT="${SM_RUN_ROOT:-/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/$RUN_ID}"
THREADS="${SM_LPS_THREADS:-12}"
PARALLEL="${SM_PARALLEL_SAMPLES:-2}"
STATUS="$RUN_ROOT/run_status.tsv"

[[ -x "$PY" ]] || { echo "ERROR Python: $PY" >&2; exit 2; }
[[ -x "$LPS" ]] || { echo "ERROR LongPhase-S: $LPS" >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "ERROR manifest: $MANIFEST" >&2; exit 2; }
[[ ! -e "$RUN_ROOT" ]] || { echo "ERROR immutable output exists: $RUN_ROOT" >&2; exit 2; }
mkdir -p "$RUN_ROOT/samples"
cp "$MANIFEST" "$RUN_ROOT/input_manifest.json"
printf 'timestamp\tsample\tstage\tstatus\tdetail\n' > "$STATUS"

record_status() {
    printf '%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "$1" "$2" "$3" "${4:-}" >> "$STATUS"
}

value() {
    jq -r --arg sample "$1" --arg field "$2" '.samples[] | select(.sample==$sample) | .[$field] // ""' "$MANIFEST"
}

run_sample() {
    local sample="$1" wd="$RUN_ROOT/samples/$1"
    local germline normal tumor reference caller_pass expected prefix log
    germline="$(value "$sample" germline_phased_vcf)"
    normal="$(value "$sample" normal_bam)"
    tumor="$(value "$sample" tumor_bam)"
    reference="$(value "$sample" reference)"
    caller_pass="$(value "$sample" caller_pass_vcf)"
    expected="$(value "$sample" caller_pass_records)"
    mkdir -p "$wd"
    for path in "$germline" "$normal" "$tumor" "$reference" "$caller_pass"; do
        [[ -f "$path" ]] || { record_status "$sample" preflight FAIL "missing=$path"; return 1; }
    done
    prefix="$wd/${sample}_production"
    log="$wd/longphase_s.production.log"
    mkfifo "${prefix}.bam"
    printf '%q ' "$LPS" somatic_haplotag -s "$germline" -b "$normal" --tumor-snv-file "$caller_pass" \
        --tumor-bam-file "$tumor" -r "$reference" -t "$THREADS" --tagSupplementary -q 20 \
        --output-somatic-vcf -o "$prefix" > "$wd/command.sh.txt"
    printf '\n' >> "$wd/command.sh.txt"
    {
        printf 'role\tpath\tsize_bytes\tmtime\n'
        for spec in "germline:$germline" "normal_bam:$normal" "tumor_bam:$tumor" \
                    "reference:$reference" "caller_pass:$caller_pass"; do
            local role="${spec%%:*}" path="${spec#*:}"
            printf '%s\t%s\t%s\t%s\n' "$role" "$path" "$(stat -c %s "$path")" "$(stat -c %y "$path")"
        done
    } > "$wd/input_files.tsv"
    sha256sum "$caller_pass" "${caller_pass}.csi" "$germline" "${germline}.csi" \
        "${normal}.bai" "${tumor}.bai" > "$wd/input.sha256"
    record_status "$sample" production_tagging START "expected_records=$expected;truth_flags=absent"
    "$PY" "$SCD/capture_longphase_tagged_bam_sidecar.py" --input-bam "${prefix}.bam" \
        --output "$wd/${sample}.read_tags.tsv" --summary "$wd/stream_capture_summary.json" \
        > "$wd/stream_capture.log" 2>&1 &
    local capture_pid="$!" rc capture_rc
    set -o pipefail
    set +e
    /usr/bin/time -v "$LPS" somatic_haplotag -s "$germline" -b "$normal" \
        --tumor-snv-file "$caller_pass" --tumor-bam-file "$tumor" -r "$reference" -t "$THREADS" \
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
    "$PY" "$SCD/validate_streamed_longphase_sidecar.py" --sidecar "$wd/${sample}.read_tags.tsv" \
        --capture-summary "$wd/stream_capture_summary.json" --execution-log "$log" \
        --input-vcf "$caller_pass" --recalibrated-vcf "${prefix}_sc.vcf" \
        --output "$wd/sidecar_validation.json"
    bgzip -@ "$THREADS" -c "$wd/${sample}.read_tags.tsv" > "$wd/${sample}.read_tags.tsv.gz"
    tabix -0 -s 1 -b 2 -e 3 "$wd/${sample}.read_tags.tsv.gz"
    bgzip -@ "$THREADS" -c "${prefix}_sc.vcf" > "$wd/${sample}.longphase_s.recalibrated.all.vcf.gz"
    bcftools index -c "$wd/${sample}.longphase_s.recalibrated.all.vcf.gz"
    bcftools view -f PASS -Oz -o "$wd/${sample}.longphase_s.recalibrated.pass.vcf.gz" "${prefix}_sc.vcf"
    bcftools index -c "$wd/${sample}.longphase_s.recalibrated.pass.vcf.gz"
    sha256sum "$wd/${sample}.read_tags.tsv" "$wd/${sample}.read_tags.tsv.gz" \
        "$wd/${sample}.longphase_s.recalibrated.all.vcf.gz" \
        "$wd/${sample}.longphase_s.recalibrated.pass.vcf.gz" "$wd/sidecar_validation.json" > "$wd/output.sha256"
    record_status "$sample" production_tagging PASS "rows=$(jq -r '.capture.rows_mapped' "$wd/sidecar_validation.json")"
    echo "[$sample] PRODUCTION SIDECAR PASS -> $wd"
}

export -f record_status value run_sample
export SCD PY LPS MANIFEST RUN_ROOT THREADS STATUS
sha256sum "$MANIFEST" "$SCD/run_longphase_production_sidecars.sh" \
    "$SCD/capture_longphase_tagged_bam_sidecar.py" "$SCD/validate_streamed_longphase_sidecar.py" \
    > "$RUN_ROOT/code.sha256"
jq -n --arg run_root "$RUN_ROOT" --argjson threads "$THREADS" --argjson parallel "$PARALLEL" \
    '{run_root:$run_root,threads:$threads,parallel_samples:$parallel,truth_flags:false,mapq:20,tag_supplementary:true,output_mode:"read_tag_sidecar"}' \
    > "$RUN_ROOT/params.json"
echo "=== LONGPHASE-S PRODUCTION SIDECARS START $(date --iso-8601=seconds) ==="
echo "INPUT MANIFEST: $MANIFEST"
echo "OUTPUT ROOT: $RUN_ROOT"
mapfile -t SAMPLES < <(jq -r '.samples[].sample' "$MANIFEST")
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$PARALLEL" -I {} bash -c 'set -euo pipefail; run_sample "$1"' _ {}
jq -s '{schema_version:"1.0",dataset_count:length,n_pass:map(select(.pass==true))|length,all_pass:all(.pass==true),samples:.}' \
    "$RUN_ROOT"/samples/*/sidecar_validation.json > "$RUN_ROOT/verification_summary.json"
jq -e '.dataset_count == 7 and .n_pass == 7 and .all_pass == true' "$RUN_ROOT/verification_summary.json" >/dev/null
record_status ALL verify PASS "7/7 production sidecars"
echo "=== LONGPHASE-S PRODUCTION SIDECARS COMPLETE $(date --iso-8601=seconds) ==="
echo "OUTPUT ROOT: $RUN_ROOT"
