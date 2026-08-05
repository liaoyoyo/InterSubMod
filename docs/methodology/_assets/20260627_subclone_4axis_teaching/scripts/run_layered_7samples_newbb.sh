#!/usr/bin/env bash
# Layered reconstruction v2: 7 datasets, chr1-22, explicit manifest, immutable run root.
set -euo pipefail

SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SM_PYTHON:-/bip7_disk/liaoyoyo2001/miniconda3/bin/python3}"
MANIFEST="${SM_INPUT_MANIFEST:-$SCD/../data/layered_v2_input_manifest.json}"
SOURCE_MANIFEST="$MANIFEST"
RUN_ID="${SM_RUN_ID:-$(date '+%Y%m%d_%H%M%S')_layered_reconstruction_v2}"
RUN_ROOT="${SM_RUN_ROOT:-/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/$RUN_ID}"
PARALLEL_SAMPLES="${SM_PARALLEL_SAMPLES:-2}"
VERIFY_EVERY="${SM_VERIFY_EVERY:-1}"
ANALYSIS_TREE_CAP="${SM_ANALYSIS_TREE_CAP:-0}"
DISPLAY_TREE_CAP="${SM_DISPLAY_TREE_CAP:-32}"
TREE_CONTRACT="$(jq -r '.tree_input_contract // ""' "$MANIFEST" 2>/dev/null || true)"
DATASET_COUNT="$(jq -r '.dataset_count' "$MANIFEST" 2>/dev/null || true)"
BIOLOGICAL_COUNT="$(jq -r '.biological_sample_count' "$MANIFEST" 2>/dev/null || true)"
TIME=/usr/bin/time
PROFILE="$RUN_ROOT/profile.log"
STATUS="$RUN_ROOT/run_status.tsv"
SPLITS=("chr1,chr6,chr11,chr16,chr21:1" "chr2,chr7,chr12,chr17,chr22:2" \
        "chr3,chr8,chr13,chr18:3" "chr4,chr9,chr14,chr19:4" "chr5,chr10,chr15,chr20:5")
SOURCE_FILES=(sm_linkage_genomewide.py sm_multilocus_combinations.py tree_enumeration_solver.py
              layered_tree_reconstruction.py build_region_view.py build_ssnv_site_ledger.py
              verify_layered_v2.py validate_layered_v2_inputs.py run_layered_7samples_newbb.sh)

[[ -x "$PY" ]] || { echo "ERROR: Python not executable: $PY" >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest missing: $MANIFEST" >&2; exit 2; }
jq -e '.schema_version == "2.1"
  and .dataset_count == 7
  and .biological_sample_count == 6
  and .dataset_count == (.samples | length)
  and .tag_contract.truth_flags == false
  and .tag_contract.PS_preserved == true
  and .tag_contract.bam_identity_locked == true
  and .tag_contract.longphase_filtering_policy == "production_default_filter"
  and (.tag_contract.production_closeout | type == "string" and length > 0)
  and (.tag_contract.production_closeout_sha256 | type == "string" and length == 64)
  and (.tag_contract.production_success | type == "string" and length > 0)
  and (.tag_contract.production_success_sha256 | type == "string" and length == 64)
  and (.tag_contract.production_artifacts_manifest | type == "string" and length > 0)
  and (.tag_contract.production_artifacts_manifest_sha256 | type == "string" and length == 64)
  and ((.tree_input_contract == "longphase_recalibrated_PASS"
        and .tag_contract.tree_backbone == "LongPhase-S _sc.vcf FILTER=PASS")
       or (.task_type == "B_BACKBONE_SENSITIVITY"
           and .tree_input_contract == "clairs_PASS_input"
           and .tag_contract.tree_backbone == "ClairS PASS sensitivity"))
  and all(.samples[];
    (.read_tag_sidecar | type == "string" and length > 0)
    and (.read_tag_validation | type == "string" and length > 0)
    and (.caller_raw_vcf | type == "string" and length > 0)
    and (.longphase_input_vcf | type == "string" and length > 0)
    and (.longphase_recalibrated_all_vcf | type == "string" and length > 0)
    and (.longphase_production_closeout | type == "string" and length > 0)
    and .longphase_tagging_scope == "production genome-wide; no truth VCF/BED flags")' \
  "$MANIFEST" >/dev/null
if [[ -e "$RUN_ROOT" ]]; then
  echo "ERROR: immutable run root already exists: $RUN_ROOT" >&2
  exit 2
fi
PREFLIGHT_PENDING="${RUN_ROOT}.input_preflight.pending.$$"
SOURCE_MANIFEST_SHA="$(sha256sum "$SOURCE_MANIFEST" | cut -d' ' -f1)"
VALIDATOR_SHA="$(sha256sum "$SCD/validate_layered_v2_inputs.py" | cut -d' ' -f1)"
"$PY" "$SCD/validate_layered_v2_inputs.py" --manifest "$SOURCE_MANIFEST" \
  --output "$PREFLIGHT_PENDING"
[[ "$SOURCE_MANIFEST_SHA" == "$(sha256sum "$SOURCE_MANIFEST" | cut -d' ' -f1)" ]] || {
  echo "ERROR: source manifest changed during preflight" >&2
  exit 2
}
[[ "$VALIDATOR_SHA" == "$(sha256sum "$SCD/validate_layered_v2_inputs.py" | cut -d' ' -f1)" ]] || {
  echo "ERROR: validator changed during preflight" >&2
  exit 2
}

mkdir -p "$RUN_ROOT/samples" "$RUN_ROOT/source"
cp "$SOURCE_MANIFEST" "$RUN_ROOT/input_manifest.json"
mv "$PREFLIGHT_PENDING" "$RUN_ROOT/input_preflight.json"
for source_file in "${SOURCE_FILES[@]}"; do
  cp "$SCD/$source_file" "$RUN_ROOT/source/$source_file"
done
MANIFEST="$RUN_ROOT/input_manifest.json"
EXEC_SCD="$RUN_ROOT/source"
sha256sum "$EXEC_SCD"/* > "$RUN_ROOT/source_bundle.sha256"
printf 'timestamp\tsample\tstage\tstatus\tdetail\n' > "$STATUS"

RUN_COMPLETED=0
mark_run_end() {
  local rc="$?"
  trap - EXIT
  if [[ "$rc" -ne 0 && "$RUN_COMPLETED" -eq 0 && -d "$RUN_ROOT" && ! -e "$RUN_ROOT/_SUCCESS" ]]; then
    {
      printf 'status\tFAILED\n'
      printf 'exit_code\t%s\n' "$rc"
      printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
    } > "$RUN_ROOT/_FAILED.pending.$$"
    mv "$RUN_ROOT/_FAILED.pending.$$" "$RUN_ROOT/_FAILED"
  fi
  exit "$rc"
}
trap mark_run_end EXIT
trap 'exit 130' INT TERM

record_status() {
  local sample="$1" stage="$2" status="$3" detail="${4:-}"
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "$sample" "$stage" "$status" "$detail" >> "$STATUS"
}

manifest_value() {
  local sample="$1" field="$2"
  jq -r --arg sample "$sample" --arg field "$field" \
    '.samples[] | select(.sample == $sample) | .[$field] // ""' "$MANIFEST"
}

validate_input() {
  local sample="$1" path="$2" label="$3"
  if [[ ! -f "$path" ]]; then
    record_status "$sample" preflight FAIL "$label missing: $path"
    echo "[$sample] ERROR: $label missing: $path" >&2
    return 1
  fi
}

run_sample() {
  local sample="$1"
  local wd="$RUN_ROOT/samples/$sample"
  local splits=("chr1,chr6,chr11,chr16,chr21:1" "chr2,chr7,chr12,chr17,chr22:2" \
                "chr3,chr8,chr13,chr18:3" "chr4,chr9,chr14,chr19:4" "chr5,chr10,chr15,chr20:5")
  local tbam treevcf lpsinput rawvcf recalall readtags readtagsidx tagvalidation lpsinventory cnbed cnsource cngain cnloss integration backbone_source
  tbam="$(manifest_value "$sample" tumor_bam)"
  treevcf="$(manifest_value "$sample" somatic_vcf)"
  lpsinput="$(manifest_value "$sample" longphase_input_vcf)"
  rawvcf="$(manifest_value "$sample" caller_raw_vcf)"
  recalall="$(manifest_value "$sample" longphase_recalibrated_all_vcf)"
  readtags="$(manifest_value "$sample" read_tag_sidecar)"
  readtagsidx="$(manifest_value "$sample" read_tag_sidecar_index)"
  tagvalidation="$(manifest_value "$sample" read_tag_validation)"
  lpsinventory="$(manifest_value "$sample" longphase_input_inventory)"
  backbone_source="$(manifest_value "$sample" somatic_vcf_role)"
  cnbed="$(manifest_value "$sample" cn_bed)"
  cnsource="$(manifest_value "$sample" cn_source)"
  cngain="$(manifest_value "$sample" cn_int_gain)"
  cnloss="$(manifest_value "$sample" cn_int_loss)"
  integration="$(manifest_value "$sample" integration_json)"

  mkdir -p "$wd"
  validate_input "$sample" "$tbam" tumor_bam
  validate_input "$sample" "${tbam}.bai" tumor_bam_index
  validate_input "$sample" "$treevcf" tree_input_vcf
  validate_input "$sample" "$lpsinput" longphase_input_vcf
  validate_input "$sample" "$rawvcf" caller_raw_vcf
  validate_input "$sample" "$recalall" longphase_recalibrated_all_vcf
  validate_input "$sample" "$readtags" read_tag_sidecar
  validate_input "$sample" "$readtagsidx" read_tag_sidecar_index
  validate_input "$sample" "$tagvalidation" read_tag_validation
  validate_input "$sample" "$lpsinventory" longphase_input_inventory
  jq -e '.pass == true
    and .region == "all"
    and .duplicate_exact_alignment_rows == 0
    and .duplicate_exact_alignment_conflicts == 0
    and .checks.truth_flags_absent == true
    and .checks.parser_count_matches_input == true
    and .checks.recalibrated_preserves_all_input_keys == true' "$tagvalidation" >/dev/null
  if [[ -n "$cnbed" ]]; then validate_input "$sample" "$cnbed" cn_bed; fi
  if [[ -n "$cngain" ]]; then validate_input "$sample" "$cngain" cn_int_gain; fi
  if [[ -n "$cnloss" ]]; then validate_input "$sample" "$cnloss" cn_int_loss; fi
  record_status "$sample" preflight PASS "cn_source=$cnsource"

  {
    printf 'role\tpath\tsize_bytes\tmtime\n'
    for spec in "tumor_bam:$tbam" "tumor_bam_index:${tbam}.bai" "tree_input_vcf:$treevcf" \
                "longphase_input_vcf:$lpsinput" \
                "caller_raw_vcf:$rawvcf" "longphase_recalibrated_all_vcf:$recalall" \
                "read_tag_sidecar:$readtags" "read_tag_sidecar_index:$readtagsidx" \
                "read_tag_validation:$tagvalidation" "longphase_input_inventory:$lpsinventory" \
                "cn_bed:$cnbed" "cn_int_gain:$cngain" "cn_int_loss:$cnloss"; do
      local role="${spec%%:*}" path="${spec#*:}"
      [[ -n "$path" && -f "$path" ]] || continue
      printf '%s\t%s\t%s\t%s\n' "$role" "$path" "$(stat -c %s "$path")" "$(stat -c %y "$path")"
    done
  } > "$wd/input_files.tsv"
  {
    sha256sum "$treevcf" "$lpsinput" "$rawvcf" "$recalall" "$readtags" "$readtagsidx" \
      "$tagvalidation" "$lpsinventory" "${tbam}.bai"
    if [[ -f "${treevcf}.csi" ]]; then sha256sum "${treevcf}.csi"; fi
    if [[ -f "${treevcf}.tbi" ]]; then sha256sum "${treevcf}.tbi"; fi
    if [[ -f "${lpsinput}.csi" ]]; then sha256sum "${lpsinput}.csi"; fi
    if [[ -f "${lpsinput}.tbi" ]]; then sha256sum "${lpsinput}.tbi"; fi
    if [[ -n "$cnbed" ]]; then sha256sum "$cnbed"; fi
    if [[ -n "$cngain" ]]; then sha256sum "$cngain"; fi
    if [[ -n "$cnloss" ]]; then sha256sum "$cnloss"; fi
  } > "$wd/input.sha256"

  local start stage1_end end failed=0
  start="$(date +%s)"
  record_status "$sample" mlhp START "5 chromosome splits"
  local pids=()
  for split in "${splits[@]}"; do
    local chroms="${split%:*}" part="${split#*:}"
    env SM_WORKDIR="$wd" SM_TBAM="$tbam" SM_VD="$(dirname "$treevcf")" \
        SM_SOMATIC_VCF="$treevcf" SM_CNBED="$cnbed" SM_CN_SOURCE="$cnsource" \
        SM_READ_TAG_SIDECAR="$readtags" \
        SM_MINREAD="${SM_MINREAD:-3}" SM_MAX_SNV="${SM_MAX_SNV:-8}" \
        SM_TIER_R="${SM_TIER_R:-50000}" SM_MAPQ_MIN="${SM_MAPQ_MIN:-20}" \
        SM_BASEQ_MIN="${SM_BASEQ_MIN:-0}" \
        "$TIME" -v "$PY" "$EXEC_SCD/sm_multilocus_combinations.py" "$chroms" "$wd/mlhp_part_${part}.json" \
        > "$wd/mlhp_part_${part}.log" 2>&1 &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then failed=1; fi
  done
  if [[ "$failed" -ne 0 ]]; then
    record_status "$sample" mlhp FAIL "one or more chromosome splits failed"
    return 1
  fi
  for part in 1 2 3 4 5; do
    jq -e '.schema_version == "2.0" and .input_funnel.check_scope_conservation == true
      and .params.read_tag_source != "BAM_HP_PS"
      and .read_tag_census.sidecar_missing == 0
      and .read_tag_census.sidecar_conflicts == 0
      and .read_tag_census.sidecar_exact_matches == .read_tag_census.alignment_group_exposures' \
      "$wd/mlhp_part_${part}.json" >/dev/null
  done
  stage1_end="$(date +%s)"
  record_status "$sample" mlhp PASS "$((stage1_end-start))s"

  record_status "$sample" layered START "full V1-V7 for all non-capped units"
  env SM_ML_GLOB="$wd/mlhp_part_*.json" SM_OUT="$wd/layered_reconstruction_${sample}.json" \
      SM_VERIFY_EVERY="$VERIFY_EVERY" SM_ANALYSIS_TREE_CAP="$ANALYSIS_TREE_CAP" \
      SM_DISPLAY_TREE_CAP="$DISPLAY_TREE_CAP" SM_CN_INT_GAIN="$cngain" SM_CN_INT_LOSS="$cnloss" \
      "$TIME" -v "$PY" "$EXEC_SCD/layered_tree_reconstruction.py" > "$wd/layered.log" 2>&1
  jq -e '.schema_version == "2.0" and .L1_ssnv_algorithm.all_eligible_V1V7_pass == true' \
    "$wd/layered_reconstruction_${sample}.json" >/dev/null
  record_status "$sample" layered PASS "$(($(date +%s)-stage1_end))s"

  record_status "$sample" region_view START
  env SM_LAYERED="$wd/layered_reconstruction_${sample}.json" \
      SM_OUT="$wd/layered_region_view_${sample}.json" SM_SAMPLE="$sample" \
      SM_ML_GLOB="$wd/mlhp_part_*.json" SM_SOMATIC_VCF="$treevcf" SM_INTEGRATION="$integration" \
      SM_BACKBONE_SOURCE="$backbone_source" \
      "$PY" "$EXEC_SCD/build_region_view.py" > "$wd/region_view.log" 2>&1
  jq -e '.schema_version == "2.0" and .census.funnel.check_six_branch_conservation == true' \
    "$wd/layered_region_view_${sample}.json" >/dev/null
  end="$(date +%s)"
  record_status "$sample" region_view PASS "$((end-stage1_end))s"

  record_status "$sample" site_ledger START "all raw ClairS records"
  "$PY" "$EXEC_SCD/build_ssnv_site_ledger.py" --sample "$sample" --caller-raw-vcf "$rawvcf" \
      --longphase-input-vcf "$lpsinput" --tree-input-vcf "$treevcf" --recalibrated-vcf "$recalall" \
      --tree-contract "$TREE_CONTRACT" \
      --mlhp-glob "$wd/mlhp_part_*.json" --output-tsv-gz "$wd/ssnv_site_ledger_${sample}.tsv.gz" \
      --output-summary "$wd/ssnv_site_ledger_${sample}.summary.json" > "$wd/site_ledger.log" 2>&1
  jq -e '.pass == true' "$wd/ssnv_site_ledger_${sample}.summary.json" >/dev/null
  [[ -s "$wd/ssnv_site_ledger_${sample}.tsv.gz.tbi" ]]
  record_status "$sample" site_ledger PASS "rows=$(jq -r '.raw_clairs_records' "$wd/ssnv_site_ledger_${sample}.summary.json")"
  end="$(date +%s)"

  sha256sum "$wd"/mlhp_part_*.json "$wd/layered_reconstruction_${sample}.json" \
    "$wd/layered_region_view_${sample}.json" "$wd/ssnv_site_ledger_${sample}.tsv.gz" \
    "$wd/ssnv_site_ledger_${sample}.tsv.gz.tbi" \
    "$wd/ssnv_site_ledger_${sample}.summary.json" > "$wd/output.sha256"
  record_status "$sample" complete PASS "$((end-start))s"
  echo "[$sample] DONE total=$((end-start))s mlhp=$((stage1_end-start))s downstream=$((end-stage1_end))s"
}

export -f record_status manifest_value validate_input run_sample
export SCD EXEC_SCD PY MANIFEST RUN_ROOT VERIFY_EVERY ANALYSIS_TREE_CAP DISPLAY_TREE_CAP TREE_CONTRACT TIME PROFILE STATUS
export SM_MINREAD="${SM_MINREAD:-3}" SM_MAX_SNV="${SM_MAX_SNV:-8}" SM_TIER_R="${SM_TIER_R:-50000}"
export SM_MAPQ_MIN="${SM_MAPQ_MIN:-20}" SM_BASEQ_MIN="${SM_BASEQ_MIN:-0}"

{
  sha256sum "$MANIFEST" "$EXEC_SCD/sm_linkage_genomewide.py" "$EXEC_SCD/sm_multilocus_combinations.py" \
    "$EXEC_SCD/tree_enumeration_solver.py" "$EXEC_SCD/layered_tree_reconstruction.py" \
    "$EXEC_SCD/build_region_view.py" "$EXEC_SCD/build_ssnv_site_ledger.py" \
    "$EXEC_SCD/verify_layered_v2.py" "$EXEC_SCD/validate_layered_v2_inputs.py" \
    "$EXEC_SCD/run_layered_7samples_newbb.sh"
} > "$RUN_ROOT/code.sha256"
{
  printf 'timestamp=%s\n' "$(date --iso-8601=seconds)"
  printf 'host=%s\n' "$(hostname)"
  printf 'kernel=%s\n' "$(uname -srvmo)"
  printf 'python='; "$PY" --version 2>&1
  printf 'python_path=%s\n' "$PY"
  printf 'pysam='; "$PY" -c 'import pysam; print(pysam.__version__)'
  printf 'samtools='; samtools --version-only
  printf 'bcftools='; bcftools --version-only
  printf 'jq='; jq --version
} > "$RUN_ROOT/environment.txt"
jq -n --arg run_id "$RUN_ID" --arg run_root "$RUN_ROOT" --arg python "$PY" \
  --arg dataset_count "$DATASET_COUNT" --arg biological_count "$BIOLOGICAL_COUNT" \
  --arg verify_every "$VERIFY_EVERY" --arg analysis_cap "$ANALYSIS_TREE_CAP" \
  --arg display_cap "$DISPLAY_TREE_CAP" --arg minread "$SM_MINREAD" \
  --arg max_snv "$SM_MAX_SNV" --arg tier_r "$SM_TIER_R" --arg mapq "$SM_MAPQ_MIN" \
  --arg baseq "$SM_BASEQ_MIN" \
  --arg tree_contract "$TREE_CONTRACT" \
  '{run_id:$run_id,run_root:$run_root,python:$python,scope:"chr1-22",dataset_count:($dataset_count|tonumber),biological_sample_count:($biological_count|tonumber),
    tree_input_contract:$tree_contract,
    params:{VERIFY_EVERY:($verify_every|tonumber),ANALYSIS_TREE_CAP:($analysis_cap|tonumber),DISPLAY_TREE_CAP:($display_cap|tonumber),
            MINREAD:($minread|tonumber),MAX_SNV:($max_snv|tonumber),TIER_R:($tier_r|tonumber),MAPQ_MIN:($mapq|tonumber),BASEQ_MIN:($baseq|tonumber)}}' \
  > "$RUN_ROOT/params.json"
jq -n --arg timestamp "$(date --iso-8601=seconds)" --arg source_manifest "$SOURCE_MANIFEST" \
  --arg frozen_manifest "$MANIFEST" --arg manifest_sha "$SOURCE_MANIFEST_SHA" \
  --arg source_bundle_sha "$(sha256sum "$RUN_ROOT/source_bundle.sha256" | cut -d' ' -f1)" \
  --arg command "$0" --arg run_root "$RUN_ROOT" \
  '{timestamp:$timestamp,source_manifest:$source_manifest,frozen_manifest:$frozen_manifest,
    manifest_sha256:$manifest_sha,source_bundle_manifest_sha256:$source_bundle_sha,
    command:$command,run_root:$run_root,workers_read_frozen_manifest:true,workers_read_source_bundle:true}' \
  > "$RUN_ROOT/launch_receipt.json"

echo "=== LAYERED V2 START $(date --iso-8601=seconds) ==="
echo "INPUT MANIFEST: $MANIFEST"
echo "OUTPUT ROOT: $RUN_ROOT"
mapfile -t SAMPLES < <(jq -r '.samples[].sample' "$MANIFEST")
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$PARALLEL_SAMPLES" -I {} \
  bash -c 'set -euo pipefail; run_sample "$1"' _ {}

"$PY" "$EXEC_SCD/verify_layered_v2.py" --run-root "$RUN_ROOT" \
  --input-manifest "$MANIFEST" --output "$RUN_ROOT/verification_summary.json"
sha256sum "$RUN_ROOT/verification_summary.json" "$RUN_ROOT/verification_summary.tsv" \
  > "$RUN_ROOT/verification.sha256"
record_status ALL verify PASS "$DATASET_COUNT/$DATASET_COUNT datasets"
echo "=== LAYERED V2 COMPLETE $(date --iso-8601=seconds) ==="
echo "OUTPUT ROOT: $RUN_ROOT"
cat "$RUN_ROOT/verification_summary.tsv"
jq -n --arg timestamp "$(date --iso-8601=seconds)" \
  --arg verification_sha "$(sha256sum "$RUN_ROOT/verification_summary.json" | cut -d' ' -f1)" \
  --argjson dataset_count "$DATASET_COUNT" \
  '{status:"SUCCESS",timestamp:$timestamp,dataset_count:$dataset_count,
    verification_sha256:$verification_sha}' > "$RUN_ROOT/_SUCCESS.pending.$$"
mv "$RUN_ROOT/_SUCCESS.pending.$$" "$RUN_ROOT/_SUCCESS"
RUN_COMPLETED=1
