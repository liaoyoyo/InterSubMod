#!/usr/bin/env bash
set -euo pipefail

# Overlap strict region construction with the tail of the interrupted-extraction
# resume. A sample is released only after all 22 extraction receipt sidecars and
# pinned manifest provenance pass. Strict workers are capped at P2.

REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
RUN_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1"
SAMPLE_ROOT="$RUN_ROOT/samples"
MANIFEST="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/release_contract/input_contract/canonical_manifest.json"
MANIFEST_SHA256="16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
BUILD_SCRIPT="$REPO_ROOT/scripts/build_strict_ps_hp_regions.py"
BUILD_SCRIPT_SHA256="912721f934bae6a58ccfe66b872706f035b658d1fc1b53db4025a998916e4b4d"
SUMMARY_SCRIPT="$REPO_ROOT/scripts/summarize_strict_ps_hp_regions.py"
SUMMARY_SCRIPT_SHA256="159d60547cb632c5ab58ef173e90cde8e49001017b188eb870ffee9062bec385"
GRAPH_CORE="$REPO_ROOT/tools/strict_endpoint_graph.py"
GRAPH_CORE_SHA256="df3d6f37615b0bfa5d382a082152cdb16ce206eaa075cd3c0dd418b5c320e37e"
RESUME_PID="${RESUME_PID:?set RESUME_PID to the active extraction-resume PID}"
RUN_ID="${RESUME_RUN_ID:?set RESUME_RUN_ID to the extraction resume run id}"
SAMPLES="${STRICT_SAMPLES:-HCC1937 HCC1954}"
AUDIT_ROOT="$RUN_ROOT/resume_audit/$RUN_ID"
LOG_PATH="$AUDIT_ROOT/strict_overlap.log"
RESUME_LOG="$AUDIT_ROOT/resume.log"
AUTOSOMES="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

mkdir -p "$AUDIT_ROOT"
exec > >(tee -a "$LOG_PATH" "$RESUME_LOG") 2>&1

sha256_of() {
    sha256sum "$1" | awk '{print $1}'
}

verify_extraction_receipt() {
    local sample="$1"
    local chrom="$2"
    local output_dir="$SAMPLE_ROOT/$sample/chromosomes/$chrom/extraction"
    local receipt="$output_dir/receipt.json"
    local sidecar="$output_dir/receipt.json.sha256"
    [[ -f "$receipt" && -f "$sidecar" ]] || return 1
    jq -e \
        --arg sample "$sample" \
        --arg chrom "$chrom" \
        --arg manifest "$MANIFEST" \
        --arg manifest_sha "$MANIFEST_SHA256" \
        '.all_pass == true
         and .scope.dataset == $sample
         and .scope.chrom == $chrom
         and .provenance.manifest.path == $manifest
         and .provenance.manifest.sha256 == $manifest_sha' \
        "$receipt" >/dev/null
    (cd "$output_dir" && sha256sum -c receipt.json.sha256 >/dev/null)
}

verify_all_extractions() {
    local sample="$1"
    local index
    for index in $(seq 1 22); do
        verify_extraction_receipt "$sample" "chr$index" || return 1
    done
}

verify_strict_receipt() {
    local sample="$1"
    local chrom="$2"
    local output_dir="$SAMPLE_ROOT/$sample/strict_regions_v1/chromosomes/$chrom"
    local receipt="$output_dir/receipt.json"
    local sidecar="$output_dir/receipt.json.sha256"
    [[ -f "$receipt" && -f "$sidecar" ]] || return 1
    jq -e \
        --arg sample "$sample" \
        --arg chrom "$chrom" \
        --arg builder "$BUILD_SCRIPT" \
        --arg builder_sha "$BUILD_SCRIPT_SHA256" \
        --arg graph_core "$GRAPH_CORE" \
        --arg graph_sha "$GRAPH_CORE_SHA256" \
        '.schema_name == "intersubmod.strict_exact_ps_hp_regions"
         and .schema_version == "1.1.0"
         and .all_pass == true
         and .scope.dataset == $sample
         and .scope.chrom == $chrom
         and .parameters.primary_threshold == 3
         and .parameters.thresholds == [1,2,3,5]
         and .parameters.linkage_rule == "strict_fixed_ra_endpoint_pair"
         and .inputs.builder.path == $builder
         and .inputs.builder.sha256 == $builder_sha
         and .inputs.strict_graph_core.path == $graph_core
         and .inputs.strict_graph_core.sha256 == $graph_sha' \
        "$receipt" >/dev/null
    (cd "$output_dir" && sha256sum -c receipt.json.sha256 >/dev/null)
}

run_strict_one() {
    local sample="$1"
    local chrom="$2"
    local extraction_dir="$SAMPLE_ROOT/$sample/chromosomes/$chrom/extraction"
    local output_dir="$SAMPLE_ROOT/$sample/strict_regions_v1/chromosomes/$chrom"
    verify_extraction_receipt "$sample" "$chrom" || {
        echo "FAIL strict worker received invalid extraction sample=$sample chrom=$chrom" >&2
        return 20
    }
    if verify_strict_receipt "$sample" "$chrom"; then
        echo "STRICT_SKIP_VALID sample=$sample chrom=$chrom"
        return 0
    fi
    [[ ! -e "$output_dir" ]] || {
        echo "FAIL strict output exists without a valid receipt: $output_dir" >&2
        return 21
    }
    mkdir -p "$output_dir"
    echo "STRICT_RUN sample=$sample chrom=$chrom output=$output_dir"
    "$PYTHON" "$BUILD_SCRIPT" \
        --dataset "$sample" \
        --chrom "$chrom" \
        --site-catalog "$extraction_dir/$sample.$chrom.site_catalog.tsv.gz" \
        --molecule-calls "$extraction_dir/$sample.$chrom.molecule_sparse_calls.tsv.gz" \
        --output-dir "$output_dir" \
        --thresholds 1,2,3,5 \
        --primary-threshold 3 \
        --require-existing-empty-output-dir
    verify_strict_receipt "$sample" "$chrom" || {
        echo "FAIL strict post-run receipt sample=$sample chrom=$chrom" >&2
        return 22
    }
    echo "STRICT_PASS sample=$sample chrom=$chrom receipt=$output_dir/receipt.json"
}
export -f verify_extraction_receipt verify_strict_receipt run_strict_one
export SAMPLE_ROOT MANIFEST MANIFEST_SHA256 PYTHON BUILD_SCRIPT BUILD_SCRIPT_SHA256
export GRAPH_CORE GRAPH_CORE_SHA256

echo "STRICT_OVERLAP_START pid=$$ extraction_resume_pid=$RESUME_PID run_id=$RUN_ID"
echo "COMMAND=RESUME_RUN_ID=$RUN_ID RESUME_PID=$RESUME_PID $REPO_ROOT/research/20260723_production_exact_ps_strict_read_linkage/scripts/run_strict_overlap_after_resume.sh"
echo "SAMPLES=$SAMPLES"
echo "BUILD_SCRIPT=$BUILD_SCRIPT sha256=$BUILD_SCRIPT_SHA256"
echo "SUMMARY_SCRIPT=$SUMMARY_SCRIPT sha256=$SUMMARY_SCRIPT_SHA256"
echo "GRAPH_CORE=$GRAPH_CORE sha256=$GRAPH_CORE_SHA256"
echo "POLICY=sample_gate_22_of_22_then_strict_P2_no_partial_no_overwrite"

[[ "$(sha256_of "$MANIFEST")" == "$MANIFEST_SHA256" ]] || exit 30
[[ "$(sha256_of "$BUILD_SCRIPT")" == "$BUILD_SCRIPT_SHA256" ]] || exit 31
[[ "$(sha256_of "$SUMMARY_SCRIPT")" == "$SUMMARY_SCRIPT_SHA256" ]] || exit 32
[[ "$(sha256_of "$GRAPH_CORE")" == "$GRAPH_CORE_SHA256" ]] || exit 33

for sample in $SAMPLES; do
    while ! verify_all_extractions "$sample"; do
        receipt_count=$(
            { find "$SAMPLE_ROOT/$sample" -path '*/extraction/receipt.json' -type f 2>/dev/null || true; } \
                | wc -l
        )
        if ! ps -p "$RESUME_PID" >/dev/null 2>&1; then
            echo "FAIL extraction resume ended before $sample reached 22/22; observed=$receipt_count" >&2
            exit 34
        fi
        echo "STRICT_WAIT sample=$sample extraction_receipts=$receipt_count/22"
        sleep 30
    done
    echo "STRICT_GATE_PASS sample=$sample extraction_receipts=22/22"

    strict_root="$SAMPLE_ROOT/$sample/strict_regions_v1"
    if [[ -e "$strict_root" ]]; then
        for index in $(seq 1 22); do
            chrom="chr$index"
            output_dir="$strict_root/chromosomes/$chrom"
            if [[ -e "$output_dir" ]] && ! verify_strict_receipt "$sample" "$chrom"; then
                echo "FAIL non-valid pre-existing strict output: $output_dir" >&2
                exit 35
            fi
        done
        if [[ -e "$strict_root/summary" ]]; then
            echo "FAIL pre-existing strict summary is not permitted: $strict_root/summary" >&2
            exit 36
        fi
    fi

    echo "STRICT_P2_COMMAND sample=$sample command=seq_1_22_pipe_xargs_-P2_build_strict_ps_hp_regions"
    (
        for index in $(seq 1 22); do
            printf '%s %s\n' "$sample" "chr$index"
        done
    ) | xargs -P 2 -n 2 bash -c 'run_strict_one "$1" "$2"' _ &
    strict_xargs_pid=$!
    echo "STRICT_P2_PID sample=$sample xargs_pid=$strict_xargs_pid"
    wait "$strict_xargs_pid"

    for index in $(seq 1 22); do
        verify_strict_receipt "$sample" "chr$index" || {
            echo "FAIL final strict receipt sample=$sample chrom=chr$index" >&2
            exit 37
        }
    done
    mkdir "$strict_root/summary"
    echo "SUMMARY_RUN sample=$sample command=$PYTHON $SUMMARY_SCRIPT --dataset $sample --input-root $strict_root --output-dir $strict_root/summary --chromosomes $AUTOSOMES --primary-threshold 3"
    "$PYTHON" "$SUMMARY_SCRIPT" \
        --dataset "$sample" \
        --input-root "$strict_root" \
        --output-dir "$strict_root/summary" \
        --chromosomes "$AUTOSOMES" \
        --primary-threshold 3
    jq -e \
        --arg sample "$sample" \
        '.all_pass == true
         and .scope.dataset == $sample
         and .scope.primary_threshold == 3
         and (.checks | to_entries | all(.value == true))' \
        "$strict_root/summary/summary.json" >/dev/null
    echo "SUMMARY_PASS sample=$sample path=$strict_root/summary/summary.json"
done
echo "STRICT_OVERLAP_COMPLETE pid=$$"
