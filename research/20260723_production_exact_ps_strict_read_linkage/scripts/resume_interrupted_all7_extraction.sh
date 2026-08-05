#!/usr/bin/env bash
set -euo pipefail

# Resume only the 43 chromosome extractions left incomplete by the interrupted
# all-7 production run. Existing valid receipts are verified and skipped;
# unexpected partial outputs fail closed. The four known interrupted HCC1937
# directories are archived rather than deleted.

REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
RUN_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/all7_production_v1"
SAMPLE_ROOT="$RUN_ROOT/samples"
MANIFEST="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/release_contract/input_contract/canonical_manifest.json"
MANIFEST_SHA256="16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
EXTRACTOR="$REPO_ROOT/research/20260718_k_gt8_read_supported_segmentation/scripts/extract_lossless_read_linkage_collapsing.py"
EXTRACTOR_SHA256="2ca7ccb67c89e816fae9284f4e2ba21b186378105086c6b0128ed5445a133490"
PYTHON="/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
SAMTOOLS="/usr/local/bin/samtools"
RUN_ID="${RESUME_RUN_ID:?set RESUME_RUN_ID to a unique timestamp-like identifier}"
AUDIT_ROOT="$RUN_ROOT/resume_audit/$RUN_ID"
ARCHIVE_ROOT="$RUN_ROOT/interrupted_archive/$RUN_ID"
LOG_PATH="$AUDIT_ROOT/resume.log"

mkdir -p "$AUDIT_ROOT" "$ARCHIVE_ROOT"
exec > >(tee -a "$LOG_PATH") 2>&1

sha256_of() {
    sha256sum "$1" | awk '{print $1}'
}

verify_receipt() {
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

echo "RESUME_START run_id=$RUN_ID"
echo "SCRIPT=$REPO_ROOT/research/20260723_production_exact_ps_strict_read_linkage/scripts/resume_interrupted_all7_extraction.sh"
echo "MANIFEST=$MANIFEST"
echo "MANIFEST_SHA256=$MANIFEST_SHA256"
echo "EXTRACTOR=$EXTRACTOR"
echo "EXTRACTOR_SHA256=$EXTRACTOR_SHA256"
echo "OUTPUT_ROOT=$SAMPLE_ROOT"
echo "ARCHIVE_ROOT=$ARCHIVE_ROOT"
echo "LOG_PATH=$LOG_PATH"

[[ "$(sha256_of "$MANIFEST")" == "$MANIFEST_SHA256" ]] || {
    echo "FAIL manifest SHA-256 mismatch" >&2
    exit 10
}
[[ "$(sha256_of "$EXTRACTOR")" == "$EXTRACTOR_SHA256" ]] || {
    echo "FAIL extractor SHA-256 mismatch" >&2
    exit 11
}

# The completed HCC1937 chr4 receipt is the provenance anchor for this resume.
verify_receipt HCC1937 chr4 || {
    echo "FAIL HCC1937 chr4 provenance anchor is invalid" >&2
    exit 12
}
echo "ANCHOR_PASS sample=HCC1937 chrom=chr4"

declare -a TASKS=()
for index in 1 2 3; do
    TASKS+=("HCC1937 chr$index")
done
for index in $(seq 5 22); do
    TASKS+=("HCC1937 chr$index")
done
for index in $(seq 1 22); do
    TASKS+=("HCC1954 chr$index")
done
[[ "${#TASKS[@]}" -eq 43 ]] || {
    echo "FAIL internal task count is ${#TASKS[@]}, expected 43" >&2
    exit 13
}

# Archive only the four partial directories observed after SIGTERM.
for chrom in chr1 chr2 chr3 chr5; do
    output_dir="$SAMPLE_ROOT/HCC1937/chromosomes/$chrom/extraction"
    archive_dir="$ARCHIVE_ROOT/HCC1937/chromosomes/$chrom/extraction"
    if verify_receipt HCC1937 "$chrom"; then
        echo "SKIP_ARCHIVE_VALID sample=HCC1937 chrom=$chrom"
        continue
    fi
    if [[ -d "$output_dir" ]]; then
        [[ -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit)" ]] || {
            echo "FAIL expected interrupted directory is empty: $output_dir" >&2
            exit 14
        }
        [[ ! -e "$archive_dir" ]] || {
            echo "FAIL archive destination already exists: $archive_dir" >&2
            exit 15
        }
        mkdir -p "$(dirname "$archive_dir")"
        mv "$output_dir" "$archive_dir"
        echo "ARCHIVED_PARTIAL sample=HCC1937 chrom=$chrom from=$output_dir to=$archive_dir"
    fi
done

# Fail closed on any unexpected nonempty target without a valid receipt.
for task in "${TASKS[@]}"; do
    read -r sample chrom <<<"$task"
    output_dir="$SAMPLE_ROOT/$sample/chromosomes/$chrom/extraction"
    if verify_receipt "$sample" "$chrom"; then
        echo "PREFLIGHT_SKIP_VALID sample=$sample chrom=$chrom"
    elif [[ -e "$output_dir" ]]; then
        echo "FAIL unexpected non-valid output target: $output_dir" >&2
        exit 16
    fi
done

run_one() {
    local sample="$1"
    local chrom="$2"
    local output_dir="$SAMPLE_ROOT/$sample/chromosomes/$chrom/extraction"
    if verify_receipt "$sample" "$chrom"; then
        echo "SKIP_VALID sample=$sample chrom=$chrom"
        return 0
    fi
    [[ ! -e "$output_dir" ]] || {
        echo "FAIL worker found existing non-valid target: $output_dir" >&2
        return 17
    }
    echo "RUN sample=$sample chrom=$chrom output=$output_dir"
    "$PYTHON" "$EXTRACTOR" \
        --manifest "$MANIFEST" \
        --sample "$sample" \
        --chrom "$chrom" \
        --output-dir "$output_dir" \
        --mapq-min 20 \
        --baseq-min 20 \
        --bridge-thresholds 1,2,3,5 \
        --samtools-threads 1 \
        --samtools "$SAMTOOLS"
    verify_receipt "$sample" "$chrom" || {
        echo "FAIL post-run receipt validation sample=$sample chrom=$chrom" >&2
        return 18
    }
    echo "PASS sample=$sample chrom=$chrom receipt=$output_dir/receipt.json"
}
export -f sha256_of verify_receipt run_one
export SAMPLE_ROOT MANIFEST MANIFEST_SHA256 EXTRACTOR PYTHON SAMTOOLS

printf '%s\n' "${TASKS[@]}" \
    | xargs -P 4 -n 2 bash -c 'run_one "$1" "$2"' _

for sample in HCC1937 HCC1954; do
    count=0
    for index in $(seq 1 22); do
        verify_receipt "$sample" "chr$index" || {
            echo "FAIL final receipt validation sample=$sample chrom=chr$index" >&2
            exit 19
        }
        count=$((count + 1))
    done
    echo "FINAL_PASS sample=$sample receipts=$count/22"
done
echo "RESUME_COMPLETE run_id=$RUN_ID"
