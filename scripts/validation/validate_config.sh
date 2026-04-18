#!/bin/bash
# validate_config.sh - Validation framework configuration
#
# Centralized thresholds, paths, and sample definitions for the
# automated single-change validation framework.
# Source this from validate.sh and runner scripts.

set -euo pipefail

# ============================================================================
# Locate project root and source pipeline config
# ============================================================================

VALIDATION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${VALIDATION_DIR}/../.." && pwd)"
PIPELINE_DIR="${PROJECT_ROOT}/scripts/pipeline"

# Source the canonical pipeline config (sample defs, paths, helpers)
source "${PIPELINE_DIR}/config.sh"

# ============================================================================
# Experiment Output
# ============================================================================

EXPERIMENT_OUTPUT_ROOT="${EXPERIMENT_OUTPUT_ROOT:-${BIG7_OUTPUT_ROOT}/experiments}"
BASELINE_STORE="${VALIDATION_DIR}/baseline/baseline_store.json"
EVIDENCE_LEDGER="${PROJECT_ROOT}/research/autoresearch/evidence_ledger.jsonl"
HYPOTHESIS_QUEUE="${PROJECT_ROOT}/research/autoresearch/hypothesis_queue.json"

# ============================================================================
# Sample Definitions
# ============================================================================

# All 7 samples for full validation
ALL_SAMPLES=("HCC1395" "HCC1395_DORADO" "COLO829" "H1437" "H2009" "HCC1937" "HCC1954")

# Quick validation defaults
QUICK_SAMPLE="${QUICK_SAMPLE:-HCC1395}"
QUICK_MODE="${QUICK_MODE:-s-pure-pileup}"

# Full validation defaults
FULL_MODE="${FULL_MODE:-s-pure-pileup}"
MAX_PARALLEL_VALIDATE="${MAX_PARALLEL_VALIDATE:-4}"

# Track presets: paired vs TO
# --track paired → pipeline-mode s-pure-pileup (default)
# --track to     → pipeline-mode to-pure
declare -A TRACK_PIPELINE_MODE=(
    [paired]="s-pure-pileup"
    [to]="to-pure"
)
declare -A TRACK_CANONICAL_MODE=(
    [paired]="paired_pileup"
    [to]="to_pileup"
)

# ============================================================================
# Thresholds
# ============================================================================

# F1 delta thresholds
F1_REGRESSION_THRESHOLD="${F1_REGRESSION_THRESHOLD:-0.005}"   # Alert if any sample drops > this
F1_IMPROVEMENT_THRESHOLD="${F1_IMPROVEMENT_THRESHOLD:-0.001}" # Minimum meaningful improvement

# Statistical test parameters
BOOTSTRAP_ITERATIONS="${BOOTSTRAP_ITERATIONS:-1000}"
BOOTSTRAP_CI_LEVEL="${BOOTSTRAP_CI_LEVEL:-0.95}"
SIGNIFICANCE_ALPHA="${SIGNIFICANCE_ALPHA:-0.05}"

# ============================================================================
# Resource Limits
# ============================================================================

# Disk space (in GB) required before starting experiments
VALIDATE_MIN_FREE_GB="${VALIDATE_MIN_FREE_GB:-100}"
# Estimated output per sample (GB)
ESTIMATED_GB_PER_SAMPLE=5
# Per-sample timeout (seconds) — 30 minutes
SAMPLE_TIMEOUT="${SAMPLE_TIMEOUT:-1800}"

# ============================================================================
# Build Configuration
# ============================================================================

BUILD_DIR="${PROJECT_ROOT}/build"
BUILD_TYPE="${BUILD_TYPE:-Release}"

# ============================================================================
# Utility Functions
# ============================================================================

# Generate unique experiment ID
generate_experiment_id() {
    local date_str
    date_str="$(date +%Y%m%d_%H%M%S)"
    local git_hash
    git_hash="$(cd "${PROJECT_ROOT}" && git rev-parse --short HEAD 2>/dev/null || echo "nohash")"
    echo "EXP_${date_str}_${git_hash}"
}

# Get current git state as JSON-compatible fields
get_git_state() {
    local git_dir="${PROJECT_ROOT}"
    local commit branch dirty diff_stat

    commit="$(cd "${git_dir}" && git rev-parse --short HEAD 2>/dev/null || echo "unknown")"
    branch="$(cd "${git_dir}" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")"
    dirty="$(cd "${git_dir}" && git diff --quiet 2>/dev/null && echo "false" || echo "true")"
    diff_stat="$(cd "${git_dir}" && git diff --stat HEAD 2>/dev/null | tail -1 || echo "")"

    echo "COMMIT=${commit}"
    echo "BRANCH=${branch}"
    echo "DIRTY=${dirty}"
    echo "DIFF_STAT=${diff_stat}"
}

# Check if baseline store exists
has_baseline() {
    [[ -f "${BASELINE_STORE}" ]]
}

log_validate() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [VALIDATE] $*" >&2
}
