#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  run_tiny_public_e2e.sh [--repo-root PATH] [--work-dir NEW_PATH]
                         [--build-dir PATH | --binary PATH] [--jobs N]
  run_tiny_public_e2e.sh --source-repository URL_OR_PATH --revision REF
                         [--work-dir NEW_PATH] [--jobs N]

The source-repository mode performs clone -> detached checkout -> build -> run -> validation.
The repo-root mode validates an existing checkout with isolated build and run directories.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_repo_root="$(cd "$script_dir/../.." && pwd)"
repo_root="$default_repo_root"
repo_root_explicit=false
source_repository=""
revision=""
work_dir=""
build_dir=""
binary=""
jobs="4"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-root) repo_root="$2"; repo_root_explicit=true; shift 2 ;;
        --source-repository) source_repository="$2"; shift 2 ;;
        --revision) revision="$2"; shift 2 ;;
        --work-dir) work_dir="$2"; shift 2 ;;
        --build-dir) build_dir="$2"; shift 2 ;;
        --binary) binary="$2"; shift 2 ;;
        --jobs) jobs="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

if [[ -n "$source_repository" && -z "$revision" ]]; then
    echo "--revision is required with --source-repository" >&2
    exit 2
fi
if [[ "$source_repository" == *"://"* ]]; then
    repository_authority="${source_repository#*://}"
    repository_authority="${repository_authority%%/*}"
    if [[ "$repository_authority" == *"@"* || "$source_repository" == *"?"* || "$source_repository" == *"#"* ]]; then
        echo "--source-repository must not contain URL credentials, query parameters, or fragments" >&2
        exit 2
    fi
fi
if [[ -n "$source_repository" && "$repo_root_explicit" == true ]]; then
    echo "--repo-root cannot be combined with --source-repository" >&2
    exit 2
fi
if [[ -n "$source_repository" && ( -n "$binary" || -n "$build_dir" ) ]]; then
    echo "--source-repository owns checkout and build identity; --binary/--build-dir are forbidden" >&2
    exit 2
fi
if [[ -z "$source_repository" && -n "$revision" ]]; then
    echo "--revision is valid only with --source-repository" >&2
    exit 2
fi
if [[ -n "$binary" && -n "$build_dir" ]]; then
    echo "Use either --binary or --build-dir, not both" >&2
    exit 2
fi
if ! [[ "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "--jobs must be a positive integer" >&2
    exit 2
fi

if [[ -z "$work_dir" ]]; then
    work_dir="$(mktemp -d /tmp/intersubmod-tiny-e2e.XXXXXX)"
elif [[ -e "$work_dir" ]]; then
    echo "Refusing to reuse existing --work-dir: $work_dir" >&2
    exit 3
else
    mkdir -p "$work_dir"
fi

if [[ -n "$source_repository" ]]; then
    checkout="$work_dir/checkout"
    git clone --no-checkout "$source_repository" "$checkout"
    git -C "$checkout" checkout --detach "$revision"
    repo_root="$checkout"
    resolved_commit="$(git -C "$checkout" rev-parse HEAD)"
else
    resolved_commit="$(git -C "$repo_root" rev-parse HEAD 2>/dev/null || true)"
fi

repo_root="$(cd "$repo_root" && pwd)"
fixture_dir="$work_dir/fixture"
output_dir="$work_dir/output"
run_log="$work_dir/inter_sub_mod.log"
command_receipt="$work_dir/run_command.txt"
receipt="$work_dir/tiny_public_e2e_receipt.json"
schema="$repo_root/tests/fixtures/tiny_public/expected_significance_schema.json"

if [[ -z "$binary" ]]; then
    if [[ -z "$build_dir" ]]; then
        build_dir="$work_dir/build"
    fi
    cmake -S "$repo_root" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
    cmake --build "$build_dir" --target inter_sub_mod -j "$jobs"
    binary="$build_dir/bin/inter_sub_mod"
fi
[[ -x "$binary" ]] || { echo "InterSubMod binary is not executable: $binary" >&2; exit 4; }

"$repo_root/scripts/handoff/build_tiny_public_fixture.sh" \
    --repo-root "$repo_root" \
    --output-dir "$fixture_dir"

run_command=(
    "$binary"
    --tumor-bam "$fixture_dir/tumor.bam"
    --reference "$fixture_dir/reference.fa"
    --vcf "$fixture_dir/variants.vcf"
    --output-dir "$output_dir"
    --window-size 1000
    --threads 1
    --min-read-length 1
    --min-common-coverage 3
    --distance-metric NHD
    --nan-distance-strategy SKIP
    --expected-coverage 12
)
printf '%q ' "${run_command[@]}" > "$command_receipt"
printf '\n' >> "$command_receipt"
set +e
"${run_command[@]}" > "$run_log" 2>&1
run_exit_code=$?
set -e
printf 'exit_code=%s\n' "$run_exit_code" >> "$command_receipt"

validator_args=(
    --repo-root "$repo_root"
    --fixture-dir "$fixture_dir"
    --output-dir "$output_dir"
    --schema "$schema"
    --receipt "$receipt"
    --binary "$binary"
    --run-log "$run_log"
    --command-receipt "$command_receipt"
    --run-exit-code "$run_exit_code"
)
if [[ -n "$source_repository" ]]; then
    validator_args+=(
        --source-repository "$source_repository"
        --requested-revision "$revision"
        --resolved-commit "$resolved_commit"
    )
elif [[ -n "$resolved_commit" ]]; then
    validator_args+=(--resolved-commit "$resolved_commit")
fi
python3 "$repo_root/scripts/handoff/validate_tiny_public_e2e.py" "${validator_args[@]}"

python3 - "$receipt" <<'PY'
import json, sys
receipt = json.load(open(sys.argv[1], encoding="utf-8"))
print("TINY_E2E_RESULT", json.dumps({
    "all_pass": receipt["all_pass"],
    "clean_clone_e2e_pass": receipt["clean_clone_e2e_pass"],
    "overall_release_gate": receipt["overall_release_gate"],
    "column_count": receipt["summary"]["column_count"],
    "tree_leaf_count": receipt["tree"]["leaf_count"],
    "tree_semantics": receipt["tree"]["semantics"],
}, sort_keys=True))
PY
echo "TINY_E2E_OUTPUT work_dir=$work_dir output=$output_dir receipt=$receipt log=$run_log"
