#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 --repo-root PATH --output-dir PATH" >&2
}

repo_root=""
output_dir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-root) repo_root="$2"; shift 2 ;;
        --output-dir) output_dir="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

if [[ -z "$repo_root" || -z "$output_dir" ]]; then
    usage
    exit 2
fi
if [[ -e "$output_dir" ]]; then
    echo "Refusing to overwrite existing fixture output: $output_dir" >&2
    exit 3
fi
command -v samtools >/dev/null 2>&1 || { echo "samtools is required" >&2; exit 4; }

source_dir="$repo_root/tests/fixtures/tiny_public"
for name in reference.fa variants.vcf tumor.sam; do
    [[ -f "$source_dir/$name" ]] || { echo "Missing fixture source: $source_dir/$name" >&2; exit 5; }
done

mkdir -p "$output_dir"
cp "$source_dir/reference.fa" "$output_dir/reference.fa"
cp "$source_dir/variants.vcf" "$output_dir/variants.vcf"

samtools faidx "$output_dir/reference.fa"
samtools view -u "$source_dir/tumor.sam" | samtools sort -o "$output_dir/tumor.bam"
samtools index "$output_dir/tumor.bam"

samtools quickcheck -v "$output_dir/tumor.bam"
read_count="$(samtools view -c "$output_dir/tumor.bam")"
if [[ "$read_count" != "12" ]]; then
    echo "Expected 12 synthetic reads, observed $read_count" >&2
    exit 6
fi

echo "FIXTURE_READY output=$output_dir reads=$read_count"
