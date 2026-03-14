#!/usr/bin/env python3
"""Build standardized benchmark comparison tables for research rounds."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

from research_common import (
    compute_metrics,
    infer_caller_name,
    infer_platform,
    infer_truth_set,
    markdown_table,
    parse_variant_counts,
    read_json,
    write_tsv_rows,
)


FIELDS = [
    "sample",
    "platform",
    "truth_set",
    "method",
    "caller",
    "mode",
    "truth_total",
    "calls_total",
    "tp",
    "fp",
    "fn",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="Source run directory with metrics.json and longphase_s/")
    parser.add_argument("--context-json", default=None, help="Optional round_context.json path")
    parser.add_argument("--output-tsv", required=True, help="Output benchmark_comparison.tsv")
    parser.add_argument("--output-md", default=None, help="Optional markdown summary path")
    parser.add_argument("--bcftools", default=shutil.which("bcftools") or "bcftools", help="bcftools binary")
    parser.add_argument(
        "--baseline-counts-file",
        default=None,
        help="Optional variant_counts.txt path for baseline method (default: run_dir/longphase_s/variant_counts.txt)",
    )
    parser.add_argument(
        "--baseline-method",
        default="LongPhase-S",
        help="Label for the baseline row (default: LongPhase-S)",
    )
    parser.add_argument(
        "--no-caller-input",
        action="store_true",
        help="Skip caller-input benchmark row even if somatic_vcf/truth are available",
    )
    parser.add_argument("--caller-name", default=None, help="Override caller name for input VCF benchmark")
    return parser.parse_args()


def run_cmd(cmd: List[str]) -> None:
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def count_vcf_records(path: Path, bcftools: str) -> int:
    proc = subprocess.run(
        [bcftools, "view", "-H", str(path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return sum(1 for line in proc.stdout.splitlines() if line.strip())


def compute_caller_counts(
    caller_vcf: Path,
    truth_vcf: Path,
    truth_bed: Optional[Path],
    bcftools: str,
    truth_total_hint: int,
) -> Dict[str, float]:
    with tempfile.TemporaryDirectory(prefix="compare_benchmark_") as tmpdir:
        tmp = Path(tmpdir)
        caller_scoped = tmp / "caller.in_scope.vcf.gz"
        truth_scoped = tmp / "truth.in_scope.vcf.gz"

        caller_cmd = [bcftools, "view", "-f", "PASS"]
        if truth_bed:
            caller_cmd.extend(["-R", str(truth_bed)])
        caller_cmd.extend([str(caller_vcf), "-Oz", "-o", str(caller_scoped)])
        run_cmd(caller_cmd)
        run_cmd([bcftools, "index", "-f", str(caller_scoped)])

        truth_cmd = [bcftools, "view"]
        if truth_bed:
            truth_cmd.extend(["-R", str(truth_bed)])
        truth_cmd.extend([str(truth_vcf), "-Oz", "-o", str(truth_scoped)])
        run_cmd(truth_cmd)
        run_cmd([bcftools, "index", "-f", str(truth_scoped)])

        isec_dir = tmp / "isec"
        run_cmd([bcftools, "isec", "-p", str(isec_dir), str(caller_scoped), str(truth_scoped)])

        tp = count_vcf_records(isec_dir / "0002.vcf", bcftools)
        fp = count_vcf_records(isec_dir / "0000.vcf", bcftools)
        truth_only = count_vcf_records(isec_dir / "0001.vcf", bcftools)

        truth_total = truth_total_hint if truth_total_hint > 0 else tp + truth_only
        metrics = compute_metrics(tp=tp, fp=fp, truth_total=truth_total)
        if truth_total_hint > 0 and truth_total_hint != tp + truth_only:
            metrics["fn"] = max(truth_total_hint - tp, 0)
            metrics["truth_total"] = truth_total_hint
        return metrics


def fmt(value: float) -> str:
    return f"{value:.4f}"


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()

    metrics_path = run_dir / "metrics.json"
    if not metrics_path.exists():
        raise SystemExit(f"metrics.json not found under {run_dir}")
    metrics = read_json(metrics_path)

    context_path = Path(args.context_json).resolve() if args.context_json else run_dir / "round_context.json"
    context = read_json(context_path) if context_path.exists() else {}

    sample = context.get("sample") or metrics.get("sample") or run_dir.name
    platform = context.get("platform") or infer_platform(sample)
    truth_set = context.get("truth_set") or infer_truth_set(sample)
    mode = context.get("analysis_mode") or metrics.get("mode", "")
    truth_total = int(metrics.get("truth_total", 0))
    baseline = metrics.get("baseline", {})
    filtered = metrics.get("filtered", {})
    rule = metrics.get("rule", {})
    counts_path = Path(args.baseline_counts_file).resolve() if args.baseline_counts_file else run_dir / "longphase_s" / "variant_counts.txt"
    counts = parse_variant_counts(counts_path)
    baseline_f1 = float(baseline.get("f1", 0.0))

    rows: List[Dict[str, str]] = []

    if not args.no_caller_input:
        caller_vcf = context.get("somatic_vcf")
        truth_vcf = context.get("truth_vcf")
        truth_bed = context.get("truth_bed")
        if caller_vcf and truth_vcf:
            caller_name = args.caller_name or context.get("caller_name") or infer_caller_name(caller_vcf)
            caller_metrics = compute_caller_counts(
                caller_vcf=Path(caller_vcf),
                truth_vcf=Path(truth_vcf),
                truth_bed=Path(truth_bed) if truth_bed else None,
                bcftools=args.bcftools,
                truth_total_hint=truth_total,
            )
            rows.append(
                {
                    "sample": sample,
                    "platform": platform,
                    "truth_set": truth_set,
                    "method": caller_name,
                    "caller": caller_name,
                    "mode": mode,
                    "truth_total": str(int(caller_metrics["truth_total"])),
                    "calls_total": str(int(caller_metrics["calls_total"])),
                    "tp": str(int(caller_metrics["tp"])),
                    "fp": str(int(caller_metrics["fp"])),
                    "fn": str(int(caller_metrics["fn"])),
                    "precision": fmt(caller_metrics["precision"]),
                    "recall": fmt(caller_metrics["recall"]),
                    "f1": fmt(caller_metrics["f1"]),
                    "delta_f1_vs_baseline": fmt(caller_metrics["f1"] - baseline_f1),
                    "notes": "caller-input PASS benchmark in truth scope",
                }
            )

    baseline_calls = int(counts.get("PASS_TOTAL", int(baseline.get("tp", 0)) + int(baseline.get("fp", 0))))
    rows.append(
        {
            "sample": sample,
            "platform": platform,
            "truth_set": truth_set,
            "method": args.baseline_method,
            "caller": context.get("caller_name") or infer_caller_name(context.get("somatic_vcf", "")),
            "mode": mode,
            "truth_total": str(truth_total),
            "calls_total": str(baseline_calls),
            "tp": str(int(baseline.get("tp", 0))),
            "fp": str(int(baseline.get("fp", 0))),
            "fn": str(max(truth_total - int(baseline.get("tp", 0)), 0)),
            "precision": fmt(float(baseline.get("precision", 0.0))),
            "recall": fmt(float(baseline.get("recall", 0.0))),
            "f1": fmt(float(baseline.get("f1", 0.0))),
            "delta_f1_vs_baseline": fmt(0.0),
            "notes": f"baseline from {counts_path}",
        }
    )

    filtered_tp = int(filtered.get("tp", 0))
    filtered_fp = int(filtered.get("fp", 0))
    rows.append(
        {
            "sample": sample,
            "platform": platform,
            "truth_set": truth_set,
            "method": "InterSubMod",
            "caller": context.get("caller_name") or infer_caller_name(context.get("somatic_vcf", "")),
            "mode": mode,
            "truth_total": str(truth_total),
            "calls_total": str(filtered_tp + filtered_fp),
            "tp": str(filtered_tp),
            "fp": str(filtered_fp),
            "fn": str(max(truth_total - filtered_tp, 0)),
            "precision": fmt(float(filtered.get("precision", 0.0))),
            "recall": fmt(float(filtered.get("recall", 0.0))),
            "f1": fmt(float(filtered.get("f1", 0.0))),
            "delta_f1_vs_baseline": fmt(float(filtered.get("f1", 0.0)) - baseline_f1),
            "notes": f"rule={rule.get('rule_id', 'unknown')}",
        }
    )

    write_tsv_rows(Path(args.output_tsv), FIELDS, rows)

    if args.output_md:
        md_columns = [
            "method",
            "caller",
            "truth_total",
            "calls_total",
            "tp",
            "fp",
            "fn",
            "precision",
            "recall",
            "f1",
            "delta_f1_vs_baseline",
        ]
        content = []
        content.append(f"# Benchmark Comparison: {sample}")
        content.append("")
        content.append(f"- Platform: `{platform}`")
        content.append(f"- Truth set: `{truth_set}`")
        content.append(f"- Source run: `{run_dir}`")
        content.append("")
        content.append(markdown_table(md_columns, rows))
        content.append("")
        content.append(f"- Output TSV: [{Path(args.output_tsv).name}]({Path(args.output_tsv).name})")
        Path(args.output_md).write_text("\n".join(content) + "\n", encoding="utf-8")

    print(f"[compare_benchmark_f1] Wrote {args.output_tsv}")
    if args.output_md:
        print(f"[compare_benchmark_f1] Wrote {args.output_md}")


if __name__ == "__main__":
    main()
