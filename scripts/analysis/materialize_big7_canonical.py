#!/usr/bin/env python3
"""Materialize big7 canonical manifests and summaries from migrated archives."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from big7_output_catalog import (  # noqa: E402
    ARCHIVE_ROOTS,
    CANONICAL_ROOT,
    MODE_SPECS,
    RUNBOOK_MANIFEST_ROOT,
    SAMPLE_SPECS,
    SYNTHESIS_ROOT,
    canonical_run_id,
)
from research_common import (  # noqa: E402
    ensure_dir,
    load_tsv_rows,
    read_json,
    write_json,
    write_tsv_rows,
)


RUN_MANIFEST_FIELDS = [
    "sample",
    "platform",
    "paired_or_to",
    "canonical_mode",
    "caller_family",
    "caller_model",
    "truth_set",
    "truth_total",
    "methylation_status",
    "source_kind",
    "archive_root_label",
    "archive_run_path",
    "canonical_run_path",
    "run_id",
    "tagged_bam_ready",
    "truth_ready",
    "baseline_ready",
    "intersubmod_ready",
    "metrics_ready",
    "decision_metrics_ready",
    "completeness_state",
    "blocking_reason",
]

EXPERIMENT_MATRIX_FIELDS = [
    "sample",
    "platform",
    "canonical_mode",
    "paired_or_to",
    "caller_family",
    "caller_model",
    "truth_set",
    "truth_total",
    "methylation_status",
    "expected",
    "discovered",
    "tagged_bam_ready",
    "truth_ready",
    "baseline_ready",
    "intersubmod_ready",
    "metrics_ready",
    "decision_metrics_ready",
    "completeness_state",
    "blocking_reason",
    "canonical_run_path",
]

SYNTHESIS_FIELDS = [
    "category",
    "archive_root_label",
    "name",
    "path",
    "sample_scope",
    "notes",
]


@dataclass
class DiscoveredRun:
    sample: str
    platform: str
    canonical_mode: str
    caller_family: str
    caller_model: str
    paired_or_to: str
    truth_set: str
    truth_total: int
    methylation_status: str
    source_kind: str
    archive_root_label: str
    archive_run_path: Path
    run_id: str
    source_name: str
    context: Dict
    sample_spec: Optional[object]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-root", default=str(CANONICAL_ROOT), help="Canonical output root")
    parser.add_argument("--synthesis-root", default=str(SYNTHESIS_ROOT), help="Synthesis output root")
    parser.add_argument("--runbook-manifest-root", default=str(RUNBOOK_MANIFEST_ROOT), help="Runbook manifest export root")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=max(os.cpu_count() or 1, 1),
        help="Maximum workers for archive scanning and benchmark normalization",
    )
    parser.add_argument("--skip-benchmark-generation", action="store_true", help="Do not generate benchmark_comparison.tsv for paired runs")
    return parser.parse_args()


def read_variant_counts(path: Path) -> Dict[str, str]:
    counts: Dict[str, str] = {}
    if not path.exists():
        return counts
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            counts[key.strip()] = value.strip()
    return counts


def path_exists(path_str: str) -> bool:
    if not path_str:
        return False
    return Path(path_str).exists()


def classify_file_role(path: Path) -> str:
    name = path.name
    if name == "metrics.json":
        return "metrics"
    if name in {"benchmark_comparison.tsv", "benchmark_comparison.md"}:
        return "benchmark_summary"
    if name == "round_context.json":
        return "round_context"
    if name == "significance_summary.csv":
        return "intersubmod_summary"
    if name == "haplotag_qc.tsv":
        return "haplotag_qc"
    if name.endswith(".bam") or name.endswith(".bam.bai"):
        return "tagged_bam" if "tagged" in name else "bam"
    if name.endswith(".vcf") or name.endswith(".vcf.gz") or name.endswith(".vcf.gz.csi") or name.endswith(".vcf.gz.tbi"):
        return "vcf"
    if name.endswith(".png"):
        return "plot"
    if name.endswith(".md"):
        return "report"
    if name.endswith(".log") or name == "run.log":
        return "log"
    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        return "table"
    return "other"


def base_context_for_sample(sample: str, canonical_mode: str) -> Dict[str, object]:
    spec = SAMPLE_SPECS[sample]
    mode_spec = MODE_SPECS[canonical_mode]
    context = {
        "sample": sample,
        "platform": spec.platform,
        "analysis_mode": canonical_mode,
        "truth_set": spec.truth_set,
        "truth_vcf": spec.truth_vcf,
        "truth_bed": spec.truth_bed,
        "truth_total": spec.truth_total,
        "caller_name": mode_spec.caller_family,
        "somatic_vcf": spec.somatic_vcf_pileup if mode_spec.caller_model == "pileup" else spec.somatic_vcf_full,
        "tumor_bam": spec.tumor_bam,
        "normal_bam": spec.normal_bam,
        "methylation_status": spec.methylation_status,
    }
    if canonical_mode.startswith("to_"):
        context["somatic_vcf"] = spec.to_somatic_vcf_pileup if mode_spec.caller_model == "pileup" else spec.to_somatic_vcf_full
    return context


def discover_official_paired_runs() -> List[DiscoveredRun]:
    runs: List[DiscoveredRun] = []
    for legacy_mode, canonical_mode in (("s-pure", "paired_full"), ("s-pure-pileup", "paired_pileup")):
        root = ARCHIVE_ROOTS["bip8_output_archive"] / legacy_mode
        if not root.exists():
            continue
        for sample_dir in sorted(root.iterdir()):
            if not sample_dir.is_dir() or sample_dir.name not in SAMPLE_SPECS:
                continue
            spec = SAMPLE_SPECS[sample_dir.name]
            for run_dir in sorted(sample_dir.iterdir()):
                if not run_dir.is_dir():
                    continue
                if run_dir.name.startswith("purity_"):
                    continue
                metrics_path = run_dir / "metrics.json"
                if not metrics_path.exists():
                    continue
                mode_spec = MODE_SPECS[canonical_mode]
                runs.append(
                    DiscoveredRun(
                        sample=spec.sample,
                        platform=spec.platform,
                        canonical_mode=canonical_mode,
                        caller_family=mode_spec.caller_family,
                        caller_model=mode_spec.caller_model,
                        paired_or_to=mode_spec.paired_or_to,
                        truth_set=spec.truth_set,
                        truth_total=spec.truth_total,
                        methylation_status=spec.methylation_status,
                        source_kind="official_pipeline_run",
                        archive_root_label="bip8_output_archive",
                        archive_run_path=run_dir,
                        run_id=canonical_run_id(spec.sample, canonical_mode, mode_spec.caller_model, run_dir.name),
                        source_name=run_dir.name,
                        context=base_context_for_sample(spec.sample, canonical_mode),
                        sample_spec=spec,
                    )
                )
    return runs


def discover_to_runs() -> List[DiscoveredRun]:
    runs: List[DiscoveredRun] = []
    research_root = ARCHIVE_ROOTS["bip8_output_archive"] / "research_rounds"
    if not research_root.exists():
        return runs
    for run_dir in sorted(research_root.iterdir()):
        if not run_dir.is_dir():
            continue
        context_path = run_dir / "round_context.json"
        benchmark_path = run_dir / "benchmark_comparison.tsv"
        if not context_path.exists() or not benchmark_path.exists():
            continue
        context = read_json(context_path)
        if str(context.get("analysis_mode", "")).lower() not in {"to-pure", "to_pure"}:
            continue
        sample = str(context.get("sample", "")).strip()
        if sample not in SAMPLE_SPECS:
            continue
        spec = SAMPLE_SPECS[sample]
        canonical_mode = "to_full" if "full" in run_dir.name.lower() else "to_pileup"
        mode_spec = MODE_SPECS[canonical_mode]
        merged_context = base_context_for_sample(sample, canonical_mode)
        merged_context.update(context)
        runs.append(
            DiscoveredRun(
                sample=sample,
                platform=str(context.get("platform") or spec.platform),
                canonical_mode=canonical_mode,
                caller_family=str(context.get("caller_name") or mode_spec.caller_family),
                caller_model=mode_spec.caller_model,
                paired_or_to=mode_spec.paired_or_to,
                truth_set=str(context.get("truth_set") or spec.truth_set),
                truth_total=int(context.get("truth_total") or spec.truth_total),
                methylation_status=spec.methylation_status,
                source_kind="tumor_only_pilot",
                archive_root_label="bip8_output_archive",
                archive_run_path=run_dir,
                run_id=canonical_run_id(sample, canonical_mode, mode_spec.caller_model, run_dir.name),
                source_name=run_dir.name,
                context=merged_context,
                sample_spec=spec,
            )
        )
    return runs


def discover_synthesis_items() -> List[Dict[str, str]]:
    items: List[Dict[str, str]] = []
    roots = [
        ("big8_output_archive", ARCHIVE_ROOTS["big8_output_archive"]),
        ("bip8_output_archive", ARCHIVE_ROOTS["bip8_output_archive"]),
    ]
    for label, root in roots:
        if not root.exists():
            continue
        for category in ["research_rounds", "pure_tumor_evaluation", "three_way_comparison", "multi_sample_quick_check"]:
            category_path = root / category
            if not category_path.exists():
                continue
            if category_path.is_file():
                items.append(
                    {
                        "category": category,
                        "archive_root_label": label,
                        "name": category_path.name,
                        "path": str(category_path),
                        "sample_scope": "mixed",
                        "notes": "top-level synthesis artifact",
                    }
                )
                continue
            if category == "three_way_comparison":
                items.append(
                    {
                        "category": category,
                        "archive_root_label": label,
                        "name": category_path.name,
                        "path": str(category_path),
                        "sample_scope": "mixed",
                        "notes": "three-way diagnostic root",
                    }
                )
                continue
            for child in sorted(category_path.iterdir()):
                if not child.exists():
                    continue
                items.append(
                    {
                        "category": category,
                        "archive_root_label": label,
                        "name": child.name,
                        "path": str(child),
                        "sample_scope": infer_sample_scope(child.name),
                        "notes": "archive synthesis artifact",
                    }
                )
    return items


def infer_sample_scope(name: str) -> str:
    for sample in SAMPLE_SPECS:
        if sample.lower() in name.lower():
            return sample
    return "mixed"


def write_file_manifest(source_root: Path, output_path: Path) -> None:
    rows = []
    for path in sorted(source_root.iterdir()):
        rel = path.relative_to(source_root)
        if path.is_file():
            rows.append(
                {
                    "relative_path": str(rel),
                    "source_path": str(path),
                    "entry_type": "file",
                    "size_bytes": str(path.stat().st_size),
                    "file_count": "1",
                    "role": classify_file_role(path),
                }
            )
            continue

        immediate_files = sum(1 for child in path.iterdir() if child.is_file())
        immediate_dirs = sum(1 for child in path.iterdir() if child.is_dir())
        rows.append(
            {
                "relative_path": str(rel),
                "source_path": str(path),
                "entry_type": "dir",
                "size_bytes": "",
                "file_count": str(immediate_files),
                "role": f"top_level_dir:{immediate_dirs}_subdirs",
            }
        )

    write_tsv_rows(output_path, ["relative_path", "source_path", "entry_type", "size_bytes", "file_count", "role"], rows)


def write_upstream_dependencies(run: DiscoveredRun, output_path: Path) -> None:
    rows = []
    context = run.context
    for label in ["tumor_bam", "normal_bam", "somatic_vcf", "truth_vcf", "truth_bed"]:
        raw = str(context.get(label, "") or "")
        rows.append(
            {
                "dependency": label,
                "path": raw,
                "exists": str(path_exists(raw)).lower(),
            }
        )
    write_tsv_rows(output_path, ["dependency", "path", "exists"], rows)


def copy_text_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    return True


def generate_paired_benchmark(run: DiscoveredRun, run_root: Path, skip_generation: bool) -> None:
    metrics_dir = ensure_dir(run_root / "metrics")
    if (metrics_dir / "benchmark_comparison.tsv").exists():
        return
    if skip_generation:
        return
    context_path = run_root / "manifest" / "run_context.json"
    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "compare_benchmark_f1.py"),
        "--run-dir",
        str(run.archive_run_path),
        "--context-json",
        str(context_path),
        "--output-tsv",
        str(metrics_dir / "benchmark_comparison.tsv"),
        "--output-md",
        str(metrics_dir / "benchmark_comparison.md"),
        "--baseline-method",
        MODE_SPECS[run.canonical_mode].baseline_method,
    ]
    try:
        proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if proc.stdout:
            print(proc.stdout, end="")
    except subprocess.CalledProcessError:
        fallback_cmd = cmd + ["--no-caller-input"]
        subprocess.run(fallback_cmd, check=True)


def run_completeness(run: DiscoveredRun, run_root: Path) -> Dict[str, str]:
    tagged_bam_candidates = [
        run.archive_run_path / "longphase_s" / f"{run.sample}_tagged.bam",
        run.archive_run_path / "step03_longphase_to" / "tumor_tagged.bam",
        run.archive_run_path / "step03_longphase_to" / "tumor_tagged.haplotag.bam",
    ]
    tagged_bam_path = next((path for path in tagged_bam_candidates if path.exists()), None)
    tp_summary = next(
        (
            path
            for path in [
                run.archive_run_path / "intersubmod_tp" / "significance_summary.csv",
                run.archive_run_path / "step05_intersubmod" / "intersubmod_tp" / "significance_summary.csv",
            ]
            if path.exists()
        ),
        None,
    )
    fp_summary = next(
        (
            path
            for path in [
                run.archive_run_path / "intersubmod_fp" / "significance_summary.csv",
                run.archive_run_path / "step05_intersubmod" / "intersubmod_fp" / "significance_summary.csv",
            ]
            if path.exists()
        ),
        None,
    )

    benchmark_path = run_root / "metrics" / "benchmark_comparison.tsv"
    metrics_json_path = run_root / "metrics" / "metrics.json"
    tagged_bam_ready = tagged_bam_path is not None
    truth_ready = path_exists(str(run.context.get("truth_vcf", ""))) and (
        not run.context.get("truth_bed") or path_exists(str(run.context.get("truth_bed", "")))
    )
    baseline_ready = benchmark_path.exists() or metrics_json_path.exists()
    intersubmod_ready = tp_summary is not None and fp_summary is not None
    metrics_ready = metrics_json_path.exists() and benchmark_path.exists()
    decision_metrics_ready = intersubmod_ready

    if tagged_bam_ready and truth_ready and baseline_ready and intersubmod_ready and metrics_ready and decision_metrics_ready:
        completeness_state = "completed"
        blocking_reason = ""
    elif truth_ready and (baseline_ready or metrics_ready):
        completeness_state = "partial"
        missing = []
        if not tagged_bam_ready:
            missing.append("missing_tagged_bam")
        if not intersubmod_ready:
            missing.append("missing_intersubmod_summary")
        if not metrics_ready:
            missing.append("missing_standard_benchmark_table")
        blocking_reason = ",".join(missing)
    else:
        completeness_state = "blocked"
        missing = []
        if not truth_ready:
            missing.append("missing_truth")
        if not baseline_ready:
            missing.append("missing_baseline")
        blocking_reason = ",".join(missing) or "missing_core_artifacts"

    row = {
        "tagged_bam_ready": str(tagged_bam_ready).lower(),
        "tagged_bam_path": str(tagged_bam_path or ""),
        "truth_ready": str(truth_ready).lower(),
        "baseline_ready": str(baseline_ready).lower(),
        "intersubmod_ready": str(intersubmod_ready).lower(),
        "metrics_ready": str(metrics_ready).lower(),
        "decision_metrics_ready": str(decision_metrics_ready).lower(),
        "tp_summary_path": str(tp_summary or ""),
        "fp_summary_path": str(fp_summary or ""),
        "completeness_state": completeness_state,
        "blocking_reason": blocking_reason,
    }
    write_tsv_rows(run_root / "metrics" / "completeness.tsv", list(row.keys()), [row])
    return row


def write_run_readme(run: DiscoveredRun, run_root: Path, completeness: Dict[str, str]) -> None:
    metrics_payload = {}
    metrics_path = run_root / "metrics" / "metrics.json"
    if metrics_path.exists():
        metrics_payload = read_json(metrics_path)
    lines = [
        f"# {run.sample} {run.canonical_mode}",
        "",
        f"- Source kind: `{run.source_kind}`",
        f"- Archive source: `{run.archive_run_path}`",
        f"- Canonical mode: `{run.canonical_mode}`",
        f"- Caller family/model: `{run.caller_family}` / `{run.caller_model}`",
        f"- Truth set: `{run.truth_set}` ({run.truth_total})",
        f"- Methylation: `{run.methylation_status}`",
        f"- Completeness: `{completeness['completeness_state']}`",
    ]
    if completeness["blocking_reason"]:
        lines.append(f"- Blocking reason: `{completeness['blocking_reason']}`")
    if metrics_payload:
        baseline = metrics_payload.get("baseline", {})
        filtered = metrics_payload.get("filtered", {})
        lines.extend(
            [
                "",
                "## Metrics",
                "",
                f"- Baseline TP/FP/F1: `{baseline.get('tp', '')}` / `{baseline.get('fp', '')}` / `{baseline.get('f1', '')}`",
                f"- Filtered TP/FP/F1: `{filtered.get('tp', '')}` / `{filtered.get('fp', '')}` / `{filtered.get('f1', '')}`",
            ]
        )
    lines.extend(
        [
            "",
            "## Canonical Contents",
            "",
            "- `manifest/`: source file inventory, upstream dependencies, migration provenance, normalized run context",
            "- `metrics/`: metrics.json, benchmark comparison, completeness table",
            "- `reports/`: this README and linked report references",
        ]
    )
    ensure_dir((run_root / "reports").parent)
    (run_root / "reports" / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def materialize_run(run: DiscoveredRun, canonical_root: Path, skip_benchmark_generation: bool) -> Dict[str, str]:
    run_root = canonical_root / run.sample / run.canonical_mode / run.run_id
    for subdir in ["inputs", "baseline", "intersubmod", "metrics", "artifacts", "reports", "manifest"]:
        ensure_dir(run_root / subdir)

    run_context = dict(run.context)
    run_context.update(
        {
            "canonical_mode": run.canonical_mode,
            "source_kind": run.source_kind,
            "archive_root_label": run.archive_root_label,
            "archive_run_path": str(run.archive_run_path),
            "canonical_run_path": str(run_root),
            "run_id": run.run_id,
            "materialized_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
    )
    write_json(run_root / "manifest" / "run_context.json", run_context)
    write_file_manifest(run.archive_run_path, run_root / "manifest" / "file_manifest.tsv")
    write_upstream_dependencies(run, run_root / "manifest" / "upstream_dependency.tsv")
    write_tsv_rows(
        run_root / "manifest" / "migration_provenance.tsv",
        ["archive_root_label", "archive_run_path", "canonical_run_path", "source_kind", "run_id"],
        [
            {
                "archive_root_label": run.archive_root_label,
                "archive_run_path": str(run.archive_run_path),
                "canonical_run_path": str(run_root),
                "source_kind": run.source_kind,
                "run_id": run.run_id,
            }
        ],
    )

    copy_text_if_exists(run.archive_run_path / "metrics.json", run_root / "metrics" / "metrics.json")
    if not copy_text_if_exists(run.archive_run_path / "benchmark_comparison.tsv", run_root / "metrics" / "benchmark_comparison.tsv"):
        if run.canonical_mode.startswith("paired_"):
            generate_paired_benchmark(run, run_root, skip_benchmark_generation)
    copy_text_if_exists(run.archive_run_path / "benchmark_comparison.md", run_root / "metrics" / "benchmark_comparison.md")

    completeness = run_completeness(run, run_root)
    write_run_readme(run, run_root, completeness)

    return {
        "sample": run.sample,
        "platform": run.platform,
        "paired_or_to": run.paired_or_to,
        "canonical_mode": run.canonical_mode,
        "caller_family": run.caller_family,
        "caller_model": run.caller_model,
        "truth_set": run.truth_set,
        "truth_total": str(run.truth_total),
        "methylation_status": run.methylation_status,
        "source_kind": run.source_kind,
        "archive_root_label": run.archive_root_label,
        "archive_run_path": str(run.archive_run_path),
        "canonical_run_path": str(run_root),
        "run_id": run.run_id,
        "tagged_bam_ready": completeness["tagged_bam_ready"],
        "truth_ready": completeness["truth_ready"],
        "baseline_ready": completeness["baseline_ready"],
        "intersubmod_ready": completeness["intersubmod_ready"],
        "metrics_ready": completeness["metrics_ready"],
        "decision_metrics_ready": completeness["decision_metrics_ready"],
        "completeness_state": completeness["completeness_state"],
        "blocking_reason": completeness["blocking_reason"],
    }


def build_experiment_matrix(run_rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    best_by_key: Dict[tuple[str, str], Dict[str, str]] = {}
    for row in run_rows:
        key = (row["sample"], row["canonical_mode"])
        prev = best_by_key.get(key)
        if prev is None or row["run_id"] > prev["run_id"]:
            best_by_key[key] = row

    rows: List[Dict[str, str]] = []
    for spec in SAMPLE_SPECS.values():
        for canonical_mode, mode_spec in MODE_SPECS.items():
            key = (spec.sample, canonical_mode)
            discovered = best_by_key.get(key)
            expected = canonical_mode in spec.expected_modes
            if discovered:
                rows.append(
                    {
                        "sample": spec.sample,
                        "platform": spec.platform,
                        "canonical_mode": canonical_mode,
                        "paired_or_to": mode_spec.paired_or_to,
                        "caller_family": mode_spec.caller_family,
                        "caller_model": mode_spec.caller_model,
                        "truth_set": spec.truth_set,
                        "truth_total": str(spec.truth_total),
                        "methylation_status": spec.methylation_status,
                        "expected": str(expected).lower(),
                        "discovered": "true",
                        "tagged_bam_ready": discovered["tagged_bam_ready"],
                        "truth_ready": discovered["truth_ready"],
                        "baseline_ready": discovered["baseline_ready"],
                        "intersubmod_ready": discovered["intersubmod_ready"],
                        "metrics_ready": discovered["metrics_ready"],
                        "decision_metrics_ready": discovered["decision_metrics_ready"],
                        "completeness_state": discovered["completeness_state"],
                        "blocking_reason": discovered["blocking_reason"],
                        "canonical_run_path": discovered["canonical_run_path"],
                    }
                )
                continue

            if not expected:
                blocking = spec.mode_unavailable_reason(canonical_mode)
                rows.append(
                    {
                        "sample": spec.sample,
                        "platform": spec.platform,
                        "canonical_mode": canonical_mode,
                        "paired_or_to": mode_spec.paired_or_to,
                        "caller_family": mode_spec.caller_family,
                        "caller_model": mode_spec.caller_model,
                        "truth_set": spec.truth_set,
                        "truth_total": str(spec.truth_total),
                        "methylation_status": spec.methylation_status,
                        "expected": "false",
                        "discovered": "false",
                        "tagged_bam_ready": "false",
                        "truth_ready": str(path_exists(spec.truth_vcf) and (not spec.truth_bed or path_exists(spec.truth_bed))).lower(),
                        "baseline_ready": "false",
                        "intersubmod_ready": "false",
                        "metrics_ready": "false",
                        "decision_metrics_ready": "false",
                        "completeness_state": "not_applicable",
                        "blocking_reason": blocking,
                        "canonical_run_path": "",
                    }
                )
                continue

            if mode_spec.paired_or_to == "paired":
                blocking = "no_discovered_paired_run"
            else:
                blocking = "tumor_only_run_missing_or_partial"
            rows.append(
                {
                    "sample": spec.sample,
                    "platform": spec.platform,
                    "canonical_mode": canonical_mode,
                    "paired_or_to": mode_spec.paired_or_to,
                    "caller_family": mode_spec.caller_family,
                    "caller_model": mode_spec.caller_model,
                    "truth_set": spec.truth_set,
                    "truth_total": str(spec.truth_total),
                    "methylation_status": spec.methylation_status,
                    "expected": "true",
                    "discovered": "false",
                    "tagged_bam_ready": "false",
                    "truth_ready": str(path_exists(spec.truth_vcf) and (not spec.truth_bed or path_exists(spec.truth_bed))).lower(),
                    "baseline_ready": "false",
                    "intersubmod_ready": "false",
                    "metrics_ready": "false",
                    "decision_metrics_ready": "false",
                    "completeness_state": "blocked",
                    "blocking_reason": blocking,
                    "canonical_run_path": "",
                }
            )
    return rows


def write_master_report(synthesis_root: Path, run_rows: List[Dict[str, str]], experiment_rows: List[Dict[str, str]], synthesis_rows: List[Dict[str, str]]) -> None:
    counts = Counter((row["canonical_mode"], row["completeness_state"]) for row in run_rows)
    runnable_gaps = [row for row in experiment_rows if row["expected"] == "true" and row["completeness_state"] != "completed"]
    unavailable_rows = [row for row in experiment_rows if row["expected"] != "true"]
    lines = [
        "# Big7 Canonical Output Report",
        "",
        f"- Generated at: `{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}`",
        f"- Canonical runs discovered: `{len(run_rows)}`",
        f"- Synthesis/archive items indexed: `{len(synthesis_rows)}`",
        f"- CPU policy: `nproc / full-machine CPU allowed; correctness and space gates remain enforced`",
        "",
        "## Canonical Run Counts",
        "",
        "| Canonical Mode | Completed | Partial | Blocked |",
        "| --- | --- | --- | --- |",
    ]
    for canonical_mode in MODE_SPECS:
        lines.append(
            f"| {canonical_mode} | {counts.get((canonical_mode, 'completed'), 0)} | "
            f"{counts.get((canonical_mode, 'partial'), 0)} | {counts.get((canonical_mode, 'blocked'), 0)} |"
        )
    lines.extend(
        [
            "",
            "## Remaining Runnable Gaps",
            "",
            "| Sample | Mode | State | Blocking Reason |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in runnable_gaps:
        lines.append(
            f"| {row['sample']} | {row['canonical_mode']} | {row['completeness_state']} | {row['blocking_reason']} |"
        )
    if unavailable_rows:
        lines.extend(
            [
                "",
                "## Unavailable Modes",
                "",
                "| Sample | Mode | Reason |",
                "| --- | --- | --- |",
            ]
        )
        for row in unavailable_rows:
            lines.append(f"| {row['sample']} | {row['canonical_mode']} | {row['blocking_reason']} |")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Canonical bundles only materialize manifests and small summary files; large BAM/VCF payloads remain at archive source paths.",
            "- Symlinked source data are not recopied as physical outputs.",
            "- `expected=false` means the current source tree does not provide that caller/model path, so it is not counted as a runnable gap.",
            "- Research rounds, pure_tumor_evaluation, and three_way_comparison remain indexed under synthesis/archive manifests instead of being flattened into sample-mode baselines.",
        ]
    )
    (synthesis_root / "master_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_runbook_copies(runbook_root: Path, run_rows: List[Dict[str, str]], experiment_rows: List[Dict[str, str]], synthesis_rows: List[Dict[str, str]]) -> None:
    ensure_dir(runbook_root)
    write_tsv_rows(runbook_root / "20260314_big7_master_run_manifest.tsv", RUN_MANIFEST_FIELDS, run_rows)
    write_tsv_rows(runbook_root / "20260314_big7_experiment_matrix.tsv", EXPERIMENT_MATRIX_FIELDS, experiment_rows)
    write_tsv_rows(runbook_root / "20260314_big7_synthesis_manifest.tsv", SYNTHESIS_FIELDS, synthesis_rows)


def main() -> None:
    args = parse_args()
    canonical_root = ensure_dir(Path(args.canonical_root))
    synthesis_root = ensure_dir(Path(args.synthesis_root))
    runbook_manifest_root = ensure_dir(Path(args.runbook_manifest_root))

    discovered_runs = discover_official_paired_runs() + discover_to_runs()
    with ThreadPoolExecutor(max_workers=max(args.max_workers, 1)) as pool:
        run_rows = list(
            pool.map(
                lambda run: materialize_run(run, canonical_root, args.skip_benchmark_generation),
                discovered_runs,
            )
        )
    synthesis_rows = discover_synthesis_items()
    experiment_rows = build_experiment_matrix(run_rows)

    write_tsv_rows(synthesis_root / "master_run_manifest.tsv", RUN_MANIFEST_FIELDS, run_rows)
    write_tsv_rows(synthesis_root / "master_experiment_matrix.tsv", EXPERIMENT_MATRIX_FIELDS, experiment_rows)
    write_tsv_rows(synthesis_root / "archive_synthesis_manifest.tsv", SYNTHESIS_FIELDS, synthesis_rows)
    write_master_report(synthesis_root, run_rows, experiment_rows, synthesis_rows)
    export_runbook_copies(runbook_manifest_root, run_rows, experiment_rows, synthesis_rows)


if __name__ == "__main__":
    main()
