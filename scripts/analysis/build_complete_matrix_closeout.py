#!/usr/bin/env python3
"""Build a final closeout report and checklist for the big7 complete experiment matrix."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


PAIR_SAMPLES = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

TO_PILEUP_SAMPLES = list(PAIR_SAMPLES)

SAMPLE_TRUTH = {
    "HCC1395": (
        "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
        "/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed",
    ),
    "HCC1395_DORADO": (
        "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
        "/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed",
    ),
    "COLO829": (
        "/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz",
        "",
    ),
    "H1437": (
        "/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        "/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark.bed",
    ),
    "H2009": (
        "/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        "/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark.bed",
    ),
    "HCC1937": (
        "/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        "/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark.bed",
    ),
    "HCC1954": (
        "/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        "/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark.bed",
    ),
}

PAIR_NORMAL_BAM = {
    "HCC1395": "/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam",
    "HCC1395_DORADO": "/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam",
    "COLO829": "/big8_disk/data/COLO829/ONT_PAO/PAO33946.bam",
    "H1437": "/big8_disk/data/H1437/ONT/H1437BL.bam",
    "H2009": "/big8_disk/data/H2009/ONT/H2009BL.bam",
    "HCC1937": "/big8_disk/data/HCC1937/ONT/HCC1937BL.bam",
    "HCC1954": "/big8_disk/data/HCC1954/ONT/HCC1954BL.bam",
}


@dataclass(frozen=True)
class RunSpec:
    sample: str
    mode: str
    expected: bool
    unavailable_reason: str = ""


def build_specs() -> List[RunSpec]:
    specs: List[RunSpec] = []
    for sample in PAIR_SAMPLES:
        specs.append(RunSpec(sample=sample, mode="paired_full", expected=True))
        specs.append(RunSpec(sample=sample, mode="paired_pileup", expected=True))
        if sample in TO_PILEUP_SAMPLES:
            specs.append(RunSpec(sample=sample, mode="to_pileup", expected=True))
            specs.append(
                RunSpec(
                    sample=sample,
                    mode="to_full",
                    expected=False,
                    unavailable_reason="tumor_only_full_model_unavailable",
                )
            )
    return specs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canonical-root",
        default="/big7_disk/liaoyoyo2001/big7_disk_output/canonical",
        help="Canonical paired output root",
    )
    parser.add_argument(
        "--synthesis-root",
        default="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis",
        help="Synthesis output root",
    )
    parser.add_argument(
        "--output-dir",
        default="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/final_closeout",
        help="Closeout output directory",
    )
    parser.add_argument(
        "--batch-root",
        action="append",
        default=[],
        help="Optional batch root(s) to summarize",
    )
    return parser.parse_args()


def read_json(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def read_tsv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def canonical_candidate_score(path: Path, sample: str) -> Tuple[int, str]:
    score = 0
    if (path / "metrics.json").exists():
        score += 5
    if (path / "benchmark_comparison.tsv").exists() or (path / "metrics" / "benchmark_comparison.tsv").exists():
        score += 4
    if (path / "intersubmod_tp" / "significance_summary.csv").exists():
        score += 3
    if (path / "intersubmod_fp" / "significance_summary.csv").exists():
        score += 3
    if (path / "longphase_s" / "filtered_snv_tp.vcf.gz").exists():
        score += 2
    if (path / "longphase_s" / "filtered_snv_fp.vcf.gz").exists():
        score += 2
    if (path / "longphase_s" / f"{sample}_tagged.bam").exists():
        score += 1
    if (path / "longphase_s" / f"{sample}_tagged.bam.bai").exists():
        score += 1
    return score, path.name


def latest_canonical_run(canonical_root: Path, sample: str, mode: str) -> Optional[Path]:
    mode_root = canonical_root / sample / mode
    if not mode_root.exists():
        return None
    candidates = [path for path in mode_root.iterdir() if path.is_dir()]
    if not candidates:
        return None
    ranked = sorted(candidates, key=lambda path: canonical_candidate_score(path, sample))
    return ranked[-1]


def to_candidates(synthesis_root: Path, sample: str) -> List[Path]:
    round_root = synthesis_root / "research_rounds"
    if not round_root.exists():
        return []
    token = sample.lower()
    return sorted(
        [
            path
            for path in round_root.iterdir()
            if path.is_dir() and token in path.name.lower() and "to_pilot" in path.name.lower()
        ]
    )


def candidate_score(path: Path) -> Tuple[int, str]:
    score = 0
    if (path / "metrics.json").exists():
        score += 4
    if (path / "benchmark_comparison.tsv").exists():
        score += 3
    if (path / "step05_intersubmod" / "intersubmod_tp" / "significance_summary.csv").exists():
        score += 2
    if (path / "step03_longphase_to" / "tumor_tagged.bam").exists():
        score += 1
    return score, path.name


def latest_to_run(synthesis_root: Path, sample: str) -> Optional[Path]:
    candidates = to_candidates(synthesis_root, sample)
    if not candidates:
        return None
    ranked = sorted(candidates, key=candidate_score)
    return ranked[-1]


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def parse_benchmark_table(path: Optional[Path]) -> Dict[str, Dict[str, float]]:
    if path is None or not path.exists():
        return {}
    rows = read_tsv_rows(path)
    metrics: Dict[str, Dict[str, float]] = {}
    for row in rows:
        method = row.get("method", "").lower()
        key = "caller"
        if "intersubmod" in method:
            key = "intersubmod"
        elif "longphase" in method:
            key = "longphase"
        metrics[key] = {
            "tp": float(row.get("tp", 0) or 0),
            "fp": float(row.get("fp", 0) or 0),
            "fn": float(row.get("fn", 0) or 0),
            "f1": float(row.get("f1", 0) or 0),
        }
    return metrics


def parse_stage_metrics_json(path: Optional[Path]) -> Dict[str, Dict[str, float]]:
    if path is None or not path.exists():
        return {}
    payload = read_json(path)
    metrics: Dict[str, Dict[str, float]] = {}
    mapping = {
        "baseline_clairs_to": "caller",
        "longphase_to": "longphase",
        "intersubmod": "intersubmod",
    }
    for source_key, target_key in mapping.items():
        if source_key not in payload:
            continue
        source = payload[source_key]
        if not isinstance(source, dict):
            continue
        metrics[target_key] = {
            "tp": float(source.get("tp", 0) or 0),
            "fp": float(source.get("fp", 0) or 0),
            "fn": float(source.get("fn", 0) or 0),
            "f1": float(source.get("f1", 0) or 0),
        }
    return metrics


def parse_metrics_json(path: Optional[Path]) -> Dict[str, Dict[str, float]]:
    if path is None or not path.exists():
        return {}
    payload = read_json(path)
    metrics: Dict[str, Dict[str, float]] = {}
    for key in ("baseline", "filtered"):
        source = payload.get(key)
        if not isinstance(source, dict):
            continue
        target_key = "longphase" if key == "baseline" else "intersubmod"
        tp = float(source.get("tp", 0) or 0)
        fp = float(source.get("fp", 0) or 0)
        truth_total = float(payload.get("truth_total", 0) or 0)
        metrics[target_key] = {
            "tp": tp,
            "fp": fp,
            "fn": max(truth_total - tp, 0),
            "f1": float(source.get("f1", 0) or 0),
        }
    return metrics


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def fmt_metric(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:.4f}"


def collect_batch_status(batch_root: Path) -> Dict[str, object]:
    status_path = batch_root / "job_status.tsv"
    if not status_path.exists():
        return {
            "batch_root": str(batch_root),
            "queued": 0,
            "running": 0,
            "completed": 0,
            "failed": 0,
            "status_path": "",
        }
    rows = read_tsv_rows(status_path)
    counts = Counter(row.get("status", "unknown") for row in rows)
    return {
        "batch_root": str(batch_root),
        "queued": counts.get("queued", 0),
        "running": counts.get("running", 0),
        "completed": counts.get("completed", 0),
        "failed": counts.get("failed", 0),
        "status_path": str(status_path),
    }


def gather_run_row(spec: RunSpec, canonical_root: Path, synthesis_root: Path) -> Dict[str, object]:
    row: Dict[str, object] = {
        "sample": spec.sample,
        "mode": spec.mode,
        "expected": bool_text(spec.expected),
        "availability_reason": spec.unavailable_reason,
        "run_dir": "",
        "source_kind": "",
        "status": "",
        "blocking_reason": "",
    }

    if not spec.expected:
        row.update(
            {
                "status": "not_applicable",
                "blocking_reason": spec.unavailable_reason,
            }
        )
        return row

    if spec.mode.startswith("paired_"):
        run_dir = latest_canonical_run(canonical_root, spec.sample, spec.mode)
        source_kind = "canonical"
        tagged_bam = run_dir / "longphase_s" / f"{spec.sample}_tagged.bam" if run_dir else None
        tagged_bai = run_dir / "longphase_s" / f"{spec.sample}_tagged.bam.bai" if run_dir else None
        tp_vcf = run_dir / "longphase_s" / "filtered_snv_tp.vcf.gz" if run_dir else None
        fp_vcf = run_dir / "longphase_s" / "filtered_snv_fp.vcf.gz" if run_dir else None
        tp_summary = run_dir / "intersubmod_tp" / "significance_summary.csv" if run_dir else None
        fp_summary = run_dir / "intersubmod_fp" / "significance_summary.csv" if run_dir else None
        metrics_path = run_dir / "metrics.json" if run_dir else None
        benchmark_path = first_existing(
            [
                run_dir / "benchmark_comparison.tsv",
                run_dir / "metrics" / "benchmark_comparison.tsv",
            ]
        ) if run_dir else None
        dashboard_path = first_existing([run_dir / "round_summary.md", run_dir / "reports" / "README.md"]) if run_dir else None
        method_design_path = run_dir / "method_design_validation.tsv" if run_dir else None
        context_path = run_dir / "run_context.json" if run_dir else None
        haplotag_qc = run_dir / "longphase_s" / "haplotag_qc.tsv" if run_dir else None
        to_stage_metrics = None
        truth_vcf_raw, truth_bed_raw = SAMPLE_TRUTH[spec.sample]
    else:
        run_dir = latest_to_run(synthesis_root, spec.sample)
        source_kind = "research_round"
        tagged_bam = run_dir / "step03_longphase_to" / "tumor_tagged.bam" if run_dir else None
        tagged_bai = run_dir / "step03_longphase_to" / "tumor_tagged.bam.bai" if run_dir else None
        tp_vcf = run_dir / "step04_benchmark_longphase_to" / "filtered_snv_tp.vcf.gz" if run_dir else None
        fp_vcf = run_dir / "step04_benchmark_longphase_to" / "filtered_snv_fp.vcf.gz" if run_dir else None
        tp_summary = run_dir / "step05_intersubmod" / "intersubmod_tp" / "significance_summary.csv" if run_dir else None
        fp_summary = run_dir / "step05_intersubmod" / "intersubmod_fp" / "significance_summary.csv" if run_dir else None
        metrics_path = run_dir / "metrics.json" if run_dir else None
        benchmark_path = run_dir / "benchmark_comparison.tsv" if run_dir else None
        dashboard_path = run_dir / "round_summary.md" if run_dir else None
        method_design_path = run_dir / "method_design_validation.tsv" if run_dir else None
        context_path = run_dir / "round_context.json" if run_dir else None
        haplotag_qc = run_dir / "haplotag_qc.tsv" if run_dir else None
        to_stage_metrics = run_dir / "to_stage_metrics.json" if run_dir else None

        truth_vcf_raw = ""
        truth_bed_raw = ""
        if context_path and context_path.exists():
            context_payload = read_json(context_path)
            truth_vcf_raw = str(context_payload.get("truth_vcf", "") or "")
            truth_bed_raw = str(context_payload.get("truth_bed", "") or "")

    if run_dir is None:
        row.update(
            {
                "status": "missing",
                "blocking_reason": "run_dir_missing",
                "source_kind": source_kind,
            }
        )
        return row

    truth_vcf = Path(truth_vcf_raw) if truth_vcf_raw else None
    truth_bed = Path(truth_bed_raw) if truth_bed_raw else None

    benchmark_metrics = parse_benchmark_table(benchmark_path)
    stage_metrics = parse_stage_metrics_json(to_stage_metrics)
    metrics_json_metrics = parse_metrics_json(metrics_path)

    merged_metrics: Dict[str, Dict[str, float]] = {}
    for source in (metrics_json_metrics, stage_metrics, benchmark_metrics):
        merged_metrics.update(source)

    tagged_bam_ready = bool(tagged_bam and tagged_bam.exists())
    tagged_bai_ready = bool(tagged_bai and tagged_bai.exists())
    truth_ready = bool(truth_vcf and truth_vcf.exists() and (truth_bed is None or str(truth_bed) == "" or truth_bed.exists()))
    longphase_ready = bool(tp_vcf and tp_vcf.exists() and fp_vcf and fp_vcf.exists())
    intersubmod_ready = bool(tp_summary and tp_summary.exists() and fp_summary and fp_summary.exists())
    metrics_ready = bool(metrics_path and metrics_path.exists())
    benchmark_ready = bool(benchmark_path and benchmark_path.exists()) or bool(merged_metrics)
    dashboard_ready = bool(dashboard_path and dashboard_path.exists())
    method_design_ready = bool(method_design_path and method_design_path.exists())
    haplotag_qc_ready = bool(haplotag_qc and haplotag_qc.exists())

    missing_checks = []
    if not tagged_bam_ready:
        missing_checks.append("tagged_bam")
    if not tagged_bai_ready:
        missing_checks.append("tagged_bam_index")
    if not truth_ready:
        missing_checks.append("truth")
    if not longphase_ready:
        missing_checks.append("longphase_tp_fp")
    if not intersubmod_ready:
        missing_checks.append("intersubmod_summary")
    if not metrics_ready:
        missing_checks.append("metrics_json")
    if not benchmark_ready:
        missing_checks.append("stage_metrics")

    status = "completed" if not missing_checks else "partial"
    blocking_reason = ",".join(missing_checks)

    caller = merged_metrics.get("caller", {})
    longphase = merged_metrics.get("longphase", {})
    intersubmod = merged_metrics.get("intersubmod", {})
    delta_vs_longphase = None
    if longphase.get("f1") is not None and intersubmod.get("f1") is not None:
        delta_vs_longphase = intersubmod.get("f1", 0.0) - longphase.get("f1", 0.0)

    row.update(
        {
            "run_dir": str(run_dir),
            "source_kind": source_kind,
            "status": status,
            "blocking_reason": blocking_reason,
            "tagged_bam_ready": bool_text(tagged_bam_ready),
            "tagged_bai_ready": bool_text(tagged_bai_ready),
            "truth_ready": bool_text(truth_ready),
            "longphase_ready": bool_text(longphase_ready),
            "intersubmod_ready": bool_text(intersubmod_ready),
            "metrics_ready": bool_text(metrics_ready),
            "benchmark_ready": bool_text(benchmark_ready),
            "dashboard_ready": bool_text(dashboard_ready),
            "method_design_ready": bool_text(method_design_ready),
            "haplotag_qc_ready": bool_text(haplotag_qc_ready),
            "tagged_bam_path": str(tagged_bam) if tagged_bam else "",
            "tagged_bai_path": str(tagged_bai) if tagged_bai else "",
            "truth_vcf_path": str(truth_vcf) if truth_vcf else "",
            "truth_bed_path": str(truth_bed) if truth_bed else "",
            "tp_vcf_path": str(tp_vcf) if tp_vcf else "",
            "fp_vcf_path": str(fp_vcf) if fp_vcf else "",
            "tp_summary_path": str(tp_summary) if tp_summary else "",
            "fp_summary_path": str(fp_summary) if fp_summary else "",
            "metrics_path": str(metrics_path) if metrics_path else "",
            "benchmark_path": str(benchmark_path) if benchmark_path else "",
            "dashboard_path": str(dashboard_path) if dashboard_path else "",
            "method_design_path": str(method_design_path) if method_design_path else "",
            "haplotag_qc_path": str(haplotag_qc) if haplotag_qc else "",
            "context_path": str(context_path) if context_path else "",
            "caller_tp": int(caller.get("tp", 0)),
            "caller_fp": int(caller.get("fp", 0)),
            "caller_fn": int(caller.get("fn", 0)),
            "caller_f1": fmt_metric(caller.get("f1")),
            "longphase_tp": int(longphase.get("tp", 0)),
            "longphase_fp": int(longphase.get("fp", 0)),
            "longphase_fn": int(longphase.get("fn", 0)),
            "longphase_f1": fmt_metric(longphase.get("f1")),
            "intersubmod_tp": int(intersubmod.get("tp", 0)),
            "intersubmod_fp": int(intersubmod.get("fp", 0)),
            "intersubmod_fn": int(intersubmod.get("fn", 0)),
            "intersubmod_f1": fmt_metric(intersubmod.get("f1")),
            "delta_f1_vs_longphase": fmt_metric(delta_vs_longphase),
        }
    )
    return row


def markdown_table(headers: List[str], rows: List[List[str]]) -> List[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return lines


def checkbox(value: bool) -> str:
    return "[x]" if value else "[ ]"


def write_markdown(path: Path, rows: List[Dict[str, object]], batch_rows: List[Dict[str, object]]) -> None:
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    counts = Counter(row["status"] for row in rows)
    expected_rows = [row for row in rows if row["expected"] == "true"]
    completed_expected = [row for row in expected_rows if row["status"] == "completed"]
    pending_expected = [row for row in expected_rows if row["status"] != "completed"]
    unavailable_rows = [row for row in rows if row["status"] == "not_applicable"]

    lines: List[str] = [
        "# Big7 Complete Matrix Closeout",
        "",
        f"- Generated at: `{generated_at}`",
        f"- Expected runnable runs: `{len(expected_rows)}`",
        f"- Completed expected runs: `{len(completed_expected)}`",
        f"- Remaining expected runs: `{len(pending_expected)}`",
        f"- Not applicable modes: `{len(unavailable_rows)}`",
        "",
        "## Global Checklist",
        "",
        f"- {checkbox(len(pending_expected) == 0)} All expected runs reached `completed`.",
        f"- {checkbox(all(row.get('tagged_bam_ready') == 'true' for row in expected_rows))} All expected runs have tagged BAM.",
        f"- {checkbox(all(row.get('tagged_bai_ready') == 'true' for row in expected_rows))} All expected runs have BAM index.",
        f"- {checkbox(all(row.get('truth_ready') == 'true' for row in expected_rows))} All expected runs have truth data in place.",
        f"- {checkbox(all(row.get('longphase_ready') == 'true' for row in expected_rows))} All expected runs have LongPhase TP/FP VCFs.",
        f"- {checkbox(all(row.get('intersubmod_ready') == 'true' for row in expected_rows))} All expected runs have InterSubMod TP/FP summaries.",
        f"- {checkbox(all(row.get('metrics_ready') == 'true' for row in expected_rows))} All expected runs have metrics.json.",
        f"- {checkbox(all(row.get('benchmark_ready') == 'true' for row in expected_rows))} All expected runs have stage-level TP/FP/FN/F1 summaries.",
        "",
    ]

    if batch_rows:
        lines.extend(["## Batch Status", ""])
        batch_table = [
            [
                Path(str(row["batch_root"])).name,
                str(row["queued"]),
                str(row["running"]),
                str(row["completed"]),
                str(row["failed"]),
                str(row["status_path"]),
            ]
            for row in batch_rows
        ]
        lines.extend(markdown_table(["Batch", "Queued", "Running", "Completed", "Failed", "Status Path"], batch_table))
        lines.append("")

    if pending_expected:
        lines.extend(["## Remaining Expected Runs", ""])
        pending_table = [
            [
                str(row["sample"]),
                str(row["mode"]),
                str(row["status"]),
                str(row["blocking_reason"]),
                str(row["run_dir"]),
            ]
            for row in pending_expected
        ]
        lines.extend(markdown_table(["Sample", "Mode", "Status", "Blocking", "Run Dir"], pending_table))
        lines.append("")

    if unavailable_rows:
        lines.extend(["## Not Applicable Modes", ""])
        unavailable_table = [
            [
                str(row["sample"]),
                str(row["mode"]),
                str(row["availability_reason"]),
            ]
            for row in unavailable_rows
        ]
        lines.extend(markdown_table(["Sample", "Mode", "Reason"], unavailable_table))
        lines.append("")

    lines.extend(["## Per-Run Detailed Checklist", ""])
    for row in rows:
        lines.extend(
            [
                f"### {row['sample']} {row['mode']}",
                "",
                f"- Status: `{row['status']}`",
                f"- Run dir: `{row['run_dir']}`",
                f"- Blocking reason: `{row['blocking_reason']}`" if row["blocking_reason"] else "- Blocking reason: `none`",
                f"- {checkbox(row.get('tagged_bam_ready') == 'true')} tagged BAM",
                f"- {checkbox(row.get('tagged_bai_ready') == 'true')} tagged BAM index",
                f"- {checkbox(row.get('truth_ready') == 'true')} truth VCF / BED",
                f"- {checkbox(row.get('longphase_ready') == 'true')} LongPhase TP / FP VCF",
                f"- {checkbox(row.get('haplotag_qc_ready') == 'true')} haplotag QC",
                f"- {checkbox(row.get('intersubmod_ready') == 'true')} InterSubMod TP / FP summaries",
                f"- {checkbox(row.get('metrics_ready') == 'true')} metrics.json",
                f"- {checkbox(row.get('benchmark_ready') == 'true')} stage metrics table",
                f"- {checkbox(row.get('dashboard_ready') == 'true')} dashboard / summary doc",
                f"- {checkbox(row.get('method_design_ready') == 'true')} method design validation",
                f"- Caller F1: `{row.get('caller_f1', '')}`",
                f"- LongPhase F1: `{row.get('longphase_f1', '')}`",
                f"- InterSubMod F1: `{row.get('intersubmod_f1', '')}`",
                f"- Delta F1 vs LongPhase: `{row.get('delta_f1_vs_longphase', '')}`",
                "",
            ]
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    canonical_root = Path(args.canonical_root).resolve()
    synthesis_root = Path(args.synthesis_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    specs = build_specs()
    rows = [gather_run_row(spec, canonical_root, synthesis_root) for spec in specs]

    batch_roots = [Path(path).resolve() for path in args.batch_root]
    batch_rows = [collect_batch_status(path) for path in batch_roots]

    write_tsv(output_dir / "completion_checklist.tsv", rows)
    write_markdown(output_dir / "completion_report.md", rows, batch_rows)

    summary = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "expected_total": sum(1 for row in rows if row["expected"] == "true"),
        "expected_completed": sum(1 for row in rows if row["expected"] == "true" and row["status"] == "completed"),
        "expected_remaining": sum(1 for row in rows if row["expected"] == "true" and row["status"] != "completed"),
        "not_applicable": sum(1 for row in rows if row["status"] == "not_applicable"),
        "batch_rows": batch_rows,
    }
    (output_dir / "completion_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
