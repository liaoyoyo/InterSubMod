#!/usr/bin/env python3
"""compute_metrics.py - Extract metrics from ISM experiment output.

Walks an experiment output directory, finds all metrics.json and
significance_summary.csv files, and produces a unified metrics bundle.

Usage:
    python3 scripts/validation/metrics/compute_metrics.py \
        --experiment-dir /path/to/experiment/output \
        --output metrics_bundle.json
"""

import argparse
import json
import os
import sys
import warnings
from collections import Counter
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    UNKNOWN_CURRENT_CLASS,
    V1_CURRENT_CLASSES,
    SchemaContractError,
    read_evidence,
    select_current_view,
)


PROVENANCE_FIELDS = [
    "VerificationClass_V1_Deprecated",
    "VerificationClass_Legacy",
    "LabelFirstSupport",
    "ClusterFirstSupport",
    "WithinHPSupport",
    "DispersionWarning",
    "EvidencePath",
    "EvidenceDerivation",
]


def find_metrics_files(root_dir):
    """Find all metrics.json files recursively."""
    results = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "metrics.json" in filenames:
            results.append(os.path.join(dirpath, "metrics.json"))
    return sorted(results)


def find_significance_summaries(root_dir):
    """Find all significance_summary.csv files recursively."""
    results = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "significance_summary.csv" in filenames:
            results.append(os.path.join(dirpath, "significance_summary.csv"))
    return sorted(results)


def load_metrics_json(filepath):
    """Load and parse a metrics.json file."""
    with open(filepath, "r") as f:
        return json.load(f)


def _serialized_distribution(series):
    return {
        str(key): int(value)
        for key, value in series.astype("string").fillna("NA").value_counts().sort_index().items()
    }


def _canonical_evidence_distribution(series):
    """Serialize typed evidence booleans as the schema's true/false/NA tokens."""
    values = series.map(
        lambda value: "NA" if pd.isna(value) else ("true" if bool(value) else "false")
    )
    return _serialized_distribution(values)


def _normalize_significance_taxonomy(frame, csv_path, allow_unversioned_v1=False):
    """Validate a current taxonomy without deriving missing evidence from class labels."""
    result = frame.copy()
    if "VerificationSchemaVersion" in result.columns:
        missing = [field for field in PROVENANCE_FIELDS if field not in result.columns]
        if missing:
            raise SchemaContractError(
                f"{csv_path}: missing v2 provenance fields: {', '.join(missing)}"
            )
        current = select_current_view(result)
        evidence = read_evidence(result)
        result["VerificationClass"] = current.values
        metadata = current.metadata()
        metadata["evidence_available"] = True
        metadata["evidence_fields"] = list(PROVENANCE_FIELDS)
        return result, evidence, metadata

    if "VerificationClass" not in result.columns:
        raise SchemaContractError(f"{csv_path}: VerificationClass is missing")
    if not allow_unversioned_v1:
        raise SchemaContractError(
            f"{csv_path}: unversioned current taxonomy requires --allow-unversioned-v1"
        )
    raw = result["VerificationClass"].astype("string")
    valid = raw.isin(V1_CURRENT_CLASSES)
    unknown_counts = {
        str(key): int(value)
        for key, value in raw[~valid].fillna("<MISSING>").value_counts().sort_index().items()
    }
    result["VerificationClass"] = raw.where(valid, UNKNOWN_CURRENT_CLASS)
    message = (
        f"UNVERSIONED: {csv_path} accepted as historical v1 current taxonomy under explicit "
        "--allow-unversioned-v1; v2 evidence provenance is unavailable"
    )
    warnings.warn(message, UserWarning, stacklevel=2)
    return result, None, {
        "selection_field": "VerificationClass",
        "schema_status": "UNVERSIONED_V1",
        "categories": list(V1_CURRENT_CLASSES) + [UNKNOWN_CURRENT_CLASS],
        "unknown_counts": unknown_counts,
        "warnings": [message],
        "evidence_available": False,
        "evidence_fields": [],
    }


def extract_significance_stats(csv_path, allow_unversioned_v1=False):
    """Extract distribution statistics from significance_summary.csv.

    Returns a dict with:
      - verification_class_dist: Counter of VerificationClass values
      - quality_score_stats: mean, median, std of Quality_Score
      - hpfine_ngroups_dist: Counter of HPFineNGroups values
      - dominant_label_dist: Counter of DominantLabel values
      - region_count: total regions in CSV
    """
    try:
        frame = pd.read_csv(csv_path, keep_default_na=False, low_memory=False)
    except (OSError, pd.errors.ParserError) as exc:
        raise SchemaContractError(f"{csv_path}: cannot parse significance summary: {exc}") from exc
    frame, evidence, taxonomy = _normalize_significance_taxonomy(
        frame,
        csv_path,
        allow_unversioned_v1=allow_unversioned_v1,
    )

    verification_classes = Counter(frame["VerificationClass"].astype(str).tolist())
    quality_scores = []
    hpfine_ngroups = Counter()
    dominant_labels = Counter()
    region_count = len(frame.index)

    if "Quality_Score" in frame.columns:
        quality_scores = pd.to_numeric(frame["Quality_Score"], errors="coerce").dropna().tolist()
    if "HPFineNGroups" in frame.columns:
        hp_values = pd.to_numeric(frame["HPFineNGroups"], errors="coerce").dropna()
        hpfine_ngroups.update(int(value) for value in hp_values.tolist())
    if "DominantLabel" in frame.columns:
        dominant_labels.update(value for value in frame["DominantLabel"].astype(str) if value)

    # Compute quality score statistics
    qs_stats = {}
    if quality_scores:
        sorted_qs = sorted(quality_scores)
        n = len(sorted_qs)
        qs_mean = sum(sorted_qs) / n
        qs_median = sorted_qs[n // 2] if n % 2 else (sorted_qs[n // 2 - 1] + sorted_qs[n // 2]) / 2
        qs_std = (sum((x - qs_mean) ** 2 for x in sorted_qs) / n) ** 0.5
        qs_stats = {
            "mean": round(qs_mean, 4),
            "median": round(qs_median, 4),
            "std": round(qs_std, 4),
            "count": n,
        }

    result = {
        "region_count": region_count,
        "verification_class_dist": dict(verification_classes),
        "verification_taxonomy": taxonomy,
        "quality_score_stats": qs_stats,
        "hpfine_ngroups_dist": {str(k): v for k, v in sorted(hpfine_ngroups.items())},
        "dominant_label_dist": dict(dominant_labels),
    }
    if evidence is not None:
        result["verification_class_v1_deprecated_dist"] = _serialized_distribution(
            frame["VerificationClass_V1_Deprecated"]
        )
        result["verification_class_legacy_dist"] = _serialized_distribution(
            frame["VerificationClass_Legacy"]
        )
        boolean_fields = {
            "LabelFirstSupport", "ClusterFirstSupport", "WithinHPSupport", "DispersionWarning"
        }
        result["verification_evidence_dist"] = {
            field: (
                _canonical_evidence_distribution(evidence[field])
                if field in boolean_fields
                else _serialized_distribution(evidence[field])
            )
            for field in evidence.columns
        }
    else:
        result["verification_evidence_dist"] = {"available": False}
    return result


def infer_sample_mode(metrics_data, filepath):
    """Infer sample and mode from metrics.json data or path."""
    sample = metrics_data.get("sample", "unknown")
    mode = metrics_data.get("mode", "unknown")
    return sample, mode


def build_metrics_bundle(experiment_dir, allow_unversioned_v1=False):
    """Build a complete metrics bundle from experiment output."""
    metrics_files = find_metrics_files(experiment_dir)
    sig_files = find_significance_summaries(experiment_dir)

    bundle = {
        "experiment_dir": experiment_dir,
        "samples": {},
        "summary": {
            "total_metrics_files": len(metrics_files),
            "total_significance_files": len(sig_files),
        },
    }

    # Process metrics.json files
    for mf in metrics_files:
        data = load_metrics_json(mf)
        sample, mode = infer_sample_mode(data, mf)

        if sample not in bundle["samples"]:
            bundle["samples"][sample] = {}

        baseline = data.get("baseline", {})
        filtered = data.get("filtered", {})
        improvement = data.get("improvement", {})

        entry = {
            "source_file": mf,
            "baseline_f1": baseline.get("f1"),
            "baseline_precision": baseline.get("precision"),
            "baseline_recall": baseline.get("recall"),
            "baseline_tp": baseline.get("tp"),
            "baseline_fp": baseline.get("fp"),
            "f1": filtered.get("f1"),
            "precision": filtered.get("precision"),
            "recall": filtered.get("recall"),
            "tp": filtered.get("tp"),
            "fp": filtered.get("fp"),
            "f1_delta": improvement.get("f1_delta"),
            "fp_removed": improvement.get("fp_removed"),
            "tp_removed": improvement.get("tp_removed"),
            "truth_total": data.get("truth_total"),
            "tp_regions_analyzed": data.get("tp_regions_analyzed"),
            "fp_regions_analyzed": data.get("fp_regions_analyzed"),
        }

        bundle["samples"][sample][mode] = entry

    # Process significance_summary.csv files (match to TP/FP labels)
    taxonomy_sources = []
    for sf in sig_files:
        sig_stats = extract_significance_stats(
            sf,
            allow_unversioned_v1=allow_unversioned_v1,
        )
        taxonomy_sources.append({"source_file": sf, **sig_stats["verification_taxonomy"]})

        # Try to match to a sample/mode by path
        # Path pattern: .../SAMPLE/MODE/RUN_ID/intersubmod_tp/...
        parts = sf.split(os.sep)
        label = None
        for i, p in enumerate(parts):
            if p.startswith("intersubmod_"):
                label = p  # intersubmod_tp or intersubmod_fp
                break

        if label:
            # Attach significance stats to the matching sample entry
            for sample in bundle["samples"]:
                for mode in bundle["samples"][sample]:
                    key = f"significance_{label}"
                    if key not in bundle["samples"][sample][mode]:
                        bundle["samples"][sample][mode][key] = sig_stats
                        break

    statuses = sorted({item["schema_status"] for item in taxonomy_sources})
    if len(statuses) > 1:
        raise SchemaContractError(
            f"significance summaries mix taxonomy schema statuses: {statuses}"
        )
    bundle["verification_taxonomy"] = {
        "schema_status": statuses[0] if statuses else "NO_SIGNIFICANCE_FILES",
        "selection_field": "VerificationClass" if statuses else None,
        "source_count": len(taxonomy_sources),
        "sources": taxonomy_sources,
    }

    return bundle


def main():
    parser = argparse.ArgumentParser(description="Extract metrics from ISM experiment output")
    parser.add_argument("--experiment-dir", required=True,
                        help="Root directory of experiment output")
    parser.add_argument("--output", default=None,
                        help="Output path for metrics bundle JSON (default: stdout)")
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize historical unversioned v1 current taxonomy summaries.",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.experiment_dir):
        print(f"[ERROR] Experiment directory not found: {args.experiment_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        bundle = build_metrics_bundle(
            args.experiment_dir,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
    except SchemaContractError as exc:
        print(f"[ERROR][schema-contract] {exc}", file=sys.stderr)
        sys.exit(2)

    # Print summary
    sample_count = len(bundle["samples"])
    print(f"[Metrics] Extracted metrics for {sample_count} samples from {args.experiment_dir}",
          file=sys.stderr)
    for sample, modes in sorted(bundle["samples"].items()):
        for mode, data in sorted(modes.items()):
            f1 = data.get("f1")
            f1_str = f"{f1:.4f}" if f1 is not None else "N/A"
            print(f"  {sample}/{mode}: F1={f1_str}", file=sys.stderr)

    # Output
    output_json = json.dumps(bundle, indent=2, ensure_ascii=False)
    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w") as f:
            f.write(output_json)
        print(f"[Metrics] Written to {args.output}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
