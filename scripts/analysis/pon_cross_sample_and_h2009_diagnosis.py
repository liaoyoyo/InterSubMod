#!/usr/bin/env python3
"""PON cross-sample removal rate + H2009 root cause diagnosis.

Task 2: PON removal rate across all 7 samples from raw TO VCFs.
Task 3: H2009 negative root cause analysis using master dataset.

Output: TSV tables + diagnostic report data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    UNKNOWN_CURRENT_CLASS,
    V1_CURRENT_CLASSES,
    SchemaContractError,
    read_evidence,
    select_current_view,
)

# ── Config ──────────────────────────────────────────────────────────────

MASTER_TSV = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces"
    "/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
)

OUTPUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/independent_analyses_20260411")

# Pre-computed raw VCF stats (from shell analysis)
RAW_VCF_STATS = {
    "HCC1395":       {"total": 3187275, "nonsomatic": 3133474, "pass": 47798, "lowqual_only": 6003},
    "HCC1395_DORADO":{"total": 3167467, "nonsomatic": 3114957, "pass": 48072, "lowqual_only": 4438},
    "COLO829":       {"total": 3444386, "nonsomatic": 3383636, "pass": 50991, "lowqual_only": 9759},
    "H1437":         {"total": 3415364, "nonsomatic": 3339874, "pass": 65448, "lowqual_only": 10042},
    "H2009":         {"total": 3447798, "nonsomatic": 3280872, "pass": 153854, "lowqual_only": 13072},
    "HCC1937":       {"total": 3124945, "nonsomatic": 3087734, "pass": 29503, "lowqual_only": 7708},
    "HCC1954":       {"total": 3775749, "nonsomatic": 3690623, "pass": 74996, "lowqual_only": 10130},
}

TRUTH_TOTALS = {
    "HCC1395": 39447, "HCC1395_DORADO": 39447, "COLO829": 41427,
    "H1437": 90129, "H2009": 168529, "HCC1937": 60691, "HCC1954": 26567,
}

# ISM features for H2009 comparison
KEY_FEATURES = [
    "CramersV", "PairwiseMeanDist", "PairwiseMedianDist", "Quality_Score",
    "HPMergedDelta", "HP_Ratio", "AlleleDelta", "NumReads", "NumCpGs",
    "ClusterPermanovaF", "HPFineNGroups", "HP1FamilyN", "HP2FamilyN",
    "LabelAllelePermanovaF", "LabelHPPermanovaF",
]


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-tsv", default=str(MASTER_TSV))
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR))
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize the historical unversioned v1 current taxonomy.",
    )
    return parser.parse_args()


def _validate_current_taxonomy(
    rows: List[Dict[str, str]],
    headers: List[str],
    allow_unversioned_v1: bool,
) -> Dict[str, object]:
    """Fail closed on schema identity and normalize only unknown class values."""
    required = ["VerificationClass"]
    if "VerificationSchemaVersion" in headers:
        required += PROVENANCE_FIELDS + ["VerificationSchemaVersion"]
        missing = [field for field in required if field not in headers]
        if missing:
            raise SchemaContractError(
                "PON master dataset: missing v2 provenance fields: " + ", ".join(missing)
            )
        schema_frame = pd.DataFrame(
            ({field: row.get(field, "") for field in required} for row in rows),
            columns=required,
        )
        current = select_current_view(schema_frame)
        read_evidence(schema_frame)
        for row, value in zip(rows, current.values.tolist()):
            row["VerificationClass"] = value
        metadata = current.metadata()
        metadata["evidence_fields"] = list(PROVENANCE_FIELDS)
        return metadata

    if "VerificationClass" not in headers:
        raise SchemaContractError("PON master dataset: VerificationClass is missing")
    if not allow_unversioned_v1:
        raise SchemaContractError(
            "PON master dataset is unversioned; --allow-unversioned-v1 is required"
        )
    unknown_counts: Counter[str] = Counter()
    allowed = set(V1_CURRENT_CLASSES)
    for row in rows:
        value = row.get("VerificationClass", "")
        if value not in allowed:
            unknown_counts[value or "<MISSING>"] += 1
            row["VerificationClass"] = UNKNOWN_CURRENT_CLASS
    message = (
        "UNVERSIONED: PON master current taxonomy accepted under explicit "
        "--allow-unversioned-v1; v2 evidence provenance is unavailable"
    )
    warnings.warn(message, UserWarning, stacklevel=2)
    return {
        "selection_field": "VerificationClass",
        "schema_status": "UNVERSIONED_V1",
        "categories": list(V1_CURRENT_CLASSES) + [UNKNOWN_CURRENT_CLASS],
        "unknown_counts": dict(sorted(unknown_counts.items())),
        "warnings": [message],
        "evidence_fields": [],
    }


def load_master_dataset(
    path: Path = MASTER_TSV,
    allow_unversioned_v1: bool = False,
) -> Tuple[List[str], List[Dict[str, str]], Dict[str, object]]:
    """Load and schema-check the master TSV before returning any analytical rows."""
    print("Loading master dataset...", file=sys.stderr)
    rows = []
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = list(reader.fieldnames or [])
        rows.extend(reader)
    print(f"  Loaded {len(rows):,} rows, {len(headers)} columns", file=sys.stderr)
    metadata = _validate_current_taxonomy(rows, headers, allow_unversioned_v1)
    return headers, rows, metadata


def _taxonomy_header(metadata: Dict[str, object]) -> List[str]:
    return [
        f"Verification taxonomy status: {metadata['schema_status']}",
        f"Verification selection field: {metadata['selection_field']}",
        f"Unknown current classes bucketed: {metadata.get('unknown_counts', {})}",
        "",
    ]


def safe_float(v: str) -> float:
    try:
        return float(v)
    except (ValueError, TypeError):
        return float("nan")


# ── Task 2: PON Cross-Sample Analysis ──────────────────────────────────

def run_pon_analysis(rows: List[Dict[str, str]], metadata: Dict[str, object] | None = None) -> str:
    """Compute PON removal rates and cross-sample comparison."""
    output_lines = []
    output_lines.append("=" * 80)
    output_lines.append("Task 2: PON Cross-Sample FP Removal Rate Analysis")
    output_lines.append("=" * 80)
    output_lines.append("")
    if metadata is not None:
        output_lines.extend(_taxonomy_header(metadata))

    # Part A: Raw VCF-level PON removal
    output_lines.append("## Part A: Raw VCF-Level PON Removal Rate")
    output_lines.append("")
    output_lines.append(f"{'Sample':<20} {'Total':>10} {'NonSomatic':>12} {'PASS':>8} "
                        f"{'LQ_only':>8} {'PON_Rate':>10} {'PASS/Total':>10}")
    output_lines.append("-" * 90)

    pon_rates = {}
    for sample in ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]:
        s = RAW_VCF_STATS[sample]
        pon_rate = s["nonsomatic"] / s["total"] * 100
        pass_rate = s["pass"] / s["total"] * 100
        pon_rates[sample] = pon_rate
        output_lines.append(
            f"{sample:<20} {s['total']:>10,} {s['nonsomatic']:>12,} {s['pass']:>8,} "
            f"{s['lowqual_only']:>8,} {pon_rate:>9.2f}% {pass_rate:>9.2f}%"
        )

    output_lines.append("")
    rates = list(pon_rates.values())
    output_lines.append(f"PON Rate Range: {min(rates):.2f}% - {max(rates):.2f}%")
    output_lines.append(f"PON Rate Mean: {np.mean(rates):.2f}% (SD: {np.std(rates):.2f}%)")
    output_lines.append(f"PON Rate Median: {np.median(rates):.2f}%")
    output_lines.append("")

    # Part B: ISM-level TP/FP counts per sample (from master dataset)
    output_lines.append("## Part B: ISM-Level TP/FP Distribution (Master Dataset)")
    output_lines.append("")

    sample_stats = defaultdict(lambda: {"TP": 0, "FP": 0})
    for row in rows:
        sample = row.get("sample", "")
        mode = row.get("mode", "")
        truth = row.get("truth_label", "")
        if truth in ("TP", "FP"):
            key = f"{sample}_{mode}"
            sample_stats[key][truth] += 1

    output_lines.append(f"{'Dataset':<35} {'TP':>8} {'FP':>8} {'Total':>8} {'FP_Rate':>8}")
    output_lines.append("-" * 70)
    for key in sorted(sample_stats.keys()):
        s = sample_stats[key]
        total = s["TP"] + s["FP"]
        fp_rate = s["FP"] / total * 100 if total > 0 else 0
        output_lines.append(f"{key:<35} {s['TP']:>8,} {s['FP']:>8,} {total:>8,} {fp_rate:>7.1f}%")

    output_lines.append("")

    # Part C: PON as fraction of non-TP (refined calculation)
    output_lines.append("## Part C: PON Removal Rate of FP (Refined)")
    output_lines.append("")
    output_lines.append("  PON_FP_Rate = NonSomatic / (Total - truth_total)")
    output_lines.append("  Note: upper bound estimate (not all truth TP are in raw VCF)")
    output_lines.append("")
    output_lines.append(f"{'Sample':<20} {'Total':>10} {'Truth_Total':>12} {'Est_FP':>10} "
                        f"{'NonSomatic':>12} {'PON_FP_Rate':>12}")
    output_lines.append("-" * 80)

    for sample in ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]:
        s = RAW_VCF_STATS[sample]
        tt = TRUTH_TOTALS[sample]
        est_fp = s["total"] - tt
        pon_fp_rate = s["nonsomatic"] / est_fp * 100 if est_fp > 0 else 0
        output_lines.append(
            f"{sample:<20} {s['total']:>10,} {tt:>12,} {est_fp:>10,} "
            f"{s['nonsomatic']:>12,} {pon_fp_rate:>11.2f}%"
        )

    output_lines.append("")
    output_lines.append("## Part D: Key Observations")
    output_lines.append("")
    output_lines.append("1. H2009 has significantly lower PON rate (95.16%) vs others (97.75-98.81%)")
    output_lines.append("2. H2009 PASS count (153,854) is 2-5x higher than other samples")
    output_lines.append("3. H2009 truth_total (168,529) is highest — more somatic mutations to detect")
    output_lines.append("4. All samples PON rate > 95% — qualitative conclusion stable")
    output_lines.append("")

    return "\n".join(output_lines)


# ── Task 3: H2009 Root Cause Diagnosis ─────────────────────────────────

def run_h2009_diagnosis(rows: List[Dict[str, str]], metadata: Dict[str, object] | None = None) -> str:
    """Diagnose H2009 negative root cause."""
    output_lines = []
    output_lines.append("=" * 80)
    output_lines.append("Task 3: H2009 Negative Root Cause Diagnosis")
    output_lines.append("=" * 80)
    output_lines.append("")
    if metadata is not None:
        output_lines.extend(_taxonomy_header(metadata))

    # Separate by sample and mode
    sample_data = defaultdict(lambda: {"TP": [], "FP": []})
    for row in rows:
        sample = row.get("sample", "")
        mode = row.get("mode", "")
        truth = row.get("truth_label", "")
        if truth in ("TP", "FP") and mode == "paired":
            sample_data[sample][truth].append(row)

    # Part A: Basic stats comparison
    output_lines.append("## Part A: Sample-Level Statistics (Paired Mode)")
    output_lines.append("")
    output_lines.append(f"{'Sample':<20} {'TP':>7} {'FP':>7} {'FP_Rate':>8} "
                        f"{'LOH_TP%':>8} {'LOH_FP%':>8} {'LOH_FP/TP':>10}")
    output_lines.append("-" * 75)

    loh_stats = {}
    for sample in sorted(sample_data.keys()):
        tp_rows = sample_data[sample]["TP"]
        fp_rows = sample_data[sample]["FP"]
        n_tp = len(tp_rows)
        n_fp = len(fp_rows)
        total = n_tp + n_fp
        fp_rate = n_fp / total * 100 if total > 0 else 0

        # LOH counts
        tp_loh = sum(1 for r in tp_rows if r.get("Potential_LOH", "") == "True")
        fp_loh = sum(1 for r in fp_rows if r.get("Potential_LOH", "") == "True")
        tp_loh_pct = tp_loh / n_tp * 100 if n_tp > 0 else 0
        fp_loh_pct = fp_loh / n_fp * 100 if n_fp > 0 else 0
        loh_ratio = (fp_loh_pct / tp_loh_pct) if tp_loh_pct > 0 else float("inf")

        loh_stats[sample] = {
            "tp_loh_pct": tp_loh_pct, "fp_loh_pct": fp_loh_pct,
            "tp_loh": tp_loh, "fp_loh": fp_loh, "n_tp": n_tp, "n_fp": n_fp,
        }

        output_lines.append(
            f"{sample:<20} {n_tp:>7,} {n_fp:>7,} {fp_rate:>7.1f}% "
            f"{tp_loh_pct:>7.1f}% {fp_loh_pct:>7.1f}% {loh_ratio:>9.2f}x"
        )

    output_lines.append("")

    # Part B: AF distribution comparison
    output_lines.append("## Part B: AF Distribution — High AF (>=0.9) FP Concentration")
    output_lines.append("")
    output_lines.append(f"{'Sample':<20} {'FP_total':>8} {'FP_AF>=0.9':>10} {'FP_AF>=0.9%':>12} "
                        f"{'TP_AF>=0.9':>10} {'TP_AF>=0.9%':>12}")
    output_lines.append("-" * 75)

    for sample in sorted(sample_data.keys()):
        tp_rows = sample_data[sample]["TP"]
        fp_rows = sample_data[sample]["FP"]
        n_fp = len(fp_rows)
        n_tp = len(tp_rows)

        fp_hi_af = sum(1 for r in fp_rows if safe_float(r.get("caller_af", "")) >= 0.9)
        tp_hi_af = sum(1 for r in tp_rows if safe_float(r.get("caller_af", "")) >= 0.9)
        fp_hi_pct = fp_hi_af / n_fp * 100 if n_fp > 0 else 0
        tp_hi_pct = tp_hi_af / n_tp * 100 if n_tp > 0 else 0

        output_lines.append(
            f"{sample:<20} {n_fp:>8,} {fp_hi_af:>10,} {fp_hi_pct:>11.1f}% "
            f"{tp_hi_af:>10,} {tp_hi_pct:>11.1f}%"
        )

    output_lines.append("")

    # Part C: LOH.bed hit comparison
    output_lines.append("## Part C: LOH.bed Hit Rate (to_loh_bed_hit)")
    output_lines.append("")
    output_lines.append(f"{'Sample':<20} {'TP_LOH_bed%':>12} {'FP_LOH_bed%':>12} {'Ratio':>8}")
    output_lines.append("-" * 55)

    for sample in sorted(sample_data.keys()):
        tp_rows = sample_data[sample]["TP"]
        fp_rows = sample_data[sample]["FP"]
        n_tp = len(tp_rows)
        n_fp = len(fp_rows)

        tp_loh_bed = sum(1 for r in tp_rows if r.get("to_loh_bed_hit", "") == "True")
        fp_loh_bed = sum(1 for r in fp_rows if r.get("to_loh_bed_hit", "") == "True")
        tp_pct = tp_loh_bed / n_tp * 100 if n_tp > 0 else 0
        fp_pct = fp_loh_bed / n_fp * 100 if n_fp > 0 else 0
        ratio = fp_pct / tp_pct if tp_pct > 0 else float("inf")

        output_lines.append(f"{sample:<20} {tp_pct:>11.1f}% {fp_pct:>11.1f}% {ratio:>7.2f}x")

    output_lines.append("")

    # Part D: Key ISM feature AUC comparison (paired mode)
    output_lines.append("## Part D: Key ISM Feature Median Comparison (Paired, TP vs FP)")
    output_lines.append("")

    # Compute per-sample TP/FP medians for key features
    feature_comparison = {}
    for feat in KEY_FEATURES:
        feature_comparison[feat] = {}
        for sample in sorted(sample_data.keys()):
            tp_vals = [safe_float(r.get(feat, "")) for r in sample_data[sample]["TP"]]
            fp_vals = [safe_float(r.get(feat, "")) for r in sample_data[sample]["FP"]]
            tp_vals = [v for v in tp_vals if not np.isnan(v)]
            fp_vals = [v for v in fp_vals if not np.isnan(v)]
            tp_med = np.median(tp_vals) if tp_vals else float("nan")
            fp_med = np.median(fp_vals) if fp_vals else float("nan")
            feature_comparison[feat][sample] = {"tp_med": tp_med, "fp_med": fp_med}

    # Print top features with biggest H2009 anomaly
    output_lines.append(f"{'Feature':<25} {'H2009_TP':>10} {'H2009_FP':>10} {'H2009_Gap':>10} "
                        f"{'Others_Gap':>10} {'Anomaly':>10}")
    output_lines.append("-" * 80)

    for feat in KEY_FEATURES:
        h2009 = feature_comparison[feat].get("H2009", {})
        h2009_tp = h2009.get("tp_med", float("nan"))
        h2009_fp = h2009.get("fp_med", float("nan"))
        h2009_gap = h2009_tp - h2009_fp if not (np.isnan(h2009_tp) or np.isnan(h2009_fp)) else float("nan")

        # Average gap for other samples
        other_gaps = []
        for sample in sorted(sample_data.keys()):
            if sample == "H2009":
                continue
            s = feature_comparison[feat].get(sample, {})
            tp_m = s.get("tp_med", float("nan"))
            fp_m = s.get("fp_med", float("nan"))
            if not (np.isnan(tp_m) or np.isnan(fp_m)):
                other_gaps.append(tp_m - fp_m)
        other_gap_mean = np.mean(other_gaps) if other_gaps else float("nan")
        anomaly = h2009_gap - other_gap_mean if not (np.isnan(h2009_gap) or np.isnan(other_gap_mean)) else float("nan")

        output_lines.append(
            f"{feat:<25} {h2009_tp:>10.4f} {h2009_fp:>10.4f} {h2009_gap:>10.4f} "
            f"{other_gap_mean:>10.4f} {anomaly:>10.4f}"
        )

    output_lines.append("")

    # Part E: Verification Class distribution
    output_lines.append("## Part E: VerificationClass Distribution (Paired)")
    output_lines.append("")

    vc_dist = defaultdict(lambda: defaultdict(lambda: {"TP": 0, "FP": 0}))
    for sample in sorted(sample_data.keys()):
        for truth in ["TP", "FP"]:
            for row in sample_data[sample][truth]:
                vc = row["VerificationClass"]
                vc_dist[sample][vc][truth] += 1

    for sample in ["H2009", "HCC1395", "COLO829"]:
        output_lines.append(f"  {sample}:")
        total_tp = sum(d["TP"] for d in vc_dist[sample].values())
        total_fp = sum(d["FP"] for d in vc_dist[sample].values())
        for vc in sorted(vc_dist[sample].keys()):
            tp_n = vc_dist[sample][vc]["TP"]
            fp_n = vc_dist[sample][vc]["FP"]
            tp_pct = tp_n / total_tp * 100 if total_tp > 0 else 0
            fp_pct = fp_n / total_fp * 100 if total_fp > 0 else 0
            output_lines.append(f"    {vc:<15} TP: {tp_n:>6} ({tp_pct:>5.1f}%)  FP: {fp_n:>6} ({fp_pct:>5.1f}%)")
        output_lines.append("")

    # Part F: H2009 root cause summary
    output_lines.append("## Part F: Root Cause Summary")
    output_lines.append("")

    h2009_loh = loh_stats.get("H2009", {})
    output_lines.append("H2009 Unique Characteristics:")
    output_lines.append(f"  1. PASS variants: 153,854 (2-5x other samples)")
    output_lines.append(f"  2. Truth total: 168,529 (highest — lung cancer, high mutation burden)")
    output_lines.append(f"  3. PON rate: 95.16% (lowest — more rare germline variants escape PON)")
    output_lines.append(f"  4. ISM LOH TP%: {h2009_loh.get('tp_loh_pct', 0):.1f}%, "
                        f"FP%: {h2009_loh.get('fp_loh_pct', 0):.1f}%")
    output_lines.append(f"  5. FP in LOH zones: {h2009_loh.get('fp_loh', 0):,} "
                        f"(out of {h2009_loh.get('n_fp', 0):,})")
    output_lines.append("")
    output_lines.append("Root Cause Hypothesis:")
    output_lines.append("  H2009 is a lung cancer line with the highest mutation burden (168K truth SNVs).")
    output_lines.append("  This results in:")
    output_lines.append("  a) More germline variants escaping PON (lower PON rate)")
    output_lines.append("  b) Higher PASS count with more FP among them")
    output_lines.append("  c) FP heavily concentrated in LOH regions")
    output_lines.append("  d) Methylation features add noise in LOH-dominant FP environment")
    output_lines.append("  e) Phase 1A methyl+context degradation: features trained on")
    output_lines.append("     other samples' TP/FP distributions don't generalize to H2009's")
    output_lines.append("     LOH-dominated FP landscape")
    output_lines.append("")

    return "\n".join(output_lines)


# ── Main ────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()

    try:
        headers, rows, verification_metadata = load_master_dataset(
            Path(args.master_tsv).resolve(),
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
    except SchemaContractError as exc:
        raise SystemExit(f"[pon-diagnosis][schema-contract] {exc}") from exc
    output_dir.mkdir(parents=True, exist_ok=True)

    # Task 2
    print("\n--- Running Task 2: PON Cross-Sample Analysis ---\n", file=sys.stderr)
    pon_report = run_pon_analysis(rows, verification_metadata)
    pon_path = output_dir / "pon_cross_sample_analysis.txt"
    pon_path.write_text(pon_report)
    print(pon_report)

    print("\n\n", file=sys.stderr)

    # Task 3
    print("--- Running Task 3: H2009 Diagnosis ---\n", file=sys.stderr)
    h2009_report = run_h2009_diagnosis(rows, verification_metadata)
    h2009_path = output_dir / "h2009_root_cause_diagnosis.txt"
    h2009_path.write_text(h2009_report)
    print(h2009_report)

    # Save TSV summary for PON rates
    tsv_path = output_dir / "pon_removal_rates.tsv"
    with open(tsv_path, "w") as f:
        f.write(
            "sample\ttotal_raw\tnonsomatic\tpass\tlowqual_only\tpon_rate_pct\ttruth_total\t"
            "pass_to_total_pct\tVerificationSchemaStatus\tVerificationSchemaVersion\n"
        )
        for sample in ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]:
            s = RAW_VCF_STATS[sample]
            pon_rate = s["nonsomatic"] / s["total"] * 100
            pass_rate = s["pass"] / s["total"] * 100
            f.write(f"{sample}\t{s['total']}\t{s['nonsomatic']}\t{s['pass']}\t"
                    f"{s['lowqual_only']}\t{pon_rate:.4f}\t{TRUTH_TOTALS[sample]}\t{pass_rate:.4f}\t"
                    f"{verification_metadata['schema_status']}\t"
                    f"{2 if verification_metadata['schema_status'] == 'V2' else ''}\n")

    print(f"\nOutputs saved to {output_dir}/", file=sys.stderr)
    print(f"  - pon_cross_sample_analysis.txt", file=sys.stderr)
    print(f"  - h2009_root_cause_diagnosis.txt", file=sys.stderr)
    print(f"  - pon_removal_rates.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
