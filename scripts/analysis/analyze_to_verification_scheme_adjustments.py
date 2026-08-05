#!/usr/bin/env python3
"""Evaluate alternative verification-class schemes on tumor-only outputs."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    V1_CURRENT_CLASSES,
    read_evidence,
    select_current_view,
)


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def to_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-csv", action="append", required=True, help="significance_summary.csv for TP/FP")
    parser.add_argument("--feature-tsv", required=True, help="to_feature_per_region.tsv")
    parser.add_argument("--benchmark-comparison", required=True, help="benchmark_comparison.tsv")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def load_tsv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def prepare_verification_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate v2 identity and expose evidence plus the frozen v1 reference."""
    current_view = select_current_view(frame)
    evidence = read_evidence(frame)
    v1_field = "VerificationClass_V1_Deprecated"
    if v1_field not in frame.columns:
        raise SchemaContractError(
            f"verification scheme comparison: missing required historical reference {v1_field}"
        )
    v1_values = frame[v1_field].astype("string")
    invalid_mask = v1_values.isna() | ~v1_values.isin(V1_CURRENT_CLASSES)
    if invalid_mask.any():
        invalid = (
            v1_values[invalid_mask]
            .fillna("<MISSING>")
            .value_counts(dropna=False)
            .to_dict()
        )
        raise SchemaContractError(
            f"verification scheme comparison: invalid {v1_field} values: {invalid}"
        )

    prepared = frame.copy()
    prepared["_verification_class_v2"] = current_view.values
    prepared["_verification_class_v1_deprecated"] = v1_values
    prepared["_label_first_support"] = evidence["LabelFirstSupport"]
    prepared["_cluster_first_support"] = evidence["ClusterFirstSupport"]
    prepared["_within_hp_support"] = evidence["WithinHPSupport"]
    prepared["_dispersion_warning"] = evidence["DispersionWarning"]
    prepared["_evidence_path"] = evidence["EvidencePath"]
    prepared["_evidence_derivation"] = evidence["EvidenceDerivation"]
    prepared["_verification_selection_field"] = v1_field
    prepared["_verification_schema_status"] = current_view.schema_status
    return prepared


def load_validated_summary_rows(path: Path) -> List[Dict[str, object]]:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    return prepare_verification_frame(frame).to_dict(orient="records")


def write_tsv(path: Path, fieldnames: Iterable[str], rows: List[Dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def f1_score(tp: int, fp: int, fn: int) -> float:
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    if precision + recall == 0.0:
        return 0.0
    return 2.0 * precision * recall / (precision + recall)


def label_tier(hp_sig: bool, allele_sig: bool) -> str:
    if hp_sig and allele_sig:
        return "Strong"
    if hp_sig:
        return "Subclone"
    if allele_sig:
        return "Weak"
    return "Noise"


def classify_scheme(cluster_support: bool, hp_sig: bool, allele_sig: bool) -> str:
    tier = label_tier(hp_sig, allele_sig)
    if cluster_support:
        if tier == "Strong":
            return "Strong"
        if tier == "Subclone":
            return "Subclone"
        if tier == "Weak":
            return "Weak"
        return "Subclone"
    return tier


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_rows = {
        row["region_key"]: row for row in load_tsv_rows(Path(args.feature_tsv))
    }

    baseline_tp = baseline_fp = baseline_fn = None
    for row in load_tsv_rows(Path(args.benchmark_comparison)):
        if row.get("method") == "InterSubMod":
            baseline_tp = int(row["tp"])
            baseline_fp = int(row["fp"])
            baseline_fn = int(row["fn"])
            break
    if baseline_tp is None or baseline_fp is None or baseline_fn is None:
        raise SystemExit("InterSubMod baseline not found in benchmark comparison.")

    per_region_rows: List[Dict[str, str]] = []
    counts: Dict[str, Dict[str, int]] = {}
    rule_hits: Dict[str, Tuple[int, int]] = {}

    def bump_count(scheme: str, truth_scope: str, cls: str) -> None:
        key = f"{scheme}:{truth_scope}:{cls}"
        counts[key] = counts.get(key, {"count": 0})
        counts[key]["count"] += 1

    def bump_rule(rule_id: str, truth_scope: str) -> None:
        tp_hit, fp_hit = rule_hits.get(rule_id, (0, 0))
        if truth_scope == "tp":
            tp_hit += 1
        elif truth_scope == "fp":
            fp_hit += 1
        rule_hits[rule_id] = (tp_hit, fp_hit)

    for summary_csv in args.summary_csv:
        for row in load_validated_summary_rows(Path(summary_csv)):
            region_key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
            feature = feature_rows.get(region_key)
            if feature is None:
                continue

            truth_scope = feature["source_scope"]
            hp_sig = (
                as_bool(row.get("HPMergedSig", "false"))
                or as_bool(row.get("HPFineSig", "false"))
                or (as_bool(row.get("LabelHPPermanovaValid", "false")) and to_float(row.get("LabelHPPermanovaP", "1")) <= 0.05)
            )
            allele_sig = (
                as_bool(row.get("AlleleSig", "false"))
                or (
                    as_bool(row.get("LabelAllelePermanovaValid", "false"))
                    and to_float(row.get("LabelAllelePermanovaP", "1")) <= 0.05
                )
            )
            global_cluster_support = as_bool(row.get("PassedGating", "false")) and to_float(row.get("GlobalP", "1")) <= 0.05
            local_cluster_support = global_cluster_support and to_float(row.get("LocalBestP", "1")) <= 0.05
            permanova_cluster_support = as_bool(row.get("ClusterPermanovaValid", "false")) and to_float(
                row.get("ClusterPermanovaP", "1")
            ) <= 0.05
            hybrid_cluster_support = global_cluster_support and (
                permanova_cluster_support or to_float(row.get("LocalBestP", "1")) <= 0.05
            )

            reference_class_v1 = row["_verification_class_v1_deprecated"]
            cluster_first_support = bool(row["_cluster_first_support"])
            scheme_label_tiered = classify_scheme(global_cluster_support, hp_sig, allele_sig)
            scheme_local_hybrid = classify_scheme(local_cluster_support, hp_sig, allele_sig)
            scheme_permanova_strict = classify_scheme(permanova_cluster_support, hp_sig, allele_sig)
            scheme_hybrid = classify_scheme(hybrid_cluster_support, hp_sig, allele_sig)
            current_label_tier = label_tier(hp_sig, allele_sig)

            out_row = {
                "region_key": region_key,
                "source_scope": truth_scope,
                "reference_class_v1_deprecated": reference_class_v1,
                "verification_class_v2": row["_verification_class_v2"],
                "verification_selection_field": row["_verification_selection_field"],
                "verification_schema_status": row["_verification_schema_status"],
                "label_first_support": str(bool(row["_label_first_support"])).lower(),
                "cluster_first_support": str(cluster_first_support).lower(),
                "within_hp_support": (
                    "NA" if pd.isna(row["_within_hp_support"])
                    else str(bool(row["_within_hp_support"])).lower()
                ),
                "dispersion_warning": (
                    "NA" if pd.isna(row["_dispersion_warning"])
                    else str(bool(row["_dispersion_warning"])).lower()
                ),
                "evidence_path": row["_evidence_path"],
                "evidence_derivation": row["_evidence_derivation"],
                "label_tier": current_label_tier,
                "scheme_label_tiered": scheme_label_tiered,
                "scheme_local_hybrid": scheme_local_hybrid,
                "scheme_permanova_strict": scheme_permanova_strict,
                "scheme_hybrid": scheme_hybrid,
                "hp_sig": str(hp_sig).lower(),
                "allele_sig": str(allele_sig).lower(),
                "global_cluster_support": str(global_cluster_support).lower(),
                "local_cluster_support": str(local_cluster_support).lower(),
                "permanova_cluster_support": str(permanova_cluster_support).lower(),
                "hybrid_cluster_support": str(hybrid_cluster_support).lower(),
                "vaf": feature["vaf"],
                "allele_delta": feature["allele_delta"],
                "qual": feature["qual"],
                "gq": feature["gq"],
                "pairwise_median_dist": feature["pairwise_median_dist"],
                "low_vaf_high_ad": feature["low_vaf_high_ad"],
                "low_vaf_high_ad_cv": feature["low_vaf_high_ad_cv"],
            }
            per_region_rows.append(out_row)

            for scheme_name, cls in (
                ("v1_deprecated_reference", reference_class_v1),
                ("label_tiered", scheme_label_tiered),
                ("local_hybrid", scheme_local_hybrid),
                ("permanova_strict", scheme_permanova_strict),
                ("hybrid", scheme_hybrid),
            ):
                bump_count(scheme_name, truth_scope, cls)

            low_vaf_high_ad = feature["low_vaf_high_ad"] == "True"
            low_vaf_high_ad_cv = feature["low_vaf_high_ad_cv"] == "True"
            # Historical rule id is retained for output compatibility, but the
            # trigger now reads explicit cluster-first evidence rather than
            # reverse-engineering support from a class name.
            current_strong_allele_only = cluster_first_support and (not hp_sig) and allele_sig
            current_strong_allele_only_low = current_strong_allele_only and low_vaf_high_ad
            current_strong_allele_only_low_cv = current_strong_allele_only and low_vaf_high_ad_cv
            tiered_weak_low = scheme_label_tiered == "Weak" and low_vaf_high_ad
            hybrid_weak_low_cv = scheme_hybrid == "Weak" and low_vaf_high_ad_cv

            for rule_id, hit in (
                ("current_low_vaf_high_ad", low_vaf_high_ad),
                ("current_low_vaf_high_ad_cv", low_vaf_high_ad_cv),
                ("current_strong_allele_only", current_strong_allele_only),
                ("current_strong_allele_only_low", current_strong_allele_only_low),
                ("current_strong_allele_only_low_cv", current_strong_allele_only_low_cv),
                ("tiered_weak_low", tiered_weak_low),
                ("hybrid_weak_low_cv", hybrid_weak_low_cv),
            ):
                if hit:
                    bump_rule(rule_id, truth_scope)

    write_tsv(
        output_dir / "verification_scheme_per_region.tsv",
        [
            "region_key",
            "source_scope",
            "reference_class_v1_deprecated",
            "verification_class_v2",
            "verification_selection_field",
            "verification_schema_status",
            "label_first_support",
            "cluster_first_support",
            "within_hp_support",
            "dispersion_warning",
            "evidence_path",
            "evidence_derivation",
            "label_tier",
            "scheme_label_tiered",
            "scheme_local_hybrid",
            "scheme_permanova_strict",
            "scheme_hybrid",
            "hp_sig",
            "allele_sig",
            "global_cluster_support",
            "local_cluster_support",
            "permanova_cluster_support",
            "hybrid_cluster_support",
            "vaf",
            "allele_delta",
            "qual",
            "gq",
            "pairwise_median_dist",
            "low_vaf_high_ad",
            "low_vaf_high_ad_cv",
        ],
        per_region_rows,
    )

    class_rows: List[Dict[str, str]] = []
    for key, payload in sorted(counts.items()):
        scheme, truth_scope, cls = key.split(":", 2)
        class_rows.append(
            {
                "scheme": scheme,
                "source_scope": truth_scope,
                "class": cls,
                "count": str(payload["count"]),
            }
        )
    write_tsv(output_dir / "verification_scheme_class_counts.tsv", ["scheme", "source_scope", "class", "count"], class_rows)

    rule_rows: List[Dict[str, str]] = []
    baseline_f1 = f1_score(baseline_tp, baseline_fp, baseline_fn)
    for rule_id, (tp_hit, fp_hit) in sorted(rule_hits.items()):
        new_tp = baseline_tp - tp_hit
        new_fp = baseline_fp - fp_hit
        new_fn = baseline_fn + tp_hit
        new_f1 = f1_score(new_tp, new_fp, new_fn)
        rule_rows.append(
            {
                "rule_id": rule_id,
                "tp_removed": str(tp_hit),
                "fp_removed": str(fp_hit),
                "fp_per_tp": f"{(fp_hit / tp_hit):.6f}" if tp_hit else "",
                "new_tp": str(new_tp),
                "new_fp": str(new_fp),
                "new_fn": str(new_fn),
                "baseline_f1": f"{baseline_f1:.6f}",
                "new_f1": f"{new_f1:.6f}",
                "delta_f1": f"{(new_f1 - baseline_f1):.6f}",
            }
        )
    write_tsv(
        output_dir / "verification_scheme_rule_comparison.tsv",
        ["rule_id", "tp_removed", "fp_removed", "fp_per_tp", "new_tp", "new_fp", "new_fn", "baseline_f1", "new_f1", "delta_f1"],
        rule_rows,
    )

    summary_lines = [
        "# Verification Scheme Adjustments",
        "",
        f"- InterSubMod baseline: TP={baseline_tp}, FP={baseline_fp}, FN={baseline_fn}, F1={baseline_f1:.6f}",
        "",
        "## Key Points",
        "",
        "- `current_strong_allele_only` 為保留歷史比較而維持的 rule id；實際條件已改讀明確的 `ClusterFirstSupport=true`，且只有 allele signal、沒有 HP signal。",
        "- class A/B 的參考欄固定為 `VerificationClass_V1_Deprecated`；canonical v2 class 與 evidence provenance 逐列輸出，不做 v2→v1 折疊。",
        "- `label_tiered` 將 `Strong` 僅留給 `hp_sig + allele_sig`；只有 allele signal 的情況降為 `Weak`。",
        "- `hybrid` 在 `label_tiered` 基礎上，另外要求 cluster support 至少通過 `GlobalP + LocalBestP/ClusterPermanova` 的混合條件。",
        "",
        "## Rule Snapshot",
        "",
        "| rule_id | tp_removed | fp_removed | fp_per_tp | new_f1 | delta_f1 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rule_rows:
        summary_lines.append(
            f"| {row['rule_id']} | {row['tp_removed']} | {row['fp_removed']} | {row['fp_per_tp']} | {row['new_f1']} | {row['delta_f1']} |"
        )
    (output_dir / "verification_scheme_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
