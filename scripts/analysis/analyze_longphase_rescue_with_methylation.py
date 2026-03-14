#!/usr/bin/env python3
"""Evaluate LongPhase TP-rescue rules combining caller and methylation features."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Callable, Dict, List

import pandas as pd

from research_common import compute_metrics, write_tsv_rows


RuleFunc = Callable[[pd.Series], bool]


RULE_FIELDS = [
    "rule_id",
    "trigger_count",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_longphase",
    "baseline_f1",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rescue-variants-tsv", required=True, help="Output from analyze_longphase_rescue.py")
    parser.add_argument("--tp-summary-csv", required=True, help="InterSubMod summary for caller_lost_tp")
    parser.add_argument("--fp-summary-csv", required=True, help="InterSubMod summary for caller_removed_fp")
    parser.add_argument("--baseline-tp", type=int, required=True)
    parser.add_argument("--baseline-fp", type=int, required=True)
    parser.add_argument("--truth-total", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def build_rules() -> Dict[str, RuleFunc]:
    return {
        "gq_ge_20": lambda row: row["gq"] >= 20,
        "pairwise_ge_020": lambda row: row["PairwiseMedianDist"] >= 0.20,
        "class_strong_or_subclone": lambda row: row["VerificationClass"] in {"Strong", "Subclone"},
        "class_strong_only": lambda row: row["VerificationClass"] == "Strong",
        "class_subclone_only": lambda row: row["VerificationClass"] == "Subclone",
        "gq20_and_pairwise_ge_020": lambda row: row["gq"] >= 20 and row["PairwiseMedianDist"] >= 0.20,
        "gq20_and_strong_or_subclone": lambda row: row["gq"] >= 20 and row["VerificationClass"] in {"Strong", "Subclone"},
        "gq20_and_strong_or_subclone_and_af_ge_020": lambda row: row["gq"] >= 20
        and row["VerificationClass"] in {"Strong", "Subclone"}
        and row["af"] >= 0.20,
        "strong_or_subclone_and_af_ge_020": lambda row: row["VerificationClass"] in {"Strong", "Subclone"}
        and row["af"] >= 0.20,
        "low_gq_and_strong_or_subclone": lambda row: row["gq"] < 20 and row["VerificationClass"] in {"Strong", "Subclone"},
        "low_gq_and_strong_or_subclone_and_af_ge_020": lambda row: row["gq"] < 20
        and row["VerificationClass"] in {"Strong", "Subclone"}
        and row["af"] >= 0.20,
        "low_gq_and_pairwise_ge_020_and_allele_delta_ge_005": lambda row: row["gq"] < 20
        and row["PairwiseMedianDist"] >= 0.20
        and row["AlleleDelta"] >= 0.05,
    }


def load_joined(args: argparse.Namespace) -> pd.DataFrame:
    rescue = pd.read_csv(args.rescue_variants_tsv, sep="\t")
    tp_summary = pd.read_csv(args.tp_summary_csv)
    fp_summary = pd.read_csv(args.fp_summary_csv)

    for df in (tp_summary, fp_summary):
        df["region_key"] = df["Chr"].astype(str) + ":" + df["Pos"].astype(str) + ":" + df["Ref"].astype(str) + ":" + df["Alt"].astype(str)

    sig = pd.concat([tp_summary, fp_summary], ignore_index=True)
    columns = [
        "region_key",
        "VerificationClass",
        "DominantLabel",
        "AlleleDelta",
        "CramersV",
        "PairwiseMedianDist",
        "GlobalP",
        "NumReads",
    ]
    joined = rescue.merge(sig[columns], on="region_key", how="left")
    joined["VerificationClass"] = joined["VerificationClass"].fillna("NA")
    joined["DominantLabel"] = joined["DominantLabel"].fillna("none")
    for column in ("AlleleDelta", "CramersV", "PairwiseMedianDist", "GlobalP"):
        joined[column] = pd.to_numeric(joined[column], errors="coerce").fillna(0.0)
    joined["NumReads"] = pd.to_numeric(joined["NumReads"], errors="coerce").fillna(0).astype(int)
    joined["af"] = pd.to_numeric(joined["af"], errors="coerce").fillna(0.0)
    joined["gq"] = pd.to_numeric(joined["gq"], errors="coerce").fillna(0).astype(int)
    return joined


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    joined = load_joined(args)
    joined.to_csv(output_dir / "joined_rescue_features.tsv", sep="\t", index=False)

    tp = joined[joined["variant_group"] == "caller_lost_tp"].copy()
    fp = joined[joined["variant_group"] == "caller_removed_fp"].copy()

    baseline_metrics = compute_metrics(args.baseline_tp, args.baseline_fp, args.truth_total)
    rows: List[Dict[str, object]] = []
    best_rule_id = "baseline"
    best_f1 = baseline_metrics["f1"]
    best_delta = 0.0

    for rule_id, rule in build_rules().items():
        tp_rescued = int(tp.apply(rule, axis=1).sum())
        fp_reintroduced = int(fp.apply(rule, axis=1).sum())
        metrics = compute_metrics(args.baseline_tp + tp_rescued, args.baseline_fp + fp_reintroduced, args.truth_total)
        delta = metrics["f1"] - baseline_metrics["f1"]
        if delta > best_delta:
            best_delta = delta
            best_rule_id = rule_id
            best_f1 = metrics["f1"]
        notes: List[str] = []
        if "gq" in rule_id:
            notes.append("caller_gq")
        if "pairwise" in rule_id:
            notes.append("methyl_pairwise")
        if "strong" in rule_id or "subclone" in rule_id:
            notes.append("label_first_class")
        if "allele_delta" in rule_id:
            notes.append("methyl_allele_delta")
        rows.append(
            {
                "rule_id": rule_id,
                "trigger_count": tp_rescued + fp_reintroduced,
                "tp_rescued": tp_rescued,
                "fp_reintroduced": fp_reintroduced,
                "precision": f"{metrics['precision']:.6f}",
                "recall": f"{metrics['recall']:.6f}",
                "f1": f"{metrics['f1']:.6f}",
                "delta_f1_vs_longphase": f"{delta:.6f}",
                "baseline_f1": f"{baseline_metrics['f1']:.6f}",
                "notes": "; ".join(notes),
            }
        )

    rows = sorted(rows, key=lambda row: float(row["delta_f1_vs_longphase"]), reverse=True)
    write_tsv_rows(output_dir / "longphase_rescue_methylation_rule_comparison.tsv", RULE_FIELDS, rows)

    md_lines = [
        "# LongPhase Rescue With Methylation Summary",
        "",
        f"- baseline LongPhase-S F1: `{baseline_metrics['f1']:.6f}`",
        f"- best rule: `{best_rule_id}`",
        f"- best F1: `{best_f1:.6f}`",
        f"- delta: `{best_delta:.6f}`",
        "",
        "## Top rules",
        "",
        "| rule_id | tp_rescued | fp_reintroduced | f1 | delta_f1_vs_longphase | notes |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows[:10]:
        md_lines.append(
            f"| {row['rule_id']} | {row['tp_rescued']} | {row['fp_reintroduced']} | {row['f1']} | {row['delta_f1_vs_longphase']} | {row['notes']} |"
        )
    (output_dir / "longphase_rescue_methylation_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[analyze_longphase_rescue_with_methylation] Wrote {output_dir / 'longphase_rescue_methylation_rule_comparison.tsv'}")


if __name__ == "__main__":
    main()
