#!/usr/bin/env python3
"""Evaluate caller-first TP rescue rules with optional methylation veto/support."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Dict, List

import pandas as pd

from research_common import compute_metrics, infer_platform, write_tsv_rows


RuleFunc = Callable[[pd.Series], bool]


RULE_FIELDS = [
    "sample",
    "platform",
    "caller",
    "mode",
    "strategy",
    "rule_id",
    "trigger_count",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "baseline_f1",
    "truth_total",
    "fp_per_tp",
    "meets_safety",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-pool-tsv", required=True, help="Output from extract_borderline_rescue_candidates.py")
    parser.add_argument("--tp-summary-csv", required=True, help="InterSubMod TP significance_summary.csv")
    parser.add_argument("--fp-summary-csv", required=True, help="InterSubMod FP significance_summary.csv")
    parser.add_argument("--agreement-tsv", default="", help="Optional label_cluster_agreement.tsv")
    parser.add_argument("--baseline-tp", type=int, required=True, help="Baseline kept TP count")
    parser.add_argument("--baseline-fp", type=int, required=True, help="Baseline kept FP count")
    parser.add_argument("--truth-total", type=int, required=True, help="Truth total count")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def build_rules(mode: str) -> List[Dict[str, object]]:
    rules: List[Dict[str, object]] = []

    def add(strategy: str, rule_id: str, notes: str, func: RuleFunc) -> None:
        rules.append({"strategy": strategy, "rule_id": rule_id, "notes": notes, "func": func})

    # Caller-only baselines.
    add("caller_only", "caller_any_gq_ge_20", "no_candidate_gate;gq>=20", lambda row: row["gq"] >= 20)
    add(
        "caller_only",
        "caller_candidate_gq_ge_20",
        "candidate_eligible;gq>=20",
        lambda row: row["candidate_eligible"] and row["gq"] >= 20,
    )
    add(
        "caller_only",
        "caller_candidate_gq_ge_15",
        "candidate_eligible;gq>=15",
        lambda row: row["candidate_eligible"] and row["gq"] >= 15,
    )
    add(
        "caller_only",
        "caller_candidate_gq_ge_20_af_ge_010",
        "candidate_eligible;gq>=20;af>=0.10",
        lambda row: row["candidate_eligible"] and row["gq"] >= 20 and row["af"] >= 0.10,
    )
    add(
        "caller_only",
        "caller_candidate_qual_ge_10_or_gq_ge_20",
        "candidate_eligible;qual>=10 or gq>=20",
        lambda row: row["candidate_eligible"] and (row["qual"] >= 10.0 or row["gq"] >= 20),
    )
    if mode.startswith("to"):
        add(
            "caller_only",
            "caller_candidate_verdict_support_gq_ge_10",
            "candidate_eligible;verdict_somatic_or_subclonal;gq>=10",
            lambda row: row["candidate_eligible"] and row["gq"] >= 10 and (row["verdict_somatic"] or row["verdict_subclonal"]),
        )

    # Caller + methylation veto.
    add(
        "caller_plus_methylation_veto",
        "caller_candidate_gq_ge_20__veto_lowvaf_highadelta_lowcv",
        "candidate_eligible;gq>=20;not(lowVAF/highAlleleDelta/lowCramersV)",
        lambda row: row["candidate_eligible"] and row["gq"] >= 20 and not artifact_low_vaf_high_adelta_low_cv(row),
    )
    add(
        "caller_plus_methylation_veto",
        "caller_candidate_gq_ge_15__veto_lowvaf_highadelta_lowcv",
        "candidate_eligible;gq>=15;not(lowVAF/highAlleleDelta/lowCramersV)",
        lambda row: row["candidate_eligible"] and row["gq"] >= 15 and not artifact_low_vaf_high_adelta_low_cv(row),
    )
    add(
        "caller_plus_methylation_veto",
        "caller_candidate_gq_ge_20__veto_combined_artifact",
        "candidate_eligible;gq>=20;combined_artifact_veto",
        lambda row: row["candidate_eligible"] and row["gq"] >= 20 and not artifact_combined(row),
    )
    add(
        "caller_plus_methylation_veto",
        "caller_candidate_gq_ge_15__veto_combined_artifact",
        "candidate_eligible;gq>=15;combined_artifact_veto",
        lambda row: row["candidate_eligible"] and row["gq"] >= 15 and not artifact_combined(row),
    )

    # Caller + methylation support.
    add(
        "caller_plus_methylation_support",
        "caller_candidate_gq_ge_10__support_strong_or_subclone",
        "candidate_eligible;gq>=10;VerificationClass in Strong/Subclone",
        lambda row: row["candidate_eligible"] and row["gq"] >= 10 and support_strong_or_subclone(row),
    )
    add(
        "caller_plus_methylation_support",
        "caller_candidate_gq_ge_10__support_pairwise_ge_020",
        "candidate_eligible;gq>=10;PairwiseMedianDist>=0.20",
        lambda row: row["candidate_eligible"] and row["gq"] >= 10 and row["PairwiseMedianDist"] >= 0.20,
    )
    add(
        "caller_plus_methylation_support",
        "caller_candidate_gq_ge_10__support_strong_pairwise_ge_020",
        "candidate_eligible;gq>=10;Strong/Subclone;PairwiseMedianDist>=0.20",
        lambda row: row["candidate_eligible"]
        and row["gq"] >= 10
        and support_strong_or_subclone(row)
        and row["PairwiseMedianDist"] >= 0.20,
    )
    add(
        "caller_plus_methylation_support",
        "caller_candidate_gq_ge_10__support_agreement_positive",
        "candidate_eligible;gq>=10;agreement=label_upgrade/consistent_strong/consistent_subclone",
        lambda row: row["candidate_eligible"] and row["gq"] >= 10 and support_agreement(row),
    )
    return rules


def artifact_low_vaf_high_adelta_low_cv(row: pd.Series) -> bool:
    return row["af"] < 0.24 and row["AlleleDelta"] > 0.25 and row["CramersV"] < 0.05


def artifact_cluster_strong_label_weak(row: pd.Series) -> bool:
    return row["cluster_class"] == "Strong" and row["label_class"] == "Weak"


def artifact_low_pairwise(row: pd.Series) -> bool:
    return row["PairwiseMedianDist"] < 0.12


def artifact_low_hp_assign(row: pd.Series) -> bool:
    return row["hp_assign_rate"] < 0.50


def artifact_combined(row: pd.Series) -> bool:
    if artifact_low_vaf_high_adelta_low_cv(row):
        return True
    if artifact_cluster_strong_label_weak(row):
        return True
    if artifact_low_pairwise(row) and row["AlleleDelta"] > 0.15:
        return True
    if artifact_low_hp_assign(row) and row["cluster_class"] == "Strong":
        return True
    return False


def support_strong_or_subclone(row: pd.Series) -> bool:
    return row["VerificationClass"] in {"Strong", "Subclone"}


def support_agreement(row: pd.Series) -> bool:
    return row["agreement_type"] in {"label_upgrade", "consistent_strong", "consistent_subclone"}


def load_joined(args: argparse.Namespace) -> pd.DataFrame:
    candidate_pool = pd.read_csv(args.candidate_pool_tsv, sep="\t")
    candidate_pool = candidate_pool[
        candidate_pool["downstream_status"].isin({"caller_lost_tp", "caller_removed_fp"})
    ].copy()
    tp_summary = pd.read_csv(args.tp_summary_csv)
    fp_summary = pd.read_csv(args.fp_summary_csv)

    for df, scope in ((tp_summary, "tp"), (fp_summary, "fp")):
        df["region_key"] = (
            df["Chr"].astype(str) + ":" + df["Pos"].astype(str) + ":" + df["Ref"].astype(str) + ":" + df["Alt"].astype(str)
        )
        df["source_scope"] = scope

    sig = pd.concat([tp_summary, fp_summary], ignore_index=True)
    sig_columns = [
        "region_key",
        "source_scope",
        "PairwiseMeanDist",
        "PairwiseMedianDist",
        "AlleleDelta",
        "CramersV",
        "VerificationClass",
        "DominantLabel",
        "PassedGating",
        "GlobalP",
        "Quality_Score",
    ]
    joined = candidate_pool.merge(sig[sig_columns], on="region_key", how="left")

    if args.agreement_tsv:
        agreement = pd.read_csv(args.agreement_tsv, sep="\t")
        agreement_columns = [
            "region_key",
            "agreement_type",
            "class_shift",
            "hp_assign_rate",
            "allele_assign_rate",
            "cluster_class",
            "label_class",
        ]
        joined = joined.merge(agreement[agreement_columns], on="region_key", how="left")
    else:
        joined["agreement_type"] = ""
        joined["class_shift"] = ""
        joined["hp_assign_rate"] = 0.0
        joined["allele_assign_rate"] = 0.0
        joined["cluster_class"] = ""
        joined["label_class"] = ""

    for column in [
        "qual",
        "gq",
        "dp",
        "af",
        "ad_ref",
        "ad_alt",
        "PairwiseMeanDist",
        "PairwiseMedianDist",
        "AlleleDelta",
        "CramersV",
        "GlobalP",
        "Quality_Score",
        "hp_assign_rate",
        "allele_assign_rate",
    ]:
        joined[column] = pd.to_numeric(joined[column], errors="coerce")

    bool_columns = [
        "candidate_eligible",
        "rescue_filter_eligible",
        "meets_numeric_thresholds",
        "has_h_flag",
        "verdict_somatic",
        "verdict_subclonal",
        "verdict_germline",
        "pon_hit",
        "PassedGating",
    ]
    for column in bool_columns:
        if column in joined.columns:
            joined[column] = joined[column].map(lambda value: str(value).strip().lower() in {"true", "1", "yes"})

    fill_map = {
        "VerificationClass": "NA",
        "DominantLabel": "none",
        "agreement_type": "",
        "class_shift": "",
        "cluster_class": "",
        "label_class": "",
    }
    for column, value in fill_map.items():
        if column in joined.columns:
            joined[column] = joined[column].fillna(value)

    numeric_fill_zero = [
        "qual",
        "gq",
        "dp",
        "af",
        "ad_ref",
        "ad_alt",
        "PairwiseMeanDist",
        "PairwiseMedianDist",
        "AlleleDelta",
        "CramersV",
        "GlobalP",
        "Quality_Score",
        "hp_assign_rate",
        "allele_assign_rate",
    ]
    for column in numeric_fill_zero:
        joined[column] = joined[column].fillna(0.0)

    return joined


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    joined = load_joined(args)

    sample = str(joined["sample"].iloc[0]) if not joined.empty else "unknown"
    platform = str(joined["platform"].iloc[0]) if not joined.empty else infer_platform(sample)
    caller = str(joined["caller"].iloc[0]) if not joined.empty else "unknown"
    mode = str(joined["mode"].iloc[0]) if not joined.empty else "unknown"

    tp = joined[joined["downstream_status"] == "caller_lost_tp"].copy()
    fp = joined[joined["downstream_status"] == "caller_removed_fp"].copy()

    export_mask = joined["candidate_eligible"] | (joined["VerificationClass"] != "NA") | (joined["PairwiseMedianDist"] > 0)
    joined.loc[export_mask].to_csv(output_dir / "rescue_joined_features.tsv", sep="\t", index=False)

    baseline_metrics = compute_metrics(args.baseline_tp, args.baseline_fp, args.truth_total)
    rows: List[Dict[str, object]] = []

    for rule in build_rules(mode):
        tp_rescued = int(tp.apply(rule["func"], axis=1).sum())
        fp_reintroduced = int(fp.apply(rule["func"], axis=1).sum())
        metrics = compute_metrics(args.baseline_tp + tp_rescued, args.baseline_fp + fp_reintroduced, args.truth_total)
        delta_f1 = metrics["f1"] - baseline_metrics["f1"]
        fp_per_tp = fp_reintroduced / tp_rescued if tp_rescued else float("inf")
        rows.append(
            {
                "sample": sample,
                "platform": platform,
                "caller": caller,
                "mode": mode,
                "strategy": rule["strategy"],
                "rule_id": rule["rule_id"],
                "trigger_count": tp_rescued + fp_reintroduced,
                "tp_rescued": tp_rescued,
                "fp_reintroduced": fp_reintroduced,
                "precision": f"{metrics['precision']:.6f}",
                "recall": f"{metrics['recall']:.6f}",
                "f1": f"{metrics['f1']:.6f}",
                "delta_f1_vs_baseline": f"{delta_f1:.6f}",
                "baseline_f1": f"{baseline_metrics['f1']:.6f}",
                "truth_total": args.truth_total,
                "fp_per_tp": f"{fp_per_tp:.6f}" if tp_rescued else "inf",
                "meets_safety": tp_rescued > 0 and fp_reintroduced <= tp_rescued,
                "notes": rule["notes"],
            }
        )

    rows = sorted(
        rows,
        key=lambda row: (
            not bool(row["meets_safety"]),
            -int(row["tp_rescued"]),
            float(row["fp_per_tp"]) if row["fp_per_tp"] != "inf" else float("inf"),
            -float(row["delta_f1_vs_baseline"]),
        ),
    )
    write_tsv_rows(output_dir / "rescue_rule_comparison.tsv", RULE_FIELDS, rows)

    safe_rows = [row for row in rows if row["meets_safety"]]
    best_safe = safe_rows[0] if safe_rows else None
    md_lines = [
        f"# Rescue Rule Summary: {sample} ({mode})",
        "",
        f"- Baseline TP / FP / FN: `{args.baseline_tp}` / `{args.baseline_fp}` / `{args.truth_total - args.baseline_tp}`",
        f"- Baseline F1: `{baseline_metrics['f1']:.6f}`",
        f"- Downstream lost TP rows: `{len(tp)}`",
        f"- Downstream removed FP rows: `{len(fp)}`",
        "",
    ]
    if best_safe:
        md_lines.extend(
            [
                "## Best safe rule",
                "",
                f"- Strategy: `{best_safe['strategy']}`",
                f"- Rule: `{best_safe['rule_id']}`",
                f"- TP rescued / FP reintroduced: `{best_safe['tp_rescued']}` / `{best_safe['fp_reintroduced']}`",
                f"- F1 delta: `{best_safe['delta_f1_vs_baseline']}`",
                "",
            ]
        )
    md_lines.extend(
        [
            "## Top rules",
            "",
            "| strategy | rule_id | tp_rescued | fp_reintroduced | fp_per_tp | f1 | delta_f1_vs_baseline | meets_safety |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in rows[:12]:
        md_lines.append(
            f"| {row['strategy']} | {row['rule_id']} | {row['tp_rescued']} | {row['fp_reintroduced']} | {row['fp_per_tp']} | {row['f1']} | {row['delta_f1_vs_baseline']} | {row['meets_safety']} |"
        )
    (output_dir / "rescue_rule_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[evaluate_rescue_with_methylation] Wrote {output_dir / 'rescue_rule_comparison.tsv'}")


if __name__ == "__main__":
    main()
