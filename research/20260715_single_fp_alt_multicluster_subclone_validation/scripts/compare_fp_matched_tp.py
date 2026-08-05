#!/usr/bin/env python3
"""Compare focal-ALT methylation outcomes in truth FP sites and matched TP controls."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scipy.stats import binomtest


ENDPOINTS = {
    "evaluable": "analysis_status",
    "forced_silhouette_multigroup": "forced_k",
    "stable_null_multigroup": "stable_null_multigroup",
    "residual_unexplained_multigroup": "residual_unexplained_multigroup",
    "phase_anchored_robust_epigenetic_candidate": "phase_anchored_robust_epigenetic_candidate",
}


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def endpoint(row: dict[str, str], name: str) -> bool:
    if name == "evaluable":
        return row.get("analysis_status") == "evaluable"
    if name == "forced_silhouette_multigroup":
        value = row.get("forced_k", "")
        return value not in {"", "0", "0.0", "None", "NA"}
    return parse_bool(row.get(ENDPOINTS[name], ""))


def variant_key(sample: str, chrom: str, pos: Any, ref: str, alt: str) -> tuple[str, str, int, str, str]:
    return sample, chrom, int(pos), ref, alt


def load_site_rows(path: Path) -> dict[tuple[str, str, int, str, str], dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    indexed = {
        variant_key(row["sample"], row["chrom"], row["pos"], row["ref"], row["alt"]): row
        for row in rows
    }
    if len(indexed) != len(rows):
        raise RuntimeError(f"Duplicate site keys in {path}: rows={len(rows)} unique={len(indexed)}")
    return indexed


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> list[float | None]:
    if total == 0:
        return [None, None]
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return [max(0.0, center - half), min(1.0, center + half)]


def paired_endpoint(rows: list[dict[str, Any]], name: str) -> dict[str, Any]:
    fp = [endpoint(row["fp"], name) for row in rows]
    tp = [endpoint(row["tp"], name) for row in rows]
    n = len(rows)
    fp_yes = sum(fp)
    tp_yes = sum(tp)
    fp_only = sum(a and not b for a, b in zip(fp, tp))
    tp_only = sum(b and not a for a, b in zip(fp, tp))
    discordant = fp_only + tp_only
    p_value = binomtest(fp_only, discordant, 0.5, alternative="two-sided").pvalue if discordant else 1.0
    return {
        "n_pairs": n,
        "fp_n": fp_yes,
        "fp_fraction": fp_yes / n if n else None,
        "fp_wilson95": wilson(fp_yes, n),
        "tp_n": tp_yes,
        "tp_fraction": tp_yes / n if n else None,
        "tp_wilson95": wilson(tp_yes, n),
        "paired_risk_difference_fp_minus_tp": (fp_yes - tp_yes) / n if n else None,
        "fp_only_discordant": fp_only,
        "tp_only_discordant": tp_only,
        "exact_mcnemar_p": p_value,
    }


def summarize_subset(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {name: paired_endpoint(rows, name) for name in ENDPOINTS}


def unique_readset_summary(rows: list[dict[str, Any]], side: str) -> dict[str, Any]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for pair in rows:
        row = pair[side]
        groups[(row["sample"], row["alt_readset_sha256"])].append(row)
    evaluable = [group for group in groups.values() if any(endpoint(row, "evaluable") for row in group)]
    result: dict[str, Any] = {
        "n_unique_alt_readsets": len(groups),
        "n_unique_evaluable_alt_readsets": len(evaluable),
        "n_sites_in_repeated_alt_readsets": sum(len(group) for group in groups.values() if len(group) > 1),
    }
    for name in ENDPOINTS:
        success = sum(any(endpoint(row, name) for row in group) for group in groups.values())
        result[f"{name}_n"] = success
        result[f"{name}_fraction_all_unique_readsets"] = success / len(groups) if groups else None
    for name in ENDPOINTS:
        success = sum(any(endpoint(row, name) for row in group) for group in evaluable)
        result[f"{name}_n_among_evaluable"] = success
        result[f"{name}_fraction_evaluable_unique_readsets"] = success / len(evaluable) if evaluable else None
    return result


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--fp-sites", type=Path, required=True)
    parser.add_argument("--tp-sites", type=Path, required=True)
    parser.add_argument(
        "--matched-preflight",
        type=Path,
        default=root / "results" / "matched_tp_control_preflight.json",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--caliper", type=float, default=1.0)
    args = parser.parse_args()

    fp_sites = load_site_rows(args.fp_sites)
    tp_sites = load_site_rows(args.tp_sites)
    preflight = json.loads(args.matched_preflight.read_text(encoding="utf-8"))
    pairs: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for entry in preflight["samples"]:
        sample = entry["sample"]
        with Path(entry["matched_pair_table"]["path"]).open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                fp_key = variant_key(sample, row["fp_chrom"], row["fp_pos"], row["fp_ref"], row["fp_alt"])
                tp_key = variant_key(sample, row["tp_chrom"], row["tp_pos"], row["tp_ref"], row["tp_alt"])
                if fp_key not in fp_sites or tp_key not in tp_sites:
                    missing.append({"sample": sample, "fp_key": fp_key, "tp_key": tp_key})
                    continue
                pairs.append(
                    {
                        "sample": sample,
                        "match_distance": float(row["match_distance"]),
                        "fp": fp_sites[fp_key],
                        "tp": tp_sites[tp_key],
                    }
                )
    if missing:
        raise RuntimeError(f"Missing {len(missing)} matched site rows; first={missing[:3]}")
    expected = sum(int(entry["matched_count"]) for entry in preflight["samples"])
    if len(pairs) != expected:
        raise RuntimeError(f"Matched pair count mismatch: expected={expected}, observed={len(pairs)}")

    subsets = {
        "all_pairs": pairs,
        "caliper_le_1": [row for row in pairs if row["match_distance"] <= args.caliper],
        "both_evaluable": [
            row for row in pairs if endpoint(row["fp"], "evaluable") and endpoint(row["tp"], "evaluable")
        ],
        "caliper_le_1_both_evaluable": [
            row
            for row in pairs
            if row["match_distance"] <= args.caliper
            and endpoint(row["fp"], "evaluable")
            and endpoint(row["tp"], "evaluable")
        ],
    }
    per_sample: dict[str, Any] = {}
    for sample in sorted({row["sample"] for row in pairs}):
        sample_pairs = [row for row in pairs if row["sample"] == sample]
        per_sample[sample] = {
            "all_pairs": summarize_subset(sample_pairs),
            "caliper_le_1": summarize_subset(
                [row for row in sample_pairs if row["match_distance"] <= args.caliper]
            ),
            "both_evaluable": summarize_subset(
                [
                    row
                    for row in sample_pairs
                    if endpoint(row["fp"], "evaluable") and endpoint(row["tp"], "evaluable")
                ]
            ),
            "match_distance": {
                "n_le_caliper": sum(row["match_distance"] <= args.caliper for row in sample_pairs),
                "fraction_le_caliper": sum(row["match_distance"] <= args.caliper for row in sample_pairs)
                / len(sample_pairs),
            },
        }

    summary = {
        "schema_name": "intersubmod.focal_alt_fp_matched_tp_comparison",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "fp_sites": str(args.fp_sites),
            "tp_sites": str(args.tp_sites),
            "matched_preflight": str(args.matched_preflight),
        },
        "n_pairs": len(pairs),
        "caliper": args.caliper,
        "subset_counts": {name: len(rows) for name, rows in subsets.items()},
        "pooled_site_weighted": {name: summarize_subset(rows) for name, rows in subsets.items()},
        "per_sample": per_sample,
        "unique_readsets": {
            "fp": unique_readset_summary(pairs, "fp"),
            "tp": unique_readset_summary(pairs, "tp"),
        },
        "sample_counts": dict(Counter(row["sample"] for row in pairs)),
        "pass": len(pairs) == expected and not missing,
        "interpretation_guardrail": (
            "Matched TP controls test truth-label specificity. Residual AF/DP/GQ imbalance, especially HCC1395, "
            "precludes a causal truth-status interpretation; caliper and per-sample results are mandatory."
        ),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "fp_matched_tp_comparison_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    table_path = args.output_dir / "fp_matched_tp_per_sample.tsv"
    fields = [
        "sample",
        "subset",
        "endpoint",
        "n_pairs",
        "fp_n",
        "fp_fraction",
        "fp_wilson95",
        "tp_n",
        "tp_fraction",
        "tp_wilson95",
        "paired_risk_difference_fp_minus_tp",
        "fp_only_discordant",
        "tp_only_discordant",
        "exact_mcnemar_p",
    ]
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for sample, sample_result in per_sample.items():
            for subset in ("all_pairs", "caliper_le_1", "both_evaluable"):
                for name, values in sample_result[subset].items():
                    writer.writerow({"sample": sample, "subset": subset, "endpoint": name, **values})
    print(json.dumps({"summary": str(summary_path), "table": str(table_path), "pass": summary["pass"]}, indent=2))


if __name__ == "__main__":
    main()
