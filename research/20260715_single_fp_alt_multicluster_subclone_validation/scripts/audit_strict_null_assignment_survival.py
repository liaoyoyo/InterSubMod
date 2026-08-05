#!/usr/bin/env python3
"""Conservative survival audit for current ALT multigroup candidates."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

import focal_alt_cluster_lib as F


def flag(value: Any) -> bool:
    return str(value).lower() == "true"


def compact_result(result: dict[str, Any]) -> dict[str, Any]:
    root = result["coarse_split_trace"][0] if result["coarse_split_trace"] else {}
    stable = result["coarse_ng"] >= 2 and not result["unstable"]
    return {
        "stable": stable,
        "coarse_ng": result["coarse_ng"],
        "modal_fraction": result["modal_fraction"],
        "assignment_ari_median": result["modal_assignment_ari_median"],
        "assignment_ari_min": result["modal_assignment_ari_min"],
        "assignment_all_pairs_ari_ge_0_8": result["modal_assignment_all_pairs_ari_ge_0_8"],
        "root_empirical_p": root.get("empirical_p"),
        "root_null_threshold": root.get("null_threshold"),
        "root_n_valid_null": root.get("n_valid_null"),
        "root_failure": root.get("failure"),
    }


def process_task(task: dict[str, Any]) -> dict[str, Any]:
    region = Path(task["region_dir"])
    reads = F.load_reads(region / "reads/reads.tsv")
    distance_ids, distance = F.load_matrix(region / "distance/BERNOULLI/matrix.csv")
    methylation_ids, methylation = F.load_matrix(region / "methylation/methylation.csv")
    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    alt_ids = [
        read_id
        for read_id in distance_ids
        if read_id in reads
        and read_id in methylation_index
        and F.is_tumor(reads[read_id]["is_tumor"])
        and reads[read_id]["alt_support"] == "ALT"
    ]
    rows = [distance_index[read_id] for read_id in alt_ids]
    sub_distance = distance[np.ix_(rows, rows)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [alt_ids[index] for index in kept]
    if len(kept_ids) < 2 * F.MIN_SIZE:
        raise RuntimeError("Primary stable candidate became non-evaluable during strict audit")
    sub_distance = sub_distance[np.ix_(kept, kept)]
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]
    base_seed = F.stable_seed(task["sample"], task["chrom"], int(task["pos"])) + 500_000
    common = {
        "base_seed": base_seed,
        "seeds": task["seeds"],
        "rnull": task["rnull"],
        "empirical_alpha": task["empirical_alpha"],
    }
    column = F.analyze_phylo(
        sub_distance,
        sub_methylation,
        null_mode="column",
        **common,
    )
    row_circular = F.analyze_phylo(
        sub_distance,
        sub_methylation,
        null_mode="row_circular",
        **common,
    )
    return {
        **{key: value for key, value in task.items() if key not in {"seeds", "rnull", "empirical_alpha"}},
        "n_alt_after_peel": len(kept_ids),
        "column": compact_result(column),
        "row_circular": compact_result(row_circular),
    }


def summarize(rows: list[dict[str, Any]], n_evaluable: int, n_all: int) -> dict[str, Any]:
    primary_e4 = [row for row in rows if row["primary_e4"]]
    background = [row for row in primary_e4 if row["background_surviving"]]

    def survives(row: dict[str, Any], mode: str, require_assignment: bool = False) -> bool:
        result = row[mode]
        return result["stable"] and (
            not require_assignment or result["assignment_all_pairs_ari_ge_0_8"]
        )

    column = [row for row in rows if survives(row, "column")]
    row_circular = [row for row in rows if survives(row, "row_circular")]
    both_assignment = [
        row
        for row in primary_e4
        if row["background_surviving"]
        and survives(row, "column", require_assignment=True)
        and survives(row, "row_circular", require_assignment=True)
    ]
    return {
        "n_primary_group_count_stable_candidates": len(rows),
        "n_primary_e4_hp_covered_candidates": len(primary_e4),
        "n_e4_background_surviving": len(background),
        "n_column_empirical_p_gated_survivors": len(column),
        "n_row_circular_empirical_p_gated_survivors": len(row_circular),
        "n_e4_column_survivors": sum(survives(row, "column") for row in primary_e4),
        "n_e4_row_circular_survivors": sum(survives(row, "row_circular") for row in primary_e4),
        "n_e4_column_assignment_robust_survivors": sum(
            survives(row, "column", require_assignment=True) for row in primary_e4
        ),
        "n_e4_row_circular_assignment_robust_survivors": sum(
            survives(row, "row_circular", require_assignment=True) for row in primary_e4
        ),
        "n_strict_epigenetic_followup_candidates": len(both_assignment),
        "strict_fraction_evaluable": len(both_assignment) / n_evaluable if n_evaluable else None,
        "strict_fraction_all_sites": len(both_assignment) / n_all if n_all else None,
        "strict_sites": [f"{row['sample']}:{row['chrom']}:{row['pos']}" for row in both_assignment],
    }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-results", type=Path, required=True)
    parser.add_argument(
        "--ref-background",
        type=Path,
        default=root / "results/focal_alt_multicluster/ref_background_v2/ref_background_site_results.tsv",
    )
    parser.add_argument("--workers", type=int, default=max(1, min(42, (os.cpu_count() or 4) - 4)))
    parser.add_argument("--seeds", type=int, default=10)
    parser.add_argument("--rnull", type=int, default=199)
    parser.add_argument("--empirical-alpha", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    with args.site_results.open(encoding="utf-8") as handle:
        site_rows = list(csv.DictReader(handle, delimiter="\t"))
    with args.ref_background.open(encoding="utf-8") as handle:
        ref_rows = {
            (row["sample"], row["chrom"], int(row["pos"])): row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    n_evaluable = sum(row["analysis_status"] == "evaluable" for row in site_rows)
    tasks = []
    for row in site_rows:
        if not flag(row["stable_null_multigroup"]):
            continue
        key = (row["sample"], row["chrom"], int(row["pos"]))
        ref = ref_rows.get(key, {})
        background_surviving = (
            ref.get("ref_status") == "evaluable"
            and not flag(ref.get("ref_stable_null_multigroup"))
            and flag(ref.get("joint_stable_null_multigroup"))
            and flag(ref.get("joint_allele_axis_aligned"))
        )
        tasks.append(
            {
                "sample": row["sample"],
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "region_dir": row["region_dir"],
                "primary_e4": flag(row["phase_anchored_robust_epigenetic_candidate"]),
                "background_surviving": background_surviving,
                "seeds": args.seeds,
                "rnull": args.rnull,
                "empirical_alpha": args.empirical_alpha,
            }
        )

    rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {executor.submit(process_task, task): task for task in tasks}
        for index, future in enumerate(as_completed(futures), 1):
            task = futures[future]
            try:
                rows.append(future.result())
            except Exception as error:
                failures.append({"site": f"{task['sample']}:{task['chrom']}:{task['pos']}", "error": repr(error)})
            if index % 50 == 0 or index == len(tasks):
                print(f"processed={index}/{len(tasks)} failures={len(failures)}", flush=True)
    if failures:
        raise RuntimeError(f"Strict survival failures={len(failures)} first={failures[:3]}")
    rows.sort(key=lambda row: (row["sample"], row["chrom"], row["pos"]))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / "strict_null_assignment_site_results.tsv"
    flattened = []
    for row in rows:
        flat = {key: value for key, value in row.items() if key not in {"column", "row_circular"}}
        for mode in ("column", "row_circular"):
            flat.update({f"{mode}_{key}": value for key, value in row[mode].items()})
        flattened.append(flat)
    fields = sorted({key for row in flattened for key in row})
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(flattened)

    pooled = summarize(rows, n_evaluable, len(site_rows))
    summary = {
        "schema_name": "intersubmod.strict_null_assignment_survival",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {"site_results": str(args.site_results), "ref_background": str(args.ref_background)},
        "parameters": {
            "seeds": args.seeds,
            "rnull": args.rnull,
            "empirical_alpha": args.empirical_alpha,
            "minimum_valid_null_fraction": 0.8,
            "assignment_gate": "minimum pairwise ARI among modal-K seeds >= 0.8",
            "row_circular_null": "independent per-read cyclic shift preserving row values, missingness, and covariance",
        },
        "pooled": pooled,
        "per_sample": {
            sample: summarize(
                [row for row in rows if row["sample"] == sample],
                sum(
                    row["analysis_status"] == "evaluable" and row["sample"] == sample
                    for row in site_rows
                ),
                sum(row["sample"] == sample for row in site_rows),
            )
            for sample in sorted({row["sample"] for row in rows})
        },
        "failures": failures,
        "pass": not failures,
        "guardrail": (
            "This is a conservative survival sensitivity among primary candidates, not a calibrated genome-wide "
            "discovery procedure. It does not provide FDR, PPV, cellular-clone identity, or lineage topology."
        ),
    }
    summary_path = args.output_dir / "strict_null_assignment_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), "table": str(table_path), "pooled": pooled, "pass": summary["pass"]}, indent=2))


if __name__ == "__main__":
    main()
