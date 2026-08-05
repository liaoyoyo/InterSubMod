#!/usr/bin/env python3
"""Test whether ALT-read methylation multigroups are also present in tumor REF reads."""

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


def true_value(value: Any) -> bool:
    return str(value).lower() == "true"


def load_detail(path: Path) -> dict[tuple[str, str, int], dict[str, Any]]:
    details: dict[tuple[str, str, int], dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            details[(row["sample"], row["chrom"], int(row["pos"]))] = row
    return details


def analyze_subset(
    ids: list[str],
    distance_ids: list[str],
    distance: np.ndarray,
    methylation_ids: list[str],
    methylation: np.ndarray,
    seed: int,
) -> dict[str, Any]:
    if len(ids) < 2 * F.MIN_SIZE:
        return {"status": "insufficient_reads", "n_raw": len(ids), "n_after_peel": 0}
    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    usable = [read_id for read_id in ids if read_id in distance_index and read_id in methylation_index]
    if len(usable) < 2 * F.MIN_SIZE:
        return {"status": "insufficient_matrix_reads", "n_raw": len(ids), "n_after_peel": 0}
    distance_rows = [distance_index[read_id] for read_id in usable]
    sub_distance = distance[np.ix_(distance_rows, distance_rows)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [usable[index] for index in kept]
    if len(kept_ids) < 2 * F.MIN_SIZE:
        return {
            "status": "incomplete_distance_below_minimum",
            "n_raw": len(ids),
            "n_after_peel": len(kept_ids),
        }
    sub_distance = sub_distance[np.ix_(kept, kept)]
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]
    result = F.analyze_phylo(sub_distance, sub_methylation, base_seed=seed)
    return {
        "status": "evaluable",
        "n_raw": len(ids),
        "n_after_peel": len(kept_ids),
        "kept_ids": kept_ids,
        "coarse_ng": result["coarse_ng"],
        "modal_fraction": result["modal_fraction"],
        "unstable": result["unstable"],
        "stable_null_multigroup": result["coarse_ng"] >= 2 and not result["unstable"],
        "labels": result["coarse_labels"],
        "cluster_sizes": dict(Counter(label for label in result["coarse_labels"] if label not in {"other", "outlier"})),
        "coarse_split_trace": result["coarse_split_trace"],
    }


def process_task(task: dict[str, Any]) -> dict[str, Any]:
    region = Path(task["region_dir"])
    reads = F.load_reads(region / "reads/reads.tsv")
    distance_ids, distance = F.load_matrix(region / "distance/BERNOULLI/matrix.csv")
    methylation_ids, methylation = F.load_matrix(region / "methylation/methylation.csv")
    tumor_alt = [
        read_id
        for read_id in distance_ids
        if read_id in reads and F.is_tumor(reads[read_id]["is_tumor"]) and reads[read_id]["alt_support"] == "ALT"
    ]
    tumor_ref = [
        read_id
        for read_id in distance_ids
        if read_id in reads and F.is_tumor(reads[read_id]["is_tumor"]) and reads[read_id]["alt_support"] == "REF"
    ]
    seed = F.stable_seed(task["sample"], task["chrom"], int(task["pos"]))
    ref = analyze_subset(tumor_ref, distance_ids, distance, methylation_ids, methylation, seed + 100_000)
    output = {
        **task,
        "n_tumor_alt": len(tumor_alt),
        "n_tumor_ref": len(tumor_ref),
        "ref_status": ref["status"],
        "ref_n_after_peel": ref.get("n_after_peel", 0),
        "ref_coarse_ng": ref.get("coarse_ng"),
        "ref_modal_fraction": ref.get("modal_fraction"),
        "ref_stable_null_multigroup": ref.get("stable_null_multigroup", False),
        "ref_cluster_sizes": json.dumps(ref.get("cluster_sizes", {}), sort_keys=True, separators=(",", ":")),
        "joint_status": "not_run_non_e4",
        "joint_stable_null_multigroup": False,
        "joint_allele_axis_aligned": False,
    }
    if not task["is_e4"]:
        return output
    if len(tumor_ref) < F.MIN_SIZE:
        output["joint_status"] = "insufficient_ref_for_joint_control"
        return output

    joint = analyze_subset(
        tumor_alt + tumor_ref,
        distance_ids,
        distance,
        methylation_ids,
        methylation,
        seed + 200_000,
    )
    output.update(
        {
            "joint_status": joint["status"],
            "joint_n_after_peel": joint.get("n_after_peel", 0),
            "joint_coarse_ng": joint.get("coarse_ng"),
            "joint_modal_fraction": joint.get("modal_fraction"),
            "joint_stable_null_multigroup": joint.get("stable_null_multigroup", False),
            "joint_cluster_sizes": json.dumps(
                joint.get("cluster_sizes", {}), sort_keys=True, separators=(",", ":")
            ),
        }
    )
    if joint["status"] == "evaluable" and joint.get("stable_null_multigroup"):
        allele_by_id = {read_id: "ALT" for read_id in tumor_alt}
        allele_by_id.update({read_id: "REF" for read_id in tumor_ref})
        core = [
            (read_id, label)
            for read_id, label in zip(joint["kept_ids"], joint["labels"])
            if label not in {"other", "outlier"}
        ]
        association = F.categorical_permutation_association(
            [allele_by_id[read_id] for read_id, _ in core],
            [label for _, label in core],
            seed + 300_000,
        )
        output.update(
            {
                "joint_allele_v": association["v"],
                "joint_allele_p_perm": association["p_perm"],
                "joint_allele_n": association["n"],
                "joint_allele_axis_aligned": association["aligned"],
            }
        )
    return output


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    ref_evaluable = [row for row in rows if row["ref_status"] == "evaluable"]
    e4 = [row for row in rows if row["is_e4"]]
    e4_ref_evaluable = [row for row in e4 if row["ref_status"] == "evaluable"]
    e4_joint_evaluable = [row for row in e4 if row["joint_status"] == "evaluable"]
    e4_joint_stable = [row for row in e4_joint_evaluable if row["joint_stable_null_multigroup"]]
    return {
        "n_alt_stable_sites": len(rows),
        "n_ref_evaluable": len(ref_evaluable),
        "n_ref_stable_null_multigroup": sum(row["ref_stable_null_multigroup"] for row in ref_evaluable),
        "ref_stable_fraction_evaluable": (
            sum(row["ref_stable_null_multigroup"] for row in ref_evaluable) / len(ref_evaluable)
            if ref_evaluable
            else None
        ),
        "n_e4": len(e4),
        "n_e4_with_no_ref_reads": sum(row["n_tumor_ref"] == 0 for row in e4),
        "n_e4_with_ref_below_joint_minimum": sum(row["n_tumor_ref"] < F.MIN_SIZE for row in e4),
        "n_e4_ref_evaluable": len(e4_ref_evaluable),
        "n_e4_ref_stable_null_multigroup": sum(
            row["ref_stable_null_multigroup"] for row in e4_ref_evaluable
        ),
        "n_e4_joint_evaluable": len(e4_joint_evaluable),
        "n_e4_joint_stable_null_multigroup": len(e4_joint_stable),
        "n_e4_joint_allele_axis_aligned": sum(row["joint_allele_axis_aligned"] for row in e4_joint_stable),
        "n_e4_joint_not_allele_aligned": sum(
            not row["joint_allele_axis_aligned"] for row in e4_joint_stable
        ),
    }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--site-results",
        type=Path,
        default=root / "results/focal_alt_multicluster/latest_full_v1/latest_site_results_with_topology.tsv",
    )
    parser.add_argument(
        "--details",
        type=Path,
        default=root
        / "results/focal_alt_multicluster/latest_full_v1/latest_stable_multigroup_read_assignments.jsonl",
    )
    parser.add_argument("--workers", type=int, default=max(1, min(42, (os.cpu_count() or 4) - 4)))
    parser.add_argument(
        "--output-dir", type=Path, default=root / "results/focal_alt_multicluster/ref_background_v1"
    )
    args = parser.parse_args()

    with args.site_results.open(encoding="utf-8") as handle:
        site_rows = list(csv.DictReader(handle, delimiter="\t"))
    details = load_detail(args.details)
    tasks = []
    for row in site_rows:
        if not true_value(row["stable_null_multigroup"]):
            continue
        detail_key = (row["sample"], row["chrom"], int(row["pos"]))
        if detail_key not in details:
            raise RuntimeError(f"Missing stable-site detail: {detail_key}")
        tasks.append(
            {
                "sample": row["sample"],
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "region_dir": row["region_dir"],
                "is_e4": true_value(row["phase_anchored_robust_epigenetic_candidate"]),
                "alt_evidence_tier": row["evidence_tier"],
                "layered_topology_context": row["layered_topology_context"],
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
        raise RuntimeError(f"REF-background analysis failures={len(failures)} first={failures[:3]}")
    rows.sort(key=lambda row: (row["sample"], row["chrom"], row["pos"]))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / "ref_background_site_results.tsv"
    fields = sorted({name for row in rows for name in row})
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "schema_name": "intersubmod.focal_alt_ref_background_controls",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "All 583 latest ALT-stable FP sites; independent tumor-REF null clustering; E4 joint ALT+REF clustering",
        "method": "Same phylo-v4.1 null95/modal K=10/RNULL=40; allele association only after stable joint clustering",
        "pooled": summarize(rows),
        "per_sample": {
            sample: summarize([row for row in rows if row["sample"] == sample])
            for sample in sorted({row["sample"] for row in rows})
        },
        "normal_background_status": "NOT_TESTED_latest_InterSubMod_run_has_no_normal_bam",
        "failures": failures,
        "pass": len(rows) == 583 and not failures,
        "guardrail": (
            "REF replication or absent allele alignment argues against ALT-specificity. Lack of REF replication does "
            "not prove a subclone because REF power, CN, mapping, normal background, and orthogonal lineage evidence remain."
        ),
    }
    summary_path = args.output_dir / "ref_background_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), "table": str(table_path), "pooled": summary["pooled"], "pass": summary["pass"]}, indent=2))


if __name__ == "__main__":
    main()
