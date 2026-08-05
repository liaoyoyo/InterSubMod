#!/usr/bin/env python3
"""Select 1:1 latest LongPhase-S PASS TP controls matched to each latest FP."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pysam
from scipy.spatial import cKDTree
from scipy.stats import ks_2samp


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_value(value: Any, default: float = 0.0) -> float:
    if value is None:
        return default
    if isinstance(value, tuple):
        value = value[0] if value else default
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def variant_row(record: pysam.VariantRecord) -> dict[str, Any]:
    call = next(iter(record.samples.values())) if record.samples else None
    af = first_value(call.get("AF") if call else None)
    dp = first_value(call.get("DP") if call else None)
    gq = first_value(call.get("GQ") if call else None)
    return {
        "key": (record.chrom, int(record.pos), record.ref, record.alts[0]),
        "chrom": record.chrom,
        "pos": int(record.pos),
        "ref": record.ref,
        "alt": record.alts[0],
        "af": af,
        "dp": dp,
        "gq": gq,
        "transformed": np.asarray(
            [
                math.log(min(max(af, 1e-4), 1 - 1e-4) / (1 - min(max(af, 1e-4), 1 - 1e-4))),
                math.log1p(max(dp, 0.0)),
                math.log1p(max(gq, 0.0)),
            ],
            dtype=float,
        ),
    }


def load_rows(path: Path) -> list[dict[str, Any]]:
    with pysam.VariantFile(str(path)) as variants:
        return [variant_row(record) for record in variants]


def robust_scale(fp_rows: list[dict[str, Any]], tp_rows: list[dict[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    fp = np.vstack([row["transformed"] for row in fp_rows])
    center = np.median(fp, axis=0)
    q75, q25 = np.percentile(fp, [75, 25], axis=0)
    scale = q75 - q25
    combined_sd = np.std(np.vstack([fp, np.vstack([row["transformed"] for row in tp_rows])]), axis=0)
    scale = np.where(scale > 1e-6, scale, np.where(combined_sd > 1e-6, combined_sd, 1.0))
    return center, scale


def select_controls(
    fp_rows: list[dict[str, Any]], tp_rows: list[dict[str, Any]], seed: int
) -> tuple[list[dict[str, Any]], list[float]]:
    if len(tp_rows) < len(fp_rows):
        raise RuntimeError(f"Not enough TP controls: FP={len(fp_rows)} TP={len(tp_rows)}")
    center, scale = robust_scale(fp_rows, tp_rows)
    fp_matrix = np.vstack([(row["transformed"] - center) / scale for row in fp_rows])
    tp_matrix = np.vstack([(row["transformed"] - center) / scale for row in tp_rows])
    tree = cKDTree(tp_matrix)
    k = min(512, len(tp_rows))
    distances, indices = tree.query(fp_matrix, k=k)
    if k == 1:
        distances = distances[:, None]
        indices = indices[:, None]
    nearest = distances[:, 0]
    rng = np.random.default_rng(seed)
    order = np.arange(len(fp_rows))
    rng.shuffle(order)
    order = sorted(order.tolist(), key=lambda index: nearest[index], reverse=True)
    used: set[int] = set()
    selected: list[dict[str, Any] | None] = [None] * len(fp_rows)
    selected_distance: list[float] = [math.nan] * len(fp_rows)
    for fp_index in order:
        chosen = None
        for distance, tp_index in zip(distances[fp_index], indices[fp_index]):
            candidate = int(tp_index)
            if candidate not in used:
                chosen = candidate
                selected_distance[fp_index] = float(distance)
                break
        if chosen is None:
            remaining = np.asarray(sorted(set(range(len(tp_rows))) - used), dtype=int)
            if not remaining.size:
                raise RuntimeError("TP control assignment exhausted")
            delta = tp_matrix[remaining] - fp_matrix[fp_index]
            local = int(np.argmin(np.sum(delta * delta, axis=1)))
            chosen = int(remaining[local])
            selected_distance[fp_index] = float(np.linalg.norm(delta[local]))
        used.add(chosen)
        selected[fp_index] = tp_rows[chosen]
    return [row for row in selected if row is not None], selected_distance


def standardized_mean_difference(first: np.ndarray, second: np.ndarray) -> float:
    denominator = math.sqrt((float(np.var(first, ddof=1)) + float(np.var(second, ddof=1))) / 2)
    return (float(np.mean(first)) - float(np.mean(second))) / denominator if denominator > 0 else 0.0


def balance_metrics(fp_rows: list[dict[str, Any]], controls: list[dict[str, Any]], distances: list[float]) -> dict[str, Any]:
    output: dict[str, Any] = {
        "nearest_distance_median": float(np.median(distances)),
        "nearest_distance_p95": float(np.percentile(distances, 95)),
        "nearest_distance_max": float(np.max(distances)),
    }
    for field in ("af", "dp", "gq"):
        first = np.asarray([row[field] for row in fp_rows], dtype=float)
        second = np.asarray([row[field] for row in controls], dtype=float)
        output[field] = {
            "fp_mean": float(first.mean()),
            "control_mean": float(second.mean()),
            "fp_median": float(np.median(first)),
            "control_median": float(np.median(second)),
            "smd": standardized_mean_difference(first, second),
            "ks_statistic": float(ks_2samp(first, second).statistic),
            "ks_p": float(ks_2samp(first, second).pvalue),
        }
    return output


def write_selected_vcf(source_path: Path, selected_keys: set[tuple[str, int, str, str]], output_path: Path) -> int:
    count = 0
    with pysam.VariantFile(str(source_path)) as source, pysam.VariantFile(
        str(output_path), "wz", header=source.header
    ) as output:
        for record in source:
            key = (record.chrom, int(record.pos), record.ref, record.alts[0])
            if key in selected_keys:
                output.write(record)
                count += 1
    pysam.tabix_index(str(output_path), preset="vcf", force=True)
    return count


def merge_windows(rows: list[dict[str, Any]], flank: int = 5000) -> list[tuple[str, int, int]]:
    order = {f"chr{index}": index for index in range(1, 23)}
    intervals = sorted(
        ((row["chrom"], max(0, row["pos"] - 1 - flank), row["pos"] + flank) for row in rows),
        key=lambda value: (order[value[0]], value[1], value[2]),
    )
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", type=Path, default=root / "results" / "latest_input_preflight.json")
    parser.add_argument("--output-root", type=Path, default=Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
        "20260715_single_fp_alt_multicluster_subclone_validation/matched_tp_control"
    ))
    parser.add_argument("--summary", type=Path, default=root / "results" / "matched_tp_control_preflight.json")
    args = parser.parse_args()

    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    args.output_root.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, Any]] = []
    sample_summaries: list[dict[str, Any]] = []
    for sample_index, source_entry in enumerate(preflight["samples"]):
        sample = source_entry["sample"]
        fp_path = Path(source_entry["latest_truth_fp"]["path"])
        tp_path = Path(source_entry["latest_truth_tp"]["path"])
        fp_rows = load_rows(fp_path)
        tp_rows = load_rows(tp_path)
        controls, distances = select_controls(fp_rows, tp_rows, seed=20260715 + sample_index * 101)
        selected_keys = {row["key"] for row in controls}
        if len(selected_keys) != len(fp_rows):
            raise RuntimeError(f"{sample} controls are not unique")
        sample_dir = args.output_root / "truth_matched" / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        output_vcf = sample_dir / f"{sample}.latest_lps_pass.truth_tp.matched_to_fp.vcf.gz"
        written = write_selected_vcf(tp_path, selected_keys, output_vcf)
        if written != len(fp_rows):
            raise RuntimeError(f"{sample} selected VCF count mismatch: expected={len(fp_rows)} written={written}")
        windows = merge_windows(controls)
        window_path = sample_dir / f"{sample}.matched_tp.w5000.merged.bed"
        with window_path.open("w", encoding="ascii") as handle:
            for chrom, start, end in windows:
                handle.write(f"{chrom}\t{start}\t{end}\n")
        pair_path = sample_dir / f"{sample}.fp_to_matched_tp.tsv"
        with pair_path.open("w", encoding="ascii") as handle:
            handle.write("fp_chrom\tfp_pos\tfp_ref\tfp_alt\tfp_af\tfp_dp\tfp_gq\t"
                         "tp_chrom\ttp_pos\ttp_ref\ttp_alt\ttp_af\ttp_dp\ttp_gq\tmatch_distance\n")
            for fp, control, distance in zip(fp_rows, controls, distances):
                handle.write(
                    f"{fp['chrom']}\t{fp['pos']}\t{fp['ref']}\t{fp['alt']}\t{fp['af']}\t{fp['dp']}\t{fp['gq']}\t"
                    f"{control['chrom']}\t{control['pos']}\t{control['ref']}\t{control['alt']}\t"
                    f"{control['af']}\t{control['dp']}\t{control['gq']}\t{distance}\n"
                )
        entry = {
            "sample": sample,
            "biological_id": source_entry["biological_id"],
            "materialization_label": "latest_lps_matched_tp_w5000",
            "raw_alignment": source_entry["raw_alignment"],
            "latest_read_tag_sidecar": source_entry["latest_read_tag_sidecar"],
            "window_bed": {"path": str(window_path), "sha256": sha256(window_path)},
            "latest_truth_fp": {"path": str(output_vcf), "sha256": sha256(output_vcf)},
            "matched_tp_vcf": {"path": str(output_vcf), "sha256": sha256(output_vcf)},
            "source_latest_truth_tp": source_entry["latest_truth_tp"],
            "source_latest_truth_fp": source_entry["latest_truth_fp"],
            "matched_pair_table": {"path": str(pair_path), "sha256": sha256(pair_path)},
            "matched_count": written,
            "merged_window_count": len(windows),
            "merged_window_bp": sum(end - start for _, start, end in windows),
        }
        entries.append(entry)
        sample_summaries.append(
            {
                "sample": sample,
                "fp_count": len(fp_rows),
                "tp_pool_count": len(tp_rows),
                "matched_tp_count": written,
                "merged_window_count": len(windows),
                "merged_window_bp": entry["merged_window_bp"],
                "balance": balance_metrics(fp_rows, controls, distances),
            }
        )
        print(
            f"[{sample}] matched={written}/{len(tp_rows)} windows={len(windows)} "
            f"AF_SMD={sample_summaries[-1]['balance']['af']['smd']:.4f}",
            flush=True,
        )
    summary = {
        "schema_name": "intersubmod.latest_lps_pass_fp_matched_tp_control_preflight",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "workspace": str(args.output_root),
        "source_preflight": str(args.preflight),
        "matching": (
            "1:1 without replacement within dataset; greedy nearest neighbor on robust-scaled "
            "logit(AF), log1p(DP), log1p(GQ); deterministic seed"
        ),
        "control_scope": "sampled_matched_control_not_all_truth_TP",
        "samples": entries,
        "sample_summaries": sample_summaries,
        "totals": {
            "fp_targets": sum(row["fp_count"] for row in sample_summaries),
            "matched_tp": sum(row["matched_tp_count"] for row in sample_summaries),
            "merged_window_bp": sum(row["merged_window_bp"] for row in sample_summaries),
        },
        "pass": all(row["fp_count"] == row["matched_tp_count"] for row in sample_summaries),
    }
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(args.summary), "totals": summary["totals"], "pass": summary["pass"]}, indent=2))


if __name__ == "__main__":
    main()
