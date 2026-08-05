#!/usr/bin/env python3
"""Analyze unsupervised methylation groups within focal ALT-supporting FP reads."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pysam

import focal_alt_cluster_lib as F


DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]


def chromosome_sort_key(chrom: str) -> tuple[int, str]:
    value = chrom.removeprefix("chr")
    if value.isdigit():
        return int(value), ""
    aliases = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return aliases.get(value.upper(), 10_000), value


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def json_text(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def scalar(value: Any) -> Any:
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def sample_value(call: pysam.libcbcf.VariantRecordSample, key: str) -> Any:
    try:
        value = call.get(key)
    except (KeyError, ValueError):
        return None
    return scalar(value)


def max_homopolymer(sequence: str) -> int:
    best = current = 0
    previous = ""
    for base in sequence.upper():
        if base == previous:
            current += 1
        else:
            current = 1
            previous = base
        best = max(best, current)
    return best


def sequence_entropy(sequence: str) -> float:
    counts = Counter(sequence.upper())
    total = sum(counts.values())
    if not total:
        return 0.0
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def load_variants(path: Path, reference: Path) -> dict[tuple[str, int], dict[str, Any]]:
    variants: dict[tuple[str, int], dict[str, Any]] = {}
    with pysam.FastaFile(str(reference)) as fasta, pysam.VariantFile(str(path)) as source:
        for record in source:
            key = (record.chrom, int(record.pos))
            if key in variants:
                raise RuntimeError(f"Position-only InterSubMod output cannot represent duplicate VCF key: {key}")
            call = next(iter(record.samples.values())) if record.samples else None
            start = max(0, record.pos - 11)
            context = fasta.fetch(record.chrom, start, record.pos + 10).upper()
            variants[key] = {
                "chrom": record.chrom,
                "pos": int(record.pos),
                "ref": record.ref,
                "alt": record.alts[0] if record.alts else None,
                "qual": record.qual,
                "filter": ";".join(record.filter.keys()),
                "caller_gt": json_text(sample_value(call, "GT")) if call else None,
                "caller_gq": sample_value(call, "GQ") if call else None,
                "caller_dp": sample_value(call, "DP") if call else None,
                "caller_af": sample_value(call, "AF") if call else None,
                "caller_ad": json_text(sample_value(call, "AD")) if call else None,
                "normal_af": sample_value(call, "NAF") if call else None,
                "normal_dp": sample_value(call, "NDP") if call else None,
                "normal_ad": json_text(sample_value(call, "NAD")) if call else None,
                "context_21bp": context,
                "context_max_homopolymer": max_homopolymer(context),
                "context_entropy": sequence_entropy(context),
            }
    return variants


def load_latest_ledger_context(path: Path, target_keys: set[tuple[str, int]]) -> dict[tuple[str, int], dict[str, Any]]:
    context: dict[tuple[str, int], dict[str, Any]] = {}
    if not path.exists():
        return context
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], int(row["pos"]))
            if key not in target_keys:
                continue
            family_coverage = json.loads(row.get("family_coverage_json") or "{}")
            context[key] = {
                "caller_original_filter": row.get("caller_filter"),
                "longphase_filter_transition": row.get("longphase_filter_transition"),
                "ssnv_branch": row.get("ssnv_branch"),
                "component_id": row.get("component_id"),
                "component_size": int(row["component_size"]) if row.get("component_size") else None,
                "selected_group_id": row.get("selected_group_id"),
                "selected_group_size": int(row["selected_group_size"]) if row.get("selected_group_size") else None,
                "layered_family_coverage": json_text(family_coverage),
                "layered_alt_families": json_text(
                    sorted(
                        family
                        for family, counts in family_coverage.items()
                        if isinstance(counts, list) and len(counts) >= 2 and counts[1] > 0
                    )
                ),
            }
    return context


def parse_site_key(reads_path: Path) -> tuple[str, int]:
    site_name = reads_path.parents[2].name
    chrom, position = site_name.rsplit("_", 1)
    return chrom, int(position)


def scan_site_tasks(
    sample: str,
    output_root: Path,
    variants: dict[tuple[str, int], dict[str, Any]],
    ledger: dict[tuple[str, int], dict[str, Any]],
    canonical_keys: set[tuple[str, int]],
) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    observed: set[tuple[str, int]] = set()
    for reads_path in sorted(output_root.rglob("reads.tsv")):
        key = parse_site_key(reads_path)
        if key in observed:
            raise RuntimeError(f"Duplicate InterSubMod site output: {sample} {key}")
        observed.add(key)
        if key not in variants:
            raise RuntimeError(f"InterSubMod output not present in source VCF: {sample} {key}")
        region_dir = reads_path.parents[1]
        metadata = dict(variants[key])
        metadata.update(ledger.get(key, {}))
        metadata["overlap_canonical_fp"] = key in canonical_keys
        tasks.append(
            {
                "sample": sample,
                "reads_path": str(reads_path),
                "distance_path": str(region_dir / "distance" / "BERNOULLI" / "matrix.csv"),
                "methylation_path": str(region_dir / "methylation" / "methylation.csv"),
                "region_dir": str(region_dir),
                "metadata": metadata,
            }
        )
    missing = set(variants) - observed
    extra = observed - set(variants)
    if missing or extra:
        raise RuntimeError(
            f"{sample} VCF/output site mismatch: vcf={len(variants)} output={len(observed)} "
            f"missing={len(missing)} extra={len(extra)}"
        )
    return tasks


def readset_digest(rows: list[dict[str, str]]) -> str:
    identities = sorted(
        f"{row['read_name']}|{row['chr']}|{row['start']}|{row['end']}|{row['strand']}" for row in rows
    )
    return hashlib.sha256("\n".join(identities).encode()).hexdigest()


def flatten_association(prefix: str, result: dict[str, Any], output: dict[str, Any]) -> None:
    for key, value in result.items():
        output[f"{prefix}_{key}"] = value


def process_site(task: dict[str, Any]) -> dict[str, Any]:
    sample = task["sample"]
    metadata = task["metadata"]
    chrom = metadata["chrom"]
    position = int(metadata["pos"])
    reads = F.load_reads(Path(task["reads_path"]))
    distance_ids, distance = F.load_matrix(Path(task["distance_path"]))
    methylation_ids, methylation = F.load_matrix(Path(task["methylation_path"]))
    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    alt_rows_raw = [
        row for row in reads.values() if F.is_tumor(row["is_tumor"]) and row["alt_support"] == "ALT"
    ]
    phase_anchored_tags = {"1-1", "2-1", "1-2", "2-2"}
    explicit_second_tags = {"1-2", "2-2"}
    ambiguous_tags = {"", ".", "0", "3", "4"}
    n_phase_anchored_raw = sum(row["hp"] in phase_anchored_tags for row in alt_rows_raw)
    n_explicit_second_raw = sum(row["hp"] in explicit_second_tags for row in alt_rows_raw)
    n_ambiguous_raw = sum(row["hp"] in ambiguous_tags for row in alt_rows_raw)
    alt_ids = [
        read_id
        for read_id in distance_ids
        if read_id in reads
        and read_id in methylation_index
        and F.is_tumor(reads[read_id]["is_tumor"])
        and reads[read_id]["alt_support"] == "ALT"
    ]
    output: dict[str, Any] = {
        "sample": sample,
        **metadata,
        "region_dir": task["region_dir"],
        "n_reads_total": len(reads),
        "n_alt_raw": len(alt_rows_raw),
        "n_alt_matrix": len(alt_ids),
        "alt_readset_sha256": readset_digest(alt_rows_raw),
        "alt_hp_counts": json_text(Counter(row["hp"] or "." for row in alt_rows_raw)),
        "alt_hp_family_counts": json_text(Counter(F.hp_family(row["hp"]) for row in alt_rows_raw)),
        "alt_strand_counts": json_text(Counter(row["strand"] for row in alt_rows_raw)),
        "phase_anchored_fraction_raw": n_phase_anchored_raw / len(alt_rows_raw) if alt_rows_raw else None,
        "explicit_second_haplotype_fraction_raw": (
            n_explicit_second_raw / len(alt_rows_raw) if alt_rows_raw else None
        ),
        "ambiguous_tag_fraction_raw": n_ambiguous_raw / len(alt_rows_raw) if alt_rows_raw else None,
        "n_alt_after_peel": 0,
        "n_alt_peeled": 0,
        "analysis_status": "insufficient_alt_reads",
        "forced_k": 0,
        "forced_silhouette": None,
        "coarse_ng": 0,
        "modal_fraction": None,
        "fine_ng": 0,
        "n_other": 0,
        "n_outlier": 0,
        "unstable": None,
        "hidden_heterogeneity": None,
        "stable_null_multigroup": False,
        "hp_axis_confound": False,
        "technical_axis_confound": False,
        "residual_unexplained_multigroup": False,
        "phase_anchored_robust_epigenetic_candidate": False,
        "orthogonal_subclone_confirmed": False,
        "evidence_tier": "E0_NOT_EVALUABLE",
        "cluster_sizes": "{}",
        "seed_group_counts": "[]",
        "cluster_detail": None,
    }
    if len(alt_ids) < 2 * F.MIN_SIZE:
        return output

    selected_distance_rows = [distance_index[read_id] for read_id in alt_ids]
    sub_distance = distance[np.ix_(selected_distance_rows, selected_distance_rows)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [alt_ids[index] for index in kept]
    output["n_alt_after_peel"] = len(kept_ids)
    output["n_alt_peeled"] = len(alt_ids) - len(kept_ids)
    if len(kept_ids) < 2 * F.MIN_SIZE:
        output["analysis_status"] = "incomplete_distance_below_minimum"
        return output

    sub_distance = sub_distance[np.ix_(kept, kept)]
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]
    forced = F.forced_silhouette_split(sub_distance.copy())
    output["forced_k"] = forced["k"]
    output["forced_silhouette"] = forced["silhouette"]
    phylo = F.analyze_phylo(sub_distance, sub_methylation)
    for key in (
        "coarse_ng",
        "modal_fraction",
        "fine_ng",
        "n_other",
        "n_outlier",
        "unstable",
        "hidden_heterogeneity",
        "modal_assignment_pair_count",
        "modal_assignment_ari_median",
        "modal_assignment_ari_min",
        "modal_assignment_all_pairs_ari_ge_0_8",
    ):
        output[key] = phylo[key]
    output["seed_group_counts"] = json_text(phylo["seed_group_counts"])
    output["analysis_status"] = "evaluable"
    output["evidence_tier"] = "E1_EVALUABLE_OR_SILHOUETTE_ONLY"
    output["stable_null_multigroup"] = phylo["coarse_ng"] >= 2 and not phylo["unstable"]

    labels = phylo["coarse_labels"]
    core = [index for index, label in enumerate(labels) if label not in {"other", "outlier"}]
    core_labels = [labels[index] for index in core]
    core_ids = [kept_ids[index] for index in core]
    output["cluster_sizes"] = json_text(Counter(core_labels))
    if output["stable_null_multigroup"]:
        output["evidence_tier"] = "E2_STABLE_NULL_MULTIGROUP"
        core_rows = [reads[read_id] for read_id in core_ids]
        cpg_called = [int(np.isfinite(sub_methylation[index]).sum()) for index in core]
        seed = F.stable_seed(sample, chrom, position)
        associations = {
            "hp_exact": F.categorical_permutation_association(
                [row["hp"] or "." for row in core_rows], core_labels, seed + 1
            ),
            "hp_family": F.categorical_permutation_association(
                [F.hp_family(row["hp"]) for row in core_rows], core_labels, seed + 2
            ),
            "strand": F.categorical_permutation_association(
                [row["strand"] for row in core_rows], core_labels, seed + 3
            ),
            "start": F.continuous_permutation_association(
                [float(row["start"]) for row in core_rows], core_labels, seed + 4
            ),
            "end": F.continuous_permutation_association(
                [float(row["end"]) for row in core_rows], core_labels, seed + 5
            ),
            "length": F.continuous_permutation_association(
                [float(row["end"]) - float(row["start"]) for row in core_rows], core_labels, seed + 6
            ),
            "mapq": F.continuous_permutation_association(
                [float(row["mapq"]) for row in core_rows], core_labels, seed + 7
            ),
            "cpg_called": F.continuous_permutation_association(cpg_called, core_labels, seed + 8),
        }
        for name, association in associations.items():
            flatten_association(name, association, output)
        output["hp_axis_confound"] = associations["hp_exact"]["aligned"] or associations["hp_family"]["aligned"]
        output["technical_axis_confound"] = any(
            associations[name]["aligned"]
            for name in ("strand", "start", "end", "length", "mapq", "cpg_called")
        )
        output["residual_unexplained_multigroup"] = not output["hp_axis_confound"] and not output[
            "technical_axis_confound"
        ]
        if output["residual_unexplained_multigroup"]:
            output["evidence_tier"] = "E3_UNEXPLAINED_AFTER_MEASURED_AXES"
        core_phase_anchored_fraction = (
            sum(row["hp"] in phase_anchored_tags for row in core_rows) / len(core_rows) if core_rows else 0.0
        )
        core_sizes = Counter(core_labels)
        output["phase_anchored_fraction_core"] = core_phase_anchored_fraction
        output["phase_anchored_robust_epigenetic_candidate"] = bool(
            output["residual_unexplained_multigroup"]
            and len(core_rows) >= 10
            and phylo["modal_fraction"] == 1.0
            and core_sizes
            and min(core_sizes.values()) >= 5
            and core_phase_anchored_fraction >= 0.80
        )
        if output["phase_anchored_robust_epigenetic_candidate"]:
            output["evidence_tier"] = "E4_PHASE_ANCHORED_ROBUST_EPIGENETIC_CANDIDATE"
        output["cluster_detail"] = {
            "sample": sample,
            "chrom": chrom,
            "pos": position,
            "region_dir": task["region_dir"],
            "read_ids": core_ids,
            "read_names": [reads[read_id]["read_name"] for read_id in core_ids],
            "labels": core_labels,
            "hp": [reads[read_id]["hp"] for read_id in core_ids],
            "strand": [reads[read_id]["strand"] for read_id in core_ids],
            "mapq": [int(reads[read_id]["mapq"]) for read_id in core_ids],
            "start": [int(reads[read_id]["start"]) for read_id in core_ids],
            "end": [int(reads[read_id]["end"]) for read_id in core_ids],
            "associations": associations,
            "coarse_split_trace": phylo["coarse_split_trace"],
            "fine_split_trace": phylo["fine_split_trace"],
            "coarse_labels_all_after_peel": labels,
            "all_after_peel_read_ids": kept_ids,
        }
    return output


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> list[float | None]:
    if total <= 0:
        return [None, None]
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    radius = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return [max(0.0, center - radius), min(1.0, center + radius)]


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    total = len(rows)
    evaluable = sum(row["analysis_status"] == "evaluable" for row in rows)
    forced = sum(row["analysis_status"] == "evaluable" and row["forced_k"] >= 2 for row in rows)
    stable = sum(bool(row["stable_null_multigroup"]) for row in rows)
    residual = sum(bool(row["residual_unexplained_multigroup"]) for row in rows)
    phase_anchored_robust = sum(bool(row["phase_anchored_robust_epigenetic_candidate"]) for row in rows)
    hp_confound = sum(bool(row["stable_null_multigroup"]) and bool(row["hp_axis_confound"]) for row in rows)
    technical = sum(
        bool(row["stable_null_multigroup"]) and bool(row["technical_axis_confound"]) for row in rows
    )
    readsets: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        readsets[(row["sample"], row["alt_readset_sha256"])].append(row)
    unique_evaluable = [group for group in readsets.values() if any(r["analysis_status"] == "evaluable" for r in group)]
    unique_stable = sum(any(bool(r["stable_null_multigroup"]) for r in group) for group in unique_evaluable)
    unique_residual = sum(any(bool(r["residual_unexplained_multigroup"]) for r in group) for group in unique_evaluable)
    return {
        "n_sites": total,
        "status_counts": dict(Counter(row["analysis_status"] for row in rows)),
        "n_evaluable_alt_ge6_after_peel": evaluable,
        "evaluable_fraction_all": evaluable / total if total else None,
        "n_forced_silhouette_multigroup": forced,
        "forced_fraction_evaluable": forced / evaluable if evaluable else None,
        "forced_fraction_evaluable_wilson95": wilson(forced, evaluable),
        "n_stable_null_multigroup": stable,
        "stable_null_fraction_evaluable": stable / evaluable if evaluable else None,
        "stable_null_fraction_evaluable_wilson95": wilson(stable, evaluable),
        "n_hp_axis_confound": hp_confound,
        "n_technical_axis_confound": technical,
        "n_residual_unexplained_multigroup": residual,
        "residual_fraction_evaluable": residual / evaluable if evaluable else None,
        "residual_fraction_evaluable_wilson95": wilson(residual, evaluable),
        "residual_fraction_stable_multigroup": residual / stable if stable else None,
        "n_phase_anchored_robust_epigenetic_candidate": phase_anchored_robust,
        "phase_anchored_robust_fraction_evaluable": (
            phase_anchored_robust / evaluable if evaluable else None
        ),
        "n_orthogonal_subclone_confirmed": sum(bool(row["orthogonal_subclone_confirmed"]) for row in rows),
        "evidence_tier_counts": dict(Counter(row["evidence_tier"] for row in rows)),
        "n_unique_alt_readsets": len(readsets),
        "n_unique_evaluable_alt_readsets": len(unique_evaluable),
        "n_unique_stable_null_multigroup_readsets": unique_stable,
        "unique_stable_fraction_evaluable": unique_stable / len(unique_evaluable) if unique_evaluable else None,
        "n_unique_residual_multigroup_readsets": unique_residual,
        "unique_residual_fraction_evaluable": unique_residual / len(unique_evaluable) if unique_evaluable else None,
        "n_sites_in_repeated_alt_readsets": sum(len(group) for group in readsets.values() if len(group) > 1),
        "n_repeated_alt_readset_groups": sum(len(group) > 1 for group in readsets.values()),
        "longphase_transition_counts": dict(Counter(row.get("longphase_filter_transition") or "NA" for row in rows)),
        "ssnv_branch_counts": dict(Counter(row.get("ssnv_branch") or "NA" for row in rows)),
    }


def output_rows(rows: list[dict[str, Any]], output_dir: Path, mode: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    details_path = output_dir / f"{mode}_stable_multigroup_read_assignments.jsonl"
    serializable_rows: list[dict[str, Any]] = []
    with details_path.open("w", encoding="utf-8") as details:
        for row in rows:
            row = dict(row)
            detail = row.pop("cluster_detail", None)
            if detail is not None:
                details.write(json.dumps(detail, ensure_ascii=False) + "\n")
            serializable_rows.append(row)
    fields = sorted({key for row in serializable_rows for key in row})
    site_path = output_dir / f"{mode}_site_results.tsv"
    with site_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(serializable_rows)

    readsets: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in serializable_rows:
        readsets[(row["sample"], row["alt_readset_sha256"])].append(row)
    repeated_path = output_dir / f"{mode}_repeated_alt_readsets.tsv"
    repeated_fields = [
        "sample",
        "alt_readset_sha256",
        "n_sites",
        "sites",
        "n_alt_raw",
        "any_evaluable",
        "any_stable_null_multigroup",
        "any_residual_unexplained_multigroup",
    ]
    with repeated_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, repeated_fields, delimiter="\t")
        writer.writeheader()
        for (sample, digest), group in sorted(readsets.items()):
            if len(group) <= 1:
                continue
            writer.writerow(
                {
                    "sample": sample,
                    "alt_readset_sha256": digest,
                    "n_sites": len(group),
                    "sites": ",".join(f"{row['chrom']}:{row['pos']}" for row in group),
                    "n_alt_raw": group[0]["n_alt_raw"],
                    "any_evaluable": any(row["analysis_status"] == "evaluable" for row in group),
                    "any_stable_null_multigroup": any(bool(row["stable_null_multigroup"]) for row in group),
                    "any_residual_unexplained_multigroup": any(
                        bool(row["residual_unexplained_multigroup"]) for row in group
                    ),
                }
            )


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", type=Path, default=root / "results" / "latest_input_preflight.json")
    parser.add_argument("--mode", choices=["latest", "canonical"], default="latest")
    parser.add_argument(
        "--result-label",
        default=None,
        help="Filename prefix for this population; defaults to --mode.",
    )
    parser.add_argument(
        "--population-label",
        default="truth-labeled FP sSNVs; focal tumor reads with alt_support=ALT",
    )
    parser.add_argument("--sample", action="append", choices=DATASETS)
    parser.add_argument("--workers", type=int, default=max(1, min(42, (os.cpu_count() or 4) - 4)))
    parser.add_argument("--reference", type=Path, default=Path("/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"))
    parser.add_argument("--latest-output-root", type=Path, default=None)
    parser.add_argument("--layered-root", type=Path, default=Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5"
    ))
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "focal_alt_multicluster")
    args = parser.parse_args()

    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    selected = set(args.sample or DATASETS)
    latest_output_root = args.latest_output_root or (Path(preflight["workspace"]) / "intersubmod_latest_fp_v1")
    tasks: list[dict[str, Any]] = []
    source_inventory: list[dict[str, Any]] = []
    for entry in preflight["samples"]:
        sample = entry["sample"]
        if sample not in selected:
            continue
        if args.mode == "latest":
            vcf_path = Path(entry["latest_truth_fp"]["path"])
            output_root = latest_output_root / sample
        else:
            vcf_path = Path(entry["canonical_fp_vcf"])
            output_root = Path(entry["canonical_fp_matrix_root"])
        variants = load_variants(vcf_path, args.reference)
        canonical_fp_vcf = entry.get("canonical_fp_vcf")
        canonical_variants = (
            load_variants(Path(canonical_fp_vcf), args.reference) if canonical_fp_vcf else {}
        )
        ledger = {}
        if args.mode == "latest":
            ledger_path = args.layered_root / "samples" / sample / f"ssnv_site_ledger_{sample}.tsv.gz"
            ledger = load_latest_ledger_context(ledger_path, set(variants))
        sample_tasks = scan_site_tasks(sample, output_root, variants, ledger, set(canonical_variants))
        tasks.extend(sample_tasks)
        source_inventory.append(
            {
                "sample": sample,
                "vcf": str(vcf_path),
                "output_root": str(output_root),
                "n_vcf_sites": len(variants),
                "n_output_sites": len(sample_tasks),
                "n_ledger_matched": len(ledger),
            }
        )
    if sum(item["n_vcf_sites"] for item in source_inventory) != len(tasks):
        raise RuntimeError("Aggregate task count does not equal aggregate VCF count")

    rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {executor.submit(process_site, task): task for task in tasks}
        completed = 0
        for future in as_completed(futures):
            task = futures[future]
            try:
                rows.append(future.result())
            except Exception as error:
                metadata = task["metadata"]
                failures.append(
                    {
                        "sample": task["sample"],
                        "site": f"{metadata['chrom']}:{metadata['pos']}",
                        "error": repr(error),
                    }
                )
            completed += 1
            if completed % 250 == 0 or completed == len(tasks):
                print(f"processed={completed}/{len(tasks)} failures={len(failures)}", flush=True)
    rows.sort(
        key=lambda row: (
            DATASETS.index(row["sample"]),
            chromosome_sort_key(row["chrom"]),
            row["pos"],
        )
    )
    if failures:
        raise RuntimeError(f"Focal ALT analysis had {len(failures)} site failures; first={failures[:3]}")

    result_label = args.result_label or args.mode
    output_rows(rows, args.output_dir, result_label)
    per_sample = {
        sample: summarize_rows([row for row in rows if row["sample"] == sample])
        for sample in DATASETS
        if sample in selected
    }
    summary = {
        "schema_name": "intersubmod.focal_alt_fp_multicluster",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "mode": args.mode,
        "population": args.population_label,
        "cluster_primary": "phylo-v4.1 null95 modal K=10, RNULL=40, min_group=3",
        "cluster_sensitivity": "historical forced silhouette best-k, min_group=3, no null gate",
        "source_inventory": source_inventory,
        "per_sample": per_sample,
        "pooled_site_weighted": summarize_rows(rows),
        "failures": failures,
        "pass": not failures and len(rows) == len(tasks),
        "interpretation_guardrail": (
            "A stable unexplained methylation multigroup is an epigenetic-heterogeneity candidate, "
            "not confirmation of a cellular subclone or linear ancestry."
        ),
    }
    summary_path = args.output_dir / f"{result_label}_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), "pooled": summary["pooled_site_weighted"], "pass": summary["pass"]}, indent=2))


if __name__ == "__main__":
    main()
