#!/usr/bin/env python3
"""Materialize bounded BAMs with latest LongPhase-S HP/PS tags from exact sidecars."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pysam


DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def cigar_digest(alignment: pysam.AlignedSegment) -> str:
    return hashlib.blake2b((alignment.cigarstring or "*").encode(), digest_size=8).hexdigest()


def alignment_key(alignment: pysam.AlignedSegment) -> tuple[str, str, int, int, int, str]:
    return (
        alignment.query_name,
        alignment.reference_name,
        int(alignment.reference_start),
        int(alignment.reference_end or alignment.reference_start + 1),
        int(alignment.flag),
        cigar_digest(alignment),
    )


def alignment_digest(alignment: pysam.AlignedSegment) -> str:
    # HP/PS are replaced by the authoritative sidecar. RG is not part of the
    # coordinate_join_v1 identity and can differ between redundant merged inputs.
    immutable_payload = (
        alignment.query_sequence,
        bytes(alignment.query_qualities or []),
        int(alignment.mapping_quality),
        int(alignment.next_reference_id),
        int(alignment.next_reference_start),
        int(alignment.template_length),
        alignment.cigarstring,
        sorted((key, repr(value)) for key, value in alignment.get_tags() if key not in {"HP", "PS", "RG"}),
    )
    return hashlib.blake2b(repr(immutable_payload).encode(), digest_size=16).hexdigest()


def analysis_payload_digest(alignment: pysam.AlignedSegment, include_methylation: bool = True) -> str:
    """Hash fields consumed by InterSubMod, excluding replaceable HP/PS tags."""
    tags = dict(alignment.get_tags())
    methylation_payload = tuple(
        (tag, repr(tags.get(tag))) for tag in ("MM", "ML", "MN")
    ) if include_methylation else ()
    payload = (
        alignment.query_sequence,
        bytes(alignment.query_qualities or []),
        int(alignment.mapping_quality),
        int(alignment.next_reference_id),
        int(alignment.next_reference_start),
        int(alignment.template_length),
        alignment.cigarstring,
        methylation_payload,
    )
    return hashlib.blake2b(repr(payload).encode(), digest_size=16).hexdigest()


def has_complete_methylation_payload(alignment: pysam.AlignedSegment) -> bool:
    return alignment.has_tag("MM") and alignment.has_tag("ML")


def differing_tag_names(first: pysam.AlignedSegment, second: pysam.AlignedSegment) -> list[str]:
    first_tags = {key: repr(value) for key, value in first.get_tags()}
    second_tags = {key: repr(value) for key, value in second.get_tags()}
    return sorted(
        key for key in set(first_tags) | set(second_tags) if first_tags.get(key) != second_tags.get(key)
    )


def resolve_duplicate(
    first: pysam.AlignedSegment, second: pysam.AlignedSegment
) -> tuple[pysam.AlignedSegment, str, list[str]]:
    """Select one analytically equivalent record or fail on a material conflict."""
    differences = differing_tag_names(first, second)
    if alignment_digest(first) == alignment_digest(second):
        return first, "exact", differences
    if analysis_payload_digest(first) == analysis_payload_digest(second):
        return first, "nonanalytic_aux_tag_variation", differences
    if analysis_payload_digest(first, include_methylation=False) == analysis_payload_digest(
        second, include_methylation=False
    ):
        first_complete = has_complete_methylation_payload(first)
        second_complete = has_complete_methylation_payload(second)
        if first_complete != second_complete:
            return (
                first if first_complete else second,
                "prefer_complete_methylation_payload",
                differences,
            )
    raise RuntimeError(f"material duplicate conflict; differing tags={differences}")


def fingerprints(alignment: pysam.AlignedSegment) -> tuple[str, str, str, bool]:
    return (
        alignment_digest(alignment),
        analysis_payload_digest(alignment),
        analysis_payload_digest(alignment, include_methylation=False),
        has_complete_methylation_payload(alignment),
    )


def load_windows(path: Path) -> list[tuple[str, int, int]]:
    windows: list[tuple[str, int, int]] = []
    previous: tuple[str, int, int] | None = None
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start0, end0 = line.rstrip().split("\t")[:3]
            current = (chrom, int(start0), int(end0))
            if current[1] >= current[2]:
                raise ValueError(f"Invalid interval: {current}")
            if previous is not None and previous[0] == current[0] and current[1] <= previous[2]:
                raise ValueError(f"Intervals must be sorted and non-overlapping: {previous}, {current}")
            windows.append(current)
            previous = current
    return windows


def fetch_sidecar(
    tabix: pysam.TabixFile, chrom: str, start0: int, end0: int
) -> tuple[dict[tuple[str, str, int, int, int, str], tuple[str | None, str | None]], set[tuple[Any, ...]], int]:
    values: dict[tuple[str, str, int, int, int, str], tuple[str | None, str | None]] = {}
    conflicts: set[tuple[Any, ...]] = set()
    rows = 0
    try:
        lines = tabix.fetch(chrom, start0, end0)
    except (ValueError, OSError):
        return values, conflicts, rows
    for line in lines:
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"Malformed sidecar row with {len(fields)} fields")
        row_chrom, start, end, qname, flag, _mapq, cigar_b2, hp, ps = fields
        key = (qname, row_chrom, int(start), int(end), int(flag), cigar_b2)
        value = (None if hp == "." else hp, None if ps == "." else ps)
        if key in values and values[key] != value:
            conflicts.add(key)
        else:
            values[key] = value
        rows += 1
    return values, conflicts, rows


def set_latest_tags(alignment: pysam.AlignedSegment, hp: str | None, ps: str | None) -> None:
    if alignment.has_tag("HP"):
        alignment.set_tag("HP", None)
    if alignment.has_tag("PS"):
        alignment.set_tag("PS", None)
    if hp is not None:
        alignment.set_tag("HP", hp, value_type="Z")
    if ps is not None:
        try:
            ps_value = int(ps)
        except ValueError as error:
            raise ValueError(f"PS must be integer or '.', observed {ps!r}") from error
        alignment.set_tag("PS", ps_value, value_type="i")


def materialize_sample(
    sample_entry: dict[str, Any],
    workspace: str,
    compression_threads: int,
    force_new_suffix: str = "",
) -> dict[str, Any]:
    sample = sample_entry["sample"]
    output_dir = Path(workspace) / "latest_tagged_subset" / sample
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = f".{force_new_suffix}" if force_new_suffix else ""
    materialization_label = sample_entry.get("materialization_label", "latest_lps_fp_w5000")
    output_bam = output_dir / f"{sample}.{materialization_label}.tagged_subset{suffix}.bam"
    partial_bam = output_dir / f"{sample}.{materialization_label}.tagged_subset{suffix}.partial.bam"
    receipt_path = output_dir / f"{sample}.{materialization_label}.materialization_receipt{suffix}.json"
    if output_bam.exists() and receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("pass") is True:
            return receipt

    raw_bam = Path(sample_entry["raw_alignment"]["path"])
    raw_bai = Path(sample_entry["raw_alignment"]["index_path"])
    sidecar = Path(sample_entry["latest_read_tag_sidecar"]["path"])
    sidecar_index = Path(sample_entry["latest_read_tag_sidecar"]["index_path"])
    window_bed = Path(sample_entry["window_bed"]["path"])
    for path in (raw_bam, raw_bai, sidecar, sidecar_index, window_bed):
        if not path.exists():
            raise FileNotFoundError(path)

    windows = load_windows(window_bed)
    diagnostics: Counter[str] = Counter(
        {
            "windows": 0,
            "sidecar_rows_exposed": 0,
            "sidecar_conflicts": 0,
            "raw_alignment_exposures": 0,
            "duplicate_identity_exposures": 0,
            "duplicate_alignment_conflicts": 0,
            "duplicate_identity_collapsed": 0,
            "duplicate_exact_collapsed": 0,
            "duplicate_nonanalytic_aux_variation_collapsed": 0,
            "duplicate_methyl_payload_resolved_prefer_complete": 0,
            "sidecar_missing": 0,
            "written_unique_alignments": 0,
            "with_MM": 0,
            "with_ML": 0,
            "with_MM_ML": 0,
        }
    )
    hp_counts: Counter[str] = Counter()
    ps_counts: Counter[str] = Counter()
    duplicate_tag_difference_counts: Counter[str] = Counter()
    duplicate_resolution_examples: list[dict[str, Any]] = []
    first_signatures: dict[tuple[str, str, int, int, int, str], tuple[str, str, str, bool]] = {}
    rank: dict[str, int] = {}
    previous_order: tuple[int, int] | None = None
    sorted_ok = True
    started_at = now_utc()

    with pysam.AlignmentFile(str(raw_bam), "rb", threads=max(1, compression_threads)) as source:
        header = source.header.to_dict()
        header.setdefault("PG", []).append(
            {
                "ID": "InterSubMod_latest_sidecar_materializer",
                "PN": "materialize_latest_tagged_subset.py",
                "VN": "1.0.0",
                "CL": "bounded exact coordinate join; latest LongPhase-S HP/PS sidecar",
            }
        )
        rank = {name: index for index, name in enumerate(source.references)}
        with pysam.TabixFile(str(sidecar), index=str(sidecar_index)) as tags, pysam.AlignmentFile(
            str(partial_bam), "wb", header=header, threads=max(1, compression_threads)
        ) as output:
            for chrom, start0, end0 in windows:
                diagnostics["windows"] += 1
                sidecar_values, sidecar_conflicts, sidecar_rows = fetch_sidecar(tags, chrom, start0, end0)
                diagnostics["sidecar_rows_exposed"] += sidecar_rows
                diagnostics["sidecar_conflicts"] += len(sidecar_conflicts)
                if sidecar_conflicts:
                    raise RuntimeError(f"{sample} sidecar conflicts in {chrom}:{start0}-{end0}")
                window_records: dict[
                    tuple[str, str, int, int, int, str], pysam.AlignedSegment
                ] = {}
                for alignment in source.fetch(chrom, start0, end0):
                    if alignment.is_unmapped:
                        continue
                    diagnostics["raw_alignment_exposures"] += 1
                    key = alignment_key(alignment)
                    if key in window_records:
                        diagnostics["duplicate_identity_exposures"] += 1
                        try:
                            selected, reason, differences = resolve_duplicate(window_records[key], alignment)
                        except RuntimeError as error:
                            diagnostics["duplicate_alignment_conflicts"] += 1
                            raise RuntimeError(f"{sample} conflicting duplicate BAM identity: {key}; {error}") from error
                        window_records[key] = selected
                        diagnostics["duplicate_identity_collapsed"] += 1
                        diagnostics[
                            {
                                "exact": "duplicate_exact_collapsed",
                                "nonanalytic_aux_tag_variation": "duplicate_nonanalytic_aux_variation_collapsed",
                                "prefer_complete_methylation_payload": (
                                    "duplicate_methyl_payload_resolved_prefer_complete"
                                ),
                            }[reason]
                        ] += 1
                        duplicate_tag_difference_counts.update(differences)
                        if reason != "exact" and len(duplicate_resolution_examples) < 20:
                            duplicate_resolution_examples.append(
                                {"identity": list(key), "reason": reason, "differing_tags": differences}
                            )
                        continue
                    window_records[key] = alignment

                for key, alignment in window_records.items():
                    signature = fingerprints(alignment)
                    if key in first_signatures:
                        diagnostics["duplicate_identity_exposures"] += 1
                        previous = first_signatures[key]
                        if previous[0] == signature[0]:
                            reason = "exact"
                        elif previous[1] == signature[1]:
                            reason = "nonanalytic_aux_tag_variation"
                        elif previous[2] == signature[2] and previous[3] and not signature[3]:
                            reason = "prefer_complete_methylation_payload"
                        else:
                            diagnostics["duplicate_alignment_conflicts"] += 1
                            raise RuntimeError(f"{sample} conflicting duplicate BAM identity: {key}")
                        diagnostics["duplicate_identity_collapsed"] += 1
                        diagnostics[
                            {
                                "exact": "duplicate_exact_collapsed",
                                "nonanalytic_aux_tag_variation": "duplicate_nonanalytic_aux_variation_collapsed",
                                "prefer_complete_methylation_payload": (
                                    "duplicate_methyl_payload_resolved_prefer_complete"
                                ),
                            }[reason]
                        ] += 1
                        continue
                    first_signatures[key] = signature
                    if key not in sidecar_values:
                        diagnostics["sidecar_missing"] += 1
                        raise RuntimeError(f"{sample} missing latest sidecar tag identity: {key}")
                    hp, ps = sidecar_values[key]
                    set_latest_tags(alignment, hp, ps)
                    order = (rank.get(alignment.reference_name, 10**9), int(alignment.reference_start))
                    if previous_order is not None and order < previous_order:
                        sorted_ok = False
                    previous_order = order
                    output.write(alignment)
                    diagnostics["written_unique_alignments"] += 1
                    hp_counts[hp or "."] += 1
                    if ps is not None:
                        ps_counts[hp or "."] += 1
                    if alignment.has_tag("MM"):
                        diagnostics["with_MM"] += 1
                    if alignment.has_tag("ML"):
                        diagnostics["with_ML"] += 1
                    if alignment.has_tag("MM") and alignment.has_tag("ML"):
                        diagnostics["with_MM_ML"] += 1

    passed = (
        diagnostics["written_unique_alignments"] > 0
        and diagnostics["sidecar_missing"] == 0
        and diagnostics["sidecar_conflicts"] == 0
        and diagnostics["duplicate_alignment_conflicts"] == 0
        and sorted_ok
    )
    if not passed:
        raise RuntimeError(f"{sample} materialization validation failed: {dict(diagnostics)} sorted={sorted_ok}")
    if output_bam.exists():
        raise FileExistsError(f"Refusing to overwrite existing final BAM: {output_bam}")
    os.rename(partial_bam, output_bam)
    pysam.index("-@", str(max(1, compression_threads)), str(output_bam))
    output_bai = Path(str(output_bam) + ".bai")
    quickcheck = pysam.quickcheck(str(output_bam))
    passed = passed and quickcheck == "" and output_bai.exists()
    receipt = {
        "schema_name": "intersubmod.latest_tagged_subset_materialization",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "started_at_utc": started_at,
        "sample": sample,
        "materialization_label": materialization_label,
        "source_raw_bam": str(raw_bam),
        "source_raw_bai": str(raw_bai),
        "source_latest_sidecar": str(sidecar),
        "source_latest_sidecar_index": str(sidecar_index),
        "window_bed": str(window_bed),
        "output_bam": str(output_bam),
        "output_bai": str(output_bai),
        "output_size_bytes": output_bam.stat().st_size,
        "diagnostics": dict(diagnostics),
        "duplicate_tag_difference_counts": dict(duplicate_tag_difference_counts),
        "duplicate_resolution_examples": duplicate_resolution_examples,
        "HP_counts": dict(hp_counts),
        "HP_with_PS_counts": dict(ps_counts),
        "coordinate_sorted": sorted_ok,
        "quickcheck_empty_output": quickcheck == "",
        "duplicate_identity_policy": (
            "collapse_exact_or_InterSubMod-equivalent_records; prefer MM+ML-complete over missing; "
            "fail_on_sequence_quality_mapq_cigar_or_conflicting_complete_methylation_payload"
        ),
        "tag_authority": "latest_LongPhase-S_external_coordinate_sidecar",
        "pass": passed,
    }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if not passed:
        raise RuntimeError(f"{sample} final BAM quickcheck/index validation failed")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--preflight",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "latest_input_preflight.json",
    )
    parser.add_argument("--sample", action="append", choices=DATASETS)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--compression-threads", type=int, default=4)
    parser.add_argument("--force-new-suffix", default="")
    parser.add_argument("--summary", type=Path, default=None)
    args = parser.parse_args()
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    selected = set(args.sample or DATASETS)
    entries = [entry for entry in preflight["samples"] if entry["sample"] in selected]
    if len(entries) != len(selected):
        raise RuntimeError("Requested sample missing from preflight")
    workspace = preflight["workspace"]
    receipts: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {
            executor.submit(
                materialize_sample,
                entry,
                workspace,
                args.compression_threads,
                args.force_new_suffix,
            ): entry["sample"]
            for entry in entries
        }
        for future in as_completed(futures):
            sample = futures[future]
            try:
                receipt = future.result()
                receipts.append(receipt)
                print(
                    f"[{sample}] PASS written={receipt['diagnostics']['written_unique_alignments']} "
                    f"size={receipt['output_size_bytes']} MM+ML={receipt['diagnostics']['with_MM_ML']}",
                    flush=True,
                )
            except Exception as error:  # fail all-scope summary while preserving per-sample partial artifact
                failures.append({"sample": sample, "error": repr(error)})
                print(f"[{sample}] FAIL {error!r}", flush=True)
    receipts.sort(key=lambda value: DATASETS.index(value["sample"]))
    summary = {
        "schema_name": "intersubmod.latest_tagged_subset_materialization_batch",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "preflight": str(args.preflight),
        "workspace": workspace,
        "requested_samples": sorted(selected, key=DATASETS.index),
        "workers": args.workers,
        "compression_threads_per_worker": args.compression_threads,
        "receipts": receipts,
        "failures": failures,
        "totals": {
            "written_unique_alignments": sum(r["diagnostics"]["written_unique_alignments"] for r in receipts),
            "duplicate_identity_collapsed": sum(r["diagnostics"]["duplicate_identity_collapsed"] for r in receipts),
            "with_MM_ML": sum(r["diagnostics"]["with_MM_ML"] for r in receipts),
            "output_size_bytes": sum(r["output_size_bytes"] for r in receipts),
        },
        "pass": len(receipts) == len(entries) and not failures and all(r["pass"] for r in receipts),
    }
    summary_path = args.summary or (args.preflight.parent / "latest_tagged_subset_materialization.json")
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), "totals": summary["totals"], "pass": summary["pass"]}, indent=2))
    raise SystemExit(0 if summary["pass"] else 1)


if __name__ == "__main__":
    main()
