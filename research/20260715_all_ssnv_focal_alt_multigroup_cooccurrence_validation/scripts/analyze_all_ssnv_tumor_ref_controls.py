#!/usr/bin/env python3
"""Run tumor-REF and joint ALT+REF controls for every primary stable focal-ALT site."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
from collections import Counter, deque
from concurrent.futures import Future, ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np


TOPIC_ROOT = Path(__file__).resolve().parents[1]
FP_SCRIPT_ROOT = (
    TOPIC_ROOT.parent / "20260715_single_fp_alt_multicluster_subclone_validation" / "scripts"
)
if str(FP_SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(FP_SCRIPT_ROOT))

import focal_alt_cluster_lib as F  # noqa: E402


SCHEMA_NAME = "intersubmod.all_ssnv_tumor_ref_controls"
SCHEMA_VERSION = "1.0.0"
PRIMARY_ASSIGNMENT_SCHEMA = "intersubmod.all_ssnv_stable_multigroup_read_assignment"
SCREEN_CONTRACT = "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"
JOINT_ALLELE_PERMUTATIONS = 499
PRIMARY_SCREEN_THRESHOLDS = {
    "minimum_group_size": 3,
    "minimum_between_within_ratio": 1.3,
    "rnull": 40,
    "modal_confidence": 0.7,
    "hidden_heterogeneity_fraction": 0.30,
}

SiteKey = tuple[str, str, int, str, str]

CONTROL_FIELDS = [
    "primary_assignment_n_core_reads",
    "primary_assignment_labels_sha256",
    "n_tumor_alt",
    "n_tumor_ref",
    "n_tumor_alt_matrix",
    "n_tumor_ref_matrix",
    "ref_status",
    "ref_evaluable",
    "ref_not_testable_reason",
    "ref_n_after_peel",
    "ref_n_peeled",
    "ref_coarse_ng",
    "ref_modal_fraction",
    "ref_unstable",
    "ref_stable_null_multigroup",
    "ref_cluster_sizes",
    "joint_status",
    "joint_evaluable",
    "joint_not_testable_reason",
    "joint_n_after_peel",
    "joint_n_peeled",
    "joint_coarse_ng",
    "joint_modal_fraction",
    "joint_unstable",
    "joint_stable_null_multigroup",
    "joint_cluster_sizes",
    "joint_allele_testable",
    "joint_allele_v",
    "joint_allele_p_perm",
    "joint_allele_n",
    "joint_allele_axis_aligned",
    "joint_allele_not_testable_reason",
    "background_control_interpretation",
    "component_dedup_key",
    "component_representative",
    "readset_dedup_key",
    "readset_representative",
]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def json_text(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def open_text(path: Path) -> Any:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def primary_flag(value: Any) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"Invalid stable_null_multigroup value: {value!r}")


def validate_primary_screen_contract() -> None:
    observed = {
        "minimum_group_size": F.MIN_SIZE,
        "minimum_between_within_ratio": F.SEP_MIN,
        "rnull": F.RNULL,
        "modal_confidence": F.MODAL_CONFIDENCE,
        "hidden_heterogeneity_fraction": F.HIDDEN_HETEROGENEITY_FRACTION,
    }
    if observed != PRIMARY_SCREEN_THRESHOLDS:
        raise RuntimeError(
            f"Primary phylo-v4.1 threshold drift: observed={observed} "
            f"expected={PRIMARY_SCREEN_THRESHOLDS}"
        )


def exact_site_key(sample: Any, chrom: Any, pos: Any, ref: Any, alt: Any) -> SiteKey:
    values = {
        "sample": sample,
        "chrom": chrom,
        "ref": ref,
        "alt": alt,
    }
    missing = [name for name, value in values.items() if value is None or str(value) == ""]
    if missing:
        raise ValueError(f"Site identity has empty fields: {missing}")
    return str(sample), str(chrom), int(pos), str(ref), str(alt)


def site_row_key(row: Mapping[str, Any]) -> SiteKey:
    return exact_site_key(row.get("sample"), row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))


def assignment_row_key(record: Mapping[str, Any]) -> SiteKey:
    posthoc = record.get("posthoc")
    if not isinstance(posthoc, Mapping):
        raise ValueError("Stable assignment record is missing posthoc metadata")
    for allele in ("ref", "alt"):
        if allele in record and record[allele] != posthoc.get(allele):
            raise RuntimeError(
                f"Stable assignment top-level/posthoc {allele} mismatch: "
                f"{record[allele]!r} != {posthoc.get(allele)!r}"
            )
    return exact_site_key(
        record.get("sample"),
        record.get("chrom"),
        record.get("pos"),
        posthoc.get("ref"),
        posthoc.get("alt"),
    )


def _format_keys(keys: Iterable[SiteKey], limit: int = 3) -> list[str]:
    return [f"{sample}:{chrom}:{pos}:{ref}>{alt}" for sample, chrom, pos, ref, alt in list(keys)[:limit]]


def load_stable_assignments(path: Path) -> dict[SiteKey, dict[str, Any]]:
    assignments: dict[SiteKey, dict[str, Any]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            if not isinstance(record, dict):
                raise ValueError(f"Assignment line {line_number} is not a JSON object")
            schema = record.get("schema_name")
            if schema != PRIMARY_ASSIGNMENT_SCHEMA or record.get("schema_version") != "1.0.0":
                raise RuntimeError(f"Unexpected stable assignment schema on line {line_number}: {schema!r}")
            screen_contract = record.get("screen_contract")
            if screen_contract != SCREEN_CONTRACT:
                raise RuntimeError(
                    f"Stable assignment screen-contract drift on line {line_number}: "
                    f"{screen_contract!r}"
                )
            key = assignment_row_key(record)
            if key in assignments:
                raise RuntimeError(f"Duplicate stable assignment site key: {_format_keys([key])[0]}")
            region_dir = record.get("region_dir")
            if not isinstance(region_dir, str) or not region_dir:
                raise ValueError(f"Stable assignment lacks region_dir at {_format_keys([key])[0]}")
            read_ids = record.get("read_ids")
            labels = record.get("labels")
            if not isinstance(read_ids, list) or not isinstance(labels, list):
                raise ValueError(f"Stable assignment lacks read_ids/labels at {_format_keys([key])[0]}")
            if len(read_ids) != len(labels):
                raise RuntimeError(f"Stable assignment read/label length mismatch at {_format_keys([key])[0]}")
            if len(read_ids) != len(set(str(read_id) for read_id in read_ids)):
                raise RuntimeError(f"Duplicate primary assignment read ID at {_format_keys([key])[0]}")
            label_payload = [[str(read_id), str(label)] for read_id, label in zip(read_ids, labels)]
            assignments[key] = {
                "region_dir": region_dir,
                "n_core_reads": len(read_ids),
                "labels_sha256": hashlib.sha256(json_text(label_payload).encode()).hexdigest(),
                "posthoc": dict(record["posthoc"]),
            }
    return assignments


def load_primary_inputs(site_results: Path, stable_assignments: Path) -> tuple[list[str], list[dict[str, Any]]]:
    """Load the complete primary stable set and require exact assignment-key equality."""
    assignments = load_stable_assignments(stable_assignments)
    stable_rows: dict[SiteKey, dict[str, str]] = {}
    ordered_keys: list[SiteKey] = []
    with open_text(site_results) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_fields = list(reader.fieldnames or [])
        required = {
            "sample",
            "chrom",
            "pos",
            "ref",
            "alt",
            "region_dir",
            "stable_null_multigroup",
        }
        missing = required - set(source_fields)
        if missing:
            raise ValueError(f"{site_results} missing required columns: {sorted(missing)}")
        collisions = sorted(set(source_fields).intersection(CONTROL_FIELDS))
        if collisions:
            raise RuntimeError(
                "Primary ALT screen fields would be overwritten by control output: "
                f"{collisions}"
            )
        for row in reader:
            if not primary_flag(row["stable_null_multigroup"]):
                continue
            screen_contract = row.get("screen_contract")
            if screen_contract != SCREEN_CONTRACT:
                raise RuntimeError(
                    f"Primary stable site screen-contract drift: {screen_contract!r}"
                )
            key = site_row_key(row)
            if key in stable_rows:
                raise RuntimeError(f"Duplicate primary stable site key: {_format_keys([key])[0]}")
            stable_rows[key] = dict(row)
            ordered_keys.append(key)

    site_keys = set(stable_rows)
    assignment_keys = set(assignments)
    if site_keys != assignment_keys:
        missing_assignments = sorted(site_keys - assignment_keys)
        extra_assignments = sorted(assignment_keys - site_keys)
        raise RuntimeError(
            "Primary stable site/assignment key mismatch: "
            f"stable_sites={len(site_keys)} assignments={len(assignment_keys)} "
            f"missing_assignments={_format_keys(missing_assignments)} "
            f"extra_assignments={_format_keys(extra_assignments)}"
        )

    tasks: list[dict[str, Any]] = []
    posthoc_fields = (
        "biological_id",
        "truth_label",
        "ref",
        "alt",
        "ssnv_branch",
        "component_id",
        "component_size",
        "selected_group_id",
        "selected_group_size",
    )
    for key in ordered_keys:
        site = stable_rows[key]
        assignment = assignments[key]
        if site["region_dir"] != assignment["region_dir"]:
            raise RuntimeError(
                f"Primary site/assignment region_dir mismatch at {_format_keys([key])[0]}: "
                f"{site['region_dir']!r} != {assignment['region_dir']!r}"
            )
        for field in posthoc_fields:
            if field not in site or field not in assignment["posthoc"]:
                continue
            site_value = site[field] if site[field] != "" else None
            assignment_value = assignment["posthoc"][field]
            if site_value is not None and assignment_value is not None:
                if str(site_value) != str(assignment_value):
                    raise RuntimeError(
                        f"Primary site/assignment posthoc {field} mismatch at "
                        f"{_format_keys([key])[0]}: {site_value!r} != {assignment_value!r}"
                    )
        tasks.append(
            {
                "site": site,
                "primary_assignment_n_core_reads": assignment["n_core_reads"],
                "primary_assignment_labels_sha256": assignment["labels_sha256"],
            }
        )
    return source_fields, tasks


def load_reads_unique(path: Path) -> dict[str, dict[str, str]]:
    reads: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"read_id", "is_tumor", "alt_support"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing required read columns: {sorted(missing)}")
        for row in reader:
            read_id = row["read_id"]
            if read_id in reads:
                raise RuntimeError(f"Duplicate reads.tsv read_id: {path} {read_id}")
            reads[read_id] = row
    return reads


def select_tumor_read_ids(
    reads: Mapping[str, Mapping[str, Any]], allele: str, ordered_ids: Iterable[str] | None = None
) -> list[str]:
    if allele not in {"ALT", "REF"}:
        raise ValueError(f"Unsupported focal allele: {allele!r}")
    candidates = ordered_ids if ordered_ids is not None else reads.keys()
    return [
        read_id
        for read_id in candidates
        if read_id in reads
        and F.is_tumor(str(reads[read_id].get("is_tumor", "")))
        and reads[read_id].get("alt_support") == allele
    ]


def validate_matrices(
    distance_ids: Sequence[str],
    distance: np.ndarray,
    methylation_ids: Sequence[str],
    methylation: np.ndarray,
) -> None:
    if len(distance_ids) != len(set(distance_ids)):
        raise RuntimeError("Duplicate read IDs in BERNOULLI matrix")
    if len(methylation_ids) != len(set(methylation_ids)):
        raise RuntimeError("Duplicate read IDs in methylation matrix")
    if distance.shape != (len(distance_ids), len(distance_ids)):
        raise RuntimeError(
            f"Invalid BERNOULLI matrix shape: {distance.shape} for {len(distance_ids)} read IDs"
        )
    if methylation.ndim != 2 or methylation.shape[0] != len(methylation_ids):
        raise RuntimeError(
            f"Invalid methylation matrix shape: {methylation.shape} for {len(methylation_ids)} read IDs"
        )


def _not_evaluable(
    status: str, reason: str, n_raw: int, n_matrix: int, n_after_peel: int = 0
) -> dict[str, Any]:
    return {
        "status": status,
        "evaluable": False,
        "not_testable_reason": reason,
        "n_raw": n_raw,
        "n_matrix": n_matrix,
        "n_after_peel": n_after_peel,
        "n_peeled": (
            max(0, n_matrix - n_after_peel)
            if status == "incomplete_distance_below_minimum"
            else 0
        ),
        "coarse_ng": None,
        "modal_fraction": None,
        "unstable": None,
        "stable_null_multigroup": False,
        "cluster_sizes": {},
    }


def analyze_subset(
    read_ids: Sequence[str],
    distance_ids: Sequence[str],
    distance: np.ndarray,
    methylation_ids: Sequence[str],
    methylation: np.ndarray,
    seed: int,
) -> dict[str, Any]:
    """Apply the primary phylo-v4.1 column-null screen to a selected read subset."""
    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    usable = [
        read_id for read_id in read_ids if read_id in distance_index and read_id in methylation_index
    ]
    if len(read_ids) < 2 * F.MIN_SIZE:
        return _not_evaluable(
            "insufficient_reads",
            f"fewer_than_2x_MIN_SIZE_raw_reads:{len(read_ids)}<{2 * F.MIN_SIZE}",
            len(read_ids),
            len(usable),
        )
    if len(usable) < 2 * F.MIN_SIZE:
        return _not_evaluable(
            "insufficient_matrix_reads",
            f"fewer_than_2x_MIN_SIZE_matrix_reads:{len(usable)}<{2 * F.MIN_SIZE}",
            len(read_ids),
            len(usable),
        )

    distance_rows = [distance_index[read_id] for read_id in usable]
    sub_distance = distance[np.ix_(distance_rows, distance_rows)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [usable[index] for index in kept]
    if len(kept_ids) < 2 * F.MIN_SIZE:
        return _not_evaluable(
            "incomplete_distance_below_minimum",
            f"fewer_than_2x_MIN_SIZE_after_peel:{len(kept_ids)}<{2 * F.MIN_SIZE}",
            len(read_ids),
            len(usable),
            len(kept_ids),
        )

    sub_distance = sub_distance[np.ix_(kept, kept)]
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]
    phylo = F.analyze_phylo(
        sub_distance,
        sub_methylation,
        base_seed=seed,
        seeds=10,
        rnull=F.RNULL,
        null_mode="column",
        empirical_alpha=None,
    )
    labels = phylo["coarse_labels"]
    stable = phylo["coarse_ng"] >= 2 and not phylo["unstable"]
    return {
        "status": "evaluable",
        "evaluable": True,
        "not_testable_reason": None,
        "n_raw": len(read_ids),
        "n_matrix": len(usable),
        "n_after_peel": len(kept_ids),
        "n_peeled": len(usable) - len(kept_ids),
        "kept_ids": kept_ids,
        "labels": labels,
        "coarse_ng": phylo["coarse_ng"],
        "modal_fraction": phylo["modal_fraction"],
        "unstable": phylo["unstable"],
        "stable_null_multigroup": stable,
        "cluster_sizes": dict(
            Counter(label for label in labels if label not in {"other", "outlier"})
        ),
    }


def compute_joint_allele_association(
    kept_ids: Sequence[str],
    labels: Sequence[str],
    allele_by_id: Mapping[str, str],
    seed: int,
    permutations: int = JOINT_ALLELE_PERMUTATIONS,
) -> dict[str, Any]:
    if len(kept_ids) != len(labels):
        raise RuntimeError("Joint clustering kept-read/label length mismatch")
    core = [
        (read_id, label)
        for read_id, label in zip(kept_ids, labels)
        if label not in {"other", "outlier"}
    ]
    missing = [read_id for read_id, _ in core if read_id not in allele_by_id]
    if missing:
        raise RuntimeError(f"Joint core reads lack focal allele labels: {missing[:3]}")
    association = F.categorical_permutation_association(
        [allele_by_id[read_id] for read_id, _ in core],
        [label for _, label in core],
        seed,
        permutations=permutations,
    )
    testable = association["v"] is not None and association["p_perm"] is not None
    return {
        "testable": testable,
        "v": association["v"],
        "p_perm": association["p_perm"],
        "n": association["n"],
        "aligned": bool(association["aligned"]),
        "not_testable_reason": (
            None if testable else "joint_core_lacks_both_focal_alleles_or_multiple_clusters"
        ),
    }


def _copy_cluster_result(prefix: str, result: Mapping[str, Any], output: dict[str, Any]) -> None:
    output.update(
        {
            f"{prefix}_status": result["status"],
            f"{prefix}_evaluable": bool(result["evaluable"]),
            f"{prefix}_not_testable_reason": result["not_testable_reason"],
            f"{prefix}_n_after_peel": result["n_after_peel"],
            f"{prefix}_n_peeled": result["n_peeled"],
            f"{prefix}_coarse_ng": result["coarse_ng"],
            f"{prefix}_modal_fraction": result["modal_fraction"],
            f"{prefix}_unstable": result["unstable"],
            f"{prefix}_stable_null_multigroup": bool(result["stable_null_multigroup"]),
            f"{prefix}_cluster_sizes": json_text(result["cluster_sizes"]),
        }
    )


def _validate_primary_alt_counts(
    site: Mapping[str, Any], n_alt_raw: int, n_alt_matrix: int
) -> None:
    for field, observed in (("n_alt_raw", n_alt_raw), ("n_alt_matrix", n_alt_matrix)):
        expected = site.get(field)
        if expected not in {None, ""} and int(expected) != observed:
            raise RuntimeError(
                f"Primary ALT count drift for {field}: site_results={expected} region={observed}"
            )


def process_site_task(task: dict[str, Any]) -> dict[str, Any]:
    site = task["site"]
    key = site_row_key(site)
    try:
        if not primary_flag(site["stable_null_multigroup"]):
            raise RuntimeError("Worker received a non-primary-stable site")
        region_dir = Path(site["region_dir"])
        reads_path = region_dir / "reads" / "reads.tsv"
        distance_path = region_dir / "distance" / "BERNOULLI" / "matrix.csv"
        methylation_path = region_dir / "methylation" / "methylation.csv"
        for required in (reads_path, distance_path, methylation_path):
            if not required.exists() or required.stat().st_size == 0:
                raise FileNotFoundError(required)

        reads = load_reads_unique(reads_path)
        distance_ids, distance = F.load_matrix(distance_path)
        methylation_ids, methylation = F.load_matrix(methylation_path)
        validate_matrices(distance_ids, distance, methylation_ids, methylation)
        methylation_id_set = set(methylation_ids)

        tumor_alt = select_tumor_read_ids(reads, "ALT")
        tumor_ref = select_tumor_read_ids(reads, "REF")
        tumor_alt_matrix = select_tumor_read_ids(
            reads, "ALT", (read_id for read_id in distance_ids if read_id in methylation_id_set)
        )
        tumor_ref_matrix = select_tumor_read_ids(
            reads, "REF", (read_id for read_id in distance_ids if read_id in methylation_id_set)
        )
        _validate_primary_alt_counts(site, len(tumor_alt), len(tumor_alt_matrix))

        seed = F.stable_seed(key[0], key[1], key[2])
        ref_result = analyze_subset(
            tumor_ref,
            distance_ids,
            distance,
            methylation_ids,
            methylation,
            seed + 100_000,
        )
        output: dict[str, Any] = dict(site)
        output.update(
            {
                "primary_assignment_n_core_reads": task["primary_assignment_n_core_reads"],
                "primary_assignment_labels_sha256": task[
                    "primary_assignment_labels_sha256"
                ],
                "n_tumor_alt": len(tumor_alt),
                "n_tumor_ref": len(tumor_ref),
                "n_tumor_alt_matrix": len(tumor_alt_matrix),
                "n_tumor_ref_matrix": len(tumor_ref_matrix),
            }
        )
        _copy_cluster_result("ref", ref_result, output)

        joint_ids = [
            read_id
            for read_id in distance_ids
            if read_id in methylation_id_set
            and read_id in reads
            and F.is_tumor(reads[read_id]["is_tumor"])
            and reads[read_id]["alt_support"] in {"ALT", "REF"}
        ]
        if len(tumor_ref_matrix) < F.MIN_SIZE:
            joint_result = _not_evaluable(
                "insufficient_ref_for_joint_control",
                f"fewer_than_MIN_SIZE_tumor_REF_matrix_reads:{len(tumor_ref_matrix)}<{F.MIN_SIZE}",
                len(tumor_alt) + len(tumor_ref),
                len(joint_ids),
            )
        else:
            joint_result = analyze_subset(
                joint_ids,
                distance_ids,
                distance,
                methylation_ids,
                methylation,
                seed + 200_000,
            )
        _copy_cluster_result("joint", joint_result, output)

        allele_result = {
            "testable": False,
            "v": None,
            "p_perm": None,
            "n": 0,
            "aligned": False,
            "not_testable_reason": (
                joint_result["not_testable_reason"]
                if not joint_result["evaluable"]
                else "joint_clustering_not_stable_multigroup"
            ),
        }
        if joint_result["evaluable"] and joint_result["stable_null_multigroup"]:
            allele_by_id = {
                read_id: reads[read_id]["alt_support"]
                for read_id in joint_result["kept_ids"]
            }
            allele_result = compute_joint_allele_association(
                joint_result["kept_ids"],
                joint_result["labels"],
                allele_by_id,
                seed + 300_000,
            )
        output.update(
            {
                "joint_allele_testable": allele_result["testable"],
                "joint_allele_v": allele_result["v"],
                "joint_allele_p_perm": allele_result["p_perm"],
                "joint_allele_n": allele_result["n"],
                "joint_allele_axis_aligned": allele_result["aligned"],
                "joint_allele_not_testable_reason": allele_result[
                    "not_testable_reason"
                ],
                "background_control_interpretation": (
                    "REF_REPLICATION_WEAKENS_ALT_SPECIFICITY"
                    if ref_result["stable_null_multigroup"]
                    else "REF_NONREPLICATION_DOES_NOT_PROVE_SUBCLONE"
                    if ref_result["evaluable"]
                    else "REF_CONTROL_NOT_TESTABLE"
                ),
            }
        )
        return output
    except Exception as error:
        raise RuntimeError(f"{_format_keys([key])[0]} tumor-REF control failed: {error}") from error


def process_chunk(chunk: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [process_task(task) for task in chunk]


def process_task(task: dict[str, Any]) -> dict[str, Any]:
    """Compatibility entry point matching the earlier REF-control analyzer."""
    return process_site_task(task)


def chunked(values: Iterable[dict[str, Any]], size: int) -> Iterator[list[dict[str, Any]]]:
    if size < 1:
        raise ValueError("Chunk size must be positive")
    chunk: list[dict[str, Any]] = []
    for value in values:
        chunk.append(value)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def bounded_chunk_results(
    executor: Any,
    chunks: Iterable[list[dict[str, Any]]],
    max_pending: int,
) -> Iterator[list[dict[str, Any]]]:
    """Submit at most max_pending chunks; Future.result propagates every site failure."""
    if max_pending < 1:
        raise ValueError("max_pending must be positive")
    iterator = iter(chunks)
    pending: deque[Future[list[dict[str, Any]]]] = deque()

    def fill() -> None:
        while len(pending) < max_pending:
            try:
                pending.append(executor.submit(process_chunk, next(iterator)))
            except StopIteration:
                break

    fill()
    while pending:
        yield pending.popleft().result()
        fill()


def component_key(row: Mapping[str, Any]) -> tuple[str, str]:
    biological_id = str(row.get("biological_id") or row.get("sample") or "NA")
    component = row.get("component_id")
    if component in {None, ""}:
        site = site_row_key(row)
        component = f"SITE:{site[1]}:{site[2]}:{site[3]}>{site[4]}"
    return biological_id, str(component)


def readset_key(row: Mapping[str, Any]) -> tuple[str, str]:
    biological_id = str(row.get("biological_id") or row.get("sample") or "NA")
    readset = row.get("alt_readset_sha256")
    if readset in {None, ""}:
        site = site_row_key(row)
        readset = f"MISSING:{site[0]}:{site[1]}:{site[2]}:{site[3]}>{site[4]}"
    return biological_id, str(readset)


class ControlSummary:
    def __init__(self) -> None:
        self.n_sites = 0
        self.n_tumor_alt = 0
        self.n_tumor_ref = 0
        self.ref_evaluable = 0
        self.ref_stable = 0
        self.joint_evaluable = 0
        self.joint_stable = 0
        self.joint_allele_testable = 0
        self.joint_allele_aligned = 0
        self.ref_reasons: Counter[str] = Counter()
        self.joint_reasons: Counter[str] = Counter()
        self.allele_reasons: Counter[str] = Counter()

    def add(self, row: Mapping[str, Any]) -> None:
        self.n_sites += 1
        self.n_tumor_alt += int(row["n_tumor_alt"])
        self.n_tumor_ref += int(row["n_tumor_ref"])
        self.ref_evaluable += bool(row["ref_evaluable"])
        self.ref_stable += bool(row["ref_stable_null_multigroup"])
        self.joint_evaluable += bool(row["joint_evaluable"])
        self.joint_stable += bool(row["joint_stable_null_multigroup"])
        self.joint_allele_testable += bool(row["joint_allele_testable"])
        self.joint_allele_aligned += bool(row["joint_allele_axis_aligned"])
        if row.get("ref_not_testable_reason"):
            self.ref_reasons[str(row["ref_not_testable_reason"])] += 1
        if row.get("joint_not_testable_reason"):
            self.joint_reasons[str(row["joint_not_testable_reason"])] += 1
        if row.get("joint_allele_not_testable_reason"):
            self.allele_reasons[str(row["joint_allele_not_testable_reason"])] += 1

    def payload(self) -> dict[str, Any]:
        return {
            "n_sites": self.n_sites,
            "n_tumor_alt_reads": self.n_tumor_alt,
            "n_tumor_ref_reads": self.n_tumor_ref,
            "n_ref_evaluable": self.ref_evaluable,
            "ref_evaluable_fraction": self.ref_evaluable / self.n_sites if self.n_sites else None,
            "n_ref_stable_null_multigroup": self.ref_stable,
            "ref_stable_fraction_evaluable": (
                self.ref_stable / self.ref_evaluable if self.ref_evaluable else None
            ),
            "n_joint_evaluable": self.joint_evaluable,
            "joint_evaluable_fraction": self.joint_evaluable / self.n_sites if self.n_sites else None,
            "n_joint_stable_null_multigroup": self.joint_stable,
            "joint_stable_fraction_evaluable": (
                self.joint_stable / self.joint_evaluable if self.joint_evaluable else None
            ),
            "n_joint_allele_testable": self.joint_allele_testable,
            "n_joint_allele_axis_aligned": self.joint_allele_aligned,
            "joint_allele_aligned_fraction_testable": (
                self.joint_allele_aligned / self.joint_allele_testable
                if self.joint_allele_testable
                else None
            ),
            "ref_not_testable_reason_counts": dict(sorted(self.ref_reasons.items())),
            "joint_not_testable_reason_counts": dict(sorted(self.joint_reasons.items())),
            "joint_allele_not_testable_reason_counts": dict(sorted(self.allele_reasons.items())),
        }


class StratumSummary:
    def __init__(self) -> None:
        self.site_weighted = ControlSummary()
        self.component_deduplicated = ControlSummary()
        self.readset_deduplicated = ControlSummary()
        self._components: set[tuple[str, str]] = set()
        self._readsets: set[tuple[str, str]] = set()

    def add(self, row: Mapping[str, Any]) -> None:
        self.site_weighted.add(row)
        component = component_key(row)
        if component not in self._components:
            self._components.add(component)
            self.component_deduplicated.add(row)
        readset = readset_key(row)
        if readset not in self._readsets:
            self._readsets.add(readset)
            self.readset_deduplicated.add(row)

    def payload(self) -> dict[str, Any]:
        return {
            "site_weighted": self.site_weighted.payload(),
            "component_deduplicated": self.component_deduplicated.payload(),
            "readset_deduplicated": self.readset_deduplicated.payload(),
        }


class SummaryBook:
    def __init__(self) -> None:
        self.pooled = StratumSummary()
        self.per_sample: dict[str, StratumSummary] = {}
        self.per_truth: dict[str, StratumSummary] = {}
        self.per_biological_id: dict[str, StratumSummary] = {}

    @staticmethod
    def _add(grouped: dict[str, StratumSummary], key: str, row: Mapping[str, Any]) -> None:
        if key not in grouped:
            grouped[key] = StratumSummary()
        grouped[key].add(row)

    def add(self, row: Mapping[str, Any]) -> None:
        self.pooled.add(row)
        self._add(self.per_sample, str(row.get("sample") or "NA"), row)
        self._add(self.per_truth, str(row.get("truth_label") or "NA"), row)
        self._add(
            self.per_biological_id,
            str(row.get("biological_id") or row.get("sample") or "NA"),
            row,
        )

    def payload(self) -> dict[str, Any]:
        return {
            "pooled": self.pooled.payload(),
            "per_sample": {
                key: value.payload() for key, value in sorted(self.per_sample.items())
            },
            "per_truth": {
                key: value.payload() for key, value in sorted(self.per_truth.items())
            },
            "per_biological_id": {
                key: value.payload()
                for key, value in sorted(self.per_biological_id.items())
            },
        }


class DedupAnnotator:
    def __init__(self) -> None:
        self.components: set[tuple[str, str]] = set()
        self.readsets: set[tuple[str, str]] = set()

    def annotate(self, row: dict[str, Any]) -> None:
        component = component_key(row)
        readset = readset_key(row)
        row["component_dedup_key"] = json_text(component)
        row["component_representative"] = component not in self.components
        row["readset_dedup_key"] = json_text(readset)
        row["readset_representative"] = readset not in self.readsets
        self.components.add(component)
        self.readsets.add(readset)


def _tsv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list, tuple)):
        return json_text(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def create_output_dir(path: Path) -> None:
    if os.path.lexists(path):
        raise FileExistsError(f"Refusing to overwrite existing output directory: {path}")
    path.mkdir(parents=True, exist_ok=False)


def parse_args() -> argparse.Namespace:
    primary_root = TOPIC_ROOT / "results" / "all_ssnv_focal_alt_multigroup_v1"
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--site-results",
        type=Path,
        default=primary_root / "all_ssnv_site_results.tsv.gz",
    )
    parser.add_argument(
        "--stable-assignments",
        "--assignments",
        dest="stable_assignments",
        type=Path,
        default=primary_root / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=TOPIC_ROOT / "results" / "all_ssnv_tumor_ref_controls_v1",
    )
    parser.add_argument("--primary-artifact-audit-pre", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=max(1, min(32, (os.cpu_count() or 4) - 4)))
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--max-pending-chunks", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_at = now_utc()
    if args.workers < 1 or args.chunk_size < 1 or args.progress_every < 1:
        raise ValueError("workers, chunk-size, and progress-every must be positive")
    if args.max_pending_chunks < 0:
        raise ValueError("max-pending-chunks must be nonnegative")
    validate_primary_screen_contract()
    input_artifacts_before = {
        "site_results": artifact(args.site_results),
        "stable_assignments": artifact(args.stable_assignments),
        "primary_artifact_audit_pre": artifact(args.primary_artifact_audit_pre),
    }
    source_fields, tasks = load_primary_inputs(args.site_results, args.stable_assignments)
    if not tasks:
        raise RuntimeError("No primary stable_null_multigroup=true sites were found")
    create_output_dir(args.output_dir)

    site_path = args.output_dir / "all_ssnv_tumor_ref_control_site_results.tsv.gz"
    summary_path = args.output_dir / "all_ssnv_tumor_ref_control_summary.json"
    manifest_path = args.output_dir / "run_manifest.json"
    fields = [*source_fields, *CONTROL_FIELDS]
    max_pending = args.max_pending_chunks or max(1, args.workers * 2)
    processed = 0
    summaries = SummaryBook()
    dedup = DedupAnnotator()

    with gzip.open(site_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        chunks = chunked(tasks, args.chunk_size)
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for result_chunk in bounded_chunk_results(executor, chunks, max_pending):
                for row in result_chunk:
                    dedup.annotate(row)
                    writer.writerow({field: _tsv_value(row[field]) for field in fields})
                    summaries.add(row)
                    processed += 1
                    if processed % args.progress_every == 0 or processed == len(tasks):
                        print(f"processed={processed}/{len(tasks)}", flush=True)

    if processed != len(tasks):
        raise RuntimeError(f"Processed stable-site count mismatch: {processed} != {len(tasks)}")
    summary = {
        "schema_name": f"{SCHEMA_NAME}.summary",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "task_type": "B_comprehensive_validation",
        "scope": "complete_primary_stable_null_multigroup_set",
        "selection_contract": {
            "site_results_filter": "stable_null_multigroup=true",
            "site_key": ["sample", "chrom", "pos", "ref", "alt"],
            "site_assignment_key_sets_exact": True,
            "site_assignment_posthoc_metadata_exact": True,
            "selected_site_duplicates": 0,
            "assignment_duplicates": 0,
            "primary_alt_labels_rewritten": False,
        },
        "clustering_contract": {
            "screen": SCREEN_CONTRACT,
            "minimum_group_size": F.MIN_SIZE,
            "rnull": F.RNULL,
            "seeds": 10,
            "null_mode": "column",
            "ref_selection": "reads.tsv is_tumor=true and alt_support=REF",
            "joint_selection": "reads.tsv is_tumor=true and alt_support in {ALT,REF}",
            "joint_allele_association": "posthoc Cramer's V with 499 label permutations after stable joint clustering",
        },
        "denominator_contract": {
            "site_weighted": "one row per dataset site",
            "component_deduplicated": ["biological_id", "component_id"],
            "readset_deduplicated": ["biological_id", "alt_readset_sha256"],
            "missing_component_or_readset": "site-specific fail-open identity; never merged",
        },
        **summaries.payload(),
        "guardrail": {
            "ref_nonreplication": "Tumor-REF nonreplication does not prove a subclone.",
            "ref_replication": "Tumor-REF replication weakens ALT-specificity.",
        },
        "pass": processed == len(tasks),
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    input_artifacts = {
        "site_results": artifact(args.site_results),
        "stable_assignments": artifact(args.stable_assignments),
        "primary_artifact_audit_pre": artifact(args.primary_artifact_audit_pre),
    }
    if input_artifacts != input_artifacts_before:
        raise RuntimeError("Tumor-REF input artifacts changed during analysis")
    source_code_artifacts = {
        "analyzer": artifact(Path(__file__)),
        "focal_alt_cluster_lib": artifact(Path(F.__file__)),
    }
    output_artifacts = {
        "site_results": artifact(site_path),
        "summary": artifact(summary_path),
    }
    finished_at = now_utc()
    manifest = {
        "schema_name": f"{SCHEMA_NAME}.run_manifest",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "command": sys.argv,
        "inputs": input_artifacts,
        "source_code": source_code_artifacts,
        "execution": {
            "workers": args.workers,
            "chunk_size": args.chunk_size,
            "max_pending_chunks": max_pending,
            "per_site_future_submission": False,
            "site_exceptions_hard_fail": True,
        },
        "contracts": {
            "complete_primary_stable_set_only": True,
            "site_assignment_keys_exact": True,
            "site_assignment_posthoc_metadata_exact": True,
            "primary_alt_labels_rewritten": False,
            "existing_results_overwritten": False,
            "guardrail": summary["guardrail"],
        },
        "outputs": output_artifacts,
        "counts": {
            "primary_stable_sites": len(tasks),
            "processed_sites": processed,
        },
        "pass": summary["pass"],
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir.resolve()),
                "site_results": str(site_path.resolve()),
                "summary": str(summary_path.resolve()),
                "run_manifest": str(manifest_path.resolve()),
                "processed_sites": processed,
                "pass": summary["pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
