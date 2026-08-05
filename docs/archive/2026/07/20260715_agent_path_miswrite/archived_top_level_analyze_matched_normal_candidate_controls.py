#!/usr/bin/env python3
"""Analyze focal-REF matched-normal methyl backgrounds for candidate sites.

Normal allele counts are reported across every normal read, but the methylation
background is estimated only from focal REF-called normal reads. A negative
normal control is meaningful only when the focal callability and REF-only
methylation evaluability gates both pass.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
FP_SCRIPT_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_single_fp_alt_multicluster_subclone_validation"
    / "scripts"
)
if str(FP_SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(FP_SCRIPT_ROOT))

import focal_alt_cluster_lib as F  # noqa: E402


SCHEMA_VERSION = "2.0.0"
NORMAL_MIN_CALLED_READS = 10
NORMAL_MIN_REF_READS = 5
NORMAL_HP_POLICY = "PROHIBITED_NOT_USED"
NORMAL_METHYL_BACKGROUND_POPULATION = "is_tumor=0_and_focal_call=REF_only"
GUARDRAIL = (
    "A negative matched-normal focal-REF methyl control does not exclude allele-specific "
    "methylation (ASM), cell-state-dependent methylation, limited read power, contamination, "
    "copy-number effects, or mapping artifacts. It does not prove a molecular-haplotype "
    "candidate is a cellular subclone, clone count, or lineage."
)

CandidateKey = tuple[str, str, int, str, str]
AlignmentIdentity = tuple[str, str, int, int, int, str]


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


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def require_nonempty(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")
    return path


def create_output_dir(path: Path) -> Path:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing analysis output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.mkdir(parents=False, exist_ok=False)
    return path


def candidate_key(row: dict[str, Any]) -> CandidateKey:
    try:
        key = (
            str(row["sample"]).strip(),
            str(row["chrom"]).strip(),
            int(row["pos"]),
            str(row["ref"]).strip().upper(),
            str(row["alt"]).strip().upper(),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"Invalid candidate identity: {row!r}") from error
    if not key[0] or not key[1] or key[2] <= 0 or not key[3] or not key[4]:
        raise RuntimeError(f"Invalid candidate identity: {key!r}")
    return key


def parse_output_site(reads_path: Path) -> tuple[str, int]:
    site_name = reads_path.parents[2].name
    if "_" not in site_name:
        raise ValueError(f"Cannot parse paired InterSubMod site directory: {reads_path}")
    chrom, raw_position = site_name.rsplit("_", 1)
    try:
        return chrom, int(raw_position)
    except ValueError as error:
        raise ValueError(f"Cannot parse paired InterSubMod site directory: {reads_path}") from error


def scan_site_outputs(output_dir: Path) -> dict[tuple[str, int], dict[str, str]]:
    outputs: dict[tuple[str, int], dict[str, str]] = {}
    for reads_path in output_dir.rglob("reads.tsv"):
        if reads_path.parent.name != "reads":
            continue
        key = parse_output_site(reads_path)
        if key in outputs:
            raise RuntimeError(f"Duplicate paired InterSubMod output site: {key}")
        region_dir = reads_path.parents[1]
        distance_path = region_dir / "distance" / "BERNOULLI" / "matrix.csv"
        methylation_path = region_dir / "methylation" / "methylation.csv"
        for path in (reads_path, distance_path, methylation_path):
            require_nonempty(path, "paired site artifact")
        outputs[key] = {
            "region_dir": str(region_dir.resolve()),
            "reads_path": str(reads_path.resolve()),
            "distance_path": str(distance_path.resolve()),
            "methylation_path": str(methylation_path.resolve()),
        }
    return outputs


def load_paired_tasks(
    paired_output_root: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any], Path]:
    receipt_path = require_nonempty(
        paired_output_root / "run_receipt.json", "paired runner receipt"
    )
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema_name") != "intersubmod.matched_normal_candidate_control_run":
        raise RuntimeError("Unexpected paired runner receipt schema")
    if receipt.get("pass") is not True:
        raise RuntimeError("Paired runner receipt is not passing")
    sample_receipts = receipt.get("receipts")
    if not isinstance(sample_receipts, list) or not sample_receipts:
        raise RuntimeError("Paired runner receipt has no sample receipts")

    tasks: list[dict[str, Any]] = []
    observed_candidates: set[CandidateKey] = set()
    for sample_receipt in sample_receipts:
        if sample_receipt.get("pass") is not True:
            raise RuntimeError(
                f"Non-passing paired sample receipt: {sample_receipt.get('sample')}"
            )
        sample = str(sample_receipt["sample"])
        outputs = scan_site_outputs(Path(sample_receipt["outputs"]["output_dir"]))
        candidates = sample_receipt.get("candidate_sites")
        if not isinstance(candidates, list) or not candidates:
            raise RuntimeError(f"Paired sample receipt has no candidates: {sample}")
        expected_positions: set[tuple[str, int]] = set()
        for payload in candidates:
            if payload.get("sample") != sample:
                raise RuntimeError(f"Paired receipt candidate sample mismatch: {payload}")
            key = candidate_key(payload)
            if key in observed_candidates:
                raise RuntimeError(f"Duplicate paired candidate identity: {key}")
            observed_candidates.add(key)
            local_position = (key[1], key[2])
            if local_position in expected_positions:
                raise RuntimeError(f"Duplicate paired candidate position: {sample} {local_position}")
            expected_positions.add(local_position)
            if local_position not in outputs:
                raise RuntimeError(f"Candidate missing paired output: {key}")
            tasks.append({"candidate": dict(payload), **outputs[local_position]})
        if set(outputs) != expected_positions:
            raise RuntimeError(
                f"{sample} paired output/candidate set mismatch: "
                f"missing={sorted(expected_positions - set(outputs))} "
                f"extra={sorted(set(outputs) - expected_positions)}"
            )
        if sample_receipt.get("validation", {}).get("exact_artifact_counts") is not True:
            raise RuntimeError(f"{sample} paired artifact reconciliation is not exact")
    return tasks, receipt, receipt_path


def primary_assignment_key(row: dict[str, Any]) -> CandidateKey:
    posthoc = row.get("posthoc") or {}
    ref = row.get("ref", posthoc.get("ref"))
    alt = row.get("alt", posthoc.get("alt"))
    if ref is None or alt is None:
        raise RuntimeError("Primary assignment lacks REF/ALT identity")
    return (
        str(row["sample"]),
        str(row["chrom"]),
        int(row["pos"]),
        str(ref).upper(),
        str(alt).upper(),
    )


def load_primary_assignments(
    path: Path, expected_keys: set[CandidateKey]
) -> dict[CandidateKey, dict[str, Any]]:
    expected_by_position = {(key[0], key[1], key[2]): key for key in expected_keys}
    assignments: dict[CandidateKey, dict[str, Any]] = {}
    with open_text(require_nonempty(path, "primary stable assignment JSONL")) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            local_position = (str(row["sample"]), str(row["chrom"]), int(row["pos"]))
            if local_position not in expected_by_position:
                continue
            key = primary_assignment_key(row)
            if key != expected_by_position[local_position]:
                raise RuntimeError(
                    f"Primary assignment REF/ALT mismatch at line {line_number}: "
                    f"expected={expected_by_position[local_position]} observed={key}"
                )
            if key in assignments:
                raise RuntimeError(f"Duplicate primary assignment candidate: {key}")
            core_reads = row.get("core_reads")
            if not isinstance(core_reads, list) or not core_reads:
                raise RuntimeError(f"Primary assignment has no core_reads: {key}")
            assignments[key] = row
    missing = expected_keys - set(assignments)
    if missing:
        raise RuntimeError(
            f"Missing primary stable assignments for {len(missing)} candidates: "
            f"{sorted(missing)[:5]}"
        )
    return assignments


def load_reads(path: Path) -> dict[str, dict[str, str]]:
    reads: dict[str, dict[str, str]] = {}
    with require_nonempty(path, "reads.tsv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "read_id",
            "read_name",
            "chr",
            "start",
            "end",
            "mapq",
            "alt_support",
            "is_tumor",
            "strand",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"{path} missing reads columns: {sorted(missing)}")
        for row in reader:
            read_id = row["read_id"]
            if read_id in reads:
                raise RuntimeError(f"Duplicate paired reads.tsv read_id: {path} {read_id}")
            reads[read_id] = row
    return reads


def load_matrix(path: Path, square: bool) -> tuple[list[str], np.ndarray]:
    identifiers, matrix = F.load_matrix(require_nonempty(path, "matrix"))
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError(f"Duplicate matrix read_id: {path}")
    if matrix.shape[0] != len(identifiers):
        raise RuntimeError(f"Matrix row count mismatch: {path} {matrix.shape}")
    if square and matrix.shape != (len(identifiers), len(identifiers)):
        raise RuntimeError(f"Distance matrix is not square: {path} {matrix.shape}")
    return identifiers, matrix


def analyze_subset(
    ids: list[str],
    distance_ids: list[str],
    distance: np.ndarray,
    methylation_ids: list[str],
    methylation: np.ndarray,
    seed: int,
) -> dict[str, Any]:
    if len(ids) < 2 * F.MIN_SIZE:
        return {
            "status": "insufficient_reads",
            "evaluable": False,
            "n_raw": len(ids),
            "n_matrix": 0,
            "n_after_peel": 0,
            "stable_multigroup": False,
            "group_sizes": {},
        }
    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    usable = [
        read_id
        for read_id in ids
        if read_id in distance_index and read_id in methylation_index
    ]
    if len(usable) < 2 * F.MIN_SIZE:
        return {
            "status": "insufficient_matrix_reads",
            "evaluable": False,
            "n_raw": len(ids),
            "n_matrix": len(usable),
            "n_after_peel": 0,
            "stable_multigroup": False,
            "group_sizes": {},
        }
    distance_rows = [distance_index[read_id] for read_id in usable]
    sub_distance = distance[np.ix_(distance_rows, distance_rows)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [usable[index] for index in kept]
    if len(kept_ids) < 2 * F.MIN_SIZE:
        return {
            "status": "incomplete_distance_below_minimum",
            "evaluable": False,
            "n_raw": len(ids),
            "n_matrix": len(usable),
            "n_after_peel": len(kept_ids),
            "stable_multigroup": False,
            "group_sizes": {},
        }
    sub_distance = sub_distance[np.ix_(kept, kept)]
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]
    result = F.analyze_phylo(
        sub_distance,
        sub_methylation,
        base_seed=seed,
        seeds=10,
        rnull=F.RNULL,
        null_mode="column",
        empirical_alpha=None,
    )
    group_sizes = Counter(
        label for label in result["coarse_labels"] if label not in {"other", "outlier"}
    )
    stable = bool(result["coarse_ng"] >= 2 and not result["unstable"])
    return {
        "status": "evaluable",
        "evaluable": True,
        "n_raw": len(ids),
        "n_matrix": len(usable),
        "n_after_peel": len(kept_ids),
        "coarse_ng": result["coarse_ng"],
        "modal_fraction": result["modal_fraction"],
        "unstable": result["unstable"],
        "stable_multigroup": stable,
        "group_sizes": dict(group_sizes),
    }


def normalize_allele(value: Any) -> str:
    normalized = str(value).strip().upper()
    return normalized if normalized in {"ALT", "REF"} else "UNKNOWN"


def focal_allele_counts(reads: Iterable[dict[str, Any]]) -> dict[str, int]:
    counts = Counter(normalize_allele(row.get("alt_support")) for row in reads)
    return {allele: int(counts.get(allele, 0)) for allele in ("ALT", "REF", "UNKNOWN")}


def partition_read_ids(reads: dict[str, dict[str, Any]]) -> tuple[list[str], list[str]]:
    """Return normal focal-REF IDs and tumor focal-REF IDs without consulting HP."""
    normal_ref_ids: list[str] = []
    tumor_ref_ids: list[str] = []
    for read_id, row in reads.items():
        allele = normalize_allele(row["alt_support"])
        if F.is_tumor(row["is_tumor"]):
            if allele == "REF":
                tumor_ref_ids.append(read_id)
        elif allele == "REF":
            normal_ref_ids.append(read_id)
    return normal_ref_ids, tumor_ref_ids


def normal_focal_callability(counts: dict[str, int]) -> bool:
    called = counts["ALT"] + counts["REF"]
    return bool(called >= NORMAL_MIN_CALLED_READS and counts["REF"] >= NORMAL_MIN_REF_READS)


def normal_not_evaluable_reason(
    counts: dict[str, int], normal_result: dict[str, Any]
) -> str | None:
    called = counts["ALT"] + counts["REF"]
    if called < NORMAL_MIN_CALLED_READS:
        return f"NORMAL_CALLED_READS_LT_{NORMAL_MIN_CALLED_READS}"
    if counts["REF"] < NORMAL_MIN_REF_READS:
        return f"NORMAL_REF_READS_LT_{NORMAL_MIN_REF_READS}"
    if not normal_result["evaluable"]:
        return f"NORMAL_REF_METHYL_{str(normal_result['status']).upper()}"
    return None


def alignment_identity(row: dict[str, Any]) -> AlignmentIdentity:
    chrom = row["chrom"] if "chrom" in row else row.get("chr")
    if chrom is None:
        raise RuntimeError(f"Alignment identity lacks chrom/chr: {row}")
    try:
        return (
            str(row["read_name"]),
            str(chrom),
            int(row["start"]),
            int(row["end"]),
            int(row["mapq"]),
            str(row["strand"]),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"Malformed alignment identity: {row}") from error


def join_primary_group_assignments(
    core_reads: list[dict[str, Any]], paired_reads: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    paired_index: dict[AlignmentIdentity, list[tuple[str, dict[str, Any]]]] = defaultdict(list)
    for read_id, row in paired_reads.items():
        if F.is_tumor(row["is_tumor"]):
            paired_index[alignment_identity(row)].append((read_id, row))

    seen_primary: set[AlignmentIdentity] = set()
    missing: list[AlignmentIdentity] = []
    collisions: list[AlignmentIdentity] = []
    allele_mismatches: list[AlignmentIdentity] = []
    joined: list[tuple[AlignmentIdentity, str, str]] = []
    for core in core_reads:
        identity = alignment_identity(core)
        if identity in seen_primary:
            raise RuntimeError(f"Primary core identity collision: {identity}")
        seen_primary.add(identity)
        matches = paired_index.get(identity, [])
        if not matches:
            missing.append(identity)
            continue
        if len(matches) != 1:
            collisions.append(identity)
            continue
        read_id, paired = matches[0]
        if normalize_allele(paired["alt_support"]) != "ALT":
            allele_mismatches.append(identity)
            continue
        label = str(core.get("label", ""))
        if not label:
            raise RuntimeError(f"Primary core read lacks group label: {identity}")
        joined.append((identity, read_id, label))

    if collisions:
        raise RuntimeError(
            f"Paired tumor identity collision for {len(collisions)} primary core reads: "
            f"{collisions[:3]}"
        )
    if missing:
        raise RuntimeError(
            f"Missing paired tumor identities for {len(missing)} primary core reads: {missing[:3]}"
        )
    if allele_mismatches:
        raise RuntimeError(
            "Primary core identities are not ALT in paired tumor output: "
            f"{allele_mismatches[:3]}"
        )
    if len(joined) != len(core_reads):
        raise RuntimeError("Primary group assignment join coverage is not exact")

    digest = hashlib.sha256()
    for identity, read_id, label in sorted(joined):
        digest.update(
            json.dumps([*identity, read_id, label], separators=(",", ":")).encode("utf-8")
        )
        digest.update(b"\n")
    return {
        "primary_core_expected": len(core_reads),
        "primary_core_joined": len(joined),
        "primary_group_assignment_coverage": len(joined) / len(core_reads),
        "primary_identity_collision_count": 0,
        "primary_identity_missing_count": 0,
        "primary_allele_mismatch_count": 0,
        "primary_group_sizes": dict(Counter(label for _, _, label in joined)),
        "primary_join_sha256": digest.hexdigest(),
    }


def prefix_result(prefix: str, result: dict[str, Any], output: dict[str, Any]) -> None:
    for name, value in result.items():
        output[f"{prefix}_{name}"] = value


def analyze_site_from_data(
    candidate: dict[str, Any],
    reads: dict[str, dict[str, Any]],
    distance_ids: list[str],
    distance: np.ndarray,
    methylation_ids: list[str],
    methylation: np.ndarray,
    primary_assignment: dict[str, Any],
) -> dict[str, Any]:
    """Analyze one site; normal methylation receives focal REF reads only."""
    normal_ref_ids, tumor_ref_ids = partition_read_ids(reads)
    seed = F.stable_seed(candidate["sample"], candidate["chrom"], int(candidate["pos"]))
    normal_result = analyze_subset(
        normal_ref_ids, distance_ids, distance, methylation_ids, methylation, seed + 100_000
    )
    tumor_ref_result = analyze_subset(
        tumor_ref_ids, distance_ids, distance, methylation_ids, methylation, seed + 200_000
    )
    normal_rows = [row for row in reads.values() if not F.is_tumor(row["is_tumor"])]
    tumor_rows = [row for row in reads.values() if F.is_tumor(row["is_tumor"])]
    normal_counts = focal_allele_counts(normal_rows)
    tumor_counts = focal_allele_counts(tumor_rows)
    normal_called = normal_counts["ALT"] + normal_counts["REF"]
    callability_gate = normal_focal_callability(normal_counts)
    not_evaluable_reason = normal_not_evaluable_reason(normal_counts, normal_result)
    normal_control_evaluable = bool(callability_gate and normal_result["evaluable"])
    normal_nonreplication_gate = (
        not bool(normal_result["stable_multigroup"])
        if normal_control_evaluable
        else None
    )
    join = join_primary_group_assignments(primary_assignment["core_reads"], reads)

    output: dict[str, Any] = {
        "sample": str(candidate["sample"]),
        "chrom": str(candidate["chrom"]),
        "pos": int(candidate["pos"]),
        "ref": str(candidate["ref"]).upper(),
        "alt": str(candidate["alt"]).upper(),
        "n_reads_total": len(reads),
        "n_tumor_reads": len(tumor_rows),
        "n_normal_reads": len(normal_rows),
        "tumor_focal_allele_counts": tumor_counts,
        "normal_focal_allele_counts": normal_counts,
        "tumor_alt_reads": tumor_counts["ALT"],
        "tumor_ref_reads": tumor_counts["REF"],
        "tumor_unknown_reads": tumor_counts["UNKNOWN"],
        "normal_called_reads": normal_called,
        "normal_alt_reads": normal_counts["ALT"],
        "normal_ref_reads": normal_counts["REF"],
        "normal_unknown_reads": normal_counts["UNKNOWN"],
        "normal_focal_callability_min_called_reads": NORMAL_MIN_CALLED_READS,
        "normal_focal_callability_min_ref_reads": NORMAL_MIN_REF_READS,
        "normal_focal_callability_gate": callability_gate,
        "normal_genetic_alt_support_count": normal_counts["ALT"],
        "normal_genetic_alt_support_present": normal_counts["ALT"] > 0,
        "normal_genetic_alt_absence_gate": normal_counts["ALT"] == 0,
        "normal_genetic_alt_fraction_called": (
            normal_counts["ALT"] / normal_called if normal_called else None
        ),
        "normal_methyl_background_population": NORMAL_METHYL_BACKGROUND_POPULATION,
        "normal_control_evaluable": normal_control_evaluable,
        "normal_control_not_evaluable_reason": not_evaluable_reason,
        "normal_ref_methyl_nonreplication_gate": normal_nonreplication_gate,
        "normal_control_status": (
            "NOT_EVALUABLE"
            if not normal_control_evaluable
            else (
                "REF_METHYL_PARTITION_REPRODUCED"
                if normal_result["stable_multigroup"]
                else "REF_METHYL_PARTITION_NOT_REPRODUCED"
            )
        ),
        "normal_hp_used": False,
        "normal_hp_policy": NORMAL_HP_POLICY,
        "guardrail": GUARDRAIL,
        **join,
    }
    prefix_result("normal_ref_methyl", normal_result, output)
    prefix_result("tumor_ref", tumor_ref_result, output)
    return output


def process_task(task: dict[str, Any], assignment: dict[str, Any]) -> dict[str, Any]:
    reads = load_reads(Path(task["reads_path"]))
    distance_ids, distance = load_matrix(Path(task["distance_path"]), square=True)
    methylation_ids, methylation = load_matrix(Path(task["methylation_path"]), square=False)
    row = analyze_site_from_data(
        task["candidate"],
        reads,
        distance_ids,
        distance,
        methylation_ids,
        methylation,
        assignment,
    )
    row.update(
        {
            "region_dir": task["region_dir"],
            "reads_path": task["reads_path"],
            "distance_path": task["distance_path"],
            "methylation_path": task["methylation_path"],
        }
    )
    return row


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    normal_evaluable = [row for row in rows if row["normal_control_evaluable"]]
    tumor_ref_evaluable = [row for row in rows if row["tumor_ref_evaluable"]]
    return {
        "n_candidates": len(rows),
        "n_normal_focal_callability_gate": sum(
            row["normal_focal_callability_gate"] for row in rows
        ),
        "n_normal_control_evaluable": len(normal_evaluable),
        "n_normal_control_not_evaluable": len(rows) - len(normal_evaluable),
        "normal_control_not_evaluable_reason_counts": dict(
            Counter(
                row["normal_control_not_evaluable_reason"]
                for row in rows
                if row["normal_control_not_evaluable_reason"] is not None
            )
        ),
        "n_normal_ref_methyl_stable_multigroup": sum(
            row["normal_ref_methyl_stable_multigroup"] for row in normal_evaluable
        ),
        "n_normal_ref_methyl_nonreplicating": sum(
            row["normal_ref_methyl_nonreplication_gate"] is True
            for row in normal_evaluable
        ),
        "normal_ref_methyl_stable_fraction_evaluable": (
            sum(row["normal_ref_methyl_stable_multigroup"] for row in normal_evaluable)
            / len(normal_evaluable)
            if normal_evaluable
            else None
        ),
        "n_normal_genetic_alt_absence_gate": sum(
            row["normal_genetic_alt_absence_gate"] for row in rows
        ),
        "n_normal_genetic_alt_support_sites": sum(
            row["normal_genetic_alt_support_present"] for row in rows
        ),
        "normal_genetic_alt_support_reads": sum(
            row["normal_genetic_alt_support_count"] for row in rows
        ),
        "n_tumor_ref_evaluable": len(tumor_ref_evaluable),
        "n_tumor_ref_stable_multigroup": sum(
            row["tumor_ref_stable_multigroup"] for row in tumor_ref_evaluable
        ),
        "tumor_ref_stable_fraction_evaluable": (
            sum(row["tumor_ref_stable_multigroup"] for row in tumor_ref_evaluable)
            / len(tumor_ref_evaluable)
            if tumor_ref_evaluable
            else None
        ),
        "primary_core_expected": sum(row["primary_core_expected"] for row in rows),
        "primary_core_joined": sum(row["primary_core_joined"] for row in rows),
        "all_primary_group_assignments_exact": all(
            row["primary_group_assignment_coverage"] == 1.0
            and row["primary_identity_collision_count"] == 0
            and row["primary_identity_missing_count"] == 0
            for row in rows
        ),
    }


def tsv_value(value: Any) -> Any:
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    if value is None:
        return ""
    return value


def write_site_table(path: Path, rows: list[dict[str, Any]]) -> list[str]:
    identity_fields = ["sample", "chrom", "pos", "ref", "alt"]
    other_fields = sorted({name for row in rows for name in row} - set(identity_fields))
    fields = identity_fields + other_fields
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({name: tsv_value(row.get(name)) for name in fields})
    return fields


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paired-output-root", type=Path, required=True)
    parser.add_argument("--primary-assignments", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=TOPIC_ROOT / "results" / "matched_normal_candidate_controls_analysis_v2",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    started_at = now_utc()
    tasks, paired_receipt, paired_receipt_path = load_paired_tasks(args.paired_output_root)
    expected_keys = {candidate_key(task["candidate"]) for task in tasks}
    if len(expected_keys) != len(tasks):
        raise RuntimeError("Duplicate candidate identities in paired tasks")
    assignments = load_primary_assignments(args.primary_assignments, expected_keys)

    rows = [
        process_task(task, assignments[candidate_key(task["candidate"])]) for task in tasks
    ]
    rows.sort(
        key=lambda row: (row["sample"], row["chrom"], row["pos"], row["ref"], row["alt"])
    )
    if len(rows) != len(expected_keys):
        raise RuntimeError("Per-site analysis count mismatch")

    create_output_dir(args.output_dir)
    table_path = args.output_dir / "matched_normal_candidate_controls.tsv.gz"
    summary_path = args.output_dir / "matched_normal_candidate_controls_summary.json"
    receipt_path = args.output_dir / "run_receipt.json"
    fields = write_site_table(table_path, rows)
    pooled = summarize(rows)
    per_sample = {
        sample: summarize([row for row in rows if row["sample"] == sample])
        for sample in sorted({row["sample"] for row in rows})
    }
    summary = {
        "schema_name": "intersubmod.matched_normal_candidate_control_analysis",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "population": "G2/R1 candidate controls from paired tumor/matched-normal output",
        "method": {
            "normal_focal_counts": "all is_tumor=0 reads, reported as REF/ALT/UNKNOWN",
            "normal_methyl_background": NORMAL_METHYL_BACKGROUND_POPULATION,
            "normal_focal_callability_gate": {
                "minimum_called_reads": NORMAL_MIN_CALLED_READS,
                "minimum_ref_reads": NORMAL_MIN_REF_READS,
            },
            "normal_control_evaluable": (
                "normal_focal_callability_gate AND focal-REF-only methyl subset evaluable"
            ),
            "not_evaluable_is_negative_result": False,
            "tumor_ref": "is_tumor=1 and focal call=REF",
            "clustering": (
                f"focal_alt_cluster_lib MIN_SIZE={F.MIN_SIZE}; peel_complete; "
                f"analyze_phylo seeds=10 RNULL={F.RNULL} column null"
            ),
            "primary_join": (
                "exact (read_name,chrom,start,end,mapq,strand) to paired tumor ALT reads"
            ),
            "normal_hp_policy": NORMAL_HP_POLICY,
        },
        "pooled": pooled,
        "per_sample": per_sample,
        "guardrail": GUARDRAIL,
        "pass_semantics": "execution_and_identity_integrity_only_not_background_negativity",
        "pass": len(rows) == len(expected_keys) and pooled["all_primary_group_assignments_exact"],
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    receipt = {
        "schema_name": "intersubmod.matched_normal_candidate_control_analysis_run",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "started_at_utc": started_at,
        "finished_at_utc": now_utc(),
        "command": [sys.executable, str(Path(__file__).resolve()), *(argv or sys.argv[1:])],
        "inputs": {
            "paired_output_root": str(args.paired_output_root.resolve()),
            "paired_run_receipt": {
                **artifact(paired_receipt_path),
                "schema_name": paired_receipt.get("schema_name"),
            },
            "primary_assignments": artifact(args.primary_assignments),
            "analyzer_script": artifact(Path(__file__)),
        },
        "outputs": {
            "site_table": artifact(table_path),
            "summary": artifact(summary_path),
            "run_receipt": str(receipt_path.resolve()),
            "site_table_fields": fields,
        },
        "counts": pooled,
        "normal_methyl_background_population": NORMAL_METHYL_BACKGROUND_POPULATION,
        "normal_hp_policy": NORMAL_HP_POLICY,
        "identity_join_policy": "collision/missing/allele mismatch hard fail",
        "not_evaluable_is_negative_result": False,
        "guardrail": GUARDRAIL,
        "pass_semantics": summary["pass_semantics"],
        "pass": summary["pass"],
    }
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "input_paired_output_root": str(args.paired_output_root.resolve()),
                "input_primary_assignments": str(args.primary_assignments.resolve()),
                "site_table": str(table_path.resolve()),
                "summary": str(summary_path.resolve()),
                "run_receipt": str(receipt_path.resolve()),
                "counts": pooled,
                "pass": summary["pass"],
            },
            indent=2,
        )
    )
    return 0 if summary["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
