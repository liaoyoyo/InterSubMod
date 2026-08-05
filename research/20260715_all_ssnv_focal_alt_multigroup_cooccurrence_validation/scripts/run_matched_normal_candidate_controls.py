#!/usr/bin/env python3
"""Run paired tumor/matched-normal InterSubMod controls for selected candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, NamedTuple, Sequence

import pysam


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import release_source_authority as SOURCE_AUTHORITY  # noqa: E402


TOPIC_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TOPIC_ROOT.parents[1]
RESULT_ROOT = TOPIC_ROOT / "results"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
CANONICAL_PYTHON_CACHE_DIRNAME = (
    ".python_cache_m2v5_completion_v2_bound_bootstrap"
)
CANONICAL_PYTHON_CACHE_ROOT = WORKSPACE_ROOT / CANONICAL_PYTHON_CACHE_DIRNAME
CANONICAL_CANDIDATE_TABLE = (
    WORKSPACE_ROOT
    / "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
    / "methyl_ssnv_site_results.tsv.gz"
)
CANONICAL_MANIFEST = RESULT_ROOT / "all_ssnv_input_manifest.json"
CANONICAL_NORMAL_AUDIT = RESULT_ROOT / "matched_normal_methyl_tag_audit.v3_pre_candidate.json"
CANONICAL_BINARY = REPO_ROOT / "build" / "bin" / "inter_sub_mod"
CANONICAL_REFERENCE = Path("/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta")
CANONICAL_OUTPUT_ROOT = (
    WORKSPACE_ROOT / "matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
)
FORMAL_SELECTION_COLUMN = "multi_marker_molecular_haplotype_base_candidate"
PRODUCER_CODE_PATHS = {"matched_normal_runner": Path(__file__).resolve()}
DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
REQUIRED_CANDIDATE_COLUMNS = {"sample", "chrom", "pos", "ref", "alt"}
TRUE_VALUES = {"1", "true", "t", "yes", "y"}
FALSE_VALUES = {"0", "false", "f", "no", "n"}
REGION_ERROR_PATTERN = re.compile(r"\[ERROR\]\s+Region\s+.*\s+failed:")
ACCEPTED_REGION_STRATIFICATION_STATUSES = {"VALID", "INSUFFICIENT_REGIONS"}


class Candidate(NamedTuple):
    sample: str
    chrom: str
    pos: int
    ref: str
    alt: str

    def payload(self) -> dict[str, Any]:
        return self._asdict()


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def capture_source_identity() -> dict[str, dict[str, Any]]:
    return {role: source_artifact(path) for role, path in PRODUCER_CODE_PATHS.items()}


def capture_source_modes() -> dict[str, str]:
    return {
        role: oct(path.resolve(strict=True).stat().st_mode & 0o777)
        for role, path in PRODUCER_CODE_PATHS.items()
    }


def canonical_python_prefix() -> list[str]:
    return [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={CANONICAL_PYTHON_CACHE_ROOT}",
    ]


def canonical_task_b_command() -> list[str]:
    return [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        "--candidate-table",
        str(CANONICAL_CANDIDATE_TABLE.resolve()),
        "--selection-column",
        FORMAL_SELECTION_COLUMN,
        "--selection-value",
        "true",
        "--manifest",
        str(CANONICAL_MANIFEST.resolve()),
        "--normal-audit",
        str(CANONICAL_NORMAL_AUDIT.resolve()),
        "--binary",
        str(CANONICAL_BINARY.resolve()),
        "--reference",
        str(CANONICAL_REFERENCE.resolve()),
        "--output-root",
        str(CANONICAL_OUTPUT_ROOT.resolve()),
        "--threads-per-sample",
        "40",
    ]


def observed_process_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise RuntimeError("Matched-normal runner process command is unavailable or malformed")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def resolve_release_command(
    argv: Sequence[str] | None, source_authority: dict[str, Any]
) -> list[str]:
    command = [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        *(list(argv) if argv is not None else sys.argv[1:]),
    ]
    if source_authority.get("authority_id") == SOURCE_AUTHORITY.AUTHORITY_ID:
        expected = canonical_task_b_command()
        if argv is not None or command != expected or observed_process_command() != expected:
            raise RuntimeError(
                "Formal matched-normal runner is direct-CLI canonical-process only"
            )
    return command


def producer_source_lock(
    identity_before: dict[str, dict[str, Any]],
    modes_before: dict[str, str],
) -> dict[str, Any]:
    identity_after = capture_source_identity()
    modes_after = capture_source_modes()
    if identity_after != identity_before or modes_after != modes_before:
        raise RuntimeError("Matched-normal runner source identity changed during execution")
    if set(modes_after.values()) != {"0o444"}:
        raise RuntimeError(f"Matched-normal runner source is not mode 0444: {modes_after}")
    return {
        "source_identity_before": identity_before,
        "source_identity_after_compute": identity_after,
        "source_modes_before": modes_before,
        "source_modes_after_compute": modes_after,
        "all_sources_read_only_and_unchanged": True,
    }


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def parse_boolean(value: Any) -> bool | None:
    normalized = str(value).strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    return None


def selection_matches(observed: Any, expected: str | None) -> bool:
    observed_boolean = parse_boolean(observed)
    if expected is None:
        return observed_boolean is True
    expected_boolean = parse_boolean(expected)
    if expected_boolean is not None:
        return observed_boolean is not None and observed_boolean == expected_boolean
    return str(observed).strip() == expected.strip()


def read_candidate_table(
    path: Path,
    selection_column: str | None = None,
    selection_value: str | None = None,
    allowed_samples: set[str] | None = None,
) -> list[Candidate]:
    """Parse and validate selected candidate rows from TSV or TSV.gz."""
    if selection_value is not None and selection_column is None:
        raise ValueError("--selection-value requires --selection-column")
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

    candidates: list[Candidate] = []
    seen: set[Candidate] = set()
    seen_positions: dict[tuple[str, str, int], Candidate] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = REQUIRED_CANDIDATE_COLUMNS - fields
        if missing:
            raise ValueError(f"Candidate table missing columns: {sorted(missing)}")
        if selection_column is not None and selection_column not in fields:
            raise ValueError(f"Candidate table missing selection column: {selection_column}")

        for line_number, row in enumerate(reader, 2):
            if selection_column is not None and not selection_matches(
                row.get(selection_column, ""), selection_value
            ):
                continue
            try:
                sample = (row.get("sample") or "").strip()
                chrom = (row.get("chrom") or "").strip()
                pos = int((row.get("pos") or "").strip())
                ref = (row.get("ref") or "").strip().upper()
                alt = (row.get("alt") or "").strip().upper()
            except (TypeError, ValueError) as error:
                raise ValueError(f"Invalid candidate identity at {path}:{line_number}") from error
            if not sample or not chrom or pos <= 0 or not ref or not alt:
                raise ValueError(f"Invalid candidate identity at {path}:{line_number}")
            if allowed_samples is not None and sample not in allowed_samples:
                raise ValueError(f"Unknown candidate sample at {path}:{line_number}: {sample}")

            candidate = Candidate(sample, chrom, pos, ref, alt)
            if candidate in seen:
                raise ValueError(f"Duplicate candidate row at {path}:{line_number}: {candidate}")
            position_key = (sample, chrom, pos)
            if position_key in seen_positions:
                raise ValueError(
                    f"Duplicate candidate position at {path}:{line_number}: "
                    f"{seen_positions[position_key]} versus {candidate}"
                )
            seen.add(candidate)
            seen_positions[position_key] = candidate
            candidates.append(candidate)

    if not candidates and selection_column is None:
        raise ValueError("Candidate selection produced zero rows")
    return candidates


def require_nonempty(path: Path, label: str) -> Path:
    if not path.exists() or not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")
    return path


def validate_manifest(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    if manifest.get("schema_name") != "intersubmod.all_ssnv_focal_alt_input_manifest":
        raise RuntimeError("Unexpected all-sSNV input manifest schema")
    if manifest.get("pass") is not True:
        raise RuntimeError("All-sSNV input manifest is not passing")
    entries = manifest.get("samples")
    if not isinstance(entries, list):
        raise RuntimeError("All-sSNV input manifest has no sample list")
    by_sample: dict[str, dict[str, Any]] = {}
    for entry in entries:
        sample = entry.get("sample")
        if sample in by_sample:
            raise RuntimeError(f"Duplicate manifest sample: {sample}")
        by_sample[sample] = entry
    if set(by_sample) != set(DATASETS) or len(by_sample) != 7:
        raise RuntimeError("All-sSNV manifest must contain the fixed 7/7 datasets")
    return by_sample


def validate_normal_audit(
    audit: dict[str, Any], expected_samples: set[str]
) -> dict[str, dict[str, Any]]:
    """Enforce the authoritative 7/7 eligibility gate without inferring BAM paths."""
    if audit.get("schema_name") != "intersubmod.matched_normal_methyl_tag_audit":
        raise RuntimeError("Unexpected matched-normal audit schema")
    if audit.get("pass") is not True:
        raise RuntimeError("Matched-normal audit is not passing")
    rows = audit.get("samples")
    if not isinstance(rows, list):
        raise RuntimeError("Matched-normal audit has no sample list")
    by_sample: dict[str, dict[str, Any]] = {}
    for row in rows:
        sample = row.get("sample")
        if sample in by_sample:
            raise RuntimeError(f"Duplicate matched-normal audit sample: {sample}")
        by_sample[sample] = row
    totals = audit.get("totals", {})
    if (
        expected_samples != set(DATASETS)
        or set(by_sample) != expected_samples
        or len(by_sample) != 7
        or int(totals.get("n_samples", -1)) != 7
        or int(totals.get("n_normal_control_eligible", -1)) != 7
    ):
        raise RuntimeError("Matched-normal audit failed the fixed 7/7 coverage gate")
    ineligible = sorted(
        sample for sample, row in by_sample.items() if row.get("normal_control_eligible") is not True
    )
    if ineligible:
        raise RuntimeError(f"Matched-normal audit has ineligible samples: {ineligible}")
    for sample, row in by_sample.items():
        if not row.get("bam") or not row.get("index"):
            raise RuntimeError(f"Matched-normal audit lacks BAM/index path for {sample}")
    return by_sample


def record_identity(record: pysam.libcbcf.VariantRecord, sample: str) -> Candidate | None:
    if not record.alts or len(record.alts) != 1:
        return None
    return Candidate(sample, record.chrom, int(record.pos), record.ref.upper(), record.alts[0].upper())


def inspect_candidate_matches(source_vcf: Path, candidates: Iterable[Candidate]) -> None:
    """Fail closed on absent, allele-mismatched, duplicate, or non-PASS source records."""
    requested = list(candidates)
    if not requested:
        raise ValueError("No candidates supplied for source VCF inspection")
    samples = {candidate.sample for candidate in requested}
    if len(samples) != 1:
        raise ValueError("Source VCF inspection requires one sample at a time")
    sample = next(iter(samples))
    by_position = {(candidate.chrom, candidate.pos): candidate for candidate in requested}
    found: set[Candidate] = set()
    observed_alleles: dict[tuple[str, int], set[tuple[str, str]]] = defaultdict(set)

    with pysam.VariantFile(str(source_vcf)) as source:
        for record in source:
            position_key = (record.chrom, int(record.pos))
            if position_key not in by_position:
                continue
            for alt in record.alts or ():
                observed_alleles[position_key].add((record.ref.upper(), alt.upper()))
            identity = record_identity(record, sample)
            if identity != by_position[position_key]:
                continue
            if identity in found:
                raise RuntimeError(f"Duplicate source VCF candidate record: {identity}")
            if list(record.filter.keys()) != ["PASS"]:
                raise RuntimeError(f"Candidate is not PASS in source VCF: {identity}")
            found.add(identity)

    for candidate in requested:
        if candidate in found:
            continue
        position_key = (candidate.chrom, candidate.pos)
        if position_key in observed_alleles:
            alleles = sorted(observed_alleles[position_key])
            raise RuntimeError(
                f"REF/ALT mismatch for candidate {candidate}; source alleles={alleles}"
            )
        raise RuntimeError(f"Candidate does not exist in source VCF: {candidate}")


def group_candidates(candidates: Iterable[Candidate]) -> dict[str, list[Candidate]]:
    grouped: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in candidates:
        grouped[candidate.sample].append(candidate)
    return dict(grouped)


def validate_candidates_against_sources(
    grouped: dict[str, list[Candidate]], manifest_by_sample: dict[str, dict[str, Any]]
) -> None:
    for sample, candidates in grouped.items():
        artifact = manifest_by_sample[sample]["all_ssnv_vcf"]
        source_vcf = require_nonempty(Path(artifact["path"]), f"{sample} source VCF")
        expected_hash = artifact.get("sha256")
        if not expected_hash or sha256(source_vcf) != expected_hash:
            raise RuntimeError(f"{sample} frozen LongPhase-S PASS VCF hash drift")
        inspect_candidate_matches(source_vcf, candidates)


def write_candidate_vcf(
    source_vcf: Path, output_vcf: Path, candidates: Iterable[Candidate]
) -> dict[str, Any]:
    requested = list(candidates)
    inspect_candidate_matches(source_vcf, requested)
    if output_vcf.exists() or Path(str(output_vcf) + ".tbi").exists():
        raise FileExistsError(f"Refusing to overwrite candidate VCF: {output_vcf}")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    requested_set = set(requested)
    sample = requested[0].sample
    written: set[Candidate] = set()
    with pysam.VariantFile(str(source_vcf)) as source, pysam.VariantFile(
        str(output_vcf), "wz", header=source.header
    ) as destination:
        for record in source:
            identity = record_identity(record, sample)
            if identity in requested_set:
                if identity in written:
                    raise RuntimeError(f"Duplicate candidate during VCF extraction: {identity}")
                destination.write(record)
                written.add(identity)
    if written != requested_set:
        raise RuntimeError("Candidate VCF extraction key set drift")
    pysam.tabix_index(str(output_vcf), preset="vcf", force=False)
    index = require_nonempty(Path(str(output_vcf) + ".tbi"), "candidate VCF index")
    return {
        "path": str(output_vcf.resolve()),
        "sha256": sha256(output_vcf),
        "index": str(index.resolve()),
        "index_sha256": sha256(index),
        "site_count": len(written),
    }


def create_output_root(path: Path) -> Path:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output root: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.mkdir(parents=False, exist_ok=False)
    return path


def parse_output_site(reads_path: Path) -> tuple[str, int]:
    site_name = reads_path.parents[2].name
    if "_" not in site_name:
        raise ValueError(f"Cannot parse InterSubMod site directory: {reads_path}")
    chrom, position = site_name.rsplit("_", 1)
    try:
        return chrom, int(position)
    except ValueError as error:
        raise ValueError(f"Cannot parse InterSubMod site directory: {reads_path}") from error


def csv_data_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in csv.reader(handle)) - 1)


def file_set_sha256(paths: Iterable[Path], base: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda value: str(value.relative_to(base))):
        digest.update(str(path.relative_to(base)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def collect_region_errors(paths: Iterable[Path]) -> list[str]:
    errors: list[str] = []
    for path in paths:
        if not path.exists():
            continue
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if REGION_ERROR_PATTERN.search(line):
                errors.append(line)
    return errors


def validate_output(
    output_dir: Path,
    candidates: Iterable[Candidate],
    log_paths: Iterable[Path],
) -> dict[str, Any]:
    requested_positions = {(candidate.chrom, candidate.pos) for candidate in candidates}
    observed: dict[tuple[str, int], dict[str, Path]] = {}
    artifact_paths: list[Path] = []
    if output_dir.exists():
        for reads_path in output_dir.rglob("reads.tsv"):
            if reads_path.parent.name != "reads":
                continue
            key = parse_output_site(reads_path)
            if key in observed:
                raise RuntimeError(f"Duplicate paired output site: {key}")
            region_dir = reads_path.parents[1]
            methylation_path = region_dir / "methylation" / "methylation.csv"
            bernoulli_path = region_dir / "distance" / "BERNOULLI" / "matrix.csv"
            observed[key] = {
                "reads": reads_path,
                "methylation": methylation_path,
                "bernoulli": bernoulli_path,
            }

    missing_artifacts: list[str] = []
    for key, paths in observed.items():
        for kind, path in paths.items():
            if not path.exists() or path.stat().st_size == 0:
                missing_artifacts.append(f"{key[0]}:{key[1]}:{kind}")
            else:
                artifact_paths.append(path)
    region_errors = collect_region_errors(log_paths)
    observed_positions = set(observed)
    status_path = output_dir / "region_stratification_status.tsv"
    status_rows: list[dict[str, str]] = []
    if status_path.exists() and status_path.stat().st_size > 0:
        with status_path.open(newline="", encoding="utf-8") as handle:
            status_rows = list(csv.DictReader(handle, delimiter="\t"))
    summary_path = output_dir / "significance_summary.csv"
    statistics_path = output_dir / "significance_statistics.txt"
    summary_rows = csv_data_rows(summary_path)
    region_status = status_rows[0].get("status") if len(status_rows) == 1 else None
    required_run_artifacts = (summary_path, statistics_path, status_path)
    missing_run_artifacts = [
        path.name for path in required_run_artifacts if not path.exists() or path.stat().st_size == 0
    ]
    for path in required_run_artifacts:
        if path.exists() and path.stat().st_size > 0:
            artifact_paths.append(path)
    exact = (
        observed_positions == requested_positions
        and not missing_artifacts
        and len(artifact_paths)
        == (3 * len(requested_positions) + len(required_run_artifacts))
        and not missing_run_artifacts
        and summary_rows == len(requested_positions)
        and region_status in ACCEPTED_REGION_STRATIFICATION_STATUSES
    )
    return {
        "requested_sites": len(requested_positions),
        "output_sites": len(observed_positions),
        "significance_summary_rows": summary_rows,
        "reads_files": sum("reads" in paths and paths["reads"].exists() for paths in observed.values()),
        "methylation_files": sum(
            "methylation" in paths and paths["methylation"].exists() for paths in observed.values()
        ),
        "bernoulli_matrix_files": sum(
            "bernoulli" in paths and paths["bernoulli"].exists() for paths in observed.values()
        ),
        "missing_sites": sorted(
            f"{chrom}:{pos}" for chrom, pos in requested_positions - observed_positions
        ),
        "extra_sites": sorted(
            f"{chrom}:{pos}" for chrom, pos in observed_positions - requested_positions
        ),
        "missing_artifacts": missing_artifacts,
        "missing_run_artifacts": missing_run_artifacts,
        "region_error_count": len(region_errors),
        "region_errors": region_errors,
        "region_stratification_status": region_status,
        "region_stratification_status_accepted": (
            region_status in ACCEPTED_REGION_STRATIFICATION_STATUSES
        ),
        "artifact_set_sha256": file_set_sha256(artifact_paths, output_dir) if artifact_paths else None,
        "exact_artifact_counts": exact,
        "pass": exact and not region_errors,
    }


def artifact_identity(path: Path, include_hash: bool = False) -> dict[str, Any]:
    require_nonempty(path, "input artifact")
    stat = path.stat()
    payload: dict[str, Any] = {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if include_hash:
        payload["sha256"] = sha256(path)
    return payload


def write_not_applicable_receipt(
    output_root: Path,
    *,
    candidate_table: Path,
    manifest_path: Path,
    normal_audit_path: Path,
    binary: Path,
    reference: Path,
    selection_column: str,
    selection_value: str | None,
    started_at_utc: str,
    command: list[str],
    source_authority: dict[str, Any],
    code: dict[str, dict[str, Any]],
    source_lock: dict[str, Any],
) -> Path:
    """Record a structural zero-candidate result without implying a normal run."""
    if not output_root.is_dir():
        raise FileNotFoundError(f"Missing output root for N/A receipt: {output_root}")
    if any(output_root.iterdir()):
        raise FileExistsError(f"N/A output root is not empty: {output_root}")

    receipt_path = output_root / "not_applicable_receipt.json"
    finished_at_utc = now_utc()
    receipt = {
        "schema_name": "intersubmod.matched_normal_candidate_control.not_applicable",
        "schema_version": "1.0.0",
        "created_at_utc": finished_at_utc,
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "task_type": "B_comprehensive_validation",
        "status": "NOT_APPLICABLE",
        "execution_status": "NOT_APPLICABLE",
        "reason": "ZERO_SELECTED_CANDIDATES",
        "n_selected_candidates": 0,
        "selection_column": selection_column,
        "selection_value": selection_value,
        "selection_contract": {
            "selection_column": selection_column,
            "selection_value": selection_value,
            "n_selected_candidates": 0,
        },
        "counts": {"n_selected_candidates": 0},
        "command": command,
        "source_authority": source_authority,
        "code": code,
        "source_lock": source_lock,
        "inputs": {
            "candidate_table": artifact_identity(candidate_table, include_hash=True),
            "all_ssnv_manifest": artifact_identity(manifest_path, include_hash=True),
            "normal_audit": artifact_identity(normal_audit_path, include_hash=True),
            "binary": artifact_identity(binary, include_hash=True),
            "reference": artifact_identity(reference, include_hash=True),
            "runner_script": artifact_identity(Path(__file__).resolve(), include_hash=True),
        },
        "outputs": {
            "output_root": str(output_root.resolve()),
            "not_applicable_receipt": str(receipt_path.resolve()),
            "sample_outputs": [],
        },
        "sample_outputs_created": False,
        "cpp_executed": False,
        "normal_control_executed": False,
        "not_evaluable_is_negative": False,
        "interpretation": "NOT_APPLICABLE is not a negative matched-normal result.",
        "pass_semantics": (
            "receipt_integrity_only_not_normal_control_execution_or_negative_evidence"
        ),
        "pass": True,
    }
    with receipt_path.open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    receipt_path.chmod(0o444)
    return receipt_path


def run_sample(
    sample: str,
    candidates: list[Candidate],
    manifest_entry: dict[str, Any],
    normal_entry: dict[str, Any],
    candidate_vcf: dict[str, Any],
    sample_root: Path,
    binary: Path,
    reference: Path,
    threads: int,
    shared_provenance: dict[str, Any],
) -> dict[str, Any]:
    output_dir = sample_root / "output"
    execution_log = sample_root / "execution.log"
    region_error_log = sample_root / "region_errors.log"
    receipt_path = sample_root / "run_receipt.json"
    for path in (output_dir, execution_log, region_error_log, receipt_path):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite sample output: {path}")

    tumor_bam = require_nonempty(Path(manifest_entry["raw_alignment"]["path"]), "tumor BAM")
    tumor_index = require_nonempty(
        Path(manifest_entry["raw_alignment_index"]["path"]), "tumor BAM index"
    )
    normal_bam = require_nonempty(Path(normal_entry["bam"]), "matched-normal BAM")
    normal_index = require_nonempty(Path(normal_entry["index"]), "matched-normal BAM index")
    output_dir.mkdir(parents=False, exist_ok=False)
    command = [
        str(binary),
        "-t",
        str(tumor_bam),
        "-n",
        str(normal_bam),
        "-r",
        str(reference),
        "-v",
        candidate_vcf["path"],
        "-o",
        str(output_dir),
        "-w",
        "5000",
        "-j",
        str(threads),
        "--distance-metric",
        "BERNOULLI",
        "--min-common-coverage",
        "3",
        "--nan-distance-strategy",
        "SKIP",
        "--methyl-low",
        "0.2",
        "--methyl-high",
        "0.8",
        "--log-level",
        "info",
    ]
    started_at = now_utc()
    error: str | None = None
    try:
        with execution_log.open("x", encoding="utf-8") as log_handle:
            completed = subprocess.run(
                command,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        exit_code: int | None = completed.returncode
    except Exception as caught:  # pragma: no cover - OS-level launch failures are environment-specific.
        exit_code = None
        error = repr(caught)

    internal_log = output_dir / "inter_sub_mod.log"
    validation = validate_output(output_dir, candidates, (execution_log, internal_log))
    region_error_log.write_text(
        "\n".join(validation["region_errors"]) + ("\n" if validation["region_errors"] else ""),
        encoding="utf-8",
    )
    passed = exit_code == 0 and validation["pass"] and error is None
    receipt = {
        "schema_name": "intersubmod.matched_normal_candidate_sample_run",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "started_at_utc": started_at,
        "finished_at_utc": now_utc(),
        "sample": sample,
        "candidate_sites": [candidate.payload() for candidate in candidates],
        "site_counts": {
            "candidate_table": len(candidates),
            "candidate_vcf": candidate_vcf["site_count"],
            "paired_output": validation["output_sites"],
        },
        "command": command,
        "inputs": {
            **shared_provenance,
            "tumor_bam": artifact_identity(tumor_bam),
            "tumor_bam_index": artifact_identity(tumor_index),
            "normal_bam_from_audit": artifact_identity(normal_bam),
            "normal_bam_index_from_audit": artifact_identity(normal_index),
            "normal_control_eligible": normal_entry["normal_control_eligible"],
            "source_longphase_pass_vcf": manifest_entry["all_ssnv_vcf"],
            "candidate_vcf": candidate_vcf,
        },
        "outputs": {
            "sample_root": str(sample_root.resolve()),
            "output_dir": str(output_dir.resolve()),
            "execution_log": str(execution_log.resolve()),
            "region_error_log": str(region_error_log.resolve()),
            "run_receipt": str(receipt_path.resolve()),
        },
        "exit_code": exit_code,
        "error": error,
        "validation": validation,
        "pass": passed,
        "guardrails": {
            "normal_bam_source": "matched_normal_methyl_tag_audit receipt only",
            "normal_hp": "not used by runner or downstream normal-control grouping",
            "overwrite": "new output root and sample outputs only",
        },
    }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--selection-column")
    parser.add_argument("--selection-value")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=TOPIC_ROOT / "results" / "all_ssnv_input_manifest.json",
    )
    parser.add_argument(
        "--normal-audit",
        type=Path,
        default=CANONICAL_NORMAL_AUDIT,
    )
    parser.add_argument("--binary", type=Path, default=Path("build/bin/inter_sub_mod"))
    parser.add_argument(
        "--reference", type=Path, default=Path("/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta")
    )
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--threads-per-sample", type=int, default=2)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    started_at_utc = now_utc()
    source_authority_before = SOURCE_AUTHORITY.validate_release_source_authority(
        {"release_source_authority_validator", "matched_normal_runner"}
    )
    command = resolve_release_command(argv, source_authority_before)
    source_identity_before = capture_source_identity()
    source_modes_before = capture_source_modes()
    if (
        source_authority_before.get("authority_id") == SOURCE_AUTHORITY.AUTHORITY_ID
        and set(source_modes_before.values()) != {"0o444"}
    ):
        raise RuntimeError(
            f"Formal matched-normal runner source is not mode 0444: {source_modes_before}"
        )
    args = parse_args(argv)
    if args.threads_per_sample < 1:
        raise ValueError("--threads-per-sample must be positive")
    manifest_path = require_nonempty(args.manifest, "all-sSNV manifest")
    normal_audit_path = require_nonempty(args.normal_audit, "matched-normal audit")
    candidate_table = require_nonempty(args.candidate_table, "candidate table")
    binary = require_nonempty(args.binary.resolve(), "InterSubMod binary")
    reference = require_nonempty(args.reference.resolve(), "reference FASTA")
    require_nonempty(Path(str(reference) + ".fai"), "reference FASTA index")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    normal_audit = json.loads(normal_audit_path.read_text(encoding="utf-8"))
    manifest_by_sample = validate_manifest(manifest)
    candidates = read_candidate_table(
        candidate_table,
        args.selection_column,
        args.selection_value,
        allowed_samples=set(manifest_by_sample),
    )

    output_root = args.output_root or (Path(manifest["workspace"]) / "matched_normal_candidate_controls_v1")
    if not candidates:
        source_authority_after = SOURCE_AUTHORITY.validate_release_source_authority(
            {"release_source_authority_validator", "matched_normal_runner"}
        )
        if source_authority_after != source_authority_before:
            raise RuntimeError("Matched-normal source authority changed during execution")
        source_lock = producer_source_lock(source_identity_before, source_modes_before)
        create_output_root(output_root)
        receipt_path = write_not_applicable_receipt(
            output_root,
            candidate_table=candidate_table,
            manifest_path=manifest_path,
            normal_audit_path=normal_audit_path,
            binary=binary,
            reference=reference,
            selection_column=args.selection_column,
            selection_value=args.selection_value,
            started_at_utc=started_at_utc,
            command=command,
            source_authority=source_authority_before,
            code=source_identity_before,
            source_lock=source_lock,
        )
        print(
            json.dumps(
                {
                    "not_applicable_receipt": str(receipt_path),
                    "output_root": str(output_root),
                    "n_selected_candidates": 0,
                    "execution_status": "NOT_APPLICABLE",
                    "reason": "ZERO_SELECTED_CANDIDATES",
                    "pass": True,
                },
                indent=2,
            )
        )
        raise SystemExit(0)

    normal_by_sample = validate_normal_audit(normal_audit, set(manifest_by_sample))
    for sample, row in normal_by_sample.items():
        require_nonempty(Path(row["bam"]), f"{sample} matched-normal BAM")
        require_nonempty(Path(row["index"]), f"{sample} matched-normal BAM index")

    grouped = group_candidates(candidates)
    validate_candidates_against_sources(grouped, manifest_by_sample)

    create_output_root(output_root)
    shared_provenance = {
        "candidate_table": artifact_identity(candidate_table, include_hash=True),
        "all_ssnv_manifest": artifact_identity(manifest_path, include_hash=True),
        "normal_audit": {
            **artifact_identity(normal_audit_path, include_hash=True),
            "command": normal_audit.get("command"),
            "source_code": normal_audit.get("source_code"),
        },
        "binary": artifact_identity(binary, include_hash=True),
        "reference": artifact_identity(reference),
        "runner_script": artifact_identity(Path(__file__).resolve(), include_hash=True),
    }

    receipts: list[dict[str, Any]] = []
    for sample in DATASETS:
        if sample not in grouped:
            continue
        sample_root = output_root / sample
        sample_root.mkdir(parents=False, exist_ok=False)
        source_vcf = Path(manifest_by_sample[sample]["all_ssnv_vcf"]["path"])
        candidate_vcf = write_candidate_vcf(
            source_vcf,
            sample_root / f"{sample}.matched_normal_candidates.vcf.gz",
            grouped[sample],
        )
        receipt = run_sample(
            sample,
            grouped[sample],
            manifest_by_sample[sample],
            normal_by_sample[sample],
            candidate_vcf,
            sample_root,
            binary,
            reference,
            args.threads_per_sample,
            shared_provenance,
        )
        receipts.append(receipt)
        print(
            f"[{sample}] pass={receipt['pass']} sites={receipt['site_counts']['candidate_table']} "
            f"reads={receipt['validation']['reads_files']} "
            f"methylation={receipt['validation']['methylation_files']} "
            f"bernoulli={receipt['validation']['bernoulli_matrix_files']}",
            flush=True,
        )

    failures = [
        {"sample": receipt["sample"], "exit_code": receipt["exit_code"], "error": receipt["error"]}
        for receipt in receipts
        if not receipt["pass"]
    ]
    source_authority_after = SOURCE_AUTHORITY.validate_release_source_authority(
        {"release_source_authority_validator", "matched_normal_runner"}
    )
    if source_authority_after != source_authority_before:
        raise RuntimeError("Matched-normal source authority changed during execution")
    source_lock = producer_source_lock(source_identity_before, source_modes_before)
    run_receipt_path = output_root / "run_receipt.json"
    finished_at_utc = now_utc()
    run_receipt = {
        "schema_name": "intersubmod.matched_normal_candidate_control_run",
        "schema_version": "1.0.0",
        "created_at_utc": finished_at_utc,
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "command": command,
        "source_authority": source_authority_before,
        "code": source_identity_before,
        "source_lock": source_lock,
        "output_root": str(output_root.resolve()),
        "selection": {
            "column": args.selection_column,
            "value": args.selection_value,
            "selected_sites": len(candidates),
            "selected_samples": [sample for sample in DATASETS if sample in grouped],
        },
        "normal_eligibility_gate": {
            "required": "7/7",
            "eligible": 7,
            "audit": shared_provenance["normal_audit"],
            "pass": True,
        },
        "receipts": receipts,
        "totals": {
            "candidate_sites": sum(receipt["site_counts"]["candidate_table"] for receipt in receipts),
            "paired_output_sites": sum(receipt["site_counts"]["paired_output"] for receipt in receipts),
            "reads_files": sum(receipt["validation"]["reads_files"] for receipt in receipts),
            "methylation_files": sum(
                receipt["validation"]["methylation_files"] for receipt in receipts
            ),
            "bernoulli_matrix_files": sum(
                receipt["validation"]["bernoulli_matrix_files"] for receipt in receipts
            ),
            "region_errors": sum(receipt["validation"]["region_error_count"] for receipt in receipts),
        },
        "failures": failures,
        "pass": len(receipts) == len(grouped) and not failures and all(
            receipt["pass"] for receipt in receipts
        ),
    }
    run_receipt_path.write_text(
        json.dumps(run_receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    run_receipt_path.chmod(0o444)
    print(
        json.dumps(
            {
                "run_receipt": str(run_receipt_path),
                "output_root": str(output_root),
                "totals": run_receipt["totals"],
                "pass": run_receipt["pass"],
            },
            indent=2,
        )
    )
    raise SystemExit(0 if run_receipt["pass"] else 1)


if __name__ == "__main__":
    main()
