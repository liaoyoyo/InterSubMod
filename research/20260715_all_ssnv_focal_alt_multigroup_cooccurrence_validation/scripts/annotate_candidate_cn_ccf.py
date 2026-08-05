#!/usr/bin/env python3
"""Attach authority-locked SAVANA CN and fit-local PyClone-VI annotations to candidates.

This program is an annotation layer, not a C1/clone confirmation step.  It validates
the complete 20260712 authority chain before publishing a fixed-schema output.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import io
import json
import math
import os
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence, TextIO


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import release_source_authority as SOURCE_AUTHORITY  # noqa: E402


SCHEMA_NAME = "intersubmod.candidate_cn_ccf_annotation"
SCHEMA_VERSION = "1.0.0"
RECEIPT_SCHEMA_NAME = "intersubmod.candidate_cn_ccf_annotation_receipt"
RECEIPT_SCHEMA_VERSION = "1.0.0"
CLAIM_CEILING = "conditional_annotation_only_not_orthogonal_confirmation"
PASS_SEMANTICS = "execution_integrity_only_not_orthogonal_confirmation"
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))
AUTOSOME_SET = frozenset(AUTOSOMES)
BASES = frozenset("ACGT")
CANONICAL_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
CANONICAL_PYTHON_CACHE_DIRNAME = (
    ".python_cache_m2v5_completion_v2_bound_bootstrap"
)
CANONICAL_PYTHON_CACHE_ROOT = WORKSPACE_ROOT / CANONICAL_PYTHON_CACHE_DIRNAME
CANONICAL_INPUT = (
    WORKSPACE_ROOT
    / "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
    / "methyl_ssnv_site_results.tsv.gz"
)
CANONICAL_OUTPUT_DIR = (
    WORKSPACE_ROOT / "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
)
FORMAL_SELECTION_COLUMN = "multi_marker_molecular_haplotype_base_candidate"
PRODUCER_CODE_PATHS = {"cn_ccf_annotator": Path(__file__).resolve()}
BLOCKED_SAMPLES = frozenset(("COLO829", "HCC1937"))

CN_STATUS_ENUM = (
    "AVAILABLE_EXACT_SEGMENT",
    "SHARED_CN_SENSITIVITY",
    "BLOCKED_CN_MISFIT",
    "NO_CN_SEGMENT",
    "RECEIPT_FAIL",
)
PYCLONE_STATUS_ENUM = (
    "MATCHED_PRIMARY",
    "MATCHED_SENSITIVITY_ONLY",
    "NOT_IN_FIT_UNIVERSE",
    "BLOCKED_CN_MISFIT",
    "RECEIPT_FAIL",
)

OUTPUT_COLUMNS = (
    "sample",
    "chrom",
    "pos",
    "ref",
    "alt",
    "mutation_id",
    "callset_status",
    "cn_status",
    "cn_authority_sample",
    "cn_transfer_policy",
    "independent_cn",
    "cn_segment_id",
    "cn_segment_start",
    "cn_segment_end",
    "savana_total_cn_raw",
    "savana_minor_cn_raw",
    "savana_total_cn_discrete",
    "savana_major_cn_discrete",
    "savana_minor_cn_discrete",
    "purity_model_value",
    "pyclone_status",
    "pyclone_primary_bundle_id",
    "pyclone_primary_model_role",
    "pyclone_fit_local_cluster_id",
    "pyclone_vi_cellular_prevalence",
    "pyclone_vi_cellular_prevalence_std",
    "pyclone_vi_assignment_probability",
    "pyclone_sensitivity_status",
    "pyclone_sensitivity_bundle_id",
    "pyclone_sensitivity_fit_local_cluster_id",
    "pyclone_sensitivity_cellular_prevalence",
    "pyclone_sensitivity_cellular_prevalence_std",
    "pyclone_sensitivity_assignment_probability",
    "authority_config_sha256",
    "input_provenance_sha256",
    "analysis_summary_sha256",
    "cn_source_sha256",
    "pyclone_primary_metadata_sha256",
    "pyclone_primary_results_sha256",
    "pyclone_primary_status_sha256",
    "pyclone_sensitivity_metadata_sha256",
    "pyclone_sensitivity_results_sha256",
    "pyclone_sensitivity_status_sha256",
    "claim_ceiling",
)


class ContractError(RuntimeError):
    """Raised when an input or authority contract fails closed."""


@dataclass(frozen=True, order=True)
class VariantKey:
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def mutation_id(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


@dataclass(frozen=True, order=True)
class CandidateKey:
    sample: str
    variant: VariantKey


@dataclass(frozen=True)
class CandidateInput:
    rows_read_total: int
    selected: tuple[CandidateKey, ...]
    fieldnames: tuple[str, ...]


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    end: int
    segment_id: str
    total_raw: float
    minor_raw: float
    total_discrete: int
    major_discrete: int
    minor_discrete: int
    near_integer: bool


@dataclass(frozen=True)
class BundleRole:
    bundle_id: str
    samples: tuple[str, ...]
    near_integer_only: bool
    model_role: str


@dataclass(frozen=True)
class FitHit:
    metadata: Mapping[str, str]
    result: Mapping[str, str]


@dataclass(frozen=True)
class BundleEvidence:
    role: BundleRole
    hits: Mapping[CandidateKey, FitHit]
    metadata_path: Path
    metadata_sha256: str
    results_path: Path
    results_sha256: str
    status_path: Path
    status_sha256: str


@dataclass(frozen=True)
class AuthorityPaths:
    config: Path
    input_provenance: Path
    analysis_summary: Path

    @classmethod
    def canonical(cls) -> "AuthorityPaths":
        repo_root = Path(__file__).resolve().parents[3]
        root = repo_root / "research" / "20260712_vaf_ccf_subclone_crosssoftware_validation"
        return cls(
            config=root / "config" / "pyclone_validation_config.json",
            input_provenance=root / "data" / "pyclone_inputs" / "provenance.json",
            analysis_summary=root / "runs" / "pyclone_vi" / "analysis" / "analysis_summary.json",
        )


@dataclass(frozen=True)
class AuthorityContext:
    config: Mapping[str, object]
    hashes: Mapping[str, str]
    bundles: Mapping[str, BundleEvidence]
    cn_indexes: Mapping[str, "SegmentIndex"]
    cn_paths: Mapping[str, Path]
    purity_by_sample: Mapping[str, float]
    config_sha256: str
    provenance_sha256: str
    summary_sha256: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_artifact(path: Path) -> dict[str, object]:
    source = resolved(path)
    return {
        "path": str(source),
        "size_bytes": source.stat().st_size,
        "sha256": sha256_file(source),
    }


def capture_source_identity() -> dict[str, dict[str, object]]:
    return {role: source_artifact(path) for role, path in PRODUCER_CODE_PATHS.items()}


def capture_source_modes() -> dict[str, str]:
    return {
        role: oct(resolved(path).stat().st_mode & 0o777)
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
        "--input",
        str(CANONICAL_INPUT.resolve()),
        "--selection-column",
        FORMAL_SELECTION_COLUMN,
        "--selection-value",
        "true",
        "--output-dir",
        str(CANONICAL_OUTPUT_DIR.resolve()),
    ]


def observed_process_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise ContractError("CN/CCF producer process command is unavailable or malformed")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def resolve_release_command(
    argv: Sequence[str] | None, source_authority: Mapping[str, object]
) -> list[str]:
    command = [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        *(list(argv) if argv is not None else sys.argv[1:]),
    ]
    if source_authority.get("authority_id") == SOURCE_AUTHORITY.AUTHORITY_ID:
        expected = canonical_task_b_command()
        if argv is not None or command != expected or observed_process_command() != expected:
            raise ContractError("Formal CN/CCF producer is direct-CLI canonical-process only")
    return command


def resolved(path: Path | str) -> Path:
    return Path(path).expanduser().resolve()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def load_json(path: Path, label: str) -> Mapping[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"Cannot read {label} JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ContractError(f"{label} must be a JSON object: {path}")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def parse_int(text: object, label: str, minimum: int | None = None) -> int:
    try:
        value = int(str(text))
    except (TypeError, ValueError) as exc:
        raise ContractError(f"Invalid integer for {label}: {text!r}") from exc
    if minimum is not None and value < minimum:
        raise ContractError(f"{label} must be >= {minimum}: {value}")
    return value


def parse_float(text: object, label: str) -> float:
    try:
        value = float(str(text))
    except (TypeError, ValueError) as exc:
        raise ContractError(f"Invalid float for {label}: {text!r}") from exc
    if not math.isfinite(value):
        raise ContractError(f"Non-finite float for {label}: {text!r}")
    return value


def round_half_up(value: float) -> int:
    return int(math.floor(value + 0.5))


def validate_variant(chrom: str, pos_text: object, ref: str, alt: str, label: str) -> VariantKey:
    pos = parse_int(pos_text, f"{label}.pos", minimum=1)
    require(chrom in AUTOSOME_SET, f"{label}.chrom must be GRCh38 chr1-chr22: {chrom!r}")
    require(len(ref) == 1 and ref in BASES, f"{label}.ref must be one uppercase A/C/G/T base: {ref!r}")
    require(len(alt) == 1 and alt in BASES, f"{label}.alt must be one uppercase A/C/G/T base: {alt!r}")
    require(ref != alt, f"{label} REF and ALT must differ")
    return VariantKey(chrom, pos, ref, alt)


def read_candidates(
    path: Path,
    selection_column: str | None = None,
    selection_value: str | None = None,
) -> CandidateInput:
    require(path.is_file(), f"Candidate input does not exist: {path}")
    require((selection_column is None) == (selection_value is None),
            "--selection-column and --selection-value must be supplied together")
    required_fields = {"sample", "chrom", "pos", "ref", "alt"}
    selected: list[CandidateKey] = []
    seen: set[CandidateKey] = set()
    rows_read = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        missing = required_fields.difference(fields)
        require(not missing, f"Candidate input misses columns: {sorted(missing)}")
        if selection_column is not None:
            require(selection_column in fields, f"Selection column is missing: {selection_column}")
        for row_number, row in enumerate(reader, start=2):
            rows_read += 1
            sample = row["sample"]
            require(sample in CANONICAL_SAMPLES,
                    f"row {row_number}.sample is not canonical: {sample!r}")
            variant = validate_variant(
                row["chrom"], row["pos"], row["ref"], row["alt"], f"row {row_number}"
            )
            key = CandidateKey(sample, variant)
            require(key not in seen, f"Duplicate candidate key at row {row_number}: {key}")
            seen.add(key)
            if selection_column is None or row[selection_column] == selection_value:
                selected.append(key)
    return CandidateInput(rows_read, tuple(selected), fields)


class SegmentIndex:
    """Fail-closed 1-based inclusive SAVANA segment index."""

    REQUIRED_FIELDS = {
        "chromosome",
        "start",
        "end",
        "segment_id",
        "copyNumber",
        "minorAlleleCopyNumber",
    }

    def __init__(self, path: Path):
        self.path = path
        self._segments: dict[str, list[Segment]] = {}
        self._starts: dict[str, list[int]] = {}
        self._load()

    def _load(self) -> None:
        require(self.path.is_file(), f"CN authority file is missing: {self.path}")
        with open_text(self.path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = self.REQUIRED_FIELDS.difference(reader.fieldnames or ())
            require(not missing, f"CN authority {self.path} misses columns: {sorted(missing)}")
            for row_number, row in enumerate(reader, start=2):
                chrom = row["chromosome"]
                if chrom not in AUTOSOME_SET:
                    continue
                start = parse_int(row["start"], f"{self.path}:{row_number}.start", minimum=1)
                end = parse_int(row["end"], f"{self.path}:{row_number}.end", minimum=start)
                segment_id = row["segment_id"]
                require(bool(segment_id), f"Empty segment_id at {self.path}:{row_number}")
                # The authority builder excludes measured segments with missing CN;
                # annotation must reproduce that no-imputation policy exactly.
                if row["copyNumber"].strip() == "" or row["minorAlleleCopyNumber"].strip() == "":
                    continue
                total_raw = parse_float(row["copyNumber"], f"{self.path}:{row_number}.copyNumber")
                minor_raw = parse_float(
                    row["minorAlleleCopyNumber"],
                    f"{self.path}:{row_number}.minorAlleleCopyNumber",
                )
                total_discrete = round_half_up(total_raw)
                minor_nearest = round_half_up(minor_raw)
                minor_discrete = min(minor_nearest, total_discrete // 2)
                major_discrete = total_discrete - minor_discrete
                if not (
                    total_discrete >= 1
                    and minor_discrete >= 0
                    and major_discrete >= 1
                    and major_discrete >= minor_discrete
                ):
                    # Match the frozen builder: invalid measured CN is excluded,
                    # never repaired or replaced with diploid CN.
                    continue
                near_integer = (
                    abs(total_raw - total_discrete) <= 0.25
                    and abs(minor_raw - minor_nearest) <= 0.25
                )
                self._segments.setdefault(chrom, []).append(
                    Segment(
                        chrom,
                        start,
                        end,
                        segment_id,
                        total_raw,
                        minor_raw,
                        total_discrete,
                        major_discrete,
                        minor_discrete,
                        near_integer,
                    )
                )
        for chrom, segments in self._segments.items():
            segments.sort(key=lambda segment: (segment.start, segment.end, segment.segment_id))
            previous: Segment | None = None
            for segment in segments:
                if previous is not None and segment.start <= previous.end:
                    raise ContractError(
                        f"Overlapping CN segments are forbidden in {self.path}: "
                        f"{previous.segment_id} and {segment.segment_id}"
                    )
                previous = segment
            self._starts[chrom] = [segment.start for segment in segments]

    def find(self, variant: VariantKey) -> Segment | None:
        starts = self._starts.get(variant.chrom)
        if not starts:
            return None
        index = bisect.bisect_right(starts, variant.pos) - 1
        if index < 0:
            return None
        segment = self._segments[variant.chrom][index]
        return segment if segment.start <= variant.pos <= segment.end else None


def expected_bundle_roles(config: Mapping[str, object]) -> Mapping[str, BundleRole]:
    pair = config.get("hcc1395_pair")
    individuals = config.get("individual_samples")
    blocked = config.get("blocked_samples")
    require(isinstance(pair, dict), "Authority config misses hcc1395_pair")
    require(isinstance(individuals, dict), "Authority config misses individual_samples")
    require(isinstance(blocked, dict), "Authority config misses blocked_samples")
    pair_samples = tuple(pair.get("samples", ()))
    require(pair_samples == ("HCC1395", "HCC1395_DORADO"),
            f"Unexpected HCC pair samples: {pair_samples}")
    require(set(individuals) == {"H1437", "H2009", "HCC1954"},
            f"Unexpected individual sample set: {sorted(individuals)}")
    require(set(blocked) == BLOCKED_SAMPLES,
            f"Unexpected blocked sample set: {sorted(blocked)}")

    roles: dict[str, BundleRole] = {}
    for suffix, near in (("main", False), ("near_integer", True)):
        for prefix, role in (
            ("hcc1395_pair_primary_joint", "JOINT_DESCRIPTIVE_ONLY"),
            ("hcc1395_pair_pass_union_joint", "JOINT_DESCRIPTIVE_ONLY"),
        ):
            bundle_id = f"{prefix}_{suffix}"
            roles[bundle_id] = BundleRole(bundle_id, pair_samples, near, role)
        for sample in pair_samples:
            bundle_id = f"hcc1395_pair_primary_separate_{sample}_{suffix}"
            roles[bundle_id] = BundleRole(
                bundle_id,
                (sample,),
                near,
                "SEPARATE_MAIN_PRIMARY" if not near else "SEPARATE_NEAR_INTEGER_SENSITIVITY",
            )
        for sample in sorted(individuals):
            bundle_id = f"{sample}_individual_{suffix}"
            roles[bundle_id] = BundleRole(
                bundle_id,
                (sample,),
                near,
                "INDIVIDUAL_MAIN_PRIMARY" if not near else "INDIVIDUAL_NEAR_INTEGER_SENSITIVITY",
            )
    return roles


def verify_receipted_file(
    receipt: Mapping[str, object],
    label: str,
    hashes: dict[str, str],
) -> Path:
    require("path" in receipt and "sha256" in receipt, f"{label} receipt misses path/sha256")
    path = resolved(str(receipt["path"]))
    require(path.is_file(), f"{label} source is missing: {path}")
    actual_sha = sha256_file(path)
    require(actual_sha == receipt["sha256"],
            f"{label} SHA drift: expected {receipt['sha256']}, observed {actual_sha}")
    if "size_bytes" in receipt:
        require(path.stat().st_size == int(receipt["size_bytes"]), f"{label} size drift: {path}")
    previous = hashes.setdefault(str(path), actual_sha)
    require(previous == actual_sha, f"Conflicting SHA expectations for {path}")
    return path


def artifact_path_and_hash(
    artifact: Mapping[str, object],
    path_key: str,
    hash_key: str,
    label: str,
    hashes: dict[str, str],
) -> tuple[Path, str]:
    require(path_key in artifact and hash_key in artifact, f"{label} misses {path_key}/{hash_key}")
    path = resolved(str(artifact[path_key]))
    expected_sha = str(artifact[hash_key])
    require(path.is_file(), f"{label} file is missing: {path}")
    actual_sha = sha256_file(path)
    require(actual_sha == expected_sha,
            f"{label} SHA drift: expected {expected_sha}, observed {actual_sha}")
    hashes[str(path)] = actual_sha
    return path, actual_sha


def require_fields(reader: csv.DictReader, required: Iterable[str], label: str) -> None:
    missing = set(required).difference(reader.fieldnames or ())
    require(not missing, f"{label} misses columns: {sorted(missing)}")


def scan_bundle(
    role: BundleRole,
    metadata_path: Path,
    results_path: Path,
    target_keys: frozenset[CandidateKey],
    expected_mutations: int,
    expected_rows: int,
) -> Mapping[CandidateKey, FitHit]:
    metadata_by_mid: dict[str, VariantKey] = {}
    metadata_hits: dict[CandidateKey, Mapping[str, str]] = {}
    with open_text(metadata_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(reader, ("mutation_id", "chrom", "pos", "ref", "alt"), str(metadata_path))
        for row_number, row in enumerate(reader, start=2):
            variant = validate_variant(
                row["chrom"], row["pos"], row["ref"], row["alt"],
                f"{metadata_path}:{row_number}",
            )
            mid = row["mutation_id"]
            require(mid == variant.mutation_id,
                    f"Metadata mutation_id/allele mismatch at {metadata_path}:{row_number}")
            require(mid not in metadata_by_mid,
                    f"Duplicate metadata mutation_id in {role.bundle_id}: {mid}")
            metadata_by_mid[mid] = variant
            for sample in role.samples:
                candidate_key = CandidateKey(sample, variant)
                if candidate_key in target_keys:
                    metadata_hits[candidate_key] = dict(row)
    require(len(metadata_by_mid) == expected_mutations,
            f"{role.bundle_id} metadata count drift: {len(metadata_by_mid)} != {expected_mutations}")

    seen_pairs: set[tuple[str, str]] = set()
    result_hits: dict[CandidateKey, Mapping[str, str]] = {}
    with open_text(results_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(
            reader,
            (
                "mutation_id",
                "sample_id",
                "cluster_id",
                "cellular_prevalence",
                "cellular_prevalence_std",
                "cluster_assignment_prob",
            ),
            str(results_path),
        )
        for row_number, row in enumerate(reader, start=2):
            mid = row["mutation_id"]
            sample = row["sample_id"]
            require(sample in role.samples,
                    f"Unexpected sample in {role.bundle_id} results at row {row_number}: {sample}")
            require(mid in metadata_by_mid,
                    f"Result key lacks exact metadata allele in {role.bundle_id}: {sample}/{mid}")
            pair = (sample, mid)
            require(pair not in seen_pairs, f"Duplicate result key in {role.bundle_id}: {pair}")
            seen_pairs.add(pair)
            require(bool(row["cluster_id"]), f"Empty fit-local cluster_id in {role.bundle_id}: {pair}")
            cellular_prevalence = parse_float(
                row["cellular_prevalence"], f"{role.bundle_id}:{row_number}.cellular_prevalence"
            )
            prevalence_std = parse_float(
                row["cellular_prevalence_std"],
                f"{role.bundle_id}:{row_number}.cellular_prevalence_std",
            )
            assignment_probability = parse_float(
                row["cluster_assignment_prob"],
                f"{role.bundle_id}:{row_number}.cluster_assignment_prob",
            )
            require(0.0 <= cellular_prevalence <= 1.0,
                    f"Cellular prevalence out of range in {role.bundle_id}: {pair}")
            require(prevalence_std >= 0.0,
                    f"Cellular prevalence std is negative in {role.bundle_id}: {pair}")
            require(0.0 <= assignment_probability <= 1.0,
                    f"Assignment probability out of range in {role.bundle_id}: {pair}")
            candidate_key = CandidateKey(sample, metadata_by_mid[mid])
            if candidate_key in metadata_hits:
                result_hits[candidate_key] = dict(row)

    require(len(seen_pairs) == expected_rows,
            f"{role.bundle_id} result count drift: {len(seen_pairs)} != {expected_rows}")
    require(expected_rows == expected_mutations * len(role.samples),
            f"{role.bundle_id} is not a complete mutation/sample matrix")
    require(len(seen_pairs) == len(metadata_by_mid) * len(role.samples),
            f"{role.bundle_id} result matrix is incomplete")
    require(set(result_hits) == set(metadata_hits),
            f"{role.bundle_id} candidate metadata/result exact-key mismatch")
    return {
        key: FitHit(metadata_hits[key], result_hits[key])
        for key in sorted(metadata_hits)
    }


def compare_path(actual: object, expected: Path, label: str) -> None:
    require(resolved(str(actual)) == expected, f"{label} path drift: {actual!r} != {expected}")


def validate_bundle(
    role: BundleRole,
    bundle: Mapping[str, object],
    summary_receipt: Mapping[str, object],
    fit_summary: Mapping[str, object],
    target_keys: frozenset[CandidateKey],
    hashes: dict[str, str],
) -> BundleEvidence:
    require(bundle.get("bundle_id") == role.bundle_id, f"Bundle ID drift for {role.bundle_id}")
    require(bundle.get("status") == "PASS", f"Bundle authority is not PASS: {role.bundle_id}")
    require(tuple(bundle.get("samples", ())) == role.samples, f"Bundle sample drift: {role.bundle_id}")
    require(bool(bundle.get("near_integer_only")) == role.near_integer_only,
            f"Bundle model mode drift: {role.bundle_id}")
    artifacts = bundle.get("artifacts")
    counters = bundle.get("counters")
    require(isinstance(artifacts, dict) and isinstance(counters, dict),
            f"Bundle {role.bundle_id} misses artifacts/counters")
    input_path, input_sha = artifact_path_and_hash(
        artifacts, "pyclone_input", "pyclone_input_sha256", f"{role.bundle_id} input", hashes
    )
    metadata_path, metadata_sha = artifact_path_and_hash(
        artifacts, "site_metadata", "site_metadata_sha256", f"{role.bundle_id} metadata", hashes
    )
    qa_path = resolved(str(artifacts.get("qa_json", "")))
    require(qa_path.is_file(), f"Bundle QA file is missing: {qa_path}")
    qa = load_json(qa_path, f"{role.bundle_id} QA")
    qa_sha = sha256_file(qa_path)
    hashes[str(qa_path)] = qa_sha
    require(qa.get("status") == "PASS" and qa.get("bundle_id") == role.bundle_id,
            f"Bundle QA is not matching PASS: {role.bundle_id}")
    require(qa.get("counters") == counters, f"Bundle QA counter drift: {role.bundle_id}")
    qa_artifacts = qa.get("artifacts")
    require(isinstance(qa_artifacts, dict), f"Bundle QA artifacts missing: {role.bundle_id}")
    compare_path(qa_artifacts.get("pyclone_input"), input_path, f"{role.bundle_id} QA input")
    compare_path(qa_artifacts.get("site_metadata"), metadata_path, f"{role.bundle_id} QA metadata")
    require(qa_artifacts.get("pyclone_input_sha256") == input_sha,
            f"Bundle QA input SHA drift: {role.bundle_id}")
    require(qa_artifacts.get("site_metadata_sha256") == metadata_sha,
            f"Bundle QA metadata SHA drift: {role.bundle_id}")

    status_path = resolved(str(summary_receipt.get("status_path", "")))
    results_path = resolved(str(summary_receipt.get("results_path", "")))
    require(status_path.is_file() and results_path.is_file(),
            f"PyClone status/results source missing: {role.bundle_id}")
    results_sha = sha256_file(results_path)
    require(results_sha == summary_receipt.get("results_sha256"),
            f"Analysis summary results SHA drift: {role.bundle_id}")
    hashes[str(results_path)] = results_sha
    status = load_json(status_path, f"{role.bundle_id} status")
    status_sha = sha256_file(status_path)
    hashes[str(status_path)] = status_sha
    require(
        status.get("status") == "PASS"
        and status.get("bundle_id") == role.bundle_id
        and status.get("fit_exit_code") == 0
        and status.get("write_results_exit_code") == 0,
        f"PyClone fit receipt is not terminal PASS: {role.bundle_id}",
    )
    status_input = status.get("input")
    status_output = status.get("outputs")
    status_shape = status.get("result_shape")
    require(isinstance(status_input, dict) and isinstance(status_output, dict)
            and isinstance(status_shape, dict), f"Malformed fit receipt: {role.bundle_id}")
    compare_path(status_input.get("path"), input_path, f"{role.bundle_id} fit input")
    require(status_input.get("sha256") == input_sha, f"Fit input SHA drift: {role.bundle_id}")
    compare_path(status_input.get("qa_path"), qa_path, f"{role.bundle_id} fit QA")
    require(status_input.get("qa_sha256") == qa_sha, f"Fit QA SHA drift: {role.bundle_id}")
    compare_path(status_output.get("results"), results_path, f"{role.bundle_id} fit results")
    require(status_output.get("results_sha256") == results_sha,
            f"Fit results SHA drift: {role.bundle_id}")

    expected_mutations = parse_int(counters.get("included_mutations"),
                                   f"{role.bundle_id}.included_mutations", minimum=0)
    expected_rows = parse_int(counters.get("included_rows"),
                              f"{role.bundle_id}.included_rows", minimum=0)
    require(parse_int(status_shape.get("mutations"), f"{role.bundle_id}.status.mutations")
            == expected_mutations, f"Fit mutation count drift: {role.bundle_id}")
    require(parse_int(status_shape.get("rows"), f"{role.bundle_id}.status.rows")
            == expected_rows, f"Fit row count drift: {role.bundle_id}")
    require(tuple(status_shape.get("samples", ())) == role.samples,
            f"Fit sample set drift: {role.bundle_id}")
    require(fit_summary.get("bundle_id") == role.bundle_id,
            f"Analysis fit summary bundle drift: {role.bundle_id}")
    require(parse_int(fit_summary.get("mutations"), f"{role.bundle_id}.summary.mutations")
            == expected_mutations, f"Analysis mutation count drift: {role.bundle_id}")
    require(parse_int(fit_summary.get("rows"), f"{role.bundle_id}.summary.rows")
            == expected_rows, f"Analysis row count drift: {role.bundle_id}")
    require(tuple(fit_summary.get("samples", ())) == role.samples,
            f"Analysis sample set drift: {role.bundle_id}")

    hits = scan_bundle(
        role,
        metadata_path,
        results_path,
        target_keys,
        expected_mutations,
        expected_rows,
    )
    return BundleEvidence(
        role,
        hits,
        metadata_path,
        metadata_sha,
        results_path,
        results_sha,
        status_path,
        status_sha,
    )


def load_authority(paths: AuthorityPaths, target_keys: frozenset[CandidateKey]) -> AuthorityContext:
    config_path = resolved(paths.config)
    provenance_path = resolved(paths.input_provenance)
    summary_path = resolved(paths.analysis_summary)
    config = load_json(config_path, "authority config")
    provenance = load_json(provenance_path, "input provenance")
    summary = load_json(summary_path, "analysis summary")
    config_sha = sha256_file(config_path)
    provenance_sha = sha256_file(provenance_path)
    summary_sha = sha256_file(summary_path)
    hashes: dict[str, str] = {
        str(config_path): config_sha,
        str(provenance_path): provenance_sha,
        str(summary_path): summary_sha,
    }

    require(config.get("schema_name") == "intersubmod.pyclone_vi_validation_config",
            "Unexpected authority config schema")
    require(provenance.get("schema_name") == "intersubmod.pyclone_input_build_provenance",
            "Unexpected input provenance schema")
    require(summary.get("schema_name") == "intersubmod.pyclone_vi_conditional_clustering_analysis",
            "Unexpected analysis summary schema")
    require(provenance.get("config") == config, "Embedded authority config drift in provenance")

    source_receipts = provenance.get("source_receipts")
    require(isinstance(source_receipts, dict), "Input provenance misses source_receipts")
    for label, receipt in source_receipts.items():
        require(isinstance(receipt, dict), f"Malformed source receipt: {label}")
        verify_receipted_file(receipt, f"source_receipts[{label}]", hashes)
    config_receipt = source_receipts.get("config")
    require(isinstance(config_receipt, dict), "Input provenance lacks config source receipt")
    compare_path(config_receipt.get("path"), config_path, "authority config receipt")
    require(config_receipt.get("sha256") == config_sha, "Authority config SHA is not locked by provenance")

    roles = expected_bundle_roles(config)
    bundle_rows = provenance.get("bundles")
    require(isinstance(bundle_rows, list), "Input provenance bundles must be a list")
    bundle_map: dict[str, Mapping[str, object]] = {}
    for row in bundle_rows:
        require(isinstance(row, dict) and isinstance(row.get("bundle_id"), str),
                "Malformed input provenance bundle")
        bundle_id = str(row["bundle_id"])
        require(bundle_id not in bundle_map, f"Duplicate provenance bundle: {bundle_id}")
        bundle_map[bundle_id] = row
    require(set(bundle_map) == set(roles),
            f"Authority bundle set drift: missing={sorted(set(roles)-set(bundle_map))}, "
            f"extra={sorted(set(bundle_map)-set(roles))}")

    summary_receipts = summary.get("source_receipts")
    fit_summaries = summary.get("fit_summaries")
    require(isinstance(summary_receipts, dict) and isinstance(fit_summaries, dict),
            "Analysis summary misses source_receipts/fit_summaries")
    require(set(summary_receipts) == set(roles), "Analysis source receipt bundle set drift")
    require(set(fit_summaries) == set(roles), "Analysis fit summary bundle set drift")
    require(summary.get("all_ready") is True, "Analysis authority all_ready is not true")
    require(summary.get("pending_or_failed_bundles") == [],
            "Analysis authority has pending/failed bundles")
    require(parse_int(summary.get("expected_bundle_count"), "expected_bundle_count") == len(roles),
            "Analysis expected bundle count drift")
    require(parse_int(summary.get("pass_bundle_count"), "pass_bundle_count") == len(roles),
            "Analysis PASS bundle count drift")

    bundles: dict[str, BundleEvidence] = {}
    for bundle_id in sorted(roles):
        receipt = summary_receipts[bundle_id]
        fit_summary = fit_summaries[bundle_id]
        require(isinstance(receipt, dict) and isinstance(fit_summary, dict),
                f"Malformed analysis authority for {bundle_id}")
        bundles[bundle_id] = validate_bundle(
            roles[bundle_id],
            bundle_map[bundle_id],
            receipt,
            fit_summary,
            target_keys,
            hashes,
        )

    pair = config["hcc1395_pair"]
    individuals = config["individual_samples"]
    assert isinstance(pair, dict) and isinstance(individuals, dict)
    pair_cn_path = resolved(str(pair["copy_number_tsv"]))
    cn_paths = {
        "HCC1395": pair_cn_path,
        "HCC1395_DORADO": pair_cn_path,
        **{
            sample: resolved(str(sample_config["copy_number_tsv"]))
            for sample, sample_config in individuals.items()
            if isinstance(sample_config, dict)
        },
    }
    for sample, cn_path in cn_paths.items():
        require(str(cn_path) in hashes,
                f"CN source for {sample} is not hash-locked by input provenance: {cn_path}")
    cn_indexes = {
        sample: SegmentIndex(path)
        for sample, path in cn_paths.items()
        if sample != "HCC1395_DORADO"
    }
    cn_indexes["HCC1395_DORADO"] = cn_indexes["HCC1395"]
    purity_by_sample = {
        "HCC1395": parse_float(pair["purity"], "HCC1395 purity"),
        "HCC1395_DORADO": parse_float(pair["purity"], "HCC1395_DORADO purity"),
        **{
            sample: parse_float(sample_config["purity"], f"{sample} purity")
            for sample, sample_config in individuals.items()
            if isinstance(sample_config, dict)
        },
    }
    return AuthorityContext(
        config,
        hashes,
        bundles,
        cn_indexes,
        cn_paths,
        purity_by_sample,
        config_sha,
        provenance_sha,
        summary_sha,
    )


def primary_bundle_id(sample: str) -> str:
    if sample in ("HCC1395", "HCC1395_DORADO"):
        return f"hcc1395_pair_primary_separate_{sample}_main"
    return f"{sample}_individual_main"


def sensitivity_bundle_id(sample: str) -> str:
    if sample in ("HCC1395", "HCC1395_DORADO"):
        return f"hcc1395_pair_primary_separate_{sample}_near_integer"
    return f"{sample}_individual_near_integer"


def format_float(value: float) -> str:
    return f"{value:.10g}"


def fit_columns(hit: FitHit | None, prefix: str) -> Mapping[str, str]:
    if hit is None:
        return {
            f"{prefix}_fit_local_cluster_id": "",
            f"{prefix}_cellular_prevalence": "",
            f"{prefix}_cellular_prevalence_std": "",
            f"{prefix}_assignment_probability": "",
        }
    result = hit.result
    return {
        f"{prefix}_fit_local_cluster_id": result["cluster_id"],
        f"{prefix}_cellular_prevalence": result["cellular_prevalence"],
        f"{prefix}_cellular_prevalence_std": result["cellular_prevalence_std"],
        f"{prefix}_assignment_probability": result["cluster_assignment_prob"],
    }


def assert_cn_metadata_matches(segment: Segment, hit: FitHit, key: CandidateKey, bundle_id: str) -> None:
    metadata = hit.metadata
    required = {
        "segment_id",
        "segment_start",
        "segment_end",
        "total_cn_raw",
        "minor_cn_raw",
        "total_cn_discrete",
        "major_cn",
        "minor_cn",
        "near_integer",
    }
    missing = required.difference(metadata)
    require(not missing, f"{bundle_id} metadata lacks CN columns: {sorted(missing)}")
    require(
        metadata["segment_id"] == segment.segment_id
        and parse_int(metadata["segment_start"], "metadata.segment_start") == segment.start
        and parse_int(metadata["segment_end"], "metadata.segment_end") == segment.end
        and math.isclose(parse_float(metadata["total_cn_raw"], "metadata.total_cn_raw"),
                         segment.total_raw, rel_tol=0.0, abs_tol=1e-8)
        and math.isclose(parse_float(metadata["minor_cn_raw"], "metadata.minor_cn_raw"),
                         segment.minor_raw, rel_tol=0.0, abs_tol=1e-8)
        and parse_int(metadata["total_cn_discrete"], "metadata.total_cn_discrete")
        == segment.total_discrete
        and parse_int(metadata["major_cn"], "metadata.major_cn") == segment.major_discrete
        and parse_int(metadata["minor_cn"], "metadata.minor_cn") == segment.minor_discrete
        and bool(parse_int(metadata["near_integer"], "metadata.near_integer"))
        == segment.near_integer,
        f"Direct CN authority and exact metadata disagree for {key} in {bundle_id}",
    )


def build_annotation_row(key: CandidateKey, authority: AuthorityContext) -> Mapping[str, str]:
    sample = key.sample
    variant = key.variant
    row = {column: "" for column in OUTPUT_COLUMNS}
    row.update(
        {
            "sample": sample,
            "chrom": variant.chrom,
            "pos": str(variant.pos),
            "ref": variant.ref,
            "alt": variant.alt,
            "mutation_id": variant.mutation_id,
            "callset_status": "INPUT_CANDIDATE_KEY_VALIDATED_NOT_RECHECKED",
            "authority_config_sha256": authority.config_sha256,
            "input_provenance_sha256": authority.provenance_sha256,
            "analysis_summary_sha256": authority.summary_sha256,
            "claim_ceiling": CLAIM_CEILING,
        }
    )

    if sample in BLOCKED_SAMPLES:
        row.update(
            {
                "cn_status": "BLOCKED_CN_MISFIT",
                "cn_transfer_policy": "BLOCKED_CN_MISFIT_NO_CN2_IMPUTATION",
                "independent_cn": "false",
                "pyclone_status": "BLOCKED_CN_MISFIT",
                "pyclone_sensitivity_status": "BLOCKED_CN_MISFIT",
            }
        )
        return row

    cn_authority_sample = "HCC1395" if sample == "HCC1395_DORADO" else sample
    cn_path = authority.cn_paths[sample]
    segment = authority.cn_indexes[sample].find(variant)
    row.update(
        {
            "cn_authority_sample": cn_authority_sample,
            "cn_transfer_policy": (
                "SAME_CELL_LINE_SHARED_CN_SENSITIVITY"
                if sample == "HCC1395_DORADO"
                else "SAMPLE_SPECIFIC_MEASURED_CN"
            ),
            "independent_cn": "false" if sample == "HCC1395_DORADO" else "true",
            "purity_model_value": format_float(authority.purity_by_sample[sample]),
            "cn_source_sha256": authority.hashes[str(cn_path)],
        }
    )
    if segment is None:
        row["cn_status"] = "NO_CN_SEGMENT"
    else:
        row.update(
            {
                "cn_status": (
                    "SHARED_CN_SENSITIVITY"
                    if sample == "HCC1395_DORADO"
                    else "AVAILABLE_EXACT_SEGMENT"
                ),
                "cn_segment_id": segment.segment_id,
                "cn_segment_start": str(segment.start),
                "cn_segment_end": str(segment.end),
                "savana_total_cn_raw": format_float(segment.total_raw),
                "savana_minor_cn_raw": format_float(segment.minor_raw),
                "savana_total_cn_discrete": str(segment.total_discrete),
                "savana_major_cn_discrete": str(segment.major_discrete),
                "savana_minor_cn_discrete": str(segment.minor_discrete),
            }
        )

    main_id = primary_bundle_id(sample)
    sensitivity_id = sensitivity_bundle_id(sample)
    main_bundle = authority.bundles[main_id]
    sensitivity_bundle = authority.bundles[sensitivity_id]
    main_hit = main_bundle.hits.get(key)
    sensitivity_hit = sensitivity_bundle.hits.get(key)
    require(not (sensitivity_hit is not None and main_hit is None),
            f"Near-integer sensitivity hit lacks primary main hit: {key}")
    if main_hit is not None:
        require(segment is not None, f"PyClone primary hit lacks direct CN segment: {key}")
        assert segment is not None
        assert_cn_metadata_matches(segment, main_hit, key, main_id)
    if sensitivity_hit is not None:
        require(segment is not None and segment.near_integer,
                f"Near-integer PyClone hit lacks near-integer direct CN: {key}")
        assert segment is not None
        assert_cn_metadata_matches(segment, sensitivity_hit, key, sensitivity_id)

    row.update(
        {
            "pyclone_status": "MATCHED_PRIMARY" if main_hit is not None else "NOT_IN_FIT_UNIVERSE",
            "pyclone_primary_bundle_id": main_id,
            "pyclone_primary_model_role": main_bundle.role.model_role,
            "pyclone_sensitivity_status": (
                "MATCHED_SENSITIVITY_ONLY"
                if sensitivity_hit is not None
                else "NOT_IN_FIT_UNIVERSE"
            ),
            "pyclone_sensitivity_bundle_id": sensitivity_id,
            "pyclone_primary_metadata_sha256": main_bundle.metadata_sha256,
            "pyclone_primary_results_sha256": main_bundle.results_sha256,
            "pyclone_primary_status_sha256": main_bundle.status_sha256,
            "pyclone_sensitivity_metadata_sha256": sensitivity_bundle.metadata_sha256,
            "pyclone_sensitivity_results_sha256": sensitivity_bundle.results_sha256,
            "pyclone_sensitivity_status_sha256": sensitivity_bundle.status_sha256,
        }
    )
    primary_values = fit_columns(main_hit, "pyclone")
    row.update(
        {
            "pyclone_fit_local_cluster_id": primary_values["pyclone_fit_local_cluster_id"],
            "pyclone_vi_cellular_prevalence": primary_values["pyclone_cellular_prevalence"],
            "pyclone_vi_cellular_prevalence_std": primary_values[
                "pyclone_cellular_prevalence_std"
            ],
            "pyclone_vi_assignment_probability": primary_values[
                "pyclone_assignment_probability"
            ],
        }
    )
    sensitivity_values = fit_columns(sensitivity_hit, "pyclone_sensitivity")
    row.update(sensitivity_values)
    return row


def write_gzip_tsv(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as text_handle:
                writer = csv.DictWriter(
                    text_handle,
                    fieldnames=OUTPUT_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="raise",
                )
                writer.writeheader()
                writer.writerows(rows)


def status_counts(rows: Sequence[Mapping[str, str]], column: str) -> Mapping[str, int]:
    return dict(sorted(Counter(row[column] for row in rows).items()))


def annotate_candidates(
    input_path: Path,
    output_dir: Path,
    *,
    selection_column: str | None = None,
    selection_value: str | None = None,
    authority_paths: AuthorityPaths | None = None,
    command: Sequence[str] | None = None,
    source_authority: Mapping[str, object] | None = None,
) -> Mapping[str, object]:
    release_authority: Mapping[str, object] = source_authority or {
        "authority_id": "TEST_ONLY_UNSIGNED_AUTHORITY",
        "pass": True,
    }
    source_identity_before = capture_source_identity()
    source_modes_before = capture_source_modes()
    if (
        release_authority.get("authority_id") == SOURCE_AUTHORITY.AUTHORITY_ID
        and set(source_modes_before.values()) != {"0o444"}
    ):
        raise ContractError(f"Formal CN/CCF producer source is not mode 0444: {source_modes_before}")
    input_path = resolved(input_path)
    output_dir = output_dir.expanduser().absolute()
    require(not output_dir.exists(), f"Output directory must not exist: {output_dir}")
    candidates = read_candidates(input_path, selection_column, selection_value)
    target_keys = frozenset(candidates.selected)
    require(len(target_keys) == len(candidates.selected), "Selected candidate keys are not unique")
    authority = load_authority(authority_paths or AuthorityPaths.canonical(), target_keys)
    rows = [build_annotation_row(key, authority) for key in candidates.selected]
    require(len(rows) == len(candidates.selected), "Input/output row conservation failed")
    n_selected_candidates = len(candidates.selected)
    zero_selected = n_selected_candidates == 0
    receipt_status = "NOT_APPLICABLE" if zero_selected else "PASS"
    execution_status = "NOT_APPLICABLE" if zero_selected else "EXECUTION_PASS"
    reason = "ZERO_SELECTED_CANDIDATES" if zero_selected else None

    if release_authority.get("authority_id") == SOURCE_AUTHORITY.AUTHORITY_ID:
        release_authority_after = SOURCE_AUTHORITY.validate_release_source_authority(
            {"release_source_authority_validator", "cn_ccf_annotator"}
        )
        require(
            release_authority_after == release_authority,
            "CN/CCF source authority changed during execution",
        )
    source_identity_after = capture_source_identity()
    source_modes_after = capture_source_modes()
    require(
        source_identity_after == source_identity_before
        and source_modes_after == source_modes_before
        and set(source_modes_after.values()) == {"0o444"},
        "CN/CCF producer source identity or mode changed during execution",
    )
    source_lock = {
        "source_identity_before": source_identity_before,
        "source_identity_after_compute": source_identity_after,
        "source_modes_before": source_modes_before,
        "source_modes_after_compute": source_modes_after,
        "all_sources_read_only_and_unchanged": True,
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output_dir.name}.staging.", dir=output_dir.parent) as temp:
        stage = Path(temp) / "payload"
        stage.mkdir()
        output_path = stage / "candidate_cn_ccf_annotations.tsv.gz"
        receipt_path = stage / "receipt.json"
        write_gzip_tsv(output_path, rows)
        output_sha = sha256_file(output_path)
        code_path = resolved(Path(__file__))
        receipt: dict[str, object] = {
            "schema_name": RECEIPT_SCHEMA_NAME,
            "schema_version": RECEIPT_SCHEMA_VERSION,
            "pass": True,
            "status": receipt_status,
            "execution_status": execution_status,
            "reason": reason,
            "n_selected_candidates": n_selected_candidates,
            "command": list(command if command is not None else sys.argv),
            "source_authority": dict(release_authority),
            "code": source_identity_after,
            "source_lock": source_lock,
            "input": {
                "path": str(input_path),
                "sha256": sha256_file(input_path),
                "rows_read_total": candidates.rows_read_total,
                "rows_in": len(candidates.selected),
                "selection_column": selection_column,
                "selection_value": selection_value,
            },
            "output": {
                "path": str((output_dir / output_path.name).resolve()),
                "sha256": output_sha,
                "rows_out": len(rows),
                "columns": list(OUTPUT_COLUMNS),
            },
            "conservation": {
                "rows_in": len(candidates.selected),
                "rows_out": len(rows),
                "rows_in_equals_rows_out": len(candidates.selected) == len(rows),
            },
            "authority": {
                "config": {
                    "path": str(resolved((authority_paths or AuthorityPaths.canonical()).config)),
                    "sha256": authority.config_sha256,
                },
                "input_provenance": {
                    "path": str(resolved((authority_paths or AuthorityPaths.canonical()).input_provenance)),
                    "sha256": authority.provenance_sha256,
                },
                "analysis_summary": {
                    "path": str(resolved((authority_paths or AuthorityPaths.canonical()).analysis_summary)),
                    "sha256": authority.summary_sha256,
                },
                "all_source_hashes": dict(sorted(authority.hashes.items())),
            },
            "status_counts": {
                "cn_status": status_counts(rows, "cn_status"),
                "pyclone_status": status_counts(rows, "pyclone_status"),
                "pyclone_sensitivity_status": status_counts(rows, "pyclone_sensitivity_status"),
            },
            "status_enums": {
                "cn_status": list(CN_STATUS_ENUM),
                "pyclone_status": list(PYCLONE_STATUS_ENUM),
            },
            "claim_ceiling": CLAIM_CEILING,
            "pass_semantics": PASS_SEMANTICS,
            "scientific_interpretation": {
                "is_negative_result": False,
                "c1_formed": False,
                "statement": (
                    "ZERO_SELECTED_CANDIDATES is not a biological negative; C1 was not evaluated or formed."
                    if zero_selected
                    else "Conditional CN/PyClone annotation does not itself form C1."
                ),
            },
        }
        require(receipt["conservation"]["rows_in_equals_rows_out"] is True,
                "Receipt row conservation failed")
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        require(not output_dir.exists(), f"Output directory appeared during staging: {output_dir}")
        os.rename(stage, output_dir)
    (output_dir / "candidate_cn_ccf_annotations.tsv.gz").chmod(0o444)
    (output_dir / "receipt.json").chmod(0o444)
    return receipt


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Candidate TSV or TSV.GZ")
    parser.add_argument("--output-dir", type=Path, required=True, help="New output directory")
    parser.add_argument("--selection-column")
    parser.add_argument("--selection-value")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    source_authority = SOURCE_AUTHORITY.validate_release_source_authority(
        {"release_source_authority_validator", "cn_ccf_annotator"}
    )
    command = resolve_release_command(argv, source_authority)
    args = parse_args(argv)
    try:
        receipt = annotate_candidates(
            args.input,
            args.output_dir,
            selection_column=args.selection_column,
            selection_value=args.selection_value,
            command=command,
            source_authority=source_authority,
        )
    except ContractError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Input: {resolved(args.input)}")
    print(f"Output: {(args.output_dir / 'candidate_cn_ccf_annotations.tsv.gz').absolute()}")
    print(f"Receipt: {(args.output_dir / 'receipt.json').absolute()}")
    print(
        f"{receipt['status']} rows={receipt['output']['rows_out']} "
        f"sha256={receipt['output']['sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
