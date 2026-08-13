#!/usr/bin/env python3
"""Load, validate, and inspect InterSubMod site profiles."""

from __future__ import annotations

import hashlib
import argparse
import json
import os
import platform
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


PROFILE_SCHEMA_VERSION = "1.1.0"
RECEIPT_SCHEMA_VERSION = "1.1.0"
AUTOSOMES = tuple(str(index) for index in range(1, 23))

PROFILE_KEYS = {
    "REFERENCE": ("reference", "fasta", "path"),
    "PROJECT_ROOT_DEFAULT": ("project_root",),
    "SAMTOOLS": ("tools", "samtools", "path"),
    "BCFTOOLS": ("tools", "bcftools", "path"),
    "LONGPHASE_S_BIN": ("tools", "longphase_s", "path"),
    "INTERSUBMOD_BIN": ("tools", "intersubmod", "path"),
    "REFERENCE_GENOME_BUILD": ("reference", "genome_build"),
    "CONTIG_NAMING": ("reference", "contig_contract", "naming"),
    "CONTIG_SCOPE": ("reference", "contig_contract", "scope"),
}

SAMPLE_INPUT_KEYS = {
    "TUMOR_BAM": "tumor_bam",
    "TUMOR_BAM_INDEX": "tumor_bam_index",
    "NORMAL_BAM": "normal_bam",
    "NORMAL_BAM_INDEX": "normal_bam_index",
    "SOMATIC_VCF": "somatic_vcf",
    "SOMATIC_VCF_INDEX": "somatic_vcf_index",
    "SOMATIC_VCF_PILEUP": "somatic_vcf_pileup",
    "SOMATIC_VCF_PILEUP_INDEX": "somatic_vcf_pileup_index",
    "SOMATIC_VCF_INDEL": "somatic_vcf_indel",
    "SOMATIC_VCF_INDEL_INDEX": "somatic_vcf_indel_index",
    "TO_SOMATIC_VCF_PILEUP": "to_somatic_vcf_pileup",
    "TO_SOMATIC_VCF_PILEUP_INDEX": "to_somatic_vcf_pileup_index",
    "TO_SOMATIC_VCF": "to_somatic_vcf",
    "TO_SOMATIC_VCF_INDEX": "to_somatic_vcf_index",
    "TO_SOMATIC_VCF_INDEL": "to_somatic_vcf_indel",
    "TO_SOMATIC_VCF_INDEL_INDEX": "to_somatic_vcf_indel_index",
    "GERMLINE_PHASED_DIR": "germline_phased_dir",
    "TRUTH_VCF": "truth_vcf",
    "TRUTH_VCF_INDEX": "truth_vcf_index",
    "TRUTH_BED": "truth_bed",
}

SUPPORTED_TOOLS = {"samtools", "bcftools", "longphase_s", "intersubmod"}
SUPPORTED_INPUT_ROLES = {
    "tumor_bam", "tumor_bam_index", "normal_bam", "normal_bam_index",
    "somatic_vcf", "somatic_vcf_index", "somatic_vcf_pileup", "somatic_vcf_pileup_index",
    "somatic_vcf_indel", "somatic_vcf_indel_index", "to_somatic_vcf", "to_somatic_vcf_index",
    "to_somatic_vcf_pileup", "to_somatic_vcf_pileup_index", "to_somatic_vcf_indel",
    "to_somatic_vcf_indel_index", "germline_phased_dir", "truth_vcf", "truth_vcf_index",
    "truth_bed",
}
REQUIRED_WORKFLOW_ROLES = {
    "tumor_bam", "tumor_bam_index", "normal_bam", "normal_bam_index",
    "somatic_vcf", "somatic_vcf_pileup",
    "somatic_vcf_indel", "to_somatic_vcf", "to_somatic_vcf_pileup",
    "to_somatic_vcf_indel", "germline_phased_dir", "truth_vcf",
    "truth_bed",
}
RANDOM_ACCESS_VCF_INDEX_ROLES = {
    "somatic_vcf": "somatic_vcf_index",
    "somatic_vcf_pileup": "somatic_vcf_pileup_index",
    "somatic_vcf_indel": "somatic_vcf_indel_index",
    "to_somatic_vcf": "to_somatic_vcf_index",
    "to_somatic_vcf_pileup": "to_somatic_vcf_pileup_index",
    "to_somatic_vcf_indel": "to_somatic_vcf_indel_index",
    "truth_vcf": "truth_vcf_index",
}
CHECKSUM_POLICIES = {"EXPECTED_SHA256", "LOCATOR_ONLY", "SIZE_ONLY", "EXISTENCE_ONLY"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_profile(path: Path) -> dict[str, Any]:
    profile = json.loads(path.read_text(encoding="utf-8"))
    errors = validate_profile(profile)
    if errors:
        raise ValueError("; ".join(errors))
    return profile


def validate_profile(profile: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    required = {
        "schema_version",
        "site_id",
        "expected_hostname",
        "project_root",
        "data_roots",
        "reference",
        "tools",
        "datasets",
    }
    missing = sorted(required - set(profile))
    if missing:
        errors.append(f"missing top-level keys: {','.join(missing)}")
    if profile.get("schema_version") != PROFILE_SCHEMA_VERSION:
        errors.append(f"schema_version must be {PROFILE_SCHEMA_VERSION}")
    extra = sorted(set(profile) - required)
    if extra:
        errors.append(f"unknown top-level keys: {','.join(extra)}")
    for key in ("site_id", "expected_hostname"):
        if not isinstance(profile.get(key), str) or not profile.get(key):
            errors.append(f"{key} must be a non-empty string")

    def validate_path(label: str, value: Any) -> None:
        if not isinstance(value, str) or not value.startswith("/") or value.startswith("/<"):
            errors.append(f"{label} must be an explicit absolute path")

    def validate_file_spec(label: str, spec: Any, *, allow_required: bool = False) -> None:
        if not isinstance(spec, dict):
            errors.append(f"{label} must be an object")
            return
        allowed = {"path", "checksum_policy", "sha256", "checksum_uri", "size_bytes"}
        if allow_required:
            allowed.add("required")
        extra_spec = sorted(set(spec) - allowed)
        if extra_spec:
            errors.append(f"{label} has unknown keys: {','.join(extra_spec)}")
        validate_path(f"{label}.path", spec.get("path"))
        policy = spec.get("checksum_policy")
        if policy not in CHECKSUM_POLICIES:
            errors.append(f"{label}.checksum_policy is invalid")
        if policy == "EXPECTED_SHA256" and not re_full_sha256(spec.get("sha256")):
            errors.append(f"{label}.sha256 must contain 64 lowercase hex characters")
        if policy == "LOCATOR_ONLY" and not isinstance(spec.get("checksum_uri"), str):
            errors.append(f"{label}.checksum_uri is required for LOCATOR_ONLY")
        if policy != "EXPECTED_SHA256" and "sha256" in spec and not re_full_sha256(spec.get("sha256")):
            errors.append(f"{label}.sha256 must contain 64 lowercase hex characters")
        if "checksum_uri" in spec and (
            not isinstance(spec.get("checksum_uri"), str) or not spec["checksum_uri"]
        ):
            errors.append(f"{label}.checksum_uri must be a non-empty string")
        invalid_size = "size_bytes" in spec and (
            not isinstance(spec.get("size_bytes"), int)
            or isinstance(spec.get("size_bytes"), bool)
            or spec["size_bytes"] <= 0
        )
        if invalid_size or (policy == "SIZE_ONLY" and "size_bytes" not in spec):
            errors.append(f"{label}.size_bytes must be a positive integer for SIZE_ONLY")
        if allow_required and "required" not in spec:
            errors.append(f"{label}.required must be explicitly true")
        if allow_required and "required" in spec and not isinstance(spec["required"], bool):
            errors.append(f"{label}.required must be boolean")

    def validate_directory_spec(label: str, spec: Any) -> None:
        if not isinstance(spec, dict) or set(spec) != {"path"}:
            errors.append(f"{label} must contain exactly path")
            return
        validate_path(f"{label}.path", spec.get("path"))

    validate_path("project_root", profile.get("project_root"))
    data_roots = profile.get("data_roots")
    if not isinstance(data_roots, dict) or not {"primary", "output"}.issubset(data_roots):
        errors.append("data_roots must include primary and output")
    elif data_roots:
        for name, path in data_roots.items():
            validate_path(f"data_roots.{name}", path)

    reference = profile.get("reference")
    reference_keys = {"genome_build", "contig_contract", "fasta", "indexes"}
    if not isinstance(reference, dict) or set(reference) != reference_keys:
        errors.append(f"reference must contain exactly {sorted(reference_keys)}")
    else:
        if not isinstance(reference.get("genome_build"), str) or not reference["genome_build"]:
            errors.append("reference.genome_build must be a non-empty string")
        contract = reference.get("contig_contract")
        if not isinstance(contract, dict) or set(contract) != {"naming", "scope", "required_contigs"}:
            errors.append("reference.contig_contract must contain exactly naming, scope, and required_contigs")
        else:
            naming = contract.get("naming")
            if naming not in {"CHR_PREFIXED", "NO_CHR_PREFIX"}:
                errors.append("reference.contig_contract.naming is invalid")
            if contract.get("scope") != "AUTOSOMES_1_22":
                errors.append("reference.contig_contract.scope must be AUTOSOMES_1_22")
            expected_contigs = [f"chr{item}" for item in AUTOSOMES] if naming == "CHR_PREFIXED" else list(AUTOSOMES)
            if contract.get("required_contigs") != expected_contigs:
                errors.append("reference.contig_contract.required_contigs must exactly list autosomes 1-22 in order")
        validate_file_spec("reference.fasta", reference.get("fasta"))
        indexes = reference.get("indexes")
        if not isinstance(indexes, list) or not indexes:
            errors.append("reference.indexes must be a non-empty array")
        else:
            for index, spec in enumerate(indexes, start=1):
                validate_file_spec(f"reference.indexes[{index}]", spec)

    tools = profile.get("tools")
    if not isinstance(tools, dict) or len(tools) < 3:
        errors.append("tools must contain at least three tool specifications")
    elif tools:
        unknown_tools = sorted(set(tools) - SUPPORTED_TOOLS)
        if unknown_tools:
            errors.append(f"tools has unknown roles: {','.join(unknown_tools)}")
        for name, spec in tools.items():
            validate_file_spec(f"tools.{name}", spec, allow_required=True)
            if name in SUPPORTED_TOOLS and isinstance(spec, dict) and spec.get("required", True) is not True:
                errors.append(f"tools.{name}.required cannot disable a workflow-required tool")
            if isinstance(spec, dict) and spec.get("checksum_policy") == "EXISTENCE_ONLY":
                errors.append(f"tools.{name}.checksum_policy cannot be EXISTENCE_ONLY")
        missing_tools = sorted(SUPPORTED_TOOLS - set(tools))
        if missing_tools:
            errors.append(f"tools missing required roles: {','.join(missing_tools)}")
    if not isinstance(profile.get("datasets"), dict) or not profile.get("datasets"):
        errors.append("datasets must be a non-empty object")
    else:
        for dataset_id, dataset in profile["datasets"].items():
            if not isinstance(dataset, dict):
                errors.append(f"datasets.{dataset_id} must be an object")
                continue
            required_dataset = {
                "biological_id", "technical_dataset_id", "genome_build", "platform_label",
                "truth_set_label", "truth_total", "inputs",
            }
            if set(dataset) != required_dataset:
                errors.append(f"datasets.{dataset_id} must contain exactly {sorted(required_dataset)}")
            for key in ("biological_id", "technical_dataset_id", "genome_build", "platform_label", "truth_set_label"):
                if not isinstance(dataset.get(key), str) or not dataset.get(key):
                    errors.append(f"datasets.{dataset_id}.{key} must be a non-empty string")
            if isinstance(reference, dict) and dataset.get("genome_build") != reference.get("genome_build"):
                errors.append(f"datasets.{dataset_id}.genome_build must equal reference.genome_build")
            truth_total = dataset.get("truth_total")
            if not isinstance(truth_total, int) or isinstance(truth_total, bool) or truth_total <= 0:
                errors.append(f"datasets.{dataset_id}.truth_total must be a positive integer")
            inputs = dataset.get("inputs")
            if not isinstance(inputs, dict):
                errors.append(f"datasets.{dataset_id}.inputs must be an object")
            else:
                unknown_roles = sorted(set(inputs) - SUPPORTED_INPUT_ROLES)
                missing_roles = sorted(REQUIRED_WORKFLOW_ROLES - set(inputs))
                if unknown_roles:
                    errors.append(f"datasets.{dataset_id}.inputs has unknown roles: {','.join(unknown_roles)}")
                if missing_roles:
                    errors.append(f"datasets.{dataset_id}.inputs missing workflow roles: {','.join(missing_roles)}")
                for vcf_role, index_role in RANDOM_ACCESS_VCF_INDEX_ROLES.items():
                    vcf_spec = inputs.get(vcf_role)
                    path = vcf_spec.get("path") if isinstance(vcf_spec, dict) else None
                    if isinstance(path, str) and path.lower().endswith((".vcf.gz", ".vcf.bgz", ".bcf")):
                        if index_role not in inputs:
                            errors.append(
                                f"datasets.{dataset_id}.inputs.{index_role} is required for random-access {vcf_role}"
                            )
                for name, spec in inputs.items():
                    if name == "germline_phased_dir":
                        validate_directory_spec(f"datasets.{dataset_id}.inputs.{name}", spec)
                    else:
                        validate_file_spec(f"datasets.{dataset_id}.inputs.{name}", spec)
    return errors


def re_full_sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def nested(profile: dict[str, Any], keys: tuple[str, ...]) -> Any:
    value: Any = profile
    for key in keys:
        value = value[key]
    return value


def dataset_for_sample(profile: dict[str, Any], sample: str) -> dict[str, Any]:
    exact = profile["datasets"].get(sample)
    if exact:
        return exact
    candidates = [
        item
        for item in profile["datasets"].values()
        if item.get("technical_dataset_id") == sample or item.get("biological_id") == sample
    ]
    if len(candidates) != 1:
        raise ValueError(f"sample {sample!r} resolves to {len(candidates)} datasets")
    return candidates[0]


def dataset_id_for_sample(profile: dict[str, Any], sample: str) -> str:
    if sample in profile["datasets"]:
        return sample
    candidates = [
        key
        for key, item in profile["datasets"].items()
        if item.get("technical_dataset_id") == sample or item.get("biological_id") == sample
    ]
    if len(candidates) != 1:
        raise ValueError(f"sample {sample!r} resolves to {len(candidates)} datasets")
    return candidates[0]


def shell_assignments(profile: dict[str, Any], sample: str) -> str:
    # Emit every known role, including absent optional index roles. This makes
    # sourcing the output deterministic and prevents inherited machine values
    # from satisfying a missing profile field.
    values: dict[str, str] = {
        variable: "" for variable in (*PROFILE_KEYS.keys(), *SAMPLE_INPUT_KEYS.keys())
    }
    for variable, keys in PROFILE_KEYS.items():
        try:
            values[variable] = str(nested(profile, keys))
        except KeyError:
            pass
    output_root = str(profile["data_roots"]["output"])
    values["BIG7_OUTPUT_ROOT"] = output_root
    values["CANONICAL_OUTPUT_ROOT"] = str(Path(output_root) / "canonical")
    values["OUTPUT_ROOT"] = str(Path(output_root) / "canonical")
    dataset = dataset_for_sample(profile, sample)
    values["PLATFORM_LABEL"] = str(dataset["platform_label"])
    values["TRUTH_SET_LABEL"] = str(dataset["truth_set_label"])
    values["TRUTH_TOTAL"] = str(dataset["truth_total"])
    for variable, input_key in SAMPLE_INPUT_KEYS.items():
        if input_key in dataset["inputs"]:
            values[variable] = str(dataset["inputs"][input_key]["path"])
    return "\n".join(f"{key}={shlex.quote(value)}" for key, value in sorted(values.items()))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    shell_parser = subparsers.add_parser("shell", help="Emit shell-safe assignments")
    shell_parser.add_argument("--profile", required=True, type=Path)
    shell_parser.add_argument("--sample", required=True)
    preflight_parser = subparsers.add_parser(
        "preflight", help="Emit a machine-readable, side-effect-free preflight receipt"
    )
    preflight_parser.add_argument("--profile", required=True, type=Path)
    preflight_parser.add_argument("--sample", required=True)
    preflight_parser.add_argument(
        "--mode", choices=("synthetic", "real-preflight"), default="real-preflight"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        profile = load_profile(args.profile)
        if args.command == "shell":
            print(shell_assignments(profile, args.sample))
            return 0
        if args.command == "preflight":
            dataset_id = dataset_id_for_sample(profile, args.sample)
            result = inspect_profile(
                profile,
                include_real_data=args.mode == "real-preflight",
                profile_path=args.profile,
                dataset_ids={dataset_id},
                receipt_type="site_preflight",
            )
            print(json.dumps(result, ensure_ascii=False, sort_keys=True))
            return 0 if result["pass"] else 5
    except (OSError, json.JSONDecodeError, ValueError) as error:
        if args.command == "shell":
            print(f"[ERROR] invalid site profile: {error}", file=sys.stderr)
        else:
            print(
                json.dumps(
                    {
                        "schema_version": RECEIPT_SCHEMA_VERSION,
                        "receipt_type": "site_preflight_error",
                        "pass": False,
                        "error": str(error),
                    },
                    ensure_ascii=False,
                    sort_keys=True,
                )
            )
        return 3
    return 2


def mount_info(path: Path) -> dict[str, Any]:
    try:
        result = subprocess.run(
            ["findmnt", "--json", "--target", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        filesystems = json.loads(result.stdout).get("filesystems", [])
        if filesystems:
            item = filesystems[0]
            return {key: item.get(key) for key in ("target", "source", "fstype", "options")}
    except (FileNotFoundError, subprocess.CalledProcessError, json.JSONDecodeError):
        pass
    return {"target": None, "source": None, "fstype": None, "options": None}


def inspect_file(
    label: str,
    spec: dict[str, Any],
    executable: bool = False,
    always_hash: bool = False,
) -> dict[str, Any]:
    path = Path(spec["path"])
    exists = path.is_file()
    record: dict[str, Any] = {
        "label": label,
        "path": str(path),
        "exists": exists,
        "executable": bool(exists and os.access(path, os.X_OK)),
        "size_bytes": path.stat().st_size if exists else None,
        "checksum_policy": spec.get("checksum_policy"),
        "checksum_uri": spec.get("checksum_uri"),
        "expected_sha256": spec.get("sha256"),
        "expected_size_bytes": spec.get("size_bytes"),
        "actual_sha256": None,
        "status": "MISSING" if not exists else "PRESENT",
    }
    if exists and always_hash:
        record["actual_sha256"] = sha256_file(path)
    if exists and path.stat().st_size == 0:
        record["status"] = "EMPTY"
    elif exists and executable and not record["executable"]:
        record["status"] = "NOT_EXECUTABLE"
    elif exists and spec.get("checksum_policy") == "EXPECTED_SHA256":
        if record["actual_sha256"] is None:
            record["actual_sha256"] = sha256_file(path)
        record["status"] = (
            "MATCH" if record["actual_sha256"] == spec.get("sha256") else "HASH_MISMATCH"
        )
    elif exists and spec.get("checksum_policy") == "SIZE_ONLY":
        record["status"] = "SIZE_MATCH" if path.stat().st_size == spec.get("size_bytes") else "SIZE_MISMATCH"
    elif exists and spec.get("checksum_policy") == "LOCATOR_ONLY":
        locator = spec.get("checksum_uri")
        if not locator:
            record["status"] = "CHECKSUM_LOCATOR_MISSING"
        elif str(locator).startswith("/") and not Path(locator).is_file():
            record["status"] = "CHECKSUM_LOCATOR_MISSING"
        elif str(locator).startswith("/"):
            record["status"] = "CHECKSUM_LOCATOR_PRESENT"
        else:
            record["status"] = "EXTERNAL_CHECKSUM_LOCATOR_DECLARED"
    return record


def inspect_directory(label: str, path: Path) -> dict[str, Any]:
    exists = path.is_dir()
    return {
        "label": label,
        "path": str(path),
        "exists": exists,
        "status": "PRESENT" if exists else "MISSING",
        "mount": mount_info(path) if exists else {"target": None, "source": None, "fstype": None, "options": None},
    }


def inspect_reference_contigs(profile: dict[str, Any]) -> dict[str, Any]:
    required = profile["reference"]["contig_contract"]["required_contigs"]
    fai_candidates = [
        Path(spec["path"])
        for spec in profile["reference"]["indexes"]
        if str(spec["path"]).endswith(".fai")
    ]
    if len(fai_candidates) != 1 or not fai_candidates[0].is_file():
        return {
            "label": "reference:contig_contract",
            "path": str(fai_candidates[0]) if fai_candidates else None,
            "status": "CONTIG_CONTRACT_MISMATCH",
            "required_contigs": required,
            "missing_contigs": required,
            "reason": "exactly one readable FASTA .fai is required",
        }
    observed = {
        line.split("\t", 1)[0]
        for line in fai_candidates[0].read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    missing = [contig for contig in required if contig not in observed]
    return {
        "label": "reference:contig_contract",
        "path": str(fai_candidates[0]),
        "status": "CONTIG_CONTRACT_MATCH" if not missing else "CONTIG_CONTRACT_MISMATCH",
        "required_contigs": required,
        "observed_contig_count": len(observed),
        "missing_contigs": missing,
    }


def inspect_profile(
    profile: dict[str, Any],
    include_real_data: bool,
    profile_path: Path | None = None,
    dataset_ids: set[str] | None = None,
    receipt_type: str = "site_doctor",
) -> dict[str, Any]:
    hostname = platform.node()
    checks: list[dict[str, Any]] = []
    checks.append(inspect_directory("project_root", Path(profile["project_root"])))
    for name, path in sorted(profile["data_roots"].items()):
        checks.append(inspect_directory(f"data_root:{name}", Path(path)))
    for name, spec in sorted(profile["tools"].items()):
        if spec.get("required", True):
            checks.append(inspect_file(f"tool:{name}", spec, executable=True, always_hash=True))
    checks.append(inspect_file("reference:fasta", profile["reference"]["fasta"]))
    for index, spec in enumerate(profile["reference"]["indexes"], start=1):
        checks.append(inspect_file(f"reference:index:{index}", spec))
    checks.append(inspect_reference_contigs(profile))
    if include_real_data:
        for dataset_id, dataset in sorted(profile["datasets"].items()):
            if dataset_ids is not None and dataset_id not in dataset_ids:
                continue
            for name, spec in sorted(dataset["inputs"].items()):
                if name == "germline_phased_dir":
                    checks.append(inspect_directory(f"dataset:{dataset_id}:{name}", Path(spec["path"])))
                else:
                    checks.append(inspect_file(f"dataset:{dataset_id}:{name}", spec))

    failure_states = {
        "MISSING", "EMPTY", "NOT_EXECUTABLE", "HASH_MISMATCH", "SIZE_MISMATCH",
        "CHECKSUM_LOCATOR_MISSING", "CONTIG_CONTRACT_MISMATCH",
    }
    hostname_match = hostname == profile["expected_hostname"]
    resolved_profile = profile_path.resolve(strict=True) if profile_path is not None else None
    tool_hashes = {
        item["label"].removeprefix("tool:"): item["actual_sha256"]
        for item in checks
        if item["label"].startswith("tool:")
    }
    return {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "receipt_type": receipt_type,
        "verified_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "site_id": profile["site_id"],
        "profile": {
            "path": str(resolved_profile) if resolved_profile is not None else None,
            "sha256": sha256_file(resolved_profile) if resolved_profile is not None else None,
            "schema_version": profile["schema_version"],
        },
        "hostname": hostname,
        "expected_hostname": profile["expected_hostname"],
        "hostname_match": hostname_match,
        "mode": "real-data-read-only-preflight" if include_real_data else "synthetic-preflight",
        "mounts": {
            name: mount_info(Path(path)) for name, path in sorted(profile["data_roots"].items())
        },
        "tool_hashes": tool_hashes,
        "selected_datasets": sorted(dataset_ids) if dataset_ids is not None else sorted(profile["datasets"]),
        "checks": checks,
        "tally": {
            state: sum(item["status"] == state for item in checks)
            for state in sorted({item["status"] for item in checks})
        },
        "pass": hostname_match and not any(item["status"] in failure_states for item in checks),
        "claim_ceiling": (
            "Read-only path/index/tool-integrity preflight; no scientific computation and no "
            "independent-host claim unless hostname matches."
        ),
    }


if __name__ == "__main__":
    raise SystemExit(main())
