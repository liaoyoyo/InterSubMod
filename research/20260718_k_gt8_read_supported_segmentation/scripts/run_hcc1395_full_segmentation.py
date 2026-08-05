#!/usr/bin/env python3
"""Run HCC1395 chr1-22 extraction and k>8 segmentation sequentially.

The runner is deliberately fail-closed:

* the frozen manifest, consumed input authorities, and both worker scripts are
  authenticated before any chromosome starts;
* an existing run root is accepted only with ``--resume`` and an exact,
  byte-authenticated run contract;
* a chromosome stage is reused only when its stage receipt, child receipt,
  scope, commands, inputs, outputs, and SHA-256 sidecars all validate;
* incomplete or altered stage output is never deleted or overwritten;
* ``_SUCCESS`` is written only for the exact canonical HCC1395 chr1-22 scope.

All chromosomes, including chr21 where the canonical k>8 target count is zero,
run the extraction stage.  Segmentation is then recorded as
``SKIP_NO_K_GT8_TARGET`` only after the extracted site catalog independently
confirms zero k>8 components.  This keeps full extraction runtime measurable.

Scheduling priority is intentionally inherited.  Apply ``nice`` or ``ionice``
to this runner from the outer launch command when required.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import shutil
import stat
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import pysam


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MANIFEST = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/"
    "input_manifest.snapshot.json"
)
DEFAULT_MANIFEST_SHA256 = (
    "16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
)
DEFAULT_EXTRACTOR = SCRIPT_DIR / "extract_lossless_read_linkage_collapsing.py"
DEFAULT_PARTITIONER = SCRIPT_DIR / "build_k_gt8_partitions.py"
AUTOSOMES = tuple(f"chr{value}" for value in range(1, 23))
TIME_BINARY = Path("/usr/bin/time")
DATASET = "HCC1395"
SCHEMA_NAME = "intersubmod.hcc1395_full_k_gt8_segmentation"
SCHEMA_VERSION = "1.0.0"

# Independent census from the canonical frozen LongPhase-S recalibrated PASS
# tree VCF.  The runner re-derives these values from the authenticated VCF and
# then derives them again from each extracted site catalog.
CANONICAL_EXPECTED: dict[str, tuple[int, int, int, int]] = {
    "chr1": (3045, 13, 146, 18),
    "chr2": (3074, 17, 253, 49),
    "chr3": (2381, 9, 103, 23),
    "chr4": (2598, 21, 245, 23),
    "chr5": (2305, 13, 146, 17),
    "chr6": (27099, 83, 25657, 3574),
    "chr7": (2845, 46, 523, 23),
    "chr8": (3178, 59, 1205, 153),
    "chr9": (1523, 6, 80, 22),
    "chr10": (1504, 9, 106, 15),
    "chr11": (1328, 2, 28, 18),
    "chr12": (1573, 13, 160, 43),
    "chr13": (1019, 5, 52, 12),
    "chr14": (1343, 10, 156, 56),
    "chr15": (1154, 9, 123, 28),
    "chr16": (18819, 63, 18162, 2674),
    "chr17": (1093, 10, 188, 65),
    "chr18": (991, 5, 53, 13),
    "chr19": (761, 3, 27, 9),
    "chr20": (1075, 6, 59, 11),
    "chr21": (436, 0, 0, 0),
    "chr22": (543, 6, 98, 40),
}


class ContractError(RuntimeError):
    """Raised when an immutable input, receipt, or resume contract fails."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_slice(path: Path, offset: int, length: int) -> str:
    with path.open("rb") as handle:
        handle.seek(offset)
        payload = handle.read(length)
    if len(payload) != length:
        raise ContractError(
            f"short sampled chunk: {path} offset={offset} "
            f"expected={length} observed={len(payload)}"
        )
    return hashlib.sha256(payload).hexdigest()


def strict_json_load(path: Path) -> Any:
    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ContractError(f"duplicate JSON key in {path}: {key}")
            result[key] = value
        return result

    try:
        return json.loads(
            path.read_text(encoding="utf-8", errors="strict"),
            object_pairs_hook=pairs_hook,
            parse_constant=lambda value: (_ for _ in ()).throw(
                ContractError(f"non-finite JSON value in {path}: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot read strict JSON {path}: {exc}") from exc


def require_regular_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise ContractError(f"{label} is missing or not a regular file: {path}")


def file_identity(path: Path) -> dict[str, Any]:
    require_regular_file(path, "file")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def write_json_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    if temporary.exists():
        raise ContractError(f"refusing to reuse temporary path: {temporary}")
    payload = json.dumps(
        value, ensure_ascii=False, allow_nan=False, indent=2, sort_keys=True
    )
    with temporary.open("x", encoding="utf-8") as handle:
        handle.write(payload)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_sha256_sidecar(path: Path) -> Path:
    sidecar = path.with_name(f"{path.name}.sha256")
    temporary = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    with temporary.open("x", encoding="ascii") as handle:
        handle.write(f"{sha256_path(path)}  {path.name}\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, sidecar)
    return sidecar


def verify_sha256_sidecar(path: Path) -> str:
    require_regular_file(path, "receipt")
    sidecar = path.with_name(f"{path.name}.sha256")
    require_regular_file(sidecar, "receipt SHA-256 sidecar")
    fields = sidecar.read_text(encoding="ascii").strip().split()
    if len(fields) != 2 or fields[1] != path.name:
        raise ContractError(f"malformed SHA-256 sidecar: {sidecar}")
    observed = sha256_path(path)
    if fields[0] != observed:
        raise ContractError(
            f"SHA-256 sidecar mismatch: {path}; "
            f"expected={fields[0]} observed={observed}"
        )
    return observed


def verify_full_identity(spec: Mapping[str, Any], label: str) -> dict[str, Any]:
    if set(spec) != {"path", "identity"}:
        raise ContractError(f"{label} must contain exactly path and identity")
    path = Path(str(spec["path"]))
    require_regular_file(path, label)
    expected = spec["identity"]
    if not isinstance(expected, Mapping) or expected.get("policy") != "full_sha256":
        raise ContractError(f"{label} is not full_sha256")
    observed = {
        "policy": "full_sha256",
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }
    if observed != dict(expected):
        raise ContractError(
            f"{label} identity drift: expected={dict(expected)} observed={observed}"
        )
    return {"path": str(path), "identity": observed}


def verify_storage_identity(
    alignment: Mapping[str, Any],
) -> dict[str, Any]:
    required = {"path", "index_path", "storage_identity_v1"}
    if not required.issubset(alignment):
        raise ContractError(
            f"alignment_payload lacks required keys: expected={sorted(required)} "
            f"observed={sorted(alignment)}"
        )
    if alignment.get("identity_policy", "storage_identity_v1") != "storage_identity_v1":
        raise ContractError("alignment_payload identity_policy is not storage_identity_v1")
    if alignment.get("embedded_hp_ps_policy", "ignore") != "ignore":
        raise ContractError("alignment_payload embedded HP/PS policy must be ignore")
    path = Path(str(alignment["path"]))
    index_path = Path(str(alignment["index_path"]))
    expected = alignment["storage_identity_v1"]
    if not isinstance(expected, Mapping):
        raise ContractError("alignment storage_identity_v1 must be an object")
    require_regular_file(path, "alignment BAM")
    require_regular_file(index_path, "alignment BAM index")
    logical = path.lstat()
    resolved = path.resolve(strict=True)
    target = resolved.stat()
    chunks = []
    for raw in expected.get("chunks", []):
        if not isinstance(raw, Mapping):
            raise ContractError("alignment sampled chunk must be an object")
        offset = int(raw["offset"])
        length = int(raw["length"])
        chunks.append(
            {
                "label": raw["label"],
                "offset": offset,
                "length": length,
                "sha256": sha256_slice(resolved, offset, length),
            }
        )
    index_expected = expected.get("index")
    if not isinstance(index_expected, Mapping):
        raise ContractError("alignment storage identity lacks index")
    index_observed = verify_full_identity(index_expected, "alignment BAM index")
    observed: dict[str, Any] = {
        "policy": "storage_identity_v1",
        "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "is_full_content_hash": False,
        "requested_path": str(path),
        "realpath": str(resolved),
        "logical_is_symlink": stat.S_ISLNK(logical.st_mode),
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "st_dev": target.st_dev,
        "st_ino": target.st_ino,
        "size_bytes": target.st_size,
        "mtime_ns": target.st_mtime_ns,
        "ctime_ns": target.st_ctime_ns,
        "chunk_size_bytes": expected.get("chunk_size_bytes"),
        "chunks": chunks,
        "index": index_observed,
    }
    observed["identity_sha256"] = canonical_sha256(observed)
    if observed != dict(expected):
        raise ContractError("alignment_payload storage_identity_v1 drift")
    if str(index_path) != str(index_expected.get("path")):
        raise ContractError("alignment_payload index_path disagrees with storage identity")
    return {
        "path": str(path),
        "index_path": str(index_path),
        "identity_policy": alignment.get("identity_policy", "storage_identity_v1"),
        "embedded_hp_ps_policy": alignment.get("embedded_hp_ps_policy", "ignore"),
        "storage_identity_v1": observed,
    }


def sample_entry(manifest: Mapping[str, Any]) -> Mapping[str, Any]:
    samples = manifest.get("samples")
    if not isinstance(samples, list):
        raise ContractError("manifest samples must be an array")
    matches = [
        item
        for item in samples
        if isinstance(item, Mapping) and item.get("sample") == DATASET
    ]
    if len(matches) != 1:
        raise ContractError(f"manifest must contain exactly one {DATASET} entry")
    return matches[0]


def verify_consumed_authorities(entry: Mapping[str, Any]) -> dict[str, Any]:
    alignment = entry.get("alignment_payload")
    somatic = entry.get("somatic")
    tags = entry.get("read_tags")
    if not all(isinstance(value, Mapping) for value in (alignment, somatic, tags)):
        raise ContractError("manifest HCC1395 entry lacks input authority objects")
    if tags.get("duplicate_identity_policy") != (
        "collapse_redundant_rows_with_identical_HP_PS"
    ):
        raise ContractError("unsupported duplicate_identity_policy")
    tree_vcf = somatic.get("tree_vcf")
    sidecar = tags.get("sidecar")
    sidecar_index = tags.get("index")
    if not all(
        isinstance(value, Mapping)
        for value in (tree_vcf, sidecar, sidecar_index)
    ):
        raise ContractError("manifest HCC1395 entry lacks tree VCF or sidecar")
    tree_index = tree_vcf.get("index")
    if not isinstance(tree_index, Mapping):
        raise ContractError("tree VCF lacks index identity")
    tree_main = {"path": tree_vcf["path"], "identity": tree_vcf["identity"]}
    return {
        "alignment_payload": verify_storage_identity(alignment),
        "tree_vcf": {
            **verify_full_identity(tree_main, "tree VCF"),
            "index": verify_full_identity(tree_index, "tree VCF index"),
        },
        "read_tag_sidecar": verify_full_identity(sidecar, "read-tag sidecar"),
        "read_tag_sidecar_index": verify_full_identity(
            sidecar_index, "read-tag sidecar index"
        ),
        "duplicate_identity_policy": tags["duplicate_identity_policy"],
    }


def positional_inventory(
    positions: Sequence[int], legacy_gap_bp: int, max_block_size: int
) -> dict[str, Any]:
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise ContractError("sSNV positions are not strictly increasing")
    components: list[list[int]] = []
    for position in positions:
        if not components or position - components[-1][-1] > legacy_gap_bp:
            components.append([position])
        else:
            components[-1].append(position)
    targets = [component for component in components if len(component) > max_block_size]
    return {
        "ssnv_sites": len(positions),
        "positional_components_all": len(components),
        "k_gt8_components": len(targets),
        "k_gt8_sites": sum(len(component) for component in targets),
        "k_gt8_max_k": max((len(component) for component in targets), default=0),
        "positions_sha256": canonical_sha256(list(positions)),
    }


def derive_vcf_inventory(
    vcf_path: Path,
    chromosomes: Sequence[str],
    legacy_gap_bp: int,
    max_block_size: int,
) -> list[dict[str, Any]]:
    rows = []
    with pysam.VariantFile(str(vcf_path)) as source:
        for chrom in chromosomes:
            positions = []
            try:
                records: Iterable[Any] = source.fetch(chrom)
            except (ValueError, OSError) as exc:
                raise ContractError(f"cannot fetch {chrom} from tree VCF: {exc}") from exc
            for record in records:
                alts = tuple(record.alts or ())
                if len(record.ref) != 1 or len(alts) != 1 or len(alts[0]) != 1:
                    continue
                if set(record.filter.keys()) not in (set(), {"PASS"}):
                    raise ContractError(
                        f"non-PASS biallelic SNV in tree VCF: {chrom}:{record.pos}"
                    )
                positions.append(int(record.pos))
            rows.append(
                {
                    "chrom": chrom,
                    **positional_inventory(
                        positions,
                        legacy_gap_bp=legacy_gap_bp,
                        max_block_size=max_block_size,
                    ),
                }
            )
    return rows


def derive_site_catalog_inventory(
    site_catalog: Path,
    expected_chrom: str,
    legacy_gap_bp: int,
    max_block_size: int,
) -> dict[str, Any]:
    positions = []
    with gzip.open(site_catalog, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"site_index", "chrom", "pos1", "ref", "alt"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ContractError(f"site catalog missing columns: {site_catalog}")
        for expected_index, row in enumerate(reader):
            if row["chrom"] != expected_chrom:
                raise ContractError(
                    f"site catalog chromosome mismatch: {row['chrom']} != {expected_chrom}"
                )
            if int(row["site_index"]) != expected_index:
                raise ContractError("site catalog indices are not contiguous")
            positions.append(int(row["pos1"]))
    return positional_inventory(positions, legacy_gap_bp, max_block_size)


def known_inventory_check(
    inventory: Sequence[Mapping[str, Any]],
) -> None:
    if tuple(row["chrom"] for row in inventory) != AUTOSOMES:
        raise ContractError("canonical inventory check requires exact chr1-22 order")
    for row in inventory:
        expected = CANONICAL_EXPECTED[row["chrom"]]
        observed = (
            row["ssnv_sites"],
            row["k_gt8_components"],
            row["k_gt8_sites"],
            row["k_gt8_max_k"],
        )
        if observed != expected:
            raise ContractError(
                f"canonical inventory mismatch for {row['chrom']}: "
                f"expected={expected} observed={observed}"
            )


def verify_output_identities(
    outputs: Mapping[str, Any], expected_parent: Path
) -> dict[str, dict[str, Any]]:
    if not outputs:
        raise ContractError("child receipt has no outputs")
    verified: dict[str, dict[str, Any]] = {}
    expected_parent = expected_parent.resolve()
    for name, raw in sorted(outputs.items()):
        if not isinstance(raw, Mapping):
            raise ContractError(f"output identity is not an object: {name}")
        required = {"path", "size_bytes", "sha256"}
        if set(raw) != required:
            raise ContractError(f"output identity keys differ: {name}")
        path = Path(str(raw["path"]))
        require_regular_file(path, f"child output {name}")
        resolved = path.resolve()
        if resolved.parent != expected_parent:
            raise ContractError(
                f"child output escaped stage directory: {resolved} "
                f"parent={expected_parent}"
            )
        observed = file_identity(path)
        if observed != dict(raw):
            raise ContractError(
                f"child output identity mismatch: {name}; "
                f"expected={dict(raw)} observed={observed}"
            )
        verified[str(name)] = observed
    return verified


def expected_extract_command(
    args: argparse.Namespace, chrom: str, output_dir: Path
) -> list[str]:
    return [
        str(args.python),
        str(args.extractor.resolve()),
        "--manifest",
        str(args.manifest.resolve()),
        "--sample",
        DATASET,
        "--chrom",
        chrom,
        "--output-dir",
        str(output_dir.resolve()),
        "--mapq-min",
        str(args.mapq_min),
        "--baseq-min",
        str(args.baseq_min),
        "--bridge-thresholds",
        args.bridge_thresholds,
        "--samtools-threads",
        str(args.samtools_threads),
    ]


def expected_partition_command(
    args: argparse.Namespace,
    chrom: str,
    extraction_dir: Path,
    output_dir: Path,
    inventory: Mapping[str, Any],
) -> list[str]:
    prefix = f"{DATASET}.{chrom}"
    return [
        str(args.python),
        str(args.partitioner.resolve()),
        "--dataset",
        DATASET,
        "--chrom",
        chrom,
        "--site-catalog",
        str((extraction_dir / f"{prefix}.site_catalog.tsv.gz").resolve()),
        "--molecule-calls",
        str(
            (
                extraction_dir / f"{prefix}.molecule_sparse_calls.tsv.gz"
            ).resolve()
        ),
        "--output-dir",
        str(output_dir.resolve()),
        "--legacy-gap-bp",
        str(args.legacy_gap_bp),
        "--max-block-size",
        str(args.max_block_size),
        "--expected-target-components",
        str(inventory["k_gt8_components"]),
        "--expected-target-sites",
        str(inventory["k_gt8_sites"]),
    ]


def verify_extract_child(
    receipt_path: Path,
    args: argparse.Namespace,
    chrom: str,
    output_dir: Path,
    expected_inventory: Mapping[str, Any],
    expected_command: Sequence[str],
    manifest_identity: Mapping[str, Any],
    extractor_identity: Mapping[str, Any],
) -> dict[str, Any]:
    receipt_sha = verify_sha256_sidecar(receipt_path)
    receipt = strict_json_load(receipt_path)
    if not isinstance(receipt, Mapping) or receipt.get("all_pass") is not True:
        raise ContractError(f"extraction child receipt is not all_pass: {receipt_path}")
    if receipt.get("schema_name") != (
        "intersubmod.lossless_read_linkage_chromosome_receipt"
    ):
        raise ContractError("unexpected extraction child receipt schema")
    expected_scope = {
        "dataset": DATASET,
        "chrom": chrom,
        "n_sSNV": expected_inventory["ssnv_sites"],
    }
    if receipt.get("scope") != expected_scope:
        raise ContractError(
            f"extraction scope mismatch: expected={expected_scope} "
            f"observed={receipt.get('scope')}"
        )
    provenance = receipt.get("provenance")
    expected_provenance = {
        "extractor": {
            "path": str(args.extractor.resolve()),
            "sha256": extractor_identity["sha256"],
        },
        "manifest": {
            "path": str(args.manifest.resolve()),
            "sha256": manifest_identity["sha256"],
        },
    }
    if provenance != expected_provenance:
        raise ContractError("extraction provenance does not match frozen contract")
    # The child records the samtools command, not its own Python invocation.
    if not isinstance(receipt.get("command"), list):
        raise ContractError("extraction child command is missing")
    outputs = verify_output_identities(receipt.get("outputs", {}), output_dir)
    return {
        "receipt": dict(receipt),
        "receipt_identity": {
            "path": str(receipt_path.resolve()),
            "size_bytes": receipt_path.stat().st_size,
            "sha256": receipt_sha,
        },
        "outputs": outputs,
        "runner_invocation": list(expected_command),
    }


def verify_partition_child(
    receipt_path: Path,
    args: argparse.Namespace,
    chrom: str,
    output_dir: Path,
    expected_inventory: Mapping[str, Any],
    expected_command: Sequence[str],
    site_catalog: Path,
    molecule_calls: Path,
) -> dict[str, Any]:
    receipt_sha = verify_sha256_sidecar(receipt_path)
    receipt = strict_json_load(receipt_path)
    if not isinstance(receipt, Mapping) or receipt.get("all_pass") is not True:
        raise ContractError(f"partition child receipt is not all_pass: {receipt_path}")
    if receipt.get("schema_name") != "intersubmod.k_gt8_read_supported_segmentation":
        raise ContractError("unexpected partition child receipt schema")
    scope = receipt.get("scope")
    if not isinstance(scope, Mapping):
        raise ContractError("partition child scope is missing")
    expected_scope_values = {
        "dataset": DATASET,
        "chrom": chrom,
        "site_catalog_sites": expected_inventory["ssnv_sites"],
    }
    for key, value in expected_scope_values.items():
        if scope.get(key) != value:
            raise ContractError(
                f"partition scope mismatch for {key}: "
                f"expected={value} observed={scope.get(key)}"
            )
    counts = receipt.get("counts")
    if not isinstance(counts, Mapping):
        raise ContractError("partition counts are missing")
    expected_counts = {
        "target_components": expected_inventory["k_gt8_components"],
        "target_sites": expected_inventory["k_gt8_sites"],
    }
    for key, value in expected_counts.items():
        if counts.get(key, 0) != value:
            raise ContractError(
                f"partition count mismatch for {key}: "
                f"expected={value} observed={counts.get(key, 0)}"
            )
    if receipt.get("command") != list(expected_command):
        raise ContractError("partition child command differs from expected command")
    expected_inputs = {
        "site_catalog": file_identity(site_catalog),
        "molecule_calls": file_identity(molecule_calls),
    }
    if receipt.get("inputs") != expected_inputs:
        raise ContractError("partition input identities differ from extraction outputs")
    outputs = verify_output_identities(receipt.get("outputs", {}), output_dir)
    return {
        "receipt": dict(receipt),
        "receipt_identity": {
            "path": str(receipt_path.resolve()),
            "size_bytes": receipt_path.stat().st_size,
            "sha256": receipt_sha,
        },
        "outputs": outputs,
        "runner_invocation": list(expected_command),
    }


def stage_paths(chrom_root: Path, stage: str) -> dict[str, Path]:
    return {
        "output": chrom_root / stage,
        "stdout": chrom_root / f"{stage}.stdout.log",
        "stderr": chrom_root / f"{stage}.stderr.log",
        "resource": chrom_root / f"{stage}.time_v.txt",
        "stage_receipt": chrom_root / f"{stage}.stage_receipt.json",
        "failure": chrom_root / f"{stage}.failure.json",
    }


def ensure_stage_absent(paths: Mapping[str, Path]) -> None:
    occupied = [path for path in paths.values() if path.exists() or path.is_symlink()]
    if occupied:
        raise ContractError(
            "refusing to overwrite incomplete/unverified stage artifacts: "
            + ", ".join(map(str, occupied))
        )


def run_timed_stage(
    *,
    stage: str,
    chrom: str,
    command: Sequence[str],
    paths: Mapping[str, Path],
    contract_sha256: str,
    tool_identity: Mapping[str, Any],
    validator: Callable[[], dict[str, Any]],
) -> dict[str, Any]:
    ensure_stage_absent(paths)
    paths["stdout"].parent.mkdir(parents=True, exist_ok=True)
    wrapped = [
        str(TIME_BINARY),
        "-v",
        "-o",
        str(paths["resource"]),
        "--",
        *command,
    ]
    started = time.perf_counter()
    with paths["stdout"].open("xb") as stdout, paths["stderr"].open("xb") as stderr:
        completed = subprocess.run(wrapped, stdout=stdout, stderr=stderr, check=False)
    wall_seconds = time.perf_counter() - started
    if completed.returncode != 0:
        failure = {
            "schema_name": f"{SCHEMA_NAME}.stage_failure",
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": utc_now(),
            "all_pass": False,
            "sample": DATASET,
            "chrom": chrom,
            "stage": stage,
            "exit_code": completed.returncode,
            "wall_seconds": wall_seconds,
            "command": list(command),
            "timed_command": wrapped,
            "contract_sha256": contract_sha256,
            "tool": dict(tool_identity),
            "logs": {
                name: file_identity(paths[name])
                for name in ("stdout", "stderr", "resource")
                if paths[name].is_file()
            },
        }
        write_json_atomic(paths["failure"], failure)
        write_sha256_sidecar(paths["failure"])
        raise ContractError(
            f"{stage} failed for {chrom} with exit={completed.returncode}; "
            f"see {paths['failure']}"
        )
    child = validator()
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.stage_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "stage": stage,
        "status": "COMPLETED",
        "exit_code": completed.returncode,
        "wall_seconds": wall_seconds,
        "command": list(command),
        "timed_command": wrapped,
        "scheduling_policy": "inherited_from_outer_process",
        "contract_sha256": contract_sha256,
        "tool": dict(tool_identity),
        "child_receipt": child["receipt_identity"],
        "logs": {
            name: file_identity(paths[name])
            for name in ("stdout", "stderr", "resource")
        },
    }
    write_json_atomic(paths["stage_receipt"], receipt)
    receipt_sha = sha256_path(paths["stage_receipt"])
    write_sha256_sidecar(paths["stage_receipt"])
    return {
        "stage_receipt": receipt,
        "stage_receipt_identity": {
            "path": str(paths["stage_receipt"].resolve()),
            "size_bytes": paths["stage_receipt"].stat().st_size,
            "sha256": receipt_sha,
        },
        "child": child,
        "resumed": False,
    }


def verify_stage_receipt(
    *,
    stage: str,
    chrom: str,
    command: Sequence[str],
    paths: Mapping[str, Path],
    contract_sha256: str,
    tool_identity: Mapping[str, Any],
    validator: Callable[[], dict[str, Any]],
) -> dict[str, Any]:
    if not paths["stage_receipt"].is_file():
        raise ContractError(
            f"existing {stage} output lacks reusable stage receipt: "
            f"{paths['stage_receipt']}"
        )
    receipt_sha = verify_sha256_sidecar(paths["stage_receipt"])
    receipt = strict_json_load(paths["stage_receipt"])
    expected_scalars = {
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "stage": stage,
        "status": "COMPLETED",
        "exit_code": 0,
        "contract_sha256": contract_sha256,
    }
    for key, expected in expected_scalars.items():
        if receipt.get(key) != expected:
            raise ContractError(
                f"stage receipt mismatch for {stage} {chrom} {key}: "
                f"expected={expected} observed={receipt.get(key)}"
            )
    if receipt.get("command") != list(command):
        raise ContractError(f"stage command drift for {stage} {chrom}")
    if receipt.get("tool") != dict(tool_identity):
        raise ContractError(f"stage tool drift for {stage} {chrom}")
    for name in ("stdout", "stderr", "resource"):
        expected = receipt.get("logs", {}).get(name)
        if expected != file_identity(paths[name]):
            raise ContractError(f"stage log identity drift: {stage} {chrom} {name}")
    child = validator()
    if receipt.get("child_receipt") != child["receipt_identity"]:
        raise ContractError(f"child receipt identity drift for {stage} {chrom}")
    return {
        "stage_receipt": dict(receipt),
        "stage_receipt_identity": {
            "path": str(paths["stage_receipt"].resolve()),
            "size_bytes": paths["stage_receipt"].stat().st_size,
            "sha256": receipt_sha,
        },
        "child": child,
        "resumed": True,
    }


def get_or_run_stage(
    *,
    resume: bool,
    stage: str,
    chrom: str,
    command: Sequence[str],
    paths: Mapping[str, Path],
    contract_sha256: str,
    tool_identity: Mapping[str, Any],
    validator: Callable[[], dict[str, Any]],
) -> dict[str, Any]:
    any_present = any(
        path.exists() or path.is_symlink() for path in paths.values()
    )
    if any_present:
        if not resume:
            raise ContractError(f"{stage} artifacts already exist without --resume")
        if paths["failure"].exists():
            raise ContractError(
                f"{stage} has a prior failure receipt; archive it before retry: "
                f"{paths['failure']}"
            )
        return verify_stage_receipt(
            stage=stage,
            chrom=chrom,
            command=command,
            paths=paths,
            contract_sha256=contract_sha256,
            tool_identity=tool_identity,
            validator=validator,
        )
    return run_timed_stage(
        stage=stage,
        chrom=chrom,
        command=command,
        paths=paths,
        contract_sha256=contract_sha256,
        tool_identity=tool_identity,
        validator=validator,
    )


def skip_partition(
    *,
    chrom: str,
    path: Path,
    extraction_receipt_identity: Mapping[str, Any],
    contract_sha256: str,
    inventory: Mapping[str, Any],
    resume: bool,
) -> dict[str, Any]:
    expected_stable = {
        "schema_name": f"{SCHEMA_NAME}.partition_skip_receipt",
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "stage": "partition",
        "status": "SKIP_NO_K_GT8_TARGET",
        "reason": (
            "authenticated tree VCF and extracted site catalog both contain "
            "zero legacy positional components with k>8"
        ),
        "contract_sha256": contract_sha256,
        "expected_target_components": 0,
        "expected_target_sites": 0,
        "extraction_receipt": dict(extraction_receipt_identity),
        "inventory": dict(inventory),
    }
    if path.exists() or path.is_symlink():
        if not resume:
            raise ContractError(f"partition skip receipt exists without --resume: {path}")
        receipt_sha = verify_sha256_sidecar(path)
        observed = strict_json_load(path)
        stable = {key: observed.get(key) for key in expected_stable}
        if stable != expected_stable:
            raise ContractError(f"partition skip receipt drift: {path}")
        return {
            "stage_receipt": dict(observed),
            "stage_receipt_identity": {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": receipt_sha,
            },
            "child": None,
            "resumed": True,
        }
    if resume:
        partition_dir = path.parent / "partition"
        if partition_dir.exists():
            raise ContractError(
                f"unexpected partition output for zero-target chromosome: {partition_dir}"
            )
    receipt = {**expected_stable, "created_at_utc": utc_now()}
    write_json_atomic(path, receipt)
    receipt_sha = sha256_path(path)
    write_sha256_sidecar(path)
    return {
        "stage_receipt": receipt,
        "stage_receipt_identity": {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": receipt_sha,
        },
        "child": None,
        "resumed": False,
    }


def establish_summary_tsv(
    path: Path, rows: Sequence[Mapping[str, Any]], *, resume: bool
) -> None:
    fieldnames = (
        "chrom",
        "ssnv_sites",
        "k_gt8_components",
        "k_gt8_sites",
        "k_gt8_max_k",
        "extraction_status",
        "extraction_wall_seconds",
        "extraction_stage_receipt",
        "extraction_stage_receipt_sha256",
        "partition_status",
        "partition_wall_seconds",
        "partition_stage_receipt",
        "partition_stage_receipt_sha256",
    )
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    payload = buffer.getvalue().encode("utf-8")
    if path.exists() or path.is_symlink():
        if not resume:
            raise ContractError(f"summary exists without --resume: {path}")
        verify_sha256_sidecar(path)
        if path.read_bytes() != payload:
            raise ContractError(f"existing summary differs on resume: {path}")
        return
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    write_sha256_sidecar(path)


def parse_chromosomes(value: str) -> tuple[str, ...]:
    chromosomes = tuple(token.strip() for token in value.split(",") if token.strip())
    if not chromosomes or len(chromosomes) != len(set(chromosomes)):
        raise argparse.ArgumentTypeError("chromosomes must be nonempty and unique")
    invalid = [chrom for chrom in chromosomes if chrom not in AUTOSOMES]
    if invalid:
        raise argparse.ArgumentTypeError(f"invalid autosomes: {invalid}")
    return chromosomes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--expected-manifest-sha256", default=DEFAULT_MANIFEST_SHA256
    )
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--extractor", type=Path, default=DEFAULT_EXTRACTOR)
    parser.add_argument("--partitioner", type=Path, default=DEFAULT_PARTITIONER)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--mapq-min", type=int, default=20)
    parser.add_argument("--baseq-min", type=int, default=20)
    parser.add_argument("--bridge-thresholds", default="1,2,3,5")
    parser.add_argument("--samtools-threads", type=int, default=1)
    parser.add_argument("--legacy-gap-bp", type=int, default=50_000)
    parser.add_argument("--max-block-size", type=int, default=8)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--test-mode",
        action="store_true",
        help=(
            "permit a chromosome subset or noncanonical frozen manifest; writes "
            "_TEST_SUCCESS, never the comprehensive _SUCCESS marker"
        ),
    )
    parser.add_argument(
        "--chromosomes",
        type=parse_chromosomes,
        default=AUTOSOMES,
        help="comma-separated autosomes; subset is allowed only with --test-mode",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.test_mode and args.chromosomes != AUTOSOMES:
        raise ContractError("production scope must be exact chr1-22 in canonical order")
    if (
        not args.test_mode
        and args.expected_manifest_sha256 != DEFAULT_MANIFEST_SHA256
    ):
        raise ContractError("production scope requires the canonical frozen manifest SHA")
    if args.mapq_min < 0 or args.baseq_min < 0:
        raise ContractError("quality thresholds must be nonnegative")
    if args.samtools_threads < 0:
        raise ContractError("samtools threads must be nonnegative")
    if args.legacy_gap_bp < 1 or args.max_block_size < 1:
        raise ContractError("gap and block size must be positive")
    thresholds = tuple(
        sorted({int(token) for token in args.bridge_thresholds.split(",")})
    )
    if not thresholds or thresholds[0] < 1:
        raise ContractError("bridge thresholds must be positive integers")
    if args.output_root.is_symlink():
        raise ContractError("output root must not be a symlink")
    require_regular_file(args.manifest, "frozen manifest")
    require_regular_file(args.extractor, "extractor")
    require_regular_file(args.partitioner, "partitioner")
    require_regular_file(args.python, "Python executable")
    require_regular_file(TIME_BINARY, "/usr/bin/time")


def build_contract(args: argparse.Namespace) -> tuple[dict[str, Any], Mapping[str, Any]]:
    manifest_identity = file_identity(args.manifest)
    if manifest_identity["sha256"] != args.expected_manifest_sha256:
        raise ContractError(
            "frozen manifest SHA mismatch: "
            f"expected={args.expected_manifest_sha256} "
            f"observed={manifest_identity['sha256']}"
        )
    manifest = strict_json_load(args.manifest)
    if not isinstance(manifest, Mapping):
        raise ContractError("frozen manifest root must be an object")
    if manifest.get("schema_name") != "intersubmod.layered_input_manifest":
        raise ContractError("unexpected frozen manifest schema")
    entry = sample_entry(manifest)
    authorities = verify_consumed_authorities(entry)
    tree_vcf = Path(str(authorities["tree_vcf"]["path"]))
    inventory = derive_vcf_inventory(
        tree_vcf,
        args.chromosomes,
        legacy_gap_bp=args.legacy_gap_bp,
        max_block_size=args.max_block_size,
    )
    if not args.test_mode:
        known_inventory_check(inventory)
    contract = {
        "schema_name": f"{SCHEMA_NAME}.run_contract",
        "schema_version": SCHEMA_VERSION,
        "sample": DATASET,
        "scope": {
            "chromosomes": list(args.chromosomes),
            "comprehensive_chr1_22": (
                not args.test_mode and args.chromosomes == AUTOSOMES
            ),
            "test_mode": args.test_mode,
        },
        "frozen_manifest": manifest_identity,
        "manifest_id": manifest.get("manifest_id"),
        "consumed_input_authorities": authorities,
        "tools": {
            "runner": file_identity(Path(__file__).resolve()),
            "extractor": file_identity(args.extractor),
            "partitioner": file_identity(args.partitioner),
            "python": file_identity(args.python),
            "time": file_identity(TIME_BINARY),
        },
        "parameters": {
            "mapq_min": args.mapq_min,
            "baseq_min": args.baseq_min,
            "bridge_thresholds": args.bridge_thresholds,
            "samtools_threads": args.samtools_threads,
            "legacy_gap_bp": args.legacy_gap_bp,
            "max_block_size": args.max_block_size,
            "sequential_chromosomes": True,
            "scheduling_policy": "inherited; apply nice/ionice externally",
            "zero_target_policy": (
                "always extract; skip partition only after VCF and extracted "
                "site-catalog inventories both confirm k>8 target count zero"
            ),
        },
        "vcf_derived_inventory": inventory,
        "canonical_expected_inventory_enforced": not args.test_mode,
    }
    return contract, entry


def establish_run_contract(
    output_root: Path, contract: Mapping[str, Any], resume: bool
) -> tuple[Path, str]:
    contract_path = output_root / "run_contract.json"
    if output_root.exists():
        if not resume:
            raise ContractError(
                f"output root exists; use --resume only for a verified run: {output_root}"
            )
        if not output_root.is_dir() or output_root.is_symlink():
            raise ContractError(f"output root is not a real directory: {output_root}")
        observed_sha = verify_sha256_sidecar(contract_path)
        observed = strict_json_load(contract_path)
        if observed != dict(contract):
            raise ContractError("resume contract differs from current frozen contract")
        return contract_path, observed_sha
    if resume:
        raise ContractError(f"--resume output root does not exist: {output_root}")
    output_root.mkdir(parents=True, exist_ok=False)
    write_json_atomic(contract_path, contract)
    contract_sha = sha256_path(contract_path)
    write_sha256_sidecar(contract_path)
    return contract_path, contract_sha


def run(args: argparse.Namespace) -> dict[str, Any]:
    validate_args(args)
    contract, _entry = build_contract(args)
    contract_path, contract_sha = establish_run_contract(
        args.output_root, contract, args.resume
    )
    inventory_by_chrom = {
        row["chrom"]: row for row in contract["vcf_derived_inventory"]
    }
    tools = contract["tools"]
    chromosome_rows = []
    chromosome_receipts = []

    for chrom in args.chromosomes:
        inventory = inventory_by_chrom[chrom]
        chrom_root = args.output_root / "chromosomes" / chrom
        chrom_root.mkdir(parents=True, exist_ok=True)

        extract_paths = stage_paths(chrom_root, "extract")
        extract_command = expected_extract_command(
            args, chrom, extract_paths["output"]
        )
        extract_receipt_path = extract_paths["output"] / "receipt.json"
        extract_validator = lambda: verify_extract_child(
            extract_receipt_path,
            args,
            chrom,
            extract_paths["output"],
            inventory,
            extract_command,
            contract["frozen_manifest"],
            tools["extractor"],
        )
        extraction = get_or_run_stage(
            resume=args.resume,
            stage="extract",
            chrom=chrom,
            command=extract_command,
            paths=extract_paths,
            contract_sha256=contract_sha,
            tool_identity=tools["extractor"],
            validator=extract_validator,
        )

        prefix = f"{DATASET}.{chrom}"
        site_catalog = (
            extract_paths["output"] / f"{prefix}.site_catalog.tsv.gz"
        )
        molecule_calls = (
            extract_paths["output"] / f"{prefix}.molecule_sparse_calls.tsv.gz"
        )
        extracted_inventory = derive_site_catalog_inventory(
            site_catalog,
            expected_chrom=chrom,
            legacy_gap_bp=args.legacy_gap_bp,
            max_block_size=args.max_block_size,
        )
        if extracted_inventory != {
            key: inventory[key]
            for key in (
                "ssnv_sites",
                "positional_components_all",
                "k_gt8_components",
                "k_gt8_sites",
                "k_gt8_max_k",
                "positions_sha256",
            )
        }:
            raise ContractError(
                f"extracted site inventory differs from frozen VCF for {chrom}"
            )

        if inventory["k_gt8_components"] == 0:
            partition = skip_partition(
                chrom=chrom,
                path=chrom_root / "partition.skip_receipt.json",
                extraction_receipt_identity=extraction[
                    "stage_receipt_identity"
                ],
                contract_sha256=contract_sha,
                inventory=inventory,
                resume=args.resume,
            )
        else:
            partition_paths = stage_paths(chrom_root, "partition")
            partition_command = expected_partition_command(
                args,
                chrom,
                extract_paths["output"],
                partition_paths["output"],
                inventory,
            )
            partition_receipt_path = partition_paths["output"] / "receipt.json"
            partition_validator = lambda: verify_partition_child(
                partition_receipt_path,
                args,
                chrom,
                partition_paths["output"],
                inventory,
                partition_command,
                site_catalog,
                molecule_calls,
            )
            partition = get_or_run_stage(
                resume=args.resume,
                stage="partition",
                chrom=chrom,
                command=partition_command,
                paths=partition_paths,
                contract_sha256=contract_sha,
                tool_identity=tools["partitioner"],
                validator=partition_validator,
            )

        partition_status = partition["stage_receipt"]["status"]
        chromosome_row = {
            "chrom": chrom,
            "ssnv_sites": inventory["ssnv_sites"],
            "k_gt8_components": inventory["k_gt8_components"],
            "k_gt8_sites": inventory["k_gt8_sites"],
            "k_gt8_max_k": inventory["k_gt8_max_k"],
            "extraction_status": extraction["stage_receipt"]["status"],
            "extraction_wall_seconds": extraction["stage_receipt"].get(
                "wall_seconds", ""
            ),
            "extraction_stage_receipt": extraction[
                "stage_receipt_identity"
            ]["path"],
            "extraction_stage_receipt_sha256": extraction[
                "stage_receipt_identity"
            ]["sha256"],
            "partition_status": partition_status,
            "partition_wall_seconds": partition["stage_receipt"].get(
                "wall_seconds", ""
            ),
            "partition_stage_receipt": partition[
                "stage_receipt_identity"
            ]["path"],
            "partition_stage_receipt_sha256": partition[
                "stage_receipt_identity"
            ]["sha256"],
        }
        chromosome_rows.append(chromosome_row)
        chromosome_receipts.append(
            {
                "chrom": chrom,
                "inventory": inventory,
                "extraction": extraction["stage_receipt_identity"],
                "partition": partition["stage_receipt_identity"],
                "partition_status": partition_status,
            }
        )

    summary_path = args.output_root / "chromosome_summary.tsv"
    establish_summary_tsv(summary_path, chromosome_rows, resume=args.resume)
    totals = {
        "chromosomes": len(chromosome_rows),
        "ssnv_sites": sum(row["ssnv_sites"] for row in chromosome_rows),
        "k_gt8_components": sum(
            row["k_gt8_components"] for row in chromosome_rows
        ),
        "k_gt8_sites": sum(row["k_gt8_sites"] for row in chromosome_rows),
        "k_gt8_max_k": max(
            (row["k_gt8_max_k"] for row in chromosome_rows), default=0
        ),
        "partitioned_chromosomes": sum(
            row["partition_status"] == "COMPLETED"
            for row in chromosome_rows
        ),
        "zero_target_skipped_chromosomes": sum(
            row["partition_status"] == "SKIP_NO_K_GT8_TARGET"
            for row in chromosome_rows
        ),
    }
    if not args.test_mode and totals != {
        "chromosomes": 22,
        "ssnv_sites": 79687,
        "k_gt8_components": 408,
        "k_gt8_sites": 47570,
        "k_gt8_max_k": 3574,
        "partitioned_chromosomes": 21,
        "zero_target_skipped_chromosomes": 1,
    }:
        raise ContractError(f"canonical final totals mismatch: {totals}")

    receipt = {
        "schema_name": f"{SCHEMA_NAME}.run_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "comprehensive_all_pass": (
            not args.test_mode and args.chromosomes == AUTOSOMES
        ),
        "sample": DATASET,
        "scope": {
            "chromosomes": list(args.chromosomes),
            "test_mode": args.test_mode,
        },
        "contract": {
            "path": str(contract_path.resolve()),
            "size_bytes": contract_path.stat().st_size,
            "sha256": contract_sha,
        },
        "command": [str(Path(__file__).resolve()), *sys.argv[1:]],
        "scheduling_policy": (
            "sequential stages; nice/ionice inherited from outer process"
        ),
        "totals": totals,
        "chromosomes": chromosome_receipts,
        "outputs": {
            "chromosome_summary": file_identity(summary_path),
        },
        "claim_ceiling": (
            "complete engineering execution and conservation for requested "
            "scope; not proof of a unique true evolutionary tree"
        ),
    }
    receipt_path = args.output_root / "receipt.json"
    if receipt_path.exists() or receipt_path.is_symlink():
        if not args.resume:
            raise ContractError(f"final receipt exists without --resume: {receipt_path}")
        previous_sha = verify_sha256_sidecar(receipt_path)
        previous = strict_json_load(receipt_path)
        stable_keys = (
            "schema_name",
            "schema_version",
            "all_pass",
            "comprehensive_all_pass",
            "sample",
            "scope",
            "contract",
            "scheduling_policy",
            "totals",
            "chromosomes",
            "outputs",
            "claim_ceiling",
        )
        if {key: previous.get(key) for key in stable_keys} != {
            key: receipt.get(key) for key in stable_keys
        }:
            raise ContractError("existing final receipt differs on resume")
        final_receipt = dict(previous)
        final_receipt_sha = previous_sha
    else:
        write_json_atomic(receipt_path, receipt)
        final_receipt_sha = sha256_path(receipt_path)
        write_sha256_sidecar(receipt_path)
        final_receipt = receipt

    marker_name = "_TEST_SUCCESS" if args.test_mode else "_SUCCESS"
    marker_path = args.output_root / marker_name
    expected_marker = {
        "schema_name": f"{SCHEMA_NAME}.success_marker",
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "comprehensive": not args.test_mode,
        "sample": DATASET,
        "scope": list(args.chromosomes),
        "run_receipt": {
            "path": str(receipt_path.resolve()),
            "sha256": final_receipt_sha,
        },
    }
    if marker_path.exists() or marker_path.is_symlink():
        observed = strict_json_load(marker_path)
        stable = {key: observed.get(key) for key in expected_marker}
        if stable != expected_marker:
            raise ContractError(f"success marker drift: {marker_path}")
    else:
        write_json_atomic(
            marker_path, {**expected_marker, "created_at_utc": utc_now()}
        )
    print(
        json.dumps(
            {
                "all_pass": True,
                "comprehensive_all_pass": final_receipt[
                    "comprehensive_all_pass"
                ],
                "totals": totals,
                "receipt": str(receipt_path),
                "receipt_sha256": final_receipt_sha,
                "marker": str(marker_path),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return final_receipt


def main() -> None:
    try:
        args = parse_args()
        run(args)
    except (ContractError, OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"FAIL-CLOSED: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
