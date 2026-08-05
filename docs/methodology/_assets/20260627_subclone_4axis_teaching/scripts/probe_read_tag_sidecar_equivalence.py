#!/usr/bin/env python3
"""Run a bounded real-data equivalence probe for direct BAM tags versus a tag sidecar."""

import argparse
import hashlib
import importlib.util
import json
import os
import sys
from collections import Counter
from pathlib import Path

import pysam


SCRIPT_DIR = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("sm_linkage_genomewide", SCRIPT_DIR / "sm_linkage_genomewide.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def canonical_sha256(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def storage_identity(path, chunk_size=1024 * 1024):
    """Return an explicitly lower-assurance identity without scanning a full BAM."""
    stat = path.stat()
    offsets = sorted({0, max(0, stat.st_size // 2 - chunk_size // 2), max(0, stat.st_size - chunk_size)})
    digest = hashlib.sha256()
    sampled = []
    with path.open("rb") as handle:
        for offset in offsets:
            handle.seek(offset)
            payload = handle.read(chunk_size)
            digest.update(f"{offset}:{len(payload)}\n".encode("ascii"))
            digest.update(payload)
            sampled.append({"offset": offset, "length": len(payload)})
    index = Path(str(path) + ".bai")
    return {
        "policy": "storage_identity_v1_not_full_content_hash",
        "path": str(path),
        "realpath": os.path.realpath(path),
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "sampled_chunks": sampled,
        "sampled_content_sha256": digest.hexdigest(),
        "index_path": str(index),
        "index_sha256": file_sha256(index),
    }


def alignment_key_text(key):
    return "\t".join(str(value) for value in key)


def load_group(vcf_path, chrom, positions):
    wanted = set(positions)
    records = {}
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf.fetch(chrom, min(wanted) - 1, max(wanted)):
            if record.pos not in wanted:
                continue
            if len(record.ref) != 1 or len(record.alts or ()) != 1 or len(record.alts[0]) != 1:
                raise RuntimeError(f"non-biallelic SNV at {chrom}:{record.pos}")
            records[record.pos] = (record.pos, record.ref.upper(), record.alts[0].upper(), "PROBE")
    missing = sorted(wanted - records.keys())
    if missing:
        raise RuntimeError(f"probe positions missing from VCF: {missing}")
    return [records[position] for position in sorted(records)]


def build_sidecar(tagged_bam, chrom, start, end, output_gz):
    plain = Path(str(output_gz).removesuffix(".gz"))
    seen = {}
    counts = Counter()
    with pysam.AlignmentFile(str(tagged_bam), "rb") as bam, plain.open("w", encoding="utf-8") as handle:
        handle.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
        for record in bam.fetch(chrom, start - 1, end):
            if record.is_unmapped:
                continue
            key = M._alignment_key(record)
            value = (str(record.get_tag("HP")) if record.has_tag("HP") else ".",
                     str(record.get_tag("PS")) if record.has_tag("PS") else ".")
            if key in seen:
                counts["duplicate_rows"] += 1
                counts["duplicate_conflicts"] += int(seen[key] != value)
            else:
                seen[key] = value
            fields = (
                record.reference_name,
                str(record.reference_start),
                str(record.reference_end or record.reference_start + 1),
                record.query_name,
                str(record.flag),
                str(record.mapping_quality),
                M._cigar_digest(record),
                value[0],
                value[1],
            )
            handle.write("\t".join(fields) + "\n")
            counts["rows"] += 1
            counts[f"HP_{value[0]}"] += 1
    if counts["duplicate_rows"]:
        raise RuntimeError(
            f"coordinate_join_v1 is not unique: duplicates={counts['duplicate_rows']} "
            f"conflicts={counts['duplicate_conflicts']}"
        )
    pysam.tabix_compress(str(plain), str(output_gz), force=False)
    pysam.tabix_index(str(output_gz), seq_col=0, start_col=1, end_col=2, zerobased=True, force=False)
    return plain, counts


def reset_sidecar(path=""):
    if M._TAG_TABIX is not None:
        M._TAG_TABIX.close()
    M._TAG_TABIX = None
    M.READ_TAG_SIDECAR = str(path) if path else ""


def observe(bam_path, chrom, group, sidecar=None):
    reset_sidecar(sidecar)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        calls, hp, ps = M.per_read_calls(bam, chrom, group)
    calls_out = {
        alignment_key_text(key): {str(position): call for position, call in sorted(values.items())}
        for key, values in sorted(calls.items(), key=lambda item: alignment_key_text(item[0]))
    }
    hp_out = {alignment_key_text(key): str(value) for key, value in sorted(hp.items(), key=lambda item: alignment_key_text(item[0]))}
    ps_out = {alignment_key_text(key): str(value) for key, value in sorted(ps.items(), key=lambda item: alignment_key_text(item[0]))}
    family_counts = Counter(
        value.split("-", 1)[0] if value not in {".", ""} else "none" for value in hp_out.values()
    )
    payload = {"calls": calls_out, "hp": hp_out, "ps": ps_out, "family_counts": dict(family_counts)}
    return {
        "payload": payload,
        "payload_sha256": canonical_sha256(payload),
        "calls_sha256": canonical_sha256(calls_out),
        "hp_sha256": canonical_sha256(hp_out),
        "ps_sha256": canonical_sha256(ps_out),
        "diagnostics": M.last_tag_diagnostics(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tagged-bam", required=True, type=Path)
    parser.add_argument("--raw-bam", required=True, type=Path)
    parser.add_argument("--somatic-vcf", required=True, type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--positions", required=True, help="Comma-separated 1-based SNV positions")
    parser.add_argument("--mapq-min", type=int, default=20)
    parser.add_argument("--baseq-min", type=int, default=0)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.chrom != "chr1":
        raise SystemExit("this pre-registered probe currently requires --chrom chr1")
    if args.output_dir.exists():
        raise SystemExit(f"immutable output already exists: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    positions = sorted({int(value) for value in args.positions.split(",")})
    if len(positions) < 2:
        raise SystemExit("at least two positions are required")
    M.MAPQ_MIN = args.mapq_min
    M.BASEQ_MIN = args.baseq_min
    group = load_group(args.somatic_vcf, args.chrom, positions)
    sidecar = args.output_dir / "COLO829.coordinate_join_v1.read_tags.tsv.gz"
    plain, sidecar_counts = build_sidecar(
        args.tagged_bam, args.chrom, min(positions), max(positions), sidecar
    )
    direct = observe(args.tagged_bam, args.chrom, group)
    joined = observe(args.raw_bam, args.chrom, group, sidecar)
    comparisons = {
        "payload_equal": direct["payload"] == joined["payload"],
        "payload_sha256_equal": direct["payload_sha256"] == joined["payload_sha256"],
        "calls_sha256_equal": direct["calls_sha256"] == joined["calls_sha256"],
        "hp_sha256_equal": direct["hp_sha256"] == joined["hp_sha256"],
        "ps_sha256_equal": direct["ps_sha256"] == joined["ps_sha256"],
        "sidecar_exact_matches_all_exposures": (
            joined["diagnostics"]["sidecar_exact_matches"]
            == joined["diagnostics"]["alignment_group_exposures"]
        ),
        "sidecar_missing_zero": joined["diagnostics"]["sidecar_missing"] == 0,
        "sidecar_conflicts_zero": joined["diagnostics"]["sidecar_conflicts"] == 0,
    }
    output = {
        "schema_version": "1.0",
        "scope": "PARTIAL_BOUNDED_REAL_DATA_PROBE",
        "claim_limit": "join equivalence only; historical tags are not clean production-tag evidence",
        "identity_schema": "coordinate_join_v1",
        "command_argv": sys.argv,
        "inputs": {
            "tagged_bam": storage_identity(args.tagged_bam),
            "raw_bam": storage_identity(args.raw_bam),
            "somatic_vcf": {"path": str(args.somatic_vcf), "sha256": file_sha256(args.somatic_vcf)},
        },
        "region": {"chrom": args.chrom, "start": min(positions), "end": max(positions), "positions": positions},
        "params": {"MAPQ_MIN": args.mapq_min, "BASEQ_MIN": args.baseq_min},
        "sidecar": {
            "path": str(sidecar),
            "sha256": file_sha256(sidecar),
            "index_path": str(Path(str(sidecar) + ".tbi")),
            "index_sha256": file_sha256(Path(str(sidecar) + ".tbi")),
            "plain_path": str(plain),
            "counts": dict(sidecar_counts),
        },
        "direct": direct,
        "joined": joined,
        "comparisons": comparisons,
        "pass": all(comparisons.values()),
    }
    result_path = args.output_dir / "equivalence_probe.json"
    pending = args.output_dir / "equivalence_probe.json.pending"
    pending.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    pending.replace(result_path)
    reset_sidecar()
    print(f"REAL-DATA SIDECAR PROBE: {'PASS' if output['pass'] else 'FAIL'} -> {result_path}")
    print(json.dumps(comparisons, sort_keys=True))
    raise SystemExit(0 if output["pass"] else 1)


if __name__ == "__main__":
    main()
