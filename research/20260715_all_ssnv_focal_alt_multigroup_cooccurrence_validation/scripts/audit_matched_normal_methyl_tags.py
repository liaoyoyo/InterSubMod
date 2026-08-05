#!/usr/bin/env python3
"""Audit matched-normal BAM availability and MM/ML tag coverage."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pysam


NORMAL_BAMS = {
    "HCC1395": Path("/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam"),
    "HCC1395_DORADO": Path("/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam"),
    "COLO829": Path("/big8_disk/data/COLO829/ONT_PAO/PAO33946.bam"),
    "H1437": Path("/big8_disk/data/H1437/ONT/H1437BL.bam"),
    "H2009": Path("/big8_disk/data/H2009/ONT/H2009BL.bam"),
    "HCC1937": Path("/big8_disk/data/HCC1937/ONT/HCC1937BL.bam"),
    "HCC1954": Path("/big8_disk/data/HCC1954/ONT/HCC1954BL.bam"),
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def header_sha256(alignment: pysam.AlignmentFile) -> str:
    payload = alignment.text.encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def find_index(path: Path) -> Path | None:
    candidates = [Path(str(path) + ".bai"), path.with_suffix(".bai"), Path(str(path) + ".csi")]
    return next((candidate for candidate in candidates if candidate.exists()), None)


def has_any_tag(read: pysam.AlignedSegment, names: tuple[str, ...]) -> bool:
    return any(read.has_tag(name) for name in names)


def summarize_mm_tag(value: str) -> list[str]:
    return [block.split(",", 1)[0] for block in value.rstrip(";").split(";") if block]


def audit_bam(sample: str, path: Path, target_reads: int, max_scanned: int) -> dict[str, Any]:
    index = find_index(path)
    base = {
        "sample": sample,
        "bam": str(path),
        "bam_exists": path.exists(),
        "bam_size_bytes": path.stat().st_size if path.exists() else None,
        "bam_mtime_ns": path.stat().st_mtime_ns if path.exists() else None,
        "index": str(index) if index else None,
        "index_exists": index is not None,
    }
    if not path.exists() or index is None:
        return {
            **base,
            "readable": False,
            "normal_control_eligible": False,
            "error": "missing_bam_or_index",
        }

    tag_codes: Counter[str] = Counter()
    hp_values: Counter[str] = Counter()
    sampled = 0
    scanned = 0
    mm = 0
    ml = 0
    mm_ml = 0
    primary_mapped = 0
    with_sequence = 0
    try:
        with pysam.AlignmentFile(str(path), "rb") as alignment:
            references = list(alignment.references)
            preferred = [f"chr{chrom}" for chrom in (1, 5, 10, 15, 20)]
            contigs = [contig for contig in preferred if contig in references]
            if not contigs:
                contigs = references[:5]
            for contig in contigs:
                for read in alignment.fetch(contig):
                    scanned += 1
                    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                        if scanned >= max_scanned:
                            break
                        continue
                    primary_mapped += 1
                    if read.query_sequence:
                        with_sequence += 1
                    has_mm = has_any_tag(read, ("MM", "Mm"))
                    has_ml = has_any_tag(read, ("ML", "Ml"))
                    mm += int(has_mm)
                    ml += int(has_ml)
                    mm_ml += int(has_mm and has_ml)
                    if has_mm:
                        name = "MM" if read.has_tag("MM") else "Mm"
                        tag_codes.update(summarize_mm_tag(str(read.get_tag(name))))
                    if read.has_tag("HP"):
                        hp_values[str(read.get_tag("HP"))] += 1
                    sampled += 1
                    if sampled >= target_reads or scanned >= max_scanned:
                        break
                if sampled >= target_reads or scanned >= max_scanned:
                    break
            fraction = mm_ml / sampled if sampled else 0.0
            eligible = sampled >= min(100, target_reads) and fraction >= 0.95
            return {
                **base,
                "readable": True,
                "header_sha256": header_sha256(alignment),
                "references_n": len(references),
                "sampled_contigs": contigs,
                "reads_scanned": scanned,
                "primary_mapped_seen": primary_mapped,
                "primary_reads_sampled": sampled,
                "reads_with_sequence": with_sequence,
                "mm_tag_reads": mm,
                "ml_tag_reads": ml,
                "mm_ml_pair_reads": mm_ml,
                "mm_ml_pair_fraction": fraction,
                "mm_code_counts": dict(sorted(tag_codes.items())),
                "hp_value_counts": dict(sorted(hp_values.items())),
                "normal_control_eligible": eligible,
                "eligibility_rule": "primary_reads_sampled>=100 and mm_ml_pair_fraction>=0.95",
                "error": None,
            }
    except Exception as error:  # pragma: no cover - real-file failure receipt
        return {
            **base,
            "readable": False,
            "normal_control_eligible": False,
            "error": repr(error),
        }


def main() -> None:
    topic_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-reads", type=int, default=1000)
    parser.add_argument("--max-scanned", type=int, default=100_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=topic_root / "results" / "matched_normal_methyl_tag_audit.json",
    )
    args = parser.parse_args()
    if args.target_reads < 100 or args.max_scanned < args.target_reads:
        raise ValueError("Require target_reads>=100 and max_scanned>=target_reads")
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit receipt: {args.output}")

    rows = [
        audit_bam(sample, path, args.target_reads, args.max_scanned)
        for sample, path in NORMAL_BAMS.items()
    ]
    receipt = {
        "schema_name": "intersubmod.matched_normal_methyl_tag_audit",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "command": sys.argv,
        "source_code": {
            "path": str(Path(__file__).resolve()),
            "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "sampling": {
            "target_primary_reads_per_bam": args.target_reads,
            "max_records_scanned_per_bam": args.max_scanned,
            "contigs": ["chr1", "chr5", "chr10", "chr15", "chr20"],
            "exclusions": ["unmapped", "secondary", "supplementary", "duplicate"],
        },
        "samples": rows,
        "totals": {
            "n_samples": len(rows),
            "n_readable": sum(row["readable"] for row in rows),
            "n_normal_control_eligible": sum(row["normal_control_eligible"] for row in rows),
        },
        "pass": len(rows) == 7 and all(row["readable"] for row in rows),
        "interpretation": (
            "MM/ML sampling verifies that matched-normal methylation evidence can be extracted. "
            "It does not validate biological normality, phasing, or absence of tumor contamination."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {"output": str(args.output), "totals": receipt["totals"], "pass": receipt["pass"]},
            ensure_ascii=False,
            indent=2,
        )
    )
    raise SystemExit(0 if receipt["pass"] else 1)


if __name__ == "__main__":
    main()
