#!/usr/bin/env python3
"""Summarize region_samtools_snapshot outputs into one TSV/Markdown report."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

from research_common import load_tsv_rows, write_tsv_rows


FIELDS = [
    "sample",
    "region_id",
    "snapshot_dir",
    "bam",
    "region",
    "target_chrom",
    "target_pos",
    "target_ref",
    "target_alt",
    "read_count",
    "target_depth",
    "target_ref_count",
    "target_alt_count",
    "target_alt_fraction",
    "target_other_count",
    "target_mpileup_depth",
    "top_hp_tag",
    "top_hp_count",
    "top_hp_fraction",
    "hp1_collapsed_count",
    "hp2_collapsed_count",
    "other_hp_count",
    "collapsed_hp_balance_delta",
    "na_hp_count",
    "na_hp_fraction",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot-root", action="append", required=True, help="Root directory containing snapshot dirs")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def read_summary(path: Path) -> Dict[str, str]:
    payload: Dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        payload[key.strip()] = value.strip()
    return payload


def parse_region_id(region_id: str) -> Tuple[str, str, str, str]:
    parts = region_id.split("_")
    if len(parts) < 4:
        return "", "", "", ""
    return parts[0], parts[1], parts[2], parts[3]


def parse_depth_at_target(path: Path, target_pos: str) -> str:
    if not path.exists() or not target_pos:
        return ""
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.split("\t")
        if len(parts) >= 3 and parts[1] == target_pos:
            return parts[2]
    return ""


def clean_mpileup_bases(bases: str) -> str:
    cleaned: List[str] = []
    i = 0
    while i < len(bases):
        base = bases[i]
        if base == "^":
            i += 2
            continue
        if base == "$":
            i += 1
            continue
        if base in "+-":
            i += 1
            digits: List[str] = []
            while i < len(bases) and bases[i].isdigit():
                digits.append(bases[i])
                i += 1
            indel_len = int("".join(digits)) if digits else 0
            i += indel_len
            continue
        cleaned.append(base)
        i += 1
    return "".join(cleaned)


def parse_mpileup_target(path: Path, target_pos: str, target_ref: str, target_alt: str) -> Tuple[str, str, str, str, str]:
    if not path.exists() or not target_pos:
        return "", "", "", "", ""
    ref_base = target_ref.upper()
    alt_base = target_alt.upper()
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.split("\t")
        if len(parts) < 5 or parts[1] != target_pos:
            continue
        depth = parts[3]
        bases = clean_mpileup_bases(parts[4])
        ref_count = 0
        alt_count = 0
        other_count = 0
        for base in bases:
            if base in ".,":  # pileup reference match
                ref_count += 1
            elif base.upper() == alt_base:
                alt_count += 1
            elif base.upper() == ref_base:
                ref_count += 1
            elif base == "*":
                other_count += 1
            else:
                other_count += 1
        denom = ref_count + alt_count + other_count
        alt_fraction = f"{(alt_count / denom):.6f}" if denom else ""
        return depth, str(ref_count), str(alt_count), alt_fraction, str(other_count)
    return "", "", "", "", ""


def hp_stats(path: Path) -> Tuple[str, str, str, str, str, str, str, str, str]:
    rows = load_tsv_rows(path)
    if not rows:
        return "", "", "", "", "", "", "", "", ""
    total = sum(int(row.get("count", 0)) for row in rows)
    top_tag = ""
    top_count = -1
    na_count = 0
    hp1_count = 0
    hp2_count = 0
    other_count = 0
    for row in rows:
        tag = str(row.get("hp_tag", ""))
        count = int(row.get("count", 0))
        if count > top_count:
            top_tag = tag
            top_count = count
        if tag == "NA":
            na_count = count
        elif tag.startswith("1"):
            hp1_count += count
        elif tag.startswith("2"):
            hp2_count += count
        else:
            other_count += count
    top_fraction = f"{(top_count / total):.6f}" if total else ""
    na_fraction = f"{(na_count / total):.6f}" if total else ""
    hp_balance_den = hp1_count + hp2_count
    hp_balance_delta = f"{(abs(hp1_count - hp2_count) / hp_balance_den):.6f}" if hp_balance_den else ""
    return (
        top_tag,
        str(top_count if top_count >= 0 else ""),
        top_fraction,
        str(hp1_count),
        str(hp2_count),
        str(other_count),
        hp_balance_delta,
        str(na_count),
        na_fraction,
    )


def sample_from_root(root: Path) -> str:
    return root.name


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, str]] = []
    for root_str in args.snapshot_root:
        root = Path(root_str).resolve()
        sample = sample_from_root(root)
        if not root.exists():
            continue
        snapshot_dirs = sorted(summary_path.parent for summary_path in root.rglob("summary.txt"))
        for snapshot_dir in snapshot_dirs:
            summary_path = snapshot_dir / "summary.txt"
            summary = read_summary(summary_path)
            region_id = snapshot_dir.parent.name if snapshot_dir.name == "samtools_snapshot" else snapshot_dir.name
            target_chrom, target_pos, target_ref, target_alt = parse_region_id(region_id)
            (
                top_tag,
                top_count,
                top_fraction,
                hp1_count,
                hp2_count,
                other_count,
                hp_balance_delta,
                na_count,
                na_fraction,
            ) = hp_stats(snapshot_dir / "hp_tag_counts.tsv")
            target_depth = parse_depth_at_target(snapshot_dir / "depth.tsv", target_pos)
            target_mpileup_depth, target_ref_count, target_alt_count, target_alt_fraction, target_other_count = (
                parse_mpileup_target(snapshot_dir / "mpileup.txt", target_pos, target_ref, target_alt)
            )
            read_count = summary.get("read_count", "")
            notes: List[str] = []
            if na_fraction and float(na_fraction) >= 0.15:
                notes.append("high_na_hp")
            if top_fraction and float(top_fraction) >= 0.60:
                notes.append("single_hp_dominant")
            if hp_balance_delta and float(hp_balance_delta) >= 0.60:
                notes.append("haplotype_skewed")
            if target_alt_fraction and float(target_alt_fraction) <= 0.15:
                notes.append("low_alt_fraction")
            rows.append(
                {
                    "sample": sample,
                    "region_id": region_id,
                    "snapshot_dir": str(snapshot_dir),
                    "bam": summary.get("bam", ""),
                    "region": summary.get("region", ""),
                    "read_count": read_count,
                    "target_chrom": target_chrom,
                    "target_pos": target_pos,
                    "target_ref": target_ref,
                    "target_alt": target_alt,
                    "target_depth": target_depth,
                    "target_ref_count": target_ref_count,
                    "target_alt_count": target_alt_count,
                    "target_alt_fraction": target_alt_fraction,
                    "target_other_count": target_other_count,
                    "target_mpileup_depth": target_mpileup_depth,
                    "top_hp_tag": top_tag,
                    "top_hp_count": top_count,
                    "top_hp_fraction": top_fraction,
                    "hp1_collapsed_count": hp1_count,
                    "hp2_collapsed_count": hp2_count,
                    "other_hp_count": other_count,
                    "collapsed_hp_balance_delta": hp_balance_delta,
                    "na_hp_count": na_count,
                    "na_hp_fraction": na_fraction,
                    "notes": ",".join(notes),
                }
            )

    write_tsv_rows(output_dir / "snapshot_summary.tsv", FIELDS, rows)

    lines = [
        "# Snapshot Summary",
        "",
        "- 主表：`snapshot_summary.tsv`",
        "",
        "| sample | region_id | read_count | target_depth | target_alt_fraction | top_hp_tag | top_hp_fraction | collapsed_hp_balance_delta | na_hp_fraction | notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['sample']} | {row['region_id']} | {row['read_count']} | {row['target_depth']} | {row['target_alt_fraction']} | {row['top_hp_tag']} | {row['top_hp_fraction']} | {row['collapsed_hp_balance_delta']} | {row['na_hp_fraction']} | {row['notes']} |"
        )
    (output_dir / "snapshot_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[summarize_samtools_snapshots] Wrote {output_dir / 'snapshot_summary.tsv'}")


if __name__ == "__main__":
    main()
