#!/usr/bin/env python3
"""Run same-scope snapshot control by comparing full vs subset BAM snapshots."""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Tuple


REPO_ROOT = Path("/big8_disk/liaoyoyo2001/InterSubMod")
SAMTOOLS_SNAPSHOT = REPO_ROOT / "scripts" / "analysis" / "region_samtools_snapshot.sh"
SUMMARIZE_SNAPSHOTS = REPO_ROOT / "scripts" / "analysis" / "summarize_samtools_snapshots.py"

DEFAULT_DIAGNOSTIC_SUMMARY = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
    "20260311_to_support_feature_diagnostics/diagnostic_summary.tsv"
)
DEFAULT_FULL_SNAPSHOT_SUMMARY = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
    "20260311_to_support_feature_diagnostics/hcc1395_5khz_to/snapshot_summary/snapshot_summary.tsv"
)
DEFAULT_FULL_BAM = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
    "20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam"
)

NUMERIC_FIELDS = [
    "read_count",
    "target_depth",
    "target_ref_count",
    "target_alt_count",
    "target_alt_fraction",
    "target_other_count",
    "target_mpileup_depth",
    "top_hp_fraction",
    "hp1_collapsed_count",
    "hp2_collapsed_count",
    "other_hp_count",
    "collapsed_hp_balance_delta",
    "na_hp_count",
    "na_hp_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-id", default="hcc1395_5khz_to", help="Dataset id inside diagnostic_summary.tsv")
    parser.add_argument("--diagnostic-summary", default=str(DEFAULT_DIAGNOSTIC_SUMMARY))
    parser.add_argument("--full-snapshot-summary", default=str(DEFAULT_FULL_SNAPSHOT_SUMMARY))
    parser.add_argument("--full-bam", default=str(DEFAULT_FULL_BAM))
    parser.add_argument("--max-reads", type=int, default=120)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def run(cmd: List[str]) -> None:
    subprocess.run(cmd, check=True)


def load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def parse_region(region: str) -> Tuple[str, int, int]:
    chrom, rest = region.split(":")
    start_s, end_s = rest.split("-")
    return chrom, int(start_s), int(end_s)


def normalize_region_key(key: str) -> str:
    return key.replace(":", "_")


def as_float(value: str) -> float | None:
    if value in ("", "NA", None):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def make_bed_rows(rows: List[Dict[str, str]]) -> List[Tuple[str, int, int, str]]:
    dedup: Dict[Tuple[str, int, int], str] = {}
    for row in rows:
        chrom, start, end = parse_region(row["region"])
        bed_key = (chrom, start - 1, end)
        dedup[bed_key] = row["region_key"]
    return [(chrom, start0, end, key) for (chrom, start0, end), key in sorted(dedup.items())]


def build_subset_bam(full_bam: Path, bed_path: Path, subset_bam: Path) -> None:
    run(["samtools", "view", "-b", "-M", "-L", str(bed_path), str(full_bam), "-o", str(subset_bam)])
    run(["samtools", "index", "-b", str(subset_bam)])


def subset_snapshot_root(output_dir: Path, dataset_id: str) -> Path:
    return output_dir / dataset_id / "subset_snapshots"


def run_subset_snapshots(rows: List[Dict[str, str]], subset_bam: Path, output_dir: Path, max_reads: int) -> Path:
    root = subset_snapshot_root(output_dir, rows[0]["dataset_id"])
    for row in rows:
        region_id = normalize_region_key(row["region_key"])
        snap_dir = root / "diagnostics" / row["downstream_status"] / region_id / "samtools_snapshot"
        snap_dir.mkdir(parents=True, exist_ok=True)
        run(
            [
                "bash",
                str(SAMTOOLS_SNAPSHOT),
                "--bam",
                str(subset_bam),
                "--region",
                row["region"],
                "--output-dir",
                str(snap_dir),
                "--max-reads",
                str(max_reads),
            ]
        )
    summary_dir = root / "snapshot_summary"
    run(
        [
            "python3",
            str(SUMMARIZE_SNAPSHOTS),
            "--snapshot-root",
            str(root),
            "--output-dir",
            str(summary_dir),
        ]
    )
    return summary_dir / "snapshot_summary.tsv"


def compare_rows(
    rows: List[Dict[str, str]], full_summary: Dict[str, Dict[str, str]], subset_summary: Dict[str, Dict[str, str]]
) -> List[Dict[str, object]]:
    compared = []
    for row in rows:
        region_id = normalize_region_key(row["region_key"])
        full_row = full_summary.get(region_id, {})
        subset_row = subset_summary.get(region_id, {})
        out: Dict[str, object] = {
            "dataset_id": row["dataset_id"],
            "region_key": row["region_key"],
            "downstream_status": row["downstream_status"],
            "truth_status": row["truth_status"],
            "region": row["region"],
        }
        identical_all = True
        for field in NUMERIC_FIELDS:
            full_v = as_float(full_row.get(field, ""))
            subset_v = as_float(subset_row.get(field, ""))
            out[f"{field}_full"] = "" if full_v is None else full_v
            out[f"{field}_subset"] = "" if subset_v is None else subset_v
            if full_v is None or subset_v is None:
                out[f"{field}_abs_delta"] = ""
                continue
            delta = subset_v - full_v
            out[f"{field}_abs_delta"] = delta
            if abs(delta) > 1e-9:
                identical_all = False
        out["top_hp_tag_full"] = full_row.get("top_hp_tag", "")
        out["top_hp_tag_subset"] = subset_row.get("top_hp_tag", "")
        out["top_hp_tag_changed"] = str(full_row.get("top_hp_tag", "") != subset_row.get("top_hp_tag", "")).lower()
        if out["top_hp_tag_changed"] == "true":
            identical_all = False
        out["identical_all_metrics"] = str(identical_all).lower()
        compared.append(out)
    return compared


def summarize_bias(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    result = []
    identical_count = sum(1 for row in rows if row["identical_all_metrics"] == "true")
    for field in NUMERIC_FIELDS:
        deltas = [abs(float(row[f"{field}_abs_delta"])) for row in rows if row[f"{field}_abs_delta"] != ""]
        result.append(
            {
                "metric": field,
                "rows_compared": len(deltas),
                "max_abs_delta": max(deltas) if deltas else "",
                "mean_abs_delta": mean(deltas) if deltas else "",
                "all_identical": str(all(delta <= 1e-9 for delta in deltas)).lower() if deltas else "",
                "identical_rows_total": identical_count,
                "row_count": len(rows),
            }
        )
    result.append(
        {
            "metric": "top_hp_tag",
            "rows_compared": len(rows),
            "max_abs_delta": "",
            "mean_abs_delta": "",
            "all_identical": str(all(row["top_hp_tag_changed"] == "false" for row in rows)).lower(),
            "identical_rows_total": identical_count,
            "row_count": len(rows),
        }
    )
    return result


def write_summary_md(
    path: Path,
    dataset_id: str,
    subset_bam: Path,
    rows: List[Dict[str, object]],
    metric_rows: List[Dict[str, object]],
) -> None:
    identical_count = sum(1 for row in rows if row["identical_all_metrics"] == "true")
    lines = [
        "# Same-scope Snapshot Bias Summary",
        "",
        f"- dataset: `{dataset_id}`",
        f"- subset_bam: `{subset_bam}`",
        f"- compared_regions: `{len(rows)}`",
        f"- fully_identical_regions: `{identical_count}`",
        "",
        "## 結論",
        "",
    ]
    if identical_count == len(rows):
        lines.extend(
            [
                "1. `candidate-window subset BAM` 在同一組 `5kHz TO` region 上，對 snapshot 指標沒有觀察到實質偏移。",
                "2. 在本輪控制下，`subset snapshot` 對同 dataset 的 read-level ranking 與方向判讀可視為安全近似。",
                "3. 這代表先前 `5kHz TO full` vs `DORADO TO subset` 的主要限制，較可能來自 dataset/platform 差異，而不是 subset 技術本身。",
            ]
        )
    else:
        lines.extend(
            [
                "1. `candidate-window subset BAM` 對部分 snapshot 指標有可觀察偏移，不能直接當作 full BAM 的等價近似。",
                "2. 之後若要做跨 dataset read-level 比較，需先補 full-BAM control 或額外偏移校正。",
            ]
        )
    lines.extend(
        [
            "",
            "## 指標摘要",
            "",
            "| metric | rows_compared | max_abs_delta | mean_abs_delta | all_identical |",
            "| --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in metric_rows:
        lines.append(
            f"| `{row['metric']}` | {row['rows_compared']} | {row['max_abs_delta']} | {row['mean_abs_delta']} | {row['all_identical']} |"
        )
    lines.extend(
        [
            "",
            "## 判讀",
            "",
            "- 若 `read_count / target_depth / target_alt_fraction / na_hp_fraction / collapsed_hp_balance_delta` 都相同，代表 subset 只是在磁碟層縮小 BAM，不會改變相同 region 的 read-level snapshot。",
            "- 這個控制只能回答 `subset` 對同 dataset 的影響，不能直接消除 `5kHz` 與 `DORADO` 在平台、候選集合與 aggregate feature 方向上的差異。",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    diagnostic_summary = Path(args.diagnostic_summary).resolve()
    full_snapshot_summary = Path(args.full_snapshot_summary).resolve()
    full_bam = Path(args.full_bam).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = [row for row in load_tsv(diagnostic_summary) if row["dataset_id"] == args.dataset_id]
    if not rows:
        raise SystemExit(f"No rows found for dataset_id={args.dataset_id}")

    bed_rows = make_bed_rows(rows)
    bed_path = output_dir / f"{args.dataset_id}_selected_regions.bed"
    with bed_path.open("w", encoding="utf-8") as handle:
        for chrom, start0, end, key in bed_rows:
            handle.write(f"{chrom}\t{start0}\t{end}\t{key}\n")

    subset_bam = output_dir / f"{args.dataset_id}_candidate_windows_subset.bam"
    build_subset_bam(full_bam, bed_path, subset_bam)
    subset_summary_path = run_subset_snapshots(rows, subset_bam, output_dir, args.max_reads)

    full_summary = {row["region_id"]: row for row in load_tsv(full_snapshot_summary)}
    subset_summary = {row["region_id"]: row for row in load_tsv(subset_summary_path)}
    compared_rows = compare_rows(rows, full_summary, subset_summary)
    metric_rows = summarize_bias(compared_rows)

    compare_fields = [
        "dataset_id",
        "region_key",
        "downstream_status",
        "truth_status",
        "region",
        "top_hp_tag_full",
        "top_hp_tag_subset",
        "top_hp_tag_changed",
        "identical_all_metrics",
    ]
    for field in NUMERIC_FIELDS:
        compare_fields.extend([f"{field}_full", f"{field}_subset", f"{field}_abs_delta"])
    write_tsv(output_dir / "same_scope_snapshot_comparison.tsv", compare_fields, compared_rows)
    write_tsv(
        output_dir / "same_scope_metric_summary.tsv",
        ["metric", "rows_compared", "max_abs_delta", "mean_abs_delta", "all_identical", "identical_rows_total", "row_count"],
        metric_rows,
    )
    write_summary_md(output_dir / "same_scope_snapshot_bias_summary.md", args.dataset_id, subset_bam, compared_rows, metric_rows)

    print(f"[run_snapshot_scope_same_control] Wrote {output_dir / 'same_scope_snapshot_comparison.tsv'}")
    print(f"[run_snapshot_scope_same_control] Wrote {output_dir / 'same_scope_metric_summary.tsv'}")
    print(f"[run_snapshot_scope_same_control] Wrote {output_dir / 'same_scope_snapshot_bias_summary.md'}")


if __name__ == "__main__":
    main()
