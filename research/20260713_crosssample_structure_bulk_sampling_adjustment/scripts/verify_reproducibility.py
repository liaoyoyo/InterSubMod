#!/usr/bin/env python3
"""Compare two complete analysis directories by relative-file SHA-256."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def files(root: Path) -> dict[str, Path]:
    return {
        str(path.relative_to(root)): path
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not args.first.is_dir() or not args.second.is_dir():
        raise FileNotFoundError("Both comparison inputs must be directories")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    first, second = files(args.first), files(args.second)
    names = sorted(set(first) | set(second))
    rows = []
    for name in names:
        first_hash = sha256(first[name]) if name in first else ""
        second_hash = sha256(second[name]) if name in second else ""
        rows.append({
            "relative_file": name,
            "present_first": name in first,
            "present_second": name in second,
            "sha256_first": first_hash,
            "sha256_second": second_hash,
            "byte_identical": bool(first_hash and first_hash == second_hash),
        })
    status = "PASS" if rows and all(row["byte_identical"] for row in rows) else "FAIL"

    tsv_path = args.output_dir / "reproducibility_hashes.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "relative_file", "present_first", "present_second",
                "sha256_first", "sha256_second", "byte_identical",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    receipt = {
        "schema_name": "intersubmod.crosssample_structure_bulk_sampling_adjustment.reproducibility",
        "schema_version": "1.0.0",
        "analysis_date": "2026-07-13",
        "status": status,
        "first": str(args.first.resolve()),
        "second": str(args.second.resolve()),
        "files_compared": len(rows),
        "byte_identical_files": sum(row["byte_identical"] for row in rows),
        "mismatches": [row["relative_file"] for row in rows if not row["byte_identical"]],
        "comparison_table": str(tsv_path.resolve()),
    }
    (args.output_dir / "reproducibility_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    report = f"""<!--
建立時間: 2026-07-13
目標: 驗證固定 seed 的兩次完整 compositional analysis 是否 byte-identical
處理範圍: 全部分析輸出檔
關聯檔案: reproducibility_hashes.tsv
-->

# Reproducibility receipt

> **{status}：{receipt['byte_identical_files']}/{receipt['files_compared']} files byte-identical。**

- First run：`{args.first.resolve()}`
- Second run：`{args.second.resolve()}`
- 比對方式：相對檔名集合＋每檔 SHA-256；包含 stochastic summaries、JSON、report與 provenance。
- Fixed seed與replicate數由各 run的 `provenance.json` 記錄。
"""
    (args.output_dir / "reproducibility_receipt.md").write_text(report, encoding="utf-8")
    print(f"STATUS -> {status} ({receipt['byte_identical_files']}/{receipt['files_compared']} byte-identical files)")
    print(f"OUTPUT -> {args.output_dir.resolve()}")
    if status != "PASS":
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
