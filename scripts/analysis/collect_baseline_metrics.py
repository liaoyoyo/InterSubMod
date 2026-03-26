#!/usr/bin/env python3
"""
collect_baseline_metrics.py

讀取所有樣本的現有 to_stage_metrics.json，輸出跨樣本 F1 基線對照表。
不重跑任何分析，僅彙整已存在的輸出結果。

用法:
  python3 scripts/analysis/collect_baseline_metrics.py
  python3 scripts/analysis/collect_baseline_metrics.py --output baseline_summary.tsv
"""

import json
import os
import sys
import argparse

SYNTHESIS_BASE = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds"

# 每個樣本的「最新最完整 TO round」，優先使用 fastresume（完整跑完）
SAMPLE_ROUNDS = [
    ("HCC1395-5kHz",    "20260315_hcc1395_to_pilot"),
    ("HCC1395-DORADO",  "20260315_hcc1395_dorado_to_pilot"),
    ("COLO829",         "20260317_colo829_to_pilot_1"),
    ("H1437",           "20260318_h1437_to_pilot_fastresume"),
    ("H2009",           "20260318_h2009_to_pilot_fastresume"),
    ("HCC1937",         "20260318_hcc1937_to_pilot_fastresume"),
    ("HCC1954",         "20260318_hcc1954_to_pilot_fastresume"),
]


def load_metrics(round_dir: str) -> dict:
    """讀取單個 round 的 to_stage_metrics.json。"""
    path = os.path.join(SYNTHESIS_BASE, round_dir, "to_stage_metrics.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def collect_all() -> list[dict]:
    """彙整所有樣本的基線 metrics。"""
    rows = []
    for sample, round_dir in SAMPLE_ROUNDS:
        m = load_metrics(round_dir)
        baseline = m.get("baseline_clairs_to", {})
        ism = m.get("intersubmod", {})

        baseline_f1 = baseline.get("f1")
        ism_f1 = ism.get("f1")
        delta = (ism_f1 - baseline_f1) if (ism_f1 is not None and baseline_f1 is not None) else None

        rows.append({
            "sample": sample,
            "round": round_dir,
            "baseline_precision": baseline.get("precision"),
            "baseline_recall": baseline.get("recall"),
            "baseline_f1": baseline_f1,
            "ism_precision": ism.get("precision"),
            "ism_recall": ism.get("recall"),
            "ism_f1": ism_f1,
            "ism_delta": delta,
            "baseline_tp": baseline.get("tp"),
            "baseline_fp": baseline.get("fp"),
            "baseline_fn": baseline.get("fn"),
        })
    return rows


def fmt(val, spec=".4f"):
    if val is None:
        return "N/A"
    return format(val, spec)


def print_summary_table(rows: list[dict]) -> None:
    """列印對照表至 stdout。"""
    header = "{:<22} {:<10} {:<10} {:<10} {}".format(
        "樣本", "基線 F1", "ISM F1", "delta", "Round"
    )
    print(header)
    print("-" * 85)
    for r in rows:
        b = fmt(r["baseline_f1"])
        i = fmt(r["ism_f1"])
        d = fmt(r["ism_delta"], "+.6f") if r["ism_delta"] is not None else "N/A"
        line = "{:<22} {:<10} {:<10} {:<10} {}".format(
            r["sample"], b, i, d, r["round"]
        )
        print(line)

    # 統計摘要
    valid = [r for r in rows if r["ism_delta"] is not None]
    if valid:
        avg_delta = sum(r["ism_delta"] for r in valid) / len(valid)
        pos = sum(1 for r in valid if r["ism_delta"] > 0)
        print()
        print("摘要: {} 個樣本，平均 ISM delta = {:+.6f}，正增益 {}/{} 個".format(
            len(valid), avg_delta, pos, len(valid)
        ))


def write_tsv(rows: list[dict], output_path: str) -> None:
    """輸出詳細 TSV 供後續分析。"""
    columns = [
        "sample", "round",
        "baseline_f1", "baseline_precision", "baseline_recall",
        "ism_f1", "ism_precision", "ism_recall", "ism_delta",
        "baseline_tp", "baseline_fp", "baseline_fn",
    ]
    with open(output_path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            vals = []
            for c in columns:
                v = r.get(c)
                vals.append("" if v is None else str(v))
            f.write("\t".join(vals) + "\n")
    print(f"\n已輸出 TSV: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="收集跨樣本 TO ISM 基線 F1 對照表（讀取現有 metrics，不重跑）"
    )
    parser.add_argument("--output", "-o", default=None,
                        help="輸出 TSV 路徑（可選，預設僅印至 stdout）")
    args = parser.parse_args()

    rows = collect_all()
    print_summary_table(rows)

    if args.output:
        write_tsv(rows, args.output)


if __name__ == "__main__":
    main()
