#!/usr/bin/env python3
"""V3F vs V5 vs V6 BAM three-way head-to-head comparison.

Extends v6c_phaseB_v3f_vs_v5_compare.py to include V6 BAM (Layer 1.5 reverted).
Outputs full evaluation matrix per Phase B 驗收 criteria.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
BASES = {
    "V3F": REPO / "research/paired_priority_bug_audit/v6c_phaseB_runs",
    "V5": REPO / "research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs",
    "V6": REPO / "research/paired_priority_bug_audit/v6c_phaseB_v6bam_runs",
}
RUNS = ["off_tp", "off_fp", "on_tp", "on_fp"]
OUT_BASE = REPO / "research/paired_priority_bug_audit/v6c_phaseB_v6bam_runs"


def collect(run_dir: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for reads_tsv in run_dir.rglob("reads.tsv"):
        region_dir = reads_tsv.parent.parent
        region_id = region_dir.name
        hp_counter: Counter[str] = Counter()
        with reads_tsv.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                hp_idx = header.index("hp")
            except ValueError:
                continue
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= hp_idx:
                    continue
                hp_counter[cols[hp_idx]] += 1
        out[region_id] = {"hp_counter": hp_counter, "n_reads": sum(hp_counter.values())}
    return out


def compute_ng(hp_counter: Counter[str]) -> int:
    return sum(1 for k, v in hp_counter.items() if k not in {"0", "hp"} and v > 0)


def run_summary(data: dict[str, dict[str, dict]], label: str, off_run: str, on_run: str):
    out = []
    for region_id in sorted(set(data[off_run]) | set(data[on_run])):
        off = data[off_run].get(region_id)
        on = data[on_run].get(region_id)
        if not off or not on:
            continue
        ng_off = compute_ng(off["hp_counter"])
        ng_on = compute_ng(on["hp_counter"])
        n_reads = off["n_reads"]
        out.append({
            "region_id": region_id,
            "label": label,
            "ng_off": ng_off,
            "ng_on": ng_on,
            "n_reads": n_reads,
        })
    return out


def aggregate_marker(rows):
    tp_off = sum(1 for r in rows if r["label"] == "TP" and r["ng_off"] >= 3)
    fp_off = sum(1 for r in rows if r["label"] == "FP" and r["ng_off"] >= 3)
    tp_on2 = sum(1 for r in rows if r["label"] == "TP" and r["ng_on"] == 2)
    fp_on2 = sum(1 for r in rows if r["label"] == "FP" and r["ng_on"] == 2)
    return {
        "marker_off_n": tp_off + fp_off,
        "marker_off_tp": tp_off,
        "marker_off_fp": fp_off,
        "marker_off_rate": tp_off / max(tp_off + fp_off, 1),
        "ng_on2_tp": tp_on2,
        "ng_on2_fp": fp_on2,
        "ng_on2_rate": tp_on2 / max(tp_on2 + fp_on2, 1),
    }


def main() -> int:
    print("[V3F vs V5 vs V6 three-way] Loading BAMs...")
    data: dict[str, dict[str, dict[str, dict]]] = {}
    for tag, base in BASES.items():
        if not base.exists():
            print(f"  [SKIP] {tag} ({base}) not found")
            continue
        data[tag] = {run: collect(base / run) for run in RUNS}
        for run in RUNS:
            print(f"  {tag} {run}: {len(data[tag][run])} regions")

    available = list(data.keys())
    if "V6" not in available:
        print("\n[WARN] V6 BAM not yet available (haplotag may still be running)")
    print(f"\n=== Available BAMs: {available} ===\n")

    # Read-level hp distribution per BAM
    print("=== chr19 全部 reads hp distribution per BAM ===")
    print(f"{'hp_value':<10}", end="")
    for tag in available:
        print(f" {tag+' off':>11} {tag+' on':>11}", end="")
    print()

    all_hp = set()
    for tag in available:
        for run in RUNS:
            for v in data[tag][run].values():
                all_hp.update(v["hp_counter"].keys())

    for hp in sorted(all_hp):
        print(f"{hp:<10}", end="")
        for tag in available:
            off_cnt = sum(d["hp_counter"].get(hp, 0) for d in data[tag]["off_tp"].values()) + sum(
                d["hp_counter"].get(hp, 0) for d in data[tag]["off_fp"].values()
            )
            on_cnt = sum(d["hp_counter"].get(hp, 0) for d in data[tag]["on_tp"].values()) + sum(
                d["hp_counter"].get(hp, 0) for d in data[tag]["on_fp"].values()
            )
            print(f" {off_cnt:>11} {on_cnt:>11}", end="")
        print()

    # hp=1-1 vs hp=2-1 ratio per BAM
    print("\n=== hp=1-1 vs hp=2-1 ratio (priority bug feature 化指標) ===")
    print(f"{'BAM':<6} {'hp=1-1':>10} {'hp=2-1':>10} {'ratio':>10}")
    for tag in available:
        h11 = sum(d["hp_counter"].get("1-1", 0) for d in data[tag]["off_tp"].values()) + sum(
            d["hp_counter"].get("1-1", 0) for d in data[tag]["off_fp"].values()
        )
        h21 = sum(d["hp_counter"].get("2-1", 0) for d in data[tag]["off_tp"].values()) + sum(
            d["hp_counter"].get("2-1", 0) for d in data[tag]["off_fp"].values()
        )
        ratio = h11 / max(h21, 1)
        print(f"{tag:<6} {h11:>10} {h21:>10} {ratio:>10.3f}")

    # hp=3 ambiguous reads per BAM
    print("\n=== hp=3 (somatic ambiguous) reads per BAM ===")
    print(f"{'BAM':<6} {'hp=3':>10}")
    for tag in available:
        h3 = sum(d["hp_counter"].get("3", 0) for d in data[tag]["off_tp"].values()) + sum(
            d["hp_counter"].get("3", 0) for d in data[tag]["off_fp"].values()
        )
        print(f"{tag:<6} {h3:>10}")

    # Marker filter aggregate per BAM
    summary_rows: list[dict] = []
    print("\n=== Marker filter (NG_off ≥ 3) — three-way ===")
    print(f"{'BAM':<6} {'NG≥3 N':>8} {'TP':>5} {'FP':>5} {'rate':>7} | {'NG_on=2 TP':>11} {'FP':>5} {'rate':>7}")
    for tag in available:
        rows = run_summary(data[tag], "TP", "off_tp", "on_tp") + run_summary(
            data[tag], "FP", "off_fp", "on_fp"
        )
        agg = aggregate_marker(rows)
        print(
            f"{tag:<6} {agg['marker_off_n']:>8} {agg['marker_off_tp']:>5} "
            f"{agg['marker_off_fp']:>5} {agg['marker_off_rate']:>7.3f} | "
            f"{agg['ng_on2_tp']:>11} {agg['ng_on2_fp']:>5} {agg['ng_on2_rate']:>7.3f}"
        )
        summary_rows.append({"BAM": tag, **agg})

    # Per-cell TP rate cross-tab per BAM
    for tag in available:
        rows = run_summary(data[tag], "TP", "off_tp", "on_tp") + run_summary(
            data[tag], "FP", "off_fp", "on_fp"
        )
        cell_tp: dict[tuple[int, int], int] = defaultdict(int)
        cell_fp: dict[tuple[int, int], int] = defaultdict(int)
        for r in rows:
            key = (r["ng_off"], r["ng_on"])
            if r["label"] == "TP":
                cell_tp[key] += 1
            else:
                cell_fp[key] += 1
        print(f"\n=== {tag} BAM per-cell TP rate ===")
        print(f"{'NG_off':>7} {'NG_on':>6} {'TP':>5} {'FP':>5} {'rate':>7}")
        for noff, non in sorted(set(cell_tp) | set(cell_fp)):
            tp = cell_tp[(noff, non)]
            fp = cell_fp[(noff, non)]
            rate = tp / max(tp + fp, 1)
            print(f"{noff:>7} {non:>6} {tp:>5} {fp:>5} {rate:>7.3f}")

    # Per-region NG agreement (V3F vs V6, V5 vs V6)
    if "V6" in available:
        print("\n=== Per-region NG_off agreement (V6 vs others) ===")
        v6_ng = {
            (r["region_id"], r["label"]): r["ng_off"]
            for r in (
                run_summary(data["V6"], "TP", "off_tp", "on_tp")
                + run_summary(data["V6"], "FP", "off_fp", "on_fp")
            )
        }
        for tag in ["V3F", "V5"]:
            if tag not in data:
                continue
            tag_ng = {
                (r["region_id"], r["label"]): r["ng_off"]
                for r in (
                    run_summary(data[tag], "TP", "off_tp", "on_tp")
                    + run_summary(data[tag], "FP", "off_fp", "on_fp")
                )
            }
            common = set(v6_ng) & set(tag_ng)
            n_diff = sum(1 for k in common if v6_ng[k] != tag_ng[k])
            n_total = len(common)
            print(f"  V6 vs {tag}: {n_diff}/{n_total} ({100 * n_diff / max(n_total, 1):.1f}%) regions disagree")

    # Phase B 驗收標準（V6 only）
    if "V6" in available:
        print("\n=== Phase B 驗收標準（V6 是否成功消除 priority bug feature 化）===")
        v3f_h11 = sum(d["hp_counter"].get("1-1", 0) for d in data["V3F"]["off_tp"].values()) + sum(
            d["hp_counter"].get("1-1", 0) for d in data["V3F"]["off_fp"].values()
        )
        v3f_h21 = sum(d["hp_counter"].get("2-1", 0) for d in data["V3F"]["off_tp"].values()) + sum(
            d["hp_counter"].get("2-1", 0) for d in data["V3F"]["off_fp"].values()
        )
        v6_h11 = sum(d["hp_counter"].get("1-1", 0) for d in data["V6"]["off_tp"].values()) + sum(
            d["hp_counter"].get("1-1", 0) for d in data["V6"]["off_fp"].values()
        )
        v6_h21 = sum(d["hp_counter"].get("2-1", 0) for d in data["V6"]["off_tp"].values()) + sum(
            d["hp_counter"].get("2-1", 0) for d in data["V6"]["off_fp"].values()
        )
        v6_ratio = v6_h11 / max(v6_h21, 1)
        v3f_ratio = v3f_h11 / max(v3f_h21, 1)
        print(f"  hp=1-1:hp=2-1 ratio  V3F={v3f_ratio:.3f}  V6={v6_ratio:.3f}  → {'PASS ratio<1.0' if v6_ratio < 1.0 else 'FAIL ratio≥1.0'}")

        v3f_h3 = sum(d["hp_counter"].get("3", 0) for d in data["V3F"]["off_tp"].values()) + sum(
            d["hp_counter"].get("3", 0) for d in data["V3F"]["off_fp"].values()
        )
        v6_h3 = sum(d["hp_counter"].get("3", 0) for d in data["V6"]["off_tp"].values()) + sum(
            d["hp_counter"].get("3", 0) for d in data["V6"]["off_fp"].values()
        )
        print(f"  hp=3 (ambiguous)     V3F={v3f_h3}  V6={v6_h3}  → {'PASS V6≥1500' if v6_h3 >= 1500 else 'FAIL V6<1500'}")

        v3f_summary = next(s for s in summary_rows if s["BAM"] == "V3F")
        v6_summary = next(s for s in summary_rows if s["BAM"] == "V6")
        print(
            f"  marker coverage (NG≥3) V3F={v3f_summary['marker_off_n']}  V6={v6_summary['marker_off_n']}  "
            f"→ {'PASS V6≥450' if v6_summary['marker_off_n'] >= 450 else 'FAIL V6<450'}"
        )
        print(
            f"  marker rate (off)    V3F={v3f_summary['marker_off_rate']:.3f}  V6={v6_summary['marker_off_rate']:.3f}  "
            f"→ {'PASS V6≥0.94' if v6_summary['marker_off_rate'] >= 0.94 else 'FAIL V6<0.94'}"
        )
        print(
            f"  flag=on NG_on=2 rate  V3F={v3f_summary['ng_on2_rate']:.3f}  V6={v6_summary['ng_on2_rate']:.3f}  "
            f"→ {'PASS V6≥0.90' if v6_summary['ng_on2_rate'] >= 0.90 else 'FAIL V6<0.90'}"
        )

    # Save TSV
    out_tsv = OUT_BASE / "v3f_vs_v5_vs_v6_summary.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        if summary_rows:
            writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
    print(f"\n[three-way] summary → {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
