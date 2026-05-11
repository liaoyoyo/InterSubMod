#!/usr/bin/env python3
"""V6-C Phase B head-to-head comparison: V3F BAM vs V5 BAM × 2 flag (germline-hp-only on/off).

Compares per-region NG distribution + marker filter TP rate between V3F and V5 BAMs
to evaluate whether V5 Layer 1.5 design defect (germline-absent 4.19:1 偏 HP1) impacts
region-level marker behavior beyond the chr19 V3F-only test (V6C Phase B 5/10).
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
V3F_BASE = REPO / "research/paired_priority_bug_audit/v6c_phaseB_runs"
V5_BASE = REPO / "research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs"
RUNS = ["off_tp", "off_fp", "on_tp", "on_fp"]
OUT_BASE = REPO / "research/paired_priority_bug_audit/v6c_phaseB_v5bam_runs"


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


def run_summary(run_data: dict[str, dict[str, dict]], label: str, off_run: str, on_run: str):
    """Returns list of (region_id, ng_off, ng_on, n_reads, somatic_tag_pct_off)."""
    out = []
    for region_id in sorted(set(run_data[off_run]) | set(run_data[on_run])):
        off = run_data[off_run].get(region_id)
        on = run_data[on_run].get(region_id)
        if not off or not on:
            continue
        ng_off = compute_ng(off["hp_counter"])
        ng_on = compute_ng(on["hp_counter"])
        n_reads = off["n_reads"]
        n_som = sum(v for k, v in off["hp_counter"].items() if k in {"1-1", "2-1", "3", "11", "21", "33"})
        out.append({
            "region_id": region_id,
            "label": label,
            "ng_off": ng_off,
            "ng_on": ng_on,
            "n_reads": n_reads,
            "somatic_pct_off": round(100 * n_som / max(n_reads, 1), 2),
        })
    return out


def aggregate_marker(rows, ng_off_thresh: int = 3):
    """Aggregate marker filter TP rate (NG_off ≥ threshold) and corresponding flag=on cell."""
    tp_off = sum(1 for r in rows if r["label"] == "TP" and r["ng_off"] >= ng_off_thresh)
    fp_off = sum(1 for r in rows if r["label"] == "FP" and r["ng_off"] >= ng_off_thresh)
    tp_on_max = sum(
        1 for r in rows if r["label"] == "TP" and r["ng_on"] == max((rr["ng_on"] for rr in rows), default=0)
    )
    # The flag=on equivalent cell = NG_on=2 (since schema collapse caps NG_on at 2)
    tp_on2 = sum(1 for r in rows if r["label"] == "TP" and r["ng_on"] == 2)
    fp_on2 = sum(1 for r in rows if r["label"] == "FP" and r["ng_on"] == 2)
    return {
        "marker_off_tp": tp_off,
        "marker_off_fp": fp_off,
        "marker_off_rate": tp_off / max(tp_off + fp_off, 1),
        "ng_on2_tp": tp_on2,
        "ng_on2_fp": fp_on2,
        "ng_on2_rate": tp_on2 / max(tp_on2 + fp_on2, 1),
    }


def main() -> int:
    print("[V3F vs V5 BAM head-to-head] Loading runs...")
    v3f_data = {run: collect(V3F_BASE / run) for run in RUNS}
    v5_data = {run: collect(V5_BASE / run) for run in RUNS}

    for tag, data in [("V3F", v3f_data), ("V5", v5_data)]:
        for run in RUNS:
            print(f"  {tag} {run}: {len(data[run])} regions")

    v3f_rows = run_summary(v3f_data, "TP", "off_tp", "on_tp") + run_summary(v3f_data, "FP", "off_fp", "on_fp")
    v5_rows = run_summary(v5_data, "TP", "off_tp", "on_tp") + run_summary(v5_data, "FP", "off_fp", "on_fp")

    # Read-level hp distribution comparison
    print("\n=== chr19 全部 reads hp distribution (V3F vs V5) ===")
    print(f"{'hp_value':<10} {'V3F off':>10} {'V3F on':>10} {'V5 off':>10} {'V5 on':>10}")
    all_hp = set()
    for run in RUNS:
        for d in [v3f_data, v5_data]:
            for v in d[run].values():
                all_hp.update(v["hp_counter"].keys())
    all_hp = sorted(all_hp)
    for hp in all_hp:
        v3f_off = sum(d["hp_counter"].get(hp, 0) for d in v3f_data["off_tp"].values()) + sum(
            d["hp_counter"].get(hp, 0) for d in v3f_data["off_fp"].values()
        )
        v3f_on = sum(d["hp_counter"].get(hp, 0) for d in v3f_data["on_tp"].values()) + sum(
            d["hp_counter"].get(hp, 0) for d in v3f_data["on_fp"].values()
        )
        v5_off = sum(d["hp_counter"].get(hp, 0) for d in v5_data["off_tp"].values()) + sum(
            d["hp_counter"].get(hp, 0) for d in v5_data["off_fp"].values()
        )
        v5_on = sum(d["hp_counter"].get(hp, 0) for d in v5_data["on_tp"].values()) + sum(
            d["hp_counter"].get(hp, 0) for d in v5_data["on_fp"].values()
        )
        print(f"{hp:<10} {v3f_off:>10} {v3f_on:>10} {v5_off:>10} {v5_on:>10}")

    # Marker filter aggregate
    print("\n=== Marker filter (NG_off ≥ 3) — V3F vs V5 ===")
    v3f_agg = aggregate_marker(v3f_rows)
    v5_agg = aggregate_marker(v5_rows)
    print(f"{'BAM':<6} {'NG_off≥3 TP':>13} {'FP':>5} {'rate':>7} | {'NG_on=2 TP':>12} {'FP':>5} {'rate':>7}")
    for tag, agg in [("V3F", v3f_agg), ("V5", v5_agg)]:
        print(
            f"{tag:<6} {agg['marker_off_tp']:>13} {agg['marker_off_fp']:>5} "
            f"{agg['marker_off_rate']:>7.3f} | "
            f"{agg['ng_on2_tp']:>12} {agg['ng_on2_fp']:>5} {agg['ng_on2_rate']:>7.3f}"
        )

    # Per-cell TP rate cross-tab (V3F vs V5)
    for tag, rows in [("V3F", v3f_rows), ("V5", v5_rows)]:
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

    # Per-region V3F-vs-V5 NG_off difference (does V5 give DIFFERENT NG count?)
    v3f_region_ng: dict[tuple[str, str], int] = {(r["region_id"], r["label"]): r["ng_off"] for r in v3f_rows}
    v5_region_ng: dict[tuple[str, str], int] = {(r["region_id"], r["label"]): r["ng_off"] for r in v5_rows}
    common = set(v3f_region_ng) & set(v5_region_ng)
    diff_count: Counter[tuple[int, int]] = Counter()
    for k in common:
        diff_count[(v3f_region_ng[k], v5_region_ng[k])] += 1
    print(f"\n=== Per-region NG_off (V3F BAM vs V5 BAM) — flag=off ===")
    print(f"{'V3F NG':>7} {'V5 NG':>7} {'count':>7}")
    for (v3f_n, v5_n), c in sorted(diff_count.items()):
        marker = " ← DIFF" if v3f_n != v5_n else ""
        print(f"{v3f_n:>7} {v5_n:>7} {c:>7}{marker}")
    n_diff = sum(c for (v3f_n, v5_n), c in diff_count.items() if v3f_n != v5_n)
    n_total = sum(diff_count.values())
    print(f"\n  Disagreement: {n_diff}/{n_total} ({100 * n_diff / max(n_total, 1):.1f}%)")

    # Save summary TSV
    out_tsv = OUT_BASE / "v3f_vs_v5_summary.tsv"
    with out_tsv.open("w", newline="") as fh:
        fields = [
            "BAM",
            "marker_off_tp",
            "marker_off_fp",
            "marker_off_rate",
            "ng_on2_tp",
            "ng_on2_fp",
            "ng_on2_rate",
        ]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"BAM": "V3F", **v3f_agg})
        writer.writerow({"BAM": "V5", **v5_agg})
    print(f"\n[V3F vs V5] summary → {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
