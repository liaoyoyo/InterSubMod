#!/usr/bin/env python3
"""Phase C: HCC1395 全基因組 V3F vs V5 vs V6 three-way head-to-head comparison.

Reads 12 ISM run dirs from phaseC_genome_three_way/ (3 BAM × 2 flag × 2 label).
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
BASE = REPO / "research/paired_priority_bug_audit/phaseC_genome_three_way"
BAMS = ["V3F", "V5", "V6"]
RUNS = ["off_tp", "off_fp", "on_tp", "on_fp"]


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


def main() -> int:
    print("[Phase C 全基因組三向] Loading 12 runs...")
    data: dict[str, dict[str, dict[str, dict]]] = {}
    for tag in BAMS:
        data[tag] = {}
        for run in RUNS:
            run_dir = BASE / f"{tag}_{run}"
            if not run_dir.exists():
                print(f"  [SKIP] {tag}_{run}: {run_dir} missing")
                continue
            data[tag][run] = collect(run_dir)
            print(f"  {tag}_{run}: {len(data[tag][run])} regions")

    available = [t for t in BAMS if all(r in data.get(t, {}) for r in RUNS)]
    print(f"\n=== Available BAMs (complete 4 runs): {available} ===\n")

    # Read-level distribution
    print("=== 全基因組 hp distribution (3 BAMs) ===")
    print(f"{'hp':<8}", end="")
    for tag in available:
        print(f" {tag+' off':>13} {tag+' on':>13}", end="")
    print()
    all_hp = set()
    for tag in available:
        for run in RUNS:
            for v in data[tag][run].values():
                all_hp.update(v["hp_counter"].keys())
    for hp in sorted(all_hp):
        print(f"{hp:<8}", end="")
        for tag in available:
            off_c = sum(d["hp_counter"].get(hp, 0) for d in data[tag]["off_tp"].values()) + sum(
                d["hp_counter"].get(hp, 0) for d in data[tag]["off_fp"].values()
            )
            on_c = sum(d["hp_counter"].get(hp, 0) for d in data[tag]["on_tp"].values()) + sum(
                d["hp_counter"].get(hp, 0) for d in data[tag]["on_fp"].values()
            )
            print(f" {off_c:>13} {on_c:>13}", end="")
        print()

    # Priority bug feature 化指標
    print("\n=== hp=1-1 vs hp=2-1 ratio (priority bug feature 化) ===")
    print(f"{'BAM':<6} {'hp=1-1':>12} {'hp=2-1':>12} {'ratio':>10}")
    ratios: dict[str, float] = {}
    for tag in available:
        h11 = sum(d["hp_counter"].get("1-1", 0) for d in data[tag]["off_tp"].values()) + sum(
            d["hp_counter"].get("1-1", 0) for d in data[tag]["off_fp"].values()
        )
        h21 = sum(d["hp_counter"].get("2-1", 0) for d in data[tag]["off_tp"].values()) + sum(
            d["hp_counter"].get("2-1", 0) for d in data[tag]["off_fp"].values()
        )
        ratio = h11 / max(h21, 1)
        ratios[tag] = ratio
        print(f"{tag:<6} {h11:>12} {h21:>12} {ratio:>10.3f}")

    # hp=3 ambiguous
    print("\n=== hp=3 (somatic ambiguous) reads ===")
    print(f"{'BAM':<6} {'hp=3':>12}")
    h3_dict: dict[str, int] = {}
    for tag in available:
        h3 = sum(d["hp_counter"].get("3", 0) for d in data[tag]["off_tp"].values()) + sum(
            d["hp_counter"].get("3", 0) for d in data[tag]["off_fp"].values()
        )
        h3_dict[tag] = h3
        print(f"{tag:<6} {h3:>12}")

    # Marker filter aggregate
    summary_rows = []
    print("\n=== Marker filter (NG_off ≥ 3) — three-way ===")
    print(f"{'BAM':<6} {'NG≥3 N':>10} {'TP':>6} {'FP':>6} {'rate':>8} | {'NG_on=2 TP':>11} {'FP':>6} {'rate':>8}")
    for tag in available:
        rows = []
        for label, off_run, on_run in [("TP", "off_tp", "on_tp"), ("FP", "off_fp", "on_fp")]:
            for region_id in sorted(set(data[tag][off_run]) | set(data[tag][on_run])):
                off = data[tag][off_run].get(region_id)
                on = data[tag][on_run].get(region_id)
                if not off or not on:
                    continue
                rows.append({
                    "label": label,
                    "ng_off": compute_ng(off["hp_counter"]),
                    "ng_on": compute_ng(on["hp_counter"]),
                })
        tp_off = sum(1 for r in rows if r["label"] == "TP" and r["ng_off"] >= 3)
        fp_off = sum(1 for r in rows if r["label"] == "FP" and r["ng_off"] >= 3)
        tp_on2 = sum(1 for r in rows if r["label"] == "TP" and r["ng_on"] == 2)
        fp_on2 = sum(1 for r in rows if r["label"] == "FP" and r["ng_on"] == 2)
        rate_off = tp_off / max(tp_off + fp_off, 1)
        rate_on = tp_on2 / max(tp_on2 + fp_on2, 1)
        print(
            f"{tag:<6} {tp_off+fp_off:>10} {tp_off:>6} {fp_off:>6} {rate_off:>8.3f} | "
            f"{tp_on2:>11} {fp_on2:>6} {rate_on:>8.3f}"
        )
        summary_rows.append({
            "BAM": tag,
            "marker_off_n": tp_off + fp_off,
            "marker_off_tp": tp_off,
            "marker_off_fp": fp_off,
            "marker_off_rate": rate_off,
            "ng_on2_tp": tp_on2,
            "ng_on2_fp": fp_on2,
            "ng_on2_rate": rate_on,
            "hp_1_1_2_1_ratio": ratios.get(tag, 0),
            "hp_3_count": h3_dict.get(tag, 0),
        })

    # Phase C 驗收
    if "V6" in available and "V3F" in available:
        print("\n=== Phase C 驗收 ===")
        v6_summary = next(s for s in summary_rows if s["BAM"] == "V6")
        v3f_summary = next(s for s in summary_rows if s["BAM"] == "V3F")
        print(
            f"  全基因組 hp=1-1:hp=2-1 ratio  V3F={ratios['V3F']:.3f}  V6={ratios['V6']:.3f}  "
            f"→ {'PASS V6 接近 V3F' if abs(ratios['V6'] - ratios['V3F']) / max(ratios['V3F'], 0.1) < 0.3 else 'FAIL'}"
        )
        cov_ratio = v6_summary["marker_off_n"] / max(v3f_summary["marker_off_n"], 1)
        print(
            f"  marker coverage ratio (V6/V3F) = {cov_ratio:.3f}  "
            f"→ {'PASS V6 ≥ V3F × 0.95' if cov_ratio >= 0.95 else 'FAIL'}"
        )

    # Save TSV
    out_tsv = BASE / "v3f_vs_v5_vs_v6_genome_summary.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        if summary_rows:
            writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
    print(f"\n[Phase C three-way] summary → {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
