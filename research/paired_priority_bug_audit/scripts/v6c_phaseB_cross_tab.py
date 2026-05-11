#!/usr/bin/env python3
"""V6-C Phase B cross-flag NG cross-tab analysis.

Reads per-region reads.tsv from 4 ISM runs (off_tp, off_fp, on_tp, on_fp)
to compute number-of-distinct-hp-groups (NG) per region under each flag
setting, then evaluates whether HPFineNGroups marker behavior changes.

Output: research/paired_priority_bug_audit/v6c_phaseB_runs/cross_tab_summary.tsv
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
RUN_BASE = REPO / "research/paired_priority_bug_audit/v6c_phaseB_runs"
RUNS = ["off_tp", "off_fp", "on_tp", "on_fp"]


def collect_regions(run_dir: Path) -> dict[str, dict]:
    """Return {region_id: {hp_counter, n_reads, region_path}} for run."""
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
        out[region_id] = {
            "hp_counter": hp_counter,
            "n_reads": sum(hp_counter.values()),
            "region_path": str(region_dir.relative_to(REPO)),
        }
    return out


def compute_ng(hp_counter: Counter[str]) -> int:
    """NG = number of distinct hp buckets EXCLUDING hp='0' (unphased) and 'hp' (header noise)."""
    return sum(1 for k, v in hp_counter.items() if k not in {"0", "hp"} and v > 0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(RUN_BASE / "cross_tab_summary.tsv"))
    args = parser.parse_args()

    print("[V6-C Phase B cross-tab] Loading 4 runs...")
    data: dict[str, dict[str, dict]] = {}
    for run in RUNS:
        data[run] = collect_regions(RUN_BASE / run)
        print(f"  {run}: {len(data[run])} regions")

    # Region universe: union of TP/FP × off/on (should match per label)
    tp_regions = set(data["off_tp"]) | set(data["on_tp"])
    fp_regions = set(data["off_fp"]) | set(data["on_fp"])
    print(f"  TP regions: {len(tp_regions)}, FP regions: {len(fp_regions)}")

    # Build per-region cross-flag table
    rows: list[dict] = []
    for label, regions, off_run, on_run in [
        ("TP", tp_regions, "off_tp", "on_tp"),
        ("FP", fp_regions, "off_fp", "on_fp"),
    ]:
        for region_id in sorted(regions):
            off_info = data[off_run].get(region_id)
            on_info = data[on_run].get(region_id)
            if not off_info or not on_info:
                continue
            ng_off = compute_ng(off_info["hp_counter"])
            ng_on = compute_ng(on_info["hp_counter"])
            n_reads = off_info["n_reads"]
            n_somatic_off = sum(
                v for k, v in off_info["hp_counter"].items() if k in {"1-1", "2-1", "3", "11", "21", "33"}
            )
            rows.append(
                {
                    "region_id": region_id,
                    "label": label,
                    "n_reads": n_reads,
                    "ng_off": ng_off,
                    "ng_on": ng_on,
                    "ng_delta": ng_on - ng_off,
                    "n_somatic_tag_reads_off": n_somatic_off,
                    "somatic_tag_pct_off": round(100.0 * n_somatic_off / max(n_reads, 1), 2),
                    "hp_off_distribution": ",".join(
                        f"{k}:{v}" for k, v in sorted(off_info["hp_counter"].items())
                    ),
                    "hp_on_distribution": ",".join(
                        f"{k}:{v}" for k, v in sorted(on_info["hp_counter"].items())
                    ),
                }
            )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        if rows:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    print(f"\n[V6-C Phase B cross-tab] wrote {len(rows)} rows → {out_path}")

    # Aggregate: NG_off × NG_on cross-tab per label
    print("\n=== NG_off × NG_on cross-tab (TP / FP) ===")
    cross: dict[tuple[str, int, int], int] = defaultdict(int)
    n_per_label: Counter[str] = Counter()
    for r in rows:
        cross[(r["label"], r["ng_off"], r["ng_on"])] += 1
        n_per_label[r["label"]] += 1
    for label in ["TP", "FP"]:
        print(f"\n--- {label} (n={n_per_label[label]}) ---")
        ng_offs = sorted({k[1] for k in cross if k[0] == label})
        ng_ons = sorted({k[2] for k in cross if k[0] == label})
        print("NG_off \\ NG_on  " + "  ".join(f"{x:>5}" for x in ng_ons) + "    total")
        for noff in ng_offs:
            row_total = sum(cross[(label, noff, non)] for non in ng_ons)
            row_str = f"  NG_off={noff}     " + "  ".join(
                f"{cross[(label, noff, non)]:>5}" for non in ng_ons
            ) + f"    {row_total}"
            print(row_str)

    # Marker filter analysis: NG=4 in off, what becomes in on?
    print("\n=== Marker filter (NG_off ≥ 3) → NG_on distribution ===")
    for label in ["TP", "FP"]:
        marker = [r for r in rows if r["label"] == label and r["ng_off"] >= 3]
        if not marker:
            print(f"  {label}: 0 regions with NG_off ≥ 3")
            continue
        ng_on_dist = Counter(r["ng_on"] for r in marker)
        print(f"  {label} (n={len(marker)}): NG_on distribution → {dict(sorted(ng_on_dist.items()))}")

    # TP rate stratified by NG_off and NG_on
    print("\n=== TP rate per NG_off × NG_on cell ===")
    cell_tp: dict[tuple[int, int], int] = defaultdict(int)
    cell_fp: dict[tuple[int, int], int] = defaultdict(int)
    for r in rows:
        key = (r["ng_off"], r["ng_on"])
        if r["label"] == "TP":
            cell_tp[key] += 1
        else:
            cell_fp[key] += 1
    cells = sorted(set(cell_tp) | set(cell_fp))
    print(f"{'NG_off':>7} {'NG_on':>6} {'TP':>6} {'FP':>6} {'TP_rate':>9}")
    for noff, non in cells:
        tp_n = cell_tp[(noff, non)]
        fp_n = cell_fp[(noff, non)]
        rate = tp_n / max(tp_n + fp_n, 1)
        print(f"{noff:>7} {non:>6} {tp_n:>6} {fp_n:>6} {rate:>9.3f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
