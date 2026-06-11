#!/usr/bin/env python3
"""
A1: HP 6-value distribution scan across HCC1395 baseline + HCC1395 V6 + V6 5-sample extension.

Background:
  PI 4-goal report requires baseline vs V6 HP-tag distribution comparison.
  longphase-to BAM tagging encodes HP into 6 categories per Errata 01:
    HP=0  : unset/unphased (default haplotag absent)
    HP=1  : assigned to haplotype 1
    HP=2  : assigned to haplotype 2
    HP=11 : Layer 1.5 forced HP1 (V6 priority-bug-fixed assignment)
    HP=21 : Layer 1.5 forced HP2
    HP=33 : Layer 1.5 ambiguous (germline-absent, both haps tied)

Scope:
  chr8 + chr19 only (known hotspots; chr8 = LOH+HPSig 7.4x FP enrichment,
  chr19 = priority-bug-sensitive germline-rich region).
  Whole-genome scan would take >4h per BAM @ 250-470GB; chr8+chr19 ~10-15min/BAM.

Outputs:
  - TSV : A1_HP_6values_5sample.tsv
  - PNG : figures/A1_HP_distribution.png (stacked bar, HP fraction by sample x version)

Sources (frozen 2026-05-21):
  - HCC1395 baseline : /big7_disk/.../baseline/tumor_tagged.bam            (278 GB, 2026-04-03)
  - HCC1395 V6       : /big7_disk/.../v6_germline_absent_revert/tumor_tagged.bam (287 GB, 2026-05-10)
  - HCC1937 V6       : /big7_disk/.../v6_5sample_extension/HCC1937/tumor_tagged.bam (472 GB, 2026-05-11)
  - HCC1954 V6       : /big7_disk/.../v6_5sample_extension/HCC1954/tumor_tagged.bam (253 GB, 2026-05-11)
  - H1437 V6         : /big7_disk/.../v6_5sample_extension/H1437/tumor_tagged.bam   (243 GB, 2026-05-10)
  - H2009 V6         : /big7_disk/.../v6_5sample_extension/H2009/tumor_tagged.bam   (327 GB, 2026-05-11)

  HCC1395 baseline 4-sample BAMs do NOT exist (per CURRENT_FOCUS + Errata 01
  Phase 2 plan v2.0 NEGATIVE flag). Only HCC1395 baseline available for cross-version compare.
"""

from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path
from collections import Counter

import pysam
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Canonical HP value categories (per Errata 01)
HP_CATEGORIES = [0, 1, 2, 11, 21, 33]

BAMS = [
    # (sample, bam_version, path)
    ("HCC1395", "baseline",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam"),
    ("HCC1395", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam"),
    ("HCC1937", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1937/tumor_tagged.bam"),
    ("HCC1954", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1954/tumor_tagged.bam"),
    ("H1437", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H1437/tumor_tagged.bam"),
    ("H2009", "V6",
     "/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H2009/tumor_tagged.bam"),
]

SCAN_CHROMS = ["chr8", "chr19"]


def scan_bam(bam_path: str, chroms: list[str]) -> dict:
    """Scan HP tag distribution across given chromosomes. Returns counter + total.

    Skips unmapped/secondary/supplementary alignments (we want primary tagged reads only).
    """
    counts: Counter = Counter()
    total = 0
    t0 = time.time()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom in chroms:
            try:
                iterator = bam.fetch(chrom, multiple_iterators=False)
            except ValueError as e:
                print(f"  [WARN] {bam_path} fetch({chrom}) failed: {e}", file=sys.stderr)
                continue
            for read in iterator:
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                total += 1
                try:
                    hp = read.get_tag("HP")
                except KeyError:
                    hp = -1  # untagged
                counts[hp] += 1
    elapsed = time.time() - t0
    return {"counts": dict(counts), "total": total, "elapsed_sec": elapsed}


def summarize(scan: dict) -> dict:
    counts = scan["counts"]
    total = scan["total"]
    row = {f"HP_{k}_count": counts.get(k, 0) for k in HP_CATEGORIES}
    # everything not in HP_CATEGORIES and not the -1 untagged sentinel
    other = sum(v for k, v in counts.items() if k not in HP_CATEGORIES and k != -1)
    untagged = counts.get(-1, 0)
    row["HP_untagged_count"] = untagged
    row["HP_other_count"] = other
    row["total"] = total
    hp1 = counts.get(1, 0)
    hp2 = counts.get(2, 0)
    row["HP1_HP2_ratio"] = round(hp1 / hp2, 4) if hp2 > 0 else float("nan")
    row["HP33_pct"] = round(100.0 * counts.get(33, 0) / total, 4) if total > 0 else 0.0
    return row


def write_tsv(rows: list[dict], out_path: Path) -> None:
    cols = [
        "sample", "bam_version", "scan_scope",
        "HP_0_count", "HP_1_count", "HP_2_count",
        "HP_11_count", "HP_21_count", "HP_33_count",
        "HP_untagged_count", "HP_other_count",
        "total", "HP1_HP2_ratio", "HP33_pct",
        "elapsed_sec", "bam_path",
    ]
    with open(out_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def plot_stacked(rows: list[dict], out_png: Path) -> None:
    """Stacked bar chart: HP fraction per (sample, bam_version)."""
    labels = [f"{r['sample']}\n{r['bam_version']}" for r in rows]
    totals = np.array([r["total"] for r in rows], dtype=float)
    totals[totals == 0] = 1.0  # avoid div0

    cats = HP_CATEGORIES + ["untagged", "other"]
    colors = {
        0: "#cccccc",   # gray - unset
        1: "#1f77b4",   # blue - HP1
        2: "#ff7f0e",   # orange - HP2
        11: "#2ca02c",  # green - Layer 1.5 HP1
        21: "#d62728",  # red - Layer 1.5 HP2
        33: "#9467bd",  # purple - Layer 1.5 ambiguous
        "untagged": "#e0e0e0",
        "other": "#000000",
    }

    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = np.zeros(len(rows))
    x = np.arange(len(rows))
    for cat in cats:
        if cat in HP_CATEGORIES:
            vals = np.array([r[f"HP_{cat}_count"] for r in rows], dtype=float)
            label = f"HP={cat}"
        elif cat == "untagged":
            vals = np.array([r["HP_untagged_count"] for r in rows], dtype=float)
            label = "untagged"
        else:
            vals = np.array([r["HP_other_count"] for r in rows], dtype=float)
            label = "other"
        fracs = vals / totals
        ax.bar(x, fracs, bottom=bottom, color=colors[cat], label=label, edgecolor="white", linewidth=0.5)
        bottom += fracs

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fraction of primary reads")
    ax.set_ylim(0, 1.0)
    ax.set_title("A1: HP 6-value distribution per sample x BAM version\n"
                 "Scope: chr8 + chr19 primary alignments (HCC1395 baseline available, 4 V6 ext)")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


def load_state(state_path: Path) -> dict:
    if state_path.exists():
        with open(state_path) as f:
            return json.load(f)
    return {"scan_scope": "+".join(SCAN_CHROMS),
            "categories": HP_CATEGORIES,
            "rows": [], "raw": []}


def save_state(state: dict, state_path: Path) -> None:
    tmp = state_path.with_suffix(".json.tmp")
    with open(tmp, "w") as f:
        json.dump(state, f, indent=2)
    tmp.replace(state_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-png", required=True)
    parser.add_argument("--out-json", required=True,
                        help="Per-BAM raw counts JSON (incremental, supports resume)")
    parser.add_argument("--only-missing", action="store_true",
                        help="Skip BAMs already present in --out-json state")
    parser.add_argument("--render-only", action="store_true",
                        help="Skip scanning; regenerate TSV+PNG from --out-json state")
    args = parser.parse_args()

    out_tsv = Path(args.out_tsv)
    out_png = Path(args.out_png)
    out_json = Path(args.out_json)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    scope = "+".join(SCAN_CHROMS)
    state = load_state(out_json)
    done_keys = {(r["sample"], r["bam_version"]) for r in state["rows"]}
    t_start = time.time()

    if not args.render_only:
        for sample, version, bam_path in BAMS:
            key = (sample, version)
            if args.only_missing and key in done_keys:
                print(f"[skip] {sample}/{version} (already in state)", flush=True)
                continue
            print(f"[scan] {sample}/{version} {bam_path}", flush=True)
            if not Path(bam_path).exists():
                print(f"  [SKIP] BAM not found", file=sys.stderr)
                continue
            scan = scan_bam(bam_path, SCAN_CHROMS)
            summary = summarize(scan)
            summary["sample"] = sample
            summary["bam_version"] = version
            summary["scan_scope"] = scope
            summary["bam_path"] = bam_path
            summary["elapsed_sec"] = round(scan["elapsed_sec"], 1)
            # Replace any prior entry for same (sample,version), then append
            state["rows"] = [r for r in state["rows"] if (r["sample"], r["bam_version"]) != key]
            state["raw"] = [r for r in state["raw"] if (r["sample"], r["version"]) != key]
            state["rows"].append(summary)
            state["raw"].append({"sample": sample, "version": version, "bam_path": bam_path,
                                 "scan_scope": scope, "raw_counts": scan["counts"],
                                 "total": scan["total"], "elapsed_sec": scan["elapsed_sec"]})
            save_state(state, out_json)  # incremental persistence
            print(f"  -> total={scan['total']:,} HP1:HP2={summary['HP1_HP2_ratio']} "
                  f"HP33%={summary['HP33_pct']} elapsed={scan['elapsed_sec']:.1f}s",
                  flush=True)

    # Sort rows to canonical order matching BAMS for stable TSV + plot
    order = {(s, v): i for i, (s, v, _) in enumerate(BAMS)}
    state["rows"].sort(key=lambda r: order.get((r["sample"], r["bam_version"]), 999))
    save_state(state, out_json)
    write_tsv(state["rows"], out_tsv)
    plot_stacked(state["rows"], out_png)
    print(f"\n[done] TSV  -> {out_tsv}")
    print(f"[done] PNG  -> {out_png}")
    print(f"[done] JSON -> {out_json}")
    print(f"[done] total elapsed = {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
