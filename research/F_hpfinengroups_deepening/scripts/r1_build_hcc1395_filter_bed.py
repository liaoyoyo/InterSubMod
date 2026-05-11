#!/usr/bin/env python3
"""R1 — Build HCC1395 filter BED for Normal BAM region-subset extraction.

Selects HCC1395 regions that match the canonical F pilot filter:
  NG (HPFineNGroups) >= 4 AND caller_af < 0.4 AND NumReads >= 80 AND NOT LOH

Output BED is used by samtools view -L to extract only these regions from the
136GB HCC1395 normal BAM, enabling a fast Phase 2 Normal BAM pilot without
transferring the full file.

Output:
  - data/r1_hcc1395_filter_regions.bed  (chr, start, end)
  - data/r1_hcc1395_filter_summary.tsv  (per-mode count + overlap stats)
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

MASTER = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
)
OUT_DIR = Path(__file__).resolve().parents[1] / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FLANK = 5000  # half-window around Pos; conservative for Normal ASM coverage
USE_COLS = [
    "sample", "mode", "Chr", "Pos", "HPFineNGroups", "caller_af",
    "NumReads", "Potential_LOH", "to_loh_bed_hit", "tool_potential_loh",
    "core_loh_like", "truth_label",
]


def main():
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", usecols=USE_COLS, low_memory=False)
    for c in ("to_loh_bed_hit", "tool_potential_loh", "core_loh_like"):
        df[c] = df[c].astype(str).str.lower().isin({"true", "1", "t"})
    df["caller_af"] = pd.to_numeric(df["caller_af"], errors="coerce")
    df["HPFineNGroups"] = pd.to_numeric(df["HPFineNGroups"], errors="coerce")
    df["NumReads"] = pd.to_numeric(df["NumReads"], errors="coerce")

    sub = df[df["sample"] == "HCC1395"].copy()
    print(f"[R1] HCC1395 total rows: {len(sub):,}")

    filt = sub[
        (sub["HPFineNGroups"] >= 4)
        & (sub["caller_af"] < 0.4)
        & (sub["NumReads"] >= 80)
    ].copy()
    print(f"[R1] After NG>=4 & AF<0.4 & NR>=80: {len(filt):,}")

    # NonLOH: mode-aware
    nonloh_to = (filt["mode"] == "to") & (~filt["to_loh_bed_hit"])
    nonloh_pa = (filt["mode"] == "paired") & (~filt["tool_potential_loh"]) & (~filt["core_loh_like"])
    filt_nl = filt[nonloh_to | nonloh_pa].copy()
    print(f"[R1] After NonLOH: {len(filt_nl):,}")

    # Summary per mode
    summary_rows = []
    for mode, g in filt_nl.groupby("mode"):
        tp = int((g["truth_label"].str.upper() == "TP").sum())
        fp = int((g["truth_label"].str.upper() == "FP").sum())
        summary_rows.append({
            "mode": mode, "regions": len(g), "TP": tp, "FP": fp,
            "TP_rate": tp / max(tp + fp, 1),
        })
    summary = pd.DataFrame(summary_rows)
    summary_path = OUT_DIR / "r1_hcc1395_filter_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"[R1] Summary written: {summary_path}")
    print(summary.to_string(index=False))

    # Merge to unique loci (union of TO + paired)
    uniq = (filt_nl[["Chr", "Pos"]].drop_duplicates().sort_values(["Chr", "Pos"]))
    print(f"[R1] Unique loci (TO+paired): {len(uniq):,}")

    # Build BED with flank
    bed = uniq.copy()
    bed["start"] = (bed["Pos"] - FLANK).clip(lower=0).astype(int)
    bed["end"] = (bed["Pos"] + FLANK).astype(int)
    bed = bed[["Chr", "start", "end"]]

    # Sort + merge overlaps
    bed = bed.sort_values(["Chr", "start"]).reset_index(drop=True)
    merged = []
    cur = None
    for row in bed.itertuples(index=False):
        if cur and cur[0] == row.Chr and row.start <= cur[2]:
            cur = (cur[0], cur[1], max(cur[2], row.end))
        else:
            if cur:
                merged.append(cur)
            cur = (row.Chr, row.start, row.end)
    if cur:
        merged.append(cur)

    bed_path = OUT_DIR / "r1_hcc1395_filter_regions.bed"
    with bed_path.open("w") as fh:
        for chrom, start, end in merged:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    print(f"[R1] BED written: {bed_path}  merged_regions={len(merged):,}  flank=±{FLANK}")

    # Report total bp covered (approximate BAM extraction size estimate)
    total_bp = sum(e - s for _, s, e in merged)
    print(f"[R1] Total genomic span: {total_bp:,} bp  (~{total_bp/1e6:.1f} Mb)")


if __name__ == "__main__":
    main()
