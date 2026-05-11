"""Build TP-only master dataset from 14-combo KDE rerun outputs.

Reads 14 significance_summary.csv files and concatenates them with sample/mode
tags to produce a single TSV.gz compatible with downstream scripts that expect
the `all_region_rows` schema.

Output: kde_rerun_B_14combos/all_region_rows_kde_B_tp.tsv.gz
"""
import glob
import os
import sys

import pandas as pd

ROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/kde_rerun_B_14combos"
OUT = f"{ROOT}/all_region_rows_kde_B_tp.tsv.gz"

# Columns required by downstream scripts (match stale all_region_rows.tsv.gz schema subset)
# Stale schema uses Chr/Pos (no Chromosome/Start/End); significance_summary.csv uses Chr/Pos too.
COLS = [
    "RegionID", "Chr", "Pos",
    "NumReads", "Coverage_Multiple", "Diploid_Coverage_Used",
    "Coverage_Category",
    # QS + downstream feature columns (added 2026-04-21 for P1/P2/P3 verification)
    "Quality_Score", "Quality_Tier", "VerificationClass",
    "HeuristicScore", "PassedGating",
    "HPFineNGroups", "HPFineSig", "AlleleDelta",
    "Potential_LOH", "LOH_Subtype", "LOH_Bed_Overlap",
    "Subclone_ID", "Significant",
]

frames = []
for combo_dir in sorted(glob.glob(f"{ROOT}/*_tp")):
    combo = os.path.basename(combo_dir).replace("_tp", "")
    parts = combo.rsplit("_", 2)
    sample = parts[0]
    mode = "_".join(parts[1:]) if len(parts) >= 3 else "?"
    csv = os.path.join(combo_dir, "significance_summary.csv")
    if not os.path.exists(csv):
        print(f"SKIP {combo}: no significance_summary.csv", file=sys.stderr)
        continue
    df = pd.read_csv(csv)
    # Pick available columns (handle schema drift)
    keep = [c for c in COLS if c in df.columns]
    df = df[keep].copy()
    df["sample"] = sample
    df["mode"] = mode
    df["truth_label"] = "TP"  # 14 combos ran with filtered_snv_tp.vcf.gz
    frames.append(df)
    print(f"OK  {combo:40s} n={len(df):>7,}  baseline={df['Diploid_Coverage_Used'].iloc[0]:.0f}x")

if not frames:
    print("ERROR: no combos found", file=sys.stderr)
    sys.exit(1)

master = pd.concat(frames, ignore_index=True)

# Sanity: per-sample baseline uniqueness
sanity = master.groupby(["sample", "mode"])["Diploid_Coverage_Used"].nunique()
bad = sanity[sanity > 1]
if len(bad) > 0:
    print(f"WARNING: non-unique baseline per sample/mode: {bad}", file=sys.stderr)

master.to_csv(OUT, sep="\t", index=False, compression="gzip")
print(f"\nWrote {OUT}")
print(f"  rows={len(master):,}  combos={master.groupby(['sample','mode']).ngroups}")
print(f"  columns={list(master.columns)}")
