#!/usr/bin/env python3
"""Build subset VCFs for R1 HCC1395 Normal BAM pilot.

Strategy: extract filter SNV positions from master dataset (HCC1395,
NG=4+AF<0.4+NR>=80+NonLOH) per mode, build subset VCFs (TP+FP merged)
for fast Phase 2 pipeline execution.

Outputs:
  research/F_hpfinengroups_deepening/data/r1_hcc1395_filter_snvs_TO.tsv
  research/F_hpfinengroups_deepening/data/r1_hcc1395_filter_snvs_paired.tsv
"""
import gzip
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
)
OUT_DIR = ROOT / "research/F_hpfinengroups_deepening/data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COLS = [
    "sample", "mode", "Chr", "Pos", "HPFineNGroups", "caller_af",
    "NumReads", "to_loh_bed_hit", "tool_potential_loh", "core_loh_like",
    "truth_label",
]

print("Loading master dataset (HCC1395 only)...")
iter_reader = pd.read_csv(
    MASTER, sep="\t", compression="gzip", usecols=COLS, chunksize=100_000
)
parts = []
for chunk in iter_reader:
    sub = chunk[chunk["sample"] == "HCC1395"]
    if len(sub) > 0:
        parts.append(sub)

df = pd.concat(parts, ignore_index=True)
print(f"HCC1395 total: {len(df)}")

# Apply F pilot canonical filter
mask = (
    (df["HPFineNGroups"] >= 4)
    & (df["caller_af"] < 0.4)
    & (df["NumReads"] >= 80)
)
filt = df[mask].copy()
print(f"After NG>=4 + AF<0.4 + NR>=80: {len(filt)}")

# Mode-aware NonLOH
nonloh_to = (filt["mode"] == "to") & (~filt["to_loh_bed_hit"])
nonloh_pa = (
    (filt["mode"] == "paired")
    & (~filt["tool_potential_loh"])
    & (~filt["core_loh_like"])
)
filt = filt[nonloh_to | nonloh_pa].copy()
print(f"After NonLOH filter: {len(filt)}")

for mode in ["to", "paired"]:
    sub = filt[filt["mode"] == mode][["Chr", "Pos", "truth_label"]].copy()
    sub["truth_label_int"] = (sub["truth_label"].astype(str) == "TP").astype(int)
    sub = sub.sort_values(["Chr", "Pos"]).reset_index(drop=True)
    out = OUT_DIR / f"r1_hcc1395_filter_snvs_{mode}.tsv"
    sub.to_csv(out, sep="\t", index=False)
    tp = int(sub["truth_label_int"].sum())
    fp = len(sub) - tp
    print(f"  {mode}: {len(sub)} SNVs -> {out}  (TP={tp}, FP={fp})")
