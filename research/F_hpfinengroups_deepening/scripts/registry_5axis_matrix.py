#!/usr/bin/env python3
"""Generate 5-axis combinatorial matrix for Research Chain Registry.

Axes:
  mode    : {TO, paired}
  zone    : {Z1, Z2, Z3, Z4, Z5}
  LOH     : {LOH, NonLOH}
  AF band : {extreme (<0.1 OR >0.9), half (0.4-0.6), intermediate (others)}
  CN tier : {CN1 (CovM<1.3), CN2 (1.3<=CovM<2.3), CN3+ (CovM>=2.3)}

Total: 2 × 5 × 2 × 3 × 3 = 180 cells (not 450; 5×5×2×3×3 would need paired+TO+paired×TO combos, but we treat mode as a single dim)

Per cell output:
  n_rows   : total regions
  n_TP     : TP count
  n_FP     : FP count
  TP_rate  : TP / (TP+FP)
  AF_median, NGroups_median
"""
from pathlib import Path
import pandas as pd

MASTER = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
)
OUT = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/research_landscape/"
    "data/10_Registry_5axis_matrix.tsv"
)
OUT.parent.mkdir(parents=True, exist_ok=True)

COLS = [
    "sample", "mode", "Chr", "Pos", "HPFineNGroups", "caller_af",
    "NumReads", "Coverage_Multiple",
    "to_loh_bed_hit", "tool_potential_loh", "core_loh_like",
    "truth_label", "Potential_LOH",
]

print("Loading master dataset...")
parts = []
for chunk in pd.read_csv(MASTER, sep="\t", compression="gzip", usecols=COLS, chunksize=200_000):
    parts.append(chunk)
df = pd.concat(parts, ignore_index=True)
print(f"Total rows: {len(df)}")

# Axis 1: mode (already column)
# Axis 2: zone (Z1-Z5)
#   Definitions per Zone-Aware Framework:
#     Z1: AF>=0.3 AND HPFineNGroups>=3 AND Coverage_Multiple 1.3-2.3 (clean somatic main)
#     Z2: AF 0.05-0.3 AND HPFineNGroups>=2 (low-AF subclone)
#     Z3: LOH AND AF extreme (<0.1 OR >0.9) AND HPFineNGroups<=1 (amplicon suspect)
#     Z4: HPFineNGroups<=1 AND AF>=0.3 (homogeneous)
#     Z5: others
def assign_zone(row):
    af = row["caller_af"]
    ng = row["HPFineNGroups"]
    cv = row["Coverage_Multiple"]
    loh = row["to_loh_bed_hit"] if row["mode"] == "to" else (row["tool_potential_loh"] or row["core_loh_like"])
    af_ext = (af < 0.1) or (af > 0.9)
    # Z3 check first
    if loh and af_ext and ng <= 1:
        return "Z3"
    if af >= 0.3 and ng >= 3 and 1.3 <= cv < 2.3:
        return "Z1"
    if 0.05 <= af < 0.3 and ng >= 2:
        return "Z2"
    if ng <= 1 and af >= 0.3:
        return "Z4"
    return "Z5"

print("Assigning zones...")
df["zone"] = df.apply(assign_zone, axis=1)

# Axis 3: LOH status
def loh_status(row):
    if row["mode"] == "to":
        return "LOH" if row["to_loh_bed_hit"] else "NonLOH"
    return "LOH" if (row["tool_potential_loh"] or row["core_loh_like"]) else "NonLOH"

df["loh_status"] = df.apply(loh_status, axis=1)

# Axis 4: AF band
def af_band(af):
    if af < 0.1 or af > 0.9:
        return "extreme"
    if 0.4 <= af <= 0.6:
        return "half"
    return "intermediate"

df["af_band"] = df["caller_af"].apply(af_band)

# Axis 5: CN tier
def cn_tier(cv):
    if pd.isna(cv):
        return "unknown"
    if cv < 1.3:
        return "CN1"
    if cv < 2.3:
        return "CN2"
    return "CN3+"

df["cn_tier"] = df["Coverage_Multiple"].apply(cn_tier)

# Group stats
def group_stats(g):
    tp = int((g["truth_label"].astype(str) == "TP").sum())
    fp = int((g["truth_label"].astype(str) == "FP").sum())
    n = tp + fp
    return pd.Series({
        "n_rows": len(g),
        "n_TP": tp,
        "n_FP": fp,
        "TP_rate": (tp / n) if n > 0 else 0.0,
        "AF_median": g["caller_af"].median(),
        "NG_median": g["HPFineNGroups"].median(),
        "CovM_median": g["Coverage_Multiple"].median(),
    })

print("Aggregating 5-axis cells...")
grp = df.groupby(["mode", "zone", "loh_status", "af_band", "cn_tier"], observed=True)
matrix = grp.apply(group_stats).reset_index()

# Enumerate full grid so empty cells show up
full_grid = pd.MultiIndex.from_product([
    ["to", "paired"],
    ["Z1", "Z2", "Z3", "Z4", "Z5"],
    ["LOH", "NonLOH"],
    ["extreme", "half", "intermediate"],
    ["CN1", "CN2", "CN3+", "unknown"],
], names=["mode", "zone", "loh_status", "af_band", "cn_tier"])

full = pd.DataFrame(index=full_grid).reset_index()
matrix = full.merge(matrix, on=["mode", "zone", "loh_status", "af_band", "cn_tier"], how="left")
matrix.fillna({"n_rows": 0, "n_TP": 0, "n_FP": 0, "TP_rate": 0.0}, inplace=True)

# Research status classification
def research_status(row):
    n = row["n_rows"]
    if n == 0:
        return "no_data"
    tp_rate = row["TP_rate"]
    zone = row["zone"]
    loh = row["loh_status"]
    af = row["af_band"]
    # Major researched cells from F pilot / B.2 / Z3
    if zone == "Z3" and loh == "LOH" and af == "extreme":
        return "concluded_Z3_pilot"
    if zone == "Z1" and loh == "NonLOH":
        return "validated_main_signal"
    if loh == "LOH" and af == "extreme":
        return "validated_B2_batch1"
    if zone == "Z2":
        return "pilot_subclone"
    if n < 20:
        return "underpowered"
    return "unexplored"

matrix["research_status"] = matrix.apply(research_status, axis=1)

# Effect direction (heuristic)
def effect_direction(row):
    if row["n_rows"] < 20:
        return "NA"
    tp = row["TP_rate"]
    if tp >= 0.8:
        return "POS_high_TP"
    if tp >= 0.6:
        return "POS_mid"
    if tp >= 0.4:
        return "neutral"
    return "NEG_low_TP"

matrix["effect_direction"] = matrix.apply(effect_direction, axis=1)

matrix.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
print(f"[+] Matrix written: {OUT}")
print(f"[+] Total cells: {len(matrix)}")
print(f"[+] Non-empty cells: {(matrix['n_rows'] > 0).sum()}")
print("\n[+] Research status distribution:")
print(matrix["research_status"].value_counts())
print("\n[+] Effect direction distribution:")
print(matrix["effect_direction"].value_counts())

# High-priority cells summary (researched / validated / POS_high_TP)
print("\n[+] Top 20 high-density cells (by n_rows):")
top = matrix.nlargest(20, "n_rows")
print(top[["mode", "zone", "loh_status", "af_band", "cn_tier",
          "n_rows", "TP_rate", "research_status", "effect_direction"]].to_string(index=False))
