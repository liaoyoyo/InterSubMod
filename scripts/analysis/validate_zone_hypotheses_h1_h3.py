#!/usr/bin/env python3
"""
Zone-Aware Confidence Framework — H1/H3 Hypothesis Validation

H1: Zone Z1 (LOH Subclonal Active) has significantly higher TP rate than global
    Condition: LOH=True AND caller_af in [0.1,0.4]∪[0.6,0.9] AND HPFineNGroups >= 3

H3: Zone Z2 (High Somatic Heterogeneity, NGroups>=4 + NR>=80) 89.1% TP rate
    holds in TO mode

Also computes full Zone Z1-Z5 assignment statistics.

Data: all_region_rows.tsv.gz (748,391 rows, 7 samples, post-HP-fix)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# ── Config ─────────────────────────────────────────────────────────
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/" \
         "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/zone_aware_validation"
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLE_ORDER = ["H2009", "H1437", "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "COLO829"]
SAMPLE_COLORS = {
    "H2009": "#1f77b4", "H1437": "#ff7f0e", "HCC1395": "#2ca02c",
    "HCC1395_DORADO": "#d62728", "HCC1937": "#9467bd",
    "HCC1954": "#8c564b", "COLO829": "#e377c2"
}

# ── Load Data ──────────────────────────────────────────────────────
print("Loading master dataset...")
cols_needed = [
    "sample", "mode", "truth_label",
    "core_loh_like", "to_loh_bed_hit",
    "caller_af", "HPFineNGroups", "NumReads",
    "Coverage_Multiple", "effective_hp_reads"
]
df = pd.read_csv(MASTER, sep='\t', usecols=cols_needed, low_memory=False)
print(f"Loaded {len(df):,} rows")

# Type conversions
df["caller_af"] = pd.to_numeric(df["caller_af"], errors="coerce")
df["HPFineNGroups"] = pd.to_numeric(df["HPFineNGroups"], errors="coerce")
df["NumReads"] = pd.to_numeric(df["NumReads"], errors="coerce")
df["Coverage_Multiple"] = pd.to_numeric(df["Coverage_Multiple"], errors="coerce")
df["effective_hp_reads"] = pd.to_numeric(df["effective_hp_reads"], errors="coerce")

# Boolean LOH columns
for col in ["core_loh_like", "to_loh_bed_hit"]:
    df[col] = df[col].map({"True": True, "False": False, True: True, False: False}).fillna(False)

df["is_tp"] = df["truth_label"] == "TP"

# Choose LOH column by mode: paired uses core_loh_like, TO uses to_loh_bed_hit
df["loh_flag"] = np.where(df["mode"] == "paired", df["core_loh_like"], df["to_loh_bed_hit"])

# ── Zone Assignment ────────────────────────────────────────────────
print("Assigning zones...")

def is_intermediate_af(af):
    """Intermediate AF: [0.1, 0.4] union [0.6, 0.9]"""
    if pd.isna(af):
        return False
    return (0.1 <= af <= 0.4) or (0.6 <= af <= 0.9)

def is_extreme_af(af):
    """Extreme AF: [0, 0.1) union (0.9, 1.0]"""
    if pd.isna(af):
        return False
    return (af < 0.1) or (af > 0.9)

df["af_intermediate"] = df["caller_af"].apply(is_intermediate_af)
df["af_extreme"] = df["caller_af"].apply(is_extreme_af)

# Zone Z1: LOH Subclonal Active
#   LOH=True AND caller_af intermediate AND HPFineNGroups >= 3
df["zone_z1"] = (df["loh_flag"]) & (df["af_intermediate"]) & (df["HPFineNGroups"] >= 3)

# Zone Z2: High Somatic Heterogeneity
#   HPFineNGroups >= 4 AND NumReads >= 80 (regardless of LOH)
df["zone_z2"] = (df["HPFineNGroups"] >= 4) & (df["NumReads"] >= 80)

# Zone Z3: Complete LOH
#   LOH=True AND extreme AF AND HPFineNGroups <= 1
df["zone_z3"] = (df["loh_flag"]) & (df["af_extreme"]) & (df["HPFineNGroups"] <= 1)

# Zone Z5: CN Gain with Low Diversity
#   Coverage_Multiple > 1.5 AND HPFineNGroups <= 2 AND NOT LOH
df["zone_z5"] = (df["Coverage_Multiple"] > 1.5) & (df["HPFineNGroups"] <= 2) & (~df["loh_flag"])

# Zone Z4: Normal Diploid (residual — not in any other zone, CovM near 1)
#   Approximation: NOT in Z1/Z2/Z3/Z5 AND Coverage_Multiple in [0.7, 1.3]
df["zone_z4"] = (~df["zone_z1"]) & (~df["zone_z2"]) & (~df["zone_z3"]) & \
                (~df["zone_z5"]) & (df["Coverage_Multiple"].between(0.7, 1.3))

# Priority-based zone label (Z1 > Z2 > Z3 > Z5 > Z4 > Unassigned)
def assign_zone(row):
    if row["zone_z1"]: return "Z1"
    if row["zone_z2"]: return "Z2"
    if row["zone_z3"]: return "Z3"
    if row["zone_z5"]: return "Z5"
    if row["zone_z4"]: return "Z4"
    return "Unassigned"

df["zone"] = df.apply(assign_zone, axis=1)

# ══════════════════════════════════════════════════════════════════
# H1 VALIDATION: Zone Z1 TP rate vs global, per sample, both modes
# ══════════════════════════════════════════════════════════════════
print("\n=== H1: Zone Z1 TP Rate Validation ===")

h1_results = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for sample in SAMPLE_ORDER:
        ds = dm[dm["sample"] == sample]
        if len(ds) == 0:
            continue

        global_tp_rate = ds["is_tp"].mean()
        global_n = len(ds)

        z1 = ds[ds["zone_z1"]]
        z1_n = len(z1)
        z1_tp_rate = z1["is_tp"].mean() if z1_n > 0 else np.nan
        z1_tp = z1["is_tp"].sum() if z1_n > 0 else 0
        z1_fp = z1_n - z1_tp

        # Fisher's exact test: Z1 vs non-Z1
        non_z1 = ds[~ds["zone_z1"]]
        non_z1_tp = non_z1["is_tp"].sum()
        non_z1_fp = len(non_z1) - non_z1_tp

        if z1_n > 0 and len(non_z1) > 0:
            table = [[z1_tp, z1_fp], [non_z1_tp, non_z1_fp]]
            odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
        else:
            odds_ratio, p_value = np.nan, np.nan

        h1_results.append({
            "mode": mode, "sample": sample,
            "global_n": global_n, "global_tp_rate": global_tp_rate,
            "z1_n": z1_n, "z1_tp_rate": z1_tp_rate,
            "z1_tp": z1_tp, "z1_fp": z1_fp,
            "delta_tp_rate": z1_tp_rate - global_tp_rate if not np.isnan(z1_tp_rate) else np.nan,
            "odds_ratio": odds_ratio, "p_value": p_value
        })

        sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
        print(f"  {mode:6s} {sample:16s}: Z1 TP={z1_tp_rate:.3f} vs Global={global_tp_rate:.3f} "
              f"(delta={z1_tp_rate - global_tp_rate:+.3f}) N={z1_n:5d} OR={odds_ratio:.2f} {sig}")

h1_df = pd.DataFrame(h1_results)
h1_df.to_csv(os.path.join(OUT_DIR, "H1_zone_z1_tp_rate.tsv"), sep='\t', index=False)

# ── H1 Figure ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)

for idx, mode in enumerate(["paired", "to"]):
    ax = axes[idx]
    mode_data = h1_df[h1_df["mode"] == mode]

    x = np.arange(len(mode_data))
    width = 0.35

    bars_global = ax.bar(x - width/2, mode_data["global_tp_rate"], width,
                         color='lightgray', edgecolor='gray', label='Global TP Rate')
    bars_z1 = ax.bar(x + width/2, mode_data["z1_tp_rate"], width,
                     color=[SAMPLE_COLORS.get(s, '#333') for s in mode_data["sample"]],
                     edgecolor='black', label='Z1 TP Rate')

    # Add N labels and significance
    for i, (_, row) in enumerate(mode_data.iterrows()):
        sig = "***" if row["p_value"] < 0.001 else "**" if row["p_value"] < 0.01 else \
              "*" if row["p_value"] < 0.05 else "ns"
        delta = row["delta_tp_rate"]
        if not np.isnan(delta):
            ax.annotate(f'{delta:+.3f}\n{sig}\nN={int(row["z1_n"])}',
                       xy=(i + width/2, row["z1_tp_rate"]),
                       xytext=(0, 5), textcoords='offset points',
                       ha='center', va='bottom', fontsize=7,
                       fontweight='bold' if row["p_value"] < 0.05 else 'normal')

    ax.set_xticks(x)
    ax.set_xticklabels(mode_data["sample"], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel("TP Rate" if idx == 0 else "")
    ax.set_title(f"H1: Zone Z1 TP Rate — {mode.upper()} mode", fontweight='bold')
    ax.set_ylim(0.5, 1.05)
    ax.axhline(y=mode_data["global_tp_rate"].mean(), color='red', linestyle='--', alpha=0.5,
               label=f'Mean Global ({mode_data["global_tp_rate"].mean():.3f})')
    ax.legend(fontsize=8, loc='lower left')
    ax.grid(axis='y', alpha=0.3)

plt.suptitle("H1 Validation: Zone Z1 (LOH Subclonal Active) TP Rate vs Global",
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "H1_zone_z1_tp_rate_validation.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  → Figure: {FIG_DIR}/H1_zone_z1_tp_rate_validation.png")

# ══════════════════════════════════════════════════════════════════
# H3 VALIDATION: NGroups>=4 + NR>=80 TP rate in TO mode
# ══════════════════════════════════════════════════════════════════
print("\n=== H3: Zone Z2 TP Rate in TO Mode ===")

h3_results = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for sample in SAMPLE_ORDER:
        ds = dm[dm["sample"] == sample]
        if len(ds) == 0:
            continue

        global_tp_rate = ds["is_tp"].mean()
        global_n = len(ds)

        z2 = ds[ds["zone_z2"]]
        z2_n = len(z2)
        z2_tp_rate = z2["is_tp"].mean() if z2_n > 0 else np.nan
        z2_tp = z2["is_tp"].sum() if z2_n > 0 else 0
        z2_fp = z2_n - z2_tp

        # Fisher's exact
        non_z2 = ds[~ds["zone_z2"]]
        non_z2_tp = non_z2["is_tp"].sum()
        non_z2_fp = len(non_z2) - non_z2_tp

        if z2_n > 0 and len(non_z2) > 0:
            table = [[z2_tp, z2_fp], [non_z2_tp, non_z2_fp]]
            odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
        else:
            odds_ratio, p_value = np.nan, np.nan

        h3_results.append({
            "mode": mode, "sample": sample,
            "global_n": global_n, "global_tp_rate": global_tp_rate,
            "z2_n": z2_n, "z2_tp_rate": z2_tp_rate,
            "z2_tp": z2_tp, "z2_fp": z2_fp,
            "delta_tp_rate": z2_tp_rate - global_tp_rate if not np.isnan(z2_tp_rate) else np.nan,
            "odds_ratio": odds_ratio, "p_value": p_value
        })

        sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
        print(f"  {mode:6s} {sample:16s}: Z2 TP={z2_tp_rate:.3f} vs Global={global_tp_rate:.3f} "
              f"(delta={z2_tp_rate - global_tp_rate:+.3f}) N={z2_n:5d} OR={odds_ratio:.2f} {sig}")

h3_df = pd.DataFrame(h3_results)
h3_df.to_csv(os.path.join(OUT_DIR, "H3_zone_z2_tp_rate.tsv"), sep='\t', index=False)

# ── H3 Figure ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)

for idx, mode in enumerate(["paired", "to"]):
    ax = axes[idx]
    mode_data = h3_df[h3_df["mode"] == mode]

    x = np.arange(len(mode_data))
    width = 0.35

    bars_global = ax.bar(x - width/2, mode_data["global_tp_rate"], width,
                         color='lightgray', edgecolor='gray', label='Global TP Rate')
    bars_z2 = ax.bar(x + width/2, mode_data["z2_tp_rate"], width,
                     color=[SAMPLE_COLORS.get(s, '#333') for s in mode_data["sample"]],
                     edgecolor='black', label='Z2 TP Rate')

    for i, (_, row) in enumerate(mode_data.iterrows()):
        sig = "***" if row["p_value"] < 0.001 else "**" if row["p_value"] < 0.01 else \
              "*" if row["p_value"] < 0.05 else "ns"
        delta = row["delta_tp_rate"]
        if not np.isnan(delta):
            ax.annotate(f'{delta:+.3f}\n{sig}\nN={int(row["z2_n"])}',
                       xy=(i + width/2, row["z2_tp_rate"]),
                       xytext=(0, 5), textcoords='offset points',
                       ha='center', va='bottom', fontsize=7,
                       fontweight='bold' if row["p_value"] < 0.05 else 'normal')

    ax.set_xticks(x)
    ax.set_xticklabels(mode_data["sample"], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel("TP Rate" if idx == 0 else "")
    ax.set_title(f"H3: Zone Z2 TP Rate — {mode.upper()} mode", fontweight='bold')
    ax.set_ylim(0.5, 1.05)
    ax.axhline(y=0.891, color='green', linestyle='--', alpha=0.5,
               label='89.1% (prior claim)')
    ax.legend(fontsize=8, loc='lower left')
    ax.grid(axis='y', alpha=0.3)

plt.suptitle("H3 Validation: Zone Z2 (HPFineNGroups≥4 + NR≥80) TP Rate",
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "H3_zone_z2_tp_rate_validation.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  → Figure: {FIG_DIR}/H3_zone_z2_tp_rate_validation.png")

# ══════════════════════════════════════════════════════════════════
# FULL ZONE ASSIGNMENT STATISTICS
# ══════════════════════════════════════════════════════════════════
print("\n=== Full Zone Assignment Statistics ===")

zone_stats = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for zone_label in ["Z1", "Z2", "Z3", "Z4", "Z5", "Unassigned"]:
        zd = dm[dm["zone"] == zone_label]
        n = len(zd)
        tp = zd["is_tp"].sum()
        fp = n - tp
        tp_rate = tp / n if n > 0 else np.nan
        pct = n / len(dm) * 100

        zone_stats.append({
            "mode": mode, "zone": zone_label,
            "n": n, "tp": tp, "fp": fp,
            "tp_rate": tp_rate, "pct_of_mode": pct
        })
        print(f"  {mode:6s} {zone_label:12s}: N={n:7d} ({pct:5.1f}%) TP={tp:6d} FP={fp:5d} TP_rate={tp_rate:.4f}")

zone_df = pd.DataFrame(zone_stats)
zone_df.to_csv(os.path.join(OUT_DIR, "zone_assignment_statistics.tsv"), sep='\t', index=False)

# Per-sample zone breakdown
print("\n=== Per-Sample Zone Distribution ===")
sample_zone_stats = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for sample in SAMPLE_ORDER:
        ds = dm[dm["sample"] == sample]
        if len(ds) == 0:
            continue
        for zone_label in ["Z1", "Z2", "Z3", "Z4", "Z5", "Unassigned"]:
            zd = ds[ds["zone"] == zone_label]
            n = len(zd)
            tp = zd["is_tp"].sum()
            fp = n - tp
            tp_rate = tp / n if n > 0 else np.nan
            pct = n / len(ds) * 100
            sample_zone_stats.append({
                "mode": mode, "sample": sample, "zone": zone_label,
                "n": n, "tp": tp, "fp": fp,
                "tp_rate": tp_rate, "pct_of_sample": pct
            })

szdf = pd.DataFrame(sample_zone_stats)
szdf.to_csv(os.path.join(OUT_DIR, "zone_per_sample_statistics.tsv"), sep='\t', index=False)

# ── Zone Overview Figure ───────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(18, 14))

# Panel 1: Zone size distribution (paired)
ax = axes[0, 0]
paired_zones = zone_df[zone_df["mode"] == "paired"]
colors_zones = ['#e74c3c', '#3498db', '#f39c12', '#2ecc71', '#9b59b6', '#95a5a6']
ax.bar(paired_zones["zone"], paired_zones["n"], color=colors_zones, edgecolor='black')
for i, (_, row) in enumerate(paired_zones.iterrows()):
    ax.annotate(f'{int(row["n"]):,}\n({row["pct_of_mode"]:.1f}%)',
               xy=(i, row["n"]), xytext=(0, 5), textcoords='offset points',
               ha='center', va='bottom', fontsize=9)
ax.set_title("Zone Size Distribution — Paired", fontweight='bold')
ax.set_ylabel("Number of Regions")
ax.grid(axis='y', alpha=0.3)

# Panel 2: Zone size distribution (TO)
ax = axes[0, 1]
to_zones = zone_df[zone_df["mode"] == "to"]
ax.bar(to_zones["zone"], to_zones["n"], color=colors_zones, edgecolor='black')
for i, (_, row) in enumerate(to_zones.iterrows()):
    ax.annotate(f'{int(row["n"]):,}\n({row["pct_of_mode"]:.1f}%)',
               xy=(i, row["n"]), xytext=(0, 5), textcoords='offset points',
               ha='center', va='bottom', fontsize=9)
ax.set_title("Zone Size Distribution — TO", fontweight='bold')
ax.set_ylabel("Number of Regions")
ax.grid(axis='y', alpha=0.3)

# Panel 3: TP rate by zone (both modes)
ax = axes[1, 0]
x = np.arange(len(paired_zones))
width = 0.35
ax.bar(x - width/2, paired_zones["tp_rate"], width, color='steelblue',
       edgecolor='black', label='Paired')
ax.bar(x + width/2, to_zones["tp_rate"], width, color='coral',
       edgecolor='black', label='TO')
for i, (_, row) in enumerate(paired_zones.iterrows()):
    if not np.isnan(row["tp_rate"]):
        ax.annotate(f'{row["tp_rate"]:.3f}', xy=(i - width/2, row["tp_rate"]),
                   xytext=(0, 3), textcoords='offset points', ha='center', fontsize=8, color='steelblue')
for i, (_, row) in enumerate(to_zones.iterrows()):
    if not np.isnan(row["tp_rate"]):
        ax.annotate(f'{row["tp_rate"]:.3f}', xy=(i + width/2, row["tp_rate"]),
                   xytext=(0, 3), textcoords='offset points', ha='center', fontsize=8, color='coral')
ax.set_xticks(x)
ax.set_xticklabels(paired_zones["zone"])
ax.set_ylabel("TP Rate")
ax.set_title("TP Rate by Zone — Paired vs TO", fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0.4, 1.05)

# Panel 4: Per-sample Z1 TP rate heatmap
ax = axes[1, 1]
# Build pivot table: sample x mode for Z1 TP rate
z1_pivot = szdf[szdf["zone"] == "Z1"].pivot_table(
    index="sample", columns="mode", values="tp_rate")
z1_pivot = z1_pivot.reindex(SAMPLE_ORDER)

# Also add Z2 for comparison
z2_pivot = szdf[szdf["zone"] == "Z2"].pivot_table(
    index="sample", columns="mode", values="tp_rate")
z2_pivot = z2_pivot.reindex(SAMPLE_ORDER)

# Combined heatmap data
hm_data = pd.DataFrame({
    "Z1 Paired": z1_pivot.get("paired", np.nan),
    "Z1 TO": z1_pivot.get("to", np.nan),
    "Z2 Paired": z2_pivot.get("paired", np.nan),
    "Z2 TO": z2_pivot.get("to", np.nan)
}, index=SAMPLE_ORDER)

im = ax.imshow(hm_data.values, cmap='RdYlGn', aspect='auto', vmin=0.6, vmax=1.0)
ax.set_xticks(range(4))
ax.set_xticklabels(hm_data.columns, rotation=30, ha='right', fontsize=9)
ax.set_yticks(range(len(SAMPLE_ORDER)))
ax.set_yticklabels(SAMPLE_ORDER, fontsize=9)
for i in range(len(SAMPLE_ORDER)):
    for j in range(4):
        val = hm_data.values[i, j]
        if not np.isnan(val):
            ax.text(j, i, f'{val:.3f}', ha='center', va='center', fontsize=9,
                   fontweight='bold', color='white' if val < 0.75 else 'black')
ax.set_title("Z1/Z2 TP Rate by Sample × Mode", fontweight='bold')
plt.colorbar(im, ax=ax, label='TP Rate', shrink=0.8)

plt.suptitle("Zone-Aware Framework — Full Zone Assignment Overview",
             fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "zone_overview_statistics.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  → Figure: {FIG_DIR}/zone_overview_statistics.png")

# ══════════════════════════════════════════════════════════════════
# Z1 vs Z3 CONTRAST (Subclonal Active vs Complete LOH)
# ══════════════════════════════════════════════════════════════════
print("\n=== Z1 vs Z3 Contrast (Subclonal vs Complete LOH) ===")

z1z3_results = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for sample in SAMPLE_ORDER:
        ds = dm[dm["sample"] == sample]
        if len(ds) == 0:
            continue
        z1 = ds[ds["zone"] == "Z1"]
        z3 = ds[ds["zone"] == "Z3"]
        z1_tp_rate = z1["is_tp"].mean() if len(z1) > 0 else np.nan
        z3_tp_rate = z3["is_tp"].mean() if len(z3) > 0 else np.nan
        z1z3_results.append({
            "mode": mode, "sample": sample,
            "z1_n": len(z1), "z1_tp_rate": z1_tp_rate,
            "z3_n": len(z3), "z3_tp_rate": z3_tp_rate,
            "z1_minus_z3": z1_tp_rate - z3_tp_rate if not (np.isnan(z1_tp_rate) or np.isnan(z3_tp_rate)) else np.nan
        })

z1z3_df = pd.DataFrame(z1z3_results)
z1z3_df.to_csv(os.path.join(OUT_DIR, "Z1_vs_Z3_contrast.tsv"), sep='\t', index=False)

for _, r in z1z3_df.iterrows():
    delta = r["z1_minus_z3"]
    dstr = f"{delta:+.3f}" if not np.isnan(delta) else "N/A"
    print(f"  {r['mode']:6s} {r['sample']:16s}: Z1={r['z1_tp_rate']:.3f}(N={r['z1_n']}) "
          f"Z3={r['z3_tp_rate']:.3f}(N={r['z3_n']}) delta={dstr}")

# ══════════════════════════════════════════════════════════════════
# H1/H3 SUMMARY STATISTICS
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

# H1 summary
h1_paired = h1_df[h1_df["mode"] == "paired"]
h1_to = h1_df[h1_df["mode"] == "to"]

print("\nH1 (Zone Z1 TP Rate > Global):")
for mode, mdf in [("Paired", h1_paired), ("TO", h1_to)]:
    n_sig = (mdf["p_value"] < 0.05).sum()
    n_positive = (mdf["delta_tp_rate"] > 0).sum()
    mean_delta = mdf["delta_tp_rate"].mean()
    mean_z1_tp = mdf["z1_tp_rate"].mean()
    total_z1_n = mdf["z1_n"].sum()
    print(f"  {mode}: {n_positive}/{len(mdf)} positive, {n_sig}/{len(mdf)} significant (p<0.05)")
    print(f"    Mean Z1 TP rate: {mean_z1_tp:.4f}, Mean delta: {mean_delta:+.4f}")
    print(f"    Total Z1 regions: {total_z1_n:,}")

# H3 summary
h3_paired = h3_df[h3_df["mode"] == "paired"]
h3_to = h3_df[h3_df["mode"] == "to"]

print("\nH3 (Zone Z2 TP Rate, 89.1% claim):")
for mode, mdf in [("Paired", h3_paired), ("TO", h3_to)]:
    n_above_891 = (mdf["z2_tp_rate"] >= 0.891).sum()
    n_sig = (mdf["p_value"] < 0.05).sum()
    mean_z2_tp = mdf["z2_tp_rate"].mean()
    total_z2_n = mdf["z2_n"].sum()
    print(f"  {mode}: {n_above_891}/{len(mdf)} above 89.1%, {n_sig}/{len(mdf)} significant")
    print(f"    Mean Z2 TP rate: {mean_z2_tp:.4f}")
    print(f"    Total Z2 regions: {total_z2_n:,}")

# Zone distribution summary
print("\nZone Distribution Summary:")
for mode in ["paired", "to"]:
    print(f"\n  {mode.upper()}:")
    mzd = zone_df[zone_df["mode"] == mode]
    for _, row in mzd.iterrows():
        print(f"    {row['zone']:12s}: {row['n']:7,d} ({row['pct_of_mode']:5.1f}%) TP_rate={row['tp_rate']:.4f}")

print(f"\nAll outputs saved to: {OUT_DIR}/")
print("Done.")
