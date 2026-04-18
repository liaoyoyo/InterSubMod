#!/usr/bin/env python3
"""
LOH Expanded Observation Analysis (2026-03-31)

Q1: LOH vs non-LOH TP/FP baseline (enrichment of non-LOH)
Q2: Per-Tier TP/FP LOH-like rates (expanded per-sample × per-tier × per-mode)
Q3: TO-specific dimensions × LOH (to_loh_bed_hit, SuggestFilter, phase_block_status)
Q4: VerificationClass × LOH × Tier 3D cross
Q5: Methylation features in LOH-like subgroups (TP vs FP, Tier A vs A+)

Input: all_region_rows.tsv.gz (post-HP-fix master dataset)
"""

import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")

# ============================================================
# Configuration
# ============================================================
ROUND1_WS = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260327_loh_round1_cross_sample_audit"
)
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260331_loh_expanded_observation"
)
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

GENERATED_AT = datetime.now().strftime("%Y-%m-%d %H:%M:%S CST")
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_LABELS = {
    "HCC1395": "HCC1395_5kHz", "HCC1395_DORADO": "HCC1395_DORADO",
    "COLO829": "COLO829", "H1437": "H1437", "H2009": "H2009",
    "HCC1937": "HCC1937", "HCC1954": "HCC1954",
}
MODE_COLORS = {"paired": "#2166ac", "to": "#b2182b"}
TRUTH_COLORS = {"TP": "#4daf4a", "FP": "#e41a1c"}
LOH_COLORS = {True: "#d62728", False: "#1f77b4"}
TIER_ORDER = ["C0", "C", "B", "A", "A+"]
TIER_COLORS = {"C0": "#d9d9d9", "C": "#fc8d59", "B": "#fee090", "A": "#4575b4", "A+": "#313695"}
VC_ORDER = ["Noise", "Weak", "Strong", "Subclone"]


def assign_tier(n):
    if pd.isna(n): return "C0"
    n = int(n)
    if n == 0: return "C0"
    elif n < 10: return "C"
    elif n < 30: return "B"
    elif n < 50: return "A"
    else: return "A+"


def write_tsv(df, path, desc=""):
    df.to_csv(path, sep="\t", index=False)
    print(f"  {path.name} ({len(df)} rows) — {desc}")


def enrichment(df_sub, loh_col="core_loh_like", truth_col="truth_label"):
    """FP/TP enrichment: (FP_loh/FP_total) / (TP_loh/TP_total).
    >1 = FP enriched, <1 = TP enriched."""
    tp = df_sub[df_sub[truth_col] == "TP"]
    fp = df_sub[df_sub[truth_col] == "FP"]
    tp_loh = tp[loh_col].sum()
    fp_loh = fp[loh_col].sum()
    tp_total = len(tp)
    fp_total = len(fp)
    if tp_total == 0 or fp_total == 0 or tp_loh == 0:
        return np.nan
    return (fp_loh / fp_total) / (tp_loh / tp_total)


# ============================================================
# Data Loading
# ============================================================
print(f"[0/5] Loading master dataset ...")
df = pd.read_csv(ROUND1_WS / "all_region_rows.tsv.gz", sep="\t", low_memory=False)
df["tier"] = df["effective_hp_reads"].apply(assign_tier)
df["tier"] = pd.Categorical(df["tier"], categories=TIER_ORDER, ordered=True)
# Ensure boolean
df["core_loh_like"] = df["core_loh_like"].astype(bool)
df["HPMergedSig"] = df["HPMergedSig"].map({"True": True, "False": False, True: True, False: False}).fillna(False).astype(bool)
print(f"  Loaded {len(df):,} rows (TO: {(df['mode']=='to').sum():,}, paired: {(df['mode']=='paired').sum():,})")
print(f"  TP: {(df['truth_label']=='TP').sum():,}, FP: {(df['truth_label']=='FP').sum():,}")
print(f"  LOH-like: {df['core_loh_like'].sum():,}, non-LOH: {(~df['core_loh_like']).sum():,}")


# ============================================================
# Q1: LOH vs non-LOH TP/FP baseline
# ============================================================
print(f"\n[1/5] Q1 — LOH vs non-LOH TP/FP baseline ...")

# --- Q1a: Global 2x2 table per mode ---
rows_q1a = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for loh_status in [True, False]:
        sub = dm[dm["core_loh_like"] == loh_status]
        n_total = len(sub)
        n_tp = (sub["truth_label"] == "TP").sum()
        n_fp = (sub["truth_label"] == "FP").sum()
        fp_rate = n_fp / n_total if n_total > 0 else np.nan
        rows_q1a.append({
            "mode": mode,
            "loh_like": loh_status,
            "n_total": n_total,
            "n_tp": n_tp,
            "n_fp": n_fp,
            "fp_rate": round(fp_rate, 5),
            "tp_rate": round(1 - fp_rate, 5) if not np.isnan(fp_rate) else np.nan,
        })
df_q1a = pd.DataFrame(rows_q1a)
write_tsv(df_q1a, OUT_DIR / "q1a_loh_vs_nonloh_fp_rate.tsv", "LOH vs non-LOH FP rate by mode")

# --- Q1b: Per-sample LOH vs non-LOH FP rate ---
rows_q1b = []
for sample in SAMPLES:
    for mode in ["paired", "to"]:
        dm = df[(df["sample"] == sample) & (df["mode"] == mode)]
        if len(dm) == 0:
            continue
        for loh_status in [True, False]:
            sub = dm[dm["core_loh_like"] == loh_status]
            n_total = len(sub)
            n_fp = (sub["truth_label"] == "FP").sum()
            n_tp = (sub["truth_label"] == "TP").sum()
            fp_rate = n_fp / n_total if n_total > 0 else np.nan
            rows_q1b.append({
                "sample": sample, "mode": mode, "loh_like": loh_status,
                "n_total": n_total, "n_tp": n_tp, "n_fp": n_fp,
                "fp_rate": round(fp_rate, 5),
            })
df_q1b = pd.DataFrame(rows_q1b)
write_tsv(df_q1b, OUT_DIR / "q1b_per_sample_loh_vs_nonloh_fp_rate.tsv", "Per-sample LOH vs non-LOH FP rate")

# --- Q1c: Per-tier LOH vs non-LOH enrichment ---
rows_q1c = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for tier in TIER_ORDER:
        dt = dm[dm["tier"] == tier]
        if len(dt) < 10:
            continue
        for loh_status in [True, False]:
            sub = dt[dt["core_loh_like"] == loh_status]
            n_total = len(sub)
            n_fp = (sub["truth_label"] == "FP").sum()
            n_tp = (sub["truth_label"] == "TP").sum()
            fp_rate = n_fp / n_total if n_total > 0 else np.nan
            rows_q1c.append({
                "mode": mode, "tier": tier, "loh_like": loh_status,
                "n_total": n_total, "n_tp": n_tp, "n_fp": n_fp,
                "fp_rate": round(fp_rate, 5),
            })
df_q1c = pd.DataFrame(rows_q1c)
write_tsv(df_q1c, OUT_DIR / "q1c_per_tier_loh_vs_nonloh_fp_rate.tsv", "Per-tier LOH vs non-LOH FP rate")

# --- Fig01: Grouped bar chart — LOH vs non-LOH FP rate by mode ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_q1a[df_q1a["mode"] == mode]
    loh_row = sub[sub["loh_like"] == True].iloc[0]
    nonloh_row = sub[sub["loh_like"] == False].iloc[0]
    bars = ax.bar(
        [0, 1],
        [loh_row["fp_rate"], nonloh_row["fp_rate"]],
        color=[LOH_COLORS[True], LOH_COLORS[False]],
        width=0.6, edgecolor="black", linewidth=0.5,
    )
    for bar, row in zip(bars, [loh_row, nonloh_row]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f'{row["fp_rate"]:.3f}\n(n={row["n_total"]:,})',
                ha="center", va="bottom", fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["LOH-like", "non-LOH"])
    ax.set_ylabel("FP Rate" if i == 0 else "")
    ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
    ax.axhline(y=sub["n_fp"].sum() / sub["n_total"].sum(), color="gray", linestyle="--", alpha=0.5, label="overall FP rate")
    ax.legend(fontsize=8)
fig.suptitle("Q1: FP Rate — LOH-like vs non-LOH", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig01_loh_vs_nonloh_fp_rate.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig01_loh_vs_nonloh_fp_rate.png")

# --- Fig02: Per-sample grouped bar (both modes) ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_q1b[df_q1b["mode"] == mode].copy()
    samples_in = sub["sample"].unique()
    x = np.arange(len(samples_in))
    w = 0.35
    loh_rates = []
    nonloh_rates = []
    for s in samples_in:
        ss = sub[sub["sample"] == s]
        lr = ss[ss["loh_like"] == True]
        nr = ss[ss["loh_like"] == False]
        loh_rates.append(lr["fp_rate"].values[0] if len(lr) > 0 else 0)
        nonloh_rates.append(nr["fp_rate"].values[0] if len(nr) > 0 else 0)
    ax.bar(x - w/2, loh_rates, w, color=LOH_COLORS[True], label="LOH-like", edgecolor="black", linewidth=0.3)
    ax.bar(x + w/2, nonloh_rates, w, color=LOH_COLORS[False], label="non-LOH", edgecolor="black", linewidth=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels([SAMPLE_LABELS.get(s, s) for s in samples_in], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("FP Rate" if i == 0 else "")
    ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8)
fig.suptitle("Q1: Per-Sample FP Rate — LOH-like vs non-LOH", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig02_per_sample_loh_vs_nonloh_fp_rate.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig02_per_sample_loh_vs_nonloh_fp_rate.png")

# --- Fig03: Per-tier LOH vs non-LOH FP rate ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_q1c[df_q1c["mode"] == mode].copy()
    tiers = [t for t in TIER_ORDER if t in sub["tier"].values]
    x = np.arange(len(tiers))
    w = 0.35
    loh_rates = []
    nonloh_rates = []
    for t in tiers:
        ss = sub[sub["tier"] == t]
        lr = ss[ss["loh_like"] == True]
        nr = ss[ss["loh_like"] == False]
        loh_rates.append(lr["fp_rate"].values[0] if len(lr) > 0 else 0)
        nonloh_rates.append(nr["fp_rate"].values[0] if len(nr) > 0 else 0)
    ax.bar(x - w/2, loh_rates, w, color=LOH_COLORS[True], label="LOH-like", edgecolor="black", linewidth=0.3)
    ax.bar(x + w/2, nonloh_rates, w, color=LOH_COLORS[False], label="non-LOH", edgecolor="black", linewidth=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels(tiers)
    ax.set_ylabel("FP Rate" if i == 0 else "")
    ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8)
fig.suptitle("Q1: Per-Tier FP Rate — LOH-like vs non-LOH", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig03_per_tier_loh_vs_nonloh_fp_rate.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig03_per_tier_loh_vs_nonloh_fp_rate.png")


# ============================================================
# Q2: Per-Tier TP/FP LOH-like rates (expanded)
# ============================================================
print(f"\n[2/5] Q2 — Per-Tier LOH-like rates ...")

# --- Q2a: Full cross: sample × mode × tier × truth × LOH-like rate ---
rows_q2a = []
for sample in SAMPLES:
    for mode in ["paired", "to"]:
        dm = df[(df["sample"] == sample) & (df["mode"] == mode)]
        if len(dm) == 0:
            continue
        for tier in TIER_ORDER:
            dt = dm[dm["tier"] == tier]
            for truth in ["TP", "FP"]:
                sub = dt[dt["truth_label"] == truth]
                n = len(sub)
                n_loh = sub["core_loh_like"].sum()
                loh_frac = n_loh / n if n > 0 else np.nan
                rows_q2a.append({
                    "sample": sample, "mode": mode, "tier": tier, "truth": truth,
                    "n": n, "n_loh": int(n_loh), "loh_frac": round(loh_frac, 4) if not np.isnan(loh_frac) else np.nan,
                })
df_q2a = pd.DataFrame(rows_q2a)
write_tsv(df_q2a, OUT_DIR / "q2a_full_cross_loh_rate.tsv", "Full cross: sample × mode × tier × truth LOH rate")

# --- Q2b: Global mode × tier enrichment ---
rows_q2b = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for tier in TIER_ORDER:
        dt = dm[dm["tier"] == tier]
        if len(dt) < 10:
            continue
        e = enrichment(dt)
        tp_loh_frac = dt[dt["truth_label"] == "TP"]["core_loh_like"].mean()
        fp_loh_frac = dt[dt["truth_label"] == "FP"]["core_loh_like"].mean()
        rows_q2b.append({
            "mode": mode, "tier": tier,
            "n_total": len(dt),
            "tp_loh_frac": round(tp_loh_frac, 4),
            "fp_loh_frac": round(fp_loh_frac, 4),
            "enrichment": round(e, 4) if not np.isnan(e) else np.nan,
            "delta_pp": round((fp_loh_frac - tp_loh_frac) * 100, 2),
        })
df_q2b = pd.DataFrame(rows_q2b)
write_tsv(df_q2b, OUT_DIR / "q2b_mode_tier_enrichment.tsv", "Mode × Tier enrichment")

# --- Fig04: Heatmap of LOH enrichment (sample × tier, paired and TO) ---
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_q2a[(df_q2a["mode"] == mode)].copy()
    # Compute enrichment per sample × tier
    pivot_data = []
    for sample in SAMPLES:
        row = {"sample": SAMPLE_LABELS.get(sample, sample)}
        for tier in TIER_ORDER:
            tp_row = sub[(sub["sample"] == sample) & (sub["tier"] == tier) & (sub["truth"] == "TP")]
            fp_row = sub[(sub["sample"] == sample) & (sub["tier"] == tier) & (sub["truth"] == "FP")]
            if len(tp_row) > 0 and len(fp_row) > 0:
                tp_frac = tp_row.iloc[0]["loh_frac"]
                fp_frac = fp_row.iloc[0]["loh_frac"]
                if tp_frac > 0 and not np.isnan(tp_frac) and not np.isnan(fp_frac):
                    row[tier] = fp_frac / tp_frac
                else:
                    row[tier] = np.nan
            else:
                row[tier] = np.nan
        pivot_data.append(row)
    pivot_df = pd.DataFrame(pivot_data).set_index("sample")
    pivot_df = pivot_df[TIER_ORDER]

    sns.heatmap(pivot_df.astype(float), ax=ax, annot=True, fmt=".2f",
                cmap="RdBu_r", center=1.0, vmin=0.3, vmax=2.5,
                linewidths=0.5, cbar_kws={"label": "FP/TP enrichment"})
    ax.set_title(f"{mode.upper()} — LOH Enrichment by Sample × Tier", fontsize=12, fontweight="bold")
    ax.set_ylabel("")
fig.suptitle("Q2: LOH Enrichment Heatmap (red=FP enriched, blue=TP enriched)", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig04_enrichment_heatmap_sample_tier.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig04_enrichment_heatmap_sample_tier.png")

# --- Fig05: LOH-like rate by tier (TP vs FP, paired vs TO) ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_q2b[df_q2b["mode"] == mode].copy()
    tiers = [t for t in TIER_ORDER if t in sub["tier"].values]
    x = np.arange(len(tiers))
    w = 0.35
    tp_fracs = [sub[sub["tier"] == t]["tp_loh_frac"].values[0] for t in tiers]
    fp_fracs = [sub[sub["tier"] == t]["fp_loh_frac"].values[0] for t in tiers]
    ax.bar(x - w/2, tp_fracs, w, color=TRUTH_COLORS["TP"], label="TP LOH-like rate", edgecolor="black", linewidth=0.3)
    ax.bar(x + w/2, fp_fracs, w, color=TRUTH_COLORS["FP"], label="FP LOH-like rate", edgecolor="black", linewidth=0.3)
    # Add text labels
    for j, t in enumerate(tiers):
        ax.text(x[j] - w/2, tp_fracs[j] + 0.01, f'{tp_fracs[j]:.2f}', ha="center", fontsize=7)
        ax.text(x[j] + w/2, fp_fracs[j] + 0.01, f'{fp_fracs[j]:.2f}', ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(tiers)
    ax.set_ylabel("LOH-like Fraction" if i == 0 else "")
    ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8)
fig.suptitle("Q2: LOH-like Rate by Tier — TP vs FP", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig05_loh_rate_by_tier_tp_vs_fp.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig05_loh_rate_by_tier_tp_vs_fp.png")


# ============================================================
# Q3: TO-specific dimensions × LOH
# ============================================================
print(f"\n[3/5] Q3 — TO-specific dimensions × LOH ...")

df_to = df[df["mode"] == "to"].copy()

# --- Q3a: to_loh_bed_hit × LOH × truth ---
rows_q3a = []
for bed_hit in [True, False]:
    sub_bed = df_to[df_to["to_loh_bed_hit"] == bed_hit] if "to_loh_bed_hit" in df_to.columns else pd.DataFrame()
    if len(sub_bed) == 0:
        # Try string matching
        sub_bed = df_to[df_to["to_loh_bed_hit"].astype(str).str.lower() == str(bed_hit).lower()]
    for truth in ["TP", "FP"]:
        sub = sub_bed[sub_bed["truth_label"] == truth]
        n = len(sub)
        n_loh = sub["core_loh_like"].sum() if n > 0 else 0
        rows_q3a.append({
            "to_loh_bed_hit": bed_hit, "truth": truth,
            "n": n, "n_loh": int(n_loh),
            "loh_frac": round(n_loh / n, 4) if n > 0 else np.nan,
            "fp_rate_in_group": round(n / len(sub_bed), 4) if len(sub_bed) > 0 else np.nan,
        })
df_q3a = pd.DataFrame(rows_q3a)
write_tsv(df_q3a, OUT_DIR / "q3a_to_loh_bed_hit_cross.tsv", "TO loh_bed_hit × truth × LOH-like")

# --- Q3b: phase_block_status × LOH × truth ---
rows_q3b = []
for pbs in df_to["phase_block_status"].unique():
    sub_pbs = df_to[df_to["phase_block_status"] == pbs]
    for truth in ["TP", "FP"]:
        sub = sub_pbs[sub_pbs["truth_label"] == truth]
        n = len(sub)
        n_loh = sub["core_loh_like"].sum() if n > 0 else 0
        rows_q3b.append({
            "phase_block_status": pbs, "truth": truth,
            "n": n, "n_loh": int(n_loh),
            "loh_frac": round(n_loh / n, 4) if n > 0 else np.nan,
        })
df_q3b = pd.DataFrame(rows_q3b)
write_tsv(df_q3b, OUT_DIR / "q3b_phase_block_status_cross.tsv", "Phase block status × truth × LOH-like")

# --- Q3c: SuggestFilter × LOH × truth (TO only) ---
rows_q3c = []
for sf in [True, False]:
    sub_sf = df_to[df_to["SuggestFilter"] == sf] if "SuggestFilter" in df_to.columns else pd.DataFrame()
    if len(sub_sf) == 0:
        sub_sf = df_to[df_to["SuggestFilter"].astype(str).str.lower() == str(sf).lower()]
    for truth in ["TP", "FP"]:
        sub = sub_sf[sub_sf["truth_label"] == truth]
        n = len(sub)
        n_loh = sub["core_loh_like"].sum() if n > 0 else 0
        rows_q3c.append({
            "suggest_filter": sf, "truth": truth,
            "n": n, "n_loh": int(n_loh),
            "loh_frac": round(n_loh / n, 4) if n > 0 else np.nan,
        })
df_q3c = pd.DataFrame(rows_q3c)
write_tsv(df_q3c, OUT_DIR / "q3c_suggest_filter_cross.tsv", "SuggestFilter × truth × LOH-like (TO)")

# --- Fig06: TO dimensions overview (3 panels) ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: to_loh_bed_hit
ax = axes[0]
sub = df_q3a.copy()
cats = [True, False]
x = np.arange(len(cats))
w = 0.35
for j, truth in enumerate(["TP", "FP"]):
    vals = [sub[(sub["to_loh_bed_hit"] == c) & (sub["truth"] == truth)]["loh_frac"].values[0]
            if len(sub[(sub["to_loh_bed_hit"] == c) & (sub["truth"] == truth)]) > 0 else 0
            for c in cats]
    offset = -w/2 if j == 0 else w/2
    bars = ax.bar(x + offset, vals, w, color=TRUTH_COLORS[truth], label=truth, edgecolor="black", linewidth=0.3)
    for k, v in enumerate(vals):
        ax.text(x[k] + offset, v + 0.01, f'{v:.3f}', ha="center", fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(["In LOH BED", "Not in LOH BED"])
ax.set_ylabel("LOH-like Fraction")
ax.set_title("to_loh_bed_hit", fontweight="bold")
ax.legend(fontsize=8)

# Panel 2: phase_block_status
ax = axes[1]
sub = df_q3b.copy()
cats = sub["phase_block_status"].unique()
x = np.arange(len(cats))
for j, truth in enumerate(["TP", "FP"]):
    vals = [sub[(sub["phase_block_status"] == c) & (sub["truth"] == truth)]["loh_frac"].values[0]
            if len(sub[(sub["phase_block_status"] == c) & (sub["truth"] == truth)]) > 0 else 0
            for c in cats]
    offset = -w/2 if j == 0 else w/2
    bars = ax.bar(x + offset, vals, w, color=TRUTH_COLORS[truth], label=truth, edgecolor="black", linewidth=0.3)
    for k, v in enumerate(vals):
        ax.text(x[k] + offset, v + 0.01, f'{v:.3f}', ha="center", fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(cats, fontsize=8)
ax.set_ylabel("")
ax.set_title("phase_block_status", fontweight="bold")
ax.legend(fontsize=8)

# Panel 3: SuggestFilter
ax = axes[2]
sub = df_q3c.copy()
cats = [True, False]
x = np.arange(len(cats))
for j, truth in enumerate(["TP", "FP"]):
    vals = [sub[(sub["suggest_filter"] == c) & (sub["truth"] == truth)]["loh_frac"].values[0]
            if len(sub[(sub["suggest_filter"] == c) & (sub["truth"] == truth)]) > 0 else 0
            for c in cats]
    offset = -w/2 if j == 0 else w/2
    bars = ax.bar(x + offset, vals, w, color=TRUTH_COLORS[truth], label=truth, edgecolor="black", linewidth=0.3)
    for k, v in enumerate(vals):
        ax.text(x[k] + offset, v + 0.01, f'{v:.3f}', ha="center", fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(["SuggestFilter=True", "SuggestFilter=False"])
ax.set_ylabel("")
ax.set_title("SuggestFilter", fontweight="bold")
ax.legend(fontsize=8)

fig.suptitle("Q3: TO-Specific Dimensions — LOH-like Rate by Truth", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig06_to_dimensions_loh_rate.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig06_to_dimensions_loh_rate.png")


# ============================================================
# Q4: VerificationClass × LOH × Tier 3D cross
# ============================================================
print(f"\n[4/5] Q4 — VerificationClass × LOH × Tier ...")

# --- Q4a: Full 3D cross ---
rows_q4a = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for vc in VC_ORDER:
        for tier in TIER_ORDER:
            dt = dm[(dm["verification_class"] == vc) & (dm["tier"] == tier)]
            for truth in ["TP", "FP"]:
                sub = dt[dt["truth_label"] == truth]
                n = len(sub)
                n_loh = sub["core_loh_like"].sum() if n > 0 else 0
                rows_q4a.append({
                    "mode": mode, "verification_class": vc, "tier": tier, "truth": truth,
                    "n": n, "n_loh": int(n_loh),
                    "loh_frac": round(n_loh / n, 4) if n > 0 else np.nan,
                })
df_q4a = pd.DataFrame(rows_q4a)
write_tsv(df_q4a, OUT_DIR / "q4a_vc_tier_truth_loh_rate.tsv", "VerifClass × Tier × Truth LOH rate")

# --- Fig07: Faceted heatmap — FP LOH-like rate (VC × Tier), paired and TO ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_q4a[(df_q4a["mode"] == mode) & (df_q4a["truth"] == "FP")].copy()
    pivot = sub.pivot_table(index="verification_class", columns="tier", values="loh_frac")
    # Reindex
    pivot = pivot.reindex(index=VC_ORDER, columns=TIER_ORDER)
    sns.heatmap(pivot.astype(float), ax=ax, annot=True, fmt=".2f",
                cmap="YlOrRd", vmin=0, vmax=0.8,
                linewidths=0.5, cbar_kws={"label": "FP LOH-like fraction"})
    ax.set_title(f"{mode.upper()} — FP LOH-like Rate", fontsize=12, fontweight="bold")
    ax.set_ylabel("VerificationClass" if i == 0 else "")
fig.suptitle("Q4: FP LOH-like Rate — VerificationClass × Tier", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig07_vc_tier_fp_loh_rate_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig07_vc_tier_fp_loh_rate_heatmap.png")

# --- Fig08: TP LOH-like rate heatmap for comparison ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_q4a[(df_q4a["mode"] == mode) & (df_q4a["truth"] == "TP")].copy()
    pivot = sub.pivot_table(index="verification_class", columns="tier", values="loh_frac")
    pivot = pivot.reindex(index=VC_ORDER, columns=TIER_ORDER)
    sns.heatmap(pivot.astype(float), ax=ax, annot=True, fmt=".2f",
                cmap="YlGn", vmin=0, vmax=0.8,
                linewidths=0.5, cbar_kws={"label": "TP LOH-like fraction"})
    ax.set_title(f"{mode.upper()} — TP LOH-like Rate", fontsize=12, fontweight="bold")
    ax.set_ylabel("VerificationClass" if i == 0 else "")
fig.suptitle("Q4: TP LOH-like Rate — VerificationClass × Tier", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig08_vc_tier_tp_loh_rate_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig08_vc_tier_tp_loh_rate_heatmap.png")

# --- Fig09: Delta (TP - FP) LOH-like rate heatmap ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    tp_sub = df_q4a[(df_q4a["mode"] == mode) & (df_q4a["truth"] == "TP")].copy()
    fp_sub = df_q4a[(df_q4a["mode"] == mode) & (df_q4a["truth"] == "FP")].copy()
    tp_pivot = tp_sub.pivot_table(index="verification_class", columns="tier", values="loh_frac")
    fp_pivot = fp_sub.pivot_table(index="verification_class", columns="tier", values="loh_frac")
    tp_pivot = tp_pivot.reindex(index=VC_ORDER, columns=TIER_ORDER)
    fp_pivot = fp_pivot.reindex(index=VC_ORDER, columns=TIER_ORDER)
    delta = tp_pivot - fp_pivot  # positive = TP has more LOH-like (TP enriched)
    sns.heatmap(delta.astype(float), ax=ax, annot=True, fmt=".2f",
                cmap="RdBu", center=0, vmin=-0.3, vmax=0.3,
                linewidths=0.5, cbar_kws={"label": "TP - FP LOH-like rate (>0 = TP enriched)"})
    ax.set_title(f"{mode.upper()} — Delta (TP - FP) LOH-like Rate", fontsize=12, fontweight="bold")
    ax.set_ylabel("VerificationClass" if i == 0 else "")
fig.suptitle("Q4: TP vs FP LOH-like Rate Delta — VerificationClass × Tier", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig09_vc_tier_delta_loh_rate_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig09_vc_tier_delta_loh_rate_heatmap.png")


# ============================================================
# Q5: Methylation features in LOH-like subgroups
# ============================================================
print(f"\n[5/5] Q5 — Methylation features in LOH-like TP vs FP ...")

FEATURES = ["quality_score", "pairwise_median_dist", "cramers_v", "allele_delta", "caller_af", "hp0_ratio"]
FEATURE_LABELS = {
    "quality_score": "Quality Score",
    "pairwise_median_dist": "Pairwise Median Dist",
    "cramers_v": "Cramér's V",
    "allele_delta": "Allele Delta",
    "caller_af": "Caller AF",
    "hp0_ratio": "HP0 Ratio",
}

# --- Q5a: Feature statistics for LOH-like TP vs FP, by Tier A and A+ ---
rows_q5a = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for tier in ["A", "A+"]:
        dt = dm[(dm["tier"] == tier) & (dm["core_loh_like"] == True)]
        for truth in ["TP", "FP"]:
            sub = dt[dt["truth_label"] == truth]
            if len(sub) < 5:
                continue
            row = {"mode": mode, "tier": tier, "truth": truth, "n": len(sub)}
            for feat in FEATURES:
                vals = sub[feat].dropna()
                if len(vals) > 0:
                    row[f"{feat}_median"] = round(vals.median(), 4)
                    row[f"{feat}_mean"] = round(vals.mean(), 4)
                    row[f"{feat}_q25"] = round(vals.quantile(0.25), 4)
                    row[f"{feat}_q75"] = round(vals.quantile(0.75), 4)
                else:
                    row[f"{feat}_median"] = np.nan
                    row[f"{feat}_mean"] = np.nan
                    row[f"{feat}_q25"] = np.nan
                    row[f"{feat}_q75"] = np.nan
            rows_q5a.append(row)
df_q5a = pd.DataFrame(rows_q5a)
write_tsv(df_q5a, OUT_DIR / "q5a_loh_tp_vs_fp_features_by_tier.tsv", "LOH-like TP vs FP features by Tier A/A+")

# --- Q5b: Also compute for non-LOH as comparison ---
rows_q5b = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for tier in ["A", "A+"]:
        for loh_status in [True, False]:
            dt = dm[(dm["tier"] == tier) & (dm["core_loh_like"] == loh_status)]
            for truth in ["TP", "FP"]:
                sub = dt[dt["truth_label"] == truth]
                if len(sub) < 5:
                    continue
                row = {"mode": mode, "tier": tier, "loh_like": loh_status, "truth": truth, "n": len(sub)}
                for feat in FEATURES:
                    vals = sub[feat].dropna()
                    row[f"{feat}_median"] = round(vals.median(), 4) if len(vals) > 0 else np.nan
                rows_q5b.append(row)
df_q5b = pd.DataFrame(rows_q5b)
write_tsv(df_q5b, OUT_DIR / "q5b_loh_nonloh_tp_fp_features_by_tier.tsv", "LOH/non-LOH × TP/FP features by tier")

# --- Fig10: Violin plots — LOH-like TP vs FP features (Tier A and A+, paired and TO) ---
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    fig, axes = plt.subplots(2, len(FEATURES), figsize=(4 * len(FEATURES), 8))
    for row_i, tier in enumerate(["A", "A+"]):
        dt = dm[(dm["tier"] == tier) & (dm["core_loh_like"] == True)].copy()
        for col_i, feat in enumerate(FEATURES):
            ax = axes[row_i, col_i]
            plot_data = dt[[feat, "truth_label"]].dropna()
            if len(plot_data) < 10:
                ax.text(0.5, 0.5, "n < 10", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"Tier {tier}: {FEATURE_LABELS[feat]}", fontsize=9)
                continue
            parts = ax.violinplot(
                [plot_data[plot_data["truth_label"] == "TP"][feat].values,
                 plot_data[plot_data["truth_label"] == "FP"][feat].values],
                positions=[0, 1], showmedians=True, showextrema=False
            )
            for j, truth in enumerate(["TP", "FP"]):
                if j < len(parts["bodies"]):
                    parts["bodies"][j].set_facecolor(TRUTH_COLORS[truth])
                    parts["bodies"][j].set_alpha(0.7)
            n_tp = (plot_data["truth_label"] == "TP").sum()
            n_fp = (plot_data["truth_label"] == "FP").sum()
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"TP\n(n={n_tp:,})", f"FP\n(n={n_fp:,})"], fontsize=7)
            ax.set_title(f"Tier {tier}: {FEATURE_LABELS[feat]}", fontsize=9)

    fig.suptitle(f"Q5: LOH-like TP vs FP Feature Distributions — {mode.upper()} mode",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"fig10_{mode}_loh_tp_vs_fp_features.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig10_{mode}_loh_tp_vs_fp_features.png")

# --- Fig11: HPMergedSig × LOH × Truth cross-tab (special: 聯合 feature) ---
rows_hpsig = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    for tier in ["A", "A+"]:
        dt = dm[dm["tier"] == tier]
        for loh in [True, False]:
            for hpsig in [True, False]:
                for truth in ["TP", "FP"]:
                    sub = dt[(dt["core_loh_like"] == loh) & (dt["HPMergedSig"] == hpsig) & (dt["truth_label"] == truth)]
                    rows_hpsig.append({
                        "mode": mode, "tier": tier, "loh_like": loh, "HPMergedSig": hpsig,
                        "truth": truth, "n": len(sub),
                    })
df_hpsig = pd.DataFrame(rows_hpsig)
write_tsv(df_hpsig, OUT_DIR / "q5c_loh_hpsig_truth_cross.tsv", "LOH × HPMergedSig × Truth cross-tab by tier")

# Compute enrichment from this cross-tab
rows_hpsig_enrich = []
for mode in ["paired", "to"]:
    for tier in ["A", "A+"]:
        for loh in [True, False]:
            for hpsig in [True, False]:
                tp_n = df_hpsig[(df_hpsig["mode"]==mode) & (df_hpsig["tier"]==tier) &
                                (df_hpsig["loh_like"]==loh) & (df_hpsig["HPMergedSig"]==hpsig) &
                                (df_hpsig["truth"]=="TP")]["n"].values
                fp_n = df_hpsig[(df_hpsig["mode"]==mode) & (df_hpsig["tier"]==tier) &
                                (df_hpsig["loh_like"]==loh) & (df_hpsig["HPMergedSig"]==hpsig) &
                                (df_hpsig["truth"]=="FP")]["n"].values
                tp_n = tp_n[0] if len(tp_n) > 0 else 0
                fp_n = fp_n[0] if len(fp_n) > 0 else 0
                total = tp_n + fp_n
                fp_rate = fp_n / total if total > 0 else np.nan
                rows_hpsig_enrich.append({
                    "mode": mode, "tier": tier, "loh_like": loh, "HPMergedSig": hpsig,
                    "n_tp": tp_n, "n_fp": fp_n, "n_total": total,
                    "fp_rate": round(fp_rate, 4) if not np.isnan(fp_rate) else np.nan,
                })
df_hpsig_enrich = pd.DataFrame(rows_hpsig_enrich)
write_tsv(df_hpsig_enrich, OUT_DIR / "q5d_loh_hpsig_fp_rate.tsv", "LOH × HPMergedSig FP rate by tier")

# --- Fig12: LOH × HPMergedSig FP rate bar chart ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df_hpsig_enrich[df_hpsig_enrich["mode"] == mode].copy()
    # 4 categories per tier
    for j, tier in enumerate(["A", "A+"]):
        tier_sub = sub[sub["tier"] == tier]
        cats = []
        vals = []
        colors = []
        for loh in [False, True]:
            for hpsig in [False, True]:
                row = tier_sub[(tier_sub["loh_like"] == loh) & (tier_sub["HPMergedSig"] == hpsig)]
                label = f"{'LOH' if loh else 'noLOH'}+{'Sig' if hpsig else 'noSig'}"
                cats.append(label)
                vals.append(row["fp_rate"].values[0] if len(row) > 0 and not np.isnan(row["fp_rate"].values[0]) else 0)
                # Color by LOH
                colors.append(LOH_COLORS[loh] if not hpsig else ("#ff7f0e" if loh else "#9467bd"))

        x = np.arange(len(cats)) + j * (len(cats) + 1)
        bars = ax.bar(x, vals, color=colors, edgecolor="black", linewidth=0.3)
        for k, v in enumerate(vals):
            n_row = tier_sub[(tier_sub["loh_like"] == ([False, False, True, True][k])) &
                             (tier_sub["HPMergedSig"] == ([False, True, False, True][k]))]
            n_total = n_row["n_total"].values[0] if len(n_row) > 0 else 0
            ax.text(x[k], v + 0.005, f'{v:.3f}\n(n={n_total:,})', ha="center", fontsize=6)
        ax.text(x.mean(), -0.05, f"Tier {tier}", ha="center", fontsize=10, fontweight="bold",
                transform=ax.get_xaxis_transform())

    all_x = list(range(len(cats))) + list(range(len(cats) + 1, 2 * len(cats) + 1))
    ax.set_xticks(all_x)
    ax.set_xticklabels(cats * 2, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("FP Rate" if i == 0 else "")
    ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")

fig.suptitle("Q5: LOH × HPMergedSig — FP Rate by Tier", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(FIG_DIR / "fig12_loh_hpsig_fp_rate.png", dpi=150, bbox_inches="tight")
plt.close()
print("  fig12_loh_hpsig_fp_rate.png")


# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print(f"LOH Expanded Observation complete!")
print(f"  Output: {OUT_DIR}")
print(f"  TSV files: {len(list(OUT_DIR.glob('*.tsv')))}")
print(f"  Figures: {len(list(FIG_DIR.glob('*.png')))}")
print(f"  Generated at: {GENERATED_AT}")
print(f"{'='*60}")
