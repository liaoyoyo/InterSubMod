#!/usr/bin/env python3
"""
Post-HP-Fix TO LOH Deep Investigation

核心問題：
  Q1: 為什麼 TO 和 paired LOH enrichment 方向相反？
      (TO: 0.805 TP enriched, paired: 1.194 FP enriched)
  Q2: TO LOH-like 能否用於 TP rescue？
  Q3: HP fix 後 Tier 結構重塑的完整視覺化

生成時間: 2026-03-30
依賴: 20260327_loh_round1_cross_sample_audit (post-HP-fix master)
"""

import os
import sys
import json
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

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import SchemaContractError, select_current_view

warnings.filterwarnings("ignore")

# ============================================================
# Configuration
# ============================================================
ROUND1_WS = Path(
    os.environ.get(
        "LOH_ROUND1_WS",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
        "20260327_loh_round1_cross_sample_audit",
    )
)
OLD_ROUND1_WS = Path(
    os.environ.get(
        "LOH_OLD_ROUND1_WS",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
        "20260327_loh_round1_cross_sample_audit_before_hp_fix",
    )
)
OUT_DIR = Path(
    os.environ.get(
        "INVESTIGATION_OUT_DIR",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
        "20260330_post_hp_fix_to_loh_investigation",
    )
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
TIER_ORDER = ["C0", "C", "B", "A", "A+"]
TIER_COLORS = {"C0": "#d9d9d9", "C": "#fc8d59", "B": "#fee090", "A": "#4575b4", "A+": "#313695"}


def assign_tier(n):
    if pd.isna(n): return "C0"
    n = int(n)
    if n == 0: return "C0"
    elif n < 10: return "C"
    elif n < 30: return "B"
    elif n < 50: return "A"
    else: return "A+"


def write_tsv(df, path, description=""):
    df.to_csv(path, sep="\t", index=False)
    print(f"  {path.name} ({len(df)} rows, {os.path.getsize(path):,} bytes) — {description}")


def attach_current_verification_view(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate the versioned C2 view and any derived lowercase exact alias."""
    view = select_current_view(frame)
    if "verification_class" in frame.columns:
        canonical = frame["VerificationClass"].astype("string")
        derived = frame["verification_class"].astype("string")
        mismatch = canonical.fillna("<MISSING>") != derived.fillna("<MISSING>")
        if mismatch.any():
            raise SchemaContractError(
                "post-HP-fix master: verification_class is not an exact alias of "
                f"VerificationClass at {int(mismatch.sum())} rows"
            )
    prepared = frame.copy()
    prepared["_verification_class_current"] = view.values
    prepared.attrs["verification_schema_metadata"] = {
        "selection_field": view.field,
        "schema_status": view.schema_status,
        "unknown_current_class_count": sum(view.unknown_counts.values()),
        "unknown_current_class_values": "; ".join(
            f"{key}:{value}" for key, value in sorted(view.unknown_counts.items())
        ),
    }
    return prepared


# ============================================================
# Data Loading
# ============================================================
print(f"[1/7] Loading master dataset from {ROUND1_WS / 'all_region_rows.tsv.gz'} ...")
df = attach_current_verification_view(
    pd.read_csv(ROUND1_WS / "all_region_rows.tsv.gz", sep="\t", low_memory=False)
)
verification_schema_metadata = dict(df.attrs["verification_schema_metadata"])
df["tier"] = df["effective_hp_reads"].apply(assign_tier)
print(f"  Loaded {len(df):,} rows ({df['mode'].value_counts().to_dict()})")
pd.DataFrame([verification_schema_metadata]).to_csv(
    OUT_DIR / "verification_schema_provenance.tsv", sep="\t", index=False
)

# Load old master for comparison
print(f"[1b/7] Loading old (pre-fix) master for comparison ...")
df_old = pd.read_csv(OLD_ROUND1_WS / "all_region_rows.tsv.gz", sep="\t", low_memory=False)
df_old["tier"] = df_old["effective_hp_reads"].apply(assign_tier)
print(f"  Loaded {len(df_old):,} rows")

# ============================================================
# Q1: Why are TO and paired LOH enrichment directions opposite?
# ============================================================
print("\n[2/7] Q1 — Investigating opposite enrichment directions ...")

# --- Analysis 1a: AF distribution by mode × truth × LOH status ---
print("  [2a] AF distribution by mode × truth × LOH status")
rows_q1a = []
for mode in ["paired", "to"]:
    for truth in ["TP", "FP"]:
        for loh in [True, False]:
            sub = df[(df["mode"] == mode) & (df["truth_label"] == truth) & (df["core_loh_like"] == loh)]
            if len(sub) == 0:
                continue
            af = sub["caller_af"]
            rows_q1a.append({
                "mode": mode, "truth": truth, "loh_like": loh,
                "n": len(sub),
                "af_median": af.median(), "af_mean": af.mean(),
                "af_q25": af.quantile(0.25), "af_q75": af.quantile(0.75),
                "eff_hp_median": sub["effective_hp_reads"].median(),
                "hp0_ratio_median": sub["hp0_ratio"].median(),
                "hp_ratio_core_median": sub["hp_ratio_core"].median(),
                "allele_delta_median": sub["allele_delta"].median() if "allele_delta" in sub.columns else np.nan,
                "pairwise_dist_median": sub["pairwise_median_dist"].median(),
            })
q1a_df = pd.DataFrame(rows_q1a)
write_tsv(q1a_df, OUT_DIR / "q1a_af_by_mode_truth_loh.tsv", "AF/feature distributions by mode×truth×LOH")

# --- Analysis 1b: Per-sample enrichment comparison ---
print("  [2b] Per-sample enrichment: paired vs TO side-by-side")
rows_q1b = []
for sample in SAMPLES:
    for mode in ["paired", "to"]:
        sub = df[(df["sample"] == sample) & (df["mode"] == mode)]
        tp = sub[sub["truth_label"] == "TP"]
        fp = sub[sub["truth_label"] == "FP"]
        if len(tp) == 0 or len(fp) == 0:
            continue
        tp_loh = tp["core_loh_like"].sum()
        fp_loh = fp["core_loh_like"].sum()
        tp_frac = tp_loh / len(tp)
        fp_frac = fp_loh / len(fp)
        enrich = fp_frac / tp_frac if tp_frac > 0 else np.nan
        rows_q1b.append({
            "sample": sample, "mode": mode,
            "tp_total": len(tp), "fp_total": len(fp),
            "tp_loh_n": int(tp_loh), "fp_loh_n": int(fp_loh),
            "tp_loh_frac": tp_frac, "fp_loh_frac": fp_frac,
            "enrichment": enrich,
            "tp_af_median": tp["caller_af"].median(),
            "fp_af_median": fp["caller_af"].median(),
            "tp_eff_hp_median": tp["effective_hp_reads"].median(),
            "fp_eff_hp_median": fp["effective_hp_reads"].median(),
        })
q1b_df = pd.DataFrame(rows_q1b)
write_tsv(q1b_df, OUT_DIR / "q1b_enrichment_paired_vs_to.tsv", "per-sample enrichment comparison")

# --- Analysis 1c: Same-locus paired vs TO LOH concordance ---
print("  [2c] Same-locus LOH concordance (why different direction?)")
slc_path = ROUND1_WS / "same_locus_compare.tsv"
if slc_path.exists():
    slc = pd.read_csv(slc_path, sep="\t", low_memory=False)
    # Identify same-locus pairs
    rows_q1c = []
    for sample in SAMPLES:
        sub = slc[slc["sample"] == sample] if "sample" in slc.columns else pd.DataFrame()
        if len(sub) == 0:
            continue
        # Check if we have paired/to LOH columns
        p_loh_col = [c for c in sub.columns if "paired" in c.lower() and "loh" in c.lower()]
        t_loh_col = [c for c in sub.columns if "to" in c.lower() and "loh" in c.lower()]
        if p_loh_col and t_loh_col:
            p_col = p_loh_col[0]
            t_col = t_loh_col[0]
            both_loh = ((sub[p_col] == True) & (sub[t_col] == True)).sum()
            paired_only = ((sub[p_col] == True) & (sub[t_col] == False)).sum()
            to_only = ((sub[p_col] == False) & (sub[t_col] == True)).sum()
            neither = ((sub[p_col] == False) & (sub[t_col] == False)).sum()
            rows_q1c.append({
                "sample": sample, "total_loci": len(sub),
                "both_loh": both_loh, "paired_only_loh": paired_only,
                "to_only_loh": to_only, "neither_loh": neither,
                "concordance": (both_loh + neither) / len(sub) if len(sub) > 0 else np.nan,
            })
    if rows_q1c:
        q1c_df = pd.DataFrame(rows_q1c)
        write_tsv(q1c_df, OUT_DIR / "q1c_same_locus_loh_concordance.tsv", "same-locus LOH concordance")
    else:
        print("    (same_locus_compare columns not compatible, skipping)")
else:
    print("    (same_locus_compare.tsv not found, skipping)")

# --- Analysis 1d: FP origin hypothesis — germline vs somatic FP ---
print("  [2d] TO FP origin: caller_filter × LOH status")
rows_q1d = []
for mode in ["paired", "to"]:
    fp = df[(df["mode"] == mode) & (df["truth_label"] == "FP")]
    if len(fp) == 0:
        continue
    for filt_val in fp["caller_filter"].dropna().unique():
        sub = fp[fp["caller_filter"] == filt_val]
        loh_n = sub["core_loh_like"].sum()
        rows_q1d.append({
            "mode": mode, "caller_filter": filt_val,
            "n": len(sub), "loh_like_n": int(loh_n),
            "loh_like_frac": loh_n / len(sub),
            "af_median": sub["caller_af"].median(),
            "eff_hp_median": sub["effective_hp_reads"].median(),
        })
if rows_q1d:
    q1d_df = pd.DataFrame(rows_q1d).sort_values(["mode", "n"], ascending=[True, False])
    write_tsv(q1d_df, OUT_DIR / "q1d_fp_caller_filter_loh.tsv", "FP by caller_filter × LOH")

# --- Analysis 1e: HP balance as explanatory variable ---
print("  [2e] HP balance (hp_ratio_core) distribution by mode × truth")
rows_q1e = []
for mode in ["paired", "to"]:
    for truth in ["TP", "FP"]:
        sub = df[(df["mode"] == mode) & (df["truth_label"] == truth)]
        hrc = sub["hp_ratio_core"].dropna()
        if len(hrc) == 0:
            continue
        # HP ratio deviation from 0.5 (balance)
        hp_dev = (hrc - 0.5).abs()
        rows_q1e.append({
            "mode": mode, "truth": truth, "n": len(hrc),
            "hp_ratio_median": hrc.median(),
            "hp_ratio_mean": hrc.mean(),
            "hp_deviation_median": hp_dev.median(),
            "hp_deviation_mean": hp_dev.mean(),
            "pct_extreme_loh": ((hrc < 0.1) | (hrc > 0.9)).mean(),
            "pct_moderate_imbalance": (((hrc >= 0.1) & (hrc < 0.3)) | ((hrc > 0.7) & (hrc <= 0.9))).mean(),
            "pct_balanced": ((hrc >= 0.3) & (hrc <= 0.7)).mean(),
        })
q1e_df = pd.DataFrame(rows_q1e)
write_tsv(q1e_df, OUT_DIR / "q1e_hp_balance_by_mode_truth.tsv", "HP balance distribution")

# ============================================================
# Q2: Can TO LOH be used for TP rescue?
# ============================================================
print("\n[3/7] Q2 — TO LOH as TP rescue candidate ...")

# --- Analysis 2a: LOH-like fraction by VerificationClass ---
print("  [3a] LOH-like fraction by VerificationClass (TO)")
to_df = df[df["mode"] == "to"]
rows_q2a = []
for vc in to_df["_verification_class_current"].dropna().unique():
    for truth in ["TP", "FP"]:
        sub = to_df[(to_df["_verification_class_current"] == vc) & (to_df["truth_label"] == truth)]
        if len(sub) == 0:
            continue
        loh_n = sub["core_loh_like"].sum()
        rows_q2a.append({
            "verification_class": vc, "truth": truth,
            "n": len(sub), "loh_like_n": int(loh_n),
            "loh_like_frac": loh_n / len(sub),
            "af_median": sub["caller_af"].median(),
            "tier_A_pct": (sub["tier"].isin(["A", "A+"])).mean(),
        })
q2a_df = pd.DataFrame(rows_q2a)
if len(q2a_df) > 0:
    q2a_df = q2a_df.sort_values(["verification_class", "truth"])
    write_tsv(q2a_df, OUT_DIR / "q2a_loh_by_verification_class.tsv", "LOH by VerificationClass")

# --- Analysis 2b: Theoretical rescue: if LOH-like TP had been filtered, how many TP rescued ---
print("  [3b] Quality score distribution: LOH-like TP vs non-LOH-like TP (TO)")
rows_q2b = []
for sample in SAMPLES:
    to_s = to_df[to_df["sample"] == sample]
    tp = to_s[to_s["truth_label"] == "TP"]
    tp_loh = tp[tp["core_loh_like"] == True]
    tp_noloh = tp[tp["core_loh_like"] == False]
    if len(tp_loh) == 0 or len(tp_noloh) == 0:
        continue
    rows_q2b.append({
        "sample": sample,
        "tp_loh_n": len(tp_loh), "tp_noloh_n": len(tp_noloh),
        "tp_loh_af_median": tp_loh["caller_af"].median(),
        "tp_noloh_af_median": tp_noloh["caller_af"].median(),
        "tp_loh_qual_median": tp_loh["quality_score"].median(),
        "tp_noloh_qual_median": tp_noloh["quality_score"].median(),
        "tp_loh_eff_hp_median": tp_loh["effective_hp_reads"].median(),
        "tp_noloh_eff_hp_median": tp_noloh["effective_hp_reads"].median(),
        "tp_loh_hp0_median": tp_loh["hp0_ratio"].median(),
        "tp_noloh_hp0_median": tp_noloh["hp0_ratio"].median(),
        "tp_loh_tierA_pct": tp_loh["tier"].isin(["A", "A+"]).mean(),
        "tp_noloh_tierA_pct": tp_noloh["tier"].isin(["A", "A+"]).mean(),
    })
q2b_df = pd.DataFrame(rows_q2b)
write_tsv(q2b_df, OUT_DIR / "q2b_tp_loh_vs_noloh_quality.tsv", "TO TP: LOH vs non-LOH quality")

# --- Analysis 2c: LOH-like among low-AF TP (potential rescue targets) ---
print("  [3c] LOH-like among low-AF TP (caller_af < 0.2)")
rows_q2c = []
af_thresholds = [0.1, 0.15, 0.2, 0.3]
for af_thr in af_thresholds:
    for mode in ["paired", "to"]:
        sub = df[(df["mode"] == mode) & (df["truth_label"] == "TP") & (df["caller_af"] < af_thr)]
        if len(sub) == 0:
            continue
        loh_n = sub["core_loh_like"].sum()
        rows_q2c.append({
            "mode": mode, "af_threshold": af_thr,
            "tp_count": len(sub), "tp_loh_n": int(loh_n),
            "tp_loh_frac": loh_n / len(sub),
            "eff_hp_median": sub["effective_hp_reads"].median(),
        })
q2c_df = pd.DataFrame(rows_q2c)
write_tsv(q2c_df, OUT_DIR / "q2c_low_af_tp_loh_fraction.tsv", "LOH fraction among low-AF TP")

# ============================================================
# Q3: HP-fix Tier structure visualization data
# ============================================================
print("\n[4/7] Q3 — Tier structure: before vs after HP fix ...")

# --- Analysis 3a: Tier distribution before/after ---
rows_q3a = []
for label, dset in [("before_fix", df_old), ("after_fix", df)]:
    for mode in ["paired", "to"]:
        sub = dset[dset["mode"] == mode]
        total = len(sub)
        for tier in TIER_ORDER:
            n = (sub["tier"] == tier).sum()
            rows_q3a.append({
                "dataset": label, "mode": mode, "tier": tier,
                "n": n, "fraction": n / total if total > 0 else 0,
            })
q3a_df = pd.DataFrame(rows_q3a)
write_tsv(q3a_df, OUT_DIR / "q3a_tier_distribution_before_after.tsv", "Tier distribution before/after fix")

# --- Analysis 3b: Per-sample tier shift ---
rows_q3b = []
for sample in SAMPLES:
    for label, dset in [("before_fix", df_old), ("after_fix", df)]:
        to_s = dset[(dset["mode"] == "to") & (dset["sample"] == sample)]
        total = len(to_s)
        for tier in TIER_ORDER:
            n = (to_s["tier"] == tier).sum()
            rows_q3b.append({
                "dataset": label, "sample": sample, "tier": tier,
                "n": n, "fraction": n / total if total > 0 else 0,
            })
q3b_df = pd.DataFrame(rows_q3b)
write_tsv(q3b_df, OUT_DIR / "q3b_per_sample_tier_shift.tsv", "Per-sample TO tier shift")

# --- Analysis 3c: Enrichment by tier before/after ---
rows_q3c = []
for label, dset in [("before_fix", df_old), ("after_fix", df)]:
    for mode in ["paired", "to"]:
        for tier in TIER_ORDER:
            sub = dset[(dset["mode"] == mode) & (dset["tier"] == tier)]
            tp = sub[sub["truth_label"] == "TP"]
            fp = sub[sub["truth_label"] == "FP"]
            if len(tp) == 0 or len(fp) == 0:
                continue
            tp_loh = tp["core_loh_like"].sum()
            fp_loh = fp["core_loh_like"].sum()
            tp_frac = tp_loh / len(tp) if len(tp) > 0 else 0
            fp_frac = fp_loh / len(fp) if len(fp) > 0 else 0
            enrich = fp_frac / tp_frac if tp_frac > 0 else np.nan
            rows_q3c.append({
                "dataset": label, "mode": mode, "tier": tier,
                "tp_total": len(tp), "fp_total": len(fp),
                "tp_loh_frac": tp_frac, "fp_loh_frac": fp_frac,
                "enrichment": enrich,
            })
q3c_df = pd.DataFrame(rows_q3c)
write_tsv(q3c_df, OUT_DIR / "q3c_enrichment_by_tier_before_after.tsv", "Enrichment by tier before/after fix")

# ============================================================
# Figure Generation
# ============================================================
print("\n[5/7] Generating figures ...")

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10,
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
})

# --- Fig 01: Enrichment paired vs TO by sample (butterfly/diverging bar) ---
print("  fig01: Enrichment butterfly chart")
fig, ax = plt.subplots(figsize=(10, 6))
q1b_sorted = q1b_df.copy()
paired_data = q1b_sorted[q1b_sorted["mode"] == "paired"].set_index("sample")["enrichment"]
to_data = q1b_sorted[q1b_sorted["mode"] == "to"].set_index("sample")["enrichment"]
y_pos = np.arange(len(SAMPLES))
bars_p = ax.barh(y_pos + 0.18, paired_data.reindex(SAMPLES).values, 0.35,
                 color=MODE_COLORS["paired"], label="Paired", alpha=0.85)
bars_t = ax.barh(y_pos - 0.18, to_data.reindex(SAMPLES).values, 0.35,
                 color=MODE_COLORS["to"], label="TO", alpha=0.85)
ax.axvline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels([SAMPLE_LABELS.get(s, s) for s in SAMPLES])
ax.set_xlabel("FP/TP LOH-like Enrichment")
ax.set_title("Q1: LOH Enrichment — Paired vs TO\n(< 1.0 = TP enriched, > 1.0 = FP enriched)")
ax.legend(loc="lower right")
for bar, val in zip(bars_p, paired_data.reindex(SAMPLES).values):
    if not np.isnan(val):
        ax.text(val + 0.02, bar.get_y() + bar.get_height()/2, f"{val:.2f}",
                va="center", fontsize=8, color=MODE_COLORS["paired"])
for bar, val in zip(bars_t, to_data.reindex(SAMPLES).values):
    if not np.isnan(val):
        ax.text(val + 0.02, bar.get_y() + bar.get_height()/2, f"{val:.2f}",
                va="center", fontsize=8, color=MODE_COLORS["to"])
fig.savefig(FIG_DIR / "fig01_enrichment_paired_vs_to.png")
plt.close(fig)
print("    ✓ fig01_enrichment_paired_vs_to.png")

# --- Fig 02: AF distribution by mode × truth × LOH status ---
print("  fig02: AF violin by mode × truth × LOH")
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
for i, mode in enumerate(["paired", "to"]):
    ax = axes[i]
    sub = df[df["mode"] == mode].copy()
    sub["group"] = sub["truth_label"] + "\n" + sub["core_loh_like"].map({True: "LOH", False: "non-LOH"})
    order = ["TP\nnon-LOH", "TP\nLOH", "FP\nnon-LOH", "FP\nLOH"]
    colors = ["#66c2a5", "#1b7837", "#fc8d62", "#d73027"]
    # Sample for speed
    plot_sub = sub.sample(min(50000, len(sub)), random_state=42)
    sns.violinplot(data=plot_sub, x="group", y="caller_af", order=order,
                   palette=colors, ax=ax, cut=0, inner="quartile", scale="width")
    ax.set_title(f"{mode.upper()}", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Caller AF" if i == 0 else "")
    ax.set_ylim(-0.05, 1.05)
fig.suptitle("Q1: AF Distribution by Truth × LOH Status", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig02_af_violin_by_truth_loh.png")
plt.close(fig)
print("    ✓ fig02_af_violin_by_truth_loh.png")

# --- Fig 03: HP ratio core distribution by mode × truth ---
print("  fig03: HP ratio core histogram by mode × truth")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
for i, mode in enumerate(["paired", "to"]):
    for j, truth in enumerate(["TP", "FP"]):
        ax = axes[i][j]
        sub = df[(df["mode"] == mode) & (df["truth_label"] == truth)]
        hrc = sub["hp_ratio_core"].dropna()
        ax.hist(hrc, bins=100, range=(0, 1), alpha=0.7,
                color=TRUTH_COLORS[truth], edgecolor="none")
        loh_pct = ((hrc < 0.1) | (hrc > 0.9)).mean() * 100
        balanced_pct = ((hrc >= 0.3) & (hrc <= 0.7)).mean() * 100
        ax.axvline(0.1, color="red", linestyle=":", alpha=0.5)
        ax.axvline(0.9, color="red", linestyle=":", alpha=0.5)
        ax.set_title(f"{mode.upper()} {truth}  (LOH-like: {loh_pct:.1f}%, Balanced: {balanced_pct:.1f}%)")
        ax.set_xlabel("hp_ratio_core")
        ax.set_ylabel("Count")
fig.suptitle("Q1: HP Ratio Core Distribution\n(Red lines = LOH threshold 0.1/0.9)", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig03_hp_ratio_core_hist.png")
plt.close(fig)
print("    ✓ fig03_hp_ratio_core_hist.png")

# --- Fig 04: Tier structure before/after HP fix ---
print("  fig04: Tier structure before/after HP fix")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, mode in enumerate(["paired", "to"]):
    ax = axes[idx]
    for k, label in enumerate(["before_fix", "after_fix"]):
        sub = q3a_df[(q3a_df["dataset"] == label) & (q3a_df["mode"] == mode)]
        sub = sub.set_index("tier").reindex(TIER_ORDER)
        bottom = 0
        for tier in TIER_ORDER:
            val = sub.loc[tier, "fraction"] if tier in sub.index else 0
            ax.bar(k, val, bottom=bottom, color=TIER_COLORS[tier], edgecolor="white", linewidth=0.5)
            if val > 0.03:
                ax.text(k, bottom + val / 2, f"{val:.1%}", ha="center", va="center", fontsize=9,
                        color="white" if tier in ["A", "A+"] else "black", fontweight="bold")
            bottom += val
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Before Fix\n(buggy)", "After Fix\n(corrected)"])
    ax.set_ylabel("Fraction")
    ax.set_title(f"{mode.upper()}", fontweight="bold")
    ax.set_ylim(0, 1.02)

# Legend
patches = [mpatches.Patch(color=TIER_COLORS[t], label=f"Tier {t}") for t in TIER_ORDER]
fig.legend(handles=patches, loc="center right", bbox_to_anchor=(1.0, 0.5))
fig.suptitle("Q3: Tier Distribution Before vs After HP Fix", fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 0.92, 0.95])
fig.savefig(FIG_DIR / "fig04_tier_before_after.png")
plt.close(fig)
print("    ✓ fig04_tier_before_after.png")

# --- Fig 05: Enrichment by tier before/after ---
print("  fig05: Enrichment by tier before/after")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for idx, mode in enumerate(["paired", "to"]):
    ax = axes[idx]
    for label, marker, ls in [("before_fix", "o", "--"), ("after_fix", "s", "-")]:
        sub = q3c_df[(q3c_df["dataset"] == label) & (q3c_df["mode"] == mode)]
        sub = sub[sub["tier"].isin(["C", "B", "A", "A+"])]  # skip C0 (NaN enrichment)
        tier_idx = [TIER_ORDER.index(t) for t in sub["tier"]]
        ax.plot(tier_idx, sub["enrichment"], marker=marker, linestyle=ls,
                label=label.replace("_", " ").title(), markersize=8, linewidth=2)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xticks(range(len(TIER_ORDER)))
    ax.set_xticklabels(TIER_ORDER)
    ax.set_xlabel("Support Tier")
    ax.set_ylabel("FP/TP Enrichment")
    ax.set_title(f"{mode.upper()}", fontweight="bold")
    ax.legend()
fig.suptitle("Q3: LOH Enrichment by Tier — Before vs After HP Fix\n(< 1.0 = TP enriched, > 1.0 = FP enriched)",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig05_enrichment_by_tier_before_after.png")
plt.close(fig)
print("    ✓ fig05_enrichment_by_tier_before_after.png")

# --- Fig 06: Per-sample TO tier shift (grouped bar) ---
print("  fig06: Per-sample TO tier shift")
fig, axes = plt.subplots(2, 4, figsize=(18, 8))
axes_flat = axes.flatten()
for i, sample in enumerate(SAMPLES):
    ax = axes_flat[i]
    for k, label in enumerate(["before_fix", "after_fix"]):
        sub = q3b_df[(q3b_df["dataset"] == label) & (q3b_df["sample"] == sample)]
        sub = sub.set_index("tier").reindex(TIER_ORDER)
        x = np.arange(len(TIER_ORDER))
        offset = -0.2 + k * 0.4
        bars = ax.bar(x + offset, sub["fraction"].values, 0.35,
                       color=["#bdbdbd", "#2166ac"][k], alpha=0.85,
                       label=label.replace("_", " ").title() if i == 0 else "")
    ax.set_xticks(np.arange(len(TIER_ORDER)))
    ax.set_xticklabels(TIER_ORDER)
    ax.set_title(SAMPLE_LABELS.get(sample, sample), fontsize=10)
    ax.set_ylim(0, 1.0)
    if i % 4 == 0:
        ax.set_ylabel("Fraction")
# Hide last empty subplot
if len(SAMPLES) < len(axes_flat):
    for j in range(len(SAMPLES), len(axes_flat)):
        axes_flat[j].axis("off")
handles = [mpatches.Patch(color="#bdbdbd", label="Before Fix"),
           mpatches.Patch(color="#2166ac", label="After Fix")]
fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.98, 0.05))
fig.suptitle("Q3: Per-Sample TO Tier Distribution Before vs After HP Fix", fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(FIG_DIR / "fig06_per_sample_tier_shift.png")
plt.close(fig)
print("    ✓ fig06_per_sample_tier_shift.png")

# --- Fig 07: TO TP LOH vs non-LOH feature comparison ---
print("  fig07: TO TP feature comparison: LOH vs non-LOH")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
features = [
    ("caller_af", "Caller AF", (0, 1)),
    ("quality_score", "Quality Score", None),
    ("effective_hp_reads", "Effective HP Reads", (0, 200)),
    ("hp0_ratio", "HP0 Ratio", (0, 0.5)),
    ("pairwise_median_dist", "Pairwise Median Dist", (0, 0.5)),
    ("allele_delta", "Allele Delta", (0, 0.5)),
]
to_tp = to_df[to_df["truth_label"] == "TP"]
to_tp_sample = to_tp.sample(min(80000, len(to_tp)), random_state=42)
for idx, (feat, label, xlim) in enumerate(features):
    ax = axes[idx // 3][idx % 3]
    loh = to_tp_sample[to_tp_sample["core_loh_like"] == True][feat].dropna()
    noloh = to_tp_sample[to_tp_sample["core_loh_like"] == False][feat].dropna()
    if len(loh) > 0 and len(noloh) > 0:
        bins = 80
        rng = xlim if xlim else (min(loh.min(), noloh.min()), max(loh.quantile(0.99), noloh.quantile(0.99)))
        ax.hist(noloh, bins=bins, range=rng, alpha=0.6, density=True,
                color="#66c2a5", label=f"non-LOH (n={len(noloh):,})")
        ax.hist(loh, bins=bins, range=rng, alpha=0.6, density=True,
                color="#d73027", label=f"LOH (n={len(loh):,})")
        # Mann-Whitney U test
        try:
            u_stat, u_p = stats.mannwhitneyu(loh, noloh, alternative="two-sided")
            ax.text(0.98, 0.95, f"MW p={u_p:.2e}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, style="italic")
        except Exception:
            pass
    ax.set_xlabel(label)
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)
fig.suptitle("Q2: TO TP Feature Comparison — LOH-like vs non-LOH-like", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig07_to_tp_loh_vs_noloh_features.png")
plt.close(fig)
print("    ✓ fig07_to_tp_loh_vs_noloh_features.png")

# --- Fig 08: HP balance explanation — why opposite directions ---
print("  fig08: HP balance explanation scatter")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for idx, mode in enumerate(["paired", "to"]):
    ax = axes[idx]
    sub = df[df["mode"] == mode].sample(min(40000, len(df[df["mode"]==mode])), random_state=42)
    tp_sub = sub[sub["truth_label"] == "TP"]
    fp_sub = sub[sub["truth_label"] == "FP"]
    ax.scatter(tp_sub["caller_af"], tp_sub["hp_ratio_core"], alpha=0.03, s=3,
               color=TRUTH_COLORS["TP"], label="TP", rasterized=True)
    ax.scatter(fp_sub["caller_af"], fp_sub["hp_ratio_core"], alpha=0.05, s=3,
               color=TRUTH_COLORS["FP"], label="FP", rasterized=True)
    ax.axhline(0.1, color="red", linestyle=":", alpha=0.3)
    ax.axhline(0.9, color="red", linestyle=":", alpha=0.3)
    ax.set_xlabel("Caller AF")
    ax.set_ylabel("hp_ratio_core")
    ax.set_title(f"{mode.upper()}", fontweight="bold")
    ax.legend(markerscale=10, loc="upper left")
fig.suptitle("Q1: AF vs HP Ratio — TP and FP Distributions\n(Red lines = LOH threshold; LOH-like = above 0.9 or below 0.1)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig08_af_vs_hp_ratio_scatter.png")
plt.close(fig)
print("    ✓ fig08_af_vs_hp_ratio_scatter.png")

# --- Fig 09: Low-AF TP LOH rescue potential ---
print("  fig09: Low-AF TP LOH rescue potential")
fig, ax = plt.subplots(figsize=(10, 6))
for mode, color in MODE_COLORS.items():
    sub = q2c_df[q2c_df["mode"] == mode]
    ax.plot(sub["af_threshold"], sub["tp_loh_frac"], marker="o", color=color,
            linewidth=2, markersize=8, label=mode.upper())
    for _, row in sub.iterrows():
        ax.annotate(f"n={row['tp_count']:,.0f}", (row["af_threshold"], row["tp_loh_frac"]),
                    textcoords="offset points", xytext=(5, 10), fontsize=8, color=color)
ax.set_xlabel("AF Threshold (TP with AF < threshold)")
ax.set_ylabel("LOH-like Fraction among TP")
ax.set_title("Q2: LOH-like Fraction Among Low-AF TP\n(Higher = more rescue potential)",
             fontweight="bold")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig09_low_af_tp_loh_rescue.png")
plt.close(fig)
print("    ✓ fig09_low_af_tp_loh_rescue.png")

# --- Fig 10: HP0 ratio by tier (corrected data, paired vs TO) ---
print("  fig10: HP0 ratio by tier violin")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for idx, mode in enumerate(["paired", "to"]):
    ax = axes[idx]
    sub = df[df["mode"] == mode].copy()
    sub["tier_cat"] = pd.Categorical(sub["tier"], categories=TIER_ORDER, ordered=True)
    sub_sample = sub.sample(min(60000, len(sub)), random_state=42)
    sns.violinplot(data=sub_sample, x="tier_cat", y="hp0_ratio",
                   palette=TIER_COLORS, ax=ax, cut=0, inner="quartile", scale="width")
    ax.set_xlabel("Support Tier")
    ax.set_ylabel("HP0 Ratio")
    ax.set_title(f"{mode.upper()}", fontweight="bold")
    ax.set_ylim(-0.02, 1.02)
fig.suptitle("HP0 Ratio by Support Tier (Corrected Data)\n(Lower HP0 = better phasing quality)",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig10_hp0_by_tier_violin.png")
plt.close(fig)
print("    ✓ fig10_hp0_by_tier_violin.png")

# ============================================================
# Summary Report
# ============================================================
print("\n[6/7] Writing context and summary ...")

context = {
    "workspace": str(OUT_DIR),
    "generated_at": GENERATED_AT,
    "master_dataset": str(ROUND1_WS / "all_region_rows.tsv.gz"),
    "old_master_dataset": str(OLD_ROUND1_WS / "all_region_rows.tsv.gz"),
    "total_rows": int(len(df)),
    "to_rows": int((df["mode"] == "to").sum()),
    "paired_rows": int((df["mode"] == "paired").sum()),
    "global_to_enrichment": 0.805,
    "global_paired_enrichment": 1.194,
    "questions_investigated": ["Q1: opposite enrichment directions", "Q2: TP rescue", "Q3: tier structure"],
    "figures_generated": 10,
    "tsv_generated": len(list(OUT_DIR.glob("*.tsv"))),
    "verification_schema": verification_schema_metadata,
}
(OUT_DIR / "investigation_context.json").write_text(json.dumps(context, indent=2, ensure_ascii=False) + "\n")

# Key findings summary
summary_lines = [
    "# Post-HP-Fix TO LOH Deep Investigation — Key Findings",
    f"\n生成時間: {GENERATED_AT}",
    f"資料來源: {ROUND1_WS / 'all_region_rows.tsv.gz'}",
    "",
    "## Q1: 為什麼 TO 和 paired LOH enrichment 方向相反？",
    "",
]

# Compute Q1 answer from data
for mode in ["paired", "to"]:
    tp_data = q1e_df[(q1e_df["mode"] == mode) & (q1e_df["truth"] == "TP")]
    fp_data = q1e_df[(q1e_df["mode"] == mode) & (q1e_df["truth"] == "FP")]
    if len(tp_data) > 0 and len(fp_data) > 0:
        summary_lines.append(f"### {mode.upper()}")
        summary_lines.append(f"- TP extreme LOH: {tp_data['pct_extreme_loh'].values[0]:.1%}")
        summary_lines.append(f"- FP extreme LOH: {fp_data['pct_extreme_loh'].values[0]:.1%}")
        summary_lines.append(f"- TP balanced: {tp_data['pct_balanced'].values[0]:.1%}")
        summary_lines.append(f"- FP balanced: {fp_data['pct_balanced'].values[0]:.1%}")
        summary_lines.append(f"- TP AF median: {q1a_df[(q1a_df['mode']==mode)&(q1a_df['truth']=='TP')&(q1a_df['loh_like']==False)]['af_median'].values[0]:.3f}")
        summary_lines.append(f"- FP AF median: {q1a_df[(q1a_df['mode']==mode)&(q1a_df['truth']=='FP')&(q1a_df['loh_like']==False)]['af_median'].values[0]:.3f}")
        summary_lines.append("")

summary_lines.extend([
    "## Q2: TO LOH 作為 TP rescue 候選",
    "",
    "### Low-AF TP 的 LOH-like 比例",
])
for _, row in q2c_df[q2c_df["mode"] == "to"].iterrows():
    summary_lines.append(f"- AF < {row['af_threshold']}: {row['tp_loh_frac']:.1%} LOH-like (n={int(row['tp_count']):,})")
summary_lines.append("")

summary_lines.extend([
    "## Q3: Tier 結構重塑",
    "",
    "### TO Tier 分佈 (before → after fix)",
])
for tier in TIER_ORDER:
    before = q3a_df[(q3a_df["dataset"] == "before_fix") & (q3a_df["mode"] == "to") & (q3a_df["tier"] == tier)]
    after = q3a_df[(q3a_df["dataset"] == "after_fix") & (q3a_df["mode"] == "to") & (q3a_df["tier"] == tier)]
    if len(before) > 0 and len(after) > 0:
        summary_lines.append(
            f"- Tier {tier}: {before['fraction'].values[0]:.1%} → {after['fraction'].values[0]:.1%} "
            f"(n: {int(before['n'].values[0]):,} → {int(after['n'].values[0]):,})"
        )

(OUT_DIR / "investigation_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
print(f"  investigation_summary.md written")

print(f"\n[7/7] ✅ Investigation complete! Output: {OUT_DIR}")
print(f"  TSV files: {len(list(OUT_DIR.glob('*.tsv')))}")
print(f"  Figures: {len(list(FIG_DIR.glob('*.png')))}")
