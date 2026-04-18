#!/usr/bin/env python3
"""Step 5: Amplicon / High-CN Blacklist Pilot (Method A).

Context: Step 4 revealed HCC1954 Z3 FPs are NOT confined to literature focal
    amplicons (HER2 chr17:35-42M captures 0 FPs) — they are arm-level LOH:
      - chr8p (0-40Mb): 68% of chr8 FPs
      - chr17p (0-20Mb): 53% of chr17 FPs
      - chr5: distributed (HCC1954 has whole-chromosome rearrangement)

    This reflects classic HER2+ breast cancer CN architecture:
      chr8p loss + chr17p TP53 LOH + HER2 focal gain + chr5 rearrangement.

    The pilot compares three blacklist strategies on ALL 7 samples:
      S1 (literature arm-level): chr8p + chr17p + chr5 + focal HER2/MYC
      S2 (whole-chr): chr5/chr8/chr17 entirely (most aggressive)
      S3 (data-driven CovM): sample-specific Coverage_Multiple > 95th %ile

    F1 impact is computed as:
      new_TP = TP not in blacklist
      new_FP = FP not in blacklist
      new_FN = FN_original + TP_in_blacklist  (filter rejects TPs too)
      precision = new_TP / (new_TP + new_FP)
      recall    = new_TP / (new_TP + new_FN)
      F1        = 2*P*R/(P+R)

    Circular-dependency guard:
      - S1 uses literature coords independent of our data
      - S3 cutoff is computed within each sample's non-Z3 distribution
        (does not peek at Z3 TP/FP labels)
"""
import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

MASTER = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
          "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/z3_internal_feature_exploration"
FIG_DIR = os.path.join(OUT_DIR, "figures")
DATA_DIR = os.path.join(OUT_DIR, "data")

SAMPLE_ORDER = ["H2009", "H1437", "HCC1395", "HCC1395_DORADO",
                "HCC1937", "HCC1954", "COLO829"]

# ── Blacklist S1: literature arm-level (HCC1954-informed) ──────────
# Reference: refLinker Hi-C 2024; cansar.ICR; CGH array studies of HCC1954
#   chr8p loss + MYC amp at 8q24
#   chr17p TP53 LOH + HER2 focal at 17q12
#   chr5 rearrangement (full chromosome LOH/gain)
S1_LITERATURE_ARMS = [
    ("chr5",   1_000_000, 180_000_000, "chr5 (HCC1954 full rearrangement)"),
    ("chr8",           0,  45_000_000, "chr8p loss"),
    ("chr8",  120_000_000, 145_000_000, "chr8q MYC amplicon"),
    ("chr17",          0,  25_000_000, "chr17p TP53 LOH"),
    ("chr17", 35_000_000,  42_000_000, "chr17q12 HER2 focal"),
]

# ── Blacklist S2: whole-chr (most aggressive, HCC1954-targeted) ────
S2_WHOLE_CHR = [
    ("chr5",  0, 300_000_000, "chr5 whole"),
    ("chr8",  0, 300_000_000, "chr8 whole"),
    ("chr17", 0, 300_000_000, "chr17 whole"),
]


def in_blacklist(df, coords):
    """Return boolean mask of rows falling in any blacklist interval."""
    mask = np.zeros(len(df), dtype=bool)
    for chrom, start, end, _ in coords:
        hit = (df["Chr"] == chrom) & (df["Pos"] >= start) & (df["Pos"] <= end)
        mask |= hit.values
    return mask


def compute_f1(tp, fp, fn):
    if tp + fp == 0 or tp + fn == 0:
        return 0.0, 0.0, 0.0
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f1


# ── Load data ──────────────────────────────────────────────────────
COLS = ["sample", "mode", "truth_label", "Chr", "Pos",
        "core_loh_like", "to_loh_bed_hit",
        "caller_af", "HPFineNGroups", "Coverage_Multiple"]

print("Loading master dataset...")
df = pd.read_csv(MASTER, sep='\t', usecols=COLS, low_memory=False)

for c in ["caller_af", "HPFineNGroups", "Coverage_Multiple", "Pos"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
for c in ["core_loh_like", "to_loh_bed_hit"]:
    df[c] = df[c].map({"True": True, "False": False,
                       True: True, False: False}).fillna(False)

df["is_tp"] = df["truth_label"] == "TP"
df["is_fn"] = df["truth_label"] == "FN"
df["loh_flag"] = np.where(df["mode"] == "paired",
                          df["core_loh_like"], df["to_loh_bed_hit"])

df_to = df[df["mode"] == "to"].copy()
df_to["af_extreme"] = df_to["caller_af"].apply(
    lambda x: (not pd.isna(x)) and (x < 0.1 or x > 0.9))
df_to["in_z3"] = (df_to["loh_flag"] & df_to["af_extreme"]
                  & (df_to["HPFineNGroups"] <= 1))

# Sample-specific CovM 95%ile cutoff (computed on non-Z3 as baseline to
# avoid circular dependency with Z3 TP/FP labels)
covm_cutoff = {}
for sample in SAMPLE_ORDER:
    baseline = df_to[(df_to["sample"] == sample) & (~df_to["in_z3"])
                     & df_to["Coverage_Multiple"].notna()]
    covm_cutoff[sample] = baseline["Coverage_Multiple"].quantile(0.95)

print("\nSample-specific CovM 95th percentile cutoffs (non-Z3 baseline):")
for s, v in covm_cutoff.items():
    print(f"  {s:16s}: CovM > {v:.3f}")


# ── Run three strategies and compute before/after F1 ───────────────
def run_strategy(df_sample, coords=None, use_covm=False, covm_thresh=None,
                 z3_only=True):
    """Apply blacklist (optionally only to in_z3 rows) and return metrics.

    z3_only=True: filter only rejects calls that are BOTH in_z3 AND in blacklist.
    This preserves high-confidence non-Z3 TPs (the core correction vs naive
    whole-chromosome blacklist).
    """
    tp_total = int(df_sample["is_tp"].sum())
    fp_total = int((df_sample["truth_label"] == "FP").sum())
    fn_total = int(df_sample["is_fn"].sum())

    if coords is not None:
        bl_mask = in_blacklist(df_sample, coords)
    elif use_covm:
        bl_mask = (df_sample["Coverage_Multiple"] > covm_thresh).fillna(False).values
    else:
        bl_mask = np.zeros(len(df_sample), dtype=bool)

    if z3_only:
        mask = bl_mask & df_sample["in_z3"].values
    else:
        mask = bl_mask

    # Before
    p0, r0, f0 = compute_f1(tp_total, fp_total, fn_total)

    tp_kept = int((df_sample["is_tp"] & (~mask)).sum())
    fp_kept = int(((df_sample["truth_label"] == "FP") & (~mask)).sum())
    tp_lost = int((df_sample["is_tp"] & mask).sum())
    fp_removed = int(((df_sample["truth_label"] == "FP") & mask).sum())
    fn_new = fn_total + tp_lost

    p1, r1, f1_new = compute_f1(tp_kept, fp_kept, fn_new)
    return {
        "tp_total": tp_total, "fp_total": fp_total, "fn_total": fn_total,
        "f1_before": f0, "p_before": p0, "r_before": r0,
        "f1_after": f1_new, "p_after": p1, "r_after": r1,
        "f1_delta": f1_new - f0,
        "fp_removed": fp_removed, "tp_lost": tp_lost,
    }


print("\n=== Strategy 1: Literature arm-level ∩ Z3 ===")
s1_rows = []
for sample in SAMPLE_ORDER:
    ds = df_to[df_to["sample"] == sample]
    r = run_strategy(ds, coords=S1_LITERATURE_ARMS, z3_only=True)
    r["sample"] = sample
    s1_rows.append(r)
    print(f"  {sample:16s}: F1 {r['f1_before']:.4f} → {r['f1_after']:.4f} "
          f"(Δ={r['f1_delta']:+.4f}), FP-={r['fp_removed']:5d}, TP-={r['tp_lost']:4d}")

print("\n=== Strategy 2: Whole-chr5/8/17 ∩ Z3 ===")
s2_rows = []
for sample in SAMPLE_ORDER:
    ds = df_to[df_to["sample"] == sample]
    r = run_strategy(ds, coords=S2_WHOLE_CHR, z3_only=True)
    r["sample"] = sample
    s2_rows.append(r)
    print(f"  {sample:16s}: F1 {r['f1_before']:.4f} → {r['f1_after']:.4f} "
          f"(Δ={r['f1_delta']:+.4f}), FP-={r['fp_removed']:5d}, TP-={r['tp_lost']:4d}")

print("\n=== Strategy 3: CovM > 95%ile (non-Z3 baseline) ∩ Z3 ===")
s3_rows = []
for sample in SAMPLE_ORDER:
    ds = df_to[df_to["sample"] == sample]
    r = run_strategy(ds, use_covm=True, covm_thresh=covm_cutoff[sample],
                     z3_only=True)
    r["sample"] = sample
    r["covm_cutoff"] = covm_cutoff[sample]
    s3_rows.append(r)
    print(f"  {sample:16s}: CovM>{covm_cutoff[sample]:.2f}, "
          f"F1 {r['f1_before']:.4f} → {r['f1_after']:.4f} "
          f"(Δ={r['f1_delta']:+.4f}), FP-={r['fp_removed']:5d}, TP-={r['tp_lost']:4d}")

# ── Strategy 4: F pilot canonical (sanity baseline): NG=1 in Z3 reject ──
print("\n=== Strategy 4: Reject all Z3 (naive upper-bound, sanity baseline) ===")
s4_rows = []
for sample in SAMPLE_ORDER:
    ds = df_to[df_to["sample"] == sample]
    r = run_strategy(ds, coords=[("chr1", 0, 0, "")], z3_only=False)  # dummy
    # Override: reject everything in_z3
    tp_total = int(ds["is_tp"].sum())
    fp_total = int((ds["truth_label"] == "FP").sum())
    fn_total = int(ds["is_fn"].sum())
    mask = ds["in_z3"].values
    tp_kept = int((ds["is_tp"] & (~mask)).sum())
    fp_kept = int(((ds["truth_label"] == "FP") & (~mask)).sum())
    tp_lost = int((ds["is_tp"] & mask).sum())
    fp_removed = int(((ds["truth_label"] == "FP") & mask).sum())
    fn_new = fn_total + tp_lost
    p0, r0, f0 = compute_f1(tp_total, fp_total, fn_total)
    p1, rr, f1_new = compute_f1(tp_kept, fp_kept, fn_new)
    r4 = {
        "sample": sample, "tp_total": tp_total, "fp_total": fp_total,
        "fn_total": fn_total, "f1_before": f0, "f1_after": f1_new,
        "f1_delta": f1_new - f0, "fp_removed": fp_removed, "tp_lost": tp_lost,
        "p_before": p0, "r_before": r0, "p_after": p1, "r_after": rr,
    }
    s4_rows.append(r4)
    print(f"  {sample:16s}: F1 {f0:.4f} → {f1_new:.4f} (Δ={f1_new-f0:+.4f}), "
          f"FP-={fp_removed:5d}, TP-={tp_lost:4d}")

# ── Save results ───────────────────────────────────────────────────
all_rows = []
for row_list, strat in [(s1_rows, "S1_literature_capZ3"),
                         (s2_rows, "S2_whole_chr_capZ3"),
                         (s3_rows, "S3_covm_95pct_capZ3"),
                         (s4_rows, "S4_reject_all_Z3")]:
    for r in row_list:
        r2 = dict(r)
        r2["strategy"] = strat
        all_rows.append(r2)

res_df = pd.DataFrame(all_rows)
res_df.to_csv(os.path.join(DATA_DIR, "step5_blacklist_pilot_results.tsv"),
              sep='\t', index=False)

# Also save the blacklist bed files
bed_rows = []
for bed_name, coords in [("S1_literature", S1_LITERATURE_ARMS),
                          ("S2_whole_chr", S2_WHOLE_CHR)]:
    for chrom, start, end, label in coords:
        bed_rows.append({"strategy": bed_name, "chr": chrom,
                         "start": start, "end": end, "label": label})
bed_df = pd.DataFrame(bed_rows)
bed_df.to_csv(os.path.join(DATA_DIR, "step5_blacklist_bed.tsv"),
              sep='\t', index=False)

# ── Figure: Delta F1 per sample per strategy ──────────────────────
fig, axes = plt.subplots(1, 4, figsize=(19, 5), sharey=True)
strategies = [("S1 Literature Arms ∩ Z3", s1_rows),
              ("S2 Whole chr5/8/17 ∩ Z3", s2_rows),
              ("S3 CovM >95%ile ∩ Z3", s3_rows),
              ("S4 Reject all Z3 (ceiling)", s4_rows)]

for ax, (name, rows) in zip(axes, strategies):
    samples = [r["sample"] for r in rows]
    deltas = [r["f1_delta"] for r in rows]
    colors = ['#e74c3c' if d < 0 else '#2ecc71' for d in deltas]
    bars = ax.bar(range(len(samples)), deltas, color=colors)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel("ΔF1 (after − before)")
    ax.set_title(name, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    for bar, d in zip(bars, deltas):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (0.001 if d >= 0 else -0.003),
                f"{d:+.4f}",
                ha='center', va='bottom' if d >= 0 else 'top',
                fontsize=7)

fig.suptitle("Amplicon Blacklist Pilot — ΔF1 per sample per strategy (TO mode)",
             fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "step5_blacklist_delta_f1.png"),
            dpi=150, bbox_inches='tight')
plt.close()

# ── Summary Table ──────────────────────────────────────────────────
print("\n" + "=" * 80)
print("SUMMARY — Blacklist Pilot")
print("=" * 80)

summary = []
for strat_name, rows in [("S1_literature_capZ3", s1_rows),
                          ("S2_whole_chr_capZ3", s2_rows),
                          ("S3_covm_95pct_capZ3", s3_rows),
                          ("S4_reject_all_Z3", s4_rows)]:
    hcc1954 = [r for r in rows if r["sample"] == "HCC1954"][0]
    others = [r for r in rows if r["sample"] != "HCC1954"]
    other_mean_delta = np.mean([r["f1_delta"] for r in others])
    n_improved = sum(1 for r in others if r["f1_delta"] > 0)
    n_hurt = sum(1 for r in others if r["f1_delta"] < 0)
    summary.append({
        "strategy": strat_name,
        "HCC1954_dF1": hcc1954["f1_delta"],
        "HCC1954_FP_removed": hcc1954["fp_removed"],
        "HCC1954_TP_lost": hcc1954["tp_lost"],
        "other_6_mean_dF1": other_mean_delta,
        "other_improved_n": n_improved,
        "other_hurt_n": n_hurt,
    })
    print(f"\n{strat_name}:")
    print(f"  HCC1954: ΔF1={hcc1954['f1_delta']:+.4f}, "
          f"FP-={hcc1954['fp_removed']}, TP-={hcc1954['tp_lost']}")
    print(f"  Other 6: mean ΔF1={other_mean_delta:+.4f}, "
          f"improved={n_improved}/6, hurt={n_hurt}/6")

sum_df = pd.DataFrame(summary)
sum_df.to_csv(os.path.join(DATA_DIR, "step5_blacklist_pilot_summary.tsv"),
              sep='\t', index=False)

print("\nOutputs:")
print(f"  {DATA_DIR}/step5_blacklist_pilot_results.tsv")
print(f"  {DATA_DIR}/step5_blacklist_bed.tsv")
print(f"  {DATA_DIR}/step5_blacklist_pilot_summary.tsv")
print(f"  {FIG_DIR}/step5_blacklist_delta_f1.png")
