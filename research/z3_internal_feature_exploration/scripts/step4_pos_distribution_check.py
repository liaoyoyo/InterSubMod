#!/usr/bin/env python3
"""Step 4: Validate HCC1954 Z3 FP Pos distribution matches known amplicons.

Context: Step 3 showed Z3 FP concentrates on chr5/8/17 (93%). Before
    applying a literature-informed amplicon blacklist, verify the FP
    positions cluster within known HER2 (chr17q12), MYC (chr8q24), and
    chr5p amplicon coordinates reported in HCC1954 literature.

Outputs:
  - figures/step4_hcc1954_fp_pos_hist.png — per-chr Pos histograms
  - data/step4_amplicon_literature_coords.tsv — literature amplicon bed
  - data/step4_hcc1954_fp_amplicon_hit.tsv — per-chr FP enrichment in amplicon vs outside
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

# Literature-informed HCC1954 amplicons (hg38)
#   HER2/ERBB2 amplicon: chr17q12 gene at 39.7-39.8Mb; typical focal amp 37-42Mb
#   MYC amplicon: chr8q24 at 127.7Mb; typical amp region 120-140Mb
#   chr5p gain: common 1-40Mb in HER2+ breast
# References:
#   - refLinker Hi-C 2024 (biorxiv 10.1101/2024.03.02.583108)
#   - Cellosaurus CVCL_1259
#   - cansar ICR HCC1954 CNV entries
AMPLICON_COORDS = [
    ("chr5",  1_000_000,  50_000_000, "5p gain"),
    ("chr8",  120_000_000, 140_000_000, "MYC amplicon (8q24)"),
    ("chr17", 35_000_000,  42_000_000, "HER2 amplicon (17q12)"),
]

COLS = ["sample", "mode", "truth_label", "Chr", "Pos",
        "core_loh_like", "to_loh_bed_hit",
        "caller_af", "HPFineNGroups", "Coverage_Multiple"]

print("Loading master dataset...")
df = pd.read_csv(MASTER, sep='\t', usecols=COLS, low_memory=False)

df["caller_af"] = pd.to_numeric(df["caller_af"], errors="coerce")
df["HPFineNGroups"] = pd.to_numeric(df["HPFineNGroups"], errors="coerce")
df["Coverage_Multiple"] = pd.to_numeric(df["Coverage_Multiple"], errors="coerce")
df["Pos"] = pd.to_numeric(df["Pos"], errors="coerce")
for col in ["core_loh_like", "to_loh_bed_hit"]:
    df[col] = df[col].map({"True": True, "False": False,
                           True: True, False: False}).fillna(False)

df["is_tp"] = df["truth_label"] == "TP"
df["loh_flag"] = np.where(df["mode"] == "paired",
                          df["core_loh_like"], df["to_loh_bed_hit"])

df_to = df[df["mode"] == "to"].copy()
df_to["af_extreme"] = df_to["caller_af"].apply(
    lambda x: (not pd.isna(x)) and (x < 0.1 or x > 0.9))
df_to["in_z3"] = (df_to["loh_flag"] & df_to["af_extreme"]
                  & (df_to["HPFineNGroups"] <= 1))

# Save amplicon bed
amplicon_df = pd.DataFrame(AMPLICON_COORDS,
                           columns=["chr", "start", "end", "label"])
amplicon_df.to_csv(os.path.join(DATA_DIR, "step4_amplicon_literature_coords.tsv"),
                   sep='\t', index=False)

# ── HCC1954 Z3 per-chr Pos distribution ─────────────────────────────
z3_1954 = df_to[(df_to["sample"] == "HCC1954") & df_to["in_z3"]].copy()

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, chrom in zip(axes, ["chr5", "chr8", "chr17"]):
    sub = z3_1954[z3_1954["Chr"] == chrom]
    tp_pos = sub[sub["is_tp"]]["Pos"].values / 1e6
    fp_pos = sub[~sub["is_tp"]]["Pos"].values / 1e6

    # Determine axis range
    all_pos = np.concatenate([tp_pos, fp_pos]) if len(tp_pos) + len(fp_pos) else [0, 250]
    ax.hist(fp_pos, bins=60, alpha=0.6, color='#e74c3c',
            label=f"FP (n={len(fp_pos)})")
    if len(tp_pos):
        ax.hist(tp_pos, bins=60, alpha=0.6, color='#2ecc71',
                label=f"TP (n={len(tp_pos)})")

    # Highlight literature amplicon region
    for ch, start, end, lab in AMPLICON_COORDS:
        if ch == chrom:
            ax.axvspan(start / 1e6, end / 1e6, alpha=0.2,
                       color='#f39c12', label=f"lit. {lab}")
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel("count")
    ax.set_title(f"HCC1954 Z3 — {chrom}")
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)

fig.suptitle("HCC1954 Z3 FP position distribution vs literature amplicon coords",
             fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "step4_hcc1954_fp_pos_hist.png"),
            dpi=150, bbox_inches='tight')
plt.close()

# ── Compute hit rate: how many HCC1954 Z3 FPs fall in literature amplicon ──
rows = []
for chrom, start, end, label in AMPLICON_COORDS:
    sub = z3_1954[z3_1954["Chr"] == chrom]
    in_amp = (sub["Pos"] >= start) & (sub["Pos"] <= end)
    out_amp = ~in_amp

    tp_in = (in_amp & sub["is_tp"]).sum()
    fp_in = (in_amp & (~sub["is_tp"])).sum()
    tp_out = (out_amp & sub["is_tp"]).sum()
    fp_out = (out_amp & (~sub["is_tp"])).sum()

    rows.append({
        "chr": chrom, "start": start, "end": end, "label": label,
        "n_total_on_chr": len(sub),
        "n_in_amplicon": int(in_amp.sum()),
        "tp_in_amp": int(tp_in), "fp_in_amp": int(fp_in),
        "tp_out_amp": int(tp_out), "fp_out_amp": int(fp_out),
        "tp_rate_in_amp": tp_in / (tp_in + fp_in) if (tp_in + fp_in) else np.nan,
        "tp_rate_out_amp": tp_out / (tp_out + fp_out) if (tp_out + fp_out) else np.nan,
        "pct_fp_in_amp": fp_in / (fp_in + fp_out) if (fp_in + fp_out) else np.nan,
    })

hit_df = pd.DataFrame(rows)
hit_df.to_csv(os.path.join(DATA_DIR, "step4_hcc1954_fp_amplicon_hit.tsv"),
              sep='\t', index=False)

print("\n=== HCC1954 Z3 FP hit rate in literature amplicon coords ===")
print(hit_df[["chr", "label", "n_total_on_chr", "n_in_amplicon",
              "tp_rate_in_amp", "tp_rate_out_amp", "pct_fp_in_amp"]].to_string(index=False))

# Sum all Z3 FP captured by literature blacklist
total_fp_1954 = int((~z3_1954["is_tp"]).sum())
total_tp_1954 = int(z3_1954["is_tp"].sum())

in_blacklist = np.zeros(len(z3_1954), dtype=bool)
for chrom, start, end, _ in AMPLICON_COORDS:
    in_blacklist |= ((z3_1954["Chr"] == chrom)
                     & (z3_1954["Pos"] >= start)
                     & (z3_1954["Pos"] <= end))

fp_captured = int((in_blacklist & (~z3_1954["is_tp"])).sum())
tp_captured = int((in_blacklist & z3_1954["is_tp"]).sum())

print(f"\nHCC1954 Z3 totals:  TP={total_tp_1954}, FP={total_fp_1954}")
print(f"Literature blacklist captures: TP={tp_captured} ({tp_captured/max(total_tp_1954,1)*100:.1f}%), "
      f"FP={fp_captured} ({fp_captured/max(total_fp_1954,1)*100:.1f}%)")
print(f"Blacklist efficiency:  FP%/TP% = "
      f"{(fp_captured/max(total_fp_1954,1)) / max(tp_captured/max(total_tp_1954,1), 1e-9):.2f}")
