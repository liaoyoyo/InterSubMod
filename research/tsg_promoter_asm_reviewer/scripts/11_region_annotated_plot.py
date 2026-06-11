#!/usr/bin/env python3
"""Region-annotated per-CpG Δβ plot showing where hypomethylation concentrates."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
plt.rcParams["font.family"] = ["DejaVu Sans"]

BRCA2_MANE_TSS = 32315508
BRCA2_ALT_TSS = 32315086
ZAR1L_TSS = 32315363
CPG_ISLAND = (32315396, 32315763)
VARIANT = 32315128

with open(f"{BASE}/output/step4_ism_results.json") as f:
    results = json.load(f)
brca2 = next(r for r in results if r["candidate"]["gene"] == "BRCA2")

xs, deltas, b_g, b_s = [], [], [], []
for r in brca2["cpg_records"]:
    if r["HP1_n"] >= 3 and r["HP1-1_n"] >= 3 and r["HP1_beta"] is not None and r["HP1-1_beta"] is not None:
        pos = VARIANT + r["dist_to_var"]
        xs.append(pos)
        deltas.append(r["HP1-1_beta"] - r["HP1_beta"])
        b_g.append(r["HP1_beta"])
        b_s.append(r["HP1-1_beta"])
xs = np.array(xs); deltas = np.array(deltas); b_g = np.array(b_g); b_s = np.array(b_s)

fig, axes = plt.subplots(2, 1, figsize=(15, 9), sharex=True, gridspec_kw={"height_ratios":[3,2]})

# Panel 1: germline vs somatic beta
ax = axes[0]
ax.scatter(xs, b_g, s=22, c="#2563eb", label="HP1 germline β", alpha=0.7, zorder=3)
ax.scatter(xs, b_s, s=22, c="#dc2626", label="HP1-1 somatic-reconstructed β", alpha=0.7, zorder=3)
# Functional region shading
ax.axvspan(CPG_ISLAND[0], CPG_ISLAND[1], color="green", alpha=0.12, zorder=0, label="CpG island core")
ax.axvspan(CPG_ISLAND[0]-2000, CPG_ISLAND[0], color="orange", alpha=0.06, zorder=0, label="CpG island shore (upstream)")
# Landmarks
ax.axvline(VARIANT, color="purple", lw=2, ls="--", zorder=4, label=f"somatic SNV {VARIANT}")
ax.axvline(BRCA2_MANE_TSS, color="#16a34a", lw=1.5, ls=":", zorder=4, label=f"BRCA2 MANE TSS")
ax.axvline(BRCA2_ALT_TSS, color="#16a34a", lw=1, ls="-.", zorder=4, alpha=0.6, label="BRCA2 alt TSS")
ax.axvline(ZAR1L_TSS, color="#b45309", lw=1.5, ls=":", zorder=4, label="ZAR1L TSS (- strand)")
ax.set_ylabel("methylation β (5mCG)", fontsize=11)
ax.set_ylim(-0.05, 1.05)
ax.set_title("BRCA2/ZAR1L bidirectional promoter — per-CpG methylation by haplotype\n"
             "Hypomethylation concentrates in CpG-island SHORE (upstream −280 to −870bp), NOT in island core (constitutively unmethylated)",
             fontsize=11)
ax.legend(loc="center right", fontsize=8, ncol=1)
ax.grid(alpha=0.3)

# Panel 2: delta beta bars colored by region
ax = axes[1]
colors = []
for p, d in zip(xs, deltas):
    if CPG_ISLAND[0] <= p <= CPG_ISLAND[1]:
        colors.append("#16a34a")   # island core
    elif p < CPG_ISLAND[0]:
        colors.append("#dc2626")   # upstream shore
    else:
        colors.append("#9ca3af")   # downstream
ax.bar(xs, deltas, width=10, color=colors, alpha=0.7)
ax.axhline(0, color="black", lw=0.7)
ax.axhline(-0.2, color="red", lw=0.6, ls=":", label="Lister |Δβ|=0.2 ribbon")
ax.axvspan(CPG_ISLAND[0], CPG_ISLAND[1], color="green", alpha=0.12, zorder=0)
ax.axvline(VARIANT, color="purple", lw=2, ls="--")
ax.axvline(BRCA2_MANE_TSS, color="#16a34a", lw=1.5, ls=":")
ax.axvline(ZAR1L_TSS, color="#b45309", lw=1.5, ls=":")
ax.set_ylabel("Δβ (somatic − germline)", fontsize=11)
ax.set_xlabel("chr13 position (hg38)", fontsize=11)
ax.set_ylim(-1.0, 0.5)
ax.legend(loc="lower left", fontsize=9)
ax.grid(alpha=0.3)
# Annotate region means
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color="#dc2626", label="upstream shore (hypomethyl)")
green_patch = mpatches.Patch(color="#16a34a", label="CpG island core (no diff)")
grey_patch = mpatches.Patch(color="#9ca3af", label="downstream 5'UTR")
ax.legend(handles=[red_patch, green_patch, grey_patch], loc="lower right", fontsize=8)

plt.tight_layout()
plt.savefig(f"{BASE}/figures/brca2_region_annotated_delta.png", dpi=130, bbox_inches="tight")
print("Saved figures/brca2_region_annotated_delta.png")

# Quantitative summary for report
shore = [(p,d) for p,d in zip(xs,deltas) if p < CPG_ISLAND[0]]
core = [(p,d) for p,d in zip(xs,deltas) if CPG_ISLAND[0] <= p <= CPG_ISLAND[1]]
down = [(p,d) for p,d in zip(xs,deltas) if p > CPG_ISLAND[1]]
print(f"\nUpstream shore (pos < {CPG_ISLAND[0]}): n={len(shore)}, mean Δβ={np.mean([d for _,d in shore]):+.4f}")
print(f"CpG island core ({CPG_ISLAND[0]}-{CPG_ISLAND[1]}): n={len(core)}, mean Δβ={np.mean([d for _,d in core]):+.4f}")
print(f"Downstream 5'UTR (pos > {CPG_ISLAND[1]}): n={len(down)}, mean Δβ={np.mean([d for _,d in down]):+.4f}")
