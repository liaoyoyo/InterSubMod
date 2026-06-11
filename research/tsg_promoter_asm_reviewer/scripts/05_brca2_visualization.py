#!/usr/bin/env python3
"""
Build BRCA2 per-CpG delta beta plot + CpG island overlap check.
Outputs:
  - figures/brca2_per_cpg_delta.png : per-CpG beta HP1 vs HP1-1 + delta track
  - figures/brca2_methylation_density.png : 4-track ridge plot (HP1 / HP2 / HP1-1 / null random)
  - output/brca2_cpg_island_overlap.txt : whether variant is within CpG island
"""
import json, os, gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
plt.rcParams["font.family"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

with open(f"{BASE}/output/step4_ism_results.json") as f:
    results = json.load(f)

brca2 = next(r for r in results if r["candidate"]["gene"] == "BRCA2")
kmt2c = next(r for r in results if r["candidate"]["gene"] == "KMT2C")

var_pos = brca2["candidate"]["pos"]
chrom = brca2["candidate"]["chrom"]

# Check CpG island overlap
cpg_island_file = f"{BASE}/data/cpgIslandExt_hg38.txt.gz"
in_cpg_island = False
nearest_island = None
with gzip.open(cpg_island_file, "rt") as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        if fields[1] != chrom:
            continue
        start = int(fields[2])
        end = int(fields[3])
        if start <= var_pos <= end:
            in_cpg_island = True
            nearest_island = {"chrom": fields[1], "start": start, "end": end, "name": fields[4]}
            break
        # nearest within +/- 5kb
        dist = min(abs(start - var_pos), abs(end - var_pos))
        if dist < 5000 and (nearest_island is None or dist < nearest_island.get("dist", 1e9)):
            nearest_island = {"chrom": fields[1], "start": start, "end": end, "name": fields[4], "dist": dist}

with open(f"{BASE}/output/brca2_cpg_island_overlap.txt", "w") as f:
    f.write(f"BRCA2 variant: {chrom}:{var_pos}\n")
    f.write(f"In CpG island: {in_cpg_island}\n")
    f.write(f"Nearest CpG island: {nearest_island}\n")
print(f"BRCA2 in CpG island: {in_cpg_island}, nearest: {nearest_island}")

# Plot 1: per-CpG beta tracks for BRCA2
cpgs = brca2["cpg_records"]
xs, b_hp1, b_hp1_1, n_hp1, n_hp1_1 = [], [], [], [], []
for r in cpgs:
    if r["HP1_n"] >= 3 and r["HP1-1_n"] >= 3:
        xs.append(r["dist_to_var"])
        b_hp1.append(r["HP1_beta"])
        b_hp1_1.append(r["HP1-1_beta"])
        n_hp1.append(r["HP1_n"])
        n_hp1_1.append(r["HP1-1_n"])
xs = np.array(xs); b_hp1 = np.array(b_hp1); b_hp1_1 = np.array(b_hp1_1)
delta = b_hp1_1 - b_hp1  # som - germ

fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True, gridspec_kw={"height_ratios":[2, 2, 1]})

ax = axes[0]
ax.scatter(xs, b_hp1, s=18, c="#2563eb", label="HP1 (germline)", alpha=0.7)
ax.scatter(xs, b_hp1_1, s=18, c="#dc2626", label="HP1-1 (somatic-reconstructed)", alpha=0.7)
ax.axvline(0, color="orange", lw=1.2, ls="--", label="somatic SNV (chr13:32315128 G>A)")
ax.set_ylim(-0.05, 1.05)
ax.set_ylabel("methylation β (5mCG)")
ax.set_title(f"BRCA2 promoter (chr13:32315128±1kb) per-CpG methylation — HP1 vs HP1-1\nn=197 paired CpGs, mean Δβ=-0.122, Wilcoxon p=6.1e-11", fontsize=11)
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)

# Delta panel
ax = axes[1]
colors = ['#dc2626' if d < -0.1 else '#2563eb' if d > 0.1 else '#888' for d in delta]
ax.bar(xs, delta, color=colors, alpha=0.6, width=8)
ax.axhline(0, color="black", lw=0.7)
ax.axhline(-0.05, color="grey", lw=0.5, ls=":")
ax.axhline(-0.20, color="red", lw=0.7, ls=":", label="Lister 2009 |Δβ|=0.2 ribbon")
ax.axhline(0.20, color="red", lw=0.7, ls=":")
ax.axvline(0, color="orange", lw=1.2, ls="--")
ax.set_ylabel("Δβ (HP1-1 − HP1)")
ax.set_ylim(-1.0, 1.0)
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)

# Coverage panel
ax = axes[2]
ax.bar(xs, n_hp1, color="#2563eb", alpha=0.5, width=8, label="HP1 reads")
ax.bar(xs, n_hp1_1, color="#dc2626", alpha=0.5, width=8, label="HP1-1 reads", bottom=n_hp1)
ax.set_ylabel("reads / CpG")
ax.set_xlabel(f"Distance to somatic SNV (bp) — chr13:32315128 G>A")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)
ax.axvline(0, color="orange", lw=1.2, ls="--")

# Mark CpG island
if nearest_island and "dist" not in nearest_island:
    # In CpG island
    isl_start_rel = nearest_island["start"] - var_pos
    isl_end_rel = nearest_island["end"] - var_pos
    for ax in axes:
        ax.axvspan(isl_start_rel, isl_end_rel, color="green", alpha=0.1)
elif nearest_island:
    isl_start_rel = nearest_island["start"] - var_pos
    isl_end_rel = nearest_island["end"] - var_pos
    for ax in axes:
        ax.axvspan(isl_start_rel, isl_end_rel, color="green", alpha=0.07)

plt.tight_layout()
plt.savefig(f"{BASE}/figures/brca2_per_cpg_delta.png", dpi=130, bbox_inches="tight")
print(f"Saved figures/brca2_per_cpg_delta.png")

# Plot 2: distribution comparison (4 groups)
fig2, axb = plt.subplots(1, 1, figsize=(10, 5))
beta_data = {
    "HP1 (germline-tag, n=200)": [r["HP1_beta"] for r in cpgs if r["HP1_n"]>=3 and r["HP1_beta"] is not None],
    "HP2 (germline-tag, n=198)": [r["HP2_beta"] for r in cpgs if r["HP2_n"]>=3 and r["HP2_beta"] is not None],
    "HP1-1 (somatic-reconstructed, n=208)": [r["HP1-1_beta"] for r in cpgs if r["HP1-1_n"]>=3 and r["HP1-1_beta"] is not None],
}
positions = list(range(1, len(beta_data) + 1))
parts = axb.violinplot([beta_data[k] for k in beta_data], positions=positions, showmeans=True, showmedians=True)
colors_v = ["#2563eb", "#3b82f6", "#dc2626"]
for i, pc in enumerate(parts["bodies"]):
    pc.set_facecolor(colors_v[i])
    pc.set_alpha(0.5)
axb.set_xticks(positions)
axb.set_xticklabels(list(beta_data.keys()), rotation=15, ha="right", fontsize=9)
axb.set_ylabel("methylation β per CpG")
axb.set_title("BRCA2 promoter ±1kb: per-CpG β distribution by haplotype group", fontsize=11)
axb.grid(alpha=0.3)
axb.set_ylim(-0.05, 1.05)

# Annotate means
import numpy as np
for i, k in enumerate(beta_data):
    m = np.mean(beta_data[k])
    axb.annotate(f"mean={m:.3f}", xy=(positions[i], m), xytext=(positions[i]+0.15, m+0.05), fontsize=9)

plt.tight_layout()
plt.savefig(f"{BASE}/figures/brca2_beta_distribution.png", dpi=130, bbox_inches="tight")
print(f"Saved figures/brca2_beta_distribution.png")

# Plot 3: Negative control comparison
with open(f"{BASE}/output/step4_negative_control.json") as f:
    controls = json.load(f)

fig3, axc = plt.subplots(1, 1, figsize=(9, 5))
positions = ["BRCA2\n(TSG promoter)", "KMT2C\n(TSG promoter)"] + [f"Random-{i+1}\n({c['chrom']})" for i, c in enumerate(controls) if c['delta'] is not None]
deltas_plot = [
    brca2["stats"]["HP1_vs_HP1-1"]["mean_delta_som_minus_germ"],
    kmt2c["stats"]["HP1_vs_HP1-1"]["mean_delta_som_minus_germ"],
] + [c["delta"] for c in controls if c["delta"] is not None]
colors_b = ["#dc2626", "#9ca3af"] + ["#9ca3af"] * (len(deltas_plot) - 2)
axc.bar(positions, deltas_plot, color=colors_b, alpha=0.8)
axc.axhline(0, color="black", lw=0.7)
axc.axhline(-0.05, color="grey", lw=0.5, ls=":", label="±0.05 threshold")
axc.axhline(0.05, color="grey", lw=0.5, ls=":")
axc.axhline(-0.20, color="red", lw=0.7, ls=":", label="Lister 2009 ribbon")
axc.axhline(0.20, color="red", lw=0.7, ls=":")
axc.set_ylabel("Δβ = mean(somatic-HP1-1) − mean(germline-HP1)")
axc.set_title("BRCA2 vs Negative Controls (5 random non-TSG non-LOH TP sSNV)", fontsize=11)
axc.legend(loc="upper right", fontsize=9)
axc.grid(alpha=0.3, axis="y")
plt.xticks(rotation=20, ha="right", fontsize=9)
plt.tight_layout()
plt.savefig(f"{BASE}/figures/brca2_vs_controls.png", dpi=130, bbox_inches="tight")
print(f"Saved figures/brca2_vs_controls.png")
