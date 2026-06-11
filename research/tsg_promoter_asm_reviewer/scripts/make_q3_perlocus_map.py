#!/usr/bin/env python3
"""
make_q3_perlocus_map.py — Q3 per-locus consistency map.

X = ari_blind_allele (blind methylation cluster vs ALT/REF allele agreement)
Y = ari_blind_hp     (blind methylation cluster vs germline HP-tag agreement)
color = class (het-NULL / TP / FP / FN)
size  = n (read depth)
marker edge / hatch = LOH
diagonal y=x        = HP == ALT/REF (germline-allele axis collapses)

Loci ON the diagonal: HP axis == allele axis (germline ASM, Layer A).
Loci BELOW-RIGHT (allele-only): cluster maps to ALT/REF but HP is blind  -> somatic-allele divergent.
Loci ABOVE-LEFT  (hp-only):     cluster maps to HP but ALT/REF is blind  -> HP/LOH-driven divergent.

Read-only inputs: genome_survey_v2/enriched_perlocus.json
Output:           figures/q3_perlocus_map.png
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod")
try:
    from scripts.lib.plot_setup import setup_plot_style

    setup_plot_style(base_size=11, dpi=150)
except Exception:
    pass

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
SRC = os.path.join(ROOT, "genome_survey_v2", "enriched_perlocus.json")
OUT = os.path.join(ROOT, "figures", "q3_perlocus_map.png")

d = json.load(open(SRC))

# Only loci with HP tags can be placed on the 2-axis map (135/186).
pts = [x for x in d if x["ari_blind_hp"] is not None and x["ari_blind_allele"] is not None]
n_no_hp = len(d) - len(pts)

CLS_COLOR = {
    "het-NULL": "#1f77b4",  # blue  - germline reference
    "TP": "#2ca02c",        # green - true somatic
    "FP": "#d62728",        # red   - false positive
    "FN": "#ff7f0e",        # orange- false negative
}
CLS_ORDER = ["het-NULL", "TP", "FP", "FN"]

fig, ax = plt.subplots(figsize=(11.4, 8.0))
fig.subplots_adjust(right=0.74)  # reserve right margin for legends OUTSIDE the axes

# diagonal HP == ALT/REF
ax.plot([-0.15, 1.0], [-0.15, 1.0], ls="--", color="#555555", lw=1.2, zorder=1)
ax.text(0.92, 0.97, "HP $\\equiv$ ALT/REF\n(germline-allele axis)",
        fontsize=9, color="#555555", ha="right", va="top",
        transform=ax.transAxes, style="italic")

# divergence-band guides
ax.axhline(0.3, color="#bbbbbb", lw=0.6, ls=":", zorder=0)
ax.axvline(0.3, color="#bbbbbb", lw=0.6, ls=":", zorder=0)

# size scaling by read depth n — smaller range so large points do not bury small ones
def msize(n):
    return 16 + (n or 0) * 0.95

for cls in CLS_ORDER:
    sub = [x for x in pts if x["cls"] == cls]
    # non-LOH
    nl = [x for x in sub if not x["loh"]]
    if nl:
        ax.scatter([x["ari_blind_allele"] for x in nl],
                   [x["ari_blind_hp"] for x in nl],
                   s=[msize(x["n"]) for x in nl],
                   c=CLS_COLOR[cls], alpha=0.55, edgecolors="white",
                   linewidths=0.5, zorder=3, label=None)
    # LOH (heavy black edge marker)
    lh = [x for x in sub if x["loh"]]
    if lh:
        ax.scatter([x["ari_blind_allele"] for x in lh],
                   [x["ari_blind_hp"] for x in lh],
                   s=[msize(x["n"]) for x in lh],
                   c=CLS_COLOR[cls], alpha=0.55, edgecolors="black",
                   linewidths=1.6, zorder=4, label=None)

# annotate the decisive divergent loci
def find(loc):
    for x in pts:
        if x["locus"] == loc:
            return x
    return None

# annotate only the single decisive lead per regime (avoid text-on-point clutter)
ann_allele_only = ["chr16:34707014"]
ann_hp_only = ["chr5:173198872"]

for loc in ann_allele_only:
    x = find(loc)
    if x:
        ax.annotate(loc + ("\n(VAF %.2f, %s)" % (x["vaf"], x["cls"])),
                    xy=(x["ari_blind_allele"], x["ari_blind_hp"]),
                    xytext=(x["ari_blind_allele"] + 0.02, x["ari_blind_hp"] - 0.18),
                    fontsize=7.5, color="#8B0000",
                    arrowprops=dict(arrowstyle="->", color="#8B0000", lw=0.8))
for loc in ann_hp_only:
    x = find(loc)
    if x:
        ax.annotate(loc + ("\n(VAF %.2f, %s, LOH)" % (x["vaf"], x["cls"])),
                    xy=(x["ari_blind_allele"], x["ari_blind_hp"]),
                    xytext=(x["ari_blind_allele"] - 0.05, x["ari_blind_hp"] + 0.12),
                    fontsize=7.5, color="#00008B",
                    arrowprops=dict(arrowstyle="->", color="#00008B", lw=0.8))

# regime region labels
ax.text(0.72, 0.06, "ALLELE-ONLY\n(somatic divergent:\nALT/REF clusters, HP blind)",
        fontsize=8, color="#8B0000", ha="left", va="bottom", weight="bold")
ax.text(0.10, 0.80, "HP-ONLY\n(HP/LOH-driven:\nHP clusters, ALT/REF blind)",
        fontsize=8, color="#00008B", ha="left", va="center", weight="bold")

ax.set_xlabel("ari_blind_allele  (blind methylation cluster $\\leftrightarrow$ ALT/REF)", fontsize=11)
ax.set_ylabel("ari_blind_hp  (blind methylation cluster $\\leftrightarrow$ germline HP)", fontsize=11)
ax.set_title("Q3 per-locus consistency map (n=%d loci with HP tags; %d HP-absent omitted)\n"
             "On-diagonal = germline-allele axis (HP $\\equiv$ ALT/REF); off-diagonal = somatic divergence"
             % (len(pts), n_no_hp), fontsize=11)
ax.set_xlim(-0.15, 1.02)
ax.set_ylim(-0.15, 1.02)
ax.grid(True, alpha=0.18)

# legends: class color + LOH edge + size
class_handles = [Line2D([0], [0], marker="o", color="w", label=c,
                        markerfacecolor=CLS_COLOR[c], markersize=9) for c in CLS_ORDER]
loh_handles = [
    Line2D([0], [0], marker="o", color="w", label="non-LOH",
           markerfacecolor="gray", markeredgecolor="white", markersize=9),
    Line2D([0], [0], marker="o", color="w", label="LOH",
           markerfacecolor="gray", markeredgecolor="black", markeredgewidth=1.6, markersize=9),
]
size_handles = [
    Line2D([0], [0], marker="o", color="w", label="n=20",
           markerfacecolor="gray", markersize=(msize(20) ** 0.5)),
    Line2D([0], [0], marker="o", color="w", label="n=60",
           markerfacecolor="gray", markersize=(msize(60) ** 0.5)),
    Line2D([0], [0], marker="o", color="w", label="n=100",
           markerfacecolor="gray", markersize=(msize(100) ** 0.5)),
]
# all three legends placed OUTSIDE the plotting area (right margin) so no data point is ever blocked
leg1 = ax.legend(handles=class_handles, title="class", loc="upper left",
                 fontsize=9, framealpha=0.95, bbox_to_anchor=(1.02, 1.00))
ax.add_artist(leg1)
leg2 = ax.legend(handles=loh_handles, title="edge", loc="upper left",
                 fontsize=9, framealpha=0.95, bbox_to_anchor=(1.02, 0.70))
ax.add_artist(leg2)
ax.legend(handles=size_handles, title="size = read depth n", loc="upper left",
          fontsize=8, framealpha=0.95, labelspacing=1.3, bbox_to_anchor=(1.02, 0.46))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
plt.savefig(OUT, dpi=150, bbox_inches="tight")
print("wrote", OUT)

# quick provenance summary
ao = [x for x in pts if x["cls"] in ("TP", "FP", "FN")
      and x["ari_blind_allele"] > 0.3 and (x["ari_blind_allele"] - x["ari_blind_hp"]) > 0.2
      and x["ari_hp_allele"] is not None and x["ari_hp_allele"] < 0.5]
ho = [x for x in pts if x["cls"] in ("TP", "FP", "FN")
      and x["ari_blind_hp"] > 0.3 and (x["ari_blind_hp"] - x["ari_blind_allele"]) > 0.2]
print("allele-only somatic:", len(ao), "hp-only somatic (diff>0.2):", len(ho))
