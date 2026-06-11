#!/usr/bin/env python3
"""
Human-eye verification: read-level methylation matrix for top robust clear&consistent ASM loci.
Reads grouped by ALLELE at the somatic SNV (REF vs ALT — universal, works in LOH where HP1-1 absent),
with HP:Z tag shown as a left color-bar. Per-CpG 5mC dots colored by methylation state.
Generalizes 06_brca2_read_level_viz.py to a locus list.

Output: figures/human_eye/<chr>_<pos>_readlevel.png per locus.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pysam

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
OUT = f"{BASE}/figures/human_eye"
os.makedirs(OUT, exist_ok=True)
WINDOW = 1000
ML_HIGH, ML_LOW = 200, 50

# top robust clear&consistent loci (from phase2 dataverify; chrom,pos,ref,alt,direction,dHP,dALLELE,loh)
LOCI = [
    ("chr13", 32315128, "G", "A", "hypo",  "-0.122", "-0.099", "nonLOH", "BRCA2/ZAR1L (validated anchor)"),
    ("chr9",  42376881, None, None, "hyper", "+0.666", "+0.285", "nonLOH", "biggest HP-axis effect"),
    ("chr22", 34684238, None, None, "hypo",  "-0.209", "-0.601", "nonLOH", "strong hypo both axes"),
    ("chr2",  23386385, None, None, "hypo",  "-0.152", "-0.152", "LOH",    "most significant p=4e-24, symmetric"),
    ("chr17", 9577080,  None, None, "hyper", "+0.584", "+0.592", "LOH",    "biggest combined effect"),
    ("chr12", 31601630, None, None, "hyper", "+0.581", "+0.213", "nonLOH", "strong hyper"),
    ("chr16", 17774746, None, None, "hypo",  "-0.221", "-0.256", "nonLOH", "hypo both axes"),
    ("chr8",  142166567, None, None, "hypo", "-0.282", "-0.278", "LOH",    "symmetric hypo"),
]

bam = pysam.AlignmentFile(BAM, "rb")

def extract_hp(read):
    for tag, val in read.tags:
        if tag == "HP":
            return str(val)
    return None

def base_at(read, ref_pos0):
    """Return read base aligned to ref_pos0 (0-based), or None."""
    for rd, rf in read.get_aligned_pairs(matches_only=True):
        if rf == ref_pos0:
            return read.query_sequence[rd] if read.query_sequence else None
    return None

def get_meth_calls(read):
    mods = {}
    try:
        modified = read.modified_bases
    except Exception:
        return mods
    if not modified:
        return mods
    rd_to_ref = {rd: rf for rd, rf in read.get_aligned_pairs(matches_only=False)
                 if rd is not None and rf is not None}
    for key, calls in modified.items():
        if key[2] != 'm':   # 5mC only for visualization
            continue
        for read_pos, ml in calls:
            rp = rd_to_ref.get(read_pos)
            if rp is not None:
                mods[rp] = ml
    return mods

HP_COLOR = {"1": "#2563eb", "1-1": "#dc2626", "2": "#10b981", "2-1": "#f59e0b", "none": "#9ca3af"}

def plot_locus(chrom, pos, ref, alt, direction, dHP, dAL, loh, note):
    var0 = pos - 1
    start, end = var0 - WINDOW, var0 + WINDOW
    # infer ref/alt from pileup if not given
    reads = []
    for read in bam.fetch(chrom, max(0, start), end):
        if read.flag & 0x900 or read.flag & 0x400:
            continue
        mods = {p: m for p, m in get_meth_calls(read).items() if start <= p <= end}
        if len(mods) < 3:
            continue
        b = base_at(read, var0)
        reads.append({"name": read.query_name, "rs": read.reference_start,
                      "re": read.reference_end, "mods": mods,
                      "hp": extract_hp(read) or "none", "base": b})
    if not reads:
        print(f"  {chrom}:{pos} no reads"); return None
    # determine ref/alt allele by frequency if not supplied
    from collections import Counter
    bc = Counter(r["base"] for r in reads if r["base"])
    common = [b for b, _ in bc.most_common()]
    if ref is None:
        ref = common[0] if common else "?"
    if alt is None:
        alt = next((b for b in common if b != ref), "?")
    # group by allele
    grp = {"ALT (%s)" % alt: [], "REF (%s)" % ref: [], "other": []}
    for r in reads:
        if r["base"] == alt:
            grp["ALT (%s)" % alt].append(r)
        elif r["base"] == ref:
            grp["REF (%s)" % ref].append(r)
        else:
            grp["other"].append(r)

    fig, ax = plt.subplots(figsize=(14, 9))
    y = 0
    order = ["ALT (%s)" % alt, "REF (%s)" % ref, "other"]
    for gname in order:
        g = sorted(grp[gname], key=lambda r: r["rs"])
        if not g:
            continue
        ax.add_patch(Rectangle((-WINDOW - 250, y - 0.5), 180, len(g),
                               color="#64748b", alpha=0.12, zorder=0))
        ax.text(-WINDOW - 300, y + len(g) / 2, f"{gname}\n(n={len(g)})",
                ha="right", va="center", fontsize=11, fontweight="bold")
        for r in g:
            xs, xe = max(r["rs"], start) - var0, min(r["re"], end) - var0
            ax.plot([xs, xe], [y, y], color="#e2e8f0", lw=0.5, zorder=1)
            # HP color bar (left)
            ax.add_patch(Rectangle((-WINDOW - 60, y - 0.4), 40, 0.8,
                                   color=HP_COLOR.get(r["hp"], "#9ca3af"), zorder=2))
            for p, ml in r["mods"].items():
                x = p - var0
                c = "#000000" if ml >= ML_HIGH else ("#ffffff" if ml <= ML_LOW else "#888888")
                ax.plot(x, y, "o", markersize=3, markerfacecolor=c,
                        markeredgecolor="black", markeredgewidth=0.3, zorder=3)
            y += 1
        y += 4
    ax.axvline(0, color="orange", lw=2, ls="--", zorder=4,
               label=f"somatic SNV {chrom}:{pos} {ref}>{alt}")
    # HP legend
    from matplotlib.lines import Line2D
    hp_legend = [Line2D([0], [0], marker='s', color='w', markerfacecolor=HP_COLOR[k],
                        markersize=9, label=f"HP:{k}") for k in ["1", "1-1", "2", "2-1"]]
    ax.set_xlim(-WINDOW - 350, WINDOW + 100)
    ax.set_ylim(-2, y)
    ax.set_xlabel(f"Distance from somatic SNV (bp) — {chrom}:{pos}", fontsize=11)
    ax.set_ylabel("Reads (grouped by ALLELE; left bar = HP tag)", fontsize=11)
    ax.set_title(
        f"{chrom}:{pos} {note} [{loh}, {direction}]  ALLELE-axis Δβ={dAL}, HP-axis Δβ={dHP}\n"
        f"Black=5mC methylated(ML>={ML_HIGH}) White=unmeth(<={ML_LOW}) Grey=ambiguous | "
        f"ALT-allele reads (top) vs REF-allele reads = somatic ASM",
        fontsize=10)
    ax.legend(handles=hp_legend + [Line2D([0], [0], color='orange', ls='--', label='somatic SNV')],
              loc="lower right", fontsize=8, ncol=2)
    ax.set_yticks([])
    plt.tight_layout()
    fn = f"{OUT}/{chrom}_{pos}_readlevel.png"
    plt.savefig(fn, dpi=120, bbox_inches="tight")
    plt.close()
    ng = {k: len(v) for k, v in grp.items() if v}
    print(f"  {chrom}:{pos} {note} -> {os.path.basename(fn)}  groups={ng}")
    return fn

written = []
for L in LOCI:
    f = plot_locus(*L)
    if f:
        written.append(f)
print(f"\nSaved {len(written)} human-eye read-level figures to {OUT}/")
