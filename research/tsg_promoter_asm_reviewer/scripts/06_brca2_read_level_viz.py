#!/usr/bin/env python3
"""
Read-level IGV-style visualization for BRCA2 promoter ±1kb.
Each row = 1 read, dots = per-CpG methylation calls colored by methylation state.
Groups: HP1 (top) / HP1-1 (middle) / HP2 (bottom).
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pysam
from collections import defaultdict

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"

CHROM = "chr13"
VAR_POS = 32315128
WINDOW = 1000
ML_HIGH = 200
ML_LOW = 50

bam = pysam.AlignmentFile(BAM, "rb")

def extract_hp(read):
    for tag, val in read.tags:
        if tag == "HP":
            return str(val)
    return None

def get_meth_calls(read):
    mods = {}
    try:
        modified = read.modified_bases
    except Exception:
        return mods
    if not modified:
        return mods
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    rd_to_ref = {rd: rf for rd, rf in aligned_pairs if rd is not None and rf is not None}
    for key, calls in modified.items():
        base, strand, mod_code = key
        if mod_code != 'm':
            continue
        for read_pos, ml in calls:
            ref_pos = rd_to_ref.get(read_pos)
            if ref_pos is None:
                continue
            mods[ref_pos] = ml
    return mods

# Collect reads by HP group
start = VAR_POS - 1 - WINDOW
end = VAR_POS - 1 + WINDOW

groups = {"HP1": [], "HP1-1": [], "HP2": []}
for read in bam.fetch(CHROM, start, end):
    if read.flag & 0x100 or read.flag & 0x800 or read.flag & 0x400:
        continue
    hp = extract_hp(read)
    if hp not in ("1", "1-1", "2"):
        continue
    key = f"HP{hp}"
    mods = get_meth_calls(read)
    # filter to window
    in_window = {pos: ml for pos, ml in mods.items() if start <= pos <= end}
    if len(in_window) >= 3:
        groups[key].append({
            "read_name": read.query_name,
            "read_start": read.reference_start,
            "read_end": read.reference_end,
            "mods": in_window,
        })

print(f"BRCA2 reads with >=3 CpG calls in ±{WINDOW}bp:")
for k, v in groups.items():
    print(f"  {k}: {len(v)} reads")

# Build figure
fig, ax = plt.subplots(figsize=(14, 9))

color_map = {"HP1": "#2563eb", "HP1-1": "#dc2626", "HP2": "#10b981"}
y_offset = 0
group_y_starts = {}

for group_name in ["HP1", "HP1-1", "HP2"]:
    group_reads = groups[group_name]
    # Sort by read_start for visual coherence
    group_reads = sorted(group_reads, key=lambda r: r["read_start"])
    group_y_starts[group_name] = y_offset

    # Group label box
    ax.add_patch(Rectangle((start - VAR_POS - 50, y_offset - 0.5),
                            -200, len(group_reads),
                            color=color_map[group_name], alpha=0.3, zorder=0))
    ax.text(start - VAR_POS - 250, y_offset + len(group_reads)/2,
            f"{group_name}\n(n={len(group_reads)})",
            ha="right", va="center", fontsize=11, fontweight="bold", color=color_map[group_name])

    for ri, r in enumerate(group_reads):
        # Read backbone (thin gray line spanning read extent within window)
        read_x_start = max(r["read_start"], start) - VAR_POS
        read_x_end = min(r["read_end"], end) - VAR_POS
        ax.plot([read_x_start, read_x_end], [y_offset, y_offset],
                color="#cbd5e1", lw=0.5, zorder=1)
        # CpG calls
        for pos, ml in r["mods"].items():
            x = pos - VAR_POS
            if ml >= ML_HIGH:
                color = "#000000"  # methylated = black
                sz = 12
            elif ml <= ML_LOW:
                color = "#ffffff"  # unmethylated = white
                sz = 12
            else:
                color = "#888888"  # ambiguous = grey
                sz = 5
            ax.plot(x, y_offset, marker="o", markersize=sz/3,
                    markerfacecolor=color, markeredgecolor="black", markeredgewidth=0.3, zorder=2)
        y_offset += 1
    y_offset += 5  # gap between groups

# Variant marker
ax.axvline(0, color="orange", lw=2, ls="--", zorder=3, label=f"somatic SNV {CHROM}:{VAR_POS} G>A (TVAF=0.189)")

# CpG island annotation
ax.axvspan(32315396 - VAR_POS, 32315763 - VAR_POS, color="green", alpha=0.15, label="CpG island (chr13:32315396-32315763)")

ax.set_xlim(start - VAR_POS - 300, end - VAR_POS + 100)
ax.set_ylim(-2, y_offset)
ax.set_xlabel(f"Distance from somatic SNV (bp) — {CHROM}:{VAR_POS}", fontsize=11)
ax.set_ylabel("Reads (grouped by HP:Z tag)", fontsize=11)
ax.set_title(
    f"BRCA2 promoter — Read-level methylation by haplotype group\n"
    f"Black=methylated (ML>={ML_HIGH}), White=unmethylated (ML<={ML_LOW}), Grey=ambiguous, Vertical dashed=somatic SNV\n"
    f"HP1-1 (somatic-reconstructed) shows visible hypomethylation vs HP1 (germline)",
    fontsize=11)
ax.legend(loc="lower right", fontsize=9)
ax.set_yticks([])

# Add per-position summary at top
plt.tight_layout()
plt.savefig(f"{BASE}/figures/brca2_read_level_methylation.png", dpi=130, bbox_inches="tight")
print(f"Saved figures/brca2_read_level_methylation.png")

# Also save IGV session XML for manual loading
igv_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg38" hasGeneTrack="true" locus="chr13:32314500-32316000" version="8">
    <Resources>
        <Resource path="{BAM}"/>
    </Resources>
    <Panel height="500" name="Panel_Reads">
        <Track autoScale="false" clazz="org.broad.igv.sam.AlignmentTrack"
               colorOption="HP_TAG" id="{BAM}"
               name="HCC1395 paired (HP:Z colored)"
               showSpliceJunctions="false" snpThreshold="0.05" sortable="true" visible="true"/>
    </Panel>
    <PanelLayout dividerFractions="0.3"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>
"""
with open(f"{BASE}/output/brca2_igv_session.xml", "w") as f:
    f.write(igv_xml)
print(f"Saved output/brca2_igv_session.xml")
