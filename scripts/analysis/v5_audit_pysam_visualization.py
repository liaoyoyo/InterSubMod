#!/usr/bin/env python3
"""V5 vs V3-Fixed audit: pysam-based read-level HP tag visualization.

Generates IGV-like read alignment plots for 12 selected sites,
side-by-side V3-Fixed vs V5, color-coded by HP tag, with
quantitative summary annotations.

Output: docs/reports/pi_reports/2026/04/figures/igv_v5_audit/
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pysam
from collections import defaultdict

V3F = "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam"
V5  = "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam"

OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 12 selected sites + 3 self-phasing extreme = 15 sites total (D class added)
SITES = [
    # Class A: Phase 4 TP (5)
    ("A_TP01", "chr6", 145444893, "G>A", "Phase4 TP_01: allele-only, cleanest ALT cluster"),
    ("A_TP02", "chr4", 70548355, "G>A", "Phase4 TP_02: HP0 background with allele separation"),
    ("A_TP03", "chr5", 153209947, "C>A", "Phase4 TP_03: low-methyl background micro-difference"),
    ("A_TP04", "chr16", 35118902, "G>A", "Phase4 TP_04: bimodal, ALT low-methyl"),
    ("A_TP05", "chr7", 109185781, "G>T", "Phase4 TP_05: HP+Allele co-linear, high power"),
    # Class B: Phase 4 FP (4)
    ("B_FPA1", "chr8", 93565727, "C>T", "Phase4 FP_A1: low VAF (2/68), insufficient cluster"),
    ("B_FPA2", "chr9", 137953060, "T>C", "Phase4 FP_A2: high CpG false-positive"),
    ("B_FPB1", "chr7", 52087777, "A>T", "Phase4 FP_B1: HP-driven, near SEQC2 INDEL"),
    ("B_FPB2", "chr9", 75383880, "T>A", "Phase4 FP_B2: MNP, HP-driven"),
    # Class C: V3F → V5 reassignment hotspots (3)
    ("C_V5max1", "chr19", 4639528, "(reassign)", "V5 reassigned 39 reads HP:i:33 → HP:i:11"),
    ("C_V5max2", "chr19", 2235521, "(reassign)", "V5 reassigned 26 reads HP:i:33 → HP:i:11"),
    ("C_V5max3", "chr19", 7405500, "(reassign)", "V5 reassigned 18 reads HP:i:33 → HP:i:21"),
    # Class D: Self-phasing extreme bias (3)
    ("D_SP1", "chr19", 17565944, "(extreme)", "Self-phasing extreme: HP2=113, HP1=0, ratio=∞"),
    ("D_SP2", "chr19", 12452332, "(extreme)", "Self-phasing extreme: HP2=109, HP1=1, ratio=109:1"),
    ("D_SP3", "chr19", 12467180, "(extreme)", "Self-phasing extreme: HP2=108, HP1=0, ratio=108:1"),
]

HP_COLORS = {
    1:  "#2E7D32",   # germline HP1: green
    2:  "#66BB6A",   # germline HP2: light green
    11: "#1565C0",   # somatic HP1-1: blue
    21: "#42A5F5",   # somatic HP2-1: light blue
    33: "#E65100",   # somatic ambiguous: orange
    0:  "#9E9E9E",   # untagged: grey
}
HP_LABELS = {
    1: "HP1 (germ)", 2: "HP2 (germ)",
    11: "HP11 (som-1)", 21: "HP21 (som-2)",
    33: "HP33 (amb)", 0: "untagged",
}


def fetch_reads(bam_path, chrom, pos, window=200):
    """Return reads covering pos with ~window bp range, sorted by start."""
    out = []
    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam.fetch(chrom, max(0, pos - window), pos + window):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.reference_start > pos or r.reference_end < pos:
            continue
        try:
            hp = r.get_tag("HP")
        except KeyError:
            hp = 0
        out.append({
            "name":  r.query_name,
            "start": r.reference_start,
            "end":   r.reference_end,
            "strand": "-" if r.is_reverse else "+",
            "hp": hp,
        })
    bam.close()
    out.sort(key=lambda x: (x["hp"] not in (1, 2, 11, 21, 33), x["hp"], x["start"]))
    return out


def draw_reads_panel(ax, reads, pos, window, title, max_reads=80):
    """Draw read tracks colored by HP tag."""
    if not reads:
        ax.text(0.5, 0.5, "No reads in window", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="grey")
        ax.set_title(title, fontsize=11)
        return

    reads = reads[:max_reads]
    n = len(reads)
    ax.set_xlim(pos - window, pos + window)
    ax.set_ylim(-1, n + 1)

    # Group separator (between HP groups)
    last_hp = None
    for i, r in enumerate(reads):
        y = n - i
        color = HP_COLORS.get(r["hp"], "#9E9E9E")
        ax.plot([r["start"], r["end"]], [y, y], "-", color=color, lw=1.5,
                solid_capstyle="butt")
        if last_hp is not None and r["hp"] != last_hp:
            ax.axhline(y + 0.5, color="black", lw=0.3, alpha=0.4)
        last_hp = r["hp"]

    # Variant line
    ax.axvline(pos, color="red", lw=1.2, alpha=0.7, ls="--")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel(f"position (window {pos-window}-{pos+window})")
    ax.set_yticks([])
    ax.set_xticks([pos - window, pos, pos + window])
    ax.tick_params(axis="x", labelsize=8)


def render_site(case_id, chrom, pos, alt, comment):
    """Render one site as 2-panel V3-F vs V5."""
    window = 1000
    reads_v3f = fetch_reads(V3F, chrom, pos, window=window)
    reads_v5  = fetch_reads(V5, chrom, pos, window=window)

    # Quantitative summary
    cnt_v3f = defaultdict(int)
    cnt_v5  = defaultdict(int)
    for r in reads_v3f: cnt_v3f[r["hp"]] += 1
    for r in reads_v5:  cnt_v5[r["hp"]] += 1

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))
    title_v3f = f"V3-Fixed (baseline)\n{HP_summary(cnt_v3f)}"
    title_v5  = f"V5 Somatic Fallback\n{HP_summary(cnt_v5)}"
    draw_reads_panel(axes[0], reads_v3f, pos, window, title_v3f)
    draw_reads_panel(axes[1], reads_v5, pos, window, title_v5)

    # Legend (shared)
    legend_patches = [
        mpatches.Patch(color=HP_COLORS[k], label=HP_LABELS[k])
        for k in [1, 2, 11, 21, 33, 0]
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=6, fontsize=10,
               bbox_to_anchor=(0.5, -0.02))

    delta_33 = cnt_v5.get(33, 0) - cnt_v3f.get(33, 0)
    delta_dir = (cnt_v5.get(11, 0) + cnt_v5.get(21, 0)) - (cnt_v3f.get(11, 0) + cnt_v3f.get(21, 0))
    suptitle = (f"{case_id}: {chrom}:{pos:,} {alt}    "
                f"|    {comment}    |    "
                f"Δ HP:i:33 = {delta_33:+d}    Δ directional = {delta_dir:+d}")
    fig.suptitle(suptitle, fontsize=13, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    out_path = OUT_DIR / f"{case_id}_{chrom}_{pos}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path, dict(cnt_v3f), dict(cnt_v5)


def HP_summary(cnt):
    order = [1, 2, 11, 21, 33, 0]
    parts = [f"HP{k}={cnt.get(k, 0)}" for k in order if cnt.get(k, 0) > 0]
    total = sum(cnt.values())
    return f"n={total}: " + ", ".join(parts)


def main():
    rows = []
    for case_id, chrom, pos, alt, comment in SITES:
        path, t3, t5 = render_site(case_id, chrom, pos, alt, comment)
        d33 = t5.get(33, 0) - t3.get(33, 0)
        d_dir = (t5.get(11, 0) + t5.get(21, 0)) - (t3.get(11, 0) + t3.get(21, 0))
        print(f"  {case_id:<10} {chrom}:{pos:<10} V3F={dict(t3)}  V5={dict(t5)}  Δ33={d33:+d}  Δdir={d_dir:+d}")
        rows.append((case_id, chrom, pos, alt, t3, t5, d33, d_dir, comment))

    # Save TSV summary
    tsv = OUT_DIR / "v5_audit_summary.tsv"
    with open(tsv, "w") as f:
        f.write("case\tchrom\tpos\talt\tcomment\tv3f_HP1\tv3f_HP2\tv3f_HP11\tv3f_HP21\tv3f_HP33\tv3f_HP0\tv5_HP1\tv5_HP2\tv5_HP11\tv5_HP21\tv5_HP33\tv5_HP0\tdelta_HP33\tdelta_directional\n")
        for case_id, chrom, pos, alt, t3, t5, d33, d_dir, comment in rows:
            f.write(f"{case_id}\t{chrom}\t{pos}\t{alt}\t{comment}\t")
            for k in [1, 2, 11, 21, 33, 0]:
                f.write(f"{t3.get(k, 0)}\t")
            for k in [1, 2, 11, 21, 33, 0]:
                f.write(f"{t5.get(k, 0)}\t")
            f.write(f"{d33}\t{d_dir}\n")
    print(f"\nSummary saved to {tsv}")


if __name__ == "__main__":
    main()
