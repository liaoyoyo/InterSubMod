#!/usr/bin/env python3
"""
42 — A6 Fig3: LOH 3-cause read-level methylation triptych (G2 visual proof).
Read-level strip plots (more reliable than headless IGV) for the 3 LOH diagnostic
exemplars, one panel per cause class:
  (a) candidate_subclone : chr11:2146993 (cnLOH, ARI=0.892) — clear HP-som vs HP-germ separation
  (b) self_phasing_artifact : chr18:64001548 (lossLOH, ARI=-0.079) — no separation
  (c) CN_regression : chr8:134018510 (gainLOH, ARI=0.523 but placebo collider, extreme baseline)
Exemplars read from loh_diagnostic_classifier.json. Black=methylated(ML>=200) /
white=unmethylated(<=50) / grey=ambiguous-or-missing. Rows grouped by HP tag.

Single-sample HCC1395. Output: figures/fig3_loh_triptych.png
"""
import json
from pathlib import Path
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts")
from lib.plot_setup import setup_plot_style
setup_plot_style()

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
LOHC = WS / "genome_survey_v2/loh_diagnostic_classifier.json"
FIG3 = WS / "figures/fig3_loh_triptych.png"
WINDOW = 600; ML_HIGH, ML_LOW = 200, 50; MIN_CORE_CPG = 4
bam = pysam.AlignmentFile(BAM, "rb")
GREEN, GREY, ORANGE = "#15803D", "#9CA3AF", "#D97757"


def hp(r):
    for t, v in r.tags:
        if t == "HP":
            return str(v)
    return None


def mcalls(r):
    o = {}
    try:
        mod = r.modified_bases
    except Exception:
        return o
    if not mod:
        return o
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False) if a is not None and b is not None}
    for k, calls in mod.items():
        if k[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None:
                o[rf] = 1 if ml >= ML_HIGH else (0 if ml <= ML_LOW else -1)
    return o


def collect(chrom, var0, groups):
    s, e = var0 - WINDOW, var0 + WINDOW
    g = {k: [] for k in groups}
    for r in bam.fetch(chrom, max(0, s), e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h not in groups:
            continue
        m = {p: st for p, st in mcalls(r).items() if s <= p <= e and st >= -1}
        if sum(1 for v in m.values() if v >= 0) >= MIN_CORE_CPG:
            g[h].append((r.reference_start, m))
    return g


def draw_panel(ax, chrom, pos, gg, gs, title, color):
    var0 = pos - 1
    g = collect(chrom, var0, [gg, gs])
    # union of CpG positions (covered by >=20% of reads)
    allreads = [m for _, m in g[gg] + g[gs]]
    from collections import Counter
    cov = Counter()
    for m in allreads:
        for c, v in m.items():
            if v >= 0:
                cov[c] += 1
    n = max(1, len(allreads))
    cpgs = sorted([c for c, k in cov.items() if k >= max(3, 0.2 * n)])
    if not cpgs:
        ax.text(0.5, 0.5, "insufficient", ha="center"); return
    cpg_idx = {c: i for i, c in enumerate(cpgs)}
    # order reads: germline group then somatic group
    row = 0
    yticks = []; ylabels = []
    for grp, gname, gcol in [(gg, f"germline ({gg})", "#1E3A8A"), (gs, f"somatic ({gs})", ORANGE)]:
        reads = sorted(g[grp], key=lambda x: -np.mean([v for v in x[1].values() if v >= 0]) if any(v >= 0 for v in x[1].values()) else 0)
        gstart = row
        for _, m in reads:
            for c, v in m.items():
                if c not in cpg_idx:
                    continue
                x = cpg_idx[c]
                col = "black" if v == 1 else ("white" if v == 0 else GREY)
                ax.add_patch(Rectangle((x, row), 1, 1, facecolor=col, edgecolor="none"))
            row += 1
        # group label bar
        ax.add_patch(Rectangle((-2.5, gstart), 1.2, max(1, row - gstart), facecolor=gcol, edgecolor="none"))
        yticks.append((gstart + row) / 2); ylabels.append(gname)
        # separator
        ax.axhline(row, color="red", lw=0.8)
    # variant column marker
    if var0 in cpg_idx:
        ax.axvline(cpg_idx[var0] + 0.5, color=ORANGE, lw=1.0, ls="--")
    ax.set_xlim(-3, len(cpgs)); ax.set_ylim(row, 0)
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(cpgs)} CpG (±{WINDOW}bp)", fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold", color=color)


def main():
    lohc = json.load(open(LOHC))
    ex = lohc['summary']['exemplars']
    fig, ax = plt.subplots(1, 3, figsize=(18, 6.5), dpi=140)
    specs = [
        ("candidate_subclone", GREEN, "(a) candidate subclone"),
        ("self_phasing_artifact", GREY, "(b) self-phasing artifact"),
        ("CN_regression", ORANGE, "(c) CN / regression"),
    ]
    for i, (klass, col, tag) in enumerate(specs):
        e = ex.get(klass)
        if not e:
            ax[i].text(0.5, 0.5, f"no {klass} exemplar", ha="center"); continue
        chrom, pos = e['locus'].split(":"); pos = int(pos)
        gg, gs = ("1", "1-1") if e['axis'] == "HP1_vs_HP1-1" else ("2", "2-1")
        title = (f"{tag}\n{e['locus']} · {e['cn_class']} · ARI={e['ari']}\n"
                 f"n_cpg={e['n_cpg']} germβ={e['germ_beta']:.2f} Δβ={e['delta']:.3f}")
        draw_panel(ax[i], chrom, pos, gg, gs, title, col)

    # legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor="black", label="methylated (ML≥200)"),
               Patch(facecolor="white", edgecolor="grey", label="unmethylated (ML≤50)"),
               Patch(facecolor=GREY, label="ambiguous / missing")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Fig3 · G2 LOH 表觀雙 haplotype 三成因 read-level 甲基證據 (HCC1395 單樣本)\n"
                 "每行=1 read; 上=germline-tag 下=somatic-tag; 橘虛線=變異位置",
                 fontsize=13, fontweight="bold", y=1.04)
    plt.tight_layout()
    plt.savefig(FIG3, facecolor="white", bbox_inches="tight")
    print(f"Fig3 -> {FIG3} ({FIG3.stat().st_size//1024} KB)")
    for klass, _, _ in specs:
        e = ex.get(klass)
        print(f"  {klass}: {e['locus'] if e else 'NA'} ARI={e['ari'] if e else 'NA'}")


if __name__ == "__main__":
    main()
