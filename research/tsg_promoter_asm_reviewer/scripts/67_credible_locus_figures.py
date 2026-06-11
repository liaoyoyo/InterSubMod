#!/usr/bin/env python3
"""
67 - Per-credible-locus methylation observation figures (read x CpG heatmap by HP).

For each sample's credible (pass_tierA) loci, render the ACTUAL methylation
observation: a read x CpG matrix with reads grouped by haplotype (germline main HP
vs somatic sub-haplotype HPx-1 vs other HP), so the allele-specific methylation
clustering can be VISUALLY confirmed locus-by-locus (見樹; 肉眼檢視, per
feedback_visual_inspection_requirement).

Reuses validated pysam extraction (mcalls binarized ML>=200 meth / <=50 unmeth /
else intermediate; HP:Z tag). Window +/-600 (matches ARI core, script 30).

USAGE: python3 67_credible_locus_figures.py <SAMPLE> [TIER=relax|strict]
  -> genome_survey_v2/cn_confound/cross_sample/<SAMPLE>_locus_figs.json
     (per-locus base64 PNG + annotations; consumed by 68 gallery builder)
"""
import os, sys, json, base64, io
os.environ.setdefault("TMPDIR", "/big7_disk")
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"

RUNDIR = {
    "HCC1395": "20260314_HCC1395_paired_full_full_complete_matrix",
    "COLO829": "20260315_COLO829_paired_full_full_complete_matrix",
    "H1437":   "20260315_H1437_paired_full_full_complete_matrix",
    "H2009":   "20260315_H2009_paired_full_full_complete_matrix",
    "HCC1937": "20260315_HCC1937_paired_full_full_complete_matrix",
    "HCC1954": "20260315_HCC1954_paired_full_full_complete_matrix",
}
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
WINDOW = 600
ML_HIGH, ML_LOW = 200, 50
MAX_READS_PER_GRP = 45     # subsample for readable heatmap


def bam_path(s):
    return (f"/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{s}/paired_full/"
            f"{RUNDIR[s]}/longphase_s/{s}_tagged.bam")


def hp(r):
    try:
        return str(r.get_tag("HP"))
    except KeyError:
        return None


def mcalls(r, s, e):
    """ref_pos -> 1(meth)/0(unmeth)/-1(intermediate) for 5mC in [s,e]."""
    o = {}
    try:
        mod = r.modified_bases
    except Exception:
        return o
    if not mod:
        return o
    r2 = {a: b for a, b in r.get_aligned_pairs(matches_only=False)
          if a is not None and b is not None}
    for k, calls in mod.items():
        if k[2] != 'm':
            continue
        for rp, ml in calls:
            rf = r2.get(rp)
            if rf is not None and s <= rf <= e:
                o[rf] = 1 if ml >= ML_HIGH else (0 if ml <= ML_LOW else -1)
    return o


def render_locus(bam, chrom, pos, axis, ann):
    var0 = pos - 1
    s, e = max(0, var0 - WINDOW), var0 + WINDOW
    main, sub = ("1", "1-1") if axis == "HP1_vs_HP1-1" else ("2", "2-1")
    other = "2" if main == "1" else "1"
    groups = {main: [], sub: [], other: []}
    for r in bam.fetch(chrom, s, e):
        if r.flag & 0x900 or r.flag & 0x400:
            continue
        h = hp(r)
        if h not in groups:
            continue
        m = mcalls(r, s, e)
        if len(m) >= 3:
            groups[h].append(m)
    # subsample
    rng = np.random.RandomState(42)
    for g in groups:
        if len(groups[g]) > MAX_READS_PER_GRP:
            sel = rng.choice(len(groups[g]), MAX_READS_PER_GRP, replace=False)
            groups[g] = [groups[g][i] for i in sorted(sel)]
    # CpG columns = union of positions
    cols = sorted(set(c for g in groups.values() for rd in g for c in rd))
    if not cols or sum(len(g) for g in groups.values()) < 6:
        return None
    cidx = {c: i for i, c in enumerate(cols)}
    rows = []
    rowgrp = []
    grp_order = [(main, "germline HP%s" % main, "#1d4ed8"),
                 (sub, "somatic HP%s" % sub, "#16a34a"),
                 (other, "other HP%s" % other, "#9333ea")]
    for g, lab, col in grp_order:
        # sort reads within group by mean methylation (meth frac) for readability
        def meanmeth(rd):
            vals = [v for v in rd.values() if v in (0, 1)]
            return np.mean(vals) if vals else 0.5
        for rd in sorted(groups[g], key=meanmeth, reverse=True):
            row = np.full(len(cols), np.nan)
            for c, v in rd.items():
                row[cidx[c]] = v
            rows.append(row); rowgrp.append(g)
    M = np.array(rows)
    # color: nan->white(3), -1->gray(2)... map to 0,1,2,3
    Mp = np.where(np.isnan(M), 3, np.where(M == -1, 2, M)).astype(int)
    cmap = ListedColormap(["#2563eb", "#dc2626", "#d1d5db", "#ffffff"])  # 0 unmeth,1 meth,2 inter,3 missing
    norm = BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], cmap.N)

    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    ax.imshow(Mp, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    # HP group separators
    y = 0
    for g, lab, col in grp_order:
        n = sum(1 for x in rowgrp if x == g)
        if n == 0:
            continue
        ax.axhline(y - 0.5, color="black", lw=0.8)
        ax.text(-0.012 * len(cols), y + n / 2, lab, va="center", ha="right",
                fontsize=7, color=col, fontweight="bold", rotation=90)
        y += n
    ax.axhline(y - 0.5, color="black", lw=0.8)
    # mark variant column
    vcol = min(range(len(cols)), key=lambda i: abs(cols[i] - var0))
    ax.axvline(vcol, color="#f59e0b", lw=1.0, ls="--")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"{len(cols)} CpG (±{WINDOW}bp; orange=variant)", fontsize=7)
    ax.set_title(f"{ann['nearest_gene']} {chrom}:{pos}\nΔβ={ann['delta']:+.2f} ARI={ann['ari']:.2f} "
                 f"plcb={ann['placebo_ari']} nCpG={ann['n_paired_cpg']} [{ann['tier']}]",
                 fontsize=7.5)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in RUNDIR:
        print(f"usage: {sys.argv[0]} <{'|'.join(RUNDIR)}> [relax|strict]")
        sys.exit(2)
    sample = sys.argv[1]
    tier = sys.argv[2] if len(sys.argv) > 2 else "relax"
    d = json.load(open(f"{CS}/{sample}_credible_discovery.json"))
    cred = [r for r in d["credible_loci"] if r["pass_tierA"]]
    if tier == "strict":
        cred = [r for r in cred if r.get("tier") == "strict_ge100"]
    cred.sort(key=lambda r: -r["ari"])
    print(f"[67] {sample} ({CANCER[sample]}): {len(cred)} pass_tierA loci (tier={tier})")

    bam = pysam.AlignmentFile(bam_path(sample), "rb")
    figs = []
    for i, r in enumerate(cred):
        chrom, pos = r["locus"].split(":")
        fb = render_locus(bam, chrom, int(pos), r["axis"], r)
        if fb:
            figs.append(dict(locus=r["locus"], gene=r["nearest_gene"], axis=r["axis"],
                             delta=r["delta"], ari=r["ari"], placebo_ari=r.get("placebo_ari"),
                             n_paired_cpg=r["n_paired_cpg"], cpg_context=r.get("cpg_context"),
                             tier=r.get("tier"), fig=fb))
        if (i + 1) % 10 == 0:
            print(f"     ...{i+1}/{len(cred)} rendered={len(figs)}")
    bam.close()
    out = dict(sample=sample, cancer=CANCER[sample], tier_filter=tier,
               n_credible=len(cred), n_rendered=len(figs), figs=figs)
    outp = f"{CS}/{sample}_locus_figs.json"
    with open(outp, "w") as f:
        json.dump(out, f)
    print(f"[67] wrote {outp} ({os.path.getsize(outp)//1024} KB, {len(figs)} figures)")


if __name__ == "__main__":
    main()
