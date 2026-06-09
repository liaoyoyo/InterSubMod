#!/usr/bin/env python3
"""
101 - Render BOTH heatmaps (methylation read×CpG + ISM read×read NHD distance) for
EVERY one of the 30,350 ISM-analysed loci, from the full ISM scan (script 100). Same
render logic as 85. Output PNGs to figs/ + fig_index.json (key -> has_fig + which class).

Classification fields stay sourced from the frozen existence-scan summary (downstream
builder 102); this only produces the figures.
"""
import os, csv, json, glob
os.environ.setdefault("TMPDIR", "/big7_disk/liaoyoyo2001/tmp")
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
SCAN = "/big7_disk/liaoyoyo2001/ism_full_scan"
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
FIGDIR = f"{DV}/figs"
os.makedirs(FIGDIR, exist_ok=True)

HP_ORDER = ["1", "1-1", "2", "2-1"]
HP_COL = {"1": "#1d4ed8", "1-1": "#16a34a", "2": "#9333ea", "2-1": "#ca8a04"}
ALT_COL = {"ALT": "#dc2626", "REF": "#0ea5e9"}
DPI = 96


def num(r, k):
    try:
        v = float(r.get(k, "")); return None if v != v else v
    except (ValueError, TypeError):
        return None


def region_dir(cls, chrom, pos):
    for d in glob.glob(f"{SCAN}/HCC1395_{cls}/filtered_snv_{cls}/{chrom}/{chrom}_{pos}/*"):
        if os.path.exists(f"{d}/distance/NHD/matrix.csv"):
            return d
    return None


def load_reads(rd):
    reads = {}
    with open(f"{rd}/reads/reads.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            reads[int(r["read_id"])] = dict(name=r["read_name"], hp=r.get("hp", ""),
                                            alt=r.get("alt_support", ""), strand=r.get("strand", ""))
    return reads


def hp_key(hp):
    return (HP_ORDER.index(hp) if hp in HP_ORDER else len(HP_ORDER), hp)


def sidebar(ax, vals, cmap, horiz=False):
    rgb = np.array([to_rgb(cmap.get(v, "#e5e7eb")) for v in vals])
    img = rgb[np.newaxis, :, :] if horiz else rgb[:, np.newaxis, :]
    ax.imshow(img, aspect="auto", interpolation="nearest"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def render_meth(rd, reads, pos, title, outpath):
    rows = list(csv.reader(open(f"{rd}/methylation/methylation.csv")))
    cpgs = [int(x) for x in rows[0][1:]]
    data = {int(r[0]): np.array([np.nan if x == "NA" else float(x) for x in r[1:]]) for r in rows[1:]}
    rids = [i for i in data if i in reads]
    rids.sort(key=lambda i: (hp_key(reads[i]["hp"]),
                             -np.nanmean(data[i]) if np.isfinite(np.nanmean(data[i])) else 0))
    M = np.vstack([data[i] for i in rids])
    var0 = pos - 1
    fig, axes = plt.subplots(1, 3, figsize=(5.0, 3.1),
                             gridspec_kw=dict(width_ratios=[0.05, 0.05, 1], wspace=0.04))
    sidebar(axes[0], [reads[i]["hp"] for i in rids], HP_COL); axes[0].set_title("HP", fontsize=6)
    sidebar(axes[1], [reads[i]["alt"] for i in rids], ALT_COL); axes[1].set_title("Allele", fontsize=6)
    ax = axes[2]; cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#d1d5db")
    ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    y = 0
    for g in HP_ORDER + sorted(set(reads[i]["hp"] for i in rids) - set(HP_ORDER)):
        n = sum(1 for i in rids if reads[i]["hp"] == g)
        if n:
            ax.axhline(y - 0.5, color="black", lw=0.6); y += n
    if cpgs:
        ax.axvline(min(range(len(cpgs)), key=lambda k: abs(cpgs[k] - var0)), color="#f59e0b", lw=1.1, ls="--")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(f"{len(cpgs)} CpG", fontsize=6.5)
    fig.suptitle(title, fontsize=7.5, y=0.99)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight"); plt.close(fig)


def render_dist(rd, reads, title, outpath):
    rows = list(csv.reader(open(f"{rd}/distance/NHD/matrix.csv")))
    rids = [int(r[0]) for r in rows[1:]]
    idx = {rid: k for k, rid in enumerate(rids)}
    D = np.array([[float(x) for x in r[1:]] for r in rows[1:]])
    name2id = {reads[i]["name"]: i for i in reads}
    order = []
    lo = f"{rd}/clustering/leaf_order.txt"
    if os.path.exists(lo):
        for nm in open(lo).read().split():
            i = name2id.get(nm)
            if i is not None and i in idx:
                order.append(i)
    for i in rids:
        if i not in order:
            order.append(i)
    perm = [idx[i] for i in order]; Dp = D[np.ix_(perm, perm)]
    fig, axes = plt.subplots(2, 2, figsize=(3.7, 3.7),
                             gridspec_kw=dict(width_ratios=[0.06, 1], height_ratios=[0.06, 1], wspace=0.03, hspace=0.03))
    axes[0, 0].axis("off")
    sidebar(axes[0, 1], [reads[i]["hp"] for i in order], HP_COL, horiz=True)
    sidebar(axes[1, 0], [reads[i]["hp"] for i in order], HP_COL, horiz=False)
    ax = axes[1, 1]
    ax.imshow(Dp, aspect="auto", cmap="magma_r", vmin=0, vmax=max(0.6, Dp.max()), interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("reads (dendrogram order)", fontsize=6.5)
    fig.suptitle(title, fontsize=7.5, y=1.0)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight"); plt.close(fig)


def work(rec):
    cls, chrom, pos = rec
    key = f"{chrom}_{pos}"
    mp, dp = f"{FIGDIR}/{key}_meth.png", f"{FIGDIR}/{key}_dist.png"
    if os.path.exists(mp) and os.path.exists(dp):
        return key, 1   # already rendered
    rd = region_dir(cls, chrom, pos)
    if rd is None:
        return key, 0
    try:
        reads = load_reads(rd)
        if len(reads) < 6:
            return key, 0
        t = f"{chrom}:{pos}  reads={len(reads)}"
        render_meth(rd, reads, int(pos), t, mp)
        render_dist(rd, reads, t, dp)
        return key, 1
    except Exception:
        return key, 0


def main():
    loci = []
    for cls in ("tp", "fp"):
        for r in csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")):
            loci.append((cls, r["Chr"], r["Pos"]))
    print(f"[101] rendering {len(loci)} loci x2 ...")
    idx = {}
    with Pool(12) as p:
        for k, (key, ok) in enumerate(p.imap_unordered(work, loci, chunksize=24)):
            idx[key] = ok
            if (k + 1) % 2000 == 0:
                done = sum(idx.values())
                print(f"   ...{k+1}/{len(loci)}  with_fig={done}")
    with open(f"{DV}/fig_index.json", "w") as f:
        json.dump(idx, f)
    nfig = sum(idx.values())
    npng = len(glob.glob(f"{FIGDIR}/*.png"))
    print(f"[101] done: {nfig}/{len(loci)} loci with figs; {npng} PNGs; fig_index.json written")


if __name__ == "__main__":
    main()
