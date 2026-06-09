#!/usr/bin/env python3
"""
85 - Render TWO heatmaps per curated locus, both directly from the ISM per-region
output so methylation matrix, read-read distance matrix and CramersV are all from
the SAME ISM computation:

  <key>_meth.png : read x CpG methylation matrix (rows grouped by HP sub-haplotype,
                   then sorted by mean methylation; HP + ALT/REF colour sidebar)
  <key>_dist.png : read x read NHD distance matrix, reordered by the ISM dendrogram
                   leaf order so cluster block-structure is visible; HP + ALT sidebars
                   on both axes -> lets the user SEE whether reads cluster by
                   haplotype/allele (i.e. whether the CramersV number is real).

Outputs under .../cross_sample/display_v2/:
  figs/<chr>_<pos>_meth.png, figs/<chr>_<pos>_dist.png
  manifest.json  (curated fields + per-region ISM stats + which figs exist)
"""
import os, csv, json, glob
os.environ.setdefault("TMPDIR", "/big7_disk/liaoyoyo2001/tmp")
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
DV = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample/display_v2"
SCAN = "/big7_disk/liaoyoyo2001/ism_display_scan"
FIGDIR = f"{DV}/figs"
os.makedirs(FIGDIR, exist_ok=True)

HP_ORDER = ["1", "1-1", "2", "2-1"]
HP_COL = {"1": "#1d4ed8", "1-1": "#16a34a", "2": "#9333ea", "2-1": "#ca8a04"}
ALT_COL = {"ALT": "#dc2626", "REF": "#0ea5e9"}
DPI = 104


def region_dir(cls, chrom, pos):
    pat = f"{SCAN}/HCC1395_{cls}/curated_{cls}/{chrom}/{chrom}_{pos}/*"
    for d in glob.glob(pat):
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
    """Draw a thin categorical colour strip aligned to the heatmap."""
    arr = np.array([cmap.get(v, "#e5e7eb") for v in vals])
    from matplotlib.colors import to_rgb
    rgb = np.array([to_rgb(c) for c in arr])
    img = rgb[np.newaxis, :, :] if horiz else rgb[:, np.newaxis, :]
    ax.imshow(img, aspect="auto", interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def render_meth(rd, reads, chrom, pos, title, outpath):
    rows = list(csv.reader(open(f"{rd}/methylation/methylation.csv")))
    header = rows[0]
    cpgs = [int(x) for x in header[1:]]
    data = {}
    for r in rows[1:]:
        rid = int(r[0])
        vals = [np.nan if x == "NA" else float(x) for x in r[1:]]
        data[rid] = np.array(vals)
    rids = [rid for rid in data if rid in reads]
    # order: HP group, then mean methylation desc
    rids.sort(key=lambda i: (hp_key(reads[i]["hp"]),
                             -np.nanmean(data[i]) if np.isfinite(np.nanmean(data[i])) else 0))
    M = np.vstack([data[i] for i in rids])
    var0 = pos - 1
    fig, axes = plt.subplots(1, 3, figsize=(5.0, 3.1),
                             gridspec_kw=dict(width_ratios=[0.05, 0.05, 1], wspace=0.04))
    sidebar(axes[0], [reads[i]["hp"] for i in rids], HP_COL)
    axes[0].set_title("HP", fontsize=6)
    sidebar(axes[1], [reads[i]["alt"] for i in rids], ALT_COL)
    axes[1].set_title("Allele", fontsize=6)
    ax = axes[2]
    cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#d1d5db")
    ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    # HP group separators
    y = 0
    for g in HP_ORDER + sorted(set(reads[i]["hp"] for i in rids) - set(HP_ORDER)):
        n = sum(1 for i in rids if reads[i]["hp"] == g)
        if n:
            ax.axhline(y - 0.5, color="black", lw=0.6); y += n
    if cpgs:
        vcol = min(range(len(cpgs)), key=lambda k: abs(cpgs[k] - var0))
        ax.axvline(vcol, color="#f59e0b", lw=1.1, ls="--")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"{len(cpgs)} CpG", fontsize=6.5)
    fig.suptitle(title, fontsize=7.5, y=0.99)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight"); plt.close(fig)


def render_dist(rd, reads, title, outpath):
    rows = list(csv.reader(open(f"{rd}/distance/NHD/matrix.csv")))
    rids = [int(r[0]) for r in rows[1:]]
    idx = {rid: k for k, rid in enumerate(rids)}
    D = np.array([[float(x) for x in r[1:]] for r in rows[1:]])
    # reorder by dendrogram leaf order (read_name -> read_id)
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
    perm = [idx[i] for i in order]
    Dp = D[np.ix_(perm, perm)]
    fig, axes = plt.subplots(2, 2, figsize=(3.7, 3.7),
                             gridspec_kw=dict(width_ratios=[0.06, 1], height_ratios=[0.06, 1],
                                              wspace=0.03, hspace=0.03))
    axes[0, 0].axis("off")
    sidebar(axes[0, 1], [reads[i]["hp"] for i in order], HP_COL, horiz=True)
    sidebar(axes[1, 0], [reads[i]["hp"] for i in order], HP_COL, horiz=False)
    ax = axes[1, 1]
    im = ax.imshow(Dp, aspect="auto", cmap="magma_r", vmin=0, vmax=max(0.6, Dp.max()),
                   interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("reads (dendrogram order)", fontsize=6.5)
    fig.suptitle(title, fontsize=7.5, y=1.0)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight"); plt.close(fig)


def region_stats(rd):
    j = f"{rd}/clustering/significance.json"
    if not os.path.exists(j):
        return {}
    s = json.load(open(j))
    return dict(
        optimal_k=s.get("optimal_k"),
        cv_alt=round(s.get("global_alt", {}).get("cramers_v", 0), 3),
        cv_hp=round(s.get("global_hp", {}).get("cramers_v", 0), 3),
        cv_hpfine=round(s.get("global_hp_fine", {}).get("cramers_v", 0), 3),
        permanova_f=round(s.get("cluster_structure", {}).get("permanova_f", 0), 1),
        dispersion_warn=bool(s.get("cluster_structure", {}).get("dispersion_warning", False)),
    )


def work(d):
    chrom, pos, cls = d["chr"], d["pos"], d["cls"]
    key = f"{chrom}_{pos}"
    rd = region_dir(cls, chrom, pos)
    if rd is None:
        return key, dict(d, has_fig=False, note="no ISM region (low reads / not emitted)")
    try:
        reads = load_reads(rd)
        if len(reads) < 6:
            return key, dict(d, has_fig=False, note=f"only {len(reads)} reads")
        t = (f"{chrom}:{pos}  CV={d['cv_max']:.2f}  Δβ={d['db']:+.2f}\n"
             f"reads={d['reads']} {'LOH' if d['loh'] else ''} {d['category']}")
        render_meth(rd, reads, chrom, pos, t, f"{FIGDIR}/{key}_meth.png")
        render_dist(rd, reads, t, f"{FIGDIR}/{key}_dist.png")
        st = region_stats(rd)
        return key, dict(d, has_fig=True, **st)
    except Exception as e:
        return key, dict(d, has_fig=False, note=f"render error: {type(e).__name__}: {e}")


def main():
    curated = json.load(open(f"{DV}/curated_loci.json"))
    print(f"[85] rendering {len(curated)} loci x2 figs ...")
    manifest = {}
    with Pool(12) as p:
        for k, (key, rec) in enumerate(p.imap_unordered(work, curated, chunksize=8)):
            manifest[key] = rec
            if (k + 1) % 200 == 0:
                ok = sum(1 for v in manifest.values() if v.get("has_fig"))
                print(f"   ...{k+1}/{len(curated)}  rendered={ok}")
    ok = sum(1 for v in manifest.values() if v.get("has_fig"))
    with open(f"{DV}/manifest.json", "w") as f:
        json.dump(list(manifest.values()), f)
    npng = len(glob.glob(f"{FIGDIR}/*.png"))
    print(f"[85] done: {ok}/{len(curated)} loci with figs; {npng} PNGs; manifest -> {DV}/manifest.json")


if __name__ == "__main__":
    main()
