#!/usr/bin/env python3
"""Generate per-case inspection figures for manual visual verification (Part 2C).

For each selected region: methylation heatmap (read x CpG) + distance heatmap (read x read),
reads sorted into HP groups (HP1 germline / HP1-1 carrier / HP2 germline / HP2-1 carrier),
with a HP+tumor/normal sidebar, and method values in the title. Lets the eye check whether
the HP labels / clustering / Δβ claims match the visible methylation structure.

All numbers in titles are read from the per-region significance_summary (this run). No fabrication.
"""
import csv
import glob
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

CASES_JSON = sys.argv[1]
OUT_ROOT = sys.argv[2]  # /tmp/dbeta_cases/out/cases/chr1
SUMMARY = sys.argv[3]   # /tmp/dbeta_cases/out/significance_summary.csv
FIG_DIR = sys.argv[4]
os.makedirs(FIG_DIR, exist_ok=True)

# HP group order + colors
HP_ORDER = {"1": 0, "HP1": 0, "1-1": 1, "2": 2, "HP2": 2, "2-1": 3}
HP_NAME = {0: "HP1(germ)", 1: "HP1-1(carrier)", 2: "HP2(germ)", 3: "HP2-1(carrier)", 4: "unlabeled"}
HP_COLOR = {0: "#1f5fd0", 1: "#7fb0f0", 2: "#d02020", 3: "#f08080", 4: "#bbbbbb"}


def load_summary(path):
    d = {}
    for r in csv.DictReader(open(path)):
        d[r["Pos"]] = r
    return d


def region_inner(pos):
    g = glob.glob(os.path.join(OUT_ROOT, f"chr1_{pos}", "chr*"))
    return g[0] if g else None


def load_reads(inner):
    meta = {}
    with open(os.path.join(inner, "reads", "reads.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            meta[row["read_id"]] = (row["hp"], row["is_tumor"] == "1")
    return meta


def load_meth(inner):
    rows, ids = [], []
    with open(os.path.join(inner, "methylation", "methylation.csv")) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            ids.append(p[0])
            rows.append([np.nan if (v == "NA" or v == "") else float(v) for v in p[1:]])
    return ids, np.array(rows, dtype=float)


def load_dist(inner):
    dm = os.path.join(inner, "distance", "BERNOULLI", "matrix.csv")
    if not os.path.exists(dm):
        return None, None
    ids, rows = [], []
    with open(dm) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            ids.append(p[0])
            rows.append([np.nan if (v == "NA" or v == "") else float(v) for v in p[1:]])
    return ids, np.array(rows, dtype=float)


def sort_order(read_ids, meta):
    """Return indices sorting reads by HP group then tumor(0)/normal(1)."""
    def key(i):
        rid = read_ids[i]
        hp, is_t = meta.get(rid, ("", False))
        return (HP_ORDER.get(hp, 4), 0 if is_t else 1)
    return sorted(range(len(read_ids)), key=key)


TN_COLOR = {True: "#e8820e", False: "#2ca02c"}  # tumor=orange, normal=green


def hp_sidebar(ax, read_ids, order, meta, n):
    arr = np.zeros((n, 1, 3))
    for row, i in enumerate(order):
        hp, _ = meta.get(read_ids[i], ("", False))
        arr[row, 0] = matplotlib.colors.to_rgb(HP_COLOR[HP_ORDER.get(hp, 4)])
    ax.imshow(arr, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("HP", fontsize=7)


def tn_sidebar(ax, read_ids, order, meta, n):
    arr = np.zeros((n, 1, 3))
    for row, i in enumerate(order):
        _, is_t = meta.get(read_ids[i], ("", False))
        arr[row, 0] = matplotlib.colors.to_rgb(TN_COLOR[is_t])
    ax.imshow(arr, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("T/N", fontsize=7)


def gen_figure(pos, label, srow):
    inner = region_inner(pos)
    if not inner:
        return None
    meta = load_reads(inner)
    mids, M = load_meth(inner)
    dids, D = load_dist(inner)
    if M.size == 0:
        return None
    order = sort_order(mids, meta)
    Msorted = M[order, :]
    # distance reordered to meth read order (align by read_id)
    Dsorted = None
    if D is not None and dids == mids:
        idx = np.array(order)
        Dsorted = D[np.ix_(idx, idx)]

    n = len(mids)
    fig = plt.figure(figsize=(16, 6.5))
    gs = fig.add_gridspec(1, 7, width_ratios=[0.08, 0.08, 2.2, 0.18, 0.08, 0.08, 2.2], wspace=0.04)

    # methylation sidebars (HP + T/N) + heatmap
    hp_sidebar(fig.add_subplot(gs[0, 0]), mids, order, meta, n)
    tn_sidebar(fig.add_subplot(gs[0, 1]), mids, order, meta, n)
    axm = fig.add_subplot(gs[0, 2])
    cmap = plt.cm.RdYlBu_r.copy()
    cmap.set_bad("#dddddd")
    im = axm.imshow(Msorted, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    axm.set_title(f"Methylation β (read×CpG), reads grouped by HP  [n={n}]", fontsize=9)
    axm.set_xlabel("CpG (genomic order)", fontsize=8)
    axm.set_ylabel("reads (HP-grouped)", fontsize=8)
    # HP group separators
    grp = [HP_ORDER.get(meta.get(mids[i], ("", 0))[0], 4) for i in order]
    for r in range(1, n):
        if grp[r] != grp[r - 1]:
            axm.axhline(r - 0.5, color="black", lw=0.8)
            if Dsorted is not None:
                pass
    plt.colorbar(im, ax=axm, fraction=0.025, pad=0.01, label="β")

    # distance sidebars + heatmap
    if Dsorted is not None:
        hp_sidebar(fig.add_subplot(gs[0, 4]), mids, order, meta, n)
        tn_sidebar(fig.add_subplot(gs[0, 5]), mids, order, meta, n)
        axd = fig.add_subplot(gs[0, 6])
        imd = axd.imshow(Dsorted, aspect="auto", cmap="viridis", interpolation="nearest")
        axd.set_title("Read×Read distance (BERNOULLI), same HP order", fontsize=9)
        axd.set_xlabel("reads (HP-grouped)", fontsize=8)
        for r in range(1, n):
            if grp[r] != grp[r - 1]:
                axd.axhline(r - 0.5, color="white", lw=0.6)
                axd.axvline(r - 0.5, color="white", lw=0.6)
        plt.colorbar(imd, ax=axd, fraction=0.025, pad=0.01, label="dist")
    else:
        axd = fig.add_subplot(gs[0, 6])
        axd.text(0.5, 0.5, "distance matrix\nnot available", ha="center", va="center")
        axd.axis("off")

    # title with method values
    def gv(k):
        v = srow.get(k, "NA")
        return v if v not in ("", None) else "NA"
    title = (
        f"[{label}] chr1:{pos}  |  NumReads={gv('NumReads')} NReadsValid={gv('NReadsValid')}  "
        f"Tumor HP1/HP2={gv('Tumor_HP1')}/{gv('Tumor_HP2')}  Normal HP1/HP2={gv('Normal_HP1')}/{gv('Normal_HP2')}\n"
        f"CramersV={gv('CramersV')} GlobalP={gv('GlobalP')} PERMANOVA F={gv('ClusterPermanovaF')} p={gv('ClusterPermanovaP')}  "
        f"VC={gv('VerificationClass')}  |  "
        f"germΔβ={gv('GermlineAsmDbeta')}({gv('GermlineAsmDbeta_Sig')}) "
        f"subHP1={gv('SubcloneDbeta_HP1')}({gv('SubcloneDbeta_HP1_Sig')}) "
        f"subHP2={gv('SubcloneDbeta_HP2')}({gv('SubcloneDbeta_HP2_Sig')}) "
        f"somΔβ={gv('SomaticResidualDbeta')}({gv('SomaticResidualDbeta_Sig')})"
    )
    fig.suptitle(title, fontsize=8.5, y=0.99)

    # legend (HP groups + tumor/normal)
    handles = [plt.Rectangle((0, 0), 1, 1, color=HP_COLOR[k]) for k in range(5)]
    handles += [plt.Rectangle((0, 0), 1, 1, color=TN_COLOR[True]), plt.Rectangle((0, 0), 1, 1, color=TN_COLOR[False])]
    labels = [HP_NAME[k] for k in range(5)] + ["tumor", "normal"]
    fig.legend(handles, labels, loc="lower center", ncol=7, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    out = os.path.join(FIG_DIR, f"case_{label}_chr1_{pos}.png")
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return out


summ = load_summary(SUMMARY)
cases = json.load(open(CASES_JSON))["cases"]
made = []
for layer, items in cases.items():
    for c in items:
        pos = c["pos"]
        srow = summ.get(pos, {})
        out = gen_figure(pos, layer, srow)
        if out:
            made.append({"layer": layer, "pos": pos, "fig": out})

print(json.dumps({"n_figures": len(made), "figures": made}, indent=2, ensure_ascii=False))
