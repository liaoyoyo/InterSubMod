#!/usr/bin/env python3
"""U1 case visualization: MAX_DIST vs SKIP distance heatmap + methylation matrix.

For selected "obvious improvement" regions (Weak->Strong, CramersV 0->1, DominantLabel=hp),
plot read-read distance heatmaps (both strategies) + methylation matrix, reads sorted by HP tag,
to visually confirm SKIP recovers the germline-HP clustering block that the MAX_DIST 1.0
pollution blurs.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

OUTDIR = "docs/experiments/in_progress/2026/06/20260614_u1_maxdist_vs_skip/figures"
os.makedirs(OUTDIR, exist_ok=True)

CASES = [
    ("chr1_73791759", "chr1_73786759_73796759"),
    ("chr1_56791672", "chr1_56786672_56796672"),
    ("chr1_58306441", "chr1_58301441_58311441"),
    ("chr1_67579144", "chr1_67574144_67584144"),
]


def load(strat, rid, rsub):
    base = f"/tmp/u1_case_{strat}/u1_cases/chr1/{rid}/{rsub}"
    D = pd.read_csv(f"{base}/distance/BERNOULLI/matrix.csv", index_col=0).values.astype(float)
    M = pd.read_csv(f"{base}/methylation/methylation.csv", index_col=0)
    reads = pd.read_csv(f"{base}/reads/reads.tsv", sep="\t")
    return D, M, reads


def hp_sort_key(h):
    h = str(h)
    if h.startswith("1"):
        return (0, h)
    if h.startswith("2"):
        return (1, h)
    return (2, h)


# HP group color strip
HP_COLORS = {"1": "#d62728", "1-1": "#ff9896", "2": "#1f77b4", "2-1": "#aec7e8"}


def hp_color(h):
    return HP_COLORS.get(str(h), "#bbbbbb")


for rid, rsub in CASES:
    Dm, Mm, reads = load("max_dist", rid, rsub)
    Ds, _, _ = load("skip", rid, rsub)
    hp = reads["hp"].astype(str).values
    order = sorted(range(len(hp)), key=lambda i: hp_sort_key(hp[i]))
    Dm_s = Dm[np.ix_(order, order)]
    Ds_s = Ds[np.ix_(order, order)]
    M_s = Mm.values[order, :].astype(float)
    hp_s = [hp[i] for i in order]
    n = len(order)
    n_hp1 = sum(1 for h in hp_s if h.startswith("1"))
    n_hp2 = sum(1 for h in hp_s if h.startswith("2"))

    vir = plt.cm.viridis.copy()
    vir.set_bad("#dddddd")  # NaN (SKIP invalid pair) = light gray
    rdbu = plt.cm.RdBu_r.copy()
    rdbu.set_bad("#dddddd")

    fig, ax = plt.subplots(1, 3, figsize=(20, 6.5))
    fig.suptitle(f"U1 case {rid}  (n_reads={n}, HP1-family={n_hp1}, HP2-family={n_hp2})  reads sorted by HP tag",
                 fontsize=13)

    im0 = ax[0].imshow(Dm_s, cmap=vir, vmin=0, vmax=1)
    ax[0].set_title("MAX_DIST distance\n(invalid pair = 1.0 pollution)")
    plt.colorbar(im0, ax=ax[0], fraction=0.046)

    im1 = ax[1].imshow(np.ma.masked_invalid(Ds_s), cmap=vir, vmin=0, vmax=1)
    ax[1].set_title("SKIP distance\n(invalid pair = NaN, gray-masked)")
    plt.colorbar(im1, ax=ax[1], fraction=0.046)

    im2 = ax[2].imshow(np.ma.masked_invalid(M_s), cmap=rdbu, vmin=0, vmax=1, aspect="auto")
    ax[2].set_title("Methylation (read x CpG)\nHP-sorted")
    plt.colorbar(im2, ax=ax[2], fraction=0.046, label="methyl prob")

    # HP boundary line (between HP1 and HP2 family blocks)
    for a in (ax[0], ax[1]):
        a.axhline(n_hp1 - 0.5, color="yellow", lw=1.2, ls="--")
        a.axvline(n_hp1 - 0.5, color="yellow", lw=1.2, ls="--")
        a.set_xlabel("read (HP-sorted)")
        a.set_ylabel("read (HP-sorted)")
    ax[2].axhline(n_hp1 - 0.5, color="yellow", lw=1.2, ls="--")
    ax[2].set_xlabel("CpG site")
    ax[2].set_ylabel("read (HP-sorted)")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = f"{OUTDIR}/case_{rid}.png"
    plt.savefig(out, dpi=110, bbox_inches="tight")
    plt.close()
    print(f"saved {out}  n={n} HP1f={n_hp1} HP2f={n_hp2} "
          f"SKIP_nan_frac={np.isnan(Ds_s).mean():.2f}")
print("DONE")
