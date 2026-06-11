#!/usr/bin/env python3
"""
97 - L5 visual panels for the independent-validation replication: for exemplar loci,
render a 2x2 grid of read×CpG methylation restricted to EVEN CpGs / ODD CpGs /
FORWARD reads / REVERSE reads (each grouped HP1 above HP2, HP sidebar), so the user
can SEE whether the HP-methylation separation replicates across independent CpG halves
and strands. Edge cases (gate-confident but replication-failing) become visible.

Output: display_v2/figs_val/<key>_valpanel.png  (one per exemplar)
"""
import os, csv, glob
os.environ.setdefault("TMPDIR", "/big7_disk/liaoyoyo2001/tmp")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
SCAN = "/big7_disk/liaoyoyo2001/ism_display_scan"
OUT = f"{DV}/figs_val"
os.makedirs(OUT, exist_ok=True)
HP_COL = {1: "#1d4ed8", 2: "#9333ea"}

# (key, cls, label)
EXEMPLARS = [
    ("chr16_3444295", "tp", "Tier A 乾淨 — 四格全一致 (鐵證)"),
    ("chr20_22561507", "tp", "Tier A 邊緣 — 正股無分離/反股有 (複製失敗→肉眼)"),
    ("chr17_47956739", "tp", "Δβ-only — 真平均位移 (位移模式)"),
    ("chr13_32740726", "tp", "none — 無法獨立驗證 (稀疏)"),
]


def region_dir(cls, key):
    chrom, pos = key.rsplit("_", 1)
    for d in glob.glob(f"{SCAN}/HCC1395_{cls}/curated_{cls}/{chrom}/{chrom}_{pos}/*"):
        if os.path.exists(f"{d}/methylation/methylation.csv"):
            return d
    return None


def load(rd):
    rows = list(csv.reader(open(f"{rd}/methylation/methylation.csv")))
    rids = [int(r[0]) for r in rows[1:]]
    M = np.array([[np.nan if x == "NA" else float(x) for x in r[1:]] for r in rows[1:]])
    base, strand = {}, {}
    for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t"):
        h = r.get("hp", "")
        base[int(r["read_id"])] = 1 if h in ("1", "1-1") else (2 if h in ("2", "2-1") else 0)
        strand[int(r["read_id"])] = r.get("strand", "")
    lab = np.array([base.get(i, 0) for i in rids])
    st = np.array([strand.get(i, "") for i in rids])
    return M, lab, st


def Fstat(M, lab):
    """between/within mean |Δβ| ratio on per-read mean over the given CpG cols."""
    rm = np.nanmean(M, axis=1)
    a, b = rm[lab == 1], rm[lab == 2]
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    within = (np.mean(np.abs(a[:, None] - a)) + np.mean(np.abs(b[:, None] - b))) / 2
    between = np.mean(np.abs(a[:, None] - b))
    return between / within if within > 0 else np.nan


def sub_heat(ax, M, lab, title):
    keep = lab > 0
    M, lab = M[keep], lab[keep]
    order = np.concatenate([np.where(lab == 1)[0], np.where(lab == 2)[0]])
    if len(order) < 2 or M.shape[1] < 1:
        ax.text(0.5, 0.5, "n/a", ha="center", va="center", fontsize=8); ax.axis("off"); return
    Mo = M[order]; lo = lab[order]
    # HP sidebar
    side = np.array([to_rgb(HP_COL[v]) for v in lo])[:, None, :]
    ax.imshow(side, aspect="auto", extent=[-0.06 * M.shape[1], 0, len(order), 0])
    cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#d1d5db")
    ax.imshow(Mo, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest",
              extent=[0, M.shape[1], len(order), 0])
    n1 = int((lo == 1).sum())
    ax.axhline(n1, color="black", lw=0.8)
    ax.set_xlim(-0.06 * M.shape[1], M.shape[1]); ax.set_ylim(len(order), 0)
    ax.set_xticks([]); ax.set_yticks([])
    f = Fstat(M, lab)
    ax.set_title(f"{title}  F={f:.2f}" if np.isfinite(f) else f"{title}  F=n/a", fontsize=8)


def main():
    for key, cls, lab in EXEMPLARS:
        rd = region_dir(cls, key)
        if rd is None:
            print(f"[97] {key} no region"); continue
        M, hp, st = load(rd)
        ncol = M.shape[1]
        ev = np.arange(ncol) % 2 == 0
        fig, ax = plt.subplots(2, 2, figsize=(8.2, 6.2))
        sub_heat(ax[0, 0], M[:, ev], hp, "偶數 CpG (split-half A)")
        sub_heat(ax[0, 1], M[:, ~ev], hp, "奇數 CpG (split-half B)")
        sub_heat(ax[1, 0], M[st == "+"], hp[st == "+"], "正股 reads (forward)")
        sub_heat(ax[1, 1], M[st == "-"], hp[st == "-"], "反股 reads (reverse)")
        fig.suptitle(f"{key.replace('_',':')}  —  {lab}\n(HP1=藍 上 / HP2=紫 下；紅=甲基 藍=未甲基；F>1 且四格一致 = 真且可複製)", fontsize=9.5)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fig.savefig(f"{OUT}/{key}_valpanel.png", dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"[97] rendered {key}")
    print(f"[97] panels -> {OUT}")


if __name__ == "__main__":
    main()
