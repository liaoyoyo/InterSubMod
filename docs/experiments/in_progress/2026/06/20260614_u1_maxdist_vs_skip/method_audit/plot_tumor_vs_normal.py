#!/usr/bin/env python3
"""HANDOFF §4 案例圖 — tumor vs normal germline-HP 距離結構對比。
每案例 3 聯：normal 距離(HP排序) | tumor 距離(HP排序) | 甲基(normal上/tumor下, HP排序)。
視覺驗證：normal 的 HP1/HP2 對角 block 清晰（高 HP-AUC）vs tumor block 鬆散（低 HP-AUC）。
AUC 標註值取自 significance_summary 的 HP_AUC_Normal/Tumor 欄（confound_stratify.json）。
read-only 讀 /tmp/u1case2_skip per-region 輸出。
"""
import os, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = "/tmp/u1case2_skip/u1case2"
OUT = "docs/experiments/in_progress/2026/06/20260614_u1_maxdist_vs_skip/method_audit/figures"
os.makedirs(OUT, exist_ok=True)

# case -> (chr, pos, normal_auc, tumor_auc, tag)
CASES = [
    ("chr2", 44942972, 0.993, 0.401, "clean diploid · HP1=18/HP2=24"),
    ("chr20", 58822325, 0.991, 0.500, "clean diploid · HP1=39/HP2=45"),
    ("chr1", 201317938, 0.983, 0.473, "clean diploid · HP1=21/HP2=46"),
    ("chr13", 44126493, 0.951, 0.463, "clean diploid · HP1=41/HP2=20"),
    ("chr7", 46568719, 1.000, 0.390, "clean diploid · HP1=10/HP2=4"),
    ("chr10", 35976901, 1.000, 0.457, "clean diploid · HP1=5/HP2=8"),
    ("chr1", 84737071, 0.901, 0.222, "CONTRAST: CNV_Gain+LOH (confound)"),
    ("chr7", 149051472, 0.989, 0.383, "CONTRAST: CNV_Gain (confound)"),
]


def hp_key(h):
    h = str(h)
    return (0, h) if h.startswith("1") else ((1, h) if h.startswith("2") else (2, h))


def load(chrom, pos):
    region_dirs = glob.glob(f"{BASE}/{chrom}/{chrom}_{pos}/*/")
    if not region_dirs:
        return None
    d = region_dirs[0]
    D = pd.read_csv(f"{d}distance/BERNOULLI/matrix.csv", index_col=0)
    reads = pd.read_csv(f"{d}reads/reads.tsv", sep="\t")
    M = pd.read_csv(f"{d}methylation/methylation.csv", index_col=0)
    return D, reads, M


def is_tum(v):
    return str(v).strip().lower() in ("1", "true", "yes", "t")


vir = plt.cm.viridis.copy(); vir.set_bad("#dddddd")
rdbu = plt.cm.RdBu_r.copy(); rdbu.set_bad("#dddddd")

for chrom, pos, n_auc, t_auc, tag in CASES:
    r = load(chrom, pos)
    if r is None:
        print(f"SKIP {chrom}:{pos} (no output)"); continue
    D, reads, M = r
    rid2hp = dict(zip(reads["read_id"].astype(str), reads["hp"].astype(str)))
    rid2tum = dict(zip(reads["read_id"].astype(str), reads["is_tumor"].map(is_tum)))
    ids = [str(x) for x in D.index]
    Dv = D.values.astype(float)

    def subset_sorted(want_tumor):
        idx = [i for i, rid in enumerate(ids) if rid2tum.get(rid, False) == want_tumor and rid in rid2hp]
        idx = sorted(idx, key=lambda i: hp_key(rid2hp[ids[i]]))
        hps = [rid2hp[ids[i]] for i in idx]
        n1 = sum(1 for h in hps if h.startswith("1"))
        return idx, n1

    n_idx, n_n1 = subset_sorted(False)
    t_idx, t_n1 = subset_sorted(True)

    fig, ax = plt.subplots(1, 3, figsize=(19, 6))
    fig.suptitle(f"{chrom}:{pos}  ·  {tag}\nNormal HP-AUC={n_auc}  vs  Tumor HP-AUC={t_auc}  "
                 f"(germline-HP 距離區分度; tumor 越低=HP 結構越被打亂)", fontsize=12)

    for col, (idx, n1, lab, auc) in enumerate([
        (n_idx, n_n1, "Normal reads", n_auc), (t_idx, t_n1, "Tumor reads", t_auc)]):
        if len(idx) >= 2:
            sub = Dv[np.ix_(idx, idx)]
            im = ax[col].imshow(np.ma.masked_invalid(sub), cmap=vir, vmin=0, vmax=1)
            ax[col].axhline(n1 - 0.5, color="yellow", lw=1.3, ls="--")
            ax[col].axvline(n1 - 0.5, color="yellow", lw=1.3, ls="--")
            plt.colorbar(im, ax=ax[col], fraction=0.046)
        ax[col].set_title(f"{lab} 距離 (HP排序, n={len(idx)})\nHP-AUC={auc}")
        ax[col].set_xlabel("read (HP-sorted)"); ax[col].set_ylabel("read")

    # methylation: normal then tumor, HP-sorted within
    order = n_idx + t_idx
    Mv = M.values.astype(float)
    # align M rows to D ids (M index = read_id)
    mid = {str(x): k for k, x in enumerate(M.index)}
    mrows = [mid[ids[i]] for i in order if ids[i] in mid]
    if mrows:
        im2 = ax[2].imshow(np.ma.masked_invalid(Mv[mrows, :]), cmap=rdbu, vmin=0, vmax=1, aspect="auto")
        ax[2].axhline(len(n_idx) - 0.5, color="lime", lw=1.6)  # normal/tumor 分界
        plt.colorbar(im2, ax=ax[2], fraction=0.046, label="methyl prob")
    ax[2].set_title("甲基 (read×CpG)\n上=Normal 下=Tumor (綠線分界), HP排序")
    ax[2].set_xlabel("CpG site"); ax[2].set_ylabel("read")

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fn = f"{OUT}/tn_{chrom}_{pos}.png"
    plt.savefig(fn, dpi=105, bbox_inches="tight"); plt.close()
    nan_t = np.isnan(Dv[np.ix_(t_idx, t_idx)]).mean() if len(t_idx) >= 2 else float("nan")
    print(f"saved {fn}  normal_n={len(n_idx)} tumor_n={len(t_idx)} tumor_nan_frac={nan_t:.2f}")
print("DONE")
