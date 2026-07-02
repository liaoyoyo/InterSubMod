#!/usr/bin/env python3
"""[每位點 tumor-vs-normal 甲基對照圖] 走保留 region dirs, 從 methylation.csv(全 read, 含 normal) + reads.tsv(is_tumor/hp)
畫「上 tumor / 下 normal」甲基 read×CpG 熱圖 + T/N + HP 側欄(canonical H 配色)。不重跑 binary。
輸出 figs_cpp_wg_full_tn/cpp_<chr>_<start>_<end>_tn.png。補 phylo 圖只有 tumor 的缺口。"""
import os, csv, glob, json, sys
import numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
KEEP = f"{WT}/output/_phylo_wg_full"
FIGS = f"{A}/figs_cpp_wg_full_tn"; os.makedirs(FIGS, exist_ok=True)


def _mb(M, rid):
    v = np.array(M[rid]); return float(np.nanmean(v)) if not np.all(np.isnan(v)) else 0.0


def render_one(region):
    try:
        key = os.path.basename(region)  # chr_start_end
        meth = f"{region}/methylation/methylation.csv"; rt = f"{region}/reads/reads.tsv"
        if not (os.path.exists(meth) and os.path.exists(rt)): return None
        mr = open(meth).read().strip().split("\n")
        if len(mr) < 2: return None
        cpgs = [int(c) for c in mr[0].split(",")[1:]]
        M = {}
        for ln in mr[1:]:
            q = ln.split(","); M[q[0]] = [np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]]
        rows = open(rt).read().splitlines(); hdr = rows[0].split("\t")
        ti = hdr.index("is_tumor"); hi = hdr.index("hp") if "hp" in hdr else None
        tum, nor = [], []
        for r in rows[1:]:
            c = r.split("\t"); rid = c[0]
            if rid not in M: continue
            hp = c[hi] if (hi is not None and len(c) > hi) else "0"
            if c[ti] == "1": tum.append((rid, hp))
            elif c[ti] == "0": nor.append((rid, hp))
        if not tum or not nor: return None  # 需 tumor+normal 才有對照意義
        tum.sort(key=lambda x: -_mb(M, x[0])); nor.sort(key=lambda x: -_mb(M, x[0]))
        order = tum + nor; tn_split = len(tum)
        mat = np.array([M[rid] for rid, _ in order])
        tnhex = [H.TN_COL['1']] * len(tum) + [H.TN_COL['0']] * len(nor)
        hphex = [H.HP_COL.get(hp, "#dddddd") for _, hp in order]
        dtn = float(np.nanmean([M[rid] for rid, _ in tum]) - np.nanmean([M[rid] for rid, _ in nor]))
        mc, _ = H.mpl_cmaps()
        fig = plt.figure(figsize=(8.6, 5.0)); gs = fig.add_gridspec(1, 3, width_ratios=[0.05, 0.05, 1.6], wspace=0.05)
        H._sb(fig.add_subplot(gs[0, 0]), tnhex, "T/N")
        H._sb(fig.add_subplot(gs[0, 1]), hphex, "HP")
        am = fig.add_subplot(gs[0, 2]); am.imshow(mat, aspect="auto", cmap=mc, vmin=0, vmax=1, interpolation="nearest")
        am.axhline(tn_split - 0.5, color="#0a0a0a", lw=1.8)  # tumor|normal 分界
        try:
            sx = H.snv_fractional_x(cpgs, (int(key.split("_")[1]) + int(key.split("_")[2])) // 2)
            if sx is not None: am.axvline(sx, color=H.SNV_COL, lw=2)
        except Exception: pass
        am.set_xticks([]); am.set_yticks([]); am.set_title("甲基 read×CpG（上 tumor / 下 normal, RdBu_r）", fontsize=8.5)
        fig.suptitle(f"{key}  tumor {len(tum)} / normal {len(nor)}   Δβ(T−N)={dtn:+.3f}", fontsize=10, y=1.01)
        fn = f"{FIGS}/cpp_{key}_tn.png"; fig.savefig(fn, dpi=118, bbox_inches="tight"); plt.close(fig)
        return key
    except Exception:
        return None


def main():
    regions = [os.path.dirname(os.path.dirname(p)) for p in glob.glob(f"{KEEP}/**/clustering/phylo_groups.tsv", recursive=True)]
    print(f"regions: {len(regions)}", flush=True)
    with Pool(24) as p:
        out = [r for r in p.map(render_one, regions, chunksize=30) if r]
    S = {"rendered": len(out), "regions": len(regions), "skipped_no_tn": len(regions) - len(out), "figdir": FIGS}
    json.dump(S, open(f"{A}/tn_figures_stats.json", "w"), indent=2)
    print("DONE " + json.dumps(S), flush=True)


if __name__ == "__main__":
    main()
