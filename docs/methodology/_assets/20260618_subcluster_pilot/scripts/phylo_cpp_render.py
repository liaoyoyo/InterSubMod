#!/usr/bin/env python3
"""[C++輸出渲染] 從 C++ binary 原生輸出的 phylo_groups.tsv 讀標籤畫圖 — Python 只讀不重算(達成目標)。
證明: binary 一次跑出 coarse/fine 標籤 → Python renderer 只讀 tsv + 矩陣 → 畫樹(C++標籤著色)+甲基+對齊。
平行: multiprocessing.Pool(優化版, 一進程多圖, 不重複啟動 python)。"""
import os, csv, glob, json, sys
import numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from multiprocessing import Pool

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS = f"{A}/figs_cpp"; os.makedirs(FIGS, exist_ok=True)
# canonical 配色取自 ism_heatmap_std(H) 防 drift: HP1 族藍(1淺/1-1深)、HP2 族紫(2淺/2-1深)、ALT 紅/REF 琥珀、群色 teal-pink
PALETTE = H.CLUSTER_COL
HPC = H.HP_COL
alc = H.ALT_COL


def render_one(region):
    """讀 C++ phylo_groups.tsv(標籤,不重算) + 矩陣 → 畫。回傳 (key, ngroups, png)。"""
    try:
        tsv = f"{region}/clustering/phylo_groups.tsv"
        if not os.path.exists(tsv): return None
        rows = list(csv.DictReader(open(tsv), delimiter="\t"))
        if len(rows) < 6: return None
        rid = [int(r["read_id"]) for r in rows]
        clab = [r["coarse_label"] for r in rows]
        flab = [r["fine_label"] for r in rows]
        hp = [r["hp"] for r in rows]; al = [r["alt_support"] for r in rows]
        # 讀全距離矩陣 → 取 tsv read_id 子集(C++ 已過濾+peel 的順序)
        dids, D = CR.loadm(f"{region}/distance/BERNOULLI/matrix.csv")
        idx = [dids.index(str(r)) if str(r) in dids else r for r in rid]
        sub = D[np.ix_(idx, idx)]
        mr = open(f"{region}/methylation/methylation.csv").read().strip().split("\n")
        cpgs = [int(c) for c in mr[0].split(",")[1:]]
        Mall = {}
        for ln in mr[1:]:
            q = ln.split(","); Mall[q[0]] = [np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]]
        meth = np.array([Mall[str(r)] for r in rid])
        n = len(rid)
        terms = sorted(set(l for l in clab if l not in ("outlier", "other")))
        # 演化樹群色 = tree-aware(H.cluster_tree_color), 取代扁平 CLUSTER_COL; other 深灰/outlier 淺灰保留可分
        def lc(l): return "#555" if l == "other" else "#dcdcdc" if l in ("outlier", "-") else H.cluster_tree_color(l)
        Z, _ = CR.linkZ(sub)
        dn = dendrogram(Z, orientation="left", no_plot=True); order = dn["leaves"][::-1]
        mo = meth[order]; do = sub[np.ix_(order, order)].copy(); np.fill_diagonal(do, 0); do[do < 0] = np.nan
        co = [clab[i] for i in order]; hpo = [hp[i] for i in order]; alo = [al[i] for i in order]
        sb = [("C++ coarse", [lc(l) for l in co]), ("HP", [HPC.get(h, "#ddd") for h in hpo]), ("ALT", [alc.get(a, "#ddd") for a in alo])]
        # 樹分支依 C++ coarse 標籤著色(與側欄一致): 子樹單一標籤=該色, 混合/含離群=灰
        descn = {i: [i] for i in range(n)}
        for i in range(len(Z)):
            aa, bb = int(Z[i, 0]), int(Z[i, 1]); descn[n + i] = descn[aa] + descn[bb]
        nodelab = {nd: (list({clab[i] for i in descn[nd]})[0] if len({clab[i] for i in descn[nd]}) == 1 else None) for nd in range(2 * n - 1)}
        mc, dc = H.mpl_cmaps(); nsb = len(sb); wr = [1.1] + [0.045] * nsb + [1.5, 0.14] + [0.045] * nsb + [1.3]
        fig = plt.figure(figsize=(12.5, 5.2)); gs = fig.add_gridspec(1, len(wr), width_ratios=wr, wspace=0.04)
        c = 0; ax = fig.add_subplot(gs[0, c]); c += 1
        dendrogram(Z, orientation="left", no_labels=True, ax=ax,
                   link_color_func=lambda k: lc(nodelab.get(k)) if nodelab.get(k) else "#888")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title("UPGMA(分支=C++ coarse群)", fontsize=8); [ax.spines[s].set_visible(False) for s in ax.spines]
        for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
        am = fig.add_subplot(gs[0, c]); c += 1; am.imshow(mo, aspect="auto", cmap=mc, vmin=0, vmax=1, interpolation="nearest")
        key = os.path.basename(region)
        try:
            sx = H.snv_fractional_x(cpgs, (int(key.split("_")[1]) + int(key.split("_")[2])) // 2)
            if sx is not None: am.axvline(sx, color=H.SNV_COL, lw=2)
        except Exception: pass
        am.set_xticks([]); am.set_yticks([]); am.set_title("甲基 read×CpG (C++標籤,Python只讀)", fontsize=8.5)
        fig.add_subplot(gs[0, c]).axis("off"); c += 1
        for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
        ad = fig.add_subplot(gs[0, c]); vmax = max(0.5, float(np.nanmax(do)) if np.isfinite(np.nanmax(do)) else 0.5)
        ad.imshow(do, aspect="auto", cmap=dc, vmin=0, vmax=vmax, interpolation="nearest"); ad.set_xticks([]); ad.set_yticks([]); ad.set_title("距離", fontsize=8.5)
        ng = len(terms); noth = sum(1 for l in clab if l == "other")
        fig.suptitle(f"[C++ native] {key}  n={n} → coarse {ng}群 標籤{terms} other{noth}  (binary→tsv→Python只畫)", fontsize=10, y=1.02)
        fn = f"{FIGS}/cpp_{key}.png"; fig.savefig(fn, dpi=130, bbox_inches="tight"); plt.close(fig)
        return (key, ng, f"figs_cpp/cpp_{key}.png")
    except Exception as e:
        return ("ERR:" + str(e)[:60], 0, None)


if __name__ == "__main__":
    od = "/big7_disk/liaoyoyo2001/tmp/_cppval2/out"
    regions = []
    for sj in glob.glob(f"{od}/**/phylo_groups_summary.json", recursive=True):
        s = json.load(open(sj))
        if s["coarse_ng"] >= 2 and s["n"] >= 60:
            regions.append(os.path.dirname(os.path.dirname(sj)))
    regions = sorted(set(regions))[:6]
    print(f"渲染 {len(regions)} 個 C++ 多群位點 (Pool 平行)")
    with Pool(min(6, len(regions))) as p:
        out = [r for r in p.map(render_one, regions) if r and r[2]]
    json.dump(out, open(f"{A}/phylo_cpp_render.json", "w"), ensure_ascii=False, indent=1)
    for k, ng, png in out: print(f"  {k}: coarse {ng}群 → {png}")
    print(f"DONE {len(out)} 圖")
