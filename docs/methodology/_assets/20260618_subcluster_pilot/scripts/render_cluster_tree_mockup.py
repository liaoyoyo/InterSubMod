#!/usr/bin/env python3
"""[演化樹配色 mockup] 同一 region 並列: 舊(扁平 CLUSTER_COL) vs 新(tree-aware magenta-depth)。
每panel: 樹枝依該方案著色 + 側欄[coarse|HP]。HP 用 canonical(撞色檢查)。挑深樹 region 顯示階層。"""
import os, sys, csv, glob, json
import numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import cluster_tree_palette as CT
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
FIGS = "/tmp/cluster_tree_mockup"; os.makedirs(FIGS, exist_ok=True)
ODS = [f"{WT}/output/_geomdiv_TP", f"{WT}/output/_geomdiv_FP"]


def find_region(chrom, rec_pos):
    for od in ODS:
        h = glob.glob(f"{od}/**/{chrom}_{rec_pos}_*/clustering/phylo_groups.tsv", recursive=True)
        if h:
            return os.path.dirname(os.path.dirname(h[0]))
    return None


def branch_colors(Z, n, clab, colorfn):
    """每 link 若子樹單一 coarse 標籤 → 該標籤色, 否則灰。"""
    desc = {i: [i] for i in range(n)}
    for i in range(len(Z)):
        a, b = int(Z[i, 0]), int(Z[i, 1]); desc[n + i] = desc[a] + desc[b]
    def lc(k):
        s = {clab[i] for i in desc[k]}
        return colorfn(list(s)[0]) if len(s) == 1 else "#666"
    return lc


def render(chrom, rec_pos):
    region = find_region(chrom, rec_pos)
    if not region:
        return None
    rows = list(csv.DictReader(open(f"{region}/clustering/phylo_groups.tsv"), delimiter="\t"))
    if len(rows) < 8:
        return None
    rid = [int(r["read_id"]) for r in rows]; clab = [r["coarse_label"] for r in rows]; hp = [r["hp"] for r in rows]
    dids, D = CR.loadm(f"{region}/distance/BERNOULLI/matrix.csv")
    idx = [dids.index(str(x)) for x in rid if str(x) in dids]
    sub = D[np.ix_(idx, idx)]; n = len(rid)
    Z, _ = CR.linkZ(sub)
    order = dendrogram(Z, orientation="left", no_plot=True)["leaves"][::-1]
    co = [clab[i] for i in order]; hpo = [hp[i] for i in order]
    # 舊扁平: distinct label → CLUSTER_COL index
    terms = sorted(set(l for l in clab if l not in ("other", "outlier", "-")))
    oldmap = {t: H.CLUSTER_COL[i % len(H.CLUSTER_COL)] for i, t in enumerate(terms)}
    oldfn = lambda l: "#666" if l == "other" else "#ccc" if l in ("outlier", "-") else oldmap.get(l, "#888")
    newfn = CT.cluster_tree_color
    hpcol = [H.HP_COL.get(h, H.NA_GREY) for h in hpo]

    fig = plt.figure(figsize=(12.5, 5.6))
    gs = fig.add_gridspec(1, 7, width_ratios=[1.0, 0.07, 0.07, 0.18, 1.0, 0.07, 0.07], wspace=0.05)
    for col0, title, cfn in [(0, "舊 · 扁平 CLUSTER_COL", oldfn), (4, "新 · tree-aware (magenta-depth)", newfn)]:
        ax = fig.add_subplot(gs[0, col0])
        lc = branch_colors(Z, n, clab, cfn)
        dendrogram(Z, orientation="left", no_labels=True, ax=ax, link_color_func=lc)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9); [ax.spines[s].set_visible(False) for s in ax.spines]
        H._sb(fig.add_subplot(gs[0, col0 + 1]), [cfn(c) for c in co], "coarse")
        H._sb(fig.add_subplot(gs[0, col0 + 2]), hpcol, "HP")
    ng = len(terms)
    fig.suptitle(f"{chrom}:{rec_pos+5000}  coarse {ng}群 標籤{terms}  (左舊扁平 vs 右新tree-aware; HP欄不動=撞色檢查)", fontsize=10, y=1.02)
    fn = f"{FIGS}/treecolor_{chrom}_{rec_pos}.png"; fig.savefig(fn, dpi=125, bbox_inches="tight"); plt.close(fig)
    return (f"{chrom}:{rec_pos+5000}", terms, fn)


if __name__ == "__main__":
    recs = json.load(open(f"{os.path.dirname(__file__)}/../geom_records.json"))
    # 挑 coarse_ng>=3 (深樹, 標籤多層) 顯示階層
    cand = sorted([r for r in recs if r["coarse_ng"] >= 3], key=lambda r: -r["coarse_ng"])[:8]
    out = []
    for r in cand:
        res = render(r["chrom"], r["rec_pos"])
        if res:
            out.append({"title": res[0], "terms": res[1], "fig": res[2]})
            print(f"  {res[0]} labels={res[1]} -> {os.path.basename(res[2])}")
        if len(out) >= 5:
            break
    json.dump(out, open(f"{FIGS}/mockup.json", "w"), ensure_ascii=False)
    print(f"DONE {len(out)} mockup figs")
