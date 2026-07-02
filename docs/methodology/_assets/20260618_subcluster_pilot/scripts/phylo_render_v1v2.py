#!/usr/bin/env python3
"""[視覺] phylo v1│v2 並排對照圖 — 變動位點，肉眼看 double-dip 假象在 v2 消失。
每位點: 左=v1(現法)演化群著色樹+標籤側欄, 右=v2(修法)同樹同甲基但 v2 標籤。中間共用甲基 read×CpG。
重用 _ws_render 快取矩陣。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
FIGS = f"{A}/figs_phylo_compare"; os.makedirs(FIGS, exist_ok=True)
MINSZ = 3; SEP_MIN = 1.3
FAM = {"1": ["#7f1d1d", "#dc2626", "#f87171", "#fca5a5"], "2": ["#1e3a8a", "#2563eb", "#60a5fa", "#93c5fd"],
       "3": ["#14532d", "#16a34a", "#4ade80"], "4": ["#713f12", "#d97706", "#fbbf24"]}
TARGETS = ["chr21_10353822", "chr20_30274614", "chr20_21855867", "chr21_40743336"]


def labcol(lab):
    if lab is None or lab in ("outlier", "-"): return "#cfcfcf"
    fam = lab.split("-")[0]; depth = lab.count("-"); pal = FAM.get(fam, ["#555", "#888", "#aaa"]); return pal[min(depth, len(pal) - 1)]


def _bw(s, S1, S2):
    bet = s[np.ix_(S1, S2)]; bet = bet[bet >= 0]
    w1 = s[np.ix_(S1, S1)][np.triu_indices(len(S1), 1)]; w1 = w1[w1 >= 0]
    w2 = s[np.ix_(S2, S2)][np.triu_indices(len(S2), 1)]; w2 = w2[w2 >= 0]
    wm = np.concatenate([w1, w2]) if (w1.size or w2.size) else np.array([])
    if bet.size == 0 or wm.size == 0 or wm.mean() <= 1e-6: return None
    return float(bet.mean()) / float(wm.mean())


def _tree(D):
    n = D.shape[0]; Z, s = CR.linkZ(D)
    desc = {i: [i] for i in range(n)}; ch = {}
    for i in range(len(Z)):
        a, b = int(Z[i, 0]), int(Z[i, 1]); desc[n + i] = desc[a] + desc[b]; ch[n + i] = (a, b)
    return Z, s, desc, ch, n


def phylo_label(sub, P, rng, mode, RNULL):
    Z, s, desc, ch, n = _tree(sub)
    gnull = None
    if mode == "v1":
        C = P.shape[1]; gnull = []
        for _ in range(RNULL):
            Pn = P.copy()
            for cc in range(C):
                col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
            Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); gnull.append(np.maximum(Dn, Dn.T))

    def split_real(node):
        if node not in ch: return False
        c1, c2 = ch[node]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
        if min(len(S1), len(S2)) < MINSZ: return False
        r = _bw(s, S1, S2)
        if r is None or r < SEP_MIN: return False
        if mode == "v1":
            ns = [_bw(Dn, S1, S2) for Dn in gnull]
        else:
            S = np.array(desc[node]); m = len(S); Psub = P[S]; ns = []
            for _ in range(RNULL):
                Pn = Psub.copy()
                for cc in range(Pn.shape[1]):
                    col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                    if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
                Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
                _, sn, dn, cn, _ = _tree(Dn); nc1, nc2 = cn[2 * m - 2]
                ns.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
        ns = [x for x in ns if x is not None]
        return r > (np.percentile(ns, 95) if ns else 0)

    lab = [None] * n

    def rec(node, label):
        leaves = desc[node]
        if node not in ch or len(leaves) < 2 * MINSZ:
            for i in leaves: lab[i] = label; return
        if split_real(node):
            c1, c2 = ch[node]; big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            rec(big, label + "-1"); rec(small, label + "-2")
        else:
            for i in leaves: lab[i] = label
    root = 2 * n - 2
    if split_real(root):
        c1, c2 = ch[root]; big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1); rec(big, "1"); rec(small, "2")
    else:
        for i in range(n): lab[i] = "1"
    lab = [l if l else "outlier" for l in lab]
    sm = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    return [("outlier" if l in sm else l) for l in lab], Z


dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd

for key in TARGETS:
    rd = dirmap.get(key)
    if not rd: print(f"skip {key}"); continue
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs = [int(c) for c in rows[0].split(",")[1:]]
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    lab1, Z = phylo_label(sub, P, np.random.default_rng(20260622), "v1", 12)
    lab2, _ = phylo_label(sub, P, np.random.default_rng(20260622), "v2", 40)
    g1 = sorted(set(l for l in lab1 if l != "outlier")); g2 = sorted(set(l for l in lab2 if l != "outlier"))
    dn = dendrogram(Z, orientation="left", no_plot=True); order = dn["leaves"][::-1]
    meth = np.array([P[i] for i in order])
    l1o = [lab1[i] for i in order]; l2o = [lab2[i] for i in order]
    hp_o = [reads[ids[i]]["hp"] for i in order]; al_o = [reads[ids[i]]["alt_support"] for i in order]
    HPC = H.HP_COL
    alc = H.ALT_COL
    fig = plt.figure(figsize=(13, 5.4))
    gs = fig.add_gridspec(1, 8, width_ratios=[1.1, 0.06, 0.06, 0.06, 1.0, 0.06, 0.06, 1.1], wspace=0.05)
    # v1 tree
    co1 = {i: lab1[i] for i in range(n)}; _, _, desc, ch, _ = _tree(sub)
    nodelab1 = {nd: (list({lab1[i] for i in desc[nd]})[0] if len({lab1[i] for i in desc[nd]}) == 1 else None) for nd in range(2 * n - 1)}
    ax0 = fig.add_subplot(gs[0, 0])
    dendrogram(Z, orientation="left", link_color_func=lambda k: labcol(nodelab1.get(k)) if nodelab1.get(k) else "#ccc", ax=ax0, no_labels=True)
    ax0.set_xticks([]); ax0.set_yticks([]); ax0.set_title(f"v1 現法: {len(g1)} 群\n{g1}", fontsize=9, color="#b91c1c"); [ax0.spines[s].set_visible(False) for s in ax0.spines]
    # sidebars: v1 label, hp, allele
    H._sb(fig.add_subplot(gs[0, 1]), [labcol(l) for l in l1o], "v1群")
    H._sb(fig.add_subplot(gs[0, 2]), [HPC.get(h, "#ddd") for h in hp_o], "HP")
    H._sb(fig.add_subplot(gs[0, 3]), [alc.get(a, "#ddd") for a in al_o], "ALT")
    # methylation
    mc, _ = H.mpl_cmaps(); axm = fig.add_subplot(gs[0, 4]); axm.imshow(meth, aspect="auto", cmap=mc, vmin=0, vmax=1, interpolation="nearest")
    sx = H.snv_fractional_x(cpgs, int(key.split('_')[1]))
    if sx is not None: axm.axvline(sx, color=H.SNV_COL, lw=2)
    axm.set_xticks([]); axm.set_yticks([]); axm.set_title("甲基 read×CpG", fontsize=9)
    # v2 sidebars + tree
    H._sb(fig.add_subplot(gs[0, 5]), [HPC.get(h, "#ddd") for h in hp_o], "HP")
    H._sb(fig.add_subplot(gs[0, 6]), [labcol(l) for l in l2o], "v2群")
    nodelab2 = {nd: (list({lab2[i] for i in desc[nd]})[0] if len({lab2[i] for i in desc[nd]}) == 1 else None) for nd in range(2 * n - 1)}
    ax7 = fig.add_subplot(gs[0, 7])
    dendrogram(Z, orientation="right", link_color_func=lambda k: labcol(nodelab2.get(k)) if nodelab2.get(k) else "#ccc", ax=ax7, no_labels=True)
    ax7.set_xticks([]); ax7.set_yticks([]); ax7.set_title(f"v2 修法: {len(g2)} 群\n{g2}", fontsize=9, color="#15803d"); [ax7.spines[s].set_visible(False) for s in ax7.spines]
    fig.suptitle(f"{key}  n={n} CpG={len(cpgs)}  —  v1 {len(g1)}群 → v2 {len(g2)}群（double-dip 假象在 v2 消失）", fontsize=11, y=1.02)
    fn = f"{FIGS}/cmp_{key}.png"; fig.savefig(fn, dpi=110, bbox_inches="tight"); plt.close(fig)
    print(f"WROTE {fn}  v1={g1} v2={g2}", flush=True)
print("DONE")
