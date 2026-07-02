#!/usr/bin/env python3
"""[視覺] 全 30 位點 phylo-v2 切割+標籤觀察圖 — 供肉眼確認『切得對不對、標籤判得對不對』。
每位點: UPGMA 樹(分支依 v2 演化群著色) + 演化群側欄 + HP + ALT + strand + 甲基 read×CpG + 距離塊 + 各群對齊圖例。
v2 = per-subgroup 重分群 null + RNULL=40(已修 double-dip)。供 verify-workstation HTML 嵌圖。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
FIGS = f"{A}/figs_phylo_v2"; os.makedirs(FIGS, exist_ok=True)
MINSZ = 3; SEP_MIN = 1.3; RNULL = 40
FAM = {"1": ["#7f1d1d", "#dc2626", "#f87171", "#fca5a5"], "2": ["#1e3a8a", "#2563eb", "#60a5fa", "#93c5fd"],
       "3": ["#14532d", "#16a34a", "#4ade80"], "4": ["#713f12", "#d97706", "#fbbf24"]}


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


def phylo_v2(sub, P, rng):
    Z, s, desc, ch, n = _tree(sub)

    def split_real(node):
        if node not in ch: return False
        c1, c2 = ch[node]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
        if min(len(S1), len(S2)) < MINSZ: return False
        r = _bw(s, S1, S2)
        if r is None or r < SEP_MIN: return False
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
    return [("outlier" if l in sm else l) for l in lab], Z, desc, ch


dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd
items = json.load(open(f"{A}/ws_items.json")); out = []
HPC = H.HP_COL
alc = H.ALT_COL
for it in items:
    key = it["key"]; rd = dirmap.get(key)
    if not rd: continue
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs = [int(c) for c in rows[0].split(",")[1:]]
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids) < MINSZ * 2: continue
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    lab, Z, desc, ch = phylo_v2(sub, P, np.random.default_rng(20260622))
    labct = Counter(l for l in lab if l != "outlier"); nterm = len(labct); nout = sum(1 for l in lab if l == "outlier")

    def dom(idxs, field):
        c = Counter(reads[ids[i]][field] for i in idxs); return c.most_common(1)[0] if c else ("-", 0)
    galign = {}
    for L in labct:
        idxs = [i for i in range(n) if lab[i] == L]
        hpd = dom(idxs, "hp"); ald = dom(idxs, "alt_support")
        galign[L] = {"n": len(idxs), "hp": f"{hpd[0]}({100*hpd[1]//len(idxs)}%)", "allele": f"{ald[0]}({100*ald[1]//len(idxs)}%)"}
    dn = dendrogram(Z, orientation="left", no_plot=True); order = dn["leaves"][::-1]
    meth = np.array([P[i] for i in order]); dist = sub[np.ix_(order, order)].copy()
    np.fill_diagonal(dist, 0); dist[dist < 0] = np.nan; lab_o = [lab[i] for i in order]
    hp_o = [reads[ids[i]]["hp"] for i in order]; al_o = [reads[ids[i]]["alt_support"] for i in order]
    st_o = [reads[ids[i]].get("strand", "?") for i in order]
    nodelab = {nd: (list({lab[i] for i in desc[nd]})[0] if len({lab[i] for i in desc[nd]}) == 1 else None) for nd in range(2 * n - 1)}
    sb = [("演化群", [labcol(l) for l in lab_o]), ("HP", [HPC.get(h, "#ddd") for h in hp_o]),
          ("ALT", [alc.get(a, "#ddd") for a in al_o]), ("strand", ["#444" if s == "+" else "#bbb" for s in st_o])]
    mc, dc = H.mpl_cmaps(); nsb = len(sb)
    wr = [1.3] + [0.05] * nsb + [1.0, 0.16] + [0.05] * nsb + [1.0]
    fig = plt.figure(figsize=(11, 5.2)); gs = fig.add_gridspec(2, len(wr), width_ratios=wr, height_ratios=[1, 0.42], wspace=0.04, hspace=0.02)
    c = 0; axdn = fig.add_subplot(gs[0, c]); c += 1
    dendrogram(Z, orientation="left", link_color_func=lambda k: labcol(nodelab.get(k)) if nodelab.get(k) else "#ccc", ax=axdn, no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA 樹(分支=v2演化群)", fontsize=8); [axdn.spines[s2].set_visible(False) for s2 in axdn.spines]
    for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
    axm = fig.add_subplot(gs[0, c]); c += 1; axm.imshow(meth, aspect="auto", cmap=mc, vmin=0, vmax=1, interpolation="nearest")
    sx = H.snv_fractional_x(cpgs, int(key.split('_')[1]))
    if sx is not None: axm.axvline(sx, color=H.SNV_COL, lw=2)
    axm.set_xticks([]); axm.set_yticks([]); axm.set_title("甲基 read×CpG", fontsize=8.5)
    fig.add_subplot(gs[0, c]).axis("off"); c += 1
    for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
    axd = fig.add_subplot(gs[0, c]); vmax = max(0.5, float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist, aspect="auto", cmap=dc, vmin=0, vmax=vmax, interpolation="nearest")
    axd.set_xticks([]); axd.set_yticks([]); axd.set_title("距離(對角塊=群)", fontsize=8.5)
    axl = fig.add_subplot(gs[1, :]); axl.axis("off")
    txt = " · ".join(f"{L}: n{galign[L]['n']} hp={galign[L]['hp']} alle={galign[L]['allele']}" for L in sorted(labct))
    axl.text(0.0, 0.7, f"v2 切割: 穩定群 {nterm} 個 + 離群 {nout}（標籤 {sorted(labct)}）", fontsize=8.5, fontweight="bold", transform=axl.transAxes)
    axl.text(0.0, 0.25, txt, fontsize=7.5, transform=axl.transAxes, wrap=True)
    fig.suptitle(f"{key}  n={n} CpG={len(cpgs)}  fine_class={it['class']}  →  v2 演化群 {nterm} 個 離群 {nout}", fontsize=9.5, y=1.0)
    fn = f"{FIGS}/v2_{key}.png"; fig.savefig(fn, dpi=104, bbox_inches="tight"); plt.close(fig)
    out.append({"key": key, "n": n, "C": len(cpgs), "fine_class": it["class"], "v2_ngroups": nterm,
                "labels": dict(labct), "align": galign, "n_outlier": nout, "png": f"figs_phylo_v2/v2_{key}.png"})
    print(f"  {key} n={n} C={len(cpgs)} → v2 群{nterm} {sorted(labct)} 離群{nout}", flush=True)
json.dump(out, open(f"{A}/phylo_v2_render.json", "w"), indent=1)
print(f"\nDONE {len(out)} 位點 · v2 多群 {sum(1 for o in out if o['v2_ngroups']>=2)} · 圖在 figs_phylo_v2/")
