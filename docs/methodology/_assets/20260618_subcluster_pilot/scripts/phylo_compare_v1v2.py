#!/usr/bin/env python3
"""[驗證] phylo 修法前後對照 — v1(現法,double-dip) vs v2(per-subgroup 重分群 null + RNULL=40)。
v2 修法: 真實分割不變(UPGMA 節點子節點=該子樹最佳2way分割)，只把 null 改成『對該子群 read 欄內置換→重分群→取其 root 分裂比』(消 double-dip + depth≥1 用對 null)。
輸出哪些位點/標籤在 v2 下消失或變淺。重用 _ws_render 快取矩陣(免 binary)。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H  # noqa (for LABMAP parity)
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
MINSZ = 3; SEP_MIN = 1.3


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
    """mode='v1': null 沿用真實 S1/S2(double-dip). mode='v2': per-subgroup 重分群 null."""
    Z, s, desc, ch, n = _tree(sub)
    # v1 全域 once null
    global_nulls = None
    if mode == "v1":
        C = P.shape[1]; global_nulls = []
        for _ in range(RNULL):
            Pn = P.copy()
            for cc in range(C):
                col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
            Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); global_nulls.append(np.maximum(Dn, Dn.T))

    def split_real(node):
        if node not in ch: return False
        c1, c2 = ch[node]
        S1, S2 = np.array(desc[c1]), np.array(desc[c2])
        if min(len(S1), len(S2)) < MINSZ: return False
        r = _bw(s, S1, S2)
        if r is None or r < SEP_MIN: return False
        if mode == "v1":
            ns = [_bw(Dn, S1, S2) for Dn in global_nulls]
        else:
            S = np.array(desc[node]); m = len(S); Psub = P[S]; ns = []
            for _ in range(RNULL):
                Pn = Psub.copy()
                for cc in range(Pn.shape[1]):
                    col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                    if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
                Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
                _, sn, dn, cn, _ = _tree(Dn)
                nc1, nc2 = cn[2 * m - 2]
                ns.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
        ns = [x for x in ns if x is not None]
        return r > (np.percentile(ns, 95) if ns else 0)

    lab = [None] * n

    def rec(node, label):
        leaves = desc[node]
        if node not in ch or len(leaves) < 2 * MINSZ:
            for i in leaves: lab[i] = label; return
        if split_real(node):
            c1, c2 = ch[node]
            big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            rec(big, label + "-1"); rec(small, label + "-2")
        else:
            for i in leaves: lab[i] = label

    root = 2 * n - 2
    if split_real(root):
        c1, c2 = ch[root]
        big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1); rec(big, "1"); rec(small, "2")
    else:
        for i in range(n): lab[i] = "1"
    lab = [l if l else "outlier" for l in lab]
    sm = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    lab = [("outlier" if l in sm else l) for l in lab]
    return lab


# ---- load loci (same as phylo_groups.py) ----
dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd
items = json.load(open(f"{A}/ws_items.json")); rng = np.random.default_rng(20260622); out = []
for it in items:
    key = it["key"]; rd = dirmap.get(key)
    if not rd: continue
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs = rows[0].split(",")[1:]
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids) < MINSZ * 2: continue
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    nC = len(cpgs)
    lab1 = phylo_label(sub, P, np.random.default_rng(20260622), "v1", 12)
    lab2 = phylo_label(sub, P, np.random.default_rng(20260622), "v2", 40)
    g1 = dict(Counter(l for l in lab1 if l != "outlier")); g2 = dict(Counter(l for l in lab2 if l != "outlier"))
    out.append({"key": key, "n": n, "C": nC, "fine": it["fine_conf"],
                "v1_ngroups": len(g1), "v1_labels": sorted(g1), "v1_maxdepth": max((L.count('-') for L in g1), default=0),
                "v2_ngroups": len(g2), "v2_labels": sorted(g2), "v2_maxdepth": max((L.count('-') for L in g2), default=0)})
    chg = "" if (len(g1) == len(g2) and out[-1]["v1_maxdepth"] == out[-1]["v2_maxdepth"]) else "  ← 變了"
    print(f"  {key} n={n} C={nC} | v1: {len(g1)}群{sorted(g1)} | v2: {len(g2)}群{sorted(g2)}{chg}", flush=True)
json.dump(out, open(f"{A}/phylo_v1v2_compare.json", "w"), indent=1)
# summary
chg_n = sum(1 for o in out if o["v1_ngroups"] != o["v2_ngroups"])
chg_depth = sum(1 for o in out if o["v1_maxdepth"] != o["v2_maxdepth"])
m1 = [o for o in out if o["v1_ngroups"] >= 2]; m2 = [o for o in out if o["v2_ngroups"] >= 2]
print(f"\nDONE {len(out)} 位點 | v1 多群 {len(m1)} → v2 多群 {len(m2)} | 群數變 {chg_n} | 最深標籤變淺 {chg_depth}")
print(f"v1 平均群 {np.mean([o['v1_ngroups'] for o in out]):.2f} → v2 平均群 {np.mean([o['v2_ngroups'] for o in out]):.2f}")
