#!/usr/bin/env python3
"""[修法+驗證] phylo v3 = v2 + quarantine-descend 修 FM1 復發(單離群吃群位)。
v3 核心: rec 進 node 時, 先剝離 (tiny,big) caterpillar 單離群(quarantine), descend 到『兩子群皆≥MINSZ』的平衡節點再測 split(per-subgroup null)。
驗證: (A) 4 個 flagged 位點 v3 是否救回切割+對齊; (B) 全 30 位點 v2→v3 群數變化(救回幾個/6個原多群是否保留);
(C) 純噪音 FP(C=76) 確認 v3 不過切。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H  # noqa
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
from scipy.cluster.hierarchy import fcluster
WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
MINSZ = 3; SEP_MIN = 1.3; RNULL = 40
FLAGGED = ["chr22_26939195", "chr22_30454004", "chr20_42981498", "chr21_40426852"]


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


def phylo_label(sub, P, rng, mode):
    """mode='v2': 只測 root sibling(原bug). mode='v3': quarantine-descend 到平衡節點."""
    Z, s, desc, ch, n = _tree(sub)

    def subnull95(leaves):
        S = np.array(leaves); m = len(S); Psub = P[S]; ns = []
        for _ in range(RNULL):
            Pn = Psub.copy()
            for cc in range(Pn.shape[1]):
                col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
            Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
            _, sn, dn, cn, _ = _tree(Dn); nc1, nc2 = cn[2 * m - 2]
            ns.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
        ns = [x for x in ns if x is not None]
        return float(np.percentile(ns, 95)) if ns else 0

    def descend(node):
        """v3: 剝離 (tiny,big) caterpillar 至兩子群皆≥MINSZ 的平衡節點。回傳(bnode, quarantined)或(None,...)。"""
        cur = node; quar = []
        while cur in ch:
            c1, c2 = ch[cur]; s1, s2 = len(desc[c1]), len(desc[c2])
            if min(s1, s2) >= MINSZ: return cur, quar
            small, big = (c1, c2) if s1 < s2 else (c2, c1)
            quar += desc[small]; cur = big
        return None, quar

    lab = [None] * n

    def test_split(bnode):
        c1, c2 = ch[bnode]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
        r = _bw(s, S1, S2)
        if r is None or r < SEP_MIN: return False, r
        return r > subnull95(desc[bnode]), r

    def rec(node, label):
        leaves = desc[node]
        if len(leaves) < 2 * MINSZ:
            for i in leaves: lab[i] = label; return
        if mode == "v3":
            # v3.1: quarantine 物有所值 — 只在分裂成立時才剝離離群; 否則整群保留不侵蝕
            bnode, quar = descend(node)
            ok = False
            if bnode is not None:
                ok, r = test_split(bnode)
            if ok:
                for i in quar: lab[i] = "outlier"  # commit quarantine only on success
                c1, c2 = ch[bnode]; big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
                rec(big, label + "-1"); rec(small, label + "-2")
            else:
                for i in leaves: lab[i] = label  # 分裂不成立 → 整 node 保留, 不 quarantine
        else:  # v2 原法
            if node not in ch:
                for i in leaves: lab[i] = label; return
            c1, c2 = ch[node]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
            if min(len(S1), len(S2)) < MINSZ:
                for i in leaves: lab[i] = label; return
            ok, r = test_split(node)
            if ok:
                big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
                rec(big, label + "-1"); rec(small, label + "-2")
            else:
                for i in leaves: lab[i] = label

    rec(2 * n - 2, "1")
    lab = [l if l else "outlier" for l in lab]
    sm = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    lab = [("outlier" if l in sm else l) for l in lab]
    # 若只有單一非outlier標籤但全是 "1" → 正規化
    return lab


dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd


def load(key):
    rd = dirmap.get(key)
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n")
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]
    P = np.array([M[mi[x]] for x in ids]); C = M.shape[1]
    return ids, sub, P, reads, C


def align(lab, ids, reads, n):
    labct = Counter(l for l in lab if l != "outlier"); g = {}
    for L in labct:
        idx = [i for i in range(n) if lab[i] == L]
        hpc = Counter(reads[ids[i]]["hp"] for i in idx).most_common(1)[0]
        alc = Counter(reads[ids[i]]["alt_support"] for i in idx).most_common(1)[0]
        g[L] = f"n{len(idx)} hp={hpc[0]}({100*hpc[1]//len(idx)}%) alle={alc[0]}({100*alc[1]//len(idx)}%)"
    return labct, g


ws = {x["key"]: x for x in json.load(open(f"{A}/ws_items.json"))}
v2final = {x["key"]: x["v2_ngroups"] for x in json.load(open(f"{A}/phylo_v2_final.json"))}
print("=== (A) 4 個 flagged 位點: v3 是否救回 ===")
for key in FLAGGED:
    ids, sub, P, reads, C = load(key); n = len(ids)
    l2 = phylo_label(sub, P, np.random.default_rng(20260622), "v2")
    l3 = phylo_label(sub, P, np.random.default_rng(20260622), "v3")
    c2, _ = align(l2, ids, reads, n); c3, g3 = align(l3, ids, reads, n)
    print(f"\n{key} n={n} 4gate={ws[key]['class']}(k{ws[key]['fine_k']}) | v2={len(c2)}群 → v3={len(c3)}群")
    for L in sorted(g3): print(f"    群{L}: {g3[L]}")

print("\n=== (B) 全 30 位點 v2→v3 群數 ===")
rec = []
for it in ws.values():
    key = it["key"]
    if key not in dirmap: continue
    ids, sub, P, reads, C = load(key); n = len(ids)
    if n < MINSZ * 2: continue
    l3 = phylo_label(sub, P, np.random.default_rng(20260622), "v3")
    c3, g3 = align(l3, ids, reads, n)
    v2g = v2final.get(key, "?")
    rec.append({"key": key, "n": n, "v2": v2g, "v3": len(c3), "4gate": it["fine_k"] if it["class"] != "NO_CLEAR_SPLIT" else 1,
                "class": it["class"], "align": g3})
    ch = " ←救回" if (isinstance(v2g, int) and len(c3) > v2g) else (" ←減少" if isinstance(v2g, int) and len(c3) < v2g else "")
    print(f"  {key} n={n} | v2={v2g} v3={len(c3)} 4gate={rec[-1]['4gate']}{ch}")
json.dump(rec, open(f"{A}/phylo_v3_validate.json", "w"), ensure_ascii=False, indent=1)
rescued = sum(1 for r in rec if isinstance(r["v2"], int) and r["v3"] > r["v2"])
reduced = sum(1 for r in rec if isinstance(r["v2"], int) and r["v3"] < r["v2"])
v3multi = sum(1 for r in rec if r["v3"] >= 2)
print(f"\nv2 多群 {sum(1 for r in rec if isinstance(r['v2'],int) and r['v2']>=2)} → v3 多群 {v3multi} | 救回 {rescued} | 減少 {reduced} | v3 平均 {np.mean([r['v3'] for r in rec]):.2f}")

print("\n=== (C) 純噪音 FP（v3, C=76, TRIALS=120）確認不過切 ===")
def pn(n, C, rng, miss=0.25):
    rt = rng.uniform(0.1, 0.9, C); Q = (rng.random((n, C)) < rt[None, :]).astype(float)
    Q[rng.random((n, C)) < miss] = np.nan
    for cc in range(C):
        if np.sum(~np.isnan(Q[:, cc])) < 2:
            ix = rng.choice(n, 2, replace=False); Q[ix, cc] = (rng.random(2) < rt[cc]).astype(float)
    return Q
rng = np.random.default_rng(99)
for n in [20, 40, 80]:
    fp = 0; TR = 120
    for t in range(TR):
        Q = pn(n, 76, rng); Dn = CR.bernoulli_dist(Q); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
        l3 = phylo_label(Dn, Q, rng, "v3")
        if len(set(x for x in l3 if x != "outlier")) >= 2: fp += 1
    print(f"  n={n}: v3 噪音 FP = {100*fp/TR:.1f}%")
