#!/usr/bin/env python3
"""[數據解釋門檻] phylo 各門檻敏感度掃描 — 用數據展示『換個門檻會怎樣』,區分 data-derived vs convention。
掃 SEP_MIN / null 百分位 / RNULL / MINSZ → 量 純噪音FP(specificity) + 30位點多群數(sensitivity) + 4旗標群數。
回答: 每個門檻是否能用數據解釋其選值。重用快取矩陣。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
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


def label_v3(sub, P, rng, SEP_MIN=1.3, NULLPCT=95, RNULL=40, MINSZ=3):
    """參數化 v3.1 labeler。回傳 (群數, 離群數)。"""
    Z, s, desc, ch, n = _tree(sub)

    def subnull(leaves):
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
        return float(np.percentile(ns, NULLPCT)) if ns else 0

    def descend(node):
        cur = node; quar = []
        while cur in ch:
            c1, c2 = ch[cur]; s1, s2 = len(desc[c1]), len(desc[c2])
            if min(s1, s2) >= MINSZ: return cur, quar
            small, big = (c1, c2) if s1 < s2 else (c2, c1); quar += desc[small]; cur = big
        return None, quar

    def test(bnode):
        c1, c2 = ch[bnode]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
        r = _bw(s, S1, S2)
        if r is None or r < SEP_MIN: return False
        return r > subnull(desc[bnode])

    lab = [None] * n

    def rec(node, label):
        leaves = desc[node]
        if len(leaves) < 2 * MINSZ:
            for i in leaves: lab[i] = label; return
        bnode, quar = descend(node)
        if bnode is not None and test(bnode):
            for i in quar: lab[i] = "outlier"
            c1, c2 = ch[bnode]; big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            rec(big, label + "-1"); rec(small, label + "-2")
        else:
            for i in leaves: lab[i] = label
    rec(2 * n - 2, "1")
    lab = [l if l else "outlier" for l in lab]
    sm = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    lab = [("outlier" if l in sm else l) for l in lab]
    return len(set(l for l in lab if l != "outlier")), sum(1 for l in lab if l == "outlier")


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
    return sub, np.array([M[mi[x]] for x in ids])


def pn(n, C, rng, miss=0.25):
    rt = rng.uniform(0.1, 0.9, C); Q = (rng.random((n, C)) < rt[None, :]).astype(float)
    Q[rng.random((n, C)) < miss] = np.nan
    for cc in range(C):
        if np.sum(~np.isnan(Q[:, cc])) < 2:
            ix = rng.choice(n, 2, replace=False); Q[ix, cc] = (rng.random(2) < rt[cc]).astype(float)
    return Q


def noise_fp(SEP_MIN, NULLPCT, RNULL, MINSZ, n=40, C=76, TR=80):
    rng = np.random.default_rng(7); fp = 0
    for _ in range(TR):
        Q = pn(n, C, rng); D = CR.bernoulli_dist(Q); np.fill_diagonal(D, 0); D = np.maximum(D, D.T)
        g, _ = label_v3(D, Q, rng, SEP_MIN, NULLPCT, RNULL, MINSZ)
        if g >= 2: fp += 1
    return 100 * fp / TR


def multi30(SEP_MIN, NULLPCT, RNULL, MINSZ):
    ws = json.load(open(f"{A}/ws_items.json")); cnt = 0; flag = {}
    for it in ws:
        key = it["key"]
        if key not in dirmap: continue
        sub, P = load(key)
        if sub.shape[0] < 6: continue
        g, _ = label_v3(sub, P, np.random.default_rng(20260622), SEP_MIN, NULLPCT, RNULL, MINSZ)
        if g >= 2: cnt += 1
        if key in FLAGGED: flag[key] = g
    return cnt, flag


out = {}
print("=== 掃 1: SEP_MIN（固定 null95/RNULL40/MINSZ3）— 核心可調門檻 ===")
print(f"{'SEP_MIN':>8} | {'噪音FP%':>7} | {'30位點多群':>9} | 4旗標群數 [26939195,30454004,42981498,40426852]")
sep_rows = []
for sm in [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]:
    fp = noise_fp(sm, 95, 40, 3); m, fl = multi30(sm, 95, 40, 3)
    sep_rows.append({"SEP_MIN": sm, "noise_fp": fp, "multi30": m, "flagged": fl})
    print(f"{sm:>8} | {fp:>6.1f}% | {m:>9} | {[fl.get(k) for k in FLAGGED]}")
out["sep_min"] = sep_rows

print("\n=== 掃 2: null 百分位（固定 SEP_MIN1.3/RNULL40/MINSZ3）===")
print(f"{'null pct':>8} | {'噪音FP%':>7} | {'30位點多群':>9}")
pct_rows = []
for pc in [90, 95, 99]:
    fp = noise_fp(1.3, pc, 40, 3); m, fl = multi30(1.3, pc, 40, 3)
    pct_rows.append({"null_pct": pc, "noise_fp": fp, "multi30": m}); print(f"{pc:>8} | {fp:>6.1f}% | {m:>9}")
out["null_pct"] = pct_rows

print("\n=== 掃 3: RNULL（null95 估計穩定性；同位點重算 5 次看 verdict 翻不翻）===")
rn_rows = []
for key in ["chr22_26939195", "chr20_42981498"]:
    sub, P = load(key); print(f"  {key}:")
    for RN in [12, 25, 40, 80]:
        gs = [label_v3(sub, P, np.random.default_rng(1000 + i), 1.3, 95, RN, 3)[0] for i in range(5)]
        rn_rows.append({"key": key, "RNULL": RN, "groups_5runs": gs})
        print(f"    RNULL={RN:>3}: 5 次群數 = {gs} {'(穩定)' if len(set(gs))==1 else '(翻動!)'}")
out["rnull"] = rn_rows

print("\n=== 掃 4: MINSZ（固定 SEP_MIN1.3/null95/RNULL40）===")
print(f"{'MINSZ':>6} | {'噪音FP%':>7} | {'30位點多群':>9}")
ms_rows = []
for ms in [3, 4, 5]:
    fp = noise_fp(1.3, 95, 40, ms); m, fl = multi30(1.3, 95, 40, ms)
    ms_rows.append({"MINSZ": ms, "noise_fp": fp, "multi30": m}); print(f"{ms:>6} | {fp:>6.1f}% | {m:>9}")
out["minsz"] = ms_rows
json.dump(out, open(f"{A}/phylo_sensitivity.json", "w"), indent=1)
print("\nDONE → phylo_sensitivity.json")
