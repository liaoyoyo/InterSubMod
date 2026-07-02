#!/usr/bin/env python3
"""[診斷] phylo 遞迴逐步 trace — 精確定位為什麼漏掉眼睛看到的平衡切割。
對 4 位點: 印出 phylo 遞迴每一步 split_real(node) 測試(node sizes / r / sep_pass / null_pass / verdict);
並對照『平衡切割』(fcluster k=2/3 去單離群後) 的 r,看它是否 ≥SEP_MIN、phylo 遞迴是否曾測到它。
結論: bug 是 (A)caterpillar 遞迴沒到平衡節點 還是 (B)平衡節點 r<SEP_MIN 被拒。"""
import os, csv, glob, sys
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
from scipy.cluster.hierarchy import fcluster
WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
MINSZ = 3; SEP_MIN = 1.3; RNULL = 40
TARGETS = ["chr22_26939195", "chr22_30454004", "chr20_42981498", "chr21_40426852"]


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


dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd

for key in TARGETS:
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
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    Z, s, desc, ch, _ = _tree(sub)
    # 全域 null for v2-style per-subgroup (簡化: 用全域近似看 r 是否>null; 真 v2 每子群算,這裡標 r 與 SEP_MIN 即可定性)
    print(f"\n{'='*70}\n{key}  n={n}")
    # trace phylo 遞迴 (v2 per-subgroup null)
    steps = []

    def subnull95(leaves):
        S = np.array(leaves); m = len(S); Psub = P[S]; ns = []
        for _ in range(RNULL):
            Pn = Psub.copy()
            for cc in range(Pn.shape[1]):
                col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                if vi.size > 1: Pn[vi, cc] = col[np.random.default_rng(20260622 + _).permutation(vi)]
            Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
            _, sn, dn, cn, _ = _tree(Dn); nc1, nc2 = cn[2 * m - 2]
            ns.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
        ns = [x for x in ns if x is not None]
        return float(np.percentile(ns, 95)) if ns else 0

    def rec(node, label, depth):
        leaves = desc[node]
        if node not in ch or len(leaves) < 2 * MINSZ:
            steps.append((depth, label, len(leaves), "終端(太小/葉)", None, None, None)); return
        c1, c2 = ch[node]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
        if min(len(S1), len(S2)) < MINSZ:
            r = _bw(s, S1, S2)
            steps.append((depth, label, len(leaves), f"子群太小 sizes({len(S1)},{len(S2)})", round(r, 3) if r else None, "—", "拒(子群<MINSZ)"))
            # 仍會 rec 進去? 不: split_real return False → 整 node 變終端
            return
        r = _bw(s, S1, S2)
        if r is None or r < SEP_MIN:
            steps.append((depth, label, len(leaves), f"sizes({len(S1)},{len(S2)})", round(r, 3) if r else None, f"<SEP_MIN {SEP_MIN}", "拒(r<1.3)")); return
        nl = subnull95(leaves)
        if r > nl:
            big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            steps.append((depth, label, len(leaves), f"sizes({len(S1)},{len(S2)})", round(r, 3), f">null95 {round(nl,3)}", "切→遞迴"))
            rec(big, label + "-1", depth + 1); rec(small, label + "-2", depth + 1)
        else:
            steps.append((depth, label, len(leaves), f"sizes({len(S1)},{len(S2)})", round(r, 3), f"<null95 {round(nl,3)}", "拒(r<null95)"))

    root = 2 * n - 2; c1, c2 = ch[root]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
    r0 = _bw(s, S1, S2)
    if r0 is not None and r0 >= SEP_MIN and min(len(S1), len(S2)) >= MINSZ:
        nl = subnull95(desc[root])
        if r0 > nl:
            big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            steps.append((0, "root", n, f"sizes({len(S1)},{len(S2)})", round(r0, 3), f">null95 {round(nl,3)}", "切→遞迴"))
            rec(big, "1", 1); rec(small, "2", 1)
        else:
            steps.append((0, "root", n, f"sizes({len(S1)},{len(S2)})", round(r0, 3), f"<null95 {round(nl,3)}", "拒"))
    else:
        steps.append((0, "root", n, f"sizes({len(S1)},{len(S2)})", round(r0, 3) if r0 else None, "—", f"拒(r<1.3 或子群<MINSZ)"))
    print("  phylo 遞迴 trace:")
    for d, lab, sz, szinfo, r, nullinfo, verdict in steps:
        print(f"    {'  '*d}[{lab}] n={sz} {szinfo} r={r} {nullinfo or ''} → {verdict}")
    # 對照: 平衡 fcluster k=2/3 去單離群後的 r
    print("  平衡切割(fcluster, 去 size<MINSZ 群)對照:")
    for k in [2, 3]:
        if n < k * MINSZ: continue
        lab = fcluster(Z, k, "maxclust")
        groups = {g: np.where(lab == g)[0] for g in set(lab)}
        big_groups = {g: ix for g, ix in groups.items() if len(ix) >= MINSZ}
        if len(big_groups) < 2:
            print(f"    k={k}: 去離群後 <2 群"); continue
        gl = list(big_groups.values())
        # 取最大兩群算 r
        gl.sort(key=len, reverse=True)
        r_bal = _bw(s, gl[0], gl[1])
        # 該 2 群對齊
        hp = [reads[ids[i]]["hp"] for i in range(n)]; al = [reads[ids[i]]["alt_support"] for i in range(n)]
        sizes = [len(g) for g in gl]
        print(f"    k={k}: 群sizes={sizes} 最大兩群 r={round(r_bal,3) if r_bal else None} (SEP_MIN={SEP_MIN} → {'≥1.3 phylo該切但沒測到' if r_bal and r_bal>=SEP_MIN else '<1.3 phylo判不夠分離'})")
