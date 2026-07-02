#!/usr/bin/env python3
"""[驗證] phylo split_real 純噪音偽陽性率校準 — 獨立驗證 workflow 的 double-dipping 指控。
純噪音 = 每 read 對每 CpG 獨立抽 Bernoulli(per-CpG rate)，有邊際結構但零 read×read 結構。
若 phylo_label 在此資料上回傳 n_groups≥2 = 偽陽性。複製 phylo_groups.py 的 split_real/phylo_label。
對照組: split_real 改成『每個 null 重新分群取最佳 sibling-split ratio』(stab_excess 式) 看 FP 是否降。"""
import sys, numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import cluster_redesign as CR
from scipy.cluster.hierarchy import fcluster

MINSZ = 3; SEP_MIN = 1.3; RNULL = 12

def _bw(s, S1, S2):
    bet = s[np.ix_(S1, S2)]; bet = bet[bet >= 0]
    w1 = s[np.ix_(S1, S1)][np.triu_indices(len(S1), 1)]; w1 = w1[w1 >= 0]
    w2 = s[np.ix_(S2, S2)][np.triu_indices(len(S2), 1)]; w2 = w2[w2 >= 0]
    wm = np.concatenate([w1, w2]) if (w1.size or w2.size) else np.array([])
    if bet.size == 0 or wm.size == 0 or wm.mean() <= 1e-6: return None
    return float(bet.mean()) / float(wm.mean())

def make_nulls(P, rng, rnull=RNULL):
    C = P.shape[1]; nulls = []
    for _ in range(rnull):
        Pn = P.copy()
        for cc in range(C):
            col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
            if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
        Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); nulls.append(np.maximum(Dn, Dn.T))
    return nulls

def phylo_nsplit(sub, P, rng, mode="original"):
    """回傳是否在 root 判定 split_real（=至少 2 群）。mode=original(複製現法) / recluster(每 null 重分群)。"""
    n = sub.shape[0]; Z, s = CR.linkZ(sub)
    desc = {i: [i] for i in range(n)}; ch = {}
    for i in range(len(Z)):
        a, b = int(Z[i, 0]), int(Z[i, 1]); desc[n + i] = desc[a] + desc[b]; ch[n + i] = (a, b)
    nulls = make_nulls(P, rng)
    root = 2 * n - 2; c1, c2 = ch[root]
    S1, S2 = np.array(desc[c1]), np.array(desc[c2])
    if min(len(S1), len(S2)) < MINSZ: return False
    r = _bw(s, S1, S2)
    if r is None or r < SEP_MIN: return False
    if mode == "original":
        # 現法: null 沿用真實樹的 S1/S2 分割（double-dip 嫌疑）
        ns = [_bw(Dn, S1, S2) for Dn in nulls]
    else:
        # 對照: 每個 null 自己重新 UPGMA 分 2 群，取其 root sibling-split ratio（消 selection bias）
        ns = []
        for Dn in nulls:
            Zn, sn = CR.linkZ(Dn)
            dn = {i: [i] for i in range(n)}; cn = {}
            for i in range(len(Zn)):
                a, b = int(Zn[i, 0]), int(Zn[i, 1]); dn[n + i] = dn[a] + dn[b]; cn[n + i] = (a, b)
            rc1, rc2 = cn[2 * n - 2]
            ns.append(_bw(sn, np.array(dn[rc1]), np.array(dn[rc2])))
    ns = [x for x in ns if x is not None]
    return r > (np.percentile(ns, 95) if ns else 0)

def pure_noise_P(n, C, rng, miss=0.25):
    rates = rng.uniform(0.1, 0.9, C)
    P = (rng.random((n, C)) < rates[None, :]).astype(float)
    P[rng.random((n, C)) < miss] = np.nan
    # 保證每欄至少 2 個非 NaN（否則距離退化）
    for cc in range(C):
        if np.sum(~np.isnan(P[:, cc])) < 2:
            ix = rng.choice(n, 2, replace=False); P[ix, cc] = (rng.random(2) < rates[cc]).astype(float)
    return P

def phylo_label_full(sub, P, rng):
    """完整複製 phylo_groups.py phylo_label（遞迴 + 全域 once null）→ 回傳終端群數(≥2=多群)。"""
    n = sub.shape[0]; Z, s = CR.linkZ(sub)
    desc = {i: [i] for i in range(n)}; ch = {}
    for i in range(len(Z)):
        a, b = int(Z[i, 0]), int(Z[i, 1]); desc[n + i] = desc[a] + desc[b]; ch[n + i] = (a, b)
    nulls = make_nulls(P, rng)  # ← 全域只算一次（workflow 指控 (b)）
    def split_real(c1, c2):
        S1, S2 = np.array(desc[c1]), np.array(desc[c2])
        if min(len(S1), len(S2)) < MINSZ: return False
        r = _bw(s, S1, S2)
        if r is None or r < SEP_MIN: return False
        ns = [_bw(Dn, S1, S2) for Dn in nulls]; ns = [x for x in ns if x is not None]
        return r > (np.percentile(ns, 95) if ns else 0)
    lab = [None] * n
    def rec(node, label):
        leaves = desc[node]
        if node not in ch or len(leaves) < 2 * MINSZ:
            for i in leaves: lab[i] = label; return
        c1, c2 = ch[node]
        if split_real(c1, c2):
            big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            rec(big, label + "-1"); rec(small, label + "-2")
        else:
            for i in leaves: lab[i] = label
    root = 2 * n - 2; c1, c2 = ch[root]
    if split_real(c1, c2):
        big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1); rec(big, "1"); rec(small, "2")
    else:
        for i in range(n): lab[i] = "1"
    from collections import Counter
    lab = [l if l else "outlier" for l in lab]
    sm = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    lab = [("outlier" if l in sm else l) for l in lab]
    return len(set(l for l in lab if l != "outlier"))

if __name__ == "__main__":
    rng = np.random.default_rng(20260622)
    TRIALS = 200
    print("=== 測試 1: ROOT split FP（C=40 CpG, 25% 缺值）— 隔離指控(a) double-dip ===")
    print(f"{'n':>4} | {'現法 FP%':>9} | {'重分群對照 FP%':>13}")
    for n in [20, 40, 80]:
        fp_o = fp_r = 0
        for t in range(TRIALS):
            P = pure_noise_P(n, 40, rng)
            D = CR.bernoulli_dist(P); np.fill_diagonal(D, 0); D = np.maximum(D, D.T)
            if phylo_nsplit(D, P, rng, "original"): fp_o += 1
            if phylo_nsplit(D, P, rng, "recluster"): fp_r += 1
        print(f"{n:>4} | {100*fp_o/TRIALS:>8.1f}% | {100*fp_r/TRIALS:>12.1f}%")
    print("\n=== 測試 2: 完整遞迴 phylo_label FP（n_groups≥2 即 FP）— 含指控(b) 遞迴 null ===")
    print(f"{'n':>4} | {'C=10':>7} | {'C=20':>7} | {'C=40':>7} | {'C=76(中位)':>10}")
    for n in [20, 40, 80]:
        row = []
        for C in [10, 20, 40, 76]:
            fp = 0
            for t in range(TRIALS):
                P = pure_noise_P(n, C, rng)
                D = CR.bernoulli_dist(P); np.fill_diagonal(D, 0); D = np.maximum(D, D.T)
                if phylo_label_full(D, P, rng) >= 2: fp += 1
            row.append(100 * fp / TRIALS)
        print(f"{n:>4} | {row[0]:>6.1f}% | {row[1]:>6.1f}% | {row[2]:>6.1f}% | {row[3]:>9.1f}%")
