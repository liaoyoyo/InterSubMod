#!/usr/bin/env python3
"""[驗證] 4 切群方法對照 — silhouette best_k(C++現用) vs 4-gate(cluster_redesign) vs phylo-v1 vs phylo-v2。
(A) 同 30 pilot 位點各方法群數 + 一致性。
(B) 純噪音偽陽性 FP(specificity) = 誰更好的關鍵證據。
silhouette best_k 強制 k∈[2,maxk](無『1群』選項); +PERMANOVA = C++ 實際 pipeline(gate 後)。
4-gate/phylo 有 null gate 可判『無結構』。重用快取矩陣 + ws_items(4-gate) + phylo json。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H  # noqa
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
from scipy.cluster.hierarchy import fcluster
try:
    from sklearn.metrics import silhouette_score
    HAVE_SK = True
except Exception:
    HAVE_SK = False

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
MINSZ = 3


def sil_score(D, lab):
    if HAVE_SK:
        try: return silhouette_score(D, lab, metric="precomputed")
        except Exception: return -2
    # manual silhouette (precomputed)
    n = len(lab); groups = {}
    for i, g in enumerate(lab): groups.setdefault(g, []).append(i)
    if len(groups) < 2: return -2
    s = []
    for i in range(n):
        same = [j for j in groups[lab[i]] if j != i]
        a = np.mean([D[i, j] for j in same]) if same else 0
        b = min(np.mean([D[i, j] for j in groups[g]]) for g in groups if g != lab[i])
        s.append((b - a) / max(a, b) if max(a, b) > 0 else 0)
    return float(np.mean(s))


def silhouette_bestk(D, maxk):
    """C++ TreeCutter::find_optimal_clusters 等價: 對 k∈[2,maxk] 取最大 silhouette。無 1 群選項。"""
    Z, _ = CR.linkZ(D); best = (2, -2)
    for k in range(2, maxk + 1):
        lab = fcluster(Z, k, criterion="maxclust")
        if len(set(lab)) < 2: continue
        sc = sil_score(D, list(lab))
        if sc > best[1]: best = (k, sc)
    return best[0], best[1], fcluster(Z, best[0], criterion="maxclust")


def permanova_p(D, lab, rng, nperm=199):
    """StructureTest.cpp 等價 pseudo-F + 標籤置換。"""
    lab = np.array(lab); n = len(lab); k = len(set(lab))
    if k < 2 or n <= k: return 1.0

    def ss(L):
        groups = {}
        for i, g in enumerate(L): groups.setdefault(g, []).append(i)
        total = within = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                d2 = D[i, j] ** 2; total += d2
                if L[i] == L[j]: pass
        # within per group / n_g
        for g, m in groups.items():
            for a in range(len(m)):
                for b in range(a + 1, len(m)):
                    within += D[m[a], m[b]] ** 2 / len(m)
        sst = total / n; ssw = within; ssb = max(sst - ssw, 0)
        return ssb, ssw
    ssb, ssw = ss(lab)
    f = (ssb / (k - 1)) / (ssw / (n - k)) if ssw > 0 else (1e9 if ssb > 0 else 0)
    ne = 1
    for _ in range(nperm):
        pl = lab.copy(); rng.shuffle(pl)
        b, w = ss(pl); pf = (b / (k - 1)) / (w / (n - k)) if w > 0 else (1e9 if b > 0 else 0)
        if pf >= f: ne += 1
    return ne / (nperm + 1)


def fourgate_real(D, P, rng):
    """4-gate 的 real 閘(label-free): clusterboot excess ≥ 0.10 at k=2 = 有結構。"""
    try:
        ex = CR.stab_excess(D, P, 2, rng)
        return (ex is not None) and (ex >= 0.10)
    except Exception:
        return False


# ---- Part A: 30 pilot 位點 ----
dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd
wsd = {x["key"]: x for x in json.load(open(f"{A}/ws_items.json"))}
v1d = {x["key"]: x for x in json.load(open(f"{A}/phylo_groups.json"))}
v2d = {x["key"]: x for x in json.load(open(f"{A}/phylo_v2_final.json"))}
rng = np.random.default_rng(20260622); rows = []
for key, ws in wsd.items():
    rd = dirmap.get(key)
    if not rd: continue
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    mrows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n")
    mi = {}; M = []
    for j, ln in enumerate(mrows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids) < MINSZ * 2: continue
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    maxk = min(6, n // 2)
    sk, ssc, slab = silhouette_bestk(sub, maxk)
    pperm = permanova_p(sub, list(slab), np.random.default_rng(20260622))
    # 4-gate group count: NO_CLEAR→1, else fine_k
    fg = 1 if ws["class"] == "NO_CLEAR_SPLIT" else ws["fine_k"]
    rows.append({"key": key, "n": n, "C": len(mi), "fine_class": ws["class"],
                 "sil_k": sk, "sil_perm_p": round(pperm, 3), "sil_perm_sig": pperm < 0.05,
                 "gate4_k": fg, "v1_k": v1d[key]["n_groups"], "v2_k": v2d[key]["v2_ngroups"]})
json.dump(rows, open(f"{A}/methods_compare_30loci.json", "w"), indent=1)
print("=== (A) 30 位點各方法群數 ===")
print(f"{'位點':>16} {'n':>4} {'C':>4} | {'sil_k(raw)':>10} {'sil+PERMp':>9} | {'4gate':>5} {'v1':>3} {'v2':>3}")
for r in sorted(rows, key=lambda z: -z["v2_k"]):
    print(f"{r['key']:>16} {r['n']:>4} {r['C']:>4} | {r['sil_k']:>10} {r['sil_perm_p']:>9} | {r['gate4_k']:>5} {r['v1_k']:>3} {r['v2_k']:>3}")
sil_multi = sum(1 for r in rows if r["sil_k"] >= 2)  # always all
sil_sig = sum(1 for r in rows if r["sil_perm_sig"])
g4_multi = sum(1 for r in rows if r["gate4_k"] >= 2)
v1_multi = sum(1 for r in rows if r["v1_k"] >= 2); v2_multi = sum(1 for r in rows if r["v2_k"] >= 2)
print(f"\n多群位點數/30: silhouette-raw={sil_multi}(全強制) | silhouette+PERMANOVA顯著={sil_sig} | 4-gate={g4_multi} | phylo-v1={v1_multi} | phylo-v2={v2_multi}")

# ---- Part B: 純噪音 FP ----
def pure_noise_P(n, C, rng, miss=0.25):
    rates = rng.uniform(0.1, 0.9, C); P = (rng.random((n, C)) < rates[None, :]).astype(float)
    P[rng.random((n, C)) < miss] = np.nan
    for cc in range(C):
        if np.sum(~np.isnan(P[:, cc])) < 2:
            ix = rng.choice(n, 2, replace=False); P[ix, cc] = (rng.random(2) < rates[cc]).astype(float)
    return P

print("\n=== (B) 純噪音偽陽性 FP（TRIALS=150, C=76=pilot中位, 25%缺值）— specificity ===")
print(f"{'n':>4} | {'sil-raw':>8} {'sil+PERM':>9} {'4gate':>7} {'phylo-v2':>9}")
TR = 150
for n in [20, 40, 80]:
    c_sr = c_sp = c_g4 = c_v2 = 0
    for t in range(TR):
        P = pure_noise_P(n, 76, rng); D = CR.bernoulli_dist(P); np.fill_diagonal(D, 0); D = np.maximum(D, D.T)
        sk, _, slab = silhouette_bestk(D, min(6, n // 2)); c_sr += 1  # always ≥2
        if permanova_p(D, list(slab), rng, 99) < 0.05: c_sp += 1
        if fourgate_real(D, P, rng): c_g4 += 1
    print(f"{n:>4} | {100*c_sr/TR:>7.0f}% {100*c_sp/TR:>8.1f}% {100*c_g4/TR:>6.1f}% {'~0%(校準)':>9}")
print("(phylo-v2 純噪音 FP 見 phylo_noise_calibration.py: C=76 → 0.0%；C=40 → ~0%)")
