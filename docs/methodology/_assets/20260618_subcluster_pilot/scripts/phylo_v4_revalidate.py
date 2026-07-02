#!/usr/bin/env python3
"""[v4 再驗] 對抗式自檢 v4 三改善是否健全:
 A. fine層(null90) 噪音FP × CpG — fine 是否在低CpG過切(最關鍵, fine是最鬆的改善)
 B. other群是否藏真結構 — other-read-set 是否比 locus 中位更緊密(緊=可能漏掉的群,壞) vs 散(對,殘留)
 C. instability 穩健性 — K=10 確認 unstable 旗標一致(非 K=5 噪音)
純讀快取 + 模擬。輸出 phylo_v4_revalidate.json。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H  # noqa
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
from phylo_v4 import v4_label, ng_of
WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
MINSZ = 3
out = {}

# ---- A. fine(null90) vs coarse(null95) 噪音 FP × CpG ----
def pn(n, C, rng, miss=0.25):
    rt = rng.uniform(0.1, 0.9, C); Q = (rng.random((n, C)) < rt[None, :]).astype(float)
    Q[rng.random((n, C)) < miss] = np.nan
    for cc in range(C):
        if np.sum(~np.isnan(Q[:, cc])) < 2:
            ix = rng.choice(n, 2, replace=False); Q[ix, cc] = (rng.random(2) < rt[cc]).astype(float)
    return Q

print("=== A. 噪音 FP: coarse(null95) vs fine(null90) × CpG (n=40, TRIALS=120) ===")
print(f"{'C':>4} | {'coarse95 FP':>11} | {'fine90 FP':>10}")
A_rows = []
rng = np.random.default_rng(11)
for C in [10, 20, 40, 76]:
    fc = ff = 0; TR = 120
    for _ in range(TR):
        Q = pn(40, C, rng); D = CR.bernoulli_dist(Q); np.fill_diagonal(D, 0); D = np.maximum(D, D.T)
        if ng_of(v4_label(D, Q, rng, 95, want_other=False)) >= 2: fc += 1
        if ng_of(v4_label(D, Q, rng, 90, want_other=False)) >= 2: ff += 1
    A_rows.append({"C": C, "coarse95_fp": round(100 * fc / TR, 1), "fine90_fp": round(100 * ff / TR, 1)})
    print(f"{C:>4} | {100*fc/TR:>10.1f}% | {100*ff/TR:>9.1f}%")
out["A_noise_fp"] = A_rows

# ---- load pilot ----
dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd
ws = list(json.load(open(f"{A}/ws_items.json")))


def load(key):
    rd = dirmap.get(key)
    if not rd: return None
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n")
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids) < MINSZ * 2: return None
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]
    return ids, sub, np.array([M[mi[x]] for x in ids]), reads

# ---- B. other 群是否藏真結構 ----
print("\n=== B. other 群凝聚度 (other-set 內均距 / locus 中位; <1=緊密可能漏群, ≥1=散=對) ===")
B_rows = []
for it in ws:
    key = it["key"]; r = load(key)
    if not r: continue
    ids, sub, P, reads = r
    lab = v4_label(sub, P, np.random.default_rng(20260622), 95)
    oth = [i for i, l in enumerate(lab) if l == "other"]
    if len(oth) < MINSZ: continue
    n = sub.shape[0]; tri = sub[np.triu_indices(n, 1)]; med = np.median(tri[tri >= 0])
    osub = sub[np.ix_(oth, oth)]; ot = osub[np.triu_indices(len(oth), 1)]; ot = ot[ot >= 0]
    ratio = float(np.mean(ot) / med) if med > 0 and ot.size else None
    # other-set 自己跑 v4 看是否成群
    osubP = P[oth]; oth_ng = ng_of(v4_label(osub, osubP, np.random.default_rng(20260622), 95, want_other=False)) if len(oth) >= 2 * MINSZ else 0
    B_rows.append({"key": key, "n_other": len(oth), "other_tightness": round(ratio, 2) if ratio else None, "other_self_ng": oth_ng})
    print(f"  {key}: other n={len(oth)} 凝聚度={round(ratio,2) if ratio else '?'} ({'緊密⚠' if ratio and ratio<0.85 else '散=殘留✓'}) other自跑={oth_ng}群")
out["B_other_coherence"] = B_rows

# ---- C. instability K=10 穩健 ----
print("\n=== C. instability K=10 (確認 unstable 旗標一致, 非 K=5 噪音) ===")
v4p = {o["key"]: o for o in json.load(open(f"{A}/phylo_v4_pilot.json"))}
C_rows = []
for it in ws:
    key = it["key"]; r = load(key)
    if not r: continue
    ids, sub, P, reads = r
    ngs = [ng_of(v4_label(sub, P, np.random.default_rng(7000 + k * 313), 95)) for k in range(10)]
    uns10 = len(set(ngs)) > 1; uns5 = v4p.get(key, {}).get("unstable", None)
    if uns10 or uns5:
        agree = "✓一致" if (bool(uns10) == bool(uns5)) else "✗不一致"
        C_rows.append({"key": key, "k10_range": [min(ngs), max(ngs)], "k10_unstable": uns10, "k5_unstable": uns5, "agree": bool(uns10) == bool(uns5)})
        print(f"  {key}: K10 range={[min(ngs),max(ngs)]} unstable={uns10} (K5={uns5}) {agree}")
out["C_instability_k10"] = C_rows

json.dump(out, open(f"{A}/phylo_v4_revalidate.json", "w"), ensure_ascii=False, indent=1)
# 摘要
fine_safe = all(r["fine90_fp"] <= max(2.0, r["coarse95_fp"] + 1) for r in A_rows if r["C"] >= 20)
oth_ok = all((b["other_tightness"] is None or b["other_tightness"] >= 0.85) and b["other_self_ng"] == 0 for b in B_rows)
inst_ok = all(c["agree"] for c in C_rows)
print(f"\n=== 摘要 ===")
print(f"A fine層 C≥20 不顯著過切於 coarse? {fine_safe} (低CpG C=10 fine FP={A_rows[0]['fine90_fp']}% 需留意)")
print(f"B other群皆散/不自成群(真殘留)? {oth_ok}")
print(f"C instability K5↔K10 一致? {inst_ok}")
