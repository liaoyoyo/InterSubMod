#!/usr/bin/env python3
"""[驗證] phylo v2 最終結果 + 各群對齊 — per-subgroup 重分群 null + RNULL=40(已修 double-dip)。
輸出 phylo_v2_final.json: 每位點 v2 群數/標籤/各群 dominant hp+allele + 對 v1 的變化。供 doc 引用。
重用 _ws_render 快取矩陣(免 binary)。數字 L1(從快取矩陣重算)。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H  # noqa
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
MINSZ = 3; SEP_MIN = 1.3; RNULL = 40


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
        c1, c2 = ch[node]
        S1, S2 = np.array(desc[c1]), np.array(desc[c2])
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
    return [("outlier" if l in sm else l) for l in lab]


dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd
v1 = {o["key"]: o for o in json.load(open(f"{A}/phylo_groups.json"))}
items = json.load(open(f"{A}/ws_items.json")); out = []
for it in items:
    key = it["key"]; rd = dirmap.get(key)
    if not rd: continue
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); nC = len(rows[0].split(",")) - 1
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids) < MINSZ * 2: continue
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    lab = phylo_v2(sub, P, np.random.default_rng(20260622))
    labct = Counter(l for l in lab if l != "outlier"); nout = sum(1 for l in lab if l == "outlier")

    def dom(idxs, field):
        c = Counter(reads[ids[i]][field] for i in idxs); return c.most_common(1)[0] if c else ("-", 0)
    galign = {}
    for L in labct:
        idxs = [i for i in range(n) if lab[i] == L]
        hpd = dom(idxs, "hp"); ald = dom(idxs, "alt_support")
        galign[L] = {"n": len(idxs), "hp": f"{hpd[0]}({100*hpd[1]//len(idxs)}%)", "allele": f"{ald[0]}({100*ald[1]//len(idxs)}%)"}
    v1g = v1.get(key, {}).get("n_groups", "?")
    out.append({"key": key, "n": n, "C": nC, "fine": it["fine_conf"], "v1_ngroups": v1g,
                "v2_ngroups": len(labct), "labels": dict(labct), "align": galign, "n_outlier": nout,
                "changed": v1g != len(labct)})
json.dump(out, open(f"{A}/phylo_v2_final.json", "w"), indent=1)
print("=== v2 多群案例(對齊) ===")
for o in sorted(out, key=lambda z: -z["v2_ngroups"]):
    if o["v2_ngroups"] >= 2:
        parts = " || ".join(f"{L}:n{o['align'][L]['n']} hp={o['align'][L]['hp']} alle={o['align'][L]['allele']}" for L in sorted(o["align"]))
        ch = "  ←v1變" if o["changed"] else ""
        print(f"{o['key']} n={o['n']} C={o['C']} v1={o['v1_ngroups']}群→v2={o['v2_ngroups']}群{ch} | {parts}")
print(f"\n總: v2 多群 {sum(1 for o in out if o['v2_ngroups']>=2)}/{len(out)} | 變動 {sum(1 for o in out if o['changed'])} | v2 平均 {np.mean([o['v2_ngroups'] for o in out]):.2f} 群")
# 對齊型態: 每個 v2 多群位點是否各群分到不同 allele(REF/ALT)或不同 hp = germline 軸(cis-ASM) vs 同 germline 內多群(subclone候選)
print("\n=== 對齊判讀(v2 多群) ===")
for o in out:
    if o["v2_ngroups"] < 2: continue
    al = o["align"]; alleles = {al[L]["allele"].split("(")[0] for L in al}; hps = {al[L]["hp"].split("(")[0] for L in al}
    same_germ = len(alleles) == 1 and len(hps) == 1
    print(f"{o['key']}: {'🔴同germline內多群(subclone候選)' if same_germ else '✅germline軸分裂(cis-ASM): allele='+str(alleles)+' hp='+str(hps)}")
