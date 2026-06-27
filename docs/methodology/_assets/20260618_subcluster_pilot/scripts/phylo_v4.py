#!/usr/bin/env python3
"""[v4 定稿] phylo v3.1 + 3 改善(用戶回饋):
 ① other 群: 殘留 outlier ≥MINSZ → 標 'other'(residual 紀錄, 非 validated cluster)
 ② instability flag: K=5 seeds 跑 coarse, 群數不一致 → unstable + 報範圍(min-max)
 ③ coarse/fine 兩層: coarse=null95(主 verdict 嚴) / fine=null90(候選更細群, 低信心)
回傳 per-locus {coarse_ng, fine_ng, n_other, unstable, ng_range, coarse_labels...}。
本檔: 跑 30 pilot 比 v3(coarse) → 看 fine 多切幾個 / 幾個 unstable / 幾個有 other。重用快取。"""
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


def v4_label(sub, P, rng, null_pct=95, want_other=True):
    """v3.1 + null_pct 參數 + other 群。回傳 per-read labels(含 'other'/'outlier')。"""
    Z, s, desc, ch, n = _tree(sub)

    def subnull(leaves):
        S = np.array(leaves); m = len(S); Ps = P[S]; ns = []
        for _ in range(RNULL):
            Pn = Ps.copy()
            for cc in range(Pn.shape[1]):
                col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
            Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
            _, sn, dn, cn, _ = _tree(Dn); nc1, nc2 = cn[2 * m - 2]
            ns.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
        ns = [x for x in ns if x is not None]; return float(np.percentile(ns, null_pct)) if ns else 0

    def descend(node):
        cur = node; q = []
        while cur in ch:
            c1, c2 = ch[cur]; s1, s2 = len(desc[c1]), len(desc[c2])
            if min(s1, s2) >= MINSZ: return cur, q
            sm, bg = (c1, c2) if s1 < s2 else (c2, c1); q += desc[sm]; cur = bg
        return None, q

    def test(bn):
        c1, c2 = ch[bn]; r = _bw(s, np.array(desc[c1]), np.array(desc[c2]))
        if r is None or r < SEP_MIN: return False
        return r > subnull(desc[bn])

    lab = [None] * n

    def rec(node, label):
        lv = desc[node]
        if len(lv) < 2 * MINSZ:
            for i in lv: lab[i] = label; return
        bn, q = descend(node)
        if bn is not None and test(bn):
            for i in q: lab[i] = "outlier"
            c1, c2 = ch[bn]; big, sm = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            rec(big, label + "-1"); rec(sm, label + "-2")
        else:
            for i in lv: lab[i] = label
    rec(2 * n - 2, "1")
    lab = [l if l else "outlier" for l in lab]
    sm = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    lab = [("outlier" if l in sm else l) for l in lab]
    # ① other 群: 殘留 outlier ≥MINSZ → 'other'
    if want_other and sum(1 for l in lab if l == "outlier") >= MINSZ:
        lab = [("other" if l == "outlier" else l) for l in lab]
    return lab


def ng_of(lab):
    return len(set(l for l in lab if l not in ("outlier", "other")))


MODAL_CONF = 0.7  # modal_fraction 門檻: <0.7=真分裂unstable, ≥0.7=多數穩定(單seed偏離不誤報)
HIDDEN_HET = 0.30  # n_other/n > 此 → 隱藏異質性警告(coarse_ng=1 但大量 other 不可報「單群」)


def analyze_v4(sub, P, base=20260622, K=10):
    """② instability redesign: K seeds → modal_ng(robust verdict, 非單seed) + modal_fraction(信心)。
    unstable=modal_frac<0.7(真分裂); ≥0.7=多數穩定報 modal。修 seed-敏感(revalidate C + 用戶 chr20:42981498 modal=3)。"""
    cs = [ng_of(v4_label(sub, P, np.random.default_rng(base + k * 101), 95)) for k in range(K)]
    cnt = Counter(cs); modal_ng, mc = cnt.most_common(1)[0]; modal_frac = mc / K
    coarse_ng = modal_ng                          # robust modal(非單一 base seed)
    unstable = modal_frac < MODAL_CONF            # 只真分裂才 unstable, 非單seed偏離
    ng_range = (min(cs), max(cs))
    # coarse_lab: 取第一個給出 modal_ng 的 seed(代表性), 否則 base
    coarse_lab = None
    for k in range(K):
        l = v4_label(sub, P, np.random.default_rng(base + k * 101), 95)
        if ng_of(l) == modal_ng: coarse_lab = l; break
    if coarse_lab is None: coarse_lab = v4_label(sub, P, np.random.default_rng(base), 95)
    # ③ fine(null90, 候選低信心)
    fine_lab = v4_label(sub, P, np.random.default_rng(base), 90); fine_ng = ng_of(fine_lab)
    n = len(coarse_lab); n_other = sum(1 for l in coarse_lab if l == "other"); n_out = sum(1 for l in coarse_lab if l == "outlier")
    # ④ 隱藏異質性: n_other 過大 → coarse_ng 不可單獨報「單群」
    hidden_het_warn = (n_other / n > HIDDEN_HET) if n else False
    # merge-back: coarse-'other' read 在 fine(更敏感)是否成群 → 高=other 藏結構
    oth_idx = [i for i, l in enumerate(coarse_lab) if l == "other"]
    other_recovered_fine = round(sum(1 for i in oth_idx if fine_lab[i] not in ("other", "outlier")) / len(oth_idx), 2) if oth_idx else 0.0
    return {"coarse_ng": coarse_ng, "modal_frac": round(modal_frac, 2), "fine_ng": fine_ng,
            "n_other": n_other, "n_outlier": n_out, "unstable": unstable, "ng_range": ng_range, "seeds_ng": cs,
            "hidden_het_warn": hidden_het_warn, "other_recovered_fine": other_recovered_fine,
            "coarse_lab": coarse_lab, "fine_lab": fine_lab}


if __name__ == "__main__":
    dirmap = {}
    for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
        rd = mp.rsplit("/distance/", 1)[0]
        for part in rd.split("/"):
            if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd
    ws = {x["key"]: x for x in json.load(open(f"{A}/ws_items.json"))}
    v3 = {x["key"]: x["v2_ngroups"] for x in json.load(open(f"{A}/phylo_v3_render.json"))}
    out = []
    print(f"{'位點':>16} {'n':>4} | v3/coarse | fine | other | unstable(range)")
    for key in ws:
        rd = dirmap.get(key)
        if not rd: continue
        reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
        dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
        rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n")
        mi = {}; M = []
        for j, ln in enumerate(rows[1:]):
            q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
        M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
        ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
        if len(ids) < MINSZ * 2: continue
        sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
        ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
        a = analyze_v4(sub, P)
        flag = ""
        if a["fine_ng"] > a["coarse_ng"]: flag += f" fine+{a['fine_ng']-a['coarse_ng']}"
        if a["n_other"] > 0: flag += f" other={a['n_other']}"
        if a["hidden_het_warn"]: flag += " 🔴hidden-het"
        if a["unstable"]: flag += f" ⚠unstable{a['ng_range']}(modal {a['coarse_ng']}@{a['modal_frac']})"
        out.append({"key": key, "n": n, "v3": v3.get(key), **{k: a[k] for k in ["coarse_ng", "modal_frac", "fine_ng", "n_other", "n_outlier", "unstable", "ng_range", "seeds_ng", "hidden_het_warn", "other_recovered_fine"]}})
        print(f"{key:>16} {n:>4} | modal {a['coarse_ng']}@{a['modal_frac']} | fine {a['fine_ng']} | other {a['n_other']} |{flag}")
    json.dump(out, open(f"{A}/phylo_v4_pilot.json", "w"))
    nfine = sum(1 for o in out if o["fine_ng"] > o["coarse_ng"]); noth = sum(1 for o in out if o["n_other"] > 0); nuns = sum(1 for o in out if o["unstable"])
    nhh = sum(1 for o in out if o["hidden_het_warn"])
    print(f"\n30 位點: coarse(modal) 多群 {sum(1 for o in out if o['coarse_ng']>=2)} | fine>coarse {nfine} | other {noth} | hidden-het警告 {nhh} | unstable(modal<0.7真分裂) {nuns}")
    chg = [(o['key'], o['v3'], o['coarse_ng']) for o in out if o['coarse_ng'] != o['v3']]
    print(f"coarse(modal) vs v3(單seed) 變動 {len(chg)} 個: {chg} (modal 較 robust, 預期部分變)")
