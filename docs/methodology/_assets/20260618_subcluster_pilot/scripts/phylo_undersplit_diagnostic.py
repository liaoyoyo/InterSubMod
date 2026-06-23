#!/usr/bin/env python3
"""[診斷] 4 個『眼睛看到多群但 phylo 判 1 群』位點 — 為什麼拒切 + 是否該切 + v1/v2 是否變差。
每位點算: root 分裂 r=between/within vs SEP_MIN(1.3) vs null95(v1全域 / v2子群內); silhouette best_k+split;
候選 k=2/3 split 的: r / clusterboot Jaccard 穩定度(Hennig:>0.75穩) / 對齊 hp+allele CramerV / 各群甲基剖面(離散vs gradient).
判別: phylo 拒切是『gradient(對)』還是『離散但門檻太嚴(假陰性)』。渲染候選切割圖。"""
import os, csv, glob, json, sys
import numpy as np
from collections import Counter
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, fcluster

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUTD = f"{WT}/output/_ws_render"
FIGS = f"{A}/figs_undersplit"; os.makedirs(FIGS, exist_ok=True)
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


def cramerV(a, b):
    ca = sorted(set(a)); cb = sorted(set(b)); ia = {c: i for i, c in enumerate(ca)}; ib = {c: i for i, c in enumerate(cb)}
    tab = np.zeros((len(ca), len(cb)))
    for x, y in zip(a, b): tab[ia[x], ib[y]] += 1
    nn = tab.sum()
    if nn == 0: return 0.0
    row = tab.sum(1, keepdims=True); col = tab.sum(0, keepdims=True); exp = row * col / nn
    chi2 = np.nansum((tab - exp) ** 2 / np.where(exp > 0, exp, 1))
    k = min(len(ca), len(cb))
    return float(np.sqrt(chi2 / (nn * (k - 1)))) if k > 1 else 0.0


def clusterboot(D, k, rng, R=60, frac=0.8):
    """Hennig clusterboot: 重抽 frac read 重切 k 群, 各原群最佳 Jaccard 回收率均值。>0.75穩/0.6-0.75邊際/<0.6不穩。"""
    n = D.shape[0]; Z0, _ = CR.linkZ(D); base = fcluster(Z0, k, "maxclust")
    jac = []
    for _ in range(R):
        idx = np.where(rng.random(n) < frac)[0]
        if len(idx) < k + 1: continue
        Dn = D[np.ix_(idx, idx)]; Zn, _ = CR.linkZ(Dn); bl = fcluster(Zn, k, "maxclust")
        for g in set(base):
            orig = set(np.where(base == g)[0]) & set(idx)
            if not orig: continue
            best = 0
            for gg in set(bl):
                new = set(idx[np.where(bl == gg)[0]])
                u = orig | new
                if u: best = max(best, len(orig & new) / len(u))
            jac.append(best)
    return float(np.mean(jac)) if jac else 0.0


def root_split_diag(sub, P, rng):
    """root 分裂: r / null95(v1全域) / null95(v2子群內) / 各 verdict。"""
    Z, s, desc, ch, n = _tree(sub)
    root = 2 * n - 2; c1, c2 = ch[root]; S1, S2 = np.array(desc[c1]), np.array(desc[c2])
    r = _bw(s, S1, S2)
    sz = (len(S1), len(S2))
    # v1 全域 null
    C = P.shape[1]; gn = []
    for _ in range(RNULL):
        Pn = P.copy()
        for cc in range(C):
            col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
            if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
        Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); gn.append(np.maximum(Dn, Dn.T))
    n1 = [_bw(Dn, S1, S2) for Dn in gn]; n1 = [x for x in n1 if x is not None]; null95_v1 = float(np.percentile(n1, 95)) if n1 else 0
    # v2 子群內 null (root = 全體, 重分群)
    n2 = []
    for _ in range(RNULL):
        Pn = P.copy()
        for cc in range(C):
            col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
            if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
        Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
        _, sn, dn, cn, _ = _tree(Dn); nc1, nc2 = cn[2 * n - 2]
        n2.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
    n2 = [x for x in n2 if x is not None]; null95_v2 = float(np.percentile(n2, 95)) if n2 else 0
    return {"r": r, "sizes": sz, "sep_min": SEP_MIN, "null95_v1": round(null95_v1, 3), "null95_v2": round(null95_v2, 3),
            "pass_sepmin": (r is not None and r >= SEP_MIN),
            "v1_split": (r is not None and r >= SEP_MIN and r > null95_v1),
            "v2_split": (r is not None and r >= SEP_MIN and r > null95_v2), "r_round": round(r, 3) if r else None}


dirmap = {}
for mp in glob.glob(f"{OUTD}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): dirmap[part] = rd
ws = {x["key"]: x for x in json.load(open(f"{A}/ws_items.json"))}
rng = np.random.default_rng(20260622); out = []
for key in TARGETS:
    rd = dirmap.get(key)
    if not rd: print(f"skip {key}"); continue
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs = [int(c) for c in rows[0].split(",")[1:]]
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    Z, s, desc, ch, _ = _tree(sub)
    rsd = root_split_diag(sub, P, np.random.default_rng(20260622))
    # silhouette best_k + k=2/3 candidate splits
    try:
        from sklearn.metrics import silhouette_score
        sil = lambda lab: silhouette_score(sub, lab, metric="precomputed") if len(set(lab)) > 1 else -2
    except Exception:
        sil = lambda lab: -2
    cand = {}
    for k in [2, 3, 4]:
        if n < k * MINSZ: continue
        lab = fcluster(Z, k, "maxclust")
        if len(set(lab)) < 2: continue
        # split ratio for this k = mean over groups of between/within(rough): use overall pseudo via _bw on 2-coarsest? 用 k=2 root r 已有; 對 k 用 silhouette + 穩定 + 對齊
        hp = [reads[ids[i]]["hp"] for i in range(n)]; al = [reads[ids[i]]["alt_support"] for i in range(n)]
        cand[k] = {"sil": round(float(sil(list(lab))), 3),
                   "jaccard_stab": round(clusterboot(sub, k, np.random.default_rng(20260622)), 3),
                   "V_hp": round(cramerV(list(lab), hp), 3), "V_allele": round(cramerV(list(lab), al), 3),
                   "group_sizes": sorted(Counter(lab).values(), reverse=True)}
    # 離散 vs gradient: k=2 split 的兩群平均甲基剖面 + 每 read 在群心距離分佈雙峰性(用 root split 投影)
    lab2 = fcluster(Z, 2, "maxclust")
    g0 = [i for i in range(n) if lab2[i] == sorted(set(lab2))[0]]; g1 = [i for i in range(n) if lab2[i] == sorted(set(lab2))[1]]
    prof0 = np.nanmean(P[g0], axis=0); prof1 = np.nanmean(P[g1], axis=0)
    both = ~np.isnan(prof0) & ~np.isnan(prof1)
    prof_absdiff = float(np.nanmean(np.abs(prof0[both] - prof1[both]))) if both.sum() else 0
    prof_corr = float(np.corrcoef(prof0[both], prof1[both])[0, 1]) if both.sum() > 2 else 0
    rec = {"key": key, "n": n, "C": len(cpgs), "ws_class": ws[key]["class"], "ws_fine_k": ws[key]["fine_k"],
           "ws_excess": ws[key].get("excess"), "ws_align_axis": ws[key].get("align_axis"), "ws_align_V": ws[key].get("align_V"),
           "root_split": rsd, "candidates": cand,
           "k2_profile_absdiff": round(prof_absdiff, 3), "k2_profile_corr": round(prof_corr, 3)}
    out.append(rec)
    print(f"\n{key} n={n} C={len(cpgs)} ws_class={ws[key]['class']} ws_fine_k={ws[key]['fine_k']} excess={ws[key].get('excess')}", flush=True)
    print(f"  root r={rsd['r_round']} (sizes {rsd['sizes']}) vs SEP_MIN={SEP_MIN} | null95 v1={rsd['null95_v1']} v2={rsd['null95_v2']} → pass_sepmin={rsd['pass_sepmin']} v1_split={rsd['v1_split']} v2_split={rsd['v2_split']}")
    for k, c in cand.items():
        print(f"  k={k}: sil={c['sil']} jaccard_stab={c['jaccard_stab']} V_hp={c['V_hp']} V_allele={c['V_allele']} sizes={c['group_sizes']}")
    print(f"  k2 甲基剖面 |Δ|={prof_absdiff:.3f} corr={prof_corr:.3f} (高corr+低Δ=gradient; 低corr+高Δ=離散)")

    # render: heatmap ordered by tree + k=2/3 candidate colors + hp/allele
    dn = dendrogram(Z, orientation="left", no_plot=True); order = dn["leaves"][::-1]
    meth = np.array([P[i] for i in order]); dist = sub[np.ix_(order, order)].copy()
    np.fill_diagonal(dist, 0); dist[dist < 0] = np.nan
    CCOL = ["#dc2626", "#2563eb", "#16a34a", "#d97706"]
    l2 = [CCOL[fcluster(Z, 2, "maxclust")[i] % 4] for i in order]
    l3 = [CCOL[fcluster(Z, 3, "maxclust")[i] % 4] for i in order] if n >= 9 else l2
    hp_o = [reads[ids[i]]["hp"] for i in order]; al_o = [reads[ids[i]]["alt_support"] for i in order]
    HPC = H.HP_COL
    alc = H.ALT_COL
    sb = [("k=2切", l2), ("k=3切", l3), ("HP", [HPC.get(h, "#ddd") for h in hp_o]), ("ALT", [alc.get(a, "#ddd") for a in al_o])]
    mc, dc = H.mpl_cmaps(); nsb = len(sb)
    wr = [1.2] + [0.05] * nsb + [1.0, 0.14] + [0.05] * nsb + [1.0]
    fig = plt.figure(figsize=(11.5, 5.2)); gs = fig.add_gridspec(2, len(wr), width_ratios=wr, height_ratios=[1, 0.36], wspace=0.04, hspace=0.02)
    c = 0; axdn = fig.add_subplot(gs[0, c]); c += 1
    dendrogram(Z, orientation="left", ax=axdn, no_labels=True, above_threshold_color="#888")
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA 樹", fontsize=8); [axdn.spines[s2].set_visible(False) for s2 in axdn.spines]
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
    k2 = cand.get(2, {})
    txt = (f"v2/v1 判定: 1 群（root r={rsd['r_round']} < SEP_MIN {SEP_MIN}? {'是→拒切' if not rsd['pass_sepmin'] else 'r≥1.3 但 < null95 '+str(rsd['null95_v2'])+'→拒切'}）"
           f"  |  k=2: silhouette={k2.get('sil')} Jaccard穩定={k2.get('jaccard_stab')} V(hp)={k2.get('V_hp')} V(allele)={k2.get('V_allele')}"
           f"  |  4-gate: {ws[key]['class']}(fine_k={ws[key]['fine_k']}, excess={ws[key].get('excess')})"
           f"  |  k2甲基剖面 corr={prof_corr:.2f} |Δ|={prof_absdiff:.3f}")
    axl.text(0.0, 0.6, txt, fontsize=7.4, transform=axl.transAxes, wrap=True)
    fig.suptitle(f"{key} n={n} CpG={len(cpgs)} — 為什麼 phylo 判 1 群（眼睛看到多群）", fontsize=10, y=1.0)
    fn = f"{FIGS}/us_{key}.png"; fig.savefig(fn, dpi=108, bbox_inches="tight"); plt.close(fig)
    rec["png"] = f"figs_undersplit/us_{key}.png"
json.dump(out, open(f"{A}/phylo_undersplit_diagnostic.json", "w"), ensure_ascii=False, indent=1)
print(f"\nDONE {len(out)} 位點 → phylo_undersplit_diagnostic.json + figs_undersplit/")
