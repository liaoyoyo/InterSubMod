#!/usr/bin/env python3
"""build_structure_result_data.py — 整合結構結果數據(2026-07-09)。
產出 JSON:①7樣本 region 層 c 分布 ②VAF 峰值+直方圖(頻率譜) ③巢狀 clone→subclone 父子統計
④HCC1395 精選代表區(多候選樹 + 每節點 VAF 標註)。供 build_structure_result_html.py 注入(§13.0 反捏造)。
c(region層)=該區跨家族(1/2/3) distinct ALT 組合數=somatic clone 數;VAF from col_coverage_by_hp+positions 對齊。
"""
import json, glob
from collections import Counter, defaultdict

MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260709_structure_result_data.json"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
MIN_COV, MIN_ALT = 6, 2

def wd(s): return PILOT if s == "HCC1395" else f"{MSROOT}/{s}"
def rv_path(s): return f"{PILOT}/layered_region_view_HCC1395.json" if s == "HCC1395" else f"{MSROOT}/{s}/layered_region_view_{s}.json"

def mlhp_lookup(s):
    lk = {}
    for f in sorted(glob.glob(f"{wd(s)}/mlhp_part_*.json")):
        for g in json.load(open(f))["groups"]:
            if g.get("n_sSNV", 0) >= 2:
                lk[(g["chrom"], g["start"])] = g
    return lk

def peaks_from_hist(vafs, binw=0.05):
    if not vafs: return [], []
    nb = int(round(1.0 / binw))
    h = [0] * nb
    for v in vafs:
        h[min(nb - 1, int(v / binw))] += 1
    centers = [round((i + 0.5) * binw, 3) for i in range(nb)]
    sm = [(h[max(0, i-1)] + h[i] + h[min(nb-1, i+1)]) / 3 for i in range(nb)]
    thr = 0.03 * sum(h)
    pk = []
    for i in range(nb):
        l = sm[i-1] if i > 0 else -1
        r = sm[i+1] if i < nb-1 else -1
        if sm[i] >= l and sm[i] >= r and h[i] >= thr:
            pk.append({"vaf": centers[i], "n": h[i]})
    return pk, h

def tree_shape(edges):
    if not edges: return "single"
    ch = defaultdict(list); nodes = set()
    for p, c in edges:
        ch[p].append(c); nodes.add(p); nodes.add(c)
    if len([n for n in nodes if n != "ROOT"]) <= 1: return "single"
    cc = {p: len(cs) for p, cs in ch.items()}
    if any(v >= 2 for p, v in cc.items() if p != "ROOT"): return "branched"
    if cc.get("ROOT", 0) >= 2: return "star"
    return "linear"

def geno_of(node):
    """節點 genotype 字串(去 H_ 前綴);ROOT→全 R 由呼叫端補。"""
    return node[2:] if str(node).startswith("H_") else str(node)

result = {"generated": "2026-07-09", "samples": {}, "curated_regions": []}

for s in SAMPLES:
    lk = mlhp_lookup(s)
    # ① region 層 c 分布
    cdist = Counter(); unresolved = 0; maxc = 0
    # ② VAF 直方圖 + 峰值
    vafs = []
    # ③ 巢狀父子
    nested = sister = 0; anc_vafs = []; der_vafs = []; anc_gt = 0; npair = 0
    for (chrom, start), g in lk.items():
        pbh = g.get("populations_by_hp", {}) or {}
        cbh = g.get("col_coverage_by_hp", {}) or {}
        positions = g.get("positions", [])
        # region 層 c
        alt = set()
        for fam in ("1", "2", "3"):
            for gt in (pbh.get(fam, {}) or {}):
                if "A" in gt: alt.add((fam, gt))
        c = len(alt)
        if c == 0: unresolved += 1
        else: cdist[c] += 1; maxc = max(maxc, c)
        # VAF 池
        for fam in ("1", "2"):
            for pos, (nr, na) in (cbh.get(fam, {}) or {}).items():
                if na >= MIN_ALT and nr + na >= MIN_COV:
                    vafs.append(na / (nr + na))
        # 巢狀父子
        for fam in ("1", "2"):
            pops = pbh.get(fam, {}) or {}; cc = cbh.get(fam, {}) or {}
            alts = [gt for gt in pops if "A" in gt]
            if len(alts) != 2: continue
            a, b = alts
            sa = {i for i, ch in enumerate(a) if ch == "A"}
            sb = {i for i, ch in enumerate(b) if ch == "A"}
            if sa < sb or sb < sa:
                nested += 1
                parent, child = (a, b) if sa < sb else (b, a)
                anc_idx = [i for i, ch in enumerate(parent) if ch == "A"]
                der_idx = [i for i, ch in enumerate(child) if ch == "A" and parent[i] != "A"]
                def vaf_at(idx):
                    vs = []
                    for i in idx:
                        if i < len(positions):
                            p = str(positions[i])
                            if p in cc:
                                nr, na = cc[p]
                                if nr + na >= MIN_COV: vs.append(na / (nr + na))
                    return sum(vs) / len(vs) if vs else None
                av, dv = vaf_at(anc_idx), vaf_at(der_idx)
                if av is not None and dv is not None:
                    anc_vafs.append(av); der_vafs.append(dv); npair += 1
                    if av > dv: anc_gt += 1
            else:
                sister += 1
    pk, hist = peaks_from_hist(vafs)
    def med(xs): return round(sorted(xs)[len(xs)//2], 3) if xs else None
    result["samples"][s] = {
        "c_dist": {str(k): cdist.get(k, 0) for k in range(1, 10)},
        "c_ge1": sum(cdist.values()), "unresolved": unresolved, "maxC": maxc,
        "c_ge2": sum(v for k, v in cdist.items() if k >= 2),
        "vaf_hist": hist, "vaf_peaks": pk, "n_somatic_pos": len(vafs),
        "nested": nested, "sister": sister, "n_pair": npair, "anc_gt_der": anc_gt,
        "anc_med_vaf": med(anc_vafs), "der_med_vaf": med(der_vafs),
    }
    print(f"{s}: c≥1={sum(cdist.values())} c≥2={sum(v for k,v in cdist.items() if k>=2)} maxC={maxc} 峰={[p['vaf'] for p in pk]} 巢狀={nested} 祖先中位VAF={med(anc_vafs)}→衍生{med(der_vafs)}")

# ④ HCC1395 精選代表區(多候選樹 + VAF 標註)
s = "HCC1395"
lk = mlhp_lookup(s)
rv = json.load(open(rv_path(s)))
buckets = {"single": [], "nested": [], "sister": [], "complex": [], "multitree": []}
for r in rv["regions"]:
    g = lk.get((r["chrom"], r["start"]))
    if not g: continue
    positions = g.get("positions", [])
    pbh = g.get("populations_by_hp", {}) or {}
    cbh = g.get("col_coverage_by_hp", {}) or {}
    for L in r["lineages"]:
        fam = L["family"]
        if fam not in ("1", "2") or L.get("capped"): continue
        trees = L.get("trees") or []
        if not trees: continue
        pops = pbh.get(fam, {}) or {}
        cc = cbh.get(fam, {}) or {}
        n_alt = sum(1 for gt in pops if "A" in gt)
        if n_alt == 0: continue
        # per-position VAF
        posvaf = {}
        for i, p in enumerate(positions):
            ps = str(p)
            if ps in cc:
                nr, na = cc[ps]
                if nr + na >= MIN_COV: posvaf[i] = round(na / (nr + na), 3)
        if len(posvaf) < len(positions): continue  # 要全位點都有 VAF 才乾淨展示
        shapes = {tree_shape(t["edges"]) for t in trees}
        shape = next(iter(shapes)) if len(shapes) == 1 else "mixed"
        # 每棵樹的節點 VAF 標註
        def annotate(edges):
            nodes = set()
            for p, c in edges: nodes.add(p); nodes.add(c)
            node_info = {}
            for nd in nodes:
                if nd == "ROOT":
                    node_info[nd] = {"geno": "R" * len(positions), "hidden": False, "root": True, "muts": []}
                    continue
                gg = geno_of(nd)
                muts = [{"pos": positions[i], "vaf": posvaf.get(i)} for i, chh in enumerate(gg) if chh == "A"]
                node_info[nd] = {"geno": gg, "hidden": str(nd).startswith("H_"), "root": False, "muts": muts}
            return node_info
        rec = {
            "chrom": r["chrom"], "start": r["start"], "family": fam,
            "c": n_alt, "shape": shape, "n_trees": L.get("n_trees"), "n_stored": len(trees),
            "positions": positions, "posvaf": {str(i): posvaf[i] for i in posvaf},
            "populations": pops,
            "trees": [{"edges": t["edges"], "shape": tree_shape(t["edges"]),
                       "n_hidden": t.get("n_hidden", 0), "nodes": annotate(t["edges"])}
                      for t in trees[:8]],
        }
        # 分桶
        if n_alt == 1 and shape == "single": buckets["single"].append(rec)
        elif n_alt >= 3: buckets["complex"].append(rec)
        elif n_alt == 2 and shape == "linear": buckets["nested"].append(rec)
        elif n_alt == 2 and shape in ("branched", "star"): buckets["sister"].append(rec)
        if (L.get("n_trees") or 0) >= 3 and len(trees) >= 3: buckets["multitree"].append(rec)

# 每桶挑代表(nested 優先挑 VAF 梯度明顯者:祖先高衍生低)
def vaf_gradient(rec):
    vs = sorted(rec["posvaf"].values())
    return (vs[-1] - vs[0]) if len(vs) >= 2 else 0
picked = []
for cat, want in [("nested", 3), ("sister", 2), ("complex", 2), ("single", 2), ("multitree", 3)]:
    lst = buckets[cat]
    if cat == "nested":
        lst = sorted(lst, key=lambda r: -vaf_gradient(r))
    elif cat == "complex":
        lst = sorted(lst, key=lambda r: -r["c"])
    elif cat == "multitree":
        lst = sorted(lst, key=lambda r: -(r["n_trees"] or 0))
    seen = {(r["chrom"], r["start"], r["family"]) for r in picked}
    for r in lst:
        key = (r["chrom"], r["start"], r["family"])
        if key in seen: continue
        r2 = dict(r); r2["category"] = cat
        picked.append(r2); seen.add(key)
        if sum(1 for x in picked if x["category"] == cat) >= want: break
result["curated_regions"] = picked

json.dump(result, open(OUT, "w"), ensure_ascii=False)
print(f"\n精選代表區 {len(picked)} 個(nested/sister/complex/single/multitree):")
for r in picked:
    print(f"  [{r['category']:9}] {r['chrom']}:{r['start']} fam{r['family']} c={r['c']} shape={r['shape']} n_trees={r['n_trees']} VAF={sorted(r['posvaf'].values())}")
print(f"\n→ 寫出 {OUT}")
