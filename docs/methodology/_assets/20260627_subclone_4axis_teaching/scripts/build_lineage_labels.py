#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lineage 標籤產生器：把每個區域的克隆樹節點標成 HP{h}-{b1}-{b2}... (Dewey 路徑)。
定義（見 20260628_lineage_label_definition_01.md）：
  h = germline haplotype 根 (H1/H2/H3-unphased)；路徑 = somatic 樹從根往下的 lineage。
  nested(祖先→後代)=往下加層；sibling(互斥)=同層新分支；co_linked=同節點。
  分支編號 = VAF 遞減 (1=最高；CN-altered/tie 標 ?)。
產出：
  (1) genome-wide per-region 節點標籤 + situation tier + 信心旗標 → lineage_labels_regions.tsv
  (2) chr17 per-read 驗證 (應得 HP1 / HP1-1 / HP1-1-1) → 印出
用法：python3 build_lineage_labels.py
"""
import json, os, csv
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

def hp_letter(v):
    return {"1-1": "H1", "2-1": "H2", "3": "H3?"}.get(v, "H?")

# locus -> HP（從 pair lists 建）
locus_hp = {}
for c in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl", "independent"):
    for hp in ("sameHP", "diffHP"):
        p = os.path.join(DATA, "lists", f"{c}_{hp}.tsv")
        if not os.path.exists(p): continue
        for d in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"):
            locus_hp[(d["chrom"], int(d["pos_a"]))] = d["hp_a"]
            locus_hp[(d["chrom"], int(d["pos_b"]))] = d["hp_b"]

def assign_paths(nodes, nested_edges, sibling_pairs, vafd, hletter):
    """回 {node_pos: 'H{h}-b1-b2...'}。nested_edges=(ancestor,descendant)。"""
    children = defaultdict(list); has_parent = set()
    for a, dsc in nested_edges:
        children[a].append(dsc); has_parent.add(dsc)
    def vof(n): return vafd.get(n, -1.0)
    roots = sorted([n for n in nodes if n not in has_parent], key=lambda n: -vof(n))
    label = {}
    def walk(n, path):
        label[n] = f"{hletter}-" + "-".join(map(str, path))
        kids = sorted(children.get(n, []), key=lambda k: -vof(k))
        for j, k in enumerate(kids, 1):
            if k not in label:  # 防多親 DAG 重入
                walk(k, path + [j])
    for i, r in enumerate(roots, 1):
        walk(r, [i])
    return label

def region_situation(r):
    if r["has_cycle"]: return "incompatible(cycle)"
    span = r["span"]
    pops = r.get("populations") or {}
    multi = sum(1 for k in pops if "A" in k) >= 1 and len([k for k in pops if "A" in k]) >= 1 and len(pops) >= 2
    if span <= 34000 and multi: return "A_single_molecule"
    if span <= 76000: return "B_spannable_pairwise"
    return "C_chained_statistical"

# ---------- (1) genome-wide per-region ----------
RI = json.load(open(os.path.join(DATA, "sm_region_integration.json"), encoding="utf-8"))
rows = []; tier_cnt = Counter()
for r in RI["regions"]:
    if r["n_sSNV"] < 2: continue
    ne = [tuple(e) for e in r.get("nested_edges", [])]
    sp = [tuple(e) for e in r.get("sibling_pairs", [])]
    nodes = set()
    for a, b in ne: nodes.add(a); nodes.add(b)
    for a, b in sp: nodes.add(a); nodes.add(b)
    if not nodes: continue
    # haplotype（節點所屬，取眾數；同 region 樹為 same-HP）
    hps = Counter(hp_letter(locus_hp.get((p.split(":")[0], int(p.split(":")[1])))) for p in nodes)
    hletter = hps.most_common(1)[0][0]
    vafd = r.get("vaf", {})
    labels = assign_paths(nodes, ne, sp, vafd, hletter)
    tier = region_situation(r)
    tier_cnt[tier] += 1
    cn_clean = r["cn"] in ("loh", "neutral")
    rows.append({"region": r["region"], "n_sSNV": r["n_sSNV"], "span": r["span"], "cn": r["cn"],
                 "tree_shape": r["tree_shape"], "haplotypes": ",".join(sorted(hps)),
                 "n_nodes": len(nodes), "max_depth": r["max_depth"],
                 "situation_tier": tier, "vaf_order_ok": "yes" if cn_clean else "CN-altered(?)",
                 "labels": ";".join(f"{p.split(':')[1]}={labels[p]}" for p in sorted(labels))})

with open(os.path.join(DATA, "lineage_labels_regions.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ---------- (2) chr17 per-read 驗證 ----------
chr17 = json.load(open(os.path.join(DATA, "chr17_subclone_data.json"), encoding="utf-8"))
A_POS, B1, B2 = "48365089", "48362515", "48365161"  # α 祖先, β1/β2 後代(co_linked)
def read_label(geno):
    a = geno.get(A_POS) == "ALT"; b = (geno.get(B1) == "ALT") or (geno.get(B2) == "ALT")
    if a and b: return "H1-1-1"   # α+β 後代
    if a: return "H1-1"           # α-only
    return "H1"                   # germline 根
rl = Counter(read_label(rd["geno"]) for rd in chr17["reads"])

print("=== genome-wide per-region 標籤 ===")
print(f"  輸出 {len(rows)} 區 → lineage_labels_regions.tsv")
print(f"  situation tier: {dict(tier_cnt)}")
print("  範例(前 3 結構區):")
for r in rows[:3]:
    print(f"    {r['region']} [{r['tree_shape']},{r['situation_tier']},{r['haplotypes']}] {r['labels'][:80]}")
print("=== chr17 per-read 驗證 (應 H1/H1-1/H1-1-1) ===")
for k in sorted(rl): print(f"    {k}: {rl[k]} reads")
print(f"  期望: H1≈9(L0) / H1-1≈20(L1 α) / H1-1-1≈19(L2 α+β)")
print("OK wrote lineage_labels_regions.tsv")
