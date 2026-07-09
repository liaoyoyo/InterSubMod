#!/usr/bin/env python3
"""
analyze_clone_subclone_vaf.py — clone/subclone 判別(read 佔比/VAF;2026-07-09 使用者定案 read 數的正確用途)。
每 germline lineage:用 mlhp col_coverage_by_hp 的 per-position [nREF,nALT] 算 VAF(=該突變在家族細胞佔比)。
  clonal(VAF≥CLONAL_T)=突變在~全部家族細胞=主幹/founding · subclonal(SUB_LO≤VAF<CLONAL_T)=一部分細胞=subclone
  · low(<SUB_LO)=弱證據。lineage 判:全 clonal→founding_clonal;有 subclonal→has_subclone。結合拓撲(region-view)。
🔴 這是「頻率軸」(clone vs subclone),與拓撲(關係軸)、分支突變數(穩健軸)三軸互補。
用法: python3 analyze_clone_subclone_vaf.py
"""
import os, sys, json, glob
from collections import Counter, defaultdict

MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
CLONAL_T = 0.75    # VAF≥此 → clonal(家族內~全細胞;考慮 ONT 噪聲不設 1.0)
SUB_LO = 0.10      # VAF∈[此,CLONAL_T) → subclonal;<此 → 弱/雜訊
MIN_COV = 6        # 位點 total 覆蓋 ≥此才判(否則 power 不足→undetermined)

def rv_path(s):
    return f"{PILOT}/layered_region_view_HCC1395.json" if s == "HCC1395" else f"{MSROOT}/{s}/layered_region_view_{s}.json"

def mlhp_lookup(s):
    wd = PILOT if s == "HCC1395" else f"{MSROOT}/{s}"
    lk = {}
    for f in sorted(glob.glob(f"{wd}/mlhp_part_*.json")):
        for g in json.load(open(f))["groups"]:
            lk[(g["chrom"], g["start"])] = g
    return lk

def tree_shape(edges):
    if not edges:
        return "single"
    ch = defaultdict(list); nodes = set()
    for p, c in edges:
        ch[p].append(c); nodes.add(p); nodes.add(c)
    if len([n for n in nodes if n != "ROOT"]) <= 1:
        return "single"
    cc = {p: len(cs) for p, cs in ch.items()}
    if any(v >= 2 for p, v in cc.items() if p != "ROOT"):
        return "branched_internal"
    if cc.get("ROOT", 0) >= 2:
        return "star_root"
    return "linear"

def classify_lineage(colcov):
    """回 (call, n_clonal, n_sub, n_low, n_undet)。colcov={pos:[nREF,nALT]}。"""
    nc = ns = nl = nu = 0
    for pos, (nr, na) in colcov.items():
        tot = nr + na
        if tot < MIN_COV:
            nu += 1; continue
        vaf = na / tot
        if vaf >= CLONAL_T:
            nc += 1
        elif vaf >= SUB_LO:
            ns += 1
        else:
            nl += 1
    # lineage call(基於有 power 的位點)
    powered = nc + ns + nl
    if powered == 0:
        return "undetermined_lowcov", nc, ns, nl, nu
    if ns >= 1:
        return "has_subclone", nc, ns, nl, nu        # ≥1 突變在一部分細胞 = subclone 出現
    if nc >= 1 and nl == 0:
        return "founding_clonal", nc, ns, nl, nu      # 全 clonal = 主幹/founding(無 subclonal 分裂)
    return "weak_only", nc, ns, nl, nu                # 只有低 VAF 弱證據

def run(s):
    lk = mlhp_lookup(s)
    regions = json.load(open(rv_path(s)))["regions"]
    call_c = Counter(); combo = Counter(); vafdist = Counter()
    for r in regions:
        g = lk.get((r["chrom"], r["start"]))
        if not g:
            continue
        cbh = g.get("col_coverage_by_hp", {}) or {}
        for L in r["lineages"]:
            fam = L["family"]
            if fam not in ("1", "2"):
                continue
            cc = cbh.get(fam, {}) or {}
            cc = {p: v for p, v in cc.items()}
            if not cc:
                continue
            call, nclo, nsub, nlow, nu = classify_lineage(cc)
            call_c[call] += 1
            # 結合拓撲(subclone 的關係型態)
            trees = L.get("trees") or []
            shapes = {tree_shape(t["edges"]) for t in trees}
            shape = next(iter(shapes)) if len(shapes) == 1 else ("mixed" if shapes else "none")
            if call == "has_subclone":
                bio = "巢狀subclone" if shape == "linear" else ("姊妹subclone" if shape in ("branched_internal", "star_root") else ("subclone(形狀未定)" if shape == "mixed" else "subclone"))
                combo[bio] += 1
    return call_c, combo

print("=" * 92)
print("clone/subclone 判別(VAF 頻率軸): VAF≥%.2f=clonal(主幹) / [%.2f,%.2f)=subclonal(亞群) / <%.2f=弱" % (CLONAL_T, SUB_LO, CLONAL_T, SUB_LO))
print("=" * 92)
print(f"{'樣本':15}{'lineage':>8}{'founding主幹':>13}{'%':>6}{'has_subclone':>13}{'%':>6}{'弱/低cov':>10}")
print("-" * 92)
tot = Counter(); tcombo = Counter()
for s in SAMPLES:
    cc, combo = run(s)
    for k, v in cc.items():
        tot[k] += v
    for k, v in combo.items():
        tcombo[k] += v
    n = sum(cc.values()); fc = cc["founding_clonal"]; hs = cc["has_subclone"]
    weak = cc["weak_only"] + cc["undetermined_lowcov"]
    print(f"{s:15}{n:>8}{fc:>13}{f'{100*fc/n:.0f}%' if n else '-':>6}{hs:>13}{f'{100*hs/n:.0f}%' if n else '-':>6}{weak:>10}")
    if s == "HCC1395":
        print(f"    HCC1395 subclone 關係型態: {dict(combo)}")
print("-" * 92)
n = sum(tot.values()); fc = tot["founding_clonal"]; hs = tot["has_subclone"]
print(f"{'全 7 樣本':15}{n:>8}{fc:>13}{f'{100*fc/n:.0f}%':>6}{hs:>13}{f'{100*hs/n:.0f}%':>6}{tot['weak_only']+tot['undetermined_lowcov']:>10}")
print(f"\n全樣本 has_subclone 關係型態(拓撲): {dict(tcombo)}")
print(f"裁決:VAF 頻率軸 → founding主幹 {100*fc/n:.0f}% / 有subclone {100*hs/n:.0f}%;subclone 中巢狀 vs 姊妹見上。")
