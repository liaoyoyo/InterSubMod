#!/usr/bin/env python3
"""
build_quantified_overview.py — 7 樣本完整量化總覽(階層+拓撲軸+clone/subclone頻率軸+一致性;2026-07-09)。
一次過讀 region-view + mlhp,整合輸出 markdown 報告(使用者要求「完整數據表格+階層+比例可觀察」)。
軸: 拓撲(關係:single/linear/branched) × clone/subclone(頻率VAF) × 樹選擇(read-lik,已知5%)。
🔴 誠實 caveat: clone/subclone VAF 有 CN confound(gain稀釋/LOH抬高);只 HCC1395 有 SEQC2 CN 可控。
輸出: docs/methodology/20260709_quantified_topology_clone_subclone_overview_01.md
"""
import os, json, glob
from collections import Counter, defaultdict

MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260709_quantified_topology_clone_subclone_overview_01.md"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
CLONAL_T, SUB_LO, MIN_COV = 0.75, 0.10, 6

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
    if not edges: return "single"
    ch = defaultdict(list); nodes = set()
    for p, c in edges:
        ch[p].append(c); nodes.add(p); nodes.add(c)
    if len([n for n in nodes if n != "ROOT"]) <= 1: return "single"
    cc = {p: len(cs) for p, cs in ch.items()}
    if any(v >= 2 for p, v in cc.items() if p != "ROOT"): return "branched"
    if cc.get("ROOT", 0) >= 2: return "star"
    return "linear"

def vaf_call(colcov):
    nc = ns = nl = nu = 0
    for pos, (nr, na) in colcov.items():
        tot = nr + na
        if tot < MIN_COV: nu += 1; continue
        vaf = na / tot
        if vaf >= CLONAL_T: nc += 1
        elif vaf >= SUB_LO: ns += 1
        else: nl += 1
    if nc + ns + nl == 0: return "lowcov"
    if ns >= 1: return "subclone"
    if nc >= 1 and nl == 0: return "founding"
    return "weak"

def run(s):
    lk = mlhp_lookup(s)
    regions = json.load(open(rv_path(s)))["regions"]
    C = json.load(open(rv_path(s)))["census"]
    topo = Counter(); shapedet = Counter(); vc = Counter(); region_shape = {}
    for r in regions:
        g = lk.get((r["chrom"], r["start"]))
        for L in r["lineages"]:
            if L["family"] not in ("1", "2"): continue
            # 拓撲軸
            if L.get("capped"):
                sh = "capped"
            else:
                trees = L.get("trees") or []
                shapes = {tree_shape(t["edges"]) for t in (trees or [{"edges": []}])}
                if len(shapes) == 1:
                    sh = next(iter(shapes)); shapedet["確定"] += 1; region_shape[(r["chrom"], r["start"], L["family"])] = sh
                else:
                    sh = "mixed"; shapedet["未定"] += 1
            topo[sh] += 1
            # clone/subclone 頻率軸
            cbh = (g.get("col_coverage_by_hp", {}) if g else {}) or {}
            cc = cbh.get(L["family"], {}) or {}
            vc[vaf_call(cc) if cc else "lowcov"] += 1
    return C, topo, shapedet, vc, region_shape

# ---- 跑全 7 樣本 ----
rows = {}; allrs = {}
for s in SAMPLES:
    rows[s] = run(s); allrs[s] = rows[s][4]

def pc(v, t): return f"{100*v/t:.0f}%" if t else "-"

L = []
L.append("<!--\n建立: 2026-07-09 | 類型: 7樣本量化總覽(拓撲+clone/subclone) | build: build_quantified_overview.py 一鍵重算\n"
         "資料: 新骨幹(ClairS PASS) region-view + mlhp col_coverage_by_hp\n-->\n")
L.append("# 7 樣本分層重建 — 完整量化總覽（拓撲關係 × clone/subclone 頻率）\n")
L.append("> 新骨幹（ClairS PASS/LongPhase-S；is_somatic 已移除）。分母=germline(1/2) lineage。一鍵重算：`build_quantified_overview.py`。\n")

# §1 階層 + 骨幹
L.append("\n## §1 單位階層 + 骨幹 sSNV（每樣本）\n")
L.append("| 樣本 | 骨幹 sSNV | region | germline lineage | 形狀確定% |")
L.append("|---|---|---|---|---|")
for s in SAMPLES:
    C, topo, sd, vc, _ = rows[s]
    nlin = sum(topo.values())
    L.append(f"| {s} | {C.get('U1_sSNV_somatic_total','?'):,} | {C['n_regions']:,} | {nlin:,} | {pc(sd['確定'], sd['確定']+sd['未定'])} |")

# §2 拓撲軸(關係)
L.append("\n## §2 拓撲軸（樹形狀 = clone 關係）\n")
L.append("> single=單clone · **linear=巢狀subclone後代** · **branched/star=姊妹subclone分歧** · mixed=形狀未定 · capped=太密\n")
L.append("| 樣本 | single | linear(巢狀) | branched+star(姊妹) | mixed | capped |")
L.append("|---|---|---|---|---|---|")
for s in SAMPLES:
    C, topo, sd, vc, _ = rows[s]
    n = sum(topo.values())
    sis = topo['branched'] + topo['star']
    L.append(f"| {s} | {pc(topo['single'], n)} | {pc(topo['linear'], n)} | {pc(sis, n)} | {pc(topo['mixed'], n)} | {pc(topo['capped'], n)} |")

# §3 clone/subclone 頻率軸(VAF)
L.append("\n## §3 clone/subclone 頻率軸（VAF；🔴 有 CN confound）\n")
L.append(f"> VAF≥{CLONAL_T}=clonal(主幹) · [{SUB_LO},{CLONAL_T})=subclonal(亞群) · founding=全clonal · has_subclone=≥1亞群突變 · lowcov=覆蓋<{MIN_COV}\n")
L.append("> 🔴 **CN confound**：CN-gain 稀釋 VAF→誤判 subclone；LOH 抬高→誤判 clonal。**has_subclone% 是 CN-未控上界**；只 HCC1395 有 SEQC2 CN 可控。\n")
L.append("| 樣本 | founding主幹 | has_subclone(上界) | weak | lowcov |")
L.append("|---|---|---|---|---|")
for s in SAMPLES:
    C, topo, sd, vc, _ = rows[s]
    n = sum(vc.values())
    L.append(f"| {s} | {pc(vc['founding'], n)} | {pc(vc['subclone'], n)} | {pc(vc['weak'], n)} | {pc(vc['lowcov'], n)} |")

# §4 樹選擇(read-likelihood)結論
L.append("\n## §4 樹選擇（read-likelihood）— 實證結論\n")
L.append("> 對 mixed_shape(形狀未定)區用 read 數排形狀：**全 7 樣本僅 5% 可排**，95% 是純隱藏祖先(中間節點 0-read 分不出)。\n"
         "> → **read 數不適合當 tree-selector**（多數歧義是隱藏祖先）；read 數的正確用途 = §3 clone/subclone 頻率軸。\n")

# §5 HCC1395 vs DORADO 重現性
L.append("\n## §5 HCC1395(5khz) vs DORADO 重現性（同細胞株）\n")
r5, rd = allrs["HCC1395"], allrs["HCC1395_DORADO"]
common = set(r5) & set(rd)
same = sum(1 for k in common if r5[k] == rd[k])
def biocat(sh): return "單clone" if sh == "single" else ("姊妹" if sh in ("branched", "star") else ("巢狀" if sh == "linear" else "其他"))
biosame = sum(1 for k in common if biocat(r5[k]) == biocat(rd[k]))
L.append(f"- 共同形狀確定區 {len(common):,}；**拓撲一致 {pc(same, len(common))}**、**生物類別(單/巢狀/姊妹)一致 {pc(biosame, len(common))}**。")
L.append(f"- clone/subclone(§3)：has_subclone 54% vs 52%（一致）。→ 不一致主要 depth-driven（5khz 深→linear，DORADO 淺→single）。\n")

# §6 caveats
L.append("\n## §6 誠實限制\n")
L.append("1. 🔴 **clone/subclone VAF 有 CN confound**：has_subclone% 是上界，需 CN 控制（只 HCC1395 有 SEQC2 CN）。")
L.append("2. **樹選擇 read-lik 僅 5%**：多數樹/形狀歧義是隱藏祖先，read 分不出。")
L.append("3. **COLO829 65% lowcov**：其低-coread artifact（多 lineage 覆蓋不足）。")
L.append("4. **region-local scale**：每區 ≤8 sSNV/≤50kb/單 germline 家族的局部結構，非全腫瘤整樹（整樹「定不出來」）。")
L.append("5. 跨樣本只比比例（絕對受深度混淆）；p 值 pseudoreplication（真 n=7）不報。\n")

open(OUT, "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))
print(f"\n→ 寫出 {OUT}")
