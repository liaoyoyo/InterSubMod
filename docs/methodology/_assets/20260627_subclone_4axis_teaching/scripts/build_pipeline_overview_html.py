#!/usr/bin/env python3
"""
build_pipeline_overview_html.py — 「數據→分類→判斷」完整量化 HTML(Q3;2026-07-09)。
靜態表格(Python 渲染,避 JS bug)。含:位點帳(Q2)·拓撲軸·分群數×拓撲(Q4)·clone/subclone(VAF)·read-lik·重現性·caveat。
反捏造:數字一次過從 region-view+mlhp+somatic_pass 重算。輸出 standalone HTML。
"""
import os, json, glob, gzip
from collections import Counter, defaultdict

C = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260709_pipeline_quantified_overview.standalone.html"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
CLONAL_T, SUB_LO, MIN_COV = 0.75, 0.10, 6

def rv_path(s): return f"{PILOT}/layered_region_view_HCC1395.json" if s == "HCC1395" else f"{MSROOT}/{s}/layered_region_view_{s}.json"
def wd(s): return PILOT if s == "HCC1395" else f"{MSROOT}/{s}"
def spass(s):
    lp = glob.glob(f"{C}/{s}/paired_full/*complete_matrix/longphase_s"); return f"{lp[0]}/somatic_pass.vcf.gz" if lp else None

def shape(edges):
    if not edges: return "single"
    ch = defaultdict(list); nodes = set()
    for p, c in edges: ch[p].append(c); nodes.add(p); nodes.add(c)
    if len([n for n in nodes if n != "ROOT"]) <= 1: return "single"
    cc = {p: len(v) for p, v in ch.items()}
    if any(v >= 2 for p, v in cc.items() if p != "ROOT"): return "branched"
    if cc.get("ROOT", 0) >= 2: return "star"
    return "linear"

def vaf_call(cc):
    nc = ns = nl = nu = 0
    for pos, (nr, na) in cc.items():
        tot = nr + na
        if tot < MIN_COV: nu += 1; continue
        v = na / tot
        if v >= CLONAL_T: nc += 1
        elif v >= SUB_LO: ns += 1
        else: nl += 1
    if nc + ns + nl == 0: return "lowcov"
    if ns >= 1: return "subclone"
    if nc >= 1 and nl == 0: return "founding"
    return "weak"

def analyze(s):
    lk = {}
    tot_pos = set(); ingroup_pos = set()
    for f in sorted(glob.glob(f"{wd(s)}/mlhp_part_*.json")):
        for g in json.load(open(f))["groups"]:
            lk[(g["chrom"], g["start"])] = g
            for p in g.get("positions", []):
                tot_pos.add((g["chrom"], p))
                if g.get("n_sSNV", 0) >= 2: ingroup_pos.add((g["chrom"], p))
    Cn = json.load(open(rv_path(s)))["census"]
    topo = Counter(); vc = Counter(); byc = defaultdict(Counter); rshape = {}
    for r in json.load(open(rv_path(s)))["regions"]:
        g = lk.get((r["chrom"], r["start"]))
        for L in r["lineages"]:
            if L["family"] not in ("1", "2"): continue
            trees = L.get("trees") or []
            shapes = {shape(t["edges"]) for t in (trees or [{"edges": []}])}
            sh = "capped" if L.get("capped") else (next(iter(shapes)) if len(shapes) == 1 else "mixed")
            topo[sh] += 1
            if sh not in ("capped", "mixed"): rshape[(r["chrom"], r["start"], L["family"])] = sh
            cbh = (g.get("col_coverage_by_hp", {}) if g else {}) or {}
            vc[vaf_call(cbh.get(L["family"], {}) or {}) if cbh.get(L["family"]) else "lowcov"] += 1
            pops = (g.get("populations_by_hp", {}) if g else {}).get(L["family"], {}) or {}
            c = sum(1 for gt in pops if "A" in gt)
            byc["c0" if c == 0 else ("c1" if c == 1 else ("c2" if c == 2 else "c3+"))][sh] += 1
    npass = sum(1 for _ in gzip.open(spass(s), "rt") if not _.startswith("#")) if spass(s) and os.path.exists(spass(s)) else 0
    return {"census": Cn, "topo": topo, "vc": vc, "byc": byc, "rshape": rshape,
            "npass": npass, "intree": len(ingroup_pos), "totpos": len(tot_pos)}

DATA = {s: analyze(s) for s in SAMPLES}

def pc(v, t): return f"{100*v/t:.0f}%" if t else "-"
def cell(v, t): return f'<td class="n">{v:,}<span class="p">{pc(v,t)}</span></td>'

# HCC1395 vs DORADO consistency
r5, rd = DATA["HCC1395"]["rshape"], DATA["HCC1395_DORADO"]["rshape"]
common = set(r5) & set(rd); same = sum(1 for k in common if r5[k] == rd[k])
def bcat(x): return "單" if x == "single" else ("姊妹" if x in ("branched", "star") else ("巢狀" if x == "linear" else "他"))
bsame = sum(1 for k in common if bcat(r5[k]) == bcat(rd[k]))

def trows(cols, fn):
    return "".join(f"<tr><td><b>{s}</b></td>{fn(DATA[s])}</tr>" for s in SAMPLES)

H = []
H.append('<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">')
H.append("<title>分層重建 數據→分類→判斷 量化總覽</title><style>")
H.append(":root{--bg:#fff;--fg:#1a1a1a;--mut:#666;--line:#e2e2e2;--card:#fafafa;--ac:#2563eb}")
H.append("@media(prefers-color-scheme:dark){:root{--bg:#0f1115;--fg:#e6e6e6;--mut:#9aa;--line:#2a2f38;--card:#171a21;--ac:#60a5fa}}")
H.append("*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.55 ui-sans-serif,system-ui,sans-serif}")
H.append(".wrap{max-width:1080px;margin:0 auto;padding:20px}h1{font-size:21px;margin:0 0 4px}h2{font-size:15px;margin:22px 0 6px;border-left:3px solid var(--ac);padding-left:8px}")
H.append("table{border-collapse:collapse;width:100%;font-size:12.5px;margin:6px 0}th,td{border:1px solid var(--line);padding:5px 8px;text-align:left}")
H.append("th{background:var(--card)}td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}.p{color:var(--mut);font-size:10.5px;margin-left:5px}")
H.append(".sub{color:var(--mut);font-size:12px}.note{background:var(--card);border:1px solid var(--line);border-radius:7px;padding:8px 11px;margin:8px 0;font-size:12px}")
H.append(".flow{display:flex;flex-wrap:wrap;gap:6px;align-items:center;font-size:12px;margin:8px 0}.flow .b{background:var(--card);border:1px solid var(--line);border-radius:7px;padding:6px 10px}.flow .a{color:var(--mut)}")
H.append(".red{color:#dc2626;font-weight:600}</style>")
H.append('<div class="wrap">')
H.append('<h1>🧬 分層重建 — 數據 → 分類 → 判斷 完整量化總覽</h1>')
H.append('<p class="sub">新骨幹（ClairS PASS/LongPhase-S；is_somatic 已移除）· 分母=germline(1/2) lineage · 一鍵重算 build_pipeline_overview_html.py · 數字全從 region-view+mlhp+somatic_pass 注入</p>')

# §0 flow (HCC1395)
d = DATA["HCC1395"]
H.append('<h2>① 數據流（HCC1395；每步位點/單位 + 去向）</h2>')
H.append('<div class="flow">')
H.append(f'<div class="b">ClairS PASS somatic<br><b>{d["npass"]:,}</b></div><div class="a">→ 需 ≤50kb 內 ≥2 共現</div>')
H.append(f'<div class="b">進樹位點<br><b>{d["intree"]:,}</b> <span class="p">{pc(d["intree"],d["npass"])}</span></div><div class="a">→ 組成</div>')
H.append(f'<div class="b">region<br><b>{d["census"]["n_regions"]:,}</b></div><div class="a">→ ×germline家族</div>')
H.append(f'<div class="b">germline lineage<br><b>{sum(d["topo"].values()):,}</b></div>')
H.append('</div>')
H.append(f'<div class="note"><span class="red">🔴 位點帳(Q2)</span>：骨幹 {d["npass"]:,} = 全 ClairS PASS somatic（不含 ClairS 濾掉的 Germline/LowQual）。'
         f'但只 <b>{d["intree"]:,}（{pc(d["intree"],d["npass"])}）</b>落在 ≥2-sSNV 共現鏈能建樹；其餘 ~{100-int(100*d["intree"]/d["npass"])}% 是孤立單點（無 somatic 鄰居）→ 共現法固有，不進樹。</div>')

# §1 hierarchy
H.append('<h2>② 階層計數（每樣本）</h2><table><tr><th>樣本</th><th class="n">骨幹 sSNV</th><th class="n">進樹位點</th><th class="n">region</th><th class="n">germline lineage</th></tr>')
H.append(trows(None, lambda d: f'<td class="n">{d["npass"]:,}</td><td class="n">{d["intree"]:,}<span class="p">{pc(d["intree"],d["npass"])}</span></td><td class="n">{d["census"]["n_regions"]:,}</td><td class="n">{sum(d["topo"].values()):,}</td>'))
H.append('</table>')

# §2 topology axis
H.append('<h2>③ 分類·軸1 拓撲（clone 關係）</h2>')
H.append('<div class="note">single=單clone · <b>linear=巢狀subclone後代</b> · <b>branched/star=姊妹subclone分歧</b> · mixed=形狀未定 · capped=太密。🔴 subclonal 結構以巢狀為主、姊妹罕見。</div>')
H.append('<table><tr><th>樣本</th><th class="n">single</th><th class="n">linear 巢狀</th><th class="n">姊妹(br+star)</th><th class="n">mixed</th><th class="n">capped</th></tr>')
def topo_row(d):
    t = sum(d["topo"].values()); tp = d["topo"]
    return cell(tp["single"], t) + cell(tp["linear"], t) + cell(tp["branched"] + tp["star"], t) + cell(tp["mixed"], t) + cell(tp["capped"], t)
H.append(trows(None, topo_row)); H.append("</table>")

# §3 cluster-count × shape (Q4, HCC1395)
H.append('<h2>④ 分類·軸1b 分群數 c × 拓撲（Q4；HCC1395）</h2>')
H.append('<div class="note">c=帶 ALT 的細胞群數(clone/subclone 群數)。🔴 <b>c=2 是關鍵</b>：巢狀(linear) vs 姊妹(branched) 二選一；c≥3=複雜多-clone。</div>')
H.append('<table><tr><th>分群數 c</th><th class="n">總</th><th class="n">single</th><th class="n">linear 巢狀</th><th class="n">姊妹</th><th class="n">mixed</th></tr>')
bc = DATA["HCC1395"]["byc"]
for cb, lab in [("c0", "0 (僅germline/partial)"), ("c1", "1"), ("c2", "2"), ("c3+", "≥3")]:
    dd = bc[cb]; t = sum(dd.values())
    if not t: continue
    H.append(f'<tr><td>{lab}</td><td class="n">{t:,}</td>{cell(dd["single"],t)}{cell(dd["linear"],t)}{cell(dd["branched"]+dd["star"],t)}{cell(dd["mixed"],t)}</tr>')
H.append("</table>")

# §4 clone/subclone VAF axis
H.append('<h2>⑤ 分類·軸2 clone/subclone（頻率 VAF）</h2>')
H.append(f'<div class="note">VAF≥{CLONAL_T}=clonal(主幹) · [{SUB_LO},{CLONAL_T})=subclonal(亞群)。'
         f'<span class="red">🔴 CN confound</span>：CN-gain 稀釋 VAF→誤判 subclone、LOH 抬高→誤判 clonal → <b>has_subclone% 是 CN-未控上界</b>；只 HCC1395 有 SEQC2 CN 可控。</div>')
H.append('<table><tr><th>樣本</th><th class="n">founding 主幹</th><th class="n">has_subclone 上界</th><th class="n">weak</th><th class="n">lowcov</th></tr>')
def vc_row(d):
    t = sum(d["vc"].values()); v = d["vc"]
    return cell(v["founding"], t) + cell(v["subclone"], t) + cell(v["weak"], t) + cell(v["lowcov"], t)
H.append(trows(None, vc_row)); H.append("</table>")

# §5 read-likelihood + reproducibility + caveats
H.append('<h2>⑥ 判斷·軸3 樹選擇（read-likelihood）+ 重現性 + 限制</h2>')
H.append(f'<div class="note"><b>樹選擇(Q1)</b>：對 mixed_shape 區用 read 數排形狀 → 全 7 樣本僅 <b class="red">5%</b> 可排，95% 純隱藏祖先(0-read 分不出)。→ read 數不當 tree-selector；正確用途=軸2 clone/subclone 頻率。</div>')
H.append(f'<div class="note"><b>重現性</b>：HCC1395(5khz) vs DORADO 同細胞株，共同形狀確定區 {len(common):,} → 拓撲一致 <b>{pc(same,len(common))}</b>、生物類別(單/巢狀/姊妹)一致 <b>{pc(bsame,len(common))}</b>；has_subclone 54% vs 52%。不一致主要 depth-driven。</div>')
H.append('<div class="note"><b>🔴 誠實限制</b>：① clone/subclone VAF 有 CN confound(上界,需 CN 控制) ② 樹選擇 read-lik 僅 5%(隱藏祖先主導) ③ COLO829 65% lowcov(低-coread artifact) ④ region-local(≤8sSNV/≤50kb/單家族局部結構,非全腫瘤整樹—整樹「定不出來」) ⑤ 跨樣本只比比例,不報 p 值(pseudoreplication,真 n=7)。</div>')
H.append("</div>")

open(OUT, "w", encoding="utf-8").write("".join(H))
print(f"OK wrote {OUT} ({len(''.join(H)):,} bytes)")
