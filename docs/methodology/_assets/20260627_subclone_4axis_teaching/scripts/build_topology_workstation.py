#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互動式拓樸工作站 HTML v2：
  - 每區 genotype 向量(S1/S2..標籤) + 克隆樹 SVG(節點=S-mut-set·reads·%·HP) + 分支(V/H)
  - 篩選: chr / genome_ctx(telo/centro/arm) / 拓樸型 / determinacy / 含FP / 最少群數 / 搜尋
  - 排序: 座標 / 複雜度(n_sSNV) / 群數 / region
  - region badge: genome_ctx + TP/FP 組成
  - chr17 完整 worked panel(S/r/m 一致標籤 + 樹 + 16 sig CpG)
§13-A: 數字由 topology_per_region.json 注入。
用法: python3 build_topology_workstation.py
"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
MULTI_OUT = os.path.normpath(os.path.join(HERE, "..", "..", "20260629_multisample_topology_workstation.standalone.html"))
# OUT(舊單樣本 20260628)已 deprecated:build 只產多樣本 MULTI_OUT(=主結果);舊單樣本檔不再寫出

# 多樣本(2026-06-29):SM_SAMPLES="name:dir,name:dir" → 多分頁;預設納入已完成樣本(HCC1395 凍結 + multisample_subclone 下有 topology 的)。
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
def _sample_dirs():
    env = os.environ.get("SM_SAMPLES", "")
    if env:
        return [tuple(x.split(":", 1)) for x in env.split(",") if ":" in x]
    pairs = [("HCC1395", DATA)]  # 凍結主樣本
    if os.path.isdir(MSROOT):
        for s in sorted(os.listdir(MSROOT)):
            tp = os.path.join(MSROOT, s, "topology_per_region.json")
            if os.path.exists(tp):
                pairs.append((s, os.path.join(MSROOT, s)))
    return pairs

def _load_sample(dr):
    d = json.load(open(os.path.join(dr, "topology_per_region.json"), encoding="utf-8"))
    rec = {"stats": d["stats"], "detail": d["detail"], "chr17": d.get("chr17_worked"),
           "chroms": d.get("chroms", []), "provenance": d.get("provenance", {})}
    csp = os.path.join(dr, "candidate_scoring.json")
    if os.path.exists(csp):
        cs = json.load(open(csp, encoding="utf-8"))
        rec["scoring"] = {"summary": {k: cs.get(k) for k in ("n_total","n_need_confirm","score_formula","situation_dist","resolution_dist","score_buckets","needs_methyl_n")},
                          "queue": cs.get("queue", [])}
    else:
        rec["scoring"] = {"summary": {}, "queue": []}
    gp = os.path.join(dr, "region_gene_annotation.json")
    rec["gene"] = json.load(open(gp, encoding="utf-8")).get("regions", {}) if os.path.exists(gp) else {}
    accp = os.path.join(dr, "single_snv_accounting.json")
    rec["accounting"] = json.load(open(accp, encoding="utf-8")) if os.path.exists(accp) else None
    ctp = os.path.join(dr, "candidate_trees.json")  # R6/Part B: enumerate_candidate_trees 誠實版
    rec["candtrees"] = {x["region"]: x for x in json.load(open(ctp, encoding="utf-8")).get("candidate_trees", [])} if os.path.exists(ctp) else {}
    return rec

SAMPLES = {name: _load_sample(dr) for name, dr in _sample_dirs()}
SAMPLE_NAMES = list(SAMPLES.keys())
SAMPLES_JSON = json.dumps(SAMPLES, ensure_ascii=False)
# R3+R7: chr17 read×read 距離矩陣 + 分群×lineage 交叉表(HCC1395 worked;固定教學甲基展板用)
_c17t = os.path.join(DATA, "chr17_tree_data.json")
CHR17TREE_JSON = json.dumps(json.load(open(_c17t, encoding="utf-8")), ensure_ascii=False) if os.path.exists(_c17t) else "null"
FIRST = SAMPLE_NAMES[0]
# 相容:DJ/B 供舊單樣本變數(UNIVERSE_BANNER 等)用 FIRST 樣本
d = json.load(open(os.path.join(_sample_dirs()[0][1], "topology_per_region.json"), encoding="utf-8"))
acc = SAMPLES[FIRST].get("accounting")
B = acc["buckets"] if acc else {"linked": {"n": 0, "pct": 0}, "underpowered": {"n": 0, "pct": 0, "ccf_tendency": {}}, "isolated": {"n": 0, "pct": 0}}
DJ = SAMPLES_JSON
UNIVERSE_BANNER = f"""<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:11px 14px;margin:10px 0;font-size:12.5px">
<b>全 sSNV 宇宙帳本（{acc['universe_total']:,} = TP {acc['tp']:,} + FP {acc['fp']:,}；無遺漏 sum-check ✓）</b><br>
<div style="display:flex;gap:18px;flex-wrap:wrap;margin-top:6px">
<span>🟢 <b>linked {B['linked']['n']:,}（{B['linked']['pct']}%）</b> 可建樹→拓樸分析（下方）</span>
<span>🟡 <b>underpowered {B['underpowered']['n']:,}（{B['underpowered']['pct']}%）</b> 有 partner 無共讀→<b>加深覆蓋可救</b>；有 CCF（clonal {B['underpowered']['ccf_tendency'].get('high_ge0.4(clonal)',0)}/mid {B['underpowered']['ccf_tendency'].get('mid_0.1-0.4',0)}/low {B['underpowered']['ccf_tendency'].get('low_lt0.1(subclonal/noise)',0)}）</span>
<span>⚪ <b>isolated {B['isolated']['n']:,}（{B['isolated']['pct']}%）</b> read-span 內無 partner（Tier-R 樹外）；<b>有 caller VAF 可刻畫 + 可能 Tier-PS partner</b></span>
</div>
<div class="note" style="margin-top:6px">🔑 單位點<b>非「全無法處理」</b>：underpowered 有 CCF+可深覆蓋救；isolated 有 caller VAF 可放 clonal 譜 + 可能 same-PS(&gt;50kb) 連鎖(Tier-PS 未做)。真正拓樸-dead = isolated 中無 same-PS partner 者（待 Tier-PS 量化）。</div></div>"""

GLOSSARY = [
 ("sSNV / S1·S2·S3", "體細胞單核苷酸變異；S1..Sk = 區內依座標排序的 sSNV（基因型向量第 i 位 = Si）。", "癌細胞才有、正常細胞沒有的點突變。一個區域有 k 個就標 S1..Sk，順序按基因座位置。"),
 ("read 群 / r / population", "同一基因型向量(如 RAR)的 read 集合 = 一種『細胞狀態』。", "一條 read=一條 DNA 分子=一個細胞的一條染色體。攜帶相同突變組合的 read 歸為同一群（population），代表一個 lineage 節點。"),
 ("m / 甲基位點", "顯著差異的 CpG 甲基化位點（chr17: m1..m16）。", "DNA 甲基化標記；m 是 L1 vs L2 lineage 間 |Δβ| 大的 CpG。⚠ 實測對齊 cis-genotype 軸，非獨立 lineage。"),
 ("HP / H1·H2·H3", "germline 單倍型（哪條親代染色體）。H1/H2=兩條;H3?=未定相(somatic-ALT read)。", "由 longphase-S haplotag 決定。正常人 2 條 haplotype 根 H1、H2。HP tag『1-1』→H1、『2-1』→H2、『3』→H3?。"),
 ("HP{h}-path", "lineage 標籤 = HP{根}-{分支1}-{分支2}…（Dewey 路徑）。", "如 H1-1=H1 上第一個 somatic 事件;H1-1-1=其後代;H1-2=姊妹分支。分支編號=VAF 遞減。"),
 ("vertical 直系 / nested", "祖先→後代（一個細胞先後累積兩突變）。", "2×2 有 AA 格(兩突變同 read) + 一側零格 → 巢狀。樹上往下一層。"),
 ("horizontal 姊妹 / sibling", "兩突變從不共現但同 haplotype → 不同 subclone 平行分支。", "2×2 的 AA=0(兩 ALT 從不同 read) + same-HP。樹上同層分叉。"),
 ("co_linked", "兩突變完美共現(只見 AA) → 同一 lineage 事件(同節點)。", "RA=AR=0、AA≥2。兩突變永遠同進退，無法內部排序。"),
 ("mutual_excl 互斥", "兩 ALT 從不共現(AA=0)。", "diff-HP→allelic(兩條染色體各自突變,非 subclone);same-HP→sibling subclone。"),
 ("independent / four-gamete 違反", "RR/RA/AR/AA 四格全有 → 不相容單一樹。", "違反無限位點假設(回復/重複突變/CNV multiplicity/偽影)。"),
 ("perfect-phylogeny 完美系統發生樹", "二元字元(REF/ALT)的樹:每位點只突變一次。", "古典定理:二元字元『每對相容』即足以保證整棵樹存在 → pairwise 拼接合法，不需單分子整跨。"),
 ("2×2 共現 (RR/RA/AR/AA)", "對每對 sSNV 數共讀 read 在兩位點的 REF/ALT 組合。", "RR=都不帶、AA=都帶、RA/AR=只一帶。哪格為零決定關係(co_linked/nested/互斥)。"),
 ("ε=2% 噪聲底線", "cell 為真 ⟺ count > coread×2%（ONT 錯誤率）。", "保留最低 1 條(低 coread);高 coread 單讀(1≤coread×2%)判噪聲。經 FP 裁判+結構穩定+塌陷集中三路定案。"),
 ("coread 共讀", "同時覆蓋兩個位點的 read 數。", "≥6 才算 powered。決定一個零格是否可信。"),
 ("VAF / CCF", "VAF=變異等位基因頻率;CCF=癌細胞比例(clonal prevalence)。", "高 VAF=clonal(早/大);低=subclonal(晚/小)。只在 CN-clean 可信。用於分支編號 + 單位點刻畫。"),
 ("determinacy", "拓樸能否唯一辨識:A_determined(單分子向量) / A_ambiguous(缺中間群順序未定) / B_pairwise(拼接非單分子) / C_underdetermined(多樹相容) / incompatible(成環)。", "『樹存在』(99.4%相容)≠『能辨識是哪棵』(僅~11%)。"),
 ("situation tier", "A 單分子整跨 / B 可整跨pairwise / C 必鏈接(span>read)。", "≥3 位點先分 situation 再處理:有沒有一條 read 穿過全部決定證據強度。"),
 ("genome_ctx", "telomere(端粒,≤3Mb 端)/ centromere(著絲點±3Mb)/ arm(染色體臂)。", "hg38 染色體長度+centromere 近似。centromere/telomere 區偽影風險高。"),
 ("TP / FP", "SEQC2 truth set 標的真/假陽性。", "🔴 只用於觀察評估,絕不進前處理/定義(build 用 TP∪FP union + normal 比對定 somatic)。"),
 ("Tier-R / Tier-PS", "Tier-R=same-read(≤50kb,同分子);Tier-PS=same phase-set(>50kb,統計相位,未做)。", "克隆連鎖只認 Tier-R;isolated 區可能有 Tier-PS partner 待救。"),
 ("cluster-count (c, k+1 上界)", "區內 distinct population 數;perfect-phylogeny 下 ≤ k+1(非 2^k)。", "實測 99.9% n_pop≤k+1、中位 2 → 拓樸搜尋空間極小。先定 c 再縮限拓樸。"),
 ("ambiguous-parentage 缺中間群", "節點突變集跳>1(中間 population 沒觀察到)→累積順序未定。", "76 區。如 {0,3,4} 缺 {0,3} 等中間群 → 0,3,4 哪個先未定。"),
 ("linked / underpowered / isolated", "全 sSNV 三桶:可建樹 / 有 partner 無共讀(可救) / 無 partner(Tier-R 樹外)。", "61% / 15.4% / 23.5%。單位點非全無法處理:underpowered 有 CCF、isolated 有 caller VAF+可能 Tier-PS。"),
 ("cis-ASM / double-dip", "甲基隨突變的 cis 局部效應 / 用同量定群又驗群的循環。", "chr17 證甲基分群對齊突變 genotype 軸(cis)非獨立 lineage → 甲基不能當獨立驗證器。06-28 normal cis-control 已測:CROSS-HP 35.4% 可控、SAME-HP 多數區 normal 無對應 within-HP 軸=結構性無法 control(需 single-cell)。"),
 ("bounded-auxiliary 甲基定位", "甲基=corroborate 非 detect 的有界輔助(Tier-3 機率層)。", "排序:genetic 共現>HP 定根>甲基。🔴 06-28 cis-control pilot 已測定案:cis-control 只對 CROSS-HP 區有效(35.4%)、SAME-HP 59% 結構性無解→甲基不能升 resolver;最終角色=cluster-count sanity + 43 區 CROSS-HP 弱排序 PoC。"),
 ("⭐3 / single-pipeline", "單樣本 HCC1395 單一 pipeline 的證據上限。", "升 ⭐4 需 ≥5/7 樣本+COLO829+single-cell 正交確認。"),
]
GLOSSARY_CHAPTERS = [("① 基本標籤", [0,1,2,3,4]), ("② SNV 關係 / 拓樸型", [5,6,7,8,9,10]), ("③ 共現與證據量", [11,12,13,14]), ("④ 可辨識性 determinacy", [15,16,20,21]), ("⑤ 基因體與真值", [17,18,19]), ("⑥ 全 sSNV 帳本與甲基", [22,23,24,25])]
def _gterm(i):
    t, s, dd = GLOSSARY[i]
    return f'<details style="border:1px solid #f1f3f5;border-radius:6px;padding:5px 9px;font-size:12px"><summary style="cursor:pointer"><b>{t}</b></summary><div style="margin-top:4px;color:#343a40">{s}</div><div style="margin-top:3px;color:#868e96;font-size:11px">{dd}</div></details>'
GLOSSARY_HTML = ('<details style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:10px 14px;margin:10px 0"><summary style="cursor:pointer;font-weight:600;color:#1971c2">📖 名詞與概念解釋（分 6 章節；點章節→點詞展開）</summary>'
 + "".join('<details style="border:1px solid #e9ecef;border-radius:6px;padding:6px 10px;margin-top:6px;background:#f8f9fa"><summary style="cursor:pointer;font-weight:600;color:#495057">' + ch + f'（{len(idxs)} 詞）</summary><div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(290px,1fr));gap:6px;margin-top:8px">' + "".join(_gterm(i) for i in idxs) + '</div></details>' for ch, idxs in GLOSSARY_CHAPTERS)
 + '</details>')

CSS = """
*{box-sizing:border-box}body{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:#212529;background:#f8f9fa}
.wrap{max-width:1320px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:.2em 0}h3{margin:.3em 0}.sub{color:#868e96;font-size:12.5px}
.stats{display:flex;gap:12px;flex-wrap:wrap;margin:10px 0}.scard{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:9px 11px;min-width:200px}
.scard h4{margin:0 0 5px;font-size:11.5px;color:#495057}.bar{display:flex;align-items:center;gap:5px;font-size:10.5px;margin:2px 0}.bar i{height:10px;background:#1c7ed6;border-radius:2px;display:inline-block}
.scard{cursor:pointer;transition:box-shadow .12s}.scard:hover{box-shadow:0 2px 12px rgba(0,0,0,.12)}.scard h4 .more{float:right;font-size:9px;color:#adb5bd;font-weight:400}
.smbg{position:fixed;inset:0;background:rgba(0,0,0,.45);z-index:50;display:none;align-items:center;justify-content:center;padding:20px}
.smbox{background:#fff;border-radius:12px;max-width:780px;width:100%;max-height:86vh;overflow:auto;padding:18px 22px;position:relative}
.smbox .x{position:absolute;top:10px;right:15px;cursor:pointer;font-size:22px;color:#868e96;background:none;border:none;line-height:1}
details.c17{background:#fff;border:1px solid #ffd8a8;border-radius:8px;padding:10px 14px;margin:10px 0}details.c17 summary{cursor:pointer;font-weight:600;color:#d9480f}
.ctrl{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:10px 0;font-size:12.5px;background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:9px}
.ctrl select,.ctrl input{padding:3px 7px;border:1px solid #ced4da;border-radius:5px;font-size:12.5px}
.main{display:grid;grid-template-columns:400px 1fr;gap:12px}@media(max-width:860px){.main{grid-template-columns:1fr}}
.list{background:#fff;border:1px solid #dee2e6;border-radius:8px;max-height:76vh;overflow:auto}
.row{padding:6px 10px;border-bottom:1px solid #f1f3f5;cursor:pointer;font-size:12px}.row:hover{background:#e7f5ff}.row.sel{background:#d0ebff}.row b{color:#1c7ed6}
.tag{font-size:9.5px;padding:1px 6px;border-radius:9px;margin-left:3px}
.t_linear{background:#d3f9d8;color:#2b8a3e}.t_branched{background:#e5dbff;color:#5f3dc4}.t_star{background:#fff3bf;color:#b08900}.t_single{background:#f1f3f5;color:#868e96}
.ctx_telomere{background:#d0ebff;color:#1971c2}.ctx_centromere{background:#ffe3e3;color:#c92a2a}.ctx_arm{background:#f1f3f5;color:#868e96}
.facets{display:flex;gap:10px;flex-wrap:wrap;margin:8px 0}
.fcard{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:8px 11px;flex:1;min-width:230px}
.fcard .fh{display:flex;align-items:center;gap:6px;font-size:11.5px;font-weight:700;color:#495057;margin-bottom:7px;padding-bottom:5px;border-bottom:1px solid #f1f3f5}
.fcard .fh .fhd{font-weight:400;color:#adb5bd;font-size:10.5px;margin-left:auto}
.chips{display:flex;flex-wrap:wrap;gap:5px}
.chip{display:inline-flex;align-items:center;gap:5px;padding:3px 9px;border-radius:14px;border:1px solid transparent;cursor:pointer;font-size:11px;user-select:none;white-space:nowrap}
.chip:hover{filter:brightness(.95)}.chip input{margin:0;cursor:pointer}
.chip.on{box-shadow:inset 0 0 0 2px currentColor;font-weight:600}
.chip .cnt{font-size:9.5px;background:rgba(0,0,0,.09);border-radius:8px;padding:0 5px;font-weight:700}
.det_A{background:#d3f9d8;color:#2b8a3e}.det_amb{background:#fff3bf;color:#b08900}.det_B{background:#d0ebff;color:#1971c2}.det_C{background:#ffe3e3;color:#e8590c}.det_incompat{background:#ffd6e7;color:#c2255c}.det_other{background:#f1f3f5;color:#868e96}
.detail{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:15px;min-height:76vh}
.zone{margin:20px 0 6px;padding:8px 14px;border-radius:7px;font-size:14px;font-weight:700;background:#e7f5ff;color:#1971c2;border-left:5px solid #1971c2}
.detail h3{position:sticky;top:0;background:#fff;margin:-15px -15px 8px;padding:12px 15px 8px;z-index:5;border-bottom:1px solid #f1f3f5;border-radius:8px 8px 0 0}
.kv{display:flex;gap:10px;flex-wrap:wrap;margin:6px 0}.kv .b{background:#f1f3f5;border-radius:6px;padding:5px 10px;font-size:11.5px}
table{border-collapse:collapse;font-size:11.5px;margin:7px 0}th,td{border:1px solid #dee2e6;padding:3px 8px}th{background:#f1f3f5}
.note{font-size:11.5px;color:#868e96}.mono{font-family:ui-monospace,Menlo,monospace}
"""

JS = r"""
function bootWS(){
const D=window.__DATA__;
const TT={'linear(全直系)':'t_linear','branched(直系+姊妹)':'t_branched','star(全姊妹)':'t_star','single':'t_single','germline_only':'t_single'};
const el=id=>document.getElementById(id);
const PAL=['#1c7ed6','#37b24d','#f59f00','#ae3ec9','#e8590c','#1098ad','#f03e3e','#adb5bd','#7048e8'];
function bars(o,m){let vs=Object.values(o),tot=vs.reduce((a,b)=>a+b,0)||1,mx=Math.max(...vs,1);return Object.entries(o).sort((a,b)=>b[1]-a[1]).slice(0,m||9).map(([k,v],idx)=>`<div class="bar"><i style="width:${Math.max(3,78*v/mx)}px;background:${PAL[idx%PAL.length]}"></i>${k}: <b>${v}</b> (${(100*v/tot).toFixed(1)}%)</div>`).join('')+`<div class="bar" style="color:#868e96">— 合計 ${tot} (100%)</div>`}
function pie(o){let es=Object.entries(o).sort((a,b)=>b[1]-a[1]),tot=es.reduce((a,b)=>a+b[1],0)||1,cx=30,cy=30,r=26,ri=14,a0=-Math.PI/2,p='';es.forEach(([k,v],i)=>{let a1=a0+2*Math.PI*v/tot,lg=(a1-a0)>Math.PI?1:0,x0=cx+r*Math.cos(a0),y0=cy+r*Math.sin(a0),x1=cx+r*Math.cos(a1),y1=cy+r*Math.sin(a1),xi1=cx+ri*Math.cos(a1),yi1=cy+ri*Math.sin(a1),xi0=cx+ri*Math.cos(a0),yi0=cy+ri*Math.sin(a0);p+=`<path d="M${x0.toFixed(1)} ${y0.toFixed(1)} A${r} ${r} 0 ${lg} 1 ${x1.toFixed(1)} ${y1.toFixed(1)} L${xi1.toFixed(1)} ${yi1.toFixed(1)} A${ri} ${ri} 0 ${lg} 0 ${xi0.toFixed(1)} ${yi0.toFixed(1)} Z" fill="${PAL[i%PAL.length]}"><title>${k}: ${v} (${(100*v/tot).toFixed(1)}%)</title></path>`;a0=a1});return `<svg width="60" height="60" viewBox="0 0 60 60" style="flex:0 0 auto">${p}</svg>`}
function pieBars(o,m){return `<div style="display:flex;gap:9px;align-items:flex-start">${pie(o)}<div style="flex:1;min-width:0">${bars(o,m)}</div></div>`}
el('s_topo').innerHTML=pieBars(D.stats.topology_type);el('s_clust').innerHTML=pieBars(Object.fromEntries(Object.entries(D.stats.n_clusters).map(([k,v])=>['c='+k,v])));
el('s_det').innerHTML=pieBars(D.stats.determinacy);el('s_root').innerHTML=pieBars(D.stats.n_roots);
// R4: 建樹位點分布(第5卡) — 區級從 detail 算 + sSNV 位點宇宙從 accounting
(function(){let dt=D.detail,ge2=dt.length,eq2=dt.filter(r=>r.n_sSNV==2).length,eq3=dt.filter(r=>r.n_sSNV==3).length,ge4=dt.filter(r=>r.n_sSNV>=4).length,strict=dt.filter(r=>(r.topology_type||'').split('(')[0]!=='single').length;
 let rows=[['可建樹區 n_sSNV≥2',ge2,'#1c7ed6'],['　┣ 恰 2 sSNV',eq2,'#4dabf7'],['　┣ 恰 3 sSNV',eq3,'#4dabf7'],['　┗ ≥4 sSNV',ge4,'#4dabf7'],['真多節點樹 lin+br+star',strict,'#37b24d']];
 let html=rows.map(([k,v,c])=>`<div class="bar"><i style="width:${Math.max(3,78*v/(ge2||1))}px;background:${c}"></i>${k}: <b>${v}</b></div>`).join('');
 let a=D.accounting;
 if(a){html+=`<div class="bar" style="color:#868e96;margin-top:4px;border-top:1px solid #f1f3f5;padding-top:3px">全 sSNV 位點宇宙 <b>${(a.universe_total||0).toLocaleString()}</b></div><div class="bar">linked <b>${a.buckets.linked.pct}%</b> 可建樹｜單位點 <b>${a.single_pct}%</b></div>`;}
 else{html+=`<div class="note" style="margin-top:4px">單位點數待 single_snv_accounting</div>`;}
 el('s_nsnv').innerHTML=html;})();
// R4: 統計卡 popup(放大 pie+全 bin+名詞字典),onclick 在 bootWS 閉包內→永遠對應當前 D(per-sample)
const STAT_DICT={topology_type:{title:'拓樸型態',desc:'每區 read 群在系統發生樹上的形狀(只計 n_sSNV≥2 區)。',items:{'single':'單群:reads 全塌成一個基因型,無分支','linear(全直系)':'全直系鏈 germline→A→AB→…','branched(直系+姊妹)':'有姊妹分支(同層平行 subclone)','star(全姊妹)':'全姊妹:多條從 germline 各自分出','germline_only':'只有 germline'}},n_clusters:{title:'群數 c',desc:'區內 distinct population(細胞狀態)數;perfect-phylogeny 下 ≤ k+1。',items:{}},determinacy:{title:'determinacy 可辨識性',desc:'樹「存在」≠「能辨識是哪棵」。',items:{'A_determined(單分子向量)':'單分子向量唯一可辨識','A_ambiguous_order(缺中間群)':'缺中間群→累積順序未定','B_pairwise_structure':'pairwise 拼接,非單分子整跨','C_underdetermined':'多樹相容,欠定','incompatible':'四配子違反→成環,無法成單一樹','other':'單群無分支'}},n_roots:{title:'HP 根數',desc:'somatic 事件散在幾條 germline 單倍型。≥2 = 跨 HP(allelic,非 subclone)。',items:{}}};
window.openStatModal=function(which){let o=D.stats[which];if(which==='n_clusters')o=Object.fromEntries(Object.entries(o).map(([k,v])=>['c='+k,v]));let d=STAT_DICT[which]||{title:which,desc:'',items:{}};let tot=Object.values(o).reduce((a,b)=>a+b,0)||1;let dict=Object.entries(d.items||{}).filter(([k])=>o[k]!=null).map(([k,v])=>`<div style="font-size:11.5px;margin:3px 0"><b class="mono">${k}</b> — ${v}</div>`).join('');el('statmodal_body').innerHTML=`<h3 style="margin-top:0">${d.title}（${(window.__DATA__&&document.querySelector('.stab.active')?document.querySelector('.stab.active').dataset.s:'')}；合計 ${tot}）</h3><div class="note" style="margin-bottom:10px">${d.desc}</div><div style="display:flex;gap:24px;align-items:flex-start;flex-wrap:wrap"><div style="transform:scale(1.7);transform-origin:top left;margin:18px 40px 40px 8px">${pie(o)}</div><div style="flex:1;min-width:280px">${bars(o,99)}</div></div>${dict?`<div style="margin-top:12px;border-top:1px solid #f1f3f5;padding-top:9px"><b style="font-size:12.5px">類別說明</b><div style="margin-top:5px">${dict}</div></div>`:''}`;el('statmodal').style.display='flex'};
window.closeStatModal=function(){el('statmodal').style.display='none'};
[['s_topo','topology_type'],['s_clust','n_clusters'],['s_det','determinacy'],['s_root','n_roots']].forEach(([id,key])=>{let sc=el(id).closest('.scard');if(sc)sc.onclick=()=>openStatModal(key)});
(function(){let a=D.accounting,u=el('universe');if(!u)return;if(!a){u.innerHTML='';return;}let B=a.buckets||{},gg=(o,k)=>(o&&o[k]!=null?o[k]:0),nf=x=>(x||0).toLocaleString();
 u.innerHTML=`<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:11px 14px;margin:10px 0;font-size:12.5px"><b>全 sSNV 宇宙帳本（${nf(a.universe_total)} = TP ${nf(a.tp)}（${(100*a.tp/(a.universe_total||1)).toFixed(1)}%）+ FP ${nf(a.fp)}（${(100*a.fp/(a.universe_total||1)).toFixed(1)}%）；隨樣本變）</b><br><div style="display:flex;gap:18px;flex-wrap:wrap;margin-top:6px"><span>🟢 <b>linked ${nf(gg(B.linked,'n'))}（${gg(B.linked,'pct')}%）</b> 可建樹</span><span>🟡 <b>underpowered ${nf(gg(B.underpowered,'n'))}（${gg(B.underpowered,'pct')}%）</b> 有 partner 無共讀→加深覆蓋可救</span><span>⚪ <b>isolated ${nf(gg(B.isolated,'n'))}（${gg(B.isolated,'pct')}%）</b> read-span 內無 partner</span><span class="note">（三桶加總＝宇宙 ✓）</span></div>${(a.derived_from_census||a.fp_truth_sparse||a.cn_all_neutral)?`<div class="note" style="margin-top:6px;color:#a37200">⚠ ${a.derived_from_census?'帳本由 per_sSNV_census 衍生':''}${a.fp_truth_sparse?'；FP truth 稀疏（無外部 truth set→TP/FP 標籤弱，勿當判別力佐證）':''}${a.cn_all_neutral?'；CN 未併 census（by_cn 全 neutral）':''}</div>`:''}</div>`;})();
// 克隆樹判讀圖鑑(R3+R7) — 固定取 HCC1395,渲染一次,與分頁無關
function renderTeaching(){
 let TS=window.__SAMPLES__||{}, H=TS['HCC1395']||TS[Object.keys(TS).find(k=>TS[k]&&TS[k].chr17)];
 let host=el('teachbody'); if(!host)return;
 if(!H||!H.chr17){host.innerHTML='<div class="note">教學資料(HCC1395 chr17)不可用</div>';return;}
 let c=H.chr17, byreg={}; (H.detail||[]).forEach(r=>byreg[r.region]=r);
 let germ=(c.populations.find(p=>p.muts=='germline')||{}).reads||0, germVec=(c.populations.find(p=>p.muts=='germline')||{}).vec;
 let treeVecs=new Set((c.edges||[]).reduce((a,e)=>a.concat(e),[]));
 let st=c.snvs.map(s=>`<tr><td class="mono"><b>${s.S}</b></td><td class="mono">${s.pos}</td><td>${s.change}</td><td>${s.role}</td><td>VAF ${s.vaf}</td><td>${s.hp}</td><td>${s.src}</td><td>${s.somatic_confirmed?'✓':'✗ normal有ALT'}</td></tr>`).join('');
 let pp=c.populations.map(p=>{let onT=treeVecs.has(p.vec)||p.vec==germVec;return `<tr style="${onT?'':'background:#fff0f6'}"><td class="mono">${p.vec}</td><td>${p.muts}</td><td>${p.reads}</td><td>${p.pct}%</td><td class="note">${onT?'':'⚠ 噪聲丟棄(無上樹)'}</td></tr>`}).join('');
 let mm=c.sig_cpg.slice(0,8).map(x=>`<tr><td class="mono">${x.m}</td><td>${x.cpg}</td><td>${x.L1}</td><td>${x.L2}</td><td>${x.dbeta}</td></tr>`).join('');
 let hero=`<div class="kv"><div class="b">locus ${c.locus}</div><div class="b">ctx ${c.genome_ctx}</div><div class="b">拓樸 ${c.topology_type}</div><div class="b">噪聲丟 ${c.dropped_noise} reads</div></div>
  <div style="background:#e7f5ff;border:1px solid #a5d8ff;border-radius:6px;padding:7px 11px;font-size:11.5px;line-height:1.85;margin:6px 0">① <b>germline 根</b>(灰框·全R向量·無somatic·起點) ② <b>直系鏈</b> germline→RRAR(+S3 α祖先)→RAAA(+S2+S4 β後代),S2/S4 <b>co_linked</b> 同步加(順序未定) ③ <b>姊妹分支</b> germline→ARRR(+S1):S1=48357368 為 <b>ClairS FP</b>→演算法正確隔離成獨立 sibling、未塞進真鏈 ④ <b>+S</b> 對到下方 S 表 ⑤ <b>甲基 cis-ASM</b>:m 對齊 genotype 軸(L1=RRAR vs L2=RAAA)非獨立 lineage(見下方展板)</div>
  <b>克隆樹</b>${tree(c.edges,Object.fromEntries(c.populations.map(p=>[p.vec,p.reads])),c.populations.length,'H1',germ,0)}
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-top:6px"><div><b>S 位點(somatic sSNV)</b><table><tr><th>S</th><th>pos</th><th>變異</th><th>角色</th><th>VAF</th><th>HP</th><th>TP/FP</th><th>som</th></tr>${st}</table></div><div><b>細胞群(粉底=噪聲未上樹)</b><table><tr><th>向量</th><th>突變</th><th>reads</th><th>%</th><th></th></tr>${pp}</table></div></div>
  <b>甲基差異位點 m(⚠ cis-ASM:對齊 genotype 軸非獨立 lineage)</b><table><tr><th>m</th><th>CpG</th><th>L1 β</th><th>L2 β</th><th>Δβ</th></tr>${mm}</table>`;
 let card=(t,note,extra)=>`<div style="border:1px solid #dee2e6;border-radius:8px;padding:9px 11px;background:#fff"><div style="font-weight:700;font-size:12.5px;color:#1971c2">${t}</div><div class="note" style="margin:3px 0 6px">${note}</div>${extra}</div>`;
 let TR=r=>tree(r.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes);
 let cards='';
 let rA=byreg['chr1:24630300-24635403']; if(rA)cards+=card('卡A · 姊妹分支 sibling','AA=0(兩突變從不共現)+ 同 HP → 平行 subclone、樹上同層分叉(橙框)。',TR(rA));
 let rB=byreg['chr1:5816053-5816054']; if(rB)cards+=card('卡B · co_linked 完美共現','只見 AA、RA=AR=0 → 兩突變永遠同進退、綁同一事件、無法內部排序。',TR(rB));
 let rC=byreg['chr1:12724814-12732374']; if(rC){let cd=enumCandidates(rC);let cl=cd?'<div style="margin-top:6px;font-size:11px;background:#f8f0fc;border-radius:5px;padding:5px 8px">候選累積序(共 '+cd.trueCount+' 種·等機率·中間群未觀察):<br>'+cd.cands.map((cc,ix)=>(ix+1)+'. '+cc.map(a=>((rC.node_paths||{})[a.node]||a.node)+': '+a.order.join('→')).join('；')).join('<br>')+'</div>':'';cards+=card('卡C · 缺中間群 ambiguous','節點一次獲≥2 變異、中間群未觀察 → 累積順序未定(黃框)。可能順序等機率列出:',TR(rC)+cl);}
 let rD=byreg['chr2:97211773-97229016']; if(rD){let cf=fourGamete(rD.populations);cards+=card('卡D · 四配子衝突 four-gamete','RR/RA/AR/AA 四格全有 → 不相容單一樹。錨 RR=germline → AA=雙突變最遠。',cf.length?'<table><tr><th>對</th><th>RR</th><th>RA</th><th>AR</th><th>AA</th></tr>'+cf.map(x=>'<tr><td class="mono"><b>'+x.pair+'</b></td><td>'+x.g.RR+'</td><td>'+x.g.RA+'</td><td>'+x.g.AR+'</td><td>'+x.g.AA+'</td></tr>').join('')+'</table>':'<div class="note">(此區四配子未觸發)</div>');}
 let rE=byreg['chr7:100784203-100798014']; if(rE)cards+=card('卡E · CROSS-HP allelic(對比卡A)','n_roots≥2:突變散在 H1/H2 兩單倍型=各自染色體突變(allelic)非 subclone → 畫兩棵分開 HP 樹。',posTree(rE)||TR(rE));
 let mt=window.__CHR17TREE__, exhibit='';
 if(mt&&mt.cross_clu2_x_lineage){let cx=mt.cross_clu2_x_lineage,lins=['L0','L1','L2','other'];
  let xrows=Object.keys(cx).sort().map(cl=>'<tr><td><b>甲基群 '+cl+'</b></td>'+lins.map(L=>'<td style="'+((cx[cl][L]||0)>=10?'background:#ffe3e3;font-weight:700':'')+'">'+(cx[cl][L]||0)+'</td>').join('')+'</tr>').join('');
  let ax=mt.axis_sig_cpg_count||{},axmax=Math.max.apply(null,Object.values(ax).concat([1]));
  let axbars=Object.entries(ax).map(([k,v])=>'<div class="bar"><i style="width:'+Math.max(3,90*v/axmax)+'px;background:'+(k.indexOf('lineage')>=0?'#e8590c':(k.indexOf('HP')>=0?'#868e96':'#1c7ed6'))+'"></i>'+k+': <b>'+v+'</b></div>').join('');
  exhibit='<div style="background:#fff5f5;border:1px solid #ffc9c9;border-radius:8px;padding:10px 13px;margin-top:12px"><b>🧬 甲基 read 距離+分群「能/不能做什麼」(chr17 worked,'+mt.n_reads+' reads × '+mt.n_cpg+' CpG)</b><div class="note" style="margin:4px 0">把 read 用甲基距離分 k=2 群,再對遺傳 lineage 交叉:</div><table><tr><th>甲基分群＼遺傳lineage</th>'+lins.map(L=>'<th>'+L+'</th>').join('')+'</tr>'+xrows+'</table><div style="color:#c92a2a;font-weight:600;margin:6px 0;font-size:11.5px">🔴 甲基群「1」同時含 L1 與 L2 各 19 → <b>甲基分群 ≠ 遺傳 lineage</b>(cis-ASM double-dip 本體);甲基不能 recover subclone。</div><div class="note">各軸顯著 CpG 數(甲基「對齊」哪個軸):</div>'+axbars+'<div class="note" style="margin-top:5px">→ α(genotype)軸最強、lineage 軸弱、HP 軸 0 ⇒ 甲基是 <b>ASM 存在性偵測器</b>非 lineage 排序器。基因組級 PERMANOVA 740 區僅 1 testable、recover=0=統計死。<b>可用=負篩/L3弱旗標/教學;不可用=分群器/排序器/定群器</b>(06-28 cis-control 裁決)。</div></div>';}
 host.innerHTML='<div class="note" style="margin-bottom:8px">節點=一種細胞狀態(read 群);往下=直系、同層分叉=姊妹;邊上 +S=新增 somatic 變異;根=germline。<b>以下範例固定取自 HCC1395,與上方分頁樣本無關。</b></div><h4 style="margin:6px 0">① HERO:chr17:48357368-48365161 逐元素判讀（綁 production 真實區·4 sSNV）</h4>'+hero+'<h4 style="margin:16px 0 6px">② 五種 SNV 關係配套真例(各一真實區)</h4><div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:8px">'+cards+'</div>'+exhibit+'<h4 style="margin:16px 0 6px">③ 標籤圖例</h4><div class="note" style="line-height:1.9"><span style="color:#1565c0;font-weight:700">●直系 vertical</span> 後代多帶1變異 ｜ <span style="color:#d9480f;font-weight:700">●姊妹 sibling</span> 同父平行 subclone ｜ <span style="color:#5f3dc4;font-weight:700">co_linked</span> 完美共現綁同事件 ｜ <span style="color:#e8590c;font-weight:700">⚠缺中間群</span> 順序未定 ｜ <span style="color:#c2255c;font-weight:700">四配子違反</span> 成環無法成樹 ｜ <span style="color:#2b8a3e;font-weight:700">+S</span> 新增變異';
}
// R6 草圖: 無法定位群偵測總覽(per-sample;「群存在但定不出位置/順序」)
(function(){let u=el('unlocatable');if(!u)return;let sm=(D.scoring&&D.scoring.summary)||{},sd=sm.situation_dist;if(!sd){u.innerHTML='';return;}
 let tot=Object.values(sd).reduce((a,b)=>a+b,0)||1;
 let order=[['已確定','#37b24d'],['pairwise 拼接','#74c0fc'],['多樹相容(欠定)','#adb5bd'],['跨HP(兩棵樹)','#f59f00'],['順序 2-3 順位待定','#ffd43b'],['衝突(成環)','#f03e3e']];
 let seg=order.filter(([k])=>sd[k]).map(([k,c])=>`<div title="${k}: ${sd[k]}（${(100*sd[k]/tot).toFixed(1)}%）" style="width:${100*sd[k]/tot}%;background:${c}"></div>`).join('');
 let unloc=[['衝突(成環)','#f03e3e','四配子違反 → 無單一樹'],['跨HP(兩棵樹)','#f59f00','突變散在 H1/H2 → 非單一譜系、是兩棵樹'],['多樹相容(欠定)','#adb5bd','群存在但無 ≥2 ALT 群/缺對 → 無法排序(需深覆蓋)　←「某群存在卻定不出位置」正典'],['順序 2-3 順位待定','#ffd43b','缺中間群 → 累積順序未定']];
 let rows=unloc.filter(([k])=>sd[k]!=null).map(([k,c,d])=>`<div onclick="filterSituation('${k}')" style="display:flex;gap:8px;align-items:baseline;margin:3px 0;font-size:11.5px;cursor:pointer;border-radius:4px;padding:2px 6px" onmouseover="this.style.background='#fff0f6'" onmouseout="this.style.background=''"><span style="display:inline-block;width:11px;height:11px;border-radius:3px;background:${c};flex:0 0 auto"></span><b style="width:120px;flex:0 0 auto">${k}</b><b style="color:${c};width:46px;flex:0 0 auto">${sd[k]}</b><span class="note">${d}</span><span class="note" style="margin-left:auto;color:#1971c2;flex:0 0 auto;white-space:nowrap">▸ 篩選佇列</span></div>`).join('');
 window.filterSituation=function(sit){let s=el('q_sit');if(!s)return;if([...s.options].some(o=>o.value===sit)){s.value=sit;if(typeof renderQ==='function')renderQ();}let qq=el('queue');if(qq)qq.scrollIntoView({behavior:'smooth',block:'center'});};
 let q=(D.scoring.queue||[]),withp=q.filter(x=>x.parsimony_first_rank_prob!=null),hi=withp.filter(x=>x.parsimony_first_rank_prob>=0.7).length;
 let nunloc=(sd['衝突(成環)']||0)+(sd['跨HP(兩棵樹)']||0)+(sd['多樹相容(欠定)']||0)+(sd['順序 2-3 順位待定']||0);
 u.innerHTML=`<div class="zone" style="background:#fff0f6;color:#c2255c;border-left-color:#f06595">🔍 無法定位群偵測（草圖）—「群存在但定不出位置/順序」(此樣本 ${nunloc} 區)</div>
  <div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:10px 13px;font-size:12.5px">
  <div style="display:flex;height:16px;border-radius:5px;overflow:hidden;margin-bottom:7px">${seg}</div>
  <div class="note" style="margin-bottom:6px">合計 ${tot} 區｜🟢已確定 ${sd['已確定']||0}｜🔵pairwise 拼接 ${sd['pairwise 拼接']||0}(可建樹非單分子)｜<b style="color:#c2255c">↓ 4 種無法定位</b></div>${rows}
  <div style="margin-top:8px;border-top:1px solid #f1f3f5;padding-top:7px;font-size:11.5px"><b>機率(誠實兩軌)</b>：遺傳 parsimony 高信度(≥0.7)= <b>${withp.length?hi:'未回填'}</b>${withp.length?` 區/${withp.length} 有值 → ${hi?'':'連遺傳都信心不足'}`:'（此樣本待上游回填）'}　｜　🧬 甲基：SAME-HP <b>不給機率</b>(cis-ASM double-dip)、CROSS-HP 弱(~35%)、乾淨可用≈0（06-28）</div>
  <div class="note" style="margin-top:5px;color:#a37200">⚠ 「不同的樹」全列舉(ranked 替代整樹)需上游 <b>enumerate_candidate_trees</b>（未實作）→ 此面板誠實顯示「什麼定不出來」,不假裝在列舉答案。</div></div>`;})();
if(!window.__teachDone){renderTeaching();window.__teachDone=true;}
// S-label 化基因型向量
function sLabels(g){let s=[...g].map((c,i)=>c=='A'?('S'+(i+1)):null).filter(Boolean);return s.length?s.join('+'):'germline'}
function gainedS(parent,child){let g=[];for(let i=0;i<child.length;i++){if(child[i]=='A'&&(!parent||parent[i]!='A'))g.push('S'+(i+1))}return g}
// 四配子檢定:對每對 sSNV(i,j)數 read 在兩位點的 RR/RA/AR/AA。四格全有=incompatible(無限位點違反)。RR=germline(normal錨REF)、AA=雙突變(最遠)。
function fourGamete(pops){let vs=Object.keys(pops||{}).filter(v=>v);if(!vs.length)return [];let L=vs[0].length,out=[];for(let i=0;i<L;i++)for(let j=i+1;j<L;j++){let g={RR:0,RA:0,AR:0,AA:0};vs.forEach(v=>{let k=v[i]+v[j];if(g[k]!=null)g[k]+=pops[v]});if(g.RR&&g.RA&&g.AR&&g.AA)out.push({pair:'S'+(i+1)+'–S'+(j+1),i:i+1,j:j+1,g:g})}return out}
// 候選累積順序列舉:缺中間群節點(獲得≥2突變)→列舉可能累積序(等機率,中間群未觀察→無讀數可分)。前端可做;ranked 替代整樹需上游 enumerate_candidate_trees。
function fact(n){let r=1;for(let i=2;i<=n;i++)r*=i;return r}
function enumCandidates(r){let par={};(r.edges||[]).forEach(([p,c])=>{if(p!='ROOT')par[c]=p});let amb=[];Object.keys(r.populations).forEach(n=>{let g=gainedS(par[n]||null,n);if(g.length>=2)amb.push({node:n,gained:g})});if(!amb.length)return null;function perms(a){if(a.length<=1)return [a];let o=[];a.forEach((x,i)=>{perms(a.slice(0,i).concat(a.slice(i+1))).forEach(p=>o.push([x].concat(p)))});return o}let trueCount=amb.reduce((acc,a)=>acc*fact(a.gained.length),1);let perNode=amb.map(a=>({node:a.node,orders:a.gained.length<=4?perms(a.gained):[a.gained]}));let cands=[[]];perNode.forEach(pn=>{let nx=[];cands.forEach(c=>pn.orders.forEach(o=>{if(nx.length<24)nx.push(c.concat([{node:pn.node,order:o}]))}));cands=nx});return {cands:cands,trueCount:trueCount,bigNode:amb.some(a=>a.gained.length>4)}}
function candCard(cs,idx,r){let c=cs[idx],np=r.node_paths||{};let body=c.map(a=>`<div style="margin:3px 0"><b class="mono" style="color:#1971c2">${np[a.node]||a.node}</b>（${a.node}）：${a.order.map((s,i)=>`<span style="background:#ebfbee;border:1px solid #2f9e44;border-radius:9px;padding:1px 7px;color:#2b8a3e;font-weight:700;margin:0 1px">${i+1}.${s}</span>`).join('→')}</div>`).join('');return `<div style="display:flex;align-items:center;gap:10px"><button onclick="candNav(-1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px">◀</button><div style="flex:1">${body}<div class="note" style="margin-top:3px">第 ${idx+1} / ${cs.length} 個候選</div></div><button onclick="candNav(1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px">▶</button></div>`}
window.candNav=function(d){let s=window.__cand;if(!s)return;s.idx=(s.idx+d+s.cands.length)%s.cands.length;let box=document.getElementById('candbox');if(box)box.innerHTML=candCard(s.cands,s.idx,s.r)}
// Part B: 上游 enumerate 的完整候選樹 carousel(虛擬中間節點;equiprobable 誠實標)
function candTreeCard(ctr,idx,r){let c=ctr.candidate_set[idx];let th=tree(c.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes);return `<div style="display:flex;align-items:flex-start;gap:8px"><button onclick="ctNav(-1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px;margin-top:46px">◀</button><div style="flex:1;min-width:0">${th}<div class="note" style="margin-top:4px">候選樹 <b>${idx+1} / ${ctr.candidate_set.length}</b>　softmax ${(c.softmax_prob*100).toFixed(0)}%${ctr.equiprobable?'　⚖ <b>等機率</b>(無 read 證據可排)':''}　虛擬中間節點(0 reads·未觀察): ${c.virtual_nodes.length?c.virtual_nodes.join('、'):'無'}</div></div><button onclick="ctNav(1)" style="font-size:15px;padding:3px 11px;cursor:pointer;border-radius:5px;margin-top:46px">▶</button></div>`}
window.ctNav=function(d){let s=window.__ctr;if(!s)return;let n=s.cs.candidate_set.length;s.idx=(s.idx+d+n)%n;let box=document.getElementById('cttreebox');if(box)box.innerHTML=candTreeCard(s.cs,s.idx,s.r)}
function tree(edges,popcount,nc,hp,germR,np,ambig){np=np||{};ambig=ambig||0;
 let germKey=Object.keys(popcount).find(k=>k&&/^R+$/.test(k)); // germline=全-R向量;germline_reads欄位不可靠(很多區為0),改用 populations 全-R count 與下方表格一致
 let germN=germKey!=null?popcount[germKey]:(germR||0);
 let ch={},par={},all=new Set();
 (edges||[]).forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT'){all.add(p);par[c]=p}});
 if(!all.size)return `<div class="note">單群／germline-only（germline ${germN} reads），無分支樹</div>`;
 let sib={};Object.keys(ch).forEach(p=>{(ch[p]||[]).forEach(c=>{sib[c]=(ch[p].length>=2)})});
 let depth={},pos={},leaf=0,seen={};
 function lay(n,dp){if(seen[n])return pos[n]!=null?pos[n]:0;seen[n]=1;depth[n]=dp;let k=(ch[n]||[]).filter(x=>!seen[x]).sort();if(!k.length){pos[n]=leaf++;return pos[n]}let xs=k.map(x=>lay(x,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]} // seen 防護:成環邊跳過,避免無限遞迴 stack overflow
 let roots=ch['ROOT']||[];roots.forEach(r=>lay(r,1));
 let gx=roots.length?(Math.min(...roots.map(r=>pos[r]))+Math.max(...roots.map(r=>pos[r])))/2:0;
 let nodes=[...all],md=Math.max(...nodes.map(n=>depth[n]),1);
 let NW=204,NH=80,GX=50,GY=128;
 let W=Math.max(450,(leaf||1)*(NW+GX)),H=108+md*GY;
 let X=p=>34+p*(NW+GX)+NW/2,Y=dp=>56+dp*GY;
 let totR=Object.values(popcount).reduce((a,b)=>a+b,0)||1; // popcount 已含 germline 向量,勿再加 germR(否則 germline 雙算→比例與下方表格對不上)
 let relSet=new Set();
 let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" style="font-family:ui-monospace,Menlo,monospace">`;
 roots.forEach(r=>{let g=gainedS(null,r),mx=(X(gx)+X(pos[r]))/2,my=(Y(0)+Y(1))/2,w=Math.max(44,g.join('+').length*7+12);
  s+=`<line x1="${X(gx)}" y1="${Y(0)+NH/2}" x2="${X(pos[r])}" y2="${Y(1)-NH/2}" stroke="#ced4da" stroke-width="1.6"/>`;
  if(g.length)s+=`<rect x="${mx-w/2}" y="${my-10}" width="${w}" height="19" rx="9" fill="#ebfbee" stroke="#2f9e44"/><text x="${mx}" y="${my+4}" text-anchor="middle" font-size="11" fill="#2b8a3e" font-weight="700">+${g.join('+')}</text>`;}); // germline→第一代 edge 也標 +S(原缺)
 (edges||[]).forEach(([p,c])=>{if(p!='ROOT'){
  let g=gainedS(p,c),mx=(X(pos[p])+X(pos[c]))/2,my=(Y(depth[p])+Y(depth[c]))/2;
  s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+NH/2}" x2="${X(pos[c])}" y2="${Y(depth[c])-NH/2}" stroke="#ced4da" stroke-width="1.6"/>`;
  let w=Math.max(44,g.join('+').length*7+12);if(g.length)s+=`<rect x="${mx-w/2}" y="${my-10}" width="${w}" height="19" rx="9" fill="#ebfbee" stroke="#2f9e44"/><text x="${mx}" y="${my+4}" text-anchor="middle" font-size="11" fill="#2b8a3e" font-weight="700">+${g.join('+')}</text>`;
 }});
 let gp=(100*germN/totR).toFixed(0);
 s+=`<rect x="${X(gx)-NW/2}" y="${Y(0)-NH/2}" width="${NW}" height="${NH}" rx="12" fill="#f1f3f5" stroke="#868e96" stroke-width="2.5"/>`;
 s+=`<text x="${X(gx)}" y="${Y(0)-NH/2+24}" text-anchor="middle" font-size="14" font-weight="700" fill="#495057">⌂ germline（${hp||'根'}）</text>`;
 s+=`<text x="${X(gx)}" y="${Y(0)-NH/2+44}" text-anchor="middle" font-size="11" fill="#868e96">無 somatic 變異 · 起點</text>`;
 s+=`<text x="${X(gx)}" y="${Y(0)-NH/2+66}" text-anchor="middle" font-size="12.5" font-weight="600" fill="#212529">${germN} reads · ${gp}%</text>`;
 nodes.forEach(n=>{
  let cnt=popcount[n]||0,pct=(100*cnt/totR).toFixed(0),lab=np[n]||'—';
  let g=gainedS(par[n],n);
  let isSib=sib[n], isMulti=g.length>=2, isAmbig=isMulti&&ambig>0, isCo=isMulti&&!isAmbig;
  if(isSib)relSet.add('姊妹(sibling)'); else if(g.length)relSet.add('直系(vertical)');
  if(isAmbig)relSet.add('⚠缺中間群(順序未定)'); else if(isCo)relSet.add('co_linked(完美共現)');
  let multiTag=isAmbig?' · ⚠缺中間群(順序未定)':(isCo?' · co_linked(完美共現)':'');
  let rel=g.length?((isSib?'姊妹分支(sibling)':'直系(vertical)')+multiTag):'（與父同型）';
  let gtxt=g.length?('獲得 +'+g.join('+')):'';
  let x=X(pos[n]),y=Y(depth[n]);
  let fill=isSib?'#fff4e6':'#e7f5ff', stroke=isSib?'#e8590c':'#1c7ed6';
  if(isAmbig){fill='#fff9db';stroke='#f08c00';}
  s+=`<rect x="${x-NW/2}" y="${y-NH/2}" width="${NW}" height="${NH}" rx="12" fill="${fill}" stroke="${stroke}" stroke-width="2.5"/>`;
  s+=`<text x="${x}" y="${y-NH/2+22}" text-anchor="middle" font-size="16" font-weight="800" fill="#1565c0">${lab}</text>`;
  s+=`<text x="${x}" y="${y-NH/2+39}" text-anchor="middle" font-size="10" font-weight="700" fill="${isSib?'#d9480f':(isAmbig?'#e8590c':'#2b8a3e')}">${rel}</text>`;
  s+=`<text x="${x}" y="${y-NH/2+53}" text-anchor="middle" font-size="10.5" font-weight="600" fill="#2b8a3e">${gtxt}</text>`;
  s+=`<text x="${x}" y="${y-NH/2+66}" text-anchor="middle" font-size="9.5" fill="#495057">基因型 ${n}（=${sLabels(n)}）</text>`;
  s+=`<text x="${x}" y="${y-NH/2+78}" text-anchor="middle" font-size="12" font-weight="600" fill="#212529">${cnt} reads · ${pct}%</text>`;
 });
 s+='</svg>';
 let leg=`<details style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;padding:5px 11px;margin-top:5px;font-size:11px"><summary style="cursor:pointer;font-weight:600">🔖 SNV 關係圖例（點開）${relSet.size?'：此樹有 '+[...relSet].join('、'):''}</summary><div style="line-height:1.7;margin-top:5px">
 <span style="color:#1565c0;font-weight:700">●直系 vertical</span> 往下一層、後代多帶 1 變異(+S，藍框)<br>
 <span style="color:#d9480f;font-weight:700">●姊妹分支 sibling</span> 同一父不同分支、平行 subclone(橙框)<br>
 <span style="color:#5f3dc4;font-weight:700">co_linked(完美共現)</span> 一節點獲≥2 變異且區無 ambiguous → 兩變異綁同一事件<br>
 <span style="color:#e8590c;font-weight:700">⚠缺中間群(順序未定)</span> 一節點獲≥2 變異但區有 ambiguous → 跳>1突變、未觀察到中間群、累積順序未定（黃框；<b>非</b>確定 co_linked）<br>
 <span style="color:#2b8a3e;font-weight:700">+S</span> 該分支新增 somatic 變異</div></details>`;
 return s+leg;
}
// 2-root: 位置樹按 HP 分兩棵
function hasCycleEdges(edges){let ch={};edges.forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c)});let col={},cyc=false;function dfs(u){col[u]=1;for(let v of (ch[u]||[])){if(col[v]==1){cyc=true;return}if(!col[v])dfs(v);if(cyc)return}col[u]=2}for(let n of Object.keys(ch)){if(!col[n])dfs(n);if(cyc)break}return cyc}
function posTable(ns,vaf,hp,nested){let parent={};nested.forEach(([a,b])=>{(parent[b]=parent[b]||[]).push(a.split(':')[1])});let rows=[...ns].sort().map(n=>{let v=vaf[n],par=(parent[n]||[]).join('、');return `<tr><td class="mono">${n.split(':')[1]}</td><td>${v!=null?v:'?'}</td><td class="note">${par?'巢狀於 '+par:'—(根/無上游)'}</td></tr>`}).join('');return `<table style="font-size:10.5px"><tr><th>${hp} 位點</th><th>VAF</th><th>巢狀於</th></tr>${rows}</table>`}
function posTree(r){
 let byhp={};Object.entries(r.node_hp||{}).forEach(([p,h])=>{if(h=='H1'||h=='H2')(byhp[h]=byhp[h]||[]).push(p)});
 let hps=Object.keys(byhp).sort();if(hps.length<2)return '';
 return '<div style="display:flex;gap:20px;flex-wrap:wrap">'+hps.map(h=>{
  let ns=new Set(byhp[h]);let ned=(r.pos_nested||[]).filter(e=>ns.has(e[0])&&ns.has(e[1]));
  let body;
  if(hasCycleEdges(ned)){body='<div class="note" style="color:#c2255c;margin:3px 0">⚠ 此 '+h+' 位點 pairwise nested <b>成環/互指(incompatible)</b>→ 無法成單一樹,改列位點+VAF 表:</div>'+posTable(ns,r.pos_vaf||{},h,ned);}
  else{let hasp=new Set(ned.map(e=>e[1]));let edges=ned.map(e=>[e[0],e[1]]);[...ns].forEach(n=>{if(!hasp.has(n))edges.unshift(['ROOT',n])});body=posSVG(edges,ns,r.pos_vaf||{},h);}
  return '<div style="min-width:0;max-width:100%;overflow-x:auto"><b style="color:#d9480f">'+h+' 樹（'+ns.size+' 位點）</b>'+body+'</div>';
 }).join('')+'</div>';
}
function posSVG(edges,ns,vaf,hp){
 let ch={},all=new Set();edges.forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT')all.add(p)});
 if(!all.size)return '<div class="note">無結構</div>';
 let depth={},pos={},leaf=0,seen={};function lay(n,dp){if(seen[n])return pos[n]!=null?pos[n]:0;seen[n]=1;depth[n]=dp;let k=(ch[n]||[]).filter(x=>!seen[x]).sort();if(!k.length){pos[n]=leaf++;return pos[n]}let xs=k.map(x=>lay(x,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]} // seen 防護:pos_nested 成環(incompatible)時跳過,避免 stack overflow
 let roots=ch['ROOT']||[];roots.forEach(r=>lay(r,1));[...all].forEach(n=>{if(pos[n]==null)lay(n,1)});let gx=roots.length?(Math.min(...roots.map(r=>pos[r]))+Math.max(...roots.map(r=>pos[r])))/2:0;
 let nodes=[...all],md=Math.max(...nodes.map(n=>depth[n]),1),W=Math.max(200,(leaf||1)*108),H=64+md*70,X=p=>40+p*108,Y=dp=>26+dp*70;
 let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}"><circle cx="${X(gx)}" cy="${Y(0)}" r="15" fill="#fff" stroke="#495057"/><text x="${X(gx)}" y="${Y(0)+3}" text-anchor="middle" font-size="9">${hp} germ</text>`;
 roots.forEach(rt=>s+=`<line x1="${X(gx)}" y1="${Y(0)+15}" x2="${X(pos[rt])}" y2="${Y(1)-20}" stroke="#adb5bd"/>`);
 edges.forEach(([p,c])=>{if(p!='ROOT')s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+20}" x2="${X(pos[c])}" y2="${Y(depth[c])-20}" stroke="#adb5bd"/>`});
 nodes.forEach(n=>{let v=vaf[n];s+=`<rect x="${X(pos[n])-46}" y="${Y(depth[n])-18}" width="92" height="36" rx="5" fill="#fff4e6" stroke="#e8590c"/><text x="${X(pos[n])}" y="${Y(depth[n])-2}" text-anchor="middle" font-size="9" class="mono">${n.split(':')[1]}</text><text x="${X(pos[n])}" y="${Y(depth[n])+11}" text-anchor="middle" font-size="8" fill="#495057">VAF ${v!=null?v:'?'}</text>`});
 return s+'</svg>';
}
// filters + sort
let det=D.detail;
let chrsel=el('f_chr');chrsel.innerHTML='<option value="">全</option>';D.chroms.forEach(c=>{let o=document.createElement('option');o.value=c;o.textContent=c;chrsel.appendChild(o)});
const FACET_CFG={topology_type:{order:['linear(全直系)','branched(直系+姊妹)','star(全姊妹)','single','germline_only'],label:{'linear(全直系)':'linear 全直系','branched(直系+姊妹)':'branched 直系+姊妹','star(全姊妹)':'star 全姊妹','single':'single 單群','germline_only':'germline only'},cls:{'linear(全直系)':'t_linear','branched(直系+姊妹)':'t_branched','star(全姊妹)':'t_star','single':'t_single','germline_only':'t_single'}},determinacy:{order:['A_determined(單分子向量)','A_ambiguous_order(缺中間群)','B_pairwise_structure','C_underdetermined','incompatible','other'],label:{'A_determined(單分子向量)':'A 唯一·單分子','A_ambiguous_order(缺中間群)':'A 順序未定','B_pairwise_structure':'B 拼接','C_underdetermined':'C 多樹相容','incompatible':'✗ 成環','other':'— 單群無分支'},cls:{'A_determined(單分子向量)':'det_A','A_ambiguous_order(缺中間群)':'det_amb','B_pairwise_structure':'det_B','C_underdetermined':'det_C','incompatible':'det_incompat','other':'det_other'}},genome_ctx:{order:['arm','telomere','centromere'],label:{'arm':'arm 臂','telomere':'telomere 端粒','centromere':'centromere 著絲點'},cls:{'arm':'ctx_arm','telomere':'ctx_telomere','centromere':'ctx_centromere'}}};
['topology_type','determinacy','genome_ctx'].forEach(k=>{let c=el('cb_'+k);c.innerHTML='';let cfg=FACET_CFG[k],cnt={};det.forEach(r=>{cnt[r[k]]=(cnt[r[k]]||0)+1});let vals=cfg.order.filter(v=>cnt[v]!=null).concat(Object.keys(cnt).filter(v=>!cfg.order.includes(v)).sort());vals.forEach(o=>{let lab=document.createElement('label');lab.className='chip '+(cfg.cls[o]||'det_other');lab.title=o;lab.innerHTML='<input type="checkbox" value="'+o+'">'+(cfg.label[o]||o)+'<span class="cnt">'+cnt[o]+'</span>';c.appendChild(lab)});c.onchange=()=>{c.querySelectorAll('.chip').forEach(ch=>ch.classList.toggle('on',ch.querySelector('input').checked));render()}});
function cset(k){return new Set([...el('cb_'+k).querySelectorAll('input:checked')].map(x=>x.value))}
const SORT={coord:(a,b)=>a.chrom.localeCompare(b.chrom,undefined,{numeric:true})||a.start-b.start,nsnv:(a,b)=>b.n_sSNV-a.n_sSNV,nclust:(a,b)=>b.n_clusters-a.n_clusters||b.n_sSNV-a.n_sSNV,region:(a,b)=>a.region.localeCompare(b.region)};
function render(){
 let ch=el('f_chr').value,mc=+el('f_minc').value,tf=el('f_tpfp').value,loh=el('f_loh').checked,undef=el('f_undef').checked,q=el('f_q').value.trim(),so=el('f_sort').value;
 let tt=cset('topology_type'),dd=cset('determinacy'),gc=cset('genome_ctx');
 let tpfpok=r=>(tf=='all')||(tf=='tp'&&r.tp>0)||(tf=='fp'&&r.fp>0)||(tf=='both'&&r.tp>0&&r.fp>0);
 let f=det.filter(r=>(!ch||r.chrom==ch)&&(!tt.size||tt.has(r.topology_type))&&(!dd.size||dd.has(r.determinacy))&&(!gc.size||gc.has(r.genome_ctx))&&r.n_clusters>=mc&&tpfpok(r)&&(!loh||(r.cn=='loh'&&(r.haplotypes=='H1'||r.haplotypes=='H2')))&&(!undef||r.undefined)&&(!q||r.region.includes(q)));
 f.sort(SORT[so]||SORT.coord);if(el('f_sortdir').value=='rev')f.reverse();
 el('cnt').textContent=f.length+' 區';
 el('list').innerHTML=f.slice(0,700).map(r=>`<div class="row" data-i="${det.indexOf(r)}"><b>${r.region}</b> <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type.split('(')[0]}</span><span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span><br><span class="note">${r.n_sSNV}sSNV·c=${r.n_clusters}·${r.haplotypes}·${r.cn}·TP${r.tp}/FP${r.fp}${r.ambig_nodes>0?'·⚠序未定':''}</span></div>`).join('')+(f.length>700?`<div class="note" style="padding:8px">...前 700（共 ${f.length}）</div>`:'');
 el('list').querySelectorAll('.row').forEach(x=>x.onclick=()=>show(+x.dataset.i,x));
}
function show(i,row){el('list').querySelectorAll('.row').forEach(x=>x.classList.remove('sel'));if(row)row.classList.add('sel');let r=det[i];
 let popcount=r.populations;
 let np=r.node_paths||{};
 let cf=fourGamete(r.populations);
 let cand=enumCandidates(r);
 let ctr=(D.candtrees||{})[r.region];
 let pt=Object.entries(r.populations).sort((a,b)=>b[1]-a[1]).map(([g,c])=>{let tot=Object.values(r.populations).reduce((a,b)=>a+b,0);return `<tr><td class="mono" style="color:#1971c2;font-weight:600">${np[g]||(g.includes('A')?'—(未定)':'germline')}</td><td class="mono">${g}</td><td>${sLabels(g)}</td><td>${c}</td><td>${(100*c/tot).toFixed(0)}%</td></tr>`}).join('');
 el('detail').innerHTML=`<h3>${r.region} <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type}</span> <span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span></h3>
  <div class="kv"><div class="b">${r.n_sSNV} sSNV</div><div class="b">span ${(r.span/1000).toFixed(1)}kb</div><div class="b">c=${r.n_clusters} 群</div><div class="b">HP: ${r.haplotypes}</div><div class="b">CN: ${r.cn}</div><div class="b">TP ${r.tp} / FP ${r.fp}</div><div class="b">${r.determinacy}</div>${r.drop_noise_frac>0?`<div class="b">噪聲過濾 ${(r.drop_noise_frac*100).toFixed(0)}%</div>`:''}${r.ambig_nodes>0?`<div class="b" style="background:#fff3bf">⚠ 順序未定 ${r.ambig_nodes}(缺中間群)</div>`:''}${r.truncated?`<div class="b" style="background:#ffe3e3;color:#c92a2a" title="genotype 向量截斷在 8 位點(上游 GCAP=8);此區 ambig/四配子/機率偵測不完整,成環可能為截斷假象">⚠ 截斷 n_sSNV>8(偵測不完整)</div>`:''}</div>
  ${r.undefined?`<div style="background:#ffe3e3;border:1px solid #ffc9c9;border-radius:6px;padding:8px;margin:6px 0"><b>⚠ 此區有無法定義的分支（順序未定/不相容）</b>→ 下方標籤為可能位置。<br>🔴 曾標『需甲基輔助確認』,但 06-28 normal cis-control pilot 裁決:此類區甲基<b>乾淨可用≈0</b>(SAME-HP 在同一 germline HP 內分化、normal 無對應 within-HP 軸=結構性無解)→ 需 single-cell/multi-region 或加深覆蓋,<b>甲基無法解鎖此區</b>。</div>`:''}
  ${cf.length?`<div style="background:#fff0f6;border:1px solid #f783ac;border-radius:6px;padding:9px;margin:6px 0"><b>⚠ 四配子違反（incompatible）→ 無法成單一樹</b>　錨點 <b>RR=germline</b>（normal 確認 REF）→ <b>AA=雙突變（最遠）</b>；RA／AR 兩單突變並存＝累積順序未定（AA 由哪個衍生？）<table style="margin-top:5px"><tr><th>衝突對</th><th>RR<br>germ根</th><th>RA<br>僅後者</th><th>AR<br>僅前者</th><th>AA<br>最遠</th><th>讀數弱提示</th></tr>${cf.map(c=>`<tr><td class="mono"><b>${c.pair}</b></td><td>${c.g.RR}</td><td>${c.g.RA}</td><td>${c.g.AR}</td><td>${c.g.AA}</td><td class="note">${c.g.AR>c.g.RA?'S'+c.i+' 單突變較多':c.g.RA>c.g.AR?'S'+c.j+' 單突變較多':'兩單突變相當'}（弱·非定論）</td></tr>`).join('')}</table><div class="note" style="margin-top:4px">下方樹為「丟掉成環邊後的近似結構」，僅參考；真實關係非單一樹。</div><div style="margin-top:5px;padding:6px 9px;border-radius:5px;font-size:11.5px;background:${r.n_roots>=2?'#e7f5ff':'#fff5f5'};border:1px solid ${r.n_roots>=2?'#74c0fc':'#ffc9c9'}">🧬 <b>此區甲基能否分辨 AA 靠 RA／AR：</b>${r.n_roots>=2?`<b style="color:#1971c2">CROSS-HP（此區跨 H1/H2，${r.n_roots} 根）</b> → 橫跨兩單倍型的衝突對屬 <b>allelic</b>（兩突變在不同染色體），有<b>獨立 germline-ASM 甲基軸</b>，甲基可給弱獨立訊號（06-28：cross-HP ~35% 可控）；但同一 HP 內的對仍 cis-confounded。`:`<b style="color:#c92a2a">SAME-HP（此區單一 germline HP：${r.haplotypes}）</b> → 衝突對皆在同單倍型，甲基隨 genotype 在 cis 共變（cis-ASM）＝<b>結構性無法解此衝突（double-dip）</b>；normal 無對應 within-HP 軸可扣。順序只能靠<b>讀數/VAF 弱先驗</b> + <b>single-cell／multi-region 確認</b>（06-28 cis-control 裁決，L2）。`}</div></div>`:''}
  ${r.n_roots>=2?`<div style="background:#fff4e6;border:1px solid #ffd8a8;border-radius:6px;padding:8px"><b>⚠ 此區跨 H1/H2（${r.n_roots} 棵樹）→ 預設顯示分開的兩棵 HP 樹（正確）：</b>${posTree(r)}</div><details style="margin-top:6px"><summary style="cursor:pointer;color:#868e96;font-size:11.5px">▶ 混合 genotype-向量樹（跨 HP 混合，僅參考）</summary>${tree(r.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes)}</details>`:`<b>克隆樹（germline→…；節點=lineage標籤·S-mut-set·reads·%；座標=向量）</b>${tree(r.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths,r.ambig_nodes)}`}
  ${(ctr&&ctr.candidate_set&&ctr.candidate_set.length)?`<div style="background:#f8f0fc;border:1px solid #d0bfff;border-radius:6px;padding:9px;margin:8px 0"><b>🔀 替代整樹候選（左右滑動完整樹）</b> <span class="note">此區 <b>${ctr.n_candidates}</b> 棵相容候選樹(缺中間群→插虛擬中間節點);${ctr.honest_note}。<b>非在給答案</b>。</span><div id="cttreebox" style="margin-top:6px">${candTreeCard(ctr,0,r)}</div></div>`:(cand?`<div style="background:#f8f0fc;border:1px solid #d0bfff;border-radius:6px;padding:9px;margin:8px 0"><b>🔀 此區某群位置未定（左右滑動看可能排列；非在給答案）</b> <span class="note">缺中間群→中間群未觀察到，下列累積序共 ${cand.trueCount}${cand.bigNode?'+':''} 種、<b>等機率</b>（讀數無法分；甲基判定見上方）${cand.cands.length<cand.trueCount?('，顯示前 '+cand.cands.length):''}</span><div id="candbox" style="margin-top:6px">${candCard(cand.cands,0,r)}</div></div>`:'')}
  <div class="note">S1..S${r.n_sSNV}=區內排序 sSNV；直系=往下、姊妹=同層分叉；germline 根標 reads·%。tree_shape(pairwise)=${r.tree_shape}。genome_ctx 為近似(±3Mb)。</div>
  <b>細胞群(lineage 標籤 → 向量 → S 突變 → reads → 佔比)</b><table><tr><th>lineage</th><th>向量</th><th>突變(S)</th><th>reads</th><th>佔比</th></tr>${pt}</table>${geneBlock(r.region)}`;
 window.__cand=cand?{cands:cand.cands,idx:0,r:r}:null;
 window.__ctr=(ctr&&ctr.candidate_set&&ctr.candidate_set.length)?{cs:ctr,idx:0,r:r}:null;
}
function geneBlock(region){
  let g=(D.gene||{})[region]; if(!g) return '';
  let cg=Object.entries(g.cancer_genes||{}).map(([n,o])=>`<span style="background:#ffe3e3;border-radius:3px;padding:1px 4px">${n}${o.role?'('+o.role+')':''}${o.tier?' T'+o.tier:''}</span>`).join(' ');
  let dr=Object.entries(g.druggable_genes||{}).map(([n,ds])=>`<span style="background:#d3f9d8;border-radius:3px;padding:1px 4px" title="${ds.join(', ')}">${n}💊</span>`).join(' ');
  let names=(g.protein_coding||g.genes||[]).slice(0,12);
  return `<div style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;padding:8px;margin-top:8px">
    <b>🧬 基因註釋(GENCODE+DGIdb${cg?'+COSMIC':''})</b>
    <div class="note">基因(${g.n_genes}): ${names.join(', ')||'(無 protein-coding)'}${g.has_promoter?` · <b style="color:#1971c2">含啟動子</b>(${(g.promoter_genes||[]).slice(0,5).join(',')})`:''}</div>
    ${cg?`<div style="margin-top:4px">🔴 癌症基因(COSMIC): ${cg}</div>`:''}
    ${dr?`<div style="margin-top:4px">💊 可用藥(DGIdb): ${dr}</div>`:'<div class="note" style="margin-top:4px">此區無 DGIdb 可用藥基因</div>'}
  </div>`;
}
['f_chr','f_minc','f_tpfp','f_loh','f_undef','f_q','f_sort','f_sortdir'].forEach(id=>{let e=el(id);if(e){e.oninput=render;e.onchange=render}});
render();
// ===== 確認佇列(評分 + 左右判讀) =====
const SC=D.scoring;
el('scoresum').innerHTML=`需確認 <b>${SC.summary.n_need_confirm}</b>/${SC.summary.n_total} 區 · 評分桶 ${JSON.stringify(SC.summary.score_buckets)} · situation ${JSON.stringify(SC.summary.situation_dist)} · <span title="06-28 cis-control 已否決:乾淨可用≈0">曾標需甲基 ${SC.summary.needs_methyl_n}(已否決·非真可用)</span> · 公式: ${SC.summary.score_formula}`;
let Q=SC.queue;
el('q_sit').innerHTML='<option value="">全</option>';[...new Set(Q.map(q=>q.situation))].sort().forEach(s=>{let o=document.createElement('option');o.value=s;o.textContent=s;el('q_sit').appendChild(o)});
const QSORT={score:(a,b)=>a.confidence_score-b.confidence_score,scoreD:(a,b)=>b.confidence_score-a.confidence_score,coord:(a,b)=>a.chrom.localeCompare(b.chrom,undefined,{numeric:true})||a.start-b.start};
const jkey=r=>'topo_judge_'+r;
window.setJ=(r,v)=>{let cur=localStorage.getItem(jkey(r));localStorage.setItem(jkey(r),cur==v?'':v);renderQ()};
const scolor=s=>s>=80?'#2b8a3e':s>=60?'#1971c2':s>=40?'#e8590c':'#c92a2a';
function renderQ(){let sit=el('q_sit').value,mo=el('q_methyl').checked,so=el('q_sort').value;
 let f=Q.filter(q=>(!sit||q.situation==sit)&&(!mo||q.needs_methyl));f.sort(QSORT[so]||QSORT.score);
 el('qcnt').textContent=f.length+' 區';
 el('queue').innerHTML=f.slice(0,500).map(q=>{let j=localStorage.getItem(jkey(q.region))||'';
  return `<div class="row" style="display:flex;gap:7px;align-items:center;flex-wrap:wrap">
   <span style="width:140px"><b>${q.region}</b></span>
   <span style="width:58px;color:${scolor(q.confidence_score)};font-weight:700" title="confidence 0-100">▮${q.confidence_score}</span>
   <span class="tag ctx_${q.genome_ctx}">${q.genome_ctx}</span>
   <span style="width:112px;font-size:11px">${q.situation}</span>
   <span style="width:96px;font-size:10.5px" class="note">候選 ${q.n_candidates}${q.parsimony_first_rank_prob!=null?` · P1=${q.parsimony_first_rank_prob}`:''}</span>
   <span style="flex:1;min-width:170px;font-size:10.5px" class="note">${q.resolution_path}</span>
   <span style="width:100%;font-size:10px;color:#a33" class="note">🔎 為何: ${q.why_conflict||''}${q.truncated?' ⚠截斷':''}</span>
   <span style="width:100%;font-size:10px;color:#268" class="note">🧬 甲基: ${q.methyl_applicability||''}</span>
   <span style="white-space:nowrap">
     <button onclick="setJ('${q.region}','agree')" style="font-size:11px;background:${j=='agree'?'#d3f9d8':'#fff'}">✓同意rank1</button>
     <button onclick="setJ('${q.region}','alt')" style="font-size:11px;background:${j=='alt'?'#fff3bf':'#fff'}">⇄偏好其他</button>
     <button onclick="setJ('${q.region}','more')" style="font-size:11px;background:${j=='more'?'#ffe3e3':'#fff'}">?需更多資訊</button>
   </span></div>`}).join('')+(f.length>500?`<div class="note" style="padding:8px">...前 500（共 ${f.length}，可篩選縮小）</div>`:'');
}
['q_sort','q_sit','q_methyl'].forEach(id=>{el(id).onchange=renderQ;el(id).oninput=renderQ});
el('q_exp').onclick=()=>{let j={};Q.forEach(q=>{let v=localStorage.getItem(jkey(q.region));if(v)j[q.region]={judgment:v,score:q.confidence_score,situation:q.situation}});
 let b=new Blob([JSON.stringify({n:Object.keys(j).length,judgments:j},null,1)],{type:'application/json'});
 let a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='topology_judgments.json';a.click()};
renderQ();
}  // end bootWS — 每分頁切換時用新樣本資料重跑(function scope,無重宣告衝突)
function selectSample(s){
  window.__DATA__ = window.__SAMPLES__[s];
  document.querySelectorAll('.stab').forEach(t=>t.classList.toggle('active', t.dataset.s===s));
  try{ bootWS(); }catch(e){ document.getElementById('detail').innerHTML='<b style="color:#c00">分頁載入錯誤: '+e.message+'</b>'; }
}
(function(){
  let names=Object.keys(window.__SAMPLES__||{});
  let bar=document.getElementById('sampletabs');
  if(bar){ bar.innerHTML = names.map(n=>`<button class="stab" data-s="${n}" onclick="selectSample('${n}')">${n}</button>`).join(''); }
  if(names.length) selectSample(names[0]);
})();
"""

PROVENANCE_FOOTER = ('<p class="note" style="margin-top:8px;color:#888">'
                     'build_branch: research/subclonal-reconstruction-202606 · '
                     '資料 topology_per_region.json（凍結 @ feat/summary-nreadsvalid@5308d9e）· '
                     '姊妹編號 = 子樹總 read 數遞減（?-1=該 lineage 分支佔比較大，含所有子孫；?-2=較小）· '
                     '甲基 = bounded-auxiliary（見 20260628_cis_control_scope_pilot_verdict_01.md）</p>')

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>多樣本克隆樹拓樸工作站 — {len(SAMPLE_NAMES)} ONT 樣本 sSNV 重建</title><style>{CSS}
.tabs{{display:flex;gap:4px;flex-wrap:wrap;margin:8px 0;border-bottom:2px solid #dee2e6;padding-bottom:0}}
.stab{{padding:7px 14px;border:1px solid #dee2e6;border-bottom:none;border-radius:7px 7px 0 0;background:#f1f3f5;cursor:pointer;font-size:13px;font-weight:600;color:#495057}}
.stab.active{{background:#1971c2;color:#fff}}
</style></head><body><div class="wrap">
<h1>多樣本克隆樹拓樸互動工作站（cluster-first + S/r/m 標籤 + 基因註釋）</h1>
<p class="sub">{len(SAMPLE_NAMES)} ONT 樣本（{", ".join(SAMPLE_NAMES)}）· 每區 genotype 向量→拓樸(perfect-phylogeny+噪聲過濾) · 分頁切換樣本 · S=sSNV/r=read群/m=甲基位點 · 數字由 JSON 注入</p>
<div id="sampletabs" class="tabs"></div>
<div class="zone">📊 整體觀察區（隨樣本變）：全 sSNV 宇宙帳本 + 拓樸/群數/determinacy/HP 統計（圓餅+長條）</div>
<div id="universe"></div>
{GLOSSARY_HTML}
<div class="stats">
<div class="scard"><h4>拓樸型態<span class="more">▸ 點看細節</span></h4><div id="s_topo"></div></div><div class="scard"><h4>群數 c<span class="more">▸ 點看細節</span></h4><div id="s_clust"></div></div>
<div class="scard"><h4>determinacy<span class="more">▸ 點看細節</span></h4><div id="s_det"></div></div><div class="scard"><h4>HP 根數<span class="more">▸ 點看細節</span></h4><div id="s_root"></div></div>
<div class="scard" style="cursor:default"><h4>建樹位點分布<span class="more">區×sSNV</span></h4><div id="s_nsnv"></div></div>
</div>
<div class="smbg" id="statmodal" onclick="if(event.target===this)closeStatModal()"><div class="smbox"><button class="x" onclick="closeStatModal()">×</button><div id="statmodal_body"></div></div></div>
<details class="c17" id="chr17wrap"><summary>📚 克隆樹判讀圖鑑（固定教學・範例取自 HCC1395，與上方分頁樣本無關；點開）</summary><div id="teachbody"></div></details>
<div id="unlocatable"></div>
<div class="zone">🔬 樣本各區域檢視（篩選 → 點選一個區 → 看克隆樹 / 四配子衝突 / 逐區甲基判定）</div>
<div class="ctrl">
chr<select id="f_chr"><option value="">全</option></select>
排序<select id="f_sort"><option value="coord">座標</option><option value="nsnv">複雜度(sSNV)</option><option value="nclust">群數</option><option value="region">region名</option></select><select id="f_sortdir" title="正序/反序"><option value="">↓預設</option><option value="rev">↑反序</option></select>
最少群數<input id="f_minc" type="number" value="0" min="0" max="6" style="width:52px">
TP/FP<select id="f_tpfp"><option value="all">全部</option><option value="tp">只含TP</option><option value="fp">只含FP</option><option value="both">同時TP&amp;FP</option></select>
<label title="LOH 區且單一 HP 標籤(LOH-unmask 觀察)"><input id="f_loh" type="checkbox">僅 LOH 單HP</label>
<label title="分支順序未定/不相容;曾標『需甲基輔助』,06-28 cis-control 已否決(乾淨可用≈0)"><input id="f_undef" type="checkbox">僅無法定義(曾標需甲基·已否決)</label>
搜尋<input id="f_q" placeholder="chr17:" style="width:120px">
<span id="cnt" class="note"></span>
</div>
<div class="zone" style="margin-top:14px;font-size:12.5px;background:#fff9db;color:#b08900;border-left-color:#ffd43b">☑ 勾選要觀察的狀況（三組可複選；不勾＝該組全納入）</div>
<div class="facets">
 <div class="fcard"><div class="fh">🌳 拓樸型態<span class="fhd">分子累積形狀</span></div><div class="chips" id="cb_topology_type"></div></div>
 <div class="fcard"><div class="fh">🎯 determinacy<span class="fhd">能否唯一辨識</span></div><div class="chips" id="cb_determinacy"></div></div>
 <div class="fcard"><div class="fh">📍 基因體位置<span class="fhd">偽影風險</span></div><div class="chips" id="cb_genome_ctx"></div></div>
</div>
<div class="main"><div class="list" id="list"></div><div class="detail" id="detail"><div class="note">← 左側點選一個區查看克隆樹（或點上方 chr17 worked example）</div></div></div>

<h3 style="margin-top:20px">✓ 候選評分確認佇列（左右選項判讀 + 觀察評分；存瀏覽器 localStorage、可匯出）</h3>
<div id="scoresum" class="note"></div>
<div class="ctrl">
排序<select id="q_sort"><option value="score">評分(低→高,最需關注)</option><option value="scoreD">評分(高→低)</option><option value="coord">座標</option></select>
situation<select id="q_sit"><option value="">全</option></select>
<label title="06-28 cis-control 裁決:這些區甲基乾淨可用≈0,非真能用甲基解"><input id="q_methyl" type="checkbox">曾標需甲基(已否決)</label>
<button id="q_exp">匯出判讀 JSON</button><span id="qcnt" class="note"></span>
</div>
<div class="list" id="queue" style="max-height:62vh"></div>
<p class="note" style="margin-top:12px">⚠ 證據層級：A_determined=單分子向量唯一可辨識(≠對 single-cell 驗證為真)；A_ambiguous=缺中間群順序未定；B_pairwise=拼接非單分子整樹；C_underdetermined=多樹相容。TP/FP=SEQC2 僅觀察不進前處理。genome_ctx 為近似(±3Mb)。甲基不參與拓樸裁決(cis-confounded;06-28 cis-control 已測→bounded-auxiliary,非 resolver)。⭐3 單樣本·regional(≤read-span)非 genome-wide tree·分子共現≠single-cell。</p>
{PROVENANCE_FOOTER}
</div>
<script>window.__SAMPLES__={SAMPLES_JSON};window.__CHR17TREE__={CHR17TREE_JSON};</script><script>{JS}</script></body></html>"""
with open(MULTI_OUT, "w", encoding="utf-8") as f: f.write(HTML)
print(f"OK wrote {MULTI_OUT} ({len(HTML):,} bytes; samples {SAMPLE_NAMES})")
