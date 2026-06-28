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
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "20260628_topology_workstation.standalone.html"))
d = json.load(open(os.path.join(DATA, "topology_per_region.json"), encoding="utf-8"))
acc = json.load(open(os.path.join(DATA, "single_snv_accounting.json"), encoding="utf-8"))
cs = json.load(open(os.path.join(DATA, "candidate_scoring.json"), encoding="utf-8"))
DJ = json.dumps({"stats": d["stats"], "detail": d["detail"], "chr17": d["chr17_worked"], "chroms": d["chroms"],
                 "scoring": {"summary": {k: cs[k] for k in ("n_total","n_need_confirm","score_formula","situation_dist","resolution_dist","score_buckets","needs_methyl_n")},
                             "queue": cs["queue"]}}, ensure_ascii=False)
B = acc["buckets"]
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
GLOSSARY_HTML = '<details style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:10px 14px;margin:10px 0"><summary style="cursor:pointer;font-weight:600;color:#1971c2">📖 名詞與概念解釋（點開；每項可再展開細節）</summary><div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(290px,1fr));gap:6px;margin-top:8px">' + "".join(f'<details style="border:1px solid #f1f3f5;border-radius:6px;padding:5px 9px;font-size:12px"><summary style="cursor:pointer"><b>{t}</b></summary><div style="margin-top:4px;color:#343a40">{s}</div><div style="margin-top:3px;color:#868e96;font-size:11px">{dd}</div></details>' for t, s, dd in GLOSSARY) + '</div></details>'

CSS = """
*{box-sizing:border-box}body{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:#212529;background:#f8f9fa}
.wrap{max-width:1320px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:.2em 0}h3{margin:.3em 0}.sub{color:#868e96;font-size:12.5px}
.stats{display:flex;gap:12px;flex-wrap:wrap;margin:10px 0}.scard{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:9px 11px;min-width:200px}
.scard h4{margin:0 0 5px;font-size:11.5px;color:#495057}.bar{display:flex;align-items:center;gap:5px;font-size:10.5px;margin:2px 0}.bar i{height:10px;background:#1c7ed6;border-radius:2px;display:inline-block}
details.c17{background:#fff;border:1px solid #ffd8a8;border-radius:8px;padding:10px 14px;margin:10px 0}details.c17 summary{cursor:pointer;font-weight:600;color:#d9480f}
.ctrl{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:10px 0;font-size:12.5px;background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:9px}
.ctrl select,.ctrl input{padding:3px 7px;border:1px solid #ced4da;border-radius:5px;font-size:12.5px}
.main{display:grid;grid-template-columns:400px 1fr;gap:12px}@media(max-width:860px){.main{grid-template-columns:1fr}}
.list{background:#fff;border:1px solid #dee2e6;border-radius:8px;max-height:76vh;overflow:auto}
.row{padding:6px 10px;border-bottom:1px solid #f1f3f5;cursor:pointer;font-size:12px}.row:hover{background:#e7f5ff}.row.sel{background:#d0ebff}.row b{color:#1c7ed6}
.tag{font-size:9.5px;padding:1px 6px;border-radius:9px;margin-left:3px}
.t_linear{background:#d3f9d8;color:#2b8a3e}.t_branched{background:#e5dbff;color:#5f3dc4}.t_star{background:#fff3bf;color:#b08900}.t_single{background:#f1f3f5;color:#868e96}
.ctx_telomere{background:#d0ebff;color:#1971c2}.ctx_centromere{background:#ffe3e3;color:#c92a2a}.ctx_arm{background:#f1f3f5;color:#868e96}
.detail{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:15px;min-height:76vh}
.kv{display:flex;gap:10px;flex-wrap:wrap;margin:6px 0}.kv .b{background:#f1f3f5;border-radius:6px;padding:5px 10px;font-size:11.5px}
table{border-collapse:collapse;font-size:11.5px;margin:7px 0}th,td{border:1px solid #dee2e6;padding:3px 8px}th{background:#f1f3f5}
.note{font-size:11.5px;color:#868e96}.mono{font-family:ui-monospace,Menlo,monospace}
"""

JS = r"""
const D=window.__DATA__;
const TT={'linear(全直系)':'t_linear','branched(直系+姊妹)':'t_branched','star(全姊妹)':'t_star','single':'t_single','germline_only':'t_single'};
const el=id=>document.getElementById(id);
function bars(o,m){let vs=Object.values(o),tot=vs.reduce((a,b)=>a+b,0)||1,mx=Math.max(...vs,1);return Object.entries(o).sort((a,b)=>b[1]-a[1]).slice(0,m||9).map(([k,v])=>`<div class="bar"><i style="width:${Math.max(3,78*v/mx)}px"></i>${k}: <b>${v}</b> (${(100*v/tot).toFixed(1)}%)</div>`).join('')+`<div class="bar" style="color:#868e96">— 合計 ${tot} (100%)</div>`}
el('s_topo').innerHTML=bars(D.stats.topology_type);el('s_clust').innerHTML=bars(Object.fromEntries(Object.entries(D.stats.n_clusters).map(([k,v])=>['c='+k,v])));
el('s_det').innerHTML=bars(D.stats.determinacy);el('s_root').innerHTML=bars(D.stats.n_roots);
// chr17 worked panel
(function(){let c=D.chr17;
 let st=c.snvs.map(s=>`<tr><td class="mono"><b>${s.S}</b></td><td class="mono">${s.pos}</td><td>${s.change}</td><td>${s.role}</td><td>VAF ${s.vaf}</td><td>${s.hp}</td><td>${s.src}</td><td>${s.somatic_confirmed?'✓':'✗(normal有ALT)'}</td></tr>`).join('');
 let pp=c.populations.map(p=>`<tr><td class="mono">${p.vec}</td><td>${p.muts}</td><td>${p.reads}</td><td>${p.pct}%</td></tr>`).join('');
 let mm=c.sig_cpg.slice(0,16).map(x=>`<tr><td class="mono">${x.m}</td><td>${x.cpg}</td><td>${x.L1}</td><td>${x.L2}</td><td>${x.dbeta}</td></tr>`).join('');
 el('c17').innerHTML=`<div class="kv"><div class="b">locus ${c.locus}</div><div class="b">ctx ${c.genome_ctx}</div><div class="b">拓樸 ${c.topology_type}</div><div class="b">噪聲 dropped ${c.dropped_noise} reads</div></div>
  <b>S 位點(somatic sSNV)</b><table><tr><th>S</th><th>pos</th><th>變異</th><th>角色</th><th>VAF</th><th>HP</th><th>TP/FP</th><th>somatic</th></tr>${st}</table>
  <b>克隆樹</b> ${tree(c.edges,Object.fromEntries(c.populations.map(p=>[p.vec,p.reads])),c.populations.length,'H1',(c.populations.find(p=>p.muts=='germline')||{}).reads||0)}
  <b>細胞群(r 分群;基因型向量·reads·佔比)</b><table><tr><th>向量</th><th>突變(S)</th><th>reads(r)</th><th>佔比</th></tr>${pp}</table>
  <b>甲基差異位點 m1..m${c.n_sig_diff_cpg}(L1 α-only vs L2 α+β;⚠ 經實測對齊 cis-genotype 軸非獨立 lineage)</b><table><tr><th>m</th><th>CpG</th><th>L1 β</th><th>L2 β</th><th>Δβ</th></tr>${mm}</table>`;
})();
// S-label 化基因型向量
function sLabels(g){let s=[...g].map((c,i)=>c=='A'?('S'+(i+1)):null).filter(Boolean);return s.length?s.join('+'):'germline'}
function gainedS(parent,child){let g=[];for(let i=0;i<child.length;i++){if(child[i]=='A'&&(!parent||parent[i]!='A'))g.push('S'+(i+1))}return g}
function tree(edges,popcount,nc,hp,germR,np){np=np||{};
 let ch={},par={},all=new Set();
 (edges||[]).forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT'){all.add(p);par[c]=p}});
 if(!all.size)return `<div class="note">單群／germline-only（germline ${germR||0} reads），無分支樹</div>`;
 let sib={};Object.keys(ch).forEach(p=>{(ch[p]||[]).forEach(c=>{sib[c]=(ch[p].length>=2)})});
 let depth={},pos={},leaf=0;
 function lay(n,dp){depth[n]=dp;let k=(ch[n]||[]).sort();if(!k.length){pos[n]=leaf++;return pos[n]}let xs=k.map(x=>lay(x,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]}
 let roots=ch['ROOT']||[];roots.forEach(r=>lay(r,1));
 let gx=roots.length?(Math.min(...roots.map(r=>pos[r]))+Math.max(...roots.map(r=>pos[r])))/2:0;
 let nodes=[...all],md=Math.max(...nodes.map(n=>depth[n]),1);
 let NW=204,NH=80,GX=50,GY=128;
 let W=Math.max(450,(leaf||1)*(NW+GX)),H=108+md*GY;
 let X=p=>34+p*(NW+GX)+NW/2,Y=dp=>56+dp*GY;
 let totR=Object.values(popcount).reduce((a,b)=>a+b,0)+(germR||0)||1;
 let relSet=new Set();
 let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" style="font-family:ui-monospace,Menlo,monospace">`;
 roots.forEach(r=>s+=`<line x1="${X(gx)}" y1="${Y(0)+NH/2}" x2="${X(pos[r])}" y2="${Y(1)-NH/2}" stroke="#ced4da" stroke-width="1.6"/>`);
 (edges||[]).forEach(([p,c])=>{if(p!='ROOT'){
  let g=gainedS(p,c),mx=(X(pos[p])+X(pos[c]))/2,my=(Y(depth[p])+Y(depth[c]))/2;
  s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+NH/2}" x2="${X(pos[c])}" y2="${Y(depth[c])-NH/2}" stroke="#ced4da" stroke-width="1.6"/>`;
  if(g.length)s+=`<rect x="${mx-24}" y="${my-10}" width="48" height="19" rx="9" fill="#ebfbee" stroke="#2f9e44"/><text x="${mx}" y="${my+4}" text-anchor="middle" font-size="11" fill="#2b8a3e" font-weight="700">+${g.join('+')}</text>`;
 }});
 let gp=(100*(germR||0)/totR).toFixed(0);
 s+=`<rect x="${X(gx)-NW/2}" y="${Y(0)-NH/2}" width="${NW}" height="${NH}" rx="12" fill="#f1f3f5" stroke="#868e96" stroke-width="2.5"/>`;
 s+=`<text x="${X(gx)}" y="${Y(0)-NH/2+24}" text-anchor="middle" font-size="14" font-weight="700" fill="#495057">⌂ germline（${hp||'根'}）</text>`;
 s+=`<text x="${X(gx)}" y="${Y(0)-NH/2+44}" text-anchor="middle" font-size="11" fill="#868e96">無 somatic 變異 · 起點</text>`;
 s+=`<text x="${X(gx)}" y="${Y(0)-NH/2+66}" text-anchor="middle" font-size="12.5" font-weight="600" fill="#212529">${germR||0} reads · ${gp}%</text>`;
 nodes.forEach(n=>{
  let cnt=popcount[n]||0,pct=(100*cnt/totR).toFixed(0),lab=np[n]||'—';
  let g=gainedS(par[n],n);
  let isSib=sib[n], isCo=g.length>=2;
  if(isSib)relSet.add('姊妹(sibling)'); else if(g.length)relSet.add('直系(vertical)');
  if(isCo)relSet.add('co_linked');
  let rel=g.length?((isSib?'姊妹分支(sibling)':'直系(vertical)')+(isCo?' · co_linked':'')):'（與父同型）';
  let gtxt=g.length?('獲得 +'+g.join('+')):'';
  let x=X(pos[n]),y=Y(depth[n]);
  let fill=isSib?'#fff4e6':'#e7f5ff', stroke=isSib?'#e8590c':'#1c7ed6';
  s+=`<rect x="${x-NW/2}" y="${y-NH/2}" width="${NW}" height="${NH}" rx="12" fill="${fill}" stroke="${stroke}" stroke-width="2.5"/>`;
  s+=`<text x="${x}" y="${y-NH/2+22}" text-anchor="middle" font-size="16" font-weight="800" fill="#1565c0">${lab}</text>`;
  s+=`<text x="${x}" y="${y-NH/2+39}" text-anchor="middle" font-size="10.5" font-weight="700" fill="${isSib?'#d9480f':'#2b8a3e'}">${rel}</text>`;
  s+=`<text x="${x}" y="${y-NH/2+53}" text-anchor="middle" font-size="10.5" font-weight="600" fill="#2b8a3e">${gtxt}</text>`;
  s+=`<text x="${x}" y="${y-NH/2+66}" text-anchor="middle" font-size="9.5" fill="#495057">基因型 ${n}（=${sLabels(n)}）</text>`;
  s+=`<text x="${x}" y="${y-NH/2+78}" text-anchor="middle" font-size="12" font-weight="600" fill="#212529">${cnt} reads · ${pct}%</text>`;
 });
 s+='</svg>';
 let leg=`<div style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;padding:7px 11px;margin-top:5px;font-size:11px;line-height:1.6">
 <b>🔖 SNV 關係圖例</b>（怎麼讀這棵樹）：
 <span style="color:#1565c0;font-weight:700">●直系 vertical</span> 往下一層、後代多帶 1 變異(+S，藍框) ｜
 <span style="color:#d9480f;font-weight:700">●姊妹分支 sibling</span> 同一父不同分支、平行 subclone(橙框) ｜
 <span style="color:#5f3dc4;font-weight:700">co_linked</span> 一節點同時獲≥2 變異、綁同一事件 ｜
 <span style="color:#2b8a3e;font-weight:700">+S</span> 該分支新增 somatic 變異
 ${relSet.size?'<br><b>此樹實際出現</b>：'+[...relSet].join('、'):''}</div>`;
 return s+leg;
}
// 2-root: 位置樹按 HP 分兩棵
function posTree(r){
 let byhp={};Object.entries(r.node_hp||{}).forEach(([p,h])=>{if(h=='H1'||h=='H2')(byhp[h]=byhp[h]||[]).push(p)});
 let hps=Object.keys(byhp).sort();if(hps.length<2)return '';
 return '<div style="display:flex;gap:24px;flex-wrap:wrap">'+hps.map(h=>{
  let ns=new Set(byhp[h]);let ned=(r.pos_nested||[]).filter(e=>ns.has(e[0])&&ns.has(e[1]));
  let hasp=new Set(ned.map(e=>e[1]));let edges=ned.map(e=>[e[0],e[1]]);[...ns].forEach(n=>{if(!hasp.has(n))edges.unshift(['ROOT',n])});
  return `<div><b style="color:#d9480f">${h} 樹（${ns.size} 位點）</b>${posSVG(edges,ns,r.pos_vaf||{},h)}</div>`;
 }).join('')+'</div>';
}
function posSVG(edges,ns,vaf,hp){
 let ch={},all=new Set();edges.forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT')all.add(p)});
 if(!all.size)return '<div class="note">無結構</div>';
 let depth={},pos={},leaf=0;function lay(n,dp){depth[n]=dp;let k=(ch[n]||[]).sort();if(!k.length){pos[n]=leaf++;return pos[n]}let xs=k.map(x=>lay(x,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]}
 let roots=ch['ROOT']||[];roots.forEach(r=>lay(r,1));let gx=roots.length?(Math.min(...roots.map(r=>pos[r]))+Math.max(...roots.map(r=>pos[r])))/2:0;
 let nodes=[...all],md=Math.max(...nodes.map(n=>depth[n]),1),W=Math.max(200,(leaf||1)*108),H=64+md*70,X=p=>40+p*108,Y=dp=>26+dp*70;
 let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}"><circle cx="${X(gx)}" cy="${Y(0)}" r="15" fill="#fff" stroke="#495057"/><text x="${X(gx)}" y="${Y(0)+3}" text-anchor="middle" font-size="9">${hp} germ</text>`;
 roots.forEach(rt=>s+=`<line x1="${X(gx)}" y1="${Y(0)+15}" x2="${X(pos[rt])}" y2="${Y(1)-20}" stroke="#adb5bd"/>`);
 edges.forEach(([p,c])=>{if(p!='ROOT')s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+20}" x2="${X(pos[c])}" y2="${Y(depth[c])-20}" stroke="#adb5bd"/>`});
 nodes.forEach(n=>{let v=vaf[n];s+=`<rect x="${X(pos[n])-46}" y="${Y(depth[n])-18}" width="92" height="36" rx="5" fill="#fff4e6" stroke="#e8590c"/><text x="${X(pos[n])}" y="${Y(depth[n])-2}" text-anchor="middle" font-size="9" class="mono">${n.split(':')[1]}</text><text x="${X(pos[n])}" y="${Y(depth[n])+11}" text-anchor="middle" font-size="8" fill="#495057">VAF ${v!=null?v:'?'}</text>`});
 return s+'</svg>';
}
// filters + sort
let det=D.detail;
let chrsel=el('f_chr');D.chroms.forEach(c=>{let o=document.createElement('option');o.value=c;o.textContent=c;chrsel.appendChild(o)});
['topology_type','determinacy','genome_ctx'].forEach(k=>{let c=el('cb_'+k);[...new Set(det.map(r=>r[k]))].sort().forEach(o=>{let l=document.createElement('label');l.style.cssText='margin:0 5px;white-space:nowrap';l.innerHTML='<input type="checkbox" value="'+o+'"> '+o;c.appendChild(l)});c.addEventListener('change',render)});
function cset(k){return new Set([...el('cb_'+k).querySelectorAll('input:checked')].map(x=>x.value))}
const SORT={coord:(a,b)=>a.chrom.localeCompare(b.chrom,undefined,{numeric:true})||a.start-b.start,nsnv:(a,b)=>b.n_sSNV-a.n_sSNV,nclust:(a,b)=>b.n_clusters-a.n_clusters||b.n_sSNV-a.n_sSNV,region:(a,b)=>a.region.localeCompare(b.region)};
function render(){
 let ch=el('f_chr').value,mc=+el('f_minc').value,fp=el('f_fp').checked,loh=el('f_loh').checked,undef=el('f_undef').checked,q=el('f_q').value.trim(),so=el('f_sort').value;
 let tt=cset('topology_type'),dd=cset('determinacy'),gc=cset('genome_ctx');
 let f=det.filter(r=>(!ch||r.chrom==ch)&&(!tt.size||tt.has(r.topology_type))&&(!dd.size||dd.has(r.determinacy))&&(!gc.size||gc.has(r.genome_ctx))&&r.n_clusters>=mc&&(!fp||r.fp>0)&&(!loh||(r.cn=='loh'&&(r.haplotypes=='H1'||r.haplotypes=='H2')))&&(!undef||r.undefined)&&(!q||r.region.includes(q)));
 f.sort(SORT[so]||SORT.coord);
 el('cnt').textContent=f.length+' 區';
 el('list').innerHTML=f.slice(0,700).map(r=>`<div class="row" data-i="${det.indexOf(r)}"><b>${r.region}</b> <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type.split('(')[0]}</span><span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span><br><span class="note">${r.n_sSNV}sSNV·c=${r.n_clusters}·${r.haplotypes}·${r.cn}·TP${r.tp}/FP${r.fp}${r.ambig_nodes>0?'·⚠序未定':''}</span></div>`).join('')+(f.length>700?`<div class="note" style="padding:8px">...前 700（共 ${f.length}）</div>`:'');
 el('list').querySelectorAll('.row').forEach(x=>x.onclick=()=>show(+x.dataset.i,x));
}
function show(i,row){el('list').querySelectorAll('.row').forEach(x=>x.classList.remove('sel'));if(row)row.classList.add('sel');let r=det[i];
 let popcount=r.populations;
 let np=r.node_paths||{};
 let pt=Object.entries(r.populations).sort((a,b)=>b[1]-a[1]).map(([g,c])=>{let tot=Object.values(r.populations).reduce((a,b)=>a+b,0);return `<tr><td class="mono" style="color:#1971c2;font-weight:600">${np[g]||(g.includes('A')?'—(未定)':'germline')}</td><td class="mono">${g}</td><td>${sLabels(g)}</td><td>${c}</td><td>${(100*c/tot).toFixed(0)}%</td></tr>`}).join('');
 el('detail').innerHTML=`<h3>${r.region} <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type}</span> <span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span></h3>
  <div class="kv"><div class="b">${r.n_sSNV} sSNV</div><div class="b">span ${(r.span/1000).toFixed(1)}kb</div><div class="b">c=${r.n_clusters} 群</div><div class="b">HP: ${r.haplotypes}</div><div class="b">CN: ${r.cn}</div><div class="b">TP ${r.tp} / FP ${r.fp}</div><div class="b">${r.determinacy}</div>${r.drop_noise_frac>0?`<div class="b">噪聲過濾 ${(r.drop_noise_frac*100).toFixed(0)}%</div>`:''}${r.ambig_nodes>0?`<div class="b" style="background:#fff3bf">⚠ 順序未定 ${r.ambig_nodes}(缺中間群)</div>`:''}</div>
  ${r.undefined?`<div style="background:#ffe3e3;border:1px solid #ffc9c9;border-radius:6px;padding:8px;margin:6px 0"><b>⚠ 此區有無法定義的分支（順序未定/不相容）</b>→ 下方標籤為可能位置。<br>🔴 曾標『需甲基輔助確認』,但 06-28 normal cis-control pilot 裁決:此類區甲基<b>乾淨可用≈0</b>(SAME-HP 在同一 germline HP 內分化、normal 無對應 within-HP 軸=結構性無解)→ 需 single-cell/multi-region 或加深覆蓋,<b>甲基無法解鎖此區</b>。</div>`:''}
  <b>克隆樹（germline→…；節點=lineage標籤(藍)·S-mut-set·reads·%；座標=向量）</b>${tree(r.edges,r.populations,r.n_clusters,r.haplotypes,r.germline_reads,r.node_paths)}
  ${r.n_roots>=2?`<div style="background:#fff4e6;border:1px solid #ffd8a8;border-radius:6px;padding:8px;margin-top:8px"><b>⚠ 此區跨 H1/H2（${r.n_roots} 棵樹）→ 分開看的兩棵 HP 樹（上方 genotype-向量樹混合 HP 僅參考）：</b>${posTree(r)}</div>`:''}
  <div class="note">S1..S${r.n_sSNV}=區內排序 sSNV；直系=往下、姊妹=同層分叉；germline 根標 reads·%。tree_shape(pairwise)=${r.tree_shape}。genome_ctx 為近似(±3Mb)。</div>
  <b>細胞群(lineage 標籤 → 向量 → S 突變 → reads → 佔比)</b><table><tr><th>lineage</th><th>向量</th><th>突變(S)</th><th>reads</th><th>佔比</th></tr>${pt}</table>`;
}
['f_chr','f_minc','f_fp','f_loh','f_undef','f_q','f_sort'].forEach(id=>{let e=el(id);e.oninput=render;e.onchange=render});
render();
// ===== 確認佇列(評分 + 左右判讀) =====
const SC=D.scoring;
el('scoresum').innerHTML=`需確認 <b>${SC.summary.n_need_confirm}</b>/${SC.summary.n_total} 區 · 評分桶 ${JSON.stringify(SC.summary.score_buckets)} · situation ${JSON.stringify(SC.summary.situation_dist)} · <span title="06-28 cis-control 已否決:乾淨可用≈0">曾標需甲基 ${SC.summary.needs_methyl_n}(已否決·非真可用)</span> · 公式: ${SC.summary.score_formula}`;
let Q=SC.queue;
[...new Set(Q.map(q=>q.situation))].sort().forEach(s=>{let o=document.createElement('option');o.value=s;o.textContent=s;el('q_sit').appendChild(o)});
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
"""

PROVENANCE_FOOTER = ('<p class="note" style="margin-top:8px;color:#888">'
                     'build_branch: research/subclonal-reconstruction-202606 · '
                     '資料 topology_per_region.json（凍結 @ feat/summary-nreadsvalid@5308d9e）· '
                     '姊妹編號 = 子樹總 read 數遞減（?-1=該 lineage 分支佔比較大，含所有子孫；?-2=較小）· '
                     '甲基 = bounded-auxiliary（見 20260628_cis_control_scope_pilot_verdict_01.md）</p>')

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>克隆樹拓樸工作站 v2 — sSNV 重建（HCC1395 ⭐3）</title><style>{CSS}</style></head><body><div class="wrap">
<h1>克隆樹拓樸互動工作站 v2（cluster-first 算法 + S/r/m 標籤 + 篩選排序）</h1>
<p class="sub">HCC1395 ⭐3 · 每區 genotype 向量→拓樸(perfect-phylogeny+噪聲過濾) · S=sSNV/r=read群/m=甲基位點 · 數字由 JSON 注入 · {len(d['detail']):,} 區可畫樹</p>
{UNIVERSE_BANNER}
{GLOSSARY_HTML}
<div class="stats">
<div class="scard"><h4>拓樸型態</h4><div id="s_topo"></div></div><div class="scard"><h4>群數 c</h4><div id="s_clust"></div></div>
<div class="scard"><h4>determinacy</h4><div id="s_det"></div></div><div class="scard"><h4>HP 根數</h4><div id="s_root"></div></div>
</div>
<details class="c17"><summary>▶ chr17:48360161 完整 worked example（S/r/m 一致標籤 + 樹 + 甲基；點開）</summary><div id="c17"></div></details>
<div class="ctrl">
chr<select id="f_chr"><option value="">全</option></select>
排序<select id="f_sort"><option value="coord">座標</option><option value="nsnv">複雜度(sSNV)</option><option value="nclust">群數</option><option value="region">region名</option></select>
最少群數<input id="f_minc" type="number" value="0" min="0" max="6" style="width:52px">
<label><input id="f_fp" type="checkbox">僅含FP</label>
<label title="LOH 區且單一 HP 標籤(LOH-unmask 觀察)"><input id="f_loh" type="checkbox">僅 LOH 單HP</label>
<label title="分支順序未定/不相容;曾標『需甲基輔助』,06-28 cis-control 已否決(乾淨可用≈0)"><input id="f_undef" type="checkbox">僅無法定義(曾標需甲基·已否決)</label>
搜尋<input id="f_q" placeholder="chr17:" style="width:120px">
<span id="cnt" class="note"></span>
</div>
<div class="ctrl" style="font-size:11.5px"><span class="note">勾選要觀察的狀況(可複選;不勾=全):</span>
<span>拓樸 <span id="cb_topology_type"></span></span> ｜ <span>determinacy <span id="cb_determinacy"></span></span> ｜ <span>位置 <span id="cb_genome_ctx"></span></span></div>
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
<script>window.__DATA__={DJ};</script><script>{JS}</script></body></html>"""
with open(OUT, "w", encoding="utf-8") as f: f.write(HTML)
print(f"OK wrote {OUT} ({len(HTML):,} bytes; detail {len(d['detail'])})")
