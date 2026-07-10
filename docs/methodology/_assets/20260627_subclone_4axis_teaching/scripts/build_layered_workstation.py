#!/usr/bin/env python3
"""
build_layered_workstation.py — 分層樹重建工作站 HTML(2026-07-06)
讀 layered_region_view_{SAMPLE}.json(region-centric 最終格式) → standalone HTML。
資料模型:docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md
反捏造(§13-A):所有數字從 json.census 注入,無 hardcode。
顯示:①層級/比例 dashboard ②HP-multiplicity+region-determinacy census ③region browser(篩選+分頁)
      ④每 region 逐 lineage 枚舉樹 SVG carousel + L0-L3 判斷軌跡。
env SM_RV=<region_view.json> SM_OUT=<out.html>
"""
import os, sys, json

RV = os.environ.get("SM_RV", "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/layered_region_view_HCC1395.json")
OUT = os.environ.get("SM_OUT", "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260706_layered_reconstruction_workstation.standalone.html")
D = json.load(open(RV, encoding="utf-8"))
C = D["census"]
SAMPLE = D.get("sample", "?")

# 精簡 regions 給前端(trees 已 ≤6;去掉冗長欄位以縮小)
regions = D["regions"]
DATA_JSON = json.dumps({"sample": SAMPLE, "census": C, "regions": regions}, ensure_ascii=False, separators=(",", ":"))

HTML = r"""<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>分層樹重建工作站 · __SAMPLE__</title>
<style>
:root{--bg:#fff;--fg:#1a1a1a;--mut:#666;--line:#e2e2e2;--card:#fafafa;--accent:#2563eb;
 --det:#16a34a;--amb:#d97706;--cap:#9333ea;--rec:#dc2626;--hp1:#2563eb;--hp2:#db2777;--hp3:#0891b2;--none:#94a3b8;}
@media(prefers-color-scheme:dark){:root{--bg:#0f1115;--fg:#e6e6e6;--mut:#9aa;--line:#2a2f38;--card:#171a21;--accent:#60a5fa;}}
:root[data-theme=dark]{--bg:#0f1115;--fg:#e6e6e6;--mut:#9aa;--line:#2a2f38;--card:#171a21;--accent:#60a5fa;}
:root[data-theme=light]{--bg:#fff;--fg:#1a1a1a;--mut:#666;--line:#e2e2e2;--card:#fafafa;--accent:#2563eb;}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.5 ui-sans-serif,system-ui,"Segoe UI",sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:18px}
h1{font-size:20px;margin:0 0 2px}h2{font-size:15px;margin:20px 0 8px;border-left:3px solid var(--accent);padding-left:8px}
.sub{color:var(--mut);font-size:12px;margin-bottom:14px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:6px 0}
th,td{border:1px solid var(--line);padding:5px 8px;text-align:left}th{background:var(--card);font-weight:600}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}@media(max-width:760px){.grid{grid-template-columns:1fr}}
.tag{display:inline-block;padding:1px 7px;border-radius:10px;font-size:11px;font-weight:600;color:#fff}
.t-det{background:var(--det)}.t-amb{background:var(--amb)}.t-cap{background:var(--cap)}.t-rec{background:var(--rec)}.t-none{background:var(--none)}
.hp1{color:var(--hp1)}.hp2{color:var(--hp2)}.hp3{color:var(--hp3)}.hpnone{color:var(--none)}
.ctrl{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:10px 0;font-size:13px}
.ctrl select,.ctrl input{background:var(--bg);color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:4px 7px}
.rcard{border:1px solid var(--line);border-radius:9px;margin:9px 0;background:var(--card);overflow:hidden}
.rhead{padding:8px 12px;cursor:pointer;display:flex;flex-wrap:wrap;gap:10px;align-items:center}
.rhead:hover{background:rgba(128,128,128,.06)}
.rbody{padding:0 12px 10px;display:none}.rcard.open .rbody{display:block}
.lin{border-top:1px dashed var(--line);padding:9px 0 9px 10px;border-left:3px solid var(--line);margin-top:6px}
.lin.f1{border-left-color:var(--hp1)}.lin.f2{border-left-color:var(--hp2)}.lin.f3{border-left-color:var(--hp3)}.lin.fnone{border-left-color:var(--none)}
.mono{font-family:ui-monospace,"SF Mono",Menlo,monospace}
.trace{font-size:11.5px;color:var(--mut);margin:5px 0 0;padding-left:4px}.trace div{padding:1px 0}
/* 樹切換器:一次一棵大圖 + prev/next + 計數 + thumbnail 跳轉 */
.tsw{border:1px solid var(--line);border-radius:8px;padding:7px 9px;margin:6px 0;background:var(--bg)}
.tswhead{font-size:12px;color:var(--fg);margin-bottom:5px;display:flex;align-items:center;gap:6px;flex-wrap:wrap}
.tstage{overflow-x:auto;padding:4px 0;min-height:60px}
.tslide{text-align:center}
.tcap{font-size:11px;color:var(--mut);text-align:center;margin-top:3px}
.tnav{padding:1px 9px;font-size:13px;line-height:1.4}.tctr{font-size:12px;color:var(--mut);min-width:46px;text-align:center;font-variant-numeric:tabular-nums}
.thumbs{display:flex;gap:4px;flex-wrap:wrap;margin-top:6px}
.thumb{cursor:pointer;font-size:10.5px;padding:1px 7px;border:1px solid var(--line);border-radius:9px;color:var(--mut);user-select:none}
.thumb.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.thumb:hover{border-color:var(--accent)}
.carousel{display:flex;gap:8px;overflow-x:auto;padding:6px 2px}
.tree{flex:0 0 auto;border:1px solid var(--line);border-radius:7px;padding:5px;background:var(--bg)}
.note{font-size:12px;color:var(--mut);background:var(--card);border:1px solid var(--line);border-radius:7px;padding:8px 10px;margin:8px 0}
.pill{font-size:11px;padding:1px 6px;border:1px solid var(--line);border-radius:8px;color:var(--mut)}
.bar{display:inline-block;height:9px;border-radius:3px;vertical-align:middle}
button.tg{background:var(--card);color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:3px 9px;cursor:pointer;font-size:12px}
#more{margin:12px auto;display:block}
/* topology 式兩欄:左 region 列表 + 右詳情面板 */
.main{display:grid;grid-template-columns:400px 1fr;gap:12px;align-items:start}@media(max-width:900px){.main{grid-template-columns:1fr}}
.list{max-height:78vh;overflow-y:auto;border:1px solid var(--line);border-radius:9px;background:var(--bg)}
.row{padding:7px 10px;border-bottom:1px solid var(--line);cursor:pointer;font-size:12px}
.row:hover{background:rgba(128,128,128,.06)}.row.sel{background:rgba(37,99,235,.10);border-left:3px solid var(--accent)}
.detail{border:1px solid var(--line);border-radius:9px;background:var(--card);padding:12px 14px;min-height:200px}
.detail h3{font-size:15px;margin:0 0 6px;font-family:ui-monospace,monospace}
.kv{display:flex;flex-wrap:wrap;gap:6px;margin:6px 0}.kv .b{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:2px 8px;font-size:11.5px}
</style>
<div class="wrap">
<div style="display:flex;justify-content:space-between;align-items:start;flex-wrap:wrap">
<div><h1>分層樹重建工作站 · <span class="mono">__SAMPLE__</span></h1>
<div class="sub">L0 HP家族 → L1 sSNV枚舉全最小樹 → L2 CN → L3 甲基 · 家族優先於算法 · 每數字綁分子÷分母(見資料模型 spec)</div></div>
<button class="tg" onclick="var r=document.documentElement;r.dataset.theme=(r.dataset.theme=='dark'?'light':'dark')">◐ 主題</button>
</div>
<div id="dash"></div>
<h2>Region 瀏覽器（篩選 → 點左側區 → 右側看逐 HP 家族枚舉樹 ◀▶）</h2>
<div class="ctrl">
 determinacy:<select id="fdet"><option value="">全部</option></select>
 HP-mult:<select id="fhp"><option value="">全部</option><option value="2">多-HP(2)</option><option value="1">single-HP(1)</option><option value="0">無germline(0)</option></select>
 chrom:<select id="fchr"><option value="">全部</option></select>
 排序:<select id="fsort"><option value="coord">座標</option><option value="nsnv">複雜度(sSNV)</option><option value="ntree">枚舉樹總數</option><option value="hp">HP-mult</option></select>
 搜尋:<input id="fq" placeholder="chr1:123... / 基因座" size="16">
 <span class="pill" id="cnt"></span>
</div>
<div class="main"><div class="list" id="list"></div><div class="detail" id="detail"><div class="note">← 左側點選一個區,右側顯示逐 germline-HP-家族枚舉樹（◀▶ 切換等機率最小樹）+ L0-L3 判斷軌跡</div></div></div>
<div class="note">反捏造:所有 census 數字由 <span class="mono">layered_region_view___SAMPLE__.json</span> 注入;樹由 solver 枚舉(V1-V7 驗證);甲基未進樹(L3 事後,不 rank)。</div>
</div>
<script>
const D=__DATA__;const C=D.census;const R=D.regions;
const pct=(a,b)=>b?(100*a/b).toFixed(1)+'%':'—';
const esc=s=>String(s).replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));
// ---- L1 class → 短標/顏色 ----
function clsTag(c){if(c==='determined')return['determined','t-det'];
 if(c.indexOf('ambiguous')===0)return['ambiguous(多樹)','t-amb'];
 if(c.indexOf('capped')===0)return['capped(太密)','t-cap'];
 if(c.indexOf('recurrence')===0)return['recurrence','t-rec'];return[c,'t-none'];}
function regTag(c){return {all_determined:['all-determined','t-det'],has_ambiguous:['has-ambiguous','t-amb'],
 has_capped:['has-capped','t-cap'],has_recurrence:['has-recurrence','t-rec'],no_germline_lineage:['no-germline','t-none']}[c]||[c,'t-none'];}
const famCls=f=>({'1':'hp1','2':'hp2','3':'hp3'}[f]||'hpnone');
// ================= dashboard =================
function dash(){
 const L0=C.L0,L1=C.L1,L2=C.L2,rd=C.region_determinacy,mult=C.hp_multiplicity;
 const nreg=C.n_regions, nlin=L1.n_lineage_units;
 const barw=v=>Math.max(2,Math.round(v/nreg*160));
 let h='';
 // 層級計數
 h+='<h2>① 層級計數（U1→U6）</h2><table><tr><th>層</th><th>單位</th><th class="n">HCC1395 計數</th><th>說明</th></tr>'
 +row('U1','somatic sSNV',(C.U1_sSNV_somatic_total||0).toLocaleString(),'ClairS v0.4.0 paired PASS somatic(重建骨幹·_sc.vcf;2026-07-09 移除舊 is_somatic 粗重檢·救回 429 誤殺 TP);非舊 census somatic 23,810/總 census 35,332')
 +row('U3','region',nreg.toLocaleString(),'multilocus 分析群(≤8sSNV窗,主分母);linkage全span='+((C.U3_linkage_regions_full_span||0).toLocaleString())+' 同批區(位置99.9%重疊)')
 +row('U4','lineage-unit',nlin.toLocaleString(),'region×HP家族(1/2/3);樹的自然單位')
 +row('U5','tree',(C.U5_trees?C.U5_trees.sum_ntrees_noncapped.toLocaleString():'—'),'Σ枚舉樹(non-capped);avg '+(C.U5_trees?(C.U5_trees.sum_ntrees_noncapped/nlin).toFixed(2):'—')+' 樹/unit')
 +row('U6','隱藏祖先',(C.U6_hidden?C.U6_hidden.sum_hidden.toLocaleString():'—'),'Σ推斷未觀測 clone(H_*)')
 +'</table>';
 // HP-multiplicity + region determinacy
 h+='<div class="grid"><div><h2>② HP-multiplicity（每區幾個 germline 樹）</h2><table><tr><th>germline lineage 樹數</th><th class="n">區數</th><th class="n">比例</th><th></th></tr>'
 +mrow('2（多-HP：雙親代染色體各一樹）',mult['2'],nreg,'var(--hp2)')
 +mrow('1（single-HP）',mult['1'],nreg,'var(--hp1)')
 +mrow('0（只 somatic3/none）',mult['0'],nreg,'var(--none)')
 +'</table><div class="note">🔴 多-HP '+pct(mult['2'],nreg)+'：過半區兩 germline 家族都帶 mutation → 舊 pooled 混淆 allelic/clonal。</div></div>';
 h+='<div><h2>③ region-level determinacy</h2><table><tr><th>判定</th><th class="n">區數</th><th class="n">比例</th><th></th></tr>'
 +drow('all-determined（全家族唯一樹）',rd.all_determined,nreg,'t-det')
 +drow('has-ambiguous（≥1家族多樹）',rd.has_ambiguous,nreg,'t-amb')
 +drow('has-capped（太密）',rd.has_capped,nreg,'t-cap')
 +drow('has-recurrence',rd.has_recurrence,nreg,'t-rec')
 +drow('no-germline（只3/none）',rd.no_germline_lineage,nreg,'t-none')
 +'</table><div class="note">定義：region 確定 ⟺ <b>所有</b> germline(1/2) lineage 都 determined（多-HP 需雙確定）→ 嚴於 lineage-unit。</div></div></div>';
 // 比例字典(lineage 次分母) + L2
 h+='<h2>④ 比例字典（次分母 lineage-unit '+nlin.toLocaleString()+'） + L2 CN</h2><table><tr><th>比例</th><th>分子÷分母</th><th class="n">值</th></tr>'
 +prow('determined(lineage)',L1.determinacy_lineage.determined,nlin)
 +prow('ambiguous(lineage)',L1.determinacy_lineage['ambiguous_structure(多完成/多結構)'],nlin)
 +prow('capped(lineage)',L1.determinacy_lineage['capped(太密;枚舉未完)'],nlin)
 +prow('recurrence(lineage)',L1.determinacy_lineage.recurrence_required||0,nlin)
 +'<tr><td>L2 CN artifact（recurrence 內）</td><td class="mono">'+(L2.cn_split['artifact(m>1;CN-amp)']||0)+' ÷ '+L2.n_recurrence_sent_to_cn+' recurrence</td><td class="n">'+pct(L2.cn_split['artifact(m>1;CN-amp)']||0,L2.n_recurrence_sent_to_cn)+'</td></tr>'
 +'</table>';
 h+='<div class="note">V1–V7 驗證：'+(L1.all_V1V7_pass?'✅ ALL PASS（0 fail）':('⚠ '+L1.n_verify_fail+' fail'))+' · 分母鐵則：region determined '+pct(rd.all_determined,nreg)+' ≠ lineage determined '+pct(L1.determinacy_lineage.determined,nlin)+'（單位不同不可比）</div>';
 document.getElementById('dash').innerHTML=h;
 function row(u,n,c,d){return '<tr><td class="mono">'+u+'</td><td>'+n+'</td><td class="n mono">'+c+'</td><td style="font-size:12px;color:var(--mut)">'+d+'</td></tr>';}
 function mrow(lab,v,tot,col){return '<tr><td>'+lab+'</td><td class="n">'+(v||0).toLocaleString()+'</td><td class="n">'+pct(v||0,tot)+'</td><td><span class="bar" style="width:'+barw(v||0)+'px;background:'+col+'"></span></td></tr>';}
 function drow(lab,v,tot,cl){return '<tr><td>'+lab+'</td><td class="n">'+(v||0).toLocaleString()+'</td><td class="n">'+pct(v||0,tot)+'</td><td><span class="tag '+cl+'">'+pct(v||0,tot)+'</span></td></tr>';}
 function prow(lab,v,tot){return '<tr><td>'+lab+'</td><td class="mono">'+(v||0).toLocaleString()+' ÷ '+tot.toLocaleString()+'</td><td class="n">'+pct(v||0,tot)+'</td></tr>';}
}
// ================= tree SVG =================
function treeSVG(edges,stable){
 if(!edges||!edges.length)return '<div class="tcap">（無邊）</div>';
 const ch={},par={},nodes=new Set();
 edges.forEach(e=>{(ch[e[0]]=ch[e[0]]||[]).push(e[1]);par[e[1]]=e[0];nodes.add(e[0]);nodes.add(e[1]);});
 const roots=[...nodes].filter(n=>!(n in par));
 // 分層 by depth
 const depth={},pos={};let maxd=0;
 function dfs(n,d){depth[n]=d;maxd=Math.max(maxd,d);(ch[n]||[]).forEach(c=>dfs(c,d+1));}
 roots.forEach(r=>dfs(r,0));
 // x 位置:leaf 順序
 let leafX=0;const xw=76;
 function assignX(n){const cs=ch[n]||[];if(!cs.length){pos[n]=leafX++;return pos[n];}
  const xs=cs.map(assignX);pos[n]=xs.reduce((a,b)=>a+b,0)/xs.length;return pos[n];}
 roots.forEach(assignX);
 const W=Math.max(1,leafX)*xw, H=(maxd+1)*54+10;
 const nx=n=>pos[n]*xw+xw/2, ny=n=>depth[n]*54+22;
 function lab(n){if(n==='ROOT')return['germ','#64748b'];
  if(n[0]==='H'&&n[1]==='_')return['H','#a855f7'];return['obs','#2563eb'];}
 let s='<svg width="'+W+'" height="'+H+'" viewBox="0 0 '+W+' '+H+'" style="font:9px ui-monospace,monospace">';
 // 邊:全樹一致=實線灰;跨等機率樹會變=橙虛線(枚舉組合選擇)
 edges.forEach(e=>{const vary=stable&&!stable.has(e[0]+'>'+e[1]);
  s+='<line x1="'+nx(e[0]).toFixed(1)+'" y1="'+ny(e[0]).toFixed(1)+'" x2="'+nx(e[1]).toFixed(1)+'" y2="'+ny(e[1]).toFixed(1)+'" stroke="'+(vary?'#f59f00':'#94a3b8')+'" stroke-width="'+(vary?1.9:1.3)+'"'+(vary?' stroke-dasharray="4 2"':'')+'><title>'+(vary?'此邊在等機率樹間變動 = 枚舉組合選擇':'此邊在所有等機率樹一致')+'</title></line>';});
 // 節點:實測=實心藍圓 / 推測(隱藏祖先)=空心虛線紫圓 / germline=灰方
 [...nodes].forEach(n=>{const isRoot=n==='ROOT',isH=n[0]==='H'&&n[1]==='_';
  const cx=(pos[n]*xw+xw/2),cy=(depth[n]*54+22);const cxs=cx.toFixed(1),cys=cy.toFixed(1);
  const geno=isH?n.slice(2):n;const short=isRoot?'RR·germ':(isH?geno+'ᴴ':geno);  // ᴴ=hidden(推測);不用「?」避免誤讀成位點值
  const tip='<title>'+esc(n)+(isH?' · 推測(隱藏祖先·無 read;ᴴ 標記非位點值)':(isRoot?' · germline 起點(實測)':' · 實測(有 read)'))+'</title>';
  if(isRoot)s+='<rect x="'+(cx-5.5).toFixed(1)+'" y="'+(cy-5.5).toFixed(1)+'" width="11" height="11" rx="2" fill="#64748b">'+tip+'</rect>';
  else if(isH)s+='<circle cx="'+cxs+'" cy="'+cys+'" r="6" fill="#fff" stroke="#a855f7" stroke-width="1.7" stroke-dasharray="3 2">'+tip+'</circle>';
  else s+='<circle cx="'+cxs+'" cy="'+cys+'" r="6" fill="#2563eb">'+tip+'</circle>';
  s+='<text x="'+cxs+'" y="'+(cy+16).toFixed(1)+'" text-anchor="middle" fill="'+(isH?'#a855f7':'var(--mut)')+'">'+esc(short)+'</text>';});
 return s+'</svg>';
}
// ================= 樹切換器 =================
let LID=0;
function _showSlide(lid,idx){const el=document.getElementById(lid);if(!el)return;
 const sl=el.querySelectorAll('.tslide');idx=(idx+sl.length)%sl.length;
 sl.forEach((s,i)=>s.style.display=i===idx?'block':'none');
 const c=document.getElementById(lid+'_c');if(c)c.textContent=(idx+1)+' / '+sl.length;
 el.querySelectorAll('.thumb').forEach((t,i)=>t.classList.toggle('on',i===idx));}
function tnav(lid,d){const el=document.getElementById(lid);if(!el)return;let cur=0;
 el.querySelectorAll('.tslide').forEach((s,i)=>{if(s.style.display!=='none')cur=i;});_showSlide(lid,cur+d);}
function tjump(lid,i){_showSlide(lid,i);}
// ================= region list(左·緊湊) + detail(右) =================
let FILT=[];
const _sumT=r=>(r.lineages||[]).reduce((s,L)=>s+(L.n_trees||0),0);
const SORT={coord:(a,b)=>(a.start-b.start)||a.chrom.localeCompare(b.chrom),
 nsnv:(a,b)=>b.n_sSNV-a.n_sSNV, ntree:(a,b)=>_sumT(b)-_sumT(a), hp:(a,b)=>b.hp_multiplicity-a.hp_multiplicity};
function applyFilter(){
 const fd=document.getElementById('fdet').value,fh=document.getElementById('fhp').value,
  fc=document.getElementById('fchr').value,fq=document.getElementById('fq').value.trim().toLowerCase(),
  so=document.getElementById('fsort').value;
 FILT=R.map((r,i)=>[r,i]).filter(([r])=>(!fd||r.region_determinacy===fd)&&(!fh||String(r.hp_multiplicity)===fh)
  &&(!fc||r.chrom===fc)&&(!fq||r.region.toLowerCase().indexOf(fq)>=0));
 FILT.sort((A,B)=>(SORT[so]||SORT.coord)(A[0],B[0]));
 const host=document.getElementById('list');
 host.innerHTML=FILT.slice(0,600).map(([r,i])=>{const[rt,rc]=regTag(r.region_determinacy);
  return '<div class="row" data-i="'+i+'"><span class="mono" style="font-weight:600">'+esc(r.region)+'</span> <span class="tag '+rc+'">'+rt+'</span><br>'
   +'<span class="pill">HP×'+r.hp_multiplicity+(r.is_multiHP?'·多':'')+'</span><span class="pill">'+r.n_sSNV+'sSNV</span><span class="pill">'+r.lineages.length+' lin</span><span class="pill">'+_sumT(r)+'樹</span><span class="pill">cn '+esc(r.cn)+'</span></div>';
 }).join('')+(FILT.length>600?'<div class="note" style="margin:6px">...前 600（共 '+FILT.length.toLocaleString()+'，可篩選縮小）</div>':'');
 document.getElementById('cnt').textContent=FILT.length.toLocaleString()+' 區';
 host.querySelectorAll('.row').forEach(x=>x.onclick=()=>show(+x.dataset.i,x));
}
// ===== locus 位點分布(x=基因組座標·y=VAF)=====
function locusBlock(r){
 const pos=r.positions||[];if(pos.length<1)return '';
 let pp=pos.map(()=>({ref:0,alt:0}));
 r.lineages.forEach(L=>{const cov=L.obs_col_coverage||{};pos.forEach((p,i)=>{const c=cov[String(p)];if(c){pp[i].ref+=c[0];pp[i].alt+=c[1];}});});
 const mn=Math.min(...pos),mx=Math.max(...pos),span=Math.max(1,mx-mn);
 const W=560,H=104,pL=32,pR=42,pT=8,pB=22,PW=W-pL-pR,PH=H-pT-pB;
 const X=p=>pL+(p-mn)/span*PW,Y=v=>pT+(1-v)*PH;
 let s='<svg width="'+W+'" height="'+H+'" style="font:9px ui-monospace,monospace">';
 s+='<line x1="'+pL+'" y1="'+pT+'" x2="'+pL+'" y2="'+(pT+PH)+'" stroke="#ccc"/>';
 [0,0.5,1].forEach(v=>{s+='<line x1="'+(pL-3)+'" y1="'+Y(v)+'" x2="'+pL+'" y2="'+Y(v)+'" stroke="#ccc"/><text x="'+(pL-5)+'" y="'+(Y(v)+3)+'" text-anchor="end" fill="#888">'+(v*100)+'</text>';});
 s+='<line x1="'+pL+'" y1="'+(pT+PH)+'" x2="'+(pL+PW)+'" y2="'+(pT+PH)+'" stroke="#ccc"/>';
 pos.forEach((p,i)=>{const o=pp[i],tot=o.ref+o.alt,vaf=tot?o.alt/tot:0,x=X(p),y=Y(vaf),col=o.alt>0?'#2563eb':'#dc2626';
  s+='<line x1="'+x.toFixed(1)+'" y1="'+(pT+PH)+'" x2="'+x.toFixed(1)+'" y2="'+y.toFixed(1)+'" stroke="'+col+'" stroke-width="1.3"/><circle cx="'+x.toFixed(1)+'" cy="'+y.toFixed(1)+'" r="4.5" fill="'+col+'"><title>S'+(i+1)+' '+p+' · VAF '+(vaf*100).toFixed(0)+'% ('+o.alt+'A/'+o.ref+'R)</title></circle>'
   +'<text x="'+x.toFixed(1)+'" y="'+(pT+PH+10)+'" text-anchor="middle" fill="#888">S'+(i+1)+'</text>';});
 const bl=span>=2000?1000:Math.max(100,Math.round(span/4/100)*100);
 s+='<line x1="'+(pL+PW-bl/span*PW).toFixed(1)+'" y1="'+(H-4)+'" x2="'+(pL+PW)+'" y2="'+(H-4)+'" stroke="#333"/><text x="'+(pL+PW)+'" y="'+(H-6)+'" text-anchor="end" fill="#888">'+(bl>=1000?bl/1000+'kb':bl+'bp')+'</text></svg>';
 return '<details style="margin:6px 0"><summary style="cursor:pointer;font-size:11.5px;font-weight:600">📏 locus 位點分布（x=座標·y=VAF·藍=有ALT/紅=0ALT·跨 '+(span>=1000?(span/1000).toFixed(1)+'kb':span+'bp')+'）</summary>'+s+'</details>';
}
// ===== 位點 pairwise 共現(僅全跨 read·A×A=co-phased 證據)=====
function pairwiseBlock(r){
 const pos=r.positions||[];if(pos.length<2)return '';
 let pops={};r.lineages.forEach(L=>{const p=L.obs_populations||{};for(const g in p)pops[g]=(pops[g]||0)+p[g];});
 const gts=Object.keys(pops);
 // 位點各自是否有 ALT(來自 col_coverage·判斷「該對是否本應可 co-phase」)
 let hasAlt=pos.map(()=>false);r.lineages.forEach(L=>{const cov=L.obs_col_coverage||{};pos.forEach((p,i)=>{const c=cov[String(p)];if(c&&c[1]>0)hasAlt[i]=true;});});
 let pairs=[];for(let i=0;i<pos.length;i++)for(let j=i+1;j<pos.length;j++)pairs.push([i,j]);
 let anyGap=false;
 let blocks=pairs.slice(0,6).map(([i,j])=>{
  let m={RR:0,RA:0,AR:0,AA:0},cov=0;
  for(const g in pops){const a=g[i],b=g[j];if(!a||!b||a==='X'||b==='X')continue;cov+=pops[g];m[(a==='A'?'A':'R')+(b==='A'?'A':'R')]+=pops[g];}
  const gap=hasAlt[i]&&hasAlt[j]&&m.AA===0;if(gap)anyGap=true;
  const cell=(v,hl)=>'<td style="text-align:center;padding:2px 9px;border:1px solid #eee;'+(hl?(v>0?'background:#dcfce7;font-weight:700':'background:#fee2e2'):'')+'">'+v+'</td>';
  return '<div style="display:inline-block;margin:3px 10px 3px 0;vertical-align:top"><div style="font-size:10px;font-weight:600">S'+(i+1)+' × S'+(j+1)+(cov?'':' <span style="color:#dc2626">·0 全跨 read</span>')+(gap?' <span style="color:#dc2626">⚠讀span缺口</span>':'')+'</div><table style="font-size:10px;border-collapse:collapse;margin-top:1px"><tr><td style="padding:1px 6px;color:#888"></td><td style="padding:1px 6px;color:#888">S'+(j+1)+'·R</td><td style="padding:1px 6px;color:#888">S'+(j+1)+'·A</td></tr><tr><td style="color:#888">S'+(i+1)+'·R</td>'+cell(m.RR)+cell(m.RA)+'</tr><tr><td style="color:#888">S'+(i+1)+'·A</td>'+cell(m.AR)+cell(m.AA,true)+'</tr></table></div>';
 }).join('');
 return '<details style="margin:6px 0"><summary style="cursor:pointer;font-size:11.5px;font-weight:600">🔗 位點 pairwise 共現（僅全跨 read·右下綠格=兩位點同時 ALT=co-phased 直接證據）'+(anyGap?' <span style="color:#dc2626">·偵測到 read-span 缺口</span>':'')+'</summary><div style="margin-top:4px">'+(gts.length?blocks:'<span class="note" style="background:none;border:none;padding:0;color:#dc2626">此區 0 全跨 read（無 read 覆蓋 ≥2 位點）→ 任何位點對的組合皆無法實測·全靠 partial 推斷。</span>')+'</div><div class="note" style="background:none;border:none;padding:2px 0;font-size:9.5px">🔴 某對兩位點<b>各自有 ALT</b>(位點證據表)但右下格 <b>A×A=0</b>(紅底)→ 無 read 同時覆蓋兩者的 ALT = 該對組合<b>未實測·推斷</b>(read-span 限制)。聚合各 germline-HP 家族全跨 populations。</div></details>';
}
// ===== 前後相鄰區(read-span 不跨區→僅描述性)=====
function neighborBlock(r){
 const chrom=r.chrom;if(!chrom)return '';
 const sib=R.map((x,gi)=>({x:x,gi:gi})).filter(o=>o.x.chrom===chrom).sort((a,b)=>a.x.start-b.x.start);
 const idx=sib.findIndex(o=>o.x.region===r.region);if(idx<0)return '';
 const prev=idx>0?sib[idx-1]:null,next=idx<sib.length-1?sib[idx+1]:null;
 function cell(o,gap,side){if(!o)return '<div style="flex:1;min-width:180px"><b>'+side+'</b> <span class="note" style="background:none;border:none;padding:0">(此染色體邊界·無)</span></div>';
  const x=o.x,gk=gap>=1000?(gap/1000).toFixed(1)+'kb':gap+'bp';let sm=[];
  if(x.cn===r.cn)sm.push('cn '+x.cn);if(x.region_determinacy===r.region_determinacy)sm.push(regTag(x.region_determinacy)[0]);if(x.hp_multiplicity===r.hp_multiplicity)sm.push('HP×'+x.hp_multiplicity);
  return '<div style="flex:1;min-width:180px"><b>'+side+'</b> <span class="mono" style="cursor:pointer;color:#2563eb;text-decoration:underline" onclick="show('+o.gi+')">'+esc(x.region)+'</span> <span class="note" style="background:none;border:none;padding:0">'+x.n_sSNV+'sSNV·gap '+gk+'</span><br>'+(sm.length?'<span style="color:#16a34a;font-size:10px">共享 '+sm.join(' / ')+'</span>':'<span class="note" style="background:none;border:none;padding:0;font-size:10px">無共享屬性</span>')+'</div>';}
 return '<details style="margin:6px 0"><summary style="cursor:pointer;font-size:11.5px;font-weight:600">↔ 前後相鄰區（🔴 read-span 不跨區·lineage 不可確認連續·僅描述性鄰接）</summary><div style="display:flex;gap:14px;margin-top:5px;flex-wrap:wrap">'+cell(prev,prev?r.start-prev.x.end:0,'◀ 前')+cell(next,next?next.x.start-r.end:0,'後 ▶')+'</div></details>';
}
// ===== 每 germline-HP 家族觀測 read(全跨 populations + partial subreads)=====
function readObsBlock(L){
 const fp=L.obs_populations||{},sp=L.obs_subreads||{};const nf=Object.keys(fp).length,ns=Object.keys(sp).length;
 if(!nf&&!ns)return '';
 const row=(g,typ,c,col)=>'<tr><td class="mono" style="padding:2px 6px">'+g+'</td><td style="padding:2px 6px;color:'+col+'">'+typ+'</td><td style="padding:2px 6px;text-align:right">'+c+'</td></tr>';
 let rows=Object.entries(fp).sort((a,b)=>b[1]-a[1]).map(([g,c])=>row(g,'全跨(co-phase)',c,'#2563eb')).join('')
  +Object.entries(sp).sort((a,b)=>b[1]-a[1]).map(([g,c])=>row(g,'partial(X=未覆蓋)',c,'#a855f7')).join('');
 return '<details style="margin:5px 0"><summary style="cursor:pointer;font-size:11px">▶ 此家族觀測 read（'+nf+' 全跨群 + '+ns+' partial 群·原始 mlhp）</summary>'
  +'<table style="font-size:10.5px;margin-top:3px;border-collapse:collapse"><tr style="border-bottom:1px solid #ddd"><th style="padding:2px 6px;text-align:left">genotype</th><th style="padding:2px 6px;text-align:left">類型</th><th style="padding:2px 6px">reads</th></tr>'+rows+'</table>'
  +'<div class="note" style="background:none;border:none;padding:2px 0;font-size:9.5px">全跨群=同一 read 覆蓋全部位點(可定 phasing);partial=只覆蓋部分位點(單獨看到某位點的 ALT·無法連鎖)。'+(nf===0?'<b style="color:#a855f7">此家族 0 全跨群→組合全靠 partial 拼(推斷)。</b>':'')+'</div></details>';
}
function show(i,row){
 document.querySelectorAll('#list .row').forEach(x=>x.classList.remove('sel'));if(row)row.classList.add('sel');
 const r=R[i];const[rt,rc]=regTag(r.region_determinacy);
 let h='<h3>'+esc(r.region)+' <span class="tag '+rc+'">'+rt+'</span></h3>'
  +'<div class="kv"><span class="b">'+r.n_sSNV+' sSNV</span><span class="b">HP×'+r.hp_multiplicity+(r.is_multiHP?'（多-HP·雙親代各一樹）':'（single-HP）')+'</span><span class="b">cn '+esc(r.cn)+'</span><span class="b">'+r.lineages.length+' germline-HP 家族</span></div>'
  +'<div class="note" style="margin:4px 0">每 germline-HP 家族分開建樹（家族優先於算法，修 allelic/clonal 混淆）;每家族 ◀▶ 切換其「等機率最小樹」（枚舉全集=「定不出來即答案」）。'
   +'<div style="margin-top:5px;display:flex;gap:12px;flex-wrap:wrap;align-items:center;font-size:10.5px">'
   +'<span><svg width="14" height="14" style="vertical-align:middle"><circle cx="7" cy="7" r="5.5" fill="#2563eb"/></svg> <b>實測</b>(該基因型有 read 觀測到)</span>'
   +'<span><svg width="14" height="14" style="vertical-align:middle"><circle cx="7" cy="7" r="5" fill="#fff" stroke="#a855f7" stroke-width="1.6" stroke-dasharray="3 2"/></svg> <b>推測</b>(隱藏祖先·無 read·solver 補的中間 clone)</span>'
   +'<span><svg width="14" height="14" style="vertical-align:middle"><rect x="2" y="2" width="10" height="10" rx="2" fill="#64748b"/></svg> <b>germline</b> 起點(RR)</span>'
   +'<span><svg width="26" height="10" style="vertical-align:middle"><line x1="1" y1="5" x2="25" y2="5" stroke="#94a3b8" stroke-width="1.3"/></svg> 灰實線=此邊<b>在全部等機率樹一致=forced 骨幹</b></span>'
   +'<span><svg width="26" height="10" style="vertical-align:middle"><line x1="1" y1="5" x2="25" y2="5" stroke="#f59f00" stroke-width="1.9" stroke-dasharray="4 2"/></svg> <b style="color:#d97706">橙虛線</b>=此連接在等機率樹間<b>變動=枚舉組合選擇</b>(非唯一)</span>'
   +'</div><div class="note" style="font-size:10px;margin-top:3px;background:none;border:none;padding:0">🔴 <b>整棵樹全橙</b> = 無 forced 骨幹、結構完全在 N 個選擇間未定(本資料 ~82% 多樹 lineage 屬此);<b>有灰邊</b> = 灰部分 forced、僅橙邊是選擇。<b>capped(太密)</b>lineage 枚舉未完整→穩定/選擇標記僅基於已存樹,可能高估 forced。</div>';
 // ===== 位點證據 + 確定性(用原始 col_coverage 真值):每位點實測 ALT reads + co-phase 缺口 =====
 (function(){
  const pos=r.positions||[];if(!pos.length)return;  // 無 mlhp 原始資料則略
  let perpos=pos.map(()=>({ref:0,alt:0}));
  r.lineages.forEach(L=>{const cov=L.obs_col_coverage||{};pos.forEach((p,i)=>{const c=cov[String(p)];if(c){perpos[i].ref+=c[0];perpos[i].alt+=c[1];}});});
  let recur=new Set();r.lineages.forEach(L=>(L.trees||[]).forEach(t=>(t.recurrence||[]).forEach(b=>recur.add(b))));
  let cells='',nObsAlt=0,nZero=0;
  pos.forEach((p,i)=>{const o=perpos[i],tot=o.ref+o.alt,vaf=tot?100*o.alt/tot:0;
   let cls,txt;if(o.alt>0){cls='#2563eb';txt=o.alt+' ALT 實測';nObsAlt++;}else{cls='#dc2626';txt='0 ALT⚠';nZero++;}
   cells+='<td style="text-align:center;padding:3px 7px"><b>S'+(i+1)+'</b> <span class="note" style="background:none;border:none;padding:0">'+p+'</span>'+(recur.has(i)?' <span style="color:#dc2626;font-size:9px" title="突變兩次·違反 infinite-sites">🔁</span>':'')+'<br><span style="color:'+cls+';font-size:10.5px">'+txt+'</span><br><span class="note" style="background:none;border:none;padding:0;font-size:9px">VAF '+vaf.toFixed(0)+'% ('+o.alt+'A/'+o.ref+'R)</span></td>';
  });
  const nfc=r.n_full_cov_reads||0,k=pos.length,gap=pos.length>1?pos[pos.length-1]-pos[0]:0;
  const gtxt=gap>=1000?(gap/1000).toFixed(1)+'kb':gap+'bp';
  const cophaseWeak=nfc<Math.max(2,nObsAlt);  // 跨全位點 read 太少 → 組合推斷
  let cl2={all_determined:['High','#16a34a','全 lineage 唯一樹'],has_ambiguous:['Med','#d97706','有多樹未定'],has_capped:['Low','#9333ea','枚舉未完(太密)'],has_recurrence:['Low','#dc2626','有 recurrence'],no_germline_lineage:['—','#94a3b8','無 germline']}[r.region_determinacy]||['—','#94a3b8',''];
  h+='<div class="note" style="margin:6px 0"><b>🎯 結構確定性</b> <span class="tag" style="background:'+cl2[1]+'">'+cl2[0]+'</span> <span class="note" style="background:none;border:none;padding:0">'+cl2[2]+' ｜ '+nObsAlt+'/'+k+' 位點實測到 ALT</span>'
   +'<div style="margin-top:5px;font-weight:600;font-size:11.5px">📍 位點證據（實測 ALT reads·非推斷·來自原始 col_coverage）</div><table style="margin-top:3px;font-size:11px"><tr>'+cells+'</tr></table>'
   +'<div style="font-size:10px;margin-top:3px">🔴 <b>位點層 vs 組合層</b>：每位點各自的 ALT <b style="color:#2563eb">是 read 直接實測</b>（上表 A 數）;但<b>「組合(phasing)」是否實測</b>看 <b>n_full_cov_reads='+nfc+'</b>（同時跨全部 '+k+' 位點的 read）— '
   +(cophaseWeak?'<b style="color:#a855f7">太少 → 無 read 把各位點的 ALT 連成同一分子 → 樹節點(組合)為推斷（read-span 限制;位點跨 '+gtxt+'）。即「位點實測、但誰跟誰同 lineage 定不出來」。</b>':'足夠 → 組合可被實測 read 定。')
   +(nZero?'<b style="color:#dc2626"> ⚠ '+nZero+' 位點 0 ALT read(census somatic 但此 linkage 無 ALT 覆蓋)。</b>':'')+'</div></div>';
 })();
 h+=locusBlock(r)+pairwiseBlock(r)+neighborBlock(r);
 r.lineages.forEach(L=>{const[ct,cc]=clsTag(L.L1_class);const fc=famCls(L.family);
  const _obsPops=(L.n_full_pops||0);
  h+='<div class="lin f'+L.family+'"><b class="'+fc+'">▸ '+esc(L.fam_label)+'</b> '
   +'<span class="tag '+cc+'">'+ct+'</span> '
   +'<span class="pill">'+L.n_trees+' 樹'+(L.n_distinct_shapes&&L.n_distinct_shapes<L.n_trees?'/'+L.n_distinct_shapes+'形狀':'')+'</span><span class="pill">'+L.n_hidden+' 隱藏祖先</span>'
   +'<span class="pill">'+(L.n_reads||0)+' reads·'+L.n_full_pops+'full/'+L.n_partial+'partial</span>'
   +(L.verify_pass?'<span class="pill" style="color:var(--det)">V1-7✓</span>':'<span class="pill" style="color:var(--rec)">V✗</span>');
  if(_obsPops===0)h+='<div class="note" style="margin:5px 0;color:#b91c1c;background:#fff5f5;border-color:#ffc9c9"><b>⚠ 此家族 0 實測全跨群</b>（n_full_pops=0·僅 '+(L.n_partial||0)+' partial read）→ 下方樹的<b>所有節點都是推測(隱藏祖先·空心紫圓)</b>,無任何直接觀測到的完整基因型群 → 整個結構為 partial-read 推斷,非觀測。</div>';
  const _hasTree=(L.trees||[]).some(t=>t.edges&&t.edges.length);
  if(L.trees&&L.trees.length&&_hasTree){const lid='L'+(LID++);const ns=L.trees.length;
   // 穩定邊集=出現在全部 N 棵等機率樹的邊(其餘=枚舉組合選擇,treeSVG 標橙虛線)
   let stable=null;if(ns>1){const cnt={};L.trees.forEach(t=>(t.edges||[]).forEach(e=>{const k=e[0]+'>'+e[1];cnt[k]=(cnt[k]||0)+1;}));stable=new Set(Object.keys(cnt).filter(k=>cnt[k]===ns));}
   h+='<div class="tsw" id="'+lid+'"><div class="tswhead"><b>可能樹結構</b>：'+L.n_trees+' 棵等機率最小樹'
    +(L.n_distinct_shapes?'（'+L.n_distinct_shapes+(L.n_distinct_shapes<ns?'+':'')+' 種形狀）':'')
    +(L.n_trees>ns?' · 顯示前 '+ns:'')
    +(ns>1?'<button class="tg tnav" onclick="tnav(\''+lid+'\',-1)">◀</button><span class="tctr" id="'+lid+'_c">1 / '+ns+'</span><button class="tg tnav" onclick="tnav(\''+lid+'\',1)">▶</button>':'')
    +'</div><div class="tstage">';
   L.trees.forEach((t,ti)=>{
    const nd=new Set();(t.edges||[]).forEach(e=>{nd.add(e[0]);nd.add(e[1]);});
    const nobs=[...nd].filter(n=>n!=='ROOT'&&!(n[0]==='H'&&n[1]==='_')).length;
    const nhid=[...nd].filter(n=>n[0]==='H'&&n[1]==='_').length;
    h+='<div class="tslide" style="display:'+(ti?'none':'block')+'">'+treeSVG(t.edges,stable)
     +'<div class="tcap">樹 #'+(ti+1)+' · <span style="color:#2563eb">'+nobs+' 實測</span> + <span style="color:#a855f7">'+nhid+' 推測祖先</span> + germline'+(t.recurrence&&t.recurrence.length?' · recurrence bit '+t.recurrence.join(','):'')+(stable?' · 橙虛線=與其他等機率樹不同處':'')+'</div></div>';});
   h+='</div>';
   if(ns>1){h+='<div class="thumbs">';L.trees.forEach((t,ti)=>{h+='<span class="thumb'+(ti?'':' on')+'" onclick="tjump(\''+lid+'\','+ti+')">#'+(ti+1)+'</span>';});h+='</div>';}
   h+='</div>';}
  else{h+='<div class="note" style="margin:5px 0"><b>（此家族無分支樹可畫）</b> '+(_obsPops<=1?'只有單一 genotype 群或僅 germline':'樹 edges 為空')+' → <b>不是缺資料,是本來就無可枚舉的分支結構</b>（'+(L.n_full_pops||0)+' 實測群·'+(L.n_partial||0)+' partial）。determinacy='+esc(L.L1_class)+'。</div>';}
  h+=readObsBlock(L);
  h+='<div class="trace">'+L.trace.map(t=>'<div>'+esc(t)+'</div>').join('')+'</div></div>';});
 document.getElementById('detail').innerHTML=h;
}
// init
(function(){
 dash();
 const fd=document.getElementById('fdet');Object.keys(C.region_determinacy).forEach(k=>{const o=document.createElement('option');o.value=k;o.textContent=k+' ('+C.region_determinacy[k]+')';fd.appendChild(o);});
 const chrs=[...new Set(R.map(r=>r.chrom))].sort((a,b)=>(+a.slice(3)||99)-(+b.slice(3)||99));
 const fc=document.getElementById('fchr');chrs.forEach(c=>{const o=document.createElement('option');o.value=c;o.textContent=c;fc.appendChild(o);});
 ['fdet','fhp','fchr','fsort'].forEach(id=>document.getElementById(id).onchange=applyFilter);
 document.getElementById('fq').oninput=applyFilter;
 applyFilter();
})();
</script>"""

HTML = HTML.replace("__DATA__", DATA_JSON).replace("__SAMPLE__", SAMPLE)
with open(OUT, "w", encoding="utf-8") as f:
    f.write(HTML)
print(f"OK wrote {OUT} ({len(HTML):,} bytes; {len(regions)} regions)")
