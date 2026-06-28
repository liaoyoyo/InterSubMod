#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互動式拓樸工作站 HTML：每區 genotype 向量 + 克隆樹 SVG + 篩選 + 全域統計。
資料來自 topology_per_region.json（topology_analysis.py 產出，cluster-first 算法）。
§13-A：數字由 JSON 注入，前端只渲染不捏造。
用法：python3 build_topology_workstation.py
輸出：../../20260628_topology_workstation.standalone.html
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "20260628_topology_workstation.standalone.html"))

d = json.load(open(os.path.join(DATA, "topology_per_region.json"), encoding="utf-8"))
stats = d["stats"]; detail = d["detail"]
DATA_JSON = json.dumps({"stats": stats, "detail": detail}, ensure_ascii=False)

CSS = """
*{box-sizing:border-box}body{margin:0;font-family:-apple-system,"Segoe UI","Noto Sans TC","Microsoft JhengHei",sans-serif;color:#212529;background:#f8f9fa}
.wrap{max-width:1280px;margin:0 auto;padding:18px}
h1{font-size:20px;margin:.2em 0}.sub{color:#868e96;font-size:12.5px}
.stats{display:flex;gap:14px;flex-wrap:wrap;margin:12px 0}
.scard{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:10px 12px;min-width:215px}
.scard h4{margin:0 0 6px;font-size:12px;color:#495057}
.bar{display:flex;align-items:center;gap:6px;font-size:11px;margin:2px 0}
.bar i{height:11px;background:#1c7ed6;border-radius:2px;display:inline-block}
.ctrl{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin:10px 0;font-size:13px}
.ctrl select,.ctrl input{padding:4px 8px;border:1px solid #ced4da;border-radius:6px;font-size:13px}
.main{display:grid;grid-template-columns:380px 1fr;gap:14px}@media(max-width:840px){.main{grid-template-columns:1fr}}
.list{background:#fff;border:1px solid #dee2e6;border-radius:8px;max-height:74vh;overflow:auto}
.row{padding:7px 10px;border-bottom:1px solid #f1f3f5;cursor:pointer;font-size:12.5px}
.row:hover{background:#e7f5ff}.row.sel{background:#d0ebff}
.row b{color:#1c7ed6}.tag{font-size:10px;padding:1px 6px;border-radius:10px;margin-left:4px}
.t_linear{background:#d3f9d8;color:#2b8a3e}.t_branched{background:#e5dbff;color:#5f3dc4}.t_star{background:#fff3bf;color:#b08900}.t_single{background:#f1f3f5;color:#868e96}.t_inc{background:#ffe3e3;color:#c92a2a}
.detail{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:16px;min-height:74vh}
.kv{display:flex;gap:14px;flex-wrap:wrap;margin:6px 0}.kv .b{background:#f1f3f5;border-radius:6px;padding:6px 11px;font-size:12px}
table{border-collapse:collapse;font-size:12px;margin:8px 0}th,td{border:1px solid #dee2e6;padding:4px 9px}th{background:#f1f3f5}
.note{font-size:12px;color:#868e96}.mono{font-family:ui-monospace,Menlo,monospace}
.legend{font-size:11.5px;color:#495057;margin-top:6px}
"""

JS = """
const D=window.__DATA__;
const TT={'linear(全直系)':'t_linear','branched(直系+姊妹)':'t_branched','star(全姊妹)':'t_star','single':'t_single','incompatible':'t_inc','germline_only':'t_single'};
function el(id){return document.getElementById(id)}
function bars(obj,max){let m=Math.max(...Object.values(obj),1);return Object.entries(obj).sort((a,b)=>b[1]-a[1]).slice(0,max||7).map(([k,v])=>`<div class="bar"><i style="width:${Math.max(4,90*v/m)}px"></i><span>${k}: <b>${v}</b></span></div>`).join('')}
// stats
el('s_topo').innerHTML=bars(D.stats.topology_type);
el('s_clust').innerHTML=bars(Object.fromEntries(Object.entries(D.stats.n_clusters).map(([k,v])=>['c='+k,v])));
el('s_det').innerHTML=bars(D.stats.determinacy);
el('s_root').innerHTML=bars(D.stats.n_roots);
// filters
let det=D.detail;
function opts(key){return [...new Set(det.map(r=>r[key]))].sort()}
['topology_type','determinacy'].forEach(k=>{let s=el('f_'+k);opts(k).forEach(o=>{let op=document.createElement('option');op.value=o;op.textContent=o;s.appendChild(op)})});
function render(){
  let tt=el('f_topology_type').value,dd=el('f_determinacy').value,q=el('f_q').value.trim(),mc=+el('f_minc').value;
  let f=det.filter(r=>(!tt||r.topology_type==tt)&&(!dd||r.determinacy==dd)&&(!q||r.region.includes(q))&&r.n_clusters>=mc);
  el('cnt').textContent=f.length+' 區';
  el('list').innerHTML=f.slice(0,600).map((r,i)=>`<div class="row" data-i="${det.indexOf(r)}"><b>${r.region}</b> <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type.split('(')[0]}</span><br><span class="note">${r.n_sSNV}sSNV · c=${r.n_clusters} · ${r.haplotypes} · ${r.cn} · ${r.determinacy.split('(')[0].split('_')[0]}</span></div>`).join('')+(f.length>600?'<div class="note" style="padding:8px">...僅顯示前 600（共 '+f.length+'）</div>':'');
  el('list').querySelectorAll('.row').forEach(x=>x.onclick=()=>show(+x.dataset.i,x));
}
function treeSVG(r){
  // build children map from edges (ROOT 為 germline)
  let ch={},all=new Set();(r.edges||[]).forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT')all.add(p)});
  if(!all.size) return '<div class="note">此區無 genotype-向量樹（單群/germline-only）</div>';
  let depth={},pos={},leaf=0;
  function lay(n,dp){depth[n]=dp;let kids=(ch[n]||[]).sort();if(!kids.length){pos[n]=leaf++;return pos[n]}let xs=kids.map(k=>lay(k,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]}
  let roots=ch['ROOT']||[];let gx=0;roots.forEach(rt=>lay(rt,1));gx=roots.length?(Math.min(...roots.map(rt=>pos[rt]))+Math.max(...roots.map(rt=>pos[rt])))/2:0;
  let nodes=[...all];let W=Math.max(280,(leaf||1)*120),maxd=Math.max(...nodes.map(n=>depth[n]),1),H=70+maxd*78;
  let X=p=>40+p*120,Y=dp=>30+dp*78;
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}">`;
  // germline 根
  s+=`<circle cx="${X(gx)}" cy="${Y(0)}" r="16" fill="#fff" stroke="#495057"/><text x="${X(gx)}" y="${Y(0)+4}" text-anchor="middle" font-size="10">germ</text>`;
  roots.forEach(rt=>{s+=`<line x1="${X(gx)}" y1="${Y(0)+16}" x2="${X(pos[rt])}" y2="${Y(1)-22}" stroke="#adb5bd"/>`});
  (r.edges||[]).forEach(([p,c])=>{if(p!='ROOT')s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+22}" x2="${X(pos[c])}" y2="${Y(depth[c])-22}" stroke="#adb5bd"/>`});
  nodes.forEach(n=>{let cnt=r.populations[n]||'';s+=`<rect x="${X(pos[n])-44}" y="${Y(depth[n])-20}" width="88" height="40" rx="6" fill="#e7f5ff" stroke="#1c7ed6"/><text x="${X(pos[n])}" y="${Y(depth[n])-3}" text-anchor="middle" font-size="12" class="mono" font-weight="600">${n}</text><text x="${X(pos[n])}" y="${Y(depth[n])+13}" text-anchor="middle" font-size="10" fill="#495057">${cnt} reads</text>`});
  s+='</svg>';return s;
}
function show(i,row){
  el('list').querySelectorAll('.row').forEach(x=>x.classList.remove('sel'));if(row)row.classList.add('sel');
  let r=det[i];
  let pt=Object.entries(r.populations).sort((a,b)=>b[1]-a[1]).map(([g,c])=>`<tr><td class="mono">${g}</td><td>${c}</td><td class="mono">${[...g].map((ch,j)=>ch=='A'?('m'+(j+1)):'·').filter(x=>x!='·').join(',')||'germline'}</td></tr>`).join('');
  el('detail').innerHTML=`<h3>${r.region} <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type}</span></h3>
   <div class="kv"><div class="b">${r.n_sSNV} sSNV</div><div class="b">span ${(r.span/1000).toFixed(1)}kb</div><div class="b">c=${r.n_clusters} 群</div><div class="b">HP: ${r.haplotypes}</div><div class="b">CN: ${r.cn}</div><div class="b">${r.determinacy}</div>${r.drop_noise_frac>0?`<div class="b">噪聲過濾 ${(r.drop_noise_frac*100).toFixed(0)}%</div>`:''}${r.ambig_nodes>0?`<div class="b" style="background:#fff3bf">⚠ 順序未定節點 ${r.ambig_nodes}(缺中間群)</div>`:''}</div>
   <h4>克隆樹（germline→…；節點=基因型向量·read 數）</h4>${treeSVG(r)}
   <div class="legend">A=帶突變/R=不帶（位置=區內排序 sSNV）；直系=往下、姊妹=同層分叉。tree_shape(pairwise): ${r.tree_shape}</div>
   <h4>觀察到的細胞群（genotype 向量）</h4><table><tr><th>向量</th><th>reads</th><th>攜帶突變</th></tr>${pt}</table>`;
}
['f_topology_type','f_determinacy','f_q','f_minc'].forEach(id=>{let e=el(id);e.oninput=render;e.onchange=render});
render();
"""

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>克隆樹拓樸工作站 — sSNV 重建（HCC1395 ⭐3）</title><style>{CSS}</style></head><body><div class="wrap">
<h1>克隆樹拓樸互動工作站（cluster-first 算法）</h1>
<p class="sub">HCC1395 ⭐3 · 每區 genotype 向量→拓樸（perfect-phylogeny laminar 解 + population 噪聲過濾）· 數字由 topology_per_region.json 注入 · {len(detail):,} 區可畫樹</p>
<div class="stats">
<div class="scard"><h4>拓樸型態</h4><div id="s_topo"></div></div>
<div class="scard"><h4>群數 (n_clusters)</h4><div id="s_clust"></div></div>
<div class="scard"><h4>determinacy</h4><div id="s_det"></div></div>
<div class="scard"><h4>HP 根數</h4><div id="s_root"></div></div>
</div>
<div class="ctrl">
拓樸<select id="f_topology_type"><option value="">全部</option></select>
determinacy<select id="f_determinacy"><option value="">全部</option></select>
最少群數<input id="f_minc" type="number" value="0" min="0" max="6" style="width:60px">
搜尋區<input id="f_q" placeholder="chr17:..." style="width:150px">
<span id="cnt" class="note"></span>
</div>
<div class="main"><div class="list" id="list"></div><div class="detail" id="detail"><div class="note">← 左側點選一個區查看克隆樹</div></div></div>
<p class="note" style="margin-top:14px">⚠ 證據層級：A_determined=單分子向量確定；B_pairwise=pairwise 拼接（非單分子整樹）；C_underdetermined=資訊不足多樹相容。甲基不參與拓樸裁決（cis-confounded，見 methylation framework）。⭐3 單樣本上限。</p>
</div>
<script>window.__DATA__={DATA_JSON};</script><script>{JS}</script></body></html>"""

with open(OUT, "w", encoding="utf-8") as f:
    f.write(HTML)
print(f"OK wrote {OUT} ({len(HTML):,} bytes; detail {len(detail)} regions)")
