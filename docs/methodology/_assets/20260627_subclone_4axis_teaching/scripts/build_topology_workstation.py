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
DJ = json.dumps({"stats": d["stats"], "detail": d["detail"], "chr17": d["chr17_worked"], "chroms": d["chroms"]}, ensure_ascii=False)
B = acc["buckets"]
UNIVERSE_BANNER = f"""<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:11px 14px;margin:10px 0;font-size:12.5px">
<b>全 sSNV 宇宙帳本（{acc['universe_total']:,} = TP {acc['tp']:,} + FP {acc['fp']:,}；無遺漏 sum-check ✓）</b><br>
<div style="display:flex;gap:18px;flex-wrap:wrap;margin-top:6px">
<span>🟢 <b>linked {B['linked']['n']:,}（{B['linked']['pct']}%）</b> 可建樹→拓樸分析（下方）</span>
<span>🟡 <b>underpowered {B['underpowered']['n']:,}（{B['underpowered']['pct']}%）</b> 有 partner 無共讀→<b>加深覆蓋可救</b>；有 CCF（clonal {B['underpowered']['ccf_tendency'].get('high_ge0.4(clonal)',0)}/mid {B['underpowered']['ccf_tendency'].get('mid_0.1-0.4',0)}/low {B['underpowered']['ccf_tendency'].get('low_lt0.1(subclonal/noise)',0)}）</span>
<span>⚪ <b>isolated {B['isolated']['n']:,}（{B['isolated']['pct']}%）</b> read-span 內無 partner（Tier-R 樹外）；<b>有 caller VAF 可刻畫 + 可能 Tier-PS partner</b></span>
</div>
<div class="note" style="margin-top:6px">🔑 單位點<b>非「全無法處理」</b>：underpowered 有 CCF+可深覆蓋救；isolated 有 caller VAF 可放 clonal 譜 + 可能 same-PS(&gt;50kb) 連鎖(Tier-PS 未做)。真正拓樸-dead = isolated 中無 same-PS partner 者（待 Tier-PS 量化）。</div></div>"""

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
function bars(o,m){let mx=Math.max(...Object.values(o),1);return Object.entries(o).sort((a,b)=>b[1]-a[1]).slice(0,m||6).map(([k,v])=>`<div class="bar"><i style="width:${Math.max(3,84*v/mx)}px"></i>${k}: <b>${v}</b></div>`).join('')}
el('s_topo').innerHTML=bars(D.stats.topology_type);el('s_clust').innerHTML=bars(Object.fromEntries(Object.entries(D.stats.n_clusters).map(([k,v])=>['c='+k,v])));
el('s_det').innerHTML=bars(D.stats.determinacy);el('s_root').innerHTML=bars(D.stats.n_roots);
// chr17 worked panel
(function(){let c=D.chr17;
 let st=c.snvs.map(s=>`<tr><td class="mono"><b>${s.S}</b></td><td class="mono">${s.pos}</td><td>${s.change}</td><td>${s.role}</td><td>VAF ${s.vaf}</td><td>${s.hp}</td><td>${s.src}</td><td>${s.somatic_confirmed?'✓':'✗(normal有ALT)'}</td></tr>`).join('');
 let pp=c.populations.map(p=>`<tr><td class="mono">${p.vec}</td><td>${p.muts}</td><td>${p.reads}</td><td>${p.pct}%</td></tr>`).join('');
 let mm=c.sig_cpg.slice(0,16).map(x=>`<tr><td class="mono">${x.m}</td><td>${x.cpg}</td><td>${x.L1}</td><td>${x.L2}</td><td>${x.dbeta}</td></tr>`).join('');
 el('c17').innerHTML=`<div class="kv"><div class="b">locus ${c.locus}</div><div class="b">ctx ${c.genome_ctx}</div><div class="b">拓樸 ${c.topology_type}</div><div class="b">噪聲 dropped ${c.dropped_noise} reads</div></div>
  <b>S 位點(somatic sSNV)</b><table><tr><th>S</th><th>pos</th><th>變異</th><th>角色</th><th>VAF</th><th>HP</th><th>TP/FP</th><th>somatic</th></tr>${st}</table>
  <b>克隆樹</b> ${tree(c.edges,Object.fromEntries(c.populations.map(p=>[p.vec,p.reads])),c.populations.length,'H1')}
  <b>細胞群(r 分群;基因型向量·reads·佔比)</b><table><tr><th>向量</th><th>突變(S)</th><th>reads(r)</th><th>佔比</th></tr>${pp}</table>
  <b>甲基差異位點 m1..m${c.n_sig_diff_cpg}(L1 α-only vs L2 α+β;⚠ 經實測對齊 cis-genotype 軸非獨立 lineage)</b><table><tr><th>m</th><th>CpG</th><th>L1 β</th><th>L2 β</th><th>Δβ</th></tr>${mm}</table>`;
})();
// S-label 化基因型向量
function sLabels(g){let s=[...g].map((c,i)=>c=='A'?('S'+(i+1)):null).filter(Boolean);return s.length?s.join('+'):'germline'}
function tree(edges,popcount,nc,hp){
 let ch={},all=new Set();(edges||[]).forEach(([p,c])=>{(ch[p]=ch[p]||[]).push(c);all.add(c);if(p!='ROOT')all.add(p)});
 if(!all.size)return '<div class="note">單群/germline-only，無分支樹</div>';
 let depth={},pos={},leaf=0;
 function lay(n,dp){depth[n]=dp;let k=(ch[n]||[]).sort();if(!k.length){pos[n]=leaf++;return pos[n]}let xs=k.map(x=>lay(x,dp+1));pos[n]=(Math.min(...xs)+Math.max(...xs))/2;return pos[n]}
 let roots=ch['ROOT']||[];roots.forEach(r=>lay(r,1));let gx=roots.length?(Math.min(...roots.map(r=>pos[r]))+Math.max(...roots.map(r=>pos[r])))/2:0;
 let nodes=[...all],md=Math.max(...nodes.map(n=>depth[n]),1),W=Math.max(300,(leaf||1)*150),H=70+md*84;
 let X=p=>50+p*150,Y=dp=>32+dp*84,tot=Object.values(popcount).reduce((a,b)=>a+b,0)||1;
 let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}">`;
 s+=`<circle cx="${X(gx)}" cy="${Y(0)}" r="17" fill="#fff" stroke="#495057"/><text x="${X(gx)}" y="${Y(0)-2}" text-anchor="middle" font-size="9">germline</text><text x="${X(gx)}" y="${Y(0)+9}" text-anchor="middle" font-size="8" fill="#868e96">${hp||''}根</text>`;
 roots.forEach(r=>s+=`<line x1="${X(gx)}" y1="${Y(0)+17}" x2="${X(pos[r])}" y2="${Y(1)-24}" stroke="#adb5bd"/>`);
 (edges||[]).forEach(([p,c])=>{if(p!='ROOT')s+=`<line x1="${X(pos[p])}" y1="${Y(depth[p])+24}" x2="${X(pos[c])}" y2="${Y(depth[c])-24}" stroke="#adb5bd"/>`});
 nodes.forEach(n=>{let cnt=popcount[n]||0,pct=(100*cnt/tot).toFixed(0);s+=`<rect x="${X(pos[n])-52}" y="${Y(depth[n])-22}" width="104" height="44" rx="6" fill="#e7f5ff" stroke="#1c7ed6"/><text x="${X(pos[n])}" y="${Y(depth[n])-6}" text-anchor="middle" font-size="11" class="mono" font-weight="600">${sLabels(n)}</text><text x="${X(pos[n])}" y="${Y(depth[n])+8}" text-anchor="middle" font-size="9" fill="#495057">${cnt} reads · ${pct}%</text><text x="${X(pos[n])}" y="${Y(depth[n])+18}" text-anchor="middle" font-size="8" fill="#868e96">${n}</text>`});
 s+='</svg>';return s;
}
// filters + sort
let det=D.detail;
let chrsel=el('f_chr');D.chroms.forEach(c=>{let o=document.createElement('option');o.value=c;o.textContent=c;chrsel.appendChild(o)});
['topology_type','determinacy','genome_ctx'].forEach(k=>{let s=el('f_'+k);[...new Set(det.map(r=>r[k]))].sort().forEach(o=>{let op=document.createElement('option');op.value=o;op.textContent=o;s.appendChild(op)})});
const SORT={coord:(a,b)=>a.chrom.localeCompare(b.chrom,undefined,{numeric:true})||a.start-b.start,nsnv:(a,b)=>b.n_sSNV-a.n_sSNV,nclust:(a,b)=>b.n_clusters-a.n_clusters||b.n_sSNV-a.n_sSNV,region:(a,b)=>a.region.localeCompare(b.region)};
function render(){
 let ch=el('f_chr').value,tt=el('f_topology_type').value,dd=el('f_determinacy').value,gc=el('f_genome_ctx').value,
  mc=+el('f_minc').value,fp=el('f_fp').checked,q=el('f_q').value.trim(),so=el('f_sort').value;
 let f=det.filter(r=>(!ch||r.chrom==ch)&&(!tt||r.topology_type==tt)&&(!dd||r.determinacy==dd)&&(!gc||r.genome_ctx==gc)&&r.n_clusters>=mc&&(!fp||r.fp>0)&&(!q||r.region.includes(q)));
 f.sort(SORT[so]||SORT.coord);
 el('cnt').textContent=f.length+' 區';
 el('list').innerHTML=f.slice(0,700).map(r=>`<div class="row" data-i="${det.indexOf(r)}"><b>${r.region}</b> <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type.split('(')[0]}</span><span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span><br><span class="note">${r.n_sSNV}sSNV·c=${r.n_clusters}·${r.haplotypes}·${r.cn}·TP${r.tp}/FP${r.fp}${r.ambig_nodes>0?'·⚠序未定':''}</span></div>`).join('')+(f.length>700?`<div class="note" style="padding:8px">...前 700（共 ${f.length}）</div>`:'');
 el('list').querySelectorAll('.row').forEach(x=>x.onclick=()=>show(+x.dataset.i,x));
}
function show(i,row){el('list').querySelectorAll('.row').forEach(x=>x.classList.remove('sel'));if(row)row.classList.add('sel');let r=det[i];
 let popcount=r.populations;
 let pt=Object.entries(r.populations).sort((a,b)=>b[1]-a[1]).map(([g,c])=>{let tot=Object.values(r.populations).reduce((a,b)=>a+b,0);return `<tr><td class="mono">${g}</td><td>${sLabels(g)}</td><td>${c}</td><td>${(100*c/tot).toFixed(0)}%</td></tr>`}).join('');
 el('detail').innerHTML=`<h3>${r.region} <span class="tag ${TT[r.topology_type]||'t_single'}">${r.topology_type}</span> <span class="tag ctx_${r.genome_ctx}">${r.genome_ctx}</span></h3>
  <div class="kv"><div class="b">${r.n_sSNV} sSNV</div><div class="b">span ${(r.span/1000).toFixed(1)}kb</div><div class="b">c=${r.n_clusters} 群</div><div class="b">HP: ${r.haplotypes}</div><div class="b">CN: ${r.cn}</div><div class="b">TP ${r.tp} / FP ${r.fp}</div><div class="b">${r.determinacy}</div>${r.drop_noise_frac>0?`<div class="b">噪聲過濾 ${(r.drop_noise_frac*100).toFixed(0)}%</div>`:''}${r.ambig_nodes>0?`<div class="b" style="background:#fff3bf">⚠ 順序未定 ${r.ambig_nodes}(缺中間群)</div>`:''}</div>
  <b>克隆樹（germline→…；節點=S-mut-set·reads·%；座標=向量）</b>${tree(r.edges,r.populations,r.n_clusters,r.haplotypes)}
  <div class="note">S1..S${r.n_sSNV}=區內排序 sSNV；直系=往下、姊妹=同層分叉。tree_shape(pairwise)=${r.tree_shape}。genome_ctx 為近似(±3Mb)。</div>
  <b>細胞群(基因型向量 → S 突變 → reads → 佔比)</b><table><tr><th>向量</th><th>突變(S)</th><th>reads</th><th>佔比</th></tr>${pt}</table>`;
}
['f_chr','f_topology_type','f_determinacy','f_genome_ctx','f_minc','f_fp','f_q','f_sort'].forEach(id=>{let e=el(id);e.oninput=render;e.onchange=render});
render();
"""

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>克隆樹拓樸工作站 v2 — sSNV 重建（HCC1395 ⭐3）</title><style>{CSS}</style></head><body><div class="wrap">
<h1>克隆樹拓樸互動工作站 v2（cluster-first 算法 + S/r/m 標籤 + 篩選排序）</h1>
<p class="sub">HCC1395 ⭐3 · 每區 genotype 向量→拓樸(perfect-phylogeny+噪聲過濾) · S=sSNV/r=read群/m=甲基位點 · 數字由 JSON 注入 · {len(d['detail']):,} 區可畫樹</p>
{UNIVERSE_BANNER}
<div class="stats">
<div class="scard"><h4>拓樸型態</h4><div id="s_topo"></div></div><div class="scard"><h4>群數 c</h4><div id="s_clust"></div></div>
<div class="scard"><h4>determinacy</h4><div id="s_det"></div></div><div class="scard"><h4>HP 根數</h4><div id="s_root"></div></div>
</div>
<details class="c17"><summary>▶ chr17:48360161 完整 worked example（S/r/m 一致標籤 + 樹 + 甲基；點開）</summary><div id="c17"></div></details>
<div class="ctrl">
chr<select id="f_chr"><option value="">全</option></select>
拓樸<select id="f_topology_type"><option value="">全</option></select>
determinacy<select id="f_determinacy"><option value="">全</option></select>
位置<select id="f_genome_ctx"><option value="">全</option></select>
排序<select id="f_sort"><option value="coord">座標</option><option value="nsnv">複雜度(sSNV)</option><option value="nclust">群數</option><option value="region">region名</option></select>
最少群數<input id="f_minc" type="number" value="0" min="0" max="6" style="width:52px">
<label><input id="f_fp" type="checkbox">僅含FP</label>
搜尋<input id="f_q" placeholder="chr17:" style="width:120px">
<span id="cnt" class="note"></span>
</div>
<div class="main"><div class="list" id="list"></div><div class="detail" id="detail"><div class="note">← 左側點選一個區查看克隆樹（或點上方 chr17 worked example）</div></div></div>
<p class="note" style="margin-top:12px">⚠ 證據層級：A_determined=單分子向量；A_ambiguous=缺中間群順序未定；B_pairwise=拼接非單分子整樹；C_underdetermined=多樹相容。TP/FP=SEQC2 僅觀察不進前處理。genome_ctx 為近似(±3Mb)。甲基不參與拓樸裁決(cis-confounded)。⭐3 單樣本。</p>
</div>
<script>window.__DATA__={DJ};</script><script>{JS}</script></body></html>"""
with open(OUT, "w", encoding="utf-8") as f: f.write(HTML)
print(f"OK wrote {OUT} ({len(HTML):,} bytes; detail {len(d['detail'])})")
