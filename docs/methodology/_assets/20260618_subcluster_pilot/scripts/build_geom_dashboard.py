#!/usr/bin/env python3
"""[geometry-vs-null 觀察工作站 generator] geom_records.json + counts + summary → standalone HTML。
數字全由 JSON 注入(§13-A,不手打)。path-based figure lazy-load。逐位點判讀 localStorage + JSON/CSV 匯出。
banner: 觀察層,非 subclone caller / 非 TP/FP filter。⚠ geometry=純幾何無 null,本就過切。
template = 純字串 + token replace(避免 f-string 與內嵌 JS 大括號衝突)。"""
import json

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
recs = json.load(open(f"{A}/geom_records.json"))
summ = json.load(open(f"{A}/geom_summary.json"))
counts = json.load(open(f"{A}/phylo_weak_structure_counts.json"))
idx = json.load(open(f"{A}/geom_render_index.json"))

WS = {"S1_CONFIRMED_MULTI_ALIGNED": "S1·確認多群對齊", "S2_CONFIRMED_MULTI_UNALIGNED": "S2·確認多群不對齊",
      "S3_UNSTABLE_MULTI": "S3·不穩定多群", "S4_WEAK_FINE_ONLY": "S4·次閾值候選", "S6_NO_STRUCTURE_CLEAN": "S6·乾淨單群"}
D = []
for r in recs:
    fig = idx.get(f"{r['chrom']}:{r['rec_pos']}", "")
    D.append([r["chrom"], r["rec_pos"] + 5000, r["set"], WS.get(r["wstate"], r["wstate"]), r["n"],
              r["coarse_ng"], r["fine_ng"], r["geometry_ng"], 1 if r["divergence"] else 0,
              r["V_hp"], r["V_allele"], fig])

states = []
for s in ["S1_CONFIRMED_MULTI_ALIGNED", "S2_CONFIRMED_MULTI_UNALIGNED", "S3_UNSTABLE_MULTI", "S4_WEAK_FINE_ONLY", "S6_NO_STRUCTURE_CLEAN"]:
    tp = counts["TP"]["states"][s]; fp = counts["FP"]["states"][s]
    states.append([WS.get(s, s), tp["n"], tp["pct"], fp["n"], fp["pct"]])

DATA = {"D": D, "states": states, "xtab": summ["coarse_vs_geometry_xtab"], "caveats": counts["caveats"],
        "summ": {k: summ[k] for k in ["n_sampled", "n_regions_rerun", "n_divergence", "divergence_pct_of_rerun", "divergence_by_wstate", "geom_ge2_total", "geom_gt_coarse_total"]},
        "wg": {"n_total": counts["n_total"], "tp_n": counts["TP"]["n"], "fp_n": counts["FP"]["n"]}}

TPL = r"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>幾何-vs-null 弱結構觀察工作站 · HCC1395</title>
<style>
:root{--bg:#0f1115;--card:#1a1d24;--bd:#2a2f3a;--tx:#e6e8ec;--mut:#9aa3b2;--acc:#d97757}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--tx);font:14px/1.55 system-ui,'Noto Sans CJK TC',sans-serif}
.wrap{max-width:1500px;margin:0 auto;padding:18px}
h1{font-size:20px;margin:0 0 4px}.sub{color:var(--mut);font-size:13px;margin-bottom:12px}
.banner{background:#3a2418;border:1px solid var(--acc);border-radius:8px;padding:10px 14px;margin:10px 0;font-size:13px}
.banner b{color:var(--acc)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:12px 0}
.kpi{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px 12px}
.kpi .v{font-size:22px;font-weight:700}.kpi .l{color:var(--mut);font-size:12px}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin:6px 0}
th,td{border:1px solid var(--bd);padding:5px 8px;text-align:right}th{background:#20242d;text-align:center}td.l{text-align:left}
.flt{display:flex;flex-wrap:wrap;gap:8px;align-items:center;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px;margin:12px 0;position:sticky;top:0;z-index:5}
.flt label{font-size:12px;color:var(--mut)}select,input{background:#0f1115;color:var(--tx);border:1px solid var(--bd);border-radius:5px;padding:4px 7px;font-size:12.5px}
.btn{background:#20242d;color:var(--tx);border:1px solid var(--bd);border-radius:5px;padding:5px 10px;cursor:pointer;font-size:12.5px}.btn:hover{border-color:var(--acc)}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px;margin:10px 0}
.card.div{border-color:var(--acc)}
.hd{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:6px}
.bd{font-size:11px;padding:2px 7px;border-radius:10px;font-weight:600}
.b-div{background:#3a2418;color:var(--acc);border:1px solid var(--acc)}.b-tp{background:#10243a;color:#60a5fa}.b-fp{background:#3a1a1a;color:#f87171}
.b-st{background:#20242d;color:var(--mut);border:1px solid var(--bd)}
.m{color:var(--mut);font-size:12px}.m b{color:var(--tx)}
img{width:100%;border-radius:5px;background:#fff;display:block;margin:4px 0}
.jd{display:flex;gap:6px;flex-wrap:wrap;margin-top:6px}
.jb{padding:4px 9px;border-radius:5px;border:1px solid var(--bd);background:#0f1115;color:var(--tx);cursor:pointer;font-size:12px}
.jb.on-a{background:#14532d;border-color:#16a34a}.jb.on-d{background:#3a2f10;border-color:#d97706}.jb.on-r{background:#3a1414;border-color:#dc2626}
.rsn{margin-top:5px;width:100%}.prog{color:var(--mut);font-size:12px}
details{margin:8px 0}summary{cursor:pointer;color:var(--acc)}
.ft{color:var(--mut);font-size:11.5px;border-top:1px solid var(--bd);margin-top:18px;padding-top:10px}
</style></head><body><div class="wrap">
<h1>幾何-vs-null 弱結構觀察工作站 — HCC1395</h1>
<div class="sub">「scipy 0.7×max 幾何平切看得到、但 phylo null 閘擋掉」的弱結構池 — 子集重跑 @@NREG@@ 位點,逐項可觀察判讀</div>
<div class="banner">⚠ <b>這是觀察/刻畫層,不是 subclone caller、不是 TP/FP filter。</b>
<b>幾何欄(geometry_ng)= 純 scipy 0.7×max 平切,無 null,本就會過切單一 clone</b>(gap/silhouette 已證 1 clone 被判 k=3-4);它是「null 閘保守度」的診斷,<b>非 subclone 宣稱</b>。
confident-multi 在 FP 比 TP 還富集(反判別,=cis-ASM 非 subclone)。單樣本 HCC1395 探索級 ⭐2-3,無外部真值。</div>
<div id="kpis" class="grid"></div>
<details open><summary>全基因組 5-態分類 + 子集 coarse×geometry 交叉表</summary>
<div style="display:flex;gap:18px;flex-wrap:wrap"><div style="flex:1;min-width:340px"><b style="font-size:12px">全基因組 5-態(N=@@NTOTAL@@)</b><table id="st"></table></div>
<div style="flex:1;min-width:300px"><b style="font-size:12px">子集 coarse(null95)×geometry 交叉表(橘=分歧)</b><table id="xt"></table></div></div></details>
<div class="flt">
<label>set <select id="fset"><option value="">全</option><option>TP</option><option>FP</option></select></label>
<label>狀態 <select id="fst"><option value="">全</option></select></label>
<label>分歧 <select id="fdv"><option value="">全</option><option value="1">★只看分歧</option><option value="0">非分歧</option></select></label>
<label>幾何≥ <input id="fg" type="number" value="0" min="0" style="width:55px"></label>
<label>排序 <select id="fsort"><option value="geom">幾何群數↓</option><option value="div">分歧優先</option><option value="n">reads↓</option></select></label>
<span class="prog" id="prog"></span>
<button class="btn" id="ej">匯出判讀JSON</button><button class="btn" id="ec">匯出CSV</button>
</div>
<div id="list"></div>
<div class="ft" id="ft"></div>
</div>
<script>
const DATA=@@DATA@@;
const LS="geomdiv_judge_v1";
const J=DATA.D.map(r=>({chrom:r[0],pos:r[1],set:r[2],st:r[3],n:r[4],coarse:r[5],fine:r[6],geom:r[7],div:r[8],vhp:r[9],val:r[10],fig:r[11]}));
let judge=JSON.parse(localStorage.getItem(LS)||"{}");
function kid(r){return r.set+":"+r.chrom+":"+r.pos}
const s=DATA.summ;
document.getElementById('kpis').innerHTML=[
 ['重跑位點',s.n_regions_rerun],['★分歧(幾何≥2但null=1)',s.n_divergence],['分歧率',s.divergence_pct_of_rerun+'%'],
 ['幾何≥2 總',s.geom_ge2_total],['幾何>null 總',s.geom_gt_coarse_total]
].map(k=>`<div class="kpi"><div class="v">${k[1]}</div><div class="l">${k[0]}</div></div>`).join('');
document.getElementById('st').innerHTML='<tr><th>狀態</th><th>TP n</th><th>TP%</th><th>FP n</th><th>FP%</th></tr>'+
 DATA.states.map(r=>`<tr><td class="l">${r[0]}</td><td>${r[1]}</td><td>${r[2]}</td><td>${r[3]}</td><td>${r[4]}</td></tr>`).join('');
const xe=Object.entries(DATA.xtab).map(([k,v])=>{const m=k.match(/coarse(\d+)_geom(\d+)/);return{c:+m[1],g:+m[2],n:v}});
let xh='<tr><th>null95\\geom</th>';for(let g=1;g<=5;g++)xh+=`<th>geom${g}${g==5?'+':''}</th>`;xh+='</tr>';
for(let c=1;c<=3;c++){xh+=`<tr><td class="l">coarse${c}${c==3?'+':''}</td>`;for(let g=1;g<=5;g++){const e=xe.find(x=>x.c==c&&x.g==g);const dv=(c==1&&g>=2);xh+=`<td style="${dv?'color:var(--acc);font-weight:700':''}">${e?e.n:0}</td>`}xh+='</tr>'}
document.getElementById('xt').innerHTML=xh;
[...new Set(J.map(r=>r.st))].sort().forEach(st=>{const o=document.createElement('option');o.value=st;o.textContent=st;document.getElementById('fst').appendChild(o)});
const fsetEl=document.getElementById('fset'),fstEl=document.getElementById('fst'),fdvEl=document.getElementById('fdv'),fgEl=document.getElementById('fg'),fsortEl=document.getElementById('fsort'),listEl=document.getElementById('list');
function render(){
 const fset=fsetEl.value,fst=fstEl.value,fdv=fdvEl.value,fg=+fgEl.value,sort=fsortEl.value;
 let rows=J.filter(r=>(!fset||r.set==fset)&&(!fst||r.st==fst)&&(fdv===''||r.div==+fdv)&&r.geom>=fg);
 rows.sort((a,b)=>sort=='geom'?b.geom-a.geom:sort=='n'?b.n-a.n:(b.div-a.div)||(b.geom-a.geom));
 document.getElementById('prog').textContent=`顯示 ${rows.length}/${J.length} · 已判讀 ${Object.keys(judge).length}`;
 listEl.innerHTML=rows.map(r=>{const id=kid(r),jv=judge[id]||{};
  return `<div class="card ${r.div?'div':''}">
   <div class="hd">
    ${r.div?'<span class="bd b-div">★分歧</span>':''}
    <span class="bd ${r.set=='TP'?'b-tp':'b-fp'}">${r.set}</span>
    <span class="bd b-st">${r.st}</span>
    <b>${r.chrom}:${r.pos}</b>
    <span class="m">n=<b>${r.n}</b> · null95=<b>${r.coarse}</b> null90=<b>${r.fine}</b> 幾何=<b style="color:var(--acc)">${r.geom}</b> · V_hp=${r.vhp} V_al=${r.val}</span>
   </div>
   ${r.fig?`<img loading="lazy" src="${r.fig}" alt="${r.chrom}:${r.pos}">`:'<div class="m">(無圖)</div>'}
   <div class="jd">
    <button class="jb ${jv.v=='a'?'on-a':''}" onclick="setj('${id}','a')">同意=真結構</button>
    <button class="jb ${jv.v=='d'?'on-d':''}" onclick="setj('${id}','d')">存疑</button>
    <button class="jb ${jv.v=='r'?'on-r':''}" onclick="setj('${id}','r')">否定=幾何過切/cis-ASM</button>
   </div>
   <input class="rsn" placeholder="判讀理由(選填)" value="${(jv.r||'').replace(/"/g,'&quot;')}" onchange="setr('${id}',this.value)">
  </div>`}).join('')||'<div class="m">無符合條件位點</div>';
}
function setj(id,v){judge[id]=judge[id]||{};judge[id].v=(judge[id].v==v?null:v);save()}
function setr(id,r){judge[id]=judge[id]||{};judge[id].r=r;save(true)}
function save(noRender){localStorage.setItem(LS,JSON.stringify(judge));if(!noRender)render()}
[fsetEl,fstEl,fdvEl,fgEl,fsortEl].forEach(e=>e.addEventListener('input',render));
document.getElementById('ej').onclick=()=>{const b=new Blob([JSON.stringify(judge,null,1)],{type:'application/json'});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='geomdiv_judgments.json';a.click()};
document.getElementById('ec').onclick=()=>{let c='locus,set,state,n,coarse,fine,geom,divergence,verdict,reason\n';J.forEach(r=>{const j=judge[kid(r)]||{};c+=`${r.chrom}:${r.pos},${r.set},${r.st},${r.n},${r.coarse},${r.fine},${r.geom},${r.div},${j.v||''},"${(j.r||'').replace(/"/g,'""')}"\n`});const b=new Blob([c],{type:'text/csv'});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='geomdiv_judgments.csv';a.click()};
document.getElementById('ft').innerHTML=`來源: geom_records.json (子集重跑 N=${s.n_regions_rerun}) + phylo_weak_structure_counts.json (全基因組 N=${DATA.wg.n_total}=TP ${DATA.wg.tp_n}+FP ${DATA.wg.fp_n}) · 單樣本 HCC1395 探索 ⭐2-3 · 幾何=純0.7×max平切無null · `+Object.values(DATA.caveats).join(' / ');
render();
</script></body></html>"""

out = (TPL.replace("@@DATA@@", json.dumps(DATA, ensure_ascii=False))
          .replace("@@NREG@@", str(DATA["summ"]["n_regions_rerun"]))
          .replace("@@NTOTAL@@", str(DATA["wg"]["n_total"])))
with open(f"{A}/20260623_geometry_divergence_observation_dashboard.standalone.html", "w") as f:
    f.write(out)
print("wrote 20260623_geometry_divergence_observation_dashboard.standalone.html")
print(f"loci={len(D)} divergence={summ['n_divergence']} figs={len(idx)} html_bytes={len(out)}")
