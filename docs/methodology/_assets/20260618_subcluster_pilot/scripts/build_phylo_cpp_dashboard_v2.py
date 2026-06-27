#!/usr/bin/env python3
"""[全位點觀察 dashboard v2] 本 session full run + SEQC2 CN/LOH 標註。
新增: SEQC2 CN 預設標籤 / 對齊只敘述不歸類(預設) / 歸類+可能分群數量(展開內文) / 位點搜尋 /
CN狀態+CN值範圍+LOH 篩選 / 多排序 / 整體統計摺疊。圖名 JS 端決定式(省 FIGMAP)。數字全由 JSON 注入。"""
import json, os
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot".replace("/InterSubMod/docs", "/InterSubMod/.claude/worktrees/ism-review-infra/docs")
OUTD = f"{A}/obs_ws/cpp_wg"; os.makedirs(OUTD, exist_ok=True)
rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_full.json"))  # +n_tumor/n_normal/geometry_ng
cnstats = json.load(open(f"{A}/cn_stats.json"))
tng = json.load(open(f"{A}/tn_geometry_stats.json"))
summ = json.load(open(f"{A}/phylo_cpp_wg_full_summary.json"))
CHRS = [f"chr{c}" for c in range(1, 23)]; ck = {c: i for i, c in enumerate(CHRS)}
CNS = ["gain", "loss", "loh", "neutral"]; cni = {s: i for i, s in enumerate(CNS)}

# 確保 figs symlink → figs_cpp_wg_full
figlink = f"{OUTD}/figs"
if not (os.path.islink(figlink) and os.path.realpath(figlink) == f"{A}/figs_cpp_wg_full"):
    try:
        if os.path.islink(figlink) or os.path.exists(figlink): os.remove(figlink)
        os.symlink(f"{A}/figs_cpp_wg_full", figlink)
    except Exception: pass
figset = set(os.listdir(f"{A}/figs_cpp_wg_full")) if os.path.isdir(f"{A}/figs_cpp_wg_full") else set()

# 緊湊 D: [ck, ps(SNV), set(0TP/1FP), n, coarse, fine, n_other, unstable, hidden, vhp, val, cns(idx), cnv(-1=NA), loh, hasfig]
D = []
for r in rec:
    p = int(r["pos"]); snv = p + 5000
    fn = f"cpp_{r['chrom']}_{p}_{p+10000}.png"
    cv = r.get("cn_value")
    D.append([ck.get(r["chrom"], 0), snv, 0 if r["set"] == "TP" else 1, r["n"], r["coarse_ng"], r["fine_ng"],
              r["n_other"], 1 if r["unstable"] else 0, 1 if r["hidden_het"] else 0,
              r.get("V_hp", 0), r.get("V_allele", 0), cni.get(r.get("cn_state", "neutral"), 3),
              (cv if isinstance(cv, (int, float)) else -1), 1 if r.get("is_loh") else 0, 1 if fn in figset else 0,
              r.get("n_tumor") if r.get("n_tumor") is not None else -1,
              r.get("n_normal") if r.get("n_normal") is not None else -1,
              r.get("geometry_ng") if r.get("geometry_ng") is not None else -1])

# 整體統計(全 JSON 注入,不手打)
def dist(key):
    return dict(Counter(r[key] for r in rec))
STATS = {
    "n": len(rec), "tp": summ["TP"]["n"], "fp": summ["FP"]["n"],
    "structure_tp": summ["TP"]["structure_pct"], "no_structure_tp": summ["TP"]["no_structure_pct"],
    "aligned_tp": summ["TP"]["aligned_pct"], "unaligned_tp": summ["TP"]["unaligned_pct"],
    "aligned_fp": summ["FP"]["aligned_pct"], "unaligned_fp": summ["FP"]["unaligned_pct"],
    "cn_state_dist": cnstats["cn_state_dist"], "cn_state_x_structure": cnstats["cn_state_x_structure"],
    "cn_value_dist": {str(k): v for k, v in cnstats["cn_value_dist"].items()},
    "coarse_dist_tp": {str(k): v for k, v in summ["TP"]["ngroups_dist"].items()},
    "coarse_dist_fp": {str(k): v for k, v in summ["FP"]["ngroups_dist"].items()},
    "unstable_tp": summ["TP"]["unstable_pct"], "fine_multi_tp": summ["TP"]["fine_multi_pct"],
    "geom_divergence": tng.get("geometry_divergence(geom≥2 but coarse<2)", 0),
    "median_nt": tng.get("median_n_tumor"), "median_nn": tng.get("median_n_normal"),
    "geom_ng_dist": {str(k): v for k, v in tng.get("geom_ng_dist", {}).items()},
}

TPL = r"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>phylo 全位點觀察 v2 (+SEQC2 CN)</title><style>
:root{--bg:#0a0e16;--card:#111826;--fg:#e7edf5;--mut:#8b9bb4;--line:#2a3344;--ac:#D97757;--gain:#ef4444;--loss:#3b82f6;--loh:#a855f7;--neu:#64748b}
*{box-sizing:border-box}body{font-family:system-ui,"Noto Sans CJK TC",sans-serif;background:var(--bg);color:var(--fg);margin:0;line-height:1.55}
.wrap{max-width:1480px;margin:0 auto;padding:16px}
h1{font-size:19px;margin:2px 0}.sub{color:var(--mut);font-size:12.5px}
details{background:var(--card);border:1px solid var(--line);border-radius:8px;margin:8px 0;padding:4px 12px}
summary{cursor:pointer;font-weight:600;font-size:13.5px;padding:6px 0}
.statg{display:grid;grid-template-columns:repeat(auto-fill,minmax(135px,1fr));gap:7px;margin:8px 0}
.st{background:#0b1222;border:1px solid var(--line);border-radius:6px;padding:6px 9px}.st .n{font-weight:700;font-size:14px}.st .l{color:var(--mut);font-size:10.5px}
table{border-collapse:collapse;font-size:12px;margin:4px 0}td,th{border:1px solid var(--line);padding:3px 8px;text-align:right}th{background:#1a2130}td.l{text-align:left}
.bar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;background:var(--card);border:1px solid var(--line);border-radius:8px;padding:9px;margin:9px 0;position:sticky;top:0;z-index:9}
.bar label{font-size:11.5px;color:var(--mut)}input,select{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:5px;padding:4px 7px;font-size:12px}
input[type=search]{width:160px}.btn{background:#1a2130;border:1px solid var(--line);border-radius:5px;color:var(--fg);padding:5px 9px;cursor:pointer;font-size:12px}.btn:hover{border-color:var(--ac)}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(330px,1fr));gap:11px;margin:10px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden}
.hd{padding:8px 10px}.ttl{font-weight:700;font-size:13px}.badges{margin-top:3px;display:flex;flex-wrap:wrap;gap:4px}
.badge{font-size:10.5px;padding:1px 6px;border-radius:9px;border:1px solid var(--line)}
.b-tp{background:#10243a;color:#7cc4ff}.b-fp{background:#3a1a1a;color:#ffa3a3}
.b-gain{background:#3a1414;color:#ffb3b3;border-color:var(--gain)}.b-loss{background:#10233a;color:#9fc6ff;border-color:var(--loss)}.b-loh{background:#241038;color:#d6b3ff;border-color:var(--loh)}.b-neu{background:#1a2130;color:var(--mut)}
.b-al{background:#13261c;color:#9ff0c4}.b-tn{background:#0f2a24;color:#7fe3c4}.b-div{background:#3a2418;color:#ffba80;border-color:#d97706}.desc{color:var(--mut);font-size:11px;margin-top:3px}
.figs{cursor:pointer;background:#0b1222;border-top:1px solid var(--line)}.figs img{width:100%;display:block}.nofig{padding:18px;text-align:center;color:var(--mut);font-size:11px}
.exp{padding:7px 10px;border-top:1px solid var(--line);font-size:11.5px}.exp b{color:var(--fg)}.exp .k{color:var(--mut)}
.jrow{display:flex;gap:5px;padding:7px 10px;border-top:1px solid var(--line)}
.jb{flex:1;padding:4px;border:1px solid var(--line);border-radius:5px;background:#0b1222;color:var(--fg);cursor:pointer;font-size:11px}.jb.on-in{background:#14532d}.jb.on-mb{background:#3a2f10}.jb.on-ex{background:#3a1414}
.card.j-in{border-color:#16a34a}.card.j-mb{border-color:#d97706}.card.j-ex{border-color:#dc2626}
.modal{display:none;position:fixed;inset:0;background:#000a;z-index:50;overflow:auto;padding:20px}.mc{max-width:1200px;margin:0 auto;background:var(--card);border:1px solid var(--ac);border-radius:10px;padding:14px}
.mc img{width:100%;border-radius:6px;background:#fff}
.pager{text-align:center;margin:12px 0}.pager button{padding:5px 12px}
.note{background:#2a2410;border:1px solid #d97706;border-radius:7px;padding:8px 12px;font-size:12px;margin:8px 0}.note b{color:#fbbf24}
</style></head><body><div class="wrap">
<h1>phylo-v4.1 全位點觀察 v2 — +SEQC2 CN/LOH</h1>
<div class="sub">本 session full run（演化樹 tree-aware 群色）+ SEQC2 NGS-benchmark CN 真值標註。圖路徑式，需與同層 figs/ 一起開。__STAMP__</div>
<div class="note">⚠ <b>對齊（V_hp/V_allele）只敘述、不當歸類</b>：對齊 germline/allele 軸 = cis-ASM 跡象，<b>非 subclone</b>。歸類與分群數量見每卡展開。confident-multi 在 FP 比 TP 多=反判別。單樣本 HCC1395 探索 ⭐2-3。</div>

<details open><summary>📊 整體統計觀察（點擊摺疊）</summary><div id="stats"></div></details>

<div class="bar">
<label>🔍 <input type="search" id="q" placeholder="搜尋 chr1:39078224 或 chr8"></label>
<label>set <select id="fset"><option value="">全</option><option value="0">TP</option><option value="1">FP</option></select></label>
<label>CN狀態 <select id="fcn"><option value="">全</option><option value="0">gain</option><option value="1">loss</option><option value="2">loh</option><option value="3">neutral</option></select></label>
<label>CN值 <input type="number" id="cnmin" placeholder="min" style="width:60px" step="0.5"> - <input type="number" id="cnmax" placeholder="max" style="width:60px" step="0.5"></label>
<label><input type="checkbox" id="floh"> 只 LOH</label>
<label>coarse群 <input type="number" id="gmin" value="1" style="width:48px"> - <input type="number" id="gmax" value="9" style="width:48px"></label>
<label><input type="checkbox" id="funs"> unstable</label>
<label><input type="checkbox" id="ffc"> fine&gt;coarse</label>
<label><input type="checkbox" id="fdiv"> 只可疑漏切(幾何≥2嚴閘1)</label>
<label>排序 <select id="sort"><option value="1">位置</option><option value="3">n reads</option><option value="4">coarse群</option><option value="12">CN值</option><option value="vmax">對齊V</option></select> <button class="btn" id="sdir">▼</button></label>
<span class="sub" id="cnt"></span>
<button class="btn" id="ej">匯出判讀</button>
</div>
<div id="grid" class="grid"></div>
<div class="pager" id="pager"></div>
</div>
<div class="modal" id="modal" onclick="if(event.target===this)closeM()"><div class="mc" id="mc"></div></div>
<script>
const D=__D__, CHRS=__CHRS__, STATS=__STATS__, CNS=['gain','loss','loh','neutral'];
const C={ck:0,ps:1,set:2,n:3,cg:4,fg:5,oth:6,uns:7,hh:8,vhp:9,val:10,cns:11,cnv:12,loh:13,fig:14,nt:15,nn:16,geom:17};
const $=id=>document.getElementById(id);
const figName=r=>'cpp_'+CHRS[r[C.ck]]+'_'+(r[C.ps]-5000)+'_'+(r[C.ps]+5000)+'.png';
const key=r=>CHRS[r[C.ck]]+':'+r[C.ps];
let J=JSON.parse(localStorage.getItem('phylo_v2_judge')||'{}'),page=0,sdir=-1,PER=24;

// ---- 整體統計 ----
(function(){const s=STATS;
 const cnColor={gain:'var(--gain)',loss:'var(--loss)',loh:'var(--loh)',neutral:'var(--neu)'};
 let h=`<div class="statg">
  <div class="st"><div class="n">${s.n.toLocaleString()}</div><div class="l">總位點 (TP ${s.tp.toLocaleString()}+FP ${s.fp.toLocaleString()})</div></div>
  <div class="st"><div class="n">${s.structure_tp}%</div><div class="l">TP 有結構(coarse≥2)</div></div>
  <div class="st"><div class="n">${s.aligned_tp}%</div><div class="l">TP 對齊(cis-ASM跡象)</div></div>
  <div class="st"><div class="n">${s.unaligned_tp}% / FP ${s.unaligned_fp}%</div><div class="l">unaligned 候選 TP/FP</div></div>
  <div class="st"><div class="n">${s.unstable_tp}%</div><div class="l">TP unstable</div></div>
  <div class="st"><div class="n">${s.fine_multi_tp}%</div><div class="l">TP fine&gt;coarse</div></div>
  <div class="st"><div class="n">${s.geom_divergence.toLocaleString()}</div><div class="l">可疑漏切(幾何≥2但嚴閘1群)</div></div>
  <div class="st"><div class="n">${s.median_nt}/${s.median_nn}</div><div class="l">中位 tumor/normal reads</div></div></div>`;
 h+='<b style="font-size:12px">SEQC2 CN 狀態分佈</b><table><tr><th class="l">CN狀態</th><th>位點</th><th>%</th><th>有結構(multi)</th><th>單群</th></tr>';
 CNS.forEach(st=>{const tot=s.cn_state_dist[st]||0;const m=s.cn_state_x_structure[st+'_multi']||0;const sg=s.cn_state_x_structure[st+'_single']||0;
  h+=`<tr><td class="l"><span style="color:${cnColor[st]}">●</span> ${st}</td><td>${tot.toLocaleString()}</td><td>${(100*tot/s.n).toFixed(1)}</td><td>${m.toLocaleString()}</td><td>${sg.toLocaleString()}</td></tr>`;});
 h+='</table><b style="font-size:12px">coarse 群數分佈 (TP)</b><table><tr><th class="l">群數</th>'+Object.keys(s.coarse_dist_tp).sort((a,b)=>a-b).map(k=>`<th>${k}</th>`).join('')+'</tr><tr><td class="l">TP</td>'+Object.keys(s.coarse_dist_tp).sort((a,b)=>a-b).map(k=>`<td>${s.coarse_dist_tp[k]}</td>`).join('')+'</tr></table>';
 h+='<b style="font-size:12px">CN 值分佈</b>: '+Object.keys(s.cn_value_dist).sort((a,b)=>a-b).map(k=>`CN${k} <b>${s.cn_value_dist[k]}</b>`).join(' ｜ ');
 h+='<br><b style="font-size:12px">幾何切群數分佈 (0.7cut, 全位點)</b>: '+Object.keys(s.geom_ng_dist).sort((a,b)=>a-b).map(k=>`${k}群 <b>${s.geom_ng_dist[k]}</b>`).join(' ｜ ')+' <span class="sub">（幾何=純0.7×max平切，本就過切，僅供「肉眼/null 對照」非 subclone）</span>';
 $('stats').innerHTML=h;})();

// ---- 篩選 ----
function filt(){let a=D.slice();
 const q=$('q').value.trim().toLowerCase();
 if(q){a=a.filter(r=>(CHRS[r[C.ck]]+':'+r[C.ps]).toLowerCase().includes(q)||CHRS[r[C.ck]].toLowerCase()===q);}
 const fs=$('fset').value;if(fs!=='')a=a.filter(r=>r[C.set]==fs);
 const fc=$('fcn').value;if(fc!=='')a=a.filter(r=>r[C.cns]==fc);
 const cmin=$('cnmin').value,cmax=$('cnmax').value;
 if(cmin!=='')a=a.filter(r=>r[C.cnv]>=0&&r[C.cnv]>=+cmin);if(cmax!=='')a=a.filter(r=>r[C.cnv]>=0&&r[C.cnv]<=+cmax);
 if($('floh').checked)a=a.filter(r=>r[C.loh]);
 a=a.filter(r=>r[C.cg]>=+$('gmin').value&&r[C.cg]<=+$('gmax').value);
 if($('funs').checked)a=a.filter(r=>r[C.uns]);if($('ffc').checked)a=a.filter(r=>r[C.fg]>r[C.cg]);
 if($('fdiv').checked)a=a.filter(r=>r[C.geom]>=2&&r[C.cg]<2);
 const sk=$('sort').value;
 a.sort((x,y)=>{const f=v=>sk==='vmax'?Math.max(v[C.vhp],v[C.val]):(sk==1?v[C.ck]*1e10+v[C.ps]:v[+sk]);return sdir*(f(x)-f(y));});
 return a;}

function cnBadge(r){const st=CNS[r[C.cns]];const cv=r[C.cnv]>=0?(' CN'+r[C.cnv]):'';const cls={0:'b-gain',1:'b-loss',2:'b-loh',3:'b-neu'}[r[C.cns]];
 return `<span class="badge ${cls}">${st==='loh'?'LOH':st+cv}</span>`;}
function alignDesc(r){const vh=r[C.vhp],va=r[C.val];const mx=Math.max(vh,va);
 return mx>=0.3?`對齊 ${va>=vh?'allele':'hp'} 軸 (V_hp ${vh} / V_al ${va}) — cis-ASM 跡象，非 subclone`:`未對齊 germline 軸 (V_hp ${vh} / V_al ${va})`;}
function structText(r){const cg=r[C.cg],fg=r[C.fg];
 let cls=cg<2?'乾淨單群（無監督切不出≥2）':(r[C.uns]?'不穩定多群（seed 分歧，borderline）':(Math.max(r[C.vhp],r[C.val])>=0.3?'確認多群·對齊 germline（cis-ASM 候選，非 subclone）':'確認多群·不對齊（結構但無標籤對應）'));
 return cls;}
function cgBadge(r){const g=r[C.geom],cg=r[C.cg],div=(g>=2&&cg<2);return `<span class="badge ${div?'b-div':'b-al'}">切群 ${cg}/${r[C.fg]}/${g>=0?g:'?'}${div?' ⚠漏?':''}</span>`;}
function tnBadge(r){return `<span class="badge b-tn">T${r[C.nt]>=0?r[C.nt]:'?'}/N${r[C.nn]>=0?r[C.nn]:'?'}</span>`;}
function clusterText(r){const cg=r[C.cg],fg=r[C.fg],g=r[C.geom];const v=(g>=2&&cg<2)?`<span style="color:#ffba80">幾何看到 ${g} 群但嚴閘判單群 = 可疑漏切</span>`:(fg>cg?`寬閘 null90 多切（${cg}→${fg}）`:'各法一致');return `null95(嚴/主判) <b>${cg}</b> ｜ null90(寬/候選) <b>${fg}</b> ｜ 幾何0.7cut(肉眼) <b>${g>=0?g:'?'}</b> → ${v}`;}

function card(r){const k=key(r),fn=r[C.fig]?figName(r):null;const j=(J[k]||{}).c;const jc=j?(' j-'+j):'';
 const fig=fn?`<div class="figs" onclick="openM('${k}')"><img loading="lazy" src="figs/${fn}"></div>`:'<div class="nofig">無圖（reads 過少）</div>';
 return `<div class="card${jc}"><div class="hd"><div class="ttl" onclick="openM('${k}')">${k}</div>
  <div class="badges">${r[C.set]===0?'<span class="badge b-tp">TP</span>':'<span class="badge b-fp">FP</span>'}${cnBadge(r)}${tnBadge(r)}${cgBadge(r)}${r[C.loh]?'<span class="badge b-loh">LOH</span>':''}</div>
  <div class="desc">${alignDesc(r)}</div></div>
  ${fig}
  <details class="exp"><summary>歸類 + 切群數 + tumor/normal</summary>
   <div><span class="k">歸類:</span> <b>${structText(r)}</b></div>
   <div><span class="k">切群數:</span> ${clusterText(r)}</div>
   <div><span class="k">tumor/normal:</span> tumor <b>${r[C.nt]>=0?r[C.nt]:'?'}</b> / normal <b>${r[C.nn]>=0?r[C.nn]:'?'}</b> reads（phylo 用 tumor ${r[C.n]} 條切群）</div>
   <div><span class="k">SEQC2 CN:</span> ${CNS[r[C.cns]]}${r[C.cnv]>=0?' CN='+r[C.cnv]:''}${r[C.loh]?' (LOH)':''} ｜ other ${r[C.oth]} ｜ unstable ${r[C.uns]?'是':'否'} ｜ hidden-het ${r[C.hh]?'是':'否'}</div></details>
  <div class="jrow"><button class="jb${j==='in'?' on-in':''}" onclick="setJ('${k}','in')">應含</button><button class="jb${j==='mb'?' on-mb':''}" onclick="setJ('${k}','mb')">可能漏</button><button class="jb${j==='ex'?' on-ex':''}" onclick="setJ('${k}','ex')">排除</button></div></div>`;}

function draw(){const a=filt();const pages=Math.max(1,Math.ceil(a.length/PER));page=Math.min(page,pages-1);
 $('cnt').textContent=a.length.toLocaleString()+' 個';
 $('grid').innerHTML=a.slice(page*PER,page*PER+PER).map(card).join('')||'<div class="sub">無符合</div>';
 $('pager').innerHTML=`<button class="btn" ${page<=0?'disabled':''} onclick="page--;draw()">‹上頁</button> 第 ${page+1}/${pages} 頁 <button class="btn" ${page>=pages-1?'disabled':''} onclick="page++;draw()">下頁›</button>`;}
['q','fset','fcn','cnmin','cnmax','floh','gmin','gmax','funs','ffc','fdiv','sort'].forEach(id=>$(id).addEventListener('input',()=>{page=0;draw();}));
$('sdir').onclick=()=>{sdir*=-1;$('sdir').textContent=sdir===-1?'▼':'▲';page=0;draw();};
window.setJ=(k,c)=>{J[k]=J[k]||{};J[k].c=(J[k].c===c?null:c);localStorage.setItem('phylo_v2_judge',JSON.stringify(J));draw();};
window.openM=k=>{const r=D.find(x=>key(x)===k);if(!r)return;const fn=r[C.fig]?figName(r):null;
 $('mc').innerHTML=`<div style="display:flex;justify-content:space-between"><h2 style="margin:0">${k}</h2><button class="btn" onclick="closeM()">關閉 ✕</button></div>
  <div class="desc" style="margin:6px 0">${alignDesc(r)}</div>
  ${fn?`<img src="figs/${fn}">`:'<div class="nofig">無圖</div>'}
  <div class="exp" style="border:none"><div><span class="k">歸類:</span> <b>${structText(r)}</b></div>
   <div><span class="k">切群數:</span> ${clusterText(r)}</div>
   <div><span class="k">tumor/normal:</span> tumor <b>${r[C.nt]>=0?r[C.nt]:'?'}</b> / normal <b>${r[C.nn]>=0?r[C.nn]:'?'}</b>（phylo 用 tumor ${r[C.n]} 條）</div>
   <div><span class="k">SEQC2 CN:</span> ${CNS[r[C.cns]]}${r[C.cnv]>=0?' CN='+r[C.cnv]:''}${r[C.loh]?' LOH':''} ｜ set ${r[C.set]===0?'TP':'FP'}</div>
   <div><span class="k">對齊:</span> V_hp ${r[C.vhp]} / V_allele ${r[C.val]} ｜ unstable ${r[C.uns]?'是':'否'} ｜ hidden-het ${r[C.hh]?'是':'否'} ｜ other ${r[C.oth]}</div></div>
  <div class="sub" style="margin-top:8px">距離對角塊乾淨+coarse 側欄對齊 HP/ALT=真分群。整片均勻=無結構。對齊 germline=cis-ASM 非 subclone。</div>`;
 $('modal').style.display='block';};
window.closeM=()=>$('modal').style.display='none';
document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});
$('ej').onclick=()=>{const o=Object.entries(J).filter(([k,v])=>v.c).map(([k,v])=>({locus:k,choice:v.c}));
 const b=new Blob([JSON.stringify({n:o.length,judgments:o},null,1)],{type:'application/json'});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='phylo_v2_judgments.json';a.click();};
draw();
</script></body></html>"""

out = (TPL.replace("__D__", json.dumps(D, separators=(",", ":")))
          .replace("__CHRS__", json.dumps(CHRS))
          .replace("__STATS__", json.dumps(STATS, ensure_ascii=False))
          .replace("__STAMP__", "2026-06-23 · SEQC2 CNV ngs_benchmark(hg38) · single-sample ⭐2-3"))
outp = f"{OUTD}/20260623_phylo_cpp_observation_dashboard_v2.html"
open(outp, "w").write(out)
print(f"WROTE {outp} ({len(out)//1024} KB) | {len(D)} loci")
