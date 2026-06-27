#!/usr/bin/env python3
"""[全位點 FULL 儀表板 v2] C++ 原生全基因組 — 路徑式外部PNG + 豐富篩選(對齊 HCC1395 script102) + 門檻試算 + 無結構顯示.
功能: reading-guide SVG / provenance / 門檻滑桿即時試算 / 多重篩選+顯示無圖 / 逐項判讀+modal / 匯出JSON+CSV+進度 / 統計圖.
需與同目錄 figs/ 一起開. 資料 phylo_cpp_wg_records.json + summary.json + figs_cpp_wg/."""
import json, os
A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUTD = f"{A}/obs_ws/cpp_wg"; os.makedirs(OUTD, exist_ok=True)
rec = json.load(open(f"{A}/phylo_cpp_wg_full_records.json"))
summ = json.load(open(f"{A}/phylo_cpp_wg_full_summary.json")); tp, fp = summ["TP"], summ["FP"]
figlink = f"{OUTD}/figs"
if os.path.islink(figlink) or os.path.exists(figlink):
    try: os.remove(figlink)
    except Exception: pass
try: os.symlink(f"{A}/figs_cpp_wg_full", figlink)
except Exception: pass
figset = set(os.listdir(f"{A}/figs_cpp_wg_full")) if os.path.isdir(f"{A}/figs_cpp_wg_full") else set()
CHRS = [f"chr{c}" for c in range(1, 23)]; ck = {c: i for i, c in enumerate(CHRS)}


def figname(r):
    p = int(r["pos"]); fn = f"cpp_{r['chrom']}_{p}_{p+10000}.png"
    return fn if fn in figset else None


def cat(r):
    return 2 if r["coarse_ng"] < 2 else (0 if r["aligned"] else 1)


# [ck,pos(SNV),set,n,g,fine,other,unstable,aligned,hidden,hasfig,cat,vhp,valle]
D = []  # figmap 移除: figname 由 chrom+pos 決定式, JS 端即算(省 ~1.8MB)
for r in rec:
    fn = figname(r)
    D.append([ck.get(r["chrom"], 0), int(r["pos"]) + 5000, 0 if r["set"] == "TP" else 1, r["n"], r["coarse_ng"],
              r["fine_ng"], r["n_other"], 1 if r["unstable"] else 0, 1 if r["aligned"] else 0,
              1 if r["hidden_het"] else 0, 1 if fn else 0, cat(r), r.get("V_hp", 0), r.get("V_allele", 0)])
CATN = ["aligned cis-ASM", "unaligned subclone候選", "no_structure(無結構)"]
n_fig = sum(d[10] for d in D)

T = """<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>phylo-v4.1 C++ 全位點觀察分類儀表板 FULL</title><style>
:root{--bg:#0a0e16;--card:#111826;--fg:#e7edf5;--mut:#8b9bb4;--line:#2a3344;--ac:#D97757}
*{box-sizing:border-box}body{font-family:system-ui,"PingFang TC","Noto Sans CJK TC",sans-serif;background:var(--bg);color:var(--fg);margin:0;line-height:1.5}
.wrap{max-width:1340px;margin:0 auto;padding:18px 16px 90px}h1{font-size:19px;margin:4px 0}h2{font-size:14px;border-left:4px solid var(--ac);padding-left:8px;margin:20px 0 8px}
.sub{color:var(--mut);font-size:12px}.box,.panel{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:11px 14px;margin:9px 0}
.tab{display:inline-block;padding:6px 12px;border:1px solid var(--line);border-radius:7px;margin:2px;cursor:pointer;font-size:12.5px}.tab.on{background:var(--ac);color:#fff;border-color:var(--ac)}
select,input[type=number],input[type=range]{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:4px 7px;font-size:12px}
.row{display:flex;flex-wrap:wrap;gap:8px;align-items:center}.row label{font-size:12px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(310px,1fr));gap:10px;margin-top:8px}
.card{background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden}.card.j-in{border-color:#16a34a}.card.j-ex{border-color:#dc2626}.card.j-mb{border-color:#d97706}
.hd{padding:7px 10px;border-bottom:1px solid var(--line);cursor:pointer}.ttl{font-size:13px;font-weight:600}.meta{color:var(--mut);font-size:11px;margin-top:2px}
.figs{cursor:pointer;background:#0b1222}.figs img{width:100%;display:block;min-height:40px}.nofig{padding:14px;text-align:center;color:var(--mut);font-size:11.5px;background:#0b1222}
.badge{display:inline-block;padding:1px 6px;border-radius:4px;font-size:10.5px;margin-left:3px}.b-al{background:#0a2540;color:#7dd3fc}.b-un{background:#3b0764;color:#d8b4fe}.b-ns{background:#1e293b;color:#94a3b8}
.b-tp{background:#0a2a14;color:#4ade80}.b-fp{background:#2a0a0a;color:#f87171}.b-uns{background:#3a2a05;color:#fbbf24}.b-oth{background:#222;color:#aaa}.b-hh{background:#2a0a0a;color:#fca5a5}.b-fn{background:#06281e;color:#5eead4}
.jrow{display:flex;gap:4px;padding:7px 10px;border-top:1px solid var(--line)}.jb{flex:1;padding:5px;border-radius:5px;border:1px solid var(--line);background:#0b1222;color:var(--fg);cursor:pointer;font-size:11px}
.jb.on-in{background:#16a34a;color:#fff}.jb.on-ex{background:#dc2626;color:#fff}.jb.on-mb{background:#d97706;color:#fff}.jb.on-un{background:#475569;color:#fff}
button{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:6px 11px;cursor:pointer;font-size:12px}
#modal{position:fixed;inset:0;background:rgba(0,0,0,.88);display:none;z-index:99;overflow:auto;padding:20px}.mc{max-width:1080px;margin:0 auto;background:var(--card);border:1px solid var(--ac);border-radius:12px;padding:16px}.mc img{width:100%;border-radius:6px}
.statg{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:7px;margin:10px 0}.st{background:#0b1222;border:1px solid var(--line);border-radius:6px;padding:6px 9px}.st .l{color:var(--mut);font-size:10.5px}.st .n{font-weight:700;font-size:13px}
.sb{display:inline-block;background:#0b1222;border:1px solid var(--line);border-radius:7px;padding:6px 12px;margin:3px}.sb .n{font-weight:700;font-size:16px}.sb .l{color:var(--mut);font-size:10.5px}.kg{color:#4ade80}.kr{color:#f87171}.ko{color:#fbbf24}
table{border-collapse:collapse;font-size:11.5px}th,td{border:1px solid var(--line);padding:3px 7px;text-align:left}th{background:#0b1222}</style></head>
<body><div class="wrap"><h1>phylo-v4.1 C++ 原生全位點觀察分類儀表板 — FULL</h1>
<div class="sub">全 __NTOTAL__ 位點 ｜ binary 原生 phylo_groups.tsv → Python 只讀畫 ｜ __STAMP__ ｜ <b>需與同目錄 figs/ 一起開</b>（路徑式，全 __NFIG__ 位點圖含單群=可驗分錯；演化樹群色=tree-aware）</div>
<div class="box"><b>🔴 結論</b>：unaligned subclone候選 TP __TPU__% ≈ FP __FPU__%（非特異）→ 單樣本無法當 subclone。aligned cis-ASM TP __TPA__%/FP __FPA__%。</div>

<h2>⓪ 怎麼看（圖例）</h2><div class="panel"><div class="row" style="align-items:flex-start">
<svg width="430" height="150" viewBox="0 0 430 150" role="img"><text x="6" y="15" fill="#cbd5e1" font-size="12">每張圖: 左=UPGMA樹(視覺排序) 中=甲基read×CpG 右=read×read距離</text>
<rect x="6" y="26" width="13" height="13" fill="#db2777"/><text x="24" y="37" fill="#94a3b8" font-size="11">C++ coarse 群1</text>
<rect x="120" y="26" width="13" height="13" fill="#0d9488"/><text x="138" y="37" fill="#94a3b8" font-size="11">群2</text>
<rect x="190" y="26" width="13" height="13" fill="#555"/><text x="208" y="37" fill="#94a3b8" font-size="11">other殘群</text>
<rect x="280" y="26" width="13" height="13" fill="#dcdcdc"/><text x="298" y="37" fill="#94a3b8" font-size="11">離群</text>
<rect x="6" y="46" width="13" height="13" fill="#60a5fa"/><text x="24" y="57" fill="#94a3b8" font-size="11">HP1</text>
<rect x="46" y="46" width="13" height="13" fill="#1e3a8a"/><text x="64" y="57" fill="#94a3b8" font-size="11">HP1-1</text>
<rect x="104" y="46" width="13" height="13" fill="#c084fc"/><text x="122" y="57" fill="#94a3b8" font-size="11">HP2</text>
<rect x="146" y="46" width="13" height="13" fill="#6b21a8"/><text x="164" y="57" fill="#94a3b8" font-size="11">HP2-1</text>
<rect x="206" y="46" width="13" height="13" fill="#fbbf24"/><text x="224" y="57" fill="#94a3b8" font-size="11">REF</text>
<rect x="252" y="46" width="13" height="13" fill="#dc2626"/><text x="270" y="57" fill="#94a3b8" font-size="11">ALT</text>
<text x="6" y="78" fill="#cbd5e1" font-size="11">甲基色: <tspan fill="#dc2626">紅</tspan>=甲基 <tspan fill="#2563eb">藍</tspan>=未甲基 灰=未覆蓋 ｜橘線=變異位置</text>
<text x="6" y="98" fill="#cbd5e1" font-size="11">距離: 暗=近 亮=遠; <tspan fill="#86efac">對角暗塊且 coarse 側欄對齊 HP/ALT = 真分群</tspan></text>
<text x="6" y="122" fill="#cbd5e1" font-size="11">無結構: 整片均勻、無對角塊 = 單群(no_structure)</text>
<text x="6" y="142" fill="#fca5a5" font-size="11">unaligned(同 germline 多群)= subclone 候選, 但單樣本非特異</text></svg>
<div style="min-width:300px;flex:1"><table><tr><th>欄位</th><th>意義</th></tr>
<tr><td>coarse 群</td><td>modal K=10 多數決 robust 群數(主verdict)</td></tr>
<tr><td>fine 群</td><td>null90 候選更細群(低信心)</td></tr>
<tr><td>aligned</td><td>各群對齊 hp/allele(CramerV≥0.3)=cis-ASM</td></tr>
<tr><td>V_allele/V_hp</td><td>對齊 allele/hp 的 CramerV 強度</td></tr>
<tr><td>unstable</td><td>modal_frac&lt;0.7=真分裂,verdict 不穩</td></tr>
<tr><td>other</td><td>殘留離群≥3 記錄成 other 群</td></tr></table></div></div></div>

<h2>① Provenance</h2><div class="panel"><table>
<tr><th style="width:90px">樣本</th><td>HCC1395 tumor-only(配對 normal 做 HP);ONT 5khz simplex 5mCG+5hmCG;<b>單樣本</b>⭐2-3</td></tr>
<tr><th>突變來源</th><td>pileup somatic SNV;TP=SEQC2 高信度集內,FP=不在(caller誤判)</td></tr>
<tr><th>切群方法</th><td>C++ PhyloLabeler phylo-v4.1(modal K=10/fine/other);binary 原生 phylo_groups.tsv</td></tr>
<tr><th>範圍</th><td>chr1-22 全 somatic call;每位點=SNV±5000 視窗 tumor+HP-tagged read×CpG</td></tr>
<tr><th>數量</th><td>TP __TPN__ + FP __FPN__ = __NTOTAL__;全 __NFIG__ 位點圖(單群+多群,可驗分錯)</td></tr></table></div>

<h2>② 統計</h2><div id="stats"></div><div class="panel"><div class="sub" id="catstat"></div></div>

<h2>③ 即時門檻試算（拉滑桿看通過數）</h2><div class="panel"><div class="row">
 n(reads)≥<input type="range" id="tn" min="6" max="200" value="6" style="width:140px"><span id="tnv">6</span>
 coarse群≥<input type="range" id="tg" min="1" max="6" value="2" style="width:90px"><span id="tgv">2</span>
 V_allele≥<input type="range" id="tv" min="0" max="1" step="0.05" value="0" style="width:120px"><span id="tvv">0</span>
</div><div id="thlive" style="margin-top:8px"></div></div>

<h2>④ 分類觀察（點 tab；可顯示無結構）</h2><div id="cattabs"></div>
<div class="panel"><div class="row">
 set <select id="fset"><option value="">全</option><option value="0">TP</option><option value="1">FP</option></select>
 群數<input id="fgmin" type="number" value="0" min="0" max="9" style="width:46px">–<input id="fgmax" type="number" value="9" min="0" max="9" style="width:46px">
 排序<select id="sort"><option value="3">reads</option><option value="4">coarse群</option><option value="13">V_allele</option><option value="1">位置</option></select>
 <span id="sdir" style="cursor:pointer;color:var(--ac)">▼降序</span>
 <label><input type="checkbox" id="ff">僅有圖</label>
 <label><input type="checkbox" id="fu">unstable</label>
 <label><input type="checkbox" id="fh">hidden_het</label>
 <label><input type="checkbox" id="fo">有other</label>
 <label><input type="checkbox" id="ffc">fine&gt;coarse</label>
 判讀<select id="jf"><option value="all">全</option><option value="in">應包含</option><option value="mb">可能漏</option><option value="ex">排除</option><option value="un">未判</option></select>
 <span class="sub" id="gcount"></span></div></div>
<div class="grid" id="grid"></div><div id="pager" style="text-align:center;margin-top:12px"></div>

<h2>⑤ 判讀紀錄（匯出/匯入/進度，無結構也可記）</h2><div class="panel"><div class="row">
 <button onclick="exportJ()">⬇ 匯出 JSON</button><button onclick="exportCSV()">⬇ 匯出 CSV</button>
 <label class="sub">⬆ 匯入 <input type="file" id="imp" accept=".json"></label>
 <button onclick="if(confirm('清空所有判讀?'))clearJ()">✕ 清空</button><span class="sub" id="impmsg"></span></div>
 <div id="jprog" style="margin-top:10px"></div></div>
</div><div id="modal" onclick="if(event.target.id==='modal')closeM()"><div class="mc" id="mc"></div></div>
<script>
const D=__D__, CHRS=__CHRS__, CATN=__CATN__, TPN=__TPN__,FPN=__FPN__,NFIG=__NFIG__;
const C={ck:0,ps:1,set:2,n:3,g:4,fine:5,oth:6,uns:7,al:8,hh:9,fig:10,cat:11,vhp:12,va:13};
const LS='phylo_cpp_full_v2';let J={};try{J=JSON.parse(localStorage.getItem(LS)||'{}')}catch(e){}
const $=id=>document.getElementById(id);const key=r=>CHRS[r[C.ck]]+':'+r[C.ps];
const figName=r=>'cpp_'+CHRS[r[C.ck]]+'_'+(r[C.ps]-5000)+'_'+(r[C.ps]+5000)+'.png';
function setJ(k,c){J[k]=J[k]||{};J[k].c=c;J[k].t=new Date().toISOString();localStorage.setItem(LS,JSON.stringify(J));draw();prog();}
$('stats').innerHTML=[['TP',TPN.toLocaleString(),''],['FP',FPN.toLocaleString(),''],
 ['aligned cis-ASM(TP)',D.filter(r=>r[C.set]===0&&r[C.cat]===0).length.toLocaleString(),'kg'],
 ['unaligned候選(TP)',D.filter(r=>r[C.set]===0&&r[C.cat]===1).length.toLocaleString(),'kr'],
 ['no_structure(TP)',D.filter(r=>r[C.set]===0&&r[C.cat]===2).length.toLocaleString(),'ko'],
 ['圖數',NFIG.toLocaleString(),'']].map(s=>`<div class="sb"><div class="n ${s[2]}">${s[1]}</div><div class="l">${s[0]}</div></div>`).join('');
$('catstat').innerHTML='分類(全 set): '+CATN.map((n,i)=>`${n} <b>${D.filter(r=>r[C.cat]===i).length.toLocaleString()}</b>`).join(' ｜ ');
// 門檻試算
function thpass(r){return r[C.n]>=+$('tn').value&&r[C.g]>=+$('tg').value&&r[C.va]>=+$('tv').value;}
function thlive(){const tp=D.filter(r=>r[C.set]===0),fp=D.filter(r=>r[C.set]===1);
 const ptp=tp.filter(thpass).length,pfp=fp.filter(thpass).length;
 $('tnv').textContent=$('tn').value;$('tgv').textContent=$('tg').value;$('tvv').textContent=(+$('tv').value).toFixed(2);
 const enr=(ptp/tp.length)/((pfp/fp.length)||1e-9);
 $('thlive').innerHTML=[['通過 TP',ptp.toLocaleString(),'kg'],['通過 FP',pfp.toLocaleString(),'kr'],['TP富集',enr.toFixed(2)+'×',enr>1.2?'kg':enr<0.85?'kr':''],['TP通過率',(100*ptp/tp.length).toFixed(1)+'%','']].map(s=>`<div class="sb"><div class="n ${s[2]}">${s[1]}</div><div class="l">${s[0]}</div></div>`).join('');}
['tn','tg','tv'].forEach(id=>$(id).addEventListener('input',thlive));
// tabs+filter
let curCat=1,page=0,PER=24,sdir=-1;
$('cattabs').innerHTML=CATN.map((n,i)=>`<div class="tab${i===1?' on':''}" data-c="${i}">${n} <span id="cc${i}"></span></div>`).join('');
document.querySelectorAll('.tab').forEach(t=>t.onclick=()=>{curCat=+t.dataset.c;page=0;document.querySelectorAll('.tab').forEach(x=>x.classList.toggle('on',x===t));draw();});
$('sdir').onclick=()=>{sdir*=-1;$('sdir').textContent=sdir===-1?'▼降序':'▲升序';page=0;draw();};
['fset','fgmin','fgmax','sort','ff','fu','fh','fo','ffc','jf'].forEach(id=>$(id).addEventListener('input',()=>{page=0;draw();}));
function filtered(){let a=D.filter(r=>r[C.cat]===curCat);
 const fs=$('fset').value;if(fs!=='')a=a.filter(r=>r[C.set]===+fs);
 a=a.filter(r=>r[C.g]>=+$('fgmin').value&&r[C.g]<=+$('fgmax').value);
 if($('ff').checked)a=a.filter(r=>r[C.fig]);if($('fu').checked)a=a.filter(r=>r[C.uns]);if($('fh').checked)a=a.filter(r=>r[C.hh]);
 if($('fo').checked)a=a.filter(r=>r[C.oth]);if($('ffc').checked)a=a.filter(r=>r[C.fine]>r[C.g]);
 const jf=$('jf').value;if(jf!=='all')a=a.filter(r=>{const c=(J[key(r)]||{}).c||'un';return jf==='un'?c==='un':c===jf;});
 const sk=+$('sort').value;a.sort((x,y)=>sdir*((sk===1?x[C.ck]*1e10+x[C.ps]:x[sk])-(sk===1?y[C.ck]*1e10+y[C.ps]:y[sk])));
 return a;}
function badges(r){let b=['<span class="badge b-al">aligned</span>','<span class="badge b-un">unaligned</span>','<span class="badge b-ns">no-struct</span>'][r[C.cat]];
 b+=r[C.set]===0?'<span class="badge b-tp">TP</span>':'<span class="badge b-fp">FP</span>';
 if(r[C.uns])b+='<span class="badge b-uns">unstable</span>';if(r[C.oth])b+='<span class="badge b-oth">other'+r[C.oth]+'</span>';
 if(r[C.hh])b+='<span class="badge b-hh">hidden-het</span>';if(r[C.fine]>r[C.g])b+='<span class="badge b-fn">fine'+r[C.fine]+'</span>';return b;}
function jbtn(k,c,lab){const on=((J[k]||{}).c===c)?(' on-'+c):'';return `<button class="jb${on}" onclick="event.stopPropagation();setJ('${k}','${c}')">${lab}</button>`;}
function card(r){const k=key(r),fn=r[C.fig]?figName(r):null;
 const fig=fn?`<div class="figs" onclick="openM('${k}')"><img loading="lazy" src="figs/${fn}"></div>`:'<div class="nofig">無圖（reads 過少）｜ n='+r[C.n]+' coarse '+r[C.g]+'群</div>';
 const cc=(J[k]||{}).c;const cls=cc==='in'?' j-in':cc==='ex'?' j-ex':cc==='mb'?' j-mb':'';
 return `<div class="card${cls}"><div class="hd" onclick="openM('${k}')"><div class="ttl">${k} ${badges(r)}</div>
  <div class="meta">n=${r[C.n]} ｜ coarse ${r[C.g]}群 ｜ fine ${r[C.fine]} ｜ V_al ${r[C.va]} ｜ V_hp ${r[C.vhp]}</div></div>${fig}
  <div class="jrow">${jbtn(k,'in','應包含')}${jbtn(k,'mb','可能漏')}${jbtn(k,'ex','排除')}${jbtn(k,'un','無判')}</div></div>`;}
function draw(){const a=filtered();
 for(let i=0;i<3;i++)$('cc'+i).textContent='('+D.filter(r=>r[C.cat]===i).length.toLocaleString()+')';
 $('gcount').textContent=a.length.toLocaleString()+' 個 (有圖 '+a.filter(r=>r[C.fig]).length+')';
 const pages=Math.max(1,Math.ceil(a.length/PER));page=Math.min(page,pages-1);
 $('grid').innerHTML=a.slice(page*PER,page*PER+PER).map(card).join('')||'<div class="sub">無符合</div>';
 $('pager').innerHTML=`<button ${page<=0?'disabled':''} onclick="page--;draw()">‹</button> 第 ${page+1}/${pages} 頁 <button ${page>=pages-1?'disabled':''} onclick="page++;draw()">›</button>`;}
window.openM=k=>{const r=D.find(x=>key(x)===k);if(!r)return;const fn=r[C.fig]?figName(r):null;
 const fig=fn?`<img src="figs/${fn}">`:'<div class="nofig">無圖（此位點 reads 過少未渲染）</div>';
 const st=[['分類',CATN[r[C.cat]]],['set',r[C.set]===0?'TP':'FP'],['n reads',r[C.n]],['coarse群',r[C.g]],['fine群',r[C.fine]],['other',r[C.oth]],['V_allele',r[C.va]],['V_hp',r[C.vhp]],['unstable',r[C.uns]?'是':'否'],['aligned',r[C.al]?'是(cis-ASM)':'否'],['hidden-het',r[C.hh]?'是':'否']];
 $('mc').innerHTML=`<div style="display:flex;justify-content:space-between"><h2 style="border:none;margin:0">${k} ${badges(r)}</h2><button onclick="closeM()">關閉✕</button></div>
  ${fig}<div class="statg">${st.map(s=>`<div class="st"><div class="l">${s[0]}</div><div class="n">${s[1]}</div></div>`).join('')}</div>
  <div class="jrow" style="margin-top:8px">${jbtn(k,'in','應包含')}${jbtn(k,'mb','可能漏')}${jbtn(k,'ex','排除')}${jbtn(k,'un','無判')}</div>
  <div class="sub" style="margin-top:8px">距離對角塊乾淨且 coarse 側欄對齊 HP/ALT=真分群。無塊整片均勻=無結構。unaligned 同germline多群=subclone候選但單樣本非特異。</div>`;
 $('modal').style.display='block';};
window.closeM=()=>$('modal').style.display='none';document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});
function prog(){const tot=D.length;const j=Object.values(J).filter(v=>v.c&&v.c!=='un').length;
 const ci={in:0,mb:0,ex:0};Object.values(J).forEach(v=>{if(ci[v.c]!==undefined)ci[v.c]++;});
 $('jprog').innerHTML=[['已判讀',j.toLocaleString(),''],['應包含',ci.in,'kg'],['可能漏',ci.mb,'ko'],['排除',ci.ex,'kr'],['進度',(100*j/tot).toFixed(2)+'%','']].map(s=>`<div class="sb"><div class="n ${s[2]}">${s[1]}</div><div class="l">${s[0]}</div></div>`).join('');}
window.exportJ=()=>{const o={generated:new Date().toISOString(),judgments:Object.entries(J).filter(([k,v])=>v.c&&v.c!=='un').map(([k,v])=>({locus:k,choice:v.c,ts:v.t}))};dl('phylo_cpp_judgments.json',JSON.stringify(o,null,1),'application/json');};
window.exportCSV=()=>{const rows=[['locus','chr','pos','set','category','coarse_ng','fine_ng','aligned','unstable','choice']];
 D.forEach(r=>{const k=key(r),c=(J[k]||{}).c;if(c&&c!=='un')rows.push([k,CHRS[r[C.ck]],r[C.ps],r[C.set]===0?'TP':'FP',CATN[r[C.cat]],r[C.g],r[C.fine],r[C.al],r[C.uns],c]);});
 dl('phylo_cpp_judgments.csv',rows.map(x=>x.join(',')).join('\\n'),'text/csv');};
function dl(name,txt,mime){const b=new Blob([txt],{type:mime});const u=URL.createObjectURL(b);const a=document.createElement('a');a.href=u;a.download=name;a.click();}
window.clearJ=()=>{J={};localStorage.removeItem(LS);draw();prog();};
$('imp').onchange=e=>{const f=e.target.files[0];if(!f)return;const rd=new FileReader();rd.onload=()=>{try{const o=JSON.parse(rd.result);(o.judgments||[]).forEach(j=>{J[j.locus]={c:j.choice,t:j.ts}});localStorage.setItem(LS,JSON.stringify(J));$('impmsg').textContent='匯入 '+(o.judgments||[]).length+' 筆';draw();prog();}catch(x){$('impmsg').textContent='匯入失敗'}};rd.readAsText(f);};
thlive();draw();prog();
</script></body></html>"""
T = (T.replace("__D__", json.dumps(D, separators=(",", ":"))).replace("__CHRS__", json.dumps(CHRS))
     .replace("__CATN__", json.dumps(CATN, ensure_ascii=False))
     .replace("__TPN__", str(tp["n"])).replace("__FPN__", str(fp["n"])).replace("__NTOTAL__", f"{len(D):,}")
     .replace("__NFIG__", str(n_fig)).replace("__STAMP__", "2026-06-23 C++ native")
     .replace("__TPU__", str(tp["unaligned_pct"])).replace("__FPU__", str(fp["unaligned_pct"]))
     .replace("__TPA__", str(tp["aligned_pct"])).replace("__FPA__", str(fp["aligned_pct"])))
out = f"{OUTD}/20260623_phylo_cpp_observation_dashboard_FULL.html"
open(out, "w").write(T)
print(f"WROTE {out} ({os.path.getsize(out)//1024} KB) | {len(D)} loci | {n_fig} figs (含 150 無結構樣本)")
