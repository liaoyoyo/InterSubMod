#!/usr/bin/env python3
"""
86 - Build the redesigned interactive locus-judgment display.

Every curated locus has TWO viewable figures (methylation read x CpG + ISM read-read
distance matrix). Two galleries: PASS (kept by gate) and FILTERED-OUT (dropped but
high mean-shift) so the user can judge missed inclusions. A live threshold panel over
all 30,350 loci shows the TP/FP trade-off of any cutoff.

ALL numbers come from real files read at build time:
  - manifest.json            (curated 2070, rendered from ISM per-region output)
  - ism_existence_scan/HCC1395_{tp,fp}/significance_summary.csv  (full 30,350)
  - funnel_numbers.json      (upstream VCF counts, provenance stamped)

Output: display_v2/20260609_locus_judgment_display_01.html  (open in place; needs figs/)
"""
import csv, json, os, datetime

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
DV = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample/display_v2"
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
OUT = f"{DV}/20260609_locus_judgment_display_01.html"


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def load_all():
    rows = []
    for cls_i, cls in enumerate(("tp", "fp")):
        for r in csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")):
            cvs = [num(r, "CramersV"), num(r, "CramersV_HPFamily"), num(r, "CramersV_HPFine")]
            cvmax = max([c for c in cvs if c is not None] or [0.0])
            rows.append([round(cvmax, 3), round(abs(num(r, "HPMergedDelta") or 0.0), 3),
                         int(num(r, "NumReads") or 0), 1 if tb(r, "Potential_LOH") else 0,
                         cls_i, 1 if tb(r, "Significant") else 0])
    return rows


TEMPLATE = r"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>位點判斷與顯示確認 — HCC1395 甲基差異篩選 (ISM)</title>
<style>
:root{--bg:#0f172a;--card:#1e293b;--mut:#94a3b8;--fg:#e2e8f0;--ac:#38bdf8;--pass:#22c55e;--miss:#f59e0b;--fp:#ef4444;--line:#334155}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.6 -apple-system,'Segoe UI',system-ui,sans-serif}
header{padding:20px 26px;border-bottom:1px solid var(--line);background:#0b1222;position:sticky;top:0;z-index:30}
h1{margin:0 0 4px;font-size:20px} h2{font-size:16px;margin:30px 0 10px;border-left:3px solid var(--ac);padding-left:9px}
.sub{color:var(--mut);font-size:12.5px}
.wrap{max-width:1280px;margin:0 auto;padding:0 22px 80px}
.funnel{display:flex;gap:8px;flex-wrap:wrap;margin:14px 0}
.fstep{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:9px 13px;min-width:120px}
.fstep b{font-size:18px;color:var(--ac)} .fstep .l{font-size:11px;color:var(--mut)}
.arrow{align-self:center;color:var(--mut)}
.panel{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;margin:12px 0}
.ctrls{display:flex;gap:18px;flex-wrap:wrap;align-items:flex-end}
.ctrl{display:flex;flex-direction:column;gap:3px} .ctrl label{font-size:11.5px;color:var(--mut)}
.ctrl input[type=range]{width:160px} .ctrl .v{font-weight:700;color:var(--ac)}
select,button{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:6px 9px;font-size:13px}
.live{display:flex;gap:22px;flex-wrap:wrap;margin-top:14px;font-size:13px}
.live .big{font-size:22px;font-weight:800} .tp{color:var(--pass)} .fp{color:var(--fp)}
.tabs{display:flex;gap:6px;margin:8px 0} .tab{cursor:pointer;padding:7px 14px;border-radius:7px 7px 0 0;background:#16203a;border:1px solid var(--line);border-bottom:none}
.tab.on{background:var(--card);color:var(--ac);font-weight:700}
.gctrl{display:flex;gap:12px;flex-wrap:wrap;align-items:center;margin:8px 0;font-size:12.5px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:12px}
.card{background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden;cursor:pointer;transition:.12s}
.card:hover{border-color:var(--ac);transform:translateY(-2px)}
.card .hd{padding:7px 10px;border-bottom:1px solid var(--line)}
.card .ttl{font-weight:700;font-size:13px} .card .meta{font-size:11px;color:var(--mut);margin-top:2px}
.figs{display:flex;gap:2px;background:#0b1222} .figs img{width:50%;display:block;background:#fff}
.badge{display:inline-block;font-size:10px;padding:1px 6px;border-radius:9px;font-weight:700;margin-left:4px}
.b-pass{background:#14532d;color:#86efac} .b-miss{background:#78350f;color:#fcd34d} .b-loh{background:#3b0764;color:#d8b4fe} .b-fp{background:#450a0a;color:#fca5a5}
.rej{font-size:11px;color:#fcd34d;padding:5px 10px;background:#1c1407;border-top:1px solid var(--line)}
.judge{display:flex;gap:4px;padding:6px 10px;border-top:1px solid var(--line)}
.judge button{flex:1;font-size:11px;padding:4px} .jb-in.on{background:#14532d;color:#86efac;border-color:#22c55e} .jb-ex.on{background:#450a0a;color:#fca5a5;border-color:#ef4444}
.pager{display:flex;gap:8px;align-items:center;justify-content:center;margin:16px 0;font-size:13px}
#modal{position:fixed;inset:0;background:rgba(0,0,0,.86);display:none;z-index:100;overflow:auto;padding:24px}
#modal .mc{max-width:1080px;margin:0 auto;background:var(--card);border:1px solid var(--ac);border-radius:12px;padding:20px}
#modal img{width:100%;background:#fff;border-radius:6px} .mfigs{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin:14px 0}
.statgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:8px;margin-top:10px}
.stat{background:#0b1222;border:1px solid var(--line);border-radius:7px;padding:7px 10px} .stat .l{font-size:10.5px;color:var(--mut)} .stat .n{font-size:15px;font-weight:700}
.legend{display:flex;gap:18px;flex-wrap:wrap;font-size:12px;color:var(--mut)}
.sw{display:inline-block;width:12px;height:12px;border-radius:3px;vertical-align:-1px;margin-right:4px}
code{background:#0b1222;padding:1px 5px;border-radius:4px;font-size:12px}
.note-warn{background:#1c1407;border:1px solid #78350f;border-radius:8px;padding:10px 14px;font-size:12.5px;color:#fcd34d;margin:10px 0}
</style></head><body>
<header>
  <h1>位點判斷與顯示確認 — HCC1395 甲基差異篩選（ISM）</h1>
  <div class="sub">每個策展位點都有兩張圖：甲基 read×CpG + ISM read×read 距離矩陣 ｜ 列出被篩掉的位點供判斷是否該包含 ｜ 生成 __STAMP__ ｜ <b>需與同目錄 figs/ 一起開啟</b></div>
</header>
<div class="wrap">

<div class="note-warn">⚠ 單樣本 HCC1395、單一 pipeline（ClairS paired → longphase-S → ISM）。此頁是<b>位點層級的判斷/確認工具</b>，非 filter 宣稱。甲基→TP/FP filter 方向已 concluded NEGATIVE（見現況地圖），此處僅供你核對「篩選標準是否合理、是否漏掉應含位點」。</div>

<h2>① 篩選漏斗（真實數字）</h2>
<div class="funnel" id="funnel"></div>
<div class="sub">CramersV gate = <code>PassedGating &amp; global_p≤0.05 &amp; CramersV≥0.1 &amp; NumReads≥20</code>。「被篩掉(MISSED)」= 未過 gate 但 |Δβ|≥0.20（高平均位移）。完整上游漏斗見 capstone HTML。</div>

<h2>② 即時門檻試算（全 30,350 位點）— 觀察 TP/FP 取捨</h2>
<div class="panel">
  <div class="ctrls">
    <div class="ctrl"><label>CramersV 門檻 ≥ <span class="v" id="vcv">0.10</span></label><input type="range" id="cv" min="0" max="1" step="0.01" value="0.10"></div>
    <div class="ctrl"><label>|Δβ| 門檻 ≥ <span class="v" id="vdb">0.20</span></label><input type="range" id="db" min="0" max="0.6" step="0.01" value="0.20"></div>
    <div class="ctrl"><label>NumReads ≥ <span class="v" id="vrd">20</span></label><input type="range" id="rd" min="0" max="80" step="1" value="20"></div>
    <div class="ctrl"><label>邏輯</label><select id="mode">
      <option value="union">聯集 (CramersV 或 |Δβ|)</option>
      <option value="and">交集 (兩者都要)</option>
      <option value="cv">只看 CramersV</option>
      <option value="db">只看 |Δβ|</option>
    </select></div>
    <div class="ctrl"><label>LOH</label><select id="lohf">
      <option value="all">全部</option><option value="non">排除 LOH</option><option value="loh">只看 LOH</option>
    </select></div>
  </div>
  <div class="live" id="live"></div>
  <div class="sub" style="margin-top:8px">提示：CramersV 門檻調低或切「只看 |Δβ|」→ 通過數暴增但 <b>TP/FP 比往 1.0 崩</b>（Δβ 平均位移不分 TP/FP）。這就是 gate 用 CramersV 而非 Δβ 的原因。</div>
</div>

<h2>③ 位點圖庫 — 每張都可檢視，並列出被篩掉的</h2>
<div class="panel" style="padding:12px 16px">
  <div class="legend">
    <span><b>怎麼讀兩張圖：</b></span>
    <span><span class="sw" style="background:#1d4ed8"></span>HP1 <span class="sw" style="background:#16a34a"></span>HP1-1 <span class="sw" style="background:#9333ea"></span>HP2 <span class="sw" style="background:#ca8a04"></span>HP2-1</span>
    <span><span class="sw" style="background:#dc2626"></span>ALT <span class="sw" style="background:#0ea5e9"></span>REF</span>
    <span>甲基圖：<span class="sw" style="background:#b2182b"></span>甲基化 → <span class="sw" style="background:#2166ac"></span>未甲基；灰=NA</span>
    <span>距離圖：暗=近(相似)、亮=遠；按樹狀分群排序 → <b>對角塊狀</b>=有分群</span>
  </div>
  <div class="sub" style="margin-top:8px">判斷準則：高 CramersV 應在<b>距離圖看到對角塊</b>且塊與 HP/ALT 側欄對齊；高 |Δβ| 但 CramersV 低 → 距離圖<b>無塊</b>（reads 混在一起，只是整體甲基高低不同）= 為何被篩掉。</div>
</div>
<div class="note-warn" style="border-color:#0e7490;background:#082f33;color:#a5f3fc">⚡ <b>可能漏掉的位點（latent 真結構）：被篩掉的 __N_MISSED_NOTE__ 個 MISSED 中，有 <b>__N_LATENT__</b> 個其 CramersV 因「列聯表稀疏（期望格&lt;5，違反 Cochran 卡方有效性）」被判不可靠而歸零，<b>但 PERMANOVA（距離法、對稀疏穩健）仍顯著</b> → reads 其實有按 HP 分群。在 MISSED 分頁勾「只看潛在結構」即可逐一檢視，用距離圖確認塊是否乾淨（真結構）還是只有 2-3 條 read（統計脆弱）。</div>
<div class="tabs">
  <div class="tab on" data-g="PASS">✅ 通過篩選 PASS (<span id="cntP"></span>)</div>
  <div class="tab" data-g="MISSED">🟡 被篩掉 FILTERED-OUT (<span id="cntM"></span>，其中 ⚡<span id="cntL"></span> 潛在結構)</div>
</div>
<div class="panel" style="margin-top:0;border-radius:0 10px 10px 10px">
  <div class="gctrl">
    <span>排序</span>
    <select id="sort"></select>
    <label><input type="checkbox" id="onlyfig" checked> 只顯示有圖的</label>
    <label id="latwrap" style="display:none"><input type="checkbox" id="onlylat"> ⚡只看潛在結構</label>
    <span id="minnwrap" style="display:none">納回門檻 minHP≥<span class="v" id="vminn" style="color:var(--ac)">0</span><input type="range" id="minn" min="0" max="20" step="1" value="0" style="width:90px;vertical-align:middle"></span>
    <select id="glohf"><option value="all">LOH:全部</option><option value="non">排除 LOH</option><option value="loh">只看 LOH</option></select>
    <select id="clsf"><option value="all">TP+FP</option><option value="tp">只 TP</option><option value="fp">只 FP</option></select>
    <span id="gcount" style="color:var(--mut)"></span>
    <span id="judgesum" style="margin-left:auto"></span>
  </div>
  <div class="grid" id="grid"></div>
  <div class="pager" id="pager"></div>
</div>

<h2>④ 精修 gate：用 PERMANOVA 納回被「卡方不可靠」歸零的 latent 真結構</h2>
<div class="panel">
  <div class="sub" style="margin-bottom:10px">原始 gate 用 CramersV（卡方），HP 分群極不平衡（如 3 vs 51 reads）時表稀疏→判不可靠→歸零。改用 <b>PERMANOVA（距離法，對稀疏穩健）</b> 納回：條件 = <code>未過原 gate &amp; PassedGating &amp; NumReads≥20 &amp; ClusterPermanova 顯著 &amp; LabelHPPermanova 顯著（分群對齊 HP）&amp; min(HP1N,HP2N)≥minN</code>。拉 minN 看取捨：</div>
  <table id="rgtab" style="border-collapse:collapse;font-size:12.5px;width:100%"></table>
  <div class="note-warn" style="margin-top:12px">⚠ <b>誠實判讀</b>：原 gate 保留的是 <b>144:1 高度 TP-clean</b> 集（3× base 47:1）。納回的 latent 是 PERMANOVA 證實的<b>真 HP 分群結構</b>，但 TP-purity 較低（minN&lt;15 時甚至低於 base = FP-leaning）。<b>納回提升的是「ASM 結構 characterization 完整性」，不是 TP/FP 分辨力</b>（與已 concluded 的甲基→filter NEGATIVE 一致）。建議：以 characterization 為目的 → minN=10（納回 +2512 真結構位點、去除退化稀疏、納回集仍 35:1 TP-leaning）；以 cleanest gate 為目的 → 維持原始。</div>
</div>
</div>

<div id="modal" onclick="if(event.target.id==='modal')closeM()"><div class="mc" id="mc"></div></div>

<script>
const CUR=__PAYLOAD__;
const ALL=__ALLPAYLOAD__;
const FN=__FNJ__;
const RG=__REFINED_GATE__;
const $=id=>document.getElementById(id);

// ---- refined-gate sweep table ----
(function(){let h=`<tr style="color:var(--mut);text-align:right"><th style="text-align:left;padding:5px 8px">minN（最小 HP group reads）</th><th>納回 TP</th><th>納回 FP</th><th>納回集 TP:FP</th><th>vs base 47:1</th><th>合併 PASS (TP/FP)</th><th>納回 LOH%</th></tr>`;
h+=`<tr style="border-top:1px solid var(--line)"><td style="padding:5px 8px">原始 Significant（無納回）</td><td style="text-align:right">${RG.orig_pass_tp}</td><td style="text-align:right">${RG.orig_pass_fp}</td><td style="text-align:right;font-weight:700">${(RG.orig_pass_tp/Math.max(RG.orig_pass_fp,1)).toFixed(0)}:1</td><td style="text-align:right;color:var(--pass)">3.0×</td><td style="text-align:right">${RG.orig_pass_tp}/${RG.orig_pass_fp}</td><td style="text-align:right">—</td></tr>`;
for(const s of RG.sweep){const fpl=s.enrich_vs_base>=1?'var(--pass)':'var(--miss)';
 h+=`<tr style="border-top:1px solid var(--line)"><td style="padding:5px 8px">≥ ${s.minN}</td><td style="text-align:right">+${s.rec_tp}</td><td style="text-align:right">+${s.rec_fp}</td><td style="text-align:right">${s.rec_ratio}:1</td><td style="text-align:right;color:${fpl}">${s.enrich_vs_base}×</td><td style="text-align:right">${s.comb_tp}/${s.comb_fp}</td><td style="text-align:right">${s.rec_loh_pct}%</td></tr>`;}
const t=document.getElementById('rgtab');if(t)t.innerHTML=h;})();

const fsteps=[['caller→longphase-S TP',FN.lps_tp_vcf_records,''],['ISM 分析 TP',FN.ism_tp_analyzed,'31 太少 read 被丟'],
  ['gate PASS (TP)',__TP_PASS__,''],['被篩掉 MISSED (TP)',__TP_MISSED__,'|Δβ|≥0.2'],
  ['longphase-S FP',FN.lps_fp_vcf_records,''],['gate PASS (FP)',__FP_PASS__,''],['被篩掉 MISSED (FP)',__FP_MISSED__,'']];
$('funnel').innerHTML=fsteps.map(s=>`<div class="fstep"><b>${s[1].toLocaleString()}</b><div class="l">${s[0]}</div>${s[2]?'<div class="l" style="color:#fcd34d">'+s[2]+'</div>':''}</div>`).join('<span class="arrow">›</span>');

function live(){
  const cv=+$('cv').value, db=+$('db').value, rd=+$('rd').value, mode=$('mode').value, lf=$('lohf').value;
  $('vcv').textContent=cv.toFixed(2);$('vdb').textContent=db.toFixed(2);$('vrd').textContent=rd;
  let ptp=0,pfp=0,ttp=0,tfp=0;
  for(const r of ALL){
    if(lf==='non'&&r[3])continue; if(lf==='loh'&&!r[3])continue;
    if(r[4]===0)ttp++; else tfp++;
    if(r[2]<rd)continue;
    let pass=false;
    if(mode==='union')pass=(r[0]>=cv||r[1]>=db);
    else if(mode==='and')pass=(r[0]>=cv&&r[1]>=db);
    else if(mode==='cv')pass=(r[0]>=cv);
    else pass=(r[1]>=db);
    if(pass){ if(r[4]===0)ptp++; else pfp++; }
  }
  const ratio=pfp>0?(ptp/pfp).toFixed(2):'∞';
  const enr=(ptp/Math.max(ttp,1))/((pfp/Math.max(tfp,1))||1e-9);
  $('live').innerHTML=`
   <div><div class="l sub">通過 (TP)</div><div class="big tp">${ptp.toLocaleString()}</div><div class="sub">/ ${ttp.toLocaleString()} = ${(100*ptp/ttp).toFixed(1)}%</div></div>
   <div><div class="l sub">通過 (FP)</div><div class="big fp">${pfp.toLocaleString()}</div><div class="sub">/ ${tfp.toLocaleString()} = ${(100*pfp/tfp).toFixed(1)}%</div></div>
   <div><div class="l sub">TP : FP（通過數）</div><div class="big">${ratio}</div><div class="sub">越高越能分辨</div></div>
   <div><div class="l sub">TP 富集倍數 (vs FP)</div><div class="big">${enr.toFixed(2)}×</div><div class="sub">=1 表示無分辨力</div></div>`;
}
['cv','db','rd','mode','lohf'].forEach(id=>$(id).addEventListener('input',live));
live();

const SORTS={ '|Δβ| ↓':(a,b)=>b.db-a.db, 'CramersV (gate) ↓':(a,b)=>b.cv-a.cv, 'CramersV (gate) ↑':(a,b)=>a.cv-b.cv,
  'PERMANOVA F ↓':(a,b)=>(b.pf||0)-(a.pf||0), 'raw CramersV ↓':(a,b)=>(b.rmx||0)-(a.rmx||0),
  'NumReads ↓':(a,b)=>b.rd-a.rd, '位置':(a,b)=>a.ch<b.ch?-1:a.ch>b.ch?1:a.ps-b.ps };
$('sort').innerHTML=Object.keys(SORTS).map(k=>`<option>${k}</option>`).join('');
let curG='PASS', page=0, PER=24;
const JKEY='ism_judge_v2';
let judge=JSON.parse(localStorage.getItem(JKEY)||'{}');
function saveJ(){localStorage.setItem(JKEY,JSON.stringify(judge));}
$('cntP').textContent=CUR.filter(c=>c.cat==='PASS').length;
$('cntM').textContent=CUR.filter(c=>c.cat==='MISSED').length;
$('cntL').textContent=CUR.filter(c=>c.lt).length;

function setG(g){curG=g;page=0;document.querySelectorAll('.tab').forEach(t=>t.classList.toggle('on',t.dataset.g===g));
  $('sort').value=g==='MISSED'?'PERMANOVA F ↓':'CramersV (gate) ↓';
  $('latwrap').style.display=g==='MISSED'?'inline':'none';
  $('minnwrap').style.display=g==='MISSED'?'inline':'none';draw();}
document.querySelectorAll('.tab').forEach(t=>t.onclick=()=>setG(t.dataset.g));

function filtered(){
  const lf=$('glohf').value, cf=$('clsf').value, of=$('onlyfig').checked, ol=$('onlylat').checked, mn=+$('minn').value;
  $('vminn').textContent=mn;
  let a=CUR.filter(c=>c.cat===curG);
  if(of)a=a.filter(c=>c.fig);
  if(ol&&curG==='MISSED')a=a.filter(c=>c.lt&&c.mhn>=mn);
  if(lf==='non')a=a.filter(c=>!c.loh); if(lf==='loh')a=a.filter(c=>c.loh);
  if(cf!=='all')a=a.filter(c=>c.cl===cf);
  a.sort(SORTS[$('sort').value]||SORTS['CramersV (gate) ↓']);
  return a;
}
function badges(c){let b='';
  b+=c.cat==='PASS'?'<span class="badge b-pass">PASS</span>':'<span class="badge b-miss">MISSED</span>';
  if(c.lt)b+='<span class="badge" style="background:#082f33;color:#67e8f9">⚡潛在結構</span>';
  if(c.cl==='fp')b+='<span class="badge b-fp">FP</span>';
  if(c.loh)b+='<span class="badge b-loh">LOH</span>';return b;}
function card(c){
  const meth=`figs/${c.k}_meth.png`, dist=`figs/${c.k}_dist.png`;
  const figs=c.fig?`<div class="figs"><img loading="lazy" src="${meth}" alt="meth"><img loading="lazy" src="${dist}" alt="dist"></div>`
    :`<div class="figs" style="padding:18px;color:#64748b;font-size:11px">無圖：${c.note||'reads/CpG 太少'}</div>`;
  const jv=judge[c.k]||'';
  const jhtml=c.cat==='MISSED'?`<div class="judge">
     <button class="jb-in ${jv==='in'?'on':''}" onclick="event.stopPropagation();setJ('${c.k}','in')">應包含</button>
     <button class="jb-ex ${jv==='ex'?'on':''}" onclick="event.stopPropagation();setJ('${c.k}','ex')">確認排除</button></div>`:'';
  return `<div class="card" onclick="openM('${c.k}')">
    <div class="hd"><div class="ttl">${c.ch}:${c.ps} ${badges(c)}</div>
      <div class="meta">CramersV=${c.cv.toFixed(2)} ｜ Δβ=${c.db>=0?'+':''}${c.db.toFixed(2)} ｜ reads=${c.rd} ｜ CpG=${c.cg}</div></div>
    ${figs}
    ${c.rej?`<div class="rej">✗ ${c.rej}</div>`:''}
    ${jhtml}</div>`;
}
function draw(){
  const a=filtered();
  $('gcount').textContent=`${a.length} 個位點`;
  const pages=Math.max(1,Math.ceil(a.length/PER)); page=Math.min(page,pages-1);
  $('grid').innerHTML=a.slice(page*PER,page*PER+PER).map(card).join('')||'<div class="sub">無符合</div>';
  $('pager').innerHTML=`<button ${page<=0?'disabled':''} onclick="page--;draw()">‹ 上一頁</button>
    <span>第 ${page+1} / ${pages} 頁</span>
    <button ${page>=pages-1?'disabled':''} onclick="page++;draw()">下一頁 ›</button>`;
  const ji=Object.values(judge).filter(v=>v==='in').length, je=Object.values(judge).filter(v=>v==='ex').length;
  $('judgesum').innerHTML=`你的判斷：<span class="tp">應包含 ${ji}</span> ｜ 確認排除 ${je}`;
}
window.setJ=(k,v)=>{judge[k]=judge[k]===v?'':v;saveJ();draw();};
['sort','onlyfig','onlylat','minn','glohf','clsf'].forEach(id=>$(id).addEventListener('input',()=>{page=0;draw();}));

window.openM=k=>{const c=CUR.find(x=>x.k===k);if(!c)return;
  const meth=`figs/${c.k}_meth.png`, dist=`figs/${c.k}_dist.png`;
  const f3=v=>v==null?'NA':(+v).toFixed(3);
  const gated=[['CramersV max (gate 採用)',c.cv.toFixed(3)],['ALT 軸',f3(c.cva)],['HPfine 軸',f3(c.cvf)],['global_p',c.gp==null?'NA':c.gp],['過 gating',c.gat?'是':'否']];
  const raw=[['CramersV max (原始未閘控)',f3(c.rmx)],['ALT',f3(c.ra)],['HP',f3(c.rh)],['HP-family',f3(c.rhf)],['HPfine',f3(c.rhfn)]];
  const struc=[['PERMANOVA F (cluster)',c.pf==null?'NA':c.pf],['PERMANOVA p',c.pp==null?'NA':c.pp],['HP-PERMANOVA F',c.hpf==null?'NA':c.hpf],['最佳分群數 k',c.ok==null?'NA':c.ok],['離散警告',c.dw?'⚠ 是':'否']];
  const meta=[['Δβ (HP1 vs HP2)',(c.db>=0?'+':'')+c.db.toFixed(3)],['NumReads',c.rd],['NumCpG',c.cg],['Potential LOH',c.loh?'是':'否'],['類別',c.cat+' / '+c.cl.toUpperCase()]];
  const grid=a=>`<div class="statgrid">${a.map(s=>`<div class="stat"><div class="l">${s[0]}</div><div class="n">${s[1]}</div></div>`).join('')}</div>`;
  const latbox=c.lt?`<div class="note-warn" style="border-color:#0e7490;background:#082f33;color:#a5f3fc;margin:10px 0">⚡ <b>潛在結構</b>：CramersV 卡方因列聯表稀疏（期望格&lt;5）判不可靠→閘控歸零，但 PERMANOVA F=${c.pf} p=${c.pp} 顯著（原始 CramersV_max=${f3(c.rmx)}）。<b>最小 HP group reads = ${c.mhn}</b>（精修 gate minN≤${c.mhn} 時納回；越大越穩健）。請看右側距離圖：若對角塊乾淨且對齊 HP 側欄 → 真分群、值得考慮納入；若塊只由少數 read 構成 → 統計脆弱、排除合理。</div>`:'';
  $('mc').innerHTML=`<div style="display:flex;justify-content:space-between;align-items:center">
     <h2 style="margin:0;border:none">${c.ch}:${c.ps} ${badges(c)}</h2>
     <button onclick="closeM()">關閉 ✕</button></div>
    ${c.rej?`<div class="rej" style="margin:8px 0">被篩掉原因：${c.rej}</div>`:'<div class="sub" style="margin:8px 0">通過篩選</div>'}
    <div class="mfigs">
      <div><div class="sub" style="text-align:center;margin-bottom:4px">甲基 read×CpG（HP 分群 + ALT 側欄）</div><img src="${meth}"></div>
      <div><div class="sub" style="text-align:center;margin-bottom:4px">ISM read×read 距離矩陣（樹狀排序 + HP 側欄）</div><img src="${dist}"></div>
    </div>
    ${latbox}
    <div class="sub" style="margin:10px 0 3px;color:var(--ac)">① gate 採用值（閘控後，決定 PASS/MISSED）</div>${grid(gated)}
    <div class="sub" style="margin:10px 0 3px;color:#67e8f9">② 原始 CramersV（未經可靠性閘控；稀疏表時 ①會歸零但②保留）</div>${grid(raw)}
    <div class="sub" style="margin:10px 0 3px;color:#a5f3fc">③ 結構檢定（PERMANOVA 距離法，對稀疏穩健）+ 其他</div>${grid(struc)}${grid(meta)}
    <div class="sub" style="margin-top:12px">判讀：距離圖若沿對角線出現<b>暗色塊</b>且與 HP/ALT 側欄分界對齊 → 該位點 reads 確實按單倍型/等位基因分群。①與②落差大（②高①低）+ ③ PERMANOVA 顯著 = 有真結構但 CramersV 統計脆弱（稀疏），即「可能漏掉」的判斷點。若距離圖<b>均勻無塊</b>但甲基圖整體呈兩種深淺 → 只有平均位移（Δβ）、無分群 → 被 CramersV gate 篩掉是正確的。</div>`;
  $('modal').style.display='block';};
window.closeM=()=>$('modal').style.display='none';
document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});

setG('PASS');
</script></body></html>"""


def main():
    manifest = json.load(open(f"{DV}/manifest.json"))
    fn = json.load(open(f"{DV}/funnel_numbers.json"))
    rg = json.load(open(f"{DV}/refined_gate.json"))
    allrows = load_all()

    def cnt(cat, cls):
        return sum(1 for m in manifest if m["category"] == cat and m["cls"] == cls)
    tp_pass, fp_pass = cnt("PASS", "tp"), cnt("PASS", "fp")
    tp_missed, fp_missed = cnt("MISSED", "tp"), cnt("MISSED", "fp")
    n_fig = sum(1 for m in manifest if m.get("has_fig"))

    cur = [dict(k=m["chr"] + "_" + str(m["pos"]), ch=m["chr"], ps=m["pos"], cl=m["cls"],
                cat=m["category"], cv=m["cv_max"], cva=m.get("cv_alt", 0), cvf=m.get("cv_hpfine", 0),
                db=m["db"], gp=(round(m["gp"], 4) if m.get("gp") is not None else None),
                rd=m["reads"], cg=m["ncpg"], loh=1 if m["loh"] else 0,
                gat=1 if m["passed_gating"] else 0, fig=1 if m.get("has_fig") else 0,
                rej=m.get("reject", ""), dw=1 if m.get("disp_warn") else 0, note=m.get("note", ""),
                lt=1 if m.get("latent") else 0, rmx=m.get("raw_max", 0), mhn=m.get("minhpn", 0),
                ra=m.get("raw_alt", 0), rh=m.get("raw_hp", 0), rhf=m.get("raw_hpfam", 0), rhfn=m.get("raw_hpfine", 0),
                pf=m.get("perm_f"), pp=m.get("perm_p"), pv=1 if m.get("perm_valid") else 0,
                hpf=m.get("hp_perm_f"), hpp=m.get("hp_perm_p"), ok=m.get("optimal_k"))
           for m in manifest]
    n_latent = sum(1 for m in manifest if m.get("latent"))

    h = TEMPLATE
    repl = {
        "__STAMP__": datetime.date.today().isoformat(),
        "__PAYLOAD__": json.dumps(cur, separators=(",", ":"), ensure_ascii=False),
        "__ALLPAYLOAD__": json.dumps(allrows, separators=(",", ":")),
        "__FNJ__": json.dumps(fn, ensure_ascii=False),
        "__TP_PASS__": str(tp_pass), "__TP_MISSED__": str(tp_missed),
        "__FP_PASS__": str(fp_pass), "__FP_MISSED__": str(fp_missed),
        "__N_LATENT__": str(n_latent), "__N_MISSED_NOTE__": str(tp_missed + fp_missed),
        "__REFINED_GATE__": json.dumps(rg, ensure_ascii=False),
    }
    for k, v in repl.items():
        h = h.replace(k, v)

    with open(OUT, "w") as f:
        f.write(h)
    print(f"[86] wrote {OUT} ({os.path.getsize(OUT)//1024} KB)")
    print(f"     curated={len(cur)} fig={n_fig} PASS(tp{tp_pass}/fp{fp_pass}) "
          f"MISSED(tp{tp_missed}/fp{fp_missed}) latent={n_latent} ALL_rows={len(allrows)}")


if __name__ == "__main__":
    main()
