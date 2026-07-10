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
from datetime import datetime

RV = os.environ.get("SM_RV", "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/layered_region_view_HCC1395.json")
OUT = os.environ.get("SM_OUT", "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260706_layered_reconstruction_workstation.standalone.html")
D = json.load(open(RV, encoding="utf-8"))
C = D["census"]
SAMPLE = D.get("sample", "?")
GENERATED_AT = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
BACKBONE = C.get("U1_backbone_source", "ClairS paired PASS somatic backbone")
SAMPLE_NAMES = ["HCC1395", "COLO829", "H1437", "H2009", "HCC1395_DORADO", "HCC1937", "HCC1954"]
SAMPLE_OPTIONS = "".join(
    '<option value="{}.html"{}>{}</option>'.format(name, " selected" if name == SAMPLE else "", name)
    for name in SAMPLE_NAMES
)

# 精簡 regions 給前端(trees 已 ≤6;去掉冗長欄位以縮小)
regions = D["regions"]
DATA_JSON = json.dumps({"sample": SAMPLE, "census": C, "regions": regions}, ensure_ascii=False, separators=(",", ":"))

HTML = r"""<!doctype html><html lang="zh-Hant"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HP-family 分層樹重建 · __SAMPLE__</title>
<style>
:root{color-scheme:light;--bg:#f4f6fa;--paper:#fff;--fg:#172033;--mut:#566176;--line:#d8dee9;--card:#f8fafc;--accent:#145fc6;
 --det:#137a37;--amb:#985000;--cap:#7136b8;--rec:#b42318;--hp1:#145fc6;--hp2:#c51f68;--hp3:#087d8f;--none:#5d6b82;
 --info-bg:#edf5ff;--warn-bg:#fff7e8;--danger-bg:#fff1ef;--shadow:0 10px 30px rgba(23,32,51,.08)}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.55 ui-sans-serif,system-ui,"Segoe UI","Microsoft JhengHei",sans-serif}
.skip{position:absolute;left:-9999px;top:8px}.skip:focus{left:12px;z-index:10;background:#fff;padding:8px 12px;border:2px solid var(--accent);border-radius:6px}
.wrap{max-width:1240px;margin:0 auto;padding:22px}
.topbar{background:var(--paper);border:1px solid var(--line);border-top:4px solid var(--accent);border-radius:12px;padding:18px 20px;box-shadow:var(--shadow)}
.crumbs,.top-actions{display:flex;align-items:center;gap:9px;flex-wrap:wrap}.crumbs{font-size:12.5px;color:var(--mut);margin-bottom:12px}.crumbs a{color:var(--accent);font-weight:700;text-decoration:none}.crumbs a:hover{text-decoration:underline}
.top-actions{margin-left:auto}.sample-switch{background:var(--paper);border:1px solid var(--line);border-radius:7px;padding:6px 9px;color:var(--fg)}
.heading-row{display:flex;gap:18px;align-items:flex-start}.heading-row>div:first-child{flex:1;min-width:0}
h1{font-size:23px;line-height:1.25;margin:0 0 5px;letter-spacing:-.01em}h2{font-size:16px;margin:22px 0 10px;border-left:4px solid var(--accent);padding-left:9px}
.sub{color:var(--mut);font-size:13px;margin:0}
.statusline{display:flex;gap:7px;flex-wrap:wrap;margin-top:14px}.status{display:inline-flex;align-items:center;gap:5px;border:1px solid var(--line);border-radius:999px;padding:3px 9px;font-size:11.5px;color:var(--mut);background:var(--card)}
.status.current{border-color:#9bd0ad;background:#eefbf2;color:#0f6330;font-weight:800}.status.pending{border-color:#edc47c;background:var(--warn-bg);color:#744000}
.explainer{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin-top:14px}.explain-card{border:1px solid var(--line);border-radius:9px;padding:10px 12px;background:var(--card);font-size:12.5px}.explain-card b{display:block;font-size:11px;letter-spacing:.07em;text-transform:uppercase;color:var(--accent);margin-bottom:3px}
.summary{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:9px;margin:14px 0}.metric{background:var(--paper);border:1px solid var(--line);border-radius:9px;padding:10px 12px;box-shadow:0 3px 12px rgba(23,32,51,.04)}.metric .v{font:700 20px/1.1 ui-monospace,"SF Mono",Menlo,monospace}.metric .k{font-size:11.5px;color:var(--mut);margin-top:4px}.metric.warn{border-top:3px solid var(--amb)}.metric.good{border-top:3px solid var(--det)}.metric.pending{border-top:3px solid var(--cap)}
.overview{background:var(--paper);border:1px solid var(--line);border-radius:10px;margin:12px 0}.overview>summary{cursor:pointer;padding:11px 14px;font-weight:750;color:var(--accent)}.overview[open]>summary{border-bottom:1px solid var(--line)}.overview-body{padding:0 14px 14px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:6px 0;background:var(--paper)}
th,td{border:1px solid var(--line);padding:6px 8px;text-align:left}th{background:var(--card);font-weight:700}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}
.tag{display:inline-block;padding:1px 7px;border-radius:10px;font-size:11px;font-weight:600;color:#fff}
.t-det{background:var(--det)}.t-amb{background:var(--amb)}.t-cap{background:var(--cap)}.t-rec{background:var(--rec)}.t-none{background:var(--none)}
.hp1{color:var(--hp1)}.hp2{color:var(--hp2)}.hp3{color:var(--hp3)}.hpnone{color:var(--none)}
.ctrl{display:flex;flex-wrap:wrap;gap:9px;align-items:flex-end;margin:10px 0 12px;font-size:13px}.control{display:grid;gap:3px}.control>span{font-size:10.5px;font-weight:750;letter-spacing:.04em;color:var(--mut)}
.ctrl select,.ctrl input{background:var(--paper);color:var(--fg);border:1px solid var(--line);border-radius:7px;padding:7px 9px;min-height:36px}.ctrl select:focus,.ctrl input:focus,button:focus-visible,a:focus-visible{outline:3px solid rgba(20,95,198,.28);outline-offset:2px}
.rcard{border:1px solid var(--line);border-radius:9px;margin:9px 0;background:var(--card);overflow:hidden}
.rhead{padding:8px 12px;cursor:pointer;display:flex;flex-wrap:wrap;gap:10px;align-items:center}
.rhead:hover{background:rgba(128,128,128,.06)}
.rbody{padding:0 12px 10px;display:none}.rcard.open .rbody{display:block}
.lin{border-top:1px dashed var(--line);padding:9px 0 9px 10px;border-left:3px solid var(--line);margin-top:6px}
.lin.f1{border-left-color:var(--hp1)}.lin.f2{border-left-color:var(--hp2)}.lin.f3{border-left-color:var(--hp3)}.lin.fnone{border-left-color:var(--none)}
.mono{font-family:ui-monospace,"SF Mono",Menlo,monospace}
.trace{font-size:12px;color:var(--mut);margin:7px 0 0;padding:7px 9px;border-radius:6px;background:#f2f5f9}.trace div{padding:1px 0}
/* 樹切換器:一次一棵大圖 + prev/next + 計數 + thumbnail 跳轉 */
.tsw{border:1px solid var(--line);border-radius:8px;padding:7px 9px;margin:6px 0;background:var(--bg)}
.tswhead{font-size:12px;color:var(--fg);margin-bottom:5px;display:flex;align-items:center;gap:6px;flex-wrap:wrap}
.tstage{overflow-x:auto;padding:4px 0;min-height:60px}
.tslide{text-align:center}
.tcap{font-size:11px;color:var(--mut);text-align:center;margin-top:3px}
.tnav{padding:1px 9px;font-size:13px;line-height:1.4}.tctr{font-size:12px;color:var(--mut);min-width:46px;text-align:center;font-variant-numeric:tabular-nums}
.thumbs{display:flex;gap:4px;flex-wrap:wrap;margin-top:6px}
.thumb{cursor:pointer;font-size:11.5px;padding:4px 9px;border:1px solid var(--line);border-radius:9px;color:var(--mut);background:var(--paper);user-select:none;min-height:30px}
.thumb.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.thumb:hover{border-color:var(--accent)}
.carousel{display:flex;gap:8px;overflow-x:auto;padding:6px 2px}
.tree{flex:0 0 auto;border:1px solid var(--line);border-radius:7px;padding:5px;background:var(--bg)}
.note{font-size:12.5px;color:var(--mut);background:var(--card);border:1px solid var(--line);border-radius:7px;padding:8px 10px;margin:8px 0}
.pill{font-size:11.5px;padding:2px 7px;border:1px solid var(--line);border-radius:8px;color:var(--mut);background:var(--paper)}
.bar{display:inline-block;height:9px;border-radius:3px;vertical-align:middle}
button.tg,.action{background:var(--paper);color:var(--fg);border:1px solid var(--line);border-radius:7px;padding:6px 10px;cursor:pointer;font-size:12px;font-weight:650}.action.primary{background:var(--accent);border-color:var(--accent);color:#fff}.action:hover,button.tg:hover{border-color:var(--accent)}
#more{margin:12px auto;display:block}
/* topology 式兩欄:左 region 列表 + 右詳情面板 */
.main{display:grid;grid-template-columns:minmax(300px,380px) minmax(0,1fr);gap:12px;align-items:start}.main>*{min-width:0}
.list{max-height:72vh;overflow-y:auto;border:1px solid var(--line);border-radius:9px;background:var(--paper)}
.row{display:block;width:100%;padding:8px 10px;border:0;border-bottom:1px solid var(--line);cursor:pointer;font:inherit;font-size:12px;text-align:left;color:var(--fg);background:var(--paper)}
.row:hover{background:rgba(128,128,128,.06)}.row.sel{background:rgba(37,99,235,.10);border-left:3px solid var(--accent)}
.detail{border:1px solid var(--line);border-radius:9px;background:var(--paper);padding:14px;min-height:220px;overflow:hidden}
.detail h3{font-size:15px;margin:0 0 6px;font-family:ui-monospace,monospace}
.kv{display:flex;flex-wrap:wrap;gap:6px;margin:6px 0}.kv .b{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:2px 8px;font-size:11.5px}
.text-link{border:0;background:none;color:var(--accent);text-decoration:underline;font:inherit;cursor:pointer;padding:0}.evidence-state{font-size:10.5px;font-weight:800;text-transform:uppercase;letter-spacing:.04em;border-radius:4px;padding:2px 6px}.observed{background:#e8f2ff;color:#0e4d9c}.inferred{background:#f2eafe;color:#5b248e}.confounded{background:#fff0e1;color:#7d3e00}.unavailable{background:#eef1f5;color:#48556b}
svg{max-width:100%;height:auto}.footer{font-size:11.5px;color:var(--mut);margin:14px 0 4px;padding-top:10px;border-top:1px solid var(--line)}
@media(max-width:900px){.wrap{padding:14px}.heading-row{display:block}.top-actions{margin:12px 0 0}.explainer{grid-template-columns:1fr}.summary{grid-template-columns:repeat(2,minmax(0,1fr))}.grid{grid-template-columns:1fr}.main{grid-template-columns:1fr}.list{max-height:44vh}.detail{scroll-margin-top:12px}#dash table{display:block;overflow-x:auto;white-space:nowrap}#dash td:last-child{white-space:normal;min-width:240px}}
@media(max-width:560px){body{font-size:13.5px}.wrap{padding:10px}.topbar{padding:14px}.summary{grid-template-columns:1fr 1fr}.metric .v{font-size:17px}.ctrl{display:grid;grid-template-columns:1fr 1fr}.control.search{grid-column:1/-1}.ctrl select,.ctrl input{width:100%}.list{max-height:38vh}.tag,.pill{font-size:10.5px}.detail{padding:10px}.kv{gap:4px}}
</style></head><body>
<a class="skip" href="#regions">跳至 Region 瀏覽器</a>
<main class="wrap">
<header class="topbar">
 <div class="crumbs"><a href="index.html">全部樣本總覽</a><span aria-hidden="true">/</span><span class="mono">__SAMPLE__</span></div>
 <div class="heading-row"><div><h1>HP-family 分層樹重建 · <span class="mono">__SAMPLE__</span></h1>
 <p class="sub">家族優先於算法：L0 HP 家族 → L1 枚舉全最小樹 → L2 CN confound → L3 甲基（目前 pending、未參與樹排序）</p></div>
 <div class="top-actions"><label for="sample-switch" class="control"><span>切換樣本</span><select id="sample-switch" class="sample-switch" onchange="location.href=this.value">__SAMPLE_OPTIONS__</select></label><button class="action" type="button" id="copy-link">複製目前 Region 連結</button></div></div>
 <div class="statusline"><span class="status current">CURRENT / CANONICAL</span><span class="status">Backbone: __BACKBONE__</span><span class="status">Generated: __GENERATED_AT__</span><span class="status pending">L3 methylation: PENDING / 不 rank 樹</span></div>
 <div class="explainer"><div class="explain-card"><b>回答什麼</b>同一 germline-HP 家族中，哪些最小樹與目前 read 證據相容，以及哪些連接仍不唯一。</div><div class="explain-card"><b>證據怎麼看</b><span class="evidence-state observed">Observed</span> 是 read 直接觀測；<span class="evidence-state inferred">Inferred</span> 是 solver 補出的隱藏祖先或組合。</div><div class="explain-card"><b>不能宣稱什麼</b>單一 bulk、CN confound 或 partial read 下的候選，不等於已確認 subclone；ambiguous/capped 本身就是結果。</div></div>
</header>
<div class="summary" id="summary" aria-label="樣本摘要"></div>
<details class="overview"><summary>展開完整樣本統計、分母字典與驗證狀態</summary><div class="overview-body" id="dash"></div></details>
<section id="regions" aria-labelledby="regions-title"><h2 id="regions-title">Region 瀏覽器：篩選 → 選取 → 檢視各 HP 家族證據</h2>
<div class="ctrl" aria-label="Region 篩選與排序">
 <label class="control" for="fdet"><span>確定性</span><select id="fdet"><option value="">全部</option></select></label>
 <label class="control" for="fhp"><span>HP multiplicity</span><select id="fhp"><option value="">全部</option><option value="2">多-HP (2)</option><option value="1">single-HP (1)</option><option value="0">無 germline (0)</option></select></label>
 <label class="control" for="fchr"><span>染色體</span><select id="fchr"><option value="">全部</option></select></label>
 <label class="control" for="fsort"><span>排序</span><select id="fsort"><option value="coord">基因組座標</option><option value="nsnv">複雜度 (sSNV)</option><option value="ntree">枚舉樹總數</option><option value="hp">HP multiplicity</option></select></label>
 <label class="control search" for="fq"><span>搜尋 Region</span><input id="fq" placeholder="例如 chr17:430..." size="18"></label>
 <span class="pill" id="cnt" role="status" aria-live="polite"></span>
</div>
<div class="main"><div class="list" id="list" role="list" aria-label="Region 清單"></div><div class="detail" id="detail" tabindex="-1" aria-live="polite"><div class="note"><b>尚未選取 Region。</b> 從左側清單選擇一區，這裡會依序呈現判定、位點證據、家族樹、read 組合與 L0–L3 軌跡。</div></div></div>
</section>
<footer class="footer">資料來源：<span class="mono">layered_region_view___SAMPLE__.json</span> · 樹由 solver 枚舉並經 V1–V7 驗證 · L3 甲基目前僅標示 availability，未進入樹排序。</footer>
</main>
<script type="application/json" id="workstation-data">__DATA__</script>
<script>
const D=JSON.parse(document.getElementById('workstation-data').textContent);const C=D.census;const R=D.regions;
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
function renderSummary(){
 const n=C.n_regions||0,rd=C.region_determinacy||{};
 document.getElementById('summary').innerHTML=
  '<div class="metric"><div class="v">'+(C.U1_sSNV_somatic_total||0).toLocaleString()+'</div><div class="k">somatic sSNV · paired PASS backbone</div></div>'+
  '<div class="metric"><div class="v">'+n.toLocaleString()+'</div><div class="k">regions · 本頁主分母</div></div>'+
  '<div class="metric good"><div class="v">'+pct(rd.all_determined||0,n)+'</div><div class="k">all-determined · 全 HP 家族唯一樹</div></div>'+
  '<div class="metric warn"><div class="v">'+pct((rd.has_ambiguous||0)+(rd.has_capped||0),n)+'</div><div class="k">ambiguous + capped · 不確定性是結果</div></div>'+
  '<div class="metric pending"><div class="v">PENDING</div><div class="k">L3 methylation · 尚未進入樹排序</div></div>';
}
// ================= dashboard =================
function dash(){
 const L0=C.L0,L1=C.L1,L2=C.L2,rd=C.region_determinacy,mult=C.hp_multiplicity;
 const nreg=C.n_regions, nlin=L1.n_lineage_units;
 const barw=v=>Math.max(2,Math.round(v/nreg*160));
 let h='';
 // 層級計數
 h+='<h2>1. 層級計數（U1→U6）</h2><table><tr><th>層</th><th>單位</th><th class="n">'+esc(D.sample)+' 計數</th><th>說明</th></tr>'
 +row('U1','somatic sSNV',(C.U1_sSNV_somatic_total||0).toLocaleString(),esc(C.U1_backbone_source||'ClairS paired PASS somatic backbone'))
 +row('U3','region',nreg.toLocaleString(),'multilocus 分析群(≤8sSNV窗,主分母);linkage全span='+((C.U3_linkage_regions_full_span||0).toLocaleString())+' 同批區(位置99.9%重疊)')
 +row('U4','lineage-unit',nlin.toLocaleString(),'region×HP家族(1/2/3);樹的自然單位')
 +row('U5','tree',(C.U5_trees?C.U5_trees.sum_ntrees_noncapped.toLocaleString():'—'),'Σ枚舉樹(non-capped);avg '+(C.U5_trees?(C.U5_trees.sum_ntrees_noncapped/nlin).toFixed(2):'—')+' 樹/unit')
 +row('U6','隱藏祖先',(C.U6_hidden?C.U6_hidden.sum_hidden.toLocaleString():'—'),'Σ推斷未觀測 clone(H_*)')
 +'</table>';
 // HP-multiplicity + region determinacy
 h+='<div class="grid"><div><h2>2. HP-multiplicity（每區幾個 germline 樹）</h2><table><tr><th>germline lineage 樹數</th><th class="n">區數</th><th class="n">比例</th><th></th></tr>'
 +mrow('2（多-HP：雙親代染色體各一樹）',mult['2'],nreg,'var(--hp2)')
 +mrow('1（single-HP）',mult['1'],nreg,'var(--hp1)')
 +mrow('0（只 somatic3/none）',mult['0'],nreg,'var(--none)')
 +'</table><div class="note"><span class="evidence-state confounded">Model change</span> 多-HP '+pct(mult['2'],nreg)+'：兩個 germline 家族須分開建樹，避免舊 pooled 模型混淆 allelic 與 clonal 關係。</div></div>';
 h+='<div><h2>3. Region-level determinacy</h2><table><tr><th>判定</th><th class="n">區數</th><th class="n">比例</th><th></th></tr>'
 +drow('all-determined（全家族唯一樹）',rd.all_determined,nreg,'t-det')
 +drow('has-ambiguous（≥1家族多樹）',rd.has_ambiguous,nreg,'t-amb')
 +drow('has-capped（太密）',rd.has_capped,nreg,'t-cap')
 +drow('has-recurrence',rd.has_recurrence,nreg,'t-rec')
 +drow('no-germline（只3/none）',rd.no_germline_lineage,nreg,'t-none')
 +'</table><div class="note">定義：region 確定 ⟺ <b>所有</b> germline(1/2) lineage 都 determined（多-HP 需雙確定）→ 嚴於 lineage-unit。</div></div></div>';
 // 比例字典(lineage 次分母) + L2
 h+='<h2>4. 比例字典（次分母 lineage-unit '+nlin.toLocaleString()+'） + L2 CN</h2><table><tr><th>比例</th><th>分子÷分母</th><th class="n">值</th></tr>'
 +prow('determined(lineage)',L1.determinacy_lineage.determined,nlin)
 +prow('ambiguous(lineage)',L1.determinacy_lineage['ambiguous_structure(多完成/多結構)'],nlin)
 +prow('capped(lineage)',L1.determinacy_lineage['capped(太密;枚舉未完)'],nlin)
 +prow('recurrence(lineage)',L1.determinacy_lineage.recurrence_required||0,nlin)
 +'<tr><td>L2 CN artifact（recurrence 內）</td><td class="mono">'+(L2.cn_split['artifact(m>1;CN-amp)']||0)+' ÷ '+L2.n_recurrence_sent_to_cn+' recurrence</td><td class="n">'+pct(L2.cn_split['artifact(m>1;CN-amp)']||0,L2.n_recurrence_sent_to_cn)+'</td></tr>'
 +'</table>';
 h+='<div class="note">V1–V7 驗證：'+(L1.all_V1V7_pass?'<b style="color:var(--det)">ALL PASS（0 fail）</b>':('<b style="color:var(--rec)">'+L1.n_verify_fail+' fail</b>'))+' · 分母鐵則：region determined '+pct(rd.all_determined,nreg)+' ≠ lineage determined '+pct(L1.determinacy_lineage.determined,nlin)+'（單位不同不可比）</div>';
 // ⑤ CN confound 分層 + ⑥ mutation-bearing 分母 + clone 數 c(從 R client-side 算)
 (function(){
  let cn={neutral:0,gain:0,amp:0,loss:0,loh:0,unknown:0,other:0};
  let mb=0,root=0,noPos=0,detMB=0,detRoot=0,cdist={c0:0,c1:0,c2:0};
  R.forEach(r=>{cn[cn[r.cn]!==undefined?r.cn:'other']++;
   const pos=r.positions||[];let anyAlt=false;
   r.lineages.forEach(L=>{const cov=L.obs_col_coverage||{};pos.forEach(p=>{const c=cov[String(p)];if(c&&c[1]>0)anyAlt=true;});});
   if(!pos.length)noPos++;else if(anyAlt)mb++;else root++;
   if(r.region_determinacy==='all_determined'&&pos.length){if(anyAlt)detMB++;else detRoot++;}
   let g=new Set();r.lineages.forEach(L=>{const p=L.obs_populations||{};for(const kk in p)if(kk.indexOf('A')>=0)g.add(kk);});
   cdist[g.size===0?'c0':g.size===1?'c1':'c2']++;});
  const clean=cn.neutral,altered=cn.gain+cn.amp+cn.loss+cn.loh,cngain=cn.gain+cn.amp,cnloss=cn.loss+cn.loh;
  const wpos=mb+root;
  h+='<div class="grid"><div><h2>5. CN confound 分層</h2><table><tr><th>CN 類</th><th class="n">區數</th><th class="n">比例</th></tr>'
   +'<tr style="background:#f0fdf4"><td>neutral（CN-clean·可算 CCF）</td><td class="n">'+clean.toLocaleString()+'</td><td class="n">'+pct(clean,nreg)+'</td></tr>'
   +'<tr style="background:#fef2f2"><td>gain/amp（<b>read≠CCF·排序停用</b>）</td><td class="n">'+cngain.toLocaleString()+'</td><td class="n">'+pct(cngain,nreg)+'</td></tr>'
   +'<tr style="background:#fffbeb"><td>loss/loh（LOH-unmask confound）</td><td class="n">'+cnloss.toLocaleString()+'</td><td class="n">'+pct(cnloss,nreg)+'</td></tr>'
   +'<tr><td>unknown</td><td class="n">'+cn.unknown.toLocaleString()+'</td><td class="n">'+pct(cn.unknown,nreg)+'</td></tr>'
   +'</table><div class="note"><span class="evidence-state confounded">CN check</span> 本樣本 neutral '+pct(clean,nreg)+'、CN-altered '+pct(altered,nreg)+'、unknown '+pct(cn.unknown,nreg)+'。只有 neutral 可直接把 read/VAF 當作未校正 CCF 參考；其餘維持 confounded。</div></div>';
  h+='<div><h2>6. Mutation-bearing 分母 + clone 數 c</h2><table><tr><th>類</th><th class="n">區數</th><th class="n">比例</th></tr>'
   +'<tr style="background:#f0fdf4"><td>mutation-bearing（≥1 位點實測 ALT）</td><td class="n">'+mb.toLocaleString()+'</td><td class="n">'+pct(mb,wpos)+'</td></tr>'
   +'<tr style="background:#fef2f2"><td>root-only（全位點 0 ALT·無突變可重建）</td><td class="n">'+root.toLocaleString()+'</td><td class="n">'+pct(root,wpos)+'</td></tr>'
   +'<tr><td>c=0（無全跨 ALT）/ c=1 / <b>c≥2（subclonal 候選·未確認）</b></td><td class="n">'+cdist.c0.toLocaleString()+' / '+cdist.c1.toLocaleString()+' / '+cdist.c2.toLocaleString()+'</td><td class="n">'+pct(cdist.c2,nreg)+' c≥2</td></tr>'
   +'</table><div class="note"><span class="evidence-state observed">Observed denominator</span> root-only = region 有 mlhp positions、但所有位點皆為 0 ALT 實測；此樣本 '+detRoot.toLocaleString()+' 個 all-determined 屬此類。扣除後 mutation-bearing determined 為 '+detMB.toLocaleString()+'（'+pct(detMB,rd.all_determined)+'）。c≥2 只代表 co-phased ALT genotype 候選，<b>不等於已確認 subclone</b>；仍須排除 CN 與 single-bulk 限制。分母（有 positions）='+wpos.toLocaleString()+'；另有 '+noPos.toLocaleString()+' 區無 mlhp positions。</div></div></div>';
 })();
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
const chromRank=c=>{const x=String(c||'').replace(/^chr/i,'').toUpperCase();if(/^\d+$/.test(x))return +x;if(x==='X')return 23;if(x==='Y')return 24;if(x==='M'||x==='MT')return 25;return 99;};
const SORT={coord:(a,b)=>(chromRank(a.chrom)-chromRank(b.chrom))||(a.start-b.start)||(a.end-b.end),
 nsnv:(a,b)=>b.n_sSNV-a.n_sSNV, ntree:(a,b)=>_sumT(b)-_sumT(a), hp:(a,b)=>b.hp_multiplicity-a.hp_multiplicity};
function applyFilter(){
 const fd=document.getElementById('fdet').value,fh=document.getElementById('fhp').value,
  fc=document.getElementById('fchr').value,fq=document.getElementById('fq').value.trim().toLowerCase(),
  so=document.getElementById('fsort').value;
 FILT=R.map((r,i)=>[r,i]).filter(([r])=>(!fd||r.region_determinacy===fd)&&(!fh||String(r.hp_multiplicity)===fh)
  &&(!fc||r.chrom===fc)&&(!fq||r.region.toLowerCase().indexOf(fq)>=0));
 FILT.sort((A,B)=>(SORT[so]||SORT.coord)(A[0],B[0]));
 const host=document.getElementById('list');
 host.innerHTML=FILT.slice(0,400).map(([r,i])=>{const[rt,rc]=regTag(r.region_determinacy);
  return '<button type="button" class="row" data-i="'+i+'" aria-selected="false"><span class="mono" style="font-weight:700">'+esc(r.region)+'</span> <span class="tag '+rc+'">'+rt+'</span><br>'
   +'<span class="pill">HP×'+r.hp_multiplicity+(r.is_multiHP?'·多':'')+'</span><span class="pill">'+r.n_sSNV+' sSNV</span><span class="pill">'+r.lineages.length+' lineage</span><span class="pill">'+_sumT(r)+' 樹</span><span class="pill">cn '+esc(r.cn)+'</span></button>';
 }).join('')+(FILT.length>400?'<div class="note" style="margin:6px">顯示前 400 / 共 '+FILT.length.toLocaleString()+' 區；請用條件縮小。</div>':'');
 document.getElementById('cnt').textContent=FILT.length.toLocaleString()+' 區';
 host.querySelectorAll('.row').forEach(x=>x.onclick=()=>show(+x.dataset.i,x,true));
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
  return '<div style="display:inline-block;margin:3px 10px 3px 0;vertical-align:top"><div style="font-size:11px;font-weight:650">S'+(i+1)+' × S'+(j+1)+(cov?'':' <span style="color:var(--rec)">· 0 全跨 read</span>')+(gap?' <span style="color:var(--rec)">· READ-SPAN GAP</span>':'')+'</div><table style="font-size:11px;border-collapse:collapse;margin-top:1px"><tr><td style="padding:1px 6px;color:var(--mut)"></td><td style="padding:1px 6px;color:var(--mut)">S'+(j+1)+'·R</td><td style="padding:1px 6px;color:var(--mut)">S'+(j+1)+'·A</td></tr><tr><td style="color:var(--mut)">S'+(i+1)+'·R</td>'+cell(m.RR)+cell(m.RA)+'</tr><tr><td style="color:var(--mut)">S'+(i+1)+'·A</td>'+cell(m.AR)+cell(m.AA,true)+'</tr></table></div>';
 }).join('');
 return '<details style="margin:6px 0"><summary style="cursor:pointer;font-size:12px;font-weight:700">位點 pairwise 共現（僅全跨 read；A×A=co-phased 直接證據）'+(anyGap?' <span style="color:var(--rec)">· READ-SPAN GAP</span>':'')+'</summary><div style="margin-top:4px">'+(gts.length?blocks:'<span class="note" style="background:none;border:none;padding:0;color:var(--rec)">此區 0 全跨 read：任何位點對的組合皆無法直接實測，只能由 partial read 推斷。</span>')+'</div><div class="note" style="background:none;border:none;padding:2px 0"><span class="evidence-state inferred">Inferred</span> 兩位點各自有 ALT、但 A×A=0，代表沒有 read 同時連結兩個 ALT；組合層仍未被直接觀測。</div></details>';
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
  return '<div style="flex:1;min-width:180px"><b>'+side+'</b> <button type="button" class="text-link mono" onclick="show('+o.gi+',null,true)">'+esc(x.region)+'</button> <span class="note" style="background:none;border:none;padding:0">'+x.n_sSNV+'sSNV·gap '+gk+'</span><br>'+(sm.length?'<span style="color:var(--det);font-size:11px">共享 '+sm.join(' / ')+'</span>':'<span class="note" style="background:none;border:none;padding:0;font-size:11px">無共享屬性</span>')+'</div>';}
 return '<details style="margin:6px 0"><summary style="cursor:pointer;font-size:12px;font-weight:700">前後相鄰區（read-span 不跨區；僅描述性鄰接）</summary><div style="display:flex;gap:14px;margin-top:5px;flex-wrap:wrap">'+cell(prev,prev?r.start-prev.x.end:0,'← 前')+cell(next,next?next.x.start-r.end:0,'後 →')+'</div></details>';
}
// ===== 每 germline-HP 家族觀測 read(全跨 populations + partial subreads)=====
function readObsBlock(L){
 const fp=L.obs_populations||{},sp=L.obs_subreads||{};const nf=Object.keys(fp).length,ns=Object.keys(sp).length;
 if(!nf&&!ns)return '';
 const row=(g,typ,c,col)=>'<tr><td class="mono" style="padding:2px 6px">'+g+'</td><td style="padding:2px 6px;color:'+col+'">'+typ+'</td><td style="padding:2px 6px;text-align:right">'+c+'</td></tr>';
 let rows=Object.entries(fp).sort((a,b)=>b[1]-a[1]).map(([g,c])=>row(g,'全跨(co-phase)',c,'#2563eb')).join('')
  +Object.entries(sp).sort((a,b)=>b[1]-a[1]).map(([g,c])=>row(g,'partial(X=未覆蓋)',c,'#a855f7')).join('');
 return '<details style="margin:5px 0"><summary style="cursor:pointer;font-size:12px">此家族觀測 read（'+nf+' 全跨群 + '+ns+' partial 群·原始 mlhp）</summary>'
  +'<table style="font-size:10.5px;margin-top:3px;border-collapse:collapse"><tr style="border-bottom:1px solid #ddd"><th style="padding:2px 6px;text-align:left">genotype</th><th style="padding:2px 6px;text-align:left">類型</th><th style="padding:2px 6px">reads</th></tr>'+rows+'</table>'
  +'<div class="note" style="background:none;border:none;padding:2px 0;font-size:9.5px">全跨群=同一 read 覆蓋全部位點(可定 phasing);partial=只覆蓋部分位點(單獨看到某位點的 ALT·無法連鎖)。'+(nf===0?'<b style="color:#a855f7">此家族 0 全跨群→組合全靠 partial 拼(推斷)。</b>':'')+'</div></details>';
}
// ===== CN confound 徽章(區級) =====
function cnBadge(r){
 const c=r.cn;
 if(c==='neutral')return '<span class="b" style="background:#dcfce7;color:#166534">cn neutral·CN-clean</span>';
 if(c==='gain'||c==='amp')return '<span class="b" style="background:#fee2e2;color:#8a1c16" title="CN-altered：read count 不可直接視為 CCF">cn '+esc(c)+' · CN-CONFOUNDED</span>';
 if(c==='loss'||c==='loh')return '<span class="b" style="background:#fef3c7;color:#92400e" title="LOH-unmask confound 82-91%">cn '+esc(c)+'·LOH-unmask confound</span>';
 return '<span class="b" style="background:#f1f5f9;color:#475569">cn unknown</span>';
}
// ===== clone 數 c = distinct 全跨 ALT genotype(co-phased) =====
function regionC(r){let g=new Set();(r.lineages||[]).forEach(L=>{const p=L.obs_populations||{};for(const k in p)if(k.indexOf('A')>=0)g.add(k);});return g.size;}
// ===== read/VAF 排序可能結構(CCF pigeonhole·祖先突變 CCF≥衍生·07-09·純遺傳非循環)=====
function ccfBlock(r,L){
 const pos=r.positions||[],trees=L.trees||[];
 if(pos.length<2||trees.length<2)return '';
 const cov=L.obs_col_coverage||{};let ccf={};
 for(let j=0;j<pos.length;j++){const c=cov[String(pos[j])];if(!c)return '';const tot=c[0]+c[1];ccf[j]=tot?c[1]/tot:0;}
 const unlab=s=>{if(s==='ROOT')return[];const s2=(s[0]==='H'&&s[1]==='_')?s.slice(2):s;let b=[];for(let i=0;i<s2.length;i++)if(s2[i]==='A')b.push(i);return b;};
 const MARGIN=0.05,TEMP=0.05;
 function score(edges){let sc=0,vi=0;(edges||[]).forEach(e=>{const P=unlab(e[0]),Ch=unlab(e[1]),acq=Ch.filter(x=>P.indexOf(x)<0);if(acq.length!==1)return;const j=acq[0];P.forEach(anc=>{sc+=(ccf[anc]-ccf[j]);if(ccf[j]>ccf[anc]+MARGIN)vi++;});});return[sc,vi];}
 let scored=trees.map((t,i)=>{const s=score(t.edges);return{i:i,s:s[0],vi:s[1]};});
 const mx=Math.max(...scored.map(o=>o.s));let exps=scored.map(o=>Math.exp((o.s-mx)/TEMP));const tot=exps.reduce((a,b)=>a+b,0)||1;
 scored.forEach((o,k)=>o.post=exps[k]/tot);
 const ranked=[...scored].sort((a,b)=>b.post-a.post),top=ranked[0];
 const ccfrow=pos.map((p,j)=>'S'+(j+1)+' '+(ccf[j]*100).toFixed(0)+'%').join(' · ');
 if(r.cn==='gain'||r.cn==='amp')
  return '<div class="note" style="margin:5px 0;background:var(--danger-bg);border-color:#efb7b0"><span class="evidence-state confounded">Confounded</span> <b>CN-'+esc(r.cn)+'：read/VAF 排序停用</b>。CN 擴增使 read 數不等於 CCF，因此維持等機率、不用 pigeonhole 排序。CCF 僅供參考：'+ccfrow+'。</div>';
 const bars=ranked.map(o=>'樹#'+(o.i+1)+' <b style="color:'+(o.post>=0.6?'#16a34a':'#d97706')+'">'+(o.post*100).toFixed(0)+'%</b>'+(o.vi?'<span style="color:#dc2626">·'+o.vi+'違反</span>':'<span style="color:#16a34a">·0違反</span>')).join(' ｜ ');
 return '<div class="note" style="margin:5px 0"><b>Read/VAF 排序（CCF pigeonhole；祖先突變 CCF ≥ 衍生）</b> '
  +(top.post>=0.6?'<span class="tag t-det">樹 #'+(top.i+1)+' 最一致 '+(top.post*100).toFixed(0)+'%</span>':'<span class="tag t-amb">仍未破等機率 ('+(top.post*100).toFixed(0)+'%)</span>')
  +'<div style="font-size:10.5px;margin-top:3px">此 HP 家族各 sSNV CCF(=<b>該家族</b> ALT 比·含 partial read·與上方 pooled 位點證據不同軸): '+ccfrow+'</div>'
  +'<div style="font-size:10.5px;margin-top:2px">posterior(softmax·TEMP0.05): '+bars+'</div>'
  +'<div class="note" style="background:none;border:none;padding:2px 0"><span class="evidence-state inferred">Heuristic</span> read 多=祖先只是排序提示，不是證明；posterior ≥60% 才標為較一致。此區 cn='+esc(r.cn)+(r.cn==='neutral'?'（CN-clean）':'（非 neutral）')+'。</div></div>';
}
function show(i,row,scrollOnMobile){
 document.querySelectorAll('#list .row').forEach(x=>{x.classList.remove('sel');x.setAttribute('aria-selected','false');});
 if(row){row.classList.add('sel');row.setAttribute('aria-selected','true');}
 const r=R[i];const[rt,rc]=regTag(r.region_determinacy);
 let h='<h3>'+esc(r.region)+' <span class="tag '+rc+'">'+rt+'</span></h3>'
  +'<div class="kv"><span class="b">'+r.n_sSNV+' sSNV</span><span class="b">HP×'+r.hp_multiplicity+(r.is_multiHP?'（多-HP·雙親代各一樹）':'（single-HP）')+'</span>'+cnBadge(r)+(function(){const c=regionC(r);return '<span class="b" style="background:'+(c>=2?'#ede9fe;color:#6d28d9':'#f1f5f9;color:#475569')+'" title="distinct 全跨(co-phased) ALT genotype 數">c='+c+(c>=2?'·subclonal 候選(未確認)':(c===1?'·單 clone':'·無全跨ALT'))+'</span>';})()+'<span class="b">'+r.lineages.length+' germline-HP 家族</span></div>'
  +'<div class="note" style="margin:4px 0">每個 germline-HP 家族分開建樹；使用箭頭切換其等機率最小樹。枚舉後仍無法唯一決定，本身就是分析結果。'
   +'<div style="margin-top:5px;display:flex;gap:12px;flex-wrap:wrap;align-items:center;font-size:10.5px">'
   +'<span><svg width="14" height="14" style="vertical-align:middle"><circle cx="7" cy="7" r="5.5" fill="#2563eb"/></svg> <b>實測</b>(該基因型有 read 觀測到)</span>'
   +'<span><svg width="14" height="14" style="vertical-align:middle"><circle cx="7" cy="7" r="5" fill="#fff" stroke="#a855f7" stroke-width="1.6" stroke-dasharray="3 2"/></svg> <b>推測</b>(隱藏祖先·無 read·solver 補的中間 clone)</span>'
   +'<span><svg width="14" height="14" style="vertical-align:middle"><rect x="2" y="2" width="10" height="10" rx="2" fill="#64748b"/></svg> <b>germline</b> 起點(RR)</span>'
   +'<span><svg width="26" height="10" style="vertical-align:middle"><line x1="1" y1="5" x2="25" y2="5" stroke="#94a3b8" stroke-width="1.3"/></svg> 灰實線=此邊<b>在全部等機率樹一致=forced 骨幹</b></span>'
   +'<span><svg width="26" height="10" style="vertical-align:middle"><line x1="1" y1="5" x2="25" y2="5" stroke="#f59f00" stroke-width="1.9" stroke-dasharray="4 2"/></svg> <b style="color:#d97706">橙虛線</b>=此連接在等機率樹間<b>變動=枚舉組合選擇</b>(非唯一)</span>'
   +'</div><div class="note" style="margin-top:3px;background:none;border:none;padding:0"><span class="evidence-state inferred">Uncertainty</span> 整棵樹全橙表示沒有 forced 骨幹；有灰邊則只有橙色部分仍可變。capped lineage 尚未完整枚舉，forced 標記可能高估。</div>';
 // ===== 位點證據 + 確定性(用原始 col_coverage 真值):每位點實測 ALT reads + co-phase 缺口 =====
 (function(){
  const pos=r.positions||[];if(!pos.length)return;  // 無 mlhp 原始資料則略
  let perpos=pos.map(()=>({ref:0,alt:0}));
  r.lineages.forEach(L=>{const cov=L.obs_col_coverage||{};pos.forEach((p,i)=>{const c=cov[String(p)];if(c){perpos[i].ref+=c[0];perpos[i].alt+=c[1];}});});
  let recur=new Set();r.lineages.forEach(L=>(L.trees||[]).forEach(t=>(t.recurrence||[]).forEach(b=>recur.add(b))));
  let cells='',nObsAlt=0,nZero=0;
  pos.forEach((p,i)=>{const o=perpos[i],tot=o.ref+o.alt,vaf=tot?100*o.alt/tot:0;
   let cls,txt;if(o.alt>0){cls='var(--accent)';txt=o.alt+' ALT observed';nObsAlt++;}else{cls='var(--rec)';txt='0 ALT';nZero++;}
   cells+='<td style="text-align:center;padding:4px 7px"><b>S'+(i+1)+'</b> <span class="note" style="background:none;border:none;padding:0">'+p+'</span>'+(recur.has(i)?' <span style="color:var(--rec);font-size:10px" title="突變兩次·違反 infinite-sites">RECUR</span>':'')+'<br><span style="color:'+cls+';font-size:11px">'+txt+'</span><br><span class="note" style="background:none;border:none;padding:0;font-size:10px">VAF '+vaf.toFixed(0)+'% ('+o.alt+'A/'+o.ref+'R)</span></td>';
  });
  const nfc=r.n_full_cov_reads||0,k=pos.length,gap=pos.length>1?pos[pos.length-1]-pos[0]:0;
  const gtxt=gap>=1000?(gap/1000).toFixed(1)+'kb':gap+'bp';
  const cophaseWeak=nfc<Math.max(2,nObsAlt);  // 跨全位點 read 太少 → 組合推斷
  let cl2={all_determined:['High','var(--det)','全 lineage 唯一樹'],has_ambiguous:['Medium','var(--amb)','有多樹未定'],has_capped:['Low','var(--cap)','枚舉未完（太密）'],has_recurrence:['Low','var(--rec)','有 recurrence'],no_germline_lineage:['—','var(--none)','無 germline']}[r.region_determinacy]||['—','var(--none)',''];
  h+='<div class="note" style="margin:6px 0"><b>結構確定性</b> <span class="tag" style="background:'+cl2[1]+'">'+cl2[0]+'</span> <span class="note" style="background:none;border:none;padding:0">'+cl2[2]+' ｜ '+nObsAlt+'/'+k+' 位點實測到 ALT</span>'
   +'<div style="margin-top:6px;font-weight:700;font-size:12px"><span class="evidence-state observed">Observed</span> 位點證據（原始 col_coverage）</div><table style="margin-top:3px;font-size:11px"><tr>'+cells+'</tr></table>'
   +'<div style="font-size:11px;margin-top:5px"><b>位點層與組合層不同：</b>每個 ALT 可由 read 直接觀測；但 phasing 是否直接實測，取決於 <b>n_full_cov_reads='+nfc+'</b>（同時跨全部 '+k+' 位點的 read）。'
   +(cophaseWeak?'<span class="evidence-state inferred">Inferred</span> 數量不足，沒有 read 把所有 ALT 連成同一分子；位點跨 '+gtxt+'。':'<span class="evidence-state observed">Observed</span> 全跨 read 足以支持組合。')
   +(nZero?'<b style="color:var(--rec)"> 另有 '+nZero+' 個位點為 0 ALT read。</b>':'')+'</div>'
   +(nObsAlt===0?'<div class="note" style="margin-top:5px;background:var(--danger-bg);border-color:#efb7b0;color:#7b1e17"><span class="evidence-state unavailable">Unavailable</span> <b>root-only：沒有突變可重建。</b> 全部 '+k+' 位點皆為 0 ALT read；即使 determinacy='+esc(r.region_determinacy)+'，也不能解讀為真 subclonal 結構。</div>':'')
   +'</div>';
 })();
 h+=locusBlock(r)+pairwiseBlock(r)+neighborBlock(r);
 r.lineages.forEach(L=>{const[ct,cc]=clsTag(L.L1_class);const fc=famCls(L.family);
  const _obsPops=(L.n_full_pops||0);
  h+='<div class="lin f'+L.family+'"><b class="'+fc+'">HP family · '+esc(L.fam_label)+'</b> '
   +'<span class="tag '+cc+'">'+ct+'</span> '
   +'<span class="pill">'+L.n_trees+' 樹'+(L.n_distinct_shapes&&L.n_distinct_shapes<L.n_trees?'/'+L.n_distinct_shapes+'形狀':'')+'</span><span class="pill">'+L.n_hidden+' 隱藏祖先</span>'
   +'<span class="pill">'+(L.n_reads||0)+' reads·'+L.n_full_pops+'full/'+L.n_partial+'partial</span>'
   +(L.verify_pass?'<span class="pill" style="color:var(--det)">V1-7✓</span>':'<span class="pill" style="color:var(--rec)">V✗</span>');
  if(_obsPops===0)h+='<div class="note" style="margin:5px 0;color:#6b2a84;background:#f8efff;border-color:#d9b8ed"><span class="evidence-state inferred">Inferred</span> <b>此家族沒有實測全跨群。</b> n_full_pops=0、只有 '+(L.n_partial||0)+' partial read；下方節點皆為 solver 推測的隱藏祖先，不是直接觀測到的完整 genotype。</div>';
  const _hasTree=(L.trees||[]).some(t=>t.edges&&t.edges.length);
  if(L.trees&&L.trees.length&&_hasTree){const lid='L'+(LID++);const ns=L.trees.length;
   // 穩定邊集=出現在全部 N 棵等機率樹的邊(其餘=枚舉組合選擇,treeSVG 標橙虛線)
   let stable=null;if(ns>1){const cnt={};L.trees.forEach(t=>(t.edges||[]).forEach(e=>{const k=e[0]+'>'+e[1];cnt[k]=(cnt[k]||0)+1;}));stable=new Set(Object.keys(cnt).filter(k=>cnt[k]===ns));}
   h+='<div class="tsw" id="'+lid+'"><div class="tswhead"><b>可能樹結構</b>：'+L.n_trees+' 棵等機率最小樹'
    +(L.n_distinct_shapes?'（'+L.n_distinct_shapes+(L.n_distinct_shapes<ns?'+':'')+' 種形狀）':'')
    +(L.n_trees>ns?' · 顯示前 '+ns:'')
    +(ns>1?'<button type="button" class="tg tnav" aria-label="上一棵等機率樹" onclick="tnav(\''+lid+'\',-1)">←</button><span class="tctr" id="'+lid+'_c">1 / '+ns+'</span><button type="button" class="tg tnav" aria-label="下一棵等機率樹" onclick="tnav(\''+lid+'\',1)">→</button>':'')
    +'</div><div class="tstage">';
   L.trees.forEach((t,ti)=>{
    const nd=new Set();(t.edges||[]).forEach(e=>{nd.add(e[0]);nd.add(e[1]);});
    const nobs=[...nd].filter(n=>n!=='ROOT'&&!(n[0]==='H'&&n[1]==='_')).length;
    const nhid=[...nd].filter(n=>n[0]==='H'&&n[1]==='_').length;
    h+='<div class="tslide" style="display:'+(ti?'none':'block')+'">'+treeSVG(t.edges,stable)
     +'<div class="tcap">樹 #'+(ti+1)+' · <span style="color:#2563eb">'+nobs+' 實測</span> + <span style="color:#a855f7">'+nhid+' 推測祖先</span> + germline'+(t.recurrence&&t.recurrence.length?' · recurrence bit '+t.recurrence.join(','):'')+(stable?' · 橙虛線=與其他等機率樹不同處':'')+'</div></div>';});
   h+='</div>';
   if(ns>1){h+='<div class="thumbs" aria-label="選擇等機率樹">';L.trees.forEach((t,ti)=>{h+='<button type="button" class="thumb'+(ti?'':' on')+'" aria-label="顯示樹 '+(ti+1)+'" onclick="tjump(\''+lid+'\','+ti+')">#'+(ti+1)+'</button>';});h+='</div>';}
   h+='</div>';}
  else{h+='<div class="note" style="margin:5px 0"><b>（此家族無分支樹可畫）</b> '+(_obsPops<=1?'只有單一 genotype 群或僅 germline':'樹 edges 為空')+' → <b>不是缺資料,是本來就無可枚舉的分支結構</b>（'+(L.n_full_pops||0)+' 實測群·'+(L.n_partial||0)+' partial）。determinacy='+esc(L.L1_class)+'。</div>';}
  h+=ccfBlock(r,L);
  h+=readObsBlock(L);
  h+='<div class="trace">'+L.trace.map(t=>'<div>'+esc(t)+'</div>').join('')+'</div></div>';});
 const detail=document.getElementById('detail');detail.innerHTML=h;
 history.replaceState(null,'','#region='+encodeURIComponent(r.region));
 if(scrollOnMobile&&matchMedia('(max-width:900px)').matches){requestAnimationFrame(()=>{detail.scrollIntoView({block:'start'});detail.focus({preventScroll:true});});}
}
// init
(function(){
 renderSummary();
 dash();
 const detLabel={all_determined:'all-determined',has_ambiguous:'ambiguous',has_capped:'capped',has_recurrence:'recurrence',no_germline_lineage:'no germline lineage'};
 const fd=document.getElementById('fdet');Object.keys(C.region_determinacy).forEach(k=>{const o=document.createElement('option');o.value=k;o.textContent=(detLabel[k]||k)+' ('+C.region_determinacy[k]+')';fd.appendChild(o);});
 const chrs=[...new Set(R.map(r=>r.chrom))].sort((a,b)=>(+a.slice(3)||99)-(+b.slice(3)||99));
 const fc=document.getElementById('fchr');chrs.forEach(c=>{const o=document.createElement('option');o.value=c;o.textContent=c;fc.appendChild(o);});
 ['fdet','fhp','fchr','fsort'].forEach(id=>document.getElementById(id).onchange=applyFilter);
 document.getElementById('fq').oninput=applyFilter;
 document.getElementById('copy-link').onclick=function(){const ta=document.createElement('textarea');ta.value=location.href;ta.setAttribute('readonly','');ta.style.position='fixed';ta.style.opacity='0';document.body.appendChild(ta);ta.select();const ok=document.execCommand('copy');ta.remove();this.textContent=ok?'已複製 Region 連結':'請複製網址列';setTimeout(()=>this.textContent='複製目前 Region 連結',1800);};
 applyFilter();
 const hp=new URLSearchParams(location.hash.slice(1));const wanted=hp.get('region');if(wanted){const i=R.findIndex(r=>r.region===wanted);if(i>=0)show(i,null,false);}
})();
</script></body></html>"""

HTML = (HTML.replace("__DATA__", DATA_JSON)
        .replace("__SAMPLE__", SAMPLE)
        .replace("__SAMPLE_OPTIONS__", SAMPLE_OPTIONS)
        .replace("__BACKBONE__", BACKBONE)
        .replace("__GENERATED_AT__", GENERATED_AT))
with open(OUT, "w", encoding="utf-8") as f:
    f.write(HTML)
print(f"OK wrote {OUT} ({len(HTML):,} bytes; {len(regions)} regions)")
