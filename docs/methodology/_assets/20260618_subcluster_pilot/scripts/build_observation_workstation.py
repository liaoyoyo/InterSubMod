#!/usr/bin/env python3
"""觀察工作站(可持續修正+重複使用)— 完整分佈表 + 各類代表案例(逐項判讀 localStorage+匯出)+ 訊號定義 + changelog。
§13-A: 數字全從 cluster_redesign_wg_records.json + ws_items.json 注入。沿 verify-workstation 設計。"""
import json, os, html, base64
from collections import Counter
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
rec=[r for r in json.load(open(f"{A}/cluster_redesign_wg_records.json")) if "err" not in r]
items=json.load(open(f"{A}/ws_items.json"))
TP=[r for r in rec if r["set"]=="TP"]; FP=[r for r in rec if r["set"]=="FP"]; nTP,nFP=len(TP),len(FP)
def b64(p):
    fp=os.path.join(A,p); return "data:image/png;base64,"+base64.b64encode(open(fp,"rb").read()).decode() if os.path.exists(fp) else None
CLS=[("CONFIRMED","#0d9488","真實+對齊 germline = cis-ASM(somatic-相關)"),
     ("NEAR_CONFIRMED","#0891b2","excess 0.08-0.10 邊界 + 對齊"),
     ("REAL_NOVEL","#7c3aed","真實+大跳但不對齊 = subclone 候選"),
     ("REAL_DIFFUSE","#d97706","真實但無大跳+無對齊 = 散/無法歸因"),
     ("NO_CLEAR_SPLIT","#6b7280","全 k 切不出超過 null 的群 = 真無結構")]
def cnt(g,k): return sum(1 for r in g if r["fine_conf"]==k)
# 分佈
def pct(g,k): return 100*cnt(g,k)/len(g) if g else 0
sigTP=nTP-cnt(TP,"NO_CLEAR_SPLIT"); gmTP=cnt(TP,"CONFIRMED")+cnt(TP,"NEAR_CONFIRMED"); unTP=cnt(TP,"REAL_NOVEL")+cnt(TP,"REAL_DIFFUSE")
# 表1
t1=""
for k,col,desc in CLS:
    t,f=cnt(TP,k),cnt(FP,k); tp,fp=pct(TP,k),pct(FP,k); en=tp/fp if fp>0 else 0
    enc="#0d9488" if en>1.3 else "#db2777" if en<0.8 else "#888"
    t1+=(f"<tr><td><span class='dot' style='background:{col}'></span><b>{k.replace('_SPLIT','')}</b><br><span class=mut>{desc}</span></td>"
         f"<td>{t}</td><td><b>{tp:.1f}%</b></td><td>{f}</td><td>{fp:.1f}%</td><td style='color:{enc};font-weight:700'>{en:.2f}×</td></tr>")
# k 分佈
def kd(g): c=Counter(r["fine_k"] for r in g if r["fine_k"]>=2); return c
kT=kd(TP)
# 案例卡
cards=""
byc={k:[] for k,_,_ in CLS}
for it in items: byc.setdefault(it["fine_conf"],[]).append(it)   # 按重算 verdict 分組(與圖一致)
for k,col,desc in CLS:
    cards+=f"<h3><span class='dot' style='background:{col}'></span>{k.replace('_SPLIT','')} 代表案例 <span class=mut>{html.escape(desc)}</span></h3>"
    for it in byc.get(k,[]):
        img=b64(it["png"]) or ""
        V=it["align_V"]; e=it["align_e"]
        met=(f"n={it['n']} · excess={it['excess']} · gap_ratio={it['gap_ratio']}(big_gap={it['big_gap']}) · "
             f"對齊={it['align_axis']} V={V} e={e} · 群大小={it['core_sizes']} 離群={it['n_outliers']} · "
             f"per-read 核{it['perread']['core']}/邊{it['perread']['edge']}/離{it['perread']['outlier']}")
        cards+=f"""<div class="card" data-id="{html.escape(it['key'])}" data-class="{k}">
          <div class="ch"><span class="cb" style="background:{col}">{k.replace('_SPLIT','')}</span> <b>{html.escape(it['key'])}</b>
            <span class=mut>coarse k={it['coarse_k']}({it['coarse_conf']}) / fine k={it['fine_k']}({it['fine_conf']})</span></div>
          <div class="met">{met}</div>
          <img src="{it['png']}" loading="lazy">
          <div class="judge">人工判讀:
            <button class="jb" data-v="agree">✓ 切割正確</button>
            <button class="jb" data-v="doubt">? 存疑</button>
            <button class="jb" data-v="wrong">✗ 切錯</button>
            <span class="jstate"></span></div>
        </div>"""

HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>切群觀察工作站</title><style>
body{{font-family:system-ui,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif;background:#FAF9F5;color:#141413;margin:0;line-height:1.55}}
.wrap{{max-width:1080px;margin:0 auto;padding:24px 18px}} h1{{font-size:21px}} h2{{font-size:16px;border-left:4px solid #D97757;padding-left:9px;margin-top:26px}} h3{{font-size:14px;margin-top:18px}}
table{{border-collapse:collapse;width:100%;font-size:12.5px;margin:8px 0}} th,td{{border:1px solid #E3DACC;padding:5px 8px;text-align:center}} th{{background:#f0ebe2}} td:first-child{{text-align:left}}
.mut{{color:#888;font-size:11px}} .dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px;vertical-align:middle}}
.box{{border:1px solid #E3DACC;border-radius:8px;padding:11px 15px;background:#fff;margin:10px 0}} .box.key{{border-color:#D97757;background:#fdf6f1}}
.card{{border:1px solid #E3DACC;border-radius:9px;margin:12px 0;background:#fff;overflow:hidden}} .ch{{padding:7px 12px;background:#faf6ef;border-bottom:1px solid #E3DACC;font-size:13px}}
.cb{{display:inline-block;color:#fff;border-radius:5px;padding:1px 7px;font-size:11px;font-weight:600}} .met{{padding:5px 12px;font-size:11px;color:#555;background:#fbfaf7}} .card img{{max-width:100%;display:block}}
.judge{{padding:7px 12px;font-size:12px;border-top:1px solid #eee}} .jb{{margin:0 4px;padding:3px 9px;border:1px solid #ccc;border-radius:5px;background:#fff;cursor:pointer;font-size:12px}}
.jb.on{{background:#0d9488;color:#fff;border-color:#0d9488}} .jstate{{margin-left:8px;color:#888;font-size:11px}}
.bar{{position:sticky;top:0;background:#FAF9F5;border-bottom:1px solid #E3DACC;padding:8px 0;z-index:9}} .bar button{{padding:5px 12px;margin-right:6px;border:1px solid #D97757;background:#fff;border-radius:6px;cursor:pointer}}
details summary{{cursor:pointer;font-weight:600;margin-top:8px}}</style></head><body><div class="wrap">
<div class="bar"><b>切群觀察工作站</b> ·
  <button onclick="exp('json')">匯出 JSON</button><button onclick="exp('csv')">匯出 CSV</button>
  <span id="prog" class="mut"></span></div>

<h1>切群觀察工作站 — HCC1395 tumor-only（可持續修正・重複使用）</h1>
<div class="mut">C++ binary <b>穩定</b>(a678f0e, 0 未編譯 src 異動) · big7 本機全基因組 TP {nTP}+FP {nFP} · 單樣本 ⭐2-3 · 判讀存 localStorage(本機)、可匯出</div>

<div class="box key"><b>訊號定義</b>:位點附近甲基(read×CpG)→ read×read 距離 → UPGMA 切群。<b>「有甲基結構訊號」= 距離矩陣切得出「真實群」</b>:
clusterboot Jaccard(重抽 read 重切)<b>扣 within-1-group null(打散 read 間結構)≥0.10</b>(切群重現且超過無結構假穩基線)。→ NO_CLEAR = 切不出此真實群。</div>

<h2>1. 判別定義（三閘 → 五類）</h2>
<div class="box"><b>① real</b> excess(扣 null)≥0.10 · <b>② big_gap</b> 樹分支跳躍≥max(8×中位,0.4×最大) · <b>③ aligned</b> vs germline CramérV≥0.3 &amp; p&lt;.05 &amp; e≥5</div>
<table><tr><th>類</th><th>①real</th><th>②big_gap</th><th>③aligned</th><th>意義</th></tr>
<tr><td>CONFIRMED</td><td>✓</td><td>–</td><td>✓</td><td>cis-ASM(somatic-相關)</td></tr>
<tr><td>NEAR_CONFIRMED</td><td>近.08-.10</td><td>–</td><td>✓</td><td>邊界 cis-ASM</td></tr>
<tr><td>REAL_NOVEL</td><td>✓</td><td>✓</td><td>✗</td><td>subclone 候選</td></tr>
<tr><td>REAL_DIFFUSE</td><td>✓</td><td>✗</td><td>✗</td><td>真實但無法歸因</td></tr>
<tr><td>NO_CLEAR</td><td>✗(全k)</td><td>–</td><td>–</td><td>真無結構</td></tr></table>

<h2>2. 全基因組分佈（TP {nTP} / FP {nFP}）</h2>
<table><tr><th>類</th><th>TP n</th><th>TP %</th><th>FP n</th><th>FP %</th><th>富集</th></tr>{t1}</table>
<div class="box"><b>甲基訊號層級(TP)</b>:有訊號(非 NO_CLEAR) <b>{sigTP} ({100*sigTP/nTP:.1f}%)</b> ｜
其中 germline 對齊(可信 cis-ASM)<b>{gmTP} ({100*gmTP/nTP:.1f}%)</b> ｜ 真實但未對齊(候選非特異)<b>{unTP} ({100*unTP/nTP:.1f}%)</b> ｜
<b>無法判別/無結構 NO_CLEAR {cnt(TP,'NO_CLEAR_SPLIT')} ({100*cnt(TP,'NO_CLEAR_SPLIT')/nTP:.1f}%)</b>。
切群數分布: {' '.join(f'k{k}:{kT.get(k,0)}' for k in range(2,7))}。</div>

<h2>3. 各類代表案例（逐項確認切割是否正確）</h2>
<div class="mut">每張 = UPGMA 樹 + 原始甲基 read×CpG + 距離 read×read + 側欄(fine|coarse|HP|ALT|Strand),tumor-only。點按鈕記錄判讀。</div>
{cards}

<h2>4. 修正過程紀錄（changelog）</h2>
<details open><summary>方法迭代 + 2 bug(用戶觀察揪出)</summary>
<ul>
<li><b>C++ 狀態</b>: 穩定(a678f0e, binary 與 src 同步, 0 未編譯異動)。</li>
<li><b>bug1</b> stab_excess 用 maxclust(單離群)→ excess=None 誤判 → 改評估實際核心群(chr4 ALT 子群觀察揪出)。</li>
<li><b>bug2</b> NO_CLEAR 把「真實但 diffuse」誤判「無結構」→ 加 REAL_DIFFUSE(「無法切是真無法還是方法問題」揪出)。</li>
<li><b>gap-significance scale-invariant</b>: max(8×中位,0.4×最大)→ chr4 fine k6→k3,防大樹過切。</li>
<li><b>+NEAR_CONFIRMED</b>: excess 0.08-0.10+對齊,收邊界。</li>
<li><b>S2 驗證</b>: big7 跑 decisionflow S2 → 3 NO_CLEAR(正確不誤切)+1 REAL_DIFFUSE(新法較敏感)。</li>
<li>⚠ <b>邊界閾值敏感性</b>: 全基因組用 Rnull=15(提速)、本案例渲染用 Rnull=25 → <b>excess 落在 0.08-0.10 邊界(NEAR_CONFIRMED)的位點會在兩次間 flip</b>(本批 6 個抽樣 NEAR 只 1 個重算維持 NEAR,其餘變 CONFIRMED/REAL_DIFFUSE)。= 邊界本就不確定,卡片按重算 verdict(與圖一致)分組;正式統計以全基因組 records 為準。</li>
</ul></details>

<div class="mut" style="margin-top:24px">資料源 cluster_redesign_wg_records.json + ws_items.json(build a678f0e)。單樣本 characterization 非 subclone 確認。判讀存本機 localStorage,可匯出 JSON/CSV 持續修正。</div>
<script>
const LS="ws_judge_v1";
function load(){{try{{return JSON.parse(localStorage.getItem(LS)||"{{}}")}}catch(e){{return {{}}}}}}
function save(o){{localStorage.setItem(LS,JSON.stringify(o))}}
function prog(){{const o=load();const n=document.querySelectorAll('.card').length;document.getElementById('prog').textContent=`已判讀 ${{Object.keys(o).length}}/${{n}}`}}
document.querySelectorAll('.card').forEach(c=>{{
  const id=c.dataset.id; const o=load(); const cur=o[id];
  c.querySelectorAll('.jb').forEach(b=>{{
    if(cur&&cur.v===b.dataset.v) b.classList.add('on');
    b.onclick=()=>{{const o=load(); o[id]={{v:b.dataset.v,class:c.dataset.class,t:Date.now()}}; save(o);
      c.querySelectorAll('.jb').forEach(x=>x.classList.remove('on')); b.classList.add('on');
      c.querySelector('.jstate').textContent='已記錄'; prog();}};
  }});
  if(cur) c.querySelector('.jstate').textContent='已記錄';
}});
prog();
function exp(fmt){{const o=load(); const rows=Object.entries(o).map(([id,v])=>({{locus:id,class:v.class,verdict:v.v,time:new Date(v.t).toISOString()}}));
  let blob; if(fmt==='json') blob=new Blob([JSON.stringify(rows,null,1)],{{type:'application/json'}});
  else{{const h='locus,class,verdict,time\\n'+rows.map(r=>`${{r.locus}},${{r.class}},${{r.verdict}},${{r.time}}`).join('\\n');blob=new Blob([h],{{type:'text/csv'}})}}
  const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='ws_judge.'+fmt;a.click();}}
</script></div></body></html>"""
out=f"{A}/20260622_observation_workstation_01.standalone.html"; open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024}KB) · {len(items)} 案例 · 有訊號 TP {100*sigTP/nTP:.1f}% · NO_CLEAR {100*cnt(TP,'NO_CLEAR_SPLIT')/nTP:.1f}%")
