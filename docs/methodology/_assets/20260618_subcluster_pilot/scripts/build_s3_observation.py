#!/usr/bin/env python3
"""S3 (切不出但 joint 顯著、per-CpG=0) 觀察 HTML — 7 案 dual-panel 確認「modkit 漏、ISM 抓」是真結構還 dispersion 假象。
§13-A 注入 s3_index.json。location-clean(ISM主張真結構) vs dispersion(可能modkit對) 分組 + 逐項判讀。"""
import json, base64, os, subprocess, sys
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
AS=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
idx=json.load(open(f"{AS}/s3_index.json"))
BC=subprocess.run(["git","-C",WT,"rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p):
    fp=f"{AS}/{p}";
    if not os.path.exists(fp): sys.exit(f"REFUSE {fp}")
    return base64.b64encode(open(fp,"rb").read()).decode()
def card(it):
    c="#059669" if it['loc'] else "#d97706"
    tag="location-clean（ISM 主張真結構：dispP≥.05）" if it['loc'] else "dispersion（散開非位移：dispP&lt;.05，可能 modkit 對）"
    return f'''<div class="card"><div class="ch"><span class="tag" style="background:{c}">{tag}</span> <b>{it['pos'].replace('_',':')}</b></div>
<div class="cm">joint {it['axis']} 軸 F={it['F']} p={it['p']} dispP={it['dispP']} · per-CpG 顯著=<b>0</b>（modkit 一無所獲）· n={it['n']}</div>
<img loading="lazy" src="data:image/png;base64,{b64(it['png'])}"/>
<div class="obs">看點：左甲基 read×CpG 是否<b>沒有單一 CpG 欄</b>分群（=per-CpG=0 對）；右 read×read 距離依「{it['axis']}軸」（最左側欄）是否<b>真有塊結構</b>（ISM 對）還是只是散開（modkit 對）？</div>
<div class="verdict" data-key="{it['pos']}"><span>判讀：</span><button data-v="ism">✓ 有真結構(ISM對)</button><button data-v="mod">✗ 只是散開/雜訊(modkit對)</button><button data-v="unsure">? 看不清</button></div></div>'''
loc=[i for i in idx if i['loc']]; dis=[i for i in idx if not i['loc']]
CSS="""
:root{--bg:#16181d;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:14px/1.6 -apple-system,system-ui,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif}
.wrap{max-width:1060px;margin:0 auto;padding:22px}h1{font-size:19px}h2{font-size:15px;border-left:3px solid var(--acc);padding-left:8px;margin-top:24px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:9px;padding:13px 16px;margin:12px 0;font-size:13px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:9px;padding:12px;margin:12px 0}
.card img{width:100%;border-radius:5px;margin:8px 0;background:#fff;border:1px solid var(--bd)}
.ch{font-size:13.5px}.cm{font-size:11.5px;color:var(--mut);margin:3px 0;font-family:ui-monospace,monospace}
.tag{font-size:10px;padding:1px 7px;border-radius:8px;color:#111;font-weight:700}
.obs{font-size:12px;background:#231f1a;border:1px solid #4a3a26;border-radius:6px;padding:7px 10px;margin-top:6px}
.verdict{margin-top:7px;display:flex;gap:6px;align-items:center;flex-wrap:wrap}
.verdict button{cursor:pointer;border:1px solid var(--bd);background:#262a32;color:var(--txt);border-radius:5px;padding:3px 9px;font-size:12px}
.verdict button.on[data-v=ism]{background:#1a6e3a;border-color:#2ea043}.verdict button.on[data-v=mod]{background:#7a4a00}.verdict button.on[data-v=unsure]{background:#374151}
.toolbar{position:sticky;top:0;background:var(--bg);padding:8px 0;border-bottom:1px solid var(--bd);z-index:5}
.toolbar button{background:#262a32;color:var(--txt);border:1px solid var(--bd);border-radius:5px;padding:4px 10px;cursor:pointer}
.foot{color:var(--mut);font-size:11px;margin-top:22px;border-top:1px solid var(--bd);padding-top:10px}
"""
JS="""
const LS='s3_obs_v1';function load(){return JSON.parse(localStorage.getItem(LS)||'{}')}function save(o){localStorage.setItem(LS,JSON.stringify(o))}
function render(){const st=load();document.querySelectorAll('.verdict').forEach(v=>{const s=st[v.dataset.key]||{};v.querySelectorAll('button').forEach(b=>b.classList.toggle('on',s.v===b.dataset.v))});cnt()}
document.addEventListener('click',e=>{if(e.target.matches('.verdict button')){const v=e.target.closest('.verdict'),st=load();st[v.dataset.key]={v:e.target.dataset.v};save(st);render()}});
function cnt(){const st=load();let a=0,m=0;Object.values(st).forEach(s=>{if(s.v==='ism')a++;else if(s.v==='mod')m++});document.getElementById('cnt').textContent=`判讀 ISM對${a} modkit對${m}/${document.querySelectorAll('.verdict').length}`}
function expt(){const st=load();let r=[['pos','verdict']];document.querySelectorAll('.verdict').forEach(v=>r.push([v.dataset.key,(st[v.dataset.key]||{}).v||'']));const b=new Blob([r.map(x=>x.join(',')).join('\\n')],{type:'text/csv'});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='s3_verdicts.csv';a.click()}
window.onload=render;
"""
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>S3 per-CpG=0 但 joint 顯著 — modkit 漏 ISM 抓 確認</title><style>{CSS}</style></head><body><div class="wrap">
<h1>「modkit per-CpG 一無所獲、但 ISM joint 說有結構」— 7 案肉眼確認</h1>
<div class="banner">這些位點 <b>per-CpG 顯著 CpG = 0</b>（modkit 式逐點找不到任何差異位點），但 <b>ISM joint PERMANOVA 顯著</b>。
關鍵分兩種：<span style="color:#34d399">location-clean</span>（dispP≥.05，真位移＝ISM 對、modkit 漏）vs <span style="color:#fbbf24">dispersion</span>（dispP&lt;.05，只是散開＝modkit 不報是對的、ISM 過度宣稱）。
逐一判讀「右距離圖依顯著軸是否真有塊結構」。</div>
<div class="toolbar"><span id="cnt"></span> <button onclick="expt()">⬇ 匯出判讀</button></div>
<h2>location-clean（ISM 主張真結構；{len(loc)} 案）</h2>
{''.join(card(i) for i in loc)}
<h2>dispersion（可能 modkit 對；{len(dis)} 案）</h2>
{''.join(card(i) for i in dis)}
<div class="foot">build {BC} · §13-A 注入 s3_index.json · 全 567 個 S3+per-CpG=0（location-clean 472/dispersion 95）抽樣 · 左甲基RdBu_r 右距離magma · 單樣本 ⭐2-3</div>
</div><script>{JS}</script></body></html>"""
out=f"{AS}/20260620_s3_percpg0_joint_observation_01.standalone.html"
open(out,"w").write(HTML); print(f"WROTE {out} ({len(HTML)//1024} KB)")
