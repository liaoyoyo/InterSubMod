#!/usr/bin/env python3
"""build outlier-tolerant 新規則檢視工作站: 14 案 dual-panel + 新規則指標 + 逐項判讀(localStorage)。
注入 outlier_review_index.json(§13-A)。供用戶確認構想後再大規模。"""
import json, base64, subprocess
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
IDX=json.load(open(f"{A}/outlier_review_index.json"))
BC=subprocess.run(["git","-C","/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p): return base64.b64encode(open(f"{A}/{p}","rb").read()).decode()
CATC={"救回+對齊":"#2ea043","救回+不對齊":"#cc7a22","舊已可分":"#4C9EE0","chr1案":"#b05fa8"}
cards=""
for it in IDX["items"]:
    c=CATC.get(it["cat"],"#888")
    cards+=f'''<div class="card" data-cat="{it['cat']}" data-key="{it['pos']}">
<div class="ch"><span class="cat" style="background:{c}">{it['cat']}</span> <b>{it['pos'].replace('_',':')}</b> n={it['n']}</div>
<div class="cm">舊: {it['old']} → 新: valid群={it['n_valid']} + outlier={it['n_outlier']}({it['outlier_pct']}%) | {it['aligned']}</div>
<img loading="lazy" src="data:image/png;base64,{b64(it['png'])}"/>
<div class="verdict" data-key="{it['pos']}"><button data-v="agree">✓ 構想成立</button><button data-v="doubt">? 存疑</button><button data-v="reject">✗ 不該救/誤判</button>
<select class="reason"><option value="">理由…</option><option>valid群對齊乾淨、outlier確實孤立</option><option>outlier其實是真群</option><option>valid群不對齊=可能雜訊</option><option>少數群被當outlier漏掉</option></select></div></div>'''
CSS="""
:root{--bg:#16181d;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:14px/1.5 system-ui,'Noto Sans CJK TC',sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:20px}h1{font-size:19px}h2{font-size:15px;border-left:3px solid var(--acc);padding-left:8px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px 16px;margin:10px 0}
.legend{font-size:12px;color:var(--mut);background:var(--card);border:1px solid var(--bd);border-radius:6px;padding:8px;margin:8px 0}
.cards{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:10px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px}.card img{width:100%;border-radius:4px;margin-top:6px;background:#fff}
.ch{font-size:13px}.cm{font-size:11.5px;color:var(--mut);margin-top:3px;font-family:ui-monospace,monospace}
.cat{font-size:10px;padding:1px 6px;border-radius:8px;color:#111;font-weight:700}
.verdict{margin-top:7px;display:flex;gap:5px;flex-wrap:wrap}.verdict button{cursor:pointer;border:1px solid var(--bd);background:#262a32;color:var(--txt);border-radius:5px;padding:3px 8px;font-size:12px}
.verdict button.on[data-v=agree]{background:#1a6e3a;border-color:#2ea043}.verdict button.on[data-v=doubt]{background:#7a5a00}.verdict button.on[data-v=reject]{background:#7a1a1a}
.reason{background:#262a32;color:var(--txt);border:1px solid var(--bd);border-radius:5px;font-size:11px}
.toolbar{position:sticky;top:0;background:var(--bg);padding:8px 0;z-index:5;border-bottom:1px solid var(--bd);display:flex;gap:8px;align-items:center;flex-wrap:wrap}
.toolbar button,.toolbar select{background:#262a32;color:var(--txt);border:1px solid var(--bd);border-radius:5px;padding:4px 9px;cursor:pointer}
.foot{color:var(--mut);font-size:11px;margin-top:20px;border-top:1px solid var(--bd);padding-top:10px}code{color:#9ecbff}
"""
JS="""
const LS='outlier_review_v1';
function load(){return JSON.parse(localStorage.getItem(LS)||'{}')}
function save(o){localStorage.setItem(LS,JSON.stringify(o))}
function render(){const st=load();document.querySelectorAll('.verdict').forEach(v=>{const k=v.dataset.key,s=st[k]||{};v.querySelectorAll('button').forEach(b=>b.classList.toggle('on',s.v===b.dataset.v));const r=v.querySelector('.reason');if(r&&s.r!=null)r.value=s.r});cnt()}
document.addEventListener('click',e=>{if(e.target.matches('.verdict button')){const v=e.target.closest('.verdict'),k=v.dataset.key,st=load();st[k]=st[k]||{};st[k].v=e.target.dataset.v;save(st);render()}});
document.addEventListener('change',e=>{if(e.target.matches('.reason')){const v=e.target.closest('.verdict'),k=v.dataset.key,st=load();st[k]=st[k]||{};st[k].r=e.target.value;save(st)}});
function cnt(){const st=load();let a=0,d=0,r=0;Object.values(st).forEach(s=>{if(s.v==='agree')a++;else if(s.v==='doubt')d++;else if(s.v==='reject')r++});document.getElementById('cnt').textContent=`判讀 ✓${a} ?${d} ✗${r}/${document.querySelectorAll('.card').length}`}
function filt(){const c=document.getElementById('fc').value;document.querySelectorAll('.card').forEach(x=>x.style.display=(c==='all'||x.dataset.cat===c)?'':'none')}
function expt(){const st=load();let rows=[['pos','verdict','reason']];document.querySelectorAll('.card').forEach(c=>{const k=c.dataset.key,s=st[k]||{};rows.push([k,s.v||'',s.r||''])});const bl=new Blob([rows.map(r=>r.join(',')).join('\\n')],{type:'text/csv'});const a=document.createElement('a');a.href=URL.createObjectURL(bl);a.download='outlier_review_verdicts.csv';a.click()}
window.onload=render;
"""
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>outlier-tolerant 新規則 — 構想檢視</title><style>{CSS}</style></head><body><div class="wrap">
<h1>outlier-tolerant 切群新規則 — 小規模構想檢視（{IDX['n']} 案）</h1>
<div class="banner"><b>新規則</b>：把「任一群&lt;3 即拒」改成「<b>≥2 個 valid 群(≥3) 即可分；&lt;3 的小群當 outlier 記錄不為準(可能雜訊)，outlier≤20%</b>」。<br>
本頁讓你**逐項肉眼確認構想是否成立**（確認後再大規模）。看點：valid 群是否對齊生物側欄、outlier(灰)是否真的孤立、有無誤判。</div>
<div class="legend">每張：左側欄 <b>new-cut</b>(valid 群=彩色/<b>灰=outlier</b>) · <b>HP</b> · <b>ALT</b> · <b>Strand</b>；左甲基 read×CpG、右 read×read 距離(對角塊=群)。<br>
分類：<span style="color:#2ea043">救回+對齊</span>(好)·<span style="color:#cc7a22">救回+不對齊</span>(precision問題,valid群不對齊基因標籤)·<span style="color:#4C9EE0">舊已可分</span>(對照)·<span style="color:#b05fa8">chr1案</span>(你提的兩案)</div>
<div class="toolbar"><span id="cnt"></span><label>分類 <select id="fc" onchange="filt()"><option value="all">全部</option><option>救回+對齊</option><option>救回+不對齊</option><option>舊已可分</option><option>chr1案</option></select></label><button onclick="expt()">⬇ 匯出判讀</button><span class="muted">(判讀存 localStorage)</span></div>
<div class="cards">{cards}</div>
<div class="foot">build_commit {BC} · 注入自 outlier_review_index.json(§13-A) · 新規則 MIN_SZ=3 / outlier≤20% / 對齊=V&gt;null_95&V≥0.3 · pilot chr21+22+chr1案 · 熱圖 SoT dual-panel · ⭐2</div>
</div><script>{JS}</script></body></html>"""
out=f"{A}/20260620_outlier_tolerant_rule_review_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB) cards={IDX['n']}")
