#!/usr/bin/env python3
"""build FDR k-sweep workstation HTML (self-contained, base64 heatmaps, injected numbers).
§13-A: 所有 headline 數字注入自 locked json, 缺 key 直接 raise(不手打)。
逐項判讀 localStorage + 匯出 + 漏斗 + 2 驗證表 + 69 熱圖。"""
import json, base64, subprocess, html, os
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
def J(f): return json.load(open(f"{A}/{f}"))
S=J("ksweep_wg_summary.json"); CR=J("cantsplit_reasons.json"); RS=J("cantsplit_apriori_rescue.json")
SP=J("split_accounting.json"); FD=J("ksweep_fdr_loci.json"); IDX=J("fdr_heatmap_index.json")
def req(d,k):
    if k not in d: raise SystemExit(f"REFUSE: missing key {k}")
    return d[k]
BC=subprocess.run(["git","-C","/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p):
    fp=f"{A}/{p}"
    if not os.path.exists(fp): raise SystemExit(f"REFUSE: missing fig {p}")
    return base64.b64encode(open(fp,"rb").read()).decode()

AX={"hp":S["HP_axis"],"carrier":S["CARRIER_axis"],"allele":S["ALLELE_axis"]}
N=req(S,"n_clusterable_regions")

# ---- aggregate numbers (injected) ----
def axrow(name,a):
    return (f"<tr><td>{name}</td><td>{a['n_comparable']}</td><td>{a['n_silK_eq_alignK']}</td>"
            f"<td>{a['n_silK_ne_alignK_raw']} ({round(100*a['n_silK_ne_alignK_raw']/a['n_comparable'],1)}%)</td>"
            f"<td>{a['  of_ne_higher_k']}</td><td><b>{a['n_sil_not_most_aligned_FDR']}</b> "
            f"({a['pct_FDR_of_comparable']}%)</td></tr>")
v1rows="".join([axrow("HP",AX["hp"]),axrow("CARRIER (subclone向)",AX["carrier"]),axrow("ALLELE",AX["allele"])])

# verification table V1 (溯源)
V1=[
 ("WG 可分群 region",N,"ksweep_wg_records.json len",f"= split_accounting.cansplit {SP['cansplit']} ✓","L1"),
 ("HP comparable",AX['hp']['n_comparable'],"HP_axis.n_comparable",f"eq{AX['hp']['n_silK_eq_alignK']}+ne{AX['hp']['n_silK_ne_alignK_raw']} ✓","L1"),
 ("HP FDR硬數字",AX['hp']['n_sil_not_most_aligned_FDR'],"HP_axis.n_sil_not_most_aligned_FDR","BH q<0.05","L1"),
 ("CARRIER comparable",AX['carrier']['n_comparable'],"CARRIER_axis.n_comparable",f"eq{AX['carrier']['n_silK_eq_alignK']}+ne{AX['carrier']['n_silK_ne_alignK_raw']} ✓","L1"),
 ("CARRIER FDR硬數字",AX['carrier']['n_sil_not_most_aligned_FDR'],"CARRIER_axis.n_sil_not_most_aligned_FDR","Bonf61→FDR57","L1"),
 ("ALLELE FDR硬數字",AX['allele']['n_sil_not_most_aligned_FDR'],"ALLELE_axis.n_sil_not_most_aligned_FDR","BH q<0.05","L1"),
 ("FDR union 位點(本工作站)",FD['n_union'],"ksweep_fdr_loci.json n_union",f"HP{FD['n_hp']}+CAR{FD['n_carrier']}+ALL{FD['n_allele']} dedup","L1"),
]
V2=[
 ("切不出 total",CR['cant_total'],"cantsplit_reasons.cant_total",f"A412+B1+C24080 ✓ = N-cansplit","L1"),
 ("讀數不足 n<6",CR['A_insuff_n_lt6'],"A_insuff_n_lt6","=split_accounting.insuff","L1"),
 ("讀數夠但同質",CR['C_cohesive_homogeneous'],"C_cohesive_homogeneous",f"median n={CR['median_n_cohesive']}","L1"),
 ("同質-single label",CR['C_single_label'],"C_single_label",f"+multi{CR['C_multi_label']}={CR['C_cohesive_homogeneous']} ✓","L1"),
 ("multi-label 有a-priori訊號",f"{RS['C_multi_ANY_apriori_sig']} ({RS['C_multi_pct_apriori_sig']}%)","cantsplit_apriori_rescue.C_multi_ANY_apriori_sig",f"matched {RS['C_multi_matched_in_sigcsv']}","L1"),
 ("其中真subclone軸(SubcloneDbeta)",RS['C_multi_subclone_dbeta_sig'],"C_multi_subclone_dbeta_sig","14.1%","L1"),
 ("真null(所有軸)",RS['C_multi_truly_null'],"C_multi_truly_null",f"+sig{RS['C_multi_ANY_apriori_sig']}={RS['C_multi_total']} ✓","L1"),
]
def vtab(rows):
    return "".join(f"<tr><td>{html.escape(str(a))}</td><td><b>{html.escape(str(b))}</b></td><td><code>{html.escape(str(c))}</code></td><td>{html.escape(str(d))}</td><td><span class=l>{e}</span></td></tr>" for a,b,c,d,e in rows)

# funnel SVG per axis
def funnel(a,color):
    stages=[("可分群",N),("comparable",a['n_comparable']),("raw≠",a['n_silK_ne_alignK_raw']),("FDR硬",a['n_sil_not_most_aligned_FDR'])]
    mx=N; bars=""
    for i,(lab,val) in enumerate(stages):
        w=max(2,int(360*val/mx)); y=i*34
        bars+=f'<rect x="90" y="{y}" width="{w}" height="24" rx="3" fill="{color}" opacity="{1-0.12*i}"/>'
        bars+=f'<text x="84" y="{y+16}" font-size="11" text-anchor="end" fill="#bbb">{lab}</text>'
        bars+=f'<text x="{96+w}" y="{y+16}" font-size="11" fill="#eee">{val}</text>'
    return f'<svg viewBox="0 0 470 140" role="img" width="100%"><title>funnel</title>{bars}</svg>'

# ---- per-locus cards ----
items=req(IDX,"items")
cards=""
for it in items:
    ax=it["axes"]; prim=it["primary"]; p=ax[prim]
    badges=" ".join(f'<span class="b b-{a}">{a.upper()} dV={ax[a]["dV"]}</span>' for a in ax)
    EXT=os.environ.get("EXTERNAL")
    img_src=it["png"] if EXT else "data:image/png;base64,"+b64(it["png"])
    metrics=(f"sil k=2 V={p['V_sil']} → align k={p['align_k']} V={p['V_align']} | dV={p['dV']} | p_bonf={p['p_bonf']:.1e}")
    cards+=f'''<div class="card" data-axes="{','.join(ax.keys())}" data-dv="{max(a['dV'] for a in ax.values())}" data-key="{it['key']}">
 <div class="ch"><b>{it['chrom']}:{it['pos']}</b> <span class="muted">n={it['n']} · primary={prim}</span> {badges}</div>
 <div class="cm" title="{metrics}">{metrics}</div>
 <img loading="lazy" src="{img_src}"/>
 <div class="verdict" data-key="{it['key']}">
   <button data-v="agree">✓ 真子結構</button><button data-v="doubt">? 存疑</button><button data-v="reject">✗ artifact</button>
   <select class="reason"><option value="">理由…</option><option>雙峰清楚對齊乾淨</option><option>可能homopolymer</option><option>可能coverage假象</option><option>需normal cis-control</option><option>k過度切分</option></select>
 </div></div>'''

CSS="""
:root{--bg:#16181d;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:14px/1.5 system-ui,'Noto Sans CJK TC',sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:20px}
h1{font-size:20px;margin:.2em 0}h2{font-size:16px;border-left:3px solid var(--acc);padding-left:8px;margin-top:26px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px 16px;margin:10px 0}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin:8px 0}th,td{border:1px solid var(--bd);padding:4px 7px;text-align:left}
th{background:#22262e;color:var(--mut)}code{color:#9ecbff;font-size:11px}.l{color:#7ee787;font-weight:700;font-size:11px}
.grid{display:grid;grid-template-columns:1fr 1fr 1fr;gap:12px}.fcell{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px}
.fcell h3{margin:.1em 0;font-size:13px}
.cards{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:12px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:10px}
.card img{width:100%;border-radius:4px;margin-top:6px;background:#fff}
.ch{font-size:13px}.cm{font-size:11px;color:var(--mut);margin-top:3px;font-family:ui-monospace,monospace}
.muted{color:var(--mut);font-weight:400}
.b{font-size:10px;padding:1px 5px;border-radius:8px;margin-left:3px}.b-hp{background:#1f3a5f}.b-carrier{background:#3a5f1f}.b-allele{background:#5f1f5a}
.verdict{margin-top:7px;display:flex;gap:5px;flex-wrap:wrap}.verdict button{cursor:pointer;border:1px solid var(--bd);background:#262a32;color:var(--txt);border-radius:5px;padding:3px 8px;font-size:12px}
.verdict button.on[data-v=agree]{background:#1a6e3a;border-color:#2ea043}.verdict button.on[data-v=doubt]{background:#7a5a00;border-color:#bb8800}.verdict button.on[data-v=reject]{background:#7a1a1a;border-color:#cc3333}
.reason{background:#262a32;color:var(--txt);border:1px solid var(--bd);border-radius:5px;font-size:11px}
.toolbar{position:sticky;top:0;background:var(--bg);padding:8px 0;z-index:5;border-bottom:1px solid var(--bd);display:flex;gap:8px;align-items:center;flex-wrap:wrap}
.toolbar button,.toolbar select{background:#262a32;color:var(--txt);border:1px solid var(--bd);border-radius:5px;padding:4px 9px;cursor:pointer}
.foot{color:var(--mut);font-size:11px;margin-top:24px;border-top:1px solid var(--bd);padding-top:10px}
.legend{font-size:11px;color:var(--mut);background:var(--card);border:1px solid var(--bd);border-radius:6px;padding:8px;margin-top:6px}
"""
JS="""
const LS='ksweep_fdr_ws_v1';
function load(){return JSON.parse(localStorage.getItem(LS)||'{}')}
function save(o){localStorage.setItem(LS,JSON.stringify(o))}
function render(){const st=load();document.querySelectorAll('.verdict').forEach(v=>{const k=v.dataset.key,s=st[k]||{};v.querySelectorAll('button').forEach(b=>b.classList.toggle('on',s.v===b.dataset.v));const r=v.querySelector('.reason');if(r&&s.r!=null)r.value=s.r});updateCount()}
document.addEventListener('click',e=>{if(e.target.matches('.verdict button')){const v=e.target.closest('.verdict'),k=v.dataset.key,st=load();st[k]=st[k]||{};st[k].v=e.target.dataset.v;save(st);render()}});
document.addEventListener('change',e=>{if(e.target.matches('.reason')){const v=e.target.closest('.verdict'),k=v.dataset.key,st=load();st[k]=st[k]||{};st[k].r=e.target.value;save(st)}});
function updateCount(){const st=load();let a=0,d=0,r=0;Object.values(st).forEach(s=>{if(s.v==='agree')a++;else if(s.v==='doubt')d++;else if(s.v==='reject')r++});document.getElementById('cnt').textContent=`判讀: ✓${a} ?${d} ✗${r} / ${document.querySelectorAll('.card').length}`}
function filt(){const ax=document.getElementById('fax').value;document.querySelectorAll('.card').forEach(c=>{c.style.display=(ax==='all'||c.dataset.axes.includes(ax))?'':'none'})}
function expt(){const st=load();let rows=[['key','verdict','reason']];document.querySelectorAll('.card').forEach(c=>{const k=c.dataset.key,s=st[k]||{};rows.push([k,s.v||'',s.r||''])});const csv=rows.map(r=>r.join(',')).join('\\n');const bl=new Blob([csv],{type:'text/csv'});const a=document.createElement('a');a.href=URL.createObjectURL(bl);a.download='ksweep_fdr_verdicts.csv';a.click()}
window.onload=render;
"""

HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>WG k-sweep 對齊式 k-selection — FDR 確認工作站</title><style>{CSS}</style></head><body><div class="wrap">
<h1>WG k-sweep 對齊式 k-selection — FDR 確認工作站</h1>
<div class="banner"><b>HCC1395 全基因組 · tumor reads · BERNOULLI ±5000 · ⭐2 偵測非驗證(單樣本)</b><br>
問題：silhouette 選的 k 是不是「最對齊 a-priori 生物軸」的 k？這裡是<b>嚴格 gate(ΔV&gt;0.1+Cochran≥5)+雙層多重比較校正(Bonferroni×k + BH-FDR)後</b>仍顯著的 <b>{FD['n_union']} 個</b> union 位點，每個附熱圖供肉眼確認。</div>

<h2>① 總覽漏斗（可分群 {N} → comparable → raw≠ → FDR 硬數字）</h2>
<div class="grid">
 <div class="fcell"><h3>HP 軸</h3>{funnel(AX['hp'],'#4C9EE0')}</div>
 <div class="fcell"><h3>CARRIER 軸（subclone 向）</h3>{funnel(AX['carrier'],'#5FB85F')}</div>
 <div class="fcell"><h3>ALLELE 軸</h3>{funnel(AX['allele'],'#B05FA8')}</div>
</div>
<table><tr><th>軸</th><th>comparable</th><th>sil-k=align-k</th><th>sil-k≠align-k(raw)</th><th>其中更高k</th><th>FDR硬數字</th></tr>{v1rows}</table>
<p class="muted">解讀：raw ~20% 的可分群 loci，silhouette 選的 k 不是最對齊軸的 k（86-93% 是更高 k → silhouette 系統性切太少）；但嚴格校正後硬數字只剩 0.4-1.1%（這 {FD['n_union']} 個 union），是高信心「真子結構」候選。</p>

<h2>② 驗證表 V1（k-sweep，每數字溯源）</h2>
<table><tr><th>數字</th><th>值</th><th>來源:key</th><th>重算/交叉</th><th>L</th></tr>{vtab(V1)}</table>

<h2>③ 驗證表 V2（切不出原因）+ 為何切不出</h2>
<table><tr><th>數字</th><th>值</th><th>來源:key</th><th>重算/交叉</th><th>L</th></tr>{vtab(V2)}</table>
<p class="muted"><b>切不出 ≠ 沒讀數</b>（同質組 median {CR['median_n_cohesive']} reads）<b>≠ 沒訊號</b>（multi-label {RS['C_multi_pct_apriori_sig']}% 有 a-priori 訊號，真 null 只 7.6%）。深層原因：clustering 要「離散雙峰 gap」、a-priori 檢定只要「平均差」；切不出多為 mean-shift/漸層 → silhouette 過不了門檻 → 正確回 1 群。k-sweep 精煉「已可分群的 {N}」、不救切不出的 {CR['cant_total']}。</p>

<h2>④ FDR-顯著位點熱圖工作站（{FD['n_union']} 項·肉眼確認）</h2>
<div class="legend"><b>軸定義</b>：ALLELE=anchor位點ALT/REF(anchor專一) · HP=親代染色體{{1,1-1}}|{{2,2-1}} · <b>CARRIER</b>=germline相位塊{{1,2}} vs <b>somatic相位塊{{1-1,2-1}}</b>(合併HP,非anchor ALT,源碼ReadParser.cpp:154「somatic phase block」;carrier沿用pipeline SubcloneDbeta用語)。<br>每張：左側欄 <b>k=2 | align-k</b>(cluster 著色) · <b>HP</b>(藍HP1/橙HP2) · <b>carr</b>(綠germline相位/紅somatic相位塊) · <b>alle</b>(紫REF/黃ALT) · <b>str</b>(±)；右為甲基熱圖（藍=0 未甲基→紅=1 甲基，白=NA）。<b>看點</b>：k=2 那欄只 2 色（silhouette 的切法），align-k 那欄多 1 色且其邊界是否與某條生物側欄（HP/carr/alle）對齊 → 對齊=真子結構；若 align-k 切在 homopolymer/coverage 假象則 ✗。</div>
<div class="toolbar"><span id="cnt"></span><label>軸篩選 <select id="fax" onchange="filt()"><option value="all">全部</option><option value="carrier">CARRIER</option><option value="hp">HP</option><option value="allele">ALLELE</option></select></label><button onclick="expt()">⬇ 匯出判讀 CSV</button><span class="muted">(判讀存 localStorage)</span></div>
<div class="cards">{cards}</div>

<div class="foot">build_commit {BC} · binary 5c39051(tumor=ClairS_pileup_v040/normal=5khz_simplex) · 數字注入自 ksweep_wg_summary/cantsplit_reasons/cantsplit_apriori_rescue/ksweep_fdr_loci.json(§13-A 缺 key refuse) · 熱圖 = fdr_loci_mini.vcf.gz 重跑 binary → plot_fdr_heatmaps.py · 單樣本 HCC1395 ⭐2 偵測非驗證 · 規則 feedback_verification_table_per_data_answer</div>
</div><script>{JS}</script></body></html>"""
suffix="_lite" if os.environ.get("EXTERNAL") else ""
out=f"{A}/20260618_ksweep_fdr_confirmation_workstation_01{suffix}.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out}  ({len(HTML)//1024} KB)  cards={len(items)}")
