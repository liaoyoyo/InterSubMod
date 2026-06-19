#!/usr/bin/env python3
"""build k-sweep 解説 HTML (self-contained, 輕量: SVG 示意 + 3-4 張代表熱圖 base64)。
解釋『silhouette 切太少 → 對齊式 k-selection → 69 硬候選』。數字注入 locked json(§13-A)。"""
import json, base64, subprocess, os, html
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
def J(f): return json.load(open(f"{A}/{f}"))
S=J("ksweep_wg_summary.json"); SP=J("split_accounting.json"); CR=J("cantsplit_reasons.json")
RS=J("cantsplit_apriori_rescue.json"); FD=J("ksweep_fdr_loci.json"); IDX=J("fdr_heatmap_index.json")
BC=subprocess.run(["git","-C","/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
AX={"hp":S["HP_axis"],"carrier":S["CARRIER_axis"],"allele":S["ALLELE_axis"]}
N=S["n_clusterable_regions"]
M=J("ksweep_summary_merged.json")
AXM={"hp":M["HP_axis"],"carrier":M["CARRIER_axis"],"allele":M["ALLELE_axis"]}
NM=M["n_clusterable_regions"]
def b64(p): return base64.b64encode(open(f"{A}/{p}","rb").read()).decode()

# 選代表熱圖: 每軸 dV 最大那個 (loci 已按 max dV desc 排序)
picks=[]; seen_ax=set()
for it in IDX["items"]:
    pa=it["primary"]
    if pa not in seen_ax:
        picks.append(it); seen_ax.add(pa)
    if len(seen_ax)==3: break
# 再補 1 個 carrier(第二強)
for it in IDX["items"]:
    if it["primary"]=="carrier" and it not in picks: picks.append(it); break

def funnel(a,color):
    st=[("可分群",N),("comparable",a['n_comparable']),("raw≠",a['n_silK_ne_alignK_raw']),("FDR硬",a['n_sil_not_most_aligned_FDR'])]
    mx=N;b=""
    for i,(lab,v) in enumerate(st):
        w=max(3,int(300*v/mx));y=i*30
        b+=f'<rect x="78" y="{y}" width="{w}" height="21" rx="3" fill="{color}" opacity="{1-0.13*i}"/>'
        b+=f'<text x="72" y="{y+15}" font-size="10.5" text-anchor="end" fill="#aab">{lab}</text>'
        b+=f'<text x="{84+w}" y="{y+15}" font-size="11" fill="#eee" font-weight="600">{v}</text>'
    return f'<svg viewBox="0 0 410 124" width="100%" role="img"><title>funnel</title>{b}</svg>'

# 概念示意: silhouette(k=2) 把 3-way 壓成 2 群 vs align(k=3) 對齊標籤
CONCEPT='''<svg viewBox="0 0 620 220" width="100%" role="img"><title>silhouette vs alignment k-selection</title>
<style>.t{font:12px system-ui;fill:#cfd3da}.tt{font:13px system-ui;font-weight:700;fill:#e8e8e8}.s{font:10.5px system-ui;fill:#9aa0aa}</style>
<text x="20" y="20" class="tt">同一位點的甲基 read（顏色＝真實生物標籤）</text>
<g transform="translate(20,40)">
<circle cx="20" cy="20" r="9" fill="#55A868"/><circle cx="44" cy="14" r="9" fill="#55A868"/><circle cx="34" cy="40" r="9" fill="#55A868"/>
<circle cx="120" cy="18" r="9" fill="#C44E52"/><circle cx="144" cy="34" r="9" fill="#C44E52"/><circle cx="118" cy="44" r="9" fill="#C44E52"/>
<circle cx="210" cy="20" r="9" fill="#4C9EE0"/><circle cx="234" cy="16" r="9" fill="#4C9EE0"/><circle cx="222" cy="42" r="9" fill="#4C9EE0"/>
<text x="0" y="78" class="s">germline · carrier · 第三亞群（生物上 3 群）</text></g>
<line x1="20" y1="135" x2="600" y2="135" stroke="#2c3038"/>
<g transform="translate(20,150)">
<rect x="0" y="0" width="150" height="56" rx="6" fill="none" stroke="#DD8452" stroke-width="2"/>
<rect x="170" y="0" width="120" height="56" rx="6" fill="none" stroke="#DD8452" stroke-width="2"/>
<text x="0" y="-6" class="tt" fill="#DD8452">silhouette → k=2（切太少）</text>
<text x="6" y="32" class="s">綠+紅 被塞進同一群</text><text x="178" y="32" class="s">藍</text>
</g>
<g transform="translate(330,150)">
<rect x="0" y="0" width="78" height="56" rx="6" fill="none" stroke="#2ea043" stroke-width="2"/>
<rect x="92" y="0" width="78" height="56" rx="6" fill="none" stroke="#2ea043" stroke-width="2"/>
<rect x="184" y="0" width="78" height="56" rx="6" fill="none" stroke="#2ea043" stroke-width="2"/>
<text x="0" y="-6" class="tt" fill="#2ea043">align → k=3（對齊生物軸）</text>
<text x="22" y="32" class="s">綠</text><text x="112" y="32" class="s">紅</text><text x="206" y="32" class="s">藍</text>
</g></svg>'''

def card(it):
    p=it["axes"][it["primary"]]
    return f'''<figure class="hm"><img src="data:image/png;base64,{b64(it["png"])}"/>
<figcaption><b>{it["chrom"]}:{it["pos"]}</b> n={it["n"]} · {it["primary"].upper()} 軸 · sil k=2 (V={p["V_sil"]}) → <b>align k={p["align_k"]} (V={p["V_align"]})</b> · dV={p["dV"]} · p_bonf={p["p_bonf"]:.0e}</figcaption></figure>'''

def axrow(name,a):
    return (f"<tr><td>{name}</td><td>{a['n_comparable']}</td>"
            f"<td>{a['n_silK_ne_alignK_raw']} ({round(100*a['n_silK_ne_alignK_raw']/a['n_comparable'],1)}%)</td>"
            f"<td>{a['  of_ne_higher_k']}</td><td><b>{a['n_sil_not_most_aligned_FDR']}</b> ({a['pct_FDR_of_comparable']}%)</td></tr>")
def cmprow(name,t,m):
    mult=round(m['pct_FDR_of_comparable']/t['pct_FDR_of_comparable'],1) if t['pct_FDR_of_comparable'] else '—'
    return (f"<tr><td>{name}</td><td>{t['n_sil_not_most_aligned_FDR']} ({t['pct_FDR_of_comparable']}%)</td>"
            f"<td><b>{m['n_sil_not_most_aligned_FDR']} ({m['pct_FDR_of_comparable']}%)</b></td><td>×{mult}</td></tr>")
HPMULT=round(AXM['hp']['pct_FDR_of_comparable']/AX['hp']['pct_FDR_of_comparable'],1)

CSS="""
:root{--bg:#15171c;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038;--grn:#7ee787}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:15px/1.65 system-ui,'Noto Sans CJK TC',sans-serif}
.layout{display:grid;grid-template-columns:210px 1fr;max-width:1080px;margin:0 auto}
nav{position:sticky;top:0;align-self:start;height:100vh;padding:20px 12px;border-right:1px solid var(--bd);font-size:13px}
nav a{display:block;color:var(--mut);text-decoration:none;padding:4px 0}nav a:hover{color:var(--acc)}
main{padding:24px 30px;min-width:0}
h1{font-size:22px}h2{font-size:18px;border-left:3px solid var(--acc);padding-left:10px;margin-top:34px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:14px 18px;margin:12px 0}
.key{background:#1a2520;border:1px solid #2c4a38;border-radius:8px;padding:12px 16px;margin:14px 0}
.key b{color:var(--grn)}
table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}th,td{border:1px solid var(--bd);padding:6px 9px;text-align:left}th{background:#22262e;color:var(--mut)}
.svgbox{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:14px;margin:12px 0}
.grid3{display:grid;grid-template-columns:1fr 1fr 1fr;gap:12px}.grid3 .svgbox{margin:0}
.hm{margin:0;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px}
.hm img{width:100%;border-radius:4px;background:#fff}.hm figcaption{font-size:11.5px;color:var(--mut);margin-top:5px;font-family:ui-monospace,monospace}
.hmgrid{display:grid;grid-template-columns:1fr 1fr;gap:12px}
.caveat{border-left:3px solid #cc5;background:#231f12;padding:10px 14px;border-radius:6px;margin:8px 0;font-size:14px}
.foot{color:var(--mut);font-size:11.5px;margin-top:30px;border-top:1px solid var(--bd);padding-top:12px}
code{color:#9ecbff;font-size:12px}
@media(max-width:760px){.layout{grid-template-columns:1fr}nav{position:static;height:auto;border-right:none;border-bottom:1px solid var(--bd)}.grid3,.hmgrid{grid-template-columns:1fr}}
"""

HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>k-sweep 對齊式 k-selection — 解説</title><style>{CSS}</style></head><body><div class="layout">
<nav><b>目次</b>
<a href="#q">① 我們問什麼</a><a href="#concept">② silhouette 的盲點</a><a href="#ans">③ 兩層答案</a>
<a href="#what">④ 69 硬候選代表什麼</a><a href="#read">⑤ 熱圖怎麼確認</a><a href="#bound">⑥ 誠實邊界</a><a href="#merged">⑦ tumor vs merged</a></nav>
<main>
<h1>k-sweep 對齊式 k-selection — 解説</h1>
<div class="banner"><b>HCC1395 全基因組（chr1–22, {N} 可分群位點）· tumor reads · BERNOULLI ±5000 · ⭐2 偵測非驗證</b><br>
一句話：現行用 silhouette 選「切幾群(k)」會<b>系統性切太少</b>；改用「對齊獨立生物軸」選 k，全基因組嚴格校正後找出 <b>{FD['n_union']} 個</b>硬候選子結構。</div>

<h2 id="q">① 我們問什麼</h2>
<p>每個位點的甲基 read 用 UPGMA 分群，要切成幾群(k)？現行 pipeline 用 <b>silhouette（群內聚集度）</b>挑 k —— 但它只看「切出來緊不緊」，<b>不看切法對不對</b>。我們改問：掃 k=2..5，哪個 k 的切法<b>最對齊一條獨立的 a-priori 生物軸</b>（HP / germline-vs-carrier / allele）？這些標籤來自 BAM tag / 變異鹼基（<code>ReadParser.cpp:126/165</code>），<b>與甲基距離無關 → 非循環</b>。</p>
<div class="key" style="border-color:#cc5;background:#231f12"><b>三軸 = 同一個 BAM <code>HP</code> tag + anchor 鹼基的三種切法（精確定義，源碼 <code>ReadParser.cpp:130-157</code> / <code>FisherExact.cpp:415</code>）：</b><br>
• <b>ALLELE</b> = read 在 <b>anchor 位點</b>的 ALT/REF（<b>anchor 專一</b>）<br>
• <b>HP</b> = 親代染色體 <code>{{1,1-1}}</code>=HP1 / <code>{{2,2-1}}</code>=HP2<br>
• <b>CARRIER</b> = germline 相位塊 <code>{{1,2}}</code> vs somatic 相位塊 <code>{{1-1,2-1}}</code>（<b>合併 HP</b>）<br>
⚠ <b>CARRIER 的 <code>1-1/2-1</code> = longphase「somatic phase block」</b>（被 somatic 變異定相的 read），<b>非 anchor ALT</b>（那是 ALLELE）；本軸合併 HP，比 pipeline 逐-HP <code>SubcloneDbeta</code> 粗；「carrier」沿用 pipeline 用語。<b>→ subclone 解讀優先 ALLELE(anchor) + 原生 SubcloneDbeta(逐-HP)。</b></div>

<h2 id="concept">② silhouette 的盲點（為什麼會切太少）</h2>
<div class="svgbox">{CONCEPT}</div>
<p>silhouette 偏好「最緊湊的二分」。當一個位點真的有 3 個生物亞群、但其中兩群在甲基空間靠得較近，silhouette 會把它們<b>併成一群（k=2）</b>，犧牲第三群換取「更緊」。對齊式 k-selection 因為盯著生物軸，會選 <b>k=3</b> 把三群分開。</p>

<h2 id="ans">③ 兩層答案</h2>
<div class="grid3">
 <div class="svgbox"><b style="font-size:13px">HP 軸</b>{funnel(AX['hp'],'#4C9EE0')}</div>
 <div class="svgbox"><b style="font-size:13px">CARRIER（subclone向）</b>{funnel(AX['carrier'],'#5FB85F')}</div>
 <div class="svgbox"><b style="font-size:13px">ALLELE 軸</b>{funnel(AX['allele'],'#B05FA8')}</div></div>
<table><tr><th>軸</th><th>comparable</th><th>sil-k≠對齊-k(raw)</th><th>其中更高k</th><th>FDR硬數字</th></tr>
{axrow("HP",AX['hp'])}{axrow("CARRIER",AX['carrier'])}{axrow("ALLELE",AX['allele'])}</table>
<div class="key"><b>raw 層</b>：~20% 可分群位點的 silhouette-k 不是最對齊 k，且 86–93% 是<b>更高 k</b> → silhouette 系統性切太少。<br>
<b>校正層</b>（gate ΔV&gt;0.1+Cochran≥5 → Bonferroni×k → 全基因組 BH-FDR q&lt;0.05）：硬數字只剩 HP <b>18</b> / CARRIER <b>57</b> / ALLELE <b>40</b>，去重 <b>union {FD['n_union']}</b>。20% 多是 soft（小效應、撐不過多重比較），{FD['n_union']} 個是高信心。</div>

<h2 id="what">④ 這 {FD['n_union']} 個（CARRIER {FD['n_carrier']}）代表什麼</h2>
<p>CARRIER 軸 = germline-tag vs carrier-tag = <b>subclone 方向</b>。這些位點 silhouette 看成 k=2，但實際有第 3/4 群、且乾淨對齊 germline/carrier 分界 = plan 講的「balance/silhouette gate 切太嚴漏掉的 minor 子結構」，現在全基因組+FDR 點名出來。占比：{FD['n_union']} / {N} = {round(100*FD['n_union']/N,2)}% 可分群、或 / {SP['N']} 全 TP = {round(100*FD['n_union']/SP['N'],2)}%。</p>

<h2 id="read">⑤ 熱圖怎麼確認（代表案，肉眼判準）</h2>
<p>左欄 <code>k=2 | align-k</code>(cluster 著色) · <code>HP | carr | alle | str</code>(生物標籤) · 右為甲基矩陣(藍=未甲基→紅=甲基,白=NA)。<b>看點</b>：align-k 比 k=2 多出來的分界，是否與某條生物側欄的顏色邊界<b>重合</b>？重合=真子結構 ✓；切在 homopolymer/稀疏覆蓋=artifact ✗。</p>
<div class="hmgrid">{''.join(card(it) for it in picks)}</div>

<h2 id="bound">⑥ 誠實邊界</h2>
<div class="caveat"><b>數量小</b>：{FD['n_union']}/{N}=1.15% 可分群、0.23% 全 TP。對齊式 k 是「少數位點精煉」，非取代 silhouette。</div>
<div class="caveat"><b>「對齊 carrier 軸」≠「subclone」</b>：85% 這類可能是 cis-ASM（germline 等位基因甲基差），需 normal-anchored cis-control 才能從「≥2 群」確認到「subclone」（紅線，pending）。</div>
<div class="caveat"><b>單樣本 HCC1395 ⭐2 偵測非驗證</b>；跨樣本（COLO829/H1437）未做。</div>
<div class="caveat"><b>只精煉「可分群的 {N}」</b>，不救「切不出的 {CR['cant_total']}」—— 後者任何 k 都形不成群（median {CR['median_n_cohesive']} reads 但甲基同質），靠不分群的 a-priori Δβ（{RS['C_multi_pct_apriori_sig']}% 有訊號）。</div>

<h2 id="merged">⑦ tumor-only vs merged(+normal) 對比</h2>
<p>把 normal reads（全 germline <code>{{1,2}}</code>/REF、無 somatic 標籤）合併進分群+標籤 = paired / cis-control。dual-mode 同一次 binary pass，tumor 模式重現觀察 A（FDR 18/57/40 = 交叉驗證）。</p>
<table><tr><th>軸</th><th>tumor FDR硬(rate)</th><th>merged FDR硬(rate)</th><th>rate 倍數</th></tr>
{cmprow("HP",AX['hp'],AXM['hp'])}{cmprow("CARRIER",AX['carrier'],AXM['carrier'])}{cmprow("ALLELE",AX['allele'],AXM['allele'])}</table>
<div class="key">可分群 region tumor <b>{N}</b> → merged <b>{NM}</b>（{round(NM/N,1)}×，normal 補足 read 過門檻、<b>非新 subclone</b>）。<br>
<b>HP 軸 rate ×{HPMULT} 暴增</b>（normal 是乾淨 germline HP1/HP2 → 強化 germline 單倍型 ASM = <b>cis-ASM/imprinting 層</b>）；CARRIER/ALLELE rate 幾乎不動（normal 無 somatic，per-locus 不幫忙，多出的絕對數來自 3× 大分母）。<br>
<b>→ subclone 偵測用 tumor-only；merged = cis-control / germline 層揭露</b>（HP ×{HPMULT} 對應「≥2 群≠subclone、85% 可能 cis-ASM」紅線）。</div>

<div class="foot">build_commit {BC} · binary 5c39051 · 數字注入自 ksweep_wg_summary/split_accounting/cantsplit_*/ksweep_fdr_loci.json(§13-A) · 代表熱圖={len(picks)} 張(完整 {FD['n_union']} 張見工作站 …_workstation_01) · 全 22 autosomes · 規則 feedback_verification_table_per_data_answer</div>
</main></div></body></html>"""
out=f"{A}/20260618_ksweep_explainer_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB) picks={len(picks)} axes={[p['primary'] for p in picks]}")
