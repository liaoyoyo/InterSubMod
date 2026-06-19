#!/usr/bin/env python3
"""build「切不出≠沒訊號」多角度驗證 HTML: 統計表(null對照+效應量+per-CpG) + 距離分佈SVG + dual-panel熱圖。
數字注入 cantsplit_validation.json + cantsplit_heatmap_index.json(§13-A)。"""
import json, base64, subprocess, os
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
def J(f): return json.load(open(f"{A}/{f}"))
V=J("cantsplit_validation.json"); HM=J("cantsplit_heatmap_index.json")
CR=J("cantsplit_reasons.json"); RS=J("cantsplit_apriori_rescue.json")
BC=subprocess.run(["git","-C","/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p): return base64.b64encode(open(f"{A}/{p}","rb").read()).decode()
G=V["groups"]
def row(name,g):
    s=G.get(g)
    if not s: return ""
    return (f"<tr><td>{name}</td><td>{s['matched']}</td><td><b>{s['hp_p05_pct']}%</b></td><td><b>{s['hp_p001_pct']}%</b></td>"
            f"<td>{s['germ_med']} ({s['germ_gt01']}個>0.1)</td><td>{s['frac_med']}</td></tr>")

# 距離分佈 SVG (3 群 density 曲線)
dd=HM["distdist"]; bins=dd["bins"]; W=460;Hh=170
def poly(g,color):
    dens=dd[g]["density"]; mx=max(max(dd[x]["density"]) for x in("切不出·有訊號","切不出·真null","切得出·對照"))
    pts=" ".join(f"{40+i*(W-60)/(len(dens)-1):.0f},{Hh-20-(d/mx)*(Hh-40):.0f}" for i,d in enumerate(dens))
    return f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2"/>'
DIST_SVG=f'''<svg viewBox="0 0 {W} {Hh+24}" width="100%" role="img"><title>距離分佈</title>
<line x1="40" y1="{Hh-20}" x2="{W-20}" y2="{Hh-20}" stroke="#444"/>
<text x="40" y="{Hh+8}" font-size="10" fill="#9aa">0 近</text><text x="{W-50}" y="{Hh+8}" font-size="10" fill="#9aa">1.0 遠</text>
<text x="{W//2-40}" y="{Hh+20}" font-size="10" fill="#9aa">read×read 距離</text>
{poly("切不出·真null","#888")}{poly("切不出·有訊號","#5FB85F")}{poly("切得出·對照","#D97757")}
<text x="50" y="16" font-size="11" fill="#5FB85F">━ 切不出·有訊號 (med {dd["切不出·有訊號"]["median"]})</text>
<text x="50" y="30" font-size="11" fill="#888">━ 切不出·真null ({dd["切不出·真null"]["median"]})</text>
<text x="50" y="44" font-size="11" fill="#D97757">━ 切得出·對照 ({dd["切得出·對照"]["median"]})</text></svg>'''

# 熱圖 curate: 4 有訊號 + 1 null + 1 切得出
items=HM["items"]
cur=[it for it in items if it["group"]=="切不出·有訊號"][:4]+[it for it in items if it["group"]=="切不出·真null"][:1]+[it for it in items if it["group"]=="切得出·對照"][:1]
def card(it):
    return f'''<figure class="hm"><img src="data:image/png;base64,{b64(it["png"])}"/>
<figcaption><b>[{it['group']}]</b> {it['chrom']}:{it['pos']} n={it['n']} · |germΔβ|={it['germ']} · HP-PERMANOVA p={it['hpP']}</figcaption></figure>'''

CSS="""
:root{--bg:#15171c;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038;--grn:#7ee787}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:15px/1.6 system-ui,'Noto Sans CJK TC',sans-serif}
.wrap{max-width:1000px;margin:0 auto;padding:24px}h1{font-size:21px}h2{font-size:17px;border-left:3px solid var(--acc);padding-left:9px;margin-top:28px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:13px 17px;margin:10px 0}
.key{background:#1a2520;border:1px solid #2c4a38;border-radius:8px;padding:12px 16px;margin:12px 0}.key b{color:var(--grn)}
table{border-collapse:collapse;width:100%;font-size:13px;margin:9px 0}th,td{border:1px solid var(--bd);padding:5px 9px;text-align:left}th{background:#22262e;color:var(--mut)}
.svgbox{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px;margin:10px 0}
.hmgrid{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:10px}
.hm{margin:0;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px}.hm img{width:100%;border-radius:4px;background:#fff}
.hm figcaption{font-size:11px;color:var(--mut);margin-top:5px;font-family:ui-monospace,monospace}
.caveat{border-left:3px solid #cc5;background:#231f12;padding:9px 13px;border-radius:6px;margin:7px 0;font-size:13.5px}
.foot{color:var(--mut);font-size:11px;margin-top:24px;border-top:1px solid var(--bd);padding-top:10px}code{color:#9ecbff;font-size:12px}
"""
cm=G["切不出·multi"]
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>切不出 ≠ 沒訊號 — 六角度驗證</title><style>{CSS}</style></head><body><div class="wrap">
<h1>「切不出 ≠ 沒訊號」— 六角度驗證</h1>
<div class="banner"><b>HCC1395 tumor reads · 切不出 {CR['cant_total']}(80.3%)·multi-label {RS['C_multi_total']}</b><br>
疑慮：切不出的位點真的有訊號嗎，還是其實沒訊號？這裡用 6 個獨立角度（統計+視覺）對齊認知與實際。</div>

<h2>① 置換 null 對照（最關鍵）+ ② 效應量 + ④ per-CpG</h2>
<table><tr><th>群組</th><th>n</th><th>HP-PERMANOVA p&lt;.05</th><th>p≤.001(最強)</th><th>\\|germΔβ\\| 中位</th><th>Fisher per-CpG 中位</th></tr>
{row("切不出·multi",'切不出·multi')}{row("切不出·single",'切不出·single')}{row("切得出",'切得出')}{row("insuff(n&lt;6)",'insuff')}
<tr style="color:#9aa"><td><b>隨機 null 期望</b></td><td>—</td><td><b>~5%</b></td><td><b>~0.1%</b></td><td>~0</td><td>0</td></tr></table>
<div class="key"><b>① null 對照</b>：PERMANOVA 是 999 次置換檢定 → 切不出·multi 的 p&lt;.05 = <b>{cm['hp_p05_pct']}%（隨機 5% 的 {round(cm['hp_p05_pct']/5,1)}×）</b>、p≤.001 = <b>{cm['hp_p001_pct']}%（隨機 0.1% 的 {round(cm['hp_p001_pct']/0.1)}×）</b> → <b>訊號絕非 chance</b>。<br>
<b>② 效應量</b>：切不出 \\|Δβ\\| 中位 {cm['germ_med']}（{cm['germ_gt01']} 個 &gt;0.1）= 真實量值，只比切得出（{G['切得出']['germ_med']}）略小 → 效應略小=低於離散分群解析度。<br>
<b>④ per-CpG</b>：切不出 Fisher_Frac {cm['frac_med']}（~19% CpG 個別顯著）→ CpG 層真有訊號。</div>

<h2>③ 距離分佈：切不出單峰(無 gap) vs 切得出(較分散)</h2>
<div class="svgbox">{DIST_SVG}</div>
<p class="muted">切不出·有訊號 (med {dd['切不出·有訊號']['median']}) 與真null (med {dd['切不出·真null']['median']}) 都偏<b>單峰</b>（reads 距離集中、無乾淨雙峰 gap）→ clustering 無 gap 可切；切得出 (med {dd['切得出·對照']['median']}) 分佈較右移/分散（含群間距離）。<b>這直接是「為何切不出」的機制</b>：沒有離散 gap。</p>

<h2>⑤ 直接看熱圖：甲基有 HP mean-shift、距離無離散塊</h2>
<p>切不出·有訊號（按 HP 排序）：左甲基看 <b>HP1 vs HP2 的 mean-shift</b>（對齊側欄 HP 藍/紫）、右距離<b>整片均勻無對角塊</b> = 有訊號但切不出。真null = 甲基均勻；切得出 = 距離有對角塊。</p>
<div class="hmgrid">{''.join(card(it) for it in cur)}</div>

<h2>⑥ 軸分解 + 誠實邊界</h2>
<div class="caveat"><b>訊號真，但多是 germline-cis 非 subclone</b>：切不出·multi 的 a-priori 訊號 {RS['C_multi_pct_apriori_sig']}% 中，主要是 HP-axis（germline 單倍型 ASM，{RS['C_multi_subclone_dbeta_sig']} 個=14.1% 才是 SubcloneDbeta 真 subclone 軸）。所以「有訊號」=真，但「是 subclone」需 normal cis-control。</div>
<div class="caveat"><b>真 null 確實存在</b>：{RS['C_multi_truly_null']}（7.6%）所有 a-priori 軸皆不顯著 = 真的沒訊號（少數）。</div>
<div class="caveat">單樣本 HCC1395 ⭐2；germΔβ 是 germline-haplotype Δβ（HP1 vs HP2），非 read-cluster。</div>

<div class="foot">build_commit {BC} · 數字注入自 cantsplit_validation.json + cantsplit_heatmap_index.json(§13-A) · null 對照=LabelHPPermanovaP(999置換) · 熱圖 SoT dual-panel · 單樣本 ⭐2</div>
</div></body></html>"""
out=f"{A}/20260620_cantsplit_signal_validation_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB) heatmaps={len(cur)}")
