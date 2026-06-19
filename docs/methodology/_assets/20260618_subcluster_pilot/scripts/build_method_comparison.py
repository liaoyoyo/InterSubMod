#!/usr/bin/env python3
"""build Fisher+V vs PERMANOVA 方法對照 HTML。數字注入 method_comparison.json(§13-A)。
3-情況示意 SVG + Q1分軸/多組 + Q2 dispersion + Q3 Venn + Q4 Δβ + 4組 HP-fine。"""
import json, base64, subprocess, glob
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
M=json.load(open(f"{A}/method_comparison.json"))
BC=subprocess.run(["git","-C","/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p):
    g=glob.glob(f"{A}/{p}")
    return base64.b64encode(open(g[0],"rb").read()).decode() if g else ""
fv=M["fisher_v"];pr=M["permanova_raw"];dp=M["dispersion"];vr=M["venn_raw"];vc=M["venn_clean"];db=M["delta_by_class"];fg=M["fourgroup"]
hm_grad=b64("figs_cantsplit/cs_切不出·有訊號_chr12_28414696.png")
hm_disc=b64("figs_fdr/fdr_chr8_131074695.png")

# 3-情況示意
CONCEPT='''<svg viewBox="0 0 640 250" width="100%" role="img"><title>三情況</title>
<style>.tt{font:12px system-ui;font-weight:700;fill:#e8e8e8}.s{font:10.5px system-ui;fill:#9aa0aa}.v{font:11px system-ui;font-weight:700}</style>
<text x="10" y="16" class="s">read 在甲基軸 β（藍=HP1 標籤 / 橙=HP2 標籤）。問: 兩標籤能不能被「看出不同」?</text>
<g transform="translate(10,30)"><text x="0" y="10" class="tt">①離散(有gap)</text>
<line x1="0" y1="34" x2="280" y2="34" stroke="#444"/>
<circle cx="35" cy="34" r="6" fill="#4C9EE0"/><circle cx="55" cy="34" r="6" fill="#4C9EE0"/><circle cx="48" cy="34" r="6" fill="#4C9EE0"/>
<circle cx="225" cy="34" r="6" fill="#DD8452"/><circle cx="245" cy="34" r="6" fill="#DD8452"/><circle cx="235" cy="34" r="6" fill="#DD8452"/>
<text x="120" y="38" class="s">[空隙]</text>
<text x="0" y="58" class="v" fill="#7ee787">Fisher+V ✓</text><text x="90" y="58" class="v" fill="#7ee787">PERMANOVA ✓</text></g>
<g transform="translate(340,30)"><text x="0" y="10" class="tt">②漸層(無gap,平均差)</text>
<line x1="0" y1="34" x2="280" y2="34" stroke="#444"/>
<circle cx="50" cy="34" r="6" fill="#4C9EE0"/><circle cx="90" cy="34" r="6" fill="#4C9EE0"/><circle cx="130" cy="34" r="6" fill="#4C9EE0"/><circle cx="110" cy="34" r="6" fill="#DD8452"/>
<circle cx="160" cy="34" r="6" fill="#DD8452"/><circle cx="200" cy="34" r="6" fill="#DD8452"/><circle cx="150" cy="34" r="6" fill="#4C9EE0"/>
<text x="0" y="58" class="v" fill="#e06c6c">Fisher+V ✗</text><text x="90" y="58" class="v" fill="#7ee787">PERMANOVA ✓</text></g>
<g transform="translate(10,130)"><text x="0" y="10" class="tt">③無訊號(平均也一樣)</text>
<line x1="0" y1="34" x2="280" y2="34" stroke="#444"/>
<circle cx="60" cy="34" r="6" fill="#4C9EE0"/><circle cx="100" cy="34" r="6" fill="#DD8452"/><circle cx="140" cy="34" r="6" fill="#4C9EE0"/><circle cx="180" cy="34" r="6" fill="#DD8452"/><circle cx="120" cy="34" r="6" fill="#4C9EE0"/><circle cx="160" cy="34" r="6" fill="#DD8452"/>
<text x="0" y="58" class="v" fill="#e06c6c">Fisher+V ✗</text><text x="90" y="58" class="v" fill="#e06c6c">PERMANOVA ✗</text></g>
<g transform="translate(340,130)"><text x="0" y="10" class="tt">關鍵</text>
<text x="0" y="32" class="s">Fisher+V 問: 能畫一條 gap 把群切開並對齊標籤?</text>
<text x="0" y="50" class="s">PERMANOVA 問: 標籤之間距離有沒有差(含平均+散開)?</text>
<text x="0" y="72" class="s" fill="#cc5">②漸層 = 只有 PERMANOVA 抓得到 → 20% vs 85% 落差來源</text></g></svg>'''

# Venn SVG (A=Fisher+V cansplit / C=clean-location)
VENN=f'''<svg viewBox="0 0 460 200" width="100%" role="img"><title>Venn</title>
<circle cx="180" cy="100" r="78" fill="#4C9EE0" fill-opacity="0.28" stroke="#4C9EE0"/>
<circle cx="280" cy="100" r="92" fill="#5FB85F" fill-opacity="0.28" stroke="#5FB85F"/>
<text x="120" y="50" font-size="12" fill="#9cd">Fisher+V cansplit</text><text x="120" y="65" font-size="12" fill="#9cd">{fv['cansplit']} ({fv['pct']}%)</text>
<text x="330" y="40" font-size="12" fill="#9d9" text-anchor="end">PERMANOVA clean-loc</text><text x="330" y="55" font-size="12" fill="#9d9" text-anchor="end">{vc['C']} ({vc['C_pct']}%)</text>
<text x="130" y="105" font-size="13" fill="#eee" text-anchor="middle">A−C<tspan x="130" dy="16">{vc['A_only']}</tspan><tspan x="130" dy="14" font-size="10">({vc['A_only_pct']}%)</tspan></text>
<text x="232" y="100" font-size="13" font-weight="700" fill="#fff" text-anchor="middle">∩<tspan x="232" dy="16">{vc['AnC']}</tspan><tspan x="232" dy="14" font-size="10">({vc['AnC_pct']}%)</tspan></text>
<text x="330" y="105" font-size="13" fill="#eee" text-anchor="middle">C−A<tspan x="330" dy="16">{vc['C_only']}</tspan><tspan x="330" dy="14" font-size="10">({vc['C_only_pct']}%)</tspan></text>
<text x="230" y="195" font-size="12" fill="#cc5" text-anchor="middle">Jaccard={vc['jaccard']} → 數量相近但位點大不同(正交兩鏡頭)</text></svg>'''

CSS="""
:root{--bg:#15171c;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038;--grn:#7ee787}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:15px/1.6 system-ui,'Noto Sans CJK TC',sans-serif}
.wrap{max-width:1000px;margin:0 auto;padding:24px}h1{font-size:21px}h2{font-size:17px;border-left:3px solid var(--acc);padding-left:9px;margin-top:28px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:13px 17px;margin:10px 0}
.key{background:#1a2520;border:1px solid #2c4a38;border-radius:8px;padding:12px 16px;margin:12px 0}.key b{color:var(--grn)}
.warn{background:#231f12;border:1px solid #4a432c;border-left:3px solid #cc5;border-radius:6px;padding:10px 14px;margin:10px 0}
table{border-collapse:collapse;width:100%;font-size:13px;margin:9px 0}th,td{border:1px solid var(--bd);padding:5px 9px;text-align:left}th{background:#22262e;color:var(--mut)}
.svgbox{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:12px;margin:10px 0}
.hmgrid{display:grid;grid-template-columns:1fr 1fr;gap:12px}.hm{margin:0;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px}.hm img{width:100%;border-radius:4px;background:#fff}.hm figcaption{font-size:11px;color:var(--mut);margin-top:5px;font-family:ui-monospace,monospace}
.foot{color:var(--mut);font-size:11px;margin-top:24px;border-top:1px solid var(--bd);padding-top:10px}code{color:#9ecbff;font-size:12px}
"""
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Fisher+V vs PERMANOVA — 方法對照</title><style>{CSS}</style></head><body><div class="wrap">
<h1>Fisher+V vs PERMANOVA — 方法對照（驗證甲基-標籤關係的兩種鏡頭）</h1>
<div class="banner"><b>HCC1395 全基因組 N={M['N']} · tumor reads · ⭐2</b><br>
兩方法測「標籤之間甲基有沒有差」用不同問法 → 找到的數量(20% vs 85%)與位點都不同。本頁逐項對照 + 校正常見誤解。</div>

<h2>① 核心概念：兩種問法（三情況示意）</h2>
<div class="svgbox">{CONCEPT}</div>
<table><tr><th></th><th>Fisher's exact + Cramér's V</th><th>PERMANOVA</th></tr>
<tr><td>作用對象</td><td><b>無監督群 × 標籤</b> 列聯表</td><td><b>read×read 距離矩陣</b> 按標籤分組</td></tr>
<tr><td>問什麼</td><td>能畫 gap 切群、且對齊標籤?</td><td>標籤間距離有沒有差(平均+散開)?</td></tr>
<tr><td>前提</td><td>要先 clustering 出離散群</td><td>只要 a-priori 標籤(不分群)</td></tr>
<tr><td>找到</td><td><b>{fv['pct']}%</b>(離散可分群)</td><td>HP <b>{pr['hp_pct']}%</b>/Allele {pr['al_pct']}%/任一 {pr['any_pct']}%</td></tr></table>

<h2>② Q1：PERMANOVA 分軸 + 多組怎麼測（源碼）</h2>
<table><tr><th>問題</th><th>答（源碼 SignificanceAnalyzer.cpp:288-304）</th></tr>
<tr><td>分不同軸?</td><td><b>是</b>：HP 軸 + Allele 軸各一個獨立 PERMANOVA(+PERMDISP)</td></tr>
<tr><td>HP 多組?</td><td><b>合併成 2 組</b>(HP1-fam{{1,1-1}} vs HP2-fam{{2,2-1}})，非 4 組</td></tr>
<tr><td>兩兩確認?</td><td><b>否</b>：PERMANOVA omnibus(多組一次測 pseudo-F)；4組 HP-fine 改用 Fisher-FH</td></tr></table>

<h2>③ Q2：PERMANOVA 的 85% 大部分是 dispersion（散開），不是平均差 🔴</h2>
<table><tr><th>軸</th><th>PERMANOVA sig</th><th>其中 dispersion 混淆</th><th>真 location(平均差)</th></tr>
<tr><td>HP</td><td>{dp['hp_sig']} ({pr['hp_pct']}%)</td><td><b>{dp['hp_disp']} ({dp['hp_disp_pct']}%)</b></td><td><b>{dp['hp_clean']} ({dp['hp_clean_pct']}% = {dp['hp_clean_of_all']}% 全位點)</b></td></tr>
<tr><td>Allele</td><td>{pr['al']} ({pr['al_pct']}%)</td><td><b>{dp['al_disp_pct']}%</b></td><td>{round(100-dp['al_disp_pct'],1)}%</td></tr></table>
<div class="warn"><b>修正常見誤解</b>：PERMANOVA「85% 有訊號」**~72% 是 dispersion(散開度差異)**，扣掉後 **HP clean-location 只 {dp['hp_clean_of_all']}%**。Allele 軸更誇張(93% 是 dispersion)。需 **PERMDISP** 分辨 location vs dispersion。</div>

<h2>④ Q3：Fisher+V vs PERMANOVA 的 Venn（位點交集差集）</h2>
<table><tr><th>對 raw PERMANOVA</th><th>數量</th><th>%</th></tr>
<tr><td>A∩B(都有)</td><td>{vr['AnB']}</td><td>{vr['AnB_pct']}%</td></tr>
<tr><td>A−B(切群但 PERMANOVA 不顯著)</td><td>{vr['A_only']}</td><td>0.1%</td></tr>
<tr><td>B−A(PERMANOVA 有但切不出群)</td><td>{vr['B_only']}</td><td>{vr['B_only_pct']}%</td></tr>
<tr><td>neither</td><td>{vr['neither']}</td><td>{vr['neither_pct']}%</td></tr></table>
<div class="svgbox">{VENN}</div>
<div class="warn"><b>🔴 關鍵更正</b>：Fisher+V({fv['pct']}%) 與 PERMANOVA clean-location({vc['C_pct']}%) **數量相近但位點大不同 —— Jaccard 僅 {vc['jaccard']}、只重疊 {vc['AnC_pct']}%**。因為一個測「能否切離散群」、一個測「標籤平均差」，是<b>正交兩鏡頭、不可互換</b>。</div>

<h2>⑤ Q4：PERMANOVA vs「Δβ(read 平均→標籤內再平均)」</h2>
<table><tr><th></th><th>PERMANOVA</th><th>Δβ(GermlineAsmDbeta)</th></tr>
<tr><td>作用</td><td>距離矩陣(多變量 per-CpG 模式)</td><td>每 read 收成一個平均 β(純量)</td></tr>
<tr><td>測</td><td>per-CpG 模式 + 平均 + dispersion</td><td>只測整體平均差(有方向)</td></tr>
<tr><td>對 dispersion</td><td>會被影響</td><td>不受影響</td></tr></table>
<table><tr><th>PERMANOVA 分類</th><th>n</th><th>median \\|Δβ\\|</th></tr>
<tr><td>clean-location sig</td><td>{db['clean_n']}</td><td>{db['clean_med']}</td></tr>
<tr><td>dispersion-only sig</td><td>{db['disp_n']}</td><td>{db['disp_med']}</td></tr>
<tr><td>PERMANOVA 不顯著</td><td>{db['nonsig_n']}</td><td>{db['nonsig_med']}</td></tr></table>
<div class="key">三類 Δβ 幾乎一樣小(~0.045) → <b>PERMANOVA 的顯著不是整體平均差(Δβ)撐起來的</b>，它抓 per-CpG 模式+dispersion(Δβ 收成一個數看不到) → <b>兩做法測不同東西</b>。</div>

<h2>⑥ 4 組 HP-fine PERMANOVA 能處理嗎？✅ 能（實測）</h2>
<table><tr><th>標籤</th><th>testable</th><th>sig(p<.05)</th></tr>
<tr><td>2 組 HP-family</td><td>{fg['fam2']['testable']}</td><td>{fg['fam2']['sig']}</td></tr>
<tr><td><b>4 組 HP-fine</b></td><td><b>{fg['fine4']['testable']}</b></td><td><b>{fg['fine4']['sig']}</b></td></tr></table>
<p class="muted">4 組 testable 中實際有效組數(≥3 reads)：{fg['fine4']['ngroups']} → PERMANOVA omnibus 吃任意 K 能跑；但 somatic 子標籤 1-1/2-1 稀少 → 多數退化 2-3 組。<b>瓶頸是資料稀疏、非 PERMANOVA</b>。現行 C++ 合併 2 組為穩定 + Fisher-FH 處理 4 組細分。</p>

<h2>⑦ 真實熱圖對照（情況①離散 vs ②漸層）</h2>
<div class="hmgrid">
<figure class="hm"><img src="data:image/png;base64,{hm_grad}"/><figcaption><b>②漸層</b> chr12:28414696 \\|Δβ\\|=0.55 PERMANOVA p=.001 但 best_k=None(切不出)：甲基有 HP mean-shift、距離無離散塊</figcaption></figure>
<figure class="hm"><img src="data:image/png;base64,{hm_disc}"/><figcaption><b>①離散</b> chr8:131074695 align k=3：距離 3 對角塊=離散群，Fisher+V 與 PERMANOVA 都抓到</figcaption></figure></div>

<h2>⑧ 整合結論（誠實校正）</h2>
<div class="warn"><b>「切不出 ≠ 沒訊號」仍成立</b>(neither 只 {vr['neither_pct']}%)，但「訊號」性質要講清楚：<br>
• 乾淨「平均差/可分群」訊號其實只 ~20-29%(Fisher+V {fv['pct']}% / clean-location {vc['C_pct']}%，且兩者位點不同)。<br>
• 「85% 有訊號」**~72% 是 dispersion + per-CpG 模式**，非整體平均差、非離散群。<br>
• 所以切不出位點多數有「per-CpG 模式 / 散開度」層級差異，但只 ~20-29% 有乾淨平均差或可分群結構。</div>

<div class="foot">build_commit {BC} · 數字注入自 method_comparison.json + permanova_clean_4group.json(§13-A) · 源碼 SignificanceAnalyzer.cpp:288-304 / StructureTest.cpp run_permanova+check_dispersion · 熱圖 SoT dual-panel · 單樣本 ⭐2 · pilot(4組=chr21+chr22)</div>
</div></body></html>"""
out=f"{A}/20260620_method_comparison_fisher_v_vs_permanova_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB)")
