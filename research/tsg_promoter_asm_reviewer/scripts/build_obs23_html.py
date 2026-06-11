#!/usr/bin/env python3
"""Build obs23 ASM-mixture HTML with the 4-panel figure base64-inlined + regimes + hypotheses."""
import base64
from pathlib import Path

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
FIG = WS / "figures/obs23_asm_mixture_4panel.png"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/05/20260531_ASM_mixture_decomposition.standalone.html")

b64 = base64.b64encode(FIG.read_bytes()).decode()

HTML = f"""<!DOCTYPE html><html lang="zh-TW"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ASM 甲基差異 = 多 regime 混合 (2026-05-31)</title>
<style>
:root{{--ac:#D97757;--acd:#B85A3F;--tx:#141413;--txs:#6B6B66;--bg:#FAF9F5;--sf:#FFFFFF;--bd:#E3DACC;--cb:#F4EFE6;
--crit:#C2410C;--warn:#A16207;--ok:#15803D;--info:#1E3A8A;--ff:-apple-system,BlinkMacSystemFont,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif;--mo:"JetBrains Mono","DejaVu Sans Mono",monospace;}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--tx);font-family:var(--ff);font-size:16px;line-height:1.6}}
.wrap{{max-width:1080px;margin:0 auto;padding:26px 18px}}
main{{background:var(--sf);border:1px solid var(--bd);border-radius:8px;box-shadow:0 2px 8px rgba(20,20,19,.08);padding:30px}}
.hero{{border-bottom:2px solid var(--ac);padding-bottom:14px;margin-bottom:18px}}
.pretitle{{font-size:12px;text-transform:uppercase;letter-spacing:1px;color:var(--acd);font-weight:700}}
h1{{font-size:26px;margin:6px 0}}.sub{{color:var(--txs);font-size:15px}}
.verdict{{color:#fff;padding:15px;margin-top:12px;border-radius:6px;line-height:1.55;background:linear-gradient(90deg,#1E3A8A,#3B82F6)}}
.verdict strong{{font-weight:700}}
h2{{font-size:21px;margin:30px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--bd)}}
h2 .n{{font-family:var(--mo);color:var(--acd);font-size:.62em;margin-right:7px}}
h3{{font-size:16px;margin:16px 0 5px}}
table{{border-collapse:collapse;width:100%;margin:11px 0;font-size:13.5px}}
th,td{{border-bottom:1px solid var(--bd);padding:7px 9px;text-align:left;vertical-align:top}}
th{{background:var(--cb);font-weight:700;border-bottom:2px solid var(--ac)}}
td.num,th.num{{text-align:right;font-family:var(--mo)}}
tr.clean td{{background:#F0FDF4}}tr.dirty td{{background:#FEF2F0}}
.fig{{margin:16px 0;text-align:center;background:var(--bg);border:1px solid var(--bd);border-radius:6px;padding:12px}}
.fig img{{max-width:100%;height:auto;border-radius:4px}}
.cap{{font-size:13px;color:var(--txs);text-align:left;margin-top:10px;line-height:1.55}}
.cap b{{color:var(--tx)}}
.callout{{padding:11px 15px;margin:12px 0;border-radius:6px;border-left:4px solid}}
.c-ok{{background:#F0FDF4;border-color:var(--ok);color:#166534}}.c-fp{{background:#FEF2F0;border-color:var(--crit);color:#7F1D1D}}
.c-info{{background:#EFF6FF;border-color:var(--info);color:#1E40AF}}.c-warn{{background:#FEF7E6;border-color:var(--warn);color:#78350F}}
.callout strong{{display:block;margin-bottom:4px}}
.hyp{{border:1px solid var(--bd);border-left:4px solid var(--ac);border-radius:6px;padding:10px 14px;margin:10px 0;background:var(--sf)}}
.hyp .hid{{font-family:var(--mo);font-weight:700;color:var(--acd)}}
.badge{{display:inline-block;padding:1px 8px;border-radius:10px;font-size:11px;font-weight:700}}
.b-test{{background:#FEF3C7;color:#78350F}}.b-conf{{background:#DCFCE7;color:#166534}}
code{{background:var(--cb);padding:1px 6px;border-radius:3px;font-family:var(--mo);font-size:.9em}}
.prov{{background:var(--cb);border-top:2px solid var(--ac);margin-top:34px;padding:16px;border-radius:6px;font-size:12px;color:var(--txs);line-height:1.6}}
@media print{{body{{background:#fff}}.wrap{{padding:0}}main{{box-shadow:none;border:none;padding:0}}}}
</style></head><body><div class="wrap"><main>
<header class="hero">
<div class="pretitle">ASM 混合分解 · 整體 + 單獨 + 差分 + 觀察圖 + 假設 · 2026-05-31</div>
<h1>甲基差異 ≠ 單一現象 — 是「多 regime 混合」</h1>
<p class="sub">問題：我們量到的 ASM Δβ 是不是混合了多種狀況？答：是。HCC1395 paired_full · HP-axis · 顯著 ASM (p&lt;0.05, n=5,142) · 單樣本</p>
<div class="verdict"><strong>核心：ASM 顯著位點分成「乾淨 regime」與「confound regime」，FP-rate 差 3-4×（嚴格 strong tail 差 5×）。</strong>
所以「甲基差異」是真生物 ASM、LOH-CNV confound、低-coverage regression artifact 三者的<strong>混合</strong> — 不能當單一現象解讀，必須先分 regime。</div>
</header>

<h2><span class="n">§1</span>觀察圖 — 4 panel（整體看混合結構）</h2>
<div class="fig"><img src="data:image/png;base64,{b64}" alt="ASM mixture 4-panel">
<div class="cap">
<b>A（左上）Δβ vs coverage</b>：明顯<b>漏斗形</b> — n_cpg 低（左）時 Δβ 散開到 ±0.6，n_cpg 高（右，綠線=100）時收斂到 ±0.2。低 coverage = 雜訊膨脹。<br>
<b>B（右上）Δβ vs germline baseline</b>：<b>regression 形</b> — baseline 極端（紅區 &lt;0.2 或 &gt;0.8）的位點 Δβ 被強制往反方向（高 baseline→hypo、低 baseline→hyper），這是 regression-to-mean 不是生物。<br>
<b>C（左下）TP-rate heatmap（coverage × baseline）</b>：左下（高 cov + 中 baseline）= 可信；右上（低 cov + 極端 baseline）= 雜訊（FP 較多）。<br>
<b>D（右下）regime 組成</b>：4 regime 各自 n + TP-rate — 條長短 = 各 regime 大小，顏色 = 可信度。
</div></div>

<h2><span class="n">§2</span>整體 vs 單獨 — regime 分解表</h2>
<table>
<thead><tr><th>regime</th><th class="num">n</th><th class="num">FP-rate</th><th class="num">med &#124;Δβ&#124;</th><th class="num">med n_cpg</th><th class="num">LOH</th><th>性質</th></tr></thead>
<tbody>
<tr><td><b>整體 OVERALL</b></td><td class="num">5,142</td><td class="num">8.7%</td><td class="num">0.074</td><td class="num">31</td><td class="num">49%</td><td>混合（下分 4 類）</td></tr>
<tr class="clean"><td>A 乾淨（hiCov+midBase+nonLOH）</td><td class="num">78</td><td class="num">3.8%</td><td class="num">0.073</td><td class="num">143</td><td class="num">0%</td><td>可信 ASM（BRCA2 屬此）</td></tr>
<tr class="clean"><td>D intermediate（nonLOH）</td><td class="num">1,221</td><td class="num">2.7%</td><td class="num">0.112</td><td class="num">28</td><td class="num">0%</td><td>可信（最低 FP）</td></tr>
<tr class="dirty"><td>C LOH-confound</td><td class="num">1,158</td><td class="num">11.1%</td><td class="num">0.113</td><td class="num">30</td><td class="num">100%</td><td>LOH/CNV 混淆</td></tr>
<tr class="dirty"><td>B regression（loCov+extremeBase）</td><td class="num">2,685</td><td class="num">10.5%</td><td class="num">0.049</td><td class="num">32</td><td class="num">51%</td><td>雜訊膨脹（占最多）</td></tr>
</tbody></table>
<div class="callout c-info"><strong>差分（整體 − 單獨）— 混合證實</strong>
乾淨 regime（A 3.8% / D 2.7%）vs confound regime（C 11.1% / B 10.5%）<b>FP-rate 差 3-4×</b>（嚴格 Bonf+&#124;Δβ&#124;≥0.1 的 strong tail 差 5×，見 obs22）。整體 FP-rate 8.7% 是這個混合的平均，<b>掩蓋了 regime 間的巨大差異</b> — 這就是「為什麼必須見樹也見林」的具體例子：只看整體會誤判 ASM 是單一弱訊號，分 regime 才看到它是「少量可信 + 大量 artifact」的混合。</div>
<div class="callout c-warn"><strong>誠實的量級校正</strong>
A 乾淨 regime 只有 <b>78 位點（1.5% of 顯著 ASM）</b> — <b>真正可信的 ASM 是稀有的</b>；占最多的是 B regression（2,685, 52%）。所以「甲基差異」整體上<b>以 artifact 為主、可信生物為輔</b>。</div>

<h2><span class="n">§3</span>可能的假設（grounded in data）</h2>
<div class="hyp"><span class="hid">H1</span> <span class="badge b-conf">已驗證</span> <b>ASM Δβ = 多 regime 混合</b>（真 ASM + LOH-confound + regression-artifact）。證據：regime 間 FP-rate 差 3-5×、coverage/baseline/LOH 3 維 Bonferroni 確認不同（obs22）。<b>意涵</b>：任何 ASM 結論前必先分 regime，不可用整體均值。</div>
<div class="hyp"><span class="hid">H2</span> <span class="badge b-test">待驗</span> <b>乾淨 regime（A，n=78）的殘餘 ASM 是否真生物 + 可跨樣本重現？</b> 若濾掉 artifact 後 A regime 在 7 樣本一致，則「可信 ASM」是真現象（雖稀有）。測法：A regime 跨樣本一致性 + gene-context 富集。</div>
<div class="hyp"><span class="hid">H3</span> <span class="badge b-test">待驗</span> <b>LOH regime（C）的「ASM」其實是 allelic copy-number/loss，不是甲基化</b>。LOH 區物理上單 haplotype，HP1-vs-HP1-1 比較混入 CNV。測法：CNV-stratified 控制（CN=1 vs CN=2 分層）+ 與 ASCAT/CN 對照。</div>
<div class="hyp"><span class="hid">H4</span> <span class="badge b-test">待驗</span> <b>hyper-ASM 極端（Δβ&gt;+0.5, germ低→som飽和）與 hypo-ASM 是不同機制</b>，或都只是 baseline regression。測法：固定 baseline bin 後看 hyper/hypo 是否仍有方向訊號（去除 regression 後的殘差）。</div>
<div class="hyp"><span class="hid">H5</span> <span class="badge b-test">待驗</span> <b>「可信 ASM 篩選準則」（n_cpg≥100 + baseline 非極端 + nonLOH）能否把 ASM 從 anti-discriminative 救成 neutral/weak-positive？</b> 測法：只在 A regime 算 TP/FP 判別力（AUC），看是否 &gt;0.5。</div>

<h2><span class="n">§4</span>對你研究目標的結論</h2>
<ul>
<li><b>能清楚分析整體 + 單獨 + 差分</b>：✅ 整體 FP-rate 8.7% 分解成 4 regime（2.7%-11.1%），差分 3-4×。</li>
<li><b>能繪圖驗證</b>：✅ 4-panel 圖直接看出漏斗（coverage）+ regression（baseline）+ heatmap（混合結構）。</li>
<li><b>甲基差異確實混合多種狀況</b>：✅ 真 ASM（稀有 78）+ LOH-confound（1158）+ regression（2685）三者混合 — <b>不是單一現象</b>。</li>
<li><b>下一步最有價值</b>：H2（乾淨 regime 跨樣本驗證，看真 ASM 是否站得住）+ H3（LOH regime 的 CNV 控制）。</li>
</ul>

<div class="prov"><strong style="color:var(--acd)">Provenance</strong><br>
2026-05-31 · Git HEAD 274f152 · HCC1395 paired_full 單樣本 · 親自重算（per feedback_existing_artifacts_must_verify）<br>
數據 <code>genome_survey_v2/asm_dualaxis_{{tp,fp}}.tsv</code> · 腳本 <code>scripts/23_asm_mixture_decomposition.py</code>（+ obs20/21/22）· 輸出 <code>obs23_regime_summary.tsv</code> · 圖 base64 inline portable<br>
Tier [F] 全 first-hand。regime 邊界（cov≥100 / extremity&gt;0.3 / LOH）為 heuristic，非最佳化分界（H2-H5 可再 refine）。單樣本，跨樣本 NOT verified。ledger 61。
</div>
</main></div></body></html>"""

OUT.write_text(HTML)
print(f"written: {OUT} ({OUT.stat().st_size//1024} KB, figure inlined)")
