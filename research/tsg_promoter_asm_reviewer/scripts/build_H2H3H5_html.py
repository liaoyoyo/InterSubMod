#!/usr/bin/env python3
"""Build interactive H2/H3/H5 results HTML with figures base64-inlined + collapsible method audits."""
import base64
from pathlib import Path

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
FIGD = WS / "figures"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/05/20260531_ASM_H2H3H5_interactive.standalone.html")

def b64(name):
    p = FIGD / name
    return base64.b64encode(p.read_bytes()).decode() if p.exists() else ""

figs = {n: b64(n) for n in [
    "obs24_H2_within_sample_consistency.png",
    "obs25_H3_cn_stratification.png", "obs25_H3_coverage_control.png",
    "obs26_H5_roc.png", "obs26_H5_perm_null.png", "obs26_H5_confound.png", "obs26_H5_regime_composition.png",
]}

def img(name, cap):
    return f'<figure class="fig"><img src="data:image/png;base64,{figs[name]}" alt="{name}"><figcaption>{cap}</figcaption></figure>'

HTML = f"""<!DOCTYPE html><html lang="zh-TW"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ASM H2/H3/H5 深入驗證 — 互動結果 (2026-05-31)</title>
<style>
:root{{--ac:#D97757;--acd:#B85A3F;--tx:#141413;--txs:#6B6B66;--bg:#FAF9F5;--sf:#FFFFFF;--bd:#E3DACC;--cb:#F4EFE6;
--crit:#C2410C;--warn:#A16207;--ok:#15803D;--info:#1E3A8A;--ff:-apple-system,BlinkMacSystemFont,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif;--mo:"JetBrains Mono","DejaVu Sans Mono",monospace;}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--tx);font-family:var(--ff);font-size:16px;line-height:1.6}}
.layout{{max-width:1180px;margin:0 auto;display:grid;grid-template-columns:200px 1fr;gap:22px;padding:24px 16px}}
aside{{position:sticky;top:16px;align-self:start;background:var(--sf);border:1px solid var(--bd);border-radius:8px;padding:13px 11px;max-height:calc(100vh - 48px);overflow-y:auto}}
aside h3{{font-size:11px;text-transform:uppercase;letter-spacing:.5px;color:var(--txs);margin:0 0 8px}}
aside a{{display:block;padding:5px 8px;border-radius:4px;color:var(--tx);text-decoration:none;font-size:13px;border-left:3px solid transparent}}
aside a:hover{{background:var(--cb);border-left-color:var(--ac)}}aside a.active{{background:var(--cb);border-left-color:var(--acd);font-weight:600}}
main{{background:var(--sf);border:1px solid var(--bd);border-radius:8px;box-shadow:0 2px 8px rgba(20,20,19,.08);padding:30px;min-width:0}}
.hero{{border-bottom:2px solid var(--ac);padding-bottom:14px;margin-bottom:18px}}
.pretitle{{font-size:12px;text-transform:uppercase;letter-spacing:1px;color:var(--acd);font-weight:700}}
h1{{font-size:25px;margin:6px 0}}.sub{{color:var(--txs);font-size:15px}}
.verdict{{color:#fff;padding:15px;margin-top:12px;border-radius:6px;line-height:1.55;background:linear-gradient(90deg,#1E3A8A,#3B82F6)}}
.verdict strong{{font-weight:700}}
h2{{font-size:21px;margin:32px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--bd);scroll-margin-top:16px}}
h2 .n{{font-family:var(--mo);color:var(--acd);font-size:.6em;margin-right:7px}}
h3{{font-size:16px;margin:16px 0 6px}}
table{{border-collapse:collapse;width:100%;margin:11px 0;font-size:13px}}
th,td{{border-bottom:1px solid var(--bd);padding:6px 9px;text-align:left;vertical-align:top}}
th{{background:var(--cb);font-weight:700;border-bottom:2px solid var(--ac)}}
td.num,th.num{{text-align:right;font-family:var(--mo)}}
.fig{{margin:14px 0;text-align:center;background:var(--bg);border:1px solid var(--bd);border-radius:6px;padding:11px}}
.fig img{{max-width:100%;height:auto;border-radius:4px}}
.fig figcaption{{font-size:12.5px;color:var(--txs);text-align:left;margin-top:8px;line-height:1.5}}
.fig-grid{{display:grid;grid-template-columns:1fr 1fr;gap:12px}}
@media(max-width:760px){{.fig-grid{{grid-template-columns:1fr}}.layout{{grid-template-columns:1fr}}aside{{position:static;max-height:none}}}}
.callout{{padding:11px 15px;margin:12px 0;border-radius:6px;border-left:4px solid}}
.c-ok{{background:#F0FDF4;border-color:var(--ok);color:#166534}}.c-fp{{background:#FEF2F0;border-color:var(--crit);color:#7F1D1D}}
.c-info{{background:#EFF6FF;border-color:var(--info);color:#1E40AF}}.c-warn{{background:#FEF7E6;border-color:var(--warn);color:#78350F}}
.callout strong{{display:block;margin-bottom:4px}}
.vb{{display:inline-block;padding:2px 10px;border-radius:12px;font-size:12px;font-weight:700}}
.v-conf{{background:#DCFCE7;color:#166534}}.v-ref{{background:#DBEAFE;color:#1E3A8A}}.v-inc{{background:#FEF3C7;color:#78350F}}.v-block{{background:#FEE2E2;color:#991B1B}}
details{{border:1px solid var(--bd);border-radius:6px;margin:10px 0;background:var(--sf)}}
details>summary{{cursor:pointer;padding:9px 14px;background:var(--cb);font-weight:600;font-size:14px;list-style:none}}
details>summary::-webkit-details-marker{{display:none}}
details>summary::before{{content:"▸ ";color:var(--txs)}}details[open]>summary::before{{content:"▾ "}}
details .body{{padding:10px 16px;font-size:13.5px}}
code{{background:var(--cb);padding:1px 6px;border-radius:3px;font-family:var(--mo);font-size:.9em}}
.prov{{background:var(--cb);border-top:2px solid var(--ac);margin-top:34px;padding:16px;border-radius:6px;font-size:12px;color:var(--txs);line-height:1.6}}
@media print{{body{{background:#fff}}.layout{{display:block}}aside{{display:none}}main{{box-shadow:none;border:none;padding:0}}details{{border:none}}details>.body{{display:block!important}}}}
</style></head><body><div class="layout">
<aside><h3>章節</h3>
<a href="#syn">統一綜合 ★</a><a href="#h3">H3 — LOH=CNV? (positive)</a><a href="#h5">H5 — 救判別?</a><a href="#h2">H2 — 跨樣本一致?</a><a href="#prov">Provenance</a>
</aside>
<main>
<header class="hero">
<div class="pretitle">ASM H2/H3/H5 深入驗證 · workflow scout+execute · 嚴格 method audit · 2026-05-31</div>
<h1>真 ASM 到底站不站得住 — 三假設驗證結果</h1>
<p class="sub">HCC1395 paired_full 單樣本 · 每 H 經 scout 驗資料/方法 + execute 完整統計 + permutation/confound 控制。每數字親自重算（ledger 63）。</p>
<div class="verdict"><strong>核心精細化：ASM「真實，但非方向、非判別、非 CNV」。</strong>
連續 ASM 訊號的判別力 AUC=0.505 <strong>本來就中性</strong>（不是 anti-discriminative，那只是 strong 子集性質）；LOH ASM <strong>不是 CNV artifact</strong>（真甲基）；但連 credible regime 也<strong>無可重現方向</strong>。唯一 positive thread = copy-neutral LOH 的真 ASM（單樣本待跨驗）。</div>
</header>

<h2 id="syn"><span class="n">§0</span>統一綜合 — 三 H 合起來看</h2>
<table>
<thead><tr><th>H</th><th>問</th><th>Verdict</th><th>關鍵數字</th></tr></thead>
<tbody>
<tr><td><b>H3</b></td><td>LOH ASM 是 CNV?</td><td><span class="vb v-ref">REFUTED CNV</span></td><td>cnLOH(CN=2)|Δβ|=0.082 &gt; gainLOH(CN&gt;2)=0.070（反 dosage）；coverage-controlled van Elteren p=7.5e-7… 7.5e-5；TP-specific r=−0.153</td></tr>
<tr><td><b>H5</b></td><td>篩選救判別?</td><td><span class="vb v-inc">INCONCLUSIVE / neutral</span></td><td>連續|Δβ| AUC=0.505（已中性，perm-p=0.50）；credible AUC 0.57 但 perm-p=0.51（FP=8 無 power）；coverage AUC 0.54 &gt; ASM</td></tr>
<tr><td><b>H2</b></td><td>跨樣本一致?</td><td><span class="vb v-block">BLOCKED → proxy</span></td><td>真測法卡（無資料+磁碟+長計算）；proxy: credible 無淨方向（100hyper/98hypo, p=0.94）</td></tr>
</tbody></table>
<div class="callout c-info"><strong>合起來的故事（精細化 anti-discriminative）</strong>
過去説「ASM anti-discriminative」其實有兩層需要分清：(1) <b>連續 |Δβ| 訊號本來就中性</b>（AUC 0.505，H5）— 不是「反判別」；(2) 真正 <b>FP-enriched 5× 的是 strong-ASM 嚴格子集</b>（低 cov + 極端 baseline 的雜訊尾）。而 <b>coverage（n_cpg）才是真正（弱）的 TP/FP 區分軸</b>（AUC 0.54 &gt; ASM 0.505）。所以 ASM 本身既不能幫也不太會害，是「中性 + 被 coverage 調制」的量。</div>
<div class="callout c-ok"><strong>唯一 positive thread（H3）</strong>
LOH 區的 ASM <b>不是 CNV 假象</b> — copy-neutral cnLOH（CN=2，沒有劑量變化）反而有最大 ASM，且抗 coverage 控制、TP-specific。這指向 <b>copy-neutral LOH 區存在真的 allele-specific 甲基化</b>。這是目前最值得追的生物訊號（單樣本，需跨樣本驗證）。</div>

<h2 id="h3"><span class="n">§1</span>H3 — LOH ASM 是 CNV 嗎？ <span class="vb v-ref">REFUTED CNV（真甲基-leaning）</span></h2>
<p>把 LOH 區顯著 ASM 位點按 CN 類別分層（SEQC2 NGS-truth CNV BED），測 Δβ 是否隨 copy-number 系統性變化（若是 CNV artifact → 劑量越大 Δβ 越大）。</p>
<table>
<thead><tr><th>CN 類別</th><th class="num">n</th><th class="num">median |Δβ|</th><th class="num">TP-rate</th></tr></thead>
<tbody>
<tr><td><b>cnLOH（CN=2 copy-neutral）</b></td><td class="num">913</td><td class="num"><b>0.082</b></td><td class="num">0.88</td></tr>
<tr><td>gainLOH（CN&gt;2 擴增）</td><td class="num">1,580</td><td class="num">0.070</td><td class="num">0.85</td></tr>
<tr><td>lossLOH（CN&lt;2，underpowered）</td><td class="num">22</td><td class="num">0.085</td><td class="num">0.50</td></tr>
</tbody></table>
<div class="callout c-ok"><strong>為何 REFUTED CNV</strong>
若 ASM 是 CNV dosage artifact，擴增（CN&gt;2）應有最大 Δβ。但實際 <b>copy-neutral cnLOH（CN=2）最大</b>（0.082 &gt; 0.070）= <b>與 dosage 方向相反</b>。且效應抗 coverage 控制（van Elteren Z=3.96 p=7.5e-5；partial Spearman p=1.6e-3，8 test Bonferroni 全過），TP-specific（FP 翻號，TP-only r=−0.153 p=3e-9）。→ 不是 CNV，是 copy-neutral 區的真 allele-specific 甲基化。</div>
<div class="fig-grid">
{img("obs25_H3_cn_stratification.png","H3-A: CN 三類分層的 Δβ/|Δβ| 分布 — cnLOH(CN=2) 反而最大 = 反 dosage")}
{img("obs25_H3_coverage_control.png","H3-B: coverage 控制 — CN 效應在固定 coverage 下仍存活（van Elteren / partial Spearman / per-stratum 同方向）")}
</div>
<details><summary>方法學自審（用戶要求）+ 完整統計 + caveat</summary><div class="body">
<b>觀察依據</b> ✓ — gain/loh/loss BED 區真實 overlap（cnLOH 913Mb / gain∩loh 520Mb / loss∩loh 57Mb），CN 分層結構可定義。<br>
<b>分析支持</b> ✓ — coverage 是主 confound（scout 正確標），naive per-stratum MWU underpowered，但 coverage-controlled pooled（van Elteren）穩健；用全 2,515 LOH-sig（非 order-dependent 的 regime-C n=1158 子集，避免人為排除）。<br>
<b>完整統計</b>：KW 3-class p=5.4e-5；MWU cnLOH-vs-gain p=3.6e-5 r=−0.099；Spearman |Δβ| vs CN-ord rho=−0.087 p=1.2e-5；coverage-controlled 全過 Bonferroni×8（α=0.0063）。<br>
<b>caveat</b>：(1) BED 只有類別非整數 CN，「CN&gt;2/=2/&lt;2」是推論，不能測 dose-response；(2) lossLOH n=22 不可測，實為 2-bin（cnLOH vs gainLOH）；(3) |Δβ| 主要由 coverage 驅動（0.13@cov&lt;20→0.04@cov≥80），CN 差是二階調制；(4) <b>單樣本</b>；(5) 「REFUTED CNV」只排除一個 artifact 解釋，<b>不正面證明是 TSG-relevant 生物 ASM</b>。
</div></details>

<h2 id="h5"><span class="n">§2</span>H5 — credible 篩選能救判別力嗎？ <span class="vb v-inc">INCONCLUSIVE（仍中性）</span></h2>
<p>算 |Δβ| 預測該位點是 TP 還是 FP 的 AUC（per-position，permutation null + coverage-bin 控制）：全部 vs 只 credible regime。</p>
<table>
<thead><tr><th>set</th><th class="num">n (TP/FP)</th><th class="num">AUC</th><th class="num">perm-p</th><th>判讀</th></tr></thead>
<tbody>
<tr><td>全 HP-axis ASM</td><td class="num">23,444 (21633/1811)</td><td class="num">0.505</td><td class="num">0.496</td><td><b>已中性</b>（非 anti）</td></tr>
<tr><td>顯著 only</td><td class="num">5,056</td><td class="num">0.502</td><td class="num">0.899</td><td>中性</td></tr>
<tr><td>credible regime</td><td class="num">189 (181/8)</td><td class="num">0.570</td><td class="num">0.509</td><td>FP=8 無 power（null95 [0.29,0.71]）</td></tr>
<tr><td><b>coverage(n_cpg)</b></td><td class="num">—</td><td class="num"><b>0.540</b></td><td class="num">—</td><td><b>比 ASM 更會分</b></td></tr>
</tbody></table>
<div class="callout c-warn"><strong>重大 reframe</strong>
H5 問題本身需修正：連續 |Δβ| 的全域 AUC <b>本來就是 0.505（中性）</b>，不是「&lt;0.5 反判別」。「anti-discriminative 5×」是 <b>strong-ASM 嚴格子集</b>（FP-enriched 尾）的性質，不是連續分數。credible filter 把 AUC 推到 0.57 但 perm-p=0.51（FP 只剩 8 個，統計上與隨機無法區分）。<b>真正（弱）的 TP/FP 軸是 coverage（AUC 0.54）不是 ASM。</b></div>
<div class="fig-grid">
{img("obs26_H5_roc.png","H5-A: ROC — 全部/顯著/credible ASM 的 |Δβ| 預測 is_TP，全部貼近對角線（中性）")}
{img("obs26_H5_perm_null.png","H5-B: permutation null — credible AUC 0.57 落在 null95 [0.29,0.71] 內，與隨機無法區分（FP=8）")}
</div>
<details><summary>方法學自審 + confound + caveat</summary><div class="body">
<b>per-position 聚合</b>（scout 強制修正）：23,840 row → 23,444 unique position（max|Δβ|），避免雙計違反 i.i.d.；TP/FP 互斥驗證（1/23,444 衝突可忽略）。<br>
<b>confound</b>：(1) FP=8 → AUC 有效樣本由少數類決定，null95 極寬，無法區分 0.57 vs 隨機；(2) <b>無 caller_af 欄</b> → auc-confound-guard 的 AF-bin cross + within-group OLS 跑不了（caller_af 是已知最強 TP/FP proxy）；(3) coverage 是真弱判別軸（AUC 0.54 &gt; ASM 0.505），within-bin AUC 0.49-0.61 同方向但 FP 稀少。<br>
<b>scope</b>：單樣本 → 即使 weak-positive 成立也 cap ⭐2；Cynefin Complex（多 regime 混合）→ probe-first 無 deterministic 結論。<br>
<b>結論</b>：credible filter <b>未能 demonstrably 救</b> ASM 成判別性；但也確認 ASM 連續訊號本來就中性（非害）。
</div></details>

<h2 id="h2"><span class="n">§3</span>H2 — credible ASM 跨樣本一致嗎？ <span class="vb v-block">真測法 BLOCKED → within-sample proxy</span></h2>
<div class="callout c-fp"><strong>真跨樣本 BLOCKED（誠實）</strong>
6 樣本<b>沒有預先算好的 MSA 甲基 ASM per-site 數據</b>（只有 phasing significance_summary）。重跑 = 三重阻擋：(1) big7 97% / big8 98% 磁碟滿（pipeline_block_check Hard Gate）；(2) ~20-30hr（6 MSA rerun）；(3) 用戶本輪未授權長計算（模型自判長計算 = 🔴）。→ 改做 within-HCC1395 proxy（染色體對半切 + 5000× bootstrap），<b>明確標 partial，不可宣稱跨獨立樣本重現</b>。</div>
<p><b>proxy 結果</b>：credible regime（n=198, 75 sig）染色體兩半 Δβ 分布<b>一致</b>（KS p=0.28 / MWU|Δβ| p=0.24 / median signed Δβ CI overlap）→ 信號內部穩定。<b>但無淨方向</b>：mean Δβ CI [−0.0075, +0.0126] 跨 0、100 hyper/98 hypo、binomial p=0.94。即 credible ASM 是<b>對稱幅度雜訊</b>，與 anti-discriminative/中性一致。FP-arm credible-FP 僅 n=7（證實 FP power 崩潰）。</p>
{img("obs24_H2_within_sample_consistency.png","H2: within-HCC1395 proxy — 染色體對半切一致性 + bootstrap。信號內部穩定但無方向性")}
<div class="callout c-warn"><strong>意涵</strong>
proxy 通過只是「跨樣本可重現」的<b>必要非充分</b>條件。更重要：credible regime <b>無淨方向</b>（對稱 hyper/hypo）→ 即使未來補跑跨樣本，最可能確認的是「<b>無可重現方向性信號</b>」而非有信號。要正式做 H2 需先 /infra-ops 清磁碟 + 用戶授權 + 6 MSA rerun；且須標 HCC1395=SEQC2 gold-LOH anchor vs 其他=pipeline-LOH 不同質（獨立樣本實為 5）。</div>

<div class="prov" id="prov"><strong style="color:var(--acd)">Provenance</strong><br>
2026-05-31 · Git HEAD 274f152 · HCC1395 paired_full 單樣本 · workflow wdvn0uhnl (6 agents, scout+execute, method-audit) · 親自重算<br>
數據 <code>InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/obs24-26_*.json</code> · 腳本 <code>scripts/{{26,30,31}}_*.py</code> · 圖 7 張 base64 inline · ledger 63<br>
Tier：H3 REFUTED-CNV 較穩（coverage-controlled, Bonferroni×8）但單樣本 ⭐3；H5 中性 ⭐2；H2 真測法未做（proxy partial）。所有 verdict 標單樣本。CN 為類別非整數（dose-response 不可測）。H2 跨樣本 = infra + 授權 + 6 MSA rerun 後可做。
</div>
</main></div>
<script>
(function(){{var L=document.querySelectorAll('aside a'),S=Array.from(L).map(a=>document.getElementById(a.getAttribute('href').slice(1))).filter(Boolean);
function u(){{var y=scrollY+90,c=S[0];S.forEach(s=>{{if(s.offsetTop<=y)c=s}});L.forEach(a=>a.classList.remove('active'));if(c){{var e=document.querySelector('aside a[href="#'+c.id+'"]');if(e)e.classList.add('active')}}}}
var t=false;addEventListener('scroll',()=>{{if(!t){{requestAnimationFrame(()=>{{u();t=false}});t=true}}}});u();}})();
</script>
</body></html>"""

OUT.write_text(HTML)
print(f"written: {OUT} ({OUT.stat().st_size//1024} KB, 7 figures inlined)")
