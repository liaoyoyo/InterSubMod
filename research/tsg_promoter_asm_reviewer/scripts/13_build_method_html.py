#!/usr/bin/env python3
"""Build a focused method/workflow standalone HTML — WHY long-read enables ASM."""
import base64

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/in_progress/2026/05/20260528_longread_ASM_method_workflow.standalone.html"

def b64(n):
    return base64.b64encode(open(f"{BASE}/figures/{n}", "rb").read()).decode()

fig_method = b64("method_shortread_vs_longread.png")
fig_readlevel = b64("brca2_read_level_methylation.png")
fig_region = b64("brca2_region_annotated_delta.png")

HTML = f"""<!DOCTYPE html>
<html lang="zh-Hant"><head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>流程與方法 — Long-read Allele-Specific Methylation (為何現在能做、以前不行)</title>
<script src="https://cdn.tailwindcss.com"></script>
<style>
  body {{ font-family: -apple-system, "Noto Sans TC", "PingFang TC", Roboto, sans-serif; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 0.5em 0.7em; vertical-align: top; }}
  th {{ background: #f1f5f9; }}
  code {{ background: #f1f5f9; padding: 0.1em 0.4em; border-radius: 0.25em; font-size: 0.9em; }}
  .step-card {{ background: white; border: 1px solid #e5e7eb; border-radius: 0.5em; padding: 1em 1.3em; margin: 0.8em 0; box-shadow: 0 1px 2px rgba(0,0,0,0.05); }}
  .step-num {{ display:inline-block; width:2em; height:2em; line-height:2em; text-align:center; background:#2563eb; color:white; border-radius:50%; font-weight:700; margin-right:0.5em; }}
  .why-long {{ background:#dcfce7; border-left:4px solid #16a34a; padding:0.5em 0.8em; margin-top:0.5em; font-size:0.9em; }}
  .yes {{ color:#15803d; font-weight:700; }}
  .no {{ color:#b91c1c; font-weight:700; }}
  #toc {{ position:sticky; top:1em; }}
  #toc a {{ display:block; padding:0.2em 0; color:#374151; text-decoration:none; font-size:0.9em; }}
  #toc a:hover {{ color:#2563eb; }}
</style></head>
<body class="bg-slate-50">
<div class="max-w-7xl mx-auto p-6 grid grid-cols-12 gap-6">

<aside class="col-span-3">
  <nav id="toc" class="bg-white border border-slate-200 rounded-lg p-4">
    <h3 class="font-bold text-slate-700 mb-3 text-sm">目錄</h3>
    <a href="#tldr">0 · 一句話 TL;DR</a>
    <a href="#core">1 · 核心技術突破對比 ⭐</a>
    <a href="#why-short-fail">2 · 為何 short-read 做不到</a>
    <a href="#why-long-work">3 · 為何 long-read 能做到</a>
    <a href="#pipeline">4 · 完整 8-step pipeline ⭐</a>
    <a href="#evidence">5 · 實際結果視覺證據</a>
    <a href="#novelty">6 · Novelty 誠實定位</a>
    <a href="#refs">7 · 工具與文獻</a>
  </nav>
</aside>

<main class="col-span-9 space-y-6">

<header class="bg-gradient-to-br from-teal-600 to-emerald-700 text-white p-6 rounded-lg shadow">
  <h1 class="text-2xl font-bold mb-2">流程與方法 — Long-read Allele-Specific Methylation</h1>
  <p class="text-emerald-100 text-sm">為什麼現在能做到「分單倍型的甲基化分析」，以前 short-read 時代不行</p>
  <p class="text-emerald-100 text-sm mt-1">2026-05-28 · HCC1395 BRCA2/ZAR1L bidirectional promoter case</p>
</header>

<section id="tldr" class="bg-white border-l-4 border-emerald-500 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-2">0. 一句話 TL;DR</h2>
  <p class="text-base leading-relaxed">
    我們用 <strong>long-read (ONT) + 5mC 修飾偵測 + longphase-S somatic haplotagging</strong>，讓<strong>同一條 read 同時帶三樣東西</strong>：
    (1) somatic 變異「身分證」、(2) 沿路所有 CpG 甲基化、(3) 足夠長度可被 phase 成單倍型。
    這讓我們能分辨「腫瘤 allele vs 正常 allele <strong>各自</strong>的甲基化」 —
    這是 short-read（Illumina / bisulfite array）<span class="no">物理上做不到</span>的 resolution，因為短片段無法同時涵蓋身分證與足夠甲基化點。
  </p>
</section>

<section id="core" class="bg-white border border-slate-200 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-3">1. 核心技術突破對比（整個方法的靈魂）</h2>
  <img src="data:image/png;base64,{fig_method}" class="w-full rounded border" alt="short vs long read schematic">
  <p class="text-sm text-slate-600 mt-2">
    <strong>上半（short-read）</strong>：每條短片段要嘛只覆蓋 SNV、要嘛只覆蓋 CpG，無法連結 → 只能算兩 allele 平均。<br>
    <strong>下半（long-read）</strong>：每條長 read 同時帶 SNV(身分證) + 沿路 CpG + 足夠長 → 黑點(甲基)/白點(去甲基)清楚顯示 HP1-1(腫瘤,紅) 比 HP1(正常,藍) 少甲基。
  </p>
</section>

<section id="why-short-fail" class="bg-white border border-slate-200 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-3">2. 為何 short-read 以前做不到（物理限制，非分析技巧能補）</h2>
  <table class="text-sm">
    <thead><tr><th>限制</th><th>原因</th><th>後果</th></tr></thead>
    <tbody>
      <tr><td>讀長太短 (~150 bp)</td><td>一個片段無法同時涵蓋 somatic 變異 + 足夠遠的 CpG</td><td>看到甲基化卻不知來自哪條 allele</td></tr>
      <tr><td>無法 phase</td><td>跨變異與 CpG 的 linkage 斷掉</td><td>無法重建單倍型</td></tr>
      <tr><td>bisulfite 轉換破壞序列</td><td>C→T 轉換降低 mapping 與 SNP 辨識</td><td>allele 辨識更難</td></tr>
    </tbody>
  </table>
  <div class="bg-red-50 border-l-4 border-red-500 p-3 mt-3 text-sm">
    <strong>關鍵</strong>：這是 <span class="no">物理限制</span> — 文獻明説 "conventional techniques cannot resolve the allelic distribution of methylation"（Do &amp; Tycko 2020；epialleleR 方法論）。不是演算法不夠好，是片段本身太短。
  </div>
</section>

<section id="why-long-work" class="bg-white border border-slate-200 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-3">3. 為何 long-read 現在能做到</h2>
  <p class="text-sm mb-3">關鍵：<strong>一條 ONT long-read（幾 kb ~ 幾十 kb）同時攜帶三樣資訊</strong>：</p>
  <table class="text-sm">
    <thead><tr><th>read 攜帶</th><th>提供什麼</th><th>工具</th></tr></thead>
    <tbody>
      <tr><td>① 5mC MM/ML modification tags</td><td>每個 CpG 的甲基化機率（測甲基）</td><td>ONT basecaller (Dorado 5mCG)</td></tr>
      <tr><td>② somatic 變異鹼基</td><td>身分證 — 知道這條 read 屬於哪條 haplotype</td><td>read 序列本身</td></tr>
      <tr><td>③ 足夠長度</td><td>跨越變異 + 多個 CpG，可被 phase</td><td>ONT 長讀特性</td></tr>
    </tbody>
  </table>
  <div class="why-long">
    <strong>三者合一</strong> → longphase-S 用 somatic 變異把 read 標成 HP:Z:1/2 (germline) 或 HP:Z:1-1/2-1 (腫瘤重建單倍型) → 我們就能在同一位點比「腫瘤 allele vs 正常 allele 各自的甲基化」= <span class="yes">allele-specific methylation</span>。
  </div>
</section>

<section id="pipeline" class="bg-white border border-slate-200 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-3">4. 完整 8-step pipeline（我們實際跑的流程）</h2>

  <div class="step-card">
    <span class="step-num">1</span><strong>Input — longphase-S somatic-haplotagged BAM</strong>
    <p class="text-sm mt-1">HCC1395 tumor-normal paired ONT BAM，經 longphase-S 用 somatic VCF 標 HP:Z tag (5 類: 1/2/1-1/2-1/3)，且 read 帶 5mCG MM/ML tags。</p>
    <div class="why-long">為何需 long-read：HP:Z:1-1/2-1 (腫瘤重建單倍型) 需要 read 夠長橫跨 somatic 變異 + germline 變異才能 phase。</div>
  </div>

  <div class="step-card">
    <span class="step-num">2</span><strong>定義 target region</strong>
    <p class="text-sm mt-1">BRCA2/ZAR1L bidirectional promoter ±1 kb (chr13:32,314,128-32,316,128)。</p>
  </div>

  <div class="step-card">
    <span class="step-num">3</span><strong>提取 region reads + 解析 HP:Z tag</strong>
    <p class="text-sm mt-1"><code>samtools view</code> 取 region reads，按 HP:Z 分組：germline HP1/HP2 vs 腫瘤重建 HP1-1/HP2-1。本 case: HP1=63, HP1-1=47。</p>
    <div class="why-long">為何 short-read 不行：短片段無 HP tag（無法 phase），只能混在一起算平均。</div>
  </div>

  <div class="step-card">
    <span class="step-num">4</span><strong>解析 5mC → per-read per-CpG 甲基化</strong>
    <p class="text-sm mt-1"><code>pysam read.modified_bases</code> 解 MM/ML tags，閾值 ML≥200=甲基 / ML≤50=去甲基。每條 read 沿路每個 CpG 都有一個甲基化 call。</p>
    <div class="why-long">為何 long-read 必要：因為 read ① 帶甲基 + ② 帶身分證在同一條上，CpG 甲基化能直接歸給特定 haplotype。</div>
  </div>

  <div class="step-card">
    <span class="step-num">5</span><strong>按 HP group 聚合 → per-CpG per-haplotype β</strong>
    <p class="text-sm mt-1">每個 CpG 位點分別算 HP1 (germline) 與 HP1-1 (somatic) 的甲基化率 β。要求每組 ≥3 reads。得 197 個 paired CpG。</p>
  </div>

  <div class="step-card">
    <span class="step-num">6</span><strong>統計檢定 — paired Wilcoxon</strong>
    <p class="text-sm mt-1">對 197 個 paired CpG 做 Wilcoxon signed-rank（HP1 vs HP1-1）。結果 Δβ=−0.122, p=6.1×10⁻¹¹, Cohen's d=−0.39。</p>
  </div>

  <div class="step-card">
    <span class="step-num">7</span><strong>Negative control — 排除 pipeline 雜訊</strong>
    <p class="text-sm mt-1">5 個隨機非-TSG 非-LOH somatic SNV 跑同流程，max |Δ|=0.054。我們訊號 0.122 超過 2.3×。</p>
  </div>

  <div class="step-card">
    <span class="step-num">8</span><strong>Region annotation + 視覺化</strong>
    <p class="text-sm mt-1">將每個差異 CpG map 到 functional region (CpG island core / shore / 5'UTR)，+ IGV (ZAR1L_ASM session) + matplotlib read-level 圖。發現差異集中在 CpG island shore (Irizarry 2009 現象)。</p>
  </div>
</section>

<section id="evidence" class="bg-white border border-slate-200 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-3">5. 實際結果視覺證據</h2>
  <h3 class="font-semibold mb-2">5a · Read-level 甲基化 (此流程的直接產物)</h3>
  <img src="data:image/png;base64,{fig_readlevel}" class="w-full rounded border" alt="read level methylation">
  <p class="text-sm text-slate-600 mt-1">每行=1條 read。黑=甲基 白=去甲基。HP1-1 腫瘤 reads (紅 box) 白點密度明顯高於 HP1 germline (藍 box) — 這正是 step 1-5 的產出。</p>

  <h3 class="font-semibold mb-2 mt-4">5b · Region 分布 (step 8 產物)</h3>
  <img src="data:image/png;base64,{fig_region}" class="w-full rounded border" alt="region annotated">
  <p class="text-sm text-slate-600 mt-1">差異集中在 CpG island shore (上游)，island core 兩 allele 都 β≈0 無差異。</p>
</section>

<section id="novelty" class="bg-white border-l-4 border-amber-500 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-3">6. Novelty 誠實定位（不 overclaim）</h2>
  <table class="text-sm">
    <thead><tr><th>面向</th><th>誰的貢獻</th></tr></thead>
    <tbody>
      <tr><td>long-read ASM 技術本身</td><td><span class="no">不是我們發明</span> — O'Neill 2024 (Cell Genomics) 等已用於 cancer</td></tr>
      <tr><td>longphase-S somatic haplotagging</td><td>CCU Bioinformatics Lab 工具</td></tr>
      <tr><td>BRCA2/ZAR1L bidirectional promoter</td><td>Liu 2010 已 published</td></tr>
      <tr><td>CpG island shore methylation 概念</td><td>Irizarry 2009 已 established</td></tr>
      <tr><td class="font-bold">我們的新貢獻</td><td class="font-bold"><span class="yes">把這個成熟技術 apply 到此 published bidirectional promoter，解析到過去 bulk 方法看不到的 allele-specific hypomethylation pattern</span></td></tr>
    </tbody>
  </table>
  <div class="bg-amber-50 border-l-4 border-amber-500 p-3 mt-3 text-sm">
    <strong>誠實 caveat</strong>：我們有 methylation data，<strong>無配對 RNA-seq</strong>。「影響 BRCA2 表達」目前是基於機制 plausibility 的 hypothesis，需 allele-specific RNA-seq 驗證。novelty 在「應用 + 特定發現」，不在「發明技術」。
  </div>
</section>

<section id="refs" class="bg-white border border-slate-200 rounded-lg p-5 shadow">
  <h2 class="text-xl font-bold mb-3">7. 工具與文獻</h2>
  <h3 class="font-semibold mb-1 text-sm">工具</h3>
  <ul class="text-sm list-disc list-inside mb-3">
    <li><strong>longphase-S</strong> (CCU Bioinformatics Lab) — somatic haplotagging, 輸出 HP:Z tag</li>
    <li><strong>ONT Dorado</strong> — 5mCG basecalling (MM/ML tags)</li>
    <li><strong>pysam</strong> — MM/ML 解析 (read.modified_bases)</li>
    <li><strong>samtools / scipy</strong> — region 提取 + Wilcoxon</li>
    <li><strong>IGV</strong> (headless JDK21) — ZAR1L_ASM session 視覺化</li>
  </ul>
  <h3 class="font-semibold mb-1 text-sm">關鍵文獻</h3>
  <ul class="text-sm list-disc list-inside">
    <li>Do &amp; Tycko 2020 <em>Genome Biology</em> — ASM in cancer / conventional 無法 resolve allelic methylation. PMID 32594908</li>
    <li>O'Neill 2024 <em>Cell Genomics</em> — long-read phased ASM in cancer (methodology peer). PMID 39406235</li>
    <li>Liu (Misra) 2010 <em>Mol Cancer</em> — BRCA2/ZAR2 bidirectional promoter. PMID 20202217</li>
    <li>Irizarry 2009 <em>Nature Genetics</em> — CpG island shore methylation. PMID 19151715</li>
    <li>Evans 2018 <em>AJHG</em> — BRCA1 c.-107A&gt;T single-base → cis-ASM precedent. PMID 30075112</li>
  </ul>
</section>

<footer class="text-center text-xs text-slate-500 py-6">
  Method/workflow companion — 2026-05-28<br>
  主報告: <code>InterSubMod/docs/reports/in_progress/2026/05/20260527_HCC1395_TSG_promoter_ASM_reviewer_response_01.md</code>
</footer>

</main></div>
</body></html>
"""

open(OUT, "w").write(HTML)
print(f"Wrote {OUT}")
print(f"Size: {len(HTML):,} chars ({len(HTML)/1024:.1f} KB)")
