#!/usr/bin/env python3
"""Build standalone HTML companion report with embedded figures."""
import base64, json, os

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/in_progress/2026/05/20260527_HCC1395_TSG_promoter_ASM_reviewer_response_01.standalone.html"

def b64(name):
    with open(f"{BASE}/figures/{name}", "rb") as f:
        return base64.b64encode(f.read()).decode()

fig1 = b64("brca2_per_cpg_delta.png")
fig2 = b64("brca2_beta_distribution.png")
fig3 = b64("brca2_vs_controls.png")
fig4 = b64("brca2_read_level_methylation.png")
fig5_igv_full = b64("brca2_igv_zar1l_session_full.png")
fig6_igv_zoom = b64("brca2_igv_zar1l_session_zoom.png")
fig7_region = b64("brca2_region_annotated_delta.png")

with open(f"{BASE}/output/prereq_a_summary.json") as f:
    summary = json.load(f)
with open(f"{BASE}/output/step3_hp_coverage.json") as f:
    step3 = json.load(f)
with open(f"{BASE}/output/step4_ism_results.json") as f:
    step4 = json.load(f)
with open(f"{BASE}/output/step4_negative_control.json") as f:
    ctrl = json.load(f)

# Lost-to-LOH TSG groups
classic_breast = ['BRCA1', 'RASSF1', 'CDH13', 'CHEK1', 'ATM', 'ATR', 'RAD51', 'RAD51C', 'RAD51D', 'PALB2', 'BRIP1', 'MLH1', 'MLH3', 'MUTYH', 'FANCD2', 'FANCF', 'FANCG', 'FANCM', 'NBN', 'MRE11', 'RAD50']
cell_cycle = ['TP53', 'WT1', 'MEN1', 'VHL', 'APC', 'NF1']
chromatin = ['ARID1A', 'ARID1B', 'ARID2', 'ASXL2', 'EP300', 'GATA3', 'KEAP1', 'KMT2D', 'PBRM1', 'SETD2', 'SMARCA4', 'SMARCE1', 'ESR1', 'PGR', 'RARB', 'TBX3', 'TWIST1']
other = ['APAF1', 'AXIN2', 'BAP1', 'CTNNA1', 'FAT1', 'FBXW7', 'STK11']

HTML = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="UTF-8">
<title>HCC1395 TSG Promoter ASM — Reviewer Answer (2026-05-27)</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<script src="https://cdn.tailwindcss.com"></script>
<style>
  body {{ font-family: -apple-system, "Noto Sans TC", "PingFang TC", Roboto, sans-serif; }}
  details > summary {{ cursor: pointer; padding: 0.5em; background: #f8fafc; border-left: 4px solid #2563eb; }}
  details[open] > summary {{ background: #e0f2fe; }}
  .tier-strict {{ background: #16a34a; color: white; padding: 0.15em 0.5em; border-radius: 0.3em; font-size: 0.85em; font-weight: 600; }}
  .tier-relaxed {{ background: #ca8a04; color: white; padding: 0.15em 0.5em; border-radius: 0.3em; font-size: 0.85em; font-weight: 600; }}
  .tier-fail {{ background: #6b7280; color: white; padding: 0.15em 0.5em; border-radius: 0.3em; font-size: 0.85em; font-weight: 600; }}
  .stat-good {{ color: #16a34a; font-weight: 600; }}
  .stat-bad {{ color: #dc2626; font-weight: 600; }}
  .stat-neutral {{ color: #6b7280; font-weight: 600; }}
  code.pos {{ background: #f1f5f9; padding: 0.1em 0.4em; border-radius: 0.25em; font-size: 0.9em; }}
  table {{ border-collapse: collapse; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 0.4em 0.7em; }}
  th {{ background: #f1f5f9; }}
  #toc {{ position: sticky; top: 1em; max-height: calc(100vh - 2em); overflow-y: auto; }}
  #toc a {{ display: block; padding: 0.2em 0; color: #374151; text-decoration: none; font-size: 0.9em; }}
  #toc a:hover {{ color: #2563eb; }}
  #toc a.active {{ color: #dc2626; font-weight: 600; border-left: 3px solid #dc2626; padding-left: 0.5em; }}
  .anchor-target {{ scroll-margin-top: 1em; }}
  blockquote {{ background: #fef3c7; padding: 1em; border-left: 4px solid #f59e0b; margin: 1em 0; }}
  .container-card {{ background: white; border: 1px solid #e5e7eb; border-radius: 0.5em; padding: 1.2em; margin: 1em 0; box-shadow: 0 1px 2px rgba(0,0,0,0.05); }}
  .verdict-positive {{ background: linear-gradient(90deg, #16a34a 0%, #22c55e 100%); color: white; padding: 1em 1.5em; border-radius: 0.5em; font-weight: 600; font-size: 1.1em; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
  .verdict-caveat {{ background: #f59e0b; color: white; padding: 0.7em 1em; border-radius: 0.3em; margin: 0.5em 0; }}
</style>
</head>
<body class="bg-slate-50">
<div class="max-w-7xl mx-auto p-6 grid grid-cols-12 gap-6">

<!-- TOC sidebar -->
<aside class="col-span-3" id="toc-wrapper">
  <nav id="toc" class="bg-white border border-slate-200 rounded-lg p-4">
    <h3 class="font-bold text-slate-700 mb-3 text-sm">目錄</h3>
    <a href="#tldr">0 · TL;DR</a>
    <a href="#workflow">0.5 · 完整研究流程 ⭐</a>
    <a href="#reviewer-q">1 · Reviewer Q 拆解</a>
    <a href="#pipeline">2 · Pipeline 5-step</a>
    <a href="#inventory">3 · Data Sources</a>
    <a href="#yield">4 · Yield Funnel 95→1</a>
    <a href="#lost-loh">　4.1 · 51 lost-to-LOH TSG</a>
    <a href="#hp-cov">　4.2 · 6→2 HP coverage</a>
    <a href="#brca2">5 · BRCA2 Deep Dive ⭐</a>
    <a href="#brca2-stats">　5.2 · Statistics</a>
    <a href="#brca2-fig1">　5.4 · Figure 1 per-CpG</a>
    <a href="#brca2-fig2">　5.4 · Figure 2 violin</a>
    <a href="#brca2-fig3">　5.4 · Figure 3 vs control</a>
    <a href="#brca2-fig4">　5.4 · Figure 4 read-level</a>
    <a href="#brca2-region">　5.5 · CpG shore region 分布 ⭐</a>
    <a href="#brca2-igv">　5.7 · IGV ZAR1L Session ⭐</a>
    <a href="#brca2-igv">　5.8 · session XML 診斷</a>
    <a href="#negctrl">6 · Negative Control</a>
    <a href="#caveats">7 · Caveats × 6</a>
    <a href="#tradeoff">8 · Methodological Trade-off</a>
    <a href="#response">9 · Response-to-Reviewer ⭐</a>
    <a href="#literature">9b · Literature Evidence Base ⭐</a>
    <a href="#repro">10 · Reproducibility</a>
    <a href="#extensions">11 · Extensions</a>
    <a href="#glossary">12 · Glossary</a>
    <a href="#refs">13 · References</a>
  </nav>
</aside>

<!-- Main content -->
<main class="col-span-9 space-y-6">

<header class="bg-gradient-to-br from-blue-600 to-indigo-700 text-white p-6 rounded-lg shadow">
  <h1 class="text-2xl font-bold mb-2">HCC1395 TSG Promoter — Reconstructed Somatic Haplotype Hypomethylation</h1>
  <p class="text-blue-100 text-sm mb-1">Reviewer answer pipeline · 2026-05-27 · Task type B validation (whole-genome, 95 curated TSG, non-LOH)</p>
  <p class="text-blue-100 text-sm">build_date 2026-05-27 / status in_progress / verdict POSITIVE (1 strong example)</p>
</header>

<section id="tldr" class="anchor-target">
  <h2 class="text-xl font-bold mb-3">0. TL;DR — Reviewer 答案</h2>
  <div class="verdict-positive mb-3">
    ✅ <strong>Yes, we have one.</strong> BRCA2 promoter (chr13:32,315,128 G>A, TVAF=0.189) 的 reconstructed somatic haplotype (HP:Z:1-1) 在 ±1kb window 內 <strong>197 paired CpG</strong>:
    <ul class="mt-2 list-disc list-inside text-base">
      <li>mean Δβ = <strong>−0.122</strong> (somatic hypomethylated vs germline) ✓ exact reviewer direction</li>
      <li>Wilcoxon signed-rank <strong>p = 6.1 × 10⁻¹¹</strong></li>
      <li>Cohen's d = −0.39</li>
      <li>max single-CpG |Δβ| = 0.93</li>
      <li>vs 5 random non-TSG controls (max |Δ|=0.054): BRCA2 <strong>exceeds null distribution 2.3 ×</strong></li>
      <li>變異點距 BRCA2 promoter CpG island 268 bp，落 canonical 啟動子區內</li>
      <li>IGV 易 visualize (47 HP:Z:1-1 reads 全帶 5mCG MM/ML，xml session 已備)</li>
    </ul>
  </div>
  <div class="verdict-caveat">
    ⚠ <strong>核心 caveat</strong>: 95 個 curated TSG 中**唯一**通過 strict 條件的 strong example。51/95 (53.7%) TSG 在 LOH 內（含 BRCA1/TP53/CDH1/RASSF1A 等 reviewer 可能預期的 classic 例子）— 這些 ASM 框架不適用。
  </div>
</section>

<section id="workflow" class="anchor-target container-card border-l-4 border-green-500">
  <h2 class="text-xl font-bold mb-3">0.5 完整研究與實作流程 narrative（SCQA）</h2>
  <p class="text-sm text-slate-700 mb-3">給 reviewer / 自己未來追溯時的 self-contained 説明 — 從 question 提出到 answer 驗證的完整路徑。</p>

  <details open class="mb-2">
    <summary><strong>Situation (場景)</strong></summary>
    <p class="text-sm mt-2">2026-05-27 reviewer 提問：「Do you see any example with reconstructed somatic haplotype in tumor suppressor gene promoters that shows significant DNA methylation changes (hypomethylation) against the germline haplotype?」這是 paper revision review 的 follow-up — reviewer 要 <strong>single concrete example</strong>（"any example"），不是 cohort statistics。</p>
  </details>

  <details class="mb-2">
    <summary><strong>Complication (複雜性 — 為何 not trivial)</strong></summary>
    <ol class="text-sm mt-2 list-decimal list-inside space-y-1">
      <li>資料維度大 — HCC1395 paired BAM 260 GB / 21.7M reads / 315 TSG</li>
      <li>ASM 框架本質性限制 — HCC1395 hyper-diploid (ploidy 2.85)，51% genome LOH，<strong>BRCA1/TP53/CDH1 全在 LOH 內</strong></li>
      <li>Reviewer expectation vs 我們框架 mismatch — reviewer 可能想看 bulk silencing，但 ASM 是 allele-level</li>
      <li>Effect size unknown — 不知道剩餘非-LOH TSG 有沒有達 "significant" 程度</li>
      <li>CpG destruction artifact 風險 — 變異剛好破壞 CpG 會看起來像 hypomethylation 但不是真的</li>
    </ol>
  </details>

  <details class="mb-2">
    <summary><strong>Question (核心)</strong></summary>
    <p class="text-sm mt-2">在 95 個 curated TSG 的非-LOH promoter 中，有沒有一個 SEQC2 high-confidence TP somatic SNV 其 reconstructed somatic haplotype (HP:Z:1-1 / 2-1) 相對 germline (HP:Z:1 / 2) 在 ±1kb 顯示統計顯著的甲基化降低 (Δβ ≤ -0.05)，且不是 random null / CpG destruction artifact？</p>
  </details>

  <details open class="mb-3">
    <summary><strong>Answer (解答) — Yes, BRCA2 chr13:32315128</strong></summary>
    <p class="text-sm mt-2 text-green-700 font-semibold">完整 evidence chain 由 13 個 step 構成，下方表格列出每 step 的 rationale + script + output + 驗證</p>
  </details>

  <h3 class="font-semibold mt-3 mb-2">13-step pipeline（執行順序 + rationale）</h3>
  <table class="text-xs w-full">
    <thead><tr><th>#</th><th>Step</th><th>Rationale (為何這樣做)</th><th>Script</th><th>Output</th><th>驗證</th></tr></thead>
    <tbody>
      <tr><td>0a</td><td>KB query HCC1395 + SEQC2 truth 路徑</td><td>用既有 KB 避免重 inventory</td><td>KB Read</td><td>HCC1395.md + SEQC2 truth-set.md</td><td>KB 2026-04-11</td></tr>
      <tr><td>0b</td><td>Prereq B: BAM MM/ML tags 確認</td><td>沒甲基 calling 整個分析不可能</td><td>samtools view + grep</td><td>ML:B:C confirmed</td><td>✓</td></tr>
      <tr><td>0c</td><td>Lit scout 1: HCC1395 TSG epigenetic</td><td>找 priority target + theoretical anchor</td><td>researcher subagent</td><td>literature_scout.md (13 priority)</td><td>✓</td></tr>
      <tr><td>1</td><td>95 curated TSG list (hardcoded)</td><td>COSMIC 404 + OncoKB 404 → fall back curated</td><td>scripts/01_*.py</td><td>95 TSG</td><td>✓</td></tr>
      <tr><td>2</td><td>RefSeq GTF → TSS (multi-transcript outer)</td><td>TSS = promoter ±2kb 中心</td><td>scripts/01_*.py</td><td>277 raw → 157 merged</td><td>95/95 (100%)</td></tr>
      <tr><td>3</td><td>建 promoter BED (TSS ±2 kb)</td><td>Reviewer "promoter" 標準</td><td>01_*.py</td><td>701,784 bp</td><td>✓</td></tr>
      <tr class="bg-yellow-50"><td>4</td><td>bedtools subtract LOH region</td><td><strong>User explicit 要求</strong> — LOH 內 HP1+HP2 不存在</td><td>01_*.py</td><td>309,728 bp (44.1% 留, 51 TSG lost)</td><td>✓</td></tr>
      <tr><td>5</td><td>bcftools view PASS × promoter BED</td><td>TP variant 答辯基準</td><td>01_*.py</td><td><strong>6 candidates</strong></td><td>✓</td></tr>
      <tr class="bg-yellow-50"><td>6</td><td>Pause + AskUserQuestion (yield=6 < threshold 10)</td><td>邊界決策需 user 確認</td><td>互動</td><td>User: pilot 6 + HP cov check</td><td>✓</td></tr>
      <tr><td>7</td><td>Step 3: HP:Z tag count filter</td><td>沒 germline + somatic 兩側做不了 paired test</td><td>scripts/02_*.py</td><td>2 strict pass (BRCA2 + KMT2C)</td><td>✓</td></tr>
      <tr class="bg-yellow-50"><td>8</td><td>Step 4: ISM per-haplotype Δβ + Wilcoxon</td><td>核心 hypothesis test</td><td>scripts/03_*.py</td><td><strong>BRCA2 Δβ=−0.122 p=6.1e-11</strong></td><td>✓</td></tr>
      <tr><td>9</td><td>Negative control: 5 random 非-TSG</td><td>排除 pipeline systematic noise</td><td>scripts/04_*.py</td><td>max |Δ|=0.054 (BRCA2 exceed 2.3×)</td><td>✓</td></tr>
      <tr><td>10</td><td>Mechanism check: ref context destroy CpG?</td><td>排除 trivial CpG destruction artifact</td><td>samtools faidx + Ensembl REST</td><td>AGA→AAA, <strong>NOT CpG</strong></td><td>cross-verify</td></tr>
      <tr><td>11</td><td>Lit scout 2: chr13:32315128 specific 3-layer</td><td>找 cite-ready citation + mechanism plausibility</td><td>researcher subagent</td><td>literature_scout_chr13_32315128.md (11 papers)</td><td>4 anchor</td></tr>
      <tr><td>12</td><td>視覺化: matplotlib + IGV headless</td><td>Reviewer "easy to visualize on IGV"</td><td>scripts/05/06_*.py + IGV batch</td><td>6 PNG + IGV session</td><td>✓</td></tr>
      <tr><td>13</td><td>撰報告 (.md + standalone HTML)</td><td>Deliverable</td><td>scripts/07_*.py</td><td>本檔</td><td>✓</td></tr>
    </tbody>
  </table>

  <h3 class="font-semibold mt-4 mb-2">Causal chain — 為何相信 BRCA2 result 是 real</h3>
  <pre class="text-xs bg-slate-100 p-3 rounded overflow-x-auto leading-tight">
SEQC2 truth-set TP sSNV @ chr13:32315128 (62/63 lab samples PASS + PacBio orthogonal)
                              ↓ (variant quality 強)
located in BRCA2 promoter TSS −346 bp + 在 95 curated TSG list
                              ↓ (target biological relevance)
non-LOH region (HP1+HP2 reads both >20 in paired BAM)
                              ↓ (ASM 框架可用)
HP1-1 reconstructed somatic reads = 45 (longphase-s)
                              ↓ (sufficient depth)
197 paired CpG with HP1 ≥3 and HP1-1 ≥3 reads each
                              ↓ (statistical power)
Wilcoxon signed-rank p = 6.1 × 10⁻¹¹ (effect = somatic less methylated)
                              ↓ (signal extremely unlikely by chance)
5 random non-TSG non-LOH control max |Δ| = 0.054 → BRCA2 exceeds 2.3×
                              ↓ (signal > null distribution)
trinucleotide AGA→AAA NOT CpG → NOT CpG-destruction artifact
                              ↓ (mechanism plausibility passes)
in Fraile-Bethencourt 2018 PMID 29766361 experimentally-active scan window
                              ↓ (regional regulatory activity established)
analogous to Evans 2018 PMID 30075112 BRCA1 c.-107A>T cis-ASM precedent
                              ↓ (mechanism has published analog)
matches Do & Tycko 2020 PMID 32594908: 6-17% somatic mutations → cis-ASM in cancer
                              ↓ (general cancer biology framework)
✅ VERDICT: statistically + mechanistically + literarily-supported example
</pre>

  <h3 class="font-semibold mt-4 mb-2">Pre-Mortem retrospect — 3 failure modes 與 guard</h3>
  <ol class="text-sm list-decimal list-inside space-y-1">
    <li><strong>TVAF=0.189 subclonal</strong> → caveat 明示「subclone-level，非 bulk silencing」</li>
    <li><strong>HP2-1 = 0</strong> → biology (somatic 物理上單側) not bug</li>
    <li><strong>n=1 example</strong> → 明示「sample of one」+ cross-sample extension plan</li>
    <li><strong>Lister 0.2 ribbon 未達</strong> → effect 0.122 small-moderate，依 p value + max |Δ|=0.93 + null exceed 三線共同支持</li>
  </ol>
</section>

<section id="reviewer-q" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">1. Reviewer 問句拆解</h2>
  <blockquote class="text-sm">
    <strong>Original</strong>: "do you see any example with reconstructed somatic haplotype in tumor suppressor gene promoters that shows significant DNA methylation changes (hypomethylation) against the germline haplotype? Should be very easy to visualize on IGV."
  </blockquote>
  <table class="text-sm w-full mt-3">
    <thead><tr><th>條件</th><th>拆解</th><th>我們對應</th></tr></thead>
    <tbody>
      <tr><td><strong>reconstructed somatic haplotype</strong></td><td>longphase-s HP:Z:1-1 / 2-1 tagged reads</td><td>paired BAM HP:Z 5-class tag</td></tr>
      <tr><td><strong>TSG promoter</strong></td><td>Canonical TSG × TSS ±2 kb</td><td>95 curated TSG × RefSeq TSS ±2 kb</td></tr>
      <tr><td><strong>significant + hypomethylation</strong></td><td>統計顯著 + 方向 = somatic 少甲基</td><td>per-CpG Wilcoxon paired + Cohen's d + Lister 0.2</td></tr>
      <tr><td><strong>against germline haplotype</strong></td><td>same allele 比較</td><td>paired CpG-level HP1 vs HP1-1 test</td></tr>
      <tr><td><strong>IGV visualize</strong></td><td>IGV session 易展示</td><td>自動產 .xml + matplotlib read-level stripe</td></tr>
    </tbody>
  </table>
</section>

<section id="pipeline" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">2. Pipeline 5-step</h2>
  <pre class="text-xs bg-slate-100 p-4 rounded overflow-x-auto leading-tight">
DATA SOURCES → SEQC2 TP VCF (40K) · SEQC2 LOH bed (320 regions / 1490 Mb) · RefSeq GTF · UCSC CpG · 260GB paired BAM
                              ↓
Step 0:  95 curated TSG (COSMIC CGC ∪ Vogelstein 138 ∪ breast DDR)
                              ↓
Step 1:  Build TSG promoter BED  (TSS ± 2 kb, 157 merged intervals)
                              ↓
Step 2:  Subtract LOH regions   (44.1% bp retained, 51 TSG lost)
                              ↓
Step 3:  Intersect TP sSNV (HC + PASS)  →  6 candidates
                              ↓
Step 4:  HP:Z tag filter (HP1≥10 & HP2≥10 & HPn-1≥5)  →  2 (BRCA2, KMT2C)
                              ↓
Step 5:  ISM per-haplotype β (paired Wilcoxon ±1kb)
              ├─ BRCA2  Δβ=−0.122  p=6.1e-11  ✓ STRONG
              └─ KMT2C  Δβ=+0.002  p=0.98     ○ baseline unmethyl
                              ↓
Negative Control: 5 random non-TSG non-LOH → max |Δ|=0.054
                              ↓
        BRCA2 |Δ|=0.122 > 2.3 × control max → real signal
</pre>
</section>

<section id="inventory" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">3. Data Sources Inventory</h2>
  <table class="text-sm w-full">
    <thead><tr><th>資料</th><th>路徑 / 來源</th><th>用途</th><th>last verify</th></tr></thead>
    <tbody class="text-xs">
      <tr><td>SEQC2 TP sSNV</td><td><code>/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_*.vcf.gz</code></td><td>TP variant ground truth</td><td>KB 2026-04-11</td></tr>
      <tr><td>SEQC2 LOH bed</td><td><code>.../CNV/ngs_benchmark_cnvs_gain_loss_loh.bed</code></td><td>LOH exclusion</td><td>KB 2026-04-11</td></tr>
      <tr><td>HCC1395 paired tagged BAM</td><td><code>/big7_disk/.../longphase_s/HCC1395_tagged.bam</code> (260 GB)</td><td>HP:Z + MM/ML</td><td>本 report 2026-05-27</td></tr>
      <tr><td>RefSeq hg38 GTF</td><td>UCSC hgdownload (download 2026-05-27)</td><td>TSS 座標</td><td>2026-05-27</td></tr>
      <tr><td>UCSC CpG island</td><td>UCSC hgdownload (download 2026-05-27)</td><td>CpG island annotation</td><td>2026-05-27</td></tr>
      <tr><td>95 curated TSG</td><td>hardcoded (COSMIC CGC + Vogelstein 2013 + breast DDR)</td><td>Target universe</td><td>2026-05-27</td></tr>
      <tr><td>Lit scout (13 priority TSG)</td><td><code>output/literature_scout.md</code></td><td>Theoretical framework</td><td>2026-05-27</td></tr>
    </tbody>
  </table>
</section>

<section id="yield" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">4. Yield Funnel 95 → 1</h2>
  <table class="text-sm w-full">
    <thead><tr><th>階段</th><th>數量</th><th>備註</th></tr></thead>
    <tbody>
      <tr><td>Curated TSG list</td><td class="text-right">{summary['tsg_list_size']}</td><td>COSMIC CGC ∪ Vogelstein 138 ∪ breast DDR</td></tr>
      <tr><td>RefSeq GTF TSS 找到</td><td class="text-right">{summary['tsg_in_refseq']}/{summary['tsg_list_size']}</td><td>100% 命中</td></tr>
      <tr><td>Promoter intervals (raw)</td><td class="text-right">{summary['promoter_intervals_raw']}</td><td>multi-transcript TSS</td></tr>
      <tr><td>Promoter intervals (merged)</td><td class="text-right">{summary['promoter_intervals_merged']}</td><td>bedtools merge</td></tr>
      <tr><td>Promoter bp 原始</td><td class="text-right">{summary['promoter_bp_raw']:,} bp</td><td>TSS ± 2 kb</td></tr>
      <tr><td>LOH 排除後 promoter bp</td><td class="text-right stat-bad">{summary['promoter_bp_after_loh']:,} bp ({summary['promoter_bp_retention_pct']}%)</td><td>HCC1395 hyper-diploid 大 LOH 範圍</td></tr>
      <tr><td>TSG 保留</td><td class="text-right">{len(summary['tsg_retained_after_loh'])}/{summary['tsg_list_size']}</td><td>53.7% 失於 LOH</td></tr>
      <tr><td>HCC1395 TP sSNV ∩ 非-LOH TSG promoter (HC)</td><td class="text-right"><span class="tier-relaxed">{summary['tp_snv_in_tsg_promoter_nonLOH']}</span></td><td>Prereq A 邊界值 (&lt; 10 但 > 0)</td></tr>
      <tr><td>HP:Z 嚴格過濾通過</td><td class="text-right"><span class="tier-strict">2</span></td><td>BRCA2, KMT2C</td></tr>
      <tr><td>ISM 強訊號 (Δβ &lt; -0.05 + p &lt; 0.05)</td><td class="text-right"><span class="tier-strict">1</span></td><td>BRCA2 ⭐</td></tr>
    </tbody>
  </table>
</section>

<section id="lost-loh" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">4.1 失於 LOH 的 51 個 TSG（reviewer 預期看到的 classic 例子在這）</h2>
  <details>
    <summary>⭐ Classic breast / DNA repair TSG (HIGH impact on reviewer expectation) — 點開看清單</summary>
    <p class="text-sm mt-3">{', '.join(classic_breast)}</p>
    <p class="text-xs text-slate-600 mt-2">→ 全在 HCC1395 LOH 區，無法用 ASM 框架。若 reviewer 問 BRCA1 hypermethylation，誠實答「需用 bulk methylation analysis (SEQC2 EpiQC WGBS) — 不在 ASM 範疇」。</p>
  </details>
  <details>
    <summary>Cell cycle TSG (HIGH impact)</summary>
    <p class="text-sm mt-3">{', '.join(cell_cycle)}</p>
  </details>
  <details>
    <summary>Chromatin / transcription TSG (MID impact)</summary>
    <p class="text-sm mt-3">{', '.join(chromatin)}</p>
  </details>
  <details>
    <summary>Other TSG (LOW-MID impact)</summary>
    <p class="text-sm mt-3">{', '.join(other)}</p>
  </details>
</section>

<section id="hp-cov" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">4.2 · 6 → 2 candidate HP coverage</h2>
  <table class="text-sm w-full">
    <thead><tr><th>Gene</th><th>chr:pos</th><th>TVAF</th><th>HP1</th><th>HP2</th><th>HP1-1</th><th>HP2-1</th><th>HP3</th><th>untag</th><th>strict?</th><th>comment</th></tr></thead>
    <tbody>
"""
for c in step3:
    cov = c.get("coverage", {})
    pass_class = "tier-strict" if c.get("pass_strict") else "tier-fail"
    pass_text = "✅" if c.get("pass_strict") else "❌"
    comment = {
        "BRCA2": "germline + somatic 平衡",
        "AXIN1": "HP1=1 (partial LOH 邊界)",
        "NFE2L2": "HPn-1=0 (longphase-s 未重建 somatic-tag)",
        "ASXL1": "HP1=1, untag 64 — high ambig",
        "KMT2C": "balance OK 但 baseline 已 unmethyl",
        "FANCC": "HP1=0 (partial LOH 邊界)",
    }.get(c["gene"], "")
    bold = "font-bold" if c.get("pass_strict") else ""
    HTML += f"""      <tr class="{bold}"><td>{c['gene']}</td><td><code>{c['chrom']}:{c['pos']}</code></td><td class="text-right">{c['tvaf']:.3f}</td>
        <td class="text-right">{cov.get('HP:1', 0)}</td><td class="text-right">{cov.get('HP:2', 0)}</td>
        <td class="text-right">{cov.get('HP:1-1', 0)}</td><td class="text-right">{cov.get('HP:2-1', 0)}</td>
        <td class="text-right">{cov.get('HP:3', 0)}</td><td class="text-right">{cov.get('untag', 0)}</td>
        <td><span class="{pass_class}">{pass_text}</span></td><td class="text-xs">{comment}</td></tr>
"""
HTML += """    </tbody>
  </table>
  <p class="text-xs text-slate-600 mt-2">
    <strong>觀察</strong>: HPn-1 在所有候選只在單側出現 (HP1-1≠0 但 HP2-1=0)。這是 biology — somatic mutation 物理上只發生在 one haplotype，不是 phasing bug。
  </p>
</section>

<section id="brca2" class="anchor-target container-card border-l-4 border-red-500">
  <h2 class="text-xl font-bold mb-3">5. BRCA2 Deep Dive ⭐</h2>
"""

brca2_stat = next(r for r in step4 if r["candidate"]["gene"] == "BRCA2")["stats"]["HP1_vs_HP1-1"]

HTML += f"""  <h3 class="font-semibold mt-3 mb-2">5.1 變異位點生物背景（v4 final reframe 2026-05-28 — 三重 gene context）</h3>
  <table class="text-sm w-full">
    <tr><th class="text-left w-1/3">位置</th><td><code>chr13:32,315,128 (hg38)</code></td></tr>
    <tr><th class="text-left">Ref → Alt</th><td>G > A</td></tr>
    <tr><th class="text-left">TVAF / NVAF</th><td>0.189 / 0.002</td></tr>
    <tr class="bg-yellow-50"><th class="text-left">BRCA2 most-upstream TSS (Ensembl ENST00000544455)</th><td><code>chr13:32,315,086</code> (+ strand)</td></tr>
    <tr class="bg-yellow-50"><th class="text-left">Variant vs BRCA2 most-upstream TSS</th><td><strong class="stat-good">+42 bp downstream (5'UTR / proximal promoter)</strong> ⭐</td></tr>
    <tr><th class="text-left">BRCA2 NM_000059.4 (MANE Select)</th><td>~chr13:32,315,077 → variant +51 bp downstream</td></tr>
    <tr class="bg-yellow-50"><th class="text-left">ZAR1L canonical TSS (ENST00000533490, − strand)</th><td><code>chr13:32,315,363</code></td></tr>
    <tr class="bg-yellow-50"><th class="text-left">Variant vs ZAR1L canonical TSS</th><td><strong class="stat-good">−235 bp = 在 ZAR1L exon 1 / 5'UTR 內</strong> ⭐</td></tr>
    <tr><th class="text-left">ZAR1L gene body</th><td>chr13:32,303,699-32,315,363 (− strand) → variant <strong>inside gene body</strong></td></tr>
    <tr class="bg-yellow-50"><th class="text-left">Liu 2010 minimal bidirectional promoter (PMID 20202217)</th><td>BRCA2 TSS <strong>−187 to +310</strong> → variant at +42 bp <strong>inside experimentally-mapped region</strong> ⭐</td></tr>
    <tr class="bg-yellow-50"><th class="text-left">ZAR2 exon 1 (Liu 2010)</th><td>BRCA2 TSS <strong>−134 to +111</strong> → variant at +42 bp <strong>inside ZAR2 exon 1</strong> ⭐</td></tr>
    <tr><th class="text-left">在 CpG island?</th><td><strong>No</strong> — 268 bp 從最近 CpG island (chr13:32315396-32315763, "CpG: 43")</td></tr>
    <tr><th class="text-left">Trinucleotide context (ref)</th><td><code>A-G-A</code> → <code>A-A-A</code>，<strong>不是 CpG</strong>（samtools faidx + Ensembl REST 雙 verify）</td></tr>
    <tr><th class="text-left">Nearest CpG 5' / 3'</th><td>chr13:32315106 (−22 bp) / chr13:32315178 (+50 bp)</td></tr>
    <tr><th class="text-left">±500 bp CpG density</th><td>59 CpGs (mean ~17 bp 間距)</td></tr>
    <tr class="bg-blue-50"><th class="text-left">生物 relevance</th><td><strong>BRCA2</strong> = DSB HR repair, breast/ovarian cancer TSG (Vogelstein 138); <strong>ZAR1L (ZAR2)</strong> = C4-type zinc finger TF, <strong>直接 silences BRCA2 forward transcription at G0/G1</strong> (Liu 2010)</td></tr>
    <tr class="bg-blue-50"><th class="text-left">Mechanism hypothesis</th><td><strong>Cis-acting disruption of ZAR2/TF binding within bidirectional promoter → allele-specific methylation cascade</strong>（NOT CpG destruction）；analogous to BRCA1 c.-107A>T (Evans 2018)</td></tr>
    <tr class="bg-red-50"><th class="text-left">Published references to this exact variant</th><td><strong>None</strong> — ClinVar 0 / dbSNP no rsID / Ensembl variation empty → <strong>novel / undocumented private somatic event</strong>（與 TVAF=0.189 subclonal 一致）</td></tr>
    <tr class="bg-amber-50"><th class="text-left">Closest published variant in same region</th><td><strong>BRCA2 c.-26G>A SNP</strong> (Healey 2007 PMID 17945002) — Liu 2010 specifically mentions this SNP is also <strong>in ZAR2 exon 1</strong>（bidirectional-promoter SNV precedent）</td></tr>
  </table>

  <h3 id="brca2-stats" class="anchor-target font-semibold mt-4 mb-2">5.2 ±1 kb window 統計</h3>
  <table class="text-sm w-full">
    <tr><th class="text-left w-2/3">指標</th><th>數值</th></tr>
    <tr><td>總 CpG positions (window)</td><td class="text-right">412</td></tr>
    <tr><td>HP1 ≥ 3 reads qualified</td><td class="text-right">200</td></tr>
    <tr><td>HP1-1 ≥ 3 reads qualified</td><td class="text-right">208</td></tr>
    <tr><td>HP2-1 ≥ 3 reads</td><td class="text-right stat-bad">0 (by biology)</td></tr>
    <tr><td><strong>Paired CpG (HP1 + HP1-1 both ≥ 3)</strong></td><td class="text-right"><strong>197</strong></td></tr>
    <tr><td>Mean β germline-HP1</td><td class="text-right">0.215</td></tr>
    <tr><td>Mean β somatic-HP1-1</td><td class="text-right">0.093</td></tr>
    <tr><td><strong>Δβ (somatic − germline)</strong></td><td class="text-right stat-bad"><strong>−0.122</strong></td></tr>
    <tr><td>max |Δβ| single CpG</td><td class="text-right">0.931</td></tr>
    <tr><td><strong>Wilcoxon signed-rank p</strong></td><td class="text-right stat-good"><strong>6.1 × 10⁻¹¹</strong></td></tr>
    <tr><td>Mann-Whitney U p</td><td class="text-right">6.0 × 10⁻³</td></tr>
    <tr><td>Cohen's d</td><td class="text-right">−0.394 (small-to-moderate)</td></tr>
    <tr><td>Lister 2009 |Δβ| ≥ 0.2 ribbon</td><td class="text-right stat-neutral">mean 0.122 &lt; ribbon (small effect)</td></tr>
  </table>

  <details class="mt-3" open>
    <summary><strong>5.3 直白比喻（非 epigenetics 專家版）</strong></summary>
    <p class="text-sm mt-3 leading-relaxed">
      想像 BRCA2 promoter 是一間房子，<strong class="text-blue-700">germline haplotype</strong> 是「老房子」、<strong class="text-red-700">somatic-reconstructed haplotype</strong> 是「同一塊地新蓋的房子」。我們在房子前後各 1 km 範圍內裝 <strong>197 個甲基化偵測器 (CpG)</strong>。老房子的偵測器平均亮度 <strong>21.5%</strong>，新房子同一批偵測器亮度只剩 <strong>9.3%</strong> — 新房子的「甲基化漆」明顯被剝掉。Wilcoxon test 説這差異不是隨機 (<strong>p = 6 × 10⁻¹¹</strong>，6 sigma+)。Cohen's d −0.39 是「明顯但不誇張」效果。
    </p>
  </details>

  <h3 id="brca2-fig1" class="anchor-target font-semibold mt-5 mb-2">Figure 1 — per-CpG Δβ scatter (核心圖)</h3>
  <img src="data:image/png;base64,{fig1}" class="w-full rounded border" alt="BRCA2 per-CpG delta beta">
  <p class="text-xs text-slate-600 mt-1">上 panel：藍 = HP1, 紅 = HP1-1。中 panel：紅 bar = somatic 少甲基的 CpG。橘虛線 = 變異位點，綠 shading = CpG island。</p>

  <h3 id="brca2-fig2" class="anchor-target font-semibold mt-5 mb-2">Figure 2 — β 分布 violin (germ vs somatic)</h3>
  <img src="data:image/png;base64,{fig2}" class="w-full rounded border" alt="BRCA2 beta distribution">

  <h3 id="brca2-fig3" class="anchor-target font-semibold mt-5 mb-2">Figure 3 — BRCA2 vs 5 random control</h3>
  <img src="data:image/png;base64,{fig3}" class="w-full rounded border" alt="BRCA2 vs controls">
  <p class="text-xs text-slate-600 mt-1">BRCA2 (深紅) Δβ=−0.122 顯著超出 5 random control (灰 max |Δ|=0.054)。KMT2C (灰) baseline unmethyl 無差異。</p>

  <h3 id="brca2-fig4" class="anchor-target font-semibold mt-5 mb-2">Figure 4 — Read-level methylation (matplotlib IGV-style — controlled)</h3>
  <img src="data:image/png;base64,{fig4}" class="w-full rounded border" alt="BRCA2 read level">
  <p class="text-xs text-slate-600 mt-1">每行 = 1 條 read。黑 = methylated, 白 = unmethylated, 灰 = ambiguous。紅 box (HP1-1, n=47) 白點密度明顯高於藍 box (HP1, n=63)。</p>
</section>

<section id="brca2-region" class="anchor-target container-card border-l-4 border-teal-500">
  <h2 class="text-xl font-bold mb-3">5.5 甲基差異位點 region 分布 — CpG island SHORE 現象 + BRCA2 影響</h2>
  <p class="text-sm text-slate-700 mb-2"><strong>用戶提問</strong>: 甲基差異位點/區域是否影響啟動子與 BRCA2？<strong>答案</strong>: 是 — 差異<strong>不在 CpG island core，而集中在 CpG island SHORE</strong> (BRCA2 upstream −280~−870 bp = ZAR1L gene body)，經典 Irizarry 2009 (PMID 19151715) shore-methylation 現象。</p>

  <h3 class="font-semibold mt-3 mb-2">Region-level Δβ 分布 (197 paired CpG)</h3>
  <table class="text-sm w-full">
    <thead><tr><th>Functional region</th><th>n_CpG</th><th>mean Δβ</th><th>Wilcoxon p</th><th>解讀</th></tr></thead>
    <tbody>
      <tr class="bg-red-50"><td><strong>ZAR1L gene body</strong></td><td class="text-right">78</td><td class="text-right stat-bad"><strong>−0.299</strong></td><td>3.98e-10 ***</td><td>強 hypomethyl</td></tr>
      <tr class="bg-red-50"><td><strong>BRCA2 upstream promoter (MANE)</strong></td><td class="text-right">110</td><td class="text-right stat-bad"><strong>−0.212</strong></td><td>3.98e-10 ***</td><td>強 hypomethyl</td></tr>
      <tr><td>ZAR2 exon 1 (Liu 2010)</td><td class="text-right">18</td><td class="text-right">−0.133</td><td>0.021 *</td><td>中度</td></tr>
      <tr><td>Liu 2010 bidirectional promoter</td><td class="text-right">32</td><td class="text-right">−0.108</td><td>0.0076 **</td><td>中度</td></tr>
      <tr><td>BRCA2 5'UTR (MANE 下游)</td><td class="text-right">87</td><td class="text-right">−0.008</td><td>0.063 ns</td><td>無差異</td></tr>
      <tr class="bg-green-50"><td><strong>CpG island core</strong></td><td class="text-right">87</td><td class="text-right"><strong>+0.000</strong></td><td>(全 0)</td><td><strong>完全無差異</strong></td></tr>
    </tbody>
  </table>

  <h3 class="font-semibold mt-4 mb-2">三段式 spatial pattern</h3>
  <table class="text-sm w-full">
    <thead><tr><th>段</th><th>座標</th><th>n_CpG</th><th>mean Δβ</th><th>β germline → somatic</th></tr></thead>
    <tbody>
      <tr class="bg-red-50"><td><strong>Upstream shore</strong></td><td>pos &lt; 32,315,396</td><td class="text-right">80</td><td class="text-right stat-bad"><strong>−0.291</strong></td><td>~0.85 高 → ~0.05 低</td></tr>
      <tr class="bg-green-50"><td><strong>CpG island core</strong></td><td>32,315,396-763</td><td class="text-right">87</td><td class="text-right">+0.000</td><td>兩 allele 都 ~0 (active)</td></tr>
      <tr><td>Downstream 5'UTR</td><td>pos &gt; 32,315,763</td><td class="text-right">30</td><td class="text-right">−0.024</td><td>微弱</td></tr>
    </tbody>
  </table>

  <h3 class="font-semibold mt-4 mb-2">Figure 7 — Region-annotated per-CpG Δβ (核心新圖)</h3>
  <img src="data:image/png;base64,{fig7_region}" class="w-full rounded border" alt="BRCA2 region annotated delta">
  <p class="text-xs text-slate-600 mt-1">上 panel 藍點(germline)在 upstream shore 高甲基、紅點(somatic)同區低甲基；綠色 CpG island core 兩 allele 都 β≈0。下 panel 紅 bar (shore) 強 hypomethyl、綠區 (core) 無 bar。</p>

  <h3 class="font-semibold mt-4 mb-2">CpG island shore 概念 (Irizarry 2009 教學)</h3>
  <div class="bg-blue-50 border-l-4 border-blue-400 p-3 text-sm">
    <p><strong>CpG island core</strong> = 啟動子核心，正常 constitutively unmethylated（保持可轉錄）。兩 allele 在 core 都 β≈0 = 正常。</p>
    <p class="mt-1"><strong>CpG island shore</strong> = island 外緣 ≤2 kb「岸邊」。Irizarry 2009 (PMID 19151715) 發現<strong>大部分 cancer methylation 變化在 shore 而非 island core，且 shore methylation 與 gene expression 強相關</strong>。</p>
    <p class="mt-1">→ 我們差異<strong>完全落在 upstream shore (−280~−870 bp)</strong>，正是 Irizarry 框架預期 cancer ASM 出現的位置。</p>
  </div>

  <h3 class="font-semibold mt-4 mb-2">對 BRCA2 + ZAR1L 表達影響 (mechanistic hypothesis)</h3>
  <table class="text-sm w-full">
    <thead><tr><th>觀察</th><th>BRCA2 視角 (+ strand)</th><th>ZAR1L 視角 (− strand)</th></tr></thead>
    <tbody>
      <tr><td>差異區 = upstream shore</td><td>distal upstream promoter / shore</td><td>ZAR1L <strong>gene body</strong> (TSS 下游)</td></tr>
      <tr><td>Somatic allele 失去甲基</td><td>shore hypomethyl → <strong>de-repression ↑</strong></td><td>gene-body hypomethyl → <strong>↓ expression</strong></td></tr>
      <tr class="bg-blue-50"><td>推論方向</td><td><strong>BRCA2 ↑</strong></td><td><strong>ZAR1L (ZAR2) ↓</strong></td></tr>
    </tbody>
  </table>

  <pre class="text-xs bg-slate-100 p-3 rounded mt-2 overflow-x-auto">somatic allele upstream shore 失去甲基
   ├── BRCA2 視角: shore hypomethyl → BRCA2 de-repression ↑
   └── ZAR1L 視角: gene-body hypomethyl → ZAR2 (repressor) ↓
              │
              └── Liu 2010: ZAR2 silences BRCA2
                     ↓ ZAR2 ↓ → 解除對 BRCA2 抑制 → BRCA2 ↑ (與 shore 方向一致!)</pre>
  <p class="text-sm mt-2 font-semibold text-green-700">→ 兩條獨立甲基化方向 coherent 指向「somatic allele BRCA2 de-repression」— 與 Liu 2010 ZAR2-silences-BRCA2 model 吻合。</p>

  <div class="bg-amber-50 border-l-4 border-amber-500 p-3 mt-2 text-sm">
    ⚠ <strong>誠實 caveat</strong>: 我們只有 methylation data，<strong>無 expression (RNA-seq) data</strong>。表達方向是基於 Irizarry shore + Liu 2010 model 的 hypothesis，需 allele-specific RNA-seq 驗證。但兩個甲基化方向的 coherence 本身是 supporting evidence。
  </div>

  <div class="bg-green-50 border-l-4 border-green-500 p-3 mt-2 text-sm">
    <strong>對 reviewer "promoter" 要求補強</strong>: 32 個差異 CpG 落在 Liu 2010 experimentally-defined bidirectional promoter (−187/+310) 內 (mean Δβ=−0.108, p=0.0076 **)；其中 chr13:32,314,959 + 32,314,914 直接在 ZAR2 exon 1 內 (Δβ −0.76)。→ 不只「附近」，而是真有顯著 hypomethylated CpG 落在 published functional promoter 內。
  </div>
</section>

<section id="brca2-igv" class="anchor-target container-card border-l-4 border-purple-500">
  <h2 class="text-xl font-bold mb-3">5.7 真正 IGV 視覺化 — <code>ZAR1L_ASM_ex.xml</code> session 擷取</h2>
  <p class="text-sm text-slate-700 mb-2">Reviewer 説 "should be very easy to visualize on IGV" — 我們用既有 <code>/big8_disk/liaoyoyo2001/IGV_session/ZAR1L_ASM_ex.xml</code> IGV session 執行 headless IGV batch（locus 設於 chr13:32314183-32316071 正好涵蓋 BRCA2 variant），對 <code>sort BASE chr13:32315128</code> 後擷取 PNG。</p>

  <details class="mb-3">
    <summary><strong>Headless 執行命令（無 Xvfb，純 Java AWT headless）</strong></summary>
    <pre class="text-xs bg-slate-100 p-3 rounded mt-2 overflow-x-auto">JDK21_BIN=/big7_disk/liaoyoyo2001/jdk21/jdk-21.0.7+6/bin   # 系統 Java 11 不夠，需 JDK 21
PATH=$JDK21_BIN:$PATH bash /big8_disk/liaoyoyo2001/igv/build/IGV-dist/igv.sh \
    -b InterSubMod/research/tsg_promoter_asm_reviewer/scripts/08_igv_batch.txt</pre>
  </details>

  <details class="mb-3">
    <summary><strong>Batch script 內容</strong></summary>
    <pre class="text-xs bg-slate-100 p-3 rounded mt-2 overflow-x-auto">new
genome hg38
load /big8_disk/liaoyoyo2001/IGV_session/ZAR1L_ASM_ex.xml
goto chr13:32314183-32316071
sort BASE chr13:32315128         # ← 用戶 explicit 要求 sort base on variant
snapshotDirectory .../figures
maxPanelHeight 4000
snapshot brca2_igv_zar1l_session_full.png
goto chr13:32314928-32315328     # ±200bp zoom
sort BASE chr13:32315128
snapshot brca2_igv_zar1l_session_zoom.png
exit</pre>
  </details>

  <h3 class="font-semibold mt-3 mb-2">Figure 5 — IGV session full window (chr13:32314183-32316071, ~1.9 kb)</h3>
  <p class="text-xs text-slate-600 mb-1">17 tracks 渲染後 1150 × 3776 px composite。Embed 自動縮放；建議下載 PNG 全螢幕 view。</p>
  <img src="data:image/png;base64,{fig5_igv_full}" class="w-full rounded border" alt="BRCA2 IGV full window from ZAR1L_ASM_ex session" style="max-height: 1600px; object-fit: contain;">

  <h3 class="font-semibold mt-4 mb-2">Figure 6 — IGV session zoom (chr13:32314928-32315328, ±200 bp)</h3>
  <img src="data:image/png;base64,{fig6_igv_zoom}" class="w-full rounded border" alt="BRCA2 IGV zoom ±200bp from ZAR1L_ASM_ex session" style="max-height: 1600px; object-fit: contain;">

  <h3 class="font-semibold mt-4 mb-2">圖中 8 個 panel 的逐個解讀（由上至下）</h3>
  <table class="text-xs w-full">
    <thead><tr><th>Panel</th><th>attributeKey (XML)</th><th>~比例</th><th>看什麼</th></tr></thead>
    <tbody>
      <tr><td>P0</td><td>DataPanel + HCC1395BL methyl phase VCF (×2)</td><td>7%</td><td><strong>HCC1395BL germline phased methylation calls</strong>。色塊密度顯示 phase 後甲基化分布</td></tr>
      <tr><td>P1</td><td>HCC1395_Tmode_tagged_ClairS_pileup BAM</td><td>16%</td><td><strong>HCC1395 tumor BAM (HKU mod + ClairS-pileup tag)</strong>，SQUISHED，render <code>BASE_MODIFICATION_2COLOR</code> (<code>groupByOption=PHASE</code>)。紫/紅點 = 5mCG 高機率；藍/綠 = 低機率</td></tr>
      <tr><td>P2</td><td>HCC1395BL_ONT_5khz tagged BAM (normal HKU MOD)</td><td>15%</td><td><strong>HCC1395BL (normal) BAM</strong>，groupBy <code>TAG</code>/<code>HP</code>。Normal control — 應顯示 standard germline 雙倍體 HP1/HP2，甲基分布對稱</td></tr>
      <tr class="bg-yellow-50"><td>P3 ⭐ KEY</td><td>HCC1395 Tumor Dorado tagged BAM</td><td>15%</td><td><strong>核心 Panel</strong> — HCC1395 tumor Dorado tagged, <strong>groupBy <code>PHASE</code> = HP1 / HP2 / HP1-1</strong>。HP1-1 群在 ±1 kb 紫點密度 <strong>明顯低於 HP1</strong> — 正是 reviewer 要看的視覺證據</td></tr>
      <tr><td>P4</td><td>alignment-sort-hcc1395bl_GeneralMode_tagged BAM</td><td>11%</td><td>另一條 normal control（不同 longphase-s mode）— 一致性檢查</td></tr>
      <tr><td>P5</td><td>HC_SEQC2.bed + sSNV + sINDEL truth</td><td>14%</td><td><strong>SEQC2 truth annotation</strong>。chr13:32315128 處顯示 1 個 sSNV bar</td></tr>
      <tr><td>P6</td><td>junctions (small)</td><td>1%</td><td>Splice junctions (unused 此區)</td></tr>
      <tr><td>P7</td><td>Reference + ClairS + DeepSomatic + RefSeq Select + fp_cluster</td><td>21%</td><td>Reference (ACGT)、cross-caller validation、BRCA2 gene track（TSS + 5'UTR 結構）、FPStreakFinder cluster</td></tr>
    </tbody>
  </table>

  <h3 class="font-semibold mt-4 mb-2">對齊 reviewer Q 的 visual evidence 3 點</h3>
  <ol class="text-sm list-decimal list-inside space-y-1">
    <li><strong>Panel P3 HP1-1 群</strong> (Tumor Dorado, sort by base @ chr13:32315128)：中央橘虛線 (變異位置) 有 ALT base (A)，群上方 5mCG 紫點分布<strong>稀疏</strong> — 對應量化 β=0.093</li>
    <li><strong>Panel P3 HP1 群</strong> (germline-tagged) 紫點分布<strong>密</strong> — 對應量化 β=0.215</li>
    <li><strong>Panel P5 sSNV track</strong> 顯示 variant 在 SEQC2 truth set 是 HighConf — 變異本身可靠</li>
  </ol>

  <div class="bg-amber-50 border-l-4 border-amber-500 p-3 mt-3 text-sm">
    ⚠ <strong>限制</strong>: IGV PNG 是 7 panel 堆疊。embed 後在普通螢幕上 panel 細節較難看清。建議：
    <ul class="list-disc list-inside mt-1">
      <li>下載 PNG 全螢幕 view（路徑見最下方）</li>
      <li>或本地 IGV GUI 載入原 XML session 互動 view</li>
      <li>Figure 4 (matplotlib read-level) 是控制過的 IGV-style 替代，trend 觀察優先參考</li>
    </ul>
  </div>

  <h3 class="font-semibold mt-5 mb-2">5.7.1 Normal sample reads cut-off 修復（2026-05-27 22:26）</h3>
  <p class="text-sm text-slate-700 mb-2"><strong>用戶觀察</strong>: 「HCC1395BL 與 HCC1395_normal 顯示不完全，不夠高到包含」— 確認 BAM 本身 phase 均衡，問題在 IGV session XML Panel height 太小切掉 HP2 reads。</p>

  <table class="text-sm w-full mb-2">
    <thead><tr><th>驗證項</th><th>結果</th></tr></thead>
    <tbody>
      <tr><td>HCC1395BL_ONT_5khz @ chr13:32314928-32315328 BAM reads</td><td><strong>HP1=20, HP2=22</strong> (total 42 reads) ← phase 均衡</td></tr>
      <tr><td>HCC1395_Normal Dorado @ 同 region BAM reads</td><td><strong>HP1=16, HP2=17</strong> (total 33 reads) ← phase 均衡</td></tr>
      <tr class="bg-red-50"><td>原 session XML normal panel heights</td><td>Panel1751228868797 (HCC1395BL) = <strong>232 px</strong>; Panel1751390524140 (Normal Dorado) = <strong>163 px</strong> ← 太小，HP2 群被切掉</td></tr>
    </tbody>
  </table>

  <p class="text-sm font-semibold mt-2">修復方案（不動 user 原 XML，複製 + patch）:</p>
  <pre class="text-xs bg-slate-100 p-3 rounded">產物: config/ZAR1L_ASM_ex_PATCHED.xml
- Panel1751228805051 (Tumor ClairS):  246 → 700 px
- Panel1751228868797 (BL Normal 5khz): 232 → 600 px  ← 修復重點
- Panel1751390469773 (Tumor Dorado):   192 → 700 px
- Panel1751390524140 (Normal Dorado):  163 → 600 px  ← 修復重點
+ 同步重算 dividerFractions
+ 移除 groupByPos="chr1:..." (cross-chr stale bug)
+ 修 basemodFilter="m," → "m" (trailing comma)
</pre>

  <p class="text-sm mt-2"><strong>Batch script</strong> (<code>scripts/08_igv_batch_v3.txt</code>):</p>
  <pre class="text-xs bg-slate-100 p-3 rounded">new
maxPanelHeight 12000
load .../config/ZAR1L_ASM_ex_PATCHED.xml   ← 不是 user 原檔
goto chr13:32314183-32316071
sort BASE chr13:32315128
snapshot brca2_igv_zar1l_session_full.png
goto chr13:32314928-32315328
sort BASE chr13:32315128
snapshot brca2_igv_zar1l_session_zoom.png
exit</pre>

  <div class="bg-green-50 border-l-4 border-green-500 p-3 mt-3 text-sm">
    ✅ <strong>結果</strong>: Figure 5/6 都已更新為 patched 版本 — HCC1395BL normal phase 1+2 兩群完整顯示, HCC1395_Normal Dorado phase 1+2 兩群完整顯示, Tumor panels 1/1-1/2/2-1 群清楚。
  </div>
  <p class="text-xs text-slate-600 mt-2">如果要 patch user 原檔（但 read-only chenhan112 owner）：可請 owner 把 Panel height 改大（建議 600-800 px），或用 own copy at <code>InterSubMod/research/tsg_promoter_asm_reviewer/config/ZAR1L_ASM_ex_PATCHED.xml</code>。</p>

  <h3 class="font-semibold mt-5 mb-2">5.8 ZAR1L_ASM_ex.xml session 格式診斷</h3>
  <p class="text-sm text-slate-700 mb-2">用戶提問「為何有問題」— 對 2026-05-27 22:01 更新版 session 跑 5 項驗證 + headless IGV 試 load</p>

  <table class="text-sm w-full">
    <thead><tr><th>驗證項</th><th>結果</th><th>詳細</th></tr></thead>
    <tbody>
      <tr><td>1. XML well-formed</td><td class="stat-good">✅ PASS</td><td>xmllint 通過</td></tr>
      <tr><td>2. Resource 路徑存在性</td><td class="stat-good">✅ 12/12 OK</td><td>所有 Resources 都存在</td></tr>
      <tr><td>3. BAM index (.bai)</td><td class="stat-good">✅ 4/4 OK</td><td>所有 BAM 都有 .bai</td></tr>
      <tr class="bg-yellow-50"><td>4. VCF index</td><td>⚠ 3/6 MISS</td><td>3 個 uncompressed VCF 沒 .idx (載入慢 10-30s auto-build):<br><code>sSNV_HC.vcf</code> / <code>output_PASS.vcf</code> / <code>filtered_snv_tp.vcf</code></td></tr>
      <tr><td>5. Panel ⇄ dividerFractions</td><td class="stat-good">✅ PASS</td><td>7 panels = 6 dividers + 1</td></tr>
      <tr class="bg-yellow-50"><td>6. RenderOptions 異常</td><td>⚠ 3 處</td><td>見下面 A/B/C 詳述</td></tr>
    </tbody>
  </table>

  <h4 class="font-semibold mt-4 mb-2">⚠ 主要問題 3 處</h4>

  <div class="bg-yellow-50 border-l-4 border-yellow-500 p-3 mb-2 text-sm">
    <strong>A · 跨染色體 <code>groupByPos</code> 殘留</strong> 🟡 Low<br>
    <pre class="text-xs bg-white p-2 mt-1 rounded">&lt;RenderOptions ... groupByOption="PHASE" groupByPos="chr1:70037563-70037564" ...&gt;</pre>
    Session locus 在 chr13，但 <code>groupByPos="chr1:..."</code> 是先前 chr1 phasing test 殘留的 stale state。<br>
    <strong>影響</strong>: 載入時 IGV log 出現 <code>SEVERE [ReferenceFrame] Null chromosome:</code>，但會 fallback 到 PHASE grouping 不影響 render<br>
    <strong>修正</strong>: 移除這個 attribute（保留 <code>groupByOption="PHASE"</code> 即可）
  </div>

  <div class="bg-yellow-50 border-l-4 border-yellow-500 p-3 mb-2 text-sm">
    <strong>B · <code>basemodFilter</code> 尾巴多餘逗號</strong> 🟢 Trivial<br>
    <pre class="text-xs bg-white p-2 mt-1 rounded">basemodFilter="m,"   ← 3 處出現</pre>
    應該是 <code>basemodFilter="m"</code> 或 <code>basemodFilter="m,h"</code><br>
    <strong>影響</strong>: parser warning 但不影響功能<br>
    <strong>修正</strong>: <code>s/basemodFilter="m,"/basemodFilter="m"/g</code>
  </div>

  <div class="bg-yellow-50 border-l-4 border-yellow-500 p-3 mb-2 text-sm">
    <strong>C · Dangling Track ID（Junctions track 找不到主 BAM）</strong> 🟢 Trivial<br>
    2 個 Junctions track 引用 BAM 不在 <code>&lt;Resources&gt;</code>:
    <ul class="list-disc list-inside text-xs mt-1">
      <li><code>...alignment-sort-hcc1395_clairS_v041_ssrs_pileup_woFilter.bam</code></li>
      <li><code>...alignment-sort-hcc1395bl_GeneralMode_tagged.bam</code></li>
    </ul>
    Disk 上檔案存在但 Resources 沒列。<br>
    <strong>影響</strong>: IGV silently skip 該 junction track (visible="false")，不影響主 alignment<br>
    <strong>修正</strong>: 移除這兩個 Junction Track 元素，或加 Resources 註冊
  </div>

  <h4 class="font-semibold mt-4 mb-2">嚴重性評級</h4>
  <table class="text-sm w-full">
    <thead><tr><th>Issue</th><th>Severity</th><th>是否需修</th></tr></thead>
    <tbody>
      <tr><td>A (chr1 groupByPos)</td><td>🟡 Low</td><td>Cosmetic — log 變乾淨；不修也能 render</td></tr>
      <tr><td>B (basemodFilter 尾逗號)</td><td>🟢 Trivial</td><td>不修也行</td></tr>
      <tr><td>C (dangling junctions)</td><td>🟢 Trivial</td><td>IGV silently skip 無功能影響</td></tr>
      <tr><td>3 個 VCF 缺 .idx</td><td>🟡 Performance</td><td>不修載入慢 10-30s；<code>igvtools index &lt;vcf&gt;</code> 預建</td></tr>
    </tbody>
  </table>

  <div class="bg-green-50 border-l-4 border-green-500 p-3 mt-3 text-sm">
    <strong>結論</strong>: Session XML <strong>能正常 load + render</strong>（headless 已驗證），3 個 警告都是 cosmetic 不影響核心視覺輸出。若要乾淨版可 apply 修正 A+B+C + <code>igvtools index</code> 加速載入。
  </div>
</section>

<section id="negctrl" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">6. Negative Control (5 random non-TSG non-LOH TP sSNV)</h2>
  <table class="text-sm w-full">
    <thead><tr><th>Control</th><th>chr:pos</th><th>TVAF</th><th>n paired CpG</th><th>germ β</th><th>som β</th><th>Δ</th><th>Wilcoxon p</th></tr></thead>
    <tbody>
"""
for i, c in enumerate(ctrl, 1):
    if c['n_paired'] >= 5:
        HTML += f"""      <tr><td>Ctrl-{i}</td><td><code>{c['chrom']}:{c['pos']}</code></td><td class="text-right">{c['tvaf']:.3f}</td>
        <td class="text-right">{c['n_paired']}</td><td class="text-right">{c['mean_g']:.3f}</td><td class="text-right">{c['mean_s']:.3f}</td>
        <td class="text-right">{c['delta']:+.3f}</td><td class="text-right">{c['w_p']:.2e}</td></tr>
"""
    else:
        HTML += f"""      <tr class="text-slate-400"><td>Ctrl-{i}</td><td><code>{c['chrom']}:{c['pos']}</code></td><td class="text-right">{c['tvaf']:.3f}</td>
        <td class="text-right">{c['n_paired']}</td><td colspan="4" class="text-center">insufficient paired CpG</td></tr>
"""
HTML += """    </tbody>
  </table>
  <div class="bg-green-50 border-l-4 border-green-500 p-3 mt-3 text-sm">
    <strong>Summary</strong>:
    <ul class="list-disc list-inside">
      <li>5 random control max |Δ| = <strong>0.054</strong></li>
      <li>0/5 達 p<0.05 且 |Δ|≥0.05</li>
      <li><strong>BRCA2 |Δ|=0.122 超出 max control 2.3 ×</strong> → 不在 null 分布內</li>
      <li>結論：BRCA2 result 真實，不是 ASM pipeline 的 systematic noise</li>
    </ul>
  </div>
</section>

<section id="caveats" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">7. Caveats × 6（必對 reviewer 揭露）</h2>
  <div class="space-y-3">
    <div class="border-l-4 border-amber-500 pl-3 py-1"><strong>C1 · Sample of one</strong> — 1 strong example / 6 candidate / 95 TSG。Generalize 須擴 sample (relax TSS ± 5 kb 或加 HCC1937/COLO829)。</div>
    <div class="border-l-4 border-amber-500 pl-3 py-1"><strong>C2 · Effect size 0.122 &lt; Lister 0.2 ribbon</strong> — Cohen's d −0.39 是 small-to-moderate，**不到** biological-meaningful ribbon，但 max single-CpG |Δ|=0.93 + p=6e-11 表示訊號真實但分散。</div>
    <div class="border-l-4 border-amber-500 pl-3 py-1"><strong>C3 · 單側 HPn-1 (HP1-1 only)</strong> — 所有 candidate HPn-1 集中 HP1 側 (HP2-1=0)，**是 biology** (somatic mutation 物理上在單側) 不是 bug。只能做 same-allele 比較。</div>
    <div class="border-l-4 border-amber-500 pl-3 py-1"><strong>C4 · TVAF=0.189 是 subclonal</strong> — Methylation differential 是「subclonal allele vs bulk germline」，**不能 generalize 到 BRCA2 在 HCC1395 全 silencing**。實際意義: BRCA2-mutated subclone 顯示 hypomethylation。</div>
    <div class="border-l-4 border-red-500 pl-3 py-1"><strong>C5 · LOH 排除使 classic TSG (BRCA1/TP53/CDH1/RASSF1A) 全不在範圍</strong> — Methodological boundary。要分析這些須用 **bulk methylation (不分 HP)**。SEQC2 EpiQC WGBS 已有 raw data。</div>
    <div class="border-l-4 border-amber-500 pl-3 py-1"><strong>C6 · 95 TSG list 是 curated</strong> — COSMIC public 404 / OncoKB GitHub 404；fall back curated 95 TSG (Vogelstein + breast DDR + COSMIC CGC well-known)。完整 CGC ~315 TSG，long-tail 未覆蓋。</div>
  </div>
</section>

<section id="tradeoff" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">8. Methodological Trade-off Diagram</h2>
  <pre class="text-xs bg-slate-100 p-4 rounded overflow-x-auto leading-tight">
                       HCC1395 全基因組 TSG (95)
                       ┌──────────┴──────────┐
                       ▼                      ▼
                  在 LOH 區內             非 LOH 區
                  51 / 95 (53.7%)         44 / 95 (46.3%)
                       │                      │
                  HP1 + HP2 失去 het      HP1 + HP2 共存
                       │                      │
              ASM 框架  ✗ 不適用     ASM 框架  ✓ 適用
                       │                      │
              ↓ classic 例子在此              ↓ 我們的 candidate
              (BRCA1 / TP53 / CDH1 /         6 → 2 strict pass
              CDKN2A / RASSF1A)              1 strong (BRCA2)
                       │
              改用 bulk methylation
              (不分 haplotype)
              SEQC2 EpiQC WGBS 有 raw
  </pre>
  <p class="text-sm mt-3 leading-relaxed">
    對 reviewer 説：「如果你想看的是 classic TSG hypermethylation silencing (例如 BRCA1 整 cell line silencing)，那是 <strong>bulk methylation 領域</strong>。我們的 ASM 答覆是『<strong>非-LOH TSG promoter 是否存在 somatic-reconstructed haplotype 顯示 differential methylation vs germline</strong>』— 答案是 yes，BRCA2 是 statistically + biologically 實質的例子。」
  </p>
</section>

<section id="response" class="anchor-target container-card border-l-4 border-blue-600 bg-blue-50">
  <h2 class="text-xl font-bold mb-3">9. Response-to-Reviewer 草稿 v2（含 mechanism + citations，2026-05-27）</h2>
  <div class="text-sm leading-relaxed bg-white p-4 rounded">
    <p>Thank you for this insightful question. We do find an example matching your criteria.</p>
    <p class="mt-3">Using <strong>longphase-s</strong> tumor-normal paired ONT haplotag output on HCC1395, we searched the 95 curated tumor suppressor genes (intersection of COSMIC CGC, Vogelstein 2013 138-driver, and breast cancer DDR pathway), restricting to promoter regions (TSS ± 2 kb). To enable allele-specific methylation (ASM) analysis we excluded LOH regions (SEQC2 CNV/LOH benchmark) since LOH removes the heterozygosity needed for paired-allele comparison — this excluded <strong>51 / 95 TSG (53.7 %)</strong>, including BRCA1, TP53, CDH1, RASSF1A, CDKN2A, which the reviewer may expect to see; these would require bulk methylation analysis (SEQC2 EpiQC WGBS reference data exists for HCC1395).</p>
    <p class="mt-3">Of the 6 SEQC2 high-confidence somatic SNVs intersecting non-LOH TSG promoters, 2 (BRCA2, KMT2C) passed minimum haplotag coverage (HP1≥10, HP2≥10, reconstructed somatic-HPn-1≥5). Of these, <strong>BRCA2 (chr13:32,315,128 G>A, TVAF=0.189)</strong> showed strong allele-specific hypomethylation: across <strong>197 paired CpGs</strong> in the ±1 kb promoter window, mean methylation β was <strong>0.215</strong> on germline-tagged reads (HP:Z:1, n=63) vs <strong>0.093</strong> on somatic-reconstructed reads (HP:Z:1-1, n=47), giving <strong>Δβ = −0.122</strong> (Wilcoxon signed-rank <strong>p = 6.1 × 10⁻¹¹</strong>, Cohen's d = −0.39). Maximum single-CpG |Δβ| was 0.93. As negative control, 5 random non-TSG non-LOH somatic SNVs produced max |Δβ| = 0.054, confirming the BRCA2 signal <strong>exceeds the null distribution by 2.3 ×</strong>.</p>
    <p class="mt-3">KMT2C in contrast showed no differential (Δβ=+0.002, p=0.98) because the KMT2C promoter is baseline-unmethylated (β=0.003 on both germline and somatic haplotypes) — there is no methylation to lose.</p>
    <p class="mt-3"><strong>Mechanism</strong>: Reference-genome inspection (GRCh38 context T-G-A-<strong>G</strong>-A-G-A-A, verified by Ensembl REST API) confirms the variant is <strong>not at a CpG dinucleotide</strong> (trinucleotide context AGA → AAA), so the observed hypomethylation <strong>cannot be explained by direct CpG destruction</strong>. Instead, the variant sits at BRCA2 TSS −346 within the BRCA2 promoter region experimentally dissected by <strong>Fraile-Bethencourt et al. (2018, PMID 29766361)</strong>, specifically between two functionally-active microdeletion segments (µdel9: −329/−280, +111% activity; µdel12: −464/−415, 78% activity). The biological plausibility that single-nucleotide changes in the BRCA promoter can induce allele-specific methylation is established by <strong>Evans et al. (2018, PMID 30075112)</strong> for the BRCA1 c.-107A>T 5'UTR variant, where the variant segregates exclusively with the methylated, silenced allele in cis. The general principle that somatic SNVs in cancer cause de novo cis-ASM at TF/CTCF binding motifs is established by <strong>Do & Tycko et al. (2020, PMID 32594908)</strong>, reporting 6–17% of per-case somatic mutations associated with ASM gains. Long-read nanopore phasing of allele-specific methylation in cancer is a peer-reviewed established methodology (<strong>O'Neill et al. 2024, PMID 39406235</strong>). We therefore propose <strong>disruption of a TF binding motif</strong> as the most likely cis-acting mechanism; candidate motif analysis at chr13:32315128 ± 100 bp via JASPAR/HOMER is recommended follow-up.</p>
    <p class="mt-3"><strong>Publication status of this specific variant</strong>: chr13:32315128 G>A is <strong>not present in ClinVar, dbSNP, or Ensembl Variation</strong> (verified 2026-05-27), consistent with a private subclonal somatic event (NVAF=0.002).</p>
    <p class="mt-3"><strong>Caveats</strong>: (a) the BRCA2 mean Δβ=0.122 is below the Lister 2009 |Δβ|≥0.2 "biologically meaningful" ribbon although highly significant; (b) the variant is subclonal (TVAF≈0.19), so the result is "a BRCA2-mutated subclone shows hypomethylation against the germline" rather than "BRCA2 is silenced in HCC1395 bulk"; (c) the asymmetric pattern (HP2-1 reads = 0) is expected — somatic mutations physically occur on one haplotype only.</p>
    <p class="mt-3">IGV visualization: an IGV session XML (color-by-HP) and a read-level methylation strip figure are provided as supplementary materials.</p>
  </div>
</section>

<section id="literature" class="anchor-target container-card border-l-4 border-purple-600">
  <h2 class="text-xl font-bold mb-3">9b. Literature Evidence Base — 3-Layer Scout (verified 2026-05-27)</h2>

  <h3 class="font-semibold mt-2 mb-2">Layer 1 · 直接 reference (chr13:32315128) — 0 hits</h3>
  <table class="text-sm w-full">
    <thead><tr><th>Source</th><th>Result</th><th>Verification</th></tr></thead>
    <tbody>
      <tr><td>ClinVar <code>chr13:32315128[CHRPOS]</code></td><td class="stat-bad">0 entries</td><td>verified 2026-05-27</td></tr>
      <tr><td>dbSNP coordinate search</td><td class="stat-bad">No rsID</td><td>verified</td></tr>
      <tr><td>Ensembl REST overlap</td><td class="stat-bad"><code>[]</code> empty</td><td>verified</td></tr>
      <tr><td>COSMIC</td><td>TBD</td><td>需手動 query</td></tr>
      <tr><td>BRCA Exchange</td><td>TBD</td><td>需手動 query</td></tr>
      <tr><td>gnomAD</td><td>TBD (預期 absent)</td><td>需手動 query</td></tr>
    </tbody>
  </table>
  <p class="text-xs text-slate-600 mt-2">→ <strong>Conclusion</strong>: 此 variant 是 novel / undocumented private somatic event</p>

  <h3 class="font-semibold mt-4 mb-2">Layer 2 · BRCA2 Promoter Region — 6 papers verified</h3>
  <table class="text-sm w-full">
    <thead><tr><th>#</th><th>Author (Year)</th><th>Journal</th><th>PMID/DOI</th><th>Tier</th><th>Relevance</th></tr></thead>
    <tbody>
      <tr class="bg-yellow-50"><td>L2-1 ⭐</td><td><strong>Fraile-Bethencourt (2018)</strong></td><td>Breast Cancer Res Treat</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/29766361/">PMID 29766361</a></td><td>L2 key</td><td>Scanned BRCA2 −938/+312；我們 −346 位置就在 µdel9/µdel12 之間 active region</td></tr>
      <tr><td>L2-2</td><td>Burke (2018)</td><td>Human Mutation</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/30204945/">PMID 30204945</a></td><td>L2</td><td>BRCA2 c.-296C>T 破壞 PAX5 binding；test c.-407G>A / c.-395C>T</td></tr>
      <tr><td>L2-3</td><td>Däster (2024)</td><td>Breast Cancer Res Treat</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/39392573/">PMID 39392573</a></td><td>L2</td><td>BRCA2 promoter methyl 5% met BC / 33% TNBC</td></tr>
      <tr><td>L2-4</td><td>Moelans (2011)</td><td>J Pathol</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/21710692/">PMID 21710692</a></td><td>L2</td><td>BRCA2 aberrant methyl ≥50% DCIS+IDC MS-MLPA</td></tr>
      <tr class="bg-yellow-50"><td>L2-5 ⭐</td><td><strong>Evans (2018)</strong></td><td>Am J Hum Genet</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/30075112/">PMID 30075112</a></td><td>L2/L3 hybrid</td><td><strong>Gold standard analog</strong>: BRCA1 c.-107A>T 單一 base → cis ASM → silencing</td></tr>
      <tr><td>L2-6</td><td>Healey (2007)</td><td>Breast Cancer Res</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/17945002/">PMID 17945002</a></td><td>L2</td><td>BRCA2 c.-26G>A 5'UTR → 2× reporter activity</td></tr>
    </tbody>
  </table>

  <h3 class="font-semibold mt-4 mb-2">Layer 3 · Mechanism (somatic SNV → cis-methylation) — 5 papers</h3>
  <table class="text-sm w-full">
    <thead><tr><th>#</th><th>Author (Year)</th><th>Journal</th><th>PMID/DOI</th><th>Tier</th><th>Relevance</th></tr></thead>
    <tbody>
      <tr class="bg-yellow-50"><td>L3-1 ⭐⭐</td><td><strong>Do & Tycko (2020)</strong></td><td>Genome Biology</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/32594908/">PMID 32594908</a></td><td>L3 canonical</td><td><strong>THE</strong> paper to cite. &gt;15K ASM regions, ASM 5-9× in cancer; 6-17% somatic mutations → de novo ASM gains</td></tr>
      <tr><td>L3-2</td><td>Onuchic (2018)</td><td>Science</td><td>PMID 30139913 (claim)</td><td>L3</td><td>ASM map 49 WGBS; haplotype-ASM ↔ SNP disrupt TF binding</td></tr>
      <tr class="bg-yellow-50"><td>L3-3 ⭐</td><td><strong>O'Neill (2024)</strong></td><td>Cell Genomics</td><td><a class="text-blue-600" href="https://pubmed.ncbi.nlm.nih.gov/39406235/">PMID 39406235</a></td><td>L3 methodology peer</td><td>POG cohort long-read nanopore phased ASM；4.46M aDMR in tumor ~5× normal — 證 long-read ASM 是 established methodology</td></tr>
      <tr><td>L3-4</td><td>Wang (2020)</td><td>Cancer Epidemiol Biomarkers Prev</td><td><a class="text-blue-600" href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7678190/">PMC7678190</a></td><td>L3</td><td>mQTL framework；somatic SNV 為 cis-mQTL</td></tr>
      <tr class="text-slate-400"><td>L3-5</td><td>Lin (2024)</td><td>Nature Genetics</td><td><a class="text-blue-600" href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11549043/">PMC11549043</a></td><td>L3 (boundary)</td><td>CpG→TpG 5mC deamination — <strong>NOT 適用</strong> 我們 G>A non-CpG context；boundary reference 用</td></tr>
    </tbody>
  </table>

  <div class="bg-red-50 border-l-4 border-red-500 p-3 mt-4 text-sm">
    <strong>🚨 引用 AVOID</strong>:
    <ul class="list-disc list-inside">
      <li>不要把 Lin 2024 Pol-ε CpG C>T 當 mechanism 引用 — 我們 variant 不在 CpG</li>
      <li>不要 over-claim Burke 2018 PAX5 — c.-296 距 −346 有 50 bp，不是同一 motif</li>
    </ul>
  </div>

  <p class="text-xs text-slate-600 mt-3">完整 lit scout 報告: <code>InterSubMod/research/tsg_promoter_asm_reviewer/output/literature_scout_chr13_32315128.md</code></p>
</section>

<section id="repro" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">10. Reproducibility</h2>
  <table class="text-xs w-full">
    <thead><tr><th>階段</th><th>Script</th><th>Output</th></tr></thead>
    <tbody>
      <tr><td>Step 1-3</td><td><code>scripts/01_build_tsg_promoter_and_yield.py</code></td><td><code>output/tp_in_tsg_promoter_nonLOH.tsv</code> · <code>prereq_a_summary.json</code></td></tr>
      <tr><td>Step 3 HP cov</td><td><code>scripts/02_step3_hp_coverage.py</code></td><td><code>step3_hp_coverage.json</code></td></tr>
      <tr><td>Step 4 ISM</td><td><code>scripts/03_step4_ism_methylation_diff.py</code></td><td><code>step4_ism_results.json</code></td></tr>
      <tr><td>Negative control</td><td><code>scripts/04_negative_control_random_regions.py</code></td><td><code>step4_negative_control.json</code></td></tr>
      <tr><td>Visualization</td><td><code>scripts/05_brca2_visualization.py</code> + <code>06_brca2_read_level_viz.py</code></td><td><code>figures/brca2_*.png</code></td></tr>
      <tr><td>IGV session</td><td>(built-in step 06)</td><td><code>output/brca2_igv_session.xml</code></td></tr>
      <tr><td>HTML report</td><td><code>scripts/07_build_standalone_html.py</code></td><td>本檔</td></tr>
    </tbody>
  </table>
  <pre class="text-xs bg-slate-100 p-3 rounded mt-3"># Run all:
cd InterSubMod/research/tsg_promoter_asm_reviewer
for i in scripts/0*.py; do python3 $i; done</pre>
</section>

<section id="extensions" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">11. 後續延伸（依 reviewer follow-up）</h2>
  <table class="text-sm w-full">
    <thead><tr><th>Extension</th><th>工程</th><th>預期</th></tr></thead>
    <tbody>
      <tr><td>擴 random control n=5 → n=50</td><td>30 min</td><td>更 robust null distribution</td></tr>
      <tr><td>BRCA2 對 SEQC2 EpiQC bulk WGBS</td><td>1 hr</td><td>upgrade L2 → L1</td></tr>
      <tr><td>Cross-sample 驗 BRCA2 在 HCC1937 / COLO829</td><td>2-4 hr</td><td>generalize 單樣本</td></tr>
      <tr><td>Bulk methylation BRCA1/TP53 (LOH 區內)</td><td>2 hr</td><td>答 reviewer 「classic TSG 怎麼辦」</td></tr>
      <tr><td>Relax TSS ± 5 kb 重跑 Prereq A</td><td>30 min</td><td>yield 預期 15-30</td></tr>
      <tr><td>加 enhancer / DHS overlap (ENCODE)</td><td>1 hr</td><td>cis-regulatory element context</td></tr>
    </tbody>
  </table>
</section>

<section id="glossary" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">12. Glossary</h2>
  <dl class="text-sm space-y-1">
    <div><dt class="font-semibold inline">ASM</dt><dd class="inline"> — Allele-specific methylation</dd></div>
    <div><dt class="font-semibold inline">HP:Z tag (longphase-s)</dt><dd class="inline"> — 5 class: 1/2 = germline-tagged; 1-1/2-1 = reconstructed-somatic-tagged; 3 = ambiguous</dd></div>
    <div><dt class="font-semibold inline">β (beta)</dt><dd class="inline"> — 單 CpG methylation rate (0-1)</dd></div>
    <div><dt class="font-semibold inline">Δβ</dt><dd class="inline"> — (somatic) − (germline)</dd></div>
    <div><dt class="font-semibold inline">Wilcoxon signed-rank</dt><dd class="inline"> — paired non-parametric test</dd></div>
    <div><dt class="font-semibold inline">Cohen's d</dt><dd class="inline"> — effect size: |d|&lt;0.2 trivial / 0.2-0.5 small / 0.5-0.8 moderate / >0.8 large</dd></div>
    <div><dt class="font-semibold inline">Lister 2009 ribbon</dt><dd class="inline"> — |Δβ|≥0.2 視為 biologically meaningful</dd></div>
    <div><dt class="font-semibold inline">LOH</dt><dd class="inline"> — Loss of Heterozygosity</dd></div>
    <div><dt class="font-semibold inline">TVAF / NVAF</dt><dd class="inline"> — Tumor / Normal Variant Allele Frequency</dd></div>
    <div><dt class="font-semibold inline">SEQC2 EpiQC</dt><dd class="inline"> — Foox 2021 Genome Biology HCC1395 reference methylome</dd></div>
  </dl>
</section>

<section id="refs" class="anchor-target container-card">
  <h2 class="text-xl font-bold mb-3">13. References (完整 list)</h2>
  <h3 class="font-semibold mt-2 mb-1">⭐ BRCA2 Promoter & ASM Triad (核心 reviewer citations)</h3>
  <ul class="text-xs space-y-1">
    <li>⭐⭐ Do C & Tycko B. <em>Genome Biology</em> 2020 — Allele-specific DNA methylation increased in cancers. PMID 32594908 / PMC7322865</li>
    <li>⭐ Fraile-Bethencourt E. <em>Breast Cancer Res Treat</em> 2018 — Genetic dissection of BRCA2 promoter. PMID 29766361</li>
    <li>⭐ Evans D.G.R. <em>Am J Hum Genet</em> 2018 — BRCA1 c.-107A>T methylation-associated silencing. PMID 30075112</li>
    <li>⭐ O'Neill K. <em>Cell Genomics</em> 2024 — Long-read cancer ASM methodology. PMID 39406235</li>
    <li>Burke L.J. <em>Hum Mutation</em> 2018 — BRCA1/2 5' noncoding variants. PMID 30204945</li>
    <li>Däster K. <em>Breast Cancer Res Treat</em> 2024 — BRCA TNBC PARP markers. PMID 39392573</li>
    <li>Moelans C.B. <em>J Pathol</em> 2011 — BRCA2 DCIS+IDC hypermethylation. PMID 21710692</li>
    <li>Healey C.S. <em>Breast Cancer Res</em> 2007 — BRCA2 -26G>A. PMID 17945002</li>
    <li>Onuchic V. <em>Science</em> 2018 — ASM map 49 WGBS. PMID 30139913 (claim)</li>
    <li>Wang Q. <em>Cancer Epidemiol Biomarkers Prev</em> 2020 — bladder mQTL. PMC7678190</li>
  </ul>
  <h3 class="font-semibold mt-3 mb-1">SEQC2 + Methodology</h3>
  <ul class="text-xs space-y-1">
    <li>Foox J et al. <em>Genome Biology</em> 2021 — SEQC2 EpiQC HCC1395. PMID 34872606</li>
    <li>Masood D et al. <em>Genome Biology</em> 2024 — SEQC2 CNV/LOH benchmark</li>
    <li>Fang LT et al. <em>Nat Biotechnology</em> 2021 — SEQC2 somatic SNV truth set</li>
    <li>Lister R et al. <em>Nature</em> 2009 — DNA methylome at base resolution</li>
  </ul>
  <h3 class="font-semibold mt-3 mb-1">Negative-result verification (chr13:32315128 directly NOT in)</h3>
  <ul class="text-xs space-y-1">
    <li>ClinVar <code>chr13:32315128[CHRPOS]</code> — 0 entries (verified 2026-05-27)</li>
    <li>dbSNP — no matching rsID (verified)</li>
    <li>Ensembl REST <code>/overlap/region/human/13:32315127-32315129?feature=variation</code> — <code>[]</code> empty (verified)</li>
  </ul>
  <h3 class="font-semibold mt-3 mb-1">Cited but NOT applicable (boundary)</h3>
  <ul class="text-xs space-y-1 text-slate-500">
    <li>Lin S.H. <em>Nature Genetics</em> 2024 — Pol ε CpG C>T (PMC11549043). <strong>NOT 適用</strong> non-CpG context — keep as boundary ref</li>
  </ul>
</section>

<footer class="text-center text-xs text-slate-500 py-6">
  Generated by InterSubMod TSG Promoter ASM Reviewer Pipeline — 2026-05-27<br>
  Companion .md: <code>InterSubMod/docs/reports/in_progress/2026/05/20260527_HCC1395_TSG_promoter_ASM_reviewer_response_01.md</code><br>
  Source data: <code>InterSubMod/research/tsg_promoter_asm_reviewer/</code>
</footer>

</main>
</div>

<script>
// Scroll spy
const sections = document.querySelectorAll('section[id]');
const tocLinks = document.querySelectorAll('#toc a');
function updateActive() {{
  let current = '';
  sections.forEach(sec => {{
    const top = sec.getBoundingClientRect().top;
    if (top < 200) current = sec.id;
  }});
  tocLinks.forEach(a => {{
    a.classList.toggle('active', a.getAttribute('href') === '#' + current);
  }});
}}
window.addEventListener('scroll', updateActive);
updateActive();
</script>

</body>
</html>
"""

with open(OUT, "w") as f:
    f.write(HTML)
print(f"Wrote {OUT}")
print(f"Size: {len(HTML):,} chars ({len(HTML)/1024:.1f} KB)")
