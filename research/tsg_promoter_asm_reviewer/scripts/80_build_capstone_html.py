#!/usr/bin/env python3
"""
80 - CAPSTONE comprehensive HTML: full layered pipeline (BAM -> caller -> longphase-S
-> ISM -> Δβ/union) with per-layer numbers/filters/parameters + per-locus traceability
(why each locus ends up high/low) + condition->FP analysis + consensus loci +
verification-method inventory. §13 layer-A.

Output: docs/experiments/in_progress/2026/06/20260608_methylation_difference_pipeline_capstone_01.standalone.html
"""
import os, csv, json, html
import numpy as np

CS = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample")
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260608_methylation_difference_pipeline_capstone_01.standalone.html")

PF = json.load(open(f"{CS}/provenance_funnel.json"))["per_sample"]["HCC1395"]
AGG = json.load(open(f"{CS}/ism_aggregate.json"))
ITER = json.load(open(f"{CS}/iteration_history.json"))
CFP = json.load(open(f"{CS}/condition_fp_consensus.json"))


def esc(x):
    return html.escape(str(x))


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


# load HCC1395 trace loci
tp = {(r["Chr"], r["Pos"]): r for r in csv.DictReader(open(f"{EX}/HCC1395_tp/significance_summary.csv"))}
fp = {(r["Chr"], r["Pos"]): r for r in csv.DictReader(open(f"{EX}/HCC1395_fp/significance_summary.csv"))}

TRACES = [
    ("乾淨共識真 ASM", "chr1", "217943814", tp, "#15803d",
     "真 somatic 子單倍型(HP1S=48)甲基乾淨分成兩群(CramersV=1.0)+強平均差(Δβ=0.36) → 兩方法皆判高 = 真 ASM"),
    ("BRCA2/ZAR1L（分群-only）", "chr13", "32315128", tp, "#1E3A8A",
     "有分群(CramersV=0.58)但 Δβ≈0 —— somatic 差異在 HP1 vs HP1-1 內部，被 HP-merged 平均掉 → Δβ-method 漏掉，ISM 抓到（互補）"),
    ("FP-prone（Δβ-only+LOH）", "chr8", "81435384", fp, "#dc2626",
     "高Δβ(0.35)但無分群(CramersV=0)+LOH+假性小子群(HP1S=3) → Δβ-method 抓到但 ISM 正確排除；此型 9× 聚集 FP"),
    ("正確排除（無訊號）", "chr1", "1004726", tp, "#6b7280",
     "Δβ=0 + CramersV=0 + p=0.29(不顯著) → 無任何甲基差異 → 兩方法皆排除"),
]
TCOLS = [("NumReads", "覆蓋"), ("CramersV", "分群強度"), ("HPMergedDelta", "Δβ平均差"),
         ("ClusterPermanovaP", "群檢定p"), ("LabelHPPermanovaP", "HP檢定p"),
         ("Potential_LOH", "LOH"), ("HPFineN_HP1", "HP1"), ("HPFineN_HP1S", "HP1-1(somatic)"),
         ("HPFineN_HP2", "HP2"), ("Significant", "ISM判定")]

trace_html = ""
for name, ch, pos, src, col, expl in TRACES:
    r = src.get((ch, pos), {})
    cells = "".join(f'<tr><td>{esc(lab)}</td><td class="mono">{esc(r.get(c,"—"))}</td></tr>' for c, lab in TCOLS)
    trace_html += f'''<div class="trace" style="border-left-color:{col}">
      <div class="tt" style="color:{col}">{esc(name)} <span class="mono">{ch}:{pos}</span></div>
      <table class="tk">{cells}</table>
      <p class="te">{esc(expl)}</p></div>'''

# condition->FP table
cond_rows = ""
for c in CFP["A_condition_fp"]:
    hot = "hot" if (c["fp_vs_tp_OR"] or 0) >= 2 else ""
    cond_rows += (f'<tr class="{hot}"><td>{esc(c["condition"])}</td>'
                  f'<td class="num">{(c["tp_rate"] or 0)*100:.1f}%</td>'
                  f'<td class="num">{(c["fp_rate"] or 0)*100:.1f}%</td>'
                  f'<td class="num"><b>{c["fp_vs_tp_OR"]}</b></td></tr>')
combo = CFP["A_combo_highDb_noClust_LOH"]

# consensus top
cons_rows = ""
for x in CFP["B_consensus"]["top"][:20]:
    cons_rows += (f'<tr><td class="mono">{esc(x["locus"])}</td>'
                  f'<td class="num">{x["cramers"]}</td><td class="num">{x["dbeta"]:+.2f}</td>'
                  f'<td class="num">{x["n_tests"]}</td><td class="num">{x["num_reads"]}</td>'
                  f'<td>{"LOH" if x["loh"] else "—"}</td></tr>')

# layered funnel
L = PF
sig_tp = AGG["existence"]["HCC1395"]["tp"]["n_significant"]
i2 = next(i for i in ITER["iterations"] if i["id"] == "I2")

VERIF = [
    ("跨檢定一致（4 獨立 permutation）", "ClusterPermanova / LabelHP / LabelAllele / HPFine 有幾個顯著", "✅ 有效", "HCC1395 Sig 中 93% ≥3/4 檢定一致"),
    ("Dispersion artifact 旗標", "PERMANOVA 顯著是否只是 group 變異差(非位置差)", "✅ 有效（內建 FP 偵測）", "Sig 中 0 個 dispersion-warn"),
    ("Permutation FDR (BH-q)", "對 GlobalP 做 Benjamini-Hochberg 控偽發現", "✅ 有效", "I2 q≤0.10，selected 殘留 q>0.1 = 0"),
    ("條件分層（LOH/覆蓋/分群）", "高Δβ 拆 TP/FP，看哪個條件 FP-enriched", "✅ 有效（本報告 A）", "無分群 OR=8.6 / LOH OR=4.1 / 小子群 OR=5.8"),
    ("跨方法共識（CramersV × Δβ × 檢定）", "多方法同時判高 = 高信心", "✅ 有效（本報告 B）", "628 個共識真 ASM"),
    ("Placebo / collider（讀長洗牌）", "依讀長分群若也顯著 = 長度 artifact", "✅ 有效（需 ISM 重跑開啟）", "script 30/69 已驗證"),
    ("正控（imprinted DMR）", "已知 imprinted 應顯著 = sensitivity", "✅ 部分（重度 LOH 限制）", "PEG3 +0.72 驗證"),
    ("Normal-anchored cis（tumor−normal 殘差）", "甲基差是否突變造成 vs 本來就有", "✅ 有效（已做 HCC1395）", "HP_Residual cis-candidate 不判別 TP/FP"),
    ("Bootstrap 穩定度", "重採樣後是否仍顯著 = 可重現", "⚠ 需重跑（此版未開）", "Stability 欄全 0"),
    ("肉眼 read×CpG 熱圖", "視覺確認甲基分群", "✅ 基準（已有 gallery）", "100 credible locus gallery"),
]
vrows = "".join(f'<tr><td>{esc(a)}</td><td class="sub">{esc(b)}</td><td>{esc(c)}</td><td class="sub">{esc(d)}</td></tr>'
                for a, b, c, d in VERIF)

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>甲基差異位點 — 完整分層 pipeline + 逐位點 + 驗證</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1080px;margin:0 auto;padding:22px 16px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:22px 26px;margin-bottom:14px}}
header.top h1{{margin:0 0 5px;font-size:1.3rem}}.meta{{font-size:.8rem;opacity:.88}}
h2.sec{{color:var(--accent);border-bottom:2px solid var(--line);padding-bottom:5px;margin-top:26px;font-size:1.1rem}}
.layer{{background:#fff;border:1px solid var(--line);border-left:5px solid var(--accent);border-radius:10px;padding:14px 18px;margin:10px 0}}
.layer .ln{{font-weight:800;color:var(--accent)}}.layer .num{{font-size:1.25rem;font-weight:800;color:#0f172a}}
.layer .p{{font-size:.84rem;margin:.3rem 0}}.layer .par{{font-family:ui-monospace,monospace;font-size:.76rem;background:var(--bg);padding:4px 8px;border-radius:6px;color:#475569}}
.arrow{{text-align:center;color:var(--mut);font-size:1.2rem;margin:-4px 0}}
table{{width:100%;border-collapse:collapse;font-size:.83rem;margin-top:6px}}
th,td{{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left}}th{{background:#f1f5f9;font-size:.76rem}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}.mono{{font-family:ui-monospace,monospace;font-size:.76rem;color:#475569}}.sub{{color:var(--mut);font-size:.76rem}}
tr.hot{{background:#fef2f2}}
.traces{{display:grid;grid-template-columns:1fr 1fr;gap:12px}}@media(max-width:780px){{.traces{{grid-template-columns:1fr}}}}
.trace{{background:#fff;border:1px solid var(--line);border-left:4px solid;border-radius:8px;padding:10px 12px}}
.trace .tt{{font-weight:800;font-size:.86rem;margin-bottom:4px}}.tk td{{padding:2px 6px;font-size:.78rem}}.te{{font-size:.78rem;color:#374151;margin-top:6px}}
.box{{border-radius:12px;padding:14px 18px;margin:12px 0;font-size:.88rem}}
.find{{background:#fef2f2;border:1px solid #fecaca;border-left:5px solid #dc2626}}
.key{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A}}
footer{{margin-top:22px;font-size:.74rem;color:var(--mut);border-top:1px solid var(--line);padding-top:10px}}
</style></head><body><div class="wrap">

<header class="top"><h1>甲基差異位點 — 完整分層 pipeline + 逐位點技術細節 + 驗證盤點</h1>
<div class="meta">從原始 BAM → caller → longphase-S → ISM → Δβ/聯集，每層數字/篩選/參數 + 為何某位點被判甲基差異高/低 · HCC1395 worked example · 2026-06-08 · §13 layer-A</div></header>

<div class="box key">
<b>一頁總覽</b>：caller 全輸出 {L["somatic_pass"]:,} → SEQC2 benchmark 取 scored 區 + 比 truth → TP {L["tp"]:,}/FP {L["fp"]:,}/FN {L["fn"]:,} → ISM 117 欄統計 → 4 條件 gate → Significant {sig_tp} → 加 Δβ 聯集(I2) → meaningful {i2["pooled_tp_pass"]:,}(全6樣本). <b>每層都是定義明確的篩選，下方逐層 + 逐位點拆解。</b>
<br><span class="sub">⚠ 層級口徑：L2=benchmark 計數(TP {L["tp"]:,}/FP {L["fp"]:,})；L3 ISM 實際分析 ~29,723 TP/627 FP（少數位點因無覆蓋/缺 MM-ML 被 ISM 丟）；③④ 條件→FP 與共識分析以 ISM-analyzed 為分母（高Δβ TP 3413/FP 196）。</span>
</div>

<h2 class="sec">① 完整分層 pipeline（每層：數字 / 篩選法 / 參數）</h2>

<div class="layer"><span class="ln">L0 原始 BAM</span> · longphase-S haplotagged tumor BAM
  <div class="p">每條 read 帶 <b>MM/ML</b>(5mC 甲基機率) + <b>HP:Z</b> 單倍型 tag(1/2/1-1 somatic 子群)。</div>
  <div class="par">輸入：tumor reads + 甲基標籤 + 相位標籤（無篩選，原始觀測）</div></div>
<div class="arrow">▼</div>
<div class="layer"><span class="ln">L1 Caller（ClairS paired）</span> · <span class="num">{L["somatic_pass"]:,}</span> somatic SNV
  <div class="p">ClairS <b>paired 模式</b>（用 matched normal 扣 germline）呼叫全基因組 somatic SNV。</div>
  <div class="par">參數：--tumor_bam + --normal_bam（paired）· platform ont_r10_dorado_sup_5khz · PASS filter</div></div>
<div class="arrow">▼ SEQC2 benchmark（限 scored 高信賴區 + 比 truth）</div>
<div class="layer"><span class="ln">L2 longphase-S benchmark（caller vs SEQC2）</span> · TP <span class="num">{L["tp"]:,}</span> / FP {L["fp"]:,} / FN {L["fn"]:,}
  <div class="p">全輸出 {L["somatic_pass"]:,} → scored 區只剩 30,381（區外無 truth 不可評）→ 比 truth：命中=TP、未命中=FP、漏掉=FN。</div>
  <table style="margin:6px 0"><thead><tr><th>階段</th><th class="num">TP</th><th class="num">FP</th><th class="num">FN</th><th class="num">Precision</th><th class="num">Recall</th><th class="num">F1</th></tr></thead>
  <tbody>
   <tr><td>ClairS 原始</td><td class="num">29,865</td><td class="num">1,430</td><td class="num">9,582</td><td class="num">0.9543</td><td class="num">0.7571</td><td class="num">0.8443</td></tr>
   <tr><td><b>LongPhase-S</b>(用此)</td><td class="num">29,754</td><td class="num">627</td><td class="num">9,693</td><td class="num">0.9794</td><td class="num">0.7543</td><td class="num"><b>0.8522</b></td></tr>
   <tr><td>InterSubMod</td><td class="num">29,752</td><td class="num">544</td><td class="num">9,695</td><td class="num">0.9820</td><td class="num">0.7542</td><td class="num">0.8532</td></tr>
  </tbody></table>
  <div class="par">參數：truth=SEQC2 v1.2.1 HC SNV(39,447) · scored 區 only · SNV only。longphase-S <b>過濾 FP 1,430→627</b>（F1 0.8443→0.8522）；InterSubMod 再→544</div></div>
<div class="arrow">▼ ISM 逐位點 ±1kb 視窗分析</div>
<div class="layer"><span class="ln">L3 ISM 統計（117 欄/位點）→ Significant gate</span> · Significant <span class="num">{sig_tp}</span>(TP, {sig_tp/L["tp"]*100:.1f}%)
  <div class="p">每位點算 read×CpG 甲基矩陣 + read-read 距離 + 分群 + 跨 3 軸(ALT/HP/HP-family)檢定。Gate = <b>4 條件 AND</b>。</div>
  <div class="par">Gate：passed_gating AND global_p≤0.05 AND <b>CramersV≥0.1</b>(分群關聯) AND NumReads≥20 · window=1000</div></div>
<div class="arrow">▼ 加 Δβ 平均差 branch（聯集 c）+ FDR</div>
<div class="layer"><span class="ln">L4 Δβ/聯集 meaningful filter（I2）</span> · meaningful <span class="num">{i2["pooled_tp_pass"]:,}</span>(全 6 樣本 TP)
  <div class="p">聯集 = CramersV≥0.1(分群式) <b>OR</b> |Δβ|≥0.10(平均差式) + ≥2/4 檢定 + 無 dispersion + BH-q≤0.10。捕捉 ISM 漏掉的平均差位點。</div>
  <div class="par">參數：CV≥0.1 OR |HPMergedDelta|≥0.1 · ≥2/4 tests · no dispersion-warn · BH-q≤0.10 · reads≥20</div></div>

<h2 class="sec">② 逐位點：為何被判甲基差異「高 / 低 / 排除」（技術細節）</h2>
<p class="sub">4 個真實位點橫跨 CramersV(分群) × Δβ(平均差) 2×2，欄位值直接讀自 significance_summary.csv。</p>
<div class="traces">{trace_html}</div>
<div class="box key" style="margin-top:12px"><b>🔑 兩方法互補</b>：ISM CramersV 抓「離散分群」(如 BRCA2 Δβ≈0 但分群明顯)；Δβ 抓「平均偏移」(ISM 漏)。但 Δβ-only 在 LOH/無分群時 FP-prone（見 ③）。<b>聯集為覆蓋合理，但 Δβ branch 須標 FP 風險。</b><br><span class="sub">註：chr13:32315128 的「BRCA2/ZAR1L」基因名為外部註解（前 validated 工作 project_zar1l_brca2_asm 佐證，非此 CSV 欄位）。</span></div>

<h2 class="sec">③ 為何 Δβ 高的位點聚集 FP（條件→FP 統計確認）</h2>
<table>
  <thead><tr><th>高Δβ(≥0.1)位點中的條件</th><th class="num">TP 率</th><th class="num">FP 率</th><th class="num">FP/TP OR</th></tr></thead>
  <tbody>{cond_rows}</tbody></table>
<div class="box find"><b>結論</b>：高Δβ → FP <b>不是</b>低覆蓋造成(OR 0.87)，而是 <b>無分群(OR 8.6) + LOH(OR 4.1) + 假性小子群(OR 5.8)</b>。三者齊備時：<b>TP {(combo["tp_rate"]*100):.2f}% vs FP {(combo["fp_rate"]*100):.2f}%（9× 聚集 FP）</b>。機制：LOH 單倍型 + 假性 somatic tag → 平均偏移大但無真分群 → caller 也在此誤判 → Δβ 反映基因組困難度非 somatic 真實。</div>

<h2 class="sec">④ 共識高甲基差異位點（各方法都判高 + 數值高）n={CFP["B_consensus"]["n"]}</h2>
<p class="sub">CramersV≥0.1 AND |Δβ|≥0.10 AND ≥3/4 檢定 AND Significant —— 分群+平均差+多檢定皆強的真 ASM。top 20：</p>
<table>
  <thead><tr><th>locus</th><th class="num">CramersV</th><th class="num">Δβ</th><th class="num">檢定</th><th class="num">reads</th><th>LOH</th></tr></thead>
  <tbody>{cons_rows}</tbody></table>

<h2 class="sec">⑤ 驗證篩選正反結果/準確度的方法（肉眼之外，AI 已小部分確認有效性）</h2>
<table>
  <thead><tr><th>驗證方法</th><th>原理</th><th>有效性</th><th>本資料實證</th></tr></thead>
  <tbody>{vrows}</tbody></table>

<h2 class="sec">⑥ 整體判斷（讓你確認方法可行/有效/合理/可驗證）</h2>
<div class="box key">
<p>① <b>可行</b>：unmodified ISM 跑全 6 樣本×TP/FP/FN ~28 萬位點 + Δβ 聯集 post-hoc，全可重跑。</p>
<p>② <b>合理</b>：每層篩選定義明確、參數可查；逐位點可追溯「為何高/低」（見 ②）。</p>
<p>③ <b>可驗證</b>：10 種驗證方法（跨檢定/dispersion/FDR/條件分層/共識/placebo/正控/cis/bootstrap/肉眼），多數已實證有效（見 ⑤）。</p>
<p>④ <b>能判好壞</b>：好結果=共識 628(分群+平均差+多檢定+非LOH)；壞結果=Δβ-only+LOH+無分群(9× FP)。明顯現象=兩方法互補 + Δβ 的 FP-prone 條件已統計確認。</p>
</div>

<footer>數據源：provenance_funnel / ism_aggregate / iteration_history / condition_fp_consensus .json（← significance_summary.csv 源頭未改）· 腳本 70-80 · §13 layer-A · ledger 89</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[80] wrote {OUT} ({len(HTML_DOC)//1024} KB)")
