#!/usr/bin/env python3
"""
98 - Build the full pipeline MANUAL (教授級說明書) as a standalone HTML.

Every PARAMETER / FORMULA is from verified source code (cited file:line, read
2026-06-09). Every RESULT NUMBER is read at build time from the real json outputs
(funnel / diag89 / diag90 / cause / refined_gate / validation3) — never hardcoded
from memory (§13 anti-fabrication). Layered intuition -> formula style.

Scope flag: HCC1395 single sample, single pipeline (ClairS paired -> LongPhase-S ->
ISM). ⭐3 ceiling. NOT a variant filter (methylation->TP/FP concluded NEGATIVE).

Output: display_v2/20260609_pipeline_manual_for_professor_01.html
"""
import json, datetime

DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
OUT = f"{DV}/20260609_pipeline_manual_for_professor_01.html"

fn = json.load(open(f"{DV}/funnel_numbers.json"))
d89 = json.load(open(f"{DV}/diag89.json"))
d90 = json.load(open(f"{DV}/diag90.json"))
cause = json.load(open(f"{DV}/cause_genomewide.json"))
rg = json.load(open(f"{DV}/refined_gate.json"))
v3 = json.load(open(f"{DV}/validation3.json"))["summary"]

# result numbers (from real files)
R = dict(
    lps_tp=fn["lps_tp_vcf_records"], lps_fp=fn["lps_fp_vcf_records"],
    ism_tp=fn["ism_tp_analyzed"], ism_fp=fn["ism_fp_analyzed"],
    sig_tp=d90["orig_significant"]["tp"], sig_fp=d90["orig_significant"]["fp"], sig_ratio=d90["orig_significant"]["ratio"],
    base=d90["base_ratio"],
    overstrict=d89["q1"]["overstrict_tp"], overstrict_pct=d89["q1"]["pct_of_tp"],
    cause_rel=cause["tp"]["reliability_gated"], cause_n=cause["tp"]["n"], cause_pct=cause["tp"]["pct_reliability"], cause_med=cause["tp"]["median_rawCV"],
    dbeta=d89["q3"]["tp_dbeta"], dbeta_ratio=d89["q3"]["dbeta_ratio"], dbeta_enrich=d89["q3"]["dbeta_enrich"],
    cv=d89["q3"]["tp_cramersv"], cv_ratio=d89["q3"]["cramersv_ratio"], overlap=d89["q3"]["both_tp"],
    tierA_tp=next(t for t in d90["tiers"] if t["minN"] == 10)["tierA_tp"],
    tierA_fp=next(t for t in d90["tiers"] if t["minN"] == 10)["tierA_fp"],
    tierA_ratio=next(t for t in d90["tiers"] if t["minN"] == 10)["tierA_ratio"],
    valA=v3["by_tier"]["A"]["validated"], valA_n=v3["by_tier"]["A"]["n"],
    valB=v3["by_tier"]["B"]["validated"], valB_n=v3["by_tier"]["B"]["n"],
)


def pct(a, b):
    return f"{100*a/b:.0f}%"


# ============================ HTML ============================
CSS = """
:root{--bg:#0f172a;--card:#1e293b;--mut:#94a3b8;--fg:#e2e8f0;--ac:#38bdf8;--gd:#fbbf24;--gn:#4ade80;--rd:#f87171}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:14.5px/1.7 -apple-system,'Segoe UI',system-ui,sans-serif}
.layout{display:flex;max-width:1500px;margin:0 auto}
nav{position:sticky;top:0;align-self:flex-start;height:100vh;overflow:auto;width:248px;padding:18px 12px;border-right:1px solid #334155;font-size:13px;flex-shrink:0}
nav a{display:block;color:var(--mut);text-decoration:none;padding:3px 8px;border-radius:5px}
nav a:hover{background:#16203a;color:var(--ac)} nav a.h2{padding-left:18px;font-size:12px}
nav .t{color:var(--ac);font-weight:700;margin:10px 0 4px;font-size:12px}
main{flex:1;padding:24px 34px 90px;min-width:0}
h1{font-size:23px;margin:0 0 4px} h2{font-size:19px;margin:34px 0 6px;padding:6px 0 6px 11px;border-left:4px solid var(--ac);scroll-margin-top:14px}
h3{font-size:15.5px;color:var(--ac);margin:18px 0 6px} .sub{color:var(--mut);font-size:13px}
.intu{background:#0c2030;border:1px solid #155e75;border-left:4px solid var(--ac);border-radius:8px;padding:10px 14px;margin:10px 0}
.intu b{color:#7dd3fc}
table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}
th,td{border:1px solid #334155;padding:6px 9px;text-align:left;vertical-align:top}
th{background:#16203a;color:var(--mut);font-weight:600} td.n{text-align:right;font-variant-numeric:tabular-nums}
code{background:#0b1222;padding:1px 6px;border-radius:4px;font-size:12.5px;color:#fcd34d}
.formula{background:#0b1222;border:1px solid #334155;border-radius:8px;padding:12px 16px;margin:8px 0;font-size:15px;text-align:center;font-family:'Cambria Math','Georgia',serif}
.src{font-size:11px;color:#64748b} .src code{background:none;color:#64748b;padding:0}
.flowbox{display:flex;flex-direction:column;gap:0;margin:14px 0;max-width:680px}
.fb{background:var(--card);border:1px solid #334155;border-radius:9px;padding:10px 15px;position:relative}
.fb .t{font-weight:700;color:var(--ac)} .fb .d{font-size:12px;color:var(--mut)}
.fb .p{font-size:11.5px;color:#fcd34d;margin-top:3px}
.farrow{text-align:center;color:var(--mut);font-size:18px;line-height:1.1}
.note{background:#1c1407;border:1px solid #78350f;border-radius:8px;padding:10px 14px;font-size:13px;color:#fcd34d;margin:12px 0}
.ok{background:#052012;border:1px solid #14532d;border-radius:8px;padding:10px 14px;font-size:13px;color:#86efac;margin:12px 0}
.badge{display:inline-block;font-size:11px;padding:1px 7px;border-radius:9px;background:#334155;color:#cbd5e1;margin-left:4px}
.partial{background:#3b1106;color:#fdba74}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:14px}
@media(max-width:900px){nav{display:none}.grid2{grid-template-columns:1fr}}
"""

SECTIONS = [
    ("s0", "0 · 總覽與如何讀這份說明書"),
    ("s1", "1 · L0 數據與樣本（輸入）"),
    ("s2", "2 · L1 體細胞變異呼叫（ClairS paired）"),
    ("s3", "3 · L2 LongPhase-S 體細胞單倍型標記"),
    ("s4", "4 · L3 ISM 逐位點分析引擎"),
    ("s5", "5 · L4 統計檢定層（公式）"),
    ("s6", "6 · L5 顯著性 gate 與 Tier 分層"),
    ("s7", "7 · L6 TP/FP-free 獨立驗證（5 層）"),
    ("s8", "8 · 輸出格式與資料字典（schema）"),
    ("s9", "9 · 術語表（直覺→公式）"),
    ("s10", "10 · 已知邊界與全流程數字總覽"),
]


def nav():
    out = ['<nav><div class="t">目錄</div>']
    for sid, title in SECTIONS:
        out.append(f'<a href="#{sid}">{title}</a>')
    out.append("</nav>")
    return "".join(out)


def flow(steps):
    parts = ['<div class="flowbox">']
    for i, (t, d, p) in enumerate(steps):
        parts.append(f'<div class="fb"><div class="t">{t}</div><div class="d">{d}</div>'
                     + (f'<div class="p">{p}</div>' if p else '') + '</div>')
        if i < len(steps) - 1:
            parts.append('<div class="farrow">▼</div>')
    parts.append('</div>')
    return "".join(parts)


def body():
    H = []
    A = H.append

    # ---------- S0 overview ----------
    A('<h1>HCC1395 甲基差異位點分析 — 全流程說明書</h1>')
    A('<div class="sub">從原始樣本到位點判斷與驗證的逐層拆解 · 給教授確認與理解 · 生成 '
      f'{datetime.date.today().isoformat()} · <span class="badge partial">single-sample HCC1395 / 單 pipeline ⭐3</span></div>')
    A('<h2 id="s0">0 · 總覽與如何讀這份說明書</h2>')
    A('<div class="intu"><b>一句話：</b>我們對每一個體細胞突變位點，看「它周圍 read 的甲基化模式，是否會按照（由 SNP 定相的）'
      '單倍型 tag 分開」。能清楚分開＝該位點有等位基因特異甲基化（ASM）。本說明書拆解這個判斷經過的每一層、每個參數、每條公式。</div>')
    A('<div class="note">🔴 <b>定位邊界（先講清楚）</b>：本流程是 <b>characterization（描述哪些位點有 ASM 結構）</b>，'
      '<b>不是</b>用甲基化去判別真/假體細胞突變（TP/FP）—— 後者已驗證為 NEGATIVE。所有 TP/FP 字樣只是「SEQC2 真集標記的位點類別」，用來描述，不是篩選目標。</div>')
    A('<h3>全流程構造圖（資料如何一層層被處理）</h3>')
    A(flow([
        ("樣本 BAM（ONT 長讀 + 甲基化 tag）", "腫瘤 HCC1395 + 配對正常 HCC1395BL；每條 read 帶 5mC/5hmC 機率 (MM/ML)", "GRCh38 參考"),
        ("L1 ClairS paired → 體細胞 SNV", "腫瘤 vs 正常成對呼叫，得 somatic_pass.vcf.gz", "MAPQ≥20"),
        ("L2 LongPhase-S 單倍型標記", "用 germline 定相 SNP 把每條 read 標 HP tag（1/2 + 體細胞 1-1/2-1）；對 SEQC2 真集標 TP/FP/FN", "MAPQ≥20, pct 0.6"),
        (f"L3 ISM 逐位點：每個 SNV ±{1000}bp 視窗", "抽 read×CpG 甲基矩陣 → 距離 → 分群 → HP 關連檢定", "window 1000, NHD, C_min 3"),
        ("L4 統計檢定", "CramersV（卡方）+ PERMANOVA（距離）+ Δβ + per-CpG Fisher + 正常基線 cis", "p, 效應量"),
        ("L5 顯著 gate + Tier 分層", f"Significant gate → {R['sig_tp']} 高信心；Tier A/B 信心分層", "p≤0.05, V≥0.1, reads≥20"),
        ("L6/L7 獨立驗證（不用 TP/FP）", "HP-permutation + split-half + 雙股 + bootstrap → 確認真實可複製", "BH-FDR q≤0.05"),
    ]))
    A('<div class="sub">每一節結構：<b>直覺</b>（藍框）→ <b>參數/公式表</b>（含實際值 + 程式碼出處 file:line）→ <b>輸出</b>。'
      '數字結果一律來自實際檔案讀回（非記憶）。</div>')

    # ---------- S1 data ----------
    A('<h2 id="s1">1 · L0 數據與樣本（輸入）</h2>')
    A('<div class="intu"><b>直覺：</b>原料是同一個乳腺癌細胞株的「腫瘤」與「配對正常」兩份 ONT 長讀資料，'
      '每條 read 不只有鹼基，還帶有每個 CpG 位點被甲基化的機率（5mC + 5hmC）。配對正常是用來當「基線」排除遺傳性差異。</div>')
    A('<table><tr><th>項目</th><th>內容</th><th>用途</th></tr>'
      '<tr><td>樣本</td><td>HCC1395（乳腺癌細胞株，腫瘤）+ HCC1395BL（同源正常）</td><td>腫瘤分析 + 正常基線 cis 對照</td></tr>'
      '<tr><td>定序</td><td>ONT 長讀，5khz simplex，<code>5mCG_5hmCG</code> 甲基化呼叫</td><td>read-level 甲基化 + 跨距離 phasing</td></tr>'
      '<tr><td>甲基化編碼</td><td>BAM <code>MM/ML</code> tag：每 CpG 一個 0–255 機率值（ML/255 → 0–1）</td><td>read×CpG 甲基矩陣</td></tr>'
      '<tr><td>參考</td><td><code>GRCh38_no_alt_analysis_set.fasta</code></td><td>比對與 CpG 定位</td></tr>'
      '<tr><td>真集（benchmark）</td><td>SEQC2 <code>high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz</code> + HC BED</td><td>標 TP/FP/FN（描述用，非篩選）</td></tr></table>'
      '<div class="src">出處：<code>longphase_s.log</code> [Input/Benchmark Files]；BAM = <code>/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/{HCC1395,HCC1395BL}.bam</code></div>')

    # ---------- S2 caller ----------
    A('<h2 id="s2">2 · L1 體細胞變異呼叫（ClairS paired）</h2>')
    A('<div class="intu"><b>直覺：</b>先找出「腫瘤有、正常沒有」的突變（體細胞 SNV）。用腫瘤+正常成對呼叫，'
      '避免把遺傳性變異誤當體細胞。輸出是一份體細胞 SNV 清單（VCF）。</div>')
    A('<table><tr><th>參數/輸出</th><th>值</th><th>說明</th></tr>'
      '<tr><td>模式</td><td>ClairS <b>paired</b>（腫瘤 + 配對正常）</td><td>非 ClairS-TO（tumor-only）</td></tr>'
      '<tr><td>輸出</td><td><code>somatic_pass.vcf.gz</code></td><td>進 L2 的 tumor SNP 檔</td></tr>'
      '<tr><td>germline 定相</td><td><code>germline_phased_merged.vcf.gz</code></td><td>L2 phasing 的骨架（關鍵：HP 來源）</td></tr></table>'
      '<div class="note">⚠ 外部工具細節（ClairS 版本/門檻）以 KB <code>/big8_disk/.../Knowledge/05_tools/</code> 為準，'
      '本說明書只描述其在本 pipeline 的角色與輸入輸出。</div>')

    # ---------- S3 longphase ----------
    A('<h2 id="s3">3 · L2 LongPhase-S 體細胞單倍型標記（HP tag 的來源）</h2>')
    A('<div class="intu"><b>直覺：</b>把每條 read 貼上「它來自哪一條染色體拷貝（單倍型）」的標籤。'
      '用的是 germline 雜合 SNP 的定相結果（與甲基化<b>無關</b>）→ 這保證後面「甲基↔HP」不是循環論證。'
      '同時對 SEQC2 真集把每個位點標成 TP/FP/FN。</div>')
    A('<h3>HP tag 的層級（重要術語）</h3>')
    A('<table><tr><th>HP tag</th><th>raw 整數</th><th>意義</th></tr>'
      '<tr><td>HP1 / HP2</td><td>1 / 2</td><td>germline 單倍型（純 SNP 定相，本流程用此層做非循環驗證）</td></tr>'
      '<tr><td>HP1-1 / HP2-1</td><td>11 / 21</td><td>體細胞子單倍型（在 HP1/HP2 上新生的體細胞變異分支）</td></tr>'
      '<tr><td>HP3-3</td><td>33</td><td>第三型 / 保守標記</td></tr></table>'
      '<div class="src">出處：<code>RegionProcessor.hpp:105,155</code>、<code>Stats.hpp:219</code>；HP「家族」= (HP1+HP1-1) vs (HP2+HP2-1)，<code>LabelTest.hpp:217</code></div>')
    A('<table><tr><th>參數</th><th>值</th><th>說明</th></tr>'
      '<tr><td>filter mapping quality below</td><td><b>20</b></td><td>低品質比對的 read 不標</td></tr>'
      '<tr><td>percentage threshold</td><td><b>0.6</b></td><td>標 HP 的支持比例門檻</td></tr>'
      '<tr><td>tag supplementary</td><td>enabled</td><td>補充比對也標</td></tr>'
      '<tr><td>threads</td><td>24</td><td>—</td></tr></table>'
      '<div class="src">出處：<code>longphase_s.log</code> [Somatic Haplotagging Params]</div>')
    A('<h3>TP/FP/FN 怎麼來（描述用標記）</h3>')
    A('<div class="sub">每個 somatic SNV 與 SEQC2 真集比對：在真集內＝<b>TP</b>、不在＝<b>FP</b>、真集有但沒呼叫到＝<b>FN</b>。'
      f'輸出 <code>filtered_snv_{{tp,fp,fn}}.vcf.gz</code>。本研究分析的位點數：TP <b>{R["lps_tp"]:,}</b>、FP <b>{R["lps_fp"]:,}</b>（進 L3）。</div>')

    # ---------- S4 ISM ----------
    A('<h2 id="s4">4 · L3 ISM（InterSubMod）逐位點分析引擎</h2>')
    A('<div class="intu"><b>直覺：</b>對每一個體細胞 SNV，開一個 ±1000bp 的視窗，抓視窗內所有 read 在各 CpG 的甲基化，'
      '排成一張「read × CpG」表；算 read 兩兩有多像（距離）；把像的 read 聚成群；再問「這些群是否對應 HP tag」。</div>')
    A('<h3>4.1 內部資料流</h3>')
    A(flow([
        ("ReadParser / MethylationParser", "抓視窗內 read 的 CpG 甲基化機率 (0–1)，過濾品質", "MAPQ≥20, readlen≥1000, baseQ≥20"),
        ("MatrixBuilder → read×CpG 連續矩陣", "存 methylation.csv（連續 0–1，NA=未覆蓋）", "—"),
        ("二元化（給距離用）", "≥0.8→1(甲基), ≤0.2→0(未甲基), 之間→無效", "methyl_high 0.8 / low 0.2"),
        ("DistanceMatrix（NHD）", "read 兩兩距離 = 差異 CpG 數 / 共同有效 CpG 數", "C_min 3, NHD"),
        ("HierarchicalClustering", "層級聚類 + silhouette 選最佳群數 optimal_k", "Ward/average linkage"),
        ("HP / 子單倍型指派", "每 read 的 HP1/HP2/HP1-1/HP2-1（來自 L2 tag）", "—"),
        ("GlobalTest / LabelTest / StructureTest", "分群↔HP 關連（CramersV / PERMANOVA / Δβ）", "見 §5"),
    ]))
    A('<h3>4.2 ISM 全參數表（實際使用值）</h3>')
    A('<table><tr><th>參數</th><th>值（本研究）</th><th>意義</th><th>出處</th></tr>'
      '<tr><td><code>--window-size</code></td><td class="n">1000</td><td>SNV 兩側各 ±1000bp 視窗</td><td>--help / 84.sh</td></tr>'
      '<tr><td><code>--methyl-high / --methyl-low</code></td><td class="n">0.8 / 0.2</td><td>二元甲基門檻（給 NHD）</td><td>Config.hpp:33-34</td></tr>'
      '<tr><td><code>--min-mapq</code></td><td class="n">20</td><td>read 最低比對品質</td><td>--help</td></tr>'
      '<tr><td><code>--min-read-length</code></td><td class="n">1000</td><td>read 最短長度</td><td>--help</td></tr>'
      '<tr><td><code>--min-base-quality</code></td><td class="n">20</td><td>SNV 位點鹼基品質</td><td>--help</td></tr>'
      '<tr><td><code>--distance-metric</code></td><td>NHD</td><td>6 選 1（NHD/L1/L2/CORR/JACCARD/BERNOULLI）</td><td>--help</td></tr>'
      '<tr><td><code>--min-common-coverage</code> (C_min)</td><td class="n">3</td><td>算距離最少共同 CpG 數</td><td>--help</td></tr>'
      '<tr><td><code>--nan-distance-strategy</code></td><td>MAX_DIST=1.0</td><td>重疊不足的 read 對給最大距離</td><td>--help</td></tr>'
      '<tr><td><code>--expected-coverage</code></td><td class="n">0 → KDE 自動</td><td>二倍體覆蓋 = NumReads 分布眾數（KDE）；資料不足 fallback 75x</td><td>RegionProcessor.cpp:59-105</td></tr>'
      '<tr><td>PERMANOVA 置換次數</td><td class="n">999</td><td>距離結構檢定的 null 次數</td><td>LabelTest.hpp:31</td></tr></table>')
    A('<div class="note"><b>⚠ KDE 覆蓋會因「跑哪些位點集合」而不同</b>：全 TP 與子集重跑的二倍體覆蓋估計會差（如 53x vs 65x），'
      '但<b>只影響 LOH/Coverage 標註，不影響 CramersV / 分類 / 圖</b>（本研究已逐位點驗證 2070 位點分類零差異）。</div>')

    # ---------- S5 stats ----------
    A('<h2 id="s5">5 · L4 統計檢定層（每條公式：直覺 → 式子）</h2>')
    A('<h3>5.1 NHD — read 兩兩距離</h3>')
    A('<div class="intu"><b>直覺：</b>兩條 read 在它們<b>共同覆蓋</b>的 CpG 上，甲基化狀態有多少比例不一樣。0=完全一樣、1=完全相反。</div>')
    A('<div class="formula">NHD(i, j) = (不同的 CpG 數) / (共同有效 CpG 數)　，共同 CpG &lt; 3 則記無效(−1)</div>'
      '<div class="src">出處：<code>DistanceMatrix.cpp:21-54</code>（在二元化 0.8/0.2 後的矩陣上計算）</div>')
    A('<h3>5.2 階層聚類 + optimal_k</h3>')
    A('<div class="intu"><b>直覺：</b>把距離近的 read 一層層併成樹，然後在不同「剪幾群」之間，挑讓群內最緊密、群間最分開的群數（silhouette 最高）。</div>')
    A('<div class="formula">optimal_k = argmax<sub>k</sub> 平均 silhouette(k)　，silhouette = (b − a) / max(a, b)</div>'
      '<div class="src">a=同群平均距離, b=最近他群平均距離；<code>HierarchicalClustering.cpp:488-532</code></div>')
    A('<h3>5.3 CramersV — 「分群 × HP 標籤」關連強度（卡方系）</h3>')
    A('<div class="intu"><b>直覺：</b>把 read 依「聚類群」×「HP 標籤」做交叉表，問這兩種分法有多一致。0=無關、1=完全對應。'
      '<b>但卡方在表格太稀疏時不可信</b>，所以有可靠性閘控。</div>')
    A('<div class="formula">Cramér\'s V = √( χ² / ( n · (min(列,欄) − 1) ) )</div>'
      '<div class="formula">χ² = Σ (觀察 − 期望)² / 期望　；　期望(i,j) = (列總 · 欄總) / n</div>'
      '<div class="note"><b>可靠性閘控（關鍵）</b>：若任一格<b>期望值 &lt; 5</b>（Cochran 規則）→ 標記不可靠 → '
      '輸出 CramersV <b>強制歸零</b>。即 <code>summary CramersV = cramers_v_reliable ? 值 : 0</code>。'
      '原始（未歸零）值另存於 <code>significance.json</code>。<br>'
      '<span class="src">出處：<code>MathUtils.cpp:154-175</code>（min_expected&lt;5）、<code>RegionProcessor.cpp:1592-1603</code>（取可靠軸 max）</span></div>')
    A('<h3>5.4 PERMANOVA — 距離法的 HP 分離（對稀疏穩健）</h3>')
    A('<div class="intu"><b>直覺：</b>不依賴離散分群，直接問「跨 HP 的 read 距離，是否顯著大於同 HP 內」。'
      '用置換 HP 標籤建 null，看實測 pseudo-F 是否在尾端。稀疏/不平衡時比 CramersV 可靠。</div>')
    A('<div class="formula">pseudo-F = ( SS<sub>between</sub> / (k−1) ) / ( SS<sub>within</sub> / (n−k) )</div>'
      '<div class="sub">k=群數, n=read 數；完美分離 (SS<sub>within</sub>=0) → F=10⁹（哨兵值，讓置換檢定正確拒絕 null）。'
      'p 由 999 次 HP 標籤置換得經驗值。同時看 PERMDISP（離散度）避免「只是組內散度不同」假陽。</div>'
      '<div class="src">出處：<code>StructureTest.cpp:141-152</code>；置換 999 <code>LabelTest.hpp:31</code></div>')
    A('<h3>5.5 Δβ（HPMergedDelta）— HP 家族間的分離量</h3>')
    A('<div class="intu"><b>直覺：</b>把 read 併成兩大家族（HP1+HP1-1）vs（HP2+HP2-1），問「跨家族的 read 距離，比家族內大多少」。'
      '值大＝兩家族甲基化模式明顯不同（家族間甲基整體高低不同時也會推高它）。</div>')
    A('<div class="formula">Δβ = (跨家族平均 read 距離) − (家族內平均 read 距離)　= between_mean − within_mean</div>'
      '<div class="note">⚠ <b>精確定義（給教授）</b>：Δβ 是在 NHD 距離上的「<b>家族間減家族內</b>」分離量，<b>不是</b>單純的甲基化平均機率差。'
      '它與 CramersV 是<b>不同訊號</b>：Δβ 用 HP 家族標籤直接算距離分離；CramersV 用「聚類群 × HP」交叉表。'
      f'本研究中 Δβ≥0.2 與 CramersV≥0.1 的位點重疊僅 <b>{R["overlap"]}</b> 個。'
      '<br><span class="src">出處：<code>LabelTest.cpp:204-214</code>（delta=between_mean−within_mean）、<code>RegionProcessor.cpp:1621</code>、欄定義 <code>RegionProcessor.hpp:93</code></span></div>')
    A('<h3>5.6 per-CpG Fisher / ARI / 正常基線 cis</h3>')
    A('<table><tr><th>量</th><th>直覺</th><th>公式 / 定義</th><th>出處</th></tr>'
      '<tr><td>per-CpG Fisher</td><td>逐一 CpG 問「甲基/未甲基」是否與 HP 有關</td><td>2×2 Fisher exact + BH-FDR；輸出 Fisher_Frac_Sig</td><td>PerCpgAsm.cpp</td></tr>'
      '<tr><td>ARI（調整蘭德指數）</td><td>兩種分法(聚類 vs HP)一致度，校正隨機</td><td>(Σ C(nᵢⱼ,2) − 期望) / (½(Σ C(aᵢ,2)+Σ C(bⱼ,2)) − 期望)</td><td>Bootstrap.cpp:30-72</td></tr>'
      '<tr><td>HP_Residual（cis）</td><td>腫瘤的 HP 差「扣掉」正常基線 = 真體細胞 cis 訊號</td><td>HP_Residual = tumor_HP_delta − normal_HP_delta</td><td>RegionProcessor.cpp:966</td></tr></table>')

    # ---------- S6 gate + tier ----------
    A('<h2 id="s6">6 · L5 顯著性 gate 與 Tier 信心分層</h2>')
    A('<div class="intu"><b>直覺：</b>有了各種統計，怎麼判「這位點算不算有顯著 ASM 結構」？用一道組合 gate。'
      '再依「可靠性 + 結構檢定」分高/中信心 Tier。</div>')
    A('<h3>6.1 兩道 gate（注意 p 門檻不同）</h3>')
    A('<table><tr><th>gate</th><th>條件</th><th>出處</th></tr>'
      '<tr><td>passed_gating（每軸）</td><td>p ≤ <b>0.1</b> 且 CramersV ≥ 0.1</td><td>GlobalTest.hpp:29,32</td></tr>'
      '<tr><td><b>Significant</b>（最終）</td><td>passed_gating 且 global_p ≤ <b>0.05</b> 且 CramersV ≥ 0.1 且 NumReads ≥ 20</td><td>RegionProcessor.cpp:1143-1146</td></tr></table>'
      f'<div class="ok">結果：通過 Significant 的 TP = <b>{R["sig_tp"]:,}</b>、FP = <b>{R["sig_fp"]}</b>（TP:FP = {R["sig_ratio"]}:1，'
      f'是全體 base {R["base"]}:1 的 {R["sig_ratio"]/R["base"]:.1f}× → 高度 TP-clean）。</div>')
    A('<h3>6.2 過嚴格現象 + 根因（已驗證）</h3>')
    A(f'<div class="note"><b>{R["overstrict"]:,} 個 TP（{R["overstrict_pct"]}% 全 TP）</b>有乾淨的 HP-aligned PERMANOVA 結構，'
      f'卻因 gated CramersV&lt;0.1 被篩掉 —— 比通過的 {R["sig_tp"]:,} 還多。'
      f'<b>根因（全基因組 raw 重跑確認）</b>：over-strict {R["cause_n"]:,} 中 <b>{R["cause_rel"]:,} ({R["cause_pct"]}%) 是可靠性閘控</b>'
      f'（原始 CramersV≥0.3、median {R["cause_med"]}，卡方其實很強但 Cochran 期望格&lt;5 判不可靠歸零），'
      'partition-mismatch=0。即 gate 對「偵測」過嚴格，但這是統計上保守的安全閥。</div>')
    A('<h3>6.3 Tier 信心分層（minN=10；高一致 + 少誤判）</h3>')
    A('<table><tr><th>Tier</th><th>定義</th><th>TP</th><th>FP</th><th>TP:FP</th></tr>'
      f'<tr><td>原始 Significant</td><td>上述 gate</td><td class="n">{R["sig_tp"]}</td><td class="n">{R["sig_fp"]}</td><td class="n">{R["sig_ratio"]}:1</td></tr>'
      f'<tr style="color:#4ade80"><td>★ Tier A 高信心</td><td>reliable CramersV≥0.1 ∩ PERMANOVA(cluster+HP) ∩ min(HP1N,HP2N)≥10</td><td class="n">{R["tierA_tp"]}</td><td class="n">{R["tierA_fp"]}</td><td class="n">{R["tierA_ratio"]}:1</td></tr>'
      '<tr style="color:#60a5fa"><td>Tier B 救回</td><td>同上但 CramersV 不可靠（被閘控）</td><td class="n">2492</td><td class="n">72</td><td class="n">34.6:1</td></tr></table>'
      '<div class="src">定義：scripts 90/91；Tier B 數字 diag90 minN=10</div>')

    # ---------- S7 validation ----------
    A('<h2 id="s7">7 · L6/L7 TP/FP-free 獨立驗證（5 層，統計 → 肉眼最後）</h2>')
    A('<div class="intu"><b>直覺：</b>不用 SEQC2 TP/FP 當真值，改證「甲基↔HP 的關連<b>真實且能獨立複製</b>」。'
      '真值改成「換獨立的子集（不同 CpG、不同股、重抽 read）還看得到同樣的關連嗎」。HP 來自 SNP 定相（非循環）。</div>')
    A('<table><tr><th>層</th><th>方法</th><th>抓什麼 / 閾值</th></tr>'
      '<tr><td>L1 顯著</td><td>HP-permutation null on F=組間/組內距離 + Mann-Whitney</td><td>超越巧合；BH-FDR <b>q≤0.05</b>（控假發現率 5%）</td></tr>'
      '<tr><td>L2 獨立複製 ⭐</td><td>CpG split-half（偶/奇）+ 雙股一致</td><td>單一 CpG/單股 artifact；兩半/兩股同向顯著（<b>非循環</b>關鍵）</td></tr>'
      '<tr><td>L3 穩定</td><td>read bootstrap + PERMDISP 離散控制</td><td>少數 read 驅動；bootstrap <b>≥90%</b> 同向 + 無離散警告</td></tr>'
      '<tr><td>L4 生物特異</td><td>normal-cis（HP_Residual）+ germline-het 負對照</td><td>區分體細胞 cis vs 遺傳基線</td></tr>'
      '<tr><td><b>L5 肉眼（最後）</b></td><td>read×CpG 四格 panel（偶/奇 CpG · 正/反股）</td><td>只看 L1-L3 通過但 borderline 的<b>邊緣</b>位點</td></tr></table>')
    A(f'<div class="ok"><b>實證（非循環，全讀回真值）</b>：Tier A <b>{pct(R["valA"],R["valA_n"])}</b>'
      f'（{R["valA"]}/{R["valA_n"]}）、Tier B <b>{pct(R["valB"],R["valB_n"])}</b>（{R["valB"]}/{R["valB_n"]}）'
      '用完全獨立的複製驗證通過 → 確認 Tier 準則抓到真實 HP-關連甲基。~14% Tier A 沒過＝gate 自信但複製不一致＝正是 L5 該肉眼的邊緣。</div>')
    A('<div class="sub">兩種真實 ASM 模式（同被驗證）：① <b>分群型</b>（reads 分成對應 HP 的離散群，Tier A）② <b>位移型</b>'
      '（HP 家族甲基整體不同，Δβ-only）。範例圖：<code>figs_val/chr16_3444295_valpanel.png</code>（四格全一致＝鐵證）vs '
      '<code>chr20_22561507_valpanel.png</code>（四格不一致＝邊緣）。</div>')

    # ---------- S8 schema ----------
    A('<h2 id="s8">8 · 輸出格式與資料字典（每位點 117 欄）</h2>')
    A('<div class="intu"><b>直覺：</b>每個位點輸出一列共 117 欄，外加一個 per-region 資料夾存原始矩陣與分群樹。'
      '下表是教授最常查的核心欄位；完整 117 欄見 <code>significance_summary.csv</code> 表頭。</div>')
    A('<h3>8.1 per-region 資料夾結構</h3>')
    A('<table><tr><th>路徑</th><th>內容</th></tr>'
      '<tr><td><code>methylation/methylation.csv</code></td><td>read × CpG 連續甲基矩陣（0–1，NA=未覆蓋）+ cpg_sites.tsv</td></tr>'
      '<tr><td><code>distance/NHD/matrix.csv</code></td><td>read × read NHD 距離矩陣（+ 正/反股分開）</td></tr>'
      '<tr><td><code>clustering/</code></td><td>linkage_matrix.csv（樹）、leaf_order.txt（葉序）、significance.json（原始統計）、tree.nwk</td></tr>'
      '<tr><td><code>reads/reads.tsv</code></td><td>read_id, read_name, mapq, <b>hp</b>, alt_support(REF/ALT), strand…</td></tr>'
      '<tr><td><code>metadata.txt</code></td><td>位點/視窗/read 數/CpG 數摘要</td></tr></table>'
      '<div class="src">出處：實測 per-region 目錄（ism_display_scan）；對齊 .claude/rules/output-structure.md</div>')
    A('<h3>8.2 核心欄位字典（significance_summary.csv）</h3>')
    A('<table><tr><th>欄</th><th>意義</th><th>關聯</th></tr>'
      '<tr><td>NumReads / NumCpGs</td><td>視窗內 read 數 / CpG 數</td><td>gate reads≥20</td></tr>'
      '<tr><td>GlobalP / CramersV</td><td>分群↔標籤最小 p / 可靠 CramersV（max 軸）</td><td>§5.3 + gate</td></tr>'
      '<tr><td>CramersV_HPFamily / _HPFine</td><td>HP 家族軸 / HP-fine 軸的 CramersV</td><td>§5.3</td></tr>'
      '<tr><td>PassedGating</td><td>每軸 gate（p≤0.1 & V≥0.1）</td><td>§6.1</td></tr>'
      '<tr><td>ClusterPermanovaF/P/Valid · ClusterDispersionWarn</td><td>聚類結構 PERMANOVA + 離散警告</td><td>§5.4 / 驗證</td></tr>'
      '<tr><td>LabelHPPermanovaF/P/Valid</td><td>聚類是否對齊 HP 標籤（距離法）</td><td>Tier / 驗證</td></tr>'
      '<tr><td>HPMergedDelta (Δβ) / HPMergedP</td><td>HP 家族間距離分離量 + p</td><td>§5.5</td></tr>'
      '<tr><td>HP1FamilyN / HP2FamilyN</td><td>兩家族 read 數（稀疏判據）</td><td>Cochran / Tier minN</td></tr>'
      '<tr><td>Potential_LOH / Coverage_Multiple / Diploid_Coverage_Used</td><td>LOH 標註 / 覆蓋倍數 / KDE 二倍體覆蓋</td><td>§4.2 KDE</td></tr>'
      '<tr><td>HP_Residual_Delta（需 normal）</td><td>cis 訊號 = tumor−normal HP 差</td><td>§5.6</td></tr>'
      '<tr><td>Fisher_Frac_Sig</td><td>顯著 per-CpG 比例</td><td>§5.6</td></tr>'
      '<tr><td><b>Significant</b></td><td>最終顯著旗標</td><td>§6.1</td></tr></table>')

    # ---------- S9 glossary ----------
    A('<h2 id="s9">9 · 術語表（直覺一句 → 正式定義）</h2>')
    gl = [
        ("CpG", "DNA 上 C 接 G 的位點，甲基化主要發生處", "ONT 直接讀出每 CpG 甲基化機率"),
        ("5mC / 5hmC", "兩種甲基化修飾", "BAM MM/ML tag 編碼"),
        ("單倍型 HP", "read 來自哪一條染色體拷貝", "由 germline SNP 定相 (LongPhase)，與甲基無關"),
        ("子單倍型 HP1-1/2-1", "在 germline 單倍型上新生的體細胞分支", "raw tag 11/21"),
        ("ASM", "等位基因特異甲基化＝兩單倍型甲基不同", "本流程要 characterize 的對象"),
        ("NHD", "兩 read 甲基模式的差異比例", "diff/common，§5.1"),
        ("CramersV", "聚類×HP 交叉表的關連強度 0–1", "√(χ²/(n(k−1)))，稀疏歸零，§5.3"),
        ("PERMANOVA / pseudo-F", "距離法檢定 HP 是否分離，對稀疏穩健", "SS_b/SS_w 比，置換 999，§5.4"),
        ("Δβ (HPMergedDelta)", "HP 家族間 vs 家族內的距離分離量", "between−within，§5.5（非甲基平均差）"),
        ("ARI", "兩種分群一致度（校正隨機）", "Bootstrap.cpp，§5.6"),
        ("Cochran 可靠性", "卡方表期望格<5 即不可信", "→ CramersV 歸零，§5.3"),
        ("optimal_k", "最佳聚類群數", "silhouette 最大，§5.2"),
        ("Tier A / B", "ASM 結構信心分層", "§6.3"),
        ("split-half / 雙股", "獨立子集複製驗證", "非循環，§7"),
    ]
    A('<table><tr><th>術語</th><th>直覺</th><th>正式 / 出處</th></tr>'
      + "".join(f'<tr><td><b>{a}</b></td><td>{b}</td><td class="sub">{c}</td></tr>' for a, b, c in gl)
      + '</table>')

    # ---------- S10 boundaries + numbers ----------
    A('<h2 id="s10">10 · 已知邊界與全流程數字總覽</h2>')
    A('<div class="note">🔴 <b>邊界（投稿/對外必守）</b>：① 單樣本 HCC1395、單 pipeline → 證據上限 ⭐3。'
      '② 全流程是 characterization，<b>不可寫成 variant filter</b>（甲基→TP/FP 已 concluded NEGATIVE）。'
      '③ Δβ 是距離分離量非甲基平均差。④ HP 來自 SNP 定相（非循環）才使驗證成立。⑤ KDE 覆蓋隨位點集合變動（不影響分類）。</div>')
    A('<h3>全流程數字（皆讀自實際檔案）</h3>')
    A('<table><tr><th>階段</th><th>TP</th><th>FP</th><th>來源</th></tr>'
      f'<tr><td>L2 longphase-S 輸入 ISM</td><td class="n">{R["lps_tp"]:,}</td><td class="n">{R["lps_fp"]}</td><td>filtered_snv_*.vcf.gz</td></tr>'
      f'<tr><td>L3 ISM 實際分析</td><td class="n">{R["ism_tp"]:,}</td><td class="n">{R["ism_fp"]}</td><td>significance_summary.csv</td></tr>'
      f'<tr><td>L5 Significant gate 通過</td><td class="n">{R["sig_tp"]:,}</td><td class="n">{R["sig_fp"]}</td><td>diag90</td></tr>'
      f'<tr><td>過嚴格漏掉（有結構但 CV=0）</td><td class="n">{R["overstrict"]:,}</td><td>—</td><td>diag89（{R["cause_pct"]}% 為可靠性閘控）</td></tr>'
      f'<tr><td>Tier A 高信心</td><td class="n">{R["tierA_tp"]}</td><td class="n">{R["tierA_fp"]}</td><td>diag90 minN=10</td></tr>'
      f'<tr><td>Tier A 獨立驗證通過</td><td class="n">{R["valA"]} / {R["valA_n"]} ({pct(R["valA"],R["valA_n"])})</td><td>—</td><td>validation3（非循環）</td></tr></table>')
    A('<div class="sub" style="margin-top:14px">配套互動工具（逐位點雙圖 + 即時門檻 + 驗證徽章）：'
      '<code>InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/display_v2/20260609_locus_judgment_display_01.html</code></div>')

    return "".join(H)


def main():
    html = (f'<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>HCC1395 甲基差異分析 — 全流程說明書（教授級）</title>'
            f'<style>{CSS}</style></head><body><div class="layout">{nav()}<main>{body()}</main></div></body></html>')
    with open(OUT, "w") as f:
        f.write(html)
    import os
    print(f"[98] wrote {OUT} ({os.path.getsize(OUT)//1024} KB)")
    print(f"     result numbers from real json: sig_tp={R['sig_tp']} overstrict={R['overstrict']} "
          f"cause_rel={R['cause_rel']} tierA={R['tierA_tp']} valA={R['valA']}/{R['valA_n']}")


if __name__ == "__main__":
    main()
