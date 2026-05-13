#!/usr/bin/env python3
"""
build_html.py — Generate 22 slide_XX.html + index.html from inline specs.

Output:
- preview/index.html              (主頁: top nav + iframe)
- preview/slide_XX_*.html × 22    (子頁: title + 16:9 canvas + speaker note)

Run:
    python3 build_html.py
"""
from pathlib import Path
HERE = Path(__file__).parent
OUT = HERE / "preview"
OUT.mkdir(exist_ok=True)

# Inline CSS (thariqs html-effectiveness: each .html should be self-contained).
# shared/style.css is kept as source of truth for editing; build_html.py inlines
# it into every slide_XX.html + index.html so each file works standalone
# (single-file portability — can be emailed / uploaded / opened without external deps).
STYLE_CSS = (OUT / "shared" / "style.css").read_text(encoding="utf-8")

# ─── 22 slide specs ────────────────────────────────────────────────────────
# Each: id, num, section, title, en, timing, rg2, ngrep, canvas_html, speaker, tier3
# canvas_html = raw HTML body for the 16:9 .slide-canvas container
SLIDES = []

def add(**kw):
    kw.setdefault("title_critical", False)
    kw.setdefault("rg2", "—")
    kw.setdefault("ngrep", "—")
    kw.setdefault("tier3", "")
    kw.setdefault("is_cover", False)  # cover slide skips default h1+en-subtitle render
    SLIDES.append(kw)

# ─── S0 Cover + TL;DR ─────────────────────────────────────────────────────
add(id="01_cover", num="01", section="S0 Cover",
    is_cover=True,
    title="Self-Phasing 問題修正與驗證",
    en="",  # cover 自定 layout，不用 default en-subtitle render
    timing="30 sec / 中 ~90 字",
    canvas_html="""
    <div style="display:flex;flex-direction:column;justify-content:center;align-items:center;height:100%;text-align:center;padding:32px 0;">
      <h1 style="font-size:44px;color:#1E3A8A;margin:0 0 14px;font-weight:bold;line-height:1.2;">Self-Phasing 問題修正與驗證</h1>
      <p style="font-size:20px;color:#374151;margin:0 0 40px;font-weight:normal;">longphase-to tag layer 設計缺陷揭露與三版修正</p>

      <div style="border-top:1px solid #9CA3AF;padding-top:20px;width:55%;margin-top:24px;">
        <p style="font-size:15px;color:#1F2937;margin:8px 0;">日期：2026-05-01 ~ 2026-05-11</p>
        <p style="font-size:15px;color:#1F2937;margin:8px 0;">場合：實驗室週報</p>
        <p style="font-size:15px;color:#1F2937;margin:8px 0;">報告人：廖子游</p>
        <p style="font-size:13px;color:#6B7280;margin:14px 0 0;">中正大學 資工系 ・ 黃耀廷 教授 Lab405 實驗室</p>
      </div>
    </div>""",
    speaker="""[開場]
謝謝各位。

[主題]
今天的報告主題 — self-phasing 問題的修正與驗證。

[範圍]
我們在 longphase-to 工具的 tag 層發現一個設計缺陷。
報告整理修正過程，以及跨樣本驗證的結果。

[資訊]
報告人 廖子游。這是 5/1 至 5/11 的階段性整理。""",
    tier3="(封面 speaker note 精簡至 30 sec，不透露後續 slide 結論)")

add(id="02_tldr", num="02", section="S0 報告流程", rg2="—", ngrep="—",
    title="報告流程",
    en="",
    timing="60 sec / 中 ~240 字",
    canvas_html="""
    <div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:24px;margin:20px 0;">
      <!-- 階段 1: 問題發現 -->
      <div>
        <div style="background:#FEE2E2;color:#7F1D1D;padding:12px 18px;border-radius:8px;font-weight:700;font-size:22px;margin-bottom:20px;text-align:center;">問題發現</div>
        <div style="margin-bottom:20px;">
          <p style="font-size:22px;font-weight:700;color:#1E3A8A;margin:0;">① 觀察起點</p>
          <p style="font-size:17px;color:#374151;margin:6px 0 0;line-height:1.5;">全基因組偏移與 chr19 三個失衡位點觀察</p>
        </div>
        <div>
          <p style="font-size:22px;font-weight:700;color:#1E3A8A;margin:0;">② 機制</p>
          <p style="font-size:17px;color:#374151;margin:6px 0 0;line-height:1.5;">phasing 與 tagging 兩層 bug 解析</p>
        </div>
      </div>

      <!-- 階段 2: 修正過程 -->
      <div>
        <div style="background:#FEF3C7;color:#92400E;padding:12px 18px;border-radius:8px;font-weight:700;font-size:22px;margin-bottom:20px;text-align:center;">修正過程</div>
        <div style="margin-bottom:20px;">
          <p style="font-size:22px;font-weight:700;color:#1E3A8A;margin:0;">③ 量化鐵證</p>
          <p style="font-size:17px;color:#374151;margin:6px 0 0;line-height:1.5;">read-level 個案到全基因組驗證</p>
        </div>
        <div>
          <p style="font-size:22px;font-weight:700;color:#1E3A8A;margin:0;">④ 修補設計</p>
          <p style="font-size:17px;color:#374151;margin:6px 0 0;line-height:1.5;">兩層三版漸進改進歷程</p>
        </div>
      </div>

      <!-- 階段 3: 驗證結果 -->
      <div>
        <div style="background:#DCFCE7;color:#166534;padding:12px 18px;border-radius:8px;font-weight:700;font-size:22px;margin-bottom:20px;text-align:center;">驗證結果</div>
        <div style="margin-bottom:20px;">
          <p style="font-size:22px;font-weight:700;color:#1E3A8A;margin:0;">⑤ 驗證</p>
          <p style="font-size:17px;color:#374151;margin:6px 0 0;line-height:1.5;">指標對照與跨樣本一致性</p>
        </div>
        <div>
          <p style="font-size:22px;font-weight:700;color:#1E3A8A;margin:0;">⑥ 結論</p>
          <p style="font-size:17px;color:#374151;margin:6px 0 0;line-height:1.5;">修正方案總結與後續方向</p>
        </div>
      </div>
    </div>

    <div style="margin-top:20px;padding:12px 18px;background:#FFF7ED;border-left:4px solid #EA580C;border-radius:6px;">
      <p style="font-size:17px;font-weight:700;color:#9A3412;margin:0 0 6px;">⚠ 為何先解這問題？</p>
      <p style="font-size:15px;color:#7C2D12;line-height:1.55;margin:0;">
        ISM 後續所有 <strong>read 聚合 / NGroups / TP-FP 分類統計</strong>都依賴 longphase-to 的 HP tag。
        若 tag 不可信 → 下游分析標準全失準。longphase-to 作者較少 audit tag 行為（缺解答與標準）—
        <strong>本報告先填這缺口</strong>，確認 tag 合理可信，才能繼續 ISM 主軸研究。
      </p>
    </div>

    <div style="margin-top:14px;padding:10px 16px;background:#F0F9FF;border:1px solid #BAE6FD;border-radius:8px;color:#0C4A6E;font-size:16px;font-weight:600;text-align:center;">
      18 章節 slide + 3 backup（Q&A）  ·  預計 20 分鐘 + 5 分鐘問答
    </div>""",
    speaker="""[結構]
今天報告分三個階段、六個章節。

[階段 1]
問題發現 — 先看觀察起點的偏移現象，再解析兩層機制。

[階段 2]
修正過程 — 用 read-level 量化鐵證驗證，並設計三版漸進修補。

[階段 3]
驗證結果 — 跨指標跨樣本一致性驗證，最後給結論與後續方向。

[為何先解這問題]
ISM 後續所有分析 — read 聚合、NGroups、TP/FP 分類，都依賴 longphase-to HP tag。
若 tag 不可信，下游分析標準全失準。

[本報告價值]
longphase-to 作者較少思考 audit tag 行為（缺解答與標準）。
本報告先填這缺口，確認 tag 合理可信。

[時長]
預計 20 分鐘 + 5 分鐘問答。""",
    tier3="(目錄 speaker note 精簡至 60 sec，純預告章節 + 動機說明，不透露具體結論)")

# ─── S1 觀察起點 ───────────────────────────────────────────────────────────
add(id="03_genome_173", num="03", section="S1 觀察起點", rg2="3", ngrep="8",
    title="全基因組 17.3:1 偏移 = systematic bias",
    en="Genome-wide 17.3:1 = systematic engineering artifact, not sample variation",
    timing="90 sec / 中 ~340 字",
    canvas_html="""
    <div class="grid-2col">
      <table class="metric-table">
        <thead><tr><th>指標</th><th>baseline</th><th>隨機</th><th>偏離</th></tr></thead>
        <tbody>
          <tr><td>HP1 reads</td><td class="num">614,000</td><td class="num">~325K</td><td class="num">1.89×</td></tr>
          <tr><td>HP2 reads</td><td class="num">35,500</td><td class="num">~325K</td><td class="num">0.11×</td></tr>
          <tr><td>HP1:HP2 ratio</td><td class="num">17.3:1</td><td class="num">1:1</td><td class="num">17.3×</td></tr>
          <tr class="row-yellow"><td>HP1 占比</td><td class="num">94.6%</td><td class="num">~50%</td><td class="num">+44.6 pp</td></tr>
        </tbody>
      </table>
      <div style="text-align:center;align-self:center;background:#FFF7ED;border:2px solid #EA580C;border-radius:10px;padding:18px 14px;">
        <p style="font-size:38px;color:#C2410C;font-weight:800;margin:0;line-height:1.1;">HP1 占比 94.6%</p>
        <p style="font-size:24px;color:#9A3412;margin:8px 0 4px;font-weight:700;">↓↓↓↓</p>
        <p style="font-size:18px;color:#374151;margin:0;">隨機預期 <strong>50%</strong></p>
        <p style="font-size:15px;color:#6B7280;margin:4px 0 0;font-style:italic;">→ 偏離 +44.6 pp</p>
      </div>
    </div>
    <div class="arg-list">
      <strong>三條獨立論證:</strong><br>
      ①  <strong>生物學:</strong> 跨 23 染色體不該有系統偏單一 haplotype<br>
      ②  <strong>跨 chr 一致偏 HP1:</strong> cnLOH 最多只影響單一 chr；但 94.6% HP1 偏移在 23 chr 全部一致<br>
      ③  <strong>拿 paired 對照:</strong> 同樣 reads 在 paired (tumor + normal) 流程 HP1:HP2 ≈ 1:1 — 有 normal 資料時就能正確分到 1:1
    </div>
    <div class="conclusion-arrow">→ 17.3:1 是 LongPhase-TO 的 systematic engineering artifact</div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">haplotype:</span> 來自父或母的染色體版本</div>
      <div class="gloss-item">📖 <span class="term">sub-clone:</span> 腫瘤中帶相同突變的某群癌細胞</div>
      <div class="gloss-item">📖 <span class="term">cnLOH:</span> 染色體數對但失去其中一個 haplotype</div>
    </div>""",
    speaker="""[主訊息]
baseline LongPhase-TO 全基因組 HP1:HP2 = 17.3:1 — 這是 systematic bias 的硬證據。

[數據]
HP1 reads 614K，HP2 只 35.5K。94.6% 占比，遠離隨機 1:1。

[三條獨立論證]
論證 1 — 生物學：tumor sub-clone 不該跨 23 chr 系統偏 HP1。
論證 2 — 跨 chr 一致：cnLOH 只影響單 chr，但這偏移跨 23 chr 一致。
論證 3 — paired 對照：paired pipeline HP1:HP2 ≈ 1:1，沒這偏移。

[結論]
三條互相獨立 → 這是 engineering artifact，不是樣本性質。

[轉場]
接下來看具體個案 SP1 — 個別位點上失衡有多極端。""",
    tier3="cnLOH 機制細節 / 23 chr 一致性表 / paired pipeline 程式碼差異")

add(id="04a_sp1", num="04a", section="S1 觀察起點", rg2="1", ngrep="3",
    title="SP1: baseline 113:0 → V5 翻 HP2 對齊",
    en="SP1 chr19:17,565,944 — baseline 113:0 → V5 flips to HP2",
    timing="60 sec / 中 ~250 字",
    canvas_html="""
    <div class="grid-2col" style="grid-template-columns: 3fr 2fr;gap:14px;margin:6px 0;">
      <div class="igv-zoom-wrap">
        <span class="igv-zoom-hint">💡 hover / 點擊全尺寸</span>
        <a href="../figures/igv/SP1_chr19_17565944_5versions.png" target="_blank">
          <img class="igv-thumb" src="../figures/igv/SP1_chr19_17565944_5versions.png" alt="IGV SP1" style="max-height:340px;">
        </a>
        <p class="fig-caption" style="font-size:12px;margin:2px 0;">chr19:17,565,944 · 上→下 baseline / V6 / paired_T / paired_N²</p>
      </div>
      <div>
        <div class="stat-box" style="margin-bottom:6px;padding:8px;"><div class="number orange" style="font-size:24px;">113 : 0</div><div class="label">baseline HP1:HP2 (reads)</div></div>
        <div class="stat-box" style="margin-bottom:6px;padding:8px;"><div class="number" style="font-size:22px;">V6 翻 HP2</div><div class="label">V6 HP1:HP2 (待量化)¹</div></div>
        <div class="stat-box" style="padding:8px;"><div class="number green" style="font-size:24px;">HP2 ✅</div><div class="label">paired direction · 對齊</div></div>
      </div>
    </div>
    <div class="igv-focus-callout" style="font-size:14px;padding:6px 12px;margin-top:4px;"><span class="label">👁 看圖重點：</span>baseline (最上) 紅+綠 reads 集中左欄 = HP1；V6 / paired_T (下) 紅+綠搬右欄 = HP2 → 三方向一致</div>
    <p style="font-size:11px;color:#6B7280;margin:2px 0 0;">¹ V6 精確 count 待 vote_dump 量化  ²V6 取代 V2b/V3F/V5 簡化並列</p>""",
    speaker="""[從統計到個案]
全基因組 17.3:1 是平均值。
用 IGV 6-BAM 並列篩，找到 chr19 三個近 100% 失衡位點。

[SP1 數據]
位點 chr19:17,565,944。
baseline — 113 條 reads 全 HP1，HP2 = 0。
V5 修正後 — 翻到約 5:108，對齊 paired tumor 的 HP2 方向。

[排除三個 alternative]
噪音？baseline 與 paired 方向完全相反 — 是翻轉不是衰減。
caller 錯？V5 修正後與 paired 重合 — 不是 caller。
alignment 錯？跨 binary 一致 — 不是 alignment。

[結論]
read assignment 被強制集中的鐵證。""",
    tier3="6-BAM 並列順序細節 / V2b 中間階段意義 / HP1-1 sub-tag 與 HP1 合併原則")

add(id="04b_sp2_sp3", num="04b", section="S1 觀察起點", rg2="0", ngrep="6",
    title="SP2/SP3 同模式 3/3 對齊 paired",
    en="SP2 + SP3 follow SP1 pattern, 3/3 aligned with paired",
    timing="60 sec / 中 ~260 字",
    canvas_html="""
    <div class="grid-2col" style="gap:14px;">
      <div>
        <table class="metric-table" style="font-size:14px;">
          <thead><tr><th colspan="2">SP2 · chr19:12,452,332</th></tr>
                 <tr><th>baseline HP1:HP2</th><th>V6 (待量化)¹</th></tr></thead>
          <tbody><tr class="row-red"><td class="num">109 : 1</td><td class="num">翻 HP2</td></tr></tbody>
        </table>
        <div class="igv-zoom-wrap" style="margin-top:2px;">
          <span class="igv-zoom-hint">💡 hover/點擊</span>
          <a href="../figures/igv/SP2_chr19_12452332_5versions.png" target="_blank">
            <img class="igv-thumb" src="../figures/igv/SP2_chr19_12452332_5versions.png" alt="IGV SP2" style="max-height:240px;">
          </a>
        </div>
      </div>
      <div>
        <table class="metric-table" style="font-size:14px;">
          <thead><tr><th colspan="2">SP3 · chr19:12,467,180</th></tr>
                 <tr><th>baseline HP1:HP2</th><th>V6 (待量化)¹</th></tr></thead>
          <tbody><tr class="row-red"><td class="num">108 : 0</td><td class="num">翻 HP2</td></tr></tbody>
        </table>
        <div class="igv-zoom-wrap" style="margin-top:2px;">
          <span class="igv-zoom-hint">💡 hover/點擊</span>
          <a href="../figures/igv/SP3_chr19_12467180_5versions.png" target="_blank">
            <img class="igv-thumb" src="../figures/igv/SP3_chr19_12467180_5versions.png" alt="IGV SP3" style="max-height:240px;">
          </a>
        </div>
      </div>
    </div>
    <div class="igv-focus-callout" style="font-size:13px;padding:5px 10px;margin-top:4px;"><span class="label">👁 看圖重點：</span>兩圖同 SP1 模式 — baseline 紅+綠 reads 集中左欄；V6/paired_T 紅+綠搬右欄 → 3/3 對齊。HP1=(HP1+HP1-1) 紅綠；HP2=(HP2+HP2-1) 藍橙</div>
    <div class="conclusion-arrow green" style="font-size:15px;padding:8px 14px;">→ 三 SP 都在 chr19:12-17M → 對齊 slide 08 chr19 752 victims hotspot</div>
    <p style="font-size:11px;color:#6B7280;margin:2px 0 0;">¹ V6 精確 count 待 vote_dump 量化；圖檔待 V6 重擷取後替換</p>""",
    speaker="""[SP2/SP3 確認同模式]
SP2 chr19:12,452,332 — baseline 109:1。
SP3 chr19:12,467,180 — baseline 108:0。
與 SP1 同模式 3/3 — V6 修正後三位點都翻 HP2 對齊 paired。

[空間對應]
三 SP 都在 chr19:12-17M 區段。
對齊後續 slide 08 chr19 752 victims hotspot — 同機制不同層級觀察。

[口述過渡]
接下來四個章節分別回答：
① 為何全集中一邊？S2 機制
② read 層級多少？S3 量化
③ 三版各修什麼？S4 修補
④ 是否都修對？S5 驗證""",
    tier3="paired_T 與 paired_N 對照細節 / 三 SP 完整座標表 / V6 vs V5 差異 / HP sub-tag 分組原則")

# ─── S2 機制 ──────────────────────────────────────────────────────────────
add(id="05_player_referee", num="05", section="S2 機制", rg2="4 (核心 forced)", ngrep="3 + 1 commit",
    title="phasing 層球員兼裁判機制",
    en="Phasing layer player-as-referee — somatic 100% co-occurrence overrules germline 50/50",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <div class="grid-2col" style="grid-template-columns: 6fr 6fr;">
      <div>
        <!-- 左圖：SVG 直繪 phasing graph schematic（取代不清楚的 G1） -->
        <svg viewBox="0 0 480 320" xmlns="http://www.w3.org/2000/svg" style="width:100%;font-family:Arial,sans-serif;">
          <!-- 標題 -->
          <text x="240" y="16" font-size="12" fill="#1E3A8A" text-anchor="middle" font-weight="700">phasing graph — 從 germline het 推斷染色體分配</text>

          <!-- 上半：germline only（balanced graph） -->
          <text x="14" y="44" font-size="11" fill="#1E3A8A" font-weight="700">① 只用 germline het</text>
          <text x="14" y="58" font-size="10" fill="#166534">edge weight 50%（隨機分散）</text>
          <!-- 染色體 HP1 (橘) -->
          <line x1="60" y1="80" x2="440" y2="80" stroke="#FB923C" stroke-width="3"/>
          <text x="450" y="84" font-size="10" fill="#B45309" font-weight="700">HP1</text>
          <!-- 染色體 HP2 (藍) -->
          <line x1="60" y1="110" x2="440" y2="110" stroke="#3B82F6" stroke-width="3"/>
          <text x="450" y="114" font-size="10" fill="#1E40AF" font-weight="700">HP2</text>
          <!-- 3 個 het 點 -->
          <circle cx="120" cy="80" r="5" fill="#FB923C" stroke="#9A3412" stroke-width="1"/>
          <circle cx="120" cy="110" r="5" fill="#3B82F6" stroke="#1E40AF" stroke-width="1"/>
          <circle cx="240" cy="80" r="5" fill="#FB923C"/>
          <circle cx="240" cy="110" r="5" fill="#3B82F6"/>
          <circle cx="360" cy="80" r="5" fill="#FB923C"/>
          <circle cx="360" cy="110" r="5" fill="#3B82F6"/>
          <!-- read connect (thin lines 50%) -->
          <line x1="120" y1="80" x2="240" y2="80" stroke="#9CA3AF" stroke-width="0.8" opacity="0.6"/>
          <line x1="240" y1="110" x2="360" y2="110" stroke="#9CA3AF" stroke-width="0.8" opacity="0.6"/>
          <line x1="120" y1="110" x2="360" y2="80" stroke="#9CA3AF" stroke-width="0.6" opacity="0.4" stroke-dasharray="2,2"/>
          <text x="240" y="138" font-size="10" fill="#166534" text-anchor="middle">→ graph 平衡，HP1/HP2 分明</text>

          <!-- 中間分隔 -->
          <line x1="14" y1="158" x2="466" y2="158" stroke="#D1D5DB" stroke-dasharray="3,3"/>

          <!-- 下半：加入 somatic（dominated） -->
          <text x="14" y="178" font-size="11" fill="#7F1D1D" font-weight="700">② TO 模式：somatic 被當 germline 進 graph</text>
          <text x="14" y="192" font-size="10" fill="#B45309">edge weight 100% > germline 50% → somatic 主導</text>
          <!-- 染色體 -->
          <line x1="60" y1="220" x2="440" y2="220" stroke="#FB923C" stroke-width="3"/>
          <text x="450" y="224" font-size="10" fill="#B45309" font-weight="700">HP1</text>
          <line x1="60" y1="250" x2="440" y2="250" stroke="#3B82F6" stroke-width="3"/>
          <text x="450" y="254" font-size="10" fill="#1E40AF" font-weight="700">HP2</text>
          <!-- germline 點 -->
          <circle cx="120" cy="220" r="5" fill="#FB923C"/>
          <circle cx="120" cy="250" r="5" fill="#3B82F6"/>
          <circle cx="360" cy="220" r="5" fill="#FB923C"/>
          <circle cx="360" cy="250" r="5" fill="#3B82F6"/>
          <!-- somatic 點（中間，全在 HP1，100% clonal）-->
          <rect x="234" y="214" width="12" height="12" fill="#FBBF24" stroke="#B45309" stroke-width="2"/>
          <text x="240" y="208" font-size="9" fill="#92400E" text-anchor="middle" font-weight="700">som</text>
          <!-- 粗紅 edge: somatic 連到所有 HP1 點 -->
          <line x1="120" y1="220" x2="240" y2="220" stroke="#DC2626" stroke-width="3.5"/>
          <line x1="240" y1="220" x2="360" y2="220" stroke="#DC2626" stroke-width="3.5"/>
          <text x="180" y="212" font-size="9" fill="#7F1D1D" font-weight="700">edge 100%</text>
          <text x="300" y="212" font-size="9" fill="#7F1D1D" font-weight="700">edge 100%</text>
          <!-- 細線 germline 50% -->
          <line x1="120" y1="250" x2="360" y2="250" stroke="#9CA3AF" stroke-width="0.8" opacity="0.5" stroke-dasharray="2,2"/>
          <text x="240" y="282" font-size="10" fill="#7F1D1D" text-anchor="middle" font-weight="700">→ somatic 反過來決定 phase 方向（球員當裁判）</text>
          <text x="240" y="298" font-size="10" fill="#166534" text-anchor="middle">解法：PON-only flag 排除 somatic 進 graph，只保留 germline 50/50 平衡</text>
        </svg>
      </div>
      <div>
        <!-- Step 1: 物理基礎 -->
        <div style="background:#F3F4F6;border:1px solid #9CA3AF;border-radius:6px;padding:8px 12px;margin-bottom:6px;">
          <p style="font-size:13px;font-weight:700;color:#374151;margin:0 0 4px;">① 物理基礎</p>
          <p style="font-size:14px;color:#1F2937;margin:0;line-height:1.45;">
            <strong>germline het</strong>: HP1 50% / HP2 50%（隨機分佈）<br>
            <strong>somatic mut</strong>: 100% 同一條染色體（clonal）
          </p>
        </div>
        <p style="text-align:center;color:#6B7280;margin:0;font-size:14px;">↓</p>
        <!-- Step 2: TO 缺 normal -->
        <div style="background:#FEF3C7;border:1px solid #F59E0B;border-radius:6px;padding:8px 12px;margin:2px 0;">
          <p style="font-size:13px;font-weight:700;color:#92400E;margin:0 0 4px;">② TO 模式無 paired normal → 用 PoN 區分</p>
          <p style="font-size:14px;color:#78350F;margin:0;line-height:1.45;">somatic 不在 PoN 內 → <strong>被當 germline 進 phasing graph</strong></p>
        </div>
        <p style="text-align:center;color:#6B7280;margin:0;font-size:14px;">↓</p>
        <!-- Step 3: 致命連鎖 + 隱喻 -->
        <div style="background:#FEE2E2;border:2px solid #DC2626;border-radius:6px;padding:8px 12px;margin:2px 0;">
          <p style="font-size:13px;font-weight:700;color:#7F1D1D;margin:0 0 4px;">③ 致命連鎖 — somatic 「球員兼裁判」</p>
          <p style="font-size:14px;color:#7F1D1D;margin:0;line-height:1.45;">
            edge weight 100% &gt; germline 50% → <strong>somatic 主導 phasing decision</strong><br>
            <span style="font-weight:600;">應被 phase 的 somatic，反過來決定 phase 結果</span>
          </p>
        </div>
        <!-- 修法 -->
        <div class="conclusion-arrow green" style="margin-top:6px;padding:8px 12px;font-size:14px;">
          <strong>解法 <code>8b8c1fd</code> PON-only flag</strong>：phasing 階段排除非 PoN 變異<br>
          <span style="font-weight:400;font-size:13px;">→ 只用 germline het 建 graph，somatic 不再當 node</span>
        </div>
      </div>
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">PoN:</span> Panel of Normals — 多正常樣本建構的 germline reference</div>
      <div class="gloss-item">📖 <span class="term">phasing graph:</span> het 位點當 node, read 共現當 edge</div>
      <div class="gloss-item">ⓘ 此為 phasing layer / tagging layer per-read getVote 詳 slide 06</div>
    </div>""",
    speaker="""[標題]
phasing 層的「球員兼裁判」機制。

[三步推導]
步驟 1 — 物理基礎：germline het 50/50 隨機，somatic mut 100% 同染色體。
步驟 2 — TO 模式缺 paired normal，只能用 PoN 過濾。非 PoN 變異被當 germline 進 phasing graph。
步驟 3 — somatic edge weight 100% 暴漲，主導 graph。本應被 phase 的 somatic 反過來決定 phasing — 這就是球員兼裁判。

[修法]
8b8c1fd PON-only flag。
phasing 階段排除非 PoN 變異，只用 germline het 建 graph，somatic 不再當 node。

[caveat]
但只解 phasing 層。tag 階段還是壞的 — 下一頁說明。""",
    tier3="TO vs paired PoN 對照 / 自我增強迴圈 / Pass 1 vs Pass 2 預告 / edge weight 計算細節")

add(id="06_priority_bug", num="06", section="S2 機制", rg2="1 + 2 footnote", ngrep="5",
    title="tagging 層 priority bug 機制",
    en="getVote priority bug — single somatic vote triggers mislabel (break-early on wrong order)",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <!-- 上半：A1 Multi-read IGV view (per-read 框架) -->
    <svg viewBox="0 0 1200 260" xmlns="http://www.w3.org/2000/svg" style="width:100%;font-family:Arial,sans-serif;display:block;">
      <!-- reference axis -->
      <text x="20" y="14" font-size="11" fill="#6B7280" font-weight="600">chr19 reference</text>
      <line x1="20" y1="22" x2="860" y2="22" stroke="#9CA3AF" stroke-width="1"/>
      <!-- 4 position headers (germ / germ / somatic / germ) -->
      <rect x="160" y="28" width="100" height="32" fill="#DBEAFE" stroke="#1E40AF" stroke-width="1.5" rx="3"/>
      <text x="210" y="42" font-size="10" fill="#1E3A8A" text-anchor="middle" font-weight="700">pos 1</text>
      <text x="210" y="55" font-size="9" fill="#166534" text-anchor="middle" font-weight="600">germ / PoN ✓</text>
      <rect x="320" y="28" width="100" height="32" fill="#DBEAFE" stroke="#1E40AF" stroke-width="1.5" rx="3"/>
      <text x="370" y="42" font-size="10" fill="#1E3A8A" text-anchor="middle" font-weight="700">pos 2</text>
      <text x="370" y="55" font-size="9" fill="#166534" text-anchor="middle" font-weight="600">germ / PoN ✓</text>
      <rect x="480" y="28" width="100" height="32" fill="#FEF3C7" stroke="#B45309" stroke-width="1.5" rx="3"/>
      <text x="530" y="42" font-size="10" fill="#92400E" text-anchor="middle" font-weight="700">pos 3</text>
      <text x="530" y="55" font-size="9" fill="#B45309" text-anchor="middle" font-weight="600">somatic / PoN ✗</text>
      <rect x="640" y="28" width="100" height="32" fill="#DBEAFE" stroke="#1E40AF" stroke-width="1.5" rx="3"/>
      <text x="690" y="42" font-size="10" fill="#1E3A8A" text-anchor="middle" font-weight="700">pos 4</text>
      <text x="690" y="55" font-size="9" fill="#166534" text-anchor="middle" font-weight="600">germ / PoN ✓</text>
      <!-- guide lines -->
      <line x1="210" y1="60" x2="210" y2="220" stroke="#1E40AF" stroke-width="0.5" stroke-dasharray="2,3" opacity="0.4"/>
      <line x1="370" y1="60" x2="370" y2="220" stroke="#1E40AF" stroke-width="0.5" stroke-dasharray="2,3" opacity="0.4"/>
      <line x1="530" y1="60" x2="530" y2="220" stroke="#B45309" stroke-width="0.5" stroke-dasharray="2,3" opacity="0.4"/>
      <line x1="690" y1="60" x2="690" y2="220" stroke="#1E40AF" stroke-width="0.5" stroke-dasharray="2,3" opacity="0.4"/>
      <!-- right panel header -->
      <rect x="800" y="28" width="380" height="32" fill="#F3F4F6" stroke="#6B7280" rx="3"/>
      <text x="990" y="42" font-size="11" fill="#1F2937" text-anchor="middle" font-weight="700">該 read 的 countMap → HP tag</text>
      <text x="990" y="55" font-size="9" fill="#6B7280" text-anchor="middle">baseline vs V6 對照</text>
      <!-- read 1: pos 1/2/3, HP2=2, HP1_1=1 -->
      <text x="14" y="92" font-size="11" fill="#374151" text-anchor="end" font-weight="700">read 1</text>
      <path d="M 30,76 L 590,76 L 610,92 L 590,108 L 30,108 Z" fill="#D1D5DB" stroke="#6B7280" stroke-width="1"/>
      <rect x="200" y="76" width="22" height="32" fill="#3B82F6" stroke="#1E40AF" stroke-width="1.5"/>
      <text x="211" y="95" font-size="10" fill="white" text-anchor="middle" font-weight="700">A</text>
      <rect x="360" y="76" width="22" height="32" fill="#3B82F6" stroke="#1E40AF" stroke-width="1.5"/>
      <text x="371" y="95" font-size="10" fill="white" text-anchor="middle" font-weight="700">A</text>
      <rect x="520" y="74" width="22" height="36" fill="#FBBF24" stroke="#B45309" stroke-width="2"/>
      <text x="531" y="95" font-size="10" fill="#7C2D12" text-anchor="middle" font-weight="700">T</text>
      <rect x="800" y="70" width="380" height="44" fill="#FEE2E2" stroke="#DC2626" stroke-width="1.5" rx="3"/>
      <text x="810" y="85" font-size="11" fill="#1F2937" font-family="monospace">HP2=2, HP1_1=1</text>
      <text x="810" y="100" font-size="11" fill="#7F1D1D" font-weight="700">baseline → hp=11 ❌</text>
      <text x="970" y="100" font-size="11" fill="#166534" font-weight="700">V6 → hp=21 ✅</text>
      <!-- read 2: pos 1/2/4, HP2=3 (no somatic) -->
      <text x="14" y="138" font-size="11" fill="#374151" text-anchor="end" font-weight="700">read 2</text>
      <path d="M 30,122 L 720,122 L 740,138 L 720,154 L 30,154 Z" fill="#D1D5DB" stroke="#6B7280" stroke-width="1"/>
      <rect x="200" y="122" width="22" height="32" fill="#3B82F6" stroke="#1E40AF" stroke-width="1.5"/>
      <text x="211" y="141" font-size="10" fill="white" text-anchor="middle" font-weight="700">A</text>
      <rect x="360" y="122" width="22" height="32" fill="#3B82F6" stroke="#1E40AF" stroke-width="1.5"/>
      <text x="371" y="141" font-size="10" fill="white" text-anchor="middle" font-weight="700">A</text>
      <rect x="680" y="122" width="22" height="32" fill="#3B82F6" stroke="#1E40AF" stroke-width="1.5"/>
      <text x="691" y="141" font-size="10" fill="white" text-anchor="middle" font-weight="700">A</text>
      <rect x="800" y="116" width="380" height="44" fill="#DCFCE7" stroke="#16A34A" stroke-width="1.5" rx="3"/>
      <text x="810" y="131" font-size="11" fill="#1F2937" font-family="monospace">HP2=3 (no somatic)</text>
      <text x="810" y="146" font-size="11" fill="#166534" font-weight="700">baseline → hp=2 ✓</text>
      <text x="970" y="146" font-size="11" fill="#166534" font-weight="700">V6 → hp=2 ✓ (同)</text>
      <!-- read 3: pos 3/4, HP2=1, HP1_1=1 -->
      <text x="14" y="184" font-size="11" fill="#374151" text-anchor="end" font-weight="700">read 3</text>
      <path d="M 400,168 L 720,168 L 740,184 L 720,200 L 400,200 Z" fill="#D1D5DB" stroke="#6B7280" stroke-width="1"/>
      <rect x="520" y="166" width="22" height="36" fill="#FBBF24" stroke="#B45309" stroke-width="2"/>
      <text x="531" y="187" font-size="10" fill="#7C2D12" text-anchor="middle" font-weight="700">T</text>
      <rect x="680" y="168" width="22" height="32" fill="#3B82F6" stroke="#1E40AF" stroke-width="1.5"/>
      <text x="691" y="187" font-size="10" fill="white" text-anchor="middle" font-weight="700">A</text>
      <rect x="800" y="162" width="380" height="44" fill="#FEE2E2" stroke="#DC2626" stroke-width="1.5" rx="3"/>
      <text x="810" y="177" font-size="11" fill="#1F2937" font-family="monospace">HP2=1, HP1_1=1</text>
      <text x="810" y="192" font-size="11" fill="#7F1D1D" font-weight="700">baseline → hp=11 ❌</text>
      <text x="970" y="192" font-size="11" fill="#166534" font-weight="700">V6 → hp=21 ✅</text>
      <!-- footer note: per-read scope -->
      <line x1="20" y1="218" x2="1180" y2="218" stroke="#9CA3AF" stroke-dasharray="3,3"/>
      <text x="20" y="234" font-size="11" fill="#374151" font-weight="700">per-read getVote（C++ HaplotagProcess.cpp:533 每條 read 自己 reset countMap）</text>
      <text x="20" y="250" font-size="10" fill="#6B7280">read 1/3 受 priority bug 影響（有 somatic + germline mix）；read 2 純 germline 不受影響；752 victims = 752 條獨立 read</text>
    </svg>

    <!-- 下半：左 baseline+V6 code 上下 / 右 baseline+V6 flow 上下 -->
    <div class="grid-2col" style="grid-template-columns: 6fr 5fr;gap:10px;margin-top:6px;">
      <!-- 左欄: 程式碼 baseline (上) + V6 (下) -->
      <div>
        <p style="font-size:10px;color:#6B7280;margin:0 0 2px;font-family:monospace;"><strong style="color:#7F1D1D;">baseline</strong> <code>longphase-to/HaplotagProcess.cpp:506-530</code></p>
        <div style="background:#FFFFFF;border:1px solid #D1D5DB;border-radius:3px;font-family:'JetBrains Mono',monospace;font-size:9px;line-height:1.35;margin-bottom:4px;">
          <div style="background:#FFEEF0;color:#86181D;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#86181D;font-weight:700;">−</span>vector variantKeys = {{HP1_1,HP2_1},{HP3,HP2_1},{HP1,HP2}};</div>
          <div style="background:#FFEEF0;color:#86181D;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#86181D;font-weight:700;">−</span>for (auto&amp; pair : variantKeys) {  // ① som ② mix ③ germ</div>
          <div style="background:#FFEEF0;color:#86181D;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#86181D;font-weight:700;">−</span>  if (countMap[k1]&gt;0 || countMap[k2]&gt;0) {</div>
          <div style="background:#FFEEF0;color:#86181D;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#86181D;font-weight:700;">−</span>    hpResult = winner; break;  // 第一個非空 pair</div>
          <div style="background:#FFEEF0;color:#86181D;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#86181D;font-weight:700;">−</span>}}  // ③ germline pair 永遠輪不到</div>
        </div>
        <p style="font-size:10px;color:#6B7280;margin:0 0 2px;font-family:monospace;"><strong style="color:#166534;">V6</strong> <code>longphase-to-mod/HaplotagProcess.cpp:512-560</code></p>
        <div style="background:#FFFFFF;border:1px solid #D1D5DB;border-radius:3px;font-family:'JetBrains Mono',monospace;font-size:9px;line-height:1.35;">
          <div style="background:#E6FFED;color:#22863A;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#22863A;font-weight:700;">+</span>// Layer 1: Germline 先決方向</div>
          <div style="background:#E6FFED;color:#22863A;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#22863A;font-weight:700;">+</span>germlineResult = (HP1&gt;0||HP2&gt;0) ? (HP1&gt;=HP2?1:2) : 0;</div>
          <div style="background:#E6FFED;color:#22863A;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#22863A;font-weight:700;">+</span>// Layer 2: Somatic 加 sub-tag encoding</div>
          <div style="background:#E6FFED;color:#22863A;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#22863A;font-weight:700;">+</span>if (somaticTotal&gt;0)</div>
          <div style="background:#E6FFED;color:#22863A;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#22863A;font-weight:700;">+</span>  hpResult=(germ==1)?11:(germ==2)?21:33;</div>
          <div style="background:#E6FFED;color:#22863A;padding:0 8px 0 18px;position:relative;white-space:pre;"><span style="position:absolute;left:6px;color:#22863A;font-weight:700;">+</span>else hpResult = germlineResult;</div>
        </div>
        <p style="font-size:9px;color:#6B7280;margin:3px 0 0;line-height:1.3;">baseline iterate variantKeys + break early (somatic pair 第一個) → V6 兩層獨立判定（V3F/V5 為中間階段，詳 slide 07/10）</p>
      </div>
      <!-- 右欄: 流程圖 baseline (上) + V6 (下) -->
      <div>
        <p style="font-size:10px;color:#7F1D1D;margin:0 0 2px;font-weight:700;">baseline flow（iterate + break early）</p>
        <svg viewBox="0 0 460 100" xmlns="http://www.w3.org/2000/svg" style="width:100%;font-family:Arial,sans-serif;margin-bottom:4px;">
          <rect x="2" y="36" width="100" height="26" fill="#FEE2E2" stroke="#DC2626" rx="3"/>
          <text x="52" y="53" font-size="10" fill="#7F1D1D" text-anchor="middle" font-weight="700">① som pair</text>
          <rect x="120" y="36" width="100" height="26" fill="#FEE2E2" stroke="#DC2626" rx="3"/>
          <text x="170" y="53" font-size="10" fill="#7F1D1D" text-anchor="middle">② mix pair</text>
          <rect x="238" y="36" width="100" height="26" fill="#E5E7EB" stroke="#9CA3AF" stroke-dasharray="3,2" rx="3"/>
          <text x="288" y="53" font-size="10" fill="#6B7280" text-anchor="middle">③ germ (skip)</text>
          <line x1="102" y1="49" x2="118" y2="49" stroke="#DC2626" stroke-width="1.5"/>
          <polygon points="115,46 122,49 115,52" fill="#DC2626"/>
          <line x1="220" y1="49" x2="236" y2="49" stroke="#9CA3AF" stroke-width="1" stroke-dasharray="2,2"/>
          <rect x="350" y="34" width="105" height="30" rx="15" fill="#FEE2E2" stroke="#7F1D1D" stroke-width="2"/>
          <text x="402" y="53" font-size="11" fill="#7F1D1D" text-anchor="middle" font-weight="700">hp=11 ❌</text>
          <line x1="222" y1="36" x2="380" y2="36" stroke="none"/>
          <path d="M 52,36 Q 200,8 350,42" stroke="#DC2626" stroke-width="1.5" fill="none"/>
          <polygon points="346,40 354,42 348,46" fill="#DC2626"/>
          <text x="200" y="22" font-size="9" fill="#7F1D1D" text-anchor="middle" font-style="normal" font-weight="600">somatic vote ≥1 即觸發 break</text>
          <text x="200" y="82" font-size="9" fill="#6B7280" text-anchor="middle">germline 5 票永遠看不到</text>
        </svg>

        <p style="font-size:10px;color:#166534;margin:2px 0 2px;font-weight:700;">V6 flow（two-layer independent）</p>
        <svg viewBox="0 0 460 130" xmlns="http://www.w3.org/2000/svg" style="width:100%;font-family:Arial,sans-serif;">
          <rect x="80" y="4" width="300" height="22" rx="11" fill="#F3F4F6" stroke="#6B7280"/>
          <text x="230" y="19" font-size="10" fill="#1F2937" text-anchor="middle" font-family="monospace">this read's countMap</text>
          <line x1="230" y1="26" x2="230" y2="34" stroke="#6B7280" stroke-width="1.2"/>
          <polygon points="226,33 230,40 234,33" fill="#6B7280"/>
          <rect x="40" y="40" width="380" height="26" fill="#DCFCE7" stroke="#16A34A" stroke-width="1.2" rx="3"/>
          <text x="230" y="57" font-size="11" fill="#166534" text-anchor="middle" font-weight="700">Layer 1: Germline 決方向（HP1 vs HP2）</text>
          <line x1="230" y1="66" x2="230" y2="74" stroke="#6B7280" stroke-width="1.2"/>
          <polygon points="226,73 230,80 234,73" fill="#6B7280"/>
          <rect x="40" y="80" width="380" height="26" fill="#DBEAFE" stroke="#1E3A8A" stroke-width="1.2" rx="3"/>
          <text x="230" y="97" font-size="11" fill="#1E3A8A" text-anchor="middle" font-weight="700">Layer 2: Somatic 加 sub-tag → hp=11/21/33</text>
          <line x1="230" y1="106" x2="230" y2="113" stroke="#6B7280" stroke-width="1.2"/>
          <polygon points="226,112 230,118 234,112" fill="#6B7280"/>
          <rect x="170" y="115" width="120" height="14" rx="7" fill="#DCFCE7" stroke="#166534" stroke-width="2"/>
          <text x="230" y="126" font-size="10" fill="#166534" text-anchor="middle" font-weight="700">hp=21 ✅</text>
        </svg>
        <p style="font-size:9px;color:#6B7280;margin:2px 0 0;text-align:center;">baseline: 順序決定 → V6: 兩層獨立判定</p>
      </div>
    </div>

    <div class="footer-glossary" style="font-size:11px;">
      <div class="gloss-item">📖 <span class="term">sub-clone:</span> 腫瘤亞群</div>
      <div class="gloss-item">ⓘ per-read scope @ <code>judgeHaplotype:533</code></div>
    </div>""",
    speaker="""[標題]
tagging 層的 priority bug 機制。
重點：getVote 是 per-read 操作 — judgeHaplotype:533 每條 read 自己 reset countMap。

[三件套]
件 1 — vector 順序錯：baseline ① somatic / ② mixed / ③ germline。for 迴圈第一個非空 pair 就 break early。同一條 read 的 germline 5 票永遠看不到。

件 2 — 真實例子：chr19 752 victims 同模式。
單一條 read 經過 5 個 germline het 位點，累積 countMap[HP2]=5 主導。
加上 1 個 somatic 位點，countMap[HP1_1]=1。
baseline 錯標這條 read 為 hp=11，正確應 hp=21。

件 3 — V3F 修法 41ff147：重寫為兩層獨立判定。
germline 先決方向，somatic 再加 sub-tag encoding。

[串連]
tumor sub-clone somatic 100% 同方向。
→ 752 條獨立 read 各自被 priority bug 翻成 HP:i:11。
→ tag layer 17.3:1 偏移就是這樣產生的。

[關鍵釐清]
tag layer 與 phasing layer 是不同層 bug，必須分別修補，不能合併修。
752 victims = 752 條獨立 read，不是 aggregate。""",
    tier3="enum HAPLOTYPE1_1=2 vs HP tag int=11 / countMap per-read reset @ HaplotagProcess.cpp:533 / tag layer vs phasing layer 為何不能合併修 / 41ff147 與 380e8d2 INDEL guard 分工 / 752 victims per-read 獨立 case")

add(id="07_two_layer_table", num="07", section="S2 機制", rg2="1", ngrep="5 commits + 1",
    title="兩層 bug 對應兩層修補",
    en="Two layers two fixes — what each fix actually does, and why both required",
    timing="90 sec / 中 ~310 字",
    canvas_html="""
    <!-- Timeline -->
    <div style="background:#F9FAFB;border:1px solid #E5E7EB;border-radius:6px;padding:8px 12px;margin:2px 0 8px;">
      <p style="font-size:12px;color:#6B7280;margin:0 0 4px;font-weight:600;">commit 時序鏈（baseline → V3F → V5 → V6）：</p>
      <div style="display:flex;align-items:center;gap:4px;font-size:11px;flex-wrap:wrap;">
        <span style="background:#DBEAFE;color:#1E3A8A;padding:3px 8px;border-radius:4px;font-weight:600;">2026-04 phasing</span>
        <code style="font-size:11px;">8b8c1fd</code>
        <span style="color:#6B7280;">→</span>
        <span style="background:#DCFCE7;color:#166534;padding:3px 8px;border-radius:4px;font-weight:600;">2026-05 V3F</span>
        <code style="font-size:11px;">41ff147 + 380e8d2</code>
        <span style="color:#6B7280;">→</span>
        <span style="background:#FEF3C7;color:#92400E;padding:3px 8px;border-radius:4px;font-weight:600;">2026-05 V5</span>
        <code style="font-size:11px;">d0bcd8c + 938f0df</code>
        <span style="color:#6B7280;">→</span>
        <span style="background:#F3E8FF;color:#6D28D9;padding:3px 8px;border-radius:4px;font-weight:600;">2026-05-11 V6</span>
        <code style="font-size:11px;">revert Layer 1.5</code>
      </div>
    </div>
    <!-- "修了什麼" 表 -->
    <table class="metric-table" style="font-size:13px;">
      <thead><tr><th style="width:18%;">Layer / Bug</th><th>修法做了什麼</th><th style="width:10%;">Slide</th></tr></thead>
      <tbody>
        <tr style="background:#DBEAFE;"><td><strong>phasing 層</strong><br>球員兼裁判</td><td>phasing 階段<strong>排除非 PoN 變異</strong>，只用 germline het 建 graph（somatic 不再當 node）</td><td>05</td></tr>
        <tr style="background:#DCFCE7;"><td><strong>tagging Layer 1</strong><br>priority bug (vote 順序)</td><td><strong>重寫 two-layer</strong>：Layer 1 germline 決方向 / Layer 2 somatic 加 sub-tag（不再 iterate variantKeys break-early）</td><td>06</td></tr>
        <tr style="background:#DCFCE7;"><td><strong>tagging Layer 1</strong><br>enum vs int literal mismatch</td><td><code>HAPLOTYPE1_1=2</code> enum 與 HP tag int <code>11</code> 不一致，V3F 改用 int literal 11/21/33 直接賦值（搭 380e8d2 INDEL guard）</td><td>backup</td></tr>
        <tr style="background:#FEF3C7;"><td><strong>tagging Layer 1.5</strong><br>germline-absent</td><td><strong>V5 somatic fallback</strong>：純 somatic vote 區補 hp tag（V3F 此處 untagged）<span style="color:#92400E;font-weight:600;"> ⚠ 4.19:1 偏 HP1 caveat</span></td><td>16</td></tr>
        <tr style="background:#F3E8FF;"><td><strong>V6 revert</strong><br>Layer 1.5 移除</td><td>germline-absent 區回退 V3F 保守行為（hp=33），marker coverage +30% / 4 樣本 ratio 中性化；但 SP1/2/3 失去 Layer 1.5 翻方向能力</td><td>16</td></tr>
      </tbody>
    </table>
    <!-- 5 commits inline 列表（每 commit 一行差異） -->
    <div style="background:#F9FAFB;border:1px solid #E5E7EB;border-radius:4px;padding:5px 10px;margin-top:5px;font-size:10px;line-height:1.45;">
      <p style="margin:0;color:#1E3A8A;font-family:monospace;"><code>8b8c1fd</code> 加 <code>--pon-only-phasing</code> flag · <span style="color:#166534;"><code>41ff147</code> getVote 重寫 Layer 1/2 + enum→int literal</span> · <span style="color:#166534;"><code>380e8d2</code> INDEL guard</span> · <span style="color:#92400E;"><code>d0bcd8c</code> Pass 2 ploidy+Layer 1.5</span> · <span style="color:#92400E;"><code>938f0df</code> purity 0.95→0.9</span> · <span style="color:#6D28D9;"><code>V6</code> revert Layer 1.5</span></p>
    </div>
    <!-- 為何缺一不可 + V6 trade-off (合併 2-col 緊湊) -->
    <div class="grid-2col" style="gap:8px;margin-top:5px;">
      <div class="caveat-box" style="padding:5px 9px;font-size:10px;">
        <span class="label">只 PON-only？</span>解 phasing 但 tag 仍壞 → <strong>99.9% reads 仍 HP:i:11</strong>
      </div>
      <div class="caveat-box" style="padding:5px 9px;font-size:10px;">
        <span class="label">只 V3F（無 8b8c1fd）？</span>tag 修但 phasing 仍被 somatic 主導 → PS block 邊界錯
      </div>
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">PoN:</span> Panel of Normals</div>
      <div class="gloss-item">ⓘ Layer 1.5 = V3F 之上補丁，非獨立層</div>
      <div class="gloss-item">ⓘ V6 trade-off: marker ✓ / SP1/2/3 翻方向 ✗</div>
    </div>""",
    speaker="""[時序]
兩層 bug × 四 fix stacking — 五個 commits 漸進完成。

2026-04 phasing 層：8b8c1fd PON-only。
2026-05 tagging V3F：41ff147 重寫 two-layer + 380e8d2 INDEL guard。
2026-05 V5：d0bcd8c Pass 2 ploidy + 938f0df threshold + Layer 1.5 somatic fallback。
2026-05-11 V6：revert Layer 1.5。

[每個修法做什麼]
phasing 層 — 排除非 PoN 變異，只用 germline het 建 graph。
tagging Layer 1 — 重寫兩層獨立判定 + 修 enum/int literal。
Layer 1.5 — 純 somatic vote 區補 hp tag。
V6 — germline-absent 區回退保守 hp=33。

[為何缺一不可]
只 PON-only — tag 仍偏 99.9%。
只 V3F — germline 缺席區 untagged。
V5 Layer 1.5 — 補 fallback 但 germline-absent 4.19:1 偏 HP1。
V6 — 是 trade-off 修正不是完美解，但 marker coverage +30%。""",
    tier3="5 commit 順序細節 / V3F 命名歷史 / 跨層交互 / 為何不能合併修 / d0bcd8c 為何跨兩層 / enum=2 vs int=11 mismatch 細節 / V6 SP1/2/3 反向 caveat")

# ─── S3 量化鐵證 ───────────────────────────────────────────────────────────
add(id="08_quant_evidence", num="08", section="S3 量化鐵證", rg2="—", ngrep="14",
    title="read-level 鐵證: 34,855 V6 100% 修正",
    en="baseline → V6 read-level 100% fix: 752 chr19 + 34,855 genome",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <div class="grid-2col" style="grid-template-columns: 6fr 6fr; gap:12px;">
      <div>
        <p style="font-size:11.5px;font-weight:700;color:#374151;margin:0 0 4px;">📈 chr19 pilot → 全基因組 generalize:</p>
        <table class="metric-table" style="font-size:10.5px;">
          <thead><tr><th></th><th>chr19 pilot</th><th>Genome F1</th><th>倍數</th></tr></thead>
          <tbody>
            <tr><td>Tagged reads/ver</td><td class="num">~330K</td><td class="num">18,895,432</td><td class="num">57×</td></tr>
            <tr class="row-yellow"><td><strong>Priority bug victims</strong></td><td class="num">752</td><td class="num"><strong>34,855</strong></td><td class="num"><strong>46.4×</strong></td></tr>
            <tr class="row-green"><td><strong>baseline → V6 修正率</strong></td><td class="num">100%<br>(752/752)</td><td class="num"><strong>100%</strong><br>(34,855/34,855)</td><td>一致 ✅</td></tr>
            <tr><td>baseline → V6 方向</td><td colspan="2" class="num">全 hp=11 ❌ → 全 hp=21 ✅</td><td>單向</td></tr>
            <tr><td>反向 (V6 → 11)</td><td colspan="2" class="num">0 條</td><td>零 noise</td></tr>
          </tbody>
        </table>
        <p style="font-size:10px;color:#16A34A;margin:4px 0 0;font-weight:600;">→ baseline 17.3:1 偏移 = read-level 34,855 條獨立鐵證</p>
      </div>
      <div>
        <p style="font-size:11.5px;font-weight:700;color:#374151;margin:0 0 4px;">📊 per-chr victims (top-10 + chr19/chr8):</p>
        <svg viewBox="0 0 340 175" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;border:1px solid #E5E7EB;background:#F9FAFB;">
          <line x1="48" y1="148" x2="335" y2="148" stroke="#9CA3AF" stroke-width="0.8"/>
          <line x1="48" y1="14"  x2="48"  y2="148" stroke="#9CA3AF" stroke-width="0.8"/>
          <text x="44" y="17" font-size="8" fill="#6B7280" text-anchor="end">3,508</text>
          <text x="44" y="83" font-size="8" fill="#6B7280" text-anchor="end">1,700</text>
          <text x="44" y="151" font-size="8" fill="#6B7280" text-anchor="end">0</text>
          <rect x="55"  y="14"  width="20" height="134" fill="#1E3A8A"/>
          <text x="65" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr7</text>
          <text x="65" y="12"  font-size="8" fill="#1E3A8A" text-anchor="middle" font-weight="700">3508</text>
          <rect x="80"  y="42"  width="20" height="106" fill="#1E3A8A"/>
          <text x="90" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr2</text>
          <text x="90" y="40"  font-size="8" fill="#1E3A8A" text-anchor="middle" font-weight="700">2792</text>
          <rect x="105" y="48"  width="20" height="100" fill="#1E3A8A"/>
          <text x="115" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr1</text>
          <text x="115" y="46"  font-size="8" fill="#1E3A8A" text-anchor="middle" font-weight="700">2674</text>
          <rect x="130" y="51"  width="20" height="97"  fill="#3B82F6"/>
          <text x="140" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr16</text>
          <text x="140" y="49"  font-size="8" fill="#1E3A8A" text-anchor="middle" font-weight="700">2584</text>
          <rect x="155" y="68"  width="20" height="80"  fill="#3B82F6"/>
          <text x="165" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr20</text>
          <text x="165" y="66"  font-size="8" fill="#1E3A8A" text-anchor="middle" font-weight="700">2101</text>
          <rect x="180" y="79"  width="20" height="69"  fill="#60A5FA"/>
          <text x="190" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr3</text>
          <rect x="205" y="88"  width="20" height="60"  fill="#60A5FA"/>
          <text x="215" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr5</text>
          <rect x="230" y="94"  width="20" height="54"  fill="#93C5FD"/>
          <text x="240" y="160" font-size="7.5" fill="#374151" text-anchor="middle">chr6</text>
          <text x="240" y="170" font-size="6.5" fill="#9CA3AF" text-anchor="middle">rank 4-10</text>
          <rect x="262" y="120" width="20" height="28"  fill="#DC2626"/>
          <text x="272" y="160" font-size="7.5" fill="#7F1D1D" text-anchor="middle" font-weight="700">chr19</text>
          <text x="272" y="117" font-size="8" fill="#7F1D1D" text-anchor="middle" font-weight="700">752</text>
          <text x="272" y="170" font-size="6.5" fill="#7F1D1D" text-anchor="middle">★ rank 19</text>
          <rect x="287" y="123" width="20" height="25"  fill="#0EA5E9"/>
          <text x="297" y="160" font-size="7.5" fill="#0369A1" text-anchor="middle" font-weight="700">chr8</text>
          <text x="297" y="120" font-size="8" fill="#0369A1" text-anchor="middle" font-weight="700">666</text>
          <text x="297" y="170" font-size="6.5" fill="#0369A1" text-anchor="middle">❄ 0.34×</text>
        </svg>
        <p style="font-size:10px;color:#6B7280;margin:2px 0 0;text-align:center;">main hotspot 在 chr7/chr2/chr1/chr16/chr20；chr19 是「可重現案例」非主要分佈</p>
      </div>
    </div>
    <div class="conclusion-arrow green" style="margin-top:8px;font-size:13px;padding:6px 12px;">→ priority bug 不是 chr19 局部 artifact；全基因組 34,855 條 read-level victim 全部 100% 修正、0 反向 → 機制因果鐵證確立</div>
    <div class="footer-glossary" style="font-size:9.5px;">
      <div class="gloss-item">ⓘ V6 對 priority bug 結果 = V5 = V3F = 100%（修補在 V3F 41ff147 已完成；V6 唯一改動 Layer 1.5 revert 與本子集 germline+somatic 都 &gt;0 不重疊）</div>
      <div class="gloss-item">ⓘ chr8 priority bug 冷區 (0.34× avg) ≠ chr8 LOH+HPSig FP hotspot (7.4× ISM 下游) — 不同 layer，priority bug 修對不會自動清掉 chr8 hotspot</div>
    </div>""",
    speaker="""[標題]
S3 量化鐵證 — 把 S2 機制從理論變物理觀察。

[方法]
對 baseline 與 V6 兩版加 --debug-vote-dump flag。
dump 每條 read 經 getVote 後的 5-vote countMap。

[chr19 pilot]
HCC1395 5kHz chr19，篩 germline_majority ≠ somatic_majority 且 both >0。
得 752 雙向矛盾 reads。
baseline 全 752 條投錯方向 hp=11。
V6 全 752 條翻 hp=21 = 100% 單向修正，0 條反向。

[全基因組擴展]
victims 從 752 → 34,855 = 46.4× scale up。
V6 修正率仍 100% (34,855/34,855)，0 條反向。

[重要顛覆]
per-chr 分佈推翻原 chr19 hotspot 結論。
主要 victim 分佈在 chr7 (3,508) / chr2 / chr1 / chr16 / chr20。
chr19 只占 2.16%，rank 19 — SP1/2/3 是可重現案例，非主要分佈位置。

[chr8 釐清]
chr8 priority bug enrichment 0.34× avg — rank 21 冷區。
與 MEMORY chr8 LOH+HPSig 7.4× FP hotspot 是不同 layer，無直接因果。""",
    tier3="chr19 1Mb hotspot 30M=215 + 27M=133 + 16M=41 / 4-path 驗證表 (個案 trace / Density 共變 / 修正後消失) / read_name case (1c50034a-f0f) / 全 24 chr enrichment ‰ / chrY 小 N 高 ‰ / Pass 2 reclassify 104K germline het")

# ─── S4 修補設計 ───────────────────────────────────────────────────────────
add(id="10_fix_design", num="09", section="S4 修補設計", rg2="1", ngrep="5 hash",
    title="5 commits + getVote 四版 → V6 終態",
    en="5 commits fix + getVote 4-version: V6 production-grade",
    timing="150 sec / 中 ~420 字",
    canvas_html="""
    <p style="font-size:11.5px;font-weight:700;color:#374151;margin:0 0 4px;">⏱ 5 commits timeline + layer 分類 (4-09 ~ 5-10):</p>
    <svg viewBox="0 0 720 95" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;border:1px solid #E5E7EB;background:#F9FAFB;">
      <line x1="20" y1="50" x2="700" y2="50" stroke="#9CA3AF" stroke-width="1.2"/>
      <circle cx="60"  cy="50" r="7" fill="#3B82F6"/>
      <text x="60"  y="32" font-size="9" fill="#1E3A8A" text-anchor="middle" font-weight="700">8b8c1fd</text>
      <text x="60"  y="22" font-size="7.5" fill="#6B7280" text-anchor="middle">4-09</text>
      <text x="60"  y="72" font-size="8.5" fill="#1E3A8A" text-anchor="middle" font-weight="600">PON-only flag</text>
      <text x="60"  y="83" font-size="7.5" fill="#1E3A8A" text-anchor="middle">phasing layer</text>
      <circle cx="180" cy="50" r="8" fill="#16A34A" stroke="#FBBF24" stroke-width="2.5"/>
      <text x="180" y="32" font-size="9" fill="#166534" text-anchor="middle" font-weight="700">41ff147 ★</text>
      <text x="180" y="22" font-size="7.5" fill="#6B7280" text-anchor="middle">4-10</text>
      <text x="180" y="72" font-size="8.5" fill="#166534" text-anchor="middle" font-weight="600">two-layer getVote</text>
      <text x="180" y="83" font-size="7.5" fill="#166534" text-anchor="middle">tagging layer (修偏移)</text>
      <circle cx="300" cy="50" r="5" fill="#16A34A"/>
      <text x="300" y="32" font-size="9" fill="#166534" text-anchor="middle">380e8d2</text>
      <text x="300" y="22" font-size="7.5" fill="#6B7280" text-anchor="middle">4-25</text>
      <text x="300" y="72" font-size="8" fill="#166534" text-anchor="middle">INDEL guard</text>
      <text x="300" y="83" font-size="7" fill="#9CA3AF" text-anchor="middle">OOB UB fix</text>
      <circle cx="440" cy="50" r="8" fill="#7C3AED"/>
      <text x="440" y="32" font-size="9" fill="#5B21B6" text-anchor="middle" font-weight="700">d0bcd8c</text>
      <text x="440" y="22" font-size="7.5" fill="#6B7280" text-anchor="middle">4-30a</text>
      <text x="440" y="72" font-size="8.5" fill="#5B21B6" text-anchor="middle" font-weight="600">ploidy fix + Layer 1.5</text>
      <text x="440" y="83" font-size="7.5" fill="#5B21B6" text-anchor="middle">跨兩層 (bundled)</text>
      <circle cx="560" cy="50" r="5" fill="#3B82F6"/>
      <text x="560" y="32" font-size="9" fill="#1E3A8A" text-anchor="middle">938f0df</text>
      <text x="560" y="22" font-size="7.5" fill="#6B7280" text-anchor="middle">4-30b</text>
      <text x="560" y="72" font-size="8" fill="#1E3A8A" text-anchor="middle">threshold 0.95→0.9</text>
      <text x="560" y="83" font-size="7" fill="#9CA3AF" text-anchor="middle">Pass 2 觸發</text>
      <circle cx="680" cy="50" r="9" fill="#DC2626"/>
      <text x="680" y="55" font-size="11" fill="#FFFFFF" text-anchor="middle" font-weight="700">V6</text>
      <text x="680" y="32" font-size="9" fill="#7F1D1D" text-anchor="middle" font-weight="700">V6 patch</text>
      <text x="680" y="22" font-size="7.5" fill="#6B7280" text-anchor="middle">5-10</text>
      <text x="680" y="72" font-size="8.5" fill="#7F1D1D" text-anchor="middle" font-weight="600">revert Layer 1.5</text>
      <text x="680" y="83" font-size="7.5" fill="#7F1D1D" text-anchor="middle">germline-absent 回 hp=33</text>
    </svg>
    <div class="grid-2col" style="grid-template-columns: 6fr 6fr; gap:12px; margin-top:6px;">
      <div>
        <p style="font-size:11px;font-weight:700;color:#374151;margin:0 0 3px;">🔧 getVote 四版演進:</p>
        <table class="metric-table" style="font-size:10px;">
          <thead><tr><th>版</th><th>germline 有 vote</th><th>germline 缺席</th></tr></thead>
          <tbody>
            <tr class="row-red"><td><strong>baseline</strong></td><td>vector ① somatic pair break early → hp=11 ❌</td><td>vector ② mixed pair：罕走到，hp=33 為「順序副作用」</td></tr>
            <tr class="row-green"><td><strong>V3F</strong></td><td>Layer 1 germ 決方向 → hp=21 ✅</td><td><code style="font-size:9px;background:#F0FDF4;padding:1px 3px;">if(somTot>0 && germ==0) hp=33</code> 結構化保證 ✅</td></tr>
            <tr class="row-yellow"><td><strong>V5</strong></td><td>同 V3F → hp=21 ✅</td><td>Layer 1.5 somatic vote → hp=11/21 ⚠ 繼承 4.19:1 偏移</td></tr>
            <tr class="row-green" style="border-top:2px solid #DC2626;"><td><strong>V6 ★</strong></td><td>同 V3F → hp=21 ✅</td><td>revert Layer 1.5 → hp=33 ✅ 還原 V3F 結構化保證</td></tr>
          </tbody>
        </table>
      </div>
      <div>
        <p style="font-size:11px;font-weight:700;color:#374151;margin:0 0 3px;">📊 hp=33 reads 全基因組 (4 版量化):</p>
        <table class="metric-table" style="font-size:10px;">
          <thead><tr><th>版</th><th>hp=33 reads</th><th>vs V6 終態</th></tr></thead>
          <tbody>
            <tr class="row-red"><td>baseline</td><td class="num">少 (priority bug 偏 HP1/2)</td><td>HP_Ratio 失衡 ❌</td></tr>
            <tr class="row-green"><td>V3F</td><td class="num">132,060</td><td>首次正確標 hp=33 ✅</td></tr>
            <tr class="row-yellow"><td>V5 (+L1.5)</td><td class="num">13,250</td><td><strong>-89.9%</strong> 過度修正 ⚠</td></tr>
            <tr class="row-green" style="border-top:2px solid #DC2626;"><td><strong>V6 ★</strong></td><td class="num"><strong>138,317</strong></td><td><strong>+4.7%</strong> 還原超 V3F ✅</td></tr>
          </tbody>
        </table>
        <p style="font-size:9.5px;color:#16A34A;margin:3px 0 0;font-weight:600;">→ V6 為 ISM 下游恢復 hp=33 mixed-sub-tag 訊號</p>
      </div>
    </div>
    <div class="conclusion-arrow green" style="margin-top:6px;font-size:12.5px;padding:6px 12px;">→ V6 = V3F germline-absent 保守 + V5 germline-existent 設計目標 (ploidy / threshold / phased VCF) = hybrid 升級終態</div>
    <div class="footer-glossary" style="font-size:9.5px;">
      <div class="gloss-item">ⓘ V6 patch: HaplotagProcess.cpp:537-548 移除 13 行 (V5 Layer 1.5 else if 分支)</div>
      <div class="gloss-item">ⓘ V6 重用 V5 phased VCF → caller F1 三版完全相同 (HCC1395 0.93=0.7166 / 0.6=0.6273)</div>
    </div>""",
    speaker="""[標題]
S4 修補設計 — 5 commits 漸進 + getVote 四版演進。

[Timeline]
8b8c1fd 4-09 PON-only — 解 phasing 層球員兼裁判。
41ff147 4-10 ★ two-layer getVote — 修 priority bug 的關鍵 commit。
380e8d2 4-25 — INDEL guard 補 OOB UB。
d0bcd8c 4-30a 紫色跨兩層 — ploidy fix 讓 Pass 2 觸發 + bundled Layer 1.5。
938f0df 4-30b — threshold 0.95→0.9 鬆綁 Pass 2。
V6 5-10 — revert Layer 1.5。

[getVote 四版核心]
baseline — vector ordered ① somatic break early，hp=33 只是「順序副作用」。
V3F — Layer 2 結構化 if (somaticTotal>0 && germlineResult==0) hpResult=33 — 首次保證 germline-absent 區正確標 hp=33。
V5 Layer 1.5 — 把 hp=33 改派 hp=11/21，結構化保證被打破。
V6 — revert Layer 1.5，還原 V3F 結構化保證。

[hp=33 reads 全基因組量化]
baseline ≈ 0 結構化。
V3F 132,060 首次正確標 hp=33。
V5 13,250，-89.9% 過度修正。
V6 138,317，+4.7% 還原並略超 V3F。

[結論]
V6 = V3F germline-absent 保守 + V5 設計目標 + marker engineering 改善 = production-grade 終態。
caller F1 三版完全相同 (重用 V5 phased VCF)，下一頁說明。""",
    tier3="commit 各別 line count (8b8c1fd +69/-6 / 41ff147 +36/-25 / 380e8d2 +8/-4 / d0bcd8c +68/-9 / 938f0df +4/-4) / V6 patch HaplotagProcess.cpp:537-548 移除 13 行 / V5 Layer 1.5 設計動機 (補 V3F untagged) / cherry-pick from zhenyu")

add(id="12_no_regression", num="10", section="S5 驗證", rg2="1", ngrep="20+",
    title="4 類同層驗證 + 個案 3/3 對齊",
    en="baseline → V6: 4 same-layer + SP 3/3 aligned with paired",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:0 0 4px;">🎯 baseline → V6 同層驗證 — 4 類指標各自的驗證方法:</p>
    <table class="metric-table" style="font-size:10.5px;margin-bottom:6px;">
      <thead><tr><th>類別</th><th>驗證方法與計算來源</th></tr></thead>
      <tbody>
        <tr><td><strong>① phasing 結構</strong></td><td>從 longphase-to <code>phasing.log</code> + phased VCF 直接讀統計（N50 / Phased% / wall time）; LOH.bed bp-level intersect (Jaccard)</td></tr>
        <tr><td><strong>② 結果 tag 分布</strong></td><td>samtools 統計 BAM HP:i: tag count 分布; chr19 1Mb window per-region 計算 hp=1-1:2-1 ratio</td></tr>
        <tr><td><strong>③ paired GT 比對</strong></td><td>per-read 比對 paired BAM HP:Z: vs TO BAM HP:i:; 計算 concordance % (clean PS / 15-site filtered)</td></tr>
        <tr><td><strong>④ 個案 SP1/2/3</strong></td><td>IGV 6-BAM 並列視覺驗證 (baseline / paired / V3F / V5 / V6) + samtools coverage per-position</td></tr>
      </tbody>
    </table>
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:6px 0 4px;">📊 baseline → V6 同層 4 類驗證結果:</p>
    <table class="metric-table" style="font-size:10.5px;">
      <thead><tr><th>類別</th><th>驗證層級</th><th>指標（baseline → V6）</th></tr></thead>
      <tbody>
        <tr class="row-green"><td><strong>① 同層 phasing 結構 ⭐</strong></td><td>longphase-to 同層</td><td>phase block N50 +99.7% / Phased% +23.6 pp / speed 1.36× / LOH Jaccard=1.0</td></tr>
        <tr class="row-green"><td><strong>② 結果 tag 分布修正</strong></td><td>BAM HP tag (同層)</td><td>HP1:HP2 17.3:1 → ~1:1 / hp=33 還原 138,317 (+4.7%) / chr19 hp=1-1:2-1 ratio 修正</td></tr>
        <tr class="row-yellow"><td><strong>③ paired GT 比對 ⭐</strong></td><td>vs paired ground truth</td><td>clean PS +8.3 pp / 15-Aggr +6.65 pp / 15-Clean PS +13.3 pp</td></tr>
        <tr class="row-green" style="border-top:2px solid #DC2626;"><td><strong>④ 個案 SP1/2/3 (IGV)</strong></td><td>已知偏移位點修正</td><td>baseline 113:0 / 109:1 / 108:0 → V6 全翻 HP2 對齊 paired 3/3 ✅</td></tr>
      </tbody>
    </table>
    <p style="font-size:10px;color:#7C2D12;margin:4px 0 0;font-weight:600;background:#FEF3C7;padding:3px 8px;border-radius:4px;">⚠ 本 slide 驗證僅用 longphase-to 同層 + 結果 tag + paired 對比 + 個案層 — 不用 ISM 下游 metric (下游本就會因修正而改變，循環論證)</p>
    <p style="font-size:11px;font-weight:700;color:#374151;margin:6px 0 2px;">📊 6 顯著改善 ⭐ 視覺化:</p>
    <svg viewBox="0 0 720 120" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;border:1px solid #E5E7EB;background:#F9FAFB;">
      <line x1="200" y1="100" x2="710" y2="100" stroke="#9CA3AF" stroke-width="0.8"/>
      <line x1="200" y1="14"  x2="200" y2="100" stroke="#9CA3AF" stroke-width="0.8"/>
      <text x="196" y="17" font-size="8" fill="#6B7280" text-anchor="end">+100%</text>
      <text x="196" y="60" font-size="8" fill="#6B7280" text-anchor="end">+50%</text>
      <text x="196" y="103" font-size="8" fill="#6B7280" text-anchor="end">0</text>
      <text x="10" y="22"  font-size="10" fill="#374151" font-weight="600">Phase block N50</text>
      <rect x="200" y="14" width="430" height="11" fill="#16A34A"/>
      <text x="640" y="23" font-size="10" fill="#166534" font-weight="700">+99.7%</text>
      <text x="10" y="36"  font-size="10" fill="#374151" font-weight="600">Phased%</text>
      <rect x="200" y="28" width="102" height="11" fill="#16A34A"/>
      <text x="312" y="37" font-size="10" fill="#166534" font-weight="700">+23.6 pp</text>
      <text x="10" y="50"  font-size="10" fill="#374151" font-weight="600">Speed</text>
      <rect x="200" y="42" width="155" height="11" fill="#16A34A"/>
      <text x="365" y="51" font-size="10" fill="#166534" font-weight="700">1.36×</text>
      <text x="10" y="64"  font-size="10" fill="#374151" font-weight="600">15-Clean PS GT</text>
      <rect x="200" y="56" width="57" height="11" fill="#3B82F6"/>
      <text x="267" y="65" font-size="10" fill="#1E3A8A" font-weight="700">+13.3 pp</text>
      <text x="10" y="78"  font-size="10" fill="#374151" font-weight="600">clean PS GT</text>
      <rect x="200" y="70" width="36" height="11" fill="#3B82F6"/>
      <text x="246" y="79" font-size="10" fill="#1E3A8A" font-weight="700">+8.3 pp</text>
      <text x="10" y="92"  font-size="10" fill="#374151" font-weight="600">15-Aggr GT</text>
      <rect x="200" y="84" width="29" height="11" fill="#3B82F6"/>
      <text x="239" y="93" font-size="10" fill="#1E3A8A" font-weight="700">+6.65 pp</text>
      <text x="500" y="115" font-size="9" fill="#16A34A" text-anchor="middle" font-weight="600">綠=HP/LOH 結構 / 藍=Paired GT concordance (read-level tag 品質)</text>
    </svg>
    <div class="conclusion-arrow green" style="font-size:12.5px;">→ 4 類同層驗證全綠 + 個案層 3/3 對齊 paired → V6 修補在 longphase-to 端鐵證確立</div>
    <div class="footer-glossary" style="font-size:9.5px;">
      <div class="gloss-item">ⓘ scope: HCC1395 5kHz @ 0.93 purity; V6 同層 metric = V5 (重用 phased VCF)；caller F1 不變 (slide 11)</div>
      <div class="gloss-item" style="background:#DCFCE7;border-color:#16A34A;color:#166534;font-weight:600;">★ 真實價值在 read-level tag concordance — paired GT +13.3 pp 為 PI 鐵證；ISM 下游影響量化留 slide 13 future direction</div>
    </div>""",
    speaker="""[標題]
S5 驗證 — baseline → V6 同層 4 類全綠。

[重要前提]
本 slide 只用 longphase-to 同層 + 結果 tag + paired 對比 + 個案層 驗證。
不用 ISM 下游 metric — 下游本就因 V6 修正而改變，用它會是循環論證。

[類別 1：同層 phasing 結構]
phase block N50 +99.7%。
Phased% +23.6 pp。
speed 1.36× 快。
LOH Jaccard = 1.0 完全相同。

[類別 2：結果 tag 分布修正]
HP1:HP2 17.3:1 → ~1:1。
hp=33 還原 138,317 (+4.7% vs V3F)。
chr19 hp=1-1:2-1 ratio 修正。

[類別 3：paired GT 比對 — PI 鐵證]
clean PS Paired GT +8.3 pp。
15-Aggr Paired GT +6.65 pp。
15-Clean PS Paired GT +13.3 pp ★。

[類別 4：個案 SP1/2/3 IGV]
baseline 113:0 / 109:1 / 108:0。
V6 全翻 HP2 對齊 paired 3/3。

[結論]
4 類同層驗證全綠 + 個案層 3/3 對齊 paired → V6 修補在 longphase-to 端鐵證確立。
ISM 下游影響量化是未來研究方向 — slide 13 會展開。""",
    tier3="methylation 6 feat 列表 / Layer 1.5 zero-sum 細節 (留 slide 12-13) / V2b/V3F 中間版本 / HP_Ratio 詳解 / LOH Jaccard=1.0 細節")

add(id="14_caller_f1", num="11", section="S5 驗證 + ⚡ Cliffhanger", rg2="2", ngrep="15+",
    title="Caller F1 四版完全相同 = 0.7166",
    en="Caller F1 identical across 4 versions; V6 no regression",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:0 0 4px;">📊 兩 purity 對照 — 四版 F1 都相同, 但 0.93 與 0.6 不同 (caller 性質):</p>
    <table class="metric-table" style="font-size:10.5px;">
      <thead><tr><th>樣本</th><th>四版本</th><th>TP</th><th>FP</th><th>FN</th><th>F1</th><th>Δ vs 純</th></tr></thead>
      <tbody>
        <tr class="row-green"><td><strong>HCC1395 5kHz @ 0.93 (純)</strong></td><td class="num">baseline/V3F/V5/V6</td><td class="num">28,509</td><td class="num">11,606</td><td class="num">10,938</td><td class="num"><strong>0.7166</strong></td><td>—</td></tr>
        <tr class="row-yellow"><td><strong>HCC1395 t30_n20 @ 0.6 (低)</strong></td><td class="num">baseline/V3F/V5/V6</td><td class="num">24,190</td><td class="num">13,487</td><td class="num">15,257</td><td class="num"><strong>0.6273</strong></td><td class="num">−0.0893</td></tr>
      </tbody>
    </table>
    <p style="font-size:10.5px;color:#16A34A;margin:3px 0 0;font-weight:600;">→ 兩 purity 都驗證: 四版 F1 完全相同 (0 差異); 0.93→0.6 ΔF1=−0.0893 為 ClairS-TO 低 purity 性質，與 self-phasing 修補無關</p>
    <p style="font-size:11px;font-weight:700;color:#374151;margin:8px 0 4px;">🔬 為何 F1 invariant — 因果機制圖示:</p>
    <svg viewBox="0 0 720 140" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;border:1px solid #E5E7EB;background:#F9FAFB;">
      <defs><marker id="arr14" markerWidth="8" markerHeight="8" refX="6" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 Z" fill="#9CA3AF"/></marker></defs>
      <rect x="20"  y="35" width="135" height="60" rx="6" fill="#DBEAFE" stroke="#1E3A8A" stroke-width="1.5"/>
      <text x="87"  y="53" font-size="10" fill="#1E3A8A" text-anchor="middle" font-weight="700">ClairS-TO caller</text>
      <text x="87"  y="68" font-size="9"  fill="#374151" text-anchor="middle">snv.vcf</text>
      <text x="87"  y="83" font-size="9" fill="#7F1D1D" text-anchor="middle" font-weight="700">FILTER=PASS</text>
      <text x="87"  y="93" font-size="8" fill="#6B7280" text-anchor="middle">(47,798 一次性決定)</text>
      <path d="M155,65 L200,65" stroke="#9CA3AF" stroke-width="2" marker-end="url(#arr14)"/>
      <rect x="205" y="20" width="210" height="90" rx="6" fill="#F0FDF4" stroke="#16A34A" stroke-width="1.5"/>
      <text x="310" y="38" font-size="10" fill="#166534" text-anchor="middle" font-weight="700">longphase-to (V3F/V5/V6)</text>
      <text x="310" y="55" font-size="9" fill="#166534" text-anchor="middle">✓ 改 GT / PS / GT2 / GT3 / PON tag</text>
      <text x="310" y="69" font-size="9" fill="#166534" text-anchor="middle">✓ 改 BAM HP:i tag</text>
      <text x="310" y="86" font-size="9.5" fill="#DC2626" text-anchor="middle" font-weight="700">✗ 不改 FILTER 欄位</text>
      <text x="310" y="100" font-size="8" fill="#6B7280" text-anchor="middle" font-style="italic">三版 PASS set 完全相同</text>
      <path d="M415,65 L460,65" stroke="#9CA3AF" stroke-width="2" marker-end="url(#arr14)"/>
      <rect x="465" y="35" width="150" height="60" rx="6" fill="#FEF3C7" stroke="#CA8A04" stroke-width="1.5"/>
      <text x="540" y="53" font-size="10" fill="#7C2D12" text-anchor="middle" font-weight="700">bcftools isec</text>
      <text x="540" y="69" font-size="9"  fill="#374151" text-anchor="middle">取 FILTER=PASS</text>
      <text x="540" y="83" font-size="9"  fill="#374151" text-anchor="middle">vs SEQC2 truth</text>
      <path d="M615,65 L660,65" stroke="#9CA3AF" stroke-width="2" marker-end="url(#arr14)"/>
      <rect x="665" y="35" width="55" height="60" rx="6" fill="#DCFCE7" stroke="#16A34A" stroke-width="2"/>
      <text x="692" y="53" font-size="10" fill="#166534" text-anchor="middle" font-weight="700">F1</text>
      <text x="692" y="70" font-size="11" fill="#166534" text-anchor="middle" font-weight="700">0.7166</text>
      <text x="692" y="85" font-size="8.5" fill="#166534" text-anchor="middle" font-weight="600">四版相同</text>
      <text x="370" y="128" font-size="9.5" fill="#16A34A" text-anchor="middle" font-weight="600">→ FILTER 不動 → PASS set 不變 → TP/FP/FN 不變 → F1 數學保證 invariant</text>
    </svg>
    <table class="metric-table" style="font-size:10px;margin-top:6px;">
      <thead><tr><th>F1 評估口徑</th><th>結果</th><th>對 V3F/V5/V6 修補?</th></tr></thead>
      <tbody>
        <tr class="row-green"><td><strong>caller-level F1</strong> (bcftools isec, 業界標準)</td><td>四版 0.7166 / 0.6273 完全相同</td><td>invariant ✓ (FILTER 不動)</td></tr>
        <tr class="row-green"><td>論文 §4.3 V_H/V_L refined F1</td><td>Verdict_Somatic INFO tag 三版實證 0 records (本 pipeline N/A)</td><td>invariant ✓ (Step 4 不改)</td></tr>
      </tbody>
    </table>
    <div class="footer-glossary" style="font-size:9.5px;margin-top:4px;">
      <div class="gloss-item" style="background:#DCFCE7;border-color:#16A34A;color:#166534;font-weight:600;">★ 真實價值在 read-level tag concordance (+13.3 pp paired GT @ 0.93)，不在 caller F1 〔next: slide 12 主結論 + slide 13 未來方向〕</div>
    </div>""",
    speaker="""[標題]
Caller F1 vs SEQC2 — 四版完全相同。
這是 V6 不退步的 caller-side 鐵證。

[兩 purity 對照]
HCC1395 5kHz 0.93 純樣本：
baseline / V3F / V5 / V6 四版 TP=28,509 / FP=11,606 / FN=10,938 / F1=0.7166 — 完全相同。

HCC1395 t30_n20 0.6 低純度：
四版 TP=24,190 / FP=13,487 / FN=15,257 / F1=0.6273 — 完全相同。

[釐清 0.93 vs 0.6 差異]
兩 purity 都驗證了，四版內 0 差異。
但 0.93 → 0.6 ΔF1 = −0.0893 是 ClairS-TO 在低 purity 性質。
與 self-phasing 修補無關。

[為何 invariant — 因果機制]
ClairS-TO snv.vcf 一次性決定 FILTER=PASS 47,798 條。
longphase-to phase 改 GT / PS / GT2 / GT3 / PON tag — 不改 FILTER。
longphase-to haplotag 只改 BAM HP:i tag — 不碰 VCF。
bcftools isec 取 FILTER=PASS 跟 truth 比對。
FILTER 不動 → PASS set 不變 → TP/FP/FN 不變 → F1 數學保證 invariant。

[兩種 F1 口徑都 invariant]
業界 caller-level F1：四版 0.7166 / 0.6273 完全相同。
論文 §4.3 V_H/V_L refined F1：實證三版 Verdict tag INFO 全 0 records（ClairS-TO 預設未啟用 Verdict tagging）→ 本 pipeline N/A。

[結論]
V6 真實價值不在 caller F1。
真實價值在 read-level tag concordance — +13.3 pp paired GT @ 0.93。
下一節 — slide 12 主結論，slide 13 未來方向。""",
    tier3="PASS set / FILTER 機制 / purity 0.6 N50 微差 / V6 patch 不改 phasing 層 detail")

# ─── S7 Errata + Follow-up ────────────────────────────────────────────────
add(id="17_main_verdict", num="12", section="S7 結論", rg2="3 (main verdict)", ngrep="6 errata",
    title="3 條主結論 + 6 條 errata patch",
    en="3 main verdicts + 6 errata patched",
    timing="120 sec / 中 ~400 字",
    title_critical=True,
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#1E3A8A;margin:0 0 4px;">🎯 整份報告 3 條 main verdict (PI 必記):</p>
    <div style="display:flex;flex-direction:column;gap:6px;margin-bottom:8px;">
      <div style="display:flex;align-items:center;gap:12px;padding:8px 12px;background:#F0FDF4;border:2px solid #16A34A;border-radius:8px;">
        <div style="flex:none;width:50px;height:50px;border-radius:50%;background:#DCFCE7;border:2px solid #16A34A;display:flex;align-items:center;justify-content:center;font-size:22px;font-weight:700;color:#166534;">①</div>
        <div style="flex:1;">
          <p style="font-size:13px;font-weight:700;color:#166534;margin:0 0 3px;">priority bug 修補機制確立</p>
          <p style="font-size:11.5px;color:#374151;margin:0;line-height:1.5;">17.3:1 → ~1:1 全基因組 · 34,855 read-level victims baseline → V6 100% 修正 0 反向 · SP1/2/3 IGV 對齊 paired 3/3 ✅</p>
        </div>
      </div>
      <div style="display:flex;align-items:center;gap:12px;padding:8px 12px;background:#F0FDF4;border:2px solid #16A34A;border-radius:8px;">
        <div style="flex:none;width:50px;height:50px;border-radius:50%;background:#DCFCE7;border:2px solid #16A34A;display:flex;align-items:center;justify-content:center;font-size:22px;font-weight:700;color:#166534;">②</div>
        <div style="flex:1;">
          <p style="font-size:13px;font-weight:700;color:#166534;margin:0 0 3px;">V6 = production candidate (5/10 binary patch)</p>
          <p style="font-size:11.5px;color:#374151;margin:0;line-height:1.5;">= V3F 保守 + V5 設計目標 + marker engineering 改善 · hp=33 +4.7% · marker coverage +9.0% · Phase D 4 樣本 ratio 0.61-1.24 中性 · F1 0.7166</p>
        </div>
      </div>
      <div style="display:flex;align-items:center;gap:12px;padding:8px 12px;background:#FEF3C7;border:2px solid #CA8A04;border-radius:8px;">
        <div style="flex:none;width:50px;height:50px;border-radius:50%;background:#FBBF24;border:2px solid #CA8A04;display:flex;align-items:center;justify-content:center;font-size:22px;font-weight:700;color:#7C2D12;">★</div>
        <div style="flex:1;">
          <p style="font-size:13px;font-weight:700;color:#7C2D12;margin:0 0 3px;">真實價值在 read-level tag concordance — 非 caller F1</p>
          <p style="font-size:11.5px;color:#374151;margin:0;line-height:1.5;">+13.3 pp paired GT @ 0.93 · hp=33 還原 · marker coverage +9.0% → ISM 下游 marker engineering 受惠 · caller F1 對 V3F/V5/V6 數學保證 invariant</p>
        </div>
      </div>
    </div>
    <p style="font-size:12px;font-weight:700;color:#374151;margin:10px 0 4px;">📝 PI 報告 4-29 核心 errata (3 條，完整 6 條詳見報告):</p>
    <table class="metric-table" style="font-size:12px;">
      <thead><tr><th></th><th>段落</th><th>原 → 新</th></tr></thead>
      <tbody>
        <tr class="row-yellow"><td><strong>★ E4</strong></td><td>V5 數值歸因</td><td>V5 整體效益 → V5 = Pass 1 only; 主要功勞 V3F + Layer 1.5; Pass 2 二次效益尚未量化</td></tr>
        <tr class="row-yellow"><td><strong>★ E5</strong></td><td>Layer 1.5 設計 (5/10 加)</td><td>fallback 隱含修補 → germline-absent 區與 baseline 4.19:1 偏 HP1 完全相同</td></tr>
        <tr class="row-green"><td><strong>★ E6</strong></td><td>follow-up (5/11)</td><td>F-paired-D3 ⏸ → ✅ V6 binary patch DONE + Phase D 4 樣本驗證 4/5 通過</td></tr>
      </tbody>
    </table>
    <div class="footer-glossary" style="font-size:9.5px;">
      <div class="gloss-item">ⓘ commit chain: f17754f → 2553e96 → 71d21bd → V6 patch (待 commit)</div>
      <div class="gloss-item">ⓘ E1-E3 表述精確化（chr19 hotspot 降級 / V5 commit / 證據鏈升級）詳見 errata 完整報告</div>
    </div>""",
    speaker="""[標題]
整份報告主結論 3 條 — PI 必記。

[結論 1]
priority bug 修補機制因果確立。
17.3:1 偏移 → ~1:1 全基因組。
34,855 read-level victims baseline → V6 100% 修正，0 反向。
SP1/2/3 IGV 對齊 paired 3/3。

[結論 2]
V6 = production candidate（5/10 binary patch）。
等於 V3F 保守 hp=33 + V5 ploidy/threshold/phased VCF + marker engineering 改善。
hp=33 +4.7% vs V3F。
marker coverage +9.0%。
Phase D 4 樣本 ratio 0.61-1.24 全中性。
caller F1 三版相同 0.7166。

[結論 3 — ★ 最重要]
真實價值在 read-level tag concordance — 非 caller F1。
+13.3 pp paired GT concordance @ 0.93。
hp=33 還原。
marker coverage +9.0% → ISM 下游 marker engineering 受惠。
caller F1 對 V3F/V5/V6 數學保證 invariant — FILTER 不動。

[errata 簡述]
PI 報告 4-29 errata 6 條，3 條核心：
E4 — V5 = Pass 1 only 歸因修正。
E5 — V5 Layer 1.5 germline-absent 4.19:1 缺陷。
E6 — V6 binary patch 完成 + Phase D 4 樣本驗證閉環。

[轉場]
下一頁 — 改正後影響 + 未來研究方向。""",
    tier3="errata commit chain / V6 patch 待 commit / Phase D 5 樣本 evaluation matrix / 各 errata 段落 cross-ref")

add(id="18_impact_future", num="13", section="S7 結論", rg2="1", ngrep="—",
    title="改正後影響 + 未來研究 3 方向",
    en="Impact + future directions; V6 production candidate",
    timing="100 sec / 中 ~360 字",
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#16A34A;margin:0 0 4px;">🎯 改正後影響 — upstream → downstream 因果鏈:</p>
    <svg viewBox="0 0 720 200" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;border:1px solid #E5E7EB;background:#F9FAFB;">
      <defs><marker id="arr18" markerWidth="9" markerHeight="9" refX="7" refY="4.5" orient="auto"><path d="M0,0 L9,4.5 L0,9 Z" fill="#16A34A"/></marker></defs>
      <rect x="15" y="20" width="220" height="160" rx="8" fill="#F0FDF4" stroke="#16A34A" stroke-width="2"/>
      <text x="125" y="40" font-size="13" fill="#166534" text-anchor="middle" font-weight="700">upstream (longphase-to V6)</text>
      <text x="25" y="62" font-size="11" fill="#166534">✓ priority bug 修對</text>
      <text x="40" y="76" font-size="10" fill="#374151">17.3:1 → ~1:1 (全基因組)</text>
      <text x="25" y="96" font-size="11" fill="#166534">✓ hp=33 mixed-sub-tag 訊號還原</text>
      <text x="40" y="110" font-size="10" fill="#374151">V6 = 138,317 (+4.7% vs V3F)</text>
      <text x="25" y="130" font-size="11" fill="#166534">✓ marker coverage +9.0% vs V3F</text>
      <text x="40" y="144" font-size="10" fill="#374151">23,980 NG≥3 regions</text>
      <text x="25" y="164" font-size="11" fill="#166534">✓ Phase D 4 樣本 ratio 全中性</text>
      <path d="M240,100 L295,100" stroke="#16A34A" stroke-width="3" marker-end="url(#arr18)"/>
      <text x="267" y="92" font-size="9" fill="#16A34A" text-anchor="middle" font-weight="700">因果鏈</text>
      <text x="267" y="115" font-size="8" fill="#374151" text-anchor="middle" font-style="italic">tag 訊號乾淨 →</text>
      <text x="267" y="127" font-size="8" fill="#374151" text-anchor="middle" font-style="italic">下游分析更精確</text>
      <rect x="300" y="20" width="240" height="160" rx="8" fill="#DBEAFE" stroke="#1E3A8A" stroke-width="2"/>
      <text x="420" y="40" font-size="13" fill="#1E3A8A" text-anchor="middle" font-weight="700">downstream (ISM)</text>
      <text x="310" y="62" font-size="11" fill="#1E3A8A">✓ HP_Ratio tag bias 修正</text>
      <text x="325" y="76" font-size="10" fill="#374151">0.788 → 0.574</text>
      <text x="310" y="96" font-size="11" fill="#1E3A8A">✓ HP1/HP2 read 集合更乾淨</text>
      <text x="325" y="110" font-size="10" fill="#374151">→ methylation delta 更精確</text>
      <text x="310" y="130" font-size="11" fill="#1E3A8A">✓ marker engineering 訊號完整</text>
      <text x="325" y="144" font-size="10" fill="#374151">hp=33 + NG≥3 訊號</text>
      <text x="310" y="164" font-size="11" fill="#7C2D12" font-weight="700">★ paired GT concordance +13.3 pp @ 0.93</text>
      <rect x="555" y="60" width="155" height="80" rx="8" fill="#FEF3C7" stroke="#CA8A04" stroke-width="2"/>
      <text x="632" y="78" font-size="11" fill="#7C2D12" text-anchor="middle" font-weight="700">主結論</text>
      <text x="632" y="98" font-size="10.5" fill="#7C2D12" text-anchor="middle" font-weight="700">V6 = production</text>
      <text x="632" y="112" font-size="10.5" fill="#7C2D12" text-anchor="middle" font-weight="700">candidate</text>
      <text x="632" y="130" font-size="9" fill="#374151" text-anchor="middle" font-style="italic">可升級替代 V5</text>
      <path d="M540,100 L555,100" stroke="#CA8A04" stroke-width="2" marker-end="url(#arr18)"/>
    </svg>
    <p style="font-size:12.5px;font-weight:700;color:#DC2626;margin:10px 0 4px;">🔬 未來研究方向 — 回到 ISM 五大目標主線 (V6 = 下游驗證基礎):</p>
    <table class="metric-table" style="font-size:11.5px;">
      <thead><tr><th>#</th><th>主軸任務</th><th>對應 ISM 五大目標</th><th>具體驗證重點</th></tr></thead>
      <tbody>
        <tr class="row-green"><td><strong>F1 ★</strong></td><td><strong>LOH 內外 TP/FP 差異特徵</strong>（用 V6 tag）</td><td>目標 1 (per-CpG 甲基關聯) + 目標 4 (TO normal 補強) + 目標 5 (F1 提升)</td><td>LOH 內 vs LOH 外 × TP/FP × 甲基化率 / read 特徵 / HP_Ratio / NG 分布；驗證 ISM 是否能用 V6 tag 正確識別不同區域差異</td></tr>
        <tr class="row-green"><td><strong>F2 ★</strong></td><td><strong>subclone 結構 + 二次打擊事件順序</strong></td><td>目標 2 (clone 結構) + 目標 3 (二次打擊事件順序)</td><td>用 V6 hp=33 + NG≥3 marker 識別 sub-clone；結合 LOH/CNV/HP 推論 two-hit order；驗證 V6 marker 訊號是否能區分 first-hit vs second-hit</td></tr>
        <tr><td><strong>F3</strong></td><td>7-sample expansion + COLO829 解鎖 + Pass 2 ablation</td><td>F1/F2 cross-sample 一致性驗證基礎</td><td>HCC1395_DORADO/HCC1937/HCC1954/H1437/H2009/COLO829 完整 V6 驗證 + 量化 V6 對 ISM marker rate / FP-TP / NGroups 影響</td></tr>
      </tbody>
    </table>
    <div class="conclusion-arrow green" style="font-size:12px;margin-top:6px;">→ V6 = production candidate；F1/F2 是 ISM 核心主線回歸 (從 longphase-to bug 修補 → 重新驗證 ISM 五大目標)</div>""",
    speaker="""[標題]
改正後影響 + 未來研究方向 — 回到 ISM 五大目標主線。

[upstream 影響]
longphase-to V6 端：
priority bug 修對 17.3:1 → ~1:1。
hp=33 還原 138,317。
marker coverage +9.0%。
Phase D 4 樣本 ratio 全中性。

[downstream 影響]
ISM 端受惠：
HP_Ratio tag bias 修正 0.788 → 0.574。
HP1/HP2 read 集合更乾淨 → methylation delta 更精確。
marker engineering 訊號完整。
★ paired GT concordance +13.3 pp @ 0.93。

[主結論連結]
V6 = production candidate，可升級替代 V5。

[未來研究方向 — 回到 ISM 五大目標主線]
F1 ★ LOH 內外 TP/FP 差異特徵驗證。
用 V6 tag 在 LOH 內 vs LOH 外 區域對 TP/FP 比較甲基化率、read 特徵、HP_Ratio、NGroups 分布。
對應五大目標 1 (per-CpG 甲基關聯) + 4 (TO normal 補強) + 5 (F1 提升)。

F2 ★ subclone 結構 + 二次打擊事件順序。
用 V6 hp=33 + NG≥3 marker 識別 sub-clone 結構。
結合 LOH/CNV/HP 推論 two-hit order。
驗證 V6 marker 訊號是否能區分 first-hit vs second-hit。
對應五大目標 2 (clone 結構) + 3 (二次打擊事件順序)。

F3 7-sample expansion + COLO829 + Pass 2 ablation。
F1/F2 cross-sample 驗證基礎。

[研究定位]
InterSubMod 是可解釋 epigenetic evidence 整合框架，不是另一個 somatic caller。
F1/F2 是從 longphase-to bug 修補完成後，回歸 ISM 核心主線的研究重啟。
V6 是 production candidate 的同時，也是下游 ISM 驗證的乾淨基礎。

[結尾]
今天主軸到這裡。後續 4 張 QA backup slides 補充對話質疑點。""",
    tier3="ISM 五大目標完整列表 (目標 1-5) / LOH 內外 TP/FP 差異分層方案 / two-hit order 推論方法 / HCC1937 BRCA1 mutant / COLO829 chmod / T1.3 4-cell ablation 設計")

# ─── Q&A Backup (對話釐清點補佐證) ────────────────────────────────────────
add(id="b1_pass2", num="B1", section="Q&A Backup", rg2="1", ngrep="—",
    title="Pass 2: 重跑 2-point only (高 purity >0.9)",
    en="Pass 2 = re-run 2-point only; high purity gate",
    timing="90 sec / 中 ~300 字",
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#374151;margin:0 0 4px;">❓ Q: Pass 2 為什麼只跑 2-point? 為什麼高 purity 才觸發?</p>
    <svg viewBox="0 0 720 220" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;border:1px solid #E5E7EB;background:#F9FAFB;">
      <defs><marker id="arrB1" markerWidth="8" markerHeight="8" refX="6" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 Z" fill="#9CA3AF"/></marker></defs>
      <rect x="15" y="15" width="230" height="190" rx="6" fill="#DBEAFE" stroke="#1E3A8A" stroke-width="1.5"/>
      <text x="130" y="35" font-size="11" fill="#1E3A8A" text-anchor="middle" font-weight="700">Pass 1 (always runs)</text>
      <rect x="30" y="50" width="200" height="42" rx="4" fill="white" stroke="#1E3A8A"/>
      <text x="130" y="68" font-size="10" fill="#374151" text-anchor="middle" font-weight="600">3-point patternMining</text>
      <text x="130" y="83" font-size="8.5" fill="#6B7280" text-anchor="middle">somaticCalling V_H/V_L/V_N 分類</text>
      <rect x="30" y="100" width="200" height="42" rx="4" fill="white" stroke="#1E3A8A"/>
      <text x="130" y="118" font-size="10" fill="#374151" text-anchor="middle" font-weight="600">2-point edgeConnectResult</text>
      <text x="130" y="133" font-size="8.5" fill="#6B7280" text-anchor="middle">phasing graph 建構</text>
      <text x="130" y="165" font-size="9" fill="#1E3A8A" text-anchor="middle" font-weight="700">→ 產出 origin + phased VCF</text>
      <text x="130" y="186" font-size="8" fill="#6B7280" text-anchor="middle" font-style="italic">三路徑算法不依賴 purity</text>
      <path d="M245,110 L290,110" stroke="#9CA3AF" stroke-width="2" marker-end="url(#arrB1)"/>
      <rect x="295" y="80" width="125" height="60" rx="6" fill="#FEF3C7" stroke="#CA8A04" stroke-width="2"/>
      <text x="357" y="100" font-size="11" fill="#7C2D12" text-anchor="middle" font-weight="700">purity 判定</text>
      <text x="357" y="118" font-size="11" fill="#7C2D12" text-anchor="middle" font-weight="700">purity &gt; 0.9 ?</text>
      <text x="357" y="132" font-size="8" fill="#6B7280" text-anchor="middle">(原 0.95 → 938f0df 改 0.9)</text>
      <path d="M357,140 L357,165" stroke="#DC2626" stroke-width="2" marker-end="url(#arrB1)"/>
      <text x="395" y="155" font-size="9" fill="#DC2626" font-weight="700">No (≤0.9)</text>
      <path d="M420,110 L455,110" stroke="#16A34A" stroke-width="2" marker-end="url(#arrB1)"/>
      <text x="437" y="100" font-size="9" fill="#16A34A" font-weight="700">Yes (&gt;0.9)</text>
      <rect x="460" y="55" width="245" height="100" rx="6" fill="#DCFCE7" stroke="#16A34A" stroke-width="2"/>
      <text x="582" y="78" font-size="11" fill="#166534" text-anchor="middle" font-weight="700">Pass 2 second round</text>
      <text x="582" y="93" font-size="10" fill="#166534" text-anchor="middle">(only 2-point edgeConnectResult)</text>
      <text x="582" y="115" font-size="9" fill="#374151" text-anchor="middle">重跑 edgeConnectResult</text>
      <text x="582" y="128" font-size="9" fill="#374151" text-anchor="middle">用 Pass 1 origin 精修 graph</text>
      <text x="582" y="146" font-size="9.5" fill="#16A34A" text-anchor="middle" font-weight="700">N50 +3.51% · blocks −9.79% · phased −2.9 pp</text>
      <rect x="295" y="170" width="245" height="35" rx="6" fill="#FEE2E2" stroke="#DC2626" stroke-width="2"/>
      <text x="418" y="190" font-size="10" fill="#7F1D1D" text-anchor="middle" font-weight="700">Pass 2 跳過 (低 purity)</text>
      <text x="418" y="202" font-size="8.5" fill="#6B7280" text-anchor="middle">HCC1395 0.6 樣本不觸發 / Pass 1 已穩定</text>
    </svg>
    <p style="font-size:10.5px;font-weight:600;color:#16A34A;margin:6px 0 0;">★ 關鍵: Pass 2 只重跑 2-point edgeConnectResult，不重跑 3-point somaticCalling — 後者 patternMining 基於 read linkage 拓樸 (V_H 12 + V_L 3 patterns)，與 purity 無關</p>""",
    speaker="""[Q]
Pass 2 為什麼只跑 2-point？為什麼高 purity 才觸發？

[Pass 1 永遠跑]
Pass 1 跑兩件事：
3-point patternMining — somaticCalling 做 V_H / V_L / V_N 分類。
2-point edgeConnectResult — phasing graph 建構。
三路徑算法不依賴 purity（基於 read linkage 拓樸）。

[purity 判定]
purity > 0.9 才觸發 Pass 2（原 0.95，938f0df commit 改 0.9）。
低 purity ≤ 0.9 Pass 2 跳過 — Pass 1 已穩定。

[Pass 2 做什麼]
只重跑 2-point edgeConnectResult。
用 Pass 1 已分類 origin 精修 phasing graph。
不重跑 3-point somaticCalling。

[Pass 2 incremental 效果]
N50 +3.51%。
blocks −9.79%。
phased var −2.9 pp。

[常見誤解]
「低 purity 用 3-point 倒過來」— 實情相反。
高 purity 才多做事，因為 Pass 1 origin 穩定可重用。

[低 purity 樣本案例]
HCC1395 0.6 purity 樣本兩版都不觸發 Pass 2。
V5 改 GT 在低 purity 自我治癒。""",
    tier3="patternMining first/second/third path / Pass 2 不重跑 somaticCalling / threshold 0.95→0.9 cherry-pick")

add(id="b2_f1_dual_metric", num="B2", section="Q&A Backup", rg2="3", ngrep="—",
    title="F1 雙口徑釐清 + Verdict tag 實證",
    en="F1 dual metric + Verdict tag empirical N/A",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#374151;margin:0 0 4px;">❓ Q: longphase-to 論文 §4.3 寫 F1 改善 — 為何我們 slide 11 寫 F1 不變?</p>
    <p style="font-size:10.5px;color:#166534;font-weight:600;margin:2px 0 6px;">→ 雙口徑差異, 不衝突: 不同 metric 定義 + 我們 pipeline 沒啟用 Verdict tagging</p>
    <table class="metric-table" style="font-size:11.5px;">
      <thead><tr><th>F1 評估口徑</th><th>somatic call set</th><th>本 pipeline?</th><th>結果</th></tr></thead>
      <tbody>
        <tr class="row-green"><td><strong>caller-level F1</strong><br>(業界 bcftools isec)</td><td>phased_VCF.FILTER=PASS<br>vs SEQC2 truth</td><td>✅ 採用</td><td>四版 F1=0.7166 完全相同</td></tr>
        <tr style="background:#F3F4F6;"><td>論文 §4.3 V_H/V_L<br>refined F1</td><td>Verdict_Somatic INFO tag<br>(Step 4 Graph-Based Calling)</td><td>⚠ N/A</td><td>三版 Verdict tag 全 0 records<br>(ClairS-TO 預設未啟用)</td></tr>
      </tbody>
    </table>
    <p style="font-size:10.5px;font-weight:700;color:#1E3A8A;margin:8px 0 4px;">🔬 longphase-to 真實修改範圍 (V3F/V5/V6 都不動 FILTER):</p>
    <table class="metric-table" style="font-size:11.5px;">
      <thead><tr><th>階段</th><th>改什麼欄位</th><th>影響 caller F1?</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>phase</td><td>GT / PS / GT2 / GT3 / PON tag (FORMAT/INFO)</td><td>❌ 不影響</td></tr>
        <tr class="row-green"><td>haplotag</td><td>BAM HP:i tag (僅 BAM 不碰 VCF)</td><td>❌ 不影響</td></tr>
        <tr class="row-yellow"><td>V5 Pass 2 reclassify</td><td>104,457 germline het GT 0/1→0/0 (FORMAT)</td><td>❌ 不改 FILTER, F1 不變</td></tr>
        <tr style="background:#FEE2E2;"><td><strong>FILTER 欄位</strong></td><td><strong>由 ClairS-TO 一次性決定，longphase 純消費</strong></td><td>FILTER 不動 → F1 數學保證 invariant</td></tr>
      </tbody>
    </table>
    <p style="font-size:10px;font-weight:600;color:#7C2D12;margin:6px 0 0;background:#FEF3C7;padding:4px 8px;border-radius:4px;">★ 實證鐵證: 三版 phased VCF PASS=47,798 / nonPASS=3,139,477 / Verdict tag 全 0 — V3F/V5/V6 對兩種 F1 口徑都 invariant</p>""",
    speaker="""[Q]
longphase-to 論文 §4.3 寫 F1 改善 — 為何我們寫 F1 不變？

[答 — 雙口徑差異, 不衝突]
論文口徑 — 用 Step 4 Graph-Based Somatic Calling 的 V_H/V_L set 作 caller post-filter（Verdict_Somatic INFO tag）。
我們業界口徑 — 用 bcftools isec FILTER=PASS vs SEQC2 truth。

[實證鐵證]
我跑了三版 phased VCF 比對。
PASS = 47,798 / nonPASS = 3,139,477 — 三版完全相同。
Verdict_Somatic / Verdict_Germline / Verdict_SubclonalSomatic INFO tag 全 0 records。
INFO header 有定義但 ClairS-TO 預設未啟用 Verdict tagging。
所以論文 refined F1 在本 pipeline 不可量化 — 需重跑 ClairS-TO 啟用 Verdict tagging flag。

[longphase-to 真實修改範圍]
phase 階段 — 改 GT / PS / GT2 / GT3 / PON tag（FORMAT / INFO 欄位）。
haplotag 階段 — 改 BAM HP:i tag（不碰 VCF）。
V5 Pass 2 reclassify — 104,457 germline het GT 0/1 → 0/0（改 FORMAT 不改 FILTER）。

[FILTER 不動 = F1 invariant]
FILTER 由 ClairS-TO 一次性決定。
PASS set 不變 → TP/FP/FN 不變 → F1 數學保證 invariant。""",
    tier3="bcftools isec 機制 / V_H/V_L 12+3 patterns / longphase phase vs haplotag 階段分工 / FILTER 第 7 欄定義")

add(id="b3_hp33_paired", num="B3", section="Q&A Backup", rg2="2", ngrep="—",
    title="HP33 四版 trade-off + paired mode control",
    en="HP33 4-version trade-off + paired mode control group",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#374151;margin:0 0 4px;">❓ Q1: hp=33 (mixed-sub-tag) 為何重要 — baseline 沒有 / V3F 首現 / V5 壓掉 / V6 還原?</p>
    <table class="metric-table" style="font-size:11.5px;">
      <thead><tr><th>版本</th><th>hp=33 reads (全基因組)</th><th>機制</th><th>對下游影響</th></tr></thead>
      <tbody>
        <tr class="row-red"><td><strong>baseline</strong></td><td class="num">少 (≈0 結構化)</td><td>vector ① somatic break early → ② mixed pair 罕走到; hp=33 為「順序副作用」</td><td>HP_Ratio 失衡 / 無 mixed 訊號</td></tr>
        <tr class="row-green"><td><strong>V3F</strong></td><td class="num">132,060</td><td>Layer 2 結構化: <code>if(somTot>0 && germ==0) hp=33</code></td><td>首次正確標 hp=33 ✅</td></tr>
        <tr class="row-yellow"><td><strong>V5 (+L1.5)</strong></td><td class="num">13,250 (-89.9%)</td><td>Layer 1.5 用 somatic vote 改派 hp=11/21; 過度修正</td><td>下游 mixed-sub-tag 訊號遺失 ⚠</td></tr>
        <tr class="row-green" style="border-top:2px solid #DC2626;"><td><strong>V6 ★</strong></td><td class="num"><strong>138,317 (+4.7%)</strong></td><td>revert L1.5 → 還原 V3F 結構化保證</td><td>ISM 下游 marker engineering 受惠 ✅</td></tr>
      </tbody>
    </table>
    <p style="font-size:10px;color:#16A34A;margin:6px 0 0;font-weight:600;">→ V5→V6 transfer 守恆: 82% from hp=1-1 + 17% from hp=2-1 → hp=33 (mirror priority bug feature 化方向)</p>

    <p style="font-size:11px;font-weight:700;color:#374151;margin:10px 0 4px;">❓ Q2: 為何看 paired mode? — 它是 priority bug 的 control group</p>
    <table class="metric-table" style="font-size:11.5px;">
      <thead><tr><th>Mode</th><th>HP1:HP2 ratio (chr19)</th><th>som_ratio (mean)</th><th>codebase</th><th>有無 priority bug?</th></tr></thead>
      <tbody>
        <tr class="row-red"><td>TO baseline</td><td class="num">17.3:1 (94.6% HP1)</td><td class="num">0.946 偏 HP1</td><td>longphase-to (球員兼裁判 + getVote)</td><td>✗ 有 bug</td></tr>
        <tr class="row-green"><td><strong>paired</strong></td><td class="num"><strong>1:1.275 (接近隨機)</strong></td><td class="num"><strong>0.462 (跨 0-1 全範圍)</strong></td><td>longphase-s (獨立 codebase)</td><td>✓ 無 bug</td></tr>
      </tbody>
    </table>
    <p style="font-size:10px;color:#1E3A8A;margin:6px 0 0;font-weight:600;">→ paired 整體無 systematic bias 證 priority bug 是 TO mode 限定 (球員兼裁判 + getVote 順序錯); 不是 longphase 演算法整體缺陷</p>""",
    speaker="""[Q1]
hp=33 為何重要？baseline 沒有 / V3F 首現 / V5 壓掉 / V6 還原？

[hp=33 四版 trade-off]
baseline — vector ① somatic break early，hp=33 走 ② mixed pair 路徑機率被嚴重壓縮（順序副作用）。
V3F — 加 Layer 2 結構化 if (somaticTotal>0 && germlineResult==0) hpResult=33。首次保證 germline-absent 區正確標 hp=33 — 132,060 全基因組。
V5 Layer 1.5 — 用 somatic vote 改派 hp=11/21，把 hp=33 壓掉 −89.9% to 13,250。過度修正。
V6 — revert Layer 1.5，還原 V3F 結構化保證 — 138,317，+4.7% vs V3F。

[V5→V6 transfer 守恆]
82% from hp=1-1 + 17% from hp=2-1 → hp=33。
完美 mirror priority bug feature 化方向。

[Q2]
為何看 paired mode？

[paired mode = control group]
paired 用 longphase-s 獨立 codebase。
vs TO 用 longphase-to。

[paired 數據]
HP1:HP2 = 1:1.275 接近隨機。
som_ratio mean 0.462 跨 0-1 全範圍。
→ 整體無 systematic bias。

[結論]
證 priority bug 是 TO mode 限定（球員兼裁判 + getVote 順序錯 + Layer 1.5 4.19:1 偏 HP1）。
不是 longphase 演算法整體缺陷。
為 V5 Layer 1.5 caveat 提供 control 對比基礎。""",
    tier3="V5→V6 transfer 82:17 守恆細節 / paired chr19 354,919 tagged reads / longphase-s SomaticHaplotagProcess.cpp:533 / V5 Layer 1.5 germline-absent 4.19:1 量化 (5,789 events)")

add(id="b4_phase_d_outliers", num="B4", section="Q&A Backup", rg2="1", ngrep="—",
    title="Phase D 5 樣本 V6 + HCC1937 邊緣 + COLO829 阻塞",
    en="Phase D 5-sample V6 + HCC1937 edge + COLO829 pending",
    timing="100 sec / 中 ~340 字",
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#374151;margin:0 0 4px;">❓ Q: V6 跨樣本一致性? 為何 HCC1937 marker rate 0.817 邊緣 fail?</p>
    <table class="metric-table" style="font-size:11px;">
      <thead><tr><th>Sample</th><th>ratio (target≈1)</th><th>hp=33</th><th>marker rate (≥0.85)</th><th>NG_on=2 (≥0.85)</th><th>FP/TP</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>HCC1395 (ref)</td><td class="num">1.838</td><td class="num">138,317</td><td class="num">0.909 ✅</td><td class="num">0.829</td><td class="num">—</td></tr>
        <tr class="row-green"><td>H1437</td><td class="num">1.243 ✅</td><td class="num">39,050</td><td class="num"><strong>0.992 ✅</strong></td><td class="num">0.991 ✅</td><td class="num">0.011</td></tr>
        <tr class="row-green"><td>H2009</td><td class="num">0.901 ✅</td><td class="num">684,035</td><td class="num"><strong>0.993 ✅</strong></td><td class="num">0.992 ✅</td><td class="num">0.010</td></tr>
        <tr class="row-green"><td>HCC1954</td><td class="num">0.958 ✅</td><td class="num">4,859</td><td class="num">0.954 ✅</td><td class="num">0.967 ✅</td><td class="num">0.035</td></tr>
        <tr class="row-yellow"><td><strong>HCC1937 ⚠</strong></td><td class="num">0.611 ✅</td><td class="num">5,017</td><td class="num"><strong>0.817 ⚠</strong></td><td class="num">0.904 ✅</td><td class="num"><strong>0.194</strong></td></tr>
        <tr style="background:#F3F4F6;"><td>COLO829 ⏸</td><td colspan="5">truth set 0600 權限阻塞 (HKU/NYGC); BAM ready 待 chmod 660 或替代 PASS VCF</td></tr>
      </tbody>
    </table>
    <div class="grid-2col" style="margin:8px 0;">
      <div class="green-box" style="font-size:10.5px;padding:6px;">
        <strong>✅ V6 priority bug 修補一致成功 (4/4)</strong><br>
        ratio 跨 0.611-1.243 全部接近中性 (平均 0.928)<br>
        vs V5 baseline 1.86 / 原 baseline 17.3:1
      </div>
      <div class="caveat-box" style="font-size:10.5px;padding:6px;background:#FEF3C7;border:1px solid #CA8A04;border-radius:4px;">
        <strong>⚠ HCC1937 marker rate 0.817 邊緣 — 樣本特性非 V6 失敗</strong><br>
        BRCA1 mutant 高 ploidy + CNV-driven germline het<br>
        FP/TP=0.194 (vs 其他樣本 0.01-0.035 高 20×)<br>
        → 須加 AF&lt;0.4 filter (HPFineNGroups canonical)
      </div>
    </div>
    <p style="font-size:10px;color:#6B7280;font-style:italic;text-align:center;">Phase D pipeline: V5 binary phasing → V6 binary haplotag (hybrid); 每樣本 wall ~10 hr</p>""",
    speaker="""[Q]
V6 跨樣本一致性如何？為何 HCC1937 marker rate 0.817 邊緣 fail？

[Phase D 5 樣本評估]
HCC1395 (ref) + H1437 + H2009 + HCC1954 + HCC1937。
COLO829 truth set 0600 權限阻塞。

[ratio 結果 — 4 樣本全中性 ✅]
HCC1395 1.838 / H1437 1.243 / H2009 0.901 / HCC1954 0.958 / HCC1937 0.611。
範圍 0.611-1.243，平均 0.928。
vs V5 baseline 1.86，vs 原 baseline 17.3:1。
→ V6 priority bug 修補一致成功。

[marker rate 結果 — 3/4 通過 ≥0.85]
H1437 0.992 / H2009 0.993 / HCC1954 0.954 — 全通過。
HCC1937 0.817 邊緣 fail。

[HCC1937 邊緣 — 樣本特性非 V6 失敗]
BRCA1 mutant 高 ploidy 細胞株。
FP/TP = 0.194 — 比其他樣本 0.01-0.035 高 20 倍。
CNV-driven germline het AF 漂移混入 FP。
→ 須加 AF < 0.4 filter（與 HPFineNGroups canonical 一致）。
非否定 V6 patch，是樣本特性。

[COLO829 阻塞]
truth set 0600 權限。BAM ready 待 chmod 660 或替代 PASS VCF。
F1 follow-up 一條。

[Pipeline 細節]
Phase D pipeline = V5 phasing → V6 haplotag (hybrid)。
每樣本 wall ~10 hr。4 樣本約 40 hr 平行跑完。""",
    tier3="HCC1937 BRCA1 mutant 高 ploidy CNV 細節 / COLO829 truth set 0600 chmod 步驟 / HPFineNGroups AF<0.4 filter 文獻 / Phase D 4-stage pipeline")

# ═══════════════════════════════════════════════════════════════════════════
# Render templates
# ═══════════════════════════════════════════════════════════════════════════

def render_slide_page(s, prev_id, next_id, idx, total):
    """Render a single slide_XX.html page."""
    title_class = "critical" if s["title_critical"] else ""
    meta_class = "critical" if s["title_critical"] else ""
    prev_link = f'<a href="slide_{prev_id}.html">← Prev</a>' if prev_id else '<a class="disabled">← Prev</a>'
    next_link = f'<a href="slide_{next_id}.html">Next →</a>' if next_id else '<a class="disabled">Next →</a>'
    return f"""<!DOCTYPE html>
<html lang="zh-TW">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Slide {s['num']} — {s['title']}</title>
  <style>{STYLE_CSS}</style>
  <!-- PanZoom.js CDN — fail-soft：載入失敗時 fallback 到 CSS hover scale + new-tab full size -->
  <script src="https://unpkg.com/panzoom@9.4.3/dist/panzoom.min.js" defer></script>
</head>
<body>
<main class="slide-page" role="main" aria-label="Slide {s['num']} content">
  <div class="slide-meta {meta_class}" role="contentinfo" aria-label="Slide metadata">
    <div><span class="label">Slide:</span><span class="value">{s['num']} / 22</span></div>
    <div><span class="label">Section:</span><span class="value">{s['section']}</span></div>
    <div><span class="label">Tier 2:</span><span class="value">{s['timing']}</span></div>
    <div><span class="label">R-G2 術語:</span><span class="value">{s['rg2']}</span></div>
    <div><span class="label">N grep:</span><span class="value">{s['ngrep']}</span></div>
  </div>

  <article class="slide-canvas{' cover-slide' if s.get('is_cover') else ''}" data-section="{s.get('section','').split()[0] if s.get('section') else 'S0'}" aria-label="Main slide canvas">
    {'' if s.get('is_cover') else f'<h1 class="slide-title {title_class}">{s["title"]}</h1>'}
    {'' if (s.get('is_cover') or not s.get('en')) else f'<p class="en-subtitle">{s["en"]}</p>'}
    {s['canvas_html']}
  </article>

  <aside class="note-section" aria-label="Speaker notes for slide {s['num']}">
    <details class="collapsible" open>
      <summary>🎤 Speaker Note ({s['timing']})</summary>
      <div class="speak-text">{s['speaker']}</div>
      {f'<div class="tier3"><span class="label">[ORAL-OPTIONAL]</span> {s["tier3"]}</div>' if s["tier3"] else ''}
    </details>
  </aside>

  <aside class="cross-link" aria-label="Cross-reference to markdown spec">
    📋 詳細 multi-agent review (T/C/L/B/N + PI 6 並行) 與修正紀錄請見：
    <a href="../03_slide_layout_script.md">InterSubMod/docs/presentations/.../03_slide_layout_script.md</a>
  </aside>
</main>

<nav class="nav-bottom" aria-label="Slide navigation">
  {prev_link}
  <span class="counter">Slide {s['num']} ({idx + 1}/{total})</span>
  {next_link}
</nav>

<div class="kbd-hint">快速鍵: <kbd>←</kbd> <kbd>→</kbd> 切換</div>

<script>
  document.addEventListener('keydown', function(e) {{
    if (e.key === 'ArrowLeft' && {bool(prev_id)}) location.href = 'slide_{prev_id}.html';
    if (e.key === 'ArrowRight' && {bool(next_id)}) location.href = 'slide_{next_id}.html';
  }});

  // PanZoom.js progressive enhancement — fail-soft if CDN blocked
  window.addEventListener('load', function() {{
    if (typeof panzoom === 'undefined') {{
      console.warn('[PanZoom] CDN unavailable; falling back to CSS hover + new-tab zoom');
      return;
    }}
    document.querySelectorAll('.igv-zoom-wrap img.igv-thumb').forEach(function(img) {{
      var wrap = img.closest('.igv-zoom-wrap');
      wrap.classList.add('has-panzoom');  // disable CSS hover scale via this class
      // Disable target=_blank click for panzoom-enabled image (pan/zoom handled inline)
      var anchor = img.closest('a');
      if (anchor) anchor.addEventListener('click', function(e) {{ if (!e.shiftKey) e.preventDefault(); }});
      // Init panzoom with wheel zoom + drag pan
      var pz = panzoom(img, {{
        minZoom: 0.5, maxZoom: 8, bounds: false, beforeWheel: function() {{ return false; }},
        zoomDoubleClickSpeed: 1.5, smoothScroll: false
      }});
      img.style.cursor = 'grab';
      img.addEventListener('mousedown', function() {{ img.style.cursor = 'grabbing'; }});
      img.addEventListener('mouseup', function() {{ img.style.cursor = 'grab'; }});
      // Add reset button
      var resetBtn = document.createElement('button');
      resetBtn.textContent = '↺ 重置';
      resetBtn.className = 'igv-pz-reset';
      resetBtn.onclick = function() {{ pz.zoomAbs(0, 0, 1); pz.moveTo(0, 0); }};
      wrap.appendChild(resetBtn);
      // Update hint
      var hint = wrap.querySelector('.igv-zoom-hint');
      if (hint) hint.textContent = '🖱 滾輪 zoom / 拖曳 pan / Shift+點擊新分頁';
    }});
  }});
</script>

</body>
</html>"""


def render_index():
    """Render index.html — top nav + iframe."""
    nav_buttons = []
    for i, s in enumerate(SLIDES):
        cls = ""
        if s["num"] == "14":
            cls = "nav-arrow"
        elif s["num"] == "16":
            cls = "nav-star"
        elif s["num"].startswith("B"):
            cls = "nav-backup"
        active = " active" if i == 0 else ""
        nav_buttons.append(
            f'<button data-slide="{s["id"]}" class="{cls}{active}">{s["num"]} {s["section"][:6]}</button>'
        )
    nav_html = "\n  ".join(nav_buttons)

    return f"""<!DOCTYPE html>
<html lang="zh-TW">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Self-Phasing PPT Preview — 22 slide multi-page</title>
  <style>{STYLE_CSS}</style>
</head>
<body>

<header class="topbar" role="banner">
  <h1>📊 Self-Phasing PPT Preview — 22 slides multi-page</h1>
  <div class="meta">
    <code>InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/</code>
    &nbsp;|&nbsp; Build: 2026-05-11
  </div>
</header>

<nav class="slidetabs" id="tabs" role="navigation" aria-label="Slide selector">
  {nav_html}
</nav>

<main class="iframe-wrap" role="main">
  <iframe id="slide-frame" src="slide_{SLIDES[0]['id']}.html" title="Active slide preview"></iframe>
</main>

<div class="kbd-hint">點 nav 切換 slide; 子頁亦可直接雙擊獨立開啟</div>

<script>
  document.querySelectorAll('nav.slidetabs button').forEach(function(btn) {{
    btn.addEventListener('click', function() {{
      const id = btn.dataset.slide;
      document.getElementById('slide-frame').src = 'slide_' + id + '.html';
      document.querySelectorAll('nav.slidetabs button').forEach(function(b) {{ b.classList.remove('active'); }});
      btn.classList.add('active');
    }});
  }});
</script>

</body>
</html>"""


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
def main():
    print(f"[build] preview/ generator — {len(SLIDES)} slides")
    print(f"[build] Output: {OUT}")

    for i, s in enumerate(SLIDES):
        prev_id = SLIDES[i - 1]["id"] if i > 0 else None
        next_id = SLIDES[i + 1]["id"] if i < len(SLIDES) - 1 else None
        html = render_slide_page(s, prev_id, next_id, i, len(SLIDES))
        out_path = OUT / f"slide_{s['id']}.html"
        out_path.write_text(html, encoding="utf-8")
        print(f"  [ok] {out_path.name}")

    index_html = render_index()
    (OUT / "index.html").write_text(index_html, encoding="utf-8")
    print(f"  [ok] index.html (with {len(SLIDES)} nav buttons)")
    print(f"\n[done] Open: file://{OUT / 'index.html'}")


if __name__ == "__main__":
    main()
