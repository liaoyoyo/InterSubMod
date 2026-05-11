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
    speaker="謝謝各位。今天的報告主題是 self-phasing 問題的修正與驗證 — 在 longphase-to 工具的 tag 層發現一個設計缺陷，整理修正過程與跨樣本驗證的結果。報告人廖子游，這是 5/1 至 5/11 的階段性整理。",
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

    <div style="margin-top:32px;padding:14px 20px;background:#F0F9FF;border:1px solid #BAE6FD;border-radius:8px;color:#0C4A6E;font-size:18px;font-weight:600;text-align:center;">
      18 章節 slide + 3 backup（Q&A）  ·  預計 20 分鐘 + 5 分鐘問答
    </div>""",
    speaker="今天報告分三個階段六個章節。問題發現階段：先看觀察起點的偏移現象，再解析兩層機制。修正過程階段：用 read-level 量化鐵證驗證，並設計三版漸進修補。驗證結果階段：跨指標跨樣本一致性驗證，最後給結論與後續方向。預計 20 分鐘 + 5 分鐘問答。",
    tier3="(目錄 speaker note 精簡至 60 sec，純預告章節，不透露具體結論)")

# ─── S1 觀察起點 ───────────────────────────────────────────────────────────
add(id="03_genome_173", num="03", section="S1 觀察起點", rg2="3", ngrep="8",
    title="全基因組 HP1:HP2 = 17.3:1 是 systematic bias，非樣本性質",
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
      ①  生物學: tumor sub-clone 跨 23 染色體不該系統偏 HP1<br>
      ②  跨 chr 一致: cnLOH artifact 只影響單一 chr; 94.6% 跨 chr 一致<br>
      ③  paired 對照: paired tumor-normal 同 reads HP1:HP2 ≈ 1:1
    </div>
    <div class="conclusion-arrow">→ 17.3:1 是 LongPhase-TO 的 systematic engineering artifact</div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">germline het:</span> 雜合位點</div>
      <div class="gloss-item">📖 <span class="term">sub-clone:</span> 基因型同細胞群</div>
      <div class="gloss-item">📖 <span class="term">haplotype:</span> 父系/母系染色體</div>
    </div>""",
    speaker="baseline LongPhase-TO 全基因組: HP1 reads 614K vs HP2 35.5K, 比例 17.3:1 vs 隨機 1:1。94.6% 占比是 systematic bias 的硬證據, 三條獨立論證: 生物學 (tumor sub-clone 不該跨 23 chr 系統偏 HP1)、跨 chr 一致 (cnLOH 只影響單 chr 但這偏移跨 23 chr 一致)、paired 對照 (paired pipeline HP1:HP2 ≈ 1:1)。三條互相獨立 → engineering artifact。",
    tier3="cnLOH 機制細節 / 23 chr 一致性表 / paired pipeline 程式碼差異")

add(id="04a_sp1", num="04a", section="S1 觀察起點", rg2="1", ngrep="3",
    title="SP1 個案 — baseline 113:0 → V5 翻轉對齊 paired",
    en="SP1 chr19:17,565,944 — baseline 113:0 → V5 flips to HP2",
    timing="60 sec / 中 ~250 字",
    canvas_html="""
    <div class="grid-3col" style="margin-bottom:12px;">
      <div class="stat-box"><div class="number orange">113 : 0</div><div class="label">baseline HP1:HP2</div></div>
      <div class="stat-box"><div class="number">HP2 主導</div><div class="label">V5 修正後</div></div>
      <div class="stat-box"><div class="number green">✅ 對齊</div><div class="label">paired ground truth</div></div>
    </div>
    <img class="igv-thumb" src="../figures/igv/D_SP1_chr19_17565944.png" alt="IGV SP1" style="max-height:400px;">
    <p class="fig-caption" style="font-size:13px;">chr19:17,565,944 · 6-BAM 並列: baseline / V2b / V3F / V5 / paired_T / paired_N</p>
    <div class="arg-list">
      <strong>為何能排除噪音 / caller / alignment?</strong><br>
      ① baseline 與 paired 方向相反 — 翻轉而非衰減<br>
      ② V5 修正後與 paired ground truth 重合<br>
      → read assignment 強制集中的鐵證
    </div>""",
    speaker="全基因組 17.3:1 是平均值; IGV 6-BAM 並列篩到 chr19 三個近 100% 失衡位點。SP1 chr19:17,565,944: baseline 113 reads 全 HP1, HP2=0; V5 翻 HP2 與 paired tumor 一致。排除噪音/caller/alignment: baseline 與 paired 完全反向 (翻轉而非衰減), V5 與 paired 重合 → read assignment 強制集中的鐵證。",
    tier3="6-BAM 並列順序 / V2b 中間階段意義")

add(id="04b_sp2_sp3", num="04b", section="S1 觀察起點", rg2="0", ngrep="6",
    title="SP2 + SP3 同模式 — 3/3 對齊 paired",
    en="SP2 + SP3 follow SP1 pattern, 3/3 aligned with paired",
    timing="60 sec / 中 ~260 字",
    canvas_html="""
    <div class="grid-2col">
      <div>
        <table class="metric-table"><thead><tr><th>SP2 · chr19:12,452,332</th><th>baseline</th><th>V5</th></tr></thead>
        <tbody><tr class="row-red"><td><strong>HP1 : HP2</strong></td><td class="num">109 : 1</td><td>HP2 主導</td></tr></tbody></table>
        <img class="igv-thumb" src="../figures/igv/D_SP2_chr19_12452332.png" alt="IGV SP2" style="max-height:280px;">
      </div>
      <div>
        <table class="metric-table"><thead><tr><th>SP3 · chr19:12,467,180</th><th>baseline</th><th>V5</th></tr></thead>
        <tbody><tr class="row-red"><td><strong>HP1 : HP2</strong></td><td class="num">108 : 0</td><td>HP2 主導</td></tr></tbody></table>
        <img class="igv-thumb" src="../figures/igv/D_SP3_chr19_12467180.png" alt="IGV SP3" style="max-height:280px;">
      </div>
    </div>
    <div class="conclusion-arrow green">→ 三 SP 都在 chr19:12-17M 區段 → 對齊 slide 09 chr19 752 victims hotspot</div>""",
    speaker="SP2 chr19:12,452,332 baseline 109:1; SP3 chr19:12,467,180 baseline 108:0; 與 SP1 同模式: baseline 全 HP1, V5 翻 HP2 對齊 paired 3/3。三 SP 都在 chr19:12-17M, 對齊 chr19 752 victims hotspot — read-level 個案與 IGV 屬同機制不同層級。 [口述過渡] 接下來四個章節將分別回答: ① 為何全集中一邊? (S2 機制) ② read 層級? (S3 量化) ③ 三版各修? (S4 修補) ④ 是否都修對? (S5 驗證)。",
    tier3="paired_T 與 paired_N 對照細節 / 三 SP 完整座標表")

# ─── S2 機制 ──────────────────────────────────────────────────────────────
add(id="05_player_referee", num="05", section="S2 機制", rg2="4 (核心 forced)", ngrep="3 + 1 commit",
    title="phasing 層球員兼裁判 — somatic 100% 共現蓋過 germline 50/50",
    en="Phasing layer player-as-referee — somatic overrules germline",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <div class="grid-2col" style="grid-template-columns: 7fr 5fr;">
      <img class="fig-thumb" src="../figures/G1_player_as_referee.png" alt="G1">
      <div>
        <p style="font-size:12px;font-weight:700;color:#7F1D1D;margin:0 0 4px;">TO 模式致命:</p>
        <ul style="margin:4px 0 0 16px;font-size:11px;line-height:1.5;color:#374151;background:#FEE2E2;border:1px solid #FCA5A5;padding:8px 12px 8px 24px;border-radius:4px;">
          <li>TO 沒 paired normal → 用 PoN (somatic 不在 PoN → 當 germline)</li>
          <li>somatic 進 graph 後 edge weight 暴漲 (100% > 50%)</li>
          <li>自我增強: 3 somatic 共現 → 強度 ×3</li>
          <li>germline 真實訊號被 overrule</li>
        </ul>
        <div class="conclusion-arrow green" style="margin-top:8px;font-size:12px;">
          解法: PON-only flag<br>(commit <code>8b8c1fd</code>)
        </div>
      </div>
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">PoN:</span> Panel of Normals — 多正常樣本建構的 germline reference</div>
      <div class="gloss-item">📖 <span class="term">phasing graph:</span> het 位點當 node, read 共現當 edge</div>
    </div>""",
    speaker="機制兩層 bug, 先 phasing 層球員兼裁判隱喻。phasing graph 物理基礎: germline het 在 HP1/HP2 50/50 隨機分佈。TO 沒 paired normal 區分 germline/somatic, 只能用 PoN; 未在 PoN 內的位點被當 germline → somatic 進 graph 後 edge weight 暴漲, graph 自我增強。球員兼裁判: somatic 應被 phase 反過來主導 graph。修法 PON-only flag (8b8c1fd)。但只解 phasing 層, tag 還壞。",
    tier3="TO vs paired PoN 對照 / 自我增強迴圈 / Pass 1 vs Pass 2 預告")

add(id="06_priority_bug", num="06", section="S2 機制", rg2="1 + 2 footnote", ngrep="5",
    title="tagging 層 getVote priority bug — 1 票 somatic 觸發誤標",
    en="getVote priority bug — 1 somatic vote triggers mislabel",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <div class="grid-2col" style="grid-template-columns: 7fr 5fr;">
      <img class="fig-thumb" src="../figures/master/F1_priority_bug_mechanism.png" alt="F1" style="max-height:300px;">
      <div>
        <p style="font-size:11px;font-weight:700;color:#7F1D1D;margin:0 0 4px;">Real read 範例 (752 同模式):</p>
        <pre class="code-panel baseline" style="font-size:10px;">germline HP1 = 0
germline HP2 = 5  ← 主導
somatic HP1_1 = 1 ← 1 票觸發
somatic HP2_1 = 0
─────────────────
baseline: → hp=11 ❌ break
        (germline 5 票被忽略)
正確答案: hp=21
(germline HP2=5 主導 + 標 21)</pre>
      </div>
    </div>
    <div class="conclusion-arrow">→ tumor sub-clone somatic 100% 同方向 → priority bug 把所有受影響 reads 標 HP:i:11 系列 → 17.3:1 偏移在 tag layer 形成</div>
    <p style="font-size:11px;color:#6B7280;font-style:italic;margin:4px 0;">注意: tag layer 與 §slide 05 phasing layer 是不同層 bug, 必須分別修補</p>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">sub-clone:</span> 基因型同細胞群</div>
      <div class="gloss-item">ⓘ priority bug / break early</div>
    </div>""",
    speaker="tagging 層 priority bug: baseline getVote() vector 順序 ① somatic ② mixed ③ germline; for 迴圈第一個非空 break early, germline 永遠看不到。範例 read: germline HP2=5 主導, somatic HP1_1=1 觸發, baseline 標 hp=11 錯, 應 hp=21。全 752 chr19 victims 同模式。tumor sub-clone somatic 100% 同方向 → 17.3:1 形成。tag layer 與 phasing layer 不同層, 不能合併修。",
    tier3="enum HAPLOTYPE1_1=2 vs HP tag int=11 / 5-vote countMap")

add(id="07_two_layer_table", num="07", section="S2 機制", rg2="1", ngrep="5 commits + 1",
    title="兩層 bug 兩層修補對應 — phasing PON-only + tagging V3F+V5 缺一不可",
    en="Two-layer bugs: PON-only + V3F+V5 both required",
    timing="90 sec / 中 ~310 字",
    canvas_html="""
    <table class="metric-table" style="font-size:13px;">
      <thead><tr><th>Layer</th><th>Bug</th><th>修補 commit</th><th>章節</th></tr></thead>
      <tbody>
        <tr style="background:#DBEAFE;"><td><strong>phasing 層</strong></td><td>球員兼裁判</td><td><code>8b8c1fd</code> PON-only</td><td>§5.2</td></tr>
        <tr style="background:#DCFCE7;"><td><strong>tagging 層</strong></td><td>priority bug</td><td>V3F: <code>41ff147</code> + <code>380e8d2</code></td><td>§5.3</td></tr>
        <tr style="background:#DCFCE7;"><td></td><td></td><td>V5: + <code>d0bcd8c</code> + <code>938f0df</code></td><td>§5.4</td></tr>
      </tbody>
    </table>
    <div class="grid-2col">
      <div class="caveat-box">
        <span class="label">為何不能只用 PON-only?</span><br>
        解 phasing 但 tag 仍壞 → 99.9% reads 仍標 HP:i:11
      </div>
      <div class="caveat-box">
        <span class="label">為何 V3F 不夠還要 V5?</span><br>
        V3F Layer 1 only → germline 缺席區 reads 全 untagged; V5 Layer 1.5 補 fallback (slide 16 詳述 caveat)
      </div>
    </div>
    <div class="footer-glossary"><div class="gloss-item">📖 <span class="term">PoN:</span> Panel of Normals</div></div>""",
    speaker="兩層 bug 獨立: phasing 層 → 8b8c1fd PON-only; tagging 層 → V3F two-layer (41ff147 + 380e8d2 INDEL guard) → V5 (+ d0bcd8c + 938f0df)。任一單修不夠必須 stacking。PON-only 不夠: tag 仍偏 99.9% 標 11; V3F 不夠: germline 缺席區 untagged; V5 補 Layer 1.5 但 germline-absent 有 caveat (slide 16)。",
    tier3="5 commit 順序 / V3F 命名 / 跨層交互")

# ─── S3 量化鐵證 ───────────────────────────────────────────────────────────
add(id="08_chr19_752", num="08", section="S3 量化鐵證", rg2="—", ngrep="6",
    title="chr19 752 read-level victims — 100% 單向 baseline=11 → V3F/V5=21",
    en="chr19 752 victims — 100% unidirectional fix",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <div class="grid-2col">
      <div>
        <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:0 0 4px;">規模:</p>
        <pre class="code-panel" style="background:#F9FAFB;border:1px solid #E5E7EB;font-size:11px;">Dump rows:        549,206
3-way merged:     1,069,832
Priority bug:        752</pre>
      </div>
      <div>
        <p style="font-size:12px;font-weight:700;color:#16A34A;margin:0 0 4px;">修正率:</p>
        <pre class="code-panel" style="background:#DCFCE7;border:1px solid #16A34A;font-size:11px;color:#166534;">V3F 修正率: 100.00% ✅
V5  修正率: 100.00% ✅
全 752 條無一條反向</pre>
      </div>
    </div>
    <div class="grid-2col" style="grid-template-columns: 7fr 5fr;">
      <table class="metric-table" style="font-size:11px;">
        <thead><tr><th>路徑</th><th>結果</th><th>判定</th></tr></thead>
        <tbody>
          <tr class="row-green"><td>① 個案 trace</td><td>752 條</td><td>✅</td></tr>
          <tr class="row-yellow"><td>② 1Mb 區域聚集</td><td>30M+27M 46%</td><td>⚠ PARTIAL</td></tr>
          <tr><td>③ Density 共變</td><td>high≥5=0; low 觸發</td><td>🔄 反向有意義</td></tr>
          <tr class="row-green"><td>④ 修正後消失</td><td>V3F/V5 100%</td><td>✅</td></tr>
        </tbody>
      </table>
      <img class="fig-thumb" src="../figures/master/F4_chr19_752_victims_scatter.png" alt="F4" style="max-height:200px;">
    </div>
    <div class="conclusion-arrow">→ 3 PASS + 1 PARTIAL = priority bug 機制因果確立 (chr19 only scope)</div>
    <div class="footer-glossary"><div class="gloss-item">ⓘ scope: chr19 only</div></div>""",
    speaker="chr19 read-level audit: dump 549,206 rows × 3 binary versions JOIN merged 1.07M events, 篩 germline_majority ≠ somatic_majority = 752 victims。V3F+V5 修正率 100%, 全 752 條 baseline=11 → V3F=21 單向。4-path 驗證 3 PASS + 1 PARTIAL: 個案 PASS、1Mb hotspot PARTIAL、density 反向有意義 (high vote sub-clone 一致已對齊, low 觸發)、修正後消失 PASS。機制因果確立。",
    tier3="4-path detail / read_name 真實 case / SP1/2/3 對應 hotspot")

add(id="09_genome_34855", num="09", section="S3 量化鐵證", rg2="—", ngrep="14",
    title="全基因組 34,855 victims (46×) — chr19 占 2.16% rank 19",
    en="Genome-wide 34,855; main hotspots chr7/chr2/chr1",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <table class="metric-table" style="font-size:11px;">
      <thead><tr><th></th><th>chr19 pilot</th><th>Genome F1</th><th>倍數</th></tr></thead>
      <tbody>
        <tr><td>Dump rows</td><td class="num">549,206</td><td class="num">29,973,253</td><td class="num">54.6×</td></tr>
        <tr class="row-yellow"><td><strong>Priority bug victims</strong></td><td class="num">752</td><td class="num">34,855</td><td class="num">46.4×</td></tr>
        <tr class="row-green"><td>V3F / V5 修正率</td><td class="num">100% / 100%</td><td class="num">100% / 100%</td><td>一致 ✅</td></tr>
      </tbody>
    </table>
    <div class="grid-2col" style="grid-template-columns: 7fr 5fr;">
      <table class="metric-table" style="font-size:11px;">
        <thead><tr><th>chr</th><th>victims</th><th>占比</th><th>rank</th></tr></thead>
        <tbody>
          <tr><td>chr7</td><td class="num">3,508</td><td class="num">10.1%</td><td class="num">1</td></tr>
          <tr><td>chr2</td><td class="num">2,792</td><td class="num">8.0%</td><td class="num">2</td></tr>
          <tr><td>chr1</td><td class="num">2,674</td><td class="num">7.7%</td><td class="num">3</td></tr>
          <tr class="row-red"><td><strong>chr19</strong></td><td class="num">752</td><td class="num">2.16%</td><td class="num">19 ★</td></tr>
          <tr style="background:#DBEAFE;"><td>chr8</td><td class="num">666</td><td class="num">1.9%</td><td class="num">21 冷區</td></tr>
        </tbody>
      </table>
      <img class="fig-thumb" src="../figures/master/F2_priority_bug_per_chr_enrichment.png" alt="F2" style="max-height:240px;">
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">ⓘ scope: 全基因組 (T1.2-F1)</div>
      <div class="gloss-item">ⓘ chr8 LOH+HPSig 是 ISM 下游 hotspot</div>
    </div>""",
    speaker="全基因組擴展 34,855 victims (46.4×)。V3F/V5 修正率仍 100%。per-chr 分佈推翻原 chr19 結論: 主要 hotspot chr7/chr2/chr1/chr16/chr20, chr19 占 2.16% rank 19。chr8 priority bug 0.34× (rank 21 冷區), 與 chr8 LOH+HPSig hotspot 是不同 layer。",
    tier3="全 chr enrichment ‰ 表 / chrY 小 N 高 ‰ / chr8 layer 區分")

# ─── S4 修補設計 ───────────────────────────────────────────────────────────
add(id="10_5_commits", num="10", section="S4 修補設計", rg2="0 + 2 footnote", ngrep="5 hash",
    title="5 commits 兩層三版 stacking — baseline → V3F → V5",
    en="5 commits two-layer three-version stacking",
    timing="120 sec / 中 ~340 字",
    canvas_html="""
    <img class="fig-thumb" src="../figures/master/F3_binary_commit_timeline.png" alt="F3" style="max-height:280px;">
    <p class="fig-caption">F3: phasing 藍 / tagging 綠 / 跨層紫</p>
    <div class="grid-2col">
      <div class="green-box">
        <strong>V3-Fixed = baseline + 41ff147 + 380e8d2</strong><br>
        ★ 41ff147 是修偏移的關鍵 commit
      </div>
      <div class="green-box">
        <strong>V5 = V3F + d0bcd8c + 938f0df</strong><br>
        d0bcd8c 是唯一跨兩層 commit
      </div>
    </div>
    <p style="font-size:11px;color:#6B7280;text-align:center;">累計 ~155 lines tagging + ~40 lines phasing; <code>HaplotagProcess.h:66-68</code> 介面契約零變動</p>
    <div class="footer-glossary">
      <div class="gloss-item">ⓘ INDEL guard: HAPLOTYPE_UNDEFINED 檢查</div>
      <div class="gloss-item">ⓘ threshold 0.95→0.9: Pass 2 觸發</div>
    </div>""",
    speaker="5 commits 漸進完成 self-phasing 修補。8b8c1fd PON-only (藍); 41ff147 two-layer getVote (綠) ★ 修偏移關鍵; 380e8d2 INDEL guard (綠); d0bcd8c Pass 2 ploidy fix + bundled Layer 1.5 (紫 跨兩層); 938f0df threshold 0.95→0.9 (藍)。V3F = baseline + 41ff147 + 380e8d2; V5 = V3F + d0bcd8c + 938f0df。累計 155+40 行; 介面契約零變動。",
    tier3="各 commit layer / 為何不合併 / cherry-pick 自 zhenyu")

add(id="11_getvote", num="11", section="S4 修補設計", rg2="4 (核心 forced)", ngrep="—",
    title="getVote 四版演進 — baseline → V3F → V5 +Layer 1.5 → V6 (Layer 1.5 移除)",
    en="getVote 4-version evolution: V6 reverts V5 Layer 1.5",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <img class="fig-thumb" src="../figures/G3_getVote_three_layer.png" alt="G3" style="max-height:280px;">
    <div class="green-box" style="font-size:11px;font-weight:600;">
      ★ V6 (5/10 binary patch): 移除 V5 Layer 1.5 邏輯 (HaplotagProcess.cpp:537-548)
      <pre style="font-family:monospace;font-size:9px;background:white;padding:4px 8px;margin:4px 0;border-radius:3px;border:1px solid #16A34A;">- else if (somaticHP1 > 0 || somaticHP2 > 0) {<br>-     germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;<br>- }<br>+ // V6: germline absent → conservative hp=33 (V3F behavior)</pre>
      → germline-absent 退回 V3F hp=33 行為; phasing 層完全不變 (caller F1 不變)
    </div>
    <p style="font-size:11px;color:#6B7280;text-align:center;">V6 = V3F germline-absent 保守 + V5 germline-existent 設計目標 (ploidy/threshold) = hybrid 升級</p>""",
    speaker="四版 code 演進。baseline 紅底 vector 順序錯 break early。V3F 綠底 explicit Layer 1 + Layer 2。V5 綠底加黃 highlight 新增 Layer 1.5 fallback (germline 缺席用 somatic phased votes) — 但 5/9 發現 Layer 1.5 在 germline-absent 區繼承 priority bug。★ V6 (5/10 binary patch): 移除 Layer 1.5 邏輯 (HaplotagProcess.cpp:537-548 那段 else if 整段去掉), germline-absent 退回 hp=33 (V3F 行為), 但保留 V5 phasing 層 — V6 重用 V5 phased VCF 故 caller F1 完全不變。V6 = V3F germline-absent 保守 + V5 germline-existent 設計目標 (ploidy/threshold) = hybrid 升級。",
    tier3="V6 patch diff / HaplotagProcess.cpp:537-548 移除 13 行 / phasing 層完全不變 / V6 重用 V5 phased VCF")

add(id="12_sp_fixed", num="12", section="S5 驗證", rg2="1", ngrep="11",
    title="個案層 V5 修正 3/3 + 全基因組 HP1:HP2 17.3:1 → ~1:1",
    en="Site-level V5 fixes 3/3 + genome 17.3:1 → ~1:1",
    timing="90 sec / 中 ~340 字",
    canvas_html="""
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:0;">個案層: SP1/2/3 修正後對齊 paired</p>
    <table class="metric-table" style="font-size:11px;">
      <thead><tr><th>位點</th><th>baseline</th><th>V5</th><th>paired</th><th>對齊?</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>SP1 chr19:17,565,944</td><td class="num">113:0</td><td>HP2 主導</td><td>HP2</td><td>✅</td></tr>
        <tr class="row-green"><td>SP2 chr19:12,452,332</td><td class="num">109:1</td><td>HP2 主導</td><td>HP2</td><td>✅</td></tr>
        <tr class="row-green"><td>SP3 chr19:12,467,180</td><td class="num">108:0</td><td>HP2 主導</td><td>HP2</td><td>✅</td></tr>
      </tbody>
    </table>
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:8px 0 4px;">全基因組層:</p>
    <table class="metric-table" style="font-size:11px;">
      <thead><tr><th>指標</th><th>baseline</th><th>V5</th><th>Δ</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>HP1:HP2 ratio</td><td>17.3:1</td><td>~1:1</td><td>消除偏移</td></tr>
        <tr class="row-green"><td>94.6% somatic→HP1</td><td>是</td><td>~50%</td><td>balanced</td></tr>
        <tr class="row-green"><td>15-site Problem PS</td><td class="num">48.5%</td><td class="num">52.0%</td><td class="num">+3.5 pp</td></tr>
      </tbody>
    </table>
    <div class="grid-2col" style="grid-template-columns: 7fr 5fr;">
      <img class="fig-thumb" src="../figures/master/F5_layer15_zero_sum_4quadrant.png" alt="F5" style="max-height:140px;">
      <div class="green-box" style="font-size:11px;">
        F5 zero-sum:<br>
        germline=0  +560,881<br>
        germline>0  −560,881<br>
        總和         =0
      </div>
    </div>""",
    speaker="個案層: SP1/2/3 V5 翻 HP2 對齊 paired 3/3。全基因組: 17.3:1 → ~1:1; 94.6% → ~50%; 15-site Problem PS 48.5%→52.0% +3.5 pp。F5 zero-sum: germline=0 +560K / germline>0 -560K / 總和 0。Pass 2 reclassify 104K germline het 為 somatic/未 phase。",
    tier3="V2b/V3F 中間版本 / Layer 1.5 zero-sum")

add(id="13_20_metrics", num="13", section="S5 驗證", rg2="1", ngrep="20+",
    title="20 指標 no regression — 6 項 ⭐ 顯著改善 +8.3 ~ +99.7%",
    en="20 metrics no regression; 6 significant improvements",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <table class="metric-table" style="font-size:11px;">
      <thead><tr><th>類別</th><th>內容</th></tr></thead>
      <tbody>
        <tr class="row-green"><td><strong>① ISM aggregate (3)</strong></td><td>TP_rate +0.005 / HP_Ratio 0.788→0.574 / Potential_LOH +3.5 pp</td></tr>
        <tr class="row-green"><td><strong>② HP_Ratio AUC (2)</strong></td><td>All -0.005 (隨機區間) / Inner +0.002</td></tr>
        <tr class="row-green"><td><strong>③ Methylation 6 feat</strong></td><td>全 ±0.01 內持平</td></tr>
        <tr class="row-yellow"><td><strong>④ Paired GT concord. ⭐ (4)</strong></td><td>clean PS +8.3 pp / 15-Aggr +6.65 pp / 15-Clean PS +13.3 pp</td></tr>
        <tr class="row-yellow"><td><strong>⑤ HP / LOH 結構 ⭐ (5)</strong></td><td>N50 +99.7% / Phased +23.6 pp / 1.36× 快 / LOH 完全相同</td></tr>
      </tbody>
    </table>
    <div class="green-box" style="text-align:center;font-weight:700;">6 顯著改善 ⭐: N50 +99.7% / Phased +23.6 pp / 1.36× / 15-Clean PS +13.3 pp / clean PS +8.3 pp / 15-Aggr +6.65 pp</div>
    <div class="conclusion-arrow">→ 20/0 指標 no regression — V5 全面 production-ready</div>
    <p style="font-size:10px;color:#6B7280;font-style:italic;text-align:center;">(HP_Ratio 0.788→0.574 是 tag bias 修正非變差 — pre-registered metrics)</p>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">LOH:</span> 雜合性丟失</div>
      <div class="gloss-item">ⓘ scope: HCC1395 5kHz @ 0.93 (PI 報告 V5 = Pass 1 only)</div>
    </div>""",
    speaker="5 大類 20 指標全綠: ① ISM 3 項 (HP_Ratio 0.788→0.574 是 tag bias 修正非變差); ② HP_Ratio AUC 隨機區間; ③ methylation 6 持平; ④ Paired GT 4 ⭐ +6.65~+13.3 pp; ⑤ HP/LOH 5 ⭐ N50 +99.7% / Phased +23.6 pp / 1.36× / LOH Jaccard=1.0。20/0 no regression production-ready。",
    tier3="methylation 6 列表 / HP_Ratio 詳解 / LOH Jaccard=1.0")

add(id="14_caller_f1", num="14", section="S5 驗證 + ⚡ Cliffhanger", rg2="2", ngrep="15+",
    title="Caller F1 vs SEQC2 三版完全相同; purity 0.6 完整對照 0 critical regression",
    en="Caller F1 identical; purity 0.6 fully verified",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:0;">HCC1395 5kHz @ 0.93 purity:</p>
    <table class="metric-table" style="font-size:10px;">
      <thead><tr><th>版本</th><th>TP</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>A1 baseline</td><td class="num">28,509</td><td class="num">11,606</td><td class="num">10,938</td><td class="num">0.7107</td><td class="num">0.7227</td><td class="num"><strong>0.7166</strong></td></tr>
        <tr class="row-green"><td>A3 V3F</td><td class="num">28,509</td><td class="num">11,606</td><td class="num">10,938</td><td class="num">0.7107</td><td class="num">0.7227</td><td class="num"><strong>0.7166</strong></td></tr>
        <tr class="row-green"><td>A5 V5</td><td class="num">28,509</td><td class="num">11,606</td><td class="num">10,938</td><td class="num">0.7107</td><td class="num">0.7227</td><td class="num"><strong>0.7166</strong></td></tr>
      </tbody>
    </table>
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:6px 0 0;">HCC1395 t30_n20 @ 0.6 purity:</p>
    <table class="metric-table" style="font-size:10px;">
      <thead><tr><th>版本</th><th>TP</th><th>FP</th><th>FN</th><th>F1</th></tr></thead>
      <tbody><tr class="row-green"><td>B1/B3/B5 三版</td><td class="num">24,190</td><td class="num">13,487</td><td class="num">15,257</td><td class="num"><strong>0.6273</strong></td></tr></tbody>
    </table>
    <div class="arg-list" style="font-size:11px;">
      <strong>因果鏈:</strong> ClairS-TO PASS set 由 caller 決定 → V5 改 GT/PS 不改 FILTER → PASS set 不變 → TP/FP/FN 不變 → F1 不變
    </div>
    <div class="cliffhanger-box">
      → V5 (與 V6) 都不改 caller; ΔF1 (0.93→0.6) = -0.0893 為 ClairS-TO 性質<br>
      <strong>→ 5/9 paired cross-ref 揭 V5 Layer 1.5 缺陷 → 5/10 V6 binary patch 已修補...</strong> 〔next: slide 15-16〕
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">purity:</span> 樣本中腫瘤細胞占比</div>
      <div class="gloss-item">ⓘ <span class="term">PASS set:</span> ClairS-TO snv.vcf FILTER=PASS 集合</div>
    </div>""",
    speaker="Caller F1 vs SEQC2: 0.93 三版 F1=0.7166 完全相同; 0.6 三版 F1=0.6273 完全相同。因果鏈: ClairS-TO PASS set 由 caller 決定, V5 改 GT/PS/GT2/GT3 不改 FILTER → PASS set 不變 → F1 不變。V5 不改 caller, ΔF1 -0.0893 是 ClairS-TO 在低 purity 性質。turning point — 但 5/9 paired cross-ref 揭露另一面。",
    tier3="PASS set / FILTER 機制 / purity 0.6 N50 微差")

add(id="15_paired_mode", num="15", section="S6 5/9 新發現", rg2="1 + 1 footnote", ngrep="11",
    title="paired mode 整體無偏移 — HP1:HP2 = 1:1.275; som_ratio mean 0.462",
    en="paired mode no systematic bias",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <div class="grid-2col">
      <div>
        <p style="font-size:11px;font-weight:700;color:#1E3A8A;margin:0;">paired chr19 HP:Z: 分布 (354,919 tagged):</p>
        <table class="metric-table" style="font-size:10px;">
          <thead><tr><th>HP:Z:</th><th>reads</th><th>%</th></tr></thead>
          <tbody>
            <tr><td>HP:Z:2</td><td class="num">183,309</td><td class="num">51.6%</td></tr>
            <tr><td>HP:Z:1</td><td class="num">143,760</td><td class="num">40.5%</td></tr>
            <tr><td>HP:Z:2-1</td><td class="num">14,504</td><td class="num">4.1%</td></tr>
            <tr><td>HP:Z:1-1</td><td class="num">12,401</td><td class="num">3.5%</td></tr>
            <tr><td>HP:Z:3</td><td class="num">1,145</td><td class="num">0.3%</td></tr>
          </tbody>
        </table>
      </div>
      <div>
        <p style="font-size:11px;font-weight:700;color:#1E3A8A;margin:0;">paired vs TO ratio:</p>
        <table class="metric-table" style="font-size:10px;">
          <thead><tr><th></th><th>paired</th><th>TO baseline</th></tr></thead>
          <tbody>
            <tr class="row-green"><td>germline</td><td>1:1.275 ✅</td><td>17.3:1 ❌</td></tr>
            <tr class="row-green"><td>somatic</td><td>1:1.169 ✅</td><td>全偏 HP1</td></tr>
          </tbody>
        </table>
        <p style="font-size:10px;color:#6B7280;margin-top:6px;">57 chr19 1Mb windows: som_ratio mean 0.462 / median 0.494 / stdev 0.332</p>
      </div>
    </div>
    <img class="fig-thumb" src="../figures/master/F6_paired_vs_TO_HP_distribution.png" alt="F6" style="max-height:180px;">
    <p style="font-size:10px;color:#9CA3AF;text-align:center;font-style:italic;">paired = longphase-s (獨立 codebase); HP tag 用字串 HP:Z:</p>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">HP:i: vs HP:Z::</span> longphase-to=整數, longphase-s=字串</div>
      <div class="gloss-item">ⓘ som_ratio: HP1-1/(HP1-1+HP2-1)</div>
    </div>""",
    speaker="paired mode 用 longphase-s (獨立 codebase), HP:Z: 字串。chr19 paired 分布 germline 1:1.275 接近隨機, somatic 1:1.169。對比 TO baseline 17.3:1 — paired 整體無 systematic bias。57 windows som_ratio mean 0.462 跨 0-1 全範圍 = 真實 sub-clone signal。chr19:17M 對稱 0.500 = SP1 附近 paired 認雙 sub-clone (vs TO 113:0 失衡)。",
    tier3="longphase-s codebase / SomaticHaplotagProcess.cpp:533 / paired 軸對齊")

add(id="16_v5_caveat", num="16", section="S6 5/9-5/10 V5 缺陷 + V6 修補 ★", rg2="5 (核心 caveat + 修補 forced)", ngrep="20+",
    title="V5 Layer 1.5 設計缺陷 → V6 binary patch (5/10) 三向 head-to-head 修補",
    en="V5 Layer 1.5 caveat → V6 patch (5/10) three-way head-to-head fix",
    timing="180 sec / 中 ~480 字 ★ 整份最關鍵+最長",
    title_critical=True,
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#DC2626;margin:0;">問題 (5/9 發現): paired chr19 germline-absent 5,789 events</p>
    <img class="fig-thumb" src="../figures/G4_germline_absent_three_versions.png" alt="G4" style="max-height:200px;">
    <p style="font-size:11px;font-weight:700;color:#16A34A;margin:4px 0 0;">解答 (5/10 V6 patch): HaplotagProcess.cpp:537-548 移除 Layer 1.5 → germline-absent 退回 hp=33</p>
    <table class="metric-table" style="font-size:10px;">
      <caption style="text-align:left;font-weight:600;color:#1E3A8A;">全基因組三向 head-to-head (2,464,863 reads)</caption>
      <thead><tr><th>指標</th><th>V3F</th><th>V5</th><th>V6</th><th>V6 vs V3F</th></tr></thead>
      <tbody>
        <tr><td>hp=33 ambiguous reads</td><td class="num">132,060</td><td class="num">13,250</td><td class="num"><strong>138,317</strong></td><td class="num green">+4.7% ✅</td></tr>
        <tr><td>hp=1-1:hp=2-1 ratio</td><td class="num">1.138</td><td class="num">2.003</td><td class="num">1.838</td><td class="num">部分改善 ⚠</td></tr>
        <tr><td>marker coverage (NG≥3)</td><td class="num">21,997</td><td class="num">18,382</td><td class="num"><strong>23,980</strong></td><td class="num green">+9.0% ✅</td></tr>
        <tr><td>marker rate (off)</td><td class="num">0.9175</td><td class="num">0.8937</td><td class="num">0.9093</td><td class="num">介於 V3F/V5 ⚠</td></tr>
        <tr class="row-green"><td>caller F1 vs SEQC2</td><td class="num">0.7166</td><td class="num">0.7166</td><td class="num">0.7166</td><td>三版相同 ✅</td></tr>
      </tbody>
    </table>
    <div class="green-box" style="font-size:12px;font-weight:600;">
      ★ V6 = V5 設計目標保留 + V3F 保守 hp=33 策略 + marker engineering 改善 → 4/5 Phase C 通過
      <span style="display:block;font-size:10px;font-style:italic;margin-top:4px;color:#374151;">V6 reverts V5 Layer 1.5; Phase B (chr19) 3/5 + Phase C (genome) 4/5 PASS</span>
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">V6:</span> V5 + Layer 1.5 移除 (5/10 binary patch on top of 938f0df)</div>
      <div class="gloss-item">ⓘ V5→V6 transfer 82% from hp=1-1 / 17% from hp=2-1 → hp=33 (守恆)</div>
    </div>""",
    speaker="整份報告最重要的一張 slide — V5 Layer 1.5 缺陷 + V6 修補敘事完整鏈。問題: 5/9 paired audit Step D 揭露 V5 Layer 1.5 在 germline-absent 區與 baseline 4.19:1 完全相同 (5,789 events, V3F 標 hp=33 保守反而穩健)。解答: 5/10 V6 binary patch — 在 HaplotagProcess.cpp:537-548 移除 Layer 1.5 邏輯 (`else if (somaticHP1 > 0 || somaticHP2 > 0)` 整段去掉), germline-absent 退回 hp=33 (V3F 行為)。全基因組三向 head-to-head: hp=33 V6=138,317 > V3F=132,060 (+4.7%, 完全恢復 V3F 保守策略並略多), V5 只 13,250 (-89.9%) 為 priority bug feature 化副作用。Marker coverage V6 23,980 > V3F 21,997 (+9.0%) > V5 18,382 (V6 比 V5 多 +30.5%) — V6 抓到最多 NG≥3 regions。Marker rate V6 0.9093 介於 V3F 0.9175 與 V5 0.8937 之間。hp=1-1:hp=2-1 ratio V6 1.838 改善但未到 V3F 1.138 (因 V6 重用 V5 phased VCF, germline-existent 區 V5 ploidy/threshold 差異 V6 沒改)。Caller F1 三版完全相同 0.7166 (V6 不改 phasing layer)。V5→V6 transfer 82:17 比例 from hp=1-1/hp=2-1 → hp=33 守恆。Phase B chr19 3/5 + Phase C genome 4/5 通過; V6 = V5 設計目標保留 + V3F 保守策略 + marker improvement = production candidate。",
    tier3="V6 patch diff / 82:17 V5→V6 transfer / Phase B chr19 5 驗收項 / Phase C genome 5 驗收項")

# ─── S7 Errata + Follow-up ────────────────────────────────────────────────
add(id="17_errata", num="17", section="S7 結論", rg2="—", ngrep="commit chain",
    title="PI 報告 4-29 6 條 errata 已 patch; V6 patch 完成主結論升級",
    en="6 errata patched; V6 patch completes the main conclusion arc",
    timing="100 sec / 中 ~370 字",
    canvas_html="""
    <table class="metric-table" style="font-size:11px;">
      <thead><tr><th></th><th>段落</th><th>原 → 新</th></tr></thead>
      <tbody>
        <tr><td><strong>E1</strong></td><td>§3.3.3 chr19 SP1/2/3</td><td>主要 hotspot → 可重現案例 (chr19 占 2.16% rank 19)</td></tr>
        <tr><td><strong>E2</strong></td><td>§5.2 V5 commit 狀態</td><td>未 commit → ✅ 已 commit (d0bcd8c + 938f0df)</td></tr>
        <tr><td><strong>E3</strong></td><td>§5.2 priority bug 證據</td><td>commit msg + 3 IGV → + 34,855 read-level 鐵證 100%</td></tr>
        <tr class="row-yellow"><td><strong>★ E4</strong></td><td>§6.4/§6.5 V5 數值歸因</td><td>V5 整體 → V5 = Pass 1 only; 主要 V3F + Layer 1.5; Pass 2 二次效益尚未量化</td></tr>
        <tr class="row-yellow"><td><strong>★ E5</strong></td><td>§5.2 Layer 1.5 設計 (5/10 加)</td><td>fallback 隱含修補 → germline-absent 區與 baseline 4.19:1 同, priority bug feature 化</td></tr>
        <tr class="row-green"><td><strong>★ E6 NEW</strong></td><td>§9.3 follow-up (5/11 加)</td><td>F-paired-D3 ⏸ 待跑 → ✅ V6 binary patch DONE (5/10) + Phase D 跨 4 樣本驗證 4/5 通過 (5/11)</td></tr>
      </tbody>
    </table>
    <div class="conclusion-arrow green">→ E1-E3 表述精確化; E4-E5 核心 errata; ★ E6 = V6 修補敘事閉環 (主結論升級為 V6 production candidate)</div>
    <div class="footer-glossary"><div class="gloss-item">ⓘ commit chain: f17754f → 2553e96 → 71d21bd → V6 patch (待 commit)</div></div>""",
    speaker="PI 報告 4-29 6 條 errata (5/11 新加 E6): E1 chr19 hotspot 降級; E2 V5 已 commit 4-30; E3 證據鏈升級; ★ E4 V5 = Pass 1 only; ★ E5 V5 Layer 1.5 germline-absent 缺陷; ★ E6 (5/11 新加) = F-paired-D3 follow-up 已實作為 V6 binary patch (5/10) + Phase D 跨 4 樣本驗證 (5/11) 4/5 通過 — V6 = production candidate。E1-E3 為表述精確化, E4-E5 為核心 caveat, E6 為敘事閉環 (從 V5 caveat 升級為 V6 修補)。",
    tier3="errata commit chain / V6 patch 待 commit / Phase D 5 樣本 evaluation matrix")

add(id="18_followup", num="18", section="S7 結論", rg2="1", ngrep="—",
    title="整體成熟度升級 + Phase D 跨 4 樣本驗證 — V6 = production candidate",
    en="Maturity upgrade + Phase D 4-sample validation; V6 production candidate",
    timing="100 sec / 中 ~360 字",
    canvas_html="""
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:0;">整體成熟度升級: 13 ✅ / 0 ⚠️ / 2 ⏸ (15 維度)</p>
    <div class="grid-3col" style="font-size:9.5px;gap:5px;">
      <div class="green-box" style="padding:3px 6px;">✅ 機制因果</div>
      <div class="green-box" style="padding:3px 6px;">✅ 修補設計</div>
      <div class="green-box" style="padding:3px 6px;">✅ chr19 SP 對齊</div>
      <div class="green-box" style="padding:3px 6px;">✅ 全基因組 34,855</div>
      <div class="green-box" style="padding:3px 6px;">✅ 20 指標 0 regression</div>
      <div class="green-box" style="padding:3px 6px;">✅ Caller F1 三版相同</div>
      <div class="green-box" style="padding:3px 6px;">✅ purity 0.6 對照</div>
      <div class="green-box" style="padding:3px 6px;">✅ Pass 2 +3.51%</div>
      <div class="green-box" style="padding:3px 6px;">✅ Paired Step A+C+D</div>
      <div class="green-box" style="padding:3px 6px;"><strong>★ V6 patch (5/10)</strong></div>
      <div class="green-box" style="padding:3px 6px;"><strong>★ marker +9.0% vs V3F</strong></div>
      <div class="green-box" style="padding:3px 6px;"><strong>★ hp=33 恢復 +4.7%</strong></div>
      <div class="green-box" style="padding:3px 6px;"><strong>★ Phase D 4 樣本 4/5</strong></div>
      <div style="background:#F3F4F6;padding:3px 6px;border:1px solid #D1D5DB;border-radius:4px;">⏸ COLO829 (truth 權限)</div>
      <div style="background:#F3F4F6;padding:3px 6px;border:1px solid #D1D5DB;border-radius:4px;">⏸ T1.3 4-cell ablation</div>
    </div>
    <p style="font-size:11px;font-weight:700;color:#1E3A8A;margin:8px 0 0;">Phase D 跨 4 樣本 V6 evaluation (HCC1395 ref + 4 cancer cell lines):</p>
    <table class="metric-table" style="font-size:10px;">
      <thead><tr><th>Sample</th><th>hp ratio</th><th>marker rate</th><th>NG_on=2</th><th>caller F1</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>H1437</td><td class="num">1.243 ✅</td><td class="num">0.992 ✅</td><td class="num">0.991 ✅</td><td>n/a</td></tr>
        <tr class="row-green"><td>H2009</td><td class="num">0.901 ✅</td><td class="num">0.993 ✅</td><td class="num">0.992 ✅</td><td>n/a</td></tr>
        <tr class="row-green"><td>HCC1954</td><td class="num">0.958 ✅</td><td class="num">0.954 ✅</td><td class="num">0.967 ✅</td><td>n/a</td></tr>
        <tr class="row-yellow"><td>HCC1937</td><td class="num">0.611 ✅</td><td class="num">0.817 ⚠</td><td class="num">0.904 ✅</td><td>BRCA1 edge</td></tr>
        <tr style="background:#F3F4F6;"><td>COLO829</td><td colspan="4">⏸ truth set 0600 權限阻塞</td></tr>
      </tbody>
    </table>
    <div class="conclusion-arrow green">→ V6 = production candidate; ratio 4/4 中性 / marker 3/4 (HCC1937 edge); 剩 COLO829 + T1.3</div>
    <div class="footer-glossary">
      <div class="gloss-item">ⓘ HCC1937 0.817 = BRCA1 mutant FP/TP=0.194 樣本特性</div>
      <div class="gloss-item">📖 <span class="term">V6:</span> V5 + Layer 1.5 移除 (5/10 patch)</div>
    </div>""",
    speaker="整體成熟度從 12 維度升級至 15 維度: 13 ✅ + 0 ⚠️ + 2 ⏸。新加 4 維度全綠: ★ V6 patch (5/10) / ★ marker +9.0% vs V3F / ★ hp=33 恢復 +4.7% / ★ Phase D 4 樣本驗證 4/5 通過; V5 Layer 1.5 從 ⚠ 升級為 ✅ (V6 已修補)。Phase D 跨 4 樣本 (H1437/H2009/HCC1954/HCC1937, COLO829 truth set 0600 權限阻塞): hp=1-1:hp=2-1 ratio 4/4 全部接近中性 0.61-1.24 (vs V5 baseline 1.86, 原 baseline 17.3:1); marker rate 3/4 通過 ≥0.85 (HCC1937 0.817 為 BRCA1 mutant FP/TP=0.194 樣本特性 edge case, 與 HPFineNGroups canonical filter AF<0.4 一致); NG_on=2 rate 4/4 通過 ≥0.85; caller F1 V6 重用 V5 phased VCF 故不變。V6 = production candidate; 剩 2 follow-up: COLO829 待 truth set 權限解除, T1.3 4-cell ablation 量化 Pass 2 second round 獨立貢獻 — 兩者均非阻塞。",
    tier3="HCC1937 BRCA1 mutant detail / COLO829 truth set 權限 / T1.3 4-cell ablation 設計")

# ─── Q&A Backup ───────────────────────────────────────────────────────────
add(id="b1_pass2", num="B1", section="Q&A Backup", rg2="4 (backup 例外)", ngrep="—",
    title="Pass 2 = 只重跑 2-point edgeConnectResult; 高 purity 才觸發",
    en="Pass 2 = re-run 2-point only; high purity gate",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <img class="fig-thumb" src="../figures/G2_pass_two_flow.png" alt="G2" style="max-height:380px;">
    <p class="fig-caption">Pass 1 always runs 2-point + 3-point; Pass 2 only re-runs 2-point</p>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">somaticCalling vs edgeConnectResult:</span> 兩 phasing graph 演算法</div>
      <div class="gloss-item">📖 <span class="term">purity:</span> 樣本中腫瘤細胞占比</div>
    </div>""",
    speaker="Q: Pass 2 為什麼只跑 2-point? somaticCalling 用 3-point patternMining (與 purity 無關); edgeConnectResult 用 2-point 永遠跑。低 purity ≤0.9 Pass 2 跳過; 高 purity >0.9 Pass 2 只重跑 2-point (Pass 1 已產出穩定 origin 分類)。Pass 2 incremental: phased var -2.90 pp / blocks -9.79% / N50 +3.51%。常見誤解: 低 purity 用 3-point 倒過來; 高 purity 才多做事。",
    tier3="patternMining first/second/third path / Pass 2 不重跑 somaticCalling")

add(id="b2_purity06", num="B2", section="Q&A Backup", rg2="1", ngrep="20+",
    title="purity 0.6 樣本 baseline vs V5 完整對照 — 0 critical regression",
    en="purity 0.6 sample fully verified",
    timing="90 sec / 中 ~320 字",
    canvas_html="""
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;">6 Caller F1 (vs SEQC2 v1.2.1) — 三版完全相同:</p>
    <table class="metric-table" style="font-size:10px;">
      <thead><tr><th>指標</th><th>baseline 0.6</th><th>V5 0.6</th><th>Δ</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>TP / FP / FN</td><td>24,190 / 13,487 / 15,257</td><td>同</td><td>0 ✅</td></tr>
        <tr class="row-green"><td>F1 (Precision/Recall)</td><td>0.6273 (0.6420 / 0.6132)</td><td>同</td><td>0 ✅</td></tr>
      </tbody>
    </table>
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin-top:8px;">9 結構指標 — 4 改善 + 1 微差 + 1 持平:</p>
    <table class="metric-table" style="font-size:10px;">
      <thead><tr><th>指標</th><th>baseline</th><th>V5</th><th>Δ</th><th>eval</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>phased%</td><td>61.82</td><td>65.83</td><td>+4.01 pp</td><td>✅</td></tr>
        <tr class="row-green"><td>n_blocks</td><td>9,748</td><td>11,514</td><td>+18.1%</td><td>✅</td></tr>
        <tr class="row-yellow"><td>N50 (bp)</td><td>798,903</td><td>683,296</td><td>-14.5%</td><td>微差 (≥600K)</td></tr>
        <tr class="row-green"><td>HP:i:33</td><td>0</td><td>20</td><td>+20</td><td>✅ conservative</td></tr>
        <tr class="row-green"><td>purity 計算</td><td>0.607</td><td>0.634</td><td>+0.027</td><td>✅ closer to 0.6</td></tr>
      </tbody>
    </table>
    <div class="footer-glossary"><div class="gloss-item">📖 <span class="term">purity:</span> 樣本中腫瘤細胞占比</div></div>""",
    speaker="Q: 低純度樣本是否變差? purity 0.6 樣本 baseline vs V5: 6 caller F1 全相同 (V5 不改 caller); 9 結構: phased% +4.01 pp / n_blocks +18.1% / HP:i:33 +20 conservative / AMB% +3.12 pp / purity 0.607→0.634 / Pass 2 兩者都不觸發 / LOH 完全相同; 唯一 N50 -14.5% 但仍 ≥600K。0 critical regression; ploidy bug 在低 purity 自我治癒。",
    tier3="V3F vs V5 0.6 對比 / N50 為何接受")

add(id="b3_cross_sample", num="B3", section="Q&A Backup", rg2="1", ngrep="—",
    title="Phase D 跨 4 樣本 V6 詳細 evaluation matrix + COLO829 推遲",
    en="Phase D V6 4-sample evaluation + COLO829 deferred",
    timing="90 sec / 中 ~320 字",
    canvas_html="""
    <p style="font-size:11px;font-weight:700;color:#1E3A8A;margin:0;">5 樣本 V6 evaluation matrix (HCC1395 ref + 4 cancer cell lines):</p>
    <table class="metric-table" style="font-size:10px;">
      <thead><tr><th>Sample</th><th>ratio</th><th>hp=33</th><th>marker rate</th><th>NG_on=2</th><th>FP/TP</th></tr></thead>
      <tbody>
        <tr class="row-green"><td>HCC1395 (ref)</td><td class="num">1.838</td><td class="num">138,317</td><td class="num">0.909</td><td class="num">0.829</td><td class="num">—</td></tr>
        <tr class="row-green"><td>H1437</td><td class="num">1.243</td><td class="num">39,050</td><td class="num"><strong>0.992</strong></td><td class="num">0.991</td><td class="num">0.011</td></tr>
        <tr class="row-green"><td>H2009</td><td class="num">0.901</td><td class="num">684,035</td><td class="num"><strong>0.993</strong></td><td class="num">0.992</td><td class="num">0.010</td></tr>
        <tr class="row-green"><td>HCC1954</td><td class="num">0.958</td><td class="num">4,859</td><td class="num">0.954</td><td class="num">0.967</td><td class="num">0.035</td></tr>
        <tr class="row-yellow"><td>HCC1937 ⚠</td><td class="num">0.611</td><td class="num">5,017</td><td class="num">0.817</td><td class="num">0.904</td><td class="num"><strong>0.194</strong></td></tr>
        <tr style="background:#F3F4F6;"><td>COLO829</td><td colspan="5">⏸ truth set 0600 權限阻塞 (HKU/NYGC); BAM ready 待 chmod 660 或替代 PASS VCF</td></tr>
      </tbody>
    </table>
    <div class="grid-2col" style="margin:6px 0;">
      <div class="green-box" style="font-size:11px;">
        ✅ <strong>4 樣本 ratio 全部接近中性 0.611-1.243</strong> (平均 0.928, vs V5 1.86, 原 baseline 17.3:1) → V6 priority bug 修補一致成功
      </div>
      <div class="caveat-box" style="font-size:11px;">
        ⚠ <strong>HCC1937 marker rate 0.817</strong>: BRCA1 mutant FP/TP=0.194, CNV-driven germline het AF 漂移; 須加 AF&lt;0.4 filter (HPFineNGroups canonical)
      </div>
    </div>
    <p style="font-size:11px;color:#6B7280;font-style:italic;text-align:center;">Phase D pipeline: V5 binary phasing → V6 binary haplotag (hybrid); 每樣本 wall ~10 hr</p>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">BRCA1 mutant:</span> 高 ploidy CNV-driven FP edge case</div>
      <div class="gloss-item">ⓘ COLO829: 待 truth set chmod 660</div>
    </div>""",
    speaker="Q: 跨樣本 V6 一致性? Phase D 5 樣本 evaluation (HCC1395 ref + H1437/H2009/HCC1954/HCC1937, COLO829 truth set 權限阻塞)。hp=1-1:hp=2-1 ratio 4 樣本全部接近中性 0.611-1.243 (平均 0.928, vs V5 1.86, vs 原 baseline 17.3:1) — V6 priority bug 修補一致成功。hp=33 ambiguous 範圍 4,859-684,035 (V6 patch 正確標 hp=33)。Marker rate 3/4 通過 ≥0.85: H1437 0.992, H2009 0.993, HCC1954 0.954, HCC1937 0.817 邊緣 fail。HCC1937 為 BRCA1 mutant 高 ploidy 細胞株, FP/TP=0.194 (vs 其他 0.01-0.035), CNV-driven germline het AF 漂移混入 FP — 須加 AF<0.4 filter (與 HPFineNGroups canonical 一致), 非否定 V6 patch。COLO829 待用戶 chmod 660 HKU/NYGC truth set 或提供替代 PASS VCF, BAM ready。Phase D pipeline: V5 phasing → V6 haplotag (hybrid), 每樣本 wall ~10 hr。",
    tier3="HCC1937 BRCA1 mutant FP/TP=0.194 / COLO829 chmod 660 / Phase D 4-stage pipeline")


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

  <article class="slide-canvas{' cover-slide' if s.get('is_cover') else ''}" aria-label="Main slide canvas">
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
