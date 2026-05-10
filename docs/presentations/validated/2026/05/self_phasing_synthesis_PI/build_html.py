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

# ─── 22 slide specs ────────────────────────────────────────────────────────
# Each: id, num, section, title, en, timing, rg2, ngrep, canvas_html, speaker, tier3
# canvas_html = raw HTML body for the 16:9 .slide-canvas container
SLIDES = []

def add(**kw):
    kw.setdefault("title_critical", False)
    kw.setdefault("rg2", "—")
    kw.setdefault("ngrep", "—")
    kw.setdefault("tier3", "")
    SLIDES.append(kw)

# ─── S0 Cover + TL;DR ─────────────────────────────────────────────────────
add(id="01_cover", num="01", section="S0 Cover",
    title="Self-Phasing 整合觀察與 V5 Layer 1.5 設計缺陷",
    en="Self-Phasing Integration Synthesis — V5 Layer 1.5 design caveat",
    timing="30 sec / 中 ~120 字",
    canvas_html="""
    <div style="display:flex;flex-direction:column;justify-content:center;align-items:center;height:75%;text-align:center;">
      <h1 style="font-size:42px;color:#1E3A8A;margin:0;">Self-Phasing 整合觀察</h1>
      <p style="font-size:22px;color:#6B7280;font-style:italic;margin:8px 0 18px;">Self-Phasing Integration Synthesis</p>
      <p style="font-size:18px;color:#1E3A8A;">── V5 Layer 1.5 設計缺陷揭露 ──</p>
      <p style="font-size:14px;color:#6B7280;margin-top:8px;">longphase-to-mod 5 commits 修補成熟 + 5/9 paired cross-ref 新發現</p>
      <p style="font-size:13px;color:#6B7280;margin-top:24px;">2026-05-10  ·  PI / lab meeting</p>
      <p style="font-size:11px;color:#9CA3AF;">Source: 5/8 整合報告 + 5/9 errata + paired Step D</p>
    </div>""",
    speaker="今天 20 分鐘的報告主題是 self-phasing 整合觀察與 V5 Layer 1.5 設計缺陷揭露。整合 longphase-to-mod 5 commits 修補成熟度，以及 5/9 paired cross-ref 新發現。目的：協助 PI 決策 V5 是否作 production tag baseline，以及是否啟動 F-paired-D3 follow-up cycle。",
    tier3="5/8 主報告 1,211 行 + 7 figures commit hash")

add(id="02_tldr", num="02", section="S0 Cover", rg2="1", ngrep="6",
    title="TL;DR — 修補主線確立，但 V5 Layer 1.5 germline-absent 區仍待補強",
    en="Main fix line established; V5 Layer 1.5 still gaps in germline-absent regions",
    timing="90 sec / 中 ~360 字",
    canvas_html="""
    <div style="display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:8px 0;">
      <div>
        <h4 style="font-size:12px;color:#6B7280;margin:0 0 6px;">修補主線 (V3F + V5)</h4>
        <div class="grid-2col" style="gap:8px;margin:0;">
          <div class="stat-box"><div class="number">17.3:1</div><div class="label">HP1 偏 baseline</div></div>
          <div class="stat-box"><div class="number">34,855</div><div class="label">read-level victims</div></div>
        </div>
      </div>
      <div>
        <h4 style="font-size:12px;color:#6B7280;margin:0 0 6px;">關鍵驗證</h4>
        <div class="grid-2col" style="gap:8px;margin:0;">
          <div class="stat-box"><div class="number green">+13.3 pp</div><div class="label">paired GT (V5)</div></div>
          <div class="stat-box"><div class="number green">100%</div><div class="label">V3F+V5 修正率</div></div>
        </div>
      </div>
    </div>
    <div style="background:#F0F9FF;border:1px solid #BAE6FD;padding:8px 14px;border-radius:6px;text-align:center;font-weight:600;color:#0C4A6E;font-size:13px;">
      整體影響: 20/0 指標 no regression (caller F1 三版完全相同)
    </div>
    <div class="caveat-box">
      <span class="label">⚠ Caveat — 5/9 新發現:</span>
      V5 Layer 1.5 germline-absent 區與 baseline <strong>4.19:1 偏 HP1 完全相同</strong>; V3F 標 hp=33 反而更穩健 (~5% germline-absent 區，占比小不阻擋整體)
      <span class="en-note">V5 Layer 1.5 retains priority bug bias in germline-absent regions (~5% chr19 events)</span>
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">sub-clone:</span> 腫瘤內基因型相同細胞群</div>
    </div>""",
    speaker="今天 thesis 兩個焦點。正面: 修補主線確立。baseline LongPhase-TO 全基因組 HP1:HP2 = 17.3:1 systematic bias; V3F + V5 兩層修補在 read-level 對 chr19 752 + 全基因組 34,855 victims 修正率 100%; paired GT +13.3 pp; 20 指標 0 regression; caller F1 三版完全相同。反面: V5 Layer 1.5 germline-absent 區仍未修對 — 與 baseline 4.19:1 完全相同; V3F hp=33 更穩健。整體 12 維度 10 ✅ / 1 ⚠️ / 2 待跑; V5 仍可作 production baseline。",
    tier3="V3F = 41ff147 / V5 = d0bcd8c + 938f0df / paired Step D = 766ec5f")

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
      <div style="text-align:center;align-self:center;">
        <p style="font-size:24px;color:#FFA500;font-weight:700;margin:0;">HP1 占比 94.6%</p>
        <p style="font-size:14px;color:#9A3412;margin:4px 0;">↓↓↓↓↓↓↓↓↓</p>
        <p style="font-size:14px;color:#6B7280;">隨機預期 50%</p>
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
    title="SP1 chr19:17,565,944 — baseline 113:0 → V5 翻轉至 HP2 ✅ paired",
    en="SP1 — baseline 113:0 → V5 flips to HP2, aligned with paired",
    timing="60 sec / 中 ~250 字",
    canvas_html="""
    <table class="metric-table">
      <thead><tr><th>位點</th><th>baseline HP1:HP2</th><th>V5</th><th>paired GT</th><th>對齊?</th></tr></thead>
      <tbody><tr class="row-red"><td><strong>SP1 chr19:17,565,944</strong></td><td class="num">113 : 0</td><td>翻轉至 HP2 主導</td><td>HP2</td><td><strong>✅</strong></td></tr></tbody>
    </table>
    <img class="igv-thumb" src="../figures/igv/D_SP1_chr19_17565944.png" alt="IGV SP1">
    <p class="fig-caption">6-BAM 並列: baseline / V2b / V3F / V5 / paired_T / paired_N</p>
    <div class="arg-list">
      <strong>為何能排除是噪音 / caller / alignment?</strong><br>
      • baseline 與 paired 方向相反 (不是衰減而是翻轉)<br>
      • V5 修正後與 paired ground truth 重合<br>
      → read assignment 強制集中的鐵證
    </div>
    <div class="footer-glossary"><div class="gloss-item">📖 <span class="term">haplotype:</span> 父系/母系染色體</div></div>""",
    speaker="全基因組 17.3:1 是平均值; IGV 6-BAM 並列篩到 chr19 三個近 100% 失衡位點。SP1 chr19:17,565,944: baseline 113 reads 全 HP1, HP2=0; V5 翻 HP2 與 paired tumor 一致。排除噪音/caller/alignment: baseline 與 paired 完全反向 (翻轉而非衰減), V5 與 paired 重合 → read assignment 強制集中的鐵證。",
    tier3="6-BAM 並列順序 / V2b 中間階段意義")

add(id="04b_sp2_sp3", num="04b", section="S1 觀察起點", rg2="0", ngrep="6",
    title="SP2 + SP3 並列 — 同模式: baseline 全 HP1 / V5 翻 HP2 / 3/3 對齊 paired",
    en="SP2 + SP3 same pattern; 3/3 aligned with paired",
    timing="60 sec / 中 ~260 字",
    canvas_html="""
    <div class="grid-2col">
      <div>
        <table class="metric-table"><thead><tr><th>位點</th><th>baseline</th><th>V5</th></tr></thead>
        <tbody><tr class="row-red"><td><strong>SP2</strong><br>chr19:12,452,332</td><td class="num">109 : 1</td><td>HP2 主導</td></tr></tbody></table>
        <img class="igv-thumb" src="../figures/igv/D_SP2_chr19_12452332.png" alt="IGV SP2" style="max-height:200px;">
      </div>
      <div>
        <table class="metric-table"><thead><tr><th>位點</th><th>baseline</th><th>V5</th></tr></thead>
        <tbody><tr class="row-red"><td><strong>SP3</strong><br>chr19:12,467,180</td><td class="num">108 : 0</td><td>HP2 主導</td></tr></tbody></table>
        <img class="igv-thumb" src="../figures/igv/D_SP3_chr19_12467180.png" alt="IGV SP3" style="max-height:200px;">
      </div>
    </div>
    <div class="conclusion-arrow green">→ 三 SP 都在 chr19:12-17M 區段 → 對齊 slide 09 chr19 752 victims hotspot 區</div>
    <p style="font-size:12px;color:#6B7280;">引出四問: ① 為何全集中一邊? (S2 機制) ② read 層級? (S3 量化) ③ 三版各修? (S4 修補) ④ 是否都修對? (S5 驗證)</p>""",
    speaker="SP2 chr19:12,452,332 baseline 109:1; SP3 chr19:12,467,180 baseline 108:0; 與 SP1 同模式: baseline 全 HP1, V5 翻 HP2 對齊 paired 3/3。三 SP 都在 chr19:12-17M, 對齊 chr19 752 victims hotspot 區段 — read-level 個案與 IGV 屬同機制不同層級。引出四問機制/量化/修補/驗證。",
    tier3="paired_T 與 paired_N 對照細節")

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
    title="getVote 三版差異 — baseline 順序錯 → V3F 兩層 → V5 +Layer 1.5 fallback",
    en="getVote 3-version diff",
    timing="120 sec / 中 ~360 字",
    canvas_html="""
    <img class="fig-thumb" src="../figures/G3_getVote_three_layer.png" alt="G3" style="max-height:380px;">
    <div class="caveat-box">
      <span class="label">⚠ Layer 1.5 caveat:</span>
      V5 Layer 1.5 在 germline-absent 區會繼承 priority bug 偏移 → <strong>(slide 16 詳述)</strong>
    </div>""",
    speaker="三版 code: baseline 紅底 vector 順序錯 break early。V3F 綠底 explicit Layer 1 (germline only) + Layer 2 (somatic annotation), Layer 1 永不被 somatic overrule。V5 綠底加黃 highlight: 保留 Layer 1, 新增 Layer 1.5 fallback (germline 缺席用 somatic phased votes), Layer 2 同 V3F。Layer 1.5 是 caveat 預告 (slide 16 詳述)。",
    tier3="enum / V3F bonus 修")

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
      → V5 不改 caller; ΔF1 (0.93→0.6) = -0.0893 為 ClairS-TO 性質<br>
      <strong>→ 但 5/9 paired cross-ref 揭露另一面...</strong> 〔next: slide 15-16〕
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

add(id="16_v5_caveat", num="16", section="S6 5/9 新發現 ★", rg2="4 (核心 caveat forced)", ngrep="15+",
    title="V5 Layer 1.5 設計缺陷: germline-absent 區與 baseline 4.19:1 完全相同",
    en="V5 Layer 1.5 design caveat — identical to baseline in germline-absent",
    timing="150 sec / 中 ~420 字 ★ 最長 slide",
    title_critical=True,
    canvas_html="""
    <img class="fig-thumb" src="../figures/G4_germline_absent_three_versions.png" alt="G4" style="max-height:340px;">
    <div class="red-frame">
      <span class="label">★ 結論:</span>
      <strong>V5 Layer 1.5 = priority bug 的 feature 化非修補</strong> — 把 baseline buggy 行為改成 designed 行為，但<strong>該區域偏移本質沒變</strong>。V3F 標 hp=33 反而更穩健。
      <span style="display:block;font-size:10px;font-style:italic;margin-top:4px;">V5 Layer 1.5 makes priority bug a feature, not a fix.</span>
    </div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">sub-clone:</span> 基因型同細胞群</div>
      <div class="gloss-item">ⓘ Layer 1.5 fallback / phased votes</div>
      <div class="gloss-item">ⓘ scope: paired chr19 germline-absent (5,789 events ≈ 5%)</div>
    </div>""",
    speaker="整份報告最重要的一張 slide。對 paired chr19 × baseline/V3F/V5 vote dump JOIN, 篩 cnt_HP1+HP2=0 且 somatic>0 = 5,789 events。cross-tab: baseline hp=11=3,312 / hp=21=791 → 4.19:1 priority bug 次峰。V3F 全標 hp=33 保守 ✅。V5 hp=11=3,313 / hp=21=790 → 4.19:1 與 baseline 完全相同！機制詮釋: V5 Layer 1.5 用 sHP1 vs sHP2 票數決方向; sub-clone somatic 100% 共現 → graph 偏向同一 haplotype → 投票偏向同邊 → 繼承 priority bug 偏移。V5 Layer 1.5 = priority bug feature 化非修補; V3F hp=33 是更穩健選擇。V5 設計時未對 paired ground truth 做 germline-absent cross-ref, 5/9 paired audit Step D 才補上發現。",
    tier3="cross-binary axis alignment / phasing graph §3.2 / Layer 1.5 改回 V3F default / F-paired-D3 工作量")

# ─── S7 Errata + Follow-up ────────────────────────────────────────────────
add(id="17_errata", num="17", section="S7 結論", rg2="—", ngrep="commit chain",
    title="PI 報告 4-29 5 條 errata 已 patch; 主結論不撤回",
    en="5 errata patched; main conclusions retained",
    timing="90 sec / 中 ~330 字",
    canvas_html="""
    <table class="metric-table" style="font-size:11px;">
      <thead><tr><th></th><th>段落</th><th>原 → 新</th></tr></thead>
      <tbody>
        <tr><td><strong>E1</strong></td><td>§3.3.3 chr19 SP1/2/3</td><td>主要 hotspot → 可重現案例 (chr19 占 2.16% rank 19)</td></tr>
        <tr><td><strong>E2</strong></td><td>§5.2 V5 commit 狀態</td><td>未 commit → ✅ 已 commit (d0bcd8c + 938f0df)</td></tr>
        <tr><td><strong>E3</strong></td><td>§5.2 priority bug 證據</td><td>commit msg + 3 IGV → + 34,855 read-level 鐵證 100%</td></tr>
        <tr class="row-yellow"><td><strong>★ E4</strong></td><td>§6.4/§6.5 V5 數值歸因</td><td>V5 整體 → V5 = Pass 1 only; 主要 V3F + Layer 1.5; Pass 2 二次效益尚未量化</td></tr>
        <tr class="row-yellow"><td><strong>★ E5</strong></td><td>§5.2 Layer 1.5 設計 (5/10 加)</td><td>fallback 隱含修補 → germline-absent 區與 baseline 4.19:1 相同, priority bug feature 化非修補</td></tr>
      </tbody>
    </table>
    <div class="conclusion-arrow">→ E1-E3 表述精確化; E4 + E5 為核心 errata; 主結論不撤回</div>
    <div class="footer-glossary"><div class="gloss-item">ⓘ commit chain: f17754f → 2553e96 → 71d21bd</div></div>""",
    speaker="PI 報告 4-29 5 條 errata: E1 chr19 SP1/2/3 hotspot 降級為可重現案例; E2 V5 已 commit 4-30; E3 證據鏈升級含 34,855 read-level; ★ E4 V5 數值 = Pass 1 only, 主要 V3F + Layer 1.5; ★ E5 5/10 加 — V5 Layer 1.5 germline-absent 與 baseline 4.19:1 同。E1-E3 表述精確化, E4+E5 核心。",
    tier3="errata commit hash / patch 工作量")

add(id="18_followup", num="18", section="S7 結論", rg2="1", ngrep="—",
    title="整體成熟度 + 5 follow-up — V5 仍可作 production baseline",
    en="Maturity status + 5 follow-up cycles",
    timing="90 sec / 中 ~340 字",
    canvas_html="""
    <p style="font-size:12px;font-weight:700;color:#1E3A8A;margin:0;">整體成熟度: 10 ✅ / 1 ⚠️ / 2 ⏸ (12 維度)</p>
    <div class="grid-3col" style="font-size:10px;gap:6px;">
      <div class="green-box" style="padding:4px 8px;">✅ 機制因果</div>
      <div class="green-box" style="padding:4px 8px;">✅ 修補設計</div>
      <div class="green-box" style="padding:4px 8px;">✅ chr19 SP 對齊</div>
      <div class="green-box" style="padding:4px 8px;">✅ 全基因組擴展</div>
      <div class="green-box" style="padding:4px 8px;">✅ V5/V3F zero-sum</div>
      <div class="green-box" style="padding:4px 8px;">✅ 20 指標 0 regression</div>
      <div class="green-box" style="padding:4px 8px;">✅ Caller F1 三版相同</div>
      <div class="green-box" style="padding:4px 8px;">✅ purity 0.6 對照</div>
      <div class="green-box" style="padding:4px 8px;">✅ 三路徑算法</div>
      <div class="green-box" style="padding:4px 8px;">✅ Pass 2 +3.51%</div>
      <div class="green-box" style="padding:4px 8px;">✅ 版本對齊 938f0df</div>
      <div class="green-box" style="padding:4px 8px;">✅ Paired Step A+C</div>
      <div class="caveat-box" style="padding:4px 8px;margin:0;">⚠ V5 Layer 1.5 (E5)</div>
      <div style="background:#F3F4F6;padding:4px 8px;border:1px solid #D1D5DB;border-radius:4px;">⏸ Pass 2 (T1.3)</div>
      <div style="background:#F3F4F6;padding:4px 8px;border:1px solid #D1D5DB;border-radius:4px;">⏸ 跨樣本 (T3)</div>
    </div>
    <table class="metric-table" style="font-size:10px;margin-top:8px;">
      <thead><tr><th>ID</th><th>內容</th><th>預估</th></tr></thead>
      <tbody>
        <tr class="row-yellow"><td>★ F-paired-D3</td><td>Layer 1.5 改 V3F ISM 影響 (決定 V5 vs V3F default)</td><td>1-2 day</td></tr>
        <tr><td>F-paired-D1</td><td>germline-absent 全基因組擴展</td><td>0.5 day</td></tr>
        <tr><td>F-paired-D2</td><td>phase block 內 axis-aligned 分析</td><td>1 day</td></tr>
        <tr><td>T3</td><td>7 樣本跨樣本擴展</td><td>1-2 day</td></tr>
        <tr><td>T1.3</td><td>4-cell ablation</td><td>3 day</td></tr>
      </tbody>
    </table>
    <div class="conclusion-arrow">→ V5 仍可作 production baseline; F-paired-D3 量化後決定是否回 V3F default</div>
    <div class="footer-glossary"><div class="gloss-item">📖 <span class="term">LOH/cnLOH:</span> 雜合性丟失/拷貝中性</div></div>""",
    speaker="整體成熟度 12 維度: 10 ✅ + 1 ⚠️ V5 Layer 1.5 + 2 ⏸ 待跑。5 follow-up 排序按 ROI: ★ F-paired-D3 1-2 day 最重要 (V5 vs V3F default); F-paired-D1 0.5 day; F-paired-D2 1 day; T3 1-2 day; T1.3 3 day。V5 仍可作 production baseline; germline-absent ~5% 不阻擋整體; F-paired-D3 量化後決定是否回 V3F default — 待 PI 決策。",
    tier3="T1.3 4-cell / 7 樣本 / cnLOH 雙親同源")

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
    title="跨樣本擴展 (T3) + cnLOH 雙親同源待開放",
    en="Cross-sample expansion + cnLOH bi-parental",
    timing="60 sec / 中 ~250 字",
    canvas_html="""
    <div class="grid-2col">
      <div>
        <p style="font-size:12px;font-weight:700;color:#1E3A8A;">T3 跨樣本擴展:</p>
        <table class="metric-table" style="font-size:10px;">
          <thead><tr><th>樣本</th><th>狀態</th><th>預估</th></tr></thead>
          <tbody>
            <tr class="row-green"><td>HCC1395 5kHz</td><td>✅ 已驗證</td><td>—</td></tr>
            <tr><td>HCC1395_DORADO</td><td>⏸ 待跑</td><td>~3 hr</td></tr>
            <tr><td>HCC1937</td><td>⏸ 待跑</td><td>~3 hr</td></tr>
            <tr><td>HCC1954</td><td>⏸ 待跑</td><td>~3 hr</td></tr>
            <tr><td>H1437</td><td>⏸ 待跑</td><td>~3 hr</td></tr>
            <tr><td>H2009</td><td>⏸ 待跑</td><td>~3 hr</td></tr>
            <tr><td>COLO829</td><td>⏸ 待跑</td><td>~3 hr</td></tr>
          </tbody>
        </table>
      </div>
      <div>
        <p style="font-size:12px;font-weight:700;color:#1E3A8A;">cnLOH 雙親同源待開放:</p>
        <div class="caveat-box" style="font-size:11px;">
          <strong>情境:</strong> parent 1 vs parent 2 染色體在 LOH 區難 phase<br><br>
          <strong>影響:</strong> Layer 1.5 在 cnLOH 區的設計選擇待量化<br><br>
          <strong>關係:</strong> 與 V5 Layer 1.5 設計缺陷 (E5) 連動
        </div>
      </div>
    </div>
    <div class="conclusion-arrow">Follow-up 排序: F-paired-D3 (1-2 d, actionable) > T3 (1-2 d, generalizability) > T1.3 (3 d, ablation)</div>
    <div class="footer-glossary">
      <div class="gloss-item">📖 <span class="term">LOH:</span> 雜合性丟失</div>
      <div class="gloss-item">📖 <span class="term">cnLOH:</span> 拷貝中性 LOH</div>
    </div>""",
    speaker="Q: 跨樣本一致性? cnLOH 區? T3 跨樣本: HCC1395 5kHz 已驗證; 6 樣本待跑 (DORADO/1937/1954/H1437/H2009/COLO829) 每樣本 ~3 hr, 共 1-2 day。cnLOH 雙親同源: parent 1 vs parent 2 在 LOH 區同源難 phase; Layer 1.5 對 cnLOH 影響待量化, 與 E5 連動。Follow-up: F-paired-D3 actionable > T3 generalizability > T1.3 ablation。",
    tier3="7 樣本特性 / cnLOH 機制 / F-paired-D3 詳細")


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
  <link rel="stylesheet" href="shared/style.css">
</head>
<body>
<main class="slide-page">
  <div class="slide-meta {meta_class}">
    <div><span class="label">Slide:</span><span class="value">{s['num']} / 22</span></div>
    <div><span class="label">Section:</span><span class="value">{s['section']}</span></div>
    <div><span class="label">Tier 2:</span><span class="value">{s['timing']}</span></div>
    <div><span class="label">R-G2 術語:</span><span class="value">{s['rg2']}</span></div>
    <div><span class="label">N grep:</span><span class="value">{s['ngrep']}</span></div>
  </div>

  <div class="slide-canvas">
    <h2 class="slide-title {title_class}">{s['title']}</h2>
    <p class="en-subtitle">{s['en']}</p>
    {s['canvas_html']}
  </div>

  <div class="note-section">
    <details class="collapsible" open>
      <summary>🎤 Speaker Note ({s['timing']})</summary>
      <div class="speak-text">{s['speaker']}</div>
      {f'<div class="tier3"><span class="label">[ORAL-OPTIONAL]</span> {s["tier3"]}</div>' if s["tier3"] else ''}
    </details>
  </div>

  <div class="cross-link">
    📋 詳細 multi-agent review (T/C/L/B/N + PI 6 並行) 與修正紀錄請見：
    <a href="../03_slide_layout_script.md">InterSubMod/docs/presentations/.../03_slide_layout_script.md</a>
  </div>
</main>

<div class="nav-bottom">
  {prev_link}
  <span class="counter">Slide {s['num']} ({idx + 1}/{total})</span>
  {next_link}
</div>

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
  <link rel="stylesheet" href="shared/style.css">
</head>
<body>

<header class="topbar">
  <h1>📊 Self-Phasing PPT Preview — 22 slides multi-page</h1>
  <div class="meta">
    <code>InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/</code>
    &nbsp;|&nbsp; Build: 2026-05-11
  </div>
</header>

<nav class="slidetabs" id="tabs">
  {nav_html}
</nav>

<div class="iframe-wrap">
  <iframe id="slide-frame" src="slide_{SLIDES[0]['id']}.html"></iframe>
</div>

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
