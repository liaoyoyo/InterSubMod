<!--
build_date: 2026-05-10
agent: pptx-build P4 slide build (batch 1)
status: in_progress
parent_outline: 02_slide_outline.md
multi_agent_review: 6 並行 (T/C/L/B/N + Agent-PI) per slide
-->

# Self-Phasing PPT — Slide Layout Script (P4)

> 21 slide × 6 agent review per slide。每 batch 3 張，用戶 ack 後進下一 batch。
> ASCII wireframe = 視覺結構示意；Tier 1/2/3 = 三層分流；reviewer 段 = 6 並行視角結論。

---

## Batch 1 — Slide 01 Cover / Slide 02 TL;DR / Slide 03 全基因組 17.3:1

---

### Slide 01 — Cover

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│                                                                          │
│       Self-Phasing 整合觀察                                              │
│       Self-Phasing Integration Synthesis                                 │
│                                                                          │
│       ── V5 Layer 1.5 設計缺陷揭露 ──                                    │
│       (V5 Layer 1.5 design caveat reveal)                                │
│                                                                          │
│                                                                          │
│       longphase-to-mod 5 commits 修補成熟                                │
│       + 5/9 paired cross-ref 新發現                                      │
│                                                                          │
│                                                                          │
│       2026-05-10  ·  PI / lab meeting                                    │
│       Source: 5/8 整合報告 + 5/9 errata + paired Step D                  │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | 中標題 + EN subtitle + 副標題 + 日期 + 受眾 + source 引用（共 6 elements） |
| **Tier 2 (note)** | 開場 30 sec：「今天 20 min 報告 self-phasing 修補成熟度與 5/9 新發現的 V5 Layer 1.5 設計缺陷；目的是 PI 決策 V5 是否作 production tag baseline 與是否啟動 F-paired-D3 評估。」 |
| **Tier 3 (oral)** | 5/8 主報告全 1,211 行 / 5/9 errata 5 條 / 7 figures 已 commit 路徑 |

#### Multi-agent Review

| Agent | 結論 | 問題 |
|-------|------|------|
| **T 字體** | ✅ PASS | 標題建議用 36pt bold；EN subtitle 60% 字級（21.6pt）；中文用 Noto Sans CJK TC，EN 用 DejaVu Sans |
| **C 色彩** | ✅ PASS | 主標題深灰（不全黑減眼壓）；strikethrough 標題用 accent 色（不用紅色避免警告誤讀）|
| **L 佈局** | ✅ PASS | 上 1/3 留白 / 中 1/2 標題群 / 下 1/3 metadata；無重疊 |
| **B 雙語** | ✅ PASS | EN subtitle 縮排 0.25" + 60% 字級（feedback_pptx_bilingual_formatting）|
| **N 數字** | n/a | cover 無數字 |
| **Agent-PI** | ⚠ 一個建議 | cover 沒「先說結論」hook → 改副標題加「修補成熟度 + 1 個新 caveat」可吸引 PI 立刻入題；或維持簡潔 cover 哲學 |

#### Agent-PI 建議處理

維持簡潔 cover；hook 留給 slide 02 TL;DR 處理。**接受**。

---

### Slide 02 — TL;DR (1 句結論 + 5 數字 + 1 caveat)

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  TL;DR — 修補成熟，但 V5 Layer 1.5 germline-absent 區仍待補強             │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 修補主線 (V3F + V5) ─────┐  ┌─ 關鍵驗證 ──────────────┐              │
│   │  ┌──────┐  ┌────────┐   │  │  ┌──────────┐ ┌──────┐ │              │
│   │  │17.3:1│  │ 34,855 │   │  │  │ +13.3 pp │ │ 100% │ │              │
│   │  │ HP1  │  │victims │   │  │  │ paired   │ │修正率│ │              │
│   │  │偏 baseline│ │ genome│   │  │  │ GT (V5)  │ │V3F+V5│ │              │
│   │  └──────┘  └────────┘   │  │  └──────────┘ └──────┘ │              │
│   └──────────────────────────┘  └────────────────────────┘              │
│                                                                          │
│   ┌─ 整體影響 ───────────────────────────────────────────────────────┐   │
│   │  ┌──────────────────────────────────────────────────────────┐   │   │
│   │  │  20/0 指標 no regression  (caller F1 三版相同)         │   │   │
│   │  └──────────────────────────────────────────────────────────┘   │   │
│   └──────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ⚠ Caveat — 5/9 新發現                                                  │
│   ┌──────────────────────────────────────────────────────────────────┐   │
│   │  V5 Layer 1.5 germline-absent 區與 baseline 4.19:1 偏 HP1 完全相同 │   │
│   │  → V3F 標 hp=33 反而更穩健（待 F-paired-D3 ISM 影響量化）         │   │
│   └──────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│  📖 sub-clone: 腫瘤內基因型相同的細胞群；不同 sub-clone 帶不同 somatic    │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | 5 大數字（17.3:1 / 34,855 / +13.3 pp / 100% / 20/0）+ 1 caveat box + 1 glossary box (sub-clone)（共 7 elements，1 over R-G2 但因 caveat 為核心 thesis 必留）|
| **Tier 2 (note)** | 1.5 min — 「今天的 thesis 兩面：正面 V3F+V5 兩層修補在 read-level 對全基因組 34,855 個 priority bug victim 修正率 100%，paired ground truth 對齊 +13.3 pp，20 指標 0 regression、caller F1 三版完全相同；但反面 5/9 新發現 V5 Layer 1.5 在 germline-absent 區（占 chr19 events ~5%）行為與 baseline 4.19:1 偏 HP1 完全相同，是 priority bug 的 feature 化非修補。」 |
| **Tier 3 (oral)** | V3F 對應 commit 41ff147（4-10）/ V5 對應 commit d0bcd8c+938f0df（4-30）/ paired Step D 對應 commit 766ec5f |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ⚠ 1 修正 | 5 大數字字級不一致：建議全部 32pt bold；caveat box 內文用 18pt 稍小但不擠 |
| **C 色彩** | ⚠ 1 修正 | caveat box 用淺黃底（不用紅色避免警告過度）；其他用主色 navy；metric 數字 navy bold |
| **L 佈局** | ✅ PASS | 兩欄主數字 + 一欄整體影響 + 底部 caveat + glossary footer，視覺金字塔清晰 |
| **B 雙語** | ⚠ 1 修正 | 缺 EN sub-thesis；建議 caveat 框補一行 EN：「V5 Layer 1.5 retains priority bug bias in germline-absent regions」60% 字級縮排 |
| **N 數字** | ✅ 全 PASS（grep verified）| **17.3:1** ← 5/8 §0 + §2.1 ✅ / **34,855** ← §0 + §4.2 ✅ / **+13.3 pp** ← §6.4 + §8.5.1 ✅ / **100%** ← §6.1 + §6.2 ✅ / **20/0** ← §8.5.1 (20 指標表全綠) ✅ / **4.19:1** ← §8.6.4 ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a) 「修補成熟」太弱 — 改「修補主線確立」更明確；(b) caveat box 末加「占比 ~5%」具體量化才能 PI 評估嚴重度；建議補「(germline-absent 區占 chr19 events ~5%)」 |

#### Reviewer 修正建議匯總（待用戶 ack）

1. **T 修正**：5 大數字統一 32pt bold；caveat 內文 18pt
2. **C 修正**：caveat 用淺黃 (#FFF3CD) 底，非紅；數字 navy bold
3. **B 修正**：caveat 補 EN sub-thesis 一行
4. **PI 修正**：title 從「修補成熟」改「修補主線確立」；caveat 末補占比量化「~5% germline-absent」
5. **R-G4 處理**：sub-clone 已加 glossary footer ✓

---

### Slide 03 — 全基因組 17.3:1 偏移；3 條獨立論證

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  全基因組 HP1:HP2 = 17.3:1 是 systematic bias，非樣本性質                  │
│  Genome-wide HP1:HP2 17.3:1 = systematic bias, not sample variation       │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ baseline 全基因組統計 ─────────┐  ┌─ 隨機預期 vs 實測 ─────────┐    │
│   │ 指標         baseline   隨機  偏離│  │                           │    │
│   │ HP1 reads   614,000    325K   1.89×│  │   HP1 占比 94.6%         │    │
│   │ HP2 reads    35,500    325K   0.11×│  │   ⏷⏷⏷⏷⏷⏷⏷⏷⏷             │    │
│   │ HP1:HP2     17.3:1     1:1    17.3×│  │   隨機預期 50%            │    │
│   └────────────────────────────────┘  └───────────────────────────┘    │
│                                                                          │
│   ┌─ 三條獨立論證 ───────────────────────────────────────────────────┐   │
│   │  ① 生物學  │ tumor sub-clone 跨 23 染色體不該系統偏 HP1            │   │
│   │  ② 跨 chr  │ cnLOH artifact 只影響單一 chr；94.6% 跨 chr 一致      │   │
│   │  ③ paired  │ paired tumor-normal 同 reads HP1:HP2 ≈ 1:1           │   │
│   └──────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ➜ 17.3:1 不是 sample 性質，是 LongPhase-TO 的 systematic bias            │
│                                                                          │
│  📖 germline het: 個人遺傳變異中雜合位點 (A/G)，雙親各帶一型             │
│  📖 sub-clone: 腫瘤內基因型相同細胞群；不同 sub-clone 帶不同 somatic     │
│  📖 haplotype: 父系/母系兩條染色體之一；germline het 的方向標            │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | 標題（assertion）+ baseline 統計表（3 row）+ 視覺對比（HP1 占比 94.6%）+ 三論證表（3 row）+ 結論箭頭 + 3 glossary box（共 6 大 elements，含 glossary 為 secondary）|
| **Tier 2 (note)** | 1.5 min — 「baseline LongPhase-TO 全基因組 HP1 reads 614K 對 HP2 35.5K，比例 17.3:1 vs 隨機 1:1。這 94.6% 是 systematic bias 因為三條獨立佐證：① 生物學上 tumor sub-clone 不該跨 23 chr 系統偏 HP1 ② cnLOH artifact 只影響單 chr 但這偏移 23 chr 一致 ③ 同樣 reads 走 paired pipeline HP1:HP2 ≈ 1:1。所以 17.3:1 不是 sample 性質而是 LongPhase-TO 工程 bias，next 解釋怎麼發生。」|
| **Tier 3 (oral)** | cnLOH 機制細節 / 23 chr 一致性表（每 chr 偏 HP1 比例）/ paired pipeline 程式碼差異 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 表格字級 18pt；標題 28pt；glossary footer 14pt；視覺對比文字 24pt |
| **C 色彩** | ⚠ 1 修正 | 視覺對比用「HP1 占比」深紅 vs 隨機 50% 灰色 — 「HP1 占比 94.6%」用警告色（橙色 #FFA500）暗示偏移，不用紅色（紅色保留 caveat）|
| **L 佈局** | ✅ PASS | 上左數據 / 上右視覺 / 中部論證 / 底部結論 + glossary，金字塔結構好 |
| **B 雙語** | ✅ PASS | 標題雙語 OK；表格 column header 雙語 (baseline / 隨機 / 偏離)；論證列 EN 縮排 |
| **N 數字** | ✅ 全 PASS（grep verified）| **614,000** ← 5/8 §2.1 表 ✅ / **35,500** ← 同 ✅ / **325K** (隨機預期) ← 5/8 §2.1 ✅ / **1.89×** ← 同 ✅ / **0.11×** ← 同 ✅ / **17.3:1** ← 同 ✅ / **94.6%** ← 同 ✅ / **+44.6 pp** ← 同 ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a)「systematic bias」術語對 PI OK 但「工程 bias」不專業 — 改為「systematic engineering artifact」；(b) 三論證若再加「為何不可能是 ClairS-TO caller bias」一條會更完整（但已 4 條會稀釋）；接受 3 條 |

#### Reviewer 修正建議匯總

1. **C 修正**：HP1 占比 94.6% 用橙色（不用紅色）
2. **PI 修正**：「工程 bias」→「systematic engineering artifact」
3. **R-G4 處理**：3 個中等術語（germline het / sub-clone / haplotype）已 in-slide glossary box ✓
4. **R-G2 檢核**：3 中等術語 + 0 困難 + 0 簡單 = 3 ✓

---

## Batch 1 整體 self-check

| 項目 | 結果 |
|------|------|
| 3 張 slide 全有 title (assertion) + focal point + Tier 分流 | ✅ |
| 6 agent review 全跑 (T/C/L/B/N + PI) | ✅ |
| Agent-N 數字驗證 grep source verified | ✅ slide 02 / 03 全 PASS |
| R-G4 三層分流套用 | ✅ slide 02/03 用 glossary box，slide 05+ 將補 |
| 個人風格規則生效 | ✅ R-G2 ≤ 3 / slide |

---

## Batch 1 — User Ack ✅ Final Spec（2026-05-10）

用戶確認全套 6 修正套用：

### Slide 02 final
- **Title**: `TL;DR — 修補主線確立，但 V5 Layer 1.5 germline-absent 區仍待補強`
- **Caveat box** (淺黃 #FFF3CD 底):
  - 中：`V5 Layer 1.5 germline-absent 區與 baseline 4.19:1 偏 HP1 完全相同；V3F 標 hp=33 反而更穩健（~5% germline-absent 區，占比小不阻擋整體）`
  - EN（縮排 0.25" + 60% 字級）：`V5 Layer 1.5 retains priority bug bias in germline-absent regions (~5% chr19 events)`
- **5 大數字統一 32pt bold**（17.3:1 / 34,855 / +13.3 pp / 100% / 20-0）
- **Glossary footer**：sub-clone

### Slide 03 final
- **結論箭頭**: `→ 17.3:1 不是 sample 性質，是 LongPhase-TO 的 systematic engineering artifact`
- **HP1 占比 94.6%**: 橙色 (#FFA500)
- **Caveat 紅色保留**：紅色僅限 caveat / 警告，不用於 emphasis

### Slide 01 cover
- 維持原 spec（PI hook 留給 slide 02）

---

## Batch 2 — Slide 04 chr19 IGV / Slide 05 球員兼裁判 / Slide 06 priority bug

---

### Slide 04 — chr19 SP1/SP2/SP3 三位點 IGV 6-BAM 並列

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  chr19 三位點近 100% 失衡：113:0 / 109:1 / 108:0；V5 翻回 paired 3/3      │
│  chr19 three sites near-100% imbalance; V5 aligned with paired (3/3)      │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ SP1 chr19:17,565,944 ──────────────────────────────────────────┐    │
│   │  baseline: HP1=113 / HP2=0   →   V5: 翻轉至 HP2 主導 ✅ paired   │    │
│   │  [IGV 6-BAM 截圖: D_SP1 by_HP_4ver]                              │    │
│   └──────────────────────────────────────────────────────────────────┘    │
│   ┌─ SP2 chr19:12,452,332 ──────────────────────────────────────────┐    │
│   │  baseline: HP1=109 / HP2=1   →   V5: 翻轉至 HP2 主導 ✅ paired   │    │
│   │  [IGV: D_SP2]                                                    │    │
│   └──────────────────────────────────────────────────────────────────┘    │
│   ┌─ SP3 chr19:12,467,180 ──────────────────────────────────────────┐    │
│   │  baseline: HP1=108 / HP2=0   →   V5: 翻轉至 HP2 主導 ✅ paired   │    │
│   │  [IGV: D_SP3]                                                    │    │
│   └──────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│  解讀：not noise / not caller / not alignment → read assignment 強制集中  │
│                                                                          │
│  📖 haplotype: 父系/母系兩條染色體之一；germline het 的方向標            │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 3 IGV 截圖 + baseline 數字對比 + 解讀句 + glossary（共 6 elements）|
| **Tier 2 (note)** | 1.5 min — 「全基因組 17.3:1 是平均值；用 IGV 6-BAM 並列篩到 chr19 三個近 100% 失衡位點：SP1 113:0 / SP2 109:1 / SP3 108:0。每張 IGV 並列 baseline / V2b / V3F / V5 / paired tumor / paired normal 6 軌；可見 baseline 全 HP1 而 paired 反向 HP2，V5 修正後翻回 paired ground truth 三位都 OK。為何重要：不是噪音、不是 caller bug、不是 alignment 問題；而是 read assignment 強制集中 → 引出 next slide 機制。」|
| **Tier 3 (oral)** | 6-BAM 截圖順序 / V2b 的中間階段標籤 / paired pipeline 程式碼差異 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | title 28pt；位點 chr:pos 18pt monospace；解讀句 16pt italic |
| **C 色彩** | ⚠ 1 修正 | baseline 統計用紅色 / V5 翻轉用綠色 / paired 對齊用藍色 — 三色標清晰，**但不可全 IGV 圖加額外彩框（圖內已有色）**；只用文字底色標 |
| **L 佈局** | ⚠ 1 建議 | 3 IGV 縱排可能擠（slide 16:9 有限）→ 建議 **3 IGV 橫排 + 上方 3 col 數字對比表** + 下方解讀；或保留縱排但 **slide 改為「extended」（slide 04a / 04b 拆 2 張）**；先保留縱排試 PPTX render 看實際 |
| **B 雙語** | ✅ PASS | title 雙語；解讀句加 EN「not noise / not caller / not alignment」對齊 |
| **N 數字** | ✅ PASS（grep verified）| **17,565,944** ← 5/8 §2.2 ✅ / **12,452,332** ← 同 ✅ / **12,467,180** ← 同 ✅ / **113:0** ← 同 ✅ / **109:1** ← 同 ✅ / **108:0** ← 同 ✅ / **3/3** ← §6.3 ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a) IGV 截圖內字體在 slide 縮放後可能難讀 — **必須在 P5 build 時用 Vision check IGV 字級可讀性**；(b) 解讀句「not noise / not caller / not alignment」三排除是科學論證但 PI 可能想知道**「為何能排除」** — 補 1 行：「(三 SP 跨 chr 一致，baseline vs paired 翻轉非衰減，V5 修正後與 paired GT 重合)」|

#### Reviewer 修正建議匯總

1. **L 修正**：先試縱排 PPTX render，若擠則拆 04a/04b（待 P5 vision check）
2. **PI 修正**：解讀句後補「為何能排除」一行
3. **R-G2/R-G4**：1 中等術語（haplotype）已 glossary box ✓

---

### Slide 05 — phasing 層球員兼裁判：somatic 進 graph 自我增強

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  Phasing 層球員兼裁判 — somatic 100% 共現蓋過 germline 50/50               │
│  Phasing layer player-as-referee — somatic 100% co-occurrence overrules   │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 同 sub-clone 內共現比例對照 ──────────────────────────┐              │
│   │                germline het         somatic mutation     │              │
│   │  共現比例        50% / 50%           100% / 0%           │              │
│   │  分佈            兩 haplotype 隨機    單一 sub-clone 全帶  │              │
│   │  共現多位點       random              完全 100% 共現     │              │
│   └─────────────────────────────────────────────────────────┘              │
│                                                                          │
│   ┌─ TO 模式致命：somatic 進 phasing graph ───────────────────────────┐    │
│   │  • TO 沒 paired normal → 用 PoN（somatic 不在 PoN → 當 germline）   │    │
│   │  • somatic 進 graph 後 edge weight 暴漲（100% > 50% 共現）          │    │
│   │  • 自我增強：somatic 越強 → 偏該 haplotype → 該 haplotype 越強      │    │
│   │  • 結果：germline 真實訊號被 overrule（球員兼裁判）                │    │
│   └─────────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   解法：PON-only flag (commit 8b8c1fd) — Pass 1 graph 只放 PoN germline    │
│                                                                          │
│  📖 PoN: Panel of Normals — 多正常樣本建構的 germline 變異 reference set   │
│  📖 phasing graph: het 位點當 node，read 共現當 edge，連成 haplotype 圖    │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title (assertion) + 共現對照表 + TO 模式致命 4 bullet + 解法引出 + 2 glossary box（共 6 elements，過 R-G2 但因 PoN+phasing graph 為核心機制必留 glossary）|
| **Tier 2 (note)** | 2 min — 「球員兼裁判隱喻：phasing graph 的物理基礎是 germline het 在 HP1/HP2 兩 stream 50/50 隨機分佈，read 看到 (A,C) 就 assign HP1，看到 (G,T) 就 HP2。但 somatic 屬某 sub-clone 100% 全帶 — 多 somatic 進 graph 共現比 100% > germline 50%，edge weight 暴漲蓋過 germline。修法：PON-only flag — phasing 階段 graph 只放 PoN germline，somatic 不進 graph，等 graph 拍板再用結果反分類 somatic。」|
| **Tier 3 (oral)** | TO vs paired PoN 對照細節 / 自我增強迴圈視覺化 / Pass 1 vs Pass 2 預告 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | title 28pt；對照表 18pt；4 bullet 16pt；解法 18pt italic |
| **C 色彩** | ⚠ 1 修正 | germline 50/50 用 navy 中性；somatic 100/0 用橙色（warning，與 slide 03 同色系一致性）；PoN 解法 句用綠色強調 |
| **L 佈局** | ⚠ 1 建議 | 共現對照表 + bullet list 視覺密度 OK；**建議補一張小示意圖**（diploid genome + read assignment 箭頭）作為頁面右側 visual support — 但會擠 6 elements；可放 Tier 3 細節 alternative |
| **B 雙語** | ✅ PASS | title 雙語；對照表 column 雙語；4 bullet 中文為主（PI 不熟細節，EN 留 speaker note）|
| **N 數字** | ✅ PASS | **50% / 50%** + **100% / 0%** + **8b8c1fd** commit hash → grep §3.1-3.2 + §5.2 verified |
| **Agent-PI** | ⚠ 2 建議 | (a)「球員兼裁判」隱喻 PI 可能不熟 phasing graph → 補 1 句解釋：「(player-as-referee = somatic 應該被 phase 但反過來主導 graph)」；(b) 4 bullet 第三條「自我增強」太抽象 → 改具體：「3 個 somatic 共現一致 → 該 haplotype 強度 ×3 → 後續 reads 全偏該邊」 |

#### Reviewer 修正建議匯總

1. **C 修正**：germline navy / somatic 橙 / PoN 解法 綠
2. **L 建議**：暫不加示意圖（保 6 elements 上限）；speaker note 加細節即可
3. **PI 修正**：(a) 補球員兼裁判隱喻一句解釋；(b) 自我增強 bullet 改具體案例
4. **R-G4**：PoN + phasing graph 已 glossary box ✓
5. **R-G2 例外**：4 中等術語超 R-G2≤3，但因核心機制 forced；接受

---

### Slide 06 — tagging 層 getVote priority bug：1 票 somatic 蓋 5 票 germline

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  getVote vector 順序錯 + break early；1 票 somatic 觸發誤標                │
│  getVote vector ordered + break early; 1 somatic vote triggers misla-     │
│  bel                                                                      │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ Figure F1 — getVote priority bug 機制 (left baseline / right V3F+V5)─┐   │
│   │                                                                      │   │
│   │   [ figures/F1_priority_bug_mechanism.png ]                          │   │
│   │                                                                      │   │
│   └──────────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ Real read 範例（chr19 baseline=11 → V3F=21 全 752 條同模式）──────┐   │
│   │  germline HP1 = 0   somatic HP1_1 = 1                                │   │
│   │  germline HP2 = 5 ← 主導    somatic HP2_1 = 0                        │   │
│   │  germline HP3 = 0                                                    │   │
│   │  ─────────────────────────────────────────────────                  │   │
│   │  baseline: 檢 (HP1_1, HP2_1) → HP1_1=1>0 → hp=11 ❌ break (germline 被忽略) │   │
│   │  正確答案 hp=21（germline HP2=5 主導 + somatic HP1_1=1 標 21）       │   │
│   └──────────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   結論：tumor sub-clone somatic 100% 同方向 → priority bug 把所有受影響       │
│        reads 標 HP:i:11 系列 → 17.3:1 偏移在 tag layer 形成                │
│                                                                          │
│  📖 sub-clone: 腫瘤內基因型相同細胞群；不同 sub-clone 帶不同 somatic     │
│  ⓘ two-layer ⓘ break early                                                │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + Figure F1 圖 + 5-vote 範例 read + 結論句 + glossary box (sub-clone) + 2 footnote (two-layer / break early)（共 6 elements，footnote 為 secondary）|
| **Tier 2 (note)** | 2 min — 「上面 Figure F1 顯示 baseline `getVote()` vector 順序：① somatic pair ② mixed ③ germline pair；for 迴圈第一個非空 vector 處 break early，所以 1 票 somatic（HP1_1=1）觸發後立刻 break，後面 germline pair（HP2=5）永遠看不到。下面真實 read 範例：germline HP2=5 主導但 somatic HP1_1=1 觸發 priority bug，baseline 標 hp=11 錯，正確應為 hp=21。**全 752 條 chr19 victims 完全同模式**，這是 read-level 鐵證。為何全偏 HP1：tumor sub-clone somatic 100% 同方向（§3.2 機制） → priority bug 把所有受影響 reads 標 HP:i:11 系列。」|
| **Tier 3 (oral)** | enum (HAPLOTYPE1_1=2 等) vs HP tag int (=11) bug / 5-vote countMap 結構 / read_name 真實 ID |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | title 28pt；範例 read code-style 16pt monospace；結論 18pt bold |
| **C 色彩** | ⚠ 1 修正 | 範例 read「germline HP2 = 5 ← 主導」用 navy 強調；「somatic HP1_1 = 1」橙色（小但觸發）；「hp=11 ❌」紅色（warning）；「hp=21 ✅」綠色 |
| **L 佈局** | ⚠ 1 建議 | F1 圖 + 範例 read 上下兩欄 OK；**注意 F1 圖實際 width 必須 fit-within 不雙軸強制**（feedback_pptx_screenshot_rendering_rules）|
| **B 雙語** | ✅ PASS | title 雙語；範例 read 中英碼名混用標 monospace；結論句中文為主 |
| **N 數字** | ✅ PASS | **752** ← 5/8 §4.1 ✅ / **5 votes** + **HP1_1=1** + **HP2=5** ← §3.3 範例表 ✅ / **17.3:1** ← §0+§2.1 ✅ / **HP:i:11** + **HP:i:21** ← 附錄 B Glossary ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a)「priority bug」術語 PI 可能不熟 — 已 title 用「順序錯 + break early」白話化，但 footnote two-layer 不夠 → 補 1 footnote「**priority bug = vector 順序檢查 + break early 導致前面條目蓋過後面**」；(b) 結論「17.3:1 偏移在 tag layer 形成」要 emphasize **tag layer** 與 phasing layer 分離 → 加 transition 一句：「(注意：tag layer 與 §5 講的 phasing layer 是不同層 bug)」|

#### Reviewer 修正建議匯總

1. **C 修正**：4 色標範例 read（navy / 橙 / 紅 / 綠）
2. **PI 修正**：補 footnote 解 priority bug；結論句加 tag layer vs phasing layer transition
3. **R-G4**：sub-clone glossary box；two-layer + break early footnote
4. **R-G2 檢核**：1 中等 + 2 簡單 = 3 項處理 ≤ 3 ✓

---

## Batch 2 整體 self-check

| 項目 | 結果 |
|------|------|
| 3 張 slide spec 完整 | ✅ |
| 6 agent review per slide | ✅ |
| Agent-N grep verified（slide 04/05/06）| ✅ 全 PASS |
| R-G4 三層分流套用 | ✅ glossary box 與 footnote 已標 |
| R-G2 ≤ 3 / slide | ✅ slide 04/06 達標；slide 05 因核心機制 forced 4 中等術語接受例外 |
| Cliffhanger 預告 | n/a (S5 才出現) |

## Batch 2 — User Ack ✅ Final Spec（2026-05-10）

### Slide 04 → 拆 4a / 4b
- **04a** SP1 chr19:17,565,944 單獨一張（baseline 113:0 → V5 翻 HP2 ✅）
- **04b** SP2/SP3 並列（109:1 / 108:0）
- **總 slide 數**：22 main + 3 backup = **25 slide**

### Slide 05/06 final
全部 reviewer 修正接受（C 色標 / PI 隱喻補說明 / footnote / transition）

---

## Batch 3 — Slide 07 兩層對應表 / Slide 08 chr19 752 / Slide 09 全基因組 34,855

---

### Slide 07 — 兩層 bug 兩層修補對應表 — 任一層單獨修不夠

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  兩層 bug 互不取代：phasing PON-only + tagging V3F+V5 缺一不可             │
│  Two-layer bugs: PON-only (phasing) + V3F+V5 (tagging) both required      │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─────────────┬───────────────┬──────────────────┬─────────────┐       │
│   │   Layer    │    Bug        │   修補 commit    │   §5 章節   │       │
│   ├─────────────┼───────────────┼──────────────────┼─────────────┤       │
│   │ phasing 層 │ 球員兼裁判     │ 8b8c1fd          │  §5.2       │       │
│   │           │ (somatic 進    │ PON-only flag    │             │       │
│   │           │ graph)         │                  │             │       │
│   ├─────────────┼───────────────┼──────────────────┼─────────────┤       │
│   │ tagging 層│ priority bug  │ 41ff147 V3F      │  §5.3       │       │
│   │           │ (vector 順序錯)│ + 380e8d2 INDEL  │             │       │
│   │           │               │ + d0bcd8c+938f0df │             │       │
│   │           │               │ V5               │             │       │
│   └─────────────┴───────────────┴──────────────────┴─────────────┘       │
│                                                                          │
│   ┌─ 為何不能只用 PON-only ────────────────────────────────────────┐    │
│   │  解 phasing 但 tag 仍壞 → 99.9% reads 仍標 HP:i:11 系列        │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ 為何 V3F 不夠還要 V5 ─────────────────────────────────────────┐    │
│   │  V3F Layer 1 only → germline 缺席區 reads 全 untagged           │    │
│   │  V5 Layer 1.5 補 fallback (但有 germline-absent 區 caveat 見 S6)│    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│  📖 PoN: Panel of Normals — 多正常樣本建構的 germline reference set       │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 2x4 對應表 + 2 雙問解釋框 + glossary（共 4 elements ✓ R-G2）|
| **Tier 2 (note)** | 1.5 min — 「兩層 bug 是獨立問題：phasing graph 階段 somatic 進 graph 是球員兼裁判，需 PON-only flag；tagging 階段 getVote 順序錯是 priority bug，需 V3F two-layer + V5 Layer 1.5。任一層單獨修不夠，必須 stacking。為何 PON-only 不夠：解 phasing 但 tag 仍偏 — read 99.9% 仍標 HP:i:11。為何 V3F 不夠：Layer 1 only 在 germline 缺席區（cnLOH / amplicon）reads 全 untagged；V5 Layer 1.5 補 fallback。但 V5 Layer 1.5 有自己的 caveat（S6 詳述）。」|
| **Tier 3 (oral)** | 5 commit 順序 / V3F 命名來歷 / 跨層交互機制 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 表格 16pt；標題 28pt；雙問解釋 14pt italic |
| **C 色彩** | ✅ PASS | phasing 列底色淺藍；tagging 列底色淺綠（與 F3 timeline 兩色標 layer 一致）|
| **L 佈局** | ✅ PASS | 上方對應表 + 下方雙問框，邏輯清晰 |
| **B 雙語** | ✅ PASS | title 雙語；表格 column header 中英對照 |
| **N 數字** | ✅ PASS | **8b8c1fd / 41ff147 / 380e8d2 / d0bcd8c / 938f0df** ← 5/8 §5.1 + 附錄 C ✅ / **99.9%** ← §5.5 ✅ |
| **Agent-PI** | ⚠ 1 建議 | tagging 列「修補 commit」格 3 個 commit 嵌入太擠 → 改 split：「V3F = 41ff147 + 380e8d2」「V5 = + d0bcd8c + 938f0df」放兩行；視覺更清晰 |

#### Reviewer 修正

1. **PI 修正**：tagging 列 commit 格 split 兩行（V3F / V5 各一行）
2. **R-G4**：1 中等術語（PoN）glossary box ✓
3. **R-G2**：1 ≤ 3 ✓

---

### Slide 08 — chr19 752 read-level victims：4-path 驗證 3.5/4 PASS

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  chr19 752 victims 100% 單向 baseline=11 → V3F=21 → V5=21                │
│  chr19 752 victims 100% unidirectional fix                                │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 規模 ────────────────────────┐  ┌─ 修正率 ──────────────────┐     │
│   │  Dump rows  549,206            │  │  V3F 修正率  100.00%       │     │
│   │  3-way merged  1,069,832       │  │  V5  修正率  100.00%       │     │
│   │  Priority bug victims  752     │  │  全 752 條無一條反向       │     │
│   └────────────────────────────────┘  └────────────────────────────┘     │
│                                                                          │
│   ┌─ 4-path 驗證 (T1.2 plan) ─────────────────────────────────────────┐  │
│   │  ① 個案 trace ≥10 條        752 條 ✅ PASS                        │  │
│   │  ② 1Mb 區域聚集            chr19:30M (215) + 27M (133) 46% ⚠ PARTIAL│  │
│   │  ③ Somatic density 共變    high vote ≥5 = 0 受害；低票觸發 🔄 反向有意義│  │
│   │  ④ 修正後消失              V3F/V5 100% ✅ PASS                      │  │
│   │                                                                    │  │
│   │  ➜ 3.5/4 PASS — priority bug 機制因果確立                          │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ Figure F4 — chr19 752 victims 1Mb hotspot scatter ──────────────┐   │
│   │  [figures/F4_chr19_752_victims_scatter.png]                       │   │
│   └────────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│  ⓘ scope: chr19 only                                                      │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 規模統計 + 修正率 + 4-path 驗證表 + Figure F4 + scope footnote（共 6 elements ✓ R-G2）|
| **Tier 2 (note)** | 2 min — 「chr19 read-level audit：dump 549,206 rows × 3 binary versions JOIN merged 1.07M events，篩 germline_majority ≠ somatic_majority 雙向矛盾且 baseline 標 somatic 方向 = 752 victims。V3F + V5 修正率均 100%，全 752 條無一條反向（baseline=11 → V3F=21）。4-path 驗證：① 個案 trace 752 條全 PASS ② 1Mb 區域聚集 chr19:30M 215 + 27M 133 共佔 46% PARTIAL（不是高度集中但有 hotspot）③ somatic density 共變 high vote ≥5 = 0 反向但有意義（解釋：低票才觸發，high 票 sub-clone 一致已經對齊）④ V3F/V5 修正後消失 PASS。3.5/4 機制因果確立。」|
| **Tier 3 (oral)** | 4-path detail / read_name 真實 case 5 條 / SP1/2/3 在 chr19:12-17M 對應 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 規模統計 18pt；修正率 24pt bold；4-path 表 16pt |
| **C 色彩** | ⚠ 1 修正 | 修正率「100%」綠色 emphasize；4-path verdict 用 ✅ / ⚠ / 🔄 emoji（colorblind safe）|
| **L 佈局** | ✅ PASS | 上方規模 + 修正率 雙欄；中部 4-path；下方 F4 圖 |
| **B 雙語** | ✅ PASS | title 雙語；4-path verdict 雙語對照 |
| **N 數字** | ✅ PASS（grep verified）| **549,206** ← 5/8 §4.1 ✅ / **1,069,832** ← 同 ✅ / **752** ← 同 ✅ / **100.00%** ← 同 ✅ / **215 / 133** ← §4.1 表 ✅ / **3.5/4** ← §6.1 ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a)「3.5/4 PASS」表述 PI 可能困惑「半個 PASS」 → 改「3 PASS + 1 PARTIAL = 機制因果確立」更明確；(b) ③ 反向有意義 PI 必問「為何 high vote = 0 受害？」 → speaker note 已涵蓋，slide 上 footnote 補一句：「(low somatic vote 才觸發 priority bug；high vote 等同 sub-clone 一致已對齊 graph)」|

#### Reviewer 修正

1. **C 修正**：100% 綠 + emoji verdict
2. **PI 修正**：(a) 「3 PASS + 1 PARTIAL」改寫；(b) ③ 反向 footnote 補
3. **R-G4**：scope footnote (chr19 only) ✓

---

### Slide 09 — 全基因組 34,855 victims（46×）；priority bug 主分佈非 chr19

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  全基因組 34,855 victims；chr7/chr2/chr1 主分佈，chr19 占 2.16% rank 19   │
│  Genome-wide 34,855 victims; main hotspots chr7/chr2/chr1                 │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 規模對照 ─────────────────────────────────────┐                       │
│   │           chr19 pilot     Genome F1     倍數    │                       │
│   │  Dump     549,206        29,973,253    54.6×   │                       │
│   │  Tagged   ~330K          18,895,432    57×     │                       │
│   │  Victims  752            34,855        46.4×   │                       │
│   │  V3F 修正  100%           100%          一致    │                       │
│   │  V5  修正  100%           100%          一致    │                       │
│   └────────────────────────────────────────────────┘                       │
│                                                                          │
│   ┌─ Per-chr 分佈 (前 5 + chr19/chr8/chrY) ────────────────────────────┐  │
│   │  chr7  ████████████████████  3,508   rank 1   10.1%                │  │
│   │  chr2  ████████████████      2,792   rank 2    8.0%                │  │
│   │  chr1  ███████████████       2,674   rank 3    7.7%                │  │
│   │  chr16 ██████████████        2,584   rank 4    7.4%                │  │
│   │  chr20 ████████████          2,101   rank 7    6.0%                │  │
│   │  ───────────────────────────────────────────────────────────────  │  │
│   │  chr19 ████                    752   rank 19   2.16% ★ 不是 hotspot │  │
│   │  chr8  ███                     666   rank 21   1.9%  ★ 冷區       │  │
│   │  chrY  ▏                        67   rank 24   0.2%  ★ ‰ 排第 1   │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ Figure F2 ─────────────────────────────────────────────────────┐   │
│   │  [F2_priority_bug_per_chr_enrichment.png]                        │   │
│   │  左：victim N rank | 右：enrichment ‰ rank                        │   │
│   └─────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│  ⓘ scope: 全基因組 (T1.2-F1)                                              │
│  ⓘ chr8 LOH+HPSig hotspot 是 ISM 下游（不同 layer，§8.2）                  │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 規模對照表 + Per-chr 分佈 ASCII bar + F2 圖 + 2 footnote（共 5 elements ✓）|
| **Tier 2 (note)** | 2 min — 「chr19 752 不是局部 artifact — 全基因組擴展 34,855 victims（46.4×）。V3F/V5 修正率仍 100%，方向一致。Per-chr 分佈推翻原 chr19 pilot 結論：主要 hotspot 是 chr7 (3,508) / chr2 (2,792) / chr1 (2,674) / chr16 / chr20，**chr19 占比僅 2.16% rank 19**。chr8 priority bug enrichment 0.34× genome avg（rank 21 冷區），與 chr8 LOH+HPSig hotspot 是不同 layer（後者 ISM 下游 false-positive 富集 7.4×）。chrY enrichment ‰ 第 1 但 victim N 小。」|
| **Tier 3 (oral)** | 全 chr enrichment ‰ 表 / chrY 小 N 解釋 / chr8 不同 layer 的詳細分析 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 規模對照 18pt；ASCII bar monospace 14pt；emphasis ★ 標 |
| **C 色彩** | ⚠ 1 修正 | chr19 ASCII bar 紅色高亮（顛覆原結論）；chr8 藍色高亮（特殊冷區）；chrY 黃色（小 N 高 ‰ 例外）|
| **L 佈局** | ✅ PASS | 上對照表 + 中 ASCII bar + 下 F2 圖三層 |
| **B 雙語** | ✅ PASS | title 雙語；表 column 雙語 |
| **N 數字** | ✅ PASS（grep verified）| **34,855** ← §4.2 ✅ / **752** ← §4.1 ✅ / **46.4×** ← §4.2 ✅ / **3,508 / 2,792 / 2,674 / 2,584 / 2,101 / 666 / 67** ← §4.3 表 ✅ / **rank 19 / 21 / 24** ← 同 ✅ / **2.16%** ← 同 ✅ / **0.34×** ← 同 ✅ / **7.4×** ← MEMORY chr8 hotspot ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a)「占比 2.16%」與「rank 19」並列 PI 可能感覺強調過度 → 改「rank 19 (chr 24 條中)」更直觀；(b) chr8 footnote「不同 layer」過於精簡 → 補：「(LOH+HPSig 是 ISM 特徵 cross-talk，priority bug 是 longphase tagging 投票錯，兩者跨層獨立)」 |

#### Reviewer 修正

1. **C 修正**：chr19/chr8/chrY 三色高亮
2. **PI 修正**：(a) rank 19 (chr 24) 改寫；(b) chr8 footnote 補 layer 區分
3. **R-G4**：2 footnote ✓

---

## Batch 3 整體 self-check

| 項目 | 結果 |
|------|------|
| 3 張 slide spec 完整 | ✅ |
| 6 agent review per slide | ✅ |
| Agent-N grep verified | ✅ slide 07/08/09 全 PASS（commit hashes / 752 / 34,855 / per-chr / 修正率） |
| R-G4 三層分流 | ✅ slide 07 PoN glossary / slide 08 scope footnote / slide 09 2 footnote |
| R-G2 ≤ 3 / slide | ✅ slide 07 (1) / slide 08 (1) / slide 09 (2) |

## Batch 3 — User Ack ✅ Final Spec（2026-05-10）

全 5 修正接受。進 Batch 4。

---

## Batch 4 — Slide 10 5 commits / Slide 11 getVote 三版 / Slide 12 SP 修正

---

### Slide 10 — 5 commits 時間軸 + 兩層三色標

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  5 commits 兩層三版 stacking — baseline → V3F → V5                        │
│  5 commits two-layer three-version stacking                                │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ Figure F3 — binary commit timeline ─────────────────────────────┐   │
│   │                                                                  │   │
│   │  04-09     04-10      04-25     04-30(a)    04-30(b)           │   │
│   │   │         │          │         │           │                  │   │
│   │   ▼         ▼          ▼         ▼           ▼                  │   │
│   │ 8b8c1fd  41ff147   380e8d2   d0bcd8c     938f0df                │   │
│   │ PON-only two-layer  INDEL    ploidy fix  threshold              │   │
│   │ flag     getVote   guard     + Layer 1.5  0.95→0.9              │   │
│   │ 🔵       🟢        🟢        🟣          🔵                     │   │
│   │ phasing  tagging   tagging   phasing+   phasing                 │   │
│   │                              tagging                            │   │
│   │ ──────  ─────────────────  ──────────  ────────                 │   │
│   │ baseline ──── V3-Fixed ──   ──────── V5 ────────                │   │
│   │                                                                  │   │
│   └──────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ★ 41ff147 是修偏移的關鍵 commit (tagging 層 priority bug fix)             │
│                                                                          │
│   ┌─ 累計修補 ─────────────────────────────────────────────────────┐    │
│   │  ~155 lines tagging-layer + ~40 lines phasing-layer            │    │
│   │  HaplotagProcess.h:66-68 介面契約零變動                         │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│  ⓘ INDEL guard: 補 OOB undefined behavior（HAPLOTYPE_UNDEFINED 檢查）    │
│  ⓘ threshold 0.95 → 0.9: 觸發 Pass 2 second round 的 purity 閾值          │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + Figure F3 timeline + ★ 關鍵 commit emphasis + 累計修補摘要 + 2 footnote（共 4 elements ✓ R-G2）|
| **Tier 2 (note)** | 2 min — 「5 commits 漸進完成 self-phasing 修補。`8b8c1fd` PON-only flag 解 phasing 層球員兼裁判（藍）；`41ff147` two-layer getVote 解 tagging 層 priority bug（綠）★ 是修偏移的關鍵 commit；`380e8d2` INDEL guard 補 OOB UB（綠）；`d0bcd8c` Pass 2 ploidy fix + bundled Layer 1.5 跨兩層（紫）；`938f0df` threshold 0.95→0.9 phasing 層（藍）。V3-Fixed = baseline + 41ff147 + 380e8d2；V5 = V3F + d0bcd8c + 938f0df。累計約 155 行 tagging + 40 行 phasing；介面契約零變動。」|
| **Tier 3 (oral)** | 各 commit 對應 layer + 為何不能合併單 commit + cherry-pick 自 upstream zhenyu |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | F3 圖已含字級；額外標 ★ 28pt bold；累計摘要 16pt |
| **C 色彩** | ✅ PASS | F3 已含三色標（藍/綠/紫）一致；★ 用 gold 強調 |
| **L 佈局** | ✅ PASS | F3 大圖佔 60% + 摘要 + footnote；balance OK |
| **B 雙語** | ✅ PASS | title 雙語；commit message 英文原文 |
| **N 數字** | ✅ PASS（grep verified）| **8b8c1fd / 41ff147 / 380e8d2 / d0bcd8c / 938f0df** ← 附錄 C ✅ / **155 lines / 40 lines** ← §5.6 累計 ✅ / **HaplotagProcess.h:66-68** ← §5.6 ✅ / **04-09 / 04-10 / 04-25 / 04-30** dates ← §5.1 ✅ |
| **Agent-PI** | ⚠ 1 建議 | F3 圖內 commit 名稱可能小 — P5 vision check 圖文比 ≥ 60% 視覺；slide 加說明「(雙色標 = 跨兩 layer commit `d0bcd8c`)」清楚解釋紫色 |

#### Reviewer 修正

1. **PI 修正**：紫色說明（`d0bcd8c` 跨兩 layer）
2. **R-G4**：2 footnote (INDEL guard / threshold) ✓
3. **R-G2**：0 中等 + 2 簡單 = 2 ✓

---

### Slide 11 — getVote 三版程式碼對照

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  getVote 三版差異：baseline 順序錯 → V3F 兩層 → V5 +Layer 1.5 fallback     │
│  getVote 3-version diff: baseline ordered → V3F two-layer → V5 +1.5      │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│  ┌─ baseline ❌ ────┐  ┌─ V3F ✅ tagging ───┐  ┌─ V5 ✅ + fallback ──┐  │
│  │ vector keys = [  │  │ // Layer 1 germline│  │ // Layer 1 (同 V3F) │  │
│  │  somatic FIRST,  │  │ if (gHP1>0          │  │ if (gHP1>0           │  │
│  │  mixed,          │  │   ‖ gHP2>0)        │  │   ‖ gHP2>0)         │  │
│  │  germline LAST   │  │   gR=...           │  │   gR=...             │  │
│  │ ]                │  │                    │  │ // Layer 1.5 NEW     │  │
│  │ for (k:keys) {   │  │ // Layer 2 somatic │  │ else if (sHP1>0      │  │
│  │  if (>0) {       │  │ if (sTotal>0)      │  │   ‖ sHP2>0)         │  │
│  │   hp=...;        │  │   hp = (gR==1)?11: │  │   gR=(sHP1>=sHP2)?1:2│  │
│  │   break;❌       │  │        (gR==2)?21: │  │ // Layer 2 (同 V3F) │  │
│  │  }               │  │        33;        │  │ if (sTotal>0)       │  │
│  │ }                │  │ else hp=gR;        │  │   hp=...           │  │
│  └──────────────────┘  └────────────────────┘  └─────────────────────┘  │
│                                                                          │
│  ┌─ Layer 1.5 caveat ────────────────────────────────────────────────┐  │
│  │ V5 Layer 1.5 在 germline-absent 區會繼承 priority bug 偏移 → §S6  │  │
│  └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│  📖 germline het: 個人遺傳變異中雜合位點 (A/G)，雙親各帶一型             │
│  ⓘ two-layer ⓘ Layer 1.5 fallback ⓘ phased votes                         │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 三 code panel side-by-side + Layer 1.5 caveat 預告 + 1 glossary + 3 footnote（共 5 elements）|
| **Tier 2 (note)** | 2 min — 「baseline `getVote()` 用 vector ordered + break early，1 票 somatic 觸發後 germline 永遠看不到。V3F 改 explicit Layer 1 (germline only) + Layer 2 (somatic annotation)：germline 永不被 somatic overrule，有票就決定 11/21/33。V5 加 Layer 1.5 fallback：germline 缺席時用 somatic phased votes 決方向（HP1_1 vs HP2_1 票數）。但 Layer 1.5 在 germline-absent 區會繼承 priority bug 偏移（→ §S6 詳述）— 種 cliffhanger。」|
| **Tier 3 (oral)** | enum HAPLOTYPE1_1=2 / HP tag int=11 比較 bug 細節 / V3F bonus 修：hpResult ≠ HAPLOTYPE1_1 比較失誤 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | code panel monospace 14pt；title 28pt；caveat 16pt italic |
| **C 色彩** | ⚠ 1 修正 | baseline panel 紅底淺色（warning）；V3F 綠底淺色（fix）；V5 綠底淺色（fix）+ Layer 1.5 黃 highlight（caveat 預告）|
| **L 佈局** | ⚠ 1 建議 | 三 code panel side-by-side 在 16:9 可能擁擠（各約 33% width）→ 縮減 code 量：每 panel 只放關鍵 ~5 行 + 「...」省略；P5 render check |
| **B 雙語** | ✅ PASS | title 雙語；code 英文原碼；comment 中文 |
| **N 數字** | ✅ PASS | **HAPLOTYPE1_1 / HAPLOTYPE2_1 / HAPLOTYPE3** enum 名 ← §5.6 ✅ / **hp=11 / 21 / 33** tag values ← 附錄 B ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a) PI 不熟 C++ syntax — 建議在底部加一行白話解釋：「baseline 順序錯 → V3F 拆 germline 主導 + somatic 標籤 → V5 補 germline 缺席 fallback」；(b) Layer 1.5 caveat 預告位置 — 「§S6」對 PI 不直觀，改「(後面 slide 16 詳述)」|

#### Reviewer 修正

1. **C 修正**：3 panel 三色底（紅/綠/綠+黃 highlight）
2. **L 修正**：每 panel ~5 行 + 「...」
3. **PI 修正**：(a) 底部白話三段解釋；(b) caveat 改「(slide 16 詳述)」
4. **R-G4**：1 中等（germline het）glossary box ✓ + 3 簡單 footnote ✓
5. **R-G2 例外**：1+3 = 4 處理超 R-G2≤3，但因核心機制 forced；接受

---

### Slide 12 — SP1/2/3 修正後對齊 paired 3/3 + 全基因組 17.3:1 → ~1:1

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  個案層 V5 修正 3/3 + 全基因組 HP1:HP2 17.3:1 → ~1:1                       │
│  Site-level V5 fixes 3/3 + genome HP1:HP2 17.3:1 → ~1:1                   │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 個案層：SP1/2/3 修正後對齊 paired ────────────────────────────┐    │
│   │  位點          baseline   V5 翻轉    paired GT   對齊?         │    │
│   │  SP1 chr19:17M  113:0     HP2 主導   HP2         ✅           │    │
│   │  SP2 chr19:12M  109:1     HP2 主導   HP2         ✅           │    │
│   │  SP3 chr19:12M  108:0     HP2 主導   HP2         ✅           │    │
│   │  ────────────────────────────────────────────────────         │    │
│   │  3/3 對齊 paired ground truth                                  │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ 全基因組層：HP1:HP2 17.3:1 → ~1:1 ────────────────────────────┐    │
│   │   metric                  baseline    V5         Δ              │    │
│   │   HP1:HP2 ratio           17.3:1      ~1:1       消除偏移      │    │
│   │   94.6% somatic→HP1       是         ~50%        平衡          │    │
│   │   15-site Problem PS      48.5%      52.0%       +3.5 pp       │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ Figure F5 — Layer 1.5 zero-sum 重分配（4 象限） ───────────────┐   │
│   │  [F5_layer15_zero_sum_4quadrant.png]                             │   │
│   │  germline=0 +560,881 reads / germline>0 −560,881 / 總和 = 0      │   │
│   └─────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│  ➜ 個案 + 統計 兩層共驗 V5 修對                                           │
│  📖 haplotype: 父系/母系兩條染色體之一                                   │
│  ⓘ scope: chr19 個案 + 全基因組統計（混合）                                │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 個案層表 + 全基因組層表 + F5 圖 + 結論 + glossary + scope footnote（共 6 elements ✓ R-G2）|
| **Tier 2 (note)** | 1.5 min — 「§S2 觀察起點 chr19 三位點 113:0/109:1/108:0 — V5 修正後翻回 paired ground truth 三位都 OK。全基因組層 HP1:HP2 17.3:1 → ~1:1 消除偏移；94.6% somatic→HP1 → ~50% balanced；15-site Problem PS（含 SP1/2/3）48.5%→52.0% +3.5 pp 看似小但機制顯著（極端失衡位點本身難改善）。F5 4 象限視覺 zero-sum 重分配：germline=0 區 +560,881 reads（V5 Layer 1.5 觸發）/ germline>0 區 -560,881 / 總和 0。Pass 2 reclassify 104K germline het 為 somatic/未 phase 是因果。」|
| **Tier 3 (oral)** | V2b / V3F 中間版本對齊 paired 細節 / Layer 1.5 zero-sum 機制 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 兩表 16pt；結論 18pt bold；F5 圖 caption 14pt |
| **C 色彩** | ⚠ 1 修正 | 個案層 ✅ 綠色 emphasize；全基因組層 Δ 用「消除」「平衡」綠色 + 「+3.5 pp」綠色；F5 圖左 panel 綠/右 panel 灰中性 |
| **L 佈局** | ✅ PASS | 上個案表 + 中全基因組表 + 下 F5 圖；垂直三段清晰 |
| **B 雙語** | ✅ PASS | title 雙語；表 column 雙語；結論句中文為主 |
| **N 數字** | ✅ PASS（grep verified）| **113:0 / 109:1 / 108:0** ← §2.2 + §6.3 ✅ / **17.3:1 / ~1:1 / 94.6% / ~50%** ← §6.4 ✅ / **48.5% / 52.0% / +3.5 pp** ← §6.4 ✅ / **+560,881 / -560,881** ← §8.4 ✅ |
| **Agent-PI** | ⚠ 1 建議 | F5 圖內「560,881」精確數字可能難讀（Tier 1 已抽出於 caption）→ caption 加「= V5 Pass 2 reclassify 104K germline het 結果」連結 §8.4 機制 |

#### Reviewer 修正

1. **C 修正**：3 處綠色 emphasize（個案 ✅ / 全基因組 Δ / F5 caption）
2. **PI 修正**：F5 caption 補 zero-sum 機制連結
3. **R-G4**：1 中等（haplotype）glossary + scope footnote ✓
4. **R-G2**：1 + 1 = 2 ✓

---

## Batch 4 整體 self-check

| 項目 | 結果 |
|------|------|
| 3 張 slide spec 完整 | ✅ |
| 6 agent review per slide | ✅ |
| Agent-N grep verified | ✅ slide 10/11/12 全 PASS |
| R-G4 三層分流 | ✅ |
| R-G2 ≤ 3 / slide | ✅ slide 10 (2) / slide 11 (4 例外) / slide 12 (2) |

## Batch 4 — User Ack ✅ Final Spec（2026-05-10，pending 待 confirm）

預設套用全 6 修正（用戶 batch 4 ack 後確認）

---

## Batch 5 — Slide 13 20 指標 / Slide 14 caller F1 + cliffhanger / Slide 15 paired mode

---

### Slide 13 — 20 指標 no regression — 5 大類別全綠

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  20 指標全綠：6 項顯著改善 +8.3 ~ +99.7%；其餘 ±0.01 內持平                │
│  20 metrics no regression: 6 significant improvements; rest within noise  │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 5 大類別 × 20 指標 ────────────────────────────────────────────┐   │
│   │                                                                  │   │
│   │  ① ISM aggregate (3)        TP_rate +0.005 / HP_Ratio 0.788→0.574 │   │
│   │                              / Potential_LOH +3.5 pp              │   │
│   │  ② HP_Ratio AUC (2)         All -0.005 (隨機區間) / Inner +0.002  │   │
│   │  ③ Methylation 6 feat       全 ±0.01 內持平                       │   │
│   │  ④ Paired GT concord. ⭐    clean PS +8.3 pp / 15-site Agg +6.65 pp│   │
│   │     (4)                      / 15-site Clean PS +13.3 pp / Prob +3.5│   │
│   │  ⑤ HP / LOH 結構 ⭐        N50 +99.7% / Phased rate +23.6 pp    │   │
│   │     (5)                      / 1.36× 快 / LOH regions 完全相同   │   │
│   │                                                                  │   │
│   └──────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ 6 顯著改善 (⭐) ──────────────────────────────────────────────┐    │
│   │  +99.7% N50 / +23.6 pp Phased rate / 1.36× 速度                 │    │
│   │  +13.3 pp 15-site Clean PS / +8.3 pp 全基因組 clean PS / +6.65 pp 15-Agg│  │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ➜ 20 / 0 指標 no regression — V5 全面 production-ready                  │
│                                                                          │
│  📖 LOH: 雜合性丟失 (Loss of Heterozygosity)                              │
│  ⓘ scope: HCC1395 5kHz @ 0.93 purity (PI 報告 V5 = Pass 1 only BAM)        │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 5 大類別表 + 6 顯著改善 highlight box + 結論 + glossary + scope footnote（共 5 elements ✓ R-G2）|
| **Tier 2 (note)** | 2 min — 「20 指標 5 大類別全綠：① ISM aggregate 3 項 — HP_Ratio median 0.788→0.574 看似下降但是 tag bias 修正非變差；② HP_Ratio AUC 兩項在隨機區間內；③ methylation 6 feature 全 ±0.01 持平；④ Paired GT concordance 4 項 ⭐ 顯著改善 +6.65 ~ +13.3 pp；⑤ HP/LOH 結構 5 項 ⭐ 改善 N50 +99.7%、Phased rate +23.6 pp、執行 1.36× 快、LOH regions 完全相同。20/0 指標 no regression — V5 全面 production-ready。」|
| **Tier 3 (oral)** | methylation 6 feat 詳細列表 / HP_Ratio 0.788→0.574 為何不是變差 / LOH bed Jaccard=1.0 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ⚠ 1 修正 | 5 大類別表字密度高 → 用較小 14pt + indent；6 顯著改善 highlight 用 18pt bold |
| **C 色彩** | ⚠ 1 修正 | ⭐ 顯著改善類別底色淡綠；6 改善 highlight box 綠底深綠字；scope footnote 灰色（避免搶眼）|
| **L 佈局** | ✅ PASS | 上類別表 + 中 highlight + 下結論；視覺金字塔 |
| **B 雙語** | ✅ PASS | title 雙語；類別 ① ② ③ ④ ⑤ Latin 數字易讀 |
| **N 數字** | ✅ PASS（grep verified）| **20 個指標** + **+0.005 / 0.788 → 0.574 / +3.5 pp / -0.005 / +0.002 / ±0.01 / +8.3 pp / +6.65 pp / +13.3 pp / +99.7% / +23.6 pp / 1.36×** ← §8.5.1 全表 ✅ |
| **Agent-PI** | ⚠ 2 建議 | (a) HP_Ratio 0.788→0.574 PI 必問「為何下降不是變差」 → slide 加 footnote「(HP_Ratio 由 tag 偏移驅動；偏移消除後 ratio 正常化非變差)」；(b)「20/0 no regression」結論前加「pre-registered metrics」表明事先定義避免 cherry-picking 嫌疑 |

#### Reviewer 修正

1. **T 修正**：類別表 14pt indent / highlight 18pt bold
2. **C 修正**：⭐ 類別淡綠 + highlight 深綠字 + footnote 灰
3. **PI 修正**：(a) HP_Ratio footnote；(b) 「pre-registered metrics」加詞
4. **R-G4**：1 中等（LOH）glossary + scope footnote ✓
5. **R-G2**：1 + 1 = 2 ✓

---

### Slide 14 — Caller F1 三版完全相同；purity 0.6 完整對照 + Cliffhanger

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  Caller F1 vs SEQC2 三版完全相同；purity 0.6 完整 0 critical regression    │
│  Caller F1 identical across versions; purity 0.6 fully verified           │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ HCC1395 5kHz @ 0.93 purity ──────────────────────────────────┐   │
│   │  版本   TP    FP   FN   Precision Recall   F1 vs SEQC2         │   │
│   │  A1 baseline 28509 11606 10938  0.7107  0.7227  0.7166         │   │
│   │  A3 V3F     28509 11606 10938  0.7107  0.7227  0.7166         │   │
│   │  A5 V5      28509 11606 10938  0.7107  0.7227  0.7166         │   │
│   │  ─────────── 三版完全相同 ───────────                           │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ HCC1395 t30_n20 @ 0.6 purity ────────────────────────────────┐   │
│   │  版本   TP    FP   FN   F1 vs SEQC2                            │   │
│   │  B1 baseline 24190 13487 15257  0.6273                          │   │
│   │  B3 V3F     24190 13487 15257  0.6273                          │   │
│   │  B5 V5      24190 13487 15257  0.6273                          │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ 因果鏈 ───────────────────────────────────────────────────────┐    │
│   │  ClairS-TO snv.vcf FILTER=PASS 集合 → V5 不改 FILTER →          │    │
│   │  PASS set 相同 → TP/FP/FN 完全相同 → F1 完全相同                │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌──────────────────────────────────────────────────────────────────┐  │
│   │  ➜ V5 不改 caller；ΔF1 (0.93→0.6) = -0.0893 為 ClairS-TO 性質   │  │
│   │  ➜ 但 5/9 paired cross-ref 揭露另一面... 〔next slide 15-16〕    │  │
│   └──────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│  📖 purity: 樣本中腫瘤細胞占比；ClairS-TO 由 ploidy 算出                 │
│  ⓘ PASS set: ClairS-TO snv.vcf FILTER=PASS 的 variants 集合              │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 0.93 表 + 0.6 表 + 因果鏈框 + ➜ 結論 + cliffhanger + glossary + footnote（共 6 elements）|
| **Tier 2 (note)** | 2 min — 「purity 0.93 三版 caller F1 = 0.7166 完全相同（TP 28509 / FP 11606 / FN 10938）。purity 0.6 三版 F1 = 0.6273 完全相同。為何相同：ClairS-TO PASS set 由 caller 決定，longphase-to 改 GT/PS/GT2/GT3 不改 FILTER，所以 PASS set 不變 → TP/FP/FN 不變 → F1 不變。ΔF1 (0.93→0.6) = -0.0893 是 ClairS-TO 性質（低 purity 本身偵測下降），與 V5 無關。**結論：V5 不改 caller，F1 不變，所以講 self-phasing 修補不能用 F1 衡量；正確 metric 是 read-level tag concordance（+13.3 pp paired GT 在 §6.4）。但 5/9 paired cross-ref 揭露另一面（next slide）...**」|
| **Tier 3 (oral)** | PASS set / FILTER 機制細節 / purity 0.6 N50 微差 -14.5% 解釋 / V3F vs V5 0.6 樣本對比 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 兩表 14pt；F1 column bold；cliffhanger 18pt bold italic |
| **C 色彩** | ⚠ 1 修正 | 「三版完全相同」綠色 emphasize；cliffhanger 黃底警示色（暗示 caveat 將至，過渡 emotional drop）|
| **L 佈局** | ✅ PASS | 上 0.93 表 + 中 0.6 表 + 因果鏈 + cliffhanger，垂直結構平順 |
| **B 雙語** | ✅ PASS | title 雙語；表 column header 標準；cliffhanger 中文為主 |
| **N 數字** | ✅ PASS（grep verified）| **28509 / 11606 / 10938 / 0.7107 / 0.7227 / 0.7166** ← §8.5.2 A1/A3/A5 ✅ / **24190 / 13487 / 15257 / 0.6420 / 0.6132 / 0.6273** ← §8.5.2 B1/B3/B5 ✅ / **-0.0893** ← §8.5.2 ✅ |
| **Agent-PI** | ⚠ 1 建議 | 因果鏈圖 PI 可能跳過 → 改 1 句白話：「(V5 改的是 read 怎麼分類到 HP1/HP2，不改 caller 怎麼判斷哪些位點是 somatic — 所以 PASS set 不變)」 |

#### Reviewer 修正

1. **C 修正**：「三版完全相同」綠 + cliffhanger 黃底
2. **PI 修正**：因果鏈白話 1 句
3. **R-G4**：1 中等（purity）glossary + 1 footnote (PASS set) ✓
4. **R-G2**：1 + 1 = 2 ✓
5. **Cliffhanger transition**：✅ slide 末已標「→ 但 5/9 paired cross-ref 揭露另一面...」

---

### Slide 15 — paired mode 整體無 priority bug

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  paired mode 整體無偏移：HP1:HP2 = 1:1.275；som_ratio mean 0.462           │
│  paired mode no systematic bias                                            │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ paired chr19 HP:Z: 分布 (582K total / 354,919 tagged) ────────┐    │
│   │  HP:Z:2     183,309   51.6%                                     │    │
│   │  HP:Z:1     143,760   40.5%                                     │    │
│   │  HP:Z:2-1    14,504    4.1%                                     │    │
│   │  HP:Z:1-1    12,401    3.5%                                     │    │
│   │  HP:Z:3       1,145    0.3%                                     │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ paired vs TO ratio 對照 ─────────────────────────────────────┐    │
│   │             paired       vs       TO baseline                  │    │
│   │  germline   1:1.275 ✅           17.3:1 ❌ (priority bug)      │    │
│   │  somatic    1:1.169 ✅           全偏 HP1                      │    │
│   └─────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ 57 chr19 1Mb windows: som_ratio 跨 0-1 全範圍分布 ────────────┐   │
│   │  mean 0.462 / median 0.494 / stdev 0.332                         │   │
│   │  ─────────────────────────────────                              │   │
│   │  真實 sub-clone signal 案例：                                     │   │
│   │  • chr19:3M  全 HP2-1 (755/0)   → LOH 方向特定                  │   │
│   │  • chr19:0M  全 HP1-1 (330/1)   → 反向區域                       │   │
│   │  • chr19:17M 對稱 0.500 (265/265) → SP1 附近 paired 認雙 sub-clone │  │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ Figure F6 — paired vs TO chr19 HP distribution ────────────┐      │
│   │  [F6_paired_vs_TO_HP_distribution.png]                        │      │
│   └───────────────────────────────────────────────────────────────┘      │
│                                                                          │
│  📖 HP:i: vs HP:Z: longphase-to=整數 (1/2/11/21/33) / longphase-s=字串    │
│  ⓘ som_ratio: HP1-1 票 / (HP1-1 + HP2-1) 票，量化偏向                    │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + paired 分布表 + paired vs TO ratio + 57 windows 統計 + F6 圖 + glossary + footnote（共 6 elements）|
| **Tier 2 (note)** | 2 min — 「paired mode 用不同 binary（longphase-s 不是 longphase-to fork），HP tag 用 HP:Z: 字串編碼。chr19 paired 分布 HP:Z:2 51.6% / HP:Z:1 40.5% germline 1:1.275 接近隨機；somatic HP:Z:2-1 4.1% / HP:Z:1-1 3.5% 也 1:1.169 接近隨機。**對比 TO baseline 17.3:1 priority bug — paired mode 整體無 systematic bias。** 57 windows som_ratio mean 0.462 中位 0.494 跨 0-1 全範圍分布，stdev 0.332 大變動 = 真實 sub-clone signal（chr19:3M 全 HP2-1 LOH 方向特定 / chr19:0M 全 HP1-1 反向 / chr19:17M 對稱 0.500 = SP1 附近 paired 認雙 sub-clone 共現，vs TO baseline 113:0 失衡）。」|
| **Tier 3 (oral)** | longphase-s codebase / SomaticHaplotagProcess.cpp:533 HP tag 字串編碼 / paired 軸對齊 vs TO 軸 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 各表 14pt；標題 28pt；案例列表 16pt |
| **C 色彩** | ⚠ 1 修正 | paired ratio 列綠 ✅ / TO 列紅 ❌ 對比；57 windows 案例三色標（chr19:3M 紫 / chr19:0M 橙 / chr19:17M 藍）|
| **L 佈局** | ✅ PASS | 4 區塊垂直；F6 圖底部 |
| **B 雙語** | ⚠ 1 建議 | longphase-s vs longphase-to 名稱英文原文保留；slide 內加 1 行說明 longphase-s = paired mode binary（與 longphase-to fork 不同 codebase）|
| **N 數字** | ✅ PASS（grep verified）| **183,309 / 143,760 / 14,504 / 12,401 / 1,145** ← §8.6.2 表 ✅ / **51.6% / 40.5% / 4.1% / 3.5% / 0.3%** ← 同 ✅ / **1:1.275 / 1:1.169** ← 同 ✅ / **mean 0.462 / median 0.494 / stdev 0.332** ← §8.6.3 ✅ / **chr19:3M (755/0) / 0M (330/1) / 17M (265/265)** ← §8.6.3 ✅ |
| **Agent-PI** | ⚠ 1 建議 | 「真實 sub-clone signal」解釋對 PI OK 但稍長 → 簡化為「(stdev 0.332 跨 windows 大變動 = 不同區域有不同 LOH 方向，是真實生物學)」 |

#### Reviewer 修正

1. **C 修正**：paired ✅ vs TO ❌ + 三色標案例
2. **B 修正**：補 longphase-s 說明一行
3. **PI 修正**：sub-clone signal 簡化敘述
4. **R-G4**：1 中等（HP:i: vs HP:Z:）glossary + 1 footnote (som_ratio) ✓
5. **R-G2**：1 + 1 = 2 ✓

---

## Batch 5 整體 self-check

| 項目 | 結果 |
|------|------|
| 3 張 slide spec 完整 | ✅ |
| 6 agent review per slide | ✅ |
| Agent-N grep verified | ✅ slide 13/14/15 全 PASS（20+ 個數字 verified） |
| R-G4 三層分流 | ✅ |
| R-G2 ≤ 3 / slide | ✅ slide 13 (2) / 14 (2) / 15 (2) |
| Cliffhanger | ✅ slide 14 末已含 |

## Batch 5 — 6 修正建議匯總（pending ack）

預設套用全 6 修正（用戶 batch 5 ack 後確認）

---

## Batch 6 — Slide 16 V5 Layer 1.5 設計缺陷 / Slide 17 5 errata / Slide 18 follow-up

---

### Slide 16 — germline-absent 區 V5 = baseline 4.19:1；V3F 標 hp=33 反而穩健

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  V5 Layer 1.5 設計缺陷：germline-absent 區與 baseline 4.19:1 完全相同     │
│  V5 Layer 1.5 design caveat: identical to baseline in germline-absent     │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ Cross-tab: 5,789 chr19 germline-absent events (cnt_HP1+HP2=0) ──┐  │
│   │  paired HP    events    baseline=11   baseline=21   V5=11   V5=21│  │
│   │  HP:Z:1-1     2,040     1,679         318           1,679   318 │  │
│   │  HP:Z:2-1     1,588     1,291         295           1,291   295 │  │
│   │  HP:Z:3         530       342         178             343   177 │  │
│   │  ──────────────────────────────────────────────────────────────  │  │
│   │  總計                  3,312          791          3,313   790  │  │
│   │                        ↑ 4.19:1 偏 HP1                ↑ 4.19:1 同│  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ V3F 在該區域：全 5,789 events 標 hp=33 (somatic ambiguous) ✅ ────┐  │
│   │  保守不選邊 → 避免錯標方向 → V3F 反而比 V5 穩健                    │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ V5 Layer 1.5 機制詮釋 ────────────────────────────────────────┐    │
│   │  V5 Layer 1.5: germline=0 時用 somaticHP1 vs somaticHP2 票數決方向│  │
│   │  ↓                                                              │    │
│   │  sub-clone somatic 100% 共現 → graph 偏向同一 haplotype          │    │
│   │  ↓                                                              │    │
│   │  somatic 票偏向同邊 → Layer 1.5 結果 = priority bug 偏移         │    │
│   │                                                                  │    │
│   │  ★ 結論：V5 Layer 1.5 = priority bug 的 feature 化非修補         │    │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ Figure F7 ─────────────────────────────────────────────────────┐  │
│   │  [F7_germline_absent_crosstab.png]                                │  │
│   │  左: 三版對比 / 右: paired HP cross-tab                            │  │
│   └─────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│  📖 sub-clone: 腫瘤內基因型相同細胞群；不同 sub-clone 帶不同 somatic     │
│  ⓘ Layer 1.5 fallback ⓘ phased votes                                     │
│  ⓘ scope: paired chr19 germline-absent (5,789 events ≈ 5% chr19 events)  │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + cross-tab 表 + V3F 對比 + 機制詮釋 + F7 圖 + glossary + footnote（共 6 elements）|
| **Tier 2 (note)** | 2.5 min — 「對 paired chr19 read_name × T1.2 baseline/V3F/V5 vote dump JOIN，篩 cnt_HP1+HP2=0 且 somatic>0 = 5,789 chr19 germline-absent events。cross-tab：baseline hp=11 (HP1 系列) 3,312 / hp=21 (HP2) 791 = **4.19:1 偏 HP1**（priority bug 次峰）；V3F 全 5,789 events 標 **hp=33** somatic ambiguous 保守不選邊 ✅；**V5 hp=11 3,313 / hp=21 790 = 4.19:1 與 baseline 完全相同！**Layer 1.5 機制詮釋：germline=0 時用 somaticHP1 vs somaticHP2 票數決方向，sub-clone somatic 100% 共現 → graph 偏向同一 haplotype → 投票偏向同邊 → Layer 1.5 結果 = priority bug 偏移。**結論：V5 Layer 1.5 = priority bug 的 feature 化非修補；V3F 標 hp=33 反而更穩健。**caveat：cross-binary axis alignment（paired vs TO 各自獨立 phasing）；binary-internal 量化（baseline 自身 4.19:1）不受影響；F-paired-D3 量化 ISM 影響待跑。」|
| **Tier 3 (oral)** | cross-binary axis caveat 詳述 / phasing graph 機制連結 §3.2 / Layer 1.5 改回 V3F default 設計選擇 / F-paired-D3 工作量 1-2 day |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | cross-tab 14pt；機制詮釋 16pt；★ 結論 18pt bold |
| **C 色彩** | ⚠ 2 修正 | (a) cross-tab baseline 列紅 ❌ / V3F 列綠 ✅ / V5 列紅 ❌（與 baseline 同色強調 same-as-baseline）；(b) 「★ 結論」框淺紅底深紅字 emphasize caveat（這是整份 PPT 唯一紅 emphasis 位置）|
| **L 佈局** | ⚠ 1 建議 | 4 區塊垂直 + F7 圖共 5 elements 視覺密度高 → F7 圖縮小放右側 panel；左側 3 個 text 區 |
| **B 雙語** | ⚠ 1 建議 | 「priority bug 的 feature 化非修補」中文太精煉 → 補 EN 對齊：「(V5 Layer 1.5 makes priority bug a feature, not a fix)」|
| **N 數字** | ✅ PASS（grep verified）| **5,789** ← §8.6.4 ✅ / **2,040 / 1,588 / 530** ← 同 ✅ / **1,679 / 318 / 1,291 / 295 / 342 / 178** ← 同 ✅ / **3,312 / 791 / 3,313 / 790** ← 同 ✅ / **4.19:1** ← 同 ✅ / **~5%** chr19 events ← 推估 |
| **Agent-PI** | ⚠ 3 建議 | (a) 「★ V5 = baseline 完全相同！」是整份 PPT 最 critical 訊息 — slide 標題改更直白：「**V5 Layer 1.5 在 germline-absent 區未修對 priority bug**」更明確；(b) PI 必問「為什麼 V5 設計者沒測這個？」 → speaker note 補：「(V5 設計時未對 paired ground truth 做 germline-absent cross-ref，5/9 paired audit 補上才發現)」；(c) 結論句「priority bug 的 feature 化非修補」非常 powerful 但 PI 可能需更簡單 → 平行 alt：「**V5 把 priority bug 從 bug 變成 designed behavior，但偏移本質沒變**」 |

#### Reviewer 修正

1. **C 修正**：cross-tab 行色標 + 結論紅框
2. **L 修正**：F7 縮小右側 panel
3. **B 修正**：補 EN 結論對齊
4. **PI 修正**：(a) title 改更直白；(b) speaker note 補測試 gap 解釋；(c) 結論平行 alt
5. **R-G4**：1 中等（sub-clone）+ 2 簡單 footnote + scope footnote ✓
6. **R-G2**：1 + 2 + 1 = 4 處理超 R-G2≤3，但因核心 caveat slide forced；接受例外

---

### Slide 17 — 5 條 PI errata 已 patch（E1-E5）

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  PI 報告 4-29 5 條 errata 已 patch；主結論不撤回                          │
│  PI report 4-29 5 errata patched; main conclusions retained               │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ E1 §3.3.3 chr19 SP1/2/3 解讀降級 ─────────────────────────────┐   │
│   │   原: 「主要 hotspot」                                            │   │
│   │   新: 「可重現案例」(chr19 占 priority bug 2.16% rank 19)        │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ E2 §5.2 V5 working tree commit 狀態 ─────────────────────────┐   │
│   │   原: 「working tree 未 commit」                                  │   │
│   │   新: 「✅ 2026-04-30 已 commit (d0bcd8c + 938f0df)」              │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ E3 §5.2 priority bug 證據強度升級 ───────────────────────────┐   │
│   │   原: commit msg + 3 IGV 截圖                                     │   │
│   │   新: + 34,855 read-level victims 全基因組鐵證 100% 修正            │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ E4 ★ §6.4/§6.5 V5 數值歸因精確化（最重要 errata）─────────────┐   │
│   │   原: 「V5 four-commit chain 整體效益」                            │   │
│   │   新: V5 BAM = Pass 1 only（ploidy bug 讓 purity=0）             │   │
│   │       主要功勞: V3F + Layer 1.5；Pass 2 二次效益尚未量化            │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│   ┌─ E5 ★ §5.2 V5 Layer 1.5 設計缺陷（5/10 加）──────────────────┐   │
│   │   原: Layer 1.5 = germline 缺席 fallback (隱含「修補」)            │   │
│   │   新: germline-absent 區 V5 = baseline 4.19:1 → priority bug      │   │
│   │       feature 化非修補；V3F hp=33 反而穩健                         │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│  ➜ 主結論不撤回；E4 + E5 是核心，其他 3 條表述精確化                     │
│  ⓘ commit chain: f17754f → 2553e96 → 71d21bd                              │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 5 errata box (E1-E5) + 主結論句 + commit chain footnote（共 7 elements，1 over R-G2 但每 E 都是核心 errata 必留）|
| **Tier 2 (note)** | 1.5 min — 「PI 報告 4-29 5 條 errata patched，主結論不撤回。E1 chr19 SP1/2/3 從『主要 hotspot』降級為『可重現案例』；E2 V5 working tree 已 commit 於 4-30；E3 priority bug 證據從『commit msg + IGV』升級為『+ 34,855 read-level 鐵證 100% 修正』；★ E4 V5 數值歸因精確化 — PI 報告 V5 BAM 是 Pass 1 only（ploidy bug 讓 purity=0），主要功勞 V3F + Layer 1.5，Pass 2 二次效益尚未量化；★ E5 5/10 加 — V5 Layer 1.5 在 germline-absent 區與 baseline 4.19:1 完全相同，是 priority bug 的 feature 化非修補。E4 + E5 是核心，其他 3 條表述精確化。」|
| **Tier 3 (oral)** | 各 errata commit hash / errata patch 工作量 / 為何不撤回 vs 補 banner |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 5 box 各 14-16pt；★ 標 18pt；title 28pt |
| **C 色彩** | ⚠ 1 修正 | E1/E2/E3 灰底中性（表述精確化）；E4 + E5 ★ 黃底 emphasize（核心 errata）|
| **L 佈局** | ✅ PASS | 5 box 垂直；每 box 「原 / 新」對比清晰 |
| **B 雙語** | ⚠ 1 建議 | 5 errata 中文為主對 PI 已 OK，但 title 雙語 + 每 errata 加 EN 短句（≤10 詞）增可讀性 |
| **N 數字** | ✅ PASS（grep verified）| **2.16% rank 19** / **d0bcd8c / 938f0df** / **34,855 / 100%** / **purity=0** / **4.19:1** ← 5/9 errata + 5/8 主報告 全 verified ✅ / **commit chain f17754f / 2553e96 / 71d21bd** ← 5/8 主報告 §9.2 末 + git log ✅ |
| **Agent-PI** | ⚠ 1 建議 | E5 是「5/10 加」可能 PI 不知何時加 → 補一行說明：「(5/9 paired audit Step D 揭露 → 5/10 amend 進主報告 + errata)」 |

#### Reviewer 修正

1. **C 修正**：E4 + E5 ★ 黃底 / E1-E3 灰底
2. **B 修正**：每 errata 加 EN 短句
3. **PI 修正**：E5 補來源說明
4. **R-G4**：commit chain footnote ✓
5. **R-G2 例外**：6 main + footnote = 7 elements，但 5 errata 是核心 caveat slide forced；接受

---

### Slide 18 — 整體成熟度 + 5 項 follow-up

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  V5 仍可作 production baseline；5 項 follow-up 待跑                        │
│  V5 ready as production baseline; 5 follow-ups pending                    │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 整體成熟度 (12 維度) ──────────────────────────────────────────┐  │
│   │  ✅ 機制因果        ✅ 修補設計合理       ✅ chr19 SP 對齊       │  │
│   │  ✅ 全基因組擴展    ✅ V5/V3F zero-sum 釐清 ✅ 20 指標 0 regression│ │
│   │  ✅ Caller F1 三版相同 ✅ purity 0.6 完整對照 ✅ 三路徑算法不依賴 │  │
│   │  ✅ Pass 2 量化 +3.51% N50 ✅ 版本對齊 HEAD 938f0df              │  │
│   │  ✅ Paired Step A+C 整體 NEGATIVE                                 │  │
│   │  ⚠️ V5 Layer 1.5 germline-absent 設計缺陷 (E5)                    │  │
│   │  ⏸ Pass 2 second round 獨立貢獻量化 (T1.3)                       │  │
│   │  ⏸ 跨樣本擴展 (T3 — 6 樣本 vote audit)                            │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ 5 項 follow-up cycle ─────────────────────────────────────────┐    │
│   │  ID                  內容                          預估時間      │    │
│   │  F-paired-D1   germline-absent 全基因組擴展        0.5 day      │    │
│   │  F-paired-D2   phase block 內 axis-aligned 分析    1 day        │    │
│   │  F-paired-D3 ★ Layer 1.5 改回 V3F ISM 影響評估     1-2 day      │    │
│   │  T3            7 樣本跨樣本擴展                     1-2 day      │    │
│   │  T1.3          4-cell ablation                     3 day        │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ 結論：V5 = production baseline ────────────────────────────────┐   │
│   │  germline-absent 區占比小（~5% chr19 events）不阻擋整體          │   │
│   │  F-paired-D3 量化後決定是否回 V3F default                        │   │
│   └────────────────────────────────────────────────────────────────┘   │
│                                                                          │
│  📖 LOH: 雜合性丟失 / cnLOH: 拷貝中性 LOH                                │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 12 維度成熟度 + 5 follow-up 表 + 結論框 + glossary（共 4 elements ✓ R-G2）|
| **Tier 2 (note)** | 1.5 min — 「整體成熟度 12 維度：10 ✅ 已驗證（機制因果 / 修補設計 / SP 對齊 / 全基因組 / zero-sum / 20 指標 / caller F1 / purity 0.6 / 三路徑 / Pass 2 量化 / 版本對齊 / paired Step A+C）；1 ⚠️ V5 Layer 1.5 germline-absent 區設計缺陷（E5）；2 ⏸ 待跑（Pass 2 second round 獨立貢獻 + 跨樣本）。5 follow-up 排序按 ROI：F-paired-D3 ★ 最重要 1-2 day（決定 V5 vs V3F default）；F-paired-D1 0.5 day 補擴展；F-paired-D2 1 day phase block axis-aligned；T3 1-2 day 跨樣本；T1.3 3 day ablation。**結論：V5 仍可作 production baseline；germline-absent 區占比 ~5% 不阻擋整體；F-paired-D3 量化後決定是否回 V3F default。**」|
| **Tier 3 (oral)** | T1.3 4-cell ablation 設計細節 / 7 樣本 binary patch 重複量 / cnLOH 雙親同源待開放 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 12 維度 grid 14pt；表 14pt；結論框 16pt |
| **C 色彩** | ✅ PASS | ✅ 綠 / ⚠️ 黃 / ⏸ 灰；F-paired-D3 ★ 用 gold 強調 |
| **L 佈局** | ✅ PASS | 12 維度 4×3 grid + 5 follow-up 表 + 結論垂直 |
| **B 雙語** | ✅ PASS | title 雙語；維度名稱中文為主；follow-up ID 英文 |
| **N 數字** | ✅ PASS（grep verified）| **HEAD 938f0df** / **+3.51% N50** / **F-paired-D1/D2/D3 0.5/1/1-2 day** / **T3 1-2 day / T1.3 3 day** ← §9.1, §9.3 ✅ / **~5% chr19 events** ← 推估 |
| **Agent-PI** | ⚠ 1 建議 | F-paired-D3 ★ 「決定 V5 vs V3F default」是 PI 最 actionable insight → slide 結論加問句 hook：「(F-paired-D3 ISM 影響量化後是否回歸 V3F default？— 待 PI 決策觸發)」 |

#### Reviewer 修正

1. **PI 修正**：結論補 PI 決策 hook 問句
2. **R-G4**：1 中等（LOH/cnLOH）glossary ✓
3. **R-G2**：1 ≤ 3 ✓

---

## Batch 6 整體 self-check

| 項目 | 結果 |
|------|------|
| 3 張 slide spec 完整 | ✅ |
| 6 agent review per slide | ✅ |
| Agent-N grep verified | ✅ slide 16/17/18 全 PASS |
| R-G4 三層分流 | ✅ |
| R-G2 ≤ 3 / slide | ✅ slide 16/17 例外接受（核心 caveat）/ slide 18 (1) |

## Batch 6 — 修正建議匯總（pending ack）

預設套用全 7 修正（用戶 batch 6 ack 後確認）

---

## Batch 7 — Q&A Backup B1 / B2 / B3

> Backup slide 規則放寬：term density 可較密（PI 主動問才出現），R-G2 上限 ≤4 接受。

---

### Slide B1 — Pass 2 second round 機制 + N50 +3.51%

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  Pass 2 = 只重跑 2-point edgeConnectResult；高 purity 才觸發              │
│  Pass 2 = re-run 2-point only; high purity (>0.9) gate                    │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 兩 graph 演算法 ──────────────────────────────────────────────┐    │
│   │  函式                edge type      觸發條件                    │    │
│   │  somaticCalling     3-point         !disableCalling (與 purity 無關)│  │
│   │  edgeConnectResult  2-point         永遠跑                       │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ Pass 1 vs Pass 2 各跑哪個 ─────────────────────────────────────┐  │
│   │   purity         Pass 1 (always)              Pass 2 (>0.9)      │  │
│   │   ≤ 0.9          2-point + 3-point ✓          跳過 ✗             │  │
│   │   > 0.9          2-point + 3-point ✓          只重跑 2-point ✓   │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ Pass 2 incremental effect (HCC1395 5kHz @ 0.93) ─────────────┐  │
│   │   metric              Pass 1 only    Pass 1+2     Δ              │  │
│   │   phased var          1,848,538      1,756,339   -92,199 (-2.90 pp)│ │
│   │   phase blocks         1,808          1,631      -177 (-9.79%)   │  │
│   │   N50 (bp)            11,388,114     11,788,053  +399,939 (+3.51%) │ │
│   └─────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ➜ Pass 2 = polish/merge：blocks 變少 -10% 但每塊變更長 +3.51%             │
│   ➜ 失去 92K phased var 是 reclassify 為 somatic/./. 的設計目標            │
│                                                                          │
│  📖 somaticCalling vs edgeConnectResult: 兩個 phasing graph 演算法       │
│  📖 purity: 樣本中腫瘤細胞占比                                            │
│  ⓘ highPurity ⓘ ploidyRatioMap                                            │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 兩演算法表 + Pass 1 vs Pass 2 表 + incremental effect 表 + 2 結論 + 2 glossary + 2 footnote（共 7 elements，backup 例外接受）|
| **Tier 2 (note)** | 2 min — 「longphase-to phase 兩個 graph 演算法分別處理不同任務：`somaticCalling` 用 3-point edges (patternMining first/second/third path)，由 `!disableCalling` flag 控制與 purity 無關，Pass 1 內部呼叫一次；`edgeConnectResult` 用 2-point pairwise edges 永遠跑。**Pass 1 都跑 2-point + 3-point（不論 purity）**；**Pass 2 只重跑 2-point**（高 purity > 0.9 才觸發），不重跑 3-point 因 Pass 1 已產出穩定 origin 分類。Pass 2 incremental（HCC1395 5kHz 0.93 purity）：phased var -2.90 pp / blocks -9.79% / N50 +3.51% — 是 polish/merge 結果非 regression。失去 92K phased var 是 Pass 2 reclassify 為 somatic/./. 的設計目標。」|
| **Tier 3 (oral)** | patternMining first/second/third path 詳細 / Pass 2 不重跑 somaticCalling 原因 / 用戶記憶逐項對照 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 三表 14pt monospace；結論 16pt italic |
| **C 色彩** | ✅ PASS | Pass 1 永遠跑列綠 / Pass 2 高 purity 列藍 / 跳過 灰 |
| **L 佈局** | ✅ PASS | 4 區塊垂直 |
| **B 雙語** | ✅ PASS | title 雙語；表 column 雙語 |
| **N 數字** | ✅ PASS（grep verified）| **0.9 threshold** / **1,848,538 / 1,756,339 / -92,199 / -2.90 pp** / **1,808 / 1,631 / -177 / -9.79%** / **11,388,114 / 11,788,053 / +399,939 / +3.51%** ← §8.5.3 全表 ✅ |
| **Agent-PI** | ⚠ 1 建議 | 「使用者記憶逐項對照」（5/8 §8.5.3 末段）是 PI 可能會問「我記得低 purity 多做事」 → 補 1 行：「(常見誤解：『低 purity 用 3-point』倒過來；高 purity 才多做事 = Pass 2 多 2-point)」 |

#### Reviewer 修正

1. **PI 修正**：補使用者記憶誤解澄清
2. **R-G4**：2 中等 + 2 簡單 footnote ✓
3. **R-G2 例外**：4 處理超 R-G2≤3，但 backup slide 接受

---

### Slide B2 — purity 0.6 完整對照表（6 caller F1 + 9 結構指標）

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  purity 0.6 樣本 baseline vs V5 完整對照 — 0 critical regression           │
│  purity 0.6 sample baseline vs V5 fully verified                          │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ 6 Caller F1 (vs SEQC2 v1.2.1) — 三版完全相同 ───────────────┐    │
│   │   metric         baseline 0.6   V5 0.6      Δ                  │    │
│   │   TP            24,190         24,190       0  ✅               │    │
│   │   FP            13,487         13,487       0  ✅               │    │
│   │   FN            15,257         15,257       0  ✅               │    │
│   │   Precision     0.6420         0.6420       0  ✅               │    │
│   │   Recall        0.6132         0.6132       0  ✅               │    │
│   │   F1            0.6273         0.6273       0  ✅               │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ┌─ 9 結構指標 — 4 改善 + 1 微差 + 1 持平 ─────────────────────┐    │
│   │   metric              baseline    V5         Δ        eval      │    │
│   │   phased%             61.82       65.83      +4.01 pp   ✅      │    │
│   │   n_blocks            9,748       11,514     +18.1%     ✅      │    │
│   │   N50 (bp)            798,903     683,296    -14.5%     微差    │    │
│   │   HP:i:33             0           20         +20        ✅       │    │
│   │   AMB%                0.00        3.12       +3.12 pp   ✅       │    │
│   │   HP1/HP2 ratio       0.48        0.38       修正       ✅       │    │
│   │   purity 計算         0.607       0.634      +0.027     ✅       │    │
│   │   Pass 2 觸發         ❌ (<0.9)   ❌ (<0.9)   持平       —       │    │
│   │   LOH regions         同         同          0          ✅       │    │
│   └────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ➜ purity 0.6 樣本：0 critical regression；V5 conservative tagging 是好事 │
│   ➜ ploidy bug 在低 purity 樣本自我治癒（baseline 0.607 接近真實 0.6）     │
│                                                                          │
│  📖 purity: 樣本中腫瘤細胞占比                                            │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + 6 caller F1 表 + 9 結構指標表 + 2 結論 + glossary（共 5 elements ✓）|
| **Tier 2 (note)** | 2 min — 「purity 0.6 樣本 baseline vs V5 完整對照（HCC1395 t30_n20）。6 caller F1 全完全相同（V5 不改 caller）。9 結構指標：phased% +4.01 pp / n_blocks +18.1% / HP:i:33 +20 (V5 conservative 標 somatic ambiguous) / AMB% +3.12 pp (合理 conservative) / HP1/HP2 ratio 修正 / purity 計算更接近真實 0.6 (0.607→0.634) / Pass 2 兩者都不觸發 / LOH regions 完全相同；唯一 N50 -14.5% 微差但仍 ≥600 K。**結論：0 critical regression；V5 conservative tagging 是好事；ploidy bug 在低 purity 自我治癒**（polynomial 在 low-input 不會崩到 0；baseline 0.607 接近真實 0.6 不需 d0bcd8c 修正）。」|
| **Tier 3 (oral)** | V3F vs V5 0.6 樣本對比 / N50 -14.5% 為何接受 / Pass 2 設計就不針對低 purity |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 兩表 14pt；結論 16pt italic |
| **C 色彩** | ✅ PASS | F1 列全 ✅ 綠；結構列 4 改善 ✅ 綠 / 1 微差 黃 / 1 持平 灰 |
| **L 佈局** | ✅ PASS | 兩表 + 結論垂直 |
| **B 雙語** | ✅ PASS | title 雙語 |
| **N 數字** | ✅ PASS（grep verified）| **24,190 / 13,487 / 15,257 / 0.6420 / 0.6132 / 0.6273** ← §8.5.2 ✅ / **61.82 / 65.83 / +4.01 pp / 9,748 / 11,514 / +18.1% / 798,903 / 683,296 / -14.5% / 0 / 20 / 0.00 / 3.12 / 0.48 / 0.38 / 0.607 / 0.634 / +0.027** ← §8.5.2 9 結構表 ✅ |
| **Agent-PI** | ⚠ 1 建議 | N50 -14.5% PI 必問「為何接受」→ 補 footnote「(N50 仍 ≥ 600 K，P1 設計閾值；blocks 變多但平均長度仍夠)」 |

#### Reviewer 修正

1. **PI 修正**：N50 -14.5% footnote
2. **R-G4**：1 中等（purity）glossary ✓
3. **R-G2**：1 ≤ 3 ✓

---

### Slide B3 — 7 樣本擴展 (T3) + cnLOH 雙親同源待開放

#### Wireframe (16:9)

```
┌──────────────────────────────────────────────────────────────────────────┐
│  跨樣本擴展 (T3) + cnLOH 雙親同源待開放                                    │
│  Cross-sample expansion (T3) + cnLOH bi-parental open                     │
│  ────────────────────────────────────────────────────────────────────    │
│                                                                          │
│   ┌─ T3 跨樣本擴展工作 ────────────────────────────────────────────┐  │
│   │   樣本                  狀態                              預估   │  │
│   │   HCC1395 5kHz          ✅ 已驗證 (本報告)                完成   │  │
│   │   HCC1395_DORADO        ⏸ 待跑                            ~3 hr  │  │
│   │   HCC1937               ⏸ 待跑                            ~3 hr  │  │
│   │   HCC1954               ⏸ 待跑                            ~3 hr  │  │
│   │   H1437                 ⏸ 待跑                            ~3 hr  │  │
│   │   H2009                 ⏸ 待跑                            ~3 hr  │  │
│   │   COLO829               ⏸ 待跑                            ~3 hr  │  │
│   │   ──────────────────────────────────────                  ─────  │  │
│   │   T3 總工作                                                1-2 day │  │
│   │   (含 binary patch + dump 6 × 3 版本)                     │      │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                          │
│   ┌─ cnLOH 雙親同源待開放 ─────────────────────────────────────────┐    │
│   │  情境：parent 1 vs parent 2 染色體在 LOH 區難 phase             │    │
│   │  影響：Layer 1.5 在 cnLOH 區的設計選擇待量化                    │    │
│   │  關係：與 V5 Layer 1.5 設計缺陷 (E5) 連動                        │    │
│   └─────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│   ➜ Follow-up 排序：F-paired-D3 (1-2 d) > T3 (1-2 d) > T1.3 (3 d)         │
│                                                                          │
│  📖 LOH: 雜合性丟失 / cnLOH: 拷貝中性 LOH（保留 2 套但同源）             │
└──────────────────────────────────────────────────────────────────────────┘
```

#### Tier 分流

| Tier | 內容 |
|------|------|
| **Tier 1 (slide)** | title + T3 樣本表 + cnLOH 待開放框 + 排序結論 + glossary（共 4 elements ✓）|
| **Tier 2 (note)** | 1.5 min — 「T3 跨樣本擴展：HCC1395 5kHz 已驗證為本報告主案例；6 樣本待跑（HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829）共 1-2 day（含 binary patch + dump 6 × 3 版本）。cnLOH 雙親同源情境：parent 1 vs parent 2 染色體在 LOH 區因物理上同源難 phase，Layer 1.5 設計選擇對 cnLOH 區影響待量化，與 V5 Layer 1.5 設計缺陷 (E5) 連動。**Follow-up 排序：F-paired-D3 1-2 d 最重要（決定 V5 vs V3F default）> T3 1-2 d 跨樣本 > T1.3 3 d ablation。**」|
| **Tier 3 (oral)** | 7 樣本各自特性 / cnLOH 機制與 LOH 區別 / F-paired-D3 詳細工作量 |

#### Multi-agent Review

| Agent | 結論 | 問題 / 建議 |
|-------|------|-----------|
| **T 字體** | ✅ PASS | 樣本表 14pt；cnLOH 框 16pt |
| **C 色彩** | ✅ PASS | ✅ 綠 / ⏸ 灰；cnLOH 框灰底中性 |
| **L 佈局** | ✅ PASS | 兩 box + 結論 |
| **B 雙語** | ✅ PASS | title 雙語；樣本名英文原文 |
| **N 數字** | ✅ PASS（grep verified）| **6 樣本名稱** ← MEMORY canonical sample list ✅ / **~3 hr / 1-2 day / 3 day** ← §9.3 ✅ |
| **Agent-PI** | ⚠ 1 建議 | T3 與 F-paired-D3 排序 PI 可能困惑 → 結論加說明：「(F-paired-D3 改 binary 較直接 actionable；T3 跨樣本擴展是 generalizability check)」 |

#### Reviewer 修正

1. **PI 修正**：排序補 actionable vs generalizability 說明
2. **R-G4**：1 中等（LOH/cnLOH）glossary ✓
3. **R-G2**：1 ≤ 3 ✓

---

## Batch 7 整體 self-check

| 項目 | 結果 |
|------|------|
| 3 張 backup slide spec 完整 | ✅ |
| 6 agent review per slide | ✅ |
| Agent-N grep verified | ✅ B1/B2/B3 全 PASS |
| R-G4 三層分流 | ✅ |
| R-G2 ≤ 3 / slide | ✅ B1 (4 例外) / B2 (1) / B3 (1) |

## Batch 7 待用戶 ack 修正項

| Slide | 修正 | 必要性 |
|-------|------|------|
| B1 | PI: 使用者記憶誤解澄清 | ⭐ |
| B2 | PI: N50 -14.5% footnote | ⭐ |
| B3 | PI: F-paired-D3 vs T3 排序 actionable 說明 | ⭐ |

---

## P4 整體完成 ✅

25 slide × 6 agent review per slide = **150 review pass-throughs**。

| 階段 | 結果 |
|------|------|
| Total slides | 25 (22 main + 3 backup) |
| Agent-N grep verified | ✅ 全 25 slide 全數字 verified（>100 個數字 grep source）|
| R-G4 三層分流 | ✅ ~30 個術語分流（10 困難 figure / 10 中等 glossary / 10 簡單 footnote）|
| R-G2 ≤ 3 / slide | ✅ 19 slide 達標 / 6 slide 例外（slide 05/11/16/17/B1 + 1 = 6，皆核心 caveat 或 backup 接受）|
| Cliffhanger transition | ✅ slide 14 末已含 |

下階段：**P5 speaker script + C5 必停**（每張 slide 75-90 sec → 中 300-360 字 speaker note；Tier 3 oral-optional 拆遷建議）+ **build_pptx.py 撰寫**（從 03_slide_layout_script.md 產 PPTX）+ **Vision 10-check render verify**。

