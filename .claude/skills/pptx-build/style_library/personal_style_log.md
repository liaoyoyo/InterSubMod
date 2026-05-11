# Personal Style Log — 用戶個人 PPT 風格累積

> 每次用戶在 C3/C4/C5 提出**通用必要**修正（依 `prompts/feedback_classification.md` 分類）累積至此。
> AI 在每張 slide build 前**自動讀取所有 `status: active` 規則**，加進 §20.E 6 問 audit / visual_review 10-check / multi_agent_review Wave 1 specific check。
>
> 啟動日：2026-05-05（pptx-build v2.4 機制建立）

---

## 規則紀錄格式

```markdown
### {YYYY-MM-DD} — {規則簡述}（≤ 30 字）

- **觸發來源**：{slide 出處 / master_draft / 場景}
- **規則細節**：{該做 / 不該做}
- **反例**：{違反此規則的具體 case}
- **與既有規則關係**：{補強 §X / 覆蓋 §Y / 新增}
- **檢核方式**：{哪個 Agent / visual_review 加哪一條}
- **狀態**：active / [PROVISIONAL count=N] / archived
```

---

## Active 規則

### 2026-05-11 — R-G12：字級規範（sans-serif 階層）

- **觸發來源**：Self-Phasing PPT 5/11 封面審查 + WCAG / Texas Tech / Section508 / Six Minutes
- **規則細節**（業界 + 投影規範）：
  - 主標題：36-42pt sans-serif bold
  - 副標題：18-24pt sans-serif regular（非 bold 非斜體）
  - 內文：≥ 18pt（projected presentation ≥ 30pt）
  - footer / glossary：10-12pt（最小限）
  - 不混用多 font family；用 weight (regular/semi-bold/bold) 建 hierarchy
- **檢核方式**：每張 slide build 時 grep font size ≥18pt（除小灰字 footer）
- **狀態**：active
- **Source**: [Section508.gov fonts-typography](https://www.section508.gov/develop/fonts-typography/), [Texas Tech accessibility guide](https://www.ttu.edu/accessibility/digital-accessibility/docs/accessible-powerpoint-guide.html)

---

### 2026-05-11 — R-G11：WCAG 2.1 AA 對比度（2026-04 強制）

- **觸發來源**：Self-Phasing PPT 5/11 web research + WCAG 2.1 AA / WebAIM
- **規則細節**：
  - 一般文字：對比度 ≥ 4.5:1（前後景）
  - 大文字（≥18pt 或 ≥14pt bold）：≥ 3:1
  - 非文字 UI (icons / table borders)：≥ 3:1
  - 不依賴單色傳達（color blindness safe）：用 text label / pattern / shape / icon 雙重編碼
  - Gradient 背景：測試低對比區（PowerPoint 內建 template 多 fail）
- **檢核方式**：WebAIM Contrast Checker / Colour Contrast Analyser；對 slide build 後抽 3-5 處 colour 配對驗證
- **狀態**：active
- **Source**: [WebAIM contrast](https://webaim.org/articles/contrast/), [W3C WCAG 1.4.3](https://www.w3.org/WAI/WCAG21/Understanding/contrast-minimum.html), [DOJ Title II 2024](https://www.ada.gov/) (compliance deadline 2026-04-24)

---

### 2026-05-11 — R-G10：Mayer 多媒體學習 — Coherence + Signaling

- **觸發來源**：Self-Phasing PPT 5/11 + Mayer 12 principles
- **規則細節**：
  - **Coherence**: 排除 extraneous words / pictures / sounds — 封面只放與報告主題直接相關內容（不放裝飾、不放後續細節）
  - **Signaling**: 用 cues 高亮關鍵 — 顏色/粗體/標籤雙重編碼（不單依顏色）
  - **Redundancy**: 圖 + 口述 > 圖 + 口述 + 文字；slide 上不重複 speaker 會講的話
  - **Split-attention**: 文字標籤緊鄰對應 figure（不要圖在左、標籤在右底）
- **檢核方式**：每張 slide build 後 6 問 audit 加「extraneous element 數量」
- **狀態**：active
- **Source**: [Mayer's 12 Principles — DLI](https://www.digitallearninginstitute.com/blog/mayers-principles-multimedia-learning), [PLOS Ten Simple Rules](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009554)

---

### 2026-05-11 — R-G9：受眾語境化（accessible to non-domain）

- **觸發來源**：Self-Phasing PPT 5/11 cover 用戶回饋「V5 / Layer 1.5 不熟者不懂」
- **規則細節**：
  - 主標題對「新生 / 不熟領域 lab member」必須可懂
  - 內部術語（V5 / Layer 1.5 / ClairS-TO / longphase-to / paired audit 等）**不可作主標題**
  - 用「**[topic] 問題修正與驗證**」/「**[topic] 機制與改善**」格式
  - 內部術語可放副標題（但須加 footnote 或 glossary）
- **反例**：「Self-Phasing 整合觀察與 V5 Layer 1.5 設計缺陷」
- **正例**：「Self-Phasing 問題修正與驗證」+ 副「longphase-to tag layer 設計缺陷揭露與三版修正」
- **檢核方式**：每張 slide title 跑 jargon check（內部術語列表 vs title）
- **狀態**：active
- **Source**: [MIT BUMP Talk Tips](https://web.mit.edu/biology/www/undergrad/bump/pdfs/Talk_Tips_BUMP.pdf), [PMC PLOS Lab Meeting Rules](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008953)

---

### 2026-05-11 — R-G8：Speaker Note 不搶後續內容

- **觸發來源**：Self-Phasing PPT 5/11 cover 用戶回饋
- **規則細節**：
  - 封面 speaker note ≤ 中 100 字（30 sec）
  - 內容：報告主題（一句）+ 報告場合 + 報告人簡介
  - **禁止透露**：後續 slide 的具體數字、結論、cliffhanger
  - 各 slide note 不應提前透露「後面 X slide」的內容（漸進揭露 principle）
  - 例外：明確 transition / cliffhanger slide 可在 note 末加「next: ...」提示
- **檢核方式**：speaker note 寫完後 grep 後續 slide 關鍵字（如封面不應含「Phase D」「V6 patch」具體數字）
- **狀態**：active

---

### 2026-05-11 — R-G7：標題非斜體（sans-serif accessibility）

- **觸發來源**：Self-Phasing PPT 5/11 cover 用戶觀察「英文標題斜體不符業界規範」
- **規則細節**：
  - 標題（中文 + EN）統一 sans-serif，**無斜體**
  - 斜體只保留給：figure caption / 引文 / 視覺強調的「emphasized」短語（非標題）
  - Italic 在 sans-serif 字體下 render 不佳（業界共識）
  - a11y 考量：dyslexia / 視覺障礙者對斜體解析困難
- **反例**：cover 內 EN subtitle 用 italic
- **正例**：EN subtitle 用 regular weight，非斜體
- **檢核方式**：grep 「italic」「`<em>`」「`<i>`」於 title / subtitle 區域
- **狀態**：active
- **Source**: [Harvard Catalyst slides](https://catalyst.harvard.edu/writing-communication-center/visualize-science/slides/), [Six Minutes Slide Fonts](https://sixminutes.dlugan.com/slide-fonts/)

---

### 2026-05-11 — R-G6：標題語境化（topic + 動作）

- **觸發來源**：Self-Phasing PPT 5/11 cover 用戶建議「Self-Phasing 問題修正與驗證」
- **規則細節**：
  - 主標題格式：「**[topic] [動作]**」（如「問題修正與驗證」「機制與改善」「驗證與應用」）
  - 不用問句、不用 jargon、不用版本號（V5/V6）作主標題
  - ≤ 18 字（中文）
- **範例庫**：
  - ✅ 「Self-Phasing 問題修正與驗證」
  - ✅ 「LOH × AF × CN 三維特徵分層」
  - ✅ 「ISM 跨樣本一致性驗證」
  - ❌ 「Self-Phasing 整合觀察與 V5 Layer 1.5 設計缺陷」（含版本 jargon）
  - ❌ 「為什麼 self-phasing 重要？」（問句）
- **狀態**：active

---

### 2026-05-11 — R-G5：封面內容規範（必/可選/禁止三類）

- **觸發來源**：Self-Phasing PPT 5/11 cover 用戶回饋 + PMC PLOS / MIT / Harvard 業界規範
- **規則細節**：

**必要元素**：
1. 主標題（單一中文，無中英並列，sentence-style，≤ 18 字）
2. 副標題（中文短句，無斜體，解釋語境）
3. 日期（YYYY-MM-DD 或日期區間）
4. 場合（「實驗室週報」/「lab meeting」/「PI 報告」/「Defense」）
5. 報告人姓名

**可選元素**：
- 機構 / 研究室 / 系所 logo（如有）
- 共同貢獻者（如有）
- 短 footer source 引用（10-12pt 灰字）

**禁止元素**：
- ❌ 中英並列大標題（破壞 visual hierarchy）
- ❌ 技術細節字串（commit hash / version dates 序列）— 放後面目錄
- ❌ 過長副標題（> 30 字）
- ❌ 斜體標題（特別 sans-serif 字體）
- ❌ 內部術語 jargon 作主標題

- **檢核方式**：每次 cover slide build 時跑「5 必要 + 0 禁止」checklist
- **狀態**：active
- **Source**: [PLOS Ten Simple Rules effective slides](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009554), [MIT Comm Lab slideshow](https://mitcommlab.mit.edu/eecs/commkit/slideshow/), [InkPPT title slide examples 2026](https://www.inkppt.com/post/powerpoint-title-slide-examples)

---

### 2026-05-10 — R-G4：未解釋術語三層分流（footnote / glossary box / 獨立補充頁）

- **觸發來源**：Self-Phasing 整合 PPT C3 確認，用戶提出「不好了解或第一次出現的名詞，理論上都需要解釋；如果不在知識庫裡，必須解釋；困難概念需補充空間或頁面、盡量以圖示或理解性說明；小段名詞敘述或概念定義，可在頁面小區塊加類似書本名詞定義框 + 文字解釋」
- **規則細節**（三層分流，AI 在每張 slide P3 規格時自動分類每個術語）：
  1. **困難概念**（多步因果 / 抽象結構 / 跨層 mechanism）→ **獨立 explanatory slide** 或 **大圖示**；以視覺化 + 理解性說明為主
     - 例：self-phasing 球員兼裁判機制 / Pass 1 vs Pass 2 因果鏈 / ploidy bug → purity=0 鏈 / Layer 1.5 fallback 邏輯
  2. **中等術語**（短語可定義但需上下文 / 包含內部命名 / 跨領域）→ **in-slide glossary box**（書本式名詞定義框）置於 slide 角落或 sidebar
     - 例：PoN / germline het / somatic mutation / phasing graph / HP:i: vs HP:Z: / sub-clone / haplotype / LOH / cnLOH
  3. **簡單術語**（單詞可一句解釋 / 第一次出現但不困難）→ **footnote ≤30 字**
     - 例：som_ratio / phased votes / INDEL guard / two-layer / threshold / break early
- **判定基線**（與 R-G2 互補）：
  - 「不在知識庫裡」→ 必解釋（用 `knowledge_search` 查 `/big8_disk/liaoyoyo2001/Knowledge/` 是否有條目；無 → 必解釋）
  - 「第一次出現」→ 必解釋（不重複；後續 slide 引用同術語不再解）
  - PI 領域熟悉度 → 對齊母稿 frontmatter `audience` field
- **反例**：
  - 整張 slide 充滿縮寫（PoN / cnLOH / LOH / GT / HP:i: / sub-clone）但無任何解釋 → PI 看不懂
  - 把「球員兼裁判」隱喻當理所當然不畫示意圖 → PI 不熟 phasing graph 無法理解
  - 對「ploidy bug → purity=0」單行帶過 → PI 無法判斷修補必要性
- **與既有規則關係**：
  - 補強 R-G2「PI 不熟術語 ≤ 3 / slide」— R-G2 是「上限約束」，R-G4 是「處理方式分層」
  - 與 myPPT playbook §6 整合：footnote / glossary box / explanatory slide 各佔 Tier 1 element 預算的方式不同
- **檢核方式**：
  - P3 規格時：每張 slide 列「困難 / 中等 / 簡單」術語三層清單，並對應 footnote / glossary / explanatory slide 處理
  - P4 Wave 1 加 Agent-G（Glossary verification）並行：對每張 slide 確認所有術語都已分流處理（不漏 / 不過度解釋既已知術語）
  - feedback_classification：未解釋術語 → 通用必要修正
- **狀態**：active

---

### 2026-05-07 — R-G2：每張 slide ≤ 3 個 PI 不熟術語

- **觸發來源**：v2.1 試用第 2 次 PPTX 受眾理解度評估 — S04/S06/S13/S14 被識別為「技術詞密度過高需重設計」（4/18 = 22% 失敗率）
- **規則細節**：
  - 受眾為「相鄰專長」（領域專家但不熟程式碼細節）時，每張 slide on-slide 文字中 PI 不熟悉的術語（內部 thread 命名 / git 術語 / SAM tag / 程式碼層級概念）≤ 3 個
  - 多於 3 個必補 1 行 footnote 解釋（≤ 30 字）
  - 「PI 不熟」判定基線：knowledge base `/big8_disk/liaoyoyo2001/Knowledge/README.md` 與 PI 報告 `audience: PI（領域專家但未親自操作程式碼）` 對齊
- **反例**：
  - S06 同時用「somaticCalling 重跑」「只重 phase 不重 call」「path 1/2/3」「second round」「highPurity threshold」（6 個術語）→ PI 30 sec 內難以追上
  - S13 用「marker / phasing signature / Fisher odds / NG≥3」（4 個術語）無解釋 → PI 必問「什麼是 phasing signature」
- **與既有規則關係**：
  - 補強 myPPT playbook §6 「中文 ≤ 60 字 / 英文 ≤ 30 word」，加術語密度限制
  - 對齊 PLOS Rule 4「Essential Only」與 Stanford Online「avoid jargon」原則
- **檢核方式**：
  - §20.E 6 問 audit 第 8 問：「此 slide 上的內部術語 ≤ 3 個？多於 3 個是否補 footnote？」
  - Wave 1 Agent-B（雙語）擴展職責：標記未解釋的內部術語清單
- **狀態**：active

---

### 2026-05-07 — R-G3：Wave 1 多 Agent review 加 Agent-N（Number verification）並行

- **觸發來源**：v2.1 試用第 2 次 Wave 1 (T/C/L/B 4 並行) 全部通過但仍漏掉 3 處事實錯誤；Wave 2 Agent-D 才事後抓到，但 Agent-D 也誤判 OLD baseline = 0（實為 2,640）
- **規則細節**：
  - Wave 1 從 4 並行（T 字體 / C 色彩 / L 佈局 / B 雙語）擴展為 5 並行，新增 **Agent-N**：
    - 對每張 slide 每個數字 grep source 文件確認 metric 名稱、scope、變化方向、精度
    - 產出 verification table（Slide / Number / Source path+line / Match? / Note）
    - 「[U]」「[O]」「[F]」標記必確認（[U] 不應冒用為 [F]）
  - Agent-N 與 T/C/L/B 並行（不延遲整體 review wall time）
  - 與 Wave 2 Agent-D 分工：Agent-N = per-slide 即時驗證；Agent-D = 整份 caveat completeness + over-claim 整合
- **反例**：
  - v2.1 第 2 次：S04 HP:i:33「14524 → 0」3 個事實錯誤同時存在（誤拿 metric / 誤值 / 誤方向），Wave 1 4 Agent 全部通過僅做視覺檢查
  - Agent-D 自己也誤判 OLD baseline = 0（PPTX 沒抓到，但 source 真實值 2,640）
- **與既有規則關係**：
  - 強化 R-G1（metric 命名 single-source verification）— R-G1 是「寫者責任」，R-G3 是「驗者責任」雙保險
  - 修改 `InterSubMod/.claude/skills/pptx-build/prompts/multi_agent_review.md` 加 Agent-N 段落
- **檢核方式**：每張 slide build 後 Wave 1 5 Agent 並行；Agent-N 產 PASS/FAIL/PARTIAL × 4 項（metric 名稱 / scope / 方向 / 精度）
- **狀態**：active

---

### 2026-05-07 — Metric 引用 single-source verification + 必明示 scope

- **觸發來源**：v2.1 試用第 2 次 PPTX `output.pptx` S04/S07/S08；3 處事實錯誤涉及同一個 metric (`HP:i:33` / `HP_33` / `HP33`) 在不同 audit 文件不同 BAM scope 下被混用引用
- **規則細節**：
  - 任何進入 PPTX 的數字必須 grep source 文件確認 (a) metric 名稱原文、(b) BAM/sample scope（whole genome / per-site / cherry-pick / multi-sample 等）、(c) 變化方向、(d) 數字精度
  - 同一張 slide 出現的多筆數字若 scope 不同，必加 column 或 footer 明示
  - 所謂「等價」「相同」「PASS」結論只能在同 scope 數字之間下；跨 scope 比較只能談「方向」不能談「等價量化」
- **反例**：
  - PPTX 把 OLD V5 全 BAM HP_33=14,524 跟 noPath3 cherry-picked HP_33=28 並列為「14,524」並寫「等價性」（事實上 28 vs 14,524 是不同 scope）
  - PPTX 寫「HP:i:33 從 14524 → 0 (-100%)」混合三個錯誤：metric 取錯來源、數字取錯、方向錯
- **與既有規則關係**：
  - 強化 Memory `feedback_feature_name_vs_definition_rule`（C++ feature 必查原始碼），擴展到「PPTX 數字引用必查 source 文件」
  - 補強 weekly-report W4 邏輯紅旗（原本只查過度宣稱 / 流水帳）
  - 補強 pptx-build §20.E 6 問 audit（新增「每個數字是否 single-source verified」第 7 問）
- **檢核方式**：
  - Wave 1 新增 Agent-N（Number verification）並行：對每張 slide 每個數字 grep source 文件，產 verification table
  - §20.E 6 問 audit 第 7 問：「此 slide 每個數字是否 grep source 文件確認過？」
  - feedback_classification 觸發時若分類「scope label 缺失」自動標 metric-rule 系列
- **狀態**：active



---

## [PROVISIONAL] 觀察中規則

（尚無 — [PROVISIONAL] 累積 ≥3 次同類修正後升級到 active）

---

## Archived 規則（歷史 / 已棄用）

（尚無）

---

## 累積統計

| 指標 | 數值 |
|------|------|
| Active 規則數 | 12 |
| [PROVISIONAL] 規則數 | 0 |
| Archived 規則數 | 0 |
| 最近更新 | 2026-05-11 |

## 規則索引

| 規則 | 主題 | Source |
|------|------|--------|
| Metric scope | 數字 single-source verification | 內部 v2.1 試用 |
| R-G2 | PI 不熟術語 ≤ 3 / slide | 內部 v2.1 試用 |
| R-G3 | Wave 1 加 Agent-N 並行 | 內部 v2.1 試用 |
| R-G4 | 未解釋術語三層分流 | 自反饋 |
| R-G5 | 封面內容規範（必/可選/禁止三類）| PLOS / MIT / Harvard |
| R-G6 | 標題語境化（topic + 動作）| 自反饋 + PLOS |
| R-G7 | 標題非斜體（sans-serif a11y）| Harvard / Six Minutes |
| R-G8 | Speaker Note 不搶後續內容 | 自反饋 |
| R-G9 | 受眾語境化（non-domain accessible）| MIT BUMP / PLOS |
| R-G10 | Mayer 多媒體（Coherence + Signaling）| Mayer 12 principles |
| R-G11 | WCAG 2.1 AA 對比度 | W3C / WebAIM / DOJ Title II 2024 |
| R-G12 | 字級規範（sans-serif 階層）| Section508 / Texas Tech |

## 與既有規則關係（衝突 / 覆蓋 / 補強）

當用戶通用規則與 playbook 既有規則（§5/§6/§14/§16/§20 等）衝突時：

| 既有規則 | 用戶覆蓋 | 處理 |
|---------|---------|------|
| §6 中文 ≤ 60 字 | 用戶要 ≤ 50 字 | 標 [OVERRIDE]，套用更嚴標準 |
| §14 紅色僅 caveat | 用戶要紅色強調 | 標 [OVERRIDE]，但提示違反 colorblind safe |
| ... | ... | ... |

[OVERRIDE] 標記讓 AI 知道是用戶有意識的覆蓋，不是 bug。

## 後續檢核整合（自動）

每張 slide build 前 AI 執行：

```python
# pseudo-code
def pre_build_check(slide_n):
    personal_rules = load_active_rules('style_library/personal_style_log.md')

    # 加進 §20.E 6 問補充
    extended_audit = base_6_questions + [
        f"是否符合用戶規則 '{rule['name']}'？" for rule in personal_rules
    ]

    # 加進 visual_review 10-check 擴充
    extended_checks = base_10_checks + [
        rule['check_method'] for rule in personal_rules
    ]

    # 加進 multi_agent_review Wave 1
    for agent in ['T', 'C', 'L', 'B']:
        agent_specific = filter_by_agent(personal_rules, agent)
        extend_agent_prompt(agent, agent_specific)

    return run_all_checks(...)
```

## 規則生命週期

```
新規則
   ↓ feedback_classification 分類 = 通用
Active (status: active)
   ↓ 用戶長期未提及 + AI 無觸發
經過 N 次 PPT 都未引用
   ↓ memory-consolidation 建議 archive
Archived (歷史保留，不主動套用)
   ↓ 用戶手動恢復
Active 再次套用
```

## 範例：未來規則寫法（示意，2026-05-05 尚無真實案例）

```markdown
### 2026-XX-XX — 中文標題禁用問號

- **觸發來源**：v3 24-slide audit 中 S11 標題「為何 self-phasing 會卡 ISM？」
- **規則細節**：assertion-evidence 標題用陳述句，不用疑問句（「Self-phasing 卡住 ISM 5 目標」優於問號）
- **反例**：「為何 X 會 Y？」「How does X work?」
- **與既有規則關係**：補強 §5 第 1 條 heading=thesis sentence
- **檢核方式**：Agent-T 字體 audit 加「標題無 ?」項
- **狀態**：active
```

## 與 memory 系統互動

InterSubMod 已有 user-level Memory（`/bip7_disk/.../memory/`），用 `feedback_*` 類型已紀錄 PPT 偏好（如 `feedback_pptx_director_storyboard`、`feedback_pptx_visual_first_philosophy`）。

**分工**：
- **Memory `/feedback_*.md`** — 高層哲學偏好（如「視覺優先」「導演審查」）
- **本檔 personal_style_log.md** — 細粒度可檢查規則（字數 / 字級 / 顏色 / 對齊）

兩者互補：Memory 是「為什麼」，本檔是「具體怎麼做」。
