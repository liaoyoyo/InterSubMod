# weekly-report skill v2 — playbook (主方法論)

> 16 章節主方法論文件。SKILL.md 是觸發殼，本檔是細則。
> 各 § 對應 SKILL.md 的 W1-W7 階段、5 個 checkpoint 互動模板、4 主線類型、4 層分類規則、handoff。

## §1 系統理念

**核心原則：先有母稿，才有簡報。**

過去 v1 weekly-report 直接從 raw data 跳到 PPTX 產生，痛點：
1. PPTX 投影片承擔過多論述責任，文字密度過高
2. 沒有「資料佐證的論述母稿」作為 source of truth，每次重做 PPT 從零開始
3. 講稿與 slide 邊界模糊，重要 caveat 遺漏
4. 教授追問常常 surprise，因事先沒預演

v2 升級理念：
- 母稿 = 結構化文字輸出（Layer 0-4 + 17 段標籤），是 PI/教授可讀的中間產物
- PPT = 母稿的視覺呈現層，由 pptx-build sub-skill 處理
- 三角色分工：
  - **AI 整理者**：raw data 收集、自動分類、紅旗掃描、追問預測
  - **AI 追問者**：透過每輪 ≤ 5 個 AskUserQuestion 確認用戶意圖
  - **用戶判斷者**：4 主線選擇、分類修正、母稿批准（C1/C4 必停）

不應該：
- AI 不替用戶決定研究立場
- AI 不把 [I] 推論寫成 [F] 事實
- AI 不過度美化或誇大
- AI 不一次塞太多問題給用戶

## §2 W1 Raw Data 收集（→ C0）

### 2.1 預期輸入來源

| 來源 | 命令範例 | 過濾條件 |
|------|---------|---------|
| Git log | `git log --since="7 days ago" --pretty=format:"%h %s (%ad)" --date=short` | 排除 `chore:` `docs:` 自動化雜項 |
| experiments/in_progress 新增 | `find docs/experiments/in_progress -newer .claude/last_weekly_report -name "*.md"` | 取本週修改 |
| evidence_ledger 新增 | `tail -50 research/autoresearch/evidence_ledger.jsonl \| jq 'select(.date >= "...")'` | 抽 verdict + hypothesis_id |
| CURRENT_FOCUS 變動 | `git log --follow -p docs/CURRENT_FOCUS.md --since="7 days ago"` | 看主軸有無切換 |
| 上週週報 Layer 3-4 | `ls -t docs/reports/validated/YYYY/MM/ \| head -2` | 萃取「仍然未知」「下週優先行動」 |

### 2.2 結構化清單格式（W1 輸出）

```markdown
### W1 Raw Data 結構化清單

**時間軸**：
- 2026-04-29 a7cb7e0 feat(skills): weekly-report v2 升級
- 2026-04-30 ... (列至少 5 個 commit)

**產出 artifacts**：
- InterSubMod/docs/experiments/in_progress/2026/04/<file>.md
- InterSubMod/.claude/skills/weekly-report/SKILL.md (升級為 v2)

**evidence ledger 新增**：
- [ID-2026-04-30-001] HPFineNGroups Phase 2B PARTIAL POSITIVE
- ...

**CURRENT_FOCUS 變動**：
- 主軸從「Phase 1A F1 優化」切到「Phase 2 Normal Methylation」

**用戶補述**（C0 後填）：
- ...
```

### 2.3 C0 互動模板 → `prompts/raw_data_collect.md`

## §3 W2 主線類型識別（→ C1，fast-track 必停）

### 3.1 4 主線類型

| 類型 | 觸發條件 | 對應 template |
|------|---------|--------------|
| progress（進展） | commit 多 + 量化 delta + [F] 標籤多 | `templates/progress_focus.md` |
| problem（問題） | commit 少 + experiments 標 [BLOCKER]/[ANOMALY] | `templates/problem_focus.md` |
| advisor_consult（求協助） | 多 candidate path 都有部分 evidence | `templates/advisor_consult.md` |
| new_direction_explore（探索） | pilot N≤3 + experiments 標 [PILOT] | `templates/new_direction_explore.md` |

### 3.2 混合場景處理

如進展 + 問題混合：以教授最關心的為主，其他降為 sub-thread。
範例：本週主要 deliverable 是進展型（5/7 通過），但有一個關鍵 caveat（HCC1395 反向）需特別處理 → 主線=progress，sub-thread=problem 在 §5 [U] 標出。

### 3.3 main statement（≤ 30 字）撰寫規則

- 必含具體數字或具體動詞
- 必能對應 §1 母稿主線
- 範例好：「7 樣本 paired-pure ΔF1=+0.0112 重新確認」
- 範例不好：「本週做了一些研究」（無動詞、無數字）

### 3.4 C1 互動模板 → `prompts/main_thread_identify.md`

## §4 W3 內容 4 層分類 [F]/[O]/[I]/[U]（→ C2）

### 4.1 4 層定義

詳見 `references/LAYER_STRUCTURE.md` §C。

| 標籤 | 必要欄位 | 描述語氣 |
|:-:|---------|---------|
| [F] Fact | source_path + validation_status | 「已驗證」「確認為」 |
| [O] Observation | sample_count + threshold_to_promote | 「初步觀察到」 |
| [I] Inference | basis_evidence + alternatives | 「推測」「可能」 |
| [U] Unconfirmed | what_to_check + when | 「待釐清」「需要 X」 |

### 4.2 標籤決策樹

```
有具體 source？
├─ 否 → 看是「推測」還是「待釐清的疑問」
│       ├─ 推測 → [I]
│       └─ 疑問 → [U]
└─ 是 → 達 validation threshold？
        ├─ 否 → [O]
        └─ 是 → [F]
```

### 4.3 4 層 vs Tier 1/2/3 並用

不同維度：
- Tier 1/2/3 = 重要性（決定 Layer 2 呈現深度）
- [F]/[O]/[I]/[U] = 真實性（決定描述語氣）

範例：Tier 1 + [O] 表「最重要的一個觀察，但 N 不足」。

### 4.4 過度宣稱紅旗（W5 預先掃描）

- [F] 描述含「初步」「可能」 → ⚠ 應降 [O] 或 [I]
- [O] 描述含「證實」「確認」 → ⚠ 改觀察語氣
- [F] 但 sample_count < threshold → ⚠ 降 [O]

## §5 W4 重點排序 + 4 桶分流（→ C2）

### 5.1 5 維度評分表

| 維度 | 1 分 | 3 分 | 5 分 |
|------|------|------|------|
| 研究重要性 | 邊緣 | 中等 | main thesis 一部分 |
| 證據強度 | [U] | [I] | [F] |
| 教授關心程度 | 不太關心 | 中等 | 必問 |
| 影響下週計畫 | 無關 | 弱關聯 | 直接決定 |
| 適合簡報呈現 | 純文字 | 部分視覺 | 強視覺 |

### 5.2 桶分流（加總分數）

| 加總 | 桶 | 上限 | Tier |
|------|---|------|------|
| 18-25 | PPT slide | ≤ 8 | Tier 1 |
| 13-17 | 講稿 speaker note | ≤ 15 | Tier 2 |
| 8-12 | 備註 oral-optional | 不限 | Tier 3 |
| <8 | 暫存 | 不限 | (棄) |

### 5.3 上限超出處理

PPT 桶 > 8 時，5 維度最弱者降到講稿桶。若用戶堅持保留標 [ESCALATION]。

### 5.4 C2 互動模板（合併 W3+W4）→ `prompts/classify_and_sort.md`

## §6 W5 邏輯紅旗檢查（→ C3）

### 6.1 過度宣稱紅旗（6 條）

詳見 `references/LAYER_STRUCTURE.md` §C 末段。

### 6.2 流水帳紅旗（4 條）

- 同層級 ≥ 5 個 bullet 平列 → 重組
- PPT 桶 > 5 筆無因果連接詞 → 補因果鏈
- 每件事獨立段落 → ≥ 2 件合併
- 「本週做了 ABCDEFG」式列表 → 重排優先序

### 6.3 教授視角缺紅旗（3 條）

- 母稿無 §17 教授追問段 → 強制補
- §16 下週計畫不接續本週發現 → 紅旗
- 求協助型主線無 §5 候選 trade-off → 紅旗

### 6.4 AI 自動掃描方式

keyword 比對 + 結構檢查（如：母稿有無 §17 H2 標題）

## §7 W6 教授問答預測（→ C3）

### 7.1 預測 5-7 個追問依據

- 主線類型 → 對應典型問題（templates/{type}_focus.md「教授可能追問」段）
- 4 層分類中的 [I] / [U] → 教授會追問
- 過度宣稱紅旗未完全修正 → 教授會質疑
- 與過去報告的不一致 → 教授會記得

### 7.2 預備回答骨架

每個追問必須有 1 段預備回答，含：
- 結論（一句）
- evidence 引用（檔案 path / 數字）
- alternative 解釋（如有）

### 7.3 C3 互動模板（合併 W5+W6）→ `prompts/check_and_predict.md`

## §8 W7 母稿產出（→ C4，fast-track 必停）

### 8.1 17 段 → Layer 0-4 mapping

詳見 `references/LAYER_STRUCTURE.md` §B。

母稿主骨架：

```
# §1 主線（≤ 30 字）+ §2 一句話重點   ← Layer 0.1
## Layer 0.2 背景知識                   ← §13
## Layer 0.3 上週前情提要
## Layer 1 已建立知識參考
## Layer 2 本週調查 (Thread × N)
   ### 證據卡 Tier 1 → §3 [F]
   ### 證據卡 Tier 2 → §4 [O][I]
   ### 證據卡 Tier 3 → §5 [U]
## §7-§8 排序 + 報告順序                 ← Layer 2 整合
## Layer 3 整合更新                     ← §11 §12 §13 §14
## Layer 4 未來方向                     ← §16 §17
## 附錄: §6 不放 PPT + §15 暫存 + §9 §10 PPT 模板建議
```

### 8.2 frontmatter schema

詳見 `references/HANDOFF_TO_PPTX_BUILD.md` §2。

### 8.3 母稿輸出位置

`InterSubMod/docs/reports/validated/YYYY/MM/YYYYMMDD_週報_<主題>/master_draft.md`

### 8.4 C4 互動模板 → `prompts/master_draft_finalize.md`

## §9 Handoff 4 選

### 9.1 4 個選項

詳見 SKILL.md handoff 段或 `prompts/handoff_to_pptx_build.md`。

A: 立即觸發 pptx-build (--from-draft)
B: 母稿留檔，下次手動
C: 母稿即終點
D: 母稿留檔 + 加寫下週計畫

### 9.2 next_week_plan.md 自動萃取規則

選 D 時，weekly-report 自動萃取下列段到 `next_week_plan.md`：
- 母稿 §16 下一步行動 → 直接複製
- 母稿 §17 教授追問 → 轉成「下週要準備哪些 evidence」
- 母稿 §11 補資料 → 轉成「下週要產的 artifacts」
- 母稿 §15 暫存 → 標 [SHELVED] 等下週判斷
- 母稿 §5 [U] → 轉成「下週 blocker 追蹤」

## §10 範例案例

詳見 `examples/master_draft_example.md`：HPFineNGroups Phase 2B 重驗證（progress 主 + problem sub）。

## §11 反例（不該這樣寫）

### 反例 1：[I] 寫成 [F]
❌ 「驗證 LOH-constrained phasing 為新方向」（這是推論，無 7 樣本驗證）
✅ 「推測 LOH-constrained phasing 可能是新方向 [I]，依據 NG=2 same-hap 6/6 觀察」

### 反例 2：流水帳 8 件
❌ 「本週做了 A、B、C、D、E、F、G、H」
✅ 重組為 3 個 narrative 群組 + 因果連接詞

### 反例 3：母稿無 §17
❌ 跳過教授追問預測（C3 不可省）
✅ 母稿必含 §17，5-7 個追問 + 預備回答

### 反例 4：主線混雜
❌ 同時宣告「本週進展 + 本週問題 + 本週新方向」3 主線
✅ 取主敘事為主軸，其他降為 sub-thread

### 反例 5：raw data 收集遺漏
❌ AI W1 只看 git log，漏掉 evidence_ledger 新增
✅ W1 全部 5 來源都掃過

## §12 互動 protocol 模板（W1-W7 對應 prompts/）

| W → C | prompts/ 檔案 | 必停 (fast-track) |
|-------|------------|------------------|
| W1 → C0 | raw_data_collect.md | ❌ 非必停 |
| W2 → C1 | main_thread_identify.md | ✅ **必停** |
| W3+W4 → C2 | classify_and_sort.md | ❌ 非必停 |
| W5+W6 → C3 | check_and_predict.md | ❌ 非必停 |
| W7 → C4 | master_draft_finalize.md | ✅ **必停** |
| C4 後 handoff | handoff_to_pptx_build.md | ✅ **必停** |

每個 prompts/ 檔案 AskUserQuestion **≤ 5 個 question**（避免一次轟炸用戶）。

## §13 fast-track 模式

### 13.1 觸發條件

用戶明示「全自動」「auto」。

### 13.2 跳過項目

- C0 raw data: AI 預設清單，30 秒倒數可中斷
- C2 分類+排序: AI 預設標籤+評分，30 秒倒數
- C3 邏輯+追問: AI 自動修正紅旗，30 秒倒數

### 13.3 仍強制項目

- **C1 主線類型 + main statement**：必停（決定整份母稿走向）
- **C4 母稿確認**：必停（最終輸出，不能未確認就 handoff）
- handoff 4 選：必停（不能自動觸發 pptx-build）

### 13.4 fast-track 報告

完成後產「自動執行報告」，列 AI 所有自動決策（用戶可事後 audit）。

## §14 與既有 skill 互引用

| Skill | 關係 | 使用場景 |
|-------|------|---------|
| pptx-build | **下游接棒** | C4 後選 A 觸發 |
| myPPT | **上游觸發**（場景識別後委派）| myPPT 識別週報情境 → 觸發本 skill |
| review-evidence | **上游工具** | W1 raw data 收集，掃 evidence_ledger |
| provenance-tier-audit | **上游工具** | W1 sanity check tier 分佈 |
| structured-tech-report | **平行不重疊** | 單一工程改動 deep dive，非週期性 |
| confirmation-protocol | **規範來源** | C0-C4 對應 Hard Gate / Gate / Review 級別 |
| doc-standards | **規範來源** | 母稿命名規範 |
| feature-layered-observation | **互補** | W1 收集已用 Step 1-6 產出的 deep-dive 報告 |

## §15 輸出檔案管理

### 15.1 路徑規範

```
InterSubMod/docs/reports/validated/YYYY/MM/YYYYMMDD_週報_<主題>/
├── master_draft.md                  # 第 1 版
├── master_draft_v2.md (optional)    # 用戶要求修改時
├── next_week_plan.md (D 選項)
└── pptx_build/ (A 選項後 pptx-build 產出)
```

### 15.2 與舊版兼容性

舊 v1 週報路徑同樣為 `docs/reports/validated/YYYY/MM/`，命名為 `YYYYMMDD_研究週報_<period>_<主題>_NN.md`。
v2 沿用此路徑，但檔名統一為 `master_draft.md`（在以週為命名的子資料夾內）。

### 15.3 INDEX 更新

C4 完成後，自動建議用戶在 `docs/reports/INDEX.md` 加新項。

## §16 success criteria（驗收標準）

| 指標 | 通過標準 |
|------|---------|
| 母稿可讀性 | PI 不需翻舊報告即可獨立讀懂 |
| 教授追問預測準確度 | 5-7 個追問覆蓋 ≥ 80% 實際追問（事後 audit）|
| [F]/[O]/[I]/[U] 標籤準確 | 用戶 audit 通過率 ≥ 90% |
| handoff 平順 | 用戶不需在 pptx-build P1 重填 main thesis |
| 紅旗檢出率 | W5 自動掃描 ≥ 80% 過度宣稱用法 |
| 4 桶上限遵守 | PPT ≤ 8 / 講稿 ≤ 15 |
| fast-track 模式可用 | C1/C4 仍互動，其他自動推進 < 60 秒 |
| 母稿字數 | 2000-4000 字（17 段平均 100-250 字/段）|

落地後第一個試用案例：下次實際週報走完整 W1-W7 + handoff 4 選，事後對照 `examples/master_draft_example.md` 結構是否齊全。
