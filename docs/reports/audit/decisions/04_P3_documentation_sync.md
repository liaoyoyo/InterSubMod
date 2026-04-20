# 04 — P3 Documentation Sync Decisions（5 題）

> **建立日期**: 2026-04-19

> **優先等級**：P3（文件同步小任務）
>
> **依賴**：全部獨立，可多線並行；預計每題 < 1 小時

---

## 問題 P3-A: HPFineNGroups filter 版本更新

**問題描述**：
`docs/CURRENT_FOCUS.md:98-103` R4 章節仍記載舊 HPFineNGroups filter（NG≥4+NR≥80），未更新為 2026-04-18 F pilot 確認的新 filter（**NG=4+AF<0.4+NR≥80**）。

**影響範圍**：
- 受影響文件：`CURRENT_FOCUS.md:98-103`, 可能其他 research_landscape 文件
- 受影響結論：C16 的 filter 定義引用
- 若不處理的風險：執行者接手時引用舊 filter，重跑結果不一致

**現況證據**：
- Audit card: `cards/C16_HPFineNGroups.md`（已用新 filter）
- 原始 pilot: `research/F_pilot/20260418_*`

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 修 `CURRENT_FOCUS.md:98-103` + grep 全 docs 找其他舊 filter 引用 | 20 分鐘 | 低 |
| B | 僅修 CURRENT_FOCUS，audit card 已是新版無需動 | 5 分鐘 | 中 — 其他文件可能殘留 |
| C | 不改，記錄於 R4 之 TODO 下 | 0 分鐘 | 高 — 執行時誤用 |

**驗證標準**：
- `grep -rn "NG.*NR.*80" docs/ | grep -v "AF.*0\.4"` 預期僅回傳歷史報告（有日期註腳）
- `CURRENT_FOCUS.md:98-103` 顯示新 filter

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P3-B: 22 結論納入 06_結論穩定性審查

**問題描述**：
`docs/reports/research_landscape/06_結論穩定性審查.md` 穩定度評分系統建立於 2026-04-04（14 結論時點）。之後新增結論 17-22（特別是 19, 20, 21 為補充結論）未完整納入主表。00_INDEX 總表與 06_審查表格不一致。

**影響範圍**：
- 受影響文件：`06_結論穩定性審查.md`
- 受影響結論：C17, C18, C19, C20, C21, C22（格局一致性）
- 若不處理的風險：新人看 06_審查會以為只有 14 結論

**現況證據**：
- Phase 1 Agent A 觀察 I4
- 對照：`00_INDEX.md` 列 22 結論 vs `06_審查.md` 主表

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 在 `06_審查.md` 新增「2026-04-18 補充結論 17-22」章節，5-6 結論各一段 | 1 小時 | 低 |
| B | 重寫 06_審查.md 整合所有 22 結論為統一主表 | 半天 | 中 — 結構大幅改動 |
| C | 不改 06_審查.md，引用時指向 audit card | 0 分鐘 | 中 — 06_審查 地位下降 |

**驗證標準**：
- `06_審查.md` 有 17-22 結論段落
- 每結論至少含穩定度、狀態、核心發現

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P3-C: TODO 編號統一

**問題描述**：
`CURRENT_FOCUS.md` P0-3（已完成）vs `06_結論穩定性審查.md` P0-1（待驗證）使用不同編號系統指同一事項（LOH.bed 生成機制）。編號衝突導致執行者誤判優先序。

**影響範圍**：
- 受影響文件：`CURRENT_FOCUS.md`, `06_結論穩定性審查.md`
- 若不處理的風險：執行者看不同文件取得不同 TODO，誤判進度

**現況證據**：
- Phase 1 Agent A 觀察 I5

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 採用 audit card 的 P0-A/P0-B/P0-C 編號為唯一標準；修兩個檔案同步 | 30 分鐘 | 低 |
| B | 僅在 06_審查.md 加「註：同 CURRENT_FOCUS P0-3」 | 10 分鐘 | 中 — 形式統一，編號仍雙系統 |
| C | 不統一 | 0 分鐘 | 高 |

**驗證標準**：
- `grep -rn "P0-3\|P0-1" docs/` 回傳一致的說明
- audit card P0-A/B/C 編號被三方引用

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P3-D: Archive 舊數據 grep 修正

**問題描述**：
早期 archive 報告仍引用 TP=30,490 / FP=4,842（矯正前）。2026-04-04 矯正後應為 TP=28,509 / FP=11,606（paired）或 TP=28,383 / FP=11,830（TO）。

**影響範圍**：
- 受影響文件：`docs/experiments/.../archive*.md`, 可能 research_landscape 早期段落
- 若不處理的風險：論文引用時抓錯數據；內部不一致

**現況證據**：
- Phase 1 Agent A 觀察 I6

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | `grep -rn "30[,.]490\|4[,.]842" docs/` + 逐處加「（舊數據，已矯正為 X）」註記 | 1 小時 | 低 |
| B | 僅在明顯位置修正（CURRENT_FOCUS / INDEX） | 20 分鐘 | 中 — archive 仍殘留 |
| C | 不改，archive 視為歷史 | 0 分鐘 | 中 |

**驗證標準**：
- `grep -rn "30[,.]490" docs/` 所有結果附近有「舊數據」標記
- `grep -rn "4[,.]842" docs/` 同樣

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P3-E: PRE-FIX 警告標示

**問題描述**：
`bip8_output_archive/` 的 39 個歷史結論中，部分使用 PRE-FIX haplotag 數據（修正前）。這些數據仍可能被部分文件引用，但未統一標示「PRE-FIX 不得引用」。

**影響範圍**：
- 受影響目錄：`bip8_output_archive/`（39 dirs）
- 若不處理的風險：未來檢視 archive 時誤引用

**現況證據**：
- Audit card 交叉驗證：`cross_cutting/CovM_Bug_Impact.md` 提及
- Phase 2 Role 3 C-INFRA-2

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 在每個 PRE-FIX 目錄新增 `_WARNING_PRE_FIX.md` 檔案（統一模板），標示「此目錄數據為 haplotag 修正前，不得直接引用於新結論」 | 2 小時 | 低 |
| B | 建立總覽 `bip8_output_archive/README.md` 標示 PRE-FIX 目錄清單 | 30 分鐘 | 中 — 檢視時可能不看 README |
| C | 不標示 | 0 分鐘 | 高 |

**驗證標準**：
- 每個 PRE-FIX 目錄內有 `_WARNING_PRE_FIX.md`
- 或 `bip8_output_archive/README.md` 列明 PRE-FIX 目錄清單

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## P3 區段總結

**推薦總動作**（全選 A）：
1. P3-A HPFineNGroups filter 更新（20 分鐘）
2. P3-B 補 17-22 結論入 06_審查（1 小時）
3. P3-C P0 編號統一（30 分鐘）
4. P3-D archive 舊數據 grep（1 小時）
5. P3-E PRE-FIX 警告（2 小時）

**總投入**：~5 小時，全部獨立可並行。

**推薦理由**：
- 全部獨立、低成本
- 文件一致性補齊後論文撰寫/reviewer response 都受益
- 其中 P3-B 最有價值（保持 06_審查 作為權威索引）

**若全選 A 後的輸出**：
- `CURRENT_FOCUS.md` / `06_結論穩定性審查.md` 更新
- `bip8_output_archive/` 新增警告檔
- archive 舊數據註解補齊

---

## 建議執行序列

```
並行：P3-A + P3-C（編輯小檔案）
並行：P3-B + P3-D（編輯中檔案）
並行：P3-E（批次建立小檔案）
```

可全部在半日完成。

---

## 關聯文件

- EXECUTIVE_DECISION_BRIEF.md 第 4-5 節
- `docs/CURRENT_FOCUS.md`, `06_結論穩定性審查.md`
- CHECKLIST.md
