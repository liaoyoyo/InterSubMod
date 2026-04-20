# 05 — Strategic R-id Decisions（4 題）

> **建立日期**: 2026-04-19

> **優先等級**：戰略（R-03/R-04/R-05/R-06）
>
> **依賴**：R-03 獨立；R-04/R-05 需 P0-B 結果；R-06 獨立可漸進

---

## 問題 R-03: 期刊目標選擇

**問題描述**：
論文投稿目標期刊尚未定。候選包含 Genome Research（IF~9）、Genome Biology（IF~13）、Bioinformatics（IF~5）。期刊選擇影響：
- 論文 scope 與 framing（methodology-heavy vs biology-heavy）
- 需補做的實驗範圍（F3 Zone-Aware / F4 ASM 深度）
- 時程（Genome Biology 要求更嚴需更長修正）

用戶 2026-04-18 已明示「分析完整性 + 對領域貢獻 > 投稿速度」（論文定位降優先）。

**影響範圍**：
- 受影響結論：全 22 結論（framing 角度不同）
- 受影響功能：F1-F5 全部（投稿對象決定哪些 function 強調）
- 若不處理的風險：繞遠路做非必要實驗；或低估 reviewer 要求範圍

**現況證據**：
- Memory: `feedback_paper_positioning_de_prioritized.md`
- CURRENT_FOCUS.md 當前階段

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **暫不決定期刊**，先補完 Phase 2 A+D（Normal BAM + ASM）再評估適配度；3 個月後重議 | 3 個月 | 低 — 與用戶「分析完整性優先」對齊 |
| B | 現在鎖定 Genome Biology，按其要求補做 F3/F4 深度實驗 | 4-6 個月 | 中 — 可能做過頭 |
| C | 現在鎖定 Bioinformatics（較快），用現有結果投稿 | 1-2 個月 | 高 — 論文貢獻不完整 |

**驗證標準**：
- **選 A**：Phase 2 A+D 完成後，對照三期刊近兩年類似主題論文（2-3 篇 sample）做 gap analysis
- **選 B**：Genome Biology 最近 3 篇 long-read methylation 論文要求範圍盤點
- **選 C**：現有結果是否足以通過 Bioinformatics 審查（F1-F5 ≥3 顆星）

**依賴關係**：
- **前置依賴**：無（策略決策）
- **被阻塞項**：論文 Introduction / Discussion 的 framing

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 R-04: 13 NO-GO 結論是否回溯 within-group 驗證

**問題描述**：
13 個 NO-GO/NEGATIVE 結論（C02 O11, C04/C05/C06 O11/O12/O13, C07 G1-G7, C08 Read-FP, C10 TO-pure, C14 Option C, others）中，部分使用 pooled OLS / n_reads confound 檢查。若將來 reviewer 問「這些 NEGATIVE 結論是否因 pooled OLS 陷阱而被低估」，目前無回溯驗證。

**影響範圍**：
- 受影響結論：C04, C05, C06（O11/O12/O13 pre-stratify AUC）、C07（G1-G7）、C10（TO-pure）、C14（Option C 雙路）
- 受影響功能：F2/F5 的 NEGATIVE 敘事完整性
- 若不處理的風險：reviewer 要求 "Have you done within-group OLS for NEGATIVE claims?" — 若無則 NEGATIVE 結論動搖

**現況證據**：
- Audit card: 上述各卡
- `cross_cutting/Pooled_OLS_Audit.md` 已列 C04/C05/C06 已知 NEGATIVE（無需重做）

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **接受現有 NEGATIVE 結論**；用戶已在 C04/C05/C06 驗證過 n_reads confound；補一行於每 NEGATIVE 卡「Within-group OLS not re-run: n_reads confound already identified」 | 1 小時 | 低 — 最省時；審查理由充分 |
| B | 全 13 NO-GO 結論統一 within-group OLS 驗證 | 2 週 | 中 — 時間大量投入 NEGATIVE 驗證 |
| C | 僅 C10 / C14 做 within-group（這兩個最有 pooled OLS 風險） | 3-5 天 | 中 — 折中 |

**驗證標準**：
- **選 A**：每 NEGATIVE 卡 D2/D5 段新增一行說明
- **選 B**：13 卡產出 `within_group_ols_rebuttal.tsv`，各結論 within-group AUC 與 pooled AUC 差異
- **選 C**：C10 + C14 within-group AUC 寫入

**依賴關係**：
- **前置依賴**：P0-B 結果（如 C17 within-group 大幅反轉，則需擴大至 R-04 選 B）
- **被阻塞項**：reviewer response 穩固度

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 R-05: Pooled OLS 全面回溯範圍

**問題描述**：
Pooled OLS residualization 陷阱（見 `feedback_pooled_ols_residualization_trap.md`）影響範圍除了 P0-B 的 C15/C16/C17，是否還有其他結論採用 pooled OLS？需統一盤點。

**影響範圍**：
- 全專案所有使用 `residualize on AF/n_reads/CovM` 的 scripts
- 受影響功能：全 F1-F5（任何使用 residualization 的分析）
- 若不處理的風險：未知數量結論受影響；無法向 reviewer 保證「all pooled OLS identified」

**現況證據**：
- `cross_cutting/Pooled_OLS_Audit.md` 列出已知使用點
- scripts: `scripts/analysis/*residualize*.py`

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **grep 掃描** `scripts/analysis/` + `output/synthesis/research_rounds/` 找所有 residualize 使用；分類為「already within-group」/「pooled-only」；清單寫入 cross_cutting | 1 天 | 低 — 純盤點，無重跑 |
| B | 盤點 + 全部 pooled-only 使用重跑 within-group | 2-3 週 | 中 — 時間投入大 |
| C | 僅 P0-B 處理 C15/C16/C17，其他結論維持現狀 | 0 天 | 中 — 潛在盲點 |

**驗證標準**：
- **必達**：`docs/reports/audit/cross_cutting/Pooled_OLS_Global_Inventory.md`（新）列所有使用點
- **驗收指令**：`grep -rn "residualize\|OLS\|residual" scripts/ output/ docs/reports/ | grep -v "within.group\|within_group"` 分類

**依賴關係**：
- **前置依賴**：無（獨立盤點）
- **被阻塞項**：R-04 NEGATIVE 卡補註的證據基礎

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 R-06: 知識庫交叉驗證範圍

**問題描述**：
22 結論中僅部分已對照 `/big8_disk/liaoyoyo2001/knowledge/` 知識庫查證（例如 C12 ASM 對照 02_DNA_Methylation_Biology）。其餘結論（尤其 F3 Zone-Aware / F5 Variant Confidence）與文獻對照不完整。

**影響範圍**：
- 受影響結論：C01-C22（部分）
- 受影響功能：論文 Discussion / Related Work 章節完整性
- 若不處理的風險：reviewer 發現已知文獻未引用；論文原創性被質疑（實際是已有現象）

**現況證據**：
- 00_INDEX.md D6 統計：僅 ~40% 結論有明確知識庫引用
- knowledge MCP server 可用

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **漸進補齊**：優先 F1-F2 核心結論（C09/C16/C17）深度查證；其餘 22 卡依投稿需要補 | 1 週 | 低 — 與論文撰寫同步 |
| B | 全 22 結論系統化查證（每卡 2-3 文獻對照） | 2-3 週 | 中 — 時間投入大 |
| C | 僅在特定結論被 reviewer challenge 時補 | 0 天 | 中 — 被動，可能延誤 |

**驗證標準**：
- **必達**：F1-F2 核心結論（C09/C16/C17）每卡至少 3 篇知識庫文獻引用
- **觀察**：D6 覆蓋率從 ~40% 提升至 ≥70%
- **驗收指令**：`grep -rn "knowledge_get_doc\|knowledge:" docs/reports/audit/cards/` 命中覆蓋率

**依賴關係**：
- **前置依賴**：R-03 期刊選擇（決定引用深度）
- **被阻塞項**：論文 Discussion 章節

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## R 區段總結

**推薦總動作**（全選 A）：
1. R-03 暫不鎖期刊，Phase 2 完成後評估
2. R-04 接受現有 NEGATIVE；每卡加「within-group not re-run, n_reads confound known」
3. R-05 grep 盤點 pooled OLS 全部使用點
4. R-06 F1-F2 核心結論文獻深度補齊

**推薦理由**：
- R-03 與用戶「分析完整性 > 投稿速度」對齊
- R-04 選 A 最省時且審查理由充分
- R-05 選 A 先盤點，再決定是否重跑（避免過度投入）
- R-06 與論文撰寫同步，不做重複投入

**若全選 A 後的輸出**：
- `cross_cutting/Pooled_OLS_Global_Inventory.md`（新）
- 13 NEGATIVE 卡補註說明
- F1-F2 核心卡文獻引用補強
- R-03 延後 3 個月再議

---

## 關聯文件

- EXECUTIVE_DECISION_BRIEF.md 第 4-5 節
- `cross_cutting/Pooled_OLS_Audit.md`
- `feedback_paper_positioning_de_prioritized.md`
- CHECKLIST.md
