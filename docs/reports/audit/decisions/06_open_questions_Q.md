# 06 — Open Questions Q-id（10 題）

> **建立日期**: 2026-04-19

> **優先等級**：戰略探索 / 方法論 / 生物學
>
> **依賴**：多數獨立；Q-01/Q-05 依賴 P0 完成

---

## 問題 Q-01: E→A+D→B→C 優先序在 CovM bug 後是否成立？

**問題描述**：
先前研究路徑規劃 `E (infra fix) → A+D (Normal BAM + ASM) → B (Phase 1B) → C (LOH cohort)` 基於「CovM bug 修 + ISM 漸進整合」假設。但 CovM bug 修正後實際重跑結果可能顯示：
- 若 C17/C20/C22 受影響大 → 需先補數據再進 A+D
- 若影響小 → 可直接啟動 A+D

目前無法預測結果規模。

**影響範圍**：
- 受影響優先序：P0/P1/P2 任務排序
- 受影響功能：F2/F3/F4 的啟動時機
- 若不處理的風險：A+D 在汙染 baseline 上啟動，或等太久浪費時程

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **維持 E→A+D→B→C，但在 P0-A 完成後設決策點**：rank correlation >0.8 → 繼續；<0.8 → 暫停並重估 | 依決策點結果 | 低 — 柔性策略 |
| B | 立即重排為 A+D 優先（跳過 E） | 現在 | 高 — 違背 R-01 |
| C | E 完全完成（含 P0-B/P0-C + 全部重算）才進 A+D | +2 週 | 中 — 保守但慢 |

**驗證標準**：
- **選 A**：P0-A 完成後產出 rank correlation 報告
- **門檻**：>0.8 直接 A+D；0.5-0.8 部分重算後 A+D；<0.5 全量重算

**依賴關係**：
- **前置依賴**：P0-A 執行結果
- **被阻塞項**：Phase 2 A+D 啟動

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-02: Zone 若失效是否重新定義？

**問題描述**：
若 P0-A CovM 修正後，C22 Zone-Aware Z1-Z5 zone cardinality 大幅變動（例 ±50%），現有 zone 定義可能失效。是否重新定義 zone？

**影響範圍**：
- 受影響結論：C19, C20, C22（全 Zone 相關）
- 受影響功能：F3 Zone-Aware
- 若不處理的風險：Zone 框架崩潰則 F3 功能不成熟

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | P0-A 後先觀察 zone cardinality 變化；若 ±30% 內維持；±30%-50% 調整 boundary；>50% 重新定義 | 依觀察 | 低 — 柔性 |
| B | 預先定義新 zone（依 KDE-based expected_coverage） | 1 週 | 中 — 過度預備 |
| C | 不變 zone 定義，只更新 zone 成員 | 1 天 | 中 — zone 定義可能失效 |

**驗證標準**：
- zone cardinality 變化 <30% → 維持
- 30-50% → 調整 boundary（Z2/Z3 邊界）
- >50% → 重新定義（新 criteria）

**依賴關係**：
- **前置依賴**：P0-A 完成
- **被阻塞項**：F3 Zone-Aware 論文章節

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-03: ISM characterization 具體下游應用？

**問題描述**：
ISM 定位為 "read-level epigenetic characterization"，但具體下游應用尚未定：
- 臨床：subclonal dynamics monitoring（治療反應追蹤）？
- 生物學：tumor evolution mechanism study？
- 研究工具：general-purpose ASM profiler？

不同下游決定論文 Discussion 章節方向。

**影響範圍**：
- 受影響功能：F1-F5 的應用定位
- 若不處理的風險：論文 Impact 章節空泛

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **三方向並列敘述**：臨床 + 生物學 + 工具（依論文讀者群調整權重） | 1-2 週論文撰寫 | 低 |
| B | 僅聚焦臨床應用（最吸引 Genome Biology） | 2-3 週（需 cohort pilot） | 中 |
| C | 僅工具定位（最務實，適合 Bioinformatics） | 1 週 | 中 |

**驗證標準**：
- 論文 Discussion 章節涵蓋至少 2 個具體應用場景
- 每場景含 1-2 篇相關文獻引用

**依賴關係**：
- **前置依賴**：R-03 期刊選擇
- **被阻塞項**：論文 Discussion

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-04: TO 甲基化增益為負是否獨立 negative result 章節？

**問題描述**：
Phase 1A 結論：TO 模式甲基化增益為負（delta F1<0）。這是重要 NEGATIVE finding，但目前散落於 C01/C03/C10 多卡。是否彙整為獨立 negative result 章節（論文 Results subsection）？

**影響範圍**：
- 受影響結論：C01, C03, C10, C14
- 受影響功能：F5 Variant Confidence 的完整性
- 若不處理的風險：論文 NEGATIVE 敘事分散，reviewer 難理解

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 論文新增 Results subsection「TO-mode Methylation Gain: A Negative Result」整合 C01/C03/C10/C14 | 2-3 天 | 低 — NEGATIVE 彙整標準做法 |
| B | 保留散落於各卡，不整合 | 0 天 | 中 — reviewer 難讀 |
| C | 放入 Discussion 的 Limitations 章節 | 1 天 | 中 — 降低 negative result 貢獻度 |

**驗證標準**：
- Subsection 含明確機制解釋（self-phasing circular dependency）
- 含 7 樣本驗證一致性圖表

**依賴關係**：
- **前置依賴**：無（可立即規劃）
- **被阻塞項**：論文 Results 章節結構

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-05: Phase 2 A+D 成功標準？

**問題描述**：
Phase 2 A+D（Normal BAM + Sample ASM）啟動後，何謂「成功」？目前無量化門檻：
- Sample ASM map 產出且生物學合理？
- 識別 ≥N 個候選 subclonal ASM regions？
- 與 LOH.bed subclone signal 一致性（Jaccard ≥ X）？

無成功標準 → 階段結束判定模糊。

**影響範圍**：
- 受影響功能：F4 ASM Profiling, F2 Subclone Marker
- 若不處理的風險：A+D 做完不知是否成功，繼續投入或停損模糊

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 先定義 3 層標準：**M1 最低**（ASM 分佈圖產出 + 7/7 樣本無 crash）、**M2 中階**（somatic ASM 數量 > germline baseline 15-30% 上限）、**M3 理想**（與 LOH.bed subclone signal Jaccard ≥0.3） | 1 天（定義）| 低 |
| B | 不預定標準，做完再看 | 0 天 | 高 — 事後合理化 |
| C | 僅最低標準（M1） | 1 天 | 中 — 缺乏深度 |

**驗證標準**：
- M1 達成 → 繼續 B 階段（Phase 1B）
- M1+M2 達成 → 論文 F4 章節可寫
- M1+M2+M3 達成 → F4 為論文主力貢獻

**依賴關係**：
- **前置依賴**：P0-A/P0-C 完成（乾淨 baseline）
- **被阻塞項**：B/C 階段啟動時機

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-06: 13 NO-GO 結論是否值得 ~2 週 within-group validation？

**問題描述**：
R-04 已決策（選 A：接受現有）。但此 Q-06 為 R-04 的延伸討論：**是否某些特定 NEGATIVE 結論特別值得補驗證**（例：F5 核心的 C08 Read-FP）？

**影響範圍**：
- 受影響結論：C08（Read-FP）為主
- 受影響功能：F5 Variant Confidence 的 NEGATIVE 敘事穩固度

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **僅 C08 within-group 驗證**（F5 主要 NEGATIVE 證據） | 2-3 天 | 低 |
| B | 全 13 NO-GO within-group | 2 週 | 中 |
| C | 不補任何 | 0 天 | 中 |

**驗證標準**：
- C08 within-group AUC 補齊
- 若 within-group AUC 仍 <0.6 → NEGATIVE 穩固

**依賴關係**：
- **前置依賴**：P0-B 方法論成熟後
- **被阻塞項**：F5 Results 章節

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-07: HCC1395 單樣本結論是否需其他 5 樣本驗證？

**問題描述**：
多個結論主要在 HCC1395 驗證（C12 ASM, C19 Zone 細節, C17 Subclone step3）。若 reviewer 問「僅 HCC1395 驗證能否代表其他癌型？」，目前無回答。

**影響範圍**：
- 受影響結論：C12, C17, C19
- 受影響功能：F2/F3/F4
- 若不處理的風險：論文 generalizability 被質疑

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 核心結論（C17 LOH Subclone）已 7/7 驗證；其他單樣本結論補「HCC1395 為代表性 example，其他樣本模式一致（見 C17）」敘述 | 1-2 天 | 低 |
| B | 全單樣本結論在 COLO829 + H2009 再驗證 | 2 週 | 中 |
| C | 只敘述 HCC1395 結論，不做 generalizability 聲明 | 0 天 | 中 |

**驗證標準**：
- C12/C17/C19 每卡補「cross-sample context」段落

**依賴關係**：
- **前置依賴**：無
- **被阻塞項**：論文 Discussion Limitations

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-08: ASM germline/somatic 區分是否依賴 Normal BAM？

**問題描述**：
P2-C 已問 ASM germline vs somatic 區分方式。此 Q-08 進一步：**是否 Normal BAM 是唯一可靠方式**？可能替代：
- 公開 WGBS reference（ENCODE 等）
- imprinting region 先驗（chr15q11-q13, chr11p15 等）

**影響範圍**：
- 受影響結論：C12 ASM
- 受影響功能：F4 ASM
- 若不處理的風險：Phase 2 A 時若 Normal BAM 取得困難，無備案

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **Normal BAM 為主**（高精度），imprinting 先驗為補充（低 cost quick check） | 2-3 週（Normal BAM 整合） | 低 |
| B | 僅 imprinting 先驗（不啟動 Normal BAM） | 3-5 天 | 中 — 精度不足 |
| C | 僅 Normal BAM（不考慮 imprinting） | 2-3 週 | 中 — 漏補 quick check |

**驗證標準**：
- Normal BAM 整合後，germline/somatic ASM split TSV 產出
- imprinting regions 標註統計（chr15q11-q13 hit rate 等）

**依賴關係**：
- **前置依賴**：Phase 2 A Normal BAM 資源
- **被阻塞項**：F4 完整性

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-09: NG=3 非單調：生物機制 vs artifact？

**問題描述**：
P2-D 已問（選 A：補假說 + NR=120 pilot）。此 Q-09 延伸：**若 NR=120 pilot 顯示 NG=3 仍非單調**，如何 frame biology？可能機制：
- **Hypothesis 1**：2 mutation events 之間的 branching 中間態本就不穩
- **Hypothesis 2**：NG=3 對應 subclone size ~25-33%，處於 detection boundary

**影響範圍**：
- 受影響結論：C16
- 受影響功能：F2
- 若不處理的風險：NG=3 無解釋 → reviewer 懷疑整個 NGroups metric 的 biological meaning

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **Pilot 決定方向後**，若仍非單調 → 採 H1 敘事（branching dynamics）；若 pilot 後單調 → 標為 artifact | 1 週（等 pilot） | 低 — 依數據決定 |
| B | 預先採 H1 敘事（不等 pilot） | 0 天 | 中 — 數據若不支持會尷尬 |
| C | 不解釋機制，只展示 filter 有效 | 0 天 | 中 — 論文深度不足 |

**驗證標準**：
- Pilot 後於 C16 補「Biological Mechanism」段落
- 段落引用 branching evolution 文獻（2-3 篇）

**依賴關係**：
- **前置依賴**：P2-D pilot（NR=120）
- **被阻塞項**：C16 完整 narrative

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 Q-10: LOH Subclone 臨床 downstream cohort？

**問題描述**：
C17（LOH Subclone AF × Methylation POSITIVE）是最強 F2 證據，但為 cell line 數據（HCC1395 等 7 樣本）。臨床 downstream 是否需 patient cohort 驗證？

**影響範圍**：
- 受影響結論：C17
- 受影響功能：F2 Subclone Marker
- 若不處理的風險：reviewer 問 "Does this work in real patient samples?"

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **現階段不補 patient cohort**；論文敘述為「cell line proof-of-concept」+ Discussion 提及「future patient cohort validation needed」 | 0 天 | 低 — 誠實定位 |
| B | 啟動小型 patient cohort pilot（3-5 例 ONT tumor samples） | 3-6 個月 | 中 — 資源大 |
| C | 僅引用公開 patient WGBS 數據間接驗證（ENCODE cancer samples） | 1-2 週 | 中 — 非 ONT 數據，說服力有限 |

**驗證標準**：
- 論文 Discussion 含 "Limitation: cell line only" 段落
- "Future work" 段落明確 patient cohort 需求

**依賴關係**：
- **前置依賴**：R-03 期刊選擇（Genome Biology 可能要求 B）
- **被阻塞項**：F2 論文投稿時機

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## Q 區段總結

**推薦總動作**（全選 A）：
1. Q-01 柔性策略（P0-A 後看 rank correlation）
2. Q-02 zone 變化 ≤30% 維持
3. Q-03 三方向並列敘述
4. Q-04 新增 NEGATIVE subsection
5. Q-05 3 層標準（M1/M2/M3）
6. Q-06 僅 C08 within-group
7. Q-07 核心結論 7/7 驗證，其他補 cross-sample context
8. Q-08 Normal BAM 為主 + imprinting 補充
9. Q-09 pilot 後依數據決定
10. Q-10 cell line proof-of-concept + future work

**推薦理由**：
- 多數 Q 題為「做分析時的選項」，現階段不需硬決策
- 延後決策（A 選項普遍為「觀察後再決」）最符合 Opus 4.7 literal 特性 — 不預判
- 與用戶「分析完整性優先」對齊

**若全選 A 後的輸出**：
- P0 執行報告含 Q-01/Q-02 決策觸發點
- 論文大綱含 Q-03/Q-04 章節規劃
- Phase 2 A+D 啟動報告含 Q-05 M1/M2/M3 標準
- C08 within-group 執行（Q-06）

---

## 關聯文件

- EXECUTIVE_DECISION_BRIEF.md 第 4-5 節
- `cross_cutting/CovM_Bug_Impact.md`, `Characterization_Functions.md`
- CHECKLIST.md
