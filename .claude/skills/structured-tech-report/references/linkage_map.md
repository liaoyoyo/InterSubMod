# 與其他 skill 的呼叫先後圖（structured-tech-report）

避免重複造輪、確保跨 skill 一致性。

---

## 上游（在本 skill 之前可能先跑的）

| Skill | 觸發情境 | 餵給本 skill 的內容 |
|------|---------|---------------------|
| `/methodology-audit` | 修改 .cpp 前的方法學審查 | §5 根因素材 + §6 候選方案 |
| `/cpp-change` 6 步驟 PDD | 已決定動 .cpp 的執行協議 | §7.2 commit hash + file:line |
| `/results-analysis` | 跑完統計分析 | §4 量化問題、§9 驗證證據 |
| `/auc-confound-guard` | AUC > 0.58 的特徵宣告 | §9 三關驗證結果 |
| `/feature-layered-observation` | 特徵 Step 1-6 觀察 | §4 / §9 圖表與 verdict |
| `/multi-sample-consistency` | 7 樣本一致性 | §9 canonical 排序表 |

## 同層（本 skill 內呼叫）

| Skill | 在哪一步呼叫 | 用途 |
|------|------------|------|
| `/doc-standards` | Step 5 登錄前 | 確認檔案命名與 INDEX 格式 |
| `/known-pitfalls` | Step 1（涉 OLS／VCF／feature 設計時） | 避免重蹈已知錯誤 |
| `/confirmation-protocol` | Step 5（涉 Hard Gate 動作時） | 取得用戶確認 |

## 下游（本 skill 之後可能再跑的）

| Skill | 觸發情境 | 從本報告承接的內容 |
|------|---------|---------------------|
| `/conclude-research` | 報告涉假說驗證收尾 | §13 結論 → evidence_ledger |
| `/weekly-report` | 本週新增重要報告 | §0 TL;DR → 週報 Layer 2 |
| `/results-report` | 報告強調實驗結果 | §8 新舊比較 → 決策章 |
| `/report` | AI session 紀錄附件 | §7.2 工程版 → session decisions |
| `/data-audit` | 報告涉新增 figures／數據 | §9 驗證 + 附錄 A → 完整性審計 |

---

## 不應重複的內容（沿用既有 skill 規範）

| 內容 | 沿用源 |
|-----|-------|
| 檔名 `YYYYMMDD_主題_NN.md` | `/doc-standards` |
| HTML comment frontmatter 欄位 | `v5_audit_suite/00_INDEX.md` 範式 |
| 圖表命名 `figXX_*.png` 與 `figures/{section}/` | `/weekly-report` |
| 7 樣本 canonical 排序 | `/multi-sample-consistency` |
| confound 三關（within-group OLS / AF-bin / permutation） | `/auc-confound-guard` |
| C++ 修改 6 步驟 | `/cpp-change` |
| 暫停判斷矩陣（Hard Gate / Gate / Review / FYI） | `/confirmation-protocol` |
| MEMORY 條目格式 | `/memory-consolidation` |

---

## 觸發鏈典型範例

### 範例 1：「整理 V5 haplotag 修改的技術報告」

```
1. 用戶觸發 /structured-tech-report
2. Step 1 定位 → 報告類別=pipeline-change，含 .cpp 改動=Y
3. 自動呼叫 /methodology-audit（讀 InterSubMod.cpp HaplotagProcess 既有設計）
4. 自動呼叫 /known-pitfalls（檢查 self-phasing 已知陷阱）
5. Step 2 盤點 → 拉出 v5_audit_suite 19 份子報告
6. Step 3 填骨架 → §9 涉 7 樣本驗證 → 呼叫 /multi-sample-consistency
7. Step 4 自檢 → 25 項 checklist
8. Step 5 登錄 → 落 validated/2026/04/，更新 INDEX
9. 完成後提示用戶：是否要呼叫 /conclude-research 把結論寫進 evidence_ledger
```

### 範例 2：「整理 ReadParser --germline-hp-only flag 的 bug fix 報告」

```
1. 用戶觸發 /structured-tech-report
2. Step 1 定位 → 報告類別=bug-fix，含 .cpp=Y
3. 自動呼叫 /cpp-change（取得 6 步驟 commit 對應）
4. Step 3 §5 用 5 Whys；§6 候選表只列 2 個
5. Step 4 自檢
6. Step 5 落 validated/，提示用戶 MEMORY 已有 `project_readparser_germline_hp_only_phase1_negative.md`，可在附錄 A 引用
```

---

## 反例（不該做）

- ❌ 把「本週多主題進度」塞進本 skill 的 13 段 → 應改用 `/weekly-report`
- ❌ 在 §6 重新定義候選方案的圖表命名規則 → 應沿用 `/weekly-report` figures/ 慣例
- ❌ 在 §9 自創統計檢定流程 → 應呼叫 `/auc-confound-guard`
- ❌ 跳過 §7.1 直接寫工程細節 → §7 雙語為硬性
