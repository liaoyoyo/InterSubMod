---
id: ism-kb-00-governance-think-before-code
name: "實作前準則（Think Before Code）"
description: "開始實作/分析/決策前的檢核清單：假設陳述規則、2D 暫停判斷矩陣（可逆性 override）、目標驅動驗證（Step→Verify）、高影響場景清單、首 turn 規格完整度檢核。對齊 Opus 4.7 literal 特性。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: howto
verified_scope: "對齊 .claude/CLAUDE.md「實作前行為準則（Think Before Code）」章節"
related_ids:
  - ism-kb-00-governance-index
  - ism-kb-00-governance-confirmation-protocol
  - ism-kb-00-governance-query-protocol
  - ism-kb-06-workflows-cpp-change-pdd
tags: [governance, think-before-code, assumptions, verification, opus-4.7, matrix]
canonical_paths: [00_governance/10_think_before_code.md]
alias_paths: []
---

# 實作前準則（Think Before Code）

- 一句結論：開工前 4 步驟 —（1）列假設（2）查暫停矩陣（3）寫強驗證標準（4）首 turn 規格檢核。對齊 Opus 4.7 **不推斷未明講需求** 的特性
- 適用對象：AI agent 接到任何實作/分析/決策任務時
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -A 60 "實作前行為準則" /big7_disk/liaoyoyo2001/InterSubMod/.claude/CLAUDE.md | head -70
  ```

---

## ⚠️ 為何 Opus 4.7 下特別重要

> 「模型不會推斷未明講需求、不會悄悄泛化指令。模糊輸入將被直接按字面執行。」— CLAUDE.md

**後果**：模糊指令 → 模型按字面做 → 做錯事且無警告。

**屏障**：本頁的 4 個規則。

---

## 🧩 Step 1 — 假設陳述規則

開始任何實作、分析、或研究決策前，**必須先列出關鍵假設**。

### 原則
- **明確寫出假設**；不確定的標 `⚠ 待確認`，不要默默選擇
- **多重解讀時**：列出所有合理解讀 + 說明傾向與理由
- **主動提簡化替代方案**；必要時推回過度複雜需求
- 模糊指令 → 查 Step 2 暫停矩陣決定動作

### 範例（好 vs 壞）

❌ **壞**（默默假設）：
> 「好的，我來實作 HPFineNGroups filter。」
> （未說明 N 門檻、是否含 AF<0.4、是否 NonLOH）

✅ **好**（列假設）：
> 「開始實作 HPFineNGroups canonical filter。假設：
> - 條件：N≥4 + AF<0.4 + NR≥80 + NonLOH（依 09_Part_B.md canonical 定義）
> - 輸出欄位：加到 significance_summary.csv 的 P 群
> - ⚠ 待確認：是否要同時加到 master_dataset.csv（Python 衍生）？」

---

## 🚦 Step 2 — 暫停判斷矩陣（SoT：見 09）

**2D + 可逆性 override 的完整定義為 SoT 單一權威**：

👉 **[09_confirmation_protocol.md §2D 暫停判斷矩陣](09_confirmation_protocol.md#-2d-暫停判斷矩陣)**（完整矩陣 + 三維度定義 + Hard Gate 清單 + 一行告知格式）

本頁不重複矩陣內容，避免未來 SoT 漂移（例如 v0.4 review 發現 09/10 兩份矩陣已有「高影響」定義微差）。

**核心原則速查**：
- **不可逆** → 永遠 🔴（見 09 §Hard Gate 清單）
- **一行告知格式**：`[決策]（影響: X, 信心: Y, 理由: 一句）` — 詳 09

---

## 🚨 Step 3 — 高影響場景清單

以下場景**典型落在中/高影響區**，需有明確高信心理由才能自主；否則預設暫停：

| 場景 | 典型影響 | 為何危險 | 正確做法 |
|------|---------|----------|----------|
| 研究重點排序 | 高（>1h 投入） | 選錯浪費數天 | 列候選方向 + 預期收益/風險，等用戶決定 |
| 假說選擇 | 高（影響結論） | 隱含假設可能導致結論偏誤 | 寫出前提條件、可能 confound |
| 統計方法選擇 | 中-高（影響結論） | 不同方法可能相反結論 | 說明為何選此方法、替代方案、預期差異 |
| 「改進」/「優化」模糊指令 | 高（方向未定） | 改進方向無數 | 要求用戶定義成功標準，或提 2-3 方向供選 |
| 多檔案/多步驟重構 | 中-高（影響範圍廣） | 影響範圍不明 | 列受影響檔案 + 預期改動，確認後再動 |

---

## ✅ Step 4 — 目標驅動驗證（Step → Verify）

多步驟任務必須轉化為**可驗證執行計劃**。每步明確驗證標準：

```
1. [步驟描述] → 驗證: [具體可觀察的檢查方式]
2. [步驟描述] → 驗證: [具體可觀察的檢查方式]
...
```

### 強驗證標準（要求）

- **具體數值**：`delta_F1 > 0`、`AUC >= 0.60`
- **特定檔案輸出含路徑**：`output/canonical/HCC1395/paired_full/*/significance_summary.csv`
- **可重現測試命令**：`./scripts/run_vcf_all_snv.sh --mode chr19-verification`
- **命令退出碼**：`echo $?` = 0
- **預期輸出片段**：`grep "PASS"` 可找到

### 弱驗證標準（禁止）

- ❌ 「讓它能動」
- ❌ 「確認結果合理」
- ❌ 「看起來正確」
- ❌ 「double-check 一下」

> **Opus 4.7 備註**：模型自我驗證能力已提升，不需加「確認看起來正確」類軟性 scaffolding。驗證標準必須是**外部可觀察**的。

### 範例：加入 Normal BAM read 過濾邏輯

```
1. 讀取現有 ReadParser 了解 BAM 讀取流程
   → 驗證: 能指出 normal read 進入的函數位置和資料結構

2. 新增 normal BAM 路徑參數與讀取邏輯
   → 驗證: make -j$(nproc) 編譯通過，無 warning

3. 實作 normal read 過濾條件
   → 驗證: 單元測試覆蓋 normal-only / tumor-only / mixed 三種情境

4. 全量測試確認無回歸
   → 驗證: run_batch_vcf_analysis.sh 結果與修改前 F1 差異 < 0.01
```

---

## 📋 首 turn 規格完整度檢核

接到任務後，**先列出「我理解的任務規格」**：

| 欄位 | 內容 |
|------|------|
| **意圖** | 使用者想達成什麼？ |
| **約束** | 時間/資源限制？不可變項？ |
| **檔案路徑** | 需讀/寫哪些檔案（含行號）？ |
| **驗收標準** | 任務完成的可觀察判斷？ |
| **scope 邊界** | 不做什麼？ |

**規則**：
- **規格缺項 ≥2 且涉及高影響決策** → 🔴 必須回問
- 否則 → 列出預設選擇與理由後繼續

### 範例：任務規格列表

> **任務**：「把 HPFineNGroups 結論更新到 KB」
>
> **我理解的規格**：
> - 意圖：將 2026-04-22 HPFineNGroups 新 canonical filter 加入 KB
> - 約束：不重寫 docs/reports/research_landscape/09_Part_B.md（權威來源）
> - 檔案：`knowledge/07_derived_features/01_hpfinengroups.md`、`knowledge/09_conclusions/01_positive_findings.md`
> - 驗收：跑 3 驗證腳本 PASS；盲測能找到 NG=4+AF<0.4+NR≥80
> - scope 邊界：不改 HPFineNGroups 本身實作；僅文件更新
>
> **不確定**：⚠ 待確認 — `docs/reports/research_landscape/09_Part_B.md` 是否需同步更新？

---

## 🔗 決策流程總結

```
使用者給任務
    │
    ▼
(1) 列假設（模糊處明示 ⚠ 待確認）
    │
    ▼
(2) 查首 turn 規格完整度
    ├─ 規格缺項 ≥2 且高影響 → 🔴 回問
    └─ 否則 → 繼續
    │
    ▼
(3) 對每個決策點查暫停矩陣
    ├─ 🔴 立即暫停
    ├─ 🟠 節點暫停
    ├─ 🟡 列假設繼續
    └─ 🟢 一行告知繼續
    │
    ▼
(4) 寫強驗證標準（Step → Verify）
    │
    ▼
執行
```

---

## 🔧 相關工具對應

| 工具 | 對應階段 |
|------|----------|
| `/confirmation-protocol` skill | Step 2 暫停判斷 |
| `/methodology-audit` skill | Step 1 假設陳述（修 C++ 前）|
| `/cpp-change` skill | Step 4 驗證（PDD 6 steps 含明確驗收）|
| `/known-pitfalls` skill | Step 1 假設陳述（避免已知陷阱）|

---

## ❌ 反模式（禁止）

- 接到任務就開工，不列假設
- 「我覺得應該是 X」不查 memory / docs 直接做
- 驗證標準寫「looks good」「應該 OK」
- 模糊指令下默默選一種解讀
- 模型自判需長計算（>10 min）時不暫停

---

## 相關

- Confirmation protocol（暫停級別定義）：[09_confirmation_protocol.md](09_confirmation_protocol.md)
- Query protocol（查詢情境）：[05_query_protocol.md](05_query_protocol.md)
- PDD 6 steps（實作階段的明確驗收）：[../06_workflows/07_cpp_change_pdd.md](../06_workflows/07_cpp_change_pdd.md)
- CLAUDE.md 權威章節：`.claude/CLAUDE.md` §「實作前行為準則（Think Before Code）」
