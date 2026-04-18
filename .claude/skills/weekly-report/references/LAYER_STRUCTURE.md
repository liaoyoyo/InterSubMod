# 五層敘事框架 + 認知減負工具包

## Layer 0：研究脈絡與背景知識（~2-3 頁）

**0.1 宏觀問題定位**（固定骨架，每週更新數值）：
- Mermaid 脈絡圖（核心三節點 + 本週焦點）
- 核心數字表（已確認結論/NEGATIVE/懸置/開放假說/當前階段）
- 一段話研究摘要（~100 字）

**0.2 本週相關背景知識**（動態，依 Layer 2 主題從知識庫選取）：
- 最多 3 概念群組，每群組 2-3 段自包含解釋
- 必須引用知識庫來源
- 涵蓋 Layer 2 所需全部背景概念

**0.3 上週前情提要**（動態，從上週 Layer 3+4 提取）：
- 敘事橋樑：「上週確認了 [X]、關閉了 [Y]，引出本週核心問題：[焦點]。」

## Layer 1：已建立知識參考（薄參考層）

**1.1 已關閉假說參考表**（僅列與本週觸及者）
**1.2 進入本週的開放問題**

## Layer 2：本週調查（~4-8 頁，核心內容）

每個 Thread 使用統一結構：

```
### Thread [A/B/C]: [名稱]
#### 問題陳述
#### 定義區塊（最多 5 術語，必含範例值）
#### 假說與可否證條件
#### 方法
#### 證據卡（Tier 1/2/3）
#### 因果鏈圖（Mermaid）
#### 結論：判決 + 穩定度 + 影響 + 已排除替代解釋 + 重新開啟條件
```

## Layer 3：整合更新（~1-2 頁）

**3.1 結論總表更新**（本週變動粗體）
**3.2 本週新增認知**（3-5 點）
**3.3 仍然未知的**（優先序+問題+依賴+預計回答時間）

## Layer 4：未來方向（~1 頁）

**4.1 下週優先行動**
**4.2 里程碑收斂圖**（Mermaid Gantt）
**4.3 風險評估**

**附錄**：參考文件索引

---

## 證據卡三層制度

**Tier 1**（Priority 1-3，7 欄位必填）：
假說、測試、結果、可否證性、Confound 檢查、意義、圖表

**Tier 2**（Priority 4-5，4 必填 + 3 選填）：
> **假說**：[陳述] | **結果**：[觀察+數字]
> **意義**：[影響] | **可否證性**：[推翻條件]

**Tier 3**（Priority 6-7，行內標注）：
`[觀察] (Evidence: [指標]=[數值]; falsifiable if [條件])`

分配規則：Phase 1.1 Priority 對應 → 用戶可覆蓋。每 Thread 至少 1 張 Tier 1/2。

---

## 認知減負工具包

| 工具 | 用途 |
|------|------|
| 定義區塊 | Thread 開頭，最多 5 術語+範例值 |
| 因果鏈圖 | Mermaid，≤6 節點，含可否證節點 |
| Worked Example | 用真實 locus 跑完整流程，跨週沿用 |
| Before/After 表 | 修正前後比較 |
| 假說消除矩陣 | 多假說同時測試 |
| LOH 雙定義對照圖 | LOH.bed vs HP_Ratio，Kappa=0.670 |
| Assertion-Evidence 標題 | 結論句而非描述詞 |

---

## 主題→知識庫對照表

| Layer 2 主題 | 所需背景概念 | KB 來源 |
|---|---|---|
| LOH 分析 | LongPhase-TO, HP tags, LOH.bed vs HP_Ratio, self-phasing | `05_tools/longphase-to.md`, `03_file_formats/bam-format.md` |
| 甲基化特徵 | MM/ML, CpG, ISM window, PERMANOVA | `03_file_formats/bam-format.md`, `05_tools/intersubmod.md` |
| Self-phasing | 循環依賴, HP bias, PS, PON-only | `06_workflows/phasing-workflow.md` |
| 跨樣本共識 | 7 樣本特性, truth set, benchmark | `02_samples/`, `06_workflows/benchmark-workflow.md` |
| PON 過濾 | gnomAD/dbSNP/CoLoRSdb, FILTER | `04_databases/`, `03_file_formats/vcf-clairs-to.md` |
| F1/AUC 評估 | Benchmark workflow, TP/FP/FN 定義 | `06_workflows/benchmark-workflow.md` |
| Quality Score | QS 組成, LOH penalty, EHR tier | `05_tools/intersubmod.md` |

使用規則：每週最多 3 群組。PI 已知概念不重複解釋，新組合或新發現需要。

---

## 固定 vs 動態

| 固定（每週不變） | 動態（每週調整） |
|-----------------|----------------|
| 五層框架結構 | Layer 1 選取哪些假說/問題 |
| Layer 0.1 骨架 | Layer 0.2 概念群組 |
| Thread 統一模板 | Layer 0.3 前情提要 |
| 證據卡三層 Tier | Thread 數量與主題 |
| 價值排序規則 | 具體排序與 Tier 分配 |
| 圖表規範（Paired/TO 分開） | 具體圖表 |
| 7 samples 固定順序 | 敘事主軸類型 |
| 雙層檢核清單 | 重點強調項目 |
