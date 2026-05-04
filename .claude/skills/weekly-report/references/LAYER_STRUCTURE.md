# Layer 0-4 五層敘事框架 + 17 段母稿 mapping + 4 層分類規則（v2 升級）

> **2026-05-02 v2 升級**：
> - **§A**（行 1-50 原內容保留）：Layer 0-4 + Tier 1/2/3 證據卡（v1 強架構，沿用）
> - **§B**（新增）：17 段母稿 → Layer 對應 mapping 表
> - **§C**（新增）：4 層分類 [F]/[O]/[I]/[U] 規則 + 與 Tier 1/2/3 並用範例
> - **§D**（新增）：4 桶分流（PPT/講稿/備註/暫存）評分規則

---

# §A 五層敘事框架 + 認知減負工具包（v1 保留）

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

---

# §B 17 段母稿 → Layer 對應 mapping（v2 新增）

母稿主骨架沿用 Layer 0-4，**17 段為 Layer 內部標籤**（不雙重結構）。

| 段 | 標題 | 對應 Layer | 內容來源（W階段）|
|:-:|------|----------|---------------|
| §1 | 本週報告主線（≤ 30 字）| Layer 0.1 | W2 → C1 |
| §2 | 本週一句話重點 | Layer 0.1 | W2 → C1 |
| §3 | 已確認內容 [F] | Layer 2 證據卡 Tier 1/2 [F] | W3 → C2 |
| §4 | 初步觀察與合理推論 [O][I] | Layer 2 證據卡 [O][I] | W3 → C2 |
| §5 | 待確認內容 [U] | Layer 2 證據卡 [U] | W3 → C2 |
| §6 | 不建議放入 PPT | 4 桶分流: 暫存 | W4 → C2 |
| §7 | 報告重點優先順序 | Layer 2 整合 | W4 → C2 |
| §8 | 建議報告順序 | Layer 2 整合 | W4 → C2 |
| §9 | 建議 PPT 模板類型 | Layer 4 / handoff 提示 | W7 → C4（指向 pptx-build templates）|
| §10 | 建議投影片架構 | Layer 4 / handoff 提示 | W7 → C4 |
| §11 | 需要補充的資料 | Layer 3 「仍然未知的」 | W4 → C2 |
| §12 | 需要製作的圖表 | Layer 3 | W4 → C2 |
| §13 | 需要補充的定義或解釋 | Layer 0.2 | W4 → C2 |
| §14 | 可用於講稿的例子 | Layer 2 / Tier 2 speaker note | W4 → C2 |
| §15 | 暫存紀錄 | 4 桶分流: 暫存 | W4 → C2 |
| §16 | 下一步行動清單 | Layer 4.1 | W7 → C4 |
| §17 | 教授可能提問 + 回答準備 | Layer 4 / 新增區段 | W6 → C3 |

**寫法規則**：
- 母稿 .md 用 Layer 0-4 為頂層 H2 結構（## Layer 0 / ## Layer 1 / ...）
- 17 段以 H3 寫在對應 Layer 下（### §1 主線 / ### §2 一句話重點 / ...）
- 範例：母稿 §3「已確認內容 [F]」放在 Layer 2 內，當 Thread 證據卡的 Tier 1 [F] 標籤項目

---

# §C 內容 4 層分類 [F]/[O]/[I]/[U] 規則（v2 新增）

## 4 層定義

| 標籤 | 全名 | 標記條件 | 必要欄位 | 描述語氣範例 |
|:-:|------|---------|---------|-------------|
| **[F]** | Fact 確定事實 | 有具體 source（檔案 path / commit hash / output csv）+ 達 validation threshold | source_path, validation_status | 「已驗證 X」「確認為 Y」 |
| **[O]** | Observation 初步觀察 | 有結果但 N 不足或未獨立驗證 | sample_count, threshold_to_promote | 「初步觀察到 X」「需 N 樣本驗證」 |
| **[I]** | Inference 合理推論 | 根據資料推測 | basis_evidence, alternatives | 「推測」「可能」「值得進一步觀察」 |
| **[U]** | Unconfirmed 待確認 | 有疑問或不確定 | what_to_check, when | 「待釐清」「需要 X 才能確認」 |

## 標籤決策樹

```
Q1: 有具體 source（檔案 path / commit hash / output csv）？
  ├─ 否 → [I] 或 [U]（看下面）
  └─ 是 → Q2

Q2: 達 validation threshold（如 N≥7 樣本 / p<0.05 / Confound 檢查通過）？
  ├─ 否 → [O]
  └─ 是 → [F]

(從 Q1 否分支)
Q3: 是否屬「根據已知資料推測」（而非「待釐清的疑問」）？
  ├─ 是 → [I]
  └─ 否 → [U]
```

## 4 層 vs Tier 1/2/3 並用範例

範例 A：Phase 1A paired-pure ΔF1=+0.0112（7 樣本驗證）
- Tier: 1（最重要，influencing main thesis）
- 標籤: **[F]**（7 樣本驗證、有 csv 來源）
- 寫法：「Tier 1 已驗證項目：7 樣本 paired-pure ΔF1=+0.0112 [F]」

範例 B：HPFineNGroups 在 5/7 樣本顯示 phasing 信號
- Tier: 2（重要 supporting evidence）
- 標籤: **[O]**（N=5 未達全 7 樣本）
- 寫法：「Tier 2 初步觀察：HPFineNGroups 5/7 phasing 信號 [O]，需補 HCC1395 + COLO829 驗證才能升 [F]」

範例 C：自我推測「這個方向可能要轉向 LOH-constrained phasing」
- Tier: 2 或 3
- 標籤: **[I]**（推論，無直接 evidence）
- 寫法：「Tier 3 推論：可能要轉向 LOH-constrained phasing [I]，依據是 NG=2 same-hap 6/6 觀察」

範例 D：「ClairS-TO 4.5.0 是否解決 LOH penalty 問題？」
- Tier: 1（影響後續方向）
- 標籤: **[U]**（待釐清）
- 寫法：「Tier 1 待確認：ClairS-TO 4.5.0 對 LOH penalty 影響 [U]，需跑全 7 樣本 benchmark 才能確認」

## 過度宣稱紅旗（W5 → C3 自動掃描）

下列用法 = 違反分類規則：

- 把 [I] 寫成「已確認」「證實」「解決」 → 過度宣稱紅旗
- 把 [O] N=3 樣本寫成「全部樣本」「跨樣本一致」 → 過度宣稱紅旗
- 把 [U] 直接寫進 §3 已確認內容 → 標籤錯位
- 證據卡 Tier 1 內容無 [F]/[O]/[I]/[U] 標籤 → 規則違反

---

# §D 4 桶分流評分規則（v2 新增）

每筆素材打 5 維度分數（1-5），加總後分桶：

| 維度 | 1 分 | 3 分 | 5 分 |
|------|------|------|------|
| 研究重要性 | 邊緣 | 中等 | main thesis 一部分 |
| 證據強度 | [U] | [I] | [F] |
| 教授關心程度 | 不太關心 | 中等 | 必問 |
| 影響下週計畫 | 無關 | 弱關聯 | 直接決定 |
| 適合簡報呈現 | 純文字無圖 | 部分視覺化 | 強視覺敘事 |

## 加總分數 → 桶

| 加總 | 桶 | 上限 | 對應 Tier |
|------|---|------|-----------|
| 18-25 | **PPT**（slide 上）| ≤ 8 筆 | Tier 1 |
| 13-17 | **講稿**（speaker note）| ≤ 15 筆 | Tier 2 |
| 8-12 | **備註**（oral-optional）| 不限 | Tier 3 |
| <8 | **暫存**（不放本次）| 不限 | (棄用本週) |

## 桶 vs Tier vs 4 層分類三維對照

| 桶 | Tier | 典型 4 層分類 | 範例 |
|---|------|-------------|------|
| PPT | 1 | [F] 為主，少數 Tier 1 [O] | main thesis + 7 樣本 ΔF1 |
| 講稿 | 2 | [F] 細節 + [O] 觀察 | 各樣本分項 + 統計補充 |
| 備註 | 3 | [I] 推論 + [O] 邊緣 | 「如果有時間可講」內容 |
| 暫存 | (棄) | [I] / [U] 弱 | 流水帳細節 |

## 上限超出處理

- PPT 桶 > 8 筆：強制降級最弱者到講稿桶（用 5 維度評分排序）
- 講稿桶 > 15 筆：降級到備註桶
- 若用戶堅持保留，標 [ESCALATION] 並 W7 母稿提示

---

# §E Highlights / Verdict 強調機制（v2.1 新增 2026-05-03）

> **解決問題**：v2.0 母稿 Layer 0-4 結構完整但「重點/結論」未突出。教授快速翻閱時，前 30 秒應能抓到本週最關鍵 3-5 件事 + 1 個 verdict。

## §E.1 母稿頂部強制 TL;DR 區段（§0）

母稿在 §1 主線之前**強制**插入 `§0 Highlights / TL;DR`：

```markdown
## §0 Highlights (TL;DR — 教授前 30 秒能抓到的本週關鍵)

⭐⭐⭐ **This Week's Verdict** (≤ 50 字)：[本週一句決定性結論]

### Verdict Detail（v2.2 新增，optional ≤ 100 字）

[若 Verdict 50 字無法涵蓋三事實以上的複雜情況，補一段 ≤ 100 字 detail。
 此段為 optional — 簡單週報可省略；複雜 multi-finding 週報建議寫。]

### Top Findings（3-5 條）
1. ⭐⭐⭐ [F] <最重要的事實，必含具體 evidence>
...

### Top Asks（教授必須判斷的決策點）
1. **[U]** <必問教授的問題 1>
... (≤ 3 條)

### Decisive Next Step（≤ 1 條）
> **Priority 1**：<下週決定性 priority，含工時估算>
```

**v2.2 修正（2026-05-05）— Verdict 50 字限制檢討**：
- Verdict ≤ 50 字硬上限保留（一句話原則，教授第一眼焦點）
- 若需涵蓋 3 事實以上（如「commit 鏈完整化 + Pass 1 only caveat + P0 行動」），加 **Verdict Detail** ≤ 100 字補充段
- Verdict Detail 為 optional：簡單 progress 週報可省略；problem 主線 + 多 finding 場景建議寫
- 經驗法則：Verdict 想超過 50 字 = 訊號「需要拆 Verdict + Verdict Detail」

**規則**：
- §0 在 §1 之前
- Top Findings ≤ 5 條（多了表示重點過散，須再壓縮）
- Top Asks ≤ 3 條（多了表示教授決策負擔過重）
- Decisive Next Step 只能 1 條（決定性 priority）

## §E.2 重要度標註規則（每段加 ⭐）

每段內容前綴加標：

| 標記 | 含義 | 適用 |
|:-:|------|------|
| ⭐⭐⭐ | 教授必看必記 | §0 Highlights、Layer 2 Thread 結論、§16 Decisive priority |
| ⭐⭐ | 重要支持證據 | Layer 2 證據卡 Tier 1、§17 必問追問 |
| ⭐ | 細節支持 | Layer 2 證據卡 Tier 2、§17 可能追問 |
| (無) | 背景 / 過程 | Layer 0.2 背景知識、§13 定義、附錄 |

## §E.3 Layer 2 Thread 結論段強制 One-line Verdict

每個 Thread 的「結論」段第一行強制為 **black bold one-line verdict**：

```markdown
#### 結論
**⭐⭐⭐ One-line Verdict**: <一句決定性判決，≤ 30 字>

- 判決細節：...
- 穩定度：...
- 影響：...
- 已排除替代解釋：...
- 重新開啟條件：...
```

## §E.4 §17 教授追問分級

7 個追問分兩級：

```markdown
## §17 教授可能提問 + 回答準備

### ⭐⭐⭐ 必問 (Must-Answer，3 個)
1. **「<追問>」** → <預備回答骨架>
2. ...
3. ...

### ⭐⭐ 可能問 (May-Ask，4 個)
4. ...
...
7. ...
```

**規則**：必問 ≤ 3 個，多了表示優先序不夠收斂。

## §E.5 §16 下週 priority 標決定性

```markdown
## §16 下一步行動清單

### ⭐⭐⭐ Decisive（決定性，影響後續路徑）
- **Priority 1 (必做，4 hr)**: <最重要，下週若沒做後續無法推進>

### ⭐⭐ Operational（執行性，預計推進）
- **Priority 2 (4 hr)**: ...
- **Priority 3 (2 hr)**: ...

### ⭐ Maintenance（維護性，可延後）
- **Priority 4 (1 hr)**: ...
```

**規則**：Decisive 必須有 1 條（不能 0），≤ 2 條（多了表示優先序混亂）。

---

# §F 混合主線規則（v2.1 新增 2026-05-03）

> **解決問題**：實務上單純 progress / problem / advisor / explore 主線少見，多為混合（如 problem 主 + progress sub）。

## §F.1 混合主線標註語法

C1 確認主線時，可寫 `主線:子線`：

```yaml
report_type: "problem:progress"   # 主線 problem，sub-thread 為本週進展
# 或
report_type: "advisor_consult:exploration"  # 主線求協助，sub-thread 為新方向 pilot
```

## §F.2 母稿 Sub-thread 標註

§1 main statement 後**強制**加「Sub-thread」段：

```markdown
# §1 主線（≤ 30 字）
[main statement]

### Sub-thread（如為混合主線必填）
- **主線**: problem — V5 變 no-op，audit chain 需重評
- **Sub-thread**: progress — 4/28-29 V5 audit + 4/29 技術報告 + 4/30 V3F ablation 已完成（未受 5/01 發現否定，仍計入週進度）
- **教授視角優先序**: 主線追問 > Sub-thread 進度
```

## §F.3 Layer 2 Thread A/B 對應

混合主線時：
- **Thread A** = 主線（依 main statement）
- **Thread B** = Sub-thread

範例：
```
Layer 2
├── Thread A: V5 變 no-op + audit chain 重評（problem 主線）
└── Thread B: V5 audit 系列工作完成歷程（progress sub-thread）
```

## §F.4 templates/ 對應

混合主線時優先載入「主線」對應 template，並在 W7 組裝時參考 sub template 的 §X 重點段建議：

| 主線:Sub | 主 template | Sub template 影響 |
|---------|------------|----------------|
| problem:progress | problem_focus.md | progress_focus.md §3 [F] 寫法（多用「已完成」語氣）|
| advisor_consult:exploration | advisor_consult.md | new_direction_explore.md §4 pilot 結果格式 |
| progress:problem | progress_focus.md | problem_focus.md §5 [U] 紅旗區塊 |

## §F.5 一份母稿至多 1 個 sub-thread

避免「3 主線混合」造成失焦。若用戶堅持 3 主線 → C1 強制要求拆兩份週報。

