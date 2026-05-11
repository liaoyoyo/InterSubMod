---
id: ism-kb-06-workflows-pptx-and-weekly-report
name: "PPTX 與週報流程（索引）"
description: "3 個 skill 對照（weekly-report / results-report / report）；PPTX 規範摘要；7 Phase 週報流程索引；既有 manual 清單。純索引入口，不重寫 skill 或 manual 內容。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "對應 .claude/skills/{weekly-report,results-report,report} + docs/references/manual/"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-06-workflows-analysis-scripts-index
tags: [workflow, pptx, weekly-report, reporting, index]
canonical_paths: [06_workflows/09_pptx_and_weekly_report.md]
alias_paths: []
---

# PPTX 與週報流程（索引）

- 一句結論：3 個 skill 分工（週報 / 實驗報告 / session 報告）；PPTX 規範於既有 manual；本頁只做入口，不重寫
- 適用對象：每週週報、實驗結果報告、AI session 報告
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/.claude/skills/ | grep -E "report|weekly"
  ```

---

## 🎯 3 個 Skill 分工對照

| Skill | 用途 | 觸發詞 | 輸出 | 受眾 |
|-------|------|--------|------|------|
| `/weekly-report` | **每週 PI 週報** | 「週報」、「weekly report」、「整理本週」 | `.md` 週報 + `.pptx` 簡報 | PI / 指導教授 |
| `/results-report` | **單一實驗結果報告** | 「寫實驗報告」、「results report」 | `docs/experiments/*.md` | 研究者自身 + PI |
| `/report` | **AI 對話執行報告** | Stop hook、「session report」 | `docs/provenance/ai_sessions/*.md` | AI 協作追溯 |

**選擇建議**：
- 每週例行 → `weekly-report`
- 單一實驗告一段落 → `results-report`（先用 `/results-analysis` 完成統計）
- AI 會話結束前 → `report`（Stop hook 提醒）

---

## 📊 Weekly Report — 8 Phase 流程（含 3.5）

對應 `.claude/skills/weekly-report/SKILL.md`（完整版）：

```
Phase 0  : 自動收集本週進展     [FYI]
Phase 1  : 協作規劃              [Review]
Phase 2  : Layer 0 基礎層撰寫    [FYI]
Phase 3  : Layer 1-4 週報主體    [Review]
Phase 3.5: 導演 Storyboard 審查  [Review]   ← 關鍵；不要跳過
Phase 4  : PPTX 三件套 + 生成    [FYI]       ← PPTX 在此，非 Phase 5
Phase 5  : 多 Agent 驗證         [FYI → Review]
Phase 6  : 最終檢核與索引更新    [Review]
```

**互動模式**：Review 暫停點（1, 3, 3.5, 5, 6）等確認
**全自動模式**：FYI 自動過；Review 5/6 強制展示（skill 內定義）

**PPTX 主腳本**：`scripts/analysis/build_weekly_report_pptx.py`（+ `build_weekly_report_pptx_oral_v2.py` 口試版 / `build_loh_weekly_pptx.py` LOH 主題版）

> **重要澄清**：各 `docs/presentations/<run>/build_pptx.py` 為**週報實例內**由模板複製的一次性檔；**全域** PPTX 主腳本是 `scripts/analysis/build_weekly_report_pptx.py`。

---

## 🎨 PPTX 規範摘要（詳見 manual）

### 設計哲學
- **視覺優先**：物件圖示為主、文字為輔（6:4 或 7:3 比例）
- **每頁 2-3 焦點**：避免資訊密集
- **結論句非數據**：「ΔF1=+0.0112 locked」比「0.0112」更有意義

### 格式規範
- **標題 ≤15 字**
- **禁用斜體**（適合標題可粗體）
  - 理由：PPTX 斜體在多數系統字型（尤其 CJK）下可讀性差；粗體對比更高
- **雙語**：EN 縮排 0.25"、60% 字號
- **首次術語**：必配圖 + 名詞解釋區塊
- **圖示比例**：60-70% 版面

### 審查
- **Storyboard 審查**：先信任內容→再邏輯分析→最後根因
- **導演審查模式**：整場 flow 先看過再調細節

### 修改前 diff 規則
⚠️ **改 `build_pptx.py` 前必先 diff**，偵測手動編輯（避免覆蓋）
- MEMORY：`feedback_pptx_presync_workflow`

---

## 📚 既有 Manual 權威清單

| Manual | 用途 | 路徑 |
|--------|------|------|
| 週報生成 Skill 規格 | 週報 skill 設計文件 | [`docs/references/manual/20260310_週報生成Skill規格_01.md`](../../docs/references/manual/20260310_週報生成Skill規格_01.md) |
| 週報撰寫指令與 skill 草案 | 週報撰寫流程 | [`docs/references/manual/20260310_研究週報撰寫指令與skill草案_01.md`](../../docs/references/manual/20260310_研究週報撰寫指令與skill草案_01.md) |
| PPTX 客製化設定與製作手冊 | PPTX 完整規範 | [`docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md`](../../docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md) |
| 投影片產生規範完整整理 | 最新投影片規範（2026-04-22）| [`docs/references/manual/20260422_投影片產生規範完整整理_01.md`](../../docs/references/manual/) |

---

## 📂 歷史週報與範例位置

```
docs/presentations/YYYY/MM/
├── 20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.pdf
├── 20260311_研究主線週報_..._oral_v01.pdf  # 口試版變體
└── ...
```

**檢視範例**：
```bash
ls /big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/ -R | head -20
```

---

## 🧩 相關腳本

| 腳本 | 用途 |
|------|------|
| `scripts/analysis/build_weekly_report_pptx.py` | ★ **週報 PPTX 主腳本**（Phase 4 使用） |
| `scripts/analysis/build_weekly_report_pptx_oral_v2.py` | 口試版 PPTX |
| `scripts/analysis/build_loh_weekly_pptx.py` | LOH 主題週報 PPTX |
| `scripts/analysis/build_ai_f1_research_pptx.py` | AI F1 研究 PPTX |
| `scripts/analysis/generate_pi_report_figures_self_phasing.py` | PI 報告用 Self-Phasing 圖表 |
| `scripts/analysis/extract_ppt_text.py` | 從 PPTX 抽文字（回歸驗證） |
| `docs/presentations/<run>/build_pptx.py` | **每個週報實例內的一次性客製化腳本**（模板複製，非全域）|

**完整腳本索引** → [08_analysis_scripts_index.md](08_analysis_scripts_index.md)

---

## ⚠️ 常見陷阱（MEMORY 記錄）

1. **改 build_pptx.py 前不 diff**：覆蓋手動編輯（MEMORY: `feedback_pptx_presync_workflow`）
2. **PPTX 斜體**：禁用；用粗體替代（MEMORY: `feedback_pptx_text_formatting_rules`）
3. **首次術語未配圖**：讀者無法理解（MEMORY: `feedback_pptx_term_explanation_rule`）
4. **Storyboard 跳過**：直接進細節調整（MEMORY: `feedback_pptx_director_storyboard`）

---

## 📅 週報頻率與排程

- **頻率**：每週（Fri-Sat 週末整理）
- **時間範圍**：Mon-Sun（週一到週日）
- **受眾**：指導教授/PI（熟悉 ONT、cancer genomics、somatic calling 領域）
- **自包含原則**：每份週報必須獨立可讀，PI 不需翻閱前期報告

---

## 🔗 執行範例

### 情境 1：本週例行週報
```
User: 「/weekly-report」
Agent: 進入 weekly-report skill，啟動 7 Phase 流程
```

### 情境 2：某實驗告一段落，寫結果報告
```
User: 「這個 HPFineNGroups 實驗寫一下 results report」
Agent:
  1. 先跑 /results-analysis 確認統計已完成
  2. 啟動 /results-report skill
  3. 輸出 docs/experiments/YYYY-MM-DD_*.md
```

### 情境 3：AI 會話結束前記錄
```
Stop hook 觸發: 「請撰寫執行報告」
Agent: 啟動 /report skill → docs/provenance/ai_sessions/
```

---

## 相關

- Workflows 索引：[00_index.md](00_index.md)
- Evidence ledger（週報資料來源）：[../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)
- Analysis scripts（PPTX 圖表生成）：[08_analysis_scripts_index.md](08_analysis_scripts_index.md)
- 權威 skill：`.claude/skills/weekly-report/SKILL.md`、`.claude/skills/results-report/SKILL.md`、`.claude/skills/report/SKILL.md`
