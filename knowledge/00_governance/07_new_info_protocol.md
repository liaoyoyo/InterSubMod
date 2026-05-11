---
id: ism-kb-00-governance-new-info-protocol
name: "新資訊補充與驗證協議"
description: "新結論、新 feature、新參數、新樣本、方法變更進入 KB 的標準流程；含 3 層驗證（事實、一致性、覆蓋度）、SoT（Single source of truth）規則、多 agent 驗證模板。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: howto
verified_scope: "new info ingestion protocol design for KB v0.2+"
related_ids:
  - ism-kb-00-governance-index
  - ism-kb-00-governance-update-workflow
  - ism-kb-00-governance-freshness-policy
  - ism-kb-00-governance-query-protocol
  - ism-kb-00-governance-confirmation-protocol
tags: [governance, new-info, protocol, validation, sot, ingest]
canonical_paths: [00_governance/07_new_info_protocol.md]
alias_paths: []
---

# 新資訊補充與驗證協議

- 一句結論：新資訊進 KB 要經 **分類 → SoT 確認 → 3 層驗證 → 入庫 → 跨驗證**；數字類必有唯一來源，避免 22 處重複問題
- 適用對象：任何要更新 KB 內容的人或 AI agent
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 跑完整驗證（3 腳本）
  cd /big7_disk/liaoyoyo2001/InterSubMod/knowledge
  /usr/bin/python3 scripts/validate_frontmatter.py . && \
    /usr/bin/python3 scripts/check_related_ids_symmetry.py . && \
    /usr/bin/python3 scripts/check_canonical_paths.py .
  ```

---

## 🗂️ 新資訊分類（決定走哪條流程）

| 類別 | 範例 | 流程 | Freshness 敏感度 |
|------|------|------|-----------------|
| **A. 事實型** | CLI 參數異動、file:line、path、數值 | [→ 流程 A](#流程-a事實型) | 🔴 高（須精確） |
| **B. 結論型** | 新 positive finding、狀態翻轉 (positive→NEGATIVE)、方向 NO-GO | [→ 流程 B](#流程-b結論型) | 🟡 中 |
| **C. 元件型** | 新 pipeline、新 sample、新 feature、新 workflow | [→ 流程 C](#流程-c元件型) | 🟡 中 |
| **D. 狀態型** | CURRENT_FOCUS 變動、active hypothesis 增減、阻塞狀態 | [→ 流程 D](#流程-d狀態型) | 🔴 高（2 週有效） |
| **E. 結構型** | 新增目錄、改 frontmatter schema、新增腳本 | [→ 流程 E](#流程-e結構型) | 🟢 低 |

---

## 🔑 核心原則：Single Source of Truth (SoT)

每個「事實」（數字、路徑、版本、決策）**在 KB 內只有一個權威位置**。其他位置只做連結或摘要。

### SoT 對照表（避免重複）

| 事實 | 權威位置 | 其他位置的做法 |
|------|----------|---------------|
| ΔF1 locked 數字 | [03_pipelines/01_paired_full.md#f1-baseline](../03_pipelines/01_paired_full.md) | 「ΔF1=+0.0112（locked，見 01_paired_full.md）」無需複製 CI |
| CLI 參數表 | [04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md) | 連結即可，不要複製參數表 |
| Truth set 路徑 | [08_truth_and_benchmark/01_truth_set_registry.md](../08_truth_and_benchmark/01_truth_set_registry.md) | 單一 sample 頁用「見 registry」 |
| 距離度量公式 | [04_parameters/02_distance_metrics.md](../04_parameters/02_distance_metrics.md) | 只用名字引用 |
| HPFineNGroups 完整條件 | [07_derived_features/01_hpfinengroups.md](../07_derived_features/01_hpfinengroups.md) | 其他頁只寫「NG=4+AF<0.4+NR≥80 canonical filter」 |
| 研究結論詳情 | `docs/reports/research_landscape/` | KB `09_conclusions/` 只做狀態索引 |
| MEMORY 個人偏好 | `/bip7_disk/.../memory/` | KB 不複製，只引用 ID |

**違反後果**：更新時需同步 N 處，容易遺漏產生不一致（例如 2026-04-22 發現的 「ΔF1 22 處重複」問題）

---

## 🔍 3 層驗證（所有類別都必經）

### Layer 1 — 事實驗證（Facts）
- 檔案/路徑存在：`ls` / `test -f`
- 原始碼 file:line：`grep -n "<符號>" <file>` 確認
- 數值來源：引自哪份 docs / evidence_ledger / output？
- 命令可執行：非破壞性則實跑一次（例：`--help`、`ls`）

### Layer 2 — 一致性驗證（Consistency）
- **與 docs/** 是否一致（尤其 `docs/reports/research_landscape/`、`docs/CURRENT_FOCUS.md`）
- 與 **MEMORY** 是否一致（`MEMORY.md` Active/Concluded 記載）
- 與 **既有 KB 頁面** 是否一致（跨頁面的同一數字、同一結論）
- 腳本通過：`validate_frontmatter` / `check_related_ids_symmetry` / `check_canonical_paths`

### Layer 3 — 覆蓋度驗證（Coverage）
- 新概念是否有索引入口（對應 `00_index.md` 更新）
- 是否在 `README.md` 快速導航表出現（若高頻查詢）
- `related_ids` 是否雙向完整
- 是否被 `AGENT.md` 或 `05_query_protocol.md` 的 5 種情境涵蓋

---

## 流程 A：事實型（CLI 參數 / file:line / 數值）

```
1. 讀原始碼確認（grep / Read 定位精確行號）
2. 判斷 SoT 位置（通常在 04_parameters/ 或 05_data_formats/）
3. 若 SoT 已存在 → Edit 更新該處；其他頁面是否需改連結？
4. 若 SoT 不存在 → 決定放哪份；考慮是否新增檔案
5. 更新 last_verified；精確描述 verified_scope
6. 跑 Layer 1 驗證：命令實際執行一次
7. 跑 3 腳本
8. 若其他 KB 頁面引用此數值 → grep 全庫檢查是否需同步
```

**範例（Phase 6 已實施）**：
- 發現 `RegionProcessor.cpp:806` 實際在 1066 → 更新 `05_data_formats/01_significance_summary_schema.md` 一處，grep 其他頁確認無重複

---

## 流程 B：結論型（positive / NEGATIVE / characterization）

```
1. 讀 docs/reports/research_landscape/XX.md 確認權威結論
2. 讀 MEMORY 對應 project_*.md 看補充資訊
3. 判斷狀態：positive / characterization / NEGATIVE / ongoing
4. 加入對應 09_conclusions/ 清單（01/02/03 之一）
5. 若對應衍生特徵存在 → 在 07_derived_features/XX.md 加一段
6. 更新 docs/reports/research_landscape/ 若需要
7. 跑 Layer 2 驗證：KB 與 docs 結論數字完全一致？
8. 跑 3 腳本
```

**範例（Phase 6 已實施）**：
- HPFineNGroups 機制更正 → `07_derived_features/01_hpfinengroups.md` 加「⚠️ 機制更正」章節；同步 `09_conclusions/01_positive_findings.md` 與 `02_characterization_only.md`

### 狀態翻轉處理
```
positive → NEGATIVE:
  1. 從 09_conclusions/01_positive_findings.md 移除
  2. 加到 09_conclusions/03_concluded_negative.md
  3. 在 07_derived_features/XX.md 加 NEGATIVE 警告區塊（不刪頁）
  4. 更新 MEMORY 對應檔案 status
  5. CHANGELOG 記錄重大變動
```

---

## 流程 C：元件型（新 pipeline / sample / feature）

```
1. 確認該元件是否該入 KB（對照 SoT 規則：已有類似元件嗎？）
2. 若新增：
   a. 複製最類似既有頁面為模板
   b. 更新 frontmatter（新 id, name, description）
   c. 更新 body（三欄必寫 + 章節內容）
   d. 在所屬 00_index.md 加一行
   e. 建立 related_ids 雙向關係
3. 若為新樣本 → 走 06_workflows/06_adding_new_sample.md 9 步驟流程
4. 跑 Layer 1-3 全部驗證
5. 跑 3 腳本
```

---

## 流程 D：狀態型（CURRENT_FOCUS / 活躍假說 / 阻塞）

```
1. 讀最新 docs/CURRENT_FOCUS.md
2. 讀最新 research/autoresearch/{hypothesis_queue,evidence_ledger}.jsonl
3. 更新 10_research_status/ 對應檔案
4. 開頭加/更新「⚠️ 此為 YYYY-MM-DD 快照，2 週有效」警示
5. last_verified 改為當日
6. 跑 3 腳本
```

**重要**：此類檔案 freshness 敏感，建議**每 2 週** `refresh_last_verified.py --file ...` 主動更新

---

## 流程 E：結構型（新目錄 / 新 schema / 新腳本）

```
1. 對照 00_governance/ 既有規範是否需改
2. 若新增目錄：
   a. 建立 NN_groupname/ 目錄 + 00_index.md
   b. 更新 README.md「目錄總覽」表
   c. 更新 00_governance/02_naming_conventions.md
3. 若改 frontmatter schema：
   a. 修 00_governance/01_frontmatter_schema.md
   b. 修 scripts/validate_frontmatter.py
   c. 所有既有檔案跑一次 migration
4. CHANGELOG 記錄架構變動
5. 跑 3 腳本 + 邀新 agent 盲測
```

---

## 🤖 Multi-Agent 驗證模板

對於**高影響**變動（B / E 類別），建議啟動 ≥2 個 agent 平行驗證：

### Agent 模板 1：File:line 事實驗證
```
啟動 general-purpose agent，prompt：
「驗證 KB 新增內容的 file:line 引用是否精確
 1. 讀 <KB 檔案>
 2. 對每個 file:line 引用，用 Read/Grep 確認
 3. 檢查參數預設值、函式位置、腳本路徑
 4. 輸出錯誤清單（按嚴重性）」
```

### Agent 模板 2：結論一致性驗證
```
啟動 general-purpose agent，prompt：
「比對 KB <檔案> 與 docs/reports/research_landscape/ 權威文件
 1. 逐項抓取 KB 中的『數字/結論』claims
 2. 到 docs/ 權威查核
 3. 特別注意方向、精度、MEMORY ID 存在性
 4. 輸出矛盾清單」
```

### Agent 模板 3：紅隊盲測
```
啟動 general-purpose agent，prompt：
「以新進者身份從 README.md 模擬查詢
 Q1-Q6: <具體 6 個查詢問題>
 記錄跳數、路徑、迷路可能
 找遺漏主題、冗餘、可執行命令失敗處」
```

**本次 v0.2 發布已跑過這 3 模板** — 參考 `CHANGELOG.md` Phase 6 紀錄。

---

## 📅 週期性維護（Recommended）

| 動作 | 週期 | 責任 |
|------|------|------|
| `refresh_last_verified.py --list-stale` | 每月 | KB 維護者 |
| `10_research_status/` 主動更新 | 每 2 週 | 研究 PI |
| `09_conclusions/05_hypothesis_queue_snapshot.md` 同步 | 每 2 週 | 研究 PI |
| 跑 3 驗證腳本 | 每次 commit 前 | 所有貢獻者 |
| Multi-agent 盲測 | 每季或重大變動後 | KB 維護者 |
| CHANGELOG 重大變動紀錄 | 每次 schema/狀態翻轉 | 貢獻者 |

---

## 🚨 違反協議的後果

| 違反 | 後果 |
|------|------|
| 跳過 SoT 直接複製數字 | 22 處重複問題（已發生於 v0.1 → v0.2 修正）|
| 不跑 3 腳本 | related_ids 不對稱、broken canonical_paths 累積 |
| 不更新 last_verified | freshness 警示無效，stale 內容被誤引用 |
| 不更新 CHANGELOG | 架構變動難以追溯 |
| 事實型變動不跑 multi-agent | 類似 big7→big8 longphase-s 路徑錯誤滲漏 |

---

## 相關

- Governance 索引：[00_index.md](00_index.md)
- 既有更新流程：[06_update_workflow.md](06_update_workflow.md)
- Freshness：[04_freshness_policy.md](04_freshness_policy.md)
- Query protocol：[05_query_protocol.md](05_query_protocol.md)
- Maintenance scripts：`../scripts/`
