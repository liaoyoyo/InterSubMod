---
id: ism-kb-00-governance-naming-conventions
name: "命名規則"
description: "KB 內檔名、目錄、id、canonical_path 的命名規則；與專案 /doc-standards 的差異點說明。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: reference
verified_scope: "naming conventions"
related_ids:
  - ism-kb-00-governance-frontmatter-schema
tags: [governance, naming, conventions]
canonical_paths: [00_governance/02_naming_conventions.md]
alias_paths: []
---

# 命名規則

- 一句結論：KB 內一律用 `NN_lowercase_underscore.md`，id 用 `ism-kb-<group>-<topic>`；與專案 `/doc-standards` 的中文檔名不同
- 適用對象：新增 KB 文件者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 檢查命名合規
  find /big7_disk/liaoyoyo2001/InterSubMod/knowledge -name "*.md" | \
    grep -vE '/(0[0-9]|README|AGENT|CHANGELOG)' | head -10
  ```

---

## 目錄命名

```
NN_group_name/           # 兩位數字 + 底線 + 小寫英文
```

**範例**：
- ✅ `00_governance/`、`03_pipelines/`、`08_truth_and_benchmark/`
- ❌ `governance/`（缺編號）、`03-pipelines/`（用連字號）、`03_Pipelines/`（大寫）

**目錄編號意義**：
- `00_` — 後設/規範（governance）
- `01_` - `10_` — 主題分組，按讀取優先級排序
- 頂層 `99_` 保留給 `scripts/` 類非內容目錄

---

## 檔案命名

```
NN_topic_keyword.md       # 兩位數字 + 底線 + 小寫英文關鍵字
```

**範例**：
- ✅ `01_paired_full.md`、`02_distance_metrics.md`、`03_concluded_negative.md`
- ❌ `paired-full.md`（用連字號）、`01_PairedFull.md`（駝峰）、`1_paired.md`（數字未補零）

**特殊檔案**：
- `00_index.md` — 每個目錄的索引（**必有**）
- `README.md`、`AGENT.md`、`CHANGELOG.md` — 頂層特殊檔案，無編號

---

## `id` 命名

格式：`ism-kb-{group-slug}-{topic-slug}`

- 前綴固定：`ism-kb-`
- `{group-slug}`：兩位數字 + 連字號 + 小寫英文（對應目錄）
- `{topic-slug}`：主題關鍵字（連字號分隔，小寫）

**範例**：
| 檔案路徑 | id |
|----------|-----|
| `00_governance/01_frontmatter_schema.md` | `ism-kb-00-governance-frontmatter-schema` |
| `03_pipelines/01_paired_full.md` | `ism-kb-03-pipelines-paired-full` |
| `04_parameters/02_distance_metrics.md` | `ism-kb-04-parameters-distance-metrics` |
| `09_conclusions/03_concluded_negative.md` | `ism-kb-09-conclusions-concluded-negative` |
| `02_samples/01_hcc1395.md` | `ism-kb-02-samples-hcc1395` |

**特殊**：
- `00_index.md` 的 id 為 `ism-kb-{group-slug}-index`（例：`ism-kb-03-pipelines-index`）
- 頂層特殊檔案無需 frontmatter（README、AGENT、CHANGELOG）

---

## `canonical_paths` 命名

- 一律**相對於 `knowledge/` 根目錄**
- 單一文件通常只有 1 個 canonical path
- 範例：`[03_pipelines/01_paired_full.md]`

---

## `alias_paths` 命名

**只放路徑變體**，不放主題詞（這是 D9 決策，對齊既有 `/big8_disk/knowledge/`）。

✅ 應放入：
- 舊路徑（檔案搬移前）
- 大小寫變體（`HCC1395.md` vs `hcc1395.md`）
- 底線 vs 連字號變體

❌ 不應放入（改放 `tags`）：
- 主題關鍵字
- 中文主題詞
- 工具功能描述

---

## 與專案 `/doc-standards` skill 的差異

| 面向 | KB 規範 | `/doc-standards`（for `docs/`） |
|------|---------|--------------------------------|
| 檔名語言 | 英文 snake_case | 繁中描述性 |
| 日期前綴 | 無（用 last_verified） | `YYYYMMDD_主題_NN.md` |
| 目錄編號 | `NN_` 前綴（固定） | `YYYY/MM/` 年月分層 |
| Frontmatter | **必要** | 不強制 |

**原因**：KB 文件較少（~43），靜態查詢；docs/ 文件多（~208），時序性強。

---

## 邊界與歧義處理

### 需要跨分組的文件
選最相關的**單一分組**，用 `related_ids` 或 tags 連結其他分組。
- 範例：Zone-Aware Framework 放 `07_derived_features/03_zone_aware_framework.md`，用 `related_ids` 指向 `09_conclusions/`。

### 多個 topic 混合的文件
若文件含多個 topic 且難以拆分，用複合名稱：
- 範例：`08_truth_and_benchmark/` 而非 `08_truth/` + `09_benchmark/`

### 淘汰文件
不刪除，改 `status: deprecated` 並更新 `description` 說明替代文件：
```yaml
status: deprecated
description: "[DEPRECATED] 已被 03_pipelines/01_paired_full.md 取代"
```

---

**相關**：
- Frontmatter schema：[01_frontmatter_schema.md](01_frontmatter_schema.md)
- 更新流程：[06_update_workflow.md](06_update_workflow.md)
