---
id: ism-kb-00-governance-frontmatter-schema
name: "Frontmatter Schema 規範"
description: "KB 所有 .md 文件的 YAML frontmatter 欄位定義、取值範圍、使用時機；對齊 /big8_disk/knowledge/ 格式。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: reference
verified_scope: "frontmatter schema specification"
related_ids:
  - ism-kb-00-governance-naming-conventions
  - ism-kb-00-governance-cross-reference-rules
  - ism-kb-00-governance-freshness-policy
  - ism-kb-00-governance-query-protocol
  - ism-kb-00-governance-update-workflow
tags: [governance, frontmatter, schema, yaml, metadata]
canonical_paths: [00_governance/01_frontmatter_schema.md]
alias_paths: []
---

# Frontmatter Schema 規範

- 一句結論：KB 所有 `.md` 文件的 YAML frontmatter 必須遵循本 schema；`alias_paths` 僅放路徑變體，`tags` 承擔主題搜尋
- 適用對象：文件作者、`scripts/validate_frontmatter.py` 維護者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  python3 /big7_disk/liaoyoyo2001/InterSubMod/knowledge/scripts/validate_frontmatter.py \
    /big7_disk/liaoyoyo2001/InterSubMod/knowledge/
  ```

---

## 欄位總表

| 欄位 | 必要 | 值型 | 用途 | 範例 |
|------|------|------|------|------|
| `id` | ✅ | string (kebab-case) | 全庫唯一識別符 | `ism-kb-03-pipelines-paired-full` |
| `name` | ✅ | string | 人類可讀標題 | `"Paired Full Pipeline"` |
| `description` | ✅ | string (60-120 字) | 搜尋結果摘要 | 見範例 |
| `status` | ✅ | `active` / `archived` / `deprecated` | 文件狀態 | `active` |
| `last_verified` | ✅ | `YYYY-MM-DD` | 最後人工驗證日期 | `2026-04-22` |
| `content_nature` | ✅ | 見下表（7 種） | 內容性質分類（v0.5+；對齊 Dublin Core `type`） | `runtime-fact` |
| `doc_type` | ✅ | `reference` / `howto` / `explanation` / `tutorial` | Diataxis 分類 | `reference` |
| `verified_scope` | ✅ | string | `last_verified` 當日核對了什麼 | `"CLI args against include/utils/ArgParser.hpp on HEAD"` |
| `decision_refs` | ⚠️ | list[string] | 若引用專案決策，列 decision ID | `[]` |
| `related_ids` | ⚠️ | list[string] | 關聯其他 KB 文件 `id`（**雙向對稱**） | `[ism-kb-04-parameters-cli-arguments]` |
| `tags` | ✅ | list[string] | 主題搜尋關鍵字 | `[pipeline, paired, benchmark]` |
| `canonical_paths` | ✅ | list[string] | 本文件正式路徑（通常 1 個） | `[03_pipelines/01_paired_full.md]` |
| `alias_paths` | ⚠️ | list[string] | **路徑變體**（不可放主題詞） | `[]` |

✅ 必要、⚠️ 有條件必要

---

## `content_nature` 取值（7 種；v0.5 由 `source_type` rename）

| 值 | 定義 | 使用時機 |
|----|------|----------|
| `runtime-fact` | 當下 code/tool/sample 的現況事實 | CLI 參數、樣本路徑、執行命令 |
| `frozen-decision` | 已凍結的專案決策或架構選擇 | governance 規範、pipeline 定義 |
| `reference-summary` | 索引、摘要、導航類（不可直接執行） | 00_index、結論索引 |
| `reference` | 一般參考文件（結構規範、benchmark 協議等，介於 summary 與 runtime 之間） | benchmark 協議、檔案格式規範 |
| `paper-derived` | 來自論文或外部文獻 | — |
| `historical-note` | 歷史脈絡，已過期不可當現行依據 | 已淘汰的 approach、過期快照 |
| `postmortem` | 事件 / 資料陷阱 / 研究失敗的結構化記錄（Google SRE blameless postmortem 格式） | pitfalls、NEGATIVE 研究回顧 |

**業界對照**：本 schema 對齊 Google SRE postmortem culture（`postmortem`）+ Diataxis（`reference-summary` vs `reference`）+ Dublin Core Type 元素。詳見下方「Dublin Core 對照」。

---

## `doc_type` 取值（Diataxis 分類）

| 值 | 特徵 | 範例 |
|----|------|------|
| `reference` | 查詢型，表格為主，精確 | 參數表、欄位字典 |
| `howto` | 解題型，步驟清楚，可執行 | workflow、benchmark 流程 |
| `explanation` | 理解型，概念解說 | 五大目標、突破策略 |
| `tutorial` | 學習型，新手 end-to-end | （本 KB 較少用） |

---

## `id` 命名規則

格式：`ism-kb-{group}-{topic}`

- `ism-kb-` — 固定前綴（InterSubMod Knowledge Base）
- `{group}` — 頂層目錄兩位數字 + 單詞（例：`00-governance`、`03-pipelines`）
- `{topic}` — 文件主題（kebab-case）

**範例**：
- `ism-kb-03-pipelines-paired-full`（`03_pipelines/01_paired_full.md`）
- `ism-kb-04-parameters-cli-arguments`
- `ism-kb-09-conclusions-concluded-negative`

**驗證**：全庫唯一，`scripts/validate_frontmatter.py` 會檢查重複。

---

## `related_ids` 雙向對稱原則

A 文件 `related_ids` 列了 B，則 B 的 `related_ids` 也必須列 A。

**自動檢查**：`scripts/check_related_ids_symmetry.py`

**何時適用**：
- 兩個文件互為補充（pipeline → parameters）
- 一個文件是另一個的展開（conclusions index → 具體 feature）

**何時不需雙向**：
- 索引文件（00_index）指向子文件：不建雙向，因為 index 會指向太多
- 連結到 **非 KB** 的外部檔案：不用 `related_ids`（用 markdown link）

---

## 範例（正確）

```yaml
---
id: ism-kb-03-pipelines-paired-full
name: "Paired Full Pipeline"
description: "ClairS paired full VCF + LongPhase-S haplotag + ISM 的 canonical benchmark pipeline；delta F1=+0.0112 locked。包含輸入鏈路、參數、輸出結構、典型命令。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "pipeline command against scripts/run_vcf_all_snv.sh HEAD"
decision_refs: []
related_ids:
  - ism-kb-04-parameters-cli-arguments
  - ism-kb-04-parameters-distance-metrics
  - ism-kb-06-workflows-full-vcf-analysis
  - ism-kb-08-truth-and-benchmark-f1-calculation
tags: [pipeline, paired, benchmark, canonical, clairs]
canonical_paths: [03_pipelines/01_paired_full.md]
alias_paths: []
---
```

---

## 📚 業界標準對照（Dublin Core / Diataxis / PROV-O / ADR）

本 schema 為專案自訂但**語義對齊**以下業界標準：

| KB 欄位 | Dublin Core (ISO 15836) | PROV-O | ADR (Nygard) | Diataxis |
|---------|-------------------------|--------|--------------|----------|
| `id` | identifier | — | — | — |
| `name` | title | — | title | — |
| `description` | description | — | — | — |
| `tags` | subject | — | — | — |
| `doc_type` | type | — | — | **✅ 四象限分類** |
| `last_verified` | date (valid) | — | — | — |
| `content_nature`（v0.5 rename 自 `source_type`）| type（延伸：內容性質分類）| Entity type | — | — |
| `related_ids` | relation | wasDerivedFrom | — | — |
| `canonical_paths` / `alias_paths` | source / hasPart | — | — | — |
| `verified_scope` | coverage | wasGeneratedBy | — | — |
| `decision_refs` | — | — | **supersedes / amends** | — |
| `entities` (v0.2+，可選) | — | **✅ Entity** | — | — |

### 缺欄位（未來可選補強）
- `creator` — Dublin Core 核心；本 KB 為單一團隊產出暫省略
- `rights` / `license` — 未來開源時需補
- `language` — 本 KB 繁中+英混合，暫以專案慣例處理

### 命名不衝突聲明（v0.5 rename 理由）
- **v0.1-v0.4 原名 `source_type`** 與 Dublin Core `source`（本資源衍生自哪個原始源）**語義衝突** — Dublin Core 的 `source` 語義較接近本 KB 的 `canonical_paths` + `verified_scope`
- **v0.5 rename 為 `content_nature`** — 明確為「內容性質分類」，對齊 Dublin Core `type` 延伸概念，避免誤解
- KB 的 `doc_type`（Diataxis 分類）與 Dublin Core `type` 重合 — 一個專用於文件類型（Diataxis 四象限），一個專用於資源內容性質（runtime-fact / frozen-decision / etc.）

---

## 反模式（禁止）

```yaml
# ❌ 錯：id 沒有 ism-kb- 前綴
id: 03-pipelines-paired-full

# ❌ 錯：alias_paths 放主題詞（應在 tags）
alias_paths: [canonical, paired-normal]

# ❌ 錯：content_nature 不在允許列表
content_nature: experimental

# ❌ 錯：使用已 deprecated 的舊欄位名 source_type（v0.5 已改為 content_nature）
source_type: runtime-fact

# ❌ 錯：last_verified 格式錯誤
last_verified: 2026/04/22

# ❌ 錯：related_ids 指向不存在的 id
related_ids: [ism-kb-99-fake-doc]
```

---

**相關**：
- 命名規則：[02_naming_conventions.md](02_naming_conventions.md)
- 雙向對稱檢查：[03_cross_reference_rules.md](03_cross_reference_rules.md)
- Freshness：[04_freshness_policy.md](04_freshness_policy.md)
- 驗證腳本：[../scripts/validate_frontmatter.py](../scripts/validate_frontmatter.py)
