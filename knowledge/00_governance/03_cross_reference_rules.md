---
id: ism-kb-00-governance-cross-reference-rules
name: "交叉引用規則"
description: "KB 內 related_ids 雙向對稱、連結 docs/research/code 的標準；何時用 markdown link 何時用 frontmatter。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: reference
verified_scope: "cross-reference policy"
related_ids:
  - ism-kb-00-governance-frontmatter-schema
  - ism-kb-00-governance-query-protocol
tags: [governance, cross-reference, related-ids, linking]
canonical_paths: [00_governance/03_cross_reference_rules.md]
alias_paths: []
---

# 交叉引用規則

- 一句結論：KB 內部用 `related_ids`（雙向）+ markdown link；KB 到 docs/research/code 只用 markdown link（相對路徑）
- 適用對象：文件作者、`scripts/check_related_ids_symmetry.py` 維護者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  python3 /big7_disk/liaoyoyo2001/InterSubMod/knowledge/scripts/check_related_ids_symmetry.py \
    /big7_disk/liaoyoyo2001/InterSubMod/knowledge/
  ```

---

## 三種連結類型

### 1. KB 內文件之間：用 `related_ids` + markdown link

**frontmatter**：
```yaml
related_ids:
  - ism-kb-04-parameters-cli-arguments
  - ism-kb-08-truth-and-benchmark-f1-calculation
```

**body 內連結**：
```markdown
完整 CLI 參數見 [04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md)
```

**雙向對稱要求**：A 列 B 則 B 必列 A。用 `check_related_ids_symmetry.py` 驗證。

### 2. KB 到 docs/：只用 markdown link（相對路徑）

**理由**：docs/ 非 KB 成員，不登記 frontmatter；使用相對路徑便於 git 追蹤。

**範例**：
```markdown
權威研究報告見 [../docs/reports/research_landscape/09_Part_B.md](../../docs/reports/research_landscape/09_Part_B.md)

當前進度快照見 [docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
```

**路徑起點**：從目前 KB 文件位置往上跳。
- `knowledge/03_pipelines/01_paired_full.md` → `../../docs/...`
- `knowledge/09_conclusions/00_index.md` → `../../docs/...`

### 3. KB 到 source code：`file:line` 文字引用（不用 hyperlink）

**理由**：行號會隨 git HEAD 變動；寫 `file:line` 便於讀者用編輯器跳轉。

**範例**：
```markdown
`--min-mapq` 預設 20，定義於 `include/utils/ArgParser.hpp:67-68`

距離計算核心邏輯：`src/core/DistanceMatrix.cpp:20-45`（NHD）
```

---

## `related_ids` 使用時機

✅ **應建立 `related_ids`**：
- 兩文件互為補充（pipeline ↔ parameters）
- 一文件是另一文件的展開（feature ↔ conclusion）
- 讀者看完 A 後「下一步自然會想看 B」

❌ **不應建立 `related_ids`**：
- 索引文件（`00_index.md`）指向大量子文件（會變噪音）
- 連結只出現在 body 一次且不重要
- 連結目標不在 KB 內

---

## 索引文件（`00_index.md`）例外

索引文件通常**不設 `related_ids`**（或只留最重要的 2-3 個），因為：
- 索引本來就列所有子文件
- 若 index 列所有子，則所有子都要回指 index → 污染 graph

替代做法：在 index body 用 markdown link 列出子文件，簡單清楚。

---

## 連結到專案既有關鍵文件的標準路徑

以下為高頻引用的外部權威文件，**應一致使用這些路徑**（從 KB 子目錄算起，`../../` 跳到專案根）：

| 外部文件 | 從 KB 的相對路徑 |
|----------|------------------|
| `docs/CURRENT_FOCUS.md` | `../../docs/CURRENT_FOCUS.md` |
| `docs/experiments/INDEX.md` | `../../docs/experiments/INDEX.md` |
| `docs/reports/research_landscape/00_INDEX.md` | `../../docs/reports/research_landscape/00_INDEX.md` |
| `docs/concepts/2026/04/20260409_研究構想總索引_01.md` | `../../docs/concepts/2026/04/20260409_研究構想總索引_01.md` |
| `docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md` | `../../docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md` |
| `.claude/CLAUDE.md` | `../../.claude/CLAUDE.md` |
| `research/autoresearch/hypothesis_queue.json` | `../../research/autoresearch/hypothesis_queue.json` |
| `research/autoresearch/evidence_ledger.jsonl` | `../../research/autoresearch/evidence_ledger.jsonl` |

---

## 反模式（禁止）

```markdown
# ❌ 錯：硬連結絕對路徑（無法在不同機器運作）
見 [/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md]

# ❌ 錯：假連結（目標不存在）
見 [完整分析](../docs/imaginary.md)

# ❌ 錯：列在 related_ids 但 body 沒解釋為何相關
related_ids: [ism-kb-99-random]
```

---

## 驗證流程

1. **frontmatter 合規**：`validate_frontmatter.py`
2. **related_ids 雙向**：`check_related_ids_symmetry.py`
3. **canonical_paths 存在**：`check_canonical_paths.py`
4. **Markdown link 有效性**：手動或 `grep` 抓連結路徑檢查存在性

每次提交 KB 變動前，**至少跑前三項**。

---

**相關**：
- Frontmatter：[01_frontmatter_schema.md](01_frontmatter_schema.md)
- 查詢流程：[05_query_protocol.md](05_query_protocol.md)
