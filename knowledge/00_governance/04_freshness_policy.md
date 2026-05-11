---
id: ism-kb-00-governance-freshness-policy
name: "Freshness Policy"
description: "last_verified 時間戳管理、freshness 分級（verified/needs-recheck/stale）、失效重新驗證流程。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: reference
verified_scope: "freshness policy specification"
related_ids:
  - ism-kb-00-governance-frontmatter-schema
  - ism-kb-00-governance-update-workflow
  - ism-kb-00-governance-new-info-protocol
tags: [governance, freshness, last-verified, staleness]
canonical_paths: [00_governance/04_freshness_policy.md]
alias_paths: []
---

# Freshness Policy

- 一句結論：`last_verified` 決定 freshness 分級；`verified`(<90d) / `needs-recheck`(90-180d) / `stale`(>180d)；過期文件需先驗證再使用
- 適用對象：KB 讀者（判斷是否可信）、文件作者（知道何時需更新）
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 列出 >90 天未驗證的 KB 文件
  python3 /big7_disk/liaoyoyo2001/InterSubMod/knowledge/scripts/refresh_last_verified.py \
    /big7_disk/liaoyoyo2001/InterSubMod/knowledge/ --list-stale
  ```

---

## Freshness 三級

| 分級 | 條件 | 使用建議 |
|------|------|----------|
| `verified` | `today - last_verified < 90 天` | 可直接引用、信任 |
| `needs-recheck` | `90 ≤ today - last_verified < 180 天` | **執行前**先跑文件內命令驗證 |
| `stale` | `today - last_verified ≥ 180 天` | **禁止**直接引用；必須先重新驗證並更新 |

**自動計算**：由 `scripts/refresh_last_verified.py --list-stale` 產生報告。

**缺值處理**：若 `last_verified` 欄位缺失或格式錯誤 → 一律視為 `needs-recheck`。

---

## `source_type` 對 freshness 的影響

不同 `source_type` 的 freshness 敏感度不同：

| source_type | 預期 staleness | 說明 |
|-------------|---------------|------|
| `runtime-fact` | **高敏感** | code/path/command 會隨 git 變動；過期風險高 |
| `frozen-decision` | **低敏感** | 決策定案不變；過期無妨（但仍應 verify 決策還在有效） |
| `reference-summary` | **中敏感** | 索引對象可能新增/淘汰 |
| `paper-derived` | **低敏感** | 論文結論不變 |
| `historical-note` | **N/A** | 本身即為歷史，不需 freshness |

---

## 重新驗證流程

### 觸發條件
1. 自動腳本報告 `stale`（>180 天）
2. 文件作者被動發現內容錯誤
3. 相關 `docs/` 或 source code 有重大變動
4. 新一輪研究結果推翻或更新既有結論

### 執行步驟
1. **讀文件開頭「可直接執行命令」** — 若命令仍能跑且結果合理 → 更新 `last_verified`
2. **若命令過期或失敗** — 需實質修訂內容：
   - 比對 `verified_scope` 列出的來源（code/doc/data）最新版
   - 修改文件內容
   - 更新 `last_verified` 為今日
   - 更新 `verified_scope` 描述這次核對了什麼
   - 在 `CHANGELOG.md` 留記錄（若為重大變動）
3. **若文件已過時且無法修復** — 改 `status: deprecated`，在 `description` 說明替代文件

---

## 時間敏感文件的特殊處理

以下目錄內容變動快，需要更短 freshness 週期：

| 目錄 | 建議重驗週期 | 原因 |
|------|--------------|------|
| `10_research_status/` | **2 週** | CURRENT_FOCUS、active 假說頻繁變動 |
| `09_conclusions/05_hypothesis_queue_snapshot.md` | **2 週** | 假說佇列快照 |
| `07_derived_features/` | 4-8 週 | feature 結論可能更新 |
| `03_pipelines/` | 3-6 個月 | pipeline 架構較穩定 |
| `04_parameters/` | 3-6 個月 | CLI 參數異動較慢 |
| `00_governance/` | 6 個月+ | 規範應穩定 |

**高頻文件**應在 body 首段加警示：
```markdown
> ⚠️ **此為快照，有效期 2 週**。最新狀態請見 [docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)。
```

---

## `verified_scope` 使用指引

`verified_scope` 一句話描述「上次驗證時核對了什麼」，便於後續驗證時參照。

**好範例**：
- `"CLI args against include/utils/ArgParser.hpp on HEAD (commit 4dc2d73)"`
- `"truth VCF paths verified via ls on /big8_disk/data/"`
- `"F1 formula against scripts/pipeline/steps/03_filter_analysis.py:229-234"`

**壞範例**：
- `"verified"`（無資訊）
- `"checked everything"`（不精確）

---

## 實務範例

### 範例 1：參數文件 90 天後重驗
```bash
# 讀文件內命令
./build/bin/inter_sub_mod --help 2>&1 | head -40

# 若輸出與文件表格一致 → 只需更新 last_verified
# 若有參數新增/變動 → 修訂表格 + 更新 verified_scope
```

### 範例 2：研究狀態快照過期
```bash
# 比對最新 CURRENT_FOCUS
diff <(head -100 docs/CURRENT_FOCUS.md) \
     <(head -100 knowledge/10_research_status/01_current_focus_snapshot.md)

# 若差異大 → 重寫 snapshot；若小 → 更新 last_verified
```

---

**相關**：
- Frontmatter：[01_frontmatter_schema.md](01_frontmatter_schema.md)
- 更新流程：[06_update_workflow.md](06_update_workflow.md)
- 腳本：[../scripts/refresh_last_verified.py](../scripts/refresh_last_verified.py)
