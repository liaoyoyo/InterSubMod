---
id: ism-kb-00-governance-update-workflow
name: "更新工作流程"
description: "新增/修改/淘汰 KB 文件的標準流程；包含前置檢查、變更、驗證、提交。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: howto
verified_scope: "KB update workflow"
related_ids:
  - ism-kb-00-governance-frontmatter-schema
  - ism-kb-00-governance-freshness-policy
  - ism-kb-00-governance-new-info-protocol
tags: [governance, workflow, update, maintenance]
canonical_paths: [00_governance/06_update_workflow.md]
alias_paths: []
---

# 更新工作流程

- 一句結論：KB 變動必經「前置檢查 → 執行變更 → 跑 3 腳本驗證 → 更新 CHANGELOG → git commit」五步
- 適用對象：所有 KB 內容貢獻者（人類或 agent）
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  cd /big7_disk/liaoyoyo2001/InterSubMod/knowledge
  python3 scripts/validate_frontmatter.py . && \
  python3 scripts/check_related_ids_symmetry.py . && \
  python3 scripts/check_canonical_paths.py .
  ```

---

## 三種變動情境

### 情境 1：新增文件

```
1. 確認目錄分組
   → 若不確定，參考 00_governance/02_naming_conventions.md「邊界與歧義處理」

2. 決定檔名
   → NN_topic_keyword.md（見 02_naming_conventions.md）
   → 找該目錄下現有最大 NN，+1

3. 寫 frontmatter
   → 複製 01_frontmatter_schema.md 範例
   → id 為 ism-kb-<group>-<topic>
   → last_verified 設今日

4. 寫 body
   → 開頭三欄：一句結論 / 適用對象 / 可直接執行命令
   → 依 doc_type 寫對應章節（見 01_frontmatter_schema.md）

5. 建立 related_ids 雙向關係
   → 在 frontmatter 列出最相關的 2-4 個既有文件
   → 去那些文件也加回指

6. 更新所屬目錄 00_index.md
   → 新增一行（檔名、一句結論、何時讀）

7. 驗證（跑 3 個腳本）

8. 更新 CHANGELOG.md 與 README.md 導航表（若為重要文件）

9. git commit
```

### 情境 2：修改既有文件

```
1. 讀現有 frontmatter 與 body，確認要改什麼

2. 執行變更

3. 更新 frontmatter
   → last_verified 改今日
   → verified_scope 描述這次核對了什麼
   → 若 body 主題有變，考慮更新 description 與 tags

4. 檢查 related_ids 是否仍合理
   → 新增、刪除的相關連結要雙向同步

5. 驗證（跑 3 個腳本）

6. 若為重大修訂（架構變、結論反轉等）→ 更新 CHANGELOG.md

7. git commit
```

### 情境 3：淘汰文件

```
1. 確認淘汰原因
   ├── 已被新文件取代 → 改 status: deprecated，description 說明替代
   ├── 整個主題不再相關 → status: archived
   └── 從未發佈（誤建） → 可直接 git rm

2. 不刪除檔案（保留 git history）
   → status: deprecated 或 archived
   → description 開頭加 [DEPRECATED] 或 [ARCHIVED]

3. 清理 related_ids
   → 其他文件若指向此 deprecated 文件，改指向替代文件或刪除 related_ids

4. 更新所屬 00_index.md
   → 加 [DEPRECATED] 標註或移至附錄

5. 驗證（跑 3 個腳本）

6. 更新 CHANGELOG.md

7. git commit
```

---

## 驗證腳本（全部必跑）

### 1. frontmatter 合規
```bash
python3 knowledge/scripts/validate_frontmatter.py knowledge/
```
**檢查項**：欄位完整性、值型正確、id 唯一、`source_type`/`doc_type`/`status` 取值合法

### 2. related_ids 雙向對稱
```bash
python3 knowledge/scripts/check_related_ids_symmetry.py knowledge/
```
**檢查項**：A 列 B 則 B 必列 A

### 3. canonical_paths 存在性
```bash
python3 knowledge/scripts/check_canonical_paths.py knowledge/
```
**檢查項**：每份文件 `canonical_paths` 指向的檔案實存

### 4.（可選）markdown link 有效性
```bash
# 抓 docs/ 外部連結檢查
grep -rohE '\]\(\.\./\.\./docs/[^)]+\)' knowledge/ | \
  sed 's/.*(\(.*\))/\1/' | sort -u | \
  while read p; do test -e "knowledge/$p" || echo "Broken: $p"; done
```

---

## 提交前檢查清單

在 `git commit` 前對照：

- [ ] Frontmatter 完整、`last_verified` 已更新
- [ ] Body 含「一句結論 / 適用對象 / 可直接執行命令」三欄
- [ ] 所有 `related_ids` 雙向對稱
- [ ] Markdown link 路徑正確
- [ ] 所屬 `00_index.md` 已同步更新
- [ ] 3 個驗證腳本通過
- [ ] 重大變動已加入 `CHANGELOG.md`

---

## 特殊變動類型

### 研究結論狀態轉換（positive ↔ NEGATIVE）
- 涉及 `09_conclusions/` 多份文件的互相移動
- **一定要**同步更新 `docs/reports/research_landscape/` 的指向
- 建議在 commit message 明確說明：`docs(kb): move HPFineNGroups from positive to characterization`

### 新樣本加入
- 需加：`02_samples/NN_<sample>.md`、`08_truth_and_benchmark/01_truth_set_registry.md` 新增一列
- 需檢查：`02_samples/00_index.md` 的表格、`README.md` 若提及 7 樣本需更新為 8

### 新 pipeline 加入（罕見）
- 整個 `03_pipelines/` 結構可能需重組
- 影響 `04_pipeline_comparison.md`、典型命令速查
- 建議先開 RFC 討論再動手

---

## Commit message 建議格式

```
docs(kb): <變動摘要>

- <具體變動 1>
- <具體變動 2>

Updated last_verified for: <列出的文件>
Scripts validated: validate_frontmatter, check_related_ids_symmetry, check_canonical_paths
```

**範例**：
```
docs(kb): add 03_pipelines/01_paired_full.md

- New pipeline reference document with canonical command and output structure
- Updated 03_pipelines/00_index.md
- Added bidirectional related_ids to 04_parameters/01_cli_arguments.md

Updated last_verified for: 03_pipelines/00_index.md, 04_parameters/01_cli_arguments.md
Scripts validated: all 3 passed
```

---

**相關**：
- Frontmatter：[01_frontmatter_schema.md](01_frontmatter_schema.md)
- 命名：[02_naming_conventions.md](02_naming_conventions.md)
- Freshness：[04_freshness_policy.md](04_freshness_policy.md)
- 驗證腳本：[../scripts/](../scripts/)
