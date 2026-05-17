---
name: memory-consolidation
description: 記憶檔案生命週期管理。掃描 status/last_relevant 標記、合併同主題 concluded、降級過時 active、更新 MEMORY.md 索引（保持 <200 行）。USE WHEN：「整理記憶」「memory consolidation」「記憶太多」「MEMORY.md 太長」。SKIP WHEN MEMORY.md <100 行不需整理、純 code edit、單一 memory 新增（直接 Write）、純 build / commit。
allowed-tools: Read, Write, Edit, Glob, Grep, Bash
user-invocable: true
---

# Memory Consolidation

整理 memory 目錄維持索引精簡、防止重複、合併相關結論。

## 觸發時機

- 記憶檔案 ≥50 個
- MEMORY.md ≥180 行（接近 200 上限）
- 用戶明確要求

## 路徑

讀取 MEMORY.md 取得當前 project 的 memory 目錄。本專案路徑為：
`/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/`

## 執行流程

### Step 1：掃描現狀 [FYI]

```bash
MEM=/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory
wc -l $MEM/MEMORY.md
ls $MEM | grep -v MEMORY.md | wc -l
```

```
Grep(pattern="^status: active",    path=MEM, glob="*.md", output_mode="count")
Grep(pattern="^status: concluded", path=MEM, glob="*.md", output_mode="count")
Grep(pattern="^status: archived",  path=MEM, glob="*.md", output_mode="count")
```

輸出：`N active / M concluded / K archived，MEMORY.md L 行`。

### Step 2：識別可合併項目 [Review]

#### 2a. 同主題多筆 concluded（≥3 筆相同前綴/關鍵詞）

```
Grep(pattern="^name:", path=MEM, glob="project_O*.md")
```

例：O11/O12/O13 都是甲基化 NEGATIVE → 合併為 `project_methylation_features_exhaustion_summary.md`。

#### 2b. 過時 active（last_relevant > 30 天）

```bash
TODAY=$(date +%s)
for f in $MEM/*.md; do
  [ "$(basename $f)" = "MEMORY.md" ] && continue
  grep -l "^status: active" $f > /dev/null || continue
  LR=$(grep "^last_relevant:" $f | head -1 | sed 's/last_relevant: *//')
  [ -z "$LR" ] && continue
  AGE=$(( ($TODAY - $(date -d "$LR" +%s)) / 86400 ))
  [ $AGE -gt 30 ] && echo "$AGE 天 - $(basename $f)"
done
```

#### 2c. 重複內容（description 高度相似）

對照 MEMORY.md 各條目的 hook，>80% 字面重疊 → 候選合併。

展示候選清單給用戶確認。

### Step 3：執行合併 [Gate]

用戶確認後：

1. **建立合併檔案**：寫入 `project_<topic>_summary.md`，frontmatter 含 `status: concluded`、`merged_from: [原始檔名清單]`
2. **原始降級**：被合併的檔案 frontmatter 改為 `status: archived`，保留檔案不刪
3. **更新 MEMORY.md**：
   - Active Research / Pending Items：保留 active 條目
   - Methodology & Preferences：保留
   - Concluded：用合併後的單行條目取代多條原始條目
   - 確保總行數 < 200

### Step 4：驗證 [FYI]

```bash
MEM=/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory
for f in $MEM/*.md; do
  base=$(basename $f)
  [ "$base" = "MEMORY.md" ] && continue
  status=$(grep "^status:" $f | head -1 | sed 's/status: *//')
  [ "$status" = "archived" ] && continue
  grep -q "$base" $MEM/MEMORY.md || echo "MISSING IN INDEX: $base"
done
wc -l $MEM/MEMORY.md
```

## 合併規則

| 條件 | 動作 |
|------|------|
| ≥3 筆同主題 NEGATIVE/NO-GO | 合併為 1 筆 summary |
| active + last_relevant >30 天 | 建議降級 concluded |
| concluded + 已被新記憶覆蓋 | 降級 archived |
| archived ≥20 | 提醒可清理最舊（仍不刪檔） |

## Status 與 MEMORY.md 對應

| status | MEMORY.md 區段 | 呈現 |
|--------|-------------|------|
| `active` | Active Research / Pending Items | 完整單行 hook |
| `concluded` | Concluded | 單行精簡 |
| `archived` | （不列入） | — |

## 鐵律

1. **永不刪除原始記憶檔案** — 只改 status
2. **合併檔案必須列出 merged_from** — 保留 provenance
3. **MEMORY.md 是索引** — 每條 <150 字元，總長 <200 行
4. **AI 提合併建議，用戶決定** — 不可自主執行 step 3

---

## 與 /scientific-rigor 元方法論的關係

本 skill 支援 `/scientific-rigor §10.1 Spaced Repetition` 的**結論週期回想機制**:
- Concluded / NEGATIVE 條目的 last_relevant / status 標記 = `/scientific-rigor §10.1` 「7d / 30d / 90d 後 spaced check」的觸發依據
- 合併 + 降級邏輯（active → archived）對應 `/scientific-rigor §8.4 知識追溯 audit` 的「來自哪份數據 / 哪個 cycle」溯源（provenance 由本 skill 的 `merged_from` 提供）
- MEMORY.md 索引維持 <200 行對應 `/scientific-rigor §0.5 最小可用子集` 的 Cognitive Load 限制原則

**Step 2c — 理解追溯檢核**（建議在 Step 2 識別可合併項時新增）: 掃描 `project_*.md` / `feedback_*.md`，檢查：
- ☐ 包含「為什麼這結論可信」的證據鏈來源？
- ☐ 是否記錄「人類閱讀 AI 分析後修正過的假設」？
- ☐ concluded 檔是否清楚說明「下次如何應用此發現」（呼應 `/scientific-rigor §8.3 Reflexion buffer`）？

**級聯觸發**: `/scientific-rigor §10.2 程式修改前 retrieval practice` → Read MEMORY.md → 本 skill 維持索引健康度
