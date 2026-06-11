---
name: doc-standards
description: InterSubMod 文檔管理規範。檔案命名格式、元數據模板、圖片規則、多步驟目錄結構、AI session 報告規則。觸發：建立新 .md 檔案或重組文件目錄時。DO NOT USE WHEN：寫 SKILL.md（用 plugin-dev:skill-development）、寫 13 段技術報告（用 structured-tech-report）、寫週報（用 weekly-report）。
allowed-tools: Read, Write, Glob, Grep
user-invocable: true
paths: ["docs/**/*.md", "research/**/*.md", "**/*.md"]
---

# InterSubMod 文檔管理規範

## 文檔目錄結構

```
docs/
├── architecture/        # 專案主軸架構說明
├── concepts/            # 構想紀錄（研究大圖景、發展樹、理論基礎）
├── plans/               # 計劃文件（YYYY/MM 分層）
├── solutions/           # 問題解決報告（任務目標/YYYY/MM 分層）
├── experiments/         # 實驗紀錄（in_progress / validated / finalized）
├── reports/             # 正式說明文件（validated / finalized，索引用 00_INDEX.md）
│   └── {topic}/         # 多檔案說明文件
│       └── figures/     # 說明文件圖表
├── methodology/         # 方法學審查與 closeout 報告
├── standards/           # 命名規範、狀態管理、output 規範等治理文件
├── decisions/           # 重整與治理決策紀錄
├── presentations/       # 簡報 PDF
├── provenance/          # AI 對話 provenance
│   └── ai_sessions/     # AI 對話紀錄（YYYY/MM 分層）
├── data_specs/          # 數據規格
├── references/          # 參考資料（含 manual/、external/）
├── research/            # 長期研究主題工作區
│   └── {study_name}/    # 單一研究主題
│       ├── README.md    # 研究概述與索引
│       ├── figures/     # 圖表（統一命名，gitignore）
│       ├── data/        # 中間數據（gitignore）
│       ├── scripts/     # 分析腳本（git 追蹤）
│       └── reports/     # 研究報告
├── refactor_baseline/   # 重構基準數據
├── archive/             # 歸檔
│   └── deep/            # 深度歸檔（需查詢歷史）
└── trash/               # 暫存待刪除
```

## 檔案命名格式

```
{YYYYMMDD}_{中文說明目標}_{流水號}.md
```

範例：`20260111_文檔庫重整計劃_01.md`

## 多步驟研究專案目錄結構

當研究/任務包含多個 Step 或子計劃時，在對應目錄下建立**專案子資料夾**：

```
docs/plans/YYYY/MM/{YYYYMMDD}_{專案主題}/
  ├── 00_總覽與執行順序.md        # 索引：波次、優先序、依賴關係
  ├── Step1_{子計劃主題}.md
  ├── Step2_{子計劃主題}.md
  └── ...

docs/architecture/{YYYYMMDD}_{專案主題}/
  ├── {YYYYMMDD}_{架構說明}_{流水號}.md
  └── {YYYYMMDD}_{資料追蹤表}_{流水號}.md
```

**規則**：
- 資料夾名稱格式：`{YYYYMMDD}_{中文專案主題}`
- 索引檔用 `00_` 開頭
- Step 計劃用 `StepN_` 前綴（方便排序）
- 同一研究的架構文件、追蹤表集中在同一專案子目錄
- 可按日期或任務主題搜尋定位

## 檔案元數據

每個 Markdown 檔案開頭需包含：

```markdown
<!--
建立時間: YYYY-MM-DD HH:MM
目標: [本檔案的目標或用途]
處理範圍: [涵蓋的工作範圍]
關聯檔案:
  - [相關檔案路徑 1]
  - [相關檔案路徑 2]
-->
```

## ⭐ 2026-06-12 新增 — Provenance Stamp + 口徑欄位（改進 ①⑦；落地盤點稽核 P-17）

> **背景**：2026-06-12 盤點稽核發現「憑記憶寫 status、不標所在 branch」是準確度+可尋性的共同根因（5 worktree 並行讓 cross-branch 幻覺更嚴重）。對應 known-pitfalls **P-17** + §13.7。

**(A) 盤點 / audit / status / 現況整理類報告**：metadata block **必加 provenance stamp**（助手 `bash scripts/provenance_stamp.sh` 一行貼上）：
```markdown
<!--
建立時間: YYYY-MM-DD
build_branch: <git rev-parse --abbrev-ref HEAD>      # 本報告的 current 相對哪個 branch
build_commit: <git rev-parse --short HEAD>
worktree: <pwd>                                       # 5 worktree 並行時定位用
data_sources: <path>,<path>                           # 數字/事實來源檔（number_provenance_check.sh gate 用）
驗證方式: <每個 status=current 的事實旁附驗證指令；未驗證標 unverified>
-->
```

**(B) 涉及不同 build / 口徑的資料**（Tmode / Paired / TO、read-instance / unique、max-collapse / 5mC-only）：metadata 加 `build_mode:` + `build_date:` + `data_version:`，避免口徑混用（pitfalls P-12）。

**(C) status 斷言鐵則（P-17）**：任何「某檔存在 / current / stale / N 個」斷言旁**必附驗證指令 + 一行輸出**；無附 → 降 `unverified`，禁寫 `current`。盤點「可尋性」卻不去 ls = 自我矛盾。

## ⭐ 2026-05-26 新增 — Partial-Scope Ribbon 強制規則（A5 落地，5/24 incident postmortem）

**規則**：當文件涵蓋的 scope < 100%（如 3/24 chr, 1/7 sample, single-cycle subset），**必須在以下 3 個位置同步標註 partial flag**：

### 1. 元數據區（HTML comment 開頭）

```markdown
<!--
建立時間: YYYY-MM-DD HH:MM
目標: ...
處理範圍: ⚠ PARTIAL SCOPE — 3/24 chr (chr1, chr8, chr19) × 1/7 sample (HCC1395)
完整 scope 未驗證: chr2-7, chr9-23, chrX, chrY × 其他 6 樣本
補強計劃: cycle N+1 擴展至全 24 chr × 7 樣本 (預估 +30-50 min)
-->
```

### 2. 文件首段 hero / TL;DR 顯著位置

```markdown
> ⚠ **PARTIAL SCOPE**: 本報告僅涵蓋 3/24 chr × 1/7 sample。完整 scope 待 cycle N+1 驗證。
```

### 3. HTML standalone footer / .md 結尾 caveat 段

```markdown
## Scope Limitation
- 已驗證: chr1 (5,000 sites), chr8 (3,200 sites), chr19 (2,100 sites) × HCC1395
- 未驗證: 剩餘 21 chr + 6 樣本（包含 chr16, chrX phasing-weak outlier 風險）
- 推論限制: 結論僅適用已驗證 scope；不可外推到 full genome 直到 cycle N+1 完成
- 補強動作: 預定 YYYY-MM-DD 跑全量
```

### Task Type 對應 partial flag 義務（per AGENTS.md §15.3）

| Task Type | partial flag 義務 |
|-----------|------------------|
| (A) Exploratory pilot | 必標 partial flag + 後續補強計劃 |
| **(B) Comprehensive validation** | **不應 partial** — 若 partial 必明示為何（時間 / 計算成本）+ 補強 deadline |
| **(C) Production deployment** | **不應 partial** — partial 視為 NOT release-ready |
| **(D) Handoff to external** | **不應 partial** — 必先補完全 scope 才能對外 |
| (E) Hotfix | 最小可重現 OK，但必加 regression test scope 説明 |
| (F) Demo | 必首段標 `[DEMO]` + 不可作 validation evidence |

**5/24 incident reference**: HKU handoff 任務（task type D）PI .md / HTML / dashboard 均無 partial flag → 5/26 修正為 mandatory 3-locations。

---

## 圖片存放與引用規則

- 純圖片子目錄統一命名 `figures/`（不用 `images/`、`plots/`）；混合類型資源目錄可用 `assets/`
- .md 引用圖片用相對路徑：`figures/xxx.png`
- 相對路徑最多 2 層（禁止 `../../../`）；引用 `output/` 或 `research/` 圖片時允許超過 2 層
- 圖片命名：`{NN}_{英文描述}.png`

## 嚴謹度繼承（/scientific-rigor）

文件元數據規範（檔案命名 / metadata / 目錄結構）+ **嚴謹度敘述格式** 雙軌。本 skill 加入 `InterSubMod/.claude/skills/scientific-rigor/SKILL.md` 對照：

- **§2 Evidence Tier ribbon 標註位置**:
  - 文件 metadata block 加 `證據等級: L1/L2/L3/L4/L5 ⭐...`
  - 主要 claim 標題 + 內文用 `⭐⭐⭐⭐ L2` inline 標
- **§3 Effect Size 敘述標準**:
  - 數字 metric 必含 Cohen ribbon: 「+0.0112 (Cohen's d 0.08, marginal)」
  - CI 95% 必標: 「n=5, 95% CI [+0.003, +0.020]」
- **§10.2 confidence 詞彙紅線**: 禁止「clearly」「strongly」「rigorously」「significant」（無 p-value 對齊）— 用具體 evidence 敘述替代
- **§8.4 Provenance 引用**: 報告必含 commit hash / cycle_id / 數據版本（standalone 報告 footer 強制）

**檔案命名 + tier 結合範例**:
- `20260518_HP_priority_bug_L1_完全佐證_01.md`（tier 嵌入 filename，便於 grep）
- 或 metadata block：`<!-- 證據等級: L1 ⭐⭐⭐⭐⭐ / Cohen: d=2.3 large effect -->`

## AI 對話紀錄撰寫

每次 AI 對話完成重要任務後，撰寫執行報告（可使用 `/report` skill）：

1. **報告位置**：`docs/provenance/ai_sessions/YYYY/MM/`
2. **檔案格式**：`{YYYYMMDD}_{對話主題}_{流水號}.md`
3. **必要內容**：
   - 對話目標
   - 關鍵決策
   - 修改的檔案清單
   - 後續行動建議
