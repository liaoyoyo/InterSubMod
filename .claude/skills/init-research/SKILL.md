---
name: init-research
description: 初始化新研究專案。建立標準目錄結構、manifest.yaml、00_PLAN.md 骨架。觸發詞：init-research、開始新研究、新研究專案。
allowed-tools: Read, Write, Bash, Glob, Grep
user-invocable: true
---

# 研究專案初始化 (Init Research)

從標準模板建立新研究專案，確保目錄結構、metadata、計劃書一致。

## 觸發時機

`/init-research {project_name}`、`開始新研究`、`新研究專案`

## 輸入參數

| 參數 | 必須 | 說明 |
|------|------|------|
| `project_name` | ✓ | 目錄名（小寫 + 底線，如 `loh_subclone_af`） |
| `title` | ✓ | 研究標題（中文或英文） |
| `hypotheses` | 建議 | 1-5 條假說摘要 |
| `datasets` | 建議 | 使用的數據集路徑 |

若用戶未提供 `title` 或 `hypotheses`，以對話方式逐一詢問。

## 執行流程

### Step 0：前置檢查（強制）

**必須先讀取 CLAUDE.md 必讀清單**，確認：
1. 讀取 `docs/CURRENT_FOCUS.md` — 確認新研究不與當前優先事項衝突
2. 讀取 `docs/experiments/INDEX.md` — 確認未與已關閉方向重複
3. 若有疑似重複，**停下來告知用戶**，提供已有結論的連結

```python
# 檢查是否與已關閉方向重複
# 搜尋 INDEX.md 中的 ❌ NEGATIVE / NO-GO 條目
# 若新假說與這些條目相似 → 警告用戶
```

### Step 1：建立目錄結構

```bash
PROJECT_DIR="research/${project_name}"

# 確認不存在
if [ -d "$PROJECT_DIR" ]; then
    echo "ERROR: $PROJECT_DIR already exists"
    exit 1
fi

# 從模板複製
cp -r research/_template/ "$PROJECT_DIR"
```

### Step 2：填寫 manifest.yaml

用用戶提供的資訊填入 `manifest.yaml`：

```yaml
project:
  name: "{project_name}"
  title: "{title}"
  created: "{YYYY-MM-DD}"
  status: "initiated"
  conclusion: null

hypotheses:
  - id: "H1"
    statement: "{hypothesis_1}"
    # ... 其他欄位留空讓執行階段填
```

### Step 3：撰寫 00_PLAN.md

從模板填入用戶提供的資訊：
- 背景與動機（從對話中提取）
- 假說（帶前提、confound、驗證標準）
- 數據來源
- 已知風險

**計劃書品質門檻**：
- 每個假說必須有明確的 positive / negative 判定標準
- 至少列出一個已知 confound
- 數據路徑必須是實際可存取的

### Step 4：更新追蹤

1. 在 `docs/CURRENT_FOCUS.md` 加入「進行中」標記（如果用戶同意）

### Step 5：確認輸出

向用戶展示建立的結構：

```
研究專案已初始化：

research/{project_name}/
├── 00_PLAN.md              ✓ 已填寫
├── manifest.yaml           ✓ 已填寫
├── scripts/                (空)
├── data/                   (空)
├── figures/                (空)
└── reports/                (空)

下一步：開始 Step 1 分析，或調用 /research-loop 進入研究迴圈。
```

## 注意事項

1. **不執行分析**：此 skill 只做 scaffold，不跑任何數據分析
2. **不修改 INDEX.md**：INDEX 在 `/conclude-research` 時才更新
3. **project_name 命名規則**：小寫英文 + 底線，不含日期（日期在 manifest 中）
4. **模板位置**：`research/_template/`（不要修改模板本身）
