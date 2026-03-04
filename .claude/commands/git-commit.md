---
allowed-tools: Bash(git:*)
argument-hint: [type] [message]
description: 提交變更並附上 Co-Author
---

提交變更到目前分支。

## 使用方式
```
/git-commit feat "add bernoulli distance calculation"
/git-commit fix "correct methylation parsing logic"
/git-commit docs "update API documentation"
```

## 執行步驟

1. 檢查當前狀態
```bash
git status
```

2. 暫存變更
```bash
git add -A
```

3. 提交變更
```bash
git commit -m "$1: $2

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

## 提交類型
| 類型 | 說明 |
|------|------|
| feat | 新功能 |
| fix | 錯誤修復 |
| docs | 文件更新 |
| refactor | 重構 |
| test | 測試相關 |
| chore | 維護工作 |
