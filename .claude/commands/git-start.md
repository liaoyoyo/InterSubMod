---
allowed-tools: Bash(git:*)
argument-hint: [branch-type] [description]
description: 開始新任務，建立功能分支
---

開始新的開發任務。

## 使用方式
```
/git-start feature add-new-distance
/git-start fix methylation-parsing-bug
/git-start docs update-readme
```

## 執行步驟

1. 確認當前在 develop 分支
```bash
git checkout develop
git pull origin develop
```

2. 建立新分支
```bash
git checkout -b $1/$2
```

## 分支類型
- `feature` - 新功能
- `fix` - 錯誤修復
- `refactor` - 重構
- `docs` - 文件更新
