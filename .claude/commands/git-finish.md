---
allowed-tools: Bash(git:*)
description: 完成任務，合併分支回 develop
---

完成當前任務，將分支合併回 develop。

## 執行步驟

1. 確認當前分支狀態
```bash
git status
git branch
```

2. 切換到 develop 並更新
```bash
git checkout develop
git pull origin develop
```

3. 合併當前分支
```bash
git merge --no-ff <current-branch>
```

4. 刪除已合併的分支
```bash
git branch -d <current-branch>
```

## 注意事項
- 合併前確保所有測試通過
- 確認無衝突
- 合併後驗證 develop 分支狀態
