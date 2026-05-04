# outline_confirm.md — P1+P2 → C1+C2

## 使用時機
P1 main thesis 鎖定 + P2 outline 結構確認。

## P1 → C1（fast-track 必停）

```yaml
question: "Main Thesis (≤30 字) + Audience tier 鎖定？"
options:
  - 接受 AI 草稿
  - 微調（用戶提供改寫）
  - 完全重寫（AI 列 4 候選 + 自訂 free-text）
  - 改 audience tier
```

## P2 → C2

```yaml
question: "5-7 段 outline，每段 thesis sentence ≤ 30 字。批准？"
options:
  - 批准 outline
  - 調整段數（4 / 6 / 8）
  - 修改順序
  - 改 main thesis（退回 P1）
```

## 從 weekly-report 母稿觸發
若 frontmatter `main_statement` 已有 → 跳過 C1，直接顯示「main thesis 從母稿讀取：<...>」一行 FYI，進 C2。
