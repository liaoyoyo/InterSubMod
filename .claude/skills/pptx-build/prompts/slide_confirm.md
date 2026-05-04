# slide_confirm.md — P4 → C4 (×24 slide)

## 使用時機
每張 slide build + Vision 10-check 完成後逐張確認。

## AskUserQuestion

```yaml
question: "Slide <N>: <focal point>。三層分流 + 10-check 結果如下。批准？"
options:
  - 批准 slide（進下張）
  - 補 Tier 1 element
  - 削減 Tier 2（拆 Tier 3）
  - 棄用 slide（§20.E 第 6 問 fail）
```

## 顯示內容
- Tier 1 (slide) 字數 / element count
- Tier 2 (speaker note) 字數 / 預估 sec
- Tier 3 [ORAL-OPTIONAL] 標數
- Vision 10-check 表格（PASS/FAIL/PARTIAL × 10）
- §20.E 6 問 audit 結果
- §20.F 紅旗清單
