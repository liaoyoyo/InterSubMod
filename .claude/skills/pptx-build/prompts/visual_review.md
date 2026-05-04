# visual_review.md — Claude Vision 10-Checkpoint

## 使用時機
每張 slide build 後立即跑（C4 前）。

## 10-Check 表格

| # | 項目 | PASS / FAIL / PARTIAL | 理由 |
|:-:|------|---------------------|------|
| 1 | Heading = thesis sentence | | |
| 2 | 1 idea/slide | | |
| 3 | ≤ 6 elements | | |
| 4 | 文字密度（中 ≤ 60 / 英 ≤ 30） | | |
| 5 | Distracted takeaway | | |
| 6 | 視覺對比（顏色/字型/字級） | | |
| 7 | 圖表 vs 文字 ≥ 60% 視覺 | | |
| 8 | 引文 / 數據來源 | | |
| 9 | 1 minute timing check | | |
| 10 | 技術災難備援 | | |

## Pass 標準
- ≥ 8 PASS = 通過
- 1-2 FAIL = 修正後再 review
- ≥ 3 FAIL = 重新設計 slide

## 從 PNG 讀取
```python
from ppt_toolkit.claude_vision_review import review_slide_png
result = review_slide_png(png_path, focal_point="...", section_thesis="...")
print(result['pass_rate'], result['issues'])
```

## P5 整體 timing alarm
P5 結束時加總所有 slide 的 1 minute check：
- 24 張 × 1 min = 24 min vs 預設 30 min → OK
- 字數 24,412 → 60 min vs 30 min → 警告，建議拆 Tier 3
