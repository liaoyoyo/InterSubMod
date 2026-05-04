# from_draft_loader.md — 從 weekly-report 母稿讀取

## 使用時機
觸發指令：`/pptx-build --from-draft <path>`
或自動偵測：用戶說「我要做簡報」+ 提供 master_draft.md path

## 自動讀取流程

```python
import yaml
import re

def load_master_draft(path):
    with open(path) as f:
        content = f.read()
    # 解析 frontmatter
    fm_match = re.search(r'^---\n(.+?)\n---\n', content, re.DOTALL)
    fm = yaml.safe_load(fm_match.group(1))

    # 必要欄位
    assert fm.get('type') == 'weekly_master_draft'
    assert fm.get('status') == 'ready_for_handoff'

    return {
        'main_statement': fm['main_statement'],     # 跳過 P1 main thesis 鎖定
        'report_type': fm['report_type'],            # 跳過 P1 報告類型識別
        'audience': fm.get('audience', 'advisor'),
        'target_duration_min': fm.get('target_duration_min', 25),
        'suggested_pptx_template': fm.get('suggested_pptx_template'),
        'source_artifacts': fm.get('source_artifacts', []),
        'priority_buckets': fm.get('priority_buckets', {}),
        'professor_qa_count': fm.get('professor_qa_count'),
    }
```

## 跳過項目（從母稿讀取後）

- ❌ P1 main thesis 鎖定（從 main_statement 取）
- ❌ P1 報告類型識別（從 report_type / suggested_pptx_template 取）
- ❌ §20.C Tier 評分（從母稿 [F]/[O]/[I]/[U] 自動 mapping：F→S, O→A, I→B, U→C）

## 仍保留項目

- ✅ P2 Outline checkpoint（拆 5-7 段）
- ✅ P3 Section checkpoint
- ✅ P4 Slide checkpoint（含 §20.E 6 問 audit + visual_review.md 10-check）
- ✅ P5 Speaker script

## 母稿 → slide 萃取

| 母稿欄位 | slide 用途 |
|---------|----------|
| Layer 0.1 + §1 + §2 | thesis_title_bar + cover slide |
| Layer 2 Thread × N | 1-2 張 slide per Thread |
| §3 [F] / §4 [O][I] / §5 [U] | slide 內容 + 4 層分類標籤帶過 |
| Layer 3 整合 | 1 張 integration slide |
| Layer 4 §16 | future tree slide |
| Layer 4 §17 教授追問 | backup slides / Q&A slides |

## Failure modes

| 情境 | 處置 |
|------|------|
| 母稿不存在 | 退回直接觸發模式（從零開始 P1）|
| frontmatter status != ready_for_handoff | 警告 + 要求先回 weekly-report C4 確認 |
| report_type 不在 6 模板列表 | 預設 improvement_report + 提示用戶調整 |
| 母稿 frontmatter 不完整 | 標記缺欄，退回 P1 補填 |
