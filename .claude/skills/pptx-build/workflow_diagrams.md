# pptx-build workflow diagrams

## 整體流程（5 階段 + checkpoint）

```mermaid
graph TB
  Start[U: 我要做簡報] --> Trigger{From draft?}
  Trigger -->|是 --from-draft| Loader[from_draft_loader.md<br/>讀 frontmatter]
  Trigger -->|否| P1[P1 Audience & Goal]
  Loader --> P2[P2 Outline]
  P1 --> C1[C1 ★ main thesis ≤ 30 字 必停]
  C1 --> P2
  P2 --> C2[C2 5-7 段 outline]
  C2 --> P3[P3 Section batch ×N]
  P3 --> C3[C3 5-7 slide + focal point]
  C3 --> P4[P4 Slide build ×24]
  P4 --> C4[C4 三層分流 + 10-check]
  C4 --> P5[P5 Speaker script]
  P5 --> C5[C5 ★ 字數↔時長 必停]
  C5 --> Output[Output: PPTX + speaker + audit]
```

## §20 6 階段過濾

```mermaid
graph LR
  A[A. Main Thesis 鎖定] --> B[B. Slide Focal Point]
  B --> C[C. Tier S/A/B/C/D 評分]
  C --> D[D. Def/Prereq/Body/Conclusion 比例]
  D --> E[E. 6 問 audit]
  E --> F[F. 雜訊紅旗清單]
  F -->|通過| Build[slide build]
  F -->|失敗| Discard[棄用]
```

## 6 模板識別流程

```
觸發 keyword 偵測
   ├─ 「修補/修復/fix」 → improvement_report
   ├─ 「A vs B/對照」 → comparison_report
   ├─ 「月會/PI/KPI」 → executive_summary
   ├─ 「結果/實驗」 → data_showcase
   ├─ 「教學/解釋」 → concept_walkthrough
   └─ 「教授/thesis」 → academic_defense

混合用例（如 v3 = improvement + academic）
→ 取主敘事用模板 + sub template 影響特定段
```

## ppt_toolkit 調用流程

```
用戶 → SKILL.md trigger
     → playbook §1-§24 規範
     → style_library YAML（colors / objects / layouts）
     → ppt_toolkit Python helpers
     → python-pptx → PPTX
     → screenshot_all.py → wireframe PNG
     → Claude Vision 10-checkpoint review
     → 反饋 → 修正
```

## --from-draft handoff 流程

```mermaid
graph LR
  WR[weekly-report W7] --> C4[C4 母稿確認]
  C4 --> H[handoff 4 選]
  H -->|A 立即| Trigger[/pptx-build --from-draft path/]
  H -->|B 留檔| Stash[master_draft.md 存檔]
  H -->|C 終點| End[週報結束]
  H -->|D 加下週計畫| Plan[next_week_plan.md]
  Trigger --> Loader[from_draft_loader.md]
  Loader --> SkipP1[跳過 P1 main thesis 鎖定]
  SkipP1 --> P2[P2 Outline]
```
