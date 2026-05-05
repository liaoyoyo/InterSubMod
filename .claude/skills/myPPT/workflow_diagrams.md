# myPPT workflow diagrams

## 整體 pipeline（三 skill 互動）

```mermaid
graph TB
  Start[U: 我要做簡報/週報/教授報告] --> Trigger{myPPT 場景識別}
  Trigger -->|場景 A 週報| WR[weekly-report W1-W7]
  Trigger -->|場景 B 已有母稿| PB1[/pptx-build --from-draft/]
  Trigger -->|場景 C 從零做| PB2[/pptx-build P1-P5/]
  Trigger -->|場景 D 完整 pipeline| WR

  WR --> WC4[C4 母稿確認]
  WC4 --> Handoff{handoff 4 選}
  Handoff -->|A 立即| PB1
  Handoff -->|B 留檔| Stash[master_draft.md 存檔]
  Handoff -->|C 終點| End[週報結束]
  Handoff -->|D 加下週計畫| Plan[next_week_plan.md]

  PB1 --> P2[P2 outline]
  PB2 --> P1[P1 audience]
  P1 --> P2
  P2 --> P3[P3 section ×N]
  P3 --> P4[P4 slide ×24 + Vision 10-check]
  P4 --> P5[P5 speaker script]
  P5 --> Output[Output: PPTX + speaker + audit]

  style Trigger fill:#FFF3E0
  style Handoff fill:#FFF3E0
  style WR fill:#E8F5E9
  style PB1 fill:#E3F2FD
  style PB2 fill:#E3F2FD
  style Output fill:#FCE4EC
```

## 場景識別決策樹

```mermaid
graph TD
  Intent[用戶意圖] --> Q1{含「週報」/「整理本週」?}
  Q1 -->|是| Q1a{要產 PPT 嗎?}
  Q1a -->|是| D[場景 D 完整 pipeline]
  Q1a -->|否| A[場景 A 只產母稿]
  Q1 -->|否| Q2{提供 master_draft path?}
  Q2 -->|是| B[場景 B 已有母稿]
  Q2 -->|否| Q3{含「PPT」/「投影片」?}
  Q3 -->|是| C[場景 C 從零做]
  Q3 -->|否| Ask[AskUserQuestion 4 選]
```

## 三 skill 責任分工

```mermaid
graph LR
  subgraph myPPT [myPPT 總入口]
    M1[場景識別]
    M2[AskUserQuestion]
    M3[委派 sub-skill]
  end

  subgraph WR [weekly-report]
    W1[W1 raw data]
    W2[W2 主線]
    W3[W3 4 層分類]
    W7[W7 母稿產出]
  end

  subgraph PB [pptx-build]
    P1[P1 audience]
    P2[P2 outline]
    P4[P4 slide build]
    P5[P5 speaker]
  end

  M3 --> WR
  M3 --> PB
  WR -->|handoff A| PB

  style myPPT fill:#FFF3E0
  style WR fill:#E8F5E9
  style PB fill:#E3F2FD
```

## 委派指令對應

| 場景 | 委派 |
|------|------|
| A 週報 | 觸發 `weekly-report` skill |
| B 已有母稿 | 執行 `/pptx-build --from-draft <path>` |
| C 從零做 | 觸發 `pptx-build` skill |
| D 完整 | `weekly-report` → C4 handoff A 自動接 `/pptx-build --from-draft` |
