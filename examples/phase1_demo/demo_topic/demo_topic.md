---
title: "Demo topic — Phase 1 closed-loop validation"
date: 2026-05-10
status: phase1_demo
authors: [Claude Opus 4.7 (assistant)]
spec: InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md
---

# Demo Topic — Phase 1 Closed-Loop Validation

> 此 .md + 同名主題資料夾用於驗證 Phase 1 三件事：
> 1. image-gen skill 能跑通 AI 軌 (concept + data) 與 cairo 軌 (flow + icon)
> 2. image-vision-check 能對 4 張圖各打 6 維分數
> 3. 手寫極簡 index.html 能 inline 4 張 base64 PNG 並正確顯示
> 此檔不是 Phase 2 html-preview skill 的設計範例，僅閉環驗證用。

## 範例圖示需求 (4 類各 1 張)

<!-- figure-needed: concept_diagram, slug=fig1_concept_loh -->
**Figure 1**: LOH 區段如何引發 phasing 失效的概念示意。

<!-- figure-needed: flow_diagram, slug=fig2_flow_pipeline -->
**Figure 2**: Self-Phasing V5 G1-G3 簡化 pipeline 流程圖。

<!-- figure-needed: data_mockup, slug=fig3_data_mockup -->
**Figure 3**: 7 樣本 baseline AUC 分布的概念型 bar chart (mock-up，非真數據)。

<!-- figure-needed: icon, slug=fig4_icon_status -->
**Figure 4**: status_pass 綠色勾勾 icon。

## 預期結果

執行：
```bash
# 1. 生圖
python3 InterSubMod/.claude/skills/image-gen/tools/dispatch.py \
  InterSubMod/examples/phase1_demo/demo_topic/prompts \
  InterSubMod/examples/phase1_demo/demo_topic/figures

# 2. 視覺檢核 (instruction phase)
python3 InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py \
  InterSubMod/examples/phase1_demo/demo_topic/figures

# (Claude agent 處理每張圖、寫 _pending_entries/*.json)

# 3. 視覺檢核 (aggregate)
python3 InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py \
  InterSubMod/examples/phase1_demo/demo_topic/figures --aggregate

# 4. 開瀏覽器
xdg-open InterSubMod/examples/phase1_demo/demo_topic/index.html
```

開瀏覽器後預期：4 張圖 inline 顯示、每張下方有 vision check verdict + score。
