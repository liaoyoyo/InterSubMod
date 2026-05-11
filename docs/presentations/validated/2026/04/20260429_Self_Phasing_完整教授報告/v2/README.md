# Self-Phasing 完整教授報告 · v2 版本

> 26-slide 演講者模式 PPT · 主軸聚焦版 · 2026-04-29
>
> **目前定位：v2 storyboard + source package**（PPTX、speaker script、wireframe 等待 build 階段產生；見下方第 6 節「執行指令」與第 7 節「v2 完成狀態」）。
>
> **本資料夾自包含程度**：來源報告 (3 份)、所有圖片 (33 張)、背景概念 (含 v5_audit_suite 6 子報告濃縮)、重點數據 (含來源行號) 都已實體複製進來；單獨拿走 `v2/` 即可閱讀整體論證鏈、重現背景理解、做 fact-check。**完整 PPTX 重現需 build 階段執行 `scripts/build_pptx.py`**（待產生）。

## 1. 為何有 v2

v1（同層 `../`）採線性敘事 34 slides。經 4-agent review（科學敘事 / 教授認知 / 數據準確性 / PPT 視覺設計）後重新組織為 **26-slide 6 段流程**：高層次重點 → 觀察問題 → 為何重要 → 解釋與原因 → 改動驗證 → 未來規劃。

關鍵差異：S1 衝擊式 TL;DR、S2 業界家族樹 framing、S11/S12 兩張 prerequisite、S13 baseline 三 bug 根因樹合併、S18 V5max1 設為 climax、S26 take-home + next step 收尾。

## 2. 資料夾結構（自包含）

```
v2/
├── README.md                          ← 本檔（入口）
├── 00_storyboard_v2.md                ← 26-slide 完整大綱（必讀）
├── source_materials/                  ← 三份來源報告（實體複製）
│   ├── README.md                      角色說明
│   ├── 01_V5_IGV_session_visual_audit.md       (PI 報告，IGV 真截圖 15 sites × 4-BAM)
│   ├── 02_Self_Phasing_Baseline_V5_Audit.md    (validated, 一句結論審核版)
│   └── 03_longphase_TO_vs_V5_技術報告.md       (主敘事骨架, 13 段結構)
├── figures/                           ← 33 張 PNG（實體複製）
│   ├── HP tag / phasing 概念圖（fig2/3/4/9/17）
│   ├── code diff schematic（fig01a-d）
│   ├── self-phasing impact 3-tier（03_self_phasing_impact）
│   ├── IGV 4-BAM 真截圖（D_SP1/2/3, C_V5max1/2/3, A_TP01）
│   └── 量化結果（fig18 concordance/AMB/F1）
├── notes/                             ← 6 份輔助索引文件
│   ├── 00_background_context.md       背景概念（含 v5_audit_suite 6 子報告濃縮）
│   ├── glossary.md                    術語表（30+ 術語對應 slide）
│   ├── key_metrics_table.md           重點數據速查表（含來源行號）
│   ├── code_references.md             commit hash + 檔案行號對照
│   ├── industry_references.md         業界對照文獻（已寫，含 LongPhase-S 2025）
│   └── qa_11_questions.md             教授 11 條關鍵「為什麼」預備答案
├── scripts/                           ← PPT 建構腳本（dispatch agent 後產生）
│   ├── build_pptx.py                  ⏭ 待寫
│   ├── screenshot_all.py              ⏭ 待寫
│   └── PPT_RENDERING_RULES.md         ⏭ 待寫
└── output/                            ← 產出
    ├── 20260429_..._v2.pptx           ⏭ 待生成
    └── wireframes/                    ⏭ 待生成
```

## 3. 閱讀順序建議

| 角色 | 閱讀順序 |
|------|---------|
| **快速 onboarding（5 min）** | README → `notes/00_background_context.md` §1 一句總結 → `notes/key_metrics_table.md` 看數字 |
| **深度準備（30 min）** | README → `00_storyboard_v2.md` → `notes/00_background_context.md` 全文 → `notes/qa_11_questions.md` |
| **演講者準備（60 min）** | 上面 + `source_materials/03_longphase_TO_vs_V5_技術報告.md` 全文 + `notes/glossary.md` 默讀 |
| **重現 PPT（執行）** | `scripts/build_pptx.py` → `scripts/screenshot_all.py` |
| **fact-check 數字** | `notes/key_metrics_table.md` 對照 `source_materials/` |
| **改寫程式碼** | `notes/code_references.md` 對照 `/big7_disk/liaoyoyo2001/longphase-to-mod/` |

## 4. 核心訊息（一句話 + 分層精確化）

> **Self-phasing 問題鏈被分層處理**：V2b（commit `8b8c1fd`）解 **Phase scaffold**（`--pon-only-phasing` flag）；V3F（commit `41ff147`）解 **Tag 層 priority bug + enum mismatch**；INDEL guard（`380e8d2`）補 UB；V5（working tree）加 **Layer 1.5 somatic fallback**（HP33 directional reassignment）。InterSubMod 是下游消費者本 repo 無 C++ 改動；V5 真實價值在 **read-level concordance +8.3 pp（clean PS blocks，全基因組）**，SEQC2 F1 不變是預期行為（V5 不改 caller、不改 phasing graph，只改 BAM HP tag）。

> **重要區分**：V5 **不修 self-phasing 本身**（self-phasing 是 phasing graph 層問題、由 V2b 處理）；V5 解的是 V3F 留下的「過於保守把所有平手 read 標 HP33」問題，把 directional 線索從 somatic vote 取回。

## 5. 關鍵數字速查（詳見 `notes/key_metrics_table.md`）

| 維度 | 數值 |
|------|------|
| Baseline HP1 family : HP2 family | **17.3 : 1**（HCC1395 全基因組 614K vs 35.5K, 94.6% 集中於 HP1）|
| AMB% | 17.5% → **8.0%**（−9.5 pp）|
| HP:i:33 reads | 239,679 → 110,197（**−54%**）|
| Clean PS paired GT concordance | Baseline 82.2% → V5 **90.5%**（+8.3 pp）|
| 15-site clean PS concordance | Baseline 74.9% → V5 **88.2%**（+13.3 pp）|
| ClairS-TO raw F1 | **0.7166（三版本完全相同 — V5 不改 caller）**|
| V5 vs Baseline ΔF1（套 ISM 後）| **−0.0003（噪音級）**|
| Sanity check | **15/15 PASS, 0 violation** |
| 跨樣本一致性 | 7/7 樣本同方向，Cohen's d = −1.20，r = 0.001 |
| ISM 影響分類（共 85 features）| 🔴 嚴重 **29** 個 / 🟡 中度 **14** 個 / 🟢 不受影響 **42** 個 *（百分比口徑來源報告誤寫為 38%/7%/55%；以 count 為準）* |

## 6. 執行指令（待 builder agent 完成腳本後）

```bash
cd InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2
python3 scripts/build_pptx.py        # 重建 PPTX
python3 scripts/screenshot_all.py    # 重做 wireframe + 結構檢查
```

## 7. 與 v1 的關係

| 維度 | v1（同層 `../`）| v2（本資料夾）|
|------|---------------|--------------|
| Slides | 34 | 26 |
| 結構 | 線性 | 6 段（結論先行）|
| TL;DR 位置 | S29 F1 釐清 | S1 衝擊式 |
| 業界對照 | 無 | S2 + S21 兩張 |
| 五大目標銜接 | 無 | S23 |
| 程式碼 diff | 無 | S15 一張 |
| Prerequisite | 無 | S11 函數 + S12 PON |
| Climax 設計 | 無 | S18 V5max1 |
| 結尾 | TL;DR | take-home + next step 承諾 |

v1 保留作為備份（不要刪），可作 fallback；v2 為當前推薦版本。
