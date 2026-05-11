# 來源資料說明（v2 PPT source materials）

本資料夾保留產生本次 v2 PPT 的三份來源 .md 報告，供 fact-check 與深度閱讀用。

## 1. 三份來源報告

| 檔案 | 對應原始位置 | 用途 |
|------|------------|------|
| `01_V5_IGV_session_visual_audit.md` | `InterSubMod/docs/reports/pi_reports/2026/04/20260427_V5_IGV_session_visual_audit_01.md` | 15 位點 × 4 BAM 並列 IGV 真截圖、HP1↔HP2 翻轉、self-phasing extreme 視覺證據 |
| `02_Self_Phasing_Baseline_V5_Audit.md` | `InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md` | 一句結論審核版，HPFineNGroups 重詮釋、舊說法錯誤訂正、caveat 清單 |
| `03_longphase_TO_vs_V5_技術報告.md` | `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` | structured-tech-report 13 段格式，4-commit 演進、F1 變動因果鏈、5 層獨立證據鏈 |

## 2. 三份報告的角色定位

```
報告 03（技術報告）  ── 主敘事骨架（13 段結構部分對應到 PPT 26 slides）
   ├─ 引用報告 01（IGV 視覺）    ── PPT 視覺證據主源（fig01d、IGV 4-BAM）
   └─ 引用報告 02（審核）         ── PPT 「修正後說法 vs 舊說法」段落主源
```

## 3. v2 內其他配套文件（不在本子資料夾）

從 v2/ 角度導覽，以下檔案在**同層**（不在 `source_materials/`）：

| 檔案路徑（從 `v2/` 起算）| 角色 |
|------|------|
| `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/00_storyboard_v2.md` | **26-slide 主敘事大綱**（必讀；對應 v2 PPT 結構）|
| `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/notes/00_background_context.md` | 背景濃縮（含 v5_audit_suite 6 子報告）|
| `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/notes/glossary.md` | 30+ 術語對應 slide 表 |
| `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/notes/key_metrics_table.md` | 重點數據速查 + 來源行號 + fact-check 警示 |
| `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/notes/code_references.md` | commit hash + 行號 + 重現指令 |
| `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/notes/industry_references.md` | 業界對照（LongPhase-S 2025 同實驗室 paired 版）|
| `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/notes/qa_11_questions.md` | 教授 11 條關鍵「為什麼」+ 預備答案 |

## 4. 不能直接從 .md 自動轉的內容

1. **演講者口語化解釋**：報告為文字密集，PPT 需大量視覺＋少字標題；逐字稿將於 PPT build 階段產出於 `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/notes/speaker_script_v2.md`（**目前尚未生成**）
2. **逐 slide 圖片配置**：報告為線性段落、PPT 為並列佈局；對應關係在 `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/00_storyboard_v2.md`
3. **背景知識補充**：教授未必熟 LongPhase-TO 內部結構，PPT 開頭需 4-5 張定義 slides；本資料夾的 `notes/00_background_context.md` 已濃縮供 onboarding

## 5. v2 當前完成狀態

| 階段 | 狀態 |
|------|------|
| Storyboard 26 slides | ✅ 完成（`00_storyboard_v2.md`）|
| Source materials 實體複製 | ✅ 完成（本資料夾）|
| Figures 33 張實體複製 | ✅ 完成（`figures/`）|
| 6 份 onboarding 文件 | ✅ 完成（`notes/`）|
| `scripts/build_pptx.py` | ⏭ **待 build 階段** |
| `scripts/screenshot_all.py` | ⏭ **待 build 階段** |
| `output/*.pptx` | ⏭ **待 build 階段** |
| `notes/speaker_script_v2.md` | ⏭ **待 build 階段** |

> **v2 目前定位為 storyboard + source package**，PPTX 與 speaker script 將於 build 階段生成。
