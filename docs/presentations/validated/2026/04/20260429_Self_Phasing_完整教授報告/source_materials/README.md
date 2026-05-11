# 來源資料說明（PPT source materials）

本資料夾保留產生本次 PPT 的三份來源 .md 報告，供反查用。

| 檔案 | 對應原始位置 | 用途 |
|------|------------|------|
| `01_V5_IGV_session_visual_audit.md` | `InterSubMod/docs/reports/pi_reports/2026/04/20260427_V5_IGV_session_visual_audit_01.md` | 15 位點 × 4 BAM 並列 IGV 真截圖、HP1↔HP2 翻轉、self-phasing extreme 視覺證據 |
| `02_Self_Phasing_Baseline_V5_Audit.md` | `InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md` | 一句結論審核版，HPFineNGroups 重詮釋、舊說法錯誤訂正、caveat 清單 |
| `03_longphase_TO_vs_V5_技術報告.md` | `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` | structured-tech-report 13 段格式，4-commit 演進、F1 變動因果鏈、5 層獨立證據鏈 |

## 三份報告的角色定位

```
報告 03（技術報告）  ── 主敘事骨架（段落結構照搬到 PPT）
   ├─ 引用報告 01（IGV 視覺）    ── PPT 視覺證據主源（fig01d、IGV 4-BAM）
   └─ 引用報告 02（審核）         ── PPT 「修正後說法 vs 舊說法」段落主源
```

## 不能直接從 .md 自動轉的內容

1. **演講者口語化解釋**：報告為文字密集，PPT 需大量視覺＋少字標題；逐字稿在 `notes/speaker_script.md`
2. **逐 slide 圖片配置**：報告為線性段落、PPT 為並列佈局；對應關係在 `00_storyboard.md`
3. **背景知識補充**：教授未必熟 LongPhase-TO 內部結構，PPT 開頭需 4-5 張定義 slides，原報告假設讀者已熟悉
