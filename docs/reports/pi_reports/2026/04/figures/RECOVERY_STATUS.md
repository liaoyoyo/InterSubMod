# pi_reports figures 復原狀態（2026-05-11）

此目錄含 IGV 截圖系列（igv_v5_audit/by_HP/, by_chr/）已從 disk 消失。

## 復原策略

IGV 截圖無法由 .py 自動重產 — 需要：
1. 開 IGV
2. 載入對應 BAM (HCC1395 / COLO829 等)
3. 跳到對應 chr/position
4. 截圖

部分 IGV 截圖在備份位置可用：
- `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/figures/`
  - `A_TP01_chr6_145444893.png` ✓
  - `C_V5max1_chr19_4639528.png`, `C_V5max2`, `C_V5max3` ✓
  - `D_SP1_chr19_17565944.png`, `D_SP2`, `D_SP3` ✓

## 自動重產 script

`InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/scripts/screenshot_all.py`
可批次重產 IGV 截圖（如果 IGV + BAM 環境還在）。
