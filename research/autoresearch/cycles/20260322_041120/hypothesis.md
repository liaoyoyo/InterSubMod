假設 ID: H004_modified (新注入)
假設內容: Quality_Score >= 50 作為獨立 rescue 規則（無需 GQ gate）
  - 移除對 CramersV 的依賴（TO 候選中 CramersV 中位數 = 0，無判別力）
  - 降低 QS 門檻從 60 到 50（更多 TP 被救回）
  - 移除 GQ 前置條件（讓 QS 訊號獨立運作）
來源: OBSERVE 分析發現 QS>=50 alone 預測 delta_F1 = +0.008556 vs GQ+QS>=60 的 +0.005551
Pipeline Track: TO (Tumor-Only)
Scale: S1 pilot (HCC1395_5kHz_TO)
修改層級: Tier 1 (Python 分析腳本)
預期方向: 正增益，預測 delta_F1 = +0.008556

OBSERVE 關鍵發現:
1. HPP 在 TO 中無判別力（TP/FP 分佈幾乎相同，各約 52% HPP<0.05）
   → H001/H002 預測 FAIL，本輪跳過
2. CramersV 在 TO 候選中中位數 = 0，QS+CramersV 組合僅捕獲 42/556 TP
   → H004 原版 (QS+CramersV) 無效
3. QS>=50 alone（無 GQ gate）預測 delta_F1 = +0.008556
   → 優於已知最佳 GQ+QS>=60 的 +0.005551
4. CramersV>0.3 在 TO: 48 candidates, 41 TP (85.4% precision) — 範圍太小
