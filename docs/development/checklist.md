# 待確認項目 (Checklist)

## High Priority: 架構與邏輯決策 (Blockers)

- [ ] **上游 VCF 的 INDEL/Multi-allelic 過濾**
  - **問題**: 系統目前假設 Biallelic SNV，若讀入 INDEL 會導致 Array Index 錯亂。
  - **建議解法**: 在 `S0 (LoadVCF)` 階段嚴格檢查 `REF` 和 `ALT` 的長度，若 `length > 1` 或 `ALT` 含逗號，直接跳過並 Log Warning。

- [ ] **Phase Set (PS) 跨越處理**
  - **問題**: 2000bp window 可能跨越 Phase Block，導致 HP 定義不一致。
  - **建議解法**: 第一版先實作「主流決策 (Majority Vote)」：若 Window 內 90% reads 屬於同一 PS，則分析；否則標記 `AMBIGUOUS_PHASE` 並跳過，避免過度工程化。

- [ ] **聚類演算法的選擇**
  - **問題**: 需確認 C++ Library 支援度。
  - **建議解法**: 優先使用 `fastcluster` (C++ interface)，它支援距離矩陣輸入且效能極佳 ($O(N^2)$)。若依賴過重，可考慮抽取其原始碼入專案。

## Medium Priority: 實作細節 (Implementation)

- [ ] **距離計算缺失值策略**
  - **問題**: 當 `Common CpG < C_min` 時，距離如何定義？
  - **建議解法**: 定案為 `MAX_DIST (1.0)`。這能確保無重疊的 Reads 在聚類樹中被推得最遠，符合直覺。

- [ ] **AltSupport 的 Base Quality 閾值**
  - **問題**: `min_base_quality` 具體數值？
  - **建議解法**: 設定預設值 `20` (Phred Score，即 1% 錯誤率)。低於此值視為 `UNKNOWN`。

## Low Priority: 優化與擴充 (Future Work)

- [ ] **加權距離 (Weighted Distance)**
  - **想法**: 依據 CpG 的 Genomic Context (e.g., PMD vs HMM) 給予不同權重，而非直接 Gating 排除。
  - **狀態**: 先完成 Gating 版本，行有餘力再做。

- [ ] **多執行緒 I/O 優化**
  - **想法**: Producer-Consumer 模型分離 I/O 與計算。
  - **狀態**: 目前 OpenMP Per-Region 已足夠快，此為後期效能調優選項。
