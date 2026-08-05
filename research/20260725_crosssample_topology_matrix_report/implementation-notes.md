<!--
建立時間：2026-07-25
目標：追蹤七資料集跨樣本矩陣與 HTML 報告的設計決定、偏離、折衷與未決事項
處理範圍：chr1–22；7 technical datasets；6 biological IDs；21 technical pairs
關聯檔案：research/20260725_crosssample_topology_matrix_report/pre-decision-audit.md
-->

# 跨樣本 sSNV／區域拓撲矩陣 implementation notes

## 設計決定

- Task Type B comprehensive validation；服務 G4（多樣本一致性）與 G5（可外部驗證）。
- 將「相似」拆為 candidate sites、active sites、W sites、primary undirected edges、exact components、marginal raw-VAF distribution、exact-PS topology profile 與 HCC coordinate-matched topology。
- 同時保留 7 technical datasets × 21 pairs 與 6 biological IDs 去重分析；同癌種推論以後者為主。
- HTML 使用 canonical Data Analytics artifact 與 portable builder；每張圖旁附 denominator、解釋與 claim ceiling。

## 偏離

- 不把 SEQC2 Figure 4 當成可直接擬合的 ground truth；本輪只回答跨樣本相似與目前 topology resolution。
- 不把 phase-invariant undirected endpoint edge 稱為 evolutionary parent→child edge。

## 折衷

- marginal raw-VAF distribution 可比較不同位點集合，但只能表示頻譜相似；不能作 clone identity。
- exact-PS profile 可覆蓋 7 datasets，但不含 locus identity；coordinate-matched tree comparison 只對 HCC same-ID pair 有合理共同母群。
- 同癌種只有 4 個 biological-ID pairs，任何 p-value 均標為探索性，不作一般化證明。

## 未決

- 是否新增 matched downsampling／PS opportunity-matching ablation。
- 是否取得第二組同 biological ID、不同 basecaller／library 的 replicate。
- 是否以單細胞或 SEQC2-compatible CCF/CN model 作 global fusion-tree 正交驗證。

## 交付與驗證

- 產出 standalone HTML：`research/20260725_crosssample_topology_matrix_report/20260725_七資料集_sSNV區域拓撲比對矩陣_01.html`。
- Canonical compact artifact validation：18 datasets、27 sources、status ready。
- Official structural verification：24 blocks、5 charts、4 tables，embedded artifact exact match。
- Local Playwright QA：桌機 1440×1000 與手機 390×844 共 11/11 checks pass；無外部請求、page error、console error 或 document-level overflow。
- 21-pair heatmap 在窄螢幕保留完整 21 欄，改由圖內水平捲動；未刪減矩陣資料。
- Official static-SVG extraction／browser verifier 在本環境逾時；因此以 structural verifier 加獨立 Playwright QA 補足，限制已寫入 `final_delivery_receipt.json`。
