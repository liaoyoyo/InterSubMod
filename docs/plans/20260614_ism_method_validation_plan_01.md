---
title: "ISM 距離→分群→顯著性 方法驗證計劃（第二/三部分）"
date: 2026-06-14
type: plan / methodology-validation
status: DRAFT — 待用戶 review 後逐步執行
scope: 單樣本 HCC1395 paired（後續多樣本 COLO829 等）
binary_baseline: develop 7c341d2（SKIP fix + U1 + gap1）
related: docs/experiments/in_progress/2026/06/20260614_u1_maxdist_vs_skip/{README,gap2,method_audit}
---

# ISM 方法驗證計劃 — 距離→分群→顯著性判定

> **工作原則**：先整理計劃 → 逐項驗證 → 每項「確認現狀→思考改進→驗證」。從最上游可穩定確認往下游走。
> **驗證地基**：HP tag = germline-haplotype 真值 ground truth（normal-only 最乾淨）；HP-AUC 為距離有效性尺。
> ⚠ **共同 caveat**：HP=germline 真值，HP-AUC 只證 germline-HP 捕捉；somatic subclone **無 ground truth**，須間接判別。

## 已確立的事實（本 session 已驗，作為計劃前提）
| 項 | 結論 | 證據 |
|----|------|------|
| 距離 metric | BERNOULLI≈NHD（HP-AUC 0.641 vs 0.636），中間帶處理對 germline-HP 無關 | method_audit/hp_auc.json |
| SKIP vs MAX_DIST | **SKIP 去污染 distance HP-AUC 0.64→0.835**；應改預設 | hp_auc.json |
| clustering 鏈 | SKIP distance 強(0.835)→CramersV 1.0 有基礎，非過度擬合 | hp_auc + U1 |
| tumor HP 結構 | **67.7% region tumor both-HP / 32.3% 單-HP**；germline 主導因 tumor HP 差小非無 HP | 全基因組 tumor 狀況 |
| tumor 複雜度 | LOH 8.5% + CNV(Gain 8.7%/Low 8%/ampl 9.6%) — 須分層 | 同上 |
| HP somatic 判定缺陷 | hp_residual_sig=tumor sig（未檢定殘差）| RegionProcessor.cpp:978 |

## 待驗證項目（優先序：上游穩定→下游/方法改進）

### P1 — 確認最終判別「基本有效」的地基（最上游，先做）
- **P1a 全基因組 HP-AUC**：把 4-region HP-AUC 擴到全基因組/分層（normal-only / tumor both-HP / 單-HP），確認「距離抓 germline-HP」普遍成立 + SKIP 0.835 優勢穩定。
  - 方法：對 N 個 region 批次算 HP-AUC（需 distance matrix 輸出，重跑子集）。可行性：中（需重跑開 matrix，限子集如 chr1）。
- **P1b SKIP 改預設**（cpp-change）：`--nan-distance-strategy` 預設 MAX_DIST→SKIP；MAX_DIST 標實驗。regression golden 需同步（golden 改用 SKIP 或保留 MAX_DIST 對照）。
  - 風險：改預設影響所有下游；須 regression 重建 + 確認無破壞。

### P2 — somatic subclone 是否存在（你的核心問題）
- **P2a tumor both-HP region 的 tumor-specific HP 訊號**：67.7% both-HP region，算 tumor-only HP-AUC（tumor 自己 HP1 vs HP2 距離分離）。
  - 若 tumor-only AUC 顯著 >0.5 且 > 對應 normal → tumor-specific（somatic 候選）。
  - caveat：須扣 germline baseline（normal 同 region HP-AUC）。
- **P2b 修正 hp_residual 顯著性**（cpp-change，對應方法缺陷 #1）：對「殘差 = tumor HP delta − normal HP delta」做**正規 permutation 檢定**（非借 tumor p）→ 才能可信數 somatic 比例。
- **P2c LOH/CNV 分層**：subclone 候選須排除 LOH/CNV confound（CNV 改變 depth→偽 HP 不平衡）。用 Coverage_Category + Potential_LOH 分層。

### P3 — 「單結構兩標籤 / 一標籤兩結構」偵測（你列的 subclone 核心）
- **P3a 一個甲基 cluster 內混兩種 HP/allele 標籤**：cluster vs label 交叉表，找「單一 cluster 內 label 混合」= 甲基結構不分 HP（subclone 跨 HP？或甲基同質）。
- **P3b 一個 HP 標籤內兩個甲基 cluster**：HP1-only reads 是否再分兩甲基群 = haplotype 內 subclone（tumor subclone within germline HP）。
  - 這是最接近「somatic subclone」的訊號：germline HP 內出現 tumor-specific 甲基亞群。
  - 方法：對單一 HP 的 reads 再跑 sub-clustering，測顯著性。可行性：中（需 per-HP 子分群）。

### P4 — 甲基位點與結構的關聯（你列的）
- **P4a 哪些 CpG 驅動分群**：per-CpG 對 cluster/HP 的判別力（已有 PerCpgAsm？Fisher）。確認「是哪些位點」。
- **P4b 甲基-結構組合**：cluster × (HP/allele/sample/LOH) 的關聯模式。

### P5 — tumor/normal 同位點甲基差異驗證輸出（你列的）
- **P5a 同 CpG 的 tumor vs normal 甲基狀況比對**：對每個 region 的每個 CpG，算 tumor 甲基率 vs normal 甲基率 → 差異 = somatic 甲基改變（DMR-like）。
  - 這直接做「位點層 tumor/normal 差異」，不經 HP（補充 HP-axis 的另一視角）。
  - 輸出：per-CpG tumor/normal 甲基差異表 + 顯著位點。
  - 可行性：高（methylation.csv 已有 tumor/normal reads，按 is_tumor 分組算 per-CpG 均值差）。
- **P5b normal-anchored cis-test 對齊**：與既有 NormalBaseline（normal-anchored cis-test，memory reference_msa_vs_ism_tool）整合，判真 cis。

## 方法學可行性總表
| 項 | 資料來源 | 重算需求 | 可行性 |
|----|---------|---------|:---:|
| P1a HP-AUC 全基因組 | distance matrix | 重跑開 matrix（子集）| 中 |
| P1b SKIP 預設 | — | cpp-change + regression | 中 |
| P2a tumor-only AUC | distance matrix + reads | 重跑子集 | 中 |
| P2b 殘差檢定 | — | cpp-change | 中 |
| P3 結構×標籤 | cluster labels + reads | summary + 子分群 | 中 |
| P4 per-CpG | methylation.csv + PerCpgAsm | 已有 | 高 |
| P5 tumor/normal 位點 | methylation.csv + reads(is_tumor) | 已有 | **高** |

## 建議執行順序（逐步，每步確認後下一步）
1. **P1a**（全基因組 HP-AUC 確認地基）→ 確立「距離+SKIP 對 germline-HP 有效」普遍成立
2. **P5a**（tumor/normal 同位點甲基差異）→ 最高可行性、直接的 somatic 甲基視角、不依賴 HP
3. **P3b**（HP 內子分群）→ 最接近 somatic subclone 的偵測
4. **P2b/P1b**（方法修正：殘差檢定 + SKIP 預設）→ cpp-change
5. P2a/P2c/P4/P5b → 補充

## 未決 / 需用戶確認
- 執行順序是否照 P1a→P5a→P3b→修正？
- P1b（SKIP 改預設）何時做（影響 regression baseline）？
- scope：先 chr1 還是直接全基因組？多樣本何時納入？
