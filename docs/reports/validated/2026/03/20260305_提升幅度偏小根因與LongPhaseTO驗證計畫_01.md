<!--
建立時間: 2026-03-05 18:00
目標: 回答「為何全樣本提升幅度很小」並提出 LongPhase-TO 純樣本的可驗證計畫
處理範圍:
  - 純 tumor（s-pure-pileup）現有 7 樣本結果之根因診斷
  - 現有特徵可分性與門檻敏感度判讀
  - LongPhase-TO 純樣本可行性與驗證步驟
關聯檔案:
  - scripts/analysis/small_gain_diagnosis.py
  - scripts/analysis/run_small_gain_diagnosis.sh
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418/data/small_gain_diagnosis_*.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418/data/paired_vs_to_caller_partial_benchmark.tsv
-->

# 提升幅度偏小根因與 LongPhase-TO 驗證計畫

## 1. 問題定義

你提出的核心疑問：
1. 為何全樣本放寬後提升仍很小？是門檻太嚴格，還是特徵本身不夠明顯？
2. 這個方向是否本質上困難？
3. 在 LongPhase-TO（tumor-only）純樣本下，是否可能有更好結果？如何驗證？

---

## 2. 已完成的定量診斷

### 2.1 產出檔案

診斷表位於：

`/big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418/data/`

- `small_gain_diagnosis_headroom.tsv`
- `small_gain_diagnosis_current_rule_precision.tsv`
- `small_gain_diagnosis_feature_auc.tsv`
- `small_gain_diagnosis_global_auc.tsv`
- `paired_vs_to_caller_partial_benchmark.tsv`（目前 4 個樣本完成）

### 2.2 根因 A：多數樣本可提升天花板本來就很低

根據 `small_gain_diagnosis_headroom.tsv`，若「理想化」只移除 FP 且不誤刪 TP，多數樣本的 F1 理論上限增幅不高：

- `H2009`: max delta `+0.001009`
- `H1437`: max delta `+0.001681`
- `HCC1954`: max delta `+0.002706`

代表即便是完美過濾器，某些樣本也不可能得到大幅 F1 上升。

### 2.3 根因 B：特徵可分性跨樣本不一致，HCC1954 特別弱

根據 `small_gain_diagnosis_feature_auc.tsv`：

- `HCC1395`：`VAF` oriented AUC `0.927`（可分性高）
- `HCC1954`：`VAF` oriented AUC `0.562`、`AlleleDelta_abs` oriented AUC `0.517`（接近隨機）

這說明同一規則在不同樣本的判別力差異很大，無法用單一全局門檻兼顧。

### 2.4 根因 C：現行規則在 HCC1954 為「過度激進」

根據 `small_gain_diagnosis_current_rule_precision.tsv`：

- `HCC1954`：triggered 93 筆中，FP only 4 筆（FP precision `0.043`），TP 誤刪 89 → F1 delta `-0.002179`
- `HCC1395`：triggered 182 筆中，FP 178 筆（FP precision `0.978`），TP 誤刪 4 → F1 delta `+0.002081`

結論：不是「所有樣本門檻都太嚴格」，而是「同一門檻在不同樣本落在不同分布區」，在 HCC1954 變成誤刪 TP。

---

## 3. 這個方向是否本質困難？

### 3.1 單一全局硬門檻：困難且上限有限

綜合上面三點，可判定：
1. 多數樣本 headroom 小
2. 可分性高度 sample-dependent
3. 單一門檻容易在弱可分樣本（如 HCC1954）產生負效益

因此「全樣本共用單一硬門檻且期待明顯 F1 提升」本質上困難。

### 3.2 但方向不是無效：在高可分樣本仍有提升

`HCC1395` 顯示同類特徵在部分樣本可有效工作，重點是要：
1. sample-aware（或分群）門檻
2. 使用組合訊號（VAF + 甲基特徵），而非單一條件

---

## 4. LongPhase-TO 純樣本是否可能更好？

### 4.1 可行性前提（已檢查）

1. 純 tumor BAM 的 MM/ML：7 樣本皆可抽到 MM/ML（可做甲基分析）
2. `ClairS_TO_v0_3_0/snv.vcf.gz`：7 樣本皆存在
3. 現有 `tmp/phasing_output` 多有 phased VCF，但 `phased_bam_output` 多為空，代表後續仍需補 `haplotag`

### 4.2 快速 caller 層比較（目前完成 4 樣本）

`paired_vs_to_caller_partial_benchmark.tsv`：

| sample | paired_f1 | to_f1 | delta(to-paired) | to_fp / paired_fp |
|---|---:|---:|---:|---:|
| COLO829 | 0.8698 | 0.7068 | -0.1629 | 6.32x |
| H1437 | 0.9143 | 0.8213 | -0.0929 | 68.13x |
| HCC1395 | 0.8443 | 0.7166 | -0.1277 | 8.12x |
| HCC1395_DORADO | 0.8565 | 0.7226 | -0.1339 | 19.69x |

解讀：
1. TO 並非天然更好，caller 基線 F1 明顯較低（主要因 FP 大增）
2. 但 TO 也帶來更大的「可被過濾掉的 FP 空間」（headroom 變大）
3. 是否能「反超 paired」取決於後段過濾能否高精度移除大量 FP 且少誤刪 TP

---

## 5. 驗證計畫（可直接執行）

## Phase A：固定現況診斷（已可重跑）

1. 執行：
   `scripts/analysis/run_small_gain_diagnosis.sh /big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418`
2. 判讀：
   - headroom 是否先天受限
   - 規則 trigger precision 是否低於可接受下限
   - AUC 是否低於可用門檻（例如 < 0.62）

## Phase B：LongPhase-TO 單樣本打通（先 HCC1395_DORADO）

1. 以 `ClairS_TO_v0_3_0/snv.vcf.gz` + tumor BAM 執行：
   - LongPhase-TO `phase`（含正確 `--caller` 與 PON）
   - LongPhase-TO `haplotag`（產生 tagged BAM）
2. 以 truth 做 isec，產生 TO 版 `filtered_snv_tp.vcf.gz` / `filtered_snv_fp.vcf.gz`
3. 套用既有 InterSubMod + 過濾分析流程，輸出與 LongPhase-S 同結構報表
4. 比較：
   - baseline F1、filtered F1、delta F1
   - TP 誤刪率、FP 去除率、trigger precision
   - 特徵 AUC（是否較 S 流程更高）

## Phase C：跨樣本擴展（7 純 tumor）

1. 完整跑 7 樣本 TO
2. 與 S 流程做配對比較（per sample）
3. 建立決策分流：
   - 若某樣本 TO 可分性佳且增益穩定，納入 TO 分支
   - 否則維持 S 為主、TO 作補充

---

## 6. 成功判準與決策門檻

建議驗證門檻：
1. `trigger_precision_fp >= 0.70`（避免 HCC1954 類型誤刪）
2. `tp_removed_rate <= 0.2%`
3. `filtered_f1 - baseline_f1 >= +0.001`（至少達可觀察改善）
4. 跨樣本不允許出現顯著負效益（例如 < -0.0015）

若未達門檻，應優先改成：
1. sample-aware 門檻
2. 以 `VAF + AD` 為主、`CV` 降為輔助訊號
3. 引入 TO 專屬特徵（如 PON 命中、Verdict、H/LOH）再做分類

---

## 7. 原始目標可實現性判斷

1. **可實現**：建立完整可重跑診斷流程、定位失效樣本與影響因子（已完成）
2. **條件式可實現**：在部分樣本取得穩定小幅提升（已見 HCC1395）
3. **高難度**：單一全局硬門檻在全部樣本都顯著提升（現階段不現實）
4. **可驗證但未證實**：LongPhase-TO 純樣本是否可帶來更好終局 F1（需完成 Phase B/C）

---

## 8. 知識庫依據

1. 根據 `Knowledge/05_tools/longphase_to.md`：LongPhase-TO 為 tumor-only phasing/LOH 工具，`phase` 需正確 caller 與 PON 設定。
2. 根據 `Knowledge/06_workflows/phasing_workflow.md`：TO 標準流程為 `ClairS-TO -> longphase-to phase -> longphase-to haplotag`。
3. 根據 `Knowledge/02_samples/HCC1395.md`：HCC1395 與 HCC1395_DORADO 的甲基標記型態不同（5mCG+5hmCG vs 5mCG），可造成特徵分布差異。
4. 根據 `Knowledge/03_file_formats/vcf_clairs_to.md`：TO 的 VCF 具 PON/Verdict/H 等欄位，可作為後續模型特徵擴充候選。

