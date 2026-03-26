<!--
建立時間: 2026-03-26
目標: 重構後 chr19 功能等效性驗證
處理範圍: chr19:29283968 (HCC1395 TP site)，重構分支 refactor/phase1-safety
關聯檔案:
  - docs/refactor_baseline/BASELINE_SUMMARY.md
  - docs/refactor_baseline/post_phase3/refactoring_final_report.md
-->

# chr19 功能驗證報告

**驗證時間**：2026-03-26
**分支**：`refactor/phase1-safety`（最終 commit: `7d4f425`）
**目的**：確認 Phase 1~3 + Task 1.4 全部修正後，核心統計輸出維持正確

---

## 執行設定

| 參數 | 值 |
|------|----|
| VCF | HCC1395 pileup filtered_snv_tp.vcf.gz |
| Site | chr19:29283968 (C→T, TP) |
| BAM | HCC1395 tumor + normal |
| Reference | hg38 |
| Distance metric | BERNOULLI |
| Window size | ±1000 bp |
| Binary | `/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod` |

---

## 執行結果

**Status：✅ 正常完成（exit 0），1.07 秒**

### 關鍵統計指標

| 指標 | 數值 | 說明 |
|------|------|------|
| NumReads | 85 | 40 forward + 45 reverse |
| NumCpGs | 11 | 覆蓋到的 CpG 位點數 |
| ClusterPermanovaF | **166.23** | 聚類結構顯著性（偽 F 統計量） |
| ClusterPermanovaP | 0.01 | 顯著（p ≤ 0.01）|
| ClusterPermanovaValid | true | PERMANOVA 有效執行 |
| HPMergedDelta | 0.1010 | HP1/HP2 甲基化差異 |
| HPMergedP | 0.001 | HP 差異顯著 |
| HPMergedSig | true | |
| LabelHPPermanovaF | 27.18 | Label HP PERMANOVA F |
| LabelHPPermanovaP | 0.01 | 顯著 |
| AlleleDelta | -0.003 | Allele 差異（無 LOH 信號）|
| AlleleP | 1.0 | 無 Allele 顯著差異 |
| VerificationClass | **Strong** | |
| Quality_Score | 85.0 | High 品質 |
| Quality_Tier | High | |

### Distance Matrix 統計

| 指標 | 數值 |
|------|------|
| Valid read pairs | 3,462 |
| Invalid pairs（coverage 不足）| 108 |
| Valid pair ratio | 97.0% |
| Avg common CpG coverage | 7.94 |

---

## 結論

chr19:29283968 在重構後分支執行**正常**：
- 無 crash、無 NaN/Inf 輸出
- PERMANOVA pseudo-F = 166.23（遠大於 0，Task 1.4 修正有效）
- VerificationClass = Strong，與 HCC1395 TP 位點預期吻合
- 統計計算路徑（FisherExact log-sum-exp、PERMANOVA 直接公式）均正常輸出

**重構分支已可提交 PR 合併至 main。**
