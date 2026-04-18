<!--
建立時間: 2026-04-11 04:30
目標: VCF 來源錯誤完整影響評估與結論盤點
處理範圍: 所有使用 pileup symlink 的 ISM runs、O1-O13 觀察、已發表結論
關聯檔案:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_*/significance_summary.csv
  - /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/to_pileup/
-->

# VCF 來源錯誤完整影響評估報告

## 1. 錯誤根因

```
/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup
  → /big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup
```

- **ClairS_ssrs** = ClairS **paired** 模式 (tumor + normal BAM)，ssrs = synthetic + real samples model
- **ClairS-TO 沒有 pileup 模型**（因為無 normal BAM）
- Symlink 在專案初期建立時指向 ClairS paired 結果，後來開始 TO 分析時未更新

### VCF Header 差異

| 項目 | ClairS paired (錯) | ClairS-TO (對) |
|------|-------------------|----------------|
| `##source` | `ClairS` | `ClairS-TO` |
| 版本 | `v0.4.1` | `v0.3.0` |
| FORMAT | GT,AF,DP,**NAF,NDP,NAD** | GT,AF,DP（無 N 系列）|
| FILTER | PASS,LowQual,RefCall,**Germline** | PASS,**NonSomatic**,LowQual,+9 others |
| TP (vs SEQC2 HC) | 30,490 | 28,509 |
| FP (vs SEQC2 HC) | 4,842 | 11,606 |
| FP 性質 | **mapping artifacts** | **germline variants** |

## 2. 影響範圍界定

### 2.1 受影響的 ISM Runs（longphase-to-mod/output/）

| ISM Run | BAM | TP | FP | VCF | 判定 |
|---------|-----|-----|-----|-----|------|
| `ism_baseline_{tp,fp}` | Baseline | 30,476 | 4,822 | ClairS paired | **錯誤** |
| `ism_pononly_v2b_{tp,fp}` | V2b | 30,476 | 4,822 | ClairS paired | **錯誤** |
| `ism_v3_fixed_{tp,fp}` | V3-Fixed | 30,476 | 4,822 | ClairS paired | **錯誤** |
| `ism_pononly_{tp,fp}` | PON-Only v1 | 30,476 | 4,822 | ClairS paired | **錯誤** |

### 2.2 不受影響的數據（2026-04-11 完整驗證）

| 數據來源 | 理由 | 驗證方式 |
|---------|------|---------|
| **Master dataset (748K rows)** | 來源=canonical outputs，HCC1395 TO TP=28,495 FP=11,601 | `all_region_rows.tsv.gz` groupby 直接確認 |
| **7 樣本 canonical TO data** | 全部使用各自的 ClairS-TO VCF | HCC1395 run_context.json 確認 caller=ClairS-TO |
| **O1-O10 觀察** | 基於 master dataset | 同上 |
| **O11-O13 觀察** | 使用 master dataset | 同上 |
| **G1-G7 Germline FP NO-GO** | 使用 master dataset | 同上 |
| **TO FP Provenance** | 使用正確 ClairS-TO 11,598 FPs | 報告內數字一致 |

**Master dataset 7 樣本 TO 行數驗證（2026-04-11 確認）**：

| Sample | TO TP | TO FP | 狀態 |
|--------|-------|-------|------|
| HCC1395 | 28,495 | 11,601 | ✓ 正確 ClairS-TO |
| HCC1395_DORADO | 28,856 | 11,572 | ✓ |
| COLO829 | 33,089 | 17,528 | ✓ |
| H1437 | 45,473 | 13,442 | ✓ |
| H2009 | 125,706 | 11,989 | ✓ |
| HCC1937 | 12,623 | 12,032 | ✓ |
| HCC1954 | 17,068 | 50,218 | ✓ |

### 2.3 已修正的 ISM Runs

| ISM Run | BAM | TP | FP | FP Catch | TP Loss | ISM F1 | SEQC2 F1 | SEQC2 ΔF1 |
|---------|-----|-----|-----|----------|---------|--------|---------|-----------|
| `ism_baseline_clairsto` | Baseline | 28,383 | 11,830 | 0.77% | 0.37% | 0.0151 | 0.7117 | -0.0009 |
| `ism_pononly_v2b_clairsto` | V2b | 28,383 | 11,830 | 0.89% | 0.44% | 0.0174 | 0.7115 | -0.0011 |
| `ism_v3fixed_clairsto` | V3-Fixed | 28,495 | 11,601 | 0.63% | 0.39% | 0.0124 | 0.7153 | -0.0012 |
| `ism_longphase_to_pass` | V3-Fixed | 28,495 | 19,236 | 0.48% | 0.39% | 0.0096 | 0.6527 | -0.0010 |

### 2.4 錯誤 VCF 的虛假數據（留作對照）

| ISM Run | TP | FP | FP Catch | TP Loss | ISM F1 |
|---------|-----|-----|----------|---------|--------|
| `ism_baseline` (WRONG) | 30,476 | 4,822 | **5.18%** | 0.37% | **0.0965** |
| `ism_pononly_v2b` (WRONG) | 30,476 | 4,822 | **4.38%** | 0.45% | **0.0816** |
| `ism_v3_fixed` (WRONG) | 30,476 | 4,822 | **5.08%** | 0.39% | **0.0945** |

## 3. 受影響結論盤點

### 已無效結論

| 結論 | 出處 | 無效原因 |
|------|------|---------|
| ~~V3-Fixed ISM F1 = 0.0945 首次超越 Paired F1 = 0.0909~~ | V3-Fixed 驗證報告 | 使用 paired FP (4,822) 而非 TO FP (11,601)；正確 F1 = 0.0124 |
| ~~ISM F1 +15.8% 改善~~ | 同上 | 虛假改善；正確數據顯示改善微乎其微 |
| ~~ISM FP Catch Rate 5.08%~~ | ism_v3_fixed | 正確值為 0.63%（paired FP = mapping artifacts 容易捕獲；TO FP = germline 幾乎不可能） |

### 仍然有效的結論

| 結論 | 理由 |
|------|------|
| HP tag 分布改善 (Balanced +9.5pp) | 不依賴 VCF 來源 |
| HP:i:33 正確產出 | 不依賴 VCF 來源 |
| TO QS AUC = 0.507 (隨機) | 2026-04-05 用正確 ClairS-TO VCF 確認 |
| ISM 對 TO germline FP 捕獲率 < 1% | 全部正確 VCF runs 一致確認 |
| TO FP = germline variants，ISM 無法區分 | 60+ 特徵全 AUC < 0.64 |
| O1-O10 全部觀察結論 | 來源 = canonical master dataset (不受影響) |
| O11-O13 觀察結論 | 同上 |
| G1-G7 Germline FP NO-GO | 同上 |
| Self-phasing causal chain | 不涉及 VCF TP/FP 分類 |
| LOH penalty 是 QS 失效根因 | 同上 |

## 4. 核心發現：ISM 為何對 TO FP 無效

### Paired FP vs TO FP 本質差異

| 維度 | Paired FP (4,822) | TO FP (11,601) |
|------|-------------------|----------------|
| 來源 | Mapping artifacts, technical noise | Germline variants (PoN 未涵蓋) |
| 甲基化模式 | 異常（artifact 讀段無生物意義） | **正常**（與真 variant 相同）|
| HP tag 模式 | Balanced（artifact 均勻分布）| **也 balanced**（germline = biallelic）|
| ISM 可偵測性 | 可（pattern deviation） | **不可**（pattern 與 TP 相同）|

**這解釋了 ISM F1 差距**：
- 用 paired FP → ISM F1=0.0945（能抓到 mapping artifacts）
- 用 TO FP → ISM F1=0.0124（幾乎抓不到 germline variants）

## 5. BAM 版本位置

| 版本 | 路徑 | 大小 | 建立日期 | 修正內容 |
|------|------|------|---------|---------|
| Baseline | `output/baseline/tumor_tagged.bam` | 260G | 2026-04-03 | 原始 longphase-to |
| V2b | `output/pononly_v2b/tumor_tagged.bam` | 260G | 2026-04-03 | PON-only phasing v2b |
| V3-Fixed | `output/pononly_v3_fixed/tumor_tagged.bam` | 268G | 2026-04-10 | +getVote 兩層 +countINDEL guard |
| V4 Alt Guard | `output/pononly_v4_alt_guard/tumor_tagged.bam` | (building) | 2026-04-11 | +countSNP altHaplotype guard |

根目錄：`/big7_disk/liaoyoyo2001/longphase-to-mod/output/`

### Binary 版本

| Binary | 路徑 | 內容 |
|--------|------|------|
| 原始 baseline | `longphase-to-baseline` | 無修改 |
| 最新版 (V4) | `longphase-to` | getVote 兩層 + countINDEL guard + countSNP altHaplotype guard |

## 6. countSNPHaplotype altHaplotype Guard 修正分析

### 修正內容

`countSNPHaplotype()` 的 alt 分支缺少 `HAPLOTYPE_UNDEFINED` guard（countINDELHaplotype 已有）。

**修正前**：`countMap[haplotypeBase.altHaplotype]++` 當 altHaplotype=-1 → **陣列越界存取 (UB)**
**修正後**：加入 `if(haplotypeBase.altHaplotype != HAPLOTYPE_UNDEFINED)` guard

### 影響範圍

Phased VCF 中 1,338,737 個 unphased variants（42%）：
- `0/1`：932,433 個 — 程式碼仍正確處理（`gtValue[0]='0', gtValue[2]='1'` → altHaplotype=HP2）
- **`1/1`：406,304 個 — altHaplotype 保持 UNDEFINED**

只有 406,304 個 `1/1` homozygous variants 受影響。讀段匹配 ALT base 時觸發 `countMap[-1]++`（UB）。

### HP Tag 分布比較（chr1:1-10M 抽樣）

| HP tag | Baseline | V3-Fixed | V3→V3 變化 |
|--------|----------|----------|-----------|
| HP:i:1 | 26,126 | 26,206 | +80 |
| HP:i:2 | 13,598 | 13,682 | +84 |
| HP:i:11 | 1,357 | 1,158 | -199 |
| HP:i:21 | 1,256 | 1,070 | -186 |
| HP:i:33 | 0 | **385** | +385 (new) |

### V3-Fixed vs V4 比較結果（2026-04-11 完成）

**Haplotag 全局統計**：
| 指標 | V3-Fixed | V4 Alt Guard | 差異 |
|------|---------|-------------|------|
| Total alignment | 40,859,727 | 40,859,727 | 0 |
| Tagged | 18,805,977 | **18,805,977** | **0** |
| Untagged | 22,053,750 | 22,053,750 | 0 |
| HP:i:33 (chr1:1-10M) | 385 | 385 | 0 |

**ISM 結果**：
| 指標 | V3-Fixed | V4 Alt Guard |
|------|---------|-------------|
| TP SF | 112 | **112** |
| FP SF | 73 | **73** |
| ISM F1 | 0.0124 | **0.0124** |
| SEQC2 F1 | 0.7153 | **0.7153** |

**結論**：V3-Fixed 和 V4 結果**完全相同**。`countSNPHaplotype` 的 `altHaplotype=-1` UB 在 V3-Fixed binary 中實際上沒有影響 haplotype voting 結果。guard 修正是防禦性修正（消除 UB），但不改變行為。

**建議**：V4 BAM 可以移除，V3-Fixed BAM 即為最佳版本。

## 7. LongPhase-TO 版本完整比較

| Config | TP | FP | TP SF | FP SF | FP Catch | TP Loss | ISM F1 | SEQC2 F1 | ΔF1 |
|--------|-----|-----|-------|-------|---------|---------|--------|---------|------|
| Baseline (orig) | 28,383 | 11,830 | 105 | 91 | 0.77% | 0.37% | 0.0151 | 0.7117 | -0.0009 |
| V2b (PON-only) | 28,383 | 11,830 | 125 | 105 | 0.89% | 0.44% | 0.0174 | 0.7115 | -0.0011 |
| **V3-Fixed** | 28,495 | 11,601 | 112 | 73 | 0.63% | 0.39% | 0.0124 | **0.7153** | -0.0012 |
| V4 (=V3-Fixed) | 28,495 | 11,601 | 112 | 73 | 0.63% | 0.39% | 0.0124 | 0.7153 | -0.0012 |

### 評估結論

1. **SEQC2 Calling F1**: V3-Fixed (0.7153) > Baseline (0.7117) > V2b (0.7115)
   - V3-Fixed 較好是因為 phased VCF 品質提升導致更多 TP 落入 HighConf (28,495 vs 28,383)
2. **ISM SuggestFilter**: 所有版本對 TO germline FP 的捕獲率 < 1%，差異不具實質意義
3. **ISM 對 SEQC2 F1 影響**: 全部略微負面（ΔF1 -0.0009 to -0.0012）
4. **V4 = V3-Fixed**：altHaplotype guard 是純防禦性修正，不改變結果

## 8. 後續行動

1. **V4 BAM 完成後** → 執行 ISM (ClairS-TO VCFs) → 比對 V3 vs V4
2. **評估** baseline vs V3-Fixed vs V4 哪版 BAM 較好
3. **等用戶確認後** 移除不需要的 BAM 版本
4. **Scripts 已更新** — `run_vcf_all_snv.sh` 和 `run_batch_vcf_analysis.sh` 已加入 `--caller-mode to|paired` 參數
