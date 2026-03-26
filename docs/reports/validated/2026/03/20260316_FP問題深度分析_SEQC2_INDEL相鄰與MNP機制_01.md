<!--
建立時間: 2026-03-16
目標: 詳細分析 Phase 4 Case Study 發現的兩類 FP 問題的規模、機制、定義與可修正性
處理範圍: FP-B1（SEQC2 INDEL 鄰近型）、FP-B2（MNP 拆分型）全面量化分析
關聯檔案:
  - docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md
  - docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md
  - scripts/analysis/screen_mnp_adjacent_fp.py
-->

# FP 問題深度分析：SEQC2 INDEL 相鄰型與 MNP 拆分型

**日期**：2026-03-16
**資料集**：HCC1395 5kHz simplex 5mCG_5hmCG（ClairS ssrs paired mode）
**分析動機**：Phase 4 Case Study 觀察到兩類系統性 FP 機制（FP-B1 chr7:52087777 與 FP-B2 chr9:75383880），需量化其規模、確認機制、判斷可修正性

> **2026-03-17 校正註記**
>
> 本文件中關於 **FP-B2 = 乾淨 MNP 拆分型** 的早期假說，已被後續正式驗證修正。  
> 根據 `20260316_IGV固定template_AI視覺初篩與正式驗證_01.md` 與 `20260316_igv_case_validation_01.tsv`，FP-B2 在主 tumor BAM 下更像 **`target SNV + -1 deletion-like gap / local alignment anomaly`**，不是乾淨 MNP。  
> 因此閱讀本文件時，請以後續修正文件為準，將本文件視為「早期假說形成過程」的紀錄。

---

## 一、基礎數字概覽

| 資料集 | 總數 |
|--------|------|
| ClairS FP SNV | **4,842** |
| ClairS TP SNV | **30,490** |
| SEQC2 sINDEL truth (PASS+HighConf) | **1,625** |
| SEQC2 sSNV truth (PASS) | **39,447** |
| ClairS INDEL PASS | **7,940** |
| ClairS INDEL PASS + LowQual | **12,077** |

---

## 二、問題 A：SEQC2 INDEL 鄰近型 FP SNV（Benchmark Gap）

### 2.1 定義與機制

**觸發條件**：某個 SEQC2 HighConf PASS 的 sINDEL 位點（已知真實體細胞 INDEL）旁邊 ±2bp 處，ClairS 呼叫了一個 SNV，但該 SNV 不在 SEQC2 sSNV truth set 中。

**生物機制**：
```
Reads 攜帶 INDEL allele 在 alignment 時，
相鄰位置的 base 可能因 soft-clip 或 alignment 偏移
→ 造成 ClairS 在相鄰位置感知到少量「SNV-like」錯配
→ ClairS 呼叫假性 SNV
```

**Benchmark gap 本質**：SEQC2 SNV truth set 設計上**刻意不評估 INDEL ±2bp 內的 SNV**，因複雜變異評估困難且容易產生假陽性。因此這類 FP 屬於「評估未定義」區域，非 ClairS 真正的生物學錯誤。

### 2.2 數量統計

| 條件 | 數量 | 佔 FP 比例 |
|------|------|-----------|
| FP SNV 在 SEQC2 sINDEL ±2bp 內 | **3** | 0.06% |
| TP SNV 在 ClairS INDEL ±2bp 內 | 5 | 0.02% |
| SEQC2 sSNV 在 SEQC2 sINDEL ±2bp 內 | 6 | 0.015% |

### 2.3 三個案例詳情

#### 案例 A-1：chr7:52087777（FP-B1，Phase 4 重點案例）

| 屬性 | 值 |
|------|-----|
| ClairS FP SNV | chr7:52087777 A>T, GQ=0.67, AF=0.1068 (11/103 reads) |
| 相鄰 SEQC2 INDEL | chr7:52087776 TA>T (+1 offset), HighConf PASS |
| SEQC2 INDEL TVAF | 0.184，63 samples 支持，bwaClassification=Strong |
| ClairS INDEL call | chr7:52087776 TA>T，**LowQual** (GQ=7, AF=0.1359, 14/103 reads) |
| SEQC2 SNV 查詢 | chr7:52087777 **無 entry**（benchmark gap 確認） |

**分析**：ClairS 確實偵測到了 INDEL（但為 LowQual），同時在 +1 位置呼叫了假性 SNV。ClairS INDEL 的 AF（13.6%）與 SNV 的 AF（10.7%）相近，顯示兩者來自相同讀段。這是 INDEL alignment artifact 造成的 SNV 假陽性。

**類型判斷**：**定義 + 處理雙重問題**
- Benchmark 不評估此位置（定義層）
- ClairS 應優先呼叫 INDEL 而非在相鄰位置呼叫 SNV（處理層）

---

#### 案例 A-2：chr20:22472592

| 屬性 | 值 |
|------|-----|
| ClairS FP SNV | chr20:22472592 C>G, GQ=0.57, AF=0.0606 (4/66 reads) |
| 相鄰 SEQC2 INDEL | chr20:22472591 GC>G (+1 offset), HighConf PASS |
| SEQC2 INDEL TVAF | **0.557**（高 VAF），63 samples 支持，bwaClassification=Strong |
| ClairS INDEL call | chr20:22472591 GC>G，**PASS** (GQ=15.4, AF=0.5577, 29/52 reads) |
| SEQC2 SNV 查詢 | chr20:22472592 **無 entry** |

**分析**：ClairS 正確以 PASS 呼叫了這個 INDEL（AF=55.8%，高 VAF）。SNV call（AF=6.1%，4 reads）是極低支持度的 artifact，來自攜帶 INDEL 的讀段在相鄰位置的殘留錯配。

**特點**：這是最乾淨的 INDEL-adjacent SNV FP 案例——INDEL 已被 ClairS 正確識別，SNV 是副產物。
**可修正性**：在 ClairS INDEL PASS 的相鄰 ±2bp 位置 post-filter SNV calls，可移除此 FP。

---

#### 案例 A-3：chr7:54585815

| 屬性 | 值 |
|------|-----|
| ClairS FP SNV | chr7:54585815 C>A, GQ=0.66, AF=0.0941 (8/85 reads) |
| 相鄰 SEQC2 INDEL | chr7:54585814 AC>A (+1 offset), HighConf PASS |
| SEQC2 INDEL TVAF | **0.245**，62 samples 支持（全部 60 個 bwa/bowtie/novo PASS） |
| ClairS INDEL call | **無（未呼叫）** |
| SEQC2 SNV 查詢 | chr7:54585815 **無 entry** |

**分析**：ClairS 在此完全未偵測到 INDEL（TVAF=0.245，理論上應可 call），卻呼叫了相鄰 SNV（AF=9.4%）。
**類型判斷**：主要是 ClairS 處理問題（INDEL missed + adjacent SNV artifact）。
**可修正性**：無法用「相鄰 INDEL PASS filter」救回，因 ClairS INDEL 未呼叫。

---

### 2.4 A 型問題小結

| 案例 | SEQC2 INDEL | ClairS INDEL | SNV AF | 主因 | 可後處理修正？ |
|------|-------------|--------------|--------|------|--------------|
| chr7:52087777 | HighConf PASS | LowQual | 10.7% | 雙重 | 否（INDEL LowQual） |
| chr20:22472592 | HighConf PASS | **PASS** | 6.1% | Alignment artifact | **是** |
| chr7:54585815 | HighConf PASS | 未呼叫 | 9.4% | ClairS INDEL miss | 否 |

**可 post-filter 修正的 Type A FP：1 / 3**
**總 Type A 規模：3 / 4842 = 0.06% of FP**（對 F1 影響 < 0.001）

---

### 2.5 SEQC2 Benchmark 設計的一致性確認

SEQC2 sSNV truth set 中：**6 / 39,447 sSNV 在 sINDEL ±2bp 內（0.015%）**

這 6 個案例的 ClairS 狀態：

| 位置 | 類型 | ClairS 狀態 | 備註 |
|------|------|------------|------|
| chr1:188602967 C>T | HighConf | TP | 正確 call |
| chr3:28490117 A>T | MedConf | FN | 未 call |
| chr4:35608608 C>A | MedConf | TP | 正確 call |
| chr11:107495093 C>G | MedConf | TP | 正確 call |
| chr18:30042532 T>C | MedConf | TP | 正確 call |
| chr19:24297950 C>A | HighConf | FN | 未 call |

**結論**：SEQC2 benchmark 設計非常保守，INDEL 鄰近的 SNV 極少（0.015%），但 ClairS 能正確 call 其中 4/6（67%）。此設計保守度支持「這類 SNV 評估困難」的假設。

---

## 三、問題 B：MNP 拆分型 FP SNV（相鄰 SNV 被誤判）

### 3.1 定義與機制

**觸發條件**：兩個緊鄰的體細胞 SNV 存在於同一條單倍型上（= MNP, Multi-Nucleotide Polymorphism），SNV caller 只 call 其中一個，相鄰位置被 ClairS 錯誤解讀為獨立 FP SNV 或 INDEL。

**生物機制（以 FP-B2 chr9:75383880 為例）**：
```
真實情況：
  HP1 上同時存在 C>A (pos 75383879) + T>A (pos 75383880)
  → MNP: CTA → AAA（連續兩個 SNV 共 phasing 於 HP1）

ClairS 呼叫結果：
  SNV: T>A at 75383880 → PASS（part of MNP，被 call 出來）
  INDEL: AC>A at 75383878 → LowQual（= C 被解讀為 deletion）
         （實際是 C>A SNV 被誤解為 deletion artifact）

SEQC2 truth set:
  兩個位置都不在 SNV 或 INDEL truth set 中
  → ClairS 的 T>A at 75383880 被評估為 FP
```

### 3.2 數量統計

| 條件 | 數量 | 佔 FP 比例 |
|------|------|-----------|
| FP SNV 在 ClairS INDEL ±2bp 內（含重疊） | **7** | 0.14% |
| 其中同時在 SEQC2 INDEL 附近（Type A 重疊） | 2 | — |
| 純 Type B（ClairS INDEL 相鄰，無 SEQC2 INDEL） | **5** | 0.10% |

### 3.3 五個純 Type B 案例

| FP SNV | ClairS INDEL | INDEL Filter | offset | SEQC2? | 備註 |
|--------|-------------|-------------|--------|--------|------|
| chr9:75383880 T>A | chr9:75383878 AC>A | LowQual | +2 | 無 | FP-B2，C>A SNV 被誤解為 del |
| chr13:21328786 T>C | chr13:21328788 大缺失 | **PASS** | -2 | 無 | 大缺失（32bp）鄰近 SNV |
| chr4:126681681 C>A | chr4:126681683 GAGGGT>G | LowQual | -2 | 無 | 缺失鄰近 SNV |
| chr8:94149032 G>T | chr8:94149031 T>TATATTTA | LowQual | +1 | 無 | 短串聯重複插入鄰近 SNV |
| chr9:99663637 A>T | chr9:99663635 A>AT | LowQual | +2 | 無 | 插入鄰近 SNV |

**Type B 機制分類**：

1. **MNP 型**（chr9:75383880）：兩個相鄰 SNV，一個被 call 為 SNV FP，另一個被誤解為 INDEL
2. **INDEL 鄰近 artifact 型**（chr13:21328786, chr4:126681681, chr8:94149032, chr9:99663637）：ClairS 偵測到 INDEL（部分 LowQual），在相鄰位置產生 alignment artifact SNV

### 3.4 FP-B2 深度驗證（BAM 證據）

**直接 BAM 分析結果**（來自 Phase 4 Case Study）：

| 讀段組別 | ALT 讀段（T>A at 75383880） | REF 讀段（ref T） |
|---------|---------------------------|--------------------|
| pos 75383879（ref=C） | **全部 15/15 reads 顯示 A** | 保持 C |
| pos 75383880（ref=T） | 顯示 A（SNV call） | 保持 T |
| CIGAR | 無 indel 符號 | 無 indel 符號 |
| 正常 BAM | 無 A at 75383879 | — |

**結論**：chr9:75383879 的 C>A 是真實的 SNV（非 alignment artifact），與 75383880 的 T>A 構成 MNP。SEQC2 不含此 MNP，可能是：
(a) SEQC2 短讀取平台無法可靠偵測相鄰 SNV
(b) 此 MNP 在 HCC1395 短讀取數據中被遺漏

---

## 四、兩類問題的對比分析

### 4.1 性質比較

| 維度 | Type A（SEQC2 INDEL 相鄰） | Type B（MNP 拆分） |
|------|---------------------------|--------------------|
| 數量 | 3 / 4842 (0.06%) | 5 / 4842 (0.10%) |
| SEQC2 中有 INDEL？ | ✅ 是（HighConf PASS） | ❌ 否 |
| ClairS 有 INDEL call？ | 2/3 有（1 PASS + 1 LowQual） | 5/5 有（全部 LowQual 或 PASS） |
| 機制 | INDEL alignment artifact → 假性 SNV | 相鄰 SNV 之一或 INDEL 鄰近 artifact |
| 是定義問題？ | ✅ 部分（benchmark 排除 INDEL ±2bp） | ❌ 否（SEQC2 完全沒有這些位置） |
| 是處理問題？ | ✅ 部分（ClairS 應優先選擇 INDEL） | ✅ 是（MNP 無法拆分 call） |
| 可 post-filter 修正？ | 1/3 案例（ClairS INDEL PASS 相鄰） | 部分可（ClairS INDEL LowQual 相鄰） |

### 4.2 整體規模評估

| 指標 | 值 |
|------|-----|
| 兩類 FP 合計 | **8 / 4842 = 0.17%** |
| 對 F1 影響（估計） | **< 0.001**（統計上可忽略） |
| 是否系統性？ | ✅ 是（可預測，非隨機） |
| 對生物學解釋影響 | 輕微（定性理解有幫助，定量無顯著影響） |

---

## 五、是否可修正：詳細評估

### 5.1 Benchmark 層修正

**方案**：將 ClairS SNV FP 評估中，位於 SEQC2 sINDEL ±2bp 的位置標記為「評估未定義（Undefined）」，從 FP 計數中排除。

**影響**：
- Type A FP 減少：3 個
- TP 不受影響（SEQC2 SNV 中僅 6 個在 sINDEL ±2bp，ClairS 已正確 call 4 個為 TP）
- F1 變化：< 0.001（忽略不計）
- 適合做法：針對特定案例的深度分析時有意義，不影響整體 benchmark 結論

**結論**：**可修正，但價值有限**（3 個 FP 恢復）

---

### 5.2 ClairS Post-processing 層修正

**方案 P1**：在 ClairS SNV 後處理中，移除位於 ClairS INDEL PASS ±2bp 內且低 AF（< 0.15）的 SNV calls

**影響估計**：
- FP 移除：1-2 個（chr20:22472592 確定，其他需驗證）
- TP 誤傷：需確認（5 個 TP SNV 在 ClairS INDEL ±2bp 內，但 AF 通常不低）
- 整體影響：極小

**方案 P2**：實作 MNP merger — 識別 ±2bp 內相鄰兩個 SNV 在同一讀段上 co-occurring → 合併為 MNP entry

**難度**：高（需要讀段級別 phasing 分析）
**影響**：直接移除 Type B FP，但需確保不誤傷真實 TP SNV 旁的另一個真實 TP SNV

---

### 5.3 InterSubMod 層修正

本分析發現的兩類 FP 問題與 InterSubMod 甲基化分析無直接關聯：

- 這些 FP 的判定來自 SEQC2 benchmark，與甲基化特徵無關
- InterSubMod 的改進方向（CramersV in VerificationClass, HP-Allele 判別等）針對不同的 FP 類型
- 若 INDEL-adjacent 或 MNP FP 進入 InterSubMod 分析，其甲基化特徵可能正常或異常，需個案判斷

**建議**：InterSubMod 可以選擇性地標記 `ClairS INDEL PASS ±2bp` 或 `ClairS FP INDEL LowQual ±2bp` 的位置，在 metadata 中加入 `near_indel` 標記，供分析時參考。

---

## 六、最終結論

### 問題本質

| 問題 | 結論 |
|------|------|
| Type A 是定義問題還是處理問題？ | **雙重問題**：Benchmark 設計排除 INDEL ±2bp（定義層），ClairS 在 INDEL 鄰近產生假性 SNV（處理層） |
| Type B 是定義問題還是處理問題？ | **主要是處理問題**：ClairS 無法合併相鄰 SNV 為 MNP，導致其中一個無法在 SEQC2 中找到對應 |
| 判斷出錯是正常還是可修正？ | **正常且可修正**：此類複雜變異（MNP、INDEL-adjacent）是短讀取 benchmark 的已知挑戰；ONT 長讀取有能力修正，但需要 MNP caller |

### 規模結論

**兩類 FP 合計 8 個（0.17% of 4842 FP）**，對 F1 影響可忽略，但具有重要的生物學意義：
- 揭示了 ONT 5kHz 在偵測 MNP 時的系統性缺口
- 確認了 SEQC2 SNV benchmark 對 INDEL 鄰近區域的刻意排除
- 提供了 ClairS → INDEL caller 整合的改進方向

### 可操作建議

1. **短期**：記錄 8 個案例為 `INDEL_adjacent` 或 `MNP_candidate` 類型，不計入核心 benchmark 評估
2. **中期**：對所有 4,842 FP 執行 `screen_mnp_adjacent_fp.py`（已建立），估計 MNP 候選總數
3. **長期**：評估 ClairS MNP calling 支援或 post-processing MNP merger 的可行性

---

## 附錄：關鍵 VCF 查詢指令

```bash
# FP-B1：chr7:52087777
tabix /big8_disk/.../filtered_snv_fp.vcf.gz chr7:52087770-52087785
tabix /big8_disk/liaoyoyo2001/ClairS_ssrs_output/indel.vcf.gz chr7:52087770-52087785
tabix /big8_disk/data/HCC1395/SEQC2/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz chr7:52087770-52087785

# FP-B2：chr9:75383880
tabix /big8_disk/.../filtered_snv_fp.vcf.gz chr9:75383875-75383885
tabix /big8_disk/liaoyoyo2001/ClairS_ssrs_output/indel.vcf.gz chr9:75383875-75383885
tabix /big8_disk/data/HCC1395/SEQC2/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz chr9:75383875-75383885

# Type A chr20:22472592
tabix /big8_disk/.../filtered_snv_fp.vcf.gz chr20:22472586-22472598
tabix /big8_disk/liaoyoyo2001/ClairS_ssrs_output/indel.vcf.gz chr20:22472586-22472598
```

---

*分析工具：Python 3 + pysam + tabix*
*ClairS 版本：ssrs（paired mode）*
*SEQC2 truth set：v1.2.1 (high-confidence_sSNV/sINDEL_in_HC_regions)*
