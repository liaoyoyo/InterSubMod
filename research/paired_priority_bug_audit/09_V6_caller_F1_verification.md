<!--
build_date: 2026-05-12
agent: V6 caller F1 vs SEQC2 truth — complete verification
status: validated
report_class: pipeline-invariance-verification
audience: PI / lab member / 任何質疑 V6 是否破壞 F1 的人
parent_reports:
  - InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md (caller F1 6-version 直接實證)
  - InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md §8.5.2 (F1 三版相同總結)
  - InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md (V6 完整說明)
inputs:
  - 6 phased VCFs (baseline / v3f_no_pononly / V5 @ 0.93 + 0.6 purity)
  - SEQC2 truth set v1.2.1: /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
outputs:
  - 本檔（caller F1 不變完整鐵證）
verdict: caller F1 vs SEQC2 truth 從 ClairS-TO → baseline → V3F → V5 → V6 五階段完全不變（TP/FP/FN 每個小數位都相同）；F1 = 0.7166 @ 0.93 purity, 0.6273 @ 0.6 purity；機制證明：longphase-to phase 階段不改 VCF 的 FILTER 欄位，只動 GT/PS 等 phasing 註記；V6 重用 V5 phased VCF 故 F1 由檔案 identity 保證一致
last_verified: 2026-05-12
-->

# V6 Caller F1 vs SEQC2 — 完整不變性驗證

## 0. TL;DR

> **caller F1 在 ClairS-TO → baseline → V3F → V5 → V6 五階段中完全不變**。@ 0.93 purity F1 = 0.7166（TP 28,509 / FP 11,606 / FN 10,938）；@ 0.6 purity F1 = 0.6273。三層獨立證據：(1) **直接實證**: 4/30 已用同一 SEQC2 truth 對 6 個 phased VCF 跑 F1，6 個版本 TP/FP/FN 每個小數位完全相同；(2) **VCF FILTER 計數**: 3 版本 phased VCF 都有 47,798 PASS variants + 3,187,275 total + 完全相同 FILTER 分布（檔案層級證明 phase 階段不動 FILTER）；(3) **機制證明**: longphase-to phase 階段只動 GT/PS/GT2/GT3（phasing 註記），FILTER 來自 ClairS-TO caller。V6 特有的證據是「**重用 V5 phased VCF**」(literally same file path) → V6 F1 由檔案 identity 數學保證 = V5 F1 = baseline F1 = caller F1。

---

## 1. F1 計算 pipeline 完整流程

```
ClairS-TO caller (run_clairs_to)
    ↓ 產生：raw VCF (含 PASS + LowQual + NonSomatic + ... FILTER 標籤)
    ↓ 例 HCC1395 5kHz: /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz
    ↓ 3,187,275 total variants, 47,798 PASS, 46,793 LowQual;NonSomatic, ...
    
longphase-to phase (binary baseline / V2b / V3F / V5)
    ↓ 讀入 raw VCF + tumor BAM + PoN + reference
    ↓ 對 PASS variants 建 phasing graph + 標 PS (phase set) / GT (genotype) / GT2 / GT3
    ↓ **不改 FILTER 欄位**（ClairS-TO 已標定）
    ↓ 產生：tumor_phased.vcf (含原始 FILTER + 新增 PS/GT/GT2/GT3 註記)

longphase-to haplotag (binary baseline / V3F / V5 / V6)
    ↓ 讀入 phased VCF + tumor BAM
    ↓ 標 BAM reads 的 HP:i: 標籤 (HP1, HP2, HP1-1, HP2-1, HP3, 33)
    ↓ **完全不動 VCF**（haplotag 是 BAM-level 操作）
    ↓ 產生：tumor_tagged.bam

F1 vs SEQC2 truth
    ↓ 從 phased VCF 取 PASS variants
    ↓ 對 SEQC2 truth (high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz)
    ↓ 計算 TP/FP/FN/Precision/Recall/F1
```

→ **F1 只依賴 phased VCF 中的 PASS variants**。

→ **PASS variants 由 ClairS-TO caller 決定**，longphase-to phase 不改、haplotag 完全不碰 VCF。

→ **任何 longphase-to 版本（baseline / V2b / V3F / V5 / V6）都會輸出相同 PASS variants 集合 → F1 數學保證相同**。

---

## 2. 五階段 F1 完全相同 — 三層獨立證據

### 2.1 證據 ① — 4/30 直接實證 6 個 phased VCF F1（最強）

來源報告：[`InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md`](../../docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md) §3a

| Label | Purity | Binary version | TP | FP | FN | Precision | Recall | **F1** |
|---|:---:|---|---:|---:|---:|---:|---:|---:|
| ClairS-TO @ 0.93 raw | 0.93 | (caller only) | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| **A1 baseline @ 0.93** | 0.93 | baseline LongPhase-TO | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| **A3 v3f_no_pononly @ 0.93** | 0.93 | V3F binary（無 PON-only flag）| 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| **A5 V5 @ 0.93** | 0.93 | V5 binary | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| ClairS-TO @ 0.6 raw | 0.6 | (caller only) | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| **B1 baseline @ 0.6** | 0.6 | baseline | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| **B3 v3f_no_pononly @ 0.6** | 0.6 | V3F binary | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| **B5 V5 @ 0.6** | 0.6 | V5 binary | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |

來源 TSV：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phased_vcf_f1.tsv`

→ **每個小數位 TP/FP/FN/Precision/Recall/F1 完全相同**（6 個版本 × 2 purity）。

### 2.2 證據 ② — Phased VCF FILTER 計數（檔案層級鐵證）

直接對 3 個 phased VCF 計數 PASS variants 與 FILTER 分布：

**指令**：
```bash
for v in baseline pononly_v2b threshold_compare/v5_flag; do
    f=/big7_disk/liaoyoyo2001/longphase-to-mod/output/$v/tumor_phased.vcf
    n=$(grep -v "^#" "$f" | awk -F'\t' '$7=="PASS"' | wc -l)
    total=$(grep -v "^#" "$f" | wc -l)
    echo "$v: PASS=$n / total=$total"
done
```

**結果**：

| Phased VCF | PASS variants | Total variants | FILTER 分布相同？ |
|---|---:|---:|---|
| `baseline/tumor_phased.vcf` | **47,798** | **3,187,275** | ✓ |
| `pononly_v2b/tumor_phased.vcf` | **47,798** | **3,187,275** | ✓（與 baseline 每個 FILTER 類別行數完全一致）|
| `threshold_compare/v5_flag/tumor_phased.vcf` | **47,798** | **3,187,275** | ✓（同上）|

FILTER 分布每個類別逐項相同（例：LowQual=550, LowQual;LowAltBQ=88, LowQual;LowAltBQ;MultiHap=2, ...）。

→ **3 版本 phased VCF 在 FILTER 欄位**檔案層級**完全相同**。

### 2.3 證據 ③ — 機制證明：longphase-to phase 不動 FILTER

從 V5 phased VCF header 看 longphase-to 加的新 INFO/FORMAT 欄位：

```
##INFO=<ID=PS,Number=1,Type=Integer,Description="Phase set ID">
##FORMAT=<ID=GT2,Number=1,Type=String,Description="Secondary genotype (somatic phased)">
##FORMAT=<ID=GT3,Number=1,Type=String,Description="Tertiary genotype">
##INFO=<ID=PON,Number=0,Type=Flag,Description="Tagged as non-somatic by PON">
```

longphase-to phase 階段對每個 variant 僅修改：
- GT (FORMAT 欄位 — 加 phasing info, e.g. 0|1)
- PS (FORMAT/INFO — phase set ID)
- GT2 / GT3 (FORMAT — secondary/tertiary genotype 給 somatic ALT)
- 加 INFO PON tag（若該 variant 命中 --pon-file 資料庫）

**FILTER 欄位（第 7 欄）完全不動** — 來自 ClairS-TO 原始輸出（PASS / LowQual / NonSomatic 等）。

源碼確認：`PhasingProcess.cpp` 與 `Phasing.cpp` 中對 vcf record 的修改只動 GT/PS/INFO，FILTER 完全 passthrough（見 `vcfParser.cpp` 中 write_vcf_record）。

---

## 3. V6 特殊性 — 重用 V5 phased VCF（檔案 identity 證明）

V6 與 V3F / V5 不同之處：

| 版本 | Phasing 階段使用的 binary | 產出的 phased VCF |
|---|---|---|
| baseline | baseline longphase-to | `output/baseline/tumor_phased.vcf` |
| V2b | V2b longphase-to | `output/pononly_v2b/tumor_phased.vcf` |
| V3F | V2b longphase-to（V3F 重用 V2b phased VCF）| 同上 V2b |
| V5 | V5 longphase-to（新 phasing 邏輯）| `output/threshold_compare/v5_flag/tumor_phased.vcf` |
| **V6** | **不新跑 phasing — 直接重用 V5 phased VCF** | **= V5 phased VCF**（檔案層級相同）|

V6 命令範本（從 `run_v6_haplotag.sh`）：

```bash
./longphase-to-v6 haplotag \
    -s /big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_phased.vcf \
    -b /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam \
    -r /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
    --log -t 24 \
    -o output/v6_germline_absent_revert/tumor_tagged
```

V6 haplotag 只動 BAM（產 `tumor_tagged.bam` + `.bai`），**完全不碰 VCF**。

→ **V6 F1 = V5 F1**（檔案 identity 數學保證）= **0.7166 @ 0.93 / 0.6273 @ 0.6**

---

## 4. 完整 F1 不變性證據鏈

```
ClairS-TO caller
   ↓ F1 = 0.7166 @ 0.93 (TP=28,509, FP=11,606, FN=10,938)
   ↓ [longphase-to phase: 只加 GT/PS/GT2/GT3 註記，不動 FILTER]
   ↓
baseline phased VCF
   ↓ F1 = 0.7166 ✓（A1 直接實證 + PASS=47,798 計數證明）
   ↓ [V3F haplotag: 不動 VCF，只改 BAM HP tag]
   ↓
V3F BAM
   ↓ F1 = 0.7166 ✓（A3 v3f_no_pononly @ 0.93 直接實證）
   ↓ [V5 phasing: 新邏輯但仍只動 GT/PS/GT2/GT3，不動 FILTER]
   ↓
V5 phased VCF
   ↓ F1 = 0.7166 ✓（A5 V5 @ 0.93 直接實證 + PASS=47,798 計數證明）
   ↓ [V6 haplotag: 重用 V5 phased VCF + 不動 VCF]
   ↓
V6 BAM
   ↓ F1 = 0.7166 ✓（檔案 identity 數學保證 = V5 F1）
```

跨 purity 一致：@ 0.6 purity 同樣 5 階段全 F1 = 0.6273。

---

## 5. 為什麼 phasing 改動不影響 F1

longphase-to phase 階段做的事：
1. 讀入 ClairS-TO PASS variants（已固定的 somatic 候選）
2. 用 BAM reads 連接相鄰 PASS variants → 建 haplotype block
3. 對每個 block 給 PS (phase set) ID
4. 對 somatic ALT variants 給 GT2 (secondary genotype) 與 GT3 (tertiary genotype, ambiguous)
5. （V5+）加 PON tag 給命中 --pon-file 的 variants（INFO 層級，不改 FILTER）

→ **這些動作都不改 FILTER 欄位**。FILTER 由 ClairS-TO caller 一次性決定，longphase-to 純消費。

→ V5 vs baseline 的差異（Pass 2 reclassify germline het, ploidy fix, threshold）只影響 **PS / GT block 結構**，不影響哪些 variants 是 PASS。

→ F1 只看 PASS variants vs SEQC2 truth → F1 不變。

---

## 6. V6 在跨樣本（Phase D）的 F1 caveat

Phase D 4 樣本（H1437/H2009/HCC1954/HCC1937）的 caller F1 vs 各自 truth set **未直接計算**，但同樣邏輯保證 F1 不變：

| 樣本 | 來源 caller | longphase-to phase | F1 vs truth |
|---|---|---|---|
| HCC1395 0.93 / 0.6 | ClairS-TO | V5 binary | **實測** 0.7166 / 0.6273 |
| H1437 | ClairS-TO | V5 binary | 沒實測，但**機制保證** = ClairS-TO caller F1 |
| H2009 | ClairS-TO | V5 binary | 同上 |
| HCC1954 | ClairS-TO | V5 binary | 同上 |
| HCC1937 | ClairS-TO | V5 binary | 同上 |

→ 對任何樣本：**caller F1 = longphase-to phased VCF F1 = ISM 下游 ground-truth F1**（前提：使用同 caller VCF 與 truth set）

→ 4 樣本 V6 BAM 的 ISM downstream marker rate 與 caller F1 是兩個不同概念：
- **caller F1** = ClairS-TO 對 truth 的準確度（不受 longphase-to 版本影響）
- **ISM downstream marker rate** = NG≥3 region 內 TP/FP ratio（受 BAM haplotag 結果影響，這才是 V6 改善的部分）

---

## 7. 操作驗證（用戶可重現）

### 7.1 確認 phased VCF PASS variants 相同

```bash
for v in baseline pononly_v2b threshold_compare/v5_flag; do
    f=/big7_disk/liaoyoyo2001/longphase-to-mod/output/$v/tumor_phased.vcf
    n=$(grep -v "^#" "$f" | awk -F'\t' '$7=="PASS"' | wc -l)
    echo "$v: PASS=$n"
done
# Expected output:
#   baseline: PASS=47798
#   pononly_v2b: PASS=47798
#   threshold_compare/v5_flag: PASS=47798
```

### 7.2 確認 FILTER 分布完全相同

```bash
for v in baseline pononly_v2b threshold_compare/v5_flag; do
    f=/big7_disk/liaoyoyo2001/longphase-to-mod/output/$v/tumor_phased.vcf
    echo "--- $v ---"
    grep -v "^#" "$f" | awk -F'\t' '{print $7}' | sort | uniq -c | head -10
done
# Expected: 同樣 LowQual=550, LowQual;LowAltBQ=88, ... 等
```

### 7.3 確認 V6 重用 V5 phased VCF

```bash
grep -E "tumor_phased.vcf" /big7_disk/liaoyoyo2001/longphase-to-mod/run_v6_haplotag.sh
# Expected: 包含 -s /big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_phased.vcf
```

### 7.4 重跑 F1（若需要再次實證）

```bash
# 已實證的 F1 計算腳本（4/30 完成）
TRUTH=/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
REGION=/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed

for PHASED in baseline pononly_v2b threshold_compare/v5_flag; do
    VCF=/big7_disk/liaoyoyo2001/longphase-to-mod/output/$PHASED/tumor_phased.vcf
    # bcftools / hap.py / vcfeval to compute F1
    # (existing pipeline used: bcftools isec for TP/FP/FN)
done
```

實證結果已存：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/{caller_f1,phased_vcf_f1}.tsv`

---

## 8. 結論

| 問題 | 答案 |
|---|---|
| longphase-to 會輸出 VCF 供 F1 比對嗎？ | ✅ 是（`tumor_phased.vcf`，含原 FILTER + 新 GT/PS）|
| caller F1 從 ClairS-TO 到 baseline 是否變動？ | ✅ 不變（longphase-to phase 不動 FILTER）|
| baseline → V3F F1 變動？ | ✅ 不變（直接實證：A1 = A3 F1 = 0.7166 @ 0.93）|
| V3F → V5 F1 變動？ | ✅ 不變（直接實證：A3 = A5 F1 = 0.7166；PASS 47,798 同；V5 phasing layer 只改 PS/GT 結構）|
| V5 → V6 F1 變動？ | ✅ **檔案 identity 數學保證不變**（V6 重用 V5 phased VCF）|
| 跨 purity（0.93/0.6）一致？ | ✅ 一致（B1=B3=B5 同樣 0.6273）|
| 4 樣本（Phase D）F1 是否驗證？ | ⚠️ 未直接計算，但**機制相同保證 = ClairS-TO caller F1**（V5 binary 跑 phasing, V6 binary 跑 haplotag, 同 V3F/V5 邏輯）|
| TP / FP / FN / Precision / Recall 是否可算？ | ✅ 可算（已有 4/30 報告完整數據 + bcftools / hap.py / vcfeval 標準 pipeline）|
| V6 修改是否會影響 F1？ | ❌ **不會**（V6 不改 phasing 不改 VCF，只改 BAM HP tag）|

**結論**：caller F1 vs SEQC2 truth 在 ClairS-TO → baseline → V3F → V5 → V6 五階段**完全不變**（直接實證 6 個版本 × 2 purity + 檔案層級 FILTER 計數 + 機制證明三層保證）。

V6 patch 只動 BAM haplotag 層，對 caller F1 數學上不可能有影響。Phase D 4 樣本未直接計算 F1 但機制保證同樣不變。

---

## 9. 引用文件

- 4/30 V3F ablation F1 報告：`InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md` §3a
- 5/8 主報告 §8.5.2 F1 三版相同：`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`
- V6 完整說明：`InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md` §8 Verdict
- F1 數據 TSV：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phased_vcf_f1.tsv`
- SEQC2 truth set：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- V6 haplotag 命令：`/big7_disk/liaoyoyo2001/longphase-to-mod/run_v6_haplotag.sh`
