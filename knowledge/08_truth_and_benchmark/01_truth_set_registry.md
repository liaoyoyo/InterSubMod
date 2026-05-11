---
id: ism-kb-08-truth-and-benchmark-truth-set-registry
name: "Truth Set Registry（7 樣本）"
description: "7 樣本 truth set 完整表：VCF 路徑、HC BED 可用性、TRUTH_TOTAL、來源、platform；F1 計算必備。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "truth set paths verified via ls on /big8_disk/data/"
related_ids:
  - ism-kb-08-truth-and-benchmark-index
  - ism-kb-08-truth-and-benchmark-f1-calculation
  - ism-kb-02-samples-index
  - ism-kb-05-data-formats-vcf-sources
  - ism-kb-02-samples-hcc1395
  - ism-kb-02-samples-colo829
  - ism-kb-02-samples-h1437
  - ism-kb-02-samples-h2009
  - ism-kb-02-samples-hcc1937
  - ism-kb-06-workflows-f1-benchmark
  - ism-kb-06-workflows-adding-new-sample
tags: [truth, benchmark, registry, seqc2, nygc, truth-vcf]
canonical_paths: [08_truth_and_benchmark/01_truth_set_registry.md]
alias_paths: []
---

# Truth Set Registry（7 樣本）

- 一句結論：HCC1395=SEQC2、COLO829=NYGC、H1437/H2009/HCC1937/HCC1954=orthogonal-tools；各有 TRUTH_TOTAL 用於 recall 分母
- 適用對象：F1 計算、論文 methods 章節
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  bcftools view /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz | \
    grep -v "^#" | wc -l
  ```

---

## 7 樣本 Truth Set 完整表

| 樣本 | Platform | Truth VCF 路徑 | HC BED | TRUTH_TOTAL | 來源 |
|------|---------|---------------|--------|-------------|------|
| **HCC1395** | ONT 5kHz | `/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz` | ✅ | **39,447** | SEQC2 consortium |
| **HCC1395_DORADO** | ONT Dorado | 同上 | ✅ | 39,447 | 同上 |
| **COLO829** | ONT PAO (R10) | `/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz` | ❌ | **41,427** | NYGC（Illumina 驗證） |
| **H1437** | ONT | `/big8_disk/data/H1437/orthogonal-tools-benchmark/...somatic-only.vcf.gz` | ✅ | **90,129** | Orthogonal tools benchmark |
| **H2009** | ONT | 同上格式 | ✅ | **168,529** | 同上 |
| **HCC1937** | ONT | 同上格式 | ✅ | **60,691** | 同上 |
| **HCC1954** | ONT | 同上格式 | ✅ | **26,567** | 同上 |

**總計**：7 樣本（HCC1395_DORADO 共用 HCC1395 truth）

---

## Truth Set 來源細節

### SEQC2（HCC1395）
- Consortium：The SEQC2 Somatic Variant Call Consortium
- 版本：high-confidence_sSNV_in_HC_regions_v1.2.1
- 特色：多 caller 共識 + 人工審核
- HC BED：`/big8_disk/data/HCC1395/SEQC2/high-confidence_regions_v1.2.bed`
- 文獻：Fang et al. 2021 Nature Biotechnology

### NYGC（COLO829）
- Producer：New York Genome Center
- 版本：NovaSeq vs NovaSeqBL 配對，v6 annotated
- 特色：無 HC BED；需用 `num_callers ≥ 2` 過濾

**替代過濾**：
```bash
bcftools view -f PASS -i 'INFO/num_callers >= 2' <NYGC.vcf.gz>
```

### Orthogonal Tools（H1437/H2009/HCC1937/HCC1954）
- 方法：多 caller orthogonal validation
- 特色：每樣本有 HC BED
- 變異數偏多（H2009 特別大：168K）

---

## 驗證路徑存在性

```bash
# 檢查所有 truth VCF 存在
for sample in HCC1395 COLO829 H1437 H2009 HCC1937 HCC1954; do
  # 以 HCC1395 為例的路徑模板
  f=$(find /big8_disk/data/$sample -name "*.vcf.gz" -path "*truth*" -o -name "*somatic-only*" -o -name "*high-confidence*" 2>/dev/null | head -1)
  echo "$sample: $f"
done
```

---

## TRUTH_TOTAL 重新計算

```bash
# HCC1395
bcftools view -f PASS /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz | \
  grep -v "^#" | wc -l

# COLO829
bcftools view -f PASS -i 'INFO/num_callers >= 2' \
  /big8_disk/data/COLO829/NYGC/...vcf.gz | grep -v "^#" | wc -l
```

**注意**：`-f PASS` 篩選；若 truth VCF 沒 PASS filter，直接數全部

---

## HC BED 的意義

**HC BED** = High-confidence BED regions

**用途**：限制 F1 計算的 region 範圍到「truth set 信心度高」的區域

**效果**：
- 僅在 HC BED 內的 caller 輸出被納入 TP/FP
- HC BED 外的 caller 輸出被忽略（不算 FP）

**影響**：
- F1 通常比全基因體的 F1 **高**（因為排除 truth 無信心的區域）
- 論文標準做法

---

## COLO829 無 HC BED 的特殊處理

**情況**：NYGC 提供的 truth 無 HC BED

**解法**：用 `num_callers ≥ 2`（多 caller 共識）作為 confidence proxy
```bash
bcftools view -i 'INFO/num_callers >= 2' <NYGC.vcf.gz> > confident.vcf.gz
```

**影響**：COLO829 F1 可能較其他樣本低（因為沒有 HC BED 排除低信心區）

---

## 權威文件

👉 **`docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md`** 為最新權威文件

本 KB 只提供結構化摘要；所有異動應同步更新權威文件。

---

## 常見誤區

- ❌ **混用 VCF 與 HC BED**：HC BED 是 region 過濾，不是 VCF 過濾
- ❌ **用過濾後的 caller VCF 算 TRUTH_TOTAL**：TRUTH_TOTAL 永遠是原始 truth VCF 的變異數
- ❌ **COLO829 用 HCC1395 truth**：各樣本 truth 不通用

---

## 相關

- F1 計算：[02_f1_calculation.md](02_f1_calculation.md)
- Benchmark 協議：[03_benchmark_protocols.md](03_benchmark_protocols.md)
- 樣本詳情：[../02_samples/](../02_samples/)
- 權威：[../../docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md](../../docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md)
