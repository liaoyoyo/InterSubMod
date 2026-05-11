---
id: ism-kb-03-pipelines-paired-pileup
name: "Paired Pileup Pipeline"
description: "ClairS pileup VCF + LongPhase-S haplotag + ISM；pileup 模型評估用途。與 paired_full 的差異主要在 VCF 模型複雜度（變異數通常較少）。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "pipeline command against scripts/run_vcf_all_snv.sh"
related_ids:
  - ism-kb-03-pipelines-index
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-pipeline-comparison
  - ism-kb-04-parameters-cli-arguments
tags: [pipeline, paired, pileup, clairs, model-evaluation]
canonical_paths: [03_pipelines/02_paired_pileup.md]
alias_paths: []
---

# Paired Pileup Pipeline

- 一句結論：與 paired_full 同架構，差別在 VCF 來源為 ClairS pileup model（變異數較少）；用於評估 pileup 模型本身的 ISM 效能
- 適用對象：比較 pileup vs full model 效能時；pileup 模型調參時
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big8_disk/liaoyoyo2001/longphase-s/output/<sample>/paired_v5/ | grep pileup
  ```

---

## 與 paired_full 的差異

| 面向 | paired_full | paired_pileup |
|------|-------------|---------------|
| VCF 來源 | ClairS paired full (pileup + full 複合) | ClairS pileup only |
| 典型變異數 | 多（full 加了 pileup 未 call 到的） | 少 |
| 檔名 | `filtered_snv_tp.vcf.gz` | `filtered_snv_pileup_tp.vcf.gz` |
| 模型複雜度 | 高 | 低 |
| ISM ΔF1 | +0.0112（locked） | 類似 paired_full（小差異） |

**其他部分**（Normal BAM、LongPhase-S haplotag、ISM 參數、輸出結構）與 paired_full 完全相同。

---

## 典型命令

```bash
./scripts/run_vcf_all_snv.sh \
  --caller-mode paired \
  --mode all-with-w5000 \
  --vcf-source pileup \
  -o output/paired_pileup_run
```

或直接指定 VCF：
```bash
./build/bin/inter_sub_mod \
  -t <tumor_haplotagged.bam> \
  -n <normal.bam> \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  -v <filtered_snv_pileup_tp.vcf.gz> \
  --loh-bed <phased.bed> \
  -w 5000 -j 120 \
  --distance-metric BERNOULLI \
  -o output/canonical/<sample>/paired_pileup/<timestamp>
```

---

## 輸出結構

同 [01_paired_full.md#輸出結構](01_paired_full.md)，只是目錄名改為 `paired_pileup/`。

---

## 何時使用

✅ **適合**：
- 評估 ClairS pileup 模型本身的 ISM 效能
- 探索 pileup vs full model 在 F1 上的差異
- Reproducibility 研究（pileup 模型較穩定）

❌ **不適合**：
- 論文 benchmark 主表（應用 paired_full）
- 與 TO pipeline 直接比較（應分別比 paired_full）

---

## 已知 canonical run 位置

所有 7 樣本在 `output/canonical/<sample>/paired_pileup/` 下皆有 canonical run（副次角色）。詳見 [../02_samples/](../02_samples/)。

---

## 相關

- 索引：[00_index.md](00_index.md)
- 對照表：[04_pipeline_comparison.md](04_pipeline_comparison.md)
- Paired Full 權威文件：[01_paired_full.md](01_paired_full.md)
