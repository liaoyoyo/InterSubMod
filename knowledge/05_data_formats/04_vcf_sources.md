---
id: ism-kb-05-data-formats-vcf-sources
name: "VCF Sources 與符號連結對照"
description: "7 樣本 × 3 pipeline 的 VCF 來源路徑、符號連結關係、歷史錯誤矯正紀錄（2026-04-04 TO VCF 來源錯誤矯正）。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "VCF source paths by sample and mode; symlink traversal"
related_ids:
  - ism-kb-05-data-formats-index
  - ism-kb-08-truth-and-benchmark-truth-set-registry
  - ism-kb-03-pipelines-tumor-only
tags: [data-format, vcf, sources, symlink, caller]
canonical_paths: [05_data_formats/04_vcf_sources.md]
alias_paths: []
---

# VCF Sources 與符號連結對照

- 一句結論：paired 用 ClairS paired output，TO 用 ClairS-TO output；2026-04-04 前 TO 曾誤用 symlink 指向 paired pileup，已矯正
- 適用對象：確認 VCF 來源正確性、重跑 benchmark 前
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 追蹤 symlink 實際指向
  readlink -f /big8_disk/liaoyoyo2001/longphase-to-mod/output/clairsto_v3fixed_work/clairsto_tp.vcf.gz
  ```

---

## Caller 類型對應

| Pipeline | Caller | 原始 VCF 路徑 |
|----------|--------|--------------|
| paired_full | ClairS (paired, full model) | `/big8_disk/liaoyoyo2001/longphase-s/output/<sample>/paired_v5/filtered_snv_tp.vcf.gz` |
| paired_pileup | ClairS (paired, pileup only) | 同上但檔名為 `filtered_snv_pileup_tp.vcf.gz` |
| TO | ClairS-TO | `/big8_disk/liaoyoyo2001/longphase-to-mod/output/clairsto_v3fixed_work/clairsto_tp.vcf.gz` |

---

## 7 樣本 VCF 路徑速查

| 樣本 | Paired VCF base 目錄 | TO VCF 位置 |
|------|---------------------|-------------|
| HCC1395 | `/big8_disk/liaoyoyo2001/longphase-s/output/HCC1395/paired_v5/` | `/big8_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/hcc1395/` |
| HCC1395_DORADO | `.../HCC1395_DORADO/paired_v5/` | `.../hcc1395_dorado/` |
| COLO829 | `.../COLO829/paired_v5/` | `.../colo829/` |
| H1437 | `.../H1437/paired_v5/` | `.../h1437/` |
| H2009 | `.../H2009/paired_v5/` | `.../h2009/` |
| HCC1937 | `.../HCC1937/paired_v5/` | `.../hcc1937/` |
| HCC1954 | `.../HCC1954/paired_v5/` | `.../hcc1954/` |

**驗證**：`ls <path>` 確認檔案存在；歷史路徑可能變動

---

## 歷史錯誤矯正（2026-04-04）

**錯誤**：TO canonical run 的 VCF symlink 誤指向 paired pileup output（非 ClairS-TO）
**發現時間**：2026-04-04
**影響**：所有 2026-04-04 前的 TO canonical run 實際上用的是 paired pileup VCF
**矯正方式**：重建 symlink 指向正確 ClairS-TO 輸出；重跑 affected runs
**記錄**：MEMORY: `project_vcf_source_error_correction`

**如何驗證沒中招**：
```bash
# 讀 run 目錄的 source_vcf_file 記錄
cat output/canonical/<sample>/TO/<run>/master_manifest.yaml  # 若有

# 或從 master_dataset 查 source_vcf_file 欄位
python3 -c "
import pandas as pd
df = pd.read_csv('.../all_region_rows.tsv.gz', sep='\t', low_memory=False)
to_df = df[df['mode'] == 'TO']
print(to_df[['sample', 'source_vcf_file']].drop_duplicates())
"
```

---

## 篩選 TP/FP 流程

原始 VCF（含所有 caller 輸出）→ 用 truth VCF 分類：

```bash
# 典型流程（scripts/pipeline/steps/02_filter_analysis.py 類）
python3 scripts/pipeline/steps/03_filter_analysis.py \
  --caller-vcf <caller_output.vcf.gz> \
  --truth-vcf <truth.vcf.gz> \
  --hc-bed <truth_hc.bed> \
  --output-tp <filtered_snv_tp.vcf.gz> \
  --output-fp <filtered_snv_fp.vcf.gz>
```

**TP 定義**：caller 輸出 ∩ truth VCF（適用 HC BED 區）
**FP 定義**：caller 輸出 ∖ truth VCF

---

## 常見陷阱

- ❌ **誤用 symlink 不知指向**：跑 benchmark 前先 `readlink -f` 確認
- ❌ **混用 paired_full 的 VCF 跑 TO pipeline**：結果不一致
- ❌ **忽略 HC BED 過濾**：部分樣本（COLO829）無 HC BED，需用 `num_callers≥2` 邏輯替代

---

## 相關

- Truth set 全表：[../08_truth_and_benchmark/01_truth_set_registry.md](../08_truth_and_benchmark/01_truth_set_registry.md)
- TO pipeline 警告：[../03_pipelines/03_tumor_only.md](../03_pipelines/03_tumor_only.md)
- F1 計算：[../08_truth_and_benchmark/02_f1_calculation.md](../08_truth_and_benchmark/02_f1_calculation.md)
