---
id: ism-kb-06-workflows-adding-new-sample
name: "Adding a New Sample"
description: "新樣本入庫流程：truth set 登記、VCF 符號連結、canonical run 跑法、更新 KB 文件。"
status: active
last_verified: 2026-04-22
content_nature: reference
doc_type: howto
verified_scope: "sample addition workflow derived from existing 7 samples structure"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-02-samples-index
  - ism-kb-08-truth-and-benchmark-truth-set-registry
tags: [workflow, new-sample, onboarding, truth-set]
canonical_paths: [06_workflows/06_adding_new_sample.md]
alias_paths: []
---

# Adding a New Sample

- 一句結論：新樣本 = truth set 登記 + VCF symlink 建立 + 3 pipeline canonical 跑 + KB 文件新增
- 適用對象：引入新 cell line 或新 platform
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 查看既有樣本清單
  ls /big7_disk/liaoyoyo2001/big7_disk_output/canonical/
  ```

---

## 9 步驟流程

### Step 1: 確認資料齊全
```bash
# 必備檔案
ls /big8_disk/data/<NEW_SAMPLE>/bam/tumor.bam           # Tumor BAM
ls /big8_disk/data/<NEW_SAMPLE>/bam/normal.bam          # Normal BAM (paired)
ls /big8_disk/data/<NEW_SAMPLE>/truth/*.vcf.gz          # Truth set
ls /big8_disk/data/<NEW_SAMPLE>/truth/*.bed             # HC BED (若有)
```

### Step 2: 登記 Truth Set
編輯 `docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md`：
```markdown
| <NEW_SAMPLE> | ONT X | /big8_disk/data/<NEW_SAMPLE>/truth/*.vcf.gz | ✓/✗ | <TRUTH_TOTAL> |
```

**TRUTH_TOTAL 計算**：
```bash
bcftools view -f PASS <truth.vcf.gz> | grep -v "^#" | wc -l
```

### Step 3: 執行 ClairS paired 產出 VCF
（由 ClairS pipeline 負責，不屬於 ISM）

### Step 4: 建立 VCF symlink
```bash
cd /big8_disk/liaoyoyo2001/longphase-s/output/
mkdir -p <NEW_SAMPLE>/paired_v5/
cd <NEW_SAMPLE>/paired_v5/
ln -s /path/to/clairs/output/<sample>/filtered_snv_tp.vcf.gz .
ln -s /path/to/clairs/output/<sample>/filtered_snv_tp.vcf.gz.tbi .
ln -s /path/to/clairs/output/<sample>/filtered_snv_fp.vcf.gz .
```

### Step 5: 執行 LongPhase-S haplotag
（LongPhase-S pipeline；不屬於 ISM）
產出 `tumor_tagged.bam` 與 `phased.bed`

### Step 6: 跑 3 條 ISM canonical pipeline
```bash
# paired_full
./scripts/run_batch_vcf_analysis.sh \
  --caller-mode paired --mode all-with-w5000 \
  --sample <NEW_SAMPLE> \
  -o output/canonical/<NEW_SAMPLE>/paired_full/$(date +%Y%m%d)

# paired_pileup
./scripts/run_batch_vcf_analysis.sh \
  --caller-mode paired --mode all-with-w5000 \
  --vcf-source pileup \
  --sample <NEW_SAMPLE> \
  -o output/canonical/<NEW_SAMPLE>/paired_pileup/$(date +%Y%m%d)

# TO（需先跑 ClairS-TO + LongPhase-TO）
./scripts/run_batch_vcf_analysis.sh \
  --caller-mode to --mode all-with-w5000 \
  --sample <NEW_SAMPLE> \
  -o output/canonical/<NEW_SAMPLE>/TO/$(date +%Y%m%d)
```

### Step 7: 更新 Master Dataset 聚合
```bash
python3 scripts/analysis/build_loh_round1_cross_sample_audit.py \
  --include-sample <NEW_SAMPLE>
```

### Step 8: 新增 KB 文件

建立 `knowledge/02_samples/0X_<new_sample>.md`：
- 複製既有樣本文件（如 `01_hcc1395.md`）為模板
- 更新 frontmatter、platform、canonical run 路徑、F1 baseline、sample-specific caveats

更新：
- `knowledge/02_samples/00_index.md` 新增一行
- `knowledge/08_truth_and_benchmark/01_truth_set_registry.md` 新增一列
- `knowledge/README.md`「關鍵術語」或樣本總覽區（若提及 7 樣本需改為 8）

### Step 9: 跑 KB 驗證腳本
```bash
python3 knowledge/scripts/validate_frontmatter.py knowledge/
python3 knowledge/scripts/check_related_ids_symmetry.py knowledge/
python3 knowledge/scripts/check_canonical_paths.py knowledge/
```

---

## Checklist

- [ ] Truth set + HC BED 登記
- [ ] TRUTH_TOTAL 確認
- [ ] VCF symlink 建立
- [ ] tumor_tagged.bam + phased.bed 產出
- [ ] 3 條 pipeline canonical run 完成
- [ ] Master dataset 聚合含新樣本
- [ ] KB 新樣本文件撰寫
- [ ] 00_index 與 truth_set_registry 更新
- [ ] KB 驗證腳本全過
- [ ] git commit 說明新樣本加入

---

## 跨樣本一致性驗證

新樣本加入後，跑跨樣本一致性：
```bash
python3 scripts/analysis/observe_fisher_signal_stability.py \
  --include-samples HCC1395 COLO829 ... <NEW_SAMPLE>
```

**預期**：新樣本的 direction consistency 應與既有 7 樣本相容（參考 HPFineNGroups、LOH×AF 等 positive findings）

---

## 常見 Pitfall

- ❌ **忘記登記 truth_total**：F1 計算會錯
- ❌ **VCF symlink 指向錯的 caller**：歷史 TO 錯誤案例（2026-04-04）
- ❌ **Platform 混淆**：ONT 5kHz vs Dorado vs R10 各有 caveat
- ❌ **沒建 HC BED**（如 COLO829）：需用 `num_callers≥2` 邏輯替代

---

## 相關

- 樣本索引：[../02_samples/00_index.md](../02_samples/00_index.md)
- Truth set 登記：[../08_truth_and_benchmark/01_truth_set_registry.md](../08_truth_and_benchmark/01_truth_set_registry.md)
- Canonical 結構：[../05_data_formats/05_output_directory_layout.md](../05_data_formats/05_output_directory_layout.md)
