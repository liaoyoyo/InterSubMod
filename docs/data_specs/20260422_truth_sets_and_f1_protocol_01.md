<!--
建立時間: 2026-04-22
目標: 集中 7 樣本 truth set 路徑 + F1 計算口徑，讓 AI session 經 CLAUDE.md 速查表即可定位
處理範圍: 7 樣本 benchmark truth VCF + F1 計算協議 + canonical 比對工具
關聯檔案:
  - scripts/pipeline/config.sh (sample 配置與 TRUTH_TOTAL)
  - scripts/pipeline/steps/03_filter_analysis.py (F1 公式實作)
  - scripts/pipeline/run_benchmark.sh (管線入口)
  - scripts/analysis/compare_benchmark_f1.py (跨 run 聚合)
  - /big8_disk/liaoyoyo2001/Knowledge/02_samples/colo829.md
  - /big8_disk/liaoyoyo2001/Knowledge/06_workflows/benchmark-workflow.md
狀態: active
最後驗證: 2026-04-22
doc_type: reference
-->

# Truth Set 答案與 F1 比對協議（7 樣本）

> **此文件用途**：讓任何 AI session 一眼定位「哪個樣本用哪個 truth、怎麼篩、F1 怎麼算」。
> **不重複知識庫與原始碼內容**，僅聚合與導航 — 權威事實在引用位置。

---

## 一、7 樣本 Truth Set 速查表

| Sample | Platform | Truth Set | TRUTH_VCF 路徑（`/big8_disk/data/...`） | HC BED | `TRUTH_TOTAL` |
|--------|----------|-----------|------------------------------------|--------|--------------|
| **HCC1395** | ONT_5kHz | SEQC2 v1.2.1 HC SNV | `HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz` | `HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed` | 39447 |
| HCC1395_DORADO | ONT_Dorado | SEQC2 v1.2.1 HC SNV | （同上） | （同上） | 39447 |
| **COLO829** | ONT_PAO | **NYGC** | `COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz` | **無公開 BED** | 41427 |
| H1437 | ONT | orthogonal-tools | `H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark_somatic-only.vcf.gz` | `H1437_orthogonal-tools-benchmark.bed` | 90129 |
| H2009 | ONT | orthogonal-tools | `H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark_somatic-only.vcf.gz` | `H2009_orthogonal-tools-benchmark.bed` | 168529 |
| HCC1937 | ONT | orthogonal-tools | `HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark_somatic-only.vcf.gz` | `HCC1937_orthogonal-tools-benchmark.bed` | 60691 |
| HCC1954 | ONT | orthogonal-tools | `HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark_somatic-only.vcf.gz` | `HCC1954_orthogonal-tools-benchmark.bed` | 26567 |

- 上表權威來源：`scripts/pipeline/config.sh:50-191`（由 `get_sample_config()` 輸出）
- `.tbi` 索引：所有 `TRUTH_VCF` 皆已備妥（COLO829 簡化版 `.tbi` 於 2026-03-30 補建，D2）
- `TRUTH_TOTAL` 定義：**整個 truth VCF 的變異數**，為 recall 分母；若改用更嚴格子集需重設此常數

---

## 二、COLO829 NYGC 專區（特殊：無 HC BED）

### 解答檔案位置

| 用途 | 路徑 | 備註 |
|------|------|------|
| **Canonical benchmark** | `/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz`（+ `.tbi`） | 標準 caller 流程 truth |
| 外部驗證簡化版 | `/big8_disk/data/colo829_nygc/colo829_truth.vcf.gz`（+ `.tbi`） | D5 非標準樹，僅做 short-read 對照 |
| NYGC Illumina tumor BAM | `/big8_disk/data/colo829_nygc/COLO829_NYGC_50x.bam` | 50× NovaSeq |
| NYGC Illumina normal BAM | `/big8_disk/data/colo829_nygc/COLO829_NYGC_BL_25x.bam` | 25× NovaSeq（COLO829BL） |

> D5 決策（2026-03-01，see `06_workflows/benchmark-workflow.md:49-50`）：`colo829_nygc/` 不納入標準 caller 測試樹；標準流程用 ONT_PAO。

### 三級篩選口徑（NYGC 以 caller-consensus 取代 HC BED）

| 信度層級 | 篩選條件 | 適用場景 |
|---------|---------|---------|
| 極高信度 | `FILTER=PASS` ∧ `HighConfidence` ∧ `num_callers >= 3` | 論文表格 |
| **Gold standard（預設）** | `FILTER=PASS` ∧ `HighConfidence` ∧ `num_callers >= 2` | 標準 benchmark |
| 寬鬆 | `FILTER=PASS` | 探索性分析 |

篩選欄位來源：NYGC annotated VCF 的 INFO 區含 `num_callers`、`called_by`、`HighConfidence` flag；完整 FILTER 類型（PASS / GRM / HighNormalVAF / LowNormalDP / LowTumorDP / LowTumorVAF / PON / PartOfMNV）見 `02_samples/colo829.md`。

### 使用範例

```bash
# 以 Gold standard 子集比對（推薦）
bcftools view -f PASS \
  -i 'INFO/HighConfidence=1 && INFO/num_callers>=2' \
  /big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz \
  -Oz -o colo829_truth_gold.vcf.gz
tabix -p vcf colo829_truth_gold.vcf.gz
```

---

## 三、F1 計算協議

### 3.1 本專案內部管線公式（canonical for 自動化 F1）

**實作位置**：`scripts/pipeline/steps/03_filter_analysis.py:229-234`

```python
def compute_metrics(tp, fp, truth_total):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / truth_total if truth_total > 0 else 0.0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"tp": tp, "fp": fp, "precision": round(precision, 4),
            "recall": round(recall, 4), "f1": round(f1, 4)}
```

**參數來源**：
- `tp` / `fp`：LongPhase-S step 01 產生的 `filtered_snv_tp.vcf.gz` / `filtered_snv_fp.vcf.gz`（分類參照 `TRUTH_VCF`；由 LongPhase-S 的 `somatic_haplotag --truth-vcf` 在 haplotag 同步執行）
- `truth_total`：從 `config.sh` 讀入的 `TRUTH_TOTAL` 常數（例：COLO829=41427、HCC1395=39447）
- `recall` 分母 = 整個 truth VCF 的變異數（**不扣 HC BED 之外的變異**；若要更嚴口徑需重設 `TRUTH_TOTAL`）

**輸出**：
- `metrics.json`：baseline / filtered precision/recall/F1
- `benchmark_comparison.tsv`：跨 run 聚合（由 `scripts/analysis/compare_benchmark_f1.py` 產生，欄位見檔首 `FIELDS`）

### 3.2 Canonical 跨工具比較（論文級報告用）

依 `06_workflows/benchmark-workflow.md` D5 凍結口徑：

| 方法 | 角色 | 適用場景 |
|------|------|---------|
| **`som.py`（hap.py 框架）** | **Canonical** | 跨工具公平比較、論文表格 |
| LongPhase-S built-in `--truth-vcf` | 輔助 | haplotag + benchmark 一步流，本專案 pipeline 內部即採此 |
| ClairS Docker `compare_vcf` | 輔助 | ClairS/ClairS-TO 生態內一致性檢查 |

⚠ **som.py 限制**：indel > 100 bp 可能不納入或誤分類；涉及 SV / 大片段 indel 時，須搭配 truvari 或 `bcftools isec` + manual curation，**不要只看 som.py 報表**。

som.py 完整命令與參數說明：見 `/big8_disk/liaoyoyo2001/Knowledge/06_workflows/benchmark-workflow.md`（第 111-122、126-150 行）。

### 3.3 TP/FP/FN 定義

| 分類 | 定義 |
|------|------|
| TP | caller 輸出且位於 truth VCF（套用 BED/篩選後） |
| FP | caller 輸出但不在 truth VCF |
| FN | truth VCF 有但 caller 未輸出 |

分類工具：LongPhase-S `somatic_haplotag` 的 `--truth-vcf` 旗標（產生 `filtered_snv_tp.vcf.gz` / `_fp.vcf.gz`）；FN 需從 truth VCF 與 caller VCF 的 `bcftools isec` 推得。

---

## 四、Pipeline 呼叫範例

### 4.1 COLO829 paired_full（ClairS full output）
```bash
./scripts/pipeline/run_benchmark.sh --mode s-pure --sample COLO829
```

### 4.2 COLO829 paired_pileup（ClairS pileup VCF）
```bash
./scripts/pipeline/run_benchmark.sh --mode s-pure-pileup --sample COLO829 \
  --vcf-source pileup
```

### 4.3 COLO829 tumor-only
```bash
./scripts/pipeline/run_benchmark.sh --mode to-pure --sample COLO829 \
  --vcf-source pileup
```

### 4.4 跨 run 聚合 F1 比較表
```bash
python3 scripts/analysis/compare_benchmark_f1.py \
  --run-dir <canonical_run_dir> \
  --output-tsv benchmark_comparison.tsv \
  --output-md benchmark_comparison.md
```

- `--dry-run` 可加在任何 `run_benchmark.sh` 呼叫上預覽命令（不執行）
- Mode 名稱映射：`s-pure → paired_full`、`s-pure-pileup → paired_pileup`、`to-pure → to_pileup`、`to-full → to_full`（見 `config.sh:200-209` 的 `canonical_mode_name()`）

---

## 五、權威來源索引（不重複，只連結）

### 原始碼
- `scripts/pipeline/config.sh:50-191` — 7 樣本 sample 配置（TRUTH_VCF / TRUTH_BED / TRUTH_TOTAL）
- `scripts/pipeline/steps/03_filter_analysis.py:229-234` — F1 公式實作
- `scripts/pipeline/run_benchmark.sh` — 管線入口（4 步驟：LongPhase-S → InterSubMod → Filter Analysis → Cleanup）
- `scripts/analysis/compare_benchmark_f1.py` — 跨 run 聚合報表

### 知識庫（`/big8_disk/liaoyoyo2001/Knowledge/`）
- `02_samples/colo829.md` — NYGC truth 細節、`num_callers` / `HighConfidence` 欄位、三級篩選
- `02_samples/hcc1395.md` — SEQC2 truth 對照（第一優先 benchmark 樣本）
- `06_workflows/benchmark-workflow.md` — F1 公式、三方法比較、som.py canonical、最佳實踐
- `04_databases/seqc2-truth-set.md` — SEQC2 truth set 詳細

### 決策紀錄
- **D5**（2026-03-01）：`colo829_nygc/` 不納入標準 caller 測試樹
  → `06_workflows/benchmark-workflow.md:49-50`
- **D2**（2026-03-30）：`colo829_truth.vcf.gz.tbi` 索引補建
  → `02_samples/colo829.md:109`

---

## 六、最佳實踐提醒

1. **一致區域限制**：有 `TRUTH_BED` 的樣本（除 COLO829 外全部）必先用 `bcftools view -R <bed>` 限制比較範圍
2. **PASS filter**：僅評估 caller 輸出中 `FILTER=PASS` 的變異
3. **SNV 與 INDEL 分開**：難度與 `TRUTH_TOTAL` 不同，不混算；本專案 `TRUTH_TOTAL` 現多為 **SNV only**
4. **版本記錄**：報告時註明 caller 版本、truth set 版本、`num_callers` 門檻（COLO829）
5. **多樣本驗證**：單一樣本的改進結論需跨 4+ 樣本一致（見 `/multi-sample-consistency` skill）
