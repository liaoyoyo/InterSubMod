<!--
build_date: 2026-04-30 00:00
status: in_progress（執行中，預估 ~70 min phase + haplotag + ~10 min metrics）
audience: PI / 工程深度查證
parent_qa: InterSubMod/docs/reports/validated/2026/04/20260429_supplement_getVote_design_intent_QA_01.md
inputs:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/longphase-to-v3fixed (V3-Fixed binary, commit 41ff147)
  - /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz
  - /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam (HCC1395 5kHz 純 tumor)
  - 4 個 PON databases
outputs:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_tagged.bam
  - InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_results_01.md (執行後)
verdict: pending
last_verified: 2026-04-30
-->

# V3-Fixed getVote-only Ablation Plan — 隔離「只修 getVote 不啟用 PON-only」效果

## §0 一句結論（待執行後填）

> 待 ablation 執行完成後填入：A3 vs A4 差距是否 <3pp（PON-only 不必要）或 >5pp（PON-only 必要）。

---

## §1 動機

從 `InterSubMod/docs/reports/validated/2026/04/20260429_supplement_getVote_design_intent_QA_01.md` Q15：

> 用戶質疑：V5 改善是否「**一定需要 PON-only phasing**」？還是「**僅修補 getVote 標記邏輯**」就足以解決 17.3:1 問題？

**現有數據邊界**：所有既有 BAM（V3-Fixed、V5）都在 PON-only flag 啟用下測試 — **無「baseline mode + V3-Fixed getVote」獨立對照組**。

→ 本 ablation 補上這個對照，嚴格隔離「getVote 修補」與「PON-only flag」兩個獨立變量。

---

## §2 4-cell ablation 設計

| 實驗代號 | PON-only flag | getVote 邏輯 | BAM 路徑 | 狀態 |
|---|:---:|:---:|---|:---:|
| **A1** baseline | OFF | 原始（priority bug + enum mismatch）| `output/baseline/tumor_tagged.bam` | ✅ 已有 |
| **A2** PON-only only | ON | 原始（V2b commit `8b8c1fd`）| `output/pononly_v2b/tumor_tagged.bam` | ✅ 已有 |
| **A3** getVote-only fix | **OFF** | **V3-Fixed (`41ff147`)** | `output/v3f_no_pononly/tumor_tagged.bam` | ⏳ **本實驗目標** |
| **A4** V3-Fixed + PON-only | ON | V3-Fixed | `output/pononly_v3_fixed/tumor_tagged.bam` | ✅ 已有 |
| A5 V5 完整 | ON | V3-Fixed + Layer 1.5 | `output/pononly_v5_somatic_fallback/tumor_tagged.bam` | ✅ 已有 |

### 隔離效果

| 比對 | 隔離 |
|------|------|
| A1 → A2 | PON-only flag 對 phase 階段的單獨效果 |
| **A1 → A3** | **getVote 修補對 tag 階段的單獨效果（本實驗目標）** |
| A1 → A4 | 兩者組合效果 |
| **A3 → A4** | **PON-only 是否在 getVote 修補下提供額外貢獻？（本實驗目標）** |
| A4 → A5 | Layer 1.5 增量效果 |

---

## §3 執行步驟

### Step 1: phase（不啟用 PON-only，故意省略 `--pon-only-phasing`）

```bash
cd /big7_disk/liaoyoyo2001/longphase-to-mod/

./longphase-to-v3fixed phase \
  -s /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz \
  -b /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam \
  -r /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
  -t 48 \
  --caller clairs_to_ssrs \
  --pon-file /big7_disk/liaoyoyo2001/data/PON/1000g-pon.sites.vcf.gz,/big7_disk/liaoyoyo2001/data/PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz \
  --strict-pon-file /big7_disk/liaoyoyo2001/data/PON/dbsnp.b138.non-somatic.sites.vcf.gz,/big7_disk/liaoyoyo2001/data/PON/gnomad.r2.1.af-ge-0.001.sites.vcf.gz \
  --ont \
  --loh \
  -o /big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_phased \
  > /big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/phase.log 2>&1
```

⚠ **故意不加 `--pon-only-phasing` flag**（這是本實驗的關鍵 ablation）。

預估時間：~45 min（與 baseline 2693s 接近）。

### Step 2: haplotag

```bash
./longphase-to-v3fixed haplotag \
  -s /big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_phased.vcf \
  -b /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam \
  -r /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
  -t 48 \
  -o /big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_tagged \
  --tagSupplementary \
  --log \
  > /big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/haplotag.log 2>&1
```

預估時間：~24 min。

### Step 3: 計算 metrics（5 BAM 對比）

對 A1/A2/A3/A4/A5 同時計算：

| Metric | 命令／來源 |
|---|---|
| HP family count（HP:i:1+11 vs 2+21）| `samtools view ... -d HP` 統計 |
| HP:i:33 count | 同上 |
| AMB%（HP:i:33 / total tagged）| 同上 |
| Phase block N50 | 從 phased VCF 統計 PS field |
| Phased rate | 從 phased VCF |
| Sanity check（4 項硬性檢查）| `InterSubMod/scripts/analysis/v5_sanity_paired_check.py` |
| Paired GT concordance（15 sites）| 同上（如有 paired BAM）|

---

## §4 預期結果（hypothesis）

基於機制推論（QA Q15 §G）：

| Metric | A1 baseline | A3 getVote-only | A4 V3F + PON-only | 預期判定 |
|--------|:---:|:---:|:---:|---|
| HP1:HP2 ratio | 17.3:1 | **3-5:1**（推估）| ~1:1 | A3 部分改善但不完全 |
| HP:i:33 出現 | 0 | **>0**（M3 修正）| 240K | A3 enum mismatch 修復 |
| AMB% | N/A | **5-10%**（推估）| 17.5% | A3 可能比 A4 略低（無 PON-only 過度保守）|
| Phase block N50 | 4,061 | ≈ 4,061 | 8,109 | A3 ≈ A1（PON-only flag 不啟用 → phase 結構不變）|
| Phased rate | 54.9% | ≈ 54.9% | 78.5% | 同上 |

### 結論判定門檻

| A3 vs A4 差距 | 結論 |
|---|---|
| HP1:HP2 / paired GT 差距 < 3 pp | **PON-only 不必要**，僅修 getVote 即足夠 |
| 差距 3-5 pp | PON-only 提供額外貢獻但非必要 |
| **差距 > 5 pp** | **PON-only 確實是必要前提**，符合機制推論 |

---

## §5 風險與限制

| 風險 | 緩解 |
|------|------|
| 磁碟空間不足（4.9T 可用，每 BAM ~150GB）| 監控；必要時清理舊實驗 BAM |
| Phase 執行時間超出預估 | background 執行；不阻塞主對話 |
| BAM 路徑不一致（big7 vs big8）| 用 big7（與 V3F 既有測試一致）|
| Working tree V5 修改影響編譯 | 不需重編，用既有 V3-Fixed binary |
| Paired GT 比對需 paired tumor BAM | 從 `/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam` 取 |

---

## §6 後續延伸（純樣本驗證後）

1. **0.6 purity 同矩陣**：用 t30_n20 BAM 跑同樣 4-cell ablation
2. **跨樣本驗證**：擴展到 HCC1937 / HCC1954 / COLO829（如必要）
3. **論文支持比對**：本 ablation 結果與 Lin 2022 / 陳鎮宇 2025 §4.3 對比

---

## §7 落點

| 檔案 | 路徑 |
|------|------|
| 本計畫 | `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_plan_01.md` |
| 結果報告 | `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_results_01.md` |
| BAM 輸出 | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/` |
| 對比表／圖 | `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/` |

---

## §8 變更歷史

| 日期 | 動作 |
|-----|------|
| 2026-04-30 00:00 | 計畫建立；環境驗證通過（V3-Fixed binary 已存在）|
| ⏳ 待執行 | Step 1 phase（~45 min）|
| ⏳ 待執行 | Step 2 haplotag（~24 min）|
| ⏳ 待執行 | Step 3 metrics + 結果報告 |
