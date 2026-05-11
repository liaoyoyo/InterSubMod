<!--
build_date: 2026-04-30 04:30
status: in_progress
audience: PI / 工程深度查證
parent_plan: InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_plan_01.md
parent_qa: InterSubMod/docs/reports/validated/2026/04/20260429_supplement_getVote_design_intent_QA_01.md
verdict: pending
-->

# V3-Fixed Ablation @ 0.6 Purity — 確認 PON 方法是否在低純度仍最佳 + F1 評估

## §0 一句結論（待 B3 完成填）

> **caller F1 已實證**：0.93 = **0.7166**；0.6 = **0.6273**（首次發布；ClairS-TO 自身性能隨 purity 下降 −12%）。
>
> 5-BAM ablation @ 0.6 待 B3 完成。

---

## §1 動機

- **0.93 純樣本 ablation 已完成**：A1-A5 5-BAM 對比；證實 PON-only flag 是 ratio 翻轉主因（4.3× 差距）
- **0.6 purity 同矩陣未驗證**：09_purity06_simulation 只測 baseline + V5，缺 V2b / V3F+PONonly / V3F-no-pononly
- **F1 從未在 0.6 計算**（09_purity06_simulation §6 標 caveat：F1 vs SEQC2 truth 仍未計算）
- 用戶要求：**確認 PON 的方法在 0.6 下仍比較好，特別 F1**

---

## §2 5-BAM 矩陣 @ 0.6 purity（HCC1395 t30_n20）

| 實驗代號 | binary | --pon-only-phasing flag | BAM 路徑 | 狀態 |
|---|---|:---:|---|:---:|
| **B1** baseline @ 0.6 | longphase-to-baseline | OFF | `output/purity_06_simulation/baseline_06/tumor_tagged.bam` | ✅ 已有 |
| **B2** V2b @ 0.6 | longphase-to-baseline | ON | （未做，本實驗暫不跑）| ⏸ |
| **B3** v3f_no_pononly @ 0.6 ★ | longphase-to-v3fixed | OFF | `output/purity_06_simulation/v3f_no_pononly_06/tumor_tagged.bam` | ⏳ **執行中**（btty29dzu）|
| **B4** V3F + PON-only @ 0.6 | longphase-to-v3fixed | ON | （未做，本實驗暫不跑）| ⏸ |
| **B5** V5 @ 0.6 | longphase-to (V5) | ON | `output/purity_06_simulation/v5_06/tumor_tagged.bam` | ✅ 已有 |

⭐ **本實驗只跑 B3**（最關鍵，與 0.93 A3 對應）。**3-BAM 對比 B1 + B3 + B5** 即可隔離 PON-only 效果。

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/longphase-to-mod/

./longphase-to-v3fixed phase \
  -s /big8_disk/data/HCC1395/ONT/subsample/t30_n20/ClairS_TO_v0_3_0/snv.vcf.gz \
  -b /big8_disk/data/HCC1395/ONT/subsample/t30_n20/HCC1395_t30_n20.bam \
  -r /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
  -t 48 \
  --caller clairs_to_ssrs \
  --pon-file ... \
  --strict-pon-file ... \
  --ont \
  --loh \
  -o output/purity_06_simulation/v3f_no_pononly_06/tumor_phased
  # ⚠ 故意不加 --pon-only-phasing flag
```

---

## §3 F1 評估

### §3.1 caller F1（不依賴 phasing/tagging）— 已完成

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/caller_f1.tsv`

| 樣本 | TP | FP | FN | Precision | Recall | **F1** |
|------|---:|---:|---:|----------:|-------:|-------:|
| ClairS-TO @ 0.93 raw | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| **ClairS-TO @ 0.6 raw** | **24,190** | **13,487** | **15,257** | **0.6420** | **0.6132** | **0.6273** |

**Δ F1 (0.93 → 0.6) = -0.0893**（caller 自身在低 purity 下顯著下降）

**驗證**：0.93 F1=0.7166 與 V3-Fixed memory + PI 報告 4 §3.4 完全一致。

### §3.2 V5/baseline 對 caller F1 的影響

從 PI 報告 4 §3.4：
- V5 不改 ClairS-TO caller → **Raw F1 對所有版本完全相同**（0.93=0.7166）
- 預期 0.6 也類似：**所有版本 Raw F1 = 0.6273**（caller 已固定）

### §3.3 ISM SuggestFilter F1（屬另案，本實驗暫不跑）

- 需跑 ISM 流程 5 BAM × ~30 min = ~2.5 hr
- 從 0.93 數據（PI 報告 4 §3.4）：
  - Baseline+ISM = 0.7157
  - V5+ISM = 0.7154（vs Baseline+ISM = -0.0003 噪音）
- 預期 0.6 也類似：ISM SuggestFilter 對 F1 影響在噪音範圍

→ 本實驗暫不跑 ISM SuggestFilter F1（已知結論：F1 不衡量 V5 品質，QA Q12）。

---

## §4 預期結果（hypothesis）

基於 0.93 ablation 推論：

| Metric | B1 baseline_06 | **B3 v3f_no_pononly_06** | B5 v5_06 |
|--------|:---:|:---:|:---:|
| HP1:HP2 ratio @ chr19-22 | **1:1.14**（已知，從 09 章）| **預期 ≈ 1:1.14**（PON-only OFF）| 1:1.78（已知）|
| HP:i:33 reads | 0（enum bug）| **預期 >0**（M3 修正）| 7,141 |
| AMB% | 2.0% | 預期 ~7-8% | 12.4% |
| Phase block N50 | 接近 baseline | 預期 ≈ B1 | V2b 結構 |
| caller F1 | 0.6273 | 0.6273 | 0.6273 |
| ISM filter F1 | 0.6273 ± 0.001 | 0.6273 ± 0.001 | 0.6273 ± 0.001（ISM 不改 caller F1）|

### 預期判定

- **B3 vs B1 ratio 接近** → 「PON-only OFF 下 0.6 purity 自然就接近平衡，getVote 修補單獨無法翻轉」
- **B3 vs B5 ratio 差距小** → 「0.6 下 V5 修補幅度有限（baseline 已自然平衡）」
- **F1 全部接近 0.6273** → 「PON-only 不影響 caller F1（與 0.93 結論一致）」

---

## §5 限制

| 限制 | 說明 |
|------|------|
| 本實驗只跑 B3 | B2/B4 未跑（工作量太大）；3-BAM 對比已可隔離 PON-only flag 效果 |
| ISM SuggestFilter F1 未跑 | 屬另案；既有 0.93 數據已證 F1 不衡量 V5 品質 |
| 僅 HCC1395 一樣本 | 跨樣本驗證屬另案 |

---

## §6 後續

| 行動 | 優先 |
|------|---|
| 等 B3 phase + haplotag 完成 | **進行中** |
| 跑 3-BAM HP metrics（B1/B3/B5）| 高 |
| 寫 0.6 結果報告（合併 0.93 ablation 結論）| 高 |
| 補 B2/B4（如需完整 5-BAM）| 中 |
| 跑 ISM SuggestFilter F1（如需）| 低 |

---

## §7 落點

| 檔案 | 路徑 |
|------|------|
| 本計畫 | `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_plan_01.md` |
| 結果報告（待） | `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md` |
| Caller F1 數據 | `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/caller_f1.tsv` |
| B3 BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v3f_no_pononly_06/tumor_tagged.bam`（待產生）|
