<!--
建立時間: 2026-03-12 15:15
目標: 盤點 HCC1395 5kHz / HCC1395_DORADO 在 tumor-only 下的雙模型可得性、benchmark scope 與 snapshot scope，確認目前哪些結論可下、哪些只是資料不存在
處理範圍:
  - TO pileup / full model source tree inventory
  - benchmark-test availability
  - full tagged BAM vs candidate-window subset tagged BAM scope 差異
  - Pool B FN 與主線整合缺口盤點
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md
-->

# TO 雙模型與 snapshot scope 盤點

## 1. 破題結論

這一輪盤點的最重要結論不是新的 F1，而是把目前 `TO` 主線中 3 個常被混在一起的問題拆清楚：

1. **TO 的 `pileup vs full model` 目前不是「只差 benchmark」，而是現有 source tree 中尚未看到 full model 產物。**
2. **`5kHz TO` 與 `DORADO TO` 的 snapshot scope 差異已被證實，不是推測，因此 read-level 圖目前只能做方向比較。**
3. **`Pool B FN` 已有足夠完整的 caller-side rescue 結論，但還沒有正式接回四象限整合與 annotation 總結。**

因此現階段最合理的口徑是：

- `TO dual-model comparison`: **unavailable in current source tree**
- `TO read-level comparison`: **directional-only**
- `Pool B FN`: **有結論，但尚未整合收尾**

---

## 2. 本輪研究問題

### Q1. 目前 `HCC1395 5kHz TO` 與 `HCC1395_DORADO TO` 是否有可直接比較的 `pileup / full model`？

### Q2. `5kHz TO full tagged BAM` 與 `DORADO TO candidate-window subset tagged BAM` 的差異會不會影響 read-level diagnostics 解讀？

### Q3. `Pool B FN` 是否已經足夠成熟，可視為 caller-side rescue 的已成立結論？

---

## 3. 盤點方法

### 3.1 固定檢查路徑

#### TO caller 來源樹

- `/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0`
- `/big8_disk/data/HCC1395/ONT/ClairS_TO_ss_v0_3_0`
- `/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0`
- `/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_ss_v0_3_0`

#### paired raw caller 對照樹

- `/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0`

#### TO diagnostics 與 candidate-specific 輸出

- `/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1`
- `/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue`

### 3.2 固定檢查內容

1. `run_clairs_to.log` 或等價 log 中是否出現：
   - `PILEUP`
   - `FULL-ALIGNMENT`
   - `fa_tensor_can`
   - `full_alignment.vcf`
2. `benchmark-test/` 是否存在
3. snapshot 用 BAM 是 `full tagged BAM` 還是 `candidate-window subset tagged BAM`
4. 與主線文件中的既有敘述是否一致

---

## 4. TO 模型可得性盤點

### 4.1 source tree 盤點結果

| Dataset | Source dir | model family | log 證據 | full-model 證據 | benchmark-test |
| --- | --- | --- | --- | --- | --- |
| `HCC1395 ONT TO v0_3_0` | `/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0` | pileup-only | `SNV PILEUP AFFIRMATIVE MODEL PATH` | 無 | 有 |
| `HCC1395 ONT TO ss_v0_3_0` | `/big8_disk/data/HCC1395/ONT/ClairS_TO_ss_v0_3_0` | pileup-only | `SNV PILEUP AFFIRMATIVE MODEL PATH` | 無 | 無 |
| `HCC1395_DORADO TO v0_3_0` | `/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0` | pileup-only | `SNV PILEUP AFFIRMATIVE MODEL PATH` | 無 | 無 |
| `HCC1395_DORADO TO ss_v0_3_0` | `/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_ss_v0_3_0` | pileup-only | `SNV PILEUP AFFIRMATIVE MODEL PATH` | 無 | 無 |

### 4.2 paired raw caller 對照

paired raw caller 在 `/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/run_clairs.log` 中可明確看到：

- `PILEUP MODEL PATH`
- `FULL-ALIGNMENT MODEL PATH`
- `fa_tensor_can`
- `full_alignment.vcf`
- `full_alignment_filter.vcf`

這表示：

- paired raw caller 的 `pileup/full` 對照是 **真實存在**
- TO 目前沒有看到等價的 full-model 產物

### 4.3 盤點推論

目前最保守且正確的結論是：

1. `TO pileup vs full model` 在現有 source tree 中 **尚不可直接比較**
2. `v0_3_0` 與 `ss_v0_3_0` 應先視為兩種 **pileup model family**
3. 若未來要做 TO 雙模型 benchmark，前提是：
   - 找到或重跑真正的 TO full model 產物
   - 補 benchmark split

---

## 5. Benchmark scope 盤點

| Dataset | benchmark-test exists | 可直接做 baseline compare | 備註 |
| --- | --- | --- | --- |
| `HCC1395 ONT TO v0_3_0` | 是 | 可以 | 有 `tp.vcf / fp.vcf / fn.vcf` |
| `HCC1395 ONT TO ss_v0_3_0` | 否 | 不行 | 缺 compare step |
| `HCC1395_DORADO TO v0_3_0` | 否 | 不行 | 缺 compare step |
| `HCC1395_DORADO TO ss_v0_3_0` | 否 | 不行 | 缺 compare step |

這張表的意義是：

- 現在即使把 TO model availability 問題放開，也不能直接做四條 TO 路線的 benchmark 比較。
- `benchmark-test` 本身就是第二層可比性門檻。

---

## 6. Snapshot scope 盤點

### 6.1 已證實的 BAM 母體差異

| Dataset | snapshot BAM | scope | 用途 |
| --- | --- | --- | --- |
| `HCC1395 5kHz TO` | `/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam` | full tagged BAM | 主 pilot read-level diagnostics |
| `HCC1395_DORADO TO` | `/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam` | candidate-window subset tagged BAM | candidate-specific low-cost diagnostics |

### 6.2 這代表什麼

目前 read-level 圖的比較方式應明確分成兩層：

#### 可以做

- 同 dataset 內的 TP vs FP 方向比較
- 同 dataset 內不同特徵的代表性現象比較
- `5kHz TO` 與 `DORADO TO` 的現象類型比較

#### 不可以做

- `5kHz TO` 與 `DORADO TO` 的 `depth / hp_assign_rate / target_alt_fraction` 絕對值硬比
- 直接把兩者放在同一個 read-level 表格裡，假設分布母體相同

### 6.3 盤點推論

snapshot scope 差異 **已辨識，但尚未消除**。  
若要正式關閉這個限制，下一步應做：

1. `5kHz TO` 建一份 `candidate-window subset snapshot`
2. 在相同 selected regions 上比較 full vs subset
3. 若偏移小，再判定 subset 是否可作 DORADO TO 的近似 read-level 方案

---

## 7. Pool B FN 盤點

### 7.1 已有定論

來自 [20260308_ClairS邊緣FN探勘與甲基救援_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md) 的已成立結論：

| 項目 | 數值 / 結論 |
| --- | --- |
| `Pool B FN` | `840` |
| 甲基覆蓋率 | `99.9% (839/840)` |
| 最佳 caller-only | `no_varcluster` |
| 最乾淨 caller-only | `no_varcluster_and_gq15` |
| 最佳 combined | `gq15_and_allele_delta_low`，但只比 caller-only 多 `+0.000386` |
| 甲基唯一規則 | 全數負效益 |
| 最終口徑 | caller-side rescue 成立；甲基不適合作獨立主規則 |

### 7.2 目前還缺什麼

缺的不是新分析，而是整合收尾：

1. 尚未接回四象限整合總表
2. 尚未接回 phase 2 / annotation layer 總結
3. 尚未在主線週報中以 closing narrative 明確寫出

---

## 8. 本輪可直接下的結論

### 已成立

1. `TO full model` 在現有 source tree 中 **尚未確認存在**
2. `v0_3_0` 與 `ss_v0_3_0` 目前不應被寫成 `pileup vs full model`
3. `5kHz TO` 與 `DORADO TO` 的 snapshot scope 差異是真實存在的
4. `Pool B FN` 目前已足以支持 caller-side rescue 空間成立

### 尚不能下的結論

1. `TO pileup vs full model` 哪個更好
2. `5kHz TO` 與 `DORADO TO` read-level 指標的絕對值誰更高、差多少
3. `Pool B FN` 是否已與主線總報告完全閉環

---

## 9. 立即下一步

### 第一優先

把目前盤點結果正式升級成 closeout：

- `TO dual-model availability = unavailable in current source tree`

### 第二優先

做 `5kHz TO` 的 same-scope subset snapshot，驗證 full vs subset 的偏移。

### 第三優先

把 `Pool B FN` 接回：

- 四象限整合表
- phase 2 / annotation layer 總結
- 主線週報 closing section

---

## 10. 對 closing plan 的影響

本輪盤點後，`未完項 closing plan` 的優先順序可以正式定案為：

1. `TO dual-model availability closeout`
2. `snapshot scope same-scope control`
3. `Pool B FN integration closeout`

對應計畫見：

- [docs/plans/2026/03/20260312_未完項closing_plan_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md)

