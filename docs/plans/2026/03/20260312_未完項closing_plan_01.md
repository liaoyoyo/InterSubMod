<!--
建立時間: 2026-03-12 15:10
目標: 將目前尚未關閉的 TO 雙模型、snapshot scope 與 Pool B FN 整合問題整理成可執行的 closing plan
處理範圍:
  - TO pileup vs full model 可得性與可比性
  - 5kHz TO / DORADO TO snapshot scope 母體差異
  - Pool B FN caller-side rescue 與甲基訊號整合收尾
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md
-->

# 未完項 Closing Plan

> **更新註記（2026-03-12 16:20）**：本文件中的三條 closing workstream 已全部完成。本文目前保留作為規劃與完成脈絡紀錄；正式 closeout 結論請優先看：
> - [docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md)
> - [docs/reports/validated/2026/03/20260312_TO_snapshot_scope_same_scope_control_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO_snapshot_scope_same_scope_control_01.md)
> - [docs/reports/validated/2026/03/20260312_PoolB_FN_integration_closeout_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_PoolB_FN_integration_closeout_01.md)

## 破題結論

本文件原本定義的 3 條未完項，現在都已完成 closeout：

1. `TO pileup vs full model` 已 closeout 為：`unavailable in current source tree`，目前 TO caller 模型證據應統一標示為 `pileup-route only`。
2. `5kHz TO` 與 `DORADO TO` 的 snapshot scope 差異已補做 `same-scope control`；目前不可做跨 dataset 絕對值硬比的主因是 dataset/platform 差異，而不是 subset 技術本身。
3. `Pool B FN` 已完成 integration closeout；正式定案為 `caller-first` 主導、甲基只提供弱補強，且 `AlleleDelta` 不可跨母體全域化。

因此這份 closing plan 目前的用途，已從「待執行任務表」轉為「本輪三條收尾主線的完成記錄」。

若要延續後續工作，重點應改為：

1. 把 closeout 結論回接到更上游的 annotation / dashboard / deck narrative
2. 繼續做新的 feature 驗證，而不是再把這三條議題當成 blocker

---

## 1. 本次研究問題

### Q1. TO 下是否真的存在可直接比較的 `pileup vs full model`？

- 若存在，才能做直接 benchmark 對照。
- 若不存在，應正式記錄為資料缺口，後續需重跑或另找來源。

### Q2. `5kHz TO` 與 `DORADO TO` 的 snapshot scope 差異會不會影響 read-level 解讀？

- 若差異足以改變 `alt_fraction / hp_assign_rate / HP skew` 的絕對值解讀，就不能再把兩者放在同一張 read-level 表裡做硬比較。

### Q3. `Pool B FN` 的結論是否已足夠成熟，可接回主線報告？

- 目前已有 `840` 個 Pool B FN、caller-side rescue 規則、甲基訊號限制。
- 但尚未整合進四象限總表、annotation layer 與 closing narrative。

---

## 2. 三條未完項總覽

| Workstream | 現況 | 主要 blocker | 成功條件 | 預期輸出 |
| --- | --- | --- | --- | --- |
| A. TO dual-model availability closeout | ✅ 已完成 | 現有 TO source tree 僅確認 pileup-route only | 正式回答現有資料沒有 full model 產物 | `20260312_TO雙模型可得性closeout_01.md` |
| B. Snapshot scope same-scope control | ✅ 已完成 | `5kHz TO` same-scope control 結果 `18/18 identical` | 已區分 subset 技術偏移與 dataset/platform 差異 | `20260312_TO_snapshot_scope_same_scope_control_01.md` |
| C. Pool B FN integration closeout | ✅ 已完成 | authoritative hierarchy 與主線 narrative 已收斂 | caller-side rescue 與甲基限制已寫回主線整合表 | `20260312_PoolB_FN_integration_closeout_01.md` |

---

## 3. 固定研究模板

### 3.1 本次假設

1. `TO full model` 在目前資料樹中大機率不存在；`v0_3_0` 與 `ss_v0_3_0` 應是兩種 pileup model family，而不是 pileup/full 的對照組。
2. snapshot 母體差異原本可能影響 `depth`、`target_alt_fraction`、`hp_assign_rate` 等 read-level 指標的絕對值；目前已完成 `5kHz TO` same-scope 控制，用來分離 `subset` 技術偏移與 dataset/platform 差異。
3. `Pool B FN` 的主救援訊號仍是 caller-side（特別是 `VariantCluster` / `GQ`），甲基訊號不適合作獨立主規則，但這個結論值得正式接回主線整合文件。

### 3.2 成功條件

1. 能用文件與資料樹證明 `TO full model` 是否存在，而不是停留在推測。
2. 能把 `5kHz TO` / `DORADO TO` 的 read-level diagnostics 限制轉成一個已驗證的 `same-scope` 結論，並回寫主線文件。
3. 能在主線總表中明確寫出：
   - `Pool B FN = 840`
   - 最佳 caller-side rescue
   - 甲基訊號不適合作獨立主規則
   - 為何這件事不應再停留在支線文件

### 3.3 失敗條件

1. 無法證明 TO 是否有 full model 產物，仍用模糊語言寫成「可能有」。
2. 知道 snapshot 母體不同，但沒有完成 same-scope 控制或沒有把結果寫回主線文件。
3. `Pool B FN` 仍只存在單一 in-progress 文件，未接回主線整合表。

### 3.4 評估指標

- `availability_confirmed`: 是否已明確確認 TO full model 產物存在與否
- `comparability_status`: `direct / directional-only / unavailable`
- `same_scope_ready`: 是否已有可執行 same-scope 對照設計
- `same_scope_closed`: 是否已完成 same-scope 控制並產出 closeout 結論
- `pool_b_integrated`: Pool B FN 是否已接回四象限整合 / 週報 / annotation narrative

---

## 4. Workstream A：TO dual-model availability closeout

### 4.1 現有證據

根據目前已查到的資料樹與 log：

- `HCC1395 ONT TO v0_3_0`
- `HCC1395 ONT TO ss_v0_3_0`
- `HCC1395_DORADO TO v0_3_0`
- `HCC1395_DORADO TO ss_v0_3_0`

都只明確出現：

- `SNV PILEUP AFFIRMATIVE MODEL PATH`
- `INDEL PILEUP AFFIRMATIVE MODEL PATH`

尚未看到 paired raw caller 中會出現的：

- `FULL-ALIGNMENT MODEL PATH`
- `fa_tensor_can`
- `full_alignment.vcf`
- `full_alignment_filter.vcf`

### 4.2 這條主線的正確問題

不是「TO pileup vs full model 要怎麼比」，而是：

> **現有 TO 來源到底有沒有 full model 可比產物？**

### 4.3 執行任務

1. 固定盤點 4 個 TO 來源目錄的 log 與中間產物。
2. 固定比對 paired raw caller 的 full-model 證據欄位。
3. 補齊 benchmark 可得性表：
   - 哪些 TO 路線已有 `benchmark-test`
   - 哪些沒有
4. 若無 full model：
   - 正式 closeout 為 `unavailable in current source tree`
   - 不再把 TO 雙模型比較列為「差一步就能跑」的任務
5. 若有 full model：
   - 才進入 phase 2 benchmark 對照

### 4.4 預計輸出

- `docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md`
- `docs/experiments/in_progress/2026/03/assets/20260312_to_model_scope_inventory_01.tsv`
- 後續若 closeout 完成，再補：
  - `docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md`

### 4.5 預期結論形式

- `TO dual-model comparison unavailable in current source tree`
- `v0_3_0 vs ss_v0_3_0 = different pileup model family, not pileup/full contrast`

---

## 5. Workstream B：Snapshot scope same-scope control

### 5.1 現有證據

- `5kHz TO` diagnostics 使用 `full tagged BAM`
- `DORADO TO` diagnostics 使用 `candidate-window subset tagged BAM`

因此在控制前只能說：

- read-level 圖可做方向比較
- 不可做絕對值硬比

### 5.2 需要消除的限制

目前最關鍵的不確定性不是「圖看起來像不像」，而是：

- `subset tagged BAM` 是否會系統性提高或降低 `target_alt_fraction`
- 是否會改變 `hp_assign_rate`
- 是否會讓 `depth / HP skew / NA HP fraction` 的分布與 full tagged BAM 不同

### 5.3 執行任務

優先採低成本路徑：

1. 對 `5kHz TO` 建立 `candidate-window subset snapshot`
2. 在相同的 18 個代表性 region 上比較：
   - `full tagged BAM`
   - `candidate-window subset tagged BAM`
3. 比對下列指標是否系統性偏移：
   - `depth`
   - `target_alt_fraction`
   - `na_hp_fraction`
   - `collapsed_hp_balance_delta`
4. 若偏移小且方向一致：
   - 正式 closeout 為 `subset-safe-for-within-dataset-ranking`
5. 若偏移大：
   - `DORADO TO` 必須補 full tagged BAM 或至少 full-BAM snapshot 驗證

### 5.4 預計輸出

- `same_scope_snapshot_comparison.tsv`
- `same_scope_snapshot_bias_summary.md`
- `docs/reports/validated/2026/03/20260312_TO_snapshot_scope_same_scope_control_01.md`
- 更新既有：
  - `snapshot_scope_bias_assessment.tsv`
  - `TO support 特徵 read-level diagnostics` 的限制說明

### 5.5 預期結論形式

- `subset-safe-for-within-dataset-ranking`
- `unsafe-for-cross-dataset-absolute-comparison`
- `no-observed-subset-bias-in-5khz-control`

---

## 6. Workstream C：Pool B FN integration closeout

### 6.1 現有證據

`Pool B FN` 已有明確數據：

- `Pool B FN = 840`
- 甲基覆蓋率補跑後為 `99.9%`
- caller-only 最佳規則：
  - `no_varcluster`
  - `no_varcluster_and_gq15`
- 甲基唯一規則為負
- caller+methyl 最佳只比 caller-only 多 `+0.000386` 的 F1 delta

### 6.2 這條主線還差什麼

不是再做一輪 rescue，而是把已成立的結論接回：

1. 四象限整合總表
2. phase 2 / annotation layer narrative
3. 主線週報 / closing summary

### 6.3 執行任務

1. 把 Pool B FN 統一整理成：
   - baseline 母數
   - candidate-specific coverage
   - caller-only 最佳規則
   - combined 規則最終定位
2. 寫進四象限整合報告的限制與補充章節。
3. 寫進 phase 2 / annotation layer 總結：
   - `Pool B FN` 再次支持 `caller-first`
   - 甲基不適合作獨立主 rescue 規則
4. 補入 closing summary 與 next-step matrix。

### 6.4 預計輸出

- `pool_b_fn_integration_summary.tsv`
- `docs/reports/validated/2026/03/20260312_PoolB_FN整合收尾報告_01.md`
- 更新：
  - `20260311_HCC1395_四象限甲基rescue整合報告_01.md`
  - `20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md`

### 6.5 預期結論形式

- `Pool B FN confirms caller-side rescue space`
- `Methylation not suitable as independent primary rescue rule in Pool B`

---

## 7. 執行順序與依賴關係

### 第一優先：A. TO dual-model availability closeout

理由：

- 這會直接決定後面是補資料、補 benchmark，還是正式關閉「TO 雙模型比較」。

### 第二優先：B. Snapshot scope same-scope control / closeout

理由：

- 目前 `5kHz TO` same-scope control 已完成，下一步是把結果回寫主線，避免文件仍停留在舊限制。

### 第三優先：C. Pool B FN integration closeout

理由：

- 本條已完成；其作用是補齊主線 narrative，而不是再另開新探索。

---

## 8. 本輪已開始的研究動作

本輪先完成 `round 0 inventory`，已確認：

1. TO 現有來源目錄中，尚未找到 full model 產物。
2. `v0_3_0` 與 `ss_v0_3_0` 目前應視為不同 pileup family，而非 pileup/full 對照。
3. `5kHz TO full tagged BAM` 與 `DORADO TO subset tagged BAM` 的 scope 差異已被文件證實，不是推測；而 `5kHz TO` same-scope control 已證明 subset 技術本身不造成指標偏移。
4. `Pool B FN` 的方向已夠清楚，未來工作應以整合收尾為主。

詳細盤點結果見：

- [docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260312_TO雙模型與snapshot_scope盤點_01.md)

---

## 9. 下一步建議

### 立即要做

1. 把三份 closeout 文件維持與 `CURRENT_FOCUS / INDEX / 四象限整合報告` 一致
2. 若後續還要深化，優先接到 annotation / dashboard / deck narrative
3. 新研究 round 不再把這三條議題當作 blocker

### 先不要做

1. 不要在沒有 full model 證據的前提下繼續寫 TO `pileup vs full` 差異結論
2. 不要在 scope 不同的前提下把 `5kHz TO` 與 `DORADO TO` 的 read-level 絕對值放在同一張比較表硬比
3. 不要再把 `Pool B FN` 當作「尚未有方向」的開放問題
