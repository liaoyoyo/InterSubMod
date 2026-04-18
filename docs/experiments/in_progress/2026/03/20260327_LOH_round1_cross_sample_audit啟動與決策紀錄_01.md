<!--
建立時間: 2026-03-27 00:30
狀態: in_progress
目標: 固定 LOH Round 1 cross-sample audit 的樣本範圍、決策、輸入優先順序、可追溯輸出與中途需回報的判斷
處理範圍:
  - 7 個具甲基資料樣本的 paired / tumor-only Round 1 盤點範圍
  - LOH-like 分箱優先、非單一固定閾值的初始判讀策略
  - PS 缺口的暫行處理方式
  - same-locus paired-vs-TO compare 與 case panel 的 Round 1 納入方式
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260326_LOH盤點執行規格_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_methodology/2026/03/20260326_InterSubMod五目標研究願景與LOH先導觀察策略_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/output/synthesis/master_run_manifest.tsv
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
-->

# LOH Round 1 Cross-Sample Audit 啟動與決策紀錄

## 1. 文件定位

本文件是 `LOH Round 1` 的專用啟動與決策紀錄。

角色是把使用者已確認的方向正式寫死，避免之後：

1. 用錯樣本範圍
2. 把 exploratory bins 誤寫成固定結論閾值
3. 忘記 paired / TO availability 不完全對稱
4. 忘記 same-locus compare 是本輪必做
5. 忘記哪些地方只能先作 `summary audit`，不能過度解釋

通用規格仍以 [20260326_LOH盤點執行規格_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260326_LOH盤點執行規格_01.md) 為主；本文件只負責 Round 1 的具體口徑。

---

## 2. Round 1 研究問題

Round 1 要優先回答的不是「LOH 是否已能提升 F1」，而是：

1. 在所有具甲基資料樣本中，`HP family imbalance / LOH-like / Potential_LOH / PS-related limitation` 的分佈到底長什麼樣子
2. 這些現象在 paired 與 TO 是否一致、互斥、互補，或只存在於特定 sample / mode
3. 哪些觀察可以形成 `evidence panel`
4. 哪些觀察只夠形成 `diagnostic / risk note`
5. 哪些地方真的卡在 `PS` 缺口，值得後續補 `PS export`

---

## 3. Round 1 已確認的總決策

### 3.1 樣本範圍

Round 1 直接做全樣本 cross-sample audit，不採最小閉環起手。

本輪預計納入 7 個有甲基資料的樣本：

1. `HCC1395_HKU_5kHZ`
2. `HCC1395_DORADO`
3. `HCC1937`
4. `HCC1954`
5. `H1437`
6. `H2009`
7. `COLO829`

**名稱對齊規則**

1. `HCC1395_HKU_5kHZ` 在現有 manifest 與 canonical 路徑中，正式 sample key 以 `HCC1395` 表示，平台欄位為 `ONT_5kHz`
2. 其餘 sample 以現有 canonical / research round 的 sample key 為主

### 3.2 閾值策略

Round 1 **不先假設單一最終 LOH-like 閾值已成立**。

本輪固定採：

1. 先輸出完整分箱與連續分佈
2. `Potential_LOH <0.1 or >0.9` 只保留為 legacy reference / comparison baseline
3. 本輪最終重點是看分佈、分層與可判別性，不是先把 `<0.1 / >0.9` 升級成最終結論

### 3.3 PS 缺口策略

Round 1 先做 `summary audit`，不先改 C++ 主輸出。

但若同位點或 case 分析需要 `phase block` 判讀，允許：

1. 依 `phased VCF` 的 `PS` 欄位做 locus-level 補充查詢
2. 在 case report 中把該資訊明確標註為 `PS recovered from phased VCF`
3. 只有當這類需求大量出現且明顯卡住解釋時，才把 `PS export to reads.tsv` 升級成後續工程任務

### 3.4 case panel 強度

Round 1 case panel 固定只做 `4–8` 個代表案例。

代表案例需覆蓋：

1. paired TP LOH-like
2. paired FP LOH-like
3. TO TP LOH-like
4. TO FP LOH-like
5. paired / TO 同位點結論一致 case
6. paired / TO 同位點結論不一致 case
7. `HP0 / HP3` 高比例但不宜強解釋 case

### 3.5 same-locus compare 優先度

Round 1 **直接納入** same-locus compare，不延後到 Round 1.5。

本輪要同時盤：

1. 同樣本中 `paired` 與 `TO` 有交集的位點
2. 同樣本中只出現在 `paired` 的位點
3. 同樣本中只出現在 `TO` 的位點
4. 上述交集 / 非交集位點內，`TP / FP` 的分佈與判讀差異

---

## 4. 資料來源優先順序

Round 1 的資料來源不要求所有 sample 都先 canonical 化為完全對稱結構，但必須記錄清楚優先順序與 availability。

### 4.1 paired 優先順序

1. `canonical paired_full`
2. 若 `paired_full` 缺必要欄位或缺 LOH 相關 sidecar，再看 `paired_pileup`
3. 舊 rerun 若只是歷史版本，不作主分析，只作 provenance 備註

### 4.2 tumor-only 優先順序

1. 若 sample 已有 `canonical to_pileup`，優先使用 canonical
2. 若尚未 canonical 化，使用對應 `research_round TO pilot` 作為 Round 1 TO 來源
3. `input_manifest.tsv` 必須明確註記：
   - 來源是 `canonical` 還是 `research_round`
   - 哪些 sample 缺 paired / TO 任一側
   - 哪些 sample 雖有 TO round，但 sidecar 或 tagged BAM 不完整

### 4.3 paired / TO availability 原則

1. Round 1 的 paired 軸以 7 樣本全納入為目標
2. TO 軸以「所有現有具甲基 TO pilot / canonical 輸出」全納入為目標
3. 若某 sample / mode 最終無法納入，不直接從 summary 消失，而是要在 manifest 與 `open_questions.md` 記錄缺口原因

---

## 5. Round 1 固定分析層次

### 5.1 Cross-sample distribution audit

固定要做：

1. 全樣本 `hp_ratio_core` 分佈
2. 全樣本 `effective_hp_reads` 與 `hp_assign_rate`
3. 全樣本 `hp0_ratio / hp3_ratio`
4. 全樣本 `Potential_LOH` 與 exploratory bins
5. 全樣本 `VerificationClass × LOH_Subtype`

### 5.2 Cross-mode paired-vs-TO compare

固定要做：

1. 同 sample paired vs TO 的 LOH-like 占比
2. 同 sample paired vs TO 的 feature distribution 對照
3. 同 sample paired vs TO 的 same-locus overlap / non-overlap 對照
4. paired tag.bam 與 TO 結論在同位點上的一致 / 不一致情況

### 5.3 TP / FP 判斷層

固定要拆：

1. `TP / FP`
2. `paired / TO`
3. `same-locus overlap / paired-only / TO-only`
4. `LOH-like / non-LOH-like`

目的是避免只看總體比例而看不見：

1. 其實只是某個 sample 拉動
2. 其實只是某個 mode 拉動
3. 其實只是 overlap / non-overlap 子集合在作怪

### 5.4 PS 補充層

若某些結論需要判讀 `phase block`：

1. 先從 phased VCF 補查 `PS`
2. 再決定是否在 case report 升級為 `PS-aware unresolved / partially resolved`
3. 不因為補查到少量 `PS` 就把整個 summary 假裝成 fully phase-aware

---

## 6. Round 1 輸出要求

除通用 `LOH` 規格的固定輸出外，本輪另外強制：

1. `input_manifest.tsv` 要多一欄 `availability_status`
2. `input_manifest.tsv` 要多一欄 `source_priority`
3. `same_locus_compare.tsv` 要多一欄 `locus_overlap_class`
4. `same_locus_compare.tsv` 要多一欄 `paired_to_concordance`
5. `decision_ledger.tsv` 要多一欄 `decision_scope_round1`
6. `open_questions.md` 要分成：
   - `blocked by data availability`
   - `blocked by PS`
   - `blocked by unclear threshold`
   - `needs user confirmation`

---

## 7. Round 1 預設輸出落點

本輪建議輸出到 observation workspace，而不是直接寫成 final research round closeout：

```text
output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/
```

理由：

1. 這輪重點是 cross-sample diagnostics 與分佈盤點
2. 尚未在方法層做 intervention
3. 尚未到可以直接 claim 成效結論的階段

---

## 8. Round 1 尚未關閉的模糊點

以下問題先不阻止 Round 1 開始，但執行中若明顯影響判讀，必須回報：

1. `paired_full` 與 `paired_pileup` 是否會對 HP / LOH proxy 給出明顯不同分佈
2. 各 TO pilot 的 sidecar 是否足以支撐完全一致的 table schema
3. same-locus key 是否只用 `chr:pos:ref:alt`，或需額外處理 caller normalization 差異
4. 若 `PS` 補查比例過低，case panel 是否仍值得做 phase-aware 描述

這些問題若在中途變成真 blocker，應先寫進 `decision_ledger.tsv` 與 `open_questions.md`，再決定是否需要你進一步拍板。

---

## 9. 本輪成功條件

Round 1 視為成功，不以 `F1` 是否提升判定，而以以下條件為準：

1. 已完成 7 樣本 paired 軸盤點
2. 已完成所有目前可用 TO sample 的盤點與 availability 註記
3. 已有完整分箱與連續分佈輸出
4. 已完成 same-locus paired-vs-TO compare
5. 已完成 4–8 個代表案例
6. 已明確指出哪些結論是直接事實、哪些只是推論、哪些仍卡在 `PS` 或資料缺口
