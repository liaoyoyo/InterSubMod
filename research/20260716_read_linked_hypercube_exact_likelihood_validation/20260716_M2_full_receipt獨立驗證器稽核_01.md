<!--
建立時間：2026-07-16
目標：為 M2 chr1–22 × 7 technical datasets 建立與 production aggregation code 解耦的 full receipt verifier
處理範圍：extraction schema 1.2、ranking schema 2.0、canonical candidate table；不讀 BAM；由持久化 pattern/state/π 獨立重算 profile likelihood 與 certificate
關聯檔案：InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py
-->

# M2 full receipt 獨立驗證器稽核

## TL;DR

已完成一個**不匯入 production runner、production aggregator、ranker 或 solver**的獨立驗證器。它會由 154 份 extraction child receipts、154 份 ranking child receipts 與其表格重新建立 task index、加總核心數字、驗證單位守恆，逐列由 child candidate tables 重建 canonical candidate table，並由 R/A/X、fixed BQ、count、candidate states 與 persisted π 獨立重算 profile likelihood、simplex residual、Frank–Wolfe/KKT gap、relative weight 與 winner/tie partition。合成測試目前 23/23 PASS；因 full M2 尚未完成，本文件只能判定「驗證器實作 PASS、full154 實際資料驗證 PENDING」，不能先宣稱最後數字通過。

服務目標：**G4 多樣本一致性與 reproducibility、G5 可被外部驗證**。

## 關鍵假設與限制

- Task Type：**B Comprehensive validation**。
- 正式 scope 固定為 7 technical datasets × chr1–22＝154 tasks；HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 technical datasets。
- 驗證器只讀已產生的 JSON/TSV.GZ 與 SHA-256 sidecar；**不讀 BAM**，因此不會與 extraction I/O 工作競爭。
- 驗證器驗證的是「持久化結果是否完整且自洽」，不替代 BAM→table 的演算法正確性測試。
- Profile audit 是在 producer 保存的 π 上獨立重算 concave mixture objective 與 global-gap certificate；不重跑 SLSQP，也不表示已證明真實 clone 數或 parent edge。
- Production code 與 verifier 不共用 aggregation helper；schema 欄位名稱雖必須依相同公開契約重述，但計算迴圈與守恆判定是另一份實作。

## Step → Verify

1. 由 full extraction receipt 取得正式 scope，再從固定的 `samples/<dataset>/<chr>/receipt.json` 逐份重讀 child receipt；驗證 task key 必須是完整且不重複的 Cartesian product，每份 schema＝1.2.0、`all_pass=true`、所有 checks 為布林 true、receipt sidecar 與所有 output SHA/bytes 均吻合。
2. 從 154 份 extraction child receipts 重加總 counts、PS counts、component k distributions 與 max k；驗證重算 aggregate 必須與 `full_extraction_receipt.json/aggregate` 完全相等，max k 使用 `max`，不錯誤相加。
3. 由磁碟逐份讀取 154 份 ranking child receipts，並重新核對其 extraction input identities；驗證 schema＝2.0.0、dataset/chrom/PS-aware scope 正確，四個 extraction inputs 的 path/bytes/SHA 與相對應 extraction child receipt 完全相同。
4. 對 primary 與 minread sensitivity 的 global、per-dataset、basis×threshold cells 獨立重加總；驗證 numeric fields、categorical counters、partial-pattern funnels 與 full ranking aggregate 逐欄相等，每個 cell 另驗證 12 項單位守恆。
5. 從每份 child `m2_compressed_vertex_set_candidates.tsv.gz` 串流重建 canonical output row；驗證與 `m2_ps_aware_candidate_table.tsv.gz` 逐列完全相等，並重算 physical SHA、semantic SHA、row/unit 數、排序、candidate ID、winner partition 與 solver-complete unit coverage。
6. 將 child `m2_symbolic_pattern_counts.tsv.gz` 與 candidate table 依 primary unit lockstep streaming join；直接重算每個 candidate 的 BQ emission、profile LL、simplex residual、KKT gap、relative likelihood weight 與 winner/tie。Full audit 要求 profile units/candidates 分別等於 canonical candidate-table units/rows，且154份 child summaries 全部存在。

## 獨立守恆條件

每個 ranking aggregate cell 都重新檢查：

1. selection status 合計＝units；
2. unique＋tied＋abstain＝units；
3. solver complete＋incomplete/not-run＝units；
4. projections＝informative＋all-X；
5. informative＝structural retained＋below-minread；
6. BQ scoring molecules＝informative molecules；
7. raw tree candidates T ≥ distinct vertex sets V；
8. coarse topology unique＋multiple＝topology evaluated；
9. partial covered＋unsatisfied＝partial denominator；
10. partial unsatisfied＝0；
11. k-route counts 合計＝units；
12. component k＝structural observed-ALT k＋非 structural ALT-active k。

## 輸入、執行命令與輸出

### 本次合成測試

輸入：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_full_m2_receipts.py`

本次 identity：verifier SHA-256 `4859598d74486f4eba6e4af6fa2dec2b4c0eb5c4e8ed86feac82483e5a7f32d8`；test SHA-256 `ae35510a44235ebeed4ffe6e84b385850cb2a29912146f9b283112820bfafb59`。

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
/usr/bin/time -v python3 -m unittest -v \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_full_m2_receipts.py
```

實際輸出片段：

```text
Ran 23 tests in 0.042s
OK
Elapsed wall time: 0.33 s
Maximum resident set size: 36,424 kB
Exit status: 0
```

### full154 完成後的正式命令

輸入路徑：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_schema_1_2/`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_ranking_schema_2_0/`

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
OPENBLAS_NUM_THREADS=1 /usr/bin/time -v python3 \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py \
  --extraction-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_schema_1_2 \
  --ranking-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_ranking_schema_2_0 \
  --output /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_independent_verification.json
```

輸出路徑：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_independent_verification.json`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_independent_verification.json.sha256`

預期成功輸出片段（**尚未執行 full154，因此不是實際結果**）：

```json
{"all_pass": true, "failure": null}
```

## 合成測試涵蓋的失敗模式

- task 數量正確但同一 dataset×chrom 重複，另一 task 遺漏；
- extraction filter funnel 或 molecule count 不守恆；
- component counts 應加總、max k 卻被錯誤相加；
- ranking full aggregate 的核心欄位被竄改；
- selection status 與 units 不守恆；
- full candidate table 的 metadata、physical SHA 與 semantic SHA 都自洽，但某一列與 child source 不一致。
- method contract 漂移，或 actual ranker source bytes 與 receipt 綁定不一致；
- 額外 LL／類 VAF additive score、BQ、π、state、KKT gap、winner/tie 任一被篡改；
- profile checks 被設為 false 但 full receipt 偽造 `all_pass=true`；
- 不同 state sets 的真正 likelihood tie（positive control）被錯分成唯一第一。

## 當前 verdict

| 層級 | 結果 | 可下的結論 |
|---|---:|---|
| 程式語法 | PASS | `py_compile` 成功 |
| 獨立單元／合成測試 | 23/23 PASS | 已能捕捉聚合、task、守恆、候選漏改、method identity與profile likelihood篡改風險 |
| Full154 真實 extraction/ranking | PENDING | full receipts 尚未完成，沒有最後數字可驗證 |
| 最終對教授報告 | NOT READY | 必須以正式 verifier receipt 的 `all_pass=true` 作為 final gate 之一 |
