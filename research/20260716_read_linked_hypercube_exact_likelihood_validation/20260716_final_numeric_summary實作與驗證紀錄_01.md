<!--
建立時間：2026-07-16
目標：紀錄 presentation-stage final numeric summary 的可推導欄位、證據綁定與測試結果
處理範圍：M2 full extraction/ranking/final independent verifier；不修改 exact-11 release 或 HTML builder
關聯檔案：InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py
-->

# Final numeric summary 實作與驗證紀錄

## TL;DR

已建立獨立的 presentation-stage 數字彙整器。它只接受通過且具有正式 release binding 的 final independent-verification receipt，重新驗證 terminal receipts、全部 child receipts、candidate table 與 runtime files 的 SHA-256，才輸出每 dataset 的 component/k、partial funnel、h*、T/V、unique/tied/abstain、tied×coarse-Topo、HP1/HP2、runtime 與明確分母比例。

本工作服務 G4（多樣本一致性與 reproducibility）與 G5（可由外部稽核的證據鏈）。

## Step → Verify

1. 驗證 final、extraction、ranking terminal receipts
   → 驗證：schema/version、`all_pass`、所有 checks、external SHA sidecar、正式 release eligibility 均通過，且 final receipt 的 path/SHA exact bind 兩份 terminal receipts。
2. 重新讀取 154 extraction child receipts
   → 驗證：每份 child sidecar、embedded child equality、final task-index SHA、alignment/count funnel、component k distribution 均守恆；加總後 exact 等於 terminal 與 final-verifier aggregate。
3. 重新讀取 154 ranking child receipts
   → 驗證：每份 child sidecar、final task-index SHA、partial funnel、unit/solver/outcome/T/V/k/topology partition 均守恆；加總後 exact 等於 terminal 與 final-verifier recomputation。
4. 串流 canonical candidate table
   → 驗證：physical SHA、semantic SHA、row/unit count、排序、candidate completeness、每 unit 的 h* 一致、`parent_choice_count >= 1`、unique/tied/optimizer-abstain partition，以及 candidate T/V/solver-complete/topology coverage exact 對回 ranking aggregate。
5. 串流 child runtime TSV
   → 驗證：每檔 SHA/size、欄位與 scope、segment ≤ total、child exact-nearest-rank summary、全資料 exact 等於 terminal 與 final verifier。
6. 產生 presentation JSON
   → 驗證：exclusive create、JSON `allow_nan=False`、SHA sidecar、0444 read-only；無可靠來源的欄位輸出 `null + reason`，不以 0 代替未知。

## 可可靠推導與不可推導

| 欄位 | 判定 | 證據與定義 |
|---|---|---|
| 每 dataset component/k | 可 | 逐 chromosome child 的 `component_summary_by_linkage_basis` 加總；保留 basis×threshold |
| partial funnel | 可 | ranking child 的 primary partial-pattern funnel；同時計數 covered/unsatisfied |
| T | 可 | solver-complete units 的 feasible parent-edge assignment 總數 |
| V | 可 | solver-complete units 的 distinct optimal vertex-set 總數；candidate table rows 應 exact 相等 |
| h* | 可 | candidate `vertex_roles` 中不含 `root` 或 `full-observed` 的 state 數；不是 hidden clone 數 |
| unique/tied/abstain | 可 | ranking aggregate；candidate table另驗證 solver-complete 子集 |
| tied×coarse-Topo | 可 | 所有 tied winner coarse-class union 大小為 1＝consistent，>1＝inconsistent |
| exact tied parent-edge topology | 不可 | canonical candidate table 只有 `parent_choice_count`，未展開全部 edges；輸出 null+reason |
| HP1/HP2 cellular clone pairing | 不可 | HP-tagged reads 無 cell-level homolog pairing；輸出 null+reason |
| 跨 thresholds 去重 biological regions | 不可 | 缺少 authenticated cross-threshold region identity table；輸出 null+reason |
| 另加 same-read VAF confirmation | 不可作獨立項 | conditional read-pattern likelihood 已含同一份 R/A evidence，再加會 double count |
| 正式 full-run peak RSS/CPU/disk envelope | 不可由科學 receipt 推導 | receipt 有每 unit monotonic runtime，沒有正式全程 process/resource telemetry；輸出 null+reason |

### Red-team 修正（P1）

初版只把 candidate row 數核對到 V，尚未把每列 `parent_choice_count` 加總核對到 raw T，卻使用了 T/V binding conserved 的敘述。現已補成逐 dataset × HP basis × threshold 的 exact T conservation；同時由 winning rows 重新推導 coarse topology partition、topology class counts 與 single-winning-parent-edge uniqueness，逐格對回 ranking aggregate。若任何一格不同，summary fail closed，不輸出報告數字。

「並列第一且 coarse topology 一致」採保守定義：所有 tied winning vertex sets 所含 coarse classes 的聯集只能有一類。即使每個 winner 都有同一組多類 ambiguity，仍不稱為 topology 一致；而且這個分類不主張 exact parent-edge topology 已唯一。

## 輸入、命令與輸出

正式執行前提是完整 7 technical datasets × chr1–chr22（154 tasks/stage）已完成，且 final verifier receipt 通過。

輸入路徑：

- `InterSubMod/<extraction-root>/full_extraction_receipt.json`
- `InterSubMod/<ranking-root>/full_ranking_receipt.json`
- `InterSubMod/<final-verification-receipt>.json`

命令：

```bash
python InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py \
  --extraction-root /absolute/path/to/full_extraction_root \
  --ranking-root /absolute/path/to/full_ranking_root \
  --final-verification-receipt /absolute/path/to/final_verification.json \
  --output /absolute/path/to/final_numeric_summary.json
```

輸出：

- `/absolute/path/to/final_numeric_summary.json`
- `/absolute/path/to/final_numeric_summary.json.sha256`

若正式 154-task evidence 尚未存在，本工具不會填入或推測真實數字；測試 fixture 的 `SYNTH` 數字只能證明程式契約，不可作研究結果。

## 驗證結果

執行命令：

```bash
python -m unittest \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_full_m2_receipts.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_ranking.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_extraction.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_final_numeric_summary.py -v
```

實際結果：`Ran 102 tests in 18.099s — OK`；`/usr/bin/time -v` 顯示 wall time 18.63 s、peak RSS 70,804 kB、exit status 0。

新增測試涵蓋：正常 bundle、candidate / extraction child / ranking child / runtime 四類 byte tamper、`parent_choice_count=0`、未知 h* role、component k conservation failure、正式 scope gate、exclusive output/sidecar、tied coarse topology consistent/inconsistent。四類 byte tamper 均 fail closed。

runtime TSV 逐列讀取，不保留完整 row；為了 exact nearest-rank p50/p95/p99，會保留 compact float64 arrays，因此是 O(N values)，不是 O(1) memory。此限制已寫入輸出 `memory_model`。

exact-11 release 的 11 個 SHA-256 在實作前後逐一使用 `sha256sum -c` 驗證，結果為 `11/11 OK`；本任務未修改 `build_validated_html_report.py`。
