<!--
建立時間: 2026-07-16 08:21 CST
目標: 獨立稽核 M0 likelihood census 的欄位契約、分母、分類守恆、receipt/TSV 雜湊與 canonical-v5 對應
處理範圍: COLO829 m0_colo829_v3 pilot 稽核 + chr1–22 × 7 datasets m0_full_v4 Task B 最終 red-team
關聯檔案: InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_receipt.json；InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz；InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json
服務目標: G3 / G4 / G5
-->

# M0 likelihood census 獨立稽核（COLO829 v3 pilot + full_v4）

用 SCQA：現有 M0 要把候選 mutation-state vertex sets 依 read-pattern likelihood 分類；但 receipt 自身的守恆旗標不足以證明每列、分母與 canonical population 都一致；因此本稽核不引用執行腳本的彙總函式，從 TSV 與 canonical-v5 獨立重算並逐 unit 對應。

> [!NOTE]
> **文件分層**：§1–§8 保留 COLO829 v3 pilot 稽核與當時的執行紀錄；§9 是後續完成的 chr1–22 × 7 datasets full_v4 最終 red-team，應以 §9 作為目前結論。

**TL;DR：full_v4 通過 7/7 datasets 獨立稽核 — 64,973 個 eligible HP-lineage units 全數與 canonical-v5 對應並深度重建，receipt、TSV、T/V 與排序分類皆守恆（影響: 高，信心: 高）。** 這些結果可作 current-v5 的 M0 engineering baseline；不能升格為真實 clone 數、唯一 parent edge、最終拓撲或 calibrated posterior。完整結果與 claim ceiling 見 §9。

## 1. 稽核結論

整體判定：**PASS with provenance caveats**。

| 檢查 | 結果 | 實際證據 |
|---|---:|---|
| receipt 宣告的 TSV SHA/大小 | PASS | SHA-256 `b3e964011d9ac46fe1a1e118d2694399e83431fa8a98767d9bb0512003df8fad`；2,378,201 bytes |
| TSV schema 26 欄、型別與值域 | PASS | 11,401/11,401 rows 可解析；無重複 `(dataset, region, family)` key |
| TSV aggregate 對 receipt | PASS | HP units、regions、ΣT、ΣV、status counts、fractions、by-dataset 全相等 |
| TSV eligible units 對 canonical-v5 | PASS | missing=0；extra=0 |
| canonical eligible＋excluded 強守恆 | PASS | 11,401 eligible + 751 capped = 12,152 primary mutation HP units |
| 所有候選 vertex-set classes 深度重建 | PASS | 11,401/11,401；stored 10,633、frozen-solver rebuilt 768 |
| T/V 分割 | PASS | T=1/V=1 4,290；T>1/V>1 7,111；不可能的 T=1/V>1 為 0 |
| 最終 selection-status 分割 | PASS | 4,290 + 6,572 + 438 + 101 = 11,401 |
| optimizer fail-closed | PASS | 101 個未收斂 unit 全數標為 `RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE` |

## 2. 資料單位與分母

### 2.1 三種不同單位

1. **Region**：canonical `region`；同一 region 可有 HP1、HP2 兩個 mutation-bearing lineage。
2. **HP-lineage unit（M0 主分母）**：`dataset × region × family(HP1/HP2)`；每個 TSV row 就是一個 unit。
3. **Candidate count**：一個 HP-lineage unit 內的候選數，不能和 region/HP unit 數量相加比較。

因此「11,401」是 eligible HP-lineage units，不是 region 數；「85,013」是把每個 unit 的候選 T 加總，也不是 85,013 個 region。

### 2.2 Population funnel

| Population | 數量 | 分母 | 比例 | 定義 |
|---|---:|---:|---:|---|
| Primary mutation regions | 7,878 | 7,878 | 100.000% | 至少一個 primary、mutation-bearing lineage 的 region |
| Regions with any eligible lineage | 7,631 | 7,878 regions | 96.865% | 至少一個 lineage 通過 M0 filter |
| Fully M0-eligible regions | 7,128 | 7,878 regions | 90.480% | region 內所有 primary mutation lineages 都 eligible |
| Primary mutation HP-lineage units | 12,152 | 12,152 units | 100.000% | 進入 eligibility 判斷的 HP units |
| Excluded: capped | 751 | 12,152 units | 6.180% | `capped=true`，不可進 M0 exact candidate census |
| Eligible M0 HP-lineage units | 11,401 | 12,152 units | 93.820% | noncapped、candidate complete、solver verification full PASS、具有 retained patterns |

強守恆式：

```text
12,152 primary mutation HP units
= 11,401 eligible output rows
+    751 excluded capped units
```

## 3. T、V、E 與排序分類的精確定義

| 欄位 | 單位 | 定義 | 不代表什麼 |
|---|---|---|---|
| `T` | candidates / HP-lineage unit | canonical raw tree candidate 數；receipt 的 `ΣT` 是跨 row 加總 | 不是 region 數、clone 數 |
| `V` | distinct vertex sets / unit | 把具有相同 mutation-state vertex set 的 edge variants 合併後數量 | 不保證每個 V 是真實 cell clone 組合 |
| `E` | edge candidates / unit | 此 schema 記錄 raw edge candidate 數，契約要求 `E=T` | snapshot likelihood 不能據此辨識 parent edge |
| `best_vertex_sets` | vertex sets / unit | log-likelihood 與第一名差距 ≤ `1e-6` 的並列第一 V 數 | 不是 posterior credible set |
| `top_relative_likelihood_weight` | 0–1 | 在枚舉 V 中用 relative likelihood 正規化後，單一最高候選的 weight | 不是校準 posterior probability |

COLO829 的特定觀察是 `ΣT=85,013` 且 `ΣV=85,013`，逐列也全部 `T=V`。其中 `V` 是 distinct vertex sets 的數量；這表示**這份輸出沒有出現「相同 vertex set N、只有 edge 不同」的候選重複**。它是資料結果，不是方法上的普遍保證。

## 4. COLO829 分類結果與雙分母比例

### 4.1 結構候選層

| 分類 | HP units | 占 eligible 11,401 | 占 T>1 7,111 |
|---|---:|---:|---:|
| T=1 且 V=1 | 4,290 | 37.628% | — |
| T>1 且 V=1 | 0 | 0.000% | 0.000% |
| T>1 且 V>1 | 7,111 | 62.372% | 100.000% |
| T=1 且 V>1（不可能） | 0 | 0.000% | — |

### 4.2 最終互斥狀態層

| 最終狀態 | HP units | 占 eligible 11,401 | 占 T>1 7,111 | 判讀 |
|---|---:|---:|---:|---|
| `T1_CANDIDATE_UNIQUE` | 4,290 | 37.628% | — | 原始候選只有一個；沒有做候選間排序 |
| `LIKELIHOOD_TIED_VERTEX_SETS` | 6,572 | 57.644% | 92.420% | 多個 V 位於 `1e-6` tie tolerance 內 |
| likelihood unique vertex set | 438 | 3.842% | 6.159% | 唯一第一；本資料全部同時是 single-edge candidate |
| `RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE` | 101 | 0.886% | 1.420% | SLSQP 至少一個 candidate fit 未收斂；安全棄權 |
| **合計** | **11,401** | **100.000%** | — | 互斥且完整 |

T>1 的獨立守恆式：

```text
7,111 T>1/V>1 units
= 6,572 likelihood tied
+    438 likelihood unique
+    101 optimizer nonconverged / abstain
```

重要：receipt 的 `all_pass=true` 與 `aggregate.all_fits_converged=false` 可以同時成立。前者表示檔案/守恆與 fail-closed 規則通過；後者明確表示存在 101 個 optimizer 未收斂 unit。報告時不能把這 101 個併入 unique 或 tied 結論。

## 5. 逐欄與逐列驗證內容

獨立 verifier 沒有 import `run_m0_likelihood_census.py`，而是另外實作下列檢查：

1. TSV header 必須精確等於 26 欄 schema。
2. 字串、整數、布林、nullable likelihood、JSON SHA-list 全部做型別與值域檢查。
3. 每列驗證 `1 ≤ V ≤ T`、`E=T`、`1 ≤ best_V ≤ V`、`best_V ≤ top_E ≤ T`。
4. `V=1` 必須沒有 best/second/delta likelihood，weight=1、rank=0；`V>1` 必須有有限數值，並驗證 `delta=best-second`。
5. 依 `T/V/best_V/top_E/converged` 獨立導出 status，和 TSV status 逐列比較。
6. `vertex_set_ids` 數量必須等於 V；top IDs 必須等於 best_V 且為全集子集。
7. 重新讀 canonical-v5，依 eligibility 規則建立預期 `(dataset, region, family)` keys，逐 row 比對座標、k、read/pattern counts、T/E。
8. 對 stored trees 自行解析 vertex sets；對 display-incomplete units 呼叫已雜湊的 canonical frozen solver 重建完整可行集合；11,401 個 unit 的 V、完整 ID 序列及 top edge count 全部吻合。
9. 從 TSV 重新計算 receipt aggregate、fractions、by-dataset 與所有互斥分割。

數值 likelihood 本身沒有為 11,401 rows 全部重跑 SLSQP；本稽核驗證的是已輸出數值的有限性、排序關係、delta、tie/unique/nonconvergence status 契約。若要驗證 optimizer bit-level reproducibility，必須另外保存並鎖定執行環境後重跑 scoring，不能只靠現有 receipt 1.0.0 證明。

## 6. 輸入、命令、輸出與實際片段

### 6.1 輸入

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_colo829_v3/m0_receipt.json`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_colo829_v3/m0_hp_lineage_likelihood_census.tsv.gz`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/samples/COLO829/layered_region_view_COLO829.json`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/source_bundle/files/imported/003_tree_enumeration_solver.py`

### 6.2 執行命令

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/usr/bin/time -v python3 \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_m0_likelihood_census.py \
  --output-dir research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_colo829_v3 \
  --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 \
  --candidate-mode deep \
  --json-output research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_colo829_v3/independent_verification_v2.json
```

測試命令：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest -v \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_m0_likelihood_census.py
```

### 6.3 輸出

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_colo829_v3/independent_verification_v2.json`
- 本稽核文件：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M0_likelihood_census獨立稽核_01.md`

實際輸出片段：

```json
{
  "verdict": "PASS",
  "checks": {
    "n_errors": 0,
    "receipt_tsv_hash_matches": true,
    "receipt_aggregate_matches_tsv": true,
    "tsv_units_match_canonical_eligible_population": true,
    "strong_eligible_excluded_conservation": true,
    "selection_partition_conserves": true,
    "T_V_partition_conserves": true,
    "required_files_are_exact_byte_hashed": true,
    "verification_runtime_versions_recorded": true
  },
  "canonical_reconciliation": {
    "n_candidate_units_deep_checked": 11401,
    "missing_tsv_units": 0,
    "extra_tsv_units": 0
  }
}
```

實際測試片段：

```text
Ran 4 tests in 0.010s
OK
```

verification 1.1.0 綁定的 runtime：CPython 3.9.12、NumPy 1.23.5、SciPy 1.13.1；`m0_rows_gzip` 的 hash scope 明定為「包含 gzip container 的 exact compressed file bytes」。它同時保存 M0 receipt、rows、census script、scoring utility、independent verifier 五個檔案的 path、size 與 SHA-256。

## 7. SHA-256 證據鏈

| Artifact | SHA-256 |
|---|---|
| COLO829 receipt | `88f85734eae301911b207a6f676a40133429685ee456094e572f145509b54882` |
| COLO829 TSV.gz | `b3e964011d9ac46fe1a1e118d2694399e83431fa8a98767d9bb0512003df8fad` |
| Independent verification JSON 1.1.0 | `fb3b4e9ff05b3be5c89e951b2fdf2b5b6cfb863ed20316c4e075922b01a7ce2f` |
| Running M0 census script（稽核當下） | `e2a3742e645a16d7254f950dee61de420662e0c6d13f956dc8338b5779535acd` |
| Independent verifier | `153ea85e19a70195cb820f7d4af01dfed278ca66a99388f7f6de57f94c195c9f` |
| Independent verifier test | `2123eccc4917f04c9925daceb92dbf6a98c7b8ef9c7d855b29cf8c2e511d3381` |
| Current scoring utility `hypercube_exact.py` | `db85c94f5cd1fda5bb1d4eca9ca64d5df6680829b1b005ac02be9d3c835500bc` |
| canonical frozen solver | `36727f4e1d8d7ce8abf869606211c93d8c1a0506dd7d142e855863c412ca0d61` |
| canonical COLO829 JSON | `438950cf287e18802065bd67ce12e0ea5b6fc3b77eaac96e9fbd71543b54cd32` |

## 8. 問題分級與後續要求

### 無 blocker

未發現會使 COLO829 v3 計數失效的 blocker；不需中止或修改正在執行的 full_v4 census。

### Caveats

1. **中度／高信心 — 原始 M0 receipt 的歷史執行 provenance 不完整**：receipt 1.0.0 有 input、frozen solver、output hashes，但沒有記錄 `run_m0_likelihood_census.py`、`hypercube_exact.py` 的執行時 SHA，亦沒有 Python/NumPy/SciPy 版本。verification 1.1.0 已綁定驗證當下的 source/runtime；這可精確識別驗證環境，但不能反向證明先前 census process 載入的歷史 bytes/runtime 完全相同。
2. **中度／高信心 — `all_pass` 不是「所有 optimizer 收斂」**：101 個未收斂 unit 已 fail-closed，資料契約通過，但生物結論必須排除或另列。
3. **低度／高信心 — receipt naming/semantic weakness**：`eligible_plus_excluded_hp_units_positive` 的原始實作只檢查輸出非空與 excluded 非負，沒有檢查等式。verifier 只重現它來確認 legacy receipt 自洽，**不把它當成真實 checked invariant**；真正使用的強 invariant 是另算的 `eligible + excluded == canonical primary units`，本次結果為 true。
4. **中度／高信心 — M0 claim ceiling**：本資料是 thresholded alignment-exposure engineering census，不是 lossless molecule likelihood；不能升格為 clone count、parent-edge identification 或 calibrated posterior。

### full_v4 驗證命令與完成狀態

歷史執行註記：2026-07-16 08:22 CST 時 full_v4 尚在執行；TSV/receipt 後於 08:40 完成，independent verification 於 08:45 完成且 PASS。下列是實際使用的一鍵驗證命令；§9 再以唯讀 red-team 獨立重算最終表格。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_m0_likelihood_census.py \
  --output-dir research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4 \
  --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 \
  --candidate-mode deep \
  --json-output research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json
```

full_v4 的下列條件已全部成立：

- `verdict=PASS`、`n_errors=0`；
- `selected_datasets` 精確為 7 datasets 且 `full_task_b_scope=true`；
- receipt/TSV SHA 相符；
- missing/extra canonical units 都為 0；
- `n_candidate_units_deep_checked == n_hp_lineage_units`；
- T/V 與 final status 兩個 partition 都守恆；
- nonconverged 保持獨立棄權分母，不混入 unique/tied。

## 9. full_v4 最終 red-team（chr1–22 × 7 datasets）

用 Claim–Evidence–Boundary：先確認哪些數字由 exact-byte artifacts 與 canonical reconciliation 支持，再明確限制能下到哪一層結論，避免把 HP-lineage engineering census 誤寫成 region、clone 或真實拓撲結果。

**TL;DR：full_v4 verification 1.1.0 為 PASS、errors=0；64,973/64,973 units 深度重建，missing=0、extra=0。** 最安全的正式敘述是：「在 current-v5 M0 eligible HP-lineage units 中，26,225 為 T=1；38,748 為 T>1，其中 8,751 有唯一第一名 mutation-state vertex set、27,270 並列、2,626 optimizer 棄權、101 只有 edge variants。」這不是 64,973 個 regions，也不是已確認 8,751 個真實拓撲。

### 9.1 最終驗證判定與證據身份

| Gate | 結果 | 實際證據 |
|---|---:|---|
| Scope | PASS | 7/7 technical datasets；TSV 只含 chr1–chr22 |
| Gzip integrity / row count | PASS | `gzip -t` PASS；64,974 lines = 1 header + 64,973 rows |
| Exact-byte receipt / rows binding | PASS | receipt SHA `eba081…`; gzip rows SHA `9df74d…` |
| TSV aggregate 對 receipt | PASS | `receipt_aggregate_matches_tsv=true` |
| TSV eligible keys 對 canonical-v5 | PASS | missing=0；extra=0 |
| Candidate classes 深度重建 | PASS | 64,973/64,973 units；stored=61,702、rebuilt=3,271 |
| T/V 與 final-status partition | PASS | 兩個 partition 均守恆 |
| Verification-time runtime/source binding | PASS | CPython 3.9.12、NumPy 1.23.5、SciPy 1.13.1；census/scoring/verifier SHA 已保存 |
| Numeric optimizer 全重跑 | **未執行** | verifier 檢查 persisted ordering/status，但沒有重新執行全部 SLSQP fits |

Source/runtime binding 描述的是 verification-time environment；由於原始 M0 receipt 1.0.0 沒有保存這些欄位，它不能單獨反向證明 census process 當時載入的 historical bytes/runtime 完全相同。這是 provenance caveat，不影響 receipt、rows、canonical keys 與分類守恆的 exact-byte 驗證。

輸入 artifact：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_receipt.json`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json`

SHA-256：

| Artifact | Size | SHA-256 |
|---|---:|---|
| full_v4 receipt | 9,177 bytes | `eba081a70f16c008f70cd97c85ec3bcbce41d3982eb6a55c5915e89149197699` |
| full_v4 TSV.gz exact compressed bytes | 12,292,887 bytes | `9df74db30bc930fee4e0a6941f371bfe12b069808423bb53be5eaf3fc77c1a6c` |
| Independent verification 1.1.0 | 13,257 bytes | `912760104cb3b7dca4e18d56ac429f0cbc901c81800719425ebaf228c2ae735c` |

### 9.2 Region 與 HP-lineage 分母不可混用

M0 的分類主單位是 `dataset × region × family(HP1/HP2)`，即 **HP-lineage unit**。一個 region 可以貢獻一個或兩個 eligible HP-lineage units，而且兩個 HP 可以落在不同分類。因此不能把 64,973、26,225、38,748、8,751 或 27,270 改稱為 region 數。

| Funnel 層 | 數量 | 分母 | 比例 | 正確意義 |
|---|---:|---:|---:|---|
| Primary mutation regions | 50,215 | 50,215 regions | 100.00% | 至少一個 primary mutation-bearing lineage 的 current-v5 region |
| Fully M0-eligible regions | 42,240 | 50,215 regions | 84.12% | region 內所有 primary mutation lineages 都通過 M0 eligibility |
| Regions with any eligible HP lineage | 46,385 | 50,215 regions | 92.37% | 至少一個 HP lineage eligible；M0 rows 涵蓋的 region union |
| Regions with no eligible lineage | 3,830 | 50,215 regions | 7.63% | 沒有任何 lineage 進入 M0 rows |
| Mixed-eligibility regions | 4,145 | 50,215 regions | 8.25% | 有 eligible lineage，但並非所有 lineage eligible |
| Primary mutation HP-lineage units | 72,994 | 72,994 units | 100.00% | eligibility 判斷主母體 |
| Excluded capped HP-lineage units | 8,021 | 72,994 units | 10.99% | 這是 HP-unit 數，不是 capped region 數 |
| Eligible M0 HP-lineage units | 64,973 | 72,994 units | 89.01% | TSV rows 與後續分類主分母 |

強守恆：

```text
72,994 primary mutation HP-lineage units
= 64,973 eligible M0 rows
+  8,021 capped/excluded HP-lineage units
```

64,973 eligible units 分布在 46,385 regions，平均 1.4007 eligible HP units/region。這個平均只用來說明單位差異，不代表平均 clone 數。

### 9.3 T/V 的最終定義與守恆

| 符號 | 精確單位與定義 | full_v4 加總 | Claim boundary |
|---|---|---:|---|
| `T` | 每 HP-lineage unit 的 raw edge-tree candidates；`ΣT` 是跨 unit 加總 | 444,007 | 不是 regions 或 clones |
| `V` | 把相同 mutation-state vertex set 的 edge variants 合併後，每 unit 的 distinct vertex sets；`ΣV` 跨 unit 加總 | 443,745 | 不等於真實 cellular populations |
| `T−V` | 被 V-collapse 合併的額外 edge candidate instances | 262 | 不是 262 個 units；實際 T>V units 是 153 |

逐 HP-unit 的互斥分類：

| T/V class | HP units | 占 eligible 64,973 | 占 T>1 38,748 | 解釋 |
|---|---:|---:|---:|---|
| T=1、V=1 | 26,225 | 40.36% | — | enumeration 只有一個 candidate，未做候選間 likelihood 比較 |
| T>1、V=1 | 101 | 0.16% | 0.26% | 多個 edge variants 共享同一 vertex set；mutation-state likelihood 無法辨識 edge |
| T>1、V>1 | 38,647 | 59.48% | 99.74% | 存在多個可比較 mutation-state vertex sets |
| T=1、V>1 | 0 | 0.00% | — | 數學上不允許，實際也為 0 |
| **合計** | **64,973** | **100.00%** | — | partition 完整 |

`ΣT−ΣV=262` 只占 raw T 的 0.059%，但出現在 153 個 units。當中 101 個是 T>1/V=1；其餘 units 是 V>1、但至少一個 vertex set 含多個 edge variants。不能以 `262` 取代「edge-unresolved units」數量。

### 9.4 T>1 units 的 likelihood 結果

M0 固定參數為 error rate=`0.01`、log-likelihood tie tolerance=`1e-6`；scoring grain 是 distinct mutation-state vertex set N。相同 N 的 edge variants 依政策取得相同分數；`V` 僅表示不同 N 的數量。

| 最終互斥狀態 | HP units | 占 eligible 64,973 | 占 T>1 38,748 | 可以說什麼 |
|---|---:|---:|---:|---|
| T>1、V=1 edge-equivalent | 101 | 0.16% | 0.26% | mutation states 唯一，但 parent edges 未辨識 |
| Likelihood unique V | 8,751 | 13.47% | 22.58% | 在枚舉 candidates 與 M0 模型內有唯一第一名 V |
| Likelihood tied V | 27,270 | 41.97% | 70.38% | 至少兩個 V 位於 `1e-6` tie tolerance 內 |
| Optimizer nonconverged / abstain | 2,626 | 4.04% | 6.78% | fail-closed；不可併入 unique 或 tied |
| **T>1 合計** | **38,748** | **59.64%** | **100.00%** | partition 完整 |

8,751 個 unique-V 又分成 8,737 個 top single-edge candidate 與 14 個 top edge-unresolved。因此正式文字應寫「唯一第一名 mutation-state vertex set」，不能直接縮寫成「唯一真實拓撲」。

另一個合法但不同的分母是「已收斂且 V>1」的 36,021 units：unique V=8,751（24.29%）、tied V=27,270（75.71%）。這組比例排除了 101 個 V=1 edge-only units 與 2,626 個 optimizer abstain；不得與「占全部 eligible」或「占 T>1」比例混寫。

### 9.5 每 technical dataset 的結果與雙分母

> HCC1395 與 HCC1395_DORADO 是同一生物樣本的兩套技術資料／處理結果。因此下表是 **7 technical datasets、6 biological samples**，不能把兩者視為獨立腫瘤個體來做生物統計。

下表 T1/T>1 百分比使用各 dataset 的 eligible HP units；edge/unique/tied/abstain 百分比改用該 dataset 的 T>1 units。表頭明示分母，避免 denominator shift。

| Dataset | Regions with any eligible HP | Eligible HP units | T=1 n（占全部） | T>1 n（占全部） | V=1 edge-only（占 T>1） | Unique V（占 T>1） | Tied V（占 T>1） | Abstain（占 T>1） |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| COLO829 | 7,631 | 11,401 | 4,290（37.63%） | 7,111（62.37%） | 0（0.00%） | 438（6.16%） | 6,572（92.42%） | 101（1.42%） |
| H1437 | 8,196 | 11,628 | 4,393（37.78%） | 7,235（62.22%） | 13（0.18%） | 1,005（13.89%） | 5,630（77.82%） | 587（8.11%） |
| H2009 | 7,684 | 10,599 | 3,117（29.41%） | 7,482（70.59%） | 31（0.41%） | 1,026（13.71%） | 5,569（74.43%） | 856（11.44%） |
| HCC1395 | 7,505 | 10,120 | 4,217（41.67%） | 5,903（58.33%） | 21（0.36%） | 2,256（38.22%） | 3,268（55.36%） | 358（6.06%） |
| HCC1395_DORADO | 7,360 | 10,028 | 4,033（40.22%） | 5,995（59.78%） | 2（0.03%） | 2,125（35.45%） | 3,329（55.53%） | 539（8.99%） |
| HCC1937 | 3,491 | 4,447 | 2,490（55.99%） | 1,957（44.01%） | 28（1.43%） | 874（44.66%） | 918（46.91%） | 137（7.00%） |
| HCC1954 | 4,518 | 6,750 | 3,685（54.59%） | 3,065（45.41%） | 6（0.20%） | 1,027（33.51%） | 1,984（64.73%） | 48（1.57%） |
| **全部 7 datasets** | **46,385** | **64,973** | **26,225（40.36%）** | **38,748（59.64%）** | **101（0.26%）** | **8,751（22.58%）** | **27,270（70.38%）** | **2,626（6.78%）** |

每列 T>1 內的四分類都加總為 100%。例如 H2009 的 nonconvergence 最高（11.44%），COLO829 的 tied 比例最高（92.42%），HCC1937 的 unique-V 比例最高（44.66%）；這些是方法行為差異，尚不能歸因為腫瘤演化差異。

在只看「已收斂且 V>1」時，各 dataset 的 unique/tied 比例如下：

| Dataset | Converged V>1 denominator | Unique V | Tied V |
|---|---:|---:|---:|
| COLO829 | 7,010 | 438（6.25%） | 6,572（93.75%） |
| H1437 | 6,635 | 1,005（15.15%） | 5,630（84.85%） |
| H2009 | 6,595 | 1,026（15.56%） | 5,569（84.44%） |
| HCC1395 | 5,524 | 2,256（40.84%） | 3,268（59.16%） |
| HCC1395_DORADO | 5,454 | 2,125（38.96%） | 3,329（61.04%） |
| HCC1937 | 1,792 | 874（48.77%） | 918（51.23%） |
| HCC1954 | 3,011 | 1,027（34.11%） | 1,984（65.89%） |
| **全部** | **36,021** | **8,751（24.29%）** | **27,270（75.71%）** |

### 9.6 其他 M0 engineering 指標

| 指標 | 數量／比例 | 判讀限制 |
|---|---:|---|
| Stored-complete candidate units | 61,702 / 64,973（94.97%） | 表示輸出儲存路徑，不是證據品質較高 |
| Frozen-solver reconstructed units | 3,271 / 64,973（5.03%） | 重建完整 analytical candidate set；不是只分析 stored prefix |
| Scoring alignment exposures | 7,086,491 | 跨 HP units 加總；不是 unique molecules |
| Reported alignment exposures | 7,341,571 | 跨 HP units 加總；不是 tumor read depth |
| M0 scoring exposure retention | 96.53% | 省略 255,080（3.47%）thresholded-out exposures；不是 lossless likelihood coverage |
| All fits monotone | true | 不代表全部收斂；仍有 2,626 abstain |

### 9.7 哪些數字可進正式報告

#### A. 可作正式主流程／資料品質數字

可直接放入正式報告，但要保留 current-v5、chr1–22、7 technical datasets 的 scope：

- 50,215 primary mutation regions；
- 72,994 primary mutation HP-lineage units；
- 42,240 fully M0-eligible regions、46,385 regions with any eligible lineage；
- 8,021 capped/excluded HP-lineage units、64,973 eligible HP-lineage units；
- independent verification PASS、missing=0、extra=0、64,973/64,973 candidate units deep checked；
- receipt/TSV/verification 的 exact SHA-256 證據鏈。

#### B. 可放正式方法結果，但標題必須寫「M0 engineering baseline」

- T/V partition：26,225 T=1/V=1、101 T>1/V=1、38,647 T>1/V>1；
- T>1 的 M0 排序：8,751 unique V、27,270 tied V、2,626 nonconverged abstain、101 edge-equivalent；
- 各 technical dataset 的相同分類與雙分母比例；
- `ΣT=444,007`、`ΣV=443,745`、153 units 有 T>V、`Σ(T−V)=262`；
- stored/rebuilt 與 alignment-exposure retention 指標。

這些數字可以說明「目前候選空間有多大、M0 能區分多少、多少仍並列或棄權」，不能作最終生物結論。

#### C. 不可用 M0 數字宣稱的結論

| 不允許的說法 | 原因 | 安全改寫 |
|---|---|---|
| 「64,973 個 regions」 | 64,973 是 HP-lineage units | 「64,973 eligible HP-lineage units，分布於 46,385 regions」 |
| 「8,751 個唯一真實拓撲」 | unique 的是 M0 枚舉內第一名 V；14 個 top edge 仍不唯一 | 「8,751 個 T>1 units 有唯一第一名 mutation-state vertex set」 |
| 「26,225 個經 VAF/likelihood 確認的樹」 | T=1 沒有候選間 ranking；M0 也不是獨立 VAF 矯正 | 「26,225 個 units 的原始 candidate set 為 T=1」 |
| 「27,270 個 biologically ambiguous regions」 | 單位是 HP lineage，且 ambiguity 依固定 error/tie/model | 「27,270 個 M0 HP-lineage units 的 top V 並列」 |
| 「已辨識 parent–child edge」 | likelihood 觀察 state，不觀察 transition；相同 vertex set N 的 edges 同分 | 「parent edge 在此模型下通常不可識別」 |
| 「unique V 是 calibrated posterior probability」 | relative likelihood weight 不是 posterior | 「在枚舉 candidates 內依 M0 likelihood 排名」 |
| 「7 個獨立癌症樣本比較」 | HCC1395 與 DORADO 是同一生物樣本 | 「7 technical datasets、6 biological samples」 |
| 「M0 已完成 lossless partial-read likelihood」 | M0 使用 thresholded alignment exposures | 「M0 是 engineering census；lossless molecule/PS-aware 模型需看後續 M2」 |

### 9.8 Claim ceiling

full_v4 能支持的最高層級結論是：

> 在 current-v5 chr1–22 的 7 technical datasets 中，M0 對 64,973 個 eligible HP-lineage units 完成 mutation-state candidate census。38,748 個 units 原始 T>1；在固定 error rate 0.01、tie tolerance `1e-6` 的 thresholded read-pattern mixture 下，8,751 個有唯一第一名 vertex set、27,270 個並列、2,626 個 optimizer 棄權，另有 101 個只有 edge variants。所有數字均為 engineering ambiguity baseline。

不能越過的界線：**不是真實 clone count、不識別唯一 parent edge、不等於已確認癌症演化拓撲、不提供 calibrated posterior，也不是最終 publication likelihood。**

### 9.9 唯讀 red-team 命令與實際片段

完整輸入路徑、hash 與 verification gate：

```bash
sha256sum \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_receipt.json \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json

gzip -t \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz

jq -e '
  .verdict == "PASS" and
  .scope.full_task_b_scope == true and
  .checks.n_errors == 0 and
  .canonical_reconciliation.n_candidate_units_deep_checked == 64973 and
  .canonical_reconciliation.missing_tsv_units == 0 and
  .canonical_reconciliation.extra_tsv_units == 0
' research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json
```

從 TSV 獨立重算的核心 AWK 邏輯：

```bash
gzip -dc \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz |
awk -F '\t' '
  NR == 1 { next }
  {
    dataset = $1
    key = dataset SUBSEP $2
    if (!(key in seen_region)) { seen_region[key] = 1; regions[dataset]++ }
    units[dataset]++
    sumT[dataset] += $11
    sumV[dataset] += $12
    if ($11 == 1 && $12 == 1) t1[dataset]++
    else if ($11 > 1 && $12 == 1) edge_only[dataset]++
    else if ($11 > 1 && $12 > 1) multiV[dataset]++
    status[dataset, $21]++
  }
  END {
    for (dataset in units)
      print dataset, regions[dataset], units[dataset], t1[dataset],
            edge_only[dataset], multiV[dataset], sumT[dataset], sumV[dataset]
  }
'
```

實際輸出片段：

```text
gzip_test=PASS
64974                         # 1 header + 64,973 rows
verification_gate=PASS

ALL regions=46385 eligible_units=64973
T1=26225 T>1=38748
T>1/V=1=101 uniqueV=8751 tiedV=27270 abstain=2626
sumT=444007 sumV=443745
```

最終 red-team 判定：**Ready to share as a clearly labeled M0 engineering baseline；not ready to share as final biological topology evidence。**
