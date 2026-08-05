<!--
建立時間: 2026-07-18 00:00 +08:00
目標: 記錄 prepared-base 與 complete-only structural cache 的修改、輸入、命令、結果、獨立審查與限制
處理範圍: bounded engineering PROBE；不含真實雙pilot或全154-task執行
關聯檔案:
  - InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/pre-decision-audit.md
  - InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/implementation-notes.md
  - InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/results/verification_summary.json
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260718_M2_exact_preserving修補正式執行Runbook_03.md
-->

# Hypercube exact-preserving修補與驗證結果

用 PREP：**修補與測試均通過，可進入v3雙pilot；不可直接宣稱全量scaling已解決。**
（影響：中；exact-equivalence信心：高；真實scale信心：待pilot）

## 1. 修改重點

1. **Prepared base MILP**
   - universe、vertices、bounds、predecessor rows與reduced group rows每次enumeration只建一次。
   - objective equality與no-good rows仍依原探索順序追加。
   - `time_limit_seconds`仍是每個SciPy `milp()` call的上限。
2. **Complete-only structural cache**
   - 生命週期只限單一ranking child process。
   - full frozen tuple是判等權威；SHA只作diagnostic。
   - 只存`complete=true`且重驗root/full mandatory、objective、connectivity、
     partial coverage與vertex-set digests的結果。
   - store與hit各deep-copy；incomplete／limit／numerical／MAX_SETS一律不存。
3. **Likelihood邊界不變**
   - BQ mixture、fixed-error grid與bootstrap每個unit重新fit。
   - 不加入VAF第二分數，不快取likelihood，不評parent edges。
4. **Receipt diagnostics**
   - run-level記錄lookup/hit/miss/store/reject/solver invocation與跨minread hits。
   - 守恆：`lookups=hits+misses`、`misses=stores+rejected`、
     `solver_invocations=misses`、`entries=stores-evictions`。
   - 科學TSV與runtime TSV既有schema不變。
5. **Release身分**
   - frozen-v2 Runbook與失敗root保持不變。
   - 新增v3 Runbook，綁定本次新ranker與solver bytes；雙pilot先於full。

## 2. 輸入、命令與輸出

### 2.1 修改輸入

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py`
- 對應solver/ranker/runbook tests。

### 2.2 修改前基準命令

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest \
  tests/test_hypercube_exact.py \
  tests/test_build_m2_patterns_and_rank.py -v
```

輸出片段：`Ran 51 tests in 2.638s`、`OK`、exit 0。

24-optima probe輸出：

```json
{"build_problem_calls":25,"complete":true,"objective":3,"vertex_set_count":24}
```

### 2.3 修改後完整命令

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest \
  discover -s tests -p 'test_*.py' -v
```

輸出片段：`Ran 341 tests in 136.479s`、`OK`、exit 0。

輸出紀錄：

- `InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/results/verification_summary.json`
- 本文件。

## 3. 數值驗證結果

| Gate | 修改前 | 修改後 | 判定 |
|---|---:|---:|---|
| `AAAA,k=4`候選數 | 24 | 24 | PASS |
| objective `h*` | 3 | 3 | PASS |
| ordered candidate-ID digest | `a1ec0d…c93a86` | `a1ec0d…c93a86` | PASS |
| base problem builds | 25 | 1 | PASS |
| MILP solves | 25 | 25 | 誠實保留 |
| final full tests | — | 341/341 | PASS |
| cache repeated lookups | — | 3 | — |
| cache miss／hit | — | 1／2 | PASS |
| solver invocations | — | 1 | PASS |
| incomplete stored | — | 0 | PASS |
| likelihood cached | — | false | PASS |

Wall-time probe由0.266119秒到0.173771秒，但只屬單次本機diagnostic，
**不可寫成穩定加速倍數**；正式效能必須由v3雙pilot量測。

## 4. Exact與fail-closed驗證

- small-k exhaustive reduction/no-good oracle：PASS。
- 獨立legacy shadow 51/51：`h*`、`complete`、ordered IDs共0 mismatch。
- first/next limit、next infeasible、max_sets=0/1/2、h*=0：PASS。
- objective/no-good全零sparse rows與duplicate exclusions：PASS。
- cache digest強制碰撞：不同full tuple仍不會false hit。
- cache回傳值mutation：不污染後續hit。
- same structure／different counts：只重用enumeration；likelihood與unit semantic SHA各自不同。
- cache on/off：四個科學輸出semantic SHA及aggregate完全相同。
- v3 exact-11 operational closure：3/3 PASS。

## 5. 三方獨立審查

| 視角 | 結果 | 主要證據 |
|---|---|---|
| solver exactness | PASS | 60/60正式tests、51/51 legacy shadow、7/7 status branches |
| cache／resource contract | PASS | 101/101 ranker＋full verifier tests、無科學schema drift |
| scale red-team | PASS_FOR_PROBE | v2歷史保留、v3新source binding、full 341/341；仍要求雙pilot |

## 6. 結論與限制

本次已證明：

- Python端重複base建模可安全移除；
- 相同完整結構問題可安全重用已證明完整的candidate catalog；
- 科學答案、candidate順序、likelihood與現有TSV schema沒有被改動。

本次尚未證明：

- 33個真實長尾units或chr6一定能在時限內完成；
- SciPy/HiGHS model與presolve被重用；
- full 154可啟動。

因此最終狀態是 **`PASS_FOR_PROBE`**。下一個合法動作是依v3 Runbook建立新frozen
contract並跑HCC1395_DORADO chr6＋H2009 chr2雙pilot；不是resume v2，也不是直接full。
